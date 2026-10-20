function [best, info] = mwspm_calib_fit(pre, C, opts)
%MWSPM_CALIB_FIT  Bounded global-ish fit of the calibration objective.
%
%   [best, info] = mwspm_calib_fit(pre, C)
%   [best, info] = mwspm_calib_fit(pre, C, struct('n_screen', 400, 'seeds', [13 47]))
%
%   C fields
%     free         cellstr, parameters to optimise (must NOT contain d_weeks)
%     fixed        struct, values for every other canonical parameter
%     box          struct of [lo hi] (mwspm_calib_theta('box'))
%     family       likelihood family passed to mwspm_calib_nll
%     family_opts  struct passed through (e.g. struct('theta', 5))
%     d_list       candidate integer-valued d_weeks (default 0:3)
%     sub_weeks    scoring subset of weeks (default all)
%     phase_penalty, target_peak, pen
%     kernel_opts  explicit propagation options; absent preserves Euler replay
%
%   opts fields
%     warm         cell array of theta structs (kept as guaranteed starts)
%     n_screen     uniform draws PER SEED evaluated before the local stage
%     seeds        random seeds for stage-1 sampling (default opts.seed)
%     n_seed       best draws promoted to Nelder-Mead (default 6)
%     restarts     Nelder-Mead restarts from each local solution (default 1)
%     seed, maxFunEvals, maxIter, TolX, TolFun, verbose
%
%   Why two stages, and why several seeds
%   -------------------------------------
%   The likelihood surface of this problem is FLAT and RIDGED: the spread of
%   the negative log-likelihood over the whole admissible rho range is of the
%   same order as the convergence noise of a single Nelder-Mead run.  A single
%   local fit therefore reports "an optimum", not "the optimum".  The objective
%   costs ~0.4 ms, so drawing a few thousand uniform points and promoting only
%   the best few to local refinement is cheap.
%
%   Stage-1 screening is exhaustive only in the weak sense of a random sample
%   and is labelled as such.  Because the sample is what stage 2 refines,
%   stage 1 is drawn from SEVERAL seeds and the best value attained by each
%   seed is reported separately (info.screen_best_by_seed).  That spread is a
%   SCREEN-NOISE diagnostic, NOT a confidence interval and NOT a measure of
%   optimiser reliability: a screen of a few hundred uniform points in 7
%   dimensions is simply far from the optimum, so the spread is large for any
%   problem with a peaked likelihood.  What bounds optimiser reliability is
%   instead (i) info.final_best_by_start, the spread of the FINAL values reached
%   from the different stage-1 starts, (ii) info.n_stage2_within_1nll, how many
%   independent local runs land within 1 nll unit of the reported optimum, and
%   (iii) info.best_start_J, the best value among the caller's guaranteed
%   starts - a final optimum WORSE than this would prove the refinement had
%   been counter-productive.
%
%   d_weeks is NOT continuous: the kernel rounds it, so the objective is a
%   staircase in d.  The fit therefore LOOPS over d in d_list and optimises the
%   remaining coordinates for each d.

if nargin < 3 || isempty(opts)
    opts = struct();
end

names = cellstr(string(C.free));
assert(~any(strcmp(names, 'd_weeks')), 'mwspm:calib', ...
    'd_weeks must not be a continuous free parameter - it is looped over');
if ~isfield(C, 'd_list') || isempty(C.d_list)
    C.d_list = 0:3;
end
if ~isfield(C, 'family'), C.family = 'poisson'; end
if ~isfield(C, 'family_opts'), C.family_opts = struct(); end
if ~isfield(C, 'sub_weeks'), C.sub_weeks = []; end
if ~isfield(C, 'y_obs'), C.y_obs = []; end
if ~isfield(C, 'phase_penalty'), C.phase_penalty = false; end
if ~isfield(C, 'target_peak'), C.target_peak = 37; end
if ~isfield(C, 'pen'), C.pen = 300.0; end

% scoring options shared by every evaluation (sub-week block + observation
% vector).  Built once so the fast path cannot drift from the recorded one.
kOpts = struct();
if isfield(C, 'kernel_opts'), kOpts = C.kernel_opts; end
kOpts.sub_weeks = C.sub_weeks;
if ~isempty(C.y_obs)
    kOpts.y_obs = C.y_obs;
end

if ~isfield(opts, 'warm'), opts.warm = {}; end
if ~isfield(opts, 'n_screen'), opts.n_screen = 400; end
if ~isfield(opts, 'n_seed'), opts.n_seed = 6; end
if ~isfield(opts, 'restarts'), opts.restarts = 1; end
if ~isfield(opts, 'seed'), opts.seed = 13; end
if ~isfield(opts, 'seeds') || isempty(opts.seeds), opts.seeds = opts.seed; end
if ~isfield(opts, 'maxFunEvals'), opts.maxFunEvals = 6000; end
if ~isfield(opts, 'maxIter'), opts.maxIter = 6000; end
if ~isfield(opts, 'TolX'), opts.TolX = 1e-8; end
if ~isfield(opts, 'TolFun'), opts.TolFun = 1e-10; end
if ~isfield(opts, 'verbose'), opts.verbose = false; end

seedList = double(opts.seeds(:))';
nFree = numel(names);

lo = zeros(1, nFree);
hi = zeros(1, nFree);
for k = 1:nFree
    b = C.box.(names{k});
    lo(k) = b(1);
    hi(k) = b(2);
end

% ---------------------------------------------------------------- stage 1
% guaranteed starts: the caller's warm thetas, plus the box centre
uSeed = zeros(0, nFree);
for k = 1:numel(opts.warm)
    th = mergeTheta(opts.warm{k}, C.fixed);
    uSeed(end + 1, :) = cell2mat(cellfun(@(f) double(th.(f)), names, 'UniformOutput', false))'; %#ok<AGROW>
end
uSeed(end + 1, :) = (lo + hi) / 2;
nGuaranteed = size(uSeed, 1);

Uall = uSeed;
owner = zeros(0, 1);
for s = seedList
    rng(s, 'mt19937ar');
    Uall = [Uall; lo + (hi - lo) .* rand(opts.n_screen, nFree)]; %#ok<AGROW>
    owner = [owner; repmat(s, opts.n_screen, 1)];                %#ok<AGROW>
end

% NOTE: nEval is shared with the nested functions below, so it MUST be
% initialised before any nested call.  Stage-1 screening is counted too, and
% that is deliberate: info.n_eval then reports the true objective count.
nEval = 0;
Jscreen = zeros(size(Uall, 1), 1);
for i = 1:size(Uall, 1)
    Jscreen(i) = evalU(Uall(i, :), C.d_list(1));
end

startsPerSeed = zeros(numel(seedList), 1);
for k = 1:numel(seedList)
    sel = owner == seedList(k);
    if any(sel)
        startsPerSeed(k) = min(Jscreen(sel));
    end
end
screenBestBySeed = [seedList(:), startsPerSeed];

[~, order] = sort(Jscreen);
keep = order(1:min(nGuaranteed + opts.n_seed, numel(order)));
keep = unique([1:nGuaranteed, keep(:)'], 'stable');
zStarts = cell(size(keep));
for i = 1:numel(keep)
    zStarts{i} = packU(Uall(keep(i), :), lo, hi);
end

% ---------------------------------------------------------------- stage 2
trace = struct('d_weeks', {}, 'start', {}, 'stage', {}, 'J', {}, 'nll', {}, ...
    'exitflag', {}, 'funcCount', {}, 'iterations', {});
best = struct('J', Inf);

fmsOpt = optimset('MaxFunEvals', opts.maxFunEvals, 'MaxIter', opts.maxIter, ...
                  'TolX', opts.TolX, 'TolFun', opts.TolFun, 'Display', 'off');

for di = 1:numel(C.d_list)
    d = C.d_list(di);
    for si = 1:numel(zStarts)
        [zOpt, J, flag, localInfo] = fminsearch(@(z) costOf(z, d), zStarts{si}, fmsOpt);
        localEvals = localInfo.funcCount;
        for r = 1:opts.restarts
            [zOpt, J, flag, localInfo] = fminsearch(@(z) costOf(z, d), zOpt, fmsOpt);
            localEvals = localEvals + localInfo.funcCount;
        end
        th = buildTheta(zOpt, d);
        L = score(th);
        trace(end + 1) = struct('d_weeks', d, 'start', si, 'stage', 'nm', ...
                                'J', J, 'nll', L.nll, 'exitflag', flag, ...
                                'funcCount', localEvals, 'iterations', localInfo.iterations); %#ok<AGROW>
        if opts.verbose
            fprintf('    [fit] d=%d start %2d/%2d  J=%.6f  m0=%.4f rho=%.4f\n', ...
                d, si, numel(zStarts), J, th.m0, th.rho);
        end
        if J < best.J
            best = finaliseRecord(th, d, J, L, zOpt);
            best.exitflag = flag;
            best.funcCount = localEvals;
        end
    end
end

% Search-reliability diagnostics.  The distinction matters for interpretation:
% screen_best_by_seed summarises the RANDOM SCREEN, so its spread measures how
% badly uniform sampling performs, NOT how reliable the optimiser is.  What
% bounds optimiser reliability is the spread of the FINAL values reached from
% the different stage-1 starts, and how many independent local runs land on the
% reported optimum.
allStart = [trace.start];
allJ     = [trace.J];
uStart   = unique(allStart);
finalByStart = zeros(numel(uStart), 1);
for j = 1:numel(uStart)
    finalByStart(j) = min(allJ(allStart == uStart(j)));
end

info = struct();
info.free = names;
info.fixed = C.fixed;
info.family = C.family;
info.family_opts = C.family_opts;
info.d_list = C.d_list;
info.n_screen = opts.n_screen;
info.seeds = seedList;
info.n_screen_total = size(Uall, 1);
info.n_starts = numel(zStarts);
info.n_eval = nEval;
info.seed = opts.seed;
info.trace = trace;
info.sub_weeks = C.sub_weeks;
info.phase_penalty = C.phase_penalty;
info.kernel_opts = kOpts;
info.screen_best_J = min(Jscreen);
info.screen_best_by_seed = screenBestBySeed;
info.screen_is_sampled = true;   % stage 1 is a random sample, not exhaustive
info.best_start_J = min(Jscreen(1:nGuaranteed));
info.final_best_J = min(allJ);
info.final_worst_J = max(allJ);
info.final_best_by_start = [uStart(:), finalByStart];
info.n_stage2_runs = numel(allJ);
info.n_stage2_within_1nll = sum(allJ <= min(allJ) + 1.0);
info.n_warm = numel(opts.warm);

    % ---------------------------------------------------------------- nested
    function J = costOf(z, d)
        nEval = nEval + 1;
        th = buildTheta(z, d);
        % One kernel evaluation serves BOTH the likelihood and the phase
        % penalty.  Asking for them separately would double the cost of every
        % optimiser step, and the objective is evaluated ~10^4 times per fit.
        % NOTE: the local is named `kern`, not `k`: a variable that appears
        % inside a nested function becomes SHARED with the parent, and `k` is
        % used as a loop counter in the parent body.
        kern = mwspm_calib_kernel(pre, th, kOpts);
        L = mwspm_calib_nll(kern.Y_obs_sub, kern.Yhat_sub, C.family, C.family_opts);
        J = L.nll;
        if C.phase_penalty
            dev = max(abs(kern.peak_week - C.target_peak), 0);
            J = J + C.pen * dev^2;
        end
        if ~isfinite(J)
            J = 1e12;
        end
    end

    function J = evalU(u, d)
        z = packU(u, lo, hi);
        J = costOf(z, d);
    end

    function th = buildTheta(z, d)
        base = C.fixed;
        base.d_weeks = d;
        th = mwspm_calib_theta('unpack', z, names, C.box, base);
    end

    function L = score(th)
        kern = mwspm_calib_kernel(pre, th, kOpts);
        L = mwspm_calib_nll(kern.Y_obs_sub, kern.Yhat_sub, C.family, C.family_opts);
        L.kernel = kern;
    end
end

% ---------------------------------------------------------------------------
function z = packU(u, lo, hi)
% u (natural units, inside the box) -> z (unbounded simplex coordinates)
w = (u(:)' - lo) ./ (hi - lo);
w = min(max(w, 1e-9), 1 - 1e-9);
z = log(w ./ (1 - w));
z = z(:);
end

% ---------------------------------------------------------------------------
function th = mergeTheta(th, fixed)
fn = fieldnames(fixed);
for k = 1:numel(fn)
    th.(fn{k}) = fixed.(fn{k});
end
end

% ---------------------------------------------------------------------------
function rec = finaliseRecord(th, d, J, L, z)
rec = struct();
rec.theta = th;
rec.d_weeks = d;
rec.J = J;
rec.nll = L.nll;
rec.nll_full = L.nll_full;
rec.pearson = L.pearson;
rec.dispersion = L.dispersion;
rec.df = L.df;
rec.family = L.family;
rec.theta_nb = L.theta;
rec.scale = L.scale;
rec.z = z;
rec.Yhat = L.kernel.Yhat;
rec.Yhat_sub = L.kernel.Yhat_sub;
rec.Y_obs_sub = L.kernel.Y_obs_sub;
rec.dC = L.kernel.dC;
rec.peak_week = L.kernel.peak_week;
rec.true_total = L.kernel.true_total;
rec.reported_total = L.kernel.reported_total;
end
