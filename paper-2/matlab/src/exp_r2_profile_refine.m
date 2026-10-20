function payload = exp_r2_profile_refine(ctx)
%EXP_R2_PROFILE_REFINE  Multi-start refinement at anchors of the R2 profile.
%
%   The R2 sweep used a single continuation start per grid point.  A candidate
%   minimum is only credible if it is reached from starts that do not depend on
%   the sweep order, so this experiment refits a small number of named anchors
%   from three independent starts with a larger evaluation budget:
%     1. the continuation solution recorded by the profile,
%     2. the R1 conditional solution at rho = 0.8,
%     3. the parameter-box centre.
%
%   Anchors per lag: the lowest-nll grid points (default 3) plus the rho values
%   listed in cfg.extra_rhos when they are on the profile grid.
%
%   Output is a local-reproducibility statement, not a confidence set.

cfg = ctx.cfg;
assert(~cfg.phase_penalty, 'R2c:penalty', 'the revised refinement must not use a peak penalty');

if isfield(ctx, 'dryRun') && ctx.dryRun
    payload = struct('outputs', {{}}); payload.dry_run = true; return
end

diary(fullfile(ctx.runDir, 'logs', 'stdout.txt'));
diaryGuard = onCleanup(@() diary('off')); %#ok<NASGU>
snap = fullfile(ctx.dataDir, 'code_snapshot');
mkdir(snap);
copyfile(fullfile(ctx.root, 'src', 'matlab', '*.m'), snap);

D    = mwspm_load_data(ctx.root, ctx.modelCard);
pre  = mwspm_calib_prepare(D, ctx.modelCard, struct('mu_variant', cfg.mu_variant));
pRef = mwspm_params(ctx.modelCard);
box  = mwspm_calib_theta('box');
warmRun = mwspm_readjson(fullfile(ctx.root, 'results', 'runs', cfg.warm_run, 'data', 'stable_fit.json'));

profPath = fullfile(ctx.root, 'results', 'runs', cfg.profile_run, 'data', 'rho_profile_pure.csv');
assert(exist(profPath, 'file') == 2, 'R2c:noProfile', 'missing %s', profPath);
P = mwspm_readcsv(profPath);
lagP = r2c_col(P, 'lag'); rhoP = r2c_col(P, 'rho'); nllP = r2c_col(P, 'nll');

rows = table();
for d = unique(lagP(:))'
    sel = find(lagP == d);
    [~, ord] = sort(nllP(sel));
    pick = sel(ord(1:min(cfg.anchors_per_lag, numel(sel))));
    for r = cfg.extra_rhos(:)'
        j = sel(abs(rhoP(sel) - r) < 1e-9);
        if ~isempty(j), pick = unique([pick; j]); end
    end
    key = sprintf('lag%d', d);
    assert(isfield(warmRun, key), 'R2c:noWarm', 'missing warm fit for %s', key);
    for i = pick(:)'
        rho = rhoP(i);
        th0 = struct('m0', r2c_val(P, 'm0', i), 'rho', rho, 'd_weeks', d, ...
            'onset_doy', r2c_val(P, 'onset_doy', i), 'offset_doy', r2c_val(P, 'offset_doy', i), ...
            'mu_s', r2c_val(P, 'mu_s', i), 'gamma', r2c_val(P, 'gamma', i), ...
            'beta', r2c_val(P, 'beta', i));

        C = struct('free', {cellstr(string(cfg.free))}, 'fixed', struct('rho', rho), ...
            'box', box, 'family', 'poisson', 'phase_penalty', false, ...
            'kernel_opts', cfg.solver, 'd_list', d);
        opt = cfg.optimiser;
        opt.warm = {th0, warmRun.(key).best.theta};
        [best, info] = mwspm_calib_fit(pre, C, opt); %#ok<ASGLU>
        assert(abs(best.J - best.nll) < 1e-10, 'R2c:penaltyLeak', 'objective is not the pure NLL');

        tr = info.trace;
        [tauW, ~, T0at, cumAt] = r2c_tau_star(pre, D, pRef, best.theta, rho, cfg.solver);

        row = table(d, rho, nllP(i), best.nll, nllP(i) - best.nll, ...
            sum([tr.exitflag] == 1), numel(tr), info.n_eval, best.exitflag, ...
            best.theta.m0, best.theta.mu_s, best.theta.onset_doy, best.theta.offset_doy, ...
            best.theta.gamma, best.theta.beta, best.peak_week, best.reported_total, ...
            best.true_total, tauW, T0at, cumAt, ...
            'VariableNames', {'lag','rho','nll_profile','nll_refined','improvement', ...
            'n_starts_converged','n_starts','n_eval','exitflag','m0','mu_s','onset_doy', ...
            'offset_doy','gamma','beta','peak_week','reported_total','true_total', ...
            'tau_star','T0_at_tau','cum_reports_at_tau'});
        rows = [rows; row]; %#ok<AGROW>

        fprintf(['R2c lag%d rho=%.2f nll %.9f -> %.9f (improve %.3g) converged %d/%d ' ...
            'tau*=%s\n'], d, rho, nllP(i), best.nll, nllP(i) - best.nll, ...
            sum([tr.exitflag] == 1), numel(tr), r2c_num(tauW));
    end
end

writetable(rows, fullfile(ctx.dataDir, 'refine_anchors.csv'));
summary = struct();
summary.interpretation = ['local reproducibility of named profile anchors; ' ...
    'no confidence set and no global-optimum claim'];
summary.n_anchors = height(rows);
summary.n_improved = sum(rows.improvement > 1e-6);
summary.all_starts_converged = all(rows.n_starts_converged == rows.n_starts);
summary.max_improvement = max(rows.improvement);
mwspm_writejson(fullfile(ctx.dataDir, 'refine_summary.json'), summary);

payload = summary;
payload.outputs = cfg.outputs;
end

% ---------------------------------------------------------------------------
function v = r2c_col(S, name)
v = NaN(S.n_rows, 1);
idx = find(strcmp(S.names, name), 1);
assert(~isempty(idx), 'R2c:col', 'column "%s" not found', name);
for k = 1:S.n_rows
    c = strtrim(S.rows{k}{idx});
    if ~isempty(c) && ~strcmpi(c, 'na'), v(k) = str2double(c); end
end
end

function v = r2c_val(S, name, row)
c = strtrim(S.rows{row}{find(strcmp(S.names, name), 1)});
if isempty(c) || strcmpi(c, 'na'), v = NaN; else, v = str2double(c); end
end

function [tauStar, r0StartWeek, T0atTau, cumAtTau] = r2c_tau_star(pre, D, pRef, th, rhoUsed, solverOpts)
nW = pre.n_weeks;
sp = struct('onset_doy', th.onset_doy, 'offset_doy', th.offset_doy, ...
    'kappa_days', pre.kappa_days, 'gamma', th.gamma, 'beta', th.beta, ...
    's_floor', pre.s_floor, 'min_window_days', pre.min_window_days, ...
    'normalize', true, 'norm_factor', NaN);
s = mwspm_season(pre.doy, sp);
R0daily = mwspm_r0(D.T_daily, th.m0 * s, th.mu_s * pre.mu_raw_daily, pRef.gamma_h);
wk = floor(((1:numel(R0daily))' - 1) / 7) + 1;
wk = min(max(wk, 1), nW);
gt1 = find(R0daily > 1, 1);
if isempty(gt1), r0StartWeek = NaN; else, r0StartWeek = wk(gt1); end
R0w = accumarray(wk, R0daily, [nW 1], @mean, NaN);
kern = mwspm_calib_kernel(pre, th, solverOpts);
cumRep = cumsum(kern.Yhat(:));
T0w = nan(nW, 1); ok = R0w > 1;
T0w(ok) = rhoUsed ./ log(R0w(ok));
cross = find(isfinite(T0w) & cumRep >= T0w, 1);
if isempty(cross)
    tauStar = NaN; T0atTau = NaN; cumAtTau = NaN;
else
    tauStar = cross; T0atTau = T0w(cross); cumAtTau = cumRep(cross);
end
end

function s = r2c_num(v)
if isfinite(v), s = sprintf('%.6g', v); else, s = 'NA'; end
end
