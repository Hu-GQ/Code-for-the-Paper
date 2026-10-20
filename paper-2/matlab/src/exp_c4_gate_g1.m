function out = exp_c4_gate_g1(ctx)
%EXP_C4_GATE_G1  T13.01 -- gate G1: is the C4 "alarm week" stable in rho?
%
%   Experiment entry point for config/experiments/exp_c4_gate_g1.json.
%
%   The gate compares four weeks on the fitted series:
%     (i)   current-caliber alarm week  = first week whose CUMULATIVE OBSERVED
%           reported count reaches the threshold (100 by default);
%     (ii)  tau*                        = min{ w : cum_w >= T0(w) }, where
%           cum_w is the cumulative DETECTED (model) count and
%           T0(w) = rho / log R0(w)   (Hartfield & Alizon 2013, with the
%           2017 Correction, rho = the detection fraction p);
%     (iii) R0 > 1 window start week;
%     (iv)  observation start week      = first week with a non-zero observed
%           local case count.
%   Verdict: GO iff the range of (ii) over rho in [0.30, 1.00] is <= 4 weeks.
%   Failing the gate triggers stop condition S1: C4 drops to a "threshold
%   sensitivity" section and T04 is reopened.  The verdict is not negotiable
%   by wording.
%
%   STRUCTURAL NOTE THAT SHAPES THE READING (theorems.tex rem:rho-decouple)
%   --------------------------------------------------------------------
%   rho appears in cum_w (as a multiplicative factor) AND in the numerator of
%   T0.  For a FIXED trajectory the crossing condition
%       rho * C_lag(w) >= rho / log R0(w)
%   is therefore EXACTLY INDEPENDENT of rho: it reduces to
%       C_lag(w) >= 1 / log R0(w).
%   So any spread of tau* over the rho grid does NOT come from rho directly; it
%   comes entirely from the fact that, at each rho, the other parameters are
%   RE-OPTIMISED, which changes the trajectory and hence R0(w).  That is what
%   the gate measures, and the report must say so (CL-22).
%
%   This experiment is POST-PROCESSING of a parent run: it does not fit
%   anything and does not re-integrate the ODE.

logFid = ctx.logFid;
root   = ctx.root;
cfg    = ctx.cfg;
modelCard = ctx.modelCard;

fprintf(logFid, '=== exp_c4_gate_g1 (T13.01, gate G1) ===\n');
fprintf(logFid, 'claims: %s\n', strjoin(cellstr(string(getdef(cfg, 'claims', {}))), ', '));

D   = mwspm_load_data(root, modelCard);
pre = mwspm_calib_prepare(D, modelCard, struct('mu_variant', 'day'));
pRef = mwspm_params(modelCard);

% ------------------------------------------------ locate the parent run ------
cfgParent = getdef(cfg, 'parent_run_id', []);
parent = gate_find_parent(root, cfgParent);
assert(~isempty(parent), 'gate:noParent', ...
    ['no run of exp_c4_rho_fit was found under results/runs/. The gate must consume a ' ...
     'real run; it may not fall back to placeholder numbers.']);
fprintf(logFid, 'parent run       : %s\n', parent.runId);
fprintf(logFid, 'parent fit_best  : %s\n', parent.fitPath);
fprintf(logFid, 'parent profile   : %s\n', parent.profilePath);

fit = mwspm_readjson(parent.fitPath);
assert(isfield(fit, 'provenance') && isfield(fit.provenance, 'parent_run_id'), ...
    'gate:badParent', 'parent fit_best.json carries no provenance.parent_run_id');
fprintf(logFid, 'parent provenance: %s\n', char(string(fit.provenance.parent_run_id)));

prof = gate_read_profile(parent.profilePath);

if isfield(ctx, 'dryRun') && ctx.dryRun
    out = struct();
    out.dry_run = true;
    out.parent_run_id = parent.runId;
    out.validated = struct('n_profile_rows', numel(prof.branch_lag), ...
        'branches', unique(prof.branch_lag(:)'));
    fprintf(logFid, '[dry-run] parent and profile readable; no gate computation performed\n');
    return
end

ACC = cfg.acceptance_criterion;
gridDef = cfg.parameter_domain.rho_gate_grid;
rhoGate = double(gridDef.from) : double(gridDef.step) : double(gridDef.to);
rhoGate = round(rhoGate, 6);
rangeMax = double(ACC.range_weeks.threshold_weeks);

fprintf(logFid, 'rho gate grid    : %g : %g : %g  (%d points)\n', ...
    rhoGate(1), double(gridDef.step), rhoGate(end), numel(rhoGate));
fprintf(logFid, 'acceptance       : range of tau* over the grid <= %g weeks\n', rangeMax);

% --------------------------------------------------------- (i) and (iv) -----
Yobs = pre.Y(:);
cumObs = cumsum(Yobs);
thresholds = double(cfg.parameter_domain.threshold_drift_cases(:)');
n = numel(Yobs);

ivObsStart = find(Yobs > 0, 1);
if isempty(ivObsStart), ivObsStart = NaN; end

driftRows = {};
for k = 1:numel(thresholds)
    th = thresholds(k);
    wk = find(cumObs >= th, 1);
    if isempty(wk), wk = NaN; end
    driftRows{end+1, 1} = {th, wk}; %#ok<AGROW>
    fprintf(logFid, '  current caliber: cumulative observed >= %-6g first at week %s\n', ...
        th, num2strOrNA(wk));
end
iPeak100 = find(cumObs >= 100, 1);
if isempty(iPeak100), iPeak100 = NaN; end
fprintf(logFid, '(i) current caliber alarm week (100 cases) = %s\n', num2strOrNA(iPeak100));
fprintf(logFid, '(iv) observation start week = %s\n', num2strOrNA(ivObsStart));

fid = fopen(fullfile(ctx.dataDir, 'threshold_drift.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, 'threshold_cumulative_cases,first_week_reaching\n');
for k = 1:numel(driftRows)
    rw = driftRows{k};
    fprintf(fid, '%g,%d\n', rw{1}, rw{2});
end
fclose(fid);

% ------------------------------------------------ per-branch gate sweep -----
branches = unique(prof.branch_lag(:)');
weeksRows = {};
verdict = struct();
verdict.branches = struct();
allRanges = [];

for bi = 1:numel(branches)
    dd = branches(bi);
    key = sprintf('lag%d', dd);
    sel = prof.branch_lag(:) == dd;
    gvals = prof.grid_value(sel);
    fprintf(logFid, '\n--- branch %s : tau* over the rho gate grid ---\n', key);

    tau = nan(numel(rhoGate), 1);
    r0start = nan(numel(rhoGate), 1);
    T0at = nan(numel(rhoGate), 1);
    cumAt = nan(numel(rhoGate), 1);

    for i = 1:numel(rhoGate)
        rg = rhoGate(i);
        [~, j] = min(abs(gvals - rg));
        if abs(gvals(j) - rg) > 1e-6
            % The gate grid must be a subset of the profile grid; if it is not,
            % the parent run was configured differently and the gate refuses to
            % interpolate a refit that was never computed.
            error('gate:gridMismatch', ...
                'rho=%g is not on the parent profile grid (closest %g). Refusing to interpolate.', ...
                rg, gvals(j));
        end
        rowIdx = find(sel);
        k = rowIdx(j);
        th = struct();
        th.m0         = prof.m0(k);
        th.rho        = prof.rho(k);
        th.onset_doy  = prof.onset_doy(k);
        th.offset_doy = prof.offset_doy(k);
        th.mu_s       = prof.mu_s(k);
        th.gamma      = prof.gamma(k);
        th.beta       = prof.beta(k);
        th.d_weeks    = prof.d_weeks(k);

        [tauW, r0StartW, T0w, cumw] = gate_tau_star(pre, D, pRef, th, prof.grid_value(k));
        tau(i) = tauW;
        r0start(i) = r0StartW;
        T0at(i) = T0w;
        cumAt(i) = cumw;

        weeksRows{end+1, 1} = {dd, prof.grid_value(k), tauW, T0w, cumw, r0StartW, ...
            ivObsStart, iPeak100}; %#ok<AGROW>
    end

    fin = isfinite(tau);
    if ~any(fin)
        v = 'NO_GO';
        rng = NaN;
        tauMin = NaN; tauMax = NaN;
    else
        tauMin = min(tau(fin));
        tauMax = max(tau(fin));
        rng = tauMax - tauMin;
        v = tern(rng <= rangeMax + 1e-12, 'GO', 'NO_GO');
    end
    allRanges(end+1) = rng; %#ok<AGROW>

    verdict.branches.(key) = struct('d_weeks', dd, 'tau_star_min', tauMin, ...
        'tau_star_max', tauMax, 'tau_star_range_weeks', rng, ...
        'verdict', v, 'threshold_weeks', rangeMax, ...
        'n_rho_points', numel(rhoGate), 'n_undefined', sum(~fin), ...
        'tau_star_by_rho', tau(:)', 'rho_grid', rhoGate(:)', ...
        'R0_gt1_start_week_min', min(r0start(isfinite(r0start))), ...
        'R0_gt1_start_week_max', max(r0start(isfinite(r0start))), ...
        'T0_at_tau_star_min', min(T0at(fin)), 'T0_at_tau_star_max', max(T0at(fin)));

    fprintf(logFid, '  tau* in [%s, %s] ; range = %s weeks ; verdict = %s\n', ...
        num2strOrNA(tauMin), num2strOrNA(tauMax), num2strOrNA(rng), v);
    fprintf(logFid, '  R0>1 window start week over the grid: [%s, %s]\n', ...
        num2strOrNA(verdict.branches.(key).R0_gt1_start_week_min), ...
        num2strOrNA(verdict.branches.(key).R0_gt1_start_week_max));
    fprintf(logFid, '  |tau* - (iv)| = %s weeks (reported, no hard cut-off in the plan)\n', ...
        num2strOrNA(tern(isfinite(tauMin) && isfinite(ivObsStart), abs(tauMin - ivObsStart), NaN)));
end

% ------------------------------------------------------------ overall verdict
% The gate is passed only if EVERY profiled branch passes: choosing the branch
% that happens to pass would be selecting the convenient answer.
vAll = 'GO';
for bi = 1:numel(branches)
    key = sprintf('lag%d', branches(bi));
    if strcmp(verdict.branches.(key).verdict, 'NO_GO')
        vAll = 'NO_GO';
    end
end
verdict.verdict = vAll;
verdict.criterion = sprintf('range of tau* over rho in [%g, %g] <= %g weeks, for EVERY profiled lag branch', ...
    rhoGate(1), rhoGate(end), rangeMax);
verdict.structural_note = ['rho cancels EXACTLY between cum_w and T0(w) for a fixed ' ...
    'trajectory (theorems.tex rem:rho-decouple). The observed spread of tau* is therefore ' ...
    'produced by the per-rho RE-OPTIMISATION of the other parameters, not by rho itself.'];
verdict.consequence_if_no_go = 'stop condition S1: C4 is downgraded to a threshold-sensitivity section and T04 is reopened.';
fprintf(logFid, '\n=== GATE G1 VERDICT: %s ===\n', vAll);
mwspm_writejson(fullfile(ctx.dataDir, 'gate_g1_verdict.json'), verdict);

fid = fopen(fullfile(ctx.dataDir, 'gate_g1_weeks.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, 'branch_lag,rho,tau_star_week,T0_at_tau_star,cum_detected_at_tau_star,R0_gt1_start_week,observation_start_week,current_caliber_100_week\n');
for i = 1:numel(weeksRows)
    rw = weeksRows{i};
    fprintf(fid, '%d,%.17g,%s,%s,%s,%s,%s,%s\n', rw{1}, rw{2}, numOrNA(rw{3}), numOrNA(rw{4}), ...
        numOrNA(rw{5}), numOrNA(rw{6}), numOrNA(rw{7}), numOrNA(rw{8}));
end
fclose(fid);

summary = struct();
summary.experiment_id = ctx.expId;
summary.task = char(string(getdef(cfg, 'task', '-')));
summary.parent_run_id = parent.runId;
summary.parent_provenance = fit.provenance;
summary.verdict = vAll;
summary.branches = verdict.branches;
summary.current_caliber = struct('week_at_100_cumulative_cases', iPeak100, ...
    'observation_start_week', ivObsStart, 'threshold_drift', cfg.parameter_domain.threshold_drift_cases);
summary.main_report_point = fit.main_report_point;
summary.honesty = ['The verdict is GO / NO_GO and is not softened by wording. If NO_GO, ' ...
    'S1 applies: C4 is downgraded and T04 is reopened.'];
mwspm_writejson(fullfile(ctx.dataDir, 'summary.json'), summary);

out = struct();
out.status = 'SUCCESS';
out.verdict = vAll;
out.parent_run_id = parent.runId;
out.ranges = allRanges;
out.tau_star_by_branch = verdict.branches;
out.outputs = {'gate_g1_weeks.csv', 'gate_g1_verdict.json', 'threshold_drift.csv', 'summary.json'};
out.deterministic = true;
out.rng = 'none (post-processing of a deterministic parent run)';
out.note = ['Gate G1. A NO_GO triggers S1: C4 becomes a threshold-sensitivity section and ' ...
    'T04 is reopened; this run does not attempt to rescue the claim by rewording.'];
end

% ---------------------------------------------------------------------------
function [tauStar, r0StartWeek, T0atTau, cumAtTau] = gate_tau_star(pre, D, pRef, th, rhoUsed)
% Weekly R0 and cumulative DETECTED incidence for one refit theta, then tau*.
nW = pre.n_weeks;

sp = struct('onset_doy', th.onset_doy, 'offset_doy', th.offset_doy, ...
    'kappa_days', pre.kappa_days, 'gamma', th.gamma, 'beta', th.beta, ...
    's_floor', pre.s_floor, 'min_window_days', pre.min_window_days, ...
    'normalize', true, 'norm_factor', NaN);
s = mwspm_season(pre.doy, sp);
trRaw = mwspm_thermal_traits(D.T_daily, 'raw');
muEff = th.mu_s * trRaw.mu;
R0daily = mwspm_r0(D.T_daily, th.m0 * s, muEff, pRef.gamma_h);

% Day -> week index (1..47) with the same convention as mwspm_r0_timeseries.
wk = floor(((1:numel(R0daily))' - 1) / 7) + 1;
wk = min(max(wk, 1), nW);

% (iii): first week containing a day with R0 > 1
gt1 = find(R0daily > 1, 1);
if isempty(gt1)
    r0StartWeek = NaN;
else
    r0StartWeek = wk(gt1);
end

% Weekly R0 = MEAN of the daily values inside the week.  Declared explicitly:
% the daily R0(t) is sawtoothed (only 96 of the 126 days inside the annually
% favourable window have R0 > 1), so any weekly aggregation is a choice.
R0w = accumarray(wk, R0daily, [nW 1], @mean, NaN);

kern = mwspm_calib_kernel(pre, th, struct());
cumDet = cumsum(kern.Yhat(:));           % cumulative DETECTED (model) cases

T0w = nan(nW, 1);
ok = R0w > 1;
T0w(ok) = rhoUsed ./ log(R0w(ok));

cross = find(isfinite(T0w) & cumDet >= T0w, 1);
if isempty(cross)
    tauStar = NaN; T0atTau = NaN; cumAtTau = NaN;
else
    tauStar = cross; T0atTau = T0w(cross); cumAtTau = cumDet(cross);
end
end

% ---------------------------------------------------------------------------
function P = gate_find_parent(root, wantId)
% Latest results/runs/<run_id>/data/fit_best.json produced by exp_c4_rho_fit.
runsDir = fullfile(root, 'results', 'runs');
d = dir(runsDir);
best = [];
for k = 1:numel(d)
    if ~d(k).isdir || strcmp(d(k).name, '.') || strcmp(d(k).name, '..')
        continue
    end
    f = fullfile(runsDir, d(k).name, 'data', 'fit_best.json');
    if exist(f, 'file') ~= 2
        continue
    end
    if ~isempty(wantId) && ischar(wantId) && ~strcmp(char(wantId), d(k).name)
        continue
    end
    try
        j = mwspm_readjson(f);
    catch
        continue
    end
    if ~isfield(j, 'provenance') || ~isfield(j.provenance, 'experiment_id')
        continue
    end
    if ~strcmp(char(string(j.provenance.experiment_id)), 'exp_c4_rho_fit')
        continue
    end
    % also require the run itself to be SUCCESS
    manF = fullfile(runsDir, d(k).name, 'manifest.json');
    if exist(manF, 'file') == 2
        man = mwspm_readjson(manF);
        if isfield(man, 'status') && ~strcmp(char(string(man.status)), 'SUCCESS')
            continue
        end
    end
    if isempty(best) || strcmp(d(k).name, best.runId) > 0
        best = struct('runId', d(k).name, 'fitPath', f, ...
            'profilePath', fullfile(runsDir, d(k).name, 'data', 'rho_profile.csv'));
    end
end
P = best;
end

% ---------------------------------------------------------------------------
function prof = gate_read_profile(path)
T = mwspm_readcsv(path);
prof = struct();
cols = {'branch_lag', 'grid_value', 'J', 'nll', 'nll_full', 'dnll', 'm0', 'rho', ...
        'onset_doy', 'offset_doy', 'mu_s', 'gamma', 'beta', 'd_weeks', ...
        'peak_week', 'true_total', 'reported_total', 'dispersion', 'n_eval'};
for k = 1:numel(cols)
    prof.(cols{k}) = mwspm_csvcol(T, cols{k}, 'numeric');
end
end

% ---------------------------------------------------------------------------
function s = num2strOrNA(v)
if isfinite(v), s = sprintf('%g', v); else s = 'NA'; end
end
function s = numOrNA(v)
if isfinite(v), s = sprintf('%.17g', v); else s = 'NA'; end
end
function v = tern(c, a, b)
if c, v = a; else v = b; end
end
function v = getdef(s, f, dflt)
if isfield(s, f) && ~isempty(s.(f))
    v = s.(f);
else
    v = dflt;
end
end
