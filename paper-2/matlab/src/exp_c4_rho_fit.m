function out = exp_c4_rho_fit(ctx)
%EXP_C4_RHO_FIT  T12.01 -- free rho fit + profile likelihood (C4 main report point).
%
%   Experiment entry point for config/experiments/exp_c4_rho_fit.json.
%
%   This is the run that turns "rho is a fitted reporting ratio" into a number
%   WITH an interval.  Unlocked by the G0 ruling S3b (decision
%   protocols/decisions/2026-10-06_G0_novelty_preemption_ruling.md).
%
%   Honest constraints carried in from the T12 audit
%   ------------------------------------------------
%   1. The observation-layer lag d enters the likelihood ONLY through round(d),
%      so the objective is a STAIRCASE in d and its sensitivity column norm is
%      exactly zero.  d must therefore be looped over integers, never optimised
%      as a continuous coordinate (that is precisely the defect that made the
%      legacy profile meaningless: its 11 rows mixed lag 1 and lag 2).
%   2. Lag 1 and lag 2 are statistically indistinguishable at this sample size
%      (they differed by at most 4.23 nll units in the audit).  Both are
%      therefore fitted AND profiled, and both are reported.  Choosing one and
%      hiding the other would be reporting only the flattering branch.
%   3. rho touches the parametrisation bound rho <= 0.80.  That is recorded as
%      `at_bound`, and it bounds the language that may be used later
%      (obligation O1): nothing here may be called "the estimated reporting
%      rate", only "the highest reporting rate the data cannot exclude".
%
%   No synthetic data, no random replicate: this is frequentist profile
%   likelihood on deterministic trajectories.  Intervals come from the profile,
%   never from the spread of optimiser seeds.

logFid = ctx.logFid;
root      = ctx.root;
modelCard = ctx.modelCard;
cfg       = ctx.cfg;

fprintf(logFid, '=== exp_c4_rho_fit (T12.01, C4 main report point) ===\n');
fprintf(logFid, 'claims: %s\n', strjoin(cellstr(string(getdef(cfg, 'claims', {}))), ', '));
fprintf(logFid, 'G0 ruling S3b (2026-10-06) unlocked this experiment.\n');

D   = mwspm_load_data(root, modelCard);
pre = mwspm_calib_prepare(D, modelCard, struct('mu_variant', 'day'));
box   = mwspm_calib_theta('box');
thRef = mwspm_calib_theta('reference');
crit  = mwspm_chisq1_quantile(0.95);

PD = cfg.parameter_domain;
rhoGrid = double(PD.rho_profile_grid.from) : double(PD.rho_profile_grid.step) : double(PD.rho_profile_grid.to);
rhoGrid = round(rhoGrid, 6);
dGrid   = double(PD.d_weeks_integer_grid(:)');
profLags = double(PD.profile_lags(:)');

if isfield(ctx, 'dryRun') && ctx.dryRun
    out = struct();
    out.dry_run = true;
    out.validated = struct('n_weeks', pre.n_weeks, 'observed_total', sum(pre.Y), ...
        'n_rho_grid', numel(rhoGrid), 'd_grid', dGrid, 'profile_lags', profLags);
    fprintf(logFid, '[dry-run] %d rho grid points, d in [%s], profile lags [%s]; no optimisation performed\n', ...
        numel(rhoGrid), mat2str(dGrid), mat2str(profLags));
    return
end

o = cfg.optimiser;
fitOpt = struct('n_screen', double(o.n_screen), 'seeds', double(o.seeds(:)'), ...
                'n_seed', double(o.n_seed), 'restarts', double(o.restarts), ...
                'maxFunEvals', double(o.maxFunEvals), 'maxIter', double(o.maxIter), ...
                'TolX', 1e-8, 'TolFun', 1e-9, 'verbose', false);

freeAll    = {'m0', 'rho', 'onset_doy', 'offset_doy', 'mu_s', 'gamma', 'beta'};
freeStruct = {'m0', 'onset_doy', 'offset_doy', 'mu_s', 'gamma', 'beta'};
baseC = struct('free', {freeAll}, 'fixed', struct(), 'box', box, ...
    'family', 'poisson', 'family_opts', struct(), 'sub_weeks', [], 'y_obs', [], ...
    'phase_penalty', true, 'target_peak', 37.0, 'pen', 300.0);

nWeeks = pre.n_weeks;
kParam = 8;                       % m0, rho, d, onset, offset, mu_s, gamma, beta
warm0  = {thRef};

fprintf(logFid, 'objective: poisson nll + phase soft penalty (pi=%g, target peak week %g)\n', ...
    double(cfg.objective.phase_penalty.pi), double(cfg.objective.phase_penalty.target_peak_week));
fprintf(logFid, 'rho box = [%.4g, %.4g] (source hard constraint)\n', box.rho(1), box.rho(2));

% =============================================================== 1. free fit
fprintf(logFid, '\n--- 1. free fit (rho free, integer lag looped) ---\n');
fits = struct();
tableRows = {};
bestJ = Inf; bestOverall = []; bestD = NaN;

for di = 1:numel(dGrid)
    dd = dGrid(di);
    Cd = baseC;
    Cd.d_list = dd;
    te = tic;
    [b, info] = mwspm_calib_fit(pre, Cd, mergeOpts(fitOpt, struct('warm', {warm0})));
    el = toc(te);

    aic = 2 * b.nll_full + 2 * kParam;
    bic = 2 * b.nll_full + kParam * log(nWeeks);
    pen = b.J - b.nll;
    atBound = (b.theta.rho >= box.rho(2) - 1e-6) || (b.theta.rho <= box.rho(1) + 1e-6);

    key = sprintf('lag%d', dd);
    fits.(key) = struct('d_weeks', dd, 'theta', b.theta, 'J', b.J, 'nll', b.nll, ...
        'nll_full', b.nll_full, 'AIC', aic, 'BIC', bic, ...
        'phase_penalty_contribution', pen, 'rho_at_bound', atBound, ...
        'pearson_dispersion', b.dispersion, 'peak_week', b.peak_week, ...
        'true_total', b.true_total, 'reported_total', b.reported_total, ...
        'n_stage2_within_1nll', info.n_stage2_within_1nll, ...
        'final_best_by_start', info.final_best_by_start, 'n_eval', info.n_eval, ...
        'elapsed_seconds', el);

    fprintf(logFid, '  d=%d : J=%.4f nll=%.4f nll_full=%.4f AIC=%.3f BIC=%.3f | rho=%.6f%s m0=%.6f mu_s=%.4f onset=%.3f offset=%.3f gamma=%.4f beta=%.4f | peak=%d pen=%.6f within1nll=%d (%.1f s)\n', ...
        dd, b.J, b.nll, b.nll_full, aic, bic, b.theta.rho, bndMark(atBound), ...
        b.theta.m0, b.theta.mu_s, b.theta.onset_doy, b.theta.offset_doy, ...
        b.theta.gamma, b.theta.beta, b.peak_week, pen, info.n_stage2_within_1nll, el);

    tableRows{end+1, 1} = {sprintf('free_fit_%s', key), dd, b.theta.m0, b.theta.rho, ...
        b.theta.onset_doy, b.theta.offset_doy, b.theta.mu_s, b.theta.gamma, b.theta.beta, ...
        b.J, b.nll, b.nll_full, aic, bic, b.dispersion, b.peak_week, b.true_total, ...
        b.reported_total, double(atBound)}; %#ok<AGROW>

    if b.J < bestJ
        bestJ = b.J; bestOverall = b; bestD = dd;
    end
end

fprintf(logFid, 'best free fit: d=%d at J=%.4f (nll=%.4f)\n', bestD, bestJ, bestOverall.nll);

% =========================================================== 2. rho profiles
fprintf(logFid, '\n--- 2. profile likelihood of rho (per lag branch) ---\n');
profiles = struct();
profHeader = {'branch_lag', 'grid_value', 'J', 'nll', 'nll_full', 'dnll', ...
              'm0', 'rho', 'onset_doy', 'offset_doy', 'mu_s', 'gamma', 'beta', ...
              'd_weeks', 'peak_week', 'true_total', 'reported_total', 'dispersion', 'n_eval'};
profRowsOut = {};

for di = 1:numel(profLags)
    dd = profLags(di);
    Cg = baseC;
    Cg.free = freeStruct;
    Cg.d_list = dd;
    key = sprintf('lag%d', dd);
    fprintf(logFid, '  [%s] grid %g:%g:%g (%d points)\n', key, rhoGrid(1), ...
        double(PD.rho_profile_grid.step), rhoGrid(end), numel(rhoGrid));
    te = tic;
    P = mwspm_calib_profile(pre, Cg, 'rho', rhoGrid, ...
        mergeOpts(fitOpt, struct('warm', {warm0}, 'anchor', {thRef}, 'crit', crit)));
    el = toc(te);

    % Anchor contract: no grid point may report worse than what the anchor
    % theta itself achieves at that rho.  Measured, not asserted.
    anchorDif = P.nll - anchorNllAt(pre, Cg, thRef, rhoGrid);
    anchorMax = max(anchorDif);
    penMax = max([P.rows.J] - [P.rows.nll]);

    fprintf(logFid, '  [%s] profile min %.4f at rho=%.4f ; Wilks 95%% band [%.4f, %.4f] connected=%d\n', ...
        key, P.nll_min, P.min_at, P.ci_low, P.ci_high, P.ci_is_connected);
    fprintf(logFid, '  [%s] ANCHOR CONTRACT max(profile - anchor) = %+.4f (<=0 required)\n', key, anchorMax);
    fprintf(logFid, '  [%s] phase penalty over grid max(J-nll) = %+.6f (0 = pure likelihood band)\n', key, penMax);
    fprintf(logFid, '  [%s] elapsed %.1f s\n', key, el);

    for i = 1:numel(rhoGrid)
        th = P.rows(i).theta;
        profRowsOut{end+1, 1} = {dd, rhoGrid(i), P.rows(i).J, P.rows(i).nll, P.rows(i).nll_full, ...
            P.dnll(i), th.m0, th.rho, th.onset_doy, th.offset_doy, th.mu_s, th.gamma, th.beta, ...
            P.rows(i).d_weeks, P.rows(i).peak_week, P.rows(i).true_total, P.rows(i).reported_total, ...
            P.rows(i).dispersion, P.rows(i).n_eval}; %#ok<AGROW>
    end

    profiles.(key) = struct('d_weeks', dd, 'grid', rhoGrid(:)', 'J', P.J(:)', 'nll', P.nll(:)', ...
        'dnll', P.dnll(:)', 'nll_min', P.nll_min, 'min_at', P.min_at, ...
        'ci_low', P.ci_low, 'ci_high', P.ci_high, 'ci_is_connected', P.ci_is_connected, ...
        'crit', P.crit, 'anchor_contract_max', anchorMax, 'phase_penalty_max', penMax, ...
        'at_BIC_best_rho', struct(), 'elapsed_seconds', el, 'n_eval_total', sum([P.rows.n_eval]));
end

% ======================================================= 3. main report point
% EXPERIMENT_PLAN Sec.6: the main report point is the BIC-optimum.  BIC uses
% nll_full (the density normalisation), NOT J (which carries the phase penalty
% and is therefore not a log-likelihood).
allBIC = zeros(1, numel(dGrid));
dKeys  = cell(1, numel(dGrid));
for di = 1:numel(dGrid)
    key = sprintf('lag%d', dGrid(di));
    dKeys{di} = key;
    allBIC(di) = fits.(key).BIC;
end
[minBIC, iBIC] = min(allBIC);
mainKey = dKeys{iBIC};
mainLag = dGrid(iBIC);
main = fits.(mainKey);

% The source parametrisation cannot express rho > 0.80, so the reported point
% must additionally be admissible there.
admissOK = main.theta.rho <= box.rho(2) + 1e-9;

fprintf(logFid, '\n--- 3. main report point (BIC optimum) ---\n');
fprintf(logFid, 'BIC over the free fits: %s\n', mat2str(allBIC, 4));
fprintf(logFid, 'main: %s  rho=%.6f  BIC=%.4f  admissible under rho<=%.2f: %d\n', ...
    mainKey, main.theta.rho, minBIC, box.rho(2), double(admissOK));

% =================================================================== outputs
bestJson = struct();
bestJson.provenance = struct('parent_run_id', ctx.runId, ...
    'experiment_id', ctx.expId, 'config_sha256', ctx.cfgHash, ...
    'note', 'exp_c4_gate_g1 (T13.01) must consume THIS run; it exits with an explicit error otherwise.');
bestJson.selection_rule = 'minimum BIC over the integer-lag free fits; BIC computed from nll_full, not from J';
bestJson.free_fits = fits;
bestJson.main_report_point = struct('d_weeks', mainLag, 'theta', main.theta, 'nll', main.nll, ...
    'nll_full', main.nll_full, 'AIC', main.AIC, 'BIC', main.BIC, ...
    'rho_at_bound', main.rho_at_bound, 'admissible_under_source_box', admissOK, ...
    'phase_penalty_contribution', main.phase_penalty_contribution, ...
    'peak_week', main.peak_week, 'true_total', main.true_total, ...
    'reported_total', main.reported_total);
bestJson.profiles = profiles;
bestJson.k_parameters = kParam;
bestJson.n_weeks = nWeeks;
bestJson.criterion = crit;
bestJson.legacy_reference_theta = thRef;
bestJson.caution = ['Both lag branches are reported because they are statistically ' ...
    'indistinguishable (audit: max 4.23 nll apart). Selecting one and suppressing the ' ...
    'other would be reporting only the convenient branch.'];
mwspm_writejson(fullfile(ctx.dataDir, 'fit_best.json'), bestJson);

% profile table
fid = fopen(fullfile(ctx.dataDir, 'rho_profile.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '%s\n', strjoin(profHeader, ','));
for i = 1:numel(profRowsOut)
    rw = profRowsOut{i};
    fprintf(fid, '%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%d,%d,%.17g,%.17g,%.17g,%d\n', rw{:});
end
fclose(fid);

% weekly fit at the main report point
kernMain = mwspm_calib_kernel(pre, main.theta, struct());
mwspm_writecsv(fullfile(ctx.dataDir, 'weekly_fit.csv'), ...
    {'week_index', 'observed_local', 'model_reported', 'model_true_incidence', 'model_true_incidence_lagged'}, ...
    [(1:nWeeks)', pre.Y(:), kernMain.Yhat(:), kernMain.dC(:), kernMain.dC_delayed(:)]);

% parameter table
fid = fopen(fullfile(ctx.dataDir, 'parameter_table.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, 'point_id,d_weeks,m0,rho,onset_doy,offset_doy,mu_s,gamma,beta,J,nll,nll_full,AIC,BIC,dispersion,peak_week,true_total,reported_total,rho_at_bound\n');
for i = 1:numel(tableRows)
    rw = tableRows{i};
    fprintf(fid, '%s,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%d,%.17g,%.17g,%d\n', rw{:});
end
fclose(fid);

summary = struct();
summary.experiment_id = ctx.expId;
summary.task = char(string(getdef(cfg, 'task', '-')));
summary.gate_satisfied = 'G0 = S3b (2026-10-06)';
summary.main_report_point = bestJson.main_report_point;
summary.n_free_fits = numel(dGrid);
summary.profile_lags = profLags;
summary.n_rho_grid = numel(rhoGrid);
summary.wilks95 = struct('lag1', [], 'lag2', []);
for k = 1:numel(profLags)
    key = sprintf('lag%d', profLags(k));
    if isfield(profiles, key)
        summary.wilks95.(key) = struct('low', profiles.(key).ci_low, ...
            'high', profiles.(key).ci_high, 'connected', profiles.(key).ci_is_connected, ...
            'min_at', profiles.(key).min_at);
    end
end
summary.what_intervals_mean = ['Wilks bands come from the profile likelihood of a DETERMINISTIC ' ...
    'trajectory; they quantify parameter uncertainty, NOT stochastic replicate scatter. ' ...
    'Time points of a single deterministic trajectory are NOT independent replicates and no ' ...
    'CI is computed from them.'];
mwspm_writejson(fullfile(ctx.dataDir, 'summary.json'), summary);

out = struct();
out.status = 'SUCCESS';
out.main_report_point = bestJson.main_report_point;
out.profiles = profiles;
out.outputs = {'fit_best.json', 'rho_profile.csv', 'weekly_fit.csv', 'parameter_table.csv', 'summary.json'};
out.deterministic = true;
out.parent_run_id = ctx.runId;
out.rng = 'fixed optimiser seeds for reproducible starts; the spread of start values is NOT an interval';
out.note = ['T12.01 main report point. rho touches the source bound, so per obligation O1 it may ' ...
    'only be described as the highest reporting rate the data cannot exclude.'];
end

% ---------------------------------------------------------------------------
function nllAt = anchorNllAt(pre, C, thRef, grid)
% Objective that the reference theta itself achieves when rho is forced to each
% grid value.  Comparing the profile against it is a measurable contract.
nllAt = nan(numel(grid), 1);
for i = 1:numel(grid)
    th = thRef;
    th.rho = grid(i);
    th.d_weeks = C.d_list(1);
    k = mwspm_calib_kernel(pre, th, struct());
    L = mwspm_calib_nll(pre.Y, k.Yhat, 'poisson', struct());
    nllAt(i) = L.nll;
end
end

% ---------------------------------------------------------------------------
function o = mergeOpts(base, over)
o = base;
fn = fieldnames(over);
for k = 1:numel(fn)
    o.(fn{k}) = over.(fn{k});
end
end

% ---------------------------------------------------------------------------
function s = bndMark(tf)
if tf, s = '*'; else s = ' '; end
end

% ---------------------------------------------------------------------------
function v = getdef(s, f, dflt)
if isfield(s, f) && ~isempty(s.(f))
    v = s.(f);
else
    v = dflt;
end
end
