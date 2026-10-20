function out = exp_e0_smoke_baseline(ctx)
%EXP_E0_SMOKE_BASELINE  T10 minimal run check (pipeline self-test).
%
%   Experiment entry point for config/experiments/exp_e0_smoke_baseline.json.
%
%   This is NOT a claim-bearing experiment.  It exercises the whole chain
%     config -> data loader -> forcing -> RHS -> stiffness/solver choice ->
%     run directory + manifest + hashes -> observation mapping -> metrics
%   with two solvers on the same model, and compares the result with the
%   legacy Python refit values recorded in the model card.  A mismatch is a
%   FINDING to be reported, not something to be tuned away.
%
%   Separation of concerns: this file owns experiment logic only; the model
%   equations live in mwspm_rhs.m and the integration machinery in
%   mwspm_solve.m / mwspm_integrate_fixed.m.

logFid = ctx.logFid;
fprintf(logFid, '=== exp_e0_smoke_baseline ===\n');

root      = ctx.root;
modelCard = ctx.modelCard;
cfg       = ctx.cfg;

p = mwspm_params(modelCard);
fprintf(logFid, 'model %s / %s / spec %s\n', p.model_id, p.model_version, p.spec_version);
fprintf(logFid, 'm0=%.6g  mu_s=%.6g  rho=%.4g  d=%d  window %s .. %s (%d days / %d weeks)\n', ...
    p.m0, p.mu_s, p.rho, round(p.d_weeks), p.t0_date, p.t_end_date, p.n_days, p.n_weeks);

D = mwspm_load_data(root, modelCard);
fprintf(logFid, 'data: T in [%.2f, %.2f] degC, import total %.0f, observed local total %.0f, obs peak week %d\n', ...
    D.diagnostics.temperature_min_c, D.diagnostics.temperature_max_c, ...
    D.diagnostics.import_total, D.diagnostics.observed_local_total, ...
    D.diagnostics.observed_peak_week);

reportGrid = (0:p.n_days)';

% ---- forcing for both mortality variants --------------------------------
Fraw = mwspm_make_forcing(D, p, struct('interp', 'zoh', 'mu_variant', 'raw'));
Fday = mwspm_make_forcing(D, p, struct('interp', 'zoh', 'mu_variant', 'day'));
stiffRaw = mwspm_stiffness(Fraw, p, 1.0);
stiffDay = mwspm_stiffness(Fday, p, 1.0);

fprintf(logFid, 'stiffness: ratio=%.4g (fast %.4g / slow %.4g per day), stable dt limit %.4g d\n', ...
    stiffRaw.stiffness_ratio, stiffRaw.rate_fast_per_day, stiffRaw.rate_slow_per_day, ...
    stiffRaw.dt_euler_stable_max);
fprintf(logFid, 'euler stability @dt=1: raw max|1-mu_eff*dt|=%.4f (%d/329 days outside), day variant %.4f (%d/329)\n', ...
    stiffRaw.max_abs_1_minus_mu_eff_dt_raw, stiffRaw.n_days_outside_stability_raw, ...
    stiffDay.max_abs_1_minus_mu_eff_dt_day, stiffDay.n_days_outside_stability_day);

if isfield(ctx, 'dryRun') && ctx.dryRun
    out = struct();
    out.dry_run = true;
    out.validated = struct( ...
        'n_days', D.n_days, 'n_weeks', D.n_weeks, ...
        'observed_local_total', D.diagnostics.observed_local_total, ...
        'stiffness_ratio', stiffRaw.stiffness_ratio, ...
        'n_solvers', numel(cfg.solvers));
    fprintf(logFid, '[dry-run] config/data validated, no integration performed\n');
    return
end

% ---- run every configured solver ----------------------------------------
nSol = numel(cfg.solvers);
runs = cell(nSol, 1);
for k = 1:nSol
    spec = cfg.solvers{k};
    if isfield(spec, 'mu_variant') && strcmp(spec.mu_variant, 'day')
        F = Fday;
        stiff = stiffDay;
    else
        F = Fraw;
        stiff = stiffRaw;
    end
    rhs = @(t, x) mwspm_rhs(t, x, p, F);
    x0  = mwspm_initial_state(D, p, F);

    fprintf(logFid, '\n[solver %d/%d] id=%s kind=%s method=%s\n', k, nSol, spec.id, spec.kind, spec.method);
    S = mwspm_solve(struct('rhs', rhs, 't0', 0, 'tf', p.n_days, 'x0', x0, ...
        'spec', spec, 'reportGrid', reportGrid, 'events', cfg.termination.events, ...
        'stiff', stiff, 'pickCfg', getfielddef(cfg, 'auto_solver', struct()), 'logFid', logFid));

    O = mwspm_observe(S.t, S.x(:, 8), p, D);

    fprintf(logFid, '  method=%s n_steps=%d step in [%.4g, %.4g] d  terminated=%s\n', ...
        S.method, S.n_steps, S.min_internal_step, S.max_internal_step, S.termination_reason);
    fprintf(logFid, '  R2_log=%.6f R2_raw=%.6f RMSE_log=%.4f peak week=%d (obs %d) reported_total=%.2f true_total=%.2f\n', ...
        O.metrics.R2_log, O.metrics.R2_raw, O.metrics.RMSE_log, ...
        O.metrics.peak_week_model, O.metrics.peak_week_observed, ...
        O.metrics.reported_total_model, O.metrics.true_total);

    r = struct();
    r.id = spec.id;
    r.solver = struct('kind', S.kind, 'method', S.method, 'mu_variant', S.mu_variant, ...
        'relTol', getfielddef(spec, 'RelTol', NaN), 'absTol', getfielddef(spec, 'AbsTol', NaN), ...
        'dt_days', getfielddef(spec, 'dt_days', NaN), ...
        'clamp_nonnegative', getfielddef(spec, 'clamp_nonnegative', false));
    r.n_steps = S.n_steps;
    r.min_internal_step = S.min_internal_step;
    r.max_internal_step = S.max_internal_step;
    r.termination_reason = S.termination_reason;
    r.solver_realisation = S.solver_realisation;
    r.selection = S.selection;
    r.metrics = O.metrics;
    r.O = O;
    r.S = S;
    runs{k} = r;
end

% ---- R0(t) ---------------------------------------------------------------
R0ts = mwspm_r0_timeseries(D, p, struct('caliber', 'Nv', 'mu_variant', 'raw'));
fprintf(logFid, '\nR0(t): max %.4f at week %d ; R0>1 in weeks %d..%d (%d days)\n', ...
    R0ts.max_R0, R0ts.peak_week, R0ts.gt1_first_week, R0ts.gt1_last_week, R0ts.n_days_above1);

% ---- outputs -------------------------------------------------------------
daily = [(0:p.n_days)', runs{1}.S.x];
mwspm_writecsv(fullfile(ctx.dataDir, 'daily_timeseries.csv'), ...
    {'t_day', 'S_h', 'E_h', 'I_h', 'R_h', 'S_v', 'E_v', 'I_v', 'Cum'}, daily);

nW = p.n_weeks;
weekly = [(1:nW)', D.Y_weekly(:), runs{1}.O.Yhat(:), runs{1}.O.dC(:), runs{2}.O.Yhat(:)];
mwspm_writecsv(fullfile(ctx.dataDir, 'weekly_fit.csv'), ...
    {'week_index', 'observed_local', 'model_reported_ode15s_raw', 'model_true_incidence', 'model_reported_euler_dt1_day'}, ...
    weekly);

Td = D.T_daily; trR = mwspm_thermal_traits(Td, 'raw');
r0daily = [(1:D.n_days)', Td, R0ts.m(:), p.mu_s * trR.mu, R0ts.R0(:)];
mwspm_writecsv(fullfile(ctx.dataDir, 'r0_daily.csv'), ...
    {'day', 'temperature_c', 'm_Nv_caliber', 'mu_eff_raw', 'R0'}, r0daily);

solverDiag = struct();
solverDiag.reference_metrics = runs{1}.metrics;
solverDiag.comparison = struct();
for k = 1:numel(runs)
    solverDiag.comparison.(matlab.lang.makeValidName(char(runs{k}.id))) = struct( ...
        'method', runs{k}.solver.method, 'kind', runs{k}.solver.kind, ...
        'mu_variant', runs{k}.solver.mu_variant, ...
        'n_steps', runs{k}.n_steps, ...
        'min_internal_step', runs{k}.min_internal_step, ...
        'max_internal_step', runs{k}.max_internal_step, ...
        'termination_reason', runs{k}.termination_reason, ...
        'solver_realisation', runs{k}.solver_realisation, ...
        'selection_rule', runs{k}.selection.rule, ...
        'R2_log', runs{k}.metrics.R2_log, ...
        'R2_raw', runs{k}.metrics.R2_raw, ...
        'RMSE_log', runs{k}.metrics.RMSE_log, ...
        'peak_week_model', runs{k}.metrics.peak_week_model, ...
        'reported_total_model', runs{k}.metrics.reported_total_model, ...
        'true_total_model', runs{k}.metrics.true_total, ...
        'reported_total_observed', runs{k}.metrics.reported_total_observed);
end
solverDiag.stiffness_raw = stiffRaw;
solverDiag.stiffness_day = stiffDay;
mwspm_writejson(fullfile(ctx.dataDir, 'solver_diagnostics.json'), solverDiag);

% ---- summary + consistency check against the model card ------------------
ref = modelCard.reference_metrics;
summary = struct();
summary.experiment_id = ctx.expId;
summary.model_version = p.model_version;
summary.n_days = D.n_days;
summary.n_weeks = D.n_weeks;
summary.n_rhs_state = 8;
summary.report_grid = sprintf('daily integers 0..%d', p.n_days);
summary.internal_step = 'adaptive for ode15s (solver-chosen), dt=1 day for the source-faithful Euler row';
summary.stiffness = struct('ratio', stiffRaw.stiffness_ratio, ...
    'dt_euler_stable_max', stiffRaw.dt_euler_stable_max, ...
    'max_abs_1_minus_mu_eff_dt_raw', stiffRaw.max_abs_1_minus_mu_eff_dt_raw, ...
    'n_days_outside_stability_raw', stiffRaw.n_days_outside_stability_raw, ...
    'n_days_outside_stability_day', stiffDay.n_days_outside_stability_day);
summary.R0 = struct('max', R0ts.max_R0, 'peak_week', R0ts.peak_week, ...
    'gt1_first_week', R0ts.gt1_first_week, 'gt1_last_week', R0ts.gt1_last_week, ...
    'caliber', R0ts.caliber, 'mu_variant', R0ts.mu_variant);
summary.legacy_reference_metrics = ref;
summary.consistency_checks = struct();
for k = 1:numel(runs)
    nm = matlab.lang.makeValidName(char(runs{k}.id));
    m = runs{k}.metrics;
    summary.consistency_checks.(nm) = struct( ...
        'R2_log_model', m.R2_log, 'R2_log_legacy', ref.R2_log, ...
        'R2_log_abs_diff', abs(m.R2_log - ref.R2_log), ...
        'R2_raw_model', m.R2_raw, 'R2_raw_legacy', ref.R2_raw, ...
        'R2_raw_abs_diff', abs(m.R2_raw - ref.R2_raw), ...
        'peak_week_model', m.peak_week_model, 'peak_week_legacy', ref.peak_week_model, ...
        'true_total_model', m.true_total, 'true_total_legacy', ref.total_infections_model, ...
        'true_total_rel_diff', abs(m.true_total - ref.total_infections_model) / ref.total_infections_model);
end
mwspm_writejson(fullfile(ctx.dataDir, 'summary.json'), summary);

% ---- payload for manifest.json ------------------------------------------
out = struct();
out.status = 'SUCCESS';
out.report_grid = sprintf('ints 0..%d', p.n_days);
out.internal_step_policy = 'adaptive solver-chosen mesh; fixed-step dt = configured dt_days; report grid never defines the internal step';
out.stiffness = summary.stiffness;
out.solvers = cell(numel(runs), 1);
for k = 1:numel(runs)
    out.solvers{k} = rmfield(runs{k}, {'O', 'S'});
end
out.R0 = summary.R0;
out.data_diagnostics = D.diagnostics;
out.outputs = {'summary.json', 'daily_timeseries.csv', 'weekly_fit.csv', 'r0_daily.csv', 'solver_diagnostics.json'};
out.deterministic = true;
out.rng = 'none (deterministic ODE; no seed-based uncertainty is reported)';
out.consistency_vs_legacy_model_card = summary.consistency_checks;
out.note = ['Pipeline self-test only. It carries no CLAIMS.csv entry and must not be cited as ' ...
    'evidence for CL-06/CL-07; the C5 grid is exp_c5_integrator_grid (T11.01).'];
end

% ---------------------------------------------------------------------------
function v = getfielddef(s, f, dflt)
if isfield(s, f) && ~isempty(s.(f))
    v = s.(f);
else
    v = dflt;
end
end
