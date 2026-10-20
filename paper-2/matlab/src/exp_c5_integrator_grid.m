function out = exp_c5_integrator_grid(ctx)
%EXP_C5_INTEGRATOR_GRID  T11.01 -- C5 integrator diagnostic grid.
%
%   Experiment entry point for config/experiments/exp_c5_integrator_grid.json.
%
%   What is measured
%   ----------------
%   The same deterministic initial-value problem is integrated with
%     dt x method x model x mortality-variant,
%   and every result is compared with an ADAPTIVE reference solution obtained
%   from ode15s at RelTol 1e-10 / AbsTol 1e-12 on the SAME model and the SAME
%   mortality variant.  Nothing here is a "good" or "bad" model: the question
%   is only how much of the difference is produced by the integrator itself.
%
%   This experiment carries CL-06 and CL-07.  It deliberately draws NO verdict
%   about whether the source project's headline numbers are wrong; it supplies
%   the measured grid from which the discussion of C5 is written.
%
%   Separation of concerns: model equations are in mwspm_rhs.m, the integrators
%   in mwspm_solve.m / mwspm_integrate_fixed.m, the observation mapping in
%   mwspm_observe.m.  This file only drives the grid and records the numbers.
%
%   Every row is reported, including rows whose trajectory diverges or fails.
%   Nothing is filtered for being ugly (project rule: do not delete runs that
%   disagree with the expectation, and do not report only the nicest curve).

logFid = ctx.logFid;
root      = ctx.root;
modelCard = ctx.modelCard;
cfg       = ctx.cfg;

fprintf(logFid, '=== exp_c5_integrator_grid (T11.01, C5) ===\n');
fprintf(logFid, 'claims: %s\n', strjoin(cellstr(string(getdef(cfg, 'claims', {}))), ', '));

D = mwspm_load_data(root, modelCard);
fprintf(logFid, 'data: %d days / %d weeks ; observed local total %.0f ; obs peak week %d\n', ...
    D.n_days, D.n_weeks, D.diagnostics.observed_local_total, D.diagnostics.observed_peak_week);

G   = cfg.grid;
nWe = (0:D.n_weeks) * 7;                      % weekly reporting edges

models   = G.models;
if isstruct(models)
    % jsondecode turns a JSON array of uniform objects into a STRUCT ARRAY.
    % Normalise to a cell array here so the indexing below cannot depend on it.
    models = num2cell(models);
end
variants = cellstr(string(G.mu_variants));
dts      = double(G.dts_days(:)');
methods  = cellstr(string(G.methods));
refIds   = {};
solversList = cfg.solvers;
if isstruct(solversList)
    solversList = num2cell(solversList);
end
for k = 1:numel(solversList)
    if strcmpi(char(string(solversList{k}.kind)), 'adaptive')
        refIds{end+1} = char(string(solversList{k}.id)); %#ok<AGROW>
    end
end

if isfield(ctx, 'dryRun') && ctx.dryRun
    out = struct();
    out.dry_run = true;
    out.validated = struct('n_models', numel(models), 'n_mu_variants', numel(variants), ...
        'n_dt', numel(dts), 'n_fixed_methods', sum(~strcmpi(methods, 'ode15s')), ...
        'reference_solvers', {refIds});
    fprintf(logFid, '[dry-run] grid validated: %d models x %d variants x (%d dt x %d fixed methods + reference)\n', ...
        numel(models), numel(variants), numel(dts), sum(~strcmpi(methods, 'ode15s')));
    fprintf(logFid, '[dry-run] no integration performed\n');
    return
end

% ---------------------------------------------------------------- run grid
rows = {};
refStore = struct();
nIntegrations = 0;

for mi = 1:numel(models)
    M = models{mi};
    for vi = 1:numel(variants)
        variant = variants{vi};

        % Season-window fields are OPTIONAL per model entry (P1, 2026-10-07):
        % before this, a model could only say season='window', which silently
        % reused the MODEL CARD's onset/offset/gamma/beta.  Re-running with the
        % T12.01 main report point therefore produced a MIXED parameterisation
        % (new m0/mu_s + old season window), which is exactly what P1 forbids.
        % Each field present in the config now overrides the card.
        ov = struct('m0', double(M.m0), 'mu_s', double(M.mu_s), ...
                    'season_mode', char(M.season));
        optFields = {'onset_doy', 'offset_doy', 'gamma', 'beta', 'rho'};
        for fi = 1:numel(optFields)
            fld = optFields{fi};
            if isfield(M, fld) && ~isempty(M.(fld))
                ov.(fld) = double(M.(fld));
            end
        end
        p = mwspm_params(modelCard, ov);
        F  = mwspm_make_forcing(D, p, struct('interp', 'zoh', 'mu_variant', variant));
        x0 = mwspm_initial_state(D, p, F);
        rhs = @(t, x) mwspm_rhs(t, x, p, F);
        stiffFull = mwspm_stiffness(F, p, 1.0);
        keyName = sprintf('model_%s__variant_%s', matlab.lang.makeValidName(char(string(M.id))), variant);

        fprintf(logFid, ['\n--- model %s (m0=%.6g, mu_s=%.6g, season=%s, onset=%.4f, ' ...
            'offset=%.4f, gamma=%.5f, beta=%.5f, rho=%.6f) | mu_variant=%s ---\n'], ...
            char(string(M.id)), p.m0, p.mu_s, char(string(M.season)), p.onset_doy, ...
            p.offset_doy, p.gamma, p.beta, p.rho, variant);
        fprintf(logFid, 'stiffness ratio=%.4g (fast %.4g / slow %.4g per day); euler stable dt <= %.4g d\n', ...
            stiffFull.stiffness_ratio, stiffFull.rate_fast_per_day, stiffFull.rate_slow_per_day, ...
            stiffFull.dt_euler_stable_max);

        % ---------------- reference solution (adaptive, tight tolerance) ----
        % MaxStep is capped below the forcing grid: both T(t) and the
        % importation Lambda(t) are zero-order held on a DAILY grid, so a step
        % larger than one day can step across a jump.  Without this cap the
        % baseline-model reference came back with 18 steps and zero incidence
        % (see mwspm_solve), which silently turned "reference" into "nothing".
        refSpec = struct('id', 'ode15s_ref', 'kind', 'adaptive', 'method', 'ode15s', ...
            'RelTol', 1e-10, 'AbsTol', 1e-12, 'MaxStep', cfg.reference_solution.max_step_days, ...
            'NonNegative', [1 2 3 4 5 6 7 8], 'mu_variant', variant);
        tRef = tic;
        SR = mwspm_solve(struct('rhs', rhs, 't0', 0, 'tf', p.n_days, 'x0', x0, ...
            'spec', refSpec, 'reportGrid', nWe(:), 'events', cfg.termination.events, ...
            'stiff', stiffFull, 'logFid', logFid));
        nIntegrations = nIntegrations + 1;
        [refDC, refOK] = weekly_incidence(SR, nWe);
        ref = struct();
        ref.model = char(M.id);
        ref.mu_variant = variant;
        ref.method = SR.method;
        ref.n_steps = SR.n_steps;
        ref.min_internal_step = SR.min_internal_step;
        ref.max_internal_step = SR.max_internal_step;
        ref.termination_reason = SR.termination_reason;
        ref.solver_realisation = SR.solver_realisation;
        ref.complete = refOK;
        ref.dC = refDC;
        % The observation layer applies a d-week lag, so the reported total is
        % rho times the SUM OF THE FIRST (n_weeks - d) weekly increments, not
        % rho times the full-window sum: the last d weeks are shifted out of
        % the 47-week window and never appear as a prediction.
        refTrueOnly = refDC(isfinite(refDC));
        ref.reported_total = p.rho * sum(refTrueOnly(1:max(numel(refTrueOnly) - round(p.d_weeks), 0)));
        ref.true_total = sum(refTrueOnly);
        ref.elapsed_seconds = toc(tRef);
        refStore.(keyName) = ref;

        fprintf(logFid, '  REF ode15s: n_steps=%d step in [%.4g, %.4g] d, terminated=%s, true_total=%.3f, elapsed %.2f s\n', ...
            SR.n_steps, SR.min_internal_step, SR.max_internal_step, SR.termination_reason, ...
            ref.true_total, ref.elapsed_seconds);

        % ---------------- fixed-step methods -------------------------------
        for mm = 1:numel(methods)
            method = methods{mm};
            if strcmpi(method, 'ode15s')
                continue
            end
            for di = 1:numel(dts)
                dt = dts(di);
                spec = struct('id', sprintf('%s_dt%g', method, dt), 'kind', 'fixed_step', ...
                    'method', method, 'dt_days', dt, 'clamp_nonnegative', true, ...
                    'mu_variant', variant);
                te = tic;
                S = mwspm_solve(struct('rhs', rhs, 't0', 0, 'tf', p.n_days, 'x0', x0, ...
                    'spec', spec, 'reportGrid', nWe(:), 'events', cfg.termination.events, ...
                    'stiff', stiffFull, 'logFid', logFid));
                nIntegrations = nIntegrations + 1;
                elapsed = toc(te);

                stiff = mwspm_stiffness(F, p, dt);
                if strcmp(variant, 'day')
                    maxMuDt = stiff.max_abs_1_minus_mu_eff_dt_day;
                    nOut = stiff.n_days_outside_stability_day;
                else
                    maxMuDt = stiff.max_abs_1_minus_mu_eff_dt_raw;
                    nOut = stiff.n_days_outside_stability_raw;
                end

                r = score_row(S, nWe, p, D, struct( ...
                    'model', char(string(M.id)), 'mu_variant', variant, 'method', method, ...
                    'dt_days', dt, 'ref', ref, 'elapsed_seconds', elapsed, ...
                    'max_mu_eff_dt', maxMuDt, 'n_days_outside_stability', nOut, ...
                    'stiffness_ratio', stiffFull.stiffness_ratio, ...
                    'dt_euler_stable_max', stiffFull.dt_euler_stable_max));
                rows{end+1, 1} = r; %#ok<AGROW>

                fprintf(logFid, '  %-8s dt=%-6g steps=%-7d term=%-22s true_tot=%-14.3f relerr=%-12.4g R2raw=%-10.4f pk=%-4d outside=%d/%d\n', ...
                    method, dt, S.n_steps, S.termination_reason, r.true_total_model, ...
                    r.true_total_rel_error, r.R2_raw, r.peak_week_model, nOut, D.n_days);
            end
        end
    end
end

% ------------------------------------------------------------------ outputs
nRow = numel(rows);
hdr = {'model', 'mu_variant', 'method', 'dt_days', 'n_steps', 'min_internal_step', ...
       'max_internal_step', 'termination_reason', 'complete', 'n_weeks_scored', ...
       'true_total_model', 'true_total_abs_error', 'true_total_rel_error', ...
       'reported_total_model', 'reported_total_rel_error', 'reference_degenerate', ...
       'R2_log', 'R2_raw', 'RMSE_log', 'peak_week_model', 'peak_week_shift_vs_reference', ...
       'max_mu_eff_dt', 'n_days_outside_euler_stability', 'stiffness_ratio', ...
       'dt_euler_stable_max', 'elapsed_seconds'};

flat = nan(nRow, numel(hdr));
cells = cell(nRow, 4);
for i = 1:nRow
    r = rows{i};
    cells(i, :) = {r.model, r.mu_variant, r.method, r.termination_reason};
    flat(i, :) = [NaN, NaN, NaN, r.dt_days, r.n_steps, r.min_internal_step, ...
        r.max_internal_step, NaN, double(r.complete), r.n_weeks_scored, ...
        r.true_total_model, r.true_total_abs_error, r.true_total_rel_error, ...
        r.reported_total_model, r.reported_total_rel_error, double(r.reference_degenerate), ...
        r.R2_log, r.R2_raw, r.RMSE_log, r.peak_week_model, r.peak_week_shift, ...
        r.max_mu_eff_dt, r.n_days_outside_stability, r.stiffness_ratio, ...
        r.dt_euler_stable_max, r.elapsed_seconds];
    flat(i, 1:3) = [mi_code(r.model), mv_code(r.mu_variant), me_code(r.method)];
end
% NOTE (MATLAB trap 3 revisited): writecsv takes a numeric matrix, so the three
% leading character columns are encoded numerically and the mapping is written
% into the summary; lossy encodings are not acceptable for a machine-readable
% result.  A second, fully explicit long-format CSV is written below.
numCsv = flat;
mwspm_writecsv(fullfile(ctx.dataDir, 'grid_results.csv'), hdr, numCsv);

longRows = {};
for i = 1:nRow
    r = rows{i};
    longRows{end+1, 1} = sprintf('%s,%s,%s,%g,%s,%d,%.17g,%.17g,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%d,%.17g,%.17g,%.17g,%d,%d,%.17g,%d,%.17g,%.17g,%.4f', ...
        r.model, r.mu_variant, r.method, r.dt_days, r.termination_reason, r.n_steps, ...
        r.min_internal_step, r.max_internal_step, r.complete, r.n_weeks_scored, ...
        r.true_total_model, r.true_total_abs_error, r.true_total_rel_error, ...
        r.reported_total_model, r.reported_total_rel_error, double(r.reference_degenerate), ...
        r.R2_log, r.R2_raw, r.RMSE_log, ...
        r.peak_week_model, r.peak_week_shift, r.max_mu_eff_dt, ...
        r.n_days_outside_stability, r.stiffness_ratio, r.dt_euler_stable_max, r.elapsed_seconds); %#ok<AGROW>
end
fid = fopen(fullfile(ctx.dataDir, 'grid_results_long.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '%s\n', strjoin(hdr, ','));
for i = 1:numel(longRows)
    fprintf(fid, '%s\n', longRows{i});
end
fclose(fid);

% reference weekly incidence, one column per (model, variant).
% NOTE: there are 48 weekly EDGES but only 47 weekly increments, so every
% incidence column is padded with one trailing NaN to keep the row count
% aligned with the edge column.
refCols = nWe(:);
refHdr = {'week_edge_day'};
fn = fieldnames(refStore);
for k = 1:numel(fn)
    R = refStore.(fn{k});
    dcol = [double(R.dC(:)); NaN];
    assert(numel(dcol) == numel(nWe), 'c5:grid', 'column length mismatch');
    refCols = [refCols, dcol]; %#ok<AGROW>
    refHdr{end+1} = sprintf('true_incidence_%s', fn{k}); %#ok<AGROW>
end
mwspm_writecsv(fullfile(ctx.dataDir, 'reference_weekly.csv'), refHdr, refCols);

diagOut = struct();
diagOut.reference = refStore;
diagOut.reference_definition = cfg.reference_solution;
diagOut.n_integrations = nIntegrations;
diagOut.grid = G;
diagOut.row_note = 'complete=0 means the trajectory did not reach t_end (divergence, NaN/Inf or a model thinning event); such rows keep their partial totals and get NaN for every R2.';
mwspm_writejson(fullfile(ctx.dataDir, 'solver_diagnostics.json'), diagOut);

% ------------------------------------------------------------------ summary
mxDev = max(cellfun(@(r) absOrZero(r.true_total_rel_error), rows));
summary = struct();
summary.experiment_id = ctx.expId;
summary.task = char(string(getdef(cfg, 'task', '-')));
summary.n_rows = nRow;
summary.n_integrations = nIntegrations;
summary.categorical_codes = struct('model', 'see grid_results_long.csv; models are named there', ...
    'mu_variant', 'raw | day', 'method', 'euler | rk4', ...
    'note', ['grid_results.csv encodes model/mu_variant/method numerically because ' ...
             'mwspm_writecsv writes a numeric matrix; grid_results_long.csv carries the ' ...
             'same rows with every categorical column spelled out and is the authoritative file.']);
summary.max_abs_true_total_rel_error = mxDev;
summary.n_rows_not_complete = sum(cellfun(@(r) ~r.complete, rows));
summary.rows_with_days_outside_stability = sum(cellfun(@(r) r.n_days_outside_stability > 0, rows));
summary.headline_rows = pick_headline(rows);
summary.honesty = ['This is a measurement, not a verdict. No row is suppressed. Rows whose ' ...
    'trajectory diverged keep their partial totals and get NaN goodness-of-fit values.'];
mwspm_writejson(fullfile(ctx.dataDir, 'summary.json'), summary);

out = struct();
out.status = 'SUCCESS';
out.n_rows = nRow;
out.n_integrations = nIntegrations;
out.max_abs_true_total_rel_error = mxDev;
out.n_rows_not_complete = summary.n_rows_not_complete;
out.headline_rows = summary.headline_rows;
out.outputs = {'grid_results.csv', 'grid_results_long.csv', 'reference_weekly.csv', ...
    'solver_diagnostics.json', 'summary.json'};
out.deterministic = true;
out.rng = 'none (deterministic ODE; grid over numerical methods, no stochastic replicate)';
out.note = ['T11.01 supplies the C5 measurement grid. It does not itself decide whether the ' ...
    'source project total is wrong; see the honesty block of the config before quoting numbers.'];
end

% ---------------------------------------------------------------------------
function [dC, complete] = weekly_incidence(S, edges)
% Weekly true incidence by differencing the cumulative accumulator exactly at
% the weekly edges (quadrature-exact; see mwspm_observe).  Works on the
% INTERNAL mesh so that a run which stopped early is not silently extrapolated.
tI = double(S.t_internal(:));
cI = double(S.x_internal(:, 8));
m  = isfinite(tI) & isfinite(cI);
tI = tI(m);
cI = cI(m);
edges = double(edges(:));
dC = nan(numel(edges) - 1, 1);
if numel(tI) < 2
    complete = false;
    return
end
Ce = nan(numel(edges), 1);
ok = edges >= tI(1) - 1e-12 & edges <= tI(end) + 1e-12;
Ce(ok) = interp1(tI, cI, edges(ok), 'linear', NaN);
dC = diff(Ce);
complete = all(isfinite(dC));
end

% ---------------------------------------------------------------------------
function r = score_row(S, edges, p, D, a)
[dC, complete] = weekly_incidence(S, edges);
fin = isfinite(dC);
nScored = sum(fin);

trueTotal = sum(dC(fin));
refTrue   = a.ref.true_total;
refRep    = a.ref.reported_total;

d = round(p.d_weeks);
d = min(max(d, 0), p.n_weeks);
dClag = [zeros(d, 1); dC(1:end - d)];
Yhat = p.rho * dClag;
repTotal = sum(Yhat(isfinite(Yhat)));

if complete
    m = mwspm_metrics(D.Y_weekly(:), Yhat);
    R2log = m.R2_log; R2raw = m.R2_raw; rmse = m.RMSE_log;
    pk = m.peak_week_model;
    refPk = find(a.ref.dC == max(a.ref.dC), 1);
    if isempty(refPk), refPk = NaN; end
    pkShift = pk - refPk;
else
    % Comparing an incomplete trajectory against the full 47-week observation
    % series would manufacture a number that looks like a fit quality but is
    % not one.  Report NaN and say why.
    R2log = NaN; R2raw = NaN; rmse = NaN; pk = NaN; pkShift = NaN;
end

r = struct();
r.model            = a.model;
r.mu_variant       = a.mu_variant;
r.method           = a.method;
r.dt_days          = a.dt_days;
r.n_steps          = S.n_steps;
r.min_internal_step = S.min_internal_step;
r.max_internal_step = S.max_internal_step;
r.termination_reason = S.termination_reason;
r.solver_realisation = S.solver_realisation;
r.complete         = complete;
r.n_weeks_scored   = nScored;
r.true_total_model = trueTotal;
r.true_total_abs_error = trueTotal - refTrue;
r.true_total_rel_error = pickFinite((trueTotal - refTrue) / max(abs(refTrue), 1e-30), abs(refTrue) >= 1e-9);
r.reported_total_model = repTotal;
% Dividing by a reference that accumulated (almost) nothing produces numbers
% like 3.7e32, which look like a result but encode "the reference is useless
% here".  Report NaN with an explicit flag instead.
if abs(refRep) < 1e-9
    r.reported_total_rel_error = NaN;
else
    r.reported_total_rel_error = (repTotal - refRep) / abs(refRep);
end
r.reference_degenerate = (abs(refTrue) < 1e-9) || (abs(refRep) < 1e-9);
r.R2_log           = R2log;
r.R2_raw           = R2raw;
r.RMSE_log         = rmse;
r.peak_week_model  = pk;
r.peak_week_shift  = pkShift;
r.max_mu_eff_dt    = a.max_mu_eff_dt;
r.n_days_outside_stability = a.n_days_outside_stability;
r.stiffness_ratio  = a.stiffness_ratio;
r.dt_euler_stable_max = a.dt_euler_stable_max;
r.elapsed_seconds  = a.elapsed_seconds;
end

% ---------------------------------------------------------------------------
function H = pick_headline(rows)
% The rows an author might be tempted to quote alone.  They are returned so
% that the summary states them explicitly rather than leaving them implicit.
H = struct();
for i = 1:numel(rows)
    r = rows{i};
    if strcmp(r.model, 'mwspm_gen3') && strcmp(r.mu_variant, 'day') && ...
            strcmp(r.method, 'euler') && abs(r.dt_days - 1.0) < 1e-12
        H.source_faithful_row = r;
    end
    if strcmp(r.model, 'mwspm_gen3') && strcmp(r.mu_variant, 'day') && ...
            strcmp(r.method, 'rk4') && abs(r.dt_days - 0.01) < 1e-12
        H.finest_fixed_step_row = r;
    end
end
end

% ---------------------------------------------------------------------------
function v = mi_code(model)
switch model
    case 'baseline_m0.249', v = 1;
    case 'mwspm_gen3',      v = 2;
    otherwise,              v = 0;
end
end
function v = mv_code(vr)
switch vr
    case 'raw', v = 1;
    case 'day', v = 2;
    otherwise,  v = 0;
end
end
function v = me_code(me)
switch me
    case 'euler', v = 1;
    case 'rk4',   v = 2;
    otherwise,    v = 0;
end
end
function v = pickFinite(x, keep)
if keep && isfinite(x), v = x; else v = NaN; end
end
function v = absOrZero(x)
if isfinite(x), v = abs(x); else v = 0; end
end
function v = getdef(s, f, dflt)
if isfield(s, f) && ~isempty(s.(f))
    v = s.(f);
else
    v = dflt;
end
end
