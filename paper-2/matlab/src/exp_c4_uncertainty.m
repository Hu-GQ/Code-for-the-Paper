function out = exp_c4_uncertainty(ctx)
%EXP_C4_UNCERTAINTY  T14.01 -- propagate the rho profile to the warning week.
%
%   Experiment entry point for config/experiments/exp_c4_uncertainty.json.
%
%   Three propagation channels, run in parallel (config "design"):
%     (1) rho continuous profile scan on the FROZEN grid [0.30, 1.00] step 0.01,
%         reported as TWO branches separately:
%           (a) profile_refit      : at each rho the other parameters are the
%                                    T12.01 refit (rho_profile.csv); this is the
%                                    channel the gate G1 measured, so the frozen
%                                    grid is READ from gate_g1_weeks.csv and
%                                    CROSS-CHECKED by an identical recomputation;
%           (b) main_report_point  : the main report point trajectory is held
%                                    fixed and rho is varied; by rem:rho-decouple
%                                    tau* is then structurally rho-invariant -
%                                    the measured spread MUST be ~0 and is
%                                    reported, not assumed.
%         Supplement (honesty): the T12.01 Wilks 95% band for rho was MEASURED
%         to cover the whole [0.02, 1.00] (lag1), so the interval over the full
%         band is computed as well; the frozen grid [0.30, 1.00] is prior-
%         specified, NOT the estimated band (G0 ruling, 2026-10-06).
%     (2) external-evidence extreme rho = 0.13 (Undurraga EF = 7.6) as a single
%         scenario, on the profile-refit theta (both lag branches) and on the
%         fixed main-report-point trajectory.
%     (3) one-at-a-time perturbations of the main report point: m0 +/-10%,
%         mu_s +/-10% (four edge scenarios) and gamma_h in {1/4, 1/8} per day,
%         i.e. infectious periods of 4 and 8 days - the endpoints of the Cheng
%         2016 range, whose coverage is an obligation (gamma_h = 1/5 is only one
%         interior point and may NOT be reported alone).
%
%   Warning-week definition (identical to T13.01 / gate G1):
%     tau* = min{ w : cum_w >= T0(w) },  T0(w) = rho / log R0(w)  (R0(w) > 1),
%     cum_w = cumulative DETECTED model cases, weekly R0 = MEAN of daily R0,
%     trajectory = source-faithful Euler dt = 1 day, mu_variant = 'day'
%     (mwspm_calib_kernel), R0 uses the CONTINUOUS (raw) mortality variant.
%
%   HONESTY (config "honesty", CL-03/CL-12 prohibitions):
%     * rho boundary values are NOT "the true reporting rate";
%     * OAT perturbations are NOT a joint confidence region - wording must be
%       "sensitivity under one-at-a-time perturbation";
%     * the propagated interval is a TRIVIAL consequence of the rho-m0 ridge:
%       on the refit branch the detected trajectory is pinned to the observed
%       cumulative curve, so the interval is narrow by construction, NOT
%       because rho is well identified (it is completely unidentifiable);
%     * if the interval is incompatible with the observation start week (19),
%       that is reported as-is (CL-03 failure handling).

logFid = ctx.logFid;
root   = ctx.root;
cfg    = ctx.cfg;
modelCard = ctx.modelCard;

fprintf(logFid, '=== exp_c4_uncertainty (T14.01, rho profile -> warning week) ===\n');
fprintf(logFid, 'claims: %s\n', strjoin(cellstr(string(getdef(cfg, 'claims', {}))), ', '));

D   = mwspm_load_data(root, modelCard);
pre = mwspm_calib_prepare(D, modelCard, struct('mu_variant', 'day'));
pRef = mwspm_params(modelCard);

% ---------------------------------------------------------- parent runs -----
fitRun = unc_find_parent(root, getdef(cfg, 'parent_run_id', []), 'exp_c4_rho_fit', 'data/fit_best.json');
assert(~isempty(fitRun), 'uncertainty:noParentFit', ...
    ['no SUCCESS run of exp_c4_rho_fit with data/fit_best.json found under results/runs/. ' ...
     'T14.01 must consume a real fit; placeholder numbers are forbidden.']);
profPath = fullfile(root, 'results', 'runs', fitRun.runId, 'data', 'rho_profile.csv');
assert(exist(profPath, 'file') == 2, 'uncertainty:noParentProfile', ...
    'parent fit run %s has no data/rho_profile.csv', fitRun.runId);

gateRun = unc_find_parent(root, [], 'exp_c4_gate_g1', 'data/gate_g1_weeks.csv');
assert(~isempty(gateRun), 'uncertainty:noParentGate', ...
    ['no SUCCESS run of exp_c4_gate_g1 with data/gate_g1_weeks.csv found under results/runs/. ' ...
     'The frozen-grid profile channel is read from the gate; placeholder numbers are forbidden.']);

fprintf(logFid, 'parent fit run   : %s\n', fitRun.runId);
fprintf(logFid, 'parent gate run  : %s\n', gateRun.runId);

fit  = mwspm_readjson(fitRun.path);
prof = unc_read_profile(profPath);
gate = unc_read_gate_weeks(fullfile(root, 'results', 'runs', gateRun.runId, 'data', 'gate_g1_weeks.csv'));

% Observation start week: recomputed from the observations (definition), then
% cross-checked against the gate file (which carries the same quantity).
obsStart = find(pre.Y > 0, 1);
assert(~isempty(obsStart), 'uncertainty:noCases', 'observed series has no positive week');
gateObs = unique(gate.observation_start_week(isfinite(gate.observation_start_week)));
assert(isscalar(gateObs) && gateObs == obsStart, 'uncertainty:obsStartMismatch', ...
    'recomputed observation start week (%d) disagrees with the gate file (%s)', ...
    obsStart, mat2str(gateObs));
fprintf(logFid, 'observation start week (recomputed, matches gate) = %d\n', obsStart);

% -------------------------------------------------------------- grids -------
PD = cfg.parameter_domain;
rp = PD.rho_profile;
rhoFrozen = double(rp.from) : double(rp.step) : double(rp.to);
rhoFrozen = round(rhoFrozen, 6);
rhoExt    = double(PD.rho_external(:)');
m0Rel     = double(PD.m0_perturbation_rel(:)');
musRel    = double(PD.mu_s_perturbation_rel(:)');
gammaDays = double(PD.gamma_h_sensitivity_days(:)');

if isfield(ctx, 'dryRun') && ctx.dryRun
    out = struct();
    out.dry_run = true;
    out.validated = struct('n_frozen_grid', numel(rhoFrozen), 'rho_external', rhoExt, ...
        'n_profile_rows', numel(prof.branch_lag), 'gate_rows', numel(gate.branch_lag), ...
        'parent_fit_run_id', fitRun.runId, 'parent_gate_run_id', gateRun.runId, ...
        'observation_start_week', obsStart);
    fprintf(logFid, '[dry-run] grids and parents readable; nothing recomputed\n');
    return
end

% --------------------------------------------------- baseline (main point) --
th0 = unc_theta_from_fit(fit.main_report_point.theta);
gammaH0 = pRef.gamma_h;
fprintf(logFid, 'main report point: lag %d, m0=%.6f rho=%.6f mu_s=%.6f onset=%.3f offset=%.3f\n', ...
    th0.d_weeks, th0.m0, th0.rho, th0.mu_s, th0.onset_doy, th0.offset_doy);
fprintf(logFid, 'baseline gamma_h = 1/%.4g per day (infectious period %.4g d)\n', 1 / gammaH0, 1 / gammaH0);

[tauBase, r0sBase, T0Base, cumBase] = unc_tau_star(pre, D, th0, th0.rho, gammaH0);
fprintf(logFid, 'baseline warning week (fixed rho = %.6f, no propagation) = %s\n', ...
    th0.rho, num2strOrNA(tauBase));

rows = {};   % warning_week_vs_rho.csv rows

% ===================================== channel (1b): main report point ======
fprintf(logFid, '\n--- channel (1b): fixed main-report-point trajectory, rho varies ---\n');
rhoMain = [rhoExt, rhoFrozen];
tauMain = nan(size(rhoMain));
for i = 1:numel(rhoMain)
    rg = rhoMain(i);
    thE = th0; thE.rho = rg;
    [tw, r0s, T0a, ca] = unc_tau_star(pre, D, thE, rg, gammaH0);
    tauMain(i) = tw;
    rows{end + 1, 1} = {'main_report_point', th0.d_weeks, rg, double(rg >= rhoFrozen(1) - 1e-9), ...
        tw, T0a, ca, r0s, 'main report point theta (fixed trajectory)', ...
        'structurally rho-invariant (rem:rho-decouple); measured'}; %#ok<AGROW>
end
finM = isfinite(tauMain);
spreadMain = tern(any(finM), max(tauMain(finM)) - min(tauMain(finM)), NaN);
fprintf(logFid, '  tau* over %d rho values: [%s, %s], spread = %s weeks (rem:rho-decouple => must be ~0)\n', ...
    numel(rhoMain), num2strOrNA(min(tauMain(finM))), num2strOrNA(max(tauMain(finM))), num2strOrNA(spreadMain));
if isfinite(spreadMain) && spreadMain > 0
    fprintf(logFid, ['  WARNING: nonzero spread on a FIXED trajectory contradicts the exact ' ...
        'rho-cancellation; it can only be a numerical artifact and is reported as measured.\n']);
end

% ====================================== channel (1a): profile refit =========
fprintf(logFid, '\n--- channel (1a): profile refit branches (full profile grid, frozen grid cross-checked vs gate) ---\n');
lags = unique(prof.branch_lag(:)');
refit = struct();
nCrossChecked = 0;
for li = 1:numel(lags)
    dd = lags(li);
    key = sprintf('lag%d', dd);
    sel = find(prof.branch_lag == dd);
    gvals = prof.grid_value(sel);
    nG = numel(sel);
    tauL = nan(nG, 1);
    T0L = nan(nG, 1); cumL = nan(nG, 1); r0sL = nan(nG, 1);
    for k = 1:nG
        th = struct('m0', prof.m0(sel(k)), 'rho', prof.rho(sel(k)), ...
            'onset_doy', prof.onset_doy(sel(k)), 'offset_doy', prof.offset_doy(sel(k)), ...
            'mu_s', prof.mu_s(sel(k)), 'gamma', prof.gamma(sel(k)), ...
            'beta', prof.beta(sel(k)), 'd_weeks', prof.d_weeks(sel(k)));
        [tauL(k), r0sL(k), T0L(k), cumL(k)] = unc_tau_star(pre, D, th, gvals(k), gammaH0);
        rows{end + 1, 1} = {'profile_refit', dd, gvals(k), double(gvals(k) >= rhoFrozen(1) - 1e-9), ...
            tauL(k), T0L(k), cumL(k), r0sL(k), 'rho_profile.csv refit theta (recomputed)', ...
            'matches gate_g1_weeks.csv on the frozen grid (asserted)'}; %#ok<AGROW>
    end
    % ---- cross-check against the gate file on the frozen grid -------------
    gsel = find(gate.branch_lag == dd);
    for i = 1:numel(rhoFrozen)
        rg = rhoFrozen(i);
        [~, jp] = min(abs(gvals - rg));
        assert(abs(gvals(jp) - rg) <= 1e-6, 'uncertainty:gridMismatch', ...
            'rho=%g is not on the parent profile grid (closest %g)', rg, gvals(jp));
        [~, jg] = min(abs(gate.rho(gsel) - rg));
        assert(abs(gate.rho(gsel(jg)) - rg) <= 1e-6, 'uncertainty:gateGridMismatch', ...
            'rho=%g is not in gate_g1_weeks.csv (closest %g)', rg, gate.rho(gsel(jg)));
        tauGate = gate.tau_star_week(gsel(jg));
        assert(isequaln(tauL(jp), tauGate), 'uncertainty:caliberDrift', ...
            ['recomputed tau* (%g) at lag %d rho %.2f differs from the gate value (%g). ' ...
             'The T14.01 caliber must be IDENTICAL to T13.01; refusing to mix calibers.'], ...
            tauL(jp), dd, rg, tauGate);
        nCrossChecked = nCrossChecked + 1;
    end
    finF = isfinite(tauL) & (gvals >= rhoFrozen(1) - 1e-9);
    finAll = isfinite(tauL);
    refit.(key) = struct('grid', gvals(:)', 'tau', tauL(:)', ...
        'frozen_low', tern(any(finF), min(tauL(finF)), NaN), ...
        'frozen_high', tern(any(finF), max(tauL(finF)), NaN), ...
        'full_low', tern(any(finAll), min(tauL(finAll)), NaN), ...
        'full_high', tern(any(finAll), max(tauL(finAll)), NaN));
    fprintf(logFid, '  %s: frozen grid [%g,%g] -> tau* in [%s, %s]; full profile grid [%g,%g] -> tau* in [%s, %s]\n', ...
        key, rhoFrozen(1), rhoFrozen(end), num2strOrNA(refit.(key).frozen_low), num2strOrNA(refit.(key).frozen_high), ...
        gvals(1), gvals(end), num2strOrNA(refit.(key).full_low), num2strOrNA(refit.(key).full_high));
end
fprintf(logFid, '  cross-check: %d (lag, rho) recomputations match gate_g1_weeks.csv EXACTLY\n', nCrossChecked);

% ========================================= channel (2): external rho=0.13 ===
fprintf(logFid, '\n--- channel (2): external extreme rho = %s ---\n', mat2str(rhoExt));
tauExt = struct();
for i = 1:numel(rhoExt)
    rg = rhoExt(i);
    for li = 1:numel(lags)
        dd = lags(li);
        key = sprintf('lag%d', dd);
        sel = find(prof.branch_lag == dd);
        [~, jp] = min(abs(prof.grid_value(sel) - rg));
        assert(abs(prof.grid_value(sel(jp)) - rg) <= 1e-6, 'uncertainty:extOffGrid', ...
            'external rho=%g is not on the parent profile grid', rg);
        g = prof.grid_value(sel(jp));
        tauExt.(key) = refit.(key).tau(prof.grid_value(sel) == g);
        tauExt.(key) = tauExt.(key)(1);
        fprintf(logFid, '  %s refit at rho=%.2f : tau* = %s\n', key, rg, num2strOrNA(tauExt.(key)));
    end
    [~, jm] = min(abs(rhoMain - rg));
    tauExt.main_report_point = tauMain(jm);
    fprintf(logFid, '  fixed main-report-point trajectory at rho=%.2f : tau* = %s (same constant as every rho on a fixed trajectory)\n', ...
        rg, num2strOrNA(tauExt.main_report_point));
end

% ========================================= channel (3): OAT perturbations ===
fprintf(logFid, '\n--- channel (3): one-at-a-time perturbations of the main report point ---\n');
oatRows = {};
oatRows{end + 1, 1} = {'baseline', 'main_report_point', '-', th0.rho, th0.rho, ...
    tauBase, 0, T0Base, cumBase, r0sBase, ...
    'main report point, fixed rho=0.800, no propagation (= current caliber)'}; %#ok<AGROW>

scen = {};
for i = 1:numel(m0Rel)
    thE = th0; thE.m0 = th0.m0 * (1 + m0Rel(i));
    scen{end + 1, 1} = {tern(m0Rel(i) < 0, 'm0_minus10', 'm0_plus10'), 'm0', 'ratio to human density', ...
        th0.m0, thE.m0, thE, gammaH0, sprintf('m0 x %.2f, everything else at the main report point', 1 + m0Rel(i))}; %#ok<AGROW>
end
for i = 1:numel(musRel)
    thE = th0; thE.mu_s = th0.mu_s * (1 + musRel(i));
    scen{end + 1, 1} = {tern(musRel(i) < 0, 'mu_s_minus10', 'mu_s_plus10'), 'mu_s', 'mortality scale', ...
        th0.mu_s, thE.mu_s, thE, gammaH0, sprintf('mu_s x %.2f, everything else at the main report point', 1 + musRel(i))}; %#ok<AGROW>
end
for i = 1:numel(gammaDays)
    gH = 1 / gammaDays(i);
    scen{end + 1, 1} = {sprintf('gamma_h_%dd', round(gammaDays(i))), 'gamma_h', 'infectious period, days', ...
        1 / gammaH0, gammaDays(i), th0, gH, ...
        sprintf('gamma_h = 1/%.4g per day (endpoint of the Cheng 2016 4-8 d range), everything else at the main report point', gammaDays(i))}; %#ok<AGROW>
end

oatTau = struct('baseline', tauBase);
for i = 1:numel(scen)
    sc = scen{i};
    [tw, r0s, T0a, ca] = unc_tau_star(pre, D, sc{6}, sc{6}.rho, sc{7});
    oatTau.(sc{1}) = tw;
    oatRows{end + 1, 1} = {sc{1}, sc{2}, sc{3}, sc{4}, sc{5}, tw, ...
        tern(isfinite(tw) && isfinite(tauBase), tw - tauBase, NaN), T0a, ca, r0s, sc{8}}; %#ok<AGROW>
    fprintf(logFid, '  %-14s : tau* = %s (delta vs baseline %s weeks)\n', sc{1}, ...
        num2strOrNA(tw), num2strOrNA(tern(isfinite(tw) && isfinite(tauBase), tw - tauBase, NaN)));
end

% ------------------------------------------------------------ metrics -------
lag1 = refit.lag1;
lo = lag1.frozen_low; hi = lag1.frozen_high;
compatible = isfinite(lo) && isfinite(hi) && (lo <= obsStart) && (obsStart <= hi);

metrics = struct();
metrics.warning_week_interval_low = lo;
metrics.warning_week_interval_high = hi;
metrics.warning_week_interval_width_weeks = tern(isfinite(lo) && isfinite(hi), hi - lo, NaN);
metrics.warning_week_at_rho_0_13 = tauExt.lag1;      % config name: warning_week_at_rho_0.13
metrics.warning_week_m0_minus10 = oatTau.m0_minus10;
metrics.warning_week_m0_plus10 = oatTau.m0_plus10;
metrics.warning_week_mu_s_minus10 = oatTau.mu_s_minus10;
metrics.warning_week_mu_s_plus10 = oatTau.mu_s_plus10;
metrics.warning_week_gamma_h_4d = oatTau.gamma_h_4d;
metrics.warning_week_gamma_h_8d = oatTau.gamma_h_8d;
metrics.current_caliber_peak_week = tauBase;
metrics.compatible_with_observation_start_week = compatible;

fprintf(logFid, '\n--- metrics ---\n');
fprintf(logFid, '  warning-week interval (profile refit lag1, frozen grid [%g,%g]) = [%s, %s], width %s weeks\n', ...
    rhoFrozen(1), rhoFrozen(end), num2strOrNA(lo), num2strOrNA(hi), num2strOrNA(metrics.warning_week_interval_width_weeks));
fprintf(logFid, '  compatible with observation start week %d : %s\n', obsStart, yesNo(compatible));

supplements = struct();
supplements.interval_lag2_frozen = [refit.lag2.frozen_low, refit.lag2.frozen_high];
supplements.interval_union_lags_frozen = [min([refit.lag1.frozen_low, refit.lag2.frozen_low]), ...
    max([refit.lag1.frozen_high, refit.lag2.frozen_high])];
supplements.interval_lag1_full_profile_grid = [refit.lag1.full_low, refit.lag1.full_high];
supplements.interval_lag2_full_profile_grid = [refit.lag2.full_low, refit.lag2.full_high];
supplements.main_report_point_tau_star_spread_over_rho = spreadMain;
supplements.warning_week_at_rho_0_13_lag2 = tauExt.lag2;
supplements.warning_week_at_rho_0_13_fixed_trajectory = tauExt.main_report_point;
supplements.oat_max_abs_delta_weeks = max(abs([oatTau.m0_minus10, oatTau.m0_plus10, ...
    oatTau.mu_s_minus10, oatTau.mu_s_plus10, oatTau.gamma_h_4d, oatTau.gamma_h_8d] - tauBase));
supplements.all_channels_min_max = [min([lo, refit.lag2.frozen_low, tauBase, ...
    oatTau.m0_minus10, oatTau.m0_plus10, oatTau.mu_s_minus10, oatTau.mu_s_plus10, ...
    oatTau.gamma_h_4d, oatTau.gamma_h_8d]), max([hi, refit.lag2.frozen_high, tauBase, ...
    oatTau.m0_minus10, oatTau.m0_plus10, oatTau.mu_s_minus10, oatTau.mu_s_plus10, ...
    oatTau.gamma_h_4d, oatTau.gamma_h_8d])];
supplements.compatible_by_branch = struct( ...
    'lag1_frozen', compatible, ...
    'lag2_frozen', (refit.lag2.frozen_low <= obsStart) && (obsStart <= refit.lag2.frozen_high), ...
    'lag1_full_profile_grid', (refit.lag1.full_low <= obsStart) && (obsStart <= refit.lag1.full_high), ...
    'lag2_full_profile_grid', (refit.lag2.full_low <= obsStart) && (obsStart <= refit.lag2.full_high));

% ------------------------------------------------------------- outputs ------
fid = fopen(fullfile(ctx.dataDir, 'warning_week_vs_rho.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, ['channel,branch_lag,rho,in_frozen_grid,tau_star_week,T0_at_tau_star,', ...
    'cum_detected_at_tau_star,R0_gt1_start_week,theta_source,cross_check\n']);
for i = 1:numel(rows)
    rw = rows{i};
    fprintf(fid, '%s,%d,%.17g,%d,%s,%s,%s,%s,"%s","%s"\n', rw{1}, rw{2}, rw{3}, rw{4}, ...
        numOrNA(rw{5}), numOrNA(rw{6}), numOrNA(rw{7}), numOrNA(rw{8}), rw{9}, rw{10});
end
fclose(fid);

fid = fopen(fullfile(ctx.dataDir, 'warning_week_oat.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, ['scenario,parameter,unit,baseline_value,scenario_value,tau_star_week,', ...
    'delta_vs_baseline_weeks,T0_at_tau_star,cum_detected_at_tau_star,R0_gt1_start_week,note\n']);
for i = 1:numel(oatRows)
    rw = oatRows{i};
    fprintf(fid, '%s,%s,"%s",%.17g,%.17g,%s,%s,%s,%s,%s,"%s"\n', rw{1}, rw{2}, rw{3}, rw{4}, rw{5}, ...
        numOrNA(rw{6}), numOrNA(rw{7}), numOrNA(rw{8}), numOrNA(rw{9}), numOrNA(rw{10}), rw{11});
end
fclose(fid);

usum = struct();
usum.experiment_id = ctx.expId;
usum.task = char(string(getdef(cfg, 'task', '-')));
usum.provenance = struct('parent_fit_run_id', fitRun.runId, 'parent_gate_run_id', gateRun.runId, ...
    'fit_best_json', fitRun.path, 'rho_profile_csv', profPath, ...
    'gate_g1_weeks_csv', fullfile(root, 'results', 'runs', gateRun.runId, 'data', 'gate_g1_weeks.csv'));
usum.definition = struct('warning_week', 'tau* = min{ w : cum_w >= T0(w) }, T0(w) = rho/log R0(w) for R0(w) > 1', ...
    'cum_w', 'cumulative DETECTED model cases (mwspm_calib_kernel, Euler dt=1 day, mu_variant=day)', ...
    'weekly_R0', 'MEAN of daily R0 inside the week; R0 uses the CONTINUOUS (raw) mortality variant', ...
    'observation_start_week', obsStart, ...
    'current_caliber_peak_week_definition', ['warning week at the main report point with fixed rho=0.800 and NO ' ...
        'propagation (the current manuscript''s single-number answer). Named peak_week because the current ' ...
        'manuscript phrases that week as the establishment-probability peak week (CL-09); under the T0 ' ...
        'caliber it is tau*. It is NOT the incidence peak week (fitted incidence peaks at week 37).']);
usum.metrics = metrics;
usum.metric_name_map = struct('warning_week_at_rho_0_13', 'config metric warning_week_at_rho_0.13 (struct fields cannot contain a dot)');
usum.value_types = struct('warning_week_interval_low', 'parameter_uncertainty_interval', ...
    'warning_week_interval_high', 'parameter_uncertainty_interval', ...
    'warning_week_interval_width_weeks', 'parameter_uncertainty_interval', ...
    'warning_week_at_rho_0_13', 'deterministic_point', ...
    'warning_week_m0_minus10', 'deterministic_point', 'warning_week_m0_plus10', 'deterministic_point', ...
    'warning_week_mu_s_minus10', 'deterministic_point', 'warning_week_mu_s_plus10', 'deterministic_point', ...
    'warning_week_gamma_h_4d', 'deterministic_point', 'warning_week_gamma_h_8d', 'deterministic_point', ...
    'current_caliber_peak_week', 'deterministic_point', ...
    'compatible_with_observation_start_week', 'deterministic_point');
usum.supplements = supplements;
usum.frozen_grid = struct('rho_from', rhoFrozen(1), 'rho_to', rhoFrozen(end), 'step', double(rp.step), ...
    'note', 'prior-specified gate grid; the MEASURED T12.01 Wilks 95% band is [0.02, 1.00] (lag1)');
usum.cross_check = struct('rule', 'recomputed tau* on the frozen grid must equal gate_g1_weeks.csv exactly', ...
    'n_points_checked', nCrossChecked, 'result', 'PASS (no caliber drift between T13.01 and T14.01)');
usum.gate_controls = struct('current_caliber_100_week', unique(gate.current_caliber_100_week(isfinite(gate.current_caliber_100_week))), ...
    'observation_start_week', obsStart);
usum.honesty = { ...
    'the propagated interval is NOT evidence that rho is well identified: the T12.01 Wilks 95% band covers the whole [0.02, 1.00] (rho completely unidentifiable) and on the refit branch the detected trajectory is pinned to the observed cumulative curve, so the narrow interval is TRIVIAL by construction (CL-03 wording constraint)'; ...
    'rho boundary values are not interpreted as the true reporting rate (CL-12 prohibition)'; ...
    'OAT perturbations are one-at-a-time sensitivities, NOT a joint confidence region'; ...
    'gamma_h = 1/4 and 1/8 per day (4 d / 8 d) cover the endpoints of the Cheng 2016 4-8 d range (obligation); the baseline 5 d is only an interior point'; ...
    'if compatible_with_observation_start_week were false it would be reported as-is (CL-03 failure handling)'};
mwspm_writejson(fullfile(ctx.dataDir, 'uncertainty_summary.json'), usum);

summary = struct();
summary.experiment_id = ctx.expId;
summary.task = char(string(getdef(cfg, 'task', '-')));
summary.parent_fit_run_id = fitRun.runId;
summary.parent_gate_run_id = gateRun.runId;
summary.deterministic = true;
summary.rng = 'none (deterministic post-processing + recomputation of deterministic trajectories)';
summary.metrics = metrics;
summary.supplements = supplements;
summary.observation_start_week = obsStart;
summary.main_report_point = fit.main_report_point;
summary.honesty = ['Interval is trivial by construction (rho-m0 ridge pins the detected trajectory); ' ...
    'OAT is not a joint region; gamma_h endpoints 4 d / 8 d covered as obligated.'];
mwspm_writejson(fullfile(ctx.dataDir, 'summary.json'), summary);

out = struct();
out.status = 'SUCCESS';
out.parent_run_id = fitRun.runId;
out.parent_gate_run_id = gateRun.runId;
out.deterministic = true;
out.rng = 'none';
out.metrics = metrics;
out.outputs = {'warning_week_vs_rho.csv', 'warning_week_oat.csv', 'uncertainty_summary.json', 'summary.json'};
out.note = ['T14.01. The warning-week interval is a parameter-uncertainty interval (value type 2); ' ...
    'the OAT weeks and the current-caliber week are deterministic points; the 47 weekly samples of ' ...
    'one trajectory are never treated as replicates.'];
end

% ---------------------------------------------------------------------------
function [tauStar, r0StartWeek, T0atTau, cumAtTau] = unc_tau_star(pre, D, th, rhoUsed, gammaH)
% Identical caliber to exp_c4_gate_g1's gate_tau_star, with gamma_h made
% explicit so the 4 d / 8 d endpoint scenarios can perturb it consistently in
% BOTH the trajectory (kernel) and R0.  The cross-check against
% gate_g1_weeks.csv asserts that with gammaH = pRef.gamma_h this reproduces
% the gate exactly.
nW = pre.n_weeks;
sp = struct('onset_doy', th.onset_doy, 'offset_doy', th.offset_doy, ...
    'kappa_days', pre.kappa_days, 'gamma', th.gamma, 'beta', th.beta, ...
    's_floor', pre.s_floor, 'min_window_days', pre.min_window_days, ...
    'normalize', true, 'norm_factor', NaN);
s = mwspm_season(pre.doy, sp);
trRaw = mwspm_thermal_traits(D.T_daily, 'raw');
muEff = th.mu_s * trRaw.mu;
R0daily = mwspm_r0(D.T_daily, th.m0 * s, muEff, gammaH);
wk = floor(((1:numel(R0daily))' - 1) / 7) + 1;
wk = min(max(wk, 1), nW);
gt1 = find(R0daily > 1, 1);
if isempty(gt1)
    r0StartWeek = NaN;
else
    r0StartWeek = wk(gt1);
end
R0w = accumarray(wk, R0daily, [nW 1], @mean, NaN);
preK = pre;
preK.gamma_h = gammaH;
kern = mwspm_calib_kernel(preK, th, struct());
cumDet = cumsum(kern.Yhat(:));
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
function th = unc_theta_from_fit(tf)
th = struct('m0', double(tf.m0), 'rho', double(tf.rho), 'onset_doy', double(tf.onset_doy), ...
    'offset_doy', double(tf.offset_doy), 'mu_s', double(tf.mu_s), 'gamma', double(tf.gamma), ...
    'beta', double(tf.beta), 'd_weeks', double(tf.d_weeks));
end

% ---------------------------------------------------------------------------
function P = unc_find_parent(root, wantId, expId, relPath)
% Latest results/runs/<run_id> whose manifest is SUCCESS for experiment <expId>
% and that carries <relPath>.  Returns [] when there is none (callers assert).
runsDir = fullfile(root, 'results', 'runs');
d = dir(runsDir);
P = [];
for k = 1:numel(d)
    if ~d(k).isdir || strcmp(d(k).name, '.') || strcmp(d(k).name, '..')
        continue
    end
    if ~isempty(wantId) && ischar(wantId) && ~strcmp(char(wantId), d(k).name)
        continue
    end
    f = fullfile(runsDir, d(k).name, relPath);
    if exist(f, 'file') ~= 2
        continue
    end
    manF = fullfile(runsDir, d(k).name, 'manifest.json');
    if exist(manF, 'file') == 2
        man = mwspm_readjson(manF);
        if isfield(man, 'status') && ~strcmp(char(string(man.status)), 'SUCCESS')
            continue
        end
        if isfield(man, 'experiment_id') && ~strcmp(char(string(man.experiment_id)), expId)
            continue
        end
    end
    if isempty(P) || strcmp(d(k).name, P.runId) > 0
        P = struct('runId', d(k).name, 'path', f);
    end
end
end

% ---------------------------------------------------------------------------
function prof = unc_read_profile(path)
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
function gate = unc_read_gate_weeks(path)
T = mwspm_readcsv(path);
gate = struct();
cols = {'branch_lag', 'rho', 'tau_star_week', 'T0_at_tau_star', 'cum_detected_at_tau_star', ...
        'R0_gt1_start_week', 'observation_start_week', 'current_caliber_100_week'};
for k = 1:numel(cols)
    gate.(cols{k}) = mwspm_csvcol(T, cols{k}, 'numeric');
end
end

% ---------------------------------------------------------------------------
function s = num2strOrNA(v)
if isfinite(v), s = sprintf('%g', v); else s = 'NA'; end
end
function s = numOrNA(v)
if isfinite(v), s = sprintf('%.17g', v); else s = 'NA'; end
end
function s = yesNo(c)
if c, s = 'YES'; else s = 'NO'; end
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
