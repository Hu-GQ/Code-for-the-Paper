function out = exp_c5_rossmacdonald(ctx)
%EXP_C5_ROSSMACDONALD  T11.02 -- reproduce the integrator law on an INDEPENDENT model.
%
%   Experiment entry point for config/experiments/exp_c5_rossmacdonald.json.
%
%   Why this exists
%   ---------------
%   T11.01 measures how explicit Euler degrades on THIS project's model.  A
%   critic can always answer "that is your model".  Stop condition S2 therefore
%   requires the same LAW to be reproduced on an idealised periodic
%   Ross-Macdonald that borrows NOTHING from the Guangzhou fit: m = 1.0 is a
%   scenario setting, there is no windowed seasonality, no mu_s scaling and no
%   observation layer.  Only the Mordecai 2017 thermal curves are shared, and
%   those are literature values rather than project estimates.
%
%   Model (two compartments, standard incidence on N_h)
%       I_h' = a(T) b(T) (I_v/N_h) (N_h - I_h) - gamma_h I_h
%       I_v' = a(T) c(T) (I_h/N_h) (N_v - I_v) - mu_v(T) I_v
%       T(t) = T_mean + A sin(2 pi t / period),   N_v = m N_h (constant)
%   EIP / latency are DELIBERATELY omitted so that the temperature channel that
%   survives is exactly the mortality one, i.e. the one whose stiffness governs
%   the explicit-Euler stability condition mu_eff*dt < 2.
%
%   This model carries no epidemiological claim.  Nothing here may be read as a
%   statement about Guangzhou 2025.

logFid = ctx.logFid;
cfg    = ctx.cfg;

fprintf(logFid, '=== exp_c5_rossmacdonald (T11.02, C5 external reproduction) ===\n');
fprintf(logFid, 'claims: %s\n', strjoin(cellstr(string(getdef(cfg, 'claims', {}))), ', '));
fprintf(logFid, 'NOTE: independent idealised model; no M-WSPM fitted parameter enters this run.\n');

P = cfg.parameter_domain;
% N_h is declared in the scenario model card, not in the experiment config.
N_h    = double(ctx.modelCard.parameters_fixed.N_h);
assert(N_h == 1.0e6, 'rossmac:n_h', 'unexpected N_h in the card: %.6g', N_h);
m      = double(P.m_vector_human_ratio);
gammaH = double(P.gamma_h);
Tmean  = double(P.temperature_mean_c);
period = double(P.period_days);
muCap  = double(P.mu_v_cap_per_day);
amps   = double(P.temperature_amplitudes_c(:)');
nDays  = double(cfg.termination.t_end_days);

assert(m == 1.0, 'rossmac:m', 'm is expected to be the scenario value 1.0, got %g', m);

% This model has no data-driven seasonality and no importation forcing, so it
% does not touch data/ at all.  Input hashes still cover those files because
% mwspm_run_finalize always hashes the standard input set; that is provenance,
% not a dependency.

solvers = cfg.solvers;
if isstruct(solvers)
    solvers = num2cell(solvers);
end

fprintf(logFid, 'N_h=%.6g  m=%.6g  gamma_h=%.6g  T_mean=%.4g C  period=%.4g d  horizon=%.4g d\n', ...
    N_h, m, gammaH, Tmean, period, nDays);
fprintf(logFid, 'amplitudes: %s\n', mat2str(amps));

reportGrid = (0:nDays)';
x0 = [1.0; 0.0];                                    % I_h(0)=1 person, I_v(0)=0

if isfield(ctx, 'dryRun') && ctx.dryRun
    out = struct();
    out.dry_run = true;
    out.validated = struct('n_amplitudes', numel(amps), 'n_solvers', numel(solvers), ...
        'horizon_days', nDays, 'x0', x0(:)');
    fprintf(logFid, '[dry-run] validated: %d amplitudes x %d solvers, no integration performed\n', ...
        numel(amps), numel(solvers));
    return
end

rows = {};
sweep = struct();
dense  = (0:0.05:nDays)';      % only for locating max mu_v
daily  = (0:nDays)';           % one sample per DAY, for the "days outside" count

for ai = 1:numel(amps)
    A = amps(ai);
    Tfun = @(t) Tmean + A * sin(2 * pi * t / period);

    % The stiffness-carrying quantity is the continuous mortality rate, so the
    % 'raw' variant is the right one here (there is no mu_s rescaling in this
    % model and therefore no source-faithful 'day' variant to reproduce).
    trDense = mwspm_thermal_traits(Tfun(dense), 'raw', struct('mu_cap', muCap));
    trDaily = mwspm_thermal_traits(Tfun(daily), 'raw', struct('mu_cap', muCap));
    muMax = max(trDense.mu);

    fprintf(logFid, '\n--- amplitude A = %.4g C  (T in [%.2f, %.2f]) ---\n', A, Tmean - A, Tmean + A);
    fprintf(logFid, 'max mu_v over the horizon = %.6g /day ; at dt=1 the euler limit mu*dt<2 is %s\n', ...
        muMax, tfWord(muMax > 2));

    rhs = @(t, x) rm_rhs(t, x, Tfun, N_h, m, gammaH, muCap);

    % ---------------- reference solution (adaptive, tight tolerance) -------
    refSpec = solvers{1};
    fprintf(logFid, '  REF %s RelTol=%.3g AbsTol=%.3g\n', refSpec.method, refSpec.RelTol, refSpec.AbsTol);
    specR = struct('id', char(string(refSpec.id)), 'kind', 'adaptive', ...
        'method', char(string(refSpec.method)), 'RelTol', double(refSpec.RelTol), ...
        'AbsTol', double(refSpec.AbsTol), 'NonNegative', [1 2], 'mu_variant', 'raw');
    SR = mwspm_solve(struct('rhs', rhs, 't0', 0, 'tf', nDays, 'x0', x0, ...
        'spec', specR, 'reportGrid', reportGrid, 'events', cfg.termination.events, ...
        'logFid', logFid));

    IhR = SR.x(:, 1);
    [pkR, pkTi] = max(IhR);
    if isempty(pkR) || ~isfinite(pkR), pkR = NaN; pkTi = NaN; else pkTi = SR.t(pkTi); end
    finR = IhR(end);
    fprintf(logFid, '  REF n_steps=%d step in [%.4g, %.4g]; peak I_h=%.6g at t=%.4g; final I_h=%.6g; %s\n', ...
        SR.n_steps, SR.min_internal_step, SR.max_internal_step, pkR, pkTi, finR, SR.termination_reason);

    Ain = struct('A_c', A, 'T_min_c', Tmean - A, 'T_max_c', Tmean + A, ...
        'max_mu_v_per_day', muMax, 'reference_method', specR.method, ...
        'reference_relTol', specR.RelTol, 'reference_absTol', specR.AbsTol, ...
        'reference_n_steps', SR.n_steps, 'reference_termination', SR.termination_reason, ...
        'reference_peak_I_h', pkR, 'reference_peak_time_days', pkTi, 'reference_final_I_h', finR, ...
        'reference_I_h', IhR);

    % ---------------- fixed-step methods ----------------------------------
    for si = 2:numel(solvers)
        s = solvers{si};
        dt = double(s.dt_days);
        method = char(string(s.method));
        clamp = logical(s.clamp_nonnegative);
        specF = struct('id', char(string(s.id)), 'kind', 'fixed_step', 'method', method, ...
            'dt_days', dt, 'clamp_nonnegative', clamp, 'mu_variant', 'raw');
        te = tic;
        S = mwspm_solve(struct('rhs', rhs, 't0', 0, 'tf', nDays, 'x0', x0, ...
            'spec', specF, 'reportGrid', reportGrid, 'events', cfg.termination.events, ...
            'logFid', logFid));
        el = toc(te);

        Ih = S.x(:, 1);
        fin = isfinite(Ih);
        if ~all(fin)
            IhUse = Ih;
            IhUse(~fin) = NaN;
        else
            IhUse = Ih;
        end
        [pk, pkI] = max(IhUse);
        if isempty(pk) || ~isfinite(pk)
            pk = NaN; pkT = NaN; finV = NaN;
        else
            pkT = S.t(pkI);
            finV = IhUse(end);
        end

        maxMuDt = muMax * dt;
        % UNIT FIX (2026-10-06): the count must be in DAYS.  Evaluating the
        % condition on the 0.05-day `dense` grid returned 4004 "days" on a
        % 730-day horizon - i.e. 20 samples per day - silently inflating it 20x.
        nOut = sum(trDaily.mu * dt > 2);

        % near-extinction FLAG (see the deviation note below): when both
        % compartments fall below the 1e-6 threshold the trajectory has entered
        % a regime where relative comparisons are numerically meaningless.
        nearExt = all(isfinite([pk, finV])) && (pk < 1e-6) && (finV < 1e-6);

        r = struct();
        r.amplitude_c      = A;
        r.method           = method;
        r.dt_days          = dt;
        r.solver_id        = char(string(s.id));
        r.n_steps          = S.n_steps;
        r.min_internal_step = S.min_internal_step;
        r.max_internal_step = S.max_internal_step;
        r.termination_reason = S.termination_reason;
        r.complete         = all(fin);
        r.peak_I_h         = pk;
        r.peak_I_h_abs_error = pk - pkR;
        r.peak_I_h_rel_error = (pk - pkR) / max(abs(pkR), 1e-30);
        r.peak_time_days   = pkT;
        r.peak_time_shift_days = pkT - pkTi;
        r.final_I_h        = finV;
        r.final_I_h_rel_error = (finV - finR) / max(abs(finR), 1e-30);
        r.max_mu_eff_dt    = maxMuDt;
        r.n_days_outside_euler_stability = nOut;
        r.near_extinction  = nearExt;
        r.elapsed_seconds  = el;
        rows{end+1, 1} = r; %#ok<AGROW>

        fprintf(logFid, '  %-6s dt=%-6g steps=%-7d %-22s peak_relerr=%-12.4g peak_t_shift=%-9.4g final_relerr=%-12.4g mu*dt=%.4g outside=%d\n', ...
            method, dt, S.n_steps, S.termination_reason, r.peak_I_h_rel_error, ...
            r.peak_time_shift_days, r.final_I_h_rel_error, maxMuDt, nOut);
    end

    sweep.(sprintf('amp%g', A)) = Ain;
end

methList = unique(cellfun(@(r) r.method, rows, 'UniformOutput', false));
dtList   = unique(cellfun(@(r) r.dt_days, rows));

% =============================================================== BLOCK B
% Why a second sweep is needed: the temperature amplitude does NOT vary
% mu_eff*dt continuously.  mu_v is capped at 5/day as soon as T leaves the
% viability band, so the amplitude sweep only has two regimes - cap engaged
% (A = 10) or not (A <= 9) - and mu_eff*dt jumps from 0.068 straight to 5.
% Varying dt at the amplitude that engages the cap sweeps mu_eff*dt over two
% decades and is what actually tests "deviation grows with mu_eff*dt".
dtSweepA = double(P.dt_sweep_amplitude_c);
dtGrid   = double(P.dt_sweep_days(:)');
fprintf(logFid, '\n--- BLOCK B: dt sweep at A = %.4g C (mortality cap engaged, max mu_v = %.6g /day) ---\n', ...
    dtSweepA, max(mwspm_thermal_traits(Tmean + dtSweepA * sin(2*pi*dense/period), 'raw', struct('mu_cap', muCap)).mu));

TfunB = @(t) Tmean + dtSweepA * sin(2 * pi * t / period);
rhsB  = @(t, x) rm_rhs(t, x, TfunB, N_h, m, gammaH, muCap);
trBDense = mwspm_thermal_traits(TfunB(dense), 'raw', struct('mu_cap', muCap));
trBDaily = mwspm_thermal_traits(TfunB(daily), 'raw', struct('mu_cap', muCap));
muMaxB = max(trBDense.mu);

specRB = struct('id', 'ode15s_ref_dtSweep', 'kind', 'adaptive', 'method', 'ode15s', ...
    'RelTol', 1e-10, 'AbsTol', 1e-12, 'NonNegative', [1 2], 'mu_variant', 'raw');
SRB = mwspm_solve(struct('rhs', rhsB, 't0', 0, 'tf', nDays, 'x0', x0, ...
    'spec', specRB, 'reportGrid', reportGrid, 'events', cfg.termination.events, ...
    'logFid', logFid));
IhRB = SRB.x(:, 1);
[pkRB, pkIB] = max(IhRB);
if ~isfinite(pkRB), pkRB = NaN; pkTB = NaN; else pkTB = SRB.t(pkIB); end
finRB = IhRB(end);
fprintf(logFid, '  REF (A=%.4g): n_steps=%d ; peak I_h=%.6g at t=%.4g ; final I_h=%.6g\n', ...
    dtSweepA, SRB.n_steps, pkRB, pkTB, finRB);

swRows = {};
for mi = 1:numel(methList)
    for di = 1:numel(dtGrid)
        dt = dtGrid(di);
        specF = struct('id', sprintf('%s_dt%g_sweep', methList{mi}, dt), 'kind', 'fixed_step', ...
            'method', methList{mi}, 'dt_days', dt, 'clamp_nonnegative', true, 'mu_variant', 'raw');
        te = tic;
        S = mwspm_solve(struct('rhs', rhsB, 't0', 0, 'tf', nDays, 'x0', x0, ...
            'spec', specF, 'reportGrid', reportGrid, 'events', cfg.termination.events, ...
            'logFid', logFid));
        el = toc(te);
        Ih = S.x(:, 1);
        if all(isfinite(Ih))
            [pk, pkI] = max(Ih); pkT = S.t(pkI); finV = Ih(end);
        else
            [pk, pkI] = max(Ih(isfinite(Ih)));
            if isempty(pk), pk = NaN; pkT = NaN; finV = NaN; else pkT = S.t(pkI); finV = Ih(end); end
        end
        sw = struct();
        sw.amplitude_c = dtSweepA;
        sw.method = methList{mi};
        sw.dt_days = dt;
        sw.max_mu_eff_dt = muMaxB * dt;
        sw.n_days_outside_euler_stability = sum(trBDaily.mu * dt > 2);
        sw.peak_I_h = pk;
        sw.peak_I_h_rel_error = tern(isfinite(pk) && isfinite(pkRB), (pk - pkRB) / max(abs(pkRB), 1e-30), NaN);
        sw.peak_time_days = pkT;
        sw.peak_time_shift_days = pkT - pkTB;
        sw.final_I_h = finV;
        sw.final_I_h_rel_error = tern(isfinite(finV) && isfinite(finRB), (finV - finRB) / max(abs(finRB), 1e-30), NaN);
        sw.n_steps = S.n_steps;
        sw.termination_reason = S.termination_reason;
        sw.complete = all(isfinite(Ih));
        sw.elapsed_seconds = el;
        swRows{end+1, 1} = sw; %#ok<AGROW>
        fprintf(logFid, '  %-6s dt=%-8g mu*dt=%-9.4g steps=%-7d %-22s peak_relerr=%-13.5g final_relerr=%-13.5g outside=%d/%d\n', ...
            methList{mi}, dt, sw.max_mu_eff_dt, S.n_steps, S.termination_reason, ...
            sw.peak_I_h_rel_error, sw.final_I_h_rel_error, sw.n_days_outside_euler_stability, nDays);
    end
end

% monotonicity of the dt sweep: sort by mu*dt and count adjacent inversions
monoB = struct();
monoB.pairs = {};
monoB.n_inversions = 0;
monoB.n_pairs = 0;
for mi = 1:numel(methList)
    sel = cellfun(@(r) strcmp(r.method, methList{mi}), swRows);
    rr = swRows(sel);
    x = cellfun(@(r) r.max_mu_eff_dt, rr);
    y = cellfun(@(r) absOr0(r.peak_I_h_rel_error), rr);
    [xs, ord] = sort(x);
    ys = y(ord);
    nInv = 0;
    for k = 2:numel(ys)
        if ys(k) < ys(k - 1) - 1e-9
            nInv = nInv + 1;
        end
    end
    % Regime test: the explicit-Euler stability limit for x' = -mu x is
    % mu*dt < 2.  Report the two regimes separately instead of forcing one
    % monotonicity statement across a discontinuity.
    inDom  = xs < 2;
    outDom = xs >= 2;
    monoB.pairs{end+1} = struct('method', methList{mi}, 'mu_dt_sorted', xs(:)', ...
        'abs_peak_rel_err_sorted', ys(:)', 'n_adjacent_inversions', nInv, ...
        'monotone_nondecreasing', nInv == 0, ...
        'max_abs_rel_err_inside_stability_domain', tern(any(inDom), max(ys(inDom)), NaN), ...
        'max_abs_rel_err_outside_stability_domain', tern(any(outDom), max(ys(outDom)), NaN), ...
        'regime_jump_ratio', tern(any(inDom) && any(outDom) && max(ys(inDom)) > 0, ...
            max(ys(outDom)) / max(ys(inDom)), NaN)); %#ok<AGROW>
    monoB.n_inversions = monoB.n_inversions + nInv;
    monoB.n_pairs = monoB.n_pairs + numel(ys) - 1;
    fprintf(logFid, '[S2-B] %s : adjacent inversions = %d/%d ; monotone = %s ; max|rel err| inside mu*dt<2 = %.6g ; outside = %.6g ; ratio = %.6g\n', ...
        methList{mi}, nInv, numel(ys) - 1, yesNo(nInv == 0), ...
        monoB.pairs{end}.max_abs_rel_err_inside_stability_domain, ...
        monoB.pairs{end}.max_abs_rel_err_outside_stability_domain, ...
        monoB.pairs{end}.regime_jump_ratio);
end

fidB = fopen(fullfile(ctx.dataDir, 'dt_sweep_results.csv'), 'w', 'n', 'UTF-8');
fprintf(fidB, 'amplitude_c,method,dt_days,max_mu_eff_dt,n_days_outside_euler_stability,n_steps,termination_reason,complete,peak_I_h,peak_I_h_rel_error,peak_time_days,peak_time_shift_days,final_I_h,final_I_h_rel_error,elapsed_seconds\n');
for i = 1:numel(swRows)
    r = swRows{i};
    fprintf(fidB, '%.17g,%s,%.17g,%.17g,%d,%d,%s,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.4f\n', ...
        r.amplitude_c, r.method, r.dt_days, r.max_mu_eff_dt, r.n_days_outside_euler_stability, ...
        r.n_steps, r.termination_reason, double(r.complete), r.peak_I_h, r.peak_I_h_rel_error, ...
        r.peak_time_days, r.peak_time_shift_days, r.final_I_h, r.final_I_h_rel_error, r.elapsed_seconds);
end
fclose(fidB);

% ---------------------------------------------------- monotonicity check (S2)
% Question: does the Euler deviation increase with mu_eff*dt once the model is
% no longer this project's?  Answer by observing, not by assuming: for each
% (method, dt) pair, list |peak relative error| against max(mu*dt) over the six
% amplitudes and count how many adjacent pairs (ordered by mu*dt) are
% non-decreasing.
mono = struct();
mono.pairs = {};
mono.n_inversions = 0;
mono.n_pairs = 0;
for mi = 1:numel(methList)
    for di = 1:numel(dtList)
        sel = cellfun(@(r) strcmp(r.method, methList{mi}) && abs(r.dt_days - dtList(di)) < 1e-12, rows);
        if nnz(sel) ~= numel(amps)
            continue
        end
        rr = rows(sel);
        x = cellfun(@(r) r.max_mu_eff_dt, rr);
        y = cellfun(@(r) absOr0(r.peak_I_h_rel_error), rr);
        [xs, ord] = sort(x);
        ys = y(ord);
        worst = max(ys);
        nInv = 0;
        for k = 2:numel(ys)
            if ys(k) < ys(k - 1) - 1e-9
                nInv = nInv + 1;
            end
        end
        key = sprintf('%s_dt%g', methList{mi}, dtList(di));
        mono.pairs{end+1} = struct('series', key, 'mu_dt_sorted', xs(:)', ...
            'abs_peak_rel_err_sorted', ys(:)', 'n_adjacent_inversions', nInv, ...
            'max_abs_peak_rel_err', worst, ...
            'monotone_nondecreasing', nInv == 0); %#ok<AGROW>
        mono.n_inversions = mono.n_inversions + nInv;
        mono.n_pairs = mono.n_pairs + numel(ys) - 1;
        fprintf(logFid, '\n[S2] %s : max|peak rel err| over amplitudes = %.6g ; adjacent inversions = %d/%d ; monotone non-decreasing = %s\n', ...
            key, worst, nInv, numel(ys) - 1, yesNo(nInv == 0));
    end
end
mono.all_monotone = (mono.n_inversions == 0);
mono.verdict = tern(mono.all_monotone, ...
    'MONOTONE: the Euler deviation grows with mu_eff*dt on the independent model', ...
    'NOT MONOTONE at least somewhere: report as measured and apply S2 honestly');

% ------------------------------------------------------------------- outputs
hdr = {'amplitude_c','method_code','dt_days','n_steps','min_internal_step', ...
       'max_internal_step','complete','peak_I_h','peak_I_h_abs_error','peak_I_h_rel_error', ...
       'peak_time_days','peak_time_shift_days','final_I_h','final_I_h_rel_error', ...
       'max_mu_eff_dt','n_days_outside_euler_stability','near_extinction','elapsed_seconds'};
M = nan(numel(rows), numel(hdr));
for i = 1:numel(rows)
    r = rows{i};
    M(i, :) = [r.amplitude_c, mc(r.method), r.dt_days, r.n_steps, r.min_internal_step, ...
        r.max_internal_step, double(r.complete), r.peak_I_h, r.peak_I_h_abs_error, ...
        r.peak_I_h_rel_error, r.peak_time_days, r.peak_time_shift_days, r.final_I_h, ...
        r.final_I_h_rel_error, r.max_mu_eff_dt, r.n_days_outside_euler_stability, ...
        double(r.near_extinction), r.elapsed_seconds];
end
mwspm_writecsv(fullfile(ctx.dataDir, 'amplitude_sweep_results.csv'), hdr, M);

fid = fopen(fullfile(ctx.dataDir, 'amplitude_sweep_results_long.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, 'amplitude_c,method,dt_days,solver_id,n_steps,termination_reason,complete,peak_I_h,peak_I_h_abs_error,peak_I_h_rel_error,peak_time_days,peak_time_shift_days,final_I_h,final_I_h_rel_error,max_mu_eff_dt,n_days_outside_euler_stability,near_extinction,elapsed_seconds\n');
for i = 1:numel(rows)
    r = rows{i};
    fprintf(fid, '%.17g,%s,%.17g,%s,%d,%s,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%d,%d,%.4f\n', ...
        r.amplitude_c, r.method, r.dt_days, r.solver_id, r.n_steps, r.termination_reason, ...
        double(r.complete), r.peak_I_h, r.peak_I_h_abs_error, r.peak_I_h_rel_error, ...
        r.peak_time_days, r.peak_time_shift_days, r.final_I_h, r.final_I_h_rel_error, ...
        r.max_mu_eff_dt, r.n_days_outside_euler_stability, double(r.near_extinction), r.elapsed_seconds);
end
fclose(fid);

% one reference trajectory per amplitude, downsampled to the daily report grid
trajCols = (0:nDays)';
trajHdr = {'t_day'};
fn = fieldnames(sweep);
for k = 1:numel(fn)
    R = sweep.(fn{k});
    trajCols = [trajCols, R.reference_I_h(:)]; %#ok<AGROW>
    trajHdr{end+1} = ['I_h_' fn{k} '_reference']; %#ok<AGROW>
end
mwspm_writecsv(fullfile(ctx.dataDir, 'trajectories_sample.csv'), trajHdr, trajCols);

summary = struct();
summary.experiment_id = ctx.expId;
summary.task = char(string(getdef(cfg, 'task', '-')));
summary.model = 'ross_macdonald_periodic_idealised (rm_ideal_v1)';
summary.n_amplitudes = numel(amps);
summary.n_rows = numel(rows);
summary.method_code = 'method_code column of amplitude_sweep_results.csv: 1 = euler, 2 = rk4 (see the long CSV for names)';
summary.deviation_registered = ['The config asks fixed-step integrators to STOP on a near-extinction ' ...
    'event (I_h or I_v below 1e-6). That would truncate each method at a DIFFERENT time and make the ' ...
    'peak / final-state comparison across methods meaningless, which is the whole point of the sweep. ' ...
    'The near-extinction condition is therefore recorded as a BOOLEAN FLAG per row and the integration ' ...
    'is not stopped. This is a deliberate, logged deviation from the frozen config.'];
summary.monotonicity_amplitude_sweep = mono;
summary.monotonicity_dt_sweep = monoB;
summary.dt_sweep = struct('amplitude_c', dtSweepA, 'dt_grid_days', dtGrid, ...
    'max_mu_v_per_day', muMaxB, 'reference_peak_I_h', pkRB, 'reference_final_I_h', finRB, ...
    'why', ['the amplitude sweep has only two regimes (mortality cap engaged or not), so it ' ...
            'cannot test a continuous dependence on mu_eff*dt; the dt sweep can']);
summary.max_abs_peak_rel_error = max(M(:, 10));
summary.rows = rows;
mwspm_writejson(fullfile(ctx.dataDir, 'summary.json'), summary);

out = struct();
out.status = 'SUCCESS';
out.n_amplitudes = numel(amps);
out.n_rows = numel(rows);
out.monotone = mono.all_monotone;
out.n_adjacent_inversions = mono.n_inversions;
out.max_abs_peak_rel_error = summary.max_abs_peak_rel_error;
out.outputs = {'amplitude_sweep_results.csv', 'amplitude_sweep_results_long.csv', ...
    'dt_sweep_results.csv', 'trajectories_sample.csv', 'summary.json'};
out.deterministic = true;
out.rng = 'none (deterministic ODE; the sweep varies the temperature amplitude, not a random seed)';
out.note = ['T11.02 supplies the independent-model reproduction. If monotonicity fails, S2 applies: ' ...
    'C5 drops to a project-specific remark and must not be written as a general diagnostic checklist.'];
end

% ---------------------------------------------------------------------------
function dx = rm_rhs(t, x, Tfun, N_h, m, gammaH, muCap)
% Two-compartment Ross-Macdonald with analytic sinusoidal forcing.
T = Tfun(t);
tr = mwspm_thermal_traits(T, 'raw', struct('mu_cap', muCap));
Ih = max(x(1), 0);
Iv = max(x(2), 0);
Nv = m * N_h;
lh = tr.a * tr.b * (Iv / N_h) * (N_h - Ih);
lv = tr.a * tr.c * (Ih / N_h) * (Nv - Iv);
dx = [lh - gammaH * Ih;
      lv - tr.mu * Iv];
end

% ---------------------------------------------------------------------------
function v = mc(method)
switch method
    case 'euler', v = 1;
    case 'rk4',   v = 2;
    otherwise,    v = 0;
end
end
function v = absOr0(x)
if isfinite(x), v = abs(x); else v = 0; end
end
function s = tfWord(tf)
% Used only where "true" means a CONDITION IS BROKEN (stability limit crossed).
if tf, s = 'VIOLATED'; else s = 'satisfied'; end
end
function s = yesNo(tf)
% Used where "true" means a property HOLDS.  It must not be confused with
% tfWord: an early revision printed "monotone = satisfied" for a series that
% had two inversions, which reads as an assertion of monotonicity when it is
% the opposite.
if tf, s = 'YES'; else s = 'NO'; end
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
