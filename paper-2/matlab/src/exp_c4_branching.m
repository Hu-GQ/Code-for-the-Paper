function out = exp_c4_branching(ctx)
%EXP_C4_BRANCHING  the analysis -- stochastic two-type branching process + rho.
%
%   Experiment entry point for config/experiments/exp_c4_branching.json.
%   Unlocked by gate G1 (the analysis) = GO, 2026-10-06.
%
%   WHAT IS COMPUTED
%   ----------------
%   A time-inhomogeneous 4-type branching process (E_h, I_h, E_v, I_v) that is
%   the linearisation of M-WSPM gen3_rhohat about the disease-free state, with
%   rates taken from the FITTED daily forcing and importation set to ZERO
%   (a single introduction cannot be studied while Lambda(t) keeps seeding).
%
%   Estimands, kept strictly apart:
%     (a) p_surv(tau)        P(lineage still alive at t_end | 1 imported I_h).
%                            rho-independent.  EXACT (backward ODE).
%     (b) p_est(tau; rho)    P(the lineage produces at least one DETECTED
%                            autochthonous human infection inside the window).
%                            EXACT (backward ODE) and SIMULATED (N_sim = 1e4).
%     (c) p_carmona(tau)     Carmona & Gandon 2020 Eq.(2) on a single-type
%                            birth-death reduction matched on (R0, r).
%                            CONTROL ONLY - cited, never claimed as ours.
%     (d) *_static           the same with the rates FROZEN at tau.  This is
%                            the quasi-static approximation whose validity is
%                            obligation O8; the gap to the time-varying value
%                            is measured here rather than argued about.
%
%   THREE KINDS OF QUANTITY, NEVER MIXED
%   ------------------------------------
%     deterministic point       : (a)-(d) evaluated once, no interval;
%     parameter uncertainty     : NOT produced here - rho is a GRID, not an
%                                 estimated interval (the Wilks band from
%                                 the analysis covers the whole [0.02, 1.00]);
%     replicate distribution    : the empirical distribution of the first
%                                 detection week over N_sim replicates, with
%                                 the Monte-Carlo standard error reported.
%   Weekly time samples of one epidemic are NEVER treated as replicates.

logFid = ctx.logFid;
root   = ctx.root;
cfg    = ctx.cfg;
modelCard = ctx.modelCard;

fprintf(logFid, '=== exp_c4_branching (the analysis, C4 stochastic control) ===\n');
fprintf(logFid, 'claims: %s\n', strjoin(cellstr(string(getdef(cfg, 'claims', {}))), ', '));
fprintf(logFid, 'unlocked by gate G1 (the analysis) = GO, 2026-10-06.\n');

D = mwspm_load_data(root, modelCard);
pBase = mwspm_params(modelCard);

PD = cfg.parameter_domain;
rhoGrid = double(PD.rho_grid.values(:)');
tauAll  = double(PD.tau_weeks.from) : double(PD.tau_weeks.step) : double(PD.tau_weeks.to);
tauAll  = tauAll(:);
tauSub  = double(PD.tau_sim_subset_weeks(:));
baseSeed = double(cfg.rng.base_seed);
sv = findSolver(cfg, 'stoch_branching');
nSim    = double(getdef(sv, 'n_sim', 10000));
pTarget = double(getdef(sv, 'p_target', 0.05));
variants = cellstr(string(getdef(cfg, 'variants', {'day', 'raw'})));

% ------------------------------------------------------------- parents ------
fitRun = find_parent(root, getdef(cfg, 'parent_run_id', []), 'exp_c4_rho_fit', 'data/fit_best.json');
assert(~isempty(fitRun), 'branching:noParent', ...
    'no SUCCESS run of exp_c4_rho_fit found; the branching rates must come from a real fit.');
fit = mwspm_readjson(fitRun.path);
th = struct();
tf = fit.main_report_point.theta;
th.m0 = double(tf.m0); th.rho = double(tf.rho);
th.onset_doy = double(tf.onset_doy); th.offset_doy = double(tf.offset_doy);
th.mu_s = double(tf.mu_s); th.gamma = double(tf.gamma); th.beta = double(tf.beta);
th.d_weeks = double(tf.d_weeks);
fprintf(logFid, 'parent fit run : %s (lag %d)\n', fitRun.runId, th.d_weeks);
fprintf(logFid, 'theta          : m0=%.6f rho=%.6f mu_s=%.6f onset=%.3f offset=%.3f\n', ...
    th.m0, th.rho, th.mu_s, th.onset_doy, th.offset_doy);

% Observation start week is recomputed HERE from the observations rather than
% read back from the gate run: the definition is "first week with a non-zero
% observed local count" and the data are already loaded, so a second-hand copy
% would only add a failure mode (the gate verdict JSON does not carry it).
Yo = D.Y_weekly(:);
obsStart = find(Yo > 0, 1);
if isempty(obsStart)
    obsStart = NaN;
end

gateRun = find_parent(root, [], 'exp_c4_gate_g1', 'data/gate_g1_verdict.json');
tauStar = NaN;
tauStarMax = NaN;
tauStarByBranch = struct();
if ~isempty(gateRun)
    gv = mwspm_readjson(gateRun.path);
    if isfield(gv, 'branches')
        bs = fieldnames(gv.branches);
        mn = nan(numel(bs), 1);
        mx = nan(numel(bs), 1);
        for b = 1:numel(bs)
            bb = gv.branches.(bs{b});
            mn(b) = double(bb.tau_star_min);
            mx(b) = double(bb.tau_star_max);
            tauStarByBranch.(bs{b}) = struct('min', mn(b), 'max', mx(b));
        end
        tauStar = min(mn);
        tauStarMax = max(mx);
    end
    fprintf(logFid, 'gate G1 run    : %s ; tau* over branches in [%s, %s]\n', ...
        gateRun.runId, num2strOrNA(tauStar), num2strOrNA(tauStarMax));
else
    fprintf(logFid, 'gate G1 run    : NOT FOUND (tau* comparison columns will be NA)\n');
end
fprintf(logFid, 'observation start week (recomputed from the observed series) = %s\n', ...
    num2strOrNA(obsStart));

if isfield(ctx, 'dryRun') && ctx.dryRun
    out = struct();
    out.dry_run = true;
    out.validated = struct('n_tau', numel(tauAll), 'n_rho', numel(rhoGrid), ...
        'n_sim', nSim, 'tau_subset', tauSub(:)', 'variants', {variants}, ...
        'parent_run_id', fitRun.runId);
    fprintf(logFid, '[dry-run] grid and parents readable; nothing simulated\n');
    return
end

tEval = 7 * (tauAll - 1);

% =========================================================== 1. analytic ====
fprintf(logFid, '\n--- 1. exact backward solve (deterministic, no interval) ---\n');
A = struct();
for vi = 1:numel(variants)
    vv = lower(char(variants{vi}));
    R = mwspm_branching_rates(D, pBase, th, vv);
    fprintf(logFid, 'variant %s : mu in [%.4g, %.4g], R0 in [%.4g, %.4g], r in [%.4g, %.4g], %d/%d days with R0>1\n', ...
        vv, min(R.mu_daily), max(R.mu_daily), min(R.R0_daily), max(R.R0_daily), ...
        min(R.r_daily), max(R.r_daily), nnz(R.R0_daily > 1), R.n_days);

    bw = mwspm_branching_backward(R, struct('rho', rhoGrid, 't_eval', tEval));
    fprintf(logFid, '  time-varying : %d solve(s), max |violation of [0,1]| = %.3e\n', ...
        bw.n_solves, bw.max_violation_of_unit_interval);
    bwS = mwspm_branching_backward(R, struct('rho', rhoGrid, 't_eval', tEval, ...
        'static', true, 'max_step', 5));
    fprintf(logFid, '  quasi-static : %d solve(s), max |violation of [0,1]| = %.3e\n', ...
        bwS.n_solves, bwS.max_violation_of_unit_interval);

    pc = carmona_control(R, tEval);

    A.(vv) = struct('R', R, 'bw', bw, 'bwS', bwS, 'pc', pc);
end

% ========================================================= 2. stochastic ====
fprintf(logFid, '\n--- 2. stochastic replicates (N_sim = %d per grid point, BOTH variants) ---\n', nSim);
% The config variants_note says day and raw must BOTH be reported ("不得只报
% 其中一个"); the tau-leap subset therefore validates the exact backward solve
% on each variant separately, using the SAME seed rule per variant (registered
% in raw_choices.csv).
simRows = {};
simDiag = struct('n_points', 0, 'max_abs_diff_in_mc_se', 0, 'sum_overshoot', 0, 'sum_cap', 0, ...
    'by_variant', struct());

for vi = 1:numel(variants)
    vv = lower(char(variants{vi}));
    Rsim = A.(vv).R;
    bwy = A.(vv).bw;
    vDiag = struct('n_points', 0, 'max_abs_diff_in_mc_se', 0, 'sum_overshoot', 0, 'sum_cap', 0);
    fprintf(logFid, '  variant %<local path>, vv);
    for ti = 1:numel(tauSub)
        tw = tauSub(ti);
        for ri = 1:numel(rhoGrid)
            rr = rhoGrid(ri);
            seed = baseSeed + 1000 * round(rr * 100) + tw;
            s = mwspm_branching_sim(Rsim, tw, rr, nSim, struct('seed', seed, 'p_target', pTarget));
            iTau = find(tauAll == tw, 1);
            pAna = 1 - bwy.u(iTau, 2, ri);
            diff = s.p_est - pAna;
            inSe = tern(s.mc_se > 0, abs(diff) / s.mc_se, NaN);
            simRows{end + 1, 1} = {vv, tw, rr, nSim, s.p_est, s.mc_se, pAna, diff, inSe, ...
                s.n_extinct, s.n_alive_at_end, s.n_cap_forced, ...
                s.det_week_q05, s.det_week_q25, s.det_week_q50, s.det_week_q75, s.det_week_q95, ...
                s.n_det_week_observed, s.n_overshoot, seed, pTarget, s.elapsed_seconds}; %#ok<AGROW>
            vDiag.n_points = vDiag.n_points + 1;
            vDiag.sum_overshoot = vDiag.sum_overshoot + s.n_overshoot;
            vDiag.sum_cap = vDiag.sum_cap + s.n_cap_forced;
            if isfinite(inSe)
                vDiag.max_abs_diff_in_mc_se = max(vDiag.max_abs_diff_in_mc_se, inSe);
            end
            fprintf(logFid, '    tau=%2d rho=%.2f : p_sim=%.4f (se %.4f) vs p_exact=%.4f | diff=%+.4f = %.2f MC se | ext=%d cap=%d | detwk med=%s | %.1f s\n', ...
                tw, rr, s.p_est, s.mc_se, pAna, diff, inSe, s.n_extinct, s.n_cap_forced, ...
                num2strOrNA(s.det_week_q50), s.elapsed_seconds);
        end
    end
    simDiag.by_variant.(vv) = vDiag;
    simDiag.n_points = simDiag.n_points + vDiag.n_points;
    simDiag.sum_overshoot = simDiag.sum_overshoot + vDiag.sum_overshoot;
    simDiag.sum_cap = simDiag.sum_cap + vDiag.sum_cap;
    simDiag.max_abs_diff_in_mc_se = max(simDiag.max_abs_diff_in_mc_se, vDiag.max_abs_diff_in_mc_se);
end

% ------------------------------------------- 2b. discretisation refinement --
fprintf(logFid, '\n--- 2b. tau-leap discretisation refinement (p_target/4, BOTH variants) ---\n');
refRows = {};
refTau = double(getdef(cfg, 'dt_refinement_tau_weeks', [19, 21, 33]));
refRho = double(getdef(cfg, 'dt_refinement_rho', [0.13, 0.80]));
for vi = 1:numel(variants)
    vv = lower(char(variants{vi}));
    Rsim = A.(vv).R;
    bwy = A.(vv).bw;
    for ti = 1:numel(refTau)
        tw = refTau(ti);
        for ri = 1:numel(refRho)
            rr = refRho(ri);
            iTau = find(tauAll == tw, 1);
            jRho = find(rhoGrid == rr, 1);
            pAna = 1 - bwy.u(iTau, 2, jRho);
            seedC = baseSeed + 1000 * round(rr * 100) + tw;
            sC = mwspm_branching_sim(Rsim, tw, rr, nSim, struct('seed', seedC, 'p_target', pTarget));
            seedF = seedC + 777000;
            sF = mwspm_branching_sim(Rsim, tw, rr, nSim, struct('seed', seedF, 'p_target', pTarget / 4));
            refRows{end + 1, 1} = {vv, tw, rr, pAna, sC.p_est, sC.mc_se, sF.p_est, sF.mc_se, ...
                sF.p_est - sC.p_est, sF.n_overshoot, sC.n_overshoot, ...
                sC.elapsed_seconds, sF.elapsed_seconds}; %#ok<AGROW>
            fprintf(logFid, '  %s tau=%2d rho=%.2f : exact=%.4f | coarse=%.4f | fine=%.4f | fine-coarse=%+.4f | overshoot %d -> %d\n', ...
                vv, tw, rr, pAna, sC.p_est, sF.p_est, sF.p_est - sC.p_est, sC.n_overshoot, sF.n_overshoot);
        end
    end
end

% ===================================================== 3. profile shape =====
fprintf(logFid, '\n--- 3. establishment profile shape (the statement test) ---\n');
profRows = {};
for vi = 1:numel(variants)
    vv = lower(char(variants{vi}));
    bw = A.(vv).bw;
    for ri = 1:numel(rhoGrid)
        rr = rhoGrid(ri);
        p = 1 - bw.u(:, 2, ri);
        [iMax, lm] = shapeOf(p);
        single = numel(lm) <= 1;
        pAtTau = atWeek(p, tauStar);
        pAtObs = atWeek(p, obsStart);
        profRows{end + 1, 1} = {vv, rr, tauAll(iMax), max(p), numel(lm), double(single), ...
            strjoinLocal(lm), tauStar, pAtTau, obsStart, pAtObs, tauAll(iMax) - obsStart}; %#ok<AGROW>
        fprintf(logFid, '  %s rho=%.2f : peak week %d (p=%.4f), interior local maxima = %d %s, single-peaked = %s; p(tau*=%s)=%s, p(obs start=%s)=%s\n', ...
            vv, rr, tauAll(iMax), max(p), numel(lm), bracket(lm), yesNo(single), ...
            num2strOrNA(tauStar), num2strOrNA(pAtTau), num2strOrNA(obsStart), num2strOrNA(pAtObs));
    end
    ps = 1 - bw.q(:, 2);
    [iMaxS, lmS] = shapeOf(ps);
    profRows{end + 1, 1} = {vv, NaN, tauAll(iMaxS), max(ps), numel(lmS), double(numel(lmS) <= 1), ...
        strjoinLocal(lmS), tauStar, atWeek(ps, tauStar), ...
        obsStart, atWeek(ps, obsStart), tauAll(iMaxS) - obsStart}; %#ok<AGROW>
    fprintf(logFid, '  %s p_surv   : peak week %d (p=%.6f), interior local maxima = %d %s\n', ...
        vv, tauAll(iMaxS), max(ps), numel(lmS), bracket(lmS));
end

% =========================================================== 4. outputs =====
writeAnalyticCsv(fullfile(ctx.dataDir, 'branching_analytic_grid.csv'), A, variants, tauAll, rhoGrid);

hdr = {'variant', 'tau_week', 'rho', 'n_sim', 'p_est_sim', 'mc_se', 'p_est_exact', 'diff', ...
    'abs_diff_in_mc_se', 'n_extinct', 'n_alive_at_end', 'n_cap_forced', ...
    'det_week_q05', 'det_week_q25', 'det_week_q50', 'det_week_q75', 'det_week_q95', ...
    'n_det_week_observed', 'n_overshoot', 'seed', 'p_target', 'elapsed_seconds'};
fid = fopen(fullfile(ctx.dataDir, 'branching_sim_subset.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '%s\n', strjoin(hdr, ','));
for i = 1:numel(simRows)
    rw = simRows{i};
    fprintf(fid, '%s', rw{1});
    for c = 2:numel(hdr)
        fprintf(fid, ',%.17g', double(rw{c}));
    end
    fprintf(fid, '\n');
end
fclose(fid);

hdr2 = {'variant', 'tau_week', 'rho', 'p_est_exact', 'p_est_coarse', 'mc_se_coarse', 'p_est_fine', ...
    'mc_se_fine', 'fine_minus_coarse', 'n_overshoot_fine', 'n_overshoot_coarse', ...
    'elapsed_coarse_s', 'elapsed_fine_s'};
fid = fopen(fullfile(ctx.dataDir, 'branching_dt_refinement.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '%s\n', strjoin(hdr2, ','));
for i = 1:numel(refRows)
    rw = refRows{i};
    fprintf(fid, '%s', rw{1});
    for c = 2:numel(hdr2)
        fprintf(fid, ',%.17g', double(rw{c}));
    end
    fprintf(fid, '\n');
end
fclose(fid);

writeProfileCsv(fullfile(ctx.dataDir, 'p_establishment_by_week.csv'), profRows);

% ---- declared modelling choices (a file, so that they cannot be forgotten) --
ch = {'seed_type', 'I_h (one imported infectious human); the seed itself is not counted as detected';
    'importation', 'Lambda(t) = 0 during the branching process; the fitted model seeds through Lambda and a single introduction is impossible while it runs';
    'linearisation', 'S_h = N_h and S_v = N_v(t) = m0*N_h*s(t); the branching process is the DFE linearisation, not the full nonlinear model';
    'detection_rule', 'each newly created human infection is detected independently with probability rho';
    'estimand_primary', 'p_est = P(at least one DETECTED autochthonous human infection inside the window)';
    'estimand_secondary', 'p_surv = P(lineage alive at t_end); EXACT only, not simulated (surviving lineages grow without bound)';
    'estimand_control', 'p_carmona = single-type birth-death reduction matched on (R0, r), Carmona & Gandon 2020 Eq.(2), horizon truncated at t_end';
    'estimand_static', 'rates frozen at the introduction day; the gap to the time-varying value is the measured size of obligation O8';
    'weekly_R0_aggregation', 'MEAN of the daily values inside the week (declared; daily R0 is sawtoothed)';
    'forcing_sampling', 'zero-order hold on the daily grid, identical to the fitted model';
    'stochastic_scheme', 'tau-leap with fixed sub-steps inside each day; sub-step count chosen so that max_rate*dt <= p_target';
    'rng', sprintf('one mt19937ar RandStream per (variant, tau, rho), seed = %d + 1000*round(100*rho) + tau; no time-based seed; the SAME seed rule is used for both variants, so the day and raw streams coincide at equal (tau,rho) - registered here because the config seed policy has no variant term', baseSeed);
    'poisson_sampler', 'mwspm_poissrnd: Knuth product-of-uniforms for lambda <= 25, normal approximation above (declared, counted)';
    'clamp_policy', 'the tau-leap draws competing exits independently, so their sum can exceed the compartment; clamps are COUNTED (n_overshoot)';
    'cap_policy', sprintf('live population > %d or cumulative human infections > %d stops the replicate and declares it established; (1-rho)^cap is < 1e-24 on this rho grid', 400, 5000);
    'what_is_not_here', 'no parameter-uncertainty interval is produced: rho is a grid, and the the analysis Wilks band for rho already covers the whole [0.02, 1.00]';
    'variants_both_reported', 'the tau-leap subset and the dt refinement run on BOTH mu variants (day and raw), because the config variants_note forbids reporting only one; the exact backward solve covers both as well'};
fid = fopen(fullfile(ctx.dataDir, 'raw_choices.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, 'choice,value\n');
for i = 1:size(ch, 1)
    fprintf(fid, '%s,"%s"\n', ch{i, 1}, ch{i, 2});
end
fclose(fid);

% ------------------------------------------------------------------ summary -
summary = struct();
summary.experiment_id = ctx.expId;
summary.task = char(string(getdef(cfg, 'task', '-')));
summary.parent_fit_run_id = fitRun.runId;
summary.gate_g1_run_id = tern(isempty(gateRun), '', gateRun.runId);
summary.main_report_point = fit.main_report_point;
summary.deterministic = false;
summary.rng = struct('generator', 'mt19937ar', 'base_seed', baseSeed, ...
    'seed_rule', 'seed = base_seed + 1000*round(100*rho) + tau; one stream per (variant,tau,rho); same seed rule on both variants (registered)');
summary.n_sim = nSim;
summary.p_target = pTarget;
summary.estimands = {'p_surv (exact, rho-free)', 'p_est (exact + simulated)', ...
    'p_carmona (cited control)', 'p_*_static (quasi-static, obligation )'};
summary.sim_diagnostics = simDiag;
summary.max_violation_of_unit_interval = max(cellfun(@(v) A.(lower(char(v))).bw.max_violation_of_unit_interval, variants));
summary.profile_shape_note = 'see p_establishment_by_week.csv (peak week, number of interior local maxima, single-peaked flag)';
summary.peak_week_by_rho = struct();
for vi = 1:numel(variants)
    vv = lower(char(variants{vi}));
    bw = A.(vv).bw;
    pk = nan(numel(rhoGrid), 1);
    for ri = 1:numel(rhoGrid)
        pv = 1 - bw.u(:, 2, ri);
        [im, ~] = shapeOf(pv);
        pk(ri) = tauAll(im);
    end
    summary.peak_week_by_rho.(vv) = struct('rho', rhoGrid, 'peak_week', pk(:)');
end
summary.comparison = struct('tau_star_deterministic', tauStar, 'observation_start_week', obsStart);
summary.honesty = ['p_est peaks are MEASURED, not asserted; the profile is reported as ' ...
    'single-peaked only if the measured number of interior local maxima is <= 1. ' ...
    'Carmona & Gandon 2020 Eq.(2) is a cited control and is never presented as a ' ...
    'contribution of this paper (the statement prohibition).'];
mwspm_writejson(fullfile(ctx.dataDir, 'summary.json'), summary);

out = struct();
out.status = 'SUCCESS';
out.parent_run_id = fitRun.runId;
out.n_sim_points = simDiag.n_points;
out.max_abs_diff_in_mc_se = simDiag.max_abs_diff_in_mc_se;
out.deterministic = false;
out.rng = sprintf('mt19937ar, base_seed %d, one stream per (tau,rho)', baseSeed);
out.outputs = {'branching_analytic_grid.csv', 'branching_sim_subset.csv', ...
    'branching_dt_refinement.csv', 'p_establishment_by_week.csv', 'raw_choices.csv', 'summary.json'};
out.note = ['the analysis. The deterministic profile, the exact branching probabilities and the ' ...
    'replicate distribution are three different objects and are reported as such; ' ...
    'the 47 weekly samples of one trajectory are never used as replicates.'];
end

% ---------------------------------------------------------------------------
function writeAnalyticCsv(path, A, variants, tauAll, rhoGrid)
fid = fopen(path, 'w', 'n', 'UTF-8');
fprintf(fid, ['variant,tau_week,rho,R0_week_mean,r_week_mean,p_surv,p_est,', ...
    'p_surv_static,p_est_static,p_carmona\n']);
nW = numel(tauAll);
for vi = 1:numel(variants)
    vv = lower(char(variants{vi}));
    S = A.(vv);
    R = S.R; bw = S.bw; bwS = S.bwS; pc = S.pc;
    wk = floor(((1:R.n_days)' - 1) / 7) + 1;
    wk = min(max(wk, 1), nW);
    R0w = accumarray(wk, R.R0_daily, [nW 1], @mean, NaN);
    rw  = accumarray(wk, R.r_daily, [nW 1], @mean, NaN);
    for i = 1:nW
        for j = 1:numel(rhoGrid)
            fprintf(fid, '%s,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n', vv, ...
                tauAll(i), rhoGrid(j), R0w(i), rw(i), 1 - bw.q(i, 2), 1 - bw.u(i, 2, j), ...
                1 - bwS.q(i, 2), 1 - bwS.u(i, 2, j), pc(i));
        end
    end
end
fclose(fid);
end

% ---------------------------------------------------------------------------
function writeProfileCsv(path, profRows)
fid = fopen(path, 'w', 'n', 'UTF-8');
fprintf(fid, ['variant,rho,peak_week,peak_value,n_interior_local_maxima,single_peaked,', ...
    'local_max_weeks,tau_star_deterministic,p_at_tau_star,observation_start_week,', ...
    'p_at_observation_start,peak_minus_observation_start\n']);
for i = 1:numel(profRows)
    rw = profRows{i};
    fprintf(fid, '%s,%s,%d,%.17g,%d,%d,%s,%s,%s,%s,%s,%s\n', rw{1}, numOrNA(rw{2}), ...
        rw{3}, rw{4}, rw{5}, rw{6}, rw{7}, numOrNA(rw{8}), numOrNA(rw{9}), ...
        numOrNA(rw{10}), numOrNA(rw{11}), numOrNA(rw{12}));
end
fclose(fid);
end

% ---------------------------------------------------------------------------
function v = atWeek(p, w)
% Index a weekly profile at week w, returning NaN when w is undefined.  A
% guarded function is used instead of an inline ternary because MATLAB
% evaluates both branches of a ternary call, so p(NaN) would raise.
if isfinite(double(w)) && double(w) >= 1 && double(w) <= numel(p)
    v = p(round(double(w)));
else
    v = NaN;
end
end

% ---------------------------------------------------------------------------
function [iMax, lm] = shapeOf(p)
p = p(:);
[~, iMax] = max(p);
n = numel(p);
lm = [];
for i = 2:(n - 1)
    if p(i) >= p(i - 1) && p(i) > p(i + 1)
        lm(end + 1) = i; %#ok<AGROW>
    end
end
end

% ---------------------------------------------------------------------------
function [p, info] = carmona_control(R, tEval)
% Carmona & Gandon 2020 PLoS Comput Biol 16(7):e1007954 Eq.(2):
%     p_e(t0) = 1 / (1 + int_{t0}^{inf} mu(t) exp(-(phi(t) - phi(t0))) dt)
% with the horizon TRUNCATED at t_end (declared: with a finite upper limit the
% same expression is exactly the survival-to-t_end probability of the
% single-type birth-death process).
% The reduction matches the single-type birth-death on (R0, r):
%     r  = dominant eigenvalue of the frozen 4x4 mean matrix  (Malthusian)
%     R0 = next-generation ratio of the same frozen system
%     b - d = r,  b/d = R0  =>  d = r/(R0-1),  b = R0*d
% This is a DECLARED reduction of a two-host cycle to one type; it is a
% control, and its deviation from the exact 4-type answer is measured.
nd = R.n_days;
r = R.r_daily(:);
R0 = R.R0_daily(:);
d = nan(nd, 1);
ok = abs(R0 - 1) > 1e-9;
d(ok) = r(ok) ./ (R0(ok) - 1);
bad = ok & ~(d > 0);
d(bad) = NaN;
info = struct('n_guard', nnz(~ok) + nnz(bad), ...
    'min_d', tern(any(isfinite(d)), min(d(isfinite(d))), NaN), ...
    'max_d', tern(any(isfinite(d)), max(d(isfinite(d))), NaN));

Phi = [0; cumsum(r)];
p = nan(numel(tEval), 1);
for i = 1:numel(tEval)
    k0 = min(max(floor(tEval(i)) + 1, 1), nd);
    acc = 0;
    for k = k0:nd
        if ~isfinite(d(k))
            continue
        end
        base = Phi(k) - Phi(k0);
        if abs(r(k)) < 1e-12
            w = 1;
        else
            w = (1 - exp(-r(k))) / r(k);
        end
        acc = acc + d(k) * exp(-base) * w;
    end
    p(i) = 1 / (1 + acc);
end
end

% ---------------------------------------------------------------------------
function P = find_parent(root, wantId, expId, relPath)
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
    P = struct('runId', d(k).name, 'path', f);
end
end

% ---------------------------------------------------------------------------
function sv = findSolver(cfg, id)
% jsondecode returns a struct array only when all array elements have the same
% fields; the solver entries deliberately do not, so it returns a CELL array.
% Both shapes are accepted here.
sv = struct();
S = cfg.solvers;
for k = 1:numel(S)
    if iscell(S)
        item = S{k};
    else
        item = S(k);
    end
    if strcmp(char(string(getdef(item, 'id', ''))), id)
        sv = item;
        return
    end
end
error('branching:noSolver', 'no solver with id "%s" in the experiment config', id);
end

% ---------------------------------------------------------------------------
function s = strjoinLocal(v)
if isempty(v)
    s = 'none';
else
    s = strjoin(cellstr(string(v(:)')), ' ');
end
end
function s = bracket(v)
if isempty(v)
    s = '[]';
else
    s = sprintf('[%s]', strjoin(cellstr(string(v(:)')), ' '));
end
end
function s = numOrNA(v)
if isfinite(double(v)), s = sprintf('%.17g', double(v)); else s = 'NA'; end
end
function s = num2strOrNA(v)
if isfinite(double(v)), s = sprintf('%.6g', double(v)); else s = 'NA'; end
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
