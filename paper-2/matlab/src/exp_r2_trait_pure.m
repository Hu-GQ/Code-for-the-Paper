function payload = exp_r2_trait_pure(ctx)
%EXP_R2_TRAIT_PURE  Trait-substitution refits under the pure-Poisson stable spec.
%
%   Three trait sets are refitted on the specification of
%   exp_r1_stable_calibration: daily zero-order-held coefficients, unclamped RK4
%   at 0.025 d, pure Poisson negative log-likelihood, no peak-week penalty.
%
%   Every set keeps the source reporting-ratio bound rho in [0.02, 0.80] and the
%   same seven free parameters, so the Delta-BIC between sets is computed on one
%   objective and one parameter count.  tau* is recomputed with the same stable
%   kernel and the set-specific continuous-rate traits.
%
%   Sets (definitions copied from exp_t14_trait_source_swap[_c1].json):
%     albopictus_tableA_baseline  Mordecai Table A, Ae. albopictus (project default)
%     aegypti_tableB_swap         Mordecai Table B, Ae. aegypti   (same paper)
%     hybrid_LH2014_a_sigma_mu    a/sigma_v/mu from Liu-Helmersson 2014, b/c kept
%
%   The third set is a partial cross-paper mixture, not a complete alternative
%   source.  None of the three is a validation cohort.

cfg = ctx.cfg;
assert(~cfg.phase_penalty, 'R2b:penalty', 'the revised refit must not use a peak penalty');
assert(isfield(cfg, 'trait_set_order') && ~isempty(cfg.trait_set_order), 'R2b:order', ...
    'trait_set_order must be explicit: the first set is the baseline');

if isfield(ctx, 'dryRun') && ctx.dryRun
    payload = struct('outputs', {{}});
    payload.dry_run = true;
    return
end

diary(fullfile(ctx.runDir, 'logs', 'stdout.txt'));
diaryGuard = onCleanup(@() diary('off')); %#ok<NASGU>
snap = fullfile(ctx.dataDir, 'code_snapshot');
mkdir(snap);
copyfile(fullfile(ctx.root, 'src', 'matlab', '*.m'), snap);
copyfile(fullfile(ctx.root, 'scripts', 'reproduce.m'), snap);

D    = mwspm_load_data(ctx.root, ctx.modelCard);
pRef = mwspm_params(ctx.modelCard);
box  = mwspm_calib_theta('box');
warmRun = mwspm_readjson(fullfile(ctx.root, 'results', 'runs', cfg.warm_run, 'data', 'stable_fit.json'));
warm0 = warmRun.(sprintf('lag%d', cfg.lag)).best.theta;

names = cellstr(string(cfg.trait_set_order));
rows = table();

for i = 1:numel(names)
    nm  = names{i};
    assert(isfield(cfg.trait_sets, nm), 'R2b:badSet', 'trait set "%s" is not defined', nm);
    S = cfg.trait_sets.(nm);

    prepOpts = struct('mu_variant', cfg.mu_variant);
    trOpts   = struct();
    if isfield(S, 'coef') && ~isempty(S.coef), prepOpts.coef = S.coef; trOpts.coef = S.coef; end
    if isfield(S, 'alt')  && ~isempty(S.alt)
        prepOpts.alt = S.alt;
        trOpts.alt   = S.alt;
    end
    pre = mwspm_calib_prepare(D, ctx.modelCard, prepOpts);

    freeNow = cellstr(string(cfg.free));
    fixedNow = struct();
    if isfield(cfg, 'fixed_rho') && ~isempty(cfg.fixed_rho)
        % Trait comparison at a COMMON reporting ratio.  With a nearly flat
        % rho profile, letting rho float makes each trait set settle at a
        % different point of the ridge, so the onset-week difference would
        % measure the rho drift rather than the trait curves.
        freeNow = setdiff(freeNow, {'rho'}, 'stable');
        fixedNow.rho = double(cfg.fixed_rho);
    end
    C = struct('free', {freeNow}, 'fixed', fixedNow, 'box', box, ...
        'family', 'poisson', 'phase_penalty', false, 'kernel_opts', cfg.solver, ...
        'd_list', cfg.lag);
    opt = cfg.optimiser;
    opt.warm = {warm0};
    % Optional extra start from an earlier run of this experiment (the same
    % trait set), so a long refinement is not re-explored from scratch.
    if isfield(cfg, 'warm_from_run') && ~isempty(cfg.warm_from_run)
        fPrev = fullfile(ctx.root, 'results', 'runs', cfg.warm_from_run, ...
            'data', sprintf('refit_%s.json', nm));
        if exist(fPrev, 'file') == 2
            opt.warm{end + 1} = mwspm_readjson(fPrev).theta;
        end
    end
    [best, info] = mwspm_calib_fit(pre, C, opt); %#ok<ASGLU>
    assert(abs(best.J - best.nll) < 1e-10, 'R2b:penaltyLeak', 'objective is not the pure NLL');

    [tauW, r0Start, T0at, cumAt] = r2b_tau_star(pre, D, pRef, best.theta, best.theta.rho, ...
        cfg.solver, trOpts);
    kB = mwspm_calib_kernel(pre, best.theta, struct('solver', 'rk4', 'dt', cfg.solver.dt / 2));
    kO = mwspm_calib_kernel(pre, best.theta, struct('solver', 'ode15s'));
    solverErr = max([max(abs(best.Yhat - kB.Yhat)), max(abs(best.Yhat - kO.Yhat))]);

    row = table({nm}, numel(freeNow), best.nll, best.nll_full, best.exitflag, info.n_eval, ...
        best.theta.m0, best.theta.rho, best.theta.mu_s, best.theta.onset_doy, ...
        best.theta.offset_doy, best.theta.gamma, best.theta.beta, best.peak_week, ...
        best.reported_total, best.true_total, tauW, T0at, cumAt, r0Start, solverErr, ...
        'VariableNames', {'trait_set','free_count','nll','nll_full','exitflag','funcCount','m0','rho', ...
        'mu_s','onset_doy','offset_doy','gamma','beta','peak_week','reported_total', ...
        'true_total','tau_star','T0_at_tau','cum_reports_at_tau','R0_gt1_start_week','solver_error'});
    rows = [rows; row]; %#ok<AGROW>

    fprintf(['R2b %-28s nll=%.6f flag=%d evals=%d rho=%.4f peak=%d ' ...
        'tau*=%s T0=%s cum=%s solverErr=%.2e\n'], ...
        nm, best.nll, best.exitflag, info.n_eval, best.theta.rho, best.peak_week, ...
        r2b_num(tauW), r2b_num(T0at), r2b_num(cumAt), solverErr);
    mwspm_writejson(fullfile(ctx.dataDir, sprintf('refit_%s.json', nm)), best);
end

writetable(rows, fullfile(ctx.dataDir, 'trait_refit_pure.csv'));

% Delta-BIC against the baseline set, same objective and same parameter count.
nObs = numel(pre.Y);
kPar = rows.free_count(1);
bic = 2 * rows.nll_full + kPar * log(nObs);
d = struct();
d.interpretation = ['Descriptive conditional refits under one specification; ' ...
    'Delta-BIC is not evidence that the fit discriminates trait sets.'];
d.solver = cfg.solver;
d.mu_variant = cfg.mu_variant;
d.n_obs = nObs;
d.n_parameters = kPar;
d.bic = bic(:)';
d.delta_bic_vs_baseline = (bic - bic(1))';
d.tau_star = rows.tau_star(:)';
d.tau_star_shift_vs_baseline = (rows.tau_star - rows.tau_star(1))';
d.sets = rows.trait_set(:)';
d.all_converged = all(rows.exitflag == 1);
d.reported_total = rows.reported_total(:)';
d.true_total = rows.true_total(:)';
d.rho_hat = rows.rho(:)';
mwspm_writejson(fullfile(ctx.dataDir, 'trait_summary.json'), d);

payload = d;
payload.outputs = cfg.outputs;
end

% ---------------------------------------------------------------------------
function [tauStar, r0StartWeek, T0atTau, cumAtTau] = r2b_tau_star(pre, D, pRef, th, rhoUsed, solverOpts, trOpts)
nW = pre.n_weeks;

sp = struct('onset_doy', th.onset_doy, 'offset_doy', th.offset_doy, ...
    'kappa_days', pre.kappa_days, 'gamma', th.gamma, 'beta', th.beta, ...
    's_floor', pre.s_floor, 'min_window_days', pre.min_window_days, ...
    'normalize', true, 'norm_factor', NaN);
s = mwspm_season(pre.doy, sp);

R0daily = mwspm_r0(D.T_daily, th.m0 * s, th.mu_s * pre.mu_raw_daily, pRef.gamma_h, trOpts);

wk = floor(((1:numel(R0daily))' - 1) / 7) + 1;
wk = min(max(wk, 1), nW);

gt1 = find(R0daily > 1, 1);
if isempty(gt1)
    r0StartWeek = NaN;
else
    r0StartWeek = wk(gt1);
end

R0w = accumarray(wk, R0daily, [nW 1], @mean, NaN);
kern = mwspm_calib_kernel(pre, th, solverOpts);
cumRep = cumsum(kern.Yhat(:));

T0w = nan(nW, 1);
ok = R0w > 1;
T0w(ok) = rhoUsed ./ log(R0w(ok));

cross = find(isfinite(T0w) & cumRep >= T0w, 1);
if isempty(cross)
    tauStar = NaN; T0atTau = NaN; cumAtTau = NaN;
else
    tauStar = cross; T0atTau = T0w(cross); cumAtTau = cumRep(cross);
end
end

% ---------------------------------------------------------------------------
function s = r2b_num(v)
if isfinite(v)
    s = sprintf('%.6g', v);
else
    s = 'NA';
end
end
