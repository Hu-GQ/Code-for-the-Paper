function payload = exp_r2_pure_profile(ctx)
%EXP_R2_PURE_PROFILE  Pure-Poisson reporting-ratio profile under the stable solver.
%
%   Replaces the archived penalized Euler profile with a profile computed on the
%   same specification as exp_r1_stable_calibration: daily zero-order-held
%   coefficients, unclamped RK4 at 0.025 d, pure Poisson negative
%   log-likelihood and NO peak-week penalty.
%
%   For every rho on the grid the six remaining continuous parameters are
%   re-optimised from the previous grid point (continuation).  Because each fit
%   is conditional on a fixed rho, the result is a DESCRIPTIVE candidate curve.
%   It is not a confidence set, not a rho estimate and not a global optimum.
%   All exit flags and evaluation counts are recorded so nonconvergence stays
%   visible.
%
%   tau* is recomputed at every grid point with the same stable kernel, using
%   T0(w) = rho / log(R0_w) and weekly R0 from the continuous-rate traits.

cfg = ctx.cfg;
assert(~cfg.phase_penalty, 'R2:penalty', 'the revised profile must not use a peak penalty');

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

D   = mwspm_load_data(ctx.root, ctx.modelCard);
pre = mwspm_calib_prepare(D, ctx.modelCard, struct('mu_variant', cfg.mu_variant));
pRef = mwspm_params(ctx.modelCard);

warmRun = mwspm_readjson(fullfile(ctx.root, 'results', 'runs', cfg.warm_run, 'data', 'stable_fit.json'));

csvPath = fullfile(ctx.dataDir, 'rho_profile_pure.csv');
fid = fopen(csvPath, 'w', 'n', 'UTF-8');
fprintf(fid, ['lag,rho,nll,nll_full,exitflag,funcCount,m0,mu_s,onset_doy,offset_doy,' ...
    'gamma,beta,peak_week,reported_total,true_total,tau_star,T0_at_tau,cum_reports_at_tau,' ...
    'R0_gt1_start_week,dispersion\n']);
fclose(fid);

box = mwspm_calib_theta('box');
allRows = table();
summary = struct('interpretation', ...
    'descriptive conditional candidate profile; fixed rho; pure Poisson; no inference acceptance', ...
    'solver', cfg.solver, 'mu_variant', cfg.mu_variant, 'phase_penalty', false, ...
    'lags', cfg.lags, 'branches', struct(), 'rho_grid_by_lag', struct());

for d = double(cfg.lags(:))'
    key = sprintf('lag%d', d);
    assert(isfield(warmRun, key), 'R2:noWarm', 'warm start for %s missing', key);
    prev = warmRun.(key).best.theta;

    if isfield(cfg, 'rho_grid_override') && ~isempty(cfg.rho_grid_override)
        rhoGrid = double(cfg.rho_grid_override(:))';
    elseif isfield(cfg, 'rho_grid_by_lag') && isfield(cfg.rho_grid_by_lag, key)
        g = cfg.rho_grid_by_lag.(key);
        rhoGrid = double(g.from) : double(g.step) : double(g.to);
    else
        rhoGrid = double(cfg.rho_grid.from) : double(cfg.rho_grid.step) : double(cfg.rho_grid.to);
    end
    rhoGrid = round(rhoGrid, 6);
    summary.rho_grid_by_lag.(key) = rhoGrid;

    % Continuation order: walk outwards from the rho of the warm start so every
    % new grid point starts from an adjacent solved point.
    [~, ord] = sort(abs(rhoGrid - prev.rho));
    fprintf('R2 %s: %d grid points, warm rho = %.4f\n', key, numel(rhoGrid), prev.rho);

    rows = table();
    for idx = ord(:)'
        rho = rhoGrid(idx);
        C = struct('free', {cellstr(string(cfg.free))}, 'fixed', struct('rho', rho), ...
            'box', box, 'family', 'poisson', 'phase_penalty', false, ...
            'kernel_opts', cfg.solver, 'd_list', d);
        opt = cfg.optimiser;
        opt.warm = {prev};
        [best, info] = mwspm_calib_fit(pre, C, opt); %#ok<ASGLU>
        assert(abs(best.J - best.nll) < 1e-10, 'R2:penaltyLeak', 'objective is not the pure NLL');

        [tauW, r0Start, T0at, cumAt] = r2_tau_star(pre, D, pRef, best.theta, rho, ...
            cfg.solver, struct());

        row = table(d, rho, best.nll, best.nll_full, best.exitflag, info.n_eval, ...
            best.theta.m0, best.theta.mu_s, best.theta.onset_doy, best.theta.offset_doy, ...
            best.theta.gamma, best.theta.beta, best.peak_week, best.reported_total, ...
            best.true_total, tauW, T0at, cumAt, r0Start, best.dispersion, ...
            'VariableNames', {'lag','rho','nll','nll_full','exitflag','funcCount','m0','mu_s', ...
            'onset_doy','offset_doy','gamma','beta','peak_week','reported_total', ...
            'true_total','tau_star','T0_at_tau','cum_reports_at_tau','R0_gt1_start_week','dispersion'});
        rows = [rows; row]; %#ok<AGROW>
        allRows = [allRows; row]; %#ok<AGROW>

        fid = fopen(csvPath, 'a', 'n', 'UTF-8');
        fprintf(fid, '%d,%.10g,%.10g,%.10g,%d,%d,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%d,%.10g,%.10g,%s,%s,%s,%s,%.10g\n', ...
            d, rho, best.nll, best.nll_full, best.exitflag, info.n_eval, ...
            best.theta.m0, best.theta.mu_s, best.theta.onset_doy, best.theta.offset_doy, ...
            best.theta.gamma, best.theta.beta, best.peak_week, best.reported_total, ...
            best.true_total, r2_num(tauW), r2_num(T0at), r2_num(cumAt), r2_num(r0Start), ...
            best.dispersion);
        fclose(fid);

        fprintf(['R2 %s rho=%.2f nll=%.6f flag=%d evals=%d peak=%d ' ...
            'tau*=%s T0=%s cum=%s\n'], ...
            key, rho, best.nll, best.exitflag, info.n_eval, best.peak_week, ...
            r2_num(tauW), r2_num(T0at), r2_num(cumAt));
        prev = best.theta;
    end

    fin = isfinite(rows.tau_star);
    b = struct();
    b.d_weeks = d;
    b.n_points = height(rows);
    b.n_converged = sum(rows.exitflag == 1);
    b.n_at_limit = sum(rows.exitflag == 0);
    b.tau_star_min = min(rows.tau_star(fin));
    b.tau_star_max = max(rows.tau_star(fin));
    b.tau_star_unique = unique(rows.tau_star(fin))';
    b.nll_min = min(rows.nll);
    b.rho_at_nll_min = rows.rho(rows.nll == min(rows.nll));
    b.nll_at_rho_080 = rows.nll(abs(rows.rho - 0.8) < 1e-9);
    b.peak_week_at_nll_min = rows.peak_week(rows.nll == min(rows.nll));
    b.tau_star_at_nll_min = rows.tau_star(rows.nll == min(rows.nll));
    b.T0_min = min(rows.T0_at_tau(fin));
    b.T0_max = max(rows.T0_at_tau(fin));
    b.cum_min = min(rows.cum_reports_at_tau(fin));
    b.cum_max = max(rows.cum_reports_at_tau(fin));
    b.reported_total_range = [min(rows.reported_total), max(rows.reported_total)];
    b.true_total_range = [min(rows.true_total), max(rows.true_total)];
    summary.branches.(key) = b;

    fprintf('R2 %s summary: tau* in [%s, %s]; nll min %.6f at rho=%s; converged %d/%d\n', ...
        key, r2_num(b.tau_star_min), r2_num(b.tau_star_max), b.nll_min, ...
        mat2str(b.rho_at_nll_min(:)'), b.n_converged, b.n_points);
end

% Solver agreement at the lowest-nll point of the primary branch.
[~, iBest] = min(allRows.nll(allRows.lag == cfg.lags(1)));
prim = allRows(allRows.lag == cfg.lags(1), :);
th = struct('m0', prim.m0(iBest), 'rho', prim.rho(iBest), 'd_weeks', prim.lag(iBest), ...
    'onset_doy', prim.onset_doy(iBest), 'offset_doy', prim.offset_doy(iBest), ...
    'mu_s', prim.mu_s(iBest), 'gamma', prim.gamma(iBest), 'beta', prim.beta(iBest));
kA  = mwspm_calib_kernel(pre, th, cfg.solver);
kH  = mwspm_calib_kernel(pre, th, struct('solver', 'rk4', 'dt', cfg.solver.dt / 2));
kO  = mwspm_calib_kernel(pre, th, struct('solver', 'ode15s'));
solverCheck = struct('at_lag', th.d_weeks, 'at_rho', th.rho, ...
    'step_halving_max_weekly', max(abs(kA.Yhat - kH.Yhat)), ...
    'ode15s_max_weekly', max(abs(kA.Yhat - kO.Yhat)), ...
    'ode15s_nll_difference', abs(sum(kA.Yhat - pre.Y .* log(max(kA.Yhat, 1e-12))) - ...
                                 sum(kO.Yhat - pre.Y .* log(max(kO.Yhat, 1e-12)))));
mwspm_writejson(fullfile(ctx.dataDir, 'solver_check.json'), solverCheck);
mwspm_writejson(fullfile(ctx.dataDir, 'profile_summary.json'), summary);

payload = summary;
payload.outputs = cfg.outputs;
end

% ---------------------------------------------------------------------------
function [tauStar, r0StartWeek, T0atTau, cumAtTau] = r2_tau_star(pre, D, pRef, th, rhoUsed, solverOpts, trOpts)
% Weekly R0 (continuous-rate traits) and cumulative expected reports, then tau*.
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
function s = r2_num(v)
if isfinite(v)
    s = sprintf('%.6g', v);
else
    s = 'NA';
end
end
