function out = exp_t14_identifiability(ctx)
%EXP_T14_IDENTIFIABILITY  T14.03 -- are onset_doy and gamma SEPARATELY
%   identifiable? (EXPERIMENT_PLAN.md Sec.5-5, the C3 legacy question.)
%
%   Four independent pieces of evidence, all at the fitted integer lag d = 1:
%     (A) 7-D local sensitivity / collinearity (mwspm_calib_sensitivity):
%         SVD, Belsley condition indices, collinear sets, and the correlation
%         matrix of the least-squares ESTIMATES;
%     (B) multi-start scatter: 64 explicit starts on an 8x8 (onset, gamma)
%         grid, each followed by one Nelder-Mead run with NO random screening
%         (T12 showed that screen spread measures sampling distance, not
%         optimiser reliability), giving one local optimum per start;
%     (C) 1-D profiles of onset and of gamma: fix one, re-optimise the other
%         six, and take the Wilks interval at chi2_1(0.95)/2 = 1.9207;
%     (D) 2-D joint profile: fix BOTH, re-optimise the remaining five, on a
%         9x9 grid, and map the joint 95% region at chi2_2(0.95)/2 = 2.9957.
%
%   The verdict rule is PRE-REGISTERED in the config (acceptance_criterion):
%   an interval that touches the grid boundary is NOT closed inside the grid,
%   and that alone makes the parameter not identifiable on this grid - the
%   grid is not quietly widened afterwards.
%
%   Failure modes this experiment is built to avoid:
%     * random-screening failure masquerading as non-identifiability (all fits
%       use the reference point as a warm start with n_screen = 0);
%     * reporting a 1-D interval on a ridged likelihood (the 2-D region is
%       always reported next to it);
%     * treating the multi-start SPREAD as a confidence region (only the
%       near-optimal set with dJ <= Wilks is used).

logFid = ctx.logFid;
root   = ctx.root;
cfg    = ctx.cfg;
modelCard = ctx.modelCard;

fprintf(logFid, '=== exp_t14_identifiability (T14.03, Sec.5-5) ===\n');
fprintf(logFid, 'claims: %s\n', strjoin(cellstr(string(getdef(cfg, 'claims', {}))), ', '));

mc  = mwspm_readjson(fullfile(root, 'config', 'experiments', 'model_mwspm_gen3_rhohat.json'));
D   = mwspm_load_data(root, modelCard);
pre = mwspm_calib_prepare(D, modelCard, struct('mu_variant', 'day'));
box = mwspm_calib_theta('box');

fitRun = t1403_find_parent(root, getdef(cfg, 'parent_run_id', []), 'exp_c4_rho_fit', 'data/fit_best.json');
assert(~isempty(fitRun), 't14:noParent', ...
    'no SUCCESS run of exp_c4_rho_fit with data/fit_best.json found under results/runs/');
fit = mwspm_readjson(fitRun.path);
th  = t1403_theta_from_fit(fit.main_report_point.theta);
Jref = double(fit.main_report_point.nll);
fprintf(logFid, 'parent fit run : %s\n', fitRun.runId);
fprintf(logFid, 'reference point: m0=%.6f rho=%.6f onset=%.4f offset=%.4f mu_s=%.5f gamma=%.5f beta=%.5f d=%d (nll=%.4f)\n', ...
    th.m0, th.rho, th.onset_doy, th.offset_doy, th.mu_s, th.gamma, th.beta, th.d_weeks, Jref);

PD = cfg.parameter_domain;

% PARTS selector (P4, 2026-10-07).  A refinement run may skip the expensive
% multi-start scatter (B, ~4 min) when only the joint grid needs a finer
% resolution.  A always runs (0.1 s); C and D are MANDATORY because both the
% pre-registered verdict and the projection control are built on them.
% Skipping B never widens or coarsens a grid: the grids live in the config.
parts = cellstr(string(getdef(PD, 'parts', {'B', 'C', 'D'})));
assert(all(ismember({'C', 'D'}, parts)), 't14:parts', ...
    'parameter_domain.parts must contain C and D (the verdict and the pre-registered projection control depend on them)');
doB = ismember('B', parts);
fprintf(logFid, 'parts          : %s (multi-start scatter %s)\n', ...
    strjoin(parts, '+'), yesNo(doB));

if isfield(ctx, 'dryRun') && ctx.dryRun
    out = struct();
    out.dry_run = true;
    out.parent_run_id = fitRun.runId;
    out.validated = struct('n_multistart', PD.multistart_grid.onset_doy.n * PD.multistart_grid.gamma.n, ...
        'n_onset_profile', PD.profile_onset.n, 'n_gamma_profile', PD.profile_gamma.n, ...
        'n_joint', PD.joint_onset.n * PD.joint_gamma.n, 'fixed_lag', PD.fixed_lag);
    fprintf(logFid, '[dry-run] parent and grid definitions readable; no optimisation performed\n');
    return
end

freeAll  = cellstr(string(PD.free_full));
free6on  = cellstr(string(PD.refit_6_when_1d));
free6ga  = cellstr(string(PD.refit_6_when_1d_gamma));
free5    = cellstr(string(PD.refit_5_when_2d));
dFix     = double(PD.fixed_lag);

% Guard: the profiled parameter must be FIXED and must NOT be free.  The first
% implementation of this experiment put gamma in both lists, so the kernel
% received a theta without onset_doy and failed; the failed run is retained.
assert(isequal(sort(free6on), sort(setdiff(freeAll, 'onset_doy'))), 't14:badConfig', ...
    'parameter_domain.refit_6_when_1d must equal free_full minus onset_doy');
assert(isequal(sort(free6ga), sort(setdiff(freeAll, 'gamma'))), 't14:badConfig', ...
    'parameter_domain.refit_6_when_1d_gamma must equal free_full minus gamma');
assert(isequal(sort(free5), sort(setdiff(freeAll, {'onset_doy', 'gamma'}))), 't14:badConfig', ...
    'parameter_domain.refit_5_when_2d must equal free_full minus {onset_doy, gamma}');

oCfg = cfg.solvers;
sv = oCfg(1);
baseOpt = struct('n_screen', 0, 'n_seed', 0, 'restarts', 0, ...
    'maxFunEvals', double(sv.maxFunEvals), 'maxIter', double(sv.maxIter), ...
    'TolX', double(sv.TolX), 'TolFun', double(sv.TolFun), 'verbose', false);

mkC = @(free, fixed) struct('free', {free}, 'fixed', fixed, 'box', box, ...
    'family', 'poisson', 'family_opts', struct(), 'sub_weeks', [], 'y_obs', [], ...
    'phase_penalty', true, 'target_peak', 37.0, 'pen', 300.0, 'd_list', dFix);

% ---------------------------------------------------------------- control ----
% the reference point must be reproducible as a local optimum of this objective
C0 = mkC(freeAll, struct());
t0 = tic;
[bRef, iRef] = mwspm_calib_fit(pre, C0, mergeOpt(baseOpt, struct('warm', {{th}})));
tRef = toc(t0);
dRef = bRef.J - Jref;
fprintf(logFid, '\n--- control: NM from the reference point ---\n');
fprintf(logFid, '  J=%.6f (parent nll=%.6f, dJ=%+.3e) in %.2f s | onset=%.4f gamma=%.5f | n_eval=%d\n', ...
    bRef.J, Jref, dRef, tRef, bRef.theta.onset_doy, bRef.theta.gamma, iRef.n_eval);
assert(dRef <= 1e-2, 't14:controlFailed', ...
    ['NM from the reference point returned J=%.6f, worse than the parent nll %.6f by more than ' ...
     '1e-2 - the objective or the bounded map must have changed'], bRef.J, Jref);
if dRef > 1e-6
    fprintf(logFid, ['  NOTE: NM stopped %.3e nll units ABOVE the parent optimum (tolerance 1e-2; ' ...
        'a failure of this control would be an O(100) discrepancy, not O(1e-4)).\n'], dRef);
end

% ------------------------------------------------------- (A) sensitivity -----
fprintf(logFid, '\n--- (A) 7-D sensitivity / collinearity at the reference point ---\n');
S = mwspm_calib_sensitivity(pre, th, freeAll, struct('rel_step', 1e-4));
iOn = find(strcmp(freeAll, 'onset_doy'));
iGa = find(strcmp(freeAll, 'gamma'));
corrOG = S.corr_estimates(iOn, iGa);
fprintf(logFid, '  svd            : %s\n', mat2str(S.svd(:)', 5));
fprintf(logFid, '  cond number    : %.4f (practical rank %d/%d)\n', S.cond_number, S.practical_rank, numel(freeAll));
fprintf(logFid, '  colnorm        : %s\n', mat2str(S.colnorm(:)', 5));
fprintf(logFid, '  corr(onset,gam): %.6f  (from (Sh''Sh)^-1)\n', corrOG);
fprintf(logFid, '  collinear sets : %s\n', t1403_cellstr(S.collinear_sets));

% --------------------------------------------------- (B) multi-start scatter -
if doB
gOn = linspace(th.onset_doy - PD.multistart_grid.onset_doy.halfwidth, ...
               th.onset_doy + PD.multistart_grid.onset_doy.halfwidth, PD.multistart_grid.onset_doy.n);
gGa = linspace(th.gamma - PD.multistart_grid.gamma.halfwidth, ...
               th.gamma + PD.multistart_grid.gamma.halfwidth, PD.multistart_grid.gamma.n);
msRows = {};
fprintf(logFid, '\n--- (B) multi-start scatter: %d x %d = %d explicit starts ---\n', ...
    numel(gOn), numel(gGa), numel(gOn) * numel(gGa));
tB = tic;
for a = 1:numel(gOn)
    for g = 1:numel(gGa)
        st = th; st.onset_doy = gOn(a); st.gamma = gGa(g);
        [b, ~] = mwspm_calib_fit(pre, C0, mergeOpt(baseOpt, struct('warm', {{st}})));
        msRows{end+1, 1} = {gOn(a), gGa(g), b.J, b.J - bRef.J, b.theta.onset_doy, ...
            b.theta.gamma, b.theta.m0, b.theta.rho, b.theta.mu_s, b.theta.offset_doy, b.theta.beta}; %#ok<AGROW>
    end
end
tB = toc(tB);
msM = t1403_num(msRows);
dJ  = msM(:, 4);
near = dJ <= double(cfg.acceptance_criterion.wilks_1d_dnll);
nNear = sum(near);
onNear = msM(near, 5); gaNear = msM(near, 6);
if nNear >= 3
    % Pearson correlation computed by hand: `corr` lives in the Statistics and
    % Machine Learning Toolbox, which is NOT installed on this installation.
    rNear = t1403_pearson(onNear, gaNear);
    Cv = cov(onNear, gaNear);
    [eV, ~] = eig(Cv);
    dirVec = eV(:, end);                       % principal axis of the near-optimal set
    dirAngle = atan2(dirVec(2), dirVec(1)) * 180 / pi;
else
    rNear = NaN; dirVec = [NaN; NaN]; dirAngle = NaN;
end
fprintf(logFid, '  %d/%d starts reached dJ <= %.4f (%.1f s total, %.2f s per start)\n', ...
    nNear, size(msM, 1), double(cfg.acceptance_criterion.wilks_1d_dnll), tB, tB / size(msM, 1));
fprintf(logFid, '  near-optimal onset range = [%.3f, %.3f] (width %.3f d)\n', min(onNear), max(onNear), max(onNear) - min(onNear));
fprintf(logFid, '  near-optimal gamma range = [%.5f, %.5f] (width %.5f)\n', min(gaNear), max(gaNear), max(gaNear) - min(gaNear));
fprintf(logFid, '  corr(onset,gamma) over near-optimal set = %.6f | principal axis angle = %.2f deg\n', rNear, dirAngle);
fprintf(logFid, '  distinct optima (onset rounded to 1 d, gamma to 0.01) = %d\n', ...
    numel(unique(round(onNear) * 1000 + round(gaNear * 100))));
else
    msM = nan(0, 11); nNear = 0; onNear = []; gaNear = [];
    rNear = NaN; dirVec = [NaN; NaN]; dirAngle = NaN; tB = NaN;
    fprintf(logFid, ['\n--- (B) multi-start scatter SKIPPED (parts=%s) ---\n' ...
        '  The scatter is not computed in this run; the near-optimal-set correlation must be\n' ...
        '  taken from the parent T14.03 run (see summary.parent_runs) and is reported as NaN here.\n'], ...
        strjoin(parts, '+'));
end

% ---------------------------------------------------------- (C) 1-D profiles -
pOnGrid = th.onset_doy + (-PD.profile_onset.halfwidth : PD.profile_onset.step : PD.profile_onset.halfwidth);
pGaGrid = th.gamma + (-PD.profile_gamma.halfwidth : PD.profile_gamma.step : PD.profile_gamma.halfwidth);
[pOn, onRows] = t1403_profile1d(pre, mkC, baseOpt, free6on, th, 'onset_doy', pOnGrid, dFix);
[pGa, gaRows] = t1403_profile1d(pre, mkC, baseOpt, free6ga, th, 'gamma', pGaGrid, dFix);
fprintf(logFid, '\n--- (C) 1-D profiles (Wilks dJ <= 1.9207) ---\n');
fprintf(logFid, '  onset_doy : min at %.4f (dJ=%.4f) | 95%% interval = [%s] | bounded inside grid = %s\n', ...
    pOn.minAt, pOn.minDJ, t1403_fmtInterval(pOn.lo, pOn.hi), yesNo(pOn.bounded));
fprintf(logFid, '  gamma     : min at %.5f (dJ=%.4f) | 95%% interval = [%s] | bounded inside grid = %s\n', ...
    pGa.minAt, pGa.minDJ, t1403_fmtInterval(pGa.lo, pGa.hi), yesNo(pGa.bounded));

% ------------------------------------------------------- (D) 2-D joint profile
jOn = th.onset_doy + (-PD.joint_onset.halfwidth : PD.joint_onset.step : PD.joint_onset.halfwidth);
jGa = th.gamma + (-PD.joint_gamma.halfwidth : PD.joint_gamma.step : PD.joint_gamma.halfwidth);
dJ2 = nan(numel(jOn), numel(jGa));
jRows = {};
fprintf(logFid, '\n--- (D) 2-D joint profile: %d x %d = %d fits (5 free) ---\n', ...
    numel(jOn), numel(jGa), numel(jOn) * numel(jGa));
tD = tic;
for a = 1:numel(jOn)
    for g = 1:numel(jGa)
        fixed = struct('onset_doy', jOn(a), 'gamma', jGa(g));
        st = th; st.onset_doy = jOn(a); st.gamma = jGa(g);
        [b, ~] = mwspm_calib_fit(pre, mkC(free5, fixed), mergeOpt(baseOpt, struct('warm', {{st}})));
        dJ2(a, g) = b.J - bRef.J;
        jRows{end+1, 1} = {jOn(a), jGa(g), b.J, b.J - bRef.J, b.theta.m0, b.theta.rho, ...
            b.theta.offset_doy, b.theta.mu_s, b.theta.beta}; %#ok<AGROW>
    end
end
tD = toc(tD);
cr2 = double(cfg.acceptance_criterion.wilks_2d_dnll);
inside = dJ2 <= cr2;
nIn = sum(inside(:));
touch = false;
if nIn > 0
    [ia, ig] = find(inside);
    touch = any(ia == 1) || any(ia == numel(jOn)) || any(ig == 1) || any(ig == numel(jGa));
    onExt = [min(jOn(ia)), max(jOn(ia))];
    gaExt = [min(jGa(ig)), max(jGa(ig))];
else
    onExt = [NaN, NaN]; gaExt = [NaN, NaN];
end
fprintf(logFid, '  %d/%d grid points inside the 95%% joint region (dJ <= %.4f), %.1f s\n', ...
    nIn, numel(dJ2), cr2, tD);
fprintf(logFid, '  joint region extent: onset [%s] | gamma [%s] | touches grid boundary = %s\n', ...
    t1403_fmtInterval(onExt(1), onExt(2)), t1403_fmtInterval(gaExt(1), gaExt(2)), yesNo(touch));

% PRE-REGISTERED CONTROL (config acceptance_criterion.must_report / controls):
% the 2-D region's projection onto an axis must NOT be narrower than that
% parameter's 1-D Wilks interval.  If it is, the joint grid is too COARSE to
% resolve that direction and the 2-D "does not touch the boundary" statement is
% not informative there.  The grid is NOT changed after seeing this - the
% violation is reported.
projNarrower = false;
projDetail = 'not evaluated (empty joint region)';
if nIn > 0 && isfinite(pOn.lo) && isfinite(pOn.hi)
    projNarrower = (onExt(1) > pOn.lo + 1e-9) || (onExt(2) < pOn.hi - 1e-9);
    projDetail = sprintf('2-D onset projection [%s] vs 1-D interval [%s]', ...
        t1403_fmtInterval(onExt(1), onExt(2)), t1403_fmtInterval(pOn.lo, pOn.hi));
    fprintf(logFid, '  CONTROL projection vs 1-D : %s -> %s\n', projDetail, ...
        t1403_yn(projNarrower, 'VIOLATED (joint grid too coarse in the onset direction; the 2-D boundedness claim is not informative for onset)', 'consistent'));
end

% ---------------------------------------------------------------- verdict ----
if ~pOn.bounded || ~pGa.bounded || touch
    verdict = 'NOT_IDENTIFIABLE';
elseif pOn.bounded && pGa.bounded && ~touch
    verdict = 'IDENTIFIABLE';
else
    verdict = 'PARTIALLY_IDENTIFIABLE';
end
fprintf(logFid, '\nVERDICT = %s\n', verdict);
if ~doB
    jointTxt = ['近优集合的 onset-gamma 相关性在本运行中未计算（parts 不含 B）⇒ ' ...
        '联合约束的相关系数须取自父 T14.03 运行，不得用本运行的 NaN 作任何推论'];
elseif abs(rNear) >= double(cfg.acceptance_criterion.corr_threshold)
    jointTxt = sprintf(['onset 与 gamma 在近优集合上 |r|=%.4f >= %.2f ⇒ 不可分别识别；' ...
        '只能给脊方向上的联合约束（主轴 %.2f 度）'], abs(rNear), ...
        double(cfg.acceptance_criterion.corr_threshold), dirAngle);
else
    jointTxt = sprintf(['onset 与 gamma 在近优集合上 |r|=%.4f < %.2f ⇒ 无单一二参数退化；' ...
        '两参数须分别报告区间，并同时给出联合域'], abs(rNear), ...
        double(cfg.acceptance_criterion.corr_threshold));
end
fprintf(logFid, 'JOINT CONSTRAINT: %s\n', jointTxt);

% ---------------------------------------------------------------- outputs ----
mwspm_writejson(fullfile(ctx.dataDir, 'sensitivity.json'), S);

hdr = {'start_onset_doy', 'start_gamma', 'J', 'dJ_vs_reference', 'opt_onset_doy', ...
    'opt_gamma', 'opt_m0', 'opt_rho', 'opt_mu_s', 'opt_offset_doy', 'opt_beta'};
if doB
    t1403_writecsv(fullfile(ctx.dataDir, 'multistart_scatter.csv'), hdr, msM);
end

hdr1 = {'onset_doy', 'J', 'dJ', 'opt_gamma', 'opt_m0', 'opt_rho', 'opt_offset_doy', 'opt_mu_s', 'opt_beta'};
t1403_writecsv(fullfile(ctx.dataDir, 'profile_1d_onset.csv'), hdr1, t1403_num(onRows));
hdr2 = {'gamma', 'J', 'dJ', 'opt_onset_doy', 'opt_m0', 'opt_rho', 'opt_offset_doy', 'opt_mu_s', 'opt_beta'};
t1403_writecsv(fullfile(ctx.dataDir, 'profile_1d_gamma.csv'), hdr2, t1403_num(gaRows));
hdr3 = {'onset_doy', 'gamma', 'J', 'dJ', 'opt_m0', 'opt_rho', 'opt_offset_doy', 'opt_mu_s', 'opt_beta'};
t1403_writecsv(fullfile(ctx.dataDir, 'joint_profile_2d.csv'), hdr3, t1403_num(jRows));

verdictOut = struct();
verdictOut.verdict = verdict;
verdictOut.rule = cfg.acceptance_criterion.verdict_rule;
verdictOut.onset_profile = struct('min_at', pOn.minAt, 'min_dJ', pOn.minDJ, ...
    'ci95_low', pOn.lo, 'ci95_high', pOn.hi, 'bounded_inside_grid', pOn.bounded, ...
    'grid', [min(pOnGrid), max(pOnGrid)]);
verdictOut.gamma_profile = struct('min_at', pGa.minAt, 'min_dJ', pGa.minDJ, ...
    'ci95_low', pGa.lo, 'ci95_high', pGa.hi, 'bounded_inside_grid', pGa.bounded, ...
    'grid', [min(pGaGrid), max(pGaGrid)]);
verdictOut.joint_2d = struct('n_inside', nIn, 'n_grid', numel(dJ2), 'dJ_threshold', cr2, ...
    'onset_extent', onExt, 'gamma_extent', gaExt, 'touches_boundary', touch, ...
    'grid_onset', [min(jOn), max(jOn)], 'grid_gamma', [min(jGa), max(jGa)], ...
    'projection_control_violated', projNarrower, 'projection_control_detail', projDetail);
verdictOut.multistart = struct('computed', doB, 'n_starts', size(msM, 1), 'n_near_optimal', nNear, ...
    'dJ_threshold', double(cfg.acceptance_criterion.wilks_1d_dnll), ...
    'onset_range', [min(onNear), max(onNear)], 'gamma_range', [min(gaNear), max(gaNear)], ...
    'corr_onset_gamma', rNear, 'principal_direction', dirVec(:)', 'principal_angle_deg', dirAngle);
verdictOut.parts = parts;
verdictOut.sensitivity = struct('condition_number', S.cond_number, 'svd', S.svd(:)', ...
    'corr_estimates_onset_gamma', corrOG, 'colnorm', S.colnorm(:)', ...
    'practical_rank', S.practical_rank);
verdictOut.sensitivity.collinear_sets = S.collinear_sets;   % assigned field-wise: a cell
                                                            % array in struct(...) would
                                                            % expand into a struct array
verdictOut.joint_constraint_text = jointTxt;
verdictOut.note = ['T12 audit reported the 7-D sensitivity at the LEGACY reference point ' ...
    '(condition number 233.9, collinear set {m0, rho, onset_doy, gamma, beta}); this run ' ...
    'recomputes it at the main report point, so the numbers differ by construction.'];
mwspm_writejson(fullfile(ctx.dataDir, 'verdict.json'), verdictOut);

summary = struct();
summary.experiment_id = ctx.expId;
summary.task = char(string(getdef(cfg, 'task', '-')));
summary.parent_fit_run_id = fitRun.runId;
summary.main_report_point = fit.main_report_point;
summary.deterministic = true;
summary.rng = 'none (deterministic optimisation; no random screening, no random source)';
summary.parts = parts;
summary.control = struct('J_reference_from_NM', bRef.J, 'J_parent_nll', Jref, 'dJ', dRef, ...
    'elapsed_seconds', tRef);
summary.verdict = verdict;
summary.onset_profile = verdictOut.onset_profile;
summary.gamma_profile = verdictOut.gamma_profile;
summary.joint_2d = verdictOut.joint_2d;
summary.multistart = verdictOut.multistart;
summary.sensitivity = verdictOut.sensitivity;
summary.elapsed = struct('multistart_s', tB, 'joint_2d_s', tD);
summary.honesty = cfg.honesty;
mwspm_writejson(fullfile(ctx.dataDir, 'summary.json'), summary);

out = struct();
out.status = 'SUCCESS';
out.parent_run_id = fitRun.runId;
out.verdict = verdict;
out.onset_ci95 = [pOn.lo, pOn.hi];
out.gamma_ci95 = [pGa.lo, pGa.hi];
out.joint_touches_boundary = touch;
out.projection_control_violated = projNarrower;
out.multistart_corr = rNear;
out.deterministic = true;
out.parts = parts;
out.outputs = {'sensitivity.json', 'profile_1d_onset.csv', 'profile_1d_gamma.csv', ...
    'joint_profile_2d.csv', 'verdict.json', 'summary.json'};
if doB
    out.outputs = [out.outputs(1), {'multistart_scatter.csv'}, out.outputs(2:end)];
end
out.note = ['T14.03. Identifiability of onset_doy and gamma at the main report point, lag d=1. ' ...
    'The verdict rule is pre-registered; the grid is never widened after seeing the result.'];
end

% ---------------------------------------------------------------------------
function [P, rowsOut] = t1403_profile1d(pre, mkC, baseOpt, free6, th, fname, grid, dFix)
rowsOut = {};
J = nan(numel(grid), 1);
opt = nan(numel(grid), 6);
for i = 1:numel(grid)
    fixed = struct(fname, grid(i));
    st = th; st.(fname) = grid(i);
    [b, ~] = mwspm_calib_fit(pre, mkC(free6, fixed), mergeOpt(baseOpt, struct('warm', {{st}})));
    J(i) = b.J;
    opt(i, :) = [b.theta.m0, b.theta.rho, b.theta.offset_doy, b.theta.mu_s, b.theta.beta, ...
        b.theta.(t1403_other(fname))];
    rowsOut{end+1, 1} = {grid(i), b.J, NaN, b.theta.(t1403_other(fname)), b.theta.m0, ...
        b.theta.rho, b.theta.offset_doy, b.theta.mu_s, b.theta.beta}; %#ok<AGROW>
end
dJv = J - min(J);
for i = 1:numel(grid)
    rowsOut{i}{3} = dJv(i);
end
P = struct();
P.grid = grid(:);
P.J = J;
P.dJ = dJv;
P.minAt = grid(find(J == min(J), 1));
P.minDJ = min(dJv);
P.opt = opt;
P.lo = NaN; P.hi = NaN; P.bounded = false;
ok = dJv <= 1.9207;
if any(ok)
    P.lo = min(grid(ok));
    P.hi = max(grid(ok));
    P.bounded = (min(grid(ok)) > min(grid) + 1e-9) && (max(grid(ok)) < max(grid) - 1e-9);
end
end

% ---------------------------------------------------------------------------
function f = t1403_other(fname)
if strcmp(fname, 'onset_doy')
    f = 'gamma';
else
    f = 'onset_doy';
end
end

% ---------------------------------------------------------------------------
function s = t1403_fmtInterval(lo, hi)
if ~isfinite(lo) || ~isfinite(hi)
    s = 'NA';
else
    s = sprintf('%.4f, %.4f', lo, hi);
end
end

% ---------------------------------------------------------------------------
function s = t1403_cellstr(c)
if isempty(c)
    s = '(none)';
    return
end
parts = cell(1, numel(c));
for i = 1:numel(c)
    if iscell(c{i})
        parts{i} = ['{' strjoin(cellstr(string(c{i})), ', ') '}'];
    else
        parts{i} = char(string(c{i}));
    end
end
s = strjoin(parts, ' ; ');
end

% ---------------------------------------------------------------------------
function r = t1403_pearson(x, y)
% Pearson correlation without the Statistics and Machine Learning Toolbox.
x = double(x(:));
y = double(y(:));
x = x - mean(x);
y = y - mean(y);
d = sqrt(sum(x .^ 2) * sum(y .^ 2));
if d <= 0
    r = NaN;
else
    r = sum(x .* y) / d;
end
end

% ---------------------------------------------------------------------------
function M = t1403_num(c)
% cell-of-1xk-cells -> numeric matrix.  cell2mat(c) alone fails here: c is a
% cell array whose elements are themselves cell arrays.
M = cell2mat(vertcat(c{:}));
end

% ---------------------------------------------------------------------------
function t1403_writecsv(path, hdr, M)
fid = fopen(path, 'w', 'n', 'UTF-8');
assert(fid > 0, 't14:io', 'cannot write %s', path);
fprintf(fid, '%s\n', strjoin(hdr, ','));
fmt = repmat('%.17g,', 1, size(M, 2));
fmt = [fmt(1:end-1) '\n'];
for i = 1:size(M, 1)
    fprintf(fid, fmt, M(i, :));
end
fclose(fid);
end

% ---------------------------------------------------------------------------
function o = mergeOpt(base, over)
o = base;
fn = fieldnames(over);
for i = 1:numel(fn)
    o.(fn{i}) = over.(fn{i});
end
end

% ---------------------------------------------------------------------------
function th = t1403_theta_from_fit(tf)
th = struct('m0', double(tf.m0), 'rho', double(tf.rho), 'onset_doy', double(tf.onset_doy), ...
    'offset_doy', double(tf.offset_doy), 'mu_s', double(tf.mu_s), 'gamma', double(tf.gamma), ...
    'beta', double(tf.beta), 'd_weeks', double(tf.d_weeks));
end

% ---------------------------------------------------------------------------
function P = t1403_find_parent(root, wantId, expId, relPath)
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
    try
        j = mwspm_readjson(f);
    catch
        continue
    end
    if ~isfield(j, 'provenance') || ~isfield(j.provenance, 'experiment_id') ...
            || ~strcmp(char(string(j.provenance.experiment_id)), expId)
        continue
    end
    if isempty(P) || strcmp(d(k).name, P.runId) > 0
        P = struct('runId', d(k).name, 'path', f);
    end
end
end

% ---------------------------------------------------------------------------
function s = yesNo(c)
if logical(c)
    s = 'YES';
else
    s = 'NO';
end
end

function s = t1403_yn(c, whenTrue, whenFalse)
if logical(c)
    s = whenTrue;
else
    s = whenFalse;
end
end

function v = getdef(s, f, dflt)
if isfield(s, f) && ~isempty(s.(f))
    v = s.(f);
else
    v = dflt;
end
end
