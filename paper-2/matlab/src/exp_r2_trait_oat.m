function payload = exp_r2_trait_oat(ctx)
%EXP_R2_TRAIT_OAT  Thirty one-at-a-time trait-coefficient perturbations on the
%   revised specification.
%
%   Replaces the initial-Euler numbers (reported totals 284 to 110 553) with the
%   same experimental design recomputed at a revised calibration point: the
%   fitted parameters are held fixed, only one trait coefficient changes at a
%   time (five traits x {c +/-10%, T0 +/-1 C, Tm +/-1 C} = 30 variants), and the
%   forward map is the stable kernel (daily ZOH, unclamped RK4 at 0.025 d) used
%   by the rest of the revision.
%
%   Reported totals are the kernel's expected weekly reports summed over the 47
%   weeks at the fixed reporting ratio.  tau* is recomputed with the perturbed
%   traits, because both R0(w) and the latent trajectory change.
%
%   One baseline variant (no perturbation) is included and must reproduce the
%   unperturbed fit; that is asserted, not assumed.

cfg = ctx.cfg;
if isfield(ctx, 'dryRun') && ctx.dryRun
    payload = struct('outputs', {{}}); payload.dry_run = true; return
end

diary(fullfile(ctx.runDir, 'logs', 'stdout.txt'));
diaryGuard = onCleanup(@() diary('off')); %#ok<NASGU>
snap = fullfile(ctx.dataDir, 'code_snapshot');
mkdir(snap);
copyfile(fullfile(ctx.root, 'src', 'matlab', '*.m'), snap);

D    = mwspm_load_data(ctx.root, ctx.modelCard);
pRef = mwspm_params(ctx.modelCard);
th0  = local_theta(ctx.root, cfg);

coef0 = mwspm_thermal_traits(20, 'raw').coef;
traits = cellstr(string(cfg.perturbations.traits));
relC   = double(cfg.perturbations.c_relative);
dT0    = double(cfg.perturbations.T0_shift_C);
dTm    = double(cfg.perturbations.Tm_shift_C);

variants = {'baseline', '-', '-', 0, struct()};
for i = 1:numel(traits)
    f = traits{i};
    tri = coef0.(f);
    for sgn = [1, -1]
        v = tri; v(1) = tri(1) * (1 + sgn * relC);
        variants(end+1, :) = {sprintf('%s_c_%+g%%', f, 100 * sgn * relC), f, 'c_relative', sgn * relC, local_ov(f, v)}; %#ok<AGROW>
    end
    for sgn = [1, -1]
        v = tri; v(2) = tri(2) + sgn * dT0;
        variants(end+1, :) = {sprintf('%s_T0_%+gC', f, sgn * dT0), f, 'T0_shift_C', sgn * dT0, local_ov(f, v)}; %#ok<AGROW>
    end
    for sgn = [1, -1]
        v = tri; v(3) = tri(3) + sgn * dTm;
        variants(end+1, :) = {sprintf('%s_Tm_%+gC', f, sgn * dTm), f, 'Tm_shift_C', sgn * dTm, local_ov(f, v)}; %#ok<AGROW>
    end
end

rows = table();
base = struct();

for i = 1:size(variants, 1)
    label = variants{i, 1};
    ov    = variants{i, 5};

    prepOpts = struct('mu_variant', cfg.mu_variant);
    trOpts   = struct();
    if ~isempty(fieldnames(ov))
        prepOpts.coef = ov;
        trOpts.coef   = ov;
    end
    pre = mwspm_calib_prepare(D, ctx.modelCard, prepOpts);
    tr  = mwspm_thermal_traits(D.T_daily, 'raw', trOpts);

    kern = mwspm_calib_kernel(pre, th0, cfg.solver);
    [tauW, T0at, cumAt] = local_tau(pre, D, pRef, th0, th0.rho, cfg.solver, trOpts, tr);

    R0daily = mwspm_r0(D.T_daily, th0.m0 * local_season(pre, th0), ...
        th0.mu_s * tr.mu_raw, pRef.gamma_h, trOpts);
    R0sup = max(R0daily);

    row = table({label}, {variants{i,2}}, {variants{i,3}}, variants{i,4}, ...
        kern.reported_total, kern.true_total, tauW, T0at, cumAt, R0sup, ...
        kern.peak_week, ...
        'VariableNames', {'variant','trait','kind','delta','reported_total', ...
        'true_total','tau_star','T0_at_tau','cum_reports_at_tau','sup_R0','peak_week'});
    rows = [rows; row]; %#ok<AGROW>

    if strcmp(label, 'baseline')
        base = struct('reported_total', kern.reported_total, ...
            'true_total', kern.true_total, 'tau_star', tauW, 'sup_R0', R0sup);
    end
    fprintf('%-16s %-8s reported=%12.4f true=%12.4f tau*=%s supR0=%.6f\n', ...
        label, variants{i,2}, kern.reported_total, kern.true_total, ...
        local_num(tauW), R0sup);
end

b = rows(strcmp(rows.variant, 'baseline'), :);
assert(abs(b.reported_total - base.reported_total) < 1e-9, 'R2oat:baseline', ...
    'baseline row does not reproduce the unperturbed fit');

res = struct();
res.point = th0;
res.solver = cfg.solver;
res.mu_variant = cfg.mu_variant;
res.baseline = base;
res.n_variants = height(rows) - 1;
res.reported_total_min = min(rows.reported_total);
res.reported_total_max = max(rows.reported_total);
res.true_total_min = min(rows.true_total);
res.true_total_max = max(rows.true_total);
res.tau_star_unique = unique(rows.tau_star(isfinite(rows.tau_star)))';
res.sup_R0_min = min(rows.sup_R0);
res.sup_R0_max = max(rows.sup_R0);
res.peak_week_unique = unique(rows.peak_week)';
mwspm_writejson(fullfile(ctx.dataDir, 'oat_summary.json'), res);
writetable(rows, fullfile(ctx.dataDir, 'oat_variants.csv'));

fprintf(['OAT: reported total %.1f - %.1f (baseline %.1f); latent %.1f - %.1f; ' ...
    'tau* values %s\n'], res.reported_total_min, res.reported_total_max, ...
    base.reported_total, res.true_total_min, res.true_total_max, ...
    mat2str(res.tau_star_unique));

payload = res;
payload.outputs = cfg.outputs;
end

% ---------------------------------------------------------------------------
function th = local_theta(root, cfg)
J = mwspm_readjson(fullfile(root, 'results', 'runs', cfg.point_run, 'data', 'stable_fit.json'));
th = J.(sprintf('lag%d', cfg.lag)).best.theta;
end

% ---------------------------------------------------------------------------
function ov = local_ov(field, triple)
ov = struct();
ov.(field) = triple(:)';
end

% ---------------------------------------------------------------------------
function s = local_season(pre, th)
sp = struct('onset_doy', th.onset_doy, 'offset_doy', th.offset_doy, ...
    'kappa_days', pre.kappa_days, 'gamma', th.gamma, 'beta', th.beta, ...
    's_floor', pre.s_floor, 'min_window_days', pre.min_window_days, ...
    'normalize', true, 'norm_factor', NaN);
s = mwspm_season(pre.doy, sp);
end

% ---------------------------------------------------------------------------
function [tauStar, T0atTau, cumAtTau] = local_tau(pre, D, pRef, th, rhoUsed, solverOpts, trOpts, tr)
nW = pre.n_weeks;
s = local_season(pre, th);
R0daily = mwspm_r0(D.T_daily, th.m0 * s, th.mu_s * tr.mu_raw, pRef.gamma_h, trOpts);
wk = floor(((1:numel(R0daily))' - 1) / 7) + 1;
wk = min(max(wk, 1), nW);
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
function s = local_num(v)
if isfinite(v), s = sprintf('%d', v); else, s = 'NA'; end
end
