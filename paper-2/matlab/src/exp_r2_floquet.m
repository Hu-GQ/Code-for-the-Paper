function payload = exp_r2_floquet(ctx)
%EXP_R2_FLOQUET  Three periodic criteria at the revised calibration points.
%
%   The archived cor:gzhou table belongs to the initial Euler fit.  This run
%   recomputes the three criteria (pointwise, envelope, Floquet) at the revised
%   conditional calibration point (fixed rho = 0.8, lag 1) and at the lowest
%   objective candidate of the revised profile (rho = 0.45, lag 1), using the
%   same 329-day zero-order-held artificial periodic repetition that the theory
%   source now declares.
%
%   A regression block first re-derives the archived T12.01 numbers from
%   theory/derivations/gzhou_t1201_daily_series.csv and asserts that
%   mwspm_floquet_metrics reproduces them; if the shared implementation ever
%   drifts, this run fails instead of publishing new numbers.

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

% ---------------------------------------------------------------- regression
serPath = fullfile(ctx.root, cfg.regression.series);
tab = readtable(serPath);
S0 = struct('m', tab.m, 'a', tab.a, 'b', tab.b, 'c', tab.c, ...
    'mu', tab.mu_eff, 'sv', tab.sigma_v);
if ismember('R0', tab.Properties.VariableNames)
    S0.R0_exported = tab.R0;
end
reg = mwspm_floquet_metrics(S0, struct('doLinear', true));
exp = cfg.regression.expected;
tol = cfg.regression.relative_tolerance;
chk = @(name, got, want) assert(abs(got - want) <= tol * max(1, abs(want)), ...
    'R2floquet:regression', '%s drifted: got %.12g want %.12g', name, got, want);
chk('sup_R0', reg.sup_R0, exp.sup_R0);
chk('R0_envelope', reg.R0_envelope, exp.R0_envelope);
chk('zoh_rho', reg.zoh_rho, exp.zoh_rho);
chk('zoh_kappa_F', reg.zoh_kappa_F, exp.zoh_kappa_F);
fprintf('regression OK: supR0=%.9g env=%.9g zoh_rho=%.6g kappaF=%.9g\n', ...
    reg.sup_R0, reg.R0_envelope, reg.zoh_rho, reg.zoh_kappa_F);

% ------------------------------------------------------------------- points
names = fieldnames(cfg.points);
rows = table();
res = struct();
res.regression = reg;
res.points = struct();

for i = 1:numel(names)
    nm = names{i};
    spec = cfg.points.(nm);
    th = local_theta(ctx.root, spec);

    sp = struct('onset_doy', th.onset_doy, 'offset_doy', th.offset_doy, ...
        'kappa_days', pRef.kappa_days, 'gamma', th.gamma, 'beta', th.beta, ...
        's_floor', pRef.s_floor, 'min_window_days', pRef.min_window_days, ...
        'normalize', true, 'norm_factor', NaN);
    s = mwspm_season(D.doy, sp);
    tr = mwspm_thermal_traits(D.T_daily, cfg.mu_variant);

    S = struct('m', th.m0 * s, 'a', tr.a, 'b', tr.b, 'c', tr.c, ...
        'mu', th.mu_s * tr.mu_raw, 'sv', tr.sigma_v);
    m = mwspm_floquet_metrics(S);

    res.points.(nm) = m;
    res.points.(nm).theta = th;
    res.points.(nm).supply_caliber = 'N_v(t) = m0 N_h s(t); mu_eff = mu_s mu_v_raw(T)';

    rows = [rows; table({nm}, th.rho, m.sup_R0, m.peak_week, m.R0_envelope, ...
        m.kappa_star, m.zoh_rho, m.zoh_kappa_F, m.safety_factor, ...
        'VariableNames', {'point','rho','sup_R0','peak_week','R0_envelope', ...
        'kappa_star','zoh_rho','zoh_kappa_F','safety_factor'})]; %#ok<AGROW>

    fprintf(['point %-14s rho=%.2f supR0=%.9g (wk %d) env=%.9g kappa*=%.6g ' ...
        'zoh_rho=%.6g kappa_F=%.9g\n'], nm, th.rho, m.sup_R0, m.peak_week, ...
        m.R0_envelope, m.kappa_star, m.zoh_rho, m.zoh_kappa_F);
end

writetable(rows, fullfile(ctx.dataDir, 'floquet_revised.csv'));
mwspm_writejson(fullfile(ctx.dataDir, 'floquet_revised.json'), res);

payload = struct();
payload.summary = rows;
payload.outputs = cfg.outputs;
end

% ---------------------------------------------------------------------------
function th = local_theta(root, spec)
runDir = fullfile(root, 'results', 'runs', spec.run);
if isfield(spec, 'key')
    J = mwspm_readjson(fullfile(runDir, spec.source));
    th = J.(spec.key).best.theta;
    return
end
% profile CSV row
S = mwspm_readcsv(fullfile(runDir, spec.source));
getcol = @(nm) cellfun(@(r) str2double(r{strcmp(S.names, nm)}), S.rows);
lag = getcol('lag');
rho = getcol('rho');
k = find(lag == spec.lag & abs(rho - spec.rho) < 1e-9, 1);
assert(~isempty(k), 'R2floquet:row', 'no profile row for lag %d rho %g', spec.lag, spec.rho);
th = struct();
th.rho = rho(k);
th.d_weeks = lag(k);
for f = {'m0','mu_s','onset_doy','offset_doy','gamma','beta'}
    th.(f{1}) = getcol(f{1});
    th.(f{1}) = th.(f{1})(k);
end
end
