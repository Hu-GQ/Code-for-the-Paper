function out = exp_t14_trait_robustness(ctx)
%EXP_T14_TRAIT_ROBUSTNESS  T14.02 -- robustness of R0(t) and tau* to the thermal
%   TRAIT PARAMETERISATION (EXPERIMENT_PLAN.md Sec.5-4).
%
%   WHAT IS MEASURED, AND WHAT IS NOT
%   ---------------------------------
%   Sec.5-4 asks for "Mordecai curves vs another source".  There is NO second
%   temperature-dependent source in this project: the initial version
%   (archive/text1/matlab/thermal_traits.m) uses the SAME Mordecai 2017 Table A
%   coefficients, and Andraud 2012 Table 2 gives temperature-INDEPENDENT ranges
%   only.  The source-swap comparison is therefore NOT PERFORMED and is
%   registered as such in summary.json and in the run README.
%
%   What IS measured here, with a single trait definition (mwspm_thermal_traits,
%   opt-in coefficient override):
%     (1) one-at-a-time perturbation of the coefficient triples [c,T0,Tm] of
%         a, b, c, sigma_v, lf: c scaled by +/-10%, T0 and Tm shifted by +/-1 C
%         (5 traits x 6 = 30 variants + baseline), with the FITTED parameters
%         held fixed (no refit);
%     (2) for every variant, the weekly R0(t) diagnostics, tau* recomputed with
%         exactly the gate-G1 caliber (trait-modified kernel -> T0(w)=rho/logR0),
%         and the reported total;
%     (3) the R0(T) anchors implied by this implementation at the fitted mu_s,
%         listed against the anchors PUBLISHED by Mordecai 2017 for
%         Ae. albopictus (peak 26.4 C [25.2,27.4]; zero limits 16.2 [13.2,19.9]
%         and 31.6 [29.4,33.7]);
%     (4) an external range screen of the effective constants against Andraud
%         2012 Table 2 (adult lifespan 4-50 d, EIP 8-12 d, biting 0.3-1 /day,
%         beta_h 0.1-0.75, beta_v 0.5-1).
%
%   The anchors in (3) are reported as DEVIATIONS, not as a pass/fail test: our
%   R0 carries the mu_s scaling and our closed form differs from the source's
%   (theirs contains mu^3 and the aquatic traits EFD/pEA/MDR).

logFid = ctx.logFid;
root   = ctx.root;
D      = mwspm_load_data(root, ctx.modelCard);
cfg    = ctx.cfg;
modelCard = ctx.modelCard;
pRef   = mwspm_params(modelCard);

fprintf(logFid, '=== exp_t14_trait_robustness (T14.02, Sec.5-4) ===\n');
fprintf(logFid, 'claims: %s\n', strjoin(cellstr(string(getdef(cfg, 'claims', {}))), ', '));

% ------------------------------------------------------------- parent fit ---
fitRun = t14_find_parent(root, getdef(cfg, 'parent_run_id', []), 'exp_c4_rho_fit', 'data/fit_best.json');
assert(~isempty(fitRun), 't14:noParent', ...
    ['no SUCCESS run of exp_c4_rho_fit with data/fit_best.json found under results/runs/. ' ...
     'This experiment freezes the fitted main report point and may not invent one.']);
fit = mwspm_readjson(fitRun.path);
th  = t14_theta_from_fit(fit.main_report_point.theta);
fprintf(logFid, 'parent fit run : %s\n', fitRun.runId);
fprintf(logFid, 'theta          : m0=%.6f rho=%.6f mu_s=%.6f onset=%.3f offset=%.3f gamma=%.4f beta=%.4f d=%d\n', ...
    th.m0, th.rho, th.mu_s, th.onset_doy, th.offset_doy, th.gamma, th.beta, th.d_weeks);

% ------------------------------------------------------------- dry-run ------
PD = cfg.parameter_domain;
if isfield(ctx, 'dryRun') && ctx.dryRun
    out = struct();
    out.dry_run = true;
    out.parent_run_id = fitRun.runId;
    out.validated = struct('n_traits', numel(PD.perturbations.traits), ...
        'n_variants', numel(PD.perturbations.traits) * 6 + 1, ...
        'r0T_grid_points', numel(PD.r0T_grid.from:PD.r0T_grid.step:PD.r0T_grid.to));
    fprintf(logFid, '[dry-run] parent theta and grids readable; no perturbation computed\n');
    return
end

% ------------------------------------------------------------- grids --------
Tgrid = (double(PD.r0T_grid.from) : double(PD.r0T_grid.step) : double(PD.r0T_grid.to))';
sp = struct('onset_doy', th.onset_doy, 'offset_doy', th.offset_doy, ...
    'kappa_days', pRef.kappa_days, 'gamma', th.gamma, 'beta', th.beta, ...
    's_floor', pRef.s_floor, 'min_window_days', pRef.min_window_days, ...
    'normalize', true, 'norm_factor', NaN);
s = mwspm_season(D.doy, sp);
nW = ceil(numel(D.T_daily) / 7);
wk = floor(((1:numel(D.T_daily))' - 1) / 7) + 1;
wk = min(max(wk, 1), nW);

traits = cellstr(string(PD.perturbations.traits));
relC   = double(PD.perturbations.c_relative);
dT0    = double(PD.perturbations.T0_shift_C);
dTm    = double(PD.perturbations.Tm_shift_C);

% ------------------------------------------------------------- variants -----
% Each row: label, trait, kind, delta, override-coef struct ([] = baseline)
coef0 = mwspm_thermal_traits(20, 'raw').coef;
variants = {};
variants(end+1, :) = {'baseline', '-', '-', 0, []};
for i = 1:numel(traits)
    f = traits{i};
    tri = coef0.(f);
    % c scaled by (1 +/- rel): keeps the sign of the published coefficient
    for sgn = [1, -1]
        if sgn > 0, sgnTxt = '+'; else, sgnTxt = '-'; end
        v = tri; v(1) = tri(1) * (1 + sgn * relC);
        variants(end+1, :) = {sprintf('%s_c_%s%g%%', f, sgnTxt, 100 * relC), ...
            f, 'c_relative', sgn * relC, t14_one(f, v)}; %#ok<AGROW>
    end
    for sgn = [1, -1]
        v = tri; v(2) = tri(2) + sgn * dT0;
        variants(end+1, :) = {sprintf('%s_T0_%+gC', f, sgn * dT0), f, 'T0_shift_C', sgn * dT0, t14_one(f, v)}; %#ok<AGROW>
    end
    for sgn = [1, -1]
        v = tri; v(3) = tri(3) + sgn * dTm;
        variants(end+1, :) = {sprintf('%s_Tm_%+gC', f, sgn * dTm), f, 'Tm_shift_C', sgn * dTm, t14_one(f, v)}; %#ok<AGROW>
    end
end
nVar = size(variants, 1);
fprintf(logFid, '\n--- %d variants (1 baseline + %d OAT perturbations) ---\n', nVar, nVar - 1);

% ------------------------------------------------------------- run ----------
rows = {};
curves = nan(numel(Tgrid), nVar + 1);
curves(:, 1) = Tgrid;
anchorRows = {};
base = struct();
% numeric accumulators (rows{} holds formatted strings for the CSV, so the
% statistics must NOT be taken from it: 'NA' is a char and breaks cell2mat)
tauNum  = nan(nVar, 1);
dTauNum = nan(nVar, 1);
dPeakNum = nan(nVar, 1);
relR0Num = nan(nVar, 1);
R0maxNum = nan(nVar, 1);
peakWkNum = nan(nVar, 1);

for i = 1:nVar
    label = variants{i, 1};
    coefOv = variants{i, 5};
    trOpts = struct();
    if ~isempty(coefOv)
        trOpts.coef = coefOv;
    end

    tr    = mwspm_thermal_traits(D.T_daily, 'raw', trOpts);
    muEff = th.mu_s * tr.mu_raw;
    R0d   = mwspm_r0(D.T_daily, th.m0 * s, muEff, pRef.gamma_h, trOpts);
    R0w   = accumarray(wk, R0d, [nW 1], @mean, NaN);

    gt1 = find(R0d > 1);
    if isempty(gt1)
        w1 = NaN; w2 = NaN;
    else
        w1 = wk(gt1(1)); w2 = wk(gt1(end));
    end
    [mx, iMx] = max(R0d);

    preOpts = struct('mu_variant', 'day');
    if ~isempty(coefOv)
        preOpts.coef = coefOv;
    end
    pre  = mwspm_calib_prepare(D, modelCard, preOpts);
    kern = mwspm_calib_kernel(pre, th, struct());
    cumDet = cumsum(kern.Yhat(:));
    T0w = nan(nW, 1);
    ok = R0w > 1;
    T0w(ok) = th.rho ./ log(R0w(ok));
    cross = find(isfinite(T0w) & cumDet >= T0w, 1);
    if isempty(cross)
        tauStar = NaN; cumTau = NaN; T0Tau = NaN;
    else
        tauStar = cross; cumTau = cumDet(cross); T0Tau = T0w(cross);
    end

    trT = mwspm_thermal_traits(Tgrid, 'raw', trOpts);
    R0T = mwspm_r0(Tgrid, th.m0, th.mu_s * trT.mu_raw, pRef.gamma_h, trOpts);
    pos = find(R0T > 0);
    if isempty(pos)
        loT = NaN; hiT = NaN; pkT = NaN; pkV = NaN;
    else
        loT = Tgrid(pos(1)); hiT = Tgrid(pos(end));
        [pkV, ip] = max(R0T); pkT = Tgrid(ip);
    end
    curves(:, i + 1) = R0T;

    if i == 1
        base = struct('tau_star', tauStar, 'R0_max', mx, 'peak_week', wk(iMx), ...
            'R0T_peak_T', pkT, 'R0T_lo_T', loT, 'R0T_hi_T', hiT, 'R0_curve', R0w);
    end

    rows{end+1, 1} = {label, variants{i, 2}, variants{i, 3}, variants{i, 4}, ...
        mx, wk(iMx), w1, w2, numel(gt1), tauStar, T0Tau, cumTau, ...
        kern.reported_total, kern.peak_week, ...
        tauStar - base.tau_star, wk(iMx) - base.peak_week, ...
        max(abs(R0w - base.R0_curve)) / max(base.R0_curve)}; %#ok<AGROW>
    tauNum(i)   = tauStar;
    dTauNum(i)  = tauStar - base.tau_star;
    dPeakNum(i) = wk(iMx) - base.peak_week;
    relR0Num(i) = max(abs(R0w - base.R0_curve)) / max(base.R0_curve);
    R0maxNum(i) = mx;
    peakWkNum(i) = wk(iMx);
    anchorRows{end+1, 1} = {label, pkT, pkV, loT, hiT, ...
        pkT - 26.4, loT - 16.2, hiT - 31.6}; %#ok<AGROW>

    fprintf(logFid, '  %-14s R0max=%.4f peak_wk=%2d gt1=[%s,%s] days>1=%3d tau*=%s cum=%s rep_total=%.3f (dtau=%+.0f wk, dpeak=%+.0f wk, maxrel dR0w=%.4f)\n', ...
        label, mx, wk(iMx), num2strOrNA(w1), num2strOrNA(w2), numel(gt1), ...
        num2strOrNA(tauStar), numOrNA(cumTau), kern.reported_total, ...
        tauStar - base.tau_star, wk(iMx) - base.peak_week, ...
        max(abs(R0w - base.R0_curve)) / max(base.R0_curve));
end

% ------------------------------------------------- external range screen ----
Ts = mean(D.T_daily(s >= 0.5));
nScr = sum(s >= 0.5);
trS = mwspm_thermal_traits(Ts, 'raw');
lifeEff = 1 / max(th.mu_s * trS.mu_raw, realmin);
EIPeff  = 1 / max(trS.sigma_v, realmin);
Ranges = cfg.andraud2012_ranges;
scr = {};
scr(end+1, :) = {'adult_lifespan_days', lifeEff, Ranges.adult_lifespan_days(1), Ranges.adult_lifespan_days(2)};
scr(end+1, :) = {'EIP_days', EIPeff, Ranges.EIP_days(1), Ranges.EIP_days(2)};
scr(end+1, :) = {'biting_rate_per_day', trS.a, Ranges.biting_rate_per_day(1), Ranges.biting_rate_per_day(2)};
scr(end+1, :) = {'beta_mosquito_to_human', trS.b, Ranges.beta_mosquito_to_human(1), Ranges.beta_mosquito_to_human(2)};
scr(end+1, :) = {'beta_human_to_mosquito', trS.c, Ranges.beta_human_to_mosquito(1), Ranges.beta_human_to_mosquito(2)};
fprintf(logFid, '\n--- Andraud 2012 Table 2 range screen at T=%.3f C (mean over %d days with s>=0.5) ---\n', Ts, nScr);
for i = 1:size(scr, 1)
    inR = (scr{i, 2} >= scr{i, 3}) && (scr{i, 2} <= scr{i, 4});
    fprintf(logFid, '  %-26s %10.5f  in [%.4g, %.4g] : %s\n', scr{i, 1}, scr{i, 2}, scr{i, 3}, scr{i, 4}, yesNo(inR));
end

% ------------------------------------------------- R0(T) anchor comparison --
fprintf(logFid, '\n--- R0(T) anchors of THIS implementation vs Mordecai 2017 published (Ae. albopictus, constant T) ---\n');
fprintf(logFid, '  ours    : peak %.2f C | lower zero %.2f C | upper zero %.2f C (at mu_s = %.4f)\n', ...
    base.R0T_peak_T, base.R0T_lo_T, base.R0T_hi_T, th.mu_s);
A = cfg.mordecai2017_published_anchors;
fprintf(logFid, '  published: peak %.1f [%.1f, %.1f] | lower %.1f [%.1f, %.1f] | upper %.1f [%.1f, %.1f]\n', ...
    A.peak_temperature_C.mean, A.peak_temperature_C.ci95(1), A.peak_temperature_C.ci95(2), ...
    A.lower_zero_temperature_C.mean, A.lower_zero_temperature_C.ci95(1), A.lower_zero_temperature_C.ci95(2), ...
    A.upper_zero_temperature_C.mean, A.upper_zero_temperature_C.ci95(1), A.upper_zero_temperature_C.ci95(2));

% --------------------------------------------------------------- tau* spread
taus = tauNum(isfinite(tauNum));
tauSpread = max(taus) - min(taus);
dtau  = dTauNum(2:end);
dpeak = dPeakNum(2:end);
relR0 = relR0Num(2:end);

% --------------------------------------------------------------- outputs ----
hdr = {'variant', 'trait', 'kind', 'delta', 'R0_max', 'R0_peak_week', ...
    'R0_gt1_first_week', 'R0_gt1_last_week', 'R0_days_gt1', 'tau_star', 'T0_at_tau', ...
    'cum_detected_at_tau', 'reported_total', 'kernel_peak_week', ...
    'd_tau_vs_baseline_weeks', 'd_R0_peak_week_vs_baseline', 'max_rel_change_R0_weekly'};
fid = fopen(fullfile(ctx.dataDir, 'trait_variants.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '%s\n', strjoin(hdr, ','));
for i = 1:numel(rows)
    rw = rows{i};
    fprintf(fid, '%s,%s,%s,%.17g,%.17g,%d,%s,%s,%d,%s,%s,%s,%.17g,%d,%s,%s,%.17g\n', ...
        rw{1}, rw{2}, rw{3}, rw{4}, rw{5}, rw{6}, num2strOrNA(rw{7}), num2strOrNA(rw{8}), rw{9}, ...
        num2strOrNA(rw{10}), num2strOrNA(rw{11}), num2strOrNA(rw{12}), rw{13}, rw{14}, ...
        num2strOrNA(rw{15}), num2strOrNA(rw{16}), rw{17});
end
fclose(fid);

fid = fopen(fullfile(ctx.dataDir, 'r0T_anchors.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, 'variant,R0T_peak_T_C,R0T_peak_value,R0T_lower_zero_T_C,R0T_upper_zero_T_C,d_peak_vs_published_C,d_lower_vs_published_C,d_upper_vs_published_C\n');
for i = 1:size(anchorRows, 1)
    a = anchorRows{i};
    fprintf(fid, '%s,%.17g,%.17g,%.17g,%.17g,%+.17g,%+.17g,%+.17g\n', ...
        a{1}, a{2}, a{3}, a{4}, a{5}, a{6}, a{7}, a{8});
end
fclose(fid);

mwspm_writecsv(fullfile(ctx.dataDir, 'r0T_curves.csv'), ...
    [{'T_C'}, cellfun(@(x) ['R0T_' x], variants(:, 1)', 'UniformOutput', false)], curves);

fid = fopen(fullfile(ctx.dataDir, 'andraud_range_screen.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, 'quantity,value,range_low,range_high,inside_range\n');
for i = 1:size(scr, 1)
    inR = (scr{i, 2} >= scr{i, 3}) && (scr{i, 2} <= scr{i, 4});
    fprintf(fid, '%s,%.17g,%.17g,%.17g,%d\n', scr{i, 1}, scr{i, 2}, scr{i, 3}, scr{i, 4}, double(inR));
end
fprintf(fid, '#screen_temperature_C,%.17g,#days_s_ge_0.5,%d\n', Ts, nScr);
fclose(fid);

summary = struct();
summary.experiment_id = ctx.expId;
summary.task = char(string(getdef(cfg, 'task', '-')));
summary.parent_fit_run_id = fitRun.runId;
summary.main_report_point = fit.main_report_point;
summary.deterministic = true;
summary.rng = 'none (deterministic recomputation; no random source)';
summary.n_variants = nVar;
summary.perturbations = PD.perturbations;
rw0 = rows{1};
summary.baseline = struct('tau_star', base.tau_star, 'R0_max', base.R0_max, ...
    'R0_peak_week', base.peak_week, 'R0T_peak_T', base.R0T_peak_T, ...
    'R0T_lower_zero_T', base.R0T_lo_T, 'R0T_upper_zero_T', base.R0T_hi_T, ...
    'reported_total', rw0{13}, 'kernel_peak_week', rw0{14});
summary.tau_star_spread_over_variants_weeks = tauSpread;
summary.tau_star_min = min(taus);
summary.tau_star_max = max(taus);
summary.max_abs_d_tau_weeks = max(abs(dtau));
summary.n_variants_changing_tau = sum(dtau ~= 0);
summary.max_abs_d_R0_peak_week = max(abs(dpeak));
summary.max_rel_change_R0_weekly = max(relR0);
summary.andraud_screen = struct('temperature_C', Ts, 'n_days', nScr, ...
    'quantity', {scr(:, 1)'}, 'value', cell2mat(scr(:, 2)), ...
    'range_low', cell2mat(scr(:, 3)), 'range_high', cell2mat(scr(:, 4)), ...
    'inside_range', cell2mat(scr(:, 2)) >= cell2mat(scr(:, 3)) & cell2mat(scr(:, 2)) <= cell2mat(scr(:, 4)), ...
    'all_inside', all(cell2mat(scr(:, 2)) >= cell2mat(scr(:, 3)) & cell2mat(scr(:, 2)) <= cell2mat(scr(:, 4))));
summary.mordecai_anchor_deviation = struct('peak_T_ours', base.R0T_peak_T, ...
    'peak_T_published_mean', A.peak_temperature_C.mean, ...
    'peak_T_deviation', base.R0T_peak_T - A.peak_temperature_C.mean, ...
    'lower_zero_ours', base.R0T_lo_T, 'lower_zero_published_mean', A.lower_zero_temperature_C.mean, ...
    'lower_zero_deviation', base.R0T_lo_T - A.lower_zero_temperature_C.mean, ...
    'upper_zero_ours', base.R0T_hi_T, 'upper_zero_published_mean', A.upper_zero_temperature_C.mean, ...
    'upper_zero_deviation', base.R0T_hi_T - A.upper_zero_temperature_C.mean, ...
    'note', ['Deviations are reported as MEASURED, not as a pass/fail test: our R0 carries the ' ...
             'mu_s scaling and our closed form omits the aquatic traits (EFD/pEA/MDR) and the ' ...
             'mu^3 term of the source R0.']);
summary.source_swap_status = cfg.source_swap_status;
summary.honesty = cfg.honesty;
summary.value_type_note = ['Every number here is a deterministic_point: the same fitted parameters ' ...
    'were re-evaluated under a changed trait curve. No interval is produced: this is an OAT ' ...
    'perturbation grid, NOT a joint uncertainty region and NOT a replicate distribution.'];
mwspm_writejson(fullfile(ctx.dataDir, 'summary.json'), summary);

out = struct();
out.status = 'SUCCESS';
out.parent_run_id = fitRun.runId;
out.n_variants = nVar;
out.tau_star_spread_over_variants_weeks = tauSpread;
out.max_abs_d_tau_weeks = max(abs(dtau));
out.source_swap_status = cfg.source_swap_status.verdict;
out.deterministic = true;
out.outputs = {'trait_variants.csv', 'r0T_curves.csv', 'r0T_anchors.csv', 'andraud_range_screen.csv', 'summary.json'};
out.note = ['T14.02. Coefficient perturbation, NOT a source swap (registered). ' ...
    'tau* is recomputed with the gate-G1 caliber after changing only the trait curves; ' ...
    'the fitted parameters are never re-optimised.'];
end

% ---------------------------------------------------------------------------
function ov = t14_one(field, triple)
ov = struct();
ov.(field) = triple(:)';
end

% ---------------------------------------------------------------------------
function th = t14_theta_from_fit(tf)
th = struct('m0', double(tf.m0), 'rho', double(tf.rho), 'onset_doy', double(tf.onset_doy), ...
    'offset_doy', double(tf.offset_doy), 'mu_s', double(tf.mu_s), 'gamma', double(tf.gamma), ...
    'beta', double(tf.beta), 'd_weeks', double(tf.d_weeks));
end

% ---------------------------------------------------------------------------
function P = t14_find_parent(root, wantId, expId, relPath)
% Latest results/runs/<run_id>/<relPath> whose run is SUCCESS for <expId>.
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
function s = num2strOrNA(v)
if isfinite(double(v))
    s = sprintf('%.6g', double(v));
else
    s = 'NA';
end
end

function s = numOrNA(v)
if isfinite(double(v))
    s = sprintf('%.17g', double(v));
else
    s = 'NA';
end
end

function s = yesNo(c)
if logical(c)
    s = 'YES';
else
    s = 'NO';
end
end

function v = getdef(s, f, dflt)
if isfield(s, f) && ~isempty(s.(f))
    v = s.(f);
else
    v = dflt;
end
end
