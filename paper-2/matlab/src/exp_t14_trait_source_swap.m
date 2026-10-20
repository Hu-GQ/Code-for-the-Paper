function out = exp_t14_trait_source_swap(ctx)
%EXP_T14_TRAIT_SOURCE_SWAP  T14.04 -- swap the WHOLE temperature-dependent trait
%   set, within one source: Mordecai 2017 S2 Text Table A (Ae. albopictus, the
%   set this project currently uses) vs Table B (Ae. aegypti).
%
%   WHY THIS EXISTS
%   ---------------
%   EXPERIMENT_PLAN Sec.5-4 asks for "Mordecai curves vs another source".  A
%   second PUBLICATION with trait-level coefficients has not been obtained, so
%   that part stays NOT PERFORMED.  What CAN be done with the source we already
%   have is a species swap: both species' fitted triples [c,T0,Tm] are published
%   in the same S2 Text, the two sets differ substantially, and both species
%   co-occur in Guangzhou.  This is registered as a SAME-SOURCE, DIFFERENT-SPECIES
%   swap and must never be described as "another source".
%
%   TWO CALIBERS PER TRAIT SET (both reported)
%   ------------------------------------------
%     (i)  FIXED parameters: the T12.01 main report point is held fixed and only
%          the trait curves change => isolates the effect of the curves;
%     (ii) REFIT: the seven parameters are re-optimised under the swapped curves
%          => answers "does the fit still stand?".
%   For both, tau* is recomputed with the gate-G1 caliber (trait-modified kernel,
%   T0(w) = rho/log R0(w)) and the same 47-week window.
%
%   ANALYTIC CROSS-CHECK (asserted in code)
%   ---------------------------------------
%   Table A's five triples must equal the values this project already uses
%   (data/metadata/parameters.csv = mwspm_thermal_traits defaults).  If they do
%   not, the SI parse or the mapping is wrong and the run must fail rather than
%   produce numbers.

logFid = ctx.logFid;
root   = ctx.root;
cfg    = ctx.cfg;
modelCard = ctx.modelCard;

fprintf(logFid, '=== exp_t14_trait_source_swap (T14.04, Sec.5-4 species swap) ===\n');
fprintf(logFid, 'claims: %s\n', strjoin(cellstr(string(getdef(cfg, 'claims', {}))), ', '));

D    = mwspm_load_data(root, modelCard);
pRef = mwspm_params(modelCard);

fitRun = t1404_find_parent(root, getdef(cfg, 'parent_run_id', []), 'exp_c4_rho_fit', 'data/fit_best.json');
assert(~isempty(fitRun), 't14:noParent', ...
    'no SUCCESS run of exp_c4_rho_fit with data/fit_best.json found under results/runs/');
fit = mwspm_readjson(fitRun.path);
th0 = t1404_theta_from_fit(fit.main_report_point.theta);
fprintf(logFid, 'parent fit run : %s\n', fitRun.runId);
fprintf(logFid, 'reference point: m0=%.6f rho=%.6f onset=%.4f offset=%.4f mu_s=%.5f gamma=%.5f beta=%.5f d=%d\n', ...
    th0.m0, th0.rho, th0.onset_doy, th0.offset_doy, th0.mu_s, th0.gamma, th0.beta, th0.d_weeks);

PD = cfg.parameter_domain;
sets = PD.trait_sets;
if isfield(PD, 'trait_set_order') && ~isempty(PD.trait_set_order)
    setNames = cellstr(string(PD.trait_set_order(:)));   % EXPLICIT order: baseline first
else
    setNames = fieldnames(sets);                          % fallback (alphabetical!)
end
for i = 1:numel(setNames)
    assert(isfield(sets, setNames{i}), 't14:badConfig', 'trait set "%s" is not defined', setNames{i});
end

% ---- analytic cross-check: the FIRST set must reproduce the project defaults -
% (shape-agnostic: mwspm_readjson may return a column where the default is a row)
coefDefault = mwspm_thermal_traits(20, 'raw').coef;
traitFields = {'a', 'b', 'c', 'sigma_v', 'lf'};
sA = sets.(setNames{1});
for i = 1:numel(traitFields)
    f = traitFields{i};
    assert(isfield(sA, f), 't14:badConfig', 'trait set 1 is missing "%s"', f);
    fromFile = double(sA.(f));
    fromFile = fromFile(:)';
    fromCode = double(coefDefault.(f));
    fromCode = fromCode(:)';
    d = max(abs(fromFile - fromCode));
    assert(d < 1e-12, 't14:crossCheck', ...
        ['trait set 1 ("%s") triple for "%s" does not match the project default [%s] ' ...
         '(max abs diff %.3e); the SI parse or the mapping is wrong, so the run must fail ' ...
         'instead of producing numbers'], setNames{1}, f, mat2str(fromCode), d);
end
fprintf(logFid, 'cross-check    : trait set 1 ("%s") == project defaults (bitwise) for %s\n', ...
    setNames{1}, strjoin(traitFields, ', '));

if isfield(ctx, 'dryRun') && ctx.dryRun
    out = struct();
    out.dry_run = true;
    out.parent_run_id = fitRun.runId;
    out.validated = struct('n_trait_sets', numel(setNames), 'sets', {setNames'}, ...
        'fixed_lag', PD.fixed_lag, 'r0T_points', numel(PD.r0T_grid.from:PD.r0T_grid.step:PD.r0T_grid.to));
    fprintf(logFid, '[dry-run] parent, config and cross-check readable; no computation performed\n');
    return
end

Tgrid = (double(PD.r0T_grid.from) : double(PD.r0T_grid.step) : double(PD.r0T_grid.to))';
nW = ceil(numel(D.T_daily) / 7);
wk = floor(((1:numel(D.T_daily))' - 1) / 7) + 1;
wk = min(max(wk, 1), nW);
sp = struct('onset_doy', th0.onset_doy, 'offset_doy', th0.offset_doy, ...
    'kappa_days', pRef.kappa_days, 'gamma', th0.gamma, 'beta', th0.beta, ...
    's_floor', pRef.s_floor, 'min_window_days', pRef.min_window_days, ...
    'normalize', true, 'norm_factor', NaN);
s = mwspm_season(D.doy, sp);

o = cfg.optimiser;
fitOpt = struct('n_screen', double(o.n_screen), 'seeds', double(o.seeds(:)'), ...
    'n_seed', double(o.n_seed), 'restarts', double(o.restarts), ...
    'maxFunEvals', double(o.maxFunEvals), 'maxIter', double(o.maxIter), ...
    'TolX', double(o.TolX), 'TolFun', double(o.TolFun), 'verbose', false);
freeAll = cellstr(string(PD.free_full));
C = struct('free', {freeAll}, 'fixed', struct(), 'box', mwspm_calib_theta('box'), ...
    'family', 'poisson', 'family_opts', struct(), 'sub_weeks', [], 'y_obs', [], ...
    'phase_penalty', true, 'target_peak', 37.0, 'pen', 300.0, 'd_list', double(PD.fixed_lag));

rowsFixed = {};
rowsRefit = {};
anchorRows = {};
scrRows = {};
base = struct();

for si = 1:numel(setNames)
    nm  = setNames{si};
    set = sets.(nm);
    coef = struct();
    for i = 1:numel(traitFields)
        if ~isfield(set, traitFields{i}), continue; end
        v = double(set.(traitFields{i}));      % temp: `x.(f) = y(:)'` trips
        coef.(traitFields{i}) = reshape(v, 1, []);   % MIR:paren_index_must_be_last
    end
    alt = struct();
    if isfield(set, 'alt') && ~isempty(set.alt)
        alt = set.alt;                          % C1: named alternative forms
    end
    trOpts = struct('coef', coef);
    if ~isempty(fieldnames(alt))
        trOpts.alt = alt;
    end
    preOpts = struct('mu_variant', 'day', 'coef', coef);
    if ~isempty(fieldnames(alt))
        preOpts.alt = alt;
    end
    % coverage: every trait must be declared, either as a triple or as a form
    for i = 1:numel(traitFields)
        f = traitFields{i};
        covered = isfield(coef, f) || isfield(alt, f);
        if strcmp(f, 'lf')
            covered = covered || isfield(alt, 'mu');   % direct mortality replaces lifespan
        end
        assert(covered, 't14:badConfig', ...
            'trait set "%s" declares neither a [c T0 Tm] triple nor an alt form for "%s"', nm, f);
    end
    assert(isfield(coef, 'lf') || isfield(alt, 'mu'), 't14:badConfig', ...
        'trait set "%s" must give either lf (lifespan) or alt.mu (direct mortality)', nm);
    label = char(string(getdef(set, 'label', nm)));

    tr    = mwspm_thermal_traits(D.T_daily, 'raw', trOpts);
    muEff = th0.mu_s * tr.mu_raw;
    R0d   = mwspm_r0(D.T_daily, th0.m0 * s, muEff, pRef.gamma_h, trOpts);
    R0w   = accumarray(wk, R0d, [nW 1], @mean, NaN);
    gt1   = find(R0d > 1);
    if isempty(gt1), w1 = NaN; w2 = NaN; else, w1 = wk(gt1(1)); w2 = wk(gt1(end)); end
    [mx, iMx] = max(R0d);

    % ---- (i) fixed parameters ------------------------------------------------
    pre = mwspm_calib_prepare(D, modelCard, preOpts);
    [tauF, cumF, T0F] = t1404_tau_star(pre, th0, R0w);
    kernF = mwspm_calib_kernel(pre, th0, struct());

    % ---- R0(T) anchors + Andraud screen -------------------------------------
    trT = mwspm_thermal_traits(Tgrid, 'raw', trOpts);
    R0T = mwspm_r0(Tgrid, th0.m0, th0.mu_s * trT.mu_raw, pRef.gamma_h, trOpts);
    pos = find(R0T > 0);
    if isempty(pos), pkT = NaN; loT = NaN; hiT = NaN; pkV = NaN;
    else, loT = Tgrid(pos(1)); hiT = Tgrid(pos(end)); [pkV, ip] = max(R0T); pkT = Tgrid(ip); end
    if isfield(set, 'anchor_key') && ~isempty(set.anchor_key)
        anchorKey = char(string(set.anchor_key));
    else
        anchorKey = t1404_anchor_key(nm);
    end
    A = cfg.mordecai2017_published_anchors.(anchorKey);
    anchorRows{end+1, 1} = {nm, label, pkT, pkV, loT, hiT, ...
        pkT - A.peak_temperature_C.mean, loT - A.lower_zero_temperature_C.mean, hiT - A.upper_zero_temperature_C.mean}; %#ok<AGROW>

    Ts = mean(D.T_daily(s >= 0.5));
    trS = mwspm_thermal_traits(Ts, 'raw', trOpts);
    Rg = cfg.andraud2012_ranges;
    vals = [1 / max(th0.mu_s * trS.mu_raw, realmin), 1 / max(trS.sigma_v, realmin), ...
            trS.a, trS.b, trS.c];
    los  = [Rg.adult_lifespan_days(1), Rg.EIP_days(1), Rg.biting_rate_per_day(1), ...
            Rg.beta_mosquito_to_human(1), Rg.beta_human_to_mosquito(1)];
    his  = [Rg.adult_lifespan_days(2), Rg.EIP_days(2), Rg.biting_rate_per_day(2), ...
            Rg.beta_mosquito_to_human(2), Rg.beta_human_to_mosquito(2)];
    for k = 1:numel(vals)
        scrRows{end+1, 1} = {nm, sprintf('q%d', k), vals(k), los(k), his(k), double(vals(k) >= los(k) && vals(k) <= his(k))}; %#ok<AGROW>
    end
    nInside = sum(vals >= los & vals <= his);

    fprintf(logFid, '\n--- trait set %d/%d: %s ---\n', si, numel(setNames), label);
    fprintf(logFid, '  FIXED params : R0max=%.4f (week %d) | R0>1 weeks [%s,%s] | max|1-mu*dt|=%.4f | tau*=%.0f | reported total=%.3f\n', ...
        mx, wk(iMx), num2strOrNA(w1), num2strOrNA(w2), max(abs(1 - th0.mu_s * tr.mu_raw)), tauF, kernF.reported_total);
    fprintf(logFid, '  R0(T) anchors: peak %.2f C | lower zero %.2f C | upper zero %.2f C  (published: %.1f | %.1f | %.1f)\n', ...
        pkT, loT, hiT, A.peak_temperature_C.mean, A.lower_zero_temperature_C.mean, A.upper_zero_temperature_C.mean);
    fprintf(logFid, '  Andraud 2012 : %d/5 inside at T=%.3f C\n', nInside, Ts);
    tr28 = mwspm_thermal_traits(28, 'raw', trOpts);
    fprintf(logFid, ['  forms        : mu_source=%s | form_override=%d | @28.0 C: a=%.6f ' ...
        'sigma_v=%.6f (EIP %.4f d) mu_raw=%.6f (lf-equiv %.3f d) b=%.6f c=%.6f\n'], ...
        tr.mu_source, tr.form_override, tr28.a, tr28.sigma_v, 1 / tr28.sigma_v, ...
        tr28.mu_raw, 1 / tr28.mu_raw, tr28.b, tr28.c);

    % ---- (ii) refit under the swapped curves --------------------------------
    te = tic;
    [b, info] = mwspm_calib_fit(pre, C, mergeOpt(fitOpt, struct('warm', {{th0}})));
    el = toc(te);
    R0wB = accumarray(wk, mwspm_r0(D.T_daily, b.theta.m0 * s, b.theta.mu_s * tr.mu_raw, pRef.gamma_h, trOpts), [nW 1], @mean, NaN);
    preB = mwspm_calib_prepare(D, modelCard, preOpts);
    [tauB, cumB, T0B] = t1404_tau_star(preB, b.theta, R0wB);
    nw = pRef.n_weeks;
    kParam = 8;
    aic = 2 * b.nll_full + 2 * kParam;
    bic = 2 * b.nll_full + kParam * log(nw);
    fprintf(logFid, '  REFIT        : nll=%.4f nll_full=%.4f AIC=%.3f BIC=%.3f | m0=%.6f rho=%.6f mu_s=%.4f onset=%.3f offset=%.3f gamma=%.4f beta=%.4f | peak week=%d | reported total=%.3f | tau*=%.0f | %.1f s\n', ...
        b.nll, b.nll_full, aic, bic, b.theta.m0, b.theta.rho, b.theta.mu_s, b.theta.onset_doy, ...
        b.theta.offset_doy, b.theta.gamma, b.theta.beta, b.peak_week, b.reported_total, tauB, el);

    rowsFixed{end+1, 1} = {nm, label, mx, wk(iMx), w1, w2, numel(gt1), tauF, T0F, cumF, ...
        kernF.reported_total, pkT, loT, hiT}; %#ok<AGROW>
    rowsRefit{end+1, 1} = {nm, label, double(b.nll), double(b.nll_full), aic, bic, ...
        double(b.theta.m0), double(b.theta.rho), double(b.theta.mu_s), double(b.theta.onset_doy), ...
        double(b.theta.offset_doy), double(b.theta.gamma), double(b.theta.beta), ...
        double(b.peak_week), double(b.reported_total), tauB, el, double(info.n_stage2_within_1nll)}; %#ok<AGROW>
    if si == 1
        base = struct('tau_fixed', tauF, 'tau_refit', tauB, 'R0max', mx, 'nll', b.nll, 'BIC', bic);
    end
end

% ---------------------------------------------------------------- outputs ----
writeNum(fullfile(ctx.dataDir, 'trait_set_fixed_params.csv'), ...
    {'trait_set', 'label', 'R0_max', 'R0_peak_week', 'R0_gt1_first_week', 'R0_gt1_last_week', ...
     'R0_days_gt1', 'tau_star', 'T0_at_tau', 'cum_detected_at_tau', 'reported_total', ...
     'R0T_peak_T', 'R0T_lower_zero_T', 'R0T_upper_zero_T'}, rowsFixed);
writeNum(fullfile(ctx.dataDir, 'trait_set_refits.csv'), ...
    {'trait_set', 'label', 'nll', 'nll_full', 'AIC', 'BIC', 'm0', 'rho', 'mu_s', 'onset_doy', ...
     'offset_doy', 'gamma', 'beta', 'peak_week', 'reported_total', 'tau_star_refit', ...
     'elapsed_seconds', 'n_stage2_within_1nll'}, rowsRefit);
writeNum(fullfile(ctx.dataDir, 'r0T_anchors_by_set.csv'), ...
    {'trait_set', 'label', 'R0T_peak_T', 'R0T_peak_value', 'R0T_lower_zero_T', 'R0T_upper_zero_T', ...
     'd_peak_vs_published', 'd_lower_vs_published', 'd_upper_vs_published'}, anchorRows);
writeNum(fullfile(ctx.dataDir, 'andraud_range_screen_by_set.csv'), ...
    {'trait_set', 'quantity', 'value', 'range_low', 'range_high', 'inside_range'}, scrRows);

summary = struct();
summary.experiment_id = ctx.expId;
summary.task = char(string(getdef(cfg, 'task', '-')));
summary.parent_fit_run_id = fitRun.runId;
summary.main_report_point = fit.main_report_point;
summary.deterministic = true;
summary.rng = 'none (deterministic recomputation and deterministic optimisation)';
summary.trait_sets = sets;
summary.source_swap_scope = cfg.source_swap_scope;
summary.cross_check = 'Table A triples matched the project defaults bitwise (asserted in code)';
summary.mean_abs_mu_eff = 'see run_log per set';
summary.second_publication_status = struct('verdict', 'NOT PERFORMED', 'reason', ...
    ['no second PUBLICATION with trait-level thermal coefficients has been obtained; ' ...
     'this experiment swaps species within the SAME source (Mordecai 2017 Table A vs Table B). ' ...
     'It must never be reported as "another source".']);
summary.honesty = cfg.honesty;
summary.value_type_note = ['deterministic_point for every row: the fitted parameters were either held ' ...
    'fixed or re-optimised under a different trait set. No interval, no replicate distribution.'];
mwspm_writejson(fullfile(ctx.dataDir, 'summary.json'), summary);

out = struct();
out.status = 'SUCCESS';
out.parent_run_id = fitRun.runId;
out.n_trait_sets = numel(setNames);
out.deterministic = true;
out.second_publication_status = 'NOT PERFORMED';
out.outputs = {'trait_set_fixed_params.csv', 'trait_set_refits.csv', 'r0T_anchors_by_set.csv', ...
    'andraud_range_screen_by_set.csv', 'summary.json'};
out.note = ['T14.04. SAME-SOURCE species swap (Mordecai 2017 Table A albopictus vs Table B aegypti), ' ...
    'reported for fixed parameters AND after a refit. tau* is recomputed with the gate-G1 caliber.'];
end

% ---------------------------------------------------------------------------
function [tauStar, cumAtTau, T0atTau] = t1404_tau_star(pre, th, R0w)
T0w = nan(numel(R0w), 1);
ok = R0w > 1;
T0w(ok) = th.rho ./ log(R0w(ok));
kern = mwspm_calib_kernel(pre, th, struct());
cumDet = cumsum(kern.Yhat(:));
cross = find(isfinite(T0w) & cumDet >= T0w, 1);
if isempty(cross)
    tauStar = NaN; cumAtTau = NaN; T0atTau = NaN;
else
    tauStar = cross; cumAtTau = cumDet(cross); T0atTau = T0w(cross);
end
end

% ---------------------------------------------------------------------------
function writeNum(path, hdr, rows)
fid = fopen(path, 'w', 'n', 'UTF-8');
assert(fid > 0, 't14:io', 'cannot write %s', path);
fprintf(fid, '%s\n', strjoin(hdr, ','));
for i = 1:numel(rows)
    rw = rows{i};
    cellTxt = false(1, numel(rw));
    for k = 1:numel(rw)
        cellTxt(k) = ischar(rw{k}) || (isstring(rw{k}) && ~isscalar(rw{k}));
    end
    parts = cell(1, numel(rw));
    for k = 1:numel(rw)
        if cellTxt(k)
            parts{k} = char(string(rw{k}));
        elseif isfinite(double(rw{k}))
            parts{k} = sprintf('%.17g', double(rw{k}));
        else
            parts{k} = 'NA';
        end
    end
    fprintf(fid, '%s\n', strjoin(parts, ','));
end
fclose(fid);
end

% ---------------------------------------------------------------------------
function key = t1404_anchor_key(setName)
if contains(lower(setName), 'aegypti')
    key = 'aegypti';
else
    key = 'albopictus';
end
end

% ---------------------------------------------------------------------------
function th = t1404_theta_from_fit(tf)
th = struct('m0', double(tf.m0), 'rho', double(tf.rho), 'onset_doy', double(tf.onset_doy), ...
    'offset_doy', double(tf.offset_doy), 'mu_s', double(tf.mu_s), 'gamma', double(tf.gamma), ...
    'beta', double(tf.beta), 'd_weeks', double(tf.d_weeks));
end

% ---------------------------------------------------------------------------
function P = t1404_find_parent(root, wantId, expId, relPath)
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
function o = mergeOpt(base, over)
o = base;
fn = fieldnames(over);
for i = 1:numel(fn)
    o.(fn{i}) = over.(fn{i});
end
end

function s = num2strOrNA(v)
if isfinite(double(v))
    s = sprintf('%.6g', double(v));
else
    s = 'NA';
end
end

function v = getdef(s, f, dflt)
if isfield(s, f) && ~isempty(s.(f))
    v = s.(f);
else
    v = dflt;
end
end
