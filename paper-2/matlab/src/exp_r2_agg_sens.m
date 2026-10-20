function payload = exp_r2_agg_sens(ctx)
%EXP_R2_AGG_SENS  Weekly R0 aggregation sensitivity of tau*, plus the
%   rho-cancelled (fixed-trajectory) onset week.
%
%   The onset index is defined on a weekly object R0(w), but the model only
%   produces a daily R0(t) that is sawtoothed inside the favourable window.
%   The profile run aggregates it to the week by the arithmetic mean of the
%   seven daily values.  That convention was never stated in the manuscript,
%   so the reviewers asked what depends on it.  This experiment recomputes
%   tau* at two conditional points under four stated aggregations:
%
%     mean      arithmetic mean of the daily R0 inside the week  (profile rule)
%     median    median of the daily R0 inside the week
%     midpoint  R0 evaluated on the middle day of the week only
%     max       maximum of the daily R0 inside the week
%
%   Nothing is refitted: each point keeps the (theta, rho) of the profile
%   row, so the expected-report curve is identical across aggregations and
%   only T0(w) moves.  This isolates the definitional choice from the fit.
%
%   The same solve also reports the rho-cancelled value of Remark 3.4
%   (fixed latent trajectory, rho divided out), which the manuscript cites
%   only as a dotted line in Fig. 1(a).
%
%   Output is a definitional sensitivity table, not a confidence set.

cfg = ctx.cfg;
assert(~cfg.phase_penalty, 'R2agg:penalty', 'this experiment must not use a peak penalty');

if isfield(ctx, 'dryRun') && ctx.dryRun
    payload = struct('outputs', {{}}); payload.dry_run = true; return
end

diary(fullfile(ctx.runDir, 'logs', 'stdout.txt'));
diaryGuard = onCleanup(@() diary('off')); %#ok<NASGU>
snap = fullfile(ctx.dataDir, 'code_snapshot');
mkdir(snap);
copyfile(fullfile(ctx.root, 'src', 'matlab', '*.m'), snap);

D    = mwspm_load_data(ctx.root, ctx.modelCard);
pre  = mwspm_calib_prepare(D, ctx.modelCard, struct('mu_variant', cfg.mu_variant));
pRef = mwspm_params(ctx.modelCard);
nW   = pre.n_weeks;

profPath = fullfile(ctx.root, 'results', 'runs', cfg.profile_run, 'data', 'rho_profile_pure.csv');
assert(exist(profPath, 'file') == 2, 'R2agg:noProfile', 'missing %s', profPath);
P = mwspm_readcsv(profPath);
lagP = agg_col(P, 'lag'); rhoP = agg_col(P, 'rho');

aggs = {'mean', 'median', 'midpoint', 'max'};
rows = table();

pts = cfg.points;
if iscell(pts), pts = [pts{:}]; end
for k = 1:numel(pts)
    pt = pts(k);
    idx = find(lagP == double(pt.lag) & abs(rhoP - double(pt.rho)) < 1e-9, 1);
    assert(~isempty(idx), 'R2agg:row', 'no profile row for lag %d rho %g', pt.lag, pt.rho);

    th = struct();
    th.rho = rhoP(idx);
    th.d_weeks = lagP(idx);
    for f = {'m0', 'mu_s', 'onset_doy', 'offset_doy', 'gamma', 'beta'}
        col = agg_col(P, f{1});
        th.(f{1}) = col(idx);
    end

    % Daily R0 exactly as the profile computed it (continuous-rate traits,
    % day-variant supply) and the expected reported series at this theta.
    sp = struct('onset_doy', th.onset_doy, 'offset_doy', th.offset_doy, ...
        'kappa_days', pre.kappa_days, 'gamma', th.gamma, 'beta', th.beta, ...
        's_floor', pre.s_floor, 'min_window_days', pre.min_window_days, ...
        'normalize', true, 'norm_factor', NaN);
    s = mwspm_season(pre.doy, sp);
    R0d = mwspm_r0(D.T_daily, th.m0 * s, th.mu_s * pre.mu_raw_daily, pRef.gamma_h, struct());

    wk = floor(((1:numel(R0d))' - 1) / 7) + 1;
    wk = min(max(wk, 1), nW);

    kern   = mwspm_calib_kernel(pre, th, cfg.solver);
    cumRep = cumsum(kern.Yhat(:));          % expected reports (contains rho)
    cumLat = cumRep / th.rho;               % latent cumulative C_{w-d} (rho-free)

    for a = 1:numel(aggs)
        R0w = aggregate(R0d, wk, nW, aggs{a});
        [tauS, T0a, cumA, R0a] = crossing(R0w, cumRep, th.rho, true);
        rows = [rows; table(lagP(idx), rhoP(idx), {aggs{a}}, tauS, T0a, cumA, R0a, ...
            sum(isfinite(R0w) & R0w > 1), ...
            'VariableNames', {'lag', 'rho', 'aggregation', 'tau_star', ...
            'T0_at_tau', 'cum_reports_at_tau', 'R0w_at_tau', 'n_weeks_R0gt1'})]; %#ok<AGROW>
    end

    % Fixed-trajectory (rho-cancelled) value: Eq. (12) of the manuscript,
    % cum_lat(w) >= 1 / log R0(w).  Reported once per point, on the profile
    % aggregation, because the cancellation itself does not depend on rho.
    R0w = aggregate(R0d, wk, nW, 'mean');
    [tauS, T0a, cumA, R0a] = crossing(R0w, cumLat, 1.0, true);
    rows = [rows; table(lagP(idx), rhoP(idx), {'fixed-latent'}, tauS, T0a, cumA, R0a, ...
        sum(isfinite(R0w) & R0w > 1), ...
        'VariableNames', {'lag', 'rho', 'aggregation', 'tau_star', ...
        'T0_at_tau', 'cum_reports_at_tau', 'R0w_at_tau', 'n_weeks_R0gt1'})]; %#ok<AGROW>

    fprintf('R2agg lag%d rho=%.2f: %s\n', th.d_weeks, th.rho, ...
        strjoin(cellstr(string(rows.tau_star(end-numel(aggs):end)))', ' / '));
end

csvPath = fullfile(ctx.dataDir, 'agg_sensitivity.csv');
fid = fopen(csvPath, 'w', 'n', 'UTF-8');
fprintf(fid, 'lag,rho,aggregation,tau_star,T0_at_tau,cum_reports_at_tau,R0w_at_tau,n_weeks_R0gt1\n');
fclose(fid);
writetable(rows, csvPath, 'WriteMode', 'append', 'WriteVariableNames', false);

summary = struct();
summary.question = ['does tau* depend on the unstated weekly aggregation of R0(t)?'];
summary.rule_used_in_profile = 'arithmetic mean of the seven daily R0 values inside the week';
summary.alternatives = aggs(2:end);
summary.fixed_latent = ['first week with the rho-free latent cumulative reaching ' ...
    '1/log R0(w); equivalent to Eq. (12) and independent of rho by construction'];
summary.no_refit = true;
summary.interpretation = ['definitional sensitivity only; the expected-report curve and ' ...
    'the fitted parameters are held at the profile row'];
summary.rows = rows;

mwspm_writejson(fullfile(ctx.dataDir, 'agg_summary.json'), summary);
payload = struct('outputs', {{'agg_sensitivity.csv', 'agg_summary.json'}});
end

% ---------------------------------------------------------------------------
function v = aggregate(R0d, wk, nW, how)
switch how
    case 'mean'
        v = accumarray(wk, R0d, [nW 1], @mean, NaN);
    case 'median'
        v = accumarray(wk, R0d, [nW 1], @median, NaN);
    case 'max'
        v = accumarray(wk, R0d, [nW 1], @max, NaN);
    case 'midpoint'
        v = nan(nW, 1);
        for w = 1:nW
            d = 7 * (w - 1) + 4;              % middle day of the week
            if d <= numel(R0d), v(w) = R0d(d); end
        end
    otherwise
        error('R2agg:agg', 'unknown aggregation "%s"', how);
end
end

% ---------------------------------------------------------------------------
function [tauStar, T0at, cumAt, R0at] = crossing(R0w, cum, rhoUsed, useRho)
T0w = nan(size(R0w));
ok = isfinite(R0w) & R0w > 1;
if useRho
    T0w(ok) = rhoUsed ./ log(R0w(ok));
else
    T0w(ok) = 1 ./ log(R0w(ok));
end
c = find(isfinite(T0w) & cum(:) >= T0w, 1);
if isempty(c)
    tauStar = NaN; T0at = NaN; cumAt = NaN; R0at = NaN;
else
    tauStar = c; T0at = T0w(c); cumAt = cum(c); R0at = R0w(c);
end
end

% ---------------------------------------------------------------------------
function v = agg_col(P, name)
j = find(strcmp(P.names, name), 1);
assert(~isempty(j), 'R2agg:col', 'missing column "%s"', name);
v = cellfun(@(r) str2double(r{j}), P.rows);
end
