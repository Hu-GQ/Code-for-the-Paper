function fig5_risk_window()
%FIG5_RISK_WINDOW  Fig 5 - the seasonal R0(t) window against the 2025 outbreak.
%
%   Deterministic export (explicitly allowed by the workflow, step 17):
%   R0(t) is computed from the SAVED main-report-point theta (the analysis
%   fit_best.json) and the SAVED temperature input via the project's
%   validated mwspm chain (mwspm_season / mwspm_thermal_traits / mwspm_r0,
%   N_v supply caliber, continuous-mortality 'raw' variant).  No ODE solve,
%   no re-fit, no parameter change.
%
%   Figure-level statement: the R0(t)>1 window (first week 15-17) contains the
%   outbreak take-off; the warning week (18-20) falls inside it, coinciding
%   with the first detected week - no lead time.
%
%   Read-only inputs: the analysis fit_best.json; the analysis gate_g1_weeks.csv;
%   model card config/experiments/model_mwspm_gen3_rhohat.json.
%   Claims: the statement, the statement, the statement.

pp_util('style');
d = pp_util('dirs');
pal = pp_util('palette');
root = pp_util('root');

runFit  = pp_util('runDir', '20261006T012933_c2b62344');
runGate = pp_util('runDir', '20261006T020359_3295cf28');
fFit  = fullfile(runFit,  'data', 'fit_best.json');
fGate = fullfile(runGate, 'data', 'gate_g1_weeks.csv');
fMC   = fullfile(root, 'config', 'experiments', 'model_mwspm_gen3_rhohat.json');
snap = pp_util('snap', {fFit, fGate, fMC});

mc  = mwspm_readjson(fMC);
fit = mwspm_readjson(fFit);
th  = fit.main_report_point.theta;

D    = mwspm_load_data(root, mc);
pre  = mwspm_calib_prepare(D, mc, struct('mu_variant', 'day'));
pRef = mwspm_params(mc);

% ---- deterministic export of R0(t) (same chain as the analysis's unc_tau_star) --
sp = struct('onset_doy', th.onset_doy, 'offset_doy', th.offset_doy, ...
    'kappa_days', pre.kappa_days, 'gamma', th.gamma, 'beta', th.beta, ...
    's_floor', pre.s_floor, 'min_window_days', pre.min_window_days, ...
    'normalize', true, 'norm_factor', NaN);
s = mwspm_season(pre.doy, sp);
trRaw = mwspm_thermal_traits(D.T_daily, 'raw');
muEff = th.mu_s * trRaw.mu;
R0daily = mwspm_r0(D.T_daily, th.m0 * s, muEff, pRef.gamma_h);

% sanity print vs registered values (reported, not hard-asserted).
% NOTE (T17, 2026-10-08): the "107-232 / 96 of 126" figure quoted in
% exp_c4_gate_g1.m and main_results.md §4 is the LEGACY source-faithful
% caliber (reproduced by the E0 r0_daily.csv: first 107 / last 232 / 96 days).
% The MAIN REPORT POINT theta gives [114, 239] / 96 of 126 - the window
% shifts by one week; both calibers must be quoted explicitly.
gt = find(R0daily > 1);
fprintf('fig5: main report point: R0>1 days: %d of %d; first day %d (week %.2f), last day %d (week %.2f)\n', ...
    numel(gt), numel(R0daily), gt(1), (gt(1) - 1) / 7 + 1, gt(end), (gt(end) - 1) / 7 + 1);
inWin = sum(R0daily(gt(1):gt(end)) > 1);
fprintf('fig5: span [%d, %d] = %d days; %d of them have R0>1 (legacy caliber: [107,232], 96 of 126)\n', ...
    gt(1), gt(end), gt(end) - gt(1) + 1, inWin);

day = (1:numel(R0daily))';
nW  = pre.n_weeks;
wk  = min(max(floor((day - 1) / 7) + 1, 1), nW);
R0w = accumarray(wk, R0daily, [nW 1], @mean, NaN);
obs = D.Y_weekly(:);

fig = figure('Units', 'inches', 'Position', [0.5 0.5 7.10 6.00]);
tl = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, 'Seasonal R_0(t) window and the 2025 local outbreak (primary fit)', ...
    'FontWeight', 'bold', 'FontSize', 9.5);

% ===== panel a: temperature =====
ax1 = nexttile(tl); hold(ax1, 'on');
yline(ax1, 31.51, '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 0.9, 'HandleVisibility', 'off');
yline(ax1, 13.41, '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 0.9, 'HandleVisibility', 'off');
plot(ax1, day, D.T_daily, '-', 'Color', [0.70 0.25 0.25], 'LineWidth', 0.9);
xlim(ax1, [0 330]); ylim(ax1, [11 34]);
xlabel(ax1, 'Day index (1 = 2025-02-03)');
ylabel(ax1, 'Temperature (deg C)');
text(ax1, 5, 32.1, 'lifespan limits of the baseline trait set: 31.51 / 13.41 degC', ...
    'FontSize', 6.8, 'Color', [0.45 0.45 0.45]);
title(ax1, '(a) Daily mean temperature', 'FontSize', 8.5, 'FontWeight', 'normal');

% ===== panel b: R0(t) + outbreak (left axis = R0, right axis = cases) =====
ax2 = nexttile(tl); hold(ax2, 'on');
% shade R0 > 1 segments (daily)
on = R0daily > 1;
edges = diff([false; on; false]);
segStart = day(edges(1:end-1) == 1);
segEnd   = day(edges(2:end) == -1);
for k = 1:numel(segStart)
    patch(ax2, [segStart(k) segEnd(k) segEnd(k) segStart(k)], [0 0 4.2 4.2], ...
        [0.90 0.93 0.97], 'EdgeColor', 'none', 'HandleVisibility', 'off');
end
plot(ax2, day, R0daily, '-', 'Color', pal.blue_l, 'LineWidth', 0.7, ...
    'DisplayName', 'R_0(t), daily');
plot(ax2, (1:nW)' * 7 - 3.5, R0w, '-o', 'Color', pal.blue, 'MarkerSize', 3, ...
    'LineWidth', 1.3, 'DisplayName', 'R_0 weekly mean');
yline(ax2, 1, '-', 'Color', [0.3 0.3 0.3], 'LineWidth', 0.8, 'HandleVisibility', 'off');
ylim(ax2, [0 4.2]);
ylabel(ax2, 'R_0(t)');
text(ax2, 20, 3.78, 'onset index on the pre-specified grid: weeks 19-20 on both lags', ...
    'FontSize', 6.9, 'Color', [0.35 0.35 0.35]);
text(ax2, 20, 3.42, sprintf('R_0 > 1 on %d of the %d days from day %d to %d', ...
    numel(gt), gt(end) - gt(1) + 1, gt(1), gt(end)), 'FontSize', 6.9, 'Color', [0.30 0.30 0.55]);
xline(ax2, 127, ':', 'Color', [0.35 0.35 0.35], 'LineWidth', 1.0, 'HandleVisibility', 'off');

yyaxis(ax2, 'right');
bar(ax2, (1:nW)' * 7 - 3.5, obs, 6, 'FaceColor', pal.orange, 'EdgeColor', 'none', ...
    'LineWidth', 0.3, 'FaceAlpha', 0.22, 'DisplayName', 'observed local (weekly)');
ylim(ax2, [0 240]);
ylabel(ax2, 'Weekly local cases');
ax2.YColor = [0.55 0.40 0.15];

xlim(ax2, [0 330]);
xlabel(ax2, 'Day index (1 = 2025-02-03)');
lg = legend(ax2, 'Location', 'northeast', 'FontSize', 7); lg.ItemTokenSize = 12;
title(ax2, '(b) R_0(t) with the R_0>1 window and the outbreak', 'FontSize', 8.5, 'FontWeight', 'normal');

% ---- source data ----
sd = fullfile(d.src, 'fig5_risk_window');
if ~exist(sd, 'dir'); mkdir(sd); end
wsave(fullfile(sd, 'r0_daily_main_report_point.csv'), ...
    {'day', 'T_c', 'mu_eff', 'season_s', 'R0'}, {day, D.T_daily, muEff, s, R0daily});
wsave(fullfile(sd, 'weekly_r0_and_observed.csv'), ...
    {'week_index', 'R0_week_mean', 'observed_local'}, {(1:nW)', R0w, obs});

pp_util('assert', snap);
pp_util('savefig', fig, 'fig5_risk_window');
fprintf('fig5: done.\n');
end

% =========================================================================
function wsave(path, names, cols)
fid = fopen(path, 'w', 'n', 'UTF-8');
assert(fid > 0, 'fig5:write', 'cannot write %s', path);
fprintf(fid, '%s\n', strjoin(names, ','));
for k = 1:numel(cols{1})
    for j = 1:numel(cols)
        v = cols{j}(k);
        if isnan(v)
            fprintf(fid, 'NA');
        elseif abs(v - round(v)) < 1e-9
            fprintf(fid, '%d', round(v));
        else
            fprintf(fid, '%.10g', v);
        end
        if j < numel(cols); fprintf(fid, ','); end
    end
    fprintf(fid, '\n');
end
fclose(fid);
end
