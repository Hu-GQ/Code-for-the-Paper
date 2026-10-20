function fig2_fit_main()
%FIG2_FIT_MAIN  Fig 2 - fit vs observation, raw and log1p scales.
%
%   Figure-level statement: the calibrated M-WSPM reproduces the observed season
%   on both scales (R2_raw 0.950, R2_log 0.976); reported total 767.2 vs 767 observed.
%
%   Read-only data: the analysis  20261006T012933_c2b62344/data/weekly_fit.csv
%   Claims: the statement, the statement.

pp_util('style');
d = pp_util('dirs');
pal = pp_util('palette');

runA = pp_util('runDir', '20261006T012933_c2b62344');
fA = fullfile(runA, 'data', 'weekly_fit.csv');
snap = pp_util('snap', {fA});

S = mwspm_readcsv(fA);
wk   = pp_util('col', S, 'week_index');
obs  = pp_util('col', S, 'observed_local');
rep  = pp_util('col', S, 'model_reported');
true_ = pp_util('col', S, 'model_true_incidence');

fig = figure('Units', 'inches', 'Position', [0.5 0.5 4.90 6.10]);
tl = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% ===== panel a: raw scale =====
ax1 = nexttile(tl); hold(ax1, 'on');
xline(ax1, 19, '--', 'Color', [0.72 0.72 0.72], 'LineWidth', 0.9, 'HandleVisibility', 'off');
pp = plot(ax1, wk, obs, 'o', 'Color', [0.15 0.15 0.15], 'MarkerFaceColor', [0.15 0.15 0.15], ...
    'MarkerSize', 3.6, 'DisplayName', 'observed (local)');
plot(ax1, wk, rep, '-', 'Color', pal.blue, 'LineWidth', 1.3, 'DisplayName', 'model: reported');
plot(ax1, wk, true_, '-', 'Color', pal.orange, 'LineWidth', 1.3, 'DisplayName', 'model: true incidence');
xlim(ax1, [0 48]); ylim(ax1, [0 85]);
xlabel(ax1, 'Analysis week (1-47)');
ylabel(ax1, 'Weekly cases');
lg = legend(ax1, 'Location', 'northwest', 'FontSize', 7); lg.ItemTokenSize = 12;
text(ax1, 2.2, 61, 'R^2_{raw} = 0.950', 'FontSize', 7.5, 'Color', pal.blue, 'FontWeight', 'bold');
text(ax1, 2.2, 54, 'reported total 767.2 / observed 767', 'FontSize', 7.2, 'Color', [0.35 0.35 0.35]);
text(ax1, 18.3, 5, 'first detected cases', 'FontSize', 6.6, 'Color', [0.45 0.45 0.45], ...
    'Rotation', 90);
title(ax1, '(a) Fit vs observation (raw scale)', 'FontSize', 8.5, 'FontWeight', 'normal');

% ===== panel b: log1p scale =====
ax2 = nexttile(tl); hold(ax2, 'on');
xline(ax2, 19, '--', 'Color', [0.72 0.72 0.72], 'LineWidth', 0.9, 'HandleVisibility', 'off');
plot(ax2, wk, log1p(obs), 'o', 'Color', [0.15 0.15 0.15], 'MarkerFaceColor', [0.15 0.15 0.15], 'MarkerSize', 3.6);
plot(ax2, wk, log1p(rep), '-', 'Color', pal.blue, 'LineWidth', 1.3);
plot(ax2, wk, log1p(true_), '-', 'Color', pal.orange, 'LineWidth', 1.3);
xlim(ax2, [0 48]); ylim(ax2, [0 5.0]);
xlabel(ax2, 'Analysis week (1-47)');
ylabel(ax2, 'log1p(weekly cases)');
text(ax2, 2.2, 4.62, 'R^2_{log} = 0.976', 'FontSize', 7.5, 'Color', pal.blue, 'FontWeight', 'bold');
text(ax2, 2.2, 4.18, 'shape (log) and level (raw) are reproduced simultaneously', ...
    'FontSize', 7.0, 'Color', [0.35 0.35 0.35]);
title(ax2, '(b) Same data on log1p scale', 'FontSize', 8.5, 'FontWeight', 'normal');

% ---- source data ----
sd = fullfile(d.src, 'fig2_fit_main');
if ~exist(sd, 'dir'); mkdir(sd); end
copyfile(fA, fullfile(sd, 'weekly_fit.csv'));
fprintf('fig2: source data copied.\n');

pp_util('assert', snap);
pp_util('savefig', fig, 'fig2_fit_main');
fprintf('fig2: done.\n');
end
