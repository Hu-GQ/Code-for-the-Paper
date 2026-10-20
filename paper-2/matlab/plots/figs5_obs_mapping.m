function figs5_obs_mapping()
%FIGS5_OBS_MAPPING  Fig S5 - observation mapping (reported vs true split).
%
%   Read-only: data/processed/weekly_analysis_input.csv (frozen input layer,
%   explicitly annotated historical save) + the analysis 20261006T012933_c2b62344/
%   data/weekly_fit.csv.
%   Claims: the statement / the statement (observation caliber; reported-vs-true separation).
%
%   Panel (a): weekly local / imported reported cases (stacked), week_index.
%   Panel (b): observed vs model-reported vs model true incidence (log).

pp_util('style');
d = pp_util('dirs');
pal = pp_util('palette');

fIn = fullfile(pp_util('root'), 'data', 'processed', 'weekly_analysis_input.csv');
runA = pp_util('runDir', '20261006T012933_c2b62344');
fFit = fullfile(runA, 'data', 'weekly_fit.csv');
snap = pp_util('snap', {fIn, fFit});

% ===== panel a: reported composition =====
I = mwspm_readcsv(fIn);
wk = pp_util('col', I, 'week_index');
lo = pp_util('col', I, 'local_reported');
im = pp_util('col', I, 'imported_reported');

fig = figure('Units', 'inches', 'Position', [0.3 0.3 7.20 2.85]);
ax1 = subplot(1, 2, 1); hold(ax1, 'on');
b = bar(ax1, wk, [lo, im], 0.85, 'stacked', 'EdgeColor', 'none');
b(1).FaceColor = pal.blue;
b(2).FaceColor = pal.orange;
xline(ax1, 19, '--', 'Color', [0.35 0.35 0.35], 'LineWidth', 0.9);
text(ax1, 17.5, 15, 'first detected cases (w19)', 'FontSize', 6.6, ...
    'Color', [0.35 0.35 0.35], 'Rotation', 90);
xlim(ax1, [0 48]); ylim(ax1, [0 85]);
xlabel(ax1, 'week\_index (1-47)');
ylabel(ax1, 'Reported cases');
text(ax1, 1.5, 79, 'local 767 total', 'FontSize', 6.8, 'Color', pal.blue, 'FontWeight', 'bold');
text(ax1, 1.5, 73.5, 'imported 78 total', 'FontSize', 6.8, 'Color', [0.80 0.45 0.10], 'FontWeight', 'bold');
text(ax1, 1.5, 68, 'GD CDC iso\_week = index + 5', 'FontSize', 6.4, 'Color', [0.40 0.40 0.40]);
title(ax1, '(a) Weekly reported series (frozen input)', 'FontSize', 8.2, 'FontWeight', 'bold');

% ===== panel b: observed / reported / true =====
F = mwspm_readcsv(fFit);
fw = pp_util('col', F, 'week_index');
fol = pp_util('col', F, 'observed_local');
fmr = pp_util('col', F, 'model_reported');
fmt = pp_util('col', F, 'model_true_incidence');

ax2 = subplot(1, 2, 2); hold(ax2, 'on');
plot(ax2, fw, max(fol, 0.7), 'o', 'Color', [0.2 0.2 0.2], 'MarkerSize', 2.6, ...
    'MarkerFaceColor', [0.2 0.2 0.2]);
plot(ax2, fw, max(fmt, 0.7), '-', 'Color', pal.orange, 'LineWidth', 1.2);
plot(ax2, fw, max(fmr, 0.7), '--', 'Color', pal.blue, 'LineWidth', 1.1);
set(ax2, 'YScale', 'log');
xlim(ax2, [0 48]); ylim(ax2, [0.5 1e4]);
xlabel(ax2, 'week\_index (1-47)');
ylabel(ax2, 'Cases per week (log_{10})');
text(ax2, 20.5, 4200, 'true incidence (model)', 'FontSize', 6.8, 'Color', pal.orange);
text(ax2, 20.5, 4.3, 'model reported = E[obs]', 'FontSize', 6.8, 'Color', pal.blue);
text(ax2, 20.5, 0.85, 'observed local', 'FontSize', 6.8, 'Color', [0.2 0.2 0.2]);
text(ax2, 1.5, 2500, 'reported = \rho x true (shift d = 1 wk)', 'FontSize', 6.6, 'Color', [0.35 0.35 0.35]);
title(ax2, '(b) Reported vs true split (primary report point)', 'FontSize', 8.2, 'FontWeight', 'bold');

% ===== source data =====
sd = fullfile(d.src, 'figs5_obs_mapping');
if ~exist(sd, 'dir'); mkdir(sd); end
fid = fopen(fullfile(sd, 'reported_composition.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, 'week_index,local_reported,imported_reported\n');
for k = 1:numel(wk); fprintf(fid, '%g,%g,%g\n', wk(k), lo(k), im(k)); end
fclose(fid);
fid = fopen(fullfile(sd, 'reported_vs_true.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, 'week_index,observed_local,model_reported,model_true_incidence\n');
for k = 1:numel(fw); fprintf(fid, '%g,%g,%g,%.6g\n', fw(k), fol(k), fmr(k), fmt(k)); end
fclose(fid);

fprintf('s5 totals: local %g, imported %g\n', sum(lo), sum(im));

pp_util('assert', snap);
pp_util('savefig_supp', fig, 'figs5_obs_mapping');
fprintf('s5: done.\n');
end
