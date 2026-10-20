function fig1_c4_warning_week()
%FIG1_C4_WARNING_WEEK  Fig 1 - threshold-free warning week (C4 core figure).
%
%   Figure-level claim: the T0-crossing warning week tau*(rho) drifts by only
%   1-2 weeks across the frozen gate grid rho in [0.30, 1.00] (both data lags),
%   while the cumulative-case caliber drifts 16 weeks across thresholds 1/100/500.
%   Honesty panel: cumulative detections at tau* are O(1-3) cases => no lead time.
%
%   Read-only data (saved runs only; no re-solve, no re-fit):
%     T14.01  20261006T104119_0570febb/data/warning_week_vs_rho.csv
%     T13.01  20261006T020359_3295cf28/data/threshold_drift.csv
%   Claims: CL-01, CL-02, CL-03, CL-22.

pp_util('style');
d = pp_util('dirs');
pal = pp_util('palette');

runA = pp_util('runDir', '20261006T104119_0570febb');
runB = pp_util('runDir', '20261006T020359_3295cf28');
fA = fullfile(runA, 'data', 'warning_week_vs_rho.csv');
fB = fullfile(runB, 'data', 'threshold_drift.csv');
snap = pp_util('snap', {fA, fB});

S = mwspm_readcsv(fA);
ch  = pp_util('strcol', S, 'channel');
lag = pp_util('col', S, 'branch_lag');
rho = pp_util('col', S, 'rho');
tau = pp_util('col', S, 'tau_star_week');
cum = pp_util('col', S, 'cum_detected_at_tau_star');
infz = pp_util('col', S, 'in_frozen_grid');

ref1 = strcmp(ch, 'profile_refit') & lag == 1;   ix1 = find(ref1);
ref2 = strcmp(ch, 'profile_refit') & lag == 2;   ix2 = find(ref2);
fix1 = strcmp(ch, 'main_report_point') & lag == 1; ix3 = find(fix1);

T = mwspm_readcsv(fB);
thr = pp_util('col', T, 'threshold_cumulative_cases');   % 1 / 100 / 500
wks = pp_util('col', T, 'first_week_reaching');          % 19 / 26 / 35

% ---- figure ----
fig = figure('Units', 'inches', 'Position', [0.5 0.5 7.10 7.00]);
tl = tiledlayout(fig, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, ['Threshold-free warning week: \tau^*(\rho) vs cumulative-case caliber'], ...
    'FontWeight', 'bold', 'FontSize', 9.5);

% ===== panel a =====
ax1 = nexttile(tl); hold(ax1, 'on');
xr = patch(ax1, [0.30 1.00 1.00 0.30], [14 14 22 22], [0.90 0.90 0.90], ...
    'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(ax1, rho(ix1), tau(ix1), '-o', 'Color', pal.blue, 'MarkerSize', 2.6, ...
    'MarkerFaceColor', pal.blue, 'LineWidth', 1.0, 'DisplayName', 'lag 1 (refit channel)');
plot(ax1, rho(ix2), tau(ix2), '--s', 'Color', pal.red, 'MarkerSize', 2.6, ...
    'MarkerFaceColor', pal.red, 'LineWidth', 1.0, 'DisplayName', 'lag 2 (refit channel)');
plot(ax1, rho(ix3), tau(ix3), ':', 'Color', pal.gray, 'LineWidth', 1.0, ...
    'DisplayName', 'fixed trajectory (\tau^*=20, lag 1)');
ylim(ax1, [15 21]); yticks(ax1, 15:21); xlim(ax1, [0 1.02]);
xlabel(ax1, 'Detection ratio \rho');
ylabel(ax1, 'Warning week \tau^*');
lg = legend(ax1, 'Location', 'northwest', 'FontSize', 7); lg.ItemTokenSize = 12;
lg.Color = 'w';                       % opaque background: the fixed-trajectory
lg.EdgeColor = [0.85 0.85 0.85];      % line runs horizontally through this area
lg.Box = 'on';
text(ax1, 0.655, 15.25, 'frozen gate grid [0.30, 1.00]', 'FontSize', 7, 'Color', [0.35 0.35 0.35]);
text(ax1, 0.985, 17.15, 'frozen grid range: 1 wk (lag 1) / 2 wk (lag 2)', ...
    'FontSize', 7.2, 'Color', pal.blue, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
text(ax1, 0.985, 16.55, 'full-grid range [16, 20]: width depends on the grid (report both)', ...
    'FontSize', 6.8, 'Color', [0.35 0.35 0.35], 'HorizontalAlignment', 'right');
title(ax1, '(a) T_0-crossing week vs detection ratio (profile-refit channel)', 'FontSize', 8.5, 'FontWeight', 'normal');

% ===== panel b =====
ax2 = nexttile(tl); hold(ax2, 'on');
% T0 caliber intervals (frozen grid): lag1 [19,20], lag2 [18,20]
plot(ax2, [1 1], [19 20], '-', 'Color', pal.blue, 'LineWidth', 3);
plot(ax2, 1, 19.5, 'o', 'Color', pal.blue, 'MarkerFaceColor', pal.blue, 'MarkerSize', 5);
plot(ax2, [2 2], [18 20], '-', 'Color', pal.blue_l, 'LineWidth', 3);
plot(ax2, 2, 19, 'o', 'Color', pal.blue_l, 'MarkerFaceColor', pal.blue_l, 'MarkerSize', 5);
% cumulative-case caliber
plot(ax2, 3:5, wks, '-D', 'Color', pal.red, 'MarkerFaceColor', pal.red, ...
    'LineWidth', 1.0, 'MarkerSize', 4.5);
for k = 1:3
    text(ax2, k + 2, wks(k) + 0.9, sprintf('%d', wks(k)), 'Color', pal.red, ...
        'FontSize', 7.2, 'HorizontalAlignment', 'center');
end
xlim(ax2, [0.5 5.5]); ylim(ax2, [15 37]); yticks(ax2, 15:2:37);
xticks(ax2, 1:5);
xticklabels(ax2, [{'T0 lag 1', 'T0 lag 2'}, ...
    cellfun(@(v) sprintf('%d case%s', v, repmat('s', 1, double(v ~= 1))), ...
            num2cell(thr(:)'), 'UniformOutput', false)]);
xlabel(ax2, 'Caliber');
ylabel(ax2, 'Warning week');
text(ax2, 1.5, 16.0, 'T_0 crossing (this work)', 'FontSize', 7, 'Color', pal.blue, ...
    'HorizontalAlignment', 'center');
text(ax2, 4.0, 16.0, 'cumulative-case threshold', 'FontSize', 7, 'Color', pal.red, ...
    'HorizontalAlignment', 'center');
plot(ax2, [2.6 2.6], [15 37], ':', 'Color', pal.lgray, 'HandleVisibility', 'off');
text(ax2, 4.0, 31.5, 'drift 16 weeks', 'FontSize', 7.2, 'Color', pal.red, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center');
text(ax2, 1.5, 31.5, 'drift 1-2 weeks', 'FontSize', 7.2, 'Color', pal.blue, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center');
title(ax2, '(b) Threshold sensitivity: T_0 crossing vs cumulative-case caliber', 'FontSize', 8.5, 'FontWeight', 'normal');

% ===== panel c =====
ax3 = nexttile(tl); hold(ax3, 'on');
plot(ax3, rho(ix1), cum(ix1), '-o', 'Color', pal.blue, 'MarkerSize', 2.4, 'LineWidth', 1.0);
plot(ax3, rho(ix2), cum(ix2), '--s', 'Color', pal.red, 'MarkerSize', 2.4, 'LineWidth', 1.0);
ylim(ax3, [0 3.4]); xlim(ax3, [0 1.02]);
xlabel(ax3, 'Detection ratio \rho');
ylabel(ax3, 'Cumulative detections at \tau^* (cases)');
text(ax3, 0.02, 3.05, '\leq 2.7 cumulative detections at all triggers: no lead time', ...
    'FontSize', 7.2, 'Color', [0.35 0.35 0.35]);
title(ax3, '(c) Honesty panel: observed cases at the warning week', 'FontSize', 8.5, 'FontWeight', 'normal');

% ---- source data export ----
sd = fullfile(d.src, 'fig1_c4_warning_week');
if ~exist(sd, 'dir'); mkdir(sd); end
wcsv(fullfile(sd, 'panel_a_warning_week_vs_rho.csv'), ...
    {'channel', 'branch_lag', 'rho', 'tau_star_week', 'in_frozen_grid'}, ...
    {ch, lag, rho, tau, infz});
wcsv(fullfile(sd, 'panel_b_threshold_comparison.csv'), ...
    {'caliber', 'label', 'week_lo', 'week_hi'}, ...
    { {'T0', 'T0', 'cum', 'cum', 'cum'}, ...
      {'lag 1', 'lag 2', '1 case', '100 cases', '500 cases'}, ...
      [19; 18; wks], ...
      [20; 20; wks] });
wcsv(fullfile(sd, 'panel_c_cum_detected_at_tau.csv'), ...
    {'branch_lag', 'rho', 'cum_detected_at_tau_star'}, ...
    {lag([ix1; ix2]), rho([ix1; ix2]), cum([ix1; ix2])});

% ---- read-only proof ----
pp_util('assert', snap);
pp_util('savefig', fig, 'fig1_c4_warning_week');
fprintf('fig1: done.\n');
end

% =========================================================================
function wcsv(path, names, cols)
fid = fopen(path, 'w', 'n', 'UTF-8');
assert(fid > 0, 'fig1:write', 'cannot write %s', path);
fprintf(fid, '%s', strjoin(names, ','));
fprintf(fid, '\n');
n = numel(cols{1});
for k = 1:n
    for j = 1:numel(cols)
        c = cols{j};
        if iscell(c)
            v = c{k};
            fprintf(fid, '%s', v);
        else
            v = c(k);
            if isnan(v)
                fprintf(fid, 'NA');
            elseif abs(v - round(v)) < 1e-9
                fprintf(fid, '%d', round(v));
            else
                fprintf(fid, '%.10g', v);
            end
        end
        if j < numel(cols); fprintf(fid, ','); end
    end
    fprintf(fid, '\n');
end
fclose(fid);
end
