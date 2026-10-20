function fig_r2_onset_profile()
%FIG_R2_ONSET_PROFILE  Revised Fig. 4 and Fig. 5 from the pure-Poisson profile.
%
%   Replaces the archived Euler versions of
%     fig1_c4_warning_week.pdf  (manuscript Fig. 4, onset index)
%     fig4_rho_profile.pdf      (manuscript Fig. 5, reporting-ratio profile)
%   with curves computed on the revised specification: daily zero-order-held
%   coefficients, unclamped RK4 at 0.025 d, pure Poisson likelihood and no
%   peak-week penalty (run recorded in figs source_data/r2_profile_run.txt).
%
%   Claims: the statement, the statement, the statement, the statement.  All numbers are descriptive
%   conditional candidates; no interval in these figures is a confidence set.
%
%   Read-only data:
%     <r2run>/data/rho_profile_pure.csv          (revised profile)
%     20261006T012933_c2b62344/data/weekly_fit.csv   (observed weekly counts)
%     20261006T020359_3295cf28/data/threshold_drift.csv (1/100/500 case weeks)

pp_util('style');
dr  = pp_util('dirs');
pal = pp_util('palette');
root = pp_util('root');

runTxt = fullfile(dr.src, 'r2_profile_run.txt');
fid = fopen(runTxt, 'r');
assert(fid > 0, 'r2fig:noRun', 'missing %s', runTxt);
runId = strtrim(fgetl(fid));
fclose(fid);

fProf = fullfile(pp_util('runDir', runId), 'data', 'rho_profile_pure.csv');
fWeek = fullfile(pp_util('runDir', '20261006T012933_c2b62344'), 'data', 'weekly_fit.csv');
fThr  = fullfile(pp_util('runDir', '20261006T020359_3295cf28'), 'data', 'threshold_drift.csv');
assert(exist(fProf, 'file') == 2, 'r2fig:noProfile', 'missing %s', fProf);
snap = pp_util('snap', {fProf, fWeek, fThr});

P = mwspm_readcsv(fProf);
lag  = pp_util('col', P, 'lag');
rho  = pp_util('col', P, 'rho');
nll  = pp_util('col', P, 'nll');
cum  = pp_util('col', P, 'cum_reports_at_tau');
tau  = pp_util('col', P, 'tau_star');
rep  = pp_util('col', P, 'reported_total');
tru  = pp_util('col', P, 'true_total');

i1 = lag == 1; i2 = lag == 2;
dnll = nll - min(nll(i1));
[~, k1] = min(nll(i1));
rhoBest = rho(i1);
rhoBest = rhoBest(k1);
tauBest = tau(i1);
tauBest = tauBest(k1);

W = mwspm_readcsv(fWeek);
wObs = pp_util('col', W, 'observed_local');
cObs = cumsum(wObs);
wk   = (1:numel(wObs))';

T = mwspm_readcsv(fThr);
thr = pp_util('col', T, 'threshold_cumulative_cases');
twk = pp_util('col', T, 'first_week_reaching');

firstLocal = find(wObs > 0, 1);

% =========================================================================
% Fig. 4 - retrospective onset index
% =========================================================================
fig = figure('Units', 'inches', 'Position', [0.5 0.5 7.00 7.50]);
tl = tiledlayout(fig, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile(tl); hold(ax1, 'on');
patch(ax1, [0.30 1.00 1.00 0.30], [14 14 22 22], [0.92 0.92 0.92], ...
    'EdgeColor', 'none', 'HandleVisibility', 'off');
patch(ax1, 'XData', [0.02 1.00 1.00 0.02], 'YData', [14 14 22 22], 'FaceColor', 'none', ...
    'EdgeColor', [0.80 0.80 0.80], 'LineStyle', '--', 'LineWidth', 0.8, ...
    'HandleVisibility', 'off');
plot(ax1, rho(i1), tau(i1), '-o', 'Color', pal.blue, 'MarkerSize', 2.6, ...
    'MarkerFaceColor', pal.blue, 'LineWidth', 1.0, 'DisplayName', 'lag 1');
plot(ax1, rho(i2), tau(i2), '--s', 'Color', pal.red, 'MarkerSize', 2.6, ...
    'MarkerFaceColor', pal.red, 'LineWidth', 1.0, 'DisplayName', 'lag 2');
plot(ax1, rho(i1), repmat(tauBest, sum(i1), 1), ':', 'Color', pal.gray, ...
    'LineWidth', 1.0, 'DisplayName', sprintf('fixed trajectory (\\tau^* = %d)', tauBest));
yline(ax1, firstLocal, '-.', 'Color', pal.green, 'LineWidth', 0.9, ...
    'Label', sprintf('first local report (week %d)', firstLocal), ...
    'LabelHorizontalAlignment', 'right', 'FontSize', 7, 'HandleVisibility', 'off');
ylim(ax1, [15 21]); yticks(ax1, 15:21); xlim(ax1, [0 1.02]);
xlabel(ax1, 'Reporting ratio \rho');
ylabel(ax1, 'Onset week \tau^*');
lg = legend(ax1, 'Location', 'se', 'FontSize', 7); lg.ItemTokenSize = 12; lg.Color = 'w';
text(ax1, 0.985, 20.55, sprintf(['pre-specified grid [0.30,1.00] (shaded): weeks %d--%d on both branches'], ...
    min(tau(rho >= 0.30)), max(tau(rho >= 0.30))), ...
    'FontSize', 7, 'Color', [0.25 0.25 0.25], 'HorizontalAlignment', 'right');
title(ax1, '(a) Retrospective onset week vs reporting ratio', 'FontSize', 8.5, 'FontWeight', 'normal');

ax2 = nexttile(tl); hold(ax2, 'on');
plot(ax2, wk, cObs, '-', 'Color', [0.55 0.55 0.55], 'LineWidth', 1.0);
for j = 1:numel(thr)
    plot(ax2, twk(j), thr(j), 'o', 'MarkerSize', 5, 'MarkerFaceColor', pal.orange, ...
        'MarkerEdgeColor', 'k', 'LineWidth', 0.4);
    if thr(j) == 1
        lbl = sprintf('first local report: week %d', twk(j));
    else
        lbl = sprintf('%g cases: week %d', thr(j), twk(j));
    end
    text(ax2, twk(j) + 0.8, thr(j) + 18, lbl, 'FontSize', 7, 'VerticalAlignment', 'bottom');
end
xlim(ax2, [0 48]); ylim(ax2, [0 820]);
xlabel(ax2, 'Analysis week');
ylabel(ax2, 'Cumulative local reports');
text(ax2, 1.0, 760, sprintf(['three different case rules span %d weeks; ' ...
    'an illustration of definition dependence, not a performance benchmark'], ...
    max(twk) - min(twk)), 'FontSize', 7, 'Color', [0.25 0.25 0.25]);
title(ax2, '(b) Cumulative case thresholds applied to the observed series', ...
    'FontSize', 8.5, 'FontWeight', 'normal');

ax3 = nexttile(tl); hold(ax3, 'on');
plot(ax3, rho(i1), cum(i1), '-o', 'Color', pal.blue, 'MarkerSize', 2.4, ...
    'LineWidth', 1.0, 'DisplayName', 'lag 1');
plot(ax3, rho(i2), cum(i2), '--s', 'Color', pal.red, 'MarkerSize', 2.4, ...
    'LineWidth', 1.0, 'DisplayName', 'lag 2');
xlim(ax3, [0 1.02]);
xlabel(ax3, 'Reporting ratio \rho');
ylabel(ax3, 'Expected reports at \tau^* (cases)');
lg3 = legend(ax3, 'Location', 'nw', 'FontSize', 7); lg3.ItemTokenSize = 12;
text(ax3, 0.985, min(cum) + 0.1 * (max(cum) - min(cum)), ...
    sprintf('range %.2f--%.2f cases: order one, reached with the whole season''s fit', ...
    min(cum), max(cum)), 'FontSize', 7, 'Color', [0.25 0.25 0.25], ...
    'HorizontalAlignment', 'right');
title(ax3, '(c) Cumulative expected reports at the crossing', 'FontSize', 8.5, 'FontWeight', 'normal');

pp_util('assert', snap);
pp_util('savefig', fig, 'fig1_c4_warning_week');

% =========================================================================
% Fig. 5 - reporting-ratio profile
% =========================================================================
fig2 = figure('Units', 'inches', 'Position', [0.5 0.5 5.30 6.60]);
tl2 = tiledlayout(fig2, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax4 = nexttile(tl2); hold(ax4, 'on');
yline(ax4, 1.9207, '--', 'Color', pal.red, 'LineWidth', 1.0, 'HandleVisibility', 'off');
plot(ax4, rho(i1), dnll(i1), '-o', 'Color', pal.blue, 'MarkerSize', 2.4, ...
    'LineWidth', 1.1, 'DisplayName', sprintf('lag 1 (%d grid points)', sum(i1)));
plot(ax4, rho(i2), dnll(i2), '--s', 'Color', pal.red, 'MarkerSize', 2.4, ...
    'LineWidth', 1.1, 'DisplayName', sprintf('lag 2 (%d grid points)', sum(i2)));
plot(ax4, rhoBest, 0, 'v', 'Color', 'k', 'MarkerSize', 5, 'MarkerFaceColor', 'w', ...
    'DisplayName', sprintf('lowest lag-1 candidate (\\rho = %.2f)', rhoBest));
xlim(ax4, [0 1.01]); ylim(ax4, [0 6.0]);
xlabel(ax4, 'Reporting ratio \rho');
ylabel(ax4, '\Delta nll from the lowest lag-1 candidate');
lg4 = legend(ax4, 'Location', 'east', 'FontSize', 7); lg4.ItemTokenSize = 12; lg4.Color = 'w';
text(ax4, 0.02, 5.85, {'Wilks cutoff 1.92: orientation only, not an interval', ...
    sprintf('lag 2 exceeds 17 at \\rho = 0.02 (off scale)')}, 'FontSize', 6.8, ...
    'Color', pal.red, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
title(ax4, '(a) Pure-Poisson profile of the reporting ratio', 'FontSize', 8.5, 'FontWeight', 'normal');

ax5 = nexttile(tl2); hold(ax5, 'on');
yyaxis(ax5, 'left');
plot(ax5, rho(i1), rep(i1), '-o', 'Color', pal.blue, 'MarkerSize', 2.4, 'LineWidth', 1.0);
ylabel(ax5, 'Expected reports (cases)');
ax5.YColor = pal.blue;
yyaxis(ax5, 'right');
plot(ax5, rho(i1), tru(i1), '-s', 'Color', pal.red, 'MarkerSize', 2.4, 'LineWidth', 1.0);
set(ax5, 'YScale', 'log');
ylabel(ax5, 'Latent local infections (cases, log scale)');
ax5.YColor = pal.red;
xlim(ax5, [0 1.01]);
xlabel(ax5, 'Reporting ratio \rho');
text(ax5, 0.03, 0.90, {'reported series pinned near the observed total', ...
    'latent burden spans two orders of magnitude'}, ...
    'Units', 'normalized', 'FontSize', 6.8, 'Color', [0.25 0.25 0.25], 'VerticalAlignment', 'top');
title(ax5, '(b) The \rho-m_0 ridge: reported versus latent totals', 'FontSize', 8.5, 'FontWeight', 'normal');

pp_util('savefig', fig2, 'fig4_rho_profile');

% ---- machine-readable source data ---------------------------------------
sd = fullfile(dr.src, 'fig1_c4_warning_week');
if ~exist(sd, 'dir'); mkdir(sd); end
wsave(fullfile(sd, 'panel_a_onset_vs_rho.csv'), ...
    {'lag', 'rho', 'tau_star', 'cum_reports_at_tau'}, {lag, rho, tau, cum});
wsave(fullfile(sd, 'panel_b_case_rules.csv'), ...
    {'threshold_cases', 'first_week_reaching'}, {thr, twk});
sd2 = fullfile(dr.src, 'fig4_rho_profile');
if ~exist(sd2, 'dir'); mkdir(sd2); end
wsave(fullfile(sd2, 'panel_a_dnll.csv'), {'lag', 'rho', 'nll', 'dnll'}, {lag, rho, nll, dnll});
wsave(fullfile(sd2, 'panel_b_ridge.csv'), ...
    {'rho', 'reported_total', 'true_total'}, {rho(i1), rep(i1), tru(i1)});

fprintf('fig_r2_onset_profile: done (run %s).\n', runId);
end

% =========================================================================
function wsave(path, names, cols)
fid = fopen(path, 'w', 'n', 'UTF-8');
assert(fid > 0, 'r2fig:write', 'cannot write %s', path);
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
