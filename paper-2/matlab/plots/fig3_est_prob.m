function fig3_est_prob()
%FIG3_EST_PROB  Fig 3 - probability of at least one local detection before window end.
%
%   Figure-level claim (frozen): the p_est profile is NOT single-peaked
%   (4 interior local maxima) with peak weeks 20/21 - the earlier single-peak
%   claim (CL-09) was refuted by the higher-intensity stochastic control.
%   Estimator comparison: p_carmona / p_surv are of the same order as each
%   other but a different order than p_est.
%
%   Read-only data: T13.02  20261006T105433_a29b2ba7/data/
%       branching_analytic_grid.csv, p_establishment_by_week.csv
%   Claims: CL-09 (updated caliber), CL-01 (mechanism context).

pp_util('style');
d = pp_util('dirs');
pal = pp_util('palette');

runA = pp_util('runDir', '20261006T105433_a29b2ba7');
fA = fullfile(runA, 'data', 'branching_analytic_grid.csv');
fB = fullfile(runA, 'data', 'p_establishment_by_week.csv');
snap = pp_util('snap', {fA, fB});

G = mwspm_readcsv(fA);
var  = pp_util('strcol', G, 'variant');
wk   = pp_util('col', G, 'tau_week');
rho  = pp_util('col', G, 'rho');
pest = pp_util('col', G, 'p_est');
psur = pp_util('col', G, 'p_surv');
pcar = pp_util('col', G, 'p_carmona');

B = mwspm_readcsv(fB);
bvar = pp_util('strcol', B, 'variant');
brho = pp_util('col', B, 'rho');
bpeak = pp_util('col', B, 'peak_week');
blmax = pp_util('strcol', B, 'local_max_weeks');

isday = strcmp(var, 'day');
rhos = [0.13; 0.30; 0.50; 0.80; 1.00];
nR = numel(rhos);

% sequential light-to-dark blue for the rho family
c1 = [0.78 0.86 0.96]; c5 = [0.05 0.20 0.55];
cols = interp1([1; nR], [c1; c5], (1:nR)', 'linear');

fig = figure('Units', 'inches', 'Position', [0.5 0.5 7.10 6.40]);
tl = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, 'Probability of at least one local detection before window end', ...
    'FontWeight', 'bold', 'FontSize', 9.5);

% ===== panel a =====
ax1 = nexttile(tl); hold(ax1, 'on');
h = gobjects(nR, 1);
for k = 1:nR
    sel = isday & abs(rho - rhos(k)) < 1e-6;
    h(k) = plot(ax1, wk(sel), pest(sel), '-', 'Color', cols(k,:), 'LineWidth', 1.2, ...
        'DisplayName', sprintf('\\rho = %.2f', rhos(k)));
end
% mark the four interior local maxima of the rho = 0.80 profile
i80 = find(strcmp(bvar, 'day') & abs(brho - 0.80) < 1e-6, 1);
lmax = sscanf(blmax{i80}, '%d')';
sel80 = isday & abs(rho - 0.80) < 1e-6;
for w = lmax
    yv = pest(sel80 & wk == w);
    plot(ax1, w, yv, 'v', 'Color', cols(4,:), 'MarkerFaceColor', 'w', ...
        'MarkerSize', 5, 'LineWidth', 0.9, 'HandleVisibility', 'off');
end
xline(ax1, 19, '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 0.9, 'HandleVisibility', 'off');
text(ax1, 19.6, 0.02, 'week 19: first detected cases', 'FontSize', 6.8, ...
    'Color', [0.45 0.45 0.45]);
xlim(ax1, [0 48]); ylim(ax1, [0 0.88]);
xlabel(ax1, 'Introduction week (analysis week 1-47)');
ylabel(ax1, 'p_{det}(w)');
lg = legend(ax1, h, 'Location', 'northwest', 'FontSize', 7, 'NumColumns', 2);
lg.ItemTokenSize = 10;
text(ax1, 24.5, 0.80, 'not single-peaked: 4 interior local maxima', ...
    'FontSize', 7.5, 'Color', cols(4,:), 'FontWeight', 'bold');
text(ax1, 24.5, 0.735, 'peak weeks 20/21 (markers: local maxima of \rho = 0.80)', ...
    'FontSize', 7.0, 'Color', [0.35 0.35 0.35]);
title(ax1, '(a) p_{det}(w) by introduction week, five reporting ratios', 'FontSize', 8.5, 'FontWeight', 'normal');

% ===== panel b: estimator comparison (rho = 0.80) =====
ax2 = nexttile(tl); hold(ax2, 'on');
sel80 = isday & abs(rho - 0.80) < 1e-6;
wkv = wk(sel80);
psv = psur(sel80);
pcv = pcar(sel80);
plot(ax2, wkv, max(pest(sel80), 1e-12), '-', 'Color', pal.blue, 'LineWidth', 1.3, ...
    'DisplayName', 'p_{det} (this work, \rho = 0.80)');
plot(ax2, wkv, max(psv, 1e-12), '--', 'Color', pal.green, 'LineWidth', 1.2, ...
    'DisplayName', 'p_{surv} (survival-type)');
plot(ax2, wkv, max(pcv, 1e-12), ':', 'Color', pal.orange, 'LineWidth', 1.4, ...
    'DisplayName', 'p_{carmona} (Carmona & Gandon 2020)');
xlim(ax2, [0 48]); ylim(ax2, [1e-12 2]);
set(ax2, 'YScale', 'log');
xlabel(ax2, 'Introduction week (analysis week 1-47)');
ylabel(ax2, 'Probability (log scale)');
lg2 = legend(ax2, 'Location', 'northwest', 'FontSize', 7); lg2.ItemTokenSize = 12;
text(ax2, 10, 2.2e-11, {'p_{surv}: survival-type comparator', ...
    'p_{carmona}: Carmona & Gandon (2020) single-type formula', ...
    'neither matches the endpoint of p_{det}'}, ...
    'FontSize', 6.8, 'Color', [0.35 0.35 0.35]);
title(ax2, '(b) Different endpoints give different quantities (log scale)', 'FontSize', 8.5, 'FontWeight', 'normal');

% ---- source data ----
sd = fullfile(d.src, 'fig3_est_prob');
if ~exist(sd, 'dir'); mkdir(sd); end
sel = isday;
wcsv(fullfile(sd, 'panel_a_est_prob_day.csv'), ...
    {'tau_week', 'rho', 'p_est'}, {wk(sel), rho(sel), pest(sel)});
wcsv(fullfile(sd, 'panel_b_estimators_rho080.csv'), ...
    {'tau_week', 'p_est', 'p_surv', 'p_carmona'}, ...
    {wkv, pest(sel80), psv, pcv});
copyfile(fB, fullfile(sd, 'p_establishment_by_week.csv'));

pp_util('assert', snap);
pp_util('savefig', fig, 'fig3_est_prob');
fprintf('fig3: done.\n');
end

% =========================================================================
function wcsv(path, names, cols)
fid = fopen(path, 'w', 'n', 'UTF-8');
assert(fid > 0, 'fig3:write', 'cannot write %s', path);
fprintf(fid, '%s\n', strjoin(names, ','));
n = numel(cols{1});
for k = 1:n
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
