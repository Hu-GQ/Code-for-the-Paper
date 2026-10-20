function figs2_model_race()
%FIGS2_MODEL_RACE  Fig S2 - model-family BIC ladder (10 structures).
%
%   Data source: output/results/step8_candidates.json - a HISTORICAL save of
%   the early Python pipeline (step8_candidates.py). Read as-is and labelled
%   as such on the figure; NOT re-derived in MATLAB (the MATLAB-side family
%   comparison lives in calibration.md with the same ordering).
%   Claims: CL-10.

pp_util('style');
d = pp_util('dirs');
pal = pp_util('palette');

fJ = fullfile(pp_util('root'), 'output', 'results', 'step8_candidates.json');
snap = pp_util('snap', {fJ});

fh = fopen(fJ, 'r', 'n', 'UTF-8');
raw = fread(fh, '*char').';
fclose(fh);
J = jsondecode(raw);
rows = J.rows;
if ~iscell(rows)
    rows = arrayfun(@(x) x, rows, 'UniformOutput', false);   % struct array -> cell
end
ispo = cellfun(@(r) strcmp(r.objective, 'poisson'), rows);
po = rows(ispo);

tag = cellfun(@(r) string(r.tag), po);
bic = cellfun(@(r) r.BIC, po);
kk  = cellfun(@(r) r.k, po);
r2  = cellfun(@(r) r.R2_log, po);

[~, ord] = sort(bic, 'descend');           % y=1: worst -> top: best
tag = tag(ord); bic = bic(ord); kk = kk(ord); r2 = r2(ord);
dbic = bic - min(bic);                     % >= 0, best = 0
n = numel(bic);

fig = figure('Units', 'inches', 'Position', [0.5 0.5 5.60 4.20]);
ax = axes(fig); hold(ax, 'on');

for k = 1:n
    if strcmp(tag(k), 'M-WSPM')
        col = pal.blue; 
    else
        col = [0.78 0.78 0.78];
    end
    barh(ax, k, dbic(k), 0.62, 'FaceColor', col, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    if dbic(k) > 0
        text(ax, dbic(k) + 25, k, sprintf('\\DeltaBIC %d (k=%d, R^2_{log}=%.3f)', ...
            round(dbic(k)), kk(k), r2(k)), 'FontSize', 6.6, 'Color', [0.30 0.30 0.30], ...
            'VerticalAlignment', 'middle');
    else
        text(ax, dbic(k) + 25, k, sprintf('best (k=%d, R^2_{log}=%.3f)', kk(k), r2(k)), ...
            'FontSize', 6.8, 'Color', pal.blue, 'FontWeight', 'bold', 'VerticalAlignment', 'middle');
    end
end

xline(ax, 10, '--', 'Color', pal.red, 'LineWidth', 0.9);
text(ax, 15, 0.50, '\DeltaBIC = 10', 'FontSize', 6.4, 'Color', pal.red);
set(ax, 'YTick', 1:n, 'YTickLabel', cellstr(tag), 'FontSize', 7.5);
xlim(ax, [-60, 3900]); ylim(ax, [0.15, n + 0.6]);
xlabel(ax, '\DeltaBIC vs best (smaller = better; BIC of best = -4109.0)');
title(ax, 'Model-race BIC ladder (10 structures)', 'FontSize', 9, 'FontWeight', 'bold');
text(ax, 1650, 6.45, 'source: historical Python save', ...
    'FontSize', 6.8, 'Color', [0.40 0.40 0.40]);
text(ax, 1650, 5.90, '(step8\_candidates.json); not re-derived', ...
    'FontSize', 6.8, 'Color', [0.40 0.40 0.40]);
text(ax, 1650, 5.40, 'on the MATLAB side', ...
    'FontSize', 6.8, 'Color', [0.40 0.40 0.40]);

sd = fullfile(d.src, 'figs2_model_race');
if ~exist(sd, 'dir'); mkdir(sd); end
fid = fopen(fullfile(sd, 'bic_ladder.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, 'tag,BIC,delta_BIC,k,R2_log\n');
for k = n:-1:1
    fprintf(fid, '%s,%.3f,%.3f,%d,%.4f\n', tag(k), bic(k), dbic(k), kk(k), r2(k));
end
fclose(fid);

pp_util('assert', snap);
pp_util('savefig_supp', fig, 'figs2_model_race');
fprintf('s2: done.\n');
end
