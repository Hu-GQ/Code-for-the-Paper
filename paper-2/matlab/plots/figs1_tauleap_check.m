function figs1_tauleap_check()
%FIGS1_TAULEAP_CHECK  Fig S1 - tau-leap simulation vs analytic p_est (T13.02).
%
%   Read-only: T13.02 20261006T105433_a29b2ba7/data/{branching_sim_subset.csv,
%   branching_dt_refinement.csv}
%   Claims: CL-09 (simulation consistency of the frozen detection-probability curve).
%
%   Panel (a): p_est simulation vs exact (log-log, 1:1 line).
%   Panel (b): dt-refinement change |fine - coarse| in units of the coarse-run
%              MC se for the 6 (tau, rho) nodes x 2 variants.

pp_util('style');
d = pp_util('dirs');
pal = pp_util('palette');

runA = pp_util('runDir', '20261006T105433_a29b2ba7');
fA = fullfile(runA, 'data', 'branching_sim_subset.csv');
fB = fullfile(runA, 'data', 'branching_dt_refinement.csv');
snap = pp_util('snap', {fA, fB});

S = mwspm_readcsv(fA);
var = pp_util('strcol', S, 'variant');
ex  = pp_util('col', S, 'p_est_exact');
sim = pp_util('col', S, 'p_est_sim');
ds  = pp_util('col', S, 'abs_diff_in_mc_se');
isday = strcmp(var, 'day');
ok = sim > 0 & ex > 0;

fig = figure('Units', 'inches', 'Position', [0.5 0.5 4.90 6.55]);

% ===== panel a: simulation vs exact =====
ax = subplot(2, 1, 1); hold(ax, 'on');
lo = 1e-9; hi = 1.2;
plot(ax, [lo hi], [lo hi], '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 0.9, 'HandleVisibility', 'off');
plot(ax, ex(ok & isday),  sim(ok & isday),  'o', 'Color', pal.blue, 'MarkerFaceColor', pal.blue, ...
    'MarkerSize', 3.4, 'DisplayName', 'day variant');
plot(ax, ex(ok & ~isday), sim(ok & ~isday), 's', 'Color', pal.orange, 'MarkerFaceColor', pal.orange, ...
    'MarkerSize', 3.4, 'DisplayName', 'raw variant');
set(ax, 'XScale', 'log', 'YScale', 'log');
xlim(ax, [lo hi]); ylim(ax, [lo hi]);
xlabel(ax, 'p_{est} analytic (exact)');
ylabel(ax, 'p_{est} tau-leap (N = 10^4)');
lg = legend(ax, 'Location', 'southeast', 'FontSize', 7); lg.ItemTokenSize = 12;
maxd = max(ds(~isnan(ds)));
nok = sum(ok); ntot = numel(var);
text(ax, 2e-9, 0.085, sprintf('max |sim - exact| = %.2f MC se', maxd), ...
    'FontSize', 7, 'Color', [0.35 0.35 0.35]);
text(ax, 2e-9, 0.020, sprintf('%d of %d points with \\geq 1 event shown', nok, ntot), ...
    'FontSize', 6.8, 'Color', [0.45 0.45 0.45]);
title(ax, '(a) tau-leap vs analytic detection probability', 'FontSize', 8.5, 'FontWeight', 'bold');

% ===== panel b: dt-refinement change =====
T = mwspm_readcsv(fB);
var2 = pp_util('strcol', T, 'variant');
t2 = pp_util('col', T, 'tau_week');
r2 = pp_util('col', T, 'rho');
pc = pp_util('col', T, 'p_est_coarse');
pf = pp_util('col', T, 'p_est_fine');
sc = pp_util('col', T, 'mc_se_coarse');
node = [t2 r2];
u = unique(node, 'rows');                 % 6 (tau, rho) nodes, sorted
ax2 = subplot(2, 1, 2); hold(ax2, 'on');
w = 0.34;
for k = 1:size(u, 1)
    for v = 1:2
        if v == 1; sel = strcmp(var2, 'day'); col = pal.blue; else; sel = strcmp(var2, 'raw'); col = pal.orange; end
        m = sel & abs(node(:, 1) - u(k, 1)) < 1e-9 & abs(node(:, 2) - u(k, 2)) < 1e-9;
        yv = abs(pf(m) - pc(m)) / sc(m);
        xb = k + (v - 1.5) * w;
        bar(ax2, xb, yv, w * 0.92, 'FaceColor', col, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    end
end
yline(ax2, 1, '--', 'Color', pal.red, 'LineWidth', 1.0);
text(ax2, 6.45, 1.10, '1 MC se', 'FontSize', 6.8, 'Color', pal.red, 'HorizontalAlignment', 'right');
lbl = cell(1, size(u, 1));
for k = 1:size(u, 1)
    lbl{k} = sprintf('%g / %.2f', u(k, 1), u(k, 2));
end
set(ax2, 'XTick', 1:size(u, 1), 'XTickLabel', lbl, 'FontSize', 7);
xlabel(ax2, 'tau week / rho of the refinement nodes');
ylabel(ax2, '|fine - coarse| / MC se');
xlim(ax2, [0.4, size(u, 1) + 0.6]); ylim(ax2, [0, 1.7]);
text(ax2, 0.55, 1.62, 'dt-refinement change stays at the MC-error scale', ...
    'FontSize', 7, 'Color', [0.35 0.35 0.35]);
text(ax2, 0.55, 1.49, '(coarse dt is already adequate for this check)', ...
    'FontSize', 6.8, 'Color', [0.45 0.45 0.45]);
title(ax2, '(b) dt-refinement: fine vs coarse p_{est}', 'FontSize', 8.5, 'FontWeight', 'bold');

sd = fullfile(d.src, 'figs1_tauleap_check');
if ~exist(sd, 'dir'); mkdir(sd); end
wsave(fullfile(sd, 'sim_vs_exact.csv'), {'variant', 'tau_week', 'rho', 'p_est_exact', 'p_est_sim', 'abs_diff_in_mc_se'}, ...
    {var, pp_util('col', S, 'tau_week'), pp_util('col', S, 'rho'), ex, sim, ds});
wsave(fullfile(sd, 'dt_refinement.csv'), {'variant', 'tau_week', 'rho', 'p_est_coarse', 'p_est_fine', 'mc_se_coarse'}, ...
    {var2, t2, r2, pc, pf, sc});

pp_util('assert', snap);
pp_util('savefig_supp', fig, 'figs1_tauleap_check');
fprintf('s1: done.\n');
end

function wsave(path, names, cols)
fid = fopen(path, 'w', 'n', 'UTF-8');
fprintf(fid, '%s\n', strjoin(names, ','));
for k = 1:numel(cols{1})
    for j = 1:numel(cols)
        c = cols{j};
        if iscell(c)
            fprintf(fid, '%s', c{k});
        else
            v = c(k);
            if isnan(v); fprintf(fid, 'NA');
            elseif abs(v - round(v)) < 1e-9; fprintf(fid, '%d', round(v));
            else; fprintf(fid, '%.10g', v); end
        end
        if j < numel(cols); fprintf(fid, ','); end
    end
    fprintf(fid, '\n');
end
fclose(fid);
end
