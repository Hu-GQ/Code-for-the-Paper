function figs3_integrator()
%FIGS3_INTEGRATOR  Fig S3 - integrator diagnostics (C5, case-study caliber).
%
%   Read-only:
%     - the analysis 20261006T011600_2a8ab00f/data/reference_weekly.csv
%     - the analysis 20261006T011600_2a8ab00f/data/grid_results_long.csv
%       (authoritative spelled-out file; grid_results.csv is the numeric-code
%        twin produced by mwspm_writecsv - noted in the run summary.json)
%     - the analysis 20261006T012855_1c1e8c7b/data/dt_sweep_results.csv
%   Claims: the statement (refuted boolean criterion; case study after S2/S4).
%
%   Panel (a): 4 ode15s reference weekly true-incidence curves.
%   Panel (b): Euler total-error vs max|1-mu_eff*dt| (Euler stability metric), 2 models x 2 variants,
%              three annotated case-study points.
%   Panel (c): dt-sweep peak error, Euler vs RK4 (convergence check).

pp_util('style');
d = pp_util('dirs');
pal = pp_util('palette');

runA = pp_util('runDir', '20261006T011600_2a8ab00f');
runB = pp_util('runDir', '20261006T012855_1c1e8c7b');
fRef = fullfile(runA, 'data', 'reference_weekly.csv');
fGr  = fullfile(runA, 'data', 'grid_results_long.csv');
fDs  = fullfile(runB, 'data', 'dt_sweep_results.csv');
snap = pp_util('snap', {fRef, fGr, fDs});

fig = figure('Units', 'inches', 'Position', [0.3 0.3 7.20 2.75]);

% ===== panel a: reference curves =====
W = mwspm_readcsv(fRef);
wed = pp_util('col', W, 'week_edge_day');
cr = pp_util('col', W, 'true_incidence_model_baseline_m0_249__variant_raw');
cd_ = pp_util('col', W, 'true_incidence_model_baseline_m0_249__variant_day');
gr = pp_util('col', W, 'true_incidence_model_mwspm_gen3__variant_raw');
gd = pp_util('col', W, 'true_incidence_model_mwspm_gen3__variant_day');

ax = subplot(1, 3, 1); hold(ax, 'on');
plot(ax, wed, gd, '-', 'Color', pal.blue, 'LineWidth', 1.2);
plot(ax, wed, gr, '-', 'Color', pal.orange, 'LineWidth', 1.0);
plot(ax, wed, cd_, '--', 'Color', [0.45 0.55 0.75], 'LineWidth', 0.9);
plot(ax, wed, cr, '--', 'Color', [0.75 0.60 0.40], 'LineWidth', 0.9);
xlim(ax, [0 330]); xlabel(ax, 'Day index'); ylabel(ax, 'True weekly incidence');
text(ax, 8, 92, 'gen3 / day', 'FontSize', 6.6, 'Color', pal.blue);
text(ax, 8, 80, 'gen3 / raw', 'FontSize', 6.6, 'Color', pal.orange);
text(ax, 8, 68, 'baseline / day', 'FontSize', 6.6, 'Color', [0.45 0.55 0.75]);
text(ax, 8, 56, 'baseline / raw', 'FontSize', 6.6, 'Color', [0.75 0.60 0.40]);
title(ax, '(a) ode15s reference solutions', 'FontSize', 8, 'FontWeight', 'bold');

% ===== panel b: euler error vs max|1-mu_eff*dt| (Euler stability metric) =====
G = mwspm_readcsv(fGr);
gmod = pp_util('strcol', G, 'model');
gvar = pp_util('strcol', G, 'mu_variant');
gme  = pp_util('col', G, 'max_mu_eff_dt');
gre  = pp_util('col', G, 'true_total_rel_error');
ism  = strcmp(pp_util('strcol', G, 'method'), 'euler');

ax2 = subplot(1, 3, 2); hold(ax2, 'on');
xline(ax2, 1, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.9);
sty = {{'mwspm_gen3', 'day', 'o', pal.blue, '-'}, {'mwspm_gen3', 'raw', 's', pal.orange, '-'}, ...
       {'baseline_m0.249', 'day', '^', [0.45 0.55 0.75], '--'}, {'baseline_m0.249', 'raw', 'd', [0.75 0.60 0.40], '--'}};
for k = 1:4
    m = ism & strcmp(gmod, sty{k}{1}) & strcmp(gvar, sty{k}{2});
    [xm, o] = sort(gme(m));
    ym = abs(gre(m)); ym = ym(o);
    plot(ax2, xm, ym, [sty{k}{5} sty{k}{3}], 'Color', sty{k}{4}, 'MarkerSize', 3.2, ...
        'MarkerFaceColor', sty{k}{4}, 'LineWidth', 1.0);
end
set(ax2, 'XScale', 'log', 'YScale', 'log');
xlim(ax2, [0.7 100]); ylim(ax2, [0.001 1000]);
xlabel(ax2, 'max|1 - \mu_{eff} \Delta t|  (Euler)');
ylabel(ax2, '|total rel. error|');
text(ax2, 1.25, 460, 'gen3 / raw: +869%', 'FontSize', 6.4, 'Color', pal.orange);
text(ax2, 1.25, 120, 'gen3 / day: -8.1%', 'FontSize', 6.4, 'Color', pal.blue);
text(ax2, 1.25, 31, 'baseline / day: -39.3%', 'FontSize', 6.4, 'Color', [0.45 0.55 0.75]);
text(ax2, 1.25, 8.0, '(0 days outside stability!)', 'FontSize', 6.2, 'Color', [0.45 0.55 0.75]);
text(ax2, 2.6, 0.0035, 'dotted: |1 - \mu_{eff}\Delta t| = 1', 'FontSize', 6.3, 'Color', [0.5 0.5 0.5]);
title(ax2, '(b) Euler dt scan (case study)', 'FontSize', 8, 'FontWeight', 'bold');

% ===== panel c: dt-sweep convergence =====
D = mwspm_readcsv(fDs);
dm = pp_util('strcol', D, 'method');
dd = pp_util('col', D, 'dt_days');
dp = abs(pp_util('col', D, 'peak_I_h_rel_error'));
ax3 = subplot(1, 3, 3); hold(ax3, 'on');
for k = 1:2
    if k == 1; m = strcmp(dm, 'euler'); mk = 'o'; col = pal.blue; else; m = strcmp(dm, 'rk4'); mk = 's'; col = pal.orange; end
    [xm, o] = sort(dd(m));
    ym = dp(m); ym = ym(o);
    plot(ax3, xm, ym, ['-' mk], 'Color', col, 'MarkerSize', 3.2, 'MarkerFaceColor', col, 'LineWidth', 1.0);
end
set(ax3, 'XScale', 'log', 'YScale', 'log');
xlim(ax3, [0.008 6]); ylim(ax3, [1e-7 1e2]);
xlabel(ax3, 'dt (days)'); ylabel(ax3, '|peak I_h rel. error|');
text(ax3, 0.0105, 0.45, 'Euler', 'FontSize', 6.6, 'Color', pal.blue);
text(ax3, 0.055, 2.2e-6, 'RK4', 'FontSize', 6.6, 'Color', pal.orange);
text(ax3, 0.013, 28, 'orders ~1 / ~4 (the analysis)', 'FontSize', 6.4, 'Color', [0.35 0.35 0.35]);
title(ax3, '(c) dt-sweep peak error', 'FontSize', 8, 'FontWeight', 'bold');

sd = fullfile(d.src, 'figs3_integrator');
if ~exist(sd, 'dir'); mkdir(sd); end
gdt = pp_util('col', G, 'dt_days');
gnd = pp_util('col', G, 'n_days_outside_euler_stability');
fid = fopen(fullfile(sd, 'euler_grid.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, 'model,mu_variant,dt_days,max_mu_eff_dt,true_total_rel_error,n_days_outside_euler_stability\n');
for k = 1:numel(gme)
    if ism(k)
        fprintf(fid, '%s,%s,%g,%.4f,%.6f,%g\n', gmod{k}, gvar{k}, gdt(k), gme(k), gre(k), gnd(k));
    end
end
fclose(fid);

pp_util('assert', snap);
pp_util('savefig_supp', fig, 'figs3_integrator');
fprintf('s3: done.\n');
end
