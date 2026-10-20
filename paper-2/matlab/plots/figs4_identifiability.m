function figs4_identifiability()
%FIGS4_IDENTIFIABILITY  Fig S4 - onset-gamma identifiability (T14.03 + P4).
%
%   Read-only: fine run 20261007T195150_93456b81/data/{joint_profile_2d.csv,
%   profile_1d_onset.csv, profile_1d_gamma.csv}; coarse run
%   20261007T162039_d1f013a8/data/profile_1d_onset.csv.
%   Claims: T14.03 (verdict IDENTIFIABLE; projection control still VIOLATED
%   after P4, now as a node-phase effect).
%
%   Panel (a): joint (onset, gamma) profile, log10 dJ, fine grid.
%   Panel (b): 1-D onset profiles, fine (1 d) vs coarse (4 d, superseded).
%   Panel (c): 1-D gamma profile (fine run).

pp_util('style');
d = pp_util('dirs');
pal = pp_util('palette');

rF = pp_util('runDir', '20261007T195150_93456b81');
rC = pp_util('runDir', '20261007T162039_d1f013a8');
fJ  = fullfile(rF, 'data', 'joint_profile_2d.csv');
f1o = fullfile(rF, 'data', 'profile_1d_onset.csv');
f1g = fullfile(rF, 'data', 'profile_1d_gamma.csv');
fCo = fullfile(rC, 'data', 'profile_1d_onset.csv');
snap = pp_util('snap', {fJ, f1o, f1g, fCo});

THR = 1.9207;   % Wilks chi2_1(0.95)/2 (frozen rule)

% ===== reconstruct joint grid (rows gamma-fast, onset outer) =====
J = mwspm_readcsv(fJ);
jo = pp_util('col', J, 'onset_doy');
jg = pp_util('col', J, 'gamma');
jd = pp_util('col', J, 'dJ');
uo = unique(jo); ug = unique(jg);
Z = NaN(numel(ug), numel(uo));
for k = 1:numel(jo)
    i = find(abs(uo - jo(k)) < 1e-9, 1);
    j = find(abs(ug - jg(k)) < 1e-9, 1);
    Z(j, i) = jd(k);
end

fig = figure('Units', 'inches', 'Position', [0.3 0.3 7.20 2.80]);

% ===== panel a: joint contour =====
ax1 = subplot(1, 3, 1); hold(ax1, 'on');
L = -4:0.3:4.6;
contourf(ax1, uo, ug, log10(max(Z, 1e-4)), L, 'LineColor', 'none');
contour(ax1, uo, ug, log10(max(Z, 1e-4)), [log10(THR) log10(2.9957)], ...
    'LineColor', 'k', 'LineWidth', 1.0);
plot(ax1, 183.0445, -1.599377, 'p', 'MarkerSize', 7, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k');
colormap(ax1, parula);
cb = colorbar(ax1); cb.Label.String = 'log_{10} \DeltaJ'; cb.Label.FontSize = 6.8; cb.FontSize = 6.5;
xlabel(ax1, 'onset doy'); ylabel(ax1, '\gamma');
title(ax1, '(a) joint profile (fine grid)', 'FontSize', 8, 'FontWeight', 'bold');
text(ax1, 154, -2.00, 'joint domain spans [177.0, 185.0] in onset', 'FontSize', 6.4, 'Color', [1 1 1]);
text(ax1, 154, -2.06, '19/341 points inside 95%', 'FontSize', 6.4, 'Color', [1 1 1]);

% ===== panel b: onset 1-D profiles =====
Of = mwspm_readcsv(f1o);
of = pp_util('col', Of, 'onset_doy');
df = pp_util('col', Of, 'dJ');
Oc = mwspm_readcsv(fCo);
oc = pp_util('col', Oc, 'onset_doy');
dc = pp_util('col', Oc, 'dJ');

ax2 = subplot(1, 3, 2); hold(ax2, 'on');
yl = [0.02 2000];
patch(ax2, [177.0445 186.0445 186.0445 177.0445], [yl(1) yl(1) yl(2) yl(2)], ...
    [0.88 0.91 0.97], 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(ax2, [143 223], [THR THR], '-', 'Color', pal.red, 'LineWidth', 1.0);
plot(ax2, of, max(df, yl(1)), '-o', 'Color', pal.blue, 'MarkerSize', 1.8, 'LineWidth', 1.1, ...
    'MarkerFaceColor', pal.blue);
plot(ax2, oc, max(dc, yl(1)), '--s', 'Color', [0.45 0.45 0.45], 'MarkerSize', 3.2, 'LineWidth', 0.9);
set(ax2, 'YScale', 'log');
ylim(ax2, yl); xlim(ax2, [143 223]);
xlabel(ax2, 'onset doy'); ylabel(ax2, '\DeltaJ');
text(ax2, 144.5, 0.62, 'fine: [177.0, 186.0]', 'FontSize', 6.4, 'Color', pal.blue);
text(ax2, 144.5, 0.15, 'coarse 4d: [179, 183]', 'FontSize', 6.4, 'Color', [0.45 0.45 0.45]);
text(ax2, 144.5, 0.045, '(superseded grid)', 'FontSize', 6.2, 'Color', [0.45 0.45 0.45]);
text(ax2, 144.5, 2.6, '\DeltaJ = 1.9207', 'FontSize', 6.4, 'Color', pal.red);
title(ax2, '(b) onset profile: coarse vs fine', 'FontSize', 8, 'FontWeight', 'bold');

% ===== panel c: gamma 1-D profile =====
Gp = mwspm_readcsv(f1g);
ag = pp_util('col', Gp, 'gamma');
dg = pp_util('col', Gp, 'dJ');

ax3 = subplot(1, 3, 3); hold(ax3, 'on');
patch(ax3, [-1.7594 -1.3594 -1.3594 -1.7594], [yl(1) yl(1) yl(2) yl(2)], ...
    [0.99 0.92 0.88], 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(ax3, [-2.45 -1.05], [THR THR], '-', 'Color', pal.red, 'LineWidth', 1.0);
plot(ax3, ag, max(dg, yl(1)), '-o', 'Color', pal.orange, 'MarkerSize', 3.0, 'LineWidth', 1.1, ...
    'MarkerFaceColor', pal.orange);
set(ax3, 'YScale', 'log');
ylim(ax3, yl); xlim(ax3, [-2.45 -1.05]);
xlabel(ax3, '\gamma'); ylabel(ax3, '\DeltaJ');
text(ax3, -2.43, 0.25, '95% = [-1.76, -1.36]', 'FontSize', 6.4, 'Color', [0.80 0.45 0.10]);
text(ax3, -2.43, 0.062, 'unchanged vs coarse', 'FontSize', 6.2, 'Color', [0.55 0.42 0.30]);
title(ax3, '(c) gamma profile (fine)', 'FontSize', 8, 'FontWeight', 'bold');

% ===== source data =====
sd = fullfile(d.src, 'figs4_identifiability');
if ~exist(sd, 'dir'); mkdir(sd); end
fid = fopen(fullfile(sd, 'onset_profiles.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, 'grid,onset_doy,dJ\n');
for k = 1:numel(of); fprintf(fid, 'fine,%.6f,%.6f\n', of(k), df(k)); end
for k = 1:numel(oc); fprintf(fid, 'coarse,%.6f,%.6f\n', oc(k), dc(k)); end
fclose(fid);
fid = fopen(fullfile(sd, 'gamma_profile.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, 'gamma,dJ\n');
for k = 1:numel(ag); fprintf(fid, '%.6f,%.6f\n', ag(k), dg(k)); end
fclose(fid);

pp_util('assert', snap);
pp_util('savefig_supp', fig, 'figs4_identifiability');
fprintf('s4: done.\n');
end
