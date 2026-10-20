function fig4_rho_profile()
%FIG4_RHO_PROFILE  Fig 4 - rho profile likelihood: boundary and (un)identifiability.
%
%   Figure-level claim: the Wilks 95% band for rho covers the WHOLE grid
%   ([0.02, 1.00] at lag 1) - rho is completely unidentifiable; the reported
%   total is pinned (766.6-767.5) while the true burden spreads (767 -> 38373).
%
%   Read-only data: T12.01  20261006T012933_c2b62344/data/rho_profile.csv
%   Claims: CL-12, CL-03.

pp_util('style');
d = pp_util('dirs');
pal = pp_util('palette');

runA = pp_util('runDir', '20261006T012933_c2b62344');
fA = fullfile(runA, 'data', 'rho_profile.csv');
snap = pp_util('snap', {fA});

S = mwspm_readcsv(fA);
lag  = pp_util('col', S, 'branch_lag');
rho  = pp_util('col', S, 'rho');
dnll = pp_util('col', S, 'dnll');
rep  = pp_util('col', S, 'reported_total');
tru  = pp_util('col', S, 'true_total');

i1 = lag == 1; i2 = lag == 2;

fig = figure('Units', 'inches', 'Position', [0.5 0.5 5.30 6.60]);
tl = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% ===== panel a: profile likelihood =====
ax1 = nexttile(tl); hold(ax1, 'on');
yline(ax1, 1.9207, '--', 'Color', pal.red, 'LineWidth', 1.0, 'HandleVisibility', 'off');
plot(ax1, rho(i1), dnll(i1), '-o', 'Color', pal.blue, 'MarkerSize', 2.4, ...
    'LineWidth', 1.1, 'DisplayName', 'lag 1 (n = 99)');
plot(ax1, rho(i2), dnll(i2), '--s', 'Color', pal.red, 'MarkerSize', 2.4, ...
    'LineWidth', 1.1, 'DisplayName', 'lag 2 (n = 99)');
xlim(ax1, [0 1.01]); ylim(ax1, [0 5.2]);
xlabel(ax1, 'Detection ratio \rho');
ylabel(ax1, '\Delta nll = J(\rho) - min J');
lg = legend(ax1, 'Location', 'northwest', 'FontSize', 7); lg.ItemTokenSize = 12;
text(ax1, 0.97, 2.10, 'Wilks 95% threshold = 1.92', 'FontSize', 7.2, 'Color', pal.red, ...
    'HorizontalAlignment', 'right');
text(ax1, 0.97, 4.75, 'profile below the Wilks threshold (lag 1)', ...
    'FontSize', 7.5, 'Color', pal.blue, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
text(ax1, 0.97, 4.35, '=> 95% band = whole grid [0.02, 1.00]', ...
    'FontSize', 7.0, 'Color', [0.35 0.35 0.35], 'HorizontalAlignment', 'right');
text(ax1, 0.97, 3.95, '\rho is unidentifiable (ridge with m_0)', ...
    'FontSize', 7.0, 'Color', [0.35 0.35 0.35], 'HorizontalAlignment', 'right');
title(ax1, '(a) Profile likelihood: no interior optimum for \rho', 'FontSize', 8.5, 'FontWeight', 'normal');

% ===== panel b: reported pinned vs true spreading =====
ax2 = nexttile(tl); hold(ax2, 'on');
yyaxis(ax2, 'left');
plot(ax2, rho(i1), rep(i1), '-o', 'Color', pal.blue, 'MarkerSize', 2.4, 'LineWidth', 1.1);
ylabel(ax2, 'Reported total (cases)');
ax2.YColor = pal.blue;
ylim(ax2, [760 775]);
text(ax2, 0.30, 772.0, 'reported total pinned at 766.6-767.5 (the \rho-m_0 ridge)', ...
    'FontSize', 7.2, 'Color', pal.blue);
yyaxis(ax2, 'right');
plot(ax2, rho(i1), tru(i1), '-s', 'Color', pal.orange, 'MarkerSize', 2.4, 'LineWidth', 1.1);
set(ax2, 'YScale', 'log');
ylabel(ax2, 'True infections (log)');
ax2.YColor = pal.orange;
ylim(ax2, [500 200000]);
xlim(ax2, [0 1.01]);
xlabel(ax2, 'Detection ratio \rho');
text(ax2, 0.30, 8e2, 'true burden spreads 767 -> 38373 (lag 1)', ...
    'FontSize', 7.2, 'Color', pal.orange);
title(ax2, '(b) Boundary value is a reporting-caliber artifact, not information', ...
    'FontSize', 8.5, 'FontWeight', 'normal');

% ---- source data ----
sd = fullfile(d.src, 'fig4_rho_profile');
if ~exist(sd, 'dir'); mkdir(sd); end
wsave(fullfile(sd, 'rho_profile_lag1.csv'), ...
    {'rho', 'dnll', 'reported_total', 'true_total'}, {rho(i1), dnll(i1), rep(i1), tru(i1)});

pp_util('assert', snap);
pp_util('savefig', fig, 'fig4_rho_profile');
fprintf('fig4: done.\n');
end

% =========================================================================
function wsave(path, names, cols)
fid = fopen(path, 'w', 'n', 'UTF-8');
assert(fid > 0, 'fig4:write', 'cannot write %s', path);
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
