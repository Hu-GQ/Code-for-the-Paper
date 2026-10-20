function fig6_traits_robustness()
%FIG6_TRAITS_ROBUSTNESS  Fig 6 - trait-source robustness of the warning week.
%
%   Figure-level statement: tau* shifts by 2 weeks under trait-set replacement
%   (20 -> 18), reproduced INDEPENDENTLY by two swaps (same-source species
%   swap the analysis; cross-source hybrid the analysis); the fit cannot discriminate
%   between trait sets (dBIC 0.82 / 4.26).
%
%   Read-only data:
%     the analysis  20261007T155558_86f11a38/data/r0T_curves.csv
%     the analysis  20261007T205217_442e99f5/data/trait_set_refits.csv
%     the analysis  20261007T210949_913918bd/data/trait_set_refits.csv
%   Claims: the statement (limits 3/4), the statement (limit).

pp_util('style');
d = pp_util('dirs');
pal = pp_util('palette');

runB = pp_util('runDir', '20261007T155558_86f11a38');
runC = pp_util('runDir', '20261007T205217_442e99f5');
runD = pp_util('runDir', '20261007T210949_913918bd');
fB = fullfile(runB, 'data', 'r0T_curves.csv');
snap = pp_util('snap', {fB});

% ===== panel a: R0(T) family under coefficient perturbation (the analysis) =====
C = mwspm_readcsv(fB);
TC = pp_util('col', C, 'T_C');
ncol = numel(C.names) - 1;
M = NaN(C.n_rows, ncol);
for k = 1:ncol
    M(:, k) = pp_util('col', C, C.names{k + 1});
end
base = M(:, 1);
env_lo = min(M, [], 2);
env_hi = max(M, [], 2);

% ===== tau* under trait-set replacement, revised specification =====
% Two fixed reporting ratios (0.8 and 0.45) so that the onset-week difference
% measures the trait curves rather than drift along the flat rho ridge.
J8 = mwspm_readjson(fullfile(pp_util('runDir', '20261018T175739_67747f71'), 'data', 'trait_summary.json'));
J45 = mwspm_readjson(fullfile(pp_util('runDir', '20261018T175933_d8109d33'), 'data', 'trait_summary.json'));
t8 = double(J8.tau_star(:)); b8 = double(J8.delta_bic_vs_baseline(:));
t45 = double(J45.tau_star(:)); b45 = double(J45.delta_bic_vs_baseline(:));

fig = figure('Units', 'inches', 'Position', [0.5 0.5 7.10 3.40]);
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, 'Onset week under trait-set replacement (revised specification)', ...
    'FontWeight', 'bold', 'FontSize', 9.5);

% panel a
ax1 = nexttile(tl); hold(ax1, 'on');
fill(ax1, [TC; flipud(TC)], [env_lo; flipud(env_hi)], [0.85 0.88 0.93], ...
    'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(ax1, TC, base, '-', 'Color', pal.blue, 'LineWidth', 1.4);
xlim(ax1, [10 37]); ylim(ax1, [0 4.1]);
xlabel(ax1, 'Temperature (deg C)');
ylabel(ax1, 'R_0(T)');
text(ax1, 11, 3.72, 'baseline: Ae. albopictus (Mordecai 2017 Table A)', 'FontSize', 7.0, 'Color', pal.blue);
text(ax1, 11, 3.38, 'envelope: 30 single-coefficient perturbations', ...
    'FontSize', 6.8, 'Color', [0.4 0.4 0.4]);
text(ax1, 11, 3.04, '(coefficient +/-10%; T0/Tm +/-1C; 5 traits)', ...
    'FontSize', 6.8, 'Color', [0.4 0.4 0.4]);
title(ax1, '(a) R_0(T) under coefficient perturbation', 'FontSize', 8.5, 'FontWeight', 'normal');

% panel b: tau* at two fixed reporting ratios (dot plot, no zero baseline)
ax2 = nexttile(tl); hold(ax2, 'on');
labs = {'baseline', 'aegypti swap', 'hybrid source'};
cols3 = [pal.blue; pal.orange; pal.green];
for k = 1:3
    plot(ax2, [k-0.16 k+0.16], [t8(k) t45(k)], '-', 'Color', [0.75 0.75 0.75], ...
        'LineWidth', 0.8, 'HandleVisibility', 'off');
    plot(ax2, k-0.16, t8(k), 'o', 'MarkerSize', 6.5, 'MarkerFaceColor', cols3(k, :), ...
        'MarkerEdgeColor', 'k', 'LineWidth', 0.3, 'HandleVisibility', 'off');
    plot(ax2, k+0.16, t45(k), 's', 'MarkerSize', 6, 'MarkerFaceColor', 'w', ...
        'MarkerEdgeColor', cols3(k, :), 'LineWidth', 1.1, 'HandleVisibility', 'off');
    text(ax2, k, 20.45, sprintf('%d / %d', t8(k), t45(k)), 'HorizontalAlignment', 'center', ...
        'FontSize', 8, 'FontWeight', 'bold', 'Color', cols3(k, :));
end
h1 = plot(ax2, NaN, NaN, 'o', 'MarkerSize', 6.5, 'MarkerFaceColor', pal.blue, ...
    'MarkerEdgeColor', 'k', 'LineWidth', 0.3);
h2 = plot(ax2, NaN, NaN, 's', 'MarkerSize', 6, 'MarkerFaceColor', 'w', ...
    'MarkerEdgeColor', pal.blue, 'LineWidth', 1.1);
lg = legend(ax2, [h1 h2], {'\rho = 0.8 (fixed)', '\rho = 0.45 (fixed)'}, ...
    'Location', 'north', 'Orientation', 'horizontal', 'FontSize', 7);
lg.ItemTokenSize = 10;
set(ax2, 'XTick', 1:3, 'XTickLabel', labs, 'FontSize', 7.5);
xlim(ax2, [0.4 3.6]); ylim(ax2, [16.2 21.6]);
ylabel(ax2, 'Onset week \tau^* (traits refitted)');
text(ax2, 2.0, 16.35, [char(916) 'BIC vs baseline ' char(8804) ' ' ...
    num2str(max(abs([b8(2:end); b45(2:end)])), '%.1f') ...
    ': the fit cannot separate the sets'], ...
    'FontSize', 6.9, 'Color', [0.35 0.35 0.35], 'HorizontalAlignment', 'center');
title(ax2, '(b) \tau^* by trait set at two fixed reporting ratios', 'FontSize', 8.5, 'FontWeight', 'normal');

% ---- source data ----
sd = fullfile(d.src, 'fig6_traits_robustness');
if ~exist(sd, 'dir'); mkdir(sd); end
wsave(fullfile(sd, 'panel_a_r0T_envelope.csv'), {'T_C', 'R0T_baseline', 'R0T_env_lo', 'R0T_env_hi'}, ...
    {TC, base, env_lo, env_hi});
wsave(fullfile(sd, 'panel_b_tau_star.csv'), ...
    {'set_index', 'tau_star_rho08', 'tau_star_rho045', 'delta_bic_rho08', 'delta_bic_rho045'}, ...
    {[1; 2; 3], t8, t45, b8, b45});

pp_util('assert', snap);
pp_util('savefig', fig, 'fig6_traits_robustness');
fprintf('fig6: done.\n');
end

% =========================================================================
function wsave(path, names, cols)
fid = fopen(path, 'w', 'n', 'UTF-8');
assert(fid > 0, 'fig6:write', 'cannot write %s', path);
fprintf(fid, '%s\n', strjoin(names, ','));
for k = 1:numel(cols{1})
    for j = 1:numel(cols)
        v = cols{j}(k);
        if iscell(v)
            fprintf(fid, '%s', v{1});
        elseif isnan(v)
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
