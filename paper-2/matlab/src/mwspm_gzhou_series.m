function mwspm_gzhou_series(root)
%MWSPM_GZHOU_SERIES  Export the the analysis-caliber daily series for cor:gzhou.
%
%   mwspm_gzhou_series()          auto-locate the project root
%   mwspm_gzhou_series(root)      explicit root
%
%   PURPOSE
%   -------
%   cor:gzhou (theory/theorems.tex) reports three criteria computed from
%   code/step20_nonautonomous.py, which hard-codes output/results/refit_v2.json
%   (legacy refit, m0=1.62154, mu_mult=12.42018, rho=1, d=2) AND, on top of that,
%   uses the daily-corrected mortality variant thermal_traits_day instead of the
%   continuous 'raw' variant that R0(t) / the envelope / Floquet must use.
%
%   This helper recomputes the exact daily series under the the analysis MAIN REPORT
%   POINT (d=1) using the authoritative mwspm_season / mwspm_thermal_traits('raw')
%   / mwspm_r0 caliber, and writes it to
%       theory/derivations/gzhou_t1201_daily_series.csv
%   for step20_gzhou_t1201.py to consume.  It changes NO production code and
%   touches no experiment run; it is a read-only derivation aid.
%
%   The recomputed max R0 must equal the the analysis baseline R0_max = 3.5569
%   (run 20261007T155558_86f11a38, parent fit 20261006T012933_c2b62344).

if nargin < 1 || isempty(root)
    root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
end
assert(exist(fullfile(root, 'config', 'project.json'), 'file') == 2, ...
    'mwspm:gzhou', 'project root not found at %s', root);
addpath(fullfile(root, 'src', 'matlab'), '-begin');

modelCard = mwspm_readjson(fullfile(root, 'config', 'experiments', ...
    'model_mwspm_gen3_rhohat.json'));
fitPath   = fullfile(root, 'results', 'runs', '20261006T012933_c2b62344', ...
    'data', 'fit_best.json');
fit = mwspm_readjson(fitPath);
th  = fit.main_report_point.theta;
assert(double(th.d_weeks) == 1, 'mwspm:gzhou', ...
    'expected the the analysis main report point (d=1), got d=%d', double(th.d_weeks));

D = mwspm_load_data(root, modelCard);
p = mwspm_params(modelCard);

sp = struct('onset_doy', double(th.onset_doy), 'offset_doy', double(th.offset_doy), ...
    'kappa_days', p.kappa_days, 'gamma', double(th.gamma), 'beta', double(th.beta), ...
    's_floor', p.s_floor, 'min_window_days', p.min_window_days, ...
    'normalize', true, 'norm_factor', NaN);
s = mwspm_season(D.doy, sp);

tr     = mwspm_thermal_traits(D.T_daily, 'raw');
muEff  = double(th.mu_s) * tr.mu_raw;
R0     = mwspm_r0(D.T_daily, double(th.m0) * s, muEff, p.gamma_h);

T = table((1:D.n_days)', D.doy, D.T_daily, s, double(th.m0) * s, ...
    tr.mu_raw, muEff, tr.a, tr.b, tr.c, tr.sigma_v, R0, ...
    'VariableNames', {'day', 'doy', 'T_daily', 's', 'm', 'mu_raw', 'mu_eff', ...
                      'a', 'b', 'c', 'sigma_v', 'R0'});

outDir  = fullfile(root, 'theory', 'derivations');
if exist(outDir, 'dir') ~= 7, mkdir(outDir); end
outPath = fullfile(outDir, 'gzhou_t1201_daily_series.csv');
writetable(T, outPath);

[mxR0, iMx] = max(R0);
fprintf('MWSPM_GZHOU_SERIES wrote %s\n', outPath);
fprintf('  n_days = %d\n', D.n_days);
fprintf('  theta  : m0=%.12f mu_s=%.12f onset=%.9f offset=%.9f gamma=%.9f beta=%.9f d=%d\n', ...
    double(th.m0), double(th.mu_s), double(th.onset_doy), double(th.offset_doy), ...
    double(th.gamma), double(th.beta), double(th.d_weeks));
fprintf('  max R0 = %.10f at day %d (week_index %d)\n', mxR0, iMx, floor((iMx - 1) / 7) + 1);
fprintf('  mbar = max m = %.6f\n', max(double(th.m0) * s));
end
