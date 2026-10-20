function D = mwspm_load_data(root, modelCard)
%MWSPM_LOAD_DATA  Load the analysis-window forcing and observations.
%
%   D = mwspm_load_data(root, modelCard)
%
%   Reads ONLY the project's own raw/processed data (all read-only):
%     data/raw/upstream_cleaned/guangzhou_2025_daily_weather_openmeteo.csv
%     data/raw/upstream_cleaned/guangzhou_2025_daily_cases_reconstructed.csv
%     data/raw/upstream_cleaned/guangzhou_population.csv
%     data/processed/weekly_analysis_input.csv
%
%   Contract (data/metadata/DATA_DICTIONARY.md):
%     * week_index 1..47 starts at 2025-02-03; week_index = iso_week - 5.
%     * the DAILY case series is a uniform within-week reallocation of the
%       weekly counts, NOT real daily data.  It is used only as a forcing of
%       I_h (and for monthly aggregation / plotting); no daily-scale
%       inference may be based on it.
%     * a missing value means UNAVAILABLE, not zero.  This loader therefore
%       ASSERTS that the analysis window contains no missing temperature or
%       importation value, instead of silently filling zeros.
%
%   Returned struct D:
%     dates (329x1 datetime), t_days (329x1, 0..328 offset in days),
%     T_daily (329x1), Lambda_daily (329x1, imported cases),
%     Y_weekly (47x1 reported local), Yimp_weekly (47x1), week_start,
%     N_h, n_days, n_weeks, diagnostics

t0 = datetime(modelCard.window.t0_date);
t1 = datetime(modelCard.window.t_end_date);

% ---- daily temperature ---------------------------------------------------
W = mwspm_readcsv(fullfile(root, modelCard.forcing.T_daily.file));
wdate = datetime(mwspm_csvcol(W, 'date'));
Tmean = mwspm_csvcol(W, 'temperature_mean_c', 'numeric');

% ---- daily importation (reconstructed) -----------------------------------
C = mwspm_readcsv(fullfile(root, modelCard.forcing.Lambda_daily.file));
cdate = datetime(mwspm_csvcol(C, 'date'));
imp   = mwspm_csvcol(C, 'imported_reported_reconstructed', 'numeric');

% ---- population ----------------------------------------------------------
P = mwspm_readcsv(fullfile(root, 'data/raw/upstream_cleaned/guangzhou_population.csv'));
pop = mwspm_csvcol(P, 'population', 'numeric');
N_h = pop(1);
assert(isfinite(N_h) && N_h > 0, 'mwspm:load', 'population could not be read');

% ---- weekly observations -------------------------------------------------
O  = mwspm_readcsv(fullfile(root, modelCard.forcing.observations_weekly.file));
wk = mwspm_csvcol(O, 'week_index', 'numeric');
Y  = mwspm_csvcol(O, 'local_reported', 'numeric');
Yi = mwspm_csvcol(O, 'imported_reported', 'numeric');
ws = mwspm_csvcol(O, 'week_start');

% ---- window selection on the daily grid ----------------------------------
mask = (wdate >= t0) & (wdate <= t1);
dates   = wdate(mask);
T_daily = Tmean(mask);

[ok, loc] = ismember(dates, cdate);
assert(all(ok), 'mwspm:load', 'daily case file does not cover all window dates');
Lambda_daily = imp(loc);

assert(numel(dates) == modelCard.window.n_days, 'mwspm:load', ...
    'window has %d days, model card declares %d', numel(dates), modelCard.window.n_days);
assert(~any(isnan(T_daily)), 'mwspm:load', '%d window days have missing temperature', sum(isnan(T_daily)));
assert(~any(isnan(Lambda_daily)), 'mwspm:load', ...
    '%d window days have missing importation (missing means UNAVAILABLE, not zero)', sum(isnan(Lambda_daily)));

assert(numel(wk) == modelCard.window.n_weeks, 'mwspm:load', ...
    'weekly observation file has %d rows, model card declares %d', numel(wk), modelCard.window.n_weeks);
assert(isequal(wk(:)', 1:modelCard.window.n_weeks), 'mwspm:load', 'week_index is not 1..47');
assert(~any(isnan(Y)) && ~any(isnan(Yi)), 'mwspm:load', 'missing weekly observation');

D = struct();
D.t0_date      = modelCard.window.t0_date;
D.t_end_date   = modelCard.window.t_end_date;
D.dates        = dates;
D.t_days       = (0:numel(dates)-1)';                 % 0 .. 328
D.doy          = day(dates, 'dayofyear');
D.T_daily      = T_daily;
D.Lambda_daily = Lambda_daily;
D.N_h          = N_h;
D.n_days       = numel(dates);
D.n_weeks      = numel(wk);
D.week_index   = wk(:)';
D.week_start   = ws;
D.Y_weekly     = Y(:);
D.Yimp_weekly  = Yi(:);

D.diagnostics = struct();
D.diagnostics.temperature_min_c      = min(T_daily);
D.diagnostics.temperature_max_c      = max(T_daily);
D.diagnostics.import_total           = sum(Lambda_daily);
D.diagnostics.import_total_weekly    = sum(Yi);
D.diagnostics.observed_local_total   = sum(Y);
D.diagnostics.observed_peak_week     = find(Y == max(Y), 1);
D.diagnostics.n_missing_temperature  = 0;
D.diagnostics.n_missing_import       = 0;
D.diagnostics.daily_cases_are_reconstructed = true;
D.diagnostics.weekly_index_convention = 'week_index = iso_week - 5; week_index 1 starts 2025-02-03';

assert(abs(D.diagnostics.import_total - D.diagnostics.import_total_weekly) < 1e-9, ...
    'mwspm:load', 'daily and weekly importation totals disagree');
end
