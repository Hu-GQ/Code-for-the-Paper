# Data

These are the inputs the calibration actually reads. Nothing here is hand
edited after download; `metadata/raw_file_hashes.csv` records the SHA-256 of the
files as fetched, and `metadata/processed_checks.json` records the consistency
checks applied when building `processed/weekly_analysis_input.csv`.

| File | Role |
|---|---|
| `guangzhou_2025_weekly_cases_official.csv` | weekly local and imported reports, 47 weeks |
| `guangzhou_2025_weekly_model_data.csv` | the weekly table in model order (week_index, local, imported) |
| `guangzhou_2025_daily_weather_openmeteo.csv` | daily mean temperature, ERA5 reanalysis, nearest grid point |
| `guangzhou_2025_daily_cases_reconstructed.csv` | **derived**, weekly totals spread uniformly within the week; display only |
| `guangzhou_2025_monthly_cases_observed.csv` | monthly totals, used only for the sanity check in the supplement |
| `guangzhou_population.csv` | N_h = 19 101 000 |
| `processed/weekly_analysis_input.csv` | the 47-row analysis table (Y_w, temperature, week_index) |

## Week numbering

`week_index` = ISO week − 5. The first non-zero local week is `week_index` 19
(ISO 24); the observed peak is `week_index` 37 (ISO 42). The paper always uses
`week_index`.

## Daily reconstruction caveat

`guangzhou_2025_daily_cases_reconstructed.csv` is not observed data. It exists
so that figures can show a daily trace. No likelihood in the paper uses it:
every fit is weekly.
