# Paper 2 — A retrospective onset index for seasonal dengue

Code, inputs and derived results for

> **A retrospective onset index for seasonal dengue: identifiability and numerical
> limits in Guangzhou, China (2025)**
> Guoqing Hu, *International Journal of Biomathematics* (submitted).

The paper defines a retrospective onset week τ\* — the first analysis week in
which cumulative expected reports reach the reference level
`T0(w) = rho / log R0(w)` — fits a seasonally forced human–mosquito model to 47
weeks of Guangzhou surveillance (767 local, 78 imported cases, 2025), and
reports what the index can and cannot support: no prospective lead time, a
practically unidentifiable reporting ratio, a box-dependent profile minimum,
and an absolute burden that is not robust to trait perturbation.

Everything reported in the manuscript is produced by the code in this
directory. Numbers that belong to the superseded daily-Euler analysis are
labelled as such in the paper and are reproduced by the same entry points.

## Layout

```
paper-2/
├── README.md                  this file
├── requirements.txt           Python side (build + verification utilities)
├── matlab/
│   ├── reproduce.m            single entry point / dispatcher
│   ├── src/                   model, calibration, inference, verification
│   └── plots/                 figure scripts, incl. pp_util.m
├── python/
│   └── verify_deposit.py      checks every file against DEPOSIT.sha256
├── data/                      surveillance and weather inputs actually used
├── results/                   derived tables cited in the paper, one file per run
├── figures/
│   ├── main/                  vector figures as placed in the manuscript
│   └── source_data/           the exact numbers behind every panel
└── DEPOSIT.sha256             checksum of every file above
```

## Requirements

* **MATLAB R2025b** — the model, calibration and figure code. No toolbox beyond
  base MATLAB is required: the scripts implement their own optimiser wrapper,
  Poisson deviance, bootstrap and `makedist`-free quantiles. `ode15s` is used
  only as an independent cross-check of the fixed-step integrator.
* **Python ≥ 3.11** with `numpy`, `scipy`, `matplotlib` — only for the auxiliary
  utilities; `python/verify_deposit.py` needs nothing but the standard library.
  See `requirements.txt`.

## Reproducing the numbers

From the directory that contains `matlab/`:

```matlab
addpath('matlab'); addpath('matlab/src'); addpath('matlab/plots');
reproduce('exp_r2_pure_profile')     % reporting-ratio profile, both lags
reproduce('exp_r2_agg_sens')         % weekly-aggregation sensitivity of tau*
reproduce('exp_r2_floquet')          % Floquet criteria at the revised points
reproduce('exp_r2_trait_oat')        % 30 one-at-a-time trait perturbations
```

Each call creates `results/runs/<run_id>/` with the frozen configuration, a
manifest, a log and the output tables. `reproduce.m` with no argument lists the
available experiments. Runs are deterministic given the seed in the
configuration; the stochastic branching-process check uses one `mt19937ar`
stream per (variant, week, rho) grid point, with the seed rule recorded in the
run manifest.

## Result files and where they appear in the paper

| Deposit file | Run | Manuscript |
|---|---|---|
| `..._rho_profile_pure.csv`, `..._profile_summary.json`, `..._solver_check.json` | `20261018T171227_3dc71f2f` | Table 3, Table 5, Fig. 4, Fig. 5 |
| `..._refine_anchors.csv`, `..._refine_summary.json` | `20261018T173554_7df7f312` | §4.3 (multi-start refinement of profile anchors) |
| `..._trait_summary.json` (ρ = 0.8, ρ = 0.45) | `20261018T175739_67747f71`, `20261018T175933_d8109d33` | Table 5, Fig. 6 |
| `..._floquet_revised.csv/.json` | `20261018T201317_eb0d5375` | `cor:gzhou` in the theory appendix |
| `..._oat_variants.csv`, `..._oat_summary.json` | `20261018T201429_dd6a4d52` | Table 5, §4.4 |
| `..._agg_sensitivity.csv`, `..._agg_summary.json` | `20261018T210504_3ffd906b` | Eq. (3.2), Table 5, §4.2 |
| `..._stable_fit.json` | `20261018T162026_4d17af05` | Table 4 (conditional fit at ρ = 0.8) |

`results/RUN_INDEX.json` records the run id of every deposited file.

## Data

All inputs are public and are included as delivered (no hand editing):

| File | Source | Licence / note |
|---|---|---|
| `data/guangzhou_2025_weekly_cases_official.csv` | Guangdong Provincial CDC weekly reports | public government releases; compiled, not redistributed verbatim |
| `data/guangzhou_2025_daily_weather_openmeteo.csv` | Open-Meteo historical API (ERA5 reanalysis), grid point nearest Guangzhou | CC-BY-4.0 (Open-Meteo terms) |
| `data/guangzhou_population.csv` | Guangzhou Municipal Bureau of Statistics | public statistics |
| `data/guangzhou_2025_daily_cases_reconstructed.csv` | **derived** | weekly totals spread uniformly inside the week; used for display only, never for fitting |
| `data/processed/weekly_analysis_input.csv` | **derived** | the 47-week table the calibration actually reads |

`data/metadata/` carries the data dictionary, the processing record, the source
list and the SHA-256 of every raw file.

Two cautions stated in the paper and repeated here:

1. The **daily** case series is a reconstruction from weekly totals. All
   likelihood-based fitting uses the weekly series.
2. Two week-numbering conventions coexist in the project's history: the
   `week_index` used in the paper is the ISO week minus 5 (peak at analysis
   week 37 = ISO week 42; first non-zero at analysis week 19 = ISO week 24).
   Every "week N" in the paper uses the `week_index` convention.

## Not included

* The raw scraped HTML/JSON snapshots of the source pages (≈1.3 MB) are kept
  with the authors; the derived CSVs and their checksums are sufficient to
  reproduce every analysis.
* The full run archive (≈20 MB of manifests, logs and code snapshots for
  superseded and failed runs) is retained locally and available on request.
  Only the runs cited above are deposited.
* Journal-format figure files (EPS/TIFF) and the submission package itself.

## Known numerical limitations (these are results, not defects)

* The reporting ratio is practically unidentifiable; the profile minimum sits
  at ρ = 0.45 on the primary branch but the fitted mortality multiplier is at
  its box bound (15.0) at **all 198** profile grid points, so the location of
  that minimum is box-dependent.
* τ\* is defined on a weekly R0 that is the arithmetic mean of seven daily
  values; median or mid-week aggregation leaves the answer unchanged, weekly
  maximum advances it by one week at ρ = 0.8.
* The branching-process probabilities in §4.5 are conditional on the initial
  daily-Euler primary report point; they are not recomputed under the revised
  calibration.

## Citation

See `CITATION.cff` in the repository root.
