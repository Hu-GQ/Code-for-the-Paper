# SVIR prevention compliance model

Code and simulated data for **An SVIR epidemic model with an evolutionary game of prevention compliance under vaccination**, by Guoqing Hu, Zunhai Gao and Yakui Xue.

## Reproduce

The tested Python environment is recorded in `ENVIRONMENT.json`; pinned dependencies are in `requirements.txt`. Run these commands from `paper-1` or use the script paths from another working directory:

```bash
python src/verify_stability.py
python src/check_equilibrium_count.py
python src/robustness.py
python src/plot_manuscript_figures.py
```

The symbolic checks verify the quartic characteristic polynomial, both Hurwitz identities and the reduced cubic factorisation. The same script checks 783 feasible pinned equilibria. The equilibrium-count script verifies the attained bound of ten and 1000 additional parameter sets.

The robustness script runs 55 SVIR integrations with initial-condition and solver controls. The figure script reads archived baseline curves and efficacy data, then integrates two reduced VIR cases to 20000 days. Those two cases are additional to the 55 SVIR cases. Generated plots are written under `figures/`. A full run can take several minutes.

The behavioural derivatives use the exact logit identity `dz/dt = gain(x)*(I-theta)`. Probabilities and their complements are evaluated separately to preserve small tails. The archived largest infected-level error in the SVIR verification set is 1.80334e-8; the largest physical-equation residual is 2.80626e-8.

## Files

| Location | Contents |
|---|---|
| `src/` | Model, analytic Jacobian, symbolic checks, equilibrium enumeration, simulations and plotting |
| `data/original/` | 18 baseline trajectory and threshold-sweep CSV files used by the figures |
| `results/` | 16 result files supporting the efficacy figures, equilibrium counts and numerical checks |
| `proofs/` | Exported symbolic identities |
| `MANIFEST_SHA256.json` | File hashes, excluding the manifest itself |

## Figure mapping

| Paper item | Data |
|---|---|
| Figure 1 | `data/original/rev_P51.csv` |
| Figure 2 | `rev_fig5a_*.csv`, `rev_fig5b_*.csv`, `results/vir_long_time_theta_*.csv` |
| Figure 3 | `data/original/bif_thetas.csv` |
| Figure 4 | Analytic boundaries; `results/fig7_checkpoints_corrected.csv` |
| Figure 5 | `results/efficacy_branches.csv`, `fig8_timeseries_theta_*.csv` |
| Supplementary S1 | `data/original/rev_P110.csv` |
| Supplementary S2 | `data/original/rev_P112.csv` |
| Supplementary S3 | `rev_fig4a_*.csv`, `rev_fig4b_*.csv` |

Baseline data can be regenerated with `julia src/recompute_all.jl` and `julia src/bifurcation_scan.jl`. The recorded baseline environment is Julia 1.10 with DifferentialEquations.jl 7.16. These Julia programs were not rerun for the submission preparation; terminal equilibria were checked independently with the Python implementation. Some generator outputs are diagnostics that are not needed in the deposited dataset.

All data are simulated. Parameters illustrate this model under contact reassignment and are not fitted to observations. The stability results concern local equilibria on the stated open branches; numerical sampling does not establish global attraction. MIT licence applies to the code; generated data have no additional restrictions.
