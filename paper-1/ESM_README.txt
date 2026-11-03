Reproducibility archive

An SVIR epidemic model with an evolutionary game of prevention compliance under vaccination
Version: bmb-reproducibility-2026-11-03

This archive contains the model code, simulated data, symbolic checks,
solver controls and figure-to-data mappings. No empirical data were used.
Numerical data and computational scripts retain their original contents.
The README provides reproduction commands; ENVIRONMENT.json records the
software versions. MANIFEST_SHA256.json records the extracted file hashes.

Figure 3 uses the archived 48-point scan at day 2000. The vertical axis is
I(2000), and the markers describe finite-time terminal samples. Some samples
near collisions may remain transient. Run python BMB_Figure3.py from the
extracted directory to reproduce the current figure. The original plotting
script remains available for provenance.

The six regime-diagram checkpoints use an initial horizon of 60000 days.
See README.md for the remaining reproduction details and licence information.
