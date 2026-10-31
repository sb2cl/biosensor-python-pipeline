# Scientific migration gates (JBE FdeR; no new model or refit)

## Gate A — Current v0.1: exported pre-calibrated timecourses

- Verified: synthetic Excel reading, explicit well/column/dose mapping, FdeR MATLAB row-index window, pLux time window, exported replicates and provenance JSON.
- **Not yet verified**: archived actual Cytation Excel formatting, instrument calibration to MEFL/cell, experimental batch assignments, dose correction log. The actual source Excel is absent.
- Feature fitting must remain in **BioSensor Feature Studio v1.0.0** until parity established; do not accidentally use the much older `biosensorStats2_smoothed.m` as final Python method.

## Gate B — Frozen posterior, for figures immediately

- Archived `Data_pFdeR_final_batchaware_VBMC_20260921_214836.mat` physical `Xs` 11×20,000 MATLAB order; read with `h5py` as 20,000×11 Python.
- Archived PPC `P41`, `P31`, `P42`, `P44` 5×300 MATLAB order; read as 300×5 with feature order Leak, DR, EC50, Smax, Hill.
- CSV `final_batchaware_PPC_features_20260921_225249.csv` contains experimental medians, 95% intervals and feature role tags. The real data have **not** been included in the public code repository.
- Notes: R4R4 DR marginal interval is nominal RBS q. Independently propagated reporter-q uncertainty was handled in a separate MATLAB audit, not in the frozen nominal PPC `.mat`. Do not present nominal bands as full hierarchical prediction intervals.

## Gate C — Direct port of the FINAL ODE, not the six-parameter predecessor

Original implementation: `logPosterior_primary_batchaware_meanwindow.m > local_model` (in archived scripts).

- State order for Ncell=1: `[FdeR, Nar, GFP, particle_number, extracellular_Nar]`.
- At each step: measured experimental growth-rate interpolation `mu(t)` (not posterior nuisance mean), consistent initial stage with `ode23t`, followed by induction `ode23s`; use `scipy.integrate.solve_ivp` stiff methods `BDF` or `Radau` with matched tolerances and observation-grid interpolation.
- Final primary parameters: `etaTF, kG_R3, kdNar_shared, beta, n_nar, KI_TF_rel, gamma_TF, rho0_TF, B1, B2, cK2`. `n_fdeO=mI_TF=1` fixed; kdfO gauge and units from archived config.
- R3/R4 q coordinates and capacity term must match original expressions; `TFcapacity` multiplies **only the induced Hact component**, not the basal beta component.
- Observation is **mean over MATLAB indices 72:108 inclusive**, after interpolation to the experimental time grid, plus B1/B2 **exactly once**.
- Gate: replay parameters saved as `cfg.seed_physical` and `final_vbmc.posterior_median` and compare EACH design architecture dose curve to original outputs. Quantify max error, relative error, features. Solver differences can generate numerical shifts; do not select tolerance solely to improve a prediction.

## Gate D — Likelihood and prior

- Exp1 6-dimensional and Exp2 3-dimensional correlated Gaussian blocks, using target mean, Cholesky, log determinant, and covariance regularization from archived `cfg.target`.
- Uniform priors in physical bounds, with log Jacobian in z-space; gamma linear; B1/B2 `log1p` transforms.
- Gate: reproduce archived log-posterior seed=22.7784658708 and posterior median=27.8026486916, with comparable Mahalanobis block values, **without** refit. This is a separate validation from just matching ODE trajectories.

## Gate E — Optional optimization/inference

- `PyBADS`: point optimum / warm start; do not call it a full posterior.
- `PyVBMC`: re-estimate an approximate **new** posterior. The archived MATLAB VBMC `vp` is not a Python PyVBMC posterior object; do not claim that it transfers in place. Physical draws `Xs` and forward predictions do transfer.
- Require several independent seeds, runtime and stability assessment, bounds, ELBO/uncertainty, posterior correlation, held-out PPC, and reporting differences from the 2026 MATLAB fit.
- Do not silently replace the original frozen model in JBE.

## Gate C updated — 2026-10-31: frozen 11D forward port PASSED at two archived points

- Original MATLAB `local_model` and 3-stage `local_sim_mean_window` translated into `biosensor_pipeline/models/final_fder.py`.
- Private original MATLAB v7.3 MAT read with h5py (11D parameters, p0, g1/g2, 6D/3D target).
- Seed + posterior median: **24 MATLAB-vs-Python primary features** (4 contexts x 3 features x 2 points) max relative error **6.6782e-6**.
- Exact 6D/3D correlated likelihood + Jacobian: **two** log-posterior differences max absolute **0.00049581**.
- Not tested: full posterior 20,000-draw rerun, MATLAB code's bitwise solver identity, a new fit or a PyVBMC re-inference.

## Gate C-mini — 2026-10-31: three-state educational model PASSED (qualitative only)

- qTF R4R3 -> R3R3 shifts **illustrative operational EC50 ~384 -> ~104 uM**.
- qOUT R4R3 -> R4R4 scales output without shifting EC50.
- The simplified activating model gets the **observed TF-EC50 direction**, but **not the reduced R3R3 dynamic range**; that discrepancy motivates, but does not license silently adding, the later TFcapacity term.
- Recovered 2025 candidate is documented independently; its promoter law is de-repression, and historical mapping/predictions are not archived.

## Gate D updated — 2026-10-31: 6D/3D likelihood first cross-language comparison PASSED

- Physical-coordinate bounds, matrix-Cholesky transpose from MATLAB v7.3 storage, experimental blocks and Jacobian are ported.
- Recorded absolute logp discrepancies < 0.00050 for original seed and median.
- No inference with PyBADS/PyVBMC has been attempted; do not use a new posterior for the frozen JBE paper.
