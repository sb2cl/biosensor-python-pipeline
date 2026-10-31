# Model lineage and numerical-validation boundaries

## Why there are three separate Python models

The project contains **three scientifically distinct pieces of code** that must
not be interchanged. The uploaded teaching notebook
`Simulation_Gene_Expression_pyLab.ipynb` (2024; Yadira Boada and Alejandro
Vignoni) is a *coding-style reference*: define states and parameters, write an
ODE function, integrate, and plot. Its constitutive / generic TF-promoter
examples are **not FdeR model evidence**.

| File | Status | Mechanistic law | Data claim |
|---|---|---|---|
| `models/minimal_fder.py` | **New 2026 simplified demonstration** | Naringenin activates FdeR, which activates PfdeO; 3 states | Qualitative RBS TF → EC50 direction; illustrative parameters **not fitted** |
| `models/historical_2025.py` | **Translation of recovered candidate 2025 lineage** | Naringenin releases promoter from TF-dependent repression; 5 states | Archival RBS mapping/calibration and original prediction unavailable; **not recovered original design-run** |
| `models/final_fder.py` | **Forward Python translation of 2026 frozen MATLAB likelihood** | Ligand-bound TF activates promoter; capacity, basal, compression and batch effects; 5 states, 11 inferred parameters | Quantitative comparison to frozen MATLAB **without refit** |

## 1. Minimal activating three-state model

State order in `minimal_fder.py`: `TF_nM, Nar_intracellular_uM, GFP_MEFL_per_cell`.
Extracellular dose in uM is an input, not integrated. A constant µ dilutes
all molecular species; no experimental growth curve. The reference abundance
of TF depends on the TF-RBS efficiency ratio, while output translation depends
independently on the reporter-RBS ratio.

- `qTF = kappa0_TF / kappa0_R3`.
- `qOUT = kappa0_OUT / kappa0_R3`.
- Example kappa0: R3/B0034 = 0.1458; R4/J61100 = 0.0049; R4/R3 = 0.03360768.
- TF production scales with `qTF`, GFP production scales with `qOUT`.
- Intracellular ligand undergoes first-order transport towards a *clamped*
  external dose plus dilution and degradation.
- Promoter: `fNar = I^n/(Knar^n+I^n)`;
  `Hact = (TF*fNar)^m/(Kprom^m+(TF*fNar)^m)`;
  `P = beta+(1-beta)*Hact`.
- **Absent:** gammaTF, rho0, TFcapacity, B1/B2, cK2, measured µ(t),
  explicit extracellular state, physical molecule-volume conversion,
  VBMC posterior, and experimental calibration.

The default illustrative run gives approximate finite-range operational
midpoints **R4R3 384 uM**, **R3R3 104 uM**, and **R4R4 384 uM** (calculated
on a 0–1400 uM finite grid). These values were selected to demonstrate the
observed *direction*; they are neither model-estimated paper EC50s nor a fit to
Exp1. The minimal model predicts a **higher** DR for R3R3 than R4R3,
whereas the frozen Exp1 data show the **opposite**. Therefore it cannot claim
full reproduction of leak/DR/EC50 of both biosensors without extra mechanisms.

## 2. Recovered 2025 candidate historical model

Recovered reference: the uploaded file `model_exp_biosensor.txt`, labelled
"Updated 2025 Alejandro Vignoni, Yadira Boada", and historical supporting
functions `parameters.txt`, `build_g_from_evol.txt`, `simulateDoseResponse_single.txt`,
`logPosterior_features_global_multi.txt`. The archive's
`initial_model/02_initial_model_equations.md` and
`initial_model/08_missing_archival_items.csv` confirm these are candidate
historical/development files, NOT authenticated original design-run snapshots.

Historical state order `[FdeR, Nar_in, GFP, Cells, Nar_external]`.
Its PfdeO regulation is

`P = [(kdfO*CN)^n * (kdNar+Nar)^n] / ([(kdfO*CN)^n*(kdNar+Nar)^n] + [(kdNar*FdeR)^n])`

and output `beta + (1-beta)*P`. At fixed ligand, increasing FdeR **reduces**
this promoter activity. The hypothetical RBS rescaling in notebook 05 is
therefore not equivalent to an activator model and does not reproduce the
observed leftward EC50 trend just by applying qTF.

**Do not back-assign** modern RBS qTF values as if the historical code
had used that mapping: no original archived RBS-to-pR/pg record, design-run
script/forecast, `all_together.mat`, `results_1303261822.mat`, or
`exp_Samples.mat` was recovered in the archival audit. The notebook provides
only an honest, runnable candidate implementation with synthetic starting
states and constant growth.

## 3. Final 2026 frozen MATLAB model: exact 11D structure

Original mathematical source in the private archive `scripts_final.zip`:
`additive_back_gate/logPosterior_primary_batchaware_meanwindow.m`, functions
`local_model`, `local_sim_mean_window`, `local_mvn_logkernel`, and `local_decode`.
Configuration is **read directly** from the original
`Data_pFdeR_final_batchaware_VBMC_20260921_214836.mat`. No paper data are
embedded in the public repo. The companion
`additive_back_gate/model_exp_biosensor.m` differs from the *actual* 11D
likelihood's `local_model` because it does not include the inferred `rho0_TF`.

State order (one modeled cell class):

`[FdeR (molecules), Nar_in (molecules), GFP (effective MEFL/cell units), Cells (experimental particles), Nar_external (molecules)]`.

Growth is interpolated in min^-1 from the saved Exp1 / Exp2 measured arrays.
ODE phases are numerical pre-equilibration, preinduction, and postinduction,
matching changes of the extracellular volume, measured growth inputs,
cell count initial conditions, external doses, and initial GFP.
Time is in **minutes**.

Physical parameters, in exact order:

`etaTF, kG_R3, kdNar_shared, beta, n_nar, KI_TF_rel, gamma_TF, rho0_TF, B1, B2, cK2`.

`qTF` for R4R3 = `qR4`; for R3R3 = `qR4**(1-gamma_TF)`
(**not** 1 in the frozen likelihood); for R4R4 = `qR4`.
`qOUT` = `1` for reporter R3 or `qR4` for reporter R4.
Exp2 `kdNar = kdNar_shared * cK2`.
`TFcapacity` multiplies only the inducible `Hact` component; `beta` remains
unmultiplied. `rho0_TF` contributes basal ligand-independent active TF.
B1 or B2 is added **once after** integrating each dose and calculating the
mean of MATLAB **1-based indices 72:108 inclusive**, before feature fitting.

### Statistical target

Exp1 correlated Gaussian (6 dimensions):
`ln(Leak_R4R3), ln(DR_R4R3), ln(EC50_R4R3),`
`ln(Leak_R3R3-Leak_R4R3), ln(DR_R3R3/DR_R4R3), ln(EC50_R3R3/EC50_R4R3)`.

Exp2 correlated Gaussian (3 dimensions):
`ln(Leak_R4R3), ln(EC50_R4R3), ln(Leak_R4R3-Leak_R4R4)`.

The final `scipy.optimize.least_squares` approximates MATLAB's `lsqcurvefit`
for Hill features. It retains MATLAB's **1000-point finite-range EC50
midpoint**. Correlated log kernels use the saved target means and the saved
lower Cholesky `L` (transposed from MATLAB v7.3 HDF5 storage). Hard bounds and
physical-uniform prior + log Jacobian match the archived code's structure.
No VBMC posterior fitting is run.

### Actual numerical parity result (private validation, 31 Oct 2026)

Private input: the archived original 11D VBMC MAT. Tested identical 11D
*initial seed* and *posterior median*, each with **four** construct-context
scenarios. Total **26** comparisons:

- 2 theta points × 4 contexts × 3 features (basal, DR, Hill EC50) = 24;
- plus 2 log-posterior evaluations (6D+3D + prior Jacobian).
- **Maximum relative feature error: 6.6782e-06 (0.000668%)**.
- **Maximum absolute log-posterior error: 0.00049581**.
- Median log-posterior: MATLAB **27.802648691589**, Python **27.802391186331**.
- No posterior refit, no implicit dose modification, no experimental data
  restaging, and no change to the original files.

**Important limits:** SciPy `BDF` is not identical to MATLAB `ode23t/ode23s`;
SciPy Hill fitting is not identical to `lsqcurvefit`; the full 20,000-draw
posterior is **not** replayed here. We can now simulate model predictions
from the frozen samples using Python, but **new inference** with PyBADS/PyVBMC
is a *separate* future project and must not replace the JBE fit silently.

### Private reproduction command

```bash
python -m scripts.validate_final_fder \
    --mat private_data/Data_pFdeR_final_batchaware_VBMC_20260921_214836.mat \
    --out outputs/private_parity --with-seed
```

This checks the original `.mat` and creates a CSV and JSON inside
`outputs/private_parity/`. Keep those outputs private until authorization.

## 4. Path to a simple biosensor fit

It is reasonable to calibrate the minimal 3-state activating model to the
*first* R4R3 data, and then predict the stronger-TF R3R3 construct with the
independent qTF ratio **without refitting**. However, doing so would be a **new
simplified-model analysis**, not a recovery of the original 2025 forecaster.
To proceed reproducibly requires the first experimental calibrated MEFL/cell
Excel or the corresponding processed per-dose data, measured growth context
if used, and an explicit selection of which data were available before the
original design decision. The public repository deliberately includes none.

Don't force the minimal model to fit the reduced R3R3 DR by silently adding a
new capacity penalty; that would defeat the stated objective of testing what
the simplest activation model can explain.

### Optional calibration script, not part of frozen inference

`scripts/calibrate_minimal_fder.py` reads a private Feature Studio normalized
CSV with `dose,signal` (and optionally `series,system`). It fits ONLY
`nar_half_uM`, `gfp_full_scale_mefl_per_cell`, and `basal_fraction` by
minimizing log1p response residuals for the R4R3 reference construct. The
TF reference abundance and PfdeO threshold are fixed, to reduce but not
eliminate non-identifiability. It then predicts R3R3 **without any refit**
by applying the qTF ratio. This analysis is a NEW model exercise and must be
labeled as such; it is not a reconstruction of the original 2025 design fit.
The result is meant to test whether a minimal mechanism explains the first
biosensor, and to identify where high-TF capacity terms might actually be
necessary, not to hide residuals.
