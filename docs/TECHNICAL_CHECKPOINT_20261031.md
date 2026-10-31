# Python FdeR models — technical checkpoint (31 Oct 2026)

**Context:** Port of MATLAB into browser-executable research notebooks without re-estimating the frozen JBE posterior.

## Verified today

- Notebook style and structure follows user-provided 2024 teaching notebook `Simulation_Gene_Expression_pyLab.ipynb` by Boada/Vignoni, with `solve_ivp` replacing `odeint`.
- Three model tracks implemented: **minimal three-state ligand-activating**, **literal candidate 2025 de-repressor**, and **direct 2026 final 11D five-state port**.
- The 2026 MATLAB code mathematically transcribed from `logPosterior_primary_batchaware_meanwindow.m`/`local_model`, with stored growth context directly read from original MAT. Matplotlib figures, no new inference.
- Private parity against frozen 2026-09-21 MATLAB posterior: original seed and posterior median, four dose-response scenarios each. **24 primary feature comparisons** (basal, DR, Hill EC50); max relative discrepancy **6.678198728097201e-06**. **Two** complete 6D/3D+Jacobian log-posteriors; max absolute discrepancy **0.0004958088670257155**. No refit.
- Minimal 3-state no-penalty example: `R4R3 EC50* 384.26 µM` → `R3R3 EC50* 103.78 µM`, with `R4R4 EC50* ~384.26 µM` and lower reporter amplitude. Asterisks denote a finite-range operational midpoint, **not** the paper's fitted Hill EC50.
- All three new notebooks passed execution in a Python environment using **synthetic inputs only**.
- Historical 2025 candidate has a different de-repression promoter law and lacks proven archived design-run RBS mappings; do not claim its simulations are the historical R4R3→R3R3 forecast.
- A new optional reference-fit script is present but no actual experimental reference data has yet been fitted. Fitting it later must be called a NEW minimalist analysis, not a reproduction of the archived 2025 design output.

## Remaining scope before full replacement of MATLAB

1. First-experiment original Excel, calibration, well-to-construct mapping, stage times, and growth/particle data to validate the Excel extraction fully.
2. New minimum-model fit to FIRST R4R3 data only and a genuine R3R3 prediction without refit; quantify held-out EC50, DR and basal mismatches. Do not add resource/capacity terms unless scientifically justified.
3. Extend final-model numerical parity to a carefully chosen diverse subset of saved 11D posterior draws, then compare response curves and features across parameter extremes; test numerical robustness of SciPy solvers.
4. Export actual paper figures (svg/pdf) from **frozen posterior results** and preserve published calibration context. Do not overwrite MATLAB posterior.
5. Later, optional PyBADS/PyVBMC benchmarks **as new independent inference**, never as a renamed version of the frozen fit.

## Public data discipline

No private `.mat`, experimental `.xlsx`, raw flow cytometry or unpublished research CSV is included in the public GitHub bundle. Public notebook 04 defaults to a clearly fabricated growth context if the original MAT is absent; the private parity test is automatically skipped on public CI.

## Citation separation

BioSensor Feature Studio software DOI `10.5281/zenodo.23461067` identifies a separate application; the Python model repository has no DOI at this point.
