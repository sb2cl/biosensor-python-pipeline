# BioSensor Python Pipeline — v0.2 model notebooks

**Open-source research tools for SB2CL / ai2 / UPV — MIT.**

This repository removes MATLAB as a bottleneck for calibrated Excel extraction,
BioSensor Feature Studio interchange, frozen posterior/PPC plotting, and
mechanistic FdeR **forward simulation**. It follows the pedagogical structure
of the lab's `Simulation_Gene_Expression_pyLab.ipynb` and rigorously separates
**three different mechanistic model lineages**. It neither fits new posterior
parameters nor replaces the archived JBE inference.

## Three notebook pathways

| Notebook | Scientific role | Public Binder usable? |
|---|---|---|
| `01_excel_to_feature_studio.ipynb` | Synthetic Excel -> normalized per-dose replicate CSV for the existing [Feature Studio](https://sb2cl.github.io/biosensor-feature-studio/) | Yes |
| `02_frozen_matlab_ppc_figures.ipynb` | Frozen 11D posterior/PPC plotting, no refit | Yes, but skips private-data plots until the original files are supplied |
| **`03_minimal_FdeR_RBS.ipynb`** | 3-state activating model, qTF/qOUT sweeps, EC50 shift; *illustrative, not fitted* | **Yes** |
| **`04_final_FdeR_MATLAB_port.ipynb`** | Direct 5-state 11D forward model with `solve_ivp(BDF)` and correlated 6D/3D likelihood | **Yes** on a fabricated synthetic context; requires private `.mat` for the real experiment |
| **`05_historical_2025_candidate.ipynb`** | Historical 2025 de-repression candidate, without falsely attributing design forecasts | **Yes**, hypothetical data and context only |

### Why the 2025 model and the minimal model are distinct

The recovered 2025 MATLAB model used a **de-repression** PfdeO response, in
which a stronger TF increases the naringenin dose needed to release the
promoter. The later FdeR ligand-**activation** core can reproduce the
observed *lower* EC50 when the TF RBS is stronger. These are **not
interchangeable models**. The original 2025 forecast files and RBS-to-parameter
mapping have not been recovered; we will not retroactively call a new model
the original design-run. See [model lineage and numerical validation](docs/MODEL_LINEAGE_AND_VALIDATION.md).

## Frozen final 11D MATLAB -> Python: validation result

For the exact archived input
`Data_pFdeR_final_batchaware_VBMC_20260921_214836.mat` we ran an independent
Python port with the **same posterior median and initial seed**. Four design
contexts each were simulated. The **24** archived basal/DR/EC50 feature values
had a maximum relative discrepancy **6.68×10⁻⁶**. The two archived
log-posteriors differed at most **4.96×10⁻⁴**. No refit was performed.

That is strong forward-model numerical agreement, not bitwise identity nor a
new posterior. Actual private MAT and all archived experimental data must
stay outside the public repo. The shipped synthetic demonstration is not
experimental data.

## Run without MATLAB

```bash
python -m pip install -r requirements.txt
python -m pytest tests -q
python -m scripts.plot_fder_models --out outputs/demo_models
```

In JupyterLab, open `notebooks/03_minimal_FdeR_RBS.ipynb` first. Experiment
with `MinimalParameters`, sweep the RBS of the TF and the output separately,
then use `04_final_FdeR_MATLAB_port.ipynb` to see the exact 2026 model
structure. Running notebook 04 without private data selects a **synthetic
context**, labelled clearly in the plots and text.

### Run frozen MATLAB final model locally

Place the original file at
`private_data/Data_pFdeR_final_batchaware_VBMC_20260921_214836.mat` (excluded
from git), then execute notebook 04, or run the private numerical audit:

```bash
python -m scripts.validate_final_fder \
    --mat private_data/Data_pFdeR_final_batchaware_VBMC_20260921_214836.mat \
    --out outputs/private_parity --with-seed
python -m scripts.plot_fder_final_parity \
    --csv outputs/private_parity/fder_final_matlab_python_parity.csv \
    --out outputs/private_parity
```

The archived MAT contains 11D parameter samples, experimental growth
contexts, observation window and covariance targets. The existing Python
`load_final_model` reads them directly; this avoids silently substituting
fabricated experimental conditions.

## Cloud options

**GitHub Codespaces (preferred)**: GitHub repository -> **Code -> Codespaces ->
Create codespace**, then open the three notebooks. The `.devcontainer`
installs requirements. Codespaces can be used with a private repo.

**MyBinder (public demo only)**: after publishing the public repo, open

`https://mybinder.org/v2/gh/sb2cl/biosensor-python-pipeline/HEAD?urlpath=lab/tree/notebooks/03_minimal_FdeR_RBS.ipynb`

Change the organization/repo if you chose different names. Binder is
short-lived; use Codespaces, a lab workstation or a cluster for computational
inference or private experimental data.

**GitHub Actions** runs tests on synthetic examples after each push.
Do not confuse GitHub Pages (static site), GitHub Actions (automation), and
Codespaces (interactive remote Python environment).

## Experimental Excel and Feature Studio

The user-provided Excel source is still required for true parsing/calibration
validation. `01_excel_to_feature_studio.ipynb` and `scripts.extract_excel`
consume **already calibrated MEFL/cell** time courses from a synthetic Excel
workbook and create Feature Studio-compatible normalized CSV. They do not
apply a raw Cytation calibration or alter historically recorded doses. The
frozen FdeR mean window `[72:108]` (MATLAB indices) and pLux 360–480 min
window are different; see [migration gates](MIGRATION_GATES.md).

Feature extraction/bootstrapping is already supplied by the separate
[BioSensor Feature Studio](https://sb2cl.github.io/biosensor-feature-studio/).
We intentionally do not fork that implementation in this repository.

## Open-science and provenance

The code is MIT-licensed; copyright (c) 2026 Universitat Politècnica de
València. The real `.mat`, `.xlsx`, experimental growth and unpublished
posterior outputs are not part of this source bundle. The current
[BioSensor Feature Studio software DOI](https://doi.org/10.5281/zenodo.23461067)
archives that **other tool**, not this Python-model repository. The new repo
has **no DOI or GitHub release yet**.

Detailed limitations and validation: [docs/MODEL_LINEAGE_AND_VALIDATION.md](docs/MODEL_LINEAGE_AND_VALIDATION.md).

### Optional: fit only the reference biosensor (NEW model analysis)

If an experimental `feature_studio_input.csv` containing a **single** R4R3
reference condition becomes available locally, this repo offers a constrained
**three-parameter exploratory fit** to the reference (ligand half-activation,
GFP scale and basal fraction). It subsequently predicts R3R3 using the
independent TF-RBS ratio **without another fit**:

```bash
python -m scripts.calibrate_minimal_fder \
    --csv private_data/first_experiment_feature_studio_input.csv \
    --series R4R3_Exp1 --tf R4 --out outputs/private_minimal_fit
```

It produces a dose-by-dose comparison CSV plus a provenance JSON. This is
**not** a recovery of the historical forecast, and it does not fit the final
VBMC parameters. The fixed TF abundance/promoter threshold are non-identifiable
with the fitted output amplitude in general. The fitting code is a *tool*;
no private experimental file or fitted parameters are included in this repo.
