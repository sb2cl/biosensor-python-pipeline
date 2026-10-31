from pathlib import Path
import nbformat as nbf
BASE=Path('/mnt/data/bfs_models_v02/notebooks')
HEADER='''from pathlib import Path
import sys
ROOT = Path.cwd()
if not (ROOT / 'biosensor_pipeline').is_dir():
    ROOT = ROOT.parent
sys.path.insert(0, str(ROOT))
print('Project root:', ROOT)
'''

def save(name, cells):
 n=nbf.v4.new_notebook(cells=cells,metadata={'kernelspec':{'display_name':'Python 3','language':'python','name':'python3'},'language_info':{'name':'python'}})
 nbf.write(n,BASE/name)
 print(BASE/name,len(cells),'cells')

save('03_minimal_FdeR_RBS.ipynb',[
 nbf.v4.new_markdown_cell(r'''# 03 — Minimal FdeR model: TF RBS shifts EC50

**Status: mechanistic *illustration*, not a calibration, historical design-run reconstruction or the VBMC 11D model.**

This notebook follows the style of the supplied `Simulation_Gene_Expression_pyLab.ipynb`: define a small ODE function, list parameters and initial conditions, integrate with SciPy, sweep a genetic parameter, then plot dose-response and time courses.

The three states are **FdeR [nM]**, **intracellular naringenin [µM]**, **GFP [illustrative MEFL/cell]**. The extracellular naringenin dose and genetic RBS multipliers `qTF`, `qOUT` are **external inputs**, not states. Units of time: **minutes**.

It is specifically the *ligand-activated FdeR* core of the later mechanistic model, without the eventual additions `rho0_TF`, `gamma_TF`, `TFcapacity`, batch offsets `B1/B2`, or context-specific `cK2`.

## Equations

$$\frac{dR}{dt}=k_R q_{TF}-(d_R+\mu)R,\qquad
\frac{dI}{dt}=k_{in}(N_e-I)-(d_I+\mu)I$$

$$ f_N=\frac{I^n}{K_N^n+I^n},\quad
P_{fdeO}=\beta+(1-\beta)\frac{R f_N}{K_P+R f_N} $$

$$\frac{dG}{dt}=k_G q_{OUT}P_{fdeO}-(d_G+\mu)G.$$

The code parameterizes $k_R$ using the reference R3 steady TF concentration, and $k_G$ using the reference GFP scale. These are illustrative, not fitted. The RBS ratios come from the independent biopart measurements. Note that *finite-dose operational EC50* in this teaching example is **not** the final paper Hill estimator.'''),
 nbf.v4.new_code_cell(HEADER),
 nbf.v4.new_markdown_cell(r'''## 1. Define parameters and state derivative function

The reusable implementation is in `biosensor_pipeline/models/minimal_fder.py`. This is the same design pattern as the provided PyLab example, using the contemporary `solve_ivp` signature `f(t, state, ...)` rather than legacy `odeint(f, state0, t, ...)`.'''),
 nbf.v4.new_code_cell('''import numpy as np
import matplotlib.pyplot as plt
from biosensor_pipeline.models.minimal_fder import (
    MinimalParameters, minimal_rhs, promoter_activity,
    simulate_minimal_one, simulate_minimal_dose_response, normalized_half_response)
p = MinimalParameters()
print('Reference TF at R3:', p.tf_R3_steady_nM, 'nM')
print('qTF (R4/R3):', p.q_tf('R4'))
print('GFP reporter scaling (R4/R3):', p.q_out('R4'))
print('PfdeO activity at zero naringenin:', promoter_activity(400, 0, p))'''),
 nbf.v4.new_markdown_cell(r'''## 2. Simulate three genetic architectures

The naming order is **TF RBS / reporter RBS**. Strong TF: R3 = B0034. Weak TF: R4 = J61100. R3R3 and R4R3 differ in *TF* translation; R4R3 and R4R4 differ in *reporter* translation.'''),
 nbf.v4.new_code_cell('''doses = np.linspace(0, 1400, 41)
fig, ax = plt.subplots(figsize=(8.8, 5.0))
rows=[]
for name in ('R4R3','R3R3','R4R4'):
    x,y = simulate_minimal_dose_response(doses,name[:2],name[2:],p)
    e = normalized_half_response(x,y)
    rows.append((name,round(e,1),round(y[0],1),round(np.ptp(y),1)))
    ax.plot(x,y,label=f'{name} | EC50* ≈ {e:.0f} µM',linewidth=2)
    ax.axvline(e,color=ax.lines[-1].get_color(),lw=.8,ls=':')
ax.set(xlabel='External naringenin (µM)',ylabel='Illustrative GFP (MEFL/cell)',
       title='Minimal activation-only FdeR model (not fitted)')
ax.legend(loc='upper right');ax.grid(alpha=.2)
plt.show()
print('Construct | EC50* [µM] | Basal [MEFL/cell] | DR [MEFL/cell]')
for row in rows:print(row)
print('* EC50 is a finite-range midpoint here, not the published Hill EC50.')'''),
 nbf.v4.new_markdown_cell(r'''## 3. Time courses at the same inducer concentration

This confirms which state changes when the TF RBS changes, independently of the reporter RBS.'''),
 nbf.v4.new_code_cell('''fig, axs = plt.subplots(1, 3, figsize=(13, 3.7), layout='constrained')
for name in ('R4R3','R3R3'):
    sol = simulate_minimal_one(200,tf_rbs=name[:2],out_rbs=name[2:],p=p)
    for i,ax in enumerate(axs):ax.plot(sol.t/60,sol.y[i],label=name)
for ax,label in zip(axs,['FdeR (nM)','Nar_in (µM)','GFP (MEFL/cell)']):
    ax.set(xlabel='Time after induction (h)',ylabel=label);ax.grid(alpha=.25)
axs[0].legend();plt.show()'''),
 nbf.v4.new_markdown_cell(r'''## 4. What this model can and cannot explain

The stronger TF RBS **decreases the effective dose threshold** in this activating model. Weakening the reporter RBS **scales GFP output** without shifting the EC50. However, without additional mechanisms, the strong TF condition also has a greater inducible capacity; this is **opposite to the archived experimental DR trend**. The model is therefore **not a quantitative reproduction of both biosensors**.

For actual calibration, a private normalized experimental CSV must be supplied and the fit would need to be labelled as a *new minimal-model fit*, not the historical forecast. Do not place unpublished experiments in the public repository.'''),
 nbf.v4.new_code_cell("print('Next: 04_final_FdeR_MATLAB_port.ipynb for frozen 11D simulation; 05_historical_2025_candidate.ipynb for the distinct older de-repression law.')")
])

save('04_final_FdeR_MATLAB_port.ipynb',[
 nbf.v4.new_markdown_cell(r'''# 04 — FdeR final 11-parameter MATLAB model → SciPy

**Forward simulation and likelihood validation only. NO refit, PyBADS or PyVBMC.**

The source is the private original `logPosterior_primary_batchaware_meanwindow.m`: functions `local_model` and `local_sim_mean_window` (21 September 2026). The original v7.3 MAT supplies the **11-parameter posterior**, **growth contexts**, **physical constants**, **doses**, **statistical targets** and **observational window**. The MAT file is intentionally NOT included in this public GitHub/Binder repository.

State order matches MATLAB: `[FdeR, Nar, GFP, Cells, Nar_external]`. Four numerical phases of work are: numerical preequilibration (960 min), measured-growth preinduction, extracellular-dose injection, and measured-growth postinduction. The observable is the **mean of GFP at MATLAB 1-based indices 72:108**, plus the relevant `B1` or `B2` **once**. MATLAB `ode23t/ode23s` is approximated using SciPy `solve_ivp(BDF)`; differing solvers may produce small numerical discrepancies.'''),
 nbf.v4.new_code_cell(HEADER),
 nbf.v4.new_markdown_cell(r'''## 1. Load original frozen configuration, if available

Put the original `Data_pFdeR_final_batchaware_VBMC_20260921_214836.mat` under `private_data/` in a **private Codespace or local workstation**, not in public Binder. Without it, this notebook uses a **fabricated synthetic context** to demonstrate that the code runs, but does not claim correspondence with the historical experiment.'''),
 nbf.v4.new_code_cell('''import numpy as np
import matplotlib.pyplot as plt
from biosensor_pipeline.models.final_fder import (
    load_final_model, simulate_final_dose_response, final_curve_features,
    construct_rbs_coordinates, logposterior_final_physical)
from biosensor_pipeline.models.synthetic_final_fixture import make_synthetic_final_model
mat_path=ROOT/'private_data/Data_pFdeR_final_batchaware_VBMC_20260921_214836.mat'
real_archive=mat_path.is_file()
model=load_final_model(mat_path) if real_archive else make_synthetic_final_model()
print('Source:', 'ORIGINAL private 2026-09-21 frozen MAT' if real_archive else 'SYNTHETIC fabricated context (not experimental)')
print('11 parameter names:', list(__import__('biosensor_pipeline.models.final_fder', fromlist=['THETA_NAMES']).THETA_NAMES))
print('Observation window [MATLAB 1-based]:',model.obs_window_matlab)
print('Construct q coordinates: R4R3=',construct_rbs_coordinates(model.theta_median,model.qR4_nominal,'R4R3'),
      'R3R3=',construct_rbs_coordinates(model.theta_median,model.qR4_nominal,'R3R3'))'''),
 nbf.v4.new_markdown_cell(r'''## 2. Simulate contexts and RBS interventions (MATLAB experiment semantics)

- `E1`: R4R3 and R3R3, with common background B1.
- `E2`: R4R3 and R4R4, with common background B2 and `kdNar_shared × cK2`.
- `gamma_TF` changes the effective TF coordinate for R3R3; it is **not** a simple raw qTF=1 in the frozen likelihood.'''),
 nbf.v4.new_code_cell('''fig,axs=plt.subplots(1,2,figsize=(11.5,4),layout='constrained')
summary=[]
for ax,(exp,constructs) in zip(axs,[('E1',('R4R3','R3R3')),('E2',('R4R3','R4R4'))]):
    for construct in constructs:
        doses,y=simulate_final_dose_response(model,experiment=exp,construct=construct)
        stats=final_curve_features(doses,y)
        summary.append((exp,construct,stats['Leakiness'],stats['DynamicRange'],stats['EC50']))
        ax.plot(doses,y,'o-',label=f'{construct}: EC50 fit ≈ {stats["EC50"]:.0f} µM')
    ax.set(title=exp+' (archived context)' if real_archive else exp+' (synthetic context)',
           xlabel='External naringenin (µM)',ylabel='GFP [MEFL/cell-equivalent]')
    ax.grid(alpha=.2);ax.legend()
plt.show()
for row in summary:print(row)'''),
 nbf.v4.new_markdown_cell(r'''## 3. Exact archived likelihood structure: 6D + 3D

The archived MATLAB code uses **correlated Gaussian blocks**, not the older diagonal five-feature likelihood. The physical-coordinate uniform prior contributes a Jacobian in transformed VBMC coordinates. A matching log posterior is a stringent numerical check, but **we do not run a new inference**.'''),
 nbf.v4.new_code_cell('''if real_archive:
    lp=logposterior_final_physical(model)
    print('Python posterior-median logp:',lp)
    print('Frozen MATLAB logp:',model.reference_logpost_median)
    print('Absolute numerical discrepancy:',abs(lp-model.reference_logpost_median))
else:
    print('Skipped logp comparison: synthetic fixture has no frozen experimental Gaussian targets.')'''),
 nbf.v4.new_markdown_cell(r'''## 4. Full private validation against original MATLAB (outside this notebook)

With the original MAT securely available, run:

```bash
python -m scripts.validate_final_fder \\
  --mat private_data/Data_pFdeR_final_batchaware_VBMC_20260921_214836.mat \\
  --out outputs/private_parity --with-seed
```

The script compares **Exp1 R4R3/R3R3** and **Exp2 R4R3/R4R4**, including the likelihood. An independent reference run on 31 October 2026 compared 26 archived summaries (initial seed and posterior median) with a maximum relative feature difference **6.68×10⁻⁶** and a maximum absolute log-posterior discrepancy **4.96×10⁻⁴**. These numerical errors are small, but this does not establish bitwise equivalence or the correctness of a future new posterior fit.''')
])

save('05_historical_2025_candidate.ipynb',[
 nbf.v4.new_markdown_cell(r'''# 05 — Recovered 2025 FdeR candidate: historical de-repression model

This implements the archival MATLAB formula in `model_exp_biosensor.txt` (2025) using the **five-state layout** `[FdeR, intracellular Nar, GFP, cell count, extracellular Nar]`. The later refined model is **not** used in this notebook.

**Critical provenance:** the precise RBS-to-`pR` and RBS-to-`pg` mappings, the exact design-run `all_together.mat`, original posterior `results_1303261822.mat`, and the original forecast output have not been recovered. This is only an executable mathematical translation of a **candidate historical lineage**. It must not be described as the numerical model that originally predicted R3R3.

The 2025 promoter law is:
$$P_{fdeO} = \frac{(k_{dfO}CN)^n(k_{dNar}+I)^n}{(k_{dfO}CN)^n(k_{dNar}+I)^n+(k_{dNar}R)^n}.$$

This law decreases with increasing FdeR at fixed internal inducer; it is a **de-repression** law and is mechanistically distinct from the final ligand-activated FdeR law. Simply strengthening FdeR translation does **not** yield the observed leftward EC50 shift under the same conditions. Do not silently change the promoter equation to force a narrative.'''),
 nbf.v4.new_code_cell(HEADER),
 nbf.v4.new_markdown_cell(r'''## 1. Integrate the recovered 2025 ODE with synthetic conditions only'''),
 nbf.v4.new_code_cell('''import numpy as np
import matplotlib.pyplot as plt
from biosensor_pipeline.models.historical_2025 import (
    HistoricalParameters,historical_promoter,historical_rhs,simulate_historical_one)
p=HistoricalParameters()
print('Archival candidate default parameter object:',p)
print('PfdeO(TF=1000, Nar=50):',historical_promoter(1000,50,p))
print('PfdeO(TF=50, Nar=50):',historical_promoter(50,50,p))'''),
 nbf.v4.new_markdown_cell(r'''## 2. Demonstrate why the two model lineages cannot be confused

The hypothetical weaker TF multiplier below is **not an archival RBS mapping**. It is a controlled mathematical intervention to inspect the direction of the old promoter law.'''),
 nbf.v4.new_code_cell('''doses=np.linspace(0,1400,18)
fig,ax=plt.subplots(figsize=(8.5,4.5))
for factor,label in [(1.0,'Nominal TF production'),(.0049/.1458,'TF translation × 0.0336 (hypothetical)')]:
    GFP=[simulate_historical_one(float(d),tf_factor=factor)[2] for d in doses]
    ax.plot(doses,GFP,'o-',label=label)
ax.set(title='Historical 2025 de-repression law (not original design forecast)',
       xlabel='External naringenin (µM)',ylabel='GFP state (model units)')
ax.legend();ax.grid(alpha=.2);plt.show()'''),
 nbf.v4.new_markdown_cell(r'''## 3. Research interpretation

An earlier model-guided design choice R4R3→R3R3 is author-confirmed; its original quantitative prediction is not recovered. This 2025 code is **not proven** to be the exact revision that generated that prediction. Its de-repression law and the later activation law are inconsistent in TF-direction under otherwise fixed assumptions. Testing mechanistic identifiability is therefore a scientific task, **not** a task of simply translating MATLAB to Python. Use notebook 03 for a deliberately minimal *activation-only* model with the desired qualitative EC50 shift; use notebook 04 for the numerically validated final 11D model.''')
])
