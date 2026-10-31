"""Literal mathematical transcription of candidate 2025-era model_exp_biosensor.txt.

Recovered CANDIDATE historical lineage, not proven exact version used for RBS design.
This promoter uses DE-REPRESSION, and is distinct from minimal_fder.py's
later activation-only core; it must never be treated as the archived design fit.
"""
from dataclasses import dataclass
import numpy as np
from scipy.integrate import solve_ivp

@dataclass(frozen=True)
class HistoricalParameters:
    CN:float=35.;Vcell_L:float=2e-15;Vext_L:float=200e-6
    D_min_inv:float=2.0; dmR:float=np.log(2)/3; dR:float=np.log(2)/(20*60)
    kR:float=4.85; pR:float=2.016
    dmg:float=np.log(2)/3; dg:float=np.log(2)/(10*60)
    kg:float=1.7014;pg:float=2.5
    dA:float=0.0004;dAe:float=0.0000481
    kdfO:float=1.75;kdNar:float=50.; beta:float=0.0135;n_nar:float=2.
    mu:float=0.0004
    avogadro:float=6.023e23


def historical_promoter(FdeR_molecules,Nar_molecules,p:HistoricalParameters):
    c=(p.kdfO*p.CN)**p.n_nar*(p.kdNar+max(Nar_molecules,0.0))**p.n_nar
    r=(p.kdNar*max(FdeR_molecules,0.0))**p.n_nar
    return c/(c+r)


def historical_rhs(t,x,p:HistoricalParameters):
    """State order [FdeR, Nar_in, GFP, Cells, Nar_external], MATLAB legacy."""
    R,I,G,Cells,Ie=map(float,x)
    mu=p.mu
    dR=p.pR*p.CN*p.kR/(p.dmR+mu)-(p.dR+mu)*max(R,0.0)
    dI=p.D_min_inv*p.Vcell_L/p.Vext_L*max(Ie,0.0)-(p.D_min_inv+p.dA+mu)*max(I,0.0)
    cG=p.pg*p.CN*p.kg/(p.dmg+mu)
    dG=cG*(p.beta+(1-p.beta)*historical_promoter(R,I,p))-(p.dg+mu)*max(G,0.0)
    dC=mu*Cells
    dIe=-p.D_min_inv*Cells*p.Vcell_L/p.Vext_L*max(Ie,0.0)+p.D_min_inv*Cells*max(I,0.0)-p.dAe*max(Ie,0.0)
    return np.array([dR,dI,dG,dC,dIe])


def simulate_historical_one(dose_uM,tf_factor=1.0,out_factor=1.0,p=None,duration_min=960):
    """Illustrative forward run, not a reproduction of the missing historical context.

    Changes pR and pg as hypothetical RBS multipliers; their historical RBS
    mappings were not recovered. The underlying archival model was not refitted.
    """
    from dataclasses import replace
    p=p or HistoricalParameters()
    p=replace(p,pR=p.pR*tf_factor,pg=p.pg*out_factor)
    # Steady TF, basal GFP at zero, inoculum not archived: DEMO only.
    R0=p.pR*p.CN*p.kR/((p.dmR+p.mu)*(p.dR+p.mu))
    C0=1e7
    Ie0=float(dose_uM)*1e-6*p.avogadro*p.Vext_L
    y0=[R0,0.,0.,C0,Ie0]
    sol=solve_ivp(lambda t,x:historical_rhs(t,x,p),(0,duration_min),y0,
                  method='BDF',rtol=1e-7,atol=1e-8,t_eval=[duration_min])
    if not sol.success:raise RuntimeError(sol.message)
    return sol.y[:,-1]
