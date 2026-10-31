"""Minimal FdeR/inducer/GFP activating model for pedagogical RBS design.

NOT the original 2025 archival model; NOT the 11-parameter final VBMC model.
This uses the activation-only mechanistic core of the later model, without
TF compression, TF-dependent basal activation, output-capacity penalties,
batch offsets or experimental growth traces. Values below are illustrative.

States in order: TF [nM], intracellular naringenin [uM], GFP [MEFL/cell].
Time in minutes; external dose [uM] is a constant control input.
"""
from __future__ import annotations
from dataclasses import dataclass, replace
import numpy as np
from scipy.integrate import solve_ivp


@dataclass(frozen=True)
class MinimalParameters:
    tf_R3_steady_nM: float = 600.0
    promoter_half_nM: float = 18.0
    nar_half_uM: float = 600.0
    nar_hill: float = 2.0
    promoter_hill: float = 1.0
    basal_fraction: float = 0.003
    gfp_full_scale_mefl_per_cell: float = 300_000.0
    growth_rate_min_inv: float = 0.0004
    tf_degradation_min_inv: float = np.log(2) / 1200
    gfp_degradation_min_inv: float = np.log(2) / 600
    transport_in_min_inv: float = 0.10
    nar_degradation_min_inv: float = 0.0004
    kappa_R3: float = 0.1458
    kappa_R4: float = 0.0049

    def q_tf(self, rbs: str) -> float:
        if rbs not in {'R3', 'R4'}:
            raise ValueError('TF RBS must be R3 or R4')
        return {'R3': self.kappa_R3, 'R4': self.kappa_R4}[rbs]/self.kappa_R3

    def q_out(self, rbs: str) -> float:
        if rbs not in {'R3','R4'}:
            raise ValueError('Output RBS must be R3 or R4')
        return {'R3': self.kappa_R3, 'R4': self.kappa_R4}[rbs]/self.kappa_R3


def promoter_activity(tf_nM: float, nar_uM: float, p: MinimalParameters) -> float:
    """Activation-only FdeR-Nar-PfdeO law from later mechanistic lineage."""
    nar=max(0.0, float(nar_uM)); tf=max(0.0,float(tf_nM))
    f_nar=nar**p.nar_hill/(p.nar_half_uM**p.nar_hill+nar**p.nar_hill)
    tf_active=tf*f_nar
    h=tf_active**p.promoter_hill/(p.promoter_half_nM**p.promoter_hill+tf_active**p.promoter_hill)
    return p.basal_fraction+(1-p.basal_fraction)*h


def minimal_rhs(t, state, external_nar_uM: float, tf_q:float, out_q:float, p:MinimalParameters):
    """d[TF, Nar_in, GFP]/dt. Constant exterior is not an integrated state."""
    tf, nar, gfp=map(float,state)
    d_tf=p.tf_degradation_min_inv+p.growth_rate_min_inv
    d_g=p.gfp_degradation_min_inv+p.growth_rate_min_inv
    # Reference R3 steady abundance is used only to parameterize TF synthesis.
    tf_synthesis=p.tf_R3_steady_nM*d_tf*tf_q
    d_tf_state=tf_synthesis-d_tf*tf
    d_nar=p.transport_in_min_inv*(external_nar_uM-max(nar,0.0))-(p.nar_degradation_min_inv+p.growth_rate_min_inv)*max(nar,0.0)
    gfp_production=p.gfp_full_scale_mefl_per_cell*d_g*out_q*promoter_activity(tf,nar,p)
    d_gfp=gfp_production-d_g*max(gfp,0.0)
    return np.array([d_tf_state,d_nar,d_gfp],dtype=float)


def simulate_minimal_one(dose_uM:float,tf_rbs='R3',out_rbs='R3',p:MinimalParameters|None=None,
                         duration_min:float=960., times=None):
    p=p or MinimalParameters()
    if dose_uM<0: raise ValueError('Nonnegative inducer concentration required')
    if times is None: times=np.linspace(0,duration_min,161)
    times=np.asarray(times,dtype=float)
    tf_q=p.q_tf(tf_rbs);out_q=p.q_out(out_rbs)
    # Preinduced TF equilibrium; GFP basal initially at its steady baseline.
    y0=[p.tf_R3_steady_nM*tf_q, 0.0,
        p.gfp_full_scale_mefl_per_cell*out_q*p.basal_fraction]
    sol=solve_ivp(minimal_rhs, (float(times[0]),float(times[-1])),y0,
        args=(float(dose_uM),tf_q,out_q,p),t_eval=times,method='BDF',rtol=1e-8,atol=1e-10)
    if not sol.success: raise RuntimeError(sol.message)
    return sol


def simulate_minimal_dose_response(doses_uM,tf_rbs='R3',out_rbs='R3',p:MinimalParameters|None=None,
                                  duration_min:float=960.):
    """Return (dose_uM, end GFP MEFL/cell); no Hill fitting or experimental calibration."""
    ds=np.asarray(doses_uM,dtype=float).ravel()
    if len(ds)<2 or np.any(ds<0): raise ValueError('Provide at least two nonnegative doses')
    ys=[]
    for d in ds:
        sol=simulate_minimal_one(float(d),tf_rbs,out_rbs,p,duration_min,
                                 times=[0,float(duration_min)])
        ys.append(float(sol.y[2,-1]))
    return ds,np.array(ys)


def normalized_half_response(doses_uM, mefl_per_cell):
    """Observed finite-grid interpolation, explicitly NOT the paper Hill EC50."""
    x=np.asarray(doses_uM,dtype=float);y=np.asarray(mefl_per_cell,dtype=float)
    if np.any(np.diff(x)<=0): raise ValueError('Sort strictly increasing doses')
    y=np.maximum.accumulate(y)
    midpoint=.5*(y[0]+y[-1])
    if y[-1]==y[0]:return np.nan
    # Works with nondecreasing plateaus (np.interp uses last plateau tie if any).
    keep=np.r_[True,np.diff(y)>1e-10]
    return float(np.interp(midpoint,y[keep],x[keep]))
