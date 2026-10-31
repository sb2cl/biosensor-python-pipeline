"""Preliminary Python forward port of the frozen 11D final FdeR MATLAB model.

DIRECTLY TRANSLATES the mathematical RHS and staged protocol of
  logPosterior_primary_batchaware_meanwindow.m (2026-09-21)
  local_model() and local_sim_mean_window().

This is a DIFFERENT code path than the pedagogical minimal model or 2025
historical candidate. No fitting, VBMC, calibration, or posterior regeneration.
Requires a locally supplied original v7.3 MAT file not stored in the repo.

Data retained as in MATLAB: x = [FdeR_molecules, Nar_molecules, GFP_MEFLeq,
Cells_count, Nar_external_molecules]. Observation = mean at 1-based indices
[72:108] of merged Exp1+Exp2 measurement time grid, THEN add B_e once.
"""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import Optional
import numpy as np
import h5py
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares
from scipy.linalg import solve_triangular

STATE_NAMES=('FdeR','Nar_in','GFP','Cells','Nar_external')
THETA_NAMES=('etaTF','kG_R3','kdNar_shared','beta','n_nar','KI_TF_rel',
             'gamma_TF','rho0_TF','B1','B2','cK2')

@dataclass
class GrowthContext:
    label: str
    t_pre: np.ndarray
    mu_pre: np.ndarray
    t1: np.ndarray
    mu1: np.ndarray
    t2: np.ndarray
    mu2: np.ndarray
    initial: np.ndarray
    time_full: np.ndarray
    doses_uM: np.ndarray
    mu_null: float
    parti_initial: float
    parti_induction: float
    gfp_initial: float
    gfp_induction: float
    GFPindex: int  # zero-based Python indices
    ODindex: int
    Nare_index: int
    rtol: float
    atol: float

@dataclass
class FinalModelData:
    parameters: dict[str,float]
    contexts: dict[str,GrowthContext]
    theta_median: np.ndarray
    theta_seed: np.ndarray
    hard_low: np.ndarray
    hard_high: np.ndarray
    qR4_nominal: float
    obs_window_matlab: tuple[int,int]
    targets:dict
    reference_logpost_median: float
    reference_logpost_seed: float


def _v(gr,key):return np.asarray(gr[key][:],dtype=float).ravel()
def _s(gr,key):return float(np.asarray(gr[key][:],dtype=float).ravel()[0])

def _read_context(gr,label):
    opts=gr['options']
    return GrowthContext(label=label, t_pre=_v(gr,'Texp0'),mu_pre=_v(gr,'Gexp0'),
        t1=_v(gr,'Texp1'),mu1=_v(gr,'Gexp1'),t2=_v(gr,'Texp2'),mu2=_v(gr,'Gexp2'),
        initial=_v(gr,'Initial'),time_full=_v(gr,'timeEXP'),
        doses_uM=_v(gr,'induction'),mu_null=_s(gr,'mu_null'),
        parti_initial=_s(gr,'parti_initial'),parti_induction=_s(gr,'parti_induction'),
        gfp_initial=_s(gr,'gfp_initial'),gfp_induction=_s(gr,'gfp_induction'),
        GFPindex=int(_s(gr,'GFPindex'))-1,ODindex=int(_s(gr,'ODindex'))-1,
        Nare_index=int(_s(gr,'Nare_index'))-1,rtol=_s(opts,'RelTol'),atol=_s(opts,'AbsTol'))


def load_final_model(mat_path) -> FinalModelData:
    """Load cfg/p0/g1/g2/target and posterior from the UNMODIFIED v7.3 .mat.

    The file stays local: do not copy any .mat into a public GitHub repo.
    """
    with h5py.File(mat_path,'r') as f:
        cfg=f['cfg'];p0=cfg['p0']
        p={k:_s(p0,k) for k,ds in p0.items() if isinstance(ds,h5py.Dataset) and np.issubdtype(ds.dtype,np.number) and ds.size==1}
        contexts={'E1':_read_context(cfg['g1'],'E1'),'E2':_read_context(cfg['g2'],'E2')}
        targets={}
        for key in ('exp1','exp2'):
            g=cfg['target'][key]
            targets[key]={'mu':_v(g,'mu'),'L':g['L'][:].T.copy(),'logdet':_s(g,'logdet')}
        return FinalModelData(parameters=p,contexts=contexts,
             theta_median=_v(f['final_vbmc'],'posterior_median'),
             theta_seed=_v(cfg,'seed_physical'),
             hard_low=_v(cfg,'physical_LB'),hard_high=_v(cfg,'physical_UB'),
             qR4_nominal=_s(cfg,'qR4_nominal'),
             obs_window_matlab=tuple(map(int,_v(cfg,'obsWindow'))),targets=targets,
             reference_logpost_median=_s(f['final_vbmc'],'logpost_median'),
             reference_logpost_seed=_s(f['final_vbmc'],'logpost_seed'))


def final_rhs(t,x,mu_times,mu_values,p,theta,tf_q,out_q,ext_volume_l,kdNar):
    """Exact algebraic structure of local_model(), one cell compartment.

    Switching scipy/BDF for MATLAB's ode23t/ode23s means numerical solutions
    are NOT guaranteed bitwise equal to MATLAB; validate experimentally.
    """
    R,N,G,C,E=map(float,x)
    eta,kG,_,beta,nNar,KI,gamma,rho0,B1,B2,cK2=map(float,theta)
    mu=float(np.interp(t,mu_times,mu_values))
    Rpos,Npos,Gpos,Epos=max(R,0),max(N,0),max(G,0),max(E,0)
    kTL_R3=eta*p['kdfO_gauge']
    Rprod=(kTL_R3*tf_q*p['CN']*p['omega_J23106']/(p['dmR']+mu))
    dR=Rprod-(p['dR']+mu)*Rpos
    dN=p['D']*p['Vcell']/ext_volume_l*E - (p['D']+p['dA']+mu)*Npos
    N_uM=Npos*p['molec_to_uM']; R_nM=Rpos*p['molec_to_nM']
    f_ind=N_uM**nNar/(kdNar**nNar+N_uM**nNar+np.finfo(float).eps)
    rho=min(max(rho0,0.0),0.999999)
    f_nar=rho+(1-rho)*f_ind
    R_active_nM=R_nM*f_nar
    n_fdeO=p['n_fdeO']
    Hact=R_active_nM**n_fdeO/(p['kdfO_gauge']**n_fdeO+R_active_nM**n_fdeO+np.finfo(float).eps)
    Rrel=R_nM/p['kdfO_gauge']
    capacity=1/(1+(Rrel/KI)**p['mI_TF'])
    PfdeO=beta+(1-beta)*Hact*capacity
    Gprod=(kG*out_q*p['CN']/(p['dmg']+mu))*PfdeO
    dG=Gprod-(p['dg']+mu)*Gpos
    dC=mu*C
    dE=-p['D']*C*p['Vcell']/ext_volume_l*Epos + p['D']*C*Npos-p['dAe']*Epos
    return np.array([dR,dN,dG,dC,dE],dtype=float)


def _stage(tspan,y0,mu_t,mu_val,p,theta,tf_q,out_q,volume,kdNar,method,rtol,atol, t_eval=None):
    fun=lambda t,x:final_rhs(t,x,mu_t,mu_val,p,theta,tf_q,out_q,volume,kdNar)
    sol=solve_ivp(fun,(float(tspan[0]),float(tspan[-1])),y0,
        method=method,rtol=rtol,atol=atol,t_eval=t_eval)
    if not sol.success:
        raise RuntimeError(f'{method} integration failed: {sol.message}')
    return sol


def construct_rbs_coordinates(theta: np.ndarray, qR4:float, name:str):
    """Exact exp-specific qTF_eff compression in 2026 likelihood.

    Unlike physical design coordinate qTF=1 for R3, the final LIKELIHOOD
    uses qR4**(1-gamma_TF) for R3R3 Exp1.
    """
    gamma=float(theta[6])
    if name=='R4R3':return qR4,1.0
    if name=='R3R3':return qR4**(1-gamma),1.0
    if name=='R4R4':return qR4,qR4
    raise ValueError(f'Construct {name!r} not one of frozen 2026 model cases')


def simulate_final_dose_response(model:FinalModelData, theta=None,experiment='E1',
                                 construct='R4R3',doses_uM=None, method='BDF',
                                 rtol=1e-8,atol=1e-9,window_only=True):
    """Reproduce MATLAB observation protocol. Return (dose, mean-GFP+B_e).

    `window_only` accelerates integration by stopping at observation index 108
    rather than the full time course; the observed mean should be unchanged.
    The rtol/atol args are SciPy tolerances, not claims of ode23t identity.
    """
    theta=model.theta_median if theta is None else np.asarray(theta,dtype=float)
    assert len(theta)==11
    ctx=model.contexts[experiment];p=model.parameters
    doses=ctx.doses_uM if doses_uM is None else np.asarray(doses_uM,dtype=float)
    tf_q,out_q=construct_rbs_coordinates(theta,model.qR4_nominal,construct)
    kdNar=float(theta[2])*(float(theta[10]) if experiment=='E2' else 1.0)
    B=float(theta[8] if experiment=='E1' else theta[9])
    # MATLAB does not recompute p.nar_to_molec after switching p.Vext
    y0=ctx.initial.copy();y0[ctx.ODindex]=1.0
    pre=_stage(ctx.t_pre,y0,ctx.t_pre,ctx.mu_pre,p,theta,tf_q,out_q,p['Vext_0'],kdNar,
               method,rtol,atol)
    y1=pre.y[:,-1].copy()
    y1[ctx.ODindex]=ctx.parti_initial
    y1[ctx.GFPindex]=ctx.gfp_initial
    first=_stage(ctx.t1,y1,ctx.t1,ctx.mu1,p,theta,tf_q,out_q,p['Vext_plate'],kdNar,
                 method,rtol,atol)
    y2=first.y[:,-1].copy()
    y2[ctx.ODindex]=ctx.parti_induction
    y2[ctx.GFPindex]=ctx.gfp_induction
    a,b=model.obs_window_matlab
    # One-based MATLAB [72:108] -> zero-based [71:108].
    obs_times=ctx.time_full[a-1:b]
    if min(obs_times)<ctx.t2[0]:raise ValueError('Observation sample lies before induction')
    ys=[]
    for dose in doses:
        y2j=y2.copy();y2j[ctx.Nare_index]=float(dose)*p['nar_to_molec']
        if window_only:
            outtimes=obs_times; tspan=(ctx.t2[0],float(obs_times[-1]))
        else:
            outtimes=ctx.t2; tspan=ctx.t2
        sim=_stage(tspan,y2j,ctx.t2,ctx.mu2,p,theta,tf_q,out_q,p['Vext_plate'],kdNar,
                   method,rtol,atol,t_eval=outtimes)
        if window_only:
            gvalues=sim.y[ctx.GFPindex]
        else:
            # measured times are a subset of t2 except MATLAB stored 60,64.433, etc
            gvalues=np.interp(obs_times,sim.t,sim.y[ctx.GFPindex])
        ys.append(float(np.mean(gvalues)+B))
    return np.asarray(doses,dtype=float),np.asarray(ys,dtype=float)


def hill_model(pars,x):
    b,n,k,K=pars
    x=np.asarray(x,dtype=float)
    num=np.power(x,n)
    return K*(b+(1-b)*num/(np.power(k,n)+num+np.finfo(float).eps))


def final_curve_features(doses, values):
    """SciPy analog of biosensorFeaturesDeterministic.m; not identical solver.

    Hill fit uses SciPy least_squares and may differ from MATLAB lsqcurvefit.
    """
    d=np.asarray(doses,dtype=float);v=np.asarray(values,dtype=float)
    idx=np.argsort(d);d=d[idx];v=v[idx]
    p0=np.array([0.1,2.0,float(np.median(d)),float(np.max(v))])
    result=least_squares(lambda pars:hill_model(pars,d)-v, p0,
                         bounds=([0,0,0,0],[1,np.inf,np.inf,np.inf]),
                         x_scale='jac',ftol=1e-11,xtol=1e-11,gtol=1e-11,max_nfev=1200)
    if not result.success and not np.isfinite(result.cost): raise ValueError('Failed Hill fit')
    xf=np.linspace(d.min(),d.max(),1000)
    yf=hill_model(result.x,xf)
    ec50=float(xf[np.argmin(abs(yf-(yf.min()+yf.max())/2))])
    return dict(Leakiness=float(v[0]),DynamicRange=float(v.max()-v.min()),
                EC50=ec50,Hill_coefficient=float(result.x[1]),hill_fit=result.x)


def logposterior_final_physical(model:FinalModelData,theta=None,method='BDF',rtol=1e-8,atol=1e-9,features=None):
    """Preliminary full 6D+3D log posterior port. Refit prohibited.

    Same normal blocks and physical-coordinate Jacobian. Fitted EC50 can
    differ from lsqcurvefit; compare to frozen diag before using for inference.
    """
    theta=model.theta_median if theta is None else np.asarray(theta,dtype=float)
    if np.any(theta<model.hard_low) or np.any(theta>model.hard_high):return -1e12
    if features is None:
        features={}
        for exp,construct in [('E1','R4R3'),('E1','R3R3'),('E2','R4R3'),('E2','R4R4')]:
            x,y=simulate_final_dose_response(model,theta,exp,construct,method=method,rtol=rtol,atol=atol)
            features[exp,construct]=final_curve_features(x,y)
    s41=features['E1','R4R3'];s31=features['E1','R3R3']
    s42=features['E2','R4R3'];s44=features['E2','R4R4']
    ld1=s31['Leakiness']-s41['Leakiness'];ld2=s42['Leakiness']-s44['Leakiness']
    if min(ld1,ld2,s41['Leakiness'],s41['DynamicRange'],s41['EC50'],
           s31['DynamicRange'],s31['EC50'],s42['Leakiness'],s42['EC50']) <=0:return -1e12
    z1=np.log([s41['Leakiness'],s41['DynamicRange'],s41['EC50'],
             ld1,s31['DynamicRange']/s41['DynamicRange'],s31['EC50']/s41['EC50']])
    z2=np.log([s42['Leakiness'],s42['EC50'],ld2])
    lls=[]
    for key,vec in [('exp1',z1),('exp2',z2)]:
        t=model.targets[key];r=vec-t['mu']
        q=solve_triangular(t['L'],r,lower=True)
        lls.append(float(-.5*(q@q+t['logdet']+len(vec)*np.log(2*np.pi))))
    eta,kG,kd,beta,nar,KI,gamma,rho,B1,B2,cK2=theta
    logjac=float(np.log(eta)+np.log(kG)+np.log(kd)+np.log(beta)+np.log(nar)+
                 np.log(KI)+np.log(rho)+np.log1p(B1)+np.log1p(B2)+np.log(cK2))
    return sum(lls)+logjac
