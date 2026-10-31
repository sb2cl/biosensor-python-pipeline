"""Public synthetic final-model demonstration context, NOT experimental data."""
import numpy as np
from .final_fder import GrowthContext, FinalModelData


def make_synthetic_final_model():
    tpre=np.arange(0,961,10,dtype=float)
    t1=np.arange(0,61,5,dtype=float)
    t2=np.arange(60,601,5,dtype=float)
    t=np.r_[t1[:-1],t2]
    av=6.023e23;cell_L=2e-15;plate_L=200e-6
    p=dict(kdfO_gauge=1e4,CN=35,omega_J23106=.0614,dmR=np.log(2)/3,
        dR=np.log(2)/1200,dA=.0004,dAe=.0000481,D=2.,
        Vcell=cell_L,Vext_0=3e-3,Vext_plate=plate_L,
        dmg=np.log(2)/3,dg=np.log(2)/600,
        molec_to_uM=1e6/(cell_L*av),molec_to_nM=1e9/(cell_L*av),
        nar_to_molec=plate_L*av*1e-6,n_fdeO=1.,mI_TF=1.)
    cs={}
    for name in ['E1','E2']:
        mu=0.0005 if name=='E1' else 0.0004
        ctx=GrowthContext(label=name,t_pre=tpre,mu_pre=np.full(len(tpre),mu),
                t1=t1,mu1=np.full(len(t1),mu),t2=t2,mu2=np.full(len(t2),mu),
                initial=np.zeros(5),time_full=t.copy(),
                doses_uM=np.array([0.,30.,100.,200.,400.,700.,1400.]),
                mu_null=mu,parti_initial=7.5e6,parti_induction=2e7,
                gfp_initial=1500.,gfp_induction=800.,GFPindex=2,ODindex=3,
                Nare_index=4,rtol=1e-6,atol=1e-7)
        cs[name]=ctx
    theta=np.array([.055,12.,1150.,.0005,1.7,20.,.85,.0001,400.,850.,.65])
    return FinalModelData(parameters=p,contexts=cs,theta_median=theta.copy(),theta_seed=theta.copy(),
             hard_low=np.full(11,-np.inf),hard_high=np.full(11,np.inf),
             qR4_nominal=.0049/.1458,obs_window_matlab=(72,108),targets={},
             reference_logpost_median=float('nan'),reference_logpost_seed=float('nan'))
