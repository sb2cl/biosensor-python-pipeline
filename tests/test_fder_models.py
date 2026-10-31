"""Public synthetic tests for three clearly separated model lineages."""
from pathlib import Path
import os
import numpy as np
import pytest

from biosensor_pipeline.models.minimal_fder import (
    MinimalParameters, promoter_activity,simulate_minimal_dose_response,
    normalized_half_response,simulate_minimal_one
)
from biosensor_pipeline.models.historical_2025 import (
    HistoricalParameters,historical_promoter,historical_rhs,simulate_historical_one
)
from biosensor_pipeline.models.final_fder import (
    final_rhs,construct_rbs_coordinates,final_curve_features,
    load_final_model,simulate_final_dose_response,logposterior_final_physical
)


def test_minimal_response_monotone_both_inputs():
    p=MinimalParameters()
    assert promoter_activity(300,100,p)>promoter_activity(50,100,p)
    assert promoter_activity(300,800,p)>promoter_activity(300,100,p)
    assert promoter_activity(300,0,p)==p.basal_fraction


def test_minimal_tf_rbs_moves_ec50_lower_for_stronger_tf():
    p=MinimalParameters();doses=np.linspace(0,1400,31)
    d,a=simulate_minimal_dose_response(doses,'R3','R3',p)
    _,b=simulate_minimal_dose_response(doses,'R4','R3',p)
    e_strong=normalized_half_response(d,a)
    e_weak=normalized_half_response(d,b)
    assert 30<e_strong<200
    assert 250<e_weak<700
    assert e_strong<e_weak


def test_minimal_reporter_rbs_scales_gfp_not_ec50():
    p=MinimalParameters();doses=np.linspace(0,1400,20)
    _,strong=simulate_minimal_dose_response(doses,'R4','R3',p)
    _,weak=simulate_minimal_dose_response(doses,'R4','R4',p)
    np.testing.assert_allclose(weak,strong*p.q_out('R4'),rtol=3e-5)
    assert np.isclose(normalized_half_response(doses,weak),
                      normalized_half_response(doses,strong),atol=0.1)


def test_minimal_timecourse_is_nonnegative_and_finite():
    a=simulate_minimal_one(100,tf_rbs='R3')
    assert a.y.shape==(3,len(a.t))
    assert np.all(np.isfinite(a.y))
    assert np.all(a.y>=0)


def test_historical_promoter_is_derepression_not_final_activator():
    p=HistoricalParameters()
    assert historical_promoter(1000,50,p)<historical_promoter(50,50,p)
    assert historical_promoter(1000,500,p)>historical_promoter(1000,50,p)


def test_historical_rhs_is_five_states():
    dx=historical_rhs(0,np.array([100,10,100,1e7,1e11]),HistoricalParameters())
    assert dx.shape==(5,)
    assert np.all(np.isfinite(dx))


def test_historical_demo_rbs_weakening_increases_output_at_low_dose():
    # This demonstrates the historical desrepression law acts in a different
    # direction from the later ligand-activated TF minimal model.
    ystr=simulate_historical_one(100,tf_factor=1)
    yweak=simulate_historical_one(100,tf_factor=.03360768)
    assert yweak[2]>ystr[2]


def test_final_qtf_compression_from_archived_matlab():
    q=.0049/.1458
    x=np.zeros(11);x[6]=.864
    assert construct_rbs_coordinates(x,q,'R4R3')==(q,1)
    assert construct_rbs_coordinates(x,q,'R4R4')==(q,q)
    q_eff=construct_rbs_coordinates(x,q,'R3R3')[0]
    assert q<q_eff<1


def test_final_rhs_finite_and_five_states():
    p=dict(kdfO_gauge=10000.,CN=35,omega_J23106=.0614,dmR=np.log(2)/3,
            dR=np.log(2)/1200,dA=.0004,dAe=.0000481,D=2,
            Vcell=2e-15,dmg=np.log(2)/3,dg=np.log(2)/600,
            molec_to_uM=.00083026953,molec_to_nM=.83026953,n_fdeO=1.,mI_TF=1.)
    theta=np.array([.0674,13.42,1232.,.00049,1.63,22.3,.864,.0001,428.,894.,.606])
    v=final_rhs(5,np.array([2e6,1000,20000,1e7,1e12]),
                np.array([0,10]),np.array([.0004,.0004]),p,theta,.6,1,200e-6,1232.)
    assert v.shape==(5,) and np.all(np.isfinite(v))


def test_no_private_data_committed_here():
    base=Path(__file__).resolve().parents[1]
    assert not list(base.rglob('*.mat'))

@pytest.mark.skipif('FDER_FROZEN_MAT' not in os.environ,reason='Private MATLAB source not on public CI')
def test_optional_frozen_matlab_parity():
    m=load_final_model(os.environ['FDER_FROZEN_MAT'])
    for exp,construct in [('E1','R4R3'),('E1','R3R3'),('E2','R4R3'),('E2','R4R4')]:
        d,y=simulate_final_dose_response(m,experiment=exp,construct=construct)
        feats=final_curve_features(d,y)
        assert 0<feats['EC50']<max(d)
    assert abs(logposterior_final_physical(m)-m.reference_logpost_median)<0.05


def test_synthetic_final_no_private_mat_needed():
    from biosensor_pipeline.models.synthetic_final_fixture import make_synthetic_final_model
    m=make_synthetic_final_model()
    d,y=simulate_final_dose_response(m,experiment='E1',construct='R4R3')
    assert len(d)==7 and np.all(np.isfinite(y))
    assert y[-1]>y[0]>0


def test_new_minimal_calibration_does_not_modify_full_posterior():
    from dataclasses import replace
    import pandas as pd
    from scripts.calibrate_minimal_fder import fit_reference_dataframe
    p=replace(MinimalParameters(),nar_half_uM=500.,gfp_full_scale_mefl_per_cell=220_000.,basal_fraction=.005)
    doses=np.linspace(0,1400,8)
    d,y=simulate_minimal_dose_response(doses,'R4','R3',p)
    frame=pd.DataFrame({'dose':d,'signal':y})
    result, curves=fit_reference_dataframe(frame,max_nfev=12)
    assert result['n_doses']==8
    assert set(curves.columns)=={'dose_uM','experimental_median_mefl','fitted_reference_mefl','predicted_R3_no_refit_mefl'}
    assert result['EC50_predicted_R3_operational']<result['EC50_reference_operational']
    assert result['fit_status'].startswith('NEW exploratory')
