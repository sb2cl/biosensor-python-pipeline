from pathlib import Path
import numpy as np
import h5py
import pytest
from biosensor_pipeline.frozen_mat import read_frozen_posterior,read_frozen_ppc,summarize_ppc,FEATURES,BLOCKS


def test_archival_numeric_arrays_are_transposed_from_matlab_layout(tmp_path):
    f=tmp_path/'x.mat'
    with h5py.File(f,'w') as h:
        h['Xs']=np.arange(11*3,dtype=float).reshape(11,3)
        for k in BLOCKS:h[k]=np.arange(5*3,dtype=float).reshape(5,3)
    post=read_frozen_posterior(f)
    assert post.shape==(3,11)
    assert post.etaTF.tolist()==[0.,1.,2.]
    ppc=read_frozen_ppc(f)
    assert len(ppc)==4*3*5
    assert ppc[ppc.Block=='R4R3_Exp1'].loc[ppc.Feature=='DynamicRange','ModelValue'].tolist()==[3.,4.,5.]
    summary=summarize_ppc(f)
    assert len(summary)==20
    assert (summary.CI2p5<=summary.Median).all()


def test_fail_on_wrong_11d_posterior_layout(tmp_path):
    f=tmp_path/'wrong.mat'
    with h5py.File(f,'w') as h:h['Xs']=np.zeros((10,30))
    with pytest.raises(ValueError,match='11'):
        read_frozen_posterior(f)
