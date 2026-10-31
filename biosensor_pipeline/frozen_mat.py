"""Read-only *specific* frozen 2026 FdeR MATLAB v7.3 outputs.

This module does not load the MATLAB VBMC variational posterior object `vp`,
which is MATLAB-class-specific. It reads numerical posterior draws and PPC
arrays already archived, WITHOUT re-fitting or simulating.
"""
from pathlib import Path
import h5py
import numpy as np
import pandas as pd

PHYSICAL_NAMES=('etaTF','kG_R3','kdNar_shared','beta','n_nar',
                'KI_TF_rel','gamma_TF','rho0_TF','B1','B2','cK2')
FEATURES=('Leakiness','DynamicRange','EC50','Sensitivity_max','Hill_coefficient')
BLOCKS={'P41':'R4R3_Exp1','P31':'R3R3_Exp1',
        'P42':'R4R3_Exp2','P44':'R4R4_Exp2'}


def read_frozen_posterior(path:str|Path)->pd.DataFrame:
    with h5py.File(path,'r') as h:
        if 'Xs' not in h:raise ValueError('Missing Xs in frozen posterior file')
        x=np.asarray(h['Xs'],dtype=float).T  # MATLAB rows/columns are reversed in HDF5
    if x.ndim!=2 or x.shape[1]!=11 or not np.isfinite(x).all():
        raise ValueError(f'Expected valid N x 11 physical samples, got {x.shape}')
    return pd.DataFrame(x,columns=PHYSICAL_NAMES)


def read_frozen_ppc(path:str|Path)->pd.DataFrame:
    rows=[]
    with h5py.File(path,'r') as h:
        for key,name in BLOCKS.items():
            if key not in h:raise ValueError(f'Missing frozen PPC block: {key}')
            x=np.asarray(h[key],dtype=float).T
            if x.ndim!=2 or x.shape[1]!=5:
                raise ValueError(f'Unexpected {key} array {x.shape}')
            for index in range(len(x)):
                for j,feature in enumerate(FEATURES):
                    rows.append((name,index,feature,float(x[index,j])))
    return pd.DataFrame(rows,columns=['Block','draw','Feature','ModelValue'])


def summarize_ppc(path:str|Path)->pd.DataFrame:
    df=read_frozen_ppc(path)
    q=df.groupby(['Block','Feature']).ModelValue.quantile([.025,.5,.975]).unstack(-1)
    q.columns=['CI2p5','Median','CI97p5']
    return q.reset_index()
