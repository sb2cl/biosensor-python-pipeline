"""OPTIONAL new fit of pedagogical minimal model to an actual reference curve.

This is NOT an archival reconstruction and NEVER alters the 11D frozen fit.
The public repository includes no real experimental observations.

Input CSV: per-dose already calibrated replicate rows, typically the
BioSensor Feature Studio normalized schema system,series,dose,replicate,signal.
Only three parameters fit: Knar, output scale, beta. TF and promoter abundance
remain fixed to illustrate a limited and potentially non-identifiable model.
"""
from __future__ import annotations
from argparse import ArgumentParser
from dataclasses import replace,asdict
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.optimize import least_squares
from biosensor_pipeline.models.minimal_fder import (
    MinimalParameters,simulate_minimal_dose_response,normalized_half_response)


def fit_reference_dataframe(df:pd.DataFrame,p:MinimalParameters|None=None,
                            tf_rbs='R4',out_rbs='R3', max_nfev=100):
    p=p or MinimalParameters()
    if not {'dose','signal'}.issubset(df):
        raise ValueError('Input must have dose and already-calibrated MEFL/cell signal')
    obs=(df.groupby('dose',as_index=False).signal.median()
            .sort_values('dose'))
    if len(obs)<5 or obs.dose.min()!=0:raise ValueError('At least five doses including 0 are required')
    doses=obs.dose.to_numpy(dtype=float);data=obs.signal.to_numpy(dtype=float)
    if np.any(~np.isfinite(data)) or np.any(data<0):raise ValueError('Invalid signals')
    def unpack(z):
        kn,sc,beta=np.exp(z)
        return replace(p,nar_half_uM=float(kn),gfp_full_scale_mefl_per_cell=float(sc),
                       basal_fraction=float(beta))
    def residual(z):
        _,pred=simulate_minimal_dose_response(doses,tf_rbs,out_rbs,unpack(z))
        return np.log1p(pred)-np.log1p(data)
    lo=np.log([20.,1e3,1e-7]);hi=np.log([3000.,3e6,.15])
    start=np.log([p.nar_half_uM,p.gfp_full_scale_mefl_per_cell,p.basal_fraction])
    result=least_squares(residual,start,bounds=(lo,hi),max_nfev=max_nfev,
                         ftol=1e-7,xtol=1e-7,gtol=1e-7)
    fitted=unpack(result.x)
    _,pred=simulate_minimal_dose_response(doses,tf_rbs,out_rbs,fitted)
    _,other=simulate_minimal_dose_response(doses,'R3',out_rbs,fitted)
    summary={'fit_status':'NEW exploratory minimal-model fit; NOT archival design forecast',
             'optimizer_success':bool(result.success),'optimizer_message':str(result.message),
             'fitted_to':tf_rbs+out_rbs,'held_out_predict': 'R3'+out_rbs,
             'n_doses':len(doses),'residual_log1p_norm':float(np.linalg.norm(result.fun)),
             'parameters':asdict(fitted),
             'EC50_reference_operational':normalized_half_response(doses,pred),
             'EC50_predicted_R3_operational':normalized_half_response(doses,other)}
    curves=pd.DataFrame({'dose_uM':doses,'experimental_median_mefl':data,
                         'fitted_reference_mefl':pred,'predicted_R3_no_refit_mefl':other})
    return summary,curves


def main():
    ap=ArgumentParser(description=__doc__)
    ap.add_argument('--csv',type=Path,required=True)
    ap.add_argument('--system',default='FdeR')
    ap.add_argument('--series',default=None,help='Required if source CSV contains multiple series')
    ap.add_argument('--tf',choices=['R3','R4'],default='R4')
    ap.add_argument('--out',type=Path,default=Path('outputs/private_minimal_fit'))
    args=ap.parse_args()
    df=pd.read_csv(args.csv)
    if 'system' in df:df=df[df.system==args.system]
    if args.series:
        if 'series' not in df:raise ValueError('Requested series but CSV lacks series column')
        df=df[df.series==args.series]
    elif 'series' in df and df.series.nunique()>1:
        raise ValueError('CSV has several series; specify --series exactly')
    if df.empty:raise ValueError('No matching rows')
    summary,curves=fit_reference_dataframe(df,tf_rbs=args.tf)
    args.out.mkdir(parents=True,exist_ok=True)
    (args.out/'minimal_fit_metadata.json').write_text(json.dumps(summary,indent=2)+'\n')
    curves.to_csv(args.out/'minimal_fit_and_prediction.csv',index=False)
    print(json.dumps({k:v for k,v in summary.items() if k!='parameters'},indent=2))

if __name__=='__main__':main()
