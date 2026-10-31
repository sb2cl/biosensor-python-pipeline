"""Private parity report comparing Python port with original frozen MATLAB .mat.

Never publish the original MAT, generated observations, or this report to a
public repo without explicit data authorization. No fitting or VBMC is run.
"""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import h5py
import numpy as np
from biosensor_pipeline.models.final_fder import (
    load_final_model,simulate_final_dose_response,final_curve_features,
    logposterior_final_physical
)


def run_report(mat_path: Path, out_dir: Path, check_seed: bool=True):
    out_dir.mkdir(parents=True,exist_ok=True)
    m=load_final_model(mat_path)
    rows=[]
    scenarios=[('E1','R4R3','R4R3_E1'),('E1','R3R3','R3R3_E1'),
               ('E2','R4R3','R4R3_E2'),('E2','R4R4','R4R4_E2')]
    with h5py.File(mat_path,'r') as f:
        for name,theta,g in [('median',m.theta_median,f['final_vbmc']['median_diag'])] + ([('seed',m.theta_seed,f['final_vbmc']['seed_diag'])] if check_seed else []):
            feature_cache={}
            for exp,construct,matkey in scenarios:
                d,y=simulate_final_dose_response(m,theta,exp,construct)
                feats=final_curve_features(d,y)
                feature_cache[exp,construct]=feats
                for feat in ['Leakiness','DynamicRange','EC50']:
                    original=float(np.asarray(g['stats'][matkey][feat]).ravel()[0])
                    estimated=float(feats[feat]);difference=estimated-original
                    rows.append(dict(point=name,experiment=exp,construct=construct,feature=feat,
                        matlab=original,python=estimated,absolute_difference=difference,
                        relative_error=abs(difference)/max(abs(original),1.0)))
            lp=logposterior_final_physical(m,theta,features=feature_cache)
            original=m.reference_logpost_median if name=='median' else m.reference_logpost_seed
            rows.append(dict(point=name,experiment='both',construct='ALL',feature='LogPosterior',
                             matlab=original,python=lp,absolute_difference=lp-original,
                             relative_error=abs(lp-original)/max(abs(original),1.0)))
    with (out_dir/'fder_final_matlab_python_parity.csv').open('w',newline='',encoding='utf-8') as h:
        wr=csv.DictWriter(h,fieldnames=list(rows[0]));wr.writeheader();wr.writerows(rows)
    maxabs_log=max(abs(r['absolute_difference']) for r in rows if r['feature']=='LogPosterior')
    maxrelative=max(r['relative_error'] for r in rows if r['feature']!='LogPosterior')
    summary=dict(reference_mat=mat_path.name,source='Frozen original MATLAB 11D VBMC posterior',
        parameters_fitted=False,science_engine='Independent Python ode23t/ode23s -> solve_ivp(BDF) port',
        simulation_comparisons=len(rows), max_absolute_logposterior_difference=maxabs_log,
        max_relative_feature_difference=maxrelative,
        caveat='SciPy BDF and Hill least_squares are not numerical identity with MATLAB ode23t/ode23s and lsqcurvefit',
        report_scope=('median and seed' if check_seed else 'median only')+'; no refit; full 20k posterior not replayed')
    (out_dir/'fder_final_parity_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    print('PARITY SUMMARY',json.dumps(summary,indent=2))
    return summary,rows


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mat',type=Path,required=True,help='Frozen original v7.3 VBMC .mat (PRIVATE)')
    parser.add_argument('--out',type=Path,default=Path('outputs/private_parity'))
    parser.add_argument('--with-seed',action='store_true',help='Also test initial 11-parameter seed')
    args=parser.parse_args()
    run_report(args.mat,args.out,args.with_seed)

if __name__=='__main__':main()
