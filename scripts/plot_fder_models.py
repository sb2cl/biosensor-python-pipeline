"""Make public synthetic comparison plots for model teaching and code checks."""
from pathlib import Path
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from biosensor_pipeline.models.minimal_fder import (MinimalParameters,
    simulate_minimal_dose_response,simulate_minimal_one,normalized_half_response)
from biosensor_pipeline.models.historical_2025 import simulate_historical_one


def make_plots(out_dir):
    out_dir=Path(out_dir);out_dir.mkdir(parents=True,exist_ok=True)
    doses=np.linspace(0,1400,61)
    p=MinimalParameters()
    fig,ax=plt.subplots(figsize=(9.0,5.4))
    features={}
    for label,color,style in [('R4R3','tab:blue','-'),('R3R3','tab:orange','-'),('R4R4','tab:green','--')]:
        d,y=simulate_minimal_dose_response(doses,label[:2],label[2:],p)
        ec50=normalized_half_response(d,y);features[label]=(float(ec50),float(y[0]),float(max(y)-min(y)))
        ax.plot(d,y,label=f'{label}  EC50* ≈ {ec50:.0f} µM',color=color,linestyle=style,lw=2)
        ax.axvline(ec50,ls=':',lw=0.9,color=color,alpha=0.6)
    ax.set(xlabel='Naringenin supplied (µM)',ylabel='GFP output (illustrative MEFL/cell)',
        title='Minimal FdeR activating model · NO posterior fit')
    ax.grid(alpha=0.2);ax.legend(loc='upper left',fontsize=9)
    fig.subplots_adjust(bottom=.22)
    fig.text(.12,.05,'* EC50: operational finite-range midpoint (not the paper Hill estimator).\n'
            'Without TFcapacity, the minimal model cannot reproduce the reduced R3R3 dynamic range.',
            fontsize=9,va='bottom',color='dimgray')
    fig.savefig(out_dir/'minimal_FdeR_RBS_dose_response.png',dpi=180)
    fig.savefig(out_dir/'minimal_FdeR_RBS_dose_response.svg')
    plt.close(fig)
    # Original 2025 candidate: TF weak response vs strong. No RBS map assumed.
    doses_hist=np.linspace(0,1400,19)
    fig,ax=plt.subplots(figsize=(9.0,5.4))
    for fac,label in [(1.0,'Nominal TF production (1.0x)'),(p.q_tf('R4'),'TF production × 0.0336 (hypothetical)')]:
        yy=np.array([simulate_historical_one(float(dose),tf_factor=fac)[2] for dose in doses_hist])
        ax.plot(doses_hist,yy,lw=2,marker='o',markersize=2,label=label)
    ax.set(xlabel='Naringenin supplied (µM)',ylabel='GFP state (legacy model units)',
            title='Candidate 2025 DE-REPRESSION law · not archived design forecast')
    ax.grid(alpha=.25);ax.legend()
    ax.text(.02,.02,'The original RBS-to-pR mapping and experimental growth contexts were not recovered.\n'
            'These curves are purely hypothetical, not the historical predictions.',
            transform=ax.transAxes,fontsize=8,color='dimgray',va='bottom')
    fig.savefig(out_dir/'historical_2025_candidate_derepression.png',dpi=180)
    plt.close(fig)
    return features


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--out',type=Path,default=Path('outputs/demo_models'))
    args=parser.parse_args()
    print(make_plots(args.out))

if __name__=='__main__':main()
