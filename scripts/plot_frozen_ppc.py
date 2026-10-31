"""Create paper-sized, editable SVG/PNG diagnostic plots from FROZEN MATLAB data.

No refit and no new ODE simulation. Roles are read directly from the archived
feature summary table, retaining held-out/PPC distinctions.
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from biosensor_pipeline.frozen_mat import summarize_ppc


def make_figure(ppc_path, experimental_summary_path, output_dir):
    output_dir=Path(output_dir);output_dir.mkdir(parents=True,exist_ok=True)
    experimental=pd.read_csv(experimental_summary_path)
    expected=['Block','Feature','Role','ExperimentalMedian','ExperimentalCI2p5',
              'ExperimentalCI97p5','PosteriorModelMedian','PosteriorModelCI2p5','PosteriorModelCI97p5']
    if not set(expected).issubset(experimental.columns):
        raise ValueError('Unexpected archived FdeR feature CSV schema')
    frozen=summarize_ppc(ppc_path)
    # MATLAB .mat exact posterior quantiles and the archived CSV should be close;
    # print max drift as a QC report (do not silently accept wrong file combos).
    verify=experimental.merge(frozen,on=['Block','Feature'],validate='one_to_one')
    relative=np.abs(verify.Median-verify.PosteriorModelMedian) / np.maximum(np.abs(verify.PosteriorModelMedian),1e-10)
    if relative.max()>0.04:
        raise ValueError(f'Frozen PPC and CSV disagree at >4% (max={relative.max():.3g}); verify files')
    # Plot ORIGINAL MATLAB-archived quantiles; Python interpolation convention
    # can differ slightly for upper/lower percentiles.
    verify[['Block','Feature','Role','Median','CI2p5','CI97p5']].to_csv(output_dir/'frozen_PPC_python_quantiles.csv',index=False)
    layout=['R4R3_Exp1','R3R3_Exp1','R4R3_Exp2','R4R4_Exp2']
    fig,axes=plt.subplots(1,2,figsize=(11,4.4),layout='constrained')
    for ax,feature in zip(axes, ['DynamicRange','EC50']):
        tab=verify[verify.Feature==feature].set_index('Block').loc[layout]
        loc=np.arange(len(layout))
        for values,low,high,shift,label,marker in [
            ('ExperimentalMedian','ExperimentalCI2p5','ExperimentalCI97p5',-0.13,'Experiment','o'),
            ('PosteriorModelMedian','PosteriorModelCI2p5','PosteriorModelCI97p5',0.13,'Frozen model PPC','s')]:
            med=tab[values].to_numpy(float)
            ymin=med-tab[low].to_numpy(float);ymax=tab[high].to_numpy(float)-med
            ax.errorbar(loc+shift,med,yerr=[ymin,ymax],marker=marker,linestyle='none',
                        capsize=3,markersize=5,label=label)
        ax.set_xticks(loc,layout,rotation=23,ha='right')
        ax.set_title('Dynamic range' if feature=='DynamicRange' else r'EC50 ($\mu$M)')
        if feature=='DynamicRange':ax.set_yscale('log');ax.set_ylabel('MEFL / cell (log scale)')
        else:ax.set_ylabel(r'Concentration ($\mu$M)')
        ax.grid(axis='y',alpha=.25)
        ax.legend(fontsize=8)
        if feature=='DynamicRange':
            # last point is held out in Exp2; annotation from historical data.
            ax.text(0.02,0.04,'R4R4 Exp2 DR: held-out PPC',transform=ax.transAxes,fontsize=8)
    fig.suptitle('FdeR: archived experimental vs posterior-predictive features',fontsize=12)
    svg=output_dir/'FdeR_frozen_PPC_features.svg'
    png=output_dir/'FdeR_frozen_PPC_features.png'
    fig.savefig(svg,bbox_inches='tight')
    fig.savefig(png,dpi=240,bbox_inches='tight')
    plt.close(fig)
    return svg,png,output_dir/'frozen_PPC_python_quantiles.csv'

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--ppc',required=True,help='Archived Data_pFdeR_final_batchaware_PPC_core*.mat')
    p.add_argument('--experimental',required=True,help='Archived final_batchaware_PPC_features*.csv')
    p.add_argument('--out',default='outputs/figures')
    args=p.parse_args()
    for path in make_figure(args.ppc,args.experimental,args.out):print('Saved:',path)
