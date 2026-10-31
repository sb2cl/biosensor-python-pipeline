"""Plot private MATLAB-vs-Python results; only local, never bundle by default."""
from pathlib import Path
import argparse
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

def plot_report(csv_path,out_dir):
    df=pd.read_csv(csv_path)
    out_dir=Path(out_dir);out_dir.mkdir(parents=True,exist_ok=True)
    fr=df[(df.feature!='LogPosterior') & (df.point=='median')].copy()
    fig,axes=plt.subplots(1,3,figsize=(11.8,4.2),layout='constrained')
    for ax,key in zip(axes,['Leakiness','DynamicRange','EC50']):
        sub=fr[fr.feature==key]
        xx=np.arange(len(sub))
        labels=[x.replace('R','R',1).replace('_',' ') for x in sub.construct+' '+sub.experiment]
        ax.plot(xx,sub.matlab,'o',label='Frozen MATLAB',ms=6)
        ax.plot(xx,sub.python,'x',label='Python BDF',ms=7)
        ax.set_title(key);ax.set_xticks(xx,labels,rotation=45,ha='right');ax.grid(alpha=.2)
        ax.ticklabel_format(axis='y',style='plain',useOffset=False)
    axes[0].legend()
    fig.suptitle('Frozen 11D FdeR median | MATLAB vs independent Python port')
    fig.savefig(out_dir/'FdeR_final_MATLAB_vs_Python.png',dpi=180)
    fig.savefig(out_dir/'FdeR_final_MATLAB_vs_Python.svg')
    plt.close(fig)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--csv',type=Path,required=True);ap.add_argument('--out',type=Path,default=Path('outputs/private_parity'));a=ap.parse_args();plot_report(a.csv,a.out)
