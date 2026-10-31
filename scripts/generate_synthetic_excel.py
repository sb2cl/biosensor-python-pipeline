"""Create reproducible fictional, ALREADY CALIBRATED MEFL/cell timecourses."""
from pathlib import Path
import numpy as np
from openpyxl import Workbook


def make_example(path: str | Path) -> Path:
    out=Path(path)
    out.parent.mkdir(parents=True,exist_ok=True)
    w=Workbook();w.remove(w.active)
    rng=np.random.default_rng(20261031)
    for sheet_name,doses,times,scale,K,n,basal in [
        ('FdeR_R3R3',[0,30,75,100,145,250,400,700,1400],np.arange(0,630,5),112000.,115.,1.7,1250.),
        ('pLux_B0032',[0,25,50,75,125,200,400],np.arange(0,610,10),210000.,110.,2.2,1100.)
    ]:
        sh=w.create_sheet(sheet_name)
        unit='dose_uM' if sheet_name.startswith('FdeR') else 'dose_nM'
        sh.append(['time_min',unit,'replicate','MEFL_per_cell'])
        for dose in doses:
            activated=(dose**n)/(K**n+dose**n) if dose>0 else 0.
            for r in range(1,4):
                fold=rng.normal(1,0.035)
                for t in times:
                    dynamics=0.2+0.8*(1-np.exp(-max(0,t-80)/115))
                    signal=max(0.,basal+scale*activated*dynamics*fold+rng.normal(0,360))
                    sh.append([int(t),int(dose),f'rep{r:02d}',float(round(signal,3))])
    w.save(out)
    return out

if __name__=='__main__':
    p=make_example('examples/synthetic_MEFL_timecourses.xlsx')
    print('Wrote synthetic-only workbook:',p)
