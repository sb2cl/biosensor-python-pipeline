"""Schema-driven, non-destructive experimental Excel -> normalized timecourses.

Research scope: calibrated MEFL/cell inputs only. No plate-reader calibration,
no automatic strain/dose assignment, no automated rejection of suspicious doses.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from openpyxl import load_workbook


class InputMappingError(ValueError):
    """Input does not meet the declared, explicit mapping contract."""


def sha256_file(path: str | Path) -> str:
    h=hashlib.sha256()
    with Path(path).open('rb') as fh:
        for block in iter(lambda: fh.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def inspect_excel(path: str | Path, *, n_preview: int = 5) -> list[dict[str, Any]]:
    """Read-only sheet inventory. Does not assume vendor-specific layouts."""
    wb=load_workbook(path, read_only=True, data_only=True)
    output=[]
    try:
        for sheet in wb.worksheets:
            sheet.reset_dimensions()  # avoids silently truncated workbook metadata
            it=sheet.iter_rows(values_only=True)
            rows=[]
            for ix,row in enumerate(it):
                if ix >= n_preview:
                    break
                rows.append([None if v is None else str(v)[:120] for v in row[:24]])
            output.append({'sheet':sheet.title,'preview_rows':rows,
                           'note':'First rows only. Original file is read-only.'})
    finally:
        wb.close()
    return output


def _require_columns(frame: pd.DataFrame, mapping: dict[str,str], fields: set[str]) -> None:
    missing=fields-set(mapping)
    if missing:
        raise InputMappingError(f'Explicit mapping missing fields: {sorted(missing)}')
    missing=[mapping[k] for k in fields if mapping[k] not in frame.columns]
    if missing:
        raise InputMappingError(f'Excel worksheet missing required headers: {missing}')


def _prepare_long(frame: pd.DataFrame, ds: dict[str, Any]) -> pd.DataFrame:
    cols=ds['columns']
    _require_columns(frame,cols,{'time','dose','replicate','signal'})
    out=frame[[cols[k] for k in ('time','dose','replicate','signal')]].copy()
    out.columns=['time_min','dose','replicate','signal']
    return out


def _prepare_wide(frame: pd.DataFrame, ds: dict[str, Any]) -> pd.DataFrame:
    """Wide Excel: time column + explicitly mapped well columns, no guesswork."""
    tcol=ds['time_column']
    if tcol not in frame.columns:
        raise InputMappingError(f'Missing time column: {tcol}')
    rows=[]
    for col in ds['wells']:
        name=col['column']
        if name not in frame.columns:
            raise InputMappingError(f'Missing mapped fluorescence well: {name}')
        rows.append(pd.DataFrame({'time_min':frame[tcol], 'dose':col['dose'],
                                  'replicate':str(col['replicate']),
                                  'signal':frame[name]}))
    if not rows:
        raise InputMappingError('Wide-timecourse configuration has no wells')
    return pd.concat(rows,ignore_index=True)


def _window_indices(g:pd.DataFrame, window:dict[str,Any]) -> pd.Index:
    method=window['kind']
    if method=='sample_index_1based':
        start=int(window['start']); stop=int(window['stop'])
        if start < 1 or stop < start:
            raise InputMappingError('Invalid MATLAB 1-based observation indices')
        if stop > len(g):
            raise InputMappingError(f'Window [{start}:{stop}] exceeds {len(g)} observed rows')
        # MATLAB indexing [72:108] is inclusive, not [72,108] minutes.
        return g.iloc[start-1:stop].index
    elif method=='time_minutes':
        start=float(window['start']); stop=float(window['stop'])
        if start<0 or stop<start:
            raise InputMappingError('Invalid time-minute window')
        return g.loc[g.time_min.between(start,stop,inclusive='both')].index
    else:
        raise InputMappingError(f'Unsupported window kind: {method}')


def extract_dataset(workbook: str | Path, ds: dict[str, Any]) -> tuple[pd.DataFrame,pd.DataFrame]:
    """Return Feature Studio rows and time-point QC for one declared dataset.

    Each replicate receives the time mean, not median, for the selected window.
    We do not resample, fit Hill curves or recalibrate any fluorescence here.
    """
    if ds.get('signal_unit') != 'MEFL/cell':
        raise InputMappingError('Only already calibrated MEFL/cell supported, not RFU')
    if ds.get('system') not in {'FdeR','pLux'}:
        raise InputMappingError('system must be exactly FdeR or pLux')
    window=ds['window']
    frame=pd.read_excel(workbook, sheet_name=ds['sheet'], header=int(ds.get('header_row',1))-1,
                        engine='openpyxl')
    fmt=ds.get('format','long_timecourse')
    if fmt=='long_timecourse':
        data=_prepare_long(frame,ds)
    elif fmt=='wide_timecourse':
        data=_prepare_wide(frame,ds)
    else:
        raise InputMappingError(f'Unrecognized Excel format: {fmt}')

    for key in ['time_min','dose','signal']:
        try:
            data[key]=pd.to_numeric(data[key],errors='raise')
        except (ValueError,TypeError) as e:
            raise InputMappingError(f'Non-numeric {key}: {e}') from e
    if data[['time_min','dose','signal']].isna().any().any():
        raise InputMappingError('Missing numeric values: resolve explicitly, do not interpolate')
    if not np.isfinite(data[['time_min','dose','signal']].to_numpy(dtype=float)).all():
        raise InputMappingError('Non-finite values are not allowed')
    if (data.dose<0).any() or (data.signal<0).any():
        raise InputMappingError('Negative dose or calibrated MEFL/cell values')
    if data.replicate.isna().any() or data.replicate.astype(str).str.strip().eq('').any():
        raise InputMappingError('Every well/replicate needs a non-empty identifier')
    data.replicate=data.replicate.astype(str)
    data=data.sort_values(['dose','replicate','time_min'],kind='stable')
    if data.duplicated(['dose','replicate','time_min']).any():
        raise InputMappingError('Duplicate dose/replicate/time entries: cannot choose silently')
    if not (data.dose==0).any():
        raise InputMappingError('Zero dose is required for this export schema')
    if data.dose.nunique()<4:
        raise InputMappingError('At least four unique dose levels required')

    res=[]; checks=[]
    for (dose,rep),g in data.groupby(['dose','replicate'], sort=True):
        g=g.sort_values('time_min',kind='stable')
        ids=_window_indices(g,window)
        segment=g.loc[ids]
        if segment.empty:
            raise InputMappingError(f'No measurements in selected window: dose={dose}, rep={rep}')
        sig=float(segment.signal.mean())
        st=float(segment.signal.std(ddof=1)) if len(segment)>1 else 0.
        wc=f"{window['kind']} [{window['start']},{window['stop']}]"
        res.append({'system':ds['system'],'series':ds['series'],
                    'source':'experiment','dose':float(dose),'replicate':str(rep),
                    'signal':sig,'dose_unit':ds['dose_unit'],
                    'signal_unit':'MEFL/cell','window':wc,
                    'qc_annotation':'','origin':ds.get('origin','user_excel')})
        checks.append({'system':ds['system'],'series':ds['series'],
                       'dose':float(dose),'replicate':str(rep),
                       'n_rows_total':len(g),'n_rows_window':len(segment),
                       'time_min_start':float(segment.time_min.iloc[0]),
                       'time_min_end':float(segment.time_min.iloc[-1]),
                       'window_mean_MEFLcell':sig,'window_SD_MEFLcell':st,
                       'window_CV':st/sig if sig>0 else np.nan})
    return pd.DataFrame(res),pd.DataFrame(checks)


def extract_workbook(workbook: str | Path, config: str | Path,
                     output_dir: str | Path) -> tuple[Path,Path,Path]:
    workbook=Path(workbook); config=Path(config); output_dir=Path(output_dir)
    specification=json.loads(config.read_text(encoding='utf-8'))
    datasets=specification.get('datasets',[])
    if not datasets:
        raise InputMappingError('Config must declare at least one dataset')
    tables=[]; qc=[]
    for ds in datasets:
        part,checks=extract_dataset(workbook,ds)
        tables.append(part);qc.append(checks)
    norm=pd.concat(tables,ignore_index=True)
    if norm.duplicated(['system','series','dose','replicate']).any():
        raise InputMappingError('Duplicate output wells: make series/replicate names unique')
    qc=pd.concat(qc,ignore_index=True)
    output_dir.mkdir(parents=True,exist_ok=True)
    csv=output_dir/'feature_studio_input.csv'
    qc_path=output_dir/'window_extraction_QC.csv'
    provenance_path=output_dir/'extraction_provenance.json'
    norm.to_csv(csv,index=False)
    qc.to_csv(qc_path,index=False)
    provenance={'source_workbook_sha256':sha256_file(workbook),
                'config_sha256':sha256_file(config),
                'datasets':datasets,'output_rows':len(norm),
                'method':'arithmetic mean over indicated window PER replicate',
                'MATLAB_index_convention':'1-based, inclusive [start:stop]',
                'raw_calibration_performed':False,
                'calibration_required_upstream':True,
                'auto_outlier_rejection':False,
                'scientific_warning':'Verify dose plate map, calibration and window alignment from laboratory records.'}
    provenance_path.write_text(json.dumps(provenance,indent=2,ensure_ascii=False),encoding='utf-8')
    return csv,qc_path,provenance_path
