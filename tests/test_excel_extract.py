from copy import deepcopy
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from openpyxl import Workbook

from biosensor_pipeline.excel_extract import (
    InputMappingError,extract_dataset,extract_workbook,inspect_excel,sha256_file
)

ROOT=Path(__file__).resolve().parents[1]
W=ROOT/'examples/synthetic_MEFL_timecourses.xlsx'
C=ROOT/'configs/example_extraction.json'


def test_inspection_non_destructive():
    prior=sha256_file(W)
    info=inspect_excel(W,n_preview=3)
    assert [x['sheet'] for x in info]==['FdeR_R3R3','pLux_B0032']
    assert info[0]['preview_rows'][0][0]=='time_min'
    assert sha256_file(W)==prior


def test_extract_synthetic_end_to_end(tmp_path):
    prior=sha256_file(W)
    csv,qc,prov=extract_workbook(W,C,tmp_path)
    obs=pd.read_csv(csv)
    assert len(obs)==(9+7)*3
    assert set(obs.system)=={'FdeR','pLux'}
    assert set(obs.source)=={'experiment'}
    assert obs.dose.min()==0
    assert obs.signal.gt(0).all()
    assert (pd.read_csv(qc).n_rows_window>=10).all()
    assert json.loads(prov.read_text())['source_workbook_sha256']==prior
    assert sha256_file(W)==prior


def test_fder_72_to_108_inclusive_indices(tmp_path):
    ds=json.loads(C.read_text())['datasets'][0]
    data, qc=extract_dataset(W,ds)
    row=data[(data.dose==30)&(data.replicate=='rep01')].iloc[0]
    raw=pd.read_excel(W,sheet_name='FdeR_R3R3')
    group=raw[(raw.dose_uM==30)&(raw.replicate=='rep01')].sort_values('time_min')
    expected=group.iloc[71:108].MEFL_per_cell.mean()
    assert np.isclose(row.signal,expected,rtol=1e-12)
    m=qc[(qc.dose==30)&(qc.replicate=='rep01')].iloc[0]
    assert m.n_rows_window==37
    assert m.time_min_start==355 and m.time_min_end==535


def test_plux_360_to_480_minutes(tmp_path):
    ds=json.loads(C.read_text())['datasets'][1]
    data, qc=extract_dataset(W,ds)
    row=data[(data.dose==50)&(data.replicate=='rep02')].iloc[0]
    raw=pd.read_excel(W,sheet_name='pLux_B0032')
    group=raw[(raw.dose_nM==50)&(raw.replicate=='rep02')]
    expected=group[group.time_min.between(360,480)].MEFL_per_cell.mean()
    assert np.isclose(row.signal,expected,rtol=1e-12)
    assert len(data)==7*3
    assert qc.n_rows_window.eq(13).all()


def test_reject_rf_units():
    ds=json.loads(C.read_text())['datasets'][0]
    ds['signal_unit']='RFU'
    with pytest.raises(InputMappingError,match='MEFL'):
        extract_dataset(W,ds)


def test_reject_nonexistent_sheet():
    ds=json.loads(C.read_text())['datasets'][0]
    ds['sheet']='NOT_A_SHEET'
    with pytest.raises(ValueError):extract_dataset(W,ds)


def test_reject_off_grid():
    ds=json.loads(C.read_text())['datasets'][0]
    ds['window']['stop']=3000
    with pytest.raises(InputMappingError,match='exceeds'):
        extract_dataset(W,ds)


def test_explicit_wide_layout_can_extract(tmp_path):
    book=Workbook();sh=book.active;sh.title='Plate'
    sh.append(['time_min','W1','W2','W3','W4'])
    for t in range(12): sh.append([t,10+t,20+t,30+t,40+t])
    xlsx=tmp_path/'wide.xlsx';book.save(xlsx)
    ds={'system':'FdeR','series':'SYNTH_WIDE','dose_unit':'uM','signal_unit':'MEFL/cell',
        'format':'wide_timecourse','sheet':'Plate','time_column':'time_min',
        'wells':[{'column':w,'dose':d,'replicate':'rep01'} for w,d in
                 [('W1',0),('W2',30),('W3',75),('W4',100)]],
        'window':{'kind':'sample_index_1based','start':2,'stop':8}}
    frame,qc=extract_dataset(xlsx,ds)
    assert len(frame)==4
    assert frame[frame.dose==30].signal.iloc[0]==np.mean([21,22,23,24,25,26,27])
    assert set(qc.n_rows_window)=={7}


def test_duplicate_time_rejected(tmp_path):
    book=Workbook();sh=book.active;sh.title='x'
    sh.append(['time_min','dose_uM','replicate','MEFL_per_cell'])
    for d in (0,30,75,100):
        for t in (0,0,5,10):sh.append([t,d,'r1',30.])
    xlsx=tmp_path/'duplicate.xlsx';book.save(xlsx)
    ds=json.loads(C.read_text())['datasets'][0]
    ds['sheet']='x'; ds['window']={'kind':'sample_index_1based','start':1,'stop':2}
    with pytest.raises(InputMappingError,match='Duplicate'):
        extract_dataset(xlsx,ds)
