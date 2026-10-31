"""Generate Feature Studio CSV + time-window QC + SHA256 provenance."""
import argparse
from biosensor_pipeline.excel_extract import extract_workbook

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--workbook',required=True)
    p.add_argument('--config',required=True)
    p.add_argument('--out',default='outputs')
    args=p.parse_args()
    for f in extract_workbook(args.workbook,args.config,args.out):
        print('Saved:',f)
