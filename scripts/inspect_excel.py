"""Print sheet/header inventory; never change the source Excel file."""
import argparse
import json
from biosensor_pipeline.excel_extract import inspect_excel

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('workbook')
    p.add_argument('--rows',type=int,default=5)
    args=p.parse_args()
    print(json.dumps(inspect_excel(args.workbook,n_preview=args.rows),indent=2,ensure_ascii=False))
