"""List archived notebook evidence; never imports or executes notebook code."""
from pathlib import Path
import argparse,csv,json
ROOT=Path(__file__).resolve().parents[1]
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--catalog',action='store_true');args=ap.parse_args()
    if args.catalog:
        with (ROOT/'docs/notebook_catalog.csv').open(encoding='utf-8-sig',newline='') as f:
            for row in csv.DictReader(f):
                if row.get('GitHub_candidate')=='NO':continue
                print(json.dumps({k:row[k] for k in ['Step','Release_file','Kernel','Python_metadata_version','JSON_structure_QA','Run_all_verified']},ensure_ascii=True))
    else:ap.print_help()
