from pathlib import Path
import argparse
import concurrent.futures
import hashlib
import json
import sys

import run

suite=Path(__file__).resolve().parent
repo=suite.parents[1]
parser=argparse.ArgumentParser(description='Check degree-ten double-output weights using the full native FE driver.')
parser.add_argument('--jobs',type=int,default=2)
parser.add_argument('--2d-only',action='store_true',dest='two_dimensional')
parser.add_argument('--folder',default='wide_double_output')
parser.add_argument('--summary',default='wide_unit_double_output.json')
args=parser.parse_args()
run.folder=suite/'results'/args.folder
run.folder.mkdir(parents=True,exist_ok=True)
manifest={str(p.relative_to(repo)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((repo/'src/quadrature/cutfem_ls').glob('*.hpp'))}
cases=[]
for element in range(1,3 if args.two_dimensional else 5):
    for field,mapping in [('plane','identity'),('graph','identity'),('sphere','identity'),('plane','stretch')]:
        cases.append({'id':f'e{element}_q10_{field}_{mapping}','element':element,'degree':10,'field':field,'mapping':mapping,'scale':1.,'tolerance':1e-12 if field=='plane' else 1e-10,'validation':1,'timeout':180})
with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
    records=list(pool.map(run.run,cases))
after={str(p.relative_to(repo)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((repo/'src/quadrature/cutfem_ls').glob('*.hpp'))}
summary={'accepted':all(not r['assessment']['failures'] for r in records) and manifest==after,'cases':len(records),'checks':sum(r['assessment']['assertions'] for r in records),'records':[{'case':r['case'],'seconds':r['runtime_seconds'],'assessment':r['assessment']} for r in records],'source_changed':manifest!=after,'source_sha256':manifest}
(suite/'review_results'/args.summary).write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps({k:v for k,v in summary.items() if k not in ['records','source_sha256']},indent=2),flush=True)
for record in records:
    if record['assessment']['failures']:print(record['case']['id'],record['assessment']['failures'][:3],flush=True)
sys.exit(0 if summary['accepted'] else 1)
