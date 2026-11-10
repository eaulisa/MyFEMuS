from pathlib import Path
import argparse
import concurrent.futures
import csv
import hashlib
import io
import json
import math
import os
import subprocess
import sys
import time

suite=Path(__file__).resolve().parent
repo=suite.parents[1]
parser=argparse.ArgumentParser(description='Check fitted volume and physical-surface weights on refined QUAD/TRI/HEX/TET meshes.')
parser.add_argument('--smoke',action='store_true')
parser.add_argument('--jobs',type=int,default=3)
parser.add_argument('--element',type=int,choices=range(4))
parser.add_argument('--degree',type=int,choices=[2,5,10])
parser.add_argument('--folder',default='wide_unit_mesh')
parser.add_argument('--summary')
args=parser.parse_args()
folder=suite/'results'/args.folder
folder.mkdir(parents=True,exist_ok=True)
manifest={str(p.relative_to(repo)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((repo/'src/quadrature/cutfem_ls').glob('*.hpp'))}


def run(item):
    element,degree=item
    maximum=(8 if element<2 else 2) if args.smoke else (64 if element<2 else 8)
    start=time.monotonic()
    try:
        output=folder/f'e{element}_q{degree}.csv'
        with output.open('w') as sink:
            p=subprocess.run([str(suite/'build/wide_unit_mesh'),str(element),str(degree),str(maximum)],stdout=sink,stderr=subprocess.PIPE,text=True,timeout=1800,env=dict(os.environ,OPENBLAS_NUM_THREADS='1'))
        content=output.read_text()
        (folder/f'e{element}_q{degree}.stderr').write_text(p.stderr)
        rows=list(csv.DictReader(io.StringIO(content)))
        expected=int(math.log2(maximum))
        failures=[]
        if p.returncode or len(rows)!=expected or any(row['status']!='PASS' for row in rows):
            failures.append({'exit':p.returncode,'rows':len(rows),'expected':expected,'stderr':p.stderr[-2000:]})
        for row in rows:
            dimension=2 if element<2 else 3
            cells=int(row['divisions'])**dimension*(1 if element%2==0 else (2 if dimension==2 else 6))
            if int(row['cells'])!=cells or int(row['cut_cells'])==0:
                failures.append({'kind':'mesh_coverage','row':row})
        result={'element':element,'degree':degree,'rows':rows,'failures':failures,'seconds':time.monotonic()-start}
    except Exception as error:
        result={'element':element,'degree':degree,'rows':[],'failures':[str(error)],'seconds':time.monotonic()-start}
    print(['QUAD','TRI','HEX','TET'][element],degree,'PASS' if not result['failures'] else 'FAIL','levels',len(result['rows']),'seconds',round(result['seconds'],3),flush=True)
    if result['failures']:print(result['failures'],flush=True)
    return result


elements=[args.element] if args.element is not None else range(4)
degrees=[args.degree] if args.degree is not None else ([2] if args.smoke else [2,5,10])
with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
    records=list(pool.map(run,[(e,q) for e in elements for q in degrees]))
after={str(p.relative_to(repo)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((repo/'src/quadrature/cutfem_ls').glob('*.hpp'))}
rows=[row for r in records for row in r['rows']]
summary={'accepted':all(not r['failures'] for r in records) and manifest==after,'groups':len(records),'mesh_levels':len(rows),'checks':4*len(rows),'maximum_requested_degree':10,'maximum_errors':{name:max((float(row[name]) for row in rows),default=0) for name in ['volume_error','surface_error','volume_x2_error','surface_x2_error']},'records':records,'source_sha256':manifest,'source_changed':manifest!=after}
name=args.summary or ('wide_unit_mesh_smoke.json' if args.smoke else 'wide_unit_mesh_summary.json')
(suite/'review_results'/name).write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps({k:v for k,v in summary.items() if k not in ['records','source_sha256']},indent=2),flush=True)
sys.exit(0 if summary['accepted'] else 1)
