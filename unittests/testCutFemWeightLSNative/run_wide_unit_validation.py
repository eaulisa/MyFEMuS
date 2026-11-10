from pathlib import Path
import argparse
import concurrent.futures
import hashlib
import json
import math
import os
import subprocess
import sys
import time

import numpy as np
from wide_unit_reference import reference, full_moment

suite=Path(__file__).resolve().parent
repo=suite.parents[1]
parser=argparse.ArgumentParser(description='Validate complete unit-basis moment integration and fitted weights through degree ten.')
parser.add_argument('--smoke',action='store_true')
parser.add_argument('--jobs',type=int,default=3)
parser.add_argument('--element',type=int,choices=range(4))
parser.add_argument('--degree',type=int,choices=[2,5,10])
parser.add_argument('--binary',type=Path,default=suite/'build/wide_unit_validation')
parser.add_argument('--folder',default='wide_unit')
parser.add_argument('--summary')
args=parser.parse_args()
folder=suite/'results'/args.folder
folder.mkdir(parents=True,exist_ok=True)
summary_path=suite/'review_results'/(args.summary or ('wide_unit_smoke.json' if args.smoke else 'wide_unit_summary.json'))
manifest={str(p.relative_to(repo)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((repo/'src/quadrature/cutfem_ls').glob('*.hpp'))}
summary_path.write_text(json.dumps({'accepted':False,'status':'running','source_sha256':manifest},indent=2)+'\n')


def requests(element,degree):
    cases=[]
    for field in ['full','empty','boundary_low','boundary_high','plane','oblique','slab','ball','tiny_ball','touch_ball','ellipse','bubble']:
        cases.append((field,1.,1e-10,int(field in ['ball','ellipse','bubble']),1,1,0))
    if args.smoke:
        return cases
    for field in (['ball','slab','bubble'] if element<2 else ['ball','slab']):
        for scale in [-1.,2.**-20,2.**20]:
            cases.append((field,scale,1e-10,0,0,1,0))
    for field in (['ball','ellipse','bubble'] if element<2 else ['ball','ellipse']):
        for tolerance in [1e-8,1e-12]:
            cases.append((field,1.,tolerance,1,1,1,0))
    cases.append(('oblique',1.,1e-10,0,1,2,0))
    if degree==2:
        for direction in range(1,[2,3,6,4][element]+1):
            cases.append(('plane',1.,1e-10,0,1,2,direction))
    return cases


def assess(element,degree,case,data):
    field,scale,tolerance,geometry,fitting,validation,direction=case
    base,oracle_error=reference(element,field,degree)
    failures=[]
    checks=0
    worst=0.
    maximum_error=0.
    seen=set()
    if data.get('error'):
        return {'checks':0,'failures':[{'kind':'exception','message':data['error']}],'maximum_error_ratio':0.,'maximum_absolute_error':0.,'oracle_error':oracle_error}
    for row in data['moments']:
        power=tuple(row[:3])
        if power in seen or min(power)<0 or sum(power)>degree:
            failures.append({'kind':'invalid_moment','power':power})
            continue
        seen.add(power)
        expected=base[power].copy()
        if scale<0:
            expected[0]=full_moment(element,power)-expected[0]
        expected[1]/=abs(scale)
        tests=[(3,expected[0]),(4,expected[1]),(5,expected[2])]
        if geometry:
            tests.append((6,expected[3]))
        if fitting:
            tests.extend([(7,expected[0]),(8,expected[1]),(9,expected[2])])
            if geometry:
                tests.append((10,expected[3]))
                if sum(power)<degree:
                    raised=(power[0]+1,power[1],power[2])
                    extra=base[raised][0]
                    if scale<0:
                        extra=full_moment(element,raised)-extra
                    tests.append((11,expected[0]+.2*extra))
        for column,target in tests:
            actual=row[column]
            error=abs(actual-target) if actual is not None else math.inf
            reference_error=oracle_error/(abs(scale) if column in [4,8] else 1.)
            limit=2*tolerance*(1+abs(target))+8*reference_error+2e-12
            ratio=error/limit
            checks+=1
            maximum_error=max(maximum_error,error)
            worst=max(worst,ratio)
            if not math.isfinite(error) or ratio>1:
                if len(failures)<8:
                    failures.append({'kind':'moment','power':power,'column':column,'actual':actual,'expected':target,'error':error,'limit':limit})
    expected_count=math.comb(degree+(2 if element<2 else 3),2 if element<2 else 3)
    if len(seen)!=expected_count:
        failures.append({'kind':'missing_moments','actual':len(seen),'expected':expected_count})
    if data.get('fallback_cells',0):
        failures.append({'kind':'unresolved_graph_cells','count':data['fallback_cells']})
    return {'checks':checks,'failures':failures,'maximum_error_ratio':worst,'maximum_absolute_error':maximum_error,'oracle_error':oracle_error}


def group(item):
    element,degree=item
    cases=requests(element,degree)
    packet=''.join(' '.join(map(str,c))+'\n' for c in cases)
    started=time.monotonic()
    output=folder/f'e{element}_q{degree}.jsonl'
    try:
        with output.open('w') as sink:
            process=subprocess.run([str(args.binary),str(element),str(degree)],input=packet,text=True,stdout=sink,stderr=subprocess.PIPE,timeout=900,env=dict(os.environ,OPENBLAS_NUM_THREADS='1',ASAN_OPTIONS=os.environ.get('ASAN_OPTIONS','detect_leaks=0')))
        rows=[json.loads(line) for line in output.read_text().splitlines()]
        records=[]
        for case,row in zip(cases,rows):
            assessment=assess(element,degree,case,row)
            records.append({'case':case,'assessment':assessment,'seconds':row.get('seconds'),'setup_seconds':row.get('setup_seconds',0),'splits':{k:row.get(k,0) for k in ['quad_splits','triangle_splits','cell_splits','quad_calls','tri_calls','mapped_panels']}})
        errors=[]
        if process.returncode or len(rows)!=len(cases):
            errors.append({'exit':process.returncode,'rows':len(rows),'expected_rows':len(cases),'stderr':process.stderr[-2000:]})
    except Exception as error:
        records=[]
        errors=[{'error':str(error)}]
    failed=[r for r in records if r['assessment']['failures']]
    report={'element':element,'degree':degree,'cases':len(cases),'records':records,'failures':errors+failed,'seconds':time.monotonic()-started}
    (folder/f'e{element}_q{degree}_assessment.json').write_text(json.dumps(report,indent=2)+'\n')
    print(['QUAD','TRI','HEX','TET'][element],degree,'PASS' if not report['failures'] else 'FAIL','cases',len(cases),'checks',sum(r['assessment']['checks'] for r in records),'seconds',round(report['seconds'],3),flush=True)
    if report['failures']:
        print(json.dumps(report['failures'][:2]),flush=True)
    return report


elements=[args.element] if args.element is not None else range(4)
degrees=[args.degree] if args.degree is not None else ([2] if args.smoke else [2,5,10])
with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
    records=list(pool.map(group,[(e,q) for e in elements for q in degrees]))
after={str(p.relative_to(repo)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((repo/'src/quadrature/cutfem_ls').glob('*.hpp'))}
summary={'accepted':all(not r['failures'] for r in records) and manifest==after,'maximum_degree':10,'groups':len(records),'cases':sum(r['cases'] for r in records),'checks':sum(v['assessment']['checks'] for r in records for v in r['records']),'maximum_error_ratio':max((v['assessment']['maximum_error_ratio'] for r in records for v in r['records']),default=0),'failures':[{'element':r['element'],'degree':r['degree'],'failures':r['failures']} for r in records if r['failures']],'source_changed':manifest!=after,'source_sha256':manifest,'split_totals':{k:sum(v['splits'][k] for r in records for v in r['records']) for k in ['quad_splits','triangle_splits','cell_splits','quad_calls','tri_calls','mapped_panels']}}
summary_path.write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps({k:v for k,v in summary.items() if k not in ['failures','source_sha256']},indent=2),flush=True)
sys.exit(0 if summary['accepted'] else 1)
