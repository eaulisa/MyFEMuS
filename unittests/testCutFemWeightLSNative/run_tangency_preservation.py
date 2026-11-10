from pathlib import Path
import argparse
import csv
import hashlib
import json
import os
import statistics
import subprocess
import sys
import time

suite=Path(__file__).resolve().parent
repo=suite.parents[1]
parser=argparse.ArgumentParser(description='Check unchanged graph outputs, sustained 2D graph timing, and unchanged 3D integration.')
parser.add_argument('--graph-rounds',type=int,default=3)
parser.add_argument('--skip-3d',action='store_true')
parser.add_argument('--skip-build',action='store_true')
args=parser.parse_args()
folder=suite/'results/tangency_preservation'
folder.mkdir(parents=True,exist_ok=True)
manifest={str(p.relative_to(repo)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((repo/'src/quadrature/cutfem_ls').glob('*.hpp'))}
baseline=suite/'results/tangency_strips/baseline_source'
source=(suite/'benchmark_graph_dimensions.cpp').read_text()
old='for(unsigned element = 0; element < 4; ++element) DimensionMeasure'
assert source.count(old)==1
source=source.replace(old,'for(unsigned element = 0; element < 2; ++element) DimensionMeasure')
probe=folder/'graph_source.cpp'
probe.write_text(source)
common=json.loads((suite/'build/tangency_strip_benchmark_command.json').read_text())
binaries={label:suite/'build'/('tangent_graph_'+label) for label in ['baseline','candidate']}
for label,root in [('baseline',baseline),('candidate',repo)]:
    if args.skip_build:continue
    command=list(common)
    command[command.index('-O2')]='-O3'
    command[1:1]=['-DFEMUS_UNIT_COORDINATE_MOMENTS','-DNDEBUG','-march=x86-64','-mtune=generic']
    command=[('-I'+str(root/'src/quadrature/cutfem_ls')) if p=='-I'+str(repo/'src/quadrature/cutfem_ls') else (str(probe) if p==str(suite/'tangency_strip_benchmark.cpp') else p) for p in command]
    binaries[label]=suite/'build'/('tangent_graph_'+label)
    command[-1]=str(binaries[label])
    with (folder/(label+'_build.log')).open('w') as log:subprocess.run(command,check=True,stdout=log,stderr=subprocess.STDOUT,timeout=300)
    print('compiled graph probe',label,flush=True)
    if label=='baseline' and not args.skip_3d:
        wide=json.loads((suite/'build/wide_unit_validation_command.json').read_text())
        wide=[('-I'+str(root/'src/quadrature/cutfem_ls')) if p=='-I'+str(repo/'src/quadrature/cutfem_ls') else p for p in wide]
        wide[-1]=str(suite/'build/wide_unit_baseline')
        with (folder/'baseline_wide_build.log').open('w') as log:subprocess.run(wide,check=True,stdout=log,stderr=subprocess.STDOUT,timeout=300)
records={}
cpu=str(min(os.sched_getaffinity(0)))
for round_index in range(args.graph_rounds):
    for q in [2,5,10]:
        for label in (['baseline','candidate'] if round_index%2==0 else ['candidate','baseline']):
            text=subprocess.check_output(['taskset','-c',cpu,str(binaries[label]),str(q)],text=True,timeout=120)
            (folder/f'{label}_q{q}_{round_index}.csv').write_text(text)
            rows=list(csv.DictReader(text.splitlines()))
            assert len(rows)==20
            records[label,q,round_index]=rows
            print('graph',label,q,round_index+1,'PASS',flush=True)
graphs=[]
failures=[]
for q in [2,5,10]:
    for element in ['QUAD','TRI']:
        values={label:statistics.mean(float(row['seconds']) for r in range(args.graph_rounds) for row in records[label,q,r] if row['element']==element) for label in ['baseline','candidate']}
        ratios=[sum(float(row['seconds']) for row in records['candidate',q,r] if row['element']==element)/sum(float(row['seconds']) for row in records['baseline',q,r] if row['element']==element) for r in range(args.graph_rounds)]
        graphs.append({'degree':q,'element':element,'baseline_ms':1000*values['baseline'],'candidate_ms':1000*values['candidate'],'candidate_over_baseline':values['candidate']/values['baseline'],'round_ratios':ratios})
        if values['candidate']>1.08*values['baseline'] and min(ratios)>1.05:failures.append({'kind':'sustained_graph_timing_regression','result':graphs[-1]})
unchanged=[]
checks=0
if not args.skip_3d:
    packet=''.join(f'{field} 1 1e-10 0 0 1 0\n' for field in ['plane','slab','ball','tiny_ball'])
    for element in [2,3]:
        for q in [2,10]:
            results={}
            for label,binary in [('baseline','wide_unit_baseline'),('candidate','wide_unit_validation')]:
                p=subprocess.run([str(suite/'build'/binary),str(element),str(q)],input=packet,capture_output=True,text=True,timeout=300)
                (folder/f'3d_{label}_e{element}_q{q}.jsonl').write_text(p.stdout)
                assert p.returncode==0
                results[label]=[json.loads(line,parse_float=str) for line in p.stdout.splitlines()]
            assert len(results['baseline'])==len(results['candidate'])==4
            for old,new in zip(results['baseline'],results['candidate']):
                good=not old.get('error') and not new.get('error') and old['moments']==new['moments'] and new.get('tangent_attempts',0)==0
                checks+=sum(3 for row in new.get('moments',[]))
                unchanged.append({'element':element,'degree':q,'field':new['field'],'unchanged':good})
                if not good:failures.append({'kind':'3d_changed','case':unchanged[-1]})
            print('3D',element,q,'checked',flush=True)
after={str(p.relative_to(repo)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((repo/'src/quadrature/cutfem_ls').glob('*.hpp'))}
result={'accepted':not failures and manifest==after,'graph_timing':graphs,'graph_rounds':args.graph_rounds,'graphs_per_element':10,'cpu':int(cpu),'three_dimensional_cases':unchanged,'three_dimensional_checks':checks,'failures':failures,'source_changed':manifest!=after,'source_sha256':manifest}
(suite/'review_results/tangency_preservation.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in ['source_sha256','three_dimensional_cases']},indent=2),flush=True)
sys.exit(0 if result['accepted'] else 1)
