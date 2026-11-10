from pathlib import Path
import argparse
import csv
import hashlib
import io
import json
import math
import os
import statistics
import subprocess
import sys
import tarfile
import time

from tangency_strip_reference import reference

suite=Path(__file__).resolve().parent
repo=suite.parents[1]
parser=argparse.ArgumentParser(description='Compare the isolated 2D tangency-strip experiment with its pinned baseline.')
parser.add_argument('--quick',action='store_true')
parser.add_argument('--skip-build',action='store_true')
parser.add_argument('--rounds',type=int,default=3)
parser.add_argument('--measure',action='store_true')
args=parser.parse_args()
baseline='1a394fdc97bedc5a6cff9b1e6c68f16261ba57e9'
folder=suite/'results/tangency_strips'
folder.mkdir(parents=True,exist_ok=True)
summary_path=suite/'review_results/tangency_strip_experiment.json'
manifest={str(p.relative_to(repo)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((repo/'src/quadrature/cutfem_ls').glob('*.hpp'))}
summary_path.write_text(json.dumps({'accepted':False,'status':'running','source_sha256':manifest},indent=2)+'\n')
binaries={name:suite/'build'/('tangent_'+name) for name in ['baseline','candidate']}
if not args.skip_build:
    original=folder/'baseline_source'
    original.mkdir(exist_ok=True)
    archive=subprocess.check_output(['git','archive',baseline,'src/quadrature/cutfem_ls'],cwd=repo)
    with tarfile.open(fileobj=io.BytesIO(archive)) as stream:stream.extractall(original,filter='data')
    for name in binaries:
        command=json.loads((suite/'build/tangency_strip_benchmark_command.json').read_text())
        command[command.index('-O2')]='-O3'
        command[1:1]=['-DNDEBUG','-march=x86-64','-mtune=generic']
        if name=='baseline':command=[('-I'+str(original/'src/quadrature/cutfem_ls')) if p=='-I'+str(repo/'src/quadrature/cutfem_ls') else p for p in command]
        command[-1]=str(binaries[name])
        (folder/(name+'_command.json')).write_text(json.dumps(command,indent=2)+'\n')
        with (folder/(name+'_build.log')).open('w') as log:subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=300)
        print('compiled',name,flush=True)


def parse(text):
    timings={}
    moments={}
    coefficients={}
    for row in csv.reader(text.splitlines()):
        if row[0]=='T':
            key=(int(row[1]),int(row[2]),row[3],int(row[4]))
            timings.setdefault(key,[]).append({'name':row[5],'round':int(row[6]),'seconds':float(row[7]),'status':row[8],'attempts':int(row[9]),'completed':int(row[10]),'rejected':int(row[11]),'analytic_panels':int(row[12]),'model_refinements':int(row[13]),'bulk_samples':int(row[14]),'geometry_builds':int(row[15]),'geometry_hits':int(row[16]),'minimum_delta':float(row[17]),'maximum_delta':float(row[18]),'quad_splits':int(row[19]),'triangle_splits':int(row[20])})
        elif row[0]=='M':
            key=(int(row[1]),int(row[2]),row[3],int(row[4]))
            power=(int(row[5]),int(row[6]))
            if power in moments.setdefault(key,{}):raise RuntimeError('Duplicated moment.')
            moments[key][power]=(row[7],row[8])
        elif row[0]=='C':coefficients[int(row[1]),int(row[2])]=tuple(row[3:])
    return timings,moments,coefficients


records=[]
failures=[]
checks=0
graph_cases=0
maximum_ratio=0.
cpu=str(min(os.sched_getaffinity(0)))
degrees=[2] if args.quick else [2,5,10]
modes=['surface'] if args.quick else ['raw','surface','physical']
for q in degrees:
    for mode_index,mode in enumerate(modes):
        outputs={}
        for name in (['baseline','candidate'] if (q+mode_index)%2==0 else ['candidate','baseline']):
            start=time.monotonic()
            p=subprocess.run(['taskset','-c',cpu,str(binaries[name]),str(q),mode,str(args.rounds),'1e-10',str(int(args.measure))],capture_output=True,text=True,timeout=600,env=dict(os.environ,OPENBLAS_NUM_THREADS='1'))
            (folder/f'{name}_q{q}_{mode}.csv').write_text(p.stdout)
            (folder/f'{name}_q{q}_{mode}.stderr').write_text(p.stderr)
            outputs[name]=parse(p.stdout)
            print(name,q,mode,'exit',p.returncode,'seconds',round(time.monotonic()-start,3),flush=True)
            if name=='candidate' and p.returncode:failures.append({'kind':'process','degree':q,'mode':mode,'exit':p.returncode,'stderr':p.stderr[-3000:]})
        candidate,actual,coefficients=outputs['candidate']
        old,previous,old_coefficients=outputs['baseline']
        if len(candidate)!=45 or len(old)!=45:failures.append({'kind':'missing_cases','degree':q,'mode':mode,'candidate':len(candidate),'baseline':len(old)})
        for key,rows in candidate.items():
            element,degree,kind,field=key
            record={'element':element,'degree':degree,'mode':kind,'field':field,'name':rows[0]['name'],'candidate_ms':1000*statistics.median(r['seconds'] for r in rows),'stats':rows[0],'failures':[]}
            if len(rows)!=args.rounds or any(r['status']!='PASS' for r in rows):record['failures'].append({'kind':'case_failure'})
            if key in old and all(r['status']=='PASS' for r in old[key]):
                record['baseline_ms']=1000*statistics.median(r['seconds'] for r in old[key])
                record['speedup']=record['baseline_ms']/record['candidate_ms']
            else:record['baseline_failure']=True
            expected_count=math.comb(q+2,2)
            if len(actual.get(key,{}))!=expected_count:record['failures'].append({'kind':'missing_moments','count':len(actual.get(key,{})),'expected':expected_count})
            if field<10:
                graph_cases+=1
                if any(r['attempts'] for r in rows):record['failures'].append({'kind':'graph_entered_new_path'})
                if actual.get(key)!=previous.get(key):record['failures'].append({'kind':'graph_not_bitwise_equal'})
            targets,discrepancy=reference(element,coefficients[element,field],degree)
            record['reference_discrepancy']=discrepancy
            error_ratio=0.
            for power,values in actual.get(key,{}).items():
                targets_pair=targets[power][[0,{'raw':1,'surface':2,'physical':3}[kind]]]
                for value,target in zip(values,targets_pair):
                    error=abs(float(value)-float(target))
                    limit=4e-10*(1+abs(float(target)))+8*discrepancy+3e-12
                    ratio=error/limit
                    checks+=1
                    error_ratio=max(error_ratio,ratio)
                    if not math.isfinite(error) or ratio>1:
                        if len(record['failures'])<5:record['failures'].append({'kind':'moment','power':power,'actual':value,'expected':float(target),'error':error,'limit':limit})
            record['maximum_error_ratio']=error_ratio
            maximum_ratio=max(maximum_ratio,error_ratio)
            records.append(record)
            if record['failures']:failures.append(record)
        print('assessed',q,mode,'failures',len(failures),flush=True)
after={str(p.relative_to(repo)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((repo/'src/quadrature/cutfem_ls').glob('*.hpp'))}
summary={'accepted':not failures and manifest==after,'baseline_sha':baseline,'maximum_requested_degree':max(degrees),'cases':len(records),'checks':checks,'graph_cases_bitwise':graph_cases,'maximum_error_ratio':maximum_ratio,'failures':failures,'records':records,'source_changed':manifest!=after,'source_sha256':manifest,'cpu':int(cpu),'sustained_measurement':args.measure,'timing':'Median of five batches lasting at least 8ms each, alternating nearby level sets to force geometry reconstruction; fixed CPU; identical O3 flags; setup and output formatting excluded.' if args.measure else 'Median of repeated full-cell nodal integrations; fixed CPU; identical O3 flags; no fitted-weight setup or output formatting in timed region.'}
summary_path.write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps({k:v for k,v in summary.items() if k not in ['records','source_sha256','failures']},indent=2),flush=True)
if failures:print(json.dumps(failures[:3],indent=2),flush=True)
sys.exit(0 if summary['accepted'] else 1)
