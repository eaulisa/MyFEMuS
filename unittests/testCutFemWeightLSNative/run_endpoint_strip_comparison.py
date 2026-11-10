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
import numpy as np

from tangency_strip_reference import reference,femus_reference

suite=Path(__file__).resolve().parent
repo=suite.parents[1]
parser=argparse.ArgumentParser(description='Compare the isolated 2D endpoint-strip quadrature experiment with its pinned baseline.')
parser.add_argument('--quick',action='store_true')
parser.add_argument('--skip-build',action='store_true')
parser.add_argument('--rounds',type=int,default=3)
parser.add_argument('--measure',action='store_true')
args=parser.parse_args()
baseline='1a394fdc97bedc5a6cff9b1e6c68f16261ba57e9'
folder=suite/'results/endpoint_strips'
folder.mkdir(parents=True,exist_ok=True)
summary_path=suite/'review_results/endpoint_strip_comparison.json'
manifest={str(p.relative_to(repo)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((repo/'src/quadrature/cutfem_ls').glob('*.hpp'))}
summary_path.write_text(json.dumps({'accepted':False,'status':'running','source_sha256':manifest},indent=2)+'\n')
binaries={name:suite/'build'/('endpoint_'+name) for name in ['baseline','candidate']}
if not args.skip_build:
    original=folder/'baseline_source'
    original.mkdir(exist_ok=True)
    archive=subprocess.check_output(['git','archive',baseline,'src/quadrature/cutfem_ls'],cwd=repo)
    with tarfile.open(fileobj=io.BytesIO(archive)) as stream:stream.extractall(original,filter='data')
    for name in binaries:
        command=json.loads((suite/'build/endpoint_strip_benchmark_command.json').read_text())
        if '-O2' in command:command[command.index('-O2')]='-O3'
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
            timings.setdefault(key,[]).append({'name':row[5],'round':int(row[6]),'seconds':float(row[7]),'status':row[8],'attempts':int(row[9]),'completed':int(row[10]),'rejected':int(row[11]),'strip_samples':int(row[12]),'strip_bisections':int(row[13]),'bulk_samples':int(row[14]),'geometry_builds':int(row[15]),'geometry_hits':int(row[16]),'minimum_delta':float(row[17]),'maximum_delta':float(row[18]),'bulk_bisections':int(row[19]),'strip_depth':int(row[20]),'bulk_depth':int(row[21]),'rule15':int(row[22]),'rule31':int(row[23]),'quad_splits':int(row[24]),'triangle_splits':int(row[25])})
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
            maximum_errors=[0.,0.]
            maximum_difference=[0.,0.]
            for power,values in actual.get(key,{}).items():
                targets_pair=targets[power][[0,{'raw':1,'surface':2,'physical':3}[kind]]]
                for component,(value,target) in enumerate(zip(values,targets_pair)):
                    error=abs(float(value)-float(target))
                    limit=1e-10*(1+abs(float(target)))+8*discrepancy+3e-12
                    ratio=error/limit
                    maximum_errors[component]=max(maximum_errors[component],error)
                    if key in previous and power in previous[key]:maximum_difference[component]=max(maximum_difference[component],abs(float(value)-float(previous[key][power][component])))
                    checks+=1
                    error_ratio=max(error_ratio,ratio)
                    if not math.isfinite(error) or ratio>1:
                        if len(record['failures'])<5:record['failures'].append({'kind':'moment','power':power,'actual':value,'expected':float(target),'error':error,'limit':limit})
            record['maximum_error_ratio']=error_ratio
            record['maximum_absolute_error']=maximum_errors
            record['maximum_difference_from_baseline']=maximum_difference
            if element==0 and key in actual:
                signed_targets,signed_discrepancy=femus_reference(element,coefficients[element,field],degree)
                signed_maximum=[0.,0.]
                signed_ratio=0.
                for power in actual[key]:
                    i,j=power
                    target_pair=signed_targets[power][[0,{'raw':1,'surface':2,'physical':3}[kind]]]
                    for component,target in enumerate(target_pair):
                        value=np.longdouble(0)
                        for a in range(i+1):
                            for b in range(j+1):
                                value+=np.longdouble(math.comb(i,a)*math.comb(j,b)*2**(a+b)*(-1)**(i+j-a-b))*np.longdouble(actual[key][a,b][component])
                        value*=4 if component==0 or kind=='raw' else (2 if kind=='surface' else 1)
                        error=float(abs(value-target))
                        limit=1e-10*(1+abs(float(target)))+8*signed_discrepancy+3e-12
                        signed_ratio=max(signed_ratio,error/limit)
                        signed_maximum[component]=max(signed_maximum[component],error)
                        checks+=1
                        if not math.isfinite(error) or error>limit:
                            if len(record['failures'])<5:record['failures'].append({'kind':'signed_reference_moment','power':power,'component':component,'actual':float(value),'expected':float(target),'error':error,'limit':limit})
                record['signed_reference_maximum_error']=signed_maximum
                record['signed_reference_error_ratio']=signed_ratio
                maximum_ratio=max(maximum_ratio,signed_ratio)
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
