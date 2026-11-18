import argparse
import concurrent.futures
import importlib.util
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time

here=Path(__file__).resolve().parent
repo=here.parents[3]
suite=repo/'unittests/testCutFemWeightLSNative'
sys.path.insert(0,str(suite))
import cases
spec=importlib.util.spec_from_file_location('native_original_assessor',suite/'run.py')
assessor=importlib.util.module_from_spec(spec)
spec.loader.exec_module(assessor)
parser=argparse.ArgumentParser()
parser.add_argument('--build',type=Path,required=True)
parser.add_argument('--output',type=Path,required=True)
parser.add_argument('--jobs',type=int,default=2)
parser.add_argument('--resume',action='store_true')
parser.add_argument('--selection',default='original')
args=parser.parse_args()
args.output.mkdir(parents=True,exist_ok=True)
selected=assessor.portable_cases(args.selection)

def execute(case):
    output=args.output/(case['id']+'.json')
    if args.resume and output.exists(): return json.loads(output.read_text())
    record={'case':case,'results':{}}
    for backend in ['old','prepared']:
        command=[str(args.build/backend/'native_suite'),str(case['element']),str(case['degree']),case['field'],case['mapping'],str(case['scale']),str(case['tolerance']),str(case['validation'])]
        started=time.monotonic()
        try:
            result=subprocess.run(command,capture_output=True,text=True,timeout=max(120,case['timeout']),env=dict(os.environ,OPENBLAS_NUM_THREADS='1'))
            data=json.loads(result.stdout) if result.returncode==0 else None
            if data is not None: assessment=assessor.assess(case,data)
            else: assessment={'assertions':0,'failures':[{'kind':'exit','code':result.returncode,'message':result.stderr[-3000:]}]}
            timing={}
            for line in result.stderr.splitlines():
                if line.startswith('EX47_NATIVE_TIMING '): timing=json.loads(line[len('EX47_NATIVE_TIMING '):])
            entry={'exit':result.returncode,'wall_seconds':time.monotonic()-started,'assessment':assessment,'data':data,'timing':timing,'stderr':result.stderr}
        except Exception as error:
            entry={'exit':124 if isinstance(error,subprocess.TimeoutExpired) else 1,'wall_seconds':time.monotonic()-started,'assessment':{'assertions':0,'failures':[{'kind':type(error).__name__,'message':str(error)}]}}
        record['results'][backend]=entry
    record['maximum_class_difference']=0
    left=record['results']['old'].get('data'); right=record['results']['prepared'].get('data')
    if left and right:
        for a,b in zip(left['moments'],right['moments']):
            if a['powers']!=b['powers']: raise RuntimeError('Moment ordering changed')
            for x,y in zip(a['values'],b['values']):
                if x is not None and y is not None: record['maximum_class_difference']=max(record['maximum_class_difference'],abs(x-y)/max(1,abs(x),abs(y)))
    output.write_text(json.dumps(record)+'\n')
    flagged=[backend for backend,entry in record['results'].items() if entry['assessment']['failures']]
    if flagged: print('FLAG',case['id'],','.join(flagged),json.dumps({backend:record['results'][backend]['assessment']['failures'][:1] for backend in flagged}),flush=True)
    return record

records=[]
with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
    futures={pool.submit(execute,case):case for case in selected}
    for future in concurrent.futures.as_completed(futures):
        records.append(future.result())
        if len(records)%20==0 or len(records)==len(selected): print('Completed',len(records),'/',len(selected),flush=True)
summary={'cases':len(records),'backends':{backend:{'flagged':sum(bool(r['results'][backend]['assessment']['failures']) for r in records),'assertions':sum(r['results'][backend]['assessment']['assertions'] for r in records)} for backend in ['old','prepared']},'maximum_class_difference':max(r['maximum_class_difference'] for r in records)}
(args.output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary),flush=True)
