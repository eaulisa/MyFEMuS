import argparse
import importlib.util
import json
import math
import os
from pathlib import Path
import statistics
import subprocess
import sys

here=Path(__file__).resolve().parent
suite=here.parents[3]/'unittests/testCutFemWeightLSNative'
sys.path.insert(0,str(suite))
spec=importlib.util.spec_from_file_location('native_assessor',suite/'run.py')
assessor=importlib.util.module_from_spec(spec)
spec.loader.exec_module(assessor)
parser=argparse.ArgumentParser(description='Time every native case with interleaved batches and unchanged accuracy assertions.')
parser.add_argument('--build',type=Path,required=True)
parser.add_argument('--accuracy',type=Path,required=True)
parser.add_argument('--output',type=Path,required=True)
parser.add_argument('--samples',type=int,default=3)
parser.add_argument('--resume',action='store_true')
args=parser.parse_args()
args.output.mkdir(parents=True,exist_ok=True)
if hasattr(os,'sched_getaffinity'):
    cpu=min(os.sched_getaffinity(0))
    os.sched_setaffinity(0,{cpu})
else: cpu=None
(args.output/'environment.json').write_text(json.dumps(dict(cpu=cpu,samples=args.samples,target_batch_seconds=.01,maximum_repetitions=10000,compiler=subprocess.check_output(['c++','--version'],text=True).splitlines()[0]),indent=2)+'\n')
records=[]
for index,path in enumerate(sorted(args.accuracy.glob('e*.json'))):
    previous=json.loads(path.read_text()); case=previous['case']; output=args.output/path.name
    if args.resume and output.exists(): records.append(json.loads(output.read_text())); continue
    means=[timing['seconds']/timing['calls'] for entry in previous['results'].values() for method,timing in entry.get('timing',{}).items() if method in ['weights','physical'] and timing['calls'] and not timing['failures']]
    repetitions=max(1,min(10000,math.ceil(.01/max(means)))) if means else 1
    result=dict(case=case,repetitions=repetitions,samples={'old':[],'prepared':[]})
    for sample in range(args.samples):
        for backend in (['old','prepared'] if sample%2==0 else ['prepared','old']):
            command=[str(args.build/backend/'native_suite_speed'),str(case['element']),str(case['degree']),case['field'],case['mapping'],str(case['scale']),str(case['tolerance']),str(case['validation']),str(repetitions)]
            try:
                process=subprocess.run(command,text=True,capture_output=True,timeout=max(120,case['timeout']),env=dict(os.environ,OPENBLAS_NUM_THREADS='1'))
                data=json.loads(process.stdout) if process.returncode==0 else None
                assessment=assessor.assess(case,data) if data else dict(failures=[dict(kind='exit',code=process.returncode)])
                timing={}
                for line in process.stderr.splitlines():
                    if line.startswith('EX47_NATIVE_TIMING '): timing=json.loads(line[len('EX47_NATIVE_TIMING '):])
                result['samples'][backend].append(dict(exit=process.returncode,assessment=assessment,timing=timing,stderr=process.stderr if process.returncode or assessment['failures'] else ''))
            except Exception as error:
                result['samples'][backend].append(dict(exit=124,assessment=dict(failures=[dict(kind=type(error).__name__,message=str(error))])))
    result['timing']={}
    for backend,samples in result['samples'].items():
        result['timing'][backend]={}
        methods=set().union(*(sample.get('timing',{}).keys() for sample in samples))
        for method in methods:
            timings=[sample['timing'][method] for sample in samples if method in sample.get('timing',{})]
            result['timing'][backend][method]=dict(seconds_per_call=statistics.median(t['seconds']/t['calls'] for t in timings),minimum_seconds_per_call=min(t['seconds']/t['calls'] for t in timings),maximum_seconds_per_call=max(t['seconds']/t['calls'] for t in timings),work_per_call={key:statistics.median(t[key]/t['calls'] for t in timings) for key in timings[0] if key not in ['seconds','calls','failures']})
    result['passed']={backend:all(not sample['exit'] and not sample['assessment']['failures'] for sample in samples) for backend,samples in result['samples'].items()}
    output.write_text(json.dumps(result)+'\n'); records.append(result)
    if (index+1)%20==0: print('Timed',index+1,flush=True)
summary=dict(cases=len(records),backends={backend:dict(passed=sum(r['passed'][backend] for r in records),failed=[r['case']['id'] for r in records if not r['passed'][backend]]) for backend in ['old','prepared']})
(args.output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary),flush=True)
raise SystemExit(any(not all(r['passed'].values()) for r in records))
