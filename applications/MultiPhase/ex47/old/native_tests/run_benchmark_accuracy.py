import argparse
import csv
import io
import json
from pathlib import Path
import subprocess
import time

parser=argparse.ArgumentParser(description='Compare every moment in the graph and degree-13 surface benchmark fixtures, outside their timing loops.')
parser.add_argument('--build',type=Path,required=True)
parser.add_argument('--output',type=Path,required=True)
args=parser.parse_args()
args.output.mkdir(parents=True,exist_ok=True)
records=[]
for name,arguments,limit in [('benchmark_graph_accuracy',['2'],4e-10),('benchmark_graph_accuracy',['10'],4e-10),('benchmark_surface_accuracy',[],1.6e-7)]:
    entry=dict(test=name,arguments=arguments,limit=limit,backends={})
    parsed={}
    for backend in ['old','prepared']:
        start=time.monotonic()
        result=subprocess.run([str(args.build/backend/name),*arguments],capture_output=True,text=True,timeout=600)
        identifier=name+('_'+'_'.join(arguments) if arguments else '')+'_'+backend
        (args.output/(identifier+'.csv')).write_text(result.stdout)
        (args.output/(identifier+'.stderr')).write_text(result.stderr)
        entry['backends'][backend]=dict(exit=result.returncode,seconds=time.monotonic()-start)
        rows={}
        for row in csv.reader(io.StringIO(result.stdout)):
            if name=='benchmark_graph_accuracy': key=tuple(row[:2]); values=list(map(float,row[2:]))
            elif row and row[0]=='M': key=tuple(row[1:6]); values=list(map(float,row[6:]))
            else: continue
            if key in rows: raise RuntimeError('Duplicated benchmark moment key')
            rows[key]=values
        parsed[backend]=rows
    checks=0;worst=0.;failed=[]
    if parsed['old'].keys()!=parsed['prepared'].keys(): failed.append(dict(kind='missing_cases'))
    for key in parsed['old'].keys()&parsed['prepared'].keys():
        a=parsed['old'][key];b=parsed['prepared'][key]
        if len(a)!=len(b): failed.append(dict(kind='missing_moments',case=key));continue
        for n,(x,y) in enumerate(zip(a,b)):
            error=abs(x-y)/max(1,abs(x),abs(y));checks+=1;worst=max(worst,error)
            if error>limit and len(failed)<8: failed.append(dict(case=key,moment=n,old=x,prepared=y,error=error))
    expected=640 if arguments==['2'] else (14080 if arguments==['10'] else 44800)
    if checks!=expected: failed.append(dict(kind='coverage',actual=checks,expected=expected))
    entry.update(checks=checks,maximum_scaled_difference=worst,failures=failed,passed=all(v['exit']==0 for v in entry['backends'].values()) and not failed)
    records.append(entry)
    print(json.dumps(entry),flush=True)
summary=dict(passed=all(r['passed'] for r in records),checks=sum(r['checks'] for r in records),records=records)
(args.output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
raise SystemExit(not summary['passed'])
