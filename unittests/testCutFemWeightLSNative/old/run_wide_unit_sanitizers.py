from pathlib import Path
import argparse
import hashlib
import json
import math
import os
import subprocess
import sys

suite=Path(__file__).resolve().parent
repo=suite.parents[1]
parser=argparse.ArgumentParser(description='Run unit-basis fitting, non-graph endpoint integration, and physical geometry with ASan and UBSan.')
parser.add_argument('--skip-build',action='store_true')
parser.add_argument('--leak-detection',action='store_true')
args=parser.parse_args()
manifest={str(p.relative_to(repo)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((repo/'src/quadrature/cutfem_ls').glob('*.hpp'))}
folder=suite/'results/wide_unit_sanitizers'
folder.mkdir(parents=True,exist_ok=True)
summary_path=suite/'review_results/wide_unit_sanitizers.json'
summary_path.write_text(json.dumps({'accepted':False,'status':'running'})+'\n')
command=json.loads((suite/'build/wide_unit_validation_command.json').read_text())
command[command.index('-O2')]='-O1'
command[1:1]=['-g1','-fsanitize=address,undefined','-fno-omit-frame-pointer','-fno-sanitize-recover=all']
binary=suite/'build/wide_unit_sanitized'
command[-1]=str(binary)
(folder/'build_command.json').write_text(json.dumps(command,indent=2)+'\n')
if not args.skip_build:
    with (folder/'build.log').open('w') as log:
        subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=600)
records=[]
for element,degree in [(0,10),(1,10),(2,2),(3,2),(2,10),(3,10)]:
    if element<2:
        packet='bubble 1 1e-12 1 1 1 0\nslab 1 1e-10 0 1 1 0\ntouch_ball 1 1e-10 0 1 1 0\n'
    elif degree==2:
        packet='tiny_ball 1 1e-10 1 1 1 0\nellipse 1 1e-10 1 1 1 0\n'
        if element==3:packet+='plane 1 1e-10 0 1 2 4\n'
    else:
        packet='plane 1 1e-10 0 1 1 0\nslab 1 1e-10 0 1 1 0\n'
    process=subprocess.run([str(binary),str(element),str(degree)],input=packet,capture_output=True,text=True,timeout=600,env=dict(os.environ,ASAN_OPTIONS=f'detect_leaks={int(args.leak_detection)}:halt_on_error=1',UBSAN_OPTIONS='halt_on_error=1:print_stacktrace=1',OPENBLAS_NUM_THREADS='1'))
    (folder/f'e{element}_q{degree}.jsonl').write_text(process.stdout)
    (folder/f'e{element}_q{degree}.stderr').write_text(process.stderr)
    try:
        rows=[json.loads(line) for line in process.stdout.splitlines()]
        count=math.comb(degree+(2 if element<2 else 3),2 if element<2 else 3)
        good=process.returncode==0 and not process.stderr and len(rows)==len(packet.splitlines()) and all(not row.get('error') and len(row.get('moments',[]))==count for row in rows)
        records.append({'element':element,'degree':degree,'cases':len(rows),'moments':sum(len(row.get('moments',[])) for row in rows),'exit':process.returncode,'accepted':good})
    except Exception as error:
        records.append({'element':element,'degree':degree,'exit':process.returncode,'accepted':False,'error':str(error)})
    print(element,degree,'PASS' if records[-1]['accepted'] else 'FAIL',flush=True)
after={str(p.relative_to(repo)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((repo/'src/quadrature/cutfem_ls').glob('*.hpp'))}
summary={'accepted':all(r['accepted'] for r in records) and manifest==after,'cases':sum(r.get('cases',0) for r in records),'moments':sum(r.get('moments',0) for r in records),'records':records,'asan_ubsan_recovery':False,'leak_detection':args.leak_detection,'leak_detection_reason':'Enabled.' if args.leak_detection else 'Disabled for ptrace-compatible execution; enable with --leak-detection.','source_changed':manifest!=after,'source_sha256':manifest}
summary_path.write_text(json.dumps(summary,indent=2)+'\n')
sys.exit(0 if summary['accepted'] else 1)
