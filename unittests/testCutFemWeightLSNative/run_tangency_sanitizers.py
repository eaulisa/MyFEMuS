from pathlib import Path
import argparse
import hashlib
import json
import os
import subprocess
import sys

suite=Path(__file__).resolve().parent
repo=suite.parents[1]
parser=argparse.ArgumentParser(description='Check the 2D tangent-strip kernels with ASan and UBSan.')
parser.add_argument('--leak-detection',action='store_true')
args=parser.parse_args()
folder=suite/'results/tangency_sanitizers'
folder.mkdir(parents=True,exist_ok=True)
manifest={str(p.relative_to(repo)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((repo/'src/quadrature/cutfem_ls').glob('*.hpp'))}
command=json.loads((suite/'build/tangency_strip_benchmark_command.json').read_text())
command[command.index('-O2')]='-O1'
command[1:1]=['-g1','-fsanitize=address,undefined','-fno-sanitize-recover=all','-fno-omit-frame-pointer']
command[-1]=str(suite/'build/tangent_asan')
with (folder/'build.log').open('w') as log:subprocess.run(command,check=True,stdout=log,stderr=subprocess.STDOUT,timeout=600)
records=[]
for degree,mode in [(2,'raw'),(2,'surface'),(10,'physical')]:
    process=subprocess.run([command[-1],str(degree),mode,'1','1e-10'],capture_output=True,text=True,timeout=300,env=dict(os.environ,ASAN_OPTIONS=f'detect_leaks={int(args.leak_detection)}:halt_on_error=1',UBSAN_OPTIONS='halt_on_error=1:print_stacktrace=1'))
    (folder/f'q{degree}_{mode}.csv').write_text(process.stdout)
    (folder/f'q{degree}_{mode}.stderr').write_text(process.stderr)
    records.append({'degree':degree,'mode':mode,'exit':process.returncode,'cases':sum(line.startswith('T,') for line in process.stdout.splitlines()),'stderr':process.stderr[-3000:]})
    print(degree,mode,'PASS' if process.returncode==0 and not process.stderr else 'FAIL',flush=True)
after={str(p.relative_to(repo)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((repo/'src/quadrature/cutfem_ls').glob('*.hpp'))}
summary={'accepted':all(r['exit']==0 and r['cases']==45 and not r['stderr'] for r in records) and manifest==after,'cases':sum(r['cases'] for r in records),'records':records,'leak_detection':args.leak_detection,'source_changed':manifest!=after,'source_sha256':manifest}
(suite/'review_results/tangency_sanitizers.json').write_text(json.dumps(summary,indent=2)+'\n')
sys.exit(0 if summary['accepted'] else 1)
