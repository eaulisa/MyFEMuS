from pathlib import Path
import hashlib
import json
import os
import subprocess
import sys

suite = Path(__file__).resolve().parent
repo = suite.parents[1]
folder = suite/'results/endpoint_strips/sanitizers'
folder.mkdir(parents=True, exist_ok=True)
manifest = {str(p.relative_to(repo)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((repo/'src/quadrature/cutfem_ls').glob('*.hpp'))}
command = json.loads((suite/'build/endpoint_strip_benchmark_command.json').read_text())
command = [p for p in command if p not in ['-O2','-O3','-DNDEBUG','-march=x86-64','-mtune=generic']]
command[1:1] = ['-O1','-g','-fsanitize=address,undefined','-fno-omit-frame-pointer']
command[-1] = str(suite/'build/endpoint_strip_sanitized')
with (folder/'build.log').open('w') as log:
    subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True, timeout=300)
records = []
for degree in [2,5,10]:
    for mode in ['raw','surface','physical']:
        command = [str(suite/'build/endpoint_strip_sanitized'), str(degree), mode, '1', '1e-10', '0']
        process = subprocess.run(command, capture_output=True, text=True, timeout=120, env=dict(os.environ, ASAN_OPTIONS='detect_leaks=0', UBSAN_OPTIONS='halt_on_error=1'))
        (folder/f'q{degree}_{mode}.csv').write_text(process.stdout)
        (folder/f'q{degree}_{mode}.stderr').write_text(process.stderr)
        cases = sum(line.startswith('T,') for line in process.stdout.splitlines())
        records.append({'degree':degree,'mode':mode,'exit':process.returncode,'cases':cases,'stderr':process.stderr[-3000:]})
        print(degree, mode, 'exit', process.returncode, 'cases', cases, flush=True)
after = {str(p.relative_to(repo)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((repo/'src/quadrature/cutfem_ls').glob('*.hpp'))}
summary = {'accepted':manifest == after and all(r['exit'] == 0 and r['cases'] == 45 and not r['stderr'] for r in records),'cases':sum(r['cases'] for r in records),'records':records,'source_sha256':manifest,'source_changed':manifest != after,'leak_detection':'Disabled: LeakSanitizer cannot run under ptrace in this environment.'}
(suite/'review_results/endpoint_strip_sanitizers.json').write_text(json.dumps(summary, indent=2)+'\n')
sys.exit(0 if summary['accepted'] else 1)
