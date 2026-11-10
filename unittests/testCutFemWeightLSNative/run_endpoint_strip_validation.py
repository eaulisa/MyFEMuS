from pathlib import Path
import argparse
import concurrent.futures
import hashlib
import json
import subprocess
import sys

suite = Path(__file__).resolve().parent
repo = suite.parents[1]
parser = argparse.ArgumentParser(description='Validate and compare regularized endpoint quadrature for QUAD and TRI through degree ten.')
parser.add_argument('--collect-only', action='store_true')
parser.add_argument('--skip-build', action='store_true')
parser.add_argument('--boost')
parser.add_argument('--eigen')
parser.add_argument('--rounds', type=int, default=3)
args = parser.parse_args()


def step(script, *options):
    subprocess.run([sys.executable, str(suite/script), *options], cwd=repo, check=True)


if not args.collect_only:
    if not args.skip_build:
        options = ['--jobs','2','--tests','endpoint_strip_benchmark','wide_unit_validation','wide_unit_mesh','native_suite','test_contract','test_ls_multiprecision','test_ls_output_range','test_order_ls_first','test_order_planar_first','line_integration','quad_tri_standard_moments']
        for name in ['boost','eigen']:
            value = getattr(args,name)
            if value:
                options.extend(['--'+name,value])
        step('build.py', *options)
    step('run_endpoint_strip_comparison.py', '--measure', '--rounds', str(args.rounds))
    for element,name in [(0,'quad'),(1,'tri')]:
        step('run_wide_unit_validation.py', '--element', str(element), '--jobs', '2', '--folder', 'endpoint_wide_'+name, '--summary', 'endpoint_wide_'+name+'.json')
        step('run_wide_unit_mesh.py', '--element', str(element), '--jobs', '2', '--folder', 'endpoint_mesh_'+name, '--summary', 'endpoint_mesh_'+name+'.json')
    step('run_wide_double_output.py', '--2d-only', '--folder', 'endpoint_double', '--summary', 'endpoint_double.json')
    step('run_endpoint_strip_sanitizers.py')
    import run
    run.folder = suite/'results/endpoint_native'
    run.folder.mkdir(parents=True,exist_ok=True)
    cases = [c for c in run.portable_cases('original') if c['element'] in [1,2] and c['degree'] <= 10]
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        records = list(pool.map(run.run,cases))
    native = {'accepted':all(not r['assessment']['failures'] for r in records),'cases':len(records),'checks':sum(r['assessment']['assertions'] for r in records),'failures':[{'case':r['case'],'assessment':r['assessment']} for r in records if r['assessment']['failures']]}
    (suite/'review_results/endpoint_native.json').write_text(json.dumps(native,indent=2)+'\n')
    contracts = []
    for name in ['test_contract','test_ls_multiprecision','test_ls_output_range','line_integration','quad_tri_standard_moments']:
        result = subprocess.run([str(suite/'build'/name)],capture_output=True,text=True,timeout=120)
        contracts.append({'test':name,'exit':result.returncode,'output':(result.stdout+result.stderr)[-2000:]})
    (suite/'review_results/endpoint_contracts.json').write_text(json.dumps({'accepted':all(r['exit'] == 0 for r in contracts),'records':contracts,'syntax_checks':['test_order_ls_first','test_order_planar_first']},indent=2)+'\n')

names = ['endpoint_strip_comparison','endpoint_wide_quad','endpoint_wide_tri','endpoint_mesh_quad','endpoint_mesh_tri','endpoint_double','endpoint_native','endpoint_strip_sanitizers','endpoint_contracts']
reports = {name:json.loads((suite/('endpoint_only_results' if name == 'endpoint_strip_comparison' else 'review_results')/(name+'.json')).read_text()) for name in names}
manifest = {str(p.relative_to(repo)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((repo/'src/quadrature/cutfem_ls').glob('*.hpp'))}
failures = [name for name,r in reports.items() if not r.get('accepted') or r.get('source_changed') or (r.get('source_sha256') and r['source_sha256'] != manifest)]
comparison = reports['endpoint_strip_comparison']
if comparison.get('cases') != 405 or comparison.get('graph_cases') != 180:
    failures.append('comparison coverage')
stats = {name:sum(r['stats'][name] for r in comparison.get('records',[])) for name in ['attempts','completed','rejected']}
if stats['attempts'] == 0 or stats['completed'] != stats['attempts'] or stats['rejected'] != 0:
    failures.append('endpoint completion coverage')
summary = {'accepted':not failures,'failures':failures,'maximum_requested_degree':10,'comparison_cases':comparison.get('cases'),'comparison_checks':comparison.get('checks'),'graph_cases':comparison.get('graph_cases'),'endpoint_stats':stats,'wide_cases':sum(reports[n]['cases'] for n in ['endpoint_wide_quad','endpoint_wide_tri']),'mesh_configurations':sum(reports[n]['mesh_levels'] for n in ['endpoint_mesh_quad','endpoint_mesh_tri']),'total_reported_checks':sum(r.get('checks',0) for r in reports.values()),'reports':names,'source_sha256':manifest}
(suite/'review_results/endpoint_strip_acceptance.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps({k:v for k,v in summary.items() if k != 'source_sha256'},indent=2))
sys.exit(0 if summary['accepted'] else 1)
