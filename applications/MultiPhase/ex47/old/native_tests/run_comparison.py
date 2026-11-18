import argparse
import ast
import concurrent.futures
import csv
import hashlib
import importlib.util
import io
import itertools
import json
import math
import os
from pathlib import Path
import statistics
import subprocess
import sys
import time
from types import SimpleNamespace

here = Path(__file__).resolve().parent
repo = here.parents[3]
suite = repo / 'unittests/testCutFemWeightLSNative'
sys.path.insert(0, str(suite))
parser = argparse.ArgumentParser(description='Run unchanged native fixtures against original and prepared ex47 integration.')
parser.add_argument('--build', type=Path, required=True)
parser.add_argument('--output', type=Path, required=True)
parser.add_argument('--group', choices=['focused','parameterized','wide','mesh','enriched','planar','all'], default='all')
parser.add_argument('--samples', type=int, default=1)
parser.add_argument('--jobs', type=int, default=1)
parser.add_argument('--resume', action='store_true')
parser.add_argument('--backend', choices=['both','old','prepared'], default='both')
parser.add_argument('--tests', nargs='*')
args = parser.parse_args()
args.output.mkdir(parents=True, exist_ok=True)
backends = ['old','prepared'] if args.backend == 'both' else [args.backend]


def original_functions(filename, names, globals):
    tree = ast.parse((suite / filename).read_text())
    tree.body = [node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name in names]
    exec(compile(tree, str(suite / filename), 'exec'), globals)
    return globals


wide = original_functions('run_wide_unit_validation.py', ['requests','assess'], dict(math=math,args=SimpleNamespace(smoke=False)))
enriched = original_functions('validate_enriched_3d.py', ['reference'], dict(math=math))
outer_spec = importlib.util.spec_from_file_location('original_outer_benchmark',suite/'run_outer_interval_benchmark.py')
outer = importlib.util.module_from_spec(outer_spec)
outer_spec.loader.exec_module(outer)


def item(name, arguments=None, group='focused', timeout=600, input='', **extra):
    return dict(test=name, arguments=list(map(str, arguments or [])), group=group, timeout=timeout, input=input, **extra)


def manifest():
    required = {'endpoint_comparison','endpoint_line_comparison','enriched_3d','hex_tet_validation','native_suite','outer_interval_benchmark','planar_suite','wide_unit_mesh','wide_unit_validation'}
    rows = [item(p.stem) for p in sorted(suite.glob('*.cpp')) if p.stem not in required]
    rows += [item('benchmark_graph_dimensions',[10]),item('single_graph_root',['benchmark']),item('test_hex_native',[10,'1e-12']),item('endpoint_surface_regressions',[2]),item('endpoint_surface_regressions',[10]),item('endpoint_nested_contacts',[10,'1e-10'])]
    rows.append(item('nodal_cut_detection',configuration='depth12',environment={'EX47_NATIVE_DETECTION_DEPTH':'12'}))
    rows.append(item('curved_public_regression'))
    rows.append(item('review_regressions',[suite/'review_results/independent_reference.txt']))
    for degree in [2,10]:
        for element in range(5):
            for field in range(6): rows.append(item('endpoint_comparison',[element,degree,field,'1e-10',1],'parameterized'))
    for element in range(5):
        for tolerance in ['1e-8','1e-12']: rows.append(item('endpoint_comparison',[element,10,2,tolerance,1],'parameterized',oracle='ball'))
    for graph in [0,1]: rows.append(item('endpoint_line_comparison',[graph,1000],'parameterized'))
    for element in range(4):
        requests=list(outer.cases(element))
        packet=''.join(' '.join(map(str,case))+'\n' for case in requests)
        for degree in [2,10]: rows.append(item('outer_interval_benchmark',[element,degree],'parameterized',1800,packet,requests=requests))
    for element in ['HEX','TET']:
        rows.append(item('hex_tet_validation',['mesh',element,2,4,'1e-9'],'mesh',600))
        for degree in [2,6,10]:
            for field in ['plane','near_plane','graph','sphere','small_sphere','tiny']: rows.append(item('hex_tet_validation',['case',element,degree,field,'1e-10'],'parameterized'))
            for divisions in [2,4,8]: rows.append(item('hex_tet_validation',['mesh',element,degree,divisions,'1e-11'],'mesh',1800))
        for tolerance in ['1e-8','1e-12']: rows.append(item('hex_tet_validation',['case',element,6,'sphere',tolerance],'parameterized'))
        for degree in [2,6,10]:
            for mode in ['raw','surface','negative_raw','negative_surface']: rows.append(item('enriched_3d',[element,degree,mode,'1e-10'],'enriched',600))
    for element in range(4):
        for degree in [2,5,10]:
            requests = wide['requests'](element,degree)
            packet = ''.join(' '.join(map(str,case))+'\n' for case in requests)
            rows.append(item('wide_unit_validation',[element,degree],'wide',1800,packet,requests=requests))
            rows.append(item('wide_unit_mesh',[element,degree,64 if element<2 else 8],'mesh',2400))
    from cases import all_cases
    cases = [case for case in all_cases() if case['mapping']=='identity' and case['scale']==1 and case['validation']==0 and case['field'] in ['plane','diagonal','boundary_low','boundary_high','tiny']]
    for case in cases:
        for direction in ([0,1,2] if case['element']==4 else [0]): rows.append(item('planar_suite',[case['element'],case['degree'],case['field'],case['scale'],direction],'planar',60))
    for row in rows:
        row['id'] = row['test'] + ('_'+'_'.join(Path(value).name if '/' in value else value for value in row['arguments']) if row['arguments'] else '')
        if row.get('configuration'): row['id'] += '_'+row['configuration']
    return rows


items = manifest()
all_sources = {p.stem for p in suite.glob('*.cpp')}
assert all_sources - {r['test'] for r in items} == {'native_suite'}
(args.output/'manifest.json').write_text(json.dumps(items, indent=2)+'\n')
items = [row for row in items if (args.group=='all' or args.group==row['group']) and (not args.tests or row['test'] in args.tests)]


def assess(row, text):
    failures = []
    result = dict(checks=0, failures=failures)
    name = row['test']
    if name=='wide_unit_validation':
        import numpy as np
        from wide_unit_reference import reference, full_moment
        wide.update(np=np,reference=reference,full_moment=full_moment)
        data = [json.loads(line) for line in text.splitlines()]
        if len(data)!=len(row['requests']): failures.append(dict(kind='missing_cases',actual=len(data),expected=len(row['requests'])))
        records = [wide['assess'](int(row['arguments'][0]),int(row['arguments'][1]),case,entry) for case,entry in zip(row['requests'],data)]
        result.update(checks=sum(r['checks'] for r in records),cases=records,failures=failures+[r for r in records if r['failures']])
    elif name=='wide_unit_mesh':
        data = list(csv.DictReader(io.StringIO(text)))
        expected = int(math.log2(int(row['arguments'][2])))
        if len(data)!=expected or any(r['status']!='PASS' for r in data): failures.append(dict(kind='mesh',rows=data,expected=expected))
        result.update(checks=4*len(data),rows=data)
    elif name=='outer_interval_benchmark':
        data=[json.loads(line) for line in text.splitlines()]; element,degree=map(int,row['arguments'])
        if len(data)!=len(row['requests']): failures.append(dict(kind='missing_cases',actual=len(data),expected=len(row['requests'])))
        worst=0.
        for entry in data:
            if entry['error']: failures.append(dict(field=entry['field'],metric=entry['metric'],error=entry['error'])); continue
            reference,discrepancy=outer.build_reference(element,entry['field'],degree)
            values=entry['values']; count=len(values)//2
            for power in reference:
                n=outer.index(2 if element<2 else 3,power)
                for pair,expected in enumerate([reference[power][0],reference[power][entry['metric']+1]]):
                    error=abs(values[pair*count+n]-expected); worst=max(worst,error); result['checks']+=1
                    if error>entry['tolerance']*(1+abs(expected))+4*discrepancy+3e-14: failures.append(dict(field=entry['field'],metric=entry['metric'],power=power,pair=pair,error=error))
        result.update(maximum_absolute_error=worst,rows=data)
    elif name=='enriched_3d':
        import numpy as np
        enriched['np'] = np
        element,degree,mode,_ = row['arguments']
        reference = enriched['reference'](element,160)
        seen = set()
        worst = 0.
        for line in text.splitlines():
            values=line.split(); power=tuple(map(int,values[:3])); actual=np.array(list(map(float,values[3:])))
            if power in seen or len(power)!=3 or sum(power)>int(degree): failures.append(dict(kind='power',power=power))
            seen.add(power)
            expected=reference[power][[0,2 if mode.endswith('surface') else 1]].copy()
            if mode.startswith('negative'):
                i,j,k=power
                full=math.factorial(i)*math.factorial(j)*math.factorial(k)/math.factorial(i+j+k+3) if element=='TET' else (0 if i%2 or j%2 or k%2 else 8/((i+1)*(j+1)*(k+1)))
                expected[0]=full-expected[0]
            error=np.abs(actual-expected); limit=4e-10+4e-10*np.abs(expected)
            worst=max(worst,float(max(error))); result['checks']+=2
            if not np.isfinite(actual).all() or (error>limit).any(): failures.append(dict(power=power,actual=actual.tolist(),expected=expected.tolist()))
        if len(seen)!=math.comb(int(degree)+3,3): failures.append(dict(kind='missing_moments'))
        result['maximum_absolute_error']=worst
    elif name=='endpoint_comparison' and row.get('oracle')=='ball':
        import mpmath as mp
        mp.mp.dps=60
        oracle=original_functions('run_endpoint_analytic_convergence.py',['powers','exact'],dict(mp=mp,itertools=itertools,math=math))
        element=int(row['arguments'][0]); dimension=1 if element==0 else (2 if element<3 else 3); simplex=element in [2,4]
        radius=mp.mpf('0.125' if simplex else '0.6'); center=mp.mpf('0.25' if simplex else '0')
        reference=[float(oracle['exact'](power,radius,center,surface)) for surface in [False,True] for power in oracle['powers'](dimension,10)]
        values=json.loads(text)['moments']; tolerance=float(row['arguments'][3])
        for n,(actual,expected) in enumerate(zip(values,reference)):
            if abs(actual-expected)>8*tolerance*(1+abs(expected)): failures.append(dict(moment=n,actual=actual,expected=expected))
        if len(values)!=len(reference): failures.append(dict(kind='missing_moments'))
        result.update(checks=len(values),maximum_absolute_error=max(abs(a-b) for a,b in zip(values,reference)))
    else:
        if any(repr_marker in text for repr_marker in [',FAIL','"status":"FAIL"']): failures.append(dict(kind='reported_failure'))
        result['checks']='embedded'
    return result


def execute(row):
    path=args.output/(row['id']+'.json')
    if args.resume and path.exists():
        previous=json.loads(path.read_text())
        if set(previous['results'])!=set(backends): raise ValueError('Resume backend mismatch: '+row['id'])
        return previous
    result=dict(case=row,results={backend:[] for backend in backends})
    for sample in range(args.samples):
        for backend in (backends if sample%2==0 else list(reversed(backends))):
            start=time.monotonic()
            try:
                process=subprocess.run([str(args.build/backend/row['test']),*row['arguments']],input=row['input'],capture_output=True,text=True,timeout=row['timeout'],env=dict(os.environ,OPENBLAS_NUM_THREADS='1',**row.get('environment',{})))
                elapsed=time.monotonic()-start
                timing={}
                for line in process.stderr.splitlines():
                    if line.startswith('EX47_NATIVE_TIMING '): timing=json.loads(line[len('EX47_NATIVE_TIMING '):])
                entry=dict(exit=process.returncode,wall_seconds=elapsed,timing=timing,stdout=process.stdout,stderr=process.stderr,assessment=assess(row,process.stdout) if process.returncode==0 else dict(failures=[dict(kind='exit',code=process.returncode)]))
            except Exception as error:
                entry=dict(exit=124 if isinstance(error,subprocess.TimeoutExpired) else 1,wall_seconds=time.monotonic()-start,assessment=dict(failures=[dict(kind=type(error).__name__,message=str(error))]))
            result['results'][backend].append(entry)
    result['passed']={backend:all(not e['exit'] and not e['assessment']['failures'] for e in values) for backend,values in result['results'].items()}
    if len(backends)==2 and all(result['passed'].values()) and row['test']=='endpoint_comparison':
        old=json.loads(result['results']['old'][0]['stdout'])['moments']; prepared=json.loads(result['results']['prepared'][0]['stdout'])['moments']
        tolerance=float(row['arguments'][3]); difference=max(abs(a-b)/max(1,abs(a),abs(b)) for a,b in zip(old,prepared))
        result['paired_accuracy']=dict(maximum_scaled_difference=difference,limit=16*tolerance,checks=len(old))
        if len(old)!=len(prepared) or difference>16*tolerance: result['passed']['prepared']=False
    if len(backends)==2 and all(result['passed'].values()) and row['test']=='benchmark_quad_tri':
        functions=original_functions('compare_quad_tri_benchmark.py',['parse','compare_moments'],dict(Path=Path,csv=csv,math=math))
        parsed={}
        for backend in ['old','prepared']:
            output=args.output/(row['id']+'_'+backend+'.csv'); output.write_text(result['results'][backend][0]['stdout']); parsed[backend]=functions['parse'](output)
        try:
            checks,fraction=functions['compare_moments'](parsed['old'][2],parsed['prepared'][2],1e-10,1e-10)
            result['paired_accuracy']=dict(checks=checks,maximum_fraction_of_limit=fraction,timing_cases=len(parsed['old'][0]),stress_cases=len(parsed['old'][1]))
        except ValueError as error:
            result['paired_accuracy']=dict(failure=str(error)); result['passed']['prepared']=False
    result['timing']={}
    for backend,entries in result['results'].items():
        methods=set().union(*(entry.get('timing',{}).keys() for entry in entries))
        result['timing'][backend]={method:{key:statistics.median(e.get('timing',{}).get(method,{}).get(key,0) for e in entries) for key in entries[0].get('timing',{}).get(method,{})} for method in methods}
        result['timing'][backend]['wall_seconds']=statistics.median(e['wall_seconds'] for e in entries)
    path.write_text(json.dumps(result)+'\n')
    print(row['id'],json.dumps(result['passed']),flush=True)
    if not all(result['passed'].values()):
        for backend,entries in result['results'].items():
            for entry in entries:
                if entry['exit'] or entry['assessment']['failures']: print(backend,entry.get('stderr','')[-1200:],str(entry['assessment']['failures'])[:1200],flush=True)
    return result


with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
    records=list(pool.map(execute,items))
summary=dict(cases=len(records),sources=len({r['case']['test'] for r in records}),backends={backend:dict(passed=sum(r['passed'][backend] for r in records),failed=[r['case']['id'] for r in records if not r['passed'][backend]]) for backend in backends})
(args.output/('summary_'+args.group+'.json')).write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary),flush=True)
raise SystemExit(any(not all(r['passed'].values()) for r in records))
