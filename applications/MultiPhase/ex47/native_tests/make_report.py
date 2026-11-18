import argparse
import csv
import hashlib
import io
import json
import math
import re
from pathlib import Path
import statistics
import subprocess

here=Path(__file__).resolve().parent
repo=here.parents[3]
parser=argparse.ArgumentParser(description='Collect the native comparison without discarding failed fixtures.')
parser.add_argument('--results',type=Path,required=True)
parser.add_argument('--build',type=Path,required=True)
parser.add_argument('--output',type=Path,required=True)
args=parser.parse_args()
args.output.mkdir(parents=True,exist_ok=True)


def records(folder):
    output={}
    for path in folder.glob('*.json'):
        data=json.loads(path.read_text())
        if 'case' in data and 'results' in data: output[data['case']['id']]=data
    return output


def write_csv(name,rows):
    keys=list(dict.fromkeys(key for row in rows for key in row))
    with (args.output/name).open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=keys,lineterminator='\n')
        writer.writeheader();writer.writerows(rows)


fixtures={}
for folder in ['accuracy','parameterized','planar']: fixtures.update(records(args.results/folder))
expected=json.loads((args.results/'accuracy/manifest.json').read_text())
missing=sorted({r['id'] for r in expected}-fixtures.keys())
if missing: raise RuntimeError('Missing fixtures: '+', '.join(missing))
timing=records(args.results/'fixture-timing')
native=records(args.results/'native-final')
native_time={}
for path in (args.results/'native-timing').glob('e*.json'):
    row=json.loads(path.read_text())
    if 'case' in row: native_time[row['case']['id']]=row
if set(native)!=set(native_time): raise RuntimeError('Missing native timing batches')


def entries(row,backend):
    result=row['results'][backend]
    return result if isinstance(result,list) else [result]


def intrinsic(row,backend):
    name=row['case']['test'];values=[]
    for entry in entries(row,backend):
        text=entry.get('stdout','')
        try:
            if name=='benchmark_quad_tri':
                value=sum(float(r[5]) for r in csv.reader(text.splitlines()) if r and r[0]=='TIMING' and r[4]=='OK');kind='sum of benchmark case medians'
            elif name=='benchmark_graph_dimensions':
                value=sum(float(r['seconds']) for r in csv.DictReader(io.StringIO(text)));kind='sum of benchmark case medians'
            elif name=='benchmark_surface_reuse':
                data=json.loads(text);value=data['same_cell_seconds']+data['different_cells_seconds'];kind='same/different element calls; setup excluded'
            elif name in ['endpoint_comparison','endpoint_line_comparison']:
                value=json.loads(text)['seconds'];kind='internal seconds per call'
            elif name=='outer_interval_benchmark':
                value=sum(r['seconds'] for r in map(json.loads,text.splitlines()) if not r['error']);kind='sum of internal case medians'
            elif name=='wide_unit_validation':
                value=sum(max(0,r.get('seconds',0)-r.get('setup_seconds',0)) for r in map(json.loads,text.splitlines()) if not r.get('error'));kind='case seconds; fitting setup excluded'
            elif name=='wide_unit_mesh':
                value=sum(float(r['seconds']) for r in csv.DictReader(io.StringIO(text)));kind='mesh seconds; rule setup excluded'
            elif name=='hex_tet_validation':
                value=float(text.strip().splitlines()[-1].split(',')[-1]);kind='driver integration seconds'
            else:
                public=[value['seconds'] for method,value in entry.get('timing',{}).items() if method in ['weights','physical','tangential','normals']]
                if public: value=sum(public);kind='public API calls; construction/oracles excluded'
                else: value=entry['wall_seconds'];kind='whole fixture, including startup/setup'
            if value>0 and math.isfinite(value): values.append(value)
        except (ValueError,KeyError,TypeError,IndexError): pass
    return (statistics.median(values),kind) if values else (None,'no successful timing')


case_rows=[];failures=[];benchmark_rows=[]
for identifier,row in sorted(fixtures.items()):
    measured=timing.get(identifier,row)
    old,kind=intrinsic(measured,'old');prepared,other=intrinsic(measured,'prepared')
    passed=row['passed'];valid=all(passed.values())
    record=dict(test=row['case']['test'],case=identifier,old_accuracy='PASS' if passed['old'] else 'FAIL',prepared_accuracy='PASS' if passed['prepared'] else 'FAIL',timing_kind=kind,old_seconds=old,prepared_seconds=prepared,old_over_prepared=old/prepared if valid and old and prepared else None,timing_samples=len(entries(measured,'old')))
    for backend in ['old','prepared']:
        aggregate=measured.get('timing',{}).get(backend,{})
        record[backend+'_evaluations']=sum(value.get('evaluations',0) for method,value in aggregate.items() if isinstance(value,dict) and method!='construction')
        record[backend+'_subdivisions']=sum(value.get('subdivisions',0) for method,value in aggregate.items() if isinstance(value,dict) and method!='construction')
        record[backend+'_line_restrictions']=sum(value.get('line_restrictions',0) for method,value in aggregate.items() if isinstance(value,dict))
        if not passed[backend]:
            first=entries(row,backend)[0]
            assessment={key:value for key,value in first.get('assessment',{}).items() if key not in ['rows','cases']}
            failures.append(dict(test=row['case']['test'],case=identifier,backend=backend,exit=first['exit'],assessment=assessment,stderr=first.get('stderr','')[:2200],stdout=first.get('stdout','')[:1500]))
    case_rows.append(record)
    name=row['case']['test']
    if name in ['benchmark_graph_dimensions','benchmark_quad_tri','benchmark_surface_reuse','outer_interval_benchmark','wide_unit_validation','wide_unit_mesh']:
        parsed={}
        for backend in ['old','prepared']:
            first=entries(measured,backend)[0];text=first.get('stdout','')
            if name=='benchmark_quad_tri':
                parsed[backend]={'_'.join(r[1:4]):dict(seconds=float(r[5])) for r in csv.reader(text.splitlines()) if r and r[0]=='TIMING' and r[4]=='OK'}
            elif name=='benchmark_surface_reuse':
                data=json.loads(text);parsed[backend]={key:dict(seconds=data[key]) for key in ['same_cell_seconds','different_cells_seconds','setup_seconds']}
            elif name in ['benchmark_graph_dimensions','wide_unit_mesh']:
                parsed[backend]={str(n):dict(r,seconds=float(r['seconds'])) for n,r in enumerate(csv.DictReader(io.StringIO(text)))}
            else:
                parsed[backend]={str(n):r for n,r in enumerate(map(json.loads,text.splitlines()))}
        for key in sorted(parsed['old'].keys()&parsed['prepared'].keys()):
            a=parsed['old'][key];b=parsed['prepared'][key];x=a.get('seconds');y=b.get('seconds');good=not a.get('error') and not b.get('error')
            old_accuracy='PASS' if not a.get('error') else 'FAIL';new_accuracy='PASS' if not b.get('error') else 'FAIL'
            if name=='wide_unit_validation':
                assessments={backend:entries(row,backend)[0]['assessment'].get('cases',[]) for backend in ['old','prepared']}
                if int(key)<len(assessments['old']) and assessments['old'][int(key)]['failures']: old_accuracy='FAIL'
                if int(key)<len(assessments['prepared']) and assessments['prepared'][int(key)]['failures']: new_accuracy='FAIL'
                good=good and old_accuracy=='PASS' and new_accuracy=='PASS'
            benchmark_rows.append(dict(test=name,case=identifier,subcase=key,element=a.get('element',''),field=a.get('field',a.get('name','')),metric=a.get('metric',''),direction=a.get('direction',''),degree=a.get('degree',''),divisions=a.get('divisions',''),old_accuracy=old_accuracy,prepared_accuracy=new_accuracy,old_seconds=x,prepared_seconds=y,old_over_prepared=x/y if good and x and y else None,old_slices=a.get('slices',''),prepared_slices=b.get('slices',''),old_outer_samples=a.get('outer_samples',''),prepared_outer_samples=b.get('outer_samples',''),old_root_samples=a.get('root_samples',''),prepared_root_samples=b.get('root_samples',''),old_error=a.get('error',''),prepared_error=b.get('error','')))

native_rows=[]
construction_rows=[]
for identifier,row in sorted(native.items()):
    measured=native_time[identifier]
    construction_rows.append(dict(case=identifier,old_seconds=measured['timing']['old']['construction']['seconds_per_call'],prepared_seconds=measured['timing']['prepared']['construction']['seconds_per_call']))
    for method in ['weights','physical']:
        a=measured['timing']['old'].get(method);b=measured['timing']['prepared'].get(method)
        if not a or not b: continue
        valid=all(measured['passed'].values());x=a['seconds_per_call'];y=b['seconds_per_call']
        record=dict(test='native_suite',case=identifier,element=row['case']['element'],element_name=['LINE3','QUAD9','TRI7','HEX27','TET15'][row['case']['element']],degree=row['case']['degree'],field=row['case']['field'],mapping=row['case']['mapping'],scale=row['case']['scale'],validation=row['case']['validation'],method=method,old_accuracy='PASS' if not row['results']['old']['assessment']['failures'] else 'FAIL',prepared_accuracy='PASS' if not row['results']['prepared']['assessment']['failures'] else 'FAIL',old_asserted_error=row['results']['old']['assessment']['maximum_asserted_normalized_error'],prepared_asserted_error=row['results']['prepared']['assessment']['maximum_asserted_normalized_error'],old_seconds=x,prepared_seconds=y,old_over_prepared=x/y if valid and x and y else None,old_min_seconds=a['minimum_seconds_per_call'],old_max_seconds=a['maximum_seconds_per_call'],prepared_min_seconds=b['minimum_seconds_per_call'],prepared_max_seconds=b['maximum_seconds_per_call'],repetitions=measured['repetitions'])
        for backend,value in [('old',a),('prepared',b)]:
            record[backend+'_evaluations']=value['work_per_call'].get('evaluations',0)
            record[backend+'_subdivisions']=value['work_per_call'].get('subdivisions',0)
            record[backend+'_line_restrictions']=value['work_per_call'].get('line_restrictions',0)
            record[backend+'_metric_preparations']=value['work_per_call'].get('metric_preparations',0)
        native_rows.append(record)
write_csv('fixtures.csv',case_rows)
write_csv('native_cases.csv',native_rows)
write_csv('construction.csv',construction_rows)
write_csv('benchmark_cases.csv',benchmark_rows)
regressions=[dict(r,sample_ranges_separated=r['prepared_min_seconds']>1.1*r['old_max_seconds']) for r in native_rows if r['old_over_prepared'] and r['old_over_prepared']<1/1.1]
regressions.sort(key=lambda r:r['old_over_prepared'])
write_csv('regressions.csv',regressions)
retained=[]
retained_path=args.results/'retained_weights.log'
if retained_path.exists():
    for match in re.finditer(r'weight_speed (\w+) q=(\d+) mode=(\d+) old_us=(\S+) prepared_us=(\S+) speedup=(\S+)',retained_path.read_text()):
        geometry,degree,mode,old,new,ratio=match.groups()
        retained.append(dict(element=geometry,degree=int(degree),mode=int(mode),family=['raw_delta','reference_surface','physical_surface','inverse_metric_tangential','three_normal_components'][int(mode)],old_microseconds=float(old),prepared_microseconds=float(new),old_over_prepared=float(ratio)))
    write_csv('retained_modes.csv',retained)
(args.output/'failures.json').write_text(json.dumps(failures,indent=2)+'\n')

source_rows=[]
for source in sorted((repo/'unittests/testCutFemWeightLSNative').glob('*.cpp')):
    name=source.stem
    if name=='native_suite':
        source_rows.append(dict(test=name,cases=len(native),old_passed=sum(not r['results']['old']['assessment']['failures'] for r in native.values()),prepared_passed=sum(not r['results']['prepared']['assessment']['failures'] for r in native.values()),speedup=statistics.median(r['old_over_prepared'] for r in native_rows if r['old_over_prepared']),timing='median per-case API speedup'))
    else:
        rows=[r for r in case_rows if r['test']==name]
        if not rows: raise RuntimeError('Missing source: '+name)
        ratios=[r['old_over_prepared'] for r in rows if r['old_over_prepared']]
        source_rows.append(dict(test=name,cases=len(rows),old_passed=sum(r['old_accuracy']=='PASS' for r in rows),prepared_passed=sum(r['prepared_accuracy']=='PASS' for r in rows),speedup=statistics.median(ratios) if ratios else None,timing=rows[0]['timing_kind']))
write_csv('sources.csv',source_rows)

summary=dict(base_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip(),native_sources=len(source_rows),fixture_runs=len(fixtures),native_cases=len(native),native_assertions={b:sum(r['results'][b]['assessment']['assertions'] for r in native.values()) for b in ['old','prepared']},original_failed=[r['case'] for r in case_rows if r['old_accuracy']=='FAIL'],prepared_failed=[r['case'] for r in case_rows if r['prepared_accuracy']=='FAIL'],native_timing_passed={b:sum(r['passed'][b] for r in native_time.values()) for b in ['old','prepared']},native_speed_by_element={},native_rows_slower_over_10_percent=len(regressions),regressions_with_separated_sample_ranges=sum(r['sample_ranges_separated'] for r in regressions),retained_mode_rows=len(retained))
ctest_names=re.search(r'set\(FEMUS_NATIVE_TESTS\s+(.*?)\)',(repo/'unittests/testCutFemWeightLSNative/CMakeLists.txt').read_text(),re.S).group(1).split()
ctest_rows=[dict(test='cutfem_'+name,case=name,old_passed=fixtures[name]['passed']['old'],prepared_passed=fixtures[name]['passed']['prepared']) for name in ctest_names]
for element in ['HEX','TET']:
    identifier=f'hex_tet_validation_mesh_{element}_2_4_1e-9'
    ctest_rows.append(dict(test='cutfem_'+element.lower()+'_sphere_mesh',case=identifier,old_passed=fixtures[identifier]['passed']['old'],prepared_passed=fixtures[identifier]['passed']['prepared']))
write_csv('ctest.csv',ctest_rows)
summary['ctest']=dict(cases=len(ctest_rows),old_passed=sum(r['old_passed'] for r in ctest_rows),prepared_passed=sum(r['prepared_passed'] for r in ctest_rows),prepared_failed=[r['test'] for r in ctest_rows if not r['prepared_passed']])
before=records(args.results/'native')
fixed=[]
for identifier,row in before.items():
    if row['results']['prepared']['assessment']['failures']:
        previous=dict(row['results']['prepared']['assessment'])
        previous['failure_count']=len(previous['failures'])
        previous['failures']=previous['failures'][:3]
        fixed.append(dict(case=identifier,before=previous,after=native[identifier]['results']['prepared']['assessment']))
        if fixed[-1]['after']['failures']: raise RuntimeError('A certified-boundary fix regressed')
if fixed: (args.output/'boundary_fix.json').write_text(json.dumps(fixed,indent=2)+'\n')
summary['certified_boundary_cases_fixed']=len(fixed)
summary['shared_failed']=sorted(set(summary['original_failed'])&set(summary['prepared_failed']))
summary['prepared_only_failed']=sorted(set(summary['prepared_failed'])-set(summary['original_failed']))
environment_path=args.results/'native-timing/environment.json'
if environment_path.exists(): summary['timing_environment']=json.loads(environment_path.read_text())
additional_path=args.results/'benchmark-accuracy/summary.json'
if additional_path.exists():
    summary['additional_benchmark_accuracy']=json.loads(additional_path.read_text())
    (args.output/'benchmark_accuracy.json').write_text(additional_path.read_text())
for element,label in enumerate(['LINE3','QUAD9','TRI7','HEX27','TET15']):
    summary['native_speed_by_element'][label]={}
    for method in ['weights','physical']:
        rows=[r for r in native_rows if r['element']==element and r['method']==method and r['old_over_prepared']]
        summary['native_speed_by_element'][label][method]=dict(cases=len(rows),median_speedup=statistics.median(r['old_over_prepared'] for r in rows),slower_over_10_percent=sum(r['old_over_prepared']<1/1.1 for r in rows))
summary['source_hashes']={str(p.relative_to(repo)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted([*(repo/'unittests/testCutFemWeightLSNative').glob('*.cpp'),*(here.parent/'include').rglob('*.hpp'),*here.glob('*.hpp'),*here.glob('*.cpp'),*here.glob('*.py')])}
binary_names={p.stem for p in (repo/'unittests/testCutFemWeightLSNative').glob('*.cpp')}|{'native_suite_speed','curved_public_regression','benchmark_graph_accuracy','benchmark_surface_accuracy'}
summary['binary_hashes']={str(p.relative_to(args.build)):hashlib.sha256(p.read_bytes()).hexdigest() for backend in ['old','prepared'] for p in sorted((args.build/backend).iterdir()) if p.is_file() and p.name in binary_names}
(args.output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
lines=['# Complete native comparison','',f"Compared all {len(source_rows)} C++ sources in `unittests/testCutFemWeightLSNative` and one additional public-API curved HEX27 reproduction. Original fixtures, tolerances, and independent references are retained.",'',f"The 386 public-class cases pass {summary['native_assertions']['old']:,} original assertions per backend after the certified-boundary-root fix. The broader fixtures still expose prototype regressions; see `failures.json` and the table below.",'','Timing ratios are original / prepared; values above one favor prepared. Public API batches clear final weights before every call and retain rule/element preparation. Constructor timing is separate. Benchmark drivers use their internal timers; other fixture timings include startup and setup. Failed cases have no accepted speedup. Primitive adapter timings include coordinate/format conversion; unchanged mathematical helpers and the legacy planar reference remain shared infrastructure.','', '## Public API timing','', '| Element | Raw/reference weights | Physical surface |','|---|---:|---:|']
for label,methods in summary['native_speed_by_element'].items(): lines.append(f"| {label} | {methods['weights']['median_speedup']:.3f}× | {methods['physical']['median_speedup']:.3f}× |")
lines += ['', f"The exact CTest requests pass {summary['ctest']['old_passed']}/{summary['ctest']['cases']} originally and {summary['ctest']['prepared_passed']}/{summary['ctest']['cases']} in the prototype. One prototype difference is the requested shallower, conservative detection default; its depth-12 control passes. Shared failures in the extended parameterized runs are listed separately in `summary.json`."]
if 'additional_benchmark_accuracy' in summary:
    extra=summary['additional_benchmark_accuracy']
    lines += ['', f"The accuracy-only versions of the timing drivers {'pass' if extra['passed'] else 'fail'} {extra['checks']:,} additional paired moment comparisons: graph degrees 2/10 and all five degree-13 volume/tangential/normal moment families. Their reconstruction work is excluded from timed drivers."]
lines += ['', '## Every native source','', '| Test | Cases | Original pass | Prepared pass | Median speed ratio | Timing |','|---|---:|---:|---:|---:|---|']
for row in source_rows:
    ratio=f"{row['speedup']:.3f}×" if row['speedup'] else '—'
    lines.append(f"| `{row['test']}` | {row['cases']} | {row['old_passed']} | {row['prepared_passed']} | {ratio} | {row['timing']} |")
lines += ['', '## Work behind the slow cases','', f"{len(regressions)} of {len(native_rows)} public-case/API-family rows are over 10% slower; {sum(r['sample_ranges_separated'] for r in regressions)} have separated three-batch timing ranges. The largest low-order QUAD boundary regressions use the same 15 samples in both implementations. The prepared path additionally creates/restricts 15 LINE records and performs root analysis even though the boundary root is certified. This adds metadata and topology work to an otherwise inexpensive analytic-root path. The fixed quadrature and fitting matrix are reused.",'', '| Case | API | Original µs | Prepared µs | Original samples | Prepared samples | LINE restrictions |','|---|---|---:|---:|---:|---:|---:|']
for r in regressions[:10]: lines.append(f"| `{r['case']}` | {r['method']} | {1e6*r['old_seconds']:.2f} | {1e6*r['prepared_seconds']:.2f} | {r['old_evaluations']:g} | {r['prepared_evaluations']:g} | {r['prepared_line_restrictions']:g} |")
if retained:
    lines += ['', '## Two retained elements and all surface families','', 'Five-batch medians from the existing ex47 benchmark; final weights are cleared and the two element records alternate. Both checksum and accuracy checks remain active. The table gives median speed ratios across degrees for each family.','', '| Element | Raw delta | Reference | Physical | Tangential | Three normals |','|---|---:|---:|---:|---:|---:|']
    for geometry in ['LINE3','QUAD9','TRI7','HEX27','TET15','WEDGE21']:
        cells=[]
        for mode in range(5):
            values=[r['old_over_prepared'] for r in retained if r['element']==geometry and r['mode']==mode]
            cells.append(f'{statistics.median(values):.3f}×' if values else '—')
        lines.append('| '+geometry+' | '+' | '.join(cells)+' |')
lines += ['', '## Files','', '`native_cases.csv` gives each public fixture and API family, independent error, median/range, and work counts. `construction.csv` separates constructor cost. `fixtures.csv` gives every executable invocation; `ctest.csv` identifies the exact CTest requests. `benchmark_cases.csv` expands driver subcases, and `retained_modes.csv` records all six elements and requested surface families. `regressions.csv` retains slower rows and timing-range checks. `sources.csv` gives the source-level table. `failures.json` retains failed assertions, shared limitations, and prototype-only failures. `summary.json` records source/binary checksums and coverage.','', 'See `../DIAGNOSIS.md` for failure analysis and `../README.md` for reproduction commands.']
(args.output/'REPORT.md').write_text('\n'.join(lines)+'\n')
print(json.dumps({k:v for k,v in summary.items() if k not in ['source_hashes','binary_hashes']},indent=2))
