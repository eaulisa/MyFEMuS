from pathlib import Path
import argparse
import hashlib
import io
import json
import os
import statistics
import subprocess
import sys
import tarfile
import tempfile
import time

suite = Path(__file__).resolve().parent
repo = suite.parents[1]
parser = argparse.ArgumentParser(description='Validate local-graph HEX/TET integration through degree ten with analytic and mesh-refinement checks.')
parser.add_argument('--skip-benchmarks', action='store_true')
parser.add_argument('--skip-sanitizers', action='store_true')
parser.add_argument('--mesh-max', type=int, default=8, choices=[2, 4, 8])
args = parser.parse_args()
folder = suite / 'results' / 'hex_tet'
folder.mkdir(parents=True, exist_ok=True)
records = []
starting_manifest = {str(path.relative_to(repo)): hashlib.sha256(path.read_bytes()).hexdigest() for path in sorted((repo / 'src/quadrature/cutfem_ls').glob('*.hpp'))}
(suite / 'review_results/hex_tet_summary.json').write_text(json.dumps({'accepted': False, 'status': 'pending'}, indent=2) + '\n')


def run(name, command, timeout=180, required=True):
    started = time.monotonic()
    try:
        process = subprocess.run(list(map(str, command)), cwd=repo, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=timeout, env=dict(os.environ, OPENBLAS_NUM_THREADS='1'))
        record = {'name': name, 'exit': process.returncode, 'seconds': time.monotonic()-started, 'output': process.stdout, 'required': required}
    except subprocess.TimeoutExpired as error:
        output = error.stdout or b''
        record = {'name': name, 'exit': 124, 'seconds': time.monotonic()-started, 'output': output.decode() if isinstance(output, bytes) else output, 'required': required}
    records.append(record)
    (folder / (name + '.log')).write_text(record['output'])
    print(name, 'PASS' if record['exit'] == 0 else 'FAIL', f"{record['seconds']:.3f}s", flush=True)
    if record['exit']:
        print(record['output'][-3000:], flush=True)
    return record


for degree in (2, 6, 10):
    run(f'analytic_q{degree}', [suite / 'build/endpoint_surface_regressions', degree], timeout=240)
for element in ('HEX', 'TET'):
    for degree in (2, 6, 10):
        for divisions in (2, 4, 8):
            if divisions <= args.mesh_max:
                run(f'mesh_{element}_q{degree}_n{divisions}', [suite / 'build/hex_tet_validation', 'mesh', element, degree, divisions, '1e-11'], timeout=240)
    for tolerance in ('1e-8', '1e-10', '1e-12'):
        run(f'tolerance_{element}_{tolerance}', [suite / 'build/hex_tet_validation', 'case', element, 6, 'sphere', tolerance], timeout=180)
    for degree in (2, 10):
        run(f'small_{element}_q{degree}', [suite / 'build/hex_tet_validation', 'case', element, degree, 'small_sphere', '1e-11'], timeout=180)

if not args.skip_sanitizers:
    sanitizer = suite / 'build/endpoint_surface_regressions_asan'
    sanitizer_flags = [os.environ.get('CXX', 'g++'), '-std=c++14', '-O1', '-g', '-fsanitize=address,undefined', '-fno-omit-frame-pointer', '-fno-sanitize-recover=all']
    sanitizer_includes = ['-I' + str(repo / path) for path in ('src/quadrature/cutfem_ls', 'src/quadrature', 'src/enums', 'src/fe')]
    sanitizer_includes.append('-I' + os.environ.get('EIGEN3_INCLUDE_DIR', '/usr/include/eigen3'))
    if os.environ.get('BOOST_INCLUDEDIR'):
        sanitizer_includes.append('-I' + os.environ['BOOST_INCLUDEDIR'])
    sanitizer_objects = [str(suite / 'build' / (name + '.o')) for name in ('Line', 'Quadrilateral', 'Triangle', 'Hexaedron', 'Tetrahedral')]
    sanitizer_objects += [str(suite / 'build' / (path.stem + '.o')) for path in sorted((repo / 'src/quadrature').glob('*Points.cpp'))]
    compiled = run('compile_sanitizers', [*sanitizer_flags, *sanitizer_includes, suite / 'endpoint_surface_regressions.cpp', *sanitizer_objects, '-o', sanitizer], timeout=300)
    if compiled['exit'] == 0:
        run('sanitizers', [sanitizer, 2], timeout=360)
    enriched_sanitizer = suite / 'build/enriched_3d_asan'
    compiled_enriched = run('compile_enriched_sanitizers', [*sanitizer_flags, *sanitizer_includes, suite / 'enriched_3d.cpp', *sanitizer_objects, '-o', enriched_sanitizer], timeout=300)
    if compiled_enriched['exit'] == 0:
        run('sanitizers_enriched_surface', [enriched_sanitizer, 'TET', 2, 'surface', '1e-10'], timeout=240)
        run('sanitizers_enriched_negative_surface', [enriched_sanitizer, 'TET', 2, 'negative_surface', '1e-10'], timeout=240)

benchmarks = {}
if not args.skip_benchmarks:
    revisions = {'pre_3d': '9ee8a7f34993af1abbb0202382d001e55ed43ab6', 'rifat': '0ff0669494e603be248d02c44c0a2e9b81f72b4d'}
    flags = [os.environ.get('CXX', 'g++'), '-std=c++14', '-O3', '-DNDEBUG', '-march=x86-64', '-mtune=generic']
    common = ['-I' + str(repo / path) for path in ('src/quadrature', 'src/enums', 'src/fe')]
    common.append('-I' + os.environ.get('EIGEN3_INCLUDE_DIR', '/usr/include/eigen3'))
    if os.environ.get('BOOST_INCLUDEDIR'):
        common.append('-I' + os.environ['BOOST_INCLUDEDIR'])
    objects = [str(suite / 'build' / (name + '.o')) for name in ('Line', 'Quadrilateral', 'Triangle', 'Hexaedron', 'Tetrahedral')]
    objects += [str(suite / 'build' / (path.stem + '.o')) for path in sorted((repo / 'src/quadrature').glob('*Points.cpp'))]
    with tempfile.TemporaryDirectory(prefix='femus_hex_tet_benchmark_') as temporary:
        temporary = Path(temporary)
        binaries = {}
        for label, revision in {'candidate': None, **revisions}.items():
            root = repo
            if revision:
                root = temporary / label
                archive = subprocess.check_output(['git', 'archive', revision, 'src/quadrature/cutfem_ls'], cwd=repo)
                with tarfile.open(fileobj=io.BytesIO(archive)) as stream:
                    stream.extractall(root, filter='data')
            executable = temporary / ('benchmark_' + label)
            compilation = run('compile_' + label, [*flags, '-I' + str(root / 'src/quadrature/cutfem_ls'), *common, suite / 'hex_tet_validation.cpp', *objects, '-o', executable], timeout=300)
            if compilation['exit'] == 0:
                binaries[label] = executable
        timings = {}
        for element in ('HEX', 'TET'):
            for degree in (2, 6, 10):
                for field in ('plane', 'near_plane', 'graph', 'sphere', 'small_sphere', 'tiny'):
                    key = f'{element}_q{degree}_{field}'
                    timings[key] = {}
                    for label, executable in binaries.items():
                        record = run('timing_' + label + '_' + key, [executable, 'case', element, degree, field, '1e-10'], timeout=60, required=label == 'candidate')
                        if record['exit'] == 0:
                            timings[key][label] = float(record['output'].strip().splitlines()[-1].split(',')[-1])
        for label in revisions:
            benchmarks[label] = {}
            for element in ('HEX', 'TET'):
                common_rows = [row for key, row in timings.items() if key.startswith(element + '_') and label in row and 'candidate' in row]
                benchmarks[label][element] = {'common_successes': len(common_rows), 'cases': 18, 'aggregate_speedup': sum(row[label] for row in common_rows)/sum(row['candidate'] for row in common_rows) if common_rows else None, 'median_case_speedup': statistics.median(row[label]/row['candidate'] for row in common_rows) if common_rows else None}
        (folder / 'benchmark_timings.json').write_text(json.dumps(timings, indent=2) + '\n')

ending_manifest = {str(path.relative_to(repo)): hashlib.sha256(path.read_bytes()).hexdigest() for path in sorted((repo / 'src/quadrature/cutfem_ls').glob('*.hpp'))}
if ending_manifest != starting_manifest:
    records.append({'name': 'source_snapshot', 'exit': 1, 'required': True})
required = [row for row in records if row['required']]
summary = {'accepted': all(row['exit'] == 0 for row in required), 'maximum_monomial_degree': 10, 'required_checks': len(required), 'passed': sum(row['exit'] == 0 for row in required), 'failures': [row['name'] for row in required if row['exit']], 'benchmark': benchmarks, 'source_sha256': {str(path.relative_to(repo)): hashlib.sha256(path.read_bytes()).hexdigest() for path in sorted((repo / 'src/quadrature/cutfem_ls').glob('*.hpp'))}, 'integration_scope': 'Standalone native kernels and actual FEMuS FE basis; full MPI/PETSc ex44 execution excluded.'}
(folder / 'records.json').write_text(json.dumps(records, indent=2) + '\n')
(suite / 'review_results/hex_tet_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
print(json.dumps({key: value for key, value in summary.items() if key != 'source_sha256'}, indent=2), flush=True)
sys.exit(0 if summary['accepted'] else 1)
