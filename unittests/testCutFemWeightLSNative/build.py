from pathlib import Path
import argparse
import concurrent.futures
import hashlib
import json
import os
import subprocess
import sys


root = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description='Build standalone native CutFemWeight_LS checks with actual FEMuS FE and quadrature sources.')
parser.add_argument('--repo', type=Path, default=root.parents[1])
parser.add_argument('--boost', type=Path, default=os.environ.get('BOOST_INCLUDEDIR'))
parser.add_argument('--eigen', type=Path, default=os.environ.get('EIGEN3_INCLUDE_DIR', '/usr/include/eigen3'))
parser.add_argument('--cxx', default=os.environ.get('CXX', 'g++'))
parser.add_argument('--jobs', type=int, default=2)
parser.add_argument('--target', choices=['all', 'regression', 'contracts', 'geometry'], default='all')
parser.add_argument('--tests', nargs='+', help='Build only named test translation units.')
args = parser.parse_args()
repo = args.repo.resolve()
build = root / 'build'
build.mkdir(exist_ok=True)
includes = [repo / 'src/quadrature/cutfem_ls', repo / 'src/quadrature', repo / 'src/enums', repo / 'src/fe', args.eigen]
if args.boost:
    includes.append(args.boost)
common = [args.cxx, '-std=c++14', '-O2', *['-I' + str(path) for path in includes]]
starting_manifest = {str(path.relative_to(repo)): hashlib.sha256(path.read_bytes()).hexdigest() for path in sorted((repo / 'src/quadrature/cutfem_ls').glob('*.hpp'))}


def execute(name, command):
    (build / (name + '_command.json')).write_text(json.dumps(command, indent=2) + '\n')
    with (build / (name + '.log')).open('w') as log:
        result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, timeout=300)
    if result.returncode:
        raise RuntimeError(name + ':\n' + (build / (name + '.log')).read_text()[-12000:])


def compile_object(source):
    target = build / (source.stem + '.o')
    if not target.exists() or source.stat().st_mtime > target.stat().st_mtime:
        execute(source.stem, [*common, '-c', str(source), '-o', str(target)])
    return str(target)


fe_sources = [repo / 'src/fe' / (name + '.cpp') for name in ['Line', 'Quadrilateral', 'Triangle', 'Hexaedron', 'Tetrahedral']]
gauss_sources = sorted((repo / 'src/quadrature').glob('*Points.cpp'))
with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
    objects = list(pool.map(compile_object, fe_sources + gauss_sources))

groups = {
    'regression': ['native_suite', 'planar_suite'],
    'contracts': ['test_contract', 'test_ls_multiprecision', 'test_ls_output_range', 'test_order_ls_first', 'test_order_planar_first'],
    'geometry': ['review_regressions', 'line_integration', 'tri_domain_splitting', 'test_reference_coefficients', 'tet_standard_moments', 'test_tet_slice_regression', 'test_root_candidate_regression', 'quad_tri_standard_moments', 'quad_native_quadrants', 'quad_scaled_plane', 'simplex_high_degree', 'test_hex_native'],
}
names = [name for group in groups.values() for name in group] if args.target == 'all' else groups[args.target]
if args.tests:
    known = [name for group in groups.values() for name in group]
    if any(name not in known for name in args.tests):
        parser.error('Unknown test name; see the groups in build.py.')
    names = args.tests


def compile_test(name):
    command = list(common)
    if name == 'planar_suite' or name.startswith('test_order_'):
        command.append('-I' + str(repo / 'src/quadrature/cutfem'))
    command.append(str(root / (name + '.cpp')))
    syntax = name.startswith('test_order_')
    if syntax:
        command.append('-fsyntax-only')
    else:
        target = 'suite_candidate' if name == 'native_suite' else ('suite_planar' if name == 'planar_suite' else name)
        command.extend([*objects, '-o', str(build / target)])
    print('Compiling', name, flush=True)
    execute(name, command)
    if name not in ('planar_suite', 'test_order_ls_first', 'test_order_planar_first'):
        dependencies = subprocess.check_output([*common, '-MM', str(root / (name + '.cpp'))], text=True)
        (build / (name + '.deps')).write_text(dependencies)
        forbidden = ['/cutfem/', '/cutfemOld/', '/algoim/', '/GramSchmidt.hpp', '/CutFem.hpp', '/CutFemWeight.hpp', '/PlanarMoments_LS.hpp']
        for path in forbidden:
            if path in dependencies:
                raise RuntimeError(name + ': forbidden dependency ' + path)
    return {'test': name, 'exit': 0}


results = []
with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
    futures = {pool.submit(compile_test, name): name for name in names}
    for future in concurrent.futures.as_completed(futures):
        try:
            results.append(future.result())
        except Exception as error:
            results.append({'test': futures[future], 'exit': 1, 'message': str(error)})
            print(str(error), flush=True)
manifest = {str(path.relative_to(repo)): hashlib.sha256(path.read_bytes()).hexdigest() for path in sorted((repo / 'src/quadrature/cutfem_ls').glob('*.hpp'))}
if manifest != starting_manifest:
    results.append({'test': 'source_snapshot', 'exit': 1, 'message': 'Implementation headers changed during compilation; rebuild before running tests.'})
(root / ('build_results_selected.json' if args.tests else 'build_results.json')).write_text(json.dumps(results, indent=2) + '\n')
(build / 'source_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
(build / 'test_source_manifest.json').write_text(json.dumps({name + '.cpp': hashlib.sha256((root / (name + '.cpp')).read_bytes()).hexdigest() for name in names}, indent=2) + '\n')
sys.exit(any(result['exit'] for result in results))
