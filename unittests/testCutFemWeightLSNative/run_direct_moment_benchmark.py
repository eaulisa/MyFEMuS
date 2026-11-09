from pathlib import Path
import csv
import argparse
import hashlib
import io
import json
import os
import statistics
import subprocess
import tarfile
import tempfile

suite = Path(__file__).resolve().parent
repo = suite.parents[1]
parser = argparse.ArgumentParser(description='Compare pinned native graph integration implementations.')
parser.add_argument('--baseline', default='bd1ddd478433e40da0e09f57a1f11f5d2ce0f07b')
parser.add_argument('--folder', default='direct_moments')
parser.add_argument('--summary', default='direct_moment_summary.json')
parser.add_argument('--rounds', type=int, default=2)
args = parser.parse_args()
if args.rounds < 1:
    parser.error('At least one complete round is required.')
baseline = args.baseline
folder = suite / 'results' / args.folder
folder.mkdir(parents=True, exist_ok=True)
summary_path = suite / 'review_results' / args.summary
summary_path.write_text(json.dumps({'accepted': False, 'status': 'pending'}) + '\n')
manifest = {str(p.relative_to(repo)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((repo / 'src/quadrature/cutfem_ls').glob('*.hpp'))}
flags = [os.environ.get('CXX', 'g++'), '-std=c++14', '-O3', '-DNDEBUG', '-march=x86-64', '-mtune=generic']
includes = ['-I' + str(repo / p) for p in ('src/quadrature', 'src/enums', 'src/fe')]
includes.append('-I' + os.environ.get('EIGEN3_INCLUDE_DIR', '/usr/include/eigen3'))
if os.environ.get('BOOST_INCLUDEDIR'):
    includes.append('-I' + os.environ['BOOST_INCLUDEDIR'])
objects = [str(suite / 'build' / (name + '.o')) for name in ('Line', 'Quadrilateral', 'Triangle', 'Hexaedron', 'Tetrahedral')]
objects += [str(suite / 'build' / (p.stem + '.o')) for p in sorted((repo / 'src/quadrature').glob('*Points.cpp'))]
cpu = str(min(os.sched_getaffinity(0)))
records = {}
with tempfile.TemporaryDirectory(prefix='femus_direct_moments_') as temporary:
    temporary = Path(temporary)
    archive = subprocess.check_output(['git', 'archive', baseline, 'src/quadrature/cutfem_ls'], cwd=repo)
    with tarfile.open(fileobj=io.BytesIO(archive)) as stream:
        stream.extractall(temporary / 'baseline', filter='data')
    binaries = {}
    for label, root in (('baseline', temporary / 'baseline'), ('candidate', repo)):
        binaries[label] = temporary / label / 'benchmark'
        binaries[label].parent.mkdir(parents=True, exist_ok=True)
        subprocess.run([*flags, '-I' + str(root / 'src/quadrature/cutfem_ls'), *includes, str(suite / 'benchmark_graph_dimensions.cpp'), *objects, '-o', str(binaries[label])], check=True, timeout=300)
    for repeat in range(args.rounds):
        for degree in (2, 5, 10):
            for label in (('baseline', 'candidate') if repeat%2 == 0 else ('candidate', 'baseline')):
                output = subprocess.check_output(['taskset', '-c', cpu, str(binaries[label]), str(degree)], text=True, timeout=90)
                (folder / f'{label}_q{degree}_{repeat}.csv').write_text(output)
                rows = list(csv.DictReader(output.splitlines()))
                if len(rows) != 40 or len({(r['element'], r['case']) for r in rows}) != 40:
                    raise RuntimeError('The graph timing matrix is incomplete or duplicated.')
                if any(int(r['splits']) != 0 for r in rows):
                    raise RuntimeError('A graph timing case required geometric splitting.')
                records[label, degree, repeat] = rows
                print(label, 'degree', degree, 'round', repeat + 1, 'PASS', flush=True)

comparisons = []
for degree in (2, 5, 10):
    for element in ('QUAD', 'TRI', 'HEX', 'TET'):
        timing = {}
        slices = {}
        for label in ('baseline', 'candidate'):
            rows = [r for repeat in range(args.rounds) for r in records[label, degree, repeat] if r['element'] == element]
            timing[label] = statistics.mean(float(r['seconds']) for r in rows)
            slices[label] = statistics.mean(int(r['slices']) for r in rows)
        entry = {'degree': degree, 'element': element, 'baseline_ms': 1000 * timing['baseline'], 'candidate_ms': 1000 * timing['candidate'], 'speedup': timing['baseline'] / timing['candidate'], 'baseline_mean_slices': slices['baseline'], 'candidate_mean_slices': slices['candidate']}
        entry['round_speedups'] = [sum(float(r['seconds']) for r in records['baseline', degree, repeat] if r['element'] == element)/sum(float(r['seconds']) for r in records['candidate', degree, repeat] if r['element'] == element) for repeat in range(args.rounds)]
        comparisons.append(entry)
        print(json.dumps(entry, sort_keys=True), flush=True)
after = {str(p.relative_to(repo)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((repo / 'src/quadrature/cutfem_ls').glob('*.hpp'))}
if after != manifest:
    raise RuntimeError('Production headers changed during benchmarking.')
summary = {'accepted': True, 'baseline_sha': baseline, 'maximum_monomial_degree': 10, 'graphs_per_element': 10, 'degrees': [2, 5, 10], 'complete_runs_per_version': args.rounds, 'cases_per_version_per_run': 120, 'absolute_tolerance': 1e-10, 'relative_tolerance': 1e-10, 'moments': 'positive volume and reference surface', 'cpu': int(cpu), 'compiler_flags': flags, 'comparisons': comparisons, 'source_sha256': manifest}
summary_path.write_text(json.dumps(summary, indent=2) + '\n')
