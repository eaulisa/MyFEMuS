from pathlib import Path
import argparse
import csv
import hashlib
import importlib.util
import json
import os
import platform
import subprocess
import sys
import tempfile
import time

suite = Path(__file__).resolve().parent
repo = suite.parents[1]
results = suite / 'review_results'
results.mkdir(exist_ok=True)
compiler = os.environ.get('CXX', 'g++')
flags = ['-std=c++14', '-O3', '-DNDEBUG', '-march=x86-64', '-mtune=generic']


def execute(command, output=None, timeout=600, expected=0):
    print('RUN', ' '.join(map(str, command)), flush=True)
    completed = subprocess.run(list(map(str, command)), cwd=repo, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=timeout)
    if output:
        Path(output).write_text(completed.stdout)
    if completed.returncode != expected:
        print(completed.stdout[-14000:], flush=True)
        raise RuntimeError(f'Command returned {completed.returncode}, expected {expected}: {command}')
    if not output or str(output).endswith('.log'):
        print(completed.stdout[-6000:], flush=True)
    return completed.stdout


def includes(root):
    paths = [root / path for path in ['src/quadrature/cutfem_ls', 'src/quadrature', 'src/enums', 'src/fe']]
    paths.append(Path(os.environ.get('EIGEN3_INCLUDE_DIR', '/usr/include/eigen3')))
    if os.environ.get('BOOST_INCLUDEDIR'):
        paths.append(Path(os.environ['BOOST_INCLUDEDIR']))
    return ['-I' + str(path) for path in paths]


def cpu_description():
    cpu = subprocess.run(['lscpu'], text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if cpu.returncode == 0:
        return cpu.stdout
    details = Path('/proc/cpuinfo')
    if details.exists():
        return next((line for line in details.read_text().splitlines() if line.startswith('model name')), platform.machine())
    return platform.machine()


def make_references():
    import mpmath as mp
    mp.mp.dps = 60
    path = results / 'independent_reference.txt'
    with path.open('w') as stream:
        for shape in [0, 1]:
            for beta in [mp.mpf('0.25'), mp.mpf('-0.25')]:
                stream.write(f'{shape} {beta}\n')
                def root(x):
                    if shape == 0:
                        return (mp.mpf('0.25') - mp.mpf('0.125') * x*x) / (1 + beta*x*x)
                    b = 1 + beta*x*(1-x)
                    discriminant = b*b - beta*x
                    return mp.mpf('0.5') / (b + mp.sqrt(discriminant))
                def derivative(x):
                    return 1 + beta*x*x if shape == 0 else 1 + beta*x*(1-x) - 2*beta*x*root(x)
                for i in range(7):
                    for j in range(7-i):
                        if shape == 0:
                            h = mp.quad(lambda x: x**i * (1-root(x)**(j+1))/(j+1), [-1, 0, 1])
                            d = mp.quad(lambda x: x**i * root(x)**j / derivative(x), [-1, 0, 1])
                        else:
                            h = mp.quad(lambda x: x**i * ((1-x)**(j+1)-root(x)**(j+1))/(j+1), [0, mp.mpf('0.75')])
                            d = mp.quad(lambda x: x**i * root(x)**j / derivative(x), [0, mp.mpf('0.75')])
                        stream.write(f'{i} {j} {mp.nstr(h, 55)} {mp.nstr(d, 55)}\n')
    return path


def test_gate():
    directory = results / 'gate_self_tests'
    directory.mkdir(exist_ok=True)
    lines = ['TIMING,QUAD,probe,2,OK,0.00001,1,0,0,0,0', 'STRESS,TRI,0,OK,0,0,0,0']
    for shape, name, degree in [('QUAD', 'probe', 2), ('TRI', 'stress_0', 4)]:
        for i in range(degree+1):
            for j in range(degree-i+1):
                lines.append(f'MOMENT,{shape},{name},{degree},{i},{j},0.25,0.5')
    good = directory / 'good.csv'
    good.write_text('\n'.join(lines) + '\n')
    variants = {
        'changed_moment': [row.replace('MOMENT,QUAD,probe,2,0,1,0.25,0.5', 'MOMENT,QUAD,probe,2,0,1,0.2501,0.5') for row in lines],
        'missing_moment': [row for row in lines if not row.startswith('MOMENT,QUAD,probe,2,0,1,')],
        'common_failure': [row.replace(',2,OK,', ',2,FAIL,') for row in lines],
        'nonfinite': [row.replace('MOMENT,QUAD,probe,2,0,1,0.25,0.5', 'MOMENT,QUAD,probe,2,0,1,nan,0.5') for row in lines],
    }
    compare = suite / 'compare_quad_tri_benchmark.py'
    execute([sys.executable, compare, good, good, good, good, '--expected-timing', '1', '--expected-stress', '1'], directory / 'valid.log')
    for name, variant in variants.items():
        bad = directory / (name + '.csv')
        bad.write_text('\n'.join(variant) + '\n')
        baseline = bad if name == 'common_failure' else good
        execute([sys.executable, compare, baseline, bad, bad, baseline, '--expected-timing', '1', '--expected-stress', '1'], directory / (name + '.log'), expected=1)
    print('REVIEW checker rejects corrupt, missing, common-failure and nonfinite results', flush=True)


def checkout(sha, name, existing=None):
    if existing:
        location = Path(existing).resolve()
        if subprocess.check_output(['git', '-C', str(location), 'rev-parse', 'HEAD'], text=True).strip() != sha:
            raise RuntimeError(f'The {name} checkout is not at the pinned commit {sha}')
        entries = subprocess.check_output(['git', '-C', str(location), 'ls-tree', '-r', sha, '--', 'src/quadrature/cutfem_ls'], text=True).splitlines()
        expected = set()
        for entry in entries:
            metadata, path = entry.split('\t', 1)
            expected.add(path)
            content = (location / path).read_bytes()
            digest = hashlib.sha1(b'blob ' + str(len(content)).encode() + b'\0' + content).hexdigest()
            if digest != metadata.split()[2]:
                raise RuntimeError(f'The pinned {name} header was modified: {path}')
        if {str(path.relative_to(location)) for path in (location / 'src/quadrature/cutfem_ls').glob('*.hpp')} != expected:
            raise RuntimeError(f'The {name} checkout contains unexpected CutFEM headers')
        return location
    location = Path(tempfile.mkdtemp(prefix='cutfem-review-' + name + '-')) / 'source'
    if subprocess.run(['git', 'cat-file', '-e', sha + '^{commit}'], cwd=repo, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).returncode:
        execute(['git', 'fetch', '--depth=1', 'origin', sha])
    execute(['git', 'worktree', 'add', '--detach', location, sha])
    return location


def main():
    parser = argparse.ArgumentParser(description='Validate degree-10 QUAD/TRI moments and compare pinned revisions.')
    parser.add_argument('--baseline-repo', type=Path)
    parser.add_argument('--reviewed-repo', type=Path)
    args = parser.parse_args()
    (results / 'acceptance.json').write_text(json.dumps({'accepted': False, 'status': 'running'}, indent=2) + '\n')
    start = time.monotonic()
    manifest = {str(path.relative_to(repo)): hashlib.sha256(path.read_bytes()).hexdigest() for path in sorted((repo / 'src/quadrature/cutfem_ls').glob('*.hpp'))}
    (results / 'source_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    environment = {'compiler': execute([compiler, '--version']), 'flags': flags, 'platform': platform.platform(), 'cpu': cpu_description(), 'cpu_affinity': sorted(os.sched_getaffinity(0)), 'input_sha': execute(['git', 'rev-parse', 'HEAD']).strip(), 'baseline_sha': '0ff0669494e603be248d02c44c0a2e9b81f72b4d', 'reviewed_sha': '39bdedb4b93e6e5f7a21ab4a154eaf3e3538fdaa'}
    environment['production_changes'] = subprocess.check_output(['git', 'status', '--porcelain', '--untracked-files=no', '--', 'src/quadrature/cutfem_ls'], cwd=repo, text=True).strip()
    environment['production_manifest_sha256'] = hashlib.sha256(json.dumps(manifest, sort_keys=True).encode()).hexdigest()
    environment['asan_options'] = os.environ.get('ASAN_OPTIONS', 'runtime defaults')
    (results / 'environment.json').write_text(json.dumps(environment, indent=2) + '\n')
    oracle = make_references()
    binary = suite / 'build/review_regressions'
    binary.parent.mkdir(exist_ok=True)
    execute([compiler, *flags, *includes(repo), suite / 'review_regressions.cpp', '-o', binary])
    execute([binary, oracle], results / 'review_regressions.log', timeout=600)
    sanitized = suite / 'build/review_regressions_sanitized'
    execute([compiler, '-std=c++14', '-O1', '-g', '-fsanitize=address,undefined', '-fno-sanitize-recover=all', '-fno-omit-frame-pointer', *includes(repo), suite / 'review_regressions.cpp', '-o', sanitized], results / 'sanitizer_build.log')
    execute([sanitized, oracle], results / 'sanitizer_run.log', timeout=600)
    test_gate()
    mesh = suite / 'build/mesh_refinement'
    execute([compiler, *flags, *includes(repo), suite / 'mesh_refinement.cpp', '-o', mesh], results / 'mesh_refinement_build.log')
    execute(['taskset', '-c', str(environment['cpu_affinity'][0]), mesh], results / 'mesh_refinement.csv', timeout=300)
    mesh_rows = list(csv.DictReader((results / 'mesh_refinement.csv').read_text().splitlines()))
    if len(mesh_rows) != 30 or any(row['status'] != 'PASS' for row in mesh_rows):
        raise RuntimeError('The QUAD/TRI mesh-refinement matrix is incomplete or failed')
    baseline = checkout(environment['baseline_sha'], 'baseline', args.baseline_repo)
    reviewed = checkout(environment['reviewed_sha'], 'reviewed', args.reviewed_repo)
    binaries = {}
    source = suite / 'benchmark_quad_tri.cpp'
    for name, revision in [('baseline', baseline), ('reviewed', reviewed), ('candidate', repo)]:
        target = suite / 'build' / ('benchmark_' + name)
        execute([compiler, *flags, *includes(revision), source, '-o', target], results / ('build_' + name + '.log'))
        binaries[name] = target
    cpu = str(environment['cpu_affinity'][0])
    for name, round_name in [('baseline','a'), ('reviewed','a'), ('candidate','a'), ('candidate','b'), ('reviewed','b'), ('baseline','b')]:
        execute(['taskset', '-c', cpu, binaries[name]], results / f'{name}_{round_name}.csv', timeout=600)
    for reference in ['baseline', 'reviewed']:
        execute([sys.executable, suite / 'compare_quad_tri_benchmark.py', results / f'{reference}_a.csv', results / 'candidate_a.csv', results / 'candidate_b.csv', results / f'{reference}_b.csv', '--json', results / f'comparison_{reference}.json'], results / f'comparison_{reference}.log')
    final_manifest = {str(path.relative_to(repo)): hashlib.sha256(path.read_bytes()).hexdigest() for path in sorted((repo / 'src/quadrature/cutfem_ls').glob('*.hpp'))}
    if manifest != final_manifest:
        raise RuntimeError('Source changed while testing')
    status = {'accepted': True, 'input_sha': environment['input_sha'], 'production_changes': environment['production_changes'], 'production_manifest_sha256': environment['production_manifest_sha256'], 'elapsed_seconds': time.monotonic()-start, 'standalone_regressions': 'PASS', 'asan_ubsan': 'PASS', 'checker_self_tests': 'PASS', 'mesh_refinement': 'PASS', 'mesh_refinement_cases': len(mesh_rows), 'momentwise_baseline_comparison': 'PASS', 'momentwise_reviewed_comparison': 'PASS'}
    (results / 'acceptance.json').write_text(json.dumps(status, indent=2) + '\n')
    print('FINAL_REVIEW_ACCEPTANCE', json.dumps(status, sort_keys=True), flush=True)


if __name__ == '__main__':
    try:
        main()
    except Exception as error:
        (results / 'acceptance.json').write_text(json.dumps({'accepted': False, 'error': str(error)}, indent=2) + '\n')
        raise
