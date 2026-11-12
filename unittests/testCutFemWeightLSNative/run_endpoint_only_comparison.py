from pathlib import Path
import argparse
import hashlib
import json
import math
import os
import platform
import statistics
import subprocess
import time

root = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description='Compare identical degree-2/10 native cases and paired timings against two checkouts.')
parser.add_argument('--accepted-root', type=Path, required=True)
parser.add_argument('--rifat-root', type=Path, required=True)
parser.add_argument('--boost', type=Path, default=os.environ.get('BOOST_INCLUDEDIR'))
parser.add_argument('--cxx', default=os.environ.get('CXX', 'g++'))
parser.add_argument('--binary-directory', type=Path, default=root / 'build' / 'comparison')
parser.add_argument('--skip-build', action='store_true')
parser.add_argument('--batches', type=int, default=7)
parser.add_argument('--maximum-3d-repetitions', type=int, default=5)
parser.add_argument('--tolerance', type=float, default=1e-10)
parser.add_argument('--output', type=Path, default=root / 'endpoint_only_results' / 'comparison.json')
args = parser.parse_args()
if args.batches < 3 or args.maximum_3d_repetitions < 1 or not args.tolerance > 0:
    parser.error('At least three batches and a positive tolerance are required.')
args.binary_directory.mkdir(parents=True, exist_ok=True)
args.output.parent.mkdir(parents=True, exist_ok=True)
cpu = min(os.sched_getaffinity(0)) if hasattr(os, 'sched_getaffinity') else None
if cpu is not None:
    os.sched_setaffinity(0, {cpu})
checkouts = {'candidate': root.parents[1], 'accepted': args.accepted_root.resolve(), 'rifat': args.rifat_root.resolve()}
commands = {}
for name, checkout in checkouts.items():
    command = [args.cxx, '-std=c++14', '-O3', '-DNDEBUG', '-march=x86-64', '-mtune=generic']
    if name == 'candidate':
        command.append('-DFEMUS_ENDPOINT_ONLY')
    command.append('-I' + str(checkout / 'src/quadrature/cutfem_ls'))
    if args.boost:
        command.append('-I' + str(args.boost))
    command += [str(root / 'endpoint_comparison.cpp'), '-o', str(args.binary_directory / ('comparison_' + name))]
    commands[name] = command
    if not args.skip_build:
        subprocess.run(command, check=True)


def integrate(name, element, degree, field, repetitions):
    command = [str(args.binary_directory / ('comparison_' + name)), str(element), str(degree), str(field), str(args.tolerance), str(repetitions)]
    started = time.monotonic()
    try:
        result = subprocess.run(command, capture_output=True, text=True, timeout=120)
        if result.returncode:
            return {'exit': result.returncode, 'error': result.stderr.strip(), 'elapsed': time.monotonic() - started}
        return dict(json.loads(result.stdout), exit=0)
    except subprocess.TimeoutExpired:
        return {'exit': 124, 'error': '120-second process budget exhausted'}


def zeroth_reference(element, field):
    dimension = 1 if element == 0 else (2 if element < 3 else 3)
    simplex = element in (2, 4)
    if field in (0, 4):
        cut = 1 - 1 / 256 if field == 4 else (0.2 if simplex else 0.15)
        if simplex:
            return ((1-cut)**dimension / math.factorial(dimension), (1-cut)**(dimension-1) / math.factorial(dimension-1))
        return ((1-cut)*2**(dimension-1), 2**(dimension-1))
    if field == 3:
        first, second = (0.2, 0.6) if simplex else (-0.4, 0.4)
        if simplex:
            volume = (1-(1-first)**dimension+(1-second)**dimension)/math.factorial(dimension)
            surface = ((1-first)**(dimension-1)+(1-second)**(dimension-1))/math.factorial(dimension-1)
            return volume, surface
        return ((2-second+first)*2**(dimension-1), 2**dimension)
    if field in (2, 5):
        radius = 1/128 if field == 5 else (0.125 if simplex else 0.6)
        volume = math.pi**(dimension/2)*radius**dimension/math.gamma(1+dimension/2)
        return volume, dimension*volume/radius
    return None


report = {'compiler': subprocess.check_output([args.cxx, '--version'], text=True).splitlines()[0], 'platform': platform.platform(), 'logical_cpu': cpu, 'tolerance': args.tolerance, 'batches': args.batches, 'maximum_3d_repetitions': args.maximum_3d_repetitions, 'protocol': 'One untimed integration inside each process; calibrate a common repetition count to 5ms using the fastest checkout, cap at 100000 and at maximum_3d_repetitions in 3D; paired batches alternate checkout order; median per-call wall time excludes JSON output and process startup.', 'commands': commands, 'driver_sha256': hashlib.sha256((root / 'endpoint_comparison.cpp').read_bytes()).hexdigest(), 'checkouts': {}, 'rows': []}
for name, checkout in checkouts.items():
    report['checkouts'][name] = {'base_commit': subprocess.check_output(['git', '-C', str(checkout), 'rev-parse', 'HEAD'], text=True).strip(), 'header_sha256': {path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in sorted((checkout / 'src/quadrature/cutfem_ls').glob('*.hpp'))}}
report['binary_sha256'] = {name: hashlib.sha256((args.binary_directory / ('comparison_' + name)).read_bytes()).hexdigest() for name in checkouts}
for degree in (2, 10):
    for element in range(5):
        for field in range(6):
            row = {'element': element, 'degree': degree, 'field': field, 'records': {}}
            for name in checkouts:
                row['records'][name] = integrate(name, element, degree, field, 1)
            good = all(record['exit'] == 0 for record in row['records'].values())
            if good:
                repetitions = max(1, min(100000, math.ceil(0.005 / min(record['seconds'] for record in row['records'].values()))))
                if element >= 3:
                    repetitions = min(repetitions, args.maximum_3d_repetitions)
                samples = {name: [] for name in checkouts}
                for batch in range(args.batches):
                    order = list(checkouts) if batch % 2 == 0 else list(reversed(checkouts))
                    for name in order:
                        measured = integrate(name, element, degree, field, repetitions)
                        if measured['exit']:
                            row['records'][name]['timing_error'] = measured
                        else:
                            samples[name].append(measured['seconds'])
                for name, record in row['records'].items():
                    record['repetitions'] = repetitions
                    record['samples_seconds'] = samples[name]
                    if len(samples[name]) == args.batches:
                        record['median_seconds'] = statistics.median(samples[name])
                        record['min_seconds'] = min(samples[name])
                        record['max_seconds'] = max(samples[name])
                    reference = zeroth_reference(element, field)
                    count = len(record['moments']) // 2
                    if reference:
                        record['zeroth_volume_error'] = abs(record['moments'][0]-reference[0])
                        record['zeroth_surface_error'] = abs(record['moments'][count]-reference[1])
                for name in ('accepted', 'rifat'):
                    count = len(row['records']['candidate']['moments']) // 2
                    difference = [abs(a-b) for a, b in zip(row['records']['candidate']['moments'], row['records'][name]['moments'])]
                    row[name + '_volume_difference'] = max(difference[:count])
                    row[name + '_surface_difference'] = max(difference[count:])
                    if all('median_seconds' in row['records'][key] for key in ('candidate', name)):
                        row[name + '_speedup'] = row['records'][name]['median_seconds']/row['records']['candidate']['median_seconds']
            report['rows'].append(row)
            args.output.write_text(json.dumps(report, indent=2) + '\n')
            print(element, degree, field, {name: record.get('error', record.get('median_seconds')) for name, record in row['records'].items()}, flush=True)
raise SystemExit(any(record['exit'] or 'timing_error' in record for row in report['rows'] for record in row['records'].values()))
