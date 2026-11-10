from pathlib import Path
import argparse
import csv
import hashlib
import json
import math
import subprocess
import sys

from tangency_strip_reference import reference, femus_reference

suite = Path(__file__).resolve().parent
repo = suite.parents[1]
parser = argparse.ArgumentParser(description='Check the retained QUAD/TRI endpoint suite against independently converged polynomial references.')
parser.add_argument('--quick', action='store_true')
parser.add_argument('--skip-build', action='store_true')
parser.add_argument('--rounds', type=int, default=1)
parser.add_argument('--measure', action='store_true')
args = parser.parse_args()
if args.rounds < 1:
    parser.error('At least one round is required.')
folder = suite / 'results' / 'endpoint_strips'
folder.mkdir(parents=True, exist_ok=True)
summary_path = suite / 'endpoint_only_results' / 'endpoint_strip_comparison.json'
summary_path.parent.mkdir(parents=True, exist_ok=True)
manifest = {str(p.relative_to(repo)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((repo / 'src/quadrature/cutfem_ls').glob('*.hpp'))}
if not args.skip_build:
    subprocess.run([sys.executable, str(suite / 'build.py'), '--tests', 'endpoint_strip_benchmark'], check=True)
records = []
failures = []
checks = 0
graph_cases = 0
maximum_ratio = 0
for degree in ([2] if args.quick else [2, 5, 10]):
    for mode in (['surface'] if args.quick else ['raw', 'surface', 'physical']):
        command = [str(suite / 'build' / 'endpoint_strip_benchmark'), str(degree), mode, str(args.rounds), '1e-10', str(int(args.measure))]
        process = subprocess.run(command, capture_output=True, text=True, timeout=600)
        (folder / ('candidate_q' + str(degree) + '_' + mode + '.csv')).write_text(process.stdout)
        if process.returncode:
            failures.append({'degree': degree, 'mode': mode, 'exit': process.returncode, 'error': process.stderr[-3000:]})
        timings = {}
        moments = {}
        coefficients = {}
        native = {}
        for row in csv.reader(process.stdout.splitlines()):
            if row[0] == 'T':
                key = int(row[1]), int(row[4])
                timings.setdefault(key, []).append(row)
            elif row[0] == 'M':
                key = int(row[1]), int(row[4])
                power = int(row[5]), int(row[6])
                if power in moments.setdefault(key, {}):
                    raise RuntimeError('Duplicated reference moment.')
                moments[key][power] = row[7], row[8]
            elif row[0] == 'F':
                native.setdefault((int(row[1]), int(row[4])), {})[int(row[5]), int(row[6])] = row[7], row[8]
            elif row[0] == 'C':
                coefficients[int(row[1]), int(row[2])] = tuple(row[3:])
        if len(timings) != 45:
            failures.append({'degree': degree, 'mode': mode, 'missing_cases': len(timings)})
        for (element, field), rows in timings.items():
            record = {'element': element, 'degree': degree, 'field': field, 'mode': mode, 'name': rows[0][5], 'stats': {'attempts': sum(int(row[9]) for row in rows), 'completed': sum(int(row[10]) for row in rows), 'rejected': sum(int(row[11]) for row in rows)}, 'maximum_absolute_error': [0, 0], 'failures': []}
            if len(rows) != args.rounds or any(row[8] != 'PASS' for row in rows):
                record['failures'].append('Integration failed or a round is missing.')
            actual = moments.get((element, field), {})
            if len(actual) != math.comb(degree+2, 2):
                record['failures'].append('Missing moments.')
            if field < 10:
                graph_cases += 1
            if actual:
                targets, discrepancy = reference(element, coefficients[element, field], degree)
                record['reference_discrepancy'] = discrepancy
                for power, values in actual.items():
                    selected = targets[power][[0, {'raw': 1, 'surface': 2, 'physical': 3}[mode]]]
                    for component, (value, target) in enumerate(zip(values, selected)):
                        error = abs(float(value)-float(target))
                        limit = 1e-10*(1+abs(float(target)))+8*discrepancy+3e-12
                        record['maximum_absolute_error'][component] = max(record['maximum_absolute_error'][component], error)
                        maximum_ratio = max(maximum_ratio, error/limit)
                        checks += 1
                        if not math.isfinite(error) or error > limit:
                            record['failures'].append({'power': power, 'component': component, 'error': error, 'limit': limit})
                if element == 0:
                    targets, discrepancy = femus_reference(element, coefficients[element, field], degree)
                    record['femus_maximum_absolute_error'] = [0, 0]
                    for (i, j) in actual:
                        selected = targets[i, j][[0, {'raw': 1, 'surface': 2, 'physical': 3}[mode]]]
                        for component, target in enumerate(selected):
                            if len(native.get((element, field), {})) != len(actual):
                                raise RuntimeError('Missing native FEMuS reference-coordinate moments.')
                            value = float(native[element, field][i,j][component])
                            error = float(abs(value-target))
                            limit = 1e-10*(1+abs(float(target)))+8*discrepancy+3e-12
                            record['femus_maximum_absolute_error'][component] = max(record['femus_maximum_absolute_error'][component], error)
                            maximum_ratio = max(maximum_ratio, error/limit)
                            checks += 1
                            if not math.isfinite(error) or error > limit:
                                record['failures'].append({'femus_power': [i,j], 'component': component, 'error': error, 'limit': limit})
            records.append(record)
            if record['failures']:
                failures.append(record)
        print(degree, mode, 'cases', len(timings), 'failures', len(failures), flush=True)
after = {str(p.relative_to(repo)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((repo / 'src/quadrature/cutfem_ls').glob('*.hpp'))}
summary = {'accepted': not failures and manifest == after, 'maximum_requested_degree': 2 if args.quick else 10, 'cases': len(records), 'checks': checks, 'graph_cases': graph_cases, 'maximum_error_ratio': maximum_ratio, 'failures': failures, 'records': records, 'source_changed': manifest != after, 'source_sha256': manifest, 'reference': 'Independent polynomial restrictions, 80-digit topology, successively doubled Gauss rules, analytical cases, and independent FEMuS coordinate conversion; discrepancy is retained in each record. Timings against accepted and rifat are in run_endpoint_only_comparison.py.'}
summary_path.write_text(json.dumps(summary, indent=2)+'\n')
print(json.dumps({k:v for k,v in summary.items() if k not in ['records', 'source_sha256', 'failures']}, indent=2))
if failures:
    print(json.dumps(failures[:3], indent=2))
sys.exit(0 if summary['accepted'] else 1)
