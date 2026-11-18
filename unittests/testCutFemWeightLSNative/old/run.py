from pathlib import Path
import concurrent.futures
import json
import math
import os
import subprocess
import sys
import time

from cases import selected_cases
from reference import dimension, outside_degree, reference

root = Path(__file__).resolve().parent
mode = 'candidate'
selection = 'original'
folder = root / 'results' / mode
methods = ('volume', 'raw_delta', 'surface_inside', 'surface_outside', 'physical_inside', 'physical_outside', 'physical_volume')


def assess(case, data):
    expected, reference_error = reference(case)
    failures = []
    approximate = []
    maximum_errors = [0.0] * len(methods)
    assertions = 0
    maximum_asserted_normalized_error = 0.0
    maximum_error_fraction_of_limit = 0.0
    for mode_index, error in enumerate(data['errors']):
        if error:
            failures.append({'kind': 'integration_exception', 'mode': mode_index, 'message': error})
    for index, error in enumerate(data['fe_errors']):
        assertions += 1
        if not math.isfinite(error) or error > 2e-11:
            failures.append({'kind': 'independent_FE', 'component': index, 'error': error})
    assertions += 1
    if data['cache_error'] != 0.0:
        failures.append({'kind': 'cache', 'error': data['cache_error']})
    constant = expected[(0, 0, 0)]
    for moment in data['moments']:
        exponent = tuple(moment['powers'])
        target = expected[exponent]
        for method, value in enumerate(moment['values']):
            if value is None:
                mode_index = 1 if method == 2 else (2 if method == 4 else 0)
                if not data['errors'][mode_index]:
                    failures.append({'kind': 'nonfinite', 'method': methods[method], 'powers': exponent})
                continue
            error = abs(value - target[method])
            maximum_errors[method] = max(maximum_errors[method], error)
            normalized = error / max(1.0, abs(constant[method]), abs(target[method]))
            exact = True
            if method in (3, 5):
                factor_degree = outside_degree(case, physical=method == 5)
                exact = factor_degree is not None and sum(exponent) + factor_degree <= case['degree']
            if method == 6 and (case['mapping'] == 'stretch' or (case['mapping'] == 'shear' and dimension(case['element']) == 1)):
                exact = sum(exponent) + 1 <= case['degree']
            if exact:
                assertions += 1
                limit = 2.0 * case['tolerance'] if case.get('strict_cap_accuracy') else 20.0 * case['tolerance'] + 3e-11
                maximum_asserted_normalized_error = max(maximum_asserted_normalized_error, normalized)
                maximum_error_fraction_of_limit = max(maximum_error_fraction_of_limit, normalized / limit)
                if not math.isfinite(value) or normalized > limit:
                    failures.append({'kind': 'moment', 'method': methods[method], 'powers': exponent, 'actual': value, 'expected': target[method], 'error': error, 'normalized_error': normalized, 'limit': limit})
            elif exponent == (0, 0, 0):
                approximate.append({'method': methods[method], 'actual': value, 'expected': target[method], 'absolute_error': error, 'relative_error': error / max(abs(target[method]), 1e-300)})
    return {'assertions': assertions, 'maximum_asserted_normalized_error': maximum_asserted_normalized_error, 'maximum_error_fraction_of_limit': maximum_error_fraction_of_limit, 'failures': failures, 'maximum_errors': dict(zip(methods, maximum_errors)), 'approximate_formulations': approximate, 'reference_quadrature_discrepancy': reference_error}


def run(case):
    path = folder / (case['id'] + '.json')
    started = time.monotonic()
    command = [str(root / 'build' / ('suite_' + mode)), str(case['element']), str(case['degree']), case['field'], case['mapping'], str(case['scale']), str(case['tolerance']), str(case['validation'])]
    record = {'case': case, 'command': command}
    try:
        result = subprocess.run(command, capture_output=True, text=True, timeout=case['timeout'], env=dict(os.environ, OPENBLAS_NUM_THREADS='1'))
        record['runtime_seconds'] = time.monotonic() - started
        record['exit_code'] = result.returncode
        if result.returncode:
            record['assessment'] = {'assertions': 0, 'failures': [{'kind': 'process', 'message': result.stderr[-4000:], 'stdout': result.stdout[-1000:]}]}
        else:
            record['data'] = json.loads(result.stdout)
            try:
                record['assessment'] = assess(case, record['data'])
            except Exception as error:
                record['assessment'] = {'assertions': 0, 'failures': [{'kind': 'reference', 'message': str(error)}]}
    except subprocess.TimeoutExpired:
        record['runtime_seconds'] = time.monotonic() - started
        record['assessment'] = {'assertions': 0, 'failures': [{'kind': 'timeout', 'seconds': case['timeout']}]}
    except Exception as error:
        record['runtime_seconds'] = time.monotonic() - started
        record['assessment'] = {'assertions': 0, 'failures': [{'kind': 'runner', 'message': str(error)}]}
    path.write_text(json.dumps(record, indent=2) + '\n')
    return record


def portable_cases(selection):
    strict_id = 'e4_q5_tiny_identity_s1_v0'
    original = selected_cases('all')
    for case in original:
        if case['id'] == strict_id:
            case['strict_cap_accuracy'] = True
    strict = next(case for case in original if case['id'] == strict_id)
    if selection == 'stress':
        return [strict]
    if selection in ('regular', 'original'):
        return original
    if selection == 'supplemental':
        cases = [dict(strict, id=strict_id + '_relaxed_1e-10', tolerance=1e-10)]
        for mapping, scales in (('identity', (2.0 ** -20, 2.0 ** 20, -3.0)), ('affine', (1.0, 2.0 ** -20, 2.0 ** 20))):
            for scale in scales:
                identifier = f'e4_q5_tiny_{mapping}_s{scale:g}_v0'
                cases.append(dict(strict, id=identifier, mapping=mapping, scale=scale, group='strict_cap_scaling'))
        scaled_plane = next(case for case in original if case['id'] == 'e4_q5_plane_identity_s9.53674e-07_v0')
        cases.append(dict(scaled_plane, id=scaled_plane['id'] + '_tol1e-13', tolerance=1e-13, group='tolerance_regression'))
        return cases
    identifiers = {case['id'] for case in selected_cases(selection)}
    return [case for case in original if case['id'] in identifiers]


def main():
    import argparse
    parser = argparse.ArgumentParser(description='Run independent native-reference integration checks.')
    parser.add_argument('--selection', default='original', choices=['original', 'regular', 'supplemental', 'quick', 'coupled', 'edge', 'scale', 'physical', 'validation', 'closed'])
    parser.add_argument('--stress', action='store_true', help='Run only the strict TET tiny-cut case at its original tolerance.')
    parser.add_argument('--jobs', type=int, default=3)
    parser.add_argument('--element', choices=['LINE', 'QUAD', 'TRI', 'HEX', 'TET'])
    args = parser.parse_args()
    global selection
    selection = 'stress' if args.stress else args.selection
    folder.mkdir(parents=True, exist_ok=True)
    cases = portable_cases(selection)
    if args.element:
        element = ['LINE', 'QUAD', 'TRI', 'HEX', 'TET'].index(args.element)
        cases = [case for case in cases if case['element'] == element]
        selection += '_' + args.element
    if not cases:
        parser.error('The requested selection contains no cases.')
    print('Running', len(cases), mode, selection, 'cases', flush=True)
    completed = 0
    failed = 0
    checks = 0
    summary = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = [pool.submit(run, case) for case in cases]
        for future in concurrent.futures.as_completed(futures):
            record = future.result()
            completed += 1
            assessment = record['assessment']
            checks += assessment['assertions']
            if assessment['failures']:
                failed += 1
                first = assessment['failures'][0]
                print('FLAG', record['case']['id'], json.dumps(first), flush=True)
            if completed % 25 == 0 or completed == len(cases):
                print('Progress', completed, '/', len(cases), 'flagged', failed, 'assertions', checks, flush=True)
            summary.append({'id': record['case']['id'], 'runtime_seconds': record['runtime_seconds'], 'assertions': assessment['assertions'], 'failures': len(assessment['failures'])})
    (root / (mode + '_' + selection + '_summary.json')).write_text(json.dumps({'cases': len(cases), 'flagged_cases': failed, 'assertions': checks, 'results': summary}, indent=2) + '\n')
    print('Completed', mode, selection, 'flagged cases', failed, flush=True)
    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(main())
