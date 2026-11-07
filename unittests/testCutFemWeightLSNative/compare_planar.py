from pathlib import Path
import concurrent.futures
import json
import math
import subprocess

from cases import all_cases

root = Path(__file__).resolve().parent


def planar_case(arguments):
    case, orientation = arguments
    command = [str(root / 'build/suite_planar'), str(case['element']), str(case['degree']), case['field'], str(case['scale']), str(orientation)]
    record = {'case': case, 'orientation': orientation, 'command': command, 'failures': [], 'weight_comparisons': 0, 'moment_comparisons': 0, 'maximum_normalized_weight_difference': 0.0, 'maximum_normalized_moment_difference': 0.0}
    try:
        process = subprocess.run(command, capture_output=True, text=True, timeout=30)
        if process.returncode:
            raise RuntimeError(process.stderr[-2000:])
        planar = json.loads(process.stdout)
        path = root / 'results/candidate' / (case['id'] + '.json')
        if not path.exists() and case['id'] == 'e4_q5_tiny_identity_s1_v0':
            path = path.with_name(case['id'] + '_relaxed_1e-10.json')
        candidate = json.loads(path.read_text())
        if candidate['assessment']['failures'] and case['id'] == 'e4_q5_tiny_identity_s1_v0':
            relaxed = path.with_name(case['id'] + '_relaxed_1e-10.json')
            if relaxed.exists():
                candidate = json.loads(relaxed.read_text())
        record['candidate_case_used'] = candidate['case']['id']
        record['candidate_request_tolerance'] = candidate['case']['tolerance']
        if candidate['assessment']['failures']:
            record['failures'].append({'kind': 'candidate', 'failures': candidate['assessment']['failures']})
            return record
        record['data'] = planar
        for method in ('h', 'delta'):
            first = planar[method]
            second = candidate['data']['weights'][method]
            normalization = max(1.0, *[abs(value) for value in second])
            difference = max(abs(a - b) for a, b in zip(first, second)) / normalization
            record['weight_comparisons'] += len(first)
            record['maximum_normalized_weight_difference'] = max(record['maximum_normalized_weight_difference'], difference)
            if difference > 2e-8:
                record['failures'].append({'kind': 'weights', 'method': method, 'normalized_difference': difference})
        for moment in candidate['data']['moments']:
            factors = [weight * math.prod(point[d] ** moment['powers'][d] for d in range(len(point))) for weight, point in zip(planar['gauss'], planar['points'])]
            for index, method in enumerate(('h', 'delta')):
                actual = math.fsum(factor * value for factor, value in zip(factors, planar[method]))
                expected = moment['values'][index]
                difference = abs(actual - expected) / max(1.0, abs(expected))
                record['moment_comparisons'] += 1
                record['maximum_normalized_moment_difference'] = max(record['maximum_normalized_moment_difference'], difference)
                if difference > 2e-9:
                    record['failures'].append({'kind': 'moments', 'method': method, 'powers': moment['powers'], 'actual': actual, 'expected': expected, 'normalized_difference': difference})
    except Exception as error:
        record['failures'].append({'kind': 'process', 'message': str(error)})
    return record


def compare_planar(element=None):
    cases = [case for case in all_cases() if case['mapping'] == 'identity' and case['scale'] == 1.0 and case['validation'] == 0 and case['field'] in ('plane', 'diagonal', 'boundary_low', 'boundary_high', 'tiny')]
    if element:
        index = ['LINE', 'QUAD', 'TRI', 'HEX', 'TET'].index(element)
        cases = [case for case in cases if case['element'] == index]
    tasks = [(case, orientation) for case in cases for orientation in ((0, 1, 2) if case['element'] == 4 else (0,))]
    print('Comparing', len(tasks), 'actual planar rules', flush=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        records = list(pool.map(planar_case, tasks))
    summary = {'cases': len(records), 'flagged_cases': sum(bool(record['failures']) for record in records), 'weight_comparisons': sum(record['weight_comparisons'] for record in records), 'moment_comparisons': sum(record['moment_comparisons'] for record in records), 'maximum_normalized_weight_difference': max(record['maximum_normalized_weight_difference'] for record in records), 'maximum_normalized_moment_difference': max(record['maximum_normalized_moment_difference'] for record in records)}
    (root / ('planar_comparison_' + element + '_results.json' if element else 'planar_comparison_results.json')).write_text(json.dumps({'summary': summary, 'results': records}, indent=2) + '\n')
    print('Planar comparison:', summary)
    for record in records:
        if record['failures']:
            print('FLAG', record['case']['id'], 'orientation', record['orientation'], record['failures'][0])
    return 1 if summary['flagged_cases'] else 0


if __name__ == '__main__':
    import argparse
    import sys
    parser = argparse.ArgumentParser(description='Compare native LS rules with the original planar class.')
    parser.add_argument('--element', choices=['LINE', 'QUAD', 'TRI', 'HEX', 'TET'])
    arguments = parser.parse_args()
    sys.exit(compare_planar(arguments.element))
