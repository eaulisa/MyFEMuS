from pathlib import Path
import argparse
import csv
import json
import math
import statistics
import sys


def parse(path):
    timing = {}
    stress = {}
    moments = {}
    for row in csv.reader(Path(path).read_text().splitlines()):
        if not row:
            continue
        if row[0] == 'TIMING':
            key = (row[1], row[2], int(row[3]))
            if key in timing:
                raise ValueError(f'Duplicate timing case {key}')
            timing[key] = {'ok': row[4] == 'OK', 'seconds': float(row[5]) if row[4] == 'OK' else None}
        elif row[0] == 'STRESS':
            key = (row[1], 'stress_' + row[2], 4)
            if key in stress:
                raise ValueError(f'Duplicate stress case {key}')
            stress[key] = {'ok': row[3] == 'OK'}
        elif row[0] == 'MOMENT':
            key = (row[1], row[2], int(row[3]))
            exponent = (int(row[4]), int(row[5]))
            if exponent in moments.setdefault(key, {}):
                raise ValueError(f'Duplicate moment {key} {exponent}')
            values = (float(row[6]), float(row[7]))
            if not all(math.isfinite(value) for value in values):
                raise ValueError(f'Nonfinite moment {key} {exponent}')
            moments[key][exponent] = values
    for key, result in {**timing, **stress}.items():
        if not result['ok']:
            raise ValueError(f'Failed integration {key}')
        degree = key[2]
        expected = {(i, j) for i in range(degree + 1) for j in range(degree - i + 1)}
        if set(moments.get(key, {})) != expected:
            raise ValueError(f'Missing or surplus moments {key}')
        if key in timing and (not math.isfinite(result['seconds']) or result['seconds'] <= 0):
            raise ValueError(f'Invalid timing {key}')
    if set(moments) != set(timing) | set(stress):
        raise ValueError('Orphan moment records')
    return timing, stress, moments


def compare_moments(reference, candidate, absolute_tolerance, relative_tolerance):
    if set(reference) != set(candidate):
        raise ValueError('The case sets differ')
    checks = 0
    worst_fraction = 0.0
    for key in reference:
        if set(reference[key]) != set(candidate[key]):
            raise ValueError(f'The moment sets differ for {key}')
        for exponent in reference[key]:
            for pair in range(2):
                a = reference[key][exponent][pair]
                b = candidate[key][exponent][pair]
                budget = 2 * absolute_tolerance + relative_tolerance * (abs(a) + abs(b))
                fraction = abs(a - b) / budget
                worst_fraction = max(worst_fraction, fraction)
                checks += 1
                if fraction > 1.0:
                    raise ValueError(f'Moment mismatch {key} exponent={exponent} pair={pair}: {a:.18g} vs {b:.18g}; budget={budget:.3g}')
    return checks, worst_fraction


def main():
    parser = argparse.ArgumentParser(description='Strict per-moment correctness gate and paired timing report.')
    parser.add_argument('baseline_a')
    parser.add_argument('candidate_a')
    parser.add_argument('candidate_b')
    parser.add_argument('baseline_b')
    parser.add_argument('--expected-timing', type=int, default=72)
    parser.add_argument('--expected-stress', type=int, default=240)
    parser.add_argument('--absolute-tolerance', type=float, default=1e-10)
    parser.add_argument('--relative-tolerance', type=float, default=1e-10)
    parser.add_argument('--json', type=Path)
    args = parser.parse_args()
    runs = [parse(path) for path in [args.baseline_a, args.candidate_a, args.candidate_b, args.baseline_b]]
    for timing, stress, moments in runs:
        if len(timing) != args.expected_timing or len(stress) != args.expected_stress:
            raise ValueError(f'Incomplete run: timing={len(timing)}, stress={len(stress)}')
    checks = 0
    worst_fraction = 0.0
    for run in runs[1:]:
        number, fraction = compare_moments(runs[0][2], run[2], args.absolute_tolerance, args.relative_tolerance)
        checks += number
        worst_fraction = max(worst_fraction, fraction)
    rows = []
    print('element,case,degree,baseline_us,candidate_us,speedup,accuracy_status')
    for key in sorted(runs[0][0]):
        baseline = statistics.median([runs[0][0][key]['seconds'], runs[3][0][key]['seconds']])
        candidate = statistics.median([runs[1][0][key]['seconds'], runs[2][0][key]['seconds']])
        row = {'element': key[0], 'case': key[1], 'degree': key[2], 'baseline_us': baseline * 1e6, 'candidate_us': candidate * 1e6, 'speedup': baseline / candidate}
        rows.append(row)
        print(f'{key[0]},{key[1]},{key[2]},{baseline*1e6:.6f},{candidate*1e6:.6f},{baseline/candidate:.6f},PASS')
    summary = {}
    for element in sorted({row['element'] for row in rows}):
        values = [row['speedup'] for row in rows if row['element'] == element]
        summary[element] = {'geometric_mean_speedup': math.exp(statistics.mean(math.log(v) for v in values)), 'median_speedup': statistics.median(values), 'cases': len(values)}
        print(element, json.dumps(summary[element], sort_keys=True))
    print(f'PASS per_moment_comparisons={checks} maximum_difference_fraction={worst_fraction:.6g} failures=0')
    if args.json:
        args.json.write_text(json.dumps({'accepted': True, 'per_moment_comparisons': checks, 'maximum_difference_fraction': worst_fraction, 'summary': summary, 'timings': rows}, indent=2) + '\n')
    return 0


if __name__ == '__main__':
    try:
        sys.exit(main())
    except (ValueError, OSError, IndexError, KeyError) as error:
        print('FAIL:', error, file=sys.stderr)
        sys.exit(1)
