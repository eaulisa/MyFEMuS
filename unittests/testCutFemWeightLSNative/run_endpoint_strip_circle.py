from pathlib import Path
import hashlib
import json
import subprocess
import mpmath as mp

mp.mp.dps = 80
root = Path(__file__).resolve().parent
folder = root / 'results' / 'endpoint_strip_circle'
folder.mkdir(parents=True, exist_ok=True)
binary = root / 'build' / 'endpoint_strip_circle'
widths = ['0.015', '0.03', '0.04393398282201787134', '0.05', '0.075']
cases = []


def exact_moment(i, j, surface, signed=False):
    cx, cy, radius = mp.mpf('0.29'), mp.mpf('0.23'), mp.mpf('0.15')
    if signed:
        cx, cy, radius = 2*cx-1, 2*cy-1, 2*radius
    value = mp.mpf(0)
    for k in range(0, i+1, 2):
        for l in range(0, j+1, 2):
            angular = 2*mp.beta(mp.mpf(k+1)/2, mp.mpf(l+1)/2)
            centered = radius**(k+l+1)*angular if surface else radius**(k+l+2)*angular/(k+l+2)
            value += mp.binomial(i, k)*mp.binomial(j, l)*cx**(i-k)*cy**(j-l)*centered
    return value


def verify(lines, triangle):
    moments = {}
    for line in lines:
        f = line.split(',')
        if f[0] == 'M':
            moments[int(f[1]), int(f[2])] = [mp.mpf(f[3]), mp.mpf(f[4])]
    unit_error = [mp.mpf(0), mp.mpf(0)]
    reference_error = [mp.mpf(0), mp.mpf(0)]
    ratios = [mp.mpf(0), mp.mpf(0)]
    for (i, j), measured in moments.items():
        for surface in [0, 1]:
            unit_error[surface] = max(unit_error[surface], abs(measured[surface]-exact_moment(i, j, surface)))
            if triangle:
                actual = measured[surface]
            else:
                actual = sum(mp.binomial(i, k)*mp.binomial(j, l)*mp.mpf(2)**(k+l)*(-1)**(i+j-k-l)*moments[k, l][surface] for k in range(i+1) for l in range(j+1))
                actual *= 2 if surface else 4
            expected = exact_moment(i, j, surface, not triangle)
            error = abs(actual-expected)
            reference_error[surface] = max(reference_error[surface], error)
            ratios[surface] = max(ratios[surface], error/(mp.mpf('1e-10')*(1+abs(expected))))
    return {'moment_count': 2*len(moments), 'unit_max_error': [float(v) for v in unit_error], 'reference_max_error': [float(v) for v in reference_error], 'max_tolerance_ratio': [float(v) for v in ratios], 'passed': all(v <= 1 for v in ratios)}


def run(method, width, degree, triangle, repetitions=16):
    name = f'{method}_{triangle}_q{degree}_{width}.csv'
    command = [str(binary), method, width, str(degree), str(triangle), str(repetitions)]
    process = subprocess.run(command, capture_output=True, text=True, timeout=120)
    (folder / name).write_text(process.stdout)
    if process.returncode:
        raise RuntimeError(str(command)+': '+process.stderr)
    lines = process.stdout.splitlines()
    f = lines[0].split(',')
    result = {'method': method, 'element': 'TRI' if triangle else 'QUAD', 'degree': degree, 'delta': float(width), 'median_seconds': float(f[5]), 'cold_seconds': float(f[6]), 'strip_samples': int(f[7]), 'bulk_samples': int(f[8]), 'strip_bisections': int(f[9]), 'bulk_bisections': int(f[10]), 'strip_depth': int(f[11]), 'bulk_depth': int(f[12]), 'minimum_bulk_width': float(f[13]), 'rule15_calls': int(f[14]), 'rule31_calls': int(f[15]), 'geometry_builds': int(f[16]), 'rejected': int(f[17]), 'baseline_samples': int(f[18]), 'command': command, 'output': str((folder/name).relative_to(root)), **verify(lines, triangle)}
    cases.append(result)
    print(result['element'], degree, method, width, f'{result["median_seconds"]*1000:.5f} ms', 'strip/bulk samples', result['strip_samples'], result['bulk_samples'], 'splits', result['strip_bisections'], result['bulk_bisections'], 'pass', result['passed'], flush=True)


for degree in [2, 5, 10]:
    for triangle in [0, 1]:
        run('baseline', widths[2], degree, triangle)
        for width in widths:
            run('mapped', width, degree, triangle)
sources = [root/'endpoint_strip_circle.cpp', root.parents[1]/'src/quadrature/cutfem_ls/EndpointStripQuadrature.hpp']
summary = {'circle': {'center': [0.29, 0.23], 'radius': 0.15}, 'requested_reference_absolute_tolerance': 1e-10, 'requested_reference_relative_tolerance': 1e-10, 'mapped_unit_absolute_tolerance': 1e-12, 'mapped_unit_relative_tolerance': 1e-12, 'timing': 'Median of nine batches after two warm-up batches; 16 calls per batch; nodal constants alternate by 2.910383045673370361328125e-11 to force geometry rebuilds.', 'reference': '80-digit analytic disk and circle moments from beta functions; QUAD moments transformed to [-1,1]^2 with volume factor four and surface factor two.', 'source_sha256': {str(p.relative_to(root.parents[1])):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}, 'cases': cases, 'all_passed': all(c['passed'] for c in cases)}
(root/'review_results/endpoint_strip_circle.json').write_text(json.dumps(summary, indent=2)+'\n')
raise SystemExit(0 if summary['all_passed'] else 1)
