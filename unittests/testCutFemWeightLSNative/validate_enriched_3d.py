from pathlib import Path
import json
import math
import subprocess
import sys
import time

import numpy as np

suite = Path(__file__).resolve().parent
folder = suite / 'results/enriched_3d'
folder.mkdir(parents=True, exist_ok=True)
summary_path = suite / 'review_results/enriched_3d_summary.json'
summary_path.write_text(json.dumps({'accepted': False, 'status': 'pending'}) + '\n')


def reference(element, order):
    nodes, weights = np.polynomial.legendre.leggauss(order)
    theta = 0.5 * math.pi * nodes[None, :]
    theta_weights = 0.5 * math.pi * weights[None, :]
    cosine = np.cos(theta)
    if element == 'HEX':
        z = nodes[:, None]
        measure = weights[:, None] * theta_weights
        transverse = 1.0 + 0.8 * z*z
        constant = 0.35 + 0.2 * z*z
        radius = np.sqrt(constant/transverse)
        y = radius * np.sin(theta)
        denominator = transverse - 0.2 * y*y
        x = np.sqrt(constant/denominator) * cosine
        volume_measure = measure * radius * cosine
        delta_measure = measure / np.sqrt(transverse * denominator)
        gradient = np.sqrt((-2*x*denominator)**2 + (-2*y*(transverse-0.2*x*x))**2 + (0.4*z-1.6*z*(x*x+y*y))**2)
        x_powers = [np.ones_like(x)]
        y_powers = [np.ones_like(y)]
        z_powers = [np.ones_like(z)]
        for degree in range(1, 12):
            x_powers.append(x_powers[-1] * x)
            y_powers.append(y_powers[-1] * y)
            z_powers.append(z_powers[-1] * z)
        result = {}
        for i in range(11):
            for j in range(11-i):
                for k in range(11-i-j):
                    if i%2 or j%2 or k%2:
                        result[(i,j,k)] = np.zeros(3)
                        continue
                    factor = y_powers[j] * z_powers[k]
                    h = np.sum(volume_measure * factor * 2*x_powers[i+1]/(i+1))
                    d = np.sum(delta_measure * factor * x_powers[i])
                    a = np.sum(delta_measure * factor * x_powers[i] * gradient)
                    result[(i,j,k)] = np.array([h,d,a])
        return result
    constant = 1.0 / 512.0
    z_roots = np.roots([-1,3,-3,1,-27*constant])
    z_roots = sorted(root.real for root in z_roots if abs(root.imag) < 1e-12 and 0 < root.real < 1)
    assert len(z_roots) == 2
    parameter = 0.5*(1+nodes)
    span = z_roots[1]-z_roots[0]
    z = (z_roots[0]+span*np.sin(0.5*math.pi*parameter)**2)[:,None]
    measure = (0.25*math.pi*span*weights*np.sin(math.pi*parameter))[:,None] * theta_weights
    remaining = 1-z[:,0]
    roots = []
    for coordinate, length in zip(z[:,0], remaining):
        values = np.roots([1,-2*length,length*length,-4*constant/coordinate])
        assert max(abs(value.imag) for value in values) < 1e-9
        roots.append(sorted(value.real for value in values))
    roots = np.array(roots)
    middle = 0.5*(roots[:,0]+roots[:,1])[:,None]
    half = 0.5*(roots[:,1]-roots[:,0])[:,None]
    y = middle + half*np.sin(theta)
    separation = half*cosine*np.sqrt((roots[:,2,None]-y)/y)
    center = 0.5*(1-y-z)
    plus = center+0.5*separation
    minus = center-0.5*separation
    delta_measure = measure / (z*np.sqrt(y*(roots[:,2,None]-y)))
    volume_measure = measure * half*cosine*separation
    norm_plus = np.sqrt((y*z*(1-2*plus-y-z))**2+(plus*z*(1-plus-2*y-z))**2+(plus*y*(1-plus-y-2*z))**2)
    norm_minus = np.sqrt((y*z*(1-2*minus-y-z))**2+(minus*z*(1-minus-2*y-z))**2+(minus*y*(1-minus-y-2*z))**2)
    plus_powers = [np.ones_like(plus)]
    minus_powers = [np.ones_like(minus)]
    y_powers = [np.ones_like(y)]
    z_powers = [np.ones_like(z)]
    for degree in range(1, 11):
        plus_powers.append(plus_powers[-1]*plus)
        minus_powers.append(minus_powers[-1]*minus)
        y_powers.append(y_powers[-1]*y)
        z_powers.append(z_powers[-1]*z)
    result = {}
    for i in range(11):
        geometric = sum(plus_powers[i-power]*minus_powers[power] for power in range(i+1))/(i+1)
        for j in range(11-i):
            for k in range(11-i-j):
                factor = y_powers[j]*z_powers[k]
                h = np.sum(volume_measure*factor*geometric)
                d = np.sum(delta_measure*factor*(plus_powers[i]+minus_powers[i]))
                a = np.sum(delta_measure*factor*(plus_powers[i]*norm_plus+minus_powers[i]*norm_minus))
                result[(i,j,k)] = np.array([h,d,a])
    return result


records = []
checks = 0
maximum_error = 0.0
reference_discrepancy = 0.0
for element in ('HEX','TET'):
    low = reference(element, 96)
    high = reference(element, 160)
    discrepancy = max(float(np.max(np.abs(low[key]-high[key]))) for key in high)
    reference_discrepancy = max(reference_discrepancy, discrepancy)
    if discrepancy > 5e-11:
        raise RuntimeError(f'{element} independent quadrature did not converge: {discrepancy}')
    for degree in (2,6,10):
        for mode in ('raw','surface','negative_raw','negative_surface'):
            started = time.monotonic()
            command = [str(suite/'build/enriched_3d'), element, str(degree), mode, '1e-10']
            try:
                process = subprocess.run(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=240)
                record = {'element': element, 'degree': degree, 'mode': mode, 'exit': process.returncode, 'seconds': time.monotonic()-started}
                (folder/f'{element}_q{degree}_{mode}.txt').write_text(process.stdout+process.stderr)
                if process.returncode:
                    record['error'] = process.stderr[-2000:]
                else:
                    failed = []
                    seen = set()
                    for line in process.stdout.splitlines():
                        values = line.split()
                        exponent = tuple(int(value) for value in values[:3])
                        if exponent in seen or any(value < 0 for value in exponent) or sum(exponent) > degree:
                            raise RuntimeError('Invalid or duplicated enriched moment.')
                        seen.add(exponent)
                        actual = np.array([float(value) for value in values[3:]])
                        expected = high[exponent][[0,2 if mode.endswith('surface') else 1]].copy()
                        if mode.startswith('negative'):
                            i,j,k = exponent
                            full = math.factorial(i)*math.factorial(j)*math.factorial(k)/math.factorial(i+j+k+3) if element == 'TET' else (0.0 if i%2 or j%2 or k%2 else 8.0/((i+1)*(j+1)*(k+1)))
                            expected[0] = full-expected[0]
                        error = np.abs(actual-expected)
                        maximum_error = max(maximum_error,float(max(error)))
                        checks += 2
                        limit = 4e-10+4e-10*np.abs(expected)
                        if not np.isfinite(actual).all() or (error > limit).any():
                            failed.append({'powers': exponent, 'actual': actual.tolist(), 'expected': expected.tolist()})
                    if len(seen) != math.comb(degree+3,3):
                        raise RuntimeError('The enriched probe omitted requested moments.')
                    if failed:
                        record.update(exit=1, failures=failed[:8])
                records.append(record)
            except subprocess.TimeoutExpired:
                records.append({'element':element,'degree':degree,'mode':mode,'exit':124,'seconds':time.monotonic()-started})
            print(element,degree,mode,'PASS' if records[-1]['exit'] == 0 else 'FAIL',f"{records[-1]['seconds']:.3f}s",flush=True)
            if records[-1]['exit']:
                print(records[-1],flush=True)
summary = {'accepted': len(records) == 24 and checks == 6080 and all(record['exit'] == 0 for record in records), 'cases': len(records), 'checks': checks, 'maximum_requested_degree': 10, 'maximum_absolute_error': maximum_error, 'independent_reference_discrepancy': reference_discrepancy, 'failures': [record for record in records if record['exit']]}
(folder/'records.json').write_text(json.dumps(records,indent=2)+'\n')
summary_path.write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2),flush=True)
sys.exit(0 if summary['accepted'] else 1)
