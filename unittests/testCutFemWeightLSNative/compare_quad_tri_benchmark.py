from pathlib import Path
import math
import statistics
import sys

def parse(path):
    timing = {}
    stress = {}
    for raw in Path(path).read_text().splitlines():
        row = raw.split(',')
        if not row:
            continue
        if row[0] == 'TIMING':
            key = (row[1], row[2], int(row[3]))
            if row[4] == 'OK':
                timing[key] = {'ok': True, 'seconds': float(row[5]), 'calls': int(row[6]), 'checksum': float(row[7]), 'h00': float(row[8]), 'd00': float(row[9]), 'max': float(row[10])}
            else:
                timing[key] = {'ok': False, 'error': ','.join(row[11:])}
        elif row[0] == 'STRESS':
            key = (row[1], int(row[2]))
            if row[3] == 'OK':
                stress[key] = {'ok': True, 'checksum': float(row[4]), 'h00': float(row[5]), 'd00': float(row[6]), 'max': float(row[7])}
            else:
                stress[key] = {'ok': False, 'error': ','.join(row[8:])}
    return timing, stress

def representative(first, second):
    keys = set(first) | set(second)
    output = {}
    for key in keys:
        a = first.get(key)
        b = second.get(key)
        if a is None:
            output[key] = b
        elif b is None:
            output[key] = a
        elif a.get('ok') and b.get('ok'):
            output[key] = dict(a)
            output[key]['seconds'] = statistics.median([a['seconds'], b['seconds']])
        else:
            output[key] = a if not a.get('ok') else b
    return output

def close(a, b, scale):
    return abs(a-b) <= 5.0e-8*max(1.0, scale, abs(a), abs(b))

base_a_t, base_a_s = parse(sys.argv[1])
cand_a_t, cand_a_s = parse(sys.argv[2])
cand_b_t, cand_b_s = parse(sys.argv[3])
base_b_t, base_b_s = parse(sys.argv[4])
base_t = representative(base_a_t, base_b_t)
cand_t = representative(cand_a_t, cand_b_t)

print('=== QUAD/TRI TIMING COMPARISON ===')
print('element,case,degree,rifat_us,rebuild_us,speedup,accuracy_status')
groups = {'QUAD': [], 'TRI': []}
for key in sorted(set(base_t) | set(cand_t)):
    base = base_t.get(key)
    cand = cand_t.get(key)
    element, name, degree = key
    status = 'OK'
    if base is None or cand is None:
        status = 'MISSING'
    elif base.get('ok') != cand.get('ok'):
        status = 'FAILURE_MISMATCH'
    elif not base.get('ok'):
        status = 'BOTH_FAIL'
    else:
        scale = max(base['max'], cand['max'])
        if not all(close(base[field],cand[field],scale) for field in ('checksum','h00','d00','max')):
            status = 'OUTPUT_MISMATCH'
    if base and cand and base.get('ok') and cand.get('ok'):
        speedup = base['seconds']/cand['seconds']
        groups[element].append(speedup)
        print(f'{element},{name},{degree},{1e6*base["seconds"]:.3f},{1e6*cand["seconds"]:.3f},{speedup:.3f},{status}')
    else:
        print(f'{element},{name},{degree},NA,NA,NA,{status}')

print('\n=== TIMING SUMMARY ===')
for element in ('QUAD','TRI'):
    values = groups[element]
    print(f'{element}: cases={len(values)} geometric_mean_speedup={math.exp(sum(math.log(v) for v in values)/len(values)):.3f} median_speedup={statistics.median(values):.3f} min={min(values):.3f} max={max(values):.3f}')

print('\n=== RANDOMIZED STRESS COMPARISON ===')
mismatch = {'QUAD': 0, 'TRI': 0}
both_ok = {'QUAD': 0, 'TRI': 0}
base_fail = {'QUAD': 0, 'TRI': 0}
cand_fail = {'QUAD': 0, 'TRI': 0}
for key in sorted(set(base_a_s) | set(cand_a_s)):
    element, index = key
    base = base_a_s.get(key)
    cand = cand_a_s.get(key)
    if base is None or cand is None:
        mismatch[element] += 1
        continue
    if not base.get('ok'):
        base_fail[element] += 1
    if not cand.get('ok'):
        cand_fail[element] += 1
    if base.get('ok') and cand.get('ok'):
        both_ok[element] += 1
        scale = max(base['max'],cand['max'])
        if not all(close(base[field],cand[field],scale) for field in ('checksum','h00','d00','max')):
            mismatch[element] += 1
    elif base.get('ok') != cand.get('ok'):
        mismatch[element] += 1
for element in ('QUAD','TRI'):
    print(f'{element}: both_ok={both_ok[element]} rifat_fail={base_fail[element]} rebuild_fail={cand_fail[element]} mismatches={mismatch[element]}')

bad_timing = 0
for key in set(base_t) | set(cand_t):
    base = base_t.get(key)
    cand = cand_t.get(key)
    if base is None or cand is None or base.get('ok') != cand.get('ok'):
        bad_timing += 1
    elif base.get('ok'):
        scale = max(base['max'],cand['max'])
        if not all(close(base[field],cand[field],scale) for field in ('checksum','h00','d00','max')):
            bad_timing += 1
print(f'\nvalidation_mismatches={bad_timing + mismatch["QUAD"] + mismatch["TRI"]}')
