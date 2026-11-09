from pathlib import Path
import argparse
import concurrent.futures
import json
import os
import subprocess
import sys
import time

root = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description='Run API, actual FEMuS coefficient, and focused native moment checks.')
parser.add_argument('--jobs', type=int, default=2)
parser.add_argument('--extended', action='store_true', help='Also run degree-8/10 simplex probes, with a five-minute process budget.')
parser.add_argument('--tests', nargs='+', help='Run only named focused tests.')
args = parser.parse_args()
tests = [('line_integration', []), ('test_contract', []), ('test_ls_multiprecision', []), ('test_ls_output_range', []), ('test_reference_coefficients', []), ('tet_standard_moments', []), ('test_tet_slice_regression', []), ('test_root_candidate_regression', []), ('quad_tri_standard_moments', []), ('quad_native_quadrants', []), ('quad_scaled_plane', []), ('simplex_high_degree', []), ('test_hex_native', ['5', '1e-12']), ('test_hex_native', ['10', '1e-12'])]
if args.tests:
    if any(name not in [test[0] for test in tests] for name in args.tests):
        parser.error('Unknown focused test name.')
    tests = [test for test in tests if test[0] in args.tests]
elif not args.extended:
    tests = [test for test in tests if test[0] != 'simplex_high_degree']


def run(test):
    name, arguments = test
    record = {'test': name, 'arguments': arguments}
    started = time.monotonic()
    try:
        process = subprocess.run([str(root / 'build' / name), *arguments], capture_output=True, text=True, timeout=300 if name == 'simplex_high_degree' else 120, env=dict(os.environ, OPENBLAS_NUM_THREADS='1'))
        record.update(exit=process.returncode, stdout=process.stdout, stderr=process.stderr)
    except Exception as error:
        record.update(exit=1, stderr=str(error))
    record['runtime_seconds'] = time.monotonic() - started
    print(name, arguments, 'exit', record['exit'], flush=True)
    if record.get('stdout'):
        print(record['stdout'].rstrip(), flush=True)
    if record.get('stderr'):
        print(record['stderr'].rstrip(), flush=True)
    return record


with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
    records = list(pool.map(run, tests))
(root / ('selected_focused_results.json' if args.tests else ('extended_results.json' if args.extended else 'focused_results.json'))).write_text(json.dumps(records, indent=2) + '\n')
sys.exit(any(record['exit'] for record in records))
