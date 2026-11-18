import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import statistics
import subprocess

parser = argparse.ArgumentParser(description='Compare eager and deferred ex47 Bernstein restriction using identical benchmark calls.')
parser.add_argument('--eager', type=Path, required=True)
parser.add_argument('--lazy', type=Path, required=True)
parser.add_argument('--output', type=Path, required=True)
parser.add_argument('--eager-micro', type=Path)
parser.add_argument('--lazy-micro', type=Path)
parser.add_argument('--micro-cpu', type=int, default=0)
parser.add_argument('--weight-cpu', type=int, default=0)
parser.add_argument('--samples', type=int, default=5)
parser.add_argument('--collect-only', action='store_true')
args = parser.parse_args()
args.output.mkdir(parents=True, exist_ok=True)
if not args.collect_only:
    for sample in range(args.samples):
        for backend in (['eager','lazy'] if sample%2 == 0 else ['lazy','eager']):
            with (args.output / (backend+'_micro_'+str(sample)+'.log')).open('w') as stream:
                subprocess.run(['taskset','-c',str(args.micro_cpu),str((getattr(args,backend+'_micro') or getattr(args,backend)).resolve()),'--benchmark-bernstein'],stdout=stream,stderr=subprocess.STDOUT,check=True)
    for dimension in [1,2,3]:
        for backend in (['eager','lazy'] if dimension%2 else ['lazy','eager']):
            with (args.output / (backend+'_weights_'+str(dimension)+'.log')).open('w') as stream:
                subprocess.run(['taskset','-c',str(args.weight_cpu),str(getattr(args,backend).resolve()),'--benchmark-weights','--weight-dimension',str(dimension)],stdout=stream,stderr=subprocess.STDOUT,check=True)

micro = {}
weights = {}
for backend in ['eager','lazy']:
    for sample in range(args.samples):
        for row in csv.reader((args.output / (backend+'_micro_'+str(sample)+'.log')).open()):
            if row and row[0] == 'bernstein_trial':
                micro.setdefault(tuple(row[1:3]),{}).setdefault(backend,[]).append((float(row[3]),row[4]))
    for dimension in [1,2,3]:
        text = (args.output / (backend+'_weights_'+str(dimension)+'.log')).read_text()
        for row in re.findall(r'weight_speed (\w+) q=(\d+) mode=(\d+) old_us=(\S+) prepared_us=(\S+)',text):
            weights.setdefault(tuple(row[:3]),{})[backend] = (float(row[4]),float(row[3]))

micro_rows = []
for key, values in micro.items():
    assert all(len(values[backend]) == args.samples for backend in ['eager','lazy'])
    assert len({row[1] for samples in values.values() for row in samples}) == 1
    row = dict(element=key[0],operation=key[1],checksum=values['eager'][0][1])
    for backend,samples in values.items():
        times = [item[0] for item in samples]
        row.update({backend+'_median_ns':statistics.median(times),backend+'_min_ns':min(times),backend+'_max_ns':max(times)})
    row['eager_over_lazy'] = row['eager_median_ns']/row['lazy_median_ns']
    micro_rows.append(row)
weight_rows = []
for key,values in weights.items():
    assert len(values) == 2
    weight_rows.append(dict(element=key[0],degree=int(key[1]),mode=int(key[2]),eager_us=values['eager'][0],lazy_us=values['lazy'][0],eager_over_lazy=values['eager'][0]/values['lazy'][0],production_eager_us=values['eager'][1],production_lazy_us=values['lazy'][1],production_ratio=values['eager'][1]/values['lazy'][1]))
assert len(micro_rows) == 22 and len(weight_rows) == 58
for name,rows in [('micro.csv',micro_rows),('weights.csv',weight_rows)]:
    with (args.output/name).open('w',newline='') as stream:
        writer = csv.DictWriter(stream,fieldnames=list(rows[0]),lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
summary = dict(baseline_commit='4ec74de8922cc52d57437c1c8bc890e426b581d9',micro_samples=args.samples,micro_checksum_matches=True,micro_rows=len(micro_rows),weight_rows=len(weight_rows),weights_median_ratio=statistics.median(row['eager_over_lazy'] for row in weight_rows),by_element={})
for element in ['LINE3','QUAD9','TRI7','HEX27','TET15','WEDGE21']:
    ratios = [row['eager_over_lazy'] for row in weight_rows if row['element'] == element]
    summary['by_element'][element] = dict(median_ratio=statistics.median(ratios),minimum_ratio=min(ratios),maximum_ratio=max(ratios))
summary['binary_hashes'] = {backend:hashlib.sha256(getattr(args,backend).read_bytes()).hexdigest() for backend in ['eager','lazy']}
summary['micro_binary_hashes'] = {backend:hashlib.sha256((getattr(args,backend+'_micro') or getattr(args,backend)).read_bytes()).hexdigest() for backend in ['eager','lazy']}
summary['micro_cpu'] = args.micro_cpu
summary['weight_cpu'] = args.weight_cpu
summary['method'] = 'micro: five interleaved process runs, 20000 restrictions per row after warmup; weights: existing five-batch median benchmark, two alternating retained elements, final weights cleared; each benchmark pinned to its recorded CPU'
(args.output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
