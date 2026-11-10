#!/usr/bin/env python3
import argparse
import csv
import functools
import gzip
import itertools
import json
import math
import os
import platform
import hashlib
from pathlib import Path
import subprocess
import sys
import time
import numpy as np
from wide_unit_reference import cofactor, exponents, full_moment, plane_reference, reference, sphere_grid


def index(dimension, power):
    i, j, k = power
    q = i+j+k
    return q*(q+1)//2+j if dimension == 2 else q*(q+1)*(q+2)//6+(q-i)*(q-i+1)//2+k


@functools.lru_cache(None)
def gaussian(power, covariance):
    if sum(power) == 0:
        return 1.0
    if sum(power) % 2:
        return 0.0
    axis = next(i for i, p in enumerate(power) if p)
    rest = list(power)
    rest[axis] -= 1
    value = 0.0
    for other, count in enumerate(rest):
        if count:
            lower = rest.copy()
            lower[other] -= 1
            value += count*covariance[axis][other]*gaussian(tuple(lower), covariance)
    return value


def ellipsoid_reference(element, field, degree, order):
    dimension = 2 if element < 2 else 3
    radius = 1/128 if field == 'tiny_ball' else (.15 if dimension == 2 else .10)
    center = np.array([radius if field in ('touch_ball','near_tangent') else .29,.23,.19][:dimension])
    if field == 'near_tangent':
        center[0] += 1e-12
    axes = np.array([.18,.07,.05][:dimension]) if field in ('ellipse','rotated_ellipse') else np.full(dimension,radius)
    level = 1 if field in ('ellipse','rotated_ellipse') else radius**2
    rotation = np.eye(dimension)
    if field == 'rotated_ellipse':
        rotation[:2,:2] = [[.8,-.6],[.6,.8]]
    matrix = rotation@np.diag(axes)
    covariance = tuple(tuple(v for v in row) for row in matrix@matrix.T)
    volume = math.pi**(dimension/2)/math.gamma(dimension/2+1)*np.prod(axes)
    unit, weight = sphere_grid(dimension, order)
    point = center[:,None]+matrix@unit
    normal = np.prod(axes)*(rotation@(unit/axes[:,None]))
    metric = np.sqrt(np.sum(normal*normal,axis=0))
    physical = cofactor(point,normal)
    result = {}
    for power in exponents(degree,dimension):
        h = d = 0.0
        for central in itertools.product(*(range(p+1) for p in power[:dimension])):
            total = sum(central)
            coefficient = math.prod(math.comb(p,c)*x**(p-c) for p,c,x in zip(power,central,center))
            central_moment = volume*gaussian(central,covariance)/math.prod(dimension+2*j for j in range(1,total//2+1))
            h += coefficient*central_moment
            d += coefficient*central_moment*(total+dimension)/(2*level)
        monomial = math.prod(point[axis]**power[axis] for axis in range(dimension))
        result[power] = np.array([h,d,np.sum(weight*metric*monomial),np.sum(weight*physical*monomial)])
    return result


def positive_root(a,b,c):
    if a == 0:
        return -c/b
    return -2*c/(b+math.sqrt(b*b-4*a*c))


def graph_reference(element,field,degree,order):
    dimension = 2 if element < 2 else 3
    simplex = element%2 != 0
    offset = .25 if field == 'near_plane' else .1
    a = -1e-8 if field == 'near_plane' else (.65 if dimension == 2 else .2)
    b = -1e-8 if field == 'near_plane' else .1
    node,weight = np.polynomial.legendre.leggauss(order)
    xmax = positive_root(a,1,offset-1) if simplex else 1
    x = xmax*(1+node)/2
    wx = xmax*weight/2
    result = {}
    if dimension == 2:
        y = offset+a*x*x
        upper = 1-x if simplex else np.ones_like(x)
        normal = np.array([-2*a*x,np.ones_like(x)])
        metric = np.sqrt(np.sum(normal*normal,axis=0))
        physical = cofactor(np.array([x,y]),normal)
        for i,j,k in exponents(degree,dimension):
            volume = np.sum(wx*x**i*(upper**(j+1)-y**(j+1))/(j+1))
            factor = wx*x**i*y**j
            result[i,j,k] = np.array([volume,np.sum(factor),np.sum(factor*metric),np.sum(factor*physical)])
        return result
    ymax = np.array([positive_root(b,1,offset+a*t*t+t-1) for t in x]) if simplex else np.ones_like(x)
    x = np.broadcast_to(x[:,None],(order,order))
    y = ymax[:,None]*(1+node[None,:])/2
    w = wx[:,None]*ymax[:,None]*weight[None,:]/2
    z = offset+a*x*x+b*y*y
    upper = 1-x-y if simplex else np.ones_like(x)
    normal = np.array([-2*a*x,-2*b*y,np.ones_like(x)])
    metric = np.sqrt(np.sum(normal*normal,axis=0))
    physical = cofactor(np.array([x,y,z]),normal)
    for i,j,k in exponents(degree,dimension):
        factor = w*x**i*y**j
        volume = np.sum(factor*(upper**(k+1)-z**(k+1))/(k+1))
        projection = factor*z**k
        result[i,j,k] = np.array([volume,np.sum(projection),np.sum(projection*metric),np.sum(projection*physical)])
    return result


@functools.lru_cache(None)
def build_reference(element,field,degree):
    if field in ('graph','parabola','near_plane'):
        low = graph_reference(element,field,degree,64)
        high = graph_reference(element,field,degree,112)
    elif field in ('rotated_ellipse','near_tangent'):
        low = ellipsoid_reference(element,field,degree,80)
        high = ellipsoid_reference(element,field,degree,144)
    elif field == 'tiny_plane':
        dimension = 2 if element < 2 else 3
        cut = 1-1/1048576
        high = {}
        for power in exponents(degree,dimension):
            i = power[0]
            transverse = power[1:dimension]
            if element%2:
                width = 1-cut
                h = sum(math.comb(i,r)*cut**(i-r)*width**(r+sum(transverse)+dimension)*math.factorial(r)*math.prod(math.factorial(p) for p in transverse)/math.factorial(r+sum(transverse)+dimension) for r in range(i+1))
                d = cut**i*width**(sum(transverse)+dimension-1)*math.prod(math.factorial(p) for p in transverse)/math.factorial(sum(transverse)+dimension-1)
            else:
                h = -math.expm1((i+1)*math.log(cut))/(i+1)*math.prod(1/(p+1) for p in transverse)
                d = cut**i*math.prod(1/(p+1) for p in transverse)
            high[power] = np.array([h,d,d,0])
        return high,0.0
    else:
        return reference(element,field,degree)
    discrepancy = max(float(np.max(np.abs(low[p]-high[p]))) for p in high)
    if discrepancy > 1e-10:
        raise RuntimeError(f'Independent reference failed: {element} {field} {discrepancy}')
    return high,discrepancy


def cases(element):
    fields = ['full','empty','plane','oblique','slab','near_plane','graph','ball','ellipse','rotated_ellipse','tiny_ball','touch_ball','near_tangent','tiny_plane']
    for field in fields:
        for metric in range(3):
            if metric == 2 and field in ('full','empty','plane','oblique','slab','tiny_plane'):
                continue
            yield field,1e-10,metric,0,3
    for field in ('graph','ball','rotated_ellipse'):
        for direction in ([1,2] if element < 2 else ([1,6] if element == 2 else [2,4])):
            yield field,1e-10,2,direction,3
    for tolerance in (1e-6,1e-8,1e-12):
        yield 'ball',tolerance,1,0,3


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--baseline',required=True)
    parser.add_argument('--candidate',required=True)
    parser.add_argument('--output',required=True)
    parser.add_argument('--reuse-existing',action='store_true')
    args = parser.parse_args()
    output = Path(args.output)
    output.mkdir(parents=True,exist_ok=True)
    cpu = min(os.sched_getaffinity(0)) if hasattr(os,'sched_getaffinity') else None
    if cpu is not None:
        os.sched_setaffinity(0,{cpu})
    metadata = dict(cpu=cpu,platform=platform.platform(),degree=10,repetitions=3,compiler=subprocess.check_output(['c++','--version'],text=True).splitlines()[0],baseline_sha='f5f81a8a088cb990780f326c4a8692ac98265536',primary_rule='G7/K15 for every interval',headers={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__).resolve().parents[2]/'src/quadrature/cutfem_ls').glob('*.hpp')})
    (output/'metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    rows = []
    archive_lines = []
    for element in range(4):
        plan = list(cases(element))
        data = '\n'.join(' '.join(map(str,case)) for case in plan)+'\n'
        for implementation,executable in [('baseline',args.baseline),('candidate',args.candidate)]:
            print(f'{implementation} element={element} degree=10 cases={len(plan)}',flush=True)
            path = output/f'{implementation}-{element}.jsonl'
            if not args.reuse_existing or not path.exists():
                completed = subprocess.run([executable,str(element),'10'],input=data,text=True,capture_output=True,check=True)
                path.write_text(completed.stdout)
                (output/f'{implementation}-{element}.stderr').write_text(completed.stderr)
            content = path.read_text()
            archive_lines.extend('{"implementation":'+json.dumps(implementation)+',"record":'+line+'}' for line in content.splitlines())
            records = [json.loads(line) for line in content.splitlines()]
            if len(records) != len(plan):
                raise RuntimeError('Missing benchmark result')
            for record in records:
                result = {k:v for k,v in record.items() if k not in ('values','statistics','legacy')}
                result['implementation'] = implementation
                result['total_quadrature_samples'] = record['root_samples']+record['outer_samples']
                result['reference_discrepancy'] = None
                result['max_abs_H'] = result['max_abs_D'] = result['max_rel_H'] = result['max_rel_D'] = None
                result['worst_H'] = result['worst_D'] = None
                result['accepted_accuracy'] = False
                if not record['error']:
                    expected,discrepancy = build_reference(element,record['field'],10)
                    result['reference_discrepancy'] = discrepancy
                    values = record['values']
                    count = len(values)//2
                    maxima = [0.,0.]
                    relatives = [0.,0.]
                    worst = [None,None]
                    passed = True
                    for power in expected:
                        n = index(2 if element < 2 else 3,power)
                        references = [expected[power][0],expected[power][record['metric']+1]]
                        for pair,exact in enumerate(references):
                            error = abs(values[pair*count+n]-exact)
                            if error > maxima[pair]:
                                maxima[pair] = error
                                worst[pair] = str(power)
                            if abs(exact) > 1e-280:
                                relatives[pair] = max(relatives[pair],error/abs(exact))
                            passed &= error <= record['tolerance']*(1+abs(exact))+4*discrepancy+3e-14
                    result.update(max_abs_H=maxima[0],max_abs_D=maxima[1],max_rel_H=relatives[0],max_rel_D=relatives[1],worst_H=worst[0],worst_D=worst[1],accepted_accuracy=bool(passed))
                if implementation == 'baseline' and len(record['statistics']) == 1:
                    result['legacy_2d_gk15_panels'] = record['statistics'][0][1]
                    result['legacy_curve_subdivisions'] = record['statistics'][0][4]
                else:
                    for dimension,stat in enumerate(record['statistics'],2):
                        for name,value in zip(('regular','transformed','emergency7','emergency3','subdivisions','collapses','resolution_limits','samples','precision_samples','accepted_moments','skipped_moments','regular_intervals','one_singular','two_singular','higher_maps'),stat):
                            result[f'dim{dimension}_{name}'] = value
                for name,value in zip(('legacy_gk31_2d','legacy_gauss16_32_2d','legacy_gk31_3d'),record.get('legacy',[])):
                    result[name] = value
                rows.append(result)
                if record['error']:
                    print(f'  FAILED {record["field"]} metric={record["metric"]} dir={record["direction"]}: {record["error"]}',flush=True)
    keys = list(dict.fromkeys(k for row in rows for k in row))
    with (output/'comparison.csv').open('w') as stream:
        writer = csv.DictWriter(stream,fieldnames=keys,lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    summary = {}
    for implementation in ('baseline','candidate'):
        group = [row for row in rows if row['implementation'] == implementation]
        summary[implementation] = dict(cases=len(group),completed=sum(row['error'] is None for row in group),accurate=sum(row['accepted_accuracy'] for row in group),max_abs_H=max((row['max_abs_H'] or 0) for row in group),max_abs_D=max((row['max_abs_D'] or 0) for row in group),failures=[{k:row[k] for k in ('element','field','metric','direction','tolerance','error')} for row in group if row['error'] or not row['accepted_accuracy']])
    (output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    (output/'moment_records.jsonl.gz').write_bytes(gzip.compress(('\n'.join(archive_lines)+'\n').encode(),mtime=0))
    print(json.dumps(summary,indent=2))


if __name__ == '__main__':
    main()
