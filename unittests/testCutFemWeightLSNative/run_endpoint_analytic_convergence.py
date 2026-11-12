from pathlib import Path
import json,itertools,math,subprocess,mpmath as mp
mp.mp.dps=60
import argparse
suite=Path(__file__).resolve().parent
parser=argparse.ArgumentParser(description='Check degree-10 closed-interface moments against 60-digit ball and sphere formulas.')
parser.add_argument('--binary',type=Path,default=suite/'build/comparison/comparison_candidate')
parser.add_argument('--output',type=Path,default=suite/'endpoint_only_results/analytic_convergence.json')
args=parser.parse_args()
out=args.output
out.parent.mkdir(parents=True,exist_ok=True)
rows=[]
def powers(dimension,degree):
 if dimension==1:return [(i,) for i in range(degree+1)]
 if dimension==2:return [(i,j) for i in range(degree+1) for j in range(degree-i+1)]
 return [(i,j,k) for i in range(degree+1) for j in range(degree-i+1) for k in range(degree-i-j+1)]
def exact(power,radius,center,surface):
 result=mp.mpf(0);dim=len(power)
 for terms in itertools.product(*(range(n+1) for n in power)):
  if any(k%2 for k in terms):continue
  n=sum(terms)
  value=radius**(dim+n)*mp.fprod(mp.gamma(mp.mpf(k+1)/2) for k in terms)/mp.gamma(mp.mpf(dim+n)/2+1)
  if surface:value *= (dim+n)/radius
  result += value*mp.fprod(math.comb(i,k)*center**(i-k) for i,k in zip(power,terms))
 return result
for e in range(5):
 dim=1 if e==0 else (2 if e<3 else 3);simplex=e in (2,4);r=mp.mpf('0.125' if simplex else '0.6');c=mp.mpf('0.25' if simplex else '0');terms=powers(dim,10)
 references=[[exact(p,r,c,surface) for p in terms] for surface in (False,True)]
 for tol in ['1e-8','1e-10','1e-12']:
  cmd=[str(args.binary),str(e),'10','2',tol,'1']
  p=subprocess.run(cmd,capture_output=True,text=True,timeout=120)
  row={'element':e,'degree':10,'tolerance':tol,'exit':p.returncode,'command':cmd}
  if p.returncode:row['error']=p.stderr
  else:
   values=json.loads(p.stdout)['moments'];n=len(terms)
   row.update(volume_max_error=float(max(abs(mp.mpf(str(a))-b) for a,b in zip(values[:n],references[0]))),surface_max_error=float(max(abs(mp.mpf(str(a))-b) for a,b in zip(values[n:],references[1]))))
  rows.append(row);out.write_text(json.dumps(rows,indent=2)+'\n');print(row,flush=True)
