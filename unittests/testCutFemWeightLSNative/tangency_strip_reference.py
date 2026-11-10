import functools
import math

import mpmath as mp
import numpy as np


def polynomial_sum(a,b):
    c=[mp.mpf(0)]*max(len(a),len(b))
    for i,v in enumerate(a):c[i]+=v
    for i,v in enumerate(b):c[i]+=v
    return c


def polynomial_product(a,b):
    c=[mp.mpf(0)]*(len(a)+len(b)-1)
    for i,u in enumerate(a):
        for j,v in enumerate(b):c[i+j]+=u*v
    return c


def real_roots(coefficient):
    c=list(coefficient)
    while len(c)>1 and c[-1]==0:c.pop()
    if len(c)<2:return []
    if len(c)==2:return [-c[0]/c[1]]
    if len(c)==3:
        d=c[1]**2-4*c[2]*c[0]
        if d<0 and abs(d)<mp.mpf('1e-65')*max(1,abs(c[1]**2),abs(4*c[2]*c[0])):d=mp.mpf(0)
        if d<0:return []
        return [(-c[1]-mp.sqrt(d))/(2*c[2]),(-c[1]+mp.sqrt(d))/(2*c[2])]
    try:roots=mp.polyroots(list(reversed(c)),maxsteps=1000,extraprec=160)
    except Exception:
        guesses=np.roots([float(v) for v in reversed(c)])
        roots=[mp.findroot(lambda x:mp.polyval(list(reversed(c)),x),mp.mpf(float(v.real)),maxsteps=200) for v in guesses if abs(v.imag)<1e-7]
    return [mp.re(v) for v in roots if abs(mp.im(v))<mp.mpf('1e-40')]


@functools.lru_cache(maxsize=None)
def geometry(element,coefficients):
    mp.mp.dps=80
    c=[[mp.mpf(coefficients[3*i+j]) for j in range(3)] for i in range(3)]
    a,b,d=[[c[i][j] for i in range(3)] for j in [2,1,0]]
    discriminant=polynomial_sum(polynomial_product(b,b),[-4*v for v in polynomial_product(a,d)])
    upper=[mp.mpf(0)]
    for i in range(3):
        for j in range(3):
            p=[mp.mpf(0)]*i+[c[i][j]]
            for n in range(j):p=polynomial_product(p,[1,-1] if element else [1])
            upper=polynomial_sum(upper,p)
    events=[mp.mpf(0),mp.mpf(1)]
    for p in [a,b,d,discriminant,upper]:events.extend(v for v in real_roots(p) if 0<v<1)
    distinct=[]
    for value in sorted(events):
        if not distinct or abs(value-distinct[-1])>mp.mpf('1e-30'):distinct.append(value)
    return c,a,b,d,discriminant,distinct


def evaluate(c,x):
    result=np.zeros_like(x)
    for v in reversed(c):result=result*x+np.longdouble(str(v))
    return result


def divide_linear(c,root):
    q=[mp.mpf(0)]*(len(c)-1)
    q[-1]=c[-1]
    for k in range(len(q)-1,0,-1):q[k-1]=c[k]+root*q[k]
    return q


def vertical_reference(element,coefficients,degree):
    mp.mp.dps=80
    c=[mp.mpf(coefficients[3*i]) for i in range(3)]
    roots=sorted(set(v for v in real_roots(c) if 0<v<1))
    cuts=[mp.mpf(0),*roots,mp.mpf(1)]
    result={}
    for i in range(degree+1):
        for j in range(degree+1-i):
            h=mp.mpf(0)
            for left,right in zip(cuts[:-1],cuts[1:]):
                if mp.polyval(list(reversed(c)),(left+right)/2)>0:
                    h+=sum((-1)**k*math.comb(j+1,k)*(right**(i+k+1)-left**(i+k+1))/(i+k+1) for k in range(j+2))/(j+1) if element else (right**(i+1)-left**(i+1))/((i+1)*(j+1))
            surface=sum(r**i*((1-r)**(j+1) if element else 1)/(j+1) for r in roots)
            delta=sum(r**i*((1-r)**(j+1) if element else 1)/((j+1)*abs(c[1]+2*c[2]*r)) for r in roots)
            result[i,j]=np.array([float(h),float(delta),float(surface),float(surface)])
    return result


def at_order(element,coefficients,degree,order):
    c,a,b,d,discriminant,events=geometry(element,coefficients)
    nodes,weights=np.polynomial.legendre.leggauss(order)
    t=np.longdouble(math.pi)*np.asarray(nodes,dtype=np.longdouble)/4+np.longdouble(math.pi)/4
    sine,cosine=np.sin(t),np.cos(t)
    results={(i,j):np.zeros(4,dtype=np.longdouble) for i in range(degree+1) for j in range(degree+1-i)}
    for left,right in zip(events[:-1],events[1:]):
        if right-left<mp.mpf('1e-28'):continue
        width=np.longdouble(str(right-left))
        x=np.longdouble(str(left))+width*sine*sine
        jacobian=np.longdouble(math.pi)/2*width*sine*cosine*np.asarray(weights,dtype=np.longdouble)
        top=1-x if element else np.ones_like(x)
        aa,bb,cc=evaluate(a,x),evaluate(b,x),evaluate(d,x)
        reduced=list(discriminant)
        left_degree=right_degree=0
        norm=max(1,sum(abs(v) for v in reduced))
        while len(reduced)>1 and abs(mp.polyval(list(reversed(reduced)),left))<mp.mpf('1e-60')*norm:
            reduced=divide_linear(reduced,left);left_degree+=1
        while len(reduced)>1 and abs(mp.polyval(list(reversed(reduced)),right))<mp.mpf('1e-60')*norm:
            reduced=divide_linear(reduced,right);right_degree+=1
        dd=evaluate(reduced,x)*(width*sine*sine)**left_degree*(-width*cosine*cosine)**right_degree
        middle=(left+right)/2
        am,bm,cm=[mp.polyval(list(reversed(p)),middle) for p in [a,b,d]]
        roots=[]
        midpoint_roots=[]
        if am==0:
            if bm!=0:roots=[-cc/bb];midpoint_roots=[-cm/bm]
        elif bm*bm-4*am*cm>0:
            s=np.sqrt(np.maximum(dd,0))
            numerator=-np.longdouble('0.5')*(bb+np.where(bb>=0,s,-s))
            first=numerator/aa
            second=cc/numerator
            roots=[np.minimum(first,second),np.maximum(first,second)]
            midpoint_roots=sorted([(-bm-mp.sqrt(bm*bm-4*am*cm))/(2*am),(-bm+mp.sqrt(bm*bm-4*am*cm))/(2*am)])
        accepted=[]
        for root,m in zip(roots,midpoint_roots):
            if 0<m<(1-middle if element else 1):accepted.append(root)
        cuts=[np.zeros_like(x),*accepted,top]
        positive=[]
        for lower,upper in zip(cuts[:-1],cuts[1:]):
            sample=(lower+upper)/2
            index=len(sample)//2
            positive.append(bool((aa*sample*sample+bb*sample+cc)[index]>0))
        surface=[]
        physical=[]
        delta=[]
        for root in accepted:
            gx=np.zeros_like(x)
            gy=np.zeros_like(x)
            for i in range(3):
                for j in range(3):
                    value=np.longdouble(str(c[i][j]))
                    if i:gx+=i*value*x**(i-1)*root**j
                    if j:gy+=j*value*x**i*root**(j-1)
            derivative=np.abs(bb) if am==0 else np.sqrt(np.maximum(dd,0))
            delta.append(jacobian/derivative)
            surface.append(jacobian*np.sqrt(gx*gx+gy*gy)/derivative)
            physical.append(jacobian*np.sqrt((gx-.24*x*gy)**2+((1+.2*x)*gy)**2)/derivative)
        for (i,j),result in results.items():
            xp=x**i
            for lower,upper,inside in zip(cuts[:-1],cuts[1:],positive):
                if inside:
                    geometric=sum(upper**(j-k)*lower**k for k in range(j+1))/(j+1)
                    result[0]+=np.sum(jacobian*xp*(upper-lower)*geometric)
            for r,root in enumerate(accepted):
                monomial=xp*root**j
                result[1]+=np.sum(delta[r]*monomial)
                result[2]+=np.sum(surface[r]*monomial)
                result[3]+=np.sum(physical[r]*monomial)
    return results


@functools.lru_cache(maxsize=None)
def reference(element,coefficients,degree):
    if all(mp.mpf(coefficients[3*i+j])==0 for i in range(3) for j in [1,2]):return vertical_reference(element,coefficients,degree),0.
    low=at_order(element,coefficients,degree,128)
    high=at_order(element,coefficients,degree,256)
    discrepancy=max(float(np.max(np.abs(low[p]-high[p]))) for p in high)
    if not math.isfinite(discrepancy) or discrepancy>5e-10:raise RuntimeError(f'Independent tangent reference did not converge: {discrepancy}')
    return high,discrepancy
