import functools
import itertools
import math
import numpy as np


def exponents(degree, dimension):
    return [(i, j, k) for i in range(degree + 1) for j in range(degree + 1 - i) for k in range(degree + 1 - i - j if dimension == 3 else 1)]


def full_moment(element, power):
    dimension = 2 if element < 2 else 3
    power = power[:dimension]
    if element % 2:
        return math.prod(math.factorial(n) for n in power) / math.factorial(sum(power) + dimension)
    return math.prod(1 / (n + 1) for n in power)


def cofactor(point, normal):
    result = np.array(normal, copy=True)
    a = 1 + .2 * point[0]
    b = .24 * point[0]
    if len(point) == 2:
        result[0] = normal[0] - b * normal[1]
        result[1] = a * normal[1]
    else:
        c, d = .08 * point[1], .08 * point[0]
        result[0] = normal[0] - b * normal[1] + (b * d - c) * normal[2]
        result[1] = a * (normal[1] - d * normal[2])
        result[2] = a * normal[2]
    return np.sqrt(np.sum(result * result, axis=0))


def shifted_ball(power, center, axes, kind, level):
    dimension = len(center)
    result = 0.0
    for even in itertools.product(*(range(0, p + 1, 2) for p in power[:dimension])):
        total = sum(even)
        coefficient = math.prod(math.comb(p, a) * c ** (p - a) for p, a, c in zip(power, even, center))
        value = math.prod(math.gamma((a + 1) / 2) * r ** (a + 1) for a, r in zip(even, axes)) / math.gamma((total + dimension) / 2 + 1)
        if kind == 'delta':
            value *= (total + dimension) / (2 * level)
        result += coefficient * value
    return result


def sphere_grid(dimension, order):
    nodes, weights = np.polynomial.legendre.leggauss(order)
    if dimension == 2:
        theta = math.pi * (1 + nodes)
        return np.array([np.cos(theta), np.sin(theta)]), math.pi * weights
    theta_nodes, theta_weights = np.polynomial.legendre.leggauss(2 * order)
    theta = math.pi * (1 + theta_nodes)[None, :]
    z = nodes[:, None]
    radius = np.sqrt(1 - z * z)
    shape = (order, 2 * order)
    return np.array([np.broadcast_to(radius * np.cos(theta), shape), np.broadcast_to(radius * np.sin(theta), shape), np.broadcast_to(z, shape)]).reshape(3, -1), (weights[:, None] * math.pi * theta_weights[None, :]).reshape(-1)


def ball_reference(element, field, degree, order):
    dimension = 2 if element < 2 else 3
    radius = 1 / 128 if field == 'tiny_ball' else (.15 if dimension == 2 else .10)
    center = np.array([radius if field == 'touch_ball' else .29, .23, .19][:dimension])
    axes = np.array([.18, .07, .05][:dimension]) if field == 'ellipse' else np.full(dimension, radius)
    level = 1 if field == 'ellipse' else radius * radius
    normal, weight = sphere_grid(dimension, order)
    point = center[:, None] + axes[:, None] * normal
    scaled_normal = math.prod(axes) * normal / axes[:, None]
    metric = np.sqrt(np.sum(scaled_normal * scaled_normal, axis=0))
    physical_metric = cofactor(point, scaled_normal)
    powers = [[np.ones(point.shape[1])] for _ in range(dimension)]
    for d in range(dimension):
        for n in range(1, degree + 2):
            powers[d].append(powers[d][-1] * point[d])
    result = {}
    for power in exponents(degree, dimension):
        h = shifted_ball(power, center, axes, 'volume', level)
        delta = shifted_ball(power, center, axes, 'delta', level)
        monomial = math.prod(powers[d][power[d]] for d in range(dimension))
        surface = float(np.sum(weight * metric * monomial))
        physical = float(np.sum(weight * physical_metric * monomial))
        result[power] = np.array([h, delta, surface, physical])
    return result


def bubble_reference(dimension, degree, order):
    nodes, weights = np.polynomial.legendre.leggauss(order)
    parameter = .5 * (nodes + 1)
    result = {}
    if dimension == 2:
        constant = 1 / 64
        roots = sorted(r.real for r in np.roots([1, -2, 1, -4 * constant]) if abs(r.imag) < 1e-10 and 0 < r.real < 1)
        span = roots[1] - roots[0]
        x = roots[0] + span * np.sin(.5 * math.pi * parameter) ** 2
        measure = .25 * math.pi * span * weights * np.sin(math.pi * parameter)
        third = sum(np.roots([1, -2, 1, -4 * constant]).real)-sum(roots)
        sine, cosine = np.sin(.5*math.pi*parameter), np.cos(.5*math.pi*parameter)
        discriminant = span*sine*cosine*np.sqrt(x*(third-x))
        delta_weight = .5*math.pi*weights/np.sqrt(x*(third-x))
        low = 2*constant/(x*(1-x)+discriminant)
        high = 1-x-low
        point_low, point_high = np.array([x, low]), np.array([x, high])
        grad_low = np.array([low*(1-2*x-low), x*(1-x-2*low)])
        grad_high = np.array([high*(1-2*x-high), x*(1-x-2*high)])
        norm_low = np.sqrt(np.sum(grad_low*grad_low, axis=0))
        norm_high = np.sqrt(np.sum(grad_high*grad_high, axis=0))
        phys_low, phys_high = cofactor(point_low, grad_low), cofactor(point_high, grad_high)
        for i,j,k in exponents(degree, dimension):
            geometric = sum(high**(j-p)*low**p for p in range(j+1))/(j+1)
            factor = measure*x**i
            h = np.sum(factor*(discriminant/x)*geometric)
            delta = np.sum(delta_weight*x**i*(low**j+high**j))
            surface = np.sum(delta_weight*x**i*(low**j*norm_low+high**j*norm_high))
            physical = np.sum(delta_weight*x**i*(low**j*phys_low+high**j*phys_high))
            result[i,j,k] = np.array([h,delta,surface,physical])
        return result
    constant = 1/512
    roots_z = sorted(r.real for r in np.roots([-1,3,-3,1,-27*constant]) if abs(r.imag) < 1e-10 and 0 < r.real < 1)
    span = roots_z[1]-roots_z[0]
    z = (roots_z[0]+span*np.sin(.5*math.pi*parameter)**2)[:,None]
    theta = .5*math.pi*nodes[None,:]
    measure = (.25*math.pi*span*weights*np.sin(math.pi*parameter))[:,None]*(.5*math.pi*weights[None,:])
    roots = np.array([sorted(r.real for r in np.roots([1,-2*(1-t),(1-t)**2,-4*constant/t])) for t in z[:,0]])
    half = .5*(roots[:,1]-roots[:,0])[:,None]
    y = .5*(roots[:,0]+roots[:,1])[:,None]+half*np.sin(theta)
    separation = half*np.cos(theta)*np.sqrt((roots[:,2,None]-y)/y)
    high = .5*(1-y-z+separation)
    low = .5*(1-y-z-separation)
    volume_measure = measure*half*np.cos(theta)*separation
    delta_measure = measure/(z*np.sqrt(y*(roots[:,2,None]-y)))
    point_low = np.array([low,y,np.broadcast_to(z,y.shape)])
    point_high = np.array([high,y,np.broadcast_to(z,y.shape)])
    grad_low = np.array([y*z*(1-2*low-y-z),low*z*(1-low-2*y-z),low*y*(1-low-y-2*z)])
    grad_high = np.array([y*z*(1-2*high-y-z),high*z*(1-high-2*y-z),high*y*(1-high-y-2*z)])
    norm_low = np.sqrt(np.sum(grad_low*grad_low,axis=0))
    norm_high = np.sqrt(np.sum(grad_high*grad_high,axis=0))
    phys_low,phys_high = cofactor(point_low,grad_low),cofactor(point_high,grad_high)
    for i,j,k in exponents(degree,dimension):
        factor=y**j*z**k
        geometric=sum(high**(i-p)*low**p for p in range(i+1))/(i+1)
        result[i,j,k]=np.array([np.sum(volume_measure*factor*geometric),np.sum(delta_measure*factor*(low**i+high**i)),np.sum(delta_measure*factor*(low**i*norm_low+high**i*norm_high)),np.sum(delta_measure*factor*(low**i*phys_low+high**i*phys_high))])
    return result


def plane_reference(element, field, power):
    dimension = 2 if element < 2 else 3
    axis = dimension-1 if field in ('plane','boundary_low','boundary_high') else 0
    transverse = [power[d] for d in range(dimension) if d != axis]
    if field == 'oblique':
        c=.43
        h=full_moment(element,power)-c**(sum(power)+dimension)*math.prod(math.factorial(p) for p in power)/math.factorial(sum(power)+dimension)
        delta=c**(sum(power)+dimension-1)*math.prod(math.factorial(p) for p in power)/math.factorial(sum(power)+dimension-1)
        return np.array([h,delta,math.sqrt(dimension)*delta,0])
    cuts=[.375] if field == 'plane' else ([0.] if field == 'boundary_low' else ([1.] if field == 'boundary_high' else [.21,.68]))
    i=power[axis]
    n=sum(transverse)+dimension-1
    factor=math.prod(math.factorial(p) for p in transverse)/math.factorial(n) if element%2 else math.prod(1/(p+1) for p in transverse)
    def primitive(x):
        if element%2:
            return sum(np.longdouble(math.comb(n,j))*(-1)**j*np.longdouble(x)**(i+j+1)/(i+j+1) for j in range(n+1))
        return np.longdouble(x)**(i+1)/(i+1)
    projection=sum(c**i*((1-c)**n if element%2 else 1)*factor for c in cuts)
    if field in ('plane','boundary_low','boundary_high'):
        h=full_moment(element,power)-factor*primitive(cuts[0])
        if field.startswith('boundary_'):
            projection *= .5
        return np.array([h,projection,projection,0],dtype=float)
    h=full_moment(element,power)-factor*(primitive(cuts[1])-primitive(cuts[0]))
    return np.array([h,projection/(cuts[1]-cuts[0]),projection,0],dtype=float)


@functools.lru_cache(maxsize=None)
def reference(element,field,degree):
    dimension=2 if element<2 else 3
    discrepancy=0.0
    if field in ('ball','tiny_ball','touch_ball','ellipse','bubble'):
        builder=(lambda n: bubble_reference(dimension,degree,n)) if field=='bubble' else (lambda n: ball_reference(element,field,degree,n))
        low,high=builder(80),builder(144)
        discrepancy=max(float(np.max(np.abs(low[p]-high[p]))) for p in high)
        if discrepancy>5e-10:
            raise RuntimeError(f'Independent {field} reference did not converge: {discrepancy}')
    else:
        high={p:(np.array([full_moment(element,p) if field=='full' else 0,0,0,0]) if field in ('full','empty') else plane_reference(element,field,p)) for p in exponents(degree,dimension)}
    return high,discrepancy
