import functools
import itertools
import math

import numpy as np


def dimension(element):
    return 1 if element == 0 else (2 if element < 3 else 3)


def powers(degree, dim):
    result = []
    for total in range(degree + 1):
        for first in range(total, -1, -1):
            if dim == 1 and first != total:
                continue
            for second in range(total - first, -1, -1):
                third = total - first - second
                if dim < 3 and third:
                    continue
                result.append((first, second, third))
    return result


def full_moment(element, exponent):
    exponent = exponent[:dimension(element)]
    if element in (2, 4):
        return math.prod(math.factorial(value) for value in exponent) / math.factorial(sum(exponent) + len(exponent))
    if any(value % 2 for value in exponent):
        return 0.0
    return math.prod(2.0 / (value + 1) for value in exponent)


def full_physical_moment(element, exponent, mapping):
    dim = dimension(element)
    if mapping == 'affine':
        return full_moment(element, exponent) * math.prod(1.25 + 0.25 * d for d in range(dim))
    if mapping == 'stretch' or (mapping == 'shear' and dim == 1):
        raised = (exponent[0] + 1, exponent[1], exponent[2])
        return full_moment(element, exponent) + 0.1 * full_moment(element, raised)
    return full_moment(element, exponent)


@functools.lru_cache(maxsize=None)
def gauss(order):
    return np.polynomial.legendre.leggauss(order)


def interval(low, high, order):
    nodes, weights = gauss(order)
    return 0.5 * (low + high) + 0.5 * (high - low) * nodes, 0.5 * (high - low) * weights


def integral_power(low, high, exponent):
    result = np.zeros_like(np.asarray(low) + np.asarray(high), dtype=float)
    for k in range(exponent + 1):
        result = result + np.asarray(high) ** (exponent - k) * np.asarray(low) ** k
    return (np.asarray(high) - np.asarray(low)) * result / (exponent + 1)


def root_for(field, element, points):
    dim = dimension(element)
    simplex = element in (2, 4)
    x = points[:, 0]
    y = points[:, 1] if dim > 1 else np.zeros_like(x)
    if field in ('plane', 'plane_y'):
        return np.full(len(points), 0.25 if simplex else 0.125)
    if field == 'boundary_low':
        return np.full(len(points), 0.0 if simplex else -1.0)
    if field == 'boundary_high':
        return np.ones(len(points))
    if field == 'tiny':
        return np.full(len(points), 1.0 - 1.0 / 256.0)
    if field == 'diagonal':
        return 0.5 - points[:, :dim - 1].sum(axis=1)
    if field == 'zero_face':
        return np.zeros(len(points))
    if field == 'monotone' or (field == 'graph' and element == 0):
        return np.full(len(points), 0.5 / (1.0 + math.sqrt(1.1)))
    if field == 'near_affine':
        epsilon = 1.0 / 1073741824.0
        if dim == 1:
            return np.full(len(points), 0.5 / (1.0 + math.sqrt(1.0 + epsilon)))
        return 0.25 - epsilon * x * x
    if field == 'graph':
        if element == 1:
            value = 0.1 + 0.2 * x * x
            return 2.0 * value / (1.0 + np.sqrt(1.0 + 0.16 * x * x * value))
        if element == 2:
            value = 0.2 + 0.15 * x * x
            coefficient = 0.81 * x
            linear = 1.0 + coefficient * (1.0 - x)
            return 2.0 * value / (linear + np.sqrt(linear * linear - 4.0 * coefficient * value))
        if element == 3:
            value = -0.05 + 0.12 * x * x + 0.07 * y * y
            return 2.0 * value / (1.0 + np.sqrt(1.0 + 0.12 * x * x * value))
        value = 0.2 + 0.08 * x * x - 0.02 * x * y
        coefficient = 0.768 * x * y
        linear = 1.0 + coefficient * (1.0 - x - y)
        return 2.0 * value / (linear + np.sqrt(linear * linear - 4.0 * coefficient * value))
    raise ValueError(field)


def analytic_gradient(field, element, points):
    dim = dimension(element)
    x = points[:, 0]
    y = points[:, 1] if dim > 1 else np.zeros_like(x)
    z = points[:, 2] if dim > 2 else np.zeros_like(x)
    gradient = np.zeros_like(points)
    if field in ('plane', 'plane_y', 'boundary_low', 'boundary_high', 'tiny'):
        gradient[:, 1 if field == 'plane_y' and dim > 1 else 0] = 1.0
    elif field == 'diagonal':
        gradient[:] = 3.0
    elif field == 'monotone' or (field == 'graph' and element == 0):
        gradient[:, 0] = 1.0 + 0.2 * x
    elif field == 'near_affine':
        gradient[:, -1] = 1.0
        gradient[:, 0] += 2.0 * x / 1073741824.0
    elif field == 'zero_face':
        gradient[:, -1] = 1.0 + 0.125 * x
        gradient[:, 0] += 0.125 * points[:, -1]
    elif field == 'graph' and element == 1:
        gradient[:, 0] = -0.4 * x + 0.08 * x * y * y
        gradient[:, 1] = 1.0 + 0.08 * x * x * y
    elif field == 'graph' and element == 2:
        gradient[:, 0] = -0.3 * x + 0.81 * y * (1.0 - 2.0 * x - y)
        gradient[:, 1] = 1.0 + 0.81 * x * (1.0 - x - 2.0 * y)
    elif field == 'graph' and element == 3:
        gradient[:, 0] = -0.24 * x + 0.06 * x * z * z
        gradient[:, 1] = -0.14 * y
        gradient[:, 2] = 1.0 + 0.06 * x * x * z
    elif field == 'graph' and element == 4:
        gradient[:, 0] = -0.16 * x + 0.02 * y + 0.768 * y * z * (1.0 - 2.0 * x - y - z)
        gradient[:, 1] = 0.02 * x + 0.768 * x * z * (1.0 - x - 2.0 * y - z)
        gradient[:, 2] = 1.0 + 0.768 * x * y * (1.0 - x - y - 2.0 * z)
    else:
        raise ValueError((field, element))
    return gradient


def cofactor_norm(mapping, points, gradient):
    dim = points.shape[1]
    transformed = gradient.copy()
    if mapping == 'affine':
        stretch = np.array([1.25 + 0.25 * d for d in range(dim)])
        transformed *= np.prod(stretch) / stretch
    elif mapping == 'shear' and dim > 1:
        transformed[:, 0] -= 0.3 * points[:, 0] * gradient[:, 1]
    elif mapping == 'stretch' or (mapping == 'shear' and dim == 1):
        transformed[:, 1:] *= (1.0 + 0.1 * points[:, 0])[:, None]
    return np.linalg.norm(transformed, axis=1)


def graph_grid(element, field, order):
    dim = dimension(element)
    simplex = element in (2, 4)
    axis = dim - 1 if field in ('graph', 'near_affine', 'zero_face', 'diagonal') else (1 if field == 'plane_y' else 0)
    if field == 'monotone':
        axis = 0
    remaining = [d for d in range(dim) if d != axis]
    boundary = 0.5 if field in ('boundary_low', 'boundary_high') or (field == 'zero_face' and simplex) else 1.0
    if dim == 1:
        points = np.zeros((1, 1))
        weights = np.ones(1)
    elif not simplex and field == 'diagonal':
        blocks = []
        weight_blocks = []
        surface_blocks = []
        if dim == 2:
            x, wx = interval(-0.5, 1.0, order)
            points = np.column_stack((x, 0.5 - x))
            weights = wx
            return points, weights, axis, points[:, axis], np.ones(len(points)), weights
        for low, high, full in ((-1.0, 0.5, False), (0.5, 1.0, False), (0.5, 1.0, True)):
            x, wx = interval(low, high, order)
            for first, weight in zip(x, wx):
                left = 1.5 - first if full else max(-1.0, -0.5 - first)
                right = 1.0 if full else min(1.0, 1.5 - first)
                y, wy = interval(left, right, order)
                z = np.full(order, -1.0) if full else 0.5 - first - y
                blocks.append(np.column_stack((np.full(order, first), y, z)))
                weight_blocks.append(weight * wy)
                surface_blocks.append(np.zeros(order) if full else weight * wy)
        points = np.concatenate(blocks)
        weights = np.concatenate(weight_blocks)
        return points, weights, axis, points[:, axis], np.ones(len(points)), np.concatenate(surface_blocks)
    elif not simplex:
        first, first_weights = interval(-1.0, 1.0, order)
        if dim == 2:
            points = np.zeros((order, dim))
            points[:, remaining[0]] = first
            weights = first_weights
        else:
            points = np.zeros((order * order, dim))
            points[:, remaining[0]] = np.repeat(first, order)
            points[:, remaining[1]] = np.tile(first, order)
            weights = np.repeat(first_weights, order) * np.tile(first_weights, order)
    else:
        probe = np.zeros((1, dim))
        if field == 'graph':
            quadratic = 0.15 if dim == 2 else 0.08
            maximum = 1.6 / (1.0 + math.sqrt(1.0 + 3.2 * quadratic))
        elif field == 'near_affine':
            maximum = 1.5 / (1.0 + math.sqrt(1.0 - 3.0 / 1073741824.0))
        else:
            maximum = 1.0 - root_for(field, element, probe)[0]
        first, first_weights = interval(0.0, maximum, order)
        if dim == 2:
            points = np.zeros((order, dim))
            points[:, remaining[0]] = first
            weights = first_weights
        else:
            blocks = []
            weight_blocks = []
            for value, weight in zip(first, first_weights):
                if field == 'graph':
                    upper = (0.8 - value - 0.08 * value * value) / (1.0 - 0.02 * value)
                elif field == 'near_affine':
                    upper = 0.75 - value + value * value / 1073741824.0
                else:
                    upper = maximum - value
                second, second_weights = interval(0.0, upper, order)
                block = np.zeros((order, dim))
                block[:, remaining[0]] = value
                block[:, remaining[1]] = second
                blocks.append(block)
                weight_blocks.append(weight * second_weights)
            points = np.concatenate(blocks)
            weights = np.concatenate(weight_blocks)
    points[:, axis] = root_for(field, element, points)
    upper = 1.0 - points[:, remaining].sum(axis=1) if simplex else np.ones(len(points))
    return points, weights, axis, points[:, axis], upper, boundary * weights


def ball_moment(dim, exponent, center, radius, kind):
    result = 0.0
    for indices in itertools.product(*(range(0, exponent[d] + 1, 2) for d in range(dim))):
        coefficient = math.prod(math.comb(exponent[d], indices[d]) * center ** (exponent[d] - indices[d]) for d in range(dim))
        total = sum(indices)
        angular = 2.0 * math.prod(math.gamma(0.5 * (index + 1)) for index in indices) / math.gamma(0.5 * (total + dim))
        if kind == 'h':
            term = angular * radius ** (total + dim) / (total + dim)
        elif kind == 'surface':
            term = angular * radius ** (total + dim - 1)
        else:
            term = 0.5 * angular * radius ** (total + dim - 2)
        result += coefficient * term
    return result


@functools.lru_cache(maxsize=None)
def reference_unscaled(element, field, mapping, maximum_degree=6, order=64):
    dim = dimension(element)
    result = {}
    if field in ('full', 'empty', 'sphere', 'two_roots') or (field == 'diagonal' and element in (2, 4)):
        for exponent in powers(maximum_degree, dim):
            if field in ('full', 'empty'):
                h = full_moment(element, exponent) if field == 'full' else 0.0
                hp = full_physical_moment(element, exponent, mapping) if field == 'full' else 0.0
                result[exponent] = (h, 0.0, 0.0, 0.0, hp)
            elif field in ('sphere', 'two_roots'):
                if mapping != 'identity':
                    raise ValueError('Sphere analytic references use identity geometry.')
                center = (1.0 / 3.0 if dim == 2 else 0.25) if element in (2, 4) else 0.0
                radius = (0.2 if dim == 2 else 0.125) if element in (2, 4) else (0.5 if dim == 1 else 0.6)
                h = ball_moment(dim, exponent, center, radius, 'h')
                if field == 'two_roots':
                    h = full_moment(element, exponent) - h
                delta = ball_moment(dim, exponent, center, radius, 'delta')
                surface = ball_moment(dim, exponent, center, radius, 'surface')
                result[exponent] = (h, delta, surface, surface, h)
            else:
                if mapping != 'identity':
                    raise ValueError('Simplex diagonal references use identity geometry.')
                total = sum(exponent)
                h = full_moment(element, exponent) * (1.0 - 0.5 ** (total + dim))
                delta = 0.5 ** (total + dim - 1) * math.prod(math.factorial(value) for value in exponent[:dim]) / (3.0 * math.factorial(total + dim - 1))
                surface = 3.0 * math.sqrt(dim) * delta
                result[exponent] = (h, delta, surface, surface, h)
        return result
    points, weights, axis, lower, upper, surface_weights = graph_grid(element, field, order)
    gradient = analytic_gradient(field, element, points)
    delta_weights = surface_weights / np.abs(gradient[:, axis])
    surface_factor = np.linalg.norm(gradient, axis=1)
    physical_factor = cofactor_norm(mapping, points, gradient)
    for exponent in powers(maximum_degree, dim):
        other = np.ones(len(points))
        on_surface = np.ones(len(points))
        for d in range(dim):
            on_surface *= points[:, d] ** exponent[d]
            if d != axis:
                other *= points[:, d] ** exponent[d]
        volume = other * integral_power(lower, upper, exponent[axis])
        h = float(np.dot(weights, volume))
        delta = float(np.dot(delta_weights, on_surface))
        surface = float(np.dot(delta_weights, on_surface * surface_factor))
        physical = float(np.dot(delta_weights, on_surface * physical_factor))
        if mapping == 'affine':
            hp = h * math.prod(1.25 + 0.25 * d for d in range(dim))
        elif mapping == 'stretch' or (mapping == 'shear' and dim == 1):
            if axis == 0:
                hp = float(np.dot(weights, other * (integral_power(lower, upper, exponent[axis]) + 0.1 * integral_power(lower, upper, exponent[axis] + 1))))
            else:
                hp = float(np.dot(weights, volume * (1.0 + 0.1 * points[:, 0])))
        else:
            hp = h
        result[exponent] = (h, delta, surface, physical, hp)
    return result


def reference(case):
    low = reference_unscaled(case['element'], case['field'], case['mapping'], 6, 40)
    high = reference_unscaled(case['element'], case['field'], case['mapping'], 6, 72)
    discrepancy = max(abs(a - b) for key in high for a, b in zip(low[key], high[key]))
    if discrepancy > 2e-10:
        raise RuntimeError('Independent reference quadrature did not converge: ' + str(discrepancy))
    scale = case['scale']
    result = {}
    for exponent, original in high.items():
        h, delta, surface, physical, hp = original
        if scale < 0.0:
            h = full_moment(case['element'], exponent) - h
            hp = full_physical_moment(case['element'], exponent, case['mapping']) - hp
        result[exponent] = (h, delta / abs(scale), surface, surface, physical, physical, hp)
    return result, discrepancy


def outside_degree(case, physical=False):
    field = case['field']
    dim = dimension(case['element'])
    mapping = case['mapping'] if physical else 'identity'
    if field in ('full', 'empty'):
        return 0
    affine = field in ('plane', 'plane_y', 'diagonal', 'boundary_low', 'boundary_high', 'tiny')
    polynomial = field == 'monotone' or (dim == 1 and field in ('graph', 'near_affine', 'zero_face'))
    if not affine and not polynomial:
        return None
    degree = 0 if affine else 1
    if mapping in ('identity', 'affine') or dim == 1:
        return degree
    axis = 1 if field == 'plane_y' else 0
    if field == 'diagonal':
        return None
    if axis == 0:
        return degree
    if mapping == 'stretch':
        return degree + 1
    return None
