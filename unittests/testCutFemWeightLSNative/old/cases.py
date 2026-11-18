def all_cases():
    result = []
    seen = set()

    def add(element, degree, field, mapping='identity', scale=1.0, group='main', validation=0):
        key = (element, degree, field, mapping, scale, validation)
        if key in seen:
            return
        seen.add(key)
        tolerance = 1e-10 if field == 'graph' else (1e-9 if field == 'sphere' else 1e-12)
        if field == 'tiny':
            tolerance = 1e-13
        label = f'e{element}_q{degree}_{field}_{mapping}_s{scale:g}_v{validation}'
        result.append({'id': label, 'element': element, 'degree': degree, 'field': field, 'mapping': mapping, 'scale': scale, 'tolerance': tolerance, 'validation': validation, 'group': group, 'timeout': 45 if field != 'sphere' else 65})

    for element in range(5):
        for degree in range(6):
            for field in ('plane', 'diagonal', 'monotone', 'graph'):
                add(element, degree, field)
        for degree in (0, 2, 5):
            for field in ('full', 'empty', 'boundary_low', 'boundary_high', 'tiny', 'near_affine', 'zero_face'):
                add(element, degree, field, group='edge')
        for degree in (2, 5):
            for scale in (2.0 ** -20, 2.0 ** 20, -3.0):
                add(element, degree, 'plane', scale=scale, group='scale')
            for mapping in ('affine', 'shear', 'stretch'):
                for field in ('plane', 'graph'):
                    add(element, degree, field, mapping=mapping, group='physical')
                if element != 0:
                    add(element, degree, 'plane_y', mapping=mapping, group='physical')
        for scale in (2.0 ** -20, -3.0):
            add(element, 2, 'graph', scale=scale, group='scale')
        for validation in (1, 2):
            add(element, 2, 'graph', validation=validation, group='validation')
    for element in range(5):
        orders = range(6) if element in (0, 1) else (0, 2, 3)
        for degree in orders:
            add(element, degree, 'sphere', group='closed')
    for degree in range(6):
        add(0, degree, 'two_roots', group='closed')
    return result


def selected_cases(selection):
    cases = all_cases()
    if selection == 'quick':
        return [case for case in cases if case['mapping'] == 'identity' and case['scale'] == 1.0 and case['validation'] == 0 and case['field'] in ('plane', 'diagonal', 'monotone')]
    if selection == 'coupled':
        return [case for case in cases if case['field'] == 'graph']
    if selection != 'all':
        return [case for case in cases if case['group'] == selection]
    return cases
