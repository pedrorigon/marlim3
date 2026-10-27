#!/usr/bin/env python3
"""Check phase Jacobian weights against independent source finite differences."""
import ast
from pathlib import Path
import re

root = Path(__file__).resolve().parent.parent
source = (root / 'src/core/SisProdTransient.cpp').read_text()
expressions = {}
for phase in 'GLC':
    match = re.search(r'state\.cells\[state\.lastCell\]\.Dmasschoke' + phase + r' = ([^;]+);', source)
    assert match, f'Missing outlet source derivative for {phase}'
    tree = ast.parse(match[1], mode='eval')
    allowed = (ast.Expression, ast.BinOp, ast.UnaryOp, ast.Mult, ast.Div, ast.Sub,
               ast.Add, ast.USub, ast.UAdd, ast.Name, ast.Constant, ast.Load)
    assert all(isinstance(n, allowed) for n in ast.walk(tree))
    expressions[phase] = compile(tree, '<outlet derivative>', 'eval')

count = 0
for gas_fraction in [0., .25, .5, .75, 1.]:
    for complementary_fraction in [0., .4]:
        densities = [950., 1000.]
        mixture_density = sum(rho * fraction for rho, fraction in
                              zip(densities, [1. - complementary_fraction, complementary_fraction]))
        liquid_mass_fraction = (1. - complementary_fraction) * densities[0] / mixture_density
        complementary_mass_fraction = complementary_fraction * densities[1] / mixture_density
        def mass_sources(pressure):
            total_outflow = 2. * pressure - 10.
            return {'G': -gas_fraction * total_outflow,
                    'L': -(1. - gas_fraction) * liquid_mass_fraction * total_outflow,
                    'C': -(1. - gas_fraction) * complementary_mass_fraction * total_outflow}
        epsilon = 1.e-3
        low, high = mass_sources(12. - epsilon), mass_sources(12. + epsilon)
        variables = {'chokeGasMassFraction': gas_fraction, 'dmaxsup': 2.,
                     'betSup': complementary_fraction, 'rholp': densities[0],
                     'rholc': densities[1], 'rholmix': mixture_density}
        for phase, expression in expressions.items():
            actual = eval(expression, {'__builtins__': {}}, variables)
            expected = (high[phase] - low[phase]) / (2. * epsilon)
            assert abs(actual - expected) < 1.e-10, (gas_fraction, complementary_fraction, phase, actual, expected)
        count += 1
print(f'PASS {count} phase partitions, including pure-liquid and pure-gas limits')
