"""Expand current steady substeps for the existing historical token proof.

The inverse only restores documented scalar scaffolding and literal extractions.
It validates every step-field binding before expansion and keeps changed physics
or state writes visible to check-steady-decomposition.py.
"""
import pathlib
import re

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from module_dependencies_layout import verified_snapshot, candidate_text



def tokens(text):
    text = re.sub(r'//[^\n]*|/\*[\s\S]*?\*/', '', text)
    return re.findall(r'[A-Za-z_]\w*|\d+(?:\.\d*)?(?:[eE][+-]?\d+)?|\S', text)


def span(text, name):
    match = re.search(r'(?m)^[\w:<>]+\s+' + re.escape(name) + r'\s*\(', text)
    assert match, name
    opening = text.index('{', match.end())
    clean = re.sub(r'//[^\n]*|/\*[\s\S]*?\*/', lambda m: ' ' * len(m[0]), text)
    depth = 1
    end = opening + 1
    while depth:
        depth += (clean[end] == '{') - (clean[end] == '}')
        end += 1
    return match.start(), opening, end


FIELDS = (
    'meanPressure', 'meanVoidFraction', 'betmed', 'liquidDensity', 'gasDensity',
    'liquidSpecificHeat', 'gasSpecificHeat', 'meanSuperficialGasVelocity',
    'meanSuperficialLiquidVelocity', 'meanCellLength', 'meanTemperatureGradientCorrection',
    'temperatureSpatialCoefficient', 'annulusResistance', 'mixtureFluxSign',
    'sourceTemperatureGradient',
)


def once(text, old, new):
    assert text.count(old) == 1, (old, text.count(old))
    return text.replace(old, new)


def restore(text):
    text = verified_snapshot({'src/core/SisProdThermal.cpp': text})['src/core/SisProdThermal.cpp']
    internal = pathlib.Path('src/include/SisProdThermalInternal.h').read_text()
    struct = re.search(r'struct SteadyTemperatureStep\s*\{([^}]*)\}', internal)[1]
    assert tokens(struct) == tokens(''.join('double ' + f + ';' for f in FIELDS)), 'steady step scalar types'
    for march, helper, neighbour, shift in (
        ('advanceSteadyTemperature', 'advanceSteadyTemperatureSubsteps', 'leftCell', '- 1'),
        ('advanceReverseSteadyTemperature', 'advanceReverseSteadyTemperatureSubsteps', 'rightCell', '+ 1')):
        a, opening, end = span(text, helper)
        helper_body = text[opening+1:end-1]
        first = helper_body.index('    int subStepCount;')
        expected_prefix = (f'Cel &cell = state.cells[cellIndex]; Cel &{neighbour} = state.cells[cellIndex {shift}];' +
                           ''.join(f'const double {f} = step.{f};' for f in FIELDS) +
                           'double liquidConductivity; double liquidViscosity; double heatFlux;')
        assert tokens(helper_body[:first]) == tokens(expected_prefix), helper + ' input aliases'
        body = helper_body[first:]
        body = once(body, f'double subStepTemperature = subStepLength * (-(-{neighbour}.temp)',
                    f'double subStepTemperature = {neighbour}.temp;\n'
                    f'    subStepTemperature = subStepLength * (-(-{neighbour}.temp)')
        a, opening, end = span(text, march)
        part = text[a:end]
        aggregate = re.search(r'const SteadyTemperatureStep step\s*\{([^}]*)\};', part)
        assert aggregate, march + ' step'
        assert tokens(aggregate[1]) == tokens(','.join('.' + f + '=' + f for f in FIELDS)), march + ' step bindings'
        call = re.search(r'\b' + helper + r'\s*\([^;]*\);', part)
        assert call and call.start() > aggregate.end(), march + ' substep call'
        assert tokens(call[0]) == tokens(helper + '(state, cellIndex, step, heatFluxTemperatureGradient);'), helper + ' arguments'
        assert not part[aggregate.end():call.start()].strip(), march + ' intervening effects'
        part = part[:aggregate.start()] + body + part[call.end():]
        part = once(part, 'double heatFlux;\n        double annulusResistance = apply',
                    f'double relaxedHeatFlux = 0;\n'
                    f'        if (state.steadyIteration > 0)\n            relaxedHeatFlux = {neighbour}.fluxcalmed;\n'
                    '        double annulusResistance = 0.;\n        double heatFlux;\n'
                    '        double gasHeatFlux;\n        annulusResistance = apply')
        kinetic = re.search(r'\n\s*double kineticTerm = compute(?:Reverse)?SteadyKineticTerm\(', part)
        assert kinetic, march + ' kinetic call'
        part = part[:kinetic.start()] + f'\n        double temperatureGradient = (-{neighbour}.temp) / meanCellLength;\n' + part[kinetic.start():]
        part = once(part, f'(void){neighbour}.calor.transperm(0);', f'double heatFlux = {neighbour}.calor.transperm(0);') if f'(void){neighbour}' in part else once(part,
                    f'(void) {neighbour}.calor.transperm(0);', f'double heatFlux = {neighbour}.calor.transperm(0);')
        text = text[:a] + part + text[end:]

    for helper, neighbour in (('computeSteadySourceTerms', 'leftCell'), ('computeReverseSteadySourceTerms', 'rightCell')):
        a, opening, end = span(text, helper)
        part = text[a:end]
        part = once(part, f'double sourceTemperature = {neighbour}.temp;',
                    f'double gasMassSourceTerm = 0.;\n    double liquidMassSourceTerm = 0.;\n'
                    f'    double sourceTemperature = {neighbour}.temp;')
        part = once(part, 'double sourceSpecificHeatRatio;', 'double sourceSpecificHeatRatio = 0.;')
        part = once(part, 'double liquidMassSourceTerm = 0;', 'liquidMassSourceTerm = 0;')
        part = once(part, f'double gasMassSourceTerm = {neighbour}.fontemassGR', f'gasMassSourceTerm = {neighbour}.fontemassGR')
        text = text[:a] + part + text[end:]

    for helper, neighbour in (('applySteadyAnnulusCoupling', 'leftCell'), ('applyReverseSteadyAnnulusCoupling', 'rightCell')):
        a, opening, end = span(text, helper)
        part = text[a:end]
        copies = '\n'.join(f'        {neighbour}.calor.{target} = state.gasCells[j].calor.{source};' for target, source in (
            ('kextern1', 'kint'), ('cpextern1', 'cpint'), ('rhoextern1', 'rhoint'), ('viscextern1', 'viscint')))
        part = once(part, copies + '\n', '')
        # The actual velocity statement stays in the comparison, so changing
        # its literal still fails rather than being silently normalized.
        pattern = (r'(?m)^([ \t]*' + re.escape(neighbour) + r'\.calor\.Vextern1 = [^;]*;)([^\n]*)\n')
        matches = list(re.finditer(pattern, part))
        assert len(matches) == 2, helper + ' external velocity branches'
        for match in reversed(matches):
            part = part[:match.end()] + '\n'.join('    ' + line for line in copies.splitlines()) + '\n' + part[match.end():]
        text = text[:a] + part + text[end:]
    return text
