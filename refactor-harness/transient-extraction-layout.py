"""Reconstruct the transient calculation phases against the pre-extraction source.

Only checked local aliases, result bindings and helper calls are expanded.
The valve reduction is checked separately: each candidate retains its source
calculation and shared state effects, and strict comparisons retain valve order.
No simulation, reference generation, or performance measurement is performed.
"""
import pathlib
import re
import subprocess
import sys

CHECKPOINT = 'abaf57a'
SIGNATURES = {'applyChokeTimeStepControls': 'void applyChokeTimeStepControls(const TransientStepState &state)',
 'applyScheduledTimeStepEvents': 'void applyScheduledTimeStepEvents(const TransientStepState &state)',
 'bufferedOutletUsesPressureDrop': 'bool bufferedOutletUsesPressureDrop(const TransientStepState &state, double '
                                   'chokeGasMassFraction,\n'
                                   '                                    double delp, double masentrada, double '
                                   'maxSup)',
 'currentOutletCanRemainOpen': 'bool currentOutletCanRemainOpen(const TransientStepState &state, double delp)',
 'currentOutletGasMassFraction': 'double currentOutletGasMassFraction(const TransientStepState &state, double '
                                 'alfSup,\n'
                                 '                                    double massgas, double masentrada, double '
                                 'romix)',
 'currentOutletUsesPressureDrop': 'bool currentOutletUsesPressureDrop(const TransientStepState &state, double '
                                  'chokeGasMassFraction,\n'
                                  '                                   double delp, double masentrada, double '
                                  'maxSup)',
 'facePhaseDensitiesOf': 'FacePhaseDensities facePhaseDensitiesOf(const TransientStepState &state, int i,\n'
                         '                                             const ProFlu &flud, double tmed, double '
                         'betI)',
 'fillMiniTableCorner': 'void fillMiniTableCorner(const TransientStepState &state, ProFlu &fluC, ProFlu &flui,\n'
                        '                         int pressureIndex, int temperatureIndex)',
 'finalizeTransientCellState': 'void finalizeTransientCellState(const TransientSolveState &state, const ProFlu '
                               '*fluiRev)',
 'penalizeOscillatingFractions': 'void penalizeOscillatingFractions(const TransientStepState &state)',
 'prepareCouplingModelFlags': 'void prepareCouplingModelFlags(const TransientSolveState &state, int celpos)',
 'preparePigFaceFractions': 'void preparePigFaceFractions(const TransientStepState &state)',
 'refreshAccessoryFluidMiniTable': 'void refreshAccessoryFluidMiniTable(const TransientStepState &state, int i, '
                                   'ProFlu &accessoryFluid)',
 'restartCouplingFractions': 'void restartCouplingFractions(const TransientSolveState &state, int kontaAcop,\n'
                             '                              double alfRev, double betRev)',
 'restoreIntermediateCouplingState': 'void restoreIntermediateCouplingState(const TransientSolveState &state, int '
                                     'kontaAcop,\n'
                                     '                                      double gasMassSource, double '
                                     'liquidMassSource,\n'
                                     '                                      double complementaryMassSource)',
 'restrictTimeStepByTransport': 'void restrictTimeStepByTransport(const TransientStepState &state)',
 'supplyAbsentGasProperties': 'void supplyAbsentGasProperties(ProFlu &flui, const std::pair<double, int> '
                              '(&titVec)[4])',
 'supplyAbsentLiquidProperties': 'void supplyAbsentLiquidProperties(ProFlu &flui, const std::pair<double, int> '
                                 '(&titVec)[4])',
 'updateActivePigState': 'void updateActivePigState(const TransientStepState &state, int i)',
 'updateBufferedNetworkSources': 'void updateBufferedNetworkSources(const TransientStepState &state, double '
                                 'betRev)',
 'updateInteriorAuxiliaryPressure': 'double updateInteriorAuxiliaryPressure(const TransientStepState &state, int '
                                    'i)',
 'updateInteriorVolumeFlows': 'void updateInteriorVolumeFlows(const TransientStepState &state, int i, double '
                              'area)',
 'updateOutletAuxiliaryPressure': 'void updateOutletAuxiliaryPressure(const TransientStepState &state, int i)',
 'updateOutletComplementaryFraction': 'void updateOutletComplementaryFraction(const TransientStepState &state, '
                                      'int i)',
 'updateOutletVoidFraction': 'void updateOutletVoidFraction(const TransientStepState &state, int i, double '
                             'alfrev)',
 'updateOutletVolumeFlows': 'void updateOutletVolumeFlows(const TransientStepState &state, int i)',
 'updateTransientMovingMeans': 'void updateTransientMovingMeans(const TransientSolveState &state)',
 'valveRestrictedTimeStep': 'double valveRestrictedTimeStep(const TransientStepState &state, int i, double dtaux)'}
LITERAL = {
    'updateInteriorVolumeFlows', 'updateOutletAuxiliaryPressure', 'updateOutletVolumeFlows',
    'applyChokeTimeStepControls', 'restrictTimeStepByTransport', 'applyScheduledTimeStepEvents',
    'penalizeOscillatingFractions', 'updateOutletVoidFraction', 'updateOutletComplementaryFraction',
    'preparePigFaceFractions', 'updateBufferedNetworkSources', 'prepareCouplingModelFlags',
    'restartCouplingFractions', 'restoreIntermediateCouplingState', 'finalizeTransientCellState',
    'updateTransientMovingMeans', 'updateActivePigState', 'supplyAbsentGasProperties',
    'supplyAbsentLiquidProperties', 'refreshAccessoryFluidMiniTable', 'fillMiniTableCorner',
}
VALUE = {'updateInteriorAuxiliaryPressure': 'area',
         'currentOutletGasMassFraction': 'chokeGasMassFraction'}
PREDICATES = {'currentOutletUsesPressureDrop', 'currentOutletCanRemainOpen',
              'bufferedOutletUsesPressureDrop'}
USED = set()


def masked(text):
    pattern = r'//[^\n]*|/\*[\s\S]*?\*/|"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\''
    return re.sub(pattern, lambda m: re.sub(r'[^\n]', ' ', m[0]), text)


def canonical(text):
    pattern = r'//[^\n]*|/\*[\s\S]*?\*/|"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\'|[A-Za-z_]\w*|\d+(?:\.\d*)?(?:[eE][+-]?\d+)?|\S'
    return ' '.join(m[0] for m in re.finditer(pattern, text)
                    if not m[0].startswith(('//', '/*')))


def closing(text, opening, left, right):
    clean, depth = masked(text), 0
    for index in range(opening, len(clean)):
        depth += (clean[index] == left) - (clean[index] == right)
        if depth == 0:
            return index
    raise AssertionError('unbalanced ' + left)


def parts(source, name):
    matches = list(re.finditer(r'(?m)^[\w:<>]+\s+' + re.escape(name) + r'\s*\(', masked(source)))
    assert len(matches) == 1, (name, len(matches))
    start = matches[0].start()
    opening = source.index('{', matches[0].end())
    end = closing(source, opening, '{', '}')
    return source[start:opening].strip(), source[opening+1:end]


def comma_parts(text):
    result, start, stack = [], 0, []
    for index, char in enumerate(text):
        if char in '(<[':
            stack.append(char)
        elif char in ')>]':
            assert stack, text
            stack.pop()
        elif char == ',' and not stack:
            result.append(text[start:index].strip())
            start = index + 1
    result.append(text[start:].strip())
    return result


def replace_once(text, before, after):
    # Locate tokens without discarding the formatting of the remaining body.
    pattern = r'//[^\n]*|/\*[\s\S]*?\*/|"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\'|[A-Za-z_]\w*|\d+(?:\.\d*)?(?:[eE][+-]?\d+)?|\S'
    spans = [m for m in re.finditer(pattern, text) if not m[0].startswith(('//', '/*'))]
    # String literals containing spaces are single tokens, so use the same lexer.
    wanted = [m[0] for m in re.finditer(pattern, before) if not m[0].startswith(('//', '/*'))]
    values = [m[0] for m in spans]
    positions = [i for i in range(len(values)-len(wanted)+1) if values[i:i+len(wanted)] == wanted]
    assert len(positions) == 1, ('replacement', before, len(positions))
    index = positions[0]
    return text[:spans[index].start()] + after + text[spans[index+len(wanted)-1].end():]


def equal(label, expected, actual):
    expected, actual = canonical(expected), canonical(actual)
    if expected != actual:
        old, new = expected.split(), actual.split()
        position = next((i for i, pair in enumerate(zip(old, new)) if pair[0] != pair[1]),
                        min(len(old), len(new)))
        raise AssertionError(f'{label}: token {position}\n'
                             f'expected: {" ".join(old[max(0,position-5):position+8])}\n'
                             f'actual:   {" ".join(new[max(0,position-5):position+8])}')


def helper(source, name):
    signature, body = parts(source, name)
    equal(name + ' interface', SIGNATURES[name], signature)
    USED.add(name)
    return signature, body


def bind(signature, body, args):
    opening = signature.index('(')
    parameters = signature[opening+1:signature.rfind(')')]
    fields = []
    for parameter in comma_parts(parameters):
        array = re.search(r'\(&(\w+)\)', parameter)
        fields.append(array[1] if array else re.findall(r'\b\w+\b', parameter)[-1])
    assert len(fields) == len(args), (signature, args)
    bindings = dict(zip(fields, args))
    return re.sub(r'\b(?:' + '|'.join(map(re.escape, fields)) + r')\b',
                  lambda m: bindings[m[0]], body)


def calls(body, names):
    pattern = r'\b(' + '|'.join(map(re.escape, sorted(names, key=len, reverse=True))) + r')\s*\('
    for match in re.finditer(pattern, masked(body)):
        opening = body.index('(', match.start())
        finish = closing(body, opening, '(', ')')
        yield match[1], match.start(), finish+1, comma_parts(body[opening+1:finish])


def live_cell_aliases(body):
    if 'Cel &cell' not in body:
        return body
    body = replace_once(body, 'Cel &cell = state.cells[i]; Cel &upstreamCell = state.cells[i - 1];', '')
    body = re.sub(r'\bupstreamCell\b', 'state.cells[i - 1]', body)
    return re.sub(r'\bcell\b', 'state.cells[i]', body)


def corner_body(source, args):
    signature, body = helper(source, 'fillMiniTableCorner')
    body = replace_once(body,
        'const double &pressure = pressureIndex == 0 ? flui.miniTabDin.pmin : flui.miniTabDin.pmax;'
        'const double &temperature = temperatureIndex == 0 ? flui.miniTabDin.tmin : flui.miniTabDin.tmax;', '')
    assert args[3] in ('0', '1') and args[4] in ('0', '1'), args
    body = re.sub(r'\bpressure\b', 'flui.miniTabDin.' + ('pmin' if args[3] == '0' else 'pmax'), body)
    body = re.sub(r'\btemperature\b', 'flui.miniTabDin.' + ('tmin' if args[4] == '0' else 'tmax'), body)
    return bind(signature, body, args)


def expand(source, name, body=None, depth=0):
    assert depth < 12, name
    if body is None:
        _, body = parts(source, name)
    body = live_cell_aliases(body)

    # Bind both physical results, retaining their property evaluation order.
    density_pattern = (r'const FacePhaseDensities densities = facePhaseDensitiesOf\(([^;]+)\);'
                       r'\s*double rhol = densities\.liquid;\s*double gasDensity = densities\.gas;')
    match = re.search(density_pattern, body)
    if match:
        signature, density = helper(source, 'facePhaseDensitiesOf')
        density = replace_once(density, 'return {rhol, gasDensity};', '')
        density = bind(signature, density, comma_parts(match[1]))
        if name in ('updateInteriorVolumeFlows', 'updateOutletVolumeFlows'):
            density = replace_once(density, 'double rhol =', 'rhol =')
        body = body[:match.start()] + density + body[match.end():]
        return expand(source, name, body, depth+1)

    for called, value in VALUE.items():
        pattern = r'\bdouble\s+' + value + r'\s*=\s*' + called + r'\s*\('
        match = re.search(pattern, masked(body))
        if match:
            opening = body.index('(', match.start())
            finish = closing(body, opening, '(', ')')
            assert body[finish+1:].lstrip().startswith(';')
            statement_end = finish + 1 + body[finish+1:].index(';') + 1
            args = comma_parts(body[opening+1:finish])
            signature, core = helper(source, called)
            core = replace_once(core, 'return ' + value + ';', '')
            core = bind(signature, core, args)
            core = expand(source, called, core, depth+1)
            body = body[:match.start()] + core + body[statement_end:]
            return expand(source, name, body, depth+1)

    for predicate in PREDICATES:
        match = re.search(r'\bif\s*\(\s*' + predicate + r'\s*\(', masked(body))
        if match:
            outer = body.index('(', match.start())
            call_open = body.index('(', outer+1)
            call_end = closing(body, call_open, '(', ')')
            outer_end = closing(body, outer, '(', ')')
            equal('predicate call shape', ')', body[call_end+1:outer_end+1])
            signature, core = helper(source, predicate)
            assert core.count('return ') == 1
            prefix, expression = core.split('return ', 1)
            assert expression.rstrip().endswith(';')
            rebuilt = prefix + 'if (' + expression.rstrip()[:-1] + ')'
            rebuilt = bind(signature, rebuilt, comma_parts(body[call_open+1:call_end]))
            body = body[:match.start()] + rebuilt + body[outer_end+1:]
            return expand(source, name, body, depth+1)

    found = list(calls(body, LITERAL))
    if found:
        called, start, finish, args = found[0]
        tail = body[finish:]
        assert tail.lstrip().startswith(';'), (called, tail[:30])
        statement_end = finish + tail.index(';') + 1
        if called == 'fillMiniTableCorner':
            core = corner_body(source, args)
        else:
            signature, core = helper(source, called)
            if called == 'updateBufferedNetworkSources':
                core = replace_once(core, 'double sinal;', '')
            core = bind(signature, core, args)
            core = expand(source, called, core, depth+1)
        body = body[:start] + core + body[statement_end:]
        return expand(source, name, body, depth+1)
    return body


def check_valves(old, current):
    signature, body = parts(current, 'restrictTimeStepByValve')
    old_signature, baseline = parts(old, 'restrictTimeStepByValve')
    equal('valve wrapper signature', old_signature, signature)
    marker = 'for (int i = 0; i <= state.input.nvalv; i++) {\n            if (i == 0)'
    first = baseline.index(marker)
    opening = baseline.index('{', first)
    finish = closing(baseline, opening, '{', '}')
    old_candidate = baseline[opening+1:finish]
    candidate_signature, candidate = helper(current, 'valveRestrictedTimeStep')
    candidate = replace_once(candidate, 'int celpos; double valveTimeStep = state.cells[0].dt;', '')
    candidate = replace_once(candidate, 'return valveTimeStep;', '')
    candidate = re.sub(r'\bvalveTimeStep\b', 'dtvec[i]', candidate)
    equal('per-valve arithmetic and shared state effects', old_candidate, candidate)
    prefix = baseline[:baseline.index('        double dtvec[')]
    prefix = replace_once(prefix, 'int celpos;', '')
    equal('initial candidate values',
          'double dtvec[state.input.nvalv + 1];'
          'for (int i = 0; i <= state.input.nvalv; i++) dtvec[i] = state.cells[0].dt;',
          baseline[baseline.index('        double dtvec['):first])
    reduce_start = finish+1
    reduce_end = baseline.index('state.cells[0].dt = dtaux;', reduce_start) + len('state.cells[0].dt = dtaux;')
    equal('source strict ordered reduction',
          'dtaux = dtvec[0];'
          'for (int i = 1; i <= state.input.nvalv; i++) if (dtvec[i] < dtaux) dtaux = dtvec[i];'
          'state.cells[0].dt = dtaux;', baseline[reduce_start:reduce_end])
    replacement = (
        'double minimumTimeStep = valveRestrictedTimeStep(state, 0, dtaux);'
        'for (int i = 1; i <= state.input.nvalv; i++) {'
        'double candidate = valveRestrictedTimeStep(state, i, dtaux);'
        'if (candidate < minimumTimeStep) minimumTimeStep = candidate; }'
        'state.cells[0].dt = minimumTimeStep;')
    equal('constant-storage strict ordered reduction', prefix + replacement + baseline[reduce_end:], body)
    return baseline


def reference(path):
    root = pathlib.Path(__file__).resolve().parents[1]
    return subprocess.check_output(['git', '-C', str(root), 'show', CHECKPOINT + ':' + path], text=True)


def check(current, header):
    USED.clear()
    old = reference('src/core/SisProdTransient.cpp')
    old_header = reference('src/include/SisProdTransientInternal.h')
    first = 'enum class TimeStepHook'
    last = 'PressureRateOfChangePolicy>;'
    equal('policy context, hooks, actions and registry',
          old_header[old_header.index(first):old_header.index(last)+len(last)],
          header[header.index(first):header.index(last)+len(last)])
    assert canonical('struct FacePhaseDensities { double liquid; double gas; };') in canonical(header)
    assert len(re.findall(r'^namespace \{', current, re.M)) == 1
    names = re.findall(r'(?m)^[\w:<>]+\s+(\w+)\s*\(', masked(old))
    for name in names:
        old_signature, old_body = parts(old, name)
        signature, _ = parts(current, name)
        equal(name + ' signature', old_signature, signature)
        body = check_valves(old, current) if name == 'restrictTimeStepByValve' else expand(current, name)
        if name == 'applyOutletBufferCondition':
            # Both spellings multiply double flows; integer 0/1 undergo the same
            # conversion to double as the explicitly floating constants.
            expression = '(state.cells[state.lastCell - 1].MRBuf - state.cells[state.lastCell - 1].MliqiniRBuf) * 0.0 + 1.0 * state.cells[state.lastCell - 1].MliqiniRBuf'
            normalized = expression.replace('0.0', '0').replace('1.0', '1')
            old_body = canonical(old_body).replace(canonical(expression), canonical(normalized))
        equal(name + ' calculation', old_body, body)
        print(f'OK {name}: {len(canonical(body).split())} tokens')
    assert USED == set(SIGNATURES), ('unchecked helpers', set(SIGNATURES)-USED, USED-set(SIGNATURES))
    print(f'{len(names)} source functions reconstructed; {len(USED)} helper interfaces and bindings checked.')
    return 0


if __name__ == '__main__':
    root = pathlib.Path(__file__).resolve().parents[1]
    source = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else root / 'src/core/SisProdTransient.cpp'
    header = pathlib.Path(sys.argv[2]) if len(sys.argv) > 2 else root / 'src/include/SisProdTransientInternal.h'
    sys.exit(check(source.read_text(), header.read_text()))
