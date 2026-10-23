#!/usr/bin/env python3
"""Reconstruct fluid-tracking calculations against the reviewed source checkpoint.

Checks every prior definition, private type, public contract, new helper interface,
and donor binding. Only the two relocated trivial declarations, baseline scalar
balance defaults and equivalent integer-one promotion are normalized explicitly.
No simulator is executed and no numerical reference is generated.
"""
import pathlib
import re
import subprocess
import sys

CHECKPOINT = 'f71be8c'
SIGNATURES = {'computePhaseLiquidVolumeBalances': 'void computePhaseLiquidVolumeBalances(const FluidTrackingState '
                                     '&state, const PhaseFace &left, const PhaseFace &right, const '
                                     'PhaseSource &source, PhaseBalance &balance, double dx, double '
                                     'flowArea, int i, double dt)',
 'copyMixtureFaceComposition': 'void copyMixtureFaceComposition(PhaseFace &face, const ProFlu &fluid, '
                               'Vcr<double> &fractions, int ncomp)',
 'prepareBlackOilAccessorySource': 'void prepareBlackOilAccessorySource(const FluidTrackingState '
                                   '&state, BlackOilSource &source, ProFlu &fluF, int i, double '
                                   'temperatureHigh, double temperatureLow)',
 'prepareBlackOilCellSourceBasis': 'void prepareBlackOilCellSourceBasis(const FluidTrackingState '
                                   '&state, BlackOilCell &cell, BlackOilSource &source, int i, double '
                                   'temperatureHigh, double temperatureLow)',
 'prepareBlackOilFaces': 'void prepareBlackOilFaces(const FluidTrackingState &state, BlackOilFace '
                         '&left, BlackOilFace &right, int i, double temperatureHigh, double '
                         'temperatureLow)',
 'prepareBlackOilStockTankBalances': 'void prepareBlackOilStockTankBalances(const FluidTrackingState '
                                     '&state, const BlackOilFace &left, const BlackOilFace &right, '
                                     'const BlackOilCell &cell, const BlackOilSource &source, '
                                     'BlackOilBalance &balance, double dx, double flowArea, int i, '
                                     'double dt)',
 'prepareLiquidFaceTransport': 'void prepareLiquidFaceTransport(PhaseFace &face, const ProFlu &fluid, '
                               'const ProFluCol &complementaryFluid, const double &pressure, const '
                               'double &temperature, const double &liquidVolumeFlow, const double '
                               '&voidFraction, Vcr<double> &fractions, int ncomp)',
 'preparePhaseFacesAndCellWeight': 'void preparePhaseFacesAndCellWeight(const FluidTrackingState '
                                   '&state, PhaseFace &left, PhaseFace &right, PhaseCell &cell, int i, '
                                   'double temperatureHigh, double temperatureLow, int ncomp)',
 'preparePhaseMixtureSource': 'void preparePhaseMixtureSource(const FluidTrackingState &state, '
                              'PhaseSource &source, ProFlu &fluF, int i, int ncomp)',
 'prepareVapourFaceComposition': 'void prepareVapourFaceComposition(PhaseFace &face, const ProFlu '
                                 '&fluid, const double &pressure, const double &temperature, const '
                                 'double &voidFraction, Vcr<double> &fractions, int ncomp)',
 'publishBlackOilCellProperties': 'void publishBlackOilCellProperties(const FluidTrackingState &state, '
                                  'int i, Vcr<double> &rgo, Vcr<double> &dg, Vcr<double> &yco2, '
                                  'Vcr<double> &API, Vcr<double> &BSW, Vcr<double> &denag, Vcr<double> '
                                  '&VISCL, Vcr<double> &VISCH, double temperatureHigh, double '
                                  'temperatureLow)',
 'publishPhaseCellProperties': 'void publishPhaseCellProperties(const FluidTrackingState &state, ProFlu '
                               '*fluC, int i, Vcr<double> &BSW, Vcr<double> &denag, Vcr<double> &VISCL, '
                               'Vcr<double> &VISCH, double temperatureHigh, double temperatureLow)',
 'readPhaseSourceLiquidQuality': 'void readPhaseSourceLiquidQuality(const FluidTrackingState &state, '
                                 'PhaseSource &source, const ProFlu &fluF, int i)',
 'refreshPhaseSourceComposition': 'void refreshPhaseSourceComposition(const FluidTrackingState &state, '
                                  'ProFlu &fluid, const double &pressure, const double &temperature)',
 'transportPhaseLiquidProperties': 'void transportPhaseLiquidProperties(const FluidTrackingState '
                                   '&state, const PhaseFace &left, const PhaseFace &right, PhaseCell '
                                   '&cell, PhaseSource &source, ProFlu &fluF, int i, Vcr<double> &BSW, '
                                   'Vcr<double> &denag, Vcr<double> &VISCL, Vcr<double> &VISCH, double '
                                   'temperatureHigh, double temperatureLow, double dx, double flowArea, '
                                   'double dt)',
 'updateBlackOilOutletProperties': 'void updateBlackOilOutletProperties(const FluidTrackingState '
                                   '&state, const ProFlu &fluiRev, Vcr<double> &rgo, Vcr<double> &dg, '
                                   'Vcr<double> &yco2, Vcr<double> &API, Vcr<double> &BSW, Vcr<double> '
                                   '&denag, Vcr<double> &VISCL, Vcr<double> &VISCH, double '
                                   'temperatureHigh, double temperatureLow)',
 'updatePhaseMoleInventory': 'void updatePhaseMoleInventory(const FluidTrackingState &state, const '
                             'PhaseCell &cell, ProFlu *fluC, int i)',
 'updatePhaseOutletProperties': 'void updatePhaseOutletProperties(const FluidTrackingState &state, '
                                'const ProFlu &fluiRev, Vcr<double> &BSW, Vcr<double> &denag, '
                                'Vcr<double> &VISCL, Vcr<double> &VISCH, double temperatureHigh, double '
                                'temperatureLow)'}



def masked(text):
    pattern = r'//[^\n]*|/\*[\s\S]*?\*/|"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\''
    return re.sub(pattern, lambda m: re.sub(r'[^\n]', ' ', m[0]), text)


def tokens(text):
    pattern = r'//[^\n]*|/\*[\s\S]*?\*/|"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\'|[A-Za-z_]\w*|\d+(?:\.\d*)?(?:[eE][+-]?\d+)?|\S'
    return [m[0] for m in re.finditer(pattern, text) if not m[0].startswith(('//', '/*'))]


def canonical(text):
    return ' '.join(tokens(text))


def closing(text, opening, left, right):
    clean, depth = masked(text), 0
    for pos in range(opening, len(clean)):
        depth += (clean[pos] == left) - (clean[pos] == right)
        if depth == 0:
            return pos
    raise AssertionError('unbalanced ' + left)


def parts(source, name):
    matches = list(re.finditer(r'(?m)^[\w:<>]+\s+' + re.escape(name) + r'\s*\(', masked(source)))
    assert len(matches) == 1, (name, len(matches))
    first = matches[0].start()
    opening = source.index('{', matches[0].end())
    end = closing(source, opening, '{', '}')
    return source[first:opening].strip(), source[opening+1:end]


def equal(label, expected, actual):
    expected, actual = tokens(expected), tokens(actual)
    if expected != actual:
        pos = next((i for i, pair in enumerate(zip(expected, actual)) if pair[0] != pair[1]),
                   min(len(expected), len(actual)))
        raise AssertionError(f'{label}: token {pos}\n'
                             f'expected: {" ".join(expected[max(0,pos-6):pos+12])}\n'
                             f'actual:   {" ".join(actual[max(0,pos-6):pos+12])}')


def replace_once(text, before, after):
    values, wanted = tokens(text), tokens(before)
    starts = [i for i in range(len(values)-len(wanted)+1) if values[i:i+len(wanted)] == wanted]
    assert len(starts) == 1, ('replacement', before, len(starts))
    i = starts[0]
    values[i:i+len(wanted)] = tokens(after)
    return ' '.join(values)


def arguments(text):
    result, start, stack = [], 0, []
    for i, char in enumerate(text):
        if char in '([':
            stack.append(char)
        elif char in ')]':
            assert stack, text
            stack.pop()
        elif char == ',' and not stack:
            result.append(text[start:i].strip())
            start = i+1
    result.append(text[start:].strip())
    return result


def params(signature):
    inner = signature[signature.index('(')+1:signature.rindex(')')]
    return [re.findall(r'[A-Za-z_]\w*', item)[-1] for item in arguments(inner)]


USED = set()
DEFAULTS = ''.join('balance.' + field + ' = 0.;' for field in
                   ['MultPe', 'MultPd', 'residuoP', 'MultAe', 'MultAd'])


def helper(source, name):
    signature, body = parts(source, name)
    equal(name + ' interface', SIGNATURES[name], signature)
    USED.add(name)
    if name == 'prepareBlackOilFaces':
        body = replace_once(body, 'right.betI = state.cells[i].betPigD;',
                            'BlackOilFace right; right.betI = state.cells[i].betPigD;')
    elif name == 'prepareBlackOilCellSourceBasis':
        body = replace_once(body, 'source.dissolvedGas = state.cells[i].fontemassLR;',
                            'BlackOilSource source; source.dissolvedGas = state.cells[i].fontemassLR;')
    elif name == 'prepareBlackOilStockTankBalances':
        # Exactly the five baseline stores are required at the extracted boundary.
        tail = 'balance.volaguaFim = cell.liquidHoldup * (1 - cell.complementaryFraction) * cell.bsw;'
        body = replace_once(body, tail + DEFAULTS, tail)
    return signature, body


def expand(source, body):
    for _ in range(128):
        candidates = list(re.finditer(r'\b(' + '|'.join(SIGNATURES) + r')\s*\(', masked(body)))
        if not candidates:
            return body
        call = candidates[0]
        name = call[1]
        opening = body.index('(', call.start())
        end = closing(body, opening, '(', ')')
        terminator = end+1
        while body[terminator].isspace():
            terminator += 1
        assert body[terminator] == ';', name
        signature, core = helper(source, name)
        actual = arguments(body[opening+1:end])
        formal = params(signature)
        assert len(actual) == len(formal), (name, actual, formal)
        bindings = dict(zip(formal, actual))
        core = re.sub(r'\b(' + '|'.join(formal) + r')\b', lambda m: bindings[m[0]], core)
        body = body[:call.start()] + core + body[terminator+1:]
    raise AssertionError('cyclic or excessive helper expansion')


def main():
    root = pathlib.Path(__file__).resolve().parent.parent
    cpp = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else root/'src/core/SisProdFluidTracking.cpp'
    header = pathlib.Path(sys.argv[2]) if len(sys.argv) > 2 else root/'src/include/SisProdFluidTrackingInternal.h'
    public = pathlib.Path(sys.argv[3]) if len(sys.argv) > 3 else root/'src/include/SisProdFluidTracking.h'
    old = subprocess.check_output(['git', 'show', CHECKPOINT+':src/core/SisProdFluidTracking.cpp'],
                                  cwd=root, text=True)
    old_header = subprocess.check_output(['git', 'show', CHECKPOINT+':src/include/SisProdFluidTrackingInternal.h'],
                                         cwd=root, text=True)
    old_public = subprocess.check_output(['git', 'show', CHECKPOINT+':src/include/SisProdFluidTracking.h'],
                                         cwd=root, text=True)
    new = cpp.read_text()
    names = re.findall(r'(?m)^[\w:<>]+\s+(\w+)\s*\(', masked(old))
    new_names = re.findall(r'(?m)^[\w:<>]+\s+(\w+)\s*\(', masked(new))
    assert len(names) == 33
    assert sorted(new_names) == sorted(names + list(SIGNATURES)), 'definition set changed'
    for name in names:
        signature, expected = parts(old, name)
        current_signature, actual = parts(new, name)
        equal(name + ' signature', signature, current_signature)
        if name == 'transportCellBlackOilProperties':
            actual = replace_once(actual, 'BlackOilFace left; BlackOilFace right;', 'BlackOilFace left;')
            actual = replace_once(actual, 'BlackOilCell cell; BlackOilSource source;', 'BlackOilCell cell;')
        actual = expand(new, actual)
        if name == 'upwindLiquidFractions':
            # With a double right operand, integer one is represented exactly as double.
            expected = canonical(expected).replace('> 1 - 1e-3', '> 1. - 1e-3')
        equal(name, expected, actual)
    assert USED == set(SIGNATURES), ('unused helper', set(SIGNATURES)-USED)
    current_header = header.read_text()
    structures = re.findall(r'(?m)^struct (\w+) \{', old_header)
    assert len(structures) == 9
    assert re.findall(r'(?m)^struct (\w+) \{', current_header) == structures
    for name in structures:
        pattern = r'(?m)^struct ' + name + r' \{[\s\S]*?^};'
        equal(name + ' layout', re.search(pattern, old_header)[0], re.search(pattern, current_header)[0])
    for name, signature in SIGNATURES.items():
        declaration = canonical(signature + ';')
        assert canonical(current_header).count(declaration) == 1, name + ' declaration'
    equal('public header contract', old_public, public.read_text())
    # The defaults predate the cleanup and retain their original scalar values.
    baseline = subprocess.check_output(['git', 'show', 'beff16e:src/core/SisProdFluidTracking.cpp'],
                                       cwd=root, text=True)
    assert canonical(DEFAULTS) in canonical(parts(baseline, 'transportCellBlackOilProperties')[1])
    print(f'{len(names)} source functions reconstructed; {len(SIGNATURES)} helper interfaces and bindings checked.')
    print('9 private layouts and the public state, callbacks and transport contracts preserved.')


if __name__ == '__main__':
    main()
