"""Literal reconstruction of the thermal mass-transfer extractions.

The checkpoint follows the scalar/interface cleanup. Comparisons retain every
formula, state write and call in the distributed models and the thermal source.
Only checked immutable value bindings and per-face control scaffolding expand.
"""
import pathlib
import re
import subprocess

CHECKPOINT = '339df28'

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from module_dependencies_layout import verified_snapshot, candidate_text



def canonical(text):
    text = re.sub(r'//[^\n]*|/\*[\s\S]*?\*/', '', text)
    return ' '.join(re.findall(r'[A-Za-z_]\w*|\d+(?:\.\d*)?(?:[eE][+-]?\d+)?|\S', text))


def parts(text, name):
    match = re.search(r'(?m)^(?:\[\[nodiscard\]\] )?[\w:<>]+\s+' + name + r'\s*\(', text)
    assert match, name + ' definition'
    opening = text.index('{', match.end())
    clean = re.sub(r'//[^\n]*|/\*[\s\S]*?\*/', lambda m: ' ' * len(m[0]), text)
    depth, end = 1, opening + 1
    while depth:
        depth += (clean[end] == '{') - (clean[end] == '}')
        end += 1
    return canonical(text[match.start():opening]), canonical(text[opening+1:end-1])


def replace(text, old, new):
    old, new = canonical(old), canonical(new)
    assert text.count(old) == 1, (old, text.count(old))
    return text.replace(old, new)


def reference():
    root = pathlib.Path(__file__).resolve().parents[1]
    return subprocess.check_output(['git', '-C', str(root), 'show', CHECKPOINT + ':src/core/SisProdThermal.cpp'], text=True)


def equal(label, old, new):
    assert old == new, label + ' differs near ' + next((f'{a!r} / {b!r}' for a, b in zip(old.split(), new.split()) if a != b), 'end')
    print(f'OK       {label}: {len(old.split())} tokens')


CELL_REFS = 'Cel &cell = state.cells[cellIndex]; Cel &leftCell = state.cells[cellIndex - 1];'
PREVIOUS = {
    'leftCellLeftFaceMixtureLiquidDensity': 'liquidDensity',
    'leftCellLeftFaceOilVolumeFactor': 'oilVolumeFactor',
    'leftCellLeftFaceSolutionGasRatio': 'solutionGasRatio',
    'leftCellLeftFaceSolutionGasPressureDerivative': 'solutionGasPressureDerivative',
}
COEFFICIENTS = {'activeDerivative':'activeDerivative', 'spatialCoupling':'spatialCoupling', 'coefficientFlowArea':'flowArea'}


def split_aliases(body):
    first = body.index(canonical(CELL_REFS))
    prefix, body = body[:first], body[first:]
    aliases = {}
    pattern = r'const double (\w+) = ([^;]+) ;\s*'
    for match in re.finditer(pattern, prefix):
        assert match[1] not in aliases
        aliases[match[1]] = match[2].strip()
    assert canonical(prefix) == canonical(' '.join('const double ' + k + ' = ' + v + ';' for k,v in aliases.items())), 'alias effects'
    return aliases, replace(body, CELL_REFS, '')


def check_mass(current):
    current = verified_snapshot({'src/core/SisProdThermal.cpp': current})['src/core/SisProdThermal.cpp']
    old = reference()
    rebuilt = {name: parts(current, name) for name in (
        'initializeDistributedMassTransferInlet', 'prepareDistributedMassTransferProperties',
        'lastCellSolutionGasPressureDerivative', 'updateDistributedMassTransferDerivatives',
        'selectDistributedMassTransferModel', 'snapshotMassTransferState',
        'applyDistributedMassTransferModel', 'clearMassTransferDerivatives')}
    sig, donor = parts(current, 'prepareMassTransferDonors')
    equal('donor input binding', canonical('MassTransferWaterCuts prepareMassTransferDonors(const ThermalState &state, int cellIndex, bool leftCellHasNegativeLiquidMassFlowRate, ProFlu &leftCellLeftFaceFluid, ProFlu &leftFaceFluid)'), sig)
    donor = replace(donor, CELL_REFS, '')
    donor = replace(donor, 'return MassTransferWaterCuts{leftFaceWaterCut, leftCellLeftFaceWaterCut, leftCellWaterCut};', '')
    sig, body = rebuilt['prepareDistributedMassTransferProperties']
    plumbing = ('const MassTransferWaterCuts waterCuts = prepareMassTransferDonors(state, cellIndex, leftCellHasNegativeLiquidMassFlowRate, leftCellLeftFaceFluid, leftFaceFluid);'
                'const double leftFaceWaterCut = waterCuts.leftFaceWaterCut;'
                'const double leftCellLeftFaceWaterCut = waterCuts.leftCellLeftFaceWaterCut;'
                'const double leftCellWaterCut = waterCuts.leftCellWaterCut;')
    rebuilt['prepareDistributedMassTransferProperties'] = sig, replace(body, plumbing, donor)

    sig, derivatives = parts(current, 'updateEquilibriumMassTransferDerivatives')
    equal('model-0 input binding', canonical('void updateEquilibriumMassTransferDerivatives(const ThermalState &state, int cellIndex, const ProFlu &leftCellLeftFaceFluid, const ProFlu &leftFaceFluid, const DistributedMassTransferProperties &properties, const DistributedMassTransferCoefficients &coefficients, const MassTransferFaceState &previousFace)'), sig)
    aliases, derivatives = split_aliases(derivatives)
    value_aliases = {k: canonical('coefficients.' + v) for k,v in COEFFICIENTS.items()}
    value_aliases.update({k:canonical('previousFace.' + v) for k,v in PREVIOUS.items()})
    for name, expression in value_aliases.items():
        assert aliases[name] == expression, name + ' value binding'
    current_sig, application = rebuilt['applyDistributedMassTransferModel']
    equal('model input binding', canonical('void applyDistributedMassTransferModel(const ThermalState &state, int cellIndex, double leftAbsoluteSuperficialVelocity, const ProFlu &leftCellLeftFaceFluid, const ProFlu &leftFaceFluid, const DistributedMassTransferProperties &properties, const DistributedMassTransferCoefficients &coefficients, const MassTransferFaceState &previousFace)'), current_sig)
    parent_aliases, application = split_aliases(application)
    assert all(aliases.get(k) == v for k,v in parent_aliases.items()), 'model alias binding'
    application = replace(application,
        'updateEquilibriumMassTransferDerivatives(state, cellIndex, leftCellLeftFaceFluid, leftFaceFluid, properties, coefficients, previousFace);', derivatives)
    guard = canonical('if (!(leftCell.TMModel == 0 || leftCell.TMModel == 1)) return;')
    assert application.count(guard) == 1
    before, after = application.split(guard)
    application = before + canonical('if (leftCell.TMModel == 0 || leftCell.TMModel == 1) {') + ' ' + after + ' }'
    prefix = ' '.join('const double ' + k + ' = ' + v + ';' for k,v in aliases.items() if k not in value_aliases)
    rebuilt['applyDistributedMassTransferModel'] = parts(old, 'applyDistributedMassTransferModel')[0], canonical(prefix + CELL_REFS) + ' ' + application

    # Check the sweep shape independently before expanding a cell into the
    # reference loop shell. Negative lastCell executes no inlet in either shape.
    equal('sweep shape', canonical('void updateDistributedMassTransfer(const ThermalState &state) { if (state.lastCell < 0) return; MassTransferFaceState previousFace = initializeDistributedMassTransferInlet(state, 0); for (int cellIndex = 1; cellIndex <= state.lastCell; cellIndex++) { updateDistributedMassTransferCell(state, cellIndex, previousFace); } }'), canonical(' '.join(parts(current, 'updateDistributedMassTransfer')[:1]) + ' { ' + parts(current, 'updateDistributedMassTransfer')[1] + ' }'))
    sig, cell = parts(current, 'updateDistributedMassTransferCell')
    equal('sweep carry binding', canonical('void updateDistributedMassTransferCell(const ThermalState &state, int cellIndex, MassTransferFaceState &previousFace)'), sig)
    cell = replace(cell, CELL_REFS + ' Cel &rightCell = state.cells[cellIndex + 1];', '')
    sig, final = parts(current, 'finalizeDistributedMassTransferFace')
    equal('face limits input binding', canonical('void finalizeDistributedMassTransferFace(const ThermalState &state, int cellIndex, const DistributedMassTransferProperties &properties, const DistributedMassTransferCoefficients &coefficients)'), sig)
    final_prefix = (CELL_REFS + 'const double betI = properties.leftFaceComposition;'
                    'const double gasDensity = properties.gasDensity;'
                    'const double liquidDensity = properties.liquidDensity;'
                    'const double coefficientFlowArea = coefficients.flowArea;')
    final = replace(final, final_prefix, '')
    cell = replace(cell, 'finalizeDistributedMassTransferFace(state, cellIndex, properties, coefficients);', final)
    cell = replace(cell, 'double leftCellWaterCut = properties.leftCellWaterCut;',
                   'double leftCellWaterCut = properties.leftCellWaterCut; double liquidDensity = properties.liquidDensity; double gasDensity = properties.gasDensity; double betI = properties.leftFaceComposition;')
    new_call = ('applyDistributedMassTransferModel(state, cellIndex, leftAbsoluteSuperficialVelocity, leftCellLeftFaceFluid, leftFaceFluid, properties, coefficients, previousFace);')
    old_call = ('applyDistributedMassTransferModel(state, cellIndex, leftAbsoluteSuperficialVelocity, leftCellLeftFaceFluid, leftFaceFluid, properties, ' + ','.join((*COEFFICIENTS, *PREVIOUS)) + ');')
    cell = replace(cell, new_call, ' '.join('double ' + k + ' = coefficients.' + v + ';' for k,v in COEFFICIENTS.items()) + old_call)
    for variable, field in PREVIOUS.items():
        cell = cell.replace(canonical('previousFace.' + field), variable)
    old_sig, old_main = parts(old, 'updateDistributedMassTransfer')
    branch = canonical('if (cellIndex != 0 && cellIndex != state.lastCell + 1) {')
    start = old_main.index(branch) + len(branch)
    end = old_main.index(canonical('} else if (cellIndex == 0)'), start)
    # The surrounding reference loop is used only after the current public
    # coordinator and all four rolling field bindings have been checked above.
    rebuilt['updateDistributedMassTransfer'] = old_sig, old_main[:start] + ' ' + cell + ' ' + old_main[end:]
    for name, (sig, body) in rebuilt.items():
        sig = sig.replace('MassTransferFaceState', 'InletMassTransferSeed')
        body = body.replace('MassTransferFaceState', 'InletMassTransferSeed')
        expected_sig, expected_body = parts(old, name)
        equal(name + ' signature', expected_sig, sig)
        equal(name, expected_body, canonical(body))
    return 0


def check_thermal_source(current):
    current = verified_snapshot({'src/core/SisProdThermal.cpp': current})['src/core/SisProdThermal.cpp']
    old = reference()
    sig, kinetic = parts(current, 'computeThermalMassTransferKineticTerm')
    equal('thermal kinetic inputs', canonical('double computeThermalMassTransferKineticTerm(const ThermalState &state, int cellIndex, double meanSuperficialGasVelocity, double meanSuperficialLiquidVelocity)'), sig)
    kinetic = replace(kinetic, CELL_REFS, '')
    kinetic = replace(kinetic, 'return kineticTerm;', '')
    sig, body = parts(current, 'computeThermalMassTransfer')
    body = replace(body, 'const double kineticTerm = computeThermalMassTransferKineticTerm(state, cellIndex, meanSuperficialGasVelocity, meanSuperficialLiquidVelocity);', kinetic)
    equal('computeThermalMassTransfer', ' '.join(parts(old, 'computeThermalMassTransfer')), sig + ' ' + body)
    equal('computeThermalMassTransferSourceTerms', ' '.join(parts(old, 'computeThermalMassTransferSourceTerms')), ' '.join(parts(current, 'computeThermalMassTransferSourceTerms')))
    return 0
