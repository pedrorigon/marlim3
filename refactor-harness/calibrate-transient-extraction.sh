#!/usr/bin/env bash
# Calibrate source reconstruction without executing the simulator.
set -euo pipefail
task_root="$(cd "$(dirname "$0")/.." && pwd)"
python3 - "$task_root" <<'PY'
import contextlib
import importlib.util
import io
import json
import os
import pathlib
import sys

root = pathlib.Path(sys.argv[1])
spec = importlib.util.spec_from_file_location('transient_proof', root / 'refactor-harness/transient-extraction-layout.py')
proof = importlib.util.module_from_spec(spec)
spec.loader.exec_module(proof)
source = (root / 'src/core/SisProdTransient.cpp').read_text()
header = (root / 'src/include/SisProdTransientInternal.h').read_text()

def check(candidate, declarations):
    with contextlib.redirect_stdout(io.StringIO()):
        return proof.check(candidate, declarations)

assert check(source, header) == 0
assert check(source.replace('namespace {\n', 'namespace {\n// Private implementation.\n', 1), header) == 0
cases = [
    ('upstream-pressure-area', 'updateInteriorCell',
     'updateInteriorAuxiliaryPressure(state, i)', 'updateInteriorAuxiliaryPressure(state, i + 1)'),
    ('density-result-binding', 'updateInteriorVolumeFlows',
     'double rhol = densities.liquid;', 'double rhol = densities.gas;'),
    ('outlet-donor-weight', 'updateOutletVolumeFlows', '* 0.99 + 0.01 *', '* 0.01 + 0.99 *'),
    ('nominal-perturbed-choke-step', 'applyOutletPressureCondition',
     'pmon = state.outletPressure * 1.0001;', 'pmon = state.outletPressure * 1.001;'),
    ('buffered-reverse-fraction', 'applyOutletBufferCondition',
     'state.surfaceChokeOpen = 1;\n        updateBufferedNetworkSources(state, betRev);',
     'state.surfaceChokeOpen = 1;\n        updateBufferedNetworkSources(state, alfRev);'),
    ('implicit-phase-order', 'computeImplicitTimeStep',
     'applyChokeTimeStepControls(state);\n    restrictTimeStepByTransport(state);',
     'restrictTimeStepByTransport(state);\n    applyChokeTimeStepControls(state);'),
    ('valve-first-candidate', 'restrictTimeStepByValve',
     'valveRestrictedTimeStep(state, 0, dtaux)', 'valveRestrictedTimeStep(state, 1, dtaux)'),
    ('valve-strict-minimum', 'restrictTimeStepByValve',
     'if (candidate < minimumTimeStep)', 'if (candidate <= minimumTimeStep)'),
    ('valve-shared-base', 'restrictTimeStepByValve',
     'valveRestrictedTimeStep(state, i, dtaux)', 'valveRestrictedTimeStep(state, i, minimumTimeStep)'),
    ('valve-artificial-gas-order', 'valveRestrictedTimeStep',
     'state.desperationMoment += 1.5;', 'state.desperationMoment += 1.0;'),
    ('fraction-sweep-outlet', 'evolveFractions',
     'updateOutletVoidFraction(state, state.lastCell, alfrev);',
     'updateOutletVoidFraction(state, state.lastCell, betrev);'),
    ('restart-cell-reference', 'restartFractionEvolutionInitial',
     'const Cel &cell = state.cells[i];', 'const Cel &cell = state.cells[i + 1];'),
    ('restart-limit-field', 'restartFractionEvolutionInitial',
     'if (cell.dt1 < state.timeStep)', 'if (cell.dt2 < state.timeStep)'),
    ('restart-strict-minimum', 'restartFractionEvolutionInitial',
     'if (cell.dt < state.timeStep)', 'if (cell.dt <= state.timeStep)'),
    ('restart-outlet-bound', 'restartFractionEvolutionInitial',
     'for (int i = 0; i <= state.lastCell; i++)', 'for (int i = 0; i < state.lastCell; i++)'),
    ('coupling-restore-binding', 'advanceCouplingIteration',
     'restoreIntermediateCouplingState(state, kontaAcop, gasMassSource, liquidMassSource, complementaryMassSource);',
     'restoreIntermediateCouplingState(state, kontaAcop, liquidMassSource, gasMassSource, complementaryMassSource);'),
    ('accepted-step-history-order', 'solveTransientStep',
     'finalizeTransientCellState(state, dependencies, fluiRev);\n\n    updateTransientMovingMeans(state);',
     'updateTransientMovingMeans(state);\n\n    finalizeTransientCellState(state, dependencies, fluiRev);'),
    ('corner-order', 'fillMiniTableCornersAtMinPressure',
     'fillMiniTableCorner(state, fluC, flui, 0, 1);', 'fillMiniTableCorner(state, fluC, flui, 1, 0);'),
    ('corner-property', 'fillMiniTableCorner',
     'fluC.drhodp(pressure, temperature);', 'fluC.drhodt(pressure, temperature);'),
    ('absent-phase-donor', 'supplyAbsentGasProperties',
     'titVec[busca].second / 2', 'titVec[busca].second % 2'),
    ('accessory-bound', 'refreshAccessoryFluidMiniTable',
     'accessoryFluid.miniTabDin.pmin = state.cells[i].flui.miniTabDin.pmin;',
     'accessoryFluid.miniTabDin.pmin = state.cells[i].flui.miniTabDin.pmax;'),
    ('pig-receiver-loop', 'updateActivePigState',
     'state.movingPigCount--;', 'state.movingPigCount--; break;'),
]
records = []
for label, method, before, after in cases:
    signature, body = proof.parts(source, method)
    assert body.count(before) == 1, (label, body.count(before))
    start = source.index(signature)
    opening = source.index('{', start + len(signature))
    finish = proof.closing(source, opening, '{', '}') + 1
    definition = source[start:finish]
    assert definition.count(before) == 1, label
    candidate = source[:start] + definition.replace(before, after) + source[finish:]
    try:
        detected = check(candidate, header) != 0
    except (AssertionError, ValueError):
        detected = True
    assert detected, label
    records.append({'case': label, 'injected_once': True, 'detected': detected})

for label, before, after in [
    ('policy-condition', 'return context.couplingIteration == 0 && context.step.input.controleDTvalv == 1;',
     'return context.couplingIteration == 1 && context.step.input.controleDTvalv == 1;'),
    ('policy-hook', 'hook = TimeStepHook::AfterPigUpdate;', 'hook = TimeStepHook::AfterValveOpenings;'),
    ('density-layout', 'double liquid;\n    double gas;', 'double gas;\n    double liquid;'),
]:
    assert header.count(before) == 1, label
    try:
        detected = check(source, header.replace(before, after)) != 0
    except (AssertionError, ValueError):
        detected = True
    assert detected, label
    records.append({'case': label, 'injected_once': True, 'detected': detected})

assert (root / 'src/core/SisProdTransient.cpp').read_text() == source
assert (root / 'src/include/SisProdTransientInternal.h').read_text() == header
report = os.environ.get('MARLIM_TRANSIENT_STATIC_REPORT')
if report:
    pathlib.Path(report).write_text(json.dumps({
        'kind': 'static source reconstruction; no simulator execution',
        'positive_controls_passed': True, 'cases': records, 'working_source_untouched': True,
    }, indent=2) + '\n')
print(f'TRANSIENT CALIBRATION PASSED: {len(records)} mutations detected; controls pass; source unchanged.')
PY
