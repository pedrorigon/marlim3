"""Find old EN key names used as string literals in Python source."""
from pathlib import Path

RENAMES = {
    'productionDucts': 'productionPipe', 'serviceDucts': 'servicePipe',
    'crossProductionProfiles': 'crossProductionProfile', 'crossServiceProfiles': 'crossServiceProfile',
    'deadOilTemps': 'deadOilTemp', 'deadOilViscs': 'deadOilVisc',
    'layers': 'layer', 'manufacturerStages': 'manufacturerStage',
    'measuredPositions': 'measuredPosition', 'productionFluids': 'productionFluid',
    'sonicFlags': 'sonicFlag', 'sonicTimes': 'sonicTime', 'stages': 'stage',
    'temperatures': 'temperatureArray', 'times': 'timeArray',
    'injectionPressures': 'injectionPressureArray',
    'maxTimestep': 'maxDT', 'bswCut': 'bswInversion',
    'blackBoxCorrelation': 'ssFlowModel',
    'hydroGradient': 'correctionFactorHydro', 'fricGradient': 'correctionFactorFric',
    'tempGradient': 'correctionFactorTemp', 'annularColumnFlag': 'annulusColumnFlag',
    'recoveryFreq': 'recoveryFactor', 'liquidRecoveryFreq': 'liquidRecoveryFactor',
    'areaGap': 'gapArea', 'pigDischCoef': 'pigCD',
    'thermalEquilBcs': 'thermalEquilEsp', 'initialFlowGuess': 'initialFlowRateGuess',
    'tubeColumn': 'tubing', 'frictionPressure': 'frictionPressureGradient',
    'emulsionVec': 'emulsionVect', 'superficialLiquidVel': 'usl',
    'superficialGasVel': 'usg', 'liquidVelocity': 'ul', 'gasVelocity': 'ug',
    'externalTemp': 'ambientTemp', 'externalVel': 'ambientVel',
    'externalConductivity': 'ambientConductivity', 'externalSpecificHeat': 'ambientSpecificHeat',
    'externalDensity': 'ambientDensity', 'externalVisc': 'ambientVisc',
    'externalEnvironment': 'ambientEnvironment', 'bet': 'complementaryFluidFraction',
}

search_dirs = [Path('tests'), Path('marlim3')]
skip_files = {'_keys.py', 'apply_renames.py', 'find_old_keys.py'}

for d in search_dirs:
    for f in d.rglob('*.py'):
        if f.name in skip_files:
            continue
        content = f.read_text(encoding='utf-8')
        for old, new in RENAMES.items():
            for pattern in [f'"{old}"', f"'{old}'"]:
                if pattern in content:
                    lines = content.split('\n')
                    for i, line in enumerate(lines, 1):
                        if pattern in line and 'ATTR_TO_JSON' not in line:
                            print(f'{f.relative_to(".")}:{i}: {old} -> {new}: {line.strip()[:100]}')
