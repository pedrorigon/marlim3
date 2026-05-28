"""Apply batch EN key renames to translations.json and English demo files.

This script is the single place defining all renames. Run once, then delete.
"""
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

# ==========================================================================
# RENAME MAP: old EN key -> new EN key
# ==========================================================================
RENAMES = {
    # --- Duct → Pipe + plural → singular ---
    "productionDucts": "productionPipe",
    "serviceDucts": "servicePipe",

    # --- Plural → Singular (no conflicts) ---
    "crossProductionProfiles": "crossProductionProfile",
    "crossServiceProfiles": "crossServiceProfile",
    "deadOilTemps": "deadOilTemp",
    "deadOilViscs": "deadOilVisc",
    "layers": "layer",
    "manufacturerStages": "manufacturerStage",
    "measuredPositions": "measuredPosition",
    "productionFluids": "productionFluid",
    "sonicFlags": "sonicFlag",
    "sonicTimes": "sonicTime",
    "stages": "stage",

    # --- Plural → Array suffix (conflicts with existing singular) ---
    "temperatures": "temperatureArray",
    "times": "timeArray",
    "injectionPressures": "injectionPressureArray",

    # --- Specific renames from review ---
    "bet": "complementaryFluidFraction",
    "maxTimestep": "maxDT",
    "bswCut": "bswInversion",
    "blackBoxCorrelation": "ssFlowModel",
    "hydroGradient": "correctionFactorHydro",
    "fricGradient": "correctionFactorFric",
    "tempGradient": "correctionFactorTemp",
    "annularColumnFlag": "annulusColumnFlag",
    "recoveryFreq": "recoveryFactor",
    "liquidRecoveryFreq": "liquidRecoveryFactor",
    "areaGap": "gapArea",
    "pigDischCoef": "pigCD",
    "thermalEquilBcs": "thermalEquilEsp",
    "initialFlowGuess": "initialFlowRateGuess",
    "tubeColumn": "tubing",
    "frictionPressure": "frictionPressureGradient",
    "emulsionVec": "emulsionVect",

    # --- Abbreviations for velocities ---
    "superficialLiquidVel": "usl",
    "superficialGasVel": "usg",
    "liquidVelocity": "ul",
    "gasVelocity": "ug",

    # --- ambient vs external (medium properties → ambient) ---
    "externalTemp": "ambientTemp",
    "externalVel": "ambientVel",
    "externalConductivity": "ambientConductivity",
    "externalSpecificHeat": "ambientSpecificHeat",
    "externalDensity": "ambientDensity",
    "externalVisc": "ambientVisc",
    "externalEnvironment": "ambientEnvironment",
    # Keep externalHeatCoef, externalGrashof, externalNusselt,
    # externalPrandtl, externalReynolds (HT surface terms)
}

# ==========================================================================
# 1. Update translations.json
# ==========================================================================
trans_path = ROOT / "marlim3" / "translations.json"
data = json.loads(trans_path.read_text(encoding="utf-8"))

new_keys = {}
renamed_count = 0
for en, pt in data["keys"].items():
    if en in RENAMES:
        new_keys[RENAMES[en]] = pt
        renamed_count += 1
    else:
        new_keys[en] = pt

data["keys"] = new_keys
trans_path.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
print(f"translations.json: renamed {renamed_count} EN keys")

# ==========================================================================
# 2. Update English demo JSON files (text-based replacement)
# ==========================================================================
demo_dir = ROOT / "demos"
demo_files = list(demo_dir.glob("*.json"))
# Exclude pt-br subdirectory (Portuguese files use PT keys, unaffected)

total_file_changes = 0
for fpath in demo_files:
    content = fpath.read_text(encoding="utf-8")
    original = content
    for old_key, new_key in RENAMES.items():
        # Replace "oldKey" with "newKey" in JSON context (as keys or references)
        content = content.replace(f'"{old_key}"', f'"{new_key}"')
    if content != original:
        fpath.write_text(content, encoding="utf-8")
        total_file_changes += 1
        print(f"  Updated: {fpath.name}")

print(f"Demo files: {total_file_changes} files updated")

# ==========================================================================
# 3. Summary of new key names for verification
# ==========================================================================
print(f"\nTotal renames applied: {len(RENAMES)}")
print("Done. Next: cmake configure + build + test.")
