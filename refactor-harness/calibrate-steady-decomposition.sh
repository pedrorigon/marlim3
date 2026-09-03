#!/usr/bin/env bash
# Calibrates check-steady-decomposition.py: a checker that never fails is not a
# checker. Each corruption below is injected into a specific line range of one
# decomposed body, and the run CONFIRMS the injection landed before asking for a
# verdict -- a corruption that fails to apply is reported as SKIP, never counted
# as a success.
#
# Usage: calibrate-steady-decomposition.sh <pre-decomposition-SisProdThermal.cpp>
set -uo pipefail
export LC_ALL=C

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "$script_dir/.." && pwd)"
target="$project_root/src/core/SisProdThermal.cpp"
baseline="${1:-}"
[[ -f "$baseline" ]] || { echo "usage: $0 <pre-decomposition-SisProdThermal.cpp>" >&2; exit 2; }

red=$'\033[0;31m'; green=$'\033[0;32m'; yellow=$'\033[1;33m'; reset=$'\033[0m'
scratch="$(mktemp -d -t marlim3-caldecomp-XXXXXX)"
trap 'cp "$scratch/pristine.cpp" "$target"; rm -rf "$scratch"' EXIT
cp "$target" "$scratch/pristine.cpp"

# Ranges are resolved from the file itself, so the calibration does not rot when
# the module shifts.
range() {
    python3 - "$1" <<'PY'
import re, sys
lines = open('src/core/SisProdThermal.cpp', encoding='utf-8').read().split('\n')
name = sys.argv[1]
i = next(k for k, l in enumerate(lines)
         if re.match(rf'^[\w:<>*&]+[\w:<>*&\s]*\b{name}\(', l))
j = i
while '{' not in lines[j]: j += 1
d = 0; k = j
while True:
    d += lines[k].count('{') - lines[k].count('}')
    if d == 0: break
    k += 1
print(f"{i+1},{k+1}")
PY
}

detected=0; missed=0; skipped=0
probe() { # <function> <sed-body> <label>
    local span; span="$(range "$1")"
    cp "$scratch/pristine.cpp" "$target"
    sed -i "${span}{$2}" "$target"
    if diff -q "$scratch/pristine.cpp" "$target" > /dev/null; then
        printf '%s  SKIP      %s (corruption did not apply)%s\n' "$yellow" "$3" "$reset"
        skipped=$((skipped + 1)); return
    fi
    if python3 "$script_dir/check-steady-decomposition.py" "$baseline" > /dev/null 2>&1; then
        printf '%s  MISSED    %s%s\n' "$red" "$3" "$reset"; missed=$((missed + 1))
    else
        printf '%s  detected  %s%s\n' "$green" "$3" "$reset"; detected=$((detected + 1))
    fi
}

probe computeSteadySourceTerms   's/sourceSpecificHeatRatio = 1\./sourceSpecificHeatRatio = 1.0001/' 'constant in the direct source terms'
probe computeSteadySourceTerms   's/cellLength/meanCellLength/'            'data source swapped, direct sources'
probe applySteadyAnnulusCoupling 's/100\./100.0001/'                       'literal in the direct annulus coupling'
probe computeSteadyKineticTerm   's/1e-3/1e-4/'                            'tolerance in the direct kinetic term'
# The literal this used to corrupt is now the named constant kPhaseChangeFloor,
# so the probe targets the comparison itself instead. A pattern that no longer
# matches is reported as SKIP and never counted as a pass -- which is how this
# case announced that it had stopped testing.
probe computeSteadyLatentHeatTerm 's/> kPhaseChangeFloor/>= kPhaseChangeFloor/' 'guard in the direct latent term'
probe computeReverseSteadySourceTerms 's/cellIndex + 1/cellIndex - 1/'     'march direction flipped, reverse sources'
probe computeReverseSteadyLatentHeatTerm 's/latentHeatTerm = 0\./latentHeatTerm = 0.1/' 'initial value, reverse latent term'
probe applyReverseSteadyAnnulusCoupling 's/annulusResistance = 0\./annulusResistance = 1./' 'reverse annulus resistance'
probe advanceSteadyTemperature   's/temperatureSpatialCoefficient/pressureSpatialCoefficient/' 'coefficient swapped in the coordinator'

# The checker normalises away breakpoint anchors -- guards whose entire body is
# an int declaration and an assignment of zero -- because fourteen of them were
# deleted from the module and the baseline still carries them. A normaliser that
# is one character too generous stops being a normaliser and starts being a
# blind spot, so this case injects a block that LOOKS like an anchor and carries
# one live statement. It must still be detected.
probe computeSteadyKineticTerm \
    's/double kineticTerm/if (cellIndex == 3) { int debugStop; debugStop = 0; kineticTermScale = 0.; }\n    double kineticTerm/' \
    'live statement hidden inside a breakpoint-anchor shape'

# expand_aliases restores what a `Cel &rightCell = ...` declaration stands for
# before comparing. Thirty-eight functions declare one, so an expansion that is
# not scoped to the declaring function lets any of them speak for all -- and a
# declaration pointed at the wrong neighbour is repaired by a correct one
# elsewhere before the comparison sees it. That is not hypothetical: the
# file-wide first version of expand_aliases missed exactly this, and this case
# is why. The corruption is one character, in a line that is not arithmetic, and
# it redirects forty-six reads.
probe computeReverseSteadySourceTerms \
    's/Cel &rightCell = state.cells\[cellIndex + 1\];/Cel \&rightCell = state.cells[cellIndex - 1];/' \
    'alias declaration pointed at the wrong neighbour'

cp "$scratch/pristine.cpp" "$target"
printf '\n%s cases: %s detected, %s missed, %s skipped\n' \
    "$((detected + missed + skipped))" "$detected" "$missed" "$skipped"
if (( missed > 0 )); then
    printf '%sDECOMPOSITION CALIBRATION FAILED%s\n' "$red" "$reset" >&2; exit 1
fi
printf '%sDECOMPOSITION CALIBRATION PASSED -- %s cases%s\n' "$green" "$detected" "$reset"
