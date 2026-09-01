#!/usr/bin/env bash
# Run every calibration in the harness and fail if any instrument is dead.
#
# Why this exists
# ---------------
# A calibration pins fixed text -- a line to corrupt, a marker to find, a
# baseline to compare against. This programme's entire job is to move and rename
# that text. So every stage silently erodes the instruments of the stages before
# it, and nothing warns you: a calibration whose pattern no longer matches does
# not crash, it simply stops testing.
#
# That is not hypothetical. At the close of stage 5 a sweep found:
#   - calibrate-thermal.sh      23 of 24 cases could not inject at all, so the
#                               thermal harness had run since T068 with no proof
#                               it detects error -- while being the only
#                               authority over eleven functions the corpus never
#                               executes;
#   - calibrate-c0ud.sh         7 of 21 dead, recorded as stage 3's next action
#                               (E3-01) and never done, so two further stages ran
#                               on a third-empty instrument;
#   - calibrate-thermal-move.sh two decomposition controls could not locate their
#                               blocks and six mutations no longer matched;
#   - calibrate-solvers.sh      its negative control failed on a renamed function.
#
# Each of those was found by RUNNING the instruments, which is why this script
# exists and why a stage gate should call it -- not only the calibration of the
# stage being closed.
#
# A calibration that reports "not injected", "pattern absent", "not found" or a
# failed control is treated as a failure here, exactly as the individual scripts
# do: the distinction between "did not inject" and "did not detect" is the whole
# reason any of this is trustworthy.
#
# Usage: calibrate-all.sh
# Exit code: 0 only when every calibration passes.
set -uo pipefail
export LC_ALL=C

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "$script_dir/.." && pwd)"
cd "$project_root" || exit 2

red=$'\033[0;31m'; green=$'\033[0;32m'; bold=$'\033[1m'; reset=$'\033[0m'
scratch="$(mktemp -d -t marlim3-calall-XXXXXX)"
trap 'rm -rf "$scratch"' EXIT

# The steady-march decomposition compares against the module as it stood before
# T071 split it. Derived from git so the caller needs no arguments.
PRE_DECOMPOSITION_COMMIT="${MARLIM_PRE_DECOMPOSITION:-949e23b}"
pre_decomposition="$scratch/pre-decomposition.cpp"
git show "$PRE_DECOMPOSITION_COMMIT:src/core/SisProdThermal.cpp" \
    > "$pre_decomposition" 2>/dev/null || {
    printf '%scannot read %s:src/core/SisProdThermal.cpp%s\n' \
           "$red" "$PRE_DECOMPOSITION_COMMIT" "$reset" >&2
    exit 2
}

failures=0
run() {
    local name="$1"; shift
    local log="$scratch/$name.log"
    printf '%s%-26s%s ' "$bold" "$name" "$reset"
    timeout 1800 bash "$script_dir/$name.sh" "$@" > "$log" 2>&1
    local status=$?
    # A dead case never counts as a pass, whatever the exit code says.
    local dead
    # Case-SENSITIVE on the uppercase status words: a summary line reading
    # "0 missed, 0 skipped" is a pass, and matching it case-insensitively turns
    # a clean run into a false failure.
    dead=$(grep -cE 'NOT INJECTED|pattern absent|not found|control failed|MISSED|INVISIBLE TO BOTH' "$log")
    if (( status == 0 && dead == 0 )); then
        printf '%sPASS%s  %s\n' "$green" "$reset" \
            "$(grep -oiE '(CALIBRAT|PASSED)[^\x1b]*' "$log" | tail -1 | cut -c1-58)"
    else
        printf '%sFAIL%s  exit=%s dead=%s\n' "$red" "$reset" "$status" "$dead"
        grep -E 'NOT INJECTED|pattern absent|not found|control failed|MISSED|INVISIBLE TO BOTH|FAILED' \
            "$log" | head -4 | sed 's/^/      /'
        failures=$((failures + 1))
    fi
}

run calibrate-thermal
run calibrate-thermal-move
run calibrate-c0ud
run calibrate-solvers
run calibrate-solver-move
run calibrate-dispatch
run calibrate-performance-gate
run calibrate-steady-decomposition "$pre_decomposition"

printf '\n'
if (( failures > 0 )); then
    printf '%sHARNESS CALIBRATION FAILED -- %s instrument(s)%s\n' \
        "$red" "$failures" "$reset" >&2
    exit 1
fi
printf '%sHARNESS CALIBRATION PASSED -- every instrument proves it detects error%s\n' \
    "$green" "$reset"
