#!/usr/bin/env bash
# Verification gate 4: wall-clock performance against the stored baseline.
#
# Design notes, derived from measurements taken in stage 0 (T008)
# -------------------------------------------------------------
# Naively comparing "median of 3" against a 3% threshold does not work on this
# corpus. Measured run-to-run spread of the SAME binary on an idle machine:
#
#     parada-longo-Combinado-BCS-GLC-PIG-completo   19.5%
#     BCS-longo-eficMotor                           17.1%
#     MultiBCS                                      10.1%
#     2zones-2GLVs-2-Check-correcThermProf           8.1%
#
# The noise floor is up to 6x the 3% threshold of SC-005/FR-024. A gate built on
# the median would fail correct refactorings and pass real regressions.
#
# Two decisions follow:
#
#   1. Use the MINIMUM, not the median. Timing noise is one-sided: interference
#      only ever adds time, never removes it. The fastest of N runs is the least
#      contaminated estimate of the true cost, and it is far more stable than
#      the median across repetitions.
#
#   2. Skip models below MIN_MEASURABLE_SECONDS. Two models complete in under
#      one second, below timer resolution, and a percentage against a ~0 s
#      baseline is undefined.
#
# The 3% threshold itself comes from the specification and is NOT changed here.
# What changes is the estimator, so that the comparison is about the code rather
# than about scheduler luck.
#
# Usage:
#   performance-gate.sh                # every model in the baseline timings
#   performance-gate.sh <model> ...    # only the named models
#
# Environment:
#   MARLIM_BIN            binary under test (default build/Marlim3)
#   MARLIM_BASELINE       baseline root (default ~/marlim3-baseline)
#   MARLIM_PERF_REPS      repetitions per model (default 3)
#   MARLIM_PERF_THRESHOLD regression threshold in percent (default 3)
#
# Exit code: 0 when no model regresses beyond the threshold; 1 otherwise.

set -uo pipefail

# Force the C locale for all numeric formatting and parsing.
#
# Without this, a pt_BR system produces comma decimal separators in awk/printf
# output ("12,0"), which then fails to parse as a number. The failure mode is
# silent and dangerous: a 12% regression was formatted as "12,0", compared as
# the integer 12 against the threshold in a way that did not trip, and the gate
# reported PASSED. A performance gate that cannot compare numbers is worse than
# no gate, because it manufactures confidence.
export LC_ALL=C

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "$script_dir/.." && pwd)"

BASELINE_DIR="${MARLIM_BASELINE:-$HOME/marlim3-baseline}"
BASELINE_TIMINGS="$BASELINE_DIR/tempos-baseline.txt"
CURRENT_BINARY="${MARLIM_BIN:-$project_root/build/Marlim3}"
REPETITIONS="${MARLIM_PERF_REPS:-3}"
THRESHOLD="${MARLIM_PERF_THRESHOLD:-3}"
MIN_MEASURABLE_SECONDS=5

red=$'\033[0;31m'; green=$'\033[0;32m'; yellow=$'\033[1;33m'; reset=$'\033[0m'

[[ -x "$CURRENT_BINARY" ]] || { printf '%scurrent binary missing: %s%s\n' "$red" "$CURRENT_BINARY" "$reset" >&2; exit 2; }
[[ -f "$BASELINE_TIMINGS" ]] || { printf '%sbaseline timings missing: %s (run T008)%s\n' "$red" "$BASELINE_TIMINGS" "$reset" >&2; exit 2; }

cd "$project_root" || exit 2

work_dir="$(mktemp -d -t marlim3-perf-XXXXXX)"
trap 'rm -rf "$work_dir"' EXIT

locate_model() {
    local name="$1"
    [[ -f "demos/$name.mr3" ]] && { echo "demos/$name.mr3"; return 0; }
    [[ -f "demos/pt-br/$name.mr3" ]] && { echo "demos/pt-br/$name.mr3"; return 0; }
    return 1
}

baseline_seconds() {
    awk -v want="$1" '$2 == want { gsub(",", ".", $1); print $1; exit }' "$BASELINE_TIMINGS"
}

if (( $# > 0 )); then
    models=("$@")
else
    mapfile -t models < <(awk '{print $2}' "$BASELINE_TIMINGS")
fi

printf '%-52s %10s %10s %8s\n' "MODEL" "BASELINE" "CURRENT" "DELTA"
printf '%s\n' "--------------------------------------------------------------------------------"

regressions=0
skipped=0
measured=0

for model in "${models[@]}"; do
    reference="$(baseline_seconds "$model")"
    [[ -n "$reference" ]] || { printf '%-52s %10s\n' "$model" "no baseline"; continue; }

    if awk -v r="$reference" -v m="$MIN_MEASURABLE_SECONDS" 'BEGIN { exit !(r < m) }'; then
        printf '%-52s %10.2f %10s %8s   %sskipped: below %ss, timer resolution%s\n' \
               "$model" "$reference" "-" "-" "$yellow" "$MIN_MEASURABLE_SECONDS" "$reset"
        skipped=$((skipped + 1))
        continue
    fi

    input_path="$(locate_model "$model")" || {
        printf '%-52s %10s\n' "$model" "input missing"
        regressions=$((regressions + 1)); continue
    }

    fastest=""
    for (( rep = 1; rep <= REPETITIONS; rep++ )); do
        run_dir="$work_dir/$model-$rep"; mkdir -p "$run_dir"
        started=$(date +%s.%N)
        "$CURRENT_BINARY" -s TRANSIENTE -i "$input_path" -p demos/ \
                          -d "$run_dir" -o "$run_dir/$model.log" \
                          > "$run_dir/_stdout.txt" 2>&1
        finished=$(date +%s.%N)
        elapsed=$(awk -v a="$started" -v b="$finished" 'BEGIN { printf "%.3f", b - a }')
        rm -rf "$run_dir"
        if [[ -z "$fastest" ]] || awk -v e="$elapsed" -v f="$fastest" 'BEGIN { exit !(e < f) }'; then
            fastest="$elapsed"
        fi
    done

    delta=$(awk -v r="$reference" -v c="$fastest" 'BEGIN { printf "%.1f", (c - r) / r * 100 }')
    measured=$((measured + 1))

    if awk -v d="$delta" -v t="$THRESHOLD" 'BEGIN { exit !(d > t) }'; then
        printf '%-52s %10.2f %10.2f %7s%%   %sREGRESSION%s\n' \
               "$model" "$reference" "$fastest" "$delta" "$red" "$reset"
        regressions=$((regressions + 1))
    else
        printf '%-52s %10.2f %10.2f %7s%%\n' "$model" "$reference" "$fastest" "$delta"
    fi
done

echo
printf 'estimator            : minimum of %d runs (noise is one-sided)\n' "$REPETITIONS"
printf 'threshold            : %s%%\n' "$THRESHOLD"
printf 'models measured      : %d\n' "$measured"
printf 'models skipped       : %d (baseline below %ss)\n' "$skipped" "$MIN_MEASURABLE_SECONDS"

if (( regressions > 0 )); then
    printf '%sPERFORMANCE GATE FAILED -- %d model(s) beyond %s%%%s\n' "$red" "$regressions" "$THRESHOLD" "$reset" >&2
    printf '%sRe-run on an idle machine before concluding: measured noise reaches 19.5%%.%s\n' "$yellow" "$reset" >&2
    exit 1
fi

printf '%sPERFORMANCE GATE PASSED -- no model beyond %s%%%s\n' "$green" "$THRESHOLD" "$reset"
exit 0
