#!/usr/bin/env bash
# Calibrate the guards around performance-gate.sh without running the slow suite.
# This does not validate timing accuracy; it proves the gate refuses three ways
# of producing a verdict without a valid measurement.
set -uo pipefail
export LC_ALL=C

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "$script_dir/.." && pwd)"
gate="$script_dir/performance-gate.sh"
current_binary="$project_root/build/Marlim3"

[[ -x "$current_binary" ]] || {
    echo "current binary missing: $current_binary" >&2
    exit 2
}

scratch="$(mktemp -d -t marlim3-perf-calibration-XXXXXX)"
trap 'rm -rf "$scratch"' EXIT

failures=0

expect_exit_and_text() {
    local case_name="$1" expected_exit="$2" expected_text="$3"
    shift 3
    local output exit_code
    set +e
    output="$("$@" 2>&1)"
    exit_code=$?
    set -e
    if (( exit_code != expected_exit )) || ! grep -qF "$expected_text" <<< "$output"; then
        printf 'FAIL  %-24s exit=%d, expected=%d, text=%q\n' \
               "$case_name" "$exit_code" "$expected_exit" "$expected_text" >&2
        printf '%s\n' "$output" >&2
        failures=$((failures + 1))
        return
    fi
    printf 'OK    %s\n' "$case_name"
}

bash -n "$gate" || exit 2

expect_exit_and_text invalid-repetitions 2 \
    'MARLIM_PERF_REPS must be a positive integer' \
    env TMPDIR="$scratch" MARLIM_PERF_REPS=abc bash "$gate" simplifiedProduction

expect_exit_and_text failed-simulation 2 \
    'failed simulations are not timings' \
    env TMPDIR="$scratch" MARLIM_BIN=/bin/false MARLIM_PERF_REPS=2 \
        bash "$gate" simplifiedProduction

expect_exit_and_text no-measurable-model 2 \
    'PERFORMANCE GATE NOT EVALUATED' \
    env TMPDIR="$scratch" MARLIM_PERF_REPS=1 bash "$gate" simplifiedProduction

expect_exit_and_text odd-count-balanced 2 \
    'minimum of 2 balanced alternating runs per side (requested 1' \
    env TMPDIR="$scratch" MARLIM_PERF_REPS=1 bash "$gate" simplifiedProduction

if (( failures > 0 )); then
    printf 'PERFORMANCE CALIBRATION FAILED -- %d case(s)\n' "$failures" >&2
    exit 1
fi

echo 'PERFORMANCE CALIBRATION PASSED -- 4 cases'