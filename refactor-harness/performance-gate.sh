#!/usr/bin/env bash
# Gate 4 -- has the refactoring made the engine slower?
#
# HOW THE COMPARISON IS MADE, AND WHY
#
# Both sides are measured now, alternately, in one session: the baseline commit
# is built from source and its binary is run head-to-head against the working
# tree's binary, model by model. Nothing is read from a file written on an
# earlier day.
#
# The two designs this replaces both failed, for reasons worth keeping:
#
#   1. Stored timings. A baseline captured once freezes one draw of a noisy
#      estimator. When that draw is fast, every later comparison inherits the
#      error and re-running cannot fix it -- the defective side of the comparison
#      is the side already on disk. Measured: a stored minimum of 125.278 s that
#      three later sessions put at 131.32, 131.93 and 135.40, producing "5.3%
#      regression" verdicts on source nobody had touched.
#
#   2. A preserved baseline BINARY. Better, but it ages against the environment.
#      The binary kept from 10 Aug stopped reproducing its own capture after a
#      system library update pulled in libmvec: same source, same compiler, same
#      flags, 99 of 102 output files different. Comparing a binary built then
#      against one built now measures the toolchain as much as the code.
#
# Building the baseline commit at comparison time removes both. Both binaries
# come from the same compiler and the same libraries, and both meet the same
# machine load at the same moment, so what differs between them is the source.
#
# The build is cached: it is redone only when the baseline commit changes or the
# cached tree is missing.
#
# THE MACHINE IS NOT IDLE, AND CANNOT BE MADE IDLE
#
# Earlier versions instructed the operator to "re-run on an idle machine". That
# advice was unusable: this is a desktop in use, with a browser, a shell and
# long-running user processes taking ~60% of the CPU and varying on their own.
# Pairing is what makes measurement possible here -- it does not need a quiet
# machine, only that both binaries meet the same machine.
#
# Usage:
#   performance-gate.sh                 full corpus
#   performance-gate.sh <model> ...     only the named models
#
# Environment:
#   MARLIM_BIN            binary under test (default build/Marlim3)
#   MARLIM_PERF_REPS      minimum alternating repetitions per side (default 4;
#                         odd values are promoted to the next even value)
#   MARLIM_PERF_THRESHOLD regression threshold in percent (default 3)
#   MARLIM_PERF_TREE      cached baseline worktree (default ~/marlim3-perf-baseline)
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
BASELINE_OUTPUTS="$BASELINE_DIR/saidas"
CURRENT_BINARY="${MARLIM_BIN:-$project_root/build/Marlim3}"
PERF_TREE="${MARLIM_PERF_TREE:-$HOME/marlim3-perf-baseline}"
REQUESTED_REPETITIONS="${MARLIM_PERF_REPS:-4}"
THRESHOLD="${MARLIM_PERF_THRESHOLD:-3}"
MIN_MEASURABLE_SECONDS=5

red=$'\033[0;31m'; green=$'\033[0;32m'; yellow=$'\033[1;33m'; reset=$'\033[0m'

[[ "$REQUESTED_REPETITIONS" =~ ^[1-9][0-9]*$ ]] || {
    printf '%sMARLIM_PERF_REPS must be a positive integer, got %s%s\n' \
           "$red" "$REQUESTED_REPETITIONS" "$reset" >&2
    exit 2
}
REPETITIONS="$REQUESTED_REPETITIONS"
if (( REPETITIONS % 2 != 0 )); then
    REPETITIONS=$((REPETITIONS + 1))
    printf '%spromoting MARLIM_PERF_REPS=%d to %d so run order is balanced%s\n' \
           "$yellow" "$REQUESTED_REPETITIONS" "$REPETITIONS" "$reset"
fi

[[ -x "$CURRENT_BINARY" ]] || {
    printf '%scurrent binary missing: %s%s\n' "$red" "$CURRENT_BINARY" "$reset" >&2; exit 2; }

BASELINE_COMMIT="$(cat "$BASELINE_DIR/commit" 2>/dev/null)"
[[ -n "$BASELINE_COMMIT" ]] || {
    printf '%sno baseline commit recorded at %s/commit%s\n' "$red" "$BASELINE_DIR" "$reset" >&2; exit 2; }

cd "$project_root" || exit 2

work_dir="$(mktemp -d -t marlim3-perf-XXXXXX)"
trap 'rm -rf "$work_dir"' EXIT

# ------------------------------------------------ baseline build, cached ----
# Rebuilding the baseline for every gate run would cost minutes on a comparison
# that does not change between runs, so the tree is kept and reused. It is
# rebuilt when the recorded baseline commit moves, which is the only thing that
# invalidates it.
baseline_binary="$PERF_TREE/tree/build/Marlim3"
stamp="$PERF_TREE/built-from-commit"

ensure_baseline_build() {
    if [[ -x "$baseline_binary" && "$(cat "$stamp" 2>/dev/null)" == "$BASELINE_COMMIT" ]]; then
        printf 'baseline binary: cached build of %s\n' "${BASELINE_COMMIT:0:9}"
        return 0
    fi

    printf '%sbuilding baseline commit %s (cached for later runs)%s\n' \
           "$yellow" "${BASELINE_COMMIT:0:9}" "$reset"

    git -C "$project_root" worktree remove --force "$PERF_TREE/tree" > /dev/null 2>&1
    rm -rf "$PERF_TREE"
    git -C "$project_root" worktree prune
    mkdir -p "$PERF_TREE"

    git -C "$project_root" worktree add --detach "$PERF_TREE/tree" "$BASELINE_COMMIT" > /dev/null 2>&1 || {
        printf '%scould not create a worktree at %s%s\n' "$red" "$BASELINE_COMMIT" "$reset" >&2
        return 1
    }

    (
        cd "$PERF_TREE/tree" || exit 1
        cmake --preset gcc-release > "$PERF_TREE/build.log" 2>&1 &&
        cmake --build --preset gcc-release >> "$PERF_TREE/build.log" 2>&1
    ) || {
        printf '%sbaseline build failed -- see %s/build.log%s\n' "$red" "$PERF_TREE" "$reset" >&2
        return 1
    }

    [[ -x "$baseline_binary" ]] || {
        printf '%sbaseline build produced no binary at %s%s\n' "$red" "$baseline_binary" "$reset" >&2
        return 1
    }

    printf '%s' "$BASELINE_COMMIT" > "$stamp"
    printf '%sbaseline built%s\n' "$green" "$reset"
}

ensure_baseline_build || exit 2

locate_model() {
    local name="$1"
    [[ -f "demos/$name.mr3" ]] && { echo "demos/$name.mr3"; return 0; }
    [[ -f "demos/pt-br/$name.mr3" ]] && { echo "demos/pt-br/$name.mr3"; return 0; }
    return 1
}

# Time one run of one binary. Prints elapsed seconds.
#
# Both binaries are invoked from the project root against the same input files,
# so the inputs are identical by construction and only the executable differs.
time_once() {
    local binary="$1" input_path="$2" run_dir="$3"
    local started finished exit_code
    mkdir -p "$run_dir"
    started=$(date +%s.%N)
    "$binary" -s TRANSIENTE -i "$input_path" -p demos/ \
              -d "$run_dir" -o "$run_dir/run.log" \
              > "$run_dir/_stdout.txt" 2>&1
    exit_code=$?
    finished=$(date +%s.%N)
    if (( exit_code != 0 )); then
        printf '%srun failed with exit %d: %s %s; output retained at %s%s\n' \
               "$red" "$exit_code" "$binary" "$input_path" "$run_dir" "$reset" >&2
        return 1
    fi
    rm -rf "$run_dir"
    awk -v a="$started" -v b="$finished" 'BEGIN { printf "%.3f", b - a }'
}

time_or_abort() {
    local output_name="$1" binary="$2" input_path="$3" run_dir="$4" measured
    if ! measured="$(time_once "$binary" "$input_path" "$run_dir")"; then
        printf '%sPERFORMANCE GATE ABORTED -- failed simulations are not timings%s\n' \
               "$red" "$reset" >&2
        trap - EXIT
        exit 2
    fi
    printf -v "$output_name" '%s' "$measured"
}

keep_fastest() {
    local candidate="$1" incumbent="$2"
    if [[ -z "$incumbent" ]] || awk -v e="$candidate" -v f="$incumbent" 'BEGIN { exit !(e < f) }'; then
        printf '%s' "$candidate"
    else
        printf '%s' "$incumbent"
    fi
}

if (( $# > 0 )); then
    models=("$@")
else
    mapfile -t models < <(cd "$BASELINE_OUTPUTS" && ls -d */ 2>/dev/null | tr -d '/')
fi

(( ${#models[@]} > 0 )) || {
    printf '%sno models to measure%s\n' "$red" "$reset" >&2; exit 2; }

printf '\n%-52s %10s %10s %8s\n' "MODEL" "BASELINE" "CURRENT" "DELTA"
printf '%s\n' "--------------------------------------------------------------------------------"

regressions=0
skipped=0
measured=0

for model in "${models[@]}"; do
    input_path="$(locate_model "$model")" || {
        printf '%-52s %10s\n' "$model" "input missing"
        regressions=$((regressions + 1)); continue
    }

    # Alternate the two binaries rather than measuring all of one and then all
    # of the other. Whatever the machine does during the run -- a browser waking
    # up, a background job starting -- then lands on both sides instead of on
    # whichever side happened to run while it happened.
    #
    # The order is swapped every repetition (baseline first, then current first)
    # so each binary spends the same number of runs in each position. Running one
    # of them always second would hand it whatever the position is worth: warm
    # caches, a settled clock, a finished background task. That advantage would
    # be small but one-directional, and this gate only fails on positive deltas,
    # so a bias favouring the working tree is a bias toward passing -- it would
    # hide exactly the regressions the gate exists to catch. Counterbalancing
    # removes the whole class by construction instead of arguing it is small.
    best_base=""; best_curr=""
    for (( rep = 1; rep <= REPETITIONS; rep++ )); do
        if (( rep % 2 == 1 )); then
            time_or_abort elapsed "$baseline_binary" "$input_path" "$work_dir/b-$model-$rep"
            best_base="$(keep_fastest "$elapsed" "$best_base")"
            time_or_abort elapsed "$CURRENT_BINARY" "$input_path" "$work_dir/c-$model-$rep"
            best_curr="$(keep_fastest "$elapsed" "$best_curr")"
        else
            time_or_abort elapsed "$CURRENT_BINARY" "$input_path" "$work_dir/c-$model-$rep"
            best_curr="$(keep_fastest "$elapsed" "$best_curr")"
            time_or_abort elapsed "$baseline_binary" "$input_path" "$work_dir/b-$model-$rep"
            best_base="$(keep_fastest "$elapsed" "$best_base")"
        fi
    done

    # The floor is judged on the freshly measured baseline, not on a stored
    # number, so it tracks the machine the gate is actually running on. Below
    # it, 3% is smaller than process start-up variance and carries no
    # information about the code.
    if awk -v r="$best_base" -v m="$MIN_MEASURABLE_SECONDS" 'BEGIN { exit !(r < m) }'; then
        printf '%-52s %10.2f %10.2f %8s   %sskipped: below %ss, timer resolution%s\n' \
               "$model" "$best_base" "$best_curr" "-" "$yellow" "$MIN_MEASURABLE_SECONDS" "$reset"
        skipped=$((skipped + 1))
        continue
    fi

    delta=$(awk -v r="$best_base" -v c="$best_curr" 'BEGIN { printf "%.1f", (c - r) / r * 100 }')
    measured=$((measured + 1))

    if awk -v d="$delta" -v t="$THRESHOLD" 'BEGIN { exit !(d > t) }'; then
        printf '%-52s %10.2f %10.2f %7s%%   %sREGRESSION%s\n' \
               "$model" "$best_base" "$best_curr" "$delta" "$red" "$reset"
        regressions=$((regressions + 1))
    else
        printf '%-52s %10.2f %10.2f %7s%%\n' "$model" "$best_base" "$best_curr" "$delta"
    fi
done

echo
printf 'comparison           : baseline commit %s built now, alternated with the working tree\n' "${BASELINE_COMMIT:0:9}"
printf 'estimator            : minimum of %d balanced alternating runs per side (requested %d; noise is one-sided)\n' \
    "$REPETITIONS" "$REQUESTED_REPETITIONS"
printf 'threshold            : %s%%\n' "$THRESHOLD"
printf 'models measured      : %d\n' "$measured"
printf 'models skipped       : %d (baseline below %ss)\n' "$skipped" "$MIN_MEASURABLE_SECONDS"

if (( regressions > 0 )); then
    printf '%sPERFORMANCE GATE FAILED -- %d model(s) beyond %s%%%s\n' \
           "$red" "$regressions" "$THRESHOLD" "$reset" >&2
    printf '%sBoth binaries were built from the same toolchain and measured in the same\n' "$yellow"
    printf 'session, so this is not explained by a stale baseline or by machine load.\n'
    printf 'Re-run once to rule out a transient, then treat it as a real regression.%s\n' "$reset" >&2
    exit 1
fi

if (( measured == 0 )); then
    printf '%sPERFORMANCE GATE NOT EVALUATED -- no model reached the %ss measurement floor%s\n' \
           "$yellow" "$MIN_MEASURABLE_SECONDS" "$reset" >&2
    exit 2
fi

printf '%sPERFORMANCE GATE PASSED -- no model beyond %s%%%s\n' "$green" "$THRESHOLD" "$reset"
exit 0
