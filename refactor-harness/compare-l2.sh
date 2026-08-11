#!/usr/bin/env bash
# Verification layer L2: compare the CURRENT binary's output against the stored
# baseline, model by model, file by file.
#
# Criterion: byte-for-byte equality after normalize-output.sh, which strips only
# timestamps, wall-clock duration, absolute paths and the 32 decorative phrases
# enumerated from src/include/SisProd.h. There is no numeric tolerance: any
# surviving difference is a real divergence (Principle I).
#
# The comparison always runs against the stored baseline, never against the
# previous stage -- comparing against the previous stage would let a small drift
# per stage accumulate across ten stages without ever tripping a gate (FR-007).
#
# Usage:
#   compare-l2.sh                 # every model in the baseline
#   compare-l2.sh <model> ...     # only the named models (fast cycle)
#
# Exit code: 0 when everything is equivalent; 1 on the first real divergence,
# reporting model, file and the first differing line.

set -uo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "$script_dir/.." && pwd)"
normalize="$script_dir/normalize-output.sh"

BASELINE_DIR="${MARLIM_BASELINE:-$HOME/marlim3-baseline}"
BASELINE_OUTPUTS="$BASELINE_DIR/saidas"
CURRENT_BINARY="${MARLIM_BIN:-$project_root/build/Marlim3}"
WORK_DIR="${MARLIM_L2_DIR:-$(mktemp -d -t marlim3-l2-XXXXXX)}"

# Files kept out of the comparison, with the reason:
#   time.txt   harness instrumentation, not engine output
#
# _stdout.txt IS compared. It is mostly decoration -- the drawn phrase and a
# banner -- but it also carries "ARQUIVO DE LOG:" and any error the engine
# prints only to the console. Excluding it would leave a blind spot exactly
# where a behaviour change is most likely to surface first, and the normalizer
# already handles everything in it that legitimately varies between runs.
EXCLUDED='^(time\.txt)$'

red=$'\033[0;31m'; green=$'\033[0;32m'; yellow=$'\033[1;33m'; reset=$'\033[0m'

fail()   { printf '%s%s%s\n' "$red" "$*" "$reset" >&2; }
pass()   { printf '%s%s%s\n' "$green" "$*" "$reset"; }
notice() { printf '%s%s%s\n' "$yellow" "$*" "$reset"; }

[[ -x "$CURRENT_BINARY" ]] || { fail "current binary missing: $CURRENT_BINARY"; exit 2; }
[[ -d "$BASELINE_OUTPUTS" ]] || { fail "baseline missing: $BASELINE_OUTPUTS (run T007)"; exit 2; }

cd "$project_root" || exit 2

if (( $# > 0 )); then
    models=("$@")
else
    mapfile -t models < <(cd "$BASELINE_OUTPUTS" && ls -d */ | tr -d '/')
fi

locate_model() {
    local name="$1"
    [[ -f "demos/$name.mr3" ]] && { echo "demos/$name.mr3"; return 0; }
    [[ -f "demos/pt-br/$name.mr3" ]] && { echo "demos/pt-br/$name.mr3"; return 0; }
    return 1
}

diverged_models=0
equivalent_models=0

for model in "${models[@]}"; do
    input_path="$(locate_model "$model")" || {
        fail "[$model] model not found under demos/ or demos/pt-br/"
        diverged_models=$((diverged_models + 1)); continue
    }

    current_run="$WORK_DIR/$model"
    rm -rf "$current_run"; mkdir -p "$current_run"

    "$CURRENT_BINARY" -s TRANSIENTE -i "$input_path" -p demos/ \
                      -d "$current_run" -o "$current_run/$model.log" \
                      > "$current_run/_stdout.txt" 2>&1
    exit_code=$?
    if (( exit_code != 0 )); then
        fail "[$model] run exited with $exit_code -- see $current_run/_stdout.txt"
        diverged_models=$((diverged_models + 1)); continue
    fi

    reference_run="$BASELINE_OUTPUTS/$model"
    model_diverged=0

    # Union of both sides, so a missing or an extra file is caught too, not just
    # differing content.
    mapfile -t artifacts < <(
        { ls "$reference_run"; ls "$current_run"; } 2>/dev/null \
        | grep -vE "$EXCLUDED" | sort -u
    )

    for artifact in "${artifacts[@]}"; do
        if [[ ! -f "$reference_run/$artifact" ]]; then
            fail "[$model] NEW file, absent from baseline: $artifact"
            model_diverged=1; continue
        fi
        if [[ ! -f "$current_run/$artifact" ]]; then
            fail "[$model] MISSING file in current run: $artifact"
            model_diverged=1; continue
        fi

        expected="$WORK_DIR/.expected"; actual="$WORK_DIR/.actual"
        bash "$normalize" "$reference_run/$artifact" "$reference_run" > "$expected" 2>/dev/null
        bash "$normalize" "$current_run/$artifact"   "$current_run"   > "$actual"   2>/dev/null

        if ! diff -q "$expected" "$actual" > /dev/null 2>&1; then
            fail "[$model] DIVERGENCE in $artifact"
            diff "$expected" "$actual" | head -6 | sed 's/^/      /' >&2
            model_diverged=1
        fi
    done

    if (( model_diverged )); then
        diverged_models=$((diverged_models + 1))
    else
        pass "[$model] equivalent (${#artifacts[@]} files)"
        equivalent_models=$((equivalent_models + 1))
    fi
done

echo
echo "----------------------------------------------------------------"
printf 'equivalent models : %d/%d\n' "$equivalent_models" "${#models[@]}"
printf 'work directory    : %s\n' "$WORK_DIR"

if (( diverged_models > 0 )); then
    fail "L2 FAILED -- $diverged_models model(s) diverged"
    notice "Do not proceed to the next task. Locate, classify, revert if needed."
    exit 1
fi

pass "L2 PASSED -- bit-for-bit equivalence against the baseline on every model"
exit 0
