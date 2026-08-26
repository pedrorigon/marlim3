#!/usr/bin/env bash
# Verification layer L1 for the live root finder: run the corpus against an
# instrumented build and compare every SProd::zriddr call to the golden table.
#
# Why this is a separate script from verify-l1.sh
# -----------------------------------------------
# verify-l1.sh compares the synthetic sweep, which needs no model and finishes
# in seconds. This one needs the whole corpus, because zriddr's inputs are
# whatever the pressure searches happen to produce -- there is no way to
# manufacture them. Merging the two would make the fast check slow, and the fast
# check is what makes L1 useful for locating a divergence during a task.
#
# What it compares, and why it is not just the root
# -------------------------------------------------
# Six columns: x1, x2, prod, tipoCC, the returned root, and the number of
# objective evaluations the call made. The last one is the point.
# contracts/busca-raiz.md is explicit that comparing roots proves nothing --
# "dois caminhos de convergência distintos podem chegar ao mesmo valor" -- and
# the golden captured in stage 0 had only the first five columns, so the
# requirement was written down and not implemented (D2-04). The count is taken
# by wrapping SProd::multMarcha rather than anything inside zriddr, which is
# what makes it comparable across the extraction: before it, zriddr calls the
# dispatcher directly; after it, the objective lambda does.
#
# The golden also has to be recaptured to gain that column, and a recapture is
# exactly where a stale artifact slips in (E1-01). So capture mode checks that
# the first five columns of the new table reproduce the stage-0 zriddr.txt
# exactly, which proves the recapture is the same population and not a different
# run wearing its name.
#
# Usage:
#   verify-l1-zriddr.sh capture <instrumented-binary>   recapture the golden
#   verify-l1-zriddr.sh compare <instrumented-binary>   compare against it
#
# Exit code: 0 when every call matches bit for bit.

set -uo pipefail
export LC_ALL=C

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "$script_dir/.." && pwd)"

GOLDEN_DIR="${MARLIM_GOLDEN_DIR:-$project_root/specs/001-refatoracao-sisprod/golden}"
GOLDEN="$GOLDEN_DIR/zriddr-trace.txt"
LEGACY="$GOLDEN_DIR/zriddr.txt"
STRIDE="${MARLIM_GOLDEN_STRIDE:-100}"

mode="${1:?usage: verify-l1-zriddr.sh capture|compare <instrumented-binary>}"
binary="${2:?usage: verify-l1-zriddr.sh capture|compare <instrumented-binary>}"

red=$'\033[0;31m'; green=$'\033[0;32m'; yellow=$'\033[1;33m'; reset=$'\033[0m'

[[ -x "$binary" ]] || { printf '%sbinary missing: %s%s\n' "$red" "$binary" "$reset" >&2; exit 2; }

# The binary has to be newer than the source it claims to represent. verify-l1.sh
# once passed against an instrumented build eight days older than the tree, and
# reported green for a version of the code that no longer existed (E1-01).
source_time=$(stat -c %Y "$project_root/src/core/SisProd.cpp")
binary_time=$(stat -c %Y "$binary")
if (( binary_time < source_time )); then
    printf '%sthe binary is older than src/core/SisProd.cpp -- rebuild it%s\n' \
           "$red" "$reset" >&2
    printf '  binary : %s\n  source : %s\n' \
           "$(date -d "@$binary_time")" "$(date -d "@$source_time")" >&2
    exit 2
fi

work="$(mktemp -d -t marlim3-l1z-XXXXXX)"
trap 'rm -rf "$work"' EXIT
capture_dir="$work/capture"
mkdir -p "$capture_dir"

# The corpus is what the baseline captured, enumerated the same way compare-l2.sh
# enumerates it. Listing demos/*.mr3 instead finds only SEVEN of the fourteen
# models: the Portuguese-language half lives in demos/pt-br/. That is how the
# first run of this script produced exactly half the golden's calls -- and it
# would have looked like a real divergence rather than a script that had not
# looked everywhere.
BASELINE_OUTPUTS="${MARLIM_BASELINE:-$HOME/marlim3-baseline}/saidas"
[[ -d "$BASELINE_OUTPUTS" ]] || {
    printf '%sbaseline outputs missing: %s%s\n' "$red" "$BASELINE_OUTPUTS" "$reset" >&2
    exit 2
}
mapfile -t models < <(cd "$BASELINE_OUTPUTS" && ls -d */ | tr -d '/')
(( ${#models[@]} > 0 )) || { printf '%sno models in %s%s\n' "$red" "$BASELINE_OUTPUTS" "$reset" >&2; exit 2; }

locate_model() {
    [[ -f "$project_root/demos/$1.mr3" ]] && { echo "$project_root/demos/$1.mr3"; return 0; }
    [[ -f "$project_root/demos/pt-br/$1.mr3" ]] && { echo "$project_root/demos/pt-br/$1.mr3"; return 0; }
    return 1
}

printf 'running %d model(s) with stride %s\n' "${#models[@]}" "$STRIDE"
for model in "${models[@]}"; do
    deck="$(locate_model "$model")" || {
        printf '%s  %s: no deck under demos/ or demos/pt-br/%s\n' "$red" "$model" "$reset" >&2
        exit 2
    }
    printf '  %s\n' "$model"
    ( cd "$project_root" && \
      MARLIM_GOLDEN_DIR="$capture_dir" MARLIM_GOLDEN_STRIDE="$STRIDE" \
        "$binary" -s TRANSIENTE -i "$deck" -p demos/ \
                  -d "$work/out-$model" -o "$work/$model.log" > /dev/null 2>&1 )
done

current="$capture_dir/zriddr.txt"
[[ -s "$current" ]] || {
    printf '%sthe run captured no zriddr calls%s\n' "$red" "$reset" >&2
    printf '%sThe binary must be built from a tree instrumented by instrument-golden.py.%s\n' \
           "$yellow" "$reset" >&2
    exit 2
}
sort -o "$current" "$current"   # models may finish in any order

case "$mode" in
capture)
    if [[ -s "$LEGACY" ]]; then
        # The stage-0 table, sorted and reduced to its five columns, must be
        # exactly what the new capture reproduces. If it is not, the recapture
        # is measuring something else and the extra column is worthless.
        if diff -q <(cut -d' ' -f1-5 "$current") <(sort "$LEGACY") > /dev/null 2>&1; then
            printf '%sthe recapture reproduces the stage-0 zriddr.txt exactly (%s calls)%s\n' \
                   "$green" "$(wc -l < "$LEGACY")" "$reset"
        else
            printf '%sthe recapture does NOT reproduce the stage-0 zriddr.txt%s\n' \
                   "$red" "$reset" >&2
            printf '  stage-0 : %s calls\n  now     : %s calls\n' \
                   "$(wc -l < "$LEGACY")" "$(wc -l < "$current")" >&2
            diff <(cut -d' ' -f1-5 "$current") <(sort "$LEGACY") | head -6 >&2
            exit 1
        fi
    else
        printf '%sno stage-0 zriddr.txt to cross-check against%s\n' "$yellow" "$reset" >&2
    fi
    cp "$current" "$GOLDEN"
    printf '%scaptured %s call(s) to %s%s\n' \
           "$green" "$(wc -l < "$GOLDEN")" "$GOLDEN" "$reset"
    exit 0
    ;;
compare)
    [[ -s "$GOLDEN" ]] || {
        printf '%sno golden at %s -- run capture against a pristine build first%s\n' \
               "$red" "$GOLDEN" "$reset" >&2
        exit 2
    }
    golden_lines=$(wc -l < "$GOLDEN")
    current_lines=$(wc -l < "$current")
    if [[ "$golden_lines" != "$current_lines" ]]; then
        printf '%sL1 FAILED -- call count differs: golden %s, current %s%s\n' \
               "$red" "$golden_lines" "$current_lines" "$reset" >&2
        exit 1
    fi
    if diff -q "$GOLDEN" "$current" > /dev/null 2>&1; then
        printf '%sL1 PASSED -- %s zriddr call(s) identical, roots AND evaluation counts%s\n' \
               "$green" "$golden_lines" "$reset"
        awk '{ e += $6 } END { printf "  %d objective evaluations across the corpus\n", e }' "$GOLDEN"
        exit 0
    fi
    printf '%sL1 FAILED -- divergent zriddr calls%s\n' "$red" "$reset" >&2
    # Say which column moved: a root that changed is a different answer, a count
    # that changed with the same root is a different path to the same answer.
    diff <(cut -d' ' -f1-5 "$GOLDEN") <(cut -d' ' -f1-5 "$current") > /dev/null 2>&1 \
        && printf '%s  roots identical -- only the evaluation counts differ.%s\n' "$yellow" "$reset" >&2 \
        || printf '%s  the returned roots differ.%s\n' "$red" "$reset" >&2
    diff "$GOLDEN" "$current" | head -10 >&2
    printf '\n  %s divergent call(s) of %s\n' \
           "$(diff "$GOLDEN" "$current" | grep -c '^<')" "$golden_lines" >&2
    exit 1
    ;;
*)
    printf 'unknown mode %s -- use capture or compare\n' "$mode" >&2
    exit 2
    ;;
esac
