#!/usr/bin/env bash
# Verification layer L1 for the drift-flux regime aggregators.
#
# WHY THIS EXISTS SEPARATELY FROM verify-l1.sh
#
# verify-l1.sh sweeps the five correlations. It cannot reach the three
# aggregators, because those take a correlation selector that the simulator
# reads from configuration, and Leitura.cpp pins arq.CorreDisper to 1. So
# selectors 0, 4 and 5 -- including the angle blend, the most delicate
# arithmetic in the module -- never execute anywhere in the demo corpus. L2 and
# L3 compare the output of code that did not run and stay green whatever
# happens to it.
#
# Since the aggregators are pure, they can be called directly: no model file and
# no configuration are read, and the values below are the sweep's own. It covers
# every selector each regime accepts, plus values outside every accepted set,
# where c0 and ud MUST come back untouched -- the switches carry no default, and
# that silence is observable behaviour.
#
# Until stage 3 this linked DriftFluxClosure.cpp and nothing else. That stopped
# being possible when driftflux::coefficient joined the same translation unit
# (R-011 puts both namespaces in one pair): the coefficient bodies reference
# Cel, Ler and the two flow-pattern map classes, and the whole object file is
# linked whether the sweep calls into it or not. So the project objects come
# along now. What is verified is unchanged -- the aggregators are still called
# directly with the sweep's own arguments.
#
# Usage:
#   verify-aggregators.sh                 compare against the stored golden
#   verify-aggregators.sh --record        write the golden (only from verified code)
#
# Exit code: 0 when every value matches bit for bit; 1 otherwise.

set -uo pipefail
export LC_ALL=C

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "$script_dir/.." && pwd)"

GOLDEN_DIR="${MARLIM_GOLDEN_DIR:-$project_root/specs/001-refatoracao-sisprod/golden}"
GOLDEN="$GOLDEN_DIR/aggregators-baseline.txt"

red=$'\033[0;31m'; green=$'\033[0;32m'; yellow=$'\033[1;33m'; reset=$'\033[0m'

work_dir="$(mktemp -d -t marlim3-agg-XXXXXX)"
trap 'rm -rf "$work_dir"' EXIT

# -ffp-contract=off matches CMakeLists.txt:135. Without it the compiler may fuse
# a multiply and an add, which changes the result and would make this harness
# disagree with the engine for a reason that has nothing to do with the code.
objects_dir="${MARLIM_BUILD:-$project_root/build}/CMakeFiles/Marlim3.dir"
mapfile -t objects < <(find "$objects_dir/src" -name '*.o' \
                            ! -name 'Num4Main.cpp.o' \
                            ! -name 'DriftFluxClosure.cpp.o' | sort)
if (( ${#objects[@]} == 0 )); then
    printf '%sno objects under %s -- build the project first%s\n' \
           "$red" "$objects_dir" "$reset" >&2
    exit 2
fi

if ! g++ -std=c++20 -O2 -ffp-contract=off -fopenmp \
         -I"$project_root/src/include" -I"$project_root/src/thirdparty" \
         -o "$work_dir/sweep" \
         "$script_dir/sweep-aggregators.cpp" \
         "$project_root/src/core/DriftFluxClosure.cpp" \
         "${objects[@]}" -lgfortran 2> "$work_dir/build.log"; then
    printf '%sthe sweep did not compile%s\n' "$red" "$reset" >&2
    grep -v 'warning: relocation' "$work_dir/build.log" >&2
    exit 2
fi

"$work_dir/sweep" > "$work_dir/current.txt"
cases=$(wc -l < "$work_dir/current.txt")

if [[ "${1:-}" == "--record" ]]; then
    mkdir -p "$GOLDEN_DIR"
    cp "$work_dir/current.txt" "$GOLDEN"
    printf '%srecorded %s cases in %s%s\n' "$green" "$cases" "$GOLDEN" "$reset"
    printf '%sRecord only from code already proven equivalent to the baseline.%s\n' \
           "$yellow" "$reset"
    exit 0
fi

[[ -s "$GOLDEN" ]] || {
    printf '%sno golden at %s -- run --record from verified code first%s\n' \
           "$red" "$GOLDEN" "$reset" >&2
    exit 2
}

if diff -q "$GOLDEN" "$work_dir/current.txt" > /dev/null 2>&1; then
    printf '%sAGGREGATORS PASSED -- %s cases identical bit for bit%s\n' \
           "$green" "$cases" "$reset"
    awk '{ seen[$1 " selector " $2]++ } END { for (k in seen) print "  " k, seen[k] "cases" }' \
        "$GOLDEN" | sort
    exit 0
fi

printf '%sAGGREGATORS FAILED -- divergent values%s\n' "$red" "$reset" >&2
diff "$GOLDEN" "$work_dir/current.txt" | grep '^<' \
    | awk '{ if (!seen[$1 " " $2]++) printf "  first divergence: %s, selector %s\n", $1, $2 }' >&2
printf '\n  %s divergent case(s) of %s\n' \
       "$(diff "$GOLDEN" "$work_dir/current.txt" | grep -c '^<')" "$cases" >&2
exit 1
