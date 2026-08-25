#!/usr/bin/env bash
# Dedicated harness for the eight trend writers of SisProdTrendOutput.
#
# Why this exists
# ---------------
# No model in the demo corpus declares a cross-section trend, so four of the
# eight writers are never executed by L2 or L3 (evidencia/trend-diff.md, D4-02).
# They also carry four of the five anomalies preserved from the baseline, and
# T057 restructures all eight into a Template Method. Six green gates would say
# nothing about whether that restructuring changed their output.
#
# The writers take a TrendState and write files, so they can be driven without
# an SProd. This script links trend-sweep.cpp against the project objects,
# sweeps 96 configurations, and stores the produced tree.
#
#   verify-trend-writers.sh capture <dir>   run the sweep, store the tree
#   verify-trend-writers.sh compare <dir>   run again, diff against a stored tree
#
# The module under test is COMPILED FROM SOURCE by this script, not taken from
# the project build. Set MARLIM_TREND_SOURCE to compare a variant against a
# stored tree without touching the working tree. Compiling here is also what
# keeps E1-01 from recurring: the harness cannot be fooled by an object file
# left over from an older version of the source.
#
# Exit code: 0 when the trees are byte-identical; 1 on any difference.
set -uo pipefail
export LC_ALL=C

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "$script_dir/.." && pwd)"
build_dir="${MARLIM_BUILD:-$project_root/build}"
objects_dir="$build_dir/CMakeFiles/Marlim3.dir"

red=$'\033[0;31m'; green=$'\033[0;32m'; reset=$'\033[0m'

mode="${1:-}"; store="${2:-}"
[[ -n "$mode" && -n "$store" ]] || { echo "usage: $0 {capture|compare} <dir>" >&2; exit 2; }

# Num4Main.cpp is excluded: it holds main() and the two globals the writers
# reach (arqRelatorioPerfis, pathPrefixoArqSaida), which trend-sweep.cpp defines
# itself so the sweep can point them at its own directory.
module_source="${MARLIM_TREND_SOURCE:-$project_root/src/core/SisProdTrendOutput.cpp}"
[[ -f "$module_source" ]] || { echo "${red}no module source at $module_source${reset}" >&2; exit 2; }

mapfile -t objects < <(find "$objects_dir/src" -name '*.o' \
                            ! -name 'Num4Main.cpp.o' \
                            ! -name 'SisProdTrendOutput.cpp.o' | sort)
(( ${#objects[@]} > 0 )) || { echo "${red}no objects under $objects_dir -- build first${reset}" >&2; exit 2; }

scratch="$(mktemp -d -t marlim3-trend-XXXXXX)"
trap 'rm -rf "$scratch"' EXIT

g++ -std=c++20 -O2 -ffp-contract=off -funroll-loops -Wall -fopenmp \
    -I"$project_root/src/include" -I"$project_root/src/thirdparty" \
    -c "$script_dir/trend-sweep.cpp" -o "$scratch/trend-sweep.o" || exit 2
g++ -std=c++20 -O2 -ffp-contract=off -funroll-loops -Wall -fopenmp \
    -I"$project_root/src/include" -I"$project_root/src/thirdparty" \
    -c "$module_source" -o "$scratch/module.o" || exit 2
g++ -fopenmp -o "$scratch/trend-sweep" "$scratch/trend-sweep.o" "$scratch/module.o" \
    "${objects[@]}" -lgfortran || exit 2

"$scratch/trend-sweep" "$scratch/out" > "$scratch/sweep.log" || exit 2
files=$(find "$scratch/out" -type f | wc -l)

if [[ "$mode" == "capture" ]]; then
    rm -rf "$store"; mkdir -p "$(dirname "$store")"; cp -r "$scratch/out" "$store"
    printf '%scaptured %s files from %s of %s%s\n' \
           "$green" "$files" "$(cat "$scratch/sweep.log")" "$module_source" "$reset"
    exit 0
fi

[[ -d "$store" ]] || { echo "${red}no stored tree at $store${reset}" >&2; exit 2; }
if diff -r "$store" "$scratch/out" > "$scratch/diff.txt" 2>&1; then
    printf '%sTREND WRITERS IDENTICAL -- %s files, %s of %s%s\n' \
           "$green" "$files" "$(cat "$scratch/sweep.log")" "$module_source" "$reset"
    exit 0
fi

printf '%sTREND WRITERS DIVERGED%s\n' "$red" "$reset" >&2
head -40 "$scratch/diff.txt" >&2
printf '%s%s differing path(s)%s\n' "$red" "$(grep -c '^' "$scratch/diff.txt")" "$reset" >&2
exit 1
