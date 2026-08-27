#!/usr/bin/env bash
# Dedicated harness for the five CalcC0Ud* variants of stage 3.
#
# Why this exists
# ---------------
# gcov over the full corpus (stage 0, evidencia/cobertura.md) measured that only
# CalcC0Ud and CalcC0UdPerm ever execute. CalcC0UdBuf, CalcC0UdIni and
# CalcC0UdIniBuf never run anywhere in the fourteen models: 718 of the 1248
# lines and 144 of the 248 conditionals are invisible to L2 and to L3. Stage 3
# restructures all five. Six green gates would say nothing about three of them.
#
# c0ud-sweep.cpp drives the five through their PUBLIC signatures, which the
# stage preserves (FR-032). That is what makes one capture comparable across the
# whole stage: the same table is produced before the extraction, after it, and
# after the unification. A harness written against the new internal API could
# only be calibrated after the change it exists to police.
#
# WHAT IT RECORDS
#   c0, ud                                     the two outputs, in %a
#   celula[ind].arranjo .transic .transic0     the pattern bookkeeping
#   celula[ind].c0Spare .udSpare               the saved pair
#   celula[ind-1].arranjoR .perdaEstratL .perdaEstratG
#
# Recording c0 and ud alone would be worthless: every variant ends with
#   if (arq.escorregaTran == 0) { ...; c0 = 1. + 0*correcaoCo; ud = 0. + ...; }
# which discards everything above it, so with slip off the table is constant.
# The six write-only fields are seeded with sentinels, never zero, so that "the
# write did not happen" is distinguishable from "the write stored a zero" --
# the stratified block stores fatorperdaLiq == 0 on most inputs, and with zero
# seeds the sweep reported it untouched while it was in fact running.
#
# The module under test is COMPILED FROM SOURCE here, never taken from the
# project build, so a stale object cannot make the harness agree with itself
# (E1-01). Set MARLIM_C0UD_SOURCES to a space-separated list to run a variant --
# a deliberately corrupted copy, say -- without touching the working tree.
#
#   verify-c0ud.sh capture <file>   run the sweep, store the table
#   verify-c0ud.sh compare <file>   run again, diff against a stored table
#
# Exit code: 0 when the tables are identical; 1 on any difference.
set -uo pipefail
export LC_ALL=C

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "$script_dir/.." && pwd)"
build_dir="${MARLIM_BUILD:-$project_root/build}"
objects_dir="$build_dir/CMakeFiles/Marlim3.dir"

red=$'\033[0;31m'; green=$'\033[0;32m'; yellow=$'\033[1;33m'; reset=$'\033[0m'

mode="${1:-}"; store="${2:-}"
[[ -n "$mode" && -n "$store" ]] || { echo "usage: $0 {capture|compare} <file>" >&2; exit 2; }
[[ "$mode" == "capture" || "$mode" == "compare" ]] || { echo "usage: $0 {capture|compare} <file>" >&2; exit 2; }

# The sources that hold the five bodies. Both are compiled here and both are
# excluded from the object list, so whichever file the bodies currently live in
# is the one under test. DriftFluxClosure.cpp is listed from the start: before
# the extraction it holds only the correlations, and compiling it here as well
# keeps the command identical on both sides of the move.
default_sources="$project_root/src/core/SisProd.cpp $project_root/src/core/DriftFluxClosure.cpp"
read -r -a sources <<< "${MARLIM_C0UD_SOURCES:-$default_sources}"
for src in "${sources[@]}"; do
    [[ -f "$src" ]] || { printf '%sno source at %s%s\n' "$red" "$src" "$reset" >&2; exit 2; }
done

mapfile -t objects < <(find "$objects_dir/src" -name '*.o' \
                            ! -name 'Num4Main.cpp.o' \
                            ! -name 'SisProd.cpp.o' \
                            ! -name 'DriftFluxClosure.cpp.o' | sort)
(( ${#objects[@]} > 0 )) || {
    printf '%sno objects under %s -- build the project first%s\n' "$red" "$objects_dir" "$reset" >&2
    exit 2
}

scratch="$(mktemp -d -t marlim3-c0ud-XXXXXX)"
trap 'rm -rf "$scratch"' EXIT

# -ffp-contract=off matches CMakeLists.txt:135; without it the compiler may fuse
# a multiply and an add and the harness would disagree with the engine for a
# reason that has nothing to do with this stage.
flags=(-std=c++20 -O2 -ffp-contract=off -funroll-loops -fopenmp
       -I"$project_root/src/include" -I"$project_root/src/thirdparty")

compiled=()
g++ "${flags[@]}" -c "$script_dir/c0ud-sweep.cpp" -o "$scratch/c0ud-sweep.o" \
    2> "$scratch/build.log" || { printf '%sthe sweep did not compile%s\n' "$red" "$reset" >&2
                                 cat "$scratch/build.log" >&2; exit 2; }
compiled+=("$scratch/c0ud-sweep.o")
i=0
for src in "${sources[@]}"; do
    obj="$scratch/module-$i.o"; i=$((i + 1))
    g++ "${flags[@]}" -c "$src" -o "$obj" 2>> "$scratch/build.log" || {
        printf '%s%s did not compile%s\n' "$red" "$src" "$reset" >&2
        tail -30 "$scratch/build.log" >&2; exit 2; }
    compiled+=("$obj")
done

g++ -fopenmp -o "$scratch/c0ud-sweep" "${compiled[@]}" "${objects[@]}" -lgfortran \
    2>> "$scratch/build.log" || { printf '%sthe sweep did not link%s\n' "$red" "$reset" >&2
                                  grep -v 'warning: relocation' "$scratch/build.log" | tail -20 >&2; exit 2; }

"$scratch/c0ud-sweep" > "$scratch/current.txt" || {
    printf '%sthe sweep did not run to completion%s\n' "$red" "$reset" >&2; exit 2; }
rows=$(wc -l < "$scratch/current.txt")
(( rows > 0 )) || { printf '%sthe sweep produced no rows%s\n' "$red" "$reset" >&2; exit 2; }

# A table whose recorded columns are all identical would compare equal against
# anything and is not a table. This is the E4-01 guard: refuse to act as an
# authority when the sweep did not actually vary its outputs.
distinct=$(cut -d' ' -f3- "$scratch/current.txt" | sort -u | wc -l)
(( distinct > rows / 4 )) || {
    printf '%sonly %s distinct outcomes in %s rows -- the sweep is not discriminating%s\n' \
           "$red" "$distinct" "$rows" "$reset" >&2; exit 2; }

if [[ "$mode" == "capture" ]]; then
    mkdir -p "$(dirname "$store")"
    cp "$scratch/current.txt" "$store"
    printf '%scaptured %s rows (%s distinct) in %s%s\n' \
           "$green" "$rows" "$distinct" "$store" "$reset"
    printf '%sCapture only from code already proven equivalent to the baseline.%s\n' \
           "$yellow" "$reset"
    exit 0
fi

[[ -s "$store" ]] || { printf '%sno stored table at %s%s\n' "$red" "$store" "$reset" >&2; exit 2; }

if diff -q "$store" "$scratch/current.txt" > /dev/null 2>&1; then
    printf '%sC0UD IDENTICAL -- %s rows, %s distinct outcomes%s\n' \
           "$green" "$rows" "$distinct" "$reset"
    awk '{ seen[$1]++ } END { for (k in seen) printf "  %-15s %d rows\n", k, seen[k] }' \
        "$store" | sort
    exit 0
fi

printf '%sC0UD DIVERGED%s\n' "$red" "$reset" >&2
differing=$(diff "$store" "$scratch/current.txt" | grep -c '^<')
diff "$store" "$scratch/current.txt" | grep '^<' \
    | awk '{ if (!seen[$1]++) printf "  first divergence in %-15s scenario %s %s\n", $1, $2, $3 }' >&2
printf '\n  %s divergent row(s) of %s\n' "$differing" "$rows" >&2
diff "$store" "$scratch/current.txt" | grep '^<' \
    | awk '{ print $1 }' | sort | uniq -c | awk '{ printf "  %-6s rows in %s\n", $1, $2 }' >&2
exit 1
