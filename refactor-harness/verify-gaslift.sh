#!/usr/bin/env bash
# Characterization harness for the Stage 6 gas-line and gas-lift surface.
set -uo pipefail
export LC_ALL=C

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "$script_dir/.." && pwd)"
build_dir="${MARLIM_BUILD:-$project_root/build}"
objects_dir="$build_dir/CMakeFiles/Marlim3.dir"

red=$'\033[0;31m'; green=$'\033[0;32m'; reset=$'\033[0m'

mode="${1:-}"
store="${2:-}"
[[ "$mode" == "capture" || "$mode" == "compare" ]] && [[ -n "$store" ]] || {
    echo "usage: $0 {capture|compare} <file>" >&2
    exit 2
}

default_sources="$project_root/src/core/SisProd.cpp"
if [[ -f "$project_root/src/core/SisProdGasLift.cpp" ]]; then
    default_sources+=" $project_root/src/core/SisProdGasLift.cpp"
fi
read -r -a sources <<< "${MARLIM_GASLIFT_SOURCES:-$default_sources}"
for source_file in "${sources[@]}"; do
    [[ -f "$source_file" ]] || {
        printf '%sno source at %s%s\n' "$red" "$source_file" "$reset" >&2
        exit 2
    }
done

mapfile -t objects < <(find "$objects_dir/src" -name '*.o' \
    ! -name 'Num4Main.cpp.o' \
    ! -name 'SisProd.cpp.o' \
    ! -name 'SisProdGasLift.cpp.o' | sort)
(( ${#objects[@]} > 0 )) || {
    printf '%sno objects under %s -- build first%s\n' "$red" "$objects_dir" "$reset" >&2
    exit 2
}

scratch="$(mktemp -d -t marlim3-gaslift-XXXXXX)"
trap 'rm -rf "$scratch"' EXIT
flags=(-std=c++20 -O2 -Wall -ffp-contract=off -funroll-loops -fopenmp
       -I"$project_root/src/include" -I"$project_root/src/thirdparty")

g++ "${flags[@]}" -c "$script_dir/gaslift-sweep.cpp" -o "$scratch/sweep.o" \
    2> "$scratch/build.log" || {
        printf '%sthe gas-lift sweep did not compile%s\n' "$red" "$reset" >&2
        cat "$scratch/build.log" >&2
        exit 2
    }

compiled=("$scratch/sweep.o")
source_index=0
for source_file in "${sources[@]}"; do
    object_file="$scratch/module-$source_index.o"
    source_index=$((source_index + 1))
    g++ "${flags[@]}" -c "$source_file" -o "$object_file" 2>> "$scratch/build.log" || {
        printf '%s%s did not compile%s\n' "$red" "$source_file" "$reset" >&2
        tail -30 "$scratch/build.log" >&2
        exit 2
    }
    compiled+=("$object_file")
done

g++ -fopenmp -o "$scratch/gaslift-sweep" "${compiled[@]}" "${objects[@]}" -lgfortran \
    2>> "$scratch/build.log" || {
        printf '%sthe gas-lift sweep did not link%s\n' "$red" "$reset" >&2
        grep -v 'warning: relocation' "$scratch/build.log" | tail -30 >&2
        exit 2
    }

"$scratch/gaslift-sweep" > "$scratch/current.txt" || {
    printf '%sthe gas-lift sweep did not run to completion%s\n' "$red" "$reset" >&2
    exit 2
}
# The count is pinned, not just compared, because the failure this guards
# against is a row that stops being emitted -- a routine that returns early on
# seeding that drifted, printing nothing. A table that shrank would still match
# on every row it kept, and "compare" would report success.
#   40 rows: the eight routines the steady half drives, four scenarios.
#   76 rows: the six unloading routines, four scenarios, several of which
#            publish more than one row (advanceInterface prints both its
#            ordinary advance and its hand-over; the searches print the state
#            they leave behind as well as their return value).
expected_rows=116
rows=$(wc -l < "$scratch/current.txt")
(( rows == expected_rows )) || {
    printf '%sexpected %s rows, got %s%s\n' "$red" "$expected_rows" "$rows" "$reset" >&2
    exit 2
}

if [[ "$mode" == "capture" ]]; then
    mkdir -p "$(dirname "$store")"
    cp "$scratch/current.txt" "$store"
    printf '%scaptured %s gas-lift rows in %s%s\n' "$green" "$rows" "$store" "$reset"
    exit 0
fi

[[ -s "$store" ]] || {
    printf '%sno stored gas-lift table at %s%s\n' "$red" "$store" "$reset" >&2
    exit 2
}
if diff -q "$store" "$scratch/current.txt" > /dev/null; then
    printf '%sGASLIFT IDENTICAL -- %s rows%s\n' "$green" "$rows" "$reset"
    exit 0
fi

printf '%sGASLIFT DIVERGED%s\n' "$red" "$reset" >&2
diff -u "$store" "$scratch/current.txt" | head -60 >&2
exit 1
