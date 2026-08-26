#!/usr/bin/env bash
# Behavioural sweep of the three root-finding solvers.
#
# What it adds over solver-move.py, which already proves the move was literal:
#
#   * It INSTANTIATES the templates. zbrent has no call site in the project and
#     falsacorda is reachable only from inside it, so after the extraction
#     nothing in the product build instantiates either -- and an uninstantiated
#     template is parsed, not type-checked. Without this the compiler never
#     looks at 96 of the 240 moved lines.
#
#   * It becomes the reference for every step AFTER the move. solver-move.py
#     compares token streams against the pristine baseline, so it fails by
#     design the moment T037r renames a local. From there on this sweep is what
#     says the behaviour did not move with the names.
#
#   * It records the whole evaluation sequence, not the returned root. The
#     contract is explicit that comparing roots proves nothing: two convergence
#     paths can reach the same value.
#
# The reference is captured immediately after the move, while solver-move.py
# still reports exact token equality against the pristine commit -- which is
# what makes a reference captured from the current tree mean anything.
#
# Usage:
#   verify-solvers.sh capture [<file>]   record the reference
#   verify-solvers.sh compare [<file>]   compare against it
#
# Exit code: 0 when the sweep matches the reference.

set -uo pipefail
export LC_ALL=C

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "$script_dir/.." && pwd)"

GOLDEN_DIR="${MARLIM_GOLDEN_DIR:-$project_root/specs/001-refatoracao-sisprod/golden}"
mode="${1:?usage: verify-solvers.sh capture|compare [file]}"
reference="${2:-$GOLDEN_DIR/solver-sweep.txt}"

# Where to take RootFindingSolvers.h from. Overridden by calibrate-solvers.sh so
# corruptions can be injected into a copy instead of the working tree.
SOLVER_INCLUDE="${MARLIM_SOLVER_INCLUDE:-$project_root/src/include}"

red=$'\033[0;31m'; green=$'\033[0;32m'; yellow=$'\033[1;33m'; reset=$'\033[0m'

work="$(mktemp -d -t marlim3-solvers-XXXXXX)"
trap 'rm -rf "$work"' EXIT

# reportIterationLimit forwards to NumError, which lives in
# FerramentasNumericas.cpp and drags the whole engine behind it. A stub keeps
# the sweep standalone while still linking the real forwarder and the real
# sign(), which is what is under test here.
cat > "$work/numerror-stub.cpp" <<'STUB'
#include <string>
void NumError(const std::string &) {}
STUB

# -Wall and the same -ffp-contract=off the product build uses. Contraction would
# fold a*b+c into an fma and change the last bit, which is exactly the class of
# difference this whole programme exists to prevent.
if ! g++ -std=c++20 -Wall -O2 -ffp-contract=off \
        -I"$SOLVER_INCLUDE" -I"$project_root/src/include" \
        -I"$project_root/src/thirdparty" \
        "$script_dir/solver-sweep.cpp" \
        "$project_root/src/core/RootFindingSolvers.cpp" \
        "$work/numerror-stub.cpp" \
        -o "$work/solver-sweep" 2> "$work/build.log"; then
    printf '%sthe sweep driver failed to compile%s\n' "$red" "$reset" >&2
    cat "$work/build.log" >&2
    exit 2
fi

# Two diagnostics are expected, and both are inherited: the baseline build log
# has them at src/core/SisProd.cpp:25654 and :25622, inside falsacorda and
# zbrent as they stood before the move. They vanish from the product build only
# because nothing there instantiates those templates any more -- this driver is
# now the only place they still surface, which is a reason to keep them visible
# rather than silence them with -Wno-.
#
# Compared as counts PER FLAG, not by message text. The first version of this
# check listed the messages, including "'e' may be used uninitialized" -- and
# T037r renamed e to previousStep, so the diagnostic changed text without
# changing meaning and the check failed on a rename it was never meant to
# police. That is D1-10 again, at a smaller scale: an identity that anchors on
# something the work legitimately changes.
expected_warning_flags() {
    cat <<'EXPECTED'
[-Wmaybe-uninitialized] 1
[-Wunused-variable] 1
EXPECTED
}

# GCC quotes identifiers typographically under a UTF-8 locale and in ASCII under
# LC_ALL=C; normalize before comparing, as gate 1 learned to (E-04).
observed_warning_flags() {
    sed -e "s/\xe2\x80\x98/'/g; s/\xe2\x80\x99/'/g" "$work/build.log" \
        | grep "warning:" \
        | sed -E 's/^.*(\[-W[a-z-]+\]).*/\1/; t; s/.*/[no-flag]/' \
        | sort | uniq -c | awk '{ printf "%s %s\n", $2, $1 }' | sort
}

if ! diff <(expected_warning_flags | sort) <(observed_warning_flags) > /dev/null 2>&1; then
    printf '%sthe sweep driver'"'"'s diagnostics are not the two inherited ones:%s\n' \
           "$red" "$reset" >&2
    diff <(expected_warning_flags | sort) <(observed_warning_flags) | sed 's/^/  /' >&2
    printf '%sfull output:%s\n' "$yellow" "$reset" >&2
    grep "warning:" "$work/build.log" | sed 's/^/  /' >&2
    exit 2
fi

"$work/solver-sweep" "$work/sweep.txt" > "$work/run.log" 2>&1 || {
    printf '%sthe sweep driver failed to run%s\n' "$red" "$reset" >&2
    cat "$work/run.log" >&2
    exit 2
}
cat "$work/run.log"

case "$mode" in
capture)
    mkdir -p "$(dirname "$reference")"
    cp "$work/sweep.txt" "$reference"
    printf '%scaptured %s lines to %s%s\n' \
           "$green" "$(wc -l < "$reference")" "$reference" "$reset"
    exit 0
    ;;
compare)
    [[ -s "$reference" ]] || {
        printf '%sno reference at %s -- run `verify-solvers.sh capture` first%s\n' \
               "$red" "$reference" "$reset" >&2
        exit 2
    }
    if diff -q "$reference" "$work/sweep.txt" > /dev/null 2>&1; then
        printf '%sSOLVERS IDENTICAL -- %s lines, evaluation sequences included%s\n' \
               "$green" "$(wc -l < "$reference")" "$reset"
        exit 0
    fi
    printf '%sSOLVERS DIVERGED%s\n' "$red" "$reset" >&2
    diff "$reference" "$work/sweep.txt" | head -30 >&2
    printf '\n  %s differing line(s) of %s\n' \
           "$(diff "$reference" "$work/sweep.txt" | grep -c '^<')" \
           "$(wc -l < "$reference")" >&2
    exit 1
    ;;
*)
    printf 'unknown mode %s -- use capture or compare\n' "$mode" >&2
    exit 2
    ;;
esac
