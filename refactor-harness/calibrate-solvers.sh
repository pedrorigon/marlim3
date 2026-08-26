#!/usr/bin/env bash
# Calibrate verify-solvers.sh, and pin down what it CANNOT see.
#
# The sweep guards the solvers once solver-move.py stops applying -- that one
# compares token streams against the pristine commit, so it fails by design as
# soon as T037r renames a local.
#
# Two groups of cases, and the second matters as much as the first.
#
#   caught  -- corruptions that change behaviour. Some change only the path to
#              the root and not the root itself, which is why the reference
#              records every evaluation rather than the returned value.
#
#   passes  -- rewrites the sweep correctly does NOT flag, each paired with the
#              checker that does. They are not all the same kind:
#
#              Three are identical by construction. Multiplication commutes in
#              IEEE 754; `e *= k` is defined as `e = e * k`; and regrouping
#              `2.0 * x * y` is exact because scaling by a power of two commutes
#              with rounding. There is no behaviour to see, and no checker that
#              works by running the code could see it.
#
#              One is not. Reassociating (q-1)*(r-1)*(s-1) genuinely CAN change
#              the last bit -- it just does not, for any triple this sweep
#              produces. Verified rather than assumed: perturbing the same line
#              to (s-1.5) is caught in 194 lines, and perturbing `d = p / q` is
#              caught in 248, so the interpolation branch runs and q's value
#              reaches the result. The rewrite is invisible for these inputs and
#              would not be for others. That is the E4-01 question -- "for which
#              INPUTS is this difference observable?" -- with the answer written
#              down instead of guessed at.
#
#              One is the dangerous one: a dead branch whose arms are identical
#              can be collapsed without changing any result, and only a token
#              comparison notices. This is the shape an accidental "fix" of
#              A2-02 would take.
#
# Writing that second group down is the point. A harness whose calibration lists
# only its successes reads as though it covers everything.
#
# Usage: calibrate-solvers.sh
# Exit code: 0 when every case behaves as declared.

set -uo pipefail
export LC_ALL=C

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "$script_dir/.." && pwd)"
reference="$project_root/specs/001-refatoracao-sisprod/golden/solver-sweep.txt"
baseline_commit="${MARLIM_BASELINE_COMMIT:-0f3b64f}"

red=$'\033[0;31m'; green=$'\033[0;32m'; yellow=$'\033[1;33m'; bold=$'\033[1m'; reset=$'\033[0m'

work="$(mktemp -d -t marlim3-solvers-cal-XXXXXX)"
trap 'rm -rf "$work"' EXIT

[[ -s "$reference" ]] || {
    printf '%sno reference at %s -- capture it first%s\n' "$red" "$reference" "$reset" >&2
    exit 2
}
git -C "$project_root" show "$baseline_commit:src/core/SisProd.cpp" > "$work/pristine.cpp"

failures=0

# Applies $2 (a python snippet given the header path as sys.argv[1]) to a copy
# of the header, then runs the sweep. $1 label, $3 expectation: caught|passes.
attempt() {
    local label="$1" patch="$2" expect="$3"
    local inc="$work/include"
    rm -rf "$inc"; mkdir -p "$inc"
    cp "$project_root/src/include/RootFindingSolvers.h" "$inc/RootFindingSolvers.h"

    if [[ -n "$patch" ]]; then
        if ! python3 -c "$patch" "$inc/RootFindingSolvers.h"; then
            printf '%s%-52s CORRUPTION NOT INJECTED%s\n' "$red" "$label" "$reset"
            failures=$((failures + 1))
            return
        fi
    fi

    local out status
    out="$(MARLIM_SOLVER_INCLUDE="$inc" bash "$script_dir/verify-solvers.sh" \
             compare "$reference" 2>&1)"
    status=$?

    local detail
    detail="$(printf '%s\n' "$out" | grep -m1 -E 'differing line|diagnostics|failed to' || true)"

    case "$expect:$status" in
        control:0)
            printf '%s%-52s PASSES (control)%s\n' "$green" "$label" "$reset" ;;
        control:*)
            printf '%s%-52s CONTROL FAILED%s\n' "$red" "$label" "$reset"
            printf '%s\n' "$out" | sed 's/^/   /'
            failures=$((failures + 1)) ;;
        caught:0)
            printf '%s%-52s NOT CAUGHT%s\n' "$red" "$label" "$reset"
            failures=$((failures + 1)) ;;
        caught:*)
            printf '%s%-52s caught%s   %s\n' "$green" "$label" "$reset" "$detail" ;;
        passes:0)
            # Declared invisible to the sweep. Prove the token comparison sees it.
            local token_out token_status
            token_out="$(python3 "$script_dir/solver-move.py" check "$work/pristine.cpp" \
                          "$inc/RootFindingSolvers.h" \
                          "$project_root/src/core/SisProd.cpp" \
                          "$project_root/src/core/RootFindingSolvers.cpp" 2>&1)"
            token_status=$?
            if (( token_status != 0 )); then
                printf '%s%-52s invisible here, caught by tokens%s   %s\n' \
                       "$green" "$label" "$reset" \
                       "$(printf '%s\n' "$token_out" | grep -m1 '^FAIL' || true)"
            else
                printf '%s%-52s INVISIBLE TO BOTH CHECKERS%s\n' "$red" "$label" "$reset"
                failures=$((failures + 1))
            fi ;;
        passes:*)
            printf '%s%-52s FLAGGED, but declared invisible%s\n' "$red" "$label" "$reset"
            printf '%s\n' "$out" | sed 's/^/   /'
            failures=$((failures + 1)) ;;
    esac
}

# Helper: a python snippet that replaces `old` with `new` once and fails loudly
# when it matched nothing, so a case can never pass without having been injected.
swap() {
    local times="${3:-1}"
    printf 'import sys\nold, new = %s, %s\np = sys.argv[1]\nt = open(p).read()\nassert t.count(old) == %s, "matched %%d times, expected %s" %% t.count(old)\nopen(p, "w").write(t.replace(old, new))\n' \
           "$1" "$2" "$times" "$times"
}

printf '%s=== verify-solvers.sh calibration ===%s\n\n' "$bold" "$reset"

attempt "control: header untouched" "" control

printf '\n%sbehaviour changed -- the sweep must catch these%s\n' "$bold" "$reset"

attempt "zbrent: 0.5 * tol becomes 0.25 * tol" \
        "$(swap "'+ 0.5 * tol;'" "'+ 0.25 * tol;'")" caught
attempt "zbrent: q and r swapped in the interpolation" \
        "$(swap "'q = fa / fc;\n                    r = fb / fc;'" "'q = fb / fc;\n                    r = fa / fc;'")" caught
attempt "falsacorda: bracket arms swapped" \
        "$(swap "'? (b = c) : (a = c, u = w);'" "'? (a = c, u = w) : (b = c);'")" caught
attempt "falsacorda: bisection midpoint drifts" \
        "$(swap "'c = a + e;'" "'c = a + e * 1.0000001;'")" caught
attempt "zriddr: xacc 1e-5 -> 2e-5 (path, not root)" \
        "$(swap "'double xacc = 1e-5;'" "'double xacc = 2e-5;'")" caught
attempt "zbrent: iteration budget starts at 1" \
        "$(swap "'for (int iter = 0; iter < maxit; iter++)'" "'for (int iter = 1; iter < maxit; iter++)'")" caught
attempt "zriddr: A2-05 guard removed, division unreachable" \
        "$(swap "'                if(j>minit)return xmin;\n            }\n            double xnew'" "'                return xmin;\n            }\n            double xnew'")" caught
attempt "zriddr: guard 1e9 -> 1e11" \
        "$(swap "'if (fabs(fl) > 1e9 || fabs(fh) > 1e9)'" "'if (fabs(fl) > 1e11 || fabs(fh) > 1e11)'")" caught
attempt "zriddr: widening step 1.0001 -> 1.001, both arms" \
        "$(swap "'                        x1 *= 1.0001;'" "'                        x1 *= 1.001;'" 2)" caught
attempt "SIGN: copysign instead of the comparison" \
        "$(swap "'return (b >= 0 ? 1.0 : -1.0) * fabs(a);'" "'return copysign(fabs(a), b);'")" caught

printf '\n%sinvisible to the sweep by construction -- tokens must catch these%s\n' \
       "$yellow" "$reset"

attempt "multiplication commuted (bit-identical in IEEE)" \
        "$(swap "'(b - a) * (r - 1.0)'" "'(r - 1.0) * (b - a)'")" passes
attempt "2.0 * x * y regrouped (2.0 is a power of two)" \
        "$(swap "'2.0 * xm * q * (q - r)'" "'2.0 * (xm * q * (q - r))'")" passes
attempt "e *= 0.5 spelled out (identical by definition)" \
        "$(swap "'        e *= 0.5;'" "'        e = e * 0.5;'")" passes
attempt "reassociated (q-1)*(r-1)*(s-1) -- see note" \
        "$(swap "'q = (q - 1.0) * (r - 1.0) * (s - 1.0);'" "'q = (q - 1.0) * ((r - 1.0) * (s - 1.0));'")" passes
attempt "A2-02 dead branch collapsed (arms are equal)" \
        "$(swap "'    if (revPerm == 0) {\n        fl = objective(x1);\n        fh = objective(x2);\n    } else {\n        fl = objective(x1);\n        fh = objective(x2);\n    }'" "'    fl = objective(x1);\n    fh = objective(x2);'")" passes

printf '\n'
if (( failures > 0 )); then
    printf '%s%d calibration case(s) behaved unexpectedly%s\n' "$red" "$failures" "$reset" >&2
    exit 1
fi
printf '%sCALIBRATED -- 10 corruptions caught; 5 declared invisible, each caught by tokens%s\n' \
       "$green" "$reset"
exit 0
