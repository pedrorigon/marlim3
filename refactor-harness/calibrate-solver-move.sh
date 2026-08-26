#!/usr/bin/env bash
# Calibrate solver-move.py against deliberately corrupted copies.
#
# A checker exercised only on code that passes is not a checker. This one is the
# ONLY verification that reaches zbrent and falsacorda at all: measured over the
# whole tree, SProd::zbrent has no call site and falsacorda is reachable only
# from inside it, so gates 2 and 3 are green whatever happens to those 96 lines.
# If solver-move.py is decorative, the stage has no verification for 39% of what
# it moved.
#
# The corruptions are the shapes this stage can actually produce: a changed
# literal, a reassociated expression, a flipped comparison, a swapped operand,
# an "accidental fix" of a dead branch, and a monitor that records the wrong
# quantity. The last two are the dangerous ones -- both leave every gate green.
#
# Each case asserts the corruption was really injected before running the
# checker. A sed that matched nothing produces a test that passed without
# testing anything, which happened once already while calibrating a different
# verifier.
#
# Usage: calibrate-solver-move.sh
# Exit code: 0 when the control passes and every corruption is caught.

set -uo pipefail
export LC_ALL=C

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "$script_dir/.." && pwd)"

BASELINE_COMMIT="${MARLIM_BASELINE_COMMIT:-0f3b64f}"

red=$'\033[0;31m'; green=$'\033[0;32m'; bold=$'\033[1m'; reset=$'\033[0m'

work="$(mktemp -d -t marlim3-solver-cal-XXXXXX)"
trap 'rm -rf "$work"' EXIT

git -C "$project_root" show "$BASELINE_COMMIT:src/core/SisProd.cpp" > "$work/pristine.cpp" || exit 2

failures=0

# Runs the checker over a copy of the four files, with one of them patched.
#   $1 label   $2 file to patch (header|sisprod|module)   $3 sed script
#   $4 pattern that MUST appear in the patched file, proving the sed matched
attempt() {
    local label="$1" target="$2" script="$3" proof="$4"
    local dir="$work/case"
    rm -rf "$dir"; mkdir -p "$dir"
    cp "$project_root/src/include/RootFindingSolvers.h" "$dir/RootFindingSolvers.h"
    cp "$project_root/src/core/RootFindingSolvers.cpp"   "$dir/RootFindingSolvers.cpp"
    cp "$project_root/src/core/SisProd.cpp"              "$dir/SisProd.cpp"

    local file
    case "$target" in
        header)  file="$dir/RootFindingSolvers.h" ;;
        module)  file="$dir/RootFindingSolvers.cpp" ;;
        sisprod) file="$dir/SisProd.cpp" ;;
        *) printf 'unknown target %s\n' "$target" >&2; exit 2 ;;
    esac

    if [[ -n "$script" ]]; then
        sed -i "$script" "$file"
        if ! grep -qF "$proof" "$file"; then
            printf '%s%-46s CORRUPTION NOT INJECTED%s\n' "$red" "$label" "$reset"
            printf '   the sed matched nothing; this case tested nothing\n'
            failures=$((failures + 1))
            return
        fi
    fi

    local out
    out="$(python3 "$script_dir/solver-move.py" check "$work/pristine.cpp" \
             "$dir/RootFindingSolvers.h" "$dir/SisProd.cpp" "$dir/RootFindingSolvers.cpp" 2>&1)"
    local status=$?

    if [[ -z "$script" ]]; then
        if (( status == 0 )); then
            printf '%s%-46s PASSES (control)%s\n' "$green" "$label" "$reset"
        else
            printf '%s%-46s CONTROL FAILED%s\n' "$red" "$label" "$reset"
            printf '%s\n' "$out" | sed 's/^/   /'
            failures=$((failures + 1))
        fi
        return
    fi

    if (( status != 0 )); then
        printf '%s%-46s caught%s   %s\n' "$green" "$label" "$reset" \
               "$(printf '%s\n' "$out" | grep -m1 '^FAIL' || true)"
    else
        printf '%s%-46s NOT CAUGHT%s\n' "$red" "$label" "$reset"
        failures=$((failures + 1))
    fi
}

printf '%s=== solver-move.py calibration ===%s\n\n' "$bold" "$reset"

attempt "control: sources untouched" header "" ""

# --- dead code: nothing but this checker can see these -----------------------
attempt "zbrent: guard literal 1e9 -> 1e8" header \
        's/fabs(fa) > 1e9 || fabs(fb) > 1e9/fabs(fa) > 1e8 || fabs(fb) > 1e9/' \
        'fabs(fa) > 1e8'
attempt "zbrent: reassociate 2.0 * EPS * fabs(b)" header \
        's/tol1 = 2\.0 \* EPS \* fabs(b) + 0\.5 \* tol;/tol1 = 2.0 * (EPS * fabs(b)) + 0.5 * tol;/' \
        '2.0 * (EPS * fabs(b))'
attempt "zbrent: drop the falsacorda fallback" header \
        's/val = falsacorda(x1, x2, objective);/val = 0.0;/' \
        'val = 0.0;'
attempt "falsacorda: delta 0.001 -> 0.0001" header \
        's/double delta = 0\.001;/double delta = 0.0001;/' \
        'double delta = 0.0001;'

# --- live code, but invisible to the gates ----------------------------------
attempt "zriddr: 'fix' the dead revPerm branch (A2-02)" header \
        '/^    } else {$/,/^    }$/{s/        fh = objective(x2);/        fh = objective(x1);/}' \
        'fh = objective(x1);'
attempt "zriddr: j>minit -> j>=minit" header \
        '0,/if(j>minit)return xmin;/s//if(j>=minit)return xmin;/' \
        'if(j>=minit)return xmin;'
attempt "zriddr: swap fl and fh in the sqrt" header \
        's/sqrt(fm \* fm - fl \* fh)/sqrt(fm * fm - fh * fl)/' \
        'fm * fm - fh * fl'
attempt "zriddr: xacc 1e-5 -> 1e-6" header \
        's/double xacc = 1e-5;/double xacc = 1e-6;/' \
        'double xacc = 1e-6;'

# --- the binding site: moved code that does not live in the module ----------
attempt "binding: monitor records the raw residual" sisprod \
        's/monitConvPerm = fabs(normalized);/monitConvPerm = fabs(residual);/' \
        'monitConvPerm = fabs(residual);'
attempt "binding: divide by a reciprocal instead" sisprod \
        's|double normalized = residual / monitConvPermBase;|double normalized = residual * (1.0 / monitConvPermBase);|' \
        'residual * (1.0 / monitConvPermBase)'
attempt "binding: minit 10 -> 5" sisprod \
        's/if(arq.acopColAnulPermForte == 1)minit=10;/if(arq.acopColAnulPermForte == 1)minit=5;/' \
        'minit=5;'
attempt "binding: objective drops tipoCC" sisprod \
        's/return multMarcha(guess, prod, tipoCC);/return multMarcha(guess, prod, 0);/' \
        'multMarcha(guess, prod, 0)'

# --- the helpers ------------------------------------------------------------
attempt "SIGN: b >= 0 -> b > 0" header \
        's/return (b >= 0 ? 1\.0 : -1\.0) \* fabs(a);/return (b > 0 ? 1.0 : -1.0) * fabs(a);/' \
        '(b > 0 ? 1.0 : -1.0)'
attempt "sign: var <= 0. -> var < 0." module \
        's/if (var <= 0\.)/if (var < 0.)/' \
        'if (var < 0.)'

printf '\n'
if (( failures > 0 )); then
    printf '%s%d calibration case(s) failed -- the checker is not trustworthy%s\n' \
           "$red" "$failures" "$reset" >&2
    exit 1
fi
printf '%sCALIBRATED -- control passes, every corruption caught%s\n' "$green" "$reset"
exit 0
