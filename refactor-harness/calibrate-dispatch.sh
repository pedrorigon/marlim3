#!/usr/bin/env bash
# Calibrate verify-dispatch.py against deliberately mis-wired dispatch tables.
#
# The sweep it runs is the only check on C-C6. Six of the nine rows are reached
# by the corpus; the injection row and the two gas-line rows are reached by no
# model, so a row wired to the wrong march there passes every gate. If this
# checker is decorative, a third of the table is unverified.
#
# Two of the cases are not about wiring at all. The selector takes the choke
# ARRAY, not its first element, so that the subscript happens only in the branch
# that needs it -- as it did when this was a nested chain. Ler::copia_chokeSup
# leaves chokep.abertura null when parserie is not positive, and reading it on
# every dispatch would turn a conditional dereference into an unconditional one.
# No corpus model takes that path, so no gate would say anything. Both the eager
# call site and the eager selector are injected here and must be caught. The
# first of the two is not hypothetical: it is what this stage wrote first.
#
# One case is expected NOT to be caught, and is asserted to pass rather than
# quietly skipped: moving the choke-opening threshold. Both arms of that branch
# select the same march (A2-01), so no value of the threshold changes the
# result. Recording that as a known blind spot is the honest form; leaving it
# out would suggest the checker sees more than it does.
#
# Usage: calibrate-dispatch.sh
# Exit code: 0 when every case behaves as declared.

set -uo pipefail
export LC_ALL=C

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "$script_dir/.." && pwd)"

red=$'\033[0;31m'; green=$'\033[0;32m'; yellow=$'\033[1;33m'; bold=$'\033[1m'; reset=$'\033[0m'

work="$(mktemp -d -t marlim3-dispatch-cal-XXXXXX)"
trap 'rm -rf "$work"' EXIT

failures=0

#   $1 label   $2 expectation: caught|passes   $3 sed script   $4 proof pattern
attempt() {
    local label="$1" expect="$2" script="$3" proof="$4"
    local file="$work/SisProd.cpp"
    cp "$project_root/src/core/SisProd.cpp" "$file"

    if [[ -n "$script" ]]; then
        sed -i "$script" "$file"
        if ! grep -qF "$proof" "$file"; then
            printf '%s%-50s CORRUPTION NOT INJECTED%s\n' "$red" "$label" "$reset"
            failures=$((failures + 1))
            return
        fi
    fi

    local out status
    out="$(python3 "$script_dir/verify-dispatch.py" --current "$file" 2>&1)"
    status=$?

    case "$expect:$status" in
        caught:0)
            printf '%s%-50s NOT CAUGHT%s\n' "$red" "$label" "$reset"
            failures=$((failures + 1)) ;;
        caught:*)
            printf '%s%-50s caught%s   %s\n' "$green" "$label" "$reset" \
                   "$(printf '%s\n' "$out" | grep -m1 -E '^(DIFF|baseline rule|.*still refers)' || true)" ;;
        passes:0)
            printf '%s%-50s passes%s\n' "$green" "$label" "$reset" ;;
        passes:*)
            printf '%s%-50s FAILED, but should pass%s\n' "$red" "$label" "$reset"
            printf '%s\n' "$out" | sed 's/^/   /'
            failures=$((failures + 1)) ;;
    esac
}

printf '%s=== verify-dispatch.py calibration ===%s\n\n' "$bold" "$reset"

attempt "control: source untouched" passes "" ""

attempt "row swapped: gasPerm2 <-> gasPerm3" caught \
        's/return tipoCC == 0 ? SteadyMarch::marchaGasPerm2 : SteadyMarch::marchaGasPerm3;/return tipoCC == 0 ? SteadyMarch::marchaGasPerm3 : SteadyMarch::marchaGasPerm2;/' \
        'tipoCC == 0 ? SteadyMarch::marchaGasPerm3'

attempt "A2-01 'fixed': else arm goes elsewhere" caught \
        's/                                              : SteadyMarch::marchaProdPresPres2;/                                              : SteadyMarch::marchaProdPresPres1;/' \
        ': SteadyMarch::marchaProdPresPres1;'

attempt "injection guard: != 0 becomes > 0" caught \
        's/    if (injectorWell != 0)/    if (injectorWell > 0)/' \
        'if (injectorWell > 0)'

attempt "prod guard: == 1 becomes >= 1" caught \
        's/    if (prod == 1) {/    if (prod >= 1) {/' \
        'if (prod >= 1) {'

attempt "reverse march inverted in ProdPerm1" caught \
        's/        return reverseMarch == 0 ? SteadyMarch::marchaProdPerm1$/        return reverseMarch != 0 ? SteadyMarch::marchaProdPerm1/' \
        'reverseMarch != 0 ? SteadyMarch::marchaProdPerm1'

attempt "reverse march inverted in ProdPresPres1" caught \
        's/    return reverseMarch == 0 ? SteadyMarch::marchaProdPresPres1$/    return reverseMarch != 0 ? SteadyMarch::marchaProdPresPres1/' \
        'reverseMarch != 0 ? SteadyMarch::marchaProdPresPres1'

attempt "call site subscripts the choke array eagerly" caught \
        's/revPerm, arq\.chokep\.abertura))/revPerm, arq.chokep.abertura[0]))/' \
        'arq.chokep.abertura[0]))'

attempt "selector subscripts the choke array eagerly" caught \
        's|^    if (injectorWell != 0)$|    const double eager = productionChokeOpening[0]; (void)eager;\n    if (injectorWell != 0)|' \
        'const double eager = productionChokeOpening[0];'

printf '\n%sdeclared blind spot%s\n' "$yellow" "$reset"
attempt "choke threshold 1e-15 -> 1e-9 (A2-01 arms equal)" passes \
        's/productionChokeOpening\[0\] > 1e-15/productionChokeOpening[0] > 1e-9/' \
        'productionChokeOpening[0] > 1e-9'

printf '\n'
if (( failures > 0 )); then
    printf '%s%d calibration case(s) behaved unexpectedly%s\n' "$red" "$failures" "$reset" >&2
    exit 1
fi
printf '%sCALIBRATED -- eight defects caught, one declared blind spot confirmed%s\n' \
       "$green" "$reset"
exit 0
