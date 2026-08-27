#!/usr/bin/env bash
# Calibration for verify-c0ud.sh -- proves the harness catches corruption, and
# declares what it does not catch.
#
# An instrument tested only on code that passes is decoration (stage 0, E-06 and
# E-07). This script injects defects one at a time into a COPY of the five
# bodies and asserts the harness reports each one. Seven of the fourteen live in
# CalcC0UdBuf, CalcC0UdIni and CalcC0UdIniBuf, which never execute anywhere in
# the demo corpus -- those are the cases that matter, because L2 and L3 are
# green for them no matter what.
#
# Every corruption is CONFIRMED PRESENT in the copy before the harness runs. In
# stage 4 a calibration sed silently matched nothing and the test reported
# success without having tested anything (D4-06); the control below would have
# caught it, and so does the per-case check here.
#
# The control case must PASS. Without that, every "caught" line beneath it could
# be reporting a broken harness rather than a detected defect -- the silent half
# of E2-05.
#
# Usage: calibrate-c0ud.sh [golden]
set -uo pipefail
export LC_ALL=C

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "$script_dir/.." && pwd)"
golden="${1:-$project_root/specs/001-refatoracao-sisprod/golden/c0ud-baseline.txt}"

red=$'\033[0;31m'; green=$'\033[0;32m'; yellow=$'\033[1;33m'; reset=$'\033[0m'

[[ -s "$golden" ]] || { printf '%sno golden at %s%s\n' "$red" "$golden" "$reset" >&2; exit 2; }

scratch="$(mktemp -d -t marlim3-c0ud-cal-XXXXXX)"
trap 'rm -rf "$scratch"' EXIT

pristine_sisprod="$project_root/src/core/SisProd.cpp"
pristine_closure="$project_root/src/core/DriftFluxClosure.cpp"

pass=0; fail=0

run_case() {
    local name="$1" target="$2" from="$3" to="$4" expect="$5"
    local work="$scratch/$name"
    mkdir -p "$work"
    cp "$pristine_sisprod" "$work/SisProd.cpp"
    cp "$pristine_closure" "$work/DriftFluxClosure.cpp"

    if [[ -n "$from" ]]; then
        python3 - "$work/$target" "$from" "$to" <<'PY'
import sys
path, frm, to = sys.argv[1], sys.argv[2], sys.argv[3]
s = open(path).read()
n = s.count(frm)
if n == 0:
    sys.stderr.write("PATTERN NOT FOUND\n"); sys.exit(3)
open(path, "w").write(s.replace(frm, to, 1))
PY
        local rc=$?
        if (( rc != 0 )); then
            printf '  %-26s %sCORRUPTION NOT INJECTED%s -- pattern absent, case proves nothing\n' \
                   "$name" "$red" "$reset"
            fail=$((fail + 1)); return
        fi
        if cmp -s "$work/$target" "$project_root/src/core/$target"; then
            printf '  %-26s %sCORRUPTION NOT INJECTED%s -- file unchanged\n' "$name" "$red" "$reset"
            fail=$((fail + 1)); return
        fi
    fi

    local out
    out=$(MARLIM_C0UD_SOURCES="$work/SisProd.cpp $work/DriftFluxClosure.cpp" \
          bash "$script_dir/verify-c0ud.sh" compare "$golden" 2>&1)
    local verdict=$?

    if [[ "$expect" == "caught" ]]; then
        if (( verdict != 0 )); then
            printf '  %-26s %scaught%s   %s\n' "$name" "$green" "$reset" \
                   "$(grep -oE '[0-9]+ divergent row' <<< "$out" | head -1)"
            pass=$((pass + 1))
        else
            printf '  %-26s %sMISSED%s   the harness approved a corrupted body\n' \
                   "$name" "$red" "$reset"
            fail=$((fail + 1))
        fi
    else
        if (( verdict == 0 )); then
            printf '  %-26s %sinvisible%s (declared)\n' "$name" "$yellow" "$reset"
            pass=$((pass + 1))
        else
            printf '  %-26s %sunexpectedly caught%s -- update the declaration\n' \
                   "$name" "$yellow" "$reset"
            fail=$((fail + 1))
        fi
    fi
}

printf 'control -- untouched source must PASS, or nothing below means anything\n'
run_case "control" "SisProd.cpp" "" "" "invisible"
if (( fail > 0 )); then
    printf '%sthe control failed: the harness disagrees with its own golden%s\n' "$red" "$reset" >&2
    exit 1
fi

printf '\nCalcC0Ud -- executes in the corpus\n'
run_case "base/simplify-0QG+1QL" "SisProd.cpp" \
    "if (((0. * celula[ind].QG + 1 * celula[ind].QL) < 0.))" \
    "if ((celula[ind].QL < 0.))" "caught"
run_case "base/correcHor-neighbour" "SisProd.cpp" \
    "if (celula[ind - 1].acsr.tipo != 5 || celula[ind - 1].acsr.chk.AreaGarg > 1e-10) {" \
    "if (celula[ind].acsr.tipo != 5 || celula[ind].acsr.chk.AreaGarg > 1e-10) {" "caught"

printf '\nCalcC0UdBuf -- NEVER executes in the corpus\n'
run_case "buf/drop-A3-01-betI" "SisProd.cpp" \
    "        betI = celula[ind].betPigE;     // duvidabeta
        double betneg;" \
    "        double betneg;" "caught"
run_case "buf/second-phase-cond" "SisProd.cpp" \
    "            if (xarr1 != -1) {

                driftflux::correlations::C0UdDisperso(rlm, rgm, tensup1, alf0, nrey, nreyl, ug1, ul1, dia1, celula[ind].duto.rug, ang,
                             c0, ud, correcHor, celula[ind].estabCol, driftSelectors.dispersed);
                if (xarr1 == -2) {" \
    "            if (xarr1 == 1) {

                driftflux::correlations::C0UdDisperso(rlm, rgm, tensup1, alf0, nrey, nreyl, ug1, ul1, dia1, celula[ind].duto.rug, ang,
                             c0, ud, correcHor, celula[ind].estabCol, driftSelectors.dispersed);
                if (xarr1 == -2) {" "caught"
run_case "buf/pmed-weighting" "SisProd.cpp" \
    "pmed = razdx * celula[ind].presBuf + (1 - razdx) * celula[ind].presLBuf;" \
    "pmed = (1 - razdx) * celula[ind].presBuf + razdx * celula[ind].presLBuf;" "caught"

printf '\nCalcC0UdIni -- NEVER executes in the corpus\n'
run_case "ini/fix-A3-05-else" "SisProd.cpp" \
    "                            if (celula[ind].transic > 19)
                                celula[ind].transic = 0;
                        }
                    } else
                        celula[ind].transic = 0;
                    celula[ind].arranjo = xarr1 = testamapa.arr;" \
    "                            if (celula[ind].transic > 19)
                                celula[ind].transic = 0;
                        } else
                            celula[ind].transic = 0;
                    }
                    celula[ind].arranjo = xarr1 = testamapa.arr;" "caught"
run_case "ini/add-arranjoR" "SisProd.cpp" \
    "                    celula[ind].arranjo = xarr1 = testamapa.arr;
                    double c0D;" \
    "                    celula[ind].arranjo = xarr1 = testamapa.arr;
                    celula[ind - 1].arranjoR = testamapa.arr;
                    double c0D;" "caught"
run_case "ini/mapaTD-dispatch" "SisProd.cpp" \
    "                testamapa.mapaTD();
                xarr1 = testamapa.arr;
                if (xarr1 == -1) {
                    if (celula[ind].arranjo != 0) {" \
    "                testamapa.mapaTD(1);
                xarr1 = testamapa.arr;
                if (xarr1 == -1) {
                    if (celula[ind].arranjo != 0) {" "caught"

printf '\nCalcC0UdIniBuf -- NEVER executes in the corpus\n'
run_case "inibuf/add-ncel-pmed" "SisProd.cpp" \
    "        pmed = presE;
        if (ind > 0)
            pmed0 = presE;
        else
            pmed0 = presE;" \
    "        pmed = presE;
        if (ind > 0)
            pmed0 = presE;
        else
            pmed0 = presE;
        if (ind == ncel)
            pmed = celula[ind].pres;" "caught"
run_case "inibuf/buffered-source" "SisProd.cpp" \
    "        double ug1 = (celula[ind].MCBuf - celula[ind].MliqiniBuf) / rgm;
        double ul1 = celula[ind].MliqiniBuf / rlm;
        double dia1 = celula[ind].duto.a;
        if (ind > 0 && ug1 >= 0)
            dia1 = celula[ind - 1].duto.a;
        double A1 = M_PI * dia1 * dia1 / 4.;

        double rmed = hns * rlm + (1 - hns) * rgm;
        double visc = (hns * viscl1 + (1 - hns) * viscg1) / pow(10., 3.);
        double nrey = dia1 * rmed * (fabs(ug1) / A1 + fabs(ul1) / A1) / visc;
        double nreyl = dia1 * rlm * (fabs(ug1) / A1 + fabs(ul1) / A1) / (viscl1 / 1000.);

        int xarr1 = 1;
        double dtot = celula[ind].dxL + celula[ind].dx;
        double razL = celula[ind].dxL;
        double raz = celula[ind].dx;
        double ang = (raz * celula[ind].dutoL.teta + razL * celula[ind].duto.teta) / dtot;
        double atenua = 20.;" \
    "        double ug1 = (celula[ind].MC - celula[ind].Mliqini) / rgm;
        double ul1 = celula[ind].MliqiniBuf / rlm;
        double dia1 = celula[ind].duto.a;
        if (ind > 0 && ug1 >= 0)
            dia1 = celula[ind - 1].duto.a;
        double A1 = M_PI * dia1 * dia1 / 4.;

        double rmed = hns * rlm + (1 - hns) * rgm;
        double visc = (hns * viscl1 + (1 - hns) * viscg1) / pow(10., 3.);
        double nrey = dia1 * rmed * (fabs(ug1) / A1 + fabs(ul1) / A1) / visc;
        double nreyl = dia1 * rlm * (fabs(ug1) / A1 + fabs(ul1) / A1) / (viscl1 / 1000.);

        int xarr1 = 1;
        double dtot = celula[ind].dxL + celula[ind].dx;
        double razL = celula[ind].dxL;
        double raz = celula[ind].dx;
        double ang = (raz * celula[ind].dutoL.teta + razL * celula[ind].duto.teta) / dtot;
        double atenua = 20.;" "caught"

printf '\nCalcC0UdPerm -- executes in the corpus\n'
run_case "perm/fix-A3-02-ang" "SisProd.cpp" \
    "        double ang = (razL * celula[ind].dutoL.teta + raz * celula[ind].duto.teta) / dtot;
        double sinalAng = 1.;" \
    "        double ang = (raz * celula[ind].dutoL.teta + razL * celula[ind].duto.teta) / dtot;
        double sinalAng = 1.;" "caught"
run_case "perm/fix-A3-03-razdx" "SisProd.cpp" \
    "        double razdx = celula[ind].dx / (celula[ind].dx + celula[ind].dxL);" \
    "        double razdx = celula[ind].dxL / (celula[ind].dx + celula[ind].dxL);" "caught"
run_case "perm/slip-field" "SisProd.cpp" \
    "    if (arq.escorregaPerm == 0) {" "    if (arq.escorregaTran == 0) {" "caught"

printf '\ndeclared blind spots -- these MUST pass, and are covered by the token comparison\n'
# betneg feeds ul0, ul0 feeds nothing but `if (ul0 < 0.) mult0 = 0.', and mult0
# is never read. The chain betneg -> ul0 -> mult0 is dead in ALL FIVE variants
# (verified by reading every mention of the three names in each body), so no
# corruption confined to it can be observed by any behavioural harness. This is
# the same shape as A2-02 in stage 2, where the anomaly itself is what makes a
# transformation safe. Guarded by the token comparison instead.
run_case "dead-chain-0.99QG" "SisProd.cpp" \
    "if ((0.99 * celula[ind - 1].QG + 0.01 * celula[ind - 1].QL) < 0.)" \
    "if ((celula[ind - 1].QG + (celula[ind - 1].QL - celula[ind - 1].QG) * 0.01) < 0.)" "invisible"
run_case "dead-local-mult0" "SisProd.cpp" \
    "                    double mult0, mult1;
                    mult0 = 1.;
                    if (ul0 < 0.)
                        mult0 = 0.;" \
    "                    double mult0, mult1;
                    mult0 = 2.;
                    if (ul0 < 0.)
                        mult0 = 0.;" "invisible"
run_case "dead-local-timeStep" "SisProd.cpp" \
    "    int timeStep = 20;
    celula[ind].transic0 = celula[ind].transic;" \
    "    int timeStep = 21;
    celula[ind].transic0 = celula[ind].transic;" "invisible"

printf '\n'
if (( fail == 0 )); then
    printf '%sCALIBRATION PASSED -- %s case(s), 0 failure(s)%s\n' "$green" "$pass" "$reset"
    exit 0
fi
printf '%sCALIBRATION FAILED -- %s of %s case(s)%s\n' "$red" "$fail" "$((pass + fail))" "$reset" >&2
exit 1
