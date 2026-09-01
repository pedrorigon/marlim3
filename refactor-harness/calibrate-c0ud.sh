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
# THE CASES TARGET src/core/DriftFluxClosure.cpp, NOT SisProd.cpp.
#
# They used to target SisProd.cpp, because that is where the five bodies lived
# when this file was written. Re-running the calibration at the END of stage 3 --
# which is the rule, not a formality -- reported CORRUPTION NOT INJECTED fourteen
# times: the bodies had moved, folded onto shared helpers, and had 229 locals
# renamed, so not one of the patterns still matched. Nothing was silently
# approved, because a case must prove the defect is PRESENT before running the
# harness; without that guard this file would have printed fourteen greens
# having tested nothing.
#
# The lesson generalises past this stage: a calibration is written against a
# particular text, and every stage that rewrites that text invalidates it. When
# stage 5 or 7 moves this code again, these patterns break again, and the fix is
# to re-point them -- never to relax the guard that noticed.
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
run_case "control" "DriftFluxClosure.cpp" "" "" "invisible"
if (( fail > 0 )); then
    printf '%sthe control failed: the harness disagrees with its own golden%s\n' "$red" "$reset" >&2
    exit 1
fi

printf '\ninstantaneous (was CalcC0Ud) -- executes in the corpus\n'
run_case "inst/simplify-0QG+1QL" "DriftFluxClosure.cpp" \
    "if (((0. * state.cells[cellIndex].QG + 1 * state.cells[cellIndex].QL) < 0.))" \
    "if ((state.cells[cellIndex].QL < 0.))" "caught"
run_case "inst/correcHor-neighbour" "DriftFluxClosure.cpp" \
    "if (state.cells[accessoryCellIndex].acsr.tipo != 5 || state.cells[accessoryCellIndex].acsr.chk.AreaGarg > 1e-10) {" \
    "if (state.cells[cellIndex].acsr.tipo != 5 || state.cells[cellIndex].acsr.chk.AreaGarg > 1e-10) {" "caught"
run_case "inst/map-dispatch" "DriftFluxClosure.cpp" \
    "                if (state.selectors.stratified == 2)
                    stratifiedMap.mapaTD();
                else
                    stratifiedMap.mapaTD(1);" \
    "                if (state.selectors.stratified == 2)
                    stratifiedMap.mapaTD(1);
                else
                    stratifiedMap.mapaTD();" "caught"

printf '\nbuffered (was CalcC0UdBuf) -- NEVER executes in the corpus\n'
run_case "buf/drop-A3-01-betI" "DriftFluxClosure.cpp" \
    "        betI = state.cells[cellIndex].betPigE;     // duvidabeta
        double betneg;
        if (cellIndex > 0) {
            betneg = state.cells[cellIndex - 1].betL;" \
    "        double betneg;
        if (cellIndex > 0) {
            betneg = state.cells[cellIndex - 1].betL;" "caught"
run_case "buf/second-phase-cond" "DriftFluxClosure.cpp" \
    "            if (flowPattern != -1) {" \
    "            if (flowPattern == 1) {" "caught"
run_case "buf/pmed-weighting" "DriftFluxClosure.cpp" \
    "        meanPressure = lengthRatio * state.cells[cellIndex].presBuf + (1 - lengthRatio) * state.cells[cellIndex].presLBuf;" \
    "        meanPressure = (1 - lengthRatio) * state.cells[cellIndex].presBuf + lengthRatio * state.cells[cellIndex].presLBuf;" "caught"

printf '\ninitialization (was CalcC0UdIni) -- NEVER executes in the corpus\n'
run_case "ini/fix-A3-05-else" "DriftFluxClosure.cpp" \
    "                    } else
                        state.cells[cellIndex].transic = 0;
                    state.cells[cellIndex].arranjo = flowPattern = stratifiedMap.arr;" \
    "                    } else
                        state.cells[cellIndex].transic = 1;
                    state.cells[cellIndex].arranjo = flowPattern = stratifiedMap.arr;" "caught"
run_case "ini/add-arranjoR" "DriftFluxClosure.cpp" \
    "                    state.cells[cellIndex].arranjo = flowPattern = stratifiedMap.arr;
                    RegimePair pair;" \
    "                    state.cells[cellIndex].arranjo = flowPattern = stratifiedMap.arr;
                    state.cells[cellIndex - 1].arranjoR = stratifiedMap.arr;
                    RegimePair pair;" "caught"
run_case "ini/mapaTD-dispatch" "DriftFluxClosure.cpp" \
    "                stratifiedMap.mapaTD();
                flowPattern = stratifiedMap.arr;" \
    "                stratifiedMap.mapaTD(1);
                flowPattern = stratifiedMap.arr;" "caught"

printf '\nbufferedInitialization (was CalcC0UdIniBuf) -- NEVER executes\n'
run_case "inibuf/add-ncel-pmed" "DriftFluxClosure.cpp" \
    "        meanPressure = state.inletPressure;
        if (cellIndex > 0)
            upstreamMeanPressure = state.inletPressure;
        else
            upstreamMeanPressure = state.inletPressure;
        double meanTemperature = state.inletTemperature;" \
    "        meanPressure = state.inletPressure;
        if (cellIndex > 0)
            upstreamMeanPressure = state.inletPressure;
        else
            upstreamMeanPressure = state.inletPressure;
        if (cellIndex == state.lastCell)
            meanPressure = state.cells[cellIndex].pres;
        double meanTemperature = state.inletTemperature;" "caught"
run_case "inibuf/wrong-source" "DriftFluxClosure.cpp" \
    "            flowScalesOf<BufferedSource>(state, cellIndex,
                                {.liquidDensity = liquidDensity,
                                 .gasDensity = gasDensity,
                                 .liquidViscosity = liquidViscosity,
                                 .gasViscosity = gasViscosity,
                                 .noSlipLiquidHoldup = noSlipLiquidHoldup});

        int flowPattern = 1;
        double totalLength = state.cells[cellIndex].dxL + state.cells[cellIndex].dx;" \
    "            flowScalesOf<InstantaneousSource>(state, cellIndex,
                                {.liquidDensity = liquidDensity,
                                 .gasDensity = gasDensity,
                                 .liquidViscosity = liquidViscosity,
                                 .gasViscosity = gasViscosity,
                                 .noSlipLiquidHoldup = noSlipLiquidHoldup});

        int flowPattern = 1;
        double totalLength = state.cells[cellIndex].dxL + state.cells[cellIndex].dx;" "caught"

printf '\nsteadyState (was CalcC0UdPerm) -- executes in the corpus\n'
run_case "perm/fix-A3-02-ang" "DriftFluxClosure.cpp" \
    "        double inclinationAngle = (leftCellLength * state.cells[cellIndex].dutoL.teta + cellLength * state.cells[cellIndex].duto.teta) / totalLength;" \
    "        double inclinationAngle = (cellLength * state.cells[cellIndex].dutoL.teta + leftCellLength * state.cells[cellIndex].duto.teta) / totalLength;" "caught"
run_case "perm/fix-A3-03-razdx" "DriftFluxClosure.cpp" \
    "        double lengthRatio = state.cells[cellIndex].dx / (state.cells[cellIndex].dx + state.cells[cellIndex].dxL);" \
    "        double lengthRatio = state.cells[cellIndex].dxL / (state.cells[cellIndex].dx + state.cells[cellIndex].dxL);" "caught"
run_case "perm/slip-field" "DriftFluxClosure.cpp" \
    "    applyNoSlipOverride(state.cells, cellIndex, state.input.escorregaPerm, c0, ud);" \
    "    applyNoSlipOverride(state.cells, cellIndex, state.input.escorregaTran, c0, ud);" "caught"

printf '\nshared helpers -- a defect here reaches all five at once\n'
run_case "helper/blend-reassociate" "DriftFluxClosure.cpp" \
    "        c0 = ((1. - blendRatio) * pair.stratifiedC0 + blendRatio * pair.dispersedC0);" \
    "        c0 = (pair.dispersedC0 + (1. - blendRatio) * (pair.stratifiedC0 - pair.dispersedC0));" "caught"
run_case "helper/diameter-unconditional" "DriftFluxClosure.cpp" \
    "    double diameter = state.cells[cellIndex].duto.a;
    if (cellIndex > 0 && gasRate >= 0)
        diameter = state.cells[cellIndex - 1].duto.a;" \
    "    double diameter = state.cells[cellIndex].duto.a;
    if (cellIndex > 0)
        diameter = state.cells[cellIndex - 1].duto.a;" "caught"
run_case "helper/noslip-threshold" "DriftFluxClosure.cpp" \
    "        double driftCorrection = 1 - (meanSuperficialLiquidVelocity - 0.15) / 0.35;" \
    "        double driftCorrection = 1 - (meanSuperficialLiquidVelocity - 0.16) / 0.35;" "invisible"

printf '\ndeclared blind spots -- these MUST pass; the token comparison covers them\n'
# betneg feeds upstreamLiquidFlowRate, which feeds nothing but
# `if (upstreamLiquidFlowRate < 0.) mult0 = 0.', and mult0 is never read. The
# chain is dead in ALL FIVE variants, so no behavioural harness can observe a
# corruption confined to it -- including the 0.99*QG + 0.01*QL expression the
# contract requires preserved literally (A3-07).
run_case "dead-chain-0.99QG" "DriftFluxClosure.cpp" \
    "if ((0.99 * state.cells[cellIndex - 1].QG + 0.01 * state.cells[cellIndex - 1].QL) < 0.)" \
    "if ((state.cells[cellIndex - 1].QG + (state.cells[cellIndex - 1].QL - state.cells[cellIndex - 1].QG) * 0.01) < 0.)" "invisible"
run_case "dead-local-mult0" "DriftFluxClosure.cpp" \
    "    mult0 = 1.;" "    mult0 = 2.;" "invisible"
run_case "dead-local-timeStep" "DriftFluxClosure.cpp" \
    "    int timeStep = 20;
    state.cells[cellIndex].transic0 = state.cells[cellIndex].transic;" \
    "    int timeStep = 21;
    state.cells[cellIndex].transic0 = state.cells[cellIndex].transic;" "invisible"

printf '\n'
if (( fail == 0 )); then
    printf '%sCALIBRATION PASSED -- %s case(s), 0 failure(s)%s\n' "$green" "$pass" "$reset"
    exit 0
fi
printf '%sCALIBRATION FAILED -- %s of %s case(s)%s\n' "$red" "$fail" "$((pass + fail))" "$reset" >&2
exit 1
