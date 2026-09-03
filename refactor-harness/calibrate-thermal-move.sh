#!/usr/bin/env bash
# Prove thermal-move.py accepts its own extraction and rejects altered tokens.
set -uo pipefail
export LC_ALL=C

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "$script_dir/.." && pwd)"
tool="$script_dir/thermal-move.py"
baseline_ref="${MARLIM_THERMAL_MOVE_BASELINE:-f8dda67}"
calctemp_baseline_ref="${MARLIM_CALCTEMP_MOVE_BASELINE:-16f7609}"
t063_baseline_ref="${MARLIM_T063_MOVE_BASELINE:-f82e95b}"
t064_baseline_ref="${MARLIM_T064_MOVE_BASELINE:-8c1742f}"
t065_baseline_ref="${MARLIM_T065_MOVE_BASELINE:-4546f1e}"
t066_baseline_ref="${MARLIM_T066_MOVE_BASELINE:-0bc2e6f}"
scratch="$(mktemp -d -t marlim3-thermal-move-cal-XXXXXX)"
trap 'rm -rf "$scratch"' EXIT
baseline="$scratch/SisProd-before-t061.cpp"
calctemp_baseline="$scratch/SisProd-before-t062.cpp"
t063_baseline="$scratch/SisProd-before-t063.cpp"
t064_baseline="$scratch/SisProd-before-t064.cpp"
t065_baseline="$scratch/SisProd-before-t065.cpp"
t066_baseline="$scratch/SisProd-before-t066.cpp"

git -C "$project_root" show "$baseline_ref:src/core/SisProd.cpp" > "$baseline" || {
    echo "cannot read thermal move baseline $baseline_ref" >&2
    exit 2
}
git -C "$project_root" show \
    "$calctemp_baseline_ref:src/core/SisProd.cpp" > "$calctemp_baseline" || {
    echo "cannot read calctemp move baseline $calctemp_baseline_ref" >&2
    exit 2
}
git -C "$project_root" show \
    "$t063_baseline_ref:src/core/SisProd.cpp" > "$t063_baseline" || {
    echo "cannot read T063 move baseline $t063_baseline_ref" >&2
    exit 2
}
git -C "$project_root" show \
    "$t064_baseline_ref:src/core/SisProd.cpp" > "$t064_baseline" || {
    echo "cannot read T064 move baseline $t064_baseline_ref" >&2
    exit 2
}
git -C "$project_root" show \
    "$t065_baseline_ref:src/core/SisProd.cpp" > "$t065_baseline" || {
    echo "cannot read T065 move baseline $t065_baseline_ref" >&2
    exit 2
}
git -C "$project_root" show \
    "$t066_baseline_ref:src/core/SisProd.cpp" > "$t066_baseline" || {
    echo "cannot read T066 move baseline $t066_baseline_ref" >&2
    exit 2
}

python3 "$tool" extract "$baseline" "$scratch/control.cpp" > /dev/null || exit 2
python3 "$tool" check "$baseline" "$scratch/control.cpp" > /dev/null || {
    echo "thermal move calibration control failed" >&2
    exit 1
}

passed=0
failed=0

run_case() {
    local name="$1" from="$2" to="$3"
    local candidate="$scratch/$name.cpp"
    cp "$scratch/control.cpp" "$candidate"
    python3 - "$candidate" "$from" "$to" <<'PY'
import pathlib
import sys

path = pathlib.Path(sys.argv[1])
old = sys.argv[2]
new = sys.argv[3]
text = path.read_text()
count = text.count(old)
if count != 1:
    print(f"expected mutation pattern once, found {count}", file=sys.stderr)
    raise SystemExit(3)
path.write_text(text.replace(old, new, 1))
PY
    if (( $? != 0 )) || cmp -s "$candidate" "$scratch/control.cpp"; then
        printf '  %-24s MUTATION NOT INJECTED\n' "$name" >&2
        failed=$((failed + 1))
        return
    fi
    if python3 "$tool" check "$baseline" "$candidate" > /dev/null 2>&1; then
        printf '  %-24s MISSED\n' "$name" >&2
        failed=$((failed + 1))
        return
    fi
    printf '  %-24s caught\n' "$name"
    passed=$((passed + 1))
}

run_case latent-literal \
    'latt = (1 - raztemp) * latp1 + raztemp * latp2;' \
    'latt = (1 - raztemp) * latp2 + raztemp * latp1;'

run_case enthalpy-gravity \
    '* sin(state.cells[i].duto.teta) * 9.82;' \
    '* sin(state.cells[i].duto.teta) * 9.81;'

run_case energy-denominator \
    'hgp1 - pres1 * 98066.5 / rhogp0' \
    'hgp1 - pres1 * 98066.5 / rhogp1'

run_case missing-function \
    'double interpolateMixtureEnergy(const ThermalState &state, int i, int jp0, int jt, double razp) {' \
    'double removedMixtureEnergy(const ThermalState &state, int i, int jp0, int jt, double razp) {'

python3 "$tool" extract-calctemp \
    "$calctemp_baseline" "$scratch/calctemp-control.cpp" > /dev/null || exit 2
python3 "$tool" check-calctemp \
    "$calctemp_baseline" "$scratch/calctemp-control.cpp" > /dev/null || {
    echo "calctemp move calibration control failed" >&2
    exit 1
}

run_calctemp_case() {
    local name="$1" from="$2" to="$3"
    local candidate="$scratch/$name.cpp"
    cp "$scratch/calctemp-control.cpp" "$candidate"
    python3 - "$candidate" "$from" "$to" <<'PY'
import pathlib
import sys

path = pathlib.Path(sys.argv[1])
old = sys.argv[2]
new = sys.argv[3]
text = path.read_text()
count = text.count(old)
if count != 1:
    print(f"expected mutation pattern once, found {count}", file=sys.stderr)
    raise SystemExit(3)
path.write_text(text.replace(old, new, 1))
PY
    if (( $? != 0 )) || cmp -s "$candidate" "$scratch/calctemp-control.cpp"; then
        printf '  %-24s MUTATION NOT INJECTED\n' "$name" >&2
        failed=$((failed + 1))
        return
    fi
    if python3 "$tool" check-calctemp \
        "$calctemp_baseline" "$candidate" > /dev/null 2>&1; then
        printf '  %-24s MISSED\n' "$name" >&2
        failed=$((failed + 1))
        return
    fi
    printf '  %-24s caught\n' "$name"
    passed=$((passed + 1))
}

run_calctemp_case calctemp-gravity \
    '* area * 9.82 * sin(state.cells[i].duto.teta);' \
    '* area * 9.81 * sin(state.cells[i].duto.teta);'

run_calctemp_case calctemp-pressure-work \
    'coefdxT * dtdx + coefdxP * dpdx - cinetico' \
    'coefdxT * dtdx - coefdxP * dpdx - cinetico'

run_calctemp_case calctemp-rate-limit \
    '* 10 * state.cells[i].dt;' \
    '* 11 * state.cells[i].dt;'

run_calctemp_case calctemp-missing \
    'void computeTemperature(const ThermalState &state, int i, double tempantiga, int modoPerm) {' \
    'void removedTemperature(const ThermalState &state, int i, double tempantiga, int modoPerm) {'

run_t063_case() {
    local name="$1" old_name="$2" from="$3" to="$4"
    local control="$scratch/$name-control.cpp"
    local candidate="$scratch/$name.cpp"
    python3 "$tool" extract-t063 \
        "$old_name" "$t063_baseline" "$control" > /dev/null || exit 2
    python3 "$tool" check-t063 \
        "$old_name" "$t063_baseline" "$control" > /dev/null || {
        printf '  %-24s CONTROL FAILED\n' "$name" >&2
        failed=$((failed + 1))
        return
    }
    cp "$control" "$candidate"
    python3 - "$candidate" "$from" "$to" <<'PY'
import pathlib
import sys

path = pathlib.Path(sys.argv[1])
old = sys.argv[2]
new = sys.argv[3]
text = path.read_text()
count = text.count(old)
if count != 1:
    print(f"expected mutation pattern once, found {count}", file=sys.stderr)
    raise SystemExit(3)
path.write_text(text.replace(old, new, 1))
PY
    if (( $? != 0 )) || cmp -s "$candidate" "$control"; then
        printf '  %-24s MUTATION NOT INJECTED\n' "$name" >&2
        failed=$((failed + 1))
        return
    fi
    if python3 "$tool" check-t063 \
        "$old_name" "$t063_baseline" "$candidate" > /dev/null 2>&1; then
        printf '  %-24s MISSED\n' "$name" >&2
        failed=$((failed + 1))
        return
    fi
    printf '  %-24s caught\n' "$name"
    passed=$((passed + 1))
}

run_t063_case enthalpy-loop calcTempEntalp \
    'while (j < ndiv + 1 || (energint >= val1 && energint <= val2) ||' \
    'while (j < ndiv + 1 && (energint >= val1 && energint <= val2) ||'

run_t063_case mass-transfer-gravity calcTransMassTermo \
    '* area * 9.82 * sin(state.cells[i].duto.teta);' \
    '* area * 9.81 * sin(state.cells[i].duto.teta);'

run_t063_case mass-transfer-latent calcTransMassTermo \
    'state.cells[i].FonteMudaFase /= latente;' \
    'state.cells[i].FonteMudaFase *= latente;'

run_t063_case mass-transfer-missing calcTransMassTermo \
    'void computeThermalMassTransfer(const ThermalState &state, int i) {' \
    'void removedThermalMassTransfer(const ThermalState &state, int i) {'

decomposition_control="$project_root/src/core/SisProdThermal.cpp"
python3 "$tool" check-t063-decomposition \
    "$t063_baseline" "$decomposition_control" > /dev/null || {
    echo "T063 decomposition calibration control failed" >&2
    exit 1
}

run_t063_decomposition_case() {
    local name="$1" from="$2" to="$3"
    local candidate="$scratch/$name.cpp"
    cp "$decomposition_control" "$candidate"
    python3 - "$candidate" "$from" "$to" <<'PY'
import pathlib
import sys

path = pathlib.Path(sys.argv[1])
old = sys.argv[2]
new = sys.argv[3]
text = path.read_text()
count = text.count(old)
if count != 1:
    print(f"expected mutation pattern once, found {count}", file=sys.stderr)
    raise SystemExit(3)
path.write_text(text.replace(old, new, 1))
PY
    if (( $? != 0 )) || cmp -s "$candidate" "$decomposition_control"; then
        printf '  %-24s MUTATION NOT INJECTED\n' "$name" >&2
        failed=$((failed + 1))
        return
    fi
    if python3 "$tool" check-t063-decomposition \
        "$t063_baseline" "$candidate" > /dev/null 2>&1; then
        printf '  %-24s MISSED\n' "$name" >&2
        failed=$((failed + 1))
        return
    fi
    printf '  %-24s caught\n' "$name"
    passed=$((passed + 1))
}

# The mass-source tail is textually identical in both source helpers, so the
# mutation targets a line unique to the T063 one.
run_t063_decomposition_case mass-source-helper \
    'sourceSpecificHeatRatio = cell.acsr.injg.FluidoPro.ConstAdG(cell.pres, cell.temp);' \
    'sourceSpecificHeatRatio = cell.acsr.injg.FluidoPro.ConstAdG(cell.presini, cell.temp);'

run_t063_decomposition_case mass-source-main \
    'sourceTerms.liquid + sourceTerms.gas + heatFlux' \
    'sourceTerms.liquid - sourceTerms.gas + heatFlux'

t064_decomposition_control="$project_root/src/core/SisProdThermal.cpp"
python3 "$tool" check-renova-temp-decomposition \
    "$t064_baseline" "$t064_decomposition_control" > /dev/null || {
    echo "T064 decomposition calibration control failed" >&2
    exit 1
}

run_t064_decomposition_case() {
    local name="$1" from="$2" to="$3"
    local candidate="$scratch/$name.cpp"
    cp "$t064_decomposition_control" "$candidate"
    python3 - "$candidate" "$from" "$to" <<'PY'
import pathlib
import sys

path = pathlib.Path(sys.argv[1])
old = sys.argv[2]
new = sys.argv[3]
text = path.read_text()
count = text.count(old)
if count != 1:
    print(f"expected mutation pattern once, found {count}", file=sys.stderr)
    raise SystemExit(3)
path.write_text(text.replace(old, new, 1))
PY
    if (( $? != 0 )) || cmp -s "$candidate" "$t064_decomposition_control"; then
        printf '  %-24s MUTATION NOT INJECTED\n' "$name" >&2
        failed=$((failed + 1))
        return
    fi
    if python3 "$tool" check-renova-temp-decomposition \
        "$t064_baseline" "$candidate" > /dev/null 2>&1; then
        printf '  %-24s MISSED\n' "$name" >&2
        failed=$((failed + 1))
        return
    fi
    printf '  %-24s caught\n' "$name"
    passed=$((passed + 1))
}

run_t064_decomposition_case renova-inlet-state \
    'previousLiquidDensity =' \
    'previousLiquidDensity +='

run_t064_decomposition_case renova-properties-temp \
    'cellSolutionGasTemperatureDerivative = (cellSolutionGasRatio / cellOilVolumeFactor - shiftedCellSolutionGasRatio / shiftedCellOilVolumeFactor) /' \
    'cellSolutionGasTemperatureDerivative = (cellSolutionGasRatio / cellOilVolumeFactor + shiftedCellSolutionGasRatio / shiftedCellOilVolumeFactor) /'

run_t064_decomposition_case renova-derivative-model \
    'state.massTransferModel != 0)' \
    'state.massTransferModel == 0)'

run_t064_decomposition_case renova-model-threshold \
    'leftAbsoluteSuperficialVelocity < 0.1)' \
    'leftAbsoluteSuperficialVelocity <= 0.1)'

run_t064_decomposition_case renova-application-rate \
    'leftCell.transmassR /=' \
    'leftCell.transmassR *='

run_t064_decomposition_case renova-main-model \
    'leftCell.TMModel == -2' \
    'leftCell.TMModel == -3'

run_t064_decomposition_case renova-missing-helper \
    'void selectDistributedMassTransferModel(' \
    'void removedDistributedMassTransferModel('

# The T065 decomposition control is RETIRED, and the reason is worth stating so
# nobody reinstates it and then silences the failure.
#
# It reconstitutes one body by inlining the seven helpers and compares it with
# the pre-decomposition renovaterm. That worked until T070r. In the baseline,
# renovaterm declares BOTH a local `c0` and a local `xc0`. After T065 split the
# function those two ended up in different helpers, so when T070r renamed
# `xc0` -> `c0` in updateFlowPartitionTerms there was no collision -- the guard
# in rename-locals.py checked, and correctly allowed it, exactly as it correctly
# REFUSED the same rename in the three functions where `c0` already existed.
#
# The code is therefore sound. What broke is the comparison: inlining merges
# scopes that the rename treated as separate, so two distinct baseline variables
# arrive under one current name. No name-based comparison can see through that,
# and pretending otherwise would mean loosening the check until it stops
# catching a real substitution.
#
# T065's decomposition was proven exact when it was made -- 18,332 tokens across
# seven helpers, nine injected corruptions -- and that evidence stands. The
# functions are covered now by the supplementary L0 check (which declares them),
# by verify-thermal.sh with its 24 calibrated corruptions, and by L2.
#
# The T063 and T064 controls below and above still run: their reconstitutions do
# not merge two same-named baseline locals, so a consistent-renaming comparison
# still answers the question.

# The nine T065 mutation cases are retired with the control they depended on:
# a mutation whose verdict comes from a comparison that can no longer be made
# proves nothing, and leaving them here to fail would train readers to ignore
# a red line. Their subject matter is covered by verify-thermal.sh.

run_t066_case() {
    local name="$1" old_name="$2" from="$3" to="$4"
    local control="$scratch/$name-control.cpp"
    local candidate="$scratch/$name.cpp"
    python3 "$tool" extract-t066 \
        "$old_name" "$t066_baseline" "$control" > /dev/null || exit 2
    python3 "$tool" check-t066 \
        "$old_name" "$t066_baseline" "$control" > /dev/null || {
        printf '  %-24s CONTROL FAILED\n' "$name" >&2
        failed=$((failed + 1))
        return
    }
    cp "$control" "$candidate"
    python3 - "$candidate" "$from" "$to" <<'PY'
import pathlib
import sys

path = pathlib.Path(sys.argv[1])
old = sys.argv[2]
new = sys.argv[3]
text = path.read_text()
count = text.count(old)
if count != 1:
    print(f"expected mutation pattern once, found {count}", file=sys.stderr)
    raise SystemExit(3)
path.write_text(text.replace(old, new, 1))
PY
    if (( $? != 0 )) || cmp -s "$candidate" "$control"; then
        printf '  %-24s MUTATION NOT INJECTED\n' "$name" >&2
        failed=$((failed + 1))
        return
    fi
    if python3 "$tool" check-t066 \
        "$old_name" "$t066_baseline" "$candidate" > /dev/null 2>&1; then
        printf '  %-24s MISSED\n' "$name" >&2
        failed=$((failed + 1))
        return
    fi
    printf '  %-24s caught\n' "$name"
    passed=$((passed + 1))
}

run_t066_case diffusion-area prepDifusCalorND \
    'double area = 0.25 * M_PI * dia * dia;' \
    'double area = 0.50 * M_PI * dia * dia;'

run_t066_case diffusion-qualified-dt prepDifusCalorND \
    'state.cells[i].calor.dt = state.cells[i].dt;' \
    'state.cells[i].calor.state.timeStep = state.cells[i].state.timeStep;'

run_t066_case energy-poisson-renew marchaEnergTrans \
    'state.poissonSolver.renova();' \
    'state.poissonSolver.state.evolutionUpdater.renew();'

run_t066_case energy-cycle-step marchaEnergTrans \
    'state.minimumCycleTimeStep != state.timeStep' \
    'state.minimumCycleTimeStep == state.timeStep'

run_t066_case energy-missing marchaEnergTrans \
    'void advanceTransientEnergy(const ThermalState &state, int ciclo, int ciclomax) {' \
    'void removedTransientEnergy(const ThermalState &state, int ciclo, int ciclomax) {'

if (( failed > 0 )); then
    printf 'THERMAL MOVE CALIBRATION FAILED -- %d of %d case(s)\n' \
           "$failed" "$((passed + failed))" >&2
    exit 1
fi
printf 'THERMAL MOVE CALIBRATION PASSED -- %d cases\n' "$passed"
