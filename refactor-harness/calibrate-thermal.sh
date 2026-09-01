#!/usr/bin/env bash
# Calibration for verify-thermal.sh. Every corruption must be injected exactly
# once, and untouched source must agree with the golden before failures count.
set -uo pipefail
export LC_ALL=C

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "$script_dir/.." && pwd)"
golden="${1:-$project_root/specs/001-refatoracao-sisprod/golden/thermal-baseline.txt}"

red=$'\033[0;31m'; green=$'\033[0;32m'; reset=$'\033[0m'

[[ -s "$golden" ]] || {
    printf '%sno thermal golden at %s%s\n' "$red" "$golden" "$reset" >&2
    exit 2
}

scratch="$(mktemp -d -t marlim3-thermal-cal-XXXXXX)"
trap 'rm -rf "$scratch"' EXIT
pristine="$project_root/src/core/SisProd.cpp"
pristine_thermal="$project_root/src/core/SisProdThermal.cpp"
passed=0
failed=0

run_case() {
    local name="$1" target="$2" from="$3" to="$4" expectation="$5"
    local work="$scratch/$name"
    mkdir -p "$work"
    cp "$pristine" "$work/SisProd.cpp"
    cp "$pristine_thermal" "$work/SisProdThermal.cpp"

    if [[ -n "$from" ]]; then
        python3 - "$work/$target" "$from" "$to" <<'PY'
import pathlib
import sys

path = pathlib.Path(sys.argv[1])
old = sys.argv[2]
new = sys.argv[3]
text = path.read_text()
count = text.count(old)
if count != 1:
    print(f"expected corruption pattern exactly once, found {count}", file=sys.stderr)
    raise SystemExit(3)
path.write_text(text.replace(old, new, 1))
PY
        if (( $? != 0 )) || cmp -s "$project_root/src/core/$target" "$work/$target"; then
            printf '  %-28s %sCORRUPTION NOT INJECTED%s\n' "$name" "$red" "$reset"
            failed=$((failed + 1))
            return 1
        fi
    fi

    local output exit_code
    output=$(MARLIM_THERMAL_SOURCES="$work/SisProd.cpp $work/SisProdThermal.cpp" \
        bash "$script_dir/verify-thermal.sh" compare "$golden" 2>&1)
    exit_code=$?

    if [[ "$expectation" == "pass" ]]; then
        if (( exit_code == 0 )); then
            printf '  %-28s %spassed%s\n' "$name" "$green" "$reset"
            passed=$((passed + 1))
            return 0
        else
            printf '  %-28s %sFAILED CONTROL%s\n%s\n' \
                   "$name" "$red" "$reset" "$output"
            failed=$((failed + 1))
            return 1
        fi
    elif (( exit_code == 1 )); then
        printf '  %-28s %scaught%s\n' "$name" "$green" "$reset"
        passed=$((passed + 1))
        return 0
    else
        printf '  %-28s %sMISSED OR INVALID%s exit=%d\n%s\n' \
               "$name" "$red" "$reset" "$exit_code" "$output"
        failed=$((failed + 1))
        return 1
    fi
}

declare -a queued_names queued_targets queued_from queued_to queued_expectations

queue_case() {
    queued_names+=("$1")
    queued_targets+=("$2")
    queued_from+=("$3")
    queued_to+=("$4")
    queued_expectations+=("$5")
}

run_queued_cases() {
    local jobs="${MARLIM_CAL_JOBS:-4}"
    [[ "$jobs" =~ ^[1-9][0-9]*$ ]] || {
        echo "MARLIM_CAL_JOBS must be a positive integer" >&2
        return 2
    }

    local index=0 total="${#queued_names[@]}"
    while (( index < total )); do
        local -a pids=()
        local batch=0
        while (( index < total && batch < jobs )); do
            run_case "${queued_names[index]}" "${queued_targets[index]}" \
                "${queued_from[index]}" "${queued_to[index]}" \
                "${queued_expectations[index]}" &
            pids+=("$!")
            index=$((index + 1))
            batch=$((batch + 1))
        done
        local pid
        for pid in "${pids[@]}"; do
            if wait "$pid"; then
                passed=$((passed + 1))
            else
                failed=$((failed + 1))
            fi
        done
    done
}

run_case control SisProd.cpp "" "" pass
if (( failed > 0 )); then
    echo "thermal calibration control failed" >&2
    exit 1
fi

# All bodies below now live in SisProdThermal.cpp with SI English locals: the
# stage moved them (T068, T069) and renamed them (T070r). The previous patterns
# targeted SisProd.cpp with the Portuguese names and stopped matching, which the
# runner reported honestly as CORRUPTION NOT INJECTED rather than as a pass.
queue_case latent-blend-order SisProdThermal.cpp \
    'latentHeat = (1 - temperatureRatio) * latentHeatAtTemperatureIndex + temperatureRatio * latentHeatAtNextTemperatureIndex;' \
    'latentHeat = (1 - temperatureRatio) * latentHeatAtNextTemperatureIndex + temperatureRatio * latentHeatAtTemperatureIndex;' caught

queue_case mixture-energy-sign SisProdThermal.cpp \
    'return previousMixtureInternalEnergy - (enthalpyFluxDivergence' \
    'return previousMixtureInternalEnergy + (enthalpyFluxDivergence' caught

queue_case boundary-temperature-write SisProdThermal.cpp \
    $'void updateProductionTemperaturePeriphery(const ThermalState &state, int cellIndex) {\n    if (cellIndex > 0)\n        state.cells[cellIndex - 1].tempR = state.cells[cellIndex].temp;' \
    $'void updateProductionTemperaturePeriphery(const ThermalState &state, int cellIndex) {\n    if (cellIndex > 0)\n        state.cells[cellIndex - 1].tempR = state.cells[cellIndex].temp + 1.;' caught

queue_case outlet-temperature-source SisProdThermal.cpp \
    'state.surfaceTemperature = state.cells[state.lastCell - 1].temp;' \
    'state.surfaceTemperature = state.cells[state.lastCell].temp;' caught

queue_case diffusion-preparation-area SisProdThermal.cpp \
    $'void prepareNonDimensionalHeatDiffusion(const ThermalState &state, int cellIndex) {\n    double diameter = state.cells[cellIndex].duto.a;\n    double flowArea = 0.25 * M_PI * diameter * diameter;' \
    $'void prepareNonDimensionalHeatDiffusion(const ThermalState &state, int cellIndex) {\n    double diameter = state.cells[cellIndex].duto.a;\n    double flowArea = 0.50 * M_PI * diameter * diameter;' caught

queue_case tabulated-gas-density SisProdThermal.cpp \
    'double upperMixtureEnergy = meanVoidFraction * upperGasDensity * (upperGasEnthalpy - upperPressure * 98066.5 / lowerGasDensity) +' \
    'double upperMixtureEnergy = meanVoidFraction * upperGasDensity * (upperGasEnthalpy - upperPressure * 98066.5 / upperGasDensity) +' caught

queue_case enthalpy-search-condition SisProdThermal.cpp \
    'while (temperatureIndex < divisionCount + 1 || (mixtureInternalEnergy >= lowerEnergy && mixtureInternalEnergy <= upperEnergy) ||' \
    'while (temperatureIndex < divisionCount + 1 && (mixtureInternalEnergy >= lowerEnergy && mixtureInternalEnergy <= upperEnergy) ||' caught

queue_case reverse-ambient-neighbor SisProdThermal.cpp \
    'state.cells[cellIndex].temp = state.cells[cellIndex + 1].calor.Textern1;' \
    'state.cells[cellIndex].temp = state.cells[cellIndex - 1].calor.Textern1;' caught

queue_case outlet-boundary-index SisProdThermal.cpp \
    $'void updateOutletFlowPartitionTerms(const ThermalState &state) {\n\n    int cellIndex = state.lastCell;' \
    $'void updateOutletFlowPartitionTerms(const ThermalState &state) {\n\n    int cellIndex = state.lastCell - 1;' caught

queue_case inlet-boundary-propagation SisProdThermal.cpp \
    $'void updateInletFlowPartitionTerms(const ThermalState &state) {\n\n    if (state.inletMassFraction < 1) {\n        int debugStop;\n        debugStop = 0;\n    }\n\n    int cellIndex = 0;' \
    $'void updateInletFlowPartitionTerms(const ThermalState &state) {\n\n    if (state.inletMassFraction < 1) {\n        int debugStop;\n        debugStop = 0;\n    }\n\n    int cellIndex = 1;' caught

queue_case forward-hot-guard SisProdThermal.cpp \
    'if (fabs(meanSuperficialGasVelocity + meanSuperficialLiquidVelocity) > 0.05 && state.thermalSourceDisabled == 0) {' \
    'if (fabs(meanSuperficialGasVelocity + meanSuperficialLiquidVelocity) > 1.e99 && state.thermalSourceDisabled == 0) {' caught

queue_case reverse-hot-guard SisProdThermal.cpp \
    'if (fabs(meanSuperficialGasVelocity + meanSuperficialLiquidVelocity) > state.slowHeatTransferThreshold) {' \
    'if (fabs(meanSuperficialGasVelocity + meanSuperficialLiquidVelocity) > state.slowHeatTransferThreshold && state.thermalSourceDisabled == 0) {' caught

queue_case forward-ambient-neighbor SisProdThermal.cpp \
    'state.cells[cellIndex].temp = state.cells[cellIndex - 1].calor.Textern1;' \
    'state.cells[cellIndex].temp = state.cells[cellIndex + 1].calor.Textern1;' caught

queue_case reverse-signed-gas SisProdThermal.cpp \
    'meanSuperficialGasVelocity = fabs(state.cells[cellIndex + 1].QG) / flowArea;' \
    'meanSuperficialGasVelocity = state.cells[cellIndex + 1].QG / flowArea;' caught

queue_case reverse-interface-pressure SisProdThermal.cpp \
    'double interfaceMeanPressure = state.cells[cellIndex + 1].presaux - state.cells[cellIndex].dpB / 98066.5;' \
    'double interfaceMeanPressure = state.cells[cellIndex + 1].presaux + state.cells[cellIndex].dpB / 98066.5;' caught

queue_case forward-network-resistance SisProdThermal.cpp \
    'heatFlux = mixtureFluxSign * state.cells[cellIndex - 1].calor.transperm(state.cells[cellIndex - 1].resAcopRedeP);' \
    'heatFlux = mixtureFluxSign * state.cells[cellIndex - 1].calor.transperm(annulusResistance);' caught

queue_case forward-velocity-cap SisProdThermal.cpp \
    'if ((*state.globals).blackOilTemp == 1 && fabs(meanSuperficialGasVelocity) > 5)' \
    'if ((*state.globals).blackOilTemp == 2 && fabs(meanSuperficialGasVelocity) > 5)' caught

queue_case reverse-bcs-gradient SisProdThermal.cpp \
    'pressureGradient = (interfaceMeanPressure - state.cells[cellIndex + 1].pres) * 98600. / cellLength;' \
    'pressureGradient = (interfaceMeanPressure - state.cells[cellIndex + 1].pres) * 98066.5 / cellLength;' caught

# Two bodies carry this guard; the cellIndex - 1 read pins it to the forward
# steady march rather than to computeTemperature.
queue_case forward-latent-limit SisProdThermal.cpp \
    $'phaseChangeSign = state.cells[cellIndex - 1].FonteMudaFase / phaseChangeMassRate;\n    if (state.input.limTransMass < phaseChangeMassRate)' \
    $'phaseChangeSign = state.cells[cellIndex - 1].FonteMudaFase / phaseChangeMassRate;\n    if (state.input.limTransMass < phaseChangeMassRate * 0.)' caught

queue_case reverse-latent-sign SisProdThermal.cpp \
    'latentHeatTerm = -interpolateLatentHeat(state, meanPressure, meanTemperature) * state.cells[cellIndex + 1].FonteMudaFase;' \
    'latentHeatTerm = interpolateLatentHeat(state, meanPressure, meanTemperature) * state.cells[cellIndex + 1].FonteMudaFase;' caught

# Same guard appears in computeTemperature; the preceding comment pins this one
# to the forward steady march.
queue_case forward-latent-switch SisProdThermal.cpp \
    $'// mass sources, latent heat, and shaft work.\n            if (state.input.latente == 0)' \
    $'// mass sources, latent heat, and shaft work.\n            if (state.input.latente == -1)' caught

queue_case reverse-annulus-gradient SisProdThermal.cpp \
    'double externalTemperatureStep = (state.gasCells[k].temp - state.gasCells[k - 1].temp) / subStepCount;' \
    'double externalTemperatureStep = (state.gasCells[k - 1].temp - state.gasCells[k].temp) / subStepCount;' caught

queue_case forward-potential-sign SisProdThermal.cpp \
    'double hydrostaticPower = (liquidDensity * meanSuperficialLiquidVelocity + gasDensity * meanSuperficialGasVelocity) * flowArea * 9.82 * sin(state.cells[cellIndex - 1].duto.teta);' \
    'double hydrostaticPower = -(liquidDensity * meanSuperficialLiquidVelocity + gasDensity * meanSuperficialGasVelocity) * flowArea * 9.82 * sin(state.cells[cellIndex - 1].duto.teta);' caught


run_queued_cases

if (( failed > 0 )); then
    printf '%sTHERMAL CALIBRATION FAILED -- %d of %d case(s)%s\n' \
           "$red" "$failed" "$((passed + failed))" "$reset" >&2
    exit 1
fi

printf '%sTHERMAL CALIBRATION PASSED -- %d cases%s\n' \
       "$green" "$passed" "$reset"
