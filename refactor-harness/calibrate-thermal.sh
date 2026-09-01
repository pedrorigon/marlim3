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

queue_case latent-blend-order SisProdThermal.cpp \
    'latt = (1 - raztemp) * latp1 + raztemp * latp2;' \
    'latt = (1 - raztemp) * latp2 + raztemp * latp1;' caught

queue_case mixture-energy-sign SisProdThermal.cpp \
    'return energintmixT0 - (delFlux' \
    'return energintmixT0 + (delFlux' caught

queue_case boundary-temperature-write SisProd.cpp \
    $'void SProd::atualizaPeriTempProd(int i) {\n    if (i > 0)\n        celula[i - 1].tempR = celula[i].temp;' \
    $'void SProd::atualizaPeriTempProd(int i) {\n    if (i > 0)\n        celula[i - 1].tempR = celula[i].temp + 1.;' caught

queue_case outlet-temperature-source SisProd.cpp \
    'tempSup = celula[ncel - 1].temp;' \
    'tempSup = celula[ncel].temp;' caught

queue_case diffusion-preparation-area SisProdThermal.cpp \
    $'void prepareNonDimensionalHeatDiffusion(const ThermalState &state, int i) {\n    double dia = state.cells[i].duto.a;\n    double area = 0.25 * M_PI * dia * dia;' \
    $'void prepareNonDimensionalHeatDiffusion(const ThermalState &state, int i) {\n    double dia = state.cells[i].duto.a;\n    double area = 0.50 * M_PI * dia * dia;' caught

queue_case tabulated-gas-density SisProdThermal.cpp \
    'double energ1 = alfmed * rhogp1 * (hgp1 - pres1 * 98066.5 / rhogp0) +' \
    'double energ1 = alfmed * rhogp1 * (hgp1 - pres1 * 98066.5 / rhogp1) +' caught

queue_case enthalpy-search-condition SisProdThermal.cpp \
    'while (j < ndiv + 1 || (energint >= val1 && energint <= val2) ||' \
    'while (j < ndiv + 1 && (energint >= val1 && energint <= val2) ||' caught

queue_case reverse-ambient-neighbor SisProd.cpp \
    'celula[i].temp = celula[i + 1].calor.Textern1;' \
    'celula[i].temp = celula[i - 1].calor.Textern1;' caught

queue_case outlet-boundary-index SisProdThermal.cpp \
    $'void updateOutletFlowPartitionTerms(const ThermalState &state) {\n\n    int i = state.lastCell;' \
    $'void updateOutletFlowPartitionTerms(const ThermalState &state) {\n\n    int i = state.lastCell - 1;' caught

queue_case inlet-boundary-propagation SisProdThermal.cpp \
    $'void updateInletFlowPartitionTerms(const ThermalState &state) {\n\n    if (state.inletMassFraction < 1) {\n        int para;\n        para = 0;\n    }\n\n    int i = 0;' \
    $'void updateInletFlowPartitionTerms(const ThermalState &state) {\n\n    if (state.inletMassFraction < 1) {\n        int para;\n        para = 0;\n    }\n\n    int i = 1;' caught

queue_case forward-hot-guard SisProd.cpp \
    'if (fabs(ugsmed + ulsmed) > 0.05 && semTermo == 0) {' \
    'if (fabs(ugsmed + ulsmed) > 1.e99 && semTermo == 0) {' caught

queue_case reverse-hot-guard SisProd.cpp \
    'if (fabs(ugsmed + ulsmed) > trocaTermicaLenta) {' \
    'if (fabs(ugsmed + ulsmed) > trocaTermicaLenta && semTermo == 0) {' caught

queue_case forward-ambient-neighbor SisProd.cpp \
    'celula[i].temp = celula[i - 1].calor.Textern1;' \
    'celula[i].temp = celula[i + 1].calor.Textern1;' caught

queue_case reverse-signed-gas SisProd.cpp \
    'ugsmed = fabs(celula[i + 1].QG) / area;' \
    'ugsmed = celula[i + 1].QG / area;' caught

queue_case reverse-interface-pressure SisProd.cpp \
    'double pmedi = celula[i + 1].presaux - celula[i].dpB / 98066.5;' \
    'double pmedi = celula[i + 1].presaux + celula[i].dpB / 98066.5;' caught

queue_case forward-network-resistance SisProd.cpp \
    'fluxcal = sinalJ * celula[i - 1].calor.transperm(celula[i - 1].resAcopRedeP);' \
    'fluxcal = sinalJ * celula[i - 1].calor.transperm(resanul);' caught

queue_case forward-velocity-cap SisProd.cpp \
    'if ((*vg1dSP).blackOilTemp == 1 && fabs(ugsmed) > 5)' \
    'if ((*vg1dSP).blackOilTemp == 2 && fabs(ugsmed) > 5)' caught

queue_case reverse-bcs-gradient SisProd.cpp \
    'dpdx = (pmedi - celula[i + 1].pres) * 98600. / dx;' \
    'dpdx = (pmedi - celula[i + 1].pres) * 98066.5 / dx;' caught

queue_case forward-latent-limit SisProd.cpp \
    'if (arq.limTransMass < valTransMass)' \
    'if (arq.limTransMass < valTransMass * 0.)' caught

queue_case reverse-latent-sign SisProd.cpp \
    'latente = -interpolaHLatente(pmed, tmed) * celula[i + 1].FonteMudaFase;' \
    'latente = interpolaHLatente(pmed, tmed) * celula[i + 1].FonteMudaFase;' caught

queue_case forward-latent-switch SisProd.cpp \
    'if (arq.latente == 0)' \
    'if (arq.latente == -1)' caught

queue_case reverse-annulus-gradient SisProd.cpp \
    'double dtext = (celulaG[k].temp - celulaG[k - 1].temp) / npasso;' \
    'double dtext = (celulaG[k - 1].temp - celulaG[k].temp) / npasso;' caught

queue_case forward-potential-sign SisProd.cpp \
    'double hidro = (rhol * ulsmed + rhog * ugsmed) * area * 9.82 * sin(celula[i - 1].duto.teta);' \
    'double hidro = -(rhol * ulsmed + rhog * ugsmed) * area * 9.82 * sin(celula[i - 1].duto.teta);' caught

run_queued_cases

if (( failed > 0 )); then
    printf '%sTHERMAL CALIBRATION FAILED -- %d of %d case(s)%s\n' \
           "$red" "$failed" "$((passed + failed))" "$reset" >&2
    exit 1
fi

printf '%sTHERMAL CALIBRATION PASSED -- %d cases%s\n' \
       "$green" "$passed" "$reset"
