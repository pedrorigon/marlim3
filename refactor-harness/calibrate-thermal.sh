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
            return
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
        else
            printf '  %-28s %sFAILED CONTROL%s\n%s\n' \
                   "$name" "$red" "$reset" "$output"
            failed=$((failed + 1))
        fi
    elif (( exit_code == 1 )); then
        printf '  %-28s %scaught%s\n' "$name" "$green" "$reset"
        passed=$((passed + 1))
    else
        printf '  %-28s %sMISSED OR INVALID%s exit=%d\n%s\n' \
               "$name" "$red" "$reset" "$exit_code" "$output"
        failed=$((failed + 1))
    fi
}

run_case control SisProd.cpp "" "" pass
if (( failed > 0 )); then
    echo "thermal calibration control failed" >&2
    exit 1
fi

run_case latent-blend-order SisProdThermal.cpp \
    'latt = (1 - raztemp) * latp1 + raztemp * latp2;' \
    'latt = (1 - raztemp) * latp2 + raztemp * latp1;' caught

run_case mixture-energy-sign SisProdThermal.cpp \
    'return energintmixT0 - (delFlux' \
    'return energintmixT0 + (delFlux' caught

run_case boundary-temperature-write SisProd.cpp \
    $'void SProd::atualizaPeriTempProd(int i) {\n    if (i > 0)\n        celula[i - 1].tempR = celula[i].temp;' \
    $'void SProd::atualizaPeriTempProd(int i) {\n    if (i > 0)\n        celula[i - 1].tempR = celula[i].temp + 1.;' caught

run_case outlet-temperature-source SisProd.cpp \
    'tempSup = celula[ncel - 1].temp;' \
    'tempSup = celula[ncel].temp;' caught

run_case diffusion-preparation-area SisProd.cpp \
    $'void SProd::prepDifusCalorND(int i) {\n    double dia = celula[i].duto.a;\n    double area = 0.25 * M_PI * dia * dia;' \
    $'void SProd::prepDifusCalorND(int i) {\n    double dia = celula[i].duto.a;\n    double area = 0.50 * M_PI * dia * dia;' caught

run_case tabulated-gas-density SisProdThermal.cpp \
    'double energ1 = alfmed * rhogp1 * (hgp1 - pres1 * 98066.5 / rhogp0) +' \
    'double energ1 = alfmed * rhogp1 * (hgp1 - pres1 * 98066.5 / rhogp1) +' caught

run_case enthalpy-search-condition SisProd.cpp \
    'while (j < ndiv + 1 || (energint >= val1 && energint <= val2) ||' \
    'while (j < ndiv + 1 && (energint >= val1 && energint <= val2) ||' caught

run_case reverse-ambient-neighbor SisProd.cpp \
    'celula[i].temp = celula[i + 1].calor.Textern1;' \
    'celula[i].temp = celula[i - 1].calor.Textern1;' caught

run_case outlet-boundary-index SisProd.cpp \
    $'void SProd::renovatermAfluFim() {\n\n    int i = ncel;' \
    $'void SProd::renovatermAfluFim() {\n\n    int i = ncel - 1;' caught

run_case inlet-boundary-propagation SisProd.cpp \
    $'    int i = 0;\n\n    double razdx = 0.5;' \
    $'    int i = 1;\n\n    double razdx = 0.5;' caught

if (( failed > 0 )); then
    printf '%sTHERMAL CALIBRATION FAILED -- %d of %d case(s)%s\n' \
           "$red" "$failed" "$((passed + failed))" "$reset" >&2
    exit 1
fi

printf '%sTHERMAL CALIBRATION PASSED -- %d cases%s\n' \
       "$green" "$passed" "$reset"
