#!/usr/bin/env bash
# Prove thermal-move.py accepts its own extraction and rejects altered tokens.
set -uo pipefail
export LC_ALL=C

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "$script_dir/.." && pwd)"
tool="$script_dir/thermal-move.py"
baseline_ref="${MARLIM_THERMAL_MOVE_BASELINE:-f8dda67}"
calctemp_baseline_ref="${MARLIM_CALCTEMP_MOVE_BASELINE:-16f7609}"
scratch="$(mktemp -d -t marlim3-thermal-move-cal-XXXXXX)"
trap 'rm -rf "$scratch"' EXIT
baseline="$scratch/SisProd-before-t061.cpp"
calctemp_baseline="$scratch/SisProd-before-t062.cpp"

git -C "$project_root" show "$baseline_ref:src/core/SisProd.cpp" > "$baseline" || {
    echo "cannot read thermal move baseline $baseline_ref" >&2
    exit 2
}
git -C "$project_root" show \
    "$calctemp_baseline_ref:src/core/SisProd.cpp" > "$calctemp_baseline" || {
    echo "cannot read calctemp move baseline $calctemp_baseline_ref" >&2
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

if (( failed > 0 )); then
    printf 'THERMAL MOVE CALIBRATION FAILED -- %d of %d case(s)\n' \
           "$failed" "$((passed + failed))" >&2
    exit 1
fi
printf 'THERMAL MOVE CALIBRATION PASSED -- %d cases\n' "$passed"
