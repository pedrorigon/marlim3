#!/usr/bin/env bash
# Calibrates verify-gaslift.sh: a harness that never fails is not a harness.
#
# Each corruption is injected into ONE routine that the sweep drives, and the run
# CONFIRMS the injection landed before asking for a verdict -- a corruption that
# fails to apply is reported as SKIP, never counted as a success. That
# distinction is what kept the thermal harness honest when the stage-5 rename
# left 23 of its 24 cases unable to inject.
set -uo pipefail
export LC_ALL=C

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "$script_dir/.." && pwd)"
target="$project_root/src/core/SisProd.cpp"
golden="$project_root/specs/001-refatoracao-sisprod/golden/gaslift-baseline.txt"

red=$'\033[0;31m'; green=$'\033[0;32m'; yellow=$'\033[1;33m'; reset=$'\033[0m'
scratch="$(mktemp -d -t marlim3-calgl-XXXXXX)"
trap 'cp "$scratch/pristine.cpp" "$target"; rm -rf "$scratch"' EXIT
cp "$target" "$scratch/pristine.cpp"

detected=0; missed=0; skipped=0
probe() { # <label> <from> <to>
    cp "$scratch/pristine.cpp" "$target"
    python3 - "$target" "$2" "$3" <<'PY'
import pathlib, sys
path = pathlib.Path(sys.argv[1]); old, new = sys.argv[2], sys.argv[3]
text = path.read_text()
count = text.count(old)
if count != 1:
    print(f"expected the pattern exactly once, found {count}", file=sys.stderr)
    raise SystemExit(3)
path.write_text(text.replace(old, new, 1))
PY
    if (( $? != 0 )) || cmp -s "$scratch/pristine.cpp" "$target"; then
        printf '%s  SKIP      %s (corruption did not apply)%s\n' "$yellow" "$1" "$reset"
        skipped=$((skipped + 1)); return
    fi
    if bash "$script_dir/verify-gaslift.sh" compare "$golden" > /dev/null 2>&1; then
        printf '%s  MISSED    %s%s\n' "$red" "$1" "$reset"; missed=$((missed + 1))
    else
        printf '%s  detected  %s%s\n' "$green" "$1" "$reset"; detected=$((detected + 1))
    fi
}

probe 'areaValvCali: calibration blend'   'PB80 = (PB80 + 14.6959488) * (80 + 460.67)' \
                                          'PB80 = (PB80 + 14.6959488) * (80 + 460.68)'
probe 'areaValvCali: opening cap'         'if (APE > areagarg)' 'if (APE > 2. * areagarg)'
# The same assignment appears in renovaGas and renovaGasBuf; the preceding line
# differs and pins this one to renovaGas.
probe 'renovaGas: neighbour source'       'celulaG[i].pres = termolivreG[3 * i];
            celulaG[i].presL = termolivreG[3 * i - 3];
            celulaG[i].presR = termolivreG[3 * i + 3];' \
                                          'celulaG[i].pres = termolivreG[3 * i];
            celulaG[i].presL = termolivreG[3 * i - 2];
            celulaG[i].presR = termolivreG[3 * i + 3];'
probe 'prescordesc: sign'                 'return sinal * precorr;' 'return -sinal * precorr;'
probe 'prescordesc: throat area'          'pow(massica / chokeVGL[ivalv].areagarg, 2.)' \
                                          'pow(massica / chokeVGL[ivalv].areagarg, 3.)'
probe 'delpGasPerm: hydrostatic constant' 'double gradhidro = celulaG[i].dPdLHidro * (9.82 * sin(celulaG[i].duto.teta) * rhog * dx);' \
                                          'double gradhidro = celulaG[i].dPdLHidro * (9.81 * sin(celulaG[i].duto.teta) * rhog * dx);'
# delpGasPerm and delpInjPerm both convert with 98066.5; the gas one is preceded
# by its own hydrostatic term written on celulaG.
probe 'delpGasPerm: pressure unit'        'double gradhidro = celulaG[i].dPdLHidro * (9.82 * sin(celulaG[i].duto.teta) * rhog * dx);

    double difpres = (gradfric + gradhidro) / 98066.5;' \
                                          'double gradhidro = celulaG[i].dPdLHidro * (9.82 * sin(celulaG[i].duto.teta) * rhog * dx);

    double difpres = (gradfric + gradhidro) / 98600.;'
probe 'delpInjPerm: interpolation weight' 'tmed = (celula[i].dx * celula[i].temp + celula[i].dxL * celula[i - 1].temp) / (celula[i].dx + celula[i].dxL);' \
                                          'tmed = (celula[i].dxL * celula[i].temp + celula[i].dx * celula[i - 1].temp) / (celula[i].dx + celula[i].dxL);'

cp "$scratch/pristine.cpp" "$target"
printf '\n%s cases: %s detected, %s missed, %s skipped\n' \
    "$((detected + missed + skipped))" "$detected" "$missed" "$skipped"
if (( missed > 0 || skipped > 0 )); then
    printf '%sGASLIFT CALIBRATION FAILED%s\n' "$red" "$reset" >&2; exit 1
fi
printf '%sGASLIFT CALIBRATION PASSED -- %s cases%s\n' "$green" "$detected" "$reset"
