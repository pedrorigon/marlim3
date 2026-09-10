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
# Two targets, because the gas-line update no longer lives where the coefficient
# routines do: stage 6 moved renovaGas and renovaGasBuf into SisProdGasLift.cpp
# and renamed their state along the way. A probe pins a fixed string, so a probe
# left pointing at the old file stops injecting and stops testing -- which is
# exactly what happened here and what the SKIP path exists to announce.
target="$project_root/src/core/SisProd.cpp"
target_gaslift="$project_root/src/core/SisProdGasLift.cpp"
golden="$project_root/specs/001-refatoracao-sisprod/golden/gaslift-baseline.txt"

red=$'\033[0;31m'; green=$'\033[0;32m'; yellow=$'\033[1;33m'; reset=$'\033[0m'
scratch="$(mktemp -d -t marlim3-calgl-XXXXXX)"
trap 'cp "$scratch/pristine.cpp" "$target"; cp "$scratch/pristine-gaslift.cpp" "$target_gaslift"; rm -rf "$scratch"' EXIT
cp "$target" "$scratch/pristine.cpp"
cp "$target_gaslift" "$scratch/pristine-gaslift.cpp"

detected=0; missed=0; skipped=0
probe() { # <label> <from> <to> [file-to-inject-into]
    local into="${4:-$target}" pristine="$scratch/pristine.cpp"
    [[ "$into" == "$target_gaslift" ]] && pristine="$scratch/pristine-gaslift.cpp"
    cp "$scratch/pristine.cpp" "$target"
    cp "$scratch/pristine-gaslift.cpp" "$target_gaslift"
    python3 - "$into" "$2" "$3" <<'PY'
import pathlib, sys
path = pathlib.Path(sys.argv[1]); old, new = sys.argv[2], sys.argv[3]
text = path.read_text()
count = text.count(old)
if count != 1:
    print(f"expected the pattern exactly once, found {count}", file=sys.stderr)
    raise SystemExit(3)
path.write_text(text.replace(old, new, 1))
PY
    if (( $? != 0 )) || cmp -s "$pristine" "$into"; then
        printf '%s  SKIP      %s (corruption did not apply)%s\n' "$yellow" "$1" "$reset"
        skipped=$((skipped + 1)); return
    fi
    if bash "$script_dir/verify-gaslift.sh" compare "$golden" > /dev/null 2>&1; then
        printf '%s  MISSED    %s%s\n' "$red" "$1" "$reset"; missed=$((missed + 1))
    else
        printf '%s  detected  %s%s\n' "$green" "$1" "$reset"; detected=$((detected + 1))
    fi
}

probe 'areaValvCali: calibration blend'   'bellowsPressureAt80F = (bellowsPressureAt80F + 14.6959488) * (80 + 460.67) / (calibrationTemperature * 1.8 + 491.67) - 14.6959488;' \
                                          'bellowsPressureAt80F = (bellowsPressureAt80F + 14.6959488) * (80 + 460.68) / (calibrationTemperature * 1.8 + 491.67) - 14.6959488;' \
                                          "$target_gaslift"
probe 'areaValvCali: opening cap'         'if (openingArea > throatArea)' \
                                          'if (openingArea > 2. * throatArea)' \
                                          "$target_gaslift"
# updateGasLine writes this assignment twice -- once per branch of its interior
# if/else. Only the interior branch follows presL with presR, so carrying the
# next line pins the corruption to that one.
probe 'updateGasLine: neighbour source'   'state.gasCells[gasCellIndex].presL = state.gasFreeTerms[3 * gasCellIndex - 3];
            state.gasCells[gasCellIndex].presR = state.gasFreeTerms[3 * gasCellIndex + 3];' \
                                          'state.gasCells[gasCellIndex].presL = state.gasFreeTerms[3 * gasCellIndex - 2];
            state.gasCells[gasCellIndex].presR = state.gasFreeTerms[3 * gasCellIndex + 3];' \
                                          "$target_gaslift"
# ------------------------------------------------------- unloading probes ---
# The six routines below were, until the sweep was extended, executed by no
# verification layer at all: the corpus never sets condicaoInicial == 3 and the
# sweep did not drive them. Rows without probes would only move the blind spot
# from "never executed" to "executed, never checked", so each gets one.
#
# Every corruption here is a LITERAL or a comparison bound. Renaming cannot
# excuse any of them, and each pattern is pinned to text that occurs exactly
# once -- advanceBufferedGasSubStep needs a two-line window, because its body is
# a near-duplicate of advanceGasSubStep and every single line it contains also
# appears there.
probe 'HidroDescargaG: hydrostatic head' \
    'meanPressure -= rho1 * 9.81 * halfLocalLength * sin(state.gasCells[gasCellIndex].duto.teta) / 98066.52;' \
    'meanPressure -= rho1 * 9.82 * halfLocalLength * sin(state.gasCells[gasCellIndex].duto.teta) / 98066.52;' \
    "$target_gaslift"

probe 'CalcPresValvDesc: wall shear' \
    'double tens1 = frictionFactor * mixtureDensity * vel1 * fabs(vel1) / 2.;' \
    'double tens1 = frictionFactor * mixtureDensity * vel1 * fabs(vel1) / 2.5;' \
    "$target_gaslift"

probe 'BuscaPresInjDesc: pressure decay' \
    'state.gasSurfacePressure *= (1 - 0.01 * state.timeStep);' \
    'state.gasSurfacePressure *= (1 - 0.02 * state.timeStep);' \
    "$target_gaslift"

# The hand-over branch, which only the second advancInter row reaches. If that
# row were ever dropped this probe would start missing, which is the point.
probe 'avancInter: handover ratio' \
    'state.gasCells[state.interfaceCell + 1].razInter = 1.0;' \
    'state.gasCells[state.interfaceCell + 1].razInter = 0.999;' \
    "$target_gaslift"

probe 'resolveDescarga: liquid hydrostatic' \
    'double hidro1L = 1 * (9.82 * sin(state.gasCells[gasCellIndex - 1].duto.teta) * rhoL) * upstreamLiquidLength;' \
    'double hidro1L = 1 * (9.81 * sin(state.gasCells[gasCellIndex - 1].duto.teta) * rhoL) * upstreamLiquidLength;' \
    "$target_gaslift"

# Splits the gas/liquid share of every control volume. Moving the threshold
# changes which side of the interface each cell is counted on.
probe 'resolveDescarga: interface split' \
    'if (state.gasCells[gasCellIndex].razInter <= 0.5)' \
    'if (state.gasCells[gasCellIndex].razInter <= 0.6)' \
    "$target_gaslift"

probe 'subtempoGasBuf: choke opening bound' \
    '        chokeOpeningFraction = state.injectionChoke.areagarg / state.gasCells[0].duto.area;
        if (chokeOpeningFraction < 0.2) {' \
    '        chokeOpeningFraction = state.injectionChoke.areagarg / state.gasCells[0].duto.area;
        if (chokeOpeningFraction < 0.3) {' \
    "$target_gaslift"

# updateBufferedGasLine writes one field, in three byte-identical branches.
# Without a probe on the offset itself the whole routine was unguarded.
probe 'renovaGasBuf: buffer offset' \
    '        state.gasCells[gasCellIndex].VGasRBuf = state.gasFreeTerms[3 * gasCellIndex + 1];
}' \
    '        state.gasCells[gasCellIndex].VGasRBuf = state.gasFreeTerms[3 * gasCellIndex + 2];
}' \
    "$target_gaslift"

probe 'prescordesc: sign'                 'return sign * pressureCorrection;' 'return -sign * pressureCorrection;' \
                                          "$target_gaslift"
probe 'prescordesc: throat area'          'pressureCorrection = pow(massFlowRate / state.gasLiftChokes[valveIndex].areagarg, 2.) / (2. * rho0 * 98066.52);' \
                                          'pressureCorrection = pow(massFlowRate / state.gasLiftChokes[valveIndex].areagarg, 3.) / (2. * rho0 * 98066.52);' \
                                          "$target_gaslift"
probe 'delpGasPerm: hydrostatic constant' 'double hydrostaticGradient = state.gasCells[cellIndex].dPdLHidro * (9.82 * sin(state.gasCells[cellIndex].duto.teta) * gasDensity * dx);' \
                                          'double hydrostaticGradient = state.gasCells[cellIndex].dPdLHidro * (9.81 * sin(state.gasCells[cellIndex].duto.teta) * gasDensity * dx);' \
                                          "$target_gaslift"
# Both steady drops divide by 98066.5; the hydrostatic line above pins this one
# to delpGasPerm.
probe 'delpGasPerm: pressure unit'        'double hydrostaticGradient = state.gasCells[cellIndex].dPdLHidro * (9.82 * sin(state.gasCells[cellIndex].duto.teta) * gasDensity * dx);

    double pressureDrop = (frictionGradient + hydrostaticGradient) / 98066.5;' \
                                          'double hydrostaticGradient = state.gasCells[cellIndex].dPdLHidro * (9.82 * sin(state.gasCells[cellIndex].duto.teta) * gasDensity * dx);

    double pressureDrop = (frictionGradient + hydrostaticGradient) / 98066.52;' \
                                          "$target_gaslift"
probe 'delpInjPerm: interpolation weight' 'meanTemperature = (state.cells[cellIndex].dx * state.cells[cellIndex].temp + state.cells[cellIndex].dxL * state.cells[cellIndex - 1].temp) / (state.cells[cellIndex].dx + state.cells[cellIndex].dxL);' \
                                          'tmed = (state.cells[cellIndex].dxL * state.cells[cellIndex].temp + state.cells[cellIndex].dx * state.cells[cellIndex - 1].temp) / (state.cells[cellIndex].dx + state.cells[cellIndex].dxL);' \
                                          "$target_gaslift"

cp "$scratch/pristine.cpp" "$target"
cp "$scratch/pristine-gaslift.cpp" "$target_gaslift"

# ---------------------------------------------------- token-check control ----
# T081r renamed the module's locals, so gaslift-move.py check now accepts a
# consistent renaming instead of demanding an exact match. A tolerance is a
# blind spot until something proves otherwise, so this control corrupts a
# LITERAL -- which no renaming can excuse -- and requires the check to reject
# it. It also confirms the check still passes on the pristine module, so a
# check that had broken outright could not masquerade as a detection.
token_baseline="$scratch/token-baseline.cpp"
if git -C "$project_root" show "${MARLIM_T079_BASELINE:-66ed35c~1}:src/core/SisProd.cpp" \
       > "$token_baseline" 2>/dev/null; then
    if python3 "$script_dir/gaslift-move.py" check delpGasPerm \
           "$token_baseline" "$target_gaslift" > /dev/null 2>&1; then
        sed -i 's/9\.82 \* sin/9.81 * sin/' "$target_gaslift"
        if python3 "$script_dir/gaslift-move.py" check delpGasPerm \
               "$token_baseline" "$target_gaslift" > /dev/null 2>&1; then
            printf '%s  MISSED    token check: literal changed under rename tolerance%s\n' \
                   "$red" "$reset"; missed=$((missed + 1))
        else
            printf '%s  detected  token check: literal changed under rename tolerance%s\n' \
                   "$green" "$reset"; detected=$((detected + 1))
        fi
        cp "$scratch/pristine-gaslift.cpp" "$target_gaslift"
    else
        printf '%s  SKIP      token check control (pristine module does not pass)%s\n' \
               "$yellow" "$reset"; skipped=$((skipped + 1))
    fi
else
    printf '%s  SKIP      token check control (baseline commit unreadable)%s\n' \
           "$yellow" "$reset"; skipped=$((skipped + 1))
fi

printf '\n%s cases: %s detected, %s missed, %s skipped\n' \
    "$((detected + missed + skipped))" "$detected" "$missed" "$skipped"
if (( missed > 0 || skipped > 0 )); then
    printf '%sGASLIFT CALIBRATION FAILED%s\n' "$red" "$reset" >&2; exit 1
fi
printf '%sGASLIFT CALIBRATION PASSED -- %s cases%s\n' "$green" "$detected" "$reset"
