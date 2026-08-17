#!/usr/bin/env bash
# Run the six constitutional gates in one pass and archive the evidence.
#
# A stage is not complete until all six pass WITH EXECUTED OUTPUT. Reasoning,
# inspection or the expectation that results did not change satisfies none of
# them. This script exists so that "I ran the gates" means the same thing every
# time, and so the evidence lands somewhere reviewable.
#
# Gate 7 (time-step series) is specific to stage 8 and is not run here.
#
# Usage:
#   run-gates.sh <evidence-directory> [model ...]
#
# With no models, the full corpus is used -- which is what a stage boundary
# requires. Naming models gives the fast cycle for intermediate checks.
#
# Exit code: 0 only when every gate passes.

set -uo pipefail
export LC_ALL=C

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "$script_dir/.." && pwd)"

BASELINE_DIR="${MARLIM_BASELINE:-$HOME/marlim3-baseline}"

evidence_dir="${1:?usage: run-gates.sh <evidence-directory> [model ...]}"
shift || true
models=("$@")

mkdir -p "$evidence_dir"
cd "$project_root" || exit 2

red=$'\033[0;31m'; green=$'\033[0;32m'; yellow=$'\033[1;33m'; bold=$'\033[1m'; reset=$'\033[0m'
failures=0

announce() { printf '\n%s=== %s ===%s\n' "$bold" "$1" "$reset"; }
verdict() {
    if (( $2 == 0 )); then
        printf '%sGATE %s: PASS%s\n' "$green" "$1" "$reset"
    else
        printf '%sGATE %s: FAIL%s\n' "$red" "$1" "$reset"
        failures=$((failures + 1))
    fi
}

# ------------------------------------------------------------ pre-flight ----
# Before any gate runs, check whether the environment still matches the one the
# baseline came from. This is not one of the six constitutional gates; it is the
# question that has to be answered before their results mean anything.
#
# When it went unasked, a system change invalidated the baseline silently and
# the next full run failed L2 on a binary built from unmodified source. That
# failure is indistinguishable from a refactoring bug until you know to suspect
# the baseline, so the check runs first and says which it is.
announce "Pre-flight - is the baseline still valid on this machine?"
bash "$script_dir/record-provenance.sh" --check 2>&1 | tee "$evidence_dir/provenance.log"
provenance_status="${PIPESTATUS[0]}"

if (( provenance_status != 0 )); then
    printf '%sThe environment differs from the one that produced the baseline.%s\n' "$yellow" "$reset"
    printf '%sRun verify-baseline.sh before trusting any gate below.%s\n' "$yellow" "$reset"
fi

# ---------------------------------------------------------------- gate 1 ----
announce "Gate 1 - clean build, no new warnings"
cmake --build --preset gcc-release > "$evidence_dir/build.log" 2>&1
build_status=$?
errors=$(grep -c "error:" "$evidence_dir/build.log")
compiled=$(grep -c "Building CXX object" "$evidence_dir/build.log")

# Warnings are compared as a SET OF IDENTITIES (file:line:column plus warning
# type), never as a count.
#
# Counting does not work here. A parallel build interleaves compiler output, so
# the same source can produce a different number of log lines from run to run --
# the baseline recorded 472 and an identical rebuild produced 485, with zero
# warnings disappearing, which is the signature of lost output rather than of
# new diagnostics. Worse, a count says nothing useful even when accurate: as
# SisProd.cpp is decomposed its 249 warnings migrate to the new modules, so the
# total is expected to move while no warning is actually new.
#
# What "no new warnings" means is that no file:line:type appears that was not
# there before. That is what this compares.
warning_identities() {
    grep "warning:" "$1" \
        | sed -E 's/^.*(src\/[^:]+:[0-9]+:[0-9]+).*(\[-W[a-z-]+\]).*/\1 \2/' \
        | grep '^src/' | sort -u
}

baseline_build_log="${MARLIM_BASELINE_BUILD_LOG:-$BASELINE_DIR/logs/build-baseline.log}"
current_warnings="$evidence_dir/warnings-current.txt"
new_warnings="$evidence_dir/warnings-new.txt"

warning_identities "$evidence_dir/build.log" > "$current_warnings"

if [[ -f "$baseline_build_log" ]]; then
    comm -13 <(warning_identities "$baseline_build_log") "$current_warnings" > "$new_warnings"
    introduced=$(wc -l < "$new_warnings")
else
    : > "$new_warnings"
    introduced=0
    printf '%sno baseline build log at %s -- cannot check for new warnings%s\n' \
           "$yellow" "$baseline_build_log" "$reset"
fi

printf 'compiled=%s errors=%s warnings=%s new=%s\n' \
       "$compiled" "$errors" "$(wc -l < "$current_warnings")" "$introduced"

if (( introduced > 0 )); then
    printf '%snew warning identities:%s\n' "$red" "$reset"
    head -10 "$new_warnings" | sed 's/^/      /'
fi

# An incremental build with nothing to do emits zero warnings and would pass the
# comparison without having verified anything. Say so, rather than reporting a
# pass that carries no information.
if (( compiled == 0 )); then
    printf '%sno translation unit was recompiled -- this gate is vacuous.%s\n' "$yellow" "$reset"
    printf '%sRun a clean build before the stage boundary.%s\n' "$yellow" "$reset"
fi

(( build_status == 0 && errors == 0 && introduced == 0 ))
verdict 1 $?

# ---------------------------------------------------------------- gate 2 ----
announce "Gate 2 - bit-for-bit equivalence (L2)"
bash "$script_dir/compare-l2.sh" "${models[@]}" 2>&1 | tee "$evidence_dir/l2.log"
verdict 2 "${PIPESTATUS[0]}"

# ---------------------------------------------------------------- gate 3 ----
announce "Gate 3 - regression suite (L3)"
# Piped into tee rather than redirected to a file: `uv run` exits 120 and
# produces no output at all when its stdout is a regular file, while the same
# command through a pipe runs normally. Redirecting made this gate report FAIL
# on a suite that passes -- a false failure, which is the one kind of gate
# result that trains people to stop trusting the gate.
uv run pytest tests/test_regression.py -m regressao 2>&1 | tee "$evidence_dir/l3.log"
verdict 3 "${PIPESTATUS[0]}"
tail -1 "$evidence_dir/l3.log"

# ---------------------------------------------------------------- gate 4 ----
announce "Gate 4 - performance within threshold"
bash "$script_dir/performance-gate.sh" "${models[@]}" 2>&1 | tee "$evidence_dir/performance.log"
verdict 4 "${PIPESTATUS[0]}"

# ---------------------------------------------------------------- gate 5 ----
announce "Gate 5 - reference files untouched"
git status --porcelain tests/comparison/ > "$evidence_dir/references.log" 2>&1
[[ ! -s "$evidence_dir/references.log" ]]
verdict 5 $?

# ---------------------------------------------------------------- gate 6 ----
announce "Gate 6 - consumers compile"
consumers_broken=0
for consumer in Num4Main.cpp FA_Hidratos.cpp FA_Hidratos_Servico.cpp SisProdVap.cpp; do
    if grep -q "error:.*$consumer" "$evidence_dir/build.log"; then
        printf '%s%s failed to compile%s\n' "$red" "$consumer" "$reset"
        consumers_broken=1
    else
        printf '%s ok\n' "$consumer"
    fi
done
(( consumers_broken == 0 ))
verdict 6 $?

# ------------------------------------------------------------------ summary --
printf '\n%s----------------------------------------------------------------%s\n' "$bold" "$reset"
printf 'evidence archived in: %s\n' "$evidence_dir"

if (( provenance_status != 0 )); then
    printf '%snote: the environment changed since the baseline was captured.%s\n' "$yellow" "$reset"
    printf '%sIf a gate failed, confirm with verify-baseline.sh before blaming the code.%s\n' \
           "$yellow" "$reset"
fi

if (( failures > 0 )); then
    printf '%s%s of 6 gates FAILED -- the stage is not complete%s\n' "$red" "$failures" "$reset" >&2
    exit 1
fi

printf '%sall 6 gates passed%s\n' "$green" "$reset"
exit 0
