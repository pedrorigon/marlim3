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

# Warnings are compared as PER-FLAG GLOBAL COUNTS. No flag may increase.
#
# Three schemes came before this one, and each failed on the thing this program
# actually does.
#
# A raw total fails because a parallel build interleaves compiler output: the
# baseline recorded 472 where an identical rebuild produced 485, with zero
# warnings disappearing -- lost log lines, not new diagnostics (E-08).
#
# Replacing the total with file:line:column identities fixed that and broke
# under code movement. Stage 1 lifted 342 lines out of SisProd.cpp, every
# warning below them shifted line, and the gate reported 250 new diagnostics
# while the multiset was conserved exactly, 485 to 485.
#
# Keying on (file, message) fails too: a warning follows its code into the new
# module, and renaming a local rewrites the message text. Both are deliberate
# steps of every extraction stage.
#
# What survives moving code between files and renaming locals is the count per
# warning flag. A genuinely new diagnostic raises one of those counts; a
# migration does not. The blind spot is real and declared: one warning vanishing
# while another appears under the same flag would cancel out. The per-file
# breakdown is written to disk for exactly that case -- it is what made the
# stage-1 migration legible.
#
# GCC quotes identifiers typographically under a UTF-8 locale and in ASCII under
# LC_ALL=C. The baseline log was captured under the former and this build runs
# under the latter, so message text is normalized before it is compared (E-04).
normalize_quotes() { sed -e "s/\xe2\x80\x98/'/g; s/\xe2\x80\x99/'/g" "$1"; }

# Not every diagnostic carries a [-Wflag]: the rapidjson C++20 comparison
# warning has none. Bucketing those under [no-flag] rather than dropping them
# keeps a new unflagged warning from being invisible to this gate.
warning_flag_counts() {
    normalize_quotes "$1" | grep "warning:" \
        | sed -E 's/^.*(\[-W[a-z-]+\]).*/\1/; t; s/.*/[no-flag]/' \
        | sort | uniq -c | awk '{ printf "%s %s\n", $2, $1 }' | sort
}

warning_details() {
    normalize_quotes "$1" | grep "warning:" \
        | sed -E 's|^.*/(src/[^:]+):[0-9]+:[0-9]+: warning: (.*)$|\1 :: \2|' \
        | grep '^src/' | sort | uniq -c | sed 's/^ *//'
}

baseline_build_log="${MARLIM_BASELINE_BUILD_LOG:-$BASELINE_DIR/logs/build-baseline.log}"
current_flags="$evidence_dir/warning-flags-current.txt"
baseline_flags="$evidence_dir/warning-flags-baseline.txt"
new_warnings="$evidence_dir/warnings-new.txt"

warning_flag_counts "$evidence_dir/build.log" > "$current_flags"
total_warnings=$(awk '{ sum += $2 } END { print sum + 0 }' "$current_flags")

if [[ -f "$baseline_build_log" ]]; then
    warning_flag_counts "$baseline_build_log" > "$baseline_flags"
    # Any flag whose count grew, including flags absent from the baseline.
    join -a1 -a2 -e 0 -o 0,1.2,2.2 "$baseline_flags" "$current_flags" \
        | awk '$3 > $2 { printf "%s %d -> %d\n", $1, $2, $3 }' > "$new_warnings"
    introduced=$(wc -l < "$new_warnings")
    # Advisory: where warnings moved, which a per-flag count cannot show.
    diff <(warning_details "$baseline_build_log") <(warning_details "$evidence_dir/build.log") \
        > "$evidence_dir/warning-migration.txt" 2>&1 || true
else
    : > "$new_warnings"
    introduced=0
    printf '%sno baseline build log at %s -- cannot check for new warnings%s\n' \
           "$yellow" "$baseline_build_log" "$reset"
fi

printf 'compiled=%s errors=%s warnings=%s flags-increased=%s\n' \
       "$compiled" "$errors" "$total_warnings" "$introduced"

if (( introduced > 0 )); then
    printf '%swarning flags that increased:%s\n' "$red" "$reset"
    sed 's/^/      /' "$new_warnings"
    printf '%ssee %s for where each warning moved%s\n' \
           "$yellow" "$evidence_dir/warning-migration.txt" "$reset"
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

# ------------------------------------------------------ supplementary L0 ----
# Not a constitutional gate, but reported alongside them because for most of
# this file it is the only verification that exists.
#
# Coverage measurement found that the corpus exercises 29.3% of the executable
# lines of SisProd.cpp and never runs 81 of its 153 functions. For those, gates
# 2 and 3 compare the output of code that did not execute, and stay green
# whatever the extraction did to it. L0 compares the moved body token for token
# against the baseline commit, which proves the move was literal even when
# nothing runs it.
#
# It is advisory here rather than blocking, because a stage legitimately renames
# locals and decomposes functions, both of which change the token stream. The
# extraction tasks call it directly on the specific functions they move, where
# the expected answer is exact equality.
announce "Supplementary - structural comparison against the baseline commit (L0)"
baseline_source="$(mktemp -d -t marlim3-l0-XXXXXX)"
if git -C "$project_root" worktree add --detach "$baseline_source/tree" \
       "$(cat "$BASELINE_DIR/commit" 2>/dev/null || echo 0f3b64f)" > /dev/null 2>&1; then
    # Compare against SisProd.cpp AND every module this refactoring introduced,
    # so a function that moved as planned is found and compared rather than
    # reported as having vanished. Extracted modules are exactly the src/core
    # sources that do not exist in the baseline commit, which keeps this correct
    # as later stages add more of them.
    l0_targets=(--current "$project_root/src/core/SisProd.cpp")
    for candidate in "$project_root"/src/core/*.cpp "$project_root"/src/include/*.h; do
        base="$(basename "$candidate")"
        [[ -e "$baseline_source/tree/src/core/$base" || -e "$baseline_source/tree/src/include/$base" ]] \
            || l0_targets+=(--current "$candidate")
    done
    # Headers are searched too, and not for symmetry. Stage 2 moved zbrent,
    # falsacorda and zriddr into RootFindingSolvers.h, because a template
    # parameterised by the objective has to be defined where it is instantiated.
    # Looking only at src/core/*.cpp would report all three as MISSING -- the
    # exact failure this loop was written to prevent, one stage after it was
    # written.
    # --declared keeps the functions each stage restructured on purpose out of
    # the failure count, so the failures that remain are the ones nobody
    # planned. Without it the report accumulates one permanent difference per
    # decomposed function, and by the last stage a function that vanished for a
    # bad reason is indistinguishable from one that moved as designed.
    python3 "$script_dir/verify-structural.py" \
        --baseline "$baseline_source/tree/src/core/SisProd.cpp" \
        "${l0_targets[@]}" \
        --declared "$script_dir/decomposed-functions.txt" \
        --all 2>&1 | grep -vE '^OK ' | tail -40 | tee "$evidence_dir/l0.log"
    git -C "$project_root" worktree remove --force "$baseline_source/tree" > /dev/null 2>&1
else
    printf '%scould not create a worktree for the baseline commit%s\n' "$yellow" "$reset"
fi
rm -rf "$baseline_source"

# ------------------------------------------- supplementary calibration ----
# Every instrument, not only the one belonging to the stage being closed.
#
# These checks pin fixed text, and this programme moves and renames that text,
# so each stage erodes the instruments of the stages before it -- silently, since
# a pattern that no longer matches does not crash, it just stops testing. Closing
# stage 5 found calibrate-thermal.sh with 23 of 24 cases unable to inject at all,
# and calibrate-c0ud.sh with 7 of 21 dead since stage 3. Both had been reporting
# nothing wrong, because nothing was being tested.
#
# Advisory rather than blocking: a stage in progress can legitimately leave an
# instrument mid-repoint. It is reported next to the gates because that is where
# someone will read it.
announce "Supplementary - every calibration still detects error"
bash "$script_dir/calibrate-all.sh" 2>&1 | tee "$evidence_dir/calibration.log" | tail -12

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
