#!/usr/bin/env bash
# Check that the stored baseline still describes what this machine produces.
#
# Why this exists
# ---------------
# The baseline is a set of outputs captured from a binary built at a point in
# time. Every later stage rebuilds the engine and compares against that capture.
# The comparison is only meaningful while an UNMODIFIED build still reproduces
# the captured outputs.
#
# That assumption broke once already, and silently. A binary built from the
# pristine baseline commit a week after the capture differed from the captured
# binary by 1.5 MB, linked libmvec where the original did not, and produced
# different results on 99 of the 102 files of
# extended-shutdown-combined-ESP-CGL-PIG-complete. Source, compiler version and
# flags were byte-identical; the environment underneath had moved. Two builds
# made on the same day agree with each other on all 102 files, so the current
# environment is self-consistent -- it simply is not the one the baseline came
# from.
#
# Discovered by accident, that failure looks exactly like a refactoring bug: L2
# fails, nobody knows why, and the search starts in the wrong place. This script
# turns it into a specific, early answer.
#
# What it does: builds the pristine baseline commit into a scratch tree and runs
# L2 with the resulting binary. Passing means the baseline is still valid on this
# machine. Failing means the baseline must be recaptured -- it does NOT mean the
# working tree is wrong.
#
# Run it at the start of a stage, after any toolchain or system update, and
# whenever L2 fails in a way that does not correspond to a change you made.
#
# Usage: verify-baseline.sh [model ...]

set -uo pipefail
export LC_ALL=C

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "$script_dir/.." && pwd)"

BASELINE_DIR="${MARLIM_BASELINE:-$HOME/marlim3-baseline}"
BASELINE_COMMIT="${MARLIM_BASELINE_COMMIT:-$(cat "$BASELINE_DIR/commit" 2>/dev/null || echo 0f3b64f)}"

red=$'\033[0;31m'; green=$'\033[0;32m'; yellow=$'\033[1;33m'; reset=$'\033[0m'

scratch="$(mktemp -d -t marlim3-baseline-check-XXXXXX)"
trap 'rm -rf "$scratch"' EXIT

printf 'building pristine %s in %s\n' "$BASELINE_COMMIT" "$scratch"

git -C "$project_root" worktree add --detach "$scratch/tree" "$BASELINE_COMMIT" > /dev/null 2>&1 || {
    printf '%scould not create a worktree at %s%s\n' "$red" "$BASELINE_COMMIT" "$reset" >&2
    exit 2
}
trap 'git -C "$project_root" worktree remove --force "$scratch/tree" >/dev/null 2>&1; rm -rf "$scratch"' EXIT

(
    cd "$scratch/tree" || exit 2
    cmake --preset gcc-release > /dev/null 2>&1 && \
    cmake --build --preset gcc-release > "$scratch/build.log" 2>&1
) || {
    printf '%spristine build failed -- see %s%s\n' "$red" "$scratch/build.log" "$reset" >&2
    exit 2
}

printf 'comparing pristine build against the stored baseline\n\n'

if MARLIM_BIN="$scratch/tree/build/Marlim3" bash "$script_dir/compare-l2.sh" "$@"; then
    printf '\n%sBASELINE VALID -- an unmodified build still reproduces the capture%s\n' \
           "$green" "$reset"
    exit 0
fi

printf '\n%sBASELINE STALE -- an unmodified build no longer reproduces the capture%s\n' \
       "$red" "$reset" >&2
printf '%sThis is not a defect in your working tree. The environment moved out from\n' "$yellow" >&2
printf 'under the capture, so the baseline must be recaptured before any stage\n' >&2
printf 'comparison means anything.%s\n' "$reset" >&2
exit 1
