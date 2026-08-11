#!/usr/bin/env bash
# Verification layer L1: compare extracted functions against golden values
# tabulated from the baseline.
#
# L1 exists so that a divergence is attributable to a FUNCTION, in seconds,
# instead of only showing up as "the model changed" after a 37-minute L2 run.
#
# Two tables, produced by instrument-golden.py against a throwaway copy of the
# baseline:
#
#   sweep-baseline.txt   synthetic sweep over a grid of physically plausible
#                        inputs, 2025 cases per drift-flux correlation
#   <function>.txt       inputs and outputs observed while running the corpus
#
# The sweep is not a nicety. Five of the eight target functions are NEVER
# called by the demo corpus -- Choi, BhagwatGhajarMod and FrancaLahey are
# unreachable because Leitura.cpp pins arq.CorreDisper to 1, and SProd::zbrent
# has no call site at all, which leaves SProd::falsacorda unreachable too. For
# those, corpus-driven capture yields nothing and L2 stays green even if the
# extraction corrupts them. The sweep is the only verification they get.
#
# Usage:
#   verify-l1.sh <binary-with-sweep>   compare that binary's sweep to the golden
#   verify-l1.sh                       use build/Marlim3
#
# Exit code: 0 when every value matches bit for bit; 1 otherwise.

set -uo pipefail
export LC_ALL=C

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "$script_dir/.." && pwd)"

GOLDEN_DIR="${MARLIM_GOLDEN_DIR:-$project_root/specs/001-refatoracao-sisprod/golden}"
GOLDEN_SWEEP="$GOLDEN_DIR/sweep-baseline.txt"
BINARY="${1:-${MARLIM_BIN:-$project_root/build/Marlim3}}"

red=$'\033[0;31m'; green=$'\033[0;32m'; yellow=$'\033[1;33m'; reset=$'\033[0m'

[[ -x "$BINARY" ]] || { printf '%sbinary missing: %s%s\n' "$red" "$BINARY" "$reset" >&2; exit 2; }
[[ -s "$GOLDEN_SWEEP" ]] || {
    printf '%sgolden sweep missing: %s (run T012)%s\n' "$red" "$GOLDEN_SWEEP" "$reset" >&2
    exit 2
}

work_dir="$(mktemp -d -t marlim3-l1-XXXXXX)"
trap 'rm -rf "$work_dir"' EXIT
current_sweep="$work_dir/sweep-current.txt"

# The sweep needs no model: the trigger fires before any input is parsed. A
# model path is still passed because the binary requires the argument.
MARLIM_GOLDEN_SWEEP="$current_sweep" \
    "$BINARY" -s TRANSIENTE -i demos/simplifiedProduction.mr3 -p demos/ \
              -d "$work_dir/out" -o "$work_dir/out.log" > /dev/null 2>&1

if [[ ! -s "$current_sweep" ]]; then
    printf '%sthe binary produced no sweep.%s\n' "$red" "$reset" >&2
    printf '%sIt must be built from a tree instrumented by instrument-golden.py.%s\n' \
           "$yellow" "$reset" >&2
    exit 2
fi

golden_lines=$(wc -l < "$GOLDEN_SWEEP")
current_lines=$(wc -l < "$current_sweep")

if [[ "$golden_lines" != "$current_lines" ]]; then
    printf '%sL1 FAILED -- case count differs: golden %s, current %s%s\n' \
           "$red" "$golden_lines" "$current_lines" "$reset" >&2
    exit 1
fi

if diff -q "$GOLDEN_SWEEP" "$current_sweep" > /dev/null 2>&1; then
    printf '%sL1 PASSED -- %s cases identical bit for bit%s\n' \
           "$green" "$golden_lines" "$reset"
    awk '{ seen[$1]++ } END { for (fn in seen) printf "  %-20s %6d cases\n", fn, seen[fn] }' \
        "$GOLDEN_SWEEP"
    exit 0
fi

printf '%sL1 FAILED -- divergent values%s\n' "$red" "$reset" >&2
echo >&2

# Localize: report the first divergence per function, which is what turns a
# failure into an address.
diff "$GOLDEN_SWEEP" "$current_sweep" \
    | grep '^<' | awk '{ if (!seen[$1]++) print }' | head -10 \
    | while read -r _ function_name rest; do
        printf '  first divergence in %s%s%s\n' "$red" "$function_name" "$reset" >&2
      done

divergent=$(diff "$GOLDEN_SWEEP" "$current_sweep" | grep -c '^<')
printf '\n  %s divergent case(s) of %s\n' "$divergent" "$golden_lines" >&2
exit 1
