#!/usr/bin/env bash
# Compare the CURRENT tree's dt series against the one T122 captured.
#
# This is the check L2 cannot perform. L2 compares final outputs: if the time
# step drifts by a last bit, L2 says "the outputs differ" and stops there. The
# series says WHERE the temporal mesh first diverged -- and since the step feeds
# back into everything, the first different step makes every later one
# incomparable, so "where" is the only actionable information.
#
# Usage: verify-determinadt.sh [golden-file]
# Exit code: 0 only when the series is identical line for line.
set -uo pipefail
export LC_ALL=C

project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
golden="${1:-$project_root/specs/001-refatoracao-sisprod/golden/determinadt.txt}"
red=$'\033[0;31m'; green=$'\033[0;32m'; yellow=$'\033[1;33m'; reset=$'\033[0m'

[[ -s "$golden" ]] || { printf '%sgolden ausente ou vazio: %s%s\n' "$red" "$golden" "$reset" >&2; exit 2; }

current="$(mktemp -t marlim3-dtnow-XXXXXX)"
trap 'rm -f "$current"' EXIT

bash "$project_root/refactor-harness/capture-determinadt.sh" "$current" > /dev/null || {
    printf '%sa captura falhou%s\n' "$red" "$reset" >&2; exit 2; }

golden_lines="$(wc -l < "$golden")"
current_lines="$(wc -l < "$current")"

if cmp -s "$golden" "$current"; then
    printf '%sDT SERIES IDENTICAL -- %s calls, bit for bit%s\n' "$green" "$golden_lines" "$reset"
    # A series that matches because both sides are empty proves nothing, and
    # six of the fourteen corpus models never call determinaDT at all.
    models="$(awk '{print $1}' "$golden" | sort -u | wc -l)"
    printf '%scoverage: %s of 14 corpus models reach determinaDT at all%s\n' "$yellow" "$models" "$reset"
    exit 0
fi

printf '%sDT SERIES DIFFERS%s\n' "$red" "$reset" >&2
printf '  golden : %s linhas\n' "$golden_lines" >&2
printf '  atual  : %s linhas\n' "$current_lines" >&2
diff <(cat "$golden") <(cat "$current") | head -8 >&2
first="$(diff --unchanged-line-format= --old-line-format='%dn ' --new-line-format= "$golden" "$current" 2>/dev/null | awk '{print $1; exit}')"
[[ -n "$first" ]] && printf '  primeira divergencia na linha %s: %s\n' "$first" "$(sed -n "${first}p" "$golden")" >&2
exit 1
