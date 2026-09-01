#!/usr/bin/env bash
# Reproduce the T067 whitespace-only comparison of the steady thermal marches.
set -euo pipefail
export LC_ALL=C

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "$script_dir/.." && pwd)"
source_file="${1:-$project_root/src/core/SisProd.cpp}"
diff_file="${2:-}"

[[ -f "$source_file" ]] || {
    echo "no source at $source_file" >&2
    exit 2
}

scratch="$(mktemp -d -t marlim3-renovatemp-diff-XXXXXX)"
trap 'rm -rf "$scratch"' EXIT

extract_function() {
    local function_name="$1" destination="$2"
    python3 - "$source_file" "$function_name" "$destination" <<'PY'
import pathlib
import sys

source = pathlib.Path(sys.argv[1]).read_text()
function_name = sys.argv[2]
destination = pathlib.Path(sys.argv[3])
signature = f"void SProd::{function_name}("
start = source.find(signature)
if start < 0:
    raise SystemExit(f"function not found: {function_name}")

brace = source.find("{", start)
if brace < 0:
    raise SystemExit(f"opening brace not found: {function_name}")

depth = 0
state = "code"
index = brace
while index < len(source):
    char = source[index]
    next_char = source[index + 1] if index + 1 < len(source) else ""
    if state == "code":
        if char == "/" and next_char == "/":
            state = "line-comment"
            index += 1
        elif char == "/" and next_char == "*":
            state = "block-comment"
            index += 1
        elif char == '"':
            state = "string"
        elif char == "'":
            state = "character"
        elif char == "{":
            depth += 1
        elif char == "}":
            depth -= 1
            if depth == 0:
                destination.write_text(source[start:index + 1] + "\n")
                break
    elif state == "line-comment":
        if char == "\n":
            state = "code"
    elif state == "block-comment":
        if char == "*" and next_char == "/":
            state = "code"
            index += 1
    elif state in {"string", "character"}:
        delimiter = '"' if state == "string" else "'"
        if char == "\\":
            index += 1
        elif char == delimiter:
            state = "code"
    index += 1
else:
    raise SystemExit(f"closing brace not found: {function_name}")
PY
}

extract_function RenovaTempPerm "$scratch/forward.cpp"
extract_function RenovaTempPermRev "$scratch/reverse.cpp"
sed 's/[[:space:]]//g' "$scratch/forward.cpp" > "$scratch/forward.normalized.cpp"
sed 's/[[:space:]]//g' "$scratch/reverse.cpp" > "$scratch/reverse.normalized.cpp"

set +e
diff -u --label RenovaTempPerm.normalized \
    --label RenovaTempPermRev.normalized \
    "$scratch/forward.normalized.cpp" "$scratch/reverse.normalized.cpp" \
    > "$scratch/normalized.diff"
diff_status=$?
set -e
(( diff_status == 0 || diff_status == 1 )) || exit "$diff_status"

if [[ -n "$diff_file" ]]; then
    mkdir -p "$(dirname "$diff_file")"
    cp "$scratch/normalized.diff" "$diff_file"
    printf 'wrote %s normalized diff lines to %s\n' \
        "$(wc -l < "$scratch/normalized.diff")" "$diff_file"
else
    cat "$scratch/normalized.diff"
fi
