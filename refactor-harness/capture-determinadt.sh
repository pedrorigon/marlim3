#!/usr/bin/env bash
# T122 -- capture the dt series that every later task in stage 8 is measured
# against.
#
# determinaDT decides the time step. A last-bit drift there does not produce a
# small difference in the answer: it produces a DIFFERENT TEMPORAL DISCRETISATION,
# and from that step onward the run is a different simulation. L2 at the end of a
# task would catch that, but only as "the outputs differ" -- it would not say
# where the mesh first diverged. The series does.
#
# The tree is never touched: everything happens in a throwaway copy, which is the
# rule instrument-golden.py established.
#
# Values are written with %a. Decimal would drop bits and turn a bit-for-bit
# comparison into an approximate one, which is the whole thing this file exists
# to prevent.
set -euo pipefail
export LC_ALL=C

project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
out="${1:-$project_root/specs/001-refatoracao-sisprod/golden/determinadt.txt}"
work="$(mktemp -d -t marlim3-dt-XXXXXX)"
trap 'rm -rf "$work"' EXIT

# Copy the whole tree, minus what is regenerated or huge. Listing the inputs by
# hand does not work and two attempts proved it: the build reads _version.py,
# then marlim3/translations.json, and there is no reason to believe that list is
# finished. An exclusion list fails safe; an inclusion list fails silently, one
# file at a time.
rsync -a --exclude build/ --exclude .git/ --exclude '.venv*' --exclude 'docs-sisprod/' \
      --exclude '*.o' "$project_root/" "$work/"
cd "$work"

python3 - <<'PY'
import io, re
p = "src/core/SisProd.cpp"
s = io.open(p, encoding="utf-8", errors="surrogateescape").read()

# determinaDT has ZERO return statements -- single exit -- so a log call before
# the closing brace sees every call. That is the drift-flux strategy from
# instrument-golden.py; the solver strategy (rename and wrap) is not needed here
# and would be the wrong tool.
m = re.search(r'^void SProd::determinaDT\(int vexpli\)\s*\{', s, re.M)
assert m, "determinaDT not found"
depth, i = 0, m.end() - 1
while True:
    if s[i] == "{":
        depth += 1
    elif s[i] == "}":
        depth -= 1
        if depth == 0:
            break
    i += 1

probe = '''
    {
        // T122. Appends one row per call: the time step this call produced, in
        // full precision, with enough context to tell two calls apart.
        static FILE *golden = fopen(getenv("MARLIM_DT_LOG"), "a");
        if (golden) {
            fprintf(golden, "%s %d %d %a %a %a\\n",
                    getenv("MARLIM_DT_MODEL") ? getenv("MARLIM_DT_MODEL") : "?",
                    vexpli, ncel, (*vg1dSP).lixo5, celula[0].pres, dt);
            fflush(golden);
        }
    }
'''
s = s[:i] + probe + s[i:]
if "#include <stdlib.h>" not in s:
    s = s.replace('#include "SisProd.h"', '#include "SisProd.h"\n#include <stdlib.h>', 1)
io.open(p, "w", encoding="utf-8", errors="surrogateescape").write(s)
print("probe injetada")
PY

grep -q "MARLIM_DT_LOG" src/core/SisProd.cpp || { echo "PROBE NAO APLICADA" >&2; exit 1; }

cmake --preset gcc-release > /dev/null 2>&1
cmake --build --preset gcc-release > /dev/null 2>&1 || { echo "build da copia falhou" >&2; exit 1; }

: > "$out"
export MARLIM_DT_LOG="$out"
for deck in demos/*.mr3 demos/pt-br/*.mr3; do
    [ -e "$deck" ] || continue
    model="$(basename "$deck" .mr3)"
    export MARLIM_DT_MODEL="$model"
    d="$work/run/$model"; mkdir -p "$d"
    ./build/Marlim3 -s TRANSIENTE -i "$deck" -p demos/ -d "$d" -o "$d/$model.log" > /dev/null 2>&1 || true
done

printf 'linhas capturadas: %s\n' "$(wc -l < "$out")"
awk '{print $1}' "$out" | sort | uniq -c | sort -rn | head -20
