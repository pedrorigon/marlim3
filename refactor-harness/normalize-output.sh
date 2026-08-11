#!/usr/bin/env bash
# Normalize engine output for equivalence comparison (verification layer L2).
#
# This exists because the engine writes four classes of content that vary
# between runs of the SAME binary on the SAME input, none of which is a
# computation result:
#
#   1. Timestamps, in four textual variants (two languages). Note that the
#      `datahora` form has NO leading zero on seconds ("19:41:3"), so the
#      pattern must accept 1 or 2 digits -- otherwise it fails intermittently,
#      only when the second is below 10.
#   2. Wall-clock duration written to the log ("DURACAO 386 segundos").
#   3. Decorative phrases drawn by rand() -- Num4Main.cpp:596 and 8084 run
#      `srand(time(NULL)); int frase = rand() % 16;`. Since the phrases differ
#      in length, they CHANGE THE FILE SIZE: comparing by hash or by size fails
#      every time, and a careless reader would conclude the numeric state
#      diverged.
#   4. Absolute path of the output directory passed via -d.
#
# There is NO numeric tolerance here. Any difference surviving this
# normalization is a real divergence and must block the stage (Principle I).
#
# Usage:  normalize-output.sh <file> <run-directory>
# Output: normalized content on stdout.

set -euo pipefail

input_file="${1:?usage: normalize-output.sh <file> <run-directory>}"
run_dir="${2:?usage: normalize-output.sh <file> <run-directory>}"

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "$script_dir/.." && pwd)"
phrase_list="${MARLIM_PHRASE_LIST:-$script_dir/decorative-phrases.txt}"

# The phrase list is derived from src/include/SisProd.h by extract-phrases.py,
# and regenerated here when missing or older than the header, so it can never
# drift out of sync with an edit to the phrases.
if [[ ! -s "$phrase_list" || "$project_root/src/include/SisProd.h" -nt "$phrase_list" ]]; then
    python3 "$script_dir/extract-phrases.py" \
        "$project_root/src/include/SisProd.h" > "$phrase_list"
fi

# Phrase removal is by CONTENT, line by line -- never by position.
#
# A positional rule ("drop whatever sits between asterisk rulers") was tried and
# REJECTED: a long transient run writes a LogEvento.dat with 159 rulers, an odd
# number, and they delimit legitimate result sections. Dropping between them
# would remove real numeric output, i.e. it could MASK a divergence, which is
# precisely the opposite of what this harness is for.
awk -v phrase_file="$phrase_list" '
    BEGIN {
        while ((getline line < phrase_file) > 0) {
            gsub(/^[ \t]+|[ \t]+$/, "", line)
            if (line != "") decorative[line] = 1
        }
    }
    {
        stripped = $0
        gsub(/^[ \t]+|[ \t]+$/, "", stripped)
        if (stripped in decorative) next     # drawn phrase: drop
        print
    }
' "$input_file" | sed -E \
    -e "s#${run_dir}#<DIR>#g" \
    -e 's#(/[A-Za-z0-9_.+-]+)+/([^/[:space:]"]+\.(dat|log|snp|txt|json))#<DIR>/\2#g' \
    -e 's/datahora = [0-9]{1,2}\/[0-9]{1,2}\/[0-9]{4} [0-9]{1,2}:[0-9]{1,2}:[0-9]{1,2}/datahora = <TS>/g' \
    -e 's/simulation date and time [0-9]{1,2}\/[0-9]{1,2}\/[0-9]{4} [0-9:]+/simulation date and time <TS>/g' \
    -e 's/data e hora da simulacao [0-9]{1,2}\/[0-9]{1,2}\/[0-9]{4} [0-9:]+/data e hora da simulacao <TS>/g' \
    -e 's/"data"[[:space:]]*:[[:space:]]*"[^"]*"/"data": "<TS>"/g' \
    -e 's/[0-9]{1,2}\/[0-9]{1,2}\/[0-9]{4} hora: [0-9]{1,2}:[0-9]{1,2}/<TS>/g' \
    -e 's/DURACAO[[:space:]]+[0-9]+[[:space:]]+segundos/DURACAO <DUR> segundos/g'
