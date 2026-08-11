#!/usr/bin/env python3
"""Extract the 32 decorative phrases embedded in src/include/SisProd.h.

The engine picks one phrase per run (`srand(time(NULL)); rand() % 16` in
Num4Main.cpp:596 and 8084) and writes it to stdout, to the .snp header and to
the tail of LogEvento.dat. Because the phrases have different lengths, they
change the file size -- they are the engine's only source of textual
non-determinism.

This script emits the literal phrase list so that normalize-output.py removes
exactly those lines and nothing else.

Why not "drop whatever sits between asterisk rulers": a long transient run
produces a LogEvento.dat with 159 rulers, an odd number, and they delimit
legitimate result sections. A positional rule would delete real numeric output,
i.e. it could MASK a divergence -- the opposite of what the harness is for.

Usage: extract-phrases.py [header-path] > decorative-phrases.txt
"""
from __future__ import annotations

import re
import sys

DEFAULT_HEADER = "src/include/SisProd.h"
PHRASE_ARRAYS = ("saidaTextoSis", "saidaSubTextoSis")
EXPECTED_PHRASES = 32


def extract(header_path: str) -> list[str]:
    source = open(header_path, encoding="utf-8", errors="replace").read()
    phrases: list[str] = []

    for array_name in PHRASE_ARRAYS:
        start = source.find(f"{array_name}[16]")
        if start < 0:
            print(f"warning: array {array_name} not found in {header_path}",
                  file=sys.stderr)
            continue
        end = source.find("};", start)
        block = source[start:end]
        # Double-quoted literals; this header has no escaped quotes inside them.
        phrases.extend(match.group(1) for match in re.finditer(r'"([^"]*)"', block))

    return phrases


def main() -> int:
    header_path = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_HEADER
    phrases = extract(header_path)

    if len(phrases) != EXPECTED_PHRASES:
        print(f"warning: expected {EXPECTED_PHRASES} phrases, extracted "
              f"{len(phrases)}", file=sys.stderr)

    for phrase in phrases:
        # The engine pads phrases to a fixed column width, so match on the
        # stripped content to stay robust against padding changes.
        print(phrase.strip())

    return 0


if __name__ == "__main__":
    sys.exit(main())
