#!/usr/bin/env python3
"""Measure the length of every file-scope function definition (C++).

Instrument for SC-004 ("no function exceeds 200 lines").

Why not a one-line awk
----------------------
The naive detector

    awk '/^[A-Za-z_].*\\{[[:space:]]*$/{...} /^\\}/{...}'

requires the opening brace on the SAME line as the signature. In
src/core/SisProd.cpp that misses 10 definitions with multi-line parameter lists
-- among them `copiaSemJson` (283 lines), `operator=` (302) and the main
constructor (235). The naive detector reports a false pass on exactly the
functions that violate the criterion.

This script accumulates the signature for up to MAX_SIGNATURE_LINES lines until
it finds the opening brace, and only then looks for the closing brace at
column 0.

Usage:
    measure-functions.py [--limit N] [--all] FILE...

Output: one line per function above the limit, largest first.
Exit code 1 when any function exceeds the limit, so it works as a gate.
"""
from __future__ import annotations

import argparse
import re
import sys

# A line opening a file-scope definition: starts at column 0 and contains '('.
# The character class includes '=', '!' and '[]' so operator overloads match --
# `SProd &SProd::operator=(const SProd &sp) {` is 302 lines and would be missed
# by a class accepting only [\w:<>,~&*\s].
SIGNATURE = re.compile(r"^[A-Za-z_][\w:<>,~&*=!\[\]\s]*\(")

# Constructs that start at column 0 but are NOT function definitions.
NOT_A_DEFINITION = re.compile(
    r"^\s*(#|//|/\*|\*|\}|using|namespace|struct|class|typedef|template|extern|"
    r"return|else|public|private|protected|enum|union)"
)

# Function name. The first alternative covers operator overloads, whose
# identifier ends in a symbol rather than a word character.
NAME = re.compile(r"(operator\s*[^\s(]+)\s*\(|([\w~]+)\s*\(")

MAX_SIGNATURE_LINES = 12   # how far to look for the opening brace


def measure(path: str) -> list[tuple[int, str, int, int]]:
    """Return [(length, name, first_line, last_line)] for each definition."""
    with open(path, encoding="utf-8", errors="replace") as handle:
        lines = handle.readlines()

    functions: list[tuple[int, str, int, int]] = []
    index, total = 0, len(lines)

    while index < total:
        line = lines[index]
        if SIGNATURE.match(line) and not NOT_A_DEFINITION.match(line):
            # Look for the opening brace; a ';' before it means a declaration.
            probe, opened = index, False
            while probe < min(index + MAX_SIGNATURE_LINES, total):
                if "{" in lines[probe]:
                    opened = True
                    break
                if ";" in lines[probe]:
                    break
                probe += 1

            if opened:
                closing = probe + 1
                while closing < total and not lines[closing].startswith("}"):
                    closing += 1
                if closing < total:
                    match = NAME.search(line)
                    name = (match.group(1) or match.group(2)) if match else "?"
                    functions.append((closing - index + 1, name,
                                      index + 1, closing + 1))
                    index = closing + 1
                    continue
        index += 1

    return functions


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("files", nargs="+")
    parser.add_argument("--limit", type=int, default=200)
    parser.add_argument("--all", action="store_true",
                        help="list every definition, not only those above the limit")
    args = parser.parse_args()

    over_limit = 0
    for path in args.files:
        try:
            functions = measure(path)
        except FileNotFoundError:
            print(f"{path}: missing (not created yet)", file=sys.stderr)
            continue

        offenders = [f for f in functions if f[0] > args.limit]
        listed = functions if args.all else offenders
        listed.sort(key=lambda f: -f[0])
        over_limit += len(offenders)

        print(f"{path}: {len(functions)} definitions, "
              f"{len(offenders)} above {args.limit} lines")
        for length, name, first, last in listed:
            print(f"  {length:5d}  {name:30s} ({first}-{last})")

    print(f"\nTOTAL above {args.limit} lines: {over_limit}")
    return 1 if over_limit else 0


if __name__ == "__main__":
    sys.exit(main())
