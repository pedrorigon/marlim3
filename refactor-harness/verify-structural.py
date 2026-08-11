#!/usr/bin/env python3
"""Verification layer L0: prove a moved function body is token-for-token identical.

Why this layer exists
---------------------
Coverage measurement in stage 0 (T014) found that the demo corpus exercises
29.3% of the executable lines in src/core/SisProd.cpp, and that 81 of 153
functions -- 53% -- are never executed at all. Among them are the five largest
functions after RenovaMassPerm: renovaFracMol2 (571 uncovered lines),
RenovaMassPermComp (499), renovaFracMol (395), RenovaMassPermCompRev (281) and
RenovaMassPermRev (244).

For those, running the gates proves nothing. L2 compares outputs of code that
never ran; L3 covers even less. Moving 571 lines and watching six green gates
teaches you exactly nothing about whether the move was faithful.

L0 does not prove numerical equivalence -- nothing that skips execution can. It
proves the move was LITERAL, which is precisely what the extraction tasks
require: "corpo transferido token a token; nenhuma expressao em ponto flutuante
reescrita". Weaker than L2 where L2 works; far stronger than nothing where it
does not.

What is normalized away
-----------------------
Whitespace, line breaks and comments. Everything else -- every operator, every
literal, every parenthesis, the order of every operand -- must match exactly.
Reassociating (a+b)+c into a+(b+c), replacing x/y with x*(1.0/y), or
pre-computing 5*M_PI/180. all change the token stream and are caught.

Renaming a local variable also changes the token stream, and is reported. That
is intentional: FR-039 requires renaming to be a step separate from moving, so a
move that also renames is a process violation even when it is numerically
harmless. Use --allow-renames when verifying the rename step itself.

Usage:
    verify-structural.py --baseline <file> --current <file> --function <name>
    verify-structural.py --baseline <file> --current <file> --all

Exit code: 0 when every compared body matches; 1 otherwise.
"""
from __future__ import annotations

import argparse
import re
import sys

# C++ token: identifier, number (including hex float and exponent), operator or
# punctuation. Numbers are matched before identifiers so 1e-5 stays one token.
TOKEN = re.compile(r"""
      (?P<number>\.?\d[\w.]*(?:[eEpP][+-]?\d+)?[fFuUlL]*)
    | (?P<name>[A-Za-z_]\w*)
    | (?P<string>"(?:[^"\\]|\\.)*"|'(?:[^'\\]|\\.)*')
    | (?P<op>::|->|\+\+|--|<<=|>>=|<<|>>|<=|>=|==|!=|&&|\|\||[-+*/%=<>!&|^~?:;,.(){}\[\]])
""", re.VERBOSE)

COMMENT = re.compile(r"//[^\n]*|/\*.*?\*/", re.DOTALL)

SIGNATURE = re.compile(r"^[A-Za-z_][\w:<>,~&*=!\[\]\s]*\(")
NOT_A_DEFINITION = re.compile(
    r"^\s*(#|//|/\*|\*|\}|using|namespace|struct|class|typedef|template|extern|"
    r"return|else|public|private|protected|enum|union)"
)
NAME = re.compile(r"(operator\s*[^\s(]+)\s*\(|([\w~]+)\s*\(")
MAX_SIGNATURE_LINES = 12


def tokenize(text: str) -> list[str]:
    """Token stream with comments and all whitespace removed."""
    return [m.group(0) for m in TOKEN.finditer(COMMENT.sub(" ", text))]


def function_bodies(path: str) -> dict[str, str]:
    """Map function name -> source text, for every file-scope definition."""
    with open(path, encoding="utf-8", errors="replace") as handle:
        lines = handle.readlines()

    bodies: dict[str, str] = {}
    index, total = 0, len(lines)

    while index < total:
        line = lines[index]
        if SIGNATURE.match(line) and not NOT_A_DEFINITION.match(line):
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
                    name = (match.group(1) or match.group(2)) if match else None
                    if name:
                        # Body only: from the opening brace to the closing one,
                        # so a changed signature (SProd:: prefix dropped on
                        # extraction) does not register as a difference.
                        body = "".join(lines[probe:closing + 1])
                        body = body[body.index("{"):]
                        bodies[name] = body
                    index = closing + 1
                    continue
        index += 1

    return bodies


def strip_names(tokens: list[str]) -> list[str]:
    """Replace identifiers with a placeholder, keeping structure and literals."""
    return ["<id>" if re.fullmatch(r"[A-Za-z_]\w*", t) else t for t in tokens]


def first_difference(expected: list[str], actual: list[str]) -> str:
    for position, (left, right) in enumerate(zip(expected, actual)):
        if left != right:
            start = max(0, position - 6)
            return (f"token {position}: expected {left!r}, found {right!r}\n"
                    f"      baseline: ...{' '.join(expected[start:position + 6])}...\n"
                    f"      current : ...{' '.join(actual[start:position + 6])}...")
    shorter, longer = ("current", "baseline") if len(actual) < len(expected) else ("baseline", "current")
    return f"{shorter} ends after {min(len(expected), len(actual))} tokens; {longer} continues"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--baseline", required=True)
    parser.add_argument("--current", required=True)
    parser.add_argument("--function", action="append", default=[])
    parser.add_argument("--all", action="store_true")
    parser.add_argument("--allow-renames", action="store_true",
                        help="ignore identifier names; compare structure and literals only")
    args = parser.parse_args()

    baseline = function_bodies(args.baseline)
    current = function_bodies(args.current)

    if args.all:
        targets = sorted(set(baseline) & set(current))
        missing = sorted(set(baseline) - set(current))
    else:
        targets = args.function
        missing = [f for f in targets if f not in current]
        targets = [f for f in targets if f in current]

    if not targets and not missing:
        print("no functions to compare", file=sys.stderr)
        return 2

    failures = 0
    for name in missing:
        print(f"MISSING  {name}: present in baseline, absent from current", file=sys.stderr)
        failures += 1

    for name in targets:
        if name not in baseline:
            print(f"NEW      {name}: absent from baseline (introduced here)")
            continue

        expected = tokenize(baseline[name])
        actual = tokenize(current[name])
        if args.allow_renames:
            expected, actual = strip_names(expected), strip_names(actual)

        if expected == actual:
            print(f"OK       {name} ({len(expected)} tokens)")
        else:
            print(f"DIFFERS  {name}", file=sys.stderr)
            print("      " + first_difference(expected, actual), file=sys.stderr)
            failures += 1

    print(f"\ncompared {len(targets)} function(s), {failures} failure(s)")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
