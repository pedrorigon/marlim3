#!/usr/bin/env python3
"""Move the five CalcC0Ud variants out of SisProd.cpp, and prove the move literal.

The move is expressed as a substitution table rather than as retyped code, for
the reason that governs this whole stage: gcov measured that CalcC0UdBuf,
CalcC0UdIni and CalcC0UdIniBuf are never executed by any of the fourteen demo
models. 718 of the 1248 lines being moved are invisible to the artifact and
regression layers, so a transcription slip in them would survive every gate.

Because the transformation is a table, it inverts. `check` applies the inverse
to the extracted module and compares the token stream against the source the
move started from, exactly. Anything the table does not describe -- a
reassociated expression, a dropped parenthesis, a changed literal, a renamed
local -- fails to invert and is reported with the offending token.

The baseline for `check` is the commit the stage starts from, not the pristine
commit 0f3b64f. Stage 1 already requalified the correlation calls in these five
bodies and swapped arq.Corre* for the cached selectors, so they differ from
pristine by construction. That earlier transformation was re-verified on its
own before this stage began -- undoing it restores exact token equality with
0f3b64f across all five bodies, 13071 tokens -- and the chain of the two
inverses is what ties this module back to the original.

Usage:
    c0ud-move.py extract <sisprod.cpp> <output-fragment.cpp>
    c0ud-move.py check   <baseline-sisprod.cpp> <extracted.cpp>
"""
from __future__ import annotations

import re
import sys

# Old SProd method -> new free function in namespace driftflux::coefficient.
FUNCTIONS = {
    "CalcC0Ud": "instantaneous",
    "CalcC0UdBuf": "buffered",
    "CalcC0UdIni": "initialization",
    "CalcC0UdIniBuf": "bufferedInitialization",
    "CalcC0UdPerm": "steadyState",
}

# SProd member -> ClosureState field. Whole-identifier matches only, so alfE
# never matches inside alfPigE or alf0E, and ncel never inside any longer name.
MEMBERS = {
    "celula": "state.cells",
    "ncel": "state.lastCell",
    "vg1dSP": "state.globals",
    "arq": "state.input",
    "driftSelectors": "state.selectors",
    "tGSup": "state.gasSurfaceTemperature",
    "alfE": "state.inletVoidFraction",
    "betaE": "state.inletColumnFraction",
    "presE": "state.inletPressure",
    "tempE": "state.inletTemperature",
    "iterperm": "state.steadyIteration",
}
FIELD_TO_MEMBER = {v.split(".", 1)[1]: k for k, v in MEMBERS.items()}

MEMBER_RE = re.compile(r"\b(" + "|".join(sorted(MEMBERS, key=len, reverse=True)) + r")\b")
FIELD_RE = re.compile(r"\bstate\.(" + "|".join(sorted(FIELD_TO_MEMBER, key=len, reverse=True)) + r")\b")

TOKEN = re.compile(r"""
      (?P<number>\.?\d[\w.]*(?:[eEpP][+-]?\d+)?[fFuUlL]*)
    | (?P<name>[A-Za-z_]\w*)
    | (?P<string>"(?:[^"\\]|\\.)*"|'(?:[^'\\]|\\.)*')
    | (?P<op>[^\s\w])
""", re.VERBOSE)
COMMENT = re.compile(r"//[^\n]*|/\*.*?\*/", re.DOTALL)


def tokenize(text: str) -> list[str]:
    return [m.group(0) for m in TOKEN.finditer(COMMENT.sub(" ", text))]


def carve(source: str) -> dict[str, tuple[int, int, str]]:
    """Locate each SProd::CalcC0Ud* definition by signature and brace depth."""
    lines = source.split("\n")
    found: dict[str, tuple[int, int, str]] = {}
    for start, line in enumerate(lines):
        match = re.match(r"void SProd::(CalcC0Ud\w*)\(", line)
        if not match:
            continue
        depth, opened, end = 0, False, start
        for probe in range(start, len(lines)):
            depth += lines[probe].count("{") - lines[probe].count("}")
            opened = opened or "{" in lines[probe]
            if opened and depth == 0:
                end = probe
                break
        found[match.group(1)] = (start, end, "\n".join(lines[start:end + 1]))
    return found


def carve_new(source: str) -> dict[str, str]:
    lines = source.split("\n")
    found: dict[str, str] = {}
    names = "|".join(FUNCTIONS.values())
    for start, line in enumerate(lines):
        match = re.match(rf"void ({names})\(const ClosureState &state", line)
        if not match:
            continue
        depth, opened, end = 0, False, start
        for probe in range(start, len(lines)):
            depth += lines[probe].count("{") - lines[probe].count("}")
            opened = opened or "{" in lines[probe]
            if opened and depth == 0:
                end = probe
                break
        found[match.group(1)] = "\n".join(lines[start:end + 1])
    return found


def forward(name: str, body: str) -> str:
    signature = re.match(r"void SProd::\w+\((.*?)\) \{", body).group(1)
    body = re.sub(r"^void SProd::\w+\(.*?\) \{",
                  f"void {FUNCTIONS[name]}(const ClosureState &state, {signature}) {{",
                  body, count=1)
    return MEMBER_RE.sub(lambda m: MEMBERS[m.group(1)], body)


def inverse(new_name: str, body: str) -> str:
    """Apply the exact inverse, so the result must equal the baseline body."""
    old_name = {v: k for k, v in FUNCTIONS.items()}[new_name]
    body = FIELD_RE.sub(lambda m: FIELD_TO_MEMBER[m.group(1)], body)
    signature = re.match(rf"void {new_name}\(const ClosureState &state, (.*?)\) \{{",
                         body).group(1)
    return re.sub(rf"^void {new_name}\(.*?\) \{{",
                  f"void SProd::{old_name}({signature}) {{", body, count=1)


def main() -> int:
    if len(sys.argv) != 4:
        print(__doc__, file=sys.stderr)
        return 2
    mode, first, second = sys.argv[1], sys.argv[2], sys.argv[3]

    if mode == "extract":
        bodies = carve(open(first, encoding="utf-8").read())
        missing = [n for n in FUNCTIONS if n not in bodies]
        if missing:
            print(f"not found in {first}: {', '.join(missing)}", file=sys.stderr)
            return 2
        parts = [forward(name, bodies[name][2]) for name in FUNCTIONS]
        open(second, "w", encoding="utf-8").write("\n\n".join(parts) + "\n")
        print(f"extracted {len(parts)} function(s) into {second}")
        for name in FUNCTIONS:
            start, end, _ = bodies[name]
            print(f"  {name:16} {start + 1}-{end + 1}")
        return 0

    if mode == "check":
        baseline = carve(open(first, encoding="utf-8").read())
        current = carve_new(open(second, encoding="utf-8").read())
        failures = 0
        for old_name, new_name in FUNCTIONS.items():
            if new_name not in current:
                print(f"MISSING  {new_name} not found in {second}")
                failures += 1
                continue
            if old_name not in baseline:
                print(f"MISSING  {old_name} not found in {first}")
                failures += 1
                continue
            want = tokenize(baseline[old_name][2])
            got = tokenize(inverse(new_name, current[new_name]))
            if want == got:
                print(f"OK       {new_name:24} <- {old_name} ({len(want)} tokens)")
                continue
            failures += 1
            where = next((i for i, (a, b) in enumerate(zip(want, got)) if a != b),
                         min(len(want), len(got)))
            print(f"DIFFERS  {new_name:24} <- {old_name}")
            print(f"           first difference at token {where}")
            print(f"           baseline: {' '.join(want[max(0, where - 6):where + 6])}")
            print(f"           current : {' '.join(got[max(0, where - 6):where + 6])}")
        total = sum(len(tokenize(baseline[n][2])) for n in FUNCTIONS if n in baseline)
        print(f"compared {len(FUNCTIONS)} body/bodies, {total} tokens, {failures} failure(s)")
        return 1 if failures else 0

    print(__doc__, file=sys.stderr)
    return 2


if __name__ == "__main__":
    sys.exit(main())
