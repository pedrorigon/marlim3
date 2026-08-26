#!/usr/bin/env python3
"""Prove the multMarcha dispatch table selects what the original chain selected.

C-C6 of contracts/busca-raiz.md requires the new selection to produce exactly
the same march for every combination of the five selectors. No gate can check
that. multMarcha is only reachable from zriddr, and the thirteen call sites pass
just six of the nine rows; the other three -- the injection row and the two
gas-line rows -- are selected by no model in the corpus. A row wired to the
wrong march there leaves L2 and L3 green.

Nor is the check something to do by reading. The original is a conditional chain
nested four deep whose leaves are eight different methods, and one of its
branches is dead on purpose (A2-01, evidencia/anomalias.md): both arms select
marchaProdPresPres2, and a "tidied" version that quietly routes one of them to
marchaProdPresPres3 would look like an improvement and change results.

So both sides are carved out of source and compiled together:

  * the baseline side, from SProd::multMarcha in the pristine commit, with each
    `return marchaX(chute);` rewritten to `return SteadyMarch::marchaX;` and the
    three member reads turned into parameters -- a substitution table, printed
    with --show so it can be read rather than trusted;

  * the current side, selectSteadyMarch, taken verbatim from SisProd.cpp.

Then every combination of the selectors is swept and compared. The choke opening
includes values on both sides of the 1e-15 threshold, the threshold itself, and
a NaN, because `NaN > 1e-15` is false and that is a path.

Usage:
    verify-dispatch.py [--baseline-commit 0f3b64f] [--show]
Exit code: 0 when every combination agrees.
"""
from __future__ import annotations

import argparse
import re
import subprocess
import sys
import tempfile
from pathlib import Path

MARCHES = (
    "marchaInjPerm1",
    "marchaGasPerm2",
    "marchaGasPerm3",
    "marchaProdPerm1",
    "marchaProdPerm1Rev",
    "marchaProdPerm2",
    "marchaProdPresPres1",
    "marchaProdPresPres1Rev",
    "marchaProdPresPres2",
)

# How the baseline chain becomes a selector. Each rule is a whole-token rewrite;
# nothing about the structure of the conditionals is touched, so the branch
# nesting, the comparison operators and the 1e-15 literal all survive intact.
BASELINE_RULES = (
    (r"^double SProd::multMarcha\(double chute, int prod, int tipoCC\) \{",
     "SteadyMarch baselineSelect(int injectorWell, int prod, int tipoCC,\n"
     "                           int reverseMarch, double productionChokeOpening) {"),
    (r"\barq\.pocinjec\b", "injectorWell"),
    (r"\brevPerm\b", "reverseMarch"),
    (r"\barq\.chokep\.abertura\[0\]", "productionChokeOpening"),
    (r"\breturn (marcha\w+)\(chute\);", r"return SteadyMarch::\1;"),
)

# The selector values swept. Chosen to cross every comparison in the chain:
# each selector is tested equal to 0, equal to 1, negative and above 1, because
# the chain tests `== 0` and `== 1` and treats everything else as a third case.
SWEEP = {
    "injectorWell": (-1, 0, 1, 2),
    "prod": (-1, 0, 1, 2, 3),
    "tipoCC": (-1, 0, 1, 2),
    "reverseMarch": (-1, 0, 1, 2),
    "productionChokeOpening": ("-1.0", "0.0", "1e-16", "1e-15", "1.0000000000000002e-15",
                               "0.5", "std::numeric_limits<double>::quiet_NaN()"),
}

DRIVER = """
#include <cstdio>
#include <limits>

enum class SteadyMarch {{
{enumerators}
}};

static const char *nameOf(SteadyMarch march) {{
    switch (march) {{
{names}
    }}
    return "<none>";
}}

{baseline}

{current}

int main() {{
    const int injectorWells[] = {{{injectorWell}}};
    const int prods[] = {{{prod}}};
    const int tipoCCs[] = {{{tipoCC}}};
    const int reverseMarches[] = {{{reverseMarch}}};
    const double openings[] = {{{productionChokeOpening}}};

    long long compared = 0, differing = 0;
    for (int injectorWell : injectorWells)
    for (int prod : prods)
    for (int tipoCC : tipoCCs)
    for (int reverseMarch : reverseMarches)
    for (double opening : openings) {{
        SteadyMarch expected = baselineSelect(injectorWell, prod, tipoCC, reverseMarch, opening);
        SteadyMarch actual = selectSteadyMarch(injectorWell, prod, tipoCC, reverseMarch, opening);
        ++compared;
        if (expected != actual) {{
            ++differing;
            if (differing <= 10)
                std::printf("DIFF pocinjec=%d prod=%d tipoCC=%d revPerm=%d abertura=%.17g"
                            "  baseline=%s current=%s\\n",
                            injectorWell, prod, tipoCC, reverseMarch, opening,
                            nameOf(expected), nameOf(actual));
        }}
    }}
    std::printf("compared %lld combination(s), %lld differing\\n", compared, differing);
    return differing == 0 ? 0 : 1;
}}
"""


def carve(source: str, pattern: str) -> str:
    """The definition whose first line matches `pattern`, to its closing brace."""
    lines = source.split("\n")
    for start, line in enumerate(lines):
        if not re.match(pattern, line):
            continue
        depth, opened = 0, False
        for probe in range(start, len(lines)):
            depth += lines[probe].count("{") - lines[probe].count("}")
            opened = opened or "{" in lines[probe]
            if opened and depth == 0:
                return "\n".join(lines[start:probe + 1])
    raise SystemExit(f"no definition matching {pattern!r}")


def to_selector(body: str) -> str:
    """Rewrite the baseline chain into a selector, one rule at a time."""
    for pattern, replacement in BASELINE_RULES:
        body, count = re.subn(pattern, replacement, body, flags=re.MULTILINE)
        if count == 0:
            raise SystemExit(f"baseline rule matched nothing: {pattern!r}\n"
                             "The chain changed shape; the rewrite would compare "
                             "something other than the original.")
    leftover = re.search(r"\b(arq|revPerm|chute)\b", body)
    if leftover:
        raise SystemExit(f"baseline selector still refers to {leftover.group(1)!r}; "
                         "a member read was not turned into a parameter")
    return body


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline-commit", default="0f3b64f")
    parser.add_argument("--current",
                        help="file to take selectSteadyMarch from; defaults to the "
                             "working tree's SisProd.cpp. Used to calibrate this "
                             "checker against deliberately corrupted copies.")
    parser.add_argument("--show", action="store_true",
                        help="print the rewritten baseline selector and exit")
    args = parser.parse_args()

    root = Path(__file__).resolve().parent.parent
    pristine = subprocess.run(
        ["git", "-C", str(root), "show", f"{args.baseline_commit}:src/core/SisProd.cpp"],
        capture_output=True, text=True, check=True).stdout
    current_path = Path(args.current) if args.current else root / "src/core/SisProd.cpp"
    current = current_path.read_text(encoding="utf-8", errors="replace")

    baseline = to_selector(carve(pristine, r"double SProd::multMarcha\("))
    selector = carve(current, r"SteadyMarch selectSteadyMarch\(")

    if args.show:
        print(baseline)
        return 0

    driver = DRIVER.format(
        enumerators="".join(f"    {m},\n" for m in MARCHES).rstrip(),
        names="".join(f'    case SteadyMarch::{m}: return "{m}";\n' for m in MARCHES).rstrip(),
        baseline=baseline,
        current=selector,
        **{key: ", ".join(str(v) for v in values) for key, values in SWEEP.items()},
    )

    with tempfile.TemporaryDirectory(prefix="marlim3-dispatch-") as work:
        source = Path(work) / "dispatch.cpp"
        binary = Path(work) / "dispatch"
        source.write_text(driver, encoding="utf-8")
        build = subprocess.run(
            ["g++", "-std=c++20", "-Wall", "-O2", "-ffp-contract=off",
             str(source), "-o", str(binary)],
            capture_output=True, text=True)
        if build.returncode != 0:
            print(build.stderr, file=sys.stderr)
            return 2
        if build.stderr.strip():
            print("the generated driver produced warnings:", file=sys.stderr)
            print(build.stderr, file=sys.stderr)
            return 2
        run = subprocess.run([str(binary)], capture_output=True, text=True)
        print(run.stdout, end="")
        if run.returncode == 0:
            print("DISPATCH IDENTICAL -- every combination selects the same march")
        else:
            print("DISPATCH DIVERGED", file=sys.stderr)
        return run.returncode


if __name__ == "__main__":
    sys.exit(main())
