#!/usr/bin/env python3
"""Prove the multMarcha dispatch table selects what the original chain selected.

C-C6 of contracts/busca-raiz.md requires the new selection to produce exactly
the same march for every combination of the five selectors. No gate can check
that. multMarcha is only reachable from zriddr, and the corpus is far narrower
than the table: over the 360 zriddr calls it makes, only three (prod, tipoCC)
pairs ever occur -- (0,0) 346 times, (1,0) four times, (1,1) ten times. That
reaches at most FOUR of the nine rows. The injection row, the second gas-line
row and all three pressure-pressure rows are selected by no model, and a row
wired to the wrong march among those five leaves L2 and L3 green.

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

  * the current side, selectSteadyMarch, taken verbatim from
    SisProdSteadyStateSearch.cpp, where T097b moved it.

Then every combination of the selectors is swept and compared. The choke opening
is passed as the ARRAY, not the value, and the sweep includes a null one: the
original subscripts it only inside one branch, and Ler::copia_chokeSup leaves
chokep.abertura null when parserie is not positive. A selector that read it
eagerly would turn a conditional dereference into an unconditional one, which no
corpus model and therefore no gate would notice. The values include both sides of
the 1e-15 threshold, the threshold itself, and a NaN, because `NaN > 1e-15` is
false and that is a path.

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
     "                           int reverseMarch, const double *productionChokeOpening) {"),
    (r"\barq\.pocinjec\b", "injectorWell"),
    (r"\brevPerm\b", "reverseMarch"),
    (r"\barq\.chokep\.abertura", "productionChokeOpening"),
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

    const double openingValues[] = {{{productionChokeOpening}}};

    // The openings are passed as ARRAYS, and one of them is null. Ler::copia_chokeSup
    // leaves chokep.abertura null when parserie is not positive, and the original chain
    // subscripts it only inside one branch. A selector that reads it eagerly would
    // dereference null on every dispatch -- and no corpus model takes that path, so no
    // gate would have said anything. Here it crashes, which is the point.
    const double *openings[] = {{
        &openingValues[0], &openingValues[1], &openingValues[2], &openingValues[3],
        &openingValues[4], &openingValues[5], &openingValues[6], nullptr,
    }};

    long long compared = 0, differing = 0;
    for (int injectorWell : injectorWells)
    for (int prod : prods)
    for (int tipoCC : tipoCCs)
    for (int reverseMarch : reverseMarches)
    for (const double *opening : openings) {{
        SteadyMarch expected = baselineSelect(injectorWell, prod, tipoCC, reverseMarch, opening);
        SteadyMarch actual = selectSteadyMarch(injectorWell, prod, tipoCC, reverseMarch, opening);
        ++compared;
        if (expected != actual) {{
            ++differing;
            if (differing <= 10)
                std::printf("DIFF pocinjec=%d prod=%d tipoCC=%d revPerm=%d abertura=%s"
                            "  baseline=%s current=%s\\n",
                            injectorWell, prod, tipoCC, reverseMarch,
                            opening ? "value" : "NULL",
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


# What SProd::multMarcha must pass as the choke argument. The sweep below compiles
# only the selector, so an eager subscript at the CALL SITE -- passing
# arq.chokep.abertura[0] instead of arq.chokep.abertura -- would be invisible to
# it. That is the mistake this check exists for, because it is the one that was
# actually made: it turns the original's conditional dereference into an
# unconditional one, and Ler::copia_chokeSup leaves that pointer null when
# parserie is not positive.
CALL_SITE = ("selectSteadyMarch(state.march.input.pocinjec, prod, tipoCC, state.reverseSteady,\n"
             "                              state.march.input.chokep.abertura)")


def check_call_site(source: str) -> str | None:
    """None when dispatchMarch hands the selector the array; a message otherwise."""
    body = carve(source, r"double dispatchMarch\(")
    if CALL_SITE in body:
        return None
    if "input.chokep.abertura[0]" in body:
        return ("dispatchMarch subscripts the choke array before the call:\n"
                "      it passes arq.chokep.abertura[0], so the dereference happens on\n"
                "      EVERY dispatch. The original only subscripted it inside one\n"
                "      branch, and the pointer can be null. Pass the array.")
    return ("dispatchMarch does not call the selector as expected.\n"
            f"      expected: {CALL_SITE}")


def check_lazy_subscript(selector: str) -> str | None:
    """None when the selector subscripts the choke array only where it must.

    Checked structurally rather than by running with a null pointer. A null was
    tried first and did not work: the read is dead if its value is unused, and
    at -O2 the compiler deletes it, so the sweep reported no difference on a
    selector that dereferences null on every call. A checker whose negative case
    depends on the optimiser not eliding the fault is not a checker.
    """
    hits = [line.strip() for line in selector.split("\n")
            if "productionChokeOpening[" in line]
    if len(hits) != 1:
        return (f"the selector subscripts the choke array {len(hits)} time(s); "
                "expected exactly 1")
    guard = selector.find("if (tipoCC != 0) {")
    subscript = selector.find("productionChokeOpening[")
    if guard < 0 or subscript < guard:
        return ("the selector subscripts the choke array outside the branch that\n"
                "      needs it, so the dereference happens on every dispatch. The\n"
                "      original only subscripted it inside `if (tipoCC != 0)`, and\n"
                "      Ler::copia_chokeSup can leave that pointer null.")
    return None


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
                             "working tree's SisProdSteadyStateSearch.cpp. Used to calibrate this "
                             "checker against deliberately corrupted copies.")
    parser.add_argument("--show", action="store_true",
                        help="print the rewritten baseline selector and exit")
    args = parser.parse_args()

    root = Path(__file__).resolve().parent.parent
    pristine = subprocess.run(
        ["git", "-C", str(root), "show", f"{args.baseline_commit}:src/core/SisProd.cpp"],
        capture_output=True, text=True, check=True).stdout
    # T097b moved the table to the search module. The default follows it: a
    # checker left pointing at the file the code used to be in reports
    # "no definition matching ..." and keeps its exit code, which is how an
    # instrument stops testing without anyone noticing.
    current_path = (Path(args.current) if args.current
                    else root / "src/core/SisProdSteadyStateSearch.cpp")
    current = current_path.read_text(encoding="utf-8", errors="replace")

    baseline = to_selector(carve(pristine, r"double SProd::multMarcha\("))
    selector = carve(current, r"SteadyMarch selectSteadyMarch\(")

    if args.show:
        print(baseline)
        return 0

    complaint = check_call_site(current) or check_lazy_subscript(selector)
    if complaint:
        print(f"CHOKE READ WRONG -- {complaint}", file=sys.stderr)
        return 1

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
