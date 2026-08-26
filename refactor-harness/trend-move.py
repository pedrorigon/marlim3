#!/usr/bin/env python3
"""Move the eight trend writers out of SisProd.cpp, and prove the move literal.

The move is expressed as a substitution table rather than as retyped code. That
matters for two reasons.

Retranscribing 521 lines by hand introduces exactly the class of error the stage
exists to avoid, and the four cross-section writers are never executed by the
demo corpus (evidencia/trend-diff.md, D4-02) -- a transcription slip in them
would survive every gate.

Because the transformation is a table, it inverts. `--check` applies the inverse
to the extracted module and compares the token stream against the pristine
baseline, exactly. Anything the table does not describe -- a reassociated
expression, a dropped parenthesis, a changed literal, a renamed local -- fails
to invert and is reported with the offending token.

Usage:
    trend-move.py extract <sisprod.cpp> <output.cpp>
    trend-move.py check   <baseline-sisprod.cpp> <extracted.cpp>
"""
from __future__ import annotations

import re
import subprocess
import sys

# Old SProd method -> new free function in namespace trendoutput.
FUNCTIONS = {
    "ImprimeTrendPCab": "writeProductionTrendHeader",
    "ImprimeTrendP": "writeProductionTrendRows",
    "ImprimeTrendGCab": "writeServiceTrendHeader",
    "ImprimeTrendG": "writeServiceTrendRows",
    "ImprimeTrendTransPCab": "writeProductionCrossSectionTrendHeader",
    "ImprimeTrendTransP": "writeProductionCrossSectionTrendRows",
    "ImprimeTrendTransGCab": "writeServiceCrossSectionTrendHeader",
    "ImprimeTrendTransG": "writeServiceCrossSectionTrendRows",
}

# SProd member -> TrendState field. Whole-identifier matches only, so ntrend
# never matches inside ntrendB and arq never matches inside arqRelatorioPerfis.
# pathPrefixoArqSaida and arqRelatorioPerfis are file-scope globals declared in
# Leitura.h, not members, and are deliberately absent from this table.
MEMBERS = {
    "arq": "state.inputData",
    "indTramo": "state.branchIndex",
    "kimpT": "state.printPassCount",
    "vg1dSP": "state.globals",
    "MatTrendP": "state.productionBuffer",
    "ntrend": "state.productionCount",
    "ntrendB": "state.productionCountBase",
    "MatTrendG": "state.serviceBuffer",
    "ntrendg": "state.serviceCount",
    "ntrendgB": "state.serviceCountBase",
    "MatTrendTransP": "state.productionCrossSectionBuffer",
    "MatTrendTransG": "state.serviceCrossSectionBuffer",
    "ntrendtrans": "state.crossSectionCount",
    "ntrendtransB": "state.crossSectionCountBase",
}

# The translation lambda captures `this` to reach arq; with the state passed in,
# it captures the state instead. Applied before the member table.
LAMBDA_OLD = "[this](const char *pt, const char *en) {"
LAMBDA_NEW = "[&state](const char *pt, const char *en) {"
THIS_ARQ_OLD = "this->arq.idiomaSaida"
THIS_ARQ_NEW = "state.inputData.idiomaSaida"

MEMBER_RE = re.compile(r"\b(" + "|".join(sorted(MEMBERS, key=len, reverse=True)) + r")\b")
FIELD_RE = re.compile(r"\bstate\.(" + "|".join(
    sorted((v.split(".", 1)[1] for v in MEMBERS.values()), key=len, reverse=True)) + r")\b")
FIELD_TO_MEMBER = {v.split(".", 1)[1]: k for k, v in MEMBERS.items()}

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
    """Locate each SProd::ImprimeTrend* definition by signature and brace depth."""
    lines = source.split("\n")
    found: dict[str, tuple[int, int, str]] = {}
    for start, line in enumerate(lines):
        match = re.match(r"void SProd::(ImprimeTrend\w*)\(", line)
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


def forward(name: str, body: str) -> str:
    """Apply the move transformation to one function body."""
    signature = re.match(r"void SProd::\w+\((.*?)\) \{", body).group(1)
    new_signature = "const TrendState &state" + (", " + signature if signature else "")
    body = re.sub(r"^void SProd::\w+\(.*?\) \{",
                  f"void {FUNCTIONS[name]}({new_signature}) {{", body, count=1)
    body = body.replace(LAMBDA_OLD, LAMBDA_NEW).replace(THIS_ARQ_OLD, THIS_ARQ_NEW)
    return MEMBER_RE.sub(lambda m: MEMBERS[m.group(1)], body)


def inverse(new_name: str, body: str) -> str:
    """Apply the exact inverse, so the result must equal the baseline body."""
    old_name = {v: k for k, v in FUNCTIONS.items()}[new_name]
    # Mirror image of forward(): the lambda capture is restored BEFORE the field
    # table runs, or state.inputData.idiomaSaida collapses to arq.idiomaSaida and
    # the `this->` can never be put back.
    body = body.replace(LAMBDA_NEW, LAMBDA_OLD).replace(THIS_ARQ_NEW, THIS_ARQ_OLD)
    body = FIELD_RE.sub(lambda m: FIELD_TO_MEMBER[m.group(1)], body)
    signature = re.match(rf"void {new_name}\(const TrendState &state,?\s*(.*?)\) \{{",
                         body).group(1)
    return re.sub(rf"^void {new_name}\(.*?\) \{{",
                  f"void SProd::{old_name}({signature}) {{", body, count=1)


def carve_new(source: str) -> dict[str, str]:
    lines = source.split("\n")
    found: dict[str, str] = {}
    for start, line in enumerate(lines):
        match = re.match(r"void (write\w+)\(const TrendState", line)
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


HEADER = '''/*
 * SisProdTrendOutput.cpp
 *
 * Trend-file output for the production and service lines, moved out of
 * SisProd.cpp. See SisProdTrendOutput.h for why the writers take their state
 * through TrendState instead of reading it from SProd.
 *
 * The bodies below were moved token for token. The five pre-existing anomalies
 * they carry -- the service cross-section writer reading the production
 * counters, the TENDTRANSG file named TENDTRANSP on the branch path, the
 * asymmetric round() in the TENDP name, the cross-section copy that ignores the
 * base offset, and the three stray endl -- are preserved deliberately, not
 * overlooked. They are catalogued in
 * specs/001-refatoracao-sisprod/evidencia/trend-diff.md.
 */
#include "SisProdTrendOutput.h"

#include "Leitura.h"
#include "Matriz.h"
#include "OutputI18n.h"
#include "variaveisGlobais1D.h"

#include <fstream>
#include <math.h>
#include <sstream>
#include <string>

namespace trendoutput {

'''

FOOTER = "\n} // namespace trendoutput\n"


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
        open(second, "w", encoding="utf-8").write(HEADER + "\n\n".join(parts) + FOOTER)
        print(f"extracted {len(parts)} function(s) into {second}")
        for name in FUNCTIONS:
            start, end, _ = bodies[name]
            print(f"  {name:24} {start + 1}-{end + 1}")
        return 0

    if mode == "check":
        module = open(second, encoding="utf-8").read()
        # This check answers one question: was the move literal? It can only
        # answer it while the module still has the shape the move produced.
        # Later steps of the stage regrouped the state, folded the writers onto
        # shared skeletons, renamed the locals and corrected four defects --
        # after any of those, the substitution table no longer inverts, and
        # emitting a token diff would report a failure that is not one. Say so
        # instead: a verifier that cries wolf is worse than no verifier.
        # A field that existed only in the shape the move produced. The scalar
        # fields survived the later regrouping, so testing for those would let
        # the check run on a module it can no longer invert.
        if MEMBERS["MatTrendP"] not in module:
            print("This module has been restructured past the literal move, so\n"
                  "the substitution table no longer inverts. The move itself was\n"
                  "verified at the commit that performed it:\n"
                  "    git show 9167524:src/core/SisProdTrendOutput.cpp > /tmp/moved.cpp\n"
                  "    trend-move.py check <baseline SisProd.cpp> /tmp/moved.cpp\n"
                  "Everything after that commit is verified by\n"
                  "refactor-harness/verify-trend-writers.sh, which compares the\n"
                  "bytes the writers produce rather than the tokens they are\n"
                  "made of.", file=sys.stderr)
            return 2
        baseline = carve(open(first, encoding="utf-8").read())
        current = carve_new(module)
        failures = 0
        for old, new in FUNCTIONS.items():
            if new not in current:
                print(f"FAIL {new}: absent from {second}")
                failures += 1
                continue
            expected = tokenize(baseline[old][2])
            actual = tokenize(inverse(new, current[new]))
            if expected == actual:
                print(f"OK   {new:42} <- {old} ({len(expected)} tokens)")
            else:
                failures += 1
                print(f"FAIL {new:42} <- {old}")
                for pos, (left, right) in enumerate(zip(expected, actual)):
                    if left != right:
                        low = max(0, pos - 6)
                        print(f"     token {pos}: baseline {left!r}, current {right!r}")
                        print(f"     baseline: ...{' '.join(expected[low:pos + 6])}...")
                        print(f"     current : ...{' '.join(actual[low:pos + 6])}...")
                        break
                else:
                    print(f"     length differs: baseline {len(expected)}, current {len(actual)}")
        print(f"\ncompared {len(FUNCTIONS)} function(s), {failures} failure(s)")
        return 1 if failures else 0

    print(f"unknown mode {mode!r}", file=sys.stderr)
    return 2


if __name__ == "__main__":
    sys.exit(main())
