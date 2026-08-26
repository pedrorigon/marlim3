#!/usr/bin/env python3
"""Move the root-finding solvers out of SisProd.cpp, and prove the move literal.

The move is expressed as a substitution table rather than as retyped code, for
the same reason the trend writers were: retranscribing arithmetic by hand
introduces exactly the class of error the stage exists to avoid.

Here the argument is sharper than usual. Measured by call site over the whole
tree, SProd::zbrent has NO caller anywhere and SProd::falsacorda is reachable
only from inside it. Neither ever executes. L2 and L3 stay green whatever this
move does to them, so a slip in 96 of the 244 lines would survive every
constitutional gate. Worse, once they become templates that nothing
instantiates, even the compiler stops checking them -- an uninstantiated
template is only parsed, not type-checked.

Because the transformation is a table, it inverts. `check` applies the inverse
to the extracted module and compares the token stream against the pristine
baseline, exactly. Anything the table does not describe -- a reassociated
expression, a dropped parenthesis, a changed literal, a renamed local -- fails
to invert and is reported with the offending token.

`check` also verifies the binding site. Three pieces of zriddr leave the solver
and land in SProd::zriddr: the two lines that derive minit from the input deck,
and the arithmetic of the convergence monitor. Those lines are still moved code,
so they are compared against the table too, not merely eyeballed.

Usage:
    solver-move.py extract <sisprod.cpp>
    solver-move.py check   <baseline-sisprod.cpp> <RootFindingSolvers.h>
                           <sisprod.cpp> <RootFindingSolvers.cpp>
"""
from __future__ import annotations

import re
import sys

# Old SProd method -> new free function in namespace rootfinding. The names are
# kept: zbrent and zriddr are the canonical Numerical Recipes names, and
# falsacorda is fixed by contracts/busca-raiz.md. Renaming a function is not
# authorised by any task in this stage -- T037r covers local variables.
SOLVERS = ("zbrent", "falsacorda", "zriddr")

# The signature of each solver, before and after. Written out in full rather
# than derived, so that a change to any parameter has to be made here, in the
# open, instead of being absorbed by a pattern.
SIGNATURES = {
    "zbrent": (
        "double SProd::zbrent(double x1, double x2, int prod, int tipoCC, "
        "double tol, double epsn, int maxit) {",
        "template <typename Objective>\n"
        "double zbrent(double x1, double x2, Objective &&objective, "
        "double tol, double epsn, int maxit) {",
    ),
    "falsacorda": (
        "double SProd::falsacorda(double a, double b, int prod, int tipoCC) {",
        "template <typename Objective>\n"
        "double falsacorda(double a, double b, Objective &&objective) {",
    ),
    "zriddr": (
        "double SProd::zriddr(double x1, double x2, int prod, int tipoCC) {",
        "template <typename Objective, typename Monitor>\n"
        "double zriddr(double x1, double x2, Objective &&objective, Monitor &&monitor,\n"
        "              int revPerm, int minit) {",
    ),
}

# Every call to the domain dispatcher becomes a call to the objective. The
# argument is always a bare identifier, in all 19 sites across the three
# solvers, which is what makes the pattern safe to invert.
CALL_OLD = re.compile(r"\bmultMarcha\((\w+), prod, tipoCC\)")
CALL_NEW = re.compile(r"\bobjective\((\w+)\)")

# zbrent's fallback passes the objective along instead of the domain pair.
FALLBACK_OLD = "falsacorda(x1, x2, prod, tipoCC)"
FALLBACK_NEW = "falsacorda(x1, x2, objective)"

# NumError lives in FerramentasNumericas.h, which declares `using namespace
# std;` at file scope and also declares global zbrent/zriddr templates that
# would collide with these. The header reaches it through a one-line wrapper
# instead, so the include graph and the global namespace stay untouched.
ERROR_OLD = "NumError("
ERROR_NEW = "reportIterationLimit("

# Applied AFTER the call table, so the pattern is written against the already
# substituted text. The convergence monitor is domain feedback -- it scales the
# residual by a member and writes another member -- so it moves to the binding
# site as a callable, and the solver keeps only the call.
MONITOR = tuple(
    (
        f"            double {var} = objective({arg}) / monitConvPermBase;\n"
        f"            if (prod != 0) {{\n"
        f"                monitConvPerm = fabs({var});\n"
        f"            }}",
        f"            double {var} = monitor(objective({arg}));",
    )
    for var, arg in (("fm", "xm"), ("fnew", "ans"))
)

# minit is derived from the input deck. `arq` has no business inside a generic
# root finder, so the two lines move to the binding site and minit arrives as a
# parameter. The body's uses of it are untouched.
MINIT_OLD = ("    int minit=0;\n"
             "    if(arq.acopColAnulPermForte == 1)minit=10;\n")

# The two sign helpers move without any transformation beyond losing the
# SProd:: qualifier, so they are compared body against body. SIGN is inline in
# the header because zriddr calls it three times per iteration; sign has no
# caller anywhere and lives out of line in the module .cpp.
HELPERS = {
    "SIGN": ("double SProd::SIGN(double a, double b) {",
             "inline double SIGN(double a, double b) {"),
    "sign": ("int SProd::sign(double var) {",
             "int sign(double var) {"),
}

# What the binding site must contain, verbatim, for the pieces above to add
# back up to the original. Compared by token stream, so indentation and line
# breaks are free but every operator, literal and operand order is not.
BINDING = """
double SProd::zriddr(double x1, double x2, int prod, int tipoCC) {
    int minit=0;
    if(arq.acopColAnulPermForte == 1)minit=10;
    return rootfinding::zriddr(
        x1, x2,
        [&](double guess) { return multMarcha(guess, prod, tipoCC); },
        [&](double residual) {
            double normalized = residual / monitConvPermBase;
            if (prod != 0) {
                monitConvPerm = fabs(normalized);
            }
            return normalized;
        },
        revPerm, minit);
}
"""

TOKEN = re.compile(r"""
      (?P<number>\.?\d[\w.]*(?:[eEpP][+-]?\d+)?[fFuUlL]*)
    | (?P<name>[A-Za-z_]\w*)
    | (?P<string>"(?:[^"\\]|\\.)*"|'(?:[^'\\]|\\.)*')
    | (?P<op>[^\s\w])
""", re.VERBOSE)
COMMENT = re.compile(r"//[^\n]*|/\*.*?\*/", re.DOTALL)


def tokenize(text: str) -> list[str]:
    return [m.group(0) for m in TOKEN.finditer(COMMENT.sub(" ", text))]


def carve(source: str, pattern: str) -> dict[str, str]:
    """Locate each definition by signature and brace depth."""
    lines = source.split("\n")
    found: dict[str, str] = {}
    for start, line in enumerate(lines):
        match = re.match(pattern, line)
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


def carve_old(source: str) -> dict[str, str]:
    return carve(source, r"double SProd::(zbrent|falsacorda|zriddr)\(")


def carve_new(source: str) -> dict[str, str]:
    """Same, but reaching back over the `template <...>` line.

    The template head is part of the signature the table rewrites, so a body
    carved without it cannot be inverted -- and the failure would look like a
    difference in the code rather than in the carving.
    """
    lines = source.split("\n")
    found: dict[str, str] = {}
    for start, line in enumerate(lines):
        match = re.match(r"double (zbrent|falsacorda|zriddr)\(double", line)
        if not match:
            continue
        head = start
        while head > 0 and lines[head - 1].startswith("template <"):
            head -= 1
        depth, opened, end = 0, False, start
        for probe in range(start, len(lines)):
            depth += lines[probe].count("{") - lines[probe].count("}")
            opened = opened or "{" in lines[probe]
            if opened and depth == 0:
                end = probe
                break
        found[match.group(1)] = "\n".join(lines[head:end + 1])
    return found


def body_at(source: str, signature: str) -> str:
    """The definition that starts at `signature`, up to its closing brace."""
    lines = source[source.index(signature):].split("\n")
    depth, opened = 0, False
    for probe, line in enumerate(lines):
        depth += line.count("{") - line.count("}")
        opened = opened or "{" in line
        if opened and depth == 0:
            return "\n".join(lines[:probe + 1])
    raise SystemExit(f"unterminated definition at {signature!r}")


def forward(name: str, body: str) -> str:
    """Apply the move transformation to one solver body."""
    old_signature, new_signature = SIGNATURES[name]
    if old_signature not in body:
        raise SystemExit(f"{name}: signature not found as written in SIGNATURES")
    body = body.replace(old_signature, new_signature, 1)
    body = body.replace(FALLBACK_OLD, FALLBACK_NEW)
    body = CALL_OLD.sub(lambda m: f"objective({m.group(1)})", body)
    body = body.replace(ERROR_OLD, ERROR_NEW)
    for old, new in MONITOR:
        body = body.replace(old, new)
    return body.replace(MINIT_OLD, "")


def inverse(name: str, body: str) -> str:
    """Apply the exact inverse, so the result must equal the baseline body.

    Mirror image of forward(): the monitor blocks are restored BEFORE the call
    table is inverted, because their replacement text was written against the
    already substituted form. Inverting in the wrong order leaves
    `monitor(multMarcha(xm, prod, tipoCC))`, which no rule can take apart.
    """
    old_signature, new_signature = SIGNATURES[name]
    if name == "zriddr":
        body = body.replace("    if ((fl > 0.0", MINIT_OLD + "    if ((fl > 0.0", 1)
    for old, new in MONITOR:
        body = body.replace(new, old)
    body = body.replace(ERROR_NEW, ERROR_OLD)
    body = CALL_NEW.sub(lambda m: f"multMarcha({m.group(1)}, prod, tipoCC)", body)
    body = body.replace(FALLBACK_NEW, FALLBACK_OLD)
    return body.replace(new_signature, old_signature, 1)


def report(name: str, expected: list[str], actual: list[str]) -> bool:
    if expected == actual:
        print(f"OK   {name:<12} ({len(expected)} tokens)")
        return True
    print(f"FAIL {name}")
    for position, (left, right) in enumerate(zip(expected, actual)):
        if left != right:
            start = max(0, position - 6)
            print(f"     token {position}: baseline {left!r}, current {right!r}")
            print(f"     baseline: ...{' '.join(expected[start:position + 6])}...")
            print(f"     current : ...{' '.join(actual[start:position + 6])}...")
            return False
    shorter = "current" if len(actual) < len(expected) else "baseline"
    print(f"     {shorter} ends after {min(len(expected), len(actual))} tokens")
    return False


def main() -> int:
    if len(sys.argv) < 2:
        print(__doc__)
        return 2
    mode = sys.argv[1]

    if mode == "extract":
        source = open(sys.argv[2], encoding="utf-8", errors="replace").read()
        bodies = carve_old(source)
        missing = [n for n in SOLVERS if n not in bodies]
        if missing:
            print(f"not found in {sys.argv[2]}: {', '.join(missing)}", file=sys.stderr)
            return 1
        for name in SOLVERS:
            print(forward(name, bodies[name]))
            print()
        return 0

    if mode == "check":
        baseline = open(sys.argv[2], encoding="utf-8", errors="replace").read()
        extracted = open(sys.argv[3], encoding="utf-8", errors="replace").read()
        binding = open(sys.argv[4], encoding="utf-8", errors="replace").read()
        module = open(sys.argv[5], encoding="utf-8", errors="replace").read()

        original = carve_old(baseline)
        current = carve_new(extracted)
        failures = 0

        for name in SOLVERS:
            if name not in original:
                print(f"FAIL {name}: absent from the baseline")
                failures += 1
                continue
            if name not in current:
                print(f"FAIL {name}: absent from {sys.argv[3]}")
                failures += 1
                continue
            expected = tokenize(original[name])
            actual = tokenize(inverse(name, current[name]))
            if not report(name, expected, actual):
                failures += 1

        # The binding site carries moved code too: the minit derivation and the
        # monitor arithmetic. Checking the solver alone would leave both
        # unverified, and both are inside the live solver's iteration.
        for name, (old_signature, new_signature) in HELPERS.items():
            where = extracted if name == "SIGN" else module
            if old_signature not in baseline or new_signature not in where:
                print(f"FAIL {name}: signature not found as written in HELPERS")
                failures += 1
                continue
            expected = tokenize(body_at(baseline, old_signature))
            actual = tokenize(body_at(where, new_signature)
                              .replace(new_signature, old_signature, 1))
            if not report(name, expected, actual):
                failures += 1

        # The binding site carries moved code too: the minit derivation and the
        # monitor arithmetic. Checking the solver alone would leave both
        # unverified, and both are inside the live solver's iteration.
        adapter = carve(binding, r"double SProd::(zriddr)\(")
        if "zriddr" not in adapter:
            print("FAIL binding: SProd::zriddr not found")
            failures += 1
        elif not report("binding", tokenize(BINDING), tokenize(adapter["zriddr"])):
            failures += 1

        total = len(SOLVERS) + len(HELPERS) + 1
        print(f"compared {total} body/bodies, {failures} failure(s)")
        return 1 if failures else 0

    print(__doc__)
    return 2


if __name__ == "__main__":
    sys.exit(main())
