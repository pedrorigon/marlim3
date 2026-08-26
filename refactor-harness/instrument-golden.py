#!/usr/bin/env python3
"""Instrument a COPY of the baseline to capture L1 golden values.

Layer L1 of the verification harness compares an extracted function against
input/output pairs tabulated from the baseline. This script produces those
pairs by patching a throwaway copy of the tree -- the working tree is never
touched.

Values are written in C99 hexadecimal float format (%a), which round-trips
exactly. Decimal formatting would silently drop bits and turn a bit-for-bit
comparison into an approximate one, defeating the purpose.

Two instrumentation strategies, chosen per function:

  * Drift-flux correlations have a single exit (zero `return` statements) and
    return their results through `double &` parameters. A logging call is
    injected just before the closing brace.

  * Root-finding solvers have multiple `return` statements, so injecting before
    the closing brace would miss most calls. They are renamed to
    `<name>__golden_impl` and a thin wrapper under the original name logs the
    arguments and the returned value.

Usage: instrument-golden.py <path-to-copy> <output-directory>
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

# The drift-flux correlations, located by signature rather than by line number.
#
# They started out in src/core/SisProd.cpp as SProd members and moved to
# src/core/DriftFluxClosure.cpp as free functions in driftflux::correlations.
# Absolute line numbers were correct for exactly one commit; searching for the
# signature works before and after the extraction, which is also the rule the
# refactoring tasks impose on every move.
CORRELATIONS = (
    "BhagwatGhajar",
    "BhagwatGhajarMod",
    "Choi",
    "HibikiIshii",
    "FrancaLahey",
)

# Where the correlations live once extracted. Empty/absent means the tree is
# still pristine and they are in SisProd.cpp.
CLOSURE_SOURCE = "src/core/DriftFluxClosure.cpp"
CLOSURE_NAMESPACE = "driftflux::correlations"

# The captured columns, in order. The two tuples are the same twelve quantities
# under the names each file uses: SisProd.cpp still has the inherited ones, and
# the extracted module was renamed to English when it moved. The order must stay
# aligned so the recorded tables keep comparing column for column.
CORRELATION_ARGS = ("rhol", "rhog", "tensup", "alf", "reymix", "reymixL",
                    "ug1", "ul1", "dia", "rug", "tet", "correcHor")
CORRELATION_ARGS_EXTRACTED = (
    "liquidDensity", "gasDensity", "surfaceTension", "voidFraction",
    "mixtureReynolds", "liquidReynolds", "gasFlowRate", "liquidFlowRate",
    "diameter", "roughness", "inclinationAngle", "horizontalCorrection")
assert len(CORRELATION_ARGS) == len(CORRELATION_ARGS_EXTRACTED)

# Solvers are wrapped rather than injected into, because they have several
# `return` statements and a probe before the closing brace would miss most
# calls.
#
# Which of these exist depends on the tree. In the pristine baseline all three
# are SProd members. From stage 2 on, zbrent and falsacorda are free templates
# in namespace rootfinding -- they have no caller anywhere and never executed,
# so nothing is lost by not capturing them -- and SProd::zriddr survives as the
# binding between the production domain and rootfinding::zriddr, with the same
# signature it always had. Instrument whichever are present; require zriddr.
SOLVERS = {
    "zbrent":     ("double", "double x1, double x2, int prod, int tipoCC, double tol, double epsn, int maxit",
                   "x1, x2, prod, tipoCC, tol, epsn, maxit",
                   [("x1", "%a"), ("x2", "%a"), ("prod", "%d"), ("tipoCC", "%d"),
                    ("tol", "%a"), ("epsn", "%a"), ("maxit", "%d")]),
    "falsacorda": ("double", "double a, double b, int prod, int tipoCC",
                   "a, b, prod, tipoCC",
                   [("a", "%a"), ("b", "%a"), ("prod", "%d"), ("tipoCC", "%d")]),
    "zriddr":     ("double", "double x1, double x2, int prod, int tipoCC",
                   "x1, x2, prod, tipoCC",
                   [("x1", "%a"), ("x2", "%a"), ("prod", "%d"), ("tipoCC", "%d")]),
}
REQUIRED_SOLVERS = ("zriddr",)

# The objective the solvers evaluate. Counting calls to it turns the capture
# from "same root" into "same trajectory", which is the difference that matters:
# two different convergence paths can land on the same value, and comparing only
# the returned root would not tell them apart.
#
# The counter wraps SProd::multMarcha rather than anything inside a solver, and
# that is what makes it survive the extraction. Before it, zriddr calls
# multMarcha directly; after it, the objective lambda does. Both go through this
# wrapper, so the recorded column means the same thing on both sides.
COUNTED = ("multMarcha", "double", "double chute, int prod, int tipoCC", "chute, prod, tipoCC")

LOGGER = '''
// ---- L1 golden capture (temporary instrumentation, never committed) ----
//
// The drift-flux correlations run per cell, per timestep, per Newton iteration.
// An unsampled capture wrote 950 MB for a single function in seconds, so this
// records one call in every MARLIM_GOLDEN_STRIDE (default 500) and stops at
// MARLIM_GOLDEN_CAP records per function (default 20000).
//
// Striding rather than truncating matters: taking the first N calls would
// sample only the first model's first timesteps, and the point of L1 is to
// cover the range of inputs the corpus actually produces.
#include <cstdio>
#include <cstdlib>
#include <mutex>
namespace golden_capture {

struct Channel {
    const char *name;
    std::FILE *file;
    long long seen;
    long long written;
};

inline long long env_long(const char *key, long long fallback) {
    const char *raw = std::getenv(key);
    return raw ? std::atoll(raw) : fallback;
}

// Returns the sink when this call should be recorded, nullptr otherwise.
inline std::FILE *sink(const char *name) {
    static std::mutex guard;
    std::lock_guard<std::mutex> held(guard);
    static Channel channels[32];
    static int count = 0;
    static const long long stride = env_long("MARLIM_GOLDEN_STRIDE", 500);
    static const long long cap = env_long("MARLIM_GOLDEN_CAP", 20000);

    Channel *channel = nullptr;
    for (int i = 0; i < count; ++i)
        if (channels[i].name == name) { channel = &channels[i]; break; }

    if (!channel) {
        const char *dir = std::getenv("MARLIM_GOLDEN_DIR");
        if (!dir || count >= 32) return nullptr;
        char path[1024];
        std::snprintf(path, sizeof(path), "%s/%s.txt", dir, name);
        channels[count] = Channel{name, std::fopen(path, "a"), 0, 0};
        channel = &channels[count];
        ++count;
    }

    if (!channel->file) return nullptr;
    long long index = channel->seen++;
    if (channel->written >= cap) return nullptr;
    if (stride > 1 && (index % stride) != 0) return nullptr;
    channel->written++;
    return channel->file;
}

// Objective evaluations so far. A solver records the delta across its own call,
// which is its iteration count in the only unit that is comparable across the
// extraction.
inline long long &evaluations() {
    static long long total = 0;
    return total;
}

}
// ---- end L1 golden capture ----
'''


def locate_correlations(lines: list[str]) -> dict[str, int]:
    """Map correlation name -> index of its closing-brace line.

    Matches `void <name>(` with or without an `SProd::` qualifier, then scans
    for the closing brace in column zero. Every correlation has a single exit,
    so that brace is the one logging must precede.
    """
    found: dict[str, int] = {}
    for name in CORRELATIONS:
        opening = re.compile(rf"^\s*void\s+(?:SProd::)?{name}\s*\(")
        for index, line in enumerate(lines):
            if not opening.match(line):
                continue
            closing = index
            while closing < len(lines) and not lines[closing].startswith("}"):
                closing += 1
            if closing < len(lines):
                found[name] = closing
            break
    return found


def instrument_correlations(source: str, extracted: bool = False) -> tuple[str, set[str]]:
    """Inject the capture call before each correlation's closing brace.

    Returns the patched source and the set of correlations actually found, so
    the caller can tell which file they live in.
    """
    args = CORRELATION_ARGS_EXTRACTED if extracted else CORRELATION_ARGS
    lines = source.splitlines(keepends=True)
    located = locate_correlations(lines)
    # Patch from the bottom up so earlier line numbers stay valid.
    for name, closing in sorted(located.items(), key=lambda item: -item[1]):
        fields = " ".join("%a" for _ in args) + " %a %a"
        values = ", ".join(args) + ", c0, ud"
        call = (f'    {{ std::FILE *gf = golden_capture::sink("{name}");\n'
                f'      if (gf) std::fprintf(gf, "{fields}\\n", {values}); }}\n')
        lines.insert(closing, call)
    return "".join(lines), set(located)


def instrument_solvers(source: str) -> tuple[str, set[str]]:
    """Wrap each solver present, plus the objective counter.

    Returns the patched source and the set of solvers actually wrapped, so the
    caller can insist on the ones that must be there and stay quiet about the
    ones a later stage legitimately moved out of the class.
    """
    name, ret, params, args = COUNTED
    if f"{ret} SProd::{name}({params})" not in source:
        raise SystemExit(f"{name} not found as `{ret} SProd::{name}({params})` -- "
                         "the signature moved and the counter would silently "
                         "record nothing")
    source = source.replace(f"{ret} SProd::{name}({params})",
                            f"{ret} SProd::{name}__golden_impl({params})", 1)
    source += f'''
{ret} SProd::{name}({params}) {{
    ++golden_capture::evaluations();
    return {name}__golden_impl({args});
}}
'''

    wrapped: set[str] = set()
    for name, (ret, params, args, logged) in SOLVERS.items():
        definition = f"{ret} SProd::{name}({params})"
        if definition not in source:
            continue
        wrapped.add(name)
        # Rename the definition; declarations in the header get the same
        # treatment separately.
        source = source.replace(definition,
                                f"{ret} SProd::{name}__golden_impl({params})", 1)

        fields = " ".join(fmt for _, fmt in logged) + " %a %lld"
        values = ", ".join(arg for arg, _ in logged)
        wrapper = f'''
{ret} SProd::{name}({params}) {{
    long long before = golden_capture::evaluations();
    {ret} result = {name}__golden_impl({args});
    {{ std::FILE *gf = golden_capture::sink("{name}");
      if (gf) std::fprintf(gf, "{fields}\\n", {values}, result,
                           golden_capture::evaluations() - before); }}
    return result;
}}
'''
        source += wrapper
    return source, wrapped


def instrument_header(header: str, wrapped: set[str]) -> str:
    targets = [(n, SOLVERS[n][0], SOLVERS[n][1]) for n in wrapped]
    targets.append((COUNTED[0], COUNTED[1], COUNTED[2]))
    for name, ret, params in targets:
        # Declare the renamed implementation next to the original declaration.
        pattern = re.compile(rf"(\n\s*{ret}\s+{name}\s*\([^;]*\);)")
        match = pattern.search(header)
        if not match:
            raise SystemExit(f"no declaration of {name} to sit beside in SisProd.h")
        header = header.replace(
            match.group(1),
            match.group(1) + f"\n    {ret} {name}__golden_impl({params});", 1)
    return header


SWEEP = '''
// ---- L1 synthetic sweep (temporary instrumentation, never committed) ----
//
// Five of the eight target functions are NEVER called by the demo corpus:
// Choi, BhagwatGhajarMod and FrancaLahey are unreachable because Leitura.cpp
// pins arq.CorreDisper to 1 and no model overrides it; SProd::zbrent has no
// call site anywhere, and SProd::falsacorda is reached only from inside it.
//
// For those, a capture driven by the corpus produces nothing, so L2 would stay
// green even if the extraction corrupted them. This sweep calls the pure
// correlations directly over a grid of physically plausible inputs, giving L1
// something to compare against that does not depend on the corpus.
//
// Triggered by MARLIM_GOLDEN_SWEEP=<file>; the process writes the table and
// exits without running a simulation.
void SProd::goldenSweep(const char *path) {
    std::FILE *out = std::fopen(path, "w");
    if (!out) return;
__CORRELATION_SCOPE__

    const double liquidDensities[]  = {700.0, 850.0, 1000.0};
    const double gasDensities[]     = {5.0, 50.0, 200.0};
    const double surfaceTensions[]  = {0.005, 0.02, 0.06};
    const double voidFractions[]    = {0.01, 0.25, 0.5, 0.75, 0.99};
    const double mixtureReynolds[]  = {1.0e2, 1.0e4, 1.0e6};
    const double inclinations[]     = {-1.4, -0.2, 0.0, 0.2, 1.4};

    for (double rhol : liquidDensities)
    for (double rhog : gasDensities)
    for (double tensup : surfaceTensions)
    for (double alf : voidFractions)
    for (double reymix : mixtureReynolds)
    for (double tet : inclinations) {
        const double reymixL = reymix * 0.5;
        const double ug1 = 0.3, ul1 = 1.7, dia = 0.15, rug = 4.5e-5;
        const double correcHor = 1.0;
        double c0 = 0.0, ud = 0.0;

        #define GOLDEN_SWEEP_CASE(fn)                                          \\
            c0 = 0.0; ud = 0.0;                                                \\
            fn(rhol, rhog, tensup, alf, reymix, reymixL, ug1, ul1, dia, rug,   \\
               tet, c0, ud, correcHor);                                        \\
            std::fprintf(out, "%s %a %a %a %a %a %a %a %a %a %a %a %a %a %a\\n",\\
                         #fn, rhol, rhog, tensup, alf, reymix, reymixL, ug1,   \\
                         ul1, dia, rug, tet, correcHor, c0, ud);

        GOLDEN_SWEEP_CASE(BhagwatGhajar)
        GOLDEN_SWEEP_CASE(BhagwatGhajarMod)
        GOLDEN_SWEEP_CASE(Choi)
        GOLDEN_SWEEP_CASE(HibikiIshii)
        GOLDEN_SWEEP_CASE(FrancaLahey)
        #undef GOLDEN_SWEEP_CASE
    }

    std::fclose(out);
}
// ---- end L1 synthetic sweep ----
'''


TRIGGER = """
    // L1 synthetic sweep trigger (temporary instrumentation)
    if (const char *sweepPath = std::getenv("MARLIM_GOLDEN_SWEEP")) {
        goldenSweep(sweepPath);
        std::exit(0);
    }
"""


def add_sweep(source: str, header: str, extracted: bool) -> tuple[str, str]:
    # Once the correlations move into driftflux::correlations they are no
    # longer members, so the sweep needs them in scope. It must stay a using
    # directive rather than qualified calls: the table's first column comes
    # from stringifying the macro argument, and `driftflux::correlations::Choi`
    # there would diverge from the recorded golden for a non-numerical reason.
    # Drops the whole placeholder line when pristine, so instrumenting an
    # un-extracted tree still reproduces the recorded golden byte for byte.
    scope = (f"    using namespace {CLOSURE_NAMESPACE};\n" if extracted else "")
    source += SWEEP.replace("__CORRELATION_SCOPE__\n", scope)
    if extracted and '#include "DriftFluxClosure.h"' not in source:
        source = source.replace('#include "SisProd.h"',
                                '#include "SisProd.h"\n#include "DriftFluxClosure.h"', 1)
    # Fire the sweep from the first line of the main constructor body, before
    # any input parsing, so it needs no model file.
    ctor = source.find("SProd::SProd(")
    if ctor >= 0:
        body = source.find("{", ctor)
        newline = source.find("\n", body)
        source = source[:newline + 1] + TRIGGER + source[newline + 1:]
    # goldenSweep has to be declared inside the class. The old anchor was the
    # zbrent declaration, which stage 2 removed -- and a .replace() that matches
    # nothing returns the string unchanged, so the sweep would simply never be
    # declared. Try several anchors and refuse to continue if none is there,
    # rather than producing an instrumented tree that quietly captures nothing.
    anchors = (
        "    double multMarcha(double chute, int prod, int tipoCC);",
        "    double zriddr(double x1, double x2, int prod, int tipoCC);",
        "    void renovaTemp();",
    )
    for anchor in anchors:
        if anchor in header:
            header = header.replace(
                anchor, "    void goldenSweep(const char *path);\n" + anchor, 1)
            break
    else:
        raise SystemExit("no anchor left in SisProd.h to declare goldenSweep beside; "
                         "tried:\n  " + "\n  ".join(anchors))
    return source, header


def add_logger(source: str) -> str:
    """Place the capture logger after the last #include, so it sees <cstdio>."""
    last_include = source.rfind("\n#include")
    end_of_line = source.find("\n", last_include + 1)
    return source[:end_of_line + 1] + LOGGER + source[end_of_line + 1:]


def main() -> int:
    if len(sys.argv) != 3:
        print(__doc__)
        return 2

    tree = Path(sys.argv[1])
    cpp = tree / "src/core/SisProd.cpp"
    hpp = tree / "src/include/SisProd.h"
    closure = tree / CLOSURE_SOURCE

    source = cpp.read_text(encoding="utf-8", errors="replace")
    if "golden_capture" in source:
        print("already instrumented; aborting", file=sys.stderr)
        return 1

    # Before the extraction the correlations are in SisProd.cpp; after it they
    # are in DriftFluxClosure.cpp. Instrument wherever they actually are.
    source, in_sisprod = instrument_correlations(source)
    source, wrapped = instrument_solvers(source)
    missing_solvers = [n for n in REQUIRED_SOLVERS if n not in wrapped]
    if missing_solvers:
        print(f"solvers not found as SProd members: {', '.join(missing_solvers)}",
              file=sys.stderr)
        return 1

    in_closure: set[str] = set()
    if closure.exists():
        closure_text = closure.read_text(encoding="utf-8", errors="replace")
        closure_text, in_closure = instrument_correlations(closure_text, extracted=True)
        if in_closure:
            closure_text = add_logger(closure_text)
            closure.write_text(closure_text, encoding="utf-8")

    located = in_sisprod | in_closure
    missing = [name for name in CORRELATIONS if name not in located]
    if missing:
        print(f"correlations not found: {', '.join(missing)}", file=sys.stderr)
        return 1

    header_text = hpp.read_text(encoding="utf-8", errors="replace")
    source, header_text = add_sweep(source, header_text, extracted=bool(in_closure))

    header_text = instrument_header(header_text, wrapped)
    if "goldenSweep" not in header_text:
        print("goldenSweep was never declared -- refusing to write", file=sys.stderr)
        return 1

    cpp.write_text(add_logger(source), encoding="utf-8")
    hpp.write_text(header_text, encoding="utf-8")

    print(f"instrumented {cpp}")
    if in_closure:
        print(f"instrumented {closure}")
    print(f"  correlations: {', '.join(sorted(located))}"
          f"{' (extracted)' if in_closure else ''}")
    print(f"  solvers     : {', '.join(sorted(wrapped))}")
    print(f"  counter     : {COUNTED[0]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
