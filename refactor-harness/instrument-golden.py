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

# name -> (first_line, closing_brace_line) in the pristine baseline.
CORRELATIONS = {
    "BhagwatGhajar":    (4905, 4988),
    "BhagwatGhajarMod": (4990, 5072),
    "Choi":             (5074, 5088),
    "HibikiIshii":      (5090, 5105),
    "FrancaLahey":      (5107, 5121),
}

CORRELATION_ARGS = ("rhol", "rhog", "tensup", "alf", "reymix", "reymixL",
                    "ug1", "ul1", "dia", "rug", "tet", "correcHor")

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

}
// ---- end L1 golden capture ----
'''


def instrument_correlations(source: str) -> str:
    lines = source.splitlines(keepends=True)
    # Patch from the bottom up so earlier line numbers stay valid.
    for name, (_, closing) in sorted(CORRELATIONS.items(),
                                     key=lambda item: -item[1][1]):
        fields = " ".join("%a" for _ in CORRELATION_ARGS) + " %a %a"
        values = ", ".join(CORRELATION_ARGS) + ", c0, ud"
        call = (f'    {{ std::FILE *gf = golden_capture::sink("{name}");\n'
                f'      if (gf) std::fprintf(gf, "{fields}\\n", {values}); }}\n')
        lines.insert(closing - 1, call)
    return "".join(lines)


def instrument_solvers(source: str) -> str:
    for name, (ret, params, args, logged) in SOLVERS.items():
        # Rename the definition; declarations in the header get the same
        # treatment separately.
        source = source.replace(f"{ret} SProd::{name}({params})",
                                f"{ret} SProd::{name}__golden_impl({params})", 1)

        fields = " ".join(fmt for _, fmt in logged) + " %a"
        values = ", ".join(arg for arg, _ in logged)
        wrapper = f'''
{ret} SProd::{name}({params}) {{
    {ret} result = {name}__golden_impl({args});
    {{ std::FILE *gf = golden_capture::sink("{name}");
      if (gf) std::fprintf(gf, "{fields}\\n", {values}, result); }}
    return result;
}}
'''
        source += wrapper
    return source


def instrument_header(header: str) -> str:
    for name, (ret, params, _, _) in SOLVERS.items():
        # Declare the renamed implementation next to the original declaration.
        pattern = re.compile(rf"(\n\s*{ret}\s+{name}\s*\([^;]*\);)")
        match = pattern.search(header)
        if match:
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


def add_sweep(source: str, header: str) -> tuple[str, str]:
    source += SWEEP
    # Fire the sweep from the first line of the main constructor body, before
    # any input parsing, so it needs no model file.
    marker = "void SProd::montasistema"
    ctor = source.find("SProd::SProd(")
    if ctor >= 0:
        body = source.find("{", ctor)
        newline = source.find("\n", body)
        source = source[:newline + 1] + TRIGGER + source[newline + 1:]
    header = header.replace("    double zbrent(double, double, int prod, int tipoCC,",
                            "    void goldenSweep(const char *path);\n"
                            "    double zbrent(double, double, int prod, int tipoCC,", 1)
    return source, header


def main() -> int:
    if len(sys.argv) != 3:
        print(__doc__)
        return 2

    tree = Path(sys.argv[1])
    cpp = tree / "src/core/SisProd.cpp"
    hpp = tree / "src/include/SisProd.h"

    source = cpp.read_text(encoding="utf-8", errors="replace")
    if "golden_capture" in source:
        print("already instrumented; aborting", file=sys.stderr)
        return 1

    source = instrument_correlations(source)
    source = instrument_solvers(source)

    header_text = hpp.read_text(encoding="utf-8", errors="replace")
    source, header_text = add_sweep(source, header_text)

    # The logger goes after the last #include so it sees <cstdio> etc.
    last_include = source.rfind("\n#include")
    end_of_line = source.find("\n", last_include + 1)
    source = source[:end_of_line + 1] + LOGGER + source[end_of_line + 1:]

    cpp.write_text(source, encoding="utf-8")
    hpp.write_text(instrument_header(header_text), encoding="utf-8")

    print(f"instrumented {cpp}")
    print(f"  correlations: {', '.join(sorted(CORRELATIONS))}")
    print(f"  solvers     : {', '.join(sorted(SOLVERS))}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
