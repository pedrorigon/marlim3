#!/usr/bin/env python3
"""List the locals that hide a member or a global inside named SProd functions.

A cut moves text; prove-cut.py shows the text is the same. What it cannot show
is that every name in the helper still BINDS to what it bound to in the host:
a local of the host that the cut failed to pass would, in the helper, silently
bind to a member or global of the same name. That can only happen for names
some local hides, so the list of hiding locals is what makes the question
finite.

This used to be a shell pipeline, and the pipeline was broken for two tasks:
GCC prints `src/core/SisProd.cpp:LINE:COL:` without a leading slash, the sed
that should strip the path did not match, and awk compared the file name with
line numbers -- every function came out clean (T101, found in T102). Hence a
tool, and a calibration it must pass before it reports anything.

Usage:
    check-shadowing.py <function> [<function> ...]

Exit code: 0 after printing, per function, each hiding local with its line,
or "none"; 2 if the calibration fails.
"""
import io
import re
import subprocess
import sys
import pathlib

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import braces

ROOT = pathlib.Path(__file__).resolve().parent.parent
SRC = ROOT / "src/core/SisProd.cpp"
WARNING = re.compile(r'SisProd\.cpp:(\d+):\d+: warning: declaration of .(\w+). shadows (a member|a global|a previous local|a parameter)')
# The calibration: renovaFracMol declares `double dt` at its outermost level,
# hiding SProd::dt; buildGasLiftLine hides nothing. If the parsing breaks, the
# first stops being found and the tool refuses to answer.
CALIBRATION = (("renovaFracMol", "dt"), ("buildGasLiftLine", None))


def ranges(text, names):
    out = {}
    for n in names:
        m = re.search(r'^[A-Za-z][^\n;{}]*?\bSProd::%s\(' % re.escape(n), text, re.M)
        if m is None:
            raise SystemExit("not found: SProd::%s" % n)
        e = braces.match(text, braces.first_brace_after(text, m.start()))
        out[n] = (text[:m.start()].count("\n") + 1, text[:e].count("\n") + 1)
    return out


def main():
    names = sys.argv[1:]
    if not names:
        print(__doc__)
        return 2
    flags = (ROOT / "build/CMakeFiles/Marlim3.dir/flags.make").read_text()
    includes = re.search(r'^CXX_INCLUDES = (.*)$', flags, re.M).group(1).split()
    run = subprocess.run(["g++", "-std=c++20", "-fopenmp", "-fsyntax-only", "-Wshadow",
                          *includes, str(SRC)], capture_output=True, text=True)
    hits = [(int(m.group(1)), m.group(2), m.group(3)) for m in WARNING.finditer(run.stderr)]
    text = io.open(SRC, encoding="utf-8", errors="surrogateescape").read()

    calib = ranges(text, [n for n, _ in CALIBRATION])
    for n, want in CALIBRATION:
        a, b = calib[n]
        got = {name for line, name, _ in hits if a <= line <= b}
        if (want is None and got) or (want is not None and want not in got):
            print("CALIBRATION FAILED on %s: found %s" % (n, sorted(got) or "nothing"))
            return 2

    for n, (a, b) in ranges(text, names).items():
        inside = [(line, name, kind) for line, name, kind in hits if a <= line <= b]
        if not inside:
            print("%-40s none" % n)
        for line, name, kind in inside:
            print("%-40s line %d: `%s` hides %s" % (n, line, name, kind[2:] if kind.startswith("a ") else kind))
    return 0


if __name__ == "__main__":
    sys.exit(main())
