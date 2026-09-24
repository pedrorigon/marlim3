#!/usr/bin/env python3
"""Prove that the named functions never write a parameter they take by value.

A cut helper receives each crossing local by reference when the analyser saw
the range write it, and by value otherwise. If the analyser misses a write --
`stream >> x`, `f(x)` with a reference parameter, `++x` hidden in a macro --
the helper modifies its own copy and the caller reads a stale value
afterwards. The compiler says nothing, and L2 says something only if the corpus
happens to take that path. This is the silent half of the T130 defect class.

The check lets the compiler answer instead: copy the source into a scratch
tree, add `const` to every by-value parameter of the named functions -- top-level
const is not part of the function type, so declarations elsewhere still match --
and compile with -fsyntax-only. Any write, direct or through a reference, is
now an error that names the parameter.

Usage:
    check-byvalue-writes.py <source.cpp> <function> [<function> ...]

A name matches a definition `... <Class>::<function>(` or `... <function>(` at
column 0. Exit code: 0 when no by-value parameter is written; 1 otherwise;
2 when a named function was not found (a check that looked at nothing is not a
pass).
"""
import io
import os
import re
import shutil
import subprocess
import sys
import tempfile
import pathlib

ROOT = pathlib.Path(__file__).resolve().parent.parent
FLAGS_MAKE = ROOT / "build/CMakeFiles/Marlim3.dir/flags.make"
SCALAR_OR_CLASS = re.compile(r'^[A-Za-z_][\w:<>, ]*\s+[A-Za-z_]\w*$')


def constify(signature_params):
    out, count = [], 0
    for p in signature_params.split(","):
        core = p.split("=")[0].strip()
        if (not core or "&" in core or "*" in core or "[" in core
                or core.startswith("const ") or " const " in core
                or not SCALAR_OR_CLASS.match(core)):
            out.append(p)
            continue
        lead = p[:len(p) - len(p.lstrip())]
        out.append(lead + "const " + p.lstrip())
        count += 1
    return ",".join(out), count


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        return 2
    source = pathlib.Path(sys.argv[1]).resolve()
    names = sys.argv[2:]
    text = io.open(source, encoding="utf-8", errors="surrogateescape").read()

    total, spans = 0, {}
    for name in names:
        m = re.search(r'^[A-Za-z_][^\n;{}]*?\b(?:\w+::)?%s\(([^()]*)\)\s*(?:const\s*)?\{' % re.escape(name),
                      text, re.M)
        if m is None:
            print("function not found: %s" % name, file=sys.stderr)
            return 2
        params, count = constify(m.group(1))
        text = text[:m.start(1)] + params + text[m.end(1):]
        spans[name] = count
        total += count

    flags = FLAGS_MAKE.read_text()
    includes = re.search(r'^CXX_INCLUDES = (.*)$', flags, re.M).group(1).split()
    scratch = pathlib.Path(tempfile.mkdtemp(prefix="byvalue-"))
    try:
        # Same directory layout, so relative includes resolve as in the build.
        rel = source.relative_to(ROOT)
        target = scratch / rel
        target.parent.mkdir(parents=True)
        io.open(target, "w", encoding="utf-8", errors="surrogateescape").write(text)
        includes = [i.replace(str(ROOT / "src/core"), str(target.parent)) for i in includes]
        cmd = ["g++", "-std=c++20", "-fopenmp", "-fsyntax-only", *includes,
               "-I", str(source.parent), str(target)]
        run = subprocess.run(cmd, capture_output=True, text=True)
    finally:
        shutil.rmtree(scratch, ignore_errors=True)

    errors = [l for l in run.stderr.splitlines() if "error:" in l]
    for name in names:
        print("%-44s %d parametro(s) por valor" % (name, spans[name]))
    if errors:
        print("WRITES TO BY-VALUE PARAMETERS (or another error):")
        for l in errors[:20]:
            print("  " + l.split("/src/", 1)[-1])
        return 1
    print("BY-VALUE CHECK PASSED -- %d parameter(s), none written" % total)
    return 0


if __name__ == "__main__":
    sys.exit(main())
