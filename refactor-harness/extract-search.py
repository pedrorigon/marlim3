#!/usr/bin/env python3
"""Lift a line range out of a function in SisProdSteadyStateSearch.cpp.

Written once instead of eight times. The searches all have the same shape --
estimate a guess, handle the march's sentinel returns, bracket the sign change,
solve -- so the same three operations keep recurring:

  * find the boundary at CHARACTER level, because `} else if (c) {` closes and
    reopens inside one line and a per-line counter walks past it (T085);
  * work out which locals cross the boundary RESPECTING SCOPE, because a
    declaration inside `if (chute < 0) { }` is not visible after that block and
    listing it invents a parameter (measured, and wrong the first time);
  * convert `return <sentinel>;` into an abort protocol when the range is not
    the function's tail, because a helper cannot return from its caller.

Usage:
    extract-search.py <function> <first> <last> <helper> [--tail] [--doc TEXT]

--tail says the range ends the function, so its returns stay returns and the
call site becomes `return helper(...)`. Without it, each return becomes
`{ abortValue = <v>; return true; }` and the helper answers "must it stop?".
"""
import argparse
import io
import re
import sys

sys.path.insert(0, "/tmp/claude-1000/-home-pedrorigon-Desktop-marlim-public-marlim3/6bb08355-c2ca-4c65-bcab-00ebb9119233/scratchpad")
import braces
import iface

PATH = "src/core/SisProdSteadyStateSearch.cpp"
RETURN = re.compile(r'(\n[ \t]*)return\s+([^;]+);')


def function_span(text, name):
    m = re.search(rf'^(?:\[\[nodiscard\]\]\s+)?(?:void|double|int|bool)\s+{name}\(', text, re.M)
    if m is None:
        raise SystemExit("function not found: " + name)
    open_brace = braces.first_brace_after(text, m.start())
    return m.start(), braces.match(text, open_brace)


def line_bounds(text, start_char, first, last):
    """Char offsets of lines `first`..`last`, counted from the function's line 1."""
    base = braces.line_of(text, start_char)          # 0-based line of the signature
    lines = text.split("\n")
    a = sum(len(lines[k]) + 1 for k in range(base + first - 1))
    b = sum(len(lines[k]) + 1 for k in range(base + last)) - 1
    return a, b


def dedent(block, extra="    "):
    rows = block.split("\n")
    pad = len(rows[0]) - len(rows[0].lstrip(" "))
    return "\n".join((extra + r[pad:]) if r.startswith(" " * pad)
                     else (extra + r.strip() if r.strip() else "") for r in rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("function")
    ap.add_argument("first", type=int)
    ap.add_argument("last", type=int)
    ap.add_argument("helper")
    ap.add_argument("--tail", action="store_true")
    ap.add_argument("--doc", default="")
    ap.add_argument("--returns", default="double")
    ap.add_argument("--body", action="store_true",
                    help="the range's first and last lines are the braces of an else or a "
                         "loop; lift what is between them and leave the brace line in place")
    args = ap.parse_args()

    text = io.open(PATH, encoding="utf-8", errors="surrogateescape").read()
    fa, fb = function_span(text, args.function)
    first, last = args.first, args.last
    if args.body:
        # `} else {` opens a block whose body is the next line onwards. Taking
        # the brace line too produced a helper opening with a stray else, which
        # the compiler caught once already in the mass-march work.
        first, last = args.first + 1, args.last - 1
    a, b = line_bounds(text, fa, first, last)
    block = text[a:b]

    scratch = ("/tmp/claude-1000/-home-pedrorigon-Desktop-marlim-public-marlim3/"
               "6bb08355-c2ca-4c65-bcab-00ebb9119233/scratchpad/_x.cpp")
    io.open(scratch, "w", encoding="utf-8", errors="surrogateescape").write(text[fa:fb + 1])
    crossing_in, crossing_out = iface.interface(scratch, first, last)
    if crossing_out:
        raise SystemExit("locals declared inside the range and read after it: %s\n"
                         "That is an out-parameter the tool will not invent."
                         % ", ".join(n for _, n in crossing_out))

    # Every local that crosses the boundary is passed. It is tempting to move a
    # scratch variable's declaration INTO the helper when the range writes it
    # and nothing reads it afterwards -- it shortens the signature and reads
    # better. That rule is wrong, and the compiler proved it: the range
    # `while (corrigechute == 1) { ... }` writes corrigechute and nothing reads
    # it after, but it also READS it on entry, and a moved declaration would
    # have re-initialised it. Deciding this properly needs to know whether the
    # first use is a write, which is flow analysis. A wider signature is the
    # cheaper mistake.
    params, moved = [], []
    for ty, name, written, after in crossing_in:
        if name in ("state", "abortValue"):
            # state is passed first by construction; abortValue is the abort
            # protocol's own out-parameter and is appended below. Cutting a
            # helper that already has one would otherwise declare it twice.
            continue
        params.append(f"{ty} {'&' if written else ''}{name}")

    body = dedent(block)

    if args.tail:
        signature = (f"{args.returns} {args.helper}(const SteadyStateSearchState &state, "
                     + ", ".join(params) + ") {")
        call = "    return %s(state%s%s);" % (
            args.helper, ", " if params else "",
            ", ".join(p.split()[-1].lstrip("&") for p in params))
    elif not re.search(r'(?<![.\w])return\b',
                       re.sub(r'/\*.*?\*/', '', re.sub(r'//[^\n]*', '', block), flags=re.S)):
        # No early exit in the range, so no abort protocol. Adding one anyway
        # would put a parameter and a return value in the signature that say
        # nothing, which is the opposite of what naming a block is for.
        signature = (f"void {args.helper}(const SteadyStateSearchState &state, "
                     + ", ".join(params) + ") {")
        call = "    %s(state%s%s);" % (
            args.helper, ", " if params else "",
            ", ".join(p.split()[-1].lstrip("&") for p in params))
    else:
        params.append("double &abortValue")
        signature = (f"bool {args.helper}(const SteadyStateSearchState &state, "
                     + ", ".join(params) + ") {")
        count = [0]

        def sub(m):
            count[0] += 1
            pad = m.group(1)
            return "%s{%s    abortValue = %s;%s    return true;%s}" % (
                pad, pad, m.group(2).strip(), pad, pad)

        body = RETURN.sub(sub, body)
        code = re.sub(r'/\*.*?\*/', '', re.sub(r'//[^\n]*', '', body), flags=re.S)
        left = [m for m in re.finditer(r'(?<![.\w])return\b(?!\s+true;)', code)]
        assert not left, "a return survived: %r" % [code[m.start() - 40:m.end() + 20] for m in left]
        body += "\n    return false;"
        call = ("    double abortValue;\n"
                "    if (%s(state%s%s, abortValue))\n        return abortValue;" % (
                    args.helper, ", " if len(params) > 1 else "",
                    ", ".join(p.split()[-1].lstrip("&") for p in params[:-1])))

    doc = args.doc or ("Part of %s." % args.function)
    helper = "".join("/// %s\n" % line if line else "///\n" for line in doc.split("\n"))
    helper += signature + "\n" + body + "\n}\n\n"

    text = text[:a] + call + text[b:]
    fa, _ = function_span(text, args.function)
    s = text.rfind("\n", 0, fa) + 1
    while s > 0:
        prev = text.rfind("\n", 0, s - 1) + 1
        if text[prev:s].lstrip().startswith("//"):
            s = prev
        else:
            break
    text = text[:s] + helper + text[s:]
    io.open(PATH, "w", encoding="utf-8", errors="surrogateescape").write(text)
    print("%-52s %3d linhas -> auxiliar; %d parametro(s), %d local(is) movido(s)" % (
        args.helper, block.count("\n") + 1, len(params), len(moved)))


if __name__ == "__main__":
    main()
