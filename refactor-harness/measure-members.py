#!/usr/bin/env python3
"""Measure which SProd data members and methods a set of member functions uses.

A module's state struct is DERIVED from this measurement, not written from a
reading of the code: stage 7 wrote parameter names from the task list instead
of the definitions and paid for it. Stage 8 measured, but the measurement
lived in a scratch directory and is not reproducible; this is that
measurement, kept.

How
---
Members and methods come from the SProd class body in SisProd.h, split into
top-level statements (inline bodies skipped). A statement with a `(` before
any `=` declares a method; otherwise its declarators are data members, with
their type (`double *`, `int [10]`, `vector<int>`).

In each function body, an identifier counts as a member use when it names a
member and is not reached through `.`, `->` or `::` -- `celula[i].pres` uses
`celula`, `arq.ncelp` uses `arq`, `(*vg1dSP).x` uses `vg1dSP`. A name the body
declares itself shadows the member and does not count; run -Wshadow on the
host first, because this tool trusts that there is no such shadow and says so.

Usage:
    measure-members.py <function> [<function> ...]

Prints, per member: type, and which functions use it; then per method: which
functions call it. Exit 0.
"""
import io
import re
import sys
import pathlib

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import braces

ROOT = pathlib.Path(__file__).resolve().parent.parent
SRC = ROOT / "src/core/SisProd.cpp"
HDR = ROOT / "src/include/SisProd.h"
KEYWORDS = {"if", "for", "while", "switch", "return", "sizeof", "new", "delete",
            "else", "do", "case", "const", "static", "int", "double", "bool",
            "char", "void", "float", "long", "unsigned", "true", "false", "this"}


def strip_comments(s):
    return re.sub(r'//[^\n]*', ' ', re.sub(r'/\*.*?\*/', ' ', s, flags=re.S))


def class_statements(hdr):
    start = hdr.index("class SProd {")
    o = hdr.index("{", start)
    e = braces.match(hdr, o)
    body = strip_comments(hdr[o + 1:e])
    body = re.sub(r'^\s*(public|private|protected)\s*:', ' ', body, flags=re.M)
    stmts, depth, cur, i = [], 0, [], 0
    while i < len(body):
        c = body[i]
        if c == "{":
            depth += 1
            if depth == 1:
                stmts.append("".join(cur) + "{}")   # an inline body ends the statement
                cur = []
        elif c == "}":
            depth -= 1
        elif depth == 0:
            if c == ";":
                stmts.append("".join(cur))
                cur = []
            else:
                cur.append(c)
        i += 1
    return [" ".join(s.split()) for s in stmts if s.strip()]


def members_and_methods(hdr):
    data, methods = {}, set()
    for s in class_statements(hdr):
        if s.startswith(("friend ", "typedef ", "using ", "enum ", "struct ", "class ")):
            continue
        paren, eq = s.find("("), s.find("=")
        if paren >= 0 and (eq < 0 or paren < eq):
            m = re.search(r'(operator\s*\S+?|~?[A-Za-z_]\w*)\s*\($', s[:paren + 1])
            if m:
                methods.add(m.group(1))
            continue
        # data members: type followed by one or more declarators
        m = re.match(r'^((?:static\s+|const\s+|mutable\s+)*[A-Za-z_][\w:]*(?:\s*<.*?>)?)\s+(.*)$', s)
        if not m:
            continue
        base, rest = m.group(1), m.group(2)
        depth, pieces, cur = 0, [], ""
        for ch in rest:
            if ch in "(<[{":
                depth += 1
            elif ch in ")>]}":
                depth -= 1
            if ch == "," and depth == 0:
                pieces.append(cur)
                cur = ""
            else:
                cur += ch
        pieces.append(cur)
        for p in pieces:
            p = p.split("=")[0].strip()
            d = re.match(r'^([*&\s]*)([A-Za-z_]\w*)\s*((?:\[[^\]]*\])*)$', p)
            if d:
                stars = d.group(1).replace(" ", "")
                ty = base + (" " + stars if stars else "") + (" " + d.group(3) if d.group(3) else "")
                data[d.group(2)] = ty
    return data, methods


def function_body(src, name):
    m = re.search(r'^[A-Za-z][^\n;{}]*?\bSProd::%s\(' % re.escape(name), src, re.M)
    if m is None:
        raise SystemExit("not found: SProd::%s" % name)
    o = braces.first_brace_after(src, m.start())
    return strip_comments(src[m.start():braces.match(src, o) + 1])


def uses(body, data, methods):
    used_data, used_methods = set(), set()
    for m in re.finditer(r'[A-Za-z_]\w*', body):
        name = m.group(0)
        before = body[:m.start()].rstrip()
        if before.endswith((".", "->", "::")):
            continue
        after = body[m.end():].lstrip()
        if name in methods and after.startswith("("):
            used_methods.add(name)
        elif name in data:
            used_data.add(name)
    return used_data, used_methods


def main():
    names = sys.argv[1:]
    if not names:
        print(__doc__)
        return 2
    hdr = io.open(HDR, encoding="utf-8", errors="surrogateescape").read()
    src = io.open(SRC, encoding="utf-8", errors="surrogateescape").read()
    data, methods = members_and_methods(hdr)
    by_member, by_method = {}, {}
    for n in names:
        d, meth = uses(function_body(src, n), data, methods)
        for x in d:
            by_member.setdefault(x, []).append(n)
        for x in meth:
            if x != n:
                by_method.setdefault(x, []).append(n)
    print("# SProd: %d data members, %d methods; %d functions measured"
          % (len(data), len(methods), len(names)))
    print("# members used: %d" % len(by_member))
    for x in sorted(by_member, key=lambda k: (-len(by_member[k]), k)):
        print("M %-28s %-26s %s" % (x, data[x], ",".join(by_member[x])))
    print("# methods called: %d" % len(by_method))
    for x in sorted(by_method):
        print("F %-28s %s" % (x, ",".join(by_method[x])))
    return 0


if __name__ == "__main__":
    sys.exit(main())
