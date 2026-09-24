#!/usr/bin/env python3
"""Prove that a cut made by extract-sprod.py is literal: put the helper back.

L2 compares what the corpus runs. A cut on a path the corpus never takes -- a
hydrate model no deck enables, a PVTSim file no deck names -- passes L2
whatever the cut did, and T130 found two real defects exactly there. This
check does not execute anything, so coverage does not matter:

  1. the call to the helper in the host AFTER the cut passes, in order, exactly
     the helper's parameter names -- so putting the body back needs no
     renaming;
  2. replacing that call with the helper's body gives, token for token, the
     host BEFORE the cut.

Together with check-byvalue-writes.py (no by-value parameter is written) and
-Wshadow on the host (no local hides a member or a global, so a name the cut
failed to pass cannot silently bind to something else), that makes the cut a
proved literal move.

Usage:
    prove-cut.py <rev-before> <rev-after> <host> <helper> [<helper> ...]

<rev-after> may be WORKTREE. Helpers are put back in the order given; list
them outermost first when one was cut out of another. Exit 0 when the token
streams are identical.
"""
import io
import re
import subprocess
import sys
import pathlib

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import braces

ROOT = pathlib.Path(__file__).resolve().parent.parent
SRC = "src/core/SisProd.cpp"
TOKEN = re.compile(r'"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\'|[A-Za-z_]\w*'
                   r'|\d[\w.]*(?:[eE][+-]?\d+)?|\.\d+(?:[eE][+-]?\d+)?'
                   r'|->|::|<<=|>>=|<<|>>|<=|>=|==|!=|&&|\|\||\+\+|--|[+\-*/%&|^]=|\S')


def source(rev):
    if rev == "WORKTREE":
        return io.open(ROOT / SRC, encoding="utf-8", errors="surrogateescape").read()
    if pathlib.Path(rev).is_file():
        # A file path stands for an arbitrary version -- used to calibrate
        # the check against a deliberately corrupted copy.
        return io.open(rev, encoding="utf-8", errors="surrogateescape").read()
    raw = subprocess.run(["git", "-C", str(ROOT), "show", f"{rev}:{SRC}"],
                         capture_output=True, check=True).stdout
    return raw.decode("utf-8", errors="surrogateescape")


def span(text, name):
    m = re.search(r'^[A-Za-z][^\n;{}]*?\bSProd::%s\(' % re.escape(name), text, re.M)
    if m is None:
        raise SystemExit("not found: SProd::%s" % name)
    o = braces.first_brace_after(text, m.start())
    return m.start(), o, braces.match(text, o)


def tokens(code):
    code = re.sub(r'/\*.*?\*/', ' ', code, flags=re.S)
    code = re.sub(r'//[^\n]*', ' ', code)
    return TOKEN.findall(code)


def param_names(signature):
    inside = signature[signature.index("(") + 1:signature.rindex(")")]
    names = []
    for p in inside.split(","):
        p = p.split("=")[0].strip()
        if p:
            names.append(re.findall(r'[A-Za-z_]\w*', p)[-1])
    return names


def main():
    if len(sys.argv) < 5:
        print(__doc__)
        return 2
    before, after, host, helpers = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:]
    old, new = source(before), source(after)

    s, o, e = span(old, host)
    want = tokens(old[o:e + 1])
    s, o, e = span(new, host)
    host_body = new[o:e + 1]

    for helper in helpers:
        hs, ho, he = span(new, helper)
        names = param_names(new[hs:ho])
        calls = list(re.finditer(r'(?<![\w.>])%s\(([^;]*)\);' % re.escape(helper), host_body))
        if len(calls) != 1:
            print("FAILED %s: %d call(s) in %s, expected 1" % (helper, len(calls), host))
            return 1
        args = [a.strip() for a in calls[0].group(1).split(",") if a.strip()]
        if args != names:
            print("FAILED %s: call passes %s, parameters are %s" % (helper, args, names))
            return 1
        body = new[ho + 1:he]
        host_body = host_body[:calls[0].start()] + body + host_body[calls[0].end():]

    got = tokens(host_body)
    if got != want:
        k = next(i for i in range(min(len(got), len(want)) + 1)
                 if i == min(len(got), len(want)) or got[i] != want[i])
        print("FAILED %s: token streams differ at token %d" % (host, k))
        print("  before: %s" % " ".join(want[max(0, k - 8):k + 8]))
        print("  put back: %s" % " ".join(got[max(0, k - 8):k + 8]))
        return 1
    print("CUT PROVED LITERAL -- %s with %s put back: %d tokens identical to %s"
          % (host, ", ".join(helpers), len(want), before))
    return 0


if __name__ == "__main__":
    sys.exit(main())
