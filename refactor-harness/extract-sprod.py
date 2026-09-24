#!/usr/bin/env python3
"""Lift a line range out of an SProd member into a new SProd member, in place.

The in-place counterpart of extract-search.py: that one writes free functions
taking a state struct into a module; this one stays inside SisProd.cpp and
writes `void SProd::helper(...)`, which is what a function being decomposed
BEFORE it moves needs -- member access stays free, and the move tool rewrites
it later in one proved step.

Written for SolveTrans, which needs about eight cuts and whose stage requires an
L2 after every one of them.

Usage:
    extract-sprod.py <function> <first> <last> <helper> [--body] [--doc TEXT]

Line numbers are relative to the function's own first line. --body treats the
first and last lines as the braces of an arm and lifts what is between them.
The range must contain no `return`: a helper cannot return from its caller, and
this tool will not invent an abort protocol for a function that does not need
one.
"""
import argparse
import io
import pathlib
import re
import textwrap
import sys
import tempfile

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import braces
import iface

SRC = "src/core/SisProd.cpp"
HDR = "src/include/SisProd.h"


def first_sentence(doc):
    flat = " ".join(doc.split())
    m = re.search(r'^(.*?\.)(?:\s|$)', flat)
    return m.group(1) if m else flat


def function_span(text, name):
    # `SProd &` and an escaped name admit operator= (T101); the rest were the
    # only hosts before.
    m = re.search(rf'^(?:void|double|int|bool|SProd &)\s*SProd::{re.escape(name)}\(', text, re.M)
    if m is None:
        raise SystemExit("function not found: " + name)
    return m.start(), braces.match(text, braces.first_brace_after(text, m.start()))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("function")
    ap.add_argument("first", type=int)
    ap.add_argument("last", type=int)
    ap.add_argument("helper")
    ap.add_argument("--body", action="store_true")
    ap.add_argument("--doc", default="")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--extra", action="append", default=[],
                    help="a parameter the analyser cannot see, as 'type name'. "
                         "iface.py hides for-init declarations on purpose -- right for "
                         "ranges AFTER the loop, wrong for ranges INSIDE it -- so a "
                         "loop body needs its loop variable passed by hand")
    args = ap.parse_args()

    text = io.open(SRC, encoding="utf-8", errors="surrogateescape").read()
    fa, fb = function_span(text, args.function)
    fn_text = text[fa:fb + 1]
    fn_lines = fn_text.split("\n")
    first, last = args.first, args.last
    if args.body:
        first, last = first + 1, last - 1
    while first <= last and not fn_lines[first - 1].strip():
        first += 1
    while last >= first and not fn_lines[last - 1].strip():
        last -= 1
    block = fn_lines[first - 1:last]

    code = re.sub(r'/\*.*?\*/', '', re.sub(r'//[^\n]*', '', "\n".join(block)), flags=re.S)
    if re.search(r'(?<![.\w])return\b', code):
        raise SystemExit("the range contains a return; this tool does not lift those")

    scratch = pathlib.Path(tempfile.gettempdir()) / "_extract_sprod_scope.cpp"
    scratch.write_text(fn_text, encoding="utf-8", errors="surrogateescape")
    crossing_in, crossing_out = iface.interface(str(scratch), first, last)
    scratch.unlink(missing_ok=True)
    if crossing_out:
        raise SystemExit("declared inside and read after: %s -- not lifting"
                         % ", ".join(n for _, n in crossing_out))

    params = []
    extra_names = []
    for spec in args.extra:
        ty, name = spec.rsplit(" ", 1)
        params.append(f"{ty} {name}")
        extra_names.append(name)
    for ty, name, written, _after in crossing_in:
        if ty.endswith("*"):
            # A pointer goes by value -- the pointee is shared either way -- or
            # by reference when the range reseats it. `const T * &` would not
            # even bind to the caller's `T *`.
            params.append(f"{ty}{'&' if written else ''}{name}")
            continue
        if ty in ("int", "double", "bool", "float", "long", "char"):
            # A scalar goes by reference when the range writes it, by value
            # otherwise -- and check-byvalue-writes.py proves the "otherwise".
            params.append(f"{ty} {'&' if written else ''}{name}")
        else:
            # An object goes by plain reference: the helper then works on the
            # very object the original code did, whatever it does to it.
            # `const T &` for "not written" relied on the analyser seeing every
            # write, and it does not see method calls or member assignments --
            # auxMiniTab's `fluC.atualizaPropComp(...)` and
            # `flui.miniTabDin.x[i][j] = ...` came out const (T101, dry run).
            params.append(f"{ty} &{name}")
    signature = ", ".join(params)
    call_args = ", ".join(extra_names + [n for _, n, _, _ in crossing_in])

    print(f"{args.helper}: linhas {first}..{last} da funcao ({len(block)}), "
          f"parametros: {signature or '(nenhum)'}")
    print(f"  primeira: {block[0].strip()[:80]}")
    print(f"  ultima  : {block[-1].strip()[:80]}")
    if args.dry_run:
        return

    pad = len(block[0]) - len(block[0].lstrip(" "))
    body = "\n".join((("    " + r[pad:]) if r.startswith(" " * pad) else
                      ("    " + r.strip() if r.strip() else "")) for r in block)
    doc = args.doc or f"Part of {args.function}."
    comment = "".join(f"/// {line}\n" if line else "///\n" for line in doc.split("\n"))
    helper = f"{comment}void SProd::{args.helper}({signature}) {{\n{body}\n}}\n\n"

    base = braces.line_of(text, fa)
    all_lines = text.split("\n")
    a, b = base + first - 1, base + last - 1
    indent = " " * pad
    all_lines[a:b + 1] = [f"{indent}{args.helper}({call_args});"]
    k = base
    while k > 0 and all_lines[k - 1].lstrip().startswith("//"):
        k -= 1
    all_lines[k:k] = helper.rstrip("\n").split("\n") + [""]
    io.open(SRC, "w", encoding="utf-8", errors="surrogateescape").write("\n".join(all_lines))

    hdr = io.open(HDR, encoding="utf-8", errors="surrogateescape").read()
    m = re.search(rf'^(\s*)(?:void|double|int|bool|SProd &)\s*{re.escape(args.function)}\(', hdr, re.M)
    if m is None:
        raise SystemExit(f"declaration of {args.function} not found in {HDR}")
    # The header gets the doc's first SENTENCE, rewrapped. It used to get the
    # first LINE, which for a wrapped doc is half a sentence: nine T100
    # helpers were declared under comments that stopped mid-clause.
    decl = "".join(f"    /// {row}\n"
                   for row in textwrap.wrap(first_sentence(doc), 76, break_on_hyphens=False))
    decl += f"    void {args.helper}({signature});\n"
    # Above the host's own doc comment, not between it and the host: inserting
    # at the declaration left montasistema's and SolveTrans's comments sitting
    # over the first helper instead (T100, both repaired).
    at = hdr.rfind("\n", 0, m.start()) + 1
    while True:
        prev = hdr.rfind("\n", 0, at - 1) + 1
        if at > 0 and hdr[prev:at].lstrip().startswith("//"):
            at = prev
        else:
            break
    hdr = hdr[:at] + decl + hdr[at:]
    io.open(HDR, "w", encoding="utf-8", errors="surrogateescape").write(hdr)
    print("  extraido e declarado")


if __name__ == "__main__":
    main()
