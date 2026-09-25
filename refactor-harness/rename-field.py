#!/usr/bin/env python3
"""Rename a module-internal identifier reached by member access: a state field,
an updater method, or a module free function.

SC-017 (T104e) found SProd members that two modules had named differently --
`presfim` was finalPressure in three states and outletPressure in the thermal
one. Reconciling them renames identifiers that are NOT surface: fields of a
module's own state struct, methods of its updater, functions of its namespace.
rename-locals.py cannot do it -- it never rewrites after a dot, by design,
because after a dot there is usually surface.

Modes:
    rename-field.py field  <old> <new> <header> <struct> <file>[@<function>] ...
        renames the member declared in `struct <struct>` in <header>, and every
        `.old` access in the listed files (optionally only inside <function>).
    rename-field.py method <old> <new> <header> <struct> <file>[@<function>] ...
        as `field`, and also the out-of-class definition `<struct>::old(`.
    rename-field.py free   <old> <new> <file> ...
        whole-word rename of a namespace function; refuses unless every
        occurrence under src/ is in the listed files.

Refusals, before anything is written: the declaration is not found exactly
once; <new> is already a member of the struct or already occurs where it would
be written; a listed file (or function) has no occurrence to rewrite. Other
types that declare a member of the same name are LISTED -- a wrong rewrite of
one of their accesses fails to compile, and the proof below catches the rest.

Proof, after the build: every object file identical to the one before, section
by section. A field name does not reach the binary without -g. A function or
method name does, in the symbol table and in call labels -- for those, compare
the demangled disassembly with <old> substituted by <new> on the old side.
"""
import io
import re
import sys
import pathlib

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "refactor-harness"))
import braces


def read(p):
    return io.open(ROOT / p, encoding="utf-8", errors="surrogateescape").read()


def write(p, t):
    io.open(ROOT / p, "w", encoding="utf-8", errors="surrogateescape").write(t)


def struct_span(text, struct):
    m = re.search(r'\bstruct\s+%s\s*\{' % re.escape(struct), text)
    if m is None:
        raise SystemExit("REFUSED struct %s not found" % struct)
    o = text.index("{", m.start())
    return o, braces.match(text, o)


def scope(text, target):
    """(start, end) of the function body, or the whole text."""
    if "@" not in target:
        return 0, len(text)
    fn = target.split("@", 1)[1]
    # The DEFINITION, not a forward declaration: SisProd.cpp declares
    # thermalStateOf and gasLiftStateOf at the top, and taking the first match
    # scoped the rewrite to whatever brace followed line 27.
    for m in re.finditer(r'^[A-Za-z][^\n;{}]*?\b%s\s*\(' % re.escape(fn), text, re.M):
        o = braces.first_brace_after(text, m.start())
        if ";" not in text[m.end():o]:
            return m.start(), braces.match(text, o) + 1
    raise SystemExit("REFUSED function %s not defined in %s" % (fn, target))


def others_with_member(name, struct):
    found = []
    for p in list((ROOT / "src/include").glob("*.h")):
        t = io.open(p, encoding="utf-8", errors="surrogateescape").read()
        for m in re.finditer(r'\b(?:struct|class)\s+(\w+)\s*(?::[^{]*)?\{', t):
            if m.group(1) == struct:
                continue
            o = t.index("{", m.start())
            try:
                body = t[o:braces.match(t, o)]
            except Exception:
                continue
            if re.search(r'[\s*&]%s\s*(?:;|\(|\[|=)' % re.escape(name), body):
                found.append("%s (%s)" % (m.group(1), p.name))
    return found


def main():
    if len(sys.argv) < 5:
        print(__doc__)
        return 2
    mode, old, new = sys.argv[1], sys.argv[2], sys.argv[3]
    edits = {}

    if mode in ("field", "method"):
        header, struct, targets = sys.argv[4], sys.argv[5], sys.argv[6:]
        h = read(header)
        o, e = struct_span(h, struct)
        body = h[o:e]
        decl = re.compile(r'([\s*&])%s(\s*(?:;|\(|\[))' % re.escape(old))
        if len(decl.findall(body)) != 1:
            raise SystemExit("REFUSED %s: %d declaration(s) of %s in struct %s"
                             % (header, len(decl.findall(body)), old, struct))
        if re.search(r'[\s*&]%s\s*(?:;|\(|\[)' % re.escape(new), body):
            raise SystemExit("REFUSED %s is already a member of %s" % (new, struct))
        edits[header] = h[:o] + decl.sub(lambda m: m.group(1) + new + m.group(2), body) + h[e:]
        access = re.compile(r'(\.\s*)%s\b' % re.escape(old))
        qualified = re.compile(r'(\b%s::)%s\b' % (re.escape(struct), re.escape(old)))
        for target in targets:
            path = target.split("@", 1)[0]
            t = edits.get(path, read(path))
            a, b = scope(t, target)
            seg = t[a:b]
            n = len(access.findall(seg)) + (len(qualified.findall(seg)) if mode == "method" else 0)
            if n == 0:
                raise SystemExit("REFUSED %s: no `.%s` to rewrite" % (target, old))
            if re.search(r'(\.\s*)%s\b' % re.escape(new), seg):
                raise SystemExit("REFUSED %s: `.%s` already occurs there" % (target, new))
            seg = access.sub(lambda m: m.group(1) + new, seg)
            if mode == "method":
                seg = qualified.sub(lambda m: m.group(1) + new, seg)
            edits[path] = t[:a] + seg + t[b:]
            print("  %-48s .%s -> .%s  (%d)" % (target, old, new, n))
        shared = others_with_member(old, struct)
        if shared:
            print("  note: other types also declare `%s`: %s" % (old, ", ".join(shared)))
    elif mode == "free":
        files = sys.argv[4:]
        word = re.compile(r'\b%s\b' % re.escape(old))
        everywhere = [p for p in (ROOT / "src").rglob("*") if p.suffix in (".cpp", ".h")]
        outside = [str(p.relative_to(ROOT)) for p in everywhere
                   if str(p.relative_to(ROOT)) not in files
                   and word.search(io.open(p, encoding="utf-8", errors="surrogateescape").read())]
        if outside:
            raise SystemExit("REFUSED %s also occurs in %s" % (old, ", ".join(outside)))
        for path in files:
            t = read(path)
            if re.search(r'\b%s\b' % re.escape(new), t):
                raise SystemExit("REFUSED %s already occurs in %s" % (new, path))
            n = len(word.findall(t))
            if n == 0:
                raise SystemExit("REFUSED %s: no %s to rewrite" % (path, old))
            edits[path] = word.sub(new, t)
            print("  %-48s %s -> %s  (%d)" % (path, old, new, n))
    else:
        print(__doc__)
        return 2

    for path, text in edits.items():
        write(path, text)
    print("renamed %s -> %s in %d file(s)" % (old, new, len(edits)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
