#!/usr/bin/env python3
"""Rename local identifiers one at a time, refusing collisions.

FR-039 forbids mass sed for renaming: sed matches substrings and crosses
scopes, and the specific hazard is collapsing two distinct identifiers into
one -- which compiles, emits no warning, and changes results.

This tool renames ONE identifier per step, inside ONE named scope, matching
whole identifiers only, and refuses when the new name already occurs in that
scope. It also refuses to touch anything reachable from outside the module:
SProd members, public methods, and fields of celula[], celulaG[], arq or
varGlob1D are surface (FR-038), and here they are always reached through a
`state.` field or a `.` access, which this tool never rewrites.

Usage: rename-locals.py <file> <scope> <old> <new> [...]
       scope is a function name, ALL for file scope, or Function@first-last
       to restrict the rewrite to those absolute lines.

The Function@first-last form exists because function scope is too coarse when
the same name is declared twice in sibling blocks -- two `for (int j ...)` in
one function are distinct variables, and giving both the same new name is not a
collapse but does produce a name that lies about one of them. Collisions are
still checked across the WHOLE function, so narrowing the rewrite never
narrows the safety check.

Exit code: 0 when every rename applied; 1 on the first refusal.
"""
from __future__ import annotations

import re
import sys

# Never rewritten, whatever the caller asks: an identifier that follows a dot
# is a struct field, and every struct these writers touch is shared surface.
AFTER_DOT = re.compile(r"(?<=\.)\s*$")


def scope_span(lines: list[str], scope: str) -> tuple[int, int]:
    if scope == "ALL":
        return 0, len(lines)
    for start, line in enumerate(lines):
        if re.search(rf"\b{re.escape(scope)}\s*\(", line) and "{" in "".join(
                lines[start:start + 6]):
            depth, opened = 0, False
            for probe in range(start, len(lines)):
                depth += lines[probe].count("{") - lines[probe].count("}")
                opened = opened or "{" in lines[probe]
                if opened and depth == 0:
                    return start, probe + 1
    raise LookupError(f"scope {scope!r} not found")


def rename(text: str, scope: str, old: str, new: str) -> tuple[str, int]:
    lines = text.split("\n")
    limit = None
    if "@" in scope:
        scope, span = scope.split("@", 1)
        first, last = (int(part) for part in span.split("-"))
        limit = (first - 1, last)
    start, end = scope_span(lines, scope)
    if limit is not None:
        if not (start <= limit[0] and limit[1] <= end):
            raise ValueError(f"{scope}: lines {limit[0] + 1}-{limit[1]} lie "
                             f"outside the function ({start + 1}-{end})")
        guard = "\n".join(lines[start:end])
        clash = re.compile(rf"(?<![\w.]){re.escape(new)}\b")
        if any(clash.finditer(guard)):
            raise ValueError(f"{scope}: {new!r} already occurs in the function -- "
                             f"renaming {old!r} to it would collapse two identifiers")
        start, end = limit
    body = "\n".join(lines[start:end])

    # Occurrences that are neither a field access nor inside a string literal.
    pattern = re.compile(rf"(?<![\w.]){re.escape(old)}\b")
    strings = [(m.start(), m.end()) for m in
               re.finditer(r'"(?:[^"\\]|\\.)*"|\'(?:[^\'\\]|\\.)*\'', body)]

    def in_string(position: int) -> bool:
        return any(a <= position < b for a, b in strings)

    hits = [m for m in pattern.finditer(body) if not in_string(m.start())]
    if not hits:
        raise ValueError(f"{scope}: {old!r} does not occur")

    collision = re.compile(rf"(?<![\w.]){re.escape(new)}\b")
    clashes = [m for m in collision.finditer(body) if not in_string(m.start())]
    if clashes:
        line = body[:clashes[0].start()].count("\n") + start + 1
        raise ValueError(f"{scope}: {new!r} already occurs (line {line}) -- "
                         f"renaming {old!r} to it would collapse two identifiers")

    out, last = [], 0
    for hit in hits:
        out.append(body[last:hit.start()])
        out.append(new)
        last = hit.end()
    out.append(body[last:])
    lines[start:end] = "".join(out).split("\n")
    return "\n".join(lines), len(hits)


def main() -> int:
    if len(sys.argv) < 5 or (len(sys.argv) - 2) % 3:
        print(__doc__, file=sys.stderr)
        return 2
    path = sys.argv[1]
    text = open(path, encoding="utf-8").read()
    steps = [tuple(sys.argv[i:i + 3]) for i in range(2, len(sys.argv), 3)]
    for scope, old, new in steps:
        try:
            text, count = rename(text, scope, old, new)
        except (ValueError, LookupError) as error:
            print(f"REFUSED {error}", file=sys.stderr)
            return 1
        print(f"  {scope:42} {old:16} -> {new:24} ({count} occurrence(s))")
    open(path, "w", encoding="utf-8").write(text)
    print(f"applied {len(steps)} rename(s) to {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
