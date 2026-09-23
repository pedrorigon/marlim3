#!/usr/bin/env python3
"""Re-indent the call sites the extraction tools emitted at a fixed four spaces.

extract-search.py (and extract-sprod.py after it) write every call site with a
hard-coded four-space indent, whatever the depth of the block it lands in. The
code is right and the layout lies: inside a `while` two levels down, a call at
column 4 followed by the loop's next statement at column 8 reads as if the
call were outside it. GCC says so once, with -Wmisleading-indentation, at the
one site where the next line is an `if` (T130); the others mislead only the
reader.

Only leading whitespace changes. A line is touched only if it is one of the
shapes the tools emit -- `double abortValue;`, `if (helper(state...))` and the
return under it, `return helper(state...);`, `helper(state...);` -- and its
indent is not four spaces per open brace. Proof: the file with all whitespace
removed is identical before and after, and the tool refuses to write otherwise.

Usage:
    reindent-calls.py [--dry-run] <file> [<file> ...]
"""
import io
import re
import sys
import pathlib

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import braces

# The abort protocol's call site, and the five tail calls retail-search.py
# rewrote from it. A generic `helper(state...)` shape was tried first and
# matched original code in the thermal and gas-lift modules, whose layout is
# not this tool's business.
SHAPES = re.compile(r'^(double abortValue;'
                    r'|if \(\w+\(state\b.*, abortValue\)\)'
                    r'|return (?:bracketReverseRoot|bracketTertiaryPressureToPressureRoot'
                    r'|bracketInjectionRoot2|bracketInjectionRoot5|bracketSecondaryBranchRoot)\(state\b.*\);)$')
IF_SHAPE = re.compile(r'^if \(\w+\(state\b.*, abortValue\)\)$')


def depth_at_line_starts(text):
    """Brace depth at the start of every line, braces in code only."""
    opens = [0] * (text.count("\n") + 1)
    events = {}
    stack = []
    for i, c in braces._scan(text):
        # A namespace's braces do not indent: the code inside sits at column 0.
        if c == "{":
            line_start = text.rfind("\n", 0, i) + 1
            stack.append(0 if re.match(r'\s*(?:inline\s+)?namespace\b', text[line_start:i]) else 1)
            events[i] = stack[-1]
        elif c == "}":
            events[i] = -stack.pop()
    depth, line = 0, 0
    for i, ch in enumerate(text):
        if i in events:
            depth += events[i]
        if ch == "\n":
            line += 1
            opens[line] = depth
    return opens


def main():
    args = sys.argv[1:]
    dry = "--dry-run" in args
    for path in [a for a in args if a != "--dry-run"]:
        text = io.open(path, encoding="utf-8", errors="surrogateescape").read()
        depth = depth_at_line_starts(text)
        lines = text.split("\n")
        changed = []
        for k, row in enumerate(lines):
            body = row.strip()
            under_if = k > 0 and IF_SHAPE.match(lines[k - 1].strip()) and body in ("return abortValue;", "return true;")
            if not (SHAPES.match(body) or under_if):
                continue
            want = 4 * depth[k]
            # The return under `if (helper(...))` is that if's body: one level in.
            if k > 0 and IF_SHAPE.match(lines[k - 1].strip()) and body.startswith("return "):
                want += 4
            have = len(row) - len(row.lstrip(" "))
            if have != want and not row.startswith("\t"):
                lines[k] = " " * want + body
                changed.append((k + 1, have, want))
        new = "\n".join(lines)
        if re.sub(r'\s', '', new) != re.sub(r'\s', '', text):
            raise SystemExit(path + ": a change reached beyond whitespace")
        for n, have, want in changed:
            print("%s:%d  %d -> %d" % (path, n, have, want))
        print("%s: %d linha(s)" % (path, len(changed)))
        if not dry and changed:
            io.open(path, "w", encoding="utf-8", errors="surrogateescape").write(new)


if __name__ == "__main__":
    main()
