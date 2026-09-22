#!/usr/bin/env python3
"""Which locals cross into a line range of a C++ function, respecting scope.

Five wrong versions preceded this one, and each is recorded because the shape of
the mistake is the useful part:

  1. "every declaration textually above the range" -- over-reports. A `double
     rhol` inside `if (chute < 0) { }` is gone before the range starts, and
     listing it invents a parameter that does not exist.
  2. "only declarations at the function's own top level" -- the opposite, and
     worse: the missing ones do not surface as a wrong signature but as code
     that does not compile. `int corrigechute = 1;` sits inside an outer while,
     at depth two, and a range inside that same while reads it.
  3. depth instead of the brace STACK -- the two arms of a sign split both sit
     at depth 3, and `} else if (c) {` closes one and opens the other inside a
     single line, so a declaration in the first arm looked live in the second.
  4. `int i` in a for-init counted as outliving its loop, because the brace that
     scopes it opens on the same line.
  5. reference parameters missed entirely, because the `&` sits between type and
     name -- which only surfaced when a generated helper became the thing being
     cut further.

One known limitation, deliberate: strip_control_headers also hides a for-init
from ranges INSIDE that loop, where `i` genuinely is visible. That direction
fails safe -- the compiler names the missing parameter immediately -- while the
opposite direction invents parameters silently.
"""
import io
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import braces

KEYWORDS = {"return", "else", "if", "while", "for", "case", "new", "delete",
            "do", "switch", "break", "continue", "const", "static", "sizeof"}
# The `[&*]?` is not decoration: without it the pattern misses every reference
# parameter -- `double &chutelim` -- which was fine while only by-value
# signatures were read, and stopped being fine the moment a generated helper
# became the thing being cut.
DECL = re.compile(r'(?<![.\w])([A-Za-z_][A-Za-z_0-9:<>]*)\s+[&*]?([A-Za-z_][A-Za-z_0-9]*)\s*(?:=[^=]|;|,|\))')
CONTROL_HEADER = re.compile(r'(?<![.\w])(?:for|if|while|switch|catch)\s*\((?:[^()]|\([^()]*\))*\)')


def strip_comments(s):
    return re.sub(r'/\*.*?\*/', '', re.sub(r'//[^\n]*', '', s), flags=re.S)


def strip_control_headers(line):
    """Blank out `for (int i = ...; ...)` and friends."""
    return CONTROL_HEADER.sub(lambda m: " " * len(m.group(0)), line)


def _stacks(lines):
    """For each line, the stack of `{` positions open before it starts.

    Depth alone is not enough: two sibling arms sit at the same depth, and
    `} else if (c) {` closes one and opens the other inside one line.
    """
    text = "\n".join(lines)
    events = list(braces._scan(text))
    out, stack, cursor = [], [], 0
    offsets, pos = [], 0
    for line in lines:
        offsets.append(pos)
        pos += len(line) + 1
    for start in offsets:
        while cursor < len(events) and events[cursor][0] < start:
            index, char = events[cursor]
            if char == "{":
                stack.append(index)
            elif stack:
                stack.pop()
            cursor += 1
        out.append(tuple(stack))
    return out


def visible_decls(lines, at):
    """Declarations visible on line `at` (1-based), as {name: (type, line)}."""
    stacks = _stacks(lines)
    target = stacks[at - 1]
    out = {}
    for k in range(at - 2, -1, -1):
        if stacks[k] != target[:len(stacks[k])]:
            continue
        for m in DECL.finditer(strip_control_headers(strip_comments(lines[k]))):
            if m.group(1) in KEYWORDS or m.group(2) in KEYWORDS:
                continue
            out.setdefault(m.group(2), (m.group(1), k + 1))
    return out


def declared_in(lines, a, b):
    """Declarations made directly in lines a..b, in the range's own scope."""
    stacks = _stacks(lines)
    base = stacks[a - 1]
    out = {}
    for k in range(a - 1, b):
        if stacks[k] != base:
            continue
        for m in DECL.finditer(strip_control_headers(strip_comments(lines[k]))):
            if m.group(1) in KEYWORDS or m.group(2) in KEYWORDS:
                continue
            out.setdefault(m.group(2), (m.group(1), k + 1))
    return out


def _escapes(name, tail):
    """True when the tail reads THIS declaration rather than one of its own.

    Two arms of a sign split each declare `int kontaReverso = 0;`. Different
    variables, different scopes; asking only "does the name appear later" calls
    the first an out-parameter of the second.
    """
    m = re.search(rf'(?<![.\w]){name}\b', tail)
    if m is None:
        return False
    line_start = tail.rfind("\n", 0, m.start()) + 1
    return not re.search(r'(?<![.\w])[A-Za-z_][A-Za-z_0-9:<>]*\s+$', tail[line_start:m.start()])


def interface(path, a, b):
    """(crossing_in, crossing_out) for lines a..b of the function in `path`."""
    lines = io.open(path, encoding="utf-8", errors="replace").read().split("\n")
    blk = strip_comments("\n".join(lines[a - 1:b]))
    tail = strip_comments("\n".join(lines[b:]))
    crossing_in = []
    for name, (ty, ln) in visible_decls(lines, a).items():
        if not re.search(rf'(?<![.\w]){name}\b', blk):
            continue
        written = re.search(rf'(?<![.\w]){name}\s*(=[^=]|\+=|-=|\*=|/=|\+\+|--)', blk) is not None
        after = re.search(rf'(?<![.\w]){name}\b', tail) is not None
        crossing_in.append((ty, name, written, after))
    crossing_out = [(ty, name) for name, (ty, _) in declared_in(lines, a, b).items()
                    if _escapes(name, tail)]
    return crossing_in, crossing_out


def report(path, a, b, label):
    ins, outs = interface(path, a, b)
    print(f"{label} ({b - a + 1} lines)")
    print("   in  (%d): %s" % (len(ins), ", ".join(
        f"{t} {n}{'(W)' if w else ''}{'(after)' if d else ''}" for t, n, w, d in ins)))
    print("   out (%d): %s" % (len(outs), ", ".join(f"{t} {n}" for t, n in outs)))


def _calibrate():
    cases = [
        ("closed before: inner does not cross", '''void f() {
    double outer = 0.;
    if (a) {
        double inner = 1.;
    }
    outer = 2.;
    g(outer);
}
''', 6, 7, {"outer"}, set()),
        ("outer scope still open: crosses", '''void f() {
    while (x) {
        int c = 1;
        while (c == 1) {
            c = 0;
        }
    }
}
''', 4, 6, {"c"}, set()),
        ("sibling arm does not leak", '''void f() {
    if (a) {
        int k = 0;
        k++;
    } else if (b) {
        int k = 0;
        k++;
    }
}
''', 6, 7, set(), set()),
        ("for-init does not outlive its loop", '''void f() {
    for (int i = n; i > 0; i--) {
        g(i);
    }
    h(1);
    k(2);
}
''', 5, 6, set(), set()),
        ("declared inside and read after: escapes", '''void f() {
    double a = 1.;
    double b = a;
    g(a);
    h(b);
}
''', 3, 4, {"a"}, {"b"}),
        ("reference parameter is seen", '''bool h(const S &state, double amp, double &lim, double &val) {
    lim = 1.;
    val = 2.;
}
''', 2, 3, {"lim", "val"}, set()),
    ]
    bad = 0
    tmp = Path("/tmp/_iface_calibrate.cpp")
    for name, text, a, b, want_in, want_out in cases:
        tmp.write_text(text, encoding="utf-8")
        ins, outs = interface(str(tmp), a, b)
        got_in = {n for _, n, _, _ in ins}
        got_out = {n for _, n in outs}
        if got_in != want_in or got_out != want_out:
            print(f"  FAILED {name}: in={sorted(got_in)} out={sorted(got_out)}, "
                  f"expected in={sorted(want_in)} out={sorted(want_out)}")
            bad += 1
    tmp.unlink(missing_ok=True)
    print(f"iface.py calibration: {len(cases) - bad}/{len(cases)}")
    return bad


if __name__ == "__main__":
    raise SystemExit(1 if _calibrate() else 0)
