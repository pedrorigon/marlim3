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

A sixth version fixed multi-declarations. `time_point begin, end;` declares
two names and the pattern only saw the first, so `end` was invisible to every
range that read it. The compiler caught it -- but not for the reason one would
hope: with `using namespace std` in scope, the undeclared `end` resolved to
std::end, and the error was a type mismatch rather than an unknown name. A name
the analyser cannot see can bind to a DIFFERENT entity, and it only failed here
because std::end is a function. So every name in a declarator list is collected
now.

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
# The `&` is not decoration: without it the pattern misses every reference
# parameter -- `double &chutelim` -- which was fine while only by-value
# signatures were read, and stopped being fine the moment a generated helper
# became the thing being cut.
#
# The `*` is KEPT in the type. It used to be matched and dropped, the same as
# `&`, which is right for a reference -- the cut decides by-reference from the
# writes -- and wrong for a pointer: montasistema's `double *compfonte` came out
# as a parameter `double compfonte` (T100, caught in a dry run).
DECL = re.compile(r'(?<![.\w])([A-Za-z_][A-Za-z_0-9:<>]*)\s+(&|\*{1,3})?([A-Za-z_][A-Za-z_0-9]*)\s*(?:=[^=]|;|,|\))')


def _typed(base, mark):
    return base + " " + mark if mark and mark.startswith("*") else base
CONTROL_HEADER = re.compile(r'(?<![.\w])(?:for|if|while|switch|catch)\s*\((?:[^()]|\([^()]*\))*\)')


def strip_comments(s):
    return re.sub(r'/\*.*?\*/', '', re.sub(r'//[^\n]*', '', s), flags=re.S)


def strip_control_headers(line):
    """Blank out `for (int i = ...; ...)` and friends."""
    return CONTROL_HEADER.sub(lambda m: " " * len(m.group(0)), line)



def _declarations(line):
    """(type, name) for every name a line declares, including `T a, b, c;`.

    DECL finds `T a`. When the declarator list continues with a comma, the rest
    of the statement up to its `;` is split on TOP-LEVEL commas -- an
    initialiser like `f(x, y)` has commas of its own -- and each piece yields
    its declarator name.
    """
    out = []
    for m in DECL.finditer(line):
        base, first = m.group(1), m.group(3)
        if base in KEYWORDS or first in KEYWORDS:
            continue
        out.append((_typed(base, m.group(2)), first))
        rest = line[m.end(3):]
        depth, cut = 0, None
        for i, ch in enumerate(rest):
            if ch in "([{":
                depth += 1
            elif ch in ")]}":
                depth -= 1
                if depth < 0:
                    cut = i
                    break
            elif ch == ";" and depth == 0:
                cut = i
                break
        tail = rest[:cut] if cut is not None else rest
        pieces, depth, start = [], 0, 0
        for i, ch in enumerate(tail):
            if ch in "([{":
                depth += 1
            elif ch in ")]}":
                depth -= 1
            elif ch == "," and depth == 0:
                pieces.append(tail[start:i])
                start = i + 1
        pieces.append(tail[start:])
        for piece in pieces[1:]:
            # A declarator-list continuation is ONE name -- `, end` or `, y = 2.`.
            # A parameter list is `T1 a, T2 b`, where each piece carries its
            # own type: there the first identifier is the TYPE, and taking it
            # as the name turned `double alfRev, int nrede` into parameters
            # called `double` and `int`. Those pieces are separate
            # declarations that DECL already found on its own, so skip them.
            name = re.match(r'\s*(&|\*{1,3})?([A-Za-z_]\w*)\s*(?:=|\[|$)', piece)
            if name and name.group(2) not in KEYWORDS:
                # `double *p, *q` and `double *p, q`: each declarator has its own star.
                out.append((_typed(base, name.group(1)), name.group(2)))
    return out


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
        for ty, name in _declarations(strip_control_headers(strip_comments(lines[k]))):
            out.setdefault(name, (ty, k + 1))
    return out


def declared_in(lines, a, b):
    """Declarations made directly in lines a..b, in the range's own scope."""
    stacks = _stacks(lines)
    base = stacks[a - 1]
    out = {}
    for k in range(a - 1, b):
        if stacks[k] != base:
            continue
        for ty, name in _declarations(strip_control_headers(strip_comments(lines[k]))):
            out.setdefault(name, (ty, k + 1))
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
        ("multi-declaration: second name is seen", '''void f() {
    clock_t begin, end;
    end = now();
    g(end - begin);
}
''', 4, 4, {"begin", "end"}, set()),
        ("initialiser commas are not declarators", '''void f() {
    double x = h(a, b), y = 2.;
    g(x + y);
}
''', 3, 3, {"x", "y"}, set()),
        ("parameter list is not a declarator list", '''void f(double a, int n, S s) {
    g(a, n);
}
''', 2, 2, {"a", "n"}, set()),
        ("reference parameter is seen", '''bool h(const S &state, double amp, double &lim, double &val) {
    lim = 1.;
    val = 2.;
}
''', 2, 3, {"lim", "val"}, set()),
        ("pointer keeps its star", '''void f(double *comp, int *posic, int n) {
    double *a, b, **c;
    posic[0] = n;
    g(comp[0], a, b, c);
}
''', 3, 4, {"comp", "posic", "n", "a", "b", "c"}, set(),
         {"comp": "double *", "posic": "int *", "n": "int", "a": "double *", "b": "double", "c": "double **"}),
    ]
    bad = 0
    tmp = Path("/tmp/_iface_calibrate.cpp")
    for name, text, a, b, want_in, want_out, *types in cases:
        tmp.write_text(text, encoding="utf-8")
        ins, outs = interface(str(tmp), a, b)
        got_in = {n for _, n, _, _ in ins}
        got_out = {n for _, n in outs}
        got_types = {n: ty for ty, n, _, _ in ins}
        if types and any(got_types.get(n) != ty for n, ty in types[0].items()):
            print(f"  FAILED {name}: types {got_types}, expected {types[0]}")
            bad += 1
        elif got_in != want_in or got_out != want_out:
            print(f"  FAILED {name}: in={sorted(got_in)} out={sorted(got_out)}, "
                  f"expected in={sorted(want_in)} out={sorted(want_out)}")
            bad += 1
    tmp.unlink(missing_ok=True)
    print(f"iface.py calibration: {len(cases) - bad}/{len(cases)}")
    return bad


if __name__ == "__main__":
    raise SystemExit(1 if _calibrate() else 0)
