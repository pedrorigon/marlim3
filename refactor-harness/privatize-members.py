#!/usr/bin/env python3
"""Make the named SProd members private WITHOUT moving any declaration (T106).

Declaration order is load-bearing twice over: the four manual copy points
(constructors, operator=, copiaSemJson) are kept in sync by it, and GCC lays
data members out in declaration order. So nothing moves. Each maximal run of
consecutive declarations whose names are all to be privatised is wrapped in
`private:` ... `public:`, where it stands; the run's leading doc comment goes
inside with it. A declaration that declares several names, some public and
some not, is left public -- splitting it would reorder.

Proof the tool demands before writing: the class body's token stream with
every `public:` / `private:` removed is identical before and after. After the
build, every object file should be identical too: access control is checked
at compile time and changes no code and, under the Itanium ABI, no layout.

Usage:
    privatize-members.py <names.json> [--dry-run]
names.json is a list of member names (T105's "no access from outside").
"""
import io
import json
import re
import sys
import pathlib

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "refactor-harness"))
import braces

HDR = ROOT / "src/include/SisProd.h"
SPEC = re.compile(r'^\s*(public|private|protected)\s*:\s*$', re.M)


def statements(body):
    """(start, end, text) of each top-level statement of the class body.

    start is where the statement's CODE begins; comments and blank lines before
    it are attached later. Inline bodies end a statement at their `}`.
    """
    out, depth, i, start = [], 0, 0, None
    code = braces  # braces._scan skips comments and strings
    marks = {k: c for k, c in code._scan(body)}
    n = len(body)
    while i < n:
        c = body[i]
        # skip comments and access specifiers when looking for a start
        if start is None:
            if body.startswith("//", i):
                i = body.find("\n", i) if body.find("\n", i) >= 0 else n
                continue
            if body.startswith("/*", i):
                i = body.find("*/", i) + 2
                continue
            m = SPEC.match(body, body.rfind("\n", 0, i) + 1)
            if m and m.start(1) <= i < m.end():
                i = m.end()
                continue
            if not c.isspace():
                start = i
        if i in marks:
            if marks[i] == "{":
                depth += 1
            elif marks[i] == "}":
                depth -= 1
                if depth == 0 and start is not None:
                    j = i + 1
                    if body[j:j + 1] == ";":
                        j += 1
                    out.append((start, j))
                    start = None
                    i = j
                    continue
        elif c == ";" and depth == 0 and start is not None:
            out.append((start, i + 1))
            start = None
        i += 1
    return out


def names_of(stmt):
    s = " ".join(re.sub(r'//[^\n]*|/\*.*?\*/', ' ', stmt, flags=re.S).split())
    paren, eq = s.find("("), s.find("=")
    if paren >= 0 and (eq < 0 or paren < eq):
        m = re.search(r'(operator\s*\S+?|~?[A-Za-z_]\w*)\s*\($', s[:paren + 1])
        return {m.group(1)} if m else set()
    m = re.match(r'^((?:static\s+|const\s+|mutable\s+)*[A-Za-z_][\w:]*(?:\s*<.*?>)?)\s+(.*?);?$', s)
    if not m:
        return set()
    names = set()
    for piece in re.split(r',(?![^<(\[]*[>)\]])', m.group(2)):
        d = re.match(r'^\s*[*&\s]*([A-Za-z_]\w*)', piece.split("=")[0])
        if d:
            names.add(d.group(1))
    return names


def leading_comment_start(body, start):
    """Move start back over the doc comment and blank lines that precede it."""
    k = body.rfind("\n", 0, start) + 1
    while k > 0:
        prev = body.rfind("\n", 0, k - 1) + 1
        line = body[prev:k].strip()
        if line.startswith(("//", "/*", "*")) or line.endswith("*/"):
            k = prev
        else:
            break
    return k


def tokens_without_specifiers(text):
    text = re.sub(r'//[^\n]*|/\*.*?\*/', ' ', text, flags=re.S)
    text = SPEC.sub(' ', text)
    return re.findall(r'\w+|\S', text)


def main():
    targets = set(json.load(open(sys.argv[1])))
    dry = "--dry-run" in sys.argv
    h = io.open(HDR, encoding="utf-8", errors="surrogateescape").read()
    s = h.index("class SProd {")
    o = h.index("{", s)
    e = braces.match(h, o)
    body = h[o + 1:e]
    # only the public region: up to the class's own `private:`
    private_at = [m.start() for m in SPEC.finditer(body) if m.group(1) == "private"]
    limit = private_at[0] if private_at else len(body)

    stmts = [(a, b) for a, b in statements(body) if b <= limit]
    flags = []
    for a, b in stmts:
        n = names_of(body[a:b])
        flags.append(bool(n) and n <= targets)
    runs, cur = [], None
    for (a, b), f in zip(stmts, flags):
        if f:
            cur = [cur[0], b] if cur else [leading_comment_start(body, a), b]
        elif cur:
            runs.append(cur)
            cur = None
    if cur:
        runs.append(cur)

    covered = set()
    for (a, b), f in zip(stmts, flags):
        if f:
            covered |= names_of(body[a:b])
    new = body
    for a, b in reversed(runs):
        new = new[:a] + "  private:\n" + new[a:b] + "\n  public:" + new[b:]
    if tokens_without_specifiers(new) != tokens_without_specifiers(body):
        raise SystemExit("REFUSED: the class's tokens changed beyond access specifiers")
    print("runs: %d | names privatised: %d of %d requested | left public (mixed or not found): %s"
          % (len(runs), len(covered), len(targets), sorted(targets - covered)))
    if not dry:
        io.open(HDR, "w", encoding="utf-8", errors="surrogateescape").write(h[:o + 1] + new + h[e:])
    return 0


if __name__ == "__main__":
    sys.exit(main())
