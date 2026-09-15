#!/usr/bin/env python3
"""Move one block out of an SProd method into a named private helper.

Why this exists
---------------
T085 decomposes RenovaMassPerm, 1,619 lines, into helpers of at most 200. The
acceptance is not "it compiles": it is that every extraction preserves
evaluation order, scope and variable lifetime, rewrites no floating-point
expression, and keeps L2 green ON ITS OWN rather than only at the end of the
batch.

Doing that by hand seven times is how a stray brace or a dropped line gets in.
So the cut is mechanical and the result is PROVEN: the body installed in the
helper is tokenised and compared against the body that was removed. Whitespace
and indentation change, nothing else may.

The block is replaced by a call, in place, so the dispatch that selects it stays
visible in the original function. That is the point of the decomposition: the
chain of accessory kinds should read as a chain, not as 1,200 lines of bodies.

Usage:
  extract-helper.py <file> <method> <first-line> <last-line> <helper> <params> <args>

  first-line / last-line are 1-based and RELATIVE to the method's first line.
  params  a C++ parameter list, e.g. "int i, double bo, double &fwI"
  args    the matching call arguments, e.g. "i, bo, fwI"

Exit code: 0 only when the token proof passes.
"""
import re, sys, subprocess

def tokens(text):
    text = re.sub(r'/\*.*?\*/', ' ', text, flags=re.S)
    text = re.sub(r'//[^\n]*', ' ', text)
    return re.findall(r'[A-Za-z_]\w*|\d+\.?\d*(?:[eE][-+]?\d+)?|[^\s\w]', text)


def match_brace(s, i):
    """i indexes a '{'; returns the index of the '}' that closes it.

    Character level, not line level, and that is the whole point. A line like

        } else if (cond) {

    closes one block and opens another WITHIN the same line, so a line-based
    depth counter never sees zero there and happily walks past the end of the
    arm into the next one. That mistake produced a helper carrying half of one
    branch and the header of the next, three times, before this function
    existed.
    """
    depth = 0
    while i < len(s):
        if s[i] == '{':
            depth += 1
        elif s[i] == '}':
            depth -= 1
            if depth == 0:
                return i
        i += 1
    sys.exit("unbalanced braces")

def chain_arms(src, start_pattern):
    """Every arm of one if / else-if chain, as (condition_line, body_first, body_last)."""
    m = re.search(start_pattern, src, re.M)
    if not m:
        sys.exit(f"chain start not found: {start_pattern}")
    line_of = lambda idx: src[:idx].count('\n') + 1
    arms = []
    pos = m.start()
    while True:
        open_brace = src.index('{', pos)
        close_brace = match_brace(src, open_brace)
        arms.append((line_of(pos), line_of(open_brace) + 1, line_of(close_brace) - 1))
        rest = src[close_brace + 1:]
        skip = re.match(r'\s*(?://[^\n]*\n\s*)*', rest).group(0)
        if not re.match(r'else\b', rest[len(skip):]):
            return arms
        pos = close_brace + 1 + len(skip)

def method_span(src, method):
    m = re.search(r'^[A-Za-z_][\w:<>,&\* ]*\bSProd::' + re.escape(method) + r'\s*\(', src, re.M)
    if not m:
        sys.exit(f"method not found: {method}")
    start = src.rindex('\n', 0, m.start()) + 1
    depth = 0; seen = False; j = start
    while j < len(src):
        if src[j] == '{': depth += 1; seen = True
        elif src[j] == '}':
            depth -= 1
            if seen and depth == 0: break
        j += 1
    return start, j + 1

def main():
    if len(sys.argv) == 9 and sys.argv[3] == '--arm':
        # extract-helper.py <file> <method> --arm <n> <helper> <params> <args> <chain-regex>
        path, method, _, n, helper, params, args, pattern = sys.argv[1:]
        src = open(path).read()
        start, end = method_span(src, method)
        body = src[start:end]
        arms = chain_arms(body, pattern)
        idx = int(n)
        if idx < 1 or idx > len(arms):
            sys.exit(f"arm {idx} out of range: the chain has {len(arms)}")
        _, a, b = arms[idx - 1]
        print(f"    arm {idx} of {len(arms)}: body lines {a}-{b}")
    elif len(sys.argv) == 8:
        path, method, a, b, helper, params, args = sys.argv[1:]
        a, b = int(a), int(b)
        src = open(path).read()
        start, end = method_span(src, method)
    else:
        sys.exit(__doc__)
    lines = src[start:end].split('\n')

    block = '\n'.join(lines[a-1:b])
    if not block.strip():
        sys.exit("the selected block is empty")

    indent = re.match(r'[ \t]*', lines[a-1]).group(0)
    # The call keeps the indentation the block had.
    call = f"{indent}{helper}({args});"
    new_method_lines = lines[:a-1] + [call] + lines[b:]

    # Dedent the body by the common leading whitespace, so the helper reads at
    # its own level. Tokens are unaffected; this is purely cosmetic.
    body_lines = [l for l in block.split('\n')]
    widths = [len(re.match(r'[ \t]*', l).group(0)) for l in body_lines if l.strip()]
    cut = min(widths) if widths else 0
    body = '\n'.join(('    ' + l[cut:]) if l.strip() else '' for l in body_lines)

    helper_src = f"void SProd::{helper}({params}) {{\n{body}\n}}\n\n"

    # Proof: the installed body must tokenise identically to the removed block.
    if tokens(body) != tokens(block):
        sys.exit("TOKEN MISMATCH: the extracted body is not the block that was removed")

    rebuilt = '\n'.join(new_method_lines)
    out = src[:start] + helper_src + rebuilt + src[end:]
    open(path, 'w').write(out)

    print(f"OK  {helper}: {b-a+1} lines moved, {len(tokens(block))} tokens identical")
    print(f"    {method} is now {len(new_method_lines)} lines")

if __name__ == '__main__':
    main()
