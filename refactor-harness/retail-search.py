#!/usr/bin/env python3
"""Turn an abort-protocol helper that ends its caller back into a tail helper.

extract-search.py has two shapes for a range that contains returns. With
--tail the range ends the function, its returns stay returns and the call site
is `return helper(...)`. Without it each return becomes
`{ abortValue = <v>; return true; }`, the helper answers "must it stop?" and
the caller goes on after it.

Five ranges that DID end their function were cut in the second shape. Every
path through them returns, so the helper always answers true -- but the
compiler cannot see that across the call, and the caller now falls off its end
as far as it knows: five -Wreturn-type warnings that a vacuous Gate 1 let
through at T099 and the T130 gate caught.

This tool rewrites such a helper into the shape --tail would have produced:

  * `bool H(..., double &abortValue)` becomes `double H(...)`;
  * each `{ abortValue = <v>; return true; }` becomes `return <v>;`;
  * a second-hop call inside it, `if (inner(..., abortValue)) return true;`,
    becomes the first-hop form, because H no longer has the parameter;
  * the trailing `return false;` goes -- if any path reached it, the compiler
    now says so with -Wreturn-type on H itself;
  * the call site `double abortValue; if (H(...)) return abortValue;` becomes
    `return H(...);`.

Proof: the forward rewrite of extract-search.py, applied to the result,
reproduces the input character for character. The tool refuses to write
otherwise.

Usage:
    retail-search.py <helper> [<helper> ...]
"""
import io
import re
import sys
import pathlib

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import braces

PATH = "src/core/SisProdSteadyStateSearch.cpp"

ABORT_BLOCK = re.compile(r'(\n[ \t]*)\{\1    abortValue = ([^;]+);\1    return true;\1\}')
SECOND_HOP = re.compile(r'(\n[ \t]*)if \((\w+)\(state((?:, [^\n]*)?), abortValue\)\)(\n[ \t]*)return true;')
FIRST_HOP = re.compile(r'(\n[ \t]*)double abortValue;\1if \((\w+)\(state((?:, [^\n]*)?), abortValue\)\)(\n[ \t]*)return abortValue;')
RETURN = re.compile(r'(\n[ \t]*)return\s+([^;]+);')


def span(text, helper):
    m = re.search(r'^bool %s\((.*), double &abortValue\) \{$' % helper, text, re.M)
    if m is None:
        raise SystemExit("not an abort-protocol helper: " + helper)
    open_brace = braces.first_brace_after(text, m.start())
    # The body proper: after the signature's `{`, up to (not including) its `}`.
    return m, open_brace + 1, braces.match(text, open_brace)


def to_tail(body):
    body = ABORT_BLOCK.sub(lambda m: "%sreturn %s;" % (m.group(1), m.group(2)), body)
    body = SECOND_HOP.sub(lambda m: "%sdouble abortValue;%sif (%s(state%s, abortValue))%sreturn abortValue;"
                          % (m.group(1), m.group(1), m.group(2), m.group(3), m.group(4)), body)
    assert body.endswith("\n    return false;\n"), repr(body[-60:])
    return body[:-len("    return false;\n")]


def forward(body):
    """extract-search.py's rewrite, second-hop form, verbatim in effect."""
    # A first-hop call is left alone by the return rewrite: its `return
    # abortValue;` is exactly what the second-hop form replaces. Mask it,
    # rewrite the remaining returns, then put it back in second-hop form.
    hops = []

    def mask(m):
        hops.append(m)
        return "\x00%d\x00" % (len(hops) - 1)

    body = FIRST_HOP.sub(mask, body)
    body = RETURN.sub(lambda m: "%s{%s    abortValue = %s;%s    return true;%s}"
                      % (m.group(1), m.group(1), m.group(2).strip(), m.group(1), m.group(1)), body)
    body = re.sub(r'\x00(\d+)\x00', lambda k: "%sif (%s(state%s, abortValue))%sreturn true;" % (
        hops[int(k.group(1))].group(1), hops[int(k.group(1))].group(2),
        hops[int(k.group(1))].group(3), hops[int(k.group(1))].group(4)), body)
    return body + "    return false;\n"


def main():
    text = io.open(PATH, encoding="utf-8", errors="surrogateescape").read()
    for helper in sys.argv[1:]:
        m, a, b = span(text, helper)
        old_body = text[a:b]
        new_body = to_tail(old_body)
        if forward(new_body) != old_body:
            raise SystemExit("%s: the forward rewrite does not reproduce the input" % helper)
        if "abortValue" in re.sub(FIRST_HOP, "", new_body):
            raise SystemExit("%s: abortValue survives outside a first-hop call" % helper)
        signature = "double %s(%s) {" % (helper, m.group(1))
        assert text[m.start():a] == m.group(0), repr(text[m.start():a])
        text = text[:m.start()] + signature + new_body + text[b:]

        calls = list(re.finditer(r'(\n[ \t]*)double abortValue;\1if \(%s\(state((?:, [^\n]*)?), abortValue\)\)\n[ \t]*return abortValue;' % helper, text))
        if len(calls) != 1:
            raise SystemExit("%s: expected one call site, found %d" % (helper, len(calls)))
        c = calls[0]
        text = text[:c.start()] + "%sreturn %s(state%s);" % (c.group(1), helper, c.group(2)) + text[c.end():]
        print("%-40s -> cauda: %d retorno(s)" % (helper, len(RETURN.findall(new_body))))
    io.open(PATH, "w", encoding="utf-8", errors="surrogateescape").write(text)


if __name__ == "__main__":
    main()
