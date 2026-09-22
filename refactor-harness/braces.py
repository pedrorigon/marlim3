#!/usr/bin/env python3
"""Character-level brace matching over C++ text.

Line-level depth counting is wrong here and was wrong once in practice, in T085:
`} else if (cond) {` closes and reopens inside one line, so a per-line counter
never sees zero and walks past the end of the arm it was measuring. Every
boundary found with this module matches a specific opening brace to its own
closing brace, with strings and comments skipped.

Lives in refactor-harness/ and not in a scratch directory. The first version did
live in one, and it was deleted between sessions along with the calibration that
made it trustworthy.
"""


def _scan(text):
    """Yield (index, char) for braces outside strings, chars and comments."""
    i, n = 0, len(text)
    while i < n:
        c = text[i]
        if text.startswith("//", i):
            i = text.find("\n", i)
            if i < 0:
                return
            continue
        if text.startswith("/*", i):
            i = text.find("*/", i)
            if i < 0:
                return
            i += 2
            continue
        if c in "\"'":
            quote, i = c, i + 1
            while i < n and text[i] != quote:
                i += 2 if text[i] == "\\" else 1
            i += 1
            continue
        if c in "{}":
            yield i, c
        i += 1


def match(text, open_index):
    """Index of the `}` that closes the `{` at open_index."""
    assert text[open_index] == "{", text[open_index - 20:open_index + 20]
    depth = 0
    for i, c in _scan(text):
        if i < open_index:
            continue
        depth += 1 if c == "{" else -1
        if depth == 0:
            return i
    raise ValueError("unbalanced from index %d" % open_index)


def line_of(text, index):
    """0-based line number containing index."""
    return text.count("\n", 0, index)


def first_brace_after(text, index):
    """Index of the first `{` at or after index, outside strings and comments."""
    for i, c in _scan(text):
        if i >= index and c == "{":
            return i
    raise ValueError("no brace after %d" % index)


def _skip_trivia(text, i):
    while i < len(text):
        if text[i].isspace():
            i += 1
        elif text.startswith("//", i):
            j = text.find("\n", i)
            i = len(text) if j < 0 else j + 1
        elif text.startswith("/*", i):
            j = text.find("*/", i)
            i = len(text) if j < 0 else j + 2
        else:
            return i
    return i


def chain_match(text, open_index):
    """Index of the last character of a whole if / else-if / else chain.

    match() answers a different question -- it ends one ARM -- and that is the
    honest answer for a brace. A chain is a sequence of arms, so this walks
    them: close an arm, look for `else`, continue into the next braced arm or to
    the end of a braceless one.
    """
    end = match(text, open_index)
    while True:
        i = _skip_trivia(text, end + 1)
        if not text.startswith("else", i) or (i + 4 < len(text) and
                                              (text[i + 4].isalnum() or text[i + 4] == "_")):
            return end
        i = _skip_trivia(text, i + 4)
        if text.startswith("if", i):
            i = _skip_trivia(text, i + 2)
            assert text[i] == "(", text[i:i + 20]
            depth = 0
            while True:
                if text[i] == "(":
                    depth += 1
                elif text[i] == ")":
                    depth -= 1
                    if depth == 0:
                        break
                i += 1
            i = _skip_trivia(text, i + 1)
        if text[i] == "{":
            end = match(text, i)
            continue
        depth = 0
        while i < len(text):
            if text[i] == "(":
                depth += 1
            elif text[i] == ")":
                depth -= 1
            elif text[i] == ";" and depth == 0:
                return i
            i += 1
        raise ValueError("braceless arm never terminated")


def _calibrate():
    """A matcher nobody has tried to fool is not a matcher."""
    cases = [
        ("braced else", 'void f() {\n if (a) {\n x();\n } else if (b) {\n y();\n } else {\n z();\n }\n}\n', 'if (a)', 8),
        ("braceless else", 'void f() {\n if (a) {\n x();\n } else\n z();\n}\n', 'if (a)', 5),
        ("no else", 'void f() {\n if (a) {\n x();\n }\n y();\n}\n', 'if (a)', 4),
        ("else-if braceless tail", 'void f() {\n if (a) {\n x();\n } else if (b)\n y();\n}\n', 'if (a)', 5),
        ("brace inside a comment", 'void f() {\n if (a) { // not a brace: {\n x();\n } else {\n z();\n }\n}\n', 'if (a)', 6),
        ("elsewhere is not else", 'void f() {\n if (a) {\n x();\n }\n elsewhere();\n}\n', 'if (a)', 4),
    ]
    bad = 0
    for name, text, anchor, want in cases:
        o = text.index("{", text.index(anchor))
        got = line_of(text, chain_match(text, o)) + 1
        if got != want:
            print(f"  FAILED {name}: ended at line {got}, expected {want}")
            bad += 1
    print(f"braces.py calibration: {len(cases) - bad}/{len(cases)}")
    return bad


if __name__ == "__main__":
    raise SystemExit(1 if _calibrate() else 0)
