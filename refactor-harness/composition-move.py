#!/usr/bin/env python3
"""Move stage 9 composition bodies out of SisProd.cpp and prove the move literal.

Same method as transient-move.py: rewrite the signature and the SProd member
accesses through an INVERTIBLE substitution, then prove the move by applying
the inverse to the installed body and comparing tokens with the pre-move
commit. The generic machinery -- tokenizer, alpha-equivalence, comment-safe
substitution, the guard for locals that hide members -- is imported from
transient-move.py, which owns and calibrates it; only what names the stage is
here.

The member map is DERIVED, not written: it is read from the `T &field; //
member` lines of SisProdComposition.h, whose seventeen fields were themselves
derived from measure-members.py. So is the adapter (`adapter` mode prints it),
in the header's declaration order, which C++20 designated initialisers
require. Table, header and adapter cannot drift apart.

Usage:
    composition-move.py extract  <function> <sisprod.cpp> <fragment.cpp>
    composition-move.py install  <function> <sisprod.cpp> <module.cpp> <out.cpp>
    composition-move.py delegate <function> <sisprod.cpp> <rewritten.cpp>
    composition-move.py check    <function> <baseline.cpp> <installed.cpp>
    composition-move.py adapter
"""
from __future__ import annotations

import importlib.util
import pathlib
import re
import sys

HERE = pathlib.Path(__file__).resolve().parent
HEADER = HERE.parent / "src/include/SisProdComposition.h"
STATE = "CompositionState"
NAMESPACE = "sisprod::composition"
ADAPTER = "compositionStateOf"

_spec = importlib.util.spec_from_file_location("transient_move", HERE / "transient-move.py")
tm = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(tm)

FUNCTIONS = {
    # T102d -- the three short routines, already under 200 lines.
    "renovaalbetini": {"new_name": "storePreviousFractionsAndMovePigs", "arguments": ""},
    "renovaMasEsp": {"new_name": "cacheCellAndFaceDensities", "arguments": ""},
    "avaliaParafina": {"new_name": "evaluateWaxDeposition", "arguments": ""},
}

# The composition routines do not call one another -- measured with
# measure-members.py, which reports corrDeng as the only method any of them
# calls.
CALLS: dict[str, tuple[str, str]] = {}

CALLBACKS = {
    "corrDeng": "state.updaters.correctGasSpecificGravity",
}


def header_fields() -> list[tuple[str, str]]:
    """(field, member) in declaration order, read from the header."""
    text = HEADER.read_text(encoding="utf-8")
    start = text.index(f"struct {STATE} {{")
    body = text[start:text.index("};", start)]
    fields = re.findall(r'^\s+[\w:<>]+\s*[*&]+\s*(\w+);\s*//\s*(\w+)\s*$', body, re.M)
    if not fields:
        raise SystemExit(f"no `T &field; // member` lines in {STATE}")
    return fields


MEMBERS = {member: f"state.{field}" for field, member in header_fields()}
FIELD_TO_MEMBER = {field: member for member, field in MEMBERS.items()}
FIELD_RE = re.compile(
    "(?<![\\w.])(" + "|".join(re.escape(f) for f in sorted(FIELD_TO_MEMBER, key=len, reverse=True)) + r")\b")
SIGNATURE_RE = tm.SIGNATURE_RE


def forward(old_name: str, body: str) -> str:
    spec = FUNCTIONS[old_name]
    new_name = spec["new_name"]
    match = SIGNATURE_RE.match(body)
    if match is None:
        raise ValueError(f"could not parse the signature of {old_name}")
    return_type, signature = match.group(1), match.group(3)
    separator = ", " if signature.strip() else ""
    head = f"{return_type}{new_name}(const {STATE} &state{separator}{signature}) {{"
    body = SIGNATURE_RE.sub(lambda _: head, body, count=1)
    shadowed = tm._declared_in_body(body)
    members = {k: v for k, v in MEMBERS.items() if k not in shadowed}
    pattern = re.compile(r"(?<![\w.])(" + "|".join(sorted(members, key=len, reverse=True)) + r")\b")
    body = tm.substitute_outside_comments(pattern, lambda m: members[m.group(1)], body)
    for member, adapter in CALLBACKS.items():
        body = tm.substitute_outside_comments(
            re.compile(rf"(?<![\w.>]){member}\("), f"{adapter}(", body)
    return body


def inverse(old_name: str, body: str) -> str:
    new_name = FUNCTIONS[old_name]["new_name"]
    for member, adapter in CALLBACKS.items():
        body = re.sub(rf"(?<![\w.>]){re.escape(adapter)}\(", f"{member}(", body)
    body = FIELD_RE.sub(lambda m: FIELD_TO_MEMBER[m.group(1)], body)
    pattern = re.compile(
        rf"^([\w:<>*&]+(?:\s+[\w:<>*&]+)*\s+){new_name}\(const {STATE} &state(?:, (.*?))?\)\s*\{{", re.S)
    match = pattern.match(body)
    if match is None:
        raise ValueError(f"could not parse the moved signature of {new_name}")
    return_type, signature = match.group(1), match.group(2) or ""
    return pattern.sub(lambda _: f"{return_type}SProd::{old_name}({signature}) {{", body, count=1)


def require(path: str, old_name: str):
    if old_name not in FUNCTIONS:
        raise ValueError(f"unsupported function: {old_name}")
    carved = tm.carve(open(path, encoding="utf-8").read(), old_name)
    if carved is None:
        raise ValueError(f"expected one SProd::{old_name} body, found none")
    return carved


def adapter() -> str:
    rows = [f"        .{field} = system.{member}," for field, member in header_fields()]
    return "\n".join([
        f"{NAMESPACE}::{STATE} {ADAPTER}(SProd &system) {{",
        f"    return {NAMESPACE}::{STATE}{{",
        *rows,
        "        .updaters = {system},",
        "    };",
        "}",
    ])


def main() -> int:
    if len(sys.argv) == 2 and sys.argv[1] == "adapter":
        print(adapter())
        return 0
    if len(sys.argv) < 5:
        print(__doc__, file=sys.stderr)
        return 2
    mode, name = sys.argv[1], sys.argv[2]

    if mode == "extract":
        start, end, body = require(sys.argv[3], name)
        open(sys.argv[4], "w", encoding="utf-8").write(forward(name, body) + "\n")
        print(f"{name:22} {start + 1}-{end + 1}")
        return 0

    if mode == "install" and len(sys.argv) == 6:
        _, _, body = require(sys.argv[3], name)
        target = open(sys.argv[4], encoding="utf-8").read()
        closing = f"}}  // namespace {NAMESPACE}\n"
        if target.count(closing) != 1:
            raise ValueError("expected one namespace closing marker")
        open(sys.argv[5], "w", encoding="utf-8").write(
            target.replace(closing, forward(name, body) + "\n\n" + closing, 1))
        print(f"installed {name}")
        return 0

    if mode == "delegate":
        source = open(sys.argv[3], encoding="utf-8").read()
        start, end, body = require(sys.argv[3], name)
        match = SIGNATURE_RE.match(body)
        return_type, signature = match.group(1), match.group(3)
        spec = FUNCTIONS[name]
        separator = ", " if spec["arguments"] else ""
        lead = "return " if return_type.strip() != "void" else ""
        call = f"{NAMESPACE}::{spec['new_name']}({ADAPTER}(*this){separator}{spec['arguments']})"
        lines = source.split("\n")
        lines[start:end + 1] = [f"{return_type}SProd::{name}({signature}) {{", f"    {lead}{call};", "}"]
        open(sys.argv[4], "w", encoding="utf-8").write("\n".join(lines))
        print(f"delegated {name} lines {start + 1}-{end + 1}")
        return 0

    if mode == "check":
        baseline = tm.carve(open(sys.argv[3], encoding="utf-8").read(), name)
        new_name = FUNCTIONS[name]["new_name"]
        lines = open(sys.argv[4], encoding="utf-8").read().split("\n")
        start = next((i for i, line in enumerate(lines)
                      if re.match(rf"^[\w:<>*&]+[\w:<>*&\s]*\b{new_name}\(const {STATE} &state", line)), None)
        if baseline is None or start is None:
            print(f"MISSING  {new_name} <- {name}")
            return 1
        brace = start
        while "{" not in lines[brace]:
            brace += 1
        depth, end = 0, brace
        while True:
            depth += lines[end].count("{") - lines[end].count("}")
            if depth == 0:
                break
            end += 1
        expected = tm.tokenize(baseline[2])
        actual = tm.tokenize(inverse(name, "\n".join(lines[start:end + 1])))
        if expected == actual:
            print(f"OK       {new_name} <- {name} ({len(expected)} tokens)")
            return 0
        position = next((i for i, pair in enumerate(zip(expected, actual)) if pair[0] != pair[1]),
                        min(len(expected), len(actual)))
        print(f"DIFFERS  {new_name} <- {name} at token {position}")
        print(f"  baseline: {' '.join(expected[max(0, position - 5):position + 6])}")
        print(f"  current : {' '.join(actual[max(0, position - 5):position + 6])}")
        return 1

    print(__doc__, file=sys.stderr)
    return 2


if __name__ == "__main__":
    sys.exit(main())
