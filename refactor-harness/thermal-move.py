#!/usr/bin/env python3
"""Extract the first Stage 5 thermal functions and prove the move literal.

Usage:
    thermal-move.py extract <sisprod.cpp> <output-fragment.cpp>
    thermal-move.py install <sisprod.cpp> <SisProdThermal.cpp>
    thermal-move.py delegate <sisprod.cpp> <rewritten-sisprod.cpp>
    thermal-move.py check   <baseline-sisprod.cpp> <extracted.cpp>
"""
from __future__ import annotations

import re
import sys


FUNCTIONS = {
    "interpolaHLatente": "interpolateLatentHeat",
    "calcHmix": "computeMixtureEnthalpy",
    "energmix": "interpolateMixtureEnergy",
}
REVERSE_FUNCTIONS = {new: old for old, new in FUNCTIONS.items()}
DELEGATE_ARGUMENTS = {
    "interpolaHLatente": "pres, temp",
    "calcHmix": "i",
    "energmix": "i, jp0, jt, razp",
}

MEMBERS = {
    "celula": "state.cells",
    "arq": "state.input",
    "HLat": "state.latentHeatTable",
}
FIELD_TO_MEMBER = {field.split(".", 1)[1]: member for member, field in MEMBERS.items()}
MEMBER_RE = re.compile(r"\b(" + "|".join(MEMBERS) + r")\b")
FIELD_RE = re.compile(r"\bstate\.(" + "|".join(FIELD_TO_MEMBER) + r")\b")

COMMENT = re.compile(r"//[^\n]*|/\*.*?\*/", re.DOTALL)


def scan_quoted(text: str, start: int) -> tuple[str, int]:
    quote = text[start]
    cursor = start + 1
    while cursor < len(text):
        if text[cursor] == "\\":
            cursor += 2
        elif text[cursor] == quote:
            return text[start:cursor + 1], cursor + 1
        else:
            cursor += 1
    return text[start:], len(text)


def scan_identifier(text: str, start: int) -> tuple[str, int]:
    cursor = start + 1
    while cursor < len(text) and (text[cursor].isalnum() or text[cursor] == "_"):
        cursor += 1
    return text[start:cursor], cursor


def scan_number(text: str, start: int) -> tuple[str, int]:
    cursor = start + 1
    while cursor < len(text):
        current = text[cursor]
        if current.isalnum() or current in "._":
            cursor += 1
            continue
        if current in "+-" and text[cursor - 1] in "eEpP":
            cursor += 1
            continue
        break
    return text[start:cursor], cursor


def tokenize(text: str) -> list[str]:
    clean = COMMENT.sub(" ", text)
    tokens: list[str] = []
    cursor = 0
    while cursor < len(clean):
        current = clean[cursor]
        if current.isspace():
            cursor += 1
        elif current in "\"'":
            token, cursor = scan_quoted(clean, cursor)
            tokens.append(token)
        elif current.isalpha() or current == "_":
            token, cursor = scan_identifier(clean, cursor)
            tokens.append(token)
        elif current.isdigit() or (current == "." and cursor + 1 < len(clean)
                                  and clean[cursor + 1].isdigit()):
            token, cursor = scan_number(clean, cursor)
            tokens.append(token)
        else:
            tokens.append(current)
            cursor += 1
    return tokens


def carve(source: str) -> dict[str, tuple[int, int, str]]:
    lines = source.split("\n")
    found: dict[str, tuple[int, int, str]] = {}
    names = "|".join(FUNCTIONS)
    for start, line in enumerate(lines):
        match = re.match(rf"double SProd::({names})\(", line)
        if not match:
            continue
        depth = 0
        opened = False
        for end in range(start, len(lines)):
            depth += lines[end].count("{") - lines[end].count("}")
            opened = opened or "{" in lines[end]
            if opened and depth == 0:
                found[match.group(1)] = (start, end, "\n".join(lines[start:end + 1]))
                break
    return found


def carve_new(source: str) -> dict[str, str]:
    lines = source.split("\n")
    found: dict[str, str] = {}
    names = "|".join(FUNCTIONS.values())
    for start, line in enumerate(lines):
        match = re.match(rf"double ({names})\(const ThermalState &state", line)
        if not match:
            continue
        depth = 0
        opened = False
        for end in range(start, len(lines)):
            depth += lines[end].count("{") - lines[end].count("}")
            opened = opened or "{" in lines[end]
            if opened and depth == 0:
                found[match.group(1)] = "\n".join(lines[start:end + 1])
                break
    return found


def forward(old_name: str, body: str) -> str:
    signature = re.match(rf"double SProd::{old_name}\((.*?)\) \{{", body).group(1)
    body = re.sub(
        rf"^double SProd::{old_name}\(.*?\) \{{",
        f"double {FUNCTIONS[old_name]}(const ThermalState &state, {signature}) {{",
        body,
        count=1,
    )
    return MEMBER_RE.sub(lambda match: MEMBERS[match.group(1)], body)


def inverse(new_name: str, body: str) -> str:
    body = FIELD_RE.sub(lambda match: FIELD_TO_MEMBER[match.group(1)], body)
    signature = re.match(
        rf"double {new_name}\(const ThermalState &state, (.*?)\) \{{", body
    ).group(1)
    return re.sub(
        rf"^double {new_name}\(.*?\) \{{",
        f"double SProd::{REVERSE_FUNCTIONS[new_name]}({signature}) {{",
        body,
        count=1,
    )


def require_bodies(source_path: str) -> dict[str, tuple[int, int, str]]:
    bodies = carve(open(source_path, encoding="utf-8").read())
    missing = [name for name in FUNCTIONS if name not in bodies]
    if missing or len(bodies) != len(FUNCTIONS):
        raise ValueError(
            f"expected {len(FUNCTIONS)} bodies, found {len(bodies)}; missing: {missing}"
        )
    return bodies


def extract_or_install(mode: str, source_path: str, output_path: str) -> int:
    bodies = require_bodies(source_path)
    parts = [forward(name, bodies[name][2]) for name in FUNCTIONS]
    content = "\n\n".join(parts) + "\n"
    if mode == "install":
        content = (
            '#include "SisProdThermal.h"\n\n'
            '#include "Leitura.h"\n'
            '#include "celula3.h"\n\n'
            '#include <math.h>\n\n'
            'namespace sisprod::thermal {\n\n'
            + content
            + '}  // namespace sisprod::thermal\n'
        )
    open(output_path, "w", encoding="utf-8").write(content)
    for name in FUNCTIONS:
        start, end, _ = bodies[name]
        print(f"{name:22} {start + 1}-{end + 1}")
    print(f"{mode}ed {len(parts)} function(s) into {output_path}")
    return 0


def delegate(source_path: str, output_path: str) -> int:
    source = open(source_path, encoding="utf-8").read()
    bodies = require_bodies(source_path)
    lines = source.split("\n")
    ordered_names = sorted(FUNCTIONS, key=lambda name: bodies[name][0], reverse=True)
    for old_name in ordered_names:
        start, end, body = bodies[old_name]
        signature = re.match(rf"double SProd::{old_name}\((.*?)\) \{{", body).group(1)
        wrapper = [
            f"double SProd::{old_name}({signature}) {{",
            f"    return sisprod::thermal::{FUNCTIONS[old_name]}(",
            f"        thermalStateOf(*this), {DELEGATE_ARGUMENTS[old_name]});",
            "}",
        ]
        lines[start:end + 1] = wrapper
        print(f"delegated {old_name:22} lines {start + 1}-{end + 1}")
    open(output_path, "w", encoding="utf-8").write("\n".join(lines))
    return 0


def first_difference(expected: list[str], actual: list[str]) -> int:
    return next(
        (index for index, pair in enumerate(zip(expected, actual)) if pair[0] != pair[1]),
        min(len(expected), len(actual)),
    )


def check(baseline_path: str, current_path: str) -> int:
    baseline = carve(open(baseline_path, encoding="utf-8").read())
    current = carve_new(open(current_path, encoding="utf-8").read())
    failures = 0
    total = 0
    for old_name, new_name in FUNCTIONS.items():
        if old_name not in baseline or new_name not in current:
            print(f"MISSING  {new_name} <- {old_name}")
            failures += 1
            continue
        expected = tokenize(baseline[old_name][2])
        actual = tokenize(inverse(new_name, current[new_name]))
        total += len(expected)
        if expected == actual:
            print(f"OK       {new_name:26} <- {old_name} ({len(expected)} tokens)")
            continue
        failures += 1
        position = first_difference(expected, actual)
        print(f"DIFFERS  {new_name:26} <- {old_name} at token {position}")
        print(f"  baseline: {' '.join(expected[max(0, position - 5):position + 6])}")
        print(f"  current : {' '.join(actual[max(0, position - 5):position + 6])}")
    print(f"compared {len(FUNCTIONS)} bodies, {total} tokens, {failures} failure(s)")
    return 1 if failures else 0


def main() -> int:
    if len(sys.argv) != 4:
        print(__doc__, file=sys.stderr)
        return 2
    mode, first, second = sys.argv[1:]
    try:
        if mode in {"extract", "install"}:
            return extract_or_install(mode, first, second)
        if mode == "delegate":
            return delegate(first, second)
        if mode == "check":
            return check(first, second)
    except (OSError, ValueError, AttributeError) as error:
        print(error, file=sys.stderr)
        return 2
    print(__doc__, file=sys.stderr)
    return 2


if __name__ == "__main__":
    sys.exit(main())