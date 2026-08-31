#!/usr/bin/env python3
"""Extract the first Stage 5 thermal functions and prove the move literal.

Usage:
    thermal-move.py extract <sisprod.cpp> <output-fragment.cpp>
    thermal-move.py install <sisprod.cpp> <SisProdThermal.cpp>
    thermal-move.py delegate <sisprod.cpp> <rewritten-sisprod.cpp>
    thermal-move.py check   <baseline-sisprod.cpp> <extracted.cpp>
    thermal-move.py extract-calctemp <sisprod.cpp> <output-fragment.cpp>
    thermal-move.py install-calctemp <sisprod.cpp> <thermal-in.cpp> <thermal-out.cpp>
    thermal-move.py delegate-calctemp <sisprod.cpp> <rewritten-sisprod.cpp>
    thermal-move.py check-calctemp <baseline-sisprod.cpp> <extracted.cpp>
    thermal-move.py extract-t063 <function> <sisprod.cpp> <fragment.cpp>
    thermal-move.py install-t063 <function> <sisprod.cpp> <thermal-in.cpp> <thermal-out.cpp>
    thermal-move.py delegate-t063 <function> <sisprod.cpp> <rewritten-sisprod.cpp>
    thermal-move.py check-t063 <function> <baseline-sisprod.cpp> <extracted.cpp>
    thermal-move.py check-t063-decomposition <baseline-sisprod.cpp> <extracted.cpp>
    thermal-move.py extract-renova-temp <sisprod.cpp> <fragment.cpp>
    thermal-move.py install-renova-temp <sisprod.cpp> <thermal-in.cpp> <thermal-out.cpp>
    thermal-move.py delegate-renova-temp <sisprod.cpp> <rewritten-sisprod.cpp>
    thermal-move.py check-renova-temp <baseline-sisprod.cpp> <extracted.cpp>
    thermal-move.py check-renova-temp-decomposition <baseline-sisprod.cpp> <extracted.cpp>
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

CALCTEMP_MEMBERS = {
    "celula": "state.cells",
    "arq": "state.input",
    "vg1dSP": "state.globals",
    "semTermo": "state.thermalSourceDisabled",
    "verificaAcopRedeS": "state.productionNetworkCoupled",
    "SecPrimIniRedeP": "state.primarySectionStart",
    "SecPrimFimRedeP": "state.primarySectionEnd",
    "acertaIndAcop": "state.coupledCellIndices",
    "poisson3D": "state.poissonSolver",
    "ncel": "state.lastCell",
    "chokeSup": "state.surfaceChoke",
    "masChkSup": "state.surfaceChokeMassCondition",
    "noextremo": "state.networkEndpoint",
    "tGSup": "state.gasSurfaceTemperature",
    "CalcLat": "state.latentHeatEnabled",
}
CALCTEMP_FIELD_TO_MEMBER = {
    field.split(".", 1)[1]: member
    for member, field in CALCTEMP_MEMBERS.items()
}
CALCTEMP_MEMBER_RE = re.compile(
    r"\b(" + "|".join(CALCTEMP_MEMBERS) + r")\b"
)
CALCTEMP_FIELD_RE = re.compile(
    r"\bstate\.(" + "|".join(CALCTEMP_FIELD_TO_MEMBER) + r")\b"
)

T063_FUNCTIONS = {
    "calcTempEntalp": {
        "new_name": "updateTemperatureFromEnthalpy",
        "arguments": "i",
        "calls": {
            "calcHmix": "computeMixtureEnthalpy",
            "energmix": "interpolateMixtureEnergy",
        },
    },
    "calcTransMassTermo": {
        "new_name": "computeThermalMassTransfer",
        "arguments": "i",
        "calls": {
            "interpolaHLatente": "interpolateLatentHeat",
        },
    },
}

RENOVA_TEMP_MEMBERS = {
    **CALCTEMP_MEMBERS,
    "modeloCompleto": "state.completeModel",
    "TransMassModel": "state.massTransferModel",
}
RENOVA_TEMP_FIELD_TO_MEMBER = {
    field.split(".", 1)[1]: member
    for member, field in RENOVA_TEMP_MEMBERS.items()
}
RENOVA_TEMP_MEMBER_RE = re.compile(
    r"\b(" + "|".join(RENOVA_TEMP_MEMBERS) + r")\b"
)
RENOVA_TEMP_FIELD_RE = re.compile(
    r"\bstate\.(" + "|".join(RENOVA_TEMP_FIELD_TO_MEMBER) + r")\b"
)

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


def carve_calctemp(source: str) -> tuple[int, int, str] | None:
    lines = source.split("\n")
    for start, line in enumerate(lines):
        if not re.match(r"void SProd::calctemp\(", line):
            continue
        depth = 0
        opened = False
        for end in range(start, len(lines)):
            depth += lines[end].count("{") - lines[end].count("}")
            opened = opened or "{" in lines[end]
            if opened and depth == 0:
                return start, end, "\n".join(lines[start:end + 1])
    return None


def carve_new_calctemp(source: str) -> str | None:
    lines = source.split("\n")
    for start, line in enumerate(lines):
        if not re.match(r"void computeTemperature\(const ThermalState &state", line):
            continue
        depth = 0
        opened = False
        for end in range(start, len(lines)):
            depth += lines[end].count("{") - lines[end].count("}")
            opened = opened or "{" in lines[end]
            if opened and depth == 0:
                return "\n".join(lines[start:end + 1])
    return None


def carve_named(source: str, signature: str) -> tuple[int, int, str] | None:
    lines = source.split("\n")
    for start, line in enumerate(lines):
        if not re.match(signature, line):
            continue
        depth = 0
        opened = False
        for end in range(start, len(lines)):
            depth += lines[end].count("{") - lines[end].count("}")
            opened = opened or "{" in lines[end]
            if opened and depth == 0:
                return start, end, "\n".join(lines[start:end + 1])
    return None


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


def forward_calctemp(body: str) -> str:
    body = re.sub(
        r"^void SProd::calctemp\((.*?)\) \{",
        r"void computeTemperature(const ThermalState &state, \1) {",
        body,
        count=1,
    )
    body = CALCTEMP_MEMBER_RE.sub(
        lambda match: CALCTEMP_MEMBERS[match.group(1)], body
    )
    return re.sub(
        r"\binterpolaHLatente\(", "interpolateLatentHeat(state, ", body
    )


def inverse_calctemp(body: str) -> str:
    body = re.sub(
        r"\binterpolateLatentHeat\(state, ", "interpolaHLatente(", body
    )
    body = CALCTEMP_FIELD_RE.sub(
        lambda match: CALCTEMP_FIELD_TO_MEMBER[match.group(1)], body
    )
    return re.sub(
        r"^void computeTemperature\(const ThermalState &state, (.*?)\) \{",
        r"void SProd::calctemp(\1) {",
        body,
        count=1,
    )


def forward_t063(old_name: str, body: str) -> str:
    spec = T063_FUNCTIONS[old_name]
    new_name = spec["new_name"]
    body = re.sub(
        rf"^void SProd::{old_name}\((.*?)\) \{{",
        rf"void {new_name}(const ThermalState &state, \1) {{",
        body,
        count=1,
    )
    body = CALCTEMP_MEMBER_RE.sub(
        lambda match: CALCTEMP_MEMBERS[match.group(1)], body
    )
    for old_call, new_call in spec["calls"].items():
        body = re.sub(rf"\b{old_call}\(", f"{new_call}(state, ", body)
    return body


def inverse_t063(old_name: str, body: str) -> str:
    spec = T063_FUNCTIONS[old_name]
    new_name = spec["new_name"]
    for old_call, new_call in spec["calls"].items():
        body = re.sub(rf"\b{new_call}\(state, ", f"{old_call}(", body)
    body = CALCTEMP_FIELD_RE.sub(
        lambda match: CALCTEMP_FIELD_TO_MEMBER[match.group(1)], body
    )
    return re.sub(
        rf"^void {new_name}\(const ThermalState &state, (.*?)\) \{{",
        rf"void SProd::{old_name}(\1) {{",
        body,
        count=1,
    )


def forward_renova_temp(body: str) -> str:
    body = re.sub(
        r"^void SProd::renovaTemp\(\) \{",
        r"void updateDistributedMassTransfer(const ThermalState &state) {",
        body,
        count=1,
    )
    body = RENOVA_TEMP_MEMBER_RE.sub(
        lambda match: RENOVA_TEMP_MEMBERS[match.group(1)], body
    )
    return re.sub(r"\brenovaFonte\(", "state.sourceUpdater(", body)


def inverse_renova_temp(body: str) -> str:
    body = re.sub(r"\bstate\.sourceUpdater\(", "renovaFonte(", body)
    body = RENOVA_TEMP_FIELD_RE.sub(
        lambda match: RENOVA_TEMP_FIELD_TO_MEMBER[match.group(1)], body
    )
    return re.sub(
        r"^void updateDistributedMassTransfer\(const ThermalState &state\) \{",
        r"void SProd::renovaTemp() {",
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


def require_calctemp(source_path: str) -> tuple[int, int, str]:
    body = carve_calctemp(open(source_path, encoding="utf-8").read())
    if body is None:
        raise ValueError("expected one SProd::calctemp body, found none")
    return body


def extract_calctemp(source_path: str, output_path: str) -> int:
    start, end, body = require_calctemp(source_path)
    open(output_path, "w", encoding="utf-8").write(forward_calctemp(body) + "\n")
    print(f"calctemp               {start + 1}-{end + 1}")
    print(f"extracted calctemp into {output_path}")
    return 0


def install_calctemp(
    source_path: str, thermal_input_path: str, thermal_output_path: str
) -> int:
    _, _, body = require_calctemp(source_path)
    thermal = open(thermal_input_path, encoding="utf-8").read()
    closing = "}  // namespace sisprod::thermal\n"
    if thermal.count(closing) != 1:
        raise ValueError("expected one thermal namespace closing marker")
    replacement = forward_calctemp(body) + "\n\n" + closing
    open(thermal_output_path, "w", encoding="utf-8").write(
        thermal.replace(closing, replacement, 1)
    )
    print(f"installed calctemp into {thermal_output_path}")
    return 0


def delegate_calctemp(source_path: str, output_path: str) -> int:
    source = open(source_path, encoding="utf-8").read()
    start, end, body = require_calctemp(source_path)
    signature = re.match(r"void SProd::calctemp\((.*?)\) \{", body).group(1)
    wrapper = [
        f"void SProd::calctemp({signature}) {{",
        "    sisprod::thermal::computeTemperature(",
        "        thermalStateOf(*this), i, tempantiga, modoPerm);",
        "}",
    ]
    lines = source.split("\n")
    lines[start:end + 1] = wrapper
    open(output_path, "w", encoding="utf-8").write("\n".join(lines))
    print(f"delegated calctemp lines {start + 1}-{end + 1}")
    return 0


def require_t063(source_path: str, old_name: str) -> tuple[int, int, str]:
    if old_name not in T063_FUNCTIONS:
        raise ValueError(f"unsupported T063 function: {old_name}")
    source = open(source_path, encoding="utf-8").read()
    body = carve_named(source, rf"void SProd::{old_name}\(")
    if body is None:
        raise ValueError(f"expected one SProd::{old_name} body, found none")
    return body


def extract_t063(old_name: str, source_path: str, output_path: str) -> int:
    start, end, body = require_t063(source_path, old_name)
    open(output_path, "w", encoding="utf-8").write(
        forward_t063(old_name, body) + "\n"
    )
    print(f"{old_name:22} {start + 1}-{end + 1}")
    print(f"extracted {old_name} into {output_path}")
    return 0


def install_t063(
    old_name: str, source_path: str, thermal_input_path: str,
    thermal_output_path: str
) -> int:
    _, _, body = require_t063(source_path, old_name)
    thermal = open(thermal_input_path, encoding="utf-8").read()
    closing = "}  // namespace sisprod::thermal\n"
    if thermal.count(closing) != 1:
        raise ValueError("expected one thermal namespace closing marker")
    replacement = forward_t063(old_name, body) + "\n\n" + closing
    open(thermal_output_path, "w", encoding="utf-8").write(
        thermal.replace(closing, replacement, 1)
    )
    print(f"installed {old_name} into {thermal_output_path}")
    return 0


def delegate_t063(old_name: str, source_path: str, output_path: str) -> int:
    source = open(source_path, encoding="utf-8").read()
    start, end, body = require_t063(source_path, old_name)
    signature = re.match(
        rf"void SProd::{old_name}\((.*?)\) \{{", body
    ).group(1)
    spec = T063_FUNCTIONS[old_name]
    wrapper = [
        f"void SProd::{old_name}({signature}) {{",
        f"    sisprod::thermal::{spec['new_name']}(",
        f"        thermalStateOf(*this), {spec['arguments']});",
        "}",
    ]
    lines = source.split("\n")
    lines[start:end + 1] = wrapper
    open(output_path, "w", encoding="utf-8").write("\n".join(lines))
    print(f"delegated {old_name} lines {start + 1}-{end + 1}")
    return 0


def check_t063(old_name: str, baseline_path: str, current_path: str) -> int:
    baseline = carve_named(
        open(baseline_path, encoding="utf-8").read(),
        rf"void SProd::{old_name}\(",
    )
    new_name = T063_FUNCTIONS[old_name]["new_name"]
    current = carve_named(
        open(current_path, encoding="utf-8").read(),
        rf"void {new_name}\(const ThermalState &state",
    )
    if baseline is None or current is None:
        print(f"MISSING  {new_name} <- {old_name}")
        return 1
    expected = tokenize(baseline[2])
    actual = tokenize(inverse_t063(old_name, current[2]))
    if expected == actual:
        print(f"OK       {new_name} <- {old_name} ({len(expected)} tokens)")
        return 0
    position = first_difference(expected, actual)
    print(f"DIFFERS  {new_name} <- {old_name} at token {position}")
    print(f"  baseline: {' '.join(expected[max(0, position - 5):position + 6])}")
    print(f"  current : {' '.join(actual[max(0, position - 5):position + 6])}")
    return 1


def thermal_mass_source_block(body: str) -> str:
    start_marker = "    double fontemassG = 0.;"
    end_marker = "        fontemassG = 0;"
    start = body.find(start_marker)
    end = body.find(end_marker, start)
    if start < 0 or end < 0:
        raise ValueError("thermal mass-source block not found")
    return body[start:end + len(end_marker)]


def check_t063_decomposition(baseline_path: str, current_path: str) -> int:
    old_name = "calcTransMassTermo"
    baseline = require_t063(baseline_path, old_name)
    expected_main = forward_t063(old_name, baseline[2])
    expected_block = thermal_mass_source_block(expected_main)

    current_source = open(current_path, encoding="utf-8").read()
    current_main = carve_named(
        current_source,
        r"void computeThermalMassTransfer\(const ThermalState &state",
    )
    current_helper = carve_named(
        current_source,
        r"TemperatureSourceTerms computeThermalMassTransferSourceTerms\(",
    )
    if current_main is None or current_helper is None:
        print("MISSING  decomposed thermal mass-transfer function or helper")
        return 1

    actual_block = thermal_mass_source_block(current_helper[2])
    replacement = (
        "    TemperatureSourceTerms sourceTerms =\n"
        "        computeThermalMassTransferSourceTerms(state, i, dx);"
    )
    expected_main = expected_main.replace(expected_block, replacement, 1)
    expected_main = expected_main.replace(
        "fontemassL + fontemassG", "sourceTerms.liquid + sourceTerms.gas", 1
    )

    expected_main_tokens = tokenize(expected_main)
    actual_main_tokens = tokenize(current_main[2])
    expected_block_tokens = tokenize(expected_block)
    actual_block_tokens = tokenize(actual_block)

    failures = 0
    for label, expected, actual in (
        ("main", expected_main_tokens, actual_main_tokens),
        ("source helper", expected_block_tokens, actual_block_tokens),
    ):
        if expected == actual:
            print(f"OK       thermal mass-transfer {label} ({len(expected)} tokens)")
            continue
        failures += 1
        position = first_difference(expected, actual)
        print(f"DIFFERS  thermal mass-transfer {label} at token {position}")
        print(f"  expected: {' '.join(expected[max(0, position - 5):position + 6])}")
        print(f"  current : {' '.join(actual[max(0, position - 5):position + 6])}")
    return 1 if failures else 0


def require_renova_temp(source_path: str) -> tuple[int, int, str]:
    source = open(source_path, encoding="utf-8").read()
    body = carve_named(source, r"void SProd::renovaTemp\(\)")
    if body is None:
        raise ValueError("expected one SProd::renovaTemp body, found none")
    return body


def extract_renova_temp(source_path: str, output_path: str) -> int:
    start, end, body = require_renova_temp(source_path)
    open(output_path, "w", encoding="utf-8").write(
        forward_renova_temp(body) + "\n"
    )
    print(f"renovaTemp             {start + 1}-{end + 1}")
    print(f"extracted renovaTemp into {output_path}")
    return 0


def install_renova_temp(
    source_path: str, thermal_input_path: str, thermal_output_path: str
) -> int:
    _, _, body = require_renova_temp(source_path)
    thermal = open(thermal_input_path, encoding="utf-8").read()
    closing = "}  // namespace sisprod::thermal\n"
    if thermal.count(closing) != 1:
        raise ValueError("expected one thermal namespace closing marker")
    replacement = forward_renova_temp(body) + "\n\n" + closing
    open(thermal_output_path, "w", encoding="utf-8").write(
        thermal.replace(closing, replacement, 1)
    )
    print(f"installed renovaTemp into {thermal_output_path}")
    return 0


def delegate_renova_temp(source_path: str, output_path: str) -> int:
    source = open(source_path, encoding="utf-8").read()
    start, end, _ = require_renova_temp(source_path)
    wrapper = [
        "void SProd::renovaTemp() {",
        "    sisprod::thermal::updateDistributedMassTransfer(",
        "        thermalStateOf(*this));",
        "}",
    ]
    lines = source.split("\n")
    lines[start:end + 1] = wrapper
    open(output_path, "w", encoding="utf-8").write("\n".join(lines))
    print(f"delegated renovaTemp lines {start + 1}-{end + 1}")
    return 0


def check_renova_temp(baseline_path: str, current_path: str) -> int:
    baseline = carve_named(
        open(baseline_path, encoding="utf-8").read(),
        r"void SProd::renovaTemp\(\)",
    )
    current = carve_named(
        open(current_path, encoding="utf-8").read(),
        r"void updateDistributedMassTransfer\(const ThermalState &state\)",
    )
    if baseline is None or current is None:
        print("MISSING  updateDistributedMassTransfer <- renovaTemp")
        return 1
    expected = tokenize(baseline[2])
    actual = tokenize(inverse_renova_temp(current[2]))
    if expected == actual:
        print(
            "OK       updateDistributedMassTransfer <- renovaTemp "
            f"({len(expected)} tokens)"
        )
        return 0
    position = first_difference(expected, actual)
    print(
        "DIFFERS  updateDistributedMassTransfer <- renovaTemp "
        f"at token {position}"
    )
    print(f"  baseline: {' '.join(expected[max(0, position - 5):position + 6])}")
    print(f"  current : {' '.join(actual[max(0, position - 5):position + 6])}")
    return 1


def text_between(text: str, start_marker: str, end_marker: str) -> str:
    start = text.find(start_marker)
    end = text.find(end_marker, start + len(start_marker))
    if start < 0 or end < 0:
        raise ValueError(
            f"expected block between {start_marker!r} and {end_marker!r}"
        )
    return text[start:end]


def function_core(body: str) -> str:
    opening = body.find("{")
    closing = body.rfind("}")
    if opening < 0 or closing < opening:
        raise ValueError("function braces not found")
    return body[opening + 1:closing]


def compare_token_blocks(
    blocks: tuple[tuple[str, str, str], ...]
) -> int:
    failures = 0
    for label, expected_text, actual_text in blocks:
        expected = tokenize(expected_text)
        actual = tokenize(actual_text)
        if expected == actual:
            print(f"OK       renovaTemp {label} ({len(expected)} tokens)")
            continue
        failures += 1
        position = first_difference(expected, actual)
        print(f"DIFFERS  renovaTemp {label} at token {position}")
        print(
            f"  expected: "
            f"{' '.join(expected[max(0, position - 5):position + 6])}"
        )
        print(
            f"  current : "
            f"{' '.join(actual[max(0, position - 5):position + 6])}"
        )
    return failures


def check_renova_temp_decomposition(
    baseline_path: str, current_path: str
) -> int:
    baseline = require_renova_temp(baseline_path)
    expected_main = forward_renova_temp(baseline[2])
    current_source = open(current_path, encoding="utf-8").read()

    signatures = {
        "main": r"void updateDistributedMassTransfer\(const ThermalState &state\)",
        "inlet": r"void initializeDistributedMassTransferInlet\(",
        "properties": (
            r"DistributedMassTransferProperties "
            r"prepareDistributedMassTransferProperties\("
        ),
        "derivatives": (
            r"DistributedMassTransferCoefficients "
            r"updateDistributedMassTransferDerivatives\("
        ),
        "selection": r"void selectDistributedMassTransferModel\(",
        "application": r"void applyDistributedMassTransferModel\(",
    }
    bodies = {
        name: carve_named(current_source, signature)
        for name, signature in signatures.items()
    }
    missing = [name for name, body in bodies.items() if body is None]
    if missing:
        print(f"MISSING  renovaTemp decomposition: {', '.join(missing)}")
        return 1

    property_start = "            double fwd;"
    derivative_start = "            double ativa = 1.;"
    selection_start = "            tmed = state.cells[i - 1].temp;"
    tail_start = "            if (state.cells[i - 1].TMModel == -2) {"
    inlet_start = "            state.cells[0].transmassLini"
    inlet_end = "            ProFlu flutemp = state.cells[i].flui;"

    expected_properties = text_between(
        expected_main, property_start, derivative_start
    )
    expected_derivatives = text_between(
        expected_main, derivative_start, selection_start
    )
    selection_position = expected_main.find(
        selection_start, expected_main.find(derivative_start)
    )
    tail_position = expected_main.find(tail_start, selection_position)
    if selection_position < 0 or tail_position < 0:
        raise ValueError("renovaTemp selection/application block not found")
    expected_selection_and_application = expected_main[
        selection_position:tail_position
    ]
    application_start = "            state.cells[i - 1].fontedissolv = 0.;"
    expected_selection = text_between(
        expected_selection_and_application,
        selection_start,
        application_start,
    )
    expected_application = expected_selection_and_application[
        expected_selection_and_application.find(application_start):
    ]
    inlet_body_end = "\n        }\n    }\n}"
    expected_inlet = text_between(
        expected_main, inlet_start, inlet_body_end
    )

    inlet_core = function_core(bodies["inlet"][2])
    for new_name, old_name in (
        ("previousLiquidDensity", "rhol0"),
        ("previousOilVolumeFactor", "boL"),
        ("previousSolutionGasRatio", "rsL"),
        ("previousSolutionGasPressureDerivative", "DRsBoL"),
        ("cellIndex", "i"),
    ):
        inlet_core = re.sub(rf"\b{new_name}\b", old_name, inlet_core)

    properties_core = text_between(
        function_core(bodies["properties"][2]),
        "    double fwd;",
        "\n\n    return DistributedMassTransferProperties{",
    )
    derivatives_core = text_between(
        function_core(bodies["derivatives"][2]),
        "    double ativa = 1.;",
        "\n\n    return DistributedMassTransferCoefficients{",
    )
    selection_core = function_core(bodies["selection"][2])
    application_core = function_core(bodies["application"][2])
    application_core = application_core[
        application_core.find("    state.cells[i - 1].fontedissolv = 0.;"):
    ]

    property_plumbing = """            DistributedMassTransferProperties properties =
                prepareDistributedMassTransferProperties(
                    state, i, tmed, flue, flud);
            double fwd = properties.downstreamWaterFraction;
            double fwe = properties.upstreamWaterFraction;
            double fwC = properties.cellWaterFraction;
            double rl = properties.liquidDensity;
            double rg = properties.gasDensity;
            double betI = properties.downstreamComposition;
            double betL = properties.upstreamComposition;
            double rhol = properties.mixtureLiquidDensity;
            double boR = properties.downstreamOilVolumeFactor;
            double rsR = properties.downstreamSolutionGasRatio;
            double DRsBoR =
                properties.downstreamSolutionGasPressureDerivative;
            double boM = properties.cellOilVolumeFactor;
            double rsM = properties.cellSolutionGasRatio;
            double DRsBoM = properties.cellSolutionGasPressureDerivative;
            double DRsBoMT =
                properties.cellSolutionGasTemperatureDerivative;
"""
    derivative_plumbing = """            DistributedMassTransferCoefficients coefficients =
                updateDistributedMassTransferDerivatives(
                    state, i, fwC, flud, DRsBoM, DRsBoMT);
            double ativa = coefficients.activeDerivative;
            double acop = coefficients.spatialCoupling;
            double A1 = coefficients.flowArea;

"""
    application_plumbing = """            applyDistributedMassTransferModel(
                state, i, tmed, tmedL, ABSjL, flue, flud, fwd, fwe,
                fwC, betI, betL, rhol, boR, rsR, DRsBoR, boM, rsM,
                ativa, acop, A1, rhol0, boL, rsL, DRsBoL);
"""
    inlet_branch = text_between(
        expected_main,
        "        } else if (i == 0) {",
        inlet_body_end,
    ) + "\n        }"
    inlet_plumbing = """        } else if (i == 0)
            initializeDistributedMassTransferInlet(
                state, i, rhol0, boL, rsL, DRsBoL);"""

    expected_decomposed_main = expected_main.replace(
        expected_properties, property_plumbing, 1
    ).replace(
        expected_derivatives, derivative_plumbing, 1
    ).replace(
        expected_selection_and_application, application_plumbing, 1
    ).replace(
        inlet_branch, inlet_plumbing, 1
    )

    failures = compare_token_blocks((
        ("main", expected_decomposed_main, bodies["main"][2]),
        ("inlet helper", expected_inlet, inlet_core),
        ("properties helper", expected_properties, properties_core),
        ("derivatives helper", expected_derivatives, derivatives_core),
        ("selection helper", expected_selection, selection_core),
        ("application helper", expected_application, application_core),
    ))
    return 1 if failures else 0


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


def check_calctemp(baseline_path: str, current_path: str) -> int:
    baseline = carve_calctemp(open(baseline_path, encoding="utf-8").read())
    current = carve_new_calctemp(open(current_path, encoding="utf-8").read())
    if baseline is None or current is None:
        print("MISSING  computeTemperature <- calctemp")
        return 1
    expected = tokenize(baseline[2])
    actual = tokenize(inverse_calctemp(current))
    if expected == actual:
        print(f"OK       computeTemperature <- calctemp ({len(expected)} tokens)")
        return 0
    position = first_difference(expected, actual)
    print(f"DIFFERS  computeTemperature <- calctemp at token {position}")
    print(f"  baseline: {' '.join(expected[max(0, position - 5):position + 6])}")
    print(f"  current : {' '.join(actual[max(0, position - 5):position + 6])}")
    return 1


def main() -> int:
    if len(sys.argv) not in {4, 5, 6}:
        print(__doc__, file=sys.stderr)
        return 2
    mode = sys.argv[1]
    try:
        if mode == "install-t063" and len(sys.argv) == 6:
            return install_t063(
                sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5]
            )
        if mode in {"extract-t063", "delegate-t063", "check-t063"} \
                and len(sys.argv) == 5:
            old_name, first, second = sys.argv[2:]
            if mode == "extract-t063":
                return extract_t063(old_name, first, second)
            if mode == "delegate-t063":
                return delegate_t063(old_name, first, second)
            return check_t063(old_name, first, second)
        if mode == "check-t063-decomposition" and len(sys.argv) == 4:
            return check_t063_decomposition(sys.argv[2], sys.argv[3])
        if mode == "install-renova-temp" and len(sys.argv) == 5:
            return install_renova_temp(sys.argv[2], sys.argv[3], sys.argv[4])
        if mode == "install-calctemp" and len(sys.argv) == 5:
            return install_calctemp(sys.argv[2], sys.argv[3], sys.argv[4])
        if len(sys.argv) != 4:
            print(__doc__, file=sys.stderr)
            return 2
        first, second = sys.argv[2:]
        if mode in {"extract", "install"}:
            return extract_or_install(mode, first, second)
        if mode == "delegate":
            return delegate(first, second)
        if mode == "check":
            return check(first, second)
        if mode == "extract-calctemp":
            return extract_calctemp(first, second)
        if mode == "delegate-calctemp":
            return delegate_calctemp(first, second)
        if mode == "check-calctemp":
            return check_calctemp(first, second)
        if mode == "extract-renova-temp":
            return extract_renova_temp(first, second)
        if mode == "delegate-renova-temp":
            return delegate_renova_temp(first, second)
        if mode == "check-renova-temp":
            return check_renova_temp(first, second)
        if mode == "check-renova-temp-decomposition":
            return check_renova_temp_decomposition(first, second)
    except (OSError, ValueError, AttributeError) as error:
        print(error, file=sys.stderr)
        return 2
    print(__doc__, file=sys.stderr)
    return 2


if __name__ == "__main__":
    sys.exit(main())
