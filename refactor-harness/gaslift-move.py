#!/usr/bin/env python3
"""Move Stage 6 gas-line bodies out of SisProd.cpp and prove the move literal.

The method is the one thermal-move.py established: rewrite the signature and the
SProd member accesses through an INVERTIBLE substitution table, then prove the
move by applying the inverse to the installed body and comparing tokens with the
pre-move commit. Nothing is retyped, so no arithmetic can drift; and because the
table inverts, the comparison is exact rather than approximate.

That matters more here than usual: coverage measurement found the demo corpus
never executes twelve of the twenty-two routines this stage moves, so for those
the token proof and the dedicated harness are the only verification there is.

Usage:
    gaslift-move.py extract  <function> <sisprod.cpp> <fragment.cpp>
    gaslift-move.py install  <function> <sisprod.cpp> <gaslift-in.cpp> <gaslift-out.cpp>
    gaslift-move.py delegate <function> <sisprod.cpp> <rewritten-sisprod.cpp>
    gaslift-move.py check    <function> <baseline-sisprod.cpp> <installed.cpp>
"""
from __future__ import annotations

import re
import sys

FUNCTIONS = {
    "HidroDescargaG": {"new_name": "computeGasUnloadingHydrostatics", "arguments": ""},
    "renovaGas": {"new_name": "updateGasLine", "arguments": ""},
    "renovaGasBuf": {"new_name": "updateBufferedGasLine", "arguments": ""},
    # areaValvCali reads nothing but its eight arguments. It gets NO state
    # parameter: handing one to a pure function would be a lie about what it
    # touches, and the header already declares it without.
    "areaValvCali": {"new_name": "calibratedValveArea", "stateless": True,
                     "arguments": "PCal, TCal, PVO, PT, dextern, areagarg, Rvalv, Temp"},
    "prescordesc": {"new_name": "unloadingPressureCorrection",
                    "arguments": "vazmax, ivalv, fator, sinal"},
    # The first parameter is vazGarg, a flow rate compared against
    # arq.vazDescControl -- not a pressure. The table said "pres" and the
    # delegation did not compile, which is how that was found.
    "CalcPresValvDesc": {"new_name": "computeUnloadingValvePressure",
                         "arguments": "vazGarg, ivalv"},
    "resolveDescarga": {"new_name": "solveUnloading", "arguments": ""},
    "avancInter": {"new_name": "advanceInterface", "arguments": ""},
}

# Calls BETWEEN moved routines. The moved body must reach the namespace version,
# and a stateful callee needs `state` threaded through; a stateless one must not
# receive it.
CALLS = {
    "prescordesc": ("unloadingPressureCorrection", True),
    "CalcPresValvDesc": ("computeUnloadingValvePressure", True),
    "areaValvCali": ("calibratedValveArea", False),
    "HidroDescargaG": ("computeGasUnloadingHydrostatics", True),
    "renovaGas": ("updateGasLine", True),
    "renovaGasBuf": ("updateBufferedGasLine", True),
    "resolveDescarga": ("solveUnloading", True),
    "avancInter": ("advanceInterface", True),
}

# SProd member -> GasLiftState field. Longest first when the pattern is built, so
# a short name cannot claim the prefix of a longer one (ncel vs ncelGas).
# Calls the moved bodies make back into SProd, routed through an adapter the
# state carries. tempDescarga is still an SProd member; the gas line reaches the
# thermal module through it rather than growing a ThermalState of its own.
CALLBACKS = {
    "tempDescarga": "state.temperatureUpdater.dischargeTemperature",
}

MEMBERS = {
    "celulaG": "state.gasCells",
    "celula": "state.cells",
    "arq": "state.input",
    "vg1dSP": "state.globals",
    "ncelGas": "state.gasCellCount",
    "ncel": "state.lastCell",
    "termolivreG": "state.gasFreeTerms",
    # SProd::dt. Longest-first ordering matters: dtInter and dtDesc must claim
    # their names before the bare dt can.
    "dt": "state.timeStep",
    "celInter": "state.interfaceCell",
    "presiniG": "state.initialGasPressure",
    "tempiniG": "state.initialGasTemperature",
    "pGSup": "state.gasSurfacePressure",
    "posicVGLG": "state.gasValveCellIndices",
    "posicVGLP": "state.productionValveCellIndices",
    "chokeVGL": "state.gasLiftChokes",
    "chokeInj": "state.injectionChoke",
    "iterperm": "state.steadyIteration",
    "semTermo": "state.thermalSourceDisabled",
    "verificaAcop": "state.networkCoupled",
    "ColunaAnulaIni": "state.annulusTubingStart",
    "ColunaAnulaFim": "state.annulusTubingEnd",
    "AnulaColunaIni": "state.tubingAnnulusStart",
    "AnulaColunaFim": "state.tubingAnnulusEnd",
    "velInter": "state.interfaceVelocity",
    "dtInter": "state.interfaceTimeStep",
    "matglobG": "state.gasSystemMatrix",
}
FIELD_TO_MEMBER = {field.split(".", 1)[1]: member for member, field in MEMBERS.items()}
MEMBER_RE = re.compile(
    r"(?<![\w.])(" + "|".join(sorted(MEMBERS, key=len, reverse=True)) + r")\b"
)
FIELD_RE = re.compile(
    r"\bstate\.(" + "|".join(sorted(FIELD_TO_MEMBER, key=len, reverse=True)) + r")\b"
)
SIGNATURE_RE = re.compile(
    r"^([\w:<>*&]+(?:\s+[\w:<>*&]+)*\s+)SProd::(\w+)\((.*?)\)\s*\{", re.S
)

COMMENT = re.compile(r"//[^\n]*|/\*.*?\*/", re.S)


def tokenize(text: str) -> list[str]:
    """Whitespace- and comment-insensitive; everything else is a token.

    Numbers keep their exact spelling, so 1 and 1. are different tokens and a
    precomputed constant cannot slip through. Each operator character is its own
    token, so reassociation and a changed sign both show up.
    """
    clean = COMMENT.sub(" ", text)
    tokens: list[str] = []
    cursor = 0
    while cursor < len(clean):
        current = clean[cursor]
        if current.isspace():
            cursor += 1
        elif current in "\"'":
            quote = current
            end = cursor + 1
            while end < len(clean) and clean[end] != quote:
                end += 2 if clean[end] == "\\" else 1
            tokens.append(clean[cursor:end + 1])
            cursor = end + 1
        elif current.isalpha() or current == "_":
            end = cursor
            while end < len(clean) and (clean[end].isalnum() or clean[end] == "_"):
                end += 1
            tokens.append(clean[cursor:end])
            cursor = end
        elif current.isdigit() or (current == "." and cursor + 1 < len(clean)
                                   and clean[cursor + 1].isdigit()):
            end = cursor
            while end < len(clean) and (clean[end].isalnum() or clean[end] in "._"
                                        or (clean[end] in "+-" and clean[end - 1] in "eE")):
                end += 1
            tokens.append(clean[cursor:end])
            cursor = end
        else:
            tokens.append(current)
            cursor += 1
    return tokens


def substitute_outside_comments(pattern, replacement, text: str) -> str:
    """Rewrite only in code. Comments keep their words.

    Stage 5 learned this the hard way: a rename that reaches into prose turned
    "mudado para &&" into "mudado debugStop &&" in fifteen comments.
    """
    out = []
    last = 0
    for match in COMMENT.finditer(text):
        out.append(pattern.sub(replacement, text[last:match.start()]))
        out.append(match.group(0))
        last = match.end()
    out.append(pattern.sub(replacement, text[last:]))
    return "".join(out)


def carve(source: str, name: str):
    lines = source.split("\n")
    start = next((i for i, line in enumerate(lines)
                  if re.match(rf"^[\w:<>*&]+[\w:<>*&\s]*\bSProd::{name}\(", line)), None)
    if start is None:
        return None
    brace = start
    while "{" not in lines[brace]:
        brace += 1
    depth, end = 0, brace
    while True:
        depth += lines[end].count("{") - lines[end].count("}")
        if depth == 0:
            break
        end += 1
    return start, end, "\n".join(lines[start:end + 1])


def forward(old_name: str, body: str) -> str:
    spec = FUNCTIONS[old_name]
    new_name = spec["new_name"]
    stateless = spec.get("stateless", False)
    match = SIGNATURE_RE.match(body)
    if match is None:
        raise ValueError(f"could not parse the signature of {old_name}")
    return_type, _, signature = match.group(1), match.group(2), match.group(3)
    if stateless:
        head = f"{return_type}{new_name}({signature}) {{"
    else:
        separator = ", " if signature.strip() else ""
        head = (f"{return_type}{new_name}(const GasLiftState &state"
                f"{separator}{signature}) {{")
    body = SIGNATURE_RE.sub(lambda _: head, body, count=1)
    body = substitute_outside_comments(MEMBER_RE,
                                       lambda m: MEMBERS[m.group(1)], body)
    for called, (renamed, needs_state) in CALLS.items():
        if called == old_name:
            continue
        prefix = "state, " if needs_state else ""
        body = substitute_outside_comments(
            re.compile(rf"(?<![\w.>]){called}\("), f"{renamed}({prefix}", body)
    for member, adapter in CALLBACKS.items():
        body = substitute_outside_comments(
            re.compile(rf"(?<![\w.>]){member}\("), f"{adapter}(", body)
    return body


def inverse(old_name: str, body: str) -> str:
    spec = FUNCTIONS[old_name]
    new_name = spec["new_name"]
    stateless = spec.get("stateless", False)
    for called, (renamed, needs_state) in reversed(list(CALLS.items())):
        if called == old_name:
            continue
        prefix = "state, " if needs_state else ""
        body = re.sub(rf"(?<![\w.>]){re.escape(renamed)}\({re.escape(prefix)}",
                      f"{called}(", body)
    for member, adapter in CALLBACKS.items():
        body = re.sub(rf"(?<![\w.>]){re.escape(adapter)}\(", f"{member}(", body)
    body = FIELD_RE.sub(lambda m: FIELD_TO_MEMBER[m.group(1)], body)
    if stateless:
        pattern = re.compile(
            rf"^([\w:<>*&]+(?:\s+[\w:<>*&]+)*\s+){new_name}\((.*?)\)\s*\{{", re.S)
        match = pattern.match(body)
        if match is None:
            raise ValueError(f"could not parse the moved signature of {new_name}")
        return_type, signature = match.group(1), match.group(2)
    else:
        pattern = re.compile(
            rf"^([\w:<>*&]+(?:\s+[\w:<>*&]+)*\s+){new_name}"
            rf"\(const GasLiftState &state(?:, (.*?))?\)\s*\{{", re.S)
        match = pattern.match(body)
        if match is None:
            raise ValueError(f"could not parse the moved signature of {new_name}")
        return_type, signature = match.group(1), match.group(2) or ""
    return pattern.sub(
        lambda _: f"{return_type}SProd::{old_name}({signature}) {{",
        body, count=1)


def require(path: str, old_name: str):
    if old_name not in FUNCTIONS:
        raise ValueError(f"unsupported function: {old_name}")
    carved = carve(open(path, encoding="utf-8").read(), old_name)
    if carved is None:
        raise ValueError(f"expected one SProd::{old_name} body, found none")
    return carved


def main() -> int:
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
        closing = "}  // namespace sisprod::gaslift\n"
        if target.count(closing) != 1:
            raise ValueError("expected one gaslift namespace closing marker")
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
        arguments = spec["arguments"]
        lead = "return " if return_type.strip() != "void" else ""
        if spec.get("stateless", False):
            call = f"sisprod::gaslift::{spec['new_name']}({arguments})"
        else:
            separator = ", " if arguments else ""
            call = (f"sisprod::gaslift::{spec['new_name']}"
                    f"(gasLiftStateOf(*this){separator}{arguments})")
        wrapper = [
            f"{return_type}SProd::{name}({signature}) {{",
            f"    {lead}{call};",
            "}",
        ]
        lines = source.split("\n")
        lines[start:end + 1] = wrapper
        open(sys.argv[4], "w", encoding="utf-8").write("\n".join(lines))
        print(f"delegated {name} lines {start + 1}-{end + 1}")
        return 0

    if mode == "check":
        baseline = carve(open(sys.argv[3], encoding="utf-8").read(), name)
        new_name = FUNCTIONS[name]["new_name"]
        current_source = open(sys.argv[4], encoding="utf-8").read()
        lines = current_source.split("\n")
        opening = ("" if FUNCTIONS[name].get("stateless", False)
                   else r"\(const GasLiftState &state")
        start = next((i for i, line in enumerate(lines)
                      if re.match(rf"^[\w:<>*&]+[\w:<>*&\s]*\b{new_name}"
                                  rf"{opening or r'\('}", line)), None)
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
        expected = tokenize(baseline[2])
        actual = tokenize(inverse(name, "\n".join(lines[start:end + 1])))
        if expected == actual:
            print(f"OK       {new_name} <- {name} ({len(expected)} tokens)")
            return 0
        position = next((i for i, pair in enumerate(zip(expected, actual))
                         if pair[0] != pair[1]), min(len(expected), len(actual)))
        print(f"DIFFERS  {new_name} <- {name} at token {position}")
        print(f"  baseline: {' '.join(expected[max(0, position - 5):position + 6])}")
        print(f"  current : {' '.join(actual[max(0, position - 5):position + 6])}")
        return 1

    print(__doc__, file=sys.stderr)
    return 2


if __name__ == "__main__":
    sys.exit(main())
