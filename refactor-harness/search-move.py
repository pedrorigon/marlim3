#!/usr/bin/env python3
"""Move Stage 7 SEARCH bodies out of SisProd.cpp and prove the move literal.

Same method as steady-move.py -- rewrite the signature and the SProd member
accesses through an INVERTIBLE substitution table, then prove the move by
applying the inverse to the installed body and comparing tokens with the
pre-move commit.

It is a separate tool rather than a flag on the other one because the two
sides read DIFFERENT state. A march takes SteadyStateState and spells a cell
`state.cells`. A search takes SteadyStateSearchState, which COMPOSES the march
state, so the same cell is `state.march.cells`, and every march it calls has to
be handed `state.march` rather than `state`. One table cannot say both, and a
flag that switched between them would be a table with a mode -- which is how
steady-move.py inherited gaslift-move.py's CALLS and called a name that did not
exist.

Usage:
    search-move.py extract  <function> <sisprod.cpp> <fragment.cpp>
    search-move.py install  <function> <sisprod.cpp> <search-in.cpp> <search-out.cpp>
    search-move.py delegate <function> <sisprod.cpp> <rewritten-sisprod.cpp>
    search-move.py check    <function> <baseline-sisprod.cpp> <installed.cpp>
"""
from __future__ import annotations

import re
import sys

FUNCTIONS = {
    # T091 -- the four bottom-hole searches. Three of them are searches; the
    # fourth, buscaProdPfundoPerm3, calls no solver and no march and inlines a
    # march of its own. Measured at 4-5% similarity to the other three, where
    # those three sit at 63-69% to each other. See evidencia/buscaprod-diff.md.
    "buscaProdPfundoPerm": {"new_name": "searchProductionBottomHolePressure", "arguments": "chute, kontaTenta"},
    "buscaProdPfundoPermRev": {"new_name": "searchReverseProductionBottomHolePressure", "arguments": "chute"},
    "buscaProdPfundoPerm2": {"new_name": "searchProductionBottomHolePressureSecondary", "arguments": "chute, kontaTenta"},
    "buscaProdPfundoPerm3": {"new_name": "searchProductionBottomHolePressureTertiary", "arguments": "pentrada"},

    # T093 -- pressure-to-pressure searches.
    "buscaProdPresPresPerm": {"new_name": "searchProductionPressureToPressure", "arguments": "mchute, vazmax, kontaiter"},
    "buscaProdPresPresPermRev": {"new_name": "searchReverseProductionPressureToPressure", "arguments": "mchute, vazmax, kontaiter"},
    "buscaProdPresPresPerm2": {"new_name": "searchProductionPressureToPressureSecondary", "arguments": "mchute, vazmax"},
    "buscaProdPresPresPerm3": {"new_name": "searchProductionPressureToPressureTertiary", "arguments": "mchute, vazmax"},

    # T094 -- gas-line searches. A march calls these two; see section 8 of
    # evidencia/marchaprod-diff.md for what that costs the stage-7 story.
    "buscaGasPresPerm2": {"new_name": "searchGasPressureSteadySecondary", "arguments": ""},
    "buscaGasPresPerm3": {"new_name": "searchGasPressureSteadyTertiary", "arguments": ""},

    # T095 -- five injection searches, numbered because the original numbered
    # them and nothing says what separates four from five.
    "buscaInjPfundoPerm1": {"new_name": "searchInjectionBottomHolePressure1", "arguments": "chute"},
    "buscaInjPfundoPerm2": {"new_name": "searchInjectionBottomHolePressure2", "arguments": "chute"},
    "buscaInjPfundoPerm3": {"new_name": "searchInjectionBottomHolePressure3", "arguments": "chute"},
    "buscaInjPfundoPerm4": {"new_name": "searchInjectionBottomHolePressure4", "arguments": ""},
    "buscaInjPfundoPerm5": {"new_name": "searchInjectionBottomHolePressure5", "arguments": "chute"},

    # T096 -- the secondary-branch search.
    "buscaTramoSecVazPerm": {"new_name": "searchSecondaryBranchFlowRate", "arguments": "chute"},
}

# Calls a search makes. A march lives in the OTHER module and takes the march
# state, so it is handed state.march. A search in this module takes the search
# state and is handed state.
CALLS = {
    "marchaProdPerm1": ("marchProductionSteady", "state.march"),
    "marchaProdPerm1Rev": ("marchReverseProductionSteady", "state.march"),
    "marchaProdPerm2": ("marchProductionSteadySecondary", "state.march"),
    "marchaProdPresPres1": ("marchProductionPressureToPressure", "state.march"),
    "marchaProdPresPres1Rev": ("marchReverseProductionPressureToPressure", "state.march"),
    "marchaProdPresPres2": ("marchProductionPressureToPressureSecondary", "state.march"),
    "marchaProdPresPres3": ("marchProductionPressureToPressureTertiary", "state.march"),
    "marchaGasPerm1": ("marchGasSteady", "state.march"),
    "marchaGasPerm2": ("marchGasSteadySecondary", "state.march"),
    "marchaGasPerm3": ("marchGasSteadyTertiary", "state.march"),
    "marchaInjPerm1": ("marchInjectionSteady", "state.march"),
    "hidroreverso": ("reverseHydrostatic", "state.march"),
    "hidroreversoInj": ("reverseInjectionHydrostatic", "state.march"),
    "hidroTramoSecundario": ("secondaryBranchHydrostatic", "state.march"),
    "hidroLinServ": ("serviceLineHydrostatic", "state.march"),
    "atualizaProp": ("refreshProperties", "state.march"),
    "atualizaVelTermPerm": ("refreshSteadyThermalVelocities", "state.march"),
    "calcDTPseudoTrans": ("computePseudoTransientTimeStep", "state.march"),
    "atualizaPeriPmonProd": ("refreshUpstreamProductionPeriphery", "state.march"),
    "atualizaPeriPjusProd": ("refreshDownstreamProductionPeriphery", "state.march"),
    "RenovaMassPerm": ("advanceSteadyMass", "state.march"),
    "RenovaMassPermRev": ("advanceReverseSteadyMass", "state.march"),
    "RenovaMassPermComp": ("advanceCompositionalSteadyMass", "state.march"),
    "RenovaMassPermCompRev": ("advanceReverseCompositionalSteadyMass", "state.march"),
    "RenovaPresPermMon": ("advanceUpstreamSteadyPressure", "state.march"),
    "RenovaPresPermJus": ("advanceDownstreamSteadyPressure", "state.march"),
    "RenovaPresPermNcel": ("steadyPressureAtLastCell", "state.march"),
    "calcDpArea": ("areaChangePressureDrop", "state.march"),
    "RenovaTransMassPerm": ("advanceSteadyMassTransfer", "state.march"),
    "RenovaTransMassPermGas": ("advanceSteadyGasMassTransfer", "state.march"),
    "corrDeng": ("correctGasSpecificGravity", "state.march"),

    "zriddr": ("solveSteadyRoot", "state"),
    "multMarcha": ("dispatchMarch", "state"),
    "buscaProdPfundoPerm": ("searchProductionBottomHolePressure", "state"),
    "buscaProdPfundoPermRev": ("searchReverseProductionBottomHolePressure", "state"),
    "buscaProdPfundoPerm2": ("searchProductionBottomHolePressureSecondary", "state"),
    "buscaProdPfundoPerm3": ("searchProductionBottomHolePressureTertiary", "state"),
    "buscaGasPresPerm2": ("searchGasPressureSteadySecondary", "state"),
    "buscaGasPresPerm3": ("searchGasPressureSteadyTertiary", "state"),
    "buscaInjPfundoPerm1": ("searchInjectionBottomHolePressure1", "state"),
    "buscaInjPfundoPerm2": ("searchInjectionBottomHolePressure2", "state"),
    "buscaInjPfundoPerm3": ("searchInjectionBottomHolePressure3", "state"),
    "buscaInjPfundoPerm4": ("searchInjectionBottomHolePressure4", "state"),
    "buscaInjPfundoPerm5": ("searchInjectionBottomHolePressure5", "state"),
}

# SProd member -> SteadyStateState field. Longest first when the pattern is built, so
# a short name cannot claim the prefix of a longer one (ncel vs ncelGas).
# Calls the moved bodies make back into SProd, routed through an adapter the
# state carries. tempDescarga is still an SProd member; the gas line reaches the
# thermal module through it rather than growing a ThermalState of its own.
CALLBACKS = {
    "CalcC0UdPerm": "state.march.updaters.steadyDriftClosure",
    "renovaFonte": "state.march.updaters.updateSource",
    "RenovaTempPermRev": "state.march.updaters.advanceReverseSteadyTemperature",
    "RenovaTempPerm": "state.march.updaters.advanceSteadyTemperature",
    "atualizaPeriTempProd": "state.march.updaters.updateProductionTemperaturePeriphery",
    "calctempGas": "state.march.updaters.computeGasTemperature",
    "calctemp": "state.march.updaters.computeTemperature",
    "IniciaVazValvGasPerm": "state.march.updaters.initializeSteadyValveGasFlowRate",
    "IniciaconectaColunaPerm": "state.march.updaters.initializeTubingConnectionSteady",
    "RenovaPresGasPerm": "state.march.updaters.updateSteadyGasPressure",
    "RenovaTempGasPerm": "state.march.updaters.updateSteadyGasTemperature",
    "calcVazGasPerm": "state.march.updaters.computeSteadyGasFlowRate",
    "conectaColunaPerm": "state.march.updaters.connectTubingSteady",
    "conectaColuna": "state.march.updaters.connectTubing",
    "delpGasPerm": "state.march.updaters.steadyGasPressureDrop",
    "delpInjPerm": "state.march.updaters.steadyInjectionPressureDrop",
    # These two ARE searches, and a march calls them. See section 8 of
    # evidencia/marchaprod-diff.md: the call graph between the two halves of
    # stage 7 has a cycle, and it is the SProd indirection -- not an absent
    # edge -- that keeps the INCLUDE one-way.
    "buscaGasPresPerm2": "state.march.updaters.searchGasPressureSteadySecondary",
    "buscaGasPresPerm3": "state.march.updaters.searchGasPressureSteadyTertiary",
}

MEMBERS = {
    # The three the search owns. They are NOT under .march: SteadyStateSearchState
    # adds them, and putting them there would be a second spelling of a fact the
    # composed struct already states once.
    "revPerm": "state.reverseSteady",
    "chuteHol": "state.holdupGuess",
    "fluiRevRede": "state.reverseNetworkFluid",

    # Longest first when the pattern is built, so a short name cannot claim the
    # prefix of a longer one: ncelGas before ncel, monitConvPermBase before
    # monitConvPerm, tempiniG before tempRev.
    "monitConvPermBase": "state.march.baseConvergenceMonitor",
    "monitConvPerm": "state.march.convergenceMonitor",
    "trocaTermicaLenta": "state.march.slowHeatTransfer",
    "verificaAcop": "state.march.networkCoupled",
    "celulaG": "state.march.gasCells",
    "celula": "state.march.cells",
    "ncelGas": "state.march.gasCellCount",
    "ncel": "state.march.lastCell",
    "posicVGLG": "state.march.gasValveCellIndices",
    "posicVGLP": "state.march.productionValveCellIndices",
    "chokeInj": "state.march.injectionChoke",
    "chokeSup": "state.march.surfaceChoke",
    "iterperm": "state.march.steadyIteration",
    "buscaIni": "state.march.searchOrigin",
    "derivaAnel": "state.march.annulusDrift",
    "noextremo": "state.march.endNode",
    "semTermo": "state.march.thermalSourceDisabled",
    "presiniG": "state.march.initialGasPressure",
    "tempiniG": "state.march.initialGasTemperature",
    "presfim": "state.march.finalPressure",
    "pGSup": "state.march.gasSurfacePressure",
    "temperatura": "state.march.ambientTemperature",
    "tempRev": "state.march.casingTemperature",
    "nfluP": "state.march.productionFluidCount",
    "titE": "state.march.inletQuality",
    "vg1dSP": "state.march.globals",
    "arq": "state.march.input",
    "dt": "state.march.timeStep",
}
# steady-move.py splits the field on its FIRST dot, which is right when every
# field is state.<name>. Here most are state.march.<name>, and that split would
# map them all to the key "march" -- one key, twenty-eight values, an inverse
# that silently loses twenty-seven. Split on the last dot, and key the regex on
# the full path so state.march.cells and a hypothetical state.cells cannot be
# confused.
FIELD_TO_MEMBER = {field: member for member, field in MEMBERS.items()}
MEMBER_RE = re.compile(
    r"(?<![\w.])(" + "|".join(sorted(MEMBERS, key=len, reverse=True)) + r")\b"
)
FIELD_RE = re.compile(
    "(?<![\\w.])(" + "|".join(re.escape(f) for f in sorted(FIELD_TO_MEMBER, key=len, reverse=True)) + r")\b"
)
SIGNATURE_RE = re.compile(
    r"^([\w:<>*&]+(?:\s+[\w:<>*&]+)*\s+)SProd::(\w+)\((.*?)\)\s*\{", re.S
)

COMMENT = re.compile(r"//[^\n]*|/\*.*?\*/", re.S)


def compare_tokens(expected: list[str], actual: list[str]) -> str:
    """Delegates to thermal-move.py, which owns the calibrated implementation."""
    import importlib.util
    import pathlib
    spec = importlib.util.spec_from_file_location(
        "thermal_move", pathlib.Path(__file__).with_name("thermal-move.py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.compare_tokens(expected, actual)


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
        head = (f"{return_type}{new_name}(const SteadyStateSearchState &state"
                f"{separator}{signature}) {{")
    body = SIGNATURE_RE.sub(lambda _: head, body, count=1)
    body = substitute_outside_comments(MEMBER_RE,
                                       lambda m: MEMBERS[m.group(1)], body)
    for called, (renamed, state_expr) in CALLS.items():
        if called == old_name:
            continue
        # A zero-argument call must become f(state), not f(state, ).
        #
        # This used to be two passes, the empty form first. That worked only
        # while every moved function was RENAMED: with called != renamed, the
        # first pass produced text the second could not match again. T090 moves
        # six helpers that keep their names, and on those the first pass emitted
        # f(state) and the second immediately rewrote it to f(state, state).
        # The compiler would have caught it, but the token proof caught it
        # first, which is the point of having one.
        #
        # One pass now decides between the two forms, so no output of this
        # substitution is input to it.
        def add_state(match, renamed=renamed, expr=state_expr):
            return f"{renamed}({expr})" if match.group(1) else f"{renamed}({expr}, "
        body = substitute_outside_comments(
            re.compile(rf"(?<![\w.>]){called}\((\s*\))?"), add_state, body)
    for member, adapter in CALLBACKS.items():
        body = substitute_outside_comments(
            re.compile(rf"(?<![\w.>]){member}\("), f"{adapter}(", body)
    return body


def inverse(old_name: str, body: str) -> str:
    spec = FUNCTIONS[old_name]
    new_name = spec["new_name"]
    stateless = spec.get("stateless", False)
    for called, (renamed, state_expr) in reversed(list(CALLS.items())):
        if called == old_name:
            continue
        body = re.sub(rf"(?<![\w.>]){re.escape(renamed)}\({re.escape(state_expr)}\)",
                      f"{called}()", body)
        body = re.sub(rf"(?<![\w.>]){re.escape(renamed)}\({re.escape(state_expr)}, ",
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
            rf"\(const SteadyStateSearchState &state(?:, (.*?))?\)\s*\{{", re.S)
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
        closing = "}  // namespace sisprod::steady\n"
        if target.count(closing) != 1:
            raise ValueError("expected one steady namespace closing marker")
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
            call = f"sisprod::steady::{spec['new_name']}({arguments})"
        else:
            separator = ", " if arguments else ""
            call = (f"sisprod::steady::{spec['new_name']}"
                    f"(searchStateOf(*this){separator}{arguments})")
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
                   else r"\(const SteadyStateSearchState &state")
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
        # T081r renames locals, so an exact match stops being possible the
        # moment it runs. Alpha-equivalence is the same test thermal-move.py
        # already uses and calibrates: it accepts a CONSISTENT renaming and
        # still rejects a changed literal, a reordering, or one variable
        # substituted for another. Reused rather than reimplemented, so there
        # is one implementation to trust and one to calibrate.
        if compare_tokens(expected, actual) == "alpha":
            print(f"OK       {new_name} <- {name} "
                  f"({len(expected)} tokens modulo renames)")
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
