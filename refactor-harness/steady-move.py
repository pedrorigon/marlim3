#!/usr/bin/env python3
"""Move Stage 7 steady-march bodies out of SisProd.cpp and prove the move literal.

The method is the one thermal-move.py established: rewrite the signature and the
SProd member accesses through an INVERTIBLE substitution table, then prove the
move by applying the inverse to the installed body and comparing tokens with the
pre-move commit. Nothing is retyped, so no arithmetic can drift; and because the
table inverts, the comparison is exact rather than approximate.

That matters more here than usual: coverage measurement found the demo corpus
never executes twelve of the twenty-two routines this stage moves, so for those
the token proof and the dedicated harness are the only verification there is.

Usage:
    steady-move.py extract  <function> <sisprod.cpp> <fragment.cpp>
    steady-move.py install  <function> <sisprod.cpp> <gaslift-in.cpp> <gaslift-out.cpp>
    steady-move.py delegate <function> <sisprod.cpp> <rewritten-sisprod.cpp>
    steady-move.py check    <function> <baseline-sisprod.cpp> <installed.cpp>
"""
from __future__ import annotations

import re
import sys

FUNCTIONS = {
    # The four mass-march variants. T084 measured them: only Rev and CompRev are
    # close enough to share anything (88.7%), and even there 11 of 14 real
    # divergences are logic. So they move as four functions, not as one
    # parameterised core, and the Template Method is applied ONLY to the
    # accessory chain those two share.
    # corrDeng belongs to T087 and lands in this same file. It is moved FIRST
    # because the mass march calls it: the alternative was a callback that would
    # be added and removed within the same stage.
    "corrDeng": {"new_name": "correctGasSpecificGravity", "arguments": "i"},
    # T087 -- pressure march. RenovaPresPermNcel takes no arguments; the tool
    # emits f(state) rather than f(state, ) for that case.
    "RenovaPresPermMon": {"new_name": "advanceUpstreamSteadyPressure", "arguments": "i, RK"},
    "RenovaPresPermJus": {"new_name": "advanceDownstreamSteadyPressure", "arguments": "i, RK"},
    "RenovaPresPermNcel": {"new_name": "steadyPressureAtLastCell", "arguments": ""},
    "calcDpArea": {"new_name": "areaChangePressureDrop", "arguments": "i, rhomix, rey, jmix"},
    # T088 -- phase mass transfer along the march.
    "RenovaTransMassPerm": {"new_name": "advanceSteadyMassTransfer", "arguments": "i"},
    "RenovaTransMassPermGas": {"new_name": "advanceSteadyGasMassTransfer", "arguments": "i"},
    "RenovaMassPerm": {"new_name": "advanceSteadyMass", "arguments": "i"},
    # The fourteen helpers T085 carved out of RenovaMassPerm. They keep their
    # names: they were already named in English when they were created.
    "applySteadyMassWithoutSource": {"new_name": "applySteadyMassWithoutSource", "arguments": "i, mudaRGO, bo, rs, tmed, boI, baI, fwI"},
    "applySteadyMassDryGasInjection": {"new_name": "applySteadyMassDryGasInjection", "arguments": "i, mudaRGO, tL, tH, bo, ba, rs, tmed, trF, boI, baI, fwI"},
    "applySteadyMassWetGasInjection": {"new_name": "applySteadyMassWetGasInjection", "arguments": "i, mudaRGO, tL, tH, bo, ba, rs, tmed, trF, boI, baI, fwI"},
    "applySteadyMassLiquidInjection": {"new_name": "applySteadyMassLiquidInjection", "arguments": "i, mudaRGO, tL, tH, bo, ba, rs, tmed, trF, boI, baI, fwI"},
    "applySteadyMassInflowPerformance": {"new_name": "applySteadyMassInflowPerformance", "arguments": "i, mudaRGO, tL, tH, bo, ba, rs, tmed, trF, boI, baI, fwI"},
    "applySteadyMassMultipleSource": {"new_name": "applySteadyMassMultipleSource", "arguments": "i, mudaRGO, tL, tH, bo, ba, rs, tmed, trF, boI, baI, fwI"},
    "applySteadyMassLeakSource": {"new_name": "applySteadyMassLeakSource", "arguments": "i, mudaRGO, tL, tH, bo, ba, rs, tmed, trF, boI, baI, fwI"},
    "applySteadyMassRadialPorous": {"new_name": "applySteadyMassRadialPorous", "arguments": "i, mudaRGO, tL, tH, bo, ba, rs, tmed, trF, boI, baI, fwI"},
    "applySteadyMassPorous2D": {"new_name": "applySteadyMassPorous2D", "arguments": "i, mudaRGO, tL, tH, bo, ba, rs, tmed, trF, boI, baI, fwI"},
    "finalizeSteadyMassWithLiquid": {"new_name": "finalizeSteadyMassWithLiquid", "arguments": "i, fwI, tmed, bo, rs, pmed, rhog, rhol, qo, masoleo"},
    "finalizeSteadyMassNoFlow": {"new_name": "finalizeSteadyMassNoFlow", "arguments": "i"},
    "finalizeSteadyMassLiquidOnly": {"new_name": "finalizeSteadyMassLiquidOnly", "arguments": "i"},
    "finalizeSteadyMassGasOnly": {"new_name": "finalizeSteadyMassGasOnly", "arguments": "i"},
    "finalizeSteadyMassTwoPhase": {"new_name": "finalizeSteadyMassTwoPhase", "arguments": "i, rhog, rhol"},
    "RenovaMassPermRev": {"new_name": "advanceReverseSteadyMass", "arguments": "i"},
    "RenovaMassPermComp": {"new_name": "advanceCompositionalSteadyMass", "arguments": "i"},
    "RenovaMassPermCompRev": {"new_name": "advanceReverseCompositionalSteadyMass", "arguments": "i"},

    # T097, and the leaves of T094 and T096, moved BEFORE the production marches
    # of T090 rather than after. The marches call all seven, and all seven land
    # in this same module; giving them callbacks first would have meant adding
    # nine entries to SteadyStateUpdaters and deleting them again inside one
    # stage. Same reasoning, and the same decision, as corrDeng above.
    "atualizaProp": {"new_name": "refreshProperties", "arguments": ""},
    "atualizaVelTermPerm": {"new_name": "refreshSteadyThermalVelocities", "arguments": ""},
    "calcDTPseudoTrans": {"new_name": "computePseudoTransientTimeStep", "arguments": ""},
    "atualizaPeriPmonProd": {"new_name": "refreshUpstreamProductionPeriphery", "arguments": "i"},
    "atualizaPeriPjusProd": {"new_name": "refreshDownstreamProductionPeriphery", "arguments": "i"},
    "hidroLinServ": {"new_name": "serviceLineHydrostatic", "arguments": ""},
    "marchaGasPerm1": {"new_name": "marchGasSteady", "arguments": "chutemass"},

    # T090 -- the three production marches and the six pieces they were cut
    # into. seedFirstCellVoidFraction and marchGasLineAndCoupleAnnulus are
    # shared; the other four belong to one march each.
    "seedFirstCellVoidFraction": {"new_name": "seedFirstCellVoidFraction",
                                  "arguments": "pchute, alfini, betini, dryGasFlashTarget"},
    "marchGasLineAndCoupleAnnulus": {"new_name": "marchGasLineAndCoupleAnnulus", "arguments": "pchute"},
    "advanceProductionColumn": {"new_name": "advanceProductionColumn", "arguments": "pchute, i, abortValue"},
    "advanceReverseProductionColumn": {"new_name": "advanceReverseProductionColumn", "arguments": "pchute, i, abortValue"},
    "advanceProductionColumnSecondary": {"new_name": "advanceProductionColumnSecondary", "arguments": "pchute, i, abortValue"},
    "surfaceChokeMassFlow": {"new_name": "surfaceChokeMassFlow", "arguments": ""},
    "marchaProdPerm1": {"new_name": "marchProductionSteady", "arguments": "pchute"},
    "marchaProdPerm1Rev": {"new_name": "marchReverseProductionSteady", "arguments": "pchute"},
    "marchaProdPerm2": {"new_name": "marchProductionSteadySecondary", "arguments": "pchute"},

    # T092, T094, T095 and T096 -- the ten marches that remain. They are moved
    # as one batch, and BEFORE any search, for a reason the task order does not
    # state: every search reaches a march through SProd::zriddr, which calls
    # multMarcha, which dispatches to eight of them. Moving a search first would
    # mean eight callbacks that T097b then deletes. Checked, not assumed: none
    # of these ten calls a search or a solver.
    "marchaProdPresPres1": {"new_name": "marchProductionPressureToPressure", "arguments": "mchute"},
    "marchaProdPresPres1Rev": {"new_name": "marchReverseProductionPressureToPressure", "arguments": "mchute"},
    "marchaProdPresPres2": {"new_name": "marchProductionPressureToPressureSecondary", "arguments": "mchute"},
    "marchaProdPresPres3": {"new_name": "marchProductionPressureToPressureTertiary", "arguments": "mchute"},
    "marchaGasPerm2": {"new_name": "marchGasSteadySecondary", "arguments": "pchute, chutemass"},
    "marchaGasPerm3": {"new_name": "marchGasSteadyTertiary", "arguments": "pchute"},
    "marchaInjPerm1": {"new_name": "marchInjectionSteady", "arguments": "chute"},
    "hidroreverso": {"new_name": "reverseHydrostatic", "arguments": "hol, vaz, vazG"},
    "hidroreversoInj": {"new_name": "reverseInjectionHydrostatic", "arguments": "hol, vaz"},
    "hidroTramoSecundario": {"new_name": "secondaryBranchHydrostatic", "arguments": "titulo"},
}

# Calls BETWEEN moved routines. The moved body must reach the namespace version,
# and a stateful callee needs `state` threaded through; a stateless one must not
# receive it.
# Calls the moved bodies make to OTHER bodies that also live in this module.
# The rewrite inserts the state argument; the inverse removes it again, which is
# what keeps the token proof honest.
#
# This table was inherited from gaslift-move.py still carrying the gas-lift
# entries. corrDeng then survived the rewrite untouched, the module called a
# name that is not in it, and the compiler said so. Swapping FUNCTIONS is not
# enough: CALLS is a second table over the same set.
CALLS = {
    "corrDeng": ("correctGasSpecificGravity", True),
    "RenovaPresPermMon": ("advanceUpstreamSteadyPressure", True),
    "RenovaPresPermJus": ("advanceDownstreamSteadyPressure", True),
    "RenovaPresPermNcel": ("steadyPressureAtLastCell", True),
    "calcDpArea": ("areaChangePressureDrop", True),
    "RenovaTransMassPerm": ("advanceSteadyMassTransfer", True),
    "RenovaTransMassPermGas": ("advanceSteadyGasMassTransfer", True),

    "applySteadyMassWithoutSource": ("applySteadyMassWithoutSource", True),
    "applySteadyMassDryGasInjection": ("applySteadyMassDryGasInjection", True),
    "applySteadyMassWetGasInjection": ("applySteadyMassWetGasInjection", True),
    "applySteadyMassLiquidInjection": ("applySteadyMassLiquidInjection", True),
    "applySteadyMassInflowPerformance": ("applySteadyMassInflowPerformance", True),
    "applySteadyMassMultipleSource": ("applySteadyMassMultipleSource", True),
    "applySteadyMassLeakSource": ("applySteadyMassLeakSource", True),
    "applySteadyMassRadialPorous": ("applySteadyMassRadialPorous", True),
    "applySteadyMassPorous2D": ("applySteadyMassPorous2D", True),
    "finalizeSteadyMassWithLiquid": ("finalizeSteadyMassWithLiquid", True),
    "finalizeSteadyMassNoFlow": ("finalizeSteadyMassNoFlow", True),
    "finalizeSteadyMassLiquidOnly": ("finalizeSteadyMassLiquidOnly", True),
    "finalizeSteadyMassGasOnly": ("finalizeSteadyMassGasOnly", True),
    "finalizeSteadyMassTwoPhase": ("finalizeSteadyMassTwoPhase", True),
    "RenovaMassPerm": ("advanceSteadyMass", True),
    "RenovaMassPermRev": ("advanceReverseSteadyMass", True),
    "RenovaMassPermComp": ("advanceCompositionalSteadyMass", True),
    "RenovaMassPermCompRev": ("advanceReverseCompositionalSteadyMass", True),

    "atualizaProp": ("refreshProperties", True),
    "atualizaVelTermPerm": ("refreshSteadyThermalVelocities", True),
    "calcDTPseudoTrans": ("computePseudoTransientTimeStep", True),
    "atualizaPeriPmonProd": ("refreshUpstreamProductionPeriphery", True),
    "atualizaPeriPjusProd": ("refreshDownstreamProductionPeriphery", True),
    "hidroLinServ": ("serviceLineHydrostatic", True),
    "marchaGasPerm1": ("marchGasSteady", True),

    "seedFirstCellVoidFraction": ("seedFirstCellVoidFraction", True),
    "marchGasLineAndCoupleAnnulus": ("marchGasLineAndCoupleAnnulus", True),
    "advanceProductionColumn": ("advanceProductionColumn", True),
    "advanceReverseProductionColumn": ("advanceReverseProductionColumn", True),
    "advanceProductionColumnSecondary": ("advanceProductionColumnSecondary", True),
    "surfaceChokeMassFlow": ("surfaceChokeMassFlow", True),

    "marchaProdPresPres1": ("marchProductionPressureToPressure", True),
    "marchaProdPresPres1Rev": ("marchReverseProductionPressureToPressure", True),
    "marchaProdPresPres2": ("marchProductionPressureToPressureSecondary", True),
    "marchaProdPresPres3": ("marchProductionPressureToPressureTertiary", True),
    "marchaGasPerm2": ("marchGasSteadySecondary", True),
    "marchaGasPerm3": ("marchGasSteadyTertiary", True),
    "marchaInjPerm1": ("marchInjectionSteady", True),
    "hidroreverso": ("reverseHydrostatic", True),
    "hidroreversoInj": ("reverseInjectionHydrostatic", True),
    "hidroTramoSecundario": ("secondaryBranchHydrostatic", True),
}

# SProd member -> SteadyStateState field. Longest first when the pattern is built, so
# a short name cannot claim the prefix of a longer one (ncel vs ncelGas).
# Calls the moved bodies make back into SProd, routed through an adapter the
# state carries. tempDescarga is still an SProd member; the gas line reaches the
# thermal module through it rather than growing a ThermalState of its own.
CALLBACKS = {
    "CalcC0UdPerm": "state.updaters.steadyDriftClosure",
    "renovaFonte": "state.updaters.updateSource",
    "RenovaTempPermRev": "state.updaters.advanceReverseSteadyTemperature",
    "RenovaTempPerm": "state.updaters.advanceSteadyTemperature",
    "atualizaPeriTempProd": "state.updaters.updateProductionTemperaturePeriphery",
    "calctempGas": "state.updaters.computeGasTemperature",
    "calctemp": "state.updaters.computeTemperature",
    "IniciaVazValvGasPerm": "state.updaters.initializeSteadyValveGasFlowRate",
    "IniciaconectaColunaPerm": "state.updaters.initializeTubingConnectionSteady",
    "RenovaPresGasPerm": "state.updaters.updateSteadyGasPressure",
    "RenovaTempGasPerm": "state.updaters.updateSteadyGasTemperature",
    "calcVazGasPerm": "state.updaters.computeSteadyGasFlowRate",
    "conectaColunaPerm": "state.updaters.connectTubingSteady",
    "conectaColuna": "state.updaters.connectTubing",
    "delpGasPerm": "state.updaters.steadyGasPressureDrop",
    "delpInjPerm": "state.updaters.steadyInjectionPressureDrop",
    # These two ARE searches, and a march calls them. See section 8 of
    # evidencia/marchaprod-diff.md: the call graph between the two halves of
    # stage 7 has a cycle, and it is the SProd indirection -- not an absent
    # edge -- that keeps the INCLUDE one-way.
    "buscaGasPresPerm2": "state.updaters.searchGasPressureSteadySecondary",
    "buscaGasPresPerm3": "state.updaters.searchGasPressureSteadyTertiary",
}

MEMBERS = {
    # Longest first when the pattern is built, so a short name cannot claim the
    # prefix of a longer one: ncelGas before ncel, monitConvPermBase before
    # monitConvPerm, tempiniG before tempRev.
    "monitConvPermBase": "state.baseConvergenceMonitor",
    "monitConvPerm": "state.convergenceMonitor",
    "trocaTermicaLenta": "state.slowHeatTransfer",
    "verificaAcop": "state.networkCoupled",
    "celulaG": "state.gasCells",
    "celula": "state.cells",
    "ncelGas": "state.gasCellCount",
    "ncel": "state.lastCell",
    "posicVGLG": "state.gasValveCellIndices",
    "posicVGLP": "state.productionValveCellIndices",
    "chokeInj": "state.injectionChoke",
    "chokeSup": "state.surfaceChoke",
    "iterperm": "state.steadyIteration",
    "buscaIni": "state.searchOrigin",
    "derivaAnel": "state.annulusDrift",
    "noextremo": "state.endNode",
    "semTermo": "state.thermalSourceDisabled",
    "presiniG": "state.initialGasPressure",
    "tempiniG": "state.initialGasTemperature",
    "presfim": "state.finalPressure",
    "pGSup": "state.gasSurfacePressure",
    "temperatura": "state.ambientTemperature",
    "tempRev": "state.casingTemperature",
    "nfluP": "state.productionFluidCount",
    "titE": "state.inletQuality",
    "vg1dSP": "state.globals",
    "arq": "state.input",
    "dt": "state.timeStep",
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
        head = (f"{return_type}{new_name}(const SteadyStateState &state"
                f"{separator}{signature}) {{")
    body = SIGNATURE_RE.sub(lambda _: head, body, count=1)
    body = substitute_outside_comments(MEMBER_RE,
                                       lambda m: MEMBERS[m.group(1)], body)
    for called, (renamed, needs_state) in CALLS.items():
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
        if needs_state:
            def add_state(match, renamed=renamed):
                return f"{renamed}(state)" if match.group(1) else f"{renamed}(state, "
            body = substitute_outside_comments(
                re.compile(rf"(?<![\w.>]){called}\((\s*\))?"), add_state, body)
        else:
            body = substitute_outside_comments(
                re.compile(rf"(?<![\w.>]){called}\("), f"{renamed}(", body)
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
        if needs_state:
            body = re.sub(rf"(?<![\w.>]){re.escape(renamed)}\(state\)",
                          f"{called}()", body)
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
            rf"\(const SteadyStateState &state(?:, (.*?))?\)\s*\{{", re.S)
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
                    f"(steadyStateOf(*this){separator}{arguments})")
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
                   else r"\(const SteadyStateState &state")
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
