#!/usr/bin/env python3
"""Move SolveTrans and the six helpers it was cut into, and prove the move.

A fourth move tool, for the reason the second and third record: the tables ARE
the tool. SolveTrans reads TransientSolveState, which COMPOSES the step state --
so a step field is spelled state.step.X here, the 52 solve-only fields are
state.X, and every step function it calls is handed state.step. One table cannot
say both, which is exactly why search-move.py is not a mode of steady-move.py.

The step half of every table below is DERIVED from transient-move.py rather than
copied from it, so the two cannot drift.

Usage: as transient-move.py.
"""
from __future__ import annotations

import re
import sys

FUNCTIONS = {
    "SolveTrans": {"new_name": "solveTransientStep", "arguments": "titRev, alfRev, betRev, nrede, fluiRev"},
    "advanceCouplingIteration": {"new_name": "advanceCouplingIteration",
                                 "arguments": "kontaAcop, celpos, vExpli, ciclomax, titRev, alfRev, betRev"},
    "writeProgressReport": {"new_name": "writeProgressReport", "arguments": "MaxKontaImpres"},
    "writeEventLog": {"new_name": "writeEventLog", "arguments": "maxEvento"},
    "writeScreenOutput": {"new_name": "writeScreenOutput", "arguments": "begin, end"},
    "writeTrends": {"new_name": "writeTrends", "arguments": "ordemImpT, velmaxdesc, nrede"},
    "writeProfiles": {"new_name": "writeProfiles", "arguments": "nrede"},
}

CALLS = {
    "renova": ("updateCells", "state.step"),
    "renovaCelulaInterior": ("updateInteriorCell", "state.step"),
    "renovaPrimeiraCelula": ("updateFirstCell", "state.step"),
    "renovaUltimaCelula": ("updateLastCell", "state.step"),
    "renovaVaz": ("updateFlowRates", "state.step"),
    "renovaBuffer": ("updateBufferFromSolution", "state.step"),
    "renovaBufferCego": ("updateBufferFromCells", "state.step"),
    "surfaceChokeIsOpen": ("surfaceChokeIsOpen", "state.step"),
    "surfaceChokeIsShut": ("surfaceChokeIsShut", "state.step"),
    "calcCCpres": ("applyOutletPressureCondition", "state.step"),
    "calcCCBuffer": ("applyOutletBufferCondition", "state.step"),
    "computeImplicitTimeStep": ("computeImplicitTimeStep", "state.step"),
    "determinaDTExpli": ("computeExplicitTimeStep", "state.step"),
    "determinaDT": ("computeTimeStep", "state.step"),
    "atenuaDtMax": ("dampMaximumTimeStep", "state.step"),
    "avaliaVariaDpDt": ("evaluatePressureRateOfChange", "state.step"),
    "restringeDTporValv": ("restrictTimeStepByValve", "state.step"),
    "aberturaVal0": ("valveOpeningLow", "state.step"),
    "aberturaVal1": ("valveOpeningHigh", "state.step"),
    "EvoluiFrac": ("evolveFractions", "state.step"),
    "ReiniEvolFrac0": ("restartFractionEvolutionInitial", "state.step"),
    "SubReiniEvolFrac": ("restartFractionEvolutionSub", "state.step"),
    "ReiniEvolFrac": ("restartFractionEvolution", "state.step"),
    "AtualizaPig": ("updatePig", "state.step"),
    "SolveAcopPV": ("solvePressureVolumeCoupling", "state.step"),
    "atualizaMiniTab": ("refreshFluidMiniTable", "state.step"),
    "atualizaCC1": ("refreshInletCondition", "state.step"),
    "advanceCouplingIteration": ("advanceCouplingIteration", "state"),
    "writeProgressReport": ("writeProgressReport", "state"),
    "writeEventLog": ("writeEventLog", "state"),
    "writeScreenOutput": ("writeScreenOutput", "state"),
    "writeTrends": ("writeTrends", "state"),
    "writeProfiles": ("writeProfiles", "state"),
}

# SProd member -> SteadyStateState field. Longest first when the pattern is built, so
# a short name cannot claim the prefix of a longer one (ncel vs ncelGas).
# Calls the moved bodies make back into SProd, routed through an adapter the
# state carries. tempDescarga is still an SProd member; the gas line reaches the
# thermal module through it rather than growing a ThermalState of its own.
# Nothing yet: the renova group reaches no SProd method. SolveTrans will add
# entries here when T127 moves it.
CALLBACKS = {
    "geraMiniTabFlu": "state.step.updaters.generateFluidMiniTable",
    "subtempoGas": "state.step.updaters.advanceGasSubStep",
    "solveHydrateEnvelopes": "state.updaters.solveHydrateEnvelopes",
    "BuscaPresInjDesc": "state.updaters.findInjectionPressureDownstream",
    "ImprimeTrendPCab": "state.updaters.writeProductionTrendHeader",
    "ImprimeTrendP": "state.updaters.writeProductionTrendRows",
    "ImprimeTrendGCab": "state.updaters.writeGasTrendHeader",
    "ImprimeTrendG": "state.updaters.writeGasTrendRows",
    "ImprimeTrendTransPCab": "state.updaters.writeProductionCrossSectionTrendHeader",
    "ImprimeTrendTransP": "state.updaters.writeProductionCrossSectionTrendRows",
    "ImprimeTrendTransGCab": "state.updaters.writeGasCrossSectionTrendHeader",
    "ImprimeTrendTransG": "state.updaters.writeGasCrossSectionTrendRows",
    "avaliaParafina": "state.updaters.evaluateParaffin",
    "conectaColuna": "state.updaters.connectTubing",
    "marchaEnergTrans": "state.updaters.marchTransientEnergy",
    "renovaFracMol2": "state.updaters.updateMolarFractions",
    "renovaMasEsp": "state.updaters.updateDensities",
    "renovaRGOdgYco2": "state.updaters.updateGasOilRatioAndCo2",
    "renovaTemp": "state.updaters.updateTemperatures",
    "renovaalbetini": "state.updaters.updateInitialFractions",
    "renovaterm": "state.updaters.updateThermal",
    "salvaFonte": "state.updaters.saveSources",
    "solveLinGas": "state.updaters.solveGasLine",
}

MEMBERS = {
    "kontaTempoTransProfG": "state.gasCrossSectionProfileTimeCounter",
    "kontaTempoTransProf": "state.productionCrossSectionProfileTimeCounter",
    "kontaTempoCelUni": "state.unitCellTimeCounters",
    "resettrendtransg": "state.gasCrossSectionTrendResetTimers",
    "saidaSubTextoSis": "state.closingSubtitles",
    "kontaRenovaComp": "state.compositionalRefreshCounter",
    "kontaTempoProfG": "state.gasProfileTimeCounter",
    "resettrendtrans": "state.productionCrossSectionTrendResetTimers",
    "fontemassCRBuf": "state.step.bufferedCompletionMassSource",
    "fontemassGRBuf": "state.step.bufferedGasMassSource",
    "fontemassPRBuf": "state.step.bufferedLiquidMassSource",
    "modeloCompleto": "state.step.fullModel",
    "ncelperftransp": "state.step.productionCrossSectionCount",
    "MatTrendTransG": "state.gasCrossSectionTrendMatrix",
    "MatTrendTransP": "state.productionCrossSectionTrendMatrix",
    "TransMassModel": "state.massTransferModel",
    "kontaTempoProf": "state.productionProfileTimeCounter",
    "ncelperftransg": "state.gasCrossSectionCellCounts",
    "tempoabertoini": "state.initialOpenTime",
    "EstadoMaster1": "state.step.masterState",
    "kontarestriDt": "state.step.timeStepRestrictionCount",
    "momentoDesesp": "state.step.desperationMoment",
    "ntrendtransgB": "state.gasCrossSectionTrendBufferedCounts",
    "saidaTextoSis": "state.closingTitles",
    "contaMaster1": "state.step.masterCounter",
    "kontaGolfada": "state.step.slugCount",
    "vRazMastCrit": "state.step.masterCriticalRatio",
    "KontaImprime": "state.printCounter",
    "ntrendtransB": "state.productionCrossSectionTrendBufferedCounts",
    "ntrendtransg": "state.gasCrossSectionTrendCounts",
    "verificaAcop": "state.networkCoupled",
    "alteraTempo": "state.step.timeChanged",
    "mudaModoChk": "state.step.chokeModeChanged",
    "tempoaberto": "state.step.openTime",
    "termolivreP": "state.step.productionSolution",
    "celInterIni": "state.step.initialInterfaceCell",
    "velInterIni": "state.step.initialInterfaceVelocity",
    "ntrendtrans": "state.productionCrossSectionTrendCounts",
    "resettrendg": "state.gasTrendResetTimers",
    "temperatura": "state.ambientTemperature",
    "dtCFLTotal": "state.step.totalCflTimeStep",
    "dtSimTotal": "state.step.totalSimulationTimeStep",
    "dtInterIni": "state.step.initialInterfaceTimeStep",
    "dtauxFinal": "state.step.finalAuxiliaryTimeStep",
    "presMedMov": "state.step.movingMeanPressure",
    "derivaAnel": "state.annulusDrift",
    "resettrend": "state.productionTrendResetTimers",
    "abertoini": "state.step.initiallyOpen",
    "indevento": "state.step.eventIndex",
    "masChkSup": "state.step.surfaceChokeMassFlag",
    "vRazMast0": "state.step.masterRatio0",
    "vRazMast1": "state.step.masterRatio1",
    "noextremo": "state.step.endNode",
    "taxaDTMax": "state.step.maximumTimeStepRates",
    "taxaDpMax": "state.step.maximumPressureRates",
    "MatTrendG": "state.gasTrendMatrix",
    "MatTrendP": "state.productionTrendMatrix",
    "alfMedMov": "state.movingMeanVoidFraction",
    "noinicial": "state.startNode",
    "poisson3D": "state.poissonSolver3D",
    "trackDeng": "state.trackGasGravity",
    "DTMaxMed": "state.step.meanMaximumTimeStep",
    "DpMaxMed": "state.step.meanMaximumPressureChange",
    "celInter": "state.step.interfaceCell",
    "dtCFLMed": "state.step.meanCflTimeStep",
    "dtSimMed": "state.step.meanSimulationTimeStep",
    "reinicia": "state.step.restart",
    "restriDt": "state.step.timeStepRestricted",
    "velInter": "state.step.interfaceVelocity",
    "chokeSup": "state.step.surfaceChoke",
    "dtauxCFL": "state.step.auxiliaryCflTimeStep",
    "indTramo": "state.step.branchIndex",
    "matglobP": "state.step.productionMatrix",
    "nfechaM1": "state.step.masterCloseCount",
    "alfTotal": "state.totalVoidFraction",
    "chokeInj": "state.injectionChoke",
    "contaLog": "state.logCounter",
    "dtCicMin": "state.minimumCycleTimeStep",
    "ktMedMov": "state.movingMeanCounter",
    "ntrendgB": "state.gasTrendBufferedCounts",
    "pGSupIni": "state.initialGasSurfacePressure",
    "presiniG": "state.initialGasPressure",
    "tempiniG": "state.initialGasTemperature",
    "trackRGO": "state.trackGasOilRatio",
    "dtInter": "state.step.interfaceTimeStep",
    "presfim": "state.step.finalPressure",
    "celulaG": "state.step.gasCells",
    "fechaM1": "state.step.masterCloseSchedule",
    "jMedMov": "state.step.movingMeanFlux",
    "menorDx": "state.step.smallestCellLength",
    "nabreM1": "state.step.masterOpenCount",
    "ncelGas": "state.step.gasCellCount",
    "tMedMov": "state.step.movingMeanTemperature",
    "ntrendB": "state.productionTrendBufferedCounts",
    "ntrendg": "state.gasTrendCounts",
    "presVet": "state.pressureHistory",
    "aberto": "state.step.open",
    "abreM1": "state.step.masterOpenSchedule",
    "celula": "state.step.cells",
    "titRev": "state.step.reverseQuality",
    "vg1dSP": "state.step.globals",
    "alfVet": "state.voidFractionHistory",
    "jTotal": "state.totalFlux",
    "ntrend": "state.productionTrendCounts",
    "pTotal": "state.totalPressure",
    "tmpLog": "state.logBuffer",
    "betaE": "state.step.inletCompletionFraction",
    "dtCFL": "state.step.cflTimeSteps",
    "dtSim": "state.step.simulationTimeSteps",
    "flutG": "state.step.gasFreeTerms",
    "pGSup": "state.step.gasSurfacePressure",
    "presE": "state.step.inletPressure",
    "tempE": "state.step.inletTemperature",
    "kimpT": "state.printTimeCounter",
    "mult": "state.step.multiplier",
    "titE": "state.step.inletQuality",
    "flut": "state.step.productionFreeTerms",
    "ncel": "state.step.lastCell",
    "jVet": "state.fluxHistory",
    "tVet": "state.temperatureHistory",
    "cpg": "state.step.gasSpecificHeatTable",
    "arq": "state.step.input",
    "kSP": "state.step.stepIndex",
    "dt": "state.step.timeStep",
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



# `a * b` is either a multiplication or a pointer declaration and no lexical
# rule tells them apart -- `tit = alfSup * celula[ncel]...` parsed as a
# declaration of `celula`, which excluded a real member from substitution and
# broke the build in the opposite direction from the bug the guard exists for.
#
# So a declaration must START a statement: preceded by line start, `;`, `{` or
# `}`. Parameters are read from the signature separately, where commas are
# unambiguous.
# Anchored on `;`, `{` or `}` and NOT on line start: a wrapped expression's
# continuation line begins a line without beginning a statement, and
# `betSup * celula[ncel]...` on such a line read as a declaration of celula.
# The "type" is captured so keywords can be rejected. Without that,
# `; else momentoDesesp = 1;` parses as a declaration of momentoDesesp, which
# then excluded a real member from substitution -- the same class of
# over-exclusion as the `<` and `*` cases above, in a third disguise.
DECLARED = re.compile(
    r"[;{}]\s*(?:const\s+)?([A-Za-z_][\w:]*)(?:<[^<>;]*>)?(?:\s*\*)?\s+"
    r"[&*]?([A-Za-z_]\w*)\s*(?:=[^=]|;|,|\[)")

NOT_A_TYPE = {"return", "else", "do", "case", "break", "continue", "goto",
              "new", "delete", "throw", "sizeof", "static", "typedef"}

PARAM = re.compile(r"[&*]?([A-Za-z_]\w*)\s*(?:=[^=][^,)]*)?$")


def _declared_in_body(body: str) -> set:
    """Names the body declares itself: parameters and locals.

    Anything here shadows the SProd member of the same name, so the member
    substitution must leave it alone. calcCCpres takes a parameter `titRev` and
    declares locals `cpg` and `abertoini`, and SProd has members of all three.
    """
    clean = COMMENT.sub(" ", body)
    names = set()
    open_paren = clean.find("(")
    close_paren = clean.find(")", open_paren)
    if open_paren >= 0 and close_paren > open_paren:
        for part in clean[open_paren + 1:close_paren].split(","):
            match = PARAM.search(part.strip())
            if match:
                names.add(match.group(1))
    for match in DECLARED.finditer(clean):
        if match.group(1) in NOT_A_TYPE:
            continue
        names.add(match.group(2))
    return names


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
        head = (f"{return_type}{new_name}(const TransientSolveState &state"
                f"{separator}{signature}) {{")
    body = SIGNATURE_RE.sub(lambda _: head, body, count=1)
    # A local or parameter whose name equals an SProd member must NOT be
    # rewritten as that member. calcCCpres takes a parameter `titRev` and
    # declares a local `cpg`, and SProd has members of both names; without this
    # guard the tool produced `double state.reverseQuality` in the signature.
    #
    # That particular breakage is loud -- a declaration with a dotted name is a
    # syntax error, so every shadowing name is caught at its declaration, which
    # is why this was never silently wrong. It is still wrong, and the fix is to
    # not rewrite what the body owns.
    shadowed = _declared_in_body(body)
    local_members = {k: v for k, v in MEMBERS.items() if k not in shadowed}
    if local_members:
        pattern = re.compile(r"(?<![\w.])(" + "|".join(
            sorted(local_members, key=len, reverse=True)) + r")\b")
        body = substitute_outside_comments(pattern,
                                           lambda m: local_members[m.group(1)], body)
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
            rf"\(const TransientSolveState &state(?:, (.*?))?\)\s*\{{", re.S)
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
        closing = "}  // namespace sisprod::transient\n"
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
            call = f"sisprod::transient::{spec['new_name']}({arguments})"
        else:
            separator = ", " if arguments else ""
            call = (f"sisprod::transient::{spec['new_name']}"
                    f"(transientSolveStateOf(*this){separator}{arguments})")
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
                   else r"\(const TransientSolveState &state")
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
