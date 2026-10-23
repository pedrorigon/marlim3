#!/usr/bin/env python3
"""Check literal steady helper bodies, call bindings and every enclosing definition."""
import argparse
from pathlib import Path
import re
import subprocess


ROOT = Path(__file__).resolve().parents[1]
REFERENCE = "38903ef"
STATE = "src/core/SisProdSteadyState.cpp"
SEARCH = "src/core/SisProdSteadyStateSearch.cpp"
INTERFACES = (
    "src/include/SisProdSteadyState.h",
    "src/include/SisProdSteadyStateSearch.h",
    "src/include/SisProdConstants.h",
    "src/core/SisProdAdapters.cpp",
)
HELPERS = ("publishSteadyVoidFraction", "publishSteadyComplementaryFraction",
           "estimateInletMixtureVelocity")
LEXEME = re.compile(
    r'"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\'|//[^\n]*|/\*[\s\S]*?\*/'
    r'|[A-Za-z_]\w*|\d+(?:\.\d*)?(?:[eE][+-]?\d+)?|\S')

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from module_dependencies_layout import verified_snapshot, candidate_text



def tokens(text):
    return [m[0] for m in LEXEME.finditer(text) if not m[0].startswith(("//", "/*"))]


def masked(text):
    return re.sub(r'//[^\n]*|/\*[\s\S]*?\*/|"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\'',
                  lambda m: re.sub(r"[^\n]", " ", m[0]), text)


def close_brace(text, opening):
    clean = masked(text)
    depth = 0
    for pos in range(opening, len(clean)):
        depth += (clean[pos] == "{") - (clean[pos] == "}")
        if depth == 0:
            return pos
    raise AssertionError("unbalanced braces")


def function(text, name):
    clean = masked(text)
    matches = list(re.finditer(r"(?m)^[\w:<>]+\s+" + name + r"\s*\(", clean))
    assert len(matches) == 1, (name, "definition count", len(matches))
    start = matches[0].start()
    opening = clean.index("{", matches[0].end())
    return text[start:opening], text[opening + 1:close_brace(text, opening)]


def declarations(text):
    """Account for functions, structs, enums and templates, ignoring namespace placement."""
    clean = masked(text)
    spans = []
    for match in re.finditer(r"(?m)^namespace(?: [\w:]+)? \{", clean):
        opening = clean.index("{", match.start())
        end = close_brace(text, opening)
        spans.extend(((match.start(), opening + 1), (end, end + 1)))
    for start, end in sorted(spans, reverse=True):
        text = text[:start] + text[end:]
    text = re.sub(r"(?m)^#.*", "", text)
    clean = masked(text)
    result = []
    pos = 0
    while (match := re.search(r"\S", clean[pos:])):
        start = pos + match.start()
        opening = clean.index("{", start)
        end = close_brace(text, opening) + 1
        if clean[end:end + 1] == ";":
            end += 1
        result.append(text[start:end])
        pos = end
    return result


def replace_sequence(sequence, needle, replacement):
    result = []
    pos = count = 0
    while pos < len(sequence):
        if sequence[pos:pos + len(needle)] == needle:
            result.extend(replacement)
            pos += len(needle)
            count += 1
        else:
            result.append(sequence[pos])
            pos += 1
    return result, count


def verify(reference, candidate):
    candidate = verified_snapshot(candidate)
    _, no_flow = function(reference[STATE], "finalizeSteadyMassNoFlow")
    start = no_flow.index("state.cells[i].alfini =")
    end = no_flow.index("if (i < state.lastCell) {", start)
    alpha = tokens(no_flow[start:end])
    beta = tokens(no_flow[end:])
    _, bracket = function(reference[SEARCH], "bracketTertiaryPressureToPressureRoot")
    start = bracket.index("double complementaryFraction = state.march.cells[1].bet;")
    end = bracket.index("if (kontaiter > 200", start)
    velocity_text = bracket[start:end]
    velocity = tokens(velocity_text)
    expected_bodies = [alpha, beta, tokens(velocity_text.replace(
        "double velocityGuess = mchute2 *", "return standardFlowRateGuess *"))]
    signatures = [
        "void publishSteadyVoidFraction(const SteadyStateState &state, int i)",
        "void publishSteadyComplementaryFraction(const SteadyStateState &state, int i)",
        "double estimateInletMixtureVelocity(const SteadyStateSearchState &state, double standardFlowRateGuess)",
    ]
    paths = [STATE, STATE, SEARCH]
    for name, signature, body, path in zip(HELPERS, signatures, expected_bodies, paths):
        actual_signature, actual_body = function(candidate[path], name)
        assert tokens(actual_signature) == tokens(signature), (name, "interface")
        assert tokens(actual_body) == body, (name, "literal body")

    counts = [0, 0, 0]
    total = 0
    for path in (STATE, SEARCH):
        clean = masked(candidate[path])
        assert len(re.findall(r"(?m)^namespace \{", clean)) == 1, (path, "private scope")
        assert len(re.findall(r"(?m)^namespace sisprod::steady \{", clean)) == 1
        assert re.findall(r"(?m)^#include.*", reference[path]) == re.findall(
            r"(?m)^#include.*", candidate[path]), (path, "includes")
        expected = sorted(tokens(d) for d in declarations(reference[path]))
        actual = []
        for definition in declarations(candidate[path]):
            signature = tokens(definition[:masked(definition).index("{")])
            if any(name in signature for name in HELPERS):
                continue
            sequence = tokens(definition)
            for index, (name, original) in enumerate(zip(HELPERS, (alpha, beta, velocity))):
                call = ("double velocityGuess = " + name + "(state, mchute2);" if index == 2
                        else name + "(state, i);")
                sequence, count = replace_sequence(sequence, tokens(call), original)
                counts[index] += count
            actual.append(sequence)
        assert sorted(actual) == expected, (path, "enclosing definitions or call bindings")
        total += len(expected)
    assert counts == [15, 12, 4], ("extraction counts", counts)
    for path in INTERFACES:
        assert tokens(candidate[path]) == tokens(reference[path]), (path, "interface or adapter")
    assert total == 109, total


def calibrate(reference, candidate):
    cases = [
        (STATE, "state.cells[i + 1].alfL =", "state.cells[i - 1].alfL =", "alpha neighbour"),
        (STATE, "if (i < state.lastCell)", "if (i <= state.lastCell)", "boundary guard"),
        (STATE, "state.cells[i - 1].betRI = state.cells[i].bet;",
         "state.cells[i - 1].betRI = state.cells[i].alf;", "beta value"),
        (STATE, "state.cells[i].alfini = state.cells[i].alf;", "", "initial transient value"),
        (STATE, "publishSteadyVoidFraction(state, i);",
         "publishSteadyVoidFraction(state, i + 1);", "alpha call binding"),
        (STATE, "publishSteadyComplementaryFraction(state, i);",
         "publishSteadyVoidFraction(state, i);", "publication identity"),
        (STATE, "localtiny * 1e-5", "localtiny * 1e-6", "physical threshold"),
        (SEARCH, "estimateInletMixtureVelocity(state, mchute2);",
         "estimateInletMixtureVelocity(state, mchute);", "velocity call binding"),
        (SEARCH, "double standardFlowRateGuess)", "float standardFlowRateGuess)", "velocity precision"),
        (SEARCH, "rmis * state.march.cells[0].duto.area",
         "rmis * state.march.cells[1].duto.area", "velocity area index"),
        (SEARCH, "MasEspFlu(1., 20.)", "MasEspFlu(1., 15.)", "standard density temperature"),
        (SEARCH, "double rmis = (1 - voidFraction) * rmisL + voidFraction * gasDensityAtGuess;",
         "double rmis = (1 - voidFraction) * rmisL - voidFraction * gasDensityAtGuess;", "mixture formula"),
        (SEARCH, "return -kSteadyConvergenceFailure;", "return kSteadyConvergenceFailure;", "sentinel sign"),
    ]
    for path, before, after, label in cases:
        assert before in candidate[path], (label, "NOT INJECTED")
        mutant = dict(candidate)
        mutant[path] = mutant[path].replace(before, after, 1)
        try:
            verify(reference, mutant)
        except AssertionError:
            print("DETECTED", label)
        else:
            raise AssertionError((label, "MISSED"))
    positive = dict(candidate)
    positive[STATE] = "// Documentation and whitespace control.\n\n" + positive[STATE]
    verify(reference, positive)
    print("POSITIVE CONTROL PASSED")
    print(f"STEADY STATE CALIBRATION PASSED: {len(cases)} negative controls, 1 positive control")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--calibrate", action="store_true")
    args = parser.parse_args()
    paths = (STATE, SEARCH, *INTERFACES)
    reference = {path: subprocess.check_output(
        ["git", "show", f"{REFERENCE}:{path}"], cwd=ROOT, text=True) for path in paths}
    candidate = {path: candidate_text(path) for path in paths}
    verify(reference, candidate)
    print("109 definitions reconstructed; 3 helper interfaces and 31 call bindings checked.")
    if args.calibrate:
        calibrate(reference, candidate)


if __name__ == "__main__":
    main()
