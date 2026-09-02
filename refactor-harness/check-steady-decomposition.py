#!/usr/bin/env python3
"""Reconstitute the two steady marches from their helpers and compare tokens.

Usage: check-steady-decomposition.py <pre-decomposition-SisProdThermal.cpp>

The pre-decomposition file is the module as it stood at commit 949e23b, before
T071 split the two marches to satisfy SC-004. Obtain it with

    git show 949e23b:src/core/SisProdThermal.cpp > /tmp/pre-decomposition.cpp


A decomposition is exact only if inlining every helper call reproduces the
pre-decomposition body token for token. This expands the four calls in each
march, strips the scaffolding the split introduced (the helper's own
declaration of the value it returns, its return statement, and the struct
unpacking at the call site) and compares against the body as it stood before.
"""
import importlib.util, re, sys
spec = importlib.util.spec_from_file_location("tm", "refactor-harness/thermal-move.py")
tm = importlib.util.module_from_spec(spec); spec.loader.exec_module(tm)

def carve(text, name):
    lines = text.split('\n')
    i = next(k for k, l in enumerate(lines)
             if re.match(rf'^[\w:<>*&]+[\w:<>*&\s]*\b{name}\(', l))
    j = i
    while '{' not in lines[j]: j += 1
    d = 0; k = j
    while True:
        d += lines[k].count('{') - lines[k].count('}')
        if d == 0: break
        k += 1
    return i, k, '\n'.join(lines[i:k + 1])

def helper_body(text, name, drop_leading=None):
    """Helper body without signature, closing brace, or trailing return."""
    i, k, block = carve(text, name)
    lines = block.split('\n')
    b = next(x for x, l in enumerate(lines) if '{' in l)
    body = lines[b + 1:-1]
    # drop everything from the last statement-initial `return` onward
    last = max((x for x, l in enumerate(body) if l.strip().startswith('return')),
               default=None)
    if last is not None:
        body = body[:last]
    while body and not body[-1].strip():
        body.pop()
    if drop_leading and body and body[0].strip() == drop_leading:
        body.pop(0)
    return '\n'.join(body)


cur = open('src/core/SisProdThermal.cpp', encoding='utf-8').read()
base = open(sys.argv[1], encoding='utf-8').read()

# The baseline predates the removal of the fourteen breakpoint anchors -- guards
# whose whole body declares an int, assigns zero to it and stops. They emit no
# code, so removing them cannot move a number, but they are text and the
# comparison is textual. Both sides are normalised, and the count is printed so
# a silent change in the pattern is visible rather than inferred.
cur, cur_anchors = tm.strip_debug_anchors(cur)
base, base_anchors = tm.strip_debug_anchors(base)
print(f"NOTE     ancoras de breakpoint normalizadas: {base_anchors} no baseline, "
      f"{cur_anchors} no modulo atual")

SPECS = [
 ("advanceSteadyTemperature", [
   ("computeSteadySourceTerms",      "TemperatureSourceTerms steadySources", 4, None),
   ("applySteadyAnnulusCoupling",    "annulusResistance = applySteadyAnnulusCoupling", 3, "double annulusResistance = 0.;"),
   ("computeSteadyKineticTerm",      "double kineticTerm = computeSteadyKineticTerm", 3, None),
   ("computeSteadyLatentHeatTerm",   "double latentHeatTerm = computeSteadyLatentHeatTerm", 3, None),
 ]),
 ("advanceReverseSteadyTemperature", [
   ("computeReverseSteadySourceTerms",   "TemperatureSourceTerms reverseSources", 4, None),
   ("applyReverseSteadyAnnulusCoupling", "annulusResistance = applyReverseSteadyAnnulusCoupling", 3, "double annulusResistance = 0.;"),
   ("computeReverseSteadyKineticTerm",   "double kineticTerm = computeReverseSteadyKineticTerm", 3, None),
   ("computeReverseSteadyLatentHeatTerm","double latentHeatTerm = computeReverseSteadyLatentHeatTerm", 3, None),
 ]),
]

rc = 0
for march, helpers in SPECS:
    _, _, march_now = carve(cur, march)
    lines = march_now.split('\n')
    for hname, marker, ncall, drop in helpers:
        idx = next((x for x, l in enumerate(lines) if marker in l), None)
        if idx is None:
            print(f"MISSING  {march}: chamada de {hname}"); rc = 1; break
        lines[idx:idx + ncall] = helper_body(cur, hname, drop).split('\n')
    rebuilt = '\n'.join(lines)
    _, _, original = carve(base, march)
    a, b = tm.tokenize(original), tm.tokenize(rebuilt)
    # compare_tokens accepts a named constant in place of the literal it is
    # PROVEN equal to by a static_assert, and a consistent renaming -- both of
    # which this module has undergone on purpose. It still rejects a changed
    # literal, a reordering, or one variable substituted for another.
    verdict = tm.compare_tokens(a, b)
    if verdict != "differs":
        note = "" if verdict == "exact" else " modulo renames and named constants"
        print(f"OK       {march} reconstituida <- 4 helpers ({len(a)} tokens{note})")
    else:
        pos = tm.first_difference(a, b)
        print(f"DIFFERS  {march} no token {pos}")
        print(f"  antes: {' '.join(a[max(0,pos-6):pos+7])}")
        print(f"  hoje : {' '.join(b[max(0,pos-6):pos+7])}")
        rc = 1
sys.exit(rc)
