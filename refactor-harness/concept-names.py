#!/usr/bin/env python3
"""Cross the glossary against every name the refactoring's own code declares (SC-017).

T104e checked the SProd members the module states rename, and sampled the
locals (flow rate, Reynolds). The sample let a synonym through: nvalv beside
valveCount, both `nvalvgas`, found only later by accident. This makes the
locals side exhaustive.

What counts as a name the refactoring declares
----------------------------------------------
In the nine module .cpp files: every local and parameter declaration, and
every function defined. In the module headers: every struct field and every
function declared. NOT the surface -- fields of celula/celulaG/arq/ProFlu
reached through `.` or `->`, which FR-038 freezes in Portuguese.

Each name is split into words at camelCase and digit boundaries, lower-cased.

The concept table
-----------------
One row per glossary concept: its glossary word, and the other spellings the
same concept could take in English or in leftover Portuguese abbreviation.
For every row the tool prints, per spelling, the names that use it. A row
where two spellings both occur is a CANDIDATE, not a verdict: `previous` is
upstream in space in one place and the previous time level in another, and
only reading the declaration tells. The tool lists; a person decides.

Introduced, not moved
---------------------
By default only names the refactoring INTRODUCED are crossed: those that do
not occur as an identifier anywhere in the baseline SisProd.cpp / SisProd.h
(commit 0f3b64f). A moved local keeps its Portuguese name under FR-039 and is
listed by T104d; it is not a synonym. Two English names for one concept are.
--all crosses every declared name instead.

Usage:
    concept-names.py [--all] [--names]   # --names dumps the names crossed
"""
import subprocess
import collections
import io
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
MODULES = ["DriftFluxClosure", "RootFindingSolvers", "SisProdThermal", "SisProdGasLift",
           "SisProdSteadyState", "SisProdSteadyStateSearch", "SisProdTrendOutput",
           "SisProdComposition", "SisProdTransient"]
HEADERS = MODULES + ["SisProdConstants"]

# glossary concept -> spellings that would be a second name for it
CONCEPTS = [
    ("flowRate",         ["flowrate", "flow", "rate", "vaz", "vazao", "q"]),
    ("pressure",         ["pressure", "pres", "press", "p"]),
    ("temperature",      ["temperature", "temp", "t"]),
    ("density",          ["density", "rho", "masesp"]),
    ("viscosity",        ["viscosity", "visc", "mu"]),
    ("enthalpy",         ["enthalpy", "entalp", "h"]),
    ("surfaceTension",   ["tension", "sigma", "tens"]),
    ("voidFraction",     ["void", "alpha", "alf", "gasfraction"]),
    ("holdup",           ["holdup", "hol", "hold"]),
    ("quality",          ["quality", "tit", "titulo"]),
    ("inclinationAngle", ["inclination", "angle", "ang"]),
    ("diameter",         ["diameter", "diam", "dia"]),
    ("roughness",        ["roughness", "rug"]),
    ("upstream",         ["upstream", "left", "previous", "prev", "montante", "mon"]),
    ("downstream",       ["downstream", "right", "next", "jusante", "jus"]),
    ("guess",            ["guess", "estimate", "chute"]),
    ("source",           ["source", "fonte", "sink"]),
    ("valve",            ["valve", "valv", "vgl"]),
    ("injection",        ["injection", "injector", "inj", "injec"]),
    ("timeStep",         ["timestep", "step", "dt"]),
    ("ratio",            ["ratio", "raz"]),
    ("residual",         ["residual", "error", "erro", "err"]),
    ("tolerance",        ["tolerance", "tol", "epsilon", "eps"]),
    ("iteration",        ["iteration", "iter", "itera"]),
    ("convergence",      ["convergence", "converged", "conv"]),
    ("relaxation",       ["relaxation", "relax", "damping"]),
    ("coupling",         ["coupling", "acop", "coupled"]),
    ("mixture",          ["mixture", "mix", "mis"]),
    ("cellIndex",        ["index", "idx", "ind"]),
    ("steadyState",      ["steady", "perm", "permanent"]),
    ("initialize",       ["initialize", "init", "ini", "initial"]),
    ("compute",          ["compute", "calc", "calculate"]),
    ("gasLine",          ["gasline", "service", "servico"]),
    ("tubing",           ["tubing", "column", "coluna"]),
    ("unloading",        ["unloading", "discharge", "descarga"]),
    ("branch",           ["branch", "tramo"]),
    ("count",            ["count", "number", "num", "n"]),
]

DECL_TYPES = (r'(?:const\s+)?(?:unsigned\s+|long\s+|short\s+)*'
              r'(?:double|int|bool|float|char|auto|size_t|ProFlu|Cel|CelG|FullMtx|choke|'
              r'std::\w+(?:<[^;(){}]*?>)?|[A-Z]\w*(?:::\w+)*)')
DECL = re.compile(r'\b' + DECL_TYPES + r'\s*[*&]*\s+([A-Za-z_]\w*)\s*(?=[,;=\[)(]|\s*$)', re.M)
FUNC_DEF = re.compile(r'^[A-Za-z_][\w:<>,\s*&]*?\b([A-Za-z_]\w*)\s*\([^;{}]*\)\s*(?:const\s*)?\{', re.M)
KEYWORDS = {"return", "if", "else", "for", "while", "switch", "case", "sizeof", "new", "delete",
            "const", "static", "inline", "namespace", "struct", "class", "operator", "void"}


def strip(text):
    return re.sub(r'//[^\n]*|/\*.*?\*/|"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])\'', ' ', text, flags=re.S)


def words(name):
    parts = re.findall(r'[A-Z]+(?=[A-Z][a-z]|\d|\b)|[A-Z]?[a-z]+|[A-Z]+|\d+', name)
    return [p.lower() for p in parts if not p.isdigit()]


def declared():
    names = collections.defaultdict(set)   # name -> {file}
    for m in MODULES:
        path = ROOT / "src/core" / (m + ".cpp")
        code = strip(io.open(path, encoding="utf-8", errors="surrogateescape").read())
        for pat in (DECL, FUNC_DEF):
            for n in pat.findall(code):
                if n not in KEYWORDS:
                    names[n].add(m + ".cpp")
    for h in HEADERS:
        path = ROOT / "src/include" / (h + ".h")
        code = strip(io.open(path, encoding="utf-8", errors="surrogateescape").read())
        for n in DECL.findall(code):
            if n not in KEYWORDS:
                names[n].add(h + ".h")
    return names


def baseline_identifiers():
    ids = set()
    for path in ("src/core/SisProd.cpp", "src/include/SisProd.h"):
        text = subprocess.run(["git", "-C", str(ROOT), "show", "0f3b64f:" + path],
                              capture_output=True, check=True).stdout.decode("utf-8", "surrogateescape")
        ids |= set(re.findall(r'[A-Za-z_]\w*', strip(text)))
    return ids


def main():
    names = declared()
    if "--all" not in sys.argv:
        old = baseline_identifiers()
        names = {n: f for n, f in names.items() if n not in old}
    if "--names" in sys.argv:
        for n in sorted(names):
            print(n, " ".join(sorted(names[n])))
        return 0
    by_word = collections.defaultdict(set)
    for n in names:
        for w in words(n):
            by_word[w].add(n)
    print("names crossed: %d (%s)" % (len(names), "all declared" if "--all" in sys.argv else "introduced by the refactoring"))
    for concept, spellings in CONCEPTS:
        used = [(s, sorted(by_word.get(s, ()))) for s in spellings]
        used = [(s, ns) for s, ns in used if ns]
        flag = "CANDIDATE" if len(used) > 1 else "one spelling" if used else "absent"
        print("\n## %s -- %s" % (concept, flag))
        for s, ns in used:
            shown = ", ".join(ns[:14]) + (" ... (+%d)" % (len(ns) - 14) if len(ns) > 14 else "")
            print("  %-12s %3d  %s" % (s, len(ns), shown))
    return 0


if __name__ == "__main__":
    sys.exit(main())
