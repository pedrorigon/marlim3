#!/usr/bin/env python3
"""Compare actual steady seed/march bodies using controlled side-effecting dependencies.

This directed ASan/UBSan check isolates dispatch, state publication and iteration
contracts. Real fluid/porous solvers are covered separately by application-level
differentials. No performance measurements or baseline updates are performed.
"""

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import subprocess


def function(text, name):
    match = re.search(r"(?m)^\w+ " + name + r"\(", text)
    if not match:
        raise AssertionError(name)
    opening = text.index("{", match.end())
    clean = re.sub(r'//[^\n]*|/\*[\s\S]*?\*/|"(?:\\.|[^"\\])*"',
                   lambda m: " " * len(m[0]), text)
    depth = 0
    for pos in range(opening, len(clean)):
        depth += (clean[pos] == "{") - (clean[pos] == "}")
        if depth == 0:
            return text[match.start():pos + 1]
    raise AssertionError("unbalanced function " + name)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--reference-dir", type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parent.parent
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    template = (root / "refactor-harness/steady-march-tests.cpp").read_text()
    compiler = os.environ.get("CXX", "g++")
    flags = ["-std=c++20", "-O1", "-g", "-ffp-contract=off", "-fsanitize=address,undefined",
             "-fno-omit-frame-pointer", f"-I{root / 'src/include'}"]
    env = os.environ.copy()
    env.update(ASAN_OPTIONS="detect_leaks=0:halt_on_error=1", UBSAN_OPTIONS="halt_on_error=1:print_stacktrace=1")
    results = {}
    for version, source_root in [("before", args.reference_dir.resolve()), ("after", root)]:
        folder = output / version
        folder.mkdir(exist_ok=True)
        source = (source_root / "src/core/SisProdSteadyState.cpp").read_text()
        old = version == "before"
        names = [] if old else ["refreshSeedCompositionalProperties", "wetGasSeedFractions",
                                "initializeClosedSourceBoundary", "resetProductionMarchInlet"]
        helpers = "\n\n".join(function(source, name) for name in names)
        seed = function(source, "seedFirstCellVoidFraction")
        seed = seed.replace("seedFirstCellVoidFraction(", "observedSeed(", 1)
        private_header = "" if old else (source_root / "src/include/SisProdSteadyStateInternal.h").read_text().replace(
            '#include "SisProdSteadyState.h"', "")
        generated = template.replace("// @PRIVATE_HEADER@", private_header).replace(
            "// @HELPERS@", helpers + "\n\n" + seed).replace(
            "// @MARCH@", function(source, "marchProductionSteady"))
        (folder / "checks.cpp").write_text(generated)
        command = [compiler, *flags, *( ["-DMARLIM_REFERENCE_BUILD"] if old else []),
                   str(folder / "checks.cpp"), "-o", str(folder / "checks")]
        with (folder / "build.log").open("w") as log:
            log.write(shlex.join(command) + "\n"); log.flush()
            p = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
        if p.returncode:
            print((folder / "build.log").read_text()); return p.returncode
        with (folder / "checks.log").open("w") as log:
            p = subprocess.run([str(folder / "checks"), str(folder / "state.bin")],
                               stdout=log, stderr=subprocess.STDOUT, env=env)
        print(version, (folder / "checks.log").read_text(), end="")
        if p.returncode:
            return p.returncode
        data = (folder / "state.bin").read_bytes()
        results[version] = {"sha256": hashlib.sha256(data).hexdigest(), "bytes": len(data)}
    results["bitwise_equal"] = (output / "before/state.bin").read_bytes() == (output / "after/state.bin").read_bytes()
    results["scope"] = "production bodies with controlled dependencies; real solvers checked separately"
    results["performance_measured"] = False
    (output / "summary.json").write_text(json.dumps(results, indent=2) + "\n")
    print("Bitwise equivalence:", results["bitwise_equal"])
    return 0 if results["bitwise_equal"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
