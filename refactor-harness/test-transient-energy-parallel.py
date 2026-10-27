#!/usr/bin/env python3
"""Check thermal thread invariance, coverage and coupling with ASan/UBSan.

An optional saved source/object directory enables an exact before/after field
comparison. No performance measurements or reference-output updates are made.
"""

import argparse
import os
from pathlib import Path
import shlex
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("build", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--reference-dir", type=Path)
    parser.add_argument("--source", default="transient-energy-parallel-tests.cpp")
    parser.add_argument("--poisson2d", type=Path)
    parser.add_argument("--poisson3d", type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parent.parent
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    objects = sorted((args.build.resolve() / "CMakeFiles/Marlim3.dir/src").rglob("*.o"))
    if not objects:
        parser.error("build the Release executable first")
    reference = args.reference_dir.resolve() if args.reference_dir else None
    compiler = os.environ.get("CXX", "g++")
    flags = ["-std=c++20", "-O1", "-g", "-ffp-contract=off", "-fopenmp",
             "-fsanitize=address,undefined", "-fno-omit-frame-pointer"]
    includes = [f"-I{root / 'src/include'}", f"-I{root / 'src/thirdparty'}"]
    main_object = reference / "Num4Main.cpp.o" if reference else next(p for p in objects if p.name == "Num4Main.cpp.o")
    commands = [["objcopy", "--redefine-sym", "main=marlim_application_main", str(main_object), str(output / "application-main.o")]]
    definitions = []
    if reference:
        definitions = ["-DMARLIM_REFERENCE_BUILD", f'-DMARLIM_THERMAL_IMPLEMENTATION="{reference / "src/core/SisProdThermal.cpp"}"']
    commands.append([compiler, *flags, *includes, "-fno-access-control", "-Wall", *definitions,
                     "-c", str(root / "refactor-harness" / args.source), "-o", str(output / "checks.o")])
    for name in ["solverPoisson", "solver3DPoisson"]:
        source = (reference or root) / f"src/core/{name}.cpp"
        commands.append([compiler, *flags, *includes, "-c", str(source), "-o", str(output / (name + ".o"))])
    excluded = {"Num4Main.cpp.o", "SisProdThermal.cpp.o", "solverPoisson.cpp.o", "solver3DPoisson.cpp.o"}
    linked = [output / "checks.o", output / "application-main.o", output / "solverPoisson.o", output / "solver3DPoisson.o"]
    linked += [(reference / p.name if reference and (reference / p.name).exists() else p) for p in objects if p.name not in excluded]
    executable = output / "checks"
    commands.append([compiler, *flags, "-o", str(executable), *map(str, linked), "-lgfortran"])
    with (output / "build.log").open("w") as log:
        for command in commands:
            log.write(shlex.join(command) + "\n")
            log.flush()
            result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
            if result.returncode:
                print((output / "build.log").read_text(errors="replace"))
                return result.returncode
    env = os.environ.copy()
    env.update(ASAN_OPTIONS="detect_leaks=0:halt_on_error=1", UBSAN_OPTIONS="halt_on_error=1:print_stacktrace=1",
               OMP_WAIT_POLICY="passive", OMP_DYNAMIC="false")
    for key, value in [("MARLIM_POISSON2D_INPUT", args.poisson2d), ("MARLIM_POISSON3D_INPUT", args.poisson3d)]:
        if value:
            env[key] = str(value.resolve())
    with (output / "checks.log").open("w") as log:
        result = subprocess.run([str(executable), str(output / "fields.bin")], env=env, stdout=log, stderr=subprocess.STDOUT)
    print((output / "checks.log").read_text(errors="replace"), end="")
    return result.returncode


if __name__ == "__main__":
    raise SystemExit(main())
