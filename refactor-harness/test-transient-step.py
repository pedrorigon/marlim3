#!/usr/bin/env python3
"""Check transient histories and ordered finalization with ASan/UBSan.

For a pre-change comparison, supply a directory containing src/ and objects/.
All objects and headers come from that same version because owner layout can
change. This check collects no performance measurements and updates no baseline.
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
    args = parser.parse_args()
    root = Path(__file__).resolve().parent.parent
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    reference = args.reference_dir.resolve() if args.reference_dir else None
    source_root = reference or root
    object_root = reference / "objects" if reference else args.build.resolve() / "CMakeFiles/Marlim3.dir/src"
    objects = sorted(object_root.rglob("*.o"))
    if not objects:
        parser.error("build the executable or supply saved object files first")
    compiler = os.environ.get("CXX", "g++")
    flags = ["-std=c++20", "-O1", "-g", "-ffp-contract=off", "-fopenmp",
             "-fsanitize=address,undefined", "-fno-omit-frame-pointer"]
    definitions = []
    if reference:
        definitions = ["-DMARLIM_REFERENCE_BUILD",
                       f'-DMARLIM_THERMAL_IMPLEMENTATION="{reference / "src/core/SisProdThermal.cpp"}"',
                       f'-DMARLIM_TRANSIENT_IMPLEMENTATION="{reference / "src/core/SisProdTransient.cpp"}"']
    main_object = next(p for p in objects if p.name == "Num4Main.cpp.o")
    executable = output / "checks"
    commands = [
        ["objcopy", "--redefine-sym", "main=marlim_application_main", str(main_object), str(output / "application-main.o")],
        [compiler, *flags, "-fno-access-control", "-Wall", f"-I{source_root / 'src/include'}",
         f"-I{source_root / 'src/thirdparty'}", *definitions, "-c", str(root / "refactor-harness/transient-step-tests.cpp"),
         "-o", str(output / "checks.o")],
        [compiler, *flags, "-o", str(executable), str(output / "checks.o"), str(output / "application-main.o"),
         *(str(p) for p in objects if p.name not in {"Num4Main.cpp.o", "SisProdTransient.cpp.o", "SisProdThermal.cpp.o"}), "-lgfortran"],
    ]
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
    with (output / "checks.log").open("w") as log:
        result = subprocess.run([str(executable), str(output / "state.bin")], env=env, stdout=log, stderr=subprocess.STDOUT)
    print((output / "checks.log").read_text(errors="replace"), end="")
    return result.returncode


if __name__ == "__main__":
    raise SystemExit(main())
