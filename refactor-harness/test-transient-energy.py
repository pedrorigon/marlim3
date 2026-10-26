#!/usr/bin/env python3
"""Check transient thermal boundaries with the actual thermal implementation.

ASan and UBSan instrument the test translation unit, including SisProdThermal.cpp.
The remaining modules use the supplied Release build. Leak detection is disabled
because ownership of the legacy model's thermal allocations is outside these checks.
"""

import argparse
import os
from pathlib import Path
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("build", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parent.parent
    build = args.build.resolve()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    objects = sorted((build / "CMakeFiles/Marlim3.dir/src").rglob("*.o"))
    if not objects:
        parser.error("build the Release executable before running these checks")
    compiler = os.environ.get("CXX", "g++")
    application_main = next(p for p in objects if p.name == "Num4Main.cpp.o")
    renamed_main = output / "application-main.o"
    flags = ["-std=c++20", "-O1", "-g", "-ffp-contract=off", "-fopenmp",
             "-fsanitize=address,undefined", "-fno-omit-frame-pointer"]
    test_object = output / "transient-energy-tests.o"
    executable = output / "transient-energy-tests"
    with (output / "build.log").open("w") as log:
        commands = [
            ["objcopy", "--redefine-sym", "main=marlim_application_main", str(application_main), str(renamed_main)],
            [compiler, *flags, "-fno-access-control", "-Wall", f"-I{root / 'src/include'}",
             f"-I{root / 'src/thirdparty'}", "-c", str(root / "refactor-harness/transient-energy-tests.cpp"),
             "-o", str(test_object)],
            [compiler, *flags, "-o", str(executable), str(test_object), str(renamed_main),
             *(str(p) for p in objects if p.name not in ("Num4Main.cpp.o", "SisProdThermal.cpp.o")), "-lgfortran"],
        ]
        for command in commands:
            result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
            if result.returncode:
                print((output / "build.log").read_text(errors="replace"))
                return result.returncode
    env = os.environ.copy()
    env["ASAN_OPTIONS"] = "detect_leaks=0:halt_on_error=1"
    env["UBSAN_OPTIONS"] = "halt_on_error=1:print_stacktrace=1"
    with (output / "checks.log").open("w") as log:
        result = subprocess.run([str(executable)], env=env, stdout=log, stderr=subprocess.STDOUT)
    print((output / "checks.log").read_text(errors="replace"), end="")
    return result.returncode


if __name__ == "__main__":
    raise SystemExit(main())
