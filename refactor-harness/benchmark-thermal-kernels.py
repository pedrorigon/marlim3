#!/usr/bin/env python3
"""Build or explicitly measure only the transient thermal kernels.

The saved directory must contain the pre-change thermal source and application
objects. Results include paired raw samples, exact-output digests and machine
settings. This is not a whole-application performance gate.
"""

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import random
import shlex
import statistics
import subprocess


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("build", type=Path)
    parser.add_argument("saved", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--measure", action="store_true", help="requires explicit authorization for local timing")
    parser.add_argument("--samples", type=int, default=7)
    parser.add_argument("--wait-policy", choices=["active", "passive"], default="passive")
    parser.add_argument("--representative-only", action="store_true")
    parser.add_argument("--kernel", choices=["all", "integrate", "advance"], default="all")
    parser.add_argument("--thread-counts", type=int, nargs="+")
    parser.add_argument("--target-batch-ms", type=float, default=80.)
    parser.add_argument("--iterations", type=int, help="fixed call count for equal thermal trajectories across teams")
    parser.add_argument("--binaries-dir", type=Path, help="reuse already built before/after Release executables")
    args = parser.parse_args()
    root = Path(__file__).resolve().parent.parent
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    if args.measure and (output / "raw.jsonl").exists():
        parser.error("measurement evidence already exists; select a new output directory")
    saved = args.saved.resolve()
    build = args.build.resolve()
    objects = sorted((build / "CMakeFiles/Marlim3.dir/src").rglob("*.o"))
    compiler = os.environ.get("CXX", "g++")
    flags = ["-std=c++20", "-O3", "-DNDEBUG", "-ffp-contract=off", "-funroll-loops", "-fopenmp"]
    binaries = {}
    if args.binaries_dir:
        for version in ["before", "after"]:
            binary = args.binaries_dir.resolve() / version / "thermal-kernel-benchmark"
            if not binary.is_file():
                parser.error(f"missing executable: {binary}")
            binaries[version] = binary
    with (output / "build.log").open("w") as log:
        for version in ([] if args.binaries_dir else ["before", "after"]):
            folder = output / version
            folder.mkdir(exist_ok=True)
            main_object = saved / "Num4Main.cpp.o" if version == "before" else next(p for p in objects if p.name == "Num4Main.cpp.o")
            definitions = [] if version == "after" else [f'-DMARLIM_THERMAL_IMPLEMENTATION="{saved / "src/core/SisProdThermal.cpp"}"']
            linked = [saved / p.name if version == "before" and (saved / p.name).exists() else p
                      for p in objects if p.name not in {"Num4Main.cpp.o", "SisProdThermal.cpp.o"}]
            binary = folder / "thermal-kernel-benchmark"
            commands = [
                ["objcopy", "--redefine-sym", "main=marlim_application_main", str(main_object), str(folder / "application-main.o")],
                [compiler, *flags, "-fno-access-control", f"-I{root / 'src/include'}", f"-I{root / 'src/thirdparty'}",
                 *definitions, "-c", str(root / "refactor-harness/thermal-kernel-benchmark.cpp"), "-o", str(folder / "benchmark.o")],
                [compiler, *flags, "-o", str(binary), str(folder / "benchmark.o"), str(folder / "application-main.o"), *map(str, linked), "-lgfortran"],
            ]
            for command in commands:
                log.write(shlex.join(command) + "\n"); log.flush()
                subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
            binaries[version] = binary
    if not args.measure:
        print("Built the two Release thermal benchmarks; no measurements requested.")
        return
    if args.samples < 3:
        parser.error("use at least three paired samples")
    if args.target_batch_ms <= 0:
        parser.error("the batch duration must be positive")
    if args.iterations is not None and args.iterations <= 0:
        parser.error("the fixed call count must be positive")
    selected, seen = [], set()
    for cpu in sorted(os.sched_getaffinity(0)):
        topology = Path(f"/sys/devices/system/cpu/cpu{cpu}/topology")
        core = ((topology / "physical_package_id").read_text(), (topology / "core_id").read_text())
        if core not in seen:
            seen.add(core); selected.append(cpu)
    os.sched_setaffinity(0, selected)
    env = os.environ.copy()
    env.update(OMP_DYNAMIC="false", OMP_WAIT_POLICY=args.wait_policy, OMP_PROC_BIND="close",
               OMP_PLACES=",".join("{" + str(cpu) + "}" for cpu in selected), OMP_MAX_ACTIVE_LEVELS="1")
    threads = [t for t in (args.thread_counts or [1, 2, 6]) if 0 < t <= len(selected)]
    if not threads:
        parser.error("no requested thread count fits the available physical cores")
    counts = [133, 250] if args.representative_only else [133, 250, 1025, 4096]
    substeps = [1] if args.representative_only else [1, 4, 16]
    cases = [("integrate", count, team, steps, source) for count in counts
             for steps in substeps for source in [0, 1] for team in threads]
    cases += [("advance", count, team, 1, 0) for count in counts for team in threads]
    if args.kernel != "all":
        cases = [case for case in cases if case[0] == args.kernel]
    random.Random(4101).shuffle(cases)
    metadata = {
        "scope": "transient thermal kernels only; not end-to-end simulation time",
        "global_performance_gate": "DEFERRED",
        "authorization": "Pode implementar e fazer uma analise de performance exclusiva para esse trecho!",
        "cpu_model": next(line.split(":", 1)[1].strip() for line in Path("/proc/cpuinfo").read_text().splitlines() if line.startswith("model name")),
        "platform": platform.platform(), "physical_cpu_affinity": selected,
        "compiler": subprocess.check_output([compiler, "--version"], text=True).splitlines()[0],
        "flags": flags, "openmp": {k: env[k] for k in ["OMP_DYNAMIC", "OMP_WAIT_POLICY", "OMP_PROC_BIND", "OMP_PLACES", "OMP_MAX_ACTIVE_LEVELS"]},
        "samples_per_version": args.samples, "target_batch_seconds": args.target_batch_ms / 1000,
        "binary_sha256": {v: sha(p) for v, p in binaries.items()},
        "allocation_and_fixture_setup_timed": False, "kernel_validation_and_commit_timed": True,
        "warmup_calls": 12, "case_order_seed": 4101,
        "fixed_iterations": args.iterations,
        "binaries_reused_from": str(args.binaries_dir.resolve()) if args.binaries_dir else None,
    }
    (output / "metadata.json").write_text(json.dumps(metadata, indent=2, ensure_ascii=False) + "\n")

    def invoke(version, case, iterations):
        command = [str(binaries[version]), *map(str, case), str(iterations)]
        result = subprocess.run(command, cwd=root, env=env, text=True, capture_output=True, check=True)
        return json.loads(result.stdout)

    summaries = []
    with (output / "raw.jsonl").open("w") as raw:
        for case in cases:
            if args.iterations is not None:
                iterations = args.iterations
            else:
                probes = [invoke(version, case, 16) for version in binaries]
                iterations = min(100000, max(8, math.ceil(args.target_batch_ms * 1_000_000 / min(p["nanoseconds_per_call"] for p in probes))))
            measured = {v: [] for v in binaries}
            for pair in range(args.samples):
                order = ["before", "after"] if pair % 2 == 0 else ["after", "before"]
                rows = {}
                for version in order:
                    row = invoke(version, case, iterations)
                    row.update(version=version, pair=pair)
                    raw.write(json.dumps(row) + "\n"); raw.flush()
                    measured[version].append(row["nanoseconds_per_call"])
                    rows[version] = row
                assert rows["before"]["initial_digest"] == rows["after"]["initial_digest"]
                assert rows["before"]["final_digest"] == rows["after"]["final_digest"], (case, rows)
            row = dict(zip(["kernel", "cells", "threads", "substeps", "sources"], case))
            for version, values in measured.items():
                quartiles = statistics.quantiles(values, n=4)
                row[version] = {"median_ns": statistics.median(values), "min_ns": min(values),
                                "q1_ns": quartiles[0], "q3_ns": quartiles[2], "samples_ns": values}
            row["speedup_before_over_after"] = row["before"]["median_ns"] / row["after"]["median_ns"]
            row["change_percent"] = 100 * (row["after"]["median_ns"] / row["before"]["median_ns"] - 1)
            row["output_digests_equal"] = True
            summaries.append(row)
            (output / "summary.json").write_text(json.dumps(summaries, indent=2) + "\n")
            print(case, f'{row["change_percent"]:+.1f}%', flush=True)


if __name__ == "__main__":
    main()
