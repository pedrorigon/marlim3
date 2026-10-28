#!/usr/bin/env python3
"""Check production inlet/trend/search bodies with ASan/UBSan and controlled dependencies.

The unchanged domain is compared byte for byte. Independent physical expectations
cover the corrected single-source ratio and the inlet trend. Application-level
regressions exercise the real fluid solvers separately. No baselines or timings.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import subprocess


def block(text, signature):
    start = text.index(signature)
    opening = text.index("{", start)
    clean = re.sub(r'//[^\n]*|/\*[\s\S]*?\*/|"(?:\\.|[^"\\])*"',
                   lambda m: " " * len(m[0]), text)
    depth = 0
    for pos in range(opening, len(clean)):
        depth += (clean[pos] == "{") - (clean[pos] == "}")
        if depth == 0:
            return text[start:pos + 1]
    raise AssertionError("unbalanced " + signature)


def tokens(text):
    return re.sub(r"\s+", "", re.sub(r"//[^\n]*|/\*[\s\S]*?\*/", "", text))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--reference-dir", type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parent.parent
    output = args.output.resolve(); output.mkdir(parents=True, exist_ok=True)
    reference = args.reference_dir.resolve()
    paths = {"steady": "src/core/SisProdSteadyState.cpp", "search": "src/core/SisProdSteadyStateSearch.cpp",
             "trend": "src/core/Leitura.cpp"}
    sources = {version: {name: (folder/path).read_text() for name,path in paths.items()}
               for version,folder in [("before",reference),("after",root)]}
    names = {"source": "void applySteadyMassLiquidInjection(",
             "trend": "if (trendp[trend].qlstTot == 1)",
             "search": "double searchProductionPressureToPressure("}
    bodies = {version: {"source": block(text["steady"], names["source"]),
                        "trend": block(text["trend"], names["trend"]),
                        "search": block(text["search"], names["search"])} for version,text in sources.items()}
    after = sources["after"]["steady"]
    old_march = block(sources["before"]["steady"], "double marchProductionPressureToPressure(")
    observed_march = block(after, "double detail::marchProductionPressureToPressure(")
    observed_march = observed_march[observed_march.index("{"):]
    observed_march = observed_march.replace("    if (stop != nullptr)\n        *stop = {};\n", "", 1)
    observed_march, count = re.subn(
        r"if \(state.cells\[i - 1\].pres <= 1\) \{\s*if \(stop != nullptr\)\s*\*stop = \{i - 1, state.cells\[i - 1\].pres\};\s*return -kSteadyMarchLimit;\s*\}",
        "if (state.cells[i - 1].pres <= 1) return -kSteadyMarchLimit;", observed_march)
    assert count == 1
    assert tokens(old_march[old_march.index("{"):]) == tokens(observed_march)
    wrapper = block(after, "double marchProductionPressureToPressure(")
    assert tokens(wrapper[wrapper.index("{"):]) == tokens("{return detail::marchProductionPressureToPressure(state, dependencies, mchute, nullptr);}")
    old_search = bodies["before"]["search"]
    new_search = bodies["after"]["search"]
    search_unobserved = new_search[:new_search.index("{")+1] + new_search[new_search.index("    state.reverseSteady = 0;"):]
    search_unobserved = re.sub(r"evaluateMarch\((mchute2?), (?:true|false)\)",
                              r"marchProductionPressureToPressure(state.march, dependencies, \1)",search_unobserved)
    search_unobserved = re.sub(r'NumError\(failureMessage\(("[^"\n]+")\)\)',r'NumError(\1)',search_unobserved)
    assert tokens(old_search) == tokens(search_unobserved)
    source_unmodified, count = re.subn(r"\s*else if \(i == 1 &&[\s\S]+?state.cells\[i\].flui.RGO = state.cells\[0\].acsr.injl.FluidoPro.RGO;", "", bodies["after"]["source"])
    assert count == 1 and tokens(source_unmodified) == tokens(bodies["before"]["source"])
    structural = {"pressure_march_calculations_and_order_unchanged": True,
                  "search_calls_conditions_limits_and_returns_unchanged": True,
                  "source_unchanged_except_approved_fallback": True}
    (output/"structural.json").write_text(json.dumps(structural, indent=2)+"\n")
    flags = ["-std=c++20", "-O1", "-g", "-ffp-contract=off", "-fsanitize=address,undefined",
             "-fno-omit-frame-pointer", "-fno-pie", "-no-pie", f"-I{root/'src/include'}"]
    env = {**os.environ, "ASAN_OPTIONS": "detect_leaks=0:halt_on_error=1",
           "UBSAN_OPTIONS": "halt_on_error=1:print_stacktrace=1"}
    compiler = os.environ.get("CXX", "g++")
    def run(label, template, replacements, corrected, expect_failure=False):
        folder = output/label; folder.mkdir(exist_ok=True)
        text = (root/"refactor-harness"/template).read_text()
        for name, body in replacements.items(): text = text.replace("// @"+name+"@", body)
        (folder/"checks.cpp").write_text(text)
        command = [compiler,*flags,f"-DMARLIM_CORRECTED_BUILD={int(corrected)}",str(folder/"checks.cpp"),"-o",str(folder/"checks")]
        with (folder/"build.log").open("w") as stream:
            stream.write(shlex.join(command)+"\n"); stream.flush()
            subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT, check=True)
        arguments = [str(folder/"checks"), str(folder/"unchanged.bin")]
        if template == "steady-inlet-tests.cpp": arguments.append(str(folder/"corrected.bin"))
        with (folder/"checks.log").open("w") as stream:
            result = subprocess.run(arguments, env=env, stdout=stream, stderr=subprocess.STDOUT)
        log = (folder/"checks.log").read_text()
        if expect_failure:
            assert result.returncode != 0, label
        else:
            assert result.returncode == 0, (label,log)
            print(label,log.splitlines()[-1],flush=True)
        return log
    for version in ("before","after"):
        run(version+"-inlet", "steady-inlet-tests.cpp",
            {"LIQUID_SOURCE":bodies[version]["source"], "TOTAL_LIQUID_TREND":bodies[version]["trend"]}, version=="after")
        run(version+"-search", "pressure-search-tests.cpp", {"PRESSURE_SEARCH":bodies[version]["search"]}, version=="after")
    for fixture in ("inlet","search"):
        assert (output/f"before-{fixture}/unchanged.bin").read_bytes() == (output/f"after-{fixture}/unchanged.bin").read_bytes(), fixture
    # Reintroduce each known defect independently; the directed contract must reject it.
    old_ratio = run("control-old-ratio","steady-inlet-tests.cpp",
        {"LIQUID_SOURCE":bodies["before"]["source"], "TOTAL_LIQUID_TREND":bodies["after"]["trend"]},True,True)
    assert "steady inlet contract failed" in old_ratio
    old_trend = run("control-old-trend","steady-inlet-tests.cpp",
        {"LIQUID_SOURCE":bodies["after"]["source"], "TOTAL_LIQUID_TREND":bodies["before"]["trend"]},True,True)
    assert "AddressSanitizer" in old_trend or "runtime error:" in old_trend
    old_diagnostic = run("control-old-diagnostic","pressure-search-tests.cpp",
        {"PRESSURE_SEARCH":bodies["before"]["search"]},True,True)
    assert "pressure search contract failed" in old_diagnostic
    def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
    summary = {"status":"PASS", "structural":structural,"sanitizers":["ASan","UBSan"],
               "unchanged_inlet_domain_bitwise_equal":True,"search_evaluations_and_state_bitwise_equal":True,
               "negative_controls_rejected":["old-ratio","old-trend","old-diagnostic"],
               "after_logs":{name:(output/f"after-{name}/checks.log").read_text().splitlines()[-1] for name in ["inlet","search"]},
               "fixture_sha256":{name:sha(root/"refactor-harness"/name) for name in ["steady-inlet-tests.cpp","pressure-search-tests.cpp",Path(__file__).name]},
               "performance_measured":False,"scope":"Production bodies with controlled dependencies; real solvers covered by native cases"}
    (output/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    print("Directed checks: PASS; three negative controls rejected",flush=True)


if __name__ == "__main__": main()
