/*
 * sisprod2_selftest.cpp
 *
 * Standalone entry point for the new SisProd architecture test suite. It is not
 * part of the Marlim3 build; the pytest test tests/test_sisprod2.py compiles it
 * together with src/core/SisProd.cpp to validate the new implementation in
 * isolation.
 *
 * Keeping the main here (instead of guarded by a macro inside the production
 * translation unit) means SisProd.cpp declares no main and stays clean.
 */

#include "SisProd2.h"

#include <cstdio>
#include <cstring>
#include <iostream>
#include <vector>

namespace sp2 = marlim::sisprod2;

static void printComparison() {
    std::vector<sp2::ComparisonRecord> records = sp2::runReferenceComparison(1e-6);
    std::printf("\nAutomatic comparison (candidate vs reference)\n");
    for (std::size_t i = 0; i < records.size(); ++i) {
        const sp2::ComparisonRecord &record = records[i];
        std::printf("  [%s] %-38s ref=%.6e cand=%.6e err=%.2e\n",
                    record.withinTolerance ? "PASS" : "FAIL", record.name.c_str(),
                    record.reference, record.candidate, record.absoluteError);
    }
}

static void printSampleProfile() {
    sp2::SimContext context;
    sp2::OstreamDiagnostics diagnostics(std::cout);
    sp2::FluidModel fluid(800.0, 1.2, 2.0e-3);
    sp2::ProductionColumn column =
        sp2::buildVerticalWell(fluid, 12.0, 10, 200.0, 0.15, 90.0);

    sp2::SteadyStateRequest request(12.0e5, 3.0, 1e-2, 200, 0);
    sp2::TramoEngine engine(context, diagnostics);
    sp2::SolveResult solved = engine.solveBottomholePressure(column, request);

    std::printf("\nSample steady-state solve\n");
    std::printf("  active implementation: %s\n", sp2::activeSisProdImplementation());
    std::printf("  segments: %zu, total height: %.1f m\n", column.segmentCount(),
                column.totalVerticalHeight());
    std::printf("  bottomhole pressure: %.3f bar (status %s, %d iterations)\n",
                solved.value / 1e5, sp2::toString(solved.status), solved.iterations);

    std::vector<sp2::ProfilePoint> profile;
    column.marchToWellhead(solved.value, context, &profile);
    for (std::size_t i = 0; i < profile.size(); ++i) {
        std::printf("    cell %2zu  P=%8.3f bar  rho=%7.2f kg/m3  v=%6.3f m/s  GOR=%7.2f\n",
                    i, profile[i].pressure / 1e5, profile[i].mixtureDensity,
                    profile[i].velocity, profile[i].gasOilRatio);
    }
}

int main(int argc, char **argv) {
    bool verbose = true;
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--quiet") == 0)
            verbose = false;
    }

    std::printf("SisProd2 architecture test suite\n");
    std::printf("--------------------------------\n");
    const bool testsOk = sp2::runAllTests(verbose);
    printComparison();
    printSampleProfile();
    std::printf("\nRESULT: %s\n", testsOk ? "ALL TESTS PASSED" : "SOME TESTS FAILED");
    return testsOk ? 0 : 1;
}
