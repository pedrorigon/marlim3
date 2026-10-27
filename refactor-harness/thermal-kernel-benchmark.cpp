#define main thermal_fixture_main
#include "transient-energy-tests.cpp"
#undef main

#include <bit>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <omp.h>

namespace {
std::uint64_t digest(const std::vector<Cel> &cells) {
    std::uint64_t hash = 1469598103934665603ULL;
    for (const auto &cell : cells) {
        for (double value : {cell.temp, cell.dTdt, cell.dTdtIni, cell.fluxcalmed}) {
            hash ^= std::bit_cast<std::uint64_t>(value);
            hash *= 1099511628211ULL;
        }
    }
    return hash;
}
}

int main(int argc, char **argv) {
    if (argc != 7) {
        std::fprintf(stderr, "Usage: %s integrate|advance cells threads substeps sources iterations\n", argv[0]);
        return 2;
    }
    try {
        const std::string kernel = argv[1];
        const int count = std::stoi(argv[2]);
        const int threads = std::stoi(argv[3]);
        const int substeps = std::stoi(argv[4]);
        const bool sources = std::stoi(argv[5]) != 0;
        const int iterations = std::stoi(argv[6]);
        require(count >= 4 && threads > 0 && substeps > 0 && iterations > 0, "positive workload dimensions");
        require(kernel == "integrate" || kernel == "advance", "supported thermal kernel");
        omp_set_dynamic(0);
        Fixture f(0.);
        Cel seed = f.cells[1];
        f.cells.resize(count, seed);
        f.system->celula = f.cells.data();
        f.system->ncel = count - 1;
        f.system->arq.modoDifus3D = 0;
        f.system->tempE = 35.;
        f.globals.solverThreadCount = threads;
        std::vector<TransientTemperatureStep> steps(count);
        for (int i = 0; i < count; ++i) {
            Cel &cell = f.cells[i];
            cell.fluiL = cell.fluiR = &cell.flui;
            cell.temp = cell.tempini = 25. + (i % 23) * 0.2;
            cell.VTemper = 1.;
            cell.QL = cell.QLini = cell.duto.area;
            cell.MC = cell.Mliqini = cell.MliqiniR = 1000. * cell.QL;
            cell.calor.Textern1 = 20.;
            auto &step = steps[i];
            step.initialTemperature = cell.temp;
            step.timeCoefficient = 1. + i;
            step.evolves = i % 19 != 0;
            step.sourceIncrement = sources ? ((i % 7) - 3) * 0.03 : 0.;
            step.sourceDecay = sources ? (i % 5) * 0.07 : 0.;
            step.donor = i % 2 == 0 && i + 1 < count ? i + 1 : i - 1;
            step.courant = step.donor >= 0 && step.evolves ? 0.8 * (substeps - 0.25) : 0.;
        }
        const auto state = f.state();
        const auto transientState = f.system->transientStateOf();
        const auto invoke = [&]() {
            if (kernel == "integrate")
                integrateTransientTemperatureField(state, steps);
            else
                advanceTransientEnergy(state, transientState, 0, 0);
        };
        const auto initialDigest = digest(f.cells);
        for (int i = 0; i < 12; ++i)
            invoke();
        const auto start = std::chrono::steady_clock::now();
        for (int i = 0; i < iterations; ++i)
            invoke();
        const auto finish = std::chrono::steady_clock::now();
        const auto elapsed = std::chrono::duration<double, std::nano>(finish - start).count();
        for (const auto &cell : f.cells)
            require(std::isfinite(cell.temp), "benchmark leaves a finite thermal field");
        std::printf("{\"kernel\":\"%s\",\"cells\":%d,\"threads\":%d,\"substeps\":%d,\"sources\":%s,"
                    "\"iterations\":%d,\"nanoseconds_per_call\":%.9g,\"initial_digest\":\"%016llx\",\"final_digest\":\"%016llx\"}\n",
                    kernel.c_str(), count, threads, substeps, sources ? "true" : "false", iterations,
                    elapsed / iterations, static_cast<unsigned long long>(initialDigest),
                    static_cast<unsigned long long>(digest(f.cells)));
    } catch (const std::exception &error) {
        std::fprintf(stderr, "%s\n", error.what());
        return 1;
    }
    return 0;
}
