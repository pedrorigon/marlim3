#define main thermal_fixture_main
#include "transient-energy-tests.cpp"
#undef main

#include <bit>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <filesystem>
#include <omp.h>

namespace {
void exactField(const std::vector<double> &expected, const std::vector<double> &actual) {
    require(expected.size() == actual.size(), "wall result dimensions");
    for (std::size_t i = 0; i < expected.size(); ++i)
        if (std::bit_cast<std::uint64_t>(expected[i]) != std::bit_cast<std::uint64_t>(actual[i])) {
            std::fprintf(stderr, "wall value %zu: %.17g != %.17g\n", i, expected[i], actual[i]);
            throw std::runtime_error("wall state differs across thread configurations");
        }
}

std::vector<double> advanceBuriedWalls(const std::string &input, const std::vector<int> &indices,
                                      int threads, int activeLevels) {
    omp_set_max_active_levels(activeLevels);
    Fixture f(1.);
    f.globals.solverThreadCount = threads;
    f.globals.simulationTime = 1.;
    f.system->arq.modoDifus3D = 0;
    f.system->transient.poisson2DCells = indices;
    for (int i : indices) {
        auto &wall = f.cells[i].calor;
        wall.poisson2D = solverP(&f.globals, input, 10., 10., 10., 100.,
                                 f.cells[i].temp, 20., 0.10, 0.12, i);
        wall.poisson2D.inicializaTransientePoisson();
        wall.poisson2D.dados.nthrdLocal = threads;
        wall.difus2D = 1;
    }
#ifndef MARLIM_REFERENCE_BUILD
    f.system->allocateEventProfileAndTrendArrays();
    require(f.system->transient.poisson2DCells == indices,
            "production setup classifies every real 2D wall, including the last cell");
#endif
    for (int step = 0; step < 3; ++step) {
        advanceTransientEnergy(f.state(), f.system->transientStateOf(), 0, 0);
        f.globals.simulationTime += f.system->dt;
    }
    std::vector<double> result;
    for (const auto &cell : f.cells) {
        require(std::isfinite(cell.temp), "finite fluid temperature with a real 2D wall");
        for (double value : {cell.temp, cell.tempL, cell.tempR, cell.dTdt, cell.fluxcalmed})
            result.push_back(value);
    }
    for (int i : indices) {
        const auto &solver = f.cells[i].calor.poisson2D;
        require(solver.dados.nthrdLocal == threads, "nested policy preserves the configured wall team");
        for (double value : {solver.dados.qTotal, solver.dados.qAcop, solver.dados.tParede}) {
            require(std::isfinite(value), "finite buried-wall coupling output");
            result.push_back(value);
        }
        for (const auto &element : solver.malha.mlh2d) {
            for (double value : {element.cel2D.tempC, element.cel2D.tempC0, element.cel2D.tempCI,
                                 element.cel2D.gradGreenT[0], element.cel2D.gradGreenT[1], element.TL[0]})
                result.push_back(value);
            for (int j = 0; j < element.cel2D.nvert; ++j)
                result.push_back(element.tempF[j]);
        }
    }
    return result;
}

std::vector<double> advance3D(const std::string &input, int threads) {
    SimulationContext globals;
    globals.poisson3DThreadCount = threads;
    std::vector<DadosGeo> geometry(1);
    solverP3D solver(input, &globals, 1, geometry, 100., 10., 50.);
    solver.inicializaTransientePoisson();
    for (int i = 0; i < 3; ++i) {
        globals.simulationTime = 1. + i * 0.25;
        solver.transientePoisson(0.25);
    }
    require(solver.itera == 5, "3D transient solver preserves its iteration condition");
    std::vector<double> result;
    for (int i = 0; i < solver.dados.CC.nAcop; ++i) {
        require(std::isfinite(solver.dados.qTotal[i]), "finite 3D interface heat flux");
        result.push_back(solver.dados.qTotal[i]);
        result.push_back(solver.dados.tParede[i]);
    }
    for (const auto &element : solver.malha.mlh3d) {
        for (double value : {element.cel3D.tempC, element.cel3D.tempC0, element.cel3D.tempCI, element.TL[0]}) {
            require(std::isfinite(value), "finite 3D element state");
            result.push_back(value);
        }
        for (int j = 0; j < element.cel3D.dim; ++j)
            result.push_back(element.cel3D.gradGreenT[j]);
        for (int j = 0; j < element.cel3D.nvert; ++j)
            result.push_back(element.tempF[j]);
    }
    return result;
}
}

int main(int argc, char **argv) {
    try {
        if (argc == 2)
            pathPrefixoArqSaida = (std::filesystem::path(argv[1]).parent_path() / "mesh-").string();
        const char *input2D = std::getenv("MARLIM_POISSON2D_INPUT");
        const char *input3D = std::getenv("MARLIM_POISSON3D_INPUT");
        require(input2D && input3D, "supply the two real Poisson input paths");
        omp_set_dynamic(0);
        std::vector<double> observed;
        int comparisons = 0;
        for (const std::vector<int> indices : {std::vector<int>{3}, {1, 3}, {0, 1, 2, 3}}) {
            const auto reference = advanceBuriedWalls(input2D, indices, 1, 1);
            observed.insert(observed.end(), reference.begin(), reference.end());
            for (int levels : {1, 2}) {
                const auto parallel = advanceBuriedWalls(input2D, indices, 6, levels);
                exactField(reference, parallel);
                observed.insert(observed.end(), parallel.begin(), parallel.end());
                ++comparisons;
            }
        }
        std::printf("PASS %d real 2D wall comparisons, including outlet, all-wall fields and nested teams\n", comparisons);
        const auto reference3D = advance3D(input3D, 1);
#ifdef MARLIM_REFERENCE_BUILD
        const auto parallel3D = reference3D;
#else
        const auto parallel3D = advance3D(input3D, 6);
        exactField(reference3D, parallel3D);
#endif
        observed.insert(observed.end(), reference3D.begin(), reference3D.end());
        observed.insert(observed.end(), parallel3D.begin(), parallel3D.end());
        std::puts("PASS real 3D transient state; saved reference uses one thread to avoid its known race");
        if (argc == 2) {
            std::ofstream output(argv[1], std::ios::binary);
            output.write(reinterpret_cast<const char *>(observed.data()), observed.size() * sizeof(double));
            require(static_cast<bool>(output), "write wall snapshot");
        }
    } catch (const std::exception &error) {
        std::fprintf(stderr, "FAIL %s\n", error.what());
        return 1;
    }
    return 0;
}
