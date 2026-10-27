#define main thermal_boundary_fixture_main
#include "transient-energy-tests.cpp"
#undef main
#include <bit>
#include <cstdint>
#include <omp.h>
#include <fstream>

namespace {
std::vector<double> observed;
void exact(double a, double b, const char *message) {
    if (std::bit_cast<std::uint64_t>(a) != std::bit_cast<std::uint64_t>(b)) {
        std::fprintf(stderr, "%s: %.17g != %.17g\n", message, a, b);
        throw std::runtime_error(message);
    }
}

std::vector<double> integrateOnly(int count, int threads, double scale, int repetition) {
    Fixture f(0.);
    f.globals.solverThreadCount = threads;
    f.cells.resize(count);
    f.system->celula = f.cells.data();
    f.system->ncel = count - 1;
    std::vector<TransientTemperatureStep> steps(count);
    for (int i = 0; i < count; ++i) {
        Cel &c = f.cells[i];
        c.dt = 0.25;
        c.VTemper = 1.;
        c.dTdLCor = 1.;
        c.tempini = 20. + (i % 23) * 0.2 + repetition;
        auto &s = steps[i];
        s.initialTemperature = c.tempini;
        s.evolves = (i % 19) != 0;
        s.sourceIncrement = ((i % 7) - 3) * 0.03;
        s.sourceDecay = (i % 5) * 0.07;
        s.timeCoefficient = 1. + i;
        if (i % 2 == 0 && i + 1 < count) s.donor = i + 1;
        else if (i > 0) s.donor = i - 1;
        s.courant = s.donor >= 0 && s.evolves ? scale * (0.1 + (i % 11) * 0.09) : 0.;
    }
    integrateTransientTemperatureField(f.state(), steps);
    std::vector<double> result;
    for (const Cel &c : f.cells) {
        result.push_back(c.temp);
        result.push_back(c.dTdt);
        result.push_back(c.dTdtIni);
    }
    if (!omp_in_parallel())
        observed.insert(observed.end(), result.begin(), result.end());
    return result;
}

std::vector<double> prepareAndIntegrate(int count, int threads, bool radialWall) {
    Fixture f(-1.);
    Cel seed = f.cells[1];
    f.cells.resize(count, seed);
    f.system->celula = f.cells.data();
    f.system->ncel = count - 1;
    f.globals.solverThreadCount = threads;
    f.system->arq.modoDifus3D = radialWall ? 0 : 1;
    f.flux.assign(count, 0.);
    f.interfaceTemperature.assign(count, 0.);
    f.interfaceCoefficient.assign(count, 0.);
    f.system->coupling3D.poissonSolver3D.dados.qTotal = f.flux.data();
    f.system->coupling3D.poissonSolver3D.dados.tInt = f.interfaceTemperature.data();
    f.system->coupling3D.poissonSolver3D.dados.hI = f.interfaceCoefficient.data();
    f.system->arq.nacop = count;
    f.system->arq.celAcop.resize(count);
    f.system->coupling3D.coupledIndexCorrections.resize(count);
    for (int i = 0; i < count; ++i) {
        Cel &c = f.cells[i];
        c.fluiL = &c.flui;
        c.fluiR = &c.flui;
        c.temp = c.tempini = 25. + 0.1 * (i % 31);
        c.calor.Textern1 = 20.;
        c.QL = c.QLini = c.duto.area * (i % 3 == 0 ? -2. : 1.);
        c.MC = c.Mliqini = c.MliqiniR = 1000. * c.QL;
        c.dTdLCor = 1.;
        c.dt = f.system->dt;
        if (i % 17 == 3) {
            c.acsr.tipo = sisprod::kAccessoryLiquidInjection;
            c.acsr.injl.temp = 50.;
            c.acsr.injl.bet = 1.;
            c.acsr.injl.FluidoPro = fluid(&f.globals);
            c.acsr.injl.fluidocol = c.fluicol;
            c.fontemassCR = 0.2;
        }
        f.system->arq.celAcop[i].indCel = i;
        f.system->arq.celAcop[i].FE = 1.;
        f.system->coupling3D.coupledIndexCorrections[i] = i;
    }
    f.refresh();
    std::vector<double> result;
    for (const Cel &c : f.cells) {
        result.push_back(c.temp);
        result.push_back(c.VTemper);
        result.push_back(c.dTdt);
        result.push_back(c.fluxcalmed);
        for (int j = 0; j < 2; ++j) result.push_back(c.calor.Tcamada[0][j]);
    }
    observed.insert(observed.end(), result.begin(), result.end());
    return result;
}

std::vector<double> coupledInterface(int threads, bool partial) {
    Fixture f(-1.);
    f.globals.solverThreadCount = threads;
    const std::vector<int> cells = partial ? std::vector<int>{2, 0} : std::vector<int>{2, 0, 3, 1};
    f.system->arq.nacop = static_cast<int>(cells.size());
    f.system->arq.celAcop.resize(cells.size());
    f.system->coupling3D.coupledIndexCorrections = partial ? std::vector<int>{1, 0} : std::vector<int>{3, 1, 0, 2};
    for (std::size_t i = 0; i < cells.size(); ++i)
        f.system->arq.celAcop[i].indCel = cells[i];
    f.flux = {10., 20., -20., 2.};
    f.system->coupling3D.poissonSolver3D.dados.qTotal = f.flux.data();
    const auto state = f.state();
    const auto inlet = prepareTransientInletTemperature(state);
    advanceCoupledThermalField(state, 1, 1, inlet);
    std::vector<double> result;
    for (const auto &cell : f.cells) {
        result.push_back(cell.temp);
        result.push_back(cell.fluxcalmed);
        result.push_back(cell.calor.dtL);
    }
    result.insert(result.end(), f.interfaceTemperature.begin(), f.interfaceTemperature.end());
    result.insert(result.end(), f.interfaceCoefficient.begin(), f.interfaceCoefficient.end());
    observed.insert(observed.end(), result.begin(), result.end());
    return result;
}

#ifndef MARLIM_REFERENCE_BUILD
void couplingValidation() {
    Fixture f;
    auto &input = f.system->arq;
    auto &solver = f.system->coupling3D.poissonSolver3D;
    solver.dados.CC.nAcop = 4;
    solver.dados.CC.rotuloAcop = new std::string[4]{"d", "b", "a", "c"};
    for (int i = 0; i < 4; ++i)
        input.celAcop[i].rotulo = std::string(1, 'a' + i);
    require(buildThermalCouplingIndices(input, solver, 3) == std::vector<int>({2, 1, 3, 0}),
            "coupling destinations follow labels rather than cell order");
    const auto rejected = [&]() {
        try { (void)buildThermalCouplingIndices(input, solver, 3); }
        catch (const std::invalid_argument &) { return; }
        throw std::runtime_error("invalid 3D coupling was accepted");
    };
    input.celAcop[1].indCel = 0; rejected(); input.celAcop[1].indCel = 1;
    input.celAcop[1].indCel = -1; rejected(); input.celAcop[1].indCel = 1;
    input.celAcop[1].indCel = 4; rejected(); input.celAcop[1].indCel = 1;
    input.celAcop[1].rotulo = "a"; rejected(); input.celAcop[1].rotulo = "b";
    input.celAcop[1].rotulo = "missing"; rejected(); input.celAcop[1].rotulo = "b";
    solver.dados.CC.rotuloAcop[0] = "a"; rejected(); solver.dados.CC.rotuloAcop[0] = "d";
    solver.dados.CC.nAcop = 3; rejected(); solver.dados.CC.nAcop = 4;
    input.nacop = 0; rejected(); input.nacop = 4;
    std::puts("PASS reordered 3D mapping and eight invalid configurations");
}

void complete2DCoverage() {
    for (int threads : {1, 6}) {
        Fixture f(0.);
        f.globals.solverThreadCount = threads;
        f.system->thermalSourceDisabled = true;
        f.cells.back().calor.difus2D = 1;
        f.system->allocateEventProfileAndTrendArrays();
        require(f.system->transient.poisson2DCells == std::vector<int>({3}),
                "last 2D control volume belongs to the preparation partition");
        f.system->allocateEventProfileAndTrendArrays();
        require(f.system->transient.poisson2DCells == std::vector<int>({3}),
                "rebuilding the classification does not duplicate a wall");
        for (double ambient : {43., 58.}) {
            f.cells.back().calor.Textern1 = ambient;
            f.refresh();
            exact(f.cells.back().temp, ambient, "last 2D operator never reuses stale cache data");
        }
    }
    std::puts("PASS last-cell 2D classification, repeated setup and refreshed cache");
}
#endif
}
int main(int argc, char **argv) {
    try {
        omp_set_dynamic(0);
        int cases = 0;
        for (int threads : {1, 2, 3, 6, 8}) {
            int actual = 0;
#pragma omp parallel num_threads(threads)
            {
#pragma omp single
                actual = omp_get_num_threads();
            }
            require(actual == threads, "requested OpenMP team is active");
        }
        for (int count : {1, 4, 17, 133, 250, 1025}) {
            for (double scale : {0.5, 2., 13.}) {
                for (int repetition = 0; repetition < 3; ++repetition) {
                    auto reference = integrateOnly(count, 1, scale, repetition);
                    for (int threads : {2, 3, 6, 8}) {
                        auto parallel = integrateOnly(count, threads, scale, repetition);
                        require(reference.size() == parallel.size(), "field dimensions agree");
                        for (std::size_t i = 0; i < reference.size(); ++i)
                            exact(reference[i], parallel[i], "frozen thermal operator is bitwise invariant across thread counts");
                        ++cases;
                    }
                }
            }
        }
        std::printf("PASS %d integration comparisons, mixed donors, sources and substep counts\n", cases);
        int preparationCases = 0;
        for (int count : {4, 17, 133, 250}) {
            for (bool radial : {false, true}) {
                auto reference = prepareAndIntegrate(count, 1, radial);
                for (int threads : {2, 3, 6, 8}) {
                    auto parallel = prepareAndIntegrate(count, threads, radial);
                    for (std::size_t i = 0; i < reference.size(); ++i)
                        exact(reference[i], parallel[i], "properties, wall and thermal field are bitwise invariant");
                    ++preparationCases;
                }
            }
        }
        std::printf("PASS %d preparation/wall/field comparisons\n", preparationCases);
        for (bool partial : {false, true}) {
            const auto reference = coupledInterface(1, partial);
            for (int threads : {2, 6}) {
                const auto parallel = coupledInterface(threads, partial);
                for (std::size_t i = 0; i < reference.size(); ++i)
                    exact(reference[i], parallel[i], "exclusive reordered 3D interface is thread invariant");
            }
        }
        std::puts("PASS four complete/partial reordered 3D interface comparisons");
        omp_set_max_active_levels(2);
        std::vector<std::vector<double>> concurrent(2);
#pragma omp parallel num_threads(2)
        {
            const int caller = omp_get_thread_num();
            concurrent[caller] = integrateOnly(133, 3, 2., caller + 10);
        }
        for (int caller = 0; caller < 2; ++caller) {
            auto reference = integrateOnly(133, 1, 2., caller + 10);
            for (std::size_t i = 0; i < reference.size(); ++i)
                exact(reference[i], concurrent[caller][i], "caller-local caches isolate independent nested invocations");
        }
        std::puts("PASS two independent concurrent callers with nested OpenMP teams");
#ifndef MARLIM_REFERENCE_BUILD
        couplingValidation();
        complete2DCoverage();
#endif
        if (argc == 2) {
            std::ofstream output(argv[1], std::ios::binary);
            output.write(reinterpret_cast<const char *>(observed.data()), observed.size() * sizeof(double));
            require(static_cast<bool>(output), "write field snapshot");
        }
        std::puts("Correctness checks collect no performance measurements.");
    } catch (const std::exception &e) {
        std::fprintf(stderr, "FAIL %s\n", e.what());
        return 1;
    }
    return 0;
}
