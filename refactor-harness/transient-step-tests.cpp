#define main transient_thermal_fixture_main
#include "transient-energy-tests.cpp"
#undef main
#ifndef MARLIM_TRANSIENT_IMPLEMENTATION
#define MARLIM_TRANSIENT_IMPLEMENTATION "../src/core/SisProdTransient.cpp"
#endif
#include MARLIM_TRANSIENT_IMPLEMENTATION

#include <bit>
#include <cstdint>
#include <fstream>
#include <sstream>

namespace {
using namespace sisprod::transient;
std::vector<double> observed;

void exactValue(double actual, double expected, const char *message) {
    if (std::bit_cast<std::uint64_t>(actual) != std::bit_cast<std::uint64_t>(expected))
        throw std::runtime_error(message);
}

void observeHistories(const sisprod::TransientRun &run, bool full) {
    for (double value : {run.movingMeanWindow, run.accumulatedPressure, run.accumulatedFlux,
                         run.accumulatedVoidFraction, run.movingMeanPressure, run.movingMeanFlux,
                         run.movingMeanVoidFraction, static_cast<double>(run.timeSamples.size())})
        observed.push_back(value);
    if (full) {
        for (std::size_t i = 0; i < run.timeSamples.size(); ++i)
            for (double value : {run.pressureSamples[i], run.fluxSamples[i], run.voidFractionSamples[i], run.timeSamples[i]})
                observed.push_back(value);
    }
}

void movingHistorySequence(int seedCount) {
    Fixture fixture(0.);
    auto &run = fixture.system->transient;
    run.movingMeanStartTime = 60.;
    const auto state = fixture.system->transientSolveStateOf();
    for (int i = 0; i < seedCount; ++i) {
        const double dt = 0.001 + (i % 7) * 0.003;
        const double pressure = (50. + (i % 23) * 0.03) * dt;
        const double flux = ((i % 11) - 5) * 0.02 * dt;
        const double fraction = 0.4 * dt;
        run.timeSamples.push_back(dt);
        run.pressureSamples.push_back(pressure);
        run.fluxSamples.push_back(flux);
        run.voidFractionSamples.push_back(fraction);
        run.movingMeanWindow += dt;
        run.accumulatedPressure += pressure;
        run.accumulatedFlux += flux;
        run.accumulatedVoidFraction += fraction;
    }
    observeHistories(run, true);
    for (int iteration = 0; iteration < 1200; ++iteration) {
        fixture.globals.simulationTime = iteration == 0 ? 59. : iteration == 1 ? 60. :
            iteration == 2 ? std::nextafter(60., 100.) : 60. + iteration * 0.05;
        fixture.system->dt = 0.001 + (iteration % 19) * 0.007;
        fixture.system->presfim = 60. + (iteration % 29) * 0.03;
        Cel &outlet = fixture.cells[fixture.cells.size() - 2];
        outlet.pres = fixture.system->presfim;
        outlet.temp = 30. + (iteration % 7) * 0.2;
        outlet.alf = 0.1 + (iteration % 7) * 0.05;
        outlet.Mliqini = ((iteration % 11) - 5) * 0.01;
        outlet.MC = outlet.Mliqini + ((iteration % 5) - 2) * 0.002;
        const auto oldSize = run.timeSamples.size();
        updateTransientMovingMeans(state);
        const auto expectedSize = oldSize + (fixture.globals.simulationTime > 60. ? 0 : 1);
        require(run.timeSamples.size() == expectedSize, "remove exactly one sample only after the threshold");
        require(run.pressureSamples.size() == expectedSize && run.fluxSamples.size() == expectedSize &&
                run.voidFractionSamples.size() == expectedSize, "weighted series remain aligned");
        observeHistories(run, false);
    }
    observeHistories(run, true);
}

void networkFinalization(int mode, int wallCount, bool screen) {
    Fixture fixture(0.);
    auto &system = *fixture.system;
    fixture.globals.networkMode = mode;
    fixture.globals.simulationTime = screen ? 30. : 0.;
    fixture.globals.networkSimulationTime = 37.125;
    system.arq.tfinal = 1000000.;
    system.arq.saidaTela = screen ? 1 : 0;
    system.arq.logevento.resize(1);
    system.arq.logevento[0].instante = 1.e9;
    system.arq.ntela = 0;
    system.arq.calculaEnvelope = 0;
    system.transient.stepCount = 23;
    system.transient.printCounter = 7;
    system.transient.movingMeanStartTime = 60.;
    system.arq.modoDifus3D = 0;
    for (int i = 0; i < wallCount; ++i) {
        const int cellIndex = static_cast<int>(fixture.cells.size()) - 1 - i;
        system.transient.poisson2DCells.insert(system.transient.poisson2DCells.begin(), cellIndex);
        auto &wall = fixture.cells[cellIndex].calor;
        wall.difus2D = 1;
        auto &solver = wall.poisson2D;
        solver.malha.nele = 1;
        solver.malha.mlh2d.emplace_back();
        solver.malha.mlh2d[0].cel2D.nvert = 0;
        solver.malha.mlh2d[0].cel2D.tempC = 31. + cellIndex;
        solver.malha.mlh2d[0].cel2D.tempC0 = -10.;
        solver.simulationContext = &fixture.globals;
        solver.dados.acop = 1;
        solver.dados.tInt = 50. + cellIndex;
        solver.dados.tInt0 = -10.;
        solver.dados.temp.parserieImp = 1;
        solver.dados.temp.tempoImp = new double[1]{1.e9};
    }
    system.transient.artificialGasThisStep.push_back({0, 0, 0.01});
    system.SolveTrans();
    exactValue(fixture.globals.simulationTime, 37.125, "network caller retains control of time advancement");
    require(system.transient.stepCount == 24 && fixture.globals.timeStepIndex == 23,
            "published time-step index precedes accepted-step increment");
    require(system.transient.printCounter == 8, "print counter increments after common reporting");
    require(system.trends.trendOutputCount == (screen ? 0 : 1),
            "trend decision uses branch time before the network clock is synchronized");
    require(system.transient.artificialGasThisStep.empty(), "accepted gas is reported before ending the step");
    for (int cellIndex : system.transient.poisson2DCells) {
        const auto &solver = fixture.cells[cellIndex].calor.poisson2D;
        exactValue(solver.dados.tInt0, solver.dados.tInt, "wall interface advances its previous time level");
        exactValue(solver.malha.mlh2d[0].cel2D.tempC0, 31. + cellIndex, "last and other buried walls are committed");
        require(solver.kontaTempo == 0, "unscheduled wall output leaves its counter unchanged");
        observed.push_back(solver.dados.tInt0);
        observed.push_back(solver.malha.mlh2d[0].cel2D.tempC0);
    }
    for (const auto &cell : fixture.cells)
        for (double value : {cell.pres, cell.temp, cell.alf, cell.bet, cell.rpC, cell.rgC, cell.rcC,
                             cell.rpL, cell.rgL, cell.rcL, cell.MC, cell.Mliqini, cell.dTdt, cell.dTdtIni})
            observed.push_back(value);
    observeHistories(system.transient, true);
}

#ifndef MARLIM_REFERENCE_BUILD
void couplingModeSelection() {
    for (int scenario = 0; scenario < 3; ++scenario) {
        Fixture fixture(0.);
        auto state = fixture.system->transientStateOf();
        state.input.master1.posic = 0;
        state.cells[0].acsr.chk.AreaGarg = state.cells[0].duto.area;
        fixture.system->chokeSup.AreaGarg = state.cells[0].duto.area;
        state.valveAreaRatiosBeforeUpdate[0] = scenario == 1 ? 0.5 : 1.;
        state.masterCriticalRatio[0] = -1.;
        for (auto &cell : fixture.cells)
            cell.presini = scenario == 2 ? cell.pres - 100. : cell.pres;
        const int expected = scenario == 0 ? 1 : 0;
        const int cycle = prepareCouplingMode(state, 0);
        require(cycle == expected && state.input.cicloAcopTerm == expected && state.fullModel == expected,
                "coupling mode follows unchanged, changed and high-rate valve conditions");
        exactValue(state.masterCriticalRatio[0], 0.5, "critical ratios reset before the policy");
    }
}
#endif
}

int main(int argc, char **argv) {
    try {
        std::ostringstream discardedScreenOutput;
        struct RestoreOutput {
            std::streambuf *original;
            ~RestoreOutput() { std::cout.rdbuf(original); }
        } restore{std::cout.rdbuf(discardedScreenOutput.rdbuf())};
        for (int count : {0, 1, 513, 20000})
            movingHistorySequence(count);
        for (int mode : {1, 2, 3})
            for (int walls : {0, 1, 2, 4})
                for (bool screen : {false, true})
                    networkFinalization(mode, walls, screen);
#ifndef MARLIM_REFERENCE_BUILD
        couplingModeSelection();
#endif
        for (double value : observed)
            require(std::isfinite(value), "finite observed state");
        if (argc == 2) {
            std::ofstream output(argv[1], std::ios::binary);
            output.write(reinterpret_cast<const char *>(observed.data()), observed.size() * sizeof(double));
            require(static_cast<bool>(output), "write observed transient state");
        }
        std::puts("PASS four variable-step history sequences, 4,800 updates, threshold and single-removal rules");
        std::puts("PASS 24 network finalizations, old/new clocks, screen modes and zero/last/multiple/all buried walls");
        std::puts("PASS accepted artificial-gas clearing and step/output counters");
        std::printf("Observed %zu finite double values. No performance measurements.\n", observed.size());
    } catch (const std::exception &error) {
        std::fprintf(stderr, "FAIL %s\n", error.what());
        return 1;
    }
    return 0;
}
