#include "SisProd.h"
#include "SimulationContext.h"
#ifndef MARLIM_THERMAL_IMPLEMENTATION
#define MARLIM_THERMAL_IMPLEMENTATION "../src/core/SisProdThermal.cpp"
#endif
#include MARLIM_THERMAL_IMPLEMENTATION

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using namespace sisprod::thermal;

void require(bool condition, const char *message) {
    if (!condition)
        throw std::runtime_error(message);
}

void near(double actual, double expected, const char *message) {
    if (!std::isfinite(actual) || std::abs(actual - expected) > 1.e-11 * std::max(1., std::abs(expected))) {
        std::fprintf(stderr, "%s: actual %.17g, expected %.17g\n", message, actual, expected);
        throw std::runtime_error(message);
    }
}

ProFlu fluid(SimulationContext *globals) {
    return ProFlu(globals, 20, 100., 0.7, 0., 1., 20, 10, 40,
                  2, 0, 0.8, 0, 0, 0);
}

TransCal wall(SimulationContext *globals, double temperature) {
    double conductivity[] = {14.};
    double diameter[] = {0.12};
    double heatCapacity[] = {500.};
    double density[] = {7800.};
    double viscosity[] = {1.};
    double expansion[] = {1.e-5};
    int type[] = {0};
    int material[] = {1};
    DadosGeo geometry(0.10, 0., 0., 4.5e-5, 0, 1, conductivity, diameter,
                      heatCapacity, density, viscosity, expansion, type, material);
    int count[] = {1};
    double spacing[] = {0.01};
    double temperatures[] = {temperature, 20.};
    double *layers[] = {temperatures};
    return TransCal(globals, geometry, 1, count, spacing, layers, temperature,
                    20., 20., 0.5, 0.25, 1, 0.25, 0.6, 4000., 1000.,
                    0.001, 0.15, 1100., 900., 0.001, 0, 0, -100.,
                    1., 1., 1., 65., 10.);
}

struct Fixture {
    SimulationContext globals;
    std::unique_ptr<SProd> system = std::make_unique<SProd>();
    std::vector<Cel> cells = std::vector<Cel>(4);
    std::vector<double> flux = std::vector<double>(4, 0.);
    std::vector<double> interfaceTemperature = std::vector<double>(4, 0.);
    std::vector<double> interfaceCoefficient = std::vector<double>(4, 0.);

    Fixture(double velocity = -1.) {
        globals.solverThreadCount = 1;
        system->simulationContext = &globals;
        system->celula = cells.data();
        system->ncel = static_cast<int>(cells.size()) - 1;
        system->dt = 0.25;
        system->transient.minimumCycleTimeStep = system->dt;
        system->transient.fullModel = 1;
        system->thermalSourceDisabled = false;
        system->tempE = 90.;
        system->temperatura = 70.;
        system->tGSup = 40.;
        system->arq.ConContEntrada = 1;
        system->arq.CCPres = {};
        system->arq.master1 = {};
        system->arq.master1.razareaativ = 0.;
        system->arq.chkv = 0;
        system->arq.modoDifus3D = 1;
        system->arq.nacop = static_cast<int>(cells.size());
        system->arq.celAcop.resize(cells.size());
        system->arq.usaTabela = 0;
        system->arq.latente = 0;
        system->arq.condlatente = 1;
        system->arq.limTransMass = 10.;
        system->CalcLat = 0;
        system->coupling3D.poissonSolver3D.dados.qTotal = flux.data();
        system->coupling3D.poissonSolver3D.dados.tInt = interfaceTemperature.data();
        system->coupling3D.poissonSolver3D.dados.hI = interfaceCoefficient.data();
        for (int i = 0; i <= system->ncel; ++i) {
            Cel &cell = cells[i];
            cell.simulationContext = &globals;
            cell.flui = fluid(&globals);
            cell.fluicol = ProFluCol(1000., 0., 0., 1.e-6, 4000.);
            cell.fluiL = &cell.flui;
            cell.fluiR = &cell.flui;
            cell.duto.a = 0.10;
            cell.duto.dia = 0.10;
            cell.duto.area = 0.25 * std::numbers::pi * 0.10 * 0.10;
            cell.duto.teta = 0.;
            cell.duto.rug = 4.5e-5;
            cell.termRug = std::pow(cell.duto.rug / cell.duto.dia / 3.7, 1.11);
            cell.dutoL = cell.dutoR = cell.duto;
            cell.dx = cell.dxL = cell.dxR = 10.;
            cell.dt = system->dt;
            cell.pres = cell.presini = cell.presaux = cell.presL = cell.presR = 65.;
            cell.alf = cell.alfini = cell.alfL = cell.alfR = 0.;
            cell.bet = cell.betini = cell.betL = cell.betR = 1.;
            cell.temp = cell.tempini = i < 3 ? 20. + 10. * i : 40.;
            cell.tempL = cell.tempLini = cell.temp;
            cell.tempR = cell.tempRini = cell.temp;
            cell.QL = cell.QLini = cell.duto.area * velocity;
            cell.QG = cell.QGini = 0.;
            cell.MC = cell.Mliqini = 1000. * cell.QL;
            cell.MliqiniR = cell.Mliqini;
            cell.rpC = cell.rpCi = cell.flui.MasEspLiq(cell.pres, cell.temp);
            cell.rcC = cell.rcCi = 1000.;
            cell.rgC = cell.rgCi = cell.flui.MasEspGas(cell.pres, cell.temp);
            cell.mipC = cell.flui.ViscOleo(cell.pres, cell.temp);
            cell.micC = cell.fluicol.VisFlu(cell.pres, cell.temp);
            cell.migC = cell.flui.ViscGas(cell.pres, cell.temp);
            cell.acsr.tipo = sisprod::kAccessoryNone;
            cell.acsrL = nullptr;
            cell.fontemassLR = cell.fontemassCR = cell.fontemassGR = 0.;
            cell.FonteMudaFase = cell.potTermo = cell.fonteCal = 0.;
            cell.VTemper = cell.VTemperini = 123.;
            cell.dTdLCor = 1.;
            cell.calor = wall(&globals, cell.temp);
            cell.calor.Textern1 = 40.;
            system->arq.celAcop[i].indCel = i;
            system->arq.celAcop[i].FE = 1.;
            system->coupling3D.coupledIndexCorrections.push_back(i);
        }
    }

    auto state() { return system->thermalStateOf(); }

    void refresh() {
        const auto state = this->state();
        for (Cel &cell : cells)
            cell.tempini = cell.temp;
        const auto inlet = prepareTransientInletTemperature(state);
        refreshTemperatureField(state, inlet);
    }

    ~Fixture() {
        system->coupling3D.poissonSolver3D.dados.qTotal = nullptr;
        system->coupling3D.poissonSolver3D.dados.tInt = nullptr;
        system->coupling3D.poissonSolver3D.dados.hI = nullptr;
    }
};

void upwindStencils() {
    Fixture f;
    auto state = f.state();
    TemperatureBalance balance{};
    balance.meanCellLength = 10.;
    double pressureGradient = 0.;
    f.cells[1].VTemper = -1.;
    near(temperatureSpatialGradientOf(state, 1, balance, pressureGradient), 1., "reverse first interior gradient");
    f.cells[0].VTemper = -1.;
    near(temperatureSpatialGradientOf(state, 0, balance, pressureGradient), 1., "reverse inlet gradient");
    f.cells[1].VTemper = 1.;
    f.cells[0].temp = 90.;
    near(temperatureSpatialGradientOf(state, 1, balance, pressureGradient), 1., "forward donor uses saved temperature");
    f.cells[1].VTemper = -1.;
    f.cells[2].dx = 30.;
    near(temperatureSpatialGradientOf(state, 1, balance, pressureGradient), 0.5, "nonuniform right spacing");
    f.cells[1].VTemper = 0.;
    near(temperatureSpatialGradientOf(state, 1, balance, pressureGradient), 0., "stationary advection");
}

void reversedAdvectionAndExternalTemperature() {
    for (double external : {-20., 90.}) {
        Fixture f;
        f.system->tempE = external;
        f.refresh();
        near(f.cells[0].temp, 20.25, "outgoing inlet evolves its stored temperature");
        near(f.cells[1].temp, 30.25, "first interior receives right temperature");
        near(f.cells[2].temp, 40., "uniform downstream segment");
    }
}

void incomingBoundaryAndUniformField() {
    Fixture incoming(1.);
    incoming.refresh();
    near(incoming.cells[0].temp, 90., "incoming external temperature");
    near(incoming.cells[1].temp, 31.5, "incoming boundary uses the temperature supplied at the start of the step");
    for (double velocity : {-1., 1.}) {
        Fixture uniform(velocity);
        uniform.system->tempE = 30.;
        for (Cel &cell : uniform.cells) {
            cell.temp = cell.tempini = 30.;
            cell.calor.Textern1 = 30.;
        }
        uniform.refresh();
        for (const Cel &cell : uniform.cells)
            near(cell.temp, 30., "uniform advection field");
    }
}

void stationaryHeating() {
    Fixture f(0.);
    const auto state = f.state();
    const auto transport = temperatureTransportOf(state, 0);
    f.cells[0].fonteCal = transport.timeCoefficient * 5. * 0.5;
    f.refresh();
    near(f.cells[0].temp, 20.125, "stationary inlet retains nonadvective heating");
}

void knownHeatInputAndPressureWork() {
    for (double dt : {0.1, 0.05, 0.025}) {
        Fixture f(0.);
        f.system->dt = dt;
        for (Cel &cell : f.cells) {
            cell.dt = dt;
            cell.temp = cell.tempini = 20.;
        }
        const auto transport = temperatureTransportOf(f.state(), 2);
        f.cells[2].fonteCal = 20. * transport.timeCoefficient * f.cells[2].dx;
        for (int i = 0; i < std::lround(2. / dt); ++i)
            f.refresh();
        near(f.cells[2].temp, 60., "known heat input satisfies the first law without rate clipping");
    }

    Fixture pressure(0.);
    Cel &cell = pressure.cells[1];
    cell.alf = cell.alfini = 1.;
    cell.pres = cell.presaux = 64.5;
    pressure.cells[2].presaux = cell.pres;
    const auto transport = temperatureTransportOf(pressure.state(), 1);
    const double expected = cell.temp - cell.flui.CalorGasPresMod(cell.presini, cell.temp) *
        cell.rgC * transport.flow.flowArea * (cell.pres - cell.presini) *
        units::kPascalPerKgfPerCm2 / transport.timeCoefficient;
    pressure.refresh();
    near(cell.temp, expected, "pressure work is integrated once over the complete step");
}

void largeCourantAndThreadConsistency() {
    for (double velocity : {-100., 100.}) {
        std::vector<double> serial;
        for (int threads : {1, 2}) {
            Fixture f(velocity);
            f.globals.solverThreadCount = threads;
            const double boundary = velocity > 0. ? 90. : 40.;
            std::vector<double> expected{velocity > 0. ? boundary : 20., 30., 40., 40.};
            const int substeps = 4;
            const double courant = 2.5 / substeps;
            for (int n = 0; n < substeps; ++n) {
                auto next = expected;
                for (int i = 0; i < 4; ++i) {
                    const int donor = velocity > 0. ? i - 1 : i + 1;
                    if (donor >= 0 && donor < 4)
                        next[i] = (1. - courant) * expected[i] + courant * expected[donor];
                }
                expected = next;
            }
            f.refresh();
            for (int i = 0; i < 4; ++i) {
                near(f.cells[i].temp, expected[i], "large-Courant synchronized upwind integration");
                require(f.cells[i].temp >= 20. && f.cells[i].temp <= 90., "pure advection preserves donor bounds");
                if (threads == 1)
                    serial.push_back(f.cells[i].temp);
                else
                    require(f.cells[i].temp == serial[i], "thread count preserves the thermal result");
            }
        }
    }
}

void temporalConvergence() {
    double previousError = 1.;
    for (double dt : {0.2, 0.1, 0.05}) {
        Fixture f(-1.);
        f.system->dt = dt;
        for (Cel &cell : f.cells)
            cell.dt = dt;
        for (int n = 0; n < std::lround(2. / dt); ++n)
            f.refresh();
        const double expectedFirst = 40. - 22. * std::exp(-0.2);
        const double expectedSecond = 40. - 10. * std::exp(-0.2);
        const double error = std::max(std::abs(f.cells[0].temp - expectedFirst),
                                      std::abs(f.cells[1].temp - expectedSecond));
        require(error < 0.55 * previousError, "upwind time refinement converges to the independent cell-cascade solution");
        previousError = error;
    }
}

void incomingGasEnthalpy() {
    for (int accessory : {sisprod::kAccessoryGasInjection, sisprod::kAccessoryLiquidInjection}) {
        Fixture f(0.);
        Cel &cell = f.cells[1];
        cell.acsr.tipo = accessory;
        cell.acsr.injg.temp = cell.acsr.injl.temp = 70.;
        cell.acsr.injg.FluidoPro = cell.acsr.injl.FluidoPro = fluid(&f.globals);
        cell.acsr.injl.fluidocol = cell.fluicol;
        cell.acsr.injl.bet = 1.;
        cell.fontemassGR = 2.;
        const double heatCapacityRate = 2. * cell.flui.CalorGas(cell.presini, cell.temp) / cell.dx;
        const auto source = computeTemperatureSourceTerms(f.state(), 1, cell.dx);
        near(source.heatCapacityRate, heatCapacityRate, "incoming gas carries cp rather than cv");
        near(source.gas, heatCapacityRate * (70. - cell.temp), "gas source enthalpy difference");
        cell.fontemassGR = -2.;
        const auto withdrawal = computeTemperatureSourceTerms(f.state(), 1, cell.dx);
        near(withdrawal.heatCapacityRate, 0., "withdrawal carries the cell's own enthalpy");
        near(withdrawal.gas, 0., "withdrawal has no artificial temperature relaxation");
    }
}

void frozenOperatorSourceAndDomain() {
    TransientTemperatureStep step;
    step.initialTemperature = 20.;
    step.sourceIncrement = 30.;
    step.sourceDecay = 0.6;
    near(integrateTemperatureSource(step, 20., 1.), 70. - 50. * std::exp(-0.6),
         "exact finite-volume source relaxation");
    step.sourceDecay = 1.e-14;
    near(integrateTemperatureSource(step, 20., 1.), 50., "small source decay is evaluated accurately");

    Fixture f(0.);
    std::vector<TransientTemperatureStep> steps(4);
    for (int i = 0; i < 4; ++i)
        steps[i].initialTemperature = 20.;
    steps[1].evolves = true;
    steps[1].sourceIncrement = -400.;
    bool rejected = false;
    try {
        integrateTransientTemperatureField(f.state(), steps);
    } catch (const std::domain_error &) {
        rejected = true;
    }
    require(rejected, "invalid absolute temperature is reported instead of clipped");
}

void localSourceHeating() {
    Fixture f(2. / (1000. * 0.25 * std::numbers::pi * 0.10 * 0.10));
    f.system->arq.ConContEntrada = 0;
    Cel &cell = f.cells[0];
    cell.QL = cell.QLini = 0.;
    cell.acsr.tipo = sisprod::kAccessoryLiquidInjection;
    cell.acsr.injl.temp = 70.;
    cell.acsr.injl.bet = 1.;
    cell.acsr.injl.FluidoPro = fluid(&f.globals);
    cell.acsr.injl.fluidocol = cell.fluicol;
    cell.fontemassCR = 2.;
    const auto state = f.state();
    const double capacity = temperatureTransportOf(state, 0).timeCoefficient;
    const double rate = (2. / cell.dx) * 4000. / capacity;
    const double expected = 70. - 50. * std::exp(-rate * cell.dt);
    f.refresh();
    near(cell.temp, expected, "closed endpoint retains source enthalpy");
    near(cell.QL, 0., "closed external face stays closed");
}

void networkTemperatureAndWeightedDirection() {
    Fixture network(1.);
    network.globals.networkMode = 1;
    network.system->arq.ConContEntrada = 0;
    network.system->temperatura = 57.;
    auto inlet = prepareTransientInletTemperature(network.state());
    require(inlet.prescribed, "incoming network characteristic");
    near(network.cells[0].temp, 57., "network-supplied inlet temperature");

    for (int sign : {-1, 1}) {
        Fixture f;
        for (Cel &cell : f.cells) {
            cell.alf = cell.alfini = cell.alfR = 0.5;
            cell.alfL = 0.2;
            cell.rgC = cell.rgCi = 50.;
            cell.QG = cell.QGini = -sign * 0.5 * cell.duto.area;
            cell.QL = cell.QLini = sign * 0.1 * cell.duto.area;
            cell.Mliqini = 1000. * cell.QL;
            cell.MC = cell.Mliqini + cell.rgC * cell.QG;
        }
        const auto state = f.state();
        const auto transport = temperatureTransportOf(state, 0);
        const double thermalVelocity = transport.spatialCoefficient / transport.timeCoefficient;
        require(sign * thermalVelocity > 0., "energy-weighted characteristic direction");
        require(sign * (f.cells[1].QL + f.cells[1].QG) < 0., "volume direction differs from thermal direction");
        const auto inlet = prepareTransientInletTemperature(state);
        require(inlet.prescribed == (sign > 0), "countercurrent inlet decision");
        auto balance = prepareTemperatureBalance(state, 0, 0, inlet.transport);
        const Cel &left = f.cells[0], &right = f.cells[1];
        const double leftVoidFraction = left.QG > 0. ? left.alfL : left.alf;
        const double rightGasVelocity = right.QG / (right.duto.area * 0.5);
        const double rightLiquidVelocity = right.QL / (right.duto.area * 0.5);
        const double leftGasVelocity = left.QG / left.duto.area / leftVoidFraction;
        const double leftLiquidFraction = 1. - (left.QL > 0. ? left.alfL : left.alf);
        const double leftLiquidVelocity = left.QL / left.duto.area / leftLiquidFraction;
        const double expected = ((right.MC - right.Mliqini) * rightGasVelocity *
                                 (rightGasVelocity - leftGasVelocity) +
                                 right.Mliqini * rightLiquidVelocity *
                                 (rightLiquidVelocity - leftLiquidVelocity)) / left.dx;
        near(computeKineticTemperatureTerm(state, 0, balance), expected, "kinetic term uses inlet face void fraction");
    }
}

void independentPhaseVelocityDonors() {
    Fixture liquid;
    liquid.cells[0].alf = liquid.cells[0].alfini = 0.;
    liquid.cells[1].alf = liquid.cells[1].alfini = 0.;
    liquid.cells[2].alf = liquid.cells[2].alfini = 0.9965;
    for (Cel &cell : liquid.cells) {
        cell.QG = cell.QGini = 0.;
        cell.QL = cell.QLini = 2. * cell.duto.area;
        cell.MC = cell.Mliqini = 1000. * cell.QL;
    }
    auto balance = prepareTemperatureBalance(liquid.state(), 1, 0,
                                            temperatureTransportOf(liquid.state(), 1));
    near(computeKineticTemperatureTerm(liquid.state(), 1, balance), 0.,
         "zero gas flow cannot select a gas-rich donor for forward liquid motion");

    for (int sign : {-1, 1}) {
        Fixture f;
        f.cells[0].alf = f.cells[0].alfini = 0.4;
        f.cells[1].alf = f.cells[1].alfini = 0.2;
        f.cells[2].alf = f.cells[2].alfini = 0.8;
        for (Cel &cell : f.cells) {
            cell.QG = cell.QGini = sign * 0.2 * cell.duto.area;
            cell.QL = cell.QLini = -sign * 0.4 * cell.duto.area;
            cell.Mliqini = 1000. * cell.QL;
            cell.MC = cell.Mliqini + cell.rgCi * cell.QG;
        }
        const double gasRight = sign > 0 ? 0.2 / 0.2 : -0.2 / 0.8;
        const double gasLeft = sign > 0 ? 0.2 / 0.4 : -0.2 / 0.2;
        const double liquidRight = sign > 0 ? -0.4 / 0.2 : 0.4 / 0.8;
        const double liquidLeft = sign > 0 ? -0.4 / 0.8 : 0.4 / 0.6;
        const Cel &right = f.cells[2];
        const double expected = ((right.MC - right.Mliqini) * gasRight * (gasRight - gasLeft) +
                                  right.Mliqini * liquidRight * (liquidRight - liquidLeft)) / f.cells[1].dx;
        auto coefficients = prepareTemperatureBalance(f.state(), 1, 0, temperatureTransportOf(f.state(), 1));
        near(computeKineticTemperatureTerm(f.state(), 1, coefficients), expected,
             "countercurrent phases select independent kinetic-energy donors");
    }
}

void reversalAndFrozenCharacteristic() {
    Fixture reversing(1.);
    reversing.refresh();
    reversing.system->tempE = 10.;
    for (Cel &cell : reversing.cells)
        cell.QL = cell.QLini = cell.MC = cell.Mliqini = 0.;
    reversing.refresh();
    near(reversing.cells[0].temp, 90., "stationary transition retains internal temperature");
    const double oldInlet = reversing.cells[0].temp;
    const double oldRight = reversing.cells[1].temp;
    for (Cel &cell : reversing.cells) {
        cell.QL = cell.QLini = -cell.duto.area;
        cell.MC = cell.Mliqini = 1000. * cell.QL;
    }
    reversing.system->tempE = 120.;
    reversing.refresh();
    near(reversing.cells[0].temp, oldInlet + reversing.system->dt * (oldRight - oldInlet) / 10.,
         "flow reversal changes donor without imposing external temperature");

    Fixture countercurrent;
    const double oldCp = countercurrent.cells[0].flui.CalorGas(65., 20.);
    const double externalCp = countercurrent.cells[0].flui.CalorGas(65., countercurrent.system->tempE);
    require(oldCp != externalCp, "temperature-dependent gas heat capacity");
    const double gasVelocity = oldCp < externalCp ? -0.5 : 0.5;
    const double liquidVelocity = -50. * gasVelocity * (0.5 * (oldCp + externalCp)) / (1000. * 4000.);
    for (Cel &cell : countercurrent.cells) {
        cell.alf = cell.alfini = cell.alfL = cell.alfR = 0.5;
        cell.rgC = cell.rgCi = 50.;
        cell.QG = cell.QGini = gasVelocity * cell.duto.area;
        cell.QL = cell.QLini = liquidVelocity * cell.duto.area;
    }
    const auto state = countercurrent.state();
    const auto inlet = prepareTransientInletTemperature(state);
    require(inlet.prescribed, "characteristic selected before external temperature is imposed");
    require(temperatureTransportOf(state, 0).spatialCoefficient < 0., "external temperature changes the nonlinear coefficient sign");
    require(inlet.transport.spatialCoefficient > 0., "frozen characteristic uses a consistent positive transport coefficient");
    refreshTemperatureField(state, inlet);
    near(countercurrent.cells[0].temp, countercurrent.system->tempE,
         "inlet policy remains fixed during the explicit thermal cycle");
}

void chokeAndOutletBoundaries() {
    Fixture f;
    const auto state = f.state();
    TemperatureBalance balance{};
    balance.meanCellLength = 10.;
    double gradient = 0.;
    f.cells[1].acsr.tipo = sisprod::kAccessoryChoke;
    f.cells[1].acsr.chk.AreaGarg = 0.;
    f.cells[1].VTemper = -1.;
    near(temperatureSpatialGradientOf(state, 1, balance, gradient), 0., "closed choke blocks reverse temperature gradient");
    f.cells[1].acsr.tipo = sisprod::kAccessoryNone;
    f.cells[0].acsr.tipo = sisprod::kAccessoryChoke;
    f.cells[0].acsr.chk.AreaGarg = 0.;
    f.cells[1].VTemper = 1.;
    f.cells[2].tempini = std::numeric_limits<double>::quiet_NaN();
    require(std::isnan(temperatureSpatialGradientOf(state, 1, balance, gradient)), "interior choke zero product retains nonfinite behavior");
    f.cells[2].tempini = f.cells[2].temp;
    f.cells[0].acsr.tipo = sisprod::kAccessoryNone;
    f.cells[2].acsr.tipo = sisprod::kAccessoryChoke;
    f.cells[2].acsr.chk.AreaGarg = 0.;
    f.cells[3].QL = -f.cells[3].QL;
    computeTemperature(state, 3, f.cells[3].tempini);
    require(std::isfinite(f.cells[3].temp), "last cell has no right control volume");

    for (int checkValve : {0, 1}) {
        Fixture outlet;
        outlet.system->arq.chkv = checkValve;
        outlet.system->masChkSup = 1;
        Cel &last = outlet.cells.back();
        Cel &upstream = outlet.cells[outlet.cells.size() - 2];
        upstream.MliqiniR = last.Mliqini = -1.;
        upstream.MR = last.MC = -0.5;
        last.QL = upstream.QLR = std::numeric_limits<double>::quiet_NaN();
        last.QG = 0.5;
        prepareTransientOutletFlow(outlet.state());
        near(last.QL, 0., "outlet liquid-volume flow matches the prohibited liquid mass flow");
        near(upstream.QLR, 0., "outlet liquid-volume flow is published to the upstream cell");
        near(last.QG, checkValve ? 0. : 0.5, "check-valve projection keeps gas volume and mass flow consistent");
    }
}

void coupledInterfaceAndThermalDisable() {
    Fixture f;
    const auto state = f.state();
    const auto inlet = prepareTransientInletTemperature(state);
    advanceCoupledThermalField(state, 1, 1, inlet);
    near(f.interfaceTemperature[0], 20., "3D interface uses internal outgoing temperature");
    near(f.cells[0].calor.dtL, 0., "3D inlet has no left control volume");
    near(f.cells[0].temp, 20.25, "coupled thermal refresh uses reverse upwind");
    Fixture disabled;
    disabled.system->thermalSourceDisabled = true;
    disabled.refresh();
    for (const Cel &cell : disabled.cells)
        near(cell.temp, 40., "disabled thermal solution retains external-temperature rule");

    Fixture single2D(0.);
    single2D.system->transient.poisson2DCells = {2};
    single2D.cells[2].calor.difus2D = 1;
    const double capacity = temperatureTransportOf(single2D.state(), 2).timeCoefficient;
    single2D.cells[2].fonteCal = 2. * capacity * single2D.cells[2].dx;
    single2D.refresh();
    near(single2D.cells[2].temp, 40.5, "a single Poisson-2D cell participates in the energy advance");
}

void outletPressureReconstructionMode() {
    for (int simple : {0, 1}) {
      for (double liquidVelocity : {0., 20.}) {
        Fixture f(0.);
        f.system->arq.MedSimpPresFront = simple;
        Cel &outlet = f.cells[3];
        Cel &upstream = f.cells[2];
        for (Cel &cell : f.cells)
            cell.duto.peri = std::numbers::pi * cell.duto.a;
        outlet.pres = 65.;
        upstream.pres = 67.;
        upstream.MC = upstream.Mliqini = 1000. * upstream.duto.area * liquidVelocity;
        outlet.MC = outlet.Mliqini = 0.;
        outlet.term1 = 1.;
        outlet.term2 = 0.;
        sisprod::transient::updateLastCell(f.system->transientStateOf(), 3, 1);
        require(std::isfinite(outlet.presaux), "outlet reconstructed pressure is finite");
        if (simple || liquidVelocity == 0.) {
            near(outlet.presaux, 66., "last cell honors the requested arithmetic face pressure");
            near(outlet.dpresaux, 0., "simple face pressure has no friction/gravity correction");
        } else {
            require(std::abs(outlet.presaux - 66.) > 1.e-3,
                    "full pressure reconstruction retains the explicit friction correction");
        }
      }
    }
}

void wallAdvanceAndCycleRestoration() {
    for (double velocity : {-1., 1.}) {
        Fixture control(velocity), actual(velocity);
        control.system->arq.modoDifus3D = actual.system->arq.modoDifus3D = 0;
        if (velocity > 0.)
            control.cells[0].temp = control.system->tempE;
        computeTemperature(control.state(), 0, control.cells[0].tempini);
        const auto state = actual.state();
        const auto step = actual.system->transientStateOf();
        advanceTransientEnergy(state, step, 0, 0);
        require(std::isfinite(actual.cells[0].temp), "field advance has a finite temperature");
        for (int j = 0; j < 2; ++j)
            near(actual.cells[0].calor.Tcamada[0][j], control.cells[0].calor.Tcamada[0][j], "wall advances exactly once");
        near(actual.cells[1].tempL, actual.cells[0].temp, "inlet temperature published to neighbour");
        near(actual.cells[1].dTdtL, actual.cells[0].dTdt, "inlet derivative published to neighbour");
    }
    Fixture f;
    f.system->arq.modoDifus3D = 0;
    const double oldWall = f.cells[0].calor.Tcamada[0][0];
    advanceTransientEnergy(f.state(), f.system->transientStateOf(), 0, 1);
    for (int i = 0; i <= f.system->ncel; ++i) {
        near(f.cells[i].temp, i < 3 ? 20. + 10. * i : 40., "intermediate cycle restores temperature");
        near(f.cells[i].tempL, f.cells[i].tempLini, "intermediate cycle restores left temperature");
    }
    near(f.cells[0].calor.Tcamada[0][0], oldWall, "intermediate cycle restores wall temperature");
}

} // namespace

int main() {
    const std::pair<const char *, void (*)()> tests[] = {
        {"upwind stencils and temporal donors", upwindStencils},
        {"reverse advection and external temperature", reversedAdvectionAndExternalTemperature},
        {"incoming boundary and uniform field", incomingBoundaryAndUniformField},
        {"stationary nonadvective heating", stationaryHeating},
        {"known heat input and pressure work", knownHeatInputAndPressureWork},
        {"large Courant and thread consistency", largeCourantAndThreadConsistency},
        {"temporal convergence to a cell-cascade solution", temporalConvergence},
        {"gas-source enthalpy", incomingGasEnthalpy},
        {"source relaxation and invalid absolute temperature", frozenOperatorSourceAndDomain},
        {"closed endpoint with source enthalpy", localSourceHeating},
        {"network temperature and countercurrent phases", networkTemperatureAndWeightedDirection},
        {"independent phase velocity donors", independentPhaseVelocityDonors},
        {"flow reversal and frozen inlet characteristic", reversalAndFrozenCharacteristic},
        {"choke and outlet boundaries", chokeAndOutletBoundaries},
        {"coupled interface and disabled thermal model", coupledInterfaceAndThermalDisable},
        {"outlet pressure reconstruction mode", outletPressureReconstructionMode},
        {"single wall advance and cycle restoration", wallAdvanceAndCycleRestoration},
    };
    try {
        for (const auto &[name, run] : tests) {
            run();
            std::printf("PASS %s\n", name);
        }
    } catch (const std::exception &error) {
        std::fprintf(stderr, "FAIL %s\n", error.what());
        return 1;
    }
    std::puts("All transient energy checks passed.");
    return 0;
}
