#include "SisProd.h"
#include "SimulationContext.h"
#include "../src/core/SisProdThermal.cpp"

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
    near(incoming.cells[1].temp, 29.75, "incoming boundary preserves explicit donor time level");
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

void localSourceHeating() {
    Fixture f;
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
    const double expected = 20.25 + cell.dt * (2. / cell.dx) * 4000. * 50. / capacity;
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
        const double leftLiquidVelocity = left.QL / left.duto.area / (1. - leftVoidFraction);
        const double expected = ((right.MC - right.Mliqini) * rightGasVelocity *
                                 (rightGasVelocity - leftGasVelocity) +
                                 right.Mliqini * rightLiquidVelocity *
                                 (rightLiquidVelocity - leftLiquidVelocity)) / left.dx;
        near(computeKineticTemperatureTerm(state, 0, balance), expected, "kinetic term uses inlet face void fraction");
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
    require(inlet.transport.spatialCoefficient < 0., "external temperature changes local coefficient sign");
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
        near(actual.cells[0].temp, control.cells[0].temp, "field advance agrees with one cell balance");
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
        {"closed endpoint with source enthalpy", localSourceHeating},
        {"network temperature and countercurrent phases", networkTemperatureAndWeightedDirection},
        {"flow reversal and frozen inlet characteristic", reversalAndFrozenCharacteristic},
        {"choke and outlet boundaries", chokeAndOutletBoundaries},
        {"coupled interface and disabled thermal model", coupledInterfaceAndThermalDisable},
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
