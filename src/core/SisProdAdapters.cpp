/// The connection between SProd and the modules: the state view each module
/// receives (sisprod::adapters), the solve context that builds the views of a call,
/// and the callbacks through which a module reaches another module's code with them.
#include "SisProd.h"
#include "SisProdSolveContext.h"
#include "SisProdFluidTracking.h"
#include "SisProdGasLift.h"
#include "SisProdSources.h"
#include "SisProdSteadyState.h"
#include "SisProdSteadyStateSearch.h"
#include "SisProdThermal.h"
#include "SisProdTransient.h"
#include "SisProdTrendOutput.h"

void sisprod::gaslift::GasLiftTemperatureUpdater::dischargeTemperature(
    int cellIndex) const {
    sisprod::thermal::computeDischargeTemperature(context.thermal(), cellIndex);
}

void sisprod::gaslift::GasLiftTemperatureUpdater::computeGasTemperature(
    int cellIndex, double previousTemperature) const {
    sisprod::thermal::computeGasTemperature(context.thermal(), cellIndex, previousTemperature, /*steadyMode=*/0);
}

double sisprod::gaslift::GasLiftTemperatureUpdater::gasLiftDischargeTemperature(
    int valveIndex) const {
    return sisprod::thermal::computeGasLiftDischargeTemperature(context.thermal(), valveIndex);
}

// ------------------------------------------- steady-state callbacks ----
//
// The steady march reaches the drift closures, the source terms, and the gas-line
// and thermal steps through these, with the views of the solve context: its own
// view does not carry their state.
void sisprod::steady::SteadyStateUpdaters::steadyDriftClosure(int cellIndex, double &c0, double &ud) const {
    driftflux::coefficient::steadyState(context.closure(), cellIndex, c0, ud);
}
void sisprod::steady::SteadyStateUpdaters::updateSource(int cellIndex) const {
    sisprod::sources::renewSourceTerms(context.sources(), cellIndex);
}
void sisprod::steady::SteadyStateUpdaters::advanceSteadyTemperature(int cellIndex, int rungeKuttaStage) const {
    sisprod::thermal::advanceSteadyTemperature(context.thermal(), cellIndex, rungeKuttaStage);
}
void sisprod::steady::SteadyStateUpdaters::advanceReverseSteadyTemperature(int cellIndex, int rungeKuttaStage) const {
    sisprod::thermal::advanceReverseSteadyTemperature(context.thermal(), cellIndex, rungeKuttaStage);
}
void sisprod::steady::SteadyStateUpdaters::computeTemperature(int cellIndex, double previousTemperature,
                                                              int steadyMode) const {
    sisprod::thermal::computeTemperature(context.thermal(), cellIndex, previousTemperature, steadyMode);
}
void sisprod::steady::SteadyStateUpdaters::computeGasTemperature(int cellIndex, double previousTemperature,
                                                                 int steadyMode) const {
    sisprod::thermal::computeGasTemperature(context.thermal(), cellIndex, previousTemperature, steadyMode);
}
void sisprod::steady::SteadyStateUpdaters::updateProductionTemperaturePeriphery(int cellIndex) const {
    sisprod::thermal::updateProductionTemperaturePeriphery(context.thermal(), cellIndex);
}
void sisprod::steady::SteadyStateUpdaters::initializeSteadyValveGasFlowRate(int cellIndex) const {
    sisprod::gaslift::initializeSteadyValveGasFlowRate(context.gasLift(), cellIndex);
}
void sisprod::steady::SteadyStateUpdaters::initializeTubingConnectionSteady() const {
    sisprod::gaslift::initializeTubingConnectionSteady(context.gasLift());
}
void sisprod::steady::SteadyStateUpdaters::updateSteadyGasPressure(int cellIndex) const {
    sisprod::gaslift::updateSteadyGasPressure(context.gasLift(), cellIndex);
}
void sisprod::steady::SteadyStateUpdaters::updateSteadyGasTemperature(int cellIndex) const {
    sisprod::gaslift::updateSteadyGasTemperature(context.gasLift(), cellIndex);
}
void sisprod::steady::SteadyStateUpdaters::computeSteadyGasFlowRate(int cellIndex) const {
    sisprod::gaslift::computeSteadyGasFlowRate(context.gasLift(), cellIndex);
}
void sisprod::steady::SteadyStateUpdaters::connectTubing() const {
    sisprod::gaslift::connectTubing(context.gasLift());
}
void sisprod::steady::SteadyStateUpdaters::connectTubingSteady() const {
    sisprod::gaslift::connectTubingSteady(context.gasLift());
}
double sisprod::steady::SteadyStateUpdaters::steadyGasPressureDrop(int cellIndex) const {
    return sisprod::gaslift::steadyGasPressureDrop(context.gasLift(), cellIndex);
}
double sisprod::steady::SteadyStateUpdaters::steadyInjectionPressureDrop(int cellIndex) const {
    return sisprod::gaslift::steadyInjectionPressureDrop(context.gasLift(), cellIndex);
}
// The two searches a march calls; their results are discarded at the call
// site.
void sisprod::steady::SteadyStateUpdaters::searchGasPressureSteadySecondary() const {
    sisprod::steady::searchGasPressureSteadySecondary(context.search());
}
void sisprod::steady::SteadyStateUpdaters::searchGasPressureSteadyTertiary() const {
    sisprod::steady::searchGasPressureSteadyTertiary(context.search());
}

void sisprod::thermal::ThermalSourceUpdater::operator()(int cellIndex) const {
    sisprod::sources::renewSourceTerms(context.sources(), cellIndex);
}

void sisprod::thermal::ThermalClosureUpdater::instantaneous(
    int cellIndex, double &distribution, double &driftVelocity) const {
    driftflux::coefficient::instantaneous(context.closure(), cellIndex, distribution, driftVelocity);
}

void sisprod::thermal::ThermalClosureUpdater::buffered(
    int cellIndex, double &distribution, double &driftVelocity) const {
    driftflux::coefficient::buffered(context.closure(), cellIndex, distribution, driftVelocity);
}

void sisprod::thermal::ThermalClosureUpdater::initialization(
    int cellIndex, double &distribution, double &driftVelocity) const {
    driftflux::coefficient::initialization(context.closure(), cellIndex, distribution, driftVelocity);
}

void sisprod::thermal::ThermalClosureUpdater::bufferedInitialization(
    int cellIndex, double &distribution, double &driftVelocity) const {
    driftflux::coefficient::bufferedInitialization(context.closure(), cellIndex, distribution, driftVelocity);
}

void sisprod::thermal::ThermalEvolutionUpdater::solvePressureVelocityCoupling(
    int cycle) const {
    sisprod::transient::solvePressureVolumeCoupling(context.transientStep(), 0, cycle);
}

void sisprod::thermal::ThermalEvolutionUpdater::renew() const {
    sisprod::transient::updateCells(context.transientStep(), 0);
}

void sisprod::transient::TransientStepUpdaters::advanceGasSubStep() const {
    sisprod::gaslift::advanceGasSubStep(context.gasLift());
}

void sisprod::fluidtracking::FluidTrackingUpdaters::correctGasSpecificGravity(int i) const {
    sisprod::steady::correctGasSpecificGravity(context.steady(), i);
}

namespace {

/// The fluid that flows back into the last cell: the collector's, or a default fluid when the system has none.
const ProFlu &reverseFlowFluid(const ProFlu *collectorFluid) {
    static const ProFlu defaultFluid;
    return collectorFluid != nullptr ? *collectorFluid : defaultFluid;
}

}  // namespace

void sisprod::transient::TransientSolveUpdaters::solveHydrateEnvelopes() const {
    context.system().solveHydrateEnvelopes();
}
double sisprod::transient::TransientSolveUpdaters::searchUnloadingInjectionPressure() const {
    return sisprod::gaslift::searchUnloadingInjectionPressure(context.gasLift());
}
void sisprod::transient::TransientSolveUpdaters::writeProductionTrendHeader(int i, int nrede) const {
    trendoutput::writeProductionTrendHeader(context.trends(), i, nrede);
}
void sisprod::transient::TransientSolveUpdaters::writeProductionTrendRows(int i, int nrede) const {
    trendoutput::writeProductionTrendRows(context.trends(), i, nrede);
}
void sisprod::transient::TransientSolveUpdaters::writeGasTrendHeader(int i, int nrede) const {
    trendoutput::writeGasLineTrendHeader(context.trends(), i, nrede);
}
void sisprod::transient::TransientSolveUpdaters::writeGasTrendRows(int i, int nrede) const {
    trendoutput::writeGasLineTrendRows(context.trends(), i, nrede);
}
void sisprod::transient::TransientSolveUpdaters::writeProductionCrossSectionTrendHeader(int i) const {
    trendoutput::writeProductionCrossSectionTrendHeader(context.trends(), i);
}
void sisprod::transient::TransientSolveUpdaters::writeProductionCrossSectionTrendRows(int i) const {
    trendoutput::writeProductionCrossSectionTrendRows(context.trends(), i);
}
void sisprod::transient::TransientSolveUpdaters::writeGasCrossSectionTrendHeader(int i) const {
    trendoutput::writeGasLineCrossSectionTrendHeader(context.trends(), i);
}
void sisprod::transient::TransientSolveUpdaters::writeGasCrossSectionTrendRows(int i) const {
    trendoutput::writeGasLineCrossSectionTrendRows(context.trends(), i);
}
void sisprod::transient::TransientSolveUpdaters::connectTubing() const {
    sisprod::gaslift::connectTubing(context.gasLift());
}
void sisprod::transient::TransientSolveUpdaters::marchTransientEnergy(int ciclo, int ciclomax) const {
    sisprod::thermal::advanceTransientEnergy(context.thermal(), ciclo, ciclomax);
}
void sisprod::transient::TransientSolveUpdaters::updateMolarFractions(const ProFlu *fluiRev) const {
    sisprod::fluidtracking::transportPhaseMolarFractions(context.fluidTracking(), reverseFlowFluid(fluiRev));
}
void sisprod::transient::TransientSolveUpdaters::updateBlackOilProperties(const ProFlu *fluiRev) const {
    sisprod::fluidtracking::transportBlackOilProperties(context.fluidTracking(), reverseFlowFluid(fluiRev));
}
void sisprod::transient::TransientSolveUpdaters::updateTemperatures() const {
    sisprod::thermal::updateDistributedMassTransfer(context.thermal());
}
void sisprod::transient::TransientSolveUpdaters::updateThermal(int aflu) const {
    sisprod::thermal::updateFlowPartitionTerms(context.thermal(), aflu);
}
void sisprod::transient::TransientSolveUpdaters::solveGasLine() const {
    sisprod::gaslift::solveGasLine(context.gasLift());
}
