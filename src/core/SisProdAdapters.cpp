/// The connection between SProd and the modules: the state view each module
/// receives (sisprod::adapters), the solve context that builds the views of a call,
/// and the callbacks through which a module reaches another module's code with them.
#include "SisProd.h"
#include "SisProdSolveContext.h"
#include "SisProdComposition.h"
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

void sisprod::composition::CompositionUpdaters::correctGasSpecificGravity(int i) const {
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
void sisprod::transient::TransientSolveUpdaters::evaluateParaffin() const {
    sisprod::composition::evaluateWaxDeposition(context.composition());
}
void sisprod::transient::TransientSolveUpdaters::connectTubing() const {
    sisprod::gaslift::connectTubing(context.gasLift());
}
void sisprod::transient::TransientSolveUpdaters::marchTransientEnergy(int ciclo, int ciclomax) const {
    sisprod::thermal::advanceTransientEnergy(context.thermal(), ciclo, ciclomax);
}
void sisprod::transient::TransientSolveUpdaters::updateMolarFractions(const ProFlu *fluiRev) const {
    sisprod::composition::transportPhaseMolarFractions(context.composition(), reverseFlowFluid(fluiRev));
}
void sisprod::transient::TransientSolveUpdaters::updateDensities() const {
    sisprod::composition::cacheCellAndFaceDensities(context.composition());
}
void sisprod::transient::TransientSolveUpdaters::updateGasOilRatioAndCo2(const ProFlu *fluiRev) const {
    sisprod::composition::transportBlackOilProperties(context.composition(), reverseFlowFluid(fluiRev));
}
void sisprod::transient::TransientSolveUpdaters::updateTemperatures() const {
    sisprod::thermal::updateDistributedMassTransfer(context.thermal());
}
void sisprod::transient::TransientSolveUpdaters::updateInitialFractions() const {
    sisprod::composition::storePreviousFractionsAndMovePigs(context.composition());
}
void sisprod::transient::TransientSolveUpdaters::updateThermal(int aflu) const {
    sisprod::thermal::updateFlowPartitionTerms(context.thermal(), aflu);
}
void sisprod::transient::TransientSolveUpdaters::solveGasLine() const {
    sisprod::gaslift::solveGasLine(context.gasLift());
}



namespace sisprod::adapters {

// Some views bind a `const T *const &` field to a `T *` member. C++20 binds such
// a reference to the member itself; a temporary there would leave the view
// with a copy that dangles.
constexpr bool bindsToTheMember() {
    int *member = nullptr;
    const int *const &field = member;
    return &field == &member;
}
static_assert(bindsToTheMember());

sisprod::gaslift::GasLiftState gasLiftStateOf(SProd &system, SolveContext &context) {
    return sisprod::gaslift::GasLiftState{
        .gasCells = system.gasLift.cells,
        .cells = system.celula,
        .input = system.arq,
        .globals = system.vg1dSP,
        .lastGasCell = system.gasLift.cellCount,
        .lastCell = system.ncel,
        .gasLiftChokes = system.gasLift.valveChokes,
        .injectionChoke = system.gasLift.injectionChoke,
        .gasValveCellIndices = system.gasLift.gasValveCells,
        .productionValveCellIndices = system.gasLift.productionValveCells,
        .gasSystemMatrix = system.gasLift.systemMatrix,
        .gasRightHandSide = system.gasLift.rightHandSide,
        .tubingCouplingStart = system.gasLift.tubingCouplingStart,
        .tubingCouplingEnd = system.gasLift.tubingCouplingEnd,
        .annulusCouplingStart = system.gasLift.annulusCouplingStart,
        .annulusCouplingEnd = system.gasLift.annulusCouplingEnd,
        .annulusThermalCoupling = system.gasLift.annulusThermalCoupling,
        .thermalSourceDisabled = system.thermalSourceDisabled,
        .injectionPressure = system.gasLift.injectionPressure,
        .injectionTemperature = system.gasLift.injectionTemperature,
        .separatorPressure = system.gasLift.separatorPressure,
        .timeStep = system.dt,
        .interfaceCell = system.gasLift.interfaceCell,
        .interfaceVelocity = system.gasLift.interfaceVelocity,
        .interfaceTimeStep = system.gasLift.interfaceTimeStep,
        .previousInterfaceCell = system.gasLift.previousInterfaceCell,
        .previousInterfaceVelocity = system.gasLift.previousInterfaceVelocity,
        .previousInterfaceTimeStep = system.gasLift.previousInterfaceTimeStep,
        .meanUnloadingFlowRate = system.gasLift.meanUnloadingFlowRate,
        .unloadingAveragingInterval = system.gasLift.unloadingAveragingInterval,
        .maximumUnloadingFlowSamples = system.gasLift.maximumUnloadingFlowSamples,
        .unloadingTimeStepSamples = system.gasLift.unloadingTimeStepSamples,
        .unloadingAveragingWindow = system.gasLift.unloadingAveragingWindow,
        .unloadingSampleCapacity = system.gasLift.unloadingSampleCapacity,
        .temperatureUpdater = {context},
    };
}

sisprod::steady::SteadyStateState steadyStateOf(SProd &system, SolveContext &context) {
    return sisprod::steady::SteadyStateState{
        .cells = system.celula,
        .gasCells = system.gasLift.cells,
        .input = system.arq,
        .globals = system.vg1dSP,
        .lastCell = system.ncel,
        .lastGasCell = system.gasLift.cellCount,
        .injectionChoke = system.gasLift.injectionChoke,
        .surfaceChoke = system.chokeSup,
        .gasValveCellIndices = system.gasLift.gasValveCells,
        .productionValveCellIndices = system.gasLift.productionValveCells,
        .steadyIteration = system.steadyIteration,
        .searchOrigin = system.steadySearch.initialEstimateSearch,
        .convergenceMonitor = system.steadySearch.convergenceMonitor,
        .baseConvergenceMonitor = system.steadySearch.convergenceMonitorReference,
        .annulusDrift = system.derivaAnel,
        .networkCoupled = system.gasLift.annulusThermalCoupling,
        .endNode = system.noextremo,
        .thermalSourceDisabled = system.thermalSourceDisabled,
        .slowHeatTransferThreshold = system.steadySearch.slowHeatTransferThreshold,
        .gasSurfacePressure = system.gasLift.separatorPressure,
        .initialGasPressure = system.gasLift.injectionPressure,
        .initialGasTemperature = system.gasLift.injectionTemperature,
        .outletPressure = system.presfim,
        .timeStep = system.dt,
        .defaultInletTemperature = system.temperatura,
        .casingTemperature = system.steadySearch.reverseFlowTemperature,
        .inletQuality = system.titE,
        .productionFluidCount = system.steadySearch.productionFluidCount,
        .updaters = {context},
    };
}

/// The state a boundary-condition search reads.
///
/// It composes the march state rather than rebuilding it: the two halves read
/// many of the same SProd members, and two bindings of one member could
/// disagree silently.
sisprod::steady::SteadyStateSearchState searchStateOf(SProd &system, SolveContext &context) {
    return sisprod::steady::SteadyStateSearchState{
        .march = steadyStateOf(system, context),
        .holdupGuess = system.steadySearch.holdupGuess,
        .reverseNetworkFluid = system.steadySearch.reverseNetworkFluid,
        .reverseSteady = system.steadySearch.steadyReverseFlow,
    };
}

/// The state one transient step reads.
sisprod::transient::TransientStepState transientStateOf(SProd &system, SolveContext &context) {
    return sisprod::transient::TransientStepState{
        .meanMaximumTimeStep = system.transient.meanMaximumTemperatureRate,
        .meanMaximumPressureChange = system.transient.meanMaximumPressureRate,
        .masterState = system.transient.masterValveState,
        .surfaceChokeOpen = system.transient.surfaceChokeOpen,
        .initiallyOpen = system.transient.previousSurfaceChokeOpen,
        .timeChanged = system.transient.adjustTimeStep,
        .artificialGasThisStep = system.transient.artificialGasThisStep,
        .inletComplementaryFraction = system.betaE,
        .interfaceCell = system.gasLift.interfaceCell,
        .masterCounter = system.transient.masterValveCounter,
        .gasSpecificHeatTable = system.tables.gasSpecificHeatTable,
        .timeStep = system.dt,
        .meanCflTimeStep = system.transient.meanCflTimeStep,
        .totalCflTimeStep = system.transient.totalCflTimeStep,
        .interfaceTimeStep = system.gasLift.interfaceTimeStep,
        .meanSimulationTimeStep = system.transient.meanAcceptedTimeStep,
        .totalSimulationTimeStep = system.transient.totalAcceptedTimeStep,
        .bufferedComplementaryMassSource = system.transient.bufferedComplementaryMassSource,
        .bufferedGasMassSource = system.transient.bufferedGasMassSource,
        .bufferedLiquidMassSource = system.transient.bufferedLiquidMassSource,
        .eventIndex = system.transient.eventIndex,
        .slugCount = system.transient.slugCount,
        .timeStepRestrictionCount = system.transient.restrictedStepsLeft,
        .surfaceChokeMassFlag = system.masChkSup,
        .fullModel = system.transient.fullModel,
        .desperationMoment = system.transient.artificialGasMultiplier,
        .chokeModeChanged = system.transient.chokeModeChanged,
        .multiplier = system.transient.cflSafetyFactor,
        .outletPressure = system.presfim,
        .restart = system.transient.rollbackRequested,
        .timeStepRestricted = system.transient.timeStepRestricted,
        .openTime = system.transient.activeChokeDelay,
        .productionSolution = system.transient.rightHandSide,
        .inletQuality = system.titE,
        .masterRatio0 = system.transient.masterValveRatiosInactive.data(),
        .masterRatio1 = system.transient.masterValveRatiosActive.data(),
        .masterCriticalRatio = system.transient.masterValveCriticalRatios.data(),
        .interfaceVelocity = system.gasLift.interfaceVelocity,
        .masterOpenSchedule = system.transient.masterValveOpenTimes,
        .input = system.arq,
        .initialInterfaceCell = system.gasLift.previousInterfaceCell,
        .cells = system.celula,
        .gasCells = system.gasLift.cells,
        .surfaceChoke = system.chokeSup,
        .cflTimeSteps = system.transient.cflTimeSteps,
        .initialInterfaceTimeStep = system.gasLift.previousInterfaceTimeStep,
        .simulationTimeSteps = system.transient.acceptedTimeSteps,
        .auxiliaryCflTimeStep = system.transient.auxiliaryCflTimeStep,
        .finalAuxiliaryTimeStep = system.transient.auxiliaryAcceptedTimeStep,
        .masterCloseSchedule = system.transient.masterValveCloseTimes,
        .productionFreeTerms = system.flut,
        .gasFreeTerms = system.flutG,
        .branchIndex = system.indTramo,
        .movingMeanFlux = system.transient.movingMeanFlux,
        .stepIndex = system.transient.stepCount,
        .productionMatrix = system.transient.systemMatrix,
        .smallestCellLength = system.transient.smallestCellLength,
        .masterOpenCount = static_cast<int>(system.transient.masterValveOpenTimes.size()),
        .lastCell = system.ncel,
        .lastGasCell = system.gasLift.cellCount,
        .productionCrossSectionCount = system.trends.productionWallProfileCells,
        .masterCloseCount = static_cast<int>(system.transient.masterValveCloseTimes.size()),
        .endNode = system.noextremo,
        .gasSurfacePressure = system.gasLift.separatorPressure,
        .inletPressure = system.presE,
        .movingMeanPressure = system.transient.movingMeanPressure,
        .movingMeanTemperature = system.transient.movingMeanStartTime,
        .maximumTimeStepRates = system.transient.maximumTemperatureRates,
        .maximumPressureRates = system.transient.maximumPressureRates,
        .inletTemperature = system.tempE,
        .reverseQuality = system.networkCoupling.reverseQuality,
        .initialInterfaceVelocity = system.gasLift.previousInterfaceVelocity,
        .globals = system.vg1dSP,
        .updaters = {context},
    };
}

/// The composition module's view of SProd: seventeen members by reference.
sisprod::composition::CompositionState compositionStateOf(SProd &system, SolveContext &context) {
    return sisprod::composition::CompositionState{
        .cells = system.celula,
        .lastCell = system.ncel,
        .input = system.arq,
        .globals = system.vg1dSP,
        .endNode = system.noextremo,
        .inletPressure = system.presE,
        .inletTemperature = system.tempE,
        .inletQuality = system.titE,
        .inletComplementaryFraction = system.betaE,
        .trackGasOilRatio = system.transient.trackGasOilRatio,
        .trackGasGravity = system.transient.trackGasGravity,
        .compositionalRefreshCounter = system.transient.compositionalRefreshCounter,
        .surfaceChokeMassFlag = system.masChkSup,
        .movingPigCount = system.transient.movingPigCount,
        .previousMovingPigCount = system.transient.previousMovingPigCount,
        .scheduledPigCount = static_cast<int>(system.transient.pigReceiverCells.size()),
        .pigReceiverCells = system.transient.pigReceiverCells,
        .updaters = {context},
    };
}


/// The state SolveTrans reads: the step state, composed rather than rebuilt,
/// plus the members only the solve touches.
sisprod::transient::TransientSolveState transientSolveStateOf(SProd &system, SolveContext &context) {
    return sisprod::transient::TransientSolveState{
        .step = transientStateOf(system, context),
        .defaultInletTemperature = system.temperatura,
        .annulusDrift = system.derivaAnel,
        .closingSubtitles = system.saidaSubTextoSis,
        .closingTitles = system.saidaTextoSis,
        .compositionalRefreshCounter = system.transient.compositionalRefreshCounter,
        .fluxHistory = system.transient.fluxSamples,
        .gasCrossSectionCellCounts = system.trends.gasWallProfileCells,
        .gasCrossSectionProfileTimeCounter = system.trends.gasWallProfileTimeIndex,
        .gasCrossSectionTrendBufferedCounts = system.trends.gasWallTrendBufferedCounts,
        .gasCrossSectionTrendCounts = system.trends.gasWallTrendCounts,
        .gasCrossSectionTrendMatrix = system.trends.gasWallTrendMatrices,
        .gasCrossSectionTrendResetTimers = system.trends.gasWallTrendResetTimes,
        .gasProfileTimeCounter = system.trends.gasProfileTimeIndex,
        .gasTrendBufferedCounts = system.trends.gasTrendBufferedCounts,
        .gasTrendCounts = system.trends.gasTrendCounts,
        .gasTrendMatrix = system.trends.gasTrendMatrices,
        .gasTrendResetTimers = system.trends.gasTrendResetTimes,
        .initialGasPressure = system.gasLift.injectionPressure,
        .initialGasSurfacePressure = system.transient.previousSeparatorPressure,
        .initialGasTemperature = system.gasLift.injectionTemperature,
        .initialOpenTime = system.transient.previousActiveChokeDelay,
        .injectionChoke = system.gasLift.injectionChoke,
        .logBuffer = system.tmpLog,
        .logCounter = system.transient.eventLogCount,
        .massTransferModel = system.transient.massTransferModel,
        .minimumCycleTimeStep = system.transient.minimumCycleTimeStep,
        .movingMeanCounter = system.transient.movingMeanWindow,
        .movingMeanVoidFraction = system.transient.movingMeanVoidFraction,
        .networkCoupled = system.gasLift.annulusThermalCoupling,
        .poissonSolver3D = system.coupling3D.poissonSolver3D,
        .pressureHistory = system.transient.pressureSamples,
        .printCounter = system.transient.printCounter,
        .printPassCount = system.trends.trendOutputCount,
        .productionCrossSectionProfileTimeCounter = system.trends.productionWallProfileTimeIndex,
        .productionCrossSectionTrendBufferedCounts = system.trends.productionWallTrendBufferedCounts,
        .productionCrossSectionTrendCounts = system.trends.productionWallTrendCounts,
        .productionCrossSectionTrendMatrix = system.trends.productionWallTrendMatrices,
        .productionCrossSectionTrendResetTimers = system.trends.productionWallTrendResetTimes,
        .productionProfileTimeCounter = system.trends.productionProfileTimeIndex,
        .productionTrendBufferedCounts = system.trends.productionTrendBufferedCounts,
        .productionTrendCounts = system.trends.productionTrendCounts,
        .productionTrendMatrix = system.trends.productionTrendMatrices,
        .productionTrendResetTimers = system.trends.productionTrendResetTimes,
        .startNode = system.noinicial,
        .temperatureHistory = system.transient.timeSamples,
        .totalFlux = system.transient.accumulatedFlux,
        .totalPressure = system.transient.accumulatedPressure,
        .totalVoidFraction = system.transient.accumulatedVoidFraction,
        .trackGasGravity = system.transient.trackGasGravity,
        .trackGasOilRatio = system.transient.trackGasOilRatio,
        .unitCellTimeCounters = system.trends.singleCellOutputTimeIndices,
        .voidFractionHistory = system.transient.voidFractionSamples,
        .updaters = {context},
    };
}



/// The state the five drift-flux closure variants read.
driftflux::coefficient::ClosureState closureStateOf(SProd &system) {
    return driftflux::coefficient::ClosureState{
        .cells = system.celula,
        .lastCell = system.ncel,
        .globals = system.vg1dSP,
        .input = system.arq,
        .selectors = {system.arq.CorreDisper, system.arq.CorreAnular, system.arq.CorreEstrat},
        .gasSurfaceTemperature = system.tGSup,
        .inletVoidFraction = system.alfE,
        .inletComplementaryFraction = system.betaE,
        .inletPressure = system.presE,
        .inletTemperature = system.tempE,
        .steadyIteration = system.steadyIteration,
    };
}

/// The state the source terms of a cell read and write.
sisprod::sources::SourceState sourceStateOf(SProd &system) {
    return sisprod::sources::SourceState{
        .cells = system.celula,
        .input = system.arq,
        .globals = system.vg1dSP,
        .lastCell = system.ncel,
        .steadyMode = system.modoPerm,
        .parallelSecondaryBranch = system.networkCoupling.parallelSecondaryBranch,
        .parallelSecondaryBoundaryCondition = system.networkCoupling.parallelSecondaryBoundaryCondition,
        .parallelSourceCells = system.networkCoupling.parallelSourceCells,
        .parallelSourceProductionLiquid = system.networkCoupling.parallelSourceProductionLiquid,
        .parallelSourceComplementaryLiquid = system.networkCoupling.parallelSourceComplementaryLiquid,
        .parallelSourceGas = system.networkCoupling.parallelSourceGas,
    };
}

sisprod::thermal::ThermalState thermalStateOf(SProd &system, SolveContext &context) {
    return sisprod::thermal::ThermalState{
        .cells = system.celula,
        .gasCells = system.gasLift.cells,
        .input = system.arq,
        .latentHeatTable = system.tables.latentHeatTable.rowPointers,
        .globals = system.vg1dSP,
        .thermalSourceDisabled = system.thermalSourceDisabled,
        .productionNetworkCoupled = system.networkCoupling.secondaryThermalCoupling,
        .primarySectionStart = system.networkCoupling.secondaryCouplingStart,
        .primarySectionEnd = system.networkCoupling.secondaryCouplingEnd,
        .coupledCellIndices = system.coupling3D.coupledIndexCorrections,
        .poissonSolver3D = system.coupling3D.poissonSolver3D,
        .lastCell = system.ncel,
        .surfaceChoke = system.chokeSup,
        .surfaceChokeMassFlag = system.masChkSup,
        .endNode = system.noextremo,
        .gasSurfaceTemperature = system.tGSup,
        .latentHeatEnabled = system.CalcLat,
        .sourceUpdater = {.context = context},
        .fullModel = system.transient.fullModel,
        .massTransferModel = system.transient.massTransferModel,
        .closureUpdater = {.context = context},
        .inletPressure = system.presE,
        .inletTemperature = system.tempE,
        .inletQuality = system.titE,
        .inletVoidFraction = system.alfE,
        .inletComplementaryFraction = system.betaE,
        .evolutionUpdater = {.context = context},
        .surfaceChokeOpen = system.transient.surfaceChokeOpen,
        .defaultInletTemperature = system.temperatura,
        .timeStep = system.dt,
        .minimumCycleTimeStep = system.transient.minimumCycleTimeStep,
        .poisson2DCellIndices = system.transient.poisson2DCells,
        .poisson2DCellCount = static_cast<int>(system.transient.poisson2DCells.size()),
        .steadyIteration = system.steadyIteration,
        .slowHeatTransferThreshold = system.steadySearch.slowHeatTransferThreshold,
        .annulusTubingStart = system.gasLift.tubingCouplingStart,
        .annulusTubingEnd = system.gasLift.tubingCouplingEnd,
        .tubingAnnulusStart = system.gasLift.annulusCouplingStart,
        .productionNetworkHeatCoupled = system.networkCoupling.primaryThermalCoupling,
        .primaryNetworkSectionEnd = system.networkCoupling.primaryCouplingEnd,
        .primaryNetworkSectionStart = system.networkCoupling.primaryCouplingStart,
        .lastGasCell = system.gasLift.cellCount,
        .gasLiftChokes = system.gasLift.valveChokes,
        .gasSurfacePressure = system.gasLift.separatorPressure,
        .outletPressure = system.presfim,
        .surfaceTemperature = system.networkCoupling.chokeOutletTemperature,
        .networkCoupled = system.gasLift.annulusThermalCoupling,
        .tubingAnnulusEnd = system.gasLift.annulusCouplingEnd,
    };
}

}  // namespace sisprod::adapters

namespace sisprod::adapters {

/// Binds the state the trend writers are allowed to read.
///
/// Designated initialisers, not positional: the four series have identical
/// types, so a positional swap would compile in silence and hand a writer
/// another line's buffer. Naming each one makes that a compile error.
///
/// Every field is a reference, never a copy -- the caller advances the counters
/// between the header call and the row call, so a copy would be read at the
/// wrong moment.
trendoutput::TrendState trendStateOf(const SProd &system) {
    return trendoutput::TrendState{
        .input = system.arq,
        .globals = system.vg1dSP,
        .branchIndex = system.indTramo,
        .printPassCount = system.trends.trendOutputCount,
        .production = {system.trends.productionTrendMatrices, system.trends.productionTrendCounts,
                       system.trends.productionTrendBufferedCounts},
        .gasLine = {system.trends.gasTrendMatrices, system.trends.gasTrendCounts, system.trends.gasTrendBufferedCounts},
        .productionCrossSection = {system.trends.productionWallTrendMatrices, system.trends.productionWallTrendCounts,
                                   system.trends.productionWallTrendBufferedCounts},
        .gasLineCrossSection = {system.trends.gasWallTrendMatrices, system.trends.gasWallTrendCounts,
                                system.trends.gasWallTrendBufferedCounts}};
}

}  // namespace sisprod::adapters

namespace sisprod {

SolveContext::SolveContext(SProd &system)
    : system_(system),
      gasLift_(adapters::gasLiftStateOf(system, *this)),
      steady_(adapters::steadyStateOf(system, *this)),
      search_(adapters::searchStateOf(system, *this)),
      transientStep_(adapters::transientStateOf(system, *this)),
      composition_(adapters::compositionStateOf(system, *this)),
      transientSolve_(adapters::transientSolveStateOf(system, *this)),
      closure_(adapters::closureStateOf(system)),
      sources_(adapters::sourceStateOf(system)),
      thermal_(adapters::thermalStateOf(system, *this)),
      trends_(adapters::trendStateOf(system)) {}

}  // namespace sisprod
