/*
 * SisProd.cpp
 *
 *  Created on: 21 Dec 2016
 *      Author: Eduardo
 */
#include "SisProd.h"
#include "FA_Hidratos.h"
#include "FA_Hidratos_Servico.h"
#include "SisProdConstants.h"
#include "SisProdGasLift.h"
#include "SisProdSources.h"
#include "SisProdSteadyState.h"
#include "SisProdSteadyStateSearch.h"
#include "SisProdTransient.h"
#include "SisProdFluidTracking.h"
#include "SisProdSolveContext.h"
#include "SisProdThermal.h"
#include "SisProdTrendOutput.h"
#include <stdexcept>

using namespace sisprod;

SProd::SProd(string inputFile, string logFile, tipoValidacaoJson_t jsonValidation,
             tipoSimulacao_t simulationType, varGlob1D *globals, int dynamicTableOverride, int blocking,
             int reversedGeometry, double *extraSourceLengths, int *extraSourceCells, int extraSourceCount,
             int steadyBranch)
    : arq(inputFile, logFile, jsonValidation, simulationType, reversedGeometry, globals, steadyBranch),
      vg1dSP(globals) {
    if (dynamicTableOverride >= 0)
        arq.tabelaDinamica = dynamicTableOverride;
    assemble({.bloq = blocking, .previousReverseLiquidFraction = 0.}, extraSourceLengths, extraSourceCells,
             extraSourceCount);
}

SProd::SProd() = default;

SProd::SProd(Ler &parsedInput, const CarriedState &carried) : vg1dSP(parsedInput.vg1dSP) {
    arq.copiaSemJson(parsedInput);
    assemble(carried);
}

void SProd::loadSection(const SProd &other) {
    if (celula != nullptr)
        throw std::logic_error("SProd::loadSection: this system was already assembled");
    arq = other.arq;
    vg1dSP = other.vg1dSP;
    networkCoupling.takeParallelRole(other.networkCoupling);
    assemble(other.carriedState());
}

SProd::CarriedState SProd::carriedState() const {
    return {.noextremo = noextremo,
            .noinicial = noinicial,
            .derivaAnel = derivaAnel,
            .bloq = bloq,
            .previousReverseLiquidFraction = networkCoupling.previousReverseLiquidFraction,
            .previousReverseQuality = networkCoupling.previousReverseQuality,
            .minimumCycleTimeStep = transient.minimumCycleTimeStep};
}

void SProd::renova(int explicitScheme) {
    SolveContext context(*this);
    transient::updateCells(context.transientStep(), explicitScheme);
}

void SProd::renovaBuffer() {
    SolveContext context(*this);
    transient::updateBufferFromSolution(context.transientStep());
}

void SProd::renovaBufferCego() {
    SolveContext context(*this);
    transient::updateBufferFromCells(context.transientStep());
}

void SProd::renovaTemp() {
    SolveContext context(*this);
    thermal::updateDistributedMassTransfer(context.thermal());
}

void SProd::renovaterm(int inflowMode) {
    SolveContext context(*this);
    thermal::updateFlowPartitionTerms(context.thermal(), inflowMode);
}

void SProd::renovatermAfluFim() {
    SolveContext context(*this);
    thermal::updateOutletFlowPartitionTerms(context.thermal());
}

void SProd::renovatermColIni() {
    SolveContext context(*this);
    thermal::updateInletFlowPartitionTerms(context.thermal());
}

void SProd::calcCCpres(double reverseQuality, double reverseVoidFraction, double reverseLiquidFraction) {
    SolveContext context(*this);
    transient::applyOutletPressureCondition(context.transientStep(), reverseQuality, reverseVoidFraction,
                                            reverseLiquidFraction);
}

void SProd::calcCCBuffer(double reverseQuality, double reverseVoidFraction, double reverseLiquidFraction) {
    SolveContext context(*this);
    transient::applyOutletBufferCondition(context.transientStep(), reverseQuality, reverseVoidFraction,
                                          reverseLiquidFraction);
}

void SProd::determinaDT(int explicitScheme) {
    SolveContext context(*this);
    transient::computeTimeStep(context.transientStep(), explicitScheme);
}

void SProd::atenuaDtMax() {
    SolveContext context(*this);
    transient::dampMaximumTimeStep(context.transientStep());
}

void SProd::avaliaVariaDpDt(int explicitScheme) {
    SolveContext context(*this);
    transient::evaluatePressureRateOfChange(context.transientStep(), explicitScheme);
}

void SProd::aberturaVal1() {
    SolveContext context(*this);
    transient::valveOpeningHigh(context.transientStep());
}

void SProd::aberturaVal0() {
    SolveContext context(*this);
    transient::valveOpeningLow(context.transientStep());
}

void SProd::restringeDTporValv() {
    SolveContext context(*this);
    transient::restrictTimeStepByValve(context.transientStep());
}

void SProd::solveLinGas() {
    SolveContext context(*this);
    gaslift::solveGasLine(context.gasLift());
}

void SProd::EvoluiFrac(double reverseVoidFraction, double reverseLiquidFraction, int cycle) {
    SolveContext context(*this);
    transient::evolveFractions(context.transientStep(), reverseVoidFraction, reverseLiquidFraction, cycle);
}

void SProd::ReiniEvolFrac0() {
    SolveContext context(*this);
    transient::restartFractionEvolutionInitial(context.transientStep());
}

void SProd::ReiniEvolFrac() {
    SolveContext context(*this);
    transient::restartFractionEvolution(context.transientStep());
}

void SProd::AtualizaPig() {
    SolveContext context(*this);
    transient::updatePig(context.transientStep());
}

void SProd::SolveAcopPV(int explicitScheme, int cycle) {
    SolveContext context(*this);
    transient::solvePressureVolumeCoupling(context.transientStep(), explicitScheme, cycle);
}

void SProd::marchaEnergTrans(int cycle, int maximumCycle) {
    SolveContext context(*this);
    thermal::advanceTransientEnergy(context.thermal(), cycle, maximumCycle);
}

void SProd::atualizaMiniTab() {
    SolveContext context(*this);
    transient::refreshFluidMiniTable(context.transientStep());
}

void SProd::atualizaCC1() {
    SolveContext context(*this);
    transient::refreshInletCondition(context.transientStep());
}

void SProd::SolveTrans(double reverseQuality, double reverseVoidFraction, double reverseLiquidFraction, int network,
                       const ProFlu *collectorFluid) {
    SolveContext context(*this);
    transient::solveTransientStep(context.transientSolve(), reverseQuality, reverseVoidFraction, reverseLiquidFraction,
                                  network, collectorFluid);
}

void SProd::ImprimeTrendP(int trend, int network) {
    SolveContext context(*this);
    trendoutput::writeProductionTrendRows(context.trends(), trend, network);
}

void SProd::ImprimeTrendPCab(int trend, int network) {
    SolveContext context(*this);
    trendoutput::writeProductionTrendHeader(context.trends(), trend, network);
}

void SProd::ImprimeTrendG(int trend, int network) {
    SolveContext context(*this);
    trendoutput::writeGasLineTrendRows(context.trends(), trend, network);
}

void SProd::ImprimeTrendGCab(int trend, int network) {
    SolveContext context(*this);
    trendoutput::writeGasLineTrendHeader(context.trends(), trend, network);
}

double SProd::marchaProdPerm1(double pressureGuess) {
    SolveContext context(*this);
    return steady::marchProductionSteady(context.steady(), pressureGuess);
}

double SProd::buscaProdPfundoPerm(double guess, int attemptCount) {
    SolveContext context(*this);
    return steady::searchProductionBottomHolePressure(context.search(), guess, attemptCount);
}

double SProd::buscaProdPfundoPermRev(double guess) {
    SolveContext context(*this);
    return steady::searchReverseProductionBottomHolePressure(context.search(), guess);
}

double SProd::buscaProdPfundoPerm2(double guess, int attemptCount) {
    SolveContext context(*this);
    return steady::searchProductionBottomHolePressureSecondary(context.search(), guess, attemptCount);
}

double SProd::buscaProdPfundoPerm3(double inletPressure) {
    SolveContext context(*this);
    return steady::searchProductionBottomHolePressureTertiary(context.search(), inletPressure);
}

double SProd::marchaProdPresPres1(double massFlowRateGuess) {
    SolveContext context(*this);
    return steady::marchProductionPressureToPressure(context.steady(), massFlowRateGuess);
}

double SProd::buscaProdPresPresPerm(double massFlowRateGuess, double maximumFlowRate, int iterationCount) {
    SolveContext context(*this);
    return steady::searchProductionPressureToPressure(context.search(), massFlowRateGuess, maximumFlowRate,
                                                      iterationCount);
}

double SProd::buscaProdPresPresPermRev(double massFlowRateGuess, double maximumFlowRate, int iterationCount) {
    SolveContext context(*this);
    return steady::searchReverseProductionPressureToPressure(context.search(), massFlowRateGuess, maximumFlowRate,
                                                             iterationCount);
}

double SProd::buscaProdPresPresPerm2(double massFlowRateGuess, double maximumFlowRate) {
    SolveContext context(*this);
    return steady::searchProductionPressureToPressureSecondary(context.search(), massFlowRateGuess, maximumFlowRate);
}

double SProd::buscaProdPresPresPerm3(double massFlowRateGuess, double maximumFlowRate) {
    SolveContext context(*this);
    return steady::searchProductionPressureToPressureTertiary(context.search(), massFlowRateGuess, maximumFlowRate);
}

double SProd::buscaInjPfundoPerm1(double guess) {
    SolveContext context(*this);
    return steady::searchInjectionBottomHolePressure1(context.search(), guess);
}

double SProd::buscaInjPfundoPerm2(double guess) {
    SolveContext context(*this);
    return steady::searchInjectionBottomHolePressure2(context.search(), guess);
}

double SProd::buscaInjPfundoPerm3(double guess) {
    SolveContext context(*this);
    return steady::searchInjectionBottomHolePressure3(context.search(), guess);
}

double SProd::buscaInjPfundoPerm4() {
    SolveContext context(*this);
    return steady::searchInjectionBottomHolePressure4(context.search());
}

double SProd::buscaInjPfundoPerm5(double guess) {
    SolveContext context(*this);
    return steady::searchInjectionBottomHolePressure5(context.search(), guess);
}

double SProd::hidroreverso(double liquidHoldup, double liquidFlowRate, double gasFlowRate) {
    SolveContext context(*this);
    return steady::reverseHydrostatic(context.steady(), liquidHoldup, liquidFlowRate, gasFlowRate);
}

double SProd::hidroreversoInj(double liquidHoldup, double liquidFlowRate) {
    SolveContext context(*this);
    return steady::reverseInjectionHydrostatic(context.steady(), liquidHoldup, liquidFlowRate);
}

double SProd::hidroTramoSecundario(double quality) {
    SolveContext context(*this);
    return steady::secondaryBranchHydrostatic(context.steady(), quality);
}

void SProd::calcTempFim() {
    SolveContext context(*this);
    thermal::computeOutletTemperature(context.thermal());
}

void SProd::assemble(const CarriedState &carried, double *extraSourceLengths, int *extraSourceCells,
                     int extraSourceCount) {
    tables = PropertyTables();
    steadySearch = SteadySearch();
    allocateSolutionBuffers();
    applyCarriedState(carried);
    montasistema(extraSourceLengths, extraSourceCells, extraSourceCount);
}

void SProd::allocateSolutionBuffers() {
    flut = FullMtx<double>(arq.ncelp, arq.nvarprofp + kProductionProfileFixedColumns);
    flutG = FullMtx<double>(arq.ncelg, arq.nvarprofg + kGasLineProfileFixedColumns);
    transient.systemMatrix =
        BandMtx<double>(kProductionUnknownsPerCell * arq.ncelp, kProductionLeftBandwidth, kProductionRightBandwidth);
    transient.rightHandSide = Vcr<double>(kProductionUnknownsPerCell * arq.ncelp);
    gasLift.systemMatrix = BandMtx<double>(kGasLineUnknownsPerCell * arq.ncelg, kGasLineBandwidth, kGasLineBandwidth);
    gasLift.rightHandSide = Vcr<double>(kGasLineUnknownsPerCell * arq.ncelg);
}

void SProd::applyCarriedState(const CarriedState &carried) {
    noextremo = carried.noextremo;
    noinicial = carried.noinicial;
    derivaAnel = carried.derivaAnel;
    bloq = carried.bloq;
    networkCoupling.previousReverseLiquidFraction = carried.previousReverseLiquidFraction;
    networkCoupling.previousReverseQuality = carried.previousReverseQuality;
    transient.minimumCycleTimeStep = carried.minimumCycleTimeStep;
}

void SProd::HidroDescargaG() {
    SolveContext context(*this);
    gaslift::computeGasUnloadingHydrostatics(context.gasLift());
}

void SProd::HidroDescargaP() {
    SolveContext context(*this);
    gaslift::computeProductionUnloadingHydrostatics(context.gasLift());
}

bool SProd::hydrateEnvelopeDue() const {
    return arq.calculaEnvelope == 1 && vg1dSP->lixo5 <= arq.tfinal;
}

void SProd::solveHydrateEnvelopes() {
    if (!hydrateEnvelopeDue())
        return;
    FA_Hidrato(*this).solverHidrato();
    if (arq.lingas > 0)
        FA_Hidrato_Servico(*this).solverHidratoG();
}

namespace {
// Some views bind a `const T *const &` field to a `T *` member. C++20 binds such
// a reference to the member itself; a temporary there would leave the view
// with a copy that dangles.
constexpr bool bindsToTheMember() {
    int *member = nullptr;
    const int *const &field = member;
    return &field == &member;
}
static_assert(bindsToTheMember());

}  // namespace

sisprod::gaslift::GasLiftState SProd::gasLiftStateOf(sisprod::SolveContext &context) {
    return sisprod::gaslift::GasLiftState{
        .gasCells = gasLift.cells,
        .cells = celula,
        .input = arq,
        .globals = vg1dSP,
        .lastGasCell = gasLift.cellCount,
        .lastCell = ncel,
        .gasLiftChokes = gasLift.valveChokes,
        .injectionChoke = gasLift.injectionChoke,
        .gasValveCellIndices = gasLift.gasValveCells,
        .productionValveCellIndices = gasLift.productionValveCells,
        .gasSystemMatrix = gasLift.systemMatrix,
        .gasRightHandSide = gasLift.rightHandSide,
        .tubingCouplingStart = gasLift.tubingCouplingStart,
        .tubingCouplingEnd = gasLift.tubingCouplingEnd,
        .annulusCouplingStart = gasLift.annulusCouplingStart,
        .annulusCouplingEnd = gasLift.annulusCouplingEnd,
        .annulusThermalCoupling = gasLift.annulusThermalCoupling,
        .thermalSourceDisabled = thermalSourceDisabled,
        .injectionPressure = gasLift.injectionPressure,
        .injectionTemperature = gasLift.injectionTemperature,
        .separatorPressure = gasLift.separatorPressure,
        .timeStep = dt,
        .interfaceCell = gasLift.interfaceCell,
        .interfaceVelocity = gasLift.interfaceVelocity,
        .interfaceTimeStep = gasLift.interfaceTimeStep,
        .previousInterfaceCell = gasLift.previousInterfaceCell,
        .previousInterfaceVelocity = gasLift.previousInterfaceVelocity,
        .previousInterfaceTimeStep = gasLift.previousInterfaceTimeStep,
        .accumulatedUnloadingVolume = gasLift.accumulatedUnloadingVolume,
        .unloadingAveragingInterval = gasLift.unloadingAveragingInterval,
        .maximumUnloadingFlowSamples = gasLift.maximumUnloadingFlowSamples,
        .unloadingTimeStepSamples = gasLift.unloadingTimeStepSamples,
        .unloadingAveragingWindow = gasLift.unloadingAveragingWindow,
        .unloadingSampleCapacity = gasLift.unloadingSampleCapacity,
        .temperatureUpdater = {context},
    };
}

sisprod::steady::SteadyStateState SProd::steadyStateOf(sisprod::SolveContext &context) {
    return sisprod::steady::SteadyStateState{
        .cells = celula,
        .gasCells = gasLift.cells,
        .input = arq,
        .globals = vg1dSP,
        .lastCell = ncel,
        .lastGasCell = gasLift.cellCount,
        .injectionChoke = gasLift.injectionChoke,
        .surfaceChoke = chokeSup,
        .gasValveCellIndices = gasLift.gasValveCells,
        .productionValveCellIndices = gasLift.productionValveCells,
        .steadyIteration = steadyIteration,
        .searchOrigin = steadySearch.initialEstimateSearch,
        .convergenceMonitor = steadySearch.convergenceMonitor,
        .baseConvergenceMonitor = steadySearch.convergenceMonitorReference,
        .annulusDrift = derivaAnel,
        .networkCoupled = gasLift.annulusThermalCoupling,
        .endNode = noextremo,
        .thermalSourceDisabled = thermalSourceDisabled,
        .slowHeatTransferThreshold = steadySearch.slowHeatTransferThreshold,
        .gasSurfacePressure = gasLift.separatorPressure,
        .initialGasPressure = gasLift.injectionPressure,
        .initialGasTemperature = gasLift.injectionTemperature,
        .outletPressure = presfim,
        .timeStep = dt,
        .casingTemperature = steadySearch.reverseFlowTemperature,
        .inletQuality = titE,
        .productionFluidCount = steadySearch.productionFluidCount,
        .updaters = {context},
    };
}

/// The state a boundary-condition search reads.
///
/// It composes the march state rather than rebuilding it: the two halves read
/// many of the same SProd members, and two bindings of one member could
/// disagree silently.
sisprod::steady::SteadyStateSearchState SProd::searchStateOf(sisprod::SolveContext &context) {
    return sisprod::steady::SteadyStateSearchState{
        .march = steadyStateOf(context),
        .holdupGuess = steadySearch.holdupGuess,
        .reverseNetworkFluid = steadySearch.reverseNetworkFluid,
        .reverseSteady = steadySearch.steadyReverseFlow,
    };
}

/// The state one transient step reads.
sisprod::transient::TransientStepState SProd::transientStateOf(sisprod::SolveContext &context) {
    return sisprod::transient::TransientStepState{
        .meanMaximumTemperatureRate = transient.meanMaximumTemperatureRate,
        .meanMaximumPressureRate = transient.meanMaximumPressureRate,
        .masterState = transient.masterValveState,
        .surfaceChokeOpen = transient.surfaceChokeOpen,
        .initiallyOpen = transient.previousSurfaceChokeOpen,
        .timeChanged = transient.adjustTimeStep,
        .artificialGasThisStep = transient.artificialGasThisStep,
        .inletComplementaryFraction = betaE,
        .interfaceCell = gasLift.interfaceCell,
        .masterCounter = transient.masterValveCounter,
        .gasSpecificHeatTable = tables.gasSpecificHeatTable,
        .timeStep = dt,
        .meanCflTimeStep = transient.meanCflTimeStep,
        .totalCflTimeStep = transient.totalCflTimeStep,
        .interfaceTimeStep = gasLift.interfaceTimeStep,
        .meanSimulationTimeStep = transient.meanAcceptedTimeStep,
        .totalSimulationTimeStep = transient.totalAcceptedTimeStep,
        .bufferedComplementaryMassSource = transient.bufferedComplementaryMassSource,
        .bufferedGasMassSource = transient.bufferedGasMassSource,
        .bufferedLiquidMassSource = transient.bufferedLiquidMassSource,
        .eventIndex = transient.eventIndex,
        .slugCount = transient.slugCount,
        .timeStepRestrictionCount = transient.restrictedStepsLeft,
        .surfaceChokeMassFlag = masChkSup,
        .fullModel = transient.fullModel,
        .desperationMoment = transient.artificialGasMultiplier,
        .chokeModeChanged = transient.chokeModeChanged,
        .multiplier = transient.cflSafetyFactor,
        .outletPressure = presfim,
        .restart = transient.rollbackRequested,
        .timeStepRestricted = transient.timeStepRestricted,
        .openTime = transient.activeChokeDelay,
        .productionSolution = transient.rightHandSide,
        .inletQuality = titE,
        .masterRatio0 = transient.masterValveRatiosInactive.data(),
        .masterRatio1 = transient.masterValveRatiosActive.data(),
        .masterCriticalRatio = transient.masterValveCriticalRatios.data(),
        .interfaceVelocity = gasLift.interfaceVelocity,
        .movingPigCount = transient.movingPigCount,
        .previousMovingPigCount = transient.previousMovingPigCount,
        .masterOpenSchedule = transient.masterValveOpenTimes,
        .input = arq,
        .initialInterfaceCell = gasLift.previousInterfaceCell,
        .cells = celula,
        .gasCells = gasLift.cells,
        .surfaceChoke = chokeSup,
        .cflTimeSteps = transient.cflTimeSteps,
        .initialInterfaceTimeStep = gasLift.previousInterfaceTimeStep,
        .simulationTimeSteps = transient.acceptedTimeSteps,
        .auxiliaryCflTimeStep = transient.auxiliaryCflTimeStep,
        .finalAuxiliaryTimeStep = transient.auxiliaryAcceptedTimeStep,
        .masterCloseSchedule = transient.masterValveCloseTimes,
        .productionFreeTerms = flut,
        .gasFreeTerms = flutG,
        .branchIndex = indTramo,
        .movingMeanFlux = transient.movingMeanFlux,
        .stepIndex = transient.stepCount,
        .productionMatrix = transient.systemMatrix,
        .smallestCellLength = transient.smallestCellLength,
        .masterOpenCount = static_cast<int>(transient.masterValveOpenTimes.size()),
        .lastCell = ncel,
        .lastGasCell = gasLift.cellCount,
        .productionCrossSectionCount = trends.productionWallProfileCells,
        .masterCloseCount = static_cast<int>(transient.masterValveCloseTimes.size()),
        .endNode = noextremo,
        .gasSurfacePressure = gasLift.separatorPressure,
        .inletPressure = presE,
        .movingMeanPressure = transient.movingMeanPressure,
        .movingMeanStartTime = transient.movingMeanStartTime,
        .maximumTemperatureRates = transient.maximumTemperatureRates,
        .maximumPressureRates = transient.maximumPressureRates,
        .inletTemperature = tempE,
        .reverseQuality = networkCoupling.reverseQuality,
        .initialInterfaceVelocity = gasLift.previousInterfaceVelocity,
        .pigReceiverCells = transient.pigReceiverCells,
        .scheduledPigCount = static_cast<int>(transient.pigReceiverCells.size()),
        .globals = vg1dSP,
        .updaters = {context},
    };
}

/// The fluid tracking module's view of SProd: eleven members by reference.
sisprod::fluidtracking::FluidTrackingState SProd::fluidTrackingStateOf(sisprod::SolveContext &context) {
    return sisprod::fluidtracking::FluidTrackingState{
        .cells = celula,
        .lastCell = ncel,
        .input = arq,
        .globals = vg1dSP,
        .endNode = noextremo,
        .inletPressure = presE,
        .inletTemperature = tempE,
        .inletComplementaryFraction = betaE,
        .trackGasOilRatio = transient.trackGasOilRatio,
        .trackGasGravity = transient.trackGasGravity,
        .compositionalRefreshCounter = transient.compositionalRefreshCounter,
        .updaters = {context},
    };
}


/// The state SolveTrans reads: the step state, composed rather than rebuilt,
/// plus the members only the solve touches.
sisprod::transient::TransientSolveState SProd::transientSolveStateOf(sisprod::SolveContext &context) {
    return sisprod::transient::TransientSolveState{
        .step = transientStateOf(context),
        .defaultInletTemperature = temperatura,
        .annulusDrift = derivaAnel,
        .closingSubtitles = saidaSubTextoSis,
        .closingTitles = saidaTextoSis,
        .compositionalRefreshCounter = transient.compositionalRefreshCounter,
        .fluxHistory = transient.fluxSamples,
        .gasCrossSectionCellCounts = trends.gasWallProfileCells,
        .gasCrossSectionProfileTimeCounter = trends.gasWallProfileTimeIndex,
        .gasCrossSectionTrendBufferedCounts = trends.gasWallTrendBufferedCounts,
        .gasCrossSectionTrendCounts = trends.gasWallTrendCounts,
        .gasCrossSectionTrendMatrix = trends.gasWallTrendMatrices,
        .gasCrossSectionTrendResetTimers = trends.gasWallTrendResetTimes,
        .gasProfileTimeCounter = trends.gasProfileTimeIndex,
        .gasTrendBufferedCounts = trends.gasTrendBufferedCounts,
        .gasTrendCounts = trends.gasTrendCounts,
        .gasTrendMatrix = trends.gasTrendMatrices,
        .gasTrendResetTimers = trends.gasTrendResetTimes,
        .initialGasPressure = gasLift.injectionPressure,
        .initialGasSurfacePressure = transient.previousSeparatorPressure,
        .initialGasTemperature = gasLift.injectionTemperature,
        .initialOpenTime = transient.previousActiveChokeDelay,
        .injectionChoke = gasLift.injectionChoke,
        .logBuffer = tmpLog,
        .logCounter = transient.eventLogCount,
        .massTransferModel = transient.massTransferModel,
        .minimumCycleTimeStep = transient.minimumCycleTimeStep,
        .movingMeanCounter = transient.movingMeanWindow,
        .movingMeanVoidFraction = transient.movingMeanVoidFraction,
        .networkCoupled = gasLift.annulusThermalCoupling,
        .poissonSolver3D = coupling3D.poissonSolver3D,
        .pressureHistory = transient.pressureSamples,
        .printCounter = transient.printCounter,
        .printPassCount = trends.trendOutputCount,
        .productionCrossSectionProfileTimeCounter = trends.productionWallProfileTimeIndex,
        .productionCrossSectionTrendBufferedCounts = trends.productionWallTrendBufferedCounts,
        .productionCrossSectionTrendCounts = trends.productionWallTrendCounts,
        .productionCrossSectionTrendMatrix = trends.productionWallTrendMatrices,
        .productionCrossSectionTrendResetTimers = trends.productionWallTrendResetTimes,
        .productionProfileTimeCounter = trends.productionProfileTimeIndex,
        .productionTrendBufferedCounts = trends.productionTrendBufferedCounts,
        .productionTrendCounts = trends.productionTrendCounts,
        .productionTrendMatrix = trends.productionTrendMatrices,
        .productionTrendResetTimers = trends.productionTrendResetTimes,
        .startNode = noinicial,
        .timeStepHistory = transient.timeSamples,
        .totalFlux = transient.accumulatedFlux,
        .totalPressure = transient.accumulatedPressure,
        .totalVoidFraction = transient.accumulatedVoidFraction,
        .trackGasGravity = transient.trackGasGravity,
        .trackGasOilRatio = transient.trackGasOilRatio,
        .unitCellTimeCounters = trends.singleCellOutputTimeIndices,
        .voidFractionHistory = transient.voidFractionSamples,
        .updaters = {context},
    };
}



/// The state the five drift-flux closure variants read.
driftflux::coefficient::ClosureState SProd::closureStateOf() {
    return driftflux::coefficient::ClosureState{
        .cells = celula,
        .lastCell = ncel,
        .globals = vg1dSP,
        .input = arq,
        .selectors = {arq.CorreDisper, arq.CorreAnular, arq.CorreEstrat},
        .gasSurfaceTemperature = tGSup,
        .inletVoidFraction = alfE,
        .inletComplementaryFraction = betaE,
        .inletPressure = presE,
        .inletTemperature = tempE,
        .steadyIteration = steadyIteration,
    };
}

/// The state the source terms of a cell read and write.
sisprod::sources::SourceState SProd::sourceStateOf() {
    return sisprod::sources::SourceState{
        .cells = celula,
        .input = arq,
        .globals = vg1dSP,
        .lastCell = ncel,
        .steadyMode = modoPerm,
        .parallelSecondaryBranch = networkCoupling.parallelSecondaryBranch,
        .parallelSecondaryBoundaryCondition = networkCoupling.parallelSecondaryBoundaryCondition,
        .parallelSourceCells = networkCoupling.parallelSourceCells,
        .parallelSourceProductionLiquid = networkCoupling.parallelSourceProductionLiquid,
        .parallelSourceComplementaryLiquid = networkCoupling.parallelSourceComplementaryLiquid,
        .parallelSourceGas = networkCoupling.parallelSourceGas,
    };
}

sisprod::thermal::ThermalState SProd::thermalStateOf(sisprod::SolveContext &context) {
    return sisprod::thermal::ThermalState{
        .cells = celula,
        .gasCells = gasLift.cells,
        .input = arq,
        .latentHeatTable = tables.latentHeatTable.rowPointers,
        .globals = vg1dSP,
        .thermalSourceDisabled = thermalSourceDisabled,
        .productionNetworkCoupled = networkCoupling.secondaryThermalCoupling,
        .primarySectionStart = networkCoupling.secondaryCouplingStart,
        .primarySectionEnd = networkCoupling.secondaryCouplingEnd,
        .coupledCellIndices = coupling3D.coupledIndexCorrections,
        .poissonSolver3D = coupling3D.poissonSolver3D,
        .lastCell = ncel,
        .surfaceChoke = chokeSup,
        .surfaceChokeMassFlag = masChkSup,
        .endNode = noextremo,
        .gasSurfaceTemperature = tGSup,
        .latentHeatEnabled = CalcLat,
        .sourceUpdater = {.context = context},
        .fullModel = transient.fullModel,
        .massTransferModel = transient.massTransferModel,
        .closureUpdater = {.context = context},
        .inletPressure = presE,
        .inletTemperature = tempE,
        .inletQuality = titE,
        .inletVoidFraction = alfE,
        .inletComplementaryFraction = betaE,
        .evolutionUpdater = {.context = context},
        .surfaceChokeOpen = transient.surfaceChokeOpen,
        .defaultInletTemperature = temperatura,
        .timeStep = dt,
        .minimumCycleTimeStep = transient.minimumCycleTimeStep,
        .poisson2DCellIndices = transient.poisson2DCells,
        .poisson2DCellCount = static_cast<int>(transient.poisson2DCells.size()),
        .steadyIteration = steadyIteration,
        .slowHeatTransferThreshold = steadySearch.slowHeatTransferThreshold,
        .annulusTubingStart = gasLift.tubingCouplingStart,
        .annulusTubingEnd = gasLift.tubingCouplingEnd,
        .tubingAnnulusStart = gasLift.annulusCouplingStart,
        .productionNetworkHeatCoupled = networkCoupling.primaryThermalCoupling,
        .primaryNetworkSectionEnd = networkCoupling.primaryCouplingEnd,
        .primaryNetworkSectionStart = networkCoupling.primaryCouplingStart,
        .lastGasCell = gasLift.cellCount,
        .gasLiftChokes = gasLift.valveChokes,
        .gasSurfacePressure = gasLift.separatorPressure,
        .outletPressure = presfim,
        .surfaceTemperature = networkCoupling.chokeOutletTemperature,
        .networkCoupled = gasLift.annulusThermalCoupling,
        .tubingAnnulusEnd = gasLift.annulusCouplingEnd,
    };
}

/// Binds the state the trend writers are allowed to read.
///
/// Designated initialisers, not positional: the four series have identical
/// types, so a positional swap would compile in silence and hand a writer
/// another line's buffer. Naming each one makes that a compile error.
///
/// Every field is a reference, never a copy -- the caller advances the counters
/// between the header call and the row call, so a copy would be read at the
/// wrong moment.
trendoutput::TrendState SProd::trendStateOf() const {
    return trendoutput::TrendState{
        .input = arq,
        .globals = vg1dSP,
        .branchIndex = indTramo,
        .printPassCount = trends.trendOutputCount,
        .production = {trends.productionTrendMatrices, trends.productionTrendCounts,
                       trends.productionTrendBufferedCounts},
        .gasLine = {trends.gasTrendMatrices, trends.gasTrendCounts, trends.gasTrendBufferedCounts},
        .productionCrossSection = {trends.productionWallTrendMatrices, trends.productionWallTrendCounts,
                                   trends.productionWallTrendBufferedCounts},
        .gasLineCrossSection = {trends.gasWallTrendMatrices, trends.gasWallTrendCounts,
                                trends.gasWallTrendBufferedCounts}};
}

namespace sisprod {

SolveContext::SolveContext(SProd &system)
    : system_(system),
      gasLift_(system.gasLiftStateOf(*this)),
      steady_(system.steadyStateOf(*this)),
      search_(system.searchStateOf(*this)),
      transientStep_(system.transientStateOf(*this)),
      fluidTracking_(system.fluidTrackingStateOf(*this)),
      transientSolve_(system.transientSolveStateOf(*this)),
      closure_(system.closureStateOf()),
      sources_(system.sourceStateOf()),
      thermal_(system.thermalStateOf(*this)),
      trends_(system.trendStateOf()) {}

}  // namespace sisprod
