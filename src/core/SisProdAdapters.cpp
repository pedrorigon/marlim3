/// The connection between SProd and the modules: the state view each module
/// receives (sisprod::adapters) and the callbacks through which a module asks
/// SProd for what it still computes itself.
#include "SisProd.h"
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
    system.tempDescarga(cellIndex);
}

void sisprod::gaslift::GasLiftTemperatureUpdater::computeGasTemperature(
    int cellIndex, double previousTemperature, int steadyMode) const {
    system.calctempGas(cellIndex, previousTemperature, steadyMode);
}

double sisprod::gaslift::GasLiftTemperatureUpdater::gasLiftDischargeTemperature(
    int valveIndex) const {
    return system.TempDescGL(valveIndex);
}

// ------------------------------------------- steady-state callbacks ----
//
// Sixteen forwards, fourteen of which reach code in another module. They come
// back through SProd because reaching sisprod::gaslift or sisprod::thermal
// needs one of their state structs, and the adapters in this file are the
// only place that builds those.
void sisprod::steady::SteadyStateUpdaters::steadyDriftClosure(int cellIndex, double &c0, double &ud) const {
    system.CalcC0UdPerm(cellIndex, c0, ud);
}
void sisprod::steady::SteadyStateUpdaters::updateSource(int cellIndex) const {
    system.renovaFonte(cellIndex);
}
void sisprod::steady::SteadyStateUpdaters::advanceSteadyTemperature(int cellIndex, int rungeKuttaStage) const {
    system.RenovaTempPerm(cellIndex, rungeKuttaStage);
}
void sisprod::steady::SteadyStateUpdaters::advanceReverseSteadyTemperature(int cellIndex, int rungeKuttaStage) const {
    system.RenovaTempPermRev(cellIndex, rungeKuttaStage);
}
void sisprod::steady::SteadyStateUpdaters::computeTemperature(int cellIndex, double previousTemperature,
                                                              int steadyMode) const {
    system.calctemp(cellIndex, previousTemperature, steadyMode);
}
void sisprod::steady::SteadyStateUpdaters::computeGasTemperature(int cellIndex, double previousTemperature,
                                                                 int steadyMode) const {
    system.calctempGas(cellIndex, previousTemperature, steadyMode);
}
void sisprod::steady::SteadyStateUpdaters::updateProductionTemperaturePeriphery(int cellIndex) const {
    system.atualizaPeriTempProd(cellIndex);
}
void sisprod::steady::SteadyStateUpdaters::initializeSteadyValveGasFlowRate(int cellIndex) const {
    system.IniciaVazValvGasPerm(cellIndex);
}
void sisprod::steady::SteadyStateUpdaters::initializeTubingConnectionSteady() const {
    system.IniciaconectaColunaPerm();
}
void sisprod::steady::SteadyStateUpdaters::updateSteadyGasPressure(int cellIndex) const {
    system.RenovaPresGasPerm(cellIndex);
}
void sisprod::steady::SteadyStateUpdaters::updateSteadyGasTemperature(int cellIndex) const {
    system.RenovaTempGasPerm(cellIndex);
}
void sisprod::steady::SteadyStateUpdaters::computeSteadyGasFlowRate(int cellIndex) const {
    system.calcVazGasPerm(cellIndex);
}
void sisprod::steady::SteadyStateUpdaters::connectTubing() const {
    system.conectaColuna();
}
void sisprod::steady::SteadyStateUpdaters::connectTubingSteady() const {
    system.conectaColunaPerm();
}
double sisprod::steady::SteadyStateUpdaters::steadyGasPressureDrop(int cellIndex) const {
    return system.delpGasPerm(cellIndex);
}
double sisprod::steady::SteadyStateUpdaters::steadyInjectionPressureDrop(int cellIndex) const {
    return system.delpInjPerm(cellIndex);
}
// The two searches a march calls; their results are discarded at the call
// site.
void sisprod::steady::SteadyStateUpdaters::searchGasPressureSteadySecondary() const {
    system.buscaGasPresPerm2();
}
void sisprod::steady::SteadyStateUpdaters::searchGasPressureSteadyTertiary() const {
    system.buscaGasPresPerm3();
}

void sisprod::thermal::ThermalSourceUpdater::operator()(int cellIndex) const {
    system.renovaFonte(cellIndex);
}

void sisprod::thermal::ThermalClosureUpdater::instantaneous(
    int cellIndex, double &distribution, double &driftVelocity) const {
    system.CalcC0Ud(cellIndex, distribution, driftVelocity);
}

void sisprod::thermal::ThermalClosureUpdater::buffered(
    int cellIndex, double &distribution, double &driftVelocity) const {
    system.CalcC0UdBuf(cellIndex, distribution, driftVelocity);
}

void sisprod::thermal::ThermalClosureUpdater::initialization(
    int cellIndex, double &distribution, double &driftVelocity) const {
    system.CalcC0UdIni(cellIndex, distribution, driftVelocity);
}

void sisprod::thermal::ThermalClosureUpdater::bufferedInitialization(
    int cellIndex, double &distribution, double &driftVelocity) const {
    system.CalcC0UdIniBuf(cellIndex, distribution, driftVelocity);
}

void sisprod::thermal::ThermalEvolutionUpdater::solvePressureVelocityCoupling(
    int cycle) const {
    system.SolveAcopPV(cycle);
}

void sisprod::thermal::ThermalEvolutionUpdater::renew() const {
    system.renova();
}

void sisprod::transient::TransientStepUpdaters::advanceGasSubStep() const {
    system.subtempoGas();
}

void sisprod::composition::CompositionUpdaters::correctGasSpecificGravity(int i) const {
    system.corrDeng(i);
}

void sisprod::transient::TransientSolveUpdaters::solveHydrateEnvelopes() const {
    system.solveHydrateEnvelopes();
}
double sisprod::transient::TransientSolveUpdaters::searchUnloadingInjectionPressure() const {
    return system.BuscaPresInjDesc();
}
void sisprod::transient::TransientSolveUpdaters::writeProductionTrendHeader(int i, int nrede) const {
    system.ImprimeTrendPCab(i, nrede);
}
void sisprod::transient::TransientSolveUpdaters::writeProductionTrendRows(int i, int nrede) const {
    system.ImprimeTrendP(i, nrede);
}
void sisprod::transient::TransientSolveUpdaters::writeGasTrendHeader(int i, int nrede) const {
    system.ImprimeTrendGCab(i, nrede);
}
void sisprod::transient::TransientSolveUpdaters::writeGasTrendRows(int i, int nrede) const {
    system.ImprimeTrendG(i, nrede);
}
void sisprod::transient::TransientSolveUpdaters::writeProductionCrossSectionTrendHeader(int i) const {
    system.ImprimeTrendTransPCab(i);
}
void sisprod::transient::TransientSolveUpdaters::writeProductionCrossSectionTrendRows(int i) const {
    system.ImprimeTrendTransP(i);
}
void sisprod::transient::TransientSolveUpdaters::writeGasCrossSectionTrendHeader(int i) const {
    system.ImprimeTrendTransGCab(i);
}
void sisprod::transient::TransientSolveUpdaters::writeGasCrossSectionTrendRows(int i) const {
    system.ImprimeTrendTransG(i);
}
void sisprod::transient::TransientSolveUpdaters::evaluateParaffin() const {
    system.avaliaParafina();
}
void sisprod::transient::TransientSolveUpdaters::connectTubing() const {
    system.conectaColuna();
}
void sisprod::transient::TransientSolveUpdaters::marchTransientEnergy(int ciclo, int ciclomax) const {
    system.marchaEnergTrans(ciclo, ciclomax);
}
void sisprod::transient::TransientSolveUpdaters::updateMolarFractions(const ProFlu &fluiRev) const {
    system.renovaFracMol2(fluiRev);
}
void sisprod::transient::TransientSolveUpdaters::updateDensities() const {
    system.renovaMasEsp();
}
void sisprod::transient::TransientSolveUpdaters::updateGasOilRatioAndCo2(const ProFlu &fluiRev) const {
    system.renovaRGOdgYco2(fluiRev);
}
void sisprod::transient::TransientSolveUpdaters::updateTemperatures() const {
    system.renovaTemp();
}
void sisprod::transient::TransientSolveUpdaters::updateInitialFractions() const {
    system.renovaalbetini();
}
void sisprod::transient::TransientSolveUpdaters::updateThermal(int aflu) const {
    system.renovaterm(aflu);
}
void sisprod::transient::TransientSolveUpdaters::solveGasLine() const {
    system.solveLinGas();
}



namespace sisprod::adapters {

sisprod::gaslift::GasLiftState gasLiftStateOf(SProd &system) {
    return sisprod::gaslift::GasLiftState{
        .gasCells = system.celulaG,
        .cells = system.celula,
        .input = system.arq,
        .globals = system.vg1dSP,
        .gasCellCount = system.ncelGas,
        .lastCell = system.ncel,
        .gasLiftChokes = system.gasLift.chokeVGL,
        .injectionChoke = system.chokeInj,
        .gasValveCellIndices = system.gasLift.posicVGLG,
        .productionValveCellIndices = system.gasLift.posicVGLP,
        .gasSystemMatrix = system.gasLift.matglobG,
        .gasFreeTerms = system.gasLift.termolivreG,
        .annulusTubingStart = system.gasLift.ColunaAnulaIni,
        .annulusTubingEnd = system.gasLift.ColunaAnulaFim,
        .tubingAnnulusStart = system.gasLift.AnulaColunaIni,
        .tubingAnnulusEnd = system.gasLift.AnulaColunaFim,
        .steadyIteration = system.iterperm,
        .networkCoupled = system.gasLift.verificaAcop,
        .thermalSourceDisabled = system.semTermo,
        .initialGasPressure = system.presiniG,
        .initialGasTemperature = system.tempiniG,
        .gasSurfacePressure = system.pGSup,
        .timeStep = system.dt,
        .interfaceCell = system.celInter,
        .interfaceVelocity = system.velInter,
        .interfaceTimeStep = system.dtInter,
        .initialInterfaceCell = system.celInterIni,
        .initialInterfaceVelocity = system.velInterIni,
        .initialInterfaceTimeStep = system.dtInterIni,
        .meanUnloadingFlowRate = system.gasLift.vazmedDesc,
        .meanUnloadingTemperature = system.gasLift.tempmedDEsc,
        .maximumMeanUnloadingFlowRates = system.gasLift.vazmaxMedDesc,
        .unloadingTimeSteps = system.gasLift.dtDesc,
        .continuousMeanUnloadingTemperature = system.gasLift.tempMedContDesc,
        .maximumContinuousUnloadingCount = system.gasLift.maxVecContDesc,
        .temperatureUpdater = {system},
    };
}

sisprod::steady::SteadyStateState steadyStateOf(SProd &system) {
    return sisprod::steady::SteadyStateState{
        .cells = system.celula,
        .gasCells = system.celulaG,
        .input = system.arq,
        .globals = system.vg1dSP,
        .lastCell = system.ncel,
        .gasCellCount = system.ncelGas,
        .injectionChoke = system.chokeInj,
        .surfaceChoke = system.chokeSup,
        .gasValveCellIndices = system.gasLift.posicVGLG,
        .productionValveCellIndices = system.gasLift.posicVGLP,
        .steadyIteration = system.iterperm,
        .searchOrigin = system.buscaIni,
        .convergenceMonitor = system.monitConvPerm,
        .baseConvergenceMonitor = system.monitConvPermBase,
        .annulusDrift = system.derivaAnel,
        .networkCoupled = system.gasLift.verificaAcop,
        .endNode = system.noextremo,
        .thermalSourceDisabled = system.semTermo,
        .slowHeatTransferThreshold = system.trocaTermicaLenta,
        .gasSurfacePressure = system.pGSup,
        .initialGasPressure = system.presiniG,
        .initialGasTemperature = system.tempiniG,
        .outletPressure = system.presfim,
        .timeStep = system.dt,
        .defaultInletTemperature = system.temperatura,
        .casingTemperature = system.tempRev,
        .inletQuality = system.titE,
        .productionFluidCount = system.nfluP,
        .updaters = {system},
    };
}

/// The state a boundary-condition search reads.
///
/// It composes the march state rather than rebuilding it: the two halves read
/// many of the same SProd members, and two bindings of one member could
/// disagree silently.
sisprod::steady::SteadyStateSearchState searchStateOf(SProd &system) {
    return sisprod::steady::SteadyStateSearchState{
        .march = steadyStateOf(system),
        .holdupGuess = system.chuteHol,
        .reverseNetworkFluid = system.fluiRevRede,
        .reverseSteady = system.revPerm,
    };
}

/// The state one transient step reads.
sisprod::transient::TransientStepState transientStateOf(SProd &system) {
    return sisprod::transient::TransientStepState{
        .meanMaximumTimeStep = system.DTMaxMed,
        .meanMaximumPressureChange = system.DpMaxMed,
        .masterState = system.EstadoMaster1,
        .surfaceChokeOpen = system.aberto,
        .initiallyOpen = system.abertoini,
        .timeChanged = system.alteraTempo,
        .inletCompletionFraction = system.betaE,
        .interfaceCell = system.celInter,
        .masterCounter = system.contaMaster1,
        .gasSpecificHeatTable = system.tables.cpg,
        .timeStep = system.dt,
        .meanCflTimeStep = system.dtCFLMed,
        .totalCflTimeStep = system.dtCFLTotal,
        .interfaceTimeStep = system.dtInter,
        .meanSimulationTimeStep = system.dtSimMed,
        .totalSimulationTimeStep = system.dtSimTotal,
        .bufferedCompletionMassSource = system.fontemassCRBuf,
        .bufferedGasMassSource = system.fontemassGRBuf,
        .bufferedLiquidMassSource = system.fontemassPRBuf,
        .eventIndex = system.indevento,
        .slugCount = system.kontaGolfada,
        .timeStepRestrictionCount = system.kontarestriDt,
        .surfaceChokeMassFlag = system.masChkSup,
        .fullModel = system.modeloCompleto,
        .desperationMoment = system.momentoDesesp,
        .chokeModeChanged = system.mudaModoChk,
        .multiplier = system.mult,
        .outletPressure = system.presfim,
        .restart = system.reinicia,
        .timeStepRestricted = system.restriDt,
        .openTime = system.tempoaberto,
        .productionSolution = system.termolivreP,
        .inletQuality = system.titE,
        .masterRatio0 = system.vRazMast0,
        .masterRatio1 = system.vRazMast1,
        .masterCriticalRatio = system.vRazMastCrit,
        .interfaceVelocity = system.velInter,
        .masterOpenSchedule = system.abreM1,
        .input = system.arq,
        .initialInterfaceCell = system.celInterIni,
        .cells = system.celula,
        .gasCells = system.celulaG,
        .surfaceChoke = system.chokeSup,
        .cflTimeSteps = system.dtCFL,
        .initialInterfaceTimeStep = system.dtInterIni,
        .simulationTimeSteps = system.dtSim,
        .auxiliaryCflTimeStep = system.dtauxCFL,
        .finalAuxiliaryTimeStep = system.dtauxFinal,
        .masterCloseSchedule = system.fechaM1,
        .productionFreeTerms = system.flut,
        .gasFreeTerms = system.flutG,
        .branchIndex = system.indTramo,
        .movingMeanFlux = system.jMedMov,
        .stepIndex = system.kSP,
        .productionMatrix = system.matglobP,
        .smallestCellLength = system.menorDx,
        .masterOpenCount = system.nabreM1,
        .lastCell = system.ncel,
        .gasCellCount = system.ncelGas,
        .productionCrossSectionCount = system.ncelperftransp,
        .masterCloseCount = system.nfechaM1,
        .endNode = system.noextremo,
        .gasSurfacePressure = system.pGSup,
        .inletPressure = system.presE,
        .movingMeanPressure = system.presMedMov,
        .movingMeanTemperature = system.tMedMov,
        .maximumTimeStepRates = system.taxaDTMax,
        .maximumPressureRates = system.taxaDpMax,
        .inletTemperature = system.tempE,
        .reverseQuality = system.titRev,
        .initialInterfaceVelocity = system.velInterIni,
        .globals = system.vg1dSP,
        .updaters = {system},
    };
}

/// The composition module's view of SProd: seventeen members by reference.
sisprod::composition::CompositionState compositionStateOf(SProd &system) {
    return sisprod::composition::CompositionState{
        .cells = system.celula,
        .lastCell = system.ncel,
        .input = system.arq,
        .globals = system.vg1dSP,
        .endNode = system.noextremo,
        .inletPressure = system.presE,
        .inletTemperature = system.tempE,
        .inletQuality = system.titE,
        .inletCompletionFraction = system.betaE,
        .trackGasOilRatio = system.trackRGO,
        .trackGasGravity = system.trackDeng,
        .compositionalRefreshCounter = system.kontaRenovaComp,
        .surfaceChokeMassFlag = system.masChkSup,
        .movingPigCount = system.indpigP,
        .previousMovingPigCount = system.indpigPini,
        .scheduledPigCount = system.npig,
        .pigReceiverCells = system.receb,
        .updaters = {system},
    };
}


/// The state SolveTrans reads: the step state, composed rather than rebuilt,
/// plus the members only the solve touches.
sisprod::transient::TransientSolveState transientSolveStateOf(SProd &system) {
    return sisprod::transient::TransientSolveState{
        .step = transientStateOf(system),
        .defaultInletTemperature = system.temperatura,
        .annulusDrift = system.derivaAnel,
        .closingSubtitles = system.saidaSubTextoSis,
        .closingTitles = system.saidaTextoSis,
        .compositionalRefreshCounter = system.kontaRenovaComp,
        .fluxHistory = system.jVet,
        .gasCrossSectionCellCounts = system.ncelperftransg,
        .gasCrossSectionProfileTimeCounter = system.kontaTempoTransProfG,
        .gasCrossSectionTrendBufferedCounts = system.ntrendtransgB,
        .gasCrossSectionTrendCounts = system.ntrendtransg,
        .gasCrossSectionTrendMatrix = system.MatTrendTransG,
        .gasCrossSectionTrendResetTimers = system.resettrendtransg,
        .gasProfileTimeCounter = system.kontaTempoProfG,
        .gasTrendBufferedCounts = system.ntrendgB,
        .gasTrendCounts = system.ntrendg,
        .gasTrendMatrix = system.MatTrendG,
        .gasTrendResetTimers = system.resettrendg,
        .initialGasPressure = system.presiniG,
        .initialGasSurfacePressure = system.pGSupIni,
        .initialGasTemperature = system.tempiniG,
        .initialOpenTime = system.tempoabertoini,
        .injectionChoke = system.chokeInj,
        .logBuffer = system.tmpLog,
        .logCounter = system.contaLog,
        .massTransferModel = system.TransMassModel,
        .minimumCycleTimeStep = system.dtCicMin,
        .movingMeanCounter = system.ktMedMov,
        .movingMeanVoidFraction = system.alfMedMov,
        .networkCoupled = system.gasLift.verificaAcop,
        .poissonSolver3D = system.poisson3D,
        .pressureHistory = system.presVet,
        .printCounter = system.KontaImprime,
        .printPassCount = system.kimpT,
        .productionCrossSectionProfileTimeCounter = system.kontaTempoTransProf,
        .productionCrossSectionTrendBufferedCounts = system.ntrendtransB,
        .productionCrossSectionTrendCounts = system.ntrendtrans,
        .productionCrossSectionTrendMatrix = system.MatTrendTransP,
        .productionCrossSectionTrendResetTimers = system.resettrendtrans,
        .productionProfileTimeCounter = system.kontaTempoProf,
        .productionTrendBufferedCounts = system.ntrendB,
        .productionTrendCounts = system.ntrend,
        .productionTrendMatrix = system.MatTrendP,
        .productionTrendResetTimers = system.resettrend,
        .startNode = system.noinicial,
        .temperatureHistory = system.tVet,
        .totalFlux = system.jTotal,
        .totalPressure = system.pTotal,
        .totalVoidFraction = system.alfTotal,
        .trackGasGravity = system.trackDeng,
        .trackGasOilRatio = system.trackRGO,
        .unitCellTimeCounters = system.kontaTempoCelUni,
        .voidFractionHistory = system.alfVet,
        .updaters = {system},
    };
}



/// The state the five drift-flux closure variants read.
driftflux::coefficient::ClosureState closureStateOf(SProd &system) {
    return driftflux::coefficient::ClosureState{
        .cells = system.celula,
        .lastCell = system.ncel,
        .globals = system.vg1dSP,
        .input = system.arq,
        .selectors = system.driftSelectors,
        .gasSurfaceTemperature = system.tGSup,
        .inletVoidFraction = system.alfE,
        .inletCompletionFraction = system.betaE,
        .inletPressure = system.presE,
        .inletTemperature = system.tempE,
        .steadyIteration = system.iterperm,
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
        .parallelSecondaryBranch = system.redeParalelaS,
        .parallelSecondaryBoundaryCondition = system.redeParalelaCCsecundario,
        .parallelSourceCells = system.indFonteRedeParalelaIni,
        .parallelSourceProductionLiquid = system.fonteMpRedeParalelaIni,
        .parallelSourceComplementaryLiquid = system.fonteMcRedeParalelaIni,
        .parallelSourceGas = system.fonteMgRedeParalelaIni,
    };
}

sisprod::thermal::ThermalState thermalStateOf(SProd &system) {
    return sisprod::thermal::ThermalState{
        .cells = system.celula,
        .gasCells = system.celulaG,
        .input = system.arq,
        .latentHeatTable = system.tables.HLat,
        .globals = system.vg1dSP,
        .thermalSourceDisabled = system.semTermo,
        .productionNetworkCoupled = system.verificaAcopRedeS,
        .primarySectionStart = system.SecPrimIniRedeP,
        .primarySectionEnd = system.SecPrimFimRedeP,
        .coupledCellIndices = system.acertaIndAcop,
        .poissonSolver3D = system.poisson3D,
        .lastCell = system.ncel,
        .surfaceChoke = system.chokeSup,
        .surfaceChokeMassFlag = system.masChkSup,
        .endNode = system.noextremo,
        .gasSurfaceTemperature = system.tGSup,
        .latentHeatEnabled = system.CalcLat,
        .sourceUpdater = {.system = system},
        .fullModel = system.modeloCompleto,
        .massTransferModel = system.TransMassModel,
        .closureUpdater = {.system = system},
        .inletPressure = system.presE,
        .inletTemperature = system.tempE,
        .inletQuality = system.titE,
        .inletVoidFraction = system.alfE,
        .inletCompletionFraction = system.betaE,
        .evolutionUpdater = {.system = system},
        .surfaceChokeOpen = system.aberto,
        .defaultInletTemperature = system.temperatura,
        .timeStep = system.dt,
        .minimumCycleTimeStep = system.dtCicMin,
        .poisson2DCellIndices = system.indCelPoisson2D,
        .poisson2DCellCount = system.nCelulaPoisson2D,
        .steadyIteration = system.iterperm,
        .slowHeatTransferThreshold = system.trocaTermicaLenta,
        .annulusTubingStart = system.gasLift.ColunaAnulaIni,
        .annulusTubingEnd = system.gasLift.ColunaAnulaFim,
        .tubingAnnulusStart = system.gasLift.AnulaColunaIni,
        .productionNetworkHeatCoupled = system.verificaAcopRedeP,
        .primaryNetworkSectionEnd = system.PrimSecFimRedeP,
        .primaryNetworkSectionStart = system.PrimSecIniRedeP,
        .gasCellCount = system.ncelGas,
        .gasLiftChokes = system.gasLift.chokeVGL,
        .gasSurfacePressure = system.pGSup,
        .outletPressure = system.presfim,
        .surfaceTemperature = system.tempSup,
        .networkCoupled = system.gasLift.verificaAcop,
        .tubingAnnulusEnd = system.gasLift.AnulaColunaFim,
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
/// Every field is a reference or a pointer, never a copy -- the caller advances
/// the counters between the header call and the row call, so a copy would be
/// read at the wrong moment.
trendoutput::TrendState trendStateOf(const SProd &system) {
    return trendoutput::TrendState{
        .input = system.arq,
        .globals = system.vg1dSP,
        .branchIndex = system.indTramo,
        .printPassCount = system.kimpT,
        .production = {system.MatTrendP, system.ntrend, system.ntrendB},
        .gasLine = {system.MatTrendG, system.ntrendg, system.ntrendgB},
        .productionCrossSection = {system.MatTrendTransP, system.ntrendtrans, system.ntrendtransB},
        .gasLineCrossSection = {system.MatTrendTransG, system.ntrendtransg, system.ntrendtransgB}};
}

}  // namespace sisprod::adapters
