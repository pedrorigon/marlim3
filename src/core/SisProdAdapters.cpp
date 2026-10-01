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
        .gasCells = system.gasLift.celulaG,
        .cells = system.celula,
        .input = system.arq,
        .globals = system.vg1dSP,
        .gasCellCount = system.gasLift.ncelGas,
        .lastCell = system.ncel,
        .gasLiftChokes = system.gasLift.chokeVGL.data(),
        .injectionChoke = system.gasLift.chokeInj,
        .gasValveCellIndices = system.gasLift.posicVGLG.data(),
        .productionValveCellIndices = system.gasLift.posicVGLP.data(),
        .gasSystemMatrix = system.gasLift.matglobG,
        .gasFreeTerms = system.gasLift.termolivreG,
        .annulusTubingStart = system.gasLift.ColunaAnulaIni,
        .annulusTubingEnd = system.gasLift.ColunaAnulaFim,
        .tubingAnnulusStart = system.gasLift.AnulaColunaIni,
        .tubingAnnulusEnd = system.gasLift.AnulaColunaFim,
        .steadyIteration = system.iterperm,
        .networkCoupled = system.gasLift.verificaAcop,
        .thermalSourceDisabled = system.semTermo,
        .initialGasPressure = system.gasLift.presiniG,
        .initialGasTemperature = system.gasLift.tempiniG,
        .gasSurfacePressure = system.gasLift.pGSup,
        .timeStep = system.dt,
        .interfaceCell = system.gasLift.celInter,
        .interfaceVelocity = system.gasLift.velInter,
        .interfaceTimeStep = system.gasLift.dtInter,
        .initialInterfaceCell = system.gasLift.celInterIni,
        .initialInterfaceVelocity = system.gasLift.velInterIni,
        .initialInterfaceTimeStep = system.gasLift.dtInterIni,
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
        .gasCells = system.gasLift.celulaG,
        .input = system.arq,
        .globals = system.vg1dSP,
        .lastCell = system.ncel,
        .gasCellCount = system.gasLift.ncelGas,
        .injectionChoke = system.gasLift.chokeInj,
        .surfaceChoke = system.chokeSup,
        .gasValveCellIndices = system.gasLift.posicVGLG.data(),
        .productionValveCellIndices = system.gasLift.posicVGLP.data(),
        .steadyIteration = system.iterperm,
        .searchOrigin = system.buscaIni,
        .convergenceMonitor = system.monitConvPerm,
        .baseConvergenceMonitor = system.monitConvPermBase,
        .annulusDrift = system.derivaAnel,
        .networkCoupled = system.gasLift.verificaAcop,
        .endNode = system.noextremo,
        .thermalSourceDisabled = system.semTermo,
        .slowHeatTransferThreshold = system.trocaTermicaLenta,
        .gasSurfacePressure = system.gasLift.pGSup,
        .initialGasPressure = system.gasLift.presiniG,
        .initialGasTemperature = system.gasLift.tempiniG,
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
        .meanMaximumTimeStep = system.transient.DTMaxMed,
        .meanMaximumPressureChange = system.transient.DpMaxMed,
        .masterState = system.transient.EstadoMaster1,
        .surfaceChokeOpen = system.transient.aberto,
        .initiallyOpen = system.transient.abertoini,
        .timeChanged = system.transient.alteraTempo,
        .inletCompletionFraction = system.betaE,
        .interfaceCell = system.gasLift.celInter,
        .masterCounter = system.transient.contaMaster1,
        .gasSpecificHeatTable = system.tables.cpg,
        .timeStep = system.dt,
        .meanCflTimeStep = system.transient.dtCFLMed,
        .totalCflTimeStep = system.transient.dtCFLTotal,
        .interfaceTimeStep = system.gasLift.dtInter,
        .meanSimulationTimeStep = system.transient.dtSimMed,
        .totalSimulationTimeStep = system.transient.dtSimTotal,
        .bufferedCompletionMassSource = system.transient.fontemassCRBuf,
        .bufferedGasMassSource = system.transient.fontemassGRBuf,
        .bufferedLiquidMassSource = system.transient.fontemassPRBuf,
        .eventIndex = system.transient.indevento,
        .slugCount = system.transient.kontaGolfada,
        .timeStepRestrictionCount = system.transient.kontarestriDt,
        .surfaceChokeMassFlag = system.masChkSup,
        .fullModel = system.transient.modeloCompleto,
        .desperationMoment = system.transient.momentoDesesp,
        .chokeModeChanged = system.transient.mudaModoChk,
        .multiplier = system.transient.mult,
        .outletPressure = system.presfim,
        .restart = system.transient.reinicia,
        .timeStepRestricted = system.transient.restriDt,
        .openTime = system.transient.tempoaberto,
        .productionSolution = system.transient.termolivreP,
        .inletQuality = system.titE,
        .masterRatio0 = system.transient.vRazMast0,
        .masterRatio1 = system.transient.vRazMast1,
        .masterCriticalRatio = system.transient.vRazMastCrit,
        .interfaceVelocity = system.gasLift.velInter,
        .masterOpenSchedule = system.transient.abreM1.data(),
        .input = system.arq,
        .initialInterfaceCell = system.gasLift.celInterIni,
        .cells = system.celula,
        .gasCells = system.gasLift.celulaG,
        .surfaceChoke = system.chokeSup,
        .cflTimeSteps = system.transient.dtCFL,
        .initialInterfaceTimeStep = system.gasLift.dtInterIni,
        .simulationTimeSteps = system.transient.dtSim,
        .auxiliaryCflTimeStep = system.transient.dtauxCFL,
        .finalAuxiliaryTimeStep = system.transient.dtauxFinal,
        .masterCloseSchedule = system.transient.fechaM1.data(),
        .productionFreeTerms = system.flut,
        .gasFreeTerms = system.flutG,
        .branchIndex = system.indTramo,
        .movingMeanFlux = system.transient.jMedMov,
        .stepIndex = system.transient.kSP,
        .productionMatrix = system.transient.matglobP,
        .smallestCellLength = system.transient.menorDx,
        .masterOpenCount = system.transient.nabreM1,
        .lastCell = system.ncel,
        .gasCellCount = system.gasLift.ncelGas,
        .productionCrossSectionCount = system.trends.ncelperftransp.data(),
        .masterCloseCount = system.transient.nfechaM1,
        .endNode = system.noextremo,
        .gasSurfacePressure = system.gasLift.pGSup,
        .inletPressure = system.presE,
        .movingMeanPressure = system.transient.presMedMov,
        .movingMeanTemperature = system.transient.tMedMov,
        .maximumTimeStepRates = system.transient.taxaDTMax,
        .maximumPressureRates = system.transient.taxaDpMax,
        .inletTemperature = system.tempE,
        .reverseQuality = system.networkCoupling.titRev,
        .initialInterfaceVelocity = system.gasLift.velInterIni,
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
        .trackGasOilRatio = system.transient.trackRGO,
        .trackGasGravity = system.transient.trackDeng,
        .compositionalRefreshCounter = system.transient.kontaRenovaComp,
        .surfaceChokeMassFlag = system.masChkSup,
        .movingPigCount = system.transient.indpigP,
        .previousMovingPigCount = system.transient.indpigPini,
        .scheduledPigCount = system.transient.npig,
        .pigReceiverCells = system.transient.receb.data(),
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
        .compositionalRefreshCounter = system.transient.kontaRenovaComp,
        .fluxHistory = system.transient.jVet,
        .gasCrossSectionCellCounts = system.trends.ncelperftransg.data(),
        .gasCrossSectionProfileTimeCounter = system.trends.kontaTempoTransProfG,
        .gasCrossSectionTrendBufferedCounts = system.trends.ntrendtransgB,
        .gasCrossSectionTrendCounts = system.trends.ntrendtransg,
        .gasCrossSectionTrendMatrix = system.trends.MatTrendTransG,
        .gasCrossSectionTrendResetTimers = system.trends.resettrendtransg,
        .gasProfileTimeCounter = system.trends.kontaTempoProfG,
        .gasTrendBufferedCounts = system.trends.ntrendgB,
        .gasTrendCounts = system.trends.ntrendg,
        .gasTrendMatrix = system.trends.MatTrendG,
        .gasTrendResetTimers = system.trends.resettrendg,
        .initialGasPressure = system.gasLift.presiniG,
        .initialGasSurfacePressure = system.transient.pGSupIni,
        .initialGasTemperature = system.gasLift.tempiniG,
        .initialOpenTime = system.transient.tempoabertoini,
        .injectionChoke = system.gasLift.chokeInj,
        .logBuffer = system.tmpLog,
        .logCounter = system.transient.contaLog,
        .massTransferModel = system.transient.TransMassModel,
        .minimumCycleTimeStep = system.transient.dtCicMin,
        .movingMeanCounter = system.transient.ktMedMov,
        .movingMeanVoidFraction = system.transient.alfMedMov,
        .networkCoupled = system.gasLift.verificaAcop,
        .poissonSolver3D = system.poisson3D,
        .pressureHistory = system.transient.presVet,
        .printCounter = system.transient.KontaImprime,
        .printPassCount = system.trends.kimpT,
        .productionCrossSectionProfileTimeCounter = system.trends.kontaTempoTransProf,
        .productionCrossSectionTrendBufferedCounts = system.trends.ntrendtransB,
        .productionCrossSectionTrendCounts = system.trends.ntrendtrans,
        .productionCrossSectionTrendMatrix = system.trends.MatTrendTransP,
        .productionCrossSectionTrendResetTimers = system.trends.resettrendtrans,
        .productionProfileTimeCounter = system.trends.kontaTempoProf,
        .productionTrendBufferedCounts = system.trends.ntrendB,
        .productionTrendCounts = system.trends.ntrend,
        .productionTrendMatrix = system.trends.MatTrendP,
        .productionTrendResetTimers = system.trends.resettrend,
        .startNode = system.noinicial,
        .temperatureHistory = system.transient.tVet,
        .totalFlux = system.transient.jTotal,
        .totalPressure = system.transient.pTotal,
        .totalVoidFraction = system.transient.alfTotal,
        .trackGasGravity = system.transient.trackDeng,
        .trackGasOilRatio = system.transient.trackRGO,
        .unitCellTimeCounters = system.trends.kontaTempoCelUni,
        .voidFractionHistory = system.transient.alfVet,
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
        .parallelSecondaryBranch = system.networkCoupling.redeParalelaS,
        .parallelSecondaryBoundaryCondition = system.networkCoupling.redeParalelaCCsecundario,
        .parallelSourceCells = system.networkCoupling.indFonteRedeParalelaIni,
        .parallelSourceProductionLiquid = system.networkCoupling.fonteMpRedeParalelaIni,
        .parallelSourceComplementaryLiquid = system.networkCoupling.fonteMcRedeParalelaIni,
        .parallelSourceGas = system.networkCoupling.fonteMgRedeParalelaIni,
    };
}

sisprod::thermal::ThermalState thermalStateOf(SProd &system) {
    return sisprod::thermal::ThermalState{
        .cells = system.celula,
        .gasCells = system.gasLift.celulaG,
        .input = system.arq,
        .latentHeatTable = system.tables.HLat.data(),
        .globals = system.vg1dSP,
        .thermalSourceDisabled = system.semTermo,
        .productionNetworkCoupled = system.networkCoupling.verificaAcopRedeS,
        .primarySectionStart = system.networkCoupling.SecPrimIniRedeP,
        .primarySectionEnd = system.networkCoupling.SecPrimFimRedeP,
        .coupledCellIndices = system.acertaIndAcop,
        .poissonSolver3D = system.poisson3D,
        .lastCell = system.ncel,
        .surfaceChoke = system.chokeSup,
        .surfaceChokeMassFlag = system.masChkSup,
        .endNode = system.noextremo,
        .gasSurfaceTemperature = system.tGSup,
        .latentHeatEnabled = system.CalcLat,
        .sourceUpdater = {.system = system},
        .fullModel = system.transient.modeloCompleto,
        .massTransferModel = system.transient.TransMassModel,
        .closureUpdater = {.system = system},
        .inletPressure = system.presE,
        .inletTemperature = system.tempE,
        .inletQuality = system.titE,
        .inletVoidFraction = system.alfE,
        .inletCompletionFraction = system.betaE,
        .evolutionUpdater = {.system = system},
        .surfaceChokeOpen = system.transient.aberto,
        .defaultInletTemperature = system.temperatura,
        .timeStep = system.dt,
        .minimumCycleTimeStep = system.transient.dtCicMin,
        .poisson2DCellIndices = system.transient.indCelPoisson2D,
        .poisson2DCellCount = system.transient.nCelulaPoisson2D,
        .steadyIteration = system.iterperm,
        .slowHeatTransferThreshold = system.trocaTermicaLenta,
        .annulusTubingStart = system.gasLift.ColunaAnulaIni,
        .annulusTubingEnd = system.gasLift.ColunaAnulaFim,
        .tubingAnnulusStart = system.gasLift.AnulaColunaIni,
        .productionNetworkHeatCoupled = system.networkCoupling.verificaAcopRedeP,
        .primaryNetworkSectionEnd = system.networkCoupling.PrimSecFimRedeP,
        .primaryNetworkSectionStart = system.networkCoupling.PrimSecIniRedeP,
        .gasCellCount = system.gasLift.ncelGas,
        .gasLiftChokes = system.gasLift.chokeVGL.data(),
        .gasSurfacePressure = system.gasLift.pGSup,
        .outletPressure = system.presfim,
        .surfaceTemperature = system.networkCoupling.tempSup,
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
        .printPassCount = system.trends.kimpT,
        .production = {system.trends.MatTrendP, system.trends.ntrend, system.trends.ntrendB},
        .gasLine = {system.trends.MatTrendG, system.trends.ntrendg, system.trends.ntrendgB},
        .productionCrossSection = {system.trends.MatTrendTransP, system.trends.ntrendtrans, system.trends.ntrendtransB},
        .gasLineCrossSection = {system.trends.MatTrendTransG, system.trends.ntrendtransg, system.trends.ntrendtransgB}};
}

}  // namespace sisprod::adapters
