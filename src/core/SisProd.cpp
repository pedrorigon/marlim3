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
#include "SisProdComposition.h"
#include "SisProdSolveContext.h"
#include "SisProdThermal.h"
#include "SisProdTrendOutput.h"

using namespace sisprod;

void SProd::resolveDriftSelectors() {
    driftSelectors = {arq.CorreDisper, arq.CorreAnular, arq.CorreEstrat};
}

void SProd::resetRunState() {
    tables = PropertyTables();
    steadySearch = SteadySearch();

    CalcLat = 0;
    celula = nullptr;
    presE = -1;
    tempE = -1;
    titE = -1;
    betaE = -1;
    alfE = -1;
    tGSup = 0.;
    thermalSourceDisabled = 0;

    gasLift.celulaG = nullptr;
    gasLift.celInter = kNoInterfaceCell;
    gasLift.dtInter = 0.;
    gasLift.velInter = 0.;
    gasLift.celInterIni = gasLift.celInter;
    gasLift.dtInterIni = gasLift.dtInter;
    gasLift.velInterIni = gasLift.velInter;
    gasLift.tempMedContDesc = kUnloadingAveragingWindow;
    gasLift.maxVecContDesc = kUnloadingSampleCapacity;
    gasLift.vazmedDesc = 0;
    gasLift.tempmedDEsc = 0;
    gasLift.vazmaxMedDesc.clear();
    gasLift.dtDesc.clear();

    transient.trackRGO = 0;
    transient.trackDeng = 0;
    transient.fontemassPRBuf = 0.;
    transient.fontemassCRBuf = 0.;
    transient.fontemassGRBuf = 0.;
    transient.dtCFLMed = 1.;
    transient.dtSimMed = 1.;
    transient.restriDt = 0;
    transient.kontarestriDt = 0;
    transient.dtauxCFL = 0.;
    transient.dtauxFinal = 0.;
    transient.kontaGolfada = kInitialSlugCount;
    transient.mudaModoChk = 0;
    transient.momentoDesesp = 0;
    transient.modeloCompleto = 1;
    transient.DpMaxMed = 1.;
    transient.DTMaxMed = 1.;
    for (int i = 0; i < kMasterRatioHistorySize; i++) {
        transient.vRazMast0[i] = 0.;
        transient.vRazMast1[i] = 0.;
        transient.vRazMastCrit[i] = kDefaultMasterCriticalRatio;
    }
    transient.kontaRenovaComp = 0;
    transient.nCelulaPoisson2D = 0;
    transient.alteraTempo = 0;
    transient.presVet.clear();
    transient.jVet.clear();
    transient.alfVet.clear();
    transient.tVet.clear();
    transient.dtSim.clear();
    transient.dtCFL.clear();
    transient.taxaDpMax.clear();
    transient.taxaDTMax.clear();
    transient.indCelPoisson2D.clear();

    trends.gasTrendMatrices = nullptr;
    trends.gasTrendResetTimes = nullptr;
    trends.gasTrendCounts = nullptr;
    trends.gasTrendBufferedCounts = nullptr;
    trends.gasWallTrendMatrices = nullptr;
    trends.gasWallTrendResetTimes = nullptr;
    trends.gasWallTrendCounts = nullptr;
    trends.gasWallTrendBufferedCounts = nullptr;
    trends.productionTrendMatrices = nullptr;
    trends.productionTrendResetTimes = nullptr;
    trends.productionTrendCounts = nullptr;
    trends.productionTrendBufferedCounts = nullptr;
    trends.productionWallTrendMatrices = nullptr;
    trends.productionWallTrendResetTimes = nullptr;
    trends.productionWallTrendCounts = nullptr;
    trends.productionWallTrendBufferedCounts = nullptr;
    trends.trendOutputCount = 0.;

    networkCoupling.previousInletPressure = -1;
    networkCoupling.previousInletTemperature = -1;
    networkCoupling.previousInletQuality = -1;
    networkCoupling.previousInletLiquidFraction = -1;
    networkCoupling.previousInletVoidFraction = -1;
    networkCoupling.previousSeparatorTemperature = 0.;
    networkCoupling.previousChokeModeChanged = 0;
    networkCoupling.parallelSourceCells.clear();
    networkCoupling.parallelSourceProductionLiquid.clear();
    networkCoupling.parallelSourceComplementaryLiquid.clear();
    networkCoupling.parallelSourceGas.clear();

    coupling3D.coupledIndexCorrections.clear();
}

void SProd::allocateSolutionBuffers() {
    flut = FullMtx<double>(arq.ncelp, arq.nvarprofp + kProductionProfileFixedColumns);
    flutG = FullMtx<double>(arq.ncelg, arq.nvarprofg + kGasLineProfileFixedColumns);
    transient.matglobP =
        BandMtx<double>(kProductionUnknownsPerCell * arq.ncelp, kProductionLeftBandwidth, kProductionRightBandwidth);
    transient.termolivreP = Vcr<double>(kProductionUnknownsPerCell * arq.ncelp);
    gasLift.matglobG = BandMtx<double>(kGasLineUnknownsPerCell * arq.ncelg, kGasLineBandwidth, kGasLineBandwidth);
    gasLift.termolivreG = Vcr<double>(kGasLineUnknownsPerCell * arq.ncelg);
}

SProd::SProd(string inputFile, string logFile, tipoValidacaoJson_t jsonValidation,
             tipoSimulacao_t simulationType, varGlob1D *globals, int dynamicTable, int blocking, int reversedGeometry,
             double *extraSourceLengths, int *extraSourceCells, int extraSourceCount, int steadyBranch)
    : networkCoupling{.previousReverseLiquidFraction = 0},
      arq(inputFile, logFile, jsonValidation, simulationType, reversedGeometry, globals, steadyBranch),
      flutG(arq.ncelg, arq.nvarprofg + kGasLineProfileFixedColumns),
      flut(arq.ncelp, arq.nvarprofp + kProductionProfileFixedColumns),
      gasLift{.matglobG = BandMtx<double>(kGasLineUnknownsPerCell * arq.ncelg, kGasLineBandwidth, kGasLineBandwidth),
              .termolivreG = Vcr<double>(kGasLineUnknownsPerCell * arq.ncelg)},
      transient{.matglobP = BandMtx<double>(kProductionUnknownsPerCell * arq.ncelp, kProductionLeftBandwidth,
                                            kProductionRightBandwidth),
                .termolivreP = Vcr<double>(kProductionUnknownsPerCell * arq.ncelp)},
      bloq(blocking), vg1dSP(globals) {
    resolveDriftSelectors();
    resetRunState();
    if (dynamicTable >= 0)
        arq.tabelaDinamica = dynamicTable;
    montasistema(extraSourceLengths, extraSourceCells, extraSourceCount);
}

SProd::SProd() {
    resolveDriftSelectors();
}

namespace {

/// Destroys the elements of v and gives its memory back.
template <typename T>
void releaseVector(vector<T> &v) {
    vector<T>().swap(v);
}

}  // namespace

void SProd::releaseOwnedStorage() {
    releaseVector(gasLift.gasCells);
    releaseVector(gasLift.chokeVGL);
    releaseVector(gasLift.posicVGLP);
    releaseVector(gasLift.posicVGLG);
    releaseVector(transient.abreM1);
    releaseVector(transient.fechaM1);

    releaseVector(trends.productionWallProfileCells);
    releaseVector(trends.gasWallProfileCells);

    trends.productionTrendSet.release();
    trends.gasTrendSet.release();
    trends.productionWallTrendSet.release();
    trends.gasWallTrendSet.release();

    tables.latentHeatTable.release();
    releaseVector(tables.pvtSimBubblePressures);
    releaseVector(tables.pvtSimBubbleTemperatures);
    tables.liviaSolutionGasRatioTable.release();

    releaseVector(productionCells);
    releaseVector(transient.receb);
    releaseVector(tables.dynamicTables);
}

SProd &SProd::operator=(const SProd &other) {
    releaseOwnedStorage();

    arq = other.arq;
    vg1dSP = other.vg1dSP;
    resolveDriftSelectors();
    allocateSolutionBuffers();
    resetRunState();
    applyCarriedState(other.carriedState());
    networkCoupling.parallelSecondaryBoundaryCondition = other.networkCoupling.parallelSecondaryBoundaryCondition;
    networkCoupling.parallelPrimaryBranch = other.networkCoupling.parallelPrimaryBranch;
    networkCoupling.parallelSecondaryBranch = other.networkCoupling.parallelSecondaryBranch;

    montasistema();

    return *this;
}

SProd::SProd(Ler &parsedInput, const CarriedState &carried) : vg1dSP(parsedInput.vg1dSP) {
    arq.copiaSemJson(parsedInput);
    resolveDriftSelectors();
    allocateSolutionBuffers();
    resetRunState();
    applyCarriedState(carried);

    montasistema();
}

SProd::CarriedState SProd::carriedState() const {
    return {.noextremo = noextremo,
            .noinicial = noinicial,
            .derivaAnel = derivaAnel,
            .bloq = bloq,
            .previousReverseLiquidFraction = networkCoupling.previousReverseLiquidFraction,
            .previousReverseQuality = networkCoupling.previousReverseQuality,
            .dtCicMin = transient.dtCicMin};
}

void SProd::applyCarriedState(const CarriedState &carried) {
    noextremo = carried.noextremo;
    noinicial = carried.noinicial;
    derivaAnel = carried.derivaAnel;
    bloq = carried.bloq;
    networkCoupling.previousReverseLiquidFraction = carried.previousReverseLiquidFraction;
    networkCoupling.previousReverseQuality = carried.previousReverseQuality;
    transient.dtCicMin = carried.dtCicMin;
}

void SProd::HidroDescargaG() {
    SolveContext context(*this);
    gaslift::computeGasUnloadingHydrostatics(context.gasLift());
}

void SProd::HidroDescargaP() {
    SolveContext context(*this);
    gaslift::computeProductionUnloadingHydrostatics(context.gasLift());
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

void SProd::aberturaVal0() {
    SolveContext context(*this);
    transient::valveOpeningLow(context.transientStep());
}

void SProd::aberturaVal1() {
    SolveContext context(*this);
    transient::valveOpeningHigh(context.transientStep());
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

void SProd::solveHydrateEnvelopes() {
    if (!(arq.calculaEnvelope == 1 && vg1dSP->lixo5 <= arq.tfinal))
        return;

    FA_Hidrato solverHidrato(*this);
    solverHidrato.solverHidrato();

    if (arq.lingas > 0) {
        FA_Hidrato_Servico solverHidratoG(*this);
        solverHidratoG.solverHidratoG();
    }
}

void SProd::SolveTrans(double reverseQuality, double reverseVoidFraction, double reverseLiquidFraction, int network,
                       const ProFlu *collectorFluid) {
    SolveContext context(*this);
    transient::solveTransientStep(context.transientSolve(), reverseQuality, reverseVoidFraction, reverseLiquidFraction,
                                  network, collectorFluid);
}

void SProd::ImprimeTrendPCab(int trend, int network) {
    SolveContext context(*this);
    trendoutput::writeProductionTrendHeader(context.trends(), trend, network);
}

void SProd::ImprimeTrendP(int trend, int network) {
    SolveContext context(*this);
    trendoutput::writeProductionTrendRows(context.trends(), trend, network);
}

void SProd::ImprimeTrendGCab(int trend, int network) {
    SolveContext context(*this);
    trendoutput::writeGasLineTrendHeader(context.trends(), trend, network);
}

void SProd::ImprimeTrendG(int trend, int network) {
    SolveContext context(*this);
    trendoutput::writeGasLineTrendRows(context.trends(), trend, network);
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

void SProd::calcTempFim() {
    SolveContext context(*this);
    thermal::computeOutletTemperature(context.thermal());
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

