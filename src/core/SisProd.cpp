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
