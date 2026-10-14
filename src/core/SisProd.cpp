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
    tables.zdranP = nullptr;
    tables.dzdpP = nullptr;
    tables.dzdtP = nullptr;
    tables.cpg = nullptr;
    tables.cpl = nullptr;
    tables.drholdT = nullptr;
    steadySearch.nfluP = 0;
    CalcLat = 0;
    transient.trackRGO = 0;
    transient.trackDeng = 0;
    gasLift.celulaG = nullptr;
    celula = nullptr;
    gasLift.celInter = kNoInterfaceCell;
    gasLift.dtInter = 0.;
    gasLift.velInter = 0.;

    trends.MatTrendG = nullptr;
    trends.resettrendg = nullptr;
    trends.ntrendg = nullptr;
    trends.ntrendgB = nullptr;
    trends.MatTrendTransG = nullptr;
    trends.resettrendtransg = nullptr;
    trends.ntrendtransg = nullptr;
    trends.ntrendtransgB = nullptr;
    trends.MatTrendP = nullptr;
    trends.resettrend = nullptr;
    trends.ntrend = nullptr;
    trends.ntrendB = nullptr;
    trends.MatTrendTransP = nullptr;
    trends.resettrendtrans = nullptr;
    trends.ntrendtrans = nullptr;
    trends.ntrendtransB = nullptr;

    transient.fontemassPRBuf = 0.;
    transient.fontemassCRBuf = 0.;
    transient.fontemassGRBuf = 0.;

    presE = -1;
    tempE = -1;
    titE = -1;
    betaE = -1;
    alfE = -1;
    networkCoupling.presEini = -1;
    networkCoupling.tempEini = -1;
    networkCoupling.titEini = -1;
    networkCoupling.betaEini = -1;
    networkCoupling.alfEini = -1;

    gasLift.tempMedContDesc = kUnloadingAveragingWindow;
    gasLift.maxVecContDesc = kUnloadingSampleCapacity;
    gasLift.vazmedDesc = 0;
    gasLift.tempmedDEsc = 0;

    tGSup = 0.;
    networkCoupling.tGSupIni = 0.;

    transient.dtCFLMed = 1.;
    transient.dtSimMed = 1.;
    transient.restriDt = 0;
    transient.kontarestriDt = 0;
    transient.dtauxCFL = 0.;
    transient.dtauxFinal = 0.;

    trends.kimpT = 0.;

    transient.kontaGolfada = kInitialSlugCount;

    transient.mudaModoChk = 0;
    networkCoupling.mudaModoChkini = 0;

    transient.momentoDesesp = 0;

    transient.modeloCompleto = 1;

    transient.DpMaxMed = 1.;
    transient.DTMaxMed = 1.;

    steadySearch.chuteHol = -1.;

    steadySearch.buscaIni = 0;

    for (int i = 0; i < kMasterRatioHistorySize; i++) {
        transient.vRazMast0[i] = 0.;
        transient.vRazMast1[i] = 0.;
        transient.vRazMastCrit[i] = kDefaultMasterCriticalRatio;
    }

    transient.kontaRenovaComp = 0;

    steadySearch.fluiRevRede = ProFlu();
    steadySearch.tempRev = 0.;
    steadySearch.revPerm = 0;
    tables.ntabDin = 0;

    transient.nCelulaPoisson2D = 0;
    steadySearch.trocaTermicaLenta = kSlowHeatTransferThreshold;

    semTermo = 0;

    steadySearch.monitConvPerm = kInitialConvergenceMonitor;
    steadySearch.monitConvPermBase = 1.;

    transient.alteraTempo = 0;
}

SProd::SProd(string nomeArquivoEntrada, string nomeArquivoLog, tipoValidacaoJson_t validacaoJson,
             tipoSimulacao_t tipoSimulacao, varGlob1D *Vvg1dSP, int TD, int vbloq, int temporario, int reverso,
             double *compfonte, int *posicfonte, int nfontes, int redeperm)
    : arq(nomeArquivoEntrada, nomeArquivoLog, validacaoJson, tipoSimulacao, reverso, Vvg1dSP, redeperm),
      flutG(arq.ncelg, arq.nvarprofg + kGasLineProfileFixedColumns),
      flut(arq.ncelp, arq.nvarprofp + kProductionProfileFixedColumns),
      gasLift{.matglobG = BandMtx<double>(kGasLineUnknownsPerCell * arq.ncelg, kGasLineBandwidth, kGasLineBandwidth),
              .termolivreG = Vcr<double>(kGasLineUnknownsPerCell * arq.ncelg)},
      transient{.matglobP = BandMtx<double>(kProductionUnknownsPerCell * arq.ncelp, kProductionLeftBandwidth,
                                            kProductionRightBandwidth),
                .termolivreP = Vcr<double>(kProductionUnknownsPerCell * arq.ncelp)} {
    resolveDriftSelectors();
    resetRunState();
    gasLift.celInterIni = gasLift.celInter;
    gasLift.dtInterIni = gasLift.dtInter;
    gasLift.velInterIni = gasLift.velInter;
    redeTemporario = temporario;
    networkCoupling.betaRevini = 0;
    bloq = vbloq;
    vg1dSP = Vvg1dSP;
    if (TD >= 0)
        arq.tabelaDinamica = TD;
    montasistema(compfonte, posicfonte, nfontes);
}

SProd::SProd()
    : arq(), flutG(1, kGasLineProfileFixedColumns), flut(1, kProductionProfileFixedColumns),
      gasLift{.matglobG = BandMtx<double>(kGasLineUnknownsPerCell, kGasLineBandwidth, kGasLineBandwidth),
              .termolivreG = Vcr<double>(kGasLineUnknownsPerCell)},
      transient{.matglobP = BandMtx<double>(kProductionUnknownsPerCell, kProductionLeftBandwidth,
                                            kProductionRightBandwidth),
                .termolivreP = Vcr<double>(kProductionUnknownsPerCell)} {
    resolveDriftSelectors();
    resetRunState();
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

    releaseVector(trends.ncelperftransp);
    releaseVector(trends.ncelperftransg);

    trends.productionTrendSet.release();
    trends.gasTrendSet.release();
    trends.productionWallTrendSet.release();
    trends.gasWallTrendSet.release();

    tables.HLat.release();
    releaseVector(tables.PBPVTSim);
    releaseVector(tables.TBPVTSim);
    tables.RSLivia.release();

    releaseVector(productionCells);
    releaseVector(transient.receb);

    if (arq.tabelaDinamica == 1) {
        tables.tabDin.clear();
    }
}

SProd &SProd::operator=(const SProd &sp) {
    releaseOwnedStorage();

    arq = sp.arq;
    resolveDriftSelectors();
    flut = sp.flut;
    flutG = sp.flutG;
    transient.matglobP = sp.transient.matglobP;
    transient.termolivreP = sp.transient.termolivreP;
    gasLift.matglobG = sp.gasLift.matglobG;
    gasLift.termolivreG = sp.gasLift.termolivreG;
    vg1dSP = sp.vg1dSP;
    resetRunState();
    gasLift.celInterIni = gasLift.celInter;
    gasLift.dtInterIni = gasLift.dtInter;
    gasLift.velInterIni = gasLift.velInter;

    noextremo = sp.noextremo;
    noinicial = sp.noinicial;
    derivaAnel = sp.derivaAnel;

    networkCoupling.betaRevini = sp.networkCoupling.betaRevini;
    networkCoupling.titRevini = sp.networkCoupling.titRevini;

    bloq = sp.bloq;

    transient.dtCicMin = sp.transient.dtCicMin;
    networkCoupling.redeParalelaCCsecundario = sp.networkCoupling.redeParalelaCCsecundario;
    networkCoupling.redeParalelaP = sp.networkCoupling.redeParalelaP;
    networkCoupling.redeParalelaS = sp.networkCoupling.redeParalelaS;

    transient.presVet.clear();
    transient.jVet.clear();
    transient.alfVet.clear();
    transient.tVet.clear();
    transient.dtSim.clear();
    transient.dtCFL.clear();
    gasLift.vazmaxMedDesc.clear();
    gasLift.dtDesc.clear();
    transient.taxaDpMax.clear();
    transient.taxaDTMax.clear();
    tables.tabDin.clear();
    coupling3D.acertaIndAcop.clear();
    transient.indCelPoisson2D.clear();
    networkCoupling.indFonteRedeParalelaIni.clear();
    networkCoupling.fonteMpRedeParalelaIni.clear();
    networkCoupling.fonteMcRedeParalelaIni.clear();
    networkCoupling.fonteMgRedeParalelaIni.clear();

    montasistema();

    return *this;
}

SProd::SProd(Ler &parsedInput, const CarriedState &carried) : SProd() {
    releaseOwnedStorage();

    arq.copiaSemJson(parsedInput);
    resolveDriftSelectors();
    flut = FullMtx<double>(arq.ncelp, arq.nvarprofp + kProductionProfileFixedColumns);
    flutG = FullMtx<double>(arq.ncelg, arq.nvarprofg + kGasLineProfileFixedColumns);
    transient.matglobP =
        BandMtx<double>(kProductionUnknownsPerCell * arq.ncelp, kProductionLeftBandwidth, kProductionRightBandwidth);
    transient.termolivreP = Vcr<double>(kProductionUnknownsPerCell * arq.ncelp);
    gasLift.matglobG = BandMtx<double>(kGasLineUnknownsPerCell * arq.ncelg, kGasLineBandwidth, kGasLineBandwidth);
    gasLift.termolivreG = Vcr<double>(kGasLineUnknownsPerCell * arq.ncelg);
    vg1dSP = arq.vg1dSP;
    resetRunState();
    gasLift.celInterIni = gasLift.celInter;
    gasLift.dtInterIni = gasLift.dtInter;
    gasLift.velInterIni = gasLift.velInter;

    noextremo = carried.noextremo;
    noinicial = carried.noinicial;
    derivaAnel = carried.derivaAnel;

    networkCoupling.betaRevini = carried.betaRevini;
    networkCoupling.titRevini = carried.titRevini;

    bloq = carried.bloq;

    transient.dtCicMin = carried.dtCicMin;

    networkCoupling.redeParalelaCCsecundario = -1;
    networkCoupling.redeParalelaP = -1;
    networkCoupling.redeParalelaS = -1;

    montasistema();
}

SProd::CarriedState SProd::carriedState() const {
    return {.noextremo = noextremo,
            .noinicial = noinicial,
            .derivaAnel = derivaAnel,
            .bloq = bloq,
            .betaRevini = networkCoupling.betaRevini,
            .titRevini = networkCoupling.titRevini,
            .dtCicMin = transient.dtCicMin};
}

void SProd::HidroDescargaG() {
    SolveContext context(*this);
    gaslift::computeGasUnloadingHydrostatics(context.gasLift());
}

void SProd::HidroDescargaP() {
    SolveContext context(*this);
    gaslift::computeProductionUnloadingHydrostatics(context.gasLift());
}

void SProd::renova(int expli) {
    SolveContext context(*this);
    transient::updateCells(context.transientStep(), expli);
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

void SProd::renovaterm(int aflu) {
    SolveContext context(*this);
    thermal::updateFlowPartitionTerms(context.thermal(), aflu);
}

void SProd::renovatermAfluFim() {
    SolveContext context(*this);
    thermal::updateOutletFlowPartitionTerms(context.thermal());
}

void SProd::renovatermColIni() {
    SolveContext context(*this);
    thermal::updateInletFlowPartitionTerms(context.thermal());
}

void SProd::calcCCpres(double titRev, double alfRev, double betRev) {
    SolveContext context(*this);
    transient::applyOutletPressureCondition(context.transientStep(), titRev, alfRev, betRev);
}

void SProd::calcCCBuffer(double titRev, double alfRev, double betRev) {
    SolveContext context(*this);
    transient::applyOutletBufferCondition(context.transientStep(), titRev, alfRev, betRev);
}

void SProd::determinaDT(int vexpli) {
    SolveContext context(*this);
    transient::computeTimeStep(context.transientStep(), vexpli);
}

void SProd::atenuaDtMax() {
    SolveContext context(*this);
    transient::dampMaximumTimeStep(context.transientStep());
}

void SProd::avaliaVariaDpDt(int vexpli) {
    SolveContext context(*this);
    transient::evaluatePressureRateOfChange(context.transientStep(), vexpli);
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

void SProd::EvoluiFrac(double alfrev, double betrev, int ciclo) {
    SolveContext context(*this);
    transient::evolveFractions(context.transientStep(), alfrev, betrev, ciclo);
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

void SProd::SolveAcopPV(int vexpli, int ciclo) {
    SolveContext context(*this);
    transient::solvePressureVolumeCoupling(context.transientStep(), vexpli, ciclo);
}

void SProd::marchaEnergTrans(int ciclo, int ciclomax) {
    SolveContext context(*this);
    thermal::advanceTransientEnergy(context.thermal(), ciclo, ciclomax);
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
    if (arq.calculaEnvelope == 1 && vg1dSP->lixo5 <= arq.tfinal) {
        FA_Hidrato solverHidrato(*this);
        solverHidrato.solverHidrato();
    }

    if (arq.lingas > 0 && arq.calculaEnvelope == 1 && vg1dSP->lixo5 <= arq.tfinal) {
        FA_Hidrato_Servico solverHidratoG(*this);
        solverHidratoG.solverHidratoG();
    }
}

void SProd::SolveTrans(double titRev, double alfRev, double betRev, int nrede, ProFlu fluiRev) {
    SolveContext context(*this);
    transient::solveTransientStep(context.transientSolve(), titRev, alfRev, betRev, nrede, fluiRev);
}

void SProd::ImprimeTrendPCab(int i, int nrede) {
    SolveContext context(*this);
    trendoutput::writeProductionTrendHeader(context.trends(), i, nrede);
}

void SProd::ImprimeTrendP(int i, int nrede) {
    SolveContext context(*this);
    trendoutput::writeProductionTrendRows(context.trends(), i, nrede);
}

void SProd::ImprimeTrendGCab(int i, int nrede) {
    SolveContext context(*this);
    trendoutput::writeGasLineTrendHeader(context.trends(), i, nrede);
}

void SProd::ImprimeTrendG(int i, int nrede) {
    SolveContext context(*this);
    trendoutput::writeGasLineTrendRows(context.trends(), i, nrede);
}

double SProd::marchaProdPerm1(double pchute) {
    SolveContext context(*this);
    return steady::marchProductionSteady(context.steady(), pchute);
}

double SProd::buscaProdPfundoPerm(double chute, int kontaTenta) {
    SolveContext context(*this);
    return steady::searchProductionBottomHolePressure(context.search(), chute, kontaTenta);
}

double SProd::buscaProdPfundoPermRev(double chute) {
    SolveContext context(*this);
    return steady::searchReverseProductionBottomHolePressure(context.search(), chute);
}

double SProd::buscaProdPfundoPerm2(double chute, int kontaTenta) {
    SolveContext context(*this);
    return steady::searchProductionBottomHolePressureSecondary(context.search(), chute, kontaTenta);
}

double SProd::buscaProdPfundoPerm3(double pentrada) {
    SolveContext context(*this);
    return steady::searchProductionBottomHolePressureTertiary(context.search(), pentrada);
}

double SProd::marchaProdPresPres1(double mchute) {
    SolveContext context(*this);
    return steady::marchProductionPressureToPressure(context.steady(), mchute);
}

double SProd::buscaProdPresPresPerm(double mchute, double maxvaz, int kontaiter) {
    SolveContext context(*this);
    return steady::searchProductionPressureToPressure(context.search(), mchute, maxvaz, kontaiter);
}

double SProd::buscaProdPresPresPermRev(double mchute, double maxvaz, int kontaiter) {
    SolveContext context(*this);
    return steady::searchReverseProductionPressureToPressure(context.search(), mchute, maxvaz, kontaiter);
}

double SProd::buscaProdPresPresPerm2(double mchute, double maxvaz) {
    SolveContext context(*this);
    return steady::searchProductionPressureToPressureSecondary(context.search(), mchute, maxvaz);
}

double SProd::buscaProdPresPresPerm3(double mchute, double maxvaz) {
    SolveContext context(*this);
    return steady::searchProductionPressureToPressureTertiary(context.search(), mchute, maxvaz);
}

void SProd::calcTempFim() {
    SolveContext context(*this);
    thermal::computeOutletTemperature(context.thermal());
}

double SProd::buscaInjPfundoPerm1(double chute) {
    SolveContext context(*this);
    return steady::searchInjectionBottomHolePressure1(context.search(), chute);
}

double SProd::buscaInjPfundoPerm2(double chute) {
    SolveContext context(*this);
    return steady::searchInjectionBottomHolePressure2(context.search(), chute);
}

double SProd::buscaInjPfundoPerm3(double chute) {
    SolveContext context(*this);
    return steady::searchInjectionBottomHolePressure3(context.search(), chute);
}

double SProd::buscaInjPfundoPerm4() {
    SolveContext context(*this);
    return steady::searchInjectionBottomHolePressure4(context.search());
}

double SProd::buscaInjPfundoPerm5(double chute) {
    SolveContext context(*this);
    return steady::searchInjectionBottomHolePressure5(context.search(), chute);
}

double SProd::hidroreverso(double hol, double vaz, double vazG) {
    SolveContext context(*this);
    return steady::reverseHydrostatic(context.steady(), hol, vaz, vazG);
}

double SProd::hidroreversoInj(double hol, double vaz) {
    SolveContext context(*this);
    return steady::reverseInjectionHydrostatic(context.steady(), hol, vaz);
}

double SProd::hidroTramoSecundario(double titulo) {
    SolveContext context(*this);
    return steady::secondaryBranchHydrostatic(context.steady(), titulo);
}

