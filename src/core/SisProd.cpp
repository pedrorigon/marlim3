/*
 * SisProd.cpp
 *
 *  Created on: 21 Dec 2016
 *      Author: Eduardo
 */
#define _USE_MATH_DEFINES // for M_PI
#include "SisProd.h"
#include "DriftFluxClosure.h"
#include "FA_Hidratos.h"
#include "FA_Hidratos_Servico.h"
#include "OutputI18n.h"
#include "RootFindingSolvers.h"
#include "SisProdConstants.h"
#include "SisProdGasLift.h"
#include "SisProdSources.h"
#include "SisProdSteadyState.h"
#include "SisProdSteadyStateSearch.h"
#include "SisProdTransient.h"
#include "SisProdComposition.h"
#include "SisProdThermal.h"
#include "SisProdTrendOutput.h"
#include <chrono>
#include <math.h>

// The state adapters live in sisprod::adapters, where SisProd.h declares them
// so that SProd can name them as friends.
using sisprod::adapters::gasLiftStateOf;
using sisprod::adapters::steadyStateOf;
using sisprod::adapters::searchStateOf;
using sisprod::adapters::transientStateOf;
using sisprod::adapters::compositionStateOf;
using sisprod::adapters::transientSolveStateOf;
using sisprod::adapters::thermalStateOf;
using sisprod::adapters::sourceStateOf;
using sisprod::adapters::closureStateOf;
using sisprod::adapters::trendStateOf;
using sisprod::kGravity;

void SProd::resolveDriftSelectors() {
    driftSelectors = {arq.CorreDisper, arq.CorreAnular, arq.CorreEstrat};
}

/// The run state every construction and reassignment starts from, before each
/// path sets what is its own and montasistema builds the system.
void SProd::resetRunState() {
    tables.zdranP = 0;
    tables.dzdpP = 0;
    tables.dzdtP = 0;
    tables.cpg = 0;
    tables.cpl = 0;
    tables.drholdT = 0;
    nfluP = 0;
    CalcLat = 0;
    transient.trackRGO = 0;
    transient.trackDeng = 0;
    gasLift.celulaG = 0;
    celula = 0;
    gasLift.celInter = 1e7;
    gasLift.dtInter = 0.;
    gasLift.velInter = 0.;

    trends.MatTrendG = 0;
    trends.resettrendg = 0;
    trends.ntrendg = 0;
    trends.ntrendgB = 0;
    trends.MatTrendTransG = 0;
    trends.resettrendtransg = 0;
    trends.ntrendtransg = 0;
    trends.ntrendtransgB = 0;
    trends.MatTrendP = 0;
    trends.resettrend = 0;
    trends.ntrend = 0;
    trends.ntrendB = 0;
    trends.MatTrendTransP = 0;
    trends.resettrendtrans = 0;
    trends.ntrendtrans = 0;
    trends.ntrendtransB = 0;

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

    gasLift.tempMedContDesc = 10.;
    gasLift.maxVecContDesc = 1000;
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

    transient.kontaGolfada = 1000;

    transient.mudaModoChk = 0;
    networkCoupling.mudaModoChkini = 0;

    transient.momentoDesesp = 0;

    transient.modeloCompleto = 1;

    transient.DpMaxMed = 1.;
    transient.DTMaxMed = 1.;

    chuteHol = -1.;

    buscaIni = 0;

    for (int i = 0; i < 10; i++) {
        transient.vRazMast0[i] = 0.;
        transient.vRazMast1[i] = 0.;
        transient.vRazMastCrit[i] = 0.5;
    }

    transient.kontaRenovaComp = 0;

    fluiRevRede = ProFlu();
    tempRev = 0.;
    revPerm = 0;
    tables.ntabDin = 0;

    transient.nCelulaPoisson2D = 0;
    trocaTermicaLenta = 0.01;

    semTermo = 0;

    monitConvPerm = 1000.;
    monitConvPermBase = 1.;

    transient.alteraTempo = 0;
}

SProd::SProd(string nomeArquivoEntrada, string nomeArquivoLog, tipoValidacaoJson_t validacaoJson,
             tipoSimulacao_t tipoSimulacao, varGlob1D *Vvg1dSP, int TD, int vbloq, int temporario, int reverso, double *compfonte,
             int *posicfonte, int nfontes, int redeperm) : arq(nomeArquivoEntrada, nomeArquivoLog, validacaoJson, tipoSimulacao, reverso, Vvg1dSP, redeperm),
                                                           flutG(arq.ncelg, arq.nvarprofg + 2 + 1 + 1 + 1 + 1 + 1), flut(arq.ncelp, arq.nvarprofp + 2 + 1 + 1 + 1 + 1),
                                                           gasLift{.matglobG = BandMtx<double>(3 * arq.ncelg, 5, 5), .termolivreG = Vcr<double>(3 * arq.ncelg)},
                                                           transient{.matglobP = BandMtx<double>(2 * arq.ncelp, 3, 2), .termolivreP = Vcr<double>(2 * arq.ncelp)} {
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

SProd::SProd() : arq(), flutG(1, 1 + 2 + 1 + 1 + 1 + 1), flut(1, 1 + 2 + 1 + 1 + 1),
                 gasLift{.matglobG = BandMtx<double>(3 * 1, 5, 5), .termolivreG = Vcr<double>(3 * 1)},
                 transient{.matglobP = BandMtx<double>(2 * 1, 3, 2), .termolivreP = Vcr<double>(2 * 1)} {
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

/// Gives back the memory of every array this object owns, so that a rebuild starts
/// with it free.
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
    acertaIndAcop.clear();
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
    flut = FullMtx<double>(arq.ncelp, arq.nvarprofp + 2 + 1 + 1 + 1 + 1);
    flutG = FullMtx<double>(arq.ncelg, arq.nvarprofg + 2 + 1 + 1 + 1 + 1 + 1);
    transient.matglobP = BandMtx<double>(2 * arq.ncelp, 3, 2);
    transient.termolivreP = Vcr<double>(2 * arq.ncelp);
    gasLift.matglobG = BandMtx<double>(3 * arq.ncelg, 5, 5);
    gasLift.termolivreG = Vcr<double>(3 * arq.ncelg);
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
    sisprod::gaslift::computeGasUnloadingHydrostatics(gasLiftStateOf(*this));
}

void SProd::HidroDescargaP() {
    sisprod::gaslift::computeProductionUnloadingHydrostatics(gasLiftStateOf(*this));
}

double SProd::areaValvCali(double PCal, double TCal, double PVO, double PT,
                           double dextern, double areagarg, double Rvalv, double Temp) {
    return sisprod::gaslift::calibratedValveArea(PCal, TCal, PVO, PT, dextern, areagarg, Rvalv, Temp);
}

void SProd::calctempGas(int i, double tempantiga, int modoPerm) {
    sisprod::thermal::computeGasTemperature(thermalStateOf(*this), i, tempantiga, modoPerm);
}

void SProd::tempDescarga(int i) {
    sisprod::thermal::computeDischargeTemperature(thermalStateOf(*this), i);
}

double SProd::TempDescGL(int igl) {
    return sisprod::thermal::computeGasLiftDischargeTemperature(thermalStateOf(*this), igl);
}

double SProd::BuscaPresInjDesc() {
    return sisprod::gaslift::searchUnloadingInjectionPressure(gasLiftStateOf(*this));
}

void SProd::subtempoGas() {
    sisprod::gaslift::advanceGasSubStep(gasLiftStateOf(*this));
}

void SProd::conectaColuna() {
    sisprod::gaslift::connectTubing(gasLiftStateOf(*this));
}

void SProd::calctemp(int i, double tempantiga, int modoPerm) {
    sisprod::thermal::computeTemperature(
        thermalStateOf(*this), i, tempantiga, modoPerm);
}

void SProd::renovaFonte(int ind) {
    sisprod::sources::renewSourceTerms(sourceStateOf(*this), ind);
}

void SProd::renovaalbetini() {
    sisprod::composition::storePreviousFractionsAndMovePigs(compositionStateOf(*this));
}

void SProd::renovaMasEsp() {
    sisprod::composition::cacheCellAndFaceDensities(compositionStateOf(*this));
}

void SProd::CalcC0Ud(int ind, double &c0, double &ud) {
    driftflux::coefficient::instantaneous(closureStateOf(*this), ind, c0, ud);
}

void SProd::CalcC0UdBuf(int ind, double &c0, double &ud) {
    driftflux::coefficient::buffered(closureStateOf(*this), ind, c0, ud);
}

void SProd::CalcC0UdIni(int ind, double &c0, double &ud) {
    driftflux::coefficient::initialization(closureStateOf(*this), ind, c0, ud);
}

void SProd::CalcC0UdIniBuf(int ind, double &c0, double &ud) {
    driftflux::coefficient::bufferedInitialization(closureStateOf(*this), ind, c0, ud);
}

void SProd::renova(int expli) {
    sisprod::transient::updateCells(transientStateOf(*this), expli);
}

void SProd::renovaBuffer() {
    sisprod::transient::updateBufferFromSolution(transientStateOf(*this));
}

void SProd::renovaBufferCego() {
    sisprod::transient::updateBufferFromCells(transientStateOf(*this));
}

void SProd::renovaTemp() {
    sisprod::thermal::updateDistributedMassTransfer(
        thermalStateOf(*this));
}

void SProd::avaliaParafina() {
    sisprod::composition::evaluateWaxDeposition(compositionStateOf(*this));
}

void SProd::renovaRGOdgYco2(ProFlu fluiRev) {
    sisprod::composition::transportBlackOilProperties(compositionStateOf(*this), fluiRev);
}

/*** change 4 ***/

void SProd::renovaFracMol2(ProFlu fluiRev) {
    sisprod::composition::transportPhaseMolarFractions(compositionStateOf(*this), fluiRev);
}

void SProd::renovaterm(int aflu) {
    sisprod::thermal::updateFlowPartitionTerms(thermalStateOf(*this), aflu);
}

void SProd::renovatermAfluFim() {
    sisprod::thermal::updateOutletFlowPartitionTerms(thermalStateOf(*this));
}

void SProd::renovatermColIni() {
    sisprod::thermal::updateInletFlowPartitionTerms(thermalStateOf(*this));
}

void SProd::calcCCpres(double titRev, double alfRev, double betRev) {
    sisprod::transient::applyOutletPressureCondition(transientStateOf(*this), titRev, alfRev, betRev);
}

void SProd::calcCCBuffer(double titRev, double alfRev, double betRev) {
    sisprod::transient::applyOutletBufferCondition(transientStateOf(*this), titRev, alfRev, betRev);
}

void SProd::determinaDT(int vexpli) {
    sisprod::transient::computeTimeStep(transientStateOf(*this), vexpli);
}

void SProd::atenuaDtMax() {
    sisprod::transient::dampMaximumTimeStep(transientStateOf(*this));
}

void SProd::avaliaVariaDpDt(double razMast, double razMast0, int vexpli) {
    sisprod::transient::evaluatePressureRateOfChange(transientStateOf(*this), razMast, razMast0, vexpli);
}

void SProd::aberturaVal0() {
    sisprod::transient::valveOpeningLow(transientStateOf(*this));
}
void SProd::aberturaVal1() {
    sisprod::transient::valveOpeningHigh(transientStateOf(*this));
}
void SProd::restringeDTporValv() {
    sisprod::transient::restrictTimeStepByValve(transientStateOf(*this));
}

void SProd::solveLinGas() {
    sisprod::gaslift::solveGasLine(gasLiftStateOf(*this));
}

void SProd::EvoluiFrac(double alfrev, double betrev, int ciclo) {
    sisprod::transient::evolveFractions(transientStateOf(*this), alfrev, betrev, ciclo);
}

void SProd::ReiniEvolFrac0() {
    sisprod::transient::restartFractionEvolutionInitial(transientStateOf(*this));
}

void SProd::ReiniEvolFrac() {
    sisprod::transient::restartFractionEvolution(transientStateOf(*this));
}

void SProd::AtualizaPig() {
    sisprod::transient::updatePig(transientStateOf(*this));
}

void SProd::SolveAcopPV(int vexpli, int ciclo) {
    sisprod::transient::solvePressureVolumeCoupling(transientStateOf(*this), vexpli, ciclo);
}

void SProd::marchaEnergTrans(int ciclo, int ciclomax) {
    sisprod::thermal::advanceTransientEnergy(thermalStateOf(*this), ciclo, ciclomax);
}

void SProd::atualizaMiniTab() {
    sisprod::transient::refreshFluidMiniTable(transientStateOf(*this));
}

void SProd::atualizaCC1() {
    sisprod::transient::refreshInletCondition(transientStateOf(*this));
}

/// Runs the hydrate-envelope solvers for the production and gas lines.
///
/// Both solvers are constructed from the whole SProd object (*this).
/// SolveTrans reaches it through a callback, so the hydrate phase runs first in
/// the step.
void SProd::solveHydrateEnvelopes() {
    if (arq.calculaEnvelope == 1 && (*vg1dSP).lixo5 <= arq.tfinal) { // *vg1dSP).lixo5>0 && //chris - hydrates
        FA_Hidrato solverHidrato(*this);
        solverHidrato.solverHidrato();
    }

    if (arq.lingas > 0 && arq.calculaEnvelope == 1 && (*vg1dSP).lixo5 <= arq.tfinal) { // *vg1dSP).lixo5>0 && //chris - hydrates
        FA_Hidrato_Servico solverHidratoG(*this);
        solverHidratoG.solverHidratoG();
    }
}

void SProd::SolveTrans(double titRev, double alfRev, double betRev, int nrede, ProFlu fluiRev) {
    sisprod::transient::solveTransientStep(transientSolveStateOf(*this), titRev, alfRev, betRev, nrede, fluiRev);
}

// The trend writers live in SisProdTrendOutput.cpp; Num4Main.cpp calls four
// of them through these.
void SProd::ImprimeTrendPCab(int i, int nrede) {
    trendoutput::writeProductionTrendHeader(trendStateOf(*this), i, nrede);
}
void SProd::ImprimeTrendP(int i, int nrede) {
    trendoutput::writeProductionTrendRows(trendStateOf(*this), i, nrede);
}
void SProd::ImprimeTrendGCab(int i, int nrede) {
    trendoutput::writeGasLineTrendHeader(trendStateOf(*this), i, nrede);
}
void SProd::ImprimeTrendG(int i, int nrede) {
    trendoutput::writeGasLineTrendRows(trendStateOf(*this), i, nrede);
}
void SProd::ImprimeTrendTransPCab(int i) {
    trendoutput::writeProductionCrossSectionTrendHeader(trendStateOf(*this), i);
}
void SProd::ImprimeTrendTransP(int i) {
    trendoutput::writeProductionCrossSectionTrendRows(trendStateOf(*this), i);
}
void SProd::ImprimeTrendTransGCab(int i) {
    trendoutput::writeGasLineCrossSectionTrendHeader(trendStateOf(*this), i);
}
void SProd::ImprimeTrendTransG(int i) {
    trendoutput::writeGasLineCrossSectionTrendRows(trendStateOf(*this), i);
}

double SProd::marchaProdPerm1(double pchute) {
    return sisprod::steady::marchProductionSteady(steadyStateOf(*this), pchute);
}

double SProd::buscaProdPfundoPerm(double chute, int kontaTenta) {
    return sisprod::steady::searchProductionBottomHolePressure(searchStateOf(*this), chute, kontaTenta);
}

double SProd::buscaProdPfundoPermRev(double chute) {
    return sisprod::steady::searchReverseProductionBottomHolePressure(searchStateOf(*this), chute);
}

double SProd::buscaProdPfundoPerm2(double chute, int kontaTenta) {
    return sisprod::steady::searchProductionBottomHolePressureSecondary(searchStateOf(*this), chute, kontaTenta);
}

double SProd::buscaProdPfundoPerm3(double pentrada) {
    return sisprod::steady::searchProductionBottomHolePressureTertiary(searchStateOf(*this), pentrada);
}

double SProd::marchaProdPresPres1(double mchute) {
    return sisprod::steady::marchProductionPressureToPressure(steadyStateOf(*this), mchute);
}

double SProd::buscaProdPresPresPerm(double chute, double maxvaz, int kontaiter) {
    return sisprod::steady::searchProductionPressureToPressure(searchStateOf(*this), chute, maxvaz, kontaiter);
}

double SProd::buscaProdPresPresPermRev(double chute, double maxvaz, int kontaiter) {
    return sisprod::steady::searchReverseProductionPressureToPressure(searchStateOf(*this), chute, maxvaz, kontaiter);
}

double SProd::buscaProdPresPresPerm2(double chute, double maxvaz) {
    return sisprod::steady::searchProductionPressureToPressureSecondary(searchStateOf(*this), chute, maxvaz);
}

double SProd::buscaProdPresPresPerm3(double chute, double maxvaz) {
    return sisprod::steady::searchProductionPressureToPressureTertiary(searchStateOf(*this), chute, maxvaz);
}

double SProd::buscaGasPresPerm2() {
    return sisprod::steady::searchGasPressureSteadySecondary(searchStateOf(*this));
}

double SProd::buscaGasPresPerm3() {
    return sisprod::steady::searchGasPressureSteadyTertiary(searchStateOf(*this));
}

void SProd::corrDeng(int i) {
    sisprod::steady::correctGasSpecificGravity(steadyStateOf(*this), i);
}















void SProd::CalcC0UdPerm(int ind, double &c0, double &ud) {
    driftflux::coefficient::steadyState(closureStateOf(*this), ind, c0, ud);
}

void SProd::RenovaTempPerm(int i, int RK) {
    sisprod::thermal::advanceSteadyTemperature(thermalStateOf(*this), i, RK);
}

void SProd::RenovaTempPermRev(int i, int RK) {
    sisprod::thermal::advanceReverseSteadyTemperature(thermalStateOf(*this), i, RK);
}

void SProd::atualizaPeriTempProd(int i) {
    sisprod::thermal::updateProductionTemperaturePeriphery(thermalStateOf(*this), i);
}

void SProd::calcTempFim() {
    sisprod::thermal::computeOutletTemperature(thermalStateOf(*this));
}

double SProd::delpGasPerm(int i) {
    return sisprod::gaslift::steadyGasPressureDrop(gasLiftStateOf(*this), i);
}

double SProd::delpInjPerm(int i) {
    return sisprod::gaslift::steadyInjectionPressureDrop(gasLiftStateOf(*this), i);
}

void SProd::RenovaPresGasPerm(int i) {
    sisprod::gaslift::updateSteadyGasPressure(gasLiftStateOf(*this), i);
}

void SProd::calcVazGasPerm(int i) {
    sisprod::gaslift::computeSteadyGasFlowRate(gasLiftStateOf(*this), i);
}

void SProd::IniciaVazValvGasPerm(int i) {
    sisprod::gaslift::initializeSteadyValveGasFlowRate(gasLiftStateOf(*this), i);
}

void SProd::RenovaTempGasPerm(int i) {
    sisprod::gaslift::updateSteadyGasTemperature(gasLiftStateOf(*this), i);
}

void SProd::conectaColunaPerm() {
    sisprod::gaslift::connectTubingSteady(gasLiftStateOf(*this));
}

void SProd::IniciaconectaColunaPerm() {
    sisprod::gaslift::initializeTubingConnectionSteady(gasLiftStateOf(*this));
}

double SProd::buscaInjPfundoPerm1(double chute) {
    return sisprod::steady::searchInjectionBottomHolePressure1(searchStateOf(*this), chute);
}

double SProd::buscaInjPfundoPerm2(double chute) {
    return sisprod::steady::searchInjectionBottomHolePressure2(searchStateOf(*this), chute);
}

double SProd::buscaInjPfundoPerm3(double chute) {
    return sisprod::steady::searchInjectionBottomHolePressure3(searchStateOf(*this), chute);
}

double SProd::buscaInjPfundoPerm4() {
    return sisprod::steady::searchInjectionBottomHolePressure4(searchStateOf(*this));
}

double SProd::buscaInjPfundoPerm5(double chute) {
    return sisprod::steady::searchInjectionBottomHolePressure5(searchStateOf(*this), chute);
}

double SProd::hidroreverso(double hol, double vaz, double vazG) {
    return sisprod::steady::reverseHydrostatic(steadyStateOf(*this), hol, vaz, vazG);
}

double SProd::hidroreversoInj(double hol, double vaz) {
    return sisprod::steady::reverseInjectionHydrostatic(steadyStateOf(*this), hol, vaz);
}

double SProd::hidroTramoSecundario(double titulo) {
    return sisprod::steady::secondaryBranchHydrostatic(steadyStateOf(*this), titulo);
}

