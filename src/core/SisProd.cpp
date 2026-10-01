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
    celulaG = 0;
    celula = 0;
    celInter = 1e7;
    dtInter = 0.;
    velInter = 0.;

    MatTrendG = 0;
    trends.resettrendg = 0;
    trends.ntrendg = 0;
    trends.ntrendgB = 0;
    trends.MatTrendTransG = 0;
    trends.resettrendtransg = 0;
    trends.ntrendtransg = 0;
    trends.ntrendtransgB = 0;
    MatTrendP = 0;
    trends.resettrend = 0;
    ntrend = 0;
    trends.ntrendB = 0;
    trends.MatTrendTransP = 0;
    trends.resettrendtrans = 0;
    trends.ntrendtrans = 0;
    trends.ntrendtransB = 0;

    fontemassPRBuf = 0.;
    fontemassCRBuf = 0.;
    fontemassGRBuf = 0.;

    presE = -1;
    tempE = -1;
    titE = -1;
    betaE = -1;
    alfE = -1;
    presEini = -1;
    tempEini = -1;
    titEini = -1;
    betaEini = -1;
    alfEini = -1;

    gasLift.tempMedContDesc = 10.;
    gasLift.maxVecContDesc = 1000;
    gasLift.vazmedDesc = 0;
    gasLift.tempmedDEsc = 0;

    tGSup = 0.;
    tGSupIni = 0.;

    transient.dtCFLMed = 1.;
    transient.dtSimMed = 1.;
    transient.restriDt = 0;
    transient.kontarestriDt = 0;
    dtauxCFL = 0.;
    dtauxFinal = 0.;

    kimpT = 0.;

    transient.kontaGolfada = 1000;

    mudaModoChk = 0;
    mudaModoChkini = 0;

    transient.momentoDesesp = 0;

    modeloCompleto = 1;

    transient.DpMaxMed = 1.;
    transient.DTMaxMed = 1.;

    chuteHol = -1.;

    buscaIni = 0;

    for (int i = 0; i < 10; i++) {
        vRazMast0[i] = 0.;
        vRazMast1[i] = 0.;
        vRazMastCrit[i] = 0.5;
    }

    transient.kontaRenovaComp = 0;

    fluiRevRede = ProFlu();
    tempRev = 0.;
    revPerm = 0;
    ntabDin = 0;

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
                                                           transient{.matglobP = BandMtx<double>(2 * arq.ncelp, 3, 2)}, termolivreP(2 * arq.ncelp) {
    resolveDriftSelectors();
    resetRunState();
    celInterIni = celInter;
    dtInterIni = dtInter;
    velInterIni = velInter;
    redeTemporario = temporario;
    betaRevini = 0;
    bloq = vbloq;
    vg1dSP = Vvg1dSP;
    if (TD >= 0)
        arq.tabelaDinamica = TD;
    montasistema(compfonte, posicfonte, nfontes);
}

SProd::SProd() : arq(), flutG(1, 1 + 2 + 1 + 1 + 1 + 1), flut(1, 1 + 2 + 1 + 1 + 1),
                 gasLift{.matglobG = BandMtx<double>(3 * 1, 5, 5), .termolivreG = Vcr<double>(3 * 1)},
                 transient{.matglobP = BandMtx<double>(2 * 1, 3, 2)}, termolivreP(2 * 1) {
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
    releaseVector(gasCells);
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
        tabDin.clear();
    }
}

SProd &SProd::operator=(const SProd &sp) {
    releaseOwnedStorage();

    arq = sp.arq;
    resolveDriftSelectors();
    flut = sp.flut;
    flutG = sp.flutG;
    transient.matglobP = sp.transient.matglobP;
    termolivreP = sp.termolivreP;
    gasLift.matglobG = sp.gasLift.matglobG;
    gasLift.termolivreG = sp.gasLift.termolivreG;
    vg1dSP = sp.vg1dSP;
    resetRunState();
    celInterIni = celInter;
    dtInterIni = dtInter;
    velInterIni = velInter;

    noextremo = sp.noextremo;
    noinicial = sp.noinicial;
    derivaAnel = sp.derivaAnel;

    betaRevini = sp.betaRevini;
    titRevini = sp.titRevini;

    bloq = sp.bloq;

    dtCicMin = sp.dtCicMin;
    redeParalelaCCsecundario = sp.redeParalelaCCsecundario;
    redeParalelaP = sp.redeParalelaP;
    redeParalelaS = sp.redeParalelaS;

    presVet.clear();
    jVet.clear();
    alfVet.clear();
    tVet.clear();
    transient.dtSim.clear();
    transient.dtCFL.clear();
    gasLift.vazmaxMedDesc.clear();
    gasLift.dtDesc.clear();
    transient.taxaDpMax.clear();
    transient.taxaDTMax.clear();
    tabDin.clear();
    acertaIndAcop.clear();
    transient.indCelPoisson2D.clear();
    indFonteRedeParalelaIni.clear();
    fonteMpRedeParalelaIni.clear();
    fonteMcRedeParalelaIni.clear();
    fonteMgRedeParalelaIni.clear();

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
    termolivreP = Vcr<double>(2 * arq.ncelp);
    gasLift.matglobG = BandMtx<double>(3 * arq.ncelg, 5, 5);
    gasLift.termolivreG = Vcr<double>(3 * arq.ncelg);
    vg1dSP = arq.vg1dSP;
    resetRunState();
    celInterIni = celInter;
    dtInterIni = dtInter;
    velInterIni = velInter;

    noextremo = carried.noextremo;
    noinicial = carried.noinicial;
    derivaAnel = carried.derivaAnel;

    betaRevini = carried.betaRevini;
    titRevini = carried.titRevini;

    bloq = carried.bloq;

    dtCicMin = carried.dtCicMin;

    redeParalelaCCsecundario = -1;
    redeParalelaP = -1;
    redeParalelaS = -1;

    montasistema();
}

SProd::CarriedState SProd::carriedState() const {
    return {.noextremo = noextremo,
            .noinicial = noinicial,
            .derivaAnel = derivaAnel,
            .bloq = bloq,
            .betaRevini = betaRevini,
            .titRevini = titRevini,
            .dtCicMin = dtCicMin};
}

void SProd::HidroDescargaG() {
    sisprod::gaslift::computeGasUnloadingHydrostatics(gasLiftStateOf(*this));
}

void SProd::HidroDescargaP() {
    sisprod::gaslift::computeProductionUnloadingHydrostatics(gasLiftStateOf(*this));
}

void SProd::renovaGas() {
    sisprod::gaslift::updateGasLine(gasLiftStateOf(*this));
}

void SProd::renovaGasBuf() {
    sisprod::gaslift::updateBufferedGasLine(gasLiftStateOf(*this));
}

double SProd::areaValvCali(double PCal, double TCal, double PVO, double PT,
                           double dextern, double areagarg, double Rvalv, double Temp) {
    return sisprod::gaslift::calibratedValveArea(PCal, TCal, PVO, PT, dextern, areagarg, Rvalv, Temp);
}

void SProd::calctempGas(int i, double tempantiga, int modoPerm) {
    sisprod::thermal::computeGasTemperature(thermalStateOf(*this), i, tempantiga, modoPerm);
}

void SProd::resolveDescarga() {
    sisprod::gaslift::solveUnloading(gasLiftStateOf(*this));
}

void SProd::tempDescarga(int i) {
    sisprod::thermal::computeDischargeTemperature(thermalStateOf(*this), i);
}

void SProd::avancInter() {
    sisprod::gaslift::advanceInterface(gasLiftStateOf(*this));
}

double SProd::TempDescGL(int igl) {
    return sisprod::thermal::computeGasLiftDischargeTemperature(thermalStateOf(*this), igl);
}

void SProd::ValvGasTrans() {
    sisprod::gaslift::updateTransientGasValves(gasLiftStateOf(*this));
}

double SProd::prescordesc(double vazmax, int ivalv, double fator, int sinal) {
    return sisprod::gaslift::unloadingPressureCorrection(gasLiftStateOf(*this), vazmax, ivalv, fator, sinal);
}

double SProd::CalcPresValvDesc(double vazGarg, int ivalv) {
    return sisprod::gaslift::computeUnloadingValvePressure(gasLiftStateOf(*this), vazGarg, ivalv);
}

double SProd::BuscaPresInjDesc() {
    return sisprod::gaslift::searchUnloadingInjectionPressure(gasLiftStateOf(*this));
}

void SProd::subtempoGas() {
    sisprod::gaslift::advanceGasSubStep(gasLiftStateOf(*this));
}

void SProd::subtempoGasBuf() {
    sisprod::gaslift::advanceBufferedGasSubStep(gasLiftStateOf(*this));
}

void SProd::conectaColuna() {
    sisprod::gaslift::connectTubing(gasLiftStateOf(*this));
}

double SProd::interpolaHLatente(double pres, double temp) {
    return sisprod::thermal::interpolateLatentHeat(
        thermalStateOf(*this), pres, temp);
}

void SProd::calctemp(int i, double tempantiga, int modoPerm) {
    sisprod::thermal::computeTemperature(
        thermalStateOf(*this), i, tempantiga, modoPerm);
}

double SProd::calcHmix(int i) {
    return sisprod::thermal::computeMixtureEnthalpy(
        thermalStateOf(*this), i);
}

double SProd::energmix(int i, int jp0, int jt, double razp) {
    return sisprod::thermal::interpolateMixtureEnergy(
        thermalStateOf(*this), i, jp0, jt, razp);
}

void SProd::calcTempEntalp(int i) {
    sisprod::thermal::updateTemperatureFromEnthalpy(
        thermalStateOf(*this), i);
}

void SProd::calcTransMassTermo(int i) {
    sisprod::thermal::computeThermalMassTransfer(
        thermalStateOf(*this), i);
}

void SProd::FonteValv(int ind) {
    sisprod::sources::addMasterValveFlow(sourceStateOf(*this), ind);
}

void SProd::salvaFonte() {
    sisprod::transient::storePreviousSources(transientStateOf(*this));
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

void SProd::correcHidroFric(int i, double &hidro, double &fric) {

    double dx = 0.5 * celula[i].dx;
    double dia = celula[i].duto.a;
    double area = 0.25 * M_PI * dia * dia;
    double si = celula[i].duto.peri;
    double alfmed = celula[i].alf;
    double rhog = celula[i].flui.MasEspGas(celula[i].pres, celula[i].temp);
    double rhol = (1 - celula[i].bet) * celula[i].flui.MasEspLiq(celula[i].pres, celula[i].temp) + celula[i].bet * celula[i].fluicol.MasEspFlu(celula[i].pres, celula[i].temp);
    double ugsmed = (celula[i].MC - celula[i].Mliqini) / (area * rhog);
    double ulsmed = celula[i].Mliqini / (area * rhol);
    double j = ugsmed + ulsmed;

    double rhomix = alfmed * rhog + (1 - alfmed) * rhol;
    double viscmix = alfmed * celula[i].flui.ViscGas(celula[i].pres, celula[i].temp) + (1 - alfmed) * ((1 - celula[i].bet) * celula[i].flui.ViscOleo(celula[i].pres, celula[i].temp) + celula[i].bet * celula[i].fluicol.VisFlu(celula[i].pres, celula[i].temp));

    double re1;
    if (celula[i].duto.revest == 0)
        re1 = celula[i].Rey(celula[i].duto.a, j, rhomix, viscmix);
    else {
        double dhid = 4 * area / si;
        re1 = celula[i].Rey(dhid, j, rhomix, viscmix);
    }
    double f1 = celula[i].fric(re1, celula[i].duto.rug / dia);
    fric = (1 - arq.MedSimpPresFront) * 0.5 * f1 * rhomix * (fabs(j) * j) * si * dx / area;
    hidro = (1 - arq.MedSimpPresFront) * kGravity * sin(celula[i].duto.teta) * rhomix * dx;
}

void SProd::auxMiniTab(ProFlu &flui) {
    sisprod::transient::fillFluidMiniTable(transientStateOf(*this), flui);
}

void SProd::geraMiniTabFlu() {
    sisprod::transient::generateFluidMiniTables(transientStateOf(*this));
}

void SProd::renova(int expli) {
    sisprod::transient::updateCells(transientStateOf(*this), expli);
}

void SProd::renovaVaz() {
    sisprod::transient::updateFlowRates(transientStateOf(*this));
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

void SProd::determinaDTExpli() {
    sisprod::transient::computeExplicitTimeStep(transientStateOf(*this));
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

void SProd::SubReiniEvolFrac() {
    sisprod::transient::restartFractionEvolutionSub(transientStateOf(*this));
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

void SProd::prepDifusCalorND(int i) {
    sisprod::thermal::prepareNonDimensionalHeatDiffusion(thermalStateOf(*this), i);
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

double SProd::marchaProdPerm1Rev(double pchute) {
    return sisprod::steady::marchReverseProductionSteady(steadyStateOf(*this), pchute);
}

double SProd::marchaProdPerm2(double pchute) {
    return sisprod::steady::marchProductionSteadySecondary(steadyStateOf(*this), pchute);
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

double SProd::marchaProdPresPres1Rev(double mchute) {
    return sisprod::steady::marchReverseProductionPressureToPressure(steadyStateOf(*this), mchute);
}

double SProd::buscaProdPresPresPerm(double chute, double maxvaz, int kontaiter) {
    return sisprod::steady::searchProductionPressureToPressure(searchStateOf(*this), chute, maxvaz, kontaiter);
}

double SProd::buscaProdPresPresPermRev(double chute, double maxvaz, int kontaiter) {
    return sisprod::steady::searchReverseProductionPressureToPressure(searchStateOf(*this), chute, maxvaz, kontaiter);
}

double SProd::marchaProdPresPres2(double mchute) {
    return sisprod::steady::marchProductionPressureToPressureSecondary(steadyStateOf(*this), mchute);
}

double SProd::buscaProdPresPresPerm2(double chute, double maxvaz) {
    return sisprod::steady::searchProductionPressureToPressureSecondary(searchStateOf(*this), chute, maxvaz);
}

double SProd::marchaProdPresPres3(double mchute) {
    return sisprod::steady::marchProductionPressureToPressureTertiary(steadyStateOf(*this), mchute);
}

double SProd::buscaProdPresPresPerm3(double chute, double maxvaz) {
    return sisprod::steady::searchProductionPressureToPressureTertiary(searchStateOf(*this), chute, maxvaz);
}

double SProd::marchaGasPerm1(double chutemass) {
    return sisprod::steady::marchGasSteady(steadyStateOf(*this), chutemass);
}

double SProd::buscaGasPresPerm2() {
    return sisprod::steady::searchGasPressureSteadySecondary(searchStateOf(*this));
}

double SProd::buscaGasPresPerm3() {
    return sisprod::steady::searchGasPressureSteadyTertiary(searchStateOf(*this));
}

double SProd::marchaGasPerm2(double pchute, double chutemass) {
    return sisprod::steady::marchGasSteadySecondary(steadyStateOf(*this), pchute, chutemass);
}

double SProd::marchaGasPerm3(double pchute) {
    return sisprod::steady::marchGasSteadyTertiary(steadyStateOf(*this), pchute);
}

void SProd::RenovaPresPermMon(int i, int RK) {
    sisprod::steady::advanceUpstreamSteadyPressure(steadyStateOf(*this), i, RK);
}

double SProd::RenovaPresPermNcel() {
    return sisprod::steady::steadyPressureAtLastCell(steadyStateOf(*this));
}

double SProd::calcDpArea(int i, double rhomix, double rey, double jmix) {
    return sisprod::steady::areaChangePressureDrop(steadyStateOf(*this), i, rhomix, rey, jmix);
}

void SProd::RenovaPresPermJus(int i, int RK) {
    sisprod::steady::advanceDownstreamSteadyPressure(steadyStateOf(*this), i, RK);
}

void SProd::corrDeng(int i) {
    sisprod::steady::correctGasSpecificGravity(steadyStateOf(*this), i);
}















void SProd::RenovaMassPerm(int i) {
    sisprod::steady::advanceSteadyMass(steadyStateOf(*this), i);
}
void SProd::RenovaMassPermRev(int i) {
    sisprod::steady::advanceReverseSteadyMass(steadyStateOf(*this), i);
}

void SProd::RenovaMassPermComp(int i) {
    sisprod::steady::advanceCompositionalSteadyMass(steadyStateOf(*this), i);
}

void SProd::RenovaMassPermCompRev(int i) {
    sisprod::steady::advanceReverseCompositionalSteadyMass(steadyStateOf(*this), i);
}

void SProd::CalcC0UdPerm(int ind, double &c0, double &ud) {
    driftflux::coefficient::steadyState(closureStateOf(*this), ind, c0, ud);
}

void SProd::RenovaTransMassPerm(int i) {
    sisprod::steady::advanceSteadyMassTransfer(steadyStateOf(*this), i);
}

void SProd::RenovaTransMassPermGas(int i) {
    sisprod::steady::advanceSteadyGasMassTransfer(steadyStateOf(*this), i);
}

void SProd::RenovaTempPerm(int i, int RK) {
    sisprod::thermal::advanceSteadyTemperature(thermalStateOf(*this), i, RK);
}

void SProd::RenovaTempPermRev(int i, int RK) {
    sisprod::thermal::advanceReverseSteadyTemperature(thermalStateOf(*this), i, RK);
}

void SProd::atualizaPeriPmonProd(int i) {
    sisprod::steady::refreshUpstreamProductionPeriphery(steadyStateOf(*this), i);
}
void SProd::atualizaPeriPjusProd(int i) {
    sisprod::steady::refreshDownstreamProductionPeriphery(steadyStateOf(*this), i);
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

void SProd::atualizaProp() {
    sisprod::steady::refreshProperties(steadyStateOf(*this));
}

void SProd::atualizaVelTermPerm() {
    sisprod::steady::refreshSteadyThermalVelocities(steadyStateOf(*this));
}

void SProd::calcDTPseudoTrans() {
    sisprod::steady::computePseudoTransientTimeStep(steadyStateOf(*this));
}
double SProd::marchaInjPerm1(double chute) {
    return sisprod::steady::marchInjectionSteady(steadyStateOf(*this), chute);
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

double SProd::multMarcha(double chute, int prod, int tipoCC) {
    return sisprod::steady::dispatchMarch(searchStateOf(*this), chute, prod, tipoCC);
}
double SProd::zriddr(double x1, double x2, int prod, int tipoCC) {
    return sisprod::steady::solveSteadyRoot(searchStateOf(*this), x1, x2, prod, tipoCC);
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

void SProd::hidroLinServ() {
    sisprod::steady::gasLineHydrostatic(steadyStateOf(*this));
}

double SProd::buscaTramoSecVazPerm(double pPartida, int indPartida) {
    return sisprod::steady::searchSecondaryBranchFlowRate(searchStateOf(*this), pPartida, indPartida);
}
