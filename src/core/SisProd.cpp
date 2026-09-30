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
using sisprod::adapters::trendStateOf;
using enum sisprod::AccessoryKind;
using sisprod::kAirDensityAtStandardConditions;
using sisprod::kBarrelPerCubicMetre;
using sisprod::kCubicFootPerCubicMetre;
using sisprod::kGravity;
using sisprod::kPascalPerKgfPerCm2Coarse;
using sisprod::kSecondsPerDay;

void SProd::resolveDriftSelectors() {
    driftSelectors = {arq.CorreDisper, arq.CorreAnular, arq.CorreEstrat};
}

/// The run state every construction and reassignment starts from, before each
/// path sets what is its own and montasistema builds the system.
void SProd::resetRunState() {
    zdranP = 0;
    dzdpP = 0;
    dzdtP = 0;
    cpg = 0;
    cpl = 0;
    drholdT = 0;
    npontos = 0;
    nfluP = 0;
    ModelCp = 0;
    Modeljtl = 0;
    CalcLat = 0;
    trackRGO = 0;
    trackDeng = 0;
    ninjgas = 0;
    lingas = 0;
    chokeVGL = 0;
    posicVGLP = 0;
    posicVGLG = 0;
    receb = 0;
    fechaM1 = 0;
    abreM1 = 0;
    celulaG = 0;
    celula = 0;
    celInter = 1e7;
    dtInter = 0.;
    velInter = 0.;

    ncelperftransg = 0;
    TrendLengthG = 0;
    MatTrendG = 0;
    resettrendg = 0;
    ntrendg = 0;
    ntrendgB = 0;
    TrendLengthTransG = 0;
    MatTrendTransG = 0;
    resettrendtransg = 0;
    ntrendtransg = 0;
    ntrendtransgB = 0;
    ncelperftransp = 0;
    TrendLengthP = 0;
    MatTrendP = 0;
    resettrend = 0;
    ntrend = 0;
    ntrendB = 0;
    TrendLengthTransP = 0;
    MatTrendTransP = 0;
    resettrendtrans = 0;
    ntrendtrans = 0;
    ntrendtransB = 0;
    LerPB = 0;
    PBPVTSim = 0;
    TBPVTSim = 0;

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

    tempMedContDesc = 10.;
    maxVecContDesc = 1000;
    vazmedDesc = 0;
    tempmedDEsc = 0;

    tGSup = 0.;
    tGSupIni = 0.;

    dtCFLMed = 1.;
    dtSimMed = 1.;
    restriDt = 0;
    kontarestriDt = 0;
    dtauxCFL = 0.;
    dtauxFinal = 0.;

    kimpT = 0.;

    kontaGolfada = 1000;

    mudaModoChk = 0;
    mudaModoChkini = 0;

    momentoDesesp = 0;

    modeloCompleto = 1;
    modeloCompleto0 = 1;
    kontaMudaModelo = 0;
    kontarestriSegrega = 0;

    DpMaxMed = 1.;
    DTMaxMed = 1.;

    chuteHol = -1.;

    buscaIni = 0;

    for (int i = 0; i < 10; i++) {
        vRazMast0[i] = 0.;
        vRazMast1[i] = 0.;
        vRazMastCrit[i] = 0.5;
    }

    kontaRenovaComp = 0;

    fluiRevRede = ProFlu();
    tempRev = 0.;
    revPerm = 0;
    ntabDin = 0;

    nCelulaPoisson2D = 0;
    trocaTermicaLenta = 0.01;

    semTermo = 0;

    monitConvPerm = 1000.;
    monitConvPermBase = 1.;

    alteraTempo = 0;
}

SProd::SProd(string nomeArquivoEntrada, string nomeArquivoLog, tipoValidacaoJson_t validacaoJson,
             tipoSimulacao_t tipoSimulacao, varGlob1D *Vvg1dSP, int TD, int vbloq, int temporario, int reverso, double *compfonte,
             int *posicfonte, int nfontes, int redeperm) : arq(nomeArquivoEntrada, nomeArquivoLog, validacaoJson, tipoSimulacao, reverso, Vvg1dSP, redeperm),
                                                           flutG(arq.ncelg, arq.nvarprofg + 2 + 1 + 1 + 1 + 1 + 1), flut(arq.ncelp, arq.nvarprofp + 2 + 1 + 1 + 1 + 1),
                                                           matglobG(3 * arq.ncelg, 5, 5), termolivreG(3 * arq.ncelg),
                                                           matglobP(2 * arq.ncelp, 3, 2), termolivreP(2 * arq.ncelp) {
    resolveDriftSelectors();
    resetRunState();
    celInterIni = celInter;
    dtInterIni = dtInter;
    velInterIni = velInter;
    RSLivia = 0;
    lerRS = 0;
    noextremo = 1;
    noinicial = 1;
    derivaAnel = -1;
    redeTemporario = temporario;

    betaRev = 0;
    betaRevini = 0;
    titRev = 1.;
    titRevini = 1.;
    dtCFLTotal = 0.;
    dtSimTotal = 0.;

    bloq = vbloq;

    vg1dSP = Vvg1dSP;

    if (TD >= 0)
        arq.tabelaDinamica = TD;

    redeParalelaCCsecundario = -1;
    redeParalelaP = -1;
    redeParalelaS = -1;
    montasistema(compfonte, posicfonte, nfontes);
}

SProd::SProd() : arq(), flutG(1, 1 + 2 + 1 + 1 + 1 + 1), flut(1, 1 + 2 + 1 + 1 + 1),
                 matglobG(3 * 1, 5, 5), termolivreG(3 * 1),
                 matglobP(2 * 1, 3, 2), termolivreP(2 * 1) {
    resolveDriftSelectors();
    resetRunState();
    tfinal = 0;
    dtini = 0;
    contaLog = 0;

    menorDx = 0.;
    iterperm = 0.;
    kSP = 0.;
    KontaImprime = 0.;
    indevento = 0.;
    modoPerm = 0.;
    ktMedMov = 0.;
    pTotal = 0.;
    jTotal = 0.;

    alfTotal = 0.;
    dt = 0.;
    nabreM1 = 0;
    nfechaM1 = 0;
    HLat = 0;

    celInterIni = 0.;
    dtInterIni = 0.;
    velInterIni = 0.;
    injPoc = 0;

    indTramo = -1;
    ncel = 0;
    reinicia = 0;
    presfim = 0;
    presfimini = 0;

    pGSup = 0;
    pGSupIni = 0.;
    temperatura = 0;

    masSup = 0;
    tempSup = 0;

    ncelGas = 0;
    presiniG = 0;
    tempiniG = 0;
    massfonte = 0;

    mult = 0;
    presMedMov = 0;
    jMedMov = 0;
    alfMedMov = 0;
    tMedMov = 0;

    aberto = 0;
    abertoini = 0;
    tempoaberto = 0;
    tempoabertoini = 0;
    EstadoMaster1 = 0;
    contaMaster1 = 0;
    masChkSup = 0;
    masChkSupini = 0;
    TransMassModel = 0;
    indpigP = 0;
    indpigPini = indpigP;
    npig = 0;

    AnulaColunaIni = 0;
    AnulaColunaFim = 0;
    ColunaAnulaIni = 0;
    ColunaAnulaFim = 0;
    verificaAcop = 0;
    verificaAcopRedeP = 0;
    verificaAcopRedeS = 0;
    SecPrimIniRedeP = 0;
    SecPrimFimRedeP = 0;
    PrimSecIniRedeP = 0;
    PrimSecFimRedeP = 0;
    kontaTempoProf = 0;
    //kontaTempoCelUni = 0;
    kontaTempoProfG = 0;
    kontaTempoTransProf = 0;
    kontaTempoTransProfG = 0;
    RSLivia = 0;
    lerRS = 0;

    noextremo = 1;
    noinicial = 1;
    derivaAnel = -1;

    titRev = 1.;
    titRevini = 1.;
    betaRev = 0.;
    betaRevini = 1.;
    redeTemporario = 0;
    dtCFLTotal = 0.;
    dtSimTotal = 0.;

    bloq = 0;

    vg1dSP = 0;
    dtCicMin = dt;

    redeParalelaCCsecundario = -1;
    redeParalelaP = -1;
    redeParalelaS = -1;
}

namespace {

/// Frees one set of trend matrices and its bookkeeping.
void releaseTrendSet(int count, double ***matrices, int *length, double *resetTimers, int *counts,
                     int *bufferedCounts) {
    for (int i = 0; i < count && matrices && length; i++) {
        if (matrices[i]) {
            for (int j = 0; j < length[i]; j++)
                delete[] matrices[i][j];
            delete[] matrices[i];
        }
    }
    delete[] matrices;
    delete[] length;
    delete[] resetTimers;
    delete[] counts;
    delete[] bufferedCounts;
}

}  // namespace

/// Frees every array this object owns, reading its current sizes and switches.
/// An array added to the construction must be released here too.
void SProd::releaseOwnedStorage() {
    if (arq.lingas > 0)
        delete[] celulaG;
    if (chokeVGL!=0 && arq.lingas > 0)
        delete[] chokeVGL;
    if (posicVGLP!=0 && arq.lingas > 0)
        delete[] posicVGLP;
    if (posicVGLG!=0 && arq.lingas > 0)
        delete[] posicVGLG;
    if (nabreM1 > 0)
        delete[] abreM1;
    if (nfechaM1 > 0)
        delete[] fechaM1;

    if (arq.nperfistransp > 0)
        delete[] ncelperftransp;
    if (arq.nperfistransg > 0 && arq.lingas > 0)
        delete[] ncelperftransg;

    if (arq.ntendp > 0 && redeTemporario == 0)
        releaseTrendSet(arq.ntendp, MatTrendP, TrendLengthP, resettrend, ntrend, ntrendB);
    if (arq.ntendg > 0 && arq.lingas > 0 && redeTemporario == 0)
        releaseTrendSet(arq.ntendg, MatTrendG, TrendLengthG, resettrendg, ntrendg, ntrendgB);
    if (arq.ntendtransp > 0 && redeTemporario == 0)
        releaseTrendSet(arq.ntendtransp, MatTrendTransP, TrendLengthTransP, resettrendtrans, ntrendtrans,
                        ntrendtransB);
    if (arq.ntendtransg > 0 && redeTemporario == 0)
        releaseTrendSet(arq.ntendtransg, MatTrendTransG, TrendLengthTransG, resettrendtransg, ntrendtransg,
                        ntrendtransgB);

    int ndiv = arq.tabent.npont - 1;
    if (CalcLat > 0 && arq.flashCompleto == 0) {
        for (int i = 0; i < ndiv + 2; i++)
            delete[] HLat[i];
        delete[] HLat;
    }
    if (LerPB > 0) {
        delete[] PBPVTSim;
        delete[] TBPVTSim;
        if (lerRS > 0) {
            for (int i = 0; i < ndiv + 2; i++)
                delete[] RSLivia[i];
            delete[] RSLivia;
        }
    }

    if (ncel > 0)
        delete[] celula;
    if (npig > 0)
        delete[] receb;

    if (arq.tabelaDinamica == 1) {
        tabDin.clear();
    }
}

SProd::~SProd() {
    releaseOwnedStorage();
}

SProd &SProd::operator=(const SProd &sp) {
    releaseOwnedStorage();

    arq = sp.arq;
    resolveDriftSelectors();
    flut = sp.flut;
    flutG = sp.flutG;
    matglobP = sp.matglobP;
    termolivreP = sp.termolivreP;
    matglobG = sp.matglobG;
    termolivreG = sp.termolivreG;
    vg1dSP = sp.vg1dSP;
    resetRunState();
    celInterIni = celInter;
    dtInterIni = dtInter;
    velInterIni = velInter;

    noextremo = sp.noextremo;
    noinicial = sp.noinicial;
    derivaAnel = sp.derivaAnel;

    betaRev = sp.betaRev;
    betaRevini = sp.betaRevini;
    titRev = sp.titRev;
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
    dtSim.clear();
    dtCFL.clear();
    vazmaxMedDesc.clear();
    dtDesc.clear();
    taxaDpMax.clear();
    taxaDTMax.clear();
    tabDin.clear();
    acertaIndAcop.clear();
    indCelPoisson2D.clear();
    indFonteRedeParalelaIni.clear();
    fonteMpRedeParalelaIni.clear();
    fonteMcRedeParalelaIni.clear();
    fonteMgRedeParalelaIni.clear();

    montasistema();

    return *this;
}

void SProd::copiaSemJson(Ler &sp, int vnoextremo, int vnoinicial, int vderivaAnel, int vbloq, double vbetaRev,
                         double vbetaRevini, double vtitRev, double vtitRevini, double vdtCicMin) {
    releaseOwnedStorage();

    arq.copiaSemJson(sp);
    resolveDriftSelectors();
    flut = FullMtx<double>(arq.ncelp, arq.nvarprofp + 2 + 1 + 1 + 1 + 1);
    flutG = FullMtx<double>(arq.ncelg, arq.nvarprofg + 2 + 1 + 1 + 1 + 1 + 1);
    matglobP = BandMtx<double>(2 * arq.ncelp, 3, 2);
    termolivreP = Vcr<double>(2 * arq.ncelp);
    matglobG = BandMtx<double>(3 * arq.ncelg, 5, 5);
    termolivreG = Vcr<double>(3 * arq.ncelg);
    vg1dSP = arq.vg1dSP;
    resetRunState();
    celInterIni = celInter;
    dtInterIni = dtInter;
    velInterIni = velInter;

    noextremo = vnoextremo;
    noinicial = vnoinicial;
    derivaAnel = vderivaAnel;

    betaRev = vbetaRev;
    betaRevini = vbetaRevini;
    titRev = vtitRev;
    titRevini = vtitRevini;

    bloq = vbloq;

    dtCicMin = vdtCicMin;

    redeParalelaCCsecundario = -1;
    redeParalelaP = -1;
    redeParalelaS = -1;

    montasistema();
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
    if (celula[ind].acsr.chk.AreaGarg < (1e-3 + arq.master1.razareaativ) * celula[ind].duto.area && celula[ind].acsr.chk.AreaGarg > 1e-3 * celula[ind].duto.area) {
        double tE = celula[ind].temp;
        double alfE = celula[ind].alf;
        double betE = celula[ind].bet;
        double sense = 1.;

        double maxSup = 0.;

        double rholp = celula[ind].flui.MasEspLiq(celula[ind].pres, celula[ind].temp);
        double rholc = celula[ind].fluicol.MasEspFlu(celula[ind].pres, celula[ind].temp);
        double rholmix = (1 - betE) * rholp + betE * rholc;

        double alfJ = celula[ind + 1].alf;
        double betJ = celula[ind + 1].bet;
        double rholpJ = celula[ind + 1].flui.MasEspLiq(celula[ind + 1].pres, celula[ind + 1].temp);
        double rholcJ = celula[ind + 1].fluicol.MasEspFlu(celula[ind + 1].pres, celula[ind + 1].temp);
        double rholmixJ = (1 - betJ) * rholpJ + betJ * rholcJ;

        double hidroM = sin(celula[ind].duto.teta) * (0.5 * celula[ind].dx) * (rholmix * (1 - alfE) + alfE * celula[ind].flui.MasEspGas(celula[ind].pres, celula[ind].temp)) * kGravity / kPascalPerKgfPerCm2Coarse;
        double hidroJ = sin(celula[ind + 1].duto.teta) * (0.5 * celula[ind + 1].dx) * (rholmixJ * (1 - alfJ) + alfJ * celula[ind + 1].flui.MasEspGas(celula[ind + 1].pres, celula[ind + 1].temp)) * kGravity / kPascalPerKgfPerCm2Coarse;

        double masentrada = celula[ind].MC;
        double massgas = celula[ind].MC - celula[ind].Mliqini;
        double tit;
        if (fabs(masentrada < 1e-9) && fabs(massgas < 1e-9)) {
            tit = alfE * celula[ind].flui.MasEspGas(celula[ind].pres, celula[ind].temp) / (celula[ind].flui.MasEspGas(celula[ind].pres, celula[ind].temp) * alfE + rholmix * (1. - alfE));
        } else {
            if ((massgas >= 0 && celula[ind].Mliqini <= 0) || (massgas < 0 && celula[ind].Mliqini == 0))
                tit = 1.;
            else if (massgas <= 0 && celula[ind].Mliqini > 0)
                tit = 0.;
            else if (masentrada < 0)
                tit = 1.;
            else
                tit = fabs(massgas / masentrada);
            if (tit > 1)
                tit = 1;
        }

        double masChk;

        double ypres = (celula[ind + 1].pres + hidroJ) / (celula[ind].pres - hidroM);
        int check = 1;
        if (ypres < 1. || check == 1) {
            if ((*vg1dSP).lixo5 > 7059) {
                int para;
                para = 1;
            }
            masChk = celula[ind].acsr.chk.vazmassSachd(ypres, celula[ind].pres - hidroM, tE, alfE,
                                                       betE, tit, celula[ind].flui, celula[ind].fluicol);
            maxSup = celula[ind].acsr.chk.vazmaxSachd(celula[ind].pres - hidroM, tE, alfE,
                                                      betE, tit, celula[ind].flui, celula[ind].fluicol);
        } else {
            tE = celula[ind + 1].temp;
            alfE = celula[ind + 1].alf;
            betE = celula[ind + 1].bet;
            ypres = 1. / ypres;
            sense = -1;
            if ((*vg1dSP).lixo5 > 7059) {
                int para;
                para = 1;
            }
            masChk = celula[ind].acsr.chk.vazmassSachd(ypres, celula[ind + 1].pres + hidroJ, tE, alfE, betE, tit, celula[ind + 1].flui,
                                                       celula[ind + 1].fluicol);
            maxSup = celula[ind].acsr.chk.vazmaxSachd(celula[ind + 1].pres + hidroJ, tE, alfE, betE, tit, celula[ind + 1].flui,
                                                      celula[ind + 1].fluicol);
        }

        if (fabs(ypres) > fabs(celula[ind].acsr.chk.razpres))
            maxSup = masChk;

        if (celula[ind].acsr.chk.AreaGarg < (1e-3) * celula[ind].duto.area || (celula[ind].pres < celula[ind + 1].pres && check == 1))
            maxSup =
                0.;

        double masliq;
        double masgas;
        masliq = sense * maxSup * (1. - tit);
        masgas = sense * maxSup * tit;

        celula[ind].fontemassLR += (-1 * masliq * (1 - betE) * rholp / rholmix);
        celula[ind].fontemassCR += (-1 * masliq * betE * rholc / rholmix);
        celula[ind].fontemassGR += (-1 * masgas);

        celula[ind + 1].fontemassLR += (masliq * (1 - betE) * rholp / rholmix);
        celula[ind + 1].fontemassCR += (masliq * betE * rholc / rholmix);
        celula[ind + 1].fontemassGR += (masgas);
    }
}

void SProd::salvaFonte() {
    sisprod::transient::storePreviousSources(transientStateOf(*this));
}

/// With the hydrate envelope on (models 2 and 3) and past the first 0.01 s, takes
/// the water and gas that hydrate formation consumed in cell ind during the step,
/// hands them back through the two out-parameters, and lowers the cell's BSW for
/// the free water that is gone.
void SProd::consumeHydrateFormationMass(double &gas_consumido_Mg, double &agua_consumida_Mw, int ind) {
    if (arq.calculaEnvelope == 1 && arq.tipoHmodel == 2 && (*vg1dSP).lixo5 > 0.01) {

        agua_consumida_Mw = celula[ind].agua_consumida_massa_step;

        gas_consumido_Mg = celula[ind].gas_consumido_massa_step;

    // Update the BSW
    double A_cross = celula[ind].duto.area;
    double Lcel    = celula[ind].dx;
    double Vlivre  = std::max(A_cross * Lcel - celula[ind].V_h, 1e-12);

    double frac_agua = std::max((1-celula[ind].alfR)*(1-celula[ind].betR)*celula[ind].FW, 1e-12);
    double frac_oleo = std::max((1-celula[ind].alfR)*(1-celula[ind].betR)*(1-celula[ind].FW), 1e-12);

    double Vagua = frac_agua * Vlivre;
    double Voil  = frac_oleo * Vlivre;

    double rho_w = std::max(celula[ind].flui.MasEspAgua(celula[ind].pres, celula[ind].temp), 1e-12);
    double Vagua_new = Vagua - agua_consumida_Mw / rho_w;
    if (Vagua_new < 0.0) Vagua_new = 0.0;

    double BSW_old = celula[ind].flui.BSW;
    double den = Voil + Vagua_new;
    if (den > 1e-12) {
    celula[ind].flui.BSW = Vagua_new / den;
    } else {
    celula[ind].flui.BSW = BSW_old;
    }
    //celula[ind].FW=celula[ind].flui.BSW;
    if (ind==3) cout << " t [s]: " << (*vg1dSP).lixo5 << " BSW: " << BSW_old << " FW: " << celula[ind].FW << " frac_agua: " << frac_agua << " BSW atualizada apos acoplamento " << celula[ind].flui.BSW << endl;
    //if (ind==3) system("pause");

    } // hydrate change

    if (arq.calculaEnvelope==1 && arq.tipoHmodel==3 && (*vg1dSP).lixo5>0.01) { // hydrate change

    agua_consumida_Mw  = celula[ind].agua_consumida_massa_step;

    gas_consumido_Mg   = celula[ind].gas_consumido_massa_step;

    // Update the BSW
    double A_cross = celula[ind].duto.area;
    double Lcel    = celula[ind].dx;
    double Vlivre  = std::max(A_cross * Lcel - celula[ind].V_h_total, 1e-12);

    double frac_agua = std::max((1-celula[ind].alfR)*(1-celula[ind].betR)*celula[ind].FW, 1e-12);
    double frac_oleo = std::max((1-celula[ind].alfR)*(1-celula[ind].betR)*(1-celula[ind].FW), 1e-12);

    double Vagua = frac_agua * Vlivre;
    double Voil  = frac_oleo * Vlivre;

    double rho_w = std::max(celula[ind].flui.MasEspAgua(celula[ind].pres, celula[ind].temp), 1e-12);
    double Vagua_new = Vagua - agua_consumida_Mw / rho_w;
    if (Vagua_new < 0.0) Vagua_new = 0.0;

    double BSW_old = celula[ind].flui.BSW;
    double den = Voil + Vagua_new;
    if (den > 1e-12) {
    celula[ind].flui.BSW = Vagua_new / den;
    } else {
    celula[ind].flui.BSW = BSW_old;
    }

    //if (ind==3) cout << " t [s]: " << (*vg1dSP).lixo5 << " BSW: " << BSW_old << " FW: " << celula[ind].FW << " frac_agua: " << frac_agua << " BSW atualizada apos acoplamento " << celula[ind].flui.BSW << endl;

    } // hydrate change
}

/// Adds to cell ind the mass its source delivers this step when the source is a
/// choke source (type 9 -- on the first iteration of a parallel network, on its
/// primary side, the flow recorded for that connection instead), a multiple
/// source (10) or a radial or 2D porous medium (15, 16).
void SProd::refreshChokeMultipleAndPorousSources(int ind) {
    if (celula[ind].acsr.tipo == kAccessoryLeak) {
        celula[ind].acsr.fontechk.fluidoP = celula[ind].flui;
        celula[ind].acsr.fontechk.presT = celula[ind].pres;
        celula[ind].acsr.fontechk.tempT = celula[ind].temp;
        double pres = celula[ind].pres;
        double temp = celula[ind].temp;
        double alf = celula[ind].alf;
        double bet = celula[ind].bet;
        double rhog = celula[ind].flui.MasEspGas(pres, temp);
        double rhoP = celula[ind].flui.MasEspLiq(pres, temp);
        double rhoC = celula[ind].fluicol.MasEspFlu(pres, temp);
        celula[ind].acsr.fontechk.titT = alf * rhog / (alf * rhog + (1 - alf) * (bet * rhoC + (1 - bet) * rhoP));
        celula[ind].acsr.fontechk.betIST = bet;
        if ((*vg1dSP).chaveRedeParalela == 0 || (*vg1dSP).iterRede > 0 || redeParalelaS == 1 || redeParalelaCCsecundario == 1) {
            celula[ind].acsr.fontechk.VMas();
            celula[ind].fontemassLR += celula[ind].acsr.fontechk.masP;
            celula[ind].fontemassCR += celula[ind].acsr.fontechk.masC;
            celula[ind].fontemassGR += celula[ind].acsr.fontechk.masG;
        } else {
            int nfonte = indFonteRedeParalelaIni.size();
            int match = 0;
            int iconex;
            for (int ifonte = 0; ifonte < nfonte; ifonte++) {
                if (indFonteRedeParalelaIni[ifonte] == ind) {
                    match = 1;
                    iconex = ifonte;
                    break;
                }
            }
            if (match == 0) {
                celula[ind].acsr.fontechk.VMas();
                celula[ind].fontemassLR += celula[ind].acsr.fontechk.masP;
                celula[ind].fontemassCR += celula[ind].acsr.fontechk.masC;
                celula[ind].fontemassGR += celula[ind].acsr.fontechk.masG;
            } else {
                celula[ind].fontemassLR += fonteMpRedeParalelaIni[iconex];
                celula[ind].fontemassCR += fonteMcRedeParalelaIni[iconex];
                celula[ind].fontemassGR += fonteMgRedeParalelaIni[iconex];
            }
        }
    }
    if (celula[ind].acsr.tipo == kAccessoryMultipleSource) {
        if (celula[ind].acsr.injm.condTermo == 1) {
            celula[ind].fontemassLR = celula[ind].acsr.injm.MassP;
            celula[ind].fontemassCR = celula[ind].acsr.injm.MassC;
            celula[ind].fontemassGR = celula[ind].acsr.injm.MassG;
        } else {
            celula[ind].fontemassCR = celula[ind].acsr.injm.MassC;
            double masT = celula[ind].acsr.injm.MassP + celula[ind].acsr.injm.MassG;
            double tit = celula[ind].acsr.injm.FluidoPro.FracMassHidra(celula[ind].pres, celula[ind].temp);
            celula[ind].acsr.injm.MassG = celula[ind].fontemassGR = masT * tit;
            celula[ind].acsr.injm.MassP = celula[ind].fontemassLR = masT - celula[ind].fontemassGR;
        }
    }
    if (celula[ind].acsr.tipo == kAccessoryRadialPorous) {
        celula[ind].acsr.radialPoro.pW.val[0] = celula[ind].pres;
        double rs;
        double bo;
        double ba;
        if (celula[ind].flui.RGO < 1e7) {
            rs = celula[ind + 1].flui.RS(celula[ind].pres, celula[ind].temp);
            bo = celula[ind + 1].flui.BOFunc(celula[ind].pres, celula[ind].temp, rs);
            ba = celula[ind + 1].flui.BAFunc(celula[ind].pres, celula[ind].temp);
            rs = rs * kBarrelPerCubicMetre / kCubicFootPerCubicMetre;
        } else {
            bo = 1;
            rs = 0;
            ba = 0.;
        }
        // in-situ BSW of the previous cell (in the march, cell i)
        double vfw = celula[ind + 1].flui.BSW * ba / (bo + ba * celula[ind + 1].flui.BSW - celula[ind + 1].flui.BSW * bo);
        celula[ind].acsr.radialPoro.sWPoc = vfw * (1. - celula[ind].acsr.radialPoro.satAconat) + celula[ind].acsr.radialPoro.satAconat;
        celula[ind].acsr.radialPoro.Pint = celula[ind].pres;
        if (modoPerm == 0) {
            celula[ind].acsr.radialPoro.avancoPressao();
            celula[ind].fontemassLR = celula[ind].acsr.radialPoro.fluxIni + celula[ind].acsr.radialPoro.fluxIniA;
            celula[ind].fontemassCR = 0.;
            celula[ind].fontemassGR = celula[ind].acsr.radialPoro.fluxIniG;
        } else {
            celula[ind].acsr.radialPoro.pseudoTrans();
            celula[ind].fontemassLR = celula[ind].acsr.radialPoro.fluxIni + celula[ind].acsr.radialPoro.fluxIniA;
            celula[ind].fontemassCR = 0.;
            celula[ind].fontemassGR = celula[ind].acsr.radialPoro.fluxIniG;
        }
    }
    if (celula[ind].acsr.tipo == kAccessoryPorous2D) {
        celula[ind].acsr.poroso2D.dados.pW.val[0] = celula[ind].pres;
        double rs;
        double bo;
        double ba;
        if (celula[ind].flui.RGO < 1e7) {
            rs = celula[ind + 1].flui.RS(celula[ind].pres, celula[ind].temp);
            bo = celula[ind + 1].flui.BOFunc(celula[ind].pres, celula[ind].temp, rs);
            ba = celula[ind + 1].flui.BAFunc(celula[ind].pres, celula[ind].temp);
            rs = rs * kBarrelPerCubicMetre / kCubicFootPerCubicMetre;
        } else {
            bo = 1;
            rs = 0;
            ba = 0.;
        }
        // in-situ BSW of the previous cell (in the march, cell i)
        double vfw = celula[ind + 1].flui.BSW * ba / (bo + ba * celula[ind + 1].flui.BSW - celula[ind + 1].flui.BSW * bo);
        celula[ind].acsr.poroso2D.sWPoc = vfw * (1. - celula[ind].acsr.poroso2D.dados.satAconat) + celula[ind].acsr.poroso2D.dados.satAconat;
        celula[ind].acsr.poroso2D.dados.transfer.sWPoc = celula[ind].acsr.poroso2D.sWPoc;
        celula[ind].acsr.poroso2D.dados.pInt = celula[ind].pres;
        celula[ind].acsr.poroso2D.dados.transfer.Pint = celula[ind].pres;
        if (modoPerm == 0) {
            celula[ind].acsr.poroso2D.avancoPressao();
            celula[ind].fontemassLR = celula[ind].acsr.poroso2D.dados.transfer.fluxIni + celula[ind].acsr.poroso2D.dados.transfer.fluxIniA;
            celula[ind].fontemassCR = 0.;
            celula[ind].fontemassGR = celula[ind].acsr.poroso2D.dados.transfer.fluxIniG;
        } else {
            celula[ind].acsr.poroso2D.pseudoTransientePoroso();
            celula[ind].fontemassLR = celula[ind].acsr.poroso2D.dados.transfer.fluxIni + celula[ind].acsr.poroso2D.dados.transfer.fluxIniA;
            celula[ind].fontemassCR = 0.;
            celula[ind].fontemassGR = celula[ind].acsr.poroso2D.dados.transfer.fluxIniG;
        }
    }
}

void SProd::renovaFonte(int ind) {

    double pr = celula[ind].pres;
    double tr = celula[ind].temp;
    if (ind == 0) {
        celula[ind].fontemassLR = 0.;
        celula[ind].fontemassCR = 0.;
        celula[ind].fontemassGR = 0.;
    } else {
        if (celula[ind - 1].acsr.tipo != kAccessoryChoke && celula[ind - 1].acsr.tipo != kAccessoryVolumetricPump && (ind < ncel || celula[ncel].acsr.tipo == kAccessoryInflowPerformance)) {
            celula[ind].fontemassLR = 0.;
            celula[ind].fontemassCR = 0.;
            celula[ind].fontemassGR = 0.;
        }
    }

    double agua_consumida_Mw = 0.;
    double gas_consumido_Mg = 0.;

    consumeHydrateFormationMass(gas_consumido_Mg, agua_consumida_Mw, ind);
    if (celula[ind].acsr.tipo == kAccessoryGasInjection) {
        if (celula[ind].acsr.injg.tipoflu == 0) {
            double masgas = celula[ind].acsr.injg.VMas(pr, tr);
            if (fabs(masgas) < (*vg1dSP).localtiny)
                masgas = 0.;
            if (celula[ind].acsr.injg.seco == 1) {
                celula[ind].fontemassGR += masgas;
                celula[ind].fontemassLR = 0.;
                celula[ind].fontemassCR = 0.;
            } else {
                double tit;
                if (arq.flashCompleto != 2)
                    tit = celula[ind].acsr.injg.FluidoPro.FracMassHidra(1., 20.);
                else
                    tit = celula[ind].acsr.injg.FluidoPro.dStockTankVaporMassFraction;
                double masT = masgas / tit;
                tit = celula[ind].acsr.injg.FluidoPro.FracMassHidra(pr, tr);
                celula[ind].fontemassGR += masT * tit;
                celula[ind].fontemassLR += masT * (1. - tit);
                double rcomp = celula[ind].acsr.injg.fluidocol.MasEspFlu(1., 20.);
                celula[ind].fontemassCR += rcomp * celula[ind].acsr.injg.razCompGas * celula[ind].acsr.injg.QGas / kSecondsPerDay;
            }
        } else {
            celula[ind].fontemassGR += 0.;
            celula[ind].fontemassLR = 0.;
            celula[ind].fontemassCR = celula[ind].acsr.injg.QGas;
        }
    }
    if (celula[ind].acsr.tipo == kAccessoryLiquidInjection) {
        double rlcA = celula[ind].acsr.injl.fluidocol.MasEspFlu(1.001, 15.);
        celula[ind].fontemassCR += rlcA * celula[ind].acsr.injl.QLiq * celula[ind].acsr.injl.bet / kSecondsPerDay;
        double massic = celula[ind].acsr.injl.QLiq * (1. - celula[ind].acsr.injl.bet) / kSecondsPerDay;
        double Rhogs = celula[ind].acsr.injl.FluidoPro.Deng * kAirDensityAtStandardConditions; // cel[ind].acsr.injl.FluidoPro.MasEspGas(1, 15);
        double Rhols = (1000 * 141.5 / (131.5 + celula[ind].acsr.injl.FluidoPro.API)) * (1 - celula[ind].acsr.injl.FluidoPro.BSW) + 1000. * celula[ind].acsr.injl.FluidoPro.Denag * celula[ind].acsr.injl.FluidoPro.BSW;
        double multiplicador = (Rhols + celula[ind].acsr.injl.FluidoPro.RGO * Rhogs * (1 - celula[ind].acsr.injl.FluidoPro.BSW));
        massic *= multiplicador;
        double fracmasshidra = celula[ind].acsr.injl.FluidoPro.FracMassHidra(pr, tr);
        celula[ind].fontemassLR += (1. - fracmasshidra) * massic;
        celula[ind].fontemassGR += fracmasshidra * massic;
    }
    if (celula[ind].acsr.tipo == kAccessoryInflowPerformance) {
        if (pr < celula[ind].acsr.ipr.Pres) {
            celula[ind].fontemassLR += celula[ind].acsr.ipr.MasL(pr, tr);
            celula[ind].fontemassCR = 0.;
            celula[ind].fontemassGR += celula[ind].acsr.ipr.MasG(pr, tr);
        } else {
            double tit;
            tit = celula[ind].alf * celula[ind].flui.MasEspGas(pr, tr) /
                  (celula[ind].alf * celula[ind].flui.MasEspGas(pr, tr) +
                   (1. - celula[ind].alf) * (1. - celula[ind].bet) * celula[ind].flui.MasEspLiq(pr, tr) +
                   (1. - celula[ind].alf) * celula[ind].bet * celula[ind].fluicol.MasEspFlu(pr, tr));
            celula[ind].fontemassLR += (1. - celula[ind].alf) * (1. - celula[ind].bet) * celula[ind].acsr.ipr.VMas(pr, tr) * celula[ind].flui.MasEspLiq(pr, tr);
            celula[ind].fontemassCR += (1. - celula[ind].alf) * celula[ind].bet * celula[ind].acsr.ipr.VMas(pr, tr) * celula[ind].fluicol.MasEspFlu(pr, tr);
            celula[ind].fontemassGR += 1 * celula[ind].alf * celula[ind].acsr.ipr.VMas(pr, tr) * celula[ind].flui.MasEspGas(pr, tr);
            celula[ind].acsr.ipr.deriP *= (1. - celula[ind].alf) * (1. - celula[ind].bet) * celula[ind].flui.MasEspLiq(pr, tr);
            celula[ind].acsr.ipr.deriC *= (1. - celula[ind].alf) * celula[ind].bet * celula[ind].fluicol.MasEspFlu(pr, tr);
            celula[ind].acsr.ipr.deriG *= 1 * celula[ind].alf * celula[ind].flui.MasEspGas(pr, tr);
        }
    }
    if (celula[ind].acsr.tipo == kAccessoryChoke) {
        celula[ind + 1].fontemassLR = 0.;
        celula[ind + 1].fontemassCR = 0.;
        celula[ind + 1].fontemassGR = 0.;
        if (modoPerm == 0)
            FonteValv(ind);
    }
    if (celula[ind].acsr.tipo == kAccessoryVolumetricPump) {
        celula[ind + 1].fontemassLR = 0.;
        celula[ind + 1].fontemassCR = 0.;
        celula[ind + 1].fontemassGR = 0.;
        celula[ind].acsr.bvol.fluido = celula[ind].flui;
        celula[ind].acsr.bvol.fluicol = celula[ind].fluicol;
        if (fabs(celula[ind].acsr.bvol.freq) > 1) {
            double alfM = celula[ind].alf;
            double betM = celula[ind].bet;
            double presM = celula[ind].pres;
            double tempM = celula[ind].temp;
            double presM1 = celula[ind].presR;
            double tempM1 = celula[ind].tempR;
            celula[ind].acsr.bvol.vazmass(presM, tempM, presM1, tempM1, betM, alfM);
            celula[ind].fontemassLR -= celula[ind].acsr.bvol.MLiqP;
            celula[ind].fontemassCR -= celula[ind].acsr.bvol.MLiqC;
            celula[ind].fontemassGR -= celula[ind].acsr.bvol.MGas;
            celula[ind + 1].fontemassLR += celula[ind].acsr.bvol.MLiqP;
            celula[ind + 1].fontemassCR += celula[ind].acsr.bvol.MLiqC;
            celula[ind + 1].fontemassGR += celula[ind].acsr.bvol.MGas;
        }
    }
    refreshChokeMultipleAndPorousSources(ind);

    if (arq.calculaEnvelope == 1 && arq.tipoHmodel == 2 && celula[ind].flui.BSW > 1e-12 && (*vg1dSP).lixo5 > 0.01) {
        celula[ind].fontemassLR -= (agua_consumida_Mw / (*vg1dSP).lixo5);
        celula[ind].fontemassGR -= (gas_consumido_Mg / (*vg1dSP).lixo5);
    }
  	// hydrate change 3
  	if (arq.calculaEnvelope==1 && arq.tipoHmodel==3 && celula[ind].flui.BSW>1e-14 && (*vg1dSP).lixo5>0.01) { // hydrate change
  		celula[ind].fontemassLR -= (agua_consumida_Mw / (*vg1dSP).lixo5);
  		celula[ind].fontemassGR -= (gas_consumido_Mg / (*vg1dSP).lixo5);
  	}
}

void SProd::renovaalbetini() {
    sisprod::composition::storePreviousFractionsAndMovePigs(compositionStateOf(*this));
}

void SProd::renovaMasEsp() {
    sisprod::composition::cacheCellAndFaceDensities(compositionStateOf(*this));
}

void SProd::CalcC0Ud(int ind, double &c0, double &ud) {
    const driftflux::coefficient::ClosureState state{
        .cells = celula,
        .lastCell = ncel,
        .globals = vg1dSP,
        .input = arq,
        .selectors = driftSelectors,
        .gasSurfaceTemperature = tGSup,
        .inletVoidFraction = alfE,
        .inletCompletionFraction = betaE,
        .inletPressure = presE,
        .inletTemperature = tempE,
        .steadyIteration = iterperm,
    };
    driftflux::coefficient::instantaneous(state, ind, c0, ud);
}

void SProd::CalcC0UdBuf(int ind, double &c0, double &ud) {
    const driftflux::coefficient::ClosureState state{
        .cells = celula,
        .lastCell = ncel,
        .globals = vg1dSP,
        .input = arq,
        .selectors = driftSelectors,
        .gasSurfaceTemperature = tGSup,
        .inletVoidFraction = alfE,
        .inletCompletionFraction = betaE,
        .inletPressure = presE,
        .inletTemperature = tempE,
        .steadyIteration = iterperm,
    };
    driftflux::coefficient::buffered(state, ind, c0, ud);
}

void SProd::CalcC0UdIni(int ind, double &c0, double &ud) {
    const driftflux::coefficient::ClosureState state{
        .cells = celula,
        .lastCell = ncel,
        .globals = vg1dSP,
        .input = arq,
        .selectors = driftSelectors,
        .gasSurfaceTemperature = tGSup,
        .inletVoidFraction = alfE,
        .inletCompletionFraction = betaE,
        .inletPressure = presE,
        .inletTemperature = tempE,
        .steadyIteration = iterperm,
    };
    driftflux::coefficient::initialization(state, ind, c0, ud);
}

void SProd::CalcC0UdIniBuf(int ind, double &c0, double &ud) {
    const driftflux::coefficient::ClosureState state{
        .cells = celula,
        .lastCell = ncel,
        .globals = vg1dSP,
        .input = arq,
        .selectors = driftSelectors,
        .gasSurfaceTemperature = tGSup,
        .inletVoidFraction = alfE,
        .inletCompletionFraction = betaE,
        .inletPressure = presE,
        .inletTemperature = tempE,
        .steadyIteration = iterperm,
    };
    driftflux::coefficient::bufferedInitialization(state, ind, c0, ud);
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
/// A member of SProd, and it must stay one: both solvers are constructed from
/// the whole SProd object (*this), which a free function taking a state struct
/// does not have. SolveTrans reaches it through a callback, so the hydrate
/// phase runs first in the step.
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
    const driftflux::coefficient::ClosureState state{
        .cells = celula,
        .lastCell = ncel,
        .globals = vg1dSP,
        .input = arq,
        .selectors = driftSelectors,
        .gasSurfaceTemperature = tGSup,
        .inletVoidFraction = alfE,
        .inletCompletionFraction = betaE,
        .inletPressure = presE,
        .inletTemperature = tempE,
        .steadyIteration = iterperm,
    };
    driftflux::coefficient::steadyState(state, ind, c0, ud);
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
