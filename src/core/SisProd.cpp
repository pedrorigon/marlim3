/*
 * SisProd.cpp
 *
 *  Created on: 21 de dez de 2016
 *      Author: Eduardo
 */
#define _USE_MATH_DEFINES // para M_PI
#include "SisProd.h"
#include "DriftFluxClosure.h"
#include "FA_Hidratos.h"
#include "FA_Hidratos_Servico.h"
#include "OutputI18n.h"
#include "RootFindingSolvers.h"
#include "SisProdGasLift.h"
#include "SisProdSteadyState.h"
#include "SisProdSteadyStateSearch.h"
#include "SisProdTransient.h"
#include "SisProdThermal.h"
#include "SisProdTrendOutput.h"
#include <chrono>
#include <math.h>

namespace {
// Defined further down, next to the rest of the thermal plumbing. Declared here
// because the first delegates that need it appear before that definition.
sisprod::thermal::ThermalState thermalStateOf(SProd &system);
sisprod::gaslift::GasLiftState gasLiftStateOf(SProd &system);
}  // namespace

void SProd::resolveDriftSelectors() {
    driftSelectors = {arq.CorreDisper, arq.CorreAnular, arq.CorreEstrat};
}

SProd::SProd(string nomeArquivoEntrada, string nomeArquivoLog, tipoValidacaoJson_t validacaoJson,
             tipoSimulacao_t tipoSimulacao, varGlob1D *Vvg1dSP, int TD, int vbloq, int temporario, int reverso, double *compfonte,
             int *posicfonte, int nfontes, int redeperm) : arq(nomeArquivoEntrada, nomeArquivoLog, validacaoJson, tipoSimulacao, reverso, Vvg1dSP, redeperm), flut(arq.ncelp, arq.nvarprofp + 2 + 1 + 1 + 1 + 1),
                                                           flutG(arq.ncelg, arq.nvarprofg + 2 + 1 + 1 + 1 + 1 + 1), matglobP(2 * arq.ncelp, 3, 2), termolivreP(
                                                                                                                                                       2 * arq.ncelp),
                                                           matglobG(3 * arq.ncelg, 5, 5), termolivreG(3 * arq.ncelg) {
    resolveDriftSelectors();

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
    celInterIni = celInter;
    dtInterIni = dtInter;
    velInterIni = velInter;

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
    RSLivia = 0;
    lerRS = 0;
    noextremo = 1;
    noinicial = 1;
    derivaAnel = -1;

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
    redeTemporario = temporario;

    betaRev = 0;
    betaRevini = 0;
    titRev = 1.;
    titRevini = 1.;

    dtCFLMed = 1.;
    dtSimMed = 1.;
    restriDt = 0;
    kontarestriDt = 0;
    dtCFLTotal = 0.;
    dtSimTotal = 0.;
    dtauxCFL = 0;
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

    bloq = vbloq;

    vg1dSP = Vvg1dSP;

    if (TD >= 0)
        arq.tabelaDinamica = TD;

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

    redeParalelaCCsecundario = -1;
    redeParalelaP = -1;
    redeParalelaS = -1;
    montasistema(compfonte, posicfonte, nfontes);
}
SProd::SProd() : arq(), flut(1, 1 + 2 + 1 + 1 + 1), flutG(1, 1 + 2 + 1 + 1 + 1 + 1), matglobP(2 * 1, 3, 2), termolivreP(2 * 1), matglobG(
                                                                                                                                    3 * 1, 5, 5),
                 termolivreG(3 * 1) {
    resolveDriftSelectors();
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
    kontaTempoProf = 0;
    //kontaTempoCelUni = 0;
    dt = 0.;
    nabreM1 = 0;
    nfechaM1 = 0;
    zdranP = 0;
    dzdpP = 0;
    dzdtP = 0;
    cpg = 0;
    cpl = 0;
    drholdT = 0;
    HLat = 0;
    npontos = 0;
    nfluP = 0;
    chokeVGL = 0;
    posicVGLP = 0;
    posicVGLG = 0;
    fechaM1 = 0;
    abreM1 = 0;
    celulaG = 0;
    celula = 0;
    celInter = 1e7;
    dtInter = 0.;
    velInter = 0.;

    celInterIni = 0.;
    dtInterIni = 0.;
    velInterIni = 0.;

    ModelCp = 0;
    Modeljtl = 0;
    CalcLat = 0;
    trackRGO = 0;
    trackDeng = 0;
    ninjgas = 0;
    lingas = 0;
    injPoc = 0;

    indTramo = -1;
    ncel = 0;
    reinicia = 0;
    presfim = 0;
    presfimini = 0;

    pGSup = 0;
    pGSupIni = 0.;
    tGSup = 0;
    tGSupIni = 0;
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
    receb = 0;

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

    LerPB = 0;
    PBPVTSim = 0;
    TBPVTSim = 0;
    RSLivia = 0;
    lerRS = 0;

    noextremo = 1;
    noinicial = 1;
    derivaAnel = -1;

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

    titRev = 1.;
    titRevini = 1.;

    tempMedContDesc = 10.;
    maxVecContDesc = 1000;
    vazmedDesc = 0;
    tempmedDEsc = 0;
    betaRev = 0.;
    betaRevini = 1.;
    redeTemporario = 0;

    dtCFLMed = 1.;
    dtSimMed = 1.;
    restriDt = 0;
    kontarestriDt = 0;
    dtCFLTotal = 0.;
    dtSimTotal = 0.;
    dtauxCFL = 0.;
    dtauxFinal = 0.;

    kimpT = 0.;

    kontaGolfada = 1000.;

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
    ntabDin = 0;

    for (int i = 0; i < 10; i++) {
        vRazMast0[i] = 0.;
        vRazMast1[i] = 0.;
        vRazMastCrit[i] = 0.5;
    }

    kontaRenovaComp = 0;

    bloq = 0;

    fluiRevRede = ProFlu();
    tempRev = 0.;
    revPerm = 0;

    vg1dSP = 0;

    nCelulaPoisson2D = 0;

    trocaTermicaLenta = 0.01;

    semTermo = 0;
    dtCicMin = dt;

    monitConvPerm = 1000.;
    monitConvPermBase = 1.;

    alteraTempo = 0;

    redeParalelaCCsecundario = -1;
    redeParalelaP = -1;
    redeParalelaS = -1;
}

SProd &SProd::operator=(const SProd &sp) {
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

    if (arq.ntendp > 0 && redeTemporario == 0) {
        for (int i = 0; i < arq.ntendp && MatTrendP && TrendLengthP; i++) {
            if (MatTrendP[i]) {
                for (int j = 0; j < TrendLengthP[i]; j++)
                    delete[] MatTrendP[i][j];
                delete[] MatTrendP[i];
            }
        }
        if (MatTrendP!=0)
            delete[] MatTrendP;
        if (TrendLengthP!=0)
            delete[] TrendLengthP;
        if (resettrend!=0)
            delete[] resettrend;
        if (ntrend!=0)
            delete[] ntrend;
        if (ntrendB!=0)
            delete[] ntrendB;
    }


    if (arq.ntendg > 0 && arq.lingas > 0 && redeTemporario == 0) {
        for (int i = 0; i < arq.ntendg && MatTrendG && TrendLengthG; i++) {
            if (MatTrendG[i]) {
                for (int j = 0; j < TrendLengthG[i]; j++)
                    delete[] MatTrendG[i][j];
                delete[] MatTrendG[i];
            }
        }
        if (MatTrendG!=0)
            delete[] MatTrendG;
        if (TrendLengthG!=0)
            delete[] TrendLengthG;
        if (resettrendg!=0)
            delete[] resettrendg;
        if (ntrendg!=0)
            delete[] ntrendg;
        if (ntrendgB!=0)
            delete[] ntrendgB;
    }

    if (arq.ntendtransp > 0 && redeTemporario == 0) {
        for (int i = 0; i < arq.ntendtransp && MatTrendTransP && TrendLengthTransP; i++) {
            if (MatTrendTransP[i]) {
                for (int j = 0; j < TrendLengthTransP[i]; j++)
                    delete[] MatTrendTransP[i][j];
                delete[] MatTrendTransP[i];
            }
        }
        if (MatTrendTransP!=0)
            delete[] MatTrendTransP;
        if (TrendLengthTransP!=0)
            delete[] TrendLengthTransP;
        if (resettrendtrans!=0)
            delete[] resettrendtrans;
        if (ntrendtrans!=0)
            delete[] ntrendtrans;
        if (ntrendtransB!=0)
            delete[] ntrendtransB;
    }

    if (arq.ntendtransg > 0 && redeTemporario == 0) {
        for (int i = 0; i < arq.ntendtransg && MatTrendTransG && TrendLengthTransG; i++) {
            if (MatTrendTransG[i]) {
                for (int j = 0; j < TrendLengthTransG[i]; j++)
                    delete[] MatTrendTransG[i][j];
                delete[] MatTrendTransG[i];
            }
        }
        if (MatTrendTransG!=0)
            delete[] MatTrendTransG;
        if (TrendLengthTransG!=0)
            delete[] TrendLengthTransG;
        if (resettrendtransg!=0)
            delete[] resettrendtransg;
        if (ntrendtransg!=0)
            delete[] ntrendtransg;
        if (ntrendtransgB!=0)
            delete[] ntrendtransgB;
    }

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

    arq = sp.arq;
    resolveDriftSelectors();
    flut = sp.flut;
    flutG = sp.flutG;
    matglobP = sp.matglobP;
    termolivreP = sp.termolivreP;
    matglobG = sp.matglobG;
    termolivreG = sp.termolivreG;
    vg1dSP = sp.vg1dSP;
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
    LerPB = 0;
    PBPVTSim = 0;
    TBPVTSim = 0;
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
    celInterIni = celInter;
    dtInterIni = dtInter;
    velInterIni = velInter;

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

    noextremo = sp.noextremo;
    noinicial = sp.noinicial;
    derivaAnel = sp.derivaAnel;

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

    betaRev = sp.betaRev;
    betaRevini = sp.betaRevini;
    titRev = sp.titRev;
    titRevini = sp.titRevini;

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
    kontaGolfada = 1000.;

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

    bloq = sp.bloq;

    fluiRevRede = ProFlu();
    tempRev = 0.;
    revPerm = 0;

    ntabDin = 0;

    nCelulaPoisson2D = 0;

    trocaTermicaLenta = 0.01;

    semTermo = 0;

    dtCicMin = sp.dtCicMin;

    monitConvPerm = 1000.;
    monitConvPermBase = 1.;

    alteraTempo = 0;
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

    if (arq.ntendp > 0 && redeTemporario == 0) {
         for (int i = 0; i < arq.ntendp && MatTrendP && TrendLengthP; i++) {
             if (MatTrendP[i]) {
                 for (int j = 0; j < TrendLengthP[i]; j++)
                     delete[] MatTrendP[i][j];
                 delete[] MatTrendP[i];
             }
         }
         if (MatTrendP!=0)
             delete[] MatTrendP;
         if (TrendLengthP!=0)
             delete[] TrendLengthP;
         if (resettrend!=0)
             delete[] resettrend;
         if (ntrend!=0)
             delete[] ntrend;
         if (ntrendB!=0)
             delete[] ntrendB;
     }

     if (arq.ntendg > 0 && arq.lingas > 0 && redeTemporario == 0) {
         for (int i = 0; i < arq.ntendg && MatTrendG && TrendLengthG; i++) {
             if (MatTrendG[i]) {
                 for (int j = 0; j < TrendLengthG[i]; j++)
                     delete[] MatTrendG[i][j];
                 delete[] MatTrendG[i];
             }
         }
         if (MatTrendG!=0)
             delete[] MatTrendG;
         if (TrendLengthG!=0)
             delete[] TrendLengthG;
         if (resettrendg!=0)
             delete[] resettrendg;
         if (ntrendg!=0)
             delete[] ntrendg;
         if (ntrendgB!=0)
             delete[] ntrendgB;
     }

     if (arq.ntendtransp > 0 && redeTemporario == 0) {
         for (int i = 0; i < arq.ntendtransp && MatTrendTransP && TrendLengthTransP; i++) {
             if (MatTrendTransP[i]) {
                 for (int j = 0; j < TrendLengthTransP[i]; j++)
                     delete[] MatTrendTransP[i][j];
                 delete[] MatTrendTransP[i];
             }
         }
         if (MatTrendTransP)
             delete[] MatTrendTransP;
         if (TrendLengthTransP)
             delete[] TrendLengthTransP;
         if (resettrendtrans)
             delete[] resettrendtrans;
         if (ntrendtrans!=0)
             delete[] ntrendtrans;
         if (ntrendtransB!=0)
             delete[] ntrendtransB;
     }

     if (arq.ntendtransg > 0 && redeTemporario == 0) {
         for (int i = 0; i < arq.ntendtransg && MatTrendTransG && TrendLengthTransG; i++) {
             if (MatTrendTransG[i]) {
                 for (int j = 0; j < TrendLengthTransG[i]; j++)
                     delete[] MatTrendTransG[i][j];
                 delete[] MatTrendTransG[i];
             }
         }
         if (MatTrendTransG)
             delete[] MatTrendTransG;
         if (TrendLengthTransG)
             delete[] TrendLengthTransG;
         if (resettrendtransg)
             delete[] resettrendtransg;
         if (ntrendtransg!=0)
             delete[] ntrendtransg;
         if (ntrendtransgB!=0)
             delete[] ntrendtransgB;
     }

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

    arq.copiaSemJson(sp);
    resolveDriftSelectors();
    flut = FullMtx<double>(arq.ncelp, arq.nvarprofp + 2 + 1 + 1 + 1 + 1);
    flutG = FullMtx<double>(arq.ncelg, arq.nvarprofg + 2 + 1 + 1 + 1 + 1 + 1);
    matglobP = BandMtx<double>(2 * arq.ncelp, 3, 2);
    termolivreP = Vcr<double>(2 * arq.ncelp);
    matglobG = BandMtx<double>(3 * arq.ncelg, 5, 5);
    termolivreG = Vcr<double>(3 * arq.ncelg);
    vg1dSP = arq.vg1dSP;
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
    LerPB = 0;
    PBPVTSim = 0;
    TBPVTSim = 0;
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
    celInterIni = celInter;
    dtInterIni = dtInter;
    velInterIni = velInter;

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

    noextremo = vnoextremo;
    noinicial = vnoinicial;
    derivaAnel = vderivaAnel;

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

    betaRev = vbetaRev;
    betaRevini = vbetaRevini;
    titRev = vtitRev;
    titRevini = vtitRevini;

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
    kontaGolfada = 1000.;

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

    bloq = vbloq;

    fluiRevRede = ProFlu();
    tempRev = 0.;
    revPerm = 0;

    ntabDin = 0;

    nCelulaPoisson2D = 0;

    trocaTermicaLenta = 0.01;

    semTermo = 0;

    dtCicMin = vdtCicMin;

    monitConvPerm = 1000.;
    monitConvPermBase = 1.;

    alteraTempo = 0;

    redeParalelaCCsecundario = -1;
    redeParalelaP = -1;
    redeParalelaS = -1;

    montasistema();
}

void SProd::HidroDescargaG() {
    sisprod::gaslift::computeGasUnloadingHydrostatics(gasLiftStateOf(*this));
}

void SProd::HidroDescargaP() {
    celula[0].massfonteCH = 0;
    celula[0].fontemassCL = 0;
    celula[0].fontemassLL = 0;
    celula[0].fontemassGL = 0;
    celula[0].fontemassCR = 0;
    celula[0].fontemassLR = 0;
    celula[0].fontemassGR = 0;
    if (0 <= arq.celdescargaP) {
        celula[0].alf = 0.;
        celula[0].bet = 1.;
        celula[0].betI = celula[0].bet;
        celula[0].alfini = celula[0].alf;
        celula[0].betini = celula[0].bet;
        celula[0].alfPigD = celula[0].alf;
        celula[0].betPigD = celula[0].bet;
        celula[0].alfPigE = celula[0].alf;
        celula[0].betPigE = celula[0].bet;
        celula[1].alfL = celula[0].alf;
        celula[1].betL = celula[0].bet;
        celula[1].betLI = celula[0].bet;
        celula[1].alfLini = celula[0].alf;
        celula[1].betLini = celula[0].bet;
    } else {
        celula[0].alf = 1.;
        celula[0].bet = 0.;
        celula[0].betI = celula[0].bet;
        celula[0].alfini = celula[0].alf;
        celula[0].betini = celula[0].bet;
        celula[0].alfPigD = celula[0].alf;
        celula[0].betPigD = celula[0].bet;
        celula[0].alfPigE = celula[0].alf;
        celula[0].betPigE = celula[0].bet;
        celula[1].alfL = celula[0].alf;
        celula[1].betL = celula[0].bet;
        celula[1].betLI = celula[0].bet;
        celula[1].alfLini = celula[0].alf;
        celula[1].betLini = celula[0].bet;
    }
    double pmed;
    double tmed;
    pmed = celula[0].acsr.ipr.Pres;
    celula[0].presL = pmed;
    celula[1].presL = pmed;
    celula[0].pres = pmed;
    celula[0].presini = pmed;
    celula[1].presLini = pmed;
    celula[0].presauxL = pmed;
    tmed = celula[0].calor.Textern1;
    celula[0].tempL = tmed;
    celula[0].temp = tmed;
    celula[0].tempini = tmed;
    celula[1].tempL = tmed;
    double rho0;
    if (0 <= arq.celdescargaP)
        rho0 = celula[0].fluicol.MasEspFlu(pmed, tmed);
    else
        rho0 = celula[0].flui.MasEspGas(pmed, tmed);
    double rho1;
    rho1 = celula[0].fluicol.MasEspFlu(pmed, tmed);

    celula[0].FW = 0;
    celula[0].FWini = 0;
    celula[0].arranjo = 0;
    celula[0].QLL = 0;
    celula[1].QLL = celula[0].QLL;
    celula[0].QL = 0;
    celula[0].QG = 0;
    celula[0].rpL = rho1;
    celula[0].rpC = rho1;
    celula[1].rpL = celula[0].rpC;
    celula[0].rcL = rho1;
    celula[0].rcC = rho1;
    celula[1].rcL = celula[0].rcC;

    celula[0].MC = 0.;
    celula[0].ML = 0.;
    celula[1].ML = celula[0].MC;
    celula[0].Mliqini = celula[0].MC;
    celula[1].MliqiniL = celula[0].MC;
    celula[0].MliqiniL = celula[0].MC;

    celula[0].rpLi = rho1;
    celula[0].rpCi = rho1;
    celula[0].rcLi = rho1;
    celula[0].rcCi = rho1;

    for (int i = 1; i <= ncel; i++) {
        double A0 = celula[i - 1].duto.area;
        double dx0 = 0.5 * celula[i].dxL;
        double A1 = celula[i].duto.area;
        double dx1 = 0.5 * celula[i].dx;
        pmed -= rho0 * 9.81 * dx0 * sin(celula[i - 1].duto.teta) / 98066.52;
        tmed = celula[i].calor.Textern1;
        double taux = (dx0 * celula[i - 1].temp + dx1 * tmed) / (dx0 + dx1);
        celula[i].presaux = pmed;
        celula[i - 1].presauxR = celula[i].presaux;
        celula[i].presauxL = celula[i - 1].presaux;
        if (i <= arq.celdescargaP)
            rho1 = celula[i].fluicol.MasEspFlu(pmed, tmed);
        else
            rho1 = celula[i].flui.MasEspGas(pmed, tmed);
        pmed -= rho1 * 9.81 * dx1 * sin(celula[i].duto.teta) / 98066.52;
        rho0 = rho1;

        celula[i].presL = celula[i - 1].pres;
        celula[i].pres = pmed;
        celula[i - 1].presR = pmed;
        celula[i].presini = pmed;
        celula[i].presLini = celula[i - 1].presini;
        celula[i].presauxL = celula[i - 1].presaux;
        celula[i].tempL = celula[i - 1].temp;
        celula[i].temp = tmed;
        celula[i - 1].tempR = tmed;
        celula[i].tempini = tmed;

        celula[i].FW = 0;
        celula[i].FWini = 0;
        celula[i].arranjo = 0;
        celula[i].QL = 0;
        celula[i - 1].QLR = celula[i].QL;
        celula[i].QG = 0;
        celula[i].rpL = rho1;
        celula[i].rpC = rho1;
        celula[i - 1].rpR = celula[i].rpC;
        celula[i].rcL = rho1;
        celula[i].rcC = rho1;
        celula[i - 1].rcR = celula[i].rcC;

        double rhoaux;
        if (i <= arq.celdescargaP)
            rhoaux = celula[i].fluicol.MasEspFlu(celula[i].presaux, taux);
        else
            rhoaux = celula[0].flui.MasEspGas(celula[i].presaux, taux);
        celula[i].rpCi = rhoaux;
        celula[i - 1].rpRi = rhoaux;
        celula[i].rcCi = rhoaux;
        celula[i - 1].rcRi = rhoaux;

        celula[i].massfonteCH = 0;
        celula[i].fontemassCL = 0;
        celula[i].fontemassLL = 0;
        celula[i].fontemassGL = 0;
        celula[i].fontemassCR = 0;
        celula[i].fontemassLR = 0;
        celula[i].fontemassGR = 0;
        if (i <= arq.celdescargaP) {
            celula[i].alf = 0.;
            celula[i].bet = 1.;
            celula[i - 1].betR = celula[i].bet;
            celula[i - 1].betRini = celula[i].bet;
            celula[i - 1].alfR = celula[i].alf;
            celula[i - 1].alfRini = celula[i].alf;
            celula[i].betI = celula[i].bet;
            celula[i].alfini = celula[i].alf;
            celula[i].betini = celula[i].bet;
            celula[i].alfPigD = celula[i].alf;
            celula[i].betPigD = celula[i].bet;
            celula[i].alfPigE = celula[i].alf;
            celula[i].betPigE = celula[i].bet;
            celula[i - 1].alfR = celula[i].alf;
            celula[i - 1].alfRini = celula[i].alf;
            celula[i - 1].betR = celula[i].bet;
            celula[i - 1].betRini = celula[i].bet;
        } else {
            celula[i].alf = 1.;
            celula[i].bet = 0.;
            celula[i - 1].betR = celula[i].bet;
            celula[i - 1].betRini = celula[i].bet;
            celula[i - 1].alfR = celula[i].alf;
            celula[i - 1].alfRini = celula[i].alf;
            celula[i].betI = celula[i].bet;
            celula[i].alfini = celula[i].alf;
            celula[i].betini = celula[i].bet;
            celula[i].alfPigD = celula[i].alf;
            celula[i].betPigD = celula[i].bet;
            celula[i].alfPigE = celula[i].alf;
            celula[i].betPigE = celula[i].bet;
            celula[i - 1].alfR = celula[i].alf;
            celula[i - 1].alfRini = celula[i].alf;
            celula[i - 1].betR = celula[i].bet;
            celula[i - 1].betRini = celula[i].bet;
        }

        if (i < ncel) {
            celula[i + 1].tempL = tmed;
            celula[i + 1].QLL = celula[i].QL;
            celula[i + 1].rpL = celula[i].rpC;
            celula[i + 1].rcL = celula[i].rcC;
            celula[i + 1].rpLi = rhoaux;
            celula[i + 1].rcLi = rhoaux;
            celula[i + 1].alfL = celula[i].alf;
            celula[i + 1].betL = celula[i].bet;
            celula[i + 1].betLI = celula[i].bet;
            celula[i + 1].alfLini = celula[i].alf;
            celula[i + 1].betLini = celula[i].bet;
        }
    }
}

/// Points every cell fluid, and every source fluid it carries, at the bubble-point
/// tables read from the PVTSim file, and switches them to saturation model 4.
void SProd::assignPvtSimBubbleTablesToCells() {
    for (int i = 0; i <= ncel; i++) {
        celula[i].flui.PBPVTSim = PBPVTSim;
        celula[i].flui.TBPVTSim = TBPVTSim;
        celula[i].flui.corrSat = 4;
        if (celula[i].acsr.tipo == 2) {
            celula[i].acsr.injl.FluidoPro.PBPVTSim = PBPVTSim;
            celula[i].acsr.injl.FluidoPro.TBPVTSim = TBPVTSim;
            celula[i].acsr.injl.FluidoPro.corrSat = 4;
        }
        if (celula[i].acsr.tipo == 3) {
            celula[i].acsr.ipr.FluidoPro.PBPVTSim = PBPVTSim;
            celula[i].acsr.ipr.FluidoPro.TBPVTSim = TBPVTSim;
            celula[i].acsr.ipr.FluidoPro.corrSat = 4;
        }
        if (celula[i].acsr.tipo == 10) {
            celula[i].acsr.injm.FluidoPro.PBPVTSim = PBPVTSim;
            celula[i].acsr.injm.FluidoPro.TBPVTSim = TBPVTSim;
            celula[i].acsr.injm.FluidoPro.corrSat = 4;
        }
        if (celula[i].acsr.tipo == 15) {
            celula[i].acsr.radialPoro.flup.PBPVTSim = PBPVTSim;
            celula[i].acsr.radialPoro.flup.TBPVTSim = TBPVTSim;
            celula[i].acsr.radialPoro.flup.corrSat = 4;
            for (int iRP = 0; iRP < celula[i].acsr.radialPoro.ncel; iRP++) {
                celula[i].acsr.radialPoro.celula[iRP].flup.PBPVTSim = PBPVTSim;
                celula[i].acsr.radialPoro.celula[iRP].flup.TBPVTSim = TBPVTSim;
                celula[i].acsr.radialPoro.celula[iRP].flup.corrSat = 4;
            }
        }
        if (celula[i].acsr.tipo == 16) {
            celula[i].acsr.poroso2D.dados.flup.PBPVTSim = PBPVTSim;
            celula[i].acsr.poroso2D.dados.flup.TBPVTSim = TBPVTSim;
            celula[i].acsr.poroso2D.dados.flup.corrSat = 4;
            for (int iRP = 0; iRP < celula[i].acsr.poroso2D.dados.transfer.ncel; iRP++) {
                celula[i].acsr.poroso2D.dados.transfer.celula[iRP].flup.PBPVTSim = PBPVTSim;
                celula[i].acsr.poroso2D.dados.transfer.celula[iRP].flup.TBPVTSim = TBPVTSim;
                celula[i].acsr.poroso2D.dados.transfer.celula[iRP].flup.corrSat = 4;
            }
            for (int iRP = 0; iRP < celula[i].acsr.poroso2D.malha.nele; iRP++) {
                celula[i].acsr.poroso2D.malha.mlh2d[iRP].flup.PBPVTSim = PBPVTSim;
                celula[i].acsr.poroso2D.malha.mlh2d[iRP].flup.TBPVTSim = TBPVTSim;
                celula[i].acsr.poroso2D.malha.mlh2d[iRP].flup.corrSat = 4;
            }
        }
    }
}

/// Reads the bubble-point curve from the PVTSim file, points the cell fluids at it
/// and writes perfilBolha; with tabRSPB on, also reads the solution gas-oil ratio
/// table and writes perfilRSLivia.
void SProd::loadPvtSimSaturationTables() {
    LerPB = 1;
    int ndiv = arq.tabent.npont - 1;
    PBPVTSim = new double[ndiv + 1];
    TBPVTSim = new double[ndiv + 1];
    double *PresPVTSim;
    PresPVTSim = new double[ndiv + 1];

    string impfile;
    impfile = arq.pvtsimarq;
    string dadosMR = impfile;
    ifstream lendoPVTSim(dadosMR.c_str(), ios_base::in);
    string chave;
    char *tenta;
    tenta = new char[400];
    double testatok;
    char line[4000];
    lendoPVTSim.get(line, 4000);
    tenta = strtok(line, " ,()=");
    lendoPVTSim >> chave;
    while (chave != "PRESSURE") {
        lendoPVTSim >> chave;
    }
    int lacoleitura = ndiv;
    lendoPVTSim.get(line, 4000);
    tenta = strtok(line, " ,()=");
    PresPVTSim[0] = atof(tenta) / 98068.059233;
    double valor;
    for (int kontaPVT = 1; kontaPVT <= lacoleitura; kontaPVT++) {
        tenta = strtok(NULL, " ,");
        testatok = atof(tenta);
        PresPVTSim[kontaPVT] = testatok / 98068.059233;
    }
    while (chave != "BUBBLEPRESSURES") {
        lendoPVTSim >> chave;
    }
    lendoPVTSim.get(line, 4000);
    tenta = strtok(line, " ,()=");
    PBPVTSim[0] = atof(tenta) * 0.00014503773800722;
    for (int kontaPVT = 1; kontaPVT <= lacoleitura; kontaPVT++) {
        tenta = strtok(NULL, " ,");
        testatok = atof(tenta);
        PBPVTSim[kontaPVT] = testatok * 0.00014503773800722;
    }
    while (chave != "BUBBLETEMPERATURES") {
        lendoPVTSim >> chave;
    }
    lendoPVTSim.get(line, 4000);
    tenta = strtok(line, " ,()=");
    TBPVTSim[0] = atof(tenta);
    for (int kontaPVT = 1; kontaPVT <= lacoleitura; kontaPVT++) {
        tenta = strtok(NULL, " ,");
        testatok = atof(tenta);
        TBPVTSim[kontaPVT] = testatok;
    }
    assignPvtSimBubbleTablesToCells();

    FullMtx<double> BolhaTemp(ndiv + 2, 2);
    for (int i = 0; i <= ndiv; i++) {
        BolhaTemp[i][0] = TBPVTSim[i];
        BolhaTemp[i][1] = PBPVTSim[i];
    }
    ostringstream saidaBolha;
    saidaBolha << pathPrefixoArqSaida << "perfilBolha";
    string tmp = saidaBolha.str();
    ofstream escreveMass(tmp.c_str(), ios_base::out);
    escreveMass << BolhaTemp;
    escreveMass.close();

    if (arq.tabRSPB == 1) {
        arq.tabRSPB = 1;
        lerRS = 1;
        RSLivia = new double *[ndiv + 2];
        for (int i = 0; i < ndiv + 2; i++) {
            RSLivia[i] = new double[ndiv + 2];
        }
        for (int i = 1; i <= ndiv + 1; i++) {
            RSLivia[i][0] = PresPVTSim[i - 1];
            RSLivia[0][i] = TBPVTSim[i - 1];
            for (int j = 1; j <= ndiv + 1; j++) {
                if (TBPVTSim[j - 1] > -10) {
                    double a1, a2, a3;
                    double TFa = 1.8 * TBPVTSim[j - 1] + 32;

                    double A0, A1, A2, A3, A4;
                    double B0, B1, B2, B3, B4;
                    double C0, C1, C2, C3, C4;
                    double D0, D1, D2;

                    A0 = 6542.69213 * 1e-11;
                    A1 = 2.60464618;
                    A2 = 1.21334544;
                    A3 = 0.16464125;
                    A4 = 1.19382967;
                    B0 = 87 * 1e-8;
                    B1 = 4.77390016;
                    B2 = 1.77267703;
                    B3 = -0.73072038;
                    B4 = 1.58093857;
                    C0 = 3226.09 * 1e-6;
                    C1 = 0.09281075;
                    C2 = 0.13633665;
                    C3 = 0.09634381;
                    C4 = 0.53238728;
                    D0 = 1.00544053;
                    D1 = -0.00134177;
                    D2 = 0.51839397;

                    double yco2 = celula[0].flui.yco2;
                    double Deng = celula[0].flui.Deng;
                    double API = celula[0].flui.API;
                    double multCor = (D0 + D1 * yco2 * pow(TFa, D2));
                    double pbtemp = PBPVTSim[j - 1];

                    a1 = A0 * pow(Deng, A1) * pow(API, A2) * pow(TFa, A3) * pow(pbtemp, A4);
                    a2 = B0 * pow(Deng, B1) * pow(API, B2) * pow(TFa, B3) * pow(pbtemp, B4);
                    a3 = C0 * pow(Deng, C1) * pow(API, C2) * pow(TFa, C3) * pow(pbtemp, C4);

                    double pcor = (RSLivia[i][0] * 0.9678411) * 14.69595 * multCor;
                    double pr = pcor / pbtemp;

                    a1 = A0 * pow(Deng, A1) * pow(API, A2) * pow(TFa, A3) * pow(pbtemp, A4);
                    a2 = B0 * pow(Deng, B1) * pow(API, B2) * pow(TFa, B3) * pow(pbtemp, B4);
                    a3 = C0 * pow(Deng, C1) * pow(API, C2) * pow(TFa, C3) * pow(pbtemp, C4);

                    double Rsr = a1 * pow(pr, a2) + (1. - a1) * pow(pr, a3);
                    double rstemp = 0.;
                    if (rstemp < 0.)
                        rstemp = 0.;
                    else if ((RSLivia[i][0] * 0.9678411) * 14.69595 > pbtemp)
                        rstemp = 1.;
                    else
                        rstemp = Rsr;
                    RSLivia[i][j] = rstemp;
                } else
                    RSLivia[i][j] = 0;
            }
        }
        FullMtx<double> RSTemp(ndiv + 2, ndiv + 2);
        for (int i = 1; i <= ndiv + 1; i++) {
            RSTemp[i][0] = RSLivia[i][0];
            RSTemp[0][i] = RSLivia[0][i];
            for (int j = 1; j <= ndiv + 1; j++) {
                if (i > 0 || j > 0)
                    RSTemp[i][j] = RSLivia[i][j] * 6.29 / 35.31467;
            }
        }
        ostringstream saidaRS;
        saidaRS << pathPrefixoArqSaida << "perfilRSLivia";
        tmp = saidaRS.str();
        ofstream escreveRS(tmp.c_str(), ios_base::out);
        escreveRS << RSTemp;
        escreveRS.close();

        for (int i = 0; i <= ncel; i++) {
            celula[i].flui.TabRSLivia = RSLivia;
            celula[i].flui.tabRSPB = 1;
            if (celula[i].acsr.tipo == 2) {
                celula[i].acsr.injl.FluidoPro.TabRSLivia = RSLivia;
                celula[i].acsr.injl.FluidoPro.tabRSPB = 1;
            }
            if (celula[i].acsr.tipo == 3) {
                celula[i].acsr.ipr.FluidoPro.TabRSLivia = RSLivia;
                celula[i].acsr.ipr.FluidoPro.tabRSPB = 1;
            }
            if (celula[i].acsr.tipo == 10) {
                celula[i].acsr.injm.FluidoPro.TabRSLivia = RSLivia;
                celula[i].acsr.injm.FluidoPro.tabRSPB = 1;
            }
            if (celula[i].acsr.tipo == 15) {
                celula[i].acsr.radialPoro.flup.TabRSLivia = RSLivia;
                celula[i].acsr.radialPoro.flup.tabRSPB = 1;
                for (int iRP = 0; iRP < celula[i].acsr.radialPoro.ncel; iRP++) {
                    celula[i].acsr.radialPoro.celula[iRP].flup.TabRSLivia = RSLivia;
                    celula[i].acsr.radialPoro.celula[iRP].flup.tabRSPB = 1;
                }
            }
            if (celula[i].acsr.tipo == 16) {
                celula[i].acsr.poroso2D.dados.flup.TabRSLivia = RSLivia;
                celula[i].acsr.poroso2D.dados.flup.tabRSPB = 1;
                for (int iRP = 0; iRP < celula[i].acsr.poroso2D.dados.transfer.ncel; iRP++) {
                    celula[i].acsr.poroso2D.dados.transfer.celula[iRP].flup.TabRSLivia = RSLivia;
                    celula[i].acsr.poroso2D.dados.transfer.celula[iRP].flup.tabRSPB = 1;
                }
                for (int iRP = 0; iRP < celula[i].acsr.poroso2D.malha.nele; iRP++) {
                    celula[i].acsr.poroso2D.malha.mlh2d[iRP].flup.TabRSLivia = RSLivia;
                    celula[i].acsr.poroso2D.malha.mlh2d[iRP].flup.tabRSPB = 1;
                }
            }
        }
    }
    delete[] PresPVTSim;
}

/// Builds the bubble-point curve and the solution gas-oil ratio table from the
/// fluid correlations over the input table's pressure-temperature grid, writes
/// perfilBolha and perfilRSLivia, and points the cell fluids at both.
void SProd::generateSaturationTablesFromCorrelations() {
    int ndiv = arq.tabent.npont - 1;
    double pteste = arq.tabent.pmin;
    double tteste = arq.tabent.tmin;
    double dpteste = (arq.tabent.pmax - pteste) / ndiv;
    double dtteste = (arq.tabent.tmax - tteste) / ndiv;
    FullMtx<double> RSTemp(ndiv + 2, ndiv + 2);
    FullMtx<double> PBTemp(ndiv + 2, 2);
    TBPVTSim = new double[ndiv + 1];
    PBPVTSim = new double[ndiv + 1];
    RSLivia = new double *[ndiv + 2];
    for (int i = 0; i < ndiv + 2; i++) {
        RSLivia[i] = new double[ndiv + 2];
    }
    double ttestepb = tteste;
    for (int i = 0; i <= ndiv; i++) {
        TBPVTSim[i] = ttestepb;
        PBPVTSim[i] = arq.flup[0].PB(pteste, ttestepb);
        PBTemp[i][0] = TBPVTSim[i];
        PBTemp[i][1] = PBPVTSim[i];
        ttestepb += dtteste;
    }

    for (int i = 1; i <= ndiv + 1; i++) {
        RSLivia[i][0] = pteste;
        RSTemp[i][0] = RSLivia[i][0];
        ttestepb = tteste;
        for (int j = 1; j <= ndiv + 1; j++) {
            RSLivia[0][j] = ttestepb;
            RSLivia[i][j] = arq.flup[0].RS(pteste, ttestepb);
            RSTemp[0][j] = RSLivia[0][j];
            RSTemp[i][j] = RSLivia[i][j] * 6.29 / 35.31467;
            ttestepb += dtteste;
        }
        pteste += dpteste;
    }

    ostringstream saidaBolha;
    saidaBolha << pathPrefixoArqSaida << "perfilBolha";
    string tmp = saidaBolha.str();
    ofstream escreveMass(tmp.c_str(), ios_base::out);
    escreveMass << PBTemp;
    escreveMass.close();

    ostringstream saidaRS;
    saidaRS << pathPrefixoArqSaida << "perfilRSLivia";
    tmp = saidaRS.str();
    ofstream escreveRS(tmp.c_str(), ios_base::out);
    escreveRS << RSTemp;
    escreveRS.close();

    for (int i = 0; i <= ncel; i++) {
        celula[i].flui.PBPVTSim = PBPVTSim;
        celula[i].flui.TBPVTSim = TBPVTSim;
        celula[i].flui.TabRSLivia = RSLivia;
        celula[i].flui.tabRSPB = 1;

        if (celula[i].acsr.tipo == 2) {
            celula[i].acsr.injl.FluidoPro.PBPVTSim = PBPVTSim;
            celula[i].acsr.injl.FluidoPro.TBPVTSim = TBPVTSim;
            celula[i].acsr.injl.FluidoPro.TabRSLivia = RSLivia;
            celula[i].acsr.injl.FluidoPro.tabRSPB = 1;
        }
        if (celula[i].acsr.tipo == 3) {
            celula[i].acsr.ipr.FluidoPro.PBPVTSim = PBPVTSim;
            celula[i].acsr.ipr.FluidoPro.TBPVTSim = TBPVTSim;
            celula[i].acsr.ipr.FluidoPro.TabRSLivia = RSLivia;
            celula[i].acsr.ipr.FluidoPro.tabRSPB = 1;
        }
        if (celula[i].acsr.tipo == 10) {
            celula[i].acsr.injm.FluidoPro.PBPVTSim = PBPVTSim;
            celula[i].acsr.injm.FluidoPro.TBPVTSim = TBPVTSim;
            celula[i].acsr.injm.FluidoPro.TabRSLivia = RSLivia;
            celula[i].acsr.injm.FluidoPro.tabRSPB = 1;
        }
        if (celula[i].acsr.tipo == 15) {
            celula[i].acsr.radialPoro.flup.PBPVTSim = PBPVTSim;
            celula[i].acsr.radialPoro.flup.TBPVTSim = TBPVTSim;
            celula[i].acsr.radialPoro.flup.TabRSLivia = RSLivia;
            celula[i].acsr.radialPoro.flup.tabRSPB = 1;
            for (int iRP = 0; iRP < celula[i].acsr.radialPoro.ncel; iRP++) {
                celula[i].acsr.radialPoro.celula[iRP].flup.PBPVTSim = PBPVTSim;
                celula[i].acsr.radialPoro.celula[iRP].flup.TBPVTSim = TBPVTSim;
                celula[i].acsr.radialPoro.celula[iRP].flup.TabRSLivia = RSLivia;
                celula[i].acsr.radialPoro.celula[iRP].flup.tabRSPB = 1;
            }
        }
        if (celula[i].acsr.tipo == 16) {
            celula[i].acsr.poroso2D.dados.flup.PBPVTSim = PBPVTSim;
            celula[i].acsr.poroso2D.dados.flup.TBPVTSim = TBPVTSim;
            celula[i].acsr.poroso2D.dados.flup.TabRSLivia = RSLivia;
            celula[i].acsr.poroso2D.dados.flup.tabRSPB = 1;
            for (int iRP = 0; iRP < celula[i].acsr.poroso2D.dados.transfer.ncel; iRP++) {
                celula[i].acsr.poroso2D.dados.transfer.celula[iRP].flup.PBPVTSim = PBPVTSim;
                celula[i].acsr.poroso2D.dados.transfer.celula[iRP].flup.TBPVTSim = TBPVTSim;
                celula[i].acsr.poroso2D.dados.transfer.celula[iRP].flup.TabRSLivia = RSLivia;
                celula[i].acsr.poroso2D.dados.transfer.celula[iRP].flup.tabRSPB = 1;
            }
            for (int iRP = 0; iRP < celula[i].acsr.poroso2D.malha.nele; iRP++) {
                celula[i].acsr.poroso2D.malha.mlh2d[iRP].flup.PBPVTSim = PBPVTSim;
                celula[i].acsr.poroso2D.malha.mlh2d[iRP].flup.TBPVTSim = TBPVTSim;
                celula[i].acsr.poroso2D.malha.mlh2d[iRP].flup.TabRSLivia = RSLivia;
                celula[i].acsr.poroso2D.malha.mlh2d[iRP].flup.tabRSPB = 1;
            }
        }
    }

    LerPB = 1;
    lerRS = 1;
}

/// Copies the run configuration into the members, checks that an injection well
/// has the IPR its boundary condition needs, builds the Cp and JTL tables, hands
/// the fluid constants to every fluid, generates the pipe and the production
/// cells, and places every source -- including the extra gas sources at
/// compfonte, whose cell indices it records in posicfonte.
void SProd::buildProductionCells(double *compfonte, int *posicfonte, int nfontes) {
    indTramo = -1;
    KontaImprime = 0;
    tempoaberto = 0.;
    tempoabertoini = 0.;
    modoPerm = 0;
    iterperm = 0;
    masChkSup = 0;
    dt = arq.dtmax;
    dtini = arq.dtmax;
    dtInter = arq.dtmax;
    tfinal = arq.tfinal;
    TransMassModel = 0;
    indpigP = 0;
    indpigPini = indpigP;
    reinicia = 0;
    masChkSupini = 0;
    betaRev = 0.;
    trackRGO = arq.trackRGO;
    trackDeng = arq.trackDeng;
    ninjgas = arq.ninjgas;
    lingas = arq.lingas;
    nfluP = arq.nfluP;
    ModelCp = arq.modelcp;
    Modeljtl = arq.modelJTL;
    CalcLat = arq.latente;
    if (arq.flashCompleto == 1)
        arq.trackRGO = 1;
    trackRGO = arq.trackRGO;
    trackDeng = arq.trackDeng;
    ninjgas = arq.ninjgas;
    lingas = arq.lingas;
    injPoc = arq.pocinjec + arq.condpocinj.tipoFlui;
    arq.fluc.injPoc = injPoc;
    if (arq.flashCompleto == 1)
        arq.fluc.injPoc = 0;
    (*vg1dSP).localtiny = arq.mono;
    (*vg1dSP).CritCond = arq.critcond;
    titRev = -1;
    if (injPoc >= 1 && arq.nipr == 0 && arq.condpocinj.CC != 3 && arq.condpocinj.CC != 4 && arq.condpocinj.CC != 5)
        NumError(
            "O simulador esta no modo Injecao de Agua em uma condicao de contorno que pede uma IPR e nao foi incluido nenhuma IPR no sistema");
    else if (injPoc >= 1 && arq.condpocinj.CC != 3 && arq.condpocinj.CC != 4 && arq.condpocinj.CC != 5) {
        int testaIPR = 0;
        for (int i = 0; i < arq.nipr; i++)
            if (arq.IPRS[i].indcel == arq.ncelp - 1)
                testaIPR = 1;
        if (testaIPR == 0)
            NumError(
                "O simulador esta no modo Injecao de Agua em uma condicao de contorno que pede uma IPR e a ultima celula nao contem uma IPR, e necessario neste modo se ter uma IPR na ultima celula");
    }

    npontos = arq.tabent.npont;
    ModelCp = arq.modelcp;
    if (arq.modelcp > 0)
        arq.geraTabCp();
    Modeljtl = arq.modelJTL;
    if (arq.modelJTL == 1)
        arq.geraTabDrholDt();
    cpg = arq.cpg;
    cpl = arq.cpl;
    drholdT = arq.drholdT;

    zdranP = arq.zdranP;
    dzdpP = arq.dzdpP;
    dzdtP = arq.dzdtP;
    fluiRevRede = arq.flup[0];
    if (arq.flashCompleto == 2)
        fluiRevRede.atualizaPropCompStandard();
    for (int i = 0; i < arq.nfluP; i++) {
        arq.flup[i].npontos = arq.tabent.npont;
        arq.flup[i].cpg = cpg;
        arq.flup[i].cpl = cpl;
        arq.flup[i].drholdT = drholdT;
        arq.flup[i].nfluP = nfluP;
    }
    arq.flug.npontos = arq.tabent.npont;
    arq.flug.cpg = cpg;
    arq.flug.cpl = cpl;
    arq.flug.drholdT = drholdT;
    arq.flug.nfluP = nfluP;

    arq.fluc.npontos = arq.tabent.npont;
    arq.fluc.RhoInj = arq.RhoInj;
    arq.fluc.ViscInj = arq.ViscInj;
    arq.fluc.CondInj = arq.CondInj;
    arq.fluc.CpInj = arq.CpInj;
    arq.fluc.DrhoDtInj = arq.DrhoDtInj;

    if (arq.modelcp > 0 || arq.modelJTL == 1 || arq.latente > 0 || (arq.pocinjec == 1 && (arq.condpocinj.tipoFlui == 2 || arq.condpocinj.tipoFlui == 3)))
        TransMassModel = arq.transmass;
    arq.geraduto();
    ncel = arq.ncelp - 1;
    temperatura = arq.celp[0].textern;
    tempRev = arq.tempReves;
    celula = new Cel[ncel + 1];
    arq.geracelp(celula);
    if (arq.nipr > 0)
        arq.geraipr(celula);
    if (arq.nvalvgas > 0 && arq.lingas > 0)
        arq.gerafgasVGL(celula);
    if (arq.ninjgas > 0)
        arq.gerafgasFonte(celula);
    for (int i = arq.ninjgas; i < arq.ninjgas + nfontes; i++) {
        int iposp = arq.buscaIndiceMeioP(compfonte[i - arq.ninjgas]);
        posicfonte[i - arq.ninjgas] = iposp;
        InjGas injgasMRT(0, 0, arq.flug);
        celula[iposp].acsr.tipo = 1;
        celula[iposp].acsr.injg = injgasMRT;
    }
    if (arq.ninjliq > 0)
        arq.gerafliqFonte(celula);
    if (arq.ninjmass > 0)
        arq.gerafmassFonte(celula);
    if (arq.nfuro > 0)
        arq.geraFuro(celula);
    if (arq.nPoroRad > 0)
        arq.gerafPoroRadFonte(celula);
    if (arq.nPoro2D > 0)
        arq.gerafPoro2DFonte(celula);
}

/// Gives the inlet the sources its boundary condition needs (and the second cell,
/// under a blockage), then places the accessories -- pumps, volumetric pumps,
/// pressure-drop requirements, heat sources, the master valve and the other
/// valves -- and sets up the outlet pressure, the surface and injection chokes and
/// the pigs.
void SProd::configureInletSourcesAndAccessories(int nfontes) {
    if (celula[0].acsr.tipo != 0 && arq.ConContEntrada == 1) {
        // RN-302: The inlet has both a pressure boundary condition and a mass source. Report a warning.
        logger.log(LOGGER_AVISO, LOG_ERR_PARSE_BUSINESS_RULE_VALIDATION,
                   "Foi escolhida uma condicao de contorno de pressao no inicio da tubulacao e foi colocada uma fonte de massa neste inicio, isto pode causar inconsistencias",
                   "", "");
    }
    if ((arq.perm == 1 || (*vg1dSP).chaverede > 0) && celula[0].acsr.tipo == 0) {
        if ((arq.ConContEntrada == 0 && nfontes == 0 && (*vg1dSP).chaverede == 0)) {
            NumError("O simulador pede calculo de permanente com condicao de vazao na entrada, mas sem nenhuma fonte na entrada");
        }
        if ((arq.ConContEntrada == 0 && nfontes == 0 && (*vg1dSP).chaverede > 0) || arq.CCPres.tit[0] < (1 - (*vg1dSP).localtiny)) {
            if (arq.tipoFluido == 1)
                celula[0].acsr.tipo = 1;
            else if (arq.tipoFluido == 0)
                celula[0].acsr.tipo = 2;
            else
                celula[0].acsr.tipo = 10;
            InjMult injmassMRT(0, 0, 0, celula[0].temp, celula[0].flui, celula[0].fluicol);
            celula[0].acsr.injm.condTermo = 1;
            celula[0].acsr.injm = injmassMRT;
            InjLiq injliqMRT(0, 0, 0, celula[0].flui, celula[0].fluicol);
            celula[0].acsr.injl = injliqMRT;
            InjGas injgasMRT(0, 0, celula[0].flui, celula[0].fluicol);
            celula[0].acsr.injg = injgasMRT;
        } else if (nfontes == 0 || arq.ConContEntrada == 1) {
            InjGas injgasMRT(0, 0, celula[0].flui, celula[0].fluicol);
            injgasMRT.seco = 0;
            celula[0].acsr.injg = injgasMRT;
            InjMult injmassMRT(0, 0, 0, celula[0].temp, celula[0].flui, celula[0].fluicol);
            celula[0].acsr.injm.condTermo = 1;
            celula[0].acsr.injm = injmassMRT;
            InjLiq injliqMRT(0, 0, 0, celula[0].flui, celula[0].fluicol);
            celula[0].acsr.injl = injliqMRT;
            if (arq.tipoFluido == 1)
                celula[0].acsr.tipo = 1;
            else if (arq.tipoFluido == 0)
                celula[0].acsr.tipo = 2;
            else
                celula[0].acsr.tipo = 10;
        }
    }
    if (bloq == 1 && celula[1].acsr.tipo == 0) {
        if ((*vg1dSP).fluidoRede == 1)
            celula[1].acsr.tipo = 2;
        else if ((*vg1dSP).fluidoRede == 0)
            celula[1].acsr.tipo = 1;
        else
            celula[1].acsr.tipo = 10;
        InjMult injmassMRT(0, 0, 0, celula[0].temp, celula[0].flui, celula[0].fluicol);
        celula[1].acsr.injm.condTermo = 1;
        celula[1].acsr.injm = injmassMRT;
        InjLiq injliqMRT(0, 0, 0, celula[0].flui, celula[0].fluicol);
        celula[1].acsr.injl = injliqMRT;
        InjGas injgasMRT(0, 0, celula[0].flui, celula[0].fluicol);
        celula[1].acsr.injg = injgasMRT;
    } else if (bloq == 1) {
        logger.log(LOGGER_AVISO, LOG_ERR_PARSE_BUSINESS_RULE_VALIDATION,
                   "Foi escolhida uma condicao de bloqueio que afeta este tramo, para que a selecao de fontes seja"
                   "feita corretamente, e necessario que a segunda celula do tramo nao tenha fontes catalogadas",
                   "", "");
    }
    if (arq.nbcs > 0)
        arq.gerafBCS(celula);
    if (arq.nmultibcs > 0)
        arq.gerafmultiBCS(celula);
    if (arq.nbvol > 0)
        arq.gerafBVOL(celula);
    if (arq.ndpreq > 0)
        arq.geraDPReq(celula);
    if (arq.ncalor > 0)
        arq.geraFonteCalor(celula);
    if (arq.master1.posic > 0.)
        arq.geraMaster1(celula);
    else {
        int verifica = 0;
        double Lverifica = 0.5 * celula[verifica].dx;
        while (celula[verifica].acsr.tipo != 0) {
            verifica++;
            Lverifica += 0.5 * (celula[verifica].dx + celula[verifica - 1].dx);
        }
        arq.master1.posic = verifica;
        arq.master1.comp = Lverifica;
        arq.geraMaster1(celula);
    }
    if (arq.nvalv > 0)
        arq.geraValv(celula);
    arq.gerapresfim(presfim, pGSup);
    pGSupIni = pGSup;
    chokeSup = choke(1., 1.);
    chokeInj = ChokeGas();
    arq.gerachokesup(chokeSup);
    npig = arq.npig;
    if (npig > 0) {
        receb = new int[npig];
        for (int i = 0; i < npig; i++)
            receb[i] = arq.pig[i].receb;
    }
}

/// Builds the gas-lift line when there is one: its cells, the second master valve,
/// the injection choke, one gas-lift valve choke per valve with its position on
/// both lines, and each gas cell's share of the annulus above the discharge cell.
/// Without a gas line, only zeroes the gas cell count.
void SProd::buildGasLiftLine() {
    if (arq.lingas == 0)
        ncelGas = 0;
    else if (arq.lingas > 0) {
        ncelGas = arq.ncelg;
        ncelGas--;
        celInter = arq.celdescarga;
        celulaG = new CelG[ncelGas + 1];
        arq.geracelg(celulaG);
        arq.geraMaster2(celulaG);
        arq.gerachokeinj(chokeInj);
        if (arq.nvalvgas > 0) {
            chokeVGL = new ChokeGas[arq.nvalvgas];
            posicVGLP = new int[arq.nvalvgas];
            posicVGLG = new int[arq.nvalvgas];
        }
        for (int i = 0; i < arq.nvalvgas; i++) {
            double diaG = arq.valvgl[i].diagarg;
            double presEstag = celulaG[arq.valvgl[i].posicG].pres;
            double tempEstag = celulaG[arq.valvgl[i].posicG].temp;
            double presGarg = celula[arq.valvgl[i].posicP].pres;
            posicVGLP[i] = arq.valvgl[i].posicP;
            posicVGLG[i] = arq.valvgl[i].posicG;
            celulaG[posicVGLG[i]].vgl = 1;
            if (arq.valvgl[i].tipo == 2) {
                arq.valvgl[i].frec = 0.;
                arq.valvgl[i].cd = 1.;
            }
            chokeVGL[i] = ChokeGas(arq.flug, M_PI * diaG * diaG / 4.,
                                   arq.valvgl[i].diaexter, arq.valvgl[i].cd, presEstag, presGarg,
                                   tempEstag, arq.valvgl[i].frec, arq.valvgl[i].tipo,
                                   (M_PI * diaG * diaG / 4.) / arq.valvgl[i].razarea,
                                   arq.valvgl[i].pcali, arq.valvgl[i].tcali, arq.valvgl[i].cdLiq, arq.valvgl[i].frecLiq);
            celulaG[posicVGLG[i]].pEstag = chokeVGL[i].presEstag;
            celulaG[posicVGLG[i]].tEstag = chokeVGL[i].tempEstag;
            celulaG[posicVGLG[i]].pGarg = chokeVGL[i].presGarg;
            celulaG[posicVGLG[i]].tGarg = chokeVGL[i].tempGarg;
            celulaG[posicVGLG[i]].qGarg = chokeVGL[i].qGarg;
            celulaG[posicVGLG[i]].areaGarg = chokeVGL[i].areagarg;
        }
        for (int i = 0; i <= ncelGas; i++) {
            celulaG[i].celInter = &celInter;
            celulaG[i].razInter = 1.;
            celulaG[i].razInterIni = 1.;
            if (i >= celInter) {
                celulaG[i].razInter = 0.;
                celulaG[i].razInterIni = 0.;
            }
        }
    }
}

/// Rejects or warns about source and boundary-condition combinations the run
/// cannot honour (RN-300, RN-301, gas-lift discharge without an IPR, no outlet
/// pressure, accessories in the last two cells), builds the gas-lift discharge
/// hydrostatics and the event log, applies the initial state of a production
/// well, and records the surface temperature and mass flow it starts from.
void SProd::validateSetupAndApplyInitialState() {
    if (celula[0].acsr.tipo != 1 && celula[0].acsr.tipo != 2 && celula[0].acsr.tipo != 3 && celula[0].acsr.tipo != 10 && celula[0].acsr.tipo != 15 && celula[0].acsr.tipo != 16) {
        if (arq.perm == 0 && arq.ConContEntrada == 0) {
            // RN-300: No source is defined in the first production system cell. Report a warning.
            logger.log(LOGGER_AVISO, LOG_ERR_PARSE_BUSINESS_RULE_VALIDATION,
                       "Nao existe nenhuma fonte na primeira celula do sistema de producao, esta nao e a condicao mais adequada para a simulacao transiente, problemas podem ocorrer nesta situacao, embora nao seja impeditivo para o modo",
                       "", "");
        } else {
            // RN-301: No source is defined in the first production system cell
            // under steady-state conditions. Report a failure.
            NumError("Nao existe nenhuma fonte na primeira celula do sistema de producao, no modo permanente esta e uma condicao necessaria");
        }
    }

    if (arq.descarga == 1 && celula[0].acsr.tipo != 3) {
        NumError(
            "Nao existe nenhuma IPR na primeira celula do sistema de producao, no modo descarga de gas lift esta e uma condicao necessaria");
    }
    if (arq.descarga == 1) {
        HidroDescargaG();
        HidroDescargaP();
    }

    arq.geraevento(noinicial, noextremo);
    if (arq.saidaClassica == 1 && (arq.transiente == 1 || (*vg1dSP).chaveredeT == 1)) {
        cout << "------------------------------------------   Eventos no Tramo -------------------------------------------" << endl;
        int contaImp = 0;
        for (int i = 0; i < arq.logevento.size(); i++)
            if (arq.logevento[i].descricao != "Gravando Perfil Linha de Producao" && arq.logevento[i].descricao != "Gravando Perfil Linha de Gas") {
                cout << contaImp << "  momento = " << arq.logevento[i].instante << " segundos " << " tipo de evento = " << arq.logevento[i].descricao << endl;
                contaImp++;
            }
    }

    if (arq.pocinjec == 0 && arq.perm == 1 && ((*vg1dSP).chaverede == 0 || noextremo == 1) && arq.ConContEntrada != 2 && arq.psep.pres[0] < -1e6) {
        NumError("sem pressao a jusante como condição de contorno para um tramo onde esta condicao é necessaria ");
    }

    if ((celula[ncel].acsr.tipo != 0 || celula[ncel - 1].acsr.tipo != 0) && arq.transiente == 1) {
        NumError("As duas celulas finais não devem ter acessorios, estas duas celulas devem ser preservadas, principalmente para simulacoes transientes ");
    }

    (*vg1dSP).lixo5 = 0.;
    (*vg1dSP).contador = 0;
    if (arq.pocinjec == 0) {
        arq.atualiza(noinicial, noextremo, derivaAnel, chokeSup, chokeInj, celula, celulaG, pGSup, temperatura, presiniG,
                     tempiniG, presE, tempE, titE, betaE, (*vg1dSP).lixo5, dt);
        pGSupIni = pGSup;
        // presiniG,tempiniG,presE,tempE,titE,betaE,(*vg1dSP).lixo5);//alteracao7
        if (chokeSup.AreaGarg >= 0.6 * celula[ncel - 1].duto.area) {
            aberto = 1;
            abertoini = 1;
        } else {
            aberto = 0;
            abertoini = 1;
        }
    }

    tempSup = celula[ncel].temp;
    masSup = celula[ncel - 1].MC;

    tGSup = celula[ncel].calor.Textern1;
}

/// With the dynamic property table on, splits the pipe into table segments that
/// end at every source cell. Then gives each cell its left and right inclination:
/// its own, or, where its pipe is horizontal, the nearest inclined neighbour's on
/// that side.
void SProd::buildDynamicTablesAndInclinations() {
    if (arq.tabelaDinamica == 1) {
        int minNPontos = 0;
        ntabDin = 1;
        tabelaDinamica temp;
        temp.npontosP = minNPontos;
        temp.npontosT = minNPontos;
        temp.celIni = 0;
        temp.rhogF = 0;
        temp.rholF = 0;
        temp.DrhogDpF = 0;
        temp.DrhogDtF = 0;
        temp.DrholDpF = 0;
        temp.DrholDtF = 0;
        temp.valBO = 0;
        temp.HgF = 0;
        temp.HlF = 0;
        temp.cpgF = 0;
        temp.cplF = 0;
        temp.valZ = 0;
        temp.valdZdT = 0;
        temp.valdZdP = 0;
        temp.tit = 0;
        temp.rs = 0;
        temp.viscG = 0;
        temp.viscO = 0;
        temp.TBF = 0;
        temp.PBF = 0;
        tabDin.push_back(temp);
        for (int i = 1; i < ncel; i++) {
            if ((celula[i].acsr.tipo == 1) ||
                (celula[i].acsr.tipo == 2) ||
                celula[i].acsr.tipo == 3 ||
                (celula[i].acsr.tipo == 10) || celula[i].acsr.tipo == 9 || celula[i].acsr.tipo == 15 || celula[i].acsr.tipo == 16) {
                tabelaDinamica temp;
                temp.npontosP = minNPontos;
                temp.npontosT = minNPontos;
                temp.celIni = i + 1;
                tabDin[ntabDin - 1].celFim = i;
                tabDin.push_back(temp);
                ntabDin++;
            }
        }
        tabDin[ntabDin - 1].celFim = ncel;
    }

    celula[0].angEsq = celula[0].duto.teta;
    for (int i = 1; i <= ncel; i++) {
        if (fabs(celula[i].duto.teta) < 1e-10)
            celula[i].angEsq = celula[i - 1].angEsq;
        else
            celula[i].angEsq = celula[i].duto.teta;
    }
    celula[ncel].angEsq = celula[ncel].duto.teta;
    for (int i = ncel - 1; i >= 0; i--) {
        if (fabs(celula[i].duto.teta) < 1e-10)
            celula[i].angDir = celula[i + 1].angDir;
        else
            celula[i].angDir = celula[i].duto.teta;
    }
}

/// Sets the latent-heat switch from the input and, when it is on and flashCompleto
/// is 0, reads the latent-heat table from the PVTSim file into HLat and writes
/// perfilLatente.
void SProd::configureLatentHeat() {
    CalcLat = arq.latente;
    if (arq.latente > 0) {
        if (arq.flashCompleto == 0) {
            string impfile;
            impfile = arq.pvtsimarq;
            string dadosMR = impfile;
            ifstream lendoPVTSim(dadosMR.c_str(), ios_base::in);
            string chave;
            char *tenta;
            tenta = new char[400];
            double testatok;
            int ndiv = arq.tabent.npont - 1;
            Vcr<double> presPVTSim(ndiv + 1, 0.);
            Vcr<double> tempPVTSim(ndiv + 1, 0.);

            FullMtx<double> HLatTemp(ndiv + 2, ndiv + 2);
            char line[4000];
            lendoPVTSim.get(line, 4000);
            tenta = strtok(line, " ,()=");
            while (strcmp(tenta, "PHASE") != 0) {
                tenta = strtok(NULL, " ,()=");
            }
            tenta = strtok(NULL, " ,()=");
            int lacoleitura = 12;
            if (strcmp(tenta, "THREE") == 0)
                lacoleitura = 18;

            while (chave != "PRESSURE")
                lendoPVTSim >> chave;
            lendoPVTSim.get(line, 4000);
            tenta = strtok(line, " ,()=");
            presPVTSim[0] = atof(tenta) * 1.01971621e-5;
            for (int kontaPVT = 1; kontaPVT <= ndiv; kontaPVT++) {
                tenta = strtok(NULL, " ,");
                presPVTSim[kontaPVT] = atof(tenta) * 1.01971621e-5;
            }
            while (chave != "TEMPERATURE")
                lendoPVTSim >> chave;
            lendoPVTSim.get(line, 4000);
            tenta = strtok(line, " ,()=");
            tempPVTSim[0] = atof(tenta);
            for (int kontaPVT = 1; kontaPVT <= ndiv; kontaPVT++) {
                tenta = strtok(NULL, " ,");
                tempPVTSim[kontaPVT] = atof(tenta);
            }

            for (int i = 1; i <= ndiv + 1; i++) {
                HLatTemp[i][0] = presPVTSim[i - 1];
                for (int j = 1; j <= ndiv + 1; j++) {
                    HLatTemp[0][j] = tempPVTSim[j - 1];
                    while (chave != "POINT")
                        lendoPVTSim >> chave;
                    lendoPVTSim.get(line, 4000);
                    tenta = strtok(line, " ,()=");
                    for (int kontaPVT = 0; kontaPVT < lacoleitura; kontaPVT++) {
                        tenta = strtok(NULL, " ,");
                        testatok = atof(tenta);
                    }
                    tenta = strtok(NULL, " ,");
                    testatok = atof(tenta);
                    HLatTemp[i][j] = testatok;
                    tenta = strtok(NULL, " ,");
                    testatok = atof(tenta);
                    HLatTemp[i][j] -= testatok;
                    if (i == ndiv + 1 && j == ndiv + 1)
                        break;
                    while (chave != "PVTTABLE")
                        lendoPVTSim >> chave;
                }
            }
            lendoPVTSim.close();
            HLat = new double *[ndiv + 2];
            for (int i = 0; i < ndiv + 2; i++) {
                HLat[i] = new double[ndiv + 2];
                for (int j = 0; j < ndiv + 2; j++)
                    HLat[i][j] = HLatTemp[i][j];
            }
            ostringstream saidaLatente;
            saidaLatente << pathPrefixoArqSaida << "perfilLatente.dat";
            string tmp = saidaLatente.str();
            ofstream escreveMass(tmp.c_str(), ios_base::out);
            escreveMass << HLatTemp;
            escreveMass.close();
            // caso nao seja simulacao POCO_INJETOR
            if (arq.tipoSimulacao != tipoSimulacao_t::poco_injetor) {
                arqRelatorioPerfis << tmp.c_str() << endl;
                arqRelatorioPerfis.flush();
            }
        }
    }
}

/// Sets the gas-density correction factors of every fluid in every cell -- the
/// cell's, its source's and its reservoir cells' -- to 1 when the correction is
/// off, or evaluates them at local pressure and temperature when it is on; then
/// makes the first cell's fluid the fluid of its source.
void SProd::applyDensityCorrectionsAndInletFluid() {
    for (int i = 0; i <= ncel; i++) {
        if (arq.corrDeng == 0) {
            celula[i].flui.rDgD = 1.;
            celula[i].flui.rDgL = 1.;
            celula[i].flui.PCis = celula[i].flui.PC;
            celula[i].flui.TCis = celula[i].flui.TC;
            if (celula[i].acsr.tipo == 2) {
                celula[i].acsr.injl.FluidoPro.rDgD = 1.;
                celula[i].acsr.injl.FluidoPro.rDgL = 1.;
            } else if (celula[i].acsr.tipo == 3) {
                celula[i].acsr.ipr.FluidoPro.rDgD = 1.;
                celula[i].acsr.ipr.FluidoPro.rDgL = 1.;
            } else if (celula[i].acsr.tipo == 10) {
                celula[i].acsr.injm.FluidoPro.rDgD = 1.;
                celula[i].acsr.injm.FluidoPro.rDgL = 1.;
            }
            if (celula[i].acsr.tipo == 15) {
                celula[i].acsr.radialPoro.flup.rDgD = 1.;
                celula[i].acsr.radialPoro.flup.rDgD = 1.;
                for (int iRP = 0; iRP < celula[i].acsr.radialPoro.ncel; iRP++) {
                    celula[i].acsr.radialPoro.celula[iRP].flup.rDgD = 1.;
                    celula[i].acsr.radialPoro.celula[iRP].flup.rDgD = 1.;
                }
            }
            if (celula[i].acsr.tipo == 16) {
                celula[i].acsr.poroso2D.dados.flup.rDgD = 1.;
                celula[i].acsr.poroso2D.dados.flup.rDgD = 1.;
                for (int iRP = 0; iRP < celula[i].acsr.poroso2D.dados.transfer.ncel; iRP++) {
                    celula[i].acsr.poroso2D.dados.transfer.celula[iRP].flup.rDgD = 1.;
                    celula[i].acsr.poroso2D.dados.transfer.celula[iRP].flup.rDgL = 1.;
                }
                for (int iRP = 0; iRP < celula[i].acsr.poroso2D.malha.nele; iRP++) {
                    celula[i].acsr.poroso2D.malha.mlh2d[iRP].flup.rDgD = 1.;
                    celula[i].acsr.poroso2D.malha.mlh2d[iRP].flup.rDgL = 1.;
                }
            }
        } else {
            celula[i].flui.razDegD(celula[i].pres, celula[i].temp);
            celula[i].flui.rzDegL(celula[i].pres, celula[i].temp);
            celula[i].flui.PcTcIS();
            if (celula[i].acsr.tipo == 2) {
                celula[i].acsr.injl.FluidoPro.razDegD(celula[i].pres, celula[i].temp);
                celula[i].acsr.injl.FluidoPro.rzDegL(celula[i].pres, celula[i].temp);
            } else if (celula[i].acsr.tipo == 3) {
                celula[i].acsr.ipr.FluidoPro.razDegD(celula[i].pres, celula[i].temp);
                celula[i].acsr.ipr.FluidoPro.rzDegL(celula[i].pres, celula[i].temp);
            } else if (celula[i].acsr.tipo == 10) {
                celula[i].acsr.injm.FluidoPro.razDegD(celula[i].pres, celula[i].temp);
                celula[i].acsr.injm.FluidoPro.rzDegL(celula[i].pres, celula[i].temp);
            } else if (celula[i].acsr.tipo == 15) {
                celula[i].acsr.radialPoro.flup.razDegD(celula[i].pres, celula[i].temp);
                celula[i].acsr.radialPoro.flup.rzDegL(celula[i].pres, celula[i].temp);
                for (int iRP = 0; iRP < celula[i].acsr.radialPoro.ncel; iRP++) {
                    double pres = celula[i].acsr.radialPoro.celula[iRP].Pcamada;
                    double temp = celula[i].acsr.radialPoro.tRes;
                    celula[i].acsr.radialPoro.celula[iRP].flup.razDegD(pres, temp);
                    celula[i].acsr.radialPoro.celula[iRP].flup.rzDegL(pres, temp);
                }
            } else if (celula[i].acsr.tipo == 16) {
                celula[i].acsr.poroso2D.dados.flup.razDegD(celula[i].pres, celula[i].temp);
                celula[i].acsr.poroso2D.dados.flup.rzDegL(celula[i].pres, celula[i].temp);
                for (int iRP = 0; iRP < celula[i].acsr.poroso2D.dados.transfer.ncel; iRP++) {
                    double pres = celula[i].acsr.poroso2D.dados.transfer.celula[iRP].Pcamada;
                    double temp = celula[i].acsr.poroso2D.dados.transfer.tRes;
                    celula[i].acsr.poroso2D.dados.transfer.celula[iRP].flup.razDegD(pres, temp);
                    celula[i].acsr.poroso2D.dados.transfer.celula[iRP].flup.rzDegL(pres, temp);
                }
                for (int iRP = 0; iRP < celula[i].acsr.poroso2D.malha.nele; iRP++) {
                    double pres = celula[i].acsr.poroso2D.malha.mlh2d[iRP].cel2D.presC;
                    double temp = celula[i].acsr.poroso2D.malha.mlh2d[iRP].tRes;
                    celula[i].acsr.poroso2D.malha.mlh2d[iRP].flup.razDegD(pres, temp);
                    celula[i].acsr.poroso2D.malha.mlh2d[iRP].flup.rzDegL(pres, temp);
                }
            }
        }
    }

    if (celula[0].acsr.tipo == 1) {
        celula[0].flui = celula[0].acsr.injg.FluidoPro;
    } else if (celula[0].acsr.tipo == 2) {
        celula[0].flui = celula[0].acsr.injl.FluidoPro;
    } else if (celula[0].acsr.tipo == 3) {
        celula[0].flui = celula[0].acsr.ipr.FluidoPro;
    } else if (celula[0].acsr.tipo == 10) {
        celula[0].flui = celula[0].acsr.injm.FluidoPro;
        ;
    } else if (celula[0].acsr.tipo == 15) {
        celula[0].flui = celula[0].acsr.radialPoro.flup;
    } else if (celula[0].acsr.tipo == 16) {
        celula[0].flui = celula[0].acsr.poroso2D.dados.flup;
    }
}

/// Lists the cells with two-dimensional heat diffusion, copies the master-valve
/// opening and closing times, maps the transient profile positions to global
/// thermal-node indices, and -- unless the network is only provisional --
/// allocates the trend matrices of the production line, the gas line and the
/// transient trends, sized to the longest simulated time and filled with the
/// -10000 sentinel.
void SProd::allocateEventProfileAndTrendArrays() {
    for (int i = 0; i < ncel; i++) {
        if (celula[i].calor.difus2D == 1) {
            indCelPoisson2D.push_back(i);
            nCelulaPoisson2D++;
        }
    }

    nabreM1 = arq.eventoabre;
    nfechaM1 = arq.eventofecha;
    abreM1 = new double[nabreM1];
    fechaM1 = new double[nfechaM1];
    for (int i = 0; i < nabreM1; i++)
        abreM1[i] = arq.Tevento[i];
    for (int i = 0; i < nfechaM1; i++)
        fechaM1[i] = arq.Teventof[i];

    int ntempGas = 0;
    if (arq.lingas > 0)
        ntempGas = ncelGas;

    if (arq.nperfistransp > 0) {
        ncelperftransp = new int[arq.nperfistransp];
        for (int i = 0; i < arq.nperfistransp; i++) {
            int posiccel = arq.proftransp.posic[i];
            ncelperftransp[i] = celula[posiccel].calor.nglobal;
        }
    }
    if (arq.nperfistransg > 0 && arq.lingas > 0) {
        ncelperftransg = new int[arq.nperfistransg];
        for (int i = 0; i < arq.nperfistransg; i++) {
            int posiccel = arq.proftransg.posic[i];
            ncelperftransg[i] = celulaG[posiccel].calor.nglobal;
        }
    }

    if ((*vg1dSP).chaverede == 0)
        (*vg1dSP).TmaxR = arq.tfinal;
    if (redeTemporario == 0) {
        if (arq.ntendp > 0) {
            TrendLengthP = new int[arq.ntendp];
            for (int i = 0; i < arq.ntendp; i++)
                TrendLengthP[i] = 1 + 1 + ceil((*vg1dSP).TmaxR / arq.trendp[i].dt); // round(arq.tfinal / arq.trendp[i].dt);
        }
        if (arq.ntendp > 0) {
            resettrend = new double[arq.ntendp];
            ntrend = new int[arq.ntendp];
            ntrendB = new int[arq.ntendp];
            MatTrendP = new double **[arq.ntendp];
            for (int i = 0; i < arq.ntendp; i++) {
                MatTrendP[i] = new double *[TrendLengthP[i]];
                for (int j = 0; j < TrendLengthP[i]; j++) {
                    MatTrendP[i][j] = new double[arq.nvartrendp[i] + 2];
                    for (int k = 0; k <= arq.nvartrendp[i]; k++)
                        MatTrendP[i][j][k] = -10000.;
                }
                resettrend[i] = 0;
                ntrend[i] = 0;
                ntrendB[i] = 0;
            }
        }

        if (arq.ntendg > 0) {
            TrendLengthG = new int[arq.ntendg];
            for (int i = 0; i < arq.ntendg; i++)
                TrendLengthG[i] = 1 + 1 + ceil((*vg1dSP).TmaxR / arq.trendg[i].dt);
        }
        if (arq.ntendg > 0) {
            resettrendg = new double[arq.ntendg];
            ntrendg = new int[arq.ntendg];
            ntrendgB = new int[arq.ntendg];
            MatTrendG = new double **[arq.ntendg];
            for (int i = 0; i < arq.ntendg; i++) {
                MatTrendG[i] = new double *[TrendLengthG[i]];
                for (int j = 0; j < TrendLengthG[i]; j++) {
                    MatTrendG[i][j] = new double[arq.nvartrendg[i] + 2];
                    for (int k = 0; k <= arq.nvartrendg[i]; k++)
                        MatTrendG[i][j][k] = -10000.;
                }
                resettrendg[i] = 0;
                ntrendg[i] = 0;
                ntrendgB[i] = 0;
            }
        }
        if (arq.ntendtransp > 0) {
            TrendLengthTransP = new int[arq.ntendtransp];
            for (int i = 0; i < arq.ntendtransp; i++)
                TrendLengthTransP[i] = 1 + 1 + ceil((*vg1dSP).TmaxR / arq.trendtransp[i].dt);
        }
        if (arq.ntendtransp > 0) {
            resettrendtrans = new double[arq.ntendtransp];
            ntrendtrans = new int[arq.ntendtransp];
            ntrendtransB = new int[arq.ntendtransp];
            MatTrendTransP = new double **[arq.ntendtransp];
            for (int i = 0; i < arq.ntendtransp; i++) {
                MatTrendTransP[i] = new double *[TrendLengthTransP[i]];
                for (int j = 0; j < TrendLengthTransP[i]; j++)
                    MatTrendTransP[i][j] = new double[2];
                for (int j = 0; j < TrendLengthTransP[i]; j++)
                    for (int k = 0; k < 2; k++)
                        MatTrendTransP[i][j][k] = -10000.;

                resettrendtrans[i] = 0;
                ntrendtrans[i] = 0;
                ntrendtransB[i] = 0;
            }
        }
        if (arq.ntendtransg > 0 && arq.lingas > 0) {
            TrendLengthTransG = new int[arq.ntendtransg];
            for (int i = 0; i < arq.ntendtransg; i++)
                TrendLengthTransG[i] = 1 + 1 + ceil((*vg1dSP).TmaxR / arq.trendtransg[i].dt);
        }
        if (arq.ntendtransg > 0) {
            resettrendtransg = new double[arq.ntendtransg];
            ntrendtransg = new int[arq.ntendtransg];
            ntrendtransgB = new int[arq.ntendtransg];
            MatTrendTransG = new double **[arq.ntendtransg];
            for (int i = 0; i < arq.ntendtransg; i++) {
                MatTrendTransG[i] = new double *[TrendLengthTransG[i]];
                for (int j = 0; j < TrendLengthTransG[i]; j++)
                    MatTrendTransG[i][j] = new double[2];
                for (int j = 0; j < TrendLengthTransG[i]; j++)
                    for (int k = 0; k < 2; k++)
                        MatTrendTransG[i][j][k] = -10000.;

                resettrendtransg[i] = 0;
                ntrendtransg[i] = 0;
                ntrendtransgB[i] = 0;
            }
        }
    }
}

/// Resets the column-annulus and network coupling flags, locates the
/// column-annulus coupling range on the gas line, sets the steady and transient
/// profile counters and their first output times, opens the event log with the
/// events known at start, finds the smallest cell length, and zeroes the
/// moving-average and running-total state the transient loop starts from.
void SProd::resetCouplingAndOutputState() {
    verificaAcop = 0;

    if (arq.lingas > 0) {
        AnulaColunaIni = arq.anulcoluini();
        AnulaColunaFim = arq.anulcolufim();
        ColunaAnulaIni = arq.coluanulini();
        ColunaAnulaFim = arq.coluanulfim();
        if (ColunaAnulaFim == -1) {
            ColunaAnulaFim = 0;
            AnulaColunaFim--;
        }
        if (AnulaColunaIni >= 0 && AnulaColunaFim >= 0 && ColunaAnulaIni >= 0 && ColunaAnulaFim >= 0)
            verificaAcop = 1;
    }

    verificaAcopRedeP = 0;
    verificaAcopRedeS = 0;
    SecPrimIniRedeP = 0;
    SecPrimFimRedeP = 0;
    PrimSecIniRedeP = 0;
    PrimSecFimRedeP = 0;

    if(arq.AP==0){
    	if (arq.nperfisp > 0) {
    		arq.imprimeProfile(celula, flut, (*vg1dSP).lixo5, indTramo);
    	}
    	if (arq.nperfisg > 0 && arq.lingas > 0) {
    		arq.imprimeProfileG(celulaG, flutG, (*vg1dSP).lixo5, indTramo);
    	}
    	if (arq.nperfistransp > 0) {
    		arq.imprimeProfileTrans(celula, ncelperftransp, (*vg1dSP).lixo5, indTramo);
    	}
    	if (arq.nperfistransg > 0 && arq.lingas > 0) {
    		arq.imprimeProfileTransG(celulaG, ncelperftransg, (*vg1dSP).lixo5, indTramo);
    	}
    }
    kontaTempoProf = 0;
    kontaTempoProfG = 0;
    if (arq.nperfisp > 0) {
        if (arq.profp.tempo[0] <= 0 + (*vg1dSP).localtiny)
            kontaTempoProf++;
    }
    if (arq.nperfisg > 0 && arq.lingas > 0) {
        if (arq.profg.tempo[0] <= 0 + (*vg1dSP).localtiny)
            kontaTempoProfG++;
    }
    kontaTempoTransProf = 0;
    kontaTempoTransProfG = 0;
    if (arq.nperfistransp > 0) {
        if (arq.proftransp.tempo[0] <= 0 + (*vg1dSP).localtiny)
            kontaTempoTransProf++;
    }
    if (arq.nperfistransg > 0 && arq.lingas > 0) {
        if (arq.proftransg.tempo[0] <= 0 + (*vg1dSP).localtiny)
            kontaTempoTransProfG++;
    }

    saidaLog << pathPrefixoArqSaida << "LogEvento" << ".dat";
    tmpLog = saidaLog.str();
    // if it's not a simulation POCO_INJETOR
    if (arq.tipoSimulacao != tipoSimulacao_t::poco_injetor) {
        contaLog = 0;
        int nevent = arq.logevento.size();
        while (fabs(arq.logevento[contaLog].instante - (*vg1dSP).lixo5) < dt && contaLog < nevent) {
            time_t now = time(0);
            tm *ltm = localtime(&now); // Retirado de https://www.tutorialspoint.com/cplusplus/cpp_date_time.htm
            ofstream escreveIni(tmpLog.c_str(), ios_base::app);
            escreveIni << "Evento Externo = ";
            escreveIni << arq.logevento[contaLog].instante << " ; ";
            escreveIni << arq.logevento[contaLog].duracao << " ; ";
            escreveIni << arq.logevento[contaLog].descricao << " ; ";
            escreveIni << "datahora = ";
            escreveIni << ltm->tm_mday << "/";
            escreveIni << 1 + ltm->tm_mon << "/";
            escreveIni << 1900 + ltm->tm_year << " ";
            escreveIni << 0 + ltm->tm_hour << ":";
            escreveIni << 0 + ltm->tm_min << ":";
            escreveIni << 0 + ltm->tm_sec;
            escreveIni << endl;
            contaLog++;
            escreveIni.close();
        }
    }

    menorDx = 1e10;
    for (int i = 0; i <= ncel; i++) {
        if (celula[i].dx < menorDx)
            menorDx = celula[i].dx;
        if (i > 0)
            celula[i].razdxTM = celula[i - 1].dx / (celula[i - 1].dx + celula[i].dx);
        if (i > 1)
            celula[i].razdxTM0 = celula[i - 2].dx / (celula[i - 1].dx + celula[i - 2].dx);
    }

    if (arq.nCelUnit>0) {
    	for(int iCelU=0;iCelU<arq.nCelUnit;iCelU++){
    		kontaTempoCelUni.push_back(0);
    	}
    }

    kSP = 0;
    indevento = 1;
    mult = 0.8;
    presMedMov = 0.;
    jMedMov = 0.;
    tMedMov = 60.;
    ktMedMov = 0.;
    pTotal = 0.;
    jTotal = 0.;
    alfTotal = 0.;
}

void SProd::montasistema(double *compfonte, int *posicfonte, int nfontes) {
    try {
        buildProductionCells(compfonte, posicfonte, nfontes);
        configureInletSourcesAndAccessories(nfontes);
        buildGasLiftLine();
        validateSetupAndApplyInitialState();
        buildDynamicTablesAndInclinations();
        configureLatentHeat();
        if (celula[0].flui.corrSat == -4) {
            loadPvtSimSaturationTables();
        } else if (arq.tabRSPB == 1) {
            generateSaturationTablesFromCorrelations();
        }
        applyDensityCorrectionsAndInletFluid();
        allocateEventProfileAndTrendArrays();
        resetCouplingAndOutputState();
    } catch (exception &excInt) {
        cout << "EXCECAO INESPERADA: " << excInt.what() << endl;
        // incluir falha
        logger.log(LOGGER_FALHA, LOG_ERR_UNEXPECTED_EXCEPTION, "", "", excInt.what());
        // gravar arquivo de log
        logger.writeOutputLog();
        // encerrar a aplicacao
        exit(EXIT_SUCCESS);
    }
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

void sisprod::gaslift::GasLiftTemperatureUpdater::dischargeTemperature(
    int cellIndex) const {
    system.tempDescarga(cellIndex);
}

void sisprod::gaslift::GasLiftTemperatureUpdater::gasTemperature(
    int cellIndex, double previousTemperature, int steadyMode) const {
    system.calctempGas(cellIndex, previousTemperature, steadyMode);
}

double sisprod::gaslift::GasLiftTemperatureUpdater::gasLiftDischargeTemperature(
    int valveIndex) const {
    return system.TempDescGL(valveIndex);
}

// ------------------------------------ steady-state callbacks (T086) ----
//
// Sixteen forwards, fourteen of which reach code that already lives in an
// extracted module. They come back through SProd anyway, because reaching
// sisprod::gaslift or sisprod::thermal needs one of THEIR state structs, and
// the only place that knows how to build those is this file.
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
// The two searches a march calls. Their results were discarded at the call
// site in the original and they are discarded here; see the note on the
// declarations.
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

void sisprod::transient::TransientStepUpdaters::generateFluidMiniTable() const {
    system.geraMiniTabFlu();
}
void sisprod::transient::TransientStepUpdaters::advanceGasSubStep() const {
    system.subtempoGas();
}

void sisprod::transient::TransientSolveUpdaters::solveHydrateEnvelopes() const {
    system.solveHydrateEnvelopes();
}
double sisprod::transient::TransientSolveUpdaters::findInjectionPressureDownstream() const {
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
void sisprod::transient::TransientSolveUpdaters::saveSources() const {
    system.salvaFonte();
}
void sisprod::transient::TransientSolveUpdaters::solveGasLine() const {
    system.solveLinGas();
}



namespace {

sisprod::gaslift::GasLiftState gasLiftStateOf(SProd &system) {
    return sisprod::gaslift::GasLiftState{
        .gasCells = system.celulaG,
        .cells = system.celula,
        .input = system.arq,
        .globals = system.vg1dSP,
        .gasCellCount = system.ncelGas,
        .lastCell = system.ncel,
        .gasLiftChokes = system.chokeVGL,
        .injectionChoke = system.chokeInj,
        .gasValveCellIndices = system.posicVGLG,
        .productionValveCellIndices = system.posicVGLP,
        .gasSystemMatrix = system.matglobG,
        .gasFreeTerms = system.termolivreG,
        .annulusTubingStart = system.ColunaAnulaIni,
        .annulusTubingEnd = system.ColunaAnulaFim,
        .tubingAnnulusStart = system.AnulaColunaIni,
        .tubingAnnulusEnd = system.AnulaColunaFim,
        .steadyIteration = system.iterperm,
        .networkCoupled = system.verificaAcop,
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
        .meanUnloadingFlowRate = system.vazmedDesc,
        .meanUnloadingTemperature = system.tempmedDEsc,
        .maximumMeanUnloadingFlowRates = system.vazmaxMedDesc,
        .unloadingTimeSteps = system.dtDesc,
        .continuousMeanUnloadingTemperature = system.tempMedContDesc,
        .maximumContinuousUnloadingCount = system.maxVecContDesc,
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
        .gasValveCellIndices = system.posicVGLG,
        .productionValveCellIndices = system.posicVGLP,
        .steadyIteration = system.iterperm,
        .searchOrigin = system.buscaIni,
        .convergenceMonitor = system.monitConvPerm,
        .baseConvergenceMonitor = system.monitConvPermBase,
        .annulusDrift = system.derivaAnel,
        .networkCoupled = system.verificaAcop,
        .endNode = system.noextremo,
        .thermalSourceDisabled = system.semTermo,
        .slowHeatTransfer = system.trocaTermicaLenta,
        .gasSurfacePressure = system.pGSup,
        .initialGasPressure = system.presiniG,
        .initialGasTemperature = system.tempiniG,
        .finalPressure = system.presfim,
        .timeStep = system.dt,
        .ambientTemperature = system.temperatura,
        .casingTemperature = system.tempRev,
        .inletQuality = system.titE,
        .productionFluidCount = system.nfluP,
        .updaters = {system},
    };
}

/// The state a boundary-condition search reads.
///
/// It composes the march state rather than rebuilding it, which is the whole
/// reason SteadyStateSearchState has a `march` member: of the 31 SProd members
/// the two halves touch, 17 are read by both, and two spellings of one fact
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
///
/// Seventy-two fields, generated from the same measurement that generated the
/// header -- what the 22 routines of this stage actually touch -- rather than
/// written twice. The ORDER is the header's declaration order, because C++20
/// designated initialisers must follow it; the first attempt sorted them by
/// name and the compiler rejected every field after the first.
sisprod::transient::TransientStepState transientStateOf(SProd &system) {
    return sisprod::transient::TransientStepState{
        .meanMaximumTimeStep = system.DTMaxMed,
        .meanMaximumPressureChange = system.DpMaxMed,
        .masterState = system.EstadoMaster1,
        .open = system.aberto,
        .initiallyOpen = system.abertoini,
        .timeChanged = system.alteraTempo,
        .inletCompletionFraction = system.betaE,
        .interfaceCell = system.celInter,
        .masterCounter = system.contaMaster1,
        .gasSpecificHeatTable = system.cpg,
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
        .finalPressure = system.presfim,
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

/// The state SolveTrans reads: the step state, composed rather than rebuilt,
/// plus the 52 members only the solve touches. Initialiser order is the
/// header's declaration order, which C++20 requires.
sisprod::transient::TransientSolveState transientSolveStateOf(SProd &system) {
    return sisprod::transient::TransientSolveState{
        .step = transientStateOf(system),
        .ambientTemperature = system.temperatura,
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
        .networkCoupled = system.verificaAcop,
        .poissonSolver3D = system.poisson3D,
        .pressureHistory = system.presVet,
        .printCounter = system.KontaImprime,
        .printTimeCounter = system.kimpT,
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



sisprod::thermal::ThermalState thermalStateOf(SProd &system) {
    return sisprod::thermal::ThermalState{
        .cells = system.celula,
        .gasCells = system.celulaG,
        .input = system.arq,
        .latentHeatTable = system.HLat,
        .globals = system.vg1dSP,
        .thermalSourceDisabled = system.semTermo,
        .productionNetworkCoupled = system.verificaAcopRedeS,
        .primarySectionStart = system.SecPrimIniRedeP,
        .primarySectionEnd = system.SecPrimFimRedeP,
        .coupledCellIndices = system.acertaIndAcop,
        .poissonSolver = system.poisson3D,
        .lastCell = system.ncel,
        .surfaceChoke = system.chokeSup,
        .surfaceChokeMassCondition = system.masChkSup,
        .networkEndpoint = system.noextremo,
        .gasSurfaceTemperature = system.tGSup,
        .latentHeatEnabled = system.CalcLat,
        .sourceUpdater = {.system = system},
        .completeModel = system.modeloCompleto,
        .massTransferModel = system.TransMassModel,
        .closureUpdater = {.system = system},
        .inletPressure = system.presE,
        .inletTemperature = system.tempE,
        .inletMassFraction = system.titE,
        .inletVoidFraction = system.alfE,
        .inletComposition = system.betaE,
        .evolutionUpdater = {.system = system},
        .surfaceChokeOpen = system.aberto,
        .defaultInletTemperature = system.temperatura,
        .timeStep = system.dt,
        .minimumCycleTimeStep = system.dtCicMin,
        .poisson2DCellIndices = system.indCelPoisson2D,
        .poisson2DCellCount = system.nCelulaPoisson2D,
        .steadyIteration = system.iterperm,
        .slowHeatTransferThreshold = system.trocaTermicaLenta,
        .annulusTubingStart = system.ColunaAnulaIni,
        .annulusTubingEnd = system.ColunaAnulaFim,
        .tubingAnnulusStart = system.AnulaColunaIni,
        .productionNetworkHeatCoupled = system.verificaAcopRedeP,
        .primaryNetworkSectionEnd = system.PrimSecFimRedeP,
        .primaryNetworkSectionStart = system.PrimSecIniRedeP,
        .gasCellCount = system.ncelGas,
        .gasLiftChokes = system.chokeVGL,
        .gasSurfacePressure = system.pGSup,
        .outletPressure = system.presfim,
        .surfaceTemperature = system.tempSup,
        .networkCoupled = system.verificaAcop,
        .tubingAnnulusEnd = system.AnulaColunaFim,
    };
}

}  // namespace

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

        double hidroM = sin(celula[ind].duto.teta) * (0.5 * celula[ind].dx) * (rholmix * (1 - alfE) + alfE * celula[ind].flui.MasEspGas(celula[ind].pres, celula[ind].temp)) * 9.82 / 98600.;
        double hidroJ = sin(celula[ind + 1].duto.teta) * (0.5 * celula[ind + 1].dx) * (rholmixJ * (1 - alfJ) + alfJ * celula[ind + 1].flui.MasEspGas(celula[ind + 1].pres, celula[ind + 1].temp)) * 9.82 / 98600.;

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
    for (int i = 0; i <= ncel; i++) {
        celula[i].fontemassLRini = celula[i].fontemassLR;
        celula[i].fontemassCRini = celula[i].fontemassCR;
        celula[i].fontemassGRini = celula[i].fontemassGR;
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
        if (celula[ind - 1].acsr.tipo != 5 && celula[ind - 1].acsr.tipo != 8 && (ind < ncel || celula[ncel].acsr.tipo == 3)) {
            celula[ind].fontemassLR = 0.;
            celula[ind].fontemassCR = 0.;
            celula[ind].fontemassGR = 0.;
        }
    }

    double agua_consumida_Mw = 0.;
    double gas_consumido_Mg = 0.;

    if (arq.calculaEnvelope == 1 && arq.tipoHmodel == 2 && (*vg1dSP).lixo5 > 0.01) {

        agua_consumida_Mw = celula[ind].agua_consumida_massa_step;

        gas_consumido_Mg = celula[ind].gas_consumido_massa_step;

	  //Atualizar BSW
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

  } //Alteracao Hidratos
  
  if (arq.calculaEnvelope==1 && arq.tipoHmodel==3 && (*vg1dSP).lixo5>0.01) { //alteracao Hidratos

	    agua_consumida_Mw  = celula[ind].agua_consumida_massa_step;

	    gas_consumido_Mg   = celula[ind].gas_consumido_massa_step;

	  //Atualizar BSW
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

  } //Alteracao Hidratos
    if (celula[ind].acsr.tipo == 1) {
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
                celula[ind].fontemassCR += rcomp * celula[ind].acsr.injg.razCompGas * celula[ind].acsr.injg.QGas / 86400.;
            }
        } else {
            celula[ind].fontemassGR += 0.;
            celula[ind].fontemassLR = 0.;
            celula[ind].fontemassCR = celula[ind].acsr.injg.QGas;
        }
    }
    if (celula[ind].acsr.tipo == 2) {
        double rlcA = celula[ind].acsr.injl.fluidocol.MasEspFlu(1.001, 15.);
        celula[ind].fontemassCR += rlcA * celula[ind].acsr.injl.QLiq * celula[ind].acsr.injl.bet / 86400;
        double massic = celula[ind].acsr.injl.QLiq * (1. - celula[ind].acsr.injl.bet) / 86400;
        double Rhogs = celula[ind].acsr.injl.FluidoPro.Deng * 1.225; // cel[ind].acsr.injl.FluidoPro.MasEspGas(1, 15);
        double Rhols = (1000 * 141.5 / (131.5 + celula[ind].acsr.injl.FluidoPro.API)) * (1 - celula[ind].acsr.injl.FluidoPro.BSW) + 1000. * celula[ind].acsr.injl.FluidoPro.Denag * celula[ind].acsr.injl.FluidoPro.BSW;
        double multiplicador = (Rhols + celula[ind].acsr.injl.FluidoPro.RGO * Rhogs * (1 - celula[ind].acsr.injl.FluidoPro.BSW));
        massic *= multiplicador;
        double fracmasshidra = celula[ind].acsr.injl.FluidoPro.FracMassHidra(pr, tr);
        celula[ind].fontemassLR += (1. - fracmasshidra) * massic;
        celula[ind].fontemassGR += fracmasshidra * massic;
    }
    if (celula[ind].acsr.tipo == 3) {
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
    if (celula[ind].acsr.tipo == 5) {
        celula[ind + 1].fontemassLR = 0.;
        celula[ind + 1].fontemassCR = 0.;
        celula[ind + 1].fontemassGR = 0.;
        if (modoPerm == 0)
            FonteValv(ind);
    }
    if (celula[ind].acsr.tipo == 8) {
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
    if (celula[ind].acsr.tipo == 9) {
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
    if (celula[ind].acsr.tipo == 10) {
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
    if (celula[ind].acsr.tipo == 15) {
        celula[ind].acsr.radialPoro.pW.val[0] = celula[ind].pres;
        double rs;
        double bo;
        double ba;
        if (celula[ind].flui.RGO < 1e7) {
            rs = celula[ind + 1].flui.RS(celula[ind].pres, celula[ind].temp);
            bo = celula[ind + 1].flui.BOFunc(celula[ind].pres, celula[ind].temp, rs);
            ba = celula[ind + 1].flui.BAFunc(celula[ind].pres, celula[ind].temp);
            rs = rs * 6.29 / 35.31467;
        } else {
            bo = 1;
            rs = 0;
            ba = 0.;
        }
        // BSW in-situ da celula anterior, na marcha, a i-esima celula
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
    if (celula[ind].acsr.tipo == 16) {
        celula[ind].acsr.poroso2D.dados.pW.val[0] = celula[ind].pres;
        double rs;
        double bo;
        double ba;
        if (celula[ind].flui.RGO < 1e7) {
            rs = celula[ind + 1].flui.RS(celula[ind].pres, celula[ind].temp);
            bo = celula[ind + 1].flui.BOFunc(celula[ind].pres, celula[ind].temp, rs);
            ba = celula[ind + 1].flui.BAFunc(celula[ind].pres, celula[ind].temp);
            rs = rs * 6.29 / 35.31467;
        } else {
            bo = 1;
            rs = 0;
            ba = 0.;
        }
        // BSW in-situ da celula anterior, na marcha, a i-esima celula
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

    if (arq.calculaEnvelope == 1 && arq.tipoHmodel == 2 && celula[ind].flui.BSW > 1e-12 && (*vg1dSP).lixo5 > 0.01) {
        celula[ind].fontemassLR -= (agua_consumida_Mw / (*vg1dSP).lixo5);
        celula[ind].fontemassGR -= (gas_consumido_Mg / (*vg1dSP).lixo5);
    }
  	//alteracao hidrato 3
  	if (arq.calculaEnvelope==1 && arq.tipoHmodel==3 && celula[ind].flui.BSW>1e-14 && (*vg1dSP).lixo5>0.01) { //alteracao Hidratos
  		celula[ind].fontemassLR -= (agua_consumida_Mw / (*vg1dSP).lixo5);
  		celula[ind].fontemassGR -= (gas_consumido_Mg / (*vg1dSP).lixo5);
  	}
}

void SProd::renovaalbetini() {

    celula[0].alfini = celula[0].alf;
    celula[0].betini = celula[0].bet;
    celula[0].alfPigDini = celula[0].alfPigD;
    celula[0].betPigDini = celula[0].betPigD;
    celula[0].alfPigEini = celula[0].alfPigE;
    celula[0].betPigEini = celula[0].betPigE;
    celula[1].alfLini = celula[0].alfini;
    celula[1].betLini = celula[0].betini;

    if (arq.ConContEntrada == 0) {
        celula[0].alfLini = celula[0].alfini;
        celula[0].betLini = celula[0].betini;
    } else {
        celula[0].alfLini = titE;
        celula[0].betLini = betaE;
    }

    indpigPini = indpigP;
    indpigP = 0;

    for (int i = 0; i <= ncel; i++) {
        celula[i].estadoPigini = celula[i].estadoPig;
        celula[i].alfini = celula[i].alf;
        if (i > 0)
            celula[i - 1].alfRini = celula[i].alfini;
        if (i < ncel)
            celula[i + 1].alfLini = celula[i].alfini;
        celula[i].betini = celula[i].bet;
        if (i > 0)
            celula[i - 1].betRini = celula[i].betini;
        if (i < ncel)
            celula[i + 1].betLini = celula[i].betini;

        celula[i].alfPigDini = celula[i].alfPigD;
        celula[i].betPigDini = celula[i].betPigD;
        celula[i].alfPigEini = celula[i].alfPigE;
        celula[i].betPigEini = celula[i].betPigE;

        celula[i].alfPigERini = celula[i].alfPigER;
        celula[i].betIini = celula[i].betI;
        celula[i].betRIini = celula[i].betRI;
        celula[i].betLIini = celula[i].betLI;

        if (i > 0)
            celula[i - 1].alfPigER = celula[i].alfPigE;
        if (i == ncel)
            celula[i].alfPigER = celula[i].alf;
        celula[i].DelPig = 0.;
        celula[i].RazAreaPig = 0.;
        celula[i].cdpig = 1.;
        if (celula[i].estadoPig == 1) {
            indpigP++;
            int ipig = celula[i].indpig;
            celula[i].DelPig = arq.pig[ipig].delpres;
            celula[i].RazAreaPig = arq.pig[ipig].razarea;
            celula[i].cdpig = arq.pig[ipig].cdPig;
            double AC = celula[i].duto.area;
            double jL = (celula[i].QL + celula[i].QG) / AC;
            double jR;
            if (i < ncel)
                jR = (celula[i + 1].QL + celula[i + 1].QG) / AC;
            else
                jR = (celula[i].QL + celula[i].QG) / AC;
            celula[i].velPigini = celula[i].velPig;
            celula[i].velPig = jL * celula[i].razPig + jR * (1. - celula[i].razPig) - celula[i].VazaPig / AC;
            for (int j = 0; j < npig; j++) {
                if (receb[j] == i) {
                    celula[i].velPig = 0.;
                    celula[i].estadoPig = 0;
                    celula[i].razPig = 0.;
                    celula[i].razPigini = 0.;
                    celula[i].alfPigDini = celula[i].alf;
                    celula[i].betPigDini = celula[i].bet;
                    celula[i].alfPigEini = celula[i].alf;
                    celula[i].betPigEini = celula[i].bet;
                    celula[i].alfPigD = celula[i].alf;
                    celula[i].betPigD = celula[i].bet;
                    celula[i].alfPigE = celula[i].alf;
                    celula[i].betPigE = celula[i].bet;
                    if (i > 0)
                        celula[i - 1].alfPigER = celula[i].alfPigE;
                    if (i == ncel)
                        celula[i].alfPigER = celula[i].alf;

                    indpigP--;
                }
            }
        }
        celula[i].razPigini = celula[i].razPig;
    }
    celula[ncel].alfRini = celula[ncel].alfini;
    celula[ncel].betRini = celula[ncel].betini;
}

void SProd::renovaMasEsp() {

    celula[0].rpC = celula[0].flui.MasEspLiq(celula[0].pres, celula[0].temp);
    celula[0].rgC = celula[0].flui.MasEspGas(celula[0].pres, celula[0].temp /*,1*/);
    celula[0].rcC = celula[0].fluicol.MasEspFlu(celula[0].pres, celula[0].temp);
    celula[0].rpL = celula[0].rpC;
    celula[0].rgL = celula[0].rgC;
    celula[0].rcL = celula[0].rcC;

    celula[0].rpCi = celula[0].rpC;
    celula[0].rgCi = celula[0].rgC;
    celula[0].rcCi = celula[0].rcC;
    celula[0].rpLi = celula[0].rpCi;
    celula[0].rgLi = celula[0].rgCi;
    celula[0].rcLi = celula[0].rcCi;

    celula[0].mipC = celula[0].flui.ViscOleo(celula[0].pres, celula[0].temp);
    celula[0].migC = celula[0].flui.ViscGas(celula[0].pres, celula[0].temp /*,1*/);
    celula[0].micC = celula[0].fluicol.VisFlu(celula[0].pres, celula[0].temp);

#pragma omp parallel for num_threads((*vg1dSP).ntrd)
    for (int i = 1; i <= ncel; i++) {
        double p;
        double t;
        p = celula[i].pres;
        t = celula[i].temp;
        celula[i].rpC = celula[i].flui.MasEspLiq(p, t);
        celula[i].rgC = celula[i].flui.MasEspGas(p, t);
        celula[i].rcC = celula[i].fluicol.MasEspFlu(p, t);

        celula[i].mipC = celula[i].flui.ViscOleo(p, t);
        celula[i].migC = celula[i].flui.ViscGas(p, t);
        celula[i].micC = celula[i].fluicol.VisFlu(p, t);

        double tmed = celula[i - 1].temp;
        if (celula[i].VTemper < 0.)
            tmed = celula[i].temp;
        ProFlu flu;
        if (celula[i].QL < 0.)
            flu = celula[i].flui;
        else
            flu = celula[i - 1].flui;
        celula[i].rpCi = flu.MasEspLiq(celula[i].presaux, tmed);
        celula[i].rgCi = flu.MasEspGas(celula[i].presaux, tmed);
        celula[i].rcCi = celula[i].fluicol.MasEspFlu(celula[i].presaux, tmed);
    }
    for (int i = 1; i <= ncel; i++) {

        celula[i - 1].rpR = celula[i].rpC;
        celula[i - 1].rgR = celula[i].rgC;
        celula[i - 1].rcR = celula[i].rcC;
        celula[i].rpL = celula[i - 1].rpC;
        celula[i].rgL = celula[i - 1].rgC;
        celula[i].rcL = celula[i - 1].rcC;

        celula[i - 1].mipR = celula[i].mipC;
        celula[i - 1].migR = celula[i].migC;
        celula[i - 1].micR = celula[i].micC;

        celula[i - 1].rpRi = celula[i].rpCi;
        celula[i - 1].rgRi = celula[i].rgCi;
        celula[i - 1].rcRi = celula[i].rcCi;
        celula[i].rpLi = celula[i - 1].rpCi;
        celula[i].rgLi = celula[i - 1].rgCi;
        celula[i].rcLi = celula[i - 1].rcCi;
    }
    celula[ncel].rpR = celula[ncel].rpC;
    celula[ncel].rgR = celula[ncel].rgC;
    celula[ncel].rcR = celula[ncel].rcC;

    celula[ncel].mipR = celula[ncel].mipC;
    celula[ncel].migR = celula[ncel].migC;
    celula[ncel].micR = celula[ncel].micC;
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
        .inletColumnFraction = betaE,
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
        .inletColumnFraction = betaE,
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
        .inletColumnFraction = betaE,
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
        .inletColumnFraction = betaE,
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
    hidro = (1 - arq.MedSimpPresFront) * 9.82 * sin(celula[i].duto.teta) * rhomix * dx;
}

void SProd::auxMiniTab(ProFlu &flui) {
    ProFlu fluC;
    fluC = flui;
    fluC.atualizaPropCompStandard();
    if (fluC.dCalculatedBeta > 0. && fluC.dCalculatedBeta < 1.)
        fluC.atualizaPropComp(flui.miniTabDin.pmin, flui.miniTabDin.tmin,
                              fluC.dCalculatedBeta, fluC.oCalculatedLiqComposition,
                              fluC.oCalculatedVapComposition, arq.pocinjec);
    else
        fluC.atualizaPropComp(flui.miniTabDin.pmin, flui.miniTabDin.tmin, -1, NULL, NULL, arq.pocinjec);
    flui.miniTabDin.rholF[0][0] =
        fluC.MasEspoleo(flui.miniTabDin.pmin, flui.miniTabDin.tmin);
    flui.miniTabDin.rhogF[0][0] =
        fluC.MasEspGas(flui.miniTabDin.pmin, flui.miniTabDin.tmin);
    flui.miniTabDin.DrhogDpF[0][0] =
        fluC.drhodp(flui.miniTabDin.pmin, flui.miniTabDin.tmin);
    flui.miniTabDin.DrholDtF[0][0] =
        fluC.DrholDT(flui.miniTabDin.pmin, flui.miniTabDin.tmin);
    flui.miniTabDin.DrhogDtF[0][0] =
        fluC.drhodt(flui.miniTabDin.pmin, flui.miniTabDin.tmin);
    flui.miniTabDin.valBO[0][0] =
        fluC.BOFunc(flui.miniTabDin.pmin, flui.miniTabDin.tmin);
    flui.miniTabDin.valZ[0][0] =
        fluC.Zdran(flui.miniTabDin.pmin, flui.miniTabDin.tmin);
    flui.miniTabDin.valdZdT[0][0] =
        fluC.DZDT(flui.miniTabDin.pmin, flui.miniTabDin.tmin);
    flui.miniTabDin.valdZdP[0][0] =
        fluC.DZDP(flui.miniTabDin.pmin, flui.miniTabDin.tmin);
    flui.miniTabDin.tit[0][0] =
        fluC.FracMass(flui.miniTabDin.pmin, flui.miniTabDin.tmin);
    flui.miniTabDin.rs[0][0] =
        fluC.RS(flui.miniTabDin.pmin, flui.miniTabDin.tmin);
    flui.miniTabDin.cplF[0][0] =
        fluC.CalorLiq(flui.miniTabDin.pmin, flui.miniTabDin.tmin);
    flui.miniTabDin.cpgF[0][0] =
        fluC.CalorGas(flui.miniTabDin.pmin, flui.miniTabDin.tmin);
    flui.miniTabDin.HlF[0][0] =
        fluC.EntalpLiq(flui.miniTabDin.pmin, flui.miniTabDin.tmin);
    flui.miniTabDin.HgF[0][0] =
        fluC.EntalpGas(flui.miniTabDin.pmin, flui.miniTabDin.tmin);
    flui.miniTabDin.PBF[0] =
        fluC.PB(flui.miniTabDin.pmin, flui.miniTabDin.tmin);
    if (fluC.dCalculatedBeta > 0. && fluC.dCalculatedBeta < 1.)
        fluC.atualizaPropComp(flui.miniTabDin.pmin, flui.miniTabDin.tmax,
                              fluC.dCalculatedBeta, fluC.oCalculatedLiqComposition,
                              fluC.oCalculatedVapComposition, arq.pocinjec);
    else
        fluC.atualizaPropComp(flui.miniTabDin.pmin, flui.miniTabDin.tmax, -1, NULL, NULL, arq.pocinjec);
    flui.miniTabDin.rholF[0][1] =
        fluC.MasEspoleo(flui.miniTabDin.pmin, flui.miniTabDin.tmax);
    flui.miniTabDin.rhogF[0][1] =
        fluC.MasEspGas(flui.miniTabDin.pmin, flui.miniTabDin.tmax);
    flui.miniTabDin.DrhogDpF[0][1] =
        fluC.drhodp(flui.miniTabDin.pmin, flui.miniTabDin.tmax);
    flui.miniTabDin.DrholDtF[0][1] =
        fluC.DrholDT(flui.miniTabDin.pmin, flui.miniTabDin.tmax);
    flui.miniTabDin.DrhogDtF[0][1] =
        fluC.drhodt(flui.miniTabDin.pmin, flui.miniTabDin.tmax);
    flui.miniTabDin.valBO[0][1] =
        fluC.BOFunc(flui.miniTabDin.pmin, flui.miniTabDin.tmax);
    flui.miniTabDin.valZ[0][1] =
        fluC.Zdran(flui.miniTabDin.pmin, flui.miniTabDin.tmax);
    flui.miniTabDin.valdZdT[0][1] =
        fluC.DZDT(flui.miniTabDin.pmin, flui.miniTabDin.tmax);
    flui.miniTabDin.valdZdP[0][1] =
        fluC.DZDP(flui.miniTabDin.pmin, flui.miniTabDin.tmax);
    flui.miniTabDin.tit[0][1] =
        fluC.FracMass(flui.miniTabDin.pmin, flui.miniTabDin.tmax);
    flui.miniTabDin.rs[0][1] =
        fluC.RS(flui.miniTabDin.pmin, flui.miniTabDin.tmax);
    flui.miniTabDin.cplF[0][1] =
        fluC.CalorLiq(flui.miniTabDin.pmin, flui.miniTabDin.tmax);
    flui.miniTabDin.cpgF[0][1] =
        fluC.CalorGas(flui.miniTabDin.pmin, flui.miniTabDin.tmax);
    flui.miniTabDin.HlF[0][1] =
        fluC.EntalpLiq(flui.miniTabDin.pmin, flui.miniTabDin.tmax);
    flui.miniTabDin.HgF[0][1] =
        fluC.EntalpGas(flui.miniTabDin.pmin, flui.miniTabDin.tmax);
    flui.miniTabDin.PBF[1] =
        fluC.PB(flui.miniTabDin.pmin, flui.miniTabDin.tmax);
    if (fluC.dCalculatedBeta > 0. && fluC.dCalculatedBeta < 1.)
        fluC.atualizaPropComp(flui.miniTabDin.pmax, flui.miniTabDin.tmin,
                              fluC.dCalculatedBeta, fluC.oCalculatedLiqComposition,
                              fluC.oCalculatedVapComposition, arq.pocinjec);
    else
        fluC.atualizaPropComp(flui.miniTabDin.pmax, flui.miniTabDin.tmin, -1, NULL, NULL, arq.pocinjec);
    flui.miniTabDin.rholF[1][0] =
        fluC.MasEspoleo(flui.miniTabDin.pmax, flui.miniTabDin.tmin);
    flui.miniTabDin.rhogF[1][0] =
        fluC.MasEspGas(flui.miniTabDin.pmax, flui.miniTabDin.tmin);
    flui.miniTabDin.DrhogDpF[1][0] =
        fluC.drhodp(flui.miniTabDin.pmax, flui.miniTabDin.tmin);
    flui.miniTabDin.DrholDtF[1][0] =
        fluC.DrholDT(flui.miniTabDin.pmax, flui.miniTabDin.tmin);
    flui.miniTabDin.DrhogDtF[1][0] =
        fluC.drhodt(flui.miniTabDin.pmax, flui.miniTabDin.tmin);
    flui.miniTabDin.valBO[1][0] =
        fluC.BOFunc(flui.miniTabDin.pmax, flui.miniTabDin.tmin);
    flui.miniTabDin.valZ[1][0] =
        fluC.Zdran(flui.miniTabDin.pmax, flui.miniTabDin.tmin);
    flui.miniTabDin.valdZdT[1][0] =
        fluC.DZDT(flui.miniTabDin.pmax, flui.miniTabDin.tmin);
    flui.miniTabDin.valdZdP[1][0] =
        fluC.DZDP(flui.miniTabDin.pmax, flui.miniTabDin.tmin);
    flui.miniTabDin.tit[1][0] =
        fluC.FracMass(flui.miniTabDin.pmax, flui.miniTabDin.tmin);
    flui.miniTabDin.rs[1][0] =
        fluC.RS(flui.miniTabDin.pmax, flui.miniTabDin.tmin);
    flui.miniTabDin.cplF[1][0] =
        fluC.CalorLiq(flui.miniTabDin.pmax, flui.miniTabDin.tmin);
    flui.miniTabDin.cpgF[1][0] =
        fluC.CalorGas(flui.miniTabDin.pmax, flui.miniTabDin.tmin);
    flui.miniTabDin.HlF[1][0] =
        fluC.EntalpLiq(flui.miniTabDin.pmax, flui.miniTabDin.tmin);
    flui.miniTabDin.HgF[1][0] =
        fluC.EntalpGas(flui.miniTabDin.pmax, flui.miniTabDin.tmin);
    if (fluC.dCalculatedBeta > 0. && fluC.dCalculatedBeta < 1.)
        fluC.atualizaPropComp(flui.miniTabDin.pmax, flui.miniTabDin.tmax,
                              fluC.dCalculatedBeta, fluC.oCalculatedLiqComposition,
                              fluC.oCalculatedVapComposition, arq.pocinjec);
    else
        fluC.atualizaPropComp(flui.miniTabDin.pmax, flui.miniTabDin.tmax, -1, NULL, NULL, arq.pocinjec);
    flui.miniTabDin.rholF[1][1] =
        fluC.MasEspoleo(flui.miniTabDin.pmax, flui.miniTabDin.tmax);
    flui.miniTabDin.rhogF[1][1] =
        fluC.MasEspGas(flui.miniTabDin.pmax, flui.miniTabDin.tmax);
    flui.miniTabDin.DrhogDpF[1][1] =
        fluC.drhodp(flui.miniTabDin.pmax, flui.miniTabDin.tmax);
    flui.miniTabDin.DrholDtF[1][1] =
        fluC.DrholDT(flui.miniTabDin.pmax, flui.miniTabDin.tmax);
    flui.miniTabDin.DrhogDtF[1][1] =
        fluC.drhodt(flui.miniTabDin.pmax, flui.miniTabDin.tmax);
    flui.miniTabDin.valBO[1][1] =
        fluC.BOFunc(flui.miniTabDin.pmax, flui.miniTabDin.tmax);
    flui.miniTabDin.valZ[1][1] =
        fluC.Zdran(flui.miniTabDin.pmax, flui.miniTabDin.tmax);
    flui.miniTabDin.valdZdT[1][1] =
        fluC.DZDT(flui.miniTabDin.pmax, flui.miniTabDin.tmax);
    flui.miniTabDin.valdZdP[1][1] =
        fluC.DZDP(flui.miniTabDin.pmax, flui.miniTabDin.tmax);
    flui.miniTabDin.tit[1][1] =
        fluC.FracMass(flui.miniTabDin.pmax, flui.miniTabDin.tmax);
    flui.miniTabDin.rs[1][1] =
        fluC.RS(flui.miniTabDin.pmax, flui.miniTabDin.tmax);
    flui.miniTabDin.cplF[1][1] =
        fluC.CalorLiq(flui.miniTabDin.pmax, flui.miniTabDin.tmax);
    flui.miniTabDin.cpgF[1][1] =
        fluC.CalorGas(flui.miniTabDin.pmax, flui.miniTabDin.tmax);
    flui.miniTabDin.HlF[1][1] =
        fluC.EntalpLiq(flui.miniTabDin.pmax, flui.miniTabDin.tmax);
    flui.miniTabDin.HgF[1][1] =
        fluC.EntalpGas(flui.miniTabDin.pmax, flui.miniTabDin.tmax);


    std::pair<double, int> titVec[4];

    for(int j=0;j<2;j++){
    	for(int k=0;k<2;k++){
    		titVec[2*j+k]={flui.miniTabDin.tit[j][k],2*j+k};
    	}
    }
    std::sort(titVec, titVec + 4);
    if(titVec[0].first<1e-3){
    	int busca=1;
    	while(busca<4 && titVec[busca].first<1e-3)busca++;
    	if(busca<4){
    		int jtroca;
    		int ktroca;
    		if(titVec[busca].second==0){
    			jtroca=0;
    			ktroca=0;
    		}
    		else if(titVec[busca].second==1){
    			jtroca=0;
    			ktroca=1;
    		}
    		else if(titVec[busca].second==2){
    			jtroca=1;
    			ktroca=0;
    		}
    		else if(titVec[busca].second==3){
    			jtroca=1;
    			ktroca=1;
    		}
    	    for(int j=0;j<2;j++){
    	    	for(int k=0;k<2;k++){
    	    		if(flui.miniTabDin.tit[j][k]<1e-3){
    	    			flui.miniTabDin.rhogF[j][k]=flui.miniTabDin.rhogF[jtroca][ktroca];
    	    			flui.miniTabDin.DrhogDpF[j][k]=flui.miniTabDin.DrhogDpF[jtroca][ktroca];
    	    			flui.miniTabDin.DrhogDtF[j][k]=flui.miniTabDin.DrhogDtF[jtroca][ktroca];
    	    			flui.miniTabDin.valZ[j][k]=flui.miniTabDin.valZ[jtroca][ktroca];
    	    			flui.miniTabDin.valdZdT[j][k]=flui.miniTabDin.valdZdT[jtroca][ktroca];
    	    			flui.miniTabDin.valdZdP[j][k]=flui.miniTabDin.valdZdP[jtroca][ktroca];
    	    			flui.miniTabDin.cpgF[j][k]=flui.miniTabDin.cpgF[jtroca][ktroca];
    	    			flui.miniTabDin.HgF[j][k]=flui.miniTabDin.HgF[jtroca][ktroca];
    	    		}
    	    	}
    	    }
    	}
    }
    if(titVec[3].first>1.-1e-3){
    	int busca=2;
    	while(busca>=0 && titVec[busca].first>1.+1e-3)busca--;
    	if(busca>=0){
    		int jtroca;
    		int ktroca;
    		if(titVec[busca].second==0){
    			jtroca=0;
    			ktroca=0;
    		}
    		else if(titVec[busca].second==1){
    			jtroca=0;
    			ktroca=1;
    		}
    		else if(titVec[busca].second==2){
    			jtroca=1;
    			ktroca=0;
    		}
    		else if(titVec[busca].second==3){
    			jtroca=1;
    			ktroca=1;
    		}
    	    for(int j=0;j<2;j++){
    	    	for(int k=0;k<2;k++){
    	    		if(flui.miniTabDin.tit[j][k]>1.-1e-3){
    	    			flui.miniTabDin.rholF[j][k]=flui.miniTabDin.rholF[jtroca][ktroca];
    	    			flui.miniTabDin.valBO[j][k]=flui.miniTabDin.valBO[jtroca][ktroca];
    	    			flui.miniTabDin.DrholDtF[j][k]=flui.miniTabDin.DrholDtF[jtroca][ktroca];
    	    			flui.miniTabDin.rs[j][k]=flui.miniTabDin.rs[jtroca][ktroca];
    	    			flui.miniTabDin.cplF[j][k]=flui.miniTabDin.cplF[jtroca][ktroca];
    	    			flui.miniTabDin.HlF[j][k]=flui.miniTabDin.HlF[jtroca][ktroca];
    	    		}
    	    	}
    	    }
    	}
    }
}

void SProd::geraMiniTabFlu() {
    (*vg1dSP).modoTransiente = 0;
#pragma omp parallel for num_threads((*vg1dSP).ntrd)
    for (int i = 0; i <= ncel; i++) {
        double delp;
        double delt;
        delp = 0.5 * celula[i].pres;
        if (delp > arq.miniTabDp)
            delp = arq.miniTabDp;
        if (delp < 5) {
            celula[i].flui.miniTabDin.pmax = celula[i].pres + 5.;
            celula[i].flui.miniTabDin.pmin = celula[i].pres - delp;
            if (celula[i].flui.miniTabDin.pmin < 0.9)
                celula[i].flui.miniTabDin.pmin = 0.9;
        } else {
            celula[i].flui.miniTabDin.pmax = celula[i].pres + delp;
            celula[i].flui.miniTabDin.pmin = celula[i].pres - delp;
        }
        delt = arq.miniTabDt;
        celula[i].flui.miniTabDin.tmax = celula[i].temp + delt;
        celula[i].flui.miniTabDin.tmin = celula[i].temp - delt;
        if(arq.miniTabAtraso > 0)auxMiniTab(celula[i].flui);
        if (celula[i].acsr.tipo == 1) {
            celula[i].acsr.injg.FluidoPro.miniTabDin.pmax = celula[i].flui.miniTabDin.pmax;
            celula[i].acsr.injg.FluidoPro.miniTabDin.pmin = celula[i].flui.miniTabDin.pmin;
            celula[i].acsr.injg.FluidoPro.miniTabDin.tmax = celula[i].flui.miniTabDin.tmax;
            celula[i].acsr.injg.FluidoPro.miniTabDin.tmin = celula[i].flui.miniTabDin.tmin;
            celula[i].acsr.injg.FluidoPro.atualizaPropCompStandard();
            if (celula[i].acsr.injg.FluidoPro.dCalculatedBeta > 0. && celula[i].acsr.injg.FluidoPro.dCalculatedBeta < 1.)
                celula[i].acsr.injg.FluidoPro.atualizaPropComp(celula[i].pres, celula[i].temp,
                                                               celula[i].acsr.injg.FluidoPro.dCalculatedBeta, celula[i].acsr.injg.FluidoPro.oCalculatedLiqComposition,
                                                               celula[i].acsr.injg.FluidoPro.oCalculatedVapComposition, arq.pocinjec);
            else
                celula[i].acsr.injg.FluidoPro.atualizaPropComp(celula[i].pres, celula[i].temp, -1, NULL, NULL, arq.pocinjec);
            if(arq.miniTabAtraso > 0)auxMiniTab(celula[i].acsr.injg.FluidoPro);
        } else if (celula[i].acsr.tipo == 2) {
            celula[i].acsr.injl.FluidoPro.miniTabDin.pmax = celula[i].flui.miniTabDin.pmax;
            celula[i].acsr.injl.FluidoPro.miniTabDin.pmin = celula[i].flui.miniTabDin.pmin;
            celula[i].acsr.injl.FluidoPro.miniTabDin.tmax = celula[i].flui.miniTabDin.tmax;
            celula[i].acsr.injl.FluidoPro.miniTabDin.tmin = celula[i].flui.miniTabDin.tmin;
            celula[i].acsr.injl.FluidoPro.atualizaPropCompStandard();
            if (celula[i].acsr.injl.FluidoPro.dCalculatedBeta > 0. && celula[i].acsr.injl.FluidoPro.dCalculatedBeta < 1.)
                celula[i].acsr.injl.FluidoPro.atualizaPropComp(celula[i].pres, celula[i].temp,
                                                               celula[i].acsr.injl.FluidoPro.dCalculatedBeta, celula[i].acsr.injl.FluidoPro.oCalculatedLiqComposition,
                                                               celula[i].acsr.injl.FluidoPro.oCalculatedVapComposition, arq.pocinjec);
            else
                celula[i].acsr.injl.FluidoPro.atualizaPropComp(celula[i].pres, celula[i].temp, -1, NULL, NULL, arq.pocinjec);
            if(arq.miniTabAtraso > 0)auxMiniTab(celula[i].acsr.injl.FluidoPro);
        } else if (celula[i].acsr.tipo == 3) {
            celula[i].acsr.ipr.FluidoPro.miniTabDin.pmax = celula[i].flui.miniTabDin.pmax;
            celula[i].acsr.ipr.FluidoPro.miniTabDin.pmin = celula[i].flui.miniTabDin.pmin;
            celula[i].acsr.ipr.FluidoPro.miniTabDin.tmax = celula[i].flui.miniTabDin.tmax;
            celula[i].acsr.ipr.FluidoPro.miniTabDin.tmin = celula[i].flui.miniTabDin.tmin;
            celula[i].acsr.ipr.FluidoPro.atualizaPropCompStandard();
            if (celula[i].acsr.ipr.FluidoPro.dCalculatedBeta > 0. && celula[i].acsr.ipr.FluidoPro.dCalculatedBeta < 1.)
                celula[i].acsr.ipr.FluidoPro.atualizaPropComp(celula[i].pres, celula[i].temp,
                                                              celula[i].acsr.ipr.FluidoPro.dCalculatedBeta, celula[i].acsr.ipr.FluidoPro.oCalculatedLiqComposition,
                                                              celula[i].acsr.ipr.FluidoPro.oCalculatedVapComposition, arq.pocinjec);
            else
                celula[i].acsr.ipr.FluidoPro.atualizaPropComp(celula[i].pres, celula[i].temp, -1, NULL, NULL, arq.pocinjec);
            if(arq.miniTabAtraso > 0)auxMiniTab(celula[i].acsr.ipr.FluidoPro);
        } else if (celula[i].acsr.tipo == 15) {
        	if(arq.miniTabAtraso > 0)celula[i].acsr.radialPoro.geraMiniTabFlu();
        } else if (celula[i].acsr.tipo == 16) {
        	if(arq.miniTabAtraso > 0)celula[i].acsr.poroso2D.geraMiniTabFlu();
        } else if (celula[i].acsr.tipo == 9) {
            celula[i].acsr.fontechk.fluidoP.miniTabDin.pmax = celula[i].flui.miniTabDin.pmax;
            celula[i].acsr.fontechk.fluidoP.miniTabDin.pmin = celula[i].flui.miniTabDin.pmin;
            celula[i].acsr.fontechk.fluidoP.miniTabDin.tmax = celula[i].flui.miniTabDin.tmax;
            celula[i].acsr.fontechk.fluidoP.miniTabDin.tmin = celula[i].flui.miniTabDin.tmin;
            celula[i].acsr.fontechk.fluidoP.atualizaPropCompStandard();
            if (celula[i].acsr.fontechk.fluidoP.dCalculatedBeta > 0. && celula[i].acsr.fontechk.fluidoP.dCalculatedBeta < 1.)
                celula[i].acsr.fontechk.fluidoP.atualizaPropComp(celula[i].pres, celula[i].temp,
                                                                 celula[i].acsr.fontechk.fluidoP.dCalculatedBeta, celula[i].acsr.fontechk.fluidoP.oCalculatedLiqComposition,
                                                                 celula[i].acsr.fontechk.fluidoP.oCalculatedVapComposition, arq.pocinjec);
            else
                celula[i].acsr.fontechk.fluidoP.atualizaPropComp(celula[i].pres, celula[i].temp, -1, NULL, NULL, arq.pocinjec);
            if(arq.miniTabAtraso > 0)auxMiniTab(celula[i].acsr.fontechk.fluidoP);
        } else if (celula[i].acsr.tipo == 10) {
            celula[i].acsr.injm.FluidoPro.miniTabDin.pmax = celula[i].flui.miniTabDin.pmax;
            celula[i].acsr.injm.FluidoPro.miniTabDin.pmin = celula[i].flui.miniTabDin.pmin;
            celula[i].acsr.injm.FluidoPro.miniTabDin.tmax = celula[i].flui.miniTabDin.tmax;
            celula[i].acsr.injm.FluidoPro.miniTabDin.tmin = celula[i].flui.miniTabDin.tmin;
            celula[i].acsr.injm.FluidoPro.atualizaPropCompStandard();
            if (celula[i].acsr.injm.FluidoPro.dCalculatedBeta > 0. && celula[i].acsr.injm.FluidoPro.dCalculatedBeta < 1.)
                celula[i].acsr.injm.FluidoPro.atualizaPropComp(celula[i].pres, celula[i].temp,
                                                               celula[i].acsr.injm.FluidoPro.dCalculatedBeta, celula[i].acsr.injm.FluidoPro.oCalculatedLiqComposition,
                                                               celula[i].acsr.injm.FluidoPro.oCalculatedVapComposition, arq.pocinjec);
            else
                celula[i].acsr.injm.FluidoPro.atualizaPropComp(celula[i].pres, celula[i].temp, -1, NULL, NULL, arq.pocinjec);
            if(arq.miniTabAtraso > 0)auxMiniTab(celula[i].acsr.injm.FluidoPro);
        }
    }
    (*vg1dSP).modoTransiente = 1;
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
    for (int i = 0; i <= ncel; i++) {
        if (i < ncel)
            celula[i].WaxDeposition(arq.detalParafina, ncel);
        else {
            if (masChkSup == 0 || celula[ncel].Mliqini > 0) {
                double cpDep = celula[i - 1].duto.cp[0];
                double kDep = celula[i - 1].duto.cond[0];
                double rhoDep = celula[i - 1].duto.rhoC[0];
                celula[i].deltaPar = celula[i - 1].deltaPar;
                if (celula[i].parafinado == 0 && celula[i].deltaPar > 0.) {
                    celula[i].duto.atualizaCamada(celula[i].deltaPar, arq.detalParafina.rug, cpDep, kDep, rhoDep);
                    celula[i].calor.atualiza(celula[i].duto, 1);
                    celula[i].parafinado = 1;
                } else if (celula[i].deltaPar > 0.) {
                    celula[i].duto.atualizaCamada2(celula[i].deltaPar, cpDep, kDep, rhoDep);
                    celula[i].calor.atualiza2(celula[i].duto);
                }
            } else
                celula[i].WaxDeposition(arq.detalParafina, ncel);
        }
    }
}

void SProd::renovaRGOdgYco2(ProFlu fluiRev) {
    Vcr<double> rgo(ncel);
    Vcr<double> dg(ncel);
    Vcr<double> yco2(ncel);
    Vcr<double> API(ncel);
    Vcr<double> BSW(ncel);
    Vcr<double> denag(ncel);
    Vcr<double> VISCL(ncel);
    Vcr<double> VISCH(ncel);
    double dt = celula[1].dt;
    double tL;
    if (arq.flashCompleto == 1 && (arq.tabent.tmin - 0) > (*vg1dSP).localtiny)
        tL = arq.tabent.tmin + 0.1;
    else
        tL = 0.;
    double tH;
    if (arq.flashCompleto == 1 && (70 - arq.tabent.tmax) > (*vg1dSP).localtiny)
        tH = arq.tabent.tmax - 0.1;
    else
        tH = 70.;
    int imin = 1;
    if ((*vg1dSP).chaverede != 0 && (celula[0].acsr.tipo == 10 || arq.ConContEntrada > 0))
        imin = 0;
    if (celula[0].acsr.tipo == 15 && (celula[0].fontemassCR + celula[0].fontemassGR + celula[0].fontemassLR) > 0.)
        celula[0].flui.BSW = celula[0].acsr.radialPoro.BSW;
    else if (celula[0].acsr.tipo == 16 && (celula[0].fontemassCR + celula[0].fontemassGR + celula[0].fontemassLR) > 0.)
        celula[0].flui.BSW = celula[0].acsr.poroso2D.dados.transfer.BSW;
#pragma omp parallel for num_threads((*vg1dSP).ntrd)
    for (int i = imin; i < ncel; i++) {
        double MultOe;
        double MultOd;
        double at = celula[i].duto.area;
        double dx = celula[i].dx;
        double ti0;
        if (celula[i].VTemper < 0.)
            ti0 = celula[i].temp;
        else {
            if (i > 0)
                ti0 = celula[i - 1].temp;
            else if (arq.ConContEntrada == 1)
                ti0 = tempE;
            else
                ti0 = celula[i].temp;
        }
        double ti1 = celula[i].temp;
        if (celula[i + 1].VTemper < 0.)
            ti1 = celula[i + 1].temp;

        double rgo0;
        double betI0;
        double bo0;
        double ba0;
        double bsw0;
        double rs0;
        double dg0O;
        double yco20O;
        double API0;
        double BSW0;
        double denag0;
        double viscL0;
        double viscH0;
        double razdgd0;
        double razdgl0;
        if (i > 0 || arq.ConContEntrada == 0) {
            if (i > 0 && celula[i].QG >= 0.)
                betI0 = celula[i - 1].betPigD;
            else
                betI0 = celula[i].betPigE;
        } else {
            if (celula[i].QG >= 0.)
                betI0 = betaE;
            else
                betI0 = celula[i].betPigE;
        }
        if ((i > 0 || arq.ConContEntrada == 1) && celula[i].QL >= 0.) {
            double pe;
            double te;
            if (i > 0) {
                pe = celula[i - 1].pres;
                te = celula[i - 1].temp;
            } else {
                pe = presE;
                te = tempE;
            }
            rgo0 = (*celula[i].fluiL).RGO;
            if (arq.ConContEntrada == 0)
                betI0 = celula[i - 1].betPigD; // testeBeta
            else
                betI0 = betaE; // testeBeta
            rs0 = (*celula[i].fluiL).RS(pe, te);
            bo0 = (*celula[i].fluiL).BOFunc(pe, te, rs0);
            ba0 = (*celula[i].fluiL).BAFunc(pe, te);
            bsw0 = (*celula[i].fluiL).BSW * ba0 / (bo0 + ba0 * (*celula[i].fluiL).BSW - (*celula[i].fluiL).BSW * bo0);
            rs0 = rs0 * 6.29 / 35.31467;
            dg0O = (*celula[i].fluiL).Deng;
            razdgd0 = 1 / (*celula[i].fluiL).rDgD;
            razdgl0 = 1 / (*celula[i].fluiL).rDgL;
            yco20O = (*celula[i].fluiL).yco2;
            API0 = (*celula[i].fluiL).API;
            BSW0 = (*celula[i].fluiL).BSW;
            denag0 = (*celula[i].fluiL).Denag;
            viscL0 = 0 * 30 + 1 * (*celula[i].fluiL).VisOM(tL);
            viscH0 = 0 * 20 + 1 * (*celula[i].fluiL).VisOM(tH);
        } else {
            betI0 = celula[i].betPigE; // testebeta
            rgo0 = celula[i].flui.RGO;
            rs0 = celula[i].flui.RS(celula[i].pres, celula[i].temp);
            bo0 = celula[i].flui.BOFunc(celula[i].pres, celula[i].temp, rs0);
            ba0 = celula[i].flui.BAFunc(celula[i].pres, celula[i].temp);
            bsw0 = celula[i].flui.BSW * ba0 / (bo0 + ba0 * celula[i].flui.BSW - celula[i].flui.BSW * bo0);
            rs0 = rs0 * 6.29 / 35.31467;
            dg0O = celula[i].flui.Deng;
            razdgd0 = 1 / celula[i].flui.rDgD;
            razdgl0 = 1 / celula[i].flui.rDgL;
            yco20O = celula[i].flui.yco2;
            API0 = celula[i].flui.API;
            BSW0 = celula[i].flui.BSW;
            denag0 = celula[i].flui.Denag;
            viscL0 = 0 * 30 + 1 * celula[i].flui.VisOM(tL);
            viscH0 = 0 * 20 + 1 * celula[i].flui.VisOM(tH);
        }
        if (bo0 < 1e-15)
            bo0 = 1e-15;

        double rgo1 = celula[i].flui.RGO;
        double betI1 = celula[i].betPigD;
        double rs1 = celula[i].flui.RS(celula[i].pres, celula[i].temp);
        double bo1 = celula[i].flui.BOFunc(celula[i].pres, celula[i].temp, rs1);
        double ba1 = celula[i].flui.BAFunc(celula[i].pres, celula[i].temp);
        double bsw1 = celula[i].flui.BSW * ba1 / (bo1 + ba1 * celula[i].flui.BSW - celula[i].flui.BSW * bo1);
        rs1 = rs1 * 6.29 / 35.31467;
        double dg1O = celula[i].flui.Deng;
        double yco21O = celula[i].flui.yco2;
        double API1 = celula[i].flui.API;
        double BSW1 = celula[i].flui.BSW;
        double denag1 = celula[i].flui.Denag;
        double viscL1 = 0 * 30 + 1 * celula[i].flui.VisOM(tL);
        double viscH1 = 0 * 20 + 1 * celula[i].flui.VisOM(tH);
        double razdgd1 = 1 / celula[i].flui.rDgD;
        double razdgl1 = 1 / celula[i].flui.rDgL;

        // betI1 = celula[i + 1].betPigE;    //duvidabeta
        if (celula[i + 1].QL < 0.) {
            betI1 = celula[i + 1].betPigE; // testeBeta
            rgo1 = celula[i + 1].flui.RGO;
            rs1 = celula[i + 1].flui.RS(celula[i + 1].pres, celula[i + 1].temp);
            bo1 = celula[i + 1].flui.BOFunc(celula[i + 1].pres, celula[i + 1].temp, rs1);
            ba1 = celula[i + 1].flui.BAFunc(celula[i + 1].pres, celula[i + 1].temp);
            bsw1 = celula[i + 1].flui.BSW * ba1 / (bo1 + ba1 * celula[i + 1].flui.BSW - celula[i + 1].flui.BSW * bo1);
            rs1 = rs1 * 6.29 / 35.31467;
            dg1O = celula[i + 1].flui.Deng;
            razdgd1 = 1 / celula[i + 1].flui.rDgD;
            razdgl1 = 1 / celula[i + 1].flui.rDgL;
            yco21O = celula[i + 1].flui.yco2;
            API1 = celula[i + 1].flui.API;
            BSW1 = celula[i + 1].flui.BSW;
            denag1 = celula[i + 1].flui.Denag;
            viscL1 = 0 * 30 + 1 * celula[i + 1].flui.VisOM(tL);
            viscH1 = 0 * 20 + 1 * celula[i + 1].flui.VisOM(tH);
        }
        if (bo1 < 1e-15)
            bo1 = 1e-15;

        double rhog0;
        double rhogST0;
        double dg0G;
        double yco20G;
        if ((i > 0 || arq.ConContEntrada == 1) && celula[i].QG > 0) {
            double pe;
            double te;
            if (i > 0) {
                pe = celula[i - 1].pres;
                te = celula[i - 1].temp;
            } else {
                pe = presE;
                te = tempE;
            }
            rhog0 = celula[i].rgL;
            rhogST0 = (*celula[i].fluiL).Deng * 1.225;
            dg0G = (*celula[i].fluiL).Deng;
            yco20G = (*celula[i].fluiL).yco2;
        } else {
            rhog0 = celula[i].rgC;
            rhogST0 = celula[i].flui.Deng * 1.225;
            dg0G = celula[i].flui.Deng;
            yco20G = celula[i].flui.yco2;
        }

        double rhog1 = celula[i].rgC;
        double rhogST1 = celula[i].flui.Deng * 1.225;
        double dg1G = celula[i].flui.Deng;
        double yco21G = celula[i].flui.yco2;
        if (celula[i + 1].QG <= 0.) {
            rhog1 = celula[i].rgR;
            rhogST1 = celula[i + 1].flui.Deng * 1.225;
            dg1G = celula[i + 1].flui.Deng;
            yco21G = celula[i + 1].flui.yco2;
        }

        if (i == 237) {
            int para;
            para = 0;
        }

        double hol = 1. - celula[i].alf;
        double bet = celula[i].bet;
        double rholST = (1 - celula[i].flui.BSW) * (1000 * 141.5 / (131.5 + celula[i].flui.API)) + celula[i].flui.BSW * 1000 * celula[i].flui.Denag;
        double rhog = celula[i].rgC;
        double rhogST = celula[i].flui.Deng * 1.225;
        double rs = celula[i].flui.RS(celula[i].pres, celula[i].temp);
        double razdgd = 1 / celula[i].flui.rDgD;
        double razdgl = 1 / celula[i].flui.rDgL;
        double bo = celula[i].flui.BOFunc(celula[i].pres, celula[i].temp, rs);
        double ba = celula[i].flui.BAFunc(celula[i].pres, celula[i].temp);
        double bsw = celula[i].flui.BSW * ba / (bo + ba * celula[i].flui.BSW - celula[i].flui.BSW * bo);
        rs = rs * 6.29 / 35.31467;

        double dgini = celula[i].flui.Deng;
        double yco2ini = celula[i].flui.yco2;
        double rgoini = celula[i].flui.RGO;
        double fonteO = celula[i].fontemassLR;
        double fonteG = celula[i].fontemassGR;
        double fonteP = celula[i].fontemassLR;
        double fonteA = celula[i].fontemassLR;
        double APIini = celula[i].flui.API;
        double BSWini = celula[i].flui.BSW;
        double denagini = celula[i].flui.Denag;
        double APIF = APIini;
        double BSWF = BSWini;
        double denagF = denagini;
        double dgFO = dgini;
        double dgFG = dgini;
        double yco2FO = yco2ini;
        double yco2FG = yco2ini;
        double rgoFO = rgoini;
        double viscLini = 0 * 30 + 1 * celula[i].flui.VisOM(tL);
        double viscHini = 0 * 20 + 1 * celula[i].flui.VisOM(tH);
        double viscLF = viscLini;
        double viscHF = viscHini;
        double rholSTF = rholST;
        double rhogSTF = rhogST;
        double razdgdF = 1.;
        double razdglF = 1.;
        double titFonte = 0.;
        ProFlu fluF;
        if (celula[i].acsr.tipo == 1 && celula[i].acsr.injg.seco == 1) {
            if (celula[i].acsr.injg.QGas > 0.)
                fluF = celula[i].acsr.injg.FluidoPro;
            else
                fluF = celula[i].flui;
            dgFG = fluF.Deng;
            yco2FG = fluF.yco2;
            rhogSTF = fluF.Deng * 1.225;
        } else if (celula[i].acsr.tipo == 1 && celula[i].acsr.injg.seco == 0) {
            if (celula[i].acsr.injg.QGas > 0.)
                fluF = celula[i].acsr.injg.FluidoPro;
            else
                fluF = celula[i].flui;

            titFonte = fluF.dStockTankVaporMassFraction;
            dgFO = fluF.Deng;
            yco2FO = fluF.yco2;
            rgoFO = fluF.RGO;

            double rsF = fluF.RS(celula[i].pres, celula[i].temp) * (6.29 / 35.31467);
            razdgdF = 1 / fluF.rDgD;
            razdglF = 1 / fluF.rDgL;
            if (celula[i].acsr.injg.FluidoPro.BSW < 1 - (*vg1dSP).localtiny)
                rholSTF = (1 - fluF.BSW) * (1000 * 141.5 / (131.5 + fluF.API)) + fluF.BSW * 1000 * fluF.Denag + fluF.Deng * 1.225 * rgoFO * (1. - fluF.BSW);
            else
                rholSTF = fluF.BSW * 1000 * fluF.Denag;
            fonteO = (fonteO + fonteG) * razdgdF * (fluF.RS(celula[i].pres, celula[i].temp) * (6.29 / 35.31467) * (1. - fluF.BSW) / rholSTF);

            if (arq.nfluP > 1 || (*vg1dSP).chaverede != 0) {
                double bswaux = fluF.BSW;
                double contrabsw = 1. - bswaux;
                if (contrabsw < (*vg1dSP).localtiny)
                    contrabsw = 0.9 * (*vg1dSP).localtiny;
                double rhoPSTF;
                if (contrabsw > (*vg1dSP).localtiny)
                    rhoPSTF = (1000 * 141.5 / (131.5 + fluF.API)) + (bswaux / contrabsw) * 1000 * fluF.Denag + fluF.Deng * 1.225 * rsF;
                else
                    rhoPSTF = 1000 * fluF.Denag;
                if (contrabsw > (*vg1dSP).localtiny) {
                    fonteP *= (1. / rhoPSTF);
                    fonteA *= ((bswaux / contrabsw) / rhoPSTF);
                } else {
                    fonteP = 0.;
                    fonteA *= (1 / (fluF.BSW * 1000 * fluF.Denag));
                }
                APIF = fluF.API;
                BSWF = fluF.BSW;
                denagF = celula[i].acsr.injg.FluidoPro.Denag;
                viscLF = 0 * 30 + 1 * fluF.VisOM(tL);
                viscHF = 0 * 20 + 1 * fluF.VisOM(tH);
            }
        } else if (celula[i].acsr.tipo == 2) {
            if (celula[i].acsr.injl.QLiq > 0.)
                fluF = celula[i].acsr.injl.FluidoPro;
            else
                fluF = celula[i].flui;

            titFonte = fluF.dStockTankVaporMassFraction;
            dgFO = fluF.Deng;
            yco2FO = fluF.yco2;
            rgoFO = fluF.RGO;

            double rsF = fluF.RS(celula[i].pres, celula[i].temp) * (6.29 / 35.31467);
            razdgdF = 1 / fluF.rDgD;
            razdglF = 1 / fluF.rDgL;
            if (fluF.BSW < 1 - (*vg1dSP).localtiny)
                rholSTF = (1 - fluF.BSW) * (1000 * 141.5 / (131.5 + fluF.API)) + fluF.BSW * 1000 * fluF.Denag + fluF.Deng * 1.225 * rgoFO * (1. - fluF.BSW);
            else
                rholSTF = fluF.BSW * 1000 * fluF.Denag;
            fonteO = (fonteO + fonteG) * razdgdF * (fluF.RS(celula[i].pres, celula[i].temp) * (6.29 / 35.31467) * (1. - fluF.BSW) / rholSTF);

            if (arq.nfluP > 1 || (*vg1dSP).chaverede != 0) {
                double bswaux = fluF.BSW;
                double contrabsw = 1. - bswaux;
                if (contrabsw < (*vg1dSP).localtiny)
                    contrabsw = 0.9 * (*vg1dSP).localtiny;
                double rhoPSTF;
                if (contrabsw > (*vg1dSP).localtiny)
                    rhoPSTF = (1000 * 141.5 / (131.5 + fluF.API)) + (bswaux / contrabsw) * 1000 * fluF.Denag + fluF.Deng * 1.225 * rsF;
                else
                    rhoPSTF = 1000 * fluF.Denag;
                if (contrabsw > (*vg1dSP).localtiny) {
                    fonteP *= (1. / rhoPSTF);
                    fonteA *= ((bswaux / contrabsw) / rhoPSTF);
                } else {
                    fonteP = 0.;
                    fonteA *= (1 / (fluF.BSW * 1000 * fluF.Denag));
                }
                APIF = fluF.API;
                BSWF = fluF.BSW;
                denagF = fluF.Denag;
                viscLF = 0 * 30 + 1 * fluF.VisOM(tL);
                viscHF = 0 * 20 + 1 * fluF.VisOM(tH);
            }
        } else if (celula[i].acsr.tipo == 10) {
            if ((celula[i].acsr.injm.MassC + celula[i].acsr.injm.MassG + celula[i].acsr.injm.MassP) > 0.)
                fluF = celula[i].acsr.injm.FluidoPro;
            else
                fluF = celula[i].flui;

            titFonte = fluF.dStockTankVaporMassFraction;
            dgFO = fluF.Deng;
            yco2FO = fluF.yco2;
            rgoFO = fluF.RGO;

            double rsF = fluF.RS(celula[i].pres, celula[i].temp) * (6.29 / 35.31467);
            razdgdF = 1 / fluF.rDgD;
            razdglF = 1 / fluF.rDgL;
            rholSTF = (1 - fluF.BSW) * (1000 * 141.5 / (131.5 + fluF.API)) +
                      fluF.BSW * 1000 *
                          fluF.Denag +
                      fluF.Deng * 1.225 * rgoFO * (1. - fluF.BSW);
            fonteO = (fonteO + fonteG) * razdgdF * (fluF.RS(celula[i].pres, celula[i].temp) * (6.29 / 35.31467) * (1. - fluF.BSW) / rholSTF);

            if (arq.nfluP > 1 || (*vg1dSP).chaverede != 0) {
                double bswaux = fluF.BSW;
                double contrabsw = 1. - bswaux;
                if (contrabsw < (*vg1dSP).localtiny)
                    contrabsw = 0.9 * (*vg1dSP).localtiny;
                double rhoPSTF;
                if (contrabsw > (*vg1dSP).localtiny)
                    rhoPSTF = (1000 * 141.5 / (131.5 + fluF.API)) +
                              (bswaux / contrabsw) * 1000 *
                                  fluF.Denag +
                              fluF.Deng * 1.225 * rsF;
                else
                    rhoPSTF = 1000 * fluF.Denag;
                if (contrabsw > (*vg1dSP).localtiny) {
                    fonteP *= (1. / rhoPSTF);
                    fonteA *= ((bswaux / contrabsw) / rhoPSTF);
                } else {
                    fonteP = 0.;
                    fonteA *= (1 / (fluF.BSW * 1000 * fluF.Denag));
                }
                APIF = fluF.API;
                BSWF = fluF.BSW;
                denagF = fluF.Denag;
                viscLF = 0 * 30 + 1 * fluF.VisOM(tL);
                viscHF = 0 * 20 + 1 * fluF.VisOM(tH);
            }
        } else if (celula[i].acsr.tipo == 3) {
            if ((celula[i].acsr.ipr.Pres) > celula[i].pres)
                fluF = celula[i].acsr.ipr.FluidoPro;
            else
                fluF = celula[i].flui;

            titFonte = fluF.dStockTankVaporMassFraction;
            dgFO = fluF.Deng;
            yco2FO = fluF.yco2;
            rgoFO = fluF.RGO;

            rholSTF = (1 - fluF.BSW) * (1000 * 141.5 / (131.5 + fluF.API)) + fluF.BSW * 1000 * fluF.Denag + fluF.Deng * 1.225 * rgoFO * (1. - fluF.BSW);
            double rsF = fluF.RS(celula[i].pres, celula[i].temp) * (6.29 / 35.31467);
            razdgdF = 1 / fluF.rDgD;
            razdglF = 1 / fluF.rDgL;
            fonteO = (fonteO + fonteG) * razdgdF * (fluF.RS(celula[i].pres, celula[i].temp) * (6.29 / 35.31467) * (1. - fluF.BSW) / rholSTF);
            if (arq.nfluP > 1 || (*vg1dSP).chaverede != 0) {
                double bswaux = fluF.BSW;
                double contrabsw = 1. - bswaux;
                if (contrabsw <= (*vg1dSP).localtiny)
                    contrabsw = 0.9 * (*vg1dSP).localtiny;
                double rhoPSTF;
                if (contrabsw > 0)
                    rhoPSTF = (1000 * 141.5 / (131.5 + fluF.API)) + (bswaux / contrabsw) * 1000 * fluF.Denag + fluF.Deng * 1.225 * rsF;
                else
                    rhoPSTF = 1000 * fluF.Denag;
                if (contrabsw > (*vg1dSP).localtiny) {
                    fonteP *= (1. / rhoPSTF);
                    fonteA *= ((bswaux / contrabsw) / rhoPSTF);
                } else {
                    fonteP = 0.;
                    fonteA *= (1 / (fluF.BSW * 1000 * fluF.Denag));
                }
                APIF = fluF.API;
                BSWF = fluF.BSW;
                denagF = fluF.Denag;
                viscLF = 0 * 30 + 1 * fluF.VisOM(tL);
                viscHF = 0 * 20 + 1 * fluF.VisOM(tH);
            }
        } else if (celula[i].acsr.tipo == 9) {
            ProFlu fluF;
            if (celula[i].acsr.fontechk.presT > celula[i].acsr.fontechk.pamb) {
                fluF = celula[i].acsr.fontechk.fluidoP;
            } else {
                fluF = celula[i].acsr.fontechk.fluidoPamb;
            }
            dgFO = fluF.Deng;
            yco2FO = fluF.yco2;
            rgoFO = fluF.RGO;
            titFonte = fluF.dStockTankVaporMassFraction;

            if (fluF.BSW < 1 - (*vg1dSP).localtiny)
                rholSTF = (1 - fluF.BSW) * (1000 * 141.5 / (131.5 + fluF.API)) + fluF.BSW * 1000 * fluF.Denag + fluF.Deng * 1.225 * rgoFO * (1. - fluF.BSW);
            else
                rholSTF = fluF.BSW * 1000 * fluF.Denag;

            double rsF = fluF.RS(celula[i].pres, celula[i].temp) * (6.29 / 35.31467);
            razdgdF = 1 / fluF.rDgD;
            razdglF = 1 / fluF.rDgL;
            if (celula[i].acsr.fontechk.ambGas != 1 || (fonteO + fonteG) < 0.)
                fonteO = (fonteO + fonteG) * razdgdF * (fluF.RS(celula[i].pres, celula[i].temp) * (6.29 / 35.31467) * (1. - fluF.BSW) / rholSTF);
            else
                fonteO = 0.;
            if (arq.nfluP > 1 || (*vg1dSP).chaverede != 0) {
                double bswaux = fluF.BSW;
                double contrabsw = 1. - bswaux;
                if (contrabsw <= (*vg1dSP).localtiny)
                    contrabsw = 0.9 * (*vg1dSP).localtiny;
                double rhoPSTF;
                if (contrabsw > (*vg1dSP).localtiny)
                    rhoPSTF = (1000 * 141.5 / (131.5 + fluF.API)) + (bswaux / contrabsw) * 1000 * fluF.Denag + fluF.Deng * 1.225 * rsF;
                else
                    rhoPSTF = 1000 * fluF.Denag;
                if (contrabsw > (*vg1dSP).localtiny) {
                    fonteP *= (1. / rhoPSTF);
                    fonteA *= ((bswaux / contrabsw) / rhoPSTF);
                } else {
                    fonteP = 0.;
                    fonteA *= (1 / (fluF.BSW * 1000 * fluF.Denag));
                }
                APIF = fluF.API;
                BSWF = fluF.BSW;
                denagF = fluF.Denag;
                viscLF = 0 * 30 + 1 * fluF.VisOM(tL);
                viscHF = 0 * 20 + 1 * fluF.VisOM(tH);
            }
        } else if (celula[i].acsr.tipo == 15) {
            if ((celula[i].fontemassLR + celula[i].fontemassGR) > 1e-15)
                fluF = celula[i].acsr.radialPoro.flup;
            else
                fluF = celula[i].flui;

            titFonte = fluF.dStockTankVaporMassFraction;
            dgFO = fluF.Deng;
            yco2FO = fluF.yco2;
            rgoFO = fluF.RGO;

            rholSTF = (1 - celula[i].acsr.radialPoro.BSW) * (1000 * 141.5 / (131.5 + fluF.API)) + celula[i].acsr.radialPoro.BSW * 1000 * fluF.Denag + fluF.Deng * 1.225 * rgoFO * (1. - celula[i].acsr.radialPoro.BSW);
            double rsF = fluF.RS(celula[i].pres, celula[i].temp) * (6.29 / 35.31467);
            razdgdF = 1 / fluF.rDgD;
            razdglF = 1 / fluF.rDgL;
            fonteO = (fonteO + fonteG) * razdgdF * (fluF.RS(celula[i].pres, celula[i].temp) * (6.29 / 35.31467) * (1. - celula[i].acsr.radialPoro.BSW) / rholSTF);
            if (arq.nfluP > 1 || (*vg1dSP).chaverede != 0) {
                double bswaux = celula[i].acsr.radialPoro.BSW;
                double contrabsw = 1. - bswaux;
                if (contrabsw <= (*vg1dSP).localtiny)
                    contrabsw = 0.9 * (*vg1dSP).localtiny;
                double rhoPSTF;
                if (contrabsw > 0)
                    rhoPSTF = (1000 * 141.5 / (131.5 + fluF.API)) + (bswaux / contrabsw) * 1000 * fluF.Denag + fluF.Deng * 1.225 * rsF;
                else
                    rhoPSTF = 1000 * fluF.Denag;
                if (contrabsw > (*vg1dSP).localtiny) {
                    fonteP *= (1. / rhoPSTF);
                    fonteA *= ((bswaux / contrabsw) / rhoPSTF);
                } else {
                    fonteP = 0.;
                    fonteA *= (1 / (celula[i].acsr.radialPoro.BSW * 1000 * fluF.Denag));
                }
                APIF = fluF.API;
                BSWF = celula[i].acsr.radialPoro.BSW;
                denagF = fluF.Denag;
                viscLF = 0 * 30 + 1 * fluF.VisOM(tL);
                viscHF = 0 * 20 + 1 * fluF.VisOM(tH);
            }
        } else if (celula[i].acsr.tipo == 16) {
            if ((celula[i].fontemassLR + celula[i].fontemassGR) > 1e-15)
                fluF = celula[i].acsr.radialPoro.flup;
            else
                fluF = celula[i].flui;

            titFonte = fluF.dStockTankVaporMassFraction;
            dgFO = fluF.Deng;
            yco2FO = fluF.yco2;
            rgoFO = fluF.RGO;

            rholSTF = (1 - celula[i].acsr.poroso2D.dados.transfer.BSW) * (1000 * 141.5 / (131.5 + fluF.API)) + celula[i].acsr.poroso2D.dados.transfer.BSW * 1000 * fluF.Denag + fluF.Deng * 1.225 * rgoFO * (1. - celula[i].acsr.poroso2D.dados.transfer.BSW);
            double rsF = fluF.RS(celula[i].pres, celula[i].temp) * (6.29 / 35.31467);
            razdgdF = 1 / fluF.rDgD;
            razdglF = 1 / fluF.rDgL;
            fonteO = (fonteO + fonteG) * razdgdF * (fluF.RS(celula[i].pres, celula[i].temp) * (6.29 / 35.31467) * (1. - celula[i].acsr.poroso2D.dados.transfer.BSW) / rholSTF);
            if (arq.nfluP > 1 || (*vg1dSP).chaverede != 0) {
                double bswaux = celula[i].acsr.poroso2D.dados.transfer.BSW;
                double contrabsw = 1. - bswaux;
                if (contrabsw <= (*vg1dSP).localtiny)
                    contrabsw = 0.9 * (*vg1dSP).localtiny;
                double rhoPSTF;
                if (contrabsw > 0)
                    rhoPSTF = (1000 * 141.5 / (131.5 + fluF.API)) + (bswaux / contrabsw) * 1000 * fluF.Denag + fluF.Deng * 1.225 * rsF;
                else
                    rhoPSTF = 1000 * fluF.Denag;
                if (contrabsw > (*vg1dSP).localtiny) {
                    fonteP *= (1. / rhoPSTF);
                    fonteA *= ((bswaux / contrabsw) / rhoPSTF);
                } else {
                    fonteP = 0.;
                    fonteA *= (1 / (celula[i].acsr.poroso2D.dados.transfer.BSW * 1000 * fluF.Denag));
                }
                APIF = fluF.API;
                BSWF = celula[i].acsr.poroso2D.dados.transfer.BSW;
                denagF = fluF.Denag;
                viscLF = 0 * 30 + 1 * fluF.VisOM(tL);
                viscHF = 0 * 20 + 1 * fluF.VisOM(tH);
            }
        } else if ((fabs(fonteO) > (*vg1dSP).localtiny && celula[i].acsr.tipo != 2 && celula[i].acsr.tipo != 3 &&
                    celula[i].acsr.tipo != 9 && celula[i].acsr.tipo != 15 && celula[i].acsr.tipo != 16) ||
                   (fabs(fonteG) > (*vg1dSP).localtiny && celula[i].acsr.tipo != 1 && celula[i].acsr.tipo != 2 && celula[i].acsr.tipo != 3 && celula[i].acsr.tipo != 9 && celula[i].acsr.tipo != 15 && celula[i].acsr.tipo != 16)) {
            if (celula[i].acsr.tipo == 5 || celula[i].acsr.tipo == 8) {
                dgFO = celula[i].flui.Deng;
                yco2FO = celula[i].flui.yco2;
                rgoFO = celula[i].flui.RGO;
                titFonte = celula[i].flui.dStockTankVaporMassFraction;
                double rholiq = (1 - celula[i].flui.BSW) * (1000 * 141.5 / (131.5 + celula[i].flui.API)) + celula[i].flui.BSW * 1000 * celula[i].flui.Denag;
                double rhogas = celula[i].flui.Deng * 1.225;
                rholSTF = rholiq + rhogas * rgoFO * (1. - celula[i].flui.BSW);
                rhogSTF = rhogas;
                double rsF = celula[i].flui.RS(celula[i].pres, celula[i].temp) * (6.29 / 35.31467);
                razdgdF = 1 / celula[i].flui.rDgD;
                razdglF = 1 / celula[i].flui.rDgL;
                fonteO = (fonteO + fonteG) * (razdgdF * rsF * (1. - celula[i].flui.BSW) / rholSTF);
                if (arq.nfluP > 1 || (*vg1dSP).chaverede != 0) {
                    double rhoPSTF = (1 - celula[i].flui.BSW) * (1000 * 141.5 / (131.5 + celula[i].flui.API)) + celula[i].flui.BSW * 1000 * celula[i].flui.Denag;
                    fonteP *= ((1 - celula[i].flui.BSW) / rhoPSTF);
                    fonteA *= (celula[i].flui.BSW / rhoPSTF);
                    APIF = celula[i].flui.API;
                    BSWF = celula[i].flui.BSW;
                    viscLF = 0 * 30 + 1 * celula[i].flui.VisOM(tL);
                    viscHF = 0 * 20 + 1 * celula[i].flui.VisOM(tH);
                }
            } else if ((*celula[i].acsrL).tipo == 5 || (*celula[i].acsrL).tipo == 8) {
                double rholiq;
                double rhogas;

                if (i > 0) {
                    dgFO = celula[i - 1].flui.Deng;
                    yco2FO = celula[i - 1].flui.yco2;
                    rgoFO = celula[i - 1].flui.RGO;
                    rholiq = (1 - celula[i - 1].flui.BSW) * (1000 * 141.5 / (131.5 + celula[i - 1].flui.API)) +
                             celula[i - 1].flui.BSW * 1000 * celula[i - 1].flui.Denag;
                    rhogas = celula[i - 1].flui.Deng * 1.225;
                    titFonte = celula[i - 1].flui.dStockTankVaporMassFraction;
                    rholSTF = rholiq + rhogas * rgoFO * (1. - celula[i - 1].flui.BSW);
                    rhogSTF = rhogas;
                    double rsF = celula[i - 1].flui.RS(celula[i].pres, celula[i].temp) * (6.29 / 35.31467);
                    razdgdF = 1 / celula[i - 1].flui.rDgD;
                    razdglF = 1 / celula[i - 1].flui.rDgL;
                    fonteO = (fonteO + fonteG) * (razdgdF * rsF *
                                                  (1. - celula[i - 1].flui.BSW) / rholSTF);
                    if (arq.nfluP > 1 || (*vg1dSP).chaverede != 0) {
                        double rhoPSTF = (1 - celula[i - 1].flui.BSW) * (1000 * 141.5 / (131.5 + celula[i - 1].flui.API)) +
                                         celula[i - 1].flui.BSW * 1000 *
                                             celula[i - 1].flui.Denag;
                        fonteP *= ((1 - celula[i - 1].flui.BSW) / rhoPSTF);
                        fonteA *= (celula[i - 1].flui.BSW / rhoPSTF);
                        APIF = celula[i - 1].flui.API;
                        BSWF = celula[i - 1].flui.BSW;
                        denagF = celula[i - 1].flui.Denag;
                        viscLF = 0 * 30 + 1 * celula[i - 1].flui.VisOM(tL);
                        viscHF = 0 * 20 + 1 * celula[i - 1].flui.VisOM(tH);
                    }
                } else {
                    dgFO = celula[i].flui.Deng;
                    yco2FO = celula[i].flui.yco2;
                    rgoFO = celula[i].flui.RGO;
                    rholiq = (1 - celula[i].flui.BSW) * (1000 * 141.5 / (131.5 + celula[i].flui.API)) +
                             celula[i].flui.BSW * 1000 * celula[i].flui.Denag;
                    rhogas = celula[i].flui.Deng * 1.225;
                    titFonte = celula[i].flui.dStockTankVaporMassFraction;
                    rholSTF = rholiq + rhogas * rgoFO * (1. - celula[i].flui.BSW);
                    rhogSTF = rhogas;
                    double rsF = celula[i].flui.RS(celula[i].pres, celula[i].temp) * (6.29 / 35.31467);
                    razdgdF = 1 / celula[i].flui.rDgD;
                    razdglF = 1 / celula[i].flui.rDgL;
                    fonteO = (fonteO + fonteG) * (razdgdF * celula[i].flui.RS(celula[i].pres, celula[i].temp) * (6.29 / 35.31467) *
                                                  (1. - celula[i].flui.BSW) / rholSTF);
                    if (arq.nfluP > 1 || (*vg1dSP).chaverede != 0) {
                        double rhoPSTF = (1 - celula[i].flui.BSW) * (1000 * 141.5 / (131.5 + celula[i].flui.API)) +
                                         celula[i].flui.BSW * 1000 *
                                             celula[i].flui.Denag;
                        fonteP *= ((1 - celula[i].flui.BSW) / rhoPSTF);
                        fonteA *= (celula[i].flui.BSW / rhoPSTF);
                        APIF = celula[i].flui.API;
                        BSWF = celula[i].flui.BSW;
                        denagF = celula[i].flui.Denag;
                        viscLF = 0 * 30 + 1 * celula[i].flui.VisOM(tL);
                        viscHF = 0 * 20 + 1 * celula[i].flui.VisOM(tH);
                    }
                }
            }
        }
        fonteG *= (razdglF / (rhogSTF));

        if (titFonte > 1. - 1e-15) {
            fonteO = 0.;
            fonteP = 0.;
            fonteA = 0.;
        }

        MultOe = 0.;
        if (celula[i].flui.dStockTankVaporMassFraction < 1. - 1e-15)
            MultOe = celula[i].QL * (1 - betI0) * (1 - bsw0) * razdgd0 * rs0 / bo0;
        MultOd = 0.;
        if (celula[i + 1].flui.dStockTankVaporMassFraction < 1. - 1e-15)
            MultOd = celula[i + 1].QL * (1 - betI1) * (1 - bsw1) * razdgd1 * rs1 / bo1;
        double MultGe = (celula[i].MC - celula[i].Mliqini) * razdgl0 / (rhogST0);
        double MultGd = (celula[i + 1].MC - celula[i + 1].Mliqini) * razdgl1 / (rhogST1);
        double volleveFim = (((1 - hol) * rhog * razdgl / (rhogST)) + hol * (1 - bet) * (1. - bsw) * rs * razdgd / (bo));
        if (volleveFim < 1e-15)
            volleveFim = 0.;
        double residuo = (volleveFim - celula[i].VolLeveST) * at / dt + (MultOd - MultOe) / dx + (MultGd - MultGe) / dx - (fonteO / dx + fonteG / dx);
        double volpesFim = hol * (1 - bet) * (1 - bsw) / bo;
        double volaguaFim = hol * (1 - bet) * bsw; // nao deveria ser dividido por Bo???????????
        double MultPe = 0.;
        double MultPd = 0.;
        double residuoP = 0.;
        double MultAe = 0.;
        double MultAd = 0.;
        double residuoA;
        if (arq.nfluP > 1 || (*vg1dSP).chaverede != 0) {
            MultPe = 0.;
            if (celula[i].flui.dStockTankVaporMassFraction < 1. - 1e-15)
                MultPe = celula[i].QL * (1 - betI0) * (1 - bsw0) / bo0;
            MultPd = 0.;
            if (celula[i + 1].flui.dStockTankVaporMassFraction < 1. - 1e-15)
                MultPd = celula[i + 1].QL * (1 - betI1) * (1 - bsw1) / bo1;
            residuoP = (volpesFim - celula[i].VolPesaST) * at / dt + (MultPd - MultPe) / dx - fonteP / dx;
            MultAe = 0.;
            if (celula[i].flui.dStockTankVaporMassFraction < 1. - 1e-15)
                MultAe = celula[i].QL * (1 - betI0) * bsw0 / bo0;
            MultAd = 0.;
            if (celula[i + 1].flui.dStockTankVaporMassFraction < 1. - 1e-15)
                MultAd = celula[i + 1].QL * (1 - betI1) * bsw1 / bo1;
            residuoA = (volaguaFim - celula[i].VolAguaST) * at / dt + (MultAd - MultAe) / dx - fonteA / dx;
        }
        rgo[i] = (*vg1dSP).RGOMax;
        if (hol > (*vg1dSP).localtiny && bet < (1. - (*vg1dSP).localtiny) && bsw < (1. - (*vg1dSP).localtiny)) {
            rgo[i] = (volleveFim - residuo * dt / at) * bo / (hol * (1 - bet) * (1 - bsw));
            if (rgo[i] > (*vg1dSP).RGOMax)
                rgo[i] = (*vg1dSP).RGOMax;
        } else if (bet >= (1. - (*vg1dSP).localtiny) || bsw >= (1. - (*vg1dSP).localtiny))
            rgo[i] = 0.;
        else
            rgo[i] = (*vg1dSP).RGOMax;

        if (volleveFim > 1e-5 && arq.flashCompleto == 0 && ((fonteG >= 0 || fonteO > 0) || ((MultGd < 0 || MultGe > 0) || (MultOd < 0 || MultOe > 0)))) {
            dg[i] = (dt * (dgFO * fonteO / dx + dgFG * fonteG / dx + 1. * dgini * residuo - (dg1O * MultOd - dg0O * MultOe) / dx - (dg1G * MultGd - dg0G * MultGe) / dx) + dgini * celula[i].VolLeveST * at) /
                    (volleveFim * at - 0. * residuo * dt);
            yco2[i] = (dt * (yco2FO * fonteO / dx + yco2FG * fonteG / dx + 1. * yco2ini * residuo - (yco21O * MultOd - yco20O * MultOe) / dx - (yco21G * MultGd - yco20G * MultGe) / dx) + yco2ini * celula[i].VolLeveST * at) /
                      (volleveFim * at - 0. * residuo * dt);
            if (yco2[i] < 0.)
                yco2[i] = 0.;
            else if (yco2[i] > 1.)
                yco2[i] = 1.;
        } else {
            dg[i] = dgini;
            yco2[i] = yco2ini;
        }
        if ((arq.nfluP > 1 && arq.flashCompleto == 0) || (*vg1dSP).chaverede != 0) {
            if (volpesFim > 1e-3 && (fonteP > 0 || (MultPd < 0 || MultPe > 0))) {
                double denmixSTDF = 141.5 / (131.5 + APIF);
                double denmixSTDini = 141.5 / (131.5 + APIini);
                double denmixSTD1 = 141.5 / (131.5 + API1);
                double denmixSTD0 = 141.5 / (131.5 + API0);
                API[i] = (dt * (denmixSTDF * fonteP / dx + 1. * denmixSTDini * residuoP - (denmixSTD1 * MultPd - denmixSTD0 * MultPe) / dx) + denmixSTDini * celula[i].VolPesaST * at) / (volpesFim * at - 0. * residuoP * dt);
                API[i] = 141.5 / API[i] - 131.5;
                VISCL[i] = (dt * (viscLF * fonteP / dx + 1. * viscLini * residuoP - (viscL1 * MultPd - viscL0 * MultPe) / dx) + viscLini * celula[i].VolPesaST * at) / (volpesFim * at - 0. * residuoP * dt);
                VISCH[i] = (dt * (viscHF * fonteP / dx + 1. * viscHini * residuoP - (viscH1 * MultPd - viscH0 * MultPe) / dx) + viscHini * celula[i].VolPesaST * at) / (volpesFim * at - 0. * residuoP * dt);
            } else {
                API[i] = APIini;
                VISCL[i] = viscLini;
                VISCH[i] = viscHini;
            }
        }
        if (arq.nfluP > 1 || (*vg1dSP).chaverede != 0) {
            if ((volaguaFim + volpesFim) > 1e-3 && ((fonteP > 0 || fonteA > 0) ||
                                                    ((MultPd < 0 || MultPe > 0) || (MultAd < 0 || MultAe > 0)))) {
                BSW[i] = (dt * (BSWF * (fonteA + fonteP) / dx + 1. * BSWini * (residuoA + residuoP) - (BSW1 * (MultAd + MultPd) - BSW0 * (MultAe + MultPe)) / dx) + BSWini * (celula[i].VolAguaST + celula[i].VolPesaST) * at) /
                         ((volaguaFim + volpesFim) * at - 0. * (residuoA + residuoP) * dt);
                denag[i] = (dt * (denagF * (fonteA) / dx + 1. * denagini * (residuoA) - (denag1 * (MultAd)-denag0 * (MultAe)) / dx) + denagini * (celula[i].VolAguaST) * at) /
                           ((volaguaFim)*at - 0. * (residuoA)*dt);
                if (BSW[i] < 0.)
                    BSW[i] = 0.;
                else if (BSW[i] > 1.)
                    BSW[i] = 1.;
                if (denag[i] < 1.)
                    denag[i] = 1.;
            } else
                BSW[i] = BSWini;
            denag[i] = denagini;
        }
        celula[i].VolLeveST = volleveFim;
        celula[i].VolPesaST = volpesFim;
        celula[i].VolAguaST = volaguaFim;
    }
    for (int i = imin; i <= ncel - 1; i++) {
        if (trackRGO > 0 && celula[i].flui.corrSat != 4)
            celula[i].flui.RGO = rgo[i];
        if (trackDeng > 0 && arq.flashCompleto == 0) {
            celula[i].flui.Deng = dg[i];
            if ((celula[i].flui.Deng > 5 || celula[i].flui.Deng < 0) && i > 0)
                celula[i].flui.Deng = celula[i - 1].flui.Deng;
            celula[i].flui.yco2 = yco2[i];
        }
        if (arq.nfluP > 1 || (*vg1dSP).chaverede != 0) {
            celula[i].flui.BSW = BSW[i];
            celula[i].flui.Denag = denag[i];
            if (arq.flashCompleto == 0) {
                celula[i].flui.API = API[i];
                celula[i].flui.LVisL = VISCL[i];
                celula[i].flui.LVisH = VISCH[i];
                celula[i].flui.TempL = tL;
                celula[i].flui.TempH = tH;
            }
        }
        if (arq.flashCompleto == 0)
            celula[i].flui.RenovaFluido();
        corrDeng(i);
    }

    if ((*vg1dSP).chaverede == 0 || noextremo == 1 || celula[ncel].Mliqini > -(*vg1dSP).localtiny) {
        if (trackRGO > 0 && celula[ncel].flui.corrSat != 4)
            celula[ncel].flui.RGO = rgo[ncel - 1];
        if (trackDeng > 0 && arq.flashCompleto == 0) {
            celula[ncel].flui.Deng = dg[ncel - 1];
            celula[ncel].flui.yco2 = yco2[ncel - 1];
        }
        if (arq.nfluP > 1 || (*vg1dSP).chaverede != 0) {
            celula[ncel].flui.BSW = BSW[ncel - 1];
            celula[ncel].flui.Denag = denag[ncel - 1];
            if (arq.flashCompleto == 0) {
                celula[ncel].flui.API = API[ncel - 1];
                celula[ncel].flui.LVisL = VISCL[ncel - 1];
                celula[ncel].flui.LVisH = VISCH[ncel - 1];
                celula[ncel].flui.TempL = tL;
                celula[ncel].flui.TempH = tH;
            }
        }
    } else {
        celula[ncel].flui.RGO = fluiRev.RGO;
        celula[ncel].flui.Deng = fluiRev.Deng;
        celula[ncel].flui.yco2 = fluiRev.yco2;
        celula[ncel].flui.API = fluiRev.API;
        celula[ncel].flui.BSW = fluiRev.BSW;
        celula[ncel].flui.Denag = fluiRev.Denag;
        celula[ncel].flui.LVisL = fluiRev.LVisL;
        celula[ncel].flui.LVisH = fluiRev.LVisH;
        celula[ncel].flui.TempL = fluiRev.TempL;
        celula[ncel].flui.TempH = fluiRev.TempH;
    }
    if (arq.flashCompleto == 0)
        celula[ncel].flui.RenovaFluido();
    if (arq.corrDeng == 0) {
        celula[ncel].flui.rDgD = 1.;
        celula[ncel].flui.rDgL = 1.;
        celula[ncel].flui.PCis = celula[ncel].flui.PC;
        celula[ncel].flui.TCis = celula[ncel].flui.TC;
    } else {
        celula[ncel].flui.razDegD(celula[ncel].pres, celula[ncel].temp);
        celula[ncel].flui.rzDegL(celula[ncel].pres, celula[ncel].temp);
        celula[ncel].flui.PcTcIS();
    }
}

/*** alteracao4 ***/

void SProd::renovaFracMol(ProFlu fluiRev) {
    Vcr<double> BSW(ncel);
    Vcr<double> VISCL(ncel);
    Vcr<double> VISCH(ncel);
    int ncomp = arq.npseudo;
    ProFlu *fluC;
    fluC = new ProFlu[ncel];
    Vcr<double> fracMol0(ncomp);
    Vcr<double> fracMol1(ncomp);
    Vcr<double> fracMolF(ncomp);
    double pesoMol0;
    double pesoMol1;
    double pesoMolF;
    double pesoMolC;
    double dt = celula[1].dt;
    double tL = 0.;
    double tH = 70.;
    if (celula[0].acsr.tipo == 15 && (celula[0].fontemassCR + celula[0].fontemassGR + celula[0].fontemassLR) > 0.)
        celula[0].flui.BSW = celula[0].acsr.radialPoro.BSW;
    else if (celula[0].acsr.tipo == 16 && (celula[0].fontemassCR + celula[0].fontemassGR + celula[0].fontemassLR) > 0.)
        celula[0].flui.BSW = celula[0].acsr.poroso2D.dados.transfer.BSW;
    int imin = 1;
    if ((*vg1dSP).chaverede != 0 && (celula[0].acsr.tipo == 10 || arq.ConContEntrada > 0))
        imin = 0;
    for (int i = imin; i < ncel; i++) {
        double at = celula[i].duto.area;
        double dx = celula[i].dx;
        double ti0;
        if (celula[i].VTemper < 0.)
            ti0 = celula[i].temp;
        else {
            if (i > 0)
                ti0 = celula[i - 1].temp;
            else if (arq.ConContEntrada == 1)
                ti0 = tempE;
            else
                ti0 = celula[i].temp;
        }
        double ti1 = celula[i].temp;
        if (celula[i + 1].VTemper < 0.)
            ti1 = celula[i + 1].temp;

        double titF = 0.;
        ProFlu fluF;

        double boF = 1.;
        double baF = 1.;
        double fwF = 1.;
        double rhoOF = 900.;
        double rhoWF = 1000.;

        double betI0;
        double bo0;
        double ba0;
        double bsw0;
        double BSW0;
        double viscL0;
        double viscH0;
        if (i > 0 || arq.ConContEntrada == 0) {
            if (i > 0 && celula[i].QG >= 0.)
                betI0 = celula[i - 1].betPigD;
            else
                betI0 = celula[i].betPigE;
        } else {
            if (celula[i].QG >= 0.)
                betI0 = betaE;
            else
                betI0 = celula[i].betPigE;
        }
        if ((i > 0 || arq.ConContEntrada == 1) && celula[i].QL >= 0.) {
            double pe;
            double te;
            if (i > 0) {
                pe = celula[i - 1].pres;
                te = celula[i - 1].temp;
            } else {
                pe = presE;
                te = tempE;
            }
            if (arq.ConContEntrada == 0)
                betI0 = celula[i - 1].betPigD; // testeBeta
            else
                betI0 = betaE; // testeBeta
            double rs0 = (*celula[i].fluiL).RS(pe, te);
            bo0 = (*celula[i].fluiL).BOFunc(pe, te, rs0);
            ba0 = (*celula[i].fluiL).BAFunc(pe, te);
            bsw0 = (*celula[i].fluiL).BSW * ba0 / (bo0 + ba0 * (*celula[i].fluiL).BSW - (*celula[i].fluiL).BSW * bo0);
            BSW0 = (*celula[i].fluiL).BSW;
            viscL0 = 0 * 30 + 1 * (*celula[i].fluiL).VisOM(tL);
            viscH0 = 0 * 20 + 1 * (*celula[i].fluiL).VisOM(tH);
        } else {
            betI0 = celula[i].betPigE; // testebeta
            double rs0 = celula[i].flui.RS(celula[i].pres, celula[i].temp);
            bo0 = celula[i].flui.BOFunc(celula[i].pres, celula[i].temp, rs0);
            ba0 = celula[i].flui.BAFunc(celula[i].pres, celula[i].temp);
            bsw0 = celula[i].flui.BSW * ba0 / (bo0 + ba0 * celula[i].flui.BSW - celula[i].flui.BSW * bo0);
            BSW0 = celula[i].flui.BSW;
            viscL0 = 0 * 30 + 1 * celula[i].flui.VisOM(tL);
            viscH0 = 0 * 20 + 1 * celula[i].flui.VisOM(tH);
        }
        if (bo0 < 1e-15)
            bo0 = 1e-15;

        double betI1 = celula[i].betPigD;
        double rs1 = celula[i].flui.RS(celula[i].pres, celula[i].temp);
        double bo1 = celula[i].flui.BOFunc(celula[i].pres, celula[i].temp, rs1);
        double ba1 = celula[i].flui.BAFunc(celula[i].pres, celula[i].temp);
        double bsw1 = celula[i].flui.BSW * ba1 / (bo1 + ba1 * celula[i].flui.BSW - celula[i].flui.BSW * bo1);
        double BSW1 = celula[i].flui.BSW;
        double viscL1 = 0 * 30 + 1 * celula[i].flui.VisOM(tL);
        double viscH1 = 0 * 20 + 1 * celula[i].flui.VisOM(tH);

        // betI1 = celula[i + 1].betPigE;    //duvidabeta
        if (celula[i + 1].QL < 0.) {
            betI1 = celula[i + 1].betPigE; // testeBet
            rs1 = celula[i + 1].flui.RS(celula[i + 1].pres, celula[i + 1].temp);
            bo1 = celula[i + 1].flui.BOFunc(celula[i + 1].pres, celula[i + 1].temp, rs1);
            ba1 = celula[i + 1].flui.BAFunc(celula[i + 1].pres, celula[i + 1].temp);
            bsw1 = celula[i + 1].flui.BSW * ba1 / (bo1 + ba1 * celula[i + 1].flui.BSW - celula[i + 1].flui.BSW * bo1);
            BSW1 = celula[i + 1].flui.BSW;
            viscL1 = 0 * 30 + 1 * celula[i + 1].flui.VisOM(tL);
            viscH1 = 0 * 20 + 1 * celula[i + 1].flui.VisOM(tH);
        }
        if (bo1 < 1e-15)
            bo1 = 1e-15;
        double titV0;
        double titV1;
        double vazMasLiq0 = celula[i].MliqiniL;
        double vazMasLiq1 = celula[i + 1].MliqiniL;
        double vazMasGas0 = celula[i].MC - celula[i].MliqiniL;
        double vazMasGas1 = celula[i + 1].MC - celula[i + 1].MliqiniL;
        pesoMolC = 0.;
        for (int kfrac = 0; kfrac < ncomp; kfrac++) {
            pesoMolC += celula[i].flui.masMol[kfrac] * celula[i].flui.fracMol[kfrac];
        }
        if ((i > 0 || arq.ConContEntrada == 1) && celula[i].MC >= 0.) {
            double pe;
            double te;
            if (i > 0) {
                pe = celula[i - 1].pres;
                te = celula[i - 1].temp;
            } else {
                pe = presE;
                te = tempE;
            }
            double fwV = bsw0;
            double rhoOV = (*celula[i].fluiL).MasEspoleo(pe, te);
            double rhoWV = (*celula[i].fluiL).MasEspAgua(pe, te);
            titV0 = (1 - fwV) * rhoOV / ((1 - fwV) * rhoOV + fwV * rhoWV);
            vazMasLiq0 -= betI0 * celula[i].QL * celula[i].fluicol.MasEspFlu(pe, te);
            vazMasLiq0 *= titV0;
            pesoMol0 = 0;
            for (int kfrac = 0; kfrac < ncomp; kfrac++) {
                fracMol0[kfrac] = (*celula[i].fluiL).fracMol[kfrac];
                pesoMol0 += (*celula[i].fluiL).masMol[kfrac] * (*celula[i].fluiL).fracMol[kfrac];
            }
        } else {
            double fwV = bsw0;
            double rhoOV = celula[i].flui.MasEspoleo(celula[i].pres, celula[i].temp);
            double rhoWV = celula[i].flui.MasEspAgua(celula[i].pres, celula[i].temp);
            titV0 = (1 - fwV) * rhoOV / ((1 - fwV) * rhoOV + fwV * rhoWV);
            vazMasLiq0 -= betI0 * celula[i].QL * celula[i].fluicol.MasEspFlu(celula[i].pres, celula[i].temp);
            vazMasLiq0 *= titV0;
            pesoMol0 = 0;
            for (int kfrac = 0; kfrac < ncomp; kfrac++) {
                fracMol0[kfrac] = celula[i].flui.fracMol[kfrac];
                pesoMol0 += celula[i].flui.masMol[kfrac] * celula[i].flui.fracMol[kfrac];
            }
        }
        if (celula[i + 1].MC >= 0.) {
            double fwV = bsw1;
            double rhoOV = celula[i].flui.MasEspoleo(celula[i].pres, celula[i].temp);
            double rhoWV = celula[i].flui.MasEspAgua(celula[i].pres, celula[i].temp);
            titV1 = (1 - fwV) * rhoOV / ((1 - fwV) * rhoOV + fwV * rhoWV);
            vazMasLiq1 -= betI1 * celula[i].QL * celula[i].fluicol.MasEspFlu(celula[i].pres, celula[i].temp);
            vazMasLiq1 *= titV1;
            pesoMol1 = 0;
            for (int kfrac = 0; kfrac < ncomp; kfrac++) {
                fracMol1[kfrac] = celula[i].flui.fracMol[kfrac];
                pesoMol1 += celula[i].flui.masMol[kfrac] * celula[i].flui.fracMol[kfrac];
            }
        } else {
            double fwV = bsw1;
            double rhoOV = celula[i + 1].flui.MasEspoleo(celula[i + 1].pres, celula[i + 1].temp);
            double rhoWV = celula[i + 1].flui.MasEspAgua(celula[i + 1].pres, celula[i + 1].temp);
            titV1 = (1 - fwV) * rhoOV / ((1 - fwV) * rhoOV + fwV * rhoWV);
            vazMasLiq1 -= betI1 * celula[i + 1].QL * celula[i].fluicol.MasEspFlu(celula[i + 1].pres, celula[i + 1].temp);
            vazMasLiq1 *= titV1;
            pesoMol1 = 0;
            for (int kfrac = 0; kfrac < ncomp; kfrac++) {
                fracMol1[kfrac] = celula[i + 1].flui.fracMol[kfrac];
                pesoMol1 += celula[i + 1].flui.masMol[kfrac] * celula[i + 1].flui.fracMol[kfrac];
            }
        }
        double fonteO = celula[i].fontemassLR;
        double fonteG = celula[i].fontemassGR;
        double fonteP = celula[i].fontemassLR;
        double fonteA = celula[i].fontemassLR;

        if (celula[i].acsr.tipo == 1) {
            celula[i].acsr.injg.FluidoPro.atualizaPropComp(celula[i].pres, celula[i].temp, celula[i].flui.dCalculatedBeta,
                                                           celula[i].flui.oCalculatedLiqComposition, celula[i].flui.oCalculatedVapComposition, arq.pocinjec);
            if (celula[i].acsr.injg.QGas > 0.)
                fluF = celula[i].acsr.injg.FluidoPro;
            else
                fluF = celula[i].flui;
            fwF = 0.;
            titF = 1.;
        } else if (celula[i].acsr.tipo == 2) {
            celula[i].acsr.injl.FluidoPro.atualizaPropComp(celula[i].pres, celula[i].temp, celula[i].flui.dCalculatedBeta,
                                                           celula[i].flui.oCalculatedLiqComposition, celula[i].flui.oCalculatedVapComposition, arq.pocinjec);
            if (celula[i].acsr.injl.QLiq > 0.)
                fluF = celula[i].acsr.injl.FluidoPro;
            else
                fluF = celula[i].flui;
            boF = fluF.BOFunc(celula[i].pres, celula[i].temp);
            baF = fluF.BAFunc(celula[i].pres, celula[i].temp);
            fwF = fluF.BSW * baF / (boF + baF * fluF.BSW - fluF.BSW * boF);
            rhoOF = fluF.MasEspoleo(celula[i].pres, celula[i].temp);
            rhoWF = fluF.MasEspAgua(celula[i].pres, celula[i].temp);
            titF = (1 - fwF) * rhoOF / ((1 - fwF) * rhoOF + fwF * rhoWF);
        } else if (celula[i].acsr.tipo == 3) {
            celula[i].acsr.ipr.FluidoPro.atualizaPropComp(celula[i].pres, celula[i].temp, celula[i].flui.dCalculatedBeta,
                                                          celula[i].flui.oCalculatedLiqComposition, celula[i].flui.oCalculatedVapComposition, arq.pocinjec);
            if ((celula[i].acsr.ipr.Pres) > celula[i].pres)
                fluF = celula[i].acsr.ipr.FluidoPro;
            else
                fluF = celula[i].flui;
            boF = fluF.BOFunc(celula[i].pres, celula[i].temp);
            baF = fluF.BAFunc(celula[i].pres, celula[i].temp);
            fwF = fluF.BSW * baF / (boF + baF * fluF.BSW - fluF.BSW * boF);
            rhoOF = fluF.MasEspoleo(celula[i].pres, celula[i].temp);
            rhoWF = fluF.MasEspAgua(celula[i].pres, celula[i].temp);
            titF = (1 - fwF) * rhoOF / ((1 - fwF) * rhoOF + fwF * rhoWF);
        } else if (celula[i].acsr.tipo == 10) {
            celula[i].acsr.injm.FluidoPro.atualizaPropComp(celula[i].pres, celula[i].temp, celula[i].flui.dCalculatedBeta,
                                                           celula[i].flui.oCalculatedLiqComposition, celula[i].flui.oCalculatedVapComposition, arq.pocinjec);
            if ((celula[i].acsr.injm.MassC + celula[i].acsr.injm.MassG + celula[i].acsr.injm.MassP) > 0.)
                fluF = celula[i].acsr.injm.FluidoPro;
            else
                fluF = celula[i].flui;
            boF = fluF.BOFunc(celula[i].pres, celula[i].temp);
            baF = fluF.BAFunc(celula[i].pres, celula[i].temp);
            fwF = fluF.BSW * baF / (boF + baF * fluF.BSW - fluF.BSW * boF);
            rhoOF = fluF.MasEspoleo(celula[i].pres, celula[i].temp);
            rhoWF = fluF.MasEspAgua(celula[i].pres, celula[i].temp);
            titF = (1 - fwF) * rhoOF / ((1 - fwF) * rhoOF + fwF * rhoWF);
        } else if (celula[i].acsr.tipo == 9 && celula[i].acsr.fontechk.abertura > 1e-6) {
            celula[i].acsr.fontechk.fluidoP.atualizaPropComp(celula[i].pres, celula[i].temp, celula[i].flui.dCalculatedBeta,
                                                             celula[i].flui.oCalculatedLiqComposition, celula[i].flui.oCalculatedVapComposition, arq.pocinjec);
            celula[i].acsr.fontechk.fluidoPamb.atualizaPropComp(celula[i].pres, celula[i].temp, celula[i].flui.dCalculatedBeta,
                                                                celula[i].flui.oCalculatedLiqComposition, celula[i].flui.oCalculatedVapComposition, arq.pocinjec);
            if (celula[i].acsr.fontechk.presT > celula[i].acsr.fontechk.pamb) {
                fluF = celula[i].acsr.fontechk.fluidoP;
            } else {
                fluF = celula[i].acsr.fontechk.fluidoPamb;
            }
            boF = fluF.BOFunc(celula[i].pres, celula[i].temp);
            baF = fluF.BAFunc(celula[i].pres, celula[i].temp);
            fwF = fluF.BSW * baF / (boF + baF * fluF.BSW - fluF.BSW * boF);
            rhoOF = fluF.MasEspoleo(celula[i].pres, celula[i].temp);
            rhoWF = fluF.MasEspAgua(celula[i].pres, celula[i].temp);
            titF = (1 - fwF) * rhoOF / ((1 - fwF) * rhoOF + fwF * rhoWF);
        } else if (celula[i].acsr.tipo == 15) {
            double tRes = celula[i].acsr.radialPoro.tRes;
            celula[i].acsr.radialPoro.flup.atualizaPropComp(celula[i].pres, tRes, celula[i].flui.dCalculatedBeta,
                                                            celula[i].flui.oCalculatedLiqComposition, celula[i].flui.oCalculatedVapComposition, arq.pocinjec);
            if ((celula[i].fontemassLR + celula[i].fontemassGR) > 1e-15)
                fluF = celula[i].acsr.radialPoro.flup;
            else
                fluF = celula[i].flui;
            boF = fluF.BOFunc(celula[i].pres, celula[i].temp);
            baF = fluF.BAFunc(celula[i].pres, celula[i].temp);
            fwF = fluF.BSW * baF / (boF + baF * fluF.BSW - fluF.BSW * boF);
            rhoOF = fluF.MasEspoleo(celula[i].pres, celula[i].temp);
            rhoWF = fluF.MasEspAgua(celula[i].pres, celula[i].temp);
            titF = (1 - fwF) * rhoOF / ((1 - fwF) * rhoOF + fwF * rhoWF);
        } else if (celula[i].acsr.tipo == 16) {
            double tRes = celula[i].acsr.poroso2D.dados.tRes;
            celula[i].acsr.poroso2D.dados.flup.atualizaPropComp(celula[i].pres, tRes, celula[i].flui.dCalculatedBeta,
                                                                celula[i].flui.oCalculatedLiqComposition, celula[i].flui.oCalculatedVapComposition, arq.pocinjec);
            if ((celula[i].fontemassLR + celula[i].fontemassGR) > 1e-15)
                fluF = celula[i].acsr.poroso2D.dados.flup;
            else
                fluF = celula[i].flui;
            boF = fluF.BOFunc(celula[i].pres, celula[i].temp);
            baF = fluF.BAFunc(celula[i].pres, celula[i].temp);
            fwF = fluF.BSW * baF / (boF + baF * fluF.BSW - fluF.BSW * boF);
            rhoOF = fluF.MasEspoleo(celula[i].pres, celula[i].temp);
            rhoWF = fluF.MasEspAgua(celula[i].pres, celula[i].temp);
            titF = (1 - fwF) * rhoOF / ((1 - fwF) * rhoOF + fwF * rhoWF);
        }

        pesoMolF = 0;
        for (int j = 0; j < ncomp; j++)
            pesoMolF += fluF.masMol[j] * fluF.fracMol[j];
        fonteO *= titF;
        double tempMol;
        fluC[i] = celula[i].flui;
        for (int corrige = 0; corrige < 2; corrige++) {
            tempMol = (fluC[i].MasEspoleo(celula[i].pres, celula[i].temp) * (1 - celula[i].alf) +
                       fluC[i].MasEspGas(celula[i].pres, celula[i].temp) *
                           celula[i].alf) *
                      celula[i].duto.area / pesoMolC;
            for (int kfrac = 0; kfrac < ncomp; kfrac++) {
                fluC[i].fracMol[kfrac] = (celula[i].nMolIni * celula[i].flui.fracMol[kfrac] +
                                          ((fonteO + fonteG) * fracMolF[kfrac] / pesoMolF -
                                           (vazMasLiq1 + vazMasGas1) * fracMol1[kfrac] / pesoMol1 -
                                           (vazMasLiq0 + vazMasGas0) * fracMol0[kfrac] / pesoMol0) *
                                              dt / dx) /
                                         tempMol;
            }
            pesoMolC = 0;
            for (int kfrac = 0; kfrac < ncomp; kfrac++)
                pesoMolC += fluC[i].masMol[kfrac] * celula[i].flui.fracMol[kfrac];
            fluC[i].atualizaPropCompStandard();
            fluC[i].atualizaPropComp(celula[i].pres, celula[i].temp, fluC[i].dCalculatedBeta,
                                     fluC[i].oCalculatedLiqComposition,
                                     fluC[i].oCalculatedVapComposition, arq.pocinjec);
        }
        celula[i].nMol = (fluC[i].MasEspoleo(celula[i].pres, celula[i].temp) * (1 - celula[i].alf) +
                          fluC[i].MasEspGas(celula[i].pres, celula[i].temp) *
                              celula[i].alf) *
                         celula[i].duto.area / pesoMolC;
        fluC[i].Pmol = pesoMolC;
        double hol = 1. - celula[i].alf;
        double bet = celula[i].bet;

        double rs = celula[i].flui.RS(celula[i].pres, celula[i].temp);
        double bo = celula[i].flui.BOFunc(celula[i].pres, celula[i].temp, rs);
        double ba = celula[i].flui.BAFunc(celula[i].pres, celula[i].temp);
        double bsw = celula[i].flui.BSW * ba / (bo + ba * celula[i].flui.BSW - celula[i].flui.BSW * bo);

        double BSWini = celula[i].flui.BSW;
        double BSWF = 0.;
        double viscLini = 0 * 30 + 1 * celula[i].flui.VisOM(tL);
        double viscHini = 0 * 20 + 1 * celula[i].flui.VisOM(tH);
        double viscLF = viscLini;
        double viscHF = viscHini;

        if (celula[i].acsr.tipo == 2) {
            double rsF = celula[i].acsr.injl.FluidoPro.RS(celula[i].pres, celula[i].temp) * (6.29 / 35.31467);
            if (arq.nfluP > 1 || (*vg1dSP).chaverede != 0) {
                double bswaux = celula[i].acsr.injl.FluidoPro.BSW;
                double contrabsw = 1. - bswaux;
                if (contrabsw < (*vg1dSP).localtiny)
                    contrabsw = 0.9 * (*vg1dSP).localtiny;
                double rhoPSTF;
                if (contrabsw > (*vg1dSP).localtiny)
                    rhoPSTF = (1000 * 141.5 / (131.5 + celula[i].acsr.injl.FluidoPro.API)) + (bswaux / contrabsw) * 1000 * celula[i].acsr.injl.FluidoPro.Denag + celula[i].acsr.injl.FluidoPro.Deng * 1.225 * rsF;
                else
                    rhoPSTF = 1000 * celula[i].acsr.injl.FluidoPro.Denag;
                if (contrabsw > (*vg1dSP).localtiny) {
                    fonteP *= (1. / rhoPSTF);
                    fonteA *= ((bswaux / contrabsw) / rhoPSTF);
                } else {
                    fonteP = 0.;
                    fonteA *= (1 / (celula[i].acsr.injl.FluidoPro.BSW * 1000 * celula[i].acsr.injl.FluidoPro.Denag));
                }
                BSWF = celula[i].acsr.injl.FluidoPro.BSW;
                viscLF = 0 * 30 + 1 * celula[i].acsr.injl.FluidoPro.VisOM(tL);
                viscHF = 0 * 20 + 1 * celula[i].acsr.injl.FluidoPro.VisOM(tH);
            }
        } else if (celula[i].acsr.tipo == 10) {
            double rsF = celula[i].acsr.injm.FluidoPro.RS(celula[i].pres, celula[i].temp) * (6.29 / 35.31467);
            if (arq.nfluP > 1 || (*vg1dSP).chaverede != 0) {
                double bswaux = celula[i].acsr.injm.FluidoPro.BSW;
                double contrabsw = 1. - bswaux;
                if (contrabsw < (*vg1dSP).localtiny)
                    contrabsw = 0.9 * (*vg1dSP).localtiny;
                double rhoPSTF;
                if (contrabsw > (*vg1dSP).localtiny)
                    rhoPSTF = (1000 * 141.5 / (131.5 + celula[i].acsr.injm.FluidoPro.API)) +
                              (bswaux / contrabsw) * 1000 *
                                  celula[i].acsr.injm.FluidoPro.Denag +
                              celula[i].acsr.injm.FluidoPro.Deng * 1.225 * rsF;
                else
                    rhoPSTF = 1000 * celula[i].acsr.injm.FluidoPro.Denag;
                if (contrabsw > (*vg1dSP).localtiny) {
                    fonteP *= (1. / rhoPSTF);
                    fonteA *= ((bswaux / contrabsw) / rhoPSTF);
                } else {
                    fonteP = 0.;
                    fonteA *= (1 / (celula[i].acsr.injm.FluidoPro.BSW * 1000 * celula[i].acsr.injm.FluidoPro.Denag));
                }
                BSWF = celula[i].acsr.injm.FluidoPro.BSW;
                viscLF = 0 * 30 + 1 * celula[i].acsr.injm.FluidoPro.VisOM(tL);
                viscHF = 0 * 20 + 1 * celula[i].acsr.injm.FluidoPro.VisOM(tH);
            }
        } else if (celula[i].acsr.tipo == 3) {
            double rsF = celula[i].acsr.ipr.FluidoPro.RS(celula[i].pres, celula[i].temp) * (6.29 / 35.31467);
            if (arq.nfluP > 1 || (*vg1dSP).chaverede != 0) {
                double bswaux = celula[i].acsr.ipr.FluidoPro.BSW;
                double contrabsw = 1. - bswaux;
                if (contrabsw <= (*vg1dSP).localtiny)
                    contrabsw = 0.9 * (*vg1dSP).localtiny;
                double rhoPSTF;
                if (contrabsw > 0)
                    rhoPSTF = (1000 * 141.5 / (131.5 + celula[i].acsr.ipr.FluidoPro.API)) + (bswaux / contrabsw) * 1000 * celula[i].acsr.ipr.FluidoPro.Denag + celula[i].acsr.ipr.FluidoPro.Deng * 1.225 * rsF;
                else
                    rhoPSTF = 1000 * celula[i].acsr.ipr.FluidoPro.Denag;
                if (contrabsw > (*vg1dSP).localtiny) {
                    fonteP *= (1. / rhoPSTF);
                    fonteA *= ((bswaux / contrabsw) / rhoPSTF);
                } else {
                    fonteP = 0.;
                    fonteA *= (1 / (celula[i].acsr.ipr.FluidoPro.BSW * 1000 * celula[i].acsr.ipr.FluidoPro.Denag));
                }
                BSWF = celula[i].acsr.ipr.FluidoPro.BSW;
                viscLF = 0 * 30 + 1 * celula[i].acsr.ipr.FluidoPro.VisOM(tL);
                viscHF = 0 * 20 + 1 * celula[i].acsr.ipr.FluidoPro.VisOM(tH);
            }
        } else if (celula[i].acsr.tipo == 9) {

            double rsF = fluF.RS(celula[i].pres, celula[i].temp) * (6.29 / 35.31467);
            if (arq.nfluP > 1 || (*vg1dSP).chaverede != 0) {
                double bswaux = fluF.BSW;
                double contrabsw = 1. - bswaux;
                if (contrabsw <= (*vg1dSP).localtiny)
                    contrabsw = 0.9 * (*vg1dSP).localtiny;
                double rhoPSTF;
                if (contrabsw > (*vg1dSP).localtiny)
                    rhoPSTF = (1000 * 141.5 / (131.5 + fluF.API)) + (bswaux / contrabsw) * 1000 * fluF.Denag + fluF.Deng * 1.225 * rsF;
                else
                    rhoPSTF = 1000 * fluF.Denag;
                if (contrabsw > (*vg1dSP).localtiny) {
                    fonteP *= (1. / rhoPSTF);
                    fonteA *= ((bswaux / contrabsw) / rhoPSTF);
                } else {
                    fonteP = 0.;
                    fonteA *= (1 / (fluF.BSW * 1000 * fluF.Denag));
                }
                BSWF = fluF.BSW;
                viscLF = 0 * 30 + 1 * fluF.VisOM(tL);
                viscHF = 0 * 20 + 1 * fluF.VisOM(tH);
            }
        } else if (celula[i].acsr.tipo == 15) {
            double rsF = celula[i].acsr.radialPoro.flup.RS(celula[i].pres, celula[i].temp) * (6.29 / 35.31467);
            if (arq.nfluP > 1 || (*vg1dSP).chaverede != 0) {
                double bswaux = celula[i].acsr.radialPoro.BSW;
                double contrabsw = 1. - bswaux;
                if (contrabsw <= (*vg1dSP).localtiny)
                    contrabsw = 0.9 * (*vg1dSP).localtiny;
                double rhoPSTF;
                if (contrabsw > 0)
                    rhoPSTF = (1000 * 141.5 / (131.5 + celula[i].acsr.radialPoro.flup.API)) + (bswaux / contrabsw) * 1000 * celula[i].acsr.radialPoro.flup.Denag + celula[i].acsr.radialPoro.flup.Deng * 1.225 * rsF;
                else
                    rhoPSTF = 1000 * celula[i].acsr.radialPoro.flup.Denag;
                if (contrabsw > (*vg1dSP).localtiny) {
                    fonteP *= (1. / rhoPSTF);
                    fonteA *= ((bswaux / contrabsw) / rhoPSTF);
                } else {
                    fonteP = 0.;
                    fonteA *= (1 / (celula[i].acsr.radialPoro.BSW * 1000 * celula[i].acsr.radialPoro.flup.Denag));
                }
                BSWF = celula[i].acsr.radialPoro.BSW;
                viscLF = 0 * 30 + 1 * celula[i].acsr.radialPoro.flup.VisOM(tL);
                viscHF = 0 * 20 + 1 * celula[i].acsr.radialPoro.flup.VisOM(tH);
            }
        } else if (celula[i].acsr.tipo == 16) {
            double rsF = celula[i].acsr.poroso2D.dados.flup.RS(celula[i].pres, celula[i].temp) * (6.29 / 35.31467);
            if (arq.nfluP > 1 || (*vg1dSP).chaverede != 0) {
                double bswaux = celula[i].acsr.poroso2D.dados.transfer.BSW;
                double contrabsw = 1. - bswaux;
                if (contrabsw <= (*vg1dSP).localtiny)
                    contrabsw = 0.9 * (*vg1dSP).localtiny;
                double rhoPSTF;
                if (contrabsw > 0)
                    rhoPSTF = (1000 * 141.5 / (131.5 + celula[i].acsr.poroso2D.dados.flup.API)) + (bswaux / contrabsw) * 1000 * celula[i].acsr.poroso2D.dados.flup.Denag + celula[i].acsr.poroso2D.dados.flup.Deng * 1.225 * rsF;
                else
                    rhoPSTF = 1000 * celula[i].acsr.poroso2D.dados.flup.Denag;
                if (contrabsw > (*vg1dSP).localtiny) {
                    fonteP *= (1. / rhoPSTF);
                    fonteA *= ((bswaux / contrabsw) / rhoPSTF);
                } else {
                    fonteP = 0.;
                    fonteA *= (1 / (celula[i].acsr.poroso2D.dados.transfer.BSW * 1000 * celula[i].acsr.poroso2D.dados.flup.Denag));
                }
                BSWF = celula[i].acsr.poroso2D.dados.transfer.BSW;
                viscLF = 0 * 30 + 1 * celula[i].acsr.poroso2D.dados.flup.VisOM(tL);
                viscHF = 0 * 20 + 1 * celula[i].acsr.poroso2D.dados.flup.VisOM(tH);
            }
        } else if ((fabs(fonteO) > (*vg1dSP).localtiny && celula[i].acsr.tipo != 2 && celula[i].acsr.tipo != 3 &&
                    celula[i].acsr.tipo != 9 && celula[i].acsr.tipo != 15 && celula[i].acsr.tipo != 16) ||
                   (fabs(fonteG) > (*vg1dSP).localtiny && celula[i].acsr.tipo != 1 && celula[i].acsr.tipo != 2 && celula[i].acsr.tipo != 3 && celula[i].acsr.tipo != 9 && celula[i].acsr.tipo != 15 && celula[i].acsr.tipo != 16)) {
            if (celula[i].acsr.tipo == 5 || celula[i].acsr.tipo == 8) {
                double rsF = celula[i].flui.RS(celula[i].pres, celula[i].temp) * (6.29 / 35.31467);
                if (arq.nfluP > 1 || (*vg1dSP).chaverede != 0) {
                    double rhoPSTF = (1 - celula[i].flui.BSW) * (1000 * 141.5 / (131.5 + celula[i].flui.API)) + celula[i].flui.BSW * 1000 * celula[i].flui.Denag;
                    fonteP *= ((1 - celula[i].flui.BSW) / rhoPSTF);
                    fonteA *= (celula[i].flui.BSW / rhoPSTF);
                    BSWF = celula[i].flui.BSW;
                    viscLF = 0 * 30 + 1 * celula[i].flui.VisOM(tL);
                    viscHF = 0 * 20 + 1 * celula[i].flui.VisOM(tH);
                }
            } else if ((*celula[i].acsrL).tipo == 5 || (*celula[i].acsrL).tipo == 8) {
                if (i > 0) {
                    double rsF = celula[i - 1].flui.RS(celula[i].pres, celula[i].temp) * (6.29 / 35.31467);
                    if (arq.nfluP > 1 || (*vg1dSP).chaverede != 0) {
                        double rhoPSTF = (1 - celula[i - 1].flui.BSW) * (1000 * 141.5 / (131.5 + celula[i - 1].flui.API)) +
                                         celula[i - 1].flui.BSW * 1000 *
                                             celula[i - 1].flui.Denag;
                        fonteP *= ((1 - celula[i - 1].flui.BSW) / rhoPSTF);
                        fonteA *= (celula[i - 1].flui.BSW / rhoPSTF);
                        BSWF = celula[i - 1].flui.BSW;
                        viscLF = 0 * 30 + 1 * celula[i - 1].flui.VisOM(tL);
                        viscHF = 0 * 20 + 1 * celula[i - 1].flui.VisOM(tH);
                    }
                } else {
                    double rsF = celula[i].flui.RS(celula[i].pres, celula[i].temp) * (6.29 / 35.31467);
                    if (arq.nfluP > 1 || (*vg1dSP).chaverede != 0) {
                        double rhoPSTF = (1 - celula[i].flui.BSW) * (1000 * 141.5 / (131.5 + celula[i].flui.API)) +
                                         celula[i].flui.BSW * 1000 *
                                             celula[i].flui.Denag;
                        fonteP *= ((1 - celula[i].flui.BSW) / rhoPSTF);
                        fonteA *= (celula[i].flui.BSW / rhoPSTF);
                        BSWF = celula[i].flui.BSW;
                        viscLF = 0 * 30 + 1 * celula[i].flui.VisOM(tL);
                        viscHF = 0 * 20 + 1 * celula[i].flui.VisOM(tH);
                    }
                }
            }
        }
        if (titF > 1. - 1e-15) {
            fonteO = 0.;
            fonteP = 0.;
            fonteA = 0.;
        }

        double volpesFim = hol * (1 - bet) * (1 - bsw) / bo;
        double volaguaFim = hol * (1 - bet) * bsw; // nao deveria ser dividido por Bo???????????
        double MultPe;
        double MultPd = 0.;
        double residuoP;
        double MultAe = 0.;
        double MultAd;
        double residuoA;
        if (arq.nfluP > 1 || (*vg1dSP).chaverede != 0) {
            MultPe = 0.;
            if (celula[i].flui.dStockTankVaporMassFraction < 1. - 1e-15)
                MultPe = celula[i].QL * (1 - betI0) * (1 - bsw0) / bo0;
            MultPd = 0.;
            if (celula[i + 1].flui.dStockTankVaporMassFraction < 1. - 1e-15)
                MultPd = celula[i + 1].QL * (1 - betI1) * (1 - bsw1) / bo1;
            residuoP = (volpesFim - celula[i].VolPesaST) * at / dt + (MultPd - MultPe) / dx - fonteP / dx;
            MultAe = 0.;
            if (celula[i].flui.dStockTankVaporMassFraction < 1. - 1e-15)
                MultAe = celula[i].QL * (1 - betI0) * bsw0 / bo0;
            MultAd = 0.;
            if (celula[i + 1].flui.dStockTankVaporMassFraction < 1. - 1e-15)
                MultAd = celula[i + 1].QL * (1 - betI1) * bsw1 / bo1;
            residuoA = (volaguaFim - celula[i].VolAguaST) * at / dt + (MultAd - MultAe) / dx - fonteA / dx;
        }
        if ((arq.nfluP > 1 && (arq.flashCompleto == 0 || celula[i].flui.viscBlackOil == 1)) || (*vg1dSP).chaverede != 0) {
            if (volpesFim > 1e-3 && (fonteP > 0 || (MultPd < 0 || MultPe > 0))) {
                VISCL[i] = (dt * (viscLF * fonteP / dx + viscLini * residuoP - (viscL1 * MultPd - viscL0 * MultPe) / dx) + viscLini * celula[i].VolPesaST * at) / (volpesFim * at);
                VISCH[i] = (dt * (viscHF * fonteP / dx + viscHini * residuoP - (viscH1 * MultPd - viscH0 * MultPe) / dx) + viscHini * celula[i].VolPesaST * at) / (volpesFim * at);
            } else {
                VISCL[i] = viscLini;
                VISCH[i] = viscHini;
            }
        }
        if (arq.nfluP > 1 || (*vg1dSP).chaverede != 0) {
            if ((volaguaFim + volpesFim) > 1e-3 && ((fonteP > 0 || fonteA > 0) ||
                                                    ((MultPd < 0 || MultPe > 0) || (MultAd < 0 || MultAe > 0))))
                BSW[i] = (dt * (BSWF * (fonteA + fonteP) / dx + BSWini * (residuoA + residuoP) - (BSW1 * (MultAd + MultPd) - BSW0 * (MultAe + MultPe)) / dx) + BSWini * (celula[i].VolAguaST + celula[i].VolPesaST) * at) / ((volaguaFim + volpesFim) * at);
            else
                BSW[i] = BSWini;
        }
        celula[i].VolPesaST = volpesFim;
        celula[i].VolAguaST = volaguaFim;
    }
    for (int i = imin; i <= ncel - 1; i++) {
        if (arq.nfluP > 1 || (*vg1dSP).chaverede != 0) {
            celula[i].flui = fluC[i];
            celula[i].flui.BSW = BSW[i];
            if (arq.flashCompleto == 0 || celula[i].flui.viscBlackOil == 1) {
                celula[i].flui.LVisL = VISCL[i];
                celula[i].flui.LVisH = VISCH[i];
                celula[i].flui.TempL = tL;
                celula[i].flui.TempH = tH;
            }
        }
        if (arq.flashCompleto == 0)
            celula[i].flui.RenovaFluido();
        corrDeng(i);
    }

    if ((*vg1dSP).chaverede == 0 || noextremo == 1 || celula[ncel].Mliqini > (*vg1dSP).localtiny) {
        if (arq.nfluP > 1 || (*vg1dSP).chaverede != 0) {
            celula[ncel].flui = celula[ncel - 1].flui;
            celula[ncel].flui.BSW = BSW[ncel - 1];
            if (arq.flashCompleto == 0) {
                celula[ncel].flui.LVisL = VISCL[ncel - 1];
                celula[ncel].flui.LVisH = VISCH[ncel - 1];
                celula[ncel].flui.TempL = tL;
                celula[ncel].flui.TempH = tH;
            }
        }
    } else {
        celula[ncel].flui = fluiRev;
        celula[ncel].flui.BSW = fluiRev.BSW;
        celula[ncel].flui.LVisL = fluiRev.LVisL;
        celula[ncel].flui.LVisH = fluiRev.LVisH;
        celula[ncel].flui.TempL = fluiRev.TempL;
        celula[ncel].flui.TempH = fluiRev.TempH;
    }
    delete[] fluC;
}

void SProd::renovaFracMol2(ProFlu fluiRev) {
    Vcr<double> BSW(ncel);
    Vcr<double> denag(ncel);
    Vcr<double> VISCL(ncel);
    Vcr<double> VISCH(ncel);
    int ncomp = arq.npseudo;
    ProFlu *fluC;
    fluC = new ProFlu[ncel];
    Vcr<double> fracMol0(ncomp);
    Vcr<double> fracMol1(ncomp);
    Vcr<double> fracMol0O(ncomp);
    Vcr<double> fracMol1O(ncomp);
    Vcr<double> fracMol0G(ncomp);
    Vcr<double> fracMol1G(ncomp);
    Vcr<double> fracMolF(ncomp);
    if (celula[0].acsr.tipo == 15 && (celula[0].fontemassCR + celula[0].fontemassGR + celula[0].fontemassLR) > 0.)
        celula[0].flui.BSW = celula[0].acsr.radialPoro.BSW;
    else if (celula[0].acsr.tipo == 16 && (celula[0].fontemassCR + celula[0].fontemassGR + celula[0].fontemassLR) > 0.)
        celula[0].flui.BSW = celula[0].acsr.poroso2D.dados.transfer.BSW;
    double dt = celula[1].dt;
    double tL = 0.;
    double tH = 70.;
    if (kontaRenovaComp == arq.miniTabAtraso  && arq.miniTabAtraso > 0)
        (*vg1dSP).modoTransiente = 0;
    int imin = 1;
    if ((*vg1dSP).chaverede != 0 && (celula[0].acsr.tipo == 10 || arq.ConContEntrada > 0))
        imin = 0;

    auto normalizarFracoes = [](vector<double>& fracMolFase, double* fracMolOriginal, int npseudo) {
        for (int kfrac = 0; kfrac < npseudo; kfrac++) {
        	fracMolFase[kfrac] = fracMolOriginal[kfrac];
        }
    	// Encontrar o menor valor
        double menorFracFase = fracMolFase[0];
        for (int kfrac = 1; kfrac < npseudo; kfrac++) {
            if (menorFracFase > fracMolFase[kfrac]) {
                menorFracFase = fracMolFase[kfrac];
            }
        }
        menorFracFase=0.;

        // Subtrair o menor valor
        for (int kfrac = 0; kfrac < npseudo; kfrac++) {
        	fracMolFase[kfrac] -= menorFracFase;
        }

        // Calcular soma total
        double fracTotFase = 0.;
        for (int kfrac = 0; kfrac < npseudo; kfrac++) {
            fracTotFase += fracMolFase[kfrac];
        }

        // Normalizar
        if (fracTotFase > 0) {
        	for (int kfrac = 0; kfrac < npseudo; kfrac++) {
        		fracMolFase[kfrac] /= fracTotFase;
            }
        }
    };

    for (int i = imin; i < ncel; i++) {

        double pesoMol0;
        double pesoMol1;
        double pesoMol0O;
        double pesoMol1O;
        double pesoMol0G;
        double pesoMol1G;
        double pesoMolF;
        double pesoMolC;
        double at = celula[i].duto.area;
        double dx = celula[i].dx;
        double ti0;
        if (celula[i].VTemper < 0.)
            ti0 = celula[i].temp;
        else {
            if (i > 0)
                ti0 = celula[i - 1].temp;
            else if (arq.ConContEntrada == 1)
                ti0 = tempE;
            else
                ti0 = celula[i].temp;
        }
        double ti1 = celula[i].temp;
        if (celula[i + 1].VTemper < 0.)
            ti1 = celula[i + 1].temp;

        double titF = 0.;
        ProFlu fluF;

        double boF = 1.;
        double baF = 1.;
        double fwF = 1.;
        double rhoOF = 900.;
        double rhoWF = 1000.;

        double betI0;
        double bo0;
        double ba0;
        double bsw0;
        double BSW0;
        double denag0;
        double viscL0;
        double viscH0;
        if (i > 0 || arq.ConContEntrada == 0) {
            if (i > 0 && celula[i].QG >= 0.)
                betI0 = celula[i - 1].betPigD;
            else
                betI0 = celula[i].betPigE;
        } else {
            if (celula[i].QG >= 0.)
                betI0 = betaE;
            else
                betI0 = celula[i].betPigE;
        }
        if ((i > 0 || arq.ConContEntrada == 1) && celula[i].QL >= 0.) {
            double pe;
            double te;
            if (i > 0) {
                pe = celula[i - 1].pres;
                te = celula[i - 1].temp;
            } else {
                pe = presE;
                te = tempE;
            }
            if (arq.ConContEntrada == 0)
                betI0 = celula[i - 1].betPigD; // testeBeta
            else
                betI0 = betaE; // testeBeta
            double rs0 = (*celula[i].fluiL).RS(pe, te);
            bo0 = (*celula[i].fluiL).BOFunc(pe, te, rs0);
            ba0 = (*celula[i].fluiL).BAFunc(pe, te);
            bsw0 = (*celula[i].fluiL).BSW * ba0 / (bo0 + ba0 * (*celula[i].fluiL).BSW - (*celula[i].fluiL).BSW * bo0);
            BSW0 = (*celula[i].fluiL).BSW;
            denag0 = (*celula[i].fluiL).Denag;
            viscL0 = 0 * 30 + 1 * (*celula[i].fluiL).VisOM(tL);
            viscH0 = 0 * 20 + 1 * (*celula[i].fluiL).VisOM(tH);
        } else {
            betI0 = celula[i].betPigE; // testebeta
            double rs0 = celula[i].flui.RS(celula[i].pres, celula[i].temp);
            bo0 = celula[i].flui.BOFunc(celula[i].pres, celula[i].temp, rs0);
            ba0 = celula[i].flui.BAFunc(celula[i].pres, celula[i].temp);
            bsw0 = celula[i].flui.BSW * ba0 / (bo0 + ba0 * celula[i].flui.BSW - celula[i].flui.BSW * bo0);
            BSW0 = celula[i].flui.BSW;
            denag0 = celula[i].flui.Denag;
            viscL0 = 0 * 30 + 1 * celula[i].flui.VisOM(tL);
            viscH0 = 0 * 20 + 1 * celula[i].flui.VisOM(tH);
        }
        if (bo0 < 1e-15)
            bo0 = 1e-15;

        double betI1 = celula[i].betPigD;
        double rs1 = celula[i].flui.RS(celula[i].pres, celula[i].temp);
        double bo1 = celula[i].flui.BOFunc(celula[i].pres, celula[i].temp, rs1);
        double ba1 = celula[i].flui.BAFunc(celula[i].pres, celula[i].temp);
        double bsw1 = celula[i].flui.BSW * ba1 / (bo1 + ba1 * celula[i].flui.BSW - celula[i].flui.BSW * bo1);
        double BSW1 = celula[i].flui.BSW;
        double denag1 = celula[i].flui.Denag;
        double viscL1 = 0 * 30 + 1 * celula[i].flui.VisOM(tL);
        double viscH1 = 0 * 20 + 1 * celula[i].flui.VisOM(tH);

        // betI1 = celula[i + 1].betPigE;    //duvidabeta
        if (celula[i + 1].QL < 0.) {
            betI1 = celula[i + 1].betPigE; // testeBet
            rs1 = celula[i + 1].flui.RS(celula[i + 1].pres, celula[i + 1].temp);
            bo1 = celula[i + 1].flui.BOFunc(celula[i + 1].pres, celula[i + 1].temp, rs1);
            ba1 = celula[i + 1].flui.BAFunc(celula[i + 1].pres, celula[i + 1].temp);
            bsw1 = celula[i + 1].flui.BSW * ba1 / (bo1 + ba1 * celula[i + 1].flui.BSW - celula[i + 1].flui.BSW * bo1);
            BSW1 = celula[i + 1].flui.BSW;
            denag1 = celula[i + 1].flui.Denag;
            viscL1 = 0 * 30 + 1 * celula[i + 1].flui.VisOM(tL);
            viscH1 = 0 * 20 + 1 * celula[i + 1].flui.VisOM(tH);
        }
        if (bo1 < 1e-15)
            bo1 = 1e-15;
        double titV0;
        double titV1;
        double vazMasLiq0 = celula[i].Mliqini;
        double vazMasLiq1 = celula[i + 1].Mliqini;
        double vazMasGas0 = celula[i].MC - celula[i].Mliqini;
        double vazMasGas1 = celula[i + 1].MC - celula[i + 1].Mliqini;
        pesoMolC = 0.;
        for (int kfrac = 0; kfrac < ncomp; kfrac++) {
            pesoMolC += celula[i].flui.masMol[kfrac] * celula[i].flui.fracMol[kfrac];
        }

        /*double menorFracFase = 0.;
        for (int kfrac = 0; kfrac < ncomp; kfrac++) {
            if (menorFracFase > fluC[i].fracMol[kfrac]) {
                menorFracFase = fluC[i].fracMol[kfrac];
            }
        }
        double fracTotFase = 0.;
        for (int kfrac = 0; kfrac < ncomp; kfrac++)
            fluC[i].fracMol[kfrac] -= menorFracFase;
        for (int kfrac = 0; kfrac < ncomp; kfrac++)
            fracTotFase += fluC[i].fracMol[kfrac];
        for (int kfrac = 0; kfrac < ncomp; kfrac++)
            fluC[i].fracMol[kfrac] /= fracTotFase;*/

        if ((i > 0 || arq.ConContEntrada == 1) && celula[i].Mliqini >= 0.) {
            double pe;
            double te;
            if (i > 0) {
                pe = celula[i - 1].pres;
                te = celula[i - 1].temp;
            } else {
                pe = presE;
                te = tempE;
            }
            double fwV = bsw0;
            double rhoOV = (*celula[i].fluiL).MasEspoleo(pe, te);
            double rhoWV = (*celula[i].fluiL).MasEspAgua(pe, te);
            double titLocal=(*celula[i].fluiL).FracMass(pe, te);
            titV0 = (1 - fwV) * rhoOV / ((1 - fwV) * rhoOV + fwV * rhoWV);
            vazMasLiq0 -= betI0 * celula[i].QL * celula[i].fluicol.MasEspFlu(pe, te);
            vazMasLiq0 *= titV0;
            pesoMol0O = 0;
            vector<double> fracMolFase(ncomp);
            normalizarFracoes(fracMolFase, (*celula[i].fluiL).oCalculatedLiqComposition, ncomp);
            for (int kfrac = 0; kfrac < ncomp; kfrac++) {
                fracMol0O[kfrac] = fracMolFase[kfrac];
                pesoMol0O += (*celula[i].fluiL).masMol[kfrac] * fracMolFase[kfrac];
            }
            if(pesoMol0O<1e-3 || titLocal>1.-1e-3 || celula[i].alfL>1.-1e-3){
                pesoMol0O = 0;
                for (int kfrac = 0; kfrac < ncomp; kfrac++) {
                    fracMol0O[kfrac] = (*celula[i].fluiL).fracMol[kfrac];
                    pesoMol0O += (*celula[i].fluiL).masMol[kfrac] * (*celula[i].fluiL).fracMol[kfrac];
                }
            }
        } else {
            double fwV = bsw0;
            double rhoOV = celula[i].flui.MasEspoleo(celula[i].pres, celula[i].temp);
            double rhoWV = celula[i].flui.MasEspAgua(celula[i].pres, celula[i].temp);
            double titLocal=celula[i].flui.FracMass(celula[i].pres, celula[i].temp);
            titV0 = (1 - fwV) * rhoOV / ((1 - fwV) * rhoOV + fwV * rhoWV);
            vazMasLiq0 -= betI0 * celula[i].QL * celula[i].fluicol.MasEspFlu(celula[i].pres, celula[i].temp);
            vazMasLiq0 *= titV0;
            pesoMol0O = 0;
            vector<double> fracMolFase(ncomp);
            normalizarFracoes(fracMolFase, celula[i].flui.oCalculatedLiqComposition, ncomp);
            for (int kfrac = 0; kfrac < ncomp; kfrac++) {
                fracMol0O[kfrac] = fracMolFase[kfrac];
                pesoMol0O += celula[i].flui.masMol[kfrac] * fracMolFase[kfrac];
            }
            if(pesoMol0O<1e-3 || titLocal>1.-1e-3 || celula[i].alf>1.-1e-3){
                pesoMol0O = 0;
                for (int kfrac = 0; kfrac < ncomp; kfrac++) {
                    fracMol0O[kfrac] = celula[i].flui.fracMol[kfrac];
                    pesoMol0O += celula[i].flui.masMol[kfrac] * celula[i].flui.fracMol[kfrac];
                }
            }
        }
        if (celula[i + 1].Mliqini >= 0.) {
            double fwV = bsw1;
            double rhoOV = celula[i].flui.MasEspoleo(celula[i].pres, celula[i].temp);
            double rhoWV = celula[i].flui.MasEspAgua(celula[i].pres, celula[i].temp);
            double titLocal=celula[i].flui.FracMass(celula[i].pres, celula[i].temp);
            titV1 = (1 - fwV) * rhoOV / ((1 - fwV) * rhoOV + fwV * rhoWV);
            vazMasLiq1 -= betI1 * celula[i].QL * celula[i].fluicol.MasEspFlu(celula[i].pres, celula[i].temp);
            vazMasLiq1 *= titV1;
            pesoMol1O = 0;
            vector<double> fracMolFase(ncomp);
            normalizarFracoes(fracMolFase, celula[i].flui.oCalculatedLiqComposition, ncomp);
            for (int kfrac = 0; kfrac < ncomp; kfrac++) {
                fracMol1O[kfrac] = fracMolFase[kfrac];
                pesoMol1O += celula[i].flui.masMol[kfrac] * fracMolFase[kfrac];
            }
            if(pesoMol1O<1e-3 || titLocal>1.-1e-3 || celula[i].alf>1-1e-3){
                pesoMol1O = 0;
                for (int kfrac = 0; kfrac < ncomp; kfrac++) {
                    fracMol1O[kfrac] = celula[i].flui.fracMol[kfrac];
                    pesoMol1O += celula[i].flui.masMol[kfrac] * celula[i].flui.fracMol[kfrac];
                }
            }
        } else {
            double fwV = bsw1;
            double rhoOV = celula[i + 1].flui.MasEspoleo(celula[i + 1].pres, celula[i + 1].temp);
            double rhoWV = celula[i + 1].flui.MasEspAgua(celula[i + 1].pres, celula[i + 1].temp);
            double titLocal=celula[i+1].flui.FracMass(celula[i+1].pres, celula[i+1].temp);
            titV1 = (1 - fwV) * rhoOV / ((1 - fwV) * rhoOV + fwV * rhoWV);
            vazMasLiq1 -= betI1 * celula[i + 1].QL * celula[i].fluicol.MasEspFlu(celula[i + 1].pres, celula[i + 1].temp);
            vazMasLiq1 *= titV1;
            pesoMol1O = 0;
            vector<double> fracMolFase(ncomp);
            normalizarFracoes(fracMolFase, celula[i + 1].flui.oCalculatedLiqComposition, ncomp);
            for (int kfrac = 0; kfrac < ncomp; kfrac++) {
                fracMol1O[kfrac] = fracMolFase[kfrac];
                pesoMol1O += celula[i + 1].flui.masMol[kfrac] * fracMolFase[kfrac];
            }
            if(pesoMol1O<1e-3 || titLocal>1.-1e-3 || celula[i+1].alf>1-1e-3){
                pesoMol1O = 0;
                for (int kfrac = 0; kfrac < ncomp; kfrac++) {
                    fracMol1O[kfrac] = celula[i + 1].flui.fracMol[kfrac];
                    pesoMol1O += celula[i + 1].flui.masMol[kfrac] * celula[i + 1].flui.fracMol[kfrac];
                }
            }
        }

        if ((i > 0 || arq.ConContEntrada == 1) && (celula[i].MC - celula[i].Mliqini) >= 0.) {
            double pe;
            double te;
            if (i > 0) {
                pe = celula[i - 1].pres;
                te = celula[i - 1].temp;
            } else {
                pe = presE;
                te = tempE;
            }
        	pesoMol0G = 0;
            double titLocal=(*celula[i].fluiL).FracMass(pe, te);
            vector<double> fracMolFase(ncomp);
            normalizarFracoes(fracMolFase, (*celula[i].fluiL).oCalculatedVapComposition, ncomp);
            for (int kfrac = 0; kfrac < ncomp; kfrac++) {
                fracMol0G[kfrac] = fracMolFase[kfrac];
                pesoMol0G += (*celula[i].fluiL).masMol[kfrac] * fracMolFase[kfrac];
            }
            if(pesoMol0G<1e-3 || titLocal<1e-3 || celula[i].alfL<1e-3){
                pesoMol0G = 0;
                for (int kfrac = 0; kfrac < ncomp; kfrac++) {
                    fracMol0G[kfrac] = (*celula[i].fluiL).fracMol[kfrac];
                    pesoMol0G += (*celula[i].fluiL).masMol[kfrac] * (*celula[i].fluiL).fracMol[kfrac];
                }
            }
        } else {
            pesoMol0G = 0;
            double titLocal=celula[i].flui.FracMass(celula[i].pres, celula[i].temp);
            vector<double> fracMolFase(ncomp);
            normalizarFracoes(fracMolFase, celula[i].flui.oCalculatedVapComposition, ncomp);
            for (int kfrac = 0; kfrac < ncomp; kfrac++) {
                fracMol0G[kfrac] = fracMolFase[kfrac];
                pesoMol0G += celula[i].flui.masMol[kfrac] * fracMolFase[kfrac];
            }
            if(pesoMol0G<1e-3 || titLocal<1e-3 || celula[i].alf<1e-3){
                pesoMol0G = 0;
                for (int kfrac = 0; kfrac < ncomp; kfrac++) {
                    fracMol0G[kfrac] = celula[i].flui.fracMol[kfrac];
                    pesoMol0G += celula[i].flui.masMol[kfrac] * celula[i].flui.fracMol[kfrac];
                }
            }
        }
        if ((celula[i + 1].MC - celula[i + 1].Mliqini) >= 0.) {
            pesoMol1G = 0;
            double titLocal=celula[i].flui.FracMass(celula[i].pres, celula[i].temp);
            vector<double> fracMolFase(ncomp);
            normalizarFracoes(fracMolFase, celula[i].flui.oCalculatedVapComposition, ncomp);
            for (int kfrac = 0; kfrac < ncomp; kfrac++) {
                fracMol1G[kfrac] = fracMolFase[kfrac];
                pesoMol1G += celula[i].flui.masMol[kfrac] * fracMolFase[kfrac];
            }
            if(pesoMol1G<1e-3 || titLocal<1e-3 || celula[i].alf<1e-3){
                pesoMol1G = 0;
                for (int kfrac = 0; kfrac < ncomp; kfrac++) {
                    fracMol1G[kfrac] = celula[i].flui.fracMol[kfrac];
                    pesoMol1G += celula[i].flui.masMol[kfrac] * celula[i].flui.fracMol[kfrac];
                }
            }
        } else {
            pesoMol1G = 0;
            double titLocal=celula[i+1].flui.FracMass(celula[i+1].pres, celula[i+1].temp);
            vector<double> fracMolFase(ncomp);
            normalizarFracoes(fracMolFase, celula[i + 1].flui.oCalculatedVapComposition, ncomp);
            for (int kfrac = 0; kfrac < ncomp; kfrac++) {
                fracMol1G[kfrac] = fracMolFase[kfrac];
                pesoMol1G += celula[i + 1].flui.masMol[kfrac] * fracMolFase[kfrac];
            }
            if(pesoMol1G<1e-3 || titLocal<1e-3 || celula[i+1].alf<1e-3){
                pesoMol1G = 0;
                for (int kfrac = 0; kfrac < ncomp; kfrac++) {
                    fracMol1G[kfrac] = celula[i + 1].flui.fracMol[kfrac];
                    pesoMol1G += celula[i + 1].flui.masMol[kfrac] * celula[i + 1].flui.fracMol[kfrac];
                }
            }
        }

        if ((i > 0 || arq.ConContEntrada == 1) && (vazMasLiq0 + vazMasGas0) >= 0) {
            pesoMol0 = 0;
            for (int kfrac = 0; kfrac < ncomp; kfrac++) {
                fracMol0[kfrac] = (*celula[i].fluiL).fracMol[kfrac];
                pesoMol0 += (*celula[i].fluiL).masMol[kfrac] * (*celula[i].fluiL).fracMol[kfrac];
            }
        } else {
            pesoMol0 = 0;
            for (int kfrac = 0; kfrac < ncomp; kfrac++) {
                fracMol0[kfrac] = celula[i].flui.fracMol[kfrac];
                pesoMol0 += celula[i].flui.masMol[kfrac] * celula[i].flui.fracMol[kfrac];
            }
        }
        if ((vazMasLiq1 + vazMasGas1) >= 0) {
            pesoMol1 = 0;
            for (int kfrac = 0; kfrac < ncomp; kfrac++) {
                fracMol1[kfrac] = celula[i].flui.fracMol[kfrac];
                pesoMol1 += celula[i].flui.masMol[kfrac] * celula[i].flui.fracMol[kfrac];
            }
        } else {
            pesoMol1 = 0;
            for (int kfrac = 0; kfrac < ncomp; kfrac++) {
                fracMol1[kfrac] = celula[i + 1].flui.fracMol[kfrac];
                pesoMol1 += celula[i + 1].flui.masMol[kfrac] * celula[i + 1].flui.fracMol[kfrac];
            }
        }

        double fonteO = celula[i].fontemassLR;
        double fonteG = celula[i].fontemassGR;
        double fonteP = celula[i].fontemassLR;
        double fonteA = celula[i].fontemassLR;

        if (celula[i].acsr.tipo == 1) {
            if (celula[i].acsr.injg.FluidoPro.dCalculatedBeta < 0. || celula[i].acsr.injg.FluidoPro.dCalculatedBeta > 1.)
                celula[i].acsr.injg.FluidoPro.atualizaPropComp(celula[i].pres, celula[i].temp, -1, NULL, NULL, arq.pocinjec);
            else
                celula[i].acsr.injg.FluidoPro.atualizaPropComp(celula[i].pres, celula[i].temp,
                                                               celula[i].acsr.injg.FluidoPro.dCalculatedBeta, celula[i].acsr.injg.FluidoPro.oCalculatedLiqComposition,
                                                               celula[i].acsr.injg.FluidoPro.oCalculatedVapComposition, arq.pocinjec);
            if (celula[i].acsr.injg.QGas > 0.)
                fluF = celula[i].acsr.injg.FluidoPro;
            else
                fluF = celula[i].flui;
            fwF = 0.;
            titF = 1.;
        } else if (celula[i].acsr.tipo == 2) {
            if (celula[i].acsr.injl.FluidoPro.dCalculatedBeta < 0. || celula[i].acsr.injl.FluidoPro.dCalculatedBeta > 1.)
                celula[i].acsr.injl.FluidoPro.atualizaPropComp(celula[i].pres, celula[i].temp, -1, NULL, NULL, arq.pocinjec);
            else
                celula[i].acsr.injl.FluidoPro.atualizaPropComp(celula[i].pres, celula[i].temp,
                                                               celula[i].acsr.injl.FluidoPro.dCalculatedBeta, celula[i].acsr.injl.FluidoPro.oCalculatedLiqComposition,
                                                               celula[i].acsr.injl.FluidoPro.oCalculatedVapComposition, arq.pocinjec);
            if (celula[i].acsr.injl.QLiq > 0.)
                fluF = celula[i].acsr.injl.FluidoPro;
            else
                fluF = celula[i].flui;
            boF = fluF.BOFunc(celula[i].pres, celula[i].temp);
            baF = fluF.BAFunc(celula[i].pres, celula[i].temp);
            fwF = fluF.BSW * baF / (boF + baF * fluF.BSW - fluF.BSW * boF);
            rhoOF = fluF.MasEspoleo(celula[i].pres, celula[i].temp);
            rhoWF = fluF.MasEspAgua(celula[i].pres, celula[i].temp);
            titF = (1 - fwF) * rhoOF / ((1 - fwF) * rhoOF + fwF * rhoWF);
        } else if (celula[i].acsr.tipo == 3) {
            if (celula[i].acsr.ipr.FluidoPro.dCalculatedBeta < 0. || celula[i].acsr.ipr.FluidoPro.dCalculatedBeta > 1.)
                celula[i].acsr.ipr.FluidoPro.atualizaPropComp(celula[i].pres, celula[i].temp, -1, NULL, NULL, arq.pocinjec);
            else
                celula[i].acsr.ipr.FluidoPro.atualizaPropComp(celula[i].pres, celula[i].temp,
                                                              celula[i].acsr.ipr.FluidoPro.dCalculatedBeta, celula[i].acsr.ipr.FluidoPro.oCalculatedLiqComposition,
                                                              celula[i].acsr.ipr.FluidoPro.oCalculatedVapComposition, arq.pocinjec);
            if ((celula[i].acsr.ipr.Pres) > celula[i].pres)
                fluF = celula[i].acsr.ipr.FluidoPro;
            else
                fluF = celula[i].flui;
            boF = fluF.BOFunc(celula[i].pres, celula[i].temp);
            baF = fluF.BAFunc(celula[i].pres, celula[i].temp);
            fwF = fluF.BSW * baF / (boF + baF * fluF.BSW - fluF.BSW * boF);
            rhoOF = fluF.MasEspoleo(celula[i].pres, celula[i].temp);
            rhoWF = fluF.MasEspAgua(celula[i].pres, celula[i].temp);
            titF = (1 - fwF) * rhoOF / ((1 - fwF) * rhoOF + fwF * rhoWF);
        } else if (celula[i].acsr.tipo == 10) {
            if (celula[i].acsr.injm.FluidoPro.dCalculatedBeta < 0. || celula[i].acsr.injm.FluidoPro.dCalculatedBeta > 1.)
                celula[i].acsr.injm.FluidoPro.atualizaPropComp(celula[i].pres, celula[i].temp, -1, NULL, NULL, arq.pocinjec);
            else
                celula[i].acsr.injm.FluidoPro.atualizaPropComp(celula[i].pres, celula[i].temp,
                                                               celula[i].acsr.injm.FluidoPro.dCalculatedBeta, celula[i].acsr.injm.FluidoPro.oCalculatedLiqComposition,
                                                               celula[i].acsr.injm.FluidoPro.oCalculatedVapComposition, arq.pocinjec);
            if ((celula[i].acsr.injm.MassC + celula[i].acsr.injm.MassG + celula[i].acsr.injm.MassP) > 0.)
                fluF = celula[i].acsr.injm.FluidoPro;
            else
                fluF = celula[i].flui;
            boF = fluF.BOFunc(celula[i].pres, celula[i].temp);
            baF = fluF.BAFunc(celula[i].pres, celula[i].temp);
            fwF = fluF.BSW * baF / (boF + baF * fluF.BSW - fluF.BSW * boF);
            rhoOF = fluF.MasEspoleo(celula[i].pres, celula[i].temp);
            rhoWF = fluF.MasEspAgua(celula[i].pres, celula[i].temp);
            titF = (1 - fwF) * rhoOF / ((1 - fwF) * rhoOF + fwF * rhoWF);
        } else if (celula[i].acsr.tipo == 9 && celula[i].acsr.fontechk.abertura > 1e-6) {
            if (celula[i].acsr.fontechk.fluidoP.dCalculatedBeta < 0. || celula[i].acsr.fontechk.fluidoP.dCalculatedBeta > 1.)
                celula[i].acsr.fontechk.fluidoP.atualizaPropComp(celula[i].pres, celula[i].temp, -1, NULL, NULL, arq.pocinjec);
            else
                celula[i].acsr.fontechk.fluidoP.atualizaPropComp(celula[i].pres, celula[i].temp,
                                                                 celula[i].acsr.fontechk.fluidoP.dCalculatedBeta, celula[i].acsr.fontechk.fluidoP.oCalculatedLiqComposition,
                                                                 celula[i].acsr.fontechk.fluidoP.oCalculatedVapComposition, arq.pocinjec);
            if (celula[i].acsr.fontechk.fluidoPamb.dCalculatedBeta < 0. || celula[i].acsr.fontechk.fluidoPamb.dCalculatedBeta > 1.)
                celula[i].acsr.fontechk.fluidoPamb.atualizaPropComp(celula[i].pres, celula[i].temp, -1, NULL, NULL, arq.pocinjec);
            else
                celula[i].acsr.fontechk.fluidoPamb.atualizaPropComp(celula[i].pres, celula[i].temp,
                                                                    celula[i].acsr.fontechk.fluidoPamb.dCalculatedBeta, celula[i].acsr.fontechk.fluidoPamb.oCalculatedLiqComposition,
                                                                    celula[i].acsr.fontechk.fluidoPamb.oCalculatedVapComposition, arq.pocinjec);
            if (celula[i].acsr.fontechk.presT > celula[i].acsr.fontechk.pamb) {
                fluF = celula[i].acsr.fontechk.fluidoP;
            } else {
                fluF = celula[i].acsr.fontechk.fluidoPamb;
            }
            boF = fluF.BOFunc(celula[i].pres, celula[i].temp);
            baF = fluF.BAFunc(celula[i].pres, celula[i].temp);
            fwF = fluF.BSW * baF / (boF + baF * fluF.BSW - fluF.BSW * boF);
            rhoOF = fluF.MasEspoleo(celula[i].pres, celula[i].temp);
            rhoWF = fluF.MasEspAgua(celula[i].pres, celula[i].temp);
            titF = (1 - fwF) * rhoOF / ((1 - fwF) * rhoOF + fwF * rhoWF);
        } else if (celula[i].acsr.tipo == 15) {
            double tRes = celula[i].acsr.radialPoro.tRes;
            if (celula[i].acsr.radialPoro.flup.dCalculatedBeta < 0. || celula[i].acsr.radialPoro.flup.dCalculatedBeta > 1.)
                celula[i].acsr.radialPoro.flup.atualizaPropComp(celula[i].pres, tRes, -1, NULL, NULL, arq.pocinjec);
            else
                celula[i].acsr.radialPoro.flup.atualizaPropComp(celula[i].pres, tRes,
                                                                celula[i].acsr.radialPoro.flup.dCalculatedBeta, celula[i].acsr.radialPoro.flup.oCalculatedLiqComposition,
                                                                celula[i].acsr.radialPoro.flup.oCalculatedVapComposition, arq.pocinjec);
            if ((celula[i].fontemassLR + celula[i].fontemassGR) > 1e-15)
                fluF = celula[i].acsr.radialPoro.flup;
            else
                fluF = celula[i].flui;
            boF = fluF.BOFunc(celula[i].pres, celula[i].temp);
            baF = fluF.BAFunc(celula[i].pres, celula[i].temp);
            fwF = fluF.BSW * baF / (boF + baF * fluF.BSW - fluF.BSW * boF);
            rhoOF = fluF.MasEspoleo(celula[i].pres, celula[i].temp);
            rhoWF = fluF.MasEspAgua(celula[i].pres, celula[i].temp);
            titF = (1 - fwF) * rhoOF / ((1 - fwF) * rhoOF + fwF * rhoWF);
        } else if (celula[i].acsr.tipo == 16) {
            double tRes = celula[i].acsr.poroso2D.dados.tRes;
            if (celula[i].acsr.poroso2D.dados.flup.dCalculatedBeta < 0. || celula[i].acsr.poroso2D.dados.flup.dCalculatedBeta > 1.)
                celula[i].acsr.poroso2D.dados.flup.atualizaPropComp(celula[i].pres, tRes, -1, NULL, NULL, arq.pocinjec);
            else
                celula[i].acsr.poroso2D.dados.flup.atualizaPropComp(celula[i].pres, tRes,
                                                                    celula[i].acsr.poroso2D.dados.flup.dCalculatedBeta, celula[i].acsr.poroso2D.dados.flup.oCalculatedLiqComposition,
                                                                    celula[i].acsr.poroso2D.dados.flup.oCalculatedVapComposition, arq.pocinjec);
            if ((celula[i].fontemassLR + celula[i].fontemassGR) > 1e-15)
                fluF = celula[i].acsr.poroso2D.dados.flup;
            else
                fluF = celula[i].flui;
            boF = fluF.BOFunc(celula[i].pres, celula[i].temp);
            baF = fluF.BAFunc(celula[i].pres, celula[i].temp);
            fwF = fluF.BSW * baF / (boF + baF * fluF.BSW - fluF.BSW * boF);
            rhoOF = fluF.MasEspoleo(celula[i].pres, celula[i].temp);
            rhoWF = fluF.MasEspAgua(celula[i].pres, celula[i].temp);
            titF = (1 - fwF) * rhoOF / ((1 - fwF) * rhoOF + fwF * rhoWF);
        } else {
            fluF = celula[i].flui;
            titF = 0.;
        }

        pesoMolF = 0;
        for (int j = 0; j < ncomp; j++) {
            pesoMolF += fluF.masMol[j] * fluF.fracMol[j];
        }
        fonteO *= titF;
        double tempMol;
        fluC[i] = celula[i].flui;
        double betIV;
        double rsV;
        double boV;
        double baV;
        double bswV;
        double rhoOVol;
        double rhoWVol;
        double titVol;
        int limCorrige = 1;
        if (kontaRenovaComp == arq.miniTabAtraso  || arq.miniTabAtraso == 0)
            limCorrige = 2;
        for (int corrige = 0; corrige < limCorrige; corrige++) {

            betIV = celula[i].bet;
            rsV = fluC[i].RS(celula[i].pres, celula[i].temp);
            boV = fluC[i].BOFunc(celula[i].pres, celula[i].temp, rsV);
            baV = fluC[i].BAFunc(celula[i].pres, celula[i].temp);
            bswV = fluC[i].BSW * baV / (boV + baV * celula[i].flui.BSW - celula[i].flui.BSW * boV);
            rhoOVol = fluC[i].MasEspoleo(celula[i].pres, celula[i].temp);
            rhoWVol = fluC[i].MasEspAgua(celula[i].pres, celula[i].temp);
            titVol = (1 - bswV) * rhoOVol / ((1 - bswV) * rhoOVol + bswV * rhoWVol);

            tempMol = (fluC[i].MasEspLiq(celula[i].pres, celula[i].temp) * (1. - celula[i].alf) * (1. - betIV) * titVol +
                       fluC[i].MasEspGas(celula[i].pres, celula[i].temp) * celula[i].alf) *
                      celula[i].duto.area * celula[i].dx / pesoMolC;

            double menorFrac = 0.;
            for (int kfrac = 0; kfrac < ncomp; kfrac++) {

                /*double vazMol1O;
                double vazMol0O;
                double vazMol1G;
                double vazMol0G;
                vazMol1O=vazMasLiq1*fracMol1O[kfrac]/pesoMol1O;
                vazMol0O=vazMasLiq0*fracMol0O[kfrac]/pesoMol0O;
                vazMol1G=vazMasGas1*fracMol1G[kfrac]/pesoMol1G;
                vazMol0G=vazMasGas0*fracMol0G[kfrac]/pesoMol0G;
            	fluC[i].fracMol[kfrac] = (celula[i].nMolIni * celula[i].flui.fracMol[kfrac] +
                                          ((fonteO + fonteG) * fluF.fracMol[kfrac] / pesoMolF -
                                           (vazMol1O + vazMol1G) +
                                           (vazMol0O + vazMol0G)) *
                                              dt) /
                                         tempMol;*/


                fluC[i].fracMol[kfrac] = (celula[i].nMolIni * celula[i].flui.fracMol[kfrac] +
                                          ((fonteO + fonteG) * fluF.fracMol[kfrac] / pesoMolF -
                                           (vazMasLiq1 + vazMasGas1) * fracMol1[kfrac] / pesoMol1 +
                                           (vazMasLiq0 + vazMasGas0) * fracMol0[kfrac] / pesoMol0) *
                                              dt) /
                                         tempMol;
                if (menorFrac > fluC[i].fracMol[kfrac]) {
                    menorFrac = fluC[i].fracMol[kfrac];
                }
            }
            double fracTot = 0.;
            if(menorFrac<0.){
            	for (int kfrac = 0; kfrac < ncomp; kfrac++)
            		fluC[i].fracMol[kfrac] -= menorFrac;
            }
            for (int kfrac = 0; kfrac < ncomp; kfrac++)
                fracTot += fluC[i].fracMol[kfrac];
            for (int kfrac = 0; kfrac < ncomp; kfrac++)
                fluC[i].fracMol[kfrac] /= fracTot;
            pesoMolC = 0;
            for (int kfrac = 0; kfrac < ncomp; kfrac++)
                pesoMolC += fluC[i].masMol[kfrac] * celula[i].flui.fracMol[kfrac];
            if (kontaRenovaComp == arq.miniTabAtraso || arq.miniTabAtraso == 0) {
                fluC[i].atualizaPropCompStandard();
                if (fluC[i].dCalculatedBeta < 0. || fluC[i].dCalculatedBeta > 1.)
                    fluC[i].atualizaPropComp(celula[i].pres, celula[i].temp, -1, NULL, NULL, arq.pocinjec);
                else
                    fluC[i].atualizaPropComp(celula[i].pres, celula[i].temp, fluC[i].dCalculatedBeta,
                                             fluC[i].oCalculatedLiqComposition,
                                             fluC[i].oCalculatedVapComposition, arq.pocinjec);
                if (fluC[i].iIER != 0) {
                    int para;
                    para = 0;
                }
            }
        }
        betIV = celula[i].bet;
        rsV = celula[i].flui.RS(celula[i].pres, celula[i].temp);
        boV = celula[i].flui.BOFunc(celula[i].pres, celula[i].temp, rsV);
        baV = celula[i].flui.BAFunc(celula[i].pres, celula[i].temp);
        bswV = celula[i].flui.BSW * baV / (boV + baV * celula[i].flui.BSW - celula[i].flui.BSW * boV);
        rhoOVol = celula[i].flui.MasEspoleo(celula[i].pres, celula[i].temp);
        rhoWVol = celula[i].flui.MasEspAgua(celula[i].pres, celula[i].temp);
        titVol = (1 - bswV) * rhoOVol / ((1 - bswV) * rhoOVol + bswV * rhoWVol);
        celula[i].nMol = (fluC[i].MasEspLiq(celula[i].pres, celula[i].temp) * (1. - celula[i].alf) * (1. - betIV) * titVol +
                          fluC[i].MasEspGas(celula[i].pres, celula[i].temp) *
                              celula[i].alf) *
                         celula[i].duto.area * celula[i].dx / pesoMolC;
        fluC[i].Pmol = pesoMolC;
        double hol = 1. - celula[i].alf;
        double bet = celula[i].bet;

        double rs = celula[i].flui.RS(celula[i].pres, celula[i].temp);
        double bo = celula[i].flui.BOFunc(celula[i].pres, celula[i].temp, rs);
        double ba = celula[i].flui.BAFunc(celula[i].pres, celula[i].temp);
        double bsw = celula[i].flui.BSW * ba / (bo + ba * celula[i].flui.BSW - celula[i].flui.BSW * bo);

        double BSWini = celula[i].flui.BSW;
        double BSWF = celula[i].flui.BSW;
        double denagini = celula[i].flui.Denag;
        double denagF = celula[i].flui.Denag;
        double viscLini = 0 * 30 + 1 * celula[i].flui.VisOM(tL);
        double viscHini = 0 * 20 + 1 * celula[i].flui.VisOM(tH);
        double viscLF = viscLini;
        double viscHF = viscHini;

        if (celula[i].acsr.tipo == 2) {
            double rsF = celula[i].acsr.injl.FluidoPro.RS(celula[i].pres, celula[i].temp) * (6.29 / 35.31467);
            if (arq.nfluP > 0 || (*vg1dSP).chaverede != 0) {
                double bswaux = celula[i].acsr.injl.FluidoPro.BSW;
                double contrabsw = 1. - bswaux;
                if (contrabsw < (*vg1dSP).localtiny)
                    contrabsw = 0.9 * (*vg1dSP).localtiny;
                double rhoPSTF;
                if (contrabsw > (*vg1dSP).localtiny)
                    rhoPSTF = (1000 * 141.5 / (131.5 + celula[i].acsr.injl.FluidoPro.API)) + (bswaux / contrabsw) * 1000 * celula[i].acsr.injl.FluidoPro.Denag + celula[i].acsr.injl.FluidoPro.Deng * 1.225 * rsF;
                else
                    rhoPSTF = 1000 * celula[i].acsr.injl.FluidoPro.Denag;
                if (contrabsw > (*vg1dSP).localtiny) {
                    fonteP *= (1. / rhoPSTF);
                    fonteA *= ((bswaux / contrabsw) / rhoPSTF);
                } else {
                    fonteP = 0.;
                    fonteA *= (1 / (celula[i].acsr.injl.FluidoPro.BSW * 1000 * celula[i].acsr.injl.FluidoPro.Denag));
                }
                BSWF = celula[i].acsr.injl.FluidoPro.BSW;
                denagF = celula[i].acsr.injl.FluidoPro.Denag;
                viscLF = 0 * 30 + 1 * celula[i].acsr.injl.FluidoPro.VisOM(tL);
                viscHF = 0 * 20 + 1 * celula[i].acsr.injl.FluidoPro.VisOM(tH);
            }
        } else if (celula[i].acsr.tipo == 10) {
            double rsF = celula[i].acsr.injm.FluidoPro.RS(celula[i].pres, celula[i].temp) * (6.29 / 35.31467);
            if (arq.nfluP > 0 || (*vg1dSP).chaverede != 0) {
                double bswaux = celula[i].acsr.injm.FluidoPro.BSW;
                double contrabsw = 1. - bswaux;
                if (contrabsw < (*vg1dSP).localtiny)
                    contrabsw = 0.9 * (*vg1dSP).localtiny;
                double rhoPSTF;
                if (contrabsw > (*vg1dSP).localtiny)
                    rhoPSTF = (1000 * 141.5 / (131.5 + celula[i].acsr.injm.FluidoPro.API)) +
                              (bswaux / contrabsw) * 1000 *
                                  celula[i].acsr.injm.FluidoPro.Denag +
                              celula[i].acsr.injm.FluidoPro.Deng * 1.225 * rsF;
                else
                    rhoPSTF = 1000 * celula[i].acsr.injm.FluidoPro.Denag;
                if (contrabsw > (*vg1dSP).localtiny) {
                    fonteP *= (1. / rhoPSTF);
                    fonteA *= ((bswaux / contrabsw) / rhoPSTF);
                } else {
                    fonteP = 0.;
                    fonteA *= (1 / (celula[i].acsr.injm.FluidoPro.BSW * 1000 * celula[i].acsr.injm.FluidoPro.Denag));
                }
                BSWF = celula[i].acsr.injm.FluidoPro.BSW;
                denagF = celula[i].acsr.injl.FluidoPro.Denag;
                viscLF = 0 * 30 + 1 * celula[i].acsr.injm.FluidoPro.VisOM(tL);
                viscHF = 0 * 20 + 1 * celula[i].acsr.injm.FluidoPro.VisOM(tH);
            }
        } else if (celula[i].acsr.tipo == 3) {
            double rsF = celula[i].acsr.ipr.FluidoPro.RS(celula[i].pres, celula[i].temp) * (6.29 / 35.31467);
            if (arq.nfluP > 0 || (*vg1dSP).chaverede != 0) {
                double bswaux = celula[i].acsr.ipr.FluidoPro.BSW;
                double contrabsw = 1. - bswaux;
                if (contrabsw <= (*vg1dSP).localtiny)
                    contrabsw = 0.9 * (*vg1dSP).localtiny;
                double rhoPSTF;
                if (contrabsw > 0)
                    rhoPSTF = (1000 * 141.5 / (131.5 + celula[i].acsr.ipr.FluidoPro.API)) + (bswaux / contrabsw) * 1000 * celula[i].acsr.ipr.FluidoPro.Denag + celula[i].acsr.ipr.FluidoPro.Deng * 1.225 * rsF;
                else
                    rhoPSTF = 1000 * celula[i].acsr.ipr.FluidoPro.Denag;
                if (contrabsw > (*vg1dSP).localtiny) {
                    fonteP *= (1. / rhoPSTF);
                    fonteA *= ((bswaux / contrabsw) / rhoPSTF);
                } else {
                    fonteP = 0.;
                    fonteA *= (1 / (celula[i].acsr.ipr.FluidoPro.BSW * 1000 * celula[i].acsr.ipr.FluidoPro.Denag));
                }
                BSWF = celula[i].acsr.ipr.FluidoPro.BSW;
                denagF = celula[i].acsr.ipr.FluidoPro.Denag;
                viscLF = 0 * 30 + 1 * celula[i].acsr.ipr.FluidoPro.VisOM(tL);
                viscHF = 0 * 20 + 1 * celula[i].acsr.ipr.FluidoPro.VisOM(tH);
            }
        } else if (celula[i].acsr.tipo == 9) {

            double rsF = fluF.RS(celula[i].pres, celula[i].temp) * (6.29 / 35.31467);
            if (arq.nfluP > 0 || (*vg1dSP).chaverede != 0) {
                double bswaux = fluF.BSW;
                double contrabsw = 1. - bswaux;
                if (contrabsw <= (*vg1dSP).localtiny)
                    contrabsw = 0.9 * (*vg1dSP).localtiny;
                double rhoPSTF;
                if (contrabsw > (*vg1dSP).localtiny)
                    rhoPSTF = (1000 * 141.5 / (131.5 + fluF.API)) + (bswaux / contrabsw) * 1000 * fluF.Denag + fluF.Deng * 1.225 * rsF;
                else
                    rhoPSTF = 1000 * fluF.Denag;
                if (contrabsw > (*vg1dSP).localtiny) {
                    fonteP *= (1. / rhoPSTF);
                    fonteA *= ((bswaux / contrabsw) / rhoPSTF);
                } else {
                    fonteP = 0.;
                    fonteA *= (1 / (fluF.BSW * 1000 * fluF.Denag));
                }
                BSWF = fluF.BSW;
                denagF = fluF.Denag;
                viscLF = 0 * 30 + 1 * fluF.VisOM(tL);
                viscHF = 0 * 20 + 1 * fluF.VisOM(tH);
            }
        } else if (celula[i].acsr.tipo == 15) {
            double rsF = celula[i].acsr.radialPoro.flup.RS(celula[i].pres, celula[i].temp) * (6.29 / 35.31467);
            if (arq.nfluP > 0 || (*vg1dSP).chaverede != 0) {
                double bswaux = celula[i].acsr.radialPoro.BSW;
                double contrabsw = 1. - bswaux;
                if (contrabsw <= (*vg1dSP).localtiny)
                    contrabsw = 0.9 * (*vg1dSP).localtiny;
                double rhoPSTF;
                if (contrabsw > 0)
                    rhoPSTF = (1000 * 141.5 / (131.5 + celula[i].acsr.radialPoro.flup.API)) + (bswaux / contrabsw) * 1000 * celula[i].acsr.radialPoro.flup.Denag + celula[i].acsr.radialPoro.flup.Deng * 1.225 * rsF;
                else
                    rhoPSTF = 1000 * celula[i].acsr.radialPoro.flup.Denag;
                if (contrabsw > (*vg1dSP).localtiny) {
                    fonteP *= (1. / rhoPSTF);
                    fonteA *= ((bswaux / contrabsw) / rhoPSTF);
                } else {
                    fonteP = 0.;
                    fonteA *= (1 / (celula[i].acsr.radialPoro.BSW * 1000 * celula[i].acsr.radialPoro.flup.Denag));
                }
                BSWF = celula[i].acsr.radialPoro.BSW;
                denagF = celula[i].acsr.radialPoro.flup.Denag;
                viscLF = 0 * 30 + 1 * celula[i].acsr.radialPoro.flup.VisOM(tL);
                viscHF = 0 * 20 + 1 * celula[i].acsr.radialPoro.flup.VisOM(tH);
            }
        } else if (celula[i].acsr.tipo == 16) {
            double rsF = celula[i].acsr.poroso2D.dados.flup.RS(celula[i].pres, celula[i].temp) * (6.29 / 35.31467);
            if (arq.nfluP > 0 || (*vg1dSP).chaverede != 0) {
                double bswaux = celula[i].acsr.poroso2D.dados.transfer.BSW;
                double contrabsw = 1. - bswaux;
                if (contrabsw <= (*vg1dSP).localtiny)
                    contrabsw = 0.9 * (*vg1dSP).localtiny;
                double rhoPSTF;
                if (contrabsw > 0)
                    rhoPSTF = (1000 * 141.5 / (131.5 + celula[i].acsr.poroso2D.dados.flup.API)) + (bswaux / contrabsw) * 1000 * celula[i].acsr.poroso2D.dados.flup.Denag + celula[i].acsr.poroso2D.dados.flup.Deng * 1.225 * rsF;
                else
                    rhoPSTF = 1000 * celula[i].acsr.poroso2D.dados.flup.Denag;
                if (contrabsw > (*vg1dSP).localtiny) {
                    fonteP *= (1. / rhoPSTF);
                    fonteA *= ((bswaux / contrabsw) / rhoPSTF);
                } else {
                    fonteP = 0.;
                    fonteA *= (1 / (celula[i].acsr.poroso2D.dados.transfer.BSW * 1000 * celula[i].acsr.poroso2D.dados.flup.Denag));
                }
                BSWF = celula[i].acsr.poroso2D.dados.transfer.BSW;
                denagF = celula[i].acsr.poroso2D.dados.flup.Denag;
                viscLF = 0 * 30 + 1 * celula[i].acsr.poroso2D.dados.flup.VisOM(tL);
                viscHF = 0 * 20 + 1 * celula[i].acsr.poroso2D.dados.flup.VisOM(tH);
            }
        } else if ((fabs(fonteO) > (*vg1dSP).localtiny && celula[i].acsr.tipo != 2 && celula[i].acsr.tipo != 3 &&
                    celula[i].acsr.tipo != 9 && celula[i].acsr.tipo != 15 && celula[i].acsr.tipo != 16) ||
                   (fabs(fonteG) > (*vg1dSP).localtiny && celula[i].acsr.tipo != 1 && celula[i].acsr.tipo != 2 && celula[i].acsr.tipo != 3 && celula[i].acsr.tipo != 9 && celula[i].acsr.tipo != 15 && celula[i].acsr.tipo != 16)) {
            if (celula[i].acsr.tipo == 5 || celula[i].acsr.tipo == 8) {
                double rsF = celula[i].flui.RS(celula[i].pres, celula[i].temp) * (6.29 / 35.31467);
                if (arq.nfluP > 0 || (*vg1dSP).chaverede != 0) {
                    double rhoPSTF = (1 - celula[i].flui.BSW) * (1000 * 141.5 / (131.5 + celula[i].flui.API)) + celula[i].flui.BSW * 1000 * celula[i].flui.Denag;
                    fonteP *= ((1 - celula[i].flui.BSW) / rhoPSTF);
                    fonteA *= (celula[i].flui.BSW / rhoPSTF);
                    BSWF = celula[i].flui.BSW;
                    denagF = fluF.Denag;
                    viscLF = 0 * 30 + 1 * celula[i].flui.VisOM(tL);
                    viscHF = 0 * 20 + 1 * celula[i].flui.VisOM(tH);
                }
            } else if ((*celula[i].acsrL).tipo == 5 || (*celula[i].acsrL).tipo == 8) {
                if (i > 0) {
                    double rsF = celula[i - 1].flui.RS(celula[i].pres, celula[i].temp) * (6.29 / 35.31467);
                    if (arq.nfluP > 0 || (*vg1dSP).chaverede != 0) {
                        double rhoPSTF = (1 - celula[i - 1].flui.BSW) * (1000 * 141.5 / (131.5 + celula[i - 1].flui.API)) +
                                         celula[i - 1].flui.BSW * 1000 *
                                             celula[i - 1].flui.Denag;
                        fonteP *= ((1 - celula[i - 1].flui.BSW) / rhoPSTF);
                        fonteA *= (celula[i - 1].flui.BSW / rhoPSTF);
                        BSWF = celula[i - 1].flui.BSW;
                        denagF = celula[i - 1].flui.Denag;
                        viscLF = 0 * 30 + 1 * celula[i - 1].flui.VisOM(tL);
                        viscHF = 0 * 20 + 1 * celula[i - 1].flui.VisOM(tH);
                    }
                } else {
                    double rsF = celula[i].flui.RS(celula[i].pres, celula[i].temp) * (6.29 / 35.31467);
                    if (arq.nfluP > 0 || (*vg1dSP).chaverede != 0) {
                        double rhoPSTF = (1 - celula[i].flui.BSW) * (1000 * 141.5 / (131.5 + celula[i].flui.API)) +
                                         celula[i].flui.BSW * 1000 *
                                             celula[i].flui.Denag;
                        fonteP *= ((1 - celula[i].flui.BSW) / rhoPSTF);
                        fonteA *= (celula[i].flui.BSW / rhoPSTF);
                        BSWF = celula[i].flui.BSW;
                        denagF = celula[i - 1].flui.Denag;
                        viscLF = 0 * 30 + 1 * celula[i].flui.VisOM(tL);
                        viscHF = 0 * 20 + 1 * celula[i].flui.VisOM(tH);
                    }
                }
            }
        }
        if (titF > 1. - 1e-15) {
            fonteO = 0.;
            fonteP = 0.;
            fonteA = 0.;
        }

        double volpesFim = hol * (1 - bet) * (1 - bsw) / bo;
        double volaguaFim = hol * (1 - bet) * bsw / bo; // nao deveria ser dividido por Bo???????????
        double MultPe;
        double MultPd = 0.;
        double residuoP = 0.;
        double MultAe = 0.;
        double MultAd;
        double residuoA = 0.;
        if ((*vg1dSP).lixo5 < 1e-15 && arq.snaps != 1) {
            celula[i].VolPesaST = volpesFim;
            celula[i].VolAguaST = volaguaFim;
        }
        if (arq.nfluP > 1 || (*vg1dSP).chaverede != 0) {
            MultPe = 0.;
            if (celula[i].flui.dStockTankVaporMassFraction < 1. - 1e-15)
                MultPe = celula[i].QL * (1 - betI0) * (1 - bsw0) / bo0;
            MultPd = 0.;
            if (celula[i + 1].flui.dStockTankVaporMassFraction < 1. - 1e-15)
                MultPd = celula[i + 1].QL * (1 - betI1) * (1 - bsw1) / bo1;
            residuoP = (volpesFim - celula[i].VolPesaST) * at / dt + (MultPd - MultPe) / dx - fonteP / dx;
            MultAe = 0.;
            if (celula[i].flui.dStockTankVaporMassFraction < 1. - 1e-15)
                MultAe = celula[i].QL * (1 - betI0) * bsw0 / bo0;
            MultAd = 0.;
            if (celula[i + 1].flui.dStockTankVaporMassFraction < 1. - 1e-15)
                MultAd = celula[i + 1].QL * (1 - betI1) * bsw1 / bo1;
            residuoA = (volaguaFim - celula[i].VolAguaST) * at / dt + (MultAd - MultAe) / dx - fonteA / dx;
        }
        if ((arq.nfluP > 1) || (*vg1dSP).chaverede != 0) {
            if (((volpesFim > 1e-3) && (fonteP > 0 || (MultPd < 0 || MultPe > 0))) && arq.nfluP > 0) {
                VISCL[i] = (dt * (viscLF * fonteP / dx + 1. * viscLini * residuoP - (viscL1 * MultPd - viscL0 * MultPe) / dx) + viscLini * celula[i].VolPesaST * at) / (volpesFim * at - 0. * residuoP * dt);
                VISCH[i] = (dt * (viscHF * fonteP / dx + 1. * viscHini * residuoP - (viscH1 * MultPd - viscH0 * MultPe) / dx) + viscHini * celula[i].VolPesaST * at) / (volpesFim * at - 0. * residuoP * dt);
            } else {
                VISCL[i] = viscLini;
                VISCH[i] = viscHini;
            }
        }
        if (arq.nfluP > 1 || (*vg1dSP).chaverede != 0) {
            if ((
                    (volaguaFim + volpesFim) > 1e-3) &&
                ((fonteP > 0 || fonteA > 0) ||
                 ((MultPd < 0 || MultPe > 0) || (MultAd < 0 || MultAe > 0)))) {
                BSW[i] = (dt * (BSWF * (fonteA + fonteP) / dx + 1. * BSWini * (residuoA + residuoP) - (BSW1 * (MultAd + MultPd) - BSW0 * (MultAe + MultPe)) / dx) + BSWini * (celula[i].VolAguaST + celula[i].VolPesaST) * at) / ((volaguaFim + volpesFim) * at -
                                                                                                                                                                                                                                  0. * dt * (residuoA + residuoP));
                denag[i] = (dt * (denagF * (fonteA) / dx + 1. * denagini * (residuoA) - (denag1 * (MultAd)-denag0 * (MultAe)) / dx) + denagini * (celula[i].VolAguaST) * at) / ((volaguaFim)*at -
                                                                                                                                                                                0. * dt * (residuoA));
                if (BSW[i] < 0.)
                    BSW[i] = 0.;
                else if (BSW[i] > 1.)
                    BSW[i] = 1.;
                else if (isnan(BSW[i]))
                    BSW[i] = BSWini;
                if (denag[i] < 1.)
                    denag[i] = 1.;
                else if (isnan(denag[i]))
                    denag[i] = denagini;
            } else
                BSW[i] = BSWini;
            denag[i] = denagini;
        }
        celula[i].VolPesaST = volpesFim - 0. * residuoP * dt / at;
        celula[i].VolAguaST = volaguaFim - 0. * residuoA * dt / at;
    }
    for (int i = imin; i <= ncel - 1; i++) {
        if (arq.nfluP > 0 || (*vg1dSP).chaverede != 0) {
            celula[i].flui = fluC[i];
            if (arq.nfluP > 1 || (*vg1dSP).chaverede != 0) {
                celula[i].flui.BSW = BSW[i];
                celula[i].flui.Denag = denag[i];
                celula[i].flui.LVisL = VISCL[i];
                celula[i].flui.LVisH = VISCH[i];
                celula[i].flui.TempL = tL;
                celula[i].flui.TempH = tH;
            }
        }
    }

    if ((*vg1dSP).chaverede == 0 || noextremo == 1 || celula[ncel].Mliqini > -(*vg1dSP).localtiny) {
        if (arq.nfluP > 0 || (*vg1dSP).chaverede != 0) {
            celula[ncel].flui = celula[ncel - 1].flui;
            if (arq.nfluP > 1 || (*vg1dSP).chaverede != 0) {
                celula[ncel].flui.BSW = BSW[ncel - 1];
                celula[ncel].flui.Denag = denag[ncel - 1];
                celula[ncel].flui.LVisL = VISCL[ncel - 1];
                celula[ncel].flui.LVisH = VISCH[ncel - 1];
                celula[ncel].flui.TempL = tL;
                celula[ncel].flui.TempH = tH;
            }
        }
    } else {
        celula[ncel].flui = fluiRev;
        celula[ncel].flui.BSW = fluiRev.BSW;
        celula[ncel].flui.Denag = fluiRev.Denag;
        celula[ncel].flui.LVisL = fluiRev.LVisL;
        celula[ncel].flui.LVisH = fluiRev.LVisH;
        celula[ncel].flui.TempL = fluiRev.TempL;
        celula[ncel].flui.TempH = fluiRev.TempH;
    }
    delete[] fluC;
    if(arq.miniTabAtraso > 0)(*vg1dSP).modoTransiente = 1;
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
/// Kept a member of SProd, and it must stay one: both solvers are constructed from
/// the whole SProd object (*this), which a free function taking a state struct
/// does not have. SolveTrans reaches it through a callback once it moves, so the
/// hydrate phase still runs first, where T127 requires it.
void SProd::solveHydrateEnvelopes() {
    if (arq.calculaEnvelope == 1 && (*vg1dSP).lixo5 <= arq.tfinal) { //*vg1dSP).lixo5>0 && //chris - Hidratos
        FA_Hidrato solverHidrato(*this);
        solverHidrato.solverHidrato();
    }

    if (arq.lingas > 0 && arq.calculaEnvelope == 1 && (*vg1dSP).lixo5 <= arq.tfinal) { //*vg1dSP).lixo5>0 && //chris - Hidratos
        FA_Hidrato_Servico solverHidratoG(*this);
        solverHidratoG.solverHidratoG();
    }
}

void SProd::SolveTrans(double titRev, double alfRev, double betRev, int nrede, ProFlu fluiRev) {
    sisprod::transient::solveTransientStep(transientSolveStateOf(*this), titRev, alfRev, betRev, nrede, fluiRev);
}


namespace {

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
        .inputData = system.arq,
        .globals = system.vg1dSP,
        .branchIndex = system.indTramo,
        .printPassCount = system.kimpT,
        .production = {system.MatTrendP, system.ntrend, system.ntrendB},
        .service = {system.MatTrendG, system.ntrendg, system.ntrendgB},
        .productionCrossSection = {system.MatTrendTransP, system.ntrendtrans, system.ntrendtransB},
        .serviceCrossSection = {system.MatTrendTransG, system.ntrendtransg, system.ntrendtransgB}};
}

} // namespace

// The eight trend writers moved to SisProdTrendOutput.cpp. Four of them are
// called from 28 sites in Num4Main.cpp, so the public signatures stay exactly
// as they were (FR-032) and these forward. Nothing in SisProd.h changed.
void SProd::ImprimeTrendPCab(int i, int nrede) {
    trendoutput::writeProductionTrendHeader(trendStateOf(*this), i, nrede);
}
void SProd::ImprimeTrendP(int i, int nrede) {
    trendoutput::writeProductionTrendRows(trendStateOf(*this), i, nrede);
}
void SProd::ImprimeTrendGCab(int i, int nrede) {
    trendoutput::writeServiceTrendHeader(trendStateOf(*this), i, nrede);
}
void SProd::ImprimeTrendG(int i, int nrede) {
    trendoutput::writeServiceTrendRows(trendStateOf(*this), i, nrede);
}
void SProd::ImprimeTrendTransPCab(int i) {
    trendoutput::writeProductionCrossSectionTrendHeader(trendStateOf(*this), i);
}
void SProd::ImprimeTrendTransP(int i) {
    trendoutput::writeProductionCrossSectionTrendRows(trendStateOf(*this), i);
}
void SProd::ImprimeTrendTransGCab(int i) {
    trendoutput::writeServiceCrossSectionTrendHeader(trendStateOf(*this), i);
}
void SProd::ImprimeTrendTransG(int i) {
    trendoutput::writeServiceCrossSectionTrendRows(trendStateOf(*this), i);
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
        .inletColumnFraction = betaE,
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
    sisprod::steady::serviceLineHydrostatic(steadyStateOf(*this));
}

double SProd::buscaTramoSecVazPerm(double pPartida, int indPartida) {
    return sisprod::steady::searchSecondaryBranchFlowRate(searchStateOf(*this), pPartida, indPartida);
}
