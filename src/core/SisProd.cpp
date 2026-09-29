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
// so that SProd can name them as friends; these keep every call site as it was.
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
using sisprod::kAtmospherePerKgfPerCm2;
using sisprod::kBarrelPerCubicMetre;
using sisprod::kCubicFootPerCubicMetre;
using sisprod::kGravity;
using sisprod::kKgfPerCm2PerPascal;
using sisprod::kPascalPerKgfPerCm2Coarse;
using sisprod::kPascalPerKgfPerCm2PvtSim;
using sisprod::kPsiPerAtmosphere;
using sisprod::kPsiPerPascal;
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

/// Allocates one set of trend matrices: for each of the count trends, length[i]
/// rows of rowWidth(i) values whose first sentinelCount(i) hold the -10000
/// sentinel, and the reset timers and sample counts, zeroed.
template <typename RowWidth, typename SentinelCount>
void allocateTrendSet(int count, const int *length, RowWidth rowWidth, SentinelCount sentinelCount,
                      double ***&matrices, double *&resetTimers, int *&counts, int *&bufferedCounts) {
    resetTimers = new double[count];
    counts = new int[count];
    bufferedCounts = new int[count];
    matrices = new double **[count];
    for (int i = 0; i < count; i++) {
        matrices[i] = new double *[length[i]];
        for (int j = 0; j < length[i]; j++) {
            matrices[i][j] = new double[rowWidth(i)];
            for (int k = 0; k < sentinelCount(i); k++)
                matrices[i][j][k] = -10000.;
        }
        resetTimers[i] = 0;
        counts[i] = 0;
        bufferedCounts[i] = 0;
    }
}

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

namespace {

/// Calls apply on each fluid the cell's source carries: that of a liquid,
/// inflow-performance or multiple source, or, for a porous reservoir, its own
/// and those of its cells. A gas injection's fluid is not visited. An apply that
/// also takes a pressure and a temperature gets the cell's for the source's own
/// fluid and each reservoir cell's for that cell's.
template <typename Apply>
void forEachSourceFluid(Cel &cell, Apply &&apply) {
    auto visit = [&](ProFlu &fluid, double pressure, double temperature) {
        if constexpr (std::invocable<Apply &, ProFlu &>)
            apply(fluid);
        else
            apply(fluid, pressure, temperature);
    };
    acessorio &source = cell.acsr;
    if (source.tipo == kAccessoryLiquidInjection) {
        visit(source.injl.FluidoPro, cell.pres, cell.temp);
    } else if (source.tipo == kAccessoryInflowPerformance) {
        visit(source.ipr.FluidoPro, cell.pres, cell.temp);
    } else if (source.tipo == kAccessoryMultipleSource) {
        visit(source.injm.FluidoPro, cell.pres, cell.temp);
    } else if (source.tipo == kAccessoryRadialPorous) {
        visit(source.radialPoro.flup, cell.pres, cell.temp);
        for (int k = 0; k < source.radialPoro.ncel; k++)
            visit(source.radialPoro.celula[k].flup, source.radialPoro.celula[k].Pcamada, source.radialPoro.tRes);
    } else if (source.tipo == kAccessoryPorous2D) {
        visit(source.poroso2D.dados.flup, cell.pres, cell.temp);
        for (int k = 0; k < source.poroso2D.dados.transfer.ncel; k++)
            visit(source.poroso2D.dados.transfer.celula[k].flup, source.poroso2D.dados.transfer.celula[k].Pcamada,
                  source.poroso2D.dados.transfer.tRes);
        for (int k = 0; k < source.poroso2D.malha.nele; k++)
            visit(source.poroso2D.malha.mlh2d[k].flup, source.poroso2D.malha.mlh2d[k].cel2D.presC,
                  source.poroso2D.malha.mlh2d[k].tRes);
    }
}

}  // namespace

/// Points every cell fluid, and every source fluid it carries, at the bubble-point
/// tables read from the PVTSim file, and switches them to saturation model 4.
void SProd::assignPvtSimBubbleTablesToCells() {
    auto assignTables = [&](ProFlu &fluid) {
        fluid.PBPVTSim = PBPVTSim;
        fluid.TBPVTSim = TBPVTSim;
        fluid.corrSat = 4;
    };
    for (int i = 0; i <= ncel; i++) {
        assignTables(celula[i].flui);
        forEachSourceFluid(celula[i], assignTables);
    }
}

namespace {

/// Reads one row of a PVTSim table file: skips words up to `name`, then
/// converts the numbers on the rest of that line into row[0] to row[last].
template <typename Row, typename Convert>
void readPvtSimRow(ifstream &file, string &word, const char *name, Row &row, int last, Convert convert) {
    while (word != name)
        file >> word;
    char line[4000];
    file.get(line, 4000);
    char *token = strtok(line, " ,()=");
    row[0] = convert(atof(token));
    for (int k = 1; k <= last; k++) {
        token = strtok(NULL, " ,");
        row[k] = convert(atof(token));
    }
}

/// Writes a table to the file at `path`, replacing it.
void writeTable(const string &path, const FullMtx<double> &table) {
    ofstream file(path.c_str(), ios_base::out);
    file << table;
}

}  // namespace

/// Reads the bubble-point curve from the PVTSim file, points the cell fluids at it
/// and writes perfilBolha; with tabRSPB on, also reads the solution gas-oil ratio
/// table and writes perfilRSLivia.
void SProd::loadPvtSimSaturationTables() {
    LerPB = 1;
    int ndiv = arq.tabent.npont - 1;
    PBPVTSim = new double[ndiv + 1];
    TBPVTSim = new double[ndiv + 1];
    vector<double> PresPVTSim(ndiv + 1);

    string impfile;
    impfile = arq.pvtsimarq;
    string dadosMR = impfile;
    ifstream lendoPVTSim(dadosMR.c_str(), ios_base::in);
    string chave;
    char line[4000];
    lendoPVTSim.get(line, 4000);
    lendoPVTSim >> chave;
    int lacoleitura = ndiv;
    readPvtSimRow(lendoPVTSim, chave, "PRESSURE", PresPVTSim, lacoleitura,
                  [](double value) { return value / kPascalPerKgfPerCm2PvtSim; });
    readPvtSimRow(lendoPVTSim, chave, "BUBBLEPRESSURES", PBPVTSim, lacoleitura,
                  [](double value) { return value * kPsiPerPascal; });
    readPvtSimRow(lendoPVTSim, chave, "BUBBLETEMPERATURES", TBPVTSim, lacoleitura, [](double value) { return value; });
    assignPvtSimBubbleTablesToCells();

    FullMtx<double> BolhaTemp(ndiv + 2, 2);
    for (int i = 0; i <= ndiv; i++) {
        BolhaTemp[i][0] = TBPVTSim[i];
        BolhaTemp[i][1] = PBPVTSim[i];
    }
    writeTable(pathPrefixoArqSaida + "perfilBolha", BolhaTemp);

    if (arq.tabRSPB == 1) {
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
                    double TFa = 1.8 * TBPVTSim[j - 1] + 32;

                    constexpr double A0 = 6542.69213 * 1e-11;
                    constexpr double A1 = 2.60464618;
                    constexpr double A2 = 1.21334544;
                    constexpr double A3 = 0.16464125;
                    constexpr double A4 = 1.19382967;
                    constexpr double B0 = 87 * 1e-8;
                    constexpr double B1 = 4.77390016;
                    constexpr double B2 = 1.77267703;
                    constexpr double B3 = -0.73072038;
                    constexpr double B4 = 1.58093857;
                    constexpr double C0 = 3226.09 * 1e-6;
                    constexpr double C1 = 0.09281075;
                    constexpr double C2 = 0.13633665;
                    constexpr double C3 = 0.09634381;
                    constexpr double C4 = 0.53238728;
                    constexpr double D0 = 1.00544053;
                    constexpr double D1 = -0.00134177;
                    constexpr double D2 = 0.51839397;

                    double yco2 = celula[0].flui.yco2;
                    double Deng = celula[0].flui.Deng;
                    double API = celula[0].flui.API;
                    double multCor = (D0 + D1 * yco2 * pow(TFa, D2));
                    double pbtemp = PBPVTSim[j - 1];

                    double pcor = (RSLivia[i][0] * kAtmospherePerKgfPerCm2) * kPsiPerAtmosphere * multCor;
                    double pr = pcor / pbtemp;

                    double a1 = A0 * pow(Deng, A1) * pow(API, A2) * pow(TFa, A3) * pow(pbtemp, A4);
                    double a2 = B0 * pow(Deng, B1) * pow(API, B2) * pow(TFa, B3) * pow(pbtemp, B4);
                    double a3 = C0 * pow(Deng, C1) * pow(API, C2) * pow(TFa, C3) * pow(pbtemp, C4);

                    double Rsr = a1 * pow(pr, a2) + (1. - a1) * pow(pr, a3);
                    double rstemp = 0.;
                    if (rstemp < 0.)
                        rstemp = 0.;
                    else if ((RSLivia[i][0] * kAtmospherePerKgfPerCm2) * kPsiPerAtmosphere > pbtemp)
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
            for (int j = 1; j <= ndiv + 1; j++)
                RSTemp[i][j] = RSLivia[i][j] * kBarrelPerCubicMetre / kCubicFootPerCubicMetre;
        }
        writeTable(pathPrefixoArqSaida + "perfilRSLivia", RSTemp);

        auto pointAtRatioTable = [&](ProFlu &fluid) {
            fluid.TabRSLivia = RSLivia;
            fluid.tabRSPB = 1;
        };
        for (int i = 0; i <= ncel; i++) {
            pointAtRatioTable(celula[i].flui);
            forEachSourceFluid(celula[i], pointAtRatioTable);
        }
    }
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
            RSTemp[i][j] = RSLivia[i][j] * kBarrelPerCubicMetre / kCubicFootPerCubicMetre;
            ttestepb += dtteste;
        }
        pteste += dpteste;
    }

    writeTable(pathPrefixoArqSaida + "perfilBolha", PBTemp);
    writeTable(pathPrefixoArqSaida + "perfilRSLivia", RSTemp);

    auto pointAtTables = [&](ProFlu &fluid) {
        fluid.PBPVTSim = PBPVTSim;
        fluid.TBPVTSim = TBPVTSim;
        fluid.TabRSLivia = RSLivia;
        fluid.tabRSPB = 1;
    };
    for (int i = 0; i <= ncel; i++) {
        pointAtTables(celula[i].flui);
        forEachSourceFluid(celula[i], pointAtTables);
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
    if (arq.modelcp > 0)
        arq.geraTabCp();
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
        celula[iposp].acsr.tipo = kAccessoryGasInjection;
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

namespace {

/// Gives `cell` a multiple, a liquid and a gas source, all with zero flow,
/// carrying the fluids of `inlet` (and, for the multiple source, its
/// temperature).
void attachEmptySources(Cel &cell, Cel &inlet) {
    cell.acsr.injm = InjMult(0, 0, 0, inlet.temp, inlet.flui, inlet.fluicol);
    cell.acsr.injl = InjLiq(0, 0, 0, inlet.flui, inlet.fluicol);
    cell.acsr.injg = InjGas(0, 0, inlet.flui, inlet.fluicol);
}

}  // namespace

/// Gives the inlet the sources its boundary condition needs (and the second cell,
/// under a blockage), then places the accessories -- pumps, volumetric pumps,
/// pressure-drop requirements, heat sources, the master valve and the other
/// valves -- and sets up the outlet pressure, the surface and injection chokes and
/// the pigs.
void SProd::configureInletSourcesAndAccessories(int nfontes) {
    if (celula[0].acsr.tipo != kAccessoryNone && arq.ConContEntrada == 1) {
        // RN-302: The inlet has both a pressure boundary condition and a mass source. Report a warning.
        logger.log(LOGGER_AVISO, LOG_ERR_PARSE_BUSINESS_RULE_VALIDATION,
                   "Foi escolhida uma condicao de contorno de pressao no inicio da tubulacao e foi colocada uma fonte de massa neste inicio, isto pode causar inconsistencias",
                   "", "");
    }
    if ((arq.perm == 1 || (*vg1dSP).chaverede > 0) && celula[0].acsr.tipo == kAccessoryNone) {
        if ((arq.ConContEntrada == 0 && nfontes == 0 && (*vg1dSP).chaverede == 0)) {
            NumError("O simulador pede calculo de permanente com condicao de vazao na entrada, mas sem nenhuma fonte na entrada");
        }
        if ((arq.ConContEntrada == 0 && nfontes == 0 && (*vg1dSP).chaverede > 0) || arq.CCPres.tit[0] < (1 - (*vg1dSP).localtiny)) {
            if (arq.tipoFluido == 1)
                celula[0].acsr.tipo = kAccessoryGasInjection;
            else if (arq.tipoFluido == 0)
                celula[0].acsr.tipo = kAccessoryLiquidInjection;
            else
                celula[0].acsr.tipo = kAccessoryMultipleSource;
            attachEmptySources(celula[0], celula[0]);
        } else if (nfontes == 0 || arq.ConContEntrada == 1) {
            attachEmptySources(celula[0], celula[0]);
            celula[0].acsr.injg.seco = 0;
            if (arq.tipoFluido == 1)
                celula[0].acsr.tipo = kAccessoryGasInjection;
            else if (arq.tipoFluido == 0)
                celula[0].acsr.tipo = kAccessoryLiquidInjection;
            else
                celula[0].acsr.tipo = kAccessoryMultipleSource;
        }
    }
    if (bloq == 1 && celula[1].acsr.tipo == kAccessoryNone) {
        if ((*vg1dSP).fluidoRede == 1)
            celula[1].acsr.tipo = kAccessoryLiquidInjection;
        else if ((*vg1dSP).fluidoRede == 0)
            celula[1].acsr.tipo = kAccessoryGasInjection;
        else
            celula[1].acsr.tipo = kAccessoryMultipleSource;
        attachEmptySources(celula[1], celula[0]);
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
        while (celula[verifica].acsr.tipo != kAccessoryNone) {
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
            double throatArea = M_PI * diaG * diaG / 4.;
            chokeVGL[i] = ChokeGas(arq.flug, throatArea,
                                   arq.valvgl[i].diaexter, arq.valvgl[i].cd, presEstag, presGarg,
                                   tempEstag, arq.valvgl[i].frec, arq.valvgl[i].tipo,
                                   throatArea / arq.valvgl[i].razarea,
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
    if (celula[0].acsr.tipo != kAccessoryGasInjection && celula[0].acsr.tipo != kAccessoryLiquidInjection && celula[0].acsr.tipo != kAccessoryInflowPerformance && celula[0].acsr.tipo != kAccessoryMultipleSource && celula[0].acsr.tipo != kAccessoryRadialPorous && celula[0].acsr.tipo != kAccessoryPorous2D) {
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

    if (arq.descarga == 1 && celula[0].acsr.tipo != kAccessoryInflowPerformance) {
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

    if ((celula[ncel].acsr.tipo != kAccessoryNone || celula[ncel - 1].acsr.tipo != kAccessoryNone) && arq.transiente == 1) {
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
        ntabDin = 1;
        tabDin.push_back(tabelaDinamica());
        for (int i = 1; i < ncel; i++) {
            if ((celula[i].acsr.tipo == kAccessoryGasInjection) ||
                (celula[i].acsr.tipo == kAccessoryLiquidInjection) ||
                celula[i].acsr.tipo == kAccessoryInflowPerformance ||
                (celula[i].acsr.tipo == kAccessoryMultipleSource) || celula[i].acsr.tipo == kAccessoryLeak || celula[i].acsr.tipo == kAccessoryRadialPorous || celula[i].acsr.tipo == kAccessoryPorous2D) {
                tabelaDinamica temp;
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

            readPvtSimRow(lendoPVTSim, chave, "PRESSURE", presPVTSim, ndiv,
                          [](double value) { return value * kKgfPerCm2PerPascal; });
            readPvtSimRow(lendoPVTSim, chave, "TEMPERATURE", tempPVTSim, ndiv, [](double value) { return value; });

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
            string tmp = pathPrefixoArqSaida + "perfilLatente.dat";
            writeTable(tmp, HLatTemp);
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
            if (celula[i].acsr.tipo == kAccessoryLiquidInjection) {
                celula[i].acsr.injl.FluidoPro.rDgD = 1.;
                celula[i].acsr.injl.FluidoPro.rDgL = 1.;
            } else if (celula[i].acsr.tipo == kAccessoryInflowPerformance) {
                celula[i].acsr.ipr.FluidoPro.rDgD = 1.;
                celula[i].acsr.ipr.FluidoPro.rDgL = 1.;
            } else if (celula[i].acsr.tipo == kAccessoryMultipleSource) {
                celula[i].acsr.injm.FluidoPro.rDgD = 1.;
                celula[i].acsr.injm.FluidoPro.rDgL = 1.;
            }
            if (celula[i].acsr.tipo == kAccessoryRadialPorous) {
                celula[i].acsr.radialPoro.flup.rDgD = 1.;
                celula[i].acsr.radialPoro.flup.rDgD = 1.;
                for (int iRP = 0; iRP < celula[i].acsr.radialPoro.ncel; iRP++) {
                    celula[i].acsr.radialPoro.celula[iRP].flup.rDgD = 1.;
                    celula[i].acsr.radialPoro.celula[iRP].flup.rDgD = 1.;
                }
            }
            if (celula[i].acsr.tipo == kAccessoryPorous2D) {
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
            forEachSourceFluid(celula[i], [](ProFlu &fluid, double pressure, double temperature) {
                fluid.razDegD(pressure, temperature);
                fluid.rzDegL(pressure, temperature);
            });
        }
    }

    if (celula[0].acsr.tipo == kAccessoryGasInjection) {
        celula[0].flui = celula[0].acsr.injg.FluidoPro;
    } else if (celula[0].acsr.tipo == kAccessoryLiquidInjection) {
        celula[0].flui = celula[0].acsr.injl.FluidoPro;
    } else if (celula[0].acsr.tipo == kAccessoryInflowPerformance) {
        celula[0].flui = celula[0].acsr.ipr.FluidoPro;
    } else if (celula[0].acsr.tipo == kAccessoryMultipleSource) {
        celula[0].flui = celula[0].acsr.injm.FluidoPro;
    } else if (celula[0].acsr.tipo == kAccessoryRadialPorous) {
        celula[0].flui = celula[0].acsr.radialPoro.flup;
    } else if (celula[0].acsr.tipo == kAccessoryPorous2D) {
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
            allocateTrendSet(arq.ntendp, TrendLengthP, [&](int i) { return arq.nvartrendp[i] + 2; },
                             [&](int i) { return arq.nvartrendp[i] + 1; }, MatTrendP, resettrend, ntrend, ntrendB);
        }
        if (arq.ntendg > 0) {
            TrendLengthG = new int[arq.ntendg];
            for (int i = 0; i < arq.ntendg; i++)
                TrendLengthG[i] = 1 + 1 + ceil((*vg1dSP).TmaxR / arq.trendg[i].dt);
            allocateTrendSet(arq.ntendg, TrendLengthG, [&](int i) { return arq.nvartrendg[i] + 2; },
                             [&](int i) { return arq.nvartrendg[i] + 1; }, MatTrendG, resettrendg, ntrendg, ntrendgB);
        }
        if (arq.ntendtransp > 0) {
            TrendLengthTransP = new int[arq.ntendtransp];
            for (int i = 0; i < arq.ntendtransp; i++)
                TrendLengthTransP[i] = 1 + 1 + ceil((*vg1dSP).TmaxR / arq.trendtransp[i].dt);
            allocateTrendSet(arq.ntendtransp, TrendLengthTransP, [](int) { return 2; }, [](int) { return 2; },
                             MatTrendTransP, resettrendtrans, ntrendtrans, ntrendtransB);
        }
        if (arq.ntendtransg > 0 && arq.lingas > 0) {
            TrendLengthTransG = new int[arq.ntendtransg];
            for (int i = 0; i < arq.ntendtransg; i++)
                TrendLengthTransG[i] = 1 + 1 + ceil((*vg1dSP).TmaxR / arq.trendtransg[i].dt);
        }
        if (arq.ntendtransg > 0) {
            allocateTrendSet(arq.ntendtransg, TrendLengthTransG, [](int) { return 2; }, [](int) { return 2; },
                             MatTrendTransG, resettrendtransg, ntrendtransg, ntrendtransgB);
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

void sisprod::gaslift::GasLiftTemperatureUpdater::computeGasTemperature(
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
        .surfaceChokeOpen = system.aberto,
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

/// The composition module's view of SProd: seventeen members by reference,
/// generated from SisProdComposition.h by composition-move.py so the field
/// order is the header's, as designated initialisers require.
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
/// plus the 52 members only the solve touches. Initialiser order is the
/// header's declaration order, which C++20 requires.
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
        .networkCoupled = system.verificaAcop,
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

}  // namespace sisprod::adapters

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
  	//alteracao hidrato 3
  	if (arq.calculaEnvelope==1 && arq.tipoHmodel==3 && celula[ind].flui.BSW>1e-14 && (*vg1dSP).lixo5>0.01) { //alteracao Hidratos
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

/*** alteracao4 ***/

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
