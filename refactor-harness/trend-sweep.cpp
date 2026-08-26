/*
 * trend-sweep.cpp -- dedicated harness for the eight trend writers.
 *
 * Why this exists
 * ---------------
 * No model in the demo corpus declares a cross-section trend, so
 * writeProductionCrossSectionTrend{Header,Rows} and their service counterparts
 * are never executed by L2 or L3 (evidencia/trend-diff.md, D4-02). They are also
 * the writers that carry four of the five preserved anomalies, and T057
 * restructures all eight into a Template Method. Restructuring the least covered
 * code with no test covering it is the largest risk in the stage.
 *
 * The writers take a TrendState and write files. That makes them callable
 * without an SProd: this driver builds synthetic input decks and buffers, runs
 * every writer across a sweep of configurations, and leaves the produced files
 * in a directory. Two builds of the module -- before and after T057 -- must
 * produce byte-identical directories.
 *
 * The sweep varies every point of variation catalogued in T052: branch index,
 * AP sequence, print pass, output language, injector guard, caption flags, the
 * -9999 sentinel, and window sizes that differ from the buffer size.
 *
 * Usage: trend-sweep <output-directory>
 */
#include "Leitura.h"
#include "SisProd.h"
#include "estruturas.h"
#include "SisProdTrendOutput.h"
#include "variaveisGlobais1D.h"

#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <sys/stat.h>
#include <unistd.h>

// Num4Main.cpp holds main() and cannot be linked in, so the file-scope globals
// it defines are redefined here. Only arqRelatorioPerfis and
// pathPrefixoArqSaida are read by the trend writers; the rest exist purely to
// satisfy the linker for the objects that come along for the ride.
ofstream arqRelatorioPerfis;
string pathPrefixoArqSaida("");

Logger logger("");
int nthrdMatriz = 1;
string versao("harness");
string pathArqEntrada("");
string pathArqExtEntrada("");
string arqSaidaSnapShot("");
string diretorioSaida("");
string nomeRedePrincipal;
int logRede = 0;
int redeLeitura = 0;
int diaIni, horaIni, minutoIni, segundoIni;
time_t nowGlobIni, nowGlobFim;
tm *ltmGlobIni;
tm *ltmGlobFim;
detTempo tempVF;
detProp prop;
detMapProp mapprop;
detCI CI;
detCC CC;
SProd *ptrSistemaProducao;
tipoSimulacao_t tipoSimulacao = tipoSimulacao_t::transiente;
const char *saidaTexto[16] = {"", "", "", "", "", "", "", "",
                              "", "", "", "", "", "", "", ""};
const char *saidaSubTexto[16] = {"", "", "", "", "", "", "", "",
                                 "", "", "", "", "", "", "", ""};

namespace {

// Deterministic, spread over the flag space: no run of identical values, and
// every flag is exercised in both states across the three patterns.
int pattern(int seed) { return ((seed * 7 + 3) % 5 == 0) ? 1 : 0; }
int allOn(int) { return 1; }

double sample(int a, int b, int c) {
    // Values that exercise precision(19): irrational-looking, mixed magnitude,
    // and a negative. Nothing here is compared numerically -- only the text.
    const double base = 1.0 + a * 0.30000000000000004 + b * 7.125e-3 - c * 2.5e4;
    return base * 1.0000000000000002 + (a - b) * 3.14159265358979312;
}

double ***makeBuffer(int series, int rows, int columns, int salt) {
    double ***buffer = new double **[series];
    for (int s = 0; s < series; s++) {
        buffer[s] = new double *[rows];
        for (int r = 0; r < rows; r++) {
            buffer[s][r] = new double[columns];
            for (int c = 0; c < columns; c++)
                buffer[s][r][c] = sample(s + salt, r, c);
        }
    }
    // The row loop stops at the first value <= -9999 in column 0. Plant one, so
    // the sentinel path is exercised and not only the run-to-the-end path.
    if (rows > 3)
        buffer[0][rows - 2][0] = -12345.5;
    return buffer;
}

void fillProduction(detTRENDP &trend, int (*flag)(int), int seed, int comp) {
    trend.dt = 1.5;
    trend.posic = 2;
    trend.comp = comp;
    trend.rotulo = "rotulo-P";
    trend.pres = flag(seed++);
    trend.temp = flag(seed++);
    trend.hol = flag(seed++);
    trend.FVH = flag(seed++);
    trend.bet = flag(seed++);
    trend.ugs = flag(seed++);
    trend.uls = flag(seed++);
    trend.ug = flag(seed++);
    trend.ul = flag(seed++);
    trend.arra = flag(seed++);
    trend.yco2 = flag(seed++);
    trend.viscl = flag(seed++);
    trend.viscg = flag(seed++);
    trend.rhog = flag(seed++);
    trend.rhol = flag(seed++);
    trend.rhoMix = flag(seed++);
    trend.masg = flag(seed++);
    trend.masl = flag(seed++);
    trend.c0 = flag(seed++);
    trend.ud = flag(seed++);
    trend.RGO = flag(seed++);
    trend.deng = flag(seed++);
    trend.calor = flag(seed++);
    trend.masstrans = flag(seed++);
    trend.qlst = flag(seed++);
    trend.qlwst = flag(seed++);
    trend.qlstTot = flag(seed++);
    trend.qgst = flag(seed++);
    trend.api = flag(seed++);
    trend.bsw = flag(seed++);
    trend.hidro = flag(seed++);
    trend.fric = flag(seed++);
    trend.dengD = flag(seed++);
    trend.dengL = flag(seed++);
    trend.mlFonte = flag(seed++);
    trend.mgFonte = flag(seed++);
    trend.mcFonte = flag(seed++);
    trend.dpB = flag(seed++);
    trend.potB = flag(seed++);
    trend.tempChokeJus = flag(seed++);
    trend.reyi = flag(seed++);
    trend.reye = flag(seed++);
    trend.Fr = flag(seed++);
    trend.grashi = flag(seed++);
    trend.grashe = flag(seed++);
    trend.nusi = flag(seed++);
    trend.nuse = flag(seed++);
    trend.hi = flag(seed++);
    trend.he = flag(seed++);
    trend.pri = flag(seed++);
    trend.pre = flag(seed++);
    trend.Rs = flag(seed++);
    trend.Bo = flag(seed++);
    trend.volMonM1PT = flag(seed++);
    trend.volJusM1PT = flag(seed++);
    trend.volMonM1ST = flag(seed++);
    trend.volJusM1ST = flag(seed++);
    trend.autoVal = flag(seed++);
    trend.autoVel = flag(seed++);
    trend.flutuacao = flag(seed++);
    trend.diamInt = flag(seed++);
    trend.TempParede = flag(seed++);
    trend.dadosParafina = flag(seed++);
    trend.inventarioGas = flag(seed++);
    trend.inventarioLiq = flag(seed++);
    trend.subResfria = flag(seed++);
}

void fillService(detTRENDG &trend, int (*flag)(int), int seed, int comp) {
    trend.dt = 2.5;
    trend.posic = 1;
    trend.comp = comp;
    trend.rotulo = "rotulo-G";
    trend.pres = flag(seed++);
    trend.temp = flag(seed++);
    trend.ugs = flag(seed++);
    trend.ug = flag(seed++);
    trend.tens = flag(seed++);
    trend.viscg = flag(seed++);
    trend.rhog = flag(seed++);
    trend.masg = flag(seed++);
    trend.masl = flag(seed++);
    trend.hidro = flag(seed++);
    trend.fric = flag(seed++);
    trend.FVHG = flag(seed++);
    trend.calor = flag(seed++);
    trend.qgst = flag(seed++);
    trend.velgarg = flag(seed++);
    trend.pEstagVGL = flag(seed++);
    trend.tEstagVGL = flag(seed++);
    trend.pGargVGL = flag(seed++);
    trend.tGargVGL = flag(seed++);
    trend.qVGL = flag(seed++);
    trend.reyi = flag(seed++);
    trend.reye = flag(seed++);
    trend.grashi = flag(seed++);
    trend.grashe = flag(seed++);
    trend.nusi = flag(seed++);
    trend.nuse = flag(seed++);
    trend.hi = flag(seed++);
    trend.he = flag(seed++);
    trend.pri = flag(seed++);
    trend.pre = flag(seed++);
    trend.diamInt = flag(seed++);
    trend.TempParede = flag(seed++);
    trend.subResfria = flag(seed++);
}

} // namespace

int main(int argc, char **argv) {
    if (argc != 2) {
        std::cerr << "usage: trend-sweep <output-directory>\n";
        return 2;
    }
    const std::string root(argv[1]);
    mkdir(root.c_str(), 0755);
    // The writers record absolute file names into arqRelatorioPerfis. Running
    // from inside the output root keeps those names relative, so two runs in
    // two scratch directories are byte-comparable -- otherwise the harness
    // reports a difference on every single _reported.txt.
    if (chdir(root.c_str()) != 0) {
        std::cerr << "cannot enter " << root << "\n";
        return 2;
    }

    const int kSeries = 2;
    const int kRows = 6;

    int configuration = 0;
    for (int branch = 0; branch < 2; branch++)
    for (int ap = 0; ap < 2; ap++)
    for (int pass = 1; pass <= 2; pass++)
    for (int language = 0; language < 2; language++)
    for (int injector = 0; injector < 2; injector++)
    for (int captions = 0; captions < 3; captions++) {
        char stem[64];
        snprintf(stem, sizeof stem, "c%03d/", configuration++);
        const std::string directory = stem;
        mkdir(directory.substr(0, directory.size() - 1).c_str(), 0755);
        pathPrefixoArqSaida = directory;

        // Deliberately leaked. ~Ler() frees arrays keyed off the counters --
        // with lingas > 0 it frees gasinj.tempo, which the default constructor
        // never allocated. Setting the counters is the point of this harness,
        // so the deck outlives the sweep instead. The process is short-lived
        // and writes files; nothing here depends on the memory being returned.
        Ler &input = *new Ler();
        input.ntendp = kSeries;
        input.ntendg = kSeries;
        input.ntendtransp = kSeries;
        input.ntendtransg = kSeries;
        input.lingas = 1;
        input.AP = ap;
        input.idiomaSaida = language;
        input.tipoSimulacao = injector ? tipoSimulacao_t::poco_injetor
                                       : tipoSimulacao_t::transiente;
        input.ncelp = 8;
        input.ncelg = 8;
        input.celp = new detcelp[input.ncelp];
        input.celg = new detcelg[input.ncelg];
        for (int c = 0; c < input.ncelp; c++)
            input.celp[c].dx = 12.5 + c * 0.125;
        for (int c = 0; c < input.ncelg; c++)
            input.celg[c].dx = 9.75 + c * 0.0625;

        input.trendp = new detTRENDP[kSeries];
        input.trendg = new detTRENDG[kSeries];
        input.trendtransp = new detTRENDTrans[kSeries];
        input.trendtransg = new detTRENDTrans[kSeries];
        input.nvartrendp = new int[kSeries];
        input.nvartrendg = new int[kSeries];

        int (*flag)(int) = (captions == 2) ? allOn : pattern;
        for (int s = 0; s < kSeries; s++) {
            // Series 1 uses a position above 1e6 on purpose. `comp` is an int,
            // and the header writer rounds it while the row writer does not
            // (anomaly A4-03). Below 1e6 both spellings print the same text and
            // the asymmetry is invisible; at 2000000 the double prints 2e+06
            // and the int prints 2000000, so header and rows land in different
            // files -- which is what the anomaly actually does.
            const int productionPosition = s ? 2000000 : 7;
            const int servicePosition = s ? 3000000 : 10;
            fillProduction(input.trendp[s], flag, captions * 11 + s, productionPosition);
            fillService(input.trendg[s], flag, captions * 13 + s, servicePosition);
            if (captions == 0) {
                // All captions off: the header collapses to the fixed columns.
                fillProduction(input.trendp[s], allOn, 0, productionPosition);
                detTRENDP &t = input.trendp[s];
                t.pres = t.temp = t.hol = t.FVH = t.bet = t.ugs = t.uls = 0;
                t.ug = t.ul = t.arra = t.autoVal = t.autoVel = t.flutuacao = 0;
                t.dadosParafina = 0;
            }
            input.trendtransp[s].camada = 2 + s;
            input.trendtransp[s].discre = 3 + s;
            input.trendtransp[s].comp = 40 + s * 5;
            input.trendtransp[s].rotulo = "rotulo-transP";
            input.trendtransg[s].camada = 4 + s;
            input.trendtransg[s].discre = 1 + s;
            input.trendtransg[s].comp = 60 + s * 5;
            input.trendtransg[s].rotulo = "rotulo-transG";
            input.nvartrendp[s] = 4 + s;
            input.nvartrendg[s] = 3 + s;
        }

        double ***productionBuffer = makeBuffer(kSeries, kRows, 16, 1);
        double ***serviceBuffer = makeBuffer(kSeries, kRows, 16, 2);
        double ***productionCross = makeBuffer(kSeries, kRows, 2, 3);
        double ***serviceCross = makeBuffer(kSeries, kRows, 2, 4);

        // Windows deliberately narrower than the buffer, and with a non-zero
        // base, so the base offset that the cross-section writers ignore
        // (anomaly A4-04) shows up as a difference if anyone "fixes" it.
        int productionCount[kSeries] = {4, 5};
        int productionBase[kSeries] = {1, 2};
        int serviceCount[kSeries] = {5, 4};
        int serviceBase[kSeries] = {2, 1};
        int crossCount[kSeries] = {4, 3};
        int crossBase[kSeries] = {1, 1};

        varGlob1D globals;
        globals.sequenciaAP = 3 + configuration % 4;

        const int branchIndex = branch ? 5 : -1;
        const double printPassCount = pass;

        const trendoutput::TrendState state{
            .inputData = input,
            .globals = &globals,
            .branchIndex = branchIndex,
            .printPassCount = printPassCount,
            .production = {productionBuffer, productionCount, productionBase},
            .service = {serviceBuffer, serviceCount, serviceBase},
            .productionCrossSection = {productionCross, crossCount, crossBase},
            // Bound to the PRODUCTION counters, mirroring trendStateOf: the
            // service cross-section writer reading the wrong group is baseline
            // behaviour, and the harness has to reproduce it to be comparable.
            .serviceCrossSection = {serviceCross, crossCount, crossBase}};

        arqRelatorioPerfis.open((directory + "_reported.txt").c_str(), ios_base::out);
        for (int s = 0; s < kSeries; s++) {
            const int network = 2 + s;
            trendoutput::writeProductionTrendHeader(state, s, network);
            trendoutput::writeProductionTrendRows(state, s, network);
            trendoutput::writeServiceTrendHeader(state, s, network);
            trendoutput::writeServiceTrendRows(state, s, network);
            trendoutput::writeProductionCrossSectionTrendHeader(state, s);
            trendoutput::writeProductionCrossSectionTrendRows(state, s);
            trendoutput::writeServiceCrossSectionTrendHeader(state, s);
            trendoutput::writeServiceCrossSectionTrendRows(state, s);
        }
        arqRelatorioPerfis.close();

        for (double ***buffer : {productionBuffer, serviceBuffer, productionCross, serviceCross}) {
            for (int s = 0; s < kSeries; s++) {
                for (int r = 0; r < kRows; r++)
                    delete[] buffer[s][r];
                delete[] buffer[s];
            }
            delete[] buffer;
        }
    }
    std::cout << configuration << " configurations swept\n";
    return 0;
}
