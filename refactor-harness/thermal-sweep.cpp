/*
 * Dedicated characterization sweep for the Stage 5 thermal surface.
 *
 * The demo corpus never calls eleven of the nineteen methods moved by this
 * stage. This driver reaches them through SProd's public API, using the same
 * valid fallback fluid and guard-cell arrangement as c0ud-sweep.cpp. Every
 * invocation starts from freshly seeded state, and every double is printed in
 * hexadecimal notation so the table is exact for binary64.
 */
#include "Leitura.h"
#include "SisProd.h"
#include "estruturas.h"
#include "variaveisGlobais1D.h"

#include <cstdio>
#include <fstream>
#include <string>

using namespace std;

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

constexpr int kCells = 7;

struct Scenario {
    const char *name;
    double pressure;
    double temperature;
    double voidFraction;
    double gasFlow;
    double liquidFlow;
    double chokeOpening;
};

const Scenario kScenarios[] = {
    {"base", 65., 55., 0.30, 0.20, 0.35, 0.50},
    {"hot", 85., 95., 0.65, 0.45, 0.10, 0.40},
    {"reverse-gas", 45., 35., 0.15, -0.12, 0.50, 0.80},
    {"reverse-liquid", 110., 70., 0.80, 0.35, -0.18, 0.20},
};

ProFlu fallbackFluid(varGlob1D *globals) {
    return ProFlu(globals, 20, 100., 0.7, 0., 1., 20, 10, 40,
                  2, 0, 0.8, 0, 0, 0);
}

void resetCells(SProd &system, Cel *cells, const Scenario &scenario) {
    for (int index = -1; index < kCells; ++index) {
        Cel &cell = cells[index];
        const double scale = 1. + 0.04 * (index + 1);
        cell.flui = fallbackFluid(system.vg1dSP);
        cell.fluiL = &cell.flui;
        cell.duto.a = 0.10 * scale;
        cell.duto.area = M_PI * cell.duto.a * cell.duto.a / 4.;
        cell.duto.rug = 4.5e-5;
        cell.duto.teta = 0.12;
        cell.dutoL = cell.duto;
        cell.dutoR = cell.duto;
        cell.dx = 10. * scale;
        cell.dxL = 12. * scale;
        cell.dxR = 11. * scale;
        cell.dt = 0.5;
        cell.pres = scenario.pressure * scale;
        cell.presini = (scenario.pressure - 1.) * scale;
        cell.presaux = (scenario.pressure + 1.) * scale;
        cell.presauxL = scenario.pressure * scale;
        cell.presBuf = (scenario.pressure + 2.) * scale;
        cell.presRBuf = (scenario.pressure + 3.) * scale;
        cell.temp = scenario.temperature * scale;
        cell.tempini = (scenario.temperature - 2.) * scale;
        cell.tempL = (scenario.temperature - 1.) * scale;
        cell.tempR = (scenario.temperature + 1.) * scale;
        cell.alf = scenario.voidFraction;
        cell.alfini = scenario.voidFraction * 0.95;
        cell.alfL = scenario.voidFraction * 0.90;
        cell.alfR = scenario.voidFraction * 1.05;
        cell.alfPigE = scenario.voidFraction * 0.92;
        cell.alfPigD = scenario.voidFraction * 1.04;
        cell.alfPigER = scenario.voidFraction * 0.94;
        cell.bet = 0.12;
        cell.betini = 0.10;
        cell.betL = 0.11;
        cell.betR = 0.13;
        cell.betPigE = 0.14;
        cell.betPigD = 0.15;
        cell.QG = scenario.gasFlow * scale;
        cell.QL = scenario.liquidFlow * scale;
        cell.QGini = cell.QG * 0.9;
        cell.QLini = cell.QL * 0.9;
        cell.MR = 1.2 * scale;
        cell.MliqiniR = 0.7 * scale;
        cell.MliqiniBuf = scenario.liquidFlow * scale;
        cell.MCBuf = cell.MliqiniBuf + scenario.gasFlow * scale;
        cell.MliqiniLBuf = scenario.liquidFlow * 0.9 * scale;
        cell.MLBuf = cell.MliqiniLBuf + scenario.gasFlow * 0.8 * scale;
        cell.MliqiniRBuf = scenario.liquidFlow * 1.1 * scale;
        cell.MRBuf = cell.MliqiniRBuf + scenario.gasFlow * 1.2 * scale;
        cell.QLR = scenario.liquidFlow * 1.1 * scale;
        cell.rpC = cell.flui.MasEspLiq(cell.pres, cell.temp);
        cell.rcC = cell.fluicol.MasEspFlu(cell.pres, cell.temp);
        cell.rgC = cell.flui.MasEspGas(cell.pres, cell.temp);
        cell.mipC = cell.flui.ViscOleo(cell.pres, cell.temp);
        cell.micC = cell.fluicol.VisFlu(cell.pres, cell.temp);
        cell.migC = cell.flui.ViscGas(cell.pres, cell.temp);
        cell.VTemper = 0.25;
        cell.acsr.tipo = 0;
        cell.acsrL = nullptr;
        cell.fontemassLR = 0.;
        cell.fontemassCR = 0.;
        cell.fontemassGR = 0.;
        cell.fontemassLL = 0.;
        cell.fontemassCL = 0.;
        cell.fontemassGL = 0.;
        cell.FonteMudaFase = 0.;
        cell.term1 = -981.;
        cell.term2 = -982.;
        cell.term1L = -983.;
        cell.term2L = -984.;
        cell.term1R = -985.;
        cell.term2R = -986.;
        cell.c0 = 1.15;
        cell.ud = 0.25;
        cell.arranjo = 1;
        cell.calor.Tint = -991.;
        cell.calor.dtL = -992.;
        cell.calor.Vint = -993.;
        cell.calor.dt = -994.;
        cell.calor.kint = -995.;
        cell.calor.cpint = -996.;
        cell.calor.rhoint = -997.;
        cell.calor.viscint = -998.;
        cell.calor.betint = -999.;
        cell.calor.Textern1 = 20. + index + 0.1 * scenario.temperature;
    }

    system.ncel = kCells - 2;
    system.semTermo = 0;
    system.presfim = scenario.pressure - 5.;
    system.pGSup = scenario.pressure;
    system.tempSup = -901.;
    system.arq.master1.razareaativ = 0.;
    system.arq.lingas = 0;
    system.arq.escorregaTran = 1;
    system.arq.escorregaPerm = 1;
    system.arq.escorregamentoCelulaContorno = 1;
    system.arq.mapaArranjo = 0;
    system.driftSelectors.dispersed = 1;
    system.driftSelectors.annularChurn = 1;
    system.driftSelectors.stratified = 1;
    system.tGSup = scenario.temperature - 5.;
    system.presE = scenario.pressure;
    system.tempE = scenario.temperature;
    system.betaE = 0.13;
    system.titE = 0.20;
    system.alfE = scenario.voidFraction;
    system.trocaTermicaLenta = 1.e99;
    system.arq.chokep.abertura[0] = scenario.chokeOpening;
}

void initializeLatentHeatTable(SProd &system) {
    system.arq.tabent.npont = 3;
    system.HLat = new double *[4];
    for (int row = 0; row < 4; ++row)
        system.HLat[row] = new double[4];

    const double pressures[4] = {0., 20., 80., 140.};
    const double temperatures[4] = {0., 10., 70., 130.};
    for (int row = 0; row < 4; ++row) {
        system.HLat[row][0] = pressures[row];
        system.HLat[0][row] = temperatures[row];
    }
    for (int row = 1; row < 4; ++row)
        for (int column = 1; column < 4; ++column)
            system.HLat[row][column] = 100. * row + 7. * column;
}

double **allocateTable(int size) {
    double **table = new double *[size];
    for (int row = 0; row < size; ++row)
        table[row] = new double[size];
    return table;
}

void initializePropertyTables(Cel &cell, const Scenario &scenario) {
    // calcTempEntalp increments j before reading Var[0][j + 1]. Its legacy
    // `while (... || ...)` can advance one interval beyond the declared table,
    // so the harness provides two guard columns and records the behavior
    // without invoking undefined storage.
    constexpr int kTableSize = 6;
    cell.flui.npontos = 3;
    cell.flui.rhogF = allocateTable(kTableSize);
    cell.flui.rholF = allocateTable(kTableSize);
    cell.flui.HgF = allocateTable(kTableSize);
    cell.flui.HlF = allocateTable(kTableSize);

    for (int row = 0; row < kTableSize; ++row) {
        for (int column = 0; column < kTableSize; ++column) {
            const double pressure = 20. + 40. * row;
            const double temperature = 10. + 40. * column;
            cell.flui.rhogF[row][column] = 35. + 3. * row - 0.04 * temperature;
            cell.flui.rholF[row][column] = 900. + 4. * row - 0.20 * temperature;
            cell.flui.HgF[row][column] = 2.0e5 + 1500. * temperature + 80. * pressure;
            cell.flui.HlF[row][column] = 5.0e4 + 900. * temperature + 40. * pressure;
        }
        cell.flui.rhogF[row][0] = 20. + 40. * row;
        cell.flui.rholF[row][0] = 20. + 40. * row;
    }
    for (int column = 0; column < kTableSize; ++column) {
        cell.flui.rhogF[0][column] = 10. + 40. * column;
        cell.flui.rholF[0][column] = 10. + 40. * column;
    }

    // Keep the table tied to the scenario so every row is independently
    // discriminating rather than four labels over one calculation.
    cell.flui.HgF[1][2] += scenario.pressure * 13.;
    cell.flui.HlF[2][2] += scenario.temperature * 17.;
}

void runScenario(SProd &system, Cel *cells, const Scenario &scenario) {
    resetCells(system, cells, scenario);
    const double latentHeat = system.interpolaHLatente(
        scenario.pressure, scenario.temperature);
    printf("%-20s %-15s value=%a\n", "interpolaHLatente", scenario.name, latentHeat);

    resetCells(system, cells, scenario);
    const double mixtureEnthalpy = system.calcHmix(2);
    printf("%-20s %-15s value=%a\n", "calcHmix", scenario.name, mixtureEnthalpy);

    resetCells(system, cells, scenario);
    system.atualizaPeriTempProd(2);
    printf("%-20s %-15s prevR=%a nextL=%a currentIni=%a\n",
           "atualizaPeriTempProd", scenario.name,
           cells[1].tempR, cells[3].tempL, cells[2].tempini);

    resetCells(system, cells, scenario);
    system.calcTempFim();
    printf("%-20s %-15s tempSup=%a\n", "calcTempFim", scenario.name, system.tempSup);

    resetCells(system, cells, scenario);
    system.prepDifusCalorND(2);
    const TransCal &heat = cells[2].calor;
    printf("%-20s %-15s Tint=%a dtL=%a Vint=%a dt=%a k=%a cp=%a rho=%a visc=%a beta=%a\n",
           "prepDifusCalorND", scenario.name, heat.Tint, heat.dtL, heat.Vint,
           heat.dt, heat.kint, heat.cpint, heat.rhoint, heat.viscint, heat.betint);

    resetCells(system, cells, scenario);
    initializePropertyTables(cells[2], scenario);
    const double tabulatedEnergy = system.energmix(2, 1, 2, 0.35);
    printf("%-20s %-15s value=%a\n", "energmix", scenario.name, tabulatedEnergy);

    if (scenario.pressure >= 60. && scenario.pressure < 140.) {
        resetCells(system, cells, scenario);
        initializePropertyTables(cells[2], scenario);
        system.calcTempEntalp(2);
        printf("%-20s %-15s temp=%a\n", "calcTempEntalp", scenario.name, cells[2].temp);
    }

    resetCells(system, cells, scenario);
    system.semTermo = 1;
    system.RenovaTempPerm(2, 0);
    printf("%-20s %-15s temp=%a\n", "RenovaTempPerm", scenario.name, cells[2].temp);

    resetCells(system, cells, scenario);
    system.RenovaTempPermRev(2, 0);
    printf("%-20s %-15s temp=%a\n", "RenovaTempPermRev", scenario.name, cells[2].temp);

        resetCells(system, cells, scenario);
        system.renovatermAfluFim();
        const Cel &outlet = cells[system.ncel];
        printf("%-20s %-15s term1=%a term2=%a c0=%a ud=%a arr=%d left1=%a left2=%a\n",
            "renovatermAfluFim", scenario.name, outlet.term1, outlet.term2,
            outlet.c0, outlet.ud, static_cast<int>(outlet.arranjo),
            outlet.term1L, outlet.term2L);

        resetCells(system, cells, scenario);
        system.renovatermColIni();
        printf("%-20s %-15s term1=%a term2=%a c0=%a ud=%a arr=%d alfE=%a right1=%a right2=%a\n",
            "renovatermColIni", scenario.name, cells[0].term1, cells[0].term2,
            cells[0].c0, cells[0].ud, static_cast<int>(cells[0].arranjo),
            system.alfE, cells[1].term1L, cells[1].term2L);
}

}  // namespace

int main() {
    varGlob1D globals;
    SProd &system = *new SProd();
    system.vg1dSP = &globals;
    system.arq.chokep.parserie = 1;
    system.arq.chokep.abertura = new double[1];

    Cel *storage = new Cel[kCells + 1];
    Cel *cells = storage + 1;
    system.celula = cells;
    initializeLatentHeatTable(system);

    for (const Scenario &scenario : kScenarios)
        runScenario(system, cells, scenario);
    return 0;
}