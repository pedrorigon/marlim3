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

enum class PermMode {
    cold,
    energyBase,
    semTermo,
    signedFlow,
    interfaceState,
    pressureAccessory,
    kineticEnergy,
    gasSource,
    latentPositive,
    latentNegative,
    latentDisabled,
    latentNan,
    coupledNetwork,
    velocityCap,
    annulusFirst,
    annulusLater,
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

TransCal seededHeatTransfer(varGlob1D *globals, double diameter,
                            double insideTemperature,
                            double outsideTemperature) {
    double conductivity[1] = {14.};
    double outsideDiameter[1] = {diameter + 0.02};
    double heatCapacity[1] = {500.};
    double density[1] = {7800.};
    double viscosity[1] = {1.};
    double expansion[1] = {1.e-5};
    int materialType[1] = {0};
    int materialIndex[1] = {1};
    DadosGeo geometry(diameter, 0., 0.12, 4.5e-5, 0, 1,
                      conductivity, outsideDiameter, heatCapacity, density,
                      viscosity, expansion, materialType, materialIndex);
    int layerCells[1] = {1};
    double radialStep[1] = {0.01};
    double layerTemperature[2] = {insideTemperature, outsideTemperature};
    double *layerTemperatures[1] = {layerTemperature};

    return TransCal(globals, geometry, 1, layerCells, radialStep,
                    layerTemperatures, insideTemperature, outsideTemperature,
                    outsideTemperature, 0.5, 0.25, 1, 0.5, 0.6, 2200., 850.,
                    2., 0.15, 1100., 900., 1., 0, 0, -100., 1., 1., 1., 65.,
                    10.);
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
        cell.MC = 1.1 * scale;
        cell.Mliqini = 0.65 * scale;
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
        cell.fluxcalmed = -880. - index;
        cell.dpB = 0.;
        cell.potBT = 0.;
        cell.potTermo = 0.;
        cell.fonteCal = 0.;
        cell.resAcopRedeP = 0.;
        cell.term1 = -981.;
        cell.term2 = -982.;
        cell.term1L = -983.;
        cell.term2L = -984.;
        cell.term1R = -985.;
        cell.term2R = -986.;
        cell.c0 = 1.15;
        cell.ud = 0.25;
        cell.arranjo = 1;
        cell.calor = seededHeatTransfer(system.vg1dSP, cell.duto.a,
                                        cell.temp, 20. + index +
                                                       0.1 * scenario.temperature);
        cell.calor.Vextern1 = 0.25;
        cell.calor.kextern1 = 0.15;
        cell.calor.cpextern1 = 1100.;
        cell.calor.rhoextern1 = 900.;
        cell.calor.viscextern1 = 0.001;
        cell.calor.fluxIni = 0.;
        cell.calor.fluxFim = 0.;
        cell.calor.resGlob = 0.;
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
    system.iterperm = 0;
    system.CalcLat = 0;
    system.arq.flashCompleto = 0;
    system.arq.latente = 0;
    system.arq.condlatente = 1;
    system.arq.limTransMass = 10.;
    system.arq.AceleraConvergPerm = 0;
    system.arq.nCompTotalUnidadesP = 100.;
    system.verificaAcopRedeP = 0;
    system.vg1dSP->blackOilTemp = 0;
    system.arq.chokep.abertura[0] = scenario.chokeOpening;
}

void resetAnnulus(SProd &system, const Scenario &scenario) {
    for (int index = 0; index < kCells + 3; ++index) {
        CelG &cell = system.celulaG[index];
        const double scale = 1. + 0.03 * index;
        cell.flui = fallbackFluid(system.vg1dSP);
        cell.temp = scenario.temperature + 8. + 0.1 * index;
        cell.VGasR = 0.4 * scale;
        cell.u1L = 2. * scale;
        cell.fluxcal = -801. - index;
        cell.calor = seededHeatTransfer(system.vg1dSP, 0.16 * scale,
                                        cell.temp, 18. + index);
        cell.calor.Tint = cell.temp;
        cell.calor.Vint = 0.5;
        cell.calor.kint = 0.6;
        cell.calor.cpint = 2200.;
        cell.calor.rhoint = 850.;
        cell.calor.viscint = 2.;
        cell.calor.Textern1 = 18. + index;
        cell.calor.Vextern1 = 0.25;
        cell.calor.kextern1 = 0.15;
        cell.calor.cpextern1 = 1100.;
        cell.calor.rhoextern1 = 900.;
        cell.calor.viscextern1 = 0.001;
        cell.calor.fluxIni = 0.;
        cell.calor.fluxFim = 0.;
        cell.calor.resGlob = 0.;
    }
}

void printPermState(const char *method, const char *scenario, const Cel &current,
                    const Cel &source, const CelG *annulus) {
    const TransCal &heat = source.calor;
    const double gasFlux = annulus ? annulus->fluxcal : -0.;
    const double gasTint = annulus ? annulus->calor.Tint : -0.;
    const double gasVint = annulus ? annulus->calor.Vint : -0.;
    const double gasK = annulus ? annulus->calor.kint : -0.;
    const double gasCp = annulus ? annulus->calor.cpint : -0.;
    const double gasRho = annulus ? annulus->calor.rhoint : -0.;
    const double gasVisc = annulus ? annulus->calor.viscint : -0.;
    const double gasText = annulus ? annulus->calor.Textern1 : -0.;
    const double gasWall0 = annulus ? annulus->calor.Tcamada[0][0] : -0.;
    const double gasWall1 = annulus ? annulus->calor.Tcamada[0][1] : -0.;
    printf("%-20s %-15s temp=%a velocity=%a phase=%a flux=%a "
           "Tint=%a Vint=%a k=%a cp=%a rho=%a visc=%a Text=%a Vext=%a "
           "kext=%a cpext=%a rhoext=%a viscext=%a res=%a wall0=%a wall1=%a "
           "gasFlux=%a gasTint=%a gasVint=%a gasK=%a gasCp=%a gasRho=%a "
           "gasVisc=%a gasText=%a gasWall0=%a gasWall1=%a\n",
           method, scenario, current.temp, current.VTemper,
           source.FonteMudaFase, source.fluxcalmed, heat.Tint, heat.Vint,
           heat.kint, heat.cpint, heat.rhoint, heat.viscint, heat.Textern1,
           heat.Vextern1, heat.kextern1, heat.cpextern1, heat.rhoextern1,
           heat.viscextern1, heat.resGlob, heat.Tcamada[0][0],
           heat.Tcamada[0][1], gasFlux, gasTint, gasVint, gasK, gasCp,
           gasRho, gasVisc, gasText, gasWall0, gasWall1);
}

void configurePermCase(SProd &system, Cel *cells, const Scenario &scenario,
                       int cellIndex, PermMode mode) {
    if (mode == PermMode::cold) {
        system.semTermo = 1;
        return;
    }

    system.trocaTermicaLenta = 0.;
    if (mode == PermMode::semTermo) {
        system.semTermo = 1;
    } else if (mode == PermMode::signedFlow) {
        for (int index = cellIndex; index <= cellIndex + 1; ++index) {
            cells[index].QG = -fabs(cells[index].QG);
            cells[index].QL = -fabs(cells[index].QL);
        }
    } else if (mode == PermMode::interfaceState) {
        system.iterperm = 1;
        cells[cellIndex].tempL -= 9.;
        cells[cellIndex].tempR += 13.;
    } else if (mode == PermMode::pressureAccessory) {
        cells[cellIndex - 1].acsr.tipo = 4;
        cells[cellIndex + 1].acsr.tipo = 4;
        cells[cellIndex - 1].acsr.bcs.freq = 45.;
        cells[cellIndex + 1].acsr.bcs.freq = 45.;
        cells[cellIndex].dpB = 7.5;
    } else if (mode == PermMode::kineticEnergy) {
        cells[cellIndex].MC = 2.4;
        cells[cellIndex].Mliqini = 0.9;
        cells[cellIndex - 1].QG *= 0.35;
        cells[cellIndex - 1].QL *= 0.45;
        cells[cellIndex + 1].QG *= 1.8;
        cells[cellIndex + 1].QL *= 1.6;
    } else if (mode == PermMode::gasSource) {
        const ProFlu fluid = fallbackFluid(system.vg1dSP);
        cells[cellIndex - 1].acsr.injg =
            InjGas(100., scenario.temperature + 25., fluid);
        cells[cellIndex + 1].acsr.injg =
            InjGas(100., scenario.temperature + 25., fluid);
        cells[cellIndex - 1].acsr.tipo = 1;
        cells[cellIndex + 1].acsr.tipo = 1;
        cells[cellIndex - 1].fontemassGR = 0.8;
        cells[cellIndex + 1].fontemassGR = 0.8;
        cells[cellIndex].flui.dVaporMassFraction = 0.4;
        cells[cellIndex - 1].FonteMudaFase = 7.;
        cells[cellIndex + 1].FonteMudaFase = 7.;
    } else if (mode == PermMode::latentPositive ||
               mode == PermMode::latentNegative ||
               mode == PermMode::latentDisabled ||
               mode == PermMode::latentNan) {
        system.CalcLat = 1;
        system.iterperm = 1;
        system.arq.latente = mode == PermMode::latentDisabled ? 0 : 1;
        cells[cellIndex].flui.dVaporMassFraction = 0.4;
        double source = mode == PermMode::latentPositive ||
                                mode == PermMode::latentDisabled
                            ? 25.
                            : -25.;
        if (mode == PermMode::latentNan)
            source = __builtin_nan("");
        cells[cellIndex - 1].FonteMudaFase = source;
        cells[cellIndex + 1].FonteMudaFase = source;
        if (mode == PermMode::latentNegative)
            system.arq.condlatente = 0;
    } else if (mode == PermMode::coupledNetwork) {
        system.verificaAcopRedeP = 1;
        system.PrimSecFimRedeP = cellIndex - 1;
        system.PrimSecIniRedeP = cellIndex - 1;
        cells[cellIndex - 1].resAcopRedeP = 2.5;
        cells[cellIndex - 1].potTermo = 1400.;
        cells[cellIndex - 1].fonteCal = 350.;
        cells[cellIndex + 1].potBT = 1750.;
    } else if (mode == PermMode::velocityCap) {
        system.vg1dSP->blackOilTemp = 1;
        cells[cellIndex].QG = 0.9;
        cells[cellIndex].QL = 0.8;
        cells[cellIndex + 1].QG = 1.1;
        cells[cellIndex + 1].QL = 1.0;
    } else if (mode == PermMode::annulusFirst ||
               mode == PermMode::annulusLater) {
        resetAnnulus(system, scenario);
        system.arq.lingas = 1;
        system.ColunaAnulaIni = 4;
        system.ColunaAnulaFim = 0;
        system.AnulaColunaIni = 4;
        if (mode == PermMode::annulusLater) {
            system.iterperm = 1;
            cells[cellIndex - 1].dTdLCor = 1.e5;
            cells[cellIndex].dTdLCor = 1.e5;
            cells[cellIndex + 1].dTdLCor = 1.e5;
        }
    }
}

void runPermCase(SProd &system, Cel *cells, const Scenario &scenario,
                 const char *caseName, int cellIndex, int RK, PermMode mode) {
    resetCells(system, cells, scenario);
    configurePermCase(system, cells, scenario, cellIndex, mode);
    system.RenovaTempPerm(cellIndex, RK);
    const CelG *forwardGas = system.arq.lingas == 1
                                 ? &system.celulaG[system.ColunaAnulaIni +
                                                  system.AnulaColunaIni -
                                                  (cellIndex - 1)]
                                 : nullptr;
    printPermState("RenovaTempPerm", caseName, cells[cellIndex],
                   cells[cellIndex - 1], forwardGas);

    resetCells(system, cells, scenario);
    configurePermCase(system, cells, scenario, cellIndex, mode);
    system.RenovaTempPermRev(cellIndex, RK);
    const CelG *reverseGas = system.arq.lingas == 1
                                 ? &system.celulaG[system.ColunaAnulaIni +
                                                  system.AnulaColunaIni -
                                                  (cellIndex + 1)]
                                 : nullptr;
    printPermState("RenovaTempPermRev", caseName, cells[cellIndex],
                   cells[cellIndex + 1], reverseGas);
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

// tempDescarga is reached in production only through arq.descarga == 1 ->
// resolveDescarga(), a gas-lift unloading run that no demo model enables. The
// corpus therefore never executes it and the L2 layer cannot speak for it,
// which is why it is driven here directly.
//
// The observable state is split across two arrays on purpose: the routine
// writes Tint, Vint, kint and cpint into celulaG[i-1].calor but rhoint and
// viscint into celula[i-1].calor. That asymmetry is preserved from the legacy
// body and is exactly the kind of thing a move must not quietly tidy up, so
// both destinations are printed.
void runDischargeCase(SProd &system, Cel *cells, const Scenario &scenario,
                      const char *caseName, int gasIndex,
                      double currentRatio, double previousRatio) {
    resetCells(system, cells, scenario);
    resetAnnulus(system, scenario);
    for (int index = gasIndex - 1; index <= gasIndex; ++index) {
        CelG &cell = system.celulaG[index];
        const double scale = 1. + 0.05 * index;
        cell.duto.a = 0.12 * scale;
        cell.duto.area = M_PI * cell.duto.a * cell.duto.a / 4.;
        cell.dx0 = 9. * scale;
        cell.dxL = 8. * scale;
        cell.pres = scenario.pressure * scale;
        cell.VGasL = 0.55 * scale;
    }
    system.celulaG[gasIndex].razInter = currentRatio;
    system.celulaG[gasIndex - 1].razInter = previousRatio;

    system.tempDescarga(gasIndex);

    const CelG &source = system.celulaG[gasIndex - 1];
    const Cel &production = cells[gasIndex - 1];
    printf("%-20s %-15s Tint=%a Vint=%a k=%a cp=%a rho=%a visc=%a prevR=%a nextL=%a\n",
           "tempDescarga", caseName, source.calor.Tint, source.calor.Vint,
           source.calor.kint, source.calor.cpint, production.calor.rhoint,
           production.calor.viscint, source.tempR,
           system.celulaG[gasIndex].tempL);
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

    runPermCase(system, cells, scenario, scenario.name, 2, 0,
                PermMode::cold);

    if (scenario.pressure == 65.) {
        runPermCase(system, cells, scenario, "energy-base", 3, 0,
                    PermMode::energyBase);
        runPermCase(system, cells, scenario, "sem-termo", 3, 0,
                    PermMode::semTermo);
        runPermCase(system, cells, scenario, "signed-flow", 3, 0,
                    PermMode::signedFlow);
        runPermCase(system, cells, scenario, "interface-rk1", 3, 1,
                    PermMode::interfaceState);
        runPermCase(system, cells, scenario, "pressure-bcs", 3, 0,
                    PermMode::pressureAccessory);
        runPermCase(system, cells, scenario, "kinetic", 3, 0,
                    PermMode::kineticEnergy);
        runPermCase(system, cells, scenario, "gas-source", 3, 0,
                    PermMode::gasSource);
        runPermCase(system, cells, scenario, "latent-positive", 3, 0,
                    PermMode::latentPositive);
        runPermCase(system, cells, scenario, "latent-negative", 3, 0,
                    PermMode::latentNegative);
        runPermCase(system, cells, scenario, "latent-disabled", 3, 0,
                    PermMode::latentDisabled);
        runPermCase(system, cells, scenario, "latent-nan", 3, 0,
                    PermMode::latentNan);
        runPermCase(system, cells, scenario, "network-shaft", 3, 0,
                    PermMode::coupledNetwork);
        runPermCase(system, cells, scenario, "velocity-cap", 3, 0,
                    PermMode::velocityCap);
        runPermCase(system, cells, scenario, "annulus-first", 3, 0,
                    PermMode::annulusFirst);
        runPermCase(system, cells, scenario, "annulus-later", 3, 0,
                    PermMode::annulusLater);
    }

    runDischargeCase(system, cells, scenario, scenario.name, 3, 0.30, 0.60);
    if (scenario.pressure == 65.) {
        // razInter below localtiny takes the liquid-density branch of vel1;
        // the 0.8 and 0.5 pairs take the two ratio guards at and past 0.5.
        runDischargeCase(system, cells, scenario, "disch-zero", 3, 0., 0.);
        runDischargeCase(system, cells, scenario, "disch-gas", 3, 0.8, 0.8);
        runDischargeCase(system, cells, scenario, "disch-half", 3, 0.5, 0.5);
    }

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
    system.celulaG = new CelG[kCells + 3];
    initializeLatentHeatTable(system);

    for (const Scenario &scenario : kScenarios)
        runScenario(system, cells, scenario);
    return 0;
}
