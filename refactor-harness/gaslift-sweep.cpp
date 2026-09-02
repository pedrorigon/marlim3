/*
 * Dedicated characterization sweep for the Stage 6 gas-line and gas-lift
 * surface.
 *
 * Coverage measurement before the stage opened found that the demo corpus never
 * executes twelve of the twenty-two methods this stage moves -- roughly 600
 * lines -- so for those the L2 and L3 gates compare the output of code that did
 * not run and stay green whatever an extraction does to it. This driver reaches
 * them through SProd's public API, on freshly seeded state, printing every
 * double in hexadecimal so the table is exact for binary64.
 *
 * It exists for the same reason thermal-sweep.cpp does, and is built on the same
 * scaffolding: the fallback fluid, the guard-cell arrangement and the seeded
 * heat-transfer object are shared, because two harnesses that disagree about
 * what a valid cell looks like would be two harnesses nobody can compare.
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
        // The thermal sweep never needed fluicol seeded because the paths it
        // drives read it only through cached cell fields. delpInjPerm calls it
        // directly, and a default-constructed ProFluCol returns NaN -- a row
        // that is deterministic and proves nothing, since almost any corruption
        // would also produce NaN. Every argument is a plain physical constant,
        // so the row is reproducible on any machine.
        // compresP and compresT must be NON-ZERO. Left at zero, MasEspFlu
        // ignores pressure and temperature and returns a constant, which made
        // delpInjPerm produce the same value in every scenario AND survive a
        // corruption of its temperature interpolation -- a row that ran and
        // proved nothing.
        cell.fluicol = ProFluCol(1000., 4.5e-10, 7.0e-4, 1e-6, 4000., 1., 20.,
                                 1., 60., 0.5, 0., 0, 0);
        cell.duto.revest = 0;
        cell.dutoL.revest = 0;
        // Both gradient switches multiply the friction and hydrostatic terms of
        // delpInjPerm. Left at zero they annihilate the scenario, which is what
        // made that row identical in all four.
        cell.dPdLFric = 1.;
        cell.dPdLHidro = 1.;
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

// ---------------------------------------------------------------- gas line --
// The gas cells the thermal sweep seeds carry only what the annulus coupling
// needed. The gas line reads a good deal more, so the fields below are set
// explicitly rather than left at whatever the constructor produced -- an
// uninitialised read would make the table depend on the allocator, and a table
// that is not reproducible proves nothing.
constexpr int kGasCells = 5;
constexpr int kValves = 2;

void resetGasLine(SProd &system, const Scenario &scenario) {
    system.ncelGas = kGasCells;
    for (int index = 0; index <= kGasCells + 1; ++index) {
        CelG &cell = system.celulaG[index];
        const double scale = 1. + 0.04 * index;
        cell.flui = fallbackFluid(system.vg1dSP);
        cell.duto.a = 0.09 * scale;
        cell.duto.area = M_PI * cell.duto.a * cell.duto.a / 4.;
        cell.duto.peri = M_PI * cell.duto.a;
        cell.duto.rug = 4.5e-5;
        cell.duto.teta = 0.15;
        cell.duto.revest = 0;
        cell.dx0 = 11. * scale;
        cell.dxL = 10. * scale;
        cell.dt = 0.5;
        cell.pres = scenario.pressure * scale;
        cell.presini = (scenario.pressure - 1.) * scale;
        cell.presL = scenario.pressure * scale;
        cell.presR = (scenario.pressure + 1.) * scale;
        cell.temp = scenario.temperature + 6. + 0.2 * index;
        cell.tempini = scenario.temperature + 5. + 0.2 * index;
        cell.tempL = scenario.temperature + 5.5 + 0.2 * index;
        // dutoL and the two gradient switches are read by delpGasPerm on its
        // second half. Left unseeded, dutoL.a is zero, the flow area is zero and
        // the routine returns NaN -- deterministic, and unable to distinguish
        // one corruption from another.
        cell.dutoL = cell.duto;
        cell.dPdLFric = 1.;
        cell.dPdLHidro = 1.;
        // Velocities carry the scenario, not just the index. Seeding them from
        // `scale` alone made delpGasPerm and delpInjPerm return the SAME value
        // in all four scenarios: the rows were deterministic and blind.
        cell.VGasL = (0.42 + scenario.gasFlow) * scale;
        cell.VGasR = (0.44 + scenario.gasFlow) * scale;
        cell.VGasRR = (0.46 + scenario.gasFlow) * scale;
        cell.u1L = (2. + scenario.liquidFlow) * scale;
        cell.rg = 0.8 * scale;
        cell.razInter = 0.3;
        cell.fluxcal = -801. - index;
        cell.calor = seededHeatTransfer(system.vg1dSP, 0.16 * scale,
                                        cell.temp, 18. + index);
    }
    for (int index = 0; index < 3 * (kGasCells + 2); ++index)
        system.termolivreG[index] = 40. + 0.5 * index;
}

void resetValves(SProd &system) {
    for (int valve = 0; valve < kValves; ++valve) {
        system.posicVGLG[valve] = valve + 1;
        system.posicVGLP[valve] = valve + 2;
        ChokeGas &choke = system.chokeVGL[valve];
        choke.presEstag = 90. + 5. * valve;
        choke.tempEstag = 60. + 2. * valve;
        choke.presGarg = 70. + 5. * valve;
        choke.areagarg = 1.2e-4 * (1. + 0.3 * valve);
        choke.flui = fallbackFluid(system.vg1dSP);
    }
}

void printGas(const char *method, const char *scenario, const CelG &cell) {
    printf("%-28s %-15s pres=%a presL=%a presR=%a VGasL=%a VGasR=%a "
           "temp=%a tempini=%a rg=%a\n",
           method, scenario, cell.pres, cell.presL, cell.presR,
           cell.VGasL, cell.VGasR, cell.temp, cell.tempini, cell.rg);
}

void printValue(const char *method, const char *scenario, double value) {
    printf("%-28s %-15s value=%a\n", method, scenario, value);
}

void runGasScenario(SProd &system, Cel *cells, const Scenario &scenario) {
    // areaValvCali is the one routine here with no state at all: eight doubles
    // in, one out. Four rows pin both sides of its opening guard.
    // The opening fraction is very nearly a binary switch: over most of the
    // input space APE either falls below the seat area (closed) or overshoots it
    // by a wide margin and is capped to 1. A calibration that corrupts the cap
    // is then undetectable, which is exactly what the first run of this harness
    // reported. Casing pressure just above the bellows threshold puts the "base"
    // scenario in the narrow band where APE exceeds the seat area WITHOUT
    // swamping it, so the cap has something to do; the other three pin the
    // saturated and closed behaviour, which is coverage of its own.
    printValue("areaValvCali", scenario.name,
               system.areaValvCali(scenario.pressure, scenario.temperature,
                                   scenario.pressure, scenario.pressure * 0.95,
                                   0.02, 1.0e-4, 0.05, scenario.temperature + 20.));
    printValue("areaValvCali-closed", scenario.name,
               system.areaValvCali(scenario.pressure * 3., scenario.temperature,
                                   scenario.pressure * 0.2, scenario.pressure * 0.1,
                                   0.05, 1.2e-4, 0.12, scenario.temperature));

    resetCells(system, cells, scenario);
    resetGasLine(system, scenario);
    resetValves(system);
    system.renovaGas();
    printGas("renovaGas", scenario.name, system.celulaG[2]);

    resetGasLine(system, scenario);
    system.renovaGasBuf();
    printGas("renovaGasBuf", scenario.name, system.celulaG[2]);

    resetCells(system, cells, scenario);
    resetGasLine(system, scenario);
    resetValves(system);
    printValue("prescordesc", scenario.name,
               system.prescordesc(0.4, 0, 1.1, 1));

    resetCells(system, cells, scenario);
    resetGasLine(system, scenario);
    resetValves(system);
    printValue("delpGasPerm", scenario.name, system.delpGasPerm(2));

    resetCells(system, cells, scenario);
    printValue("delpInjPerm", scenario.name, system.delpInjPerm(2));

    resetCells(system, cells, scenario);
    resetGasLine(system, scenario);
    resetValves(system);
    system.RenovaPresGasPerm(2);
    printGas("RenovaPresGasPerm", scenario.name, system.celulaG[2]);

    resetCells(system, cells, scenario);
    resetGasLine(system, scenario);
    resetValves(system);
    system.RenovaTempGasPerm(2);
    printGas("RenovaTempGasPerm", scenario.name, system.celulaG[2]);
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
    system.ncel = kCells - 1;
    system.celulaG = new CelG[kGasCells + 3];
    system.termolivreG = Vcr<double>(3 * (kGasCells + 3), 0.);
    system.chokeVGL = new ChokeGas[kValves];
    system.posicVGLG = new int[kValves];
    system.posicVGLP = new int[kValves];

    for (const Scenario &scenario : kScenarios)
        runGasScenario(system, cells, scenario);
    return 0;
}
