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

// ------------------------------------------------------------- unloading ---
// The six routines driven by runUnloadingScenario are the ones NEITHER the demo
// corpus NOR the rest of this file executed. Measured with gcov, not estimated:
// the corpus runs 12 of the 22 functions, this sweep ran 8, and the union left
// six untouched -- 355 executable lines resting on token identity alone.
//
// They are all on the unloading path, which the product reaches only when
// arq.descarga == 1, and that comes from configuracaoInicial/condicaoInicial
// == 3. No model in the corpus sets it, so no model will ever reach them: the
// gap is a property of the corpus, not of how many models happen to be run.
//
// Seeding the unloading path is not seeding the steady one. The schedule is a
// state machine around the interface cell, and a run that leaves celInter at
// zero exercises the guards rather than the arithmetic -- a deterministic row
// that proves nothing, which is the trap resetCells already documents for
// fluicol. Two constraints follow:
//
//   - celInter must be at least 1. celulaG is allocated WITHOUT the guard-cell
//     offset that celula gets, so advanceInterface reading celulaG[celInter - 1]
//     at celInter == 0 walks off the front of the array.
//   - lixo5 must exceed tempoLatenciaDesc. Despite the name, lixo5 is not
//     scratch here: searchUnloadingInjectionPressure gates its entire body on
//     `lixo5 > tempoLatenciaDesc`, using it as the simulation clock.
constexpr int kInterfaceCell = 2;

// The gas line must stand ABOVE the tubing it feeds, or nothing is injected.
// resetGasLine and resetCells happen to give the valve cell and its production
// cell the same pressure, and searchUnloadingInjectionPressure gates its mass
// flow on `pmed < gasCells[valve].pres` -- so at parity every valve contributes
// exactly zero and the routine returns 0 in all four scenarios. That is a row
// that cannot fail. The factor restores the physical ordering; its value is
// arbitrary, its being greater than one is not.
constexpr double kInjectionOverPressure = 1.35;

void resetUnloading(SProd &system, const Scenario &scenario) {
    system.arq.descarga = 1;
    system.arq.controDesc = 1;
    system.arq.nvalvgas = kValves;
    system.arq.presIniDescG = scenario.pressure * 1.2;
    system.arq.presMaxDesc = 100000.;
    system.arq.tempoLatenciaDesc = 100.;
    // Grams per litre, NOT ppm. ChokeGas::MasEspFlu forms the mass fraction as
    // salin/1000, so a seawater figure in ppm (35000) makes the fraction 35,
    // the brine density negative, and massica the square root of a negative
    // number. The routine then returns NaN, `NaN > vazmax` is false, and
    // searchUnloadingInjectionPressure returns a clean 0 in every scenario --
    // a row that looks stable and tests nothing. Scenario-varying, so the four
    // rows differ; the range is ordinary completion brine.
    system.arq.salinDescarga = 35. + 40. * scenario.voidFraction;
    system.arq.gasinj.tipoCC = 0;
    system.arq.gasinj.presinj[0] = scenario.pressure * 1.1;
    (*system.vg1dSP).lixo5 = 5000.;

    system.celInter = kInterfaceCell;
    system.celInterIni = kInterfaceCell;
    system.velInter = 0.35 + scenario.gasFlow;
    system.velInterIni = 0.30 + scenario.gasFlow;
    system.dtInter = 0.4;
    system.dtInterIni = 0.45;

    system.pGSup = scenario.pressure * 1.15;
    // ABOVE the gas line, not below it. The injection choke sits upstream, and
    // ChokeGas::massica zeroes its own output when the stagnation pressure is
    // below the throat pressure. Seeded under the line, the choke delivered
    // exactly nothing, so the whole source term advanceBufferedGasSubStep
    // assembles was zero and moving the opening bound changed no digit.
    system.presiniG = scenario.pressure * kInjectionOverPressure * 1.15;
    system.tempiniG = scenario.temperature + 4.;
    system.dt = 0.5;

    // Sliding windows. advanceGasSubStep push_backs at the tail and erases the
    // front once the window passes maxVecContDesc, so both must start non-empty
    // and the bound must be small enough that the erase branch is reachable.
    system.vazmedDesc = 0.22 + scenario.gasFlow;
    system.tempmedDEsc = scenario.temperature + 3.;
    system.tempMedContDesc = scenario.temperature + 2.;
    system.maxVecContDesc = 3.;
    system.vazmaxMedDesc = {0.18, 0.20, 0.24};
    system.dtDesc = {0.5, 0.5, 0.5};

    // resetValves puts the valves at gas cells 1 and 2, which is right for the
    // steady half but leaves both at or below the interface here -- and
    // searchUnloadingInjectionPressure only looks at valves ABOVE it, so its
    // entire valve loop was skipping its body. Re-placed for this half only;
    // resetValves runs again before every measurement, so the steady rows are
    // untouched.
    for (int valve = 0; valve < kValves; ++valve) {
        system.posicVGLG[valve] = kInterfaceCell + 1 + valve;
        ChokeGas &choke = system.chokeVGL[valve];
        // frec must not be 1: it divides (1 - frec).
        choke.frec = 0.1 + 0.05 * valve;
        // areafole must not be 0: it divides areagarg.
        choke.areafole = 9.5e-4 * (1. + 0.2 * valve);
        choke.pcalib = scenario.pressure * 1.25;
        choke.tcalib = scenario.temperature;
        choke.dextern = 0.0381;
        // One calibrated valve and one orifice, so the calibratedValveArea
        // branch and the branch that skips it are both driven.
        choke.tipo = valve == 0 ? 1 : 0;
    }

    // ProFluCol::VisFlu branches on its own `descarga` flag: set, it uses the
    // completion-fluid correlation; clear, it falls into an ASTM dead-oil fit
    // whose LVisL/LVisH/TempL/TempH this harness never seeds, so it takes
    // log10 of a negative number and returns NaN. The NaN then propagates to
    // the mixture viscosity, Reynolds is NaN, the friction factor comes back 0,
    // and the wall-shear term of computeUnloadingValvePressure vanishes -- a
    // whole term of the pressure march silently absent from the table. The
    // product wires this the same way, at Leitura.cpp:12648.
    for (int index = -1; index < kCells; ++index) {
        Cel &cell = system.celula[index];
        cell.fluicol.descarga = 1;
        // duto.dia and termRug are the other half of the same hole. Cel::fric
        // reads the PRECOMPUTED termRug rather than deriving it, and the
        // harness never built one, so it held whatever new Cel[] left there --
        // NaN in practice. Computed exactly as celula3.cpp:333 does, from
        // duto.dia, which also has to exist: resetCells seeds duto.a and
        // derives the area from it, but nothing sets dia.
        cell.duto.dia = cell.duto.a;
        cell.termRug = pow(cell.duto.rug / cell.duto.dia / 3.7, 1.11);
        // And the wetted perimeter. resetGasLine seeds peri on the GAS cells but
        // resetCells never did on the production ones, so it was zero -- and the
        // wall-shear term of the pressure march reads
        // `tens1 * perimeter / flowArea`. At perimeter == 0 the entire shear
        // contribution is multiplied out of existence: setting tens1 to 12345.678
        // moved not one digit of the published table. The friction was not merely
        // mis-seeded, it was absent, and no row could ever have noticed.
        cell.duto.peri = M_PI * cell.duto.a;
    }

    // The control band computeUnloadingValvePressure corrects against. Both of
    // its branches compare the surface pressure against these, so leaving them
    // at zero pins the routine to its do-nothing path.
    system.arq.vazDescControl = 0.5;
    // BELOW the surface pressure. Above it, every branch that nudges pGSup is
    // immediately clamped back up to this floor, and the nudge -- including the
    // 1% decay -- never reaches the table.
    system.arq.presMinDesc = scenario.pressure * 0.90;
    system.arq.presMinDescG = scenario.pressure * 0.80;
    system.arq.presMaxDescG = scenario.pressure * 2.00;

    // presMaxDesc starts at 100000 inside the routine and only comes down
    // through an IPR accessory, so without one the first branch can never fire:
    // the surface pressure is never >= 100000. One reservoir accessory puts the
    // ceiling within reach.
    for (int index = 0; index < kCells; ++index)
        system.celula[index].acsr.tipo = 0;
    system.celula[1].acsr.tipo = 3;
    system.celula[1].acsr.ipr.Pres = scenario.pressure * 1.10;

    system.chokeInj.presEstag = scenario.pressure * 1.3;
    system.chokeInj.tempEstag = scenario.temperature + 8.;
    system.chokeInj.presGarg = scenario.pressure * 1.1;
    // Sized as a FRACTION of the pipe, so abertoChk lands around the 0.2 bound
    // advanceBufferedGasSubStep switches on rather than far below it, where
    // moving the bound changes nothing.
    system.chokeInj.areagarg =
        (0.15 + 0.3 * scenario.voidFraction) * system.celulaG[0].duto.area;
    system.chokeInj.flui = fallbackFluid(system.vg1dSP);

    for (int index = 0; index <= kGasCells + 1; ++index) {
        CelG &cell = system.celulaG[index];
        // celInter is a POINTER on the cell, aimed at SProd::celInter, and
        // CelG::GeraLocal dereferences it unconditionally. Unseeded it is null
        // and advanceBufferedGasSubStep segfaults -- which is how this line came
        // to be written. Same wiring the product does at SisProd.cpp:1483.
        cell.celInter = &system.celInter;
        cell.celInterini = &system.celInterIni;
        cell.posic = index;

        // The product's convention: fully gas ahead of the interface, fully
        // liquid behind it, and the interface cell itself partial. Seeding a
        // uniform ratio everywhere would be physically impossible AND would
        // leave solveUnloading's `razInter <= 0.5` and `>= 0.5` branches on one
        // side each, so half its arithmetic would never run.
        cell.razInter = index < kInterfaceCell ? 1. : 0.;
        cell.razInterIni = index < kInterfaceCell ? 1. : 0.;
        if (index == kInterfaceCell) {
            // Scenario-dependent, so the four rows differ. A fixed fraction
            // here made every scenario print the same interface arithmetic.
            // Chosen so the four scenarios STRADDLE the 0.5 threshold that
            // solveUnloading splits gas from liquid on, and so at least one
            // lands just above it: a set of ratios all on one side makes the
            // threshold itself untestable.
            cell.razInter = 0.30 + 0.45 * scenario.voidFraction;
            cell.razInterIni = 0.25 + 0.45 * scenario.voidFraction;
        }

        cell.pres *= kInjectionOverPressure;
        cell.presini *= kInjectionOverPressure;
        cell.presL *= kInjectionOverPressure;
        cell.presR *= kInjectionOverPressure;

        cell.u1R = (1.9 + scenario.liquidFlow) * (1. + 0.04 * index);
        cell.u1LL = (2.1 + scenario.liquidFlow) * (1. + 0.04 * index);
        cell.tipoCC = 0;
        cell.massfonteCH = 0.;
        cell.fonteM2 = 0.;
        cell.dTdt = 0.;
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

void printInterface(const char *method, const char *scenario, const SProd &system) {
    printf("%-28s %-15s celInter=%d velInter=%a dtInter=%a pGSup=%a presiniG=%a "
           "vazmedDesc=%a tempmedDEsc=%a\n",
           method, scenario, system.celInter, system.velInter, system.dtInter,
           system.pGSup, system.presiniG, system.vazmedDesc, system.tempmedDEsc);
}

void printUnloadingControl(const char *method, const char *scenario,
                           const SProd &system) {
    printf("%-28s %-15s presiniG=%a presMaxDesc=%a pGSup=%a\n",
           method, scenario, system.presiniG, system.arq.presMaxDesc,
           system.pGSup);
}

// printGas publishes pressures, velocities and temperatures -- none of which
// advanceInterface's hand-over touches. It moves the gas/liquid RATIO, so a row
// built from printGas alone watched the wrong fields and the hand-over could
// have been corrupted freely.
void printRatios(const char *method, const char *scenario, const CelG &cell) {
    printf("%-28s %-15s razInter=%a razInterIni=%a u1L=%a u1LL=%a\n",
           method, scenario, cell.razInter, cell.razInterIni, cell.u1L,
           cell.u1LL);
}

// Same problem, worse: everything advanceBufferedGasSubStep computes ends up in
// the choke source term and the local system it assembles, and printGas shows
// neither.
//
// What is published is the ASSEMBLY, not the solution. VGasRBuf -- the field the
// buffered update writes -- comes back NaN here, because the band system built
// from synthetic cells is singular and the elimination divides by a zero pivot.
// A NaN is deterministic and discriminates nothing, so publishing it would be a
// field that cannot fail. TL and local survive the solve as cell members, they
// are finite, and they are where GeraLocal's arithmetic and the choke source
// term actually land.
void printBuffered(const char *method, const char *scenario, const CelG &cell,
                   const CelG &inlet) {
    printf("%-28s %-15s TL0=%a TL1=%a TL2=%a inletTL0=%a massfonteCH=%a\n",
           method, scenario, cell.TL[0], cell.TL[1], cell.TL[2], inlet.TL[0],
           inlet.massfonteCH);
}

void printValve(const char *method, const char *scenario, const ChokeGas &choke) {
    printf("%-28s %-15s presEstag=%a tempEstag=%a presGarg=%a areagarg=%a\n",
           method, scenario, choke.presEstag, choke.tempEstag, choke.presGarg,
           choke.areagarg);
}

// Every call below re-seeds first. These six write through each other's state
// -- searchUnloadingInjectionPressure rewrites pGSup and every valve pressure,
// solveUnloading moves the interface -- so chaining them would make each row
// depend on the ones before it, and a single corruption would move the whole
// table at once. Independent rows localise a failure to one routine.
void runUnloadingScenario(SProd &system, Cel *cells, const Scenario &scenario) {
    auto seed = [&] {
        resetCells(system, cells, scenario);
        resetGasLine(system, scenario);
        resetValves(system);
        resetUnloading(system, scenario);
    };

    seed();
    system.HidroDescargaG();
    printGas("HidroDescargaG", scenario.name, system.celulaG[0]);
    printGas("HidroDescargaG-mid", scenario.name, system.celulaG[kInterfaceCell]);

    // computeUnloadingValvePressure returns velmax, which it sets to 0 and never
    // assigns again -- the return is a constant, in the product as much as here,
    // so a row carrying only that value can never fail. What the routine
    // actually does is move presiniG and lower arq.presMaxDesc, so those are
    // what get published.
    //
    // Two calls with throat rates on either side of vazDescControl, to drive the
    // raise branch and the lower branch rather than one of them twice.
    seed();
    system.CalcPresValvDesc(4.0 * system.arq.vazDescControl, 0);
    printUnloadingControl("CalcPresValvDesc-high", scenario.name, system);

    seed();
    system.pGSup = system.arq.presMinDesc * 0.99;
    system.CalcPresValvDesc(0.1 * system.arq.vazDescControl, 1);
    printUnloadingControl("CalcPresValvDesc-low", scenario.name, system);

    // Two rows again. The search corrects the surface pressure through the
    // valve when the tubing still carries completion fluid, and lets it decay
    // on a fixed 1% ramp when it does not -- selected by the mass content of
    // the cell above the last. Only the correcting branch was reachable with
    // resetCells' seeding, so the decay ramp was never executed at all.
    seed();
    system.celula[system.ncel - 1].MC = -1.;
    printValue("BuscaPresInjDesc-decay", scenario.name, system.BuscaPresInjDesc());
    printInterface("BuscaPresInjDesc-decay-state", scenario.name, system);

    seed();
    printValue("BuscaPresInjDesc", scenario.name, system.BuscaPresInjDesc());
    printInterface("BuscaPresInjDesc-state", scenario.name, system);
    printValve("BuscaPresInjDesc-valve", scenario.name, system.chokeVGL[0]);

    // Two rows, because advanceInterface is two routines behind one name: the
    // ordinary advance, and the hand-over that fires only when the interface is
    // in the last cell AND has filled it. Driving one leaves the other blind.
    seed();
    system.avancInter();
    printGas("avancInter", scenario.name, system.celulaG[kInterfaceCell]);

    seed();
    system.celInter = system.ncelGas - 1;
    system.celulaG[system.celInter].razInterIni = 0.995;
    system.celulaG[system.celInter].razInter = 0.995;
    system.avancInter();
    printGas("avancInter-handover", scenario.name, system.celulaG[system.ncelGas - 1]);
    printRatios("avancInter-ratios", scenario.name, system.celulaG[system.ncelGas - 1]);
    // The hand-over writes the cell AHEAD of the interface as well.
    printRatios("avancInter-ahead", scenario.name, system.celulaG[system.ncelGas]);

    seed();
    system.resolveDescarga();
    printGas("resolveDescarga", scenario.name, system.celulaG[kInterfaceCell]);
    printGas("resolveDescarga-last", scenario.name, system.celulaG[kGasCells]);

    // subtempoGasBuf is the only routine here that solves the band system, and
    // BandMtx::GaussElimPP rejects a right-hand side whose size differs from its
    // own -- through the Logger, which has no open file in this harness, so the
    // rejection arrives as a segfault rather than a message.
    //
    // The assembly loop runs 0..gasCellCount inclusive at three rows per cell,
    // so both sides must be 3 * (kGasCells + 1). termolivreG is deliberately
    // wider than that for the rest of the sweep -- resetGasLine seeds
    // 3 * (kGasCells + 2) entries and updateGasLine reads past 3 * gasCellCount
    // -- and narrowing it globally would change rows this file already
    // publishes. So it is narrowed for this call only and restored after.
    // Nothing is lost: the assembly overwrites every entry it then solves.
    seed();
    Vcr<double> wideFreeTerms = system.termolivreG;
    system.termolivreG = Vcr<double>(3 * (kGasCells + 1), 0.);
    system.subtempoGasBuf();
    printGas("subtempoGasBuf", scenario.name, system.celulaG[kInterfaceCell]);
    printBuffered("subtempoGasBuf-buffer", scenario.name,
                  system.celulaG[kInterfaceCell], system.celulaG[0]);
    // The opening bound also rewrites the injection choke itself, and that
    // write appears in no other row.
    printValve("subtempoGasBuf-choke", scenario.name, system.chokeInj);
    printInterface("subtempoGasBuf-state", scenario.name, system);
    system.termolivreG = wideFreeTerms;
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
    // VGasRBuf is the ONLY field updateBufferedGasLine writes, and printGas does
    // not carry it, so until this row existed the routine could have been
    // rewritten freely. Here the free terms come straight from resetGasLine
    // rather than from a solve, so the values are finite and discriminating.
    printf("%-28s %-15s VGasRBuf0=%a VGasRBuf2=%a VGasRBuf4=%a\n",
           "renovaGasBuf-buffer", scenario.name, system.celulaG[0].VGasRBuf,
           system.celulaG[2].VGasRBuf, system.celulaG[4].VGasRBuf);

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

    // advanceBufferedGasSubStep assembles into the band matrix, which nothing
    // else in this sweep touches, so main never sized it.
    //
    // The size is 3 * (kGasCells + 1), not 3 * (kGasCells + 3) as termolivreG
    // above: the assembly loop runs 0..gasCellCount INCLUSIVE and writes three
    // rows per cell, so it fills exactly 3 * (kGasCells + 1) of them. Sized any
    // larger, the surplus rows stay zero, the matrix is singular, and
    // GaussElimPP reports it through the Logger -- which in this harness has no
    // open file and segfaults instead. Cost of getting this wrong: a crash that
    // looks like a bug in the routine under test.
    system.matglobG = BandMtx<double>(3 * (kGasCells + 1), 5, 5);
    system.arq.gasinj.presinj = new double[2];

    for (const Scenario &scenario : kScenarios)
        runGasScenario(system, cells, scenario);
    for (const Scenario &scenario : kScenarios)
        runUnloadingScenario(system, cells, scenario);
    return 0;
}
