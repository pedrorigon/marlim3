// The runner inserts actual production bodies into this controlled-dependency fixture.
#include "SisProdConstants.h"
#include "UnitConversions.h"
#include <array>
#include <bit>
#include <cfenv>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <vector>

using std::fabs;
using std::isnan;
using namespace sisprod;

namespace {
std::vector<double> events;
void event(double id, double p = 0., double t = 0.) { events.insert(events.end(), {id, p, t}); }
void require(bool condition) { if (!condition) throw std::runtime_error("steady march contract failed"); }
}

struct ProFlu {
    double id = 1., BSW = .2, RGO = 100., rDgD = 1., Deng = .7;
    double dStockTankVaporMassFraction = .3, flashed = 0.;
    void atualizaPropComp(double p, double t, double beta = -1., double * = nullptr,
                          double * = nullptr, int gasMode = 0) {
        event(id * 100. + 1., p, t); event(gasMode, beta);
        flashed = .001 * p + .0001 * t;
    }
    double FracMassHidra(double p, double t) const { event(id * 100. + 2., p, t); return .25 + flashed; }
    double MasEspGas(double p, double t) const { event(id * 100. + 3., p, t); return 8. + .01 * p + flashed; }
    double MasEspLiq(double p, double t) const { event(id * 100. + 4., p, t); return 800. + .01 * t + flashed; }
    double RS(double p, double t) const { event(id * 100. + 5., p, t); return .01 * p; }
    double BOFunc(double p, double t) const { event(id * 100. + 6., p, t); return 1.1 + flashed; }
    double BAFunc(double p, double t) const { event(id * 100. + 7., p, t); return 1.01; }
};
struct ComplementaryFluid {
    double standardDensity() const { event(900.); return 1100.; }
    double MasEspFlu(double p, double t) const { event(901., p, t); return 1100. + .01 * p; }
};
struct GasSource {
    double temp = 61., QGas = 300., razCompGas = .1;
    int seco = 0;
    ProFlu FluidoPro{2.};
    ComplementaryFluid fluidocol;
    double VMas(double p, double t) { event(1000., p, t); return QGas * .01; }
};
struct LiquidSource {
    double temp = 62., QLiq = 300., bet = .2;
    ProFlu FluidoPro{3.};
};
struct InflowSource {
    double Tres = 63., Pres = 300., deriG = -1., deriP = -2., deriC = -3., Psat = 50.;
    ProFlu FluidoPro{4.};
    double MasG(double p, double t) { event(1100., p, t); deriG = p * .01; deriP = -p; Psat += .125; return 3.; }
    double MasL(double p, double t) { event(1101., p, t); deriG = -p; deriP = p * .02; deriC = 0.; Psat += .125; return 8.; }
};
struct PorousData {
    double tRes = 64., fluxIniG = 2., fluxIni = 3., fluxIniA = 4., BSW = .3;
    ProFlu flup{5.};
};
struct Porous2DData {
    double tRes = 65.;
    ProFlu flup{6.};
    PorousData transfer;
};
struct MultipleSource {
    double temp = 66., MassG = 3., MassP = 8., MassC = 2.;
    int condTermo = 0;
    ProFlu FluidoPro{7.};
    ComplementaryFluid fluidocol;
};
struct Accessory {
    int tipo = kAccessoryNone;
    GasSource injg;
    LiquidSource injl;
    InflowSource ipr;
    PorousData radialPoro;
    struct { Porous2DData dados; } poroso2D;
    MultipleSource injm;
};
struct Cell {
    double temp = 70., tempL = 10., tempini = 11.;
    double ML = 2., MC = 3., MliqiniL = 4., Mliqini = 5., MComp = 6., QLL = 7., QL = 8., QG = 9.;
    double presauxL = 10., presLini = 11., presL = 12., pres = 13., presini = 14., presaux = 15.;
    double alf = .1, alfini = .2, bet = .3, betini = .4, alfL = .5, alfLini = .6;
    double alfPigD = .1, alfPigDini = .2, alfPigE = .3, alfPigEini = .4;
    double betL = .5, betLini = .6, betPigD = .7, betPigDini = .8, betPigE = .9, betPigEini = .1, betI = .2, betLI = .3;
    ProFlu flui;
    Accessory acsr;
};
struct Input {
    int flashCompleto = 0, tabelaDinamica = 0, AceleraConvergPerm = 0, lingas = 0, nvalvgas = 0;
    int usaTabela = 0, AP = 0, corrigeContSep = 0;
    double CriterioConvergPerm = .001;
    struct { double pmin = .2; } tabent;
    struct External { double textern = 60.; };
    std::array<External, 1> celp;
};
struct Globals { int networkMode = 0, parallelNetworkMode = 0; double smallValueTolerance = 1e-9; };
struct GasCell { int tipoCC = 0; };

namespace sisprod::steady {
struct SteadyStateState {
    Cell *cells;
    Input &input;
    Globals *globals;
    GasCell *gasCells;
    int lastCell = 3, endNode = 1;
    bool networkCoupled = true;
    mutable int steadyIteration = 99;
    double convergenceMonitor = .5, gasSurfacePressure = 20.;
    mutable double baseConvergenceMonitor = -7.;
};
struct SteadyDependencies { int sources = 0, gasLift = 0; };
}

namespace {
using namespace sisprod::steady;
SteadyStateState *active = nullptr;
int attempts = 0, behavior = 0;
double seedAlpha = .2, seedBeta = .3;
double expectedAlpha = .2, expectedBeta = .3;
void NumError(const char *) { throw std::runtime_error("unexpected iteration-limit error"); }
}
namespace sisprod::sources {
void renewSourceTerms(int, int i) {
    event(1200., i);
    auto &cell = active->cells[i];
    cell.acsr.radialPoro.fluxIniG += .25;
    cell.acsr.radialPoro.BSW += .125;
    cell.acsr.poroso2D.dados.transfer.fluxIniA += .5;
    cell.acsr.poroso2D.dados.transfer.BSW += .125;
    cell.acsr.injm.MassG += .25;
    cell.acsr.injm.MassP -= .25;
}
}
namespace sisprod::gaslift {
void initializeTubingConnectionSteady(int) { event(1300.); }
void connectTubingSteady(int) { event(1301.); }
void initializeSteadyValveGasFlowRate(int, int i) { event(1302., i); }
}

// @PRIVATE_HEADER@

namespace sisprod::steady {
namespace {
#ifdef MARLIM_REFERENCE_BUILD
struct InitialPhaseFractions { double voidFraction, complementaryLiquidFraction; };
void seedFirstCellVoidFraction(const SteadyStateState &, const SteadyDependencies &, double p, double &a, double &b) {
    event(1400., p); a = seedAlpha; b = seedBeta;
}
#else
InitialPhaseFractions seedFirstCellVoidFraction(const SteadyStateState &, const SteadyDependencies &, double p) {
    event(1400., p); return {seedAlpha, seedBeta};
}
#endif
template <typename Condition>
void estimateValveGasPressure(const SteadyStateState &state, Condition) {
    event(1401., std::is_same_v<Condition, InjectionPressureCondition>);
    state.cells[0].temp += .125;
}
bool advanceProductionCells(const SteadyStateState &state, const SteadyDependencies &, double p, int &i, double &abort) {
    event(1402., p, i);
    require(state.cells[0].pres == p);
    require(std::bit_cast<std::uint64_t>(state.cells[0].alf) == std::bit_cast<std::uint64_t>(expectedAlpha));
    require(std::bit_cast<std::uint64_t>(state.cells[0].bet) == std::bit_cast<std::uint64_t>(expectedBeta));
    ++attempts;
    if (behavior == 1) { abort = -kSteadyMarchLimit; return true; }
    if (behavior == 2 && attempts % 2 == 1) {
        state.cells[0].pres += 1.; state.cells[0].alf = -.5; state.cells[0].bet = -.5; i = 1; return false;
    }
    state.cells[state.lastCell - 1].MC = behavior == 3 ? 0. : attempts * 2.;
    state.cells[state.lastCell].pres = 10. + attempts;
    i = state.lastCell + 1;
    return false;
}
void marchGasLineAndCoupleAnnulus(const SteadyStateState &, const SteadyDependencies &, double p) { event(1403., p); }
double steadyPressureAtLastCell(const SteadyStateState &) { event(1404.); return .125; }

// @HELPERS@

}  // namespace

// @MARCH@

}  // namespace sisprod::steady

namespace {
void record(std::FILE *out, double value) { require(std::fwrite(&value, sizeof(value), 1, out) == 1); }
void stateRecord(std::FILE *out, const SteadyStateState &state, double a, double b) {
    record(out, a); record(out, b); record(out, state.steadyIteration); record(out, state.baseConvergenceMonitor);
    for (int i = 0; i <= state.lastCell; ++i) {
        const auto &c = state.cells[i];
        for (double v : {c.temp,c.tempL,c.tempini,c.ML,c.MC,c.MliqiniL,c.Mliqini,c.MComp,c.QLL,c.QL,c.QG,
                         c.presauxL,c.presLini,c.presL,c.pres,c.presini,c.presaux,c.alf,c.alfini,c.bet,c.betini,
                         c.alfL,c.alfLini,c.alfPigD,c.alfPigDini,c.alfPigE,c.alfPigEini,c.betL,c.betLini,
                         c.betPigD,c.betPigDini,c.betPigE,c.betPigEini,c.betI,c.betLI,c.flui.BSW,c.flui.flashed,
                         c.acsr.ipr.deriG,c.acsr.ipr.deriP,c.acsr.ipr.deriC,c.acsr.ipr.Psat,
                         c.acsr.radialPoro.fluxIniG,c.acsr.radialPoro.BSW,c.acsr.poroso2D.dados.transfer.BSW,
                         c.acsr.injm.MassG,c.acsr.injm.MassP,c.acsr.injg.FluidoPro.flashed,
                         c.acsr.injl.FluidoPro.flashed,c.acsr.ipr.FluidoPro.flashed,c.acsr.radialPoro.flup.flashed,
                         c.acsr.poroso2D.dados.flup.flashed,c.acsr.injm.FluidoPro.flashed}) record(out,v);
    }
    record(out,events.size()); for (double e : events) record(out,e);
}
}

int main(int argc, char **argv) {
    require(argc == 2);
    auto *out = std::fopen(argv[1], "wb"); require(out != nullptr);
    int seedCases = 0, marchCases = 0;
    for (int kind : {0,1,2,3,10,15,16,4,-1,99})
        for (int flash : {0,1,2}) for (int table : {0,1})
            for (int dry : {0,1,2}) for (int condition : {0,1}) {
                std::array<Cell,4> cells{}; Input input; Globals globals; GasCell gas;
                SteadyStateState state{cells.data(),input,&globals,&gas}; active = &state;
                input.flashCompleto=flash; input.tabelaDinamica=table;
                cells[0].acsr.tipo=kind; cells[0].acsr.injg.seco=dry; cells[0].acsr.injm.condTermo=condition;
                events.clear(); double a=0.,b=0.;
#ifdef MARLIM_REFERENCE_BUILD
                observedSeed(state,{},65.,a,b);
#else
                auto fractions=observedSeed(state,{},65.); a=fractions.voidFraction; b=fractions.complementaryLiquidFraction;
#endif
                stateRecord(out,state,a,b); ++seedCases;
            }
    for (int accelerated : {0,1,2}) for (double criterion : {0.,.001,11.,12.,std::numeric_limits<double>::quiet_NaN()})
        for (int mode : {0,1,2}) for (int gasMode : {0,1,2}) for (int scenario : {0,1,2,3,4,5,6}) {
            std::array<Cell,4> cells{}; Input input; Globals globals; GasCell gas;
            SteadyStateState state{cells.data(),input,&globals,&gas}; active=&state;
            input.AceleraConvergPerm=accelerated; input.CriterioConvergPerm=criterion;
            input.lingas=gasMode!=0; input.nvalvgas=gasMode!=0 ? 2 : 0; gas.tipoCC=gasMode-1;
            globals.networkMode=mode; input.corrigeContSep=mode%2;
            seedAlpha=expectedAlpha=.2; seedBeta=expectedBeta=.3; behavior=scenario; attempts=0;
            double p=65.;
            if (scenario==4) p=.1;
            if (scenario==5) {input.usaTabela=1; input.tabent.pmin=66.;}
            if (scenario==6) {cells[0].acsr.tipo=kAccessoryInflowPerformance; cells[0].acsr.ipr.Pres=p;}
            events.clear(); std::feclearexcept(FE_ALL_EXCEPT);
            double residual=marchProductionSteady(state,{},p);
            const int exceptions=std::fetestexcept(FE_ALL_EXCEPT);
            require(state.steadyIteration>=0 && state.steadyIteration<=(accelerated==1 ? 1 : 2));
            stateRecord(out,state,residual,attempts); record(out,exceptions); ++marchCases;
        }
    const double cutoff=1e-6, below=std::nextafter(cutoff,0.), above=std::nextafter(cutoff,1.);
    const double nan=std::numeric_limits<double>::quiet_NaN(), infinity=std::numeric_limits<double>::infinity();
    for (const auto &values : std::array<std::array<double,4>,9>{{
             {below,.4,0.,.4},{cutoff,.4,cutoff,.4},{above,.4,above,.4},
             {.2,-below,.2,0.},{.2,-cutoff,.2,-cutoff},{-above,below,-above,0.},
             {-0.,0.,0.,0.},{nan,.3,nan,.3},{.2,infinity,.2,infinity}}}) {
        std::array<Cell,4> cells{}; Input input; Globals globals; GasCell gas;
        SteadyStateState state{cells.data(),input,&globals,&gas}; active=&state;
        input.AceleraConvergPerm=1;
        seedAlpha=values[0];seedBeta=values[1];expectedAlpha=values[2];expectedBeta=values[3];
        behavior=0;attempts=0;events.clear();
        const double residual=marchProductionSteady(state,{},65.);
        require(attempts==1 && state.steadyIteration==1);
        stateRecord(out,state,residual,attempts);++marchCases;
    }
    require(std::fclose(out)==0);
    std::printf("PASS seed cases=%d, march cases=%d; production bodies with controlled dependencies\n",seedCases,marchCases);
}
