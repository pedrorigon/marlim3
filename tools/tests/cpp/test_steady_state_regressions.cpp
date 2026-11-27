// Steady-state regressions of the production system (SProd):
//  - marchaProdPerm2 repeats the march while the outlet mass flow changes, also
//    when that flow is negative;
//  - the choked arm of marchaInjPerm1 takes the black-oil stock-tank quality
//    when the compositional flash does not run (flashCompleto == 3);
//  - buscaProdPfundoPerm raises a low guess limited by the inlet IPR by 0.1%
//    up to ten times before giving up.
// Each scenario writes its deck from the smoothPipeFriction demo into its own
// run directory, builds the system as Marlim3's main() does and calls the
// routine directly. The control scenarios pass with and without the fixes.
//
// usage: test_steady_state_regressions <smoothPipeFriction.mr3> <run root> <scenario>
// Exit status: 0 expected behaviour, 1 otherwise, 2 usage or setup error. The
// no-root scenario expects the engine to end the run with the retry-limit
// failure (NumError calls exit()); every other scenario expects the routine to
// return.
#include "Leitura.h"
#include "Log.h"
#include "SisProd.h"
#include "estruturas.h"
#include "variaveisGlobais1D.h"

#include "rapidjson/document.h"
#include "rapidjson/prettywriter.h"
#include "rapidjson/stringbuffer.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>
#include <system_error>

// Globals defined by the engine entry point.
detTempo tempVF;
detProp prop;
detMapProp mapprop;
detCI CI;
detCC CC;
int nthrdMatriz = 1;
string versao = "test_steady_state_regressions";
int diaIni;
int horaIni;
int minutoIni;
int segundoIni;
int logRede = 0;
string nomeRedePrincipal;
ofstream arqRelatorioPerfis;
Logger logger("");
string pathArqEntrada;
string pathArqExtEntrada;
string pathPrefixoArqSaida;

namespace {

using rapidjson::Document;
using rapidjson::Value;
using Allocator = Document::AllocatorType;

constexpr int kPass = 0;
constexpr int kFail = 1;
constexpr int kSetupError = 2;

enum class Routine { March2, InjectionChoked, LowGuess };

struct Scenario {
    const char *name;
    Routine routine;
    double separatorPressure; // low-guess scenarios only
    const char *expectation;
};

const Scenario kScenarios[] = {
    {"march2_positive_flow_control", Routine::March2, 0., "positive flow, outlet at 1 kgf/cm2: 2 marches"},
    {"march2_high_outlet_control", Routine::March2, 0., "negative flow, outlet at 20 kgf/cm2: 2 marches"},
    {"march2_negative_flow", Routine::March2, 0., "negative flow, outlet at 1 kgf/cm2: 2 marches"},
    {"march2_negative_flow_hot_oil", Routine::March2, 0., "same with oil at 80 C: 2 marches"},
    {"injection_choked_quality_control", Routine::InjectionChoked, 0.,
     "stock-tank vapour fraction as read: the march returns, black-oil void fraction"},
    {"injection_choked_zero_vapour_fraction", Routine::InjectionChoked, 0.,
     "stock-tank vapour fraction 0: the march returns, black-oil void fraction"},
    {"low_guess_root_below_limit_control", Routine::LowGuess, 140., "separator 140: converged pressure"},
    {"low_guess_root_near_reservoir_pressure", Routine::LowGuess, 149.9, "separator 149.9: converged pressure"},
    {"low_guess_no_root", Routine::LowGuess, 160., "separator 160: retry-limit failure after ten retries"},
};

varGlob1D lineGlobals;
SProd *sp = nullptr;
std::string logPath;

// Failure that ends buscaProdPfundoPerm once its retries are exhausted.
constexpr const char *kRetryLimitFailure =
    "Busca de valores iniciais para calculo de zero de funcao em buscaProdPfundoPerm atingiu maximo de iteracoes";

[[noreturn]] void setupError(const std::string &what) {
    std::printf("SETUP ERROR %s\n", what.c_str());
    std::fflush(stdout);
    std::_Exit(kSetupError);
}

// Verdict reported when the engine ends the process (NumError calls exit()).
int (*verdictOnEngineExit)() = nullptr;

void reportEngineExit() {
    if (verdictOnEngineExit != nullptr) {
        const int status = verdictOnEngineExit();
        std::fflush(stdout);
        std::_Exit(status);
    }
}

int unexpectedEngineExit() {
    std::printf("SETUP ERROR the engine ended the process outside the call under test\n");
    return kSetupError;
}

// Description of the last failure in the simulation log the engine wrote before ending the process, or "".
std::string lastLoggedFailure() {
    std::ifstream in(logPath, std::ios::binary);
    if (!in.is_open())
        return "";
    const std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    rapidjson::Document log;
    log.Parse(text.c_str());
    if (log.HasParseError() || !log.IsObject() || !log.HasMember("resultadoSimulacao"))
        return "";
    const rapidjson::Value &result = log["resultadoSimulacao"];
    if (!result.IsObject() || !result.HasMember("logs") || !result["logs"].IsArray())
        return "";
    std::string last;
    for (const rapidjson::Value &entry : result["logs"].GetArray())
        if (entry.IsObject() && entry.HasMember("log") && entry["log"].IsString() &&
            std::string(entry["log"].GetString()) == "FALHA" && entry.HasMember("descricao") &&
            entry["descricao"].IsString())
            last = entry["descricao"].GetString();
    return last;
}

// ---------------------------------------------------------------- decks

Value &member(Value &object, const char *key) {
    if (!object.IsObject() || !object.HasMember(key))
        setupError(std::string("demo deck without ") + key);
    return object[key];
}

Value &firstItem(Value &object, const char *key) {
    Value &array = member(object, key);
    if (!array.IsArray() || array.Empty())
        setupError(std::string("demo deck without items in ") + key);
    return array[0];
}

void put(Value &object, const char *key, Value &&value, Allocator &alloc) {
    if (object.HasMember(key))
        object[key] = value;
    else
        object.AddMember(Value(key, alloc), value, alloc);
}

template <class T>
Value single(T item, Allocator &alloc) {
    Value array(rapidjson::kArrayType);
    array.PushBack(item, alloc);
    return array;
}

void readDemo(const std::string &path, Document &deck) {
    std::ifstream in(path, std::ios::binary);
    if (!in.is_open())
        setupError("cannot open " + path);
    const std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    deck.Parse(text.c_str());
    if (deck.HasParseError() || !deck.IsObject())
        setupError("cannot parse " + path);
}

// The demo is a horizontal line of 50 cells of 100 m; all decks use the default
// drift-flux model instead of its black-box friction correlation.
void editDeck(Document &deck, const Scenario &scenario) {
    Allocator &a = deck.GetAllocator();
    member(deck, "initialConfig").RemoveMember("driftModel");
    firstItem(deck, "productionPipe").RemoveMember("ssFlowModel");
    switch (scenario.routine) {
    case Routine::March2: {
        // surface choke at 0.5: a production run takes buscaProdPfundoPerm2/marchaProdPerm2
        Value choke(rapidjson::kObjectType);
        choke.AddMember("dischargeCoefficient", 0.84, a);
        choke.AddMember("model", 0, a);
        choke.AddMember("opening", single(0.5, a), a);
        choke.AddMember("time", single(0, a), a);
        put(deck, "surfaceChoke", std::move(choke), a);
        // keeps the friction loss below 0.1 kgf/cm2, so the outlet can be at 1 kgf/cm2 with the flow reversed
        put(firstItem(deck, "liquidSource"), "liquidFlowRate", single(500.0, a), a);
        break;
    }
    case Routine::InjectionChoked: {
        put(firstItem(deck, "productionFluid"), "gor", Value(100), a);
        Value &gas = firstItem(deck, "gasSource");
        put(gas, "measuredLength", Value(0.1), a);
        put(gas, "active", Value(true), a);
        put(gas, "dry", Value(false), a);
        put(gas, "prodFluidId", Value(0), a);
        put(firstItem(deck, "liquidSource"), "active", Value(false), a);
        break;
    }
    case Routine::LowGuess: {
        put(firstItem(deck, "productionFluid"), "gor", Value(50), a);
        Value &ipr = firstItem(deck, "ipr");
        put(ipr, "active", Value(true), a);
        put(ipr, "staticPressure", single(150, a), a);
        put(ipr, "ip", single(50, a), a);
        put(ipr, "ii", single(50, a), a);
        put(firstItem(deck, "liquidSource"), "active", Value(false), a);
        put(member(deck, "separator"), "pressure", single(scenario.separatorPressure, a), a);
        break;
    }
    }
}

void writeDeck(const Document &deck, const std::string &path) {
    rapidjson::StringBuffer buffer;
    rapidjson::PrettyWriter<rapidjson::StringBuffer> writer(buffer);
    deck.Accept(writer);
    std::ofstream out(path, std::ios::binary);
    out << buffer.GetString() << '\n';
    if (!out)
        setupError("cannot write " + path);
}

// Same steps as main() for a lone line: logger, paths, profile report, SProd.
SProd *buildSystem(const std::string &deckPath, const std::string &logPath, const std::string &dir) {
    const time_t now = time(nullptr);
    const tm *local = localtime(&now);
    diaIni = local->tm_mday;
    horaIni = local->tm_hour;
    minutoIni = local->tm_min;
    segundoIni = local->tm_sec;
    logger.changeResultadoSimulacao(deckPath, logPath);
    pathArqEntrada = dir;
    pathArqExtEntrada = dir;
    pathPrefixoArqSaida = dir;
    arqRelatorioPerfis.open(dir + "resultado.log", ios_base::out);
    SProd *system = new SProd(deckPath, logPath, tipoValidacaoJson_t::none, tipoSimulacao_t::transiente, &lineGlobals);
    logger.setNomeArqEntrada(deckPath);
    system->modoPerm = 1; // as SolveTramoSolteiro
    return system;
}

// ---------------------------------------------------------------- marchaProdPerm2

// With accelerateSteadyConvergence off the march repeats at most once, and the
// second while test compares the first march with the initial values
// masfim = presteste = 1. With the first-march outlet at 1 kgf/cm2 only the flow
// test can ask for the second march.

// A criterion of 5 lets the first while test pass (ratios of 11 against the
// initial -10) and makes the second fail for |masfim| of a few kg/s, with or
// without the fix: exactly one march per call.
double outletAfterOneMarch(double inlet) {
    const double criterion = sp->arq.CriterioConvergPerm;
    sp->arq.CriterioConvergPerm = 5.;
    const double residual = sp->marchaProdPerm2(inlet);
    sp->arq.CriterioConvergPerm = criterion;
    if (residual <= -0.9e10)
        return -1e300; // the pressure ran out: below any target
    if (residual >= 0.9e10 || sp->iterperm != 1)
        setupError("tuning march failed");
    return sp->celula[sp->ncel].pres;
}

// Bisection on the inlet pressure until the first-march outlet is within
// 1e-7 kgf/cm2 of the target (the outlet grows with the inlet pressure).
double tuneInlet(double target) {
    double low = 0.2;
    double high = target + 1.;
    for (int k = 0; k < 20 && outletAfterOneMarch(high) < target; k++)
        high *= 2.;
    double inlet = high;
    double error = outletAfterOneMarch(inlet) - target;
    for (int k = 0; k < 100 && std::fabs(error) > 1e-7; k++) {
        inlet = 0.5 * (low + high);
        error = outletAfterOneMarch(inlet) - target;
        (error < 0. ? low : high) = inlet;
    }
    if (std::fabs(error) > 1e-7)
        setupError("inlet pressure tuning did not converge");
    return inlet;
}

bool march2(const std::string &name) {
    InjLiq &source = sp->celula[0].acsr.injl;
    if (sp->celula[0].acsr.tipo != 2 || sp->arq.AceleraConvergPerm != 0)
        setupError("expected an inlet liquid source and accelerateSteadyConvergence off");
    const double rate = std::fabs(source.QLiq);
    double target = 1.;
    source.QLiq = -rate;
    if (name == "march2_positive_flow_control") {
        source.QLiq = rate;
    } else if (name == "march2_high_outlet_control") {
        target = 20.;
    } else if (name == "march2_negative_flow_hot_oil") {
        // the second march starts from the boundary temperatures of the first
        source.bet = 0.;
        source.temp = 80.;
    }
    const double inlet = tuneInlet(target);
    const double masfim = sp->celula[sp->ncel - 1].MC;
    std::printf("first march: inlet %.12g, outlet %.12g kgf/cm2, outlet mass flow %.12g kg/s\n", inlet,
                sp->celula[sp->ncel].pres, masfim);
    const double residual = sp->marchaProdPerm2(inlet);
    std::printf("marchaProdPerm2, criterion %g: %d march(es), outlet %.12g kgf/cm2, residual %.12g\n",
                sp->arq.CriterioConvergPerm, sp->iterperm, sp->celula[sp->ncel].pres, residual);
    return sp->iterperm == 2;
}

// ---------------------------------------------------------------- marchaInjPerm1

double expectedVoidFraction = 0.;

// Choked arm with the black-oil quality: the total mass cancels in the ratio.
double blackOilVoidFraction(double rate) {
    InjGas &gas = sp->celula[0].acsr.injg;
    const double saved = gas.QGas;
    gas.QGas = rate;
    const double massGas = gas.VMas(sp->pGSup, gas.temp);
    gas.QGas = saved;
    double quality = gas.FluidoPro.FracMassHidra(1., 20.);
    const double massTotal = massGas / quality;
    quality = gas.FluidoPro.FracMassHidra(sp->pGSup, gas.temp);
    const double qGas = massTotal * quality / gas.FluidoPro.MasEspGas(sp->pGSup, gas.temp);
    const double qLiq = massTotal * (1. - quality) / gas.FluidoPro.MasEspLiq(sp->pGSup, gas.temp);
    return qGas / (qLiq + qGas);
}

bool voidFractionMatches() {
    const double got = sp->celula[0].alfini;
    std::printf("void fraction of the first cell %.17g (expected %.17g)\n", got, expectedVoidFraction);
    return std::isfinite(got) && std::fabs(got - expectedVoidFraction) <= 1e-12 * std::fabs(expectedVoidFraction);
}

// Both injection scenarios expect marchaInjPerm1 to return: ending the process fails them.
int injectionVerdictOnExit() {
    const std::string failure = lastLoggedFailure();
    std::printf("the engine ended the process during marchaInjPerm1: %s\n",
                failure.empty() ? "no failure logged" : failure.c_str());
    voidFractionMatches();
    std::printf("expected marchaInjPerm1 to return\nFAIL\n");
    return kFail;
}

bool injectionChoked(const std::string &name) {
    InjGas &gas = sp->celula[0].acsr.injg;
    if (sp->celula[0].acsr.tipo != 1 || gas.seco != 0)
        setupError("expected a non-dry gas source in the first cell");
    // flags of the black-oil pass that precedes a compositional injector solution
    sp->arq.flashCompleto = 3;
    sp->CalcLat = 0;
    sp->arq.tipoFluido = 0;
    for (int i = 0; i <= sp->ncel; i++)
        sp->celula[i].flui.flashCompleto = 3;
    gas.FluidoPro.flashCompleto = 3;
    sp->arq.condpocinj.CC = 1;                          // injected gas rate
    sp->chokeSup.AreaGarg = 0.3 * sp->chokeSup.AreaTub; // choked arm
    sp->pGSup = 100.;                                   // pressure upstream of the choke
    std::printf("stock-tank vapour fraction as read %.17g, black-oil FracMassHidra(1, 20) %.17g\n",
                gas.FluidoPro.dStockTankVaporMassFraction, gas.FluidoPro.FracMassHidra(1., 20.));
    if (name == "injection_choked_zero_vapour_fraction")
        gas.FluidoPro.dStockTankVaporMassFraction = 0.;
    const double rate = gas.QGas;
    expectedVoidFraction = blackOilVoidFraction(rate);
    verdictOnEngineExit = injectionVerdictOnExit;
    const double residual = sp->marchaInjPerm1(rate);
    verdictOnEngineExit = unexpectedEngineExit;
    // the march returns the outlet mass flow, or +-1e10 when it fails
    const bool returned = std::isfinite(residual) && std::fabs(residual) < 0.9e10;
    std::printf("marchaInjPerm1(%g) returned %.12g%s\n", rate, residual, returned ? "" : ", a failure value");
    const bool matches = voidFractionMatches();
    return returned && matches;
}

// ---------------------------------------------------------------- buscaProdPfundoPerm

double guessAfterTenRetries = 0.;
bool expectNoRoot = false;

// Without a root the search must end the run with the retry-limit failure after raising the limited guess ten times;
// the same message ends the search at other points, so the last guess identifies the retries.
int lowGuessVerdictOnExit() {
    const double last = sp->celula[0].pres;
    const std::string failure = lastLoggedFailure();
    std::printf("the engine ended the process during buscaProdPfundoPerm: %s\n",
                failure.empty() ? "no failure logged" : failure.c_str());
    std::printf("last guess %.17g (after ten retries %.17g)\n", last, guessAfterTenRetries);
    const bool retried = std::fabs(last - guessAfterTenRetries) <= 1e-9;
    const bool diagnosed = failure.find(kRetryLimitFailure) != std::string::npos;
    if (expectNoRoot)
        std::printf("ten retries: %s; retry-limit failure logged: %s\n", retried ? "yes" : "no", diagnosed ? "yes" : "no");
    else
        std::printf("expected buscaProdPfundoPerm to return\n");
    const bool pass = expectNoRoot && retried && diagnosed;
    std::printf("%s\n", pass ? "PASS" : "FAIL");
    return pass ? kPass : kFail;
}

bool lowGuess(const std::string &name) {
    const IPR &ipr = sp->celula[0].acsr.ipr;
    // permanenteSimples calls buscaProdPfundoPerm(chutePerm) for this deck
    if (sp->arq.pocinjec != 0 || sp->arq.ConContEntrada != 0 || sp->celula[0].acsr.tipo != 3 ||
        !(sp->arq.chokep.abertura[0] > 0.6) || sp->arq.chokep.curvaDinamic != 0 || sp->arq.lingas != 0 ||
        sp->arq.chutePerm >= 0.)
        setupError("unexpected deck for buscaProdPfundoPerm");
    // guess limited by the IPR when chute < 0, and that guess after ten retries
    double limited = -(10 / ipr.ip) + ipr.Pres;
    if (limited < 0.999 * ipr.Pres)
        limited = 0.999 * ipr.Pres;
    guessAfterTenRetries = limited;
    for (int k = 0; k < 10; k++)
        guessAfterTenRetries *= 1.001;
    expectNoRoot = name == "low_guess_no_root";
    std::printf("reservoir pressure %g, separator %g, limited guess %.17g\n", ipr.Pres, sp->pGSup, limited);
    verdictOnEngineExit = lowGuessVerdictOnExit;
    const double pressure = sp->buscaProdPfundoPerm(sp->arq.chutePerm);
    verdictOnEngineExit = unexpectedEngineExit;
    const double last = sp->celula[0].pres;
    std::printf("buscaProdPfundoPerm returned %.17g, last guess %.17g\n", pressure, last);
    if (expectNoRoot) {
        std::printf("expected the run to end with the retry-limit failure\n");
        return false;
    }
    const bool found = std::fabs(pressure) < 0.9e10 && pressure > 0. && pressure < ipr.Pres;
    const double residual = found ? sp->marchaProdPerm1(pressure) : 1e10;
    std::printf("search residual at the returned pressure %.6g kgf/cm2\n", residual);
    return found && std::fabs(residual) <= 1e-4 * sp->pGSup;
}

} // namespace

int main(int argc, char **argv) {
    std::setvbuf(stdout, nullptr, _IONBF, 0);
    if (argc != 4) {
        std::printf("usage: %s <smoothPipeFriction.mr3> <run root> <scenario>\n", argv[0]);
        return kSetupError;
    }
    const std::string name = argv[3];
    const Scenario *scenario = nullptr;
    for (const Scenario &candidate : kScenarios)
        if (name == candidate.name)
            scenario = &candidate;
    if (scenario == nullptr)
        setupError("unknown scenario " + name);
    std::printf("scenario %s: %s\n", scenario->name, scenario->expectation);

    Document deck;
    readDemo(argv[1], deck);
    editDeck(deck, *scenario);
    const std::filesystem::path runDir = std::filesystem::absolute(std::filesystem::path(argv[2]) / name);
    std::error_code error;
    std::filesystem::create_directories(runDir, error);
    if (!error)
        std::filesystem::current_path(runDir, error); // for files the engine writes to the working directory
    if (error)
        setupError("cannot use the run directory " + runDir.string());
    const std::string dir = runDir.string() + "/";
    writeDeck(deck, dir + name + ".mr3");

    std::atexit(reportEngineExit);
    verdictOnEngineExit = unexpectedEngineExit;
    logPath = dir + name + ".log";
    sp = buildSystem(dir + name + ".mr3", logPath, dir);

    bool pass = false;
    switch (scenario->routine) {
    case Routine::March2:
        pass = march2(name);
        break;
    case Routine::InjectionChoked:
        pass = injectionChoked(name);
        break;
    case Routine::LowGuess:
        pass = lowGuess(name);
        break;
    }
    verdictOnEngineExit = nullptr;
    std::printf("%s\n", pass ? "PASS" : "FAIL");
    std::fflush(stdout);
    std::_Exit(pass ? kPass : kFail); // skips the engine's static destructors
}
