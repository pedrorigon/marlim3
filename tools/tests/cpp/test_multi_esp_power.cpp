// Regression checks of the multi-stage ESP (multiBomCentSub) on the public demos/MultiESP deck: marches repeated with
// identical inputs must give the same pressure increment, power and heat to the cell, without motor heat
// (motorHeatFraction 0, the demo's value) and with it (motorHeatFraction 0.5).
// Each scenario copies the deck and its PVTSim table into its own run directory, builds the system as Marlim3's
// main() does and calls the march of the pump directly.
//
// usage: test_multi_esp_power <demos/MultiESP directory> <run root> <scenario>
// Exit status: 0 expected behaviour, 1 otherwise, 2 usage or setup error (also when the engine ends the process).
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
string versao = "test_multi_esp_power";
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

constexpr int kPass = 0;
constexpr int kFail = 1;
constexpr int kSetupError = 2;

struct Scenario {
    const char *name;
    double motorHeatFraction; // < 0: as in the demo
    const char *expectation;
};

const Scenario kScenarios[] = {
    {"repeated_marches_without_motor_heat", -1., "motorHeatFraction as in the demo (0): every march gives the same results"},
    {"repeated_marches_with_motor_heat", 0.5, "motorHeatFraction 0.5: every march gives the same results"},
};

// Inlet of the pump at the steady state of demos/MultiESP (resumoPermanente.dat and the 2000 m trend at 0 s):
// in-situ gas and liquid flow at the suction (m3/s), pressure (kgf/cm2), temperature (C), void and complementary
// liquid fractions.
constexpr double kGasFlow = 0.0253087;
constexpr double kLiquidFlow = 0.0154487;
constexpr double kPressure = 41.4434;
constexpr double kTemperature = 49.2072;
constexpr double kVoidFraction = 0.510551;
constexpr double kComplementaryFraction = 0.;
constexpr int kMarches = 3;

varGlob1D lineGlobals;

bool unexpectedExitArmed = false;

void reportEngineExit() {
    if (unexpectedExitArmed) {
        std::printf("SETUP ERROR the engine ended the process\n");
        std::fflush(stdout);
        std::_Exit(kSetupError);
    }
}

[[noreturn]] void setupError(const std::string &what) {
    std::printf("SETUP ERROR %s\n", what.c_str());
    std::fflush(stdout);
    std::_Exit(kSetupError);
}

void writeDeck(const std::filesystem::path &demo, const std::filesystem::path &dir, double motorHeatFraction) {
    std::ifstream in(demo / "MultiESP.mr3", std::ios::binary);
    if (!in.is_open())
        setupError("cannot open " + (demo / "MultiESP.mr3").string());
    const std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    rapidjson::Document deck;
    deck.Parse(text.c_str());
    if (deck.HasParseError() || !deck.IsObject() || !deck.HasMember("multiESP") || !deck["multiESP"].IsArray() ||
        deck["multiESP"].Empty())
        setupError("demo deck without multiESP");
    if (motorHeatFraction >= 0.) {
        rapidjson::Value &pump = deck["multiESP"][0];
        if (pump.HasMember("motorHeatFraction"))
            pump["motorHeatFraction"] = motorHeatFraction;
        else
            pump.AddMember("motorHeatFraction", motorHeatFraction, deck.GetAllocator());
    }
    rapidjson::StringBuffer buffer;
    rapidjson::PrettyWriter<rapidjson::StringBuffer> writer(buffer);
    deck.Accept(writer);
    std::ofstream out(dir / "MultiESP.mr3", std::ios::binary);
    out << buffer.GetString() << '\n';
    if (!out)
        setupError("cannot write the deck");
    std::error_code error;
    std::filesystem::copy_file(demo / "PVTSIM-MARLIM.tab", dir / "PVTSIM-MARLIM.tab",
                               std::filesystem::copy_options::overwrite_existing, error);
    if (error)
        setupError("cannot copy the PVTSim table");
}

// Same steps as main() for a lone line: logger, paths, profile report, SProd.
SProd *buildSystem(const std::string &dir) {
    const time_t now = time(nullptr);
    const tm *local = localtime(&now);
    diaIni = local->tm_mday;
    horaIni = local->tm_hour;
    minutoIni = local->tm_min;
    segundoIni = local->tm_sec;
    const std::string deck = dir + "MultiESP.mr3";
    const std::string log = dir + "simulacao.log";
    logger.changeResultadoSimulacao(deck, log);
    pathArqEntrada = dir;
    pathArqExtEntrada = dir;
    pathPrefixoArqSaida = dir;
    arqRelatorioPerfis.open(dir + "resultado.log", ios_base::out);
    SProd *system = new SProd(deck, log, tipoValidacaoJson_t::none, tipoSimulacao_t::transiente, &lineGlobals);
    logger.setNomeArqEntrada(deck);
    return system;
}

int pumpCell(const SProd &system) {
    for (int i = 0; i <= system.ncel; i++)
        if (system.celula[i].acsr.tipo == 17)
            return i;
    return -1;
}

struct MarchResult {
    double dpB, potB, potBT, potTermo, heat;
};

// The fluids and the heat given to the cell as Cel::GeraLocal and the steady-state march do.
MarchResult march(Cel &cell) {
    multiBomCentSub &pump = cell.acsr.multibcs;
    pump.flui = cell.flui;
    pump.fluicol = cell.fluicol;
    pump.marchaMultiBcs(kGasFlow, kLiquidFlow, kPressure, kTemperature, kVoidFraction, kComplementaryFraction);
    const double heat = pump.potTermo + pump.potBT * (1. - pump.eficM / 100.) * pump.fracTermMotorEfic;
    return {pump.dpB, pump.potB, pump.potBT, pump.potTermo, heat};
}

bool same(double a, double b) {
    return std::isfinite(a) && std::isfinite(b) && std::fabs(a - b) <= 1e-12 * std::fmax(std::fabs(a), std::fabs(b));
}

bool repeatedMarches(Cel &cell, const Scenario &scenario) {
    const multiBomCentSub &pump = cell.acsr.multibcs;
    std::printf("%d internal pumps, motor efficiency %g%%, motorHeatFraction %g\n", pump.nBCS, pump.eficM,
                pump.fracTermMotorEfic);
    if (scenario.motorHeatFraction >= 0. && pump.fracTermMotorEfic != scenario.motorHeatFraction)
        setupError("motorHeatFraction not read from the deck");
    MarchResult first{};
    bool equal = true;
    for (int k = 1; k <= kMarches; k++) {
        const MarchResult r = march(cell);
        std::printf("march %d: dpB %.12g kgf/cm2, potB %.12g W, potBT %.12g W, potTermo %.12g W, heat to the cell %.12g W\n",
                    k, r.dpB, r.potB, r.potBT, r.potTermo, r.heat);
        if (k == 1) {
            first = r;
            if (!(r.dpB > 0. && r.potB > 0.))
                setupError("the pump gives no pressure increment at the demo's operating point");
            continue;
        }
        const bool marchEqual = same(r.dpB, first.dpB) && same(r.potB, first.potB) && same(r.potBT, first.potBT) &&
                                same(r.potTermo, first.potTermo) && same(r.heat, first.heat);
        if (!marchEqual)
            std::printf("march %d differs from march 1\n", k);
        equal = equal && marchEqual;
    }
    return equal;
}

} // namespace

int main(int argc, char **argv) {
    std::setvbuf(stdout, nullptr, _IONBF, 0);
    if (argc != 4) {
        std::printf("usage: %s <demos/MultiESP directory> <run root> <scenario>\n", argv[0]);
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

    const std::filesystem::path demo = std::filesystem::absolute(argv[1]);
    const std::filesystem::path runDir = std::filesystem::absolute(std::filesystem::path(argv[2]) / name);
    std::error_code error;
    std::filesystem::create_directories(runDir, error);
    if (!error)
        std::filesystem::current_path(runDir, error); // for files the engine writes to the working directory
    if (error)
        setupError("cannot use the run directory " + runDir.string());
    writeDeck(demo, runDir, scenario->motorHeatFraction);

    std::atexit(reportEngineExit);
    unexpectedExitArmed = true;
    SProd *system = buildSystem(runDir.string() + "/");
    const int k = pumpCell(*system);
    if (k < 0)
        setupError("no multi-stage ESP cell");
    std::printf("multi-stage ESP in cell %d\n", k);
    const bool pass = repeatedMarches(system->celula[k], *scenario);
    unexpectedExitArmed = false;
    std::printf("%s\n", pass ? "PASS" : "FAIL");
    std::fflush(stdout);
    std::_Exit(pass ? kPass : kFail); // skips the engine's static destructors
}
