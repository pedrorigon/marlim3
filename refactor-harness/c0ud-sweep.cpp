/*
 * c0ud-sweep.cpp -- dedicated harness for the five CalcC0Ud* variants.
 *
 * Why this exists
 * ---------------
 * Coverage measured with gcov in stage 0 (evidencia/cobertura.md) found that
 * only CalcC0Ud and CalcC0UdPerm ever execute in the demo corpus. CalcC0UdBuf,
 * CalcC0UdIni and CalcC0UdIniBuf never run: 718 of the 1248 lines and 144 of
 * the 248 conditionals are invisible to L2 and to L3. Stage 3 restructures
 * exactly those five. Six green gates would say nothing about three of them.
 *
 * The five are public members of SProd with signatures the stage preserves
 * (FR-032), so this driver reaches them through the PUBLIC SURFACE rather than
 * through whatever internal shape the refactoring settles on. That is
 * deliberate: the same driver runs before the extraction, after the extraction
 * and after the unification, and a capture from one is comparable with a
 * capture from the others. A harness bound to the new internal API could only
 * be calibrated after the change it is supposed to police.
 *
 * WHAT IS RECORDED, AND WHY NOT ONLY c0/ud
 * ----------------------------------------
 * The epilogue of every variant is
 *
 *     if (arq.escorregaTran == 0) { ...; c0 = 1. + 0*correcaoCo; ud = 0. + 0.*ud*...; }
 *
 * which DISCARDS everything computed above it. A sweep that records only c0 and
 * ud, run with slip disabled, returns 1 and 0 for every input and would accept
 * any corruption whatsoever -- the failure mode of E4-01. So the sweep varies
 * escorregaTran/escorregaPerm over both states, and records every field the
 * five write, which is also what contract guarantee C-B5 is about:
 *
 *     celula[ind].arranjo   celula[ind].transic   celula[ind].transic0
 *     celula[ind].c0Spare   celula[ind].udSpare
 *     celula[ind-1].arranjoR  celula[ind-1].perdaEstratL  celula[ind-1].perdaEstratG
 *
 * Doubles are printed with %a, which is exact for binary64.
 *
 * celula[ind-1] AND ind == 0
 * --------------------------
 * CalcC0Ud and CalcC0UdBuf read celula[ind-1].velPig before any ind > 0 guard,
 * so ind == 0 reads one cell before the array. The driver allocates a guard
 * cell and hands SProd a pointer one past the start, making index -1 defined
 * storage. That keeps ind == 0 in the sweep -- it selects the else arms of the
 * several `if (ind > 0)` -- without the run depending on undefined behaviour.
 *
 * Usage: c0ud-sweep            writes the table to stdout
 *        c0ud-sweep --bench N runs the sweep N times and prints nothing
 *
 * --bench exists because one pass is far under a second, which is too short
 * to time. It drives the SAME code by the same path; it only suppresses the
 * printing, so a timing run and a recording run exercise the same work.
 * Note what it can and cannot say: the sweep spends most of its time inside
 * ProFlu, which is identical on both sides of the comparison, so a measured
 * delta UNDERSTATES the delta of the closure code itself. Directional, not a
 * number to quote.
 */
#include "Leitura.h"
#include "SisProd.h"
#include "estruturas.h"
#include "variaveisGlobais1D.h"

#include <cstdio>
#include <cstdlib>
#include <limits>
#include <fstream>
#include <string>

using namespace std;

// Num4Main.cpp holds main() and cannot be linked in, so the file-scope globals
// it defines are redefined here, purely to satisfy the linker for the objects
// that come along for the ride.
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

constexpr int kCells = 7;   // indices 0..6 usable, plus one guard cell below

// One point of the sweep. Every field here selects a branch catalogued in
// evidencia/c0ud-diff.md; nothing is varied that no variant reads.
struct Scenario {
    const char *name;
    int slip;           // arq.escorregaTran / escorregaPerm
    int pig;            // 0 none, 1 upstream pig, 2 local pig
    int acsrTipo;       // 4 selects the ESP guard, 5 the choke guard
    double freqnova;    // paired with acsrTipo == 4
    double areaGarg;    // paired with acsrTipo == 5
    double qgSign;      // sign multiplier on QG / MCBuf-MliqiniBuf
    double qlSign;      // sign multiplier on QL / Mliqini / MliqiniBuf
    double tiny;        // scales QG down to the localtiny*1e-5 guard
    double alf;         // void fraction, drives hol0 through the 0.01/0.99 window
    double teta;        // pipe inclination, straddles the 45 degree test
    double vtemper;     // sign selects the tmed branch
    int arranjo;        // stored flow pattern: 0, 1, -1 or -2
    int transic;        // transition counter: 0, 1, 19 or 20
    int selDisp, selAnn, selStrat;   // driftSelectors
    int acelera;        // arq.AceleraConvergPerm
    int iterperm;       // SProd::iterperm
    double presRatio;   // scales pressures, moving rgm relative to rlm
    double flowScale;   // scales every mass flow, moving the flow regime
                        // between stratified, slug and annular
    double gasScale;    // scales the GAS flow alone, so ug1/ul1 is not pinned
                        // to one ratio -- mapaTD() and mapaTD(1) disagree only
                        // where ug1 and ul1 are within a factor of about two
    int slipPerm;       // arq.escorregaPerm, separate from escorregaTran: with
                        // one field driving both, swapping them in the epilogue
                        // of CalcC0UdPerm was invisible
    int acsrCell;       // which cell carries the accessory: 9 means all of them,
                        // otherwise only cell ind+acsrCell does. The three
                        // variants read acsr at ind, ind-1 and ind-2, and with
                        // every cell identical those three are indistinguishable
    int nonFinite;      // 0 none, 1 QG = NaN, 2 QG = +inf, 3 QG = -inf.
                        // 0.*QG is NOT zero for these, which is the whole
                        // reason the contract forbids simplifying 0.*QG + 1*QL
    double angEsq;      // left and right junction angles. correcHor's two arms
    double angDir;      // differ ONLY when angDir < 0 <= angEsq: the first arm
                        // needs both negative, the second reads angDir alone.
                        // Leaving both at zero made the arms agree and hid which
                        // cell's accessory the guard reads.
};

// The scenario list. Each row differs from kBase in as few fields as possible,
// so a divergence names the branch it came from. The rows are not a factorial
// product: a full cross of 19 fields is millions of runs and most of it is
// unreachable. What matters is that every branch of every variant is entered
// and left at least once, which the coverage note in the log records.
const Scenario kBase = {"base", 1, 0, 0, 0.5, 1.0, 1.0, 1.0, 1.0, 0.40, 0.30,
                        1.0, 1, 0, 1, 1, 1, 0, 1, 1.0, 1.0, 1.0, 1, 9, 0,
                        0.20, 0.20};

Scenario with(const char *name) { Scenario s = kBase; s.name = name; return s; }

void buildScenarios(Scenario *out, int &count) {
    int n = 0;
    out[n++] = kBase;

    { Scenario s = with("noslip");            s.slip = 0;              out[n++] = s; }
    { Scenario s = with("pig-upstream");      s.pig = 1;               out[n++] = s; }
    { Scenario s = with("pig-local");         s.pig = 2;               out[n++] = s; }
    { Scenario s = with("esp-fast");          s.acsrTipo = 4; s.freqnova = 60.; out[n++] = s; }
    { Scenario s = with("esp-slow");          s.acsrTipo = 4; s.freqnova = 0.5; out[n++] = s; }
    { Scenario s = with("choke-open");        s.acsrTipo = 5; s.areaGarg = 1.0; out[n++] = s; }
    { Scenario s = with("choke-shut");        s.acsrTipo = 5; s.areaGarg = 1e-14; out[n++] = s; }
    { Scenario s = with("qg-negative");       s.qgSign = -1.;          out[n++] = s; }
    { Scenario s = with("ql-negative");       s.qlSign = -1.;          out[n++] = s; }
    { Scenario s = with("both-negative");     s.qgSign = -1.; s.qlSign = -1.; out[n++] = s; }
    { Scenario s = with("qg-tiny");           s.tiny = 1e-9;           out[n++] = s; }
    { Scenario s = with("qg-zero");           s.tiny = 0.;             out[n++] = s; }
    { Scenario s = with("alf-low");           s.alf = 0.005;           out[n++] = s; }
    { Scenario s = with("alf-high");          s.alf = 0.995;           out[n++] = s; }
    { Scenario s = with("alf-mid-low");       s.alf = 0.02;            out[n++] = s; }
    { Scenario s = with("alf-mid-high");      s.alf = 0.98;            out[n++] = s; }
    { Scenario s = with("teta-steep");        s.teta = 1.2;            out[n++] = s; }
    { Scenario s = with("teta-horizontal");   s.teta = 0.;             out[n++] = s; }
    { Scenario s = with("teta-negative");     s.teta = -0.30;          out[n++] = s; }
    { Scenario s = with("teta-near45");       s.teta = 45. * M_PI / 180.; out[n++] = s; }
    { Scenario s = with("vtemper-negative");  s.vtemper = -1.;         out[n++] = s; }
    { Scenario s = with("arranjo-zero");      s.arranjo = 0;           out[n++] = s; }
    { Scenario s = with("arranjo-strat");     s.arranjo = -1;          out[n++] = s; }
    { Scenario s = with("arranjo-annular");   s.arranjo = -2;          out[n++] = s; }
    { Scenario s = with("transic-one");       s.transic = 1;           out[n++] = s; }
    { Scenario s = with("transic-19");        s.transic = 19;          out[n++] = s; }
    { Scenario s = with("transic-20");        s.transic = 20;          out[n++] = s; }
    { Scenario s = with("strat-2");           s.selStrat = 2;          out[n++] = s; }
    { Scenario s = with("ann3-disp1");        s.selAnn = 3; s.selDisp = 1; out[n++] = s; }
    { Scenario s = with("sel-choi");          s.selDisp = 0; s.selAnn = 0; s.selStrat = 0; out[n++] = s; }
    { Scenario s = with("sel-mod");           s.selDisp = 4; s.selAnn = 4; s.selStrat = 4; out[n++] = s; }
    { Scenario s = with("sel-blend");         s.selDisp = 5; s.selAnn = 5; s.selStrat = 5; out[n++] = s; }
    { Scenario s = with("acelera-on");        s.acelera = 1;           out[n++] = s; }
    { Scenario s = with("iterperm-zero");     s.iterperm = 0;          out[n++] = s; }
    { Scenario s = with("gas-rich");          s.presRatio = 0.02;      out[n++] = s; }
    { Scenario s = with("liquid-rich");       s.presRatio = 8.0;       out[n++] = s; }
    { Scenario s = with("noslip-strat");      s.slip = 0; s.arranjo = -1; out[n++] = s; }
    { Scenario s = with("transic-strat");     s.transic = 19; s.arranjo = -1; out[n++] = s; }
    { Scenario s = with("neg-strat");         s.qgSign = -1.; s.arranjo = -1; out[n++] = s; }
    { Scenario s = with("horiz-annular");     s.teta = 0.; s.arranjo = -2; out[n++] = s; }

    // The stratified map returns -1 -- the only way into the block that writes
    // perdaEstratL/perdaEstratG and does the transic bookkeeping -- for
    // volumetric rates ug1 and ul1 between about 1e-4 and 1e-2 m3/s in a
    // HORIZONTAL or DOWNHILL pipe. Measured directly against estratificado
    // over a grid, not guessed: at teta = +0.01 the map returns -1 nowhere in
    // that range, which is why the first set of low-flow rows reached the block
    // exactly once in 1855. ug1 ~ 4.4e-2*flowScale and ul1 ~ 3.9e-3*flowScale,
    // so flowScale between 0.01 and 0.3 puts both inside the window.
    { Scenario s = with("strat-h-30");     s.teta =  0.0;  s.flowScale = 0.30; out[n++] = s; }
    { Scenario s = with("strat-h-10");     s.teta =  0.0;  s.flowScale = 0.10; out[n++] = s; }
    { Scenario s = with("strat-h-03");     s.teta =  0.0;  s.flowScale = 0.03; out[n++] = s; }
    { Scenario s = with("strat-h-01");     s.teta =  0.0;  s.flowScale = 0.01; out[n++] = s; }
    { Scenario s = with("strat-d-10");     s.teta = -0.20; s.flowScale = 0.10; out[n++] = s; }
    { Scenario s = with("strat-d-03");     s.teta = -0.20; s.flowScale = 0.03; out[n++] = s; }
    { Scenario s = with("strat-d-steep");  s.teta = -0.60; s.flowScale = 0.10; out[n++] = s; }

    // The same window crossed with the state the bookkeeping reads. arranjo and
    // transic drive the transic++ / reset arms; selStrat picks mapaTD() against
    // mapaTD(1), which only CalcC0Ud dispatches on.
    { Scenario s = with("strat-tr0");      s.teta = 0.0; s.flowScale = 0.10; s.transic = 0;  s.arranjo = 1;  out[n++] = s; }
    { Scenario s = with("strat-tr1");      s.teta = 0.0; s.flowScale = 0.10; s.transic = 1;  s.arranjo = 1;  out[n++] = s; }
    { Scenario s = with("strat-tr19");     s.teta = 0.0; s.flowScale = 0.10; s.transic = 19; s.arranjo = 1;  out[n++] = s; }
    { Scenario s = with("strat-tr20");     s.teta = 0.0; s.flowScale = 0.10; s.transic = 20; s.arranjo = 1;  out[n++] = s; }
    { Scenario s = with("strat-a0");       s.teta = 0.0; s.flowScale = 0.10; s.arranjo = 0;  out[n++] = s; }
    { Scenario s = with("strat-a0-tr19");  s.teta = 0.0; s.flowScale = 0.10; s.arranjo = 0;  s.transic = 19; out[n++] = s; }
    { Scenario s = with("strat-am1");      s.teta = 0.0; s.flowScale = 0.10; s.arranjo = -1; out[n++] = s; }
    { Scenario s = with("strat-am1-tr19"); s.teta = 0.0; s.flowScale = 0.10; s.arranjo = -1; s.transic = 19; out[n++] = s; }
    { Scenario s = with("strat-am2");      s.teta = 0.0; s.flowScale = 0.10; s.arranjo = -2; out[n++] = s; }
    { Scenario s = with("strat-map2");     s.teta = 0.0; s.flowScale = 0.10; s.selStrat = 2; out[n++] = s; }
    { Scenario s = with("strat-noslip2");  s.teta = 0.0; s.flowScale = 0.10; s.slip = 0;     out[n++] = s; }
    { Scenario s = with("strat-qgneg");    s.teta = 0.0; s.flowScale = 0.10; s.qgSign = -1.; out[n++] = s; }
    { Scenario s = with("strat-qlneg");    s.teta = 0.0; s.flowScale = 0.10; s.qlSign = -1.; out[n++] = s; }
    { Scenario s = with("strat-acel");     s.teta = 0.0; s.flowScale = 0.10; s.acelera = 1;  out[n++] = s; }
    { Scenario s = with("strat-iter0");    s.teta = 0.0; s.flowScale = 0.10; s.iterperm = 0; out[n++] = s; }
    { Scenario s = with("strat-vtemp");    s.teta = 0.0; s.flowScale = 0.10; s.vtemper = -1.; out[n++] = s; }
    { Scenario s = with("strat-blend");    s.teta = 0.0; s.flowScale = 0.10; s.selDisp = 5; s.selAnn = 5; s.selStrat = 5; out[n++] = s; }

    { Scenario s = with("high-flow");         s.flowScale = 50.; out[n++] = s; }
    { Scenario s = with("high-flow-horiz");   s.flowScale = 50.; s.teta = 0.01; out[n++] = s; }

    // --- rows added after calibration found the sweep blind to four defects ---

    // The epilogue of CalcC0UdPerm reads arq.escorregaPerm while the other four
    // read arq.escorregaTran. With one field driving both, swapping them was
    // invisible: the two states have to be independent.
    { Scenario s = with("slip-tran-only");  s.slip = 1; s.slipPerm = 0; out[n++] = s; }
    { Scenario s = with("slip-perm-only");  s.slip = 0; s.slipPerm = 1; out[n++] = s; }

    // correcHor reads celula[ind-1].acsr in CalcC0Ud and celula[ind].acsr in the
    // other three; the ang adjustment reads celula[ind-2].acsr. Giving every
    // cell the same accessory made all three read the same value, so swapping
    // the index changed nothing. Each row below puts the accessory on exactly
    // one of the three cells, in a horizontal pipe where correcHor is reached.
    for (int off = -2; off <= 0; off++) {
        { Scenario s = with(off == -2 ? "acsr-at-ind-2-chk"
                          : off == -1 ? "acsr-at-ind-1-chk" : "acsr-at-ind-chk");
          s.teta = 0.; s.acsrCell = off; s.acsrTipo = 5; s.areaGarg = 1e-14; out[n++] = s; }
        { Scenario s = with(off == -2 ? "acsr-at-ind-2-open"
                          : off == -1 ? "acsr-at-ind-1-open" : "acsr-at-ind-open");
          s.teta = 0.; s.acsrCell = off; s.acsrTipo = 5; s.areaGarg = 1.0; out[n++] = s; }
        { Scenario s = with(off == -2 ? "acsr-at-ind-2-esp"
                          : off == -1 ? "acsr-at-ind-1-esp" : "acsr-at-ind-esp");
          s.teta = 0.; s.acsrCell = off; s.acsrTipo = 4; s.freqnova = 60.; out[n++] = s; }
    }

    // 0.*QG + 1*QL equals QL exactly for every FINITE QG, so a sweep of finite
    // values cannot tell the written form from the simplified one. It is
    // precisely for NaN and infinity that they differ, and that is the reason
    // the contract forbids the simplification -- so the guard has to be
    // exercised where it bites.
    { Scenario s = with("qg-nan");        s.nonFinite = 1; out[n++] = s; }
    { Scenario s = with("qg-posinf");     s.nonFinite = 2; out[n++] = s; }
    { Scenario s = with("qg-neginf");     s.nonFinite = 3; out[n++] = s; }
    { Scenario s = with("qg-nan-horiz");  s.nonFinite = 1; s.teta = 0.; out[n++] = s; }
    { Scenario s = with("qg-nan-strat");  s.nonFinite = 1; s.teta = 0.; s.flowScale = 0.10; out[n++] = s; }
    { Scenario s = with("qg-inf-strat");  s.nonFinite = 2; s.teta = 0.; s.flowScale = 0.10; out[n++] = s; }

    // CalcC0Ud dispatches mapaTD() against mapaTD(1) on driftSelectors.stratified;
    // CalcC0UdIni and CalcC0UdPerm always call mapaTD(). Measured over 1792
    // combinations, the two disagree in 8 -- all of them where ug1 and ul1 are
    // within a factor of about two, which the fixed 11:1 ratio of the rows above
    // never produces. gasScale decouples them.
    for (int gi = 0; gi < 3; gi++) {
        const double gs[3] = {0.05, 0.09, 0.2};
        for (int fi = 0; fi < 3; fi++) {
            const double fs[3] = {0.10, 0.25, 0.50};
            for (int ai = 0; ai < 3; ai++) {
                const double av[3] = {0.20, 0.30, 0.50};
                static char names[81][24];
                static int slot = 0;
                snprintf(names[slot], sizeof(names[0]), "mapbound-%d%d%d", gi, fi, ai);
                Scenario s = with(names[slot]); slot++;
                s.teta = (ai == 0) ? -0.05 : (ai == 1 ? 0.0 : -0.20);
                s.flowScale = fs[fi]; s.gasScale = gs[gi]; s.alf = av[ai];
                out[n++] = s;
            }
        }
    }

    // correcHor: the two arms of the accessory guard produce a DIFFERENT
    // correcHor only when angDir < 0 <= angEsq -- the first arm requires both
    // angles negative and does nothing here, the second reads angDir alone and
    // yields -1. Crossed with the accessory sitting on exactly one of the three
    // cells the variants read, this is what tells celula[ind-1].acsr from
    // celula[ind].acsr.
    for (int off = -2; off <= 0; off++) {
        static char cn[6][28]; static int cs = 0;
        snprintf(cn[cs], sizeof(cn[0]), "corrhor-chk-%d", -off);
        { Scenario s = with(cn[cs++]); s.teta = 0.; s.acsrCell = off;
          s.acsrTipo = 5; s.areaGarg = 1e-14;
          s.angEsq = 0.20; s.angDir = -0.20; out[n++] = s; }
        snprintf(cn[cs], sizeof(cn[0]), "corrhor-open-%d", -off);
        { Scenario s = with(cn[cs++]); s.teta = 0.; s.acsrCell = off;
          s.acsrTipo = 5; s.areaGarg = 1.0;
          s.angEsq = 0.20; s.angDir = -0.20; out[n++] = s; }
    }
    { Scenario s = with("corrhor-both-neg"); s.teta = 0.; s.angEsq = -0.2; s.angDir = -0.2; out[n++] = s; }
    { Scenario s = with("corrhor-both-pos"); s.teta = 0.; s.angEsq =  0.2; s.angDir =  0.2; out[n++] = s; }

    // 0.*QG + 1*QL differs from QL only when QG is non-finite, AND the guard it
    // feeds only changes betI when the two disagree in SIGN. With QG = NaN the
    // comparison is false either way unless QL is negative, so the non-finite
    // rows have to be crossed with a negative liquid rate.
    for (int nf = 1; nf <= 3; nf++) {
        static char qn[6][24]; static int qs = 0;
        snprintf(qn[qs], sizeof(qn[0]), "qg-nf%d-qlneg", nf);
        { Scenario s = with(qn[qs++]); s.nonFinite = nf; s.qlSign = -1.; out[n++] = s; }
        snprintf(qn[qs], sizeof(qn[0]), "qg-nf%d-qlneg-h", nf);
        { Scenario s = with(qn[qs++]); s.nonFinite = nf; s.qlSign = -1.;
          s.teta = 0.; s.flowScale = 0.10; out[n++] = s; }
    }

    // A denser crossing of the map boundary. mapaTD() and mapaTD(1) disagree in
    // 8 of 1792 measured combinations, so the grid has to be fine enough to land
    // in one: gasScale sweeps ug1 across two decades while flowScale moves ul1.
    for (int gi = 0; gi < 7; gi++) {
        const double gs[7] = {0.02, 0.04, 0.06, 0.09, 0.13, 0.20, 0.35};
        for (int fi = 0; fi < 4; fi++) {
            const double fs[4] = {0.08, 0.15, 0.25, 0.45};
            for (int ai = 0; ai < 3; ai++) {
                const double av[3] = {0.20, 0.30, 0.50};
                const double tv[3] = {0.0, -0.05, -0.20};
                static char names[128][24]; static int slot = 0;
                snprintf(names[slot], sizeof(names[0]), "mapfine-%d%d%d", gi, fi, ai);
                Scenario s = with(names[slot]); slot++;
                s.teta = tv[ai]; s.alf = av[ai];
                s.flowScale = fs[fi]; s.gasScale = gs[gi];
                out[n++] = s;
            }
        }
    }

    count = n;
}

// Rebuilds the whole array from the scenario. Called before EVERY invocation,
// so each recorded row is a function of (scenario, ind, variant) alone and the
// table does not depend on the order the variants are called in.
void resetCells(SProd &sp, Cel *cells, const Scenario &s, int ind) {
    for (int i = -1; i < kCells; i++) {
        Cel &c = cells[i];
        const double k = 1.0 + 0.05 * (i + 1);

        // The engine's own fallback fluid (Leitura.cpp:12624), used verbatim so
        // the harness fluid is one the program itself considers valid. A
        // default-constructed ProFlu has api=1, rgo=1 and templ==temph, which
        // makes ViscOleo return NaN and MasEspLiq negative.
        c.flui = ProFlu(sp.vg1dSP, 20, 100., 0.7, 0., 1., 20, 10, 40,
                        2, 0, 0.8, 0, 0, 0);
        c.fluiL = &c.flui;
        c.duto.a = 0.1 * k;
        c.duto.area = M_PI * (0.1 * k) * (0.1 * k) / 4.;
        c.duto.rug = 4.5e-5;
        c.duto.teta = s.teta;
        c.dutoL.teta = s.teta * 0.75;
        c.dutoR.teta = s.teta * 1.25;
        c.dx = 10. * k;
        c.dxL = 12. * k;

        c.alf = s.alf;
        c.alfPigE = s.alf * 0.9;
        c.alfPigD = s.alf * 1.05 > 1. ? 1. : s.alf * 1.05;

        // MC is total mass flow, Mliqini its liquid part, so MC - Mliqini is the
        // gas flow the variants take the sign of. Keeping MC > Mliqini makes
        // that difference positive before qgSign is applied.
        c.QG = s.qgSign * 2.0 * k * s.tiny * s.flowScale * s.gasScale;
        c.QL = s.qlSign * 3.0 * k * s.flowScale;
        c.Mliqini = s.qlSign * 3.0 * k * s.flowScale;
        c.MC = c.Mliqini + s.qgSign * 2.0 * k * s.tiny * s.flowScale * s.gasScale;
        c.MliqiniBuf = s.qlSign * 3.5 * k * s.flowScale;
        c.MCBuf = c.MliqiniBuf + s.qgSign * 2.2 * k * s.tiny * s.flowScale * s.gasScale;
        c.MliqiniLBuf = s.qlSign * 3.25 * k * s.flowScale;

        // Only QG carries the non-finite value: it is the operand of the 0.*QG
        // term the contract protects. MC and Mliqini stay finite so the rest of
        // the body still reaches its branches instead of drowning in NaN.
        if (s.nonFinite == 1) c.QG = std::numeric_limits<double>::quiet_NaN();
        if (s.nonFinite == 2) c.QG = std::numeric_limits<double>::infinity();
        if (s.nonFinite == 3) c.QG = -std::numeric_limits<double>::infinity();

        c.dt = 1.0;
        // celula3.h:117 -- pressure in kgf/cm2, temperature in degrees Celsius,
        // mass flow in kg/s. Feeding SI values here puts every fluid property
        // far outside its correlation range and NaNs the Reynolds number, which
        // is how the first version of this sweep reached no branch at all.
        c.pres = 60. * s.presRatio * k;
        c.presaux = 62. * s.presRatio * k;
        c.presauxL = 61. * s.presRatio * k;
        c.presBuf = 64. * s.presRatio * k;
        c.presLBuf = 63. * s.presRatio * k;
        c.temp = 55. * k;
        c.tempL = 52. * k;
        c.VTemper = s.vtemper;

        c.betL = 0.10;
        c.bet = 0.12;
        c.betPigE = 0.14;
        c.betPigD = 0.16;

        c.rpCi = 920. * k; c.rcCi = 1000. * k; c.rgCi = 55. * k;
        c.rpLi = 915. * k; c.rcLi = 995. * k; c.rgLi = 54. * k;

        c.angEsq = s.angEsq;
        c.angDir = s.angDir;
        c.estabCol = 0;
        // arranjo and transic are READ by the bodies, so they carry the scenario
        // value; arranjoR and transic0 are write-only and carry sentinels.
        c.arranjo = s.arranjo;
        c.arranjoR = 77;
        c.transic = s.transic;
        c.transic0 = -95;
        c.c0 = 1.15;
        c.ud = 0.25;
        // Sentinels, not zeros. Every field recorded below is WRITE-ONLY in the
        // five bodies (verified: 10 arranjoR, 6 perdaEstratL, 6 perdaEstratG,
        // 6 c0Spare, 6 udSpare, 1 transic0 -- all on the left of an =), so a
        // value that survives the call means the write did not happen. Seeding
        // them with 0 made "wrote zero" and "never written" identical, and the
        // stratified block writes fatorperdaLiq == 0 on most inputs: the sweep
        // reported perdaEstratL untouched while the block was in fact running.
        c.c0Spare = -97.;
        c.udSpare = -96.;
        c.perdaEstratL = -99.;
        c.perdaEstratG = -98.;

        c.velPig = 0.;
        c.estadoPig = 0;
        if (s.pig == 1) { c.velPig = 1.; c.estadoPig = 1; }
        if (s.pig == 2) { c.velPig = -1.; c.estadoPig = 1; }

        // acsrCell == 9 gives every cell the same accessory; otherwise only
        // cell ind+acsrCell carries it and the rest are neutral, which is what
        // makes celula[ind].acsr, celula[ind-1].acsr and celula[ind-2].acsr
        // tell each other apart.
        if (s.acsrCell == 9 || i == ind + s.acsrCell) {
            c.acsr.tipo = s.acsrTipo;
            c.acsr.bcs.freqnova = s.freqnova;
            c.acsr.chk.AreaGarg = s.areaGarg;
        } else {
            c.acsr.tipo = 0;
            c.acsr.bcs.freqnova = 0.5;
            c.acsr.chk.AreaGarg = 1.0;
        }
    }

    sp.ncel = kCells - 1;
    sp.iterperm = s.iterperm;
    sp.tGSup = 25.;
    sp.alfE = s.alf * 0.95;
    sp.betaE = 0.13;
    sp.presE = 61. * s.presRatio;
    sp.tempE = 53.;
    sp.driftSelectors.dispersed = s.selDisp;
    sp.driftSelectors.annularChurn = s.selAnn;
    sp.driftSelectors.stratified = s.selStrat;
    sp.arq.escorregaTran = s.slip;
    sp.arq.escorregaPerm = s.slipPerm;
    sp.arq.AceleraConvergPerm = s.acelera;
    sp.arq.mapaArranjo = 0;
}

bool gQuiet = false;

void record(const char *variant, const Scenario &s, int ind,
            double c0, double ud, const Cel *cells) {
    if (gQuiet)
        return;
    const Cel &here = cells[ind];
    const Cel &prev = cells[ind - 1];
    printf("%-14s %-18s ind=%d c0=%a ud=%a arranjo=%d transic=%d transic0=%d "
           "c0Spare=%a udSpare=%a arranjoR=%d perdaL=%a perdaG=%a\n",
           variant, s.name, ind, c0, ud,
           static_cast<int>(here.arranjo), here.transic, here.transic0,
           here.c0Spare, here.udSpare,
           static_cast<int>(prev.arranjoR), prev.perdaEstratL, prev.perdaEstratG);
}

}  // namespace

int main(int argc, char **argv) {
    long repeats = 1;
    bool quiet = false;
    if (argc == 3 && std::string(argv[1]) == "--bench") {
        repeats = std::atol(argv[2]);
        quiet = true;
    }
    varGlob1D globals;

    // Deliberately leaked, as in trend-sweep.cpp: ~SProd frees arrays keyed off
    // counters this driver sets by hand, and the process is about to exit.
    SProd &sp = *new SProd();
    sp.vg1dSP = &globals;

    // One guard cell below index 0, so celula[-1] is defined storage.
    Cel *storage = new Cel[kCells + 1];
    Cel *cells = storage + 1;
    sp.celula = cells;

    gQuiet = quiet;
    static Scenario scenarios[512];
    int scenarioCount = 0;
    buildScenarios(scenarios, scenarioCount);

    for (long pass = 0; pass < repeats; pass++)
    for (int si = 0; si < scenarioCount; si++) {
        const Scenario &s = scenarios[si];
        for (int ind = 0; ind <= kCells - 1; ind++) {
            double c0, ud;

            resetCells(sp, cells, s, ind); c0 = -7.; ud = -7.;
            sp.CalcC0Ud(ind, c0, ud);        record("CalcC0Ud", s, ind, c0, ud, cells);

            resetCells(sp, cells, s, ind); c0 = -7.; ud = -7.;
            sp.CalcC0UdBuf(ind, c0, ud);     record("CalcC0UdBuf", s, ind, c0, ud, cells);

            resetCells(sp, cells, s, ind); c0 = -7.; ud = -7.;
            sp.CalcC0UdIni(ind, c0, ud);     record("CalcC0UdIni", s, ind, c0, ud, cells);

            resetCells(sp, cells, s, ind); c0 = -7.; ud = -7.;
            sp.CalcC0UdIniBuf(ind, c0, ud);  record("CalcC0UdIniBuf", s, ind, c0, ud, cells);

            resetCells(sp, cells, s, ind); c0 = -7.; ud = -7.;
            sp.CalcC0UdPerm(ind, c0, ud);    record("CalcC0UdPerm", s, ind, c0, ud, cells);
        }
    }
    return 0;
}
