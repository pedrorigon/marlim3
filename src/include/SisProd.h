/*
 * SisProd.h
 *
 * Created on: December 21, 2016
 *     Author: Eduardo
 */

#ifndef SISPROD_H_
#define SISPROD_H_
#define _USE_MATH_DEFINES // Enables M_PI on supported platforms

#include "Acidentes2.h"
#include "DriftFluxClosure.h"
#include "Bcsm2.h"
#include "BombaVol.h"
#include "FerramentasNumericas.h"
#include "FonteMas.h"
#include "FonteMassCHK.h"
#include "Geometria.h"
#include "GradientCorrelations.h"
#include "Leitura.h"
#include "Log.h"
#include "Matriz.h"
#include "PropFlu.h"
#include "PropFluCol.h"
#include "TrocaCalor.h"
#include "Vetor.h"
#include "acessorios.h"
#include "celula3.h"
#include "celulaGas.h"
#include "chokegas.h"
#include "criterioIntermiSevera.h"
#include "dados3DPoisson.h"
#include "estrat.h"
#include "estruturaTabDin.h"
#include "mapa.h"
#include "multiBCS.h"
#include "solver3DPoisson.h"
#include "variaveisGlobais1D.h"
#include <ctime>
#include <fstream>
#include <iostream>
#include <math.h>
#include <omp.h>
#include <sstream>
#include <stdio.h>
#include <string.h>
#include <string>
#include <vector>
#include <algorithm>
#include <utility>

using namespace std;

/// Application or simulator version string.
extern string versao;
/// Global simulation start timestamp.
extern time_t nowGlobIni;
/// Broken-down local time corresponding to nowGlobIni.
extern tm *ltmGlobIni;
/// Day component of the simulation start time.
extern int diaIni;
/// Hour component of the simulation start time.
extern int horaIni;
/// Minute component of the simulation start time.
extern int minutoIni;
/// Second component of the simulation start time.
extern int segundoIni;
/// Global simulation end timestamp.
extern time_t nowGlobFim;
/// Broken-down local time corresponding to nowGlobFim.
extern tm *ltmGlobFim;

// The state adapters and the modules' updaters read members of SProd that no
// consumer reads. They are friends, so those members can stay private. Their
// types are only declared here: the files that define them include the module
// headers, and consumers never need them.
class SProd;
namespace sisprod::composition { struct CompositionState; struct CompositionUpdaters; }
namespace sisprod::gaslift { struct GasLiftState; struct GasLiftTemperatureUpdater; }
namespace sisprod::sources { struct SourceState; }
namespace sisprod::steady { struct SteadyStateSearchState; struct SteadyStateState; struct SteadyStateUpdaters; }
namespace sisprod::thermal { struct ThermalClosureUpdater; struct ThermalSourceUpdater; struct ThermalState; }
namespace sisprod::transient { struct TransientSolveState; struct TransientSolveUpdaters; struct TransientStepState; struct TransientStepUpdaters; }
namespace trendoutput { struct TrendState; }
namespace sisprod::adapters {
sisprod::gaslift::GasLiftState gasLiftStateOf(SProd &system);
sisprod::steady::SteadyStateState steadyStateOf(SProd &system);
sisprod::steady::SteadyStateSearchState searchStateOf(SProd &system);
sisprod::transient::TransientStepState transientStateOf(SProd &system);
sisprod::composition::CompositionState compositionStateOf(SProd &system);
sisprod::transient::TransientSolveState transientSolveStateOf(SProd &system);
sisprod::thermal::ThermalState thermalStateOf(SProd &system);
sisprod::sources::SourceState sourceStateOf(SProd &system);
driftflux::coefficient::ClosureState closureStateOf(SProd &system);
trendoutput::TrendState trendStateOf(const SProd &system);
}  // namespace sisprod::adapters

namespace sisprod {

/// A square table the fluids read through a double **: its rows, and a pointer to each.
struct SquareTable {
    vector<vector<double>> rows;
    vector<double *> rowPointers;

    /// n rows of n values each.
    void allocate(int n) {
        rows.assign(n, vector<double>(n));
        rowPointers.resize(n);
        for (int i = 0; i < n; i++)
            rowPointers[i] = rows[i].data();
    }
    /// Destroys the rows and gives their memory back.
    void release() {
        vector<vector<double>>().swap(rows);
        vector<double *>().swap(rowPointers);
    }
    /// What the fluids are given: the row pointers, null when there is no table.
    double **data() { return rowPointers.data(); }
    double *operator[](int i) { return rowPointers[i]; }
};

/// The fluid-property tables the fluids of a system's cells point at. The black-oil
/// compressibility-factor, specific-heat and liquid-density-derivative tables belong to the
/// input; the latent-heat table and the Livia bubble-point and solution-gas-ratio tables are
/// built and owned by the system, and so are the dynamic tables of a compositional network.
struct PropertyTables {
    /**
     * @brief Black-oil gas-compressibility-factor table.
     */
    double **zdranP = nullptr;
    /**
     * @brief Pressure derivative of the black-oil compressibility-factor table.
     */
    double **dzdpP = nullptr;
    /**
     * @brief Temperature derivative of the black-oil compressibility-factor table.
     */
    double **dzdtP = nullptr;
    /**
     * @brief Black-oil gas specific-heat table.
     */
    double **cpg = nullptr;
    /**
     * @brief Black-oil produced-liquid specific-heat table.
     */
    double **cpl = nullptr;
    /**
     * @brief Temperature derivative of liquid density.
     */
    double **drholdT = nullptr;
    /**
     * @brief Black-oil latent-heat table.
     */
    SquareTable HLat;
    /**
     * @brief Bubble-pressure values imported from PVTSim for the Livia solution-gas-ratio correlation. This table
     * is separate from the full PVTSim fluid-property model.
     */
    vector<double> PBPVTSim;
    /**
     * @brief Bubble-temperature values imported from PVTSim for the Livia solution-gas-ratio correlation. This
     * table is separate from the full PVTSim fluid-property model.
     */
    vector<double> TBPVTSim;
    /**
     * @brief Precomputed solution-gas-ratio table for the Livia correlation, used to avoid repeating its
     * expensive calculation during the simulation.
     */
    SquareTable RSLivia;
    /**
     * @brief Dynamic fluid-property tables.
     */
    vector<tabelaDinamica> tabDin;
    /**
     * @brief Number of dynamic property tables.
     */
    int ntabDin = 0;
};

/// The gas-lift line: its cells, the injection choke, the valves and the cells they join, its
/// pressure-velocity system, its surface and initial conditions, the column-annulus thermal
/// coupling, the completion-fluid interface, and the unloading controller's samples and
/// settings.
struct GasLiftLine {
    /**
     * @brief Global pressure-velocity coupling matrix for the gas line.
     */
    BandMtx<double> matglobG;
    /**
     * @brief Right-hand side and solution vector for the gas-line pressure-velocity system.
     */
    Vcr<double> termolivreG;
    /**
     * @brief Production-line cell indices associated with gas-lift valves.
     */
    vector<int> posicVGLP;
    /**
     * @brief Service-line cell indices associated with gas-lift valves.
     */
    vector<int> posicVGLG;
    /**
     * @brief Service-line index where column-annulus thermal coupling begins.
     */
    int AnulaColunaIni = 0;
    /**
     * @brief Service-line index where column-annulus thermal coupling ends.
     */
    int AnulaColunaFim = 0;
    /**
     * @brief Production-column index aligned with AnulaColunaIni.
     */
    int ColunaAnulaIni = 0;
    /**
     * @brief Production-column index aligned with AnulaColunaFim.
     */
    int ColunaAnulaFim = 0;
    /**
     * @brief Indicates whether column-annulus thermal coupling is enabled.
     */
    int verificaAcop = 0;
    /**
     * @brief Gas-lift valves installed in the system.
     */
    vector<ChokeGas> chokeVGL;
    /**
     * @brief Time horizon used by the gas-lift unloading controller.
     */
    double tempMedContDesc = 10.;
    /**
     * @brief Maximum number of flow samples retained by the unloading PI controller.
     */
    double maxVecContDesc = 1000;
    /**
     * @brief Average maximum completion-fluid mass flow through the gas-lift valves.
     */
    double vazmedDesc = 0;
    /**
     * @brief Averaging interval used for gas-lift-valve mass flow.
     */
    double tempmedDEsc = 0;
    /**
     * @brief Maximum valve mass-flow samples used by the unloading controller.
     */
    vector<double> vazmaxMedDesc;
    /**
     * @brief Time-step samples associated with the unloading-flow average.
     */
    vector<double> dtDesc;
    /**
     * @brief Separator pressure or pressure at the inlet of the downstream network section.
     */
    double pGSup = 0;
    /**
     * @brief Number of control volumes in the gas service line.
     */
    int ncelGas = 0;
    /**
     * @brief Gas-injection pressure. Zero until an injection-pressure condition or a sensitivity
     * analysis case sets it; the steady summary and the restart file record it either way.
     */
    double presiniG = 0.;
    /**
     * @brief Gas-injection temperature in the service line. Zero until an injection condition
     * or a sensitivity analysis case sets it; the restart file records it either way.
     */
    double tempiniG = 0.;
    /**
     * @brief Gas-injection choke model.
     */
    ChokeGas chokeInj;
    /**
     * @brief Gas service-line control volumes, held by gasCells.
     */
    CelG *celulaG = nullptr;
    /**
     * @brief The gas service-line control volumes celulaG points at.
     */
    vector<CelG> gasCells;
    /**
     * @brief Current service-line cell containing the completion-fluid/gas interface.
     */
    int celInter = 1e7;
    /**
     * @brief Maximum time step that keeps the unloading interface within one gas-line cell.
     */
    double dtInter = 0.;
    /**
     * @brief Current completion-fluid/gas interface velocity.
     */
    double velInter = 0.;
    /**
     * @brief Previous-time-level interface cell.
     */
    int celInterIni = 0.;
    /**
     * @brief Previous-time-level interface-limited time step.
     */
    double dtInterIni = 0.;
    /**
     * @brief Previous-time-level interface velocity.
     */
    double velInterIni = 0.;
};

/// The state of a transient run: the time step, its history and the restrictions on it, the
/// moving averages and running totals, the Master1 valve state and schedule, the
/// production-line pressure-velocity system, the buffered sources, the event and log
/// counters, the transport switches, the pigs and the cells of the two-dimensional Poisson
/// model.
struct TransientRun {
    /**
     * @brief Enables transport equations for primitive black-oil properties, including API gravity, BSW, gas-oil
     * ratio, and light/heavy mass fractions.
     */
    int trackRGO = 0;
    /**
     * @brief Enables transport equations for gas density and the gas-phase CO2 molar fraction.
     */
    int trackDeng = 0;
    /**
     * @brief Current Master1 valve state.
     */
    int EstadoMaster1 = 0;
    /**
     * @brief Counter used while changing the Master1 state.
     */
    int contaMaster1 = 0;
    /**
     * @brief Number of pigs scheduled for launch.
     */
    int npig = 0;
    /**
     * @brief Cell indices where pigs are received.
     */
    vector<int> receb;
    /**
     * @brief Global pressure-velocity coupling matrix for the multiphase production line.
     */
    BandMtx<double> matglobP;
    /**
     * @brief Number of scheduled Master1 opening events.
     */
    int nabreM1 = 0;
    /**
     * @brief Number of scheduled Master1 closing events.
     */
    int nfechaM1 = 0;
    /**
     * @brief Times at which Master1 closes.
     */
    vector<double> fechaM1;
    /**
     * @brief Times at which Master1 opens.
     */
    vector<double> abreM1;
    /**
     * @brief Number of events written to the event log.
     */
    int contaLog = 0;
    /**
     * @brief Smallest production-line control-volume length.
     */
    double menorDx = 0.;
    /**
     * @brief Simulation time-step counter.
     */
    int kSP = 0.;
    /**
     * @brief Controls event-log output frequency.
     */
    int KontaImprime = 0.;
    /**
     * @brief Index of the next scheduled simulation event.
     */
    int indevento = 0.;
    /**
     * @brief History of recently accepted time steps.
     */
    vector<double> dtSim;
    /**
     * @brief History of time steps proposed by the CFL criterion.
     */
    vector<double> dtCFL;
    /**
     * @brief Average time step proposed by the CFL criterion.
     */
    double dtCFLMed = 1.;
    /**
     * @brief Average time step actually used by the simulation.
     */
    double dtSimMed = 1.;
    /**
     * @brief Indicates that time-step growth must remain restricted.
     */
    int restriDt = 0;
    /**
     * @brief Number of remaining steps under the current time-step restriction.
     */
    int kontarestriDt = 0;
    /**
     * @brief Accumulated CFL time steps used to compute dtCFLMed.
     */
    double dtCFLTotal = 0.;
    /**
     * @brief Accumulated accepted time steps used to compute dtSimMed.
     */
    double dtSimTotal = 0.;
    /**
     * @brief Counts alternating liquid-flow oscillations near an active surface choke.
     */
    int kontaGolfada = 1000;
    /**
     * @brief Multiplier applied to a small artificial upstream gas flow during difficult Master1 closures. It
     * mitigates a pressure blind spot when the upstream side contains only liquid and increases if the
     * problem persists.
     */
    double momentoDesesp = 0;
    /**
     * @brief Recent maximum pressure-change rates.
     */
    vector<double> taxaDpMax;
    /**
     * @brief Average maximum pressure-change rate.
     */
    double DpMaxMed = 1.;
    /**
     * @brief Recent maximum temperature-change rates.
     */
    vector<double> taxaDTMax;
    /**
     * @brief Average maximum temperature-change rate.
     */
    double DTMaxMed = 1.;
    /**
     * @brief Counter controlling compositional-property refreshes.
     */
    int kontaRenovaComp = 0;
    /**
     * @brief Production cells handled by the two-dimensional Poisson model.
     */
    vector<int> indCelPoisson2D;
    /**
     * @brief Number of cells handled by the two-dimensional Poisson model.
     */
    int nCelulaPoisson2D = 0;
    /**
     * @brief Signals that the time step must be adjusted.
     */
    int alteraTempo = 0;
    /**
     * @brief Requests rollback and time-step reevaluation when at least one control volume produces a holdup or
     * volume fraction outside the physical [0, 1] range.
     */
    int reinicia = 0;
    /**
     * @brief Recent Master1 ratio values for the active state.
     */
    double vRazMast1[10] = {};
    /**
     * @brief Recent Master1 ratio values for the inactive state.
     */
    double vRazMast0[10] = {};
    /**
     * @brief Critical Master1 ratio history used by the switching logic.
     */
    double vRazMastCrit[10] = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
    /**
     * @brief Previous-time-level value of pGSup.
     */
    double pGSupIni = 0.;
    /**
     * @brief Kept for the restart file, which records it; nothing in the
     * simulation reads it.
     */
    double massfonte = 0.;
    /**
     * @brief CFL safety factor, typically set to 0.8.
     */
    double mult = 0;
    /**
     * @brief Time-averaged pressure in the final production-line control volume. Used to decide whether the
     * surface choke behaves as a localized pressure loss or as a discharge-flow model.
     */
    double presMedMov = 0;
    /**
     * @brief Time-averaged mixture volumetric flux in the final production-line control volume. Used by the
     * surface-choke operating-mode logic.
     */
    double jMedMov = 0;
    /**
     * @brief Time-averaged void fraction in the final production-line control volume, updated once the
     * averaging window fills. Zero until then; only the restart file reads it.
     */
    double alfMedMov = 0.;
    /**
     * @brief Start time of the moving-average window.
     */
    double tMedMov = 0;
    /**
     * @brief Duration of the moving-average window.
     */
    double ktMedMov = 0.;
    /**
     * @brief Accumulated pressure used to compute presMedMov.
     */
    double pTotal = 0.;
    /**
     * @brief Accumulated mixture flux used to compute jMedMov.
     */
    double jTotal = 0.;
    /**
     * @brief Accumulated void fraction used to compute alfMedMov.
     */
    double alfTotal = 0.;
    /**
     * @brief Pressure samples used by the moving-average calculation.
     */
    vector<double> presVet;
    /**
     * @brief Mixture-flux samples used by the moving-average calculation.
     */
    vector<double> jVet;
    /**
     * @brief Void-fraction samples used by the moving-average calculation.
     */
    vector<double> alfVet;
    /**
     * @brief Time samples associated with the moving-average window.
     */
    vector<double> tVet;
    /**
     * @brief Current surface-choke open/closed state.
     */
    int aberto = 0;
    /**
     * @brief Previous-time-level surface-choke state.
     */
    int abertoini = 0;
    /**
     * @brief Counter that delays transitions out of active-choke mode.
     */
    int tempoaberto = 0;
    /**
     * @brief Previous-time-level value of tempoaberto.
     */
    int tempoabertoini = 0;
    /**
     * @brief Signals a surface-choke operating-mode transition.
     */
    int mudaModoChk = 0;
    /**
     * @brief Selects the interphase mass-transfer model: 0 = complete, 1 = fully explicit, 2 = simplified, and 3
     * = disabled.
     */
    int TransMassModel = 0;
    /**
     * @brief Number of pigs currently moving through the line.
     */
    int indpigP = 0;
    /**
     * @brief Previous-time-level value of indpigP.
     */
    int indpigPini = 0;
    /**
     * @brief Right-hand side and solution vector for the production-line pressure-velocity system.
     */
    Vcr<double> termolivreP;
    /**
     * @brief Simulation end time.
     */
    double tfinal = 0;
    /**
     * @brief Current complete-model activation state.
     */
    int modeloCompleto = 1;
    /**
     * @brief Auxiliary CFL time-step accumulator.
     */
    double dtauxCFL = 0.;
    /**
     * @brief Auxiliary accepted-time-step accumulator.
     */
    double dtauxFinal = 0.;
    /**
     * @brief Intermediate production-liquid inflow estimate for a network section.
     */
    double fontemassPRBuf = 0.;
    /**
     * @brief Intermediate complementary-liquid inflow estimate for a network section.
     */
    double fontemassCRBuf = 0.;
    /**
     * @brief Intermediate gas inflow estimate for a network section.
     */
    double fontemassGRBuf = 0.;
    /**
     * @brief Minimum time step allowed by the current cycle. Zero until a
     * transient cycle sets it.
     */
    double dtCicMin = 0.;
};

/// The storage of one family of trends: for each trend, its rows of samples, and its reset time
/// and sample counts. The solve reads it through the pointers allocate() points at it.
struct TrendSet {
    vector<vector<vector<double>>> values;
    vector<vector<double *>> rows;
    vector<double **> matrices;
    vector<double> resetTimers;
    vector<int> counts;
    vector<int> bufferedCounts;

    /// For each of the count trends, length(i) rows of rowWidth(i) values whose first
    /// sentinelCount(i) hold the -10000 sentinel, and the reset timers and sample counts,
    /// zeroed; points matrixPointer, resetTimerPointer, countPointer and bufferedCountPointer
    /// at them.
    template <typename Length, typename RowWidth, typename SentinelCount>
    void allocate(int count, Length length, RowWidth rowWidth, SentinelCount sentinelCount, double ***&matrixPointer,
                  double *&resetTimerPointer, int *&countPointer, int *&bufferedCountPointer) {
        values.resize(count);
        rows.resize(count);
        matrices.resize(count);
        for (int i = 0; i < count; i++) {
            const int n = length(i);
            values[i].assign(n, vector<double>(rowWidth(i)));
            rows[i].resize(n);
            for (int j = 0; j < n; j++) {
                rows[i][j] = values[i][j].data();
                for (int k = 0; k < sentinelCount(i); k++)
                    values[i][j][k] = -10000.;
            }
            matrices[i] = rows[i].data();
        }
        resetTimers.assign(count, 0.);
        counts.assign(count, 0);
        bufferedCounts.assign(count, 0);
        matrixPointer = matrices.data();
        resetTimerPointer = resetTimers.data();
        countPointer = counts.data();
        bufferedCountPointer = bufferedCounts.data();
    }
    /// Destroys the rows and the counts and gives their memory back.
    void release() {
        vector<vector<vector<double>>>().swap(values);
        vector<vector<double *>>().swap(rows);
        vector<double **>().swap(matrices);
        vector<double>().swap(resetTimers);
        vector<int>().swap(counts);
        vector<int>().swap(bufferedCounts);
    }
};

/// The trend and profile output a system writes: the trend matrices of every family
/// (production, gas line, and the wall temperatures of both), their sample counters and reset
/// times, the profile and single-cell output indices, the cells whose radial temperature
/// profiles are written, and the count of output passes.
struct TrendRecorder {
    /**
     * @brief Gas-line cells where radial temperature profiles are written.
     */
    vector<int> ncelperftransg;
    /**
     * @brief Times at which gas-line trend buffers are reset.
     */
    double *resettrendg = nullptr;
    /**
     * @brief Number of gas-line trend samples currently stored.
     */
    int *ntrendg = nullptr;
    /**
     * @brief Number of gas-line trend samples stored before the last flush.
     */
    int *ntrendgB = nullptr;
    /**
     * @brief Buffered gas-line wall-temperature trend data.
     */
    double ***MatTrendTransG = nullptr;
    /**
     * @brief Times at which gas-line wall-temperature buffers are reset.
     */
    double *resettrendtransg = nullptr;
    /**
     * @brief Number of gas-line wall-temperature samples currently stored.
     */
    int *ntrendtransg = nullptr;
    /**
     * @brief Number of gas-line wall-temperature samples stored before the last flush.
     */
    int *ntrendtransgB = nullptr;
    /**
     * @brief Production-line cells where radial temperature profiles are written.
     */
    vector<int> ncelperftransp;
    /**
     * @brief Times at which production-line trend buffers are reset.
     */
    double *resettrend = nullptr;
    /**
     * @brief Number of production-line trend samples stored before the last flush.
     */
    int *ntrendB = nullptr;
    /**
     * @brief Buffered production-line wall-temperature trend data.
     */
    double ***MatTrendTransP = nullptr;
    /**
     * @brief Times at which production-line wall-temperature buffers are reset.
     */
    double *resettrendtrans = nullptr;
    /**
     * @brief Number of production-line wall-temperature samples currently stored.
     */
    int *ntrendtrans = nullptr;
    /**
     * @brief Number of production-line wall-temperature samples stored before the last flush.
     */
    int *ntrendtransB = nullptr;
    /**
     * @brief Current production-profile output index.
     */
    int kontaTempoProf = 0;
    /**
     * @brief Current gas-line profile output index.
     */
    int kontaTempoProfG = 0;
    /**
     * @brief Current production-wall-temperature profile output index.
     */
    int kontaTempoTransProf = 0;
    /**
     * @brief Current gas-line wall-temperature profile output index.
     */
    int kontaTempoTransProfG = 0;
    /**
     * @brief For each single-cell output, the index of its next output time.
     */
    vector<int> kontaTempoCelUni;
    /**
     * @brief The production-line trends, read through MatTrendP, resettrend, ntrend and ntrendB.
     */
    TrendSet productionTrendSet;
    /**
     * @brief The gas-line trends, read through MatTrendG, resettrendg, ntrendg and ntrendgB.
     */
    TrendSet gasTrendSet;
    /**
     * @brief The production-line wall-temperature trends, read through MatTrendTransP,
     * resettrendtrans, ntrendtrans and ntrendtransB.
     */
    TrendSet productionWallTrendSet;
    /**
     * @brief The gas-line wall-temperature trends, read through MatTrendTransG,
     * resettrendtransg, ntrendtransg and ntrendtransgB.
     */
    TrendSet gasWallTrendSet;
    /**
     * @brief Buffered gas-line trend data.
     */
    double ***MatTrendG = nullptr;
    /**
     * @brief Buffered production-line trend data.
     */
    double ***MatTrendP = nullptr;
    /**
     * @brief Number of production-line trend samples currently stored.
     */
    int *ntrend = nullptr;
    /**
     * @brief Trend-output cycle counter; a value of one triggers header output.
     */
    double kimpT = 0.;
};

/// What a system keeps as a branch of a network: the inlet and outlet conditions of the
/// previous time level, the reverse-flow state, and the coupling with the other branch of a
/// parallel network.
struct NetworkCoupling {
    /**
     * @brief Inlet pressure stored at the previous time level.
     */
    double presEini = -1;
    /**
     * @brief Inlet temperature stored at the previous time level.
     */
    double tempEini = -1;
    /**
     * @brief Inlet gas mass fraction stored at the previous time level.
     */
    double titEini = -1;
    /**
     * @brief Inlet void fraction stored at the previous time level.
     */
    double alfEini = -1;
    /**
     * @brief Inlet complementary-liquid fraction stored at the previous time level.
     */
    double betaEini = -1;
    /**
     * @brief Previous-time-level value of presfim.
     */
    double presfimini = 0;
    /**
     * @brief Gas mass fraction imposed during reverse flow at the last control volume. Used only by transient
     * network simulations.
     */
    double titRev = 1.;
    /**
     * @brief Previous-time-level value of titRev.
     */
    double titRevini = 1.;
    /**
     * @brief Complementary-liquid fraction imposed during reverse flow at the last control volume. Used only by
     * transient network simulations.
     */
    double betaRev = 0.;
    /**
     * @brief Previous-time-level value of betaRev.
     */
    double betaRevini = 1.;
    /**
     * @brief Previous-time-level downstream or separator temperature.
     */
    double tGSupIni = 0.;
    /**
     * @brief Reserved temperature state; currently unused.
     */
    double tempSup = 0;
    /**
     * @brief Previous-time-level value of masChkSup.
     */
    int masChkSupini = 0;
    /**
     * @brief Previous-time-level value of mudaModoChk.
     */
    int mudaModoChkini = 0;
    /**
     * @brief Secondary-line index where parallel-network coupling with the primary line begins.
     */
    int SecPrimIniRedeP = 0;
    /**
     * @brief Secondary-line index where parallel-network coupling with the primary line ends.
     */
    int SecPrimFimRedeP = 0;
    /**
     * @brief Primary-line index aligned with SecPrimIniRedeP.
     */
    int PrimSecIniRedeP = 0;
    /**
     * @brief Primary-line index aligned with SecPrimFimRedeP.
     */
    int PrimSecFimRedeP = 0;
    /**
     * @brief Indicates thermal coupling on the primary branch of a parallel network.
     */
    int verificaAcopRedeP = 0;
    /**
     * @brief Indicates thermal coupling on the secondary branch of a parallel network.
     */
    int verificaAcopRedeS = 0;
    /**
     * @brief Initial source indices for parallel-network coupling.
     */
    vector<int> indFonteRedeParalelaIni;
    /**
     * @brief Initial production-liquid sources for the parallel network.
     */
    vector<double> fonteMpRedeParalelaIni;
    /**
     * @brief Initial complementary-liquid sources for the parallel network.
     */
    vector<double> fonteMcRedeParalelaIni;
    /**
     * @brief Initial gas sources for the parallel network.
     */
    vector<double> fonteMgRedeParalelaIni;
    /**
     * @brief Boundary-condition type applied to the secondary parallel-network branch.
     */
    int redeParalelaCCsecundario = -1;
    /**
     * @brief Primary branch index in a parallel network.
     */
    int redeParalelaP = -1;
    /**
     * @brief Secondary branch index in a parallel network.
     */
    int redeParalelaS = -1;
};

/// What the steady-state search carries from one iteration to the next: the convergence
/// monitor and its reference, the initial estimates, the reverse flow of a network solution
/// and its fluid and temperature, the number of production fluids and the slow thermal
/// coupling parameter.
struct SteadySearch {
    /**
     * @brief Control parameter for slowly varying thermal coupling.
     */
    double trocaTermicaLenta = 0.01;
    /**
     * @brief Number of production fluids.
     */
    int nfluP = 0;
    /**
     * @brief Initial holdup estimate used by the steady-state solver.
     */
    double chuteHol = -1.;
    /**
     * @brief Controls the search for an initial steady-state estimate.
     */
    int buscaIni = 0;
    /**
     * @brief Fluid state received from downstream during reverse network flow.
     */
    ProFlu fluiRevRede;
    /**
     * @brief Temperature associated with reverse network flow.
     */
    double tempRev = 0.;
    /**
     * @brief Indicates reverse flow in the steady-state network solution.
     */
    int revPerm = 0;
    /**
     * @brief Current steady-state convergence monitor.
     */
    double monitConvPerm = 1000.;
    /**
     * @brief Reference value for the steady-state convergence monitor.
     */
    double monitConvPermBase = 1.;
};

}  // namespace sisprod

/**
 * @brief Models and solves a one-dimensional production system.
 *
 * The class owns the production and gas-service-line state, fluid-property
 * tables, boundary conditions, network-coupling data, transient buffers,
 * and steady-state/transient solution procedures.
 */
class SProd {
    // See the note above the class: the adapters and updaters are friends so the
    // members only they read can be private.
    friend sisprod::gaslift::GasLiftState sisprod::adapters::gasLiftStateOf(SProd &);
    friend sisprod::steady::SteadyStateState sisprod::adapters::steadyStateOf(SProd &);
    friend sisprod::steady::SteadyStateSearchState sisprod::adapters::searchStateOf(SProd &);
    friend sisprod::transient::TransientStepState sisprod::adapters::transientStateOf(SProd &);
    friend sisprod::composition::CompositionState sisprod::adapters::compositionStateOf(SProd &);
    friend sisprod::transient::TransientSolveState sisprod::adapters::transientSolveStateOf(SProd &);
    friend sisprod::thermal::ThermalState sisprod::adapters::thermalStateOf(SProd &);
    friend sisprod::sources::SourceState sisprod::adapters::sourceStateOf(SProd &);
    friend driftflux::coefficient::ClosureState sisprod::adapters::closureStateOf(SProd &);
    friend trendoutput::TrendState sisprod::adapters::trendStateOf(const SProd &);
    friend struct sisprod::steady::SteadyStateUpdaters;
    friend struct sisprod::transient::TransientStepUpdaters;
    friend struct sisprod::transient::TransientSolveUpdaters;
    friend struct sisprod::gaslift::GasLiftTemperatureUpdater;
    friend struct sisprod::thermal::ThermalClosureUpdater;
    friend struct sisprod::thermal::ThermalSourceUpdater;
    friend struct sisprod::composition::CompositionUpdaters;

  public:
    /**
     * @brief Enables latent-heat calculations for black-oil simulations.
     */
    int CalcLat = 0;
  private:
  public:

    /**
     * @brief Section index when this object belongs to a pipeline network.
     */
    int indTramo = -1;
    /**
     * @brief Number of control volumes in the production line.
     */
    int ncel = 0;

  private:
  public:
    /**
     * @brief What the steady-state search carries from one iteration to the next.
     */
    sisprod::SteadySearch steadySearch;
  private:
  public:

    /**
     * @brief Current inlet pressure boundary condition.
     */
    double presE = -1;
    /**
     * @brief Current inlet temperature boundary condition.
     */
    double tempE = -1;
    /**
     * @brief Current inlet gas mass fraction.
     */
    double titE = -1;
    /**
     * @brief Computed inlet void fraction.
     */
    double alfE = -1;
    /**
     * @brief Computed inlet complementary-liquid fraction.
     */
    double betaE = -1;
    /**
     * @brief What this system keeps as a branch of a network.
     */
    sisprod::NetworkCoupling networkCoupling;

    /**
     * @brief Pressure in the last production-line control volume. When the surface choke is open, this value is
     * equal to the separator pressure.
     */
    double presfim = 0;
    /**
     * @brief Separator temperature or temperature at the inlet of the downstream section.
     */
    double tGSup = 0.;
    /**
     * @brief Production-line inlet temperature when no inlet-pressure boundary condition is imposed.
     */
    double temperatura = 0;

  private:
  public:


  private:
  public:
    /**
     * @brief Indicates whether the surface choke is active.
     */
    int masChkSup = 0;
  private:
  public:
    /**
     * @brief The fluid-property tables the fluids of this system's cells point at.
     */
    sisprod::PropertyTables tables;
  private:
  public:

    /**
     * @brief Parsed user input and simulation configuration.
     */
    Ler arq;
  private:
    /**
     * @brief Drift-flux correlation chosen for each flow regime.
     *
     * Resolved from arq once, when the object is built, so the per cell path
     * never reads configuration. Kept beside arq because that is its source.
     */
    driftflux::correlations::RegimeSelectors driftSelectors;
  public:
    /**
     * @brief Buffer used to write gas-line profiles.
     */
    FullMtx<double> flutG;
    /**
     * @brief Buffer used to write production-line profiles.
     */
    FullMtx<double> flut;
  private:
  public:
    /**
     * @brief The gas-lift line: its cells, chokes and valves, its pressure-velocity system, and
     * the unloading of its completion fluid.
     */
    sisprod::GasLiftLine gasLift;
  private:
  public:
    /**
     * @brief The state of a transient run.
     */
    sisprod::TransientRun transient;
  private:
  public:

    /**
     * @brief Current time step.
     */
    double dt = 0.;
  private:
  public:

  private:


  public:
  private:

  public:
  private:
  public:
  private:
  public:
    /**
     * @brief The trend and profile output this system writes.
     */
    sisprod::TrendRecorder trends;
  private:

  public:

  private:
  public:
    /**
     * @brief Event-log file name.
     */
    string tmpLog;
  private:

    /**
     * @brief Number of iterations used to bracket the initial steady-state root.
     */
    int iterperm = 0.;
  public:
    /**
     * @brief Steady-state mode flag: 1 while Num4Main solves a network branch at
     * steady state. The source terms read it (sisprod::sources::SourceState).
     */
    int modoPerm = 0.;
  private:
  public:
  private:
  public:

    /**
     * @brief Surface-choke model.
     */
    choke chokeSup;
  private:
  public:
    /**
     * @brief Multiphase production-line control volumes, held by productionCells.
     */
    Cel *celula = nullptr;

  private:
    /**
     * @brief The production-line control volumes celula points at.
     */
    vector<Cel> productionCells;
  public:


  private:
  public:

    /**
     * @brief Indicates that the section outlet is not connected to another network section.
     */
    int noextremo = 1;
    /**
     * @brief Indicates that the section inlet is not connected to another network section.
     */
    int noinicial = 1;
    /**
     * @brief Indicates that this section is a branch of a gas-lift ring network.
     */
    int derivaAnel = -1;


  private:
    /**
     * @brief Temporary-network flag; currently expected to remain zero.
     */
    int redeTemporario = 0;
  public:


  private:


  public:

  private:
  public:
    /**
     * @brief Section-blocking state.
     */
    int bloq = 0;

  private:
  public:
    /**
     * @brief Index corrections used by coupled thermal sections.
     */
    vector<int> acertaIndAcop;
    /**
     * @brief Shared one-dimensional simulation settings.
     */
    varGlob1D *vg1dSP = nullptr;
    /**
     * @brief Three-dimensional Poisson solver used by the thermal model.
     */
    solverP3D poisson3D;

  private:
    /**
     * @brief Indicates that the thermal source term is disabled.
     */
    int semTermo = 0;
  public:


  private:

    static constexpr const char *saidaTextoSis[16] = {"                          Post Coitum Omine Animal Triste Est                   ",
                                     "           'Ouca-me. O fim quase nunca esta longe, em nenhum momento!'          ",
                                     "      So nos curamos de um sofrimento depois de o haver suportado ate o fim.    ",
                                     "                   Infeliz e o espirito ansioso pelo futuro.                    ",
                                     "                                    Memento Mori                                ",
                                     " Somente um progresso calmo e constante, livre de precipitacao, conduz ao objetivo.",
                                     "             Paciencia, nove mulheres nao conseguem gerar uma crianca em um mes. ",
                                     "                  A necessidade e a mae da inovacao, mas a paciencia e o pai    ",
                                     "O sucesso nao e uma linha reta, e um jogo de resistencia, e cada tropeco e apenas um degrau a mais para a vitoria!",
                                     "                        Quem vive de navegar, o vento e quem lhe comanda                ",
                                     "    Uma vez me perguntaram o que achava da passagem do tempo, e eu disse: sou contra    ",
                                     "                 Nao importa o quanto voce va devagar, desde que nao pare                ",
                                     "Um simulador que resolve uma parada de producao, comeca avancando pequenos incrementos de tempo",
                                     "                            Nada e permanente, exceto a mudanca                           ",
                                     "                  Uma jornada de mil quilometros comeca com um unico passo                ",
									 "Seja paciente. Espere ate que a lama assente e a agua fique limpa. Permaneça imovel ate que a acao correta suja por si so"};
    static constexpr const char *saidaSubTextoSis[16] = {
        "                         Galeno de Pergamo do Transiente Longo                          ",
        "                     J. California Cooper depois da simulacao divergir                  ",
        "                                Marcel Proust no CrossFit                               ",
        "                              Seneca do Mindfulness                                     ",
        "                                   Zuleica da Funeraria                                 ",
        "                                       China In Box                                     ",
        "                                      Tiao do Linkedin                                  ",
        "                                      Marcao da Oficina                                 ",
        "                                    Mario Pascal do Insta                               ",
        "                          Seu Pereira na feira de artesanatos numericos                 ",
        "                        Luis Fernando Verissimo das Simulacoes Permanentes              ",
        "                             Confucio vendo a simulacao emperrar                       ",
        "                               Confucio das simulacoes sem fim                         ",
        "          Heraclito de Efeso vendo tudo mudar a cada incremento de tempo                ",
        "    Lao-Tse tomando coragem para simular um caso de parafinacao em dutos de producao    ",
	    "            Lao-Tse, vendo o incremento de tempo ficar cada vez menor                   "};
  public:

    /// Constructs and initializes a production-system simulation from the input and log files.
    SProd(string nomeArquivoEntrada, string nomeArquivoLog, tipoValidacaoJson_t validacaoJson,
          tipoSimulacao_t tipoSimulacao, varGlob1D *Vvg1dSP = 0, int TD = -1, int vbloq = 0,
          int temporario = 0,
          int reverso = 0,
          double *compfonte = 0,
          int *posicfonte = 0,
          int nfontes = 0,
          int redeperm = 1);
  private:
    /// The run state every construction and reassignment starts from.
    void resetRunState();
  public:
    /// Creates an empty production-system object.
    SProd();

    /// Not copyable: the object owns raw arrays, which operator= rebuilds
    /// instead of sharing.
    SProd(const SProd &) = delete;

  private:
    /// Gives back the memory of every array this object owns. operator= and the
    /// constructor from a parsed input call it before rebuilding, so that the arrays
    /// of the previous system are free before the input is read again.
    void releaseOwnedStorage();
  public:
    /// Rebuilds this system from sp's input file, read again from disk, and takes
    /// over sp's place in the network and reverse-flow state. Not a copy: the run
    /// state starts afresh and the system is assembled again.
    SProd &operator=(const SProd &);

    /// What a system built from an already parsed input takes from the system that
    /// input came from: its place in the network and its reverse-flow state.
    struct CarriedState {
        int noextremo;
        int noinicial;
        int derivaAnel;
        int bloq;
        double betaRevini;
        double titRevini;
        double dtCicMin;
    };

    /// The state another system built from this one's input takes over.
    CarriedState carriedState() const;

    /// Builds a system from an already parsed input, without reading the JSON file
    /// again, with the place in the network and the reverse-flow state in carried.
    SProd(Ler &parsedInput, const CarriedState &carried);

  private:
    /// Points every cell fluid, and every source fluid it carries, at the
    /// bubble-point tables read from the PVTSim file, and switches them to
    /// saturation model 4.
    void assignPvtSimBubbleTablesToCells();
    /// Reads the bubble-point curve from the PVTSim file, points the cell fluids at
    /// it and writes perfilBolha; with tabRSPB on, also reads the solution gas-oil
    /// ratio table and writes perfilRSLivia.
    void loadPvtSimSaturationTables();
    /// Builds the bubble-point curve and the solution gas-oil ratio table from the
    /// fluid correlations over the input table's pressure-temperature grid, writes
    /// perfilBolha and perfilRSLivia, and points the cell fluids at both.
    void generateSaturationTablesFromCorrelations();
    /// Copies the run configuration into the members, checks that an injection well
    /// has the IPR its boundary condition needs, builds the Cp and JTL tables,
    /// hands the fluid constants to every fluid, generates the pipe and the
    /// production cells, and places every source -- including the extra gas sources
    /// at compfonte, whose cell indices it records in posicfonte.
    void buildProductionCells(double *compfonte, int *posicfonte, int nfontes);
    /// Gives the inlet the sources its boundary condition needs (and the second
    /// cell, under a blockage), then places the accessories -- pumps, volumetric
    /// pumps, pressure-drop requirements, heat sources, the master valve and the
    /// other valves -- and sets up the outlet pressure, the surface and injection
    /// chokes and the pigs.
    void configureInletSourcesAndAccessories(int nfontes);
    /// Builds the gas-lift line when there is one: its cells, the second master
    /// valve, the injection choke, one gas-lift valve choke per valve with its
    /// position on both lines, and each gas cell's share of the annulus above the
    /// discharge cell.
    void buildGasLiftLine();
    /// Rejects or warns about source and boundary-condition combinations the run
    /// cannot honour (RN-300, RN-301, gas-lift discharge without an IPR, no outlet
    /// pressure, accessories in the last two cells), builds the gas-lift discharge
    /// hydrostatics and the event log, applies the initial state of a production
    /// well, and records the surface temperature and mass flow it starts from.
    void validateSetupAndApplyInitialState();
    /// With the dynamic property table on, splits the pipe into table segments that
    /// end at every source cell.
    void buildDynamicTablesAndInclinations();
    /// Sets the latent-heat switch from the input and, when it is on and
    /// flashCompleto is 0, reads the latent-heat table from the PVTSim file into
    /// HLat and writes perfilLatente.
    void configureLatentHeat();
    /// Sets the gas-density correction factors of every fluid in every cell -- the
    /// cell's, its source's and its reservoir cells' -- to 1 when the correction is
    /// off, or evaluates them at local pressure and temperature when it is on; then
    /// makes the first cell's fluid the fluid of its source.
    void applyDensityCorrectionsAndInletFluid();
    /// Lists the cells with two-dimensional heat diffusion, copies the master-valve
    /// opening and closing times, maps the transient profile positions to global
    /// thermal-node indices, and -- unless the network is only provisional --
    /// allocates the trend matrices of the production line, the gas line and the
    /// transient trends, sized to the longest simulated time and filled with the
    /// -10000 sentinel.
    void allocateEventProfileAndTrendArrays();
    /// Resets the column-annulus and network coupling flags, locates the
    /// column-annulus coupling range on the gas line, sets the steady and transient
    /// profile counters and their first output times, opens the event log with the
    /// events known at start, finds the smallest cell length, and zeroes the
    /// moving-average and running-total state the transient loop starts from.
    void resetCouplingAndOutputState();
    /// Builds the production section after input parsing.
    void montasistema(double *compfonte = 0,
                      int *posicfonte = 0,
                      int nfontes = 0);

    /// Initializes the gas service line for gas-lift unloading.
    void HidroDescargaG();
    /// Initializes the production line for gas-lift unloading.
    void HidroDescargaP();

    /// Controls injection and upstream-choke pressures from gas-lift-valve flow rates.
    double BuscaPresInjDesc();

    /// Calculates gas-lift-valve opening area from calibration and operating conditions.
    double areaValvCali(double PCal, double TCal, double PVO, double PT,
                        double dextern, double areagarg, double Rvalv, double Temp);
    /// Advances the temperature of one gas-line control volume.
    void calctempGas(int i, double tempantiga, int modoPerm = 0);
    /// Updates gas-line temperature in the completion-fluid region during unloading.
    void tempDescarga(int i);
    /// Calculates gas temperature across a gas-lift valve using the Joule-Thomson model.
    double TempDescGL(int igl);
    /// Advances the coupled gas-line pressure, velocity, and temperature solution.
    void subtempoGas();
    /// Exchanges heat-transfer data between the production column and annulus.
    void conectaColuna();
    /// Advances the temperature of one production-line control volume.
    void calctemp(int i, double tempantiga, int modoPerm = 0);

    /// Updates IPR, gas, liquid, leak, and gas-lift source terms.
    void renovaFonte(int ind);
    /// Stores previous void fractions and updates pig motion and reception.
    void renovaalbetini();
    /// Caches cell and face densities to avoid repeated property calculations.
    void renovaMasEsp();

    /// Selects and evaluates the slip correlation at a production-line face.
    void CalcC0Ud(int ind, double &c0, double &ud);
    /// Evaluates slip parameters for the intermediate network state.
    void CalcC0UdBuf(int ind, double &c0, double &ud);
    /// Evaluates slip parameters at the inlet of an internal network section.
    void CalcC0UdIni(int ind, double &c0, double &ud);
    /// Evaluates slip parameters at the inlet of an internal network section.
    void CalcC0UdIniBuf(int ind, double &c0, double &ud);
  public:

    /// Loads pressure-velocity results into cell and face state variables.
    void renova(int expli = 0);
  private:
  public:
    /// Updates only section-end states during intermediate network convergence.
    void renovaBuffer();
    /// Copies previous states when an active choke bypasses the intermediate network solve.
    void renovaBufferCego();

    /// Updates distributed mass-transfer terms used by void-fraction and mixture-mass equations.
    void renovaTemp();
  private:
    /// Evaluates wax deposition and its effects.
    void avaliaParafina();
    /// Transports black-oil properties such as GOR, API, gas density, and CO2 fraction.
    void renovaRGOdgYco2(ProFlu fluiRev = ProFlu());
    /// Applies the alternate compositional molar-fraction transport update.
    void renovaFracMol2(ProFlu fluiRev = ProFlu());
  public:
    /// Computes T1 and T2 used to split mixture mass flow into liquid and gas flows.
    void renovaterm(int aflu = 0);
    /// Computes T1 and T2 at the outlet of an internal network section.
    void renovatermAfluFim();
    /// Computes T1 and T2 at the inlet of an internal network section.
    void renovatermColIni();

    /// Calculates outlet-choke flow from the last-cell and separator pressures.
    void calcCCpres(double titRev = 1., double alfRev = 1., double betRev = 0.);
    /// Calculates outlet-choke flow for the intermediate network state.
    void calcCCBuffer(double titRev = 1., double alfRev = 1., double betRev = 0.);

    /// Selects a stable time step from CFL and additional model restrictions.
    void determinaDT(int vexpli = 0);
  private:
  public:
    /// Limits time-step growth when accepted steps remain well below the CFL estimate.
    void atenuaDtMax();
    /// Checks whether liquid-density variation requires the complete formulation.
    void avaliaVariaDpDt(double razMast = 0, double razMast0 = 0, int vexpli = 0);
    /// Updates the valve logic for the active-state formulation.
    void aberturaVal1();
    /// Updates the valve logic for the inactive-state formulation.
    void aberturaVal0();
    /// Restricts the time step during valve-state transitions.
    void restringeDTporValv();

    /// Advances the gas-injection system and prepares its coupling with the production line.
    void solveLinGas();
    /// Advances phase volume fractions.
    void EvoluiFrac(double alfrev = 1., double betrev = 0., int ciclo = 0);
    /// Restores the initial fraction state after an invalid update.
    void ReiniEvolFrac0();
  private:
  public:
    /// Restores previous volume fractions after a nonphysical update.
    void ReiniEvolFrac();
    /// Updates phase fractions in cells affected by a moving pig.
    void AtualizaPig();
    /// Assembles and solves the global pressure-velocity coupling system.
    void SolveAcopPV(int vexpli = 0, int ciclo = 0);

  private:
  public:
    /// Advances the transient energy equation.
    void marchaEnergTrans(int ciclo = 0, int ciclomax = 0);
    /// Refreshes local dynamic fluid-property tables.
    void atualizaMiniTab();
    /// Updates the first boundary condition.
    void atualizaCC1();

  private:
    /// Runs the hydrate-envelope solvers for the production and gas lines.
    void solveHydrateEnvelopes();
  public:
    /// Advances the complete transient production-system solution and handles output and logging.
    void SolveTrans(double titRev = 1., double alfRev = 1., double betRev = 0.,
                    int nrede = -1, ProFlu fluiRev = ProFlu());

    /// Appends one production-line trend sample to its output buffer.
    void ImprimeTrendP(int i, int nrede = -1);
    /// Appends one production-line trend sample to its output buffer.
    void ImprimeTrendPCab(int i, int nrede = -1);
    /// Appends one gas-line trend sample to its output buffer.
    void ImprimeTrendG(int i, int nrede = -1);
    /// Appends one gas-line trend sample to its output buffer.
    void ImprimeTrendGCab(int i, int nrede = -1);
  private:
    /// Appends one production-line wall-temperature trend sample.
    void ImprimeTrendTransP(int i);
    /// Appends one production-line wall-temperature trend sample.
    void ImprimeTrendTransPCab(int i);
    /// Appends one gas-line wall-temperature trend sample.
    void ImprimeTrendTransG(int i);
    /// Appends one gas-line wall-temperature trend sample.
    void ImprimeTrendTransGCab(int i);
  public:

    /// Marches the steady production solution using a bottomhole-pressure guess with outlet pressure prescribed.
    double marchaProdPerm1(double pchute);
  private:
  public:
    /// Brackets and solves the bottomhole-pressure root for marchaProdPerm1.
    double buscaProdPfundoPerm(double chute = -1., int kontaTenta = -1);
    /// Brackets and solves the reverse-flow bottomhole-pressure root.
    double buscaProdPfundoPermRev(double chute = -1.);
    /// Brackets and solves the bottomhole-pressure root for marchaProdPerm2.
    double buscaProdPfundoPerm2(double chute = -1., int kontaTenta = -1);
    /// Runs a direct march when inlet pressure and flow rate are known.
    double buscaProdPfundoPerm3(double pentrada);
    /// Marches the steady production solution using a bottomhole mass-flow guess.
    double marchaProdPresPres1(double mchute);
  private:
  public:
    /// Brackets and solves the mass-flow root for marchaProdPresPres1.
    double buscaProdPresPresPerm(double mchute, double maxvaz = 0., int kontaiter = 0);
    /// Brackets and solves the reverse-flow mass-flow root.
    double buscaProdPresPresPermRev(double mchute, double maxvaz = 0., int kontaiter = 0);
  private:
  public:
    /// Brackets and solves the mass-flow root for marchaProdPresPres2.
    double buscaProdPresPresPerm2(double mchute, double maxvaz = 0.);
  private:
  public:
    /// Brackets and solves the mass-flow root for the closed-choke case, chosen
    /// by Num4Main when arq.chokep.abertura[0] <= 1e-15.
    ///
    /// It reaches marchaProdPresPres2, not marchaProdPresPres3 -- via zriddr(.., 2,
    /// 1) and multMarcha -- which is correct by the reduction described above.
    double buscaProdPresPresPerm3(double mchute, double maxvaz = 0.);

  private:
    /// Brackets and solves the pressure root for marchaGasPerm2.
    double buscaGasPresPerm2();
    /// Brackets and solves the pressure root for marchaGasPerm3.
    double buscaGasPresPerm3();
    /// Corrects gas density in one control volume.
    void corrDeng(int i);
    /// Calculates steady-state slip parameters at a downstream face.
    void CalcC0UdPerm(int ind, double &c0, double &ud);
    /// Marches steady-state temperature from cell i-1 to cell i.
    void RenovaTempPerm(int i, int RK);
    /// Marches steady-state temperature in the reverse direction.
    void RenovaTempPermRev(int i, int RK);
    /// Synchronizes neighboring face temperatures after updating a cell-center temperature.
    void atualizaPeriTempProd(int i);

    /// Marches steady-state pressure along the gas service line.
    void RenovaPresGasPerm(int i);
    /// Estimates gas-line pressure variation without advancing the full pressure march.
    double delpGasPerm(int i);
    /// Estimates pressure variation in an injection-well system.
    double delpInjPerm(int i);
    /// Calculates gas-lift-valve flow from gas-line and production-line pressures.
    void calcVazGasPerm(int i);
    /// Initializes gas-lift-valve flow estimates before the steady-state march.
    void IniciaVazValvGasPerm(int i);
    /// Marches steady-state gas-line temperature from cell i-1 to cell i.
    void RenovaTempGasPerm(int i);

    /// Exchanges steady-state heat-transfer data between column and annulus.
    void conectaColunaPerm();
    /// Initializes estimated column-annulus heat transfer before the steady-state march.
    void IniciaconectaColunaPerm();

  public:
    /// Solves injection cases CC1 and CC3.
    double buscaInjPfundoPerm1(double chute = -1.);
    /// Solves injection case CC0.
    double buscaInjPfundoPerm2(double chute = -1.);
    /// Solves injection case CC2.
    double buscaInjPfundoPerm3(double chute = -1.);
    /// Solves injection case CC4.
    double buscaInjPfundoPerm4();
    /// Solves injection case CC5.
    double buscaInjPfundoPerm5(double chute = -1.);

  private:

  public:

    /// Estimates steady production-network node pressures from hydrostatics.
    double hidroreverso(double hol, double vaz = 0, double vazG = 0);
    /// Estimates steady injection-network node pressures from hydrostatics.
    double hidroreversoInj(double hol, double vaz = 0);
    /// Estimates secondary-branch pressures in a parallel production network.
    double hidroTramoSecundario(double titulo);
  private:

    /// Marches the steady secondary branch of a parallel network.
    double marchaTramoSecVaz(double pchute, double chutemass = -1);
  public:

    /// Calculates temperature downstream from the surface choke for network coupling.
    void calcTempFim();

  private:

    /**
     * @brief Caches the drift-flux correlation of each regime from arq.
     *
     * Called wherever arq is built or replaced, so driftSelectors never goes
     * stale.
     */
    void resolveDriftSelectors();
};

#endif /* SISPROD_H_ */
