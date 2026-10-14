/*
 * SisProd.h
 *
 * Created on: December 21, 2016
 *     Author: Eduardo
 */

#ifndef SISPROD_H_
#define SISPROD_H_

#include "Acidentes2.h"
#include "DriftFluxClosure.h"
#include "Bcsm2.h"
#include "BombaVol.h"
#include "FonteMas.h"
#include "FonteMassCHK.h"
#include "Geometria.h"
#include "GradientCorrelations.h"
#include "Leitura.h"
#include "Matriz.h"
#include "PropFlu.h"
#include "PropFluCol.h"
#include "RuntimeDiagnostics.h"
#include "SisProdConstants.h"
#include "Vetor.h"
#include "celula3.h"
#include "celulaGas.h"
#include "chokegas.h"
#include "dados3DPoisson.h"
#include "estrat.h"
#include "estruturaTabDin.h"
#include "mapa.h"
#include "solver3DPoisson.h"
#include "variaveisGlobais1D.h"
#include <ctime>
#include <math.h>
#include <string>
#include <vector>

class SProd;
namespace trendoutput { struct TrendState; }

namespace sisprod {

class SolveContext;
namespace composition { struct CompositionState; struct CompositionUpdaters; }
namespace gaslift { struct GasLiftState; struct GasLiftTemperatureUpdater; }
namespace sources { struct SourceState; }
namespace steady { struct SteadyStateSearchState; struct SteadyStateState; struct SteadyStateUpdaters; }
namespace thermal { struct ThermalClosureUpdater; struct ThermalSourceUpdater; struct ThermalState; }
namespace transient {
struct TransientSolveState;
struct TransientSolveUpdaters;
struct TransientStepState;
struct TransientStepUpdaters;
}  // namespace transient

namespace adapters {
gaslift::GasLiftState gasLiftStateOf(SProd &system, SolveContext &context);
steady::SteadyStateState steadyStateOf(SProd &system, SolveContext &context);
steady::SteadyStateSearchState searchStateOf(SProd &system, SolveContext &context);
transient::TransientStepState transientStateOf(SProd &system, SolveContext &context);
composition::CompositionState compositionStateOf(SProd &system, SolveContext &context);
transient::TransientSolveState transientSolveStateOf(SProd &system, SolveContext &context);
thermal::ThermalState thermalStateOf(SProd &system, SolveContext &context);
sources::SourceState sourceStateOf(SProd &system);
driftflux::coefficient::ClosureState closureStateOf(SProd &system);
trendoutput::TrendState trendStateOf(const SProd &system);
}  // namespace adapters

/// A pressure-temperature table the fluids read through a double **: its rows, and a pointer to each.
struct GridTable {
    std::vector<std::vector<double>> rows;
    std::vector<double *> rowPointers;

    /// rowCount rows of columnCount values each.
    void allocate(int rowCount, int columnCount) {
        rows.assign(rowCount, std::vector<double>(columnCount));
        rowPointers.resize(rowCount);
        for (int i = 0; i < rowCount; i++)
            rowPointers[i] = rows[i].data();
    }
    /// Destroys the rows and gives their memory back.
    void release() {
        std::vector<std::vector<double>>().swap(rows);
        std::vector<double *>().swap(rowPointers);
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
    /// Black-oil gas-compressibility-factor table.
    double **zdranP = nullptr;
    /// Pressure derivative of the black-oil compressibility-factor table.
    double **dzdpP = nullptr;
    /// Temperature derivative of the black-oil compressibility-factor table.
    double **dzdtP = nullptr;
    /// Black-oil gas specific-heat table.
    double **cpg = nullptr;
    /// Black-oil produced-liquid specific-heat table.
    double **cpl = nullptr;
    /// Temperature derivative of liquid density.
    double **drholdT = nullptr;
    /// Black-oil latent-heat table.
    GridTable HLat;
    /// Bubble-pressure values imported from PVTSim for the Livia solution-gas-ratio correlation. This table is
    /// separate from the full PVTSim fluid-property model.
    std::vector<double> PBPVTSim;
    /// Bubble-temperature values imported from PVTSim for the Livia solution-gas-ratio correlation. This table is
    /// separate from the full PVTSim fluid-property model.
    std::vector<double> TBPVTSim;
    /// Precomputed solution-gas-ratio table for the Livia correlation, used to avoid repeating its expensive
    /// calculation during the simulation.
    GridTable RSLivia;
    /// Dynamic fluid-property tables.
    std::vector<tabelaDinamica> tabDin;
    /// Number of dynamic property tables.
    int ntabDin = 0;
};

/// The gas-lift line: its cells, the injection choke, the valves and the cells they join, its
/// pressure-velocity system, its surface and initial conditions, the column-annulus thermal
/// coupling, the completion-fluid interface, and the unloading controller's samples and
/// settings.
struct GasLiftLine {
    /// Global pressure-velocity coupling matrix for the gas line, sized for one cell until the system is sized for its
    /// input.
    BandMtx<double> matglobG{kGasLineUnknownsPerCell, kGasLineBandwidth, kGasLineBandwidth};
    /// Right-hand side and solution vector for the gas-line pressure-velocity system, sized like matglobG.
    Vcr<double> termolivreG{kGasLineUnknownsPerCell};
    /// Production-line cell indices associated with gas-lift valves.
    std::vector<int> posicVGLP;
    /// Service-line cell indices associated with gas-lift valves.
    std::vector<int> posicVGLG;
    /// Service-line index where column-annulus thermal coupling begins.
    int AnulaColunaIni = 0;
    /// Service-line index where column-annulus thermal coupling ends.
    int AnulaColunaFim = 0;
    /// Production-column index aligned with AnulaColunaIni.
    int ColunaAnulaIni = 0;
    /// Production-column index aligned with AnulaColunaFim.
    int ColunaAnulaFim = 0;
    /// Indicates whether column-annulus thermal coupling is enabled.
    int verificaAcop = 0;
    /// Gas-lift valves installed in the system.
    std::vector<ChokeGas> chokeVGL;
    /// Time horizon used by the gas-lift unloading controller.
    double tempMedContDesc = kUnloadingAveragingWindow;
    /// Maximum number of flow samples retained by the unloading PI controller.
    double maxVecContDesc = kUnloadingSampleCapacity;
    /// Average maximum completion-fluid mass flow through the gas-lift valves.
    double vazmedDesc = 0;
    /// Averaging interval used for gas-lift-valve mass flow.
    double tempmedDEsc = 0;
    /// Maximum valve mass-flow samples used by the unloading controller.
    std::vector<double> vazmaxMedDesc;
    /// Time-step samples associated with the unloading-flow average.
    std::vector<double> dtDesc;
    /// Separator pressure or pressure at the inlet of the downstream network section.
    double pGSup = 0;
    /// Number of control volumes in the gas service line.
    int ncelGas = 0;
    /// Gas-injection pressure. Zero until an injection-pressure condition or a sensitivity analysis case sets it; the
    /// steady summary and the restart file record it either way.
    double presiniG = 0.;
    /// Gas-injection temperature in the service line. Zero until an injection condition or a sensitivity analysis
    /// case sets it; the restart file records it either way.
    double tempiniG = 0.;
    /// Gas-injection choke model.
    ChokeGas chokeInj;
    /// Gas service-line control volumes, held by gasCells.
    CelG *celulaG = nullptr;
    /// The gas service-line control volumes celulaG points at.
    std::vector<CelG> gasCells;
    /// Current service-line cell containing the completion-fluid/gas interface.
    int celInter = kNoInterfaceCell;
    /// Maximum time step that keeps the unloading interface within one gas-line cell.
    double dtInter = 0.;
    /// Current completion-fluid/gas interface velocity.
    double velInter = 0.;
    /// Previous-time-level interface cell.
    int celInterIni = 0;
    /// Previous-time-level interface-limited time step.
    double dtInterIni = 0.;
    /// Previous-time-level interface velocity.
    double velInterIni = 0.;
};

/// The state of a transient run: the time step, its history and the restrictions on it, the
/// moving averages and running totals, the Master1 valve state and schedule, the
/// production-line pressure-velocity system, the buffered sources, the event and log
/// counters, the transport switches, the pigs and the cells of the two-dimensional Poisson
/// model.
struct TransientRun {
    /// Enables transport equations for primitive black-oil properties, including API gravity, BSW, gas-oil ratio, and
    /// light/heavy mass fractions.
    int trackRGO = 0;
    /// Enables transport equations for gas density and the gas-phase CO2 molar fraction.
    int trackDeng = 0;
    /// Current Master1 valve state.
    int EstadoMaster1 = 0;
    /// Counter used while changing the Master1 state.
    int contaMaster1 = 0;
    /// Number of pigs scheduled for launch.
    int npig = 0;
    /// Cell indices where pigs are received.
    std::vector<int> receb;
    /// Global pressure-velocity coupling matrix for the multiphase production line, sized for one cell until the
    /// system is sized for its input.
    BandMtx<double> matglobP{kProductionUnknownsPerCell, kProductionLeftBandwidth, kProductionRightBandwidth};
    /// Number of scheduled Master1 opening events.
    int nabreM1 = 0;
    /// Number of scheduled Master1 closing events.
    int nfechaM1 = 0;
    /// Times at which Master1 closes.
    std::vector<double> fechaM1;
    /// Times at which Master1 opens.
    std::vector<double> abreM1;
    /// Number of events written to the event log.
    int contaLog = 0;
    /// Smallest production-line control-volume length.
    double menorDx = 0.;
    /// Simulation time-step counter.
    int kSP = 0;
    /// Controls event-log output frequency.
    int KontaImprime = 0;
    /// Index of the next scheduled simulation event.
    int indevento = 0;
    /// History of recently accepted time steps.
    std::vector<double> dtSim;
    /// History of time steps proposed by the CFL criterion.
    std::vector<double> dtCFL;
    /// Average time step proposed by the CFL criterion.
    double dtCFLMed = 1.;
    /// Average time step actually used by the simulation.
    double dtSimMed = 1.;
    /// Indicates that time-step growth must remain restricted.
    int restriDt = 0;
    /// Number of remaining steps under the current time-step restriction.
    int kontarestriDt = 0;
    /// Accumulated CFL time steps used to compute dtCFLMed.
    double dtCFLTotal = 0.;
    /// Accumulated accepted time steps used to compute dtSimMed.
    double dtSimTotal = 0.;
    /// Counts alternating liquid-flow oscillations near an active surface choke.
    int kontaGolfada = kInitialSlugCount;
    /// Multiplier applied to a small artificial upstream gas flow during difficult Master1 closures. It mitigates a
    /// pressure blind spot when the upstream side contains only liquid and increases if the problem persists.
    double momentoDesesp = 0;
    /// The artificial gas the closing-valve control added in the current time step, kept until the step is accepted,
    /// when it goes to the run's diagnostics.
    std::vector<diagnostics::ArtificialGasAddition> artificialGasThisStep;
    /// Recent maximum pressure-change rates.
    std::vector<double> taxaDpMax;
    /// Average maximum pressure-change rate.
    double DpMaxMed = 1.;
    /// Recent maximum temperature-change rates.
    std::vector<double> taxaDTMax;
    /// Average maximum temperature-change rate.
    double DTMaxMed = 1.;
    /// Counter controlling compositional-property refreshes.
    int kontaRenovaComp = 0;
    /// Production cells handled by the two-dimensional Poisson model.
    std::vector<int> indCelPoisson2D;
    /// Number of cells handled by the two-dimensional Poisson model.
    int nCelulaPoisson2D = 0;
    /// Signals that the time step must be adjusted.
    int alteraTempo = 0;
    /// Requests rollback and time-step reevaluation when at least one control volume produces a holdup or volume
    /// fraction outside the physical [0, 1] range.
    int reinicia = 0;
    /// Recent Master1 ratio values for the active state.
    double vRazMast1[kMasterRatioHistorySize] = {};
    /// Recent Master1 ratio values for the inactive state.
    double vRazMast0[kMasterRatioHistorySize] = {};
    /// Critical Master1 ratio history used by the switching logic.
    double vRazMastCrit[kMasterRatioHistorySize] = {
        kDefaultMasterCriticalRatio, kDefaultMasterCriticalRatio, kDefaultMasterCriticalRatio,
        kDefaultMasterCriticalRatio, kDefaultMasterCriticalRatio, kDefaultMasterCriticalRatio,
        kDefaultMasterCriticalRatio, kDefaultMasterCriticalRatio, kDefaultMasterCriticalRatio,
        kDefaultMasterCriticalRatio};
    /// Previous-time-level value of pGSup.
    double pGSupIni = 0.;
    /// Recorded in the restart file; nothing in the simulation reads it.
    double massfonte = 0.;
    /// CFL safety factor, typically set to 0.8.
    double mult = 0;
    /// Time-averaged pressure in the final production-line control volume. Used to decide whether the surface choke
    /// behaves as a localized pressure loss or as a discharge-flow model.
    double presMedMov = 0;
    /// Time-averaged mixture volumetric flux in the final production-line control volume. Used by the surface-choke
    /// operating-mode logic.
    double jMedMov = 0;
    /// Time-averaged void fraction in the final production-line control volume, updated once the averaging window
    /// fills. Zero until then; only the restart file reads it.
    double alfMedMov = 0.;
    /// Start time of the moving-average window.
    double tMedMov = 0;
    /// Duration of the moving-average window.
    double ktMedMov = 0.;
    /// Accumulated pressure used to compute presMedMov.
    double pTotal = 0.;
    /// Accumulated mixture flux used to compute jMedMov.
    double jTotal = 0.;
    /// Accumulated void fraction used to compute alfMedMov.
    double alfTotal = 0.;
    /// Pressure samples used by the moving-average calculation.
    std::vector<double> presVet;
    /// Mixture-flux samples used by the moving-average calculation.
    std::vector<double> jVet;
    /// Void-fraction samples used by the moving-average calculation.
    std::vector<double> alfVet;
    /// Time samples associated with the moving-average window.
    std::vector<double> tVet;
    /// Current surface-choke open/closed state.
    int aberto = 0;
    /// Previous-time-level surface-choke state.
    int abertoini = 0;
    /// Counter that delays transitions out of active-choke mode.
    int tempoaberto = 0;
    /// Previous-time-level value of tempoaberto.
    int tempoabertoini = 0;
    /// Signals a surface-choke operating-mode transition.
    int mudaModoChk = 0;
    /// Selects the interphase mass-transfer model: 0 = complete, 1 = fully explicit, 2 = simplified, and 3 =
    /// disabled.
    int TransMassModel = 0;
    /// Number of pigs currently moving through the line.
    int indpigP = 0;
    /// Previous-time-level value of indpigP.
    int indpigPini = 0;
    /// Right-hand side and solution vector for the production-line pressure-velocity system, sized like matglobP.
    Vcr<double> termolivreP{kProductionUnknownsPerCell};
    /// Simulation end time.
    double tfinal = 0;
    /// Current complete-model activation state.
    int modeloCompleto = 1;
    /// Auxiliary CFL time-step accumulator.
    double dtauxCFL = 0.;
    /// Auxiliary accepted-time-step accumulator.
    double dtauxFinal = 0.;
    /// Intermediate production-liquid inflow estimate for a network section.
    double fontemassPRBuf = 0.;
    /// Intermediate complementary-liquid inflow estimate for a network section.
    double fontemassCRBuf = 0.;
    /// Intermediate gas inflow estimate for a network section.
    double fontemassGRBuf = 0.;
    /// Minimum time step allowed by the current cycle. Zero until a transient cycle sets it.
    double dtCicMin = 0.;
};

/// The storage of one family of trends: for each trend, its rows of samples, and its reset time
/// and sample counts. The solve reads it through the pointers allocate() points at it.
struct TrendSet {
    std::vector<std::vector<std::vector<double>>> values;
    std::vector<std::vector<double *>> rows;
    std::vector<double **> matrices;
    std::vector<double> resetTimers;
    std::vector<int> counts;
    std::vector<int> bufferedCounts;

    /// For each of the count trends, length(i) rows of rowWidth(i) values whose first
    /// sentinelCount(i) hold kTrendSentinel, and the reset timers and sample counts,
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
            values[i].assign(n, std::vector<double>(rowWidth(i)));
            rows[i].resize(n);
            for (int j = 0; j < n; j++) {
                rows[i][j] = values[i][j].data();
                for (int k = 0; k < sentinelCount(i); k++)
                    values[i][j][k] = kTrendSentinel;
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
        std::vector<std::vector<std::vector<double>>>().swap(values);
        std::vector<std::vector<double *>>().swap(rows);
        std::vector<double **>().swap(matrices);
        std::vector<double>().swap(resetTimers);
        std::vector<int>().swap(counts);
        std::vector<int>().swap(bufferedCounts);
    }
};

/// The trend and profile output a system writes: the trend matrices of every family
/// (production, gas line, and the wall temperatures of both), their sample counters and reset
/// times, the profile and single-cell output indices, the cells whose radial temperature
/// profiles are written, and the count of output passes.
struct TrendRecorder {
    /// Gas-line cells where radial temperature profiles are written.
    std::vector<int> ncelperftransg;
    /// Times at which gas-line trend buffers are reset.
    double *resettrendg = nullptr;
    /// Number of gas-line trend samples currently stored.
    int *ntrendg = nullptr;
    /// Number of gas-line trend samples stored before the last flush.
    int *ntrendgB = nullptr;
    /// Buffered gas-line wall-temperature trend data.
    double ***MatTrendTransG = nullptr;
    /// Times at which gas-line wall-temperature buffers are reset.
    double *resettrendtransg = nullptr;
    /// Number of gas-line wall-temperature samples currently stored.
    int *ntrendtransg = nullptr;
    /// Number of gas-line wall-temperature samples stored before the last flush.
    int *ntrendtransgB = nullptr;
    /// Production-line cells where radial temperature profiles are written.
    std::vector<int> ncelperftransp;
    /// Times at which production-line trend buffers are reset.
    double *resettrend = nullptr;
    /// Number of production-line trend samples stored before the last flush.
    int *ntrendB = nullptr;
    /// Buffered production-line wall-temperature trend data.
    double ***MatTrendTransP = nullptr;
    /// Times at which production-line wall-temperature buffers are reset.
    double *resettrendtrans = nullptr;
    /// Number of production-line wall-temperature samples currently stored.
    int *ntrendtrans = nullptr;
    /// Number of production-line wall-temperature samples stored before the last flush.
    int *ntrendtransB = nullptr;
    /// Current production-profile output index.
    int kontaTempoProf = 0;
    /// Current gas-line profile output index.
    int kontaTempoProfG = 0;
    /// Current production-wall-temperature profile output index.
    int kontaTempoTransProf = 0;
    /// Current gas-line wall-temperature profile output index.
    int kontaTempoTransProfG = 0;
    /// For each single-cell output, the index of its next output time.
    std::vector<int> kontaTempoCelUni;
    /// The production-line trends, read through MatTrendP, resettrend, ntrend and ntrendB.
    TrendSet productionTrendSet;
    /// The gas-line trends, read through MatTrendG, resettrendg, ntrendg and ntrendgB.
    TrendSet gasTrendSet;
    /// The production-line wall-temperature trends, read through MatTrendTransP, resettrendtrans, ntrendtrans and
    /// ntrendtransB.
    TrendSet productionWallTrendSet;
    /// The gas-line wall-temperature trends, read through MatTrendTransG, resettrendtransg, ntrendtransg and
    /// ntrendtransgB.
    TrendSet gasWallTrendSet;
    /// Buffered gas-line trend data.
    double ***MatTrendG = nullptr;
    /// Buffered production-line trend data.
    double ***MatTrendP = nullptr;
    /// Number of production-line trend samples currently stored.
    int *ntrend = nullptr;
    /// Trend-output cycle counter; a value of one triggers header output.
    double kimpT = 0.;
};

/// What a system keeps as a branch of a network: the inlet and outlet conditions of the
/// previous time level, the reverse-flow state, and the coupling with the other branch of a
/// parallel network.
struct NetworkCoupling {
    /// Inlet pressure stored at the previous time level.
    double presEini = -1;
    /// Inlet temperature stored at the previous time level.
    double tempEini = -1;
    /// Inlet gas mass fraction stored at the previous time level.
    double titEini = -1;
    /// Inlet void fraction stored at the previous time level.
    double alfEini = -1;
    /// Inlet complementary-liquid fraction stored at the previous time level.
    double betaEini = -1;
    /// Previous-time-level value of presfim.
    double presfimini = 0;
    /// Gas mass fraction imposed during reverse flow at the last control volume. Used only by transient network
    /// simulations.
    double titRev = 1.;
    /// Previous-time-level value of titRev.
    double titRevini = 1.;
    /// Complementary-liquid fraction imposed during reverse flow at the last control volume. Used only by transient
    /// network simulations.
    double betaRev = 0.;
    /// Previous-time-level value of betaRev.
    double betaRevini = 1.;
    /// Previous-time-level downstream or separator temperature.
    double tGSupIni = 0.;
    /// Temperature downstream of the surface choke, delivered to the next network section.
    double tempSup = 0;
    /// Previous-time-level value of masChkSup.
    int masChkSupini = 0;
    /// Previous-time-level value of mudaModoChk.
    int mudaModoChkini = 0;
    /// Secondary-line index where parallel-network coupling with the primary line begins.
    int SecPrimIniRedeP = 0;
    /// Secondary-line index where parallel-network coupling with the primary line ends.
    int SecPrimFimRedeP = 0;
    /// Primary-line index aligned with SecPrimIniRedeP.
    int PrimSecIniRedeP = 0;
    /// Primary-line index aligned with SecPrimFimRedeP.
    int PrimSecFimRedeP = 0;
    /// Indicates thermal coupling on the primary branch of a parallel network.
    int verificaAcopRedeP = 0;
    /// Indicates thermal coupling on the secondary branch of a parallel network.
    int verificaAcopRedeS = 0;
    /// Initial source indices for parallel-network coupling.
    std::vector<int> indFonteRedeParalelaIni;
    /// Initial production-liquid sources for the parallel network.
    std::vector<double> fonteMpRedeParalelaIni;
    /// Initial complementary-liquid sources for the parallel network.
    std::vector<double> fonteMcRedeParalelaIni;
    /// Initial gas sources for the parallel network.
    std::vector<double> fonteMgRedeParalelaIni;
    /// Boundary-condition type applied to the secondary parallel-network branch.
    int redeParalelaCCsecundario = -1;
    /// Primary branch index in a parallel network.
    int redeParalelaP = -1;
    /// Secondary branch index in a parallel network.
    int redeParalelaS = -1;
};

/// What the steady-state search carries from one iteration to the next: the convergence
/// monitor and its reference, the initial estimates, the reverse flow of a network solution
/// and its fluid and temperature, the number of production fluids and the slow thermal
/// coupling parameter.
struct SteadySearch {
    /// Control parameter for slowly varying thermal coupling.
    double trocaTermicaLenta = kSlowHeatTransferThreshold;
    /// Number of production fluids.
    int nfluP = 0;
    /// Initial holdup estimate used by the steady-state solver.
    double chuteHol = -1.;
    /// Controls the search for an initial steady-state estimate.
    int buscaIni = 0;
    /// Fluid state received from downstream during reverse network flow.
    ProFlu fluiRevRede;
    /// Temperature associated with reverse network flow.
    double tempRev = 0.;
    /// Indicates reverse flow in the steady-state network solution.
    int revPerm = 0;
    /// Current steady-state convergence monitor.
    double monitConvPerm = kInitialConvergenceMonitor;
    /// Reference value for the steady-state convergence monitor.
    double monitConvPermBase = 1.;
};

/// The coupling with the three-dimensional thermal model: its Poisson solver and the index
/// corrections of the coupled sections.
struct Coupling3D {
    /// Index corrections used by coupled thermal sections.
    std::vector<int> acertaIndAcop;
    /// Three-dimensional Poisson solver used by the thermal model.
    solverP3D poisson3D;
};

}  // namespace sisprod

/// Models and solves a one-dimensional production system.
///
/// The class owns the production and gas-service-line state, fluid-property tables, boundary conditions,
/// network-coupling data, transient buffers, and steady-state/transient solution procedures.
class SProd {
    friend sisprod::gaslift::GasLiftState sisprod::adapters::gasLiftStateOf(SProd &, sisprod::SolveContext &);
    friend sisprod::steady::SteadyStateState sisprod::adapters::steadyStateOf(SProd &, sisprod::SolveContext &);
    friend sisprod::steady::SteadyStateSearchState sisprod::adapters::searchStateOf(SProd &, sisprod::SolveContext &);
    friend sisprod::transient::TransientStepState sisprod::adapters::transientStateOf(SProd &, sisprod::SolveContext &);
    friend sisprod::composition::CompositionState
    sisprod::adapters::compositionStateOf(SProd &, sisprod::SolveContext &);
    friend sisprod::transient::TransientSolveState
    sisprod::adapters::transientSolveStateOf(SProd &, sisprod::SolveContext &);
    friend sisprod::thermal::ThermalState sisprod::adapters::thermalStateOf(SProd &, sisprod::SolveContext &);
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

    /// Constructs and initializes a production-system simulation from the input and log files.
    SProd(std::string inputFile, std::string logFile, tipoValidacaoJson_t jsonValidation,
          tipoSimulacao_t simulationType, varGlob1D *globals = nullptr, int dynamicTable = -1, int blocking = 0,
          int reversedGeometry = 0,
          double *extraSourceLengths = nullptr,
          int *extraSourceCells = nullptr,
          int extraSourceCount = 0,
          int steadyBranch = 1);

    /// Creates an empty production-system object.
    SProd();

    /// Builds a system from an already parsed input, without reading the JSON file
    /// again, with the place in the network and the reverse-flow state in carried.
    SProd(Ler &parsedInput, const CarriedState &carried);

    /// Not copyable: the object owns raw arrays, which operator= rebuilds
    /// instead of sharing.
    SProd(const SProd &) = delete;

    /// Rebuilds this system from other's input file, read again from disk, and takes
    /// over other's place in the network and reverse-flow state. Not a copy: the run
    /// state starts afresh and the system is assembled again.
    SProd &operator=(const SProd &other);

    /// The state another system built from this one's input takes over.
    CarriedState carriedState() const;

    /// Loads pressure-velocity results into cell and face state variables.
    void renova(int explicitScheme = 0);
    /// Updates only section-end states during intermediate network convergence.
    void renovaBuffer();
    /// Copies previous states when an active choke bypasses the intermediate network solve.
    void renovaBufferCego();

    /// Updates distributed mass-transfer terms used by void-fraction and mixture-mass equations.
    void renovaTemp();
    /// Computes T1 and T2 used to split mixture mass flow into liquid and gas flows.
    void renovaterm(int inflowMode = 0);
    /// Computes T1 and T2 at the outlet of an internal network section.
    void renovatermAfluFim();
    /// Computes T1 and T2 at the inlet of an internal network section.
    void renovatermColIni();

    /// Calculates outlet-choke flow from the last-cell and separator pressures.
    void calcCCpres(double reverseQuality = 1., double reverseVoidFraction = 1.,
                    double reverseLiquidFraction = 0.);
    /// Calculates outlet-choke flow for the intermediate network state.
    void calcCCBuffer(double reverseQuality = 1., double reverseVoidFraction = 1.,
                      double reverseLiquidFraction = 0.);

    /// Selects a stable time step from CFL and additional model restrictions.
    void determinaDT(int explicitScheme = 0);
    /// Limits time-step growth when accepted steps remain well below the CFL estimate.
    void atenuaDtMax();
    /// Checks whether liquid-density variation requires the complete formulation.
    void avaliaVariaDpDt(int explicitScheme = 0);
    /// Updates the valve logic for the active-state formulation.
    void aberturaVal1();
    /// Updates the valve logic for the inactive-state formulation.
    void aberturaVal0();
    /// Restricts the time step during valve-state transitions.
    void restringeDTporValv();

    /// Advances the gas-injection system and prepares its coupling with the production line.
    void solveLinGas();
    /// Advances phase volume fractions.
    void EvoluiFrac(double reverseVoidFraction = 1., double reverseLiquidFraction = 0., int cycle = 0);
    /// Restores the initial fraction state after an invalid update.
    void ReiniEvolFrac0();
    /// Restores previous volume fractions after a nonphysical update.
    void ReiniEvolFrac();
    /// Updates phase fractions in cells affected by a moving pig.
    void AtualizaPig();
    /// Assembles and solves the global pressure-velocity coupling system.
    void SolveAcopPV(int explicitScheme = 0, int cycle = 0);

    /// Advances the transient energy equation.
    void marchaEnergTrans(int cycle = 0, int maximumCycle = 0);
    /// Refreshes local dynamic fluid-property tables.
    void atualizaMiniTab();
    /// Updates the first boundary condition.
    void atualizaCC1();

    /// Advances the complete transient production-system solution and handles output and logging. collectorFluid is
    /// the fluid that flows back from the collector, null when the system has none.
    void SolveTrans(double reverseQuality = 1., double reverseVoidFraction = 1., double reverseLiquidFraction = 0.,
                    int network = -1, const ProFlu *collectorFluid = nullptr);

    /// Appends the buffered rows of a production-line trend file.
    void ImprimeTrendP(int trend, int network = -1);
    /// Writes the column captions of a production-line trend file, truncating it.
    void ImprimeTrendPCab(int trend, int network = -1);
    /// Appends the buffered rows of a service-line trend file.
    void ImprimeTrendG(int trend, int network = -1);
    /// Writes the column captions of a service-line trend file, truncating it.
    void ImprimeTrendGCab(int trend, int network = -1);

    /// Marches the steady production solution using a bottomhole-pressure guess with outlet pressure prescribed.
    double marchaProdPerm1(double pressureGuess);
    /// Brackets and solves the bottomhole-pressure root for marchaProdPerm1.
    double buscaProdPfundoPerm(double guess = -1., int attemptCount = -1);
    /// Brackets and solves the reverse-flow bottomhole-pressure root.
    double buscaProdPfundoPermRev(double guess = -1.);
    /// Brackets and solves the bottomhole-pressure root for marchaProdPerm2.
    double buscaProdPfundoPerm2(double guess = -1., int attemptCount = -1);
    /// Runs a direct march when inlet pressure and flow rate are known.
    double buscaProdPfundoPerm3(double inletPressure);
    /// Marches the steady production solution using a bottomhole mass-flow guess.
    double marchaProdPresPres1(double massFlowRateGuess);
    /// Brackets and solves the mass-flow root for marchaProdPresPres1.
    double buscaProdPresPresPerm(double massFlowRateGuess, double maximumFlowRate = 0., int iterationCount = 0);
    /// Brackets and solves the reverse-flow mass-flow root.
    double buscaProdPresPresPermRev(double massFlowRateGuess, double maximumFlowRate = 0., int iterationCount = 0);
    /// Brackets and solves the mass-flow root for marchaProdPresPres2.
    double buscaProdPresPresPerm2(double massFlowRateGuess, double maximumFlowRate = 0.);
    /// Brackets and solves the mass-flow root for the closed-choke case, chosen
    /// by Num4Main when arq.chokep.abertura[0] <= 1e-15.
    double buscaProdPresPresPerm3(double massFlowRateGuess, double maximumFlowRate = 0.);

    /// Solves injection cases CC1 and CC3.
    double buscaInjPfundoPerm1(double guess = -1.);
    /// Solves injection case CC0.
    double buscaInjPfundoPerm2(double guess = -1.);
    /// Solves injection case CC2.
    double buscaInjPfundoPerm3(double guess = -1.);
    /// Solves injection case CC4.
    double buscaInjPfundoPerm4();
    /// Solves injection case CC5.
    double buscaInjPfundoPerm5(double guess = -1.);

    /// Estimates steady production-network node pressures from hydrostatics.
    double hidroreverso(double liquidHoldup, double liquidFlowRate = 0, double gasFlowRate = 0);
    /// Estimates steady injection-network node pressures from hydrostatics.
    double hidroreversoInj(double liquidHoldup, double liquidFlowRate = 0);
    /// Estimates secondary-branch pressures in a parallel production network.
    double hidroTramoSecundario(double quality);

    /// Calculates temperature downstream from the surface choke for network coupling.
    void calcTempFim();

    /// Enables latent-heat calculations for black-oil simulations.
    int CalcLat = 0;

    /// Section index when this object belongs to a pipeline network.
    int indTramo = -1;
    /// Number of control volumes in the production line.
    int ncel = 0;

    /// What the steady-state search carries from one iteration to the next.
    sisprod::SteadySearch steadySearch;

    /// Current inlet pressure boundary condition.
    double presE = -1;
    /// Current inlet temperature boundary condition.
    double tempE = -1;
    /// Current inlet gas mass fraction.
    double titE = -1;
    /// Computed inlet void fraction.
    double alfE = -1;
    /// Computed inlet complementary-liquid fraction.
    double betaE = -1;
    /// What this system keeps as a branch of a network.
    sisprod::NetworkCoupling networkCoupling;

    /// Pressure in the last production-line control volume. When the surface choke is open, this value is equal to
    /// the separator pressure.
    double presfim = 0;
    /// Separator temperature or temperature at the inlet of the downstream section.
    double tGSup = 0.;
    /// Production-line inlet temperature when no inlet-pressure boundary condition is imposed.
    double temperatura = 0;

    /// Indicates whether the surface choke is active.
    int masChkSup = 0;
    /// The fluid-property tables the fluids of this system's cells point at.
    sisprod::PropertyTables tables;

    /// Parsed user input and simulation configuration.
    Ler arq;

  private:
    /// Drift-flux correlation chosen for each flow regime.
    ///
    /// Resolved from arq once, when the object is built, so the per cell path never reads configuration.
    driftflux::correlations::RegimeSelectors driftSelectors;

  public:
    /// Buffer used to write gas-line profiles, one row until the system is sized for its input.
    FullMtx<double> flutG{1, sisprod::kGasLineProfileFixedColumns};
    /// Buffer used to write production-line profiles, one row until the system is sized for its input.
    FullMtx<double> flut{1, sisprod::kProductionProfileFixedColumns};
    /// The gas-lift line: its cells, chokes and valves, its pressure-velocity system, and the unloading of its
    /// completion fluid.
    sisprod::GasLiftLine gasLift;
    /// The state of a transient run.
    sisprod::TransientRun transient;

    /// Current time step.
    double dt = 0.;

    /// The trend and profile output this system writes.
    sisprod::TrendRecorder trends;

    /// Event-log file name.
    std::string tmpLog;

  private:
    /// Number of iterations used to bracket the initial steady-state root.
    int steadyIteration = 0;

  public:
    /// Steady-state mode flag: 1 while Num4Main solves a network branch at steady state. The source terms read it
    /// (sisprod::sources::SourceState).
    int modoPerm = 0;

    /// Surface-choke model.
    choke chokeSup;
    /// Multiphase production-line control volumes, held by productionCells.
    Cel *celula = nullptr;

  private:
    /// The production-line control volumes celula points at.
    std::vector<Cel> productionCells;

  public:
    /// Indicates that the section outlet is not connected to another network section.
    int noextremo = 1;
    /// Indicates that the section inlet is not connected to another network section.
    int noinicial = 1;
    /// Indicates that this section is a branch of a gas-lift ring network.
    int derivaAnel = -1;

    /// Section-blocking state.
    int bloq = 0;

    /// The coupling with the three-dimensional thermal model.
    sisprod::Coupling3D coupling3D;
    /// Shared one-dimensional simulation settings.
    varGlob1D *vg1dSP = nullptr;

  private:
    /// Indicates that the thermal source term is disabled.
    int thermalSourceDisabled = 0;

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

    /// Puts the run state back where a system built from an input starts: the property tables and the steady-state
    /// search as declared, the trend pointers null, the boundary conditions, the time-step statistics and the counters
    /// at their initial values, the previous-time-level interface state equal to the current one, and the histories
    /// empty.
    void resetRunState();
    /// Sizes the profile buffers and the pressure-velocity systems of both lines for arq's cells, all zero.
    void allocateSolutionBuffers();

    /// Gives back the memory of every array this object owns. operator= calls it before rebuilding, so that the
    /// arrays of the previous system are free before the input is read again.
    void releaseOwnedStorage();

    /// Takes over the place in the network and the reverse-flow state in carried.
    void applyCarriedState(const CarriedState &carried);
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
    /// at extraSourceLengths, whose cell indices it records in extraSourceCells.
    void buildProductionCells(double *extraSourceLengths, int *extraSourceCells, int extraSourceCount);
    /// Gives the inlet the sources its boundary condition needs (and the second
    /// cell, under a blockage), then places the accessories -- pumps, volumetric
    /// pumps, pressure-drop requirements, heat sources, the master valve and the
    /// other valves -- and sets up the outlet pressure, the surface and injection
    /// chokes and the pigs.
    void configureInletSourcesAndAccessories(int extraSourceCount);
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
    /// Lists the cells with two-dimensional heat diffusion, copies the master-valve opening and closing times, maps
    /// the transient profile positions to global thermal-node indices, and allocates the trend matrices of the
    /// production line, the gas line and the transient trends, sized to the longest simulated time and filled with
    /// kTrendSentinel.
    void allocateEventProfileAndTrendArrays();
    /// Resets the column-annulus and network coupling flags, locates the
    /// column-annulus coupling range on the gas line, sets the steady and transient
    /// profile counters and their first output times, opens the event log with the
    /// events known at start, finds the smallest cell length, and zeroes the
    /// moving-average and running-total state the transient loop starts from.
    void resetCouplingAndOutputState();
    /// Builds the production section after input parsing.
    void montasistema(double *extraSourceLengths = nullptr,
                      int *extraSourceCells = nullptr,
                      int extraSourceCount = 0);

    /// Initializes the gas service line for gas-lift unloading.
    void HidroDescargaG();
    /// Initializes the production line for gas-lift unloading.
    void HidroDescargaP();

    /// Runs the hydrate-envelope solvers for the production and gas lines.
    void solveHydrateEnvelopes();

    /// Caches the drift-flux correlation of each regime from arq.
    ///
    /// Called wherever arq is built or replaced, so driftSelectors never goes stale.
    void resolveDriftSelectors();
};

#endif /* SISPROD_H_ */
