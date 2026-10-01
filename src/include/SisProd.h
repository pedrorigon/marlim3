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

/// The fluid-property tables the fluids of a system's cells point at. The black-oil
/// compressibility-factor, specific-heat and liquid-density-derivative tables belong to the
/// input; the latent-heat table and the Livia bubble-point and solution-gas-ratio tables are
/// built and owned by the system.
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
    double **HLat = nullptr;
    /**
     * @brief Bubble-pressure values imported from PVTSim for the Livia solution-gas-ratio correlation. This table
     * is separate from the full PVTSim fluid-property model.
     */
    double *PBPVTSim = nullptr;
    /**
     * @brief Bubble-temperature values imported from PVTSim for the Livia solution-gas-ratio correlation. This
     * table is separate from the full PVTSim fluid-property model.
     */
    double *TBPVTSim = nullptr;
    /**
     * @brief Precomputed solution-gas-ratio table for the Livia correlation, used to avoid repeating its
     * expensive calculation during the simulation.
     */
    double **RSLivia = nullptr;
};

/// The gas-lift line's pressure-velocity system, its valves and the cells they join, the
/// column-annulus thermal coupling, and the unloading controller's samples and settings.
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
    int *posicVGLP = nullptr;
    /**
     * @brief Service-line cell indices associated with gas-lift valves.
     */
    int *posicVGLG = nullptr;
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
    ChokeGas *chokeVGL = nullptr;
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
};

/// The state of a transient run: the time-step and change-rate histories and the
/// restrictions on the step, the Master1 valve schedule, the production-line
/// pressure-velocity matrix, the event and log counters, the transport switches, the pigs
/// and the cells of the two-dimensional Poisson model.
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
    int *receb = nullptr;
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
    double *fechaM1 = nullptr;
    /**
     * @brief Times at which Master1 opens.
     */
    double *abreM1 = nullptr;
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
    /**
     * @brief Requests rollback and time-step reevaluation when at least one control volume produces a holdup or
     * volume fraction outside the physical [0, 1] range.
     */
    int reinicia = 0;

  private:
    /**
     * @brief Control parameter for slowly varying thermal coupling.
     */
    double trocaTermicaLenta = 0.01;
  public:
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
     * @brief Pressure in the last production-line control volume. When the surface choke is open, this value is
     * equal to the separator pressure.
     */
    double presfim = 0;
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
     * @brief Separator pressure or pressure at the inlet of the downstream network section.
     */
    double pGSup = 0;
    /**
     * @brief Previous-time-level value of pGSup.
     */
    double pGSupIni = 0.;
    /**
     * @brief Previous-time-level downstream or separator temperature.
     */
    double tGSupIni = 0.;
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
    /**
     * @brief Reserved temperature state; currently unused.
     */
    double tempSup = 0;
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
  private:
  public:
    /**
     * @brief Indicates whether the surface choke is active.
     */
    int masChkSup = 0;
    /**
     * @brief Previous-time-level value of masChkSup.
     */
    int masChkSupini = 0;
    /**
     * @brief Signals a surface-choke operating-mode transition.
     */
    int mudaModoChk = 0;
    /**
     * @brief Previous-time-level value of mudaModoChk.
     */
    int mudaModoChkini = 0;
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
  private:
    /**
     * @brief Number of production fluids.
     */
    int nfluP = 0;
    /**
     * @brief The fluid-property tables the fluids of this system's cells point at.
     */
    sisprod::PropertyTables tables;
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
    /**
     * @brief The gas-lift line's pressure-velocity system, valves, column-annulus coupling
     * and unloading controller.
     */
    sisprod::GasLiftLine gasLift;
    /**
     * @brief The state of a transient run.
     */
    sisprod::TransientRun transient;
  public:
    /**
     * @brief Right-hand side and solution vector for the production-line pressure-velocity system.
     */
    Vcr<double> termolivreP;

    /**
     * @brief Current time step.
     */
    double dt = 0.;
  private:
  public:
    /**
     * @brief Simulation end time.
     */
    double tfinal = 0;

  private:


    /**
     * @brief Gas-line cells where radial temperature profiles are written.
     */
    int *ncelperftransg = nullptr;
    /**
     * @brief Maximum number of samples stored for each gas-line trend.
     */
    int *TrendLengthG = nullptr;
  public:
    /**
     * @brief Buffered gas-line trend data.
     */
    double ***MatTrendG = nullptr;
  private:
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
     * @brief Maximum number of wall-temperature samples stored for each gas-line trend.
     */
    int *TrendLengthTransG = nullptr;
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
    int *ncelperftransp = nullptr;
    /**
     * @brief Maximum number of samples stored for each production-line trend.
     */
    int *TrendLengthP = nullptr;
  public:
    /**
     * @brief Buffered production-line trend data.
     */
    double ***MatTrendP = nullptr;
  private:
    /**
     * @brief Times at which production-line trend buffers are reset.
     */
    double *resettrend = nullptr;
  public:
    /**
     * @brief Number of production-line trend samples currently stored.
     */
    int *ntrend = nullptr;
  private:
    /**
     * @brief Number of production-line trend samples stored before the last flush.
     */
    int *ntrendB = nullptr;
    /**
     * @brief Maximum number of wall-temperature samples stored for each production-line trend.
     */
    int *TrendLengthTransP = nullptr;
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

  public:
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

  private:
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
    /**
     * @brief Current complete-model activation state.
     */
    int modeloCompleto = 1;
  private:
  public:
    /**
     * @brief Auxiliary CFL time-step accumulator.
     */
    double dtauxCFL = 0.;
    /**
     * @brief Auxiliary accepted-time-step accumulator.
     */
    double dtauxFinal = 0.;
  private:
  public:

    /**
     * @brief Surface-choke model.
     */
    choke chokeSup;
    /**
     * @brief Gas-injection choke model.
     */
    ChokeGas chokeInj;
  private:
  public:
    /**
     * @brief Gas service-line control volumes.
     */
    CelG *celulaG = nullptr;
    /**
     * @brief Multiphase production-line control volumes.
     */
    Cel *celula = nullptr;

  private:
    /**
     * @brief Enables reading bubble-pressure and bubble-temperature tables in black-oil mode.
     */
    int LerPB = 0;
    /**
     * @brief Enables reading a solution-gas-ratio table in black-oil mode.
     */
    int lerRS = 0;
  public:

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

  private:
    /**
     * @brief Temporary-network flag; currently expected to remain zero.
     */
    int redeTemporario = 0;
  public:

    /**
     * @brief Trend-output cycle counter; a value of one triggers header output.
     */
    double kimpT = 0.;

  private:


  public:
    /**
     * @brief Initial holdup estimate used by the steady-state solver.
     */
    double chuteHol = -1.;
    /**
     * @brief Controls the search for an initial steady-state estimate.
     */
    int buscaIni = 0;

    /**
     * @brief Dynamic fluid-property tables.
     */
    vector<tabelaDinamica> tabDin;
    /**
     * @brief Number of dynamic property tables.
     */
    int ntabDin = 0;
  private:
  public:
    /**
     * @brief Section-blocking state.
     */
    int bloq = 0;

    /**
     * @brief Fluid state received from downstream during reverse network flow.
     */
    ProFlu fluiRevRede;
    /**
     * @brief Temperature associated with reverse network flow.
     */
    double tempRev = 0.;
  private:
    /**
     * @brief Indicates reverse flow in the steady-state network solution.
     */
    int revPerm = 0;
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
    /**
     * @brief Minimum time step allowed by the current cycle. Zero until a
     * transient cycle sets it.
     */
    double dtCicMin = 0.;

  private:
    /**
     * @brief Indicates that the thermal source term is disabled.
     */
    int semTermo = 0;
    /**
     * @brief Current steady-state convergence monitor.
     */
    double monitConvPerm = 1000.;
    /**
     * @brief Reference value for the steady-state convergence monitor.
     */
    double monitConvPermBase = 1.;
  public:

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

  private:
    vector<int> kontaTempoCelUni;

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

    /// Releases dynamically allocated simulation buffers and cell arrays.
    ~SProd();

  private:
    /// Frees every array this object owns, reading its current sizes and switches.
    /// The destructor calls it, and operator= and the constructor from a parsed
    /// input call it before copying anything in.
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

    /// Estimates the gas-injection pressure correction required to avoid erosional valve velocity.
    double prescordesc(double velmax, int ivalv, double fator, int sinal);
    /// Computes the unloading injection-pressure correction for one gas-lift valve.
    double CalcPresValvDesc(double velGarg, int ivalv);
    /// Controls injection and upstream-choke pressures from gas-lift-valve flow rates.
    double BuscaPresInjDesc();
    /// Updates gas-line state after solving pressure-velocity coupling.
    void renovaGas();
    /// Updates intermediate gas-line state during network convergence.
    void renovaGasBuf();

    /// Calculates gas-lift-valve opening area from calibration and operating conditions.
    double areaValvCali(double PCal, double TCal, double PVO, double PT,
                        double dextern, double areagarg, double Rvalv, double Temp);
    /// Advances the temperature of one gas-line control volume.
    void calctempGas(int i, double tempantiga, int modoPerm = 0);
    /// Solves gas-line pressure and flow in the completion-fluid region during unloading.
    void resolveDescarga();
    /// Updates gas-line temperature in the completion-fluid region during unloading.
    void tempDescarga(int i);
    /// Advances the completion-fluid/gas interface in the service line.
    void avancInter();
    /// Calculates gas temperature across a gas-lift valve using the Joule-Thomson model.
    double TempDescGL(int igl);
    /// Maps gas-lift-valve positions to gas-line control volumes.
    void ValvGasTrans();
    /// Advances the coupled gas-line pressure, velocity, and temperature solution.
    void subtempoGas();
    /// Advances the intermediate gas-line state used by network convergence.
    void subtempoGasBuf();
    /// Exchanges heat-transfer data between the production column and annulus.
    void conectaColuna();
    /// Interpolates latent heat from enthalpy tables.
    double interpolaHLatente(double pres, double temp);
    /// Advances the temperature of one production-line control volume.
    void calctemp(int i, double tempantiga, int modoPerm = 0);
    /// Returns the mixture enthalpy helper value; currently unused.
    double calcHmix(int i);
    /// Returns the mixture-energy helper value; currently unused.
    double energmix(int i, int jp0, int jt, double razp);
    /// Updates temperature from enthalpy; currently unused.
    void calcTempEntalp(int i);
    /// Evaluates thermal mass-transfer terms; currently unused.
    void calcTransMassTermo(int i);

    /// Calculates flow through Master1 while it operates as a choke.
    void FonteValv(int ind);
    /// Stores source terms from the previous time level for possible rollback.
    /// Nothing in the product calls it.
    void salvaFonte();
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
    /// Applies hydrostatic and friction corrections; currently unused.
    void correcHidroFric(int i, double &hidro, double &fric);
    /// Prepares auxiliary data for a local fluid-property table. Nothing in the
    /// product calls it.
    void auxMiniTab(ProFlu &flu);
    /// Generates the local fluid-property table. Nothing in the product calls it.
    void geraMiniTabFlu();
  public:

    /// Loads pressure-velocity results into cell and face state variables.
    void renova(int expli = 0);
  private:
    /// Updates phase and mixture flow rates.
    void renovaVaz();
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
    /// Selects a stable time step from CFL and additional model restrictions.
    void determinaDTExpli();
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
    /// Stores only the volume fractions required by fraction rollback.
    void SubReiniEvolFrac();
  public:
    /// Restores previous volume fractions after a nonphysical update.
    void ReiniEvolFrac();
    /// Updates phase fractions in cells affected by a moving pig.
    void AtualizaPig();
    /// Assembles and solves the global pressure-velocity coupling system.
    void SolveAcopPV(int vexpli = 0, int ciclo = 0);

  private:
    /// Prepares the multidimensional heat-diffusion problem for one cell.
    void prepDifusCalorND(int i);
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
    /// Marches the steady production solution using a bottomhole-pressure guess with outlet pressure prescribed.
    double marchaProdPerm1Rev(double pchute);
    /// Marches the steady production solution using a bottomhole-pressure guess with an outlet choke.
    double marchaProdPerm2(double pchute);
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
    /// Marches the steady production solution using a bottomhole mass-flow guess.
    double marchaProdPresPres1Rev(double mchute);
  public:
    /// Brackets and solves the mass-flow root for marchaProdPresPres1.
    double buscaProdPresPresPerm(double mchute, double maxvaz = 0., int kontaiter = 0);
    /// Brackets and solves the reverse-flow mass-flow root.
    double buscaProdPresPresPermRev(double mchute, double maxvaz = 0., int kontaiter = 0);
  private:
    /// Marches the steady production solution for the second pressure-pressure boundary formulation.
    double marchaProdPresPres2(double mchute);
  public:
    /// Brackets and solves the mass-flow root for marchaProdPresPres2.
    double buscaProdPresPresPerm2(double mchute, double maxvaz = 0.);
  private:
    /// Marches the steady production solution for the third pressure-pressure
    /// boundary formulation: marchaProdPresPres2 with the choke's throughput
    /// forced to zero, which is the closed-choke case.
    ///
    /// Nothing calls it: marchaProdPresPres2 already reduces to exactly this when
    /// the choke is shut, because the throat area is `abertura[0] * area` and both
    /// vazmaxSachd and vazmassSachd scale to zero with it -- leaving the same
    /// `0. - celula[ncel - 1].MR` residual this function returns literally.
    double marchaProdPresPres3(double mchute);
  public:
    /// Brackets and solves the mass-flow root for the closed-choke case, chosen
    /// by Num4Main when arq.chokep.abertura[0] <= 1e-15.
    ///
    /// It reaches marchaProdPresPres2, not marchaProdPresPres3 -- via zriddr(.., 2,
    /// 1) and multMarcha -- which is correct by the reduction described above.
    double buscaProdPresPresPerm3(double mchute, double maxvaz = 0.);

  private:
    /// Marches the steady gas line with prescribed injection pressure and a mass-flow guess.
    double marchaGasPerm1(double chutemass = -1);
    /// Marches the steady gas line with prescribed injection flow and a pressure guess.
    double marchaGasPerm2(double pchute, double chutemass = -1);
    /// Marches the steady gas line across the injection choke using a downstream-pressure guess.
    double marchaGasPerm3(double pchute);
    /// Brackets and solves the pressure root for marchaGasPerm2.
    double buscaGasPresPerm2();
    /// Brackets and solves the pressure root for marchaGasPerm3.
    double buscaGasPresPerm3();
    /// Marches pressure from the previous cell center to the downstream face.
    void RenovaPresPermMon(int i, int RK);
    /// Marches pressure from the last cell center to the outlet face.
    double RenovaPresPermNcel();
    /// Calculates the pressure contribution caused by an area change.
    double calcDpArea(int i, double rhomix, double rey, double jmix);
    /// Marches pressure from the upstream face to the current cell center.
    void RenovaPresPermJus(int i, int RK);
    /// Corrects gas density in one control volume.
    void corrDeng(int i);
    /// Updates steady-state mass flow and fluid state after a source term.
    void RenovaMassPerm(int i);
    /// Updates steady-state mass flow for reverse flow.
    void RenovaMassPermRev(int i);
    /// Updates steady-state mass flow and pseudocomponent composition after a source term.
    void RenovaMassPermComp(int i);
    /// Updates steady-state mass flow and pseudocomponent composition after a source term.
    void RenovaMassPermCompRev(int i);
    /// Calculates steady-state slip parameters at a downstream face.
    void CalcC0UdPerm(int ind, double &c0, double &ud);
    /// Calculates steady-state interphase mass transfer.
    void RenovaTransMassPerm(int i);
    /// Calculates steady-state interphase mass transfer.
    void RenovaTransMassPermGas(int i);
    /// Marches steady-state temperature from cell i-1 to cell i.
    void RenovaTempPerm(int i, int RK);
    /// Marches steady-state temperature in the reverse direction.
    void RenovaTempPermRev(int i, int RK);
    /// Adds pump pressure gain at the upstream face of a production cell.
    void atualizaPeriPmonProd(int i);
    /// Synchronizes neighboring face pressures after updating a cell-center pressure.
    void atualizaPeriPjusProd(int i);
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
    /// Refreshes fluid properties.
    void atualizaProp();
    /// Updates steady-state velocities and thermal terms.
    void atualizaVelTermPerm();
    /// Calculates the pseudo-transient time step.
    void calcDTPseudoTrans();

    /// Marches the steady injection-well solution using a pressure or flow-rate guess.
    double marchaInjPerm1(double chute);
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
    /// Dispatches the selected production, gas-line, or injection steady-state marching method.
    double multMarcha(double chute, int prod, int tipoCC);

    /// Solves for the steady-state boundary condition by Ridders' method.
    ///
    /// Binds the production domain to rootfinding::zriddr: it captures prod and
    /// tipoCC in the objective, derives the minimum iteration count from the
    /// input deck, and carries the convergence monitor. The algorithm itself
    /// knows nothing of that.
    double zriddr(double x1, double x2, int prod, int tipoCC);
  public:

    /// Estimates steady production-network node pressures from hydrostatics.
    double hidroreverso(double hol, double vaz = 0, double vazG = 0);
    /// Estimates steady injection-network node pressures from hydrostatics.
    double hidroreversoInj(double hol, double vaz = 0);
    /// Estimates secondary-branch pressures in a parallel production network.
    double hidroTramoSecundario(double titulo);
  private:
    /// Builds a hydrostatic estimate for the gas service line.
    void hidroLinServ();

    /// Solves the secondary-branch flow rate in a steady parallel network.
    double buscaTramoSecVazPerm(double pPartida, int indPartida);
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
