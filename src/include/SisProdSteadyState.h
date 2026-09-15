#ifndef SISPRODSTEADYSTATE_H_
#define SISPRODSTEADYSTATE_H_

// Declared, not included: SteadyStateState holds only pointers and references
// to these, so this header stays free of the cell, input-deck and choke headers
// and can still be compiled on its own.
class Cel;
class CelG;
class ChokeGas;
class Ler;
class ProFlu;
class choke;
struct varGlob1D;
class SProd;

namespace sisprod::steady {

/// The callbacks the steady march needs from outside itself.
///
/// Fourteen of these already live in extracted modules -- nine in
/// sisprod::gaslift, five in sisprod::thermal -- and could in principle be
/// called directly. They cannot be: reaching them needs a GasLiftState or a
/// ThermalState, and the adapters that build those (gasLiftStateOf,
/// thermalStateOf) are internal to SisProd.cpp, where the SProd they read is.
/// So the call goes back through SProd exactly as GasLiftTemperatureUpdater
/// does, and SisProd.cpp remains the one place that knows how to assemble any
/// module's state.
///
/// The remaining two are still SProd's own: CalcC0UdPerm and renovaFonte.
///
/// Four of these signatures were WRONG when this header was first written, and
/// the compiler said so on the first attempt to define them. CalcC0UdPerm
/// returns its two coefficients through reference parameters, the two
/// temperature marches take a Runge-Kutta stage, and calctemp takes the
/// previous temperature and a steady-mode flag. They are transcribed from
/// SisProd.h now rather than assumed from the call sites.
struct SteadyStateUpdaters {
    SProd &system;

    // --- drift closure and sources -----------------------------------------
    void steadyDriftClosure(int cellIndex, double &c0, double &ud) const;
    void updateSource(int cellIndex) const;

    // --- thermal ------------------------------------------------------------
    void advanceSteadyTemperature(int cellIndex, int rungeKuttaStage) const;
    void advanceReverseSteadyTemperature(int cellIndex, int rungeKuttaStage) const;
    void computeTemperature(int cellIndex, double previousTemperature, int steadyMode) const;
    void computeGasTemperature(int cellIndex, double previousTemperature, int steadyMode) const;
    void updateProductionTemperaturePeriphery(int cellIndex) const;

    // --- gas line -----------------------------------------------------------
    void initializeSteadyValveGasFlowRate(int cellIndex) const;
    void initializeTubingConnectionSteady() const;
    void updateSteadyGasPressure(int cellIndex) const;
    void updateSteadyGasTemperature(int cellIndex) const;
    void computeSteadyGasFlowRate(int cellIndex) const;
    void connectTubing() const;
    void connectTubingSteady() const;
    [[nodiscard]] double steadyGasPressureDrop(int cellIndex) const;
    [[nodiscard]] double steadyInjectionPressureDrop(int cellIndex) const;
};

/// The state the steady march reads, and the only state it may read.
///
/// Same role as ThermalState and GasLiftState: it names in one place what this
/// domain touches, and it makes the routines callable WITHOUT an SProd, which
/// is what lets a dedicated harness drive them over synthetic cells.
///
/// The field list is not a guess. It was derived by walking all 31 march
/// bodies -- 6,492 lines -- and intersecting every identifier against SProd's
/// declared data members. Twenty-eight came back; seventeen of them are also
/// read by the searches and are the reason SteadyStateSearchState composes this
/// struct rather than restating it.
///
/// Scalars are held BY REFERENCE, not by value. Copying them into the struct
/// would read every one at construction, before the branch that decides whether
/// the original would have read it at all -- the defect the root-finding stage
/// had to undo, and the reason DriftFluxClosure's ClosureState and
/// GasLiftState carry the same warning.
///
/// const marks what the march does not write. That was measured the same way,
/// not assumed: stage 6 shipped a header promising `const Ler &input` for a
/// deck the module writes through, and only the compiler caught it.
struct SteadyStateState {
    /// Production cells -- SProd::celula. Written.
    Cel *cells;
    /// Gas-line cells -- SProd::celulaG. Written.
    CelG *gasCells;
    /// Input deck -- SProd::arq. NOT const: the march writes back into it.
    Ler &input;
    /// Shared 1D globals -- SProd::vg1dSP. Read only.
    const varGlob1D *globals;

    /// Index of the last production cell -- SProd::ncel.
    const int &lastCell;
    /// Index of the last gas-line cell -- SProd::ncelGas.
    const int &gasCellCount;

    /// Injection choke -- SProd::chokeInj. Written.
    ChokeGas &injectionChoke;
    /// Surface choke -- SProd::chokeSup. Read only.
    const choke &surfaceChoke;
    /// Gas-line and production cell index of each gas-lift valve --
    /// SProd::posicVGLG and posicVGLP.
    const int *gasValveCellIndices;
    const int *productionValveCellIndices;

    /// Steady-state iteration counter -- SProd::iterperm. Written.
    int &steadyIteration;
    /// Which boundary the search started from -- SProd::buscaIni. Written.
    int &searchOrigin;
    /// Convergence monitors -- SProd::monitConvPerm and monitConvPermBase.
    /// Both written.
    double &convergenceMonitor;
    double &baseConvergenceMonitor;

    /// Annulus drift flag -- SProd::derivaAnel. Read only.
    const int &annulusDrift;
    /// Network coupling flag -- SProd::verificaAcop. Read only.
    const int &networkCoupled;
    /// End-node flag -- SProd::noextremo. Read only.
    const int &endNode;
    /// Thermal source switch -- SProd::semTermo. Read only.
    const int &thermalSourceDisabled;
    /// Slow-heat-transfer switch -- SProd::trocaTermicaLenta. Written.
    double &slowHeatTransfer;

    /// Surface gas pressure -- SProd::pGSup. Written.
    double &gasSurfacePressure;
    /// Previous-step gas pressure and temperature -- SProd::presiniG and
    /// tempiniG. Read only here; the gas line owns them.
    const double &initialGasPressure;
    const double &initialGasTemperature;
    /// Pressure the march ends on -- SProd::presfim. Written.
    double &finalPressure;
    /// Time step -- SProd::dt. Written: the pseudo-transient step shortens it.
    double &timeStep;

    /// Ambient temperature -- SProd::temperatura. Read only.
    const double &ambientTemperature;
    /// Casing temperature -- SProd::tempRev. Read only.
    const double &casingTemperature;
    /// Inlet quality -- SProd::titE. Read only.
    const double &inletQuality;
    /// Number of production fluids -- SProd::nfluP. Read only.
    const int &productionFluidCount;

    /// Everything the march needs that is not its own.
    SteadyStateUpdaters updaters;
};

// ------------------------------------------------------------ mass march ----

/// Advances the steady mass balance one cell. Four variants: forward and
/// reverse, each with and without the compositional treatment.
void advanceSteadyMass(const SteadyStateState &state, int cellIndex);
void advanceReverseSteadyMass(const SteadyStateState &state, int cellIndex);
void advanceCompositionalSteadyMass(const SteadyStateState &state, int cellIndex);
void advanceReverseCompositionalSteadyMass(const SteadyStateState &state, int cellIndex);

/// Mass transfer between phases along the steady march.
void advanceSteadyMassTransfer(const SteadyStateState &state, int cellIndex);
void advanceSteadyGasMassTransfer(const SteadyStateState &state, int cellIndex);

// -------------------------------------------------------- pressure march ----

/// Upstream and downstream halves of the steady pressure march, and the area
/// change between them. RK selects the Runge-Kutta stage.
void advanceUpstreamSteadyPressure(const SteadyStateState &state, int cellIndex, int rungeKuttaStage);
void advanceDownstreamSteadyPressure(const SteadyStateState &state, int cellIndex, int rungeKuttaStage);
[[nodiscard]] double areaChangePressureDrop(const SteadyStateState &state, int cellIndex,
                                            double mixtureDensity, double reynolds, double mixtureFlux);

/// Pressure at the last cell of the steady march.
[[nodiscard]] double steadyPressureAtLastCell(const SteadyStateState &state);

/// Corrects the gas specific gravity along the march.
void correctGasSpecificGravity(const SteadyStateState &state, int cellIndex);

// ------------------------------------------------------- production march ----

/// Marches the production column from a pressure guess and returns the
/// residual the search drives to zero. Forward and reverse.
[[nodiscard]] double marchProductionSteady(const SteadyStateState &state, double pressureGuess);
[[nodiscard]] double marchReverseProductionSteady(const SteadyStateState &state, double pressureGuess);
[[nodiscard]] double marchProductionSteadySecondary(const SteadyStateState &state, double pressureGuess);

/// Same column, driven from a mass-flow guess between two fixed pressures.
[[nodiscard]] double marchProductionPressureToPressure(const SteadyStateState &state, double massFlowGuess);
[[nodiscard]] double marchReverseProductionPressureToPressure(const SteadyStateState &state, double massFlowGuess);
[[nodiscard]] double marchProductionPressureToPressureSecondary(const SteadyStateState &state, double massFlowGuess);
[[nodiscard]] double marchProductionPressureToPressureTertiary(const SteadyStateState &state, double massFlowGuess);

// -------------------------------------------------------------- gas march ----

/// Marches the gas line. The mass guess defaults to -1, meaning "derive it".
[[nodiscard]] double marchGasSteady(const SteadyStateState &state, double massGuess = -1);
[[nodiscard]] double marchGasSteadySecondary(const SteadyStateState &state, double pressureGuess,
                                             double massGuess = -1);
[[nodiscard]] double marchGasSteadyTertiary(const SteadyStateState &state, double pressureGuess);

// ------------------------------------------------------- injection march ----

/// Marches an injection column from a pressure guess.
[[nodiscard]] double marchInjectionSteady(const SteadyStateState &state, double guess);

// ------------------------------------------------------------ hydrostatics --

/// Hydrostatic head walked back up a column, for production, for injection,
/// for a secondary branch and for the service line.
[[nodiscard]] double reverseHydrostatic(const SteadyStateState &state, double holdup,
                                        double liquidFlow = 0, double gasFlow = 0);
[[nodiscard]] double reverseInjectionHydrostatic(const SteadyStateState &state, double holdup,
                                                 double liquidFlow = 0);
[[nodiscard]] double secondaryBranchHydrostatic(const SteadyStateState &state, double quality);
void serviceLineHydrostatic(const SteadyStateState &state);

// ------------------------------------------------------ property refresh ----

/// Pseudo-transient time step for the steady solve.
void computePseudoTransientTimeStep(const SteadyStateState &state);

/// Refreshes fluid properties and the thermal velocities between iterations.
void refreshProperties(const SteadyStateState &state);
void refreshSteadyThermalVelocities(const SteadyStateState &state);

/// Refreshes the periphery of one production cell, upstream and downstream.
void refreshUpstreamProductionPeriphery(const SteadyStateState &state, int cellIndex);
void refreshDownstreamProductionPeriphery(const SteadyStateState &state, int cellIndex);

}  // namespace sisprod::steady

#endif  // SISPRODSTEADYSTATE_H_
