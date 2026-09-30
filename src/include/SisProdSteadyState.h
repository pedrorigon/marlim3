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
/// Fourteen of these live in other modules -- nine in sisprod::gaslift, five in
/// sisprod::thermal -- but reaching them needs a GasLiftState or a
/// ThermalState, which only the adapters build from an SProd. So the call goes
/// back through SProd, as GasLiftTemperatureUpdater does, and the adapters
/// remain the one place that knows how to assemble a module's state.
///
/// The remaining two are SProd's own: CalcC0UdPerm and renovaFonte.
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

    // --- searches, called by a march ----------------------------------------
    // These two break the shape the rest of this struct has: everything above is
    // something the march needs from another domain, while these are searches.
    // marchProductionSteady and marchProductionSteadySecondary call them after the
    // column converges, to march the gas line.
    //
    // Their return value is discarded at the call site, so these are declared
    // void.
    void searchGasPressureSteadySecondary() const;
    void searchGasPressureSteadyTertiary() const;
};

/// The state the steady march reads, and the only state it may read.
///
/// Same role as ThermalState and GasLiftState: it names in one place what this
/// domain touches, and it makes the routines callable without an SProd, so they
/// can be driven over synthetic cells. The searches read seventeen of these
/// fields, which is why SteadyStateSearchState composes this struct.
///
/// Scalars are held by reference, not by value: a copy would read every one at
/// construction, before the branch that decides whether it is read at all.
///
/// const marks what the march does not write.
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
    /// Surface choke -- SProd::chokeSup. Not const: marchProductionSteadySecondary
    /// calls vazmassSachd and vazmaxSachd on it, and neither is const-qualified.
    choke &surfaceChoke;
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
    double &slowHeatTransferThreshold;

    /// Surface gas pressure -- SProd::pGSup. Written.
    double &gasSurfacePressure;
    /// Previous-step gas pressure and temperature -- SProd::presiniG and
    /// tempiniG. Read only here; the gas line owns them.
    const double &initialGasPressure;
    const double &initialGasTemperature;
    /// Pressure the march ends on -- SProd::presfim. Written.
    double &outletPressure;
    /// Time step -- SProd::dt. Written: the pseudo-transient step shortens it.
    double &timeStep;

    /// Ambient temperature -- SProd::temperatura. Read only.
    const double &defaultInletTemperature;
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
[[nodiscard]] double marchProductionPressureToPressure(const SteadyStateState &state, double massFlowRateGuess);
[[nodiscard]] double marchReverseProductionPressureToPressure(const SteadyStateState &state, double massFlowRateGuess);
[[nodiscard]] double marchProductionPressureToPressureSecondary(const SteadyStateState &state, double massFlowRateGuess);
[[nodiscard]] double marchProductionPressureToPressureTertiary(const SteadyStateState &state, double massFlowRateGuess);

// -------------------------------------------------------------- gas march ----

/// Marches the gas line. The mass guess defaults to -1, meaning "derive it".
/// Not [[nodiscard]]: the production marches call it for its effect on the gas
/// cells and drop the value.
double marchGasSteady(const SteadyStateState &state, double massGuess = -1);
[[nodiscard]] double marchGasSteadySecondary(const SteadyStateState &state, double pressureGuess,
                                             double massGuess = -1);
[[nodiscard]] double marchGasSteadyTertiary(const SteadyStateState &state, double pressureGuess);

// ------------------------------------------------------- injection march ----

/// Marches an injection column from a pressure guess.
[[nodiscard]] double marchInjectionSteady(const SteadyStateState &state, double guess);

// ------------------------------------------------------------ hydrostatics --

/// Hydrostatic head walked back up a column, for production, for injection,
/// for a secondary branch and for the service line.
[[nodiscard]] double reverseHydrostatic(const SteadyStateState &state, double liquidHoldup,
                                        double liquidFlowRate = 0, double gasFlowRate = 0);
[[nodiscard]] double reverseInjectionHydrostatic(const SteadyStateState &state, double liquidHoldup,
                                                 double liquidFlowRate = 0);
[[nodiscard]] double secondaryBranchHydrostatic(const SteadyStateState &state, double quality);
void gasLineHydrostatic(const SteadyStateState &state);

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
