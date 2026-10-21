#ifndef SISPRODGASLIFT_H_
#define SISPRODGASLIFT_H_

#include <vector>

// Declared, not included: GasLiftState holds only references to these, so this
// header stays free of the cell, input-deck, choke and matrix headers and can
// still be compiled on its own.
class Cel;
class CelG;
class ChokeGas;
class Ler;
struct varGlob1D;
namespace sisprod { class SolveContext; }
template <class T> class Vcr;
template <class T> class BandMtx;

namespace sisprod::gaslift {

/// The thermal steps the gas line calls: the discharge, gas and gas-lift discharge
/// temperatures. They take a ThermalState the gas line does not carry, which the solve
/// context provides.
struct GasLiftTemperatureUpdater {
    SolveContext &context;

    void dischargeTemperature(int cellIndex) const;
    void computeGasTemperature(int cellIndex, double previousTemperature, int steadyMode = 0) const;
    [[nodiscard]] double gasLiftDischargeTemperature(int valveIndex) const;
};

/// The state the gas line and gas-lift routines read, and the only state they
/// may read.
///
/// Scalars are held by reference, not by value: a copy would read every one at
/// construction, before the branch that decides whether it is read at all.
struct GasLiftState {
    /// Gas-line cells, written as well as read.
    CelG *const &gasCells;
    /// Production cells. The gas line reads the tubing it feeds, and writes back
    /// at the connection points.
    Cel *const &cells;
    /// The parsed input. Not const: the unloading schedule uses
    /// presMaxDesc as scratch, recomputing it as a minimum over the IPR
    /// accessories and writing it back, and clears descarga when the
    /// interface reaches the end of the gas line.
    Ler &input;
    /// Shared 1D globals.
    varGlob1D *const &globals;

    /// Index of the last gas-line cell.
    const int &lastGasCell;
    /// Index of the last production cell.
    const int &lastCell;

    /// Gas-lift valve chokes.
    std::vector<ChokeGas> &gasLiftChokes;
    /// Injection choke.
    ChokeGas &injectionChoke;
    /// Gas-line cell index of each gas-lift valve.
    const std::vector<int> &gasValveCellIndices;
    /// Production cell index of each gas-lift valve.
    const std::vector<int> &productionValveCellIndices;

    /// Band matrix and right-hand side of the gas-line pressure-velocity system.
    BandMtx<double> &gasSystemMatrix;
    Vcr<double> &gasRightHandSide;

    /// Bounds of the column-annulus thermal coupling: tubingCouplingStart and
    /// tubingCouplingEnd are production-column cells, annulusCouplingStart and
    /// annulusCouplingEnd gas-line cells.
    const int &tubingCouplingStart;
    const int &tubingCouplingEnd;
    const int &annulusCouplingStart;
    const int &annulusCouplingEnd;

    /// Whether the column-annulus thermal coupling is on.
    const bool &annulusThermalCoupling;
    /// When set, the steady gas temperature is the external one, neither marched nor coupled to the tubing.
    const bool &thermalSourceDisabled;

    /// Gas-injection pressure and temperature. The pressure is written.
    double &injectionPressure;
    const double &injectionTemperature;
    /// Separator pressure, or the inlet pressure of the downstream network section. Written.
    double &separatorPressure;

    /// Time step. Written: the gas line can shorten it.
    double &timeStep;

    /// The completion-fluid/gas interface: its cell, velocity and time step, and their values at
    /// the previous time level. All written.
    int &interfaceCell;
    double &interfaceVelocity;
    double &interfaceTimeStep;
    int &previousInterfaceCell;
    double &previousInterfaceVelocity;
    double &previousInterfaceTimeStep;

    /// The unloading controller's moving average: the accumulated flow and the averaging interval,
    /// the samples they sum, and the window and the sample capacity that bound them.
    ///
    /// The two vectors are sliding windows: searchUnloadingInjectionPressure appends each step's
    /// sample, and drops the oldest when they hold more samples than the capacity or span more time
    /// than the window.
    double &meanUnloadingFlowRate;
    double &unloadingAveragingInterval;
    std::vector<double> &maximumUnloadingFlowSamples;
    std::vector<double> &unloadingTimeStepSamples;
    const double &unloadingAveragingWindow;
    const int &unloadingSampleCapacity;

    /// The gas line's temperatures, which the thermal module computes.
    GasLiftTemperatureUpdater temperatureUpdater;
};

// Unloading the completion fluid through the gas-lift valves.

/// Gives the gas line its initial unloading state: the pressure marched cell by cell from the start pressure through
/// the hydrostatic head, of gas before the interface cell and of completion fluid from it on; each cell's external
/// temperature; no flow and no source. A controlled unloading also takes its initial gas pressure as the injection
/// pressure.
void computeGasUnloadingHydrostatics(const GasLiftState &state);

/// Gives the production column its initial unloading state: the pressure marched cell by cell from the first cell's
/// IPR reservoir pressure through the hydrostatic head, of completion fluid up to the unloading cell and of gas above
/// it; each cell's external temperature; the matching phases; no flow and no mass source.
void computeProductionUnloadingHydrostatics(const GasLiftState &state);

/// The opening of a calibrated gas-lift valve as a fraction of its throat area, at most 1; a closed valve is not a
/// special case. The correlation is imperial: the calibration pressure, the valve opening pressure (the casing's) and
/// the tubing pressure are gauge pressures in psi, and the bottom-hole temperature is in degrees Fahrenheit; the
/// calibration temperature comes in degrees Celsius and the external diameter in metres. valveRatio is the throat
/// area over the bellows area.
[[nodiscard]] double calibratedValveOpeningFraction(double calibrationPressure, double calibrationTemperature,
                                                    double valveOpeningPressure, double tubingPressure,
                                                    double externalDiameter, double throatArea, double valveRatio,
                                                    double bottomHoleTemperatureFahrenheit);

/// The pressure correction the unloading control applies for one valve, in kgf/cm^2 and with the given sign: the
/// dynamic pressure, at the valve's throat, of the completion-fluid flow that separates maximumFlowRate from factor
/// times the controlled unloading flow rate. Primes the valve's choke with its gas cell's pressure and temperature.
[[nodiscard]] double unloadingPressureCorrection(const GasLiftState &state, double maximumFlowRate, int valveIndex,
                                                 double factor, int sign);

/// Recomputes presMaxDesc along the production column and moves the injection pressure towards the unloading
/// flow band.
void computeUnloadingValvePressure(const GasLiftState &state, double throatFlowRate, int valveIndex);

/// Searches the injection pressure that satisfies the unloading schedule, writing it into the separator and
/// injection pressures and into the stagnation and throat pressures of the valves behind the interface. Returns the
/// maximum unloading flow rate the search reached.
[[nodiscard]] double searchUnloadingInjectionPressure(const GasLiftState &state);

/// Solves one unloading step from the interface cell to the end of the gas line: marches the pressure through the
/// friction and the head of each half cell's gas and completion-fluid parts, updates the discharge temperatures, and
/// sets the flows that carry the valves' sources along the line and the interface velocity.
void solveUnloading(const GasLiftState &state);

/// Advances the unloading interface over the time step. When it crosses a face of its cell, the interface time step
/// is shortened to the crossing and the interface moves to the next or the previous cell; at the end of the gas line
/// it is handed over and the unloading ends (interfaceCell becomes kNoInterfaceCell, input.descarga 0).
void advanceInterface(const GasLiftState &state);

// The gas line in the transient.

/// Advances the gas line by one sub-step in the transient: updates the valves, records on each valve's production
/// cell what the valve injects (gas ahead of the interface, completion fluid behind it), at what rate and temperature,
/// copies the valve's state onto its gas cell, then advances the sub-step. Does nothing without a gas line.
void solveGasLine(const GasLiftState &state);

/// Clears the gas cells' sources from the interface cell on, primes each gas-lift valve's choke and sets the mass
/// source it gives its gas cell (completion fluid behind the interface, gas ahead of it), then the sources of the
/// master2 valve's two cells, which only flow while it throttles.
void updateTransientGasValves(const GasLiftState &state);

/// Advances the gas line by one sub-time-step. During the unloading it first advances the interface, shortening the
/// time step of the gas line and of the production cells when the interface needs a shorter one; then, for each
/// thermal coupling cycle, it solves the pressure-flow system and the gas temperatures; last come the gas densities
/// and, during the unloading, the unloading step.
void advanceGasSubStep(const GasLiftState &state);

/// Opens the injection choke if it is throttling, assembles the band system from every gas cell and solves it.
void assembleAndSolveGasSystem(const GasLiftState &state);

/// Copies the solved pressure-flow system into the gas cells: each cell's pressure and right-face mass flow, and the
/// values its neighbours' faces hold of them.
void updateGasLine(const GasLiftState &state);

/// Couples the gas line to the production tubing in the transient: each coupled production cell takes the facing gas
/// cell's temperatures, velocity, diameters and fluid (gas ahead of the interface, completion fluid behind it) as its
/// outside, and the gas cell takes the temperature of the tubing's outermost layer.
void connectTubing(const GasLiftState &state);

// The gas line in the steady state.

/// Connects the gas line to the tubing at steady state: each coupled production cell takes the facing gas cell's
/// temperature, velocity and transport properties as its external fluid.
void connectTubingSteady(const GasLiftState &state);

/// Connects the gas line to the tubing before the first steady march iteration, from the facing production cell's
/// pressure, the deck's external temperature and a fixed gas viscosity.
void initializeTubingConnectionSteady(const GasLiftState &state);

/// The friction and hydrostatic pressure drop across a gas-line cell's two half cells at steady state, in kgf/cm^2.
[[nodiscard]] double steadyGasPressureDrop(const GasLiftState &state, int cellIndex);

/// The friction and hydrostatic pressure drop across the two half cells of an injection well's production cell at
/// steady state, in kgf/cm^2.
[[nodiscard]] double steadyInjectionPressureDrop(const GasLiftState &state, int cellIndex);

/// Marches the steady gas pressure into a cell from its left neighbour, records on the neighbour the friction and
/// hydrostatic gradients of the march's first half, and copies the pressure to the neighbours' faces.
void updateSteadyGasPressure(const GasLiftState &state, int cellIndex);

/// The steady gas flow through a cell's right face: the left face's, less what a gas-lift valve in the cell
/// delivers to the tubing. For each valve in the cell it also primes the choke and records the injection on the
/// valve's production cell; when the flow runs out before the end of the line, the valves after it get no source.
void computeSteadyGasFlowRate(const GasLiftState &state, int cellIndex);

/// Primes the steady gas mass source of the gas-lift valves that discharge into this production cell, and records
/// the injection on the cell.
void initializeSteadyValveGasFlowRate(const GasLiftState &state, int cellIndex);

/// The steady gas temperature of a cell: marched from its left neighbour, or coupled to the tubing inside the
/// column-annulus coupling.
void updateSteadyGasTemperature(const GasLiftState &state, int cellIndex);

}  // namespace sisprod::gaslift

#endif  // SISPRODGASLIFT_H_
