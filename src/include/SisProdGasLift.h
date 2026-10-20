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
    const int &gasCellCount;
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
    Vcr<double> &gasFreeTerms;

    /// Bounds of the column-annulus thermal coupling: annulusTubingStart and
    /// annulusTubingEnd are production-column cells, tubingAnnulusStart and
    /// tubingAnnulusEnd gas-line cells.
    const int &annulusTubingStart;
    const int &annulusTubingEnd;
    const int &tubingAnnulusStart;
    const int &tubingAnnulusEnd;

    /// Whether the column-annulus thermal coupling is on.
    const bool &networkCoupled;
    /// Thermal source switch.
    const bool &thermalSourceDisabled;

    /// Gas-injection pressure and temperature. The pressure is written.
    double &initialGasPressure;
    const double &initialGasTemperature;
    /// Separator pressure, or the inlet pressure of the downstream network section. Written.
    double &gasSurfacePressure;

    /// Time step. Written: the gas line can shorten it.
    double &timeStep;

    /// The completion-fluid/gas interface: its cell, velocity and time step, and their values at
    /// the previous time level. All written.
    int &interfaceCell;
    double &interfaceVelocity;
    double &interfaceTimeStep;
    int &initialInterfaceCell;
    double &initialInterfaceVelocity;
    double &initialInterfaceTimeStep;

    /// The unloading controller's moving average: the accumulated flow and the averaging interval,
    /// the samples they sum, and the window and the sample capacity that bound them.
    ///
    /// The two vectors are sliding windows: searchUnloadingInjectionPressure appends each step's
    /// sample, and drops the oldest when they hold more samples than the capacity or span more time
    /// than the window.
    double &meanUnloadingFlowRate;
    double &meanUnloadingTemperature;
    std::vector<double> &maximumMeanUnloadingFlowRates;
    std::vector<double> &unloadingTimeSteps;
    const double &continuousMeanUnloadingTemperature;
    const int &maximumContinuousUnloadingCount;

    /// Discharge temperature, computed by the thermal module.
    GasLiftTemperatureUpdater temperatureUpdater;
};

// Unloading the completion fluid through the gas-lift valves.

/// Hydrostatics of the gas column during unloading.
void computeGasUnloadingHydrostatics(const GasLiftState &state);

/// Hydrostatics of the production column during unloading.
void computeProductionUnloadingHydrostatics(const GasLiftState &state);

/// Throat area of a calibrated gas-lift valve.
[[nodiscard]] double calibratedValveArea(double calibrationPressure, double calibrationTemperature,
                                         double valveOpeningPressure, double tubingPressure,
                                         double externalDiameter, double throatArea, double valveRatio,
                                         double bottomHoleTemperatureFahrenheit);

/// Pressure correction applied to the unloading flow of one valve.
[[nodiscard]] double unloadingPressureCorrection(const GasLiftState &state, double maximumFlowRate, int valveIndex,
                                                 double factor, int sign);

/// Recomputes presMaxDesc along the production column and moves the injection pressure towards the unloading
/// flow band.
void computeUnloadingValvePressure(const GasLiftState &state, double throatFlowRate, int valveIndex);

/// Searches the injection pressure that satisfies the unloading schedule,
/// writing it into the gas surface and initial pressures and into each valve's
/// stage and throat pressures. Returns the maximum unloading flow rate the
/// search reached -- the caller keeps it as velmaxdesc.
[[nodiscard]] double searchUnloadingInjectionPressure(const GasLiftState &state);

/// Solves one unloading step.
void solveUnloading(const GasLiftState &state);

/// Advances the unloading interface.
void advanceInterface(const GasLiftState &state);

// The gas line in the transient.

/// Solves the gas line.
void solveGasLine(const GasLiftState &state);

/// Updates the transient gas-lift valves.
void updateTransientGasValves(const GasLiftState &state);

/// Advances the gas line by one sub-time-step.
void advanceGasSubStep(const GasLiftState &state);

/// Updates the gas-line cells after a sub-step.
void updateGasLine(const GasLiftState &state);

/// Solves the gas line once and keeps the solved velocities in the buffer, without advancing the cells.
void advanceBufferedGasSubStep(const GasLiftState &state);

/// Copies the solved gas velocities into the buffer.
void updateBufferedGasLine(const GasLiftState &state);

/// Connects the gas line to the tubing.
void connectTubing(const GasLiftState &state);

// The gas line in the steady state.

void connectTubingSteady(const GasLiftState &state);
void initializeTubingConnectionSteady(const GasLiftState &state);
[[nodiscard]] double steadyGasPressureDrop(const GasLiftState &state, int cellIndex);
[[nodiscard]] double steadyInjectionPressureDrop(const GasLiftState &state, int cellIndex);
void updateSteadyGasPressure(const GasLiftState &state, int cellIndex);
void computeSteadyGasFlowRate(const GasLiftState &state, int cellIndex);
void initializeSteadyValveGasFlowRate(const GasLiftState &state, int cellIndex);
void updateSteadyGasTemperature(const GasLiftState &state, int cellIndex);

}  // namespace sisprod::gaslift

#endif  // SISPRODGASLIFT_H_
