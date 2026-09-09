#ifndef SISPRODGASLIFT_H_
#define SISPRODGASLIFT_H_

#include <vector>

// Declared, not included: GasLiftState holds only pointers and references to
// these, so this header stays free of the cell, input-deck, choke and matrix
// headers and can still be compiled on its own.
class Cel;
class CelG;
class ChokeGas;
class Ler;
struct varGlob1D;
class SProd;
template <class T> class Vcr;
template <class T> class BandMtx;

namespace sisprod::gaslift {

/// Direct adapter for the discharge-temperature step owned by SProd.
///
/// SProd::tempDescarga already forwards to the thermal module, but reaching it
/// needs a ThermalState the gas line does not carry. Same shape the thermal
/// module uses for the closures it does not own.
struct GasLiftTemperatureUpdater {
    SProd &system;

    void dischargeTemperature(int cellIndex) const;
    /// steadyMode defaults to 0, matching SProd::calctempGas: callers in the
    /// moved bodies omit it.
    void gasTemperature(int cellIndex, double previousTemperature,
                        int steadyMode = 0) const;
};

/// The state the gas line and gas-lift routines read, and the only state they
/// may read.
///
/// Same role as ThermalState in the thermal module: it names in one place what
/// this domain is allowed to touch, and it makes the routines callable WITHOUT
/// an SProd. That second property is not tidiness -- it is what lets a dedicated
/// harness drive them over synthetic cells, which for this stage is the only
/// verification that reaches twelve of the twenty-two functions. The demo corpus
/// never executes those, so the artifact and regression layers stay green for
/// them whatever an extraction does.
///
/// Scalars are held BY REFERENCE, not by value. Copying them into the struct
/// would read every one at construction, before the branch that decides whether
/// the original would have read it at all -- the defect the root-finding stage
/// had to undo, and the reason DriftFluxClosure's ClosureState carries the same
/// warning.
struct GasLiftState {
    /// Gas-line cells -- SProd::celulaG. Written as well as read.
    CelG *gasCells;
    /// Production cells -- SProd::celula. The gas line reads the tubing it
    /// feeds, and writes back at the connection points.
    Cel *cells;
    /// Input deck -- SProd::arq. NOT const: the unloading schedule uses
    /// presMaxDesc as scratch, recomputing it as a minimum over the IPR
    /// accessories and writing it back. Declaring it const would have been a
    /// claim this module does not honour.
    Ler &input;
    /// Shared 1D globals -- SProd::vg1dSP.
    varGlob1D *globals;

    /// Index of the last gas-line cell -- SProd::ncelGas.
    const int &gasCellCount;
    /// Index of the last production cell -- SProd::ncel.
    const int &lastCell;

    /// Gas-lift valve chokes -- SProd::chokeVGL.
    ChokeGas *gasLiftChokes;
    /// Injection choke -- SProd::chokeInj.
    ChokeGas &injectionChoke;
    /// Gas-line cell index of each gas-lift valve -- SProd::posicVGLG.
    const int *gasValveCellIndices;
    /// Production cell index of each gas-lift valve -- SProd::posicVGLP.
    const int *productionValveCellIndices;

    /// Band matrix and free-term vector of the gas line -- SProd::matglobG and
    /// SProd::termolivreG.
    BandMtx<double> &gasSystemMatrix;
    Vcr<double> &gasFreeTerms;

    /// Annulus/tubing coupling bounds, shared with the thermal module and named
    /// as they are named there -- SProd::ColunaAnulaIni, ColunaAnulaFim,
    /// AnulaColunaIni, AnulaColunaFim.
    const int &annulusTubingStart;
    const int &annulusTubingEnd;
    const int &tubingAnnulusStart;
    const int &tubingAnnulusEnd;

    /// Steady-state iteration counter -- SProd::iterperm.
    const int &steadyIteration;
    /// Network coupling flag -- SProd::verificaAcop.
    const int &networkCoupled;
    /// Thermal source switch -- SProd::semTermo.
    const int &thermalSourceDisabled;

    /// Previous-step gas pressure and temperature -- SProd::presiniG and
    /// SProd::tempiniG. The pressure is written.
    double &initialGasPressure;
    const double &initialGasTemperature;
    /// Surface gas pressure -- SProd::pGSup. Written.
    double &gasSurfacePressure;

    /// Time step -- SProd::dt. Written: the gas line can shorten it.
    double &timeStep;

    /// Unloading interface state -- SProd::celInter, velInter, dtInter, and
    /// their initial counterparts. All written.
    int &interfaceCell;
    double &interfaceVelocity;
    double &interfaceTimeStep;
    int &initialInterfaceCell;
    double &initialInterfaceVelocity;
    double &initialInterfaceTimeStep;

    /// Unloading averages -- SProd::vazmedDesc, tempmedDEsc (spelling as in
    /// SProd), and the bounds they are compared against.
    ///
    /// The two vectors are NOT const: they are sliding windows. advanceGasSubStep
    /// push_backs the current step at the tail and erases the front once the
    /// window passes maximumContinuousUnloadingCount. Declaring them const was a
    /// claim this module does not honour.
    double &meanUnloadingFlowRate;
    double &meanUnloadingTemperature;
    std::vector<double> &maximumMeanUnloadingFlowRates;
    std::vector<double> &unloadingTimeSteps;
    const double &continuousMeanUnloadingTemperature;
    const double &maximumContinuousUnloadingCount;

    /// Discharge temperature, computed by the thermal module.
    GasLiftTemperatureUpdater temperatureUpdater;
};

/// Hydrostatics of the gas column during unloading.
void computeGasUnloadingHydrostatics(const GasLiftState &state);

/// Advances the gas line by one sub-time-step.
void advanceGasSubStep(const GasLiftState &state);

/// Same, for the buffered variant. No call site exists in the product.
void advanceBufferedGasSubStep(const GasLiftState &state);

/// Updates the gas-line cells after a sub-step.
void updateGasLine(const GasLiftState &state);

/// Same, buffered. Reachable only from advanceBufferedGasSubStep, which nothing
/// calls, so this is dead code moved for completeness.
void updateBufferedGasLine(const GasLiftState &state);

/// Throat area of a calibrated gas-lift valve.
double calibratedValveArea(double calibrationPressure, double calibrationTemperature,
                           double valveOpeningPressure, double tubingPressure,
                           double externalDiameter, double throatArea,
                           double valveRatio, double temperature);

/// Pressure correction applied to the unloading flow of one valve.
double unloadingPressureCorrection(const GasLiftState &state, double maximumFlowRate,
                                   int valveIndex, double factor, int sign);

/// Pressure at an unloading valve.
double computeUnloadingValvePressure(const GasLiftState &state,
                                     double throatFlowRate, int valveIndex);

/// Searches the injection pressure that satisfies the unloading schedule,
/// writing it into the gas surface and initial pressures and into each valve's
/// stage and throat pressures. Returns the maximum unloading flow rate the
/// search reached -- the caller keeps it as velmaxdesc.
double searchUnloadingInjectionPressure(const GasLiftState &state);

/// Solves one unloading step.
void solveUnloading(const GasLiftState &state);

/// Advances the unloading interface.
void advanceInterface(const GasLiftState &state);

/// Updates the transient gas-lift valves.
void updateTransientGasValves(const GasLiftState &state);

/// Solves the gas line.
void solveGasLine(const GasLiftState &state);

/// Connects the gas line to the tubing.
void connectTubing(const GasLiftState &state);

/// Steady-state counterparts.
void connectTubingSteady(const GasLiftState &state);
void initializeTubingConnectionSteady(const GasLiftState &state);
double steadyGasPressureDrop(const GasLiftState &state, int cellIndex);
double steadyInjectionPressureDrop(const GasLiftState &state, int cellIndex);
void updateSteadyGasPressure(const GasLiftState &state, int cellIndex);
void computeSteadyGasFlowRate(const GasLiftState &state, int cellIndex);
void initializeSteadyValveGasFlowRate(const GasLiftState &state, int cellIndex);
void updateSteadyGasTemperature(const GasLiftState &state, int cellIndex);

}  // namespace sisprod::gaslift

#endif  // SISPRODGASLIFT_H_
