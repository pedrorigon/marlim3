#ifndef SISPRODTHERMAL_H_
#define SISPRODTHERMAL_H_

#include <vector>

class Cel;
class CelG;
class ChokeGas;
class Ler;
class SProd;
class choke;
class solverP3D;
struct varGlob1D;

namespace sisprod::thermal {

/// Direct adapter for the legacy source refresh owned by SProd.
struct ThermalSourceUpdater {
    SProd &system;

    void operator()(int cellIndex) const;
};

/// Direct adapter for the four legacy drift-closure entry points owned by SProd.
struct ThermalClosureUpdater {
    SProd &system;

    void instantaneous(int cellIndex, double &distribution,
                       double &driftVelocity) const;
    void buffered(int cellIndex, double &distribution,
                  double &driftVelocity) const;
    void initialization(int cellIndex, double &distribution,
                        double &driftVelocity) const;
    void bufferedInitialization(int cellIndex, double &distribution,
                                double &driftVelocity) const;
};

/// Direct adapter for the legacy transient coupling owned by SProd.
struct ThermalEvolutionUpdater {
    SProd &system;

    void solvePressureVelocityCoupling(int cycle) const;
    void renew() const;
};

/// The slice of SProd's state this module reads, so it can be called without
/// access to all of SProd.
struct ThermalState {
    Cel *cells;
    CelG *gasCells;
    const Ler &input;
    double **latentHeatTable;
    varGlob1D *globals;
    const int &thermalSourceDisabled;
    const int &productionNetworkCoupled;
    const int &primarySectionStart;
    const int &primarySectionEnd;
    const std::vector<int> &coupledCellIndices;
    solverP3D &poissonSolver;
    const int &lastCell;
    const choke &surfaceChoke;
    const int &surfaceChokeMassCondition;
    const int &networkEndpoint;
    const double &gasSurfaceTemperature;
    const int &latentHeatEnabled;
    ThermalSourceUpdater sourceUpdater;
    const int &completeModel;
    const int &massTransferModel;
    ThermalClosureUpdater closureUpdater;
    const double &inletPressure;
    const double &inletTemperature;
    const double &inletMassFraction;
    double &inletVoidFraction;
    const double &inletComposition;
    ThermalEvolutionUpdater evolutionUpdater;
    const int &surfaceChokeOpen;
    const double &defaultInletTemperature;
    const double &timeStep;
    const double &minimumCycleTimeStep;
    const std::vector<int> &poisson2DCellIndices;
    const int &poisson2DCellCount;
    const int &steadyIteration;
    const double &slowHeatTransferThreshold;
    const int &annulusTubingStart;
    const int &annulusTubingEnd;
    const int &tubingAnnulusStart;
    const int &productionNetworkHeatCoupled;
    const int &primaryNetworkSectionEnd;
    const int &primaryNetworkSectionStart;
    const int &gasCellCount;
    ChokeGas *gasLiftChokes;
    const double &gasSurfacePressure;
    const double &outletPressure;
    /// Written, not just read -- computeOutletTemperature assigns it.
    double &surfaceTemperature;
    const int &networkCoupled;
    const int &tubingAnnulusEnd;
};

/// Interpolates latent heat in the pressure-temperature table.
double interpolateLatentHeat(const ThermalState &state, double pressure,
                             double temperature);

/// Computes the mixture enthalpy balance for one control volume.
double computeMixtureEnthalpy(const ThermalState &state, int cellIndex);

/// Interpolates mixture energy between two pressure rows of a property table.
double interpolateMixtureEnergy(const ThermalState &state, int cellIndex,
                                int pressureIndex, int temperatureIndex,
                                double pressureRatio);

/// Updates one control-volume temperature from tabulated mixture enthalpy.
void updateTemperatureFromEnthalpy(const ThermalState &state, int cellIndex);

/// Computes the legacy thermal phase-change mass source for one cell.
void computeThermalMassTransfer(const ThermalState &state, int cellIndex);

/// Updates one control-volume temperature from the thermal energy balance.
void computeTemperature(const ThermalState &state, int cellIndex,
                        double previousTemperature, int steadyStateMode = 0);

/// Renews distributed phase-change mass transfer along the production cells.
void updateDistributedMassTransfer(const ThermalState &state);

/// Computes the mixture-flow partition terms along the production cells.
void updateFlowPartitionTerms(const ThermalState &state, int inflowMode);

/// Computes mixture-flow partition terms at an internal-section outlet.
void updateOutletFlowPartitionTerms(const ThermalState &state);

/// Computes mixture-flow partition terms at an internal-section inlet.
void updateInletFlowPartitionTerms(const ThermalState &state);

/// Prepares the non-dimensional heat-diffusion properties for one cell.
void prepareNonDimensionalHeatDiffusion(const ThermalState &state,
                                        int cellIndex);

/// Advances the transient thermal-energy solution by one coupling cycle.
void advanceTransientEnergy(const ThermalState &state, int cycle,
                            int maximumCycle);

/// Advances the legacy steady-state temperature march in the direct direction.
void advanceSteadyTemperature(const ThermalState &state, int cellIndex,
                              int rungeKuttaStage);

/// Advances the distinct legacy steady-state temperature march in reverse.
void advanceReverseSteadyTemperature(const ThermalState &state, int cellIndex,
                                     int rungeKuttaStage);

/// Solves the gas-line temperature for one cell.
void computeGasTemperature(const ThermalState &state, int cellIndex,
                           double previousTemperature, int steadyMode);

/// Temperature at the discharge of an accessory.
void computeDischargeTemperature(const ThermalState &state, int cellIndex);

/// Discharge temperature of a gas-lift valve.
double computeGasLiftDischargeTemperature(const ThermalState &state,
                                          int valveIndex);

/// Propagates the cell temperature to its neighbours' interface fields.
void updateProductionTemperaturePeriphery(const ThermalState &state,
                                          int cellIndex);

/// Surface temperature at the end of the production line.
void computeOutletTemperature(const ThermalState &state);


// Values passed between the kernels inside SisProdThermal.cpp. Not part of the
// module's interface -- nothing outside the .cpp constructs or reads them.

/// The enthalpy the source term carries into the cell, and the temperature it
/// arrives at. Produced by the accessory dispatch in sourceEnthalpyOf.
struct SourceEnthalpy {
    double temperature;
    double gasEnthalpy;
    double liquidEnthalpy;
    /// Complementary-fluid heat. The legacy name is kept: three assignment
    /// sites carry a comment saying it is still to be corrected.
    double hcF;
};

/// Cell geometry and superficial velocities. They come from this cell unless the
/// upstream neighbour is a throttled choke, in which case from the downstream one.
struct CellFlowBasis {
    double flowArea;
    double voidFraction;
    double betmed;
    double gasSuperficialVelocity;
    double liquidSuperficialVelocity;
};

struct TemperatureBalance {
    double cellLength;
    double meanCellLength;
    double flowArea;
    double voidFraction;
    double bet;
    double gasSuperficialVelocity;
    double liquidSuperficialVelocity;
    double referenceMixtureVelocity;
    double rp;
    double rc;
    double liquidDensity;
    double gasDensity;
    double liquidHeatCapacity;
    double liquidIsochoricHeatCapacity;
    double gasHeatCapacity;
    double gasIsochoricHeatCapacity;
    double liquidJouleThomson;
    double gasJouleThomson;
    double hydrostaticPower;
    double heatFlux;
    double timeCoefficient;
    double pressureTimeCoefficient;
    double temperatureSpatialCoefficient;
    double pressureSpatialCoefficient;
};

struct TemperatureSourceTerms {
    double gas;
    double liquid;
};

struct DistributedMassTransferProperties {
    double downstreamWaterFraction;
    double upstreamWaterFraction;
    double cellWaterFraction;
    double liquidDensity;
    double gasDensity;
    double downstreamComposition;
    double upstreamComposition;
    double mixtureLiquidDensity;
    double downstreamOilVolumeFactor;
    double downstreamSolutionGasRatio;
    double downstreamSolutionGasPressureDerivative;
    double cellOilVolumeFactor;
    double cellSolutionGasRatio;
    double cellSolutionGasPressureDerivative;
    double cellSolutionGasTemperatureDerivative;
};

struct DistributedMassTransferCoefficients {
    double activeDerivative;
    double spatialCoupling;
    double flowArea;
};

/// The drift-flux pair of the slug regime: C0 = 1.2 and ud = 0.32*sqrt(g*D),
/// signed by the inclination, collapsing to the homogeneous limit (C0 = 1,
/// ud = 0) once the density ratio passes 0.9.
struct SlugClosure {
    double c0;
    double ud;
    /// Read after the call by the interior selector only. Returned rather than
    /// passed by reference so the other three callers do not carry a variable
    /// they never read, which would trade duplication for a warning.
    double meanDiameter;
};

/// Which face the upstream properties come from, decided by the sign of the
/// gas flow rate, for the two drift-closure selectors that share the choice.
struct UpstreamFaceBasis {
    double gasDensity;
    double flowArea;
    double noSlipLiquidHoldup;
};

}  // namespace sisprod::thermal

#endif  // SISPRODTHERMAL_H_
