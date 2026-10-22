#ifndef SISPRODTHERMAL_H_
#define SISPRODTHERMAL_H_

#include <vector>

class Cel;
class CelG;
class ChokeGas;
class Ler;
class SProd;
namespace sisprod {
class SolveContext;
}
class choke;
class solverP3D;
struct varGlob1D;

namespace sisprod::thermal {

/// The source refresh the thermal steps call, with the solve context's source view.
struct ThermalSourceUpdater {
    SolveContext &context;

    /// Refreshes mass and heat sources for the selected cell.
    void operator()(int cellIndex) const;
};

/// The four drift-closure entry points the thermal steps call, with the solve context's
/// closure view.
struct ThermalClosureUpdater {
    SolveContext &context;

    /// Evaluates the instantaneous drift closure.
    void instantaneous(int cellIndex, double &distribution,
                       double &driftVelocity) const;
    /// Evaluates the buffered drift closure.
    void buffered(int cellIndex, double &distribution,
                  double &driftVelocity) const;
    /// Evaluates the instantaneous inlet closure.
    void initialization(int cellIndex, double &distribution,
                        double &driftVelocity) const;
    /// Evaluates the buffered inlet closure.
    void bufferedInitialization(int cellIndex, double &distribution,
                                double &driftVelocity) const;
};

/// The transient coupling the thermal steps call, with the solve context's view of the
/// transient step.
struct ThermalEvolutionUpdater {
    SolveContext &context;

    /// Solves the pressure/velocity coupling for one cycle.
    void solvePressureVelocityCoupling(int cycle) const;
    /// Refreshes the transient step state after thermal coupling.
    void renew() const;
};

/// State and callbacks read and updated by the thermal and phase-transfer steps.
struct ThermalState {
    Cel *const &cells;
    CelG *const &gasCells;
    const Ler &input;
    /// The latent heat table's row pointers -- SProd::tables.latentHeatTable.
    const std::vector<double *> &latentHeatTable;
    varGlob1D *const &globals;
    const bool &thermalSourceDisabled;
    const bool &productionNetworkCoupled;
    const int &primarySectionStart;
    const int &primarySectionEnd;
    const std::vector<int> &coupledCellIndices;
    solverP3D &poissonSolver3D;
    const int &lastCell;
    const choke &surfaceChoke;
    const int &surfaceChokeMassFlag;
    const int &endNode;
    const double &gasSurfaceTemperature;
    const int &latentHeatEnabled;
    ThermalSourceUpdater sourceUpdater;
    const int &fullModel;
    const int &massTransferModel;
    ThermalClosureUpdater closureUpdater;
    const double &inletPressure;
    const double &inletTemperature;
    const double &inletQuality;
    double &inletVoidFraction;
    const double &inletComplementaryFraction;
    ThermalEvolutionUpdater evolutionUpdater;
    const int &surfaceChokeOpen;
    const double &defaultInletTemperature;
    const double &timeStep;
    const double &minimumCycleTimeStep;
    const std::vector<int> &poisson2DCellIndices;
    int poisson2DCellCount;
    const int &steadyIteration;
    const double &slowHeatTransferThreshold;
    const int &annulusTubingStart;
    const int &annulusTubingEnd;
    const int &tubingAnnulusStart;
    const bool &productionNetworkHeatCoupled;
    const int &primaryNetworkSectionEnd;
    const int &primaryNetworkSectionStart;
    const int &lastGasCell;
    std::vector<ChokeGas> &gasLiftChokes;
    const double &gasSurfacePressure;
    const double &outletPressure;
    /// Written, not just read -- computeOutletTemperature assigns it.
    double &surfaceTemperature;
    const bool &networkCoupled;
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

/// Computes the thermal phase-change mass source for one cell.
void computeThermalMassTransfer(const ThermalState &state, int cellIndex);

/// Updates one control-volume temperature from the thermal energy balance.
/// previousTemperature does not participate in the calculation.
/// Uses previous pressure, current temperature, and prior face temperatures.
void computeTemperature(const ThermalState &state, int cellIndex,
                        double previousTemperature, int steadyStateMode = 0);

/// Renews distributed phase-change mass transfer along the production cells.
void updateDistributedMassTransfer(const ThermalState &state);

/// Computes the mixture-flow partition terms along the production cells.
/// inflowMode does not affect the sweep, which uses the outlet boundary.
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

/// Advances the steady-state temperature march in the direct direction.
void advanceSteadyTemperature(const ThermalState &state, int cellIndex,
                              int rungeKuttaStage);

/// Advances the reverse steady temperature march with its directional source
/// and coupling rules.
void advanceReverseSteadyTemperature(const ThermalState &state, int cellIndex,
                                     int rungeKuttaStage);

/// Solves the gas-line temperature for one cell and publishes face temperatures.
/// previousTemperature does not participate in the calculation.
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

} // namespace sisprod::thermal

#endif // SISPRODTHERMAL_H_
