#ifndef SISPRODTHERMAL_H_
#define SISPRODTHERMAL_H_

#include <vector>

class Cel;
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

/// State read by the first thermal kernels extracted from SProd.
///
/// This context grows only when a moved function demonstrates another state
/// dependency. Keeping the cell array, input deck and latent-heat table named
/// here lets the module be called without granting it access to all of SProd.
struct ThermalState {
    Cel *cells;
    const Ler &input;
    double **latentHeatTable;
    varGlob1D *globals;
    const int &thermalSourceDisabled;
    const int &productionNetworkCoupled;
    const int &primarySectionStart;
    const int &primarySectionEnd;
    const std::vector<int> &coupledCellIndices;
    const solverP3D &poissonSolver;
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
                        double previousTemperature, int steadyStateMode);

/// Renews distributed phase-change mass transfer along the production cells.
void updateDistributedMassTransfer(const ThermalState &state);

/// Computes the mixture-flow partition terms along the production cells.
void updateFlowPartitionTerms(const ThermalState &state, int inflowMode);

/// Computes mixture-flow partition terms at an internal-section outlet.
void updateOutletFlowPartitionTerms(const ThermalState &state);

/// Computes mixture-flow partition terms at an internal-section inlet.
void updateInletFlowPartitionTerms(const ThermalState &state);

}  // namespace sisprod::thermal

#endif  // SISPRODTHERMAL_H_
