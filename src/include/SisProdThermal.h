#ifndef SISPRODTHERMAL_H_
#define SISPRODTHERMAL_H_

#include <vector>

class Cel;
class Ler;
class choke;
class solverP3D;
struct varGlob1D;

namespace sisprod::thermal {

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

}  // namespace sisprod::thermal

#endif  // SISPRODTHERMAL_H_
