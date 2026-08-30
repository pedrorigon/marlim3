#ifndef SISPRODTHERMAL_H_
#define SISPRODTHERMAL_H_

class Cel;
class Ler;

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

}  // namespace sisprod::thermal

#endif  // SISPRODTHERMAL_H_