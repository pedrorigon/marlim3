#ifndef SISPRODSOURCESINTERNAL_H_
#define SISPRODSOURCESINTERNAL_H_

#include "SisProdSources.h"

namespace sisprod::sources {
namespace {

/// Reads the water and gas consumed by hydrate formation after 0.01 s for
/// envelope models 2 and 3, and updates BSW from the remaining free-water volume.
/// Inactive models leave both output arguments unchanged.
void consumeHydrateFormationMass(const SourceState &state, double &gas_consumido_Mg,
                                double &agua_consumida_Mw, int cellIndex);

/// Evaluates the donor fluid's in-situ water fraction at the receiving cell's
/// pressure and temperature, using the receiving fluid's gas-oil-ratio criterion.
double inSituWaterFraction(const Cel &cell, Cel &fluidCell);

/// Updates source terms for leaks, multiple sources and radial or 2D porous media.
/// Leaks accumulate mass; the other models assign their phase contributions.
/// A primary parallel-network connection may supply its recorded leak rates.
void refreshChokeMultipleAndPorousSources(const SourceState &state, int cellIndex);

}  // namespace
}  // namespace sisprod::sources

#endif  // SISPRODSOURCESINTERNAL_H_
