/*
 * BlackOilProperties.h
 *
 * R04 — Black-oil fluid-property correlation declarations.
 *
 * This header re-exports the black-oil property functions declared in
 * SisProd.h for modular compilation. Include this when you need only
 * the black-oil property API without the full SisProd2 architecture.
 *
 * Functions declared here (defined in BlackOilProperties.cpp):
 *   - zFactor(), gasDensityBlackOil(), gasViscosityBlackOil()
 *   - solutionGOR(), oilFVF(), waterFVFBlackOil()
 *   - waterDensityBlackOil(), oilDensityBlackOil(), liquidDensityBlackOil()
 *   - waterViscosityBlackOil(), deadOilViscosity*(), oilViscosityBlackOil()
 *   - liquidViscosityBlackOil(), liquidSpecificHeatBlackOil()
 *   - gasSpecificHeatBlackOil(), liquidDensityDerivativeTBlackOil()
 *   - liquidThermalConductivityBlackOil(), gasThermalConductivityBlackOil()
 *   - liquidJouleThomsonBlackOil(), gasJouleThomsonBlackOil()
 *   - liquidEnthalpyBlackOil(), gasEnthalpyBlackOil()
 *
 * Migration reference: issues/sisprod-migration-plan.md, region R04.
 */

#ifndef MARLIM_BLACK_OIL_PROPERTIES_H_
#define MARLIM_BLACK_OIL_PROPERTIES_H_

#include "SisProd.h"

// All black-oil property functions are declared in SisProd.h within
// the marlim::sisprod2 namespace. This header exists for documentation
// and modular dependency tracking purposes.

#endif // MARLIM_BLACK_OIL_PROPERTIES_H_