/*
 * PressureGradientEngine.h
 *
 * R05 — Two-phase pressure-gradient correlation declarations.
 *
 * This header re-exports the pressure-gradient functions declared in
 * SisProd.h for modular compilation. Include this when you need only
 * this API without the full SisProd2 architecture.
 *
 * Functions/Classes declared here (defined in PressureGradientEngine.cpp):
 *   - computePressureGradient()
 *   - ProductionColumn::marchToWellheadPhysical()
 *
 * Migration reference: issues/sisprod-migration-plan.md, region R05.
 */

#ifndef MARLIM_PRESSURE_GRADIENT_ENGINE_H_
#define MARLIM_PRESSURE_GRADIENT_ENGINE_H_

#include "SisProd.h"

// All pressure-gradient functions are declared in SisProd.h within
// the marlim::sisprod2 namespace. This header exists for documentation
// and modular dependency tracking purposes.

#endif // MARLIM_PRESSURE_GRADIENT_ENGINE_H_