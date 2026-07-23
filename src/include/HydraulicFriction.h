/*
 * HydraulicFriction.h
 *
 * R06/R10 — Hydraulic friction and GLV equipment declarations.
 *
 * This header re-exports the hydraulic friction and GLV functions declared in
 * SisProd.h for modular compilation. Include this when you need only
 * this API without the full SisProd2 architecture.
 *
 * Functions declared here (defined in HydraulicFriction.cpp):
 *   - colebrookFrictionFactor(), colebrookFrictionSeed()
 *   - areaValvCali()
 *
 * Migration reference: issues/sisprod-migration-plan.md, regions R06/R10.
 */

#ifndef MARLIM_HYDRAULIC_FRICTION_H_
#define MARLIM_HYDRAULIC_FRICTION_H_

#include "SisProd.h"

// All hydraulic friction and GLV functions are declared in SisProd.h within
// the marlim::sisprod2 namespace. This header exists for documentation
// and modular dependency tracking purposes.

#endif // MARLIM_HYDRAULIC_FRICTION_H_