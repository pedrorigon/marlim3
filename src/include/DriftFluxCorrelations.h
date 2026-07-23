/*
 * DriftFluxCorrelations.h
 *
 * R03 — Two-phase drift-flux correlation declarations.
 *
 * This header re-exports the drift-flux correlation functions declared in
 * SisProd.h for modular compilation. Include this when you need only
 * the drift-flux correlation API without the full SisProd2 architecture.
 *
 * Functions declared here (defined in DriftFluxCorrelations.cpp):
 *   - bhagwatGhajar(), bhagwatGhajarMod()
 *   - choi(), hibikiIshii(), francaLahey()
 *   - c0UdDisperso(), c0UdAnularChurn(), c0UdEstratificado()
 *
 * Migration reference: issues/sisprod-migration-plan.md, region R03.
 */

#ifndef MARLIM_DRIFT_FLUX_CORRELATIONS_H_
#define MARLIM_DRIFT_FLUX_CORRELATIONS_H_

#include "SisProd.h"

// All drift-flux correlation functions are declared in SisProd.h within
// the marlim::sisprod2 namespace. This header exists for documentation
// and modular dependency tracking purposes.

#endif // MARLIM_DRIFT_FLUX_CORRELATIONS_H_