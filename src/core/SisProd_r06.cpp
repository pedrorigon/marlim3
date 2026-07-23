/*
 * SisProd_r06.cpp
 *
 * R06/R10 — Hydraulic friction and gas-lift valve kernels.
 *
 * Pure functions ported verbatim from SProd methods.  They carry no global
 * state and do not access celula[], arq.* or any SProd field.
 *
 * colebrookFrictionFactor
 *   Full Colebrook–White iterative algorithm (Halland seed + two refinement
 *   steps).  Identical to the per-function static helper in SisProd_r03.cpp
 *   but exported so the pressure-march and thermal kernels can share it.
 *
 * areaValvCali  (R10 — control equipment)
 *   GLV opening fraction.  Mirrors SProd::areaValvCali verbatim; every
 *   constant and formula is preserved.  The only change is that it is a free
 *   function instead of a member.
 *
 * Migration reference: issues/sisprod-migration-plan.md, regions R06 / R10.
 */

#include "SisProd2.h"

#include <cmath>

namespace marlim {
namespace sisprod2 {

// ===========================================================================
// colebrookFrictionFactor
// Darcy friction factor, Colebrook–White implicit equation solved with the
// Halland approximation as seed and two Newton-style refinement steps.
// Laminar (Re ≤ 2400): 64/Re.
// The same algorithm is used as a static helper in SisProd_r03.cpp (for
// BhagwatGhajar) and matches the Colebrook loop in the legacy SProd.
// ===========================================================================
double colebrookFrictionFactor(double reynolds, double relRoughness) {
    if (reynolds <= 0.0) return 64.0; // degenerate: return laminar limit

    if (reynolds <= 2400.0)
        return 64.0 / (reynolds + 1e-15);

    // Halland approximation — good to ~1 % absolute for any Re and ε/D.
    const double eps  = relRoughness;
    const double logArg = std::pow(eps / 3.7, 1.11) + 6.9 / (reynolds + 1e-15);
    double f = 1.0 / std::pow(-1.8 * std::log10(logArg + 1e-15), 2.0);
    if (f <= 0.0) f = 0.01; // safety floor

    // Two Colebrook refinement steps (converges to within ~1e-5 for all Re).
    for (int k = 0; k < 2; ++k) {
        const double sqrtF = std::sqrt(std::fabs(f));
        const double den   = -2.0 * std::log10(eps / 3.7 +
                                                2.51 / ((reynolds + 1e-15) *
                                                         sqrtF));
        if (std::fabs(den) < 1e-15) break;
        f = 1.0 / (den * den);
    }
    return f;
}

// ===========================================================================
// areaValvCali
// Gas-lift GLV opening fraction.
// Mirrors SProd::areaValvCali verbatim (see SisProd_old.cpp L2476).
// All physics constants and the formula structure are unchanged.
// Units: pressures in psi (legacy), diameter in metres, temperature in °F.
// ===========================================================================
double areaValvCali(double PCal, double TCal, double PVO, double PT,
                    double dextern, double areagarg, double Rvalv, double Temp) {
    using std::sqrt; using std::pow; using std::fabs;

    // Pressure of the bellows at 80 °F, corrected to the actual temperature.
    double PB80 = PCal * (1.0 - Rvalv);
    PB80 = (PB80 + 14.6959488) * (80.0 + 460.67) / (TCal * 1.8 + 491.67)
           - 14.6959488;
    const double PBT = PB80 * (1.0 + 0.00215 * (Temp - 80.0));

    // Effective port area (folha) and valve stroke.
    // XMVS > 0 iff the pressure difference forces the stem open; if XMVS ≤ 0
    // the valve is mechanically closed regardless of the static opening test.
    const double areafol = areagarg / Rvalv;
    const double BSR = (dextern * 100.0 / 2.54 > 1.1)
                       ? 500.0 * areafol
                       : 1950.0 * areafol;
    const double XMVS = ((PVO - PBT) * areafol - (PVO - PT) * areagarg) / BSR;

    // If the valve is mechanically closed, return zero.
    if (XMVS <= 0.0)
        return 0.0;

    // Effective port area from geometry (annular opening around the stem).
    const double DP  = sqrt(areagarg * 4.0 / constants::kPi);
    const double RP  = DP * 0.5;
    const double RB  = sqrt(areafol / constants::kPi);
    const double RB2mRP2 = RB * RB - RP * RP;
    const double sqrtRB2mRP2 = (RB2mRP2 > 0.0) ? sqrt(RB2mRP2) : 0.0;
    double APE = constants::kPi * RP * XMVS *
                 (XMVS + 2.0 * sqrtRB2mRP2);
    const double hyp = XMVS + sqrtRB2mRP2;
    const double denom = sqrt(hyp * hyp + RP * RP);
    APE /= (denom > 1e-15 ? denom : 1e-15);
    if (APE > areagarg) APE = areagarg;

    return APE / areagarg;
}

} // namespace sisprod2
} // namespace marlim
