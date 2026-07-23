/*
 * SisProd_r03.cpp
 *
 * R03 — Two-phase drift-flux correlations (C0 / Ud).
 *
 * Pure functions ported verbatim from SProd methods. The only structural
 * change is that they are free functions in the marlim::sisprod2 namespace
 * instead of SProd member methods, and they do not access arq.* fields.
 * The dispatcher functions (c0UdDisperso, c0UdAnularChurn, c0UdEstratificado)
 * receive the mode as an explicit integer parameter.
 *
 * Migration reference: issues/sisprod-migration-plan.md, region R03.
 */

#include "SisProd.h"

#include <cmath>
#include <limits>

namespace marlim {
namespace sisprod2 {

// ===========================================================================
// Internal helpers (anonymous namespace — not exported).
// ===========================================================================
namespace {

/// Iterative Colebrook friction factor.
/// Turbulent (Re > 2400): Halland seed + one refinement step.
/// Laminar: 64/Re.
static double colebrookFriction(double re, double relRoughness) {
    if (re <= 2400.0)
        return 64.0 / (re + 1e-15);
    // Halland approximation as seed.
    const double halland = 1.0 / std::pow(
        -1.8 * std::log10(std::pow(relRoughness / 3.7, 1.11) +
                          6.9 / (re + 1e-15)), 2.0);
    // One Colebrook step (delta < 1e-3 after one iteration for all Re).
    const double den = -2.0 * std::log10(
        relRoughness / 3.7 +
        2.51 / ((re + 1e-15) * std::sqrt(std::fabs(halland))));
    return 1.0 / (den * den);
}

/// Sign-correction for very-low-velocity cases used by all C0/Ud correlations.
static void applyLowVelocitySign(double ug1, double ul1, double area,
                                  double tet, double &ud) {
    const double vmix = (std::fabs(ug1) + std::fabs(ul1)) / area;
    if (vmix < 0.01) {
        if (tet > 0.0 && ud < 0.0)       ud =  std::fabs(ud);
        else if (tet < 0.0 && ud > 0.0)  ud = -std::fabs(ud);
    }
}

} // anonymous namespace

// ===========================================================================
// bhagwatGhajar
// Drift-flux C0/Ud for dispersed/slug/churn regime, mixture Reynolds.
// Mirrors SProd::BhagwatGhajar verbatim.
// ===========================================================================
void bhagwatGhajar(double rhol, double rhog, double tensup, double alf,
                   double reymix, double /*reymixL*/, double ug1, double ul1,
                   double dia, double rug, double tet,
                   double &c0, double &ud, double correcHor) {
    using std::fabs; using std::sqrt; using std::pow;
    using std::sin;  using std::cos;  using std::log10;

    const double A1     = constants::kPi * dia * dia / 4.0;
    const double rhomix = alf * rhog + (1.0 - alf) * rhol;
    // Legacy always uses sinal = 1 regardless of sign of tet (preserved).
    const double sinal  = 1.0;
    const double re     = (reymix < 1e-7) ? 1e-7 : reymix;
    const double fat    = colebrookFriction(re, rug / dia);

    const double Froude = sqrt(rhog / (rhol - rhog + 1e-15)) *
                          (fabs(ug1) / A1) /
                          sqrt(9.81 * dia * sinal * cos(tet) + 1e-15);
    const double massQ  = (rhog * fabs(ug1)) /
                          (rhog * fabs(ug1) + rhol * fabs(ul1) + 1e-15);
    const double Beta   = fabs(ug1) / (fabs(ug1) + fabs(ul1) + 1e-15);
    const double rgrl2  = (rhog / rhol) * (rhog / rhol);
    const double re2    = (re / 1000.0) * (re / 1000.0);
    const double t1     = (2.0 - rgrl2) / (1.0 + re2);
    const double t2     = pow((1.0 + rgrl2 * sinal * cos(tet)) /
                              (1.0 + cos(tet) + 1e-15),
                              (1.0 - alf) / 5.0) /
                          (1.0 + 1.0 / (re2 + 1e-15));
    const double C1     = 0.2; // circular/annular duct (rectangular: 0.4)
    double C01 = (C1 - C1 * sqrt(rhog / rhol)) *
                 (pow(2.6 - Beta, 0.15) - sqrt(fat)) *
                 pow(1.0 - massQ, 1.5);
    if (ug1 * ul1 < 0.0)  C01 = 0.0;
    if (tet >= -50.0 * constants::kPi / 180.0 && tet <= 0.0 && Froude <= 0.1)
        C01 = 0.0;
    c0 = t1 + t2 + C01;

    const double viscl  = dia * (fabs(ug1 / A1) + fabs(ul1 / A1)) *
                          rhomix / (re + 1e-15);
    const double t_ud1  = 0.35 * sin(tet) + 0.45 * cos(tet) * sinal;
    const double t_ud2  = sqrt(9.81 * dia * (rhol - rhog) / (rhol + 1e-15)) *
                          sqrt(1.0 - alf);
    const double C2     = (viscl / 0.001 > 10.0)
                          ? pow(0.434 / (log10(viscl / 0.001) + 1e-15), 0.15)
                          : 1.0;
    const double La     = sqrt(tensup / (9.81 * (rhol - rhog + 1e-15))) / dia;
    const double C3     = (La < 0.025) ? pow(La / 0.025, 0.90) : 1.0;
    double C4 = 1.0;
    if (tet >= -50.0 * constants::kPi / 180.0 && tet < 0.0 && Froude <= 0.1)
        C4 = -1.0;
    ud = correcHor * t_ud1 * t_ud2 * C2 * C3 * C4;
    applyLowVelocitySign(ug1, ul1, A1, tet, ud);
}

// ===========================================================================
// bhagwatGhajarMod
// Same algorithm as bhagwatGhajar but driven by the liquid-only Reynolds.
// Mirrors SProd::BhagwatGhajarMod verbatim.
// ===========================================================================
void bhagwatGhajarMod(double rhol, double rhog, double tensup, double alf,
                      double reymix, double reymixL, double ug1, double ul1,
                      double dia, double rug, double tet,
                      double &c0, double &ud, double correcHor) {
    bhagwatGhajar(rhol, rhog, tensup, alf, reymixL, reymixL,
                  ug1, ul1, dia, rug, tet, c0, ud, correcHor);
    (void)reymix;
}

// ===========================================================================
// choi
// Simplified drift-flux correlation. Mirrors SProd::Choi verbatim.
// ===========================================================================
void choi(double rhol, double rhog, double tensup, double alf,
          double reymix, double /*reymixL*/, double ug1, double ul1,
          double dia, double /*rug*/, double tet,
          double &c0, double &ud, double correcHor) {
    using std::fabs; using std::sqrt; using std::pow;
    using std::sin;  using std::cos;  using std::exp;

    const double A1    = constants::kPi * dia * dia / 4.0;
    const double sinal = (tet < 0.0) ? -1.0 : 1.0;
    ud = correcHor * sinal * 0.0246 * cos(tet) +
         1.606 * pow(9.82 * tensup * (rhol - rhog) /
                     (rhol * rhol + 1e-15), 0.25) * sin(tet);
    c0 = 2.0 / (1.0 + pow(reymix / 1000.0, 2.0)) +
         (1.2 - 0.2 * sqrt(rhog / (rhol + 1e-15)) *
          (1.0 - exp(-18.0 * alf))) /
         (1.0 + pow(1000.0 / (reymix + 1e-15), 2.0));
    applyLowVelocitySign(ug1, ul1, A1, tet, ud);
}

// ===========================================================================
// hibikiIshii
// Annular-churn drift-flux. Mirrors SProd::HibikiIshii verbatim.
// ===========================================================================
void hibikiIshii(double rhol, double rhog, double /*tensup*/, double alf,
                 double /*reymix*/, double /*reymixL*/, double ug1, double ul1,
                 double dia, double /*rug*/, double tet,
                 double &c0, double &ud, double correcHor) {
    using std::fabs; using std::sqrt; using std::sin;

    const double A1    = constants::kPi * dia * dia / 4.0;
    const double sinal = (tet < 0.0) ? -1.0 : 1.0;
    const double denom = alf + 4.0 * sqrt(rhog / (rhol + 1e-15)) + 1e-15;
    c0 = 1.0 + (1.0 - alf) / denom;
    ud = (correcHor * sinal * (1.0 - alf) / denom) *
         sqrt(9.82 * fabs(sin(tet)) * dia * (rhol - rhog) * (1.0 - alf) /
              (0.015 * rhol + 1e-15));
    applyLowVelocitySign(ug1, ul1, A1, tet, ud);
}

// ===========================================================================
// francaLahey
// Constant-parameter stratified correlation. Mirrors SProd::FrancaLahey.
// ===========================================================================
void francaLahey(double /*rhol*/, double /*rhog*/, double /*tensup*/,
                 double /*alf*/,  double /*reymix*/, double /*reymixL*/,
                 double ug1,      double ul1,
                 double dia,      double /*rug*/, double tet,
                 double &c0, double &ud, double correcHor) {
    const double A1    = constants::kPi * dia * dia / 4.0;
    const double sinal = (tet < 0.0) ? -1.0 : 1.0;
    c0 = 1.04;
    ud = correcHor * sinal * 0.466;
    applyLowVelocitySign(ug1, ul1, A1, tet, ud);
}

// ===========================================================================
// Dispatcher helper (anonymous, not exported).
// mode: 0=Choi  1=BhagwatGhajar  2=FrancaLahey  3=HibikiIshii
//       4=BhagwatGhajarMod  5=angle-blended BhagwatGhajarMod/Choi
// ===========================================================================
namespace {

static void dispatchC0Ud(int mode,
                         double rhol, double rhog, double tensup, double alf,
                         double reymix, double reymixL, double ug1, double ul1,
                         double dia, double rug, double tet,
                         double &c0, double &ud, double correcHor) {
    using std::fabs;
    switch (mode) {
    case 1:
        bhagwatGhajar(rhol, rhog, tensup, alf, reymix, reymixL,
                      ug1, ul1, dia, rug, tet, c0, ud, correcHor);
        break;
    case 2:
        francaLahey(rhol, rhog, tensup, alf, reymix, reymixL,
                    ug1, ul1, dia, rug, tet, c0, ud, correcHor);
        break;
    case 3:
        hibikiIshii(rhol, rhog, tensup, alf, reymix, reymixL,
                    ug1, ul1, dia, rug, tet, c0, ud, correcHor);
        break;
    case 4:
        bhagwatGhajarMod(rhol, rhog, tensup, alf, reymix, reymixL,
                         ug1, ul1, dia, rug, tet, c0, ud, correcHor);
        break;
    case 5: {
        const double kPi = constants::kPi;
        if (fabs(tet) < 5.0 * kPi / 180.0) {
            bhagwatGhajar(rhol, rhog, tensup, alf, reymix, reymixL,
                          ug1, ul1, dia, rug, tet, c0, ud, correcHor);
        } else if (fabs(tet) > 20.0 * kPi / 180.0) {
            choi(rhol, rhog, tensup, alf, reymix, reymixL,
                 ug1, ul1, dia, rug, tet, c0, ud, correcHor);
        } else {
            double c0m = 0.0, udm = 0.0;
            bhagwatGhajarMod(rhol, rhog, tensup, alf, reymix, reymixL,
                             ug1, ul1, dia, rug, tet, c0m, udm, correcHor);
            choi(rhol, rhog, tensup, alf, reymix, reymixL,
                 ug1, ul1, dia, rug, tet, c0, ud, correcHor);
            const double raz =
                (fabs(tet) - 5.0 * kPi / 180.0) / (15.0 * kPi / 180.0);
            c0 = raz * c0 + (1.0 - raz) * c0m;
            ud = raz * ud + (1.0 - raz) * udm;
        }
        break;
    }
    default: // 0 — Choi
        choi(rhol, rhog, tensup, alf, reymix, reymixL,
             ug1, ul1, dia, rug, tet, c0, ud, correcHor);
        break;
    }
}

} // anonymous namespace

// ===========================================================================
// Public dispatcher functions (replace arq.CorreDisper/Anular/Estrat switch).
// ===========================================================================

void c0UdDisperso(double rhol, double rhog, double tensup, double alf,
                  double reymix, double reymixL, double ug1, double ul1,
                  double dia, double rug, double tet,
                  double &c0, double &ud, double correcHor, int mode) {
    dispatchC0Ud(mode, rhol, rhog, tensup, alf, reymix, reymixL,
                 ug1, ul1, dia, rug, tet, c0, ud, correcHor);
}

void c0UdAnularChurn(double rhol, double rhog, double tensup, double alf,
                     double reymix, double reymixL, double ug1, double ul1,
                     double dia, double rug, double tet,
                     double &c0, double &ud, double correcHor, int mode) {
    dispatchC0Ud(mode, rhol, rhog, tensup, alf, reymix, reymixL,
                 ug1, ul1, dia, rug, tet, c0, ud, correcHor);
}

void c0UdEstratificado(double rhol, double rhog, double tensup, double alf,
                       double reymix, double reymixL, double ug1, double ul1,
                       double dia, double rug, double tet,
                       double &c0, double &ud, double correcHor, int mode) {
    dispatchC0Ud(mode, rhol, rhog, tensup, alf, reymix, reymixL,
                 ug1, ul1, dia, rug, tet, c0, ud, correcHor);
}

} // namespace sisprod2
} // namespace marlim
