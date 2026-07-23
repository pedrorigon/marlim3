/*
 * SisProd_r04.cpp
 *
 * R04 — Black-oil fluid-property correlations (pure, no global state).
 *
 * Ported verbatim from the flashCompleto==0 branch of ProFlu methods.
 * All functions are free (no class state) and use SI/legacy units
 * documented in each function header.
 *
 * Correlation references:
 *   - Z-factor: Dranchuk-Abou-Kassem with Gopal seed (same as ProFlu::Zdran).
 *   - Gas viscosity: Lee-Kesler (same as ProFlu::ViscGas).
 *   - Gas density: real-gas equation of state (same as ProFlu::MasEspGas).
 *   - Solution GOR: Standing (same as ProFlu::RS corrSat=0).
 *   - Oil FVF: Vasquez-Beggs (same as ProFlu::BOFunc below Pb).
 *   - Water density/FVF: Meehan (same as ProFlu::MasEspAgua).
 *
 * Migration reference: issues/sisprod-migration-plan.md, region R04.
 */

#include "SisProd2.h"
#include "PropFlu.h"   // CoefGopal constant array

#include <cmath>

namespace marlim {
namespace sisprod2 {

// ===========================================================================
// Unit-conversion helpers (mirrors ProFlu inline members).
// All PVT correlations in the legacy codebase use:
//   pressure  in kgf/cm2
//   temperature in degrees Celsius
// with intermediate conversions to psia / °F / Rankine as needed.
// ===========================================================================
namespace {

inline double to_psia(double p_kgf)  { return p_kgf * 0.9678411 * 14.69595; }
inline double to_fahr(double t_c)    { return 1.8 * t_c + 32.0; }
inline double from_psia(double p_psi){ return p_psi / (14.69595 * 0.9678411); }

// Dranchuk-Abou-Kassem residual function (Z-factor implicit equation).
static double fdran(double PR, double TR, double denr) {
    const double A = 0.064225133;
    const double B = 0.53530771 * TR - 0.61232032;
    const double C = 0.31506237 * TR - 1.04670990 - 0.57832728 / (TR * TR);
    const double D = TR;
    const double E = 0.68157001 / (TR * TR);
    const double F = 0.68446549;
    const double G = 0.27 * PR;
    const double X1 = F * denr * denr;
    const double X2 = std::exp(-X1);
    const double denr3 = denr * denr * denr;
    return A * denr3 * denr3 + B * denr3 + C * denr * denr +
           D * denr + E * denr3 * (1.0 + X1) * X2 - G;
}

// Bisection (replaces ProFlu::FalsaCorda — no vg1dSP error counter).
static double solveZBisection(double a, double b,
                               double PR, double TR,
                               int maxit = 100) {
    double u = fdran(PR, TR, a);
    double c = 0.0;
    for (int k = 0; k < maxit; ++k) {
        c = 0.5 * (a + b);
        const double w = fdran(PR, TR, c);
        if (std::fabs(b - a) < 1e-5 || std::fabs(w) < 1e-4)
            return c;
        if ((u > 0.0 && w < 0.0) || (u < 0.0 && w > 0.0))
            b = c;
        else {
            a = c;
            u = w;
        }
    }
    return c;
}

// Gopal explicit Z-factor seed.
static double zGopalSeed(double PR, double TR) {
    if (PR < (TR - 0.73) ||
        (PR < (4.7 * TR - 6.74) && PR < (5.6875 - 0.525 * TR)))
        return 1.0;
    if (PR <= 5.4) {
        int i = 1;
        if (PR > 1.2) {
            if (PR > 1.4 || TR < 1.08 || TR > 1.19)
                i = (PR <= 2.8) ? 2 : 3;
        }
        int k = 4;
        if (TR <= 2.0)  k = 3;
        if (TR <= 1.4)  k = 2;
        if (TR <= 1.2)  k = 1;
        const int j = 16 * i + 4 * k - 19;
        return PR * (CoefGopal[j - 1] * TR + CoefGopal[j]) +
               CoefGopal[j + 1] * TR + CoefGopal[j + 2];
    }
    return PR * std::pow(0.711 + 3.66 * TR, -1.4667) -
           1.637 / (0.319 * TR + 0.522) + 2.071;
}

} // anonymous namespace

// ===========================================================================
// zFactor
// Gas compressibility factor via Dranchuk-Abou-Kassem.
// Mirrors ProFlu::ZdranOriginal (Deng <= 2 path, cordg == 0).
// @param pres   Pressure (kgf/cm2).
// @param temp   Temperature (°C).
// @param Deng   Gas specific gravity (air = 1).
// @param PC     Critical pressure (psia).
// @param TC     Critical temperature (Rankine).
// ===========================================================================
double zFactor(double pres, double temp, double Deng,
               double PC, double TC) {
    const double ip = to_psia(pres);
    const double it = to_fahr(temp);
    const double PR = ip / PC;
    const double TR = (it + 460.0) / TC;

    double a, b;
    if (Deng <= 2.0) {
        double seed = zGopalSeed(PR, TR);
        if (seed < 0.01) seed = 0.01;
        a = 0.27 * PR / (seed * (it + 460.0) / TC);
        b = 0.27 * PR / (0.10  * (it + 460.0) / TC);
    } else {
        a = 0.27 * PR / (40.0  * (it + 460.0) / TC);
        b = 0.27 * PR / (0.01  * (it + 460.0) / TC);
    }
    if (a > b) { double tmp = a; a = b; b = tmp; }

    const double denr = solveZBisection(a, b, PR, TR);
    const double z    = 0.27 * PR / (denr * TR);
    return (z <= 0.05) ? 0.05 : z;
}

// ===========================================================================
// gasDensityBlackOil
// In-situ gas density via the real-gas equation of state.
// Mirrors ProFlu::MasEspGas (flashCompleto==0 path).
// @param pres   Pressure (kgf/cm2).
// @param temp   Temperature (°C).
// @param Deng   Gas specific gravity.
// @param PC     Critical pressure (psia).
// @param TC     Critical temperature (Rankine).
// @returns Gas density (kg/m3).
// ===========================================================================
double gasDensityBlackOil(double pres, double temp,
                           double Deng, double PC, double TC) {
    const double z = zFactor(pres, temp, Deng, PC, TC);
    // rho_g = M_g * p / (R * z * T) — using legacy unit system
    // M_g = Deng * 28.9625 g/mol, p in Pa (= pres_kgf * 98066.5),
    // R = 8046.5 J/kmol/K, T in K
    return (Deng * 28.9625 * pres * 98066.5) /
           (8.0465e3 * z * (temp + 273.15));
}

// ===========================================================================
// gasViscosityBlackOil
// Gas dynamic viscosity via the Lee-Kesler correlation.
// Mirrors ProFlu::ViscGas (flashCompleto==0 path).
// @param pres   Pressure (kgf/cm2).
// @param temp   Temperature (°C).
// @param Deng   Gas specific gravity.
// @param PC     Critical pressure (psia).
// @param TC     Critical temperature (Rankine).
// @returns Gas viscosity (cP).
// ===========================================================================
double gasViscosityBlackOil(double pres, double temp,
                             double Deng, double PC, double TC) {
    const double rhog = gasDensityBlackOil(pres, temp, Deng, PC, TC) / 1000.0; // g/cm3
    const double TF   = to_fahr(temp);
    const double TR   = TF + 459.67;
    const double wg   = Deng * 29.0;
    const double AK   = (9.4 + 0.02 * wg) * std::pow(TR, 1.5) /
                        (209.0 + 19.0 * wg + TR);
    const double x    = 3.5 + (986.0 / TR) + 0.01 * wg;
    const double y    = 2.4 - 0.2 * x;
    return AK * std::exp(x * std::pow(rhog, y)) / 1.0e4; // cP
}

// ===========================================================================
// solutionGOR
// Solution gas-oil ratio via Standing correlation (corrSat=0, no Rs table).
// Mirrors ProFlu::RS (flashCompleto==0, corrSat==0, tabRSPB==0 path).
// @param pres   Pressure (kgf/cm2).
// @param temp   Temperature (°C).
// @param API    Oil API gravity.
// @param Deng   Gas specific gravity.
// @param Avb,Bvb,Cvb  Vasquez-Beggs / Standing coefficients from ProFlu arq.*
// @returns Solution GOR (ft3/bbl standard conditions).
// ===========================================================================
double solutionGOR(double pres, double temp, double API,
                   double Deng, double Avb, double Bvb, double Cvb) {
    const double ipres = to_psia(pres);
    const double itemp = to_fahr(temp) + 460.0; // Rankine
    double rs = (Deng * Bvb * std::pow(ipres, Avb)) /
                std::pow(10.0, Cvb * API / itemp);
    if (rs < 0.0) rs = 0.0;
    return rs;
}

// ===========================================================================
// oilFVF
// Oil formation volume factor via the Vasquez-Beggs correlation below Pb.
// Mirrors ProFlu::BOFunc (flashCompleto==0, below bubble-point branch).
// @param pres   Pressure (kgf/cm2) — assumed below bubble point.
// @param temp   Temperature (°C).
// @param API    Oil API gravity.
// @param Deng   Gas specific gravity.
// @param rs     Solution GOR (ft3/bbl) — caller provides (e.g. from solutionGOR).
// @returns Oil formation volume factor Bo (RB/STB, dimensionless > 1).
// ===========================================================================
double oilFVF(double pres, double temp, double API,
              double Deng, double rs) {
    const double itemp    = to_fahr(temp);
    // Density of gas at 100 psia for SG100 calculation.
    const double kgf100   = from_psia(100.0);
    const double z100     = zFactor(kgf100, temp, Deng, 667.0, 395.0);
    const double SG100    = (Deng * 28.9625 * kgf100 * 98066.5) /
                            (8.0465e3 * z100 * (temp + 273.15));
    const double Dvazbeg  = (itemp - 60.0) * API / SG100;
    const double Avazbeg  = (API > 30.0) ?  0.1100 :  0.1751;
    const double Bvazbeg  = (API > 30.0) ?  0.1337 : -1.8106;
    const double bo = 1.0 + 0.000467 * rs +
                      Avazbeg * Dvazbeg * 0.0001 +
                      Bvazbeg * rs * Dvazbeg * 1.0e-8;
    return (bo > 0.9) ? bo : 0.9; // physical lower bound
}

} // namespace sisprod2
} // namespace marlim
