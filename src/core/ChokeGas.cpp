/*
 * ChokeGas.cpp
 *
 * R10 Phase 1 — Gas/liquid flow through choke and Venturi restrictions.
 *
 * Ported from src/core/chokegas.cpp into marlim::sisprod2.
 *
 * Key algorithms:
 *   - calculateMassFlow(): Main mass-flow calculation with 3 code paths:
 *       fluido == 0 → gas (tipo 0/1 → simple choked, tipo ≥ 2 → Venturi table)
 *       fluido != 0 → liquid
 *   - criticalPressureRatio(): Isentropic critical pressure ratio
 *   - calcPressureRatio(): Iterative pressure-ratio solver (Newton)
 *   - throatTemperature(): Isentropic throat temperature
 *   - fluidDensity(): Brine/saline water density
 *
 * Migration reference: issues/sisprod-migration-plan.md, R10 Phase 1.
 */

#include "ChokeGas.h"
#include "BlackOilProperties.h"
#include <cmath>
#include <cstdlib>

namespace marlim {
namespace sisprod2 {

// =========================================================================
// Public API
// =========================================================================

double ChokeGas::calculateMassFlow(int fluido, double salin) {
    double massica = 0.0;
    double sens = 1.0;

    if (presEstag < presGarg) {
        sens = 0.0;
    }

    if (fluido == 0 && sens > 0.) {
        // --- GAS PATH ---
        double kad = calcAdiabaticIndex();
        double rho0 = gasDensity(presEstag, tempEstag);   // upstream density
        double rt = presEstag * 98066.5 / (rho0 * calcZFactor());

        if (tipo == 0 || tipo == 1) {
            // Simple choked/nozzle flow (tipo 0 or 1)
            double razcrit = pow((kad + 1) / 2., -kad / (kad - 1.));

            if (presGarg / presEstag > razcrit) {
                // Sub-critical flow
                massica = sqrt(2. * kad / (kad - 1.) *
                               pow(presGarg / presEstag, 2. / kad) *
                               (1. - pow(presGarg / presEstag, (kad - 1.) / kad)));
                massica *= (areagarg * presEstag * 98066.5 / sqrt(rt));
                tempGarg = pow(presGarg / presEstag, (kad - 1) / kad) * tempEstag;
            } else {
                // Critical (sonic) flow
                double max = sqrt(2. * kad / (kad - 1.) *
                                  pow(razcrit, 2. / kad) *
                                  (1. - pow(razcrit, (kad - 1.) / kad)));
                max *= (areagarg * presEstag * 98066.5 / sqrt(rt));
                tempGarg = pow(razcrit, (kad - 1) / kad) * tempEstag;
                massica = max;
            }
        } else {
            // Venturi empirical correlation (tipo >= 2)
            double SGG = flui.Deng;
            double PBAR = 0.980665 * presEstag;
            double TC = tempEstag;
            double agmm = areagarg * (1000.0 * 1000.0);

            // Bin SGG → sgg_group index [0, 1, 2, 3, 4]
            int JDG = static_cast<int>(std::floor(10.0 * (SGG - 0.55)));
            if (JDG < 0) JDG = 0;
            if (JDG > 3) JDG = 3;

            // Bin pressure (bar) → pressure_group index [0, ..., 8]
            int JPR = static_cast<int>(std::floor(0.02 * PBAR));
            if (JPR > 7) JPR = 7;

            // Bin temperature (°C) → temperature_group index [0, ..., 8]
            int JTE = static_cast<int>(std::floor(0.05 * TC)) - 1;
            if (JTE < 0) JTE = 0;
            if (JTE > 8) JTE = 8;

            // 3D bilinear interpolation in ventCR table
            // Interpolate along temperature at [JDG][JPR] and [JDG+1][JPR]
            double CRPant = (0.05 * TC - 1 - JTE) *
                                (kVentCR[JDG][JPR][JTE + 1] - kVentCR[JDG][JPR][JTE]) +
                            kVentCR[JDG][JPR][JTE];
            double CRPpos = (0.05 * TC - 1 - JTE) *
                                (kVentCR[JDG][JPR + 1][JTE + 1] - kVentCR[JDG][JPR + 1][JTE]) +
                            kVentCR[JDG][JPR + 1][JTE];
            double CRDant = (0.02 * PBAR - JPR) * (CRPpos - CRPant) + CRPant;

            CRPant = (0.05 * TC - 1 - JTE) *
                         (kVentCR[JDG + 1][JPR][JTE + 1] - kVentCR[JDG + 1][JPR][JTE]) +
                     kVentCR[JDG + 1][JPR][JTE];
            CRPpos = (0.05 * TC - 1 - JTE) *
                         (kVentCR[JDG + 1][JPR + 1][JTE + 1] - kVentCR[JDG + 1][JPR + 1][JTE]) +
                     kVentCR[JDG + 1][JPR + 1][JTE];
            double CRDpos = (0.02 * PBAR - JPR) * (CRPpos - CRPant) + CRPant;

            double CR;
            if (JDG == 1) {
                CR = ((SGG - 0.5538) / (0.65 - 0.5538)) * (CRDpos - CRDant) + CRDant;
            } else {
                CR = (10.0 * (SGG - 0.55) - JDG) * (CRDpos - CRDant) + CRDant;
            }

            double Rhostd = gasDensity(1.03322745, 20.0);
            double Rg = (8.0465 * 1000.0 * 1e5) / ((flui.Deng * 28.9625) * 100000.0);

            double QGcrt = (8640.0 * agmm * CR * PBAR / Rhostd) /
                           sqrt(Rg * (TC + 273.15));

            double PSIC = 0.97;
            if (PSIC > 0.97) PSIC = 0.97;
            if (PSIC < 0.54) PSIC = 0.54;

            double P3P1 = presGarg / presEstag;
            int ICRIT = 0;
            if (P3P1 < PSIC) ICRIT = 1;

            double QG;
            if (ICRIT == 1) {
                // Choked flow
                QG = QGcrt;
            } else {
                // Sub-critical empirical formula
                QG = QGcrt * (1.0 + pow(P3P1 - PSIC, 3)) *
                     (1.0 - pow((P3P1 - PSIC) / (1.0 - PSIC), 2.5));
            }

            if (P3P1 < PSIC) {
                tempGarg = pow(PSIC, (kad - 1) / kad) * tempEstag;
            } else {
                tempGarg = pow(P3P1, (kad - 1) / kad) * tempEstag;
            }

            massica = QG * Rhostd / 86400.0;
        }

        double rho1 = gasDensity(presGarg, tempGarg);
        qGarg = cd * sens * massica / rho1;
        return cd * sens * massica;
    } else {
        // --- LIQUID PATH ---
        double rho0 = fluidDensity(presEstag, tempEstag, salin);
        if (presEstag > presGarg) {
            massica = areagarg * sqrt(2. * rho0 * (presEstag - presGarg) * 98066.5);
        } else {
            massica = 0.0;
        }
        tempGarg = tempEstag;
        double rho1 = fluidDensity(presGarg, tempGarg, salin);
        qGarg = cd * sens * massica / rho1;
        return cdliq * sens * massica;
    }
}

double ChokeGas::throatTemperature() {
    double kad = calcAdiabaticIndex();
    double raz = pow(presEstag / presGarg, (kad - 1) / kad);
    return (tempEstag + 273.1) / raz - 273.1;
}

double ChokeGas::criticalPressureRatio() const {
    double kad = calcAdiabaticIndex();
    return pow(1.0 + 0.5 * (kad - 1.0), -kad / (kad - 1.0));
}

double ChokeGas::fluidDensity(double pres, double temper, double salin) const {
    double tfarAmb = faren(20.0);
    double ppsiAmb = psia(1.0);
    double tfar = faren(temper);
    double ppsi = psia(pres);

    double bwAmb = 1.0 + 1.2e-04 * (tfarAmb - 60.0) +
                   1.0e-06 * pow(tfarAmb - 60.0, 2.0) -
                   3.33e-06 * ppsiAmb;
    double bw = 1.0 + 1.2e-04 * (tfar - 60.0) +
                1.0e-06 * pow(tfar - 60.0, 2.0) -
                3.33e-06 * ppsi;

    double rholw = (1000.0 / bw) * bwAmb;
    double rhosal = 2160.0;
    double x = salin / 1000.0;
    double rhomist = (1.0 - x) / rholw + x / rhosal;
    return 1.0 / rhomist;
}

// =========================================================================
// Internal pressure-ratio helpers
// =========================================================================

double ChokeGas::calcPressureRatio(double mass, double rp) const {
    double kad = calcAdiabaticIndex();
    double rho0 = gasDensity(presEstag, tempEstag);
    double rt = presEstag * 98066.5 / (rho0 * calcZFactor());
    double rpcrit = criticalPressureRatio();
    rp = razpresSimples(mass, rp);
    double aux = newtonRoot(kad, rt, mass, pow(rp, 1.0 / kad));
    aux = pow(std::fabs(aux), kad);
    rp = aux;
    if (rp > 1.) rp = 1.;
    if (rp < rpcrit) rp = rpcrit;
    return rp;
}

// =========================================================================
// Internal helpers
// =========================================================================

int ChokeGas::sggBin(double sgg) const {
    // Maps specific-gravity bins:  0.55, 0.65, 0.75, 0.85, 0.95
    return static_cast<int>(std::floor(10.0 * (sgg - 0.55)));
}

double ChokeGas::interpolateVentCR(double sgg, double pbar, double tc) const {
    // Bin indices per legacy chokegas.cpp
    int JDG = static_cast<int>(std::floor(10.0 * (sgg - 0.55)));
    if (JDG < 0) JDG = 0;
    if (JDG > 3) JDG = 3;

    int JPR = static_cast<int>(std::floor(0.02 * pbar));
    if (JPR > 7) JPR = 7;

    int JTE = static_cast<int>(std::floor(0.05 * tc)) - 1;
    if (JTE < 0) JTE = 0;
    if (JTE > 8) JTE = 8;

    return kVentCR[JDG][JPR][JTE];
}

double ChokeGas::calcAdiabaticIndex() const {
    return flui.ConstAdG(presEstag, tempEstag);
}

// --- Numerics (ported verbatim) ---

double ChokeGas::fraiz(double kad, double rt, double mass, double rp) const {
    return pow((mass * sqrt(rt) / (areagarg * presEstag * 98066.5)), 2.) -
           (2. * kad / (kad - 1.)) * rp * rp * (1. - pow(rp, kad - 1.));
}

double ChokeGas::derraiz(double kad, double rp) const {
    return -(2. * kad / (kad - 1.)) *
           (2. * rp * (1. - pow(rp, kad - 1.)) -
            (kad - 1.) * rp * rp * pow(rp, kad - 2.));
}

double ChokeGas::newtonRoot(double kad, double rt, double mass, double rp) const {
    double errox = 100.;
    double errof = fraiz(kad, rt, mass, rp);

    int konta = 0;
    while ((std::fabs(errox) > 1e-5 || std::fabs(errof) > 1e-5) && konta < 100) {
        double deri = derraiz(kad, rp);
        if (std::fabs(deri) < 1e-30) break;    // safeguard against division by zero
        errox = -errof / deri;
        rp = rp + errox;
        errof = fraiz(kad, rt, mass, rp);
        konta++;
    }
    return rp;
}

double ChokeGas::razpresSimples(double mass, double rp) const {
    double kad = calcAdiabaticIndex();
    double rho0 = gasDensity(presEstag, tempEstag);
    double rhoG = gasDensity(presGarg, tempEstag);
    double rt = presEstag * 98066.5 / (rho0 * calcZFactor());
    double cd = 0.885;
    double vgarg = mass / (rhoG * areagarg * cd);
    double vsom = sqrt((2. * kad / (kad + 1.)) * rt);
    double Mach = vgarg / vsom;
    double rpcrit = criticalPressureRatio();
    double aux = pow(1. + 0.5 * (kad - 1.) * Mach * Mach, kad / (kad - 1.));
    double rpS = 1. / aux;
    if (rpS > 1.) rpS = 1.;
    if (rpS < rpcrit) rpS = rpcrit;
    return rpS;
}

// Venturi Mach number numerics
double ChokeGas::fMachVenturi(double Mach, double kad) const {
    double razarea = areafole / areagarg;
    return razarea -
           pow((1. + 0.5 * (kad - 1.) * Mach * Mach) / (0.5 * (kad + 1.)),
               0.5 * (kad + 1.) / (kad - 1.)) /
           Mach;
}

double ChokeGas::deriFMachVenturi(double Mach, double kad) const {
    double expre1 = (1. + 0.5 * (kad - 1.) * Mach * Mach) / (0.5 * (kad + 1.));
    return (1. / (Mach * Mach)) * pow(expre1, 0.5 * (kad + 1.) / (kad - 1.)) -
           pow(expre1, (-0.5 * kad + 1.5) / (kad - 1.));
}

double ChokeGas::findMachVenturi(double kad) const {
    double errox = 100.;
    double MachChute = 0.05;
    double errof = fMachVenturi(MachChute, kad);

    int konta = 0;
    while ((std::fabs(errox) > 1e-5 || std::fabs(errof) > 1e-5) && konta < 100) {
        double deri = deriFMachVenturi(MachChute, kad);
        if (std::fabs(deri) < 1e-30) break;
        errox = -errof / deri;
        MachChute = MachChute + errox;
        errof = fMachVenturi(MachChute, kad);
        konta++;
    }
    return MachChute;
}

double ChokeGas::criticalPressureRatioVenturi() const {
    double kad = calcAdiabaticIndex();
    double Ma = findMachVenturi(kad);
    return 1. / pow(1. + 0.5 * (kad - 1.) * Ma * Ma, kad / (kad - 1.));
}

}  // namespace sisprod2
}  // namespace marlim
