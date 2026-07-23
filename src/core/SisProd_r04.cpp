/*
 * SisProd_r04.cpp
 *
 * R04 — Black-oil fluid-property correlations (pure, no global state).
 */

#include "SisProd2.h"

#include <cmath>

namespace marlim {
namespace sisprod2 {

namespace {

const double kCoefGopal[] = {
    1.6643, -2.2114, -0.3647, 1.4385, 0.5220, -0.8511,
    -0.0364, 1.0490, 0.1391, -0.2988, 0.0007, 0.9969,
    0.0295, -0.0825, 0.0009, 0.9967, -1.3570, 1.4942,
    4.6315, -4.7009, 0.1717, -0.3232, 0.5869, 0.1229,
    0.0984, -0.2053, 0.0621, 0.8580, 0.0211, -0.0527,
    0.0127, 0.9549, -0.3278, 0.4752, 1.8223, -1.9036,
    -0.2521, 0.3871, 1.6087, -1.6635, -0.0284, 0.0625,
    0.4714, -0.0011, 0.0041, 0.0039, 0.0607, 0.7927
};

inline double to_psia(double p_kgf) { return p_kgf * 0.9678411 * 14.69595; }
inline double to_fahr(double t_c) { return 1.8 * t_c + 32.0; }
inline double from_psia(double p_psi) { return p_psi / (14.69595 * 0.9678411); }

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

static double solveZBisection(double a, double b, double PR, double TR,
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
        if (TR <= 2.0) k = 3;
        if (TR <= 1.4) k = 2;
        if (TR <= 1.2) k = 1;
        const int j = 16 * i + 4 * k - 19;
        return PR * (kCoefGopal[j - 1] * TR + kCoefGopal[j]) +
               kCoefGopal[j + 1] * TR + kCoefGopal[j + 2];
    }
    return PR * std::pow(0.711 + 3.66 * TR, -1.4667) -
           1.637 / (0.319 * TR + 0.522) + 2.071;
}

} // namespace

double zFactor(double pres, double temp, double Deng, double PC, double TC) {
    const double ip = to_psia(pres);
    const double it = to_fahr(temp);
    const double PR = ip / PC;
    const double TR = (it + 460.0) / TC;

    double a, b;
    if (Deng <= 2.0) {
        double seed = zGopalSeed(PR, TR);
        if (seed < 0.01) seed = 0.01;
        a = 0.27 * PR / (seed * (it + 460.0) / TC);
        b = 0.27 * PR / (0.10 * (it + 460.0) / TC);
    } else {
        a = 0.27 * PR / (40.0 * (it + 460.0) / TC);
        b = 0.27 * PR / (0.01 * (it + 460.0) / TC);
    }
    if (a > b) {
        const double tmp = a;
        a = b;
        b = tmp;
    }

    const double denr = solveZBisection(a, b, PR, TR);
    const double z = 0.27 * PR / (denr * TR);
    return (z <= 0.05) ? 0.05 : z;
}

double gasDensityBlackOil(double pres, double temp,
                          double Deng, double PC, double TC) {
    const double z = zFactor(pres, temp, Deng, PC, TC);
    return (Deng * 28.9625 * pres * 98066.5) /
           (8.0465e3 * z * (temp + 273.15));
}

double gasViscosityBlackOil(double pres, double temp,
                            double Deng, double PC, double TC) {
    const double rhog = gasDensityBlackOil(pres, temp, Deng, PC, TC) / 1000.0;
    const double TF = to_fahr(temp);
    const double TR = TF + 459.67;
    const double wg = Deng * 29.0;
    const double AK = (9.4 + 0.02 * wg) * std::pow(TR, 1.5) /
                      (209.0 + 19.0 * wg + TR);
    const double x = 3.5 + (986.0 / TR) + 0.01 * wg;
    const double y = 2.4 - 0.2 * x;
    return AK * std::exp(x * std::pow(rhog, y)) / 1.0e4;
}

double solutionGOR(double pres, double temp, double API,
                   double Deng, double Avb, double Bvb, double Cvb) {
    const double ipres = to_psia(pres);
    const double itemp = to_fahr(temp) + 460.0;
    double rs = (Deng * Bvb * std::pow(ipres, Avb)) /
                std::pow(10.0, Cvb * API / itemp);
    if (rs < 0.0) rs = 0.0;
    return rs;
}

double oilFVF(double /*pres*/, double temp, double API,
              double Deng, double rs) {
    const double itemp = to_fahr(temp);
    const double kgf100 = from_psia(100.0);
    const double z100 = zFactor(kgf100, temp, Deng, 667.0, 395.0);
    const double SG100 = (Deng * 28.9625 * kgf100 * 98066.5) /
                         (8.0465e3 * z100 * (temp + 273.15));
    const double Dvazbeg = (itemp - 60.0) * API / SG100;
    const double Avazbeg = (API > 30.0) ? 0.1100 : 0.1751;
    const double Bvazbeg = (API > 30.0) ? 0.1337 : -1.8106;
    const double bo = 1.0 + 0.000467 * rs +
                      Avazbeg * Dvazbeg * 0.0001 +
                      Bvazbeg * rs * Dvazbeg * 1.0e-8;
    return (bo > 0.9) ? bo : 0.9;
}

double waterDensityBlackOil(double pres, double temp, double Denag) {
    const double tfarAmb = to_fahr(20.0);
    const double ppsiAmb = to_psia(1.0);
    const double tfar = to_fahr(temp);
    const double ppsi = to_psia(pres);
    const double bwAmb = 1.0 + 1.2e-4 * (tfarAmb - 60.0) +
                         1.0e-6 * std::pow(tfarAmb - 60.0, 2.0) -
                         3.33e-6 * ppsiAmb;
    const double bw = 1.0 + 1.2e-4 * (tfar - 60.0) +
                      1.0e-6 * std::pow(tfar - 60.0, 2.0) -
                      3.33e-6 * ppsi;
    const double rholw = (1000.0 / bw) * bwAmb;
    const double din = 2.16;
    const double vi = 1.0 - (Denag - din) / (1.0 - din);
    const double x = vi * din / Denag;
    const double rhomist = (1.0 - x) / rholw + x / (din * 1000.0);
    return 1.0 / rhomist;
}

double waterFVFBlackOil(double pres, double temp, double Denag) {
    return (1000.0 * Denag) / waterDensityBlackOil(pres, temp, Denag);
}

double oilDensityBlackOil(double pres, double temp, double API,
                          double Deng, double rs, double rDgD) {
    const double bo = oilFVF(pres, temp, API, Deng, rs);
    return (1000.0 * 141.5 / (131.5 + API) +
            rDgD * Deng * 1.225 * rs * 6.29 / 35.31467) / bo;
}

double liquidDensityBlackOil(double pres, double temp, double API,
                             double Deng, double BSW, double Denag,
                             double rs, double rDgD) {
    const double masoleo = (1.0 - BSW) *
        (1000.0 * 141.5 / (131.5 + API) +
         rDgD * Deng * 1.225 * rs * 6.29 / 35.31467);
    const double masagua = BSW * 1000.0 * Denag;
    double ba = 1.0;
    if (BSW > 0.0)
        ba = waterFVFBlackOil(pres, temp, Denag);
    return (masoleo + masagua) / ((1.0 - BSW) * oilFVF(pres, temp, API, Deng, rs) + BSW * ba);
}

double waterViscosityBlackOil(double temp) {
    return 2.414e-5 * std::pow(10.0, 247.8 / (temp + 133.15)) * 1000.0;
}

double deadOilViscosityBeggsRobinson(double temp, double API) {
    double tf = to_fahr(temp);
    if (tf < 50.0)
        tf = 50.0;
    const double sgo = 141.5 / (131.5 + API);
    const double y = std::exp(13.108 - 6.591 / sgo);
    const double x = y * std::pow(tf, -1.163);
    return std::pow(10.0, x) - 1.0;
}

double deadOilViscosityASTM(double temp, double API,
                            double TempL, double LVisL,
                            double TempH, double LVisH) {
    const double rhol = 141.5 / (131.5 + API);
    const double bASTM1 = std::log10(std::log10(LVisL / rhol + 0.7)) -
                          std::log10(std::log10(LVisH / rhol + 0.7));
    const double bASTM2 = std::log10(std::log10(LVisL / rhol + 0.7));
    const double bASTM = bASTM1 / std::log10((TempL + 273.0) / (TempH + 273.0));
    return rhol * (std::pow(10.0,
        std::pow(10.0, bASTM * std::log10((temp + 273.0) / (TempL + 273.0)) + bASTM2)) - 0.7);
}

double oilViscosityBlackOil(double rs, double deadOilViscosity) {
    if (rs < 0.0)
        rs = 0.0;
    const double xvisc = 10.715 / std::pow(rs + 100.0, 0.515);
    const double yvisc = 5.439 / std::pow(rs + 150.0, 0.338);
    return xvisc * std::pow(deadOilViscosity, yvisc);
}

double liquidViscosityBlackOil(double /*pres*/, double temp, double API,
                               double Deng, double BSW, double Denag,
                               double rs, double rDgD,
                               double TempL, double LVisL,
                               double TempH, double LVisH) {
    const double viso = oilViscosityBlackOil(
        rs, deadOilViscosityASTM(temp, API, TempL, LVisL, TempH, LVisH));
    const double masAgua = BSW * 1000.0 * Denag;
    const double masLiq = (1.0 - BSW) *
        (1000.0 * 141.5 / (131.5 + API) +
         rDgD * Deng * 1.225 * rs * 6.29 / 35.31467) + masAgua;
    const double xliq = masLiq > 0.0 ? masAgua / masLiq : 0.0;
    return (1.0 - xliq) * viso + waterViscosityBlackOil(temp) * xliq;
}

double liquidSpecificHeatBlackOil(double /*pres*/, double temp, double API,
                                  double Deng, double BSW, double Denag,
                                  double rs, double rDgD) {
    const double Bcp = 0.06103;
    const double tempfar = to_fahr(temp);
    const double tempK = temp + 273.16;
    const double masAgua = BSW * 1000.0 * Denag;
    const double masLiq = (1.0 - BSW) *
        (1000.0 * 141.5 / (131.5 + API) +
         rDgD * Deng * 1.225 * rs * 6.29 / 35.31467) + masAgua;
    const double xliq = masLiq > 0.0 ? masAgua / masLiq : 0.0;
    const double CPOI = 4187.0 * ((2.6948e-6 * API + 3.88402e-4) * tempfar +
                                  (0.0027665 * API + 0.366079) - Bcp);
    double CPWI;
    if (tempK < 410.0) {
        CPWI = 4185.5 * (2.13974 - 9.68137 * tempK / 1000.0 +
                         2.68536 * tempK * tempK / 100000.0 -
                         2.42139e-8 * tempK * tempK * tempK);
    } else {
        CPWI = 4185.5 * (-11.1558 + 7.96443 * tempK / 100.0 -
                         1.74799 * tempK * tempK / 10000.0 +
                         1.29156e-7 * tempK * tempK * tempK);
    }
    return xliq * CPWI + (1.0 - xliq) * CPOI;
}

double liquidThermalConductivityBlackOil(double pres, double temp,
                                         double API, double Deng,
                                         double BSW, double Denag,
                                         double rs, double rDgD) {
    const double tempfar = to_fahr(temp);
    const double tempK = temp + 273.16;
    const double XKOI = 116.8 * (1.0 - 3.0 * (tempfar - 32.0) / 10000.0) / 1000.0;

    const double masAgua = BSW * 1000.0 * Denag;
    const double masLiq = (1.0 - BSW) *
        (1000.0 * 141.5 / (131.5 + API) +
         rDgD * Deng * 1.225 * rs * 6.29 / 35.31467) + masAgua;
    const double xliq = masLiq > 0.0 ? masAgua / masLiq : 0.0;

    double XKWDI;
    if (tempK < 273.16)
        XKWDI = 418.4 * (273.778 + 3.9 * tempK) / 1000000.0;
    else if (tempK < 413.16)
        XKWDI = 418.4 * (-1390.53 + 15.1937 * tempK - 0.0190398 * tempK * tempK) / 1000000.0;
    else
        XKWDI = 418.4 * (-339.838 + 9.86669 * tempK - 0.0123045 * tempK * tempK) / 1000000.0;

    return XKOI * std::pow(XKWDI / XKOI, xliq);
    (void)pres;
}

double gasThermalConductivityBlackOil(double pres, double temp) {
    const double ppas = 98066.5 * pres;
    const double fatT = (temp + 273.15) / 191.1;
    const double XK1 = 3.04314 / 100.0 +
                       (1.3242 / 10000.0 + 1.27534 * temp / 10000000.0) * temp;
    const double RXK = 0.99783 +
                       (1.973e-8 + 7.8868e-16 * ppas) * ppas;
    if (fatT < 3.0)
        return XK1 * (1.0 + (1.0 - RXK) * (fatT - 3.0) / 1.354);
    return XK1;
}

double liquidJouleThomsonBlackOil(double pres, double temp, double API,
                                  double Deng, double BSW, double Denag,
                                  double rs, double rDgD,
                                  double liquidSimple) {
    const double bo = oilFVF(pres, temp, API, Deng, rs);
    const double roIS = oilDensityBlackOil(pres, temp, API, Deng, rs, rDgD);
    const double fw = BSW / (bo + BSW - BSW * bo);
    const double drholDt = liquidDensityDerivativeTBlackOil(pres, temp, API,
                                                            Deng, rs, rDgD);
    const double jtl = -(1.0 / roIS +
        (1.0 - liquidSimple) * (temp + 273.15) * (drholDt / (roIS * roIS)));
    const double jtw = -1.0 / (1000.0 * Denag)
                     + (temp + 273.15) * (3.0e-4) / (1000.0 * Denag);
    return (1.0 - fw) * jtl + fw * jtw;
}

double gasJouleThomsonBlackOil(double pres, double temp,
                               double Deng, double PCis, double TCis,
                               double rhog) {
    const double rho = rhog < 0.0
        ? gasDensityBlackOil(pres, temp, Deng, PCis, TCis)
        : rhog;
    const double z0 = zFactor(pres, temp, Deng, PCis, TCis);
    const double dzdt = (zFactor(pres, temp + 1e-3, Deng, PCis, TCis)
                       - zFactor(pres, temp - 1e-3, Deng, PCis, TCis)) / (2e-3);
    return (temp + 273.16) * dzdt / (z0 * rho);
}

 double gasSpecificHeatBlackOil(double pres, double temp, double Deng,
                                double PCis, double TCis,
                                double yco2, double rDgL) {
    const double RG = (8.0465 * 1000.0) / (rDgL * Deng * 28.9625);
    double tempK = temp + 273.16;
    if (tempK < 240.0) tempK = 240.0;
    if (tempK > 480.0) tempK = 480.0;

    const double CP0mM = 0.047801 * tempK + 21.714;
    const double CP0mE = 0.12532 * tempK + 15.341;
    const double CP0mP = 0.19347 * tempK + 16.371;
    const double CP0mB = 0.24799 * tempK + 25.170;
    const double CP0mCO2 = -3.7498e-5 * tempK * tempK + 6.6944e-2 * tempK + 20.538;

    double yM = (yco2 < 1.0) ? 1.7278095 + 0.26875568 * yco2 - 1.3142579 * rDgL * Deng : 0.0;
    if (yM > 1.0) yM = 1.0;
    if (yM < 0.0) yM = 0.0;
    double yE = 4.0 * (1.0 - yco2 - yM) / 7.0;
    if (yE < 0.0) yE = 0.0;
    const double yP = 0.5 * yE;
    double yB = 1.0 - yM - yE - yP - yco2;
    if (yB < 0.0) yB = 0.0;

    const double CP0 = 1000.0 * (yM * CP0mM + yE * CP0mE + yP * CP0mP +
                                 yB * CP0mB + yco2 * CP0mCO2) /
                       (rDgL * Deng * 28.9625);
    const double CV0 = CP0 - RG;

    const double zg = zFactor(pres, temp, Deng, PCis, TCis);
    const double TF = to_fahr(temp);
    const double TRan = TF + 459.67;
    const double Ppsi = to_psia(pres);
    const double PR = Ppsi / PCis;
    const double TR = TRan / TCis;
    const double rhoR = 0.27 * PR / zg / TR;
    const double dzdt = (zFactor(pres, temp + 1e-3, Deng, PCis, TCis) -
                         zFactor(pres, temp - 1e-3, Deng, PCis, TCis)) / (2e-3);
    const double dzdp = (zFactor(pres + 1e-3, temp, Deng, PCis, TCis) -
                         zFactor(pres - 1e-3, temp, Deng, PCis, TCis)) / (2e-3);

    const double CV = CV0 - 5.104577 * RG * (1.0 -
        (0.6275544 - 0.03688454 / TR + 0.2023674 / (TR * TR)) * rhoR +
        0.03612444 * TR * (rhoR * rhoR) -
        0.0015258965 * TR * std::pow(rhoR, 5.0) -
        std::exp(-0.7210 * (rhoR * rhoR)) * (1.0 + 0.3605 * (rhoR * rhoR))) /
        TR / (TR * TR);

    return CV + RG * std::pow(zg + TR * dzdt, 2.0) / (zg - PR * dzdp);
}

double liquidDensityDerivativeTBlackOil(double pres, double temp, double API,
                                        double Deng, double rs, double rDgD) {
    const double rho = oilDensityBlackOil(pres, temp, API, Deng, rs, rDgD);
    return -624.0 / rho;
}

} // namespace sisprod2
} // namespace marlim
