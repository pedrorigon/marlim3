// Standalone sweep of the three regime aggregators, every selector value.
// Links only DriftFluxClosure.cpp -- no simulator, no input files. This is
// possible only because the namespace is pure, and it covers the selectors the
// demo corpus never reaches (it pins CorreDisper to 1).
#include "DriftFluxClosure.h"
#include <cstdio>

using namespace driftflux::correlations;

int main() {
    const double liquidDensities[] = {700.0, 850.0, 1000.0};
    const double gasDensities[]    = {5.0, 50.0, 200.0};
    const double surfaceTensions[] = {0.005, 0.02, 0.06};
    const double voidFractions[]   = {0.01, 0.25, 0.5, 0.75, 0.99};
    const double reynolds[]        = {1.0e2, 1.0e4, 1.0e6};
    // Angles chosen to straddle the blend thresholds at 5 and 20 degrees.
    const double angles[] = {-1.4, -0.35, -0.2, -0.05, 0.0, 0.05, 0.2, 0.35, 1.4};
    // -1 and 9 are outside every accepted set: they must leave c0/ud untouched.
    const int selectors[] = {-1, 0, 1, 2, 3, 4, 5, 9};

    for (double rhol : liquidDensities)
    for (double rhog : gasDensities)
    for (double sigma : surfaceTensions)
    for (double alpha : voidFractions)
    for (double rey : reynolds)
    for (double angle : angles)
    for (int sel : selectors) {
        const double reyL = rey * 0.5;
        const double qg = 0.3, ql = 1.7, dia = 0.15, rug = 4.5e-5, corr = 1.0;
        for (int which = 0; which < 3; ++which) {
            // Sentinels: if the selector is not accepted, these must survive.
            double c0 = -12345.0, ud = -54321.0;
            const char *name = "";
            switch (which) {
            case 0: name = "Disperso";
                C0UdDisperso(rhol, rhog, sigma, alpha, rey, reyL, qg, ql, dia, rug, angle, c0, ud, corr, 0, sel);
                break;
            case 1: name = "AnularChurn";
                C0UdAnularChurn(rhol, rhog, sigma, alpha, rey, reyL, qg, ql, dia, rug, angle, c0, ud, corr, 0, sel);
                break;
            case 2: name = "Estratificado";
                C0UdEstratificado(rhol, rhog, sigma, alpha, rey, reyL, qg, ql, dia, rug, angle, c0, ud, corr, 0, sel);
                break;
            }
            std::printf("%s %d %a %a %a %a %a %a %a %a\n", name, sel,
                        rhol, rhog, sigma, alpha, rey, angle, c0, ud);
        }
    }
    return 0;
}
