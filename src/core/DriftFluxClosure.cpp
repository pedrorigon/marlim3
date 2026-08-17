#include "DriftFluxClosure.h"

// Matches the include SisProd.cpp used when these functions lived there.
// It is load-bearing: BhagwatGhajar and BhagwatGhajarMod call abs() on double
// operands inside the Colebrook loop, and picking up int abs(int) instead
// would truncate valHalland and the convergence delta, ending the loop early
// with badly wrong values. The static_assert below fails the build if the
// overload ever stops being the floating-point one.
#include <math.h>

#include <type_traits>

static_assert(std::is_same<decltype(abs(1.5)), double>::value,
              "abs() must resolve to the double overload; see the Colebrook "
              "loop in BhagwatGhajar and BhagwatGhajarMod");

namespace driftflux {
namespace correlations {

void BhagwatGhajar(double rhol, double rhog, double tensup, double alf, double reymix, double reymixL,
                          double ug1, double ul1, double dia, double rug, double tet, double &c0,
                          double &ud, double correcHor) {

    double Beta, Froude, MassQuality, termo1, termo2, C01, C1;
    double termoud1, termoud2, La, C2, C3, C4;
    double A1 = M_PI * dia * dia / 4.;
    double rhomix = alf * rhog + (1. - alf) * rhol;

    double sinal = 1.;
    if (tet < 0.)
        sinal = 1.;
    if (tet == 0) {
        int para;
        para = 1.;
    }

    double rmed = alf * rhog + (1. - alf) * rhol;
    Froude = sqrt(rhog / (rhol - rhog)) * ((ug1) / A1) / sqrt(9.81 * dia * sinal * cos(tet));
    MassQuality = (rhog * fabs(ug1) / A1) / ((rhog * fabs(ug1) / A1) + (rhol * fabs(ul1) / A1));
    Beta = (fabs(ug1) / A1) / ((fabs(ug1) / A1) + (fabs(ul1) / A1));

    double eps;      // rugosidade relativa.
    eps = rug / dia; // rug - rugosidade absoluta. eps - rugosidade relativa

    if (reymix < 0.0000001)
        reymix = 0.0000001;
    double fat, valHalland, den, dif;
    int III;
    if (reymix > 2400) { // regime turbulento do escoamento
        valHalland = (1 / (-18e-1 * log10(pow((eps / (3.7)), 1.11) + (69e-1 / (reymix + 1e-15)))));
        valHalland *= valHalland; // Halland.
        III = 0;
    repeat7:
        III = III + 1;
        den = -2 * log10(((eps) / 3.7) + 2.51 / ((reymix + 1e-15) * sqrt(abs(valHalland))));
        fat = 1 / (den * den); // Colebrook.
        dif = abs(fat - valHalland);
        valHalland = fat;
        if (dif >= 1e-3)
            goto repeat7;
    } else {                  // regime laminar
        fat = 64. / (reymix); // 16.
    }

    double rgrl2 = rhog / rhol;
    rgrl2 *= rgrl2;
    double reymix2 = reymix / 1000;
    reymix2 *= reymix2;
    termo1 = (2 - rgrl2) / (1 + reymix2);
    termo2 = (pow(((1 + rgrl2 * sinal * cos(tet)) / (1 + cos(tet))), (1 - alf) / 5.)) /
             (1 + 1 / reymix2);
    C1 = 0.2; // duto circular ou anular. Retangular seria 0.4.
    C01 = (C1 - C1 * sqrt(rhog / rhol)) * (pow((2.6 - Beta), 0.15) - sqrt(fat)) * pow((1 - MassQuality), 1.5);
    if (ug1 * ul1 < 0.)
        C01 = 0;
    if (tet >= -50 * M_PI / 180. && tet <= 0 && Froude <= 0.1)
        C01 = 0.0;
    c0 = termo1 + termo2 + C01; // Calculo do ParÃƒÂ¢metro de Distribuicao.

    double viscl1 = dia * (fabs(ug1 / A1) + fabs(ul1 / A1)) * rhomix / reymix;
    termoud1 = (0.35 * sin(tet) + 0.45 * cos(tet) * sinal);
    termoud2 = sqrt((9.81 * dia * (rhol - rhog) / rhol)) * sqrt(1 - alf);
    if (viscl1 / 0.001 > 10) {
        C2 = pow((0.434 / (log10(viscl1 / 0.001))), 0.15);
    } else {
        C2 = 1.0;
    }
    La = sqrt(tensup / (9.81 * (rhol - rhog))) / dia;
    if (La < 0.025) {
        C3 = pow((La / 0.025), 0.90);
    } else {
        C3 = 1.0;
        ;
    }
    C4 = 1.0;
    if (tet >= -(50 * M_PI / 180.) && tet < 0 && Froude <= 0.1)
        C4 = -1.0;
    ud = correcHor * termoud1 * termoud2 * C2 * C3 * C4; // Calculo da Velocidade de Deslizamento.
    if ((fabs(ug1 / A1) + fabs(ul1 / A1)) < 0.01 && tet > 0. && ud < 0.)
        ud = fabs(ud);
    else if ((fabs(ug1 / A1) + fabs(ul1 / A1)) < 0.01 && tet < 0. && ud > 0.)
        ud = -fabs(ud);
}

void BhagwatGhajarMod(double rhol, double rhog, double tensup, double alf, double reymix, double reymixL,
                             double ug1, double ul1, double dia, double rug, double tet, double &c0,
                             double &ud, double correcHor) {
    double Beta, Froude, MassQuality, termo1, termo2, C01, C1;
    double termoud1, termoud2, La, C2, C3, C4;
    double A1 = M_PI * dia * dia / 4.;
    double rhomix = alf * rhog + (1. - alf) * rhol;

    double sinal = 1.;
    if (tet < 0.)
        sinal = 1.;
    if (tet == 0) {
        int para;
        para = 1.;
    }

    double rmed = alf * rhog + (1. - alf) * rhol;
    Froude = sqrt(rhog / (rhol - rhog)) * ((ug1) / A1) / sqrt(9.81 * dia * sinal * cos(tet));
    MassQuality = (rhog * fabs(ug1) / A1) / ((rhog * fabs(ug1) / A1) + (rhol * fabs(ul1) / A1));
    Beta = (fabs(ug1) / A1) / ((fabs(ug1) / A1) + (fabs(ul1) / A1));

    double eps;      // rugosidade relativa.
    eps = rug / dia; // rug - rugosidade absoluta. eps - rugosidade relativa

    if (reymixL < 0.0000001)
        reymixL = 0.0000001;
    double fat, valHalland, den, dif;
    int III;
    if (reymixL > 2400) { // regime turbulento do escoamento
        valHalland = (1 / (-18e-1 * log10(pow((eps / (3.7)), 1.11) + (69e-1 / (reymixL + 1e-15)))));
        valHalland *= valHalland; // Halland.
        III = 0;
    repeat7:
        III = III + 1;
        den = -2 * log10(((eps) / 3.7) + 2.51 / ((reymixL + 1e-15) * sqrt(abs(valHalland))));
        fat = 1 / (den * den); // Colebrook.
        dif = abs(fat - valHalland);
        valHalland = fat;
        if (dif >= 1e-3)
            goto repeat7;
    } else {                   // regime laminar
        fat = 64. / (reymixL); // 16.
    }

    double rgrl2 = rhog / rhol;
    rgrl2 *= rgrl2;
    double reymix2 = reymixL / 1000;
    reymix2 *= reymix2;
    termo1 = (2 - rgrl2) / (1 + reymix2);
    termo2 = (pow(((1 + rgrl2 * sinal * cos(tet)) / (1 + cos(tet))), (1 - alf) / 5.)) /
             (1 + 1 / reymix2);
    C1 = 0.2; // duto circular ou anular. Retangular seria 0.4.
    C01 = (C1 - C1 * sqrt(rhog / rhol)) * (pow((2.6 - Beta), 0.15) - sqrt(fat)) * pow((1 - MassQuality), 1.5);
    if (ug1 * ul1 < 0.)
        C01 = 0;
    if (tet >= -50 * M_PI / 180. && tet <= 0 && Froude <= 0.1)
        C01 = 0.0;
    c0 = termo1 + termo2 + C01; // Calculo do ParÃƒÂ¢metro de Distribuicao.

    double viscl1 = dia * (fabs(ug1 / A1) + fabs(ul1 / A1)) * rhomix / reymixL;
    termoud1 = (0.35 * sin(tet) + 0.45 * cos(tet) * sinal);
    termoud2 = sqrt((9.81 * dia * (rhol - rhog) / rhol)) * sqrt(1 - alf);
    if (viscl1 / 0.001 > 10) {
        C2 = pow((0.434 / (log10(viscl1 / 0.001))), 0.15);
    } else {
        C2 = 1.0;
    }
    La = sqrt(tensup / (9.81 * (rhol - rhog))) / dia;
    if (La < 0.025) {
        C3 = pow((La / 0.025), 0.90);
    } else {
        C3 = 1.0;
        ;
    }
    C4 = 1.0;
    if (tet >= -(50 * M_PI / 180.) && tet < 0 && Froude <= 0.1)
        C4 = -1.0;
    ud = correcHor * termoud1 * termoud2 * C2 * C3 * C4; // Calculo da Velocidade de Deslizamento.
    if ((fabs(ug1 / A1) + fabs(ul1 / A1)) < 0.01 && tet > 0. && ud < 0.)
        ud = fabs(ud);
    else if ((fabs(ug1 / A1) + fabs(ul1 / A1)) < 0.01 && tet < 0. && ud > 0.)
        ud = -fabs(ud);
}

void Choi(double rhol, double rhog, double tensup, double alf, double reymix, double reymixL,
                 double ug1, double ul1, double dia, double rug, double tet, double &c0,
                 double &ud, double correcHor) {

    double sinal = 1.;
    if (tet < 0.)
        sinal = -1.;
    double A1 = M_PI * dia * dia / 4.;
    ud = correcHor * sinal * 0.0246 * cos(tet) + 1.606 * pow(9.82 * tensup * (rhol - rhog) / (rhol * rhol), 0.25) * sin(tet);
    c0 = 2. / (1 + pow(reymix / 1000., 2.)) + (1.2 - 0.2 * sqrt(rhog / rhol) * (1 - exp(-18 * alf))) / (1 + pow(1000. / reymix, 2.));
    if ((fabs(ug1 / A1) + fabs(ul1 / A1)) < 0.01 && tet > 0. && ud < 0.)
        ud = fabs(ud);
    else if ((fabs(ug1 / A1) + fabs(ul1 / A1)) < 0.01 && tet < 0. && ud > 0.)
        ud = -fabs(ud);
}

void HibikiIshii(double rhol, double rhog, double tensup, double alf, double reymix, double reymixL,
                        double ug1, double ul1, double dia, double rug, double tet, double &c0,
                        double &ud, double correcHor) {

    double sinal = 1.;
    if (tet < 0.) {
        sinal = -1.;
    }
    double A1 = M_PI * dia * dia / 4.;
    c0 = 1. + (1. - alf) / (alf + 4. * sqrt(rhog / rhol));
    ud = (correcHor * sinal * (1. - alf) / (alf + 4. * sqrt(rhog / rhol))) * sqrt(9.82 * fabs(sin(tet)) * dia * (rhol - rhog) * (1. - alf) / (0.015 * rhol));
    if ((fabs(ug1 / A1) + fabs(ul1 / A1)) < 0.01 && tet > 0. && ud < 0.)
        ud = fabs(ud);
    else if ((fabs(ug1 / A1) + fabs(ul1 / A1)) < 0.01 && tet < 0. && ud > 0.)
        ud = -fabs(ud);
}

void FrancaLahey(double rhol, double rhog, double tensup, double alf, double reymix, double reymixL,
                        double ug1, double ul1, double dia, double rug, double tet, double &c0,
                        double &ud, double correcHor) {

    double sinal = 1.;
    if (tet < 0.)
        sinal = -1.;
    c0 = 1.04;
    ud = correcHor * sinal * 0.466;
    double A1 = M_PI * dia * dia / 4.;
    if ((fabs(ug1 / A1) + fabs(ul1 / A1)) < 0.01 && tet > 0. && ud < 0.)
        ud = fabs(ud);
    else if ((fabs(ug1 / A1) + fabs(ul1 / A1)) < 0.01 && tet < 0. && ud > 0.)
        ud = -fabs(ud);
}

void C0UdDisperso(double rhol, double rhog, double tensup, double alf, double reymix, double reymixL,
                         double ug1, double ul1, double dia, double rug, double tet, double &c0,
                         double &ud, double correcHor, int estabCol, int correlationIndex) {

    // case 0:
    switch (correlationIndex) {
    case 0:
        Choi(rhol, rhog, tensup, alf, reymix, reymixL, ug1, ul1, dia, rug, tet, c0,
             ud, correcHor);
        break;
    case 1:
        BhagwatGhajar(rhol, rhog, tensup, alf, reymix, reymixL, ug1, ul1, dia, rug, tet,
                      c0, ud, correcHor);
        break;
    case 4:
        BhagwatGhajarMod(rhol, rhog, tensup, alf, reymix, reymixL, ug1, ul1, dia, rug, tet,
                         c0, ud, correcHor);
        break;
    case 5:
        if (fabs(tet) < 5 * M_PI / 180.)
            BhagwatGhajar(rhol, rhog, tensup, alf, reymix, reymixL, ug1, ul1, dia, rug, tet,
                          c0, ud, correcHor);
        else if (fabs(tet) > 20 * M_PI / 180.)
            Choi(rhol, rhog, tensup, alf, reymix, reymixL, ug1, ul1, dia, rug, tet, c0,
                 ud, correcHor);
        else {
            double raz = (fabs(tet) - 5 * M_PI / 180.) / (15 * M_PI / 180.);
            BhagwatGhajarMod(rhol, rhog, tensup, alf, reymix, reymixL, ug1, ul1, dia, rug, tet,
                             c0, ud, correcHor);
            double c0temp = c0;
            double udtemp = ud;
            Choi(rhol, rhog, tensup, alf, reymix, reymixL, ug1, ul1, dia, rug, tet, c0,
                 ud, correcHor);
            c0 = raz * c0 + (1. - raz) * c0temp;
            ud = raz * ud + (1. - raz) * udtemp;
        }
        break;
    }
}
void C0UdAnularChurn(double rhol, double rhog, double tensup, double alf, double reymix, double reymixL,
                            double ug1, double ul1, double dia, double rug, double tet, double &c0,
                            double &ud, double correcHor, int estabCol, int correlationIndex) {

    switch (correlationIndex) {
    case 3:
        HibikiIshii(rhol, rhog, tensup, alf, reymix, reymixL, ug1, ul1, dia, rug, tet,
                    c0, ud, correcHor);
        break;
    case 0:
        Choi(rhol, rhog, tensup, alf, reymix, reymixL, ug1, ul1, dia, rug, tet,
             c0, ud, correcHor);
        break;
    case 1:
        BhagwatGhajar(rhol, rhog, tensup, alf, reymix, reymixL, ug1, ul1, dia, rug, tet,
                      c0, ud, correcHor);
        break;
    case 4:
        BhagwatGhajarMod(rhol, rhog, tensup, alf, reymix, reymixL, ug1, ul1, dia, rug, tet,
                         c0, ud, correcHor);
        break;
    case 5:
        if (fabs(tet) < 5 * M_PI / 180.)
            BhagwatGhajar(rhol, rhog, tensup, alf, reymix, reymixL, ug1, ul1, dia, rug, tet,
                          c0, ud, correcHor);
        else if (fabs(tet) > 20 * M_PI / 180.)
            Choi(rhol, rhog, tensup, alf, reymix, reymixL, ug1, ul1, dia, rug, tet, c0,
                 ud, correcHor);
        else {
            double raz = (fabs(tet) - 5 * M_PI / 180.) / (15 * M_PI / 180.);
            BhagwatGhajarMod(rhol, rhog, tensup, alf, reymix, reymixL, ug1, ul1, dia, rug, tet,
                             c0, ud, correcHor);
            double c0temp = c0;
            double udtemp = ud;
            Choi(rhol, rhog, tensup, alf, reymix, reymixL, ug1, ul1, dia, rug, tet, c0,
                 ud, correcHor);
            c0 = raz * c0 + (1. - raz) * c0temp;
            ud = raz * ud + (1. - raz) * udtemp;
        }
        break;
    }
}
void C0UdEstratificado(double rhol, double rhog, double tensup, double alf, double reymix, double reymixL,
                              double ug1, double ul1, double dia, double rug, double tet, double &c0,
                              double &ud, double correcHor, int estabCol, int correlationIndex) {
    // case 0:
    switch (correlationIndex) {
    case (2):
        FrancaLahey(rhol, rhog, tensup, alf, reymix, reymixL, ug1, ul1, dia, rug, tet,
                    c0, ud, correcHor);
        break;
    case (0):
        Choi(rhol, rhog, tensup, alf, reymix, reymixL, ug1, ul1, dia, rug, tet,
             c0, ud, correcHor);
        break;
    case (1):
        BhagwatGhajar(rhol, rhog, tensup, alf, reymix, reymixL, ug1, ul1, dia, rug, tet,
                      c0, ud, correcHor);
        break;
    case (4):
        BhagwatGhajarMod(rhol, rhog, tensup, alf, reymix, reymixL, ug1, ul1, dia, rug, tet,
                         c0, ud, correcHor);
        break;
    case 5:
        if (fabs(tet) < 5 * M_PI / 180.)
            BhagwatGhajar(rhol, rhog, tensup, alf, reymix, reymixL, ug1, ul1, dia, rug, tet,
                          c0, ud, correcHor);
        else if (fabs(tet) > 20 * M_PI / 180.)
            Choi(rhol, rhog, tensup, alf, reymix, reymixL, ug1, ul1, dia, rug, tet, c0,
                 ud, correcHor);
        else {
            double raz = (fabs(tet) - 5 * M_PI / 180.) / (15 * M_PI / 180.);
            BhagwatGhajarMod(rhol, rhog, tensup, alf, reymix, reymixL, ug1, ul1, dia, rug, tet,
                             c0, ud, correcHor);
            double c0temp = c0;
            double udtemp = ud;
            Choi(rhol, rhog, tensup, alf, reymix, reymixL, ug1, ul1, dia, rug, tet, c0,
                 ud, correcHor);
            c0 = raz * c0 + (1. - raz) * c0temp;
            ud = raz * ud + (1. - raz) * udtemp;
        }
        break;
    }
}

}  // namespace correlations
}  // namespace driftflux
