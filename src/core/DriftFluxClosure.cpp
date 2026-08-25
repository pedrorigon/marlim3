#include "DriftFluxClosure.h"

// Matches the include SisProd.cpp used when these functions lived there.
// It is load-bearing: BhagwatGhajar and BhagwatGhajarMod call abs() on double
// operands inside the Colebrook loop, and picking up int abs(int) instead
// would truncate frictionFactorEstimate and the convergence delta, ending the loop early
// with badly wrong values. The static_assert below fails the build if the
// overload ever stops being the floating-point one.
#include <math.h>

#include <type_traits>

static_assert(std::is_same<decltype(abs(1.5)), double>::value,
              "abs() must resolve to the double overload; see the Colebrook "
              "loop in BhagwatGhajar and BhagwatGhajarMod");

namespace driftflux {
namespace correlations {

void BhagwatGhajar(double liquidDensity, double gasDensity, double surfaceTension, double voidFraction,
                   double mixtureReynolds, double liquidReynolds, double gasFlowRate, double liquidFlowRate,
                   double diameter, double roughness, double inclinationAngle, double &c0, double &ud,
                   double horizontalCorrection) {

    double noSlipGasFraction, froudeNumber, massQuality, distributionTerm1, distributionTerm2, ductShapeTerm, ductShapeCoefficient;
    double inclinationFactor, buoyancyVelocityScale, laplaceNumber, viscosityCorrection, laplaceCorrection, downwardFlowSign;
    double flowArea = M_PI * diameter * diameter / 4.;
    double mixtureDensity = voidFraction * gasDensity + (1. - voidFraction) * liquidDensity;

    double inclinationSign = 1.;
    if (inclinationAngle < 0.)
        inclinationSign = 1.;
    if (inclinationAngle == 0) {
        int para;
        para = 1.;
    }

    double meanDensity = voidFraction * gasDensity + (1. - voidFraction) * liquidDensity;
    froudeNumber = sqrt(gasDensity / (liquidDensity - gasDensity)) * ((gasFlowRate) / flowArea) / sqrt(9.81 * diameter * inclinationSign * cos(inclinationAngle));
    massQuality = (gasDensity * fabs(gasFlowRate) / flowArea) / ((gasDensity * fabs(gasFlowRate) / flowArea) + (liquidDensity * fabs(liquidFlowRate) / flowArea));
    noSlipGasFraction = (fabs(gasFlowRate) / flowArea) / ((fabs(gasFlowRate) / flowArea) + (fabs(liquidFlowRate) / flowArea));

    double relativeRoughness;      // rugosidade relativa.
    relativeRoughness = roughness / diameter; // roughness - rugosidade absoluta. relativeRoughness - rugosidade relativa

    if (mixtureReynolds < 0.0000001)
        mixtureReynolds = 0.0000001;
    double frictionFactor, frictionFactorEstimate, colebrookDenominator, convergenceDelta;
    int iterationCount;
    if (mixtureReynolds > 2400) { // regime turbulento do escoamento
        frictionFactorEstimate = (1 / (-18e-1 * log10(pow((relativeRoughness / (3.7)), 1.11) + (69e-1 / (mixtureReynolds + 1e-15)))));
        frictionFactorEstimate *= frictionFactorEstimate; // Halland.
        iterationCount = 0;
    repeat7:
        iterationCount = iterationCount + 1;
        colebrookDenominator = -2 * log10(((relativeRoughness) / 3.7) + 2.51 / ((mixtureReynolds + 1e-15) * sqrt(abs(frictionFactorEstimate))));
        frictionFactor = 1 / (colebrookDenominator * colebrookDenominator); // Colebrook.
        convergenceDelta = abs(frictionFactor - frictionFactorEstimate);
        frictionFactorEstimate = frictionFactor;
        if (convergenceDelta >= 1e-3)
            goto repeat7;
    } else {                  // regime laminar
        frictionFactor = 64. / (mixtureReynolds); // 16.
    }

    double densityRatioSquared = gasDensity / liquidDensity;
    densityRatioSquared *= densityRatioSquared;
    double scaledReynoldsSquared = mixtureReynolds / 1000;
    scaledReynoldsSquared *= scaledReynoldsSquared;
    distributionTerm1 = (2 - densityRatioSquared) / (1 + scaledReynoldsSquared);
    distributionTerm2 = (pow(((1 + densityRatioSquared * inclinationSign * cos(inclinationAngle)) / (1 + cos(inclinationAngle))), (1 - voidFraction) / 5.)) /
             (1 + 1 / scaledReynoldsSquared);
    ductShapeCoefficient = 0.2; // duto circular ou anular. Retangular seria 0.4.
    ductShapeTerm = (ductShapeCoefficient - ductShapeCoefficient * sqrt(gasDensity / liquidDensity)) * (pow((2.6 - noSlipGasFraction), 0.15) - sqrt(frictionFactor)) * pow((1 - massQuality), 1.5);
    if (gasFlowRate * liquidFlowRate < 0.)
        ductShapeTerm = 0;
    if (inclinationAngle >= -50 * M_PI / 180. && inclinationAngle <= 0 && froudeNumber <= 0.1)
        ductShapeTerm = 0.0;
    c0 = distributionTerm1 + distributionTerm2 + ductShapeTerm; // Calculo do ParÃƒÂ¢metro de Distribuicao.

    double mixtureViscosity = diameter * (fabs(gasFlowRate / flowArea) + fabs(liquidFlowRate / flowArea)) * mixtureDensity / mixtureReynolds;
    inclinationFactor = (0.35 * sin(inclinationAngle) + 0.45 * cos(inclinationAngle) * inclinationSign);
    buoyancyVelocityScale = sqrt((9.81 * diameter * (liquidDensity - gasDensity) / liquidDensity)) * sqrt(1 - voidFraction);
    if (mixtureViscosity / 0.001 > 10) {
        viscosityCorrection = pow((0.434 / (log10(mixtureViscosity / 0.001))), 0.15);
    } else {
        viscosityCorrection = 1.0;
    }
    laplaceNumber = sqrt(surfaceTension / (9.81 * (liquidDensity - gasDensity))) / diameter;
    if (laplaceNumber < 0.025) {
        laplaceCorrection = pow((laplaceNumber / 0.025), 0.90);
    } else {
        laplaceCorrection = 1.0;
        ;
    }
    downwardFlowSign = 1.0;
    if (inclinationAngle >= -(50 * M_PI / 180.) && inclinationAngle < 0 && froudeNumber <= 0.1)
        downwardFlowSign = -1.0;
    ud = horizontalCorrection * inclinationFactor * buoyancyVelocityScale * viscosityCorrection * laplaceCorrection * downwardFlowSign; // Calculo da Velocidade de Deslizamento.
    if ((fabs(gasFlowRate / flowArea) + fabs(liquidFlowRate / flowArea)) < 0.01 && inclinationAngle > 0. && ud < 0.)
        ud = fabs(ud);
    else if ((fabs(gasFlowRate / flowArea) + fabs(liquidFlowRate / flowArea)) < 0.01 && inclinationAngle < 0. && ud > 0.)
        ud = -fabs(ud);
}

void BhagwatGhajarMod(double liquidDensity, double gasDensity, double surfaceTension, double voidFraction,
                      double mixtureReynolds, double liquidReynolds, double gasFlowRate,
                      double liquidFlowRate, double diameter, double roughness, double inclinationAngle,
                      double &c0, double &ud, double horizontalCorrection) {
    double noSlipGasFraction, froudeNumber, massQuality, distributionTerm1, distributionTerm2, ductShapeTerm, ductShapeCoefficient;
    double inclinationFactor, buoyancyVelocityScale, laplaceNumber, viscosityCorrection, laplaceCorrection, downwardFlowSign;
    double flowArea = M_PI * diameter * diameter / 4.;
    double mixtureDensity = voidFraction * gasDensity + (1. - voidFraction) * liquidDensity;

    double inclinationSign = 1.;
    if (inclinationAngle < 0.)
        inclinationSign = 1.;
    if (inclinationAngle == 0) {
        int para;
        para = 1.;
    }

    double meanDensity = voidFraction * gasDensity + (1. - voidFraction) * liquidDensity;
    froudeNumber = sqrt(gasDensity / (liquidDensity - gasDensity)) * ((gasFlowRate) / flowArea) / sqrt(9.81 * diameter * inclinationSign * cos(inclinationAngle));
    massQuality = (gasDensity * fabs(gasFlowRate) / flowArea) / ((gasDensity * fabs(gasFlowRate) / flowArea) + (liquidDensity * fabs(liquidFlowRate) / flowArea));
    noSlipGasFraction = (fabs(gasFlowRate) / flowArea) / ((fabs(gasFlowRate) / flowArea) + (fabs(liquidFlowRate) / flowArea));

    double relativeRoughness;      // rugosidade relativa.
    relativeRoughness = roughness / diameter; // roughness - rugosidade absoluta. relativeRoughness - rugosidade relativa

    if (liquidReynolds < 0.0000001)
        liquidReynolds = 0.0000001;
    double frictionFactor, frictionFactorEstimate, colebrookDenominator, convergenceDelta;
    int iterationCount;
    if (liquidReynolds > 2400) { // regime turbulento do escoamento
        frictionFactorEstimate = (1 / (-18e-1 * log10(pow((relativeRoughness / (3.7)), 1.11) + (69e-1 / (liquidReynolds + 1e-15)))));
        frictionFactorEstimate *= frictionFactorEstimate; // Halland.
        iterationCount = 0;
    repeat7:
        iterationCount = iterationCount + 1;
        colebrookDenominator = -2 * log10(((relativeRoughness) / 3.7) + 2.51 / ((liquidReynolds + 1e-15) * sqrt(abs(frictionFactorEstimate))));
        frictionFactor = 1 / (colebrookDenominator * colebrookDenominator); // Colebrook.
        convergenceDelta = abs(frictionFactor - frictionFactorEstimate);
        frictionFactorEstimate = frictionFactor;
        if (convergenceDelta >= 1e-3)
            goto repeat7;
    } else {                   // regime laminar
        frictionFactor = 64. / (liquidReynolds); // 16.
    }

    double densityRatioSquared = gasDensity / liquidDensity;
    densityRatioSquared *= densityRatioSquared;
    double scaledReynoldsSquared = liquidReynolds / 1000;
    scaledReynoldsSquared *= scaledReynoldsSquared;
    distributionTerm1 = (2 - densityRatioSquared) / (1 + scaledReynoldsSquared);
    distributionTerm2 = (pow(((1 + densityRatioSquared * inclinationSign * cos(inclinationAngle)) / (1 + cos(inclinationAngle))), (1 - voidFraction) / 5.)) /
             (1 + 1 / scaledReynoldsSquared);
    ductShapeCoefficient = 0.2; // duto circular ou anular. Retangular seria 0.4.
    ductShapeTerm = (ductShapeCoefficient - ductShapeCoefficient * sqrt(gasDensity / liquidDensity)) * (pow((2.6 - noSlipGasFraction), 0.15) - sqrt(frictionFactor)) * pow((1 - massQuality), 1.5);
    if (gasFlowRate * liquidFlowRate < 0.)
        ductShapeTerm = 0;
    if (inclinationAngle >= -50 * M_PI / 180. && inclinationAngle <= 0 && froudeNumber <= 0.1)
        ductShapeTerm = 0.0;
    c0 = distributionTerm1 + distributionTerm2 + ductShapeTerm; // Calculo do ParÃƒÂ¢metro de Distribuicao.

    double mixtureViscosity = diameter * (fabs(gasFlowRate / flowArea) + fabs(liquidFlowRate / flowArea)) * mixtureDensity / liquidReynolds;
    inclinationFactor = (0.35 * sin(inclinationAngle) + 0.45 * cos(inclinationAngle) * inclinationSign);
    buoyancyVelocityScale = sqrt((9.81 * diameter * (liquidDensity - gasDensity) / liquidDensity)) * sqrt(1 - voidFraction);
    if (mixtureViscosity / 0.001 > 10) {
        viscosityCorrection = pow((0.434 / (log10(mixtureViscosity / 0.001))), 0.15);
    } else {
        viscosityCorrection = 1.0;
    }
    laplaceNumber = sqrt(surfaceTension / (9.81 * (liquidDensity - gasDensity))) / diameter;
    if (laplaceNumber < 0.025) {
        laplaceCorrection = pow((laplaceNumber / 0.025), 0.90);
    } else {
        laplaceCorrection = 1.0;
        ;
    }
    downwardFlowSign = 1.0;
    if (inclinationAngle >= -(50 * M_PI / 180.) && inclinationAngle < 0 && froudeNumber <= 0.1)
        downwardFlowSign = -1.0;
    ud = horizontalCorrection * inclinationFactor * buoyancyVelocityScale * viscosityCorrection * laplaceCorrection * downwardFlowSign; // Calculo da Velocidade de Deslizamento.
    if ((fabs(gasFlowRate / flowArea) + fabs(liquidFlowRate / flowArea)) < 0.01 && inclinationAngle > 0. && ud < 0.)
        ud = fabs(ud);
    else if ((fabs(gasFlowRate / flowArea) + fabs(liquidFlowRate / flowArea)) < 0.01 && inclinationAngle < 0. && ud > 0.)
        ud = -fabs(ud);
}

void Choi(double liquidDensity, double gasDensity, double surfaceTension, double voidFraction,
          double mixtureReynolds, double liquidReynolds, double gasFlowRate, double liquidFlowRate,
          double diameter, double roughness, double inclinationAngle, double &c0, double &ud,
          double horizontalCorrection) {

    double inclinationSign = 1.;
    if (inclinationAngle < 0.)
        inclinationSign = -1.;
    double flowArea = M_PI * diameter * diameter / 4.;
    ud = horizontalCorrection * inclinationSign * 0.0246 * cos(inclinationAngle) + 1.606 * pow(9.82 * surfaceTension * (liquidDensity - gasDensity) / (liquidDensity * liquidDensity), 0.25) * sin(inclinationAngle);
    c0 = 2. / (1 + pow(mixtureReynolds / 1000., 2.)) + (1.2 - 0.2 * sqrt(gasDensity / liquidDensity) * (1 - exp(-18 * voidFraction))) / (1 + pow(1000. / mixtureReynolds, 2.));
    if ((fabs(gasFlowRate / flowArea) + fabs(liquidFlowRate / flowArea)) < 0.01 && inclinationAngle > 0. && ud < 0.)
        ud = fabs(ud);
    else if ((fabs(gasFlowRate / flowArea) + fabs(liquidFlowRate / flowArea)) < 0.01 && inclinationAngle < 0. && ud > 0.)
        ud = -fabs(ud);
}

void HibikiIshii(double liquidDensity, double gasDensity, double surfaceTension, double voidFraction,
                 double mixtureReynolds, double liquidReynolds, double gasFlowRate, double liquidFlowRate,
                 double diameter, double roughness, double inclinationAngle, double &c0, double &ud,
                 double horizontalCorrection) {

    double inclinationSign = 1.;
    if (inclinationAngle < 0.) {
        inclinationSign = -1.;
    }
    double flowArea = M_PI * diameter * diameter / 4.;
    c0 = 1. + (1. - voidFraction) / (voidFraction + 4. * sqrt(gasDensity / liquidDensity));
    ud = (horizontalCorrection * inclinationSign * (1. - voidFraction) / (voidFraction + 4. * sqrt(gasDensity / liquidDensity))) * sqrt(9.82 * fabs(sin(inclinationAngle)) * diameter * (liquidDensity - gasDensity) * (1. - voidFraction) / (0.015 * liquidDensity));
    if ((fabs(gasFlowRate / flowArea) + fabs(liquidFlowRate / flowArea)) < 0.01 && inclinationAngle > 0. && ud < 0.)
        ud = fabs(ud);
    else if ((fabs(gasFlowRate / flowArea) + fabs(liquidFlowRate / flowArea)) < 0.01 && inclinationAngle < 0. && ud > 0.)
        ud = -fabs(ud);
}

void FrancaLahey(double liquidDensity, double gasDensity, double surfaceTension, double voidFraction,
                 double mixtureReynolds, double liquidReynolds, double gasFlowRate, double liquidFlowRate,
                 double diameter, double roughness, double inclinationAngle, double &c0, double &ud,
                 double horizontalCorrection) {

    double inclinationSign = 1.;
    if (inclinationAngle < 0.)
        inclinationSign = -1.;
    c0 = 1.04;
    ud = horizontalCorrection * inclinationSign * 0.466;
    double flowArea = M_PI * diameter * diameter / 4.;
    if ((fabs(gasFlowRate / flowArea) + fabs(liquidFlowRate / flowArea)) < 0.01 && inclinationAngle > 0. && ud < 0.)
        ud = fabs(ud);
    else if ((fabs(gasFlowRate / flowArea) + fabs(liquidFlowRate / flowArea)) < 0.01 && inclinationAngle < 0. && ud > 0.)
        ud = -fabs(ud);
}

namespace {

/// Selector 5: blend two correlations across the inclination band between 5 and
/// 20 degrees. Identical in all three regimes, and the most delicate arithmetic
/// in this file, which is why it lives in one place now instead of three.
///
/// Two things here are load-bearing and must not be tidied. BhagwatGhajarMod
/// runs BEFORE Choi because both write c0 and ud, so the order decides which
/// result the temporaries hold. And the interpolation keeps the form
/// blendRatio*x + (1 - blendRatio)*xMod: rearranging it to
/// xMod + blendRatio*(x - xMod) is algebraically identical and rounds
/// differently.
void blendAcrossInclination(double liquidDensity, double gasDensity, double surfaceTension, double voidFraction,
                      double mixtureReynolds, double liquidReynolds, double gasFlowRate,
                      double liquidFlowRate, double diameter, double roughness,
                      double inclinationAngle, double &c0, double &ud,
                      double horizontalCorrection) {
    if (fabs(inclinationAngle) < 5 * M_PI / 180.)
        BhagwatGhajar(liquidDensity, gasDensity, surfaceTension, voidFraction, mixtureReynolds,
                  liquidReynolds, gasFlowRate, liquidFlowRate, diameter, roughness,
                  inclinationAngle, c0, ud, horizontalCorrection);
    else if (fabs(inclinationAngle) > 20 * M_PI / 180.)
        Choi(liquidDensity, gasDensity, surfaceTension, voidFraction, mixtureReynolds,
                  liquidReynolds, gasFlowRate, liquidFlowRate, diameter, roughness,
                  inclinationAngle, c0, ud, horizontalCorrection);
    else {
        double blendRatio = (fabs(inclinationAngle) - 5 * M_PI / 180.) / (15 * M_PI / 180.);
        BhagwatGhajarMod(liquidDensity, gasDensity, surfaceTension, voidFraction, mixtureReynolds,
                  liquidReynolds, gasFlowRate, liquidFlowRate, diameter, roughness,
                  inclinationAngle, c0, ud, horizontalCorrection);
        double c0BhagwatGhajarMod = c0;
        double udBhagwatGhajarMod = ud;
        Choi(liquidDensity, gasDensity, surfaceTension, voidFraction, mixtureReynolds,
                  liquidReynolds, gasFlowRate, liquidFlowRate, diameter, roughness,
                  inclinationAngle, c0, ud, horizontalCorrection);
        c0 = blendRatio * c0 + (1. - blendRatio) * c0BhagwatGhajarMod;
        ud = blendRatio * ud + (1. - blendRatio) * udBhagwatGhajarMod;
    }
}

/// The selectors every regime accepts. Anything else falls through and leaves
/// c0 and ud exactly as the caller passed them -- the original switches carried
/// no default, and that silence is observable behaviour rather than an
/// oversight, so it is preserved deliberately.
void applyCommonCorrelation(int correlationIndex, double liquidDensity, double gasDensity, double surfaceTension, double voidFraction,
                      double mixtureReynolds, double liquidReynolds, double gasFlowRate,
                      double liquidFlowRate, double diameter, double roughness,
                      double inclinationAngle, double &c0, double &ud,
                      double horizontalCorrection) {
    switch (correlationIndex) {
    case 0:
        Choi(liquidDensity, gasDensity, surfaceTension, voidFraction, mixtureReynolds,
                  liquidReynolds, gasFlowRate, liquidFlowRate, diameter, roughness,
                  inclinationAngle, c0, ud, horizontalCorrection);
        break;
    case 1:
        BhagwatGhajar(liquidDensity, gasDensity, surfaceTension, voidFraction, mixtureReynolds,
                  liquidReynolds, gasFlowRate, liquidFlowRate, diameter, roughness,
                  inclinationAngle, c0, ud, horizontalCorrection);
        break;
    case 4:
        BhagwatGhajarMod(liquidDensity, gasDensity, surfaceTension, voidFraction, mixtureReynolds,
                  liquidReynolds, gasFlowRate, liquidFlowRate, diameter, roughness,
                  inclinationAngle, c0, ud, horizontalCorrection);
        break;
    case 5:
        blendAcrossInclination(liquidDensity, gasDensity, surfaceTension, voidFraction, mixtureReynolds,
                  liquidReynolds, gasFlowRate, liquidFlowRate, diameter, roughness,
                  inclinationAngle, c0, ud, horizontalCorrection);
        break;
    }
}

}  // namespace

void C0UdDisperso(double liquidDensity, double gasDensity, double surfaceTension, double voidFraction,
                  double mixtureReynolds, double liquidReynolds, double gasFlowRate,
                  double liquidFlowRate, double diameter, double roughness,
                  double inclinationAngle, double &c0, double &ud,
                  double horizontalCorrection, int estabCol, int correlationIndex) {
    applyCommonCorrelation(correlationIndex, liquidDensity, gasDensity, surfaceTension, voidFraction, mixtureReynolds,
                            liquidReynolds, gasFlowRate, liquidFlowRate, diameter, roughness,
                            inclinationAngle, c0, ud, horizontalCorrection);
}

void C0UdAnularChurn(double liquidDensity, double gasDensity, double surfaceTension, double voidFraction,
                  double mixtureReynolds, double liquidReynolds, double gasFlowRate,
                  double liquidFlowRate, double diameter, double roughness,
                  double inclinationAngle, double &c0, double &ud,
                  double horizontalCorrection, int estabCol, int correlationIndex) {
    // Hibiki-Ishii is accepted in this regime and in no other.
    if (correlationIndex == 3) {
        HibikiIshii(liquidDensity, gasDensity, surfaceTension, voidFraction, mixtureReynolds,
                  liquidReynolds, gasFlowRate, liquidFlowRate, diameter, roughness,
                  inclinationAngle, c0, ud, horizontalCorrection);
        return;
    }
    applyCommonCorrelation(correlationIndex, liquidDensity, gasDensity, surfaceTension, voidFraction, mixtureReynolds,
                            liquidReynolds, gasFlowRate, liquidFlowRate, diameter, roughness,
                            inclinationAngle, c0, ud, horizontalCorrection);
}

void C0UdEstratificado(double liquidDensity, double gasDensity, double surfaceTension, double voidFraction,
                  double mixtureReynolds, double liquidReynolds, double gasFlowRate,
                  double liquidFlowRate, double diameter, double roughness,
                  double inclinationAngle, double &c0, double &ud,
                  double horizontalCorrection, int estabCol, int correlationIndex) {
    // Franca-Lahey is accepted in this regime and in no other.
    if (correlationIndex == 2) {
        FrancaLahey(liquidDensity, gasDensity, surfaceTension, voidFraction, mixtureReynolds,
                  liquidReynolds, gasFlowRate, liquidFlowRate, diameter, roughness,
                  inclinationAngle, c0, ud, horizontalCorrection);
        return;
    }
    applyCommonCorrelation(correlationIndex, liquidDensity, gasDensity, surfaceTension, voidFraction, mixtureReynolds,
                            liquidReynolds, gasFlowRate, liquidFlowRate, diameter, roughness,
                            inclinationAngle, c0, ud, horizontalCorrection);
}

}  // namespace correlations
}  // namespace driftflux
