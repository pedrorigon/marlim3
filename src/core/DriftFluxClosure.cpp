#include "DriftFluxClosure.h"

// Matches the include SisProd.cpp used when these functions lived there.
// It is load-bearing: the Colebrook loop below calls abs() on double operands,
// and picking up int abs(int) instead would truncate the iterate and the
// convergence delta, ending the loop early with badly wrong values. The
// static_assert fails the build if the overload ever stops being the
// floating-point one.
#include <math.h>

#include <type_traits>

// Needed only by driftflux::coefficient below: the correlations part above
// reads nothing but its arguments. estrat.h and mapa.h bring the two
// flow-pattern map classes the variants construct on the stack.
#include "Leitura.h"
#include "celula3.h"
#include "estrat.h"
#include "mapa.h"
#include "variaveisGlobais1D.h"

static_assert(std::is_same<decltype(abs(1.5)), double>::value,
              "abs() must resolve to the double overload; see the Colebrook "
              "loop in bhagwatGhajarCore");

namespace driftflux {
namespace correlations {
namespace {

/// Sign carried by the duct inclination, negative for downward flow.
///
/// Choi, Hibiki Ishii and Franca Lahey all derive it the same way. It is a
/// plain comparison rather than copysign, which would return -1 for negative
/// zero where the original returns 1.
inline double inclinationSignOf(double inclinationAngle) {
    return inclinationAngle < 0. ? -1. : 1.;
}

/// Cross sectional area of a circular duct.
inline double ductArea(double diameter) {
    return M_PI * diameter * diameter / 4.;
}

/// Forces the drift velocity to point along the inclination when the flow is
/// nearly stagnant.
///
/// Below a combined superficial velocity of 0.01 m/s the correlations can
/// return a drift velocity whose sign contradicts the duct inclination, which
/// is unphysical: buoyancy drives gas upward. All five correlations carried an
/// identical copy of this guard, so it lives in one place now.
///
/// @param gasFlowRate      Volumetric gas flow rate.
/// @param liquidFlowRate   Volumetric liquid flow rate.
/// @param flowArea         Duct cross sectional area.
/// @param inclinationAngle Duct inclination in radians, negative downward.
/// @param ud               Drift velocity, corrected in place.
inline void alignDriftWithInclination(double gasFlowRate, double liquidFlowRate, double flowArea,
                                      double inclinationAngle, double &ud) {
    const double superficialVelocity = fabs(gasFlowRate / flowArea) + fabs(liquidFlowRate / flowArea);
    if (superficialVelocity < 0.01 && inclinationAngle > 0. && ud < 0.)
        ud = fabs(ud);
    else if (superficialVelocity < 0.01 && inclinationAngle < 0. && ud > 0.)
        ud = -fabs(ud);
}

/// Darcy friction factor, Haaland estimate refined by Colebrook iteration.
///
/// @param relativeRoughness Roughness divided by diameter.
/// @param reynoldsNumber    Reynolds number, already floored away from zero.
/// @return Darcy friction factor.
///
/// The fixed point converges in at most two iterations over a grid far wider
/// than anything the engine produces, Reynolds from 1e-9 to 1e12 and roughness
/// from smooth to 1e-2, so the bound below is insurance and never binds. It is
/// checked after the body runs, which keeps the loop a do while and leaves the
/// arithmetic untouched.
///
/// A stall was never the real hazard anyway: if the iterate turns into NaN the
/// convergence delta does too, and a NaN comparison is false, so the loop
/// exits on its own.
double darcyFrictionFactor(double relativeRoughness, double reynoldsNumber) {
    if (reynoldsNumber > 2400) { // regime turbulento do escoamento
        double frictionFactorEstimate =
            (1 / (-18e-1 * log10(pow((relativeRoughness / (3.7)), 1.11) + (69e-1 / (reynoldsNumber + 1e-15)))));
        frictionFactorEstimate *= frictionFactorEstimate; // Haaland.
        // Never reached in practice, see the note above.
        constexpr int kMaxColebrookIterations = 100;
        int iterations = 0;
        double frictionFactor;
        double convergenceDelta;
        do {
            const double colebrookDenominator =
                -2 * log10(((relativeRoughness) / 3.7) + 2.51 / ((reynoldsNumber + 1e-15) * sqrt(abs(frictionFactorEstimate))));
            frictionFactor = 1 / (colebrookDenominator * colebrookDenominator); // Colebrook.
            convergenceDelta = abs(frictionFactor - frictionFactorEstimate);
            frictionFactorEstimate = frictionFactor;
        } while (convergenceDelta >= 1e-3 && ++iterations < kMaxColebrookIterations);
        return frictionFactor;
    }
    return 64. / (reynoldsNumber); // 16.
}

/// Shared implementation of the two Bhagwat and Ghajar variants.
///
/// Bhagwat, S. M. and Ghajar, A. J. (2014), "A flow pattern independent drift
/// flux model based void fraction correlation for a wide range of gas liquid
/// two phase flow", International Journal of Multiphase Flow, vol. 59.
///
/// The published correlation builds the distribution parameter from three
/// additive terms and the drift velocity from an inclination factor, a
/// buoyancy scale and three correction factors. The two exported variants run
/// exactly this algorithm and differ only in which Reynolds number feeds the
/// friction factor and the Reynolds dependent terms, which is why they share
/// one body here instead of the eighty duplicated lines they used to be.
///
/// @param reynoldsNumber Reynolds number the variant selects, mixture or liquid.
///
/// @note This variant carries no flow direction sign, unlike Choi, Hibiki Ishii
///       and Franca Lahey. The original code declared one, then assigned 1 on
///       both branches of a test on the inclination, so it was always 1 and the
///       branch was dead. It has been folded away, which changes nothing.
///
///       Substituting -1, the reading the dead branch invites, is not a fix. It
///       drives the Froude radicand negative, and although the resulting NaN
///       never reaches c0 or ud, it makes both guards below permanently false,
///       so the duct shape term would stop being zeroed and downwardFlowSign
///       would stop flipping for downward low Froude flow. A published
///       correlation does not carry guards that can never fire, which is the
///       argument that this variant simply has no such sign.
///
///       Direction is handled here by downwardFlowSign, factor C4 of the
///       published model. See evidencia/estagio-1/achado-sinal-bhagwatghajar.md.
void bhagwatGhajarCore(double liquidDensity, double gasDensity, double surfaceTension, double voidFraction,
                       double reynoldsNumber, double gasFlowRate, double liquidFlowRate, double diameter,
                       double roughness, double inclinationAngle, double &c0, double &ud,
                       double horizontalCorrection) {
    const double flowArea = ductArea(diameter);
    const double mixtureDensity = voidFraction * gasDensity + (1. - voidFraction) * liquidDensity;


    const double froudeNumber = sqrt(gasDensity / (liquidDensity - gasDensity)) * ((gasFlowRate) / flowArea) / sqrt(9.81 * diameter * cos(inclinationAngle));
    const double massQuality = (gasDensity * fabs(gasFlowRate) / flowArea) / ((gasDensity * fabs(gasFlowRate) / flowArea) + (liquidDensity * fabs(liquidFlowRate) / flowArea));
    const double noSlipGasFraction = (fabs(gasFlowRate) / flowArea) / ((fabs(gasFlowRate) / flowArea) + (fabs(liquidFlowRate) / flowArea));

    const double relativeRoughness = roughness / diameter;

    if (reynoldsNumber < 0.0000001)
        reynoldsNumber = 0.0000001;
    const double frictionFactor = darcyFrictionFactor(relativeRoughness, reynoldsNumber);

    double densityRatioSquared = gasDensity / liquidDensity;
    densityRatioSquared *= densityRatioSquared;
    double scaledReynoldsSquared = reynoldsNumber / 1000;
    scaledReynoldsSquared *= scaledReynoldsSquared;
    const double distributionTerm1 = (2 - densityRatioSquared) / (1 + scaledReynoldsSquared);
    const double distributionTerm2 = (pow(((1 + densityRatioSquared * cos(inclinationAngle)) / (1 + cos(inclinationAngle))), (1 - voidFraction) / 5.)) /
             (1 + 1 / scaledReynoldsSquared);
    const double ductShapeCoefficient = 0.2; // duto circular ou anular. Retangular seria 0.4.
    double ductShapeTerm = (ductShapeCoefficient - ductShapeCoefficient * sqrt(gasDensity / liquidDensity)) * (pow((2.6 - noSlipGasFraction), 0.15) - sqrt(frictionFactor)) * pow((1 - massQuality), 1.5);
    if (gasFlowRate * liquidFlowRate < 0.)
        ductShapeTerm = 0;
    if (inclinationAngle >= -50 * M_PI / 180. && inclinationAngle <= 0 && froudeNumber <= 0.1)
        ductShapeTerm = 0.0;
    c0 = distributionTerm1 + distributionTerm2 + ductShapeTerm; // Calculo do Parametro de Distribuicao.

    const double mixtureViscosity = diameter * (fabs(gasFlowRate / flowArea) + fabs(liquidFlowRate / flowArea)) * mixtureDensity / reynoldsNumber;
    const double inclinationFactor = (0.35 * sin(inclinationAngle) + 0.45 * cos(inclinationAngle));
    const double buoyancyVelocityScale = sqrt((9.81 * diameter * (liquidDensity - gasDensity) / liquidDensity)) * sqrt(1 - voidFraction);
    const double viscosityCorrection =
        (mixtureViscosity / 0.001 > 10) ? pow((0.434 / (log10(mixtureViscosity / 0.001))), 0.15) : 1.0;
    const double laplaceNumber = sqrt(surfaceTension / (9.81 * (liquidDensity - gasDensity))) / diameter;
    const double laplaceCorrection =
        (laplaceNumber < 0.025) ? pow((laplaceNumber / 0.025), 0.90) : 1.0;
    const double downwardFlowSign =
        (inclinationAngle >= -(50 * M_PI / 180.) && inclinationAngle < 0 && froudeNumber <= 0.1) ? -1.0 : 1.0;
    ud = horizontalCorrection * inclinationFactor * buoyancyVelocityScale * viscosityCorrection * laplaceCorrection * downwardFlowSign; // Calculo da Velocidade de Deslizamento.
    alignDriftWithInclination(gasFlowRate, liquidFlowRate, flowArea, inclinationAngle, ud);
}

} // namespace

void BhagwatGhajar(double liquidDensity, double gasDensity, double surfaceTension, double voidFraction,
                   double mixtureReynolds, double liquidReynolds, double gasFlowRate, double liquidFlowRate,
                   double diameter, double roughness, double inclinationAngle, double &c0, double &ud,
                   double horizontalCorrection) {
    bhagwatGhajarCore(liquidDensity, gasDensity, surfaceTension, voidFraction, mixtureReynolds,
                      gasFlowRate, liquidFlowRate, diameter, roughness, inclinationAngle, c0, ud,
                      horizontalCorrection);
}

void BhagwatGhajarMod(double liquidDensity, double gasDensity, double surfaceTension, double voidFraction,
                      double mixtureReynolds, double liquidReynolds, double gasFlowRate,
                      double liquidFlowRate, double diameter, double roughness, double inclinationAngle,
                      double &c0, double &ud, double horizontalCorrection) {
    bhagwatGhajarCore(liquidDensity, gasDensity, surfaceTension, voidFraction, liquidReynolds,
                      gasFlowRate, liquidFlowRate, diameter, roughness, inclinationAngle, c0, ud,
                      horizontalCorrection);
}

void Choi(double liquidDensity, double gasDensity, double surfaceTension, double voidFraction,
          double mixtureReynolds, double liquidReynolds, double gasFlowRate, double liquidFlowRate,
          double diameter, double roughness, double inclinationAngle, double &c0, double &ud,
          double horizontalCorrection) {
    const double inclinationSign = inclinationSignOf(inclinationAngle);
    const double flowArea = ductArea(diameter);
    ud = horizontalCorrection * inclinationSign * 0.0246 * cos(inclinationAngle) + 1.606 * pow(9.82 * surfaceTension * (liquidDensity - gasDensity) / (liquidDensity * liquidDensity), 0.25) * sin(inclinationAngle);
    c0 = 2. / (1 + pow(mixtureReynolds / 1000., 2.)) + (1.2 - 0.2 * sqrt(gasDensity / liquidDensity) * (1 - exp(-18 * voidFraction))) / (1 + pow(1000. / mixtureReynolds, 2.));
    alignDriftWithInclination(gasFlowRate, liquidFlowRate, flowArea, inclinationAngle, ud);
}

void HibikiIshii(double liquidDensity, double gasDensity, double surfaceTension, double voidFraction,
                 double mixtureReynolds, double liquidReynolds, double gasFlowRate, double liquidFlowRate,
                 double diameter, double roughness, double inclinationAngle, double &c0, double &ud,
                 double horizontalCorrection) {
    const double inclinationSign = inclinationSignOf(inclinationAngle);
    const double flowArea = ductArea(diameter);
    c0 = 1. + (1. - voidFraction) / (voidFraction + 4. * sqrt(gasDensity / liquidDensity));
    ud = (horizontalCorrection * inclinationSign * (1. - voidFraction) / (voidFraction + 4. * sqrt(gasDensity / liquidDensity))) * sqrt(9.82 * fabs(sin(inclinationAngle)) * diameter * (liquidDensity - gasDensity) * (1. - voidFraction) / (0.015 * liquidDensity));
    alignDriftWithInclination(gasFlowRate, liquidFlowRate, flowArea, inclinationAngle, ud);
}

void FrancaLahey(double liquidDensity, double gasDensity, double surfaceTension, double voidFraction,
                 double mixtureReynolds, double liquidReynolds, double gasFlowRate, double liquidFlowRate,
                 double diameter, double roughness, double inclinationAngle, double &c0, double &ud,
                 double horizontalCorrection) {
    const double inclinationSign = inclinationSignOf(inclinationAngle);
    c0 = 1.04;
    ud = horizontalCorrection * inclinationSign * 0.466;
    const double flowArea = ductArea(diameter);
    alignDriftWithInclination(gasFlowRate, liquidFlowRate, flowArea, inclinationAngle, ud);
}

namespace {

/// Blends Bhagwat Ghajar and Choi across the inclination band from 5 to 20
/// degrees, which is selector 5 in all three regimes.
///
/// @warning Two properties here are load bearing. BhagwatGhajarMod runs before
///          Choi because both write c0 and ud, so the order decides which
///          result the temporaries hold. And the interpolation keeps the form
///          ratio*x + (1 - ratio)*xMod, because rearranging it to
///          xMod + ratio*(x - xMod) is algebraically identical and rounds
///          differently.
void blendAcrossInclination(double liquidDensity, double gasDensity, double surfaceTension,
                            double voidFraction, double mixtureReynolds, double liquidReynolds,
                            double gasFlowRate, double liquidFlowRate, double diameter,
                            double roughness, double inclinationAngle, double &c0, double &ud,
                            double horizontalCorrection) {
    if (fabs(inclinationAngle) < 5 * M_PI / 180.) {
        BhagwatGhajar(liquidDensity, gasDensity, surfaceTension, voidFraction, mixtureReynolds,
                      liquidReynolds, gasFlowRate, liquidFlowRate, diameter, roughness,
                      inclinationAngle, c0, ud, horizontalCorrection);
    } else if (fabs(inclinationAngle) > 20 * M_PI / 180.) {
        Choi(liquidDensity, gasDensity, surfaceTension, voidFraction, mixtureReynolds,
             liquidReynolds, gasFlowRate, liquidFlowRate, diameter, roughness, inclinationAngle,
             c0, ud, horizontalCorrection);
    } else {
        const double blendRatio = (fabs(inclinationAngle) - 5 * M_PI / 180.) / (15 * M_PI / 180.);
        BhagwatGhajarMod(liquidDensity, gasDensity, surfaceTension, voidFraction, mixtureReynolds,
                         liquidReynolds, gasFlowRate, liquidFlowRate, diameter, roughness,
                         inclinationAngle, c0, ud, horizontalCorrection);
        const double c0BhagwatGhajarMod = c0;
        const double udBhagwatGhajarMod = ud;
        Choi(liquidDensity, gasDensity, surfaceTension, voidFraction, mixtureReynolds,
             liquidReynolds, gasFlowRate, liquidFlowRate, diameter, roughness, inclinationAngle,
             c0, ud, horizontalCorrection);
        c0 = blendRatio * c0 + (1. - blendRatio) * c0BhagwatGhajarMod;
        ud = blendRatio * ud + (1. - blendRatio) * udBhagwatGhajarMod;
    }
}

/// Applies the selectors that every flow regime accepts.
///
/// Anything outside this set falls through and leaves c0 and ud exactly as the
/// caller passed them. The original switches carried no default, and that
/// silence is observable behaviour rather than an oversight, which is also why
/// the three regimes keep separate case sets instead of sharing one table.
void applyCommonCorrelation(int correlationIndex, double liquidDensity, double gasDensity,
                            double surfaceTension, double voidFraction, double mixtureReynolds,
                            double liquidReynolds, double gasFlowRate, double liquidFlowRate,
                            double diameter, double roughness, double inclinationAngle, double &c0,
                            double &ud, double horizontalCorrection) {
    switch (correlationIndex) {
    case 0:
        Choi(liquidDensity, gasDensity, surfaceTension, voidFraction, mixtureReynolds,
             liquidReynolds, gasFlowRate, liquidFlowRate, diameter, roughness, inclinationAngle,
             c0, ud, horizontalCorrection);
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
        blendAcrossInclination(liquidDensity, gasDensity, surfaceTension, voidFraction,
                               mixtureReynolds, liquidReynolds, gasFlowRate, liquidFlowRate,
                               diameter, roughness, inclinationAngle, c0, ud, horizontalCorrection);
        break;
    }
}

} // namespace

void C0UdDisperso(double liquidDensity, double gasDensity, double surfaceTension, double voidFraction,
                  double mixtureReynolds, double liquidReynolds, double gasFlowRate, double liquidFlowRate,
                  double diameter, double roughness, double inclinationAngle, double &c0, double &ud,
                  double horizontalCorrection, int estabCol, int correlationIndex) {
    applyCommonCorrelation(correlationIndex, liquidDensity, gasDensity, surfaceTension, voidFraction,
                           mixtureReynolds, liquidReynolds, gasFlowRate, liquidFlowRate, diameter,
                           roughness, inclinationAngle, c0, ud, horizontalCorrection);
}

void C0UdAnularChurn(double liquidDensity, double gasDensity, double surfaceTension, double voidFraction,
                     double mixtureReynolds, double liquidReynolds, double gasFlowRate,
                     double liquidFlowRate, double diameter, double roughness, double inclinationAngle,
                     double &c0, double &ud, double horizontalCorrection, int estabCol, int correlationIndex) {
    // Hibiki Ishii is accepted in this regime and in no other.
    if (correlationIndex == 3) {
        HibikiIshii(liquidDensity, gasDensity, surfaceTension, voidFraction, mixtureReynolds,
                    liquidReynolds, gasFlowRate, liquidFlowRate, diameter, roughness,
                    inclinationAngle, c0, ud, horizontalCorrection);
        return;
    }
    applyCommonCorrelation(correlationIndex, liquidDensity, gasDensity, surfaceTension, voidFraction,
                           mixtureReynolds, liquidReynolds, gasFlowRate, liquidFlowRate, diameter,
                           roughness, inclinationAngle, c0, ud, horizontalCorrection);
}

void C0UdEstratificado(double liquidDensity, double gasDensity, double surfaceTension, double voidFraction,
                       double mixtureReynolds, double liquidReynolds, double gasFlowRate,
                       double liquidFlowRate, double diameter, double roughness, double inclinationAngle,
                       double &c0, double &ud, double horizontalCorrection, int estabCol,
                       int correlationIndex) {
    // Franca Lahey is accepted in this regime and in no other.
    if (correlationIndex == 2) {
        FrancaLahey(liquidDensity, gasDensity, surfaceTension, voidFraction, mixtureReynolds,
                    liquidReynolds, gasFlowRate, liquidFlowRate, diameter, roughness,
                    inclinationAngle, c0, ud, horizontalCorrection);
        return;
    }
    applyCommonCorrelation(correlationIndex, liquidDensity, gasDensity, surfaceTension, voidFraction,
                           mixtureReynolds, liquidReynolds, gasFlowRate, liquidFlowRate, diameter,
                           roughness, inclinationAngle, c0, ud, horizontalCorrection);
}

}  // namespace correlations

namespace coefficient {

/*
 * The five bodies below were MOVED, token for token, from SisProd.cpp. The
 * transformation is a substitution table -- the method signature, and eleven
 * SProd members rewritten as ClosureState fields -- and refactor-harness/
 * c0ud-move.py applies its inverse and compares the token stream against the
 * commit the move started from. 13556 tokens, exact.
 *
 * Nothing here is tidied. The variants disagree about which cell's accessory
 * the horizontal correction reads, whether the flow-pattern map runs at all,
 * whether the transition counter is kept, and which of arq.escorregaTran and
 * arq.escorregaPerm ends the calculation. Six preserved anomalies are
 * catalogued as A3-01 to A3-06 in
 * specs/001-refatoracao-sisprod/evidencia/c0ud-diff.md; the dead chain
 * betneg -> ul0 -> mult0, dead in all five, is A3-07.
 */

void instantaneous(const ClosureState &state, int ind, double &c0, double &ud) {
    int timeStep = 20;
    state.cells[ind].transic0 = state.cells[ind].transic;
    if (state.cells[ind].dt < 0.8)
        timeStep *= (0.8 / state.cells[ind].dt);
    c0 = 1.;
    ud = 0.;
    if (state.cells[ind - 1].velPig > 0 && state.cells[ind - 1].estadoPig == 1) {
        c0 = 1.;
        ud = 0.;
        state.cells[ind].arranjo = 1.;
        state.cells[ind - 1].arranjoR = 1.;
        state.cells[ind - 1].perdaEstratL = 0.;
        state.cells[ind - 1].perdaEstratG = 0.;
    } else if (state.cells[ind].velPig < 0 && state.cells[ind].estadoPig == 1) {
        c0 = 1.;
        ud = 0.;
        state.cells[ind].arranjo = 1.;
        state.cells[ind - 1].arranjoR = 1.;
        state.cells[ind - 1].perdaEstratL = 0.;
        state.cells[ind - 1].perdaEstratG = 0.;
    } else if ((state.cells[ind].acsr.tipo != 4 || state.cells[ind].acsr.bcs.freqnova <= 1.)) {
        double hns;
        double razdx = state.cells[ind].dxL / (state.cells[ind].dx + state.cells[ind].dxL);
        double razdx0;
        if (ind > 0)
            razdx0 = state.cells[ind - 1].dxL / (state.cells[ind - 1].dx + state.cells[ind - 1].dxL);
        else
            razdx0 = razdx;

        if (ind > 0)
            hns = 1. - state.cells[ind - 1].alfPigD;
        else
            hns = 1. - state.cells[ind].alf;
        if (state.cells[ind].QG < 0)
            hns = 1. - state.cells[ind].alfPigE;
        if (fabs(state.cells[ind].QG) < (*state.globals).localtiny * 1e-5) {
            if (fabs(state.cells[ind].alfPigE) < (*state.globals).localtiny && fabs(1. - state.cells[ind].alfPigE) < (*state.globals).localtiny && fabs(state.cells[ind - 1].alfPigD) > (*state.globals).localtiny && fabs(1. - state.cells[ind - 1].alfPigD) > (*state.globals).localtiny)
                hns = 1. - state.cells[ind - 1].alfPigD;
            else if (fabs(state.cells[ind - 1].alfPigD) < (*state.globals).localtiny && fabs(1. - state.cells[ind - 1].alfPigD) < (*state.globals).localtiny && fabs(state.cells[ind].alfPigE) > (*state.globals).localtiny && fabs(1. - state.cells[ind].alfPigE) > (*state.globals).localtiny)
                hns = 1. - state.cells[ind].alfPigE;
            else
                hns = 0.5 * (1. - state.cells[ind].alfPigE + 1. - state.cells[ind - 1].alfPigD);
        }
        if (hns < (*state.globals).localtiny || hns > 1. - (*state.globals).localtiny)
            hns = 0.5 * (1. - state.cells[ind].alfPigE + 1. - state.cells[ind - 1].alfPigD);

        double hol0 = hns;
        double alf0 = 1 - hol0;
        double alf1;
        alf1 = state.cells[ind].alf;

        double alfneg;
        if (ind > 1)
            alfneg = state.cells[ind - 2].alf;
        else if (ind > 0)
            alfneg = state.cells[ind - 1].alf;
        else
            alfneg = state.cells[ind].alf;

        double betI = state.cells[ind].betL;
        if (ind > 0)
            betI = state.cells[ind - 1].betPigD;
        if (((0. * state.cells[ind].QG + 1 * state.cells[ind].QL) < 0.))
            betI = state.cells[ind].betPigE; // duvidabeta

        double betneg;
        if (ind > 0) {
            betneg = state.cells[ind - 1].betL;
            if (ind > 1)
                betneg = state.cells[ind - 2].betPigD;
            if ((0.99 * state.cells[ind - 1].QG + 0.01 * state.cells[ind - 1].QL) < 0.)
                betneg = state.cells[ind - 1].betPigE; // duvidabeta

        } else
            betneg = state.cells[ind].bet;

        double pmed;
        double pmed0 = 0.;

        pmed = state.cells[ind].presaux;
        if (ind > 0)
            pmed0 = state.cells[ind - 1].presaux;
        if (ind == state.lastCell)
            pmed = state.cells[ind].pres;
        else
            pmed0 = state.cells[ind].presaux;
        double tmed = razdx * state.cells[ind].temp + (1 - razdx) * state.cells[ind].tempL;
        tmed = state.cells[ind].tempL;
        if (state.cells[ind].VTemper < 0.)
            tmed = state.cells[ind].temp;
        double tmed0;
        if (ind > 0)
            tmed0 = razdx0 * state.cells[ind - 1].temp + (1 - razdx0) * state.cells[ind - 1].tempL;
        else
            tmed0 = tmed;
        if (ind < state.lastCell)
            tmed = state.cells[ind].temp;
        else
            tmed = state.gasSurfaceTemperature;

        double correcHor = 1.;
        if (fabs(state.cells[ind].duto.teta) < 1e-10) {
            if (state.cells[ind - 1].acsr.tipo != 5 || state.cells[ind - 1].acsr.chk.AreaGarg > 1e-10) {
                if (state.cells[ind].angEsq < 0 && state.cells[ind].angDir < 0)
                    correcHor = -1.;
                else if (state.cells[ind].angEsq > 0 && state.cells[ind].angDir > 0)
                    correcHor = 1.;
            } else {
                if (state.cells[ind].angDir < 0)
                    correcHor = -1.;
                else if (state.cells[ind].angDir > 0)
                    correcHor = 1.;
            }
        }

        double rlm;
        double viscl1;
        double tensup1;
        if (state.cells[ind].QL < 0.) { // testeBeta
            if (ind == 0 || ind == state.lastCell)
                rlm = (1 - betI) * state.cells[ind].flui.MasEspLiq(pmed, tmed) + betI * state.cells[ind].fluicol.MasEspFlu(pmed, tmed);
            else
                rlm = (1 - betI) * state.cells[ind].rpCi + betI * state.cells[ind].rcCi;
            viscl1 = (1 - betI) * state.cells[ind].flui.ViscOleo(pmed, tmed) + betI * state.cells[ind].fluicol.VisFlu(pmed, tmed);
            tensup1 = (1 - betI) * state.cells[ind].flui.TensSuper(pmed, tmed) + betI * state.cells[ind].fluicol.TensSuper(pmed, tmed);
        } else {
            if (ind == 0 || ind == state.lastCell)
                rlm = (1 - betI) * state.cells[ind - 1].flui.MasEspLiq(pmed, tmed) + betI * state.cells[ind - 1].fluicol.MasEspFlu(pmed, tmed);
            else
                rlm = (1 - betI) * state.cells[ind].rpCi + betI * state.cells[ind - 1].rcCi;
            viscl1 = (1 - betI) * state.cells[ind - 1].flui.ViscOleo(pmed, tmed) + betI * state.cells[ind - 1].fluicol.VisFlu(pmed, tmed);
            tensup1 = (1 - betI) * state.cells[ind - 1].flui.TensSuper(pmed, tmed) + betI * state.cells[ind - 1].fluicol.TensSuper(pmed, tmed);
        }

        double rgm;
        double viscg1;
        if (state.cells[ind].QG < 0.) {
            if (ind == 0 || ind == state.lastCell)
                rgm = state.cells[ind].flui.MasEspGas(pmed, tmed);
            else
                rgm = state.cells[ind].rgCi;
            viscg1 = state.cells[ind].flui.ViscGas(pmed, tmed);
        } else {
            if (ind == 0 || ind == state.lastCell)
                rgm = state.cells[ind - 1].flui.MasEspGas(pmed, tmed);
            else
                rgm = state.cells[ind].rgCi;
            viscg1 = state.cells[ind - 1].flui.ViscGas(pmed, tmed);
        }

        double ug1 = (state.cells[ind].MC - state.cells[ind].Mliqini) / rgm;
        double ul1 = state.cells[ind].Mliqini / rlm;
        double dia1 = state.cells[ind].duto.a;
        if (ind > 0 && ug1 >= 0)
            dia1 = state.cells[ind - 1].duto.a;

        double A1 = M_PI * dia1 * dia1 / 4.;

        double rmed = hns * rlm + (1 - hns) * rgm;
        double visc = (hns * viscl1 + (1 - hns) * viscg1) / pow(10., 3.);
        double nrey = dia1 * rmed * (fabs(ug1) / A1 + fabs(ul1) / A1) / visc;
        double nreyl = dia1 * rlm * (fabs(ug1) / A1 + fabs(ul1) / A1) / (viscl1 / 1000.);

        int xarr1 = 1;
        double dtot = state.cells[ind].dxL + state.cells[ind].dx;
        double razL = state.cells[ind].dx;
        double raz = state.cells[ind].dxL;
        double ang = (razL * state.cells[ind].dutoL.teta + raz * state.cells[ind].duto.teta) / dtot;
        if (ind >= 2) {
            if (state.cells[ind - 2].acsr.tipo == 5 && state.cells[ind - 2].acsr.chk.AreaGarg <= (1e-3)) {
                if (state.cells[ind].QG >= 0)
                    ang = state.cells[ind].duto.teta;
                else
                    ang = state.cells[ind].dutoR.teta;
            } else {
                if (state.cells[ind].QG >= 0)
                    ang = state.cells[ind].dutoL.teta;
                else
                    ang = state.cells[ind].duto.teta;
            }
        }

        double atenua = 20;

        if (nrey > 1e-30) {
            if (fabs(0 * ang + 1 * state.cells[ind].duto.teta) < 45. * M_PI / 180. && hol0 < 0.99 && hol0 > 0.01 && ind < state.lastCell - 1) {
                double ug0;
                double ul0;

                ug1 = (state.cells[ind].MC - state.cells[ind].Mliqini) / rgm;
                ul1 = state.cells[ind].Mliqini / rlm;

                ug0 = (state.cells[ind - 1].MC - state.cells[ind - 1].Mliqini) / state.cells[ind].rgLi;
                // ul0 = state.cells[ind - 1].Mliqini
                //  / ((1 - betneg) * state.cells[ind].flui.MasEspLiq(pmed0, tmed0)
                ul0 = state.cells[ind - 1].Mliqini / ((1 - betneg) * state.cells[ind].rpLi + betneg * state.cells[ind].rcLi);

                estratificado testamapa(dia1, ul1, ug1, rlm, rgm, viscl1 / pow(10., 3.), viscg1 / pow(10., 3.), hol0,
                                        state.cells[ind].duto.teta, state.cells[ind].duto.rug / dia1);

                if (state.selectors.stratified == 2)
                    testamapa.mapaTD();
                else
                    testamapa.mapaTD(1);

                xarr1 = testamapa.arr;

                if (xarr1 == -1) {
                    if (state.cells[ind].arranjo != 0) {
                        if (((state.cells[ind].arranjo != xarr1) || state.cells[ind].transic > 0)) {
                            if ((state.cells[ind].arranjo != xarr1) && state.cells[ind].transic > 0)
                                state.cells[ind].transic = 0;
                            state.cells[ind].transic++;
                            if (state.cells[ind].transic > atenua - 1)
                                state.cells[ind].transic = 0;
                        } else
                            state.cells[ind].transic = 0;
                    }
                    state.cells[ind].arranjo = xarr1 = testamapa.arr;
                    state.cells[ind - 1].arranjoR = testamapa.arr;
                    state.cells[ind - 1].perdaEstratL = testamapa.fatorperdaLiq;
                    state.cells[ind - 1].perdaEstratG = testamapa.fatorperdaGas;
                    double c0D;
                    double udD;
                    double c0E;
                    double udE;

                    driftflux::correlations::C0UdDisperso(rlm, rgm, tensup1, alf0, nrey, nreyl, ug1, ul1, dia1, state.cells[ind].duto.rug, ang,
                                 c0D, udD, correcHor, state.cells[ind].estabCol, state.selectors.dispersed);
                    driftflux::correlations::C0UdEstratificado(rlm, rgm, tensup1, alf0, nrey, nreyl, ug1, ul1, dia1, state.cells[ind].duto.rug, ang,
                                      c0E, udE, correcHor, state.cells[ind].estabCol, state.selectors.stratified);

                    double mult0, mult1;
                    mult0 = 1.;
                    if (ul0 < 0.)
                        mult0 = 0.;
                    mult1 = 0.;
                    if (ul1 < 0.)
                        mult1 = 1.;
                    double alf0E = state.cells[ind - 1].alf;

                    double jmax = 0.05;
                    double jmin = 0.005;
                    if ((fabs(ug1) + fabs(ul1)) / A1 > jmax) {
                        c0 = c0E;
                        ud = udE;
                    } else if ((fabs(ug1) + fabs(ul1)) / A1 < jmin) {
                        c0 = c0D;
                        ud = udD;
                    } else {
                        double raz = (jmax - (fabs(ug1) + fabs(ul1)) / A1) / (jmax - jmin);
                        c0 = ((1. - raz) * c0E + raz * c0D);
                        ud = ((1. - raz) * udE + raz * udD);
                    }

                    if (state.cells[ind].transic > 0) {
                        c0 = (c0 * state.cells[ind].transic + state.cells[ind].c0 * (atenua - state.cells[ind].transic)) / atenua;
                        ud = (ud * state.cells[ind].transic + state.cells[ind].ud * (atenua - state.cells[ind].transic)) / atenua;
                    }
                }
            }
            if (xarr1 == 1) {

                arranjo testamapa2(dia1, ul1 / A1, ug1 / A1, rlm, rgm, viscl1 / pow(10., 3.), viscg1 / pow(10., 3.), hol0,
                                   state.cells[ind].duto.teta, tensup1, state.input.mapaArranjo, state.globals);
                xarr1 = testamapa2.verificaArr();

                driftflux::correlations::C0UdDisperso(rlm, rgm, tensup1, alf0, nrey, nreyl, ug1, ul1, dia1, state.cells[ind].duto.rug, ang,
                             c0, ud, correcHor, state.cells[ind].estabCol, state.selectors.dispersed);

                if (xarr1 == -2) {
                    driftflux::correlations::C0UdAnularChurn(rlm, rgm, tensup1, alf0, nrey, nreyl, ug1, ul1, dia1, state.cells[ind].duto.rug, ang,
                                    c0, ud, correcHor, state.cells[ind].estabCol, state.selectors.annularChurn);
                }
                if (fabs(ug1 / state.cells[ind].duto.area) > 5. && alf0 >= 0.75) {
                    atenua = 20;
                    if (state.selectors.annularChurn == 3 && state.selectors.dispersed == 1)
                        atenua = 200;
                }

                if (state.cells[ind].arranjo != 0) {
                    if ((xarr1 != state.cells[ind].arranjo || state.cells[ind].transic > 0)) {
                        if (xarr1 != state.cells[ind].arranjo && state.cells[ind].transic > 0)
                            state.cells[ind].transic = 0;
                        state.cells[ind].transic++;
                        if (state.cells[ind].transic > atenua - 1)
                            state.cells[ind].transic = 0;
                    } else
                        state.cells[ind].transic = 0;
                }
                if (state.cells[ind].transic > 0) {
                    c0 = (c0 * state.cells[ind].transic + state.cells[ind].c0 * (atenua - state.cells[ind].transic)) / atenua;
                    ud = (ud * state.cells[ind].transic + state.cells[ind].ud * (atenua - state.cells[ind].transic)) / atenua;
                }
                state.cells[ind].arranjo = xarr1;
                state.cells[ind - 1].arranjoR = xarr1;
            }
            state.cells[ind].c0Spare = c0;
            state.cells[ind].udSpare = ud;
        }
    }
    if (state.input.escorregaTran == 0) {
        double ulsmed = state.cells[ind].QL / state.cells[ind].duto.area;
        double correcaoUd = 1 - (ulsmed - 0.15) / 0.35;
        double correcaoCo = c0 - (c0 - 1) * (ulsmed - 0.15) / 0.35;
        if (correcaoUd > 1.)
            correcaoUd = 1.;
        if (correcaoUd < 0.)
            correcaoUd = 0.;
        if (correcaoCo > c0)
            correcaoCo = c0;
        if (correcaoCo < 1)
            correcaoCo = 1;
        c0 = 1. + 0 * correcaoCo;
        ud = 0. + 0. * ud * correcaoUd;
    }
}

void buffered(const ClosureState &state, int ind, double &c0, double &ud) {
    int timeStep = 20;
    if (state.cells[ind].dt < 0.8)
        timeStep *= (0.8 / state.cells[ind].dt);
    c0 = 1.;
    ud = 0.;
    if (state.cells[ind - 1].velPig > 0 && state.cells[ind - 1].estadoPig == 1) {
        c0 = 1.;
        ud = 0.;
        state.cells[ind].arranjo = 1.;
        state.cells[ind - 1].arranjoR = 1.;
        state.cells[ind - 1].perdaEstratL = 0.;
        state.cells[ind - 1].perdaEstratG = 0.;
    } else if (state.cells[ind].velPig < 0 && state.cells[ind].estadoPig == 1) {
        c0 = 1.;
        ud = 0.;
        state.cells[ind].arranjo = 1.;
        state.cells[ind - 1].arranjoR = 1.;
        state.cells[ind - 1].perdaEstratL = 0.;
        state.cells[ind - 1].perdaEstratG = 0.;
    } else if ((state.cells[ind].acsr.tipo != 4 || state.cells[ind].acsr.bcs.freqnova <= 1.)) {
        double hns;
        double razdx = state.cells[ind].dxL / (state.cells[ind].dx + state.cells[ind].dxL);
        double razdx0;
        if (ind > 0)
            razdx0 = state.cells[ind - 1].dxL / (state.cells[ind - 1].dx + state.cells[ind - 1].dxL);
        else
            razdx0 = razdx;

        if (ind > 0)
            hns = 1. - state.cells[ind - 1].alfPigD;
        else
            hns = 1. - state.cells[ind].alf;
        if ((state.cells[ind].MCBuf - state.cells[ind].MliqiniBuf) < 0)
            hns = 1. - state.cells[ind].alfPigE;
        if (fabs((state.cells[ind].MCBuf - state.cells[ind].MliqiniBuf)) < (*state.globals).localtiny * 1e-5) {
            if (fabs(state.cells[ind].alfPigE) < (*state.globals).localtiny && fabs(1. - state.cells[ind].alfPigE) < (*state.globals).localtiny && fabs(state.cells[ind - 1].alfPigD) > (*state.globals).localtiny && fabs(1. - state.cells[ind - 1].alfPigD) > (*state.globals).localtiny)
                hns = 1. - state.cells[ind - 1].alfPigD;
            else if (fabs(state.cells[ind - 1].alfPigD) < (*state.globals).localtiny && fabs(1. - state.cells[ind - 1].alfPigD) < (*state.globals).localtiny && fabs(state.cells[ind].alfPigE) > (*state.globals).localtiny && fabs(1. - state.cells[ind].alfPigE) > (*state.globals).localtiny)
                hns = 1. - state.cells[ind].alfPigE;
            else
                hns = 0.5 * (1. - state.cells[ind].alfPigE + 1. - state.cells[ind - 1].alfPigD);
        }
        if (hns < (*state.globals).localtiny || hns > 1. - (*state.globals).localtiny)
            hns = 0.5 * (1. - state.cells[ind].alfPigE + 1. - state.cells[ind - 1].alfPigD);

        double hol0 = hns;
        double alf0 = 1 - hol0;
        double alf1;
        alf1 = state.cells[ind].alf;

        double alfneg;
        if (ind > 1)
            alfneg = state.cells[ind - 2].alf;
        else if (ind > 0)
            alfneg = state.cells[ind - 1].alf;
        else
            alfneg = state.cells[ind].alf;

        double betI = state.cells[ind].betL;
        if (ind > 0)
            betI = state.cells[ind - 1].betPigD;
        if ((state.cells[ind].MliqiniBuf) < 0.)
            betI = state.cells[ind].betPigE; // testeBeta
        betI = state.cells[ind].betPigE;     // duvidabeta
        double betneg;
        if (ind > 0) {
            betneg = state.cells[ind - 1].betL;
            if (ind > 1)
                betneg = state.cells[ind - 2].betPigD;
            if (state.cells[ind].MliqiniLBuf < 0.)
                betneg = state.cells[ind - 1].betPigE; // testeBeta
            betneg = state.cells[ind - 1].betPigE;     // duvidabeta
        } else
            betneg = state.cells[ind].bet;

        double pmed;
        double pmed0 = 0.;

        pmed = razdx * state.cells[ind].presBuf + (1 - razdx) * state.cells[ind].presLBuf;
        if (ind == state.lastCell)
            pmed = state.cells[ind].presBuf;
        pmed0 = state.cells[ind].presauxL;
        double tmed = razdx * state.cells[ind].temp + (1 - razdx) * state.cells[ind].tempL;
        tmed = state.cells[ind].tempL;
        if (state.cells[ind].VTemper < 0.) {
            if (ind < state.lastCell)
                tmed = state.cells[ind].temp;
            else
                tmed = state.gasSurfaceTemperature;
        }
        double tmed0;
        if (ind > 0)
            tmed0 = razdx0 * state.cells[ind - 1].temp + (1 - razdx0) * state.cells[ind - 1].tempL;
        else
            tmed0 = tmed;

        double correcHor = 1.;
        if (fabs(state.cells[ind].duto.teta) < 1e-10) {
            if (state.cells[ind].acsr.tipo != 5 || state.cells[ind].acsr.chk.AreaGarg > 1e-10) {
                if (state.cells[ind].angEsq < 0 && state.cells[ind].angDir < 0)
                    correcHor = -1.;
                else if (state.cells[ind].angEsq > 0 && state.cells[ind].angDir > 0)
                    correcHor = 1.;
            } else {
                if (state.cells[ind].angDir < 0)
                    correcHor = -1.;
                else if (state.cells[ind].angDir > 0)
                    correcHor = 1.;
            }
        }

        double rlm;
        double viscl1;
        double tensup1;
        if ((state.cells[ind].MliqiniBuf) < 0.) { // testeBeta
            rlm = (1 - betI) * state.cells[ind].flui.MasEspLiq(pmed, tmed) + betI * state.cells[ind].fluicol.MasEspFlu(pmed, tmed);
            viscl1 = (1 - betI) * state.cells[ind].flui.ViscOleo(pmed, tmed) + betI * state.cells[ind].fluicol.VisFlu(pmed, tmed);
            tensup1 = (1 - betI) * state.cells[ind].flui.TensSuper(pmed, tmed) + betI * state.cells[ind].fluicol.TensSuper(pmed, tmed);
        } else {
            rlm = (1 - betI) * state.cells[ind - 1].flui.MasEspLiq(pmed, tmed) + betI * state.cells[ind - 1].fluicol.MasEspFlu(pmed, tmed);
            viscl1 = (1 - betI) * state.cells[ind - 1].flui.ViscOleo(pmed, tmed) + betI * state.cells[ind - 1].fluicol.VisFlu(pmed, tmed);
            tensup1 = (1 - betI) * state.cells[ind - 1].flui.TensSuper(pmed, tmed) + betI * state.cells[ind - 1].fluicol.TensSuper(pmed, tmed);
        }

        double rgm;
        double viscg1;
        if ((state.cells[ind].MCBuf - state.cells[ind].MliqiniBuf) < 0.) {
            rgm = state.cells[ind].flui.MasEspGas(pmed, tmed);
            viscg1 = state.cells[ind].flui.ViscGas(pmed, tmed);
        } else {
            rgm = state.cells[ind - 1].flui.MasEspGas(pmed, tmed);
            viscg1 = state.cells[ind - 1].flui.ViscGas(pmed, tmed);
        }

        double ug1 = (state.cells[ind].MCBuf - state.cells[ind].MliqiniBuf) / rgm;
        double ul1 = (state.cells[ind].MliqiniBuf) / rlm;
        double dia1 = state.cells[ind].duto.a;
        if (ind > 0 && ug1 >= 0)
            dia1 = state.cells[ind - 1].duto.a;
        double A1 = M_PI * dia1 * dia1 / 4.;

        double rmed = hns * rlm + (1 - hns) * rgm;
        double visc = (hns * viscl1 + (1 - hns) * viscg1) / pow(10., 3.);
        double nrey = dia1 * rmed * (fabs(ug1) / A1 + fabs(ul1) / A1) / visc;
        double nreyl = dia1 * rlm * (fabs(ug1) / A1 + fabs(ul1) / A1) / (viscl1 / 1000.);

        int xarr1 = 1;
        double dtot = state.cells[ind].dxL + state.cells[ind].dx;
        double razL = state.cells[ind].dxL;
        double raz = state.cells[ind].dx;
        double ang = (raz * state.cells[ind].dutoL.teta + razL * state.cells[ind].duto.teta) / dtot;
        if (ind >= 2) {
            if (state.cells[ind - 2].acsr.tipo == 5 && state.cells[ind - 2].acsr.chk.AreaGarg <= (1e-3)) {
                if ((state.cells[ind].MCBuf - state.cells[ind].MliqiniBuf) >= 0)
                    ang = state.cells[ind].duto.teta;
                else
                    ang = state.cells[ind].dutoR.teta;
            } else {
                if ((state.cells[ind].MCBuf - state.cells[ind].MliqiniBuf) >= 0)
                    ang = state.cells[ind].dutoL.teta;
                else
                    ang = state.cells[ind].duto.teta;
            }
        }
        double atenua = 20;
        if (nrey > 1e-30) {
            if (fabs(0 * ang + 1 * state.cells[ind].duto.teta) < 45. * M_PI / 180. && hol0 < 0.99 && hol0 > 0.01 && ind < state.lastCell - 1) {
                double ug0;
                double ul0;

                ug1 = (state.cells[ind].MCBuf - state.cells[ind].MliqiniBuf) / rgm;
                ul1 = (state.cells[ind].MliqiniBuf) / rlm;

                ug0 = (state.cells[ind - 1].MCBuf - state.cells[ind - 1].MliqiniBuf) / state.cells[ind].flui.MasEspGas(pmed0, tmed0);
                ul0 = (state.cells[ind - 1].MliqiniBuf) / ((1 - betneg) * state.cells[ind].flui.MasEspLiq(pmed0, tmed0) + betneg * state.cells[ind].fluicol.MasEspFlu(pmed0, tmed0));

                xarr1 = state.cells[ind].arranjo;
                if (xarr1 == -1) {

                    double c0D;
                    double udD;
                    double c0E;
                    double udE;
                    driftflux::correlations::C0UdDisperso(rlm, rgm, tensup1, alf0, nrey, nreyl, ug1, ul1, dia1, state.cells[ind].duto.rug, ang,
                                 c0D, udD, correcHor, state.cells[ind].estabCol, state.selectors.dispersed);
                    driftflux::correlations::C0UdEstratificado(rlm, rgm, tensup1, alf0, nrey, nreyl, ug1, ul1, dia1, state.cells[ind].duto.rug, ang,
                                      c0E, udE, correcHor, state.cells[ind].estabCol, state.selectors.stratified);

                    double mult0, mult1;
                    mult0 = 1.;
                    if (ul0 < 0.)
                        mult0 = 0.;
                    mult1 = 0.;
                    if (ul1 < 0.)
                        mult1 = 1.;
                    double alf0E = state.cells[ind - 1].alf;

                    double jmax = 0.05;
                    double jmin = 0.005;
                    if ((fabs(ug1) + fabs(ul1)) / A1 > jmax) {
                        c0 = c0E;
                        ud = udE;
                    } else if ((fabs(ug1) + fabs(ul1)) / A1 < jmin) {
                        c0 = c0D;
                        ud = udD;
                    } else {
                        double raz = (jmax - (fabs(ug1) + fabs(ul1)) / A1) / (jmax - jmin);
                        c0 = ((1. - raz) * c0E + raz * c0D);
                        ud = ((1. - raz) * udE + raz * udD);
                    }

                    if (state.cells[ind].transic > 0) {
                        c0 = (c0 * state.cells[ind].transic + state.cells[ind].c0 * (atenua - state.cells[ind].transic)) / atenua;
                        ud = (ud * state.cells[ind].transic + state.cells[ind].ud * (atenua - state.cells[ind].transic)) / atenua;
                    }
                }
            }
            if (xarr1 != -1) {

                driftflux::correlations::C0UdDisperso(rlm, rgm, tensup1, alf0, nrey, nreyl, ug1, ul1, dia1, state.cells[ind].duto.rug, ang,
                             c0, ud, correcHor, state.cells[ind].estabCol, state.selectors.dispersed);
                if (xarr1 == -2) {
                    driftflux::correlations::C0UdAnularChurn(rlm, rgm, tensup1, alf0, nrey, nreyl, ug1, ul1, dia1, state.cells[ind].duto.rug, ang,
                                    c0, ud, correcHor, state.cells[ind].estabCol, state.selectors.annularChurn);
                }
                if (fabs(ug1 / state.cells[ind].duto.area) > 5. && alf0 >= 0.75) {
                    atenua = 20;
                    if (state.selectors.annularChurn == 3 && state.selectors.dispersed == 1)
                        atenua = 200;
                }

                if (state.cells[ind].transic > 0) {
                    c0 = (c0 * state.cells[ind].transic + state.cells[ind].c0 * (atenua - state.cells[ind].transic)) / atenua;
                    ud = (ud * state.cells[ind].transic + state.cells[ind].ud * (atenua - state.cells[ind].transic)) / atenua;
                }
                state.cells[ind].arranjo = xarr1;
                state.cells[ind - 1].arranjoR = xarr1;
            }
            state.cells[ind].c0Spare = c0;
            state.cells[ind].udSpare = ud;
        }
    }
    if (state.input.escorregaTran == 0) {
        double ulsmed = state.cells[ind].QL / state.cells[ind].duto.area;
        double correcaoUd = 1 - (ulsmed - 0.15) / 0.35;
        double correcaoCo = c0 - (c0 - 1) * (ulsmed - 0.15) / 0.35;
        if (correcaoUd > 1.)
            correcaoUd = 1.;
        if (correcaoUd < 0.)
            correcaoUd = 0.;
        if (correcaoCo > c0)
            correcaoCo = c0;
        if (correcaoCo < 1)
            correcaoCo = 1;
        c0 = 1. + 0 * correcaoCo;
        ud = 0. + 0 * ud * correcaoUd;
    }
}

void initialization(const ClosureState &state, int ind, double &c0, double &ud) {
    int timeStep = 20;
    if (state.cells[ind].dt < 0.8)
        timeStep *= (0.8 / state.cells[ind].dt);
    c0 = 1.;
    ud = 0.;
    if (state.cells[ind].velPig < 0 && state.cells[ind].estadoPig == 1) {
        c0 = 1.;
        ud = 0.;
        state.cells[ind].arranjo = 1.;

    } else if ((state.cells[ind].acsr.tipo != 4 || state.cells[ind].acsr.bcs.freqnova <= 1.)) {
        double hns;

        hns = 1 - state.inletVoidFraction;

        if (state.cells[ind].QG < 0)
            hns = 1. - state.cells[ind].alfPigE;
        if (fabs(state.cells[ind].QG) < (*state.globals).localtiny * 1e-5) {
            if (fabs(state.cells[ind].alfPigE) < (*state.globals).localtiny && fabs(1. - state.cells[ind].alfPigE) < (*state.globals).localtiny && fabs(state.inletVoidFraction) > (*state.globals).localtiny && fabs(1. - state.inletVoidFraction) > (*state.globals).localtiny)
                hns = 1 - state.inletVoidFraction;
            else if (fabs(state.inletVoidFraction) < (*state.globals).localtiny && fabs(1. - state.cells[ind - 1].alfPigD) < (*state.globals).localtiny && fabs(state.cells[ind].alfPigE) > (*state.globals).localtiny && fabs(1. - state.cells[ind].alfPigE) > (*state.globals).localtiny)
                hns = 1. - state.cells[ind].alfPigE;
            else
                hns = 0.5 * (1. - state.cells[ind].alfPigE + 1. - state.inletVoidFraction);
        }
        if (hns < (*state.globals).localtiny || hns > 1. - (*state.globals).localtiny)
            hns = 0.5 * (1. - state.cells[ind].alfPigE + 1. - state.inletVoidFraction);

        double hol0 = hns;
        double alf0 = 1 - hol0;
        double alf1;
        alf1 = state.cells[ind].alf;

        double alfneg;
        if (ind > 1)
            alfneg = state.inletVoidFraction;
        else if (ind > 0)
            alfneg = state.inletVoidFraction;
        else
            alfneg = state.cells[ind].alf;

        double betI = state.cells[ind].betL;
        if (ind > 0)
            betI = state.inletColumnFraction;
        if (state.cells[ind].QL < 0.)
            betI = state.cells[ind].betPigE; // testeBeta
        betI = state.cells[ind].betPigE;     // duvidabeta
        double betneg;
        if (ind > 0) {
            betneg = state.inletColumnFraction;

        } else
            betneg = state.cells[ind].bet;

        double pmed;
        double pmed0 = 0.;

        pmed = state.inletPressure;
        if (ind > 0)
            pmed0 = state.inletPressure;
        if (ind == state.lastCell)
            pmed = state.cells[ind].pres;
        else
            pmed0 = state.inletPressure;
        double tmed = state.inletTemperature;

        double tmed0;
        tmed0 = tmed;

        double correcHor = 1.;
        if (fabs(state.cells[ind].duto.teta) < 1e-10) {
            if (state.cells[ind].acsr.tipo != 5 || state.cells[ind].acsr.chk.AreaGarg > 1e-10) {
                if (state.cells[ind].angEsq < 0 && state.cells[ind].angDir < 0)
                    correcHor = -1.;
                else if (state.cells[ind].angEsq > 0 && state.cells[ind].angDir > 0)
                    correcHor = 1.;
            } else {
                if (state.cells[ind].angDir < 0)
                    correcHor = -1.;
                else if (state.cells[ind].angDir > 0)
                    correcHor = 1.;
            }
        }

        double rlm;
        double viscl1;
        double tensup1;
        if (state.cells[ind].QL < 0.) { // testeBeta
            rlm = (1 - betI) * state.cells[ind].flui.MasEspLiq(pmed, tmed) + betI * state.cells[ind].fluicol.MasEspFlu(pmed, tmed);
            viscl1 = (1 - betI) * state.cells[ind].flui.ViscOleo(pmed, tmed) + betI * state.cells[ind].fluicol.VisFlu(pmed, tmed);
            tensup1 = (1 - betI) * state.cells[ind].flui.TensSuper(pmed, tmed) + betI * state.cells[ind].fluicol.TensSuper(pmed, tmed);
        } else {
            rlm = (1 - betI) * (*state.cells[ind].fluiL).MasEspLiq(pmed, tmed) + betI * state.cells[ind].fluicol.MasEspFlu(pmed, tmed);
            viscl1 = (1 - betI) * (*state.cells[ind].fluiL).ViscOleo(pmed, tmed) + betI * state.cells[ind].fluicol.VisFlu(pmed, tmed);
            tensup1 = (1 - betI) * (*state.cells[ind].fluiL).TensSuper(pmed, tmed) + betI * state.cells[ind].fluicol.TensSuper(pmed, tmed);
        }

        double rgm;
        double viscg1;
        if (state.cells[ind].QG < 0.) {
            rgm = state.cells[ind].flui.MasEspGas(pmed, tmed);
            viscg1 = state.cells[ind].flui.ViscGas(pmed, tmed);
        } else {
            rgm = (*state.cells[ind].fluiL).MasEspGas(pmed, tmed);
            viscg1 = (*state.cells[ind].fluiL).ViscGas(pmed, tmed);
        }

        double ug1 = (state.cells[ind].MC - state.cells[ind].Mliqini) / rgm;
        double ul1 = state.cells[ind].Mliqini / rlm;
        double dia1 = state.cells[ind].duto.a;
        if (ind > 0 && ug1 >= 0)
            dia1 = state.cells[ind - 1].duto.a;
        double A1 = M_PI * dia1 * dia1 / 4.;

        double rmed = hns * rlm + (1 - hns) * rgm;
        double visc = (hns * viscl1 + (1 - hns) * viscg1) / pow(10., 3.);
        double nrey = dia1 * rmed * (fabs(ug1) / A1 + fabs(ul1) / A1) / visc;
        double nreyl = dia1 * rlm * (fabs(ug1) / A1 + fabs(ul1) / A1) / (viscl1 / 1000.);

        int xarr1 = 1;
        double dtot = state.cells[ind].dxL + state.cells[ind].dx;
        double razL = state.cells[ind].dxL;
        double raz = state.cells[ind].dx;
        double ang = (raz * state.cells[ind].dutoL.teta + razL * state.cells[ind].duto.teta) / dtot;
        double atenua = 20.;
        if (nrey > 1e-30) {
            if (fabs(0 * ang + 1 * state.cells[ind].duto.teta) < 45. * M_PI / 180. && hol0 < 0.99 && hol0 > 0.01 && ind < state.lastCell - 1) {
                double ug0;
                double ul0;

                ug1 = (state.cells[ind].MC - state.cells[ind].Mliqini) / rgm;
                ul1 = state.cells[ind].Mliqini / rlm;

                ug0 = (state.cells[ind].MC - state.cells[ind].Mliqini) / (*state.cells[ind].fluiL).MasEspGas(pmed0, tmed0);
                ul0 = state.cells[ind].Mliqini / ((1 - betneg) * (*state.cells[ind].fluiL).MasEspLiq(pmed0, tmed0) + betneg * state.cells[ind].fluicol.MasEspFlu(pmed0, tmed0));

                estratificado testamapa(dia1, ul1, ug1, rlm, rgm, viscl1 / pow(10., 3.), viscg1 / pow(10., 3.), hol0,
                                        state.cells[ind].duto.teta, state.cells[ind].duto.rug / dia1);

                testamapa.mapaTD();
                xarr1 = testamapa.arr;
                if (xarr1 == -1) {
                    if (state.cells[ind].arranjo != 0) {
                        if (((state.cells[ind].arranjo != xarr1) || state.cells[ind].transic > 0)) {
                            if ((state.cells[ind].arranjo != xarr1) && state.cells[ind].transic > 0)
                                state.cells[ind].transic = 0;
                            state.cells[ind].transic++;
                            if (state.cells[ind].transic > 19)
                                state.cells[ind].transic = 0;
                        }
                    } else
                        state.cells[ind].transic = 0;
                    state.cells[ind].arranjo = xarr1 = testamapa.arr;
                    double c0D;
                    double udD;
                    double c0E;
                    double udE;
                    driftflux::correlations::C0UdDisperso(rlm, rgm, tensup1, alf0, nrey, nreyl, ug1, ul1, dia1, state.cells[ind].duto.rug, ang,
                                 c0D, udD, correcHor, state.cells[ind].estabCol, state.selectors.dispersed);
                    driftflux::correlations::C0UdEstratificado(rlm, rgm, tensup1, alf0, nrey, nreyl, ug1, ul1, dia1, state.cells[ind].duto.rug, ang,
                                      c0E, udE, correcHor, state.cells[ind].estabCol, state.selectors.stratified);
                    double mult0, mult1;
                    mult0 = 1.;
                    if (ul0 < 0.)
                        mult0 = 0.;
                    mult1 = 0.;
                    if (ul1 < 0.)
                        mult1 = 1.;
                    double alf0E = state.inletVoidFraction;

                    double jmax = 0.05;
                    double jmin = 0.005;
                    if ((fabs(ug1) + fabs(ul1)) / A1 > jmax) {
                        c0 = c0E;
                        ud = udE;
                    } else if ((fabs(ug1) + fabs(ul1)) / A1 < jmin) {
                        c0 = c0D;
                        ud = udD;
                    } else {
                        double raz = (jmax - (fabs(ug1) + fabs(ul1)) / A1) / (jmax - jmin);
                        c0 = ((1. - raz) * c0E + raz * c0D);
                        ud = ((1. - raz) * udE + raz * udD);
                    }

                    if (state.cells[ind].transic > 0) {
                        c0 = (c0 * state.cells[ind].transic + state.cells[ind].c0 * (atenua - state.cells[ind].transic)) / atenua;
                        ud = (ud * state.cells[ind].transic + state.cells[ind].ud * (atenua - state.cells[ind].transic)) / atenua;
                    }
                }
            }
            if (xarr1 == 1) {

                arranjo testamapa2(dia1, ul1 / A1, ug1 / A1, rlm, rgm, viscl1 / pow(10., 3.), viscg1 / pow(10., 3.), hol0,
                                   state.cells[ind].duto.teta, tensup1, state.input.mapaArranjo, state.globals);
                xarr1 = testamapa2.verificaArr();

                driftflux::correlations::C0UdDisperso(rlm, rgm, tensup1, alf0, nrey, nreyl, ug1, ul1, dia1, state.cells[ind].duto.rug, ang,
                             c0, ud, correcHor, state.cells[ind].estabCol, state.selectors.dispersed);
                if (xarr1 == -2) {
                    driftflux::correlations::C0UdAnularChurn(rlm, rgm, tensup1, alf0, nrey, nreyl, ug1, ul1, dia1, state.cells[ind].duto.rug, ang,
                                    c0, ud, correcHor, state.cells[ind].estabCol, state.selectors.annularChurn);
                }

                if (fabs(ug1 / state.cells[ind].duto.area) > 5. && alf0 >= 0.75) {
                    atenua = 20;
                    if (state.selectors.annularChurn == 3 && state.selectors.dispersed == 1)
                        atenua = 200;
                }
                if (state.cells[ind].arranjo != 0) {
                    if ((xarr1 != state.cells[ind].arranjo || state.cells[ind].transic > 0)) {
                        if (xarr1 != state.cells[ind].arranjo && state.cells[ind].transic > 0)
                            state.cells[ind].transic = 0;
                        state.cells[ind].transic++;
                        if (state.cells[ind].transic > atenua - 1)
                            state.cells[ind].transic = 0;
                    }
                } else
                    state.cells[ind].transic = 0;
                if (state.cells[ind].transic > 0) {
                    c0 = (c0 * state.cells[ind].transic + state.cells[ind].c0 * (atenua - state.cells[ind].transic)) / atenua;
                    ud = (ud * state.cells[ind].transic + state.cells[ind].ud * (atenua - state.cells[ind].transic)) / atenua;
                }
                state.cells[ind].arranjo = xarr1;
            }
            state.cells[ind].c0Spare = c0;
            state.cells[ind].udSpare = ud;
        }
    }
    if (state.input.escorregaTran == 0) {
        double ulsmed = state.cells[ind].QL / state.cells[ind].duto.area;
        double correcaoUd = 1 - (ulsmed - 0.15) / 0.35;
        double correcaoCo = c0 - (c0 - 1) * (ulsmed - 0.15) / 0.35;
        if (correcaoUd > 1.)
            correcaoUd = 1.;
        if (correcaoUd < 0.)
            correcaoUd = 0.;
        if (correcaoCo > c0)
            correcaoCo = c0;
        if (correcaoCo < 1)
            correcaoCo = 1;
        c0 = 1. + 0 * correcaoCo;
        ud = 0. + 0 * ud * correcaoUd;
    }
}

void bufferedInitialization(const ClosureState &state, int ind, double &c0, double &ud) {
    int timeStep = 20;
    if (state.cells[ind].dt < 0.8)
        timeStep *= (0.8 / state.cells[ind].dt);
    c0 = 1.;
    ud = 0.;
    if (state.cells[ind].velPig < 0 && state.cells[ind].estadoPig == 1) {
        c0 = 1.;
        ud = 0.;
        state.cells[ind].arranjo = 1.;

    } else if ((state.cells[ind].acsr.tipo != 4 || state.cells[ind].acsr.bcs.freqnova <= 1.)) {
        double hns;

        hns = 1 - state.inletVoidFraction;

        if ((state.cells[ind].MCBuf - state.cells[ind].MliqiniBuf) < 0)
            hns = 1. - state.cells[ind].alfPigE;
        if (fabs(state.cells[ind].MCBuf - state.cells[ind].MliqiniBuf) < (*state.globals).localtiny * 1e-5) {
            if (fabs(state.cells[ind].alfPigE) < (*state.globals).localtiny && fabs(1. - state.cells[ind].alfPigE) < (*state.globals).localtiny && fabs(state.inletVoidFraction) > (*state.globals).localtiny && fabs(1. - state.inletVoidFraction) > (*state.globals).localtiny)
                hns = 1 - state.inletVoidFraction;
            else if (fabs(state.inletVoidFraction) < (*state.globals).localtiny && fabs(1. - state.cells[ind - 1].alfPigD) < (*state.globals).localtiny && fabs(state.cells[ind].alfPigE) > (*state.globals).localtiny && fabs(1. - state.cells[ind].alfPigE) > (*state.globals).localtiny)
                hns = 1. - state.cells[ind].alfPigE;
            else
                hns = 0.5 * (1. - state.cells[ind].alfPigE + 1. - state.inletVoidFraction);
        }
        if (hns < (*state.globals).localtiny || hns > 1. - (*state.globals).localtiny)
            hns = 0.5 * (1. - state.cells[ind].alfPigE + 1. - state.inletVoidFraction);

        double hol0 = hns;
        double alf0 = 1 - hol0;
        double alf1;
        alf1 = state.cells[ind].alf;

        double alfneg;
        if (ind > 1)
            alfneg = state.inletVoidFraction;
        else if (ind > 0)
            alfneg = state.inletVoidFraction;
        else
            alfneg = state.cells[ind].alf;

        double betI = state.cells[ind].betL;
        if (ind > 0)
            betI = state.inletColumnFraction;
        if (state.cells[ind].QL < 0.)
            betI = state.cells[ind].betPigE; // testeBeta
        betI = state.cells[ind].betPigE;     // duvidabeta
        double betneg;
        if (ind > 0) {
            betneg = state.inletColumnFraction;

        } else
            betneg = state.cells[ind].bet;

        double pmed;
        double pmed0 = 0.;

        pmed = state.inletPressure;
        if (ind > 0)
            pmed0 = state.inletPressure;
        else
            pmed0 = state.inletPressure;
        double tmed = state.inletTemperature;

        double tmed0;
        tmed0 = tmed;

        double correcHor = 1.;
        if (fabs(state.cells[ind].duto.teta) < 1e-10) {
            if (state.cells[ind].acsr.tipo != 5 || state.cells[ind].acsr.chk.AreaGarg > 1e-10) {
                if (state.cells[ind].angEsq < 0 && state.cells[ind].angDir < 0)
                    correcHor = -1.;
                else if (state.cells[ind].angEsq > 0 && state.cells[ind].angDir > 0)
                    correcHor = 1.;
            } else {
                if (state.cells[ind].angDir < 0)
                    correcHor = -1.;
                else if (state.cells[ind].angDir > 0)
                    correcHor = 1.;
            }
        }

        double rlm;
        double viscl1;
        double tensup1;
        if (state.cells[ind].MliqiniBuf < 0.) { // testeBeta
            rlm = (1 - betI) * state.cells[ind].flui.MasEspLiq(pmed, tmed) + betI * state.cells[ind].fluicol.MasEspFlu(pmed, tmed);
            viscl1 = (1 - betI) * state.cells[ind].flui.ViscOleo(pmed, tmed) + betI * state.cells[ind].fluicol.VisFlu(pmed, tmed);
            tensup1 = (1 - betI) * state.cells[ind].flui.TensSuper(pmed, tmed) + betI * state.cells[ind].fluicol.TensSuper(pmed, tmed);
        } else {
            rlm = (1 - betI) * (*state.cells[ind].fluiL).MasEspLiq(pmed, tmed) + betI * state.cells[ind].fluicol.MasEspFlu(pmed, tmed);
            viscl1 = (1 - betI) * (*state.cells[ind].fluiL).ViscOleo(pmed, tmed) + betI * state.cells[ind].fluicol.VisFlu(pmed, tmed);
            tensup1 = (1 - betI) * (*state.cells[ind].fluiL).TensSuper(pmed, tmed) + betI * state.cells[ind].fluicol.TensSuper(pmed, tmed);
        }

        double rgm;
        double viscg1;
        if ((state.cells[ind].MCBuf - state.cells[ind].MliqiniBuf) < 0.) {
            rgm = state.cells[ind].flui.MasEspGas(pmed, tmed);
            viscg1 = state.cells[ind].flui.ViscGas(pmed, tmed);
        } else {
            rgm = (*state.cells[ind].fluiL).MasEspGas(pmed, tmed);
            viscg1 = (*state.cells[ind].fluiL).ViscGas(pmed, tmed);
        }

        double ug1 = (state.cells[ind].MCBuf - state.cells[ind].MliqiniBuf) / rgm;
        double ul1 = state.cells[ind].MliqiniBuf / rlm;
        double dia1 = state.cells[ind].duto.a;
        if (ind > 0 && ug1 >= 0)
            dia1 = state.cells[ind - 1].duto.a;
        double A1 = M_PI * dia1 * dia1 / 4.;

        double rmed = hns * rlm + (1 - hns) * rgm;
        double visc = (hns * viscl1 + (1 - hns) * viscg1) / pow(10., 3.);
        double nrey = dia1 * rmed * (fabs(ug1) / A1 + fabs(ul1) / A1) / visc;
        double nreyl = dia1 * rlm * (fabs(ug1) / A1 + fabs(ul1) / A1) / (viscl1 / 1000.);

        int xarr1 = 1;
        double dtot = state.cells[ind].dxL + state.cells[ind].dx;
        double razL = state.cells[ind].dxL;
        double raz = state.cells[ind].dx;
        double ang = (raz * state.cells[ind].dutoL.teta + razL * state.cells[ind].duto.teta) / dtot;
        double atenua = 20.;
        if (nrey > 1e-30) {
            if (fabs(0 * ang + 1 * state.cells[ind].duto.teta) < 45. * M_PI / 180. && hol0 < 0.99 && hol0 > 0.01 && ind < state.lastCell - 1) {
                double ug0;
                double ul0;

                ug1 = (state.cells[ind].MCBuf - state.cells[ind].MliqiniBuf) / rgm;
                ul1 = state.cells[ind].MliqiniBuf / rlm;

                ug0 = (state.cells[ind].MCBuf - state.cells[ind].MliqiniBuf) / (*state.cells[ind].fluiL).MasEspGas(pmed0, tmed0);
                ul0 = state.cells[ind].MliqiniBuf / ((1 - betneg) * (*state.cells[ind].fluiL).MasEspLiq(pmed0, tmed0) + betneg * state.cells[ind].fluicol.MasEspFlu(pmed0, tmed0));

                xarr1 = state.cells[ind].arranjo;
                if (xarr1 == -1) {

                    double c0D;
                    double udD;
                    double c0E;
                    double udE;
                    driftflux::correlations::C0UdDisperso(rlm, rgm, tensup1, alf0, nrey, nreyl, ug1, ul1, dia1, state.cells[ind].duto.rug, ang,
                                 c0D, udD, correcHor, state.cells[ind].estabCol, state.selectors.dispersed);
                    driftflux::correlations::C0UdEstratificado(rlm, rgm, tensup1, alf0, nrey, nreyl, ug1, ul1, dia1, state.cells[ind].duto.rug, ang,
                                      c0E, udE, correcHor, state.cells[ind].estabCol, state.selectors.stratified);
                    double mult0, mult1;
                    mult0 = 1.;
                    if (ul0 < 0.)
                        mult0 = 0.;
                    mult1 = 0.;
                    if (ul1 < 0.)
                        mult1 = 1.;

                    double jmax = 0.05;
                    double jmin = 0.005;
                    if ((fabs(ug1) + fabs(ul1)) / A1 > jmax) {
                        c0 = c0E;
                        ud = udE;
                    } else if ((fabs(ug1) + fabs(ul1)) / A1 < jmin) {
                        c0 = c0D;
                        ud = udD;
                    } else {
                        double raz = (jmax - (fabs(ug1) + fabs(ul1)) / A1) / (jmax - jmin);
                        c0 = ((1. - raz) * c0E + raz * c0D);
                        ud = ((1. - raz) * udE + raz * udD);
                    }

                    if (state.cells[ind].transic > 0) {
                        c0 = (c0 * state.cells[ind].transic + state.cells[ind].c0 * (atenua - state.cells[ind].transic)) / atenua;
                        ud = (ud * state.cells[ind].transic + state.cells[ind].ud * (atenua - state.cells[ind].transic)) / atenua;
                    }
                }
            }
            if (xarr1 != -1) {

                driftflux::correlations::C0UdDisperso(rlm, rgm, tensup1, alf0, nrey, nreyl, ug1, ul1, dia1, state.cells[ind].duto.rug, ang,
                             c0, ud, correcHor, state.cells[ind].estabCol, state.selectors.dispersed);
                if (xarr1 == -2) {
                    driftflux::correlations::C0UdAnularChurn(rlm, rgm, tensup1, alf0, nrey, nreyl, ug1, ul1, dia1, state.cells[ind].duto.rug, ang,
                                    c0, ud, correcHor, state.cells[ind].estabCol, state.selectors.annularChurn);
                }
                if (fabs(ug1 / state.cells[ind].duto.area) > 5. && alf0 >= 0.75) {
                    atenua = 20;
                    if (state.selectors.annularChurn == 3 && state.selectors.dispersed == 1)
                        atenua = 200;
                }

                if (state.cells[ind].transic > 0) {
                    c0 = (c0 * state.cells[ind].transic + state.cells[ind].c0 * (atenua - state.cells[ind].transic)) / atenua;
                    ud = (ud * state.cells[ind].transic + state.cells[ind].ud * (atenua - state.cells[ind].transic)) / atenua;
                }
                state.cells[ind].arranjo = xarr1;
            }
            state.cells[ind].c0Spare = c0;
            state.cells[ind].udSpare = ud;
        }
    }
    if (state.input.escorregaTran == 0) {
        double ulsmed = state.cells[ind].QL / state.cells[ind].duto.area;
        double correcaoUd = 1 - (ulsmed - 0.15) / 0.35;
        double correcaoCo = c0 - (c0 - 1) * (ulsmed - 0.15) / 0.35;
        if (correcaoUd > 1.)
            correcaoUd = 1.;
        if (correcaoUd < 0.)
            correcaoUd = 0.;
        if (correcaoCo > c0)
            correcaoCo = c0;
        if (correcaoCo < 1)
            correcaoCo = 1;
        c0 = 1. + 0 * correcaoCo;
        ud = 0. + 0 * ud * correcaoUd;
    }
}

void steadyState(const ClosureState &state, int ind, double &c0, double &ud) {

    c0 = 1.;
    ud = 0.;
    if (state.cells[ind].acsr.tipo != 4 || state.cells[ind].acsr.bcs.freqnova <= 1.) {
        double hns;
        double razdx = state.cells[ind].dx / (state.cells[ind].dx + state.cells[ind].dxL);
        double razdx0;
        if (ind > 0)
            razdx0 = state.cells[ind - 1].dx / (state.cells[ind - 1].dx + state.cells[ind - 1].dxL);
        else
            razdx0 = razdx;
        hns = 1. - state.cells[ind].alf;

        double hol0 = hns;
        double alf0 = 1 - hol0;
        double alf1;
        alf1 = state.cells[ind].alf;

        double alfneg;
        if (ind > 1)
            alfneg = state.cells[ind - 2].alf;
        else if (ind > 0)
            alfneg = state.cells[ind - 1].alf;
        else
            alfneg = state.cells[ind].alf;

        double betI = state.cells[ind].betL;
        double betneg = state.cells[ind].betL;

        double pmed;
        double pmed0 = 0.;

        pmed = state.cells[ind].presaux;
        if (ind > 0)
            pmed0 = state.cells[ind - 1].presaux;
        else
            pmed0 = state.cells[ind].presaux;
        double tmed;
        if (state.steadyIteration != 0 && state.input.AceleraConvergPerm == 0)
            tmed = razdx * state.cells[ind].temp + (1 - razdx) * state.cells[ind].tempL;
        else
            tmed = state.cells[ind - 1].temp;
        double tmed0;
        if (ind > 0 && state.input.AceleraConvergPerm == 0)
            tmed0 = razdx0 * state.cells[ind - 1].temp + (1 - razdx0) * state.cells[ind - 1].tempL;
        else
            tmed0 = tmed;

        double correcHor = 1.;
        if (fabs(state.cells[ind].duto.teta) < 1e-10) {
            if (state.cells[ind].angEsq < 0 && state.cells[ind].angDir < 0)
                correcHor = -1.;
            else if (state.cells[ind].angEsq > 0 && state.cells[ind].angDir > 0)
                correcHor = 1.;
        }

        double rlm;
        double viscl1;
        double tensup1;
        if (ind > 0)
            rlm = (1 - betI) * state.cells[ind - 1].flui.MasEspLiq(pmed0, tmed0) + betI * state.cells[ind - 1].fluicol.MasEspFlu(pmed0, tmed0);
        else
            rlm = (1 - betI) * (*state.cells[ind].fluiL).MasEspLiq(pmed, tmed) + betI * state.cells[ind].fluicol.MasEspFlu(pmed, tmed);

        viscl1 = (1 - betI) * (*state.cells[ind].fluiL).ViscOleo(pmed0, tmed0) + betI * state.cells[ind].fluicol.VisFlu(pmed0, tmed0);
        tensup1 = (1 - betI) * (*state.cells[ind].fluiL).TensSuper(pmed, tmed) + betI * state.cells[ind].fluicol.TensSuper(pmed, tmed);

        double rgm;
        double viscg1;
        if (ind > 0)
            rgm = state.cells[ind - 1].flui.MasEspGas(pmed0, tmed0);
        else
            rgm = (*state.cells[ind].fluiL).MasEspGas(pmed, tmed);
        viscg1 = (*state.cells[ind].fluiL).ViscGas(pmed, tmed);

        double ug1 = fabs(state.cells[ind].MC - state.cells[ind].Mliqini) / rgm;
        double ul1 = fabs(state.cells[ind].Mliqini) / rlm;
        double dia1 = state.cells[ind].duto.a;
        if (ind > 0 && ug1 >= 0)
            dia1 = state.cells[ind - 1].duto.a;
        double A1 = M_PI * dia1 * dia1 / 4.;

        double rmed = hns * rlm + (1 - hns) * rgm;
        double visc = (hns * viscl1 + (1 - hns) * viscg1) / pow(10., 3.);
        double nrey = dia1 * rmed * (fabs(ug1) / A1 + fabs(ul1) / A1) / visc;
        double nreyl = dia1 * rlm * (fabs(ug1) / A1 + fabs(ul1) / A1) / (viscl1 / 1000.);

        int xarr1 = 1;
        double dtot = state.cells[ind].dxL + state.cells[ind].dx;
        double razL = state.cells[ind].dxL;
        double raz = state.cells[ind].dx;
        double ang = (razL * state.cells[ind].dutoL.teta + raz * state.cells[ind].duto.teta) / dtot;
        double sinalAng = 1.;
        if (fabs(state.cells[ind].MC) > 1e-15)
            sinalAng = state.cells[ind].MC / fabs(state.cells[ind].MC);
        if (rgm < 0.9 * rlm) {
            if (nrey > 1e-30) {
                if (fabs(0 * ang + sinalAng * state.cells[ind].duto.teta) < 45. * M_PI / 180. && hol0 < 0.99 && hol0 > 0.01 && ind < state.lastCell - 1) {
                    double ug0;
                    double ul0;

                    ug1 = fabs(state.cells[ind].MC - state.cells[ind].Mliqini) / rgm;
                    ul1 = fabs(state.cells[ind].Mliqini) / rlm;

                    ug0 = ug1;
                    ul0 = ul1;
                    if (ind > 0) {
                        ug0 = fabs(state.cells[ind - 1].MC - state.cells[ind - 1].Mliqini) / state.cells[ind].flui.MasEspGas(pmed0, tmed0);
                        ul0 = fabs(state.cells[ind - 1].Mliqini) / ((1 - betneg) * state.cells[ind].flui.MasEspLiq(pmed0, tmed0) + betneg * state.cells[ind].fluicol.MasEspFlu(pmed0, tmed0));
                    }

                    estratificado testamapa(dia1, ul1, ug1, rlm, rgm, viscl1 / pow(10., 3.), viscg1 / pow(10., 3.), hol0,
                                            sinalAng * state.cells[ind].duto.teta, state.cells[ind].duto.rug / dia1);

                    testamapa.mapaTD();
                    xarr1 = testamapa.arr;
                    if (xarr1 == -1) {
                        if (((state.cells[ind].arranjo != xarr1) || state.cells[ind].transic > 0)) {
                            if ((state.cells[ind].arranjo != xarr1) && state.cells[ind].transic > 0)
                                state.cells[ind].transic = 0;
                            state.cells[ind].transic++;
                            if (state.cells[ind].transic > 19)
                                state.cells[ind].transic = 0;
                        }
                        state.cells[ind].arranjo = xarr1 = testamapa.arr;
                        if (ind > 0) {
                            state.cells[ind - 1].arranjoR = testamapa.arr;
                            state.cells[ind - 1].perdaEstratL = testamapa.fatorperdaLiq;
                            state.cells[ind - 1].perdaEstratG = testamapa.fatorperdaGas;
                        }

                        double c0D;
                        double udD;
                        double c0E;
                        double udE;
                        driftflux::correlations::C0UdDisperso(rlm, rgm, tensup1, alf0, nrey, nreyl, ug1, ul1, dia1, state.cells[ind].duto.rug, ang,
                                     c0D, udD, correcHor, state.cells[ind].estabCol, state.selectors.dispersed);
                        driftflux::correlations::C0UdEstratificado(rlm, rgm, tensup1, alf0, nrey, nreyl, ug1, ul1, dia1, state.cells[ind].duto.rug, ang,
                                          c0E, udE, correcHor, state.cells[ind].estabCol, state.selectors.stratified);

                        double mult0, mult1;
                        mult0 = 1.;
                        if (ul0 < 0.)
                            mult0 = 0.;
                        mult1 = 0.;
                        if (ul1 < 0.)
                            mult1 = 1.;
                        double alf0E = state.cells[ind].alf;
                        if (ind > 0)
                            alf0E = state.cells[ind - 1].alf;

                        double jmax = 0.05;
                        double jmin = 0.005;
                        if ((fabs(ug1) + fabs(ul1)) / A1 > jmax) {
                            c0 = c0E;
                            ud = udE;
                        } else if ((fabs(ug1) + fabs(ul1)) / A1 < jmin) {
                            c0 = c0D;
                            ud = udD;
                        } else {
                            double raz = (jmax - (fabs(ug1) + fabs(ul1)) / A1) / (jmax - jmin);
                            c0 = ((1. - raz) * c0E + raz * c0D);
                            ud = ((1. - raz) * udE + raz * udD);
                        }
                    }
                }
                if (xarr1 == 1) {

                    arranjo testamapa2(dia1, ul1 / A1, ug1 / A1, rlm, rgm, viscl1 / pow(10., 3.), viscg1 / pow(10., 3.), hol0,
                                       sinalAng * state.cells[ind].duto.teta, tensup1, state.input.mapaArranjo, state.globals);
                    xarr1 = testamapa2.verificaArr();

                    driftflux::correlations::C0UdDisperso(rlm, rgm, tensup1, alf0, nrey, nreyl, ug1, ul1, dia1, state.cells[ind].duto.rug, ang,
                                 c0, ud, correcHor, state.cells[ind].estabCol, state.selectors.dispersed);
                    if (xarr1 == -2) {
                        driftflux::correlations::C0UdAnularChurn(rlm, rgm, tensup1, alf0, nrey, nreyl, ug1, ul1, dia1, state.cells[ind].duto.rug, ang,
                                        c0, ud, correcHor, state.cells[ind].estabCol, state.selectors.annularChurn);
                    }
                    state.cells[ind].arranjo = xarr1;
                    if (ind > 0)
                        state.cells[ind - 1].arranjoR = xarr1;
                }
                state.cells[ind].c0Spare = c0;
                state.cells[ind].udSpare = ud;
            }
        } else {
            c0 = 1.;
            ud = 0.;
            state.cells[ind].arranjo = 1;
            if (ind > 0)
                state.cells[ind - 1].arranjoR = 1;
            state.cells[ind].c0Spare = c0;
            state.cells[ind].udSpare = ud;
        }
    }
    if (state.input.escorregaPerm == 0) {
        double ulsmed = state.cells[ind].QL / state.cells[ind].duto.area;
        double correcaoUd = 1 - (ulsmed - 0.15) / 0.35;
        double correcaoCo = c0 - (c0 - 1) * (ulsmed - 0.15) / 0.35;
        if (correcaoUd > 1.)
            correcaoUd = 1.;
        if (correcaoUd < 0.)
            correcaoUd = 0.;
        if (correcaoCo > c0)
            correcaoCo = c0;
        if (correcaoCo < 1)
            correcaoCo = 1;
        c0 = 1. + 0 * correcaoCo;
        ud = 0. + 0. * ud * correcaoUd;
    }
}

}  // namespace coefficient
}  // namespace driftflux
