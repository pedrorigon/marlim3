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
 * betneg -> upstreamLiquidFlowRate -> mult0, dead in all five, is A3-07.
 *
 * The `// duvidabeta` and `// testeBeta` markers below are the original
 * author's, kept where they were. They are Portuguese for "beta doubt" and
 * "beta test", and they sit on the assignments of betI and betneg -- which is
 * anomaly A3-01, where a later unconditional assignment makes the selection
 * above it dead. They are evidence that someone was unsure here, so they stay.
 */

namespace {

/*
 * Data Source Policy.
 *
 * The contract named five sources, one per variant. Measurement says there are
 * THREE: CalcC0Ud and CalcC0UdIni read exactly the same fields, and so do
 * CalcC0UdBuf and CalcC0UdIniBuf -- the initialisation variants are not a
 * different source, they are different control flow over the same source, which
 * is not something a data-source policy can express. What actually varies is
 * instantaneous versus buffered versus steady state.
 *
 * The hooks are called AT THE POINT OF USE, never hoisted into a local at the
 * top of a body. A call substituted for an expression is evaluated where the
 * expression was; a value read once and reused is not, and the difference is
 * exactly how a conditional read became unconditional in the root-finding
 * stage. gasForSign and gasFlowRate are separate hooks because they are
 * separate expressions in the instantaneous variants: the sign tests read QG,
 * while the flow rate is MC - Mliqini.
 *
 * Stateless structs with static members, resolved at compile time, defined in
 * the same translation unit as their only callers: no indirect call survives.
 */

/// Instantaneous state -- SProd::CalcC0Ud and SProd::CalcC0UdIni.
struct InstantaneousSource {
    /// Value whose SIGN the variant tests. Not the flow rate: these variants
    /// branch on QG but divide MC - Mliqini by the density.
    static double gasForSign(const Cel *cells, int index) { return cells[index].QG; }
    /// Gas mass flow: total minus its liquid part.
    static double gasFlowRate(const Cel *cells, int index) {
        return cells[index].MC - cells[index].Mliqini;
    }
    /// Liquid mass flow.
    static double liquidFlowRate(const Cel *cells, int index) { return cells[index].Mliqini; }
};

/// Buffered network state -- SProd::CalcC0UdBuf and SProd::CalcC0UdIniBuf.
struct BufferedSource {
    /// Value whose SIGN the variant tests. Here it is the same expression as the
    /// flow rate, unlike the instantaneous source; both hooks are kept so the
    /// two roles stay distinguishable at the call sites.
    static double gasForSign(const Cel *cells, int index) {
        return cells[index].MCBuf - cells[index].MliqiniBuf;
    }
    /// Buffered gas mass flow: total minus its liquid part.
    static double gasFlowRate(const Cel *cells, int index) {
        return cells[index].MCBuf - cells[index].MliqiniBuf;
    }
    /// Buffered liquid mass flow.
    static double liquidFlowRate(const Cel *cells, int index) { return cells[index].MliqiniBuf; }
};

/// Steady state -- SProd::CalcC0UdPerm. It takes the magnitude of both rates and
/// never tests their sign, so it has no gasForSign.
struct SteadyStateSource {
    /// Gas mass flow, magnitude only -- the steady-state variant never branches
    /// on direction, which is why this source has no gasForSign.
    static double gasFlowRate(const Cel *cells, int index) {
        return fabs(cells[index].MC - cells[index].Mliqini);
    }
    /// Liquid mass flow, magnitude only.
    static double liquidFlowRate(const Cel *cells, int index) {
        return fabs(cells[index].Mliqini);
    }
};

/// No-slip liquid holdup at the face, for the transient variants.
///
/// Identical in CalcC0Ud and CalcC0UdBuf once the sign source is the policy's:
/// the two differed only in reading QG against MCBuf - MliqiniBuf, which is
/// exactly what gasForSign is for. The hook is called at each test, never read
/// once into a local, so a value the original consulted twice is still consulted
/// twice.
///
/// The two localtiny guards are the original's: the first replaces a holdup
/// derived from a vanishing gas rate, the second rejects one that has collapsed
/// onto either end of its range.
template <typename Source>
double transientNoSlipHoldup(const ClosureState &state, int cellIndex) {
    double noSlipLiquidHoldup;
    if (cellIndex > 0)
        noSlipLiquidHoldup = 1. - state.cells[cellIndex - 1].alfPigD;
    else
        noSlipLiquidHoldup = 1. - state.cells[cellIndex].alf;
    if (Source::gasForSign(state.cells, cellIndex) < 0)
        noSlipLiquidHoldup = 1. - state.cells[cellIndex].alfPigE;
    if (fabs(Source::gasForSign(state.cells, cellIndex)) < (*state.globals).localtiny * 1e-5) {
        if (fabs(state.cells[cellIndex].alfPigE) < (*state.globals).localtiny && fabs(1. - state.cells[cellIndex].alfPigE) < (*state.globals).localtiny && fabs(state.cells[cellIndex - 1].alfPigD) > (*state.globals).localtiny && fabs(1. - state.cells[cellIndex - 1].alfPigD) > (*state.globals).localtiny)
            noSlipLiquidHoldup = 1. - state.cells[cellIndex - 1].alfPigD;
        else if (fabs(state.cells[cellIndex - 1].alfPigD) < (*state.globals).localtiny && fabs(1. - state.cells[cellIndex - 1].alfPigD) < (*state.globals).localtiny && fabs(state.cells[cellIndex].alfPigE) > (*state.globals).localtiny && fabs(1. - state.cells[cellIndex].alfPigE) > (*state.globals).localtiny)
            noSlipLiquidHoldup = 1. - state.cells[cellIndex].alfPigE;
        else
            noSlipLiquidHoldup = 0.5 * (1. - state.cells[cellIndex].alfPigE + 1. - state.cells[cellIndex - 1].alfPigD);
    }
    if (noSlipLiquidHoldup < (*state.globals).localtiny || noSlipLiquidHoldup > 1. - (*state.globals).localtiny)
        noSlipLiquidHoldup = 0.5 * (1. - state.cells[cellIndex].alfPigE + 1. - state.cells[cellIndex - 1].alfPigD);
    return noSlipLiquidHoldup;
}

/// Duct inclination seen by the face, for the transient variants.
///
/// The length-weighted mean of the two duct inclinations, overridden two cells
/// downstream of a shut choke by whichever single duct the flow actually comes
/// from. Identical in CalcC0Ud and CalcC0UdBuf once the sign source is the
/// policy's.
///
/// The weighting reads cellLength against dutoL and leftCellLength against
/// duto -- crossed, which is the usual interpolation. CalcC0UdPerm weights it
/// the other way round and therefore does NOT use this function; that is
/// anomaly A3-02, preserved.
template <typename Source>
double transientInclinationAngle(const ClosureState &state, int cellIndex) {
    const double totalLength = state.cells[cellIndex].dxL + state.cells[cellIndex].dx;
    const double cellLength = state.cells[cellIndex].dx;
    const double leftCellLength = state.cells[cellIndex].dxL;
    double inclinationAngle = (cellLength * state.cells[cellIndex].dutoL.teta + leftCellLength * state.cells[cellIndex].duto.teta) / totalLength;
    if (cellIndex >= 2) {
        if (state.cells[cellIndex - 2].acsr.tipo == 5 && state.cells[cellIndex - 2].acsr.chk.AreaGarg <= (1e-3)) {
            if (Source::gasForSign(state.cells, cellIndex) >= 0)
                inclinationAngle = state.cells[cellIndex].duto.teta;
            else
                inclinationAngle = state.cells[cellIndex].dutoR.teta;
        } else {
            if (Source::gasForSign(state.cells, cellIndex) >= 0)
                inclinationAngle = state.cells[cellIndex].dutoL.teta;
            else
                inclinationAngle = state.cells[cellIndex].duto.teta;
        }
    }
    return inclinationAngle;
}

/// Sign correction applied to the drift term on a horizontal face.
///
/// Identical in the four transient variants except for WHICH cell's accessory
/// opens the guard: CalcC0Ud reads the upstream one, the other three read the
/// face's own. That difference is anomaly A3-01's neighbour -- site 14 of the
/// normalized comparison -- and it is preserved by making the accessory cell a
/// parameter rather than by picking one and calling the rest wrong.
///
/// The two arms differ in more than the accessory: the first requires both
/// junction angles to agree in sign, the second reads the downstream angle
/// alone. So the guard is not decorative, and a call site that passed the wrong
/// cell would change results wherever the two cells carry different accessories.
double horizontalCorrectionOf(const ClosureState &state, int cellIndex, int accessoryCellIndex) {
    double horizontalCorrection = 1.;
    if (fabs(state.cells[cellIndex].duto.teta) < 1e-10) {
        if (state.cells[accessoryCellIndex].acsr.tipo != 5 || state.cells[accessoryCellIndex].acsr.chk.AreaGarg > 1e-10) {
            if (state.cells[cellIndex].angEsq < 0 && state.cells[cellIndex].angDir < 0)
                horizontalCorrection = -1.;
            else if (state.cells[cellIndex].angEsq > 0 && state.cells[cellIndex].angDir > 0)
                horizontalCorrection = 1.;
        } else {
            if (state.cells[cellIndex].angDir < 0)
                horizontalCorrection = -1.;
            else if (state.cells[cellIndex].angDir > 0)
                horizontalCorrection = 1.;
        }
    }
    return horizontalCorrection;
}

/// No-slip liquid holdup at the face, for the initialisation variants.
///
/// Same shape as the transient one, reading the inlet void fraction instead of
/// the neighbouring cell's. Identical in CalcC0UdIni and CalcC0UdIniBuf once the
/// sign source is the policy's.
///
/// One guard mixes the two: the second condition of the else-if still tests
/// cells[cellIndex - 1].alfPigD while everything around it moved to the inlet
/// fraction. That is anomaly A3-04, an incomplete substitution in the original,
/// and it is preserved exactly -- correcting it here would change results in a
/// variant no model in the corpus executes.
template <typename Source>
double inletNoSlipHoldup(const ClosureState &state, int cellIndex) {
    double noSlipLiquidHoldup = 1 - state.inletVoidFraction;
    if (Source::gasForSign(state.cells, cellIndex) < 0)
        noSlipLiquidHoldup = 1. - state.cells[cellIndex].alfPigE;
    if (fabs(Source::gasForSign(state.cells, cellIndex)) < (*state.globals).localtiny * 1e-5) {
        if (fabs(state.cells[cellIndex].alfPigE) < (*state.globals).localtiny && fabs(1. - state.cells[cellIndex].alfPigE) < (*state.globals).localtiny && fabs(state.inletVoidFraction) > (*state.globals).localtiny && fabs(1. - state.inletVoidFraction) > (*state.globals).localtiny)
            noSlipLiquidHoldup = 1 - state.inletVoidFraction;
        else if (fabs(state.inletVoidFraction) < (*state.globals).localtiny && fabs(1. - state.cells[cellIndex - 1].alfPigD) < (*state.globals).localtiny && fabs(state.cells[cellIndex].alfPigE) > (*state.globals).localtiny && fabs(1. - state.cells[cellIndex].alfPigE) > (*state.globals).localtiny)
            noSlipLiquidHoldup = 1. - state.cells[cellIndex].alfPigE;
        else
            noSlipLiquidHoldup = 0.5 * (1. - state.cells[cellIndex].alfPigE + 1. - state.inletVoidFraction);
    }
    if (noSlipLiquidHoldup < (*state.globals).localtiny || noSlipLiquidHoldup > 1. - (*state.globals).localtiny)
        noSlipLiquidHoldup = 0.5 * (1. - state.cells[cellIndex].alfPigE + 1. - state.inletVoidFraction);
    return noSlipLiquidHoldup;
}

/// No slip while a pig occupies the face: the drift terms are forced off and the
/// flow pattern is pinned to slug.
///
/// This is the body of BOTH arms of the pig guard in the two transient variants,
/// textually identical in all four places -- the conditions differ, what they do
/// does not. Sharing it does not merge the branches: it makes the fact that they
/// agree visible, the way the root-finding stage argued a preserved branch is
/// worth keeping precisely as a marker.
void applyPigOverride(const ClosureState &state, int cellIndex, double &c0, double &ud) {
    c0 = 1.;
    ud = 0.;
    state.cells[cellIndex].arranjo = 1.;
    state.cells[cellIndex - 1].arranjoR = 1.;
    state.cells[cellIndex - 1].perdaEstratL = 0.;
    state.cells[cellIndex - 1].perdaEstratG = 0.;
}

/// The phase properties the flow scales are built from.
///
/// Same reasoning as MixtureProperties: five adjacent doubles is five chances to
/// transpose a pair with nothing to catch it. This one bit before it was caught:
/// the densities were passed gas-then-liquid while the viscosities went
/// liquid-then-gas, an asymmetry with no reason behind it and no way for the
/// compiler to notice a call site that got it wrong.
struct PhaseProperties {
    double liquidDensity;       ///< rlm
    double gasDensity;          ///< rgm
    double liquidViscosity;     ///< viscl1
    double gasViscosity;        ///< viscg1
    double noSlipLiquidHoldup;  ///< hns
};

/// The scalars the closure helpers below read, named instead of counted.
///
/// evaluateRegimePair and evaluateDispersedOrAnnular took eleven and twelve
/// doubles positionally, in the order the correlation signatures use. That order
/// is a real convention and worth keeping, but eleven adjacent doubles is also
/// eleven chances to transpose a pair silently -- every one of them is the same
/// type, so neither the compiler nor the sweep would say a word about a call
/// site that swapped two. Built once per body with designated initializers,
/// against locals of the same name, a transposition is visible on the line
/// where it happens.
///
/// Constructed once per variant, immediately before the Reynolds guard that
/// gates both helpers.
///
/// The fields are REFERENCES, for the reason ClosureState holds references:
/// several of these locals are assigned again further down, and a copy taken
/// here would freeze the value at construction rather than at use. Nothing
/// between construction and use writes them today -- but "nothing writes it
/// today" is how a read moves without anyone noticing.
struct MixtureProperties {
    const double &liquidDensity;         ///< rlm
    const double &gasDensity;            ///< rgm
    const double &surfaceTension;        ///< tensup1
    const double &voidFraction;          ///< alf0
    const double &gasFlowRate;           ///< ug1, a volumetric rate
    const double &liquidFlowRate;        ///< ul1, a volumetric rate
    const double &diameter;              ///< dia1
    const double &flowArea;              ///< A1
    const double &mixtureReynolds;       ///< nrey
    const double &liquidReynolds;        ///< nreyl
    const double &inclinationAngle;      ///< ang
    const double &horizontalCorrection;  ///< correcHor
};

/// Pressure and temperature the property model is evaluated at, for CalcC0Ud.
///
/// Used once, and named because it is one question answered over twenty lines:
/// at what conditions are the phase properties taken.
///
/// Three assignments here are immediately overwritten, and that is the
/// original's shape, not an oversight in the move: the length-weighted mean
/// temperature is replaced by the left cell's, which is replaced again by the
/// face's or by the surface temperature. The upstream pressure and temperature
/// it also computes are never read by anything -- part of the dead chain
/// catalogued as A3-06 -- and are kept because removing them would change the
/// token stream of code no model in the corpus executes.
struct MeanConditions {
    double pressure;     ///< pmed
    double temperature;  ///< tmed
};

MeanConditions instantaneousMeanConditions(const ClosureState &state, int cellIndex,
                                           double lengthRatio, double upstreamLengthRatio) {
    double meanPressure;
    double upstreamMeanPressure = 0.;

    meanPressure = state.cells[cellIndex].presaux;
    if (cellIndex > 0)
        upstreamMeanPressure = state.cells[cellIndex - 1].presaux;
    if (cellIndex == state.lastCell)
        meanPressure = state.cells[cellIndex].pres;
    else
        upstreamMeanPressure = state.cells[cellIndex].presaux;
    double meanTemperature = lengthRatio * state.cells[cellIndex].temp + (1 - lengthRatio) * state.cells[cellIndex].tempL;
    meanTemperature = state.cells[cellIndex].tempL;
    if (state.cells[cellIndex].VTemper < 0.)
        meanTemperature = state.cells[cellIndex].temp;
    double upstreamMeanTemperature;
    if (cellIndex > 0)
        upstreamMeanTemperature = upstreamLengthRatio * state.cells[cellIndex - 1].temp + (1 - upstreamLengthRatio) * state.cells[cellIndex - 1].tempL;
    else
        upstreamMeanTemperature = meanTemperature;
    if (cellIndex < state.lastCell)
        meanTemperature = state.cells[cellIndex].temp;
    else
        meanTemperature = state.gasSurfaceTemperature;
    return MeanConditions{meanPressure, meanTemperature};
}

/// Phase properties at the face, from the instantaneous state.
///
/// Used once, by CalcC0Ud, and named rather than inlined because it is one
/// question -- what are the two phases like here -- answered over thirty lines
/// in the middle of a much longer one.
///
/// It is the only variant that consults the cached densities rpCi, rcCi and
/// rgCi: at the two ends of the line it calls the property model, and in between
/// it takes the cache. The other four always call the model. That asymmetry is
/// the original's and is preserved.
///
/// The `// testeBeta` marker on the first branch is the original author's; see
/// the note on the preserved anomalies above.
PhaseProperties instantaneousPhaseProperties(const ClosureState &state, int cellIndex, double betI,
                                             double meanPressure, double meanTemperature,
                                             double noSlipLiquidHoldup, double &surfaceTension) {
    double liquidDensity;
    double liquidViscosity;
    if (state.cells[cellIndex].QL < 0.) { // testeBeta
        if (cellIndex == 0 || cellIndex == state.lastCell)
            liquidDensity = (1 - betI) * state.cells[cellIndex].flui.MasEspLiq(meanPressure, meanTemperature) + betI * state.cells[cellIndex].fluicol.MasEspFlu(meanPressure, meanTemperature);
        else
            liquidDensity = (1 - betI) * state.cells[cellIndex].rpCi + betI * state.cells[cellIndex].rcCi;
        liquidViscosity = (1 - betI) * state.cells[cellIndex].flui.ViscOleo(meanPressure, meanTemperature) + betI * state.cells[cellIndex].fluicol.VisFlu(meanPressure, meanTemperature);
        surfaceTension = (1 - betI) * state.cells[cellIndex].flui.TensSuper(meanPressure, meanTemperature) + betI * state.cells[cellIndex].fluicol.TensSuper(meanPressure, meanTemperature);
    } else {
        if (cellIndex == 0 || cellIndex == state.lastCell)
            liquidDensity = (1 - betI) * state.cells[cellIndex - 1].flui.MasEspLiq(meanPressure, meanTemperature) + betI * state.cells[cellIndex - 1].fluicol.MasEspFlu(meanPressure, meanTemperature);
        else
            liquidDensity = (1 - betI) * state.cells[cellIndex].rpCi + betI * state.cells[cellIndex - 1].rcCi;
        liquidViscosity = (1 - betI) * state.cells[cellIndex - 1].flui.ViscOleo(meanPressure, meanTemperature) + betI * state.cells[cellIndex - 1].fluicol.VisFlu(meanPressure, meanTemperature);
        surfaceTension = (1 - betI) * state.cells[cellIndex - 1].flui.TensSuper(meanPressure, meanTemperature) + betI * state.cells[cellIndex - 1].fluicol.TensSuper(meanPressure, meanTemperature);
    }

    double gasDensity;
    double gasViscosity;
    if (state.cells[cellIndex].QG < 0.) {
        if (cellIndex == 0 || cellIndex == state.lastCell)
            gasDensity = state.cells[cellIndex].flui.MasEspGas(meanPressure, meanTemperature);
        else
            gasDensity = state.cells[cellIndex].rgCi;
        gasViscosity = state.cells[cellIndex].flui.ViscGas(meanPressure, meanTemperature);
    } else {
        if (cellIndex == 0 || cellIndex == state.lastCell)
            gasDensity = state.cells[cellIndex - 1].flui.MasEspGas(meanPressure, meanTemperature);
        else
            gasDensity = state.cells[cellIndex].rgCi;
        gasViscosity = state.cells[cellIndex - 1].flui.ViscGas(meanPressure, meanTemperature);
    }
    return PhaseProperties{liquidDensity, gasDensity, liquidViscosity, gasViscosity, noSlipLiquidHoldup};
}

/// The dispersed and stratified closures, evaluated as a pair and then blended.
///
/// They travelled as four loose doubles between the two helpers, which is four
/// values of one type in a row and no way for anything to notice a swap. Kept
/// together they are named at every use, and the two helpers now agree on one
/// shape: one fills it, the other reads it.
struct RegimePair {
    double dispersedC0;   ///< c0D
    double dispersedUd;   ///< udD
    double stratifiedC0;  ///< c0E
    double stratifiedUd;  ///< udE
};

/// The flow rates and the two Reynolds numbers built from them.
///
/// The five call sites take these apart with a structured binding, so THE ORDER
/// OF THESE FIELDS IS LOAD-BEARING: reordering them silently rebinds every call
/// site to the wrong values. It replaced six lines of hand unpacking per body,
/// which had the same hazard five times over and no reason for anyone to check
/// it -- `const double diameter = scales.area;` would have compiled.
struct FlowScales {
    double gasRate;      ///< ug1
    double liquidRate;   ///< ul1
    double diameter;     ///< dia1
    double area;         ///< A1
    double mixture;      ///< nrey
    double liquid;       ///< nreyl
};

/// Flow rates, duct size and Reynolds numbers -- the one block that is both
/// identical in all five variants AND parameterised by where the rates come
/// from. 145 tokens, proven identical across the five before being shared.
///
/// The duct diameter is chosen INSIDE this function, not passed in, because the
/// original reads cells[ind - 1].duto.a only when ind > 0 && ug1 >= 0. Taking it
/// as an argument would make that read unconditional.
///
/// The four phase properties are in one order -- liquid, then gas, for the
/// densities and again for the viscosities. They were not, at first: the
/// densities went gas-then-liquid while the viscosities went liquid-then-gas,
/// and since all four are double, a call site that got a pair the wrong way
/// round would have compiled in silence and returned wrong Reynolds numbers.
/// Nothing here can catch that; only the order being unsurprising can.
template <typename Source>
FlowScales flowScalesOf(const ClosureState &state, int cellIndex, const PhaseProperties &phases) {
    double gasRate = Source::gasFlowRate(state.cells, cellIndex) / phases.gasDensity;
    double liquidRate = Source::liquidFlowRate(state.cells, cellIndex) / phases.liquidDensity;
    double diameter = state.cells[cellIndex].duto.a;
    if (cellIndex > 0 && gasRate >= 0)
        diameter = state.cells[cellIndex - 1].duto.a;
    double flowArea = M_PI * diameter * diameter / 4.;

    double mixtureDensity = phases.noSlipLiquidHoldup * phases.liquidDensity + (1 - phases.noSlipLiquidHoldup) * phases.gasDensity;
    double mixtureViscosity = (phases.noSlipLiquidHoldup * phases.liquidViscosity + (1 - phases.noSlipLiquidHoldup) * phases.gasViscosity) / pow(10., 3.);
    double mixtureReynolds = diameter * mixtureDensity * (fabs(gasRate) / flowArea + fabs(liquidRate) / flowArea) / mixtureViscosity;
    double liquidReynolds = diameter * phases.liquidDensity * (fabs(gasRate) / flowArea + fabs(liquidRate) / flowArea) / (phases.liquidViscosity / 1000.);
    return FlowScales{gasRate, liquidRate, diameter, flowArea, mixtureReynolds, liquidReynolds};
}

/// Dispersed and stratified closure, evaluated as a pair. 167 tokens, identical
/// in all five.
///
/// mult0 and mult1 are dead -- assigned from ul0 and ul1 and never read, in
/// every variant. They are kept because they are part of the block that was
/// proven identical, and because deleting them would erase the only evidence
/// that a weighting was once intended here (A3-06).
void evaluateRegimePair(const ClosureState &state, int cellIndex, const MixtureProperties &mix,
                        double upstreamLiquidFlowRate, RegimePair &pair) {
    driftflux::correlations::C0UdDisperso(mix.liquidDensity, mix.gasDensity, mix.surfaceTension, mix.voidFraction, mix.mixtureReynolds, mix.liquidReynolds, mix.gasFlowRate, mix.liquidFlowRate, mix.diameter,
                                          state.cells[cellIndex].duto.rug, mix.inclinationAngle, pair.dispersedC0, pair.dispersedUd, mix.horizontalCorrection,
                                          state.cells[cellIndex].estabCol, state.selectors.dispersed);
    driftflux::correlations::C0UdEstratificado(mix.liquidDensity, mix.gasDensity, mix.surfaceTension, mix.voidFraction, mix.mixtureReynolds, mix.liquidReynolds, mix.gasFlowRate, mix.liquidFlowRate, mix.diameter,
                                               state.cells[cellIndex].duto.rug, mix.inclinationAngle, pair.stratifiedC0, pair.stratifiedUd, mix.horizontalCorrection,
                                               state.cells[cellIndex].estabCol, state.selectors.stratified);

    double mult0, mult1;
    mult0 = 1.;
    if (upstreamLiquidFlowRate < 0.)
        mult0 = 0.;
    mult1 = 0.;
    if (mix.liquidFlowRate < 0.)
        mult1 = 1.;
}

/// Blends the dispersed and stratified results by superficial velocity, with a
/// linear ramp between jmin and jmax. 128 tokens, identical in all five.
///
/// The ramp is written (1. - raz) * c0E + raz * c0D and MUST stay that way: the
/// algebraically equal c0D + (1. - raz) * (c0E - c0D) rounds differently.
void blendBySuperficialVelocity(const MixtureProperties &mix, const RegimePair &pair,
                                double &c0, double &ud) {
    double maxSuperficialVelocity = 0.05;
    double minSuperficialVelocity = 0.005;
    if ((fabs(mix.gasFlowRate) + fabs(mix.liquidFlowRate)) / mix.flowArea > maxSuperficialVelocity) {
        c0 = pair.stratifiedC0;
        ud = pair.stratifiedUd;
    } else if ((fabs(mix.gasFlowRate) + fabs(mix.liquidFlowRate)) / mix.flowArea < minSuperficialVelocity) {
        c0 = pair.dispersedC0;
        ud = pair.dispersedUd;
    } else {
        double blendRatio = (maxSuperficialVelocity - (fabs(mix.gasFlowRate) + fabs(mix.liquidFlowRate)) / mix.flowArea) / (maxSuperficialVelocity - minSuperficialVelocity);
        c0 = ((1. - blendRatio) * pair.stratifiedC0 + blendRatio * pair.dispersedC0);
        ud = ((1. - blendRatio) * pair.stratifiedUd + blendRatio * pair.dispersedUd);
    }
}

/// Dispersed closure, upgraded to annular/churn when the pattern says so.
/// 132 tokens, identical in all five.
void evaluateDispersedOrAnnular(const ClosureState &state, int cellIndex, const MixtureProperties &mix,
                                int flowPattern, double &c0, double &ud) {
    driftflux::correlations::C0UdDisperso(mix.liquidDensity, mix.gasDensity, mix.surfaceTension, mix.voidFraction, mix.mixtureReynolds, mix.liquidReynolds, mix.gasFlowRate, mix.liquidFlowRate, mix.diameter,
                                          state.cells[cellIndex].duto.rug, mix.inclinationAngle, c0, ud, mix.horizontalCorrection,
                                          state.cells[cellIndex].estabCol, state.selectors.dispersed);
    if (flowPattern == -2) {
        driftflux::correlations::C0UdAnularChurn(mix.liquidDensity, mix.gasDensity, mix.surfaceTension, mix.voidFraction, mix.mixtureReynolds, mix.liquidReynolds, mix.gasFlowRate, mix.liquidFlowRate,
                                                 mix.diameter, state.cells[cellIndex].duto.rug, mix.inclinationAngle, c0, ud,
                                                 mix.horizontalCorrection, state.cells[cellIndex].estabCol,
                                                 state.selectors.annularChurn);
    }
}

/// Discards the computed slip when the deck disables it. 113 tokens; the five
/// bodies carried two spellings of it, differing only in which configuration
/// field is read -- escorregaTran in the four transient variants, escorregaPerm
/// in the steady-state one -- and in writing the second zero of the last line as
/// `0.` in two of them and `0` in the other three.
///
/// The field is a parameter, which is what makes the two spellings one. The
/// zero is not a second difference to preserve: `0 * ud` converts the int to
/// 0.0 before multiplying, so it IS `0. * ud`, same operation and same rounding.
///
/// Everything above the two assignments is dead: correcaoUd and correcaoCo are
/// computed, clamped, and then multiplied by zero. Preserved, not removed --
/// c0 is forced to 1 and ud to 0 whenever slip is off, and the arithmetic that
/// says so is the record of what was once intended (A3-06).
void applyNoSlipOverride(const Cel *cells, int cellIndex, int slipEnabled, double &c0, double &ud) {
    if (slipEnabled == 0) {
        double meanSuperficialLiquidVelocity = cells[cellIndex].QL / cells[cellIndex].duto.area;
        double driftCorrection = 1 - (meanSuperficialLiquidVelocity - 0.15) / 0.35;
        double distributionCorrection = c0 - (c0 - 1) * (meanSuperficialLiquidVelocity - 0.15) / 0.35;
        if (driftCorrection > 1.)
            driftCorrection = 1.;
        if (driftCorrection < 0.)
            driftCorrection = 0.;
        if (distributionCorrection > c0)
            distributionCorrection = c0;
        if (distributionCorrection < 1)
            distributionCorrection = 1;
        c0 = 1. + 0 * distributionCorrection;
        ud = 0. + 0. * ud * driftCorrection;
    }
}

}  // namespace

void instantaneous(const ClosureState &state, int cellIndex, double &c0, double &ud) {
    int timeStep = 20;
    state.cells[cellIndex].transic0 = state.cells[cellIndex].transic;
    if (state.cells[cellIndex].dt < 0.8)
        timeStep *= (0.8 / state.cells[cellIndex].dt);
    c0 = 1.;
    ud = 0.;
    if (state.cells[cellIndex - 1].velPig > 0 && state.cells[cellIndex - 1].estadoPig == 1) {
        applyPigOverride(state, cellIndex, c0, ud);
    } else if (state.cells[cellIndex].velPig < 0 && state.cells[cellIndex].estadoPig == 1) {
        applyPigOverride(state, cellIndex, c0, ud);
    } else if ((state.cells[cellIndex].acsr.tipo != 4 || state.cells[cellIndex].acsr.bcs.freqnova <= 1.)) {
        double noSlipLiquidHoldup;
        double lengthRatio = state.cells[cellIndex].dxL / (state.cells[cellIndex].dx + state.cells[cellIndex].dxL);
        double upstreamLengthRatio;
        if (cellIndex > 0)
            upstreamLengthRatio = state.cells[cellIndex - 1].dxL / (state.cells[cellIndex - 1].dx + state.cells[cellIndex - 1].dxL);
        else
            upstreamLengthRatio = lengthRatio;

        noSlipLiquidHoldup = transientNoSlipHoldup<InstantaneousSource>(state, cellIndex);

        double liquidHoldup = noSlipLiquidHoldup;
        double voidFraction = 1 - liquidHoldup;
        double cellVoidFraction;
        cellVoidFraction = state.cells[cellIndex].alf;

        double alfneg;
        if (cellIndex > 1)
            alfneg = state.cells[cellIndex - 2].alf;
        else if (cellIndex > 0)
            alfneg = state.cells[cellIndex - 1].alf;
        else
            alfneg = state.cells[cellIndex].alf;

        double betI = state.cells[cellIndex].betL;
        if (cellIndex > 0)
            betI = state.cells[cellIndex - 1].betPigD;
        if (((0. * state.cells[cellIndex].QG + 1 * state.cells[cellIndex].QL) < 0.))
            betI = state.cells[cellIndex].betPigE; // duvidabeta

        double betneg;
        if (cellIndex > 0) {
            betneg = state.cells[cellIndex - 1].betL;
            if (cellIndex > 1)
                betneg = state.cells[cellIndex - 2].betPigD;
            if ((0.99 * state.cells[cellIndex - 1].QG + 0.01 * state.cells[cellIndex - 1].QL) < 0.)
                betneg = state.cells[cellIndex - 1].betPigE; // duvidabeta

        } else
            betneg = state.cells[cellIndex].bet;

        const MeanConditions conditions =
            instantaneousMeanConditions(state, cellIndex, lengthRatio, upstreamLengthRatio);
        const double meanPressure = conditions.pressure;
        const double meanTemperature = conditions.temperature;

        const double horizontalCorrection = horizontalCorrectionOf(state, cellIndex, cellIndex - 1);

        double surfaceTension;
        const PhaseProperties phases = instantaneousPhaseProperties(
            state, cellIndex, betI, meanPressure, meanTemperature, noSlipLiquidHoldup, surfaceTension);
        const double liquidDensity = phases.liquidDensity;
        const double gasDensity = phases.gasDensity;
        const double liquidViscosity = phases.liquidViscosity;
        const double gasViscosity = phases.gasViscosity;

        auto [gasFlowRate, liquidFlowRate, diameter, flowArea, mixtureReynolds, liquidReynolds] =
            flowScalesOf<InstantaneousSource>(state, cellIndex,
                                {.liquidDensity = liquidDensity,
                                 .gasDensity = gasDensity,
                                 .liquidViscosity = liquidViscosity,
                                 .gasViscosity = gasViscosity,
                                 .noSlipLiquidHoldup = noSlipLiquidHoldup});

        int flowPattern = 1;
        const double inclinationAngle = transientInclinationAngle<InstantaneousSource>(state, cellIndex);

        double transitionWindow = 20;

        const MixtureProperties mix{
            .liquidDensity = liquidDensity,
            .gasDensity = gasDensity,
            .surfaceTension = surfaceTension,
            .voidFraction = voidFraction,
            .gasFlowRate = gasFlowRate,
            .liquidFlowRate = liquidFlowRate,
            .diameter = diameter,
            .flowArea = flowArea,
            .mixtureReynolds = mixtureReynolds,
            .liquidReynolds = liquidReynolds,
            .inclinationAngle = inclinationAngle,
            .horizontalCorrection = horizontalCorrection,
        };
        if (mixtureReynolds > 1e-30) {
            if (fabs(0 * inclinationAngle + 1 * state.cells[cellIndex].duto.teta) < 45. * M_PI / 180. && liquidHoldup < 0.99 && liquidHoldup > 0.01 && cellIndex < state.lastCell - 1) {
                double upstreamGasFlowRate;
                double upstreamLiquidFlowRate;

                gasFlowRate = (state.cells[cellIndex].MC - state.cells[cellIndex].Mliqini) / gasDensity;
                liquidFlowRate = state.cells[cellIndex].Mliqini / liquidDensity;

                upstreamGasFlowRate = (state.cells[cellIndex - 1].MC - state.cells[cellIndex - 1].Mliqini) / state.cells[cellIndex].rgLi;
                // upstreamLiquidFlowRate = state.cells[cellIndex - 1].Mliqini
                //  / ((1 - betneg) * state.cells[cellIndex].flui.MasEspLiq(upstreamMeanPressure, upstreamMeanTemperature)
                upstreamLiquidFlowRate = state.cells[cellIndex - 1].Mliqini / ((1 - betneg) * state.cells[cellIndex].rpLi + betneg * state.cells[cellIndex].rcLi);

                estratificado stratifiedMap(diameter, liquidFlowRate, gasFlowRate, liquidDensity, gasDensity, liquidViscosity / pow(10., 3.), gasViscosity / pow(10., 3.), liquidHoldup,
                                        state.cells[cellIndex].duto.teta, state.cells[cellIndex].duto.rug / diameter);

                if (state.selectors.stratified == 2)
                    stratifiedMap.mapaTD();
                else
                    stratifiedMap.mapaTD(1);

                flowPattern = stratifiedMap.arr;

                if (flowPattern == -1) {
                    if (state.cells[cellIndex].arranjo != 0) {
                        if (((state.cells[cellIndex].arranjo != flowPattern) || state.cells[cellIndex].transic > 0)) {
                            if ((state.cells[cellIndex].arranjo != flowPattern) && state.cells[cellIndex].transic > 0)
                                state.cells[cellIndex].transic = 0;
                            state.cells[cellIndex].transic++;
                            if (state.cells[cellIndex].transic > transitionWindow - 1)
                                state.cells[cellIndex].transic = 0;
                        } else
                            state.cells[cellIndex].transic = 0;
                    }
                    state.cells[cellIndex].arranjo = flowPattern = stratifiedMap.arr;
                    state.cells[cellIndex - 1].arranjoR = stratifiedMap.arr;
                    state.cells[cellIndex - 1].perdaEstratL = stratifiedMap.fatorperdaLiq;
                    state.cells[cellIndex - 1].perdaEstratG = stratifiedMap.fatorperdaGas;
                    RegimePair pair;
                    evaluateRegimePair(state, cellIndex, mix, upstreamLiquidFlowRate, pair);
                    double alf0E = state.cells[cellIndex - 1].alf;

                    blendBySuperficialVelocity(mix, pair, c0, ud);

                    if (state.cells[cellIndex].transic > 0) {
                        c0 = (c0 * state.cells[cellIndex].transic + state.cells[cellIndex].c0 * (transitionWindow - state.cells[cellIndex].transic)) / transitionWindow;
                        ud = (ud * state.cells[cellIndex].transic + state.cells[cellIndex].ud * (transitionWindow - state.cells[cellIndex].transic)) / transitionWindow;
                    }
                }
            }
            if (flowPattern == 1) {

                arranjo flowPatternMap(diameter, liquidFlowRate / flowArea, gasFlowRate / flowArea, liquidDensity, gasDensity, liquidViscosity / pow(10., 3.), gasViscosity / pow(10., 3.), liquidHoldup,
                                   state.cells[cellIndex].duto.teta, surfaceTension, state.input.mapaArranjo, state.globals);
                flowPattern = flowPatternMap.verificaArr();

                evaluateDispersedOrAnnular(state, cellIndex, mix, flowPattern, c0, ud);
                if (fabs(gasFlowRate / state.cells[cellIndex].duto.area) > 5. && voidFraction >= 0.75) {
                    transitionWindow = 20;
                    if (state.selectors.annularChurn == 3 && state.selectors.dispersed == 1)
                        transitionWindow = 200;
                }

                if (state.cells[cellIndex].arranjo != 0) {
                    if ((flowPattern != state.cells[cellIndex].arranjo || state.cells[cellIndex].transic > 0)) {
                        if (flowPattern != state.cells[cellIndex].arranjo && state.cells[cellIndex].transic > 0)
                            state.cells[cellIndex].transic = 0;
                        state.cells[cellIndex].transic++;
                        if (state.cells[cellIndex].transic > transitionWindow - 1)
                            state.cells[cellIndex].transic = 0;
                    } else
                        state.cells[cellIndex].transic = 0;
                }
                if (state.cells[cellIndex].transic > 0) {
                    c0 = (c0 * state.cells[cellIndex].transic + state.cells[cellIndex].c0 * (transitionWindow - state.cells[cellIndex].transic)) / transitionWindow;
                    ud = (ud * state.cells[cellIndex].transic + state.cells[cellIndex].ud * (transitionWindow - state.cells[cellIndex].transic)) / transitionWindow;
                }
                state.cells[cellIndex].arranjo = flowPattern;
                state.cells[cellIndex - 1].arranjoR = flowPattern;
            }
            state.cells[cellIndex].c0Spare = c0;
            state.cells[cellIndex].udSpare = ud;
        }
    }
    applyNoSlipOverride(state.cells, cellIndex, state.input.escorregaTran, c0, ud);
}

void buffered(const ClosureState &state, int cellIndex, double &c0, double &ud) {
    int timeStep = 20;
    if (state.cells[cellIndex].dt < 0.8)
        timeStep *= (0.8 / state.cells[cellIndex].dt);
    c0 = 1.;
    ud = 0.;
    if (state.cells[cellIndex - 1].velPig > 0 && state.cells[cellIndex - 1].estadoPig == 1) {
        applyPigOverride(state, cellIndex, c0, ud);
    } else if (state.cells[cellIndex].velPig < 0 && state.cells[cellIndex].estadoPig == 1) {
        applyPigOverride(state, cellIndex, c0, ud);
    } else if ((state.cells[cellIndex].acsr.tipo != 4 || state.cells[cellIndex].acsr.bcs.freqnova <= 1.)) {
        double noSlipLiquidHoldup;
        double lengthRatio = state.cells[cellIndex].dxL / (state.cells[cellIndex].dx + state.cells[cellIndex].dxL);
        double upstreamLengthRatio;
        if (cellIndex > 0)
            upstreamLengthRatio = state.cells[cellIndex - 1].dxL / (state.cells[cellIndex - 1].dx + state.cells[cellIndex - 1].dxL);
        else
            upstreamLengthRatio = lengthRatio;

        noSlipLiquidHoldup = transientNoSlipHoldup<BufferedSource>(state, cellIndex);

        double liquidHoldup = noSlipLiquidHoldup;
        double voidFraction = 1 - liquidHoldup;
        double cellVoidFraction;
        cellVoidFraction = state.cells[cellIndex].alf;

        double alfneg;
        if (cellIndex > 1)
            alfneg = state.cells[cellIndex - 2].alf;
        else if (cellIndex > 0)
            alfneg = state.cells[cellIndex - 1].alf;
        else
            alfneg = state.cells[cellIndex].alf;

        double betI = state.cells[cellIndex].betL;
        if (cellIndex > 0)
            betI = state.cells[cellIndex - 1].betPigD;
        if ((state.cells[cellIndex].MliqiniBuf) < 0.)
            betI = state.cells[cellIndex].betPigE; // testeBeta
        betI = state.cells[cellIndex].betPigE;     // duvidabeta
        double betneg;
        if (cellIndex > 0) {
            betneg = state.cells[cellIndex - 1].betL;
            if (cellIndex > 1)
                betneg = state.cells[cellIndex - 2].betPigD;
            if (state.cells[cellIndex].MliqiniLBuf < 0.)
                betneg = state.cells[cellIndex - 1].betPigE; // testeBeta
            betneg = state.cells[cellIndex - 1].betPigE;     // duvidabeta
        } else
            betneg = state.cells[cellIndex].bet;

        double meanPressure;
        double upstreamMeanPressure = 0.;

        meanPressure = lengthRatio * state.cells[cellIndex].presBuf + (1 - lengthRatio) * state.cells[cellIndex].presLBuf;
        if (cellIndex == state.lastCell)
            meanPressure = state.cells[cellIndex].presBuf;
        upstreamMeanPressure = state.cells[cellIndex].presauxL;
        double meanTemperature = lengthRatio * state.cells[cellIndex].temp + (1 - lengthRatio) * state.cells[cellIndex].tempL;
        meanTemperature = state.cells[cellIndex].tempL;
        if (state.cells[cellIndex].VTemper < 0.) {
            if (cellIndex < state.lastCell)
                meanTemperature = state.cells[cellIndex].temp;
            else
                meanTemperature = state.gasSurfaceTemperature;
        }
        double upstreamMeanTemperature;
        if (cellIndex > 0)
            upstreamMeanTemperature = upstreamLengthRatio * state.cells[cellIndex - 1].temp + (1 - upstreamLengthRatio) * state.cells[cellIndex - 1].tempL;
        else
            upstreamMeanTemperature = meanTemperature;

        const double horizontalCorrection = horizontalCorrectionOf(state, cellIndex, cellIndex);

        double liquidDensity;
        double liquidViscosity;
        double surfaceTension;
        if ((state.cells[cellIndex].MliqiniBuf) < 0.) { // testeBeta
            liquidDensity = (1 - betI) * state.cells[cellIndex].flui.MasEspLiq(meanPressure, meanTemperature) + betI * state.cells[cellIndex].fluicol.MasEspFlu(meanPressure, meanTemperature);
            liquidViscosity = (1 - betI) * state.cells[cellIndex].flui.ViscOleo(meanPressure, meanTemperature) + betI * state.cells[cellIndex].fluicol.VisFlu(meanPressure, meanTemperature);
            surfaceTension = (1 - betI) * state.cells[cellIndex].flui.TensSuper(meanPressure, meanTemperature) + betI * state.cells[cellIndex].fluicol.TensSuper(meanPressure, meanTemperature);
        } else {
            liquidDensity = (1 - betI) * state.cells[cellIndex - 1].flui.MasEspLiq(meanPressure, meanTemperature) + betI * state.cells[cellIndex - 1].fluicol.MasEspFlu(meanPressure, meanTemperature);
            liquidViscosity = (1 - betI) * state.cells[cellIndex - 1].flui.ViscOleo(meanPressure, meanTemperature) + betI * state.cells[cellIndex - 1].fluicol.VisFlu(meanPressure, meanTemperature);
            surfaceTension = (1 - betI) * state.cells[cellIndex - 1].flui.TensSuper(meanPressure, meanTemperature) + betI * state.cells[cellIndex - 1].fluicol.TensSuper(meanPressure, meanTemperature);
        }

        double gasDensity;
        double gasViscosity;
        if ((state.cells[cellIndex].MCBuf - state.cells[cellIndex].MliqiniBuf) < 0.) {
            gasDensity = state.cells[cellIndex].flui.MasEspGas(meanPressure, meanTemperature);
            gasViscosity = state.cells[cellIndex].flui.ViscGas(meanPressure, meanTemperature);
        } else {
            gasDensity = state.cells[cellIndex - 1].flui.MasEspGas(meanPressure, meanTemperature);
            gasViscosity = state.cells[cellIndex - 1].flui.ViscGas(meanPressure, meanTemperature);
        }

        auto [gasFlowRate, liquidFlowRate, diameter, flowArea, mixtureReynolds, liquidReynolds] =
            flowScalesOf<BufferedSource>(state, cellIndex,
                                {.liquidDensity = liquidDensity,
                                 .gasDensity = gasDensity,
                                 .liquidViscosity = liquidViscosity,
                                 .gasViscosity = gasViscosity,
                                 .noSlipLiquidHoldup = noSlipLiquidHoldup});

        int flowPattern = 1;
        const double inclinationAngle = transientInclinationAngle<BufferedSource>(state, cellIndex);
        double transitionWindow = 20;
        const MixtureProperties mix{
            .liquidDensity = liquidDensity,
            .gasDensity = gasDensity,
            .surfaceTension = surfaceTension,
            .voidFraction = voidFraction,
            .gasFlowRate = gasFlowRate,
            .liquidFlowRate = liquidFlowRate,
            .diameter = diameter,
            .flowArea = flowArea,
            .mixtureReynolds = mixtureReynolds,
            .liquidReynolds = liquidReynolds,
            .inclinationAngle = inclinationAngle,
            .horizontalCorrection = horizontalCorrection,
        };
        if (mixtureReynolds > 1e-30) {
            if (fabs(0 * inclinationAngle + 1 * state.cells[cellIndex].duto.teta) < 45. * M_PI / 180. && liquidHoldup < 0.99 && liquidHoldup > 0.01 && cellIndex < state.lastCell - 1) {
                double upstreamGasFlowRate;
                double upstreamLiquidFlowRate;

                gasFlowRate = (state.cells[cellIndex].MCBuf - state.cells[cellIndex].MliqiniBuf) / gasDensity;
                liquidFlowRate = (state.cells[cellIndex].MliqiniBuf) / liquidDensity;

                upstreamGasFlowRate = (state.cells[cellIndex - 1].MCBuf - state.cells[cellIndex - 1].MliqiniBuf) / state.cells[cellIndex].flui.MasEspGas(upstreamMeanPressure, upstreamMeanTemperature);
                upstreamLiquidFlowRate = (state.cells[cellIndex - 1].MliqiniBuf) / ((1 - betneg) * state.cells[cellIndex].flui.MasEspLiq(upstreamMeanPressure, upstreamMeanTemperature) + betneg * state.cells[cellIndex].fluicol.MasEspFlu(upstreamMeanPressure, upstreamMeanTemperature));

                flowPattern = state.cells[cellIndex].arranjo;
                if (flowPattern == -1) {

                    RegimePair pair;
                    evaluateRegimePair(state, cellIndex, mix, upstreamLiquidFlowRate, pair);
                    double alf0E = state.cells[cellIndex - 1].alf;

                    blendBySuperficialVelocity(mix, pair, c0, ud);

                    if (state.cells[cellIndex].transic > 0) {
                        c0 = (c0 * state.cells[cellIndex].transic + state.cells[cellIndex].c0 * (transitionWindow - state.cells[cellIndex].transic)) / transitionWindow;
                        ud = (ud * state.cells[cellIndex].transic + state.cells[cellIndex].ud * (transitionWindow - state.cells[cellIndex].transic)) / transitionWindow;
                    }
                }
            }
            if (flowPattern != -1) {

                evaluateDispersedOrAnnular(state, cellIndex, mix, flowPattern, c0, ud);
                if (fabs(gasFlowRate / state.cells[cellIndex].duto.area) > 5. && voidFraction >= 0.75) {
                    transitionWindow = 20;
                    if (state.selectors.annularChurn == 3 && state.selectors.dispersed == 1)
                        transitionWindow = 200;
                }

                if (state.cells[cellIndex].transic > 0) {
                    c0 = (c0 * state.cells[cellIndex].transic + state.cells[cellIndex].c0 * (transitionWindow - state.cells[cellIndex].transic)) / transitionWindow;
                    ud = (ud * state.cells[cellIndex].transic + state.cells[cellIndex].ud * (transitionWindow - state.cells[cellIndex].transic)) / transitionWindow;
                }
                state.cells[cellIndex].arranjo = flowPattern;
                state.cells[cellIndex - 1].arranjoR = flowPattern;
            }
            state.cells[cellIndex].c0Spare = c0;
            state.cells[cellIndex].udSpare = ud;
        }
    }
    applyNoSlipOverride(state.cells, cellIndex, state.input.escorregaTran, c0, ud);
}

void initialization(const ClosureState &state, int cellIndex, double &c0, double &ud) {
    int timeStep = 20;
    if (state.cells[cellIndex].dt < 0.8)
        timeStep *= (0.8 / state.cells[cellIndex].dt);
    c0 = 1.;
    ud = 0.;
    if (state.cells[cellIndex].velPig < 0 && state.cells[cellIndex].estadoPig == 1) {
        c0 = 1.;
        ud = 0.;
        state.cells[cellIndex].arranjo = 1.;

    } else if ((state.cells[cellIndex].acsr.tipo != 4 || state.cells[cellIndex].acsr.bcs.freqnova <= 1.)) {
        double noSlipLiquidHoldup;

        noSlipLiquidHoldup = inletNoSlipHoldup<InstantaneousSource>(state, cellIndex);

        double liquidHoldup = noSlipLiquidHoldup;
        double voidFraction = 1 - liquidHoldup;
        double cellVoidFraction;
        cellVoidFraction = state.cells[cellIndex].alf;

        double alfneg;
        if (cellIndex > 1)
            alfneg = state.inletVoidFraction;
        else if (cellIndex > 0)
            alfneg = state.inletVoidFraction;
        else
            alfneg = state.cells[cellIndex].alf;

        double betI = state.cells[cellIndex].betL;
        if (cellIndex > 0)
            betI = state.inletColumnFraction;
        if (state.cells[cellIndex].QL < 0.)
            betI = state.cells[cellIndex].betPigE; // testeBeta
        betI = state.cells[cellIndex].betPigE;     // duvidabeta
        double betneg;
        if (cellIndex > 0) {
            betneg = state.inletColumnFraction;

        } else
            betneg = state.cells[cellIndex].bet;

        double meanPressure;
        double upstreamMeanPressure = 0.;

        meanPressure = state.inletPressure;
        if (cellIndex > 0)
            upstreamMeanPressure = state.inletPressure;
        if (cellIndex == state.lastCell)
            meanPressure = state.cells[cellIndex].pres;
        else
            upstreamMeanPressure = state.inletPressure;
        double meanTemperature = state.inletTemperature;

        double upstreamMeanTemperature;
        upstreamMeanTemperature = meanTemperature;

        const double horizontalCorrection = horizontalCorrectionOf(state, cellIndex, cellIndex);

        double liquidDensity;
        double liquidViscosity;
        double surfaceTension;
        if (state.cells[cellIndex].QL < 0.) { // testeBeta
            liquidDensity = (1 - betI) * state.cells[cellIndex].flui.MasEspLiq(meanPressure, meanTemperature) + betI * state.cells[cellIndex].fluicol.MasEspFlu(meanPressure, meanTemperature);
            liquidViscosity = (1 - betI) * state.cells[cellIndex].flui.ViscOleo(meanPressure, meanTemperature) + betI * state.cells[cellIndex].fluicol.VisFlu(meanPressure, meanTemperature);
            surfaceTension = (1 - betI) * state.cells[cellIndex].flui.TensSuper(meanPressure, meanTemperature) + betI * state.cells[cellIndex].fluicol.TensSuper(meanPressure, meanTemperature);
        } else {
            liquidDensity = (1 - betI) * (*state.cells[cellIndex].fluiL).MasEspLiq(meanPressure, meanTemperature) + betI * state.cells[cellIndex].fluicol.MasEspFlu(meanPressure, meanTemperature);
            liquidViscosity = (1 - betI) * (*state.cells[cellIndex].fluiL).ViscOleo(meanPressure, meanTemperature) + betI * state.cells[cellIndex].fluicol.VisFlu(meanPressure, meanTemperature);
            surfaceTension = (1 - betI) * (*state.cells[cellIndex].fluiL).TensSuper(meanPressure, meanTemperature) + betI * state.cells[cellIndex].fluicol.TensSuper(meanPressure, meanTemperature);
        }

        double gasDensity;
        double gasViscosity;
        if (state.cells[cellIndex].QG < 0.) {
            gasDensity = state.cells[cellIndex].flui.MasEspGas(meanPressure, meanTemperature);
            gasViscosity = state.cells[cellIndex].flui.ViscGas(meanPressure, meanTemperature);
        } else {
            gasDensity = (*state.cells[cellIndex].fluiL).MasEspGas(meanPressure, meanTemperature);
            gasViscosity = (*state.cells[cellIndex].fluiL).ViscGas(meanPressure, meanTemperature);
        }

        auto [gasFlowRate, liquidFlowRate, diameter, flowArea, mixtureReynolds, liquidReynolds] =
            flowScalesOf<InstantaneousSource>(state, cellIndex,
                                {.liquidDensity = liquidDensity,
                                 .gasDensity = gasDensity,
                                 .liquidViscosity = liquidViscosity,
                                 .gasViscosity = gasViscosity,
                                 .noSlipLiquidHoldup = noSlipLiquidHoldup});

        int flowPattern = 1;
        double totalLength = state.cells[cellIndex].dxL + state.cells[cellIndex].dx;
        double leftCellLength = state.cells[cellIndex].dxL;
        double cellLength = state.cells[cellIndex].dx;
        double inclinationAngle = (cellLength * state.cells[cellIndex].dutoL.teta + leftCellLength * state.cells[cellIndex].duto.teta) / totalLength;
        double transitionWindow = 20.;
        const MixtureProperties mix{
            .liquidDensity = liquidDensity,
            .gasDensity = gasDensity,
            .surfaceTension = surfaceTension,
            .voidFraction = voidFraction,
            .gasFlowRate = gasFlowRate,
            .liquidFlowRate = liquidFlowRate,
            .diameter = diameter,
            .flowArea = flowArea,
            .mixtureReynolds = mixtureReynolds,
            .liquidReynolds = liquidReynolds,
            .inclinationAngle = inclinationAngle,
            .horizontalCorrection = horizontalCorrection,
        };
        if (mixtureReynolds > 1e-30) {
            if (fabs(0 * inclinationAngle + 1 * state.cells[cellIndex].duto.teta) < 45. * M_PI / 180. && liquidHoldup < 0.99 && liquidHoldup > 0.01 && cellIndex < state.lastCell - 1) {
                double upstreamGasFlowRate;
                double upstreamLiquidFlowRate;

                gasFlowRate = (state.cells[cellIndex].MC - state.cells[cellIndex].Mliqini) / gasDensity;
                liquidFlowRate = state.cells[cellIndex].Mliqini / liquidDensity;

                upstreamGasFlowRate = (state.cells[cellIndex].MC - state.cells[cellIndex].Mliqini) / (*state.cells[cellIndex].fluiL).MasEspGas(upstreamMeanPressure, upstreamMeanTemperature);
                upstreamLiquidFlowRate = state.cells[cellIndex].Mliqini / ((1 - betneg) * (*state.cells[cellIndex].fluiL).MasEspLiq(upstreamMeanPressure, upstreamMeanTemperature) + betneg * state.cells[cellIndex].fluicol.MasEspFlu(upstreamMeanPressure, upstreamMeanTemperature));

                estratificado stratifiedMap(diameter, liquidFlowRate, gasFlowRate, liquidDensity, gasDensity, liquidViscosity / pow(10., 3.), gasViscosity / pow(10., 3.), liquidHoldup,
                                        state.cells[cellIndex].duto.teta, state.cells[cellIndex].duto.rug / diameter);

                stratifiedMap.mapaTD();
                flowPattern = stratifiedMap.arr;
                if (flowPattern == -1) {
                    if (state.cells[cellIndex].arranjo != 0) {
                        if (((state.cells[cellIndex].arranjo != flowPattern) || state.cells[cellIndex].transic > 0)) {
                            if ((state.cells[cellIndex].arranjo != flowPattern) && state.cells[cellIndex].transic > 0)
                                state.cells[cellIndex].transic = 0;
                            state.cells[cellIndex].transic++;
                            if (state.cells[cellIndex].transic > 19)
                                state.cells[cellIndex].transic = 0;
                        }
                    } else
                        state.cells[cellIndex].transic = 0;
                    state.cells[cellIndex].arranjo = flowPattern = stratifiedMap.arr;
                    RegimePair pair;
                    evaluateRegimePair(state, cellIndex, mix, upstreamLiquidFlowRate, pair);
                    double alf0E = state.inletVoidFraction;

                    blendBySuperficialVelocity(mix, pair, c0, ud);

                    if (state.cells[cellIndex].transic > 0) {
                        c0 = (c0 * state.cells[cellIndex].transic + state.cells[cellIndex].c0 * (transitionWindow - state.cells[cellIndex].transic)) / transitionWindow;
                        ud = (ud * state.cells[cellIndex].transic + state.cells[cellIndex].ud * (transitionWindow - state.cells[cellIndex].transic)) / transitionWindow;
                    }
                }
            }
            if (flowPattern == 1) {

                arranjo flowPatternMap(diameter, liquidFlowRate / flowArea, gasFlowRate / flowArea, liquidDensity, gasDensity, liquidViscosity / pow(10., 3.), gasViscosity / pow(10., 3.), liquidHoldup,
                                   state.cells[cellIndex].duto.teta, surfaceTension, state.input.mapaArranjo, state.globals);
                flowPattern = flowPatternMap.verificaArr();

                evaluateDispersedOrAnnular(state, cellIndex, mix, flowPattern, c0, ud);

                if (fabs(gasFlowRate / state.cells[cellIndex].duto.area) > 5. && voidFraction >= 0.75) {
                    transitionWindow = 20;
                    if (state.selectors.annularChurn == 3 && state.selectors.dispersed == 1)
                        transitionWindow = 200;
                }
                if (state.cells[cellIndex].arranjo != 0) {
                    if ((flowPattern != state.cells[cellIndex].arranjo || state.cells[cellIndex].transic > 0)) {
                        if (flowPattern != state.cells[cellIndex].arranjo && state.cells[cellIndex].transic > 0)
                            state.cells[cellIndex].transic = 0;
                        state.cells[cellIndex].transic++;
                        if (state.cells[cellIndex].transic > transitionWindow - 1)
                            state.cells[cellIndex].transic = 0;
                    }
                } else
                    state.cells[cellIndex].transic = 0;
                if (state.cells[cellIndex].transic > 0) {
                    c0 = (c0 * state.cells[cellIndex].transic + state.cells[cellIndex].c0 * (transitionWindow - state.cells[cellIndex].transic)) / transitionWindow;
                    ud = (ud * state.cells[cellIndex].transic + state.cells[cellIndex].ud * (transitionWindow - state.cells[cellIndex].transic)) / transitionWindow;
                }
                state.cells[cellIndex].arranjo = flowPattern;
            }
            state.cells[cellIndex].c0Spare = c0;
            state.cells[cellIndex].udSpare = ud;
        }
    }
    applyNoSlipOverride(state.cells, cellIndex, state.input.escorregaTran, c0, ud);
}

void bufferedInitialization(const ClosureState &state, int cellIndex, double &c0, double &ud) {
    int timeStep = 20;
    if (state.cells[cellIndex].dt < 0.8)
        timeStep *= (0.8 / state.cells[cellIndex].dt);
    c0 = 1.;
    ud = 0.;
    if (state.cells[cellIndex].velPig < 0 && state.cells[cellIndex].estadoPig == 1) {
        c0 = 1.;
        ud = 0.;
        state.cells[cellIndex].arranjo = 1.;

    } else if ((state.cells[cellIndex].acsr.tipo != 4 || state.cells[cellIndex].acsr.bcs.freqnova <= 1.)) {
        double noSlipLiquidHoldup;

        noSlipLiquidHoldup = inletNoSlipHoldup<BufferedSource>(state, cellIndex);

        double liquidHoldup = noSlipLiquidHoldup;
        double voidFraction = 1 - liquidHoldup;
        double cellVoidFraction;
        cellVoidFraction = state.cells[cellIndex].alf;

        double alfneg;
        if (cellIndex > 1)
            alfneg = state.inletVoidFraction;
        else if (cellIndex > 0)
            alfneg = state.inletVoidFraction;
        else
            alfneg = state.cells[cellIndex].alf;

        double betI = state.cells[cellIndex].betL;
        if (cellIndex > 0)
            betI = state.inletColumnFraction;
        if (state.cells[cellIndex].QL < 0.)
            betI = state.cells[cellIndex].betPigE; // testeBeta
        betI = state.cells[cellIndex].betPigE;     // duvidabeta
        double betneg;
        if (cellIndex > 0) {
            betneg = state.inletColumnFraction;

        } else
            betneg = state.cells[cellIndex].bet;

        double meanPressure;
        double upstreamMeanPressure = 0.;

        meanPressure = state.inletPressure;
        if (cellIndex > 0)
            upstreamMeanPressure = state.inletPressure;
        else
            upstreamMeanPressure = state.inletPressure;
        double meanTemperature = state.inletTemperature;

        double upstreamMeanTemperature;
        upstreamMeanTemperature = meanTemperature;

        const double horizontalCorrection = horizontalCorrectionOf(state, cellIndex, cellIndex);

        double liquidDensity;
        double liquidViscosity;
        double surfaceTension;
        if (state.cells[cellIndex].MliqiniBuf < 0.) { // testeBeta
            liquidDensity = (1 - betI) * state.cells[cellIndex].flui.MasEspLiq(meanPressure, meanTemperature) + betI * state.cells[cellIndex].fluicol.MasEspFlu(meanPressure, meanTemperature);
            liquidViscosity = (1 - betI) * state.cells[cellIndex].flui.ViscOleo(meanPressure, meanTemperature) + betI * state.cells[cellIndex].fluicol.VisFlu(meanPressure, meanTemperature);
            surfaceTension = (1 - betI) * state.cells[cellIndex].flui.TensSuper(meanPressure, meanTemperature) + betI * state.cells[cellIndex].fluicol.TensSuper(meanPressure, meanTemperature);
        } else {
            liquidDensity = (1 - betI) * (*state.cells[cellIndex].fluiL).MasEspLiq(meanPressure, meanTemperature) + betI * state.cells[cellIndex].fluicol.MasEspFlu(meanPressure, meanTemperature);
            liquidViscosity = (1 - betI) * (*state.cells[cellIndex].fluiL).ViscOleo(meanPressure, meanTemperature) + betI * state.cells[cellIndex].fluicol.VisFlu(meanPressure, meanTemperature);
            surfaceTension = (1 - betI) * (*state.cells[cellIndex].fluiL).TensSuper(meanPressure, meanTemperature) + betI * state.cells[cellIndex].fluicol.TensSuper(meanPressure, meanTemperature);
        }

        double gasDensity;
        double gasViscosity;
        if ((state.cells[cellIndex].MCBuf - state.cells[cellIndex].MliqiniBuf) < 0.) {
            gasDensity = state.cells[cellIndex].flui.MasEspGas(meanPressure, meanTemperature);
            gasViscosity = state.cells[cellIndex].flui.ViscGas(meanPressure, meanTemperature);
        } else {
            gasDensity = (*state.cells[cellIndex].fluiL).MasEspGas(meanPressure, meanTemperature);
            gasViscosity = (*state.cells[cellIndex].fluiL).ViscGas(meanPressure, meanTemperature);
        }

        auto [gasFlowRate, liquidFlowRate, diameter, flowArea, mixtureReynolds, liquidReynolds] =
            flowScalesOf<BufferedSource>(state, cellIndex,
                                {.liquidDensity = liquidDensity,
                                 .gasDensity = gasDensity,
                                 .liquidViscosity = liquidViscosity,
                                 .gasViscosity = gasViscosity,
                                 .noSlipLiquidHoldup = noSlipLiquidHoldup});

        int flowPattern = 1;
        double totalLength = state.cells[cellIndex].dxL + state.cells[cellIndex].dx;
        double leftCellLength = state.cells[cellIndex].dxL;
        double cellLength = state.cells[cellIndex].dx;
        double inclinationAngle = (cellLength * state.cells[cellIndex].dutoL.teta + leftCellLength * state.cells[cellIndex].duto.teta) / totalLength;
        double transitionWindow = 20.;
        const MixtureProperties mix{
            .liquidDensity = liquidDensity,
            .gasDensity = gasDensity,
            .surfaceTension = surfaceTension,
            .voidFraction = voidFraction,
            .gasFlowRate = gasFlowRate,
            .liquidFlowRate = liquidFlowRate,
            .diameter = diameter,
            .flowArea = flowArea,
            .mixtureReynolds = mixtureReynolds,
            .liquidReynolds = liquidReynolds,
            .inclinationAngle = inclinationAngle,
            .horizontalCorrection = horizontalCorrection,
        };
        if (mixtureReynolds > 1e-30) {
            if (fabs(0 * inclinationAngle + 1 * state.cells[cellIndex].duto.teta) < 45. * M_PI / 180. && liquidHoldup < 0.99 && liquidHoldup > 0.01 && cellIndex < state.lastCell - 1) {
                double upstreamGasFlowRate;
                double upstreamLiquidFlowRate;

                gasFlowRate = (state.cells[cellIndex].MCBuf - state.cells[cellIndex].MliqiniBuf) / gasDensity;
                liquidFlowRate = state.cells[cellIndex].MliqiniBuf / liquidDensity;

                upstreamGasFlowRate = (state.cells[cellIndex].MCBuf - state.cells[cellIndex].MliqiniBuf) / (*state.cells[cellIndex].fluiL).MasEspGas(upstreamMeanPressure, upstreamMeanTemperature);
                upstreamLiquidFlowRate = state.cells[cellIndex].MliqiniBuf / ((1 - betneg) * (*state.cells[cellIndex].fluiL).MasEspLiq(upstreamMeanPressure, upstreamMeanTemperature) + betneg * state.cells[cellIndex].fluicol.MasEspFlu(upstreamMeanPressure, upstreamMeanTemperature));

                flowPattern = state.cells[cellIndex].arranjo;
                if (flowPattern == -1) {

                    RegimePair pair;
                    evaluateRegimePair(state, cellIndex, mix, upstreamLiquidFlowRate, pair);

                    blendBySuperficialVelocity(mix, pair, c0, ud);

                    if (state.cells[cellIndex].transic > 0) {
                        c0 = (c0 * state.cells[cellIndex].transic + state.cells[cellIndex].c0 * (transitionWindow - state.cells[cellIndex].transic)) / transitionWindow;
                        ud = (ud * state.cells[cellIndex].transic + state.cells[cellIndex].ud * (transitionWindow - state.cells[cellIndex].transic)) / transitionWindow;
                    }
                }
            }
            if (flowPattern != -1) {

                evaluateDispersedOrAnnular(state, cellIndex, mix, flowPattern, c0, ud);
                if (fabs(gasFlowRate / state.cells[cellIndex].duto.area) > 5. && voidFraction >= 0.75) {
                    transitionWindow = 20;
                    if (state.selectors.annularChurn == 3 && state.selectors.dispersed == 1)
                        transitionWindow = 200;
                }

                if (state.cells[cellIndex].transic > 0) {
                    c0 = (c0 * state.cells[cellIndex].transic + state.cells[cellIndex].c0 * (transitionWindow - state.cells[cellIndex].transic)) / transitionWindow;
                    ud = (ud * state.cells[cellIndex].transic + state.cells[cellIndex].ud * (transitionWindow - state.cells[cellIndex].transic)) / transitionWindow;
                }
                state.cells[cellIndex].arranjo = flowPattern;
            }
            state.cells[cellIndex].c0Spare = c0;
            state.cells[cellIndex].udSpare = ud;
        }
    }
    applyNoSlipOverride(state.cells, cellIndex, state.input.escorregaTran, c0, ud);
}

void steadyState(const ClosureState &state, int cellIndex, double &c0, double &ud) {

    c0 = 1.;
    ud = 0.;
    if (state.cells[cellIndex].acsr.tipo != 4 || state.cells[cellIndex].acsr.bcs.freqnova <= 1.) {
        double noSlipLiquidHoldup;
        double lengthRatio = state.cells[cellIndex].dx / (state.cells[cellIndex].dx + state.cells[cellIndex].dxL);
        double upstreamLengthRatio;
        if (cellIndex > 0)
            upstreamLengthRatio = state.cells[cellIndex - 1].dx / (state.cells[cellIndex - 1].dx + state.cells[cellIndex - 1].dxL);
        else
            upstreamLengthRatio = lengthRatio;
        noSlipLiquidHoldup = 1. - state.cells[cellIndex].alf;

        double liquidHoldup = noSlipLiquidHoldup;
        double voidFraction = 1 - liquidHoldup;
        double cellVoidFraction;
        cellVoidFraction = state.cells[cellIndex].alf;

        double alfneg;
        if (cellIndex > 1)
            alfneg = state.cells[cellIndex - 2].alf;
        else if (cellIndex > 0)
            alfneg = state.cells[cellIndex - 1].alf;
        else
            alfneg = state.cells[cellIndex].alf;

        double betI = state.cells[cellIndex].betL;
        double betneg = state.cells[cellIndex].betL;

        double meanPressure;
        double upstreamMeanPressure = 0.;

        meanPressure = state.cells[cellIndex].presaux;
        if (cellIndex > 0)
            upstreamMeanPressure = state.cells[cellIndex - 1].presaux;
        else
            upstreamMeanPressure = state.cells[cellIndex].presaux;
        double meanTemperature;
        if (state.steadyIteration != 0 && state.input.AceleraConvergPerm == 0)
            meanTemperature = lengthRatio * state.cells[cellIndex].temp + (1 - lengthRatio) * state.cells[cellIndex].tempL;
        else
            meanTemperature = state.cells[cellIndex - 1].temp;
        double upstreamMeanTemperature;
        if (cellIndex > 0 && state.input.AceleraConvergPerm == 0)
            upstreamMeanTemperature = upstreamLengthRatio * state.cells[cellIndex - 1].temp + (1 - upstreamLengthRatio) * state.cells[cellIndex - 1].tempL;
        else
            upstreamMeanTemperature = meanTemperature;

        double horizontalCorrection = 1.;
        if (fabs(state.cells[cellIndex].duto.teta) < 1e-10) {
            if (state.cells[cellIndex].angEsq < 0 && state.cells[cellIndex].angDir < 0)
                horizontalCorrection = -1.;
            else if (state.cells[cellIndex].angEsq > 0 && state.cells[cellIndex].angDir > 0)
                horizontalCorrection = 1.;
        }

        double liquidDensity;
        double liquidViscosity;
        double surfaceTension;
        if (cellIndex > 0)
            liquidDensity = (1 - betI) * state.cells[cellIndex - 1].flui.MasEspLiq(upstreamMeanPressure, upstreamMeanTemperature) + betI * state.cells[cellIndex - 1].fluicol.MasEspFlu(upstreamMeanPressure, upstreamMeanTemperature);
        else
            liquidDensity = (1 - betI) * (*state.cells[cellIndex].fluiL).MasEspLiq(meanPressure, meanTemperature) + betI * state.cells[cellIndex].fluicol.MasEspFlu(meanPressure, meanTemperature);

        liquidViscosity = (1 - betI) * (*state.cells[cellIndex].fluiL).ViscOleo(upstreamMeanPressure, upstreamMeanTemperature) + betI * state.cells[cellIndex].fluicol.VisFlu(upstreamMeanPressure, upstreamMeanTemperature);
        surfaceTension = (1 - betI) * (*state.cells[cellIndex].fluiL).TensSuper(meanPressure, meanTemperature) + betI * state.cells[cellIndex].fluicol.TensSuper(meanPressure, meanTemperature);

        double gasDensity;
        double gasViscosity;
        if (cellIndex > 0)
            gasDensity = state.cells[cellIndex - 1].flui.MasEspGas(upstreamMeanPressure, upstreamMeanTemperature);
        else
            gasDensity = (*state.cells[cellIndex].fluiL).MasEspGas(meanPressure, meanTemperature);
        gasViscosity = (*state.cells[cellIndex].fluiL).ViscGas(meanPressure, meanTemperature);

        auto [gasFlowRate, liquidFlowRate, diameter, flowArea, mixtureReynolds, liquidReynolds] =
            flowScalesOf<SteadyStateSource>(state, cellIndex,
                                {.liquidDensity = liquidDensity,
                                 .gasDensity = gasDensity,
                                 .liquidViscosity = liquidViscosity,
                                 .gasViscosity = gasViscosity,
                                 .noSlipLiquidHoldup = noSlipLiquidHoldup});

        int flowPattern = 1;
        double totalLength = state.cells[cellIndex].dxL + state.cells[cellIndex].dx;
        double leftCellLength = state.cells[cellIndex].dxL;
        double cellLength = state.cells[cellIndex].dx;
        double inclinationAngle = (leftCellLength * state.cells[cellIndex].dutoL.teta + cellLength * state.cells[cellIndex].duto.teta) / totalLength;
        double inclinationSign = 1.;
        if (fabs(state.cells[cellIndex].MC) > 1e-15)
            inclinationSign = state.cells[cellIndex].MC / fabs(state.cells[cellIndex].MC);
        if (gasDensity < 0.9 * liquidDensity) {
        const MixtureProperties mix{
            .liquidDensity = liquidDensity,
            .gasDensity = gasDensity,
            .surfaceTension = surfaceTension,
            .voidFraction = voidFraction,
            .gasFlowRate = gasFlowRate,
            .liquidFlowRate = liquidFlowRate,
            .diameter = diameter,
            .flowArea = flowArea,
            .mixtureReynolds = mixtureReynolds,
            .liquidReynolds = liquidReynolds,
            .inclinationAngle = inclinationAngle,
            .horizontalCorrection = horizontalCorrection,
        };
            if (mixtureReynolds > 1e-30) {
                if (fabs(0 * inclinationAngle + inclinationSign * state.cells[cellIndex].duto.teta) < 45. * M_PI / 180. && liquidHoldup < 0.99 && liquidHoldup > 0.01 && cellIndex < state.lastCell - 1) {
                    double upstreamGasFlowRate;
                    double upstreamLiquidFlowRate;

                    gasFlowRate = fabs(state.cells[cellIndex].MC - state.cells[cellIndex].Mliqini) / gasDensity;
                    liquidFlowRate = fabs(state.cells[cellIndex].Mliqini) / liquidDensity;

                    upstreamGasFlowRate = gasFlowRate;
                    upstreamLiquidFlowRate = liquidFlowRate;
                    if (cellIndex > 0) {
                        upstreamGasFlowRate = fabs(state.cells[cellIndex - 1].MC - state.cells[cellIndex - 1].Mliqini) / state.cells[cellIndex].flui.MasEspGas(upstreamMeanPressure, upstreamMeanTemperature);
                        upstreamLiquidFlowRate = fabs(state.cells[cellIndex - 1].Mliqini) / ((1 - betneg) * state.cells[cellIndex].flui.MasEspLiq(upstreamMeanPressure, upstreamMeanTemperature) + betneg * state.cells[cellIndex].fluicol.MasEspFlu(upstreamMeanPressure, upstreamMeanTemperature));
                    }

                    estratificado stratifiedMap(diameter, liquidFlowRate, gasFlowRate, liquidDensity, gasDensity, liquidViscosity / pow(10., 3.), gasViscosity / pow(10., 3.), liquidHoldup,
                                            inclinationSign * state.cells[cellIndex].duto.teta, state.cells[cellIndex].duto.rug / diameter);

                    stratifiedMap.mapaTD();
                    flowPattern = stratifiedMap.arr;
                    if (flowPattern == -1) {
                        if (((state.cells[cellIndex].arranjo != flowPattern) || state.cells[cellIndex].transic > 0)) {
                            if ((state.cells[cellIndex].arranjo != flowPattern) && state.cells[cellIndex].transic > 0)
                                state.cells[cellIndex].transic = 0;
                            state.cells[cellIndex].transic++;
                            if (state.cells[cellIndex].transic > 19)
                                state.cells[cellIndex].transic = 0;
                        }
                        state.cells[cellIndex].arranjo = flowPattern = stratifiedMap.arr;
                        if (cellIndex > 0) {
                            state.cells[cellIndex - 1].arranjoR = stratifiedMap.arr;
                            state.cells[cellIndex - 1].perdaEstratL = stratifiedMap.fatorperdaLiq;
                            state.cells[cellIndex - 1].perdaEstratG = stratifiedMap.fatorperdaGas;
                        }

                        RegimePair pair;
                        evaluateRegimePair(state, cellIndex, mix, upstreamLiquidFlowRate, pair);
                        double alf0E = state.cells[cellIndex].alf;
                        if (cellIndex > 0)
                            alf0E = state.cells[cellIndex - 1].alf;

                        blendBySuperficialVelocity(mix, pair, c0, ud);
                    }
                }
                if (flowPattern == 1) {

                    arranjo flowPatternMap(diameter, liquidFlowRate / flowArea, gasFlowRate / flowArea, liquidDensity, gasDensity, liquidViscosity / pow(10., 3.), gasViscosity / pow(10., 3.), liquidHoldup,
                                       inclinationSign * state.cells[cellIndex].duto.teta, surfaceTension, state.input.mapaArranjo, state.globals);
                    flowPattern = flowPatternMap.verificaArr();

                    evaluateDispersedOrAnnular(state, cellIndex, mix, flowPattern, c0, ud);
                    state.cells[cellIndex].arranjo = flowPattern;
                    if (cellIndex > 0)
                        state.cells[cellIndex - 1].arranjoR = flowPattern;
                }
                state.cells[cellIndex].c0Spare = c0;
                state.cells[cellIndex].udSpare = ud;
            }
        } else {
            c0 = 1.;
            ud = 0.;
            state.cells[cellIndex].arranjo = 1;
            if (cellIndex > 0)
                state.cells[cellIndex - 1].arranjoR = 1;
            state.cells[cellIndex].c0Spare = c0;
            state.cells[cellIndex].udSpare = ud;
        }
    }
    applyNoSlipOverride(state.cells, cellIndex, state.input.escorregaPerm, c0, ud);
}

}  // namespace coefficient
}  // namespace driftflux
