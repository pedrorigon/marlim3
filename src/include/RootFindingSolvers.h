#ifndef ROOTFINDINGSOLVERS_H_
#define ROOTFINDINGSOLVERS_H_

/// Generic root-finding algorithms.
///
/// These three solvers were members of SProd and took `prod` and `tipoCC`,
/// two production-domain selectors, purely to hand them back to
/// SProd::multMarcha on every evaluation. Neither means anything to Brent's
/// method. They are gone from the signatures: the caller builds a callable
/// that has already captured them, and the solver sees only `objective(x)`.
///
/// Everything here lives in namespace rootfinding, and that is load-bearing
/// rather than tidy. FerramentasNumericas.h already declares templates named
/// zbrent and zriddr at global scope, and SisProd.h includes it. Free
/// functions of those names would be ambiguous at best; at worst a call would
/// resolve silently to the other overload, which takes `const T *const` and
/// would compile.

// fabs and sqrt. This is the include SisProd.cpp used when the solvers lived
// there, so overload resolution is the one they were written against.
#include <math.h>

namespace rootfinding {

/// Magnitude of the first argument carrying the sign of the second.
///
/// inline, and defined here rather than in the .cpp, because zriddr calls it
/// three times per iteration. Out of line it would become a cross-TU call in
/// the one solver that is actually hot.
inline double SIGN(double magnitude, double signSource) {
    return (signSource >= 0 ? 1.0 : -1.0) * fabs(magnitude);
}

/// Sign of a value, as -1 or 1, with zero counting as negative.
///
/// It has no caller anywhere in the project and never had one; it is preserved
/// because removing dead code is a behaviour change this programme is not
/// authorised to make. Out of line precisely because nothing calls it.
int sign(double value);

/// Reports that a solver exhausted its iteration budget.
///
/// A one-line forward to NumError, which is declared in FerramentasNumericas.h
/// together with `using namespace std;` at file scope and the colliding
/// zbrent/zriddr templates. Reaching it through the .cpp keeps all of that out
/// of every translation unit that includes this header.
void reportIterationLimit(const char *message);

/// Finds a root by bisection -- despite the name, which says false position.
///
/// The name is Portuguese for "false chord", the regula falsi. The body is not
/// that: it takes the midpoint of the bracket, not the intercept of the secant,
/// and the rest of the loop is a textbook bisection update. The loop carries a
/// comment from the original -- "this block treats the 'falsacorda' properly",
/// with the name in quotes -- and multFC, the unused 0.5 below, is the very
/// multiplier hardcoded into the line that halves the interval. Whoever wrote
/// it knew. See A2-06 in evidencia/anomalias.md.
///
/// The name is kept because renaming a function is not in this stage's scope
/// and the stage contract fixes it. If it is ever renamed, the right name is
/// bisect, not falsePosition -- that would make the name lie with more
/// authority.
///
/// Reachable only from zbrent, which nothing calls, so it never executes. Its
/// verification is structural: refactor-harness/solver-move.py inverts the move
/// and compares the token stream against the pristine baseline.
template <typename Objective>
double falsacorda(double bracketLow, double bracketHigh, Objective &&objective) {
    double lowValue = objective(bracketLow);
    double halfWidth = bracketHigh - bracketLow;
    double midpoint;
    int maximumIterations = 100;
    double intervalTolerance = 0.001;
    double valueTolerance = 0.001;
    double multFC = 0.5;

    for (int iteration = 1; iteration <= maximumIterations; iteration++) { // this block treats the 'falsacorda' properly
        halfWidth = bracketHigh - bracketLow;
        halfWidth *= 0.5;
        midpoint = bracketLow + halfWidth;
        double midpointValue = objective(midpoint);
        if (fabs(halfWidth) < intervalTolerance || fabs(midpointValue) < valueTolerance)
            return midpoint;
        ((lowValue > 0 && midpointValue < 0) || (lowValue < 0 && midpointValue > 0)) ? (bracketHigh = midpoint) : (bracketLow = midpoint, lowValue = midpointValue);
    }
    return midpoint;
}

/// Finds a root by Brent's method: bracketing with inverse quadratic
/// interpolation, falling back to false position when the interval does not
/// bracket a sign change.
///
/// Measured over the whole tree, this has NO call site. It is moved as it
/// stands, and never runs. Two consequences worth stating where they will be
/// read: L2 and L3 cannot see a defect introduced here, and because an
/// uninstantiated template is only parsed, neither can the compiler. What
/// covers it is solver-move.py for the move and verify-solvers.sh, which
/// instantiates and exercises it, for everything after.
///
/// The declaration this replaced carried default arguments -- tol and epsn both
/// 0.00001, maxit 100. They are not reproduced, because a default on a function
/// with no caller only invites one to be written without thinking about the
/// tolerance; the values are recorded here instead, since they are the only
/// statement anyone ever made about what this solver expects.
template <typename Objective>
double zbrent(double bracketLow, double bracketHigh, Objective &&objective, double absoluteTolerance, double relativeTolerance, int maximumIterations) {
    double relativePrecision = relativeTolerance;
    double previousEstimate = bracketLow;
    double currentEstimate = bracketHigh;
    double oppositeSignPoint = bracketHigh;
    double previousValue = objective(previousEstimate);
    double currentValue = objective(currentEstimate);
    if (fabs(previousValue) > 1e9 || fabs(currentValue) > 1e9)
        return 1e10;
    double previousStep = 0.;
    double step, oppositeSignValue, stepNumerator, stepDenominator, valueRatioOpposite, valueRatioPrevious, workingTolerance, halfBracketWidth;

    if ((previousValue > 0.0 && currentValue > 0.0) || (previousValue < 0.0 && currentValue < 0.0)) {
        double fallbackRoot;
        fallbackRoot = falsacorda(bracketLow, bracketHigh, objective);
        return fallbackRoot;
    } else {
        oppositeSignValue = currentValue;
        for (int iteration = 0; iteration < maximumIterations; iteration++) {
            if ((currentValue > 0.0 && oppositeSignValue > 0.0) || (currentValue < 0.0 && oppositeSignValue < 0.0)) {
                oppositeSignPoint = previousEstimate;
                oppositeSignValue = previousValue;
                previousStep = step = currentEstimate - previousEstimate;
            }
            if (fabs(oppositeSignValue) < fabs(currentValue)) {
                previousEstimate = currentEstimate;
                currentEstimate = oppositeSignPoint;
                oppositeSignPoint = previousEstimate;
                previousValue = currentValue;
                currentValue = oppositeSignValue;
                oppositeSignValue = previousValue;
            }
            workingTolerance = 2.0 * relativePrecision * fabs(currentEstimate) + 0.5 * absoluteTolerance;
            halfBracketWidth = 0.5 * (oppositeSignPoint - currentEstimate);
            if (fabs(halfBracketWidth) <= workingTolerance || currentValue == 0.0)
                return currentEstimate;
            if (fabs(previousStep) >= workingTolerance && fabs(previousValue) > fabs(currentValue)) {
                valueRatioPrevious = currentValue / previousValue;
                if (previousEstimate == oppositeSignPoint) {
                    stepNumerator = 2.0 * halfBracketWidth * valueRatioPrevious;
                    stepDenominator = 1.0 - valueRatioPrevious;
                } else {
                    stepDenominator = previousValue / oppositeSignValue;
                    valueRatioOpposite = currentValue / oppositeSignValue;
                    stepNumerator = valueRatioPrevious * (2.0 * halfBracketWidth * stepDenominator * (stepDenominator - valueRatioOpposite) - (currentEstimate - previousEstimate) * (valueRatioOpposite - 1.0));
                    stepDenominator = (stepDenominator - 1.0) * (valueRatioOpposite - 1.0) * (valueRatioPrevious - 1.0);
                }
                if (stepNumerator > 0.0)
                    stepDenominator = -stepDenominator;
                stepNumerator = fabs(stepNumerator);
                double interpolationLimit = 3.0 * halfBracketWidth * stepDenominator - fabs(workingTolerance * stepDenominator);
                double previousStepLimit = fabs(previousStep * stepDenominator);
                if (2.0 * stepNumerator < (interpolationLimit < previousStepLimit ? interpolationLimit : previousStepLimit)) {
                    previousStep = step;
                    step = stepNumerator / stepDenominator;
                } else {
                    step = halfBracketWidth;
                    previousStep = step;
                }
            } else {
                step = halfBracketWidth;
                previousStep = step;
            }
            previousEstimate = currentEstimate;
            previousValue = currentValue;
            if (fabs(step) > workingTolerance)
                currentEstimate += step;
            else
                currentEstimate += SIGN(workingTolerance, halfBracketWidth);
            currentValue = objective(currentEstimate);
        }
        reportIterationLimit("Metodo Van Winjngaarden-Dekker-Brent para calcular zero de funcaoo atingiu maximo de iteracoes");
        return 0.0;
    }
}

/// Finds a root by Ridders' method. The only solver here that executes.
///
/// Two callables, not one, because the original evaluates the objective in two
/// different ways. Twelve of its fourteen evaluations are raw; two are divided
/// by a convergence-monitor base and recorded in a member the outer pressure
/// loops read back. That scaling and that write are domain feedback, so they
/// travel in `monitor` and the solver keeps only the composition
/// `monitor(objective(x))` -- the same order of operations the original had.
///
/// `reverseMarch` selects nothing today: all three branches that test it have
/// identical arms. They are kept verbatim rather than collapsed. Collapsing
/// would be behaviour-preserving -- reading an int member has no side effect --
/// but the branch is the only surviving evidence that someone meant to treat
/// the reverse march differently here, and erasing it would make that
/// unrecoverable from the code. See A2-02 to A2-04 in evidencia/anomalias.md.
///
/// `minimumIterations` is derived from the input deck at the binding site. It
/// also gates three early returns that would otherwise be unconditional, which
/// is how the division at A2-05 becomes reachable.
template <typename Objective, typename Monitor>
double zriddr(double bracketLow, double bracketHigh, Objective &&objective, Monitor &&monitor,
              int reverseMarch, int minimumIterations) {
    double rootAccuracy = 1e-5;
    int maximumIterations = 100;
    double bestValue;
    double bestPoint;
    double lowValue;
    double highValue;
    if (reverseMarch == 0) {
        lowValue = objective(bracketLow);
        highValue = objective(bracketHigh);
    } else {
        lowValue = objective(bracketLow);
        highValue = objective(bracketHigh);
    }
    if (fabs(lowValue) > 1e9 || fabs(highValue) > 1e9)
        return 1e10;
    if (lowValue >= 0.) {
        if (lowValue > 0.9e10) {
            bracketLow *= 0.9999;
            lowValue = objective(bracketLow);
        } else {
            int attempts = 0;
            while (attempts < 100 && lowValue > 0.) {
                if (reverseMarch == 0) {
                    if (bracketHigh < bracketLow)
                        bracketLow *= 1.0001;
                    else
                        bracketLow *= 0.999;
                } else {
                    if (bracketHigh < bracketLow)
                        bracketLow *= 1.0001;
                    else
                        bracketLow *= 0.999;
                }
                lowValue = objective(bracketLow);
                attempts++;
            }
        }
    } else if (highValue <= 0.) {
        if (highValue < -0.9e10) {
            bracketHigh *= 1.00001;
            highValue = objective(bracketHigh);
        } else {
            int attempts = 0;
            while (attempts < 100 && highValue < 0.) {
                if (reverseMarch == 0) {
                    if (bracketLow < bracketHigh)
                        bracketHigh *= 1.0001;
                    else
                        bracketHigh *= 0.999;
                } else {
                    if (bracketLow < bracketHigh)
                        bracketHigh *= 1.0001;
                    else
                        bracketHigh *= 0.999;
                }
                highValue = objective(bracketHigh);
                attempts++;
            }
        }
    }
    if (fabs(highValue) < fabs(lowValue)) {
        bestValue = highValue;
        bestPoint = bracketHigh;
    } else {
        bestValue = lowValue;
        bestPoint = bracketLow;
    }
    if ((lowValue > 0.0 && highValue < 0.0) || (lowValue < 0.0 && highValue > 0.0)) {
        double intervalLow = bracketLow;
        double intervalHigh = bracketHigh;
        double previousAnswer = -1.e20;
        for (int iteration = 0; iteration < maximumIterations; iteration++) {
            double midpoint = 0.5 * (intervalLow + intervalHigh);
            double midpointValue = monitor(objective(midpoint));
            if (fabs(midpointValue) < fabs(bestValue)) {
                bestValue = midpointValue;
                bestPoint = midpoint;
            }
            double discriminant = sqrt(midpointValue * midpointValue - lowValue * highValue);
            if (discriminant == 0.0) {
                bestValue = objective(bestPoint);
                if(iteration>minimumIterations)return bestPoint;
            }
            double nextPoint = midpoint + (midpoint - intervalLow) * ((lowValue >= highValue ? 1.0 : -1.0) * midpointValue / discriminant);
            if (fabs(nextPoint - previousAnswer) <= rootAccuracy) {
                bestValue = objective(bestPoint);
                if(iteration>minimumIterations)return bestPoint;
            }
            previousAnswer = nextPoint;
            double nextValue = monitor(objective(previousAnswer));
            if (fabs(nextValue) < fabs(bestValue)) {
                bestValue = nextValue;
                bestPoint = previousAnswer;
            }
            if (fabs(nextValue) <= rootAccuracy) {
                bestValue = objective(bestPoint);
                if(iteration>minimumIterations)return bestPoint;
            }
            if (SIGN(midpointValue, nextValue) != midpointValue) {
                intervalLow = midpoint;
                lowValue = midpointValue;
                intervalHigh = previousAnswer;
                highValue = nextValue;
            } else if (SIGN(lowValue, nextValue) != lowValue) {
                intervalHigh = previousAnswer;
                highValue = nextValue;
            } else if (SIGN(highValue, nextValue) != highValue) {
                intervalLow = previousAnswer;
                lowValue = nextValue;
            } else
                return -1.e10;
            if (fabs(intervalHigh - intervalLow) <= rootAccuracy) {
                bestValue = objective(bestPoint);
                return bestPoint;
            }
        }
        return 1.e10;
    } else {
        if (fabs(lowValue) <= rootAccuracy) {
            return bracketLow;
        }
        if (fabs(highValue) <= rootAccuracy) {
            return bracketHigh;
        }
        return -1e10;
    }
}

}  // namespace rootfinding

#endif  // ROOTFINDINGSOLVERS_H_
