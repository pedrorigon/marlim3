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

// For the two concepts below. Both are lightweight and neither drags in
// anything that could change overload resolution for fabs or sqrt.
#include <concepts>
#include <type_traits>

namespace rootfinding {

/// A residual function a solver drives to zero: position in, residual out.
///
/// Constrained rather than left as a bare `typename`, and the reason is
/// specific to this module. zbrent and bisect have no call site in the product,
/// so the only thing that ever instantiates them is the verification harness. A
/// caller that passes the wrong shape would otherwise get a page of diagnostics
/// from inside the arithmetic, naming variables it never heard of, instead of
/// one line saying the argument does not satisfy ObjectiveFunction.
template <typename Function>
concept ObjectiveFunction =
    std::invocable<Function, double> &&
    std::convertible_to<std::invoke_result_t<Function, double>, double>;

/// Post-processing a solver applies to a residual before using it.
///
/// Same shape as ObjectiveFunction and a different role, which is why it has
/// its own name: zriddr composes them as monitor(objective(x)), and reading
/// that signature should not require working out which double means what.
template <typename Function>
concept ResidualMonitor = ObjectiveFunction<Function>;

/// Magnitude of the first argument carrying the sign of the second.
///
/// inline, and defined here rather than in the .cpp, because zriddr calls it
/// three times per iteration. Out of line it would become a cross-TU call in
/// the one solver that is actually hot.
/// \param magnitude   Value whose absolute value is taken.
/// \param signSource  Value whose sign is applied; zero counts as positive.
/// \return |magnitude|, signed like signSource.
[[nodiscard]] inline double SIGN(double magnitude, double signSource) {
    return (signSource >= 0 ? 1.0 : -1.0) * fabs(magnitude);
}

/// Sign of a value, as -1 or 1, with zero counting as negative.
///
/// It has no caller anywhere in the project and never had one; it is preserved
/// because removing dead code is a behaviour change this programme is not
/// authorised to make. Out of line precisely because nothing calls it.
///
/// \param value  Value to inspect.
/// \return -1 when value <= 0, otherwise 1.
[[nodiscard]] int sign(double value);

/// Reports that a solver exhausted its iteration budget.
///
/// A one-line forward to NumError, which is declared in FerramentasNumericas.h
/// together with `using namespace std;` at file scope and the colliding
/// zbrent/zriddr templates. Reaching it through the .cpp keeps all of that out
/// of every translation unit that includes this header.
///
/// \param message  Text handed to NumError verbatim.
void reportIterationLimit(const char *message);

/// Finds a root by bisection: halve the bracket, keep the half that still
/// straddles the sign change.
///
/// Was SProd::falsacorda -- Portuguese for "false chord", the regula falsi --
/// and the body was never that. It takes the midpoint of the bracket, not the
/// intercept of the secant. Whoever wrote it knew: the loop carries the original
/// comment "this block treats the 'falsacorda' properly", with the name in
/// quotes, and multFC below is the unused 0.5 that the halving line hardcodes.
/// Renamed rather than translated, because falsePosition would have made the
/// name lie with more authority. See A2-06 in evidencia/anomalias.md.
///
/// Reachable only from zbrent, which nothing calls, so it never executes. Its
/// verification is refactor-harness/solver-move.py for the move and
/// verify-solvers.sh, which instantiates and exercises it, for everything since.
/// \tparam Objective  Residual function; see the ObjectiveFunction concept.
/// \param bracketLow   Interval endpoint the sign test treats as the low side.
/// \param bracketHigh  Interval endpoint the sign test treats as the high side.
/// \param objective    Evaluated once before the loop and once per iteration.
/// \return The midpoint reached when the interval or the residual falls below
///         its tolerance, or the last midpoint if the budget runs out.
template <ObjectiveFunction Objective>
[[nodiscard]] double bisect(double bracketLow, double bracketHigh, Objective &&objective) {
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
/// interpolation, falling back to bisection when the interval does not bracket
/// a sign change.
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
/// \tparam Objective          Residual function; see the ObjectiveFunction concept.
/// \param bracketLow          Interval endpoint.
/// \param bracketHigh         Interval endpoint.
/// \param objective           Evaluated twice up front, then once per iteration.
/// \param absoluteTolerance   Absolute half-width the bracket must reach.
/// \param relativeTolerance   Relative precision scaling the working tolerance.
/// \param maximumIterations   Iteration budget before reportIterationLimit fires.
/// \return The root, or 1e10 when either endpoint evaluates beyond 1e9, or 0.0
///         when the budget is exhausted.
template <ObjectiveFunction Objective>
[[nodiscard]] double zbrent(double bracketLow, double bracketHigh, Objective &&objective, double absoluteTolerance, double relativeTolerance, int maximumIterations) {
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
        fallbackRoot = bisect(bracketLow, bracketHigh, objective);
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
/// is how the division at A2-05 becomes reachable. Deriving it before the solve
/// rather than inside it is safe because nothing in SisProd.cpp assigns the flag
/// it comes from -- all seven assignments are in Leitura.cpp, parsing decks.
///
/// `reverseMarch` is passed once, on entry, and the original read it three times
/// during the iteration. That is only safe because those three branches have
/// identical arms today. If A2-02 is ever corrected so that they differ, this
/// has to go back to being read inside the loop -- the value can change mid
/// solve, since SisProd.cpp assigns revPerm in eighteen places.
/// \tparam Objective          Residual function; see the ObjectiveFunction concept.
/// \tparam Monitor            Residual post-processing; see ResidualMonitor.
/// \param bracketLow          Endpoint where the residual is expected negative.
/// \param bracketHigh         Endpoint where the residual is expected positive.
/// \param objective           Evaluated twelve times per solve on the raw path.
/// \param monitor             Applied to the two evaluations that feed the
///                            convergence monitor, as monitor(objective(x)).
/// \param reverseMarch        Selects nothing today; see A2-02 above.
/// \param minimumIterations   Iterations that must pass before the three early
///                            returns are honoured.
/// \return The best position found, or 1e10 / -1e10 / 1.e10 sentinels for the
///         out-of-range, no-bracket and budget-exhausted cases respectively.
template <ObjectiveFunction Objective, ResidualMonitor Monitor>
[[nodiscard]] double zriddr(double bracketLow, double bracketHigh, Objective &&objective, Monitor &&monitor,
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
