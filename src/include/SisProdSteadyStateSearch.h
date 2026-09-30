#ifndef SISPRODSTEADYSTATESEARCH_H_
#define SISPRODSTEADYSTATESEARCH_H_

#include "RootFindingSolvers.h"
#include "SisProdSteadyState.h"

// SisProdSteadyState.h does not include this header, and must not. The march
// does reach two searches -- marchProductionSteady and
// marchProductionSteadySecondary march the gas line after the column
// converges, through searchGasPressureSteadySecondary and
// searchGasPressureSteadyTertiary -- but it reaches them through
// SteadyStateUpdaters, back via SProd, the same way it reaches calctemp. The
// cycle between the two halves is in the data, not in the build, and the
// include stays one-way.

namespace sisprod::steady {

/// The state a boundary-condition search reads.
///
/// It COMPOSES the march state rather than restating it. That is not tidiness:
/// of the 31 SProd members these routines touch, 17 are read by both halves.
/// Restating them would create two spellings of one fact, and the first time
/// they disagreed the compiler would not say so.
///
/// Three members belong to the searches alone, and they are what this struct
/// adds.
struct SteadyStateSearchState {
    /// Everything the march reads. A search passes this straight through when
    /// it calls one.
    SteadyStateState march;

    /// Holdup guess carried between search attempts -- SProd::chuteHol.
    /// Read only.
    const double &holdupGuess;
    /// Reverse-flow network fluid -- SProd::fluiRevRede. Read only.
    const ProFlu &reverseNetworkFluid;
    /// Whether this steady solve runs in reverse -- SProd::revPerm. Written.
    int &reverseSteady;
};

// ------------------------------------------------------ dispatch table ----

/// Runs whichever march the pair (production flag, boundary-condition kind)
/// selects, and returns its residual.
///
/// This is the seam between the two halves: the searches know it by name, the
/// marches do not know it exists. It lives here, with the callers, and not in
/// the march module -- a dispatch table belongs on the side that dispatches.
[[nodiscard]] double dispatchMarch(const SteadyStateSearchState &state, double guess,
                                   int isProduction, int boundaryConditionKind);

/// Drives dispatchMarch to a root between two bracketing guesses.
///
/// The generic half lives in rootfinding::zriddr and knows nothing about
/// production. The domain half stays here: reading acopColAnulPermForte to
/// decide the minimum iteration count, and normalising the residual against
/// the convergence monitor.
[[nodiscard]] double solveSteadyRoot(const SteadyStateSearchState &state, double x1, double x2,
                                     int isProduction, int boundaryConditionKind);

// -------------------------------------------- production bottom-hole search --

/// Searches the bottom-hole pressure that closes the production balance.
/// A guess of -1 means "bracket it yourself"; a negative attempt count means
/// "this is the first attempt".
[[nodiscard]] double searchProductionBottomHolePressure(const SteadyStateSearchState &state,
                                                        double guess = -1., int attemptCount = -1);
[[nodiscard]] double searchReverseProductionBottomHolePressure(const SteadyStateSearchState &state,
                                                               double guess = -1.);
[[nodiscard]] double searchProductionBottomHolePressureSecondary(const SteadyStateSearchState &state,
                                                                 double guess = -1., int attemptCount = -1);
[[nodiscard]] double searchProductionBottomHolePressureTertiary(const SteadyStateSearchState &state,
                                                                double inletPressure);

// ----------------------------------------- production pressure-to-pressure --

/// Searches the mass flow that connects two fixed pressures.
[[nodiscard]] double searchProductionPressureToPressure(const SteadyStateSearchState &state, double massFlowRateGuess,
                                                        double maximumFlowRate = 0., int iterationCount = 0);
[[nodiscard]] double searchReverseProductionPressureToPressure(const SteadyStateSearchState &state,
                                                               double massFlowRateGuess, double maximumFlowRate = 0.,
                                                               int iterationCount = 0);
[[nodiscard]] double searchProductionPressureToPressureSecondary(const SteadyStateSearchState &state,
                                                                 double massFlowRateGuess, double maximumFlowRate = 0.);
[[nodiscard]] double searchProductionPressureToPressureTertiary(const SteadyStateSearchState &state,
                                                                double massFlowRateGuess, double maximumFlowRate = 0.);

// --------------------------------------------------------- gas-line search --

/// Searches the gas-line pressure that closes the injection balance. Not
/// [[nodiscard]]: the march runs them for their effect on the gas line and
/// drops the value.
double searchGasPressureSteadySecondary(const SteadyStateSearchState &state);
double searchGasPressureSteadyTertiary(const SteadyStateSearchState &state);

// -------------------------------------------- injection bottom-hole search --

/// Five variants of the injection bottom-hole search; the well's boundary
/// condition selects one. The first, second and fifth search a root, the third
/// iterates to a fixed point and the fourth is a march.
[[nodiscard]] double searchInjectionBottomHolePressure1(const SteadyStateSearchState &state, double guess = -1.);
[[nodiscard]] double searchInjectionBottomHolePressure2(const SteadyStateSearchState &state, double guess = -1.);
[[nodiscard]] double searchInjectionBottomHolePressure3(const SteadyStateSearchState &state, double guess = -1.);
[[nodiscard]] double searchInjectionBottomHolePressure4(const SteadyStateSearchState &state);
[[nodiscard]] double searchInjectionBottomHolePressure5(const SteadyStateSearchState &state, double guess = -1.);

// --------------------------------------------------- secondary-branch search --

/// Searches the flow rate through the secondary branch.
[[nodiscard]] double searchSecondaryBranchFlowRate(const SteadyStateSearchState &state, double startPressure,
                                                  int startIndex);

}  // namespace sisprod::steady

#endif  // SISPRODSTEADYSTATESEARCH_H_
