#ifndef SISPRODSTEADYSTATESEARCH_H_
#define SISPRODSTEADYSTATESEARCH_H_

#include "RootFindingSolvers.h"
#include "SisProdSteadyState.h"

// SisProdSteadyState.h deliberately does NOT include this header, and that
// stays true. The reason first written here did not.
//
// It claimed the dependency runs one way, that no march calls a search. T090
// measured otherwise: marchaProdPerm1 and marchaProdPerm2 march the gas line
// after the column converges, and that block calls buscaGasPresPerm2 and
// buscaGasPresPerm3, both of which live on this side. The call graph between
// the two halves has a cycle.
//
// What holds the cut together is therefore weaker than an absent edge, and
// worth stating as what it is: the march reaches those two through
// SteadyStateUpdaters, back via SProd, the same way it reaches calctemp. The
// cycle is an edge in the data, not an edge in the build. So the include stays
// one-way and T098b's grep still means something -- but it means "no
// compile-time cycle", not "no dependency".
//
// If this include ever needs to be reversed, the cut is wrong, not the include.
// See section 8 of evidencia/marchaprod-diff.md.

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
[[nodiscard]] double searchProductionPressureToPressure(const SteadyStateSearchState &state, double massFlowGuess,
                                                        double maximumFlowRate = 0., int iterationCount = 0);
[[nodiscard]] double searchReverseProductionPressureToPressure(const SteadyStateSearchState &state,
                                                               double massFlowGuess, double maximumFlowRate = 0.,
                                                               int iterationCount = 0);
[[nodiscard]] double searchProductionPressureToPressureSecondary(const SteadyStateSearchState &state,
                                                                 double massFlowGuess, double maximumFlowRate = 0.);
[[nodiscard]] double searchProductionPressureToPressureTertiary(const SteadyStateSearchState &state,
                                                                double massFlowGuess, double maximumFlowRate = 0.);

// --------------------------------------------------------- gas-line search --

/// Searches the gas-line pressure that closes the injection balance.
[[nodiscard]] double searchGasPressureSteadySecondary(const SteadyStateSearchState &state);
[[nodiscard]] double searchGasPressureSteadyTertiary(const SteadyStateSearchState &state);

// -------------------------------------------- injection bottom-hole search --

/// Five variants of the injection bottom-hole search. They are numbered rather
/// than named because the original numbered them and nothing in the code says
/// what distinguishes four from five; T095 compares them before any of them is
/// unified.
[[nodiscard]] double searchInjectionBottomHolePressure1(const SteadyStateSearchState &state, double guess = -1.);
[[nodiscard]] double searchInjectionBottomHolePressure2(const SteadyStateSearchState &state, double guess = -1.);
[[nodiscard]] double searchInjectionBottomHolePressure3(const SteadyStateSearchState &state, double guess = -1.);
[[nodiscard]] double searchInjectionBottomHolePressure4(const SteadyStateSearchState &state);
[[nodiscard]] double searchInjectionBottomHolePressure5(const SteadyStateSearchState &state, double guess = -1.);

}  // namespace sisprod::steady

#endif  // SISPRODSTEADYSTATESEARCH_H_
