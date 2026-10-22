#ifndef SISPRODCOMPOSITION_H_
#define SISPRODCOMPOSITION_H_

// Declared, not included: the state holds only references, so the
// definitions are needed where it is built and not here, and the header
// compiles on its own with nothing but -Isrc/include.
class Cel;
class Ler;
class ProFlu;
struct varGlob1D;
class SProd;
namespace sisprod { class SolveContext; }

namespace sisprod::composition {

/// The one call the transports make outside this module: the gas specific gravity
/// correction, which lives in sisprod::steady and takes the solve context's
/// SteadyStateState.
struct CompositionUpdaters {
    SolveContext &context;

    void correctGasSpecificGravity(int i) const;
};

/// What the composition transport reads and writes.
///
/// There is no time step here although the three transports use one: each
/// declares its own dt, from the first cell's, at its outermost level.
///
/// Scalars are held by reference, not by value: a copy would read each one at
/// construction, before the branch that decides whether it is read at all.
///
/// Nothing is marked const. Field names follow the ones the other modules give
/// the same members.
struct CompositionState {
    Cel *&cells;                          // celula
    int &lastCell;                        // ncel
    Ler &input;                           // arq
    varGlob1D *&globals;                  // vg1dSP
    int &endNode;                         // noextremo
    double &inletPressure;                // presE
    double &inletTemperature;             // tempE
    double &inletComplementaryFraction;      // betaE
    int &trackGasOilRatio;                // trackRGO
    int &trackGasGravity;                 // trackDeng
    int &compositionalRefreshCounter;     // kontaRenovaComp

    CompositionUpdaters updaters;
};

/// Transports black-oil properties -- gas-oil ratio, API gravity, gas density,
/// CO2 fraction -- along the line. fluiRev is the fluid entering through the
/// outlet under reverse flow, taken by value.
void transportBlackOilProperties(const CompositionState &state, ProFlu fluiRev);

/// Transports the molar fractions of the oil and the gas phases separately, and
/// the water density with them. The one the transient step calls.
void transportPhaseMolarFractions(const CompositionState &state, ProFlu fluiRev);

}  // namespace sisprod::composition

#endif  // SISPRODCOMPOSITION_H_
