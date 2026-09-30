#ifndef SISPRODCOMPOSITION_H_
#define SISPRODCOMPOSITION_H_

// Declared, not included: the state holds only references and pointers, so the
// definitions are needed where it is built and not here, and the header
// compiles on its own with nothing but -Isrc/include.
class Cel;
class Ler;
class ProFlu;
struct varGlob1D;
class SProd;

namespace sisprod::composition {

/// Calls the transports make back into SProd: one, the gas specific gravity
/// correction (SProd::corrDeng), which delegates to
/// sisprod::steady::correctGasSpecificGravity. Routing it through SProd keeps
/// this module depending on the class rather than on the steady module's
/// state, as the transient module does for its own callbacks.
struct CompositionUpdaters {
    SProd &system;

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
    double &inletQuality;                 // titE
    double &inletCompletionFraction;      // betaE
    int &trackGasOilRatio;                // trackRGO
    int &trackGasGravity;                 // trackDeng
    int &compositionalRefreshCounter;     // kontaRenovaComp
    int &surfaceChokeMassFlag;            // masChkSup
    int &movingPigCount;                  // indpigP
    int &previousMovingPigCount;          // indpigPini
    int &scheduledPigCount;               // npig
    int *&pigReceiverCells;               // receb

    CompositionUpdaters updaters;
};

/// Transports black-oil properties -- gas-oil ratio, API gravity, gas density,
/// CO2 fraction -- along the line. fluiRev is the fluid entering through the
/// outlet under reverse flow, taken by value.
void transportBlackOilProperties(const CompositionState &state, ProFlu fluiRev);

/// Transports the molar fractions of the oil and the gas phases separately, and
/// the water density with them. The one the transient step calls.
void transportPhaseMolarFractions(const CompositionState &state, ProFlu fluiRev);

/// Stores the void and completion fractions of the time level just finished,
/// then moves the pigs and receives those that reach their receiver cell.
void storePreviousFractionsAndMovePigs(const CompositionState &state);

/// Caches cell and face densities so the step does not recompute them.
void cacheCellAndFaceDensities(const CompositionState &state);

/// Evaluates wax deposition in every cell and applies its effect on the pipe
/// wall and its heat transfer.
void evaluateWaxDeposition(const CompositionState &state);

}  // namespace sisprod::composition

#endif  // SISPRODCOMPOSITION_H_
