#ifndef SISPRODCOMPOSITION_H_
#define SISPRODCOMPOSITION_H_

// Declared, not included, for the reason SisProdTransient.h records: the state
// holds only references and pointers, so the definitions are needed at the
// adapter in SisProd.cpp and not here, and the header compiles on its own with
// nothing but -Isrc/include.
class Cel;
class Ler;
class ProFlu;
struct varGlob1D;
class SProd;

namespace sisprod::composition {

/// Calls the moved bodies make back into SProd.
///
/// One, measured: the black-oil and molar-fraction transports correct each
/// cell's gas specific gravity through SProd::corrDeng, which is itself a
/// delegation to sisprod::steady::correctGasSpecificGravity. Routing it through
/// SProd keeps this module depending on the class rather than on the steady
/// module's state, as the transient module does for its own callbacks.
struct CompositionUpdaters {
    SProd &system;

    void correctGasSpecificGravity(int i) const;
};

/// What the composition transport reads and writes.
///
/// Seventeen fields, derived from measure-members.py over the six routines this
/// module takes: renovaRGOdgYco2, renovaFracMol, renovaFracMol2, renovaalbetini,
/// renovaMasEsp and avaliaParafina. Against 72 for TransientStepState and 28
/// for SteadyStateState, this is the narrowest state any stage has needed --
/// composition transport is a separable domain that happened to live in the
/// class, as the thermal and gas-lift ones were.
///
/// There is no time step here although the three transports use one: each
/// declares `double dt = celula[1].dt;` at its outermost level, which hides
/// SProd::dt for the whole body. The first measurement counted those as uses of
/// the member and put a timeStep field in; the compiler's -Wshadow and the
/// corrected measurement agree it was never read.
///
/// Scalars are held BY REFERENCE, not by value, for the reason every state in
/// this refactoring records: copying them in would read each one at
/// construction, before the branch that decides whether the original would
/// have read it at all.
///
/// Nothing is marked const yet. Four headers before this one promised const
/// from a reading of the code and the compiler refused the promise each time;
/// here the marking is left to the move, where the compiler decides it per
/// field.
///
/// Field names follow the ones the other modules already gave the same members
/// (SC-017); pig names are new and come from the members' own documentation.
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
/// outlet under reverse flow; by value, as the original takes it (FR-035 weighs
/// const& separately, at T110a).
void transportBlackOilProperties(const CompositionState &state, ProFlu fluiRev);

/// Transports the OVERALL compositional molar fractions. No caller anywhere, in
/// this tree or in main; kept as public surface until its owner decides, and
/// compared line by line with the one that runs in evidencia/fracmol-diff.md --
/// 71 of their 73 differences are logic, so they stay two functions.
void transportOverallMolarFractions(const CompositionState &state, ProFlu fluiRev);

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
