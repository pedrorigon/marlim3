#ifndef SISPRODTRANSIENT_H_
#define SISPRODTRANSIENT_H_

// Declared, not included. TransientStepState holds only references and
// pointers, so the definitions are needed at the adapter in SisProd.cpp and not
// here -- which is what lets this header compile on its own with nothing but
// -Isrc/include. Including Acidentes2.h for `choke` pulls in Log.h and through
// it rapidjson, and the header stops being self-contained; that is what the
// first draft did and what the acceptance check caught.
//
// The template forward declarations follow SisProdGasLift.h, which established
// the same three for the same reason.
class Cel;
class CelG;
class Ler;
class choke;
struct varGlob1D;
class SProd;
template <class T> class Vcr;
template <class T> class FullMtx;
template <class T> class BandMtx;

#include <vector>

namespace sisprod::transient {

/// The sequence contract, from T015d: the time loop is driven by Num4Main.cpp
/// and is NOT internalised here (FR-031). Nothing in this module advances time
/// on its own; it computes one step and the step size, and the caller decides
/// when to call again.

/// What one transient step reads and writes.
///
/// Seventy-two fields, against twenty-eight for SteadyStateState, and that
/// number is the most useful thing this header records.
///
/// Measured, not estimated: the twenty-two routines this stage moves touch 119
/// of SProd's members. SolveTrans alone touches 85 of them, 48 of which nothing
/// else in the stage touches; the other twenty-one routines share 72 between
/// them, and most touch fewer than twenty.
///
/// That distribution is why the state splits in two rather than being one
/// struct of 119 fields, and it is also a warning. The thermal, gas-lift and
/// steady modules were separable domains that happened to live in one class.
/// The transient step is closer to the class's main loop: it reaches 42% of
/// SProd. Extracting it is still worth doing -- the step and the time-step
/// policy are the numerically load-bearing part of the engine and deserve to be
/// readable -- but nobody should expect the seam to be as clean as stage 5's.
///
/// Scalars are held BY REFERENCE, not by value, for the reason DriftFluxClosure
/// and GasLiftState record: copying them into the struct would read every one at
/// construction, before the branch that decides whether the original would have
/// read it at all.
///
/// const marks what the step does not write. That was measured per field by
/// looking for an assignment, not assumed -- stage 6 shipped a header promising
/// const for a deck the module writes through, stage 7 did the same for a choke,
/// and in both cases the compiler is what said so.
/// The two things the transient step needs from outside itself.
///
/// Only two, and both for stated reasons rather than convenience:
///
///   * geraMiniTabFlu STAYS in SisProd.cpp. PorosoRad-Simples.cpp and
///     solverPoroso.cpp consume it too, so it is surface under FR-038, and
///     T126's acceptance requires atualizaMiniTab to keep invoking it there.
///   * subtempoGas moved to the gas-lift module in stage 6, and reaching it
///     needs a GasLiftState that only SisProd.cpp knows how to assemble -- the
///     same routing SteadyStateUpdaters uses, for the same reason.
struct TransientStepUpdaters {
    SProd &system;

    void generateFluidMiniTable() const;
    void advanceGasSubStep() const;
};

/// const marks what the step does not write, and it marks SCALARS ONLY.
///
/// A class or pointer field can be mutated two ways that do not look like an
/// assignment to the field itself: `input.valTempChokeJus = ...` writes through
/// it, and a non-const method call on `surfaceChoke` mutates it. The generator
/// looked for assignments to the NAME and so promised const for both. The
/// compiler said so -- the fourth and fifth time in this refactoring that a
/// header made that promise, after stage 6's deck and stage 7's choke.
struct TransientStepState {
    /// SProd::DTMaxMed -- escrito.
    double &meanMaximumTimeStep;
    /// SProd::DpMaxMed -- escrito.
    double &meanMaximumPressureChange;
    /// SProd::EstadoMaster1 -- escrito.
    int &masterState;
    /// SProd::aberto -- escrito.
    int &open;
    /// SProd::abertoini -- escrito.
    int &initiallyOpen;
    /// SProd::alteraTempo -- escrito.
    int &timeChanged;
    /// SProd::betaE -- escrito.
    double &inletCompletionFraction;
    /// SProd::celInter -- escrito.
    int &interfaceCell;
    /// SProd::contaMaster1 -- escrito.
    int &masterCounter;
    /// SProd::cpg -- escrito.
    double** gasSpecificHeatTable;
    /// SProd::dt -- escrito.
    double &timeStep;
    /// SProd::dtCFLMed -- escrito.
    double &meanCflTimeStep;
    /// SProd::dtCFLTotal -- escrito.
    double &totalCflTimeStep;
    /// SProd::dtInter -- escrito.
    double &interfaceTimeStep;
    /// SProd::dtSimMed -- escrito.
    double &meanSimulationTimeStep;
    /// SProd::dtSimTotal -- escrito.
    double &totalSimulationTimeStep;
    /// SProd::fontemassCRBuf -- escrito.
    double &bufferedCompletionMassSource;
    /// SProd::fontemassGRBuf -- escrito.
    double &bufferedGasMassSource;
    /// SProd::fontemassPRBuf -- escrito.
    double &bufferedLiquidMassSource;
    /// SProd::indevento -- escrito.
    int &eventIndex;
    /// SProd::kontaGolfada -- escrito.
    int &slugCount;
    /// SProd::kontarestriDt -- escrito.
    int &timeStepRestrictionCount;
    /// SProd::masChkSup -- escrito.
    int &surfaceChokeMassFlag;
    /// SProd::modeloCompleto -- escrito.
    int &fullModel;
    /// SProd::momentoDesesp -- escrito.
    double &desperationMoment;
    /// SProd::mudaModoChk -- escrito.
    int &chokeModeChanged;
    /// SProd::mult -- escrito.
    double &multiplier;
    /// SProd::presfim -- escrito.
    double &finalPressure;
    /// SProd::reinicia -- escrito.
    int &restart;
    /// SProd::restriDt -- escrito.
    int &timeStepRestricted;
    /// SProd::tempoaberto -- escrito.
    int &openTime;
    /// SProd::termolivreP -- escrito.
    Vcr<double> &productionSolution;
    /// SProd::titE -- escrito.
    double &inletQuality;
    /// SProd::vRazMast0 -- escrito.
    double *masterRatio0;
    /// SProd::vRazMast1 -- escrito.
    double *masterRatio1;
    /// SProd::vRazMastCrit -- escrito.
    double *masterCriticalRatio;
    /// SProd::velInter -- escrito.
    double &interfaceVelocity;
    /// SProd::abreM1 -- so lido.
    double* masterOpenSchedule;
    /// SProd::arq -- so lido.
    Ler &input;
    /// SProd::celInterIni -- so lido.
    const int &initialInterfaceCell;
    /// SProd::celula -- so lido.
    Cel* cells;
    /// SProd::celulaG -- so lido.
    CelG* gasCells;
    /// SProd::chokeSup -- so lido.
    choke &surfaceChoke;
    /// SProd::dtCFL -- so lido.
    std::vector<double> &cflTimeSteps;
    /// SProd::dtInterIni -- so lido.
    const double &initialInterfaceTimeStep;
    /// SProd::dtSim -- so lido.
    std::vector<double> &simulationTimeSteps;
    /// SProd::dtauxCFL -- so lido.
    const double &auxiliaryCflTimeStep;
    /// SProd::dtauxFinal -- so lido.
    const double &finalAuxiliaryTimeStep;
    /// SProd::fechaM1 -- so lido.
    double* masterCloseSchedule;
    /// SProd::flut -- so lido.
    FullMtx<double> &productionFreeTerms;
    /// SProd::flutG -- so lido.
    FullMtx<double> &gasFreeTerms;
    /// SProd::indTramo -- so lido.
    const int &branchIndex;
    /// SProd::jMedMov -- so lido.
    const double &movingMeanFlux;
    /// SProd::kSP -- so lido.
    const int &stepIndex;
    /// SProd::matglobP -- so lido.
    BandMtx<double> &productionMatrix;
    /// SProd::menorDx -- so lido.
    const double &smallestCellLength;
    /// SProd::nabreM1 -- so lido.
    const int &masterOpenCount;
    /// SProd::ncel -- so lido.
    const int &lastCell;
    /// SProd::ncelGas -- so lido.
    const int &gasCellCount;
    /// SProd::ncelperftransp -- so lido.
    int* productionCrossSectionCount;
    /// SProd::nfechaM1 -- so lido.
    const int &masterCloseCount;
    /// SProd::noextremo -- so lido.
    const int &endNode;
    /// SProd::pGSup -- so lido.
    const double &gasSurfacePressure;
    /// SProd::presE -- so lido.
    const double &inletPressure;
    /// SProd::presMedMov -- so lido.
    const double &movingMeanPressure;
    /// SProd::tMedMov -- so lido.
    const double &movingMeanTemperature;
    /// SProd::taxaDTMax -- so lido.
    std::vector<double> &maximumTimeStepRates;
    /// SProd::taxaDpMax -- so lido.
    std::vector<double> &maximumPressureRates;
    /// SProd::tempE -- so lido.
    const double &inletTemperature;
    /// SProd::titRev -- so lido.
    const double &reverseQuality;
    /// SProd::velInterIni -- so lido.
    const double &initialInterfaceVelocity;
    /// SProd::vg1dSP -- so lido.
    varGlob1D* globals;
    /// Everything the step needs that is not its own.
    TransientStepUpdaters updaters;
};

// ------------------------------------------------------------ cell update ----

/// Updates every cell after a transient solve, in index order.
///
/// The iteration order over cells is load-bearing and is why the three arms
/// below are separate functions rather than one parameterised by position: the
/// split is the loop body's own structure, so the order is untouched.
void updateCells(const TransientStepState &state, int expli);
void updateInteriorCell(const TransientStepState &state, int i, int expli);
void updateFirstCell(const TransientStepState &state, int i, int expli);
void updateLastCell(const TransientStepState &state, int i, int expli);

/// Updates the flow rates. Writes a strict subset of what updateCells writes --
/// measured in evidencia/renova-diff.md, where it is also recorded that the
/// textual similarity metric reports 37% for a pair that shares every field.
void updateFlowRates(const TransientStepState &state);

// ---------------------------------------------------------- buffer update ----

/// Fills the buffered state from the solver's free-term vector, and from the
/// cells' own current values.
///
/// TWO functions, and they stay two. Num4Main.cpp picks between them in
/// alternative branches of one if, and they differ by more than their source:
/// updateBufferFromSolution propagates state to the last cell's right face and
/// does not touch the mass sources, updateBufferFromCells does exactly the
/// opposite. Nothing in the code says whether that asymmetry is intended.
void updateBufferFromSolution(const TransientStepState &state);
void updateBufferFromCells(const TransientStepState &state);

// ------------------------------------------------- outlet boundary condition ----

/// The surface choke is open: throat area above a thousandth of the pipe's.
/// And shut: throat area BELOW that.
///
/// NOT each other's negation, and nothing should write them as if they were.
/// `!(a > b)` is `a <= b`; the second of these is `a < b`. They differ at exact
/// equality.
[[nodiscard]] bool surfaceChokeIsOpen(const TransientStepState &state);
[[nodiscard]] bool surfaceChokeIsShut(const TransientStepState &state);

/// Applies the outlet boundary condition, against the real state and against
/// the buffered state.
///
/// Two functions, 52% alike by structure and 37% by text; the gap is measured
/// in evidencia/calccc-diff.md. They differ by more than their source: the
/// pressure form computes the Joule-Thomson temperature downstream of the choke
/// and the buffer form does not; the buffer form propagates the mass sources
/// and the pressure form does not.
///
/// SProdVap has its OWN calcCCpres and calcCCBuffer, two arguments instead of
/// three -- parallel implementations in a different class, not overloads, and
/// out of scope under FR-038. Anyone reading this module for "the" boundary
/// condition is reading half of it.
void applyOutletPressureCondition(const TransientStepState &state, double titRev, double alfRev, double betRev);
void applyOutletBufferCondition(const TransientStepState &state, double titRev, double alfRev, double betRev);

// ------------------------------------------------------------- time step ----

/// Decides the time step, explicitly or implicitly.
///
/// This is the numerically load-bearing function of the whole refactoring. A
/// last-bit drift here does not make a small difference in the answer: it makes
/// a DIFFERENT TEMPORAL DISCRETISATION, and from that step onward the run is a
/// different simulation.
///
/// So it is not verified by L2 alone. T122 captured the complete series of
/// 210,206 calls in %a, and refactor-harness/verify-determinadt.sh compares the
/// current tree against it call by call. L2 would say "the outputs differ"; the
/// series says where the mesh first diverged, which is the only actionable
/// answer once the step feeds back into everything.
///
/// Coverage, stated because it limits what any of this proves: EIGHT of the
/// fourteen corpus models reach this function at all. A defect here is
/// invisible to the other six, with their L2 still green.
void computeTimeStep(const TransientStepState &state, int vexpli);
void computeExplicitTimeStep(const TransientStepState &state);
void computeImplicitTimeStep(const TransientStepState &state);

}  // namespace sisprod::transient

#endif  // SISPRODTRANSIENT_H_
