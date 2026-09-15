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
    const Ler &input;
    /// SProd::celInterIni -- so lido.
    const int &initialInterfaceCell;
    /// SProd::celula -- so lido.
    Cel* cells;
    /// SProd::celulaG -- so lido.
    CelG* gasCells;
    /// SProd::chokeSup -- so lido.
    const choke &surfaceChoke;
    /// SProd::dtCFL -- so lido.
    const std::vector<double> &cflTimeSteps;
    /// SProd::dtInterIni -- so lido.
    const double &initialInterfaceTimeStep;
    /// SProd::dtSim -- so lido.
    const std::vector<double> &simulationTimeSteps;
    /// SProd::dtauxCFL -- so lido.
    const double &auxiliaryCflTimeStep;
    /// SProd::dtauxFinal -- so lido.
    const double &finalAuxiliaryTimeStep;
    /// SProd::fechaM1 -- so lido.
    double* masterCloseSchedule;
    /// SProd::flut -- so lido.
    const FullMtx<double> &productionFreeTerms;
    /// SProd::flutG -- so lido.
    const FullMtx<double> &gasFreeTerms;
    /// SProd::indTramo -- so lido.
    const int &branchIndex;
    /// SProd::jMedMov -- so lido.
    const double &movingMeanFlux;
    /// SProd::kSP -- so lido.
    const int &stepIndex;
    /// SProd::matglobP -- so lido.
    const BandMtx<double> &productionMatrix;
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
    const std::vector<double> &maximumTimeStepRates;
    /// SProd::taxaDpMax -- so lido.
    const std::vector<double> &maximumPressureRates;
    /// SProd::tempE -- so lido.
    const double &inletTemperature;
    /// SProd::titRev -- so lido.
    const double &reverseQuality;
    /// SProd::velInterIni -- so lido.
    const double &initialInterfaceVelocity;
    /// SProd::vg1dSP -- so lido.
    varGlob1D* globals;};

}  // namespace sisprod::transient

#endif  // SISPRODTRANSIENT_H_
