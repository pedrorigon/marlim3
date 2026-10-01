#ifndef SISPRODTRANSIENT_H_
#define SISPRODTRANSIENT_H_

// Declared, not included. TransientStepState holds only references and
// pointers, so the definitions are needed where the state is built and not
// here, and this header compiles on its own with nothing but -Isrc/include.
// Including Acidentes2.h for `choke` would pull in Log.h and through it
// rapidjson, and the header would stop being self-contained.
//
// The template forward declarations are the same three SisProdGasLift.h makes,
// for the same reason.
class Cel;
class CelG;
class Ler;
class choke;
class ChokeGas;
class ProFlu;
class solverP3D;
struct varGlob1D;
class SProd;
template <class T> class Vcr;
template <class T> class FullMtx;
template <class T> class BandMtx;

#include <chrono>
#include <string>
#include <vector>

namespace sisprod::transient {

/// The time loop is driven by Num4Main.cpp, not by this module: it computes one
/// step and the step size, and the caller decides when to call again.

/// The one thing the transient step needs from outside itself: subtempoGas,
/// which lives in the gas-lift module and needs a GasLiftState that only the
/// adapters assemble -- the same routing SteadyStateUpdaters uses.
struct TransientStepUpdaters {
    SProd &system;

    void advanceGasSubStep() const;
};

/// What one transient step reads and writes.
///
/// Scalars are held by reference, not by value: a copy would read every one at
/// construction, before the branch that decides whether it is read at all.
///
/// const marks what the step does not write, and it marks scalars only. A class
/// or pointer field can be changed without an assignment to the field itself:
/// `input.valTempChokeJus = ...` writes through it, and a non-const method call
/// on `surfaceChoke` changes it. The fields SolveTrans writes, directly or by
/// reference through Ler::atualiza (pGSup, presE, tempE, titE, betaE), are not
/// const either.
struct TransientStepState {
    /// SProd::transient.DTMaxMed -- written.
    double &meanMaximumTimeStep;
    /// SProd::transient.DpMaxMed -- written.
    double &meanMaximumPressureChange;
    /// SProd::transient.EstadoMaster1 -- written.
    int &masterState;
    /// SProd::transient.aberto -- written.
    int &surfaceChokeOpen;
    /// SProd::transient.abertoini -- written.
    int &initiallyOpen;
    /// SProd::transient.alteraTempo -- written.
    int &timeChanged;
    /// SProd::betaE -- written.
    double &inletCompletionFraction;
    /// SProd::gasLift.celInter -- written.
    int &interfaceCell;
    /// SProd::transient.contaMaster1 -- written.
    int &masterCounter;
    /// SProd::tables.cpg -- written.
    double** gasSpecificHeatTable;
    /// SProd::dt -- written.
    double &timeStep;
    /// SProd::transient.dtCFLMed -- written.
    double &meanCflTimeStep;
    /// SProd::transient.dtCFLTotal -- written.
    double &totalCflTimeStep;
    /// SProd::gasLift.dtInter -- written.
    double &interfaceTimeStep;
    /// SProd::transient.dtSimMed -- written.
    double &meanSimulationTimeStep;
    /// SProd::transient.dtSimTotal -- written.
    double &totalSimulationTimeStep;
    /// SProd::transient.fontemassCRBuf -- written.
    double &bufferedCompletionMassSource;
    /// SProd::transient.fontemassGRBuf -- written.
    double &bufferedGasMassSource;
    /// SProd::transient.fontemassPRBuf -- written.
    double &bufferedLiquidMassSource;
    /// SProd::transient.indevento -- written.
    int &eventIndex;
    /// SProd::transient.kontaGolfada -- written.
    int &slugCount;
    /// SProd::transient.kontarestriDt -- written.
    int &timeStepRestrictionCount;
    /// SProd::masChkSup -- written.
    int &surfaceChokeMassFlag;
    /// SProd::transient.modeloCompleto -- written.
    int &fullModel;
    /// SProd::transient.momentoDesesp -- written.
    double &desperationMoment;
    /// SProd::transient.mudaModoChk -- written.
    int &chokeModeChanged;
    /// SProd::transient.mult -- written.
    double &multiplier;
    /// SProd::presfim -- written.
    double &outletPressure;
    /// SProd::transient.reinicia -- written.
    int &restart;
    /// SProd::transient.restriDt -- written.
    int &timeStepRestricted;
    /// SProd::transient.tempoaberto -- written.
    int &openTime;
    /// SProd::transient.termolivreP -- written.
    Vcr<double> &productionSolution;
    /// SProd::titE -- written.
    double &inletQuality;
    /// SProd::transient.vRazMast0 -- written.
    double *masterRatio0;
    /// SProd::transient.vRazMast1 -- written.
    double *masterRatio1;
    /// SProd::transient.vRazMastCrit -- written.
    double *masterCriticalRatio;
    /// SProd::gasLift.velInter -- written.
    double &interfaceVelocity;
    /// SProd::transient.abreM1 -- read only.
    double* masterOpenSchedule;
    /// SProd::arq -- read only.
    Ler &input;
    /// SProd::gasLift.celInterIni -- read only.
    const int &initialInterfaceCell;
    /// SProd::celula -- read only.
    Cel* cells;
    /// SProd::gasLift.celulaG -- read only.
    CelG* gasCells;
    /// SProd::chokeSup -- read only.
    choke &surfaceChoke;
    /// SProd::transient.dtCFL -- read only.
    std::vector<double> &cflTimeSteps;
    /// SProd::gasLift.dtInterIni -- read only.
    const double &initialInterfaceTimeStep;
    /// SProd::transient.dtSim -- read only.
    std::vector<double> &simulationTimeSteps;
    /// SProd::transient.dtauxCFL -- read only.
    double &auxiliaryCflTimeStep;
    /// SProd::transient.dtauxFinal -- read only.
    double &finalAuxiliaryTimeStep;
    /// SProd::transient.fechaM1 -- read only.
    double* masterCloseSchedule;
    /// SProd::flut -- read only.
    FullMtx<double> &productionFreeTerms;
    /// SProd::flutG -- read only.
    FullMtx<double> &gasFreeTerms;
    /// SProd::indTramo -- read only.
    const int &branchIndex;
    /// SProd::transient.jMedMov -- read only.
    double &movingMeanFlux;
    /// SProd::transient.kSP -- read only.
    int &stepIndex;
    /// SProd::transient.matglobP -- read only.
    BandMtx<double> &productionMatrix;
    /// SProd::transient.menorDx -- read only.
    const double &smallestCellLength;
    /// SProd::transient.nabreM1 -- read only.
    const int &masterOpenCount;
    /// SProd::ncel -- read only.
    const int &lastCell;
    /// SProd::gasLift.ncelGas -- read only.
    const int &gasCellCount;
    /// SProd::trends.ncelperftransp -- read only.
    int* productionCrossSectionCount;
    /// SProd::transient.nfechaM1 -- read only.
    const int &masterCloseCount;
    /// SProd::noextremo -- read only.
    const int &endNode;
    /// SProd::gasLift.pGSup -- read only.
    double &gasSurfacePressure;
    /// SProd::presE -- read only.
    double &inletPressure;
    /// SProd::transient.presMedMov -- read only.
    double &movingMeanPressure;
    /// SProd::transient.tMedMov -- read only.
    const double &movingMeanTemperature;
    /// SProd::transient.taxaDTMax -- read only.
    std::vector<double> &maximumTimeStepRates;
    /// SProd::transient.taxaDpMax -- read only.
    std::vector<double> &maximumPressureRates;
    /// SProd::tempE -- read only.
    double &inletTemperature;
    /// SProd::titRev -- read only.
    const double &reverseQuality;
    /// SProd::gasLift.velInterIni -- read only.
    const double &initialInterfaceVelocity;
    /// SProd::vg1dSP -- read only.
    varGlob1D* globals;
    /// Everything the step needs that is not its own.
    TransientStepUpdaters updaters;
};

// ------------------------------------------------------------ cell update ----

/// Updates every cell after a transient solve, in index order, which is
/// load-bearing: the first, interior and last cells each have their own
/// update, called from the one loop.
void updateCells(const TransientStepState &state, int expli = 0);
void updateInteriorCell(const TransientStepState &state, int i, int expli);
void updateFirstCell(const TransientStepState &state, int i, int expli);
void updateLastCell(const TransientStepState &state, int i, int expli);

/// Updates the flow rates. Writes a strict subset of what updateCells writes.
void updateFlowRates(const TransientStepState &state);

// ---------------------------------------------------------- buffer update ----

/// Fills the buffered state from the solver's free-term vector, and from the
/// cells' own current values.
///
/// Two functions, not one: Num4Main.cpp picks between them in alternative
/// branches of one if, and they differ by more than their source:
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
/// The two differ by more than their source: the pressure form computes the
/// Joule-Thomson temperature downstream of the choke and the buffer form does
/// not; the buffer form propagates the mass sources and the pressure form does
/// not.
///
/// SProdVap has its own calcCCpres and calcCCBuffer, taking two arguments
/// instead of three: parallel implementations in another class, not
/// overloads.
void applyOutletPressureCondition(const TransientStepState &state, double titRev, double alfRev, double betRev);
void applyOutletBufferCondition(const TransientStepState &state, double titRev, double alfRev, double betRev);

// ------------------------------------------------------------- time step ----

/// Decides the time step, explicitly or implicitly.
///
/// A last-bit change here does not make a small difference in the answer: it
/// makes a different temporal discretisation, and from that step onward the run
/// is a different simulation.
void computeTimeStep(const TransientStepState &state, int vexpli = 0);
void computeExplicitTimeStep(const TransientStepState &state);
void computeImplicitTimeStep(const TransientStepState &state);

// ------------------------------------------------------ time-step policy ----

/// The restrictions applied to the time step after computeTimeStep proposes it.
/// valveOpeningLow and valveOpeningHigh are the two ends of the valve ramp.
///
/// The step this module runs applies the first three as policies of a registry
/// (TimeStepPolicies in SisProdTransient.cpp), where a new restriction is added
/// by listing it. They are declared here because SProd's methods for the same
/// jobs delegate to them, and Num4Main.cpp calls those methods when it
/// sequences a network run itself.
void dampMaximumTimeStep(const TransientStepState &state);
void evaluatePressureRateOfChange(const TransientStepState &state, double razMast, double razMast0, int vexpli);
void restrictTimeStepByValve(const TransientStepState &state);
void valveOpeningLow(const TransientStepState &state);
void valveOpeningHigh(const TransientStepState &state);

// ---------------------------------------------------- fraction evolution ----

/// Advances the phase fractions in time, and the three ways of restarting that
/// evolution. restartFractionEvolutionInitial and restartFractionEvolution are
/// two functions because Num4Main.cpp selects between them.
void evolveFractions(const TransientStepState &state, double alfrev = 1., double betrev = 0., int ciclo = 0);
void restartFractionEvolutionInitial(const TransientStepState &state);
void restartFractionEvolutionSub(const TransientStepState &state);
void restartFractionEvolution(const TransientStepState &state);

// ------------------------------------------------------ step bookkeeping ----

/// Pig position, the pressure-volume coupling, the fluid mini-table and the
/// inlet condition.
///
/// refreshFluidMiniTable invokes generateFluidMiniTables, below.
void updatePig(const TransientStepState &state);
void solvePressureVolumeCoupling(const TransientStepState &state, int vexpli = 0, int ciclo = 0);
void refreshFluidMiniTable(const TransientStepState &state);
void refreshInletCondition(const TransientStepState &state);

/// The per-cell loops run every time step: recentring the fluid mini-tables,
/// and storing the mass sources at the previous time level.
void generateFluidMiniTables(const TransientStepState &state);
void fillFluidMiniTable(const TransientStepState &state, ProFlu &flui);
void storePreviousSources(const TransientStepState &state);

// =============================================================== the solve ====

/// What SolveTrans needs from outside itself: the SProd methods it calls.
///
/// solveHydrateEnvelopes constructs the hydrate solvers from the whole SProd,
/// which no state struct can supply.
///
/// Default arguments are carried only where a call relies on them:
/// updateThermal is called with none.
struct TransientSolveUpdaters {
    SProd &system;

    void solveHydrateEnvelopes() const;
    [[nodiscard]] double searchUnloadingInjectionPressure() const;
    void writeProductionTrendHeader(int i, int nrede) const;
    void writeProductionTrendRows(int i, int nrede) const;
    void writeGasTrendHeader(int i, int nrede) const;
    void writeGasTrendRows(int i, int nrede) const;
    void writeProductionCrossSectionTrendHeader(int i) const;
    void writeProductionCrossSectionTrendRows(int i) const;
    void writeGasCrossSectionTrendHeader(int i) const;
    void writeGasCrossSectionTrendRows(int i) const;
    void evaluateParaffin() const;
    void connectTubing() const;
    void marchTransientEnergy(int ciclo, int ciclomax) const;
    void updateMolarFractions(const ProFlu &fluiRev) const;
    void updateDensities() const;
    void updateGasOilRatioAndCo2(const ProFlu &fluiRev) const;
    void updateTemperatures() const;
    void updateInitialFractions() const;
    void updateThermal(int aflu = 0) const;
    void solveGasLine() const;
};

/// The state SolveTrans reads.
///
/// It composes the step state rather than restating it, as
/// SteadyStateSearchState does: the fields below are the members SolveTrans
/// touches that none of the step's routines do.
///
/// None of them is const: writes through a field and by reference do not look
/// like assignments to it. Pointers are held by reference (`T *&`), so the
/// state stays exact even if the solve reseats one.
struct TransientSolveState {
    /// Everything the step reads. SolveTrans hands this to every step routine.
    TransientStepState step;

    /// SProd::temperatura -- written or read by the solve; not promised const.
    double &defaultInletTemperature;
    /// SProd::derivaAnel -- written or read by the solve; not promised const.
    int &annulusDrift;
    /// SProd::saidaSubTextoSis -- read only.
    const char *const *closingSubtitles;
    /// SProd::saidaTextoSis -- read only.
    const char *const *closingTitles;
    /// SProd::transient.kontaRenovaComp -- written or read by the solve; not promised const.
    int &compositionalRefreshCounter;
    /// SProd::transient.jVet -- written or read by the solve; not promised const.
    std::vector<double> &fluxHistory;
    /// SProd::trends.ncelperftransg -- written or read by the solve; not promised const.
    int* gasCrossSectionCellCounts;
    /// SProd::trends.kontaTempoTransProfG -- written or read by the solve; not promised const.
    int &gasCrossSectionProfileTimeCounter;
    /// SProd::trends.ntrendtransgB -- written or read by the solve; not promised const.
    int* &gasCrossSectionTrendBufferedCounts;
    /// SProd::trends.ntrendtransg -- written or read by the solve; not promised const.
    int* &gasCrossSectionTrendCounts;
    /// SProd::trends.MatTrendTransG -- written or read by the solve; not promised const.
    double*** &gasCrossSectionTrendMatrix;
    /// SProd::trends.resettrendtransg -- written or read by the solve; not promised const.
    double* &gasCrossSectionTrendResetTimers;
    /// SProd::trends.kontaTempoProfG -- written or read by the solve; not promised const.
    int &gasProfileTimeCounter;
    /// SProd::trends.ntrendgB -- written or read by the solve; not promised const.
    int* &gasTrendBufferedCounts;
    /// SProd::trends.ntrendg -- written or read by the solve; not promised const.
    int* &gasTrendCounts;
    /// SProd::trends.MatTrendG -- written or read by the solve; not promised const.
    double*** &gasTrendMatrix;
    /// SProd::trends.resettrendg -- written or read by the solve; not promised const.
    double* &gasTrendResetTimers;
    /// SProd::gasLift.presiniG -- written or read by the solve; not promised const.
    double &initialGasPressure;
    /// SProd::transient.pGSupIni -- written or read by the solve; not promised const.
    double &initialGasSurfacePressure;
    /// SProd::gasLift.tempiniG -- written or read by the solve; not promised const.
    double &initialGasTemperature;
    /// SProd::transient.tempoabertoini -- written or read by the solve; not promised const.
    int &initialOpenTime;
    /// SProd::gasLift.chokeInj -- written or read by the solve; not promised const.
    ChokeGas &injectionChoke;
    /// SProd::tmpLog -- written or read by the solve; not promised const.
    std::string &logBuffer;
    /// SProd::transient.contaLog -- written or read by the solve; not promised const.
    int &logCounter;
    /// SProd::transient.TransMassModel -- written or read by the solve; not promised const.
    int &massTransferModel;
    /// SProd::transient.dtCicMin -- written or read by the solve; not promised const.
    double &minimumCycleTimeStep;
    /// SProd::transient.ktMedMov -- written or read by the solve; not promised const.
    double &movingMeanCounter;
    /// SProd::transient.alfMedMov -- written or read by the solve; not promised const.
    double &movingMeanVoidFraction;
    /// SProd::gasLift.verificaAcop -- written or read by the solve; not promised const.
    int &networkCoupled;
    /// SProd::poisson3D -- written or read by the solve; not promised const.
    solverP3D &poissonSolver3D;
    /// SProd::transient.presVet -- written or read by the solve; not promised const.
    std::vector<double> &pressureHistory;
    /// SProd::transient.KontaImprime -- written or read by the solve; not promised const.
    int &printCounter;
    /// SProd::trends.kimpT -- written or read by the solve; not promised const.
    double &printPassCount;
    /// SProd::trends.kontaTempoTransProf -- written or read by the solve; not promised const.
    int &productionCrossSectionProfileTimeCounter;
    /// SProd::trends.ntrendtransB -- written or read by the solve; not promised const.
    int* &productionCrossSectionTrendBufferedCounts;
    /// SProd::trends.ntrendtrans -- written or read by the solve; not promised const.
    int* &productionCrossSectionTrendCounts;
    /// SProd::trends.MatTrendTransP -- written or read by the solve; not promised const.
    double*** &productionCrossSectionTrendMatrix;
    /// SProd::trends.resettrendtrans -- written or read by the solve; not promised const.
    double* &productionCrossSectionTrendResetTimers;
    /// SProd::trends.kontaTempoProf -- written or read by the solve; not promised const.
    int &productionProfileTimeCounter;
    /// SProd::trends.ntrendB -- written or read by the solve; not promised const.
    int* &productionTrendBufferedCounts;
    /// SProd::trends.ntrend -- written or read by the solve; not promised const.
    int* &productionTrendCounts;
    /// SProd::trends.MatTrendP -- written or read by the solve; not promised const.
    double*** &productionTrendMatrix;
    /// SProd::trends.resettrend -- written or read by the solve; not promised const.
    double* &productionTrendResetTimers;
    /// SProd::noinicial -- written or read by the solve; not promised const.
    int &startNode;
    /// SProd::transient.tVet -- written or read by the solve; not promised const.
    std::vector<double> &temperatureHistory;
    /// SProd::transient.jTotal -- written or read by the solve; not promised const.
    double &totalFlux;
    /// SProd::transient.pTotal -- written or read by the solve; not promised const.
    double &totalPressure;
    /// SProd::transient.alfTotal -- written or read by the solve; not promised const.
    double &totalVoidFraction;
    /// SProd::transient.trackDeng -- written or read by the solve; not promised const.
    int &trackGasGravity;
    /// SProd::transient.trackRGO -- written or read by the solve; not promised const.
    int &trackGasOilRatio;
    /// SProd::trends.kontaTempoCelUni -- written or read by the solve; not promised const.
    std::vector<int> &unitCellTimeCounters;
    /// SProd::transient.alfVet -- written or read by the solve; not promised const.
    std::vector<double> &voidFractionHistory;

    /// Everything the solve needs that is not its own.
    TransientSolveUpdaters updaters;
};

/// One transient step, in six parts that run in this order: hydrates,
/// mini-table, the t = 0 trend, computeTimeStep, the coupling loop (fraction
/// evolution, pig, pressure-volume coupling), then the output phases.
void solveTransientStep(const TransientSolveState &state, double titRev, double alfRev, double betRev,
                        int nrede, ProFlu fluiRev);

}  // namespace sisprod::transient

#endif  // SISPRODTRANSIENT_H_
