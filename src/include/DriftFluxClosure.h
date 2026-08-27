#ifndef DRIFTFLUXCLOSURE_H_
#define DRIFTFLUXCLOSURE_H_

/// Drift-flux closure relations.
///
/// The functions in driftflux::correlations are pure: they read nothing but
/// their arguments and write nothing but the c0 and ud output references. That
/// property is what makes them verifiable in isolation against a tabulated
/// sweep, which matters here more than usual -- three of the five correlations
/// are never reached by the demo corpus, so running the models proves nothing
/// about them.
// Declared, not included: ClosureState below holds only pointers and references
// to these, so this header stays free of the cell, input-deck and globals
// headers and can still be compiled on its own.
class Cel;
class Ler;
struct varGlob1D;

namespace driftflux {
namespace correlations {

/// Evaluates C0 and Ud using the Bhagwat-Ghajar drift-flux correlation.
void BhagwatGhajar(double liquidDensity, double gasDensity, double surfaceTension, double voidFraction,
                   double mixtureReynolds, double liquidReynolds, double gasFlowRate, double liquidFlowRate,
                   double diameter, double roughness, double inclinationAngle, double &c0, double &ud,
                   double horizontalCorrection);

/// Evaluates C0 and Ud using the Bhagwat-Ghajar drift-flux correlation.
void BhagwatGhajarMod(double liquidDensity, double gasDensity, double surfaceTension, double voidFraction,
                      double mixtureReynolds, double liquidReynolds, double gasFlowRate,
                      double liquidFlowRate, double diameter, double roughness, double inclinationAngle,
                      double &c0, double &ud, double horizontalCorrection);

/// Evaluates C0 and Ud using the Choi drift-flux correlation.
void Choi(double liquidDensity, double gasDensity, double surfaceTension, double voidFraction,
          double mixtureReynolds, double liquidReynolds, double gasFlowRate, double liquidFlowRate,
          double diameter, double roughness, double inclinationAngle, double &c0, double &ud,
          double horizontalCorrection);

/// Evaluates C0 and Ud using the Hibiki-Ishii drift-flux correlation.
void HibikiIshii(double liquidDensity, double gasDensity, double surfaceTension, double voidFraction,
                 double mixtureReynolds, double liquidReynolds, double gasFlowRate, double liquidFlowRate,
                 double diameter, double roughness, double inclinationAngle, double &c0, double &ud,
                 double horizontalCorrection);

/// Evaluates C0 and Ud using the Franca-Lahey drift-flux correlation.
void FrancaLahey(double liquidDensity, double gasDensity, double surfaceTension, double voidFraction,
                 double mixtureReynolds, double liquidReynolds, double gasFlowRate, double liquidFlowRate,
                 double diameter, double roughness, double inclinationAngle, double &c0, double &ud,
                 double horizontalCorrection);

/// Correlation selected by the aggregators below, resolved from the integer
/// configuration field before the call.
///
/// The three aggregators accept different subsets of these values, and each
/// leaves c0 and ud untouched when the selector falls outside its own subset.
/// That is observable behaviour of the original switch statements, so the
/// subsets are kept apart rather than merged into one table.
///
///   0  Choi              accepted by all three
///   1  BhagwatGhajar     accepted by all three
///   2  FrancaLahey       stratified only
///   3  HibikiIshii       annular/churn only
///   4  BhagwatGhajarMod  accepted by all three
///   5  angle blend       accepted by all three

/// The correlation each flow regime selects, read from configuration once per
/// run rather than on every cell.
///
/// The three fields are constant for a whole simulation: the engine parses them
/// from the input file and never writes them again. Holding them together says
/// so, and lets the per cell path stop consulting configuration, which is what
/// the performance requirement asks for.
///
/// The values are the ones documented above. Each regime accepts its own subset
/// and ignores the rest, so the fields are not interchangeable.
struct RegimeSelectors {
    int dispersed;    ///< arq.CorreDisper, accepts 0, 1, 4 and 5.
    int annularChurn; ///< arq.CorreAnular, accepts 3 as well.
    int stratified;   ///< arq.CorreEstrat, accepts 2 as well.
};

/// Evaluates C0 and Ud for dispersed flow, using the already-resolved
/// correlationIndex (arq.CorreDisper at the call site) instead of reading the
/// configuration itself.
void C0UdDisperso(double liquidDensity, double gasDensity, double surfaceTension, double voidFraction,
                  double mixtureReynolds, double liquidReynolds, double gasFlowRate, double liquidFlowRate,
                  double diameter, double roughness, double inclinationAngle, double &c0, double &ud,
                  double horizontalCorrection, int estabCol, int correlationIndex);

/// Evaluates C0 and Ud for annular or churn flow, using the already-resolved
/// correlationIndex (arq.CorreAnular at the call site).
void C0UdAnularChurn(double liquidDensity, double gasDensity, double surfaceTension, double voidFraction,
                     double mixtureReynolds, double liquidReynolds, double gasFlowRate,
                     double liquidFlowRate, double diameter, double roughness, double inclinationAngle,
                     double &c0, double &ud, double horizontalCorrection, int estabCol, int correlationIndex);

/// Evaluates C0 and Ud for stratified flow, using the already-resolved
/// correlationIndex (arq.CorreEstrat at the call site).
void C0UdEstratificado(double liquidDensity, double gasDensity, double surfaceTension, double voidFraction,
                       double mixtureReynolds, double liquidReynolds, double gasFlowRate,
                       double liquidFlowRate, double diameter, double roughness, double inclinationAngle,
                       double &c0, double &ud, double horizontalCorrection, int estabCol,
                       int correlationIndex);

}  // namespace correlations

/// Distribution coefficient: the five variants of CalcC0Ud.
///
/// These are not pure. They read celula[] intensively, they write back into it
/// -- arranjo, arranjoR, transic, transic0, perdaEstratL/G, c0Spare, udSpare --
/// and they call the correlations above. What they do NOT do is call any other
/// method of SProd or touch `this`, which is what makes moving them possible at
/// all (verified: the only calls in the five bodies are fabs, pow and the two
/// local flow-pattern map objects).
///
/// The five differ far more than their names suggest. A normalized comparison
/// of the two closest non-initialisation variants found 38 divergence sites, of
/// which 21 are control flow and only 7 are the data source that the original
/// design expected to be the whole story: one builds the flow-pattern map and
/// writes the pattern back, another reads the pattern it was given; one runs the
/// transition counter, another does not. Measurement and the full table are in
/// specs/001-refatoracao-sisprod/evidencia/c0ud-diff.md.
///
/// So the five keep their own control flow, and only the blocks a normalized
/// comparison proved identical are shared. Sharing more would mean either a
/// runtime mode test inside a shared core, or a policy whose hooks have one
/// caller each -- the original function wearing a hat.
namespace coefficient {

/// The state the five read, and the only state they may read.
///
/// Same role as TrendState in the trend module: it names in one place what a
/// closure evaluation is allowed to touch, and it makes the five callable
/// WITHOUT an SProd. That second property is not tidiness -- it is what lets a
/// dedicated harness drive them over synthetic cells, which is the only
/// verification that reaches CalcC0UdBuf, CalcC0UdIni and CalcC0UdIniBuf. The
/// demo corpus never executes those three, so the artifact and regression
/// layers are green for them whatever happens.
///
/// The scalars are held BY REFERENCE, not by value. Copying them into the
/// struct would read every one of them at construction, before the branch that
/// decides whether the original would have read it at all. Nothing writes them
/// during a call, so a copy would give the same numbers here -- but "the same
/// numbers today" is how a conditional read silently becomes unconditional, and
/// that is precisely the defect the root-finding stage had to undo.
struct ClosureState {
    /// The cell array -- SProd::celula. Written as well as read.
    Cel *cells;
    /// Index of the last cell -- SProd::ncel.
    const int &lastCell;
    /// Shared 1D globals; only localtiny is read, plus the pointer handed to
    /// the flow-pattern map -- SProd::vg1dSP.
    varGlob1D *globals;
    /// Input deck -- SProd::arq. Read for mapaArranjo, escorregaTran,
    /// escorregaPerm and AceleraConvergPerm.
    const Ler &input;
    /// Correlation choice per regime, resolved once per run -- see
    /// RegimeSelectors above and SProd::driftSelectors.
    const correlations::RegimeSelectors &selectors;

    /// Gas temperature at the surface -- SProd::tGSup.
    const double &gasSurfaceTemperature;
    /// Inlet void fraction -- SProd::alfE. Read by the initialisation variants.
    const double &inletVoidFraction;
    /// Inlet column fraction -- SProd::betaE.
    const double &inletColumnFraction;
    /// Inlet pressure -- SProd::presE.
    const double &inletPressure;
    /// Inlet temperature -- SProd::tempE.
    const double &inletTemperature;
    /// Steady-state iteration counter -- SProd::iterperm.
    const int &steadyIteration;
};

/// Slip parameters at a production-line face, from the instantaneous state.
/// Was SProd::CalcC0Ud.
void instantaneous(const ClosureState &state, int ind, double &c0, double &ud);

/// Slip parameters for the intermediate network state, from the buffered
/// fields. Was SProd::CalcC0UdBuf.
void buffered(const ClosureState &state, int ind, double &c0, double &ud);

/// Slip parameters at the inlet of an internal network section.
/// Was SProd::CalcC0UdIni.
void initialization(const ClosureState &state, int ind, double &c0, double &ud);

/// Slip parameters at the inlet of an internal network section, from the
/// buffered fields. Was SProd::CalcC0UdIniBuf.
void bufferedInitialization(const ClosureState &state, int ind, double &c0, double &ud);

/// Steady-state slip parameters at a downstream face.
/// Was SProd::CalcC0UdPerm.
void steadyState(const ClosureState &state, int ind, double &c0, double &ud);

}  // namespace coefficient
}  // namespace driftflux

#endif  // DRIFTFLUXCLOSURE_H_
