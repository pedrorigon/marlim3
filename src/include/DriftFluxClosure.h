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
}  // namespace driftflux

#endif  // DRIFTFLUXCLOSURE_H_
