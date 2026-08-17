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
void BhagwatGhajar(double rhol, double rhog, double tensup, double alf, double reymix, double reymixL,
                   double ug1, double ul1, double dia, double rug, double tet, double &c0,
                   double &ud, double correcHor);

/// Evaluates C0 and Ud using the Bhagwat-Ghajar drift-flux correlation.
void BhagwatGhajarMod(double rhol, double rhog, double tensup, double alf, double reymix, double reymixL,
                      double ug1, double ul1, double dia, double rug, double tet, double &c0,
                      double &ud, double correcHor);

/// Evaluates C0 and Ud using the Choi drift-flux correlation.
void Choi(double rhol, double rhog, double tensup, double alf, double reymix, double reymixL,
          double ug1, double ul1, double dia, double rug, double tet, double &c0,
          double &ud, double correcHor);

/// Evaluates C0 and Ud using the Hibiki-Ishii drift-flux correlation.
void HibikiIshii(double rhol, double rhog, double tensup, double alf, double reymix, double reymixL,
                 double ug1, double ul1, double dia, double rug, double tet, double &c0,
                 double &ud, double correcHor);

/// Evaluates C0 and Ud using the Franca-Lahey drift-flux correlation.
void FrancaLahey(double rhol, double rhog, double tensup, double alf, double reymix, double reymixL,
                 double ug1, double ul1, double dia, double rug, double tet, double &c0,
                 double &ud, double correcHor);

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

/// Evaluates C0 and Ud for dispersed flow, using the already-resolved
/// correlationIndex (arq.CorreDisper at the call site) instead of reading the
/// configuration itself.
void C0UdDisperso(double rhol, double rhog, double tensup, double alf, double reymix, double reymixL,
                  double ug1, double ul1, double dia, double rug, double tet, double &c0,
                  double &ud, double correcHor, int estabCol, int correlationIndex);

/// Evaluates C0 and Ud for annular or churn flow, using the already-resolved
/// correlationIndex (arq.CorreAnular at the call site).
void C0UdAnularChurn(double rhol, double rhog, double tensup, double alf, double reymix, double reymixL,
                     double ug1, double ul1, double dia, double rug, double tet, double &c0,
                     double &ud, double correcHor, int estabCol, int correlationIndex);

/// Evaluates C0 and Ud for stratified flow, using the already-resolved
/// correlationIndex (arq.CorreEstrat at the call site).
void C0UdEstratificado(double rhol, double rhog, double tensup, double alf, double reymix, double reymixL,
                       double ug1, double ul1, double dia, double rug, double tet, double &c0,
                       double &ud, double correcHor, int estabCol, int correlationIndex);

}  // namespace correlations
}  // namespace driftflux

#endif  // DRIFTFLUXCLOSURE_H_
