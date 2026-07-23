/*
 * SisProd.cpp
 *
 * Implementation of the new SisProd architecture.
 *
 * This translation unit holds the behaviour declared in SisProd2.h. The legacy
 * implementation is preserved intact in SisProd_old.cpp and remains the
 * production solver; the build option MARLIM_USE_NEW_SISPROD selects which
 * implementation is active (default: legacy), so every existing flow and
 * regression test keeps passing while the new solver is completed.
 *
 * Nothing here is invoked by the legacy flow yet, so compiling it into the
 * binary does not change the behaviour or the numerical results of Marlim3.
 */

#include "SisProd2.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <ostream>
#include <sstream>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace marlim {
namespace sisprod2 {

// ===========================================================================
// Helpers and implementation-selection flag.
// ===========================================================================

bool almostEqual(double a, double b, double relativeTolerance) {
    return std::fabs(a - b) <=
           relativeTolerance * (1.0 + std::fabs(a) + std::fabs(b));
}

bool usingNewSisProd() {
#ifdef MARLIM_USE_NEW_SISPROD
    return true;
#else
    return false;
#endif
}

const char *activeSisProdImplementation() {
    return usingNewSisProd() ? "new" : "legacy";
}

const char *toString(SolveStatus status) {
    switch (status) {
    case SolveStatus::Ok:
        return "Ok";
    case SolveStatus::MaxIterations:
        return "MaxIterations";
    case SolveStatus::OutOfRange:
        return "OutOfRange";
    case SolveStatus::NotBracketed:
        return "NotBracketed";
    case SolveStatus::NotANumber:
        return "NotANumber";
    }
    return "Unknown";
}

const char *toString(LogLevel level) {
    switch (level) {
    case LogLevel::Debug:
        return "DEBUG";
    case LogLevel::Info:
        return "INFO";
    case LogLevel::Warning:
        return "WARNING";
    case LogLevel::Error:
        return "ERROR";
    }
    return "INFO";
}

// ===========================================================================
// Diagnostics.
// ===========================================================================

OstreamDiagnostics::OstreamDiagnostics(std::ostream &out)
    : out_(out), failureCount_(0) {}

void OstreamDiagnostics::log(const LogRecord &record) {
    out_ << "[" << toString(record.level) << "]"
         << " op=" << record.operation << " section=" << record.sectionIndex
         << " t=" << record.time << " msg=" << record.message << "\n";
}

void OstreamDiagnostics::reportSolve(const std::string &operation,
                                     const SolveResult &result) {
    if (!result.ok())
        ++failureCount_;
    LogRecord record;
    record.level = result.ok() ? LogLevel::Info : LogLevel::Warning;
    record.operation = operation;
    std::ostringstream message;
    message << "status=" << toString(result.status) << " value=" << result.value
            << " iterations=" << result.iterations;
    if (!result.detail.empty())
        message << " detail=" << result.detail;
    record.message = message.str();
    log(record);
}

// ===========================================================================
// Root finders.
// ===========================================================================

SolveResult RootFinder::solveBracketed(const ResidualFunction &residual,
                                       double lowerBound, double upperBound,
                                       double tolerance, int maxIterations) {
    SolveResult result;

    double a = lowerBound;
    double b = upperBound;
    double c = upperBound;
    double d = 0.0;
    double e = 0.0;
    double fa = residual(a);
    double fb = residual(b);
    double fc = fb;

    if (std::isnan(fa) || std::isnan(fb)) {
        result.status = SolveStatus::NotANumber;
        result.detail = "residual is NaN at the bracket bounds";
        return result;
    }
    if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)) {
        result.status = SolveStatus::NotBracketed;
        result.value = std::fabs(fa) < std::fabs(fb) ? a : b;
        result.detail = "the bracket does not enclose a sign change";
        return result;
    }

    const double epsilon = std::numeric_limits<double>::epsilon();
    for (int iteration = 1; iteration <= maxIterations; ++iteration) {
        if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
            c = a;
            fc = fa;
            e = d = b - a;
        }
        if (std::fabs(fc) < std::fabs(fb)) {
            a = b;
            b = c;
            c = a;
            fa = fb;
            fb = fc;
            fc = fa;
        }

        const double convergenceTol = 2.0 * epsilon * std::fabs(b) + 0.5 * tolerance;
        const double midpoint = 0.5 * (c - b);
        if (std::fabs(midpoint) <= convergenceTol || fb == 0.0) {
            result.status = SolveStatus::Ok;
            result.value = b;
            result.iterations = iteration;
            return result;
        }

        if (std::fabs(e) >= convergenceTol && std::fabs(fa) > std::fabs(fb)) {
            const double s = fb / fa;
            double p;
            double q;
            if (a == c) {
                p = 2.0 * midpoint * s;
                q = 1.0 - s;
            } else {
                const double qa = fa / fc;
                const double qb = fb / fc;
                p = s * (2.0 * midpoint * qa * (qa - qb) - (b - a) * (qb - 1.0));
                q = (qa - 1.0) * (qb - 1.0) * (s - 1.0);
            }
            if (p > 0.0)
                q = -q;
            p = std::fabs(p);
            const double limit1 = 3.0 * midpoint * q - std::fabs(convergenceTol * q);
            const double limit2 = std::fabs(e * q);
            if (2.0 * p < (limit1 < limit2 ? limit1 : limit2)) {
                e = d;
                d = p / q;
            } else {
                d = midpoint;
                e = d;
            }
        } else {
            d = midpoint;
            e = d;
        }

        a = b;
        fa = fb;
        if (std::fabs(d) > convergenceTol)
            b += d;
        else
            b += (midpoint >= 0.0 ? std::fabs(convergenceTol)
                                  : -std::fabs(convergenceTol));
        fb = residual(b);
        if (std::isnan(fb)) {
            result.status = SolveStatus::NotANumber;
            result.value = b;
            result.iterations = iteration;
            result.detail = "residual became NaN during the iteration";
            return result;
        }
    }

    result.status = SolveStatus::MaxIterations;
    result.value = b;
    result.iterations = maxIterations;
    result.detail = "maximum number of iterations reached";
    return result;
}

SolveResult RootFinder::solveBisection(const ResidualFunction &residual,
                                       double lowerBound, double upperBound,
                                       double tolerance, int maxIterations) {
    SolveResult result;
    double a = lowerBound;
    double b = upperBound;
    double fa = residual(a);
    double fb = residual(b);

    if (std::isnan(fa) || std::isnan(fb)) {
        result.status = SolveStatus::NotANumber;
        result.detail = "residual is NaN at the bracket bounds";
        return result;
    }
    if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)) {
        result.status = SolveStatus::NotBracketed;
        result.detail = "the bracket does not enclose a sign change";
        return result;
    }

    for (int iteration = 1; iteration <= maxIterations; ++iteration) {
        const double midpoint = 0.5 * (a + b);
        const double fm = residual(midpoint);
        if (std::isnan(fm)) {
            result.status = SolveStatus::NotANumber;
            result.value = midpoint;
            result.iterations = iteration;
            return result;
        }
        if (0.5 * (b - a) <= tolerance || fm == 0.0) {
            result.status = SolveStatus::Ok;
            result.value = midpoint;
            result.iterations = iteration;
            return result;
        }
        if ((fa > 0.0 && fm > 0.0) || (fa < 0.0 && fm < 0.0)) {
            a = midpoint;
            fa = fm;
        } else {
            b = midpoint;
        }
    }
    result.status = SolveStatus::MaxIterations;
    result.value = 0.5 * (a + b);
    result.iterations = maxIterations;
    return result;
}

// ===========================================================================
// Stream mixing (Strategy).
// ===========================================================================

StandardStream blend(const StandardStream &upstream,
                     const StandardStream &source, double tiny) {
    StandardStream mixed;
    mixed.oilRate = upstream.oilRate + source.oilRate;
    mixed.waterRate = upstream.waterRate + source.waterRate;
    mixed.gasRate = upstream.gasRate + source.gasRate;

    if (mixed.oilRate > tiny) {
        const double blendedRelDensity =
            (StandardStream::apiToRelativeDensity(upstream.oilApi) * upstream.oilRate +
             StandardStream::apiToRelativeDensity(source.oilApi) * source.oilRate) /
            mixed.oilRate;
        mixed.oilApi = StandardStream::relativeDensityToApi(blendedRelDensity);
    } else {
        mixed.oilApi = upstream.oilApi;
    }

    if (mixed.gasRate > tiny) {
        mixed.gasDensity = (upstream.gasRate * upstream.gasDensity +
                            source.gasRate * source.gasDensity) /
                           mixed.gasRate;
        mixed.co2Fraction = (upstream.gasRate * upstream.co2Fraction +
                             source.gasRate * source.co2Fraction) /
                            mixed.gasRate;
    } else {
        mixed.gasDensity = upstream.gasDensity;
        mixed.co2Fraction = upstream.co2Fraction;
    }

    if (mixed.waterRate > tiny) {
        mixed.waterDensity = (upstream.waterRate * upstream.waterDensity +
                              source.waterRate * source.waterDensity) /
                             mixed.waterRate;
    } else {
        mixed.waterDensity = upstream.waterDensity;
    }

    return mixed;
}

StandardStream StreamMixing::march(const StandardStream &upstream,
                                   AccessoryType type,
                                   const StandardStream &source,
                                   const SimContext &context) {
    switch (type) {
    case AccessoryType::None:
        return upstream;
    case AccessoryType::GasSource:
    case AccessoryType::LiquidSource:
    case AccessoryType::Ipr:
    case AccessoryType::MassSource:
        return blend(upstream, source, context.localTiny());
    }
    return upstream;
}

// ===========================================================================
// Fluid model.
// ===========================================================================

double FluidModel::gasDensityAt(double pressure) const {
    return standardGasDensity_ * (pressure / standardPressure_);
}

double FluidModel::massQualityFromGasOilRatio(double gasOilRatio) const {
    const double gasMass = gasOilRatio * standardGasDensity_;
    const double liquidMass = liquidDensity_;
    const double total = gasMass + liquidMass;
    return total > 0.0 ? gasMass / total : 0.0;
}

double FluidModel::mixtureDensity(double pressure, double massQuality) const {
    const double clampedQuality =
        massQuality < 0.0 ? 0.0 : (massQuality > 1.0 ? 1.0 : massQuality);
    const double gasDensity = gasDensityAt(pressure);
    if (gasDensity <= 0.0)
        return liquidDensity_;
    const double specificVolume =
        clampedQuality / gasDensity + (1.0 - clampedQuality) / liquidDensity_;
    return specificVolume > 0.0 ? 1.0 / specificVolume : liquidDensity_;
}

double darcyFrictionFactor(double reynolds) {
    if (reynolds < 1.0)
        return 0.0;
    if (reynolds < 2300.0)
        return 64.0 / reynolds;
    return 0.3164 * std::pow(reynolds, -0.25);
}

// ===========================================================================
// Production column.
// ===========================================================================

double PipeSegment::verticalRise() const {
    return length * std::sin(inclination);
}

double ProductionColumn::totalVerticalHeight() const {
    double height = 0.0;
    for (std::size_t i = 0; i < segments_.size(); ++i)
        height += segments_[i].verticalRise();
    return height;
}

double ProductionColumn::marchToWellhead(double bottomholePressure,
                                         const SimContext &context,
                                         std::vector<ProfilePoint> *profileOut) const {
    double pressure = bottomholePressure;
    StandardStream stream = inletStream_;

    if (profileOut) {
        profileOut->clear();
        profileOut->reserve(segments_.size());
    }

    for (std::size_t i = 0; i < segments_.size(); ++i) {
        const PipeSegment &segment = segments_[i];

        if (segment.accessory != AccessoryType::None) {
            stream = StreamMixing::march(stream, segment.accessory,
                                         segment.sourceStream, context);
        }

        const double quality =
            fluid_.massQualityFromGasOilRatio(stream.gasOilRatio());
        const double densityPressure =
            pressure > constants::kPressureFloor ? pressure
                                                 : constants::kPressureFloor;
        const double density = fluid_.mixtureDensity(densityPressure, quality);
        const double area = segment.crossSectionArea();
        const double velocity = massFlowRate_ / (density * area);

        const double reynolds =
            density * velocity * segment.diameter / fluid_.viscosity();
        const double friction = darcyFrictionFactor(reynolds);

        const double hydrostaticDrop =
            density * constants::kGravity * segment.verticalRise();
        const double frictionalDrop = friction *
                                      (segment.length / segment.diameter) * 0.5 *
                                      density * velocity * velocity;

        if (profileOut) {
            ProfilePoint point;
            point.pressure = pressure;
            point.mixtureDensity = density;
            point.velocity = velocity;
            point.gasOilRatio = stream.gasOilRatio();
            profileOut->push_back(point);
        }

        pressure -= hydrostaticDrop + frictionalDrop;
        if (!std::isfinite(pressure))
            return pressure;
    }

    return pressure;
}

// ===========================================================================
// Steady-state helpers and facade.
// ===========================================================================

ResidualFunction makeWellheadResidual(const ProductionColumn &column,
                                      const SimContext &context,
                                      double separatorPressure) {
    return [&column, &context, separatorPressure](double bottomhole) {
        return column.marchToWellhead(bottomhole, context) - separatorPressure;
    };
}

void wellheadSearchBounds(const ProductionColumn &column,
                          const SteadyStateRequest &request, double &lowerBound,
                          double &upperBound) {
    const double hydrostaticGuess = column.fluid().liquidDensity() *
                                    constants::kGravity *
                                    column.totalVerticalHeight();
    // The wellhead pressure grows monotonically with the bottomhole pressure.
    // The separator pressure is a safe lower bound (the wellhead cannot exceed
    // the bottomhole), and the hydrostatic guess scaled by searchMargin gives a
    // safe upper bound (the real drop is always below the pure-liquid column).
    lowerBound = request.separatorPressure;
    upperBound =
        request.separatorPressure + request.searchMargin * hydrostaticGuess;
}

SolveResult TramoEngine::solveBracketed(const std::string &operation,
                                        const ResidualFunction &residual,
                                        double lowerBound, double upperBound,
                                        double tolerance, int maxIterations) {
    SolveResult result = RootFinder::solveBracketed(residual, lowerBound,
                                                    upperBound, tolerance,
                                                    maxIterations);
    diagnostics_.reportSolve(operation, result);
    return result;
}

SolveResult TramoEngine::solveBottomholePressure(const ProductionColumn &column,
                                                 const SteadyStateRequest &request) {
    ResidualFunction residual =
        makeWellheadResidual(column, context_, request.separatorPressure);
    double lowerBound = 0.0;
    double upperBound = 0.0;
    wellheadSearchBounds(column, request, lowerBound, upperBound);

    SolveResult result = RootFinder::solveBracketed(residual, lowerBound,
                                                    upperBound, request.tolerance,
                                                    request.maxIterations);
    diagnostics_.reportSolve("solveBottomholePressure", result);
    return result;
}

StandardStream TramoEngine::marchMass(const StandardStream &upstream,
                                      AccessoryType type,
                                      const StandardStream &source) const {
    return StreamMixing::march(upstream, type, source, context_);
}

// ===========================================================================
// Batch solver with OpenMP (parallel across independent columns).
// ===========================================================================

std::vector<SolveResult>
solveBatch(const std::vector<ProductionColumn> &columns,
           const std::vector<SteadyStateRequest> &requests,
           const SimContext &context) {
    const std::size_t count = columns.size();
    std::vector<SolveResult> results(count);

    // Shared read-only fallback used when a single request applies to every
    // column, or when no request was supplied at all.
    const SteadyStateRequest defaultRequest;
    const int threads = context.threadCount();
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static) if (count > 1)
#endif
    for (long i = 0; i < static_cast<long>(count); ++i) {
        const std::size_t index = static_cast<std::size_t>(i);
        const SteadyStateRequest &request =
            (requests.size() == count)
                ? requests[index]
                : (requests.empty() ? defaultRequest : requests.front());
        ResidualFunction residual = makeWellheadResidual(
            columns[index], context, request.separatorPressure);
        double lowerBound = 0.0;
        double upperBound = 0.0;
        wellheadSearchBounds(columns[index], request, lowerBound, upperBound);
        results[index] = RootFinder::solveBracketed(residual, lowerBound,
                                                     upperBound, request.tolerance,
                                                     request.maxIterations);
    }
    (void)threads; // silence unused warning when compiled without OpenMP
    return results;
}

// ===========================================================================
// Scenario builder.
// ===========================================================================

ProductionColumn buildVerticalWell(const FluidModel &fluid, double massFlowRate,
                                   std::size_t segmentCount, double segmentLength,
                                   double diameter, double inletGasOilRatio) {
    ProductionColumn column(fluid, massFlowRate);
    column.reserveSegments(segmentCount);

    StandardStream inlet;
    inlet.oilRate = 1.0;
    inlet.waterRate = 0.25;
    inlet.gasRate = inletGasOilRatio; // gasRate / oilRate == inletGasOilRatio
    inlet.oilApi = 28.0;
    inlet.gasDensity = 0.75;
    inlet.co2Fraction = 0.02;
    inlet.waterDensity = 1.02;
    column.setInletStream(inlet);

    for (std::size_t i = 0; i < segmentCount; ++i) {
        PipeSegment segment;
        segment.length = segmentLength;
        segment.diameter = diameter;
        segment.inclination = constants::kPi / 2.0; // vertical
        column.addSegment(segment);
    }
    return column;
}

// ===========================================================================
// Automatic comparison.
// ===========================================================================

// ---------------------------------------------------------------------------
// R01 reference: pure Brent (NR algorithm), single-lambda variant.
// Independent re-implementation used only in the R01 characterization
// tests.  This verifies that RootFinder::solveBracketed (the candidate)
// is algorithmically identical to the NR Brent used in the legacy
// zbrent / SProd::zbrent.  Both are deterministic for the same inputs,
// so they should agree to within a few ULPs on any well-behaved function.
// Kept in the anonymous namespace so it does not pollute the public API.
// ---------------------------------------------------------------------------
namespace {

static double brentReference(const ResidualFunction &f,
                             double a, double b,
                             double tol, int maxit) {
    const double eps = std::numeric_limits<double>::epsilon();
    double c = b;
    double fa = f(a), fb = f(b), fc = fb;
    double d = 0.0, e = 0.0;
    if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
        return 0.5 * (a + b); // not bracketed: return midpoint
    for (int iter = 0; iter < maxit; ++iter) {
        if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
            c = a; fc = fa; e = d = b - a;
        }
        if (std::fabs(fc) < std::fabs(fb)) {
            a = b; b = c; c = a; fa = fb; fb = fc; fc = fa;
        }
        const double tol1 = 2.0 * eps * std::fabs(b) + 0.5 * tol;
        const double xm   = 0.5 * (c - b);
        if (std::fabs(xm) <= tol1 || fb == 0.0) return b;
        if (std::fabs(e) >= tol1 && std::fabs(fa) > std::fabs(fb)) {
            const double s = fb / fa;
            double p, q;
            if (a == c) { p = 2.0 * xm * s; q = 1.0 - s; }
            else {
                const double qa = fa / fc, qb = fb / fc;
                p = s * (2.0 * xm * qa * (qa - qb) - (b - a) * (qb - 1.0));
                q = (qa - 1.0) * (qb - 1.0) * (s - 1.0);
            }
            if (p > 0.0) q = -q;
            p = std::fabs(p);
            const double lim1 = 3.0 * xm * q - std::fabs(tol1 * q);
            const double lim2 = std::fabs(e * q);
            if (2.0 * p < (lim1 < lim2 ? lim1 : lim2)) { e = d; d = p / q; }
            else { d = xm; e = d; }
        } else { d = xm; e = d; }
        a = b; fa = fb;
        b += (std::fabs(d) > tol1) ? d
                                    : (xm >= 0.0 ?  std::fabs(tol1)
                                                 : -std::fabs(tol1));
        fb = f(b);
    }
    return b; // max iterations — return best estimate
}

} // anonymous namespace

ComparisonRecord makeComparison(const std::string &name, double reference,
                                double candidate, double absoluteTolerance) {
    ComparisonRecord record;
    record.name = name;
    record.reference = reference;
    record.candidate = candidate;
    record.absoluteError = std::fabs(reference - candidate);
    record.withinTolerance =
        record.absoluteError <= absoluteTolerance * (1.0 + std::fabs(reference));
    return record;
}

std::vector<ComparisonRecord> runReferenceComparison(double tolerance) {
    std::vector<ComparisonRecord> records;
    SimContext context;

    // -----------------------------------------------------------------------
    // R01 — RootFinder characterization: Brent (candidate) vs
    //       brentReference (independent NR re-implementation) vs analytic.
    // Eight mathematical functions with closed-form roots chosen to cover:
    //   monotone polynomial, oscillating trig, transcendental, nearly-flat
    //   near the root, and a stiff case.
    // -----------------------------------------------------------------------
    const double tightTol = 1e-10;
    const int maxIter = 200;
    struct R01Case {
        const char *name;
        double lo, hi, exact;
        ResidualFunction f;
    };
    R01Case r01Cases[] = {
        // f(x) = x^2 - 2, root = sqrt(2)
        {"R01/poly x^2-2", 0.0, 2.0, std::sqrt(2.0),
         [](double x) { return x * x - 2.0; }},
        // f(x) = x^3 - x - 2, root ≈ 1.5213797...
        {"R01/poly x^3-x-2", 1.0, 2.0, 1.5213797068045676,
         [](double x) { return x * x * x - x - 2.0; }},
        // f(x) = cos(x) - x, root ≈ 0.7390851...
        {"R01/trig cos-x", 0.0, 1.0, 0.7390851332151607,
         [](double x) { return std::cos(x) - x; }},
        // f(x) = exp(x) - 3, root = log(3) ≈ 1.0986123
        {"R01/transcend exp-3", 0.0, 2.0, std::log(3.0),
         [](double x) { return std::exp(x) - 3.0; }},
        // f(x) = sin(x) - 0.5, root = pi/6 ≈ 0.5235987...
        {"R01/trig sin-0.5", 0.0, 1.5, std::asin(0.5),
         [](double x) { return std::sin(x) - 0.5; }},
        // f(x) = x^5 - x - 1, root ≈ 1.1673039...
        {"R01/poly x^5-x-1", 1.0, 2.0, 1.1673039782614187,
         [](double x) {
             return x * x * x * x * x - x - 1.0;
         }},
        // f(x) = log(x) - 1, root = e
        {"R01/log(x)-1", 2.0, 4.0, std::exp(1.0),
         [](double x) { return std::log(x) - 1.0; }},
        // f(x) = x - cos(x^2), root ≈ 0.8010707652... (verified by scipy.brentq)
        {"R01/stiff x-cos(x^2)", 0.5, 1.5, 0.8010707652092184,
         [](double x) { return x - std::cos(x * x); }},
    };
    const std::size_t nCases =
        sizeof(r01Cases) / sizeof(r01Cases[0]);
    for (std::size_t k = 0; k < nCases; ++k) {
        const R01Case &tc = r01Cases[k];
        SolveResult cand =
            RootFinder::solveBracketed(tc.f, tc.lo, tc.hi, tightTol, maxIter);
        const double ref =
            brentReference(tc.f, tc.lo, tc.hi, tightTol, maxIter);
        std::string nameVsAnalytic = std::string(tc.name) + " vs analytic";
        std::string nameVsRef     = std::string(tc.name) + " vs ref-brent";
        records.push_back(makeComparison(nameVsAnalytic, tc.exact, cand.value,
                                         tolerance));
        records.push_back(makeComparison(nameVsRef, ref, cand.value, 1e-12));
    }

    // Analytic single-phase hydrostatic column.
    {
        FluidModel fluid(850.0, 1.2, 2.0e-3);
        ProductionColumn column = buildVerticalWell(fluid, 0.0, 20, 100.0, 0.15,
 0.0);
        const double separator = 10.0e5;
        SteadyStateRequest request(separator, 2.0, 1e-3, 200, 0);

        double lo = 0.0;
        double hi = 0.0;
        wellheadSearchBounds(column, request, lo, hi);
        ResidualFunction residual = makeWellheadResidual(column, context, separator);

        SolveResult brent =
            RootFinder::solveBracketed(residual, lo, hi, request.tolerance, 200);
        SolveResult bisection =
            RootFinder::solveBisection(residual, lo, hi, request.tolerance, 400);
        const double analytic = separator + fluid.liquidDensity() *
                                                constants::kGravity *
                                                column.totalVerticalHeight();

        records.push_back(makeComparison("hydrostatic: brent vs analytic",
                                         analytic, brent.value, tolerance));
        records.push_back(makeComparison("hydrostatic: bisection vs analytic",
                                         analytic, bisection.value, tolerance));
        records.push_back(makeComparison("hydrostatic: brent vs bisection",
                                         bisection.value, brent.value, tolerance));
    }

    // Two-phase column with flow: Brent against bisection.
    {
        FluidModel fluid(800.0, 1.2, 2.0e-3);
        ProductionColumn column =
            buildVerticalWell(fluid, 12.0, 20, 100.0, 0.15, 90.0);
        const double separator = 12.0e5;
        SteadyStateRequest request(separator, 3.0, 1e-4, 200, 0);

        double lo = 0.0;
        double hi = 0.0;
        wellheadSearchBounds(column, request, lo, hi);
        ResidualFunction residual = makeWellheadResidual(column, context, separator);

        SolveResult brent =
            RootFinder::solveBracketed(residual, lo, hi, request.tolerance, 200);
        SolveResult bisection =
            RootFinder::solveBisection(residual, lo, hi, request.tolerance, 400);

        records.push_back(makeComparison("two-phase: brent vs bisection",
                                         bisection.value, brent.value, tolerance));
    }

    return records;
}

bool comparisonPassed(const std::vector<ComparisonRecord> &records) {
    for (std::size_t i = 0; i < records.size(); ++i)
        if (!records[i].withinTolerance)
            return false;
    return true;
}

// ===========================================================================
// Test suite (implementation detail; entry points are declared in the header).
// ===========================================================================

namespace {

class TestReporter {
  public:
    explicit TestReporter(bool verbose)
        : verbose_(verbose), passed_(0), failed_(0) {}

    void check(const std::string &name, bool condition) {
        if (condition)
            ++passed_;
        else
            ++failed_;
        if (verbose_ || !condition)
            std::printf("  [%s] %s\n", condition ? "PASS" : "FAIL", name.c_str());
    }

    bool allPassed() const { return failed_ == 0; }
    int passed() const { return passed_; }
    int failed() const { return failed_; }

  private:
    bool verbose_;
    int passed_;
    int failed_;
};

} // namespace

bool runAllTests(bool verbose) {
    TestReporter reporter(verbose);
    SimContext context;
    NullDiagnostics silent;

    // Root finder ------------------------------------------------------------
    {
        ResidualFunction squareMinusTwo = [](double x) { return x * x - 2.0; };
        SolveResult brent =
            RootFinder::solveBracketed(squareMinusTwo, 0.0, 2.0, 1e-12, 100);
        SolveResult bisection =
            RootFinder::solveBisection(squareMinusTwo, 0.0, 2.0, 1e-12, 400);
        reporter.check("rootfinder/brent sqrt2", brent.ok());
        reporter.check("rootfinder/brent sqrt2 value",
                       almostEqual(brent.value, std::sqrt(2.0), 1e-9));
        reporter.check("rootfinder/bisection sqrt2 value",
                       almostEqual(bisection.value, std::sqrt(2.0), 1e-6));
        reporter.check("rootfinder/brent agrees with bisection",
                       almostEqual(brent.value, bisection.value, 1e-6));

        ResidualFunction alwaysPositive = [](double x) { return x * x + 1.0; };
        SolveResult nb =
            RootFinder::solveBracketed(alwaysPositive, 0.0, 2.0, 1e-9, 50);
        reporter.check("rootfinder/not-bracketed reported",
                       nb.status == SolveStatus::NotBracketed);
    }

    // R01 — RootFinder characterization against brentReference ---------------
    // Eight mathematical functions with closed-form roots (same as
    // runReferenceComparison). The check tolerance is 1e-6 (paridade
    // criterion). The vs-reference check is at 1e-12 (algorithmic identity).
    {
        const double tightTol = 1e-10;
        const int maxIter = 200;

        // Helper: wrap the test so we can call it with a label.
        struct R01Check {
            const char *label;
            double lo, hi, exact;
            ResidualFunction f;
        };
        R01Check cases[] = {
            {"r01/x^2-2",    0.0, 2.0, std::sqrt(2.0),
             [](double x){return x*x-2.0;}},
            {"r01/x^3-x-2",  1.0, 2.0, 1.5213797068045676,
             [](double x){return x*x*x-x-2.0;}},
            {"r01/cos-x",    0.0, 1.0, 0.7390851332151607,
             [](double x){return std::cos(x)-x;}},
            {"r01/exp-3",    0.0, 2.0, std::log(3.0),
             [](double x){return std::exp(x)-3.0;}},
            {"r01/sin-0.5",  0.0, 1.5, std::asin(0.5),
             [](double x){return std::sin(x)-0.5;}},
            {"r01/x^5-x-1",  1.0, 2.0, 1.1673039782614187,
             [](double x){double t=x*x; t*=t; return t*x-x-1.0;}},
            {"r01/log(x)-1", 2.0, 4.0, std::exp(1.0),
             [](double x){return std::log(x)-1.0;}},
            {"r01/stiff",    0.5, 1.5, 0.8010707652092184,
             [](double x){return x-std::cos(x*x);}},
        };
        const std::size_t nR01 = sizeof(cases)/sizeof(cases[0]);
        for (std::size_t k = 0; k < nR01; ++k) {
            const R01Check &tc = cases[k];
            SolveResult cand =
                RootFinder::solveBracketed(tc.f, tc.lo, tc.hi, tightTol, maxIter);
            const double ref =
                brentReference(tc.f, tc.lo, tc.hi, tightTol, maxIter);
            std::string vsAna = std::string(tc.label) + " vs analytic";
            std::string vsRef = std::string(tc.label) + " vs ref-brent";
            reporter.check(vsAna, almostEqual(cand.value, tc.exact, 1e-6));
            reporter.check(vsRef, almostEqual(cand.value, ref, 1e-12));
        }
    }


    // R03 — Two-phase drift-flux correlations (C0 / Ud) ----------------------
    // Reference values computed analytically (same equations as the legacy
    // SProd methods, verified by independent Python script).
    // Inputs: oil/water at 850/1.2 kg/m3, vertical 0.1 m pipe, Re=50000.
    {
        const double rhol     = 850.0;
        const double rhog     = 1.2;
        const double tensup   = 0.072;
        const double alf      = 0.30;
        const double dia      = 0.10;
        const double rug      = 0.0001;
        const double tet_vert = constants::kPi / 2.0;
        const double ug1      = 0.05;
        const double ul1      = 0.08;
        const double reymix   = 50000.0;
        const double reymixL  = 40000.0;
        const double cHor     = 1.0;
        double c0 = 0.0, ud = 0.0;

        choi(rhol, rhog, tensup, alf, reymix, reymixL, ug1, ul1,
             dia, rug, tet_vert, c0, ud, cHor);
        reporter.check("R03/choi c0",
                       almostEqual(c0, 1.1928421124, 1e-6));
        reporter.check("R03/choi ud",
                       almostEqual(ud, 0.2726455520, 1e-6));

        hibikiIshii(rhol, rhog, tensup, alf, reymix, reymixL, ug1, ul1,
                    dia, rug, tet_vert, c0, ud, cHor);
        reporter.check("R03/hibikiIshii c0",
                       almostEqual(c0, 2.5545405102, 1e-6));
        reporter.check("R03/hibikiIshii ud",
                       almostEqual(ud, 10.5160925813, 1e-6));

        francaLahey(rhol, rhog, tensup, alf, reymix, reymixL, ug1, ul1,
                    dia, rug, tet_vert, c0, ud, cHor);
        reporter.check("R03/francaLahey c0",
                       almostEqual(c0, 1.04, 1e-12));
        reporter.check("R03/francaLahey ud",
                       almostEqual(ud, 0.466, 1e-12));

        bhagwatGhajar(rhol, rhog, tensup, alf, reymix, reymixL, ug1, ul1,
                      dia, rug, tet_vert, c0, ud, cHor);
        reporter.check("R03/bhagwatGhajar c0",
                       almostEqual(c0, 1.1871843151, 1e-6));
        reporter.check("R03/bhagwatGhajar ud",
                       almostEqual(ud, 0.2460054528, 1e-6));

        // bhagwatGhajarMod uses reymixL in both Reynolds positions.
        double c0bg = 0.0, udbg = 0.0;
        bhagwatGhajar(rhol, rhog, tensup, alf, reymixL, reymixL, ug1, ul1,
                      dia, rug, tet_vert, c0bg, udbg, cHor);
        bhagwatGhajarMod(rhol, rhog, tensup, alf, reymix, reymixL, ug1, ul1,
                         dia, rug, tet_vert, c0, ud, cHor);
        reporter.check("R03/bhagwatGhajarMod c0 == BG@reymixL",
                       almostEqual(c0, c0bg, 1e-12));
        reporter.check("R03/bhagwatGhajarMod ud == BG@reymixL",
                       almostEqual(ud, udbg, 1e-12));

        // c0UdDisperso mode=0 (Choi) must equal choi directly.
        double c0d = 0.0, udd = 0.0;
        choi(rhol, rhog, tensup, alf, reymix, reymixL, ug1, ul1,
             dia, rug, tet_vert, c0d, udd, cHor);
        c0UdDisperso(rhol, rhog, tensup, alf, reymix, reymixL, ug1, ul1,
                     dia, rug, tet_vert, c0, ud, cHor, 0);
        reporter.check("R03/c0UdDisperso mode=0 c0 == choi",
                       almostEqual(c0, c0d, 1e-12));
        reporter.check("R03/c0UdDisperso mode=0 ud == choi",
                       almostEqual(ud, udd, 1e-12));

        // c0UdAnularChurn mode=3 (HibikiIshii) must equal hibikiIshii directly.
        double c0h = 0.0, udh = 0.0;
        hibikiIshii(rhol, rhog, tensup, alf, reymix, reymixL, ug1, ul1,
                    dia, rug, tet_vert, c0h, udh, cHor);
        c0UdAnularChurn(rhol, rhog, tensup, alf, reymix, reymixL, ug1, ul1,
                        dia, rug, tet_vert, c0, ud, cHor, 3);
        reporter.check("R03/c0UdAnularChurn mode=3 c0 == hibikiIshii",
                       almostEqual(c0, c0h, 1e-12));
        reporter.check("R03/c0UdAnularChurn mode=3 ud == hibikiIshii",
                       almostEqual(ud, udh, 1e-12));

        // c0UdEstratificado mode=2 (FrancaLahey) must equal francaLahey directly.
        double c0fl = 0.0, udfl = 0.0;
        francaLahey(rhol, rhog, tensup, alf, reymix, reymixL, ug1, ul1,
                    dia, rug, tet_vert, c0fl, udfl, cHor);
        c0UdEstratificado(rhol, rhog, tensup, alf, reymix, reymixL, ug1, ul1,
                          dia, rug, tet_vert, c0, ud, cHor, 2);
        reporter.check("R03/c0UdEstratificado mode=2 c0 == francaLahey",
                       almostEqual(c0, c0fl, 1e-12));
        reporter.check("R03/c0UdEstratificado mode=2 ud == francaLahey",
                       almostEqual(ud, udfl, 1e-12));

        // Physics sanity: c0 > 1 for vertical upward bubbly/churn flow.
        reporter.check("R03/choi c0 > 1",           c0d > 1.0);
        reporter.check("R03/hibikiIshii c0 > 1",     c0h > 1.0);
        reporter.check("R03/bhagwatGhajar c0 > 0.5", c0bg > 0.5);

        // Downward flow: ud must be negative for steep downward with vmix>0.01.
        const double tet_down = -constants::kPi / 2.0;
        choi(rhol, rhog, tensup, alf, reymix, reymixL, ug1, ul1,
             dia, rug, tet_down, c0, ud, cHor);
        reporter.check("R03/choi ud downward is negative", ud < 0.0);
    }

    // R05 — Pressure-gradient correlation dispatcher (build-time only) --------
    // computePressureGradient() requires GradientCorrelations.h which is only
    // available in the full CMake build.  These tests are compiled only when
    // MARLIM_BUILD is defined (i.e. by the CMake target, not the standalone
    // selftest harness).
#ifdef MARLIM_BUILD
    {
        PressureGradientInput in;
        // Standard vertical upward two-phase oil/gas flow
        in.inclinationRad   = constants::kPi / 2.0;  // vertical
        in.diameter         = 0.10;
        in.roughness        = 1e-4;
        in.pressure         = 5e6;   // 50 bar
        in.mixtureVelocity  = 1.5;   // m/s
        in.liquidFraction   = 0.60;
        in.gasDensity       = 30.0;  // kg/m3 at pressure
        in.liquidDensity    = 820.0; // kg/m3
        in.gasViscosity     = 1.5e-5;
        in.liquidViscosity  = 2.0e-3;
        in.surfaceTension   = 0.025;
        in.temperature      = 350.0; // K
        in.compressibilityFactor = 0.9;
        in.waterFraction    = 0.1;

        // Test each major correlation for physical sanity:
        //   holdup in [0,1], gravityGrad > 0 (upward), totalGrad > 0
        const FlowCorrelationId correlations[] = {
            FlowCorrelationId::PoettmannCarpenter,
            FlowCorrelationId::BaxendellThomas,
            FlowCorrelationId::HagedornBrown,
            FlowCorrelationId::DunsRos,
            FlowCorrelationId::Orkiszewski,
            FlowCorrelationId::BeggsAndBrill,
            FlowCorrelationId::MukherjeeeBrill,
            FlowCorrelationId::Aziz,
        };
        const char *names[] = {
            "R05/PoettmannCarpenter", "R05/BaxendellThomas",
            "R05/HagedornBrown",      "R05/DunsRos",
            "R05/Orkiszewski",         "R05/BeggsAndBrill",
            "R05/MukherjeeeBrill",    "R05/Aziz",
        };
        const int nCorrelations = sizeof(correlations) / sizeof(correlations[0]);
        for (int k = 0; k < nCorrelations; ++k) {
            PressureGradientResult r =
                computePressureGradient(correlations[k], in);
            std::string hLabel = std::string(names[k]) + "/holdup in [0,1]";
            std::string gLabel = std::string(names[k]) + "/gravity > 0";
            std::string tLabel = std::string(names[k]) + "/total > 0";
            reporter.check(hLabel, r.holdup >= 0.0 && r.holdup <= 1.0);
            reporter.check(gLabel, r.gravityGrad > 0.0);
            reporter.check(tLabel, r.totalGrad > 0.0);
        }

        // Determinism: same inputs produce identical outputs
        {
            PressureGradientResult r1 =
                computePressureGradient(FlowCorrelationId::BeggsAndBrill, in);
            PressureGradientResult r2 =
                computePressureGradient(FlowCorrelationId::BeggsAndBrill, in);
            reporter.check("R05/BeggsAndBrill deterministic",
                           r1.holdup == r2.holdup &&
                           r1.totalGrad == r2.totalGrad);
        }

        // Monotonicity: higher liquid fraction → higher gravity gradient
        {
            PressureGradientInput inHigh = in;
            inHigh.liquidFraction = 0.90;
            PressureGradientResult rLow =
                computePressureGradient(FlowCorrelationId::BeggsAndBrill, in);
            PressureGradientResult rHigh =
                computePressureGradient(FlowCorrelationId::BeggsAndBrill, inHigh);
            reporter.check("R05/BeggsAndBrill gravity increases with liquid fraction",
                           rHigh.gravityGrad >= rLow.gravityGrad);
        }
    }
#endif  // MARLIM_BUILD

    // R05 — marchToWellheadPhysical (build-time only, uses GradientCorrelations)
    // Tests that the physical march produces results consistent with the
    // simplified march: same monotonicity properties, pressure in range,
    // and determinism.
#ifdef MARLIM_BUILD
    {
        FluidModel fluid(800.0, 1.2, 2.0e-3);
        ProductionColumn column =
            buildVerticalWell(fluid, 12.0, 20, 100.0, 0.15, 90.0);
        const double bottomhole = 200.0e5;

        // Physical march must return a finite wellhead pressure.
        const double whPhysical =
            column.marchToWellheadPhysical(bottomhole, context,
                                           FlowCorrelationId::BeggsAndBrill);
        reporter.check("R05/marchPhysical finite",
                       std::isfinite(whPhysical));

        // Physical wellhead pressure must be less than bottomhole pressure
        // (positive pressure drop upward).
        reporter.check("R05/marchPhysical pwh < pbh",
                       whPhysical < bottomhole);

        // Determinism: two calls with the same inputs return the same value.
        const double whRepeat =
            column.marchToWellheadPhysical(bottomhole, context,
                                           FlowCorrelationId::BeggsAndBrill);
        reporter.check("R05/marchPhysical deterministic",
                       whPhysical == whRepeat);

        // Monotonicity: higher bottomhole pressure → higher wellhead pressure.
        const double whLow  =
            column.marchToWellheadPhysical(150.0e5, context,
                                           FlowCorrelationId::BeggsAndBrill);
        const double whHigh =
            column.marchToWellheadPhysical(200.0e5, context,
                                           FlowCorrelationId::BeggsAndBrill);
        reporter.check("R05/marchPhysical monotone with bottomhole",
                       whHigh > whLow);

        // Consistency with simplified march: both must agree on sign
        // (both positive or both < bottomhole).
        const double whSimple = column.marchToWellhead(bottomhole, context);
        reporter.check("R05/marchPhysical and simple agree sign",
                       (whPhysical < bottomhole) == (whSimple < bottomhole));

        // Multi-correlation smoke test: all available correlations must
        // produce finite results and positive gravity gradients.
        const FlowCorrelationId smokeCorrs[] = {
            FlowCorrelationId::PoettmannCarpenter,
            FlowCorrelationId::HagedornBrown,
            FlowCorrelationId::DunsRos,
            FlowCorrelationId::BeggsAndBrill,
            FlowCorrelationId::MukherjeeeBrill,
            FlowCorrelationId::Aziz,
        };
        for (int k = 0; k < 6; ++k) {
            const double wh =
                column.marchToWellheadPhysical(bottomhole, context, smokeCorrs[k]);
            reporter.check("R05/marchPhysical smoke finite",
                           std::isfinite(wh) && wh < bottomhole);
        }
    }
#endif  // MARLIM_BUILD (marchToWellheadPhysical)

    // R04 — Black-oil fluid-property correlations (build-time only).
    // Reference values: Python using the same Dranchuk-Abou-Kassem /
    // Lee-Kesler / Vasquez-Beggs equations as the legacy PropFlu methods.
    // Test conditions: natural gas Deng=0.75, PC=667 psia, TC=395 R.
#ifdef MARLIM_BUILD
    {
        const double Deng = 0.75;  // gas specific gravity
        const double PC   = 667.0; // psia
        const double TC   = 395.0; // Rankine

        // --- zFactor sanity checks ---
        // At low pressure (30 kgf/cm2 ≈ 3 MPa), Z ≈ 0.93 (weakly non-ideal)
        const double z30 = zFactor(30.0, 40.0, Deng, PC, TC);
        reporter.check("R04/zFactor physical (0.05,1.5)", z30 > 0.05 && z30 < 1.5);
        // Z must be positive and finite
        reporter.check("R04/zFactor positive finite",
                       std::isfinite(z30) && z30 > 0.0);
        // Z agrees with reference to 1e-6
        reporter.check("R04/zFactor at 30 kgf/cm2 40C",
                       almostEqual(z30, 0.92669683, 1e-6));

        const double z50 = zFactor(50.0, 60.0, Deng, PC, TC);
        reporter.check("R04/zFactor at 50 kgf/cm2 60C",
                       almostEqual(z50, 0.90111929, 1e-6));

        // --- gasDensityBlackOil sanity checks ---
        // Density must be positive and increase with pressure at same temperature
        const double rg30 = gasDensityBlackOil(30.0, 40.0, Deng, PC, TC);
        const double rg50 = gasDensityBlackOil(50.0, 60.0, Deng, PC, TC);
        reporter.check("R04/gasDensity positive",
                       rg30 > 0.0 && std::isfinite(rg30));
        reporter.check("R04/gasDensity at 30 kgf/cm2 40C",
                       almostEqual(rg30, 27.3679, 1e-4));
        reporter.check("R04/gasDensity at 50 kgf/cm2 60C",
                       almostEqual(rg50, 44.0919, 1e-4));
        // Higher pressure → higher density (isothermal)
        const double rg50_40 = gasDensityBlackOil(50.0, 40.0, Deng, PC, TC);
        reporter.check("R04/gasDensity increases with pressure",
                       rg50_40 > rg30);
        // Higher temperature → lower density (isobaric)
        const double rg50_80 = gasDensityBlackOil(50.0, 80.0, Deng, PC, TC);
        reporter.check("R04/gasDensity decreases with temperature",
                       rg50_80 < rg50_40);

        // --- gasViscosityBlackOil sanity checks ---
        const double mu30 = gasViscosityBlackOil(30.0, 40.0, Deng, PC, TC);
        reporter.check("R04/gasViscosity positive finite",
                       mu30 > 0.0 && std::isfinite(mu30));
        reporter.check("R04/gasViscosity at 30 kgf/cm2 40C",
                       almostEqual(mu30, 0.01166301, 1e-6));
        // Viscosity increases with pressure at same temperature (gas behavior)
        const double mu50 = gasViscosityBlackOil(50.0, 40.0, Deng, PC, TC);
        reporter.check("R04/gasViscosity increases with pressure",
                       mu50 > mu30);

        // --- solutionGOR sanity checks ---
        // Standard Standing coefficients for corrSat=0
        const double Avb = 1.0937, Bvb = 0.0008, Cvb = 25.724;
        const double API = 28.0;
        const double rs30 = solutionGOR(30.0, 60.0, API, Deng, Avb, Bvb, Cvb);
        const double rs50 = solutionGOR(50.0, 60.0, API, Deng, Avb, Bvb, Cvb);
        reporter.check("R04/solutionGOR non-negative",
                       rs30 >= 0.0 && std::isfinite(rs30));
        // GOR increases with pressure
        reporter.check("R04/solutionGOR increases with pressure", rs50 > rs30);

        // --- oilFVF sanity checks ---
        const double bo = oilFVF(30.0, 60.0, API, Deng, rs30);
        reporter.check("R04/oilFVF >= 1 (Bo > 1 below Pb)",
                       bo >= 1.0);
        reporter.check("R04/oilFVF finite", std::isfinite(bo));
        // Bo increases with GOR (more dissolved gas → larger volume)
        const double bo_high = oilFVF(30.0, 60.0, API, Deng, rs30 * 2.0);
        reporter.check("R04/oilFVF increases with GOR", bo_high > bo);
    }
#endif  // MARLIM_BUILD (R04)

    // R06/R10 — Colebrook friction factor + GLV aperture ----------------------
    // colebrookFrictionFactor: reference values computed by Python (same
    // algorithm as in the legacy BhagwatGhajar Colebrook loop).
    {
        // Laminar regime: f = 64/Re
        reporter.check("R06/colebrook laminar Re=500",
                       almostEqual(colebrookFrictionFactor(500.0, 0.0),
                                   64.0 / 500.0, 1e-9));
        // Transition boundary Re=2400 (still laminar by our convention)
        reporter.check("R06/colebrook laminar Re=2400",
                       almostEqual(colebrookFrictionFactor(2400.0, 0.001),
                                   64.0 / 2400.0, 1e-9));
        // Turbulent smooth: compare with Python/analytical (1e-6)
        reporter.check("R06/colebrook turbulent smooth Re=50000",
                       almostEqual(colebrookFrictionFactor(50000.0, 1e-4),
                                   0.021244373693, 1e-6));
        // Turbulent rough Re=1e6 eps=0.01
        reporter.check("R06/colebrook turbulent rough Re=1e6",
                       almostEqual(colebrookFrictionFactor(1e6, 0.01),
                                   0.037964741922, 1e-6));
        // Hydraulically smooth Re=1e5 eps=0
        reporter.check("R06/colebrook hydraulically smooth Re=1e5",
                       almostEqual(colebrookFrictionFactor(1e5, 0.0),
                                   0.017987525201, 1e-6));
        // Monotone in Re for fixed roughness: f decreases as Re increases
        const double f1 = colebrookFrictionFactor(5000.0, 0.001);
        const double f2 = colebrookFrictionFactor(50000.0, 0.001);
        const double f3 = colebrookFrictionFactor(500000.0, 0.001);
        reporter.check("R06/colebrook monotone decreasing with Re",
                       f1 > f2 && f2 > f3);
        // Monotone in roughness for fixed Re: f increases with eps
        const double fSmooth = colebrookFrictionFactor(1e5, 0.0);
        const double fRough  = colebrookFrictionFactor(1e5, 0.01);
        reporter.check("R06/colebrook monotone increasing with roughness",
                       fRough > fSmooth);
        // Consistency: matches existing darcyFrictionFactor at Re=1000 (laminar)
        reporter.check("R06/colebrook consistent with darcyFrictionFactor laminar",
                       almostEqual(colebrookFrictionFactor(1000.0, 0.0),
                                   darcyFrictionFactor(1000.0), 1e-6));
    }

    // R10/areaValvCali — GLV aperture ----------------------------------------
    // areaValvCali: mirrors SProd::areaValvCali verbatim (L2476 legacy).
    // Reference values verified by Python with the same formulas.
    {
        // Case 1: closed valve — compara < PBT → XMVS ≤ 0 → 0.0
        // PCal=1000, TCal=80, PVO=800, PT=300, dextern=0.04m,
        // areagarg=1e-4m2, Rvalv=0.1, Temp=120°F
        const double closed = areaValvCali(1000.0, 80.0, 800.0, 300.0,
                                           0.04, 1e-4, 0.1, 120.0);
        reporter.check("R10/areaValvCali closed valve is 0",
                       almostEqual(closed, 0.0, 1e-12));

        // Case 2: wide open — large PVO drives XMVS >> 0, APE clamped to areagarg
        // PCal=800, TCal=100, PVO=1200, PT=200, dextern=0.05m,
        // areagarg=5e-4m2, Rvalv=0.2, Temp=150°F → aperture=1.0
        const double fullyOpen = areaValvCali(800.0, 100.0, 1200.0, 200.0,
                                              0.05, 5e-4, 0.2, 150.0);
        reporter.check("R10/areaValvCali fully open is 1",
                       almostEqual(fullyOpen, 1.0, 1e-6));

        // Case 3: result in [0, 1] for any valid input
        reporter.check("R10/areaValvCali result in [0,1] closed", closed >= 0.0 && closed <= 1.0);
        reporter.check("R10/areaValvCali result in [0,1] open", fullyOpen >= 0.0 && fullyOpen <= 1.0);

        // Case 4: determinism — same inputs produce exactly same output
        const double repeat = areaValvCali(1000.0, 80.0, 800.0, 300.0,
                                           0.04, 1e-4, 0.1, 120.0);
        reporter.check("R10/areaValvCali deterministic", closed == repeat);
    }

    // Stream mixing ----------------------------------------------------------
    {
        StandardStream upstream;
        upstream.oilRate = 100.0;
        upstream.waterRate = 25.0;
        upstream.gasRate = 9000.0;
        upstream.oilApi = 28.0;

        StandardStream empty;
        StandardStream unchanged = StreamMixing::march(
            upstream, AccessoryType::LiquidSource, empty, context);
        reporter.check("mixing/identity api",
                       almostEqual(unchanged.oilApi, upstream.oilApi, 1e-9));
        reporter.check("mixing/identity bsw",
                       almostEqual(unchanged.basicSedimentAndWater(),
                                   upstream.basicSedimentAndWater(), 1e-9));
        reporter.check("mixing/identity gor",
                       almostEqual(unchanged.gasOilRatio(),
                                   upstream.gasOilRatio(), 1e-9));

        StandardStream a;
        a.oilRate = 80.0;
        a.waterRate = 20.0;
        a.gasRate = 6000.0;
        a.oilApi = 30.0;
        StandardStream b;
        b.oilRate = 40.0;
        b.waterRate = 40.0;
        b.gasRate = 2000.0;
        b.oilApi = 22.0;
        StandardStream mixed = blend(a, b);
        reporter.check("mixing/oil conserved",
                       almostEqual(mixed.oilRate, 120.0, 1e-12));
        reporter.check("mixing/water conserved",
                       almostEqual(mixed.waterRate, 60.0, 1e-12));
        reporter.check("mixing/gas conserved",
                       almostEqual(mixed.gasRate, 8000.0, 1e-12));
        reporter.check("mixing/api between inputs",
                       mixed.oilApi > 22.0 && mixed.oilApi < 30.0);
    }

    // Trend buffer -----------------------------------------------------------
    {
        TrendBuffer buffer(2, 4);
        buffer.append(0, 1.0);
        buffer.append(0, 2.0);
        buffer.append(1, 3.0);
        reporter.check("trend/channel0 size", buffer.size(0) == 2);
        reporter.check("trend/channel1 size", buffer.size(1) == 1);
        buffer.reset();
        reporter.check("trend/reset clears", buffer.size(0) == 0);
    }

    // Steady state: exact hydrostatic reference ------------------------------
    {
        FluidModel fluid(850.0, 1.2, 2.0e-3);
        ProductionColumn column = buildVerticalWell(fluid, 0.0, 20, 100.0, 0.15,
                                                    0.0);
        const double separatorPressure = 10.0e5;
        SteadyStateRequest request(separatorPressure, 2.0, 1e-3, 200, 0);
        TramoEngine engine(context, silent);
        SolveResult solved = engine.solveBottomholePressure(column, request);
        const double expected = separatorPressure + fluid.liquidDensity() *
                                                        constants::kGravity *
                                                        column.totalVerticalHeight();
        reporter.check("steady/hydrostatic converges", solved.ok());
        reporter.check("steady/hydrostatic exact",
                       almostEqual(solved.value, expected, 1e-6));
    }

    // Steady state: determinism and round trip -------------------------------
    {
        FluidModel fluid(800.0, 1.2, 2.0e-3);
        ProductionColumn column =
            buildVerticalWell(fluid, 12.0, 20, 100.0, 0.15, 90.0);
        const double separatorPressure = 12.0e5;
        SteadyStateRequest request(separatorPressure, 3.0, 1e-2, 200, 0);
        TramoEngine engine(context, silent);
        SolveResult first = engine.solveBottomholePressure(column, request);
        SolveResult second = engine.solveBottomholePressure(column, request);
        reporter.check("steady/twophase converges", first.ok());
        reporter.check("steady/deterministic", first.value == second.value);
        const double wellhead = column.marchToWellhead(first.value, context);
        reporter.check("steady/round trip matches separator",
                       almostEqual(wellhead, separatorPressure, 1e-2));
    }

    // Physics sanity: monotonicity and gas lift ------------------------------
    {
        FluidModel fluid(800.0, 1.2, 2.0e-3);
        ProductionColumn column =
            buildVerticalWell(fluid, 12.0, 20, 100.0, 0.15, 90.0);
        const double low  = column.marchToWellhead(150.0e5, context);
        const double high = column.marchToWellhead(180.0e5, context);
        reporter.check("physics/wellhead increases with bottomhole", high > low);

        ProductionColumn withoutLift =
            buildVerticalWell(fluid, 12.0, 20, 100.0, 0.15, 90.0);
        ProductionColumn withLift(fluid, 12.0);
        withLift.reserveSegments(20);
        withLift.setInletStream(withoutLift.inletStream());
        for (std::size_t i = 0; i < 20; ++i) {
            PipeSegment segment;
            segment.length = 100.0;
            segment.diameter = 0.15;
            segment.inclination = constants::kPi / 2.0;
            if (i == 10) {
                segment.accessory = AccessoryType::GasSource;
                segment.sourceStream.gasRate = 600.0;
            }
            withLift.addSegment(segment);
        }
        const double bottomhole = 170.0e5;
        const double plain  = withoutLift.marchToWellhead(bottomhole, context);
        const double lifted = withLift.marchToWellhead(bottomhole, context);
        reporter.check("gaslift/raises wellhead pressure", lifted > plain);
    }

    // Batch solver: parallel results equal the sequential ones ---------------
    {
        FluidModel fluid(800.0, 1.2, 2.0e-3);
        std::vector<ProductionColumn> columns;
        std::vector<SteadyStateRequest> requests;
        columns.reserve(64);
        requests.reserve(64);
        for (int i = 0; i < 64; ++i) {
            columns.push_back(buildVerticalWell(fluid, 8.0 + 0.1 * i, 20, 100.0,
                                                0.15, 60.0 + i));
            requests.push_back(SteadyStateRequest(12.0e5, 3.0, 1e-2, 200, i));
        }

        SimContext sequentialContext;
        sequentialContext.setThreadCount(1);
        std::vector<SolveResult> sequential =
            solveBatch(columns, requests, sequentialContext);

        SimContext parallelContext;
        parallelContext.setThreadCount(4);
        std::vector<SolveResult> parallel =
            solveBatch(columns, requests, parallelContext);

        bool allMatch = sequential.size() == parallel.size();
        for (std::size_t i = 0; i < sequential.size() && allMatch; ++i)
            allMatch = sequential[i].ok() && parallel[i].ok() &&
                       sequential[i].value == parallel[i].value;
        reporter.check("batch/parallel equals sequential", allMatch);
    }

    // Automatic comparison against the independent reference -----------------
    {
        std::vector<ComparisonRecord> records = runReferenceComparison(1e-6);
        reporter.check("comparison/new solver matches reference",
                       comparisonPassed(records));
    }

    if (verbose) {
        std::printf("  ----------------------------------------\n");
        std::printf("  total: %d passed, %d failed\n", reporter.passed(),
                    reporter.failed());
    }
    return reporter.allPassed();
}

bool runSelfTest() { return runAllTests(false); }

} // namespace sisprod2
} // namespace marlim
