// Drive the three root-finding solvers over synthetic objectives and record
// every evaluation, in order.
//
// Two things this covers that nothing else does.
//
// First, instantiation. zbrent has no call site anywhere in the project and
// falsacorda is reachable only from inside it, so after the extraction they are
// templates that nothing instantiates -- and an uninstantiated template is only
// parsed, never type-checked. Without this driver the compiler would not look
// at 96 of the 240 moved lines.
//
// Second, trajectory. The contract is explicit that a root finder comparing
// only its returned value proves nothing: two different convergence paths can
// land on the same number. So the objective logs every argument it is handed
// and every value it returns, in order, and the comparison is over that whole
// sequence.
//
// Values are written in C99 hexadecimal float (%a), which round-trips exactly.
//
// Built and run by verify-solvers.sh. Not part of the product build: it lives
// in refactor-harness/, not src/core/, so file(GLOB src/core/*.cpp) never sees
// it.

#include "RootFindingSolvers.h"

#include <cstdio>
#include <cstdlib>

namespace {

std::FILE *out = nullptr;

// The shapes an objective can take. Chosen for the paths they force rather than
// for physical meaning: the solvers branch on sign, on magnitude above 1e9, and
// on whether the interval brackets a root at all.
enum Shape {
    kLinear,        // one crossing, well behaved
    kCubic,         // one crossing, flat near the root -- stresses interpolation
    kQuadratic,     // two crossings, so some brackets contain none
    kAlwaysPositive,// no crossing: forces the widening loops and zbrent's fallback
    kAlwaysNegative,// no crossing, other side
    kHuge,          // above 1e9, forces the guard that returns 1e10
    kFlatZero,      // identically zero: forces fb == 0.0
    kTiny,          // crosses zero at 1e-200 scale, so that fm*fm and fl*fh both
                    // underflow to zero and zriddr's `s == 0.0` branch -- the one
                    // A2-05 lives in -- is actually reached. Without it the sweep
                    // never enters that branch, which the calibration found.
    kShapeCount
};

const char *nameOf(Shape shape) {
    switch (shape) {
    case kLinear:         return "linear";
    case kCubic:          return "cubic";
    case kQuadratic:      return "quadratic";
    case kAlwaysPositive: return "always-positive";
    case kAlwaysNegative: return "always-negative";
    case kHuge:           return "huge";
    case kFlatZero:       return "flat-zero";
    case kTiny:           return "tiny";
    case kShapeCount:     break;
    }
    return "?";
}

double shapeValue(Shape shape, double x) {
    switch (shape) {
    case kLinear:         return 2.5 * x - 30.0;
    case kCubic:          return (x - 12.0) * (x - 12.0) * (x - 12.0);
    case kQuadratic:      return (x - 5.0) * (x - 40.0);
    case kAlwaysPositive: return x * x + 1.0;
    case kAlwaysNegative: return -(x * x) - 1.0;
    case kHuge:           return 5e9 * x + 1.0;
    case kFlatZero:       return 0.0;
    case kTiny:           return 1e-200 * (x - 12.0);
    case kShapeCount:     break;
    }
    return 0.0;
}

// The objective handed to the solvers. Logging happens here, so the record is
// of what the solver asked for, not of what the driver assumed it would ask.
struct Objective {
    Shape shape;
    mutable long long calls = 0;

    double operator()(double x) const {
        double value = shapeValue(shape, x);
        ++calls;
        std::fprintf(out, "  eval %lld %a %a\n", calls, x, value);
        return value;
    }
};

// Mirrors the shape of the lambda SProd::zriddr builds: divide by the base,
// then record when the domain says to. `records` stands in for `prod != 0`.
struct Monitor {
    double base;
    int records;
    mutable double lastRecorded = 0.0;

    double operator()(double residual) const {
        double normalized = residual / base;
        if (records != 0) {
            lastRecorded = fabs(normalized);
        }
        return normalized;
    }
};

struct Bracket { double low, high; };

const Bracket kBrackets[] = {
    {  1.0,  30.0},   // brackets the linear and cubic roots
    { 30.0,   1.0},   // reversed, so x2 < x1
    {  6.0,  35.0},   // inside the quadratic's negative region: no crossing
    {  0.0,  50.0},   // spans both quadratic roots
    { 11.999, 12.001},// tight around the cubic root
    { -5.0,  -1.0},   // entirely left of every root
};

const double kMonitorBases[] = {1.0, 2.5};
const int kReverseMarch[] = {0, 1};
const int kMinIterations[] = {0, 10};

}  // namespace

int main(int argc, char **argv) {
    if (argc != 2) {
        std::fprintf(stderr, "usage: solver-sweep <output-file>\n");
        return 2;
    }
    out = std::fopen(argv[1], "w");
    if (!out) {
        std::fprintf(stderr, "cannot write %s\n", argv[1]);
        return 2;
    }

    long long configurations = 0;

    for (int s = 0; s < kShapeCount; ++s) {
        Shape shape = static_cast<Shape>(s);
        for (const Bracket &bracket : kBrackets) {

            // ---- zriddr: the live solver, swept over its two extra selectors.
            for (int reverseMarch : kReverseMarch)
            for (int minIterations : kMinIterations)
            for (double base : kMonitorBases)
            for (int records : kReverseMarch) {
                Objective objective{shape};
                Monitor monitor{base, records};
                std::fprintf(out, "zriddr %s %a %a rev=%d min=%d base=%a rec=%d\n",
                             nameOf(shape), bracket.low, bracket.high,
                             reverseMarch, minIterations, base, records);
                double root = rootfinding::zriddr(bracket.low, bracket.high,
                                                  objective, monitor,
                                                  reverseMarch, minIterations);
                std::fprintf(out, "  root %a evals %lld monitor %a\n",
                             root, objective.calls, monitor.lastRecorded);
                ++configurations;
            }

            // ---- zbrent, and through it falsacorda. Never reached in the
            // product; this is the only place either one runs at all.
            Objective brent{shape};
            std::fprintf(out, "zbrent %s %a %a\n",
                         nameOf(shape), bracket.low, bracket.high);
            double brentRoot = rootfinding::zbrent(bracket.low, bracket.high, brent,
                                                   0.00001, 0.00001, 100);
            std::fprintf(out, "  root %a evals %lld\n", brentRoot, brent.calls);
            ++configurations;

            // ---- zbrent again with a budget of two iterations, which is the
            // only way to reach the exhaustion path and, with it,
            // reportIterationLimit. Nothing else in the sweep converges slowly
            // enough to run out of a hundred.
            Objective starved{shape};
            std::fprintf(out, "zbrent-starved %s %a %a\n",
                         nameOf(shape), bracket.low, bracket.high);
            double starvedRoot = rootfinding::zbrent(bracket.low, bracket.high, starved,
                                                     0.00001, 0.00001, 2);
            std::fprintf(out, "  root %a evals %lld\n", starvedRoot, starved.calls);
            ++configurations;

            // ---- falsacorda on its own, so a defect in it is attributable to
            // it rather than showing up only through zbrent's fallback.
            Objective chord{shape};
            std::fprintf(out, "falsacorda %s %a %a\n",
                         nameOf(shape), bracket.low, bracket.high);
            double chordRoot = rootfinding::falsacorda(bracket.low, bracket.high, chord);
            std::fprintf(out, "  root %a evals %lld\n", chordRoot, chord.calls);
            ++configurations;
        }
    }

    // The two sign helpers, over the values where they branch. sign() has no
    // caller in the product either.
    const double kSignInputs[] = {-2.5, -0.0, 0.0, 1e-300, 2.5};
    for (double a : kSignInputs)
        for (double b : kSignInputs)
            std::fprintf(out, "SIGN %a %a %a\n", a, b, rootfinding::SIGN(a, b));
    for (double v : kSignInputs)
        std::fprintf(out, "sign %a %d\n", v, rootfinding::sign(v));

    std::fprintf(out, "configurations %lld\n", configurations);
    std::fclose(out);
    std::printf("%lld configuration(s) written to %s\n", configurations, argv[1]);
    return 0;
}
