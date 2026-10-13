#ifndef RUNTIMEDIAGNOSTICS_H_
#define RUNTIMEDIAGNOSTICS_H_

/// Warnings the solution gathers while it runs, for conditions that can recur at
/// every cell of every time step. One log entry per occurrence would flood the
/// log, so each kind keeps a count, its first occurrence and the range it
/// covered, and becomes one entry when the log is written.
///
/// Recording is thread-safe, and a caller pays for it only when the condition
/// occurs: it tests the condition first. The entries reach the logger only
/// through flushTo, which runs on the thread that writes the log, so the shared
/// logger is never called from inside a parallel region.

#include <string>

class Logger;

namespace diagnostics {

/// The compositional flash routines whose general error code, iIER, is checked.
enum class FlashRoutine {
    /// Marlim_CalculateCommonlyRequiredValuesAtStockTankConditions, for API, gas
    /// density and gas-oil ratio at the stock tank.
    StockTank,
    /// Marlim_CalculateMixtureThermodynamicCondition, at a given pressure and
    /// temperature.
    Mixture,
};

/// A flash returned a nonzero general error code, and the caller used its results
/// as returned. Pressure in kgf/cm2, temperature in C, simulation time in s.
void recordFlashError(FlashRoutine routine, int errorCode, double pressure, double temperature, double time);

/// What the caller of the stock-tank flash puts in place of a result it cannot use.
enum class StockTankDefault {
    /// No stock-tank liquid density: the API gravity is set to 50.
    ApiGravity,
    /// Vapor without a gas-oil ratio: the gas-oil ratio is set to 1e6.
    GasOilRatio,
};

/// The stock-tank flash's caller applied `kind` at simulation time `time`.
void recordStockTankDefault(StockTankDefault kind, double time);

/// A mixture flash failed with a positive code and, the mixture being light, the
/// caller replaced its results with vapor of the feed composition.
void recordLightMixtureFallback(int errorCode, double pressure, double temperature, double mixtureMolarMass,
                                double time);

/// Gas the closing-valve time-step control adds to the source of the valve's
/// cell, when a closing valve holds almost no gas. The time step that uses it
/// may still be shortened, so a step keeps its additions until it is accepted.
struct ArtificialGasAddition {
    /// 0 for the Master1 valve, k for the production valve k.
    int valve;
    int cell;
    /// kg/s.
    double massRate;
};

/// An accepted time step, `step` of `system`'s run, that started at `time` and
/// lasted `timeStep` seconds with `addition` in its sources.
void recordArtificialGas(const std::string &system, const ArtificialGasAddition &addition, int step, double time,
                         double timeStep);

/// Moves every gathered warning into `logger` and forgets it, so a second call
/// adds nothing.
void flushTo(Logger &logger);

} // namespace diagnostics

#endif
