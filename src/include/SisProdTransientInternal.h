#ifndef SISPRODTRANSIENTINTERNAL_H_
#define SISPRODTRANSIENTINTERNAL_H_

#include "SisProdTransient.h"
#include "Leitura.h"

#include <concepts>
#include <ctime>
#include <iosfwd>

namespace sisprod::transient {
namespace {

/// Ordered application points of the single-line time-step policies.
enum class TimeStepHook {
    CouplingIterationStart, // start of a pressure-volume coupling iteration
    AfterPigUpdate,         // after the pig moves, in the same iteration
    AfterValveOpenings,     // after the step re-reads the valve openings
};

/// Live step bindings, coupling iteration, and explicit-scheme selector.
struct TimeStepPolicyContext {
    const TransientStepState &step;
    int couplingIteration;
    int explicitScheme;
};

constexpr int kOutsideCouplingLoop = -1;

/// Compile-time interface required of each time-step policy.
template <typename Policy>
concept TimeStepPolicy = requires(const TimeStepPolicyContext &context) {
    { Policy::hook } -> std::convertible_to<TimeStepHook>;
    { Policy::applies(context) } -> std::same_as<bool>;
    { Policy::apply(context) } -> std::same_as<void>;
};

/// Valve restriction at the first coupling iteration when enabled.
struct RestrictTimeStepByValvePolicy {
    static constexpr TimeStepHook hook = TimeStepHook::CouplingIterationStart;
    /// Tests whether this policy acts at the current step stage.
    static bool applies(const TimeStepPolicyContext &context) {
        return context.couplingIteration == 0 && context.step.input.controleDTvalv == 1;
    }
    /// Executes the policy action using the live step bindings.
    static void apply(const TimeStepPolicyContext &context) {
        restrictTimeStepByValve(context.step);
    }
};

/// Maximum-step damping at the coupling iteration selected by the model.
struct DampMaximumTimeStepPolicy {
    static constexpr TimeStepHook hook = TimeStepHook::AfterPigUpdate;
    /// Tests whether this policy acts at the current step stage.
    static bool applies(const TimeStepPolicyContext &context) {
        return context.couplingIteration == 1 * context.step.fullModel;
    }
    /// Executes the policy action using the live step bindings.
    static void apply(const TimeStepPolicyContext &context) {
        dampMaximumTimeStep(context.step);
    }
};

/// Pressure and temperature rate assessment after valve openings are read.
struct PressureRateOfChangePolicy {
    static constexpr TimeStepHook hook = TimeStepHook::AfterValveOpenings;
    /// Tests whether this policy acts at the current step stage.
    static bool applies(const TimeStepPolicyContext &context) {
        return context.step.fullModel == 1;
    }
    /// Executes the policy action using the live step bindings.
    static void apply(const TimeStepPolicyContext &context) {
        evaluatePressureRateOfChange(context.step, context.explicitScheme);
    }
};

/// Applies policies of the requested hook in their declared order.
template <TimeStepPolicy... Policies>
struct TimeStepPolicyRegistry {
    /// Executes matching policies without indirect dispatch.
    template <TimeStepHook Hook>
    static void apply(const TimeStepPolicyContext &context) {
        (applyIfAt<Policies, Hook>(context), ...);
    }

  private:
    /// Tests and executes a policy only at its own hook.
    template <typename Policy, TimeStepHook Hook>
    static void applyIfAt(const TimeStepPolicyContext &context) {
        if constexpr (Policy::hook == Hook) {
            if (Policy::applies(context))
                Policy::apply(context);
        }
    }
};

/// Policies applied by the single-line transient step.
using TimeStepPolicies = TimeStepPolicyRegistry<
    RestrictTimeStepByValvePolicy,
    DampMaximumTimeStepPolicy,
    PressureRateOfChangePolicy>;

/// Returns the absolute gas-to-total mass-flow ratio at the outlet.
double flowQuality(double gasMassFlow, double totalMassFlow);

/// Computes the in-situ gas mass fraction from void fraction and mixture density.
double inSituGasMassFraction(const Cel &cell, double voidFraction, double pressure, double mixtureDensity);

/// Records accepted artificial gas with the final step length and clears pending additions.
void reportAcceptedArtificialGas(const TransientStepState &state);

/// Advances one ordered fraction, porous-flow, pressure-volume, and thermal coupling iteration.
void advanceCouplingIteration(const TransientSolveState &state, int kontaAcop, int celpos, int vExpli, int ciclomax, double titRev, double alfRev, double betRev);

/// Emits scheduled line, wall, and unit-cell profiles and advances their schedule indices.
void writeProfiles(const TransientSolveState &state, int nrede);

/// Samples and flushes the four trend channels using their individual initial-time rules.
void writeTrends(const TransientSolveState &state, int ordemImpT, double velmaxdesc, int nrede);

/// Writes selected screen variables and the elapsed duration of the step calculation.
void writeScreenOutput(const TransientSolveState &state, const chrono::steady_clock::time_point &begin, const chrono::steady_clock::time_point &end);

/// Appends due events using the input event sentinel and advances the event index.
void writeEventLog(const TransientSolveState &state, int maxEvento);

/// Writes progress, monitored values, and final-run information when a report is due.
void writeProgressReport(const TransientSolveState &state, int MaxKontaImpres);

/// Fills the minimum-pressure corners in temperature order, including bubble pressures.
void fillMiniTableCornersAtMinPressure(const TransientStepState &state, ProFlu &fluC, ProFlu &flui);

/// Fills the maximum-pressure corners in temperature order.
void fillMiniTableCornersAtMaxPressure(const TransientStepState &state, ProFlu &fluC, ProFlu &flui);

/// Propagates the current time step to the four cell increments in index order.
void assignCellTimeSteps(const TransientStepState &state);

/// Tests the scheduled instant with short-circuit reads of the time array and index.
bool profileIsDue(const TransientStepState &state, double *const &times, const int &index);

/// Advances a profile index, retaining the final schedule entry.
void advanceProfileIndex(int &index, const int &limit);

/// Appends the branch-qualified log path to the caller-owned stream.
void appendLogPath(std::ostringstream &saidaT, const TransientSolveState &state);

/// Writes the supplied local timestamp and its terminating flush.
void writeLogTimestamp(std::ostream &escreveIni, const std::tm *ltm);

/// Appends elapsed-run and closing information when the final time is reached.
void writeRunCompletionReport(const TransientSolveState &state, std::ostream &escreveIni);

/// Tests progress, report-count, and final-time triggers at the current clock value.
bool progressReportIsDue(const TransientSolveState &state, int MaxKontaImpres);

}  // namespace
}  // namespace sisprod::transient

#endif  // SISPRODTRANSIENTINTERNAL_H_
