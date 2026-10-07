#ifndef SISPRODCONSTANTS_H_
#define SISPRODCONSTANTS_H_

/// Physical constants, accessory kinds and gas-line inlet conditions shared by
/// SProd and its modules. The static_asserts at the end pin each value to its
/// spelling.
namespace sisprod {

// ---------------------------------------------------------------- units ----

/// Used as the ratio 6.29 / 35.31467, converting a solution gas-oil ratio from
/// scf/bbl to m3/m3. Leave it as a division: the quotient written by hand
/// rounds to a different double.
inline constexpr double kBarrelPerCubicMetre = 6.29;
inline constexpr double kCubicFootPerCubicMetre = 35.31467;

/// Density of air at standard conditions, in kg/m3. Always multiplied by a gas
/// specific gravity, giving the gas density at standard conditions.
inline constexpr double kAirDensityAtStandardConditions = 1.225;

/// Gravitational acceleration as this program uses it. Not 9.80665 --
/// "correcting" it would change every hydrostatic term.
inline constexpr double kGravity = 9.82;

/// The gas-lift unloading hydrostatics spell gravity differently. As with the
/// pressure variant above, the two are not interchangeable without changing
/// results.
inline constexpr double kGravityUnloadingVariant = 9.81;

/// Written into the inlet source rate when a pressure-to-pressure search finds no
/// flowing solution; the network reads it back to report the branch.
inline constexpr double kNoFlowSolutionMarker = -2121212121.;

/// Standard conditions as this program spells them, for gas density at
/// surface: one kgf/cm^2 and fifteen degrees Celsius.
inline constexpr double kStandardPressureKgfPerCm2 = 1.;
inline constexpr double kStandardTemperatureCelsius = 15.;

// ------------------------------------------------------------- numerics ----

/// Backward perturbation for numerical derivatives.
inline constexpr double kDerivativePerturbationFactor = 0.999;

/// Below this a phase-change mass rate is treated as zero.
inline constexpr double kPhaseChangeFloor = 1e-25;

/// The temperature range the solution is clamped to, in degrees Celsius.
inline constexpr double kMinimumTemperatureCelsius = -50.;
inline constexpr double kMaximumTemperatureCelsius = 200.;

/// The lowest cell pressure, in kgf/cm^2, a production steady march goes on
/// with; at or below it the march stops as a pressure too low.
inline constexpr double kMarchMinimumPressure = 0.1;

/// True when a cell's pressure stops a production steady march: at or below
/// kMarchMinimumPressure, or NaN, which fails every comparison.
inline constexpr bool marchPressureTooLow(double pressure) {
    return !(pressure > kMarchMinimumPressure);
}

// ------------------------------------------------------ accessory kinds ----

/// The value of `celula[i].acsr.tipo`, selecting the accessory attached to a
/// cell. Read from where the field is assigned, in Leitura.cpp and
/// LeituraVapor.cpp; the contradictory table in acessorios.h belongs to a dead
/// function. Type 7 is absent because no assignment names it.
enum AccessoryKind : int {
    kAccessoryNone = 0,             ///< no source attached
    kAccessoryGasInjection = 1,     ///< acsr.injg
    kAccessoryLiquidInjection = 2,  ///< acsr.injl
    kAccessoryInflowPerformance = 3,///< acsr.ipr
    kAccessoryPump = 4,             ///< acsr.bcs, electrical submersible pump
    kAccessoryChoke = 5,            ///< acsr.chk
    kAccessoryVolumetricPump = 8,   ///< acsr.bvol
    kAccessoryLeak = 9,             ///< acsr.fontechk
    kAccessoryMultipleSource = 10,  ///< acsr.injm
    kAccessoryRadialPorous = 15,    ///< acsr.radialPoro
    kAccessoryPorous2D = 16,        ///< acsr.poroso2D
    kAccessoryMultiPump = 17,       ///< acsr.multibcs
};

// ------------------------------------------------ gas-line inlet condition ----

/// The value of `tipoCC` on the gas line's first cell, and of the configured
/// `gasinj.tipoCC`: what drives the gas line at the injection point. From the
/// input schema (docs/schemas/branch.pt.json, gasInj.tipoCC): 0 is the
/// injection pressure, 1 the injection flow rate.
enum : int {
    kGasInletInjectionPressure = 0, ///< the injection pressure is given
    kGasInletInjectionFlowRate = 1, ///< the injection flow rate is given
};

/// The two conditions as types. Where the gas line behaves differently under
/// each, the difference is written once per condition, as an overload taking
/// one of these, in the module that owns the data, instead of the flag being
/// tested at every place that cares.
struct InjectionPressureCondition {};
struct InjectionFlowRateCondition {};

/// The one place the inlet condition is decided: calls `visitor` with the
/// condition `tipoCC` names. Zero is the injection pressure and anything else
/// the flow rate; the input does not reject a value outside {0, 1}. The one site
/// that reads the flag the other way round says so where it is.
///
/// Called where the condition matters, not once in an adapter: the first gas
/// cell exists only when there is a gas line, and every call sits behind a
/// guard that says so.
template <typename Visitor>
decltype(auto) withGasInletCondition(int tipoCC, Visitor &&visitor) {
    if (tipoCC == kGasInletInjectionPressure)
        return visitor(InjectionPressureCondition{});
    return visitor(InjectionFlowRateCondition{});
}

static_assert(kAccessoryNone == 0);
static_assert(kAccessoryGasInjection == 1);
static_assert(kAccessoryLiquidInjection == 2);
static_assert(kAccessoryInflowPerformance == 3);
static_assert(kAccessoryPump == 4);
static_assert(kAccessoryChoke == 5);
static_assert(kAccessoryVolumetricPump == 8);
static_assert(kAccessoryLeak == 9);
static_assert(kAccessoryMultipleSource == 10);
static_assert(kAccessoryRadialPorous == 15);
static_assert(kAccessoryPorous2D == 16);
static_assert(kAccessoryMultiPump == 17);
static_assert(kGasInletInjectionPressure == 0);
static_assert(kGasInletInjectionFlowRate == 1);

static_assert(kBarrelPerCubicMetre == 6.29);
static_assert(kCubicFootPerCubicMetre == 35.31467);
static_assert(kAirDensityAtStandardConditions == 1.225);
static_assert(kGravity == 9.82);
static_assert(kGravityUnloadingVariant == 9.81);
static_assert(kNoFlowSolutionMarker == -2121212121.);
static_assert(kStandardPressureKgfPerCm2 == 1.);
static_assert(kStandardTemperatureCelsius == 15.);
static_assert(kDerivativePerturbationFactor == 0.999);
static_assert(kPhaseChangeFloor == 1e-25);
static_assert(kMinimumTemperatureCelsius == -50.);
static_assert(kMaximumTemperatureCelsius == 200.);
static_assert(kMarchMinimumPressure == 0.1);
static_assert(marchPressureTooLow(0.1) && !marchPressureTooLow(0.2));

}  // namespace sisprod

#endif  // SISPRODCONSTANTS_H_
