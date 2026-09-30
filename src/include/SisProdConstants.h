#ifndef SISPRODCONSTANTS_H_
#define SISPRODCONSTANTS_H_

/// Physical constants, accessory kinds and gas-line inlet conditions shared by
/// SProd and its modules. The static_asserts at the end pin each value to its
/// spelling.
namespace sisprod {

// ---------------------------------------------------------------- units ----

/// Pascals per kgf/cm^2. Pressure is carried in kgf/cm^2 and needed in Pa.
inline constexpr double kPascalPerKgfPerCm2 = 98066.5;

/// The same conversion, spelled differently in some places. Using one value
/// for both would change results.
inline constexpr double kPascalPerKgfPerCm2Variant = 98066.52;

/// The PVTSim table reader's spelling. With it one atmosphere is 1.033211
/// kgf/cm^2, the value kAtmosphereInKgfPerCm2 keeps.
inline constexpr double kPascalPerKgfPerCm2PvtSim = 98068.059233;

/// A coarse spelling, 0.54% above 98066.5. The half-cell hydrostatic terms
/// beside a choke or valve use it, as do the interfacial work terms and the
/// pressure gradients of the gas and reverse steady temperature marches.
/// Preserved like the spellings above.
inline constexpr double kPascalPerKgfPerCm2Coarse = 98600.;

/// Kgf/cm^2 per pascal, 1 / 98066.5 written as a factor. The latent-heat
/// reader multiplies by it; a product and the division round differently.
inline constexpr double kKgfPerCm2PerPascal = 1.01971621e-5;

/// Psi per pascal, 1 / 6894.757.
inline constexpr double kPsiPerPascal = 0.00014503773800722;

/// Atmospheres per kgf/cm^2 (98066.5 / 101325) and psi per atmosphere. The
/// saturation correlation of the PVTSim tables goes from kgf/cm^2 to psi
/// through the atmosphere.
inline constexpr double kAtmospherePerKgfPerCm2 = 0.9678411;
inline constexpr double kPsiPerAtmosphere = 14.69595;

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

/// Psi per kgf/cm^2. Gas-lift valve correlations are imperial throughout.
inline constexpr double kPsiPerKgfPerCm2 = 14.223595;

/// One standard atmosphere in kgf/cm^2, subtracted to turn absolute pressure
/// into gauge before the imperial conversion.
inline constexpr double kAtmosphereInKgfPerCm2 = 1.033211;

/// One standard atmosphere in psi, added and subtracted to go between gauge
/// and absolute pressure in the gas-lift valve calibration. The PVTSim reader
/// rounds the same quantity as kPsiPerAtmosphere.
inline constexpr double kAtmosphereInPsi = 14.6959488;

/// Seconds per day. Flow rates are carried in kg/s and reported in m3/day.
inline constexpr double kSecondsPerDay = 86400.;

/// Pascal-seconds per centipoise. Viscosity correlations return cP; the
/// heat-transfer objects want SI.
inline constexpr double kPascalSecondPerCentipoise = 1.e-3;

/// Standard conditions as this program spells them, for gas density at
/// surface: one kgf/cm^2 and fifteen degrees Celsius.
inline constexpr double kStandardPressureKgfPerCm2 = 1.;
inline constexpr double kStandardTemperatureCelsius = 15.;

/// Celsius to Fahrenheit, written exactly as the sites it replaces spell it.
/// A function rather than two constants: the two factors are never used apart,
/// and naming the conversion is the point.
inline constexpr double celsiusToFahrenheit(double celsius) {
    return 1.8 * celsius + 32;
}

// ------------------------------------------------------------- numerics ----

/// Backward perturbation for numerical derivatives.
inline constexpr double kDerivativePerturbationFactor = 0.999;

/// Below this a phase-change mass rate is treated as zero.
inline constexpr double kPhaseChangeFloor = 1e-25;

/// The temperature range the solution is clamped to, in degrees Celsius.
inline constexpr double kMinimumTemperatureCelsius = -50.;
inline constexpr double kMaximumTemperatureCelsius = 200.;

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

static_assert(kPascalPerKgfPerCm2 == 98066.5);
static_assert(kPascalPerKgfPerCm2Variant == 98066.52);
static_assert(kPascalPerKgfPerCm2PvtSim == 98068.059233);
static_assert(kPascalPerKgfPerCm2Coarse == 98600.);
static_assert(kKgfPerCm2PerPascal == 1.01971621e-5);
static_assert(kPsiPerPascal == 0.00014503773800722);
static_assert(kAtmospherePerKgfPerCm2 == 0.9678411);
static_assert(kPsiPerAtmosphere == 14.69595);
static_assert(kBarrelPerCubicMetre == 6.29);
static_assert(kCubicFootPerCubicMetre == 35.31467);
static_assert(kAirDensityAtStandardConditions == 1.225);
static_assert(kGravity == 9.82);
static_assert(kGravityUnloadingVariant == 9.81);
static_assert(kPsiPerKgfPerCm2 == 14.223595);
static_assert(kAtmosphereInKgfPerCm2 == 1.033211);
static_assert(kAtmosphereInPsi == 14.6959488);
static_assert(kSecondsPerDay == 86400.);
static_assert(kPascalSecondPerCentipoise == 1.e-3);
static_assert(kStandardPressureKgfPerCm2 == 1.);
static_assert(kStandardTemperatureCelsius == 15.);
static_assert(celsiusToFahrenheit(0.) == 32.);
static_assert(celsiusToFahrenheit(100.) == 212.);
static_assert(kDerivativePerturbationFactor == 0.999);
static_assert(kPhaseChangeFloor == 1e-25);
static_assert(kMinimumTemperatureCelsius == -50.);
static_assert(kMaximumTemperatureCelsius == 200.);

}  // namespace sisprod

#endif  // SISPRODCONSTANTS_H_
