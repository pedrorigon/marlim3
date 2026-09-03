#ifndef SISPRODCONSTANTS_H_
#define SISPRODCONSTANTS_H_

/// Physical constants and accessory kinds shared by the extracted modules.
/// Each keeps the exact spelling of the literal it replaces; the static_asserts
/// at the end prove it.
namespace sisprod {

// ---------------------------------------------------------------- units ----

/// Pascals per kgf/cm^2. Pressure is carried in kgf/cm^2 and needed in Pa.
inline constexpr double kPascalPerKgfPerCm2 = 98066.5;

/// The same conversion, spelled differently elsewhere in the tree. The
/// difference is preserved: harmonising it would change results.
inline constexpr double kPascalPerKgfPerCm2Variant = 98066.52;

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
enum : int {
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

static_assert(kPascalPerKgfPerCm2 == 98066.5);
static_assert(kPascalPerKgfPerCm2Variant == 98066.52);
static_assert(kBarrelPerCubicMetre == 6.29);
static_assert(kCubicFootPerCubicMetre == 35.31467);
static_assert(kAirDensityAtStandardConditions == 1.225);
static_assert(kGravity == 9.82);
static_assert(kDerivativePerturbationFactor == 0.999);
static_assert(kPhaseChangeFloor == 1e-25);
static_assert(kMinimumTemperatureCelsius == -50.);
static_assert(kMaximumTemperatureCelsius == 200.);

}  // namespace sisprod

#endif  // SISPRODCONSTANTS_H_
