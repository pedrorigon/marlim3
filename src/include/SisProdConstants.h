#ifndef SISPRODCONSTANTS_H_
#define SISPRODCONSTANTS_H_

/// Physical constants and accessory kinds shared by the extracted modules.
///
/// This is a header rather than a block in one .cpp because the same values
/// appear in more than one place: 98066.5 is written 33 times in
/// SisProdThermal.cpp, twice in SisProdGasLift.cpp and 72 times in SisProd.cpp.
/// One definition means one place to read, and one place a mistake could hide.
///
/// Every constant keeps the EXACT spelling of the literal it replaces, and the
/// static_asserts at the end prove it. Naming a value must not move a bit --
/// that is the condition under which this refactoring is allowed to touch
/// anything at all.
namespace sisprod {

// ---------------------------------------------------------------- units ----

/// Pascals per kgf/cm^2. The input deck and the cell state carry pressure in
/// kgf/cm^2; the energy and momentum balances need Pa.
inline constexpr double kPascalPerKgfPerCm2 = 98066.5;

/// The same conversion as written in two other places. It is NOT the value
/// above, and the difference is preserved rather than harmonised: changing
/// either number changes results, and which one is intended is a question for
/// the technical owner, not for this refactoring.
inline constexpr double kPascalPerKgfPerCm2Variant = 98066.52;

/// Barrels per cubic metre and cubic feet per cubic metre. They appear as a
/// ratio, 6.29 / 35.31467, converting a solution gas-oil ratio from scf/bbl to
/// m3/m3. The division is left in place: the compiler folds it at compile time,
/// and writing the quotient by hand would round it to a different double.
inline constexpr double kBarrelPerCubicMetre = 6.29;
inline constexpr double kCubicFootPerCubicMetre = 35.31467;

/// Gravitational acceleration as this program uses it. Not 9.80665: the value
/// is 9.82 throughout, and "correcting" it would change every hydrostatic term.
inline constexpr double kGravity = 9.82;

// ------------------------------------------------------------- numerics ----

/// Backward perturbation used to take numerical derivatives: a property is
/// re-evaluated at 0.999 times the pressure or temperature and differenced.
inline constexpr double kDerivativePerturbationFactor = 0.999;

/// Below this, a phase-change mass rate is treated as zero rather than used to
/// derive a sign.
inline constexpr double kPhaseChangeFloor = 1e-25;

/// The temperature range the solution is clamped to, in degrees Celsius.
inline constexpr double kMinimumTemperatureCelsius = -50.;
inline constexpr double kMaximumTemperatureCelsius = 200.;

// ------------------------------------------------------ accessory kinds ----

/// The value of `celula[i].acsr.tipo`, which selects the accessory attached to
/// a cell and therefore which branch every source term takes.
///
/// The mapping is read from where the field is ASSIGNED -- Leitura.cpp and
/// LeituraVapor.cpp -- cross-checked against the member each branch then reads
/// (`acsr.injg`, `acsr.bcs`, ...) and the comments already in the code. Note
/// that acessorios.h carries a contradictory table mapping choke to 4; the
/// function that would apply it, carrega(const choke&), has no caller anywhere
/// in the tree, so it is dead and was not used as evidence.
///
/// Type 7 is deliberately absent: it is tested in this module and reads
/// `acsr.delp`, but no assignment names it, so what it stands for is not
/// established. Guessing a name is how a wrong one gets read as fact.
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
static_assert(kGravity == 9.82);
static_assert(kDerivativePerturbationFactor == 0.999);
static_assert(kPhaseChangeFloor == 1e-25);
static_assert(kMinimumTemperatureCelsius == -50.);
static_assert(kMaximumTemperatureCelsius == 200.);

}  // namespace sisprod

#endif  // SISPRODCONSTANTS_H_
