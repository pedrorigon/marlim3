#ifndef SISPRODGASLIFTINTERNAL_H_
#define SISPRODGASLIFTINTERNAL_H_

#include "SisProdConstants.h"
#include "SisProdGasLift.h"

/// The gas-lift module's own constants and types, which no other module uses.
namespace sisprod::gaslift {

/// The interface ratio at the middle of a cell: once the interface is past it, the cell counts as on the gas side.
inline constexpr double kMidCellInterfaceRatio = 0.5;

/// What a gas-lift valve is primed for: the fluid whose recovery factor turns the stagnation pressure into the
/// throat pressure, and the fluid whose mass flow the choke passes.
struct ValveFluids {
    InjectedFluid recovery;
    InjectedFluid flow;
};

/// Gas through the valve.
inline constexpr ValveFluids kGasThroughValve{kInjectedGas, kInjectedGas};
/// Completion fluid through the valve.
inline constexpr ValveFluids kCompletionFluidThroughValve{kInjectedCompletionFluid, kInjectedCompletionFluid};
/// The unloading search's combination: the gas recovery factor with the completion-fluid mass flow.
inline constexpr ValveFluids kUnloadingThroughValve{kInjectedGas, kInjectedCompletionFluid};

/// The transport properties of the fluid outside a production cell, as its heat transfer reads them.
struct ExternalFluidProperties {
    double conductivity;
    double specificHeat;
    double density;
    double viscosityCentipoise;
};

// Unloading the completion fluid through the gas-lift valves.

/// Gas-lift valve geometry. The bore threshold and the two spring rates are
/// imperial, as the whole valve correlation is.
inline constexpr double kLargeBoreValveInches = 1.1;
inline constexpr double kLargeBoreSpringRate = 500.0;
inline constexpr double kSmallBoreSpringRate = 1950.0;

/// The bellows of a calibrated valve: the temperature its dome charge is referred to, in degrees Fahrenheit, and
/// the fractional change of the dome pressure per degree Fahrenheit.
inline constexpr double kBellowsReferenceTemperatureFahrenheit = 80.;
inline constexpr double kBellowsPressureChangePerFahrenheit = 0.00215;

/// The unloading controller corrects the injection or the separator pressure when the completion-fluid flow
/// through the selected valve leaves the band from 1 - kUnloadingLowerBandMargin to 1 - kUnloadingUpperBandMargin
/// of the controlled flow rate.
inline constexpr double kUnloadingUpperBandMargin = 0.1;
inline constexpr double kUnloadingLowerBandMargin = 0.4;

/// The two edges of that band, as fractions of the controlled flow rate.
inline constexpr double kUnloadingUpperBandEdge = 1 - kUnloadingUpperBandMargin;
inline constexpr double kUnloadingLowerBandEdge = 1 - kUnloadingLowerBandMargin;

/// The largest pressure correction of one step, as a fraction of the corrected pressure per second of the step.
inline constexpr double kMaximumUnloadingCorrectionPerSecond = 0.01;

/// The separator pressure counts as at the unloading's maximum or minimum within these factors.
inline constexpr double kReachedMaximumPressureFactor = 0.9999999;
inline constexpr double kReachedMinimumPressureFactor = 1.0000001;

/// presMaxDesc before the minimum over the IPR accessories: above any pressure the column reaches, in kgf/cm^2.
inline constexpr double kUnboundedUnloadingPressure = 100000.;

/// The pressure the unloading hydrostatics start the gas line from when nothing else gives it, in kgf/cm^2.
inline constexpr double kDefaultUnloadingStartPressure = 10.;

/// The weight of the step's own flow against the moving average, in the flow the unloading controller corrects.
inline constexpr double kCurrentFlowWeight = 0.5;

/// The fraction of the separator pressure the unloading takes off per second while the column's last cell but
/// one holds no completion fluid.
inline constexpr double kSeparatorPressureDecayPerSecond = 0.01;

/// Before the unloading control starts, the separator pressure is brought down to the unloading's minimum by this
/// fraction of the latency, and held there after it.
inline constexpr double kSeparatorRampLatencyFraction = 0.5;

/// The factors a pressure is relaxed by per step towards its target, after the unloading.
inline constexpr double kRelaxationDownFactor = 0.95;
inline constexpr double kRelaxationUpFactor = 1.05;

/// The ratio of its cell the interface must reach in the last cell but one before it is handed over to the end of
/// the gas line.
inline constexpr double kInterfaceHandoverRatio = 0.99;

/// At or below this mixture velocity, in m/s, the unloading column counts as still and takes no wall friction.
inline constexpr double kStillMixtureVelocity = 1e-15;

// The gas line in the transient.

/// Below this opening the injection choke is throttling and becomes the mass
/// source of the first cell instead of a free boundary.
inline constexpr double kThrottlingChokeOpening = 0.2;

/// The temperature of the gas line's last cell while the flow there runs backwards, in degrees Celsius.
inline constexpr double kReverseFlowOutletTemperatureCelsius = 20.;

/// The gas line's master2 valve throttles while its throat area lies between these fractions of the smaller pipe
/// area around it.
inline constexpr double kMaster2ThrottlingAreaFraction = 0.01;
inline constexpr double kMaster2ClosedAreaFraction = 1e-5;

// The gas line in the steady state.

/// Fixed gas viscosity used to prime the steady annulus connection, in cP.
inline constexpr double kInitialAnnulusGasViscosityCentipoise = 0.16;

/// The steady gas temperature marches a cell only above this superficial gas velocity, in m/s, and between
/// diameters within this ratio of each other; elsewhere the cell takes the external temperature.
inline constexpr double kMinimumMarchGasVelocity = 1e-3;
inline constexpr double kMinimumMarchDiameterRatio = 0.5;

/// The steady gas temperature marches a cell in substeps when its length over its global thermal resistance
/// passes this value.
inline constexpr double kThermalSubstepLength = 1000.;

}  // namespace sisprod::gaslift

#endif  // SISPRODGASLIFTINTERNAL_H_
