#ifndef SISPRODTHERMALDETAIL_H_
#define SISPRODTHERMALDETAIL_H_

/// Result types for the kernels inside SisProdThermal.cpp.
///
/// These carry what one step of the thermal solve hands to the next: a
/// temperature balance, the source terms feeding it, the distributed
/// mass-transfer coefficients, and the smaller bundles the shared helpers
/// return. They are declared here rather than mid-file so the shape of the data
/// can be read without reading the arithmetic, and so the arithmetic in the
/// .cpp is uninterrupted by twenty-line type definitions.
///
/// This header is DETAIL: SisProdThermal.cpp is its only intended includer, and
/// nothing in it is part of the module's contract -- that is SisProdThermal.h.
/// The functions producing these values keep internal linkage in the .cpp; only
/// the types moved.

namespace sisprod::thermal {

/// The enthalpy the source term carries into the cell, and the temperature it
/// arrives at.
///
/// Seventy-nine lines of accessory dispatch lifted whole out of
/// computeMixtureEnthalpy, which is the only caller. The chain is NOT a switch
/// waiting to become a table: three arms test acsr.tipo, the fourth tests a
/// different member entirely (acsrL, the accessory on the left face) and then
/// branches again on ITS type and on how far the choke is closed. A dispatch
/// table keyed on tipo cannot express that, so the chain stays and only the
/// boundary moves. What the extraction buys is that the four values below are
/// now visibly the whole output of eighty lines, instead of four locals a
/// reader has to trace to be sure of.
struct SourceEnthalpy {
    double temperature;
    double gasEnthalpy;
    double liquidEnthalpy;
    /// Complementary-fluid heat. The legacy comment on three of the assignments
    /// says it is still to be corrected; the spelling is kept so that note
    /// still finds its subject.
    double hcF;
};

/// Cell geometry and the superficial velocities that both preparation steps
/// start from.
///
/// These thirty lines stood token for token identical in
/// prepareTemperatureBalance and prepareNonDimensionalHeatDiffusion -- the
/// second copy differs only by a blank line and one line break. What earns the
/// name is the branch: the velocities come from this cell unless the upstream
/// neighbour is a choke throttled below the active-area ratio, in which case
/// they come from the cell downstream.
struct CellFlowBasis {
    double flowArea;
    double voidFraction;
    double betmed;
    double gasSuperficialVelocity;
    double liquidSuperficialVelocity;
};

struct TemperatureBalance {
    double cellLength;
    double meanCellLength;
    double flowArea;
    double voidFraction;
    double bet;
    double gasSuperficialVelocity;
    double liquidSuperficialVelocity;
    double referenceMixtureVelocity;
    double rp;
    double rc;
    double liquidDensity;
    double gasDensity;
    double liquidHeatCapacity;
    double liquidIsochoricHeatCapacity;
    double gasHeatCapacity;
    double gasIsochoricHeatCapacity;
    double liquidJouleThomson;
    double gasJouleThomson;
    double hydrostaticPower;
    double heatFlux;
    double timeCoefficient;
    double pressureTimeCoefficient;
    double temperatureSpatialCoefficient;
    double pressureSpatialCoefficient;
};

struct TemperatureSourceTerms {
    double gas;
    double liquid;
};

struct DistributedMassTransferProperties {
    double downstreamWaterFraction;
    double upstreamWaterFraction;
    double cellWaterFraction;
    double liquidDensity;
    double gasDensity;
    double downstreamComposition;
    double upstreamComposition;
    double mixtureLiquidDensity;
    double downstreamOilVolumeFactor;
    double downstreamSolutionGasRatio;
    double downstreamSolutionGasPressureDerivative;
    double cellOilVolumeFactor;
    double cellSolutionGasRatio;
    double cellSolutionGasPressureDerivative;
    double cellSolutionGasTemperatureDerivative;
};

struct DistributedMassTransferCoefficients {
    double activeDerivative;
    double spatialCoupling;
    double flowArea;
};

/// The drift-flux pair of the slug regime: C0 = 1.2 and ud = 0.32*sqrt(g*D),
/// signed by the inclination, collapsing to the homogeneous limit (C0 = 1,
/// ud = 0) once the density ratio passes 0.9.
///
/// These fourteen lines stood in all four regime selectors, byte for byte. What
/// differs is the line BEFORE them, seeding meanVoidFraction from the pig field
/// or from the inlet, so the seeding stays at each call site and only the
/// proven-identical closure is shared -- the rule stage 3 settled on for the
/// five CalcC0Ud variants: keep the divergent control flow, share what a
/// normalised comparison proves equal.
struct SlugClosure {
    double c0;
    double ud;
    /// Read after the call by the interior selector only. Returned rather than
    /// passed by reference so the other three callers do not carry a variable
    /// they never read, which would trade duplication for a warning.
    double meanDiameter;
};

}  // namespace sisprod::thermal

#endif  // SISPRODTHERMALDETAIL_H_
