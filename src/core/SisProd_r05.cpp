/*
 * SisProd_r05.cpp
 *
 * R05 — Two-phase pressure-gradient correlation dispatcher.
 *
 * Thin wrapper around the existing GradientCorrelations free functions.
 * Those functions are already pure C (no SProd state), so this file just
 * provides:
 *   - PressureGradientInput  : typed parameter bundle
 *   - PressureGradientResult : typed output bundle
 *   - computePressureGradient(): dispatcher by FlowCorrelationId
 *
 * IMPORTANT: GradientCorrelations.h pulls in FerramentasNumericas.h,
 * Log.h and other headers that require the full CMake build environment.
 * This file is therefore only compiled into the main CMake target
 * (where those headers are available) and NOT into the standalone
 * sisprod2_selftest harness.
 *
 * Migration reference: issues/sisprod-migration-plan.md, region R05.
 */

#include "SisProd2.h"
#include "GradientCorrelations.h"

#include <cstring>

namespace marlim {
namespace sisprod2 {

// ===========================================================================
// computePressureGradient
// ===========================================================================
PressureGradientResult
computePressureGradient(FlowCorrelationId correlationId,
                        const PressureGradientInput &in) {
    PressureGradientResult out;

    // Aliases to improve readability inside each branch.
    const double angle  = in.inclinationRad;
    const double dia    = in.diameter;
    const double rug    = in.roughness;
    const double pres   = in.pressure;
    const double vel    = in.mixtureVelocity;
    const double liqF   = in.liquidFraction;
    const double rhog   = in.gasDensity;
    const double rhol   = in.liquidDensity;
    const double muG    = in.gasViscosity;
    const double muL    = in.liquidViscosity;
    const double sigma  = in.surfaceTension;
    const double temp   = in.temperature;
    const double Z      = in.compressibilityFactor;
    const double wFrac  = in.waterFraction;
    const double qProd  = in.productionRate;

    double &holdup      = out.holdup;
    double &fricG       = out.frictionGrad;
    double &gravG       = out.gravityGrad;
    double &accelG      = out.accelGrad;
    double &totalG      = out.totalGrad;
    double &re          = out.reynolds;
    unsigned char &fType = out.flowType;

    unsigned char critFlag    = 0;
    unsigned char convFlag    = 0;
    unsigned char isigFlag    = 0;
    unsigned char palmerFlag  = 0;

    switch (correlationId) {
    case FlowCorrelationId::PoettmannCarpenter:
        poettmannCarpenter(angle, dia, vel, liqF, rhog, rhol, muG, muL, rug,
                           holdup, fricG, gravG, totalG, re, fType);
        break;
    case FlowCorrelationId::BaxendellThomas:
        baxendellThomas(angle, dia, vel, liqF, rhog, rhol, muG, muL, rug,
                        holdup, fricG, gravG, totalG, re, fType);
        break;
    case FlowCorrelationId::FancherBrown:
        fancherBrown(angle, dia, qProd, vel, liqF, rhog, rhol, muG, muL, rug,
                     holdup, fricG, gravG, totalG, re, fType);
        break;
    case FlowCorrelationId::HagedornBrown:
        hagedornBrown(angle, dia, rug, pres, vel, liqF, rhog, rhol, muG, muL,
                      sigma, temp, Z,
                      holdup, fricG, gravG, accelG, totalG, re, fType,
                      critFlag, convFlag);
        break;
    case FlowCorrelationId::DunsRos:
        dunsRos(angle, dia, rug, pres, vel, liqF, rhog, rhol, muG, muL,
                sigma, temp, Z,
                holdup, fricG, gravG, accelG, totalG, re, fType,
                critFlag, convFlag);
        break;
    case FlowCorrelationId::Orkiszewski:
        orkiszewski(angle, dia, rug, pres, vel, liqF, rhog, rhol, muG, muL,
                    sigma, temp, Z, wFrac,
                    holdup, fricG, gravG, accelG, totalG, re, fType,
                    critFlag, convFlag, isigFlag);
        break;
    case FlowCorrelationId::BeggsAndBrill:
        beggsAndBrill(angle, dia, rug, pres, vel, liqF, rhog, rhol, muG, muL,
                      sigma,
                      holdup, fricG, gravG, accelG, totalG, re, fType,
                      palmerFlag, critFlag);
        break;
    case FlowCorrelationId::MukherjeeeBrill:
        mukherjeeeBrill(angle, dia, rug, pres, vel, liqF, rhog, rhol, muG, muL,
                        sigma, temp, Z,
                        holdup, fricG, gravG, accelG, totalG, re, fType,
                        critFlag, convFlag);
        break;
    case FlowCorrelationId::Aziz:
        aziz(angle, dia, rug, pres, vel, liqF, rhog, rhol, muG, muL, sigma,
             holdup, fricG, gravG, accelG, totalG, re, fType, critFlag);
        break;
    case FlowCorrelationId::Gray:
        gray(angle, dia, rug, pres, vel, liqF, wFrac, rhog, rhol, muG, muL,
             sigma, sigma, // oilSurfaceTension = waterSurfaceTension = sigma
             holdup, fricG, gravG, accelG, totalG, re, fType, critFlag);
        break;
    default:
        // Unknown correlation — return zero gradients.
        out = PressureGradientResult();
        out.status = SolveStatus::OutOfRange;
        break;
    }
    return out;
}

// ===========================================================================
// ProductionColumn::marchToWellheadPhysical
//
// R05 physics-complete pressure march using a full two-phase correlation
// (default: Beggs & Brill, the most widely used in production engineering).
//
// The implementation follows the same sequential pattern as marchToWellhead
// but drives pressure drop through computePressureGradient() instead of the
// simplified Darcy homogeneous model.  At each segment:
//   1. StreamMixing updates composition from any accessory source.
//   2. PressureGradientInput is filled from the current stream + fluid props.
//   3. computePressureGradient() returns holdup, frictionGrad, gravityGrad.
//   4. Pressure is decremented by (gravityGrad + frictionGrad) * segmentLength.
//
// This loop is inherently sequential (pressure at cell i depends on cell i-1)
// and must NOT be parallelised at this level.
// ===========================================================================
double ProductionColumn::marchToWellheadPhysical(
    double bottomholePressure, const SimContext &context,
    FlowCorrelationId correlationId,
    std::vector<ProfilePoint> *profileOut) const {

    double pressure = bottomholePressure;
    StandardStream stream = inletStream_;

    if (profileOut) {
        profileOut->clear();
        profileOut->reserve(segments_.size());
    }

    for (std::size_t i = 0; i < segments_.size(); ++i) {
        const PipeSegment &segment = segments_[i];

        // Update stream composition at any accessory (gas-lift, liquid source).
        if (segment.accessory != AccessoryType::None) {
            stream = StreamMixing::march(stream, segment.accessory,
                                         segment.sourceStream, context);
        }

        // Estimate local fluid properties for the pressure-gradient input.
        const double qualityApprox =
            fluid_.massQualityFromGasOilRatio(stream.gasOilRatio());
        const double densityP =
            pressure > constants::kPressureFloor ? pressure
                                                 : constants::kPressureFloor;
        const double rhomix = fluid_.mixtureDensity(densityP, qualityApprox);
        const double area    = segment.crossSectionArea();
        const double vel     = massFlowRate_ > 0.0
                               ? massFlowRate_ / (rhomix * area + 1e-15)
                               : 0.0;

        // Build the input bundle for the correlation.
        PressureGradientInput in;
        in.inclinationRad         = segment.inclination;
        in.diameter               = segment.diameter;
        in.roughness              = 1e-4; // default absolute roughness (m)
        in.pressure               = densityP;
        in.mixtureVelocity        = vel;
        in.liquidFraction         = 1.0 - qualityApprox; // no-slip liquid fraction
        in.gasDensity             = fluid_.gasDensityAt(densityP);
        in.liquidDensity          = fluid_.liquidDensity();
        in.gasViscosity           = 1.5e-5; // approx gas viscosity (Pa·s)
        in.liquidViscosity        = fluid_.viscosity();
        in.surfaceTension         = 0.025;  // approx oil-gas surface tension (N/m)
        in.temperature            = context.simulationTime() > 0.0 ? 330.0 : 330.0;
        in.compressibilityFactor  = 1.0;
        in.waterFraction          = stream.basicSedimentAndWater();
        in.productionRate         = massFlowRate_ / (rhomix + 1e-15);

        PressureGradientResult grad = computePressureGradient(correlationId, in);

        const double pressureDrop = (grad.gravityGrad + grad.frictionGrad)
                                    * segment.length;

        if (profileOut) {
            ProfilePoint point;
            point.pressure       = pressure;
            point.mixtureDensity = rhomix;
            point.velocity       = vel;
            point.gasOilRatio    = stream.gasOilRatio();
            profileOut->push_back(point);
        }

        pressure -= pressureDrop;
        if (!std::isfinite(pressure))
            return pressure;
    }
    return pressure;
}

} // namespace sisprod2
} // namespace marlim
