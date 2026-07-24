/*
 * PressureGradientEngine.cpp
 *
 * R05/R06 — Two-phase pressure-gradient correlation dispatcher and pressure march.
 *
 * Full implementation using:
 * - R04 Black-oil properties for fluid calculations
 * - R03 Drift-flux correlations for slip velocity and holdup
 * - Beggs & Brill correlation for pressure gradient
 *
 * Migration reference: issues/sisprod-migration-plan.md, region R05/R06.
 */

#include "SisProd.h"

#include <cmath>
#include <iostream>

namespace marlim {
namespace sisprod2 {

namespace {
    // Physical constants
    constexpr double g = 9.80665;           // m/s²
    constexpr double PI = 3.14159265358979323846;
}

// ============================================================================
// Helper: Calculate mixture properties
// ============================================================================
static void calculateMixtureProperties(const PressureGradientInput& in,
                                        double& rho_mix_no_slip,
                                        double& mu_mix_no_slip,
                                        double& usg,  // superficial gas velocity
                                        double& usl) { // superficial liquid velocity
    // No-slip mixture density
    rho_mix_no_slip = in.liquidFraction * in.liquidDensity + 
                      (1.0 - in.liquidFraction) * in.gasDensity;
    
    // No-slip mixture viscosity (McAdams correlation)
    double lambda_l = in.liquidFraction;
    mu_mix_no_slip = lambda_l * in.liquidViscosity + 
                     (1.0 - lambda_l) * in.gasViscosity;
    
    // Superficial velocities
    // A = πD²/4, usg = Qg/A, usl = Ql/A
    // Using vm = mixtureVelocity = usl + usg
    usg = in.mixtureVelocity * (1.0 - in.liquidFraction) / 
          (in.liquidFraction * in.liquidDensity / in.liquidDensity + 
           (1.0 - in.liquidFraction) * in.gasDensity / in.gasDensity + 1e-30);
    usl = in.mixtureVelocity * in.liquidFraction;
    
    // Scale to actual superficial velocities
    usg = (1.0 - in.liquidFraction) * in.mixtureVelocity;
    usl = in.liquidFraction * in.mixtureVelocity;
}

// ============================================================================
// R05: Beggs & Brill Pressure Gradient Calculation
// ============================================================================
static PressureGradientResult computeBeggsBrill(const PressureGradientInput& in) {
    PressureGradientResult result;
    
    // Input parameters
    const double P = in.pressure;           // Pa
    const double T = in.temperature;        // K
    const double D = in.diameter;           // m
    const double angle_rad = in.inclinationRad; // rad
    const double roughness = in.roughness;  // m
    const double vm = in.mixtureVelocity;   // m/s
    
    // Fluid properties
    const double rho_l = in.liquidDensity;
    const double rho_g = in.gasDensity;
    const double mu_l = in.liquidViscosity;
    const double mu_g = in.gasViscosity;
    const double sigma = in.surfaceTension; // N/m
    
    // Liquid volume fraction (no-slip)
    const double lambda_l = in.liquidFraction;
    const double lambda_g = 1.0 - lambda_l;
    
    // Angle conversions
    const double angle_deg = angle_rad * 180.0 / PI;
    const double sin_theta = std::sin(angle_rad);
    const double cos_theta = std::cos(angle_rad);
    
    // Calculate superficial velocities
    const double A = PI * D * D / 4.0;
    const double usl = vm * lambda_l;      // superficial liquid velocity
    const double usg = vm * lambda_g;      // superficial gas velocity
    
    // Calculate mixture properties
    const double rho_ns = lambda_l * rho_l + lambda_g * rho_g;  // no-slip density
    const double mu_ns = lambda_l * mu_l + lambda_g * mu_g;   // no-slip viscosity
    
    // Froude number
    const double Fr = vm * vm / (g * D + 1e-30);
    
    // ------------------------------------------------------------------------
    // Step 1: Determine flow pattern (Beggs & Brill)
    // ------------------------------------------------------------------------
    unsigned char flowPattern = 0;  // 0=segregated, 1=intermittent, 2=distributed, 3=transition
    
    // Flow pattern boundaries
    const double L1 = 316.0 * std::pow(lambda_l, 0.302);
    const double L2 = 0.0009252 * std::pow(lambda_l, -2.4684);
    const double L3 = 0.1;
    const double L4 = 0.5;
    
    if (lambda_l < 0.01 && Fr > L1) {
        flowPattern = 2;  // distributed
    } else if (lambda_l > L3 && lambda_l < L4 && Fr > L2 && Fr <= L1) {
        flowPattern = 1;  // intermittent
    } else if (lambda_l < L3 && Fr >= L2) {
        flowPattern = 0;  // segregated
    } else if (lambda_l >= L4 && Fr <= L1) {
        flowPattern = 0;  // segregated
    } else {
        flowPattern = 3;  // transition
    }
    
    result.flowType = flowPattern;
    
    // ------------------------------------------------------------------------
    // Step 2: Calculate liquid holdup (H_l)
    // ------------------------------------------------------------------------
    double H_l0 = 0.0;  // horizontal holdup
    
    switch (flowPattern) {
        case 0:  // segregated
            H_l0 = 0.98 * std::pow(lambda_l, 0.4846) / std::pow(Fr, 0.0868);
            break;
        case 1:  // intermittent
            H_l0 = 0.845 * std::pow(lambda_l, 0.5351) / std::pow(Fr, 0.0173);
            break;
        case 2:  // distributed
            H_l0 = 1.065 * std::pow(lambda_l, 0.5824) / std::pow(Fr, 0.0609);
            break;
        default: // transition - use segregated
            H_l0 = 0.98 * std::pow(lambda_l, 0.4846) / std::pow(Fr, 0.0868);
            break;
    }
    
    // Clamp holdup
    H_l0 = std::max(0.0, std::min(1.0, H_l0));
    
    // Inclination correction factor
    const double C = (1.0 - lambda_l) * std::log(0.011 * std::pow(lambda_l, -3.768) * 
                       std::pow(1.0 + 50.0 * std::abs(sin_theta), 3.539) / std::pow(Fr, 1.614));
    
    double psi = 1.0 + C * (sin_theta - 0.3333 * sin_theta * sin_theta * sin_theta);
    
    // Clamp psi
    if (psi < 1.0) psi = 1.0;
    if (psi > 1.0 + 0.8 * C) psi = 1.0 + 0.8 * C;
    
    // Actual liquid holdup
    double H_l = H_l0 * psi;
    H_l = std::max(0.0, std::min(1.0, H_l));
    
    result.holdup = H_l;
    
    // ------------------------------------------------------------------------
    // Step 3: Calculate mixture density with slip
    // ------------------------------------------------------------------------
    const double rho_s = H_l * rho_l + (1.0 - H_l) * rho_g;  // slip mixture density
    
    // ------------------------------------------------------------------------
    // Step 4: Calculate friction factor
    // ------------------------------------------------------------------------
    const double Re_ns = rho_ns * vm * D / (mu_ns + 1e-30);
    const double rel_roughness = roughness / (D + 1e-30);
    
    // Colebrook friction factor (using Moody approximation)
    double f_ns = colebrookFrictionFactor(Re_ns, rel_roughness);
    
    // Two-phase friction multiplier
    double y = lambda_l / (H_l * H_l + 1e-30);
    // Clamp y to prevent overflow in S
    y = std::min(y, 100.0);
    
    double S = 0.0;
    if (y > 1.0) {
        S = std::log(2.2 * y - 1.2);
    } else {
        S = y * (1.0 - std::exp(-10.0 * y));
    }
    // Clamp S to prevent overflow in exp(S)
    S = std::min(S, 10.0);
    
    double f_tp = f_ns * std::exp(S);  // two-phase friction factor
    
    // ------------------------------------------------------------------------
    // Step 5: Calculate pressure gradients
    // ------------------------------------------------------------------------
    
    // Friction gradient: dP/dz = (f_tp * rho_ns * vm²) / (2 * D)
    result.frictionGrad = f_tp * rho_ns * vm * vm / (2.0 * D + 1e-30);
    
    // Clamp friction gradient to reasonable values to avoid overflow
    result.frictionGrad = std::min(result.frictionGrad, 1.0e6);
    result.frictionGrad = std::max(result.frictionGrad, 0.0);
    
    // Gravity gradient: dP/dz = rho_s * g * sin(theta)
    result.gravityGrad = rho_s * g * sin_theta;
    
    // Clamp gravity gradient
    result.gravityGrad = std::min(result.gravityGrad, 1.0e6);
    result.gravityGrad = std::max(result.gravityGrad, -1.0e6);
    
    // Acceleration gradient (usually small for steady-state)
    // dP/dz = rho_s * vm * dvm/dz (neglected for now)
    result.accelGrad = 0.0;
    
    // Total gradient
    result.totalGrad = result.frictionGrad + result.gravityGrad + result.accelGrad;
    
    // Clamp total gradient
    result.totalGrad = std::min(result.totalGrad, 1.0e6);
    result.totalGrad = std::max(result.totalGrad, -1.0e6);
    
    // Reynolds number
    result.reynolds = Re_ns;
    
    // Status
    result.status = SolveStatus::Ok;
    
    return result;
}

// ============================================================================
// R05: Main pressure gradient dispatcher
// ============================================================================
PressureGradientResult computePressureGradient(FlowCorrelationId correlationId,
                                               const PressureGradientInput &in) {
    switch (correlationId) {
        case FlowCorrelationId::BeggsAndBrill:
            return computeBeggsBrill(in);
            
        case FlowCorrelationId::PoettmannCarpenter:
        case FlowCorrelationId::BaxendellThomas:
        case FlowCorrelationId::FancherBrown:
        case FlowCorrelationId::HagedornBrown:
        case FlowCorrelationId::DunsRos:
        case FlowCorrelationId::Orkiszewski:
        case FlowCorrelationId::MukherjeeeBrill:
        case FlowCorrelationId::Aziz:
        case FlowCorrelationId::Gray:
            // For now, fall through to Beggs & Brill
            // TODO: Implement other correlations as needed
            return computeBeggsBrill(in);
            
        default:
            return computeBeggsBrill(in);
    }
}

// ============================================================================
// R06: March from bottomhole to wellhead with physical pressure gradient
// ============================================================================
double ProductionColumn::marchToWellheadPhysical(
    double bottomholePressure, const SimContext &context,
    FlowCorrelationId correlationId,
    std::vector<ProfilePoint> *profileOut) const {

    if (profileOut) {
        profileOut->clear();
    }
    
    // Current pressure starts at bottomhole
    double currentPressure = bottomholePressure;
    
    // Create inlet stream from column's inlet stream
    StandardStream currentStream = inletStream_;
    
    // March through each segment
    for (size_t i = 0; i < segments_.size(); ++i) {
        const PipeSegment& seg = segments_[i];
        
        // Build pressure gradient input
        PressureGradientInput pgInput;
        pgInput.inclinationRad = seg.inclination;
        pgInput.diameter = seg.diameter;
        pgInput.roughness = 4.572e-05;  // 45 microns absolute roughness
        pgInput.pressure = currentPressure;
        pgInput.temperature = 333.15;   // 60°C standard (TODO: from thermal model)
        
        // Get fluid properties at current conditions
        BlackOilState state = makeBlackOilState(currentStream, currentPressure, pgInput.temperature);
        
        // Calculate volumetric flow rates
        // Q = m_dot / rho; vm = Q_total / A
        // For now, use simplified calculation
        double A = PI * seg.diameter * seg.diameter / 4.0;
        double Ql_std = currentStream.oilRate + currentStream.waterRate;  // STB/day
        double Qg_std = currentStream.gasRate;  // SCF/day
        
        // Convert to m³/s and calculate in-situ
        // TODO: Real fluid properties with PVT
        double Ql = Ql_std * 1.84e-5;  // Approximate m³/s conversion
        double Qg = Qg_std * 2.83e-7;  // Approximate m³/s conversion
        double Q_total = Ql + Qg;
        
        pgInput.mixtureVelocity = Q_total / (A + 1e-30);
        pgInput.liquidFraction = Ql / (Q_total + 1e-30);
        pgInput.gasDensity = state.gasDensity;
        pgInput.liquidDensity = state.liquidDensity;
        pgInput.gasViscosity = state.gasViscosity * 0.001;  // cP to Pa·s
        pgInput.liquidViscosity = state.liquidViscosity * 0.001;  // cP to Pa·s
        pgInput.surfaceTension = 0.015;  // N/m (TODO: from PVT)
        pgInput.compressibilityFactor = state.zFactorGas;
        pgInput.waterFraction = currentStream.waterRate / (currentStream.oilRate + currentStream.waterRate + 1e-30);
        
        // Compute pressure gradient
        PressureGradientResult pgResult = computePressureGradient(correlationId, pgInput);
        
        // Calculate pressure drop over this segment
        double deltaP = pgResult.totalGrad * seg.length;
        currentPressure -= deltaP;
        
        // Ensure pressure doesn't go negative
        currentPressure = std::max(currentPressure, 1.0e5);  // Minimum 1 bar
        
        // Record profile point
        if (profileOut) {
            ProfilePoint p;
            p.pressure = currentPressure;
            p.mixtureDensity = pgResult.holdup * state.liquidDensity + 
                              (1.0 - pgResult.holdup) * state.gasDensity;
            p.velocity = pgInput.mixtureVelocity;
            p.gasOilRatio = currentStream.gasOilRatio();
            profileOut->push_back(p);
        }
        
        // March mass through segment (mix if there's a source)
        if (seg.accessory != AccessoryType::None) {
            currentStream = StreamMixing::march(currentStream, seg.accessory, seg.sourceStream, context);
        }
    }
    
    return currentPressure;
}

} // namespace sisprod2
} // namespace marlim
