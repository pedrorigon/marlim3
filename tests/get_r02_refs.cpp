/*
 * get_r02_refs.cpp - Get reference values for R02 stream mixing functions
 */

#include "SisProd.h"
#include <cstdio>

int main() {
    using namespace marlim::sisprod2;
    
    std::puts("=== R02 Stream Mixing Reference Values ===\n");
    
    // Test blend function
    std::puts("--- blend() ---");
    StandardStream upstream;
    upstream.oilRate = 100.0;
    upstream.waterRate = 50.0;
    upstream.gasRate = 200.0;
    upstream.oilApi = 35.0;
    upstream.gasDensity = 0.65;
    upstream.waterDensity = 1.05;
    upstream.co2Fraction = 0.05;
    
    StandardStream source;
    source.oilRate = 20.0;
    source.waterRate = 10.0;
    source.gasRate = 50.0;
    source.oilApi = 30.0;
    source.gasDensity = 0.7;
    source.waterDensity = 1.02;
    source.co2Fraction = 0.08;
    
    SimContext ctx;
    const double tiny = 1e-30;
    
    StandardStream mixed = blend(upstream, source, tiny);
    std::printf("blend oilRate: %.12f\n", mixed.oilRate);
    std::printf("blend waterRate: %.12f\n", mixed.waterRate);
    std::printf("blend gasRate: %.12f\n", mixed.gasRate);
    std::printf("blend oilApi: %.12f\n", mixed.oilApi);
    std::printf("blend gasDensity: %.12f\n", mixed.gasDensity);
    std::printf("blend waterDensity: %.12f\n", mixed.waterDensity);
    std::printf("blend co2Fraction: %.12f\n", mixed.co2Fraction);
    
    // Test StreamMixing::march with different accessory types
    std::puts("\n--- StreamMixing::march() ---");
    
    StandardStream resultNone = StreamMixing::march(upstream, AccessoryType::None, source, ctx);
    std::printf("march None oilRate: %.12f\n", resultNone.oilRate);
    
    StandardStream resultGasSource = StreamMixing::march(upstream, AccessoryType::GasSource, source, ctx);
    std::printf("march GasSource oilRate: %.12f\n", resultGasSource.oilRate);
    
    // Test computePhaseTransferRate
    std::puts("\n--- computePhaseTransferRate() ---");
    PhaseTransferInput ptInput;
    ptInput.center.pressurePa = 100.0 * 98066.5;
    ptInput.center.temperatureK = 60.0 + 273.15;
    ptInput.center.pressureAuxPa = 90.0 * 98066.5;
    ptInput.center.liquidRate = 150.0;
    ptInput.center.waterCut = 0.3;
    ptInput.center.gasSpecificGravity = 0.65;
    ptInput.center.dissolvedGasDensityRatio = 1.0;
    ptInput.center.pigFraction = 0.0;
    ptInput.center.stream.oilRate = 100.0;
    ptInput.center.stream.waterRate = 50.0;
    ptInput.center.stream.gasRate = 200.0;
    ptInput.center.stream.oilApi = 35.0;
    ptInput.center.stream.gasDensity = 0.65;
    ptInput.center.stream.waterDensity = 1.05;
    ptInput.center.stream.co2Fraction = 0.05;
    
    ptInput.left = ptInput.center;
    ptInput.left.liquidRate = 140.0;
    ptInput.left.pressurePa = 95.0 * 98066.5;
    ptInput.left.pressureAuxPa = 85.0 * 98066.5;
    
    ptInput.right = ptInput.center;
    
    ptInput.cellLength = 100.0;
    ptInput.accessoryIsNone = true;
    
    double transferRate = computePhaseTransferRate(ptInput);
    std::printf("computePhaseTransferRate: %.12e kg/m3/s\n", transferRate);
    
    // Test with accessoryIsNone = false (should return 0)
    ptInput.accessoryIsNone = false;
    double transferRateZero = computePhaseTransferRate(ptInput);
    std::printf("computePhaseTransferRate (no accessory): %.12e kg/m3/s\n", transferRateZero);
    
    return 0;
}
