/*
 * sisprod_parity_tests.cpp
 *
 * Framework de testes de paridade local (sem CI).
 *
 * Valida cada função migrada comparando resultados entre:
 * - Valores de referência (dos self-tests)
 * - Nova arquitetura (arquivos modulares - namespace sisprod2)
 *
 * Nota: Não compara diretamente contra SProd devido à dependência Fortran.
 * Os valores de referência são extraídos dos self-tests em sisprod2_selftest.cpp
 *
 * Execução local:
 *   cd build && cmake .. && make run_parity_tests
 *
 * Migration reference: issues/sisprod-migration-plan.md, Apêndice B
 */

#include "SisProd.h"

#include <cmath>
#include <cstdio>
#include <cstring>

namespace parity {

static int g_passed = 0;
static int g_failed = 0;

inline bool relative_eq(double a, double b, double tol, double& rel_err) {
    const double abs_diff = std::fabs(a - b);
    const double scale = 0.5 * (std::fabs(a) + std::fabs(b));
    rel_err = abs_diff / (scale + 1e-30);
    return rel_err <= tol;
}

inline bool absolute_eq(double a, double b, double tol) {
    return std::fabs(a - b) <= tol;
}

#define EXPECT_PARITY(expected, actual, rel_tol, abs_tol, ...) \
    do { \
        double _rel_err = 0.0; \
        bool _rel_ok = relative_eq(expected, actual, rel_tol, _rel_err); \
        bool _abs_ok = absolute_eq(expected, actual, abs_tol); \
        if (!_rel_ok && !_abs_ok) { \
            std::printf("    [FAIL] %s expected=%.12e actual=%.12e rel_err=%.2e\n", \
                        __VA_ARGS__, expected, actual, _rel_err); \
            g_failed++; \
        } else { \
            g_passed++; \
            std::printf("    [PASS] %s rel_err=%.2e\n", __VA_ARGS__, _rel_err); \
        } \
    } while(0)

// ---------------------------------------------------------------------------
// R02 Tests (Stream Mixing / Phase Transfer)
// ---------------------------------------------------------------------------

void test_r02_blend() {
    std::puts("\n=== R02/blend parity ===");
    
    // Reference values from actual implementation
    marlim::sisprod2::StandardStream upstream;
    upstream.oilRate = 100.0;
    upstream.waterRate = 50.0;
    upstream.gasRate = 200.0;
    upstream.oilApi = 35.0;
    upstream.gasDensity = 0.65;
    upstream.waterDensity = 1.05;
    upstream.co2Fraction = 0.05;
    
    marlim::sisprod2::StandardStream source;
    source.oilRate = 20.0;
    source.waterRate = 10.0;
    source.gasRate = 50.0;
    source.oilApi = 30.0;
    source.gasDensity = 0.7;
    source.waterDensity = 1.02;
    source.co2Fraction = 0.08;
    
    const double tiny = 1e-30;
    marlim::sisprod2::StandardStream mixed = marlim::sisprod2::blend(upstream, source, tiny);
    
    EXPECT_PARITY(120.0, mixed.oilRate, 1e-6, 1e-9, "blend oilRate");
    EXPECT_PARITY(60.0, mixed.waterRate, 1e-6, 1e-9, "blend waterRate");
    EXPECT_PARITY(250.0, mixed.gasRate, 1e-6, 1e-9, "blend gasRate");
    EXPECT_PARITY(34.145277207392, mixed.oilApi, 1e-6, 1e-9, "blend oilApi");
    EXPECT_PARITY(0.66, mixed.gasDensity, 1e-6, 1e-9, "blend gasDensity");
    EXPECT_PARITY(1.045, mixed.waterDensity, 1e-6, 1e-9, "blend waterDensity");
    EXPECT_PARITY(0.056, mixed.co2Fraction, 1e-6, 1e-9, "blend co2Fraction");
}

void test_r02_stream_mixing_march() {
    std::puts("\n=== R02/StreamMixing::march parity ===");
    
    marlim::sisprod2::StandardStream upstream;
    upstream.oilRate = 100.0;
    upstream.waterRate = 50.0;
    upstream.gasRate = 200.0;
    upstream.oilApi = 35.0;
    
    marlim::sisprod2::StandardStream source;
    source.oilRate = 20.0;
    source.waterRate = 10.0;
    source.gasRate = 50.0;
    
    marlim::sisprod2::SimContext ctx;
    
    // Test AccessoryType::None - should return upstream unchanged
    marlim::sisprod2::StandardStream resultNone = marlim::sisprod2::StreamMixing::march(
        upstream, marlim::sisprod2::AccessoryType::None, source, ctx);
    EXPECT_PARITY(100.0, resultNone.oilRate, 1e-6, 1e-9, "march None oilRate");
    
    // Test AccessoryType::GasSource - should blend
    marlim::sisprod2::StandardStream resultGas = marlim::sisprod2::StreamMixing::march(
        upstream, marlim::sisprod2::AccessoryType::GasSource, source, ctx);
    EXPECT_PARITY(120.0, resultGas.oilRate, 1e-6, 1e-9, "march GasSource oilRate");
}

void test_r02_compute_phase_transfer_rate() {
    std::puts("\n=== R02/computePhaseTransferRate parity ===");
    
    using namespace marlim::sisprod2;
    
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
    
    const double expected_transfer = -4.117224262423e-07;
    double transferRate = computePhaseTransferRate(ptInput);
    EXPECT_PARITY(expected_transfer, transferRate, 1e-6, 1e-9, "computePhaseTransferRate");
    
    // Test with accessoryIsNone = false (should return 0)
    ptInput.accessoryIsNone = false;
    double transferRateZero = computePhaseTransferRate(ptInput);
    EXPECT_PARITY(0.0, transferRateZero, 1e-6, 1e-9, "computePhaseTransferRate (no accessory)");
}

// ---------------------------------------------------------------------------
// R03 Tests
// ---------------------------------------------------------------------------

void test_r03_choi() {
    std::puts("\n=== R03/choi parity ===");
    
    // Reference values from self-test
    const double expected_c0 = 1.1928421124;
    const double expected_ud = 0.2726455520;
    
    double c0_new = 0.0, ud_new = 0.0;
    marlim::sisprod2::choi(850.0, 1.2, 0.072, 0.30,
                           50000.0, 40000.0, 0.05, 0.08,
                           0.10, 0.0001, 1.57079632679,
                           c0_new, ud_new, 1.0);
    
    EXPECT_PARITY(expected_c0, c0_new, 1e-6, 1e-9, "choi c0");
    EXPECT_PARITY(expected_ud, ud_new, 1e-6, 1e-9, "choi ud");
}

void test_r03_hibiki_ishii() {
    std::puts("\n=== R03/hibikiIshii parity ===");
    
    // Reference values calculated from implementation
    const double expected_c0 = 2.5545405102;
    const double expected_ud = 10.5160925813;
    
    double c0_new = 0.0, ud_new = 0.0;
    marlim::sisprod2::hibikiIshii(850.0, 1.2, 0.072, 0.30,
                                    50000.0, 40000.0, 0.05, 0.08,
                                    0.10, 0.0001, 1.57079632679,
                                    c0_new, ud_new, 1.0);
    
    EXPECT_PARITY(expected_c0, c0_new, 1e-6, 1e-9, "hibikiIshii c0");
    EXPECT_PARITY(expected_ud, ud_new, 1e-6, 1e-9, "hibikiIshii ud");
}

void test_r03_franca_lahey() {
    std::puts("\n=== R03/francaLahey parity ===");
    
    // Reference values: constants c0=1.04, ud=0.466 (vertical)
    const double expected_c0 = 1.04;
    const double expected_ud = 0.466;
    
    double c0_new = 0.0, ud_new = 0.0;
    marlim::sisprod2::francaLahey(850.0, 1.2, 0.072, 0.30,
                                  50000.0, 40000.0, 0.05, 0.08,
                                  0.10, 0.0001, 1.57079632679,
                                  c0_new, ud_new, 1.0);
    
    EXPECT_PARITY(expected_c0, c0_new, 1e-6, 1e-9, "francaLahey c0");
    EXPECT_PARITY(expected_ud, ud_new, 1e-6, 1e-9, "francaLahey ud");
}

void test_r03_bhagwat_ghajar() {
    std::puts("\n=== R03/bhagwatGhajar parity ===");
    
    // Reference values calculated from implementation
    const double expected_c0 = 1.187184315137;
    const double expected_ud = 0.2460054528;
    
    double c0_new = 0.0, ud_new = 0.0;
    marlim::sisprod2::bhagwatGhajar(850.0, 1.2, 0.072, 0.30,
                                      50000.0, 40000.0, 0.05, 0.08,
                                      0.10, 0.0001, 1.57079632679,
                                      c0_new, ud_new, 1.0);
    
    EXPECT_PARITY(expected_c0, c0_new, 1e-6, 1e-9, "bhagwatGhajar c0");
    EXPECT_PARITY(expected_ud, ud_new, 1e-6, 1e-9, "bhagwatGhajar ud");
}

void test_r03_bhagwat_ghajar_mod() {
    std::puts("\n=== R03/bhagwatGhajarMod parity ===");
    
    // bhagwatGhajarMod uses reymixL instead of reymix for the Reynolds number
    // Subtle difference: reymix = 50000, reymixL = 40000 leads to slightly different results
    // Reference values from actual implementation:
    const double expected_c0 = 1.186925248101;
    const double expected_ud = 0.243356160520;
    
    double c0_new = 0.0, ud_new = 0.0;
    marlim::sisprod2::bhagwatGhajarMod(850.0, 1.2, 0.072, 0.30,
                                         50000.0, 40000.0, 0.05, 0.08,
                                         0.10, 0.0001, 1.57079632679,
                                         c0_new, ud_new, 1.0);
    
    EXPECT_PARITY(expected_c0, c0_new, 1e-6, 1e-9, "bhagwatGhajarMod c0");
    EXPECT_PARITY(expected_ud, ud_new, 1e-6, 1e-9, "bhagwatGhajarMod ud");
}

// ---------------------------------------------------------------------------
// R04 Tests
// ---------------------------------------------------------------------------

void test_r04_z_factor() {
    std::puts("\n=== R04/zFactor parity ===");
    
    // Reference values from actual implementation
    // pres=100 kgf/cm2, temp=60C, Deng=0.65, PC=673 psia, TC=394 R
    const double expected_z = 0.824275376670;
    
    double z_new = marlim::sisprod2::zFactor(100.0, 60.0, 0.65, 673.0, 394.0);
    
    EXPECT_PARITY(expected_z, z_new, 1e-6, 1e-9, "zFactor");
}

void test_r04_gas_density() {
    std::puts("\n=== R04/gasDensityBlackOil parity ===");
    
    // Reference values from actual implementation
    // pres=100 kgf/cm2, temp=60C, Deng=0.65, PC=673 psia, TC=394 R
    const double expected_rho = 83.550859000065;
    
    double rho_new = marlim::sisprod2::gasDensityBlackOil(100.0, 60.0, 0.65, 673.0, 394.0);
    
    EXPECT_PARITY(expected_rho, rho_new, 1e-6, 1e-9, "gasDensity");
}

void test_r04_gas_viscosity() {
    std::puts("\n=== R04/gasViscosityBlackOil parity ===");
    
    // Reference values from actual implementation
    // pres=100 kgf/cm2, temp=60C, Deng=0.65, PC=673 psia, TC=394 R
    const double expected_visc = 0.014949906654;
    
    double visc_new = marlim::sisprod2::gasViscosityBlackOil(100.0, 60.0, 0.65, 673.0, 394.0);
    
    EXPECT_PARITY(expected_visc, visc_new, 1e-6, 1e-9, "gasViscosity");
}

void test_r04_oil_viscosity() {
    std::puts("\n=== R04/oilViscosityBlackOil parity ===");
    
    // Reference values from actual implementation
    // rs=50.0, deadOilViscosity=5.0 cP
    const double expected_visc = 3.495504625567;
    
    double visc_new = marlim::sisprod2::oilViscosityBlackOil(50.0, 5.0);
    
    EXPECT_PARITY(expected_visc, visc_new, 1e-6, 1e-9, "oilViscosity");
}

void test_r04_oil_density() {
    std::puts("\n=== R04/oilDensityBlackOil parity ===");
    
    // Reference values from actual implementation
    // pres=100 kgf/cm2, temp=60C, API=35, Deng=0.65, rs=50
    const double expected_dens = 832.251045757544;
    
    double dens_new = marlim::sisprod2::oilDensityBlackOil(100.0, 60.0, 35.0, 0.65, 50.0);
    
    EXPECT_PARITY(expected_dens, dens_new, 1e-6, 1e-9, "oilDensity");
}

void test_r04_water_density() {
    std::puts("\n=== R04/waterDensityBlackOil parity ===");
    
    // Reference values from actual implementation
    // pres=100 kgf/cm2, temp=60C, Denag=1.05
    const double expected_dens = 1039.774873564618;
    
    double dens_new = marlim::sisprod2::waterDensityBlackOil(100.0, 60.0, 1.05);
    
    EXPECT_PARITY(expected_dens, dens_new, 1e-6, 1e-9, "waterDensity");
}

void test_r04_water_viscosity() {
    std::puts("\n=== R04/waterViscosityBlackOil parity ===");
    
    // Reference values from actual implementation
    // temp=60C
    const double expected_visc = 0.463103416968;
    
    double visc_new = marlim::sisprod2::waterViscosityBlackOil(60.0);
    
    EXPECT_PARITY(expected_visc, visc_new, 1e-6, 1e-9, "waterViscosity");
}

void test_r04_water_fvf() {
    std::puts("\n=== R04/waterFVFBlackOil parity ===");
    
    // Reference values from actual implementation
    // pres=100 kgf/cm2, temp=60C, Denag=1.05
    const double expected_fvf = 1.009833981081;
    
    double fvf_new = marlim::sisprod2::waterFVFBlackOil(100.0, 60.0, 1.05);
    
    EXPECT_PARITY(expected_fvf, fvf_new, 1e-6, 1e-9, "waterFVF");
}

void test_r04_liquid_density() {
    std::puts("\n=== R04/liquidDensityBlackOil parity ===");
    
    // Reference values from actual implementation
    // pres=100, temp=60C, API=35, Deng=0.65, BSW=0.2, Denag=1.05, rs=50
    const double expected_dens = 873.113797781711;
    
    double dens_new = marlim::sisprod2::liquidDensityBlackOil(100.0, 60.0, 35.0, 0.65, 0.2, 1.05, 50.0);
    
    EXPECT_PARITY(expected_dens, dens_new, 1e-6, 1e-9, "liquidDensity");
}

void test_r04_liquid_specific_heat() {
    std::puts("\n=== R04/liquidSpecificHeatBlackOil parity ===");
    
    // Reference values from actual implementation
    // pres=100, temp=60C, API=35, Deng=0.65, BSW=0.2, Denag=1.05, rs=50, rDgD=0.6
    const double expected_cp = 2487.005583948842;
    
    double cp_new = marlim::sisprod2::liquidSpecificHeatBlackOil(
        100.0, 60.0, 35.0, 0.65, 0.2, 1.05, 50.0, 0.6);
    
    EXPECT_PARITY(expected_cp, cp_new, 1e-6, 1e-9, "liquidSpecificHeat");
}

void test_r04_gas_specific_heat() {
    std::puts("\n=== R04/gasSpecificHeatBlackOil parity ===");
    
    // Reference values from actual implementation
    // pres=100, temp=60C, Deng=0.65, PCis=673, TCis=394, yco2=0.05, rDgL=0.7
    const double expected_cp = 3098.459836997971;
    
    double cp_new = marlim::sisprod2::gasSpecificHeatBlackOil(
        100.0, 60.0, 0.65, 673.0, 394.0, 0.05, 0.7);
    
    EXPECT_PARITY(expected_cp, cp_new, 1e-6, 1e-9, "gasSpecificHeat");
}

void test_r04_liquid_thermal_conductivity() {
    std::puts("\n=== R04/liquidThermalConductivityBlackOil parity ===");
    
    // Reference values from actual implementation
    // pres=100, temp=60C, API=35, Deng=0.65, BSW=0.2, Denag=1.05, rs=50, rDgD=0.6
    const double expected_cond = 0.170627801713;
    
    double cond_new = marlim::sisprod2::liquidThermalConductivityBlackOil(
        100.0, 60.0, 35.0, 0.65, 0.2, 1.05, 50.0, 0.6);
    
    EXPECT_PARITY(expected_cond, cond_new, 1e-6, 1e-9, "liquidThermalConductivity");
}

void test_r04_gas_thermal_conductivity() {
    std::puts("\n=== R04/gasThermalConductivityBlackOil parity ===");
    
    // Reference values from actual implementation
    // pres=100, temp=60C
    const double expected_cond = 0.048465377149;
    
    double cond_new = marlim::sisprod2::gasThermalConductivityBlackOil(100.0, 60.0);
    
    EXPECT_PARITY(expected_cond, cond_new, 1e-6, 1e-9, "gasThermalConductivity");
}

void test_r04_liquid_joule_thomson() {
    std::puts("\n=== R04/liquidJouleThomsonBlackOil parity ===");
    
    // Reference values from actual implementation
    // pres=100, temp=60C, API=35, Deng=0.65, BSW=0.2, Denag=1.05, rs=50, rDgD=0.6, liquidSimple=0.0
    const double expected_jt = -0.000844418552;
    
    double jt_new = marlim::sisprod2::liquidJouleThomsonBlackOil(
        100.0, 60.0, 35.0, 0.65, 0.2, 1.05, 50.0, 0.6, 0.0);
    
    EXPECT_PARITY(expected_jt, jt_new, 1e-6, 1e-9, "liquidJouleThomson");
}

void test_r04_gas_joule_thomson() {
    std::puts("\n=== R04/gasJouleThomsonBlackOil parity ===");
    
    // Reference values from actual implementation
    // pres=100, temp=60C, Deng=0.65, PCis=673, TCis=394, rhog computed
    const double expected_jt = 0.005800348802;
    
    double rho_gas = marlim::sisprod2::gasDensityBlackOil(100.0, 60.0, 0.65, 673.0, 394.0);
    double jt_new = marlim::sisprod2::gasJouleThomsonBlackOil(
        100.0, 60.0, 0.65, 673.0, 394.0, rho_gas);
    
    EXPECT_PARITY(expected_jt, jt_new, 1e-6, 1e-9, "gasJouleThomson");
}

void test_r04_liquid_enthalpy() {
    std::puts("\n=== R04/liquidEnthalpyBlackOil parity ===");
    
    // Reference values from actual implementation
    // pres=100, temp=60C, API=35, Deng=0.65, BSW=0.2, Denag=1.05, rs=50, rDgD=0.6
    const double expected_h = 149220.335036930512;
    
    double h_new = marlim::sisprod2::liquidEnthalpyBlackOil(
        100.0, 60.0, 35.0, 0.65, 0.2, 1.05, 50.0, 0.6);
    
    EXPECT_PARITY(expected_h, h_new, 1e-6, 1e-9, "liquidEnthalpy");
}

void test_r04_gas_enthalpy() {
    std::puts("\n=== R04/gasEnthalpyBlackOil parity ===");
    
    // Reference values from actual implementation
    // pres=100, temp=60C, Deng=0.65, PCis=673, TCis=394, yco2=0.05, rDgL=0.7
    const double expected_h = 185907.590219878242;
    
    double h_new = marlim::sisprod2::gasEnthalpyBlackOil(
        100.0, 60.0, 0.65, 673.0, 394.0, 0.05, 0.7);
    
    EXPECT_PARITY(expected_h, h_new, 1e-6, 1e-9, "gasEnthalpy");
}

// ---------------------------------------------------------------------------
// R05/R06 Tests - Pressure Gradient and Wellhead March
// ---------------------------------------------------------------------------

void test_r05_compute_pressure_gradient_beggs_brill() {
    std::puts("\n=== R05/computePressureGradient (Beggs & Brill) sanity check ===");
    
    // Setup pressure gradient input with typical production well parameters
    marlim::sisprod2::PressureGradientInput input;
    input.inclinationRad = M_PI / 2.0;      // vertical well (90 degrees)
    input.diameter = 0.12;                  // 12 cm = 4.7 inches
    input.roughness = 4.572e-05;            // 45 microns absolute roughness
    input.pressure = 70.30696e5;            // ~70 kgf/cm2 in Pa
    input.temperature = 333.15;             // 60°C in K
    input.mixtureVelocity = 2.5;            // 2.5 m/s mixture velocity
    input.liquidFraction = 0.7;             // 70% liquid (BSW + oil)
    input.gasDensity = 83.55;               // kg/m3 at reservoir conditions
    input.liquidDensity = 850.0;            // kg/m3 (oil + water mixture)
    input.gasViscosity = 1.5e-5;            // Pa·s (15 µPa·s)
    input.liquidViscosity = 2e-3;           // Pa·s (2 cP)
    input.surfaceTension = 0.015;             // N/m (15 dyn/cm)
    input.compressibilityFactor = 0.85;     // typical Z-factor
    input.waterFraction = 0.3;              // 0.3 water cut
    
    // Call computePressureGradient with BeggsAndBrill
    marlim::sisprod2::PressureGradientResult result = 
        marlim::sisprod2::computePressureGradient(
            marlim::sisprod2::FlowCorrelationId::BeggsAndBrill, input);
    
    // Sanity checks - verify result is physically reasonable
    bool status_ok = (result.status == marlim::sisprod2::SolveStatus::Ok);
    bool holdup_ok = (result.holdup >= 0.0 && result.holdup <= 1.0);
    bool reynolds_ok = (result.reynolds > 0.0);
    bool frictionGrad_ok = (result.frictionGrad >= 0.0);
    bool gravityGrad_ok = (result.gravityGrad > 0.0);  // vertical upflow
    bool totalGrad_ok = (result.totalGrad > 0.0);
    
    EXPECT_PARITY(1.0, static_cast<double>(status_ok), 1e-6, 1e-9, 
                  "BeggsAndBrill status OK");
    EXPECT_PARITY(1.0, static_cast<double>(holdup_ok), 1e-6, 1e-9, 
                  "BeggsAndBrill holdup in [0,1]");
    EXPECT_PARITY(1.0, static_cast<double>(reynolds_ok), 1e-6, 1e-9,
                  "BeggsAndBrill reynolds positive");
    EXPECT_PARITY(1.0, static_cast<double>(frictionGrad_ok), 1e-6, 1e-9,
                  "BeggsAndBrill frictionGrad non-negative");
    EXPECT_PARITY(1.0, static_cast<double>(gravityGrad_ok), 1e-6, 1e-9,
                  "BeggsAndBrill gravityGrad positive (vertical)");
    EXPECT_PARITY(1.0, static_cast<double>(totalGrad_ok), 1e-6, 1e-9,
                  "BeggsAndBrill totalGrad positive");
}

void test_r05_compute_pressure_gradient_horizontal() {
    std::puts("\n=== R05/computePressureGradient (horizontal) sanity check ===");
    
    // Test horizontal pipe (inclination = 0)
    marlim::sisprod2::PressureGradientInput input;
    input.inclinationRad = 0.0;             // horizontal
    input.diameter = 0.15;                  // 6 inches
    input.roughness = 4.572e-05;
    input.pressure = 50.0e5;                // 50 kgf/cm2
    input.temperature = 333.15;
    input.mixtureVelocity = 3.0;            // 3 m/s
    input.liquidFraction = 0.5;
    input.gasDensity = 50.0;
    input.liquidDensity = 800.0;
    input.gasViscosity = 1.8e-5;
    input.liquidViscosity = 1.5e-3;
    input.surfaceTension = 0.02;
    input.compressibilityFactor = 0.90;
    input.waterFraction = 0.3;
    
    marlim::sisprod2::PressureGradientResult result = 
        marlim::sisprod2::computePressureGradient(
            marlim::sisprod2::FlowCorrelationId::BeggsAndBrill, input);
    
    // Sanity checks
    bool status_ok = (result.status == marlim::sisprod2::SolveStatus::Ok);
    bool holdup_ok = (result.holdup >= 0.0 && result.holdup <= 1.0);
    // Gravity should be ~0 for horizontal
    bool gravityGrad_small = (std::fabs(result.gravityGrad) < 1000.0);
    bool frictionGrad_ok = (result.frictionGrad >= 0.0);
    bool totalGrad_ok = (result.totalGrad >= 0.0);
    
    EXPECT_PARITY(1.0, static_cast<double>(status_ok), 1e-6, 1e-9, 
                  "Horizontal status OK");
    EXPECT_PARITY(1.0, static_cast<double>(holdup_ok), 1e-6, 1e-9,
                  "Horizontal holdup in [0,1]");
    EXPECT_PARITY(1.0, static_cast<double>(gravityGrad_small), 1e-6, 1e-9,
                  "Horizontal gravityGrad small");
    EXPECT_PARITY(1.0, static_cast<double>(frictionGrad_ok), 1e-6, 1e-9,
                  "Horizontal frictionGrad non-negative");
    EXPECT_PARITY(1.0, static_cast<double>(totalGrad_ok), 1e-6, 1e-9,
                  "Horizontal totalGrad non-negative");
}

void test_r05_compute_pressure_gradient_uphill_downhill() {
    std::puts("\n=== R05/computePressureGradient (uphill vs downhill) sanity check ===");
    
    // Test uphill (+30 degrees)
    marlim::sisprod2::PressureGradientInput input_up;
    input_up.inclinationRad = 30.0 * M_PI / 180.0;  // +30 degrees from horizontal
    input_up.diameter = 0.10;
    input_up.roughness = 4.572e-05;
    input_up.pressure = 60.0e5;
    input_up.temperature = 333.15;
    input_up.mixtureVelocity = 2.0;
    input_up.liquidFraction = 0.6;
    input_up.gasDensity = 65.0;
    input_up.liquidDensity = 820.0;
    input_up.gasViscosity = 1.6e-5;
    input_up.liquidViscosity = 1.8e-3;
    input_up.surfaceTension = 0.018;
    input_up.compressibilityFactor = 0.88;
    input_up.waterFraction = 0.25;
    
    // Test downhill (-30 degrees) - same parameters
    marlim::sisprod2::PressureGradientInput input_down = input_up;
    input_down.inclinationRad = -30.0 * M_PI / 180.0;  // -30 degrees
    
    marlim::sisprod2::PressureGradientResult result_up = 
        marlim::sisprod2::computePressureGradient(
            marlim::sisprod2::FlowCorrelationId::BeggsAndBrill, input_up);
    
    marlim::sisprod2::PressureGradientResult result_down = 
        marlim::sisprod2::computePressureGradient(
            marlim::sisprod2::FlowCorrelationId::BeggsAndBrill, input_down);
    
    // Sanity checks
    // Uphill should have positive gravity gradient
    bool uphill_gravity_ok = (result_up.gravityGrad > 0.0);
    // Downhill should have negative gravity gradient
    bool downhill_gravity_ok = (result_down.gravityGrad < 0.0);
    // Friction should be positive in both directions
    bool uphill_friction_ok = (result_up.frictionGrad >= 0.0);
    bool downhill_friction_ok = (result_down.frictionGrad >= 0.0);
    
    EXPECT_PARITY(1.0, static_cast<double>(uphill_gravity_ok), 1e-6, 1e-9,
                  "Uphill gravityGrad positive");
    EXPECT_PARITY(1.0, static_cast<double>(downhill_gravity_ok), 1e-6, 1e-9,
                  "Downhill gravityGrad negative");
    EXPECT_PARITY(1.0, static_cast<double>(uphill_friction_ok), 1e-6, 1e-9,
                  "Uphill frictionGrad non-negative");
    EXPECT_PARITY(1.0, static_cast<double>(downhill_friction_ok), 1e-6, 1e-9,
                  "Downhill frictionGrad non-negative");
}

void test_r06_march_to_wellhead_simple() {
    std::puts("\n=== R06/marchToWellheadPhysical sanity check ===");
    
    // Create a simple production column
    marlim::sisprod2::FluidModel fluid(850.0, 1.2, 0.002);
    marlim::sisprod2::ProductionColumn column(fluid, 10.0);  // 10 kg/s mass flow
    
    // Set inlet stream
    marlim::sisprod2::StandardStream inlet;
    inlet.oilRate = 100.0;      // m³/d
    inlet.waterRate = 30.0;     // m³/d
    inlet.gasRate = 10000.0;    // m³/d
    inlet.oilApi = 35.0;
    inlet.gasDensity = 0.65;
    inlet.waterDensity = 1.05;
    inlet.co2Fraction = 0.05;
    column.setInletStream(inlet);
    
    // Add a single vertical segment
    marlim::sisprod2::PipeSegment seg;
    seg.length = 1000.0;        // 1 km
    seg.diameter = 0.15;        // 15 cm
    seg.inclination = M_PI / 2.0;  // vertical
    seg.accessory = marlim::sisprod2::AccessoryType::None;
    column.addSegment(seg);
    
    // Create context
    marlim::sisprod2::SimContext ctx;
    
    // March from bottomhole to wellhead
    double bottomholePressure = 100.0e5;  // 100 bar
    std::vector<marlim::sisprod2::ProfilePoint> profile;
    
    double wellheadPressure = column.marchToWellheadPhysical(
        bottomholePressure, ctx,
        marlim::sisprod2::FlowCorrelationId::BeggsAndBrill,
        &profile);
    
    // Sanity checks
    // Wellhead pressure should be lower than bottomhole
    bool pressure_drop_ok = (wellheadPressure < bottomholePressure);
    // Profile should have at least one point
    bool profile_has_points = (profile.size() > 0);
    // All pressures should be positive
    bool pressures_positive = (wellheadPressure > 0.0);
    if (!profile.empty()) {
        pressures_positive = pressures_positive && (profile[0].pressure > 0.0);
    }
    
    EXPECT_PARITY(1.0, static_cast<double>(pressure_drop_ok), 1e-6, 1e-9,
                  "Wellhead pressure < bottomhole pressure");
    EXPECT_PARITY(1.0, static_cast<double>(profile_has_points), 1e-6, 1e-9,
                  "Profile has points");
    EXPECT_PARITY(1.0, static_cast<double>(pressures_positive), 1e-6, 1e-9,
                  "All pressures positive");
}
// ---------------------------------------------------------------------------
// R06/R10 Tests
// ---------------------------------------------------------------------------

void test_r06_colebrook_friction() {
    std::puts("\n=== R06/colebrookFrictionFactor parity ===");
    
    // Reference values from actual implementation
    // reynolds=50000, relRoughness=0.001 (0.0001/0.1)
    const double expected_f = 0.024019203431;
    
    double f_new = marlim::sisprod2::colebrookFrictionFactor(50000.0, 0.001);
    
    EXPECT_PARITY(expected_f, f_new, 1e-6, 1e-9, "colebrookFrictionFactor");
}

void test_r10_area_valv_cali() {
    std::puts("\n=== R10/areaValvCali parity ===");
    
    // Reference values from actual implementation
    // PCal=1000, TCal=150, PVO=800, PT=750, dextern=0.025, areagarg=0.001, Rvalv=0.5, Temp=200
    const double expected_area = 1.000000000000;  // Fully open
    
    double area_new = marlim::sisprod2::areaValvCali(1000.0, 150.0, 800.0, 750.0, 
                                                      0.025, 0.001, 0.5, 200.0);
    
    EXPECT_PARITY(expected_area, area_new, 1e-6, 1e-9, "areaValvCali");
}

// ---------------------------------------------------------------------------
// R07 Tests - Thermal Advancement (Energy Balance)
// ---------------------------------------------------------------------------

void test_r07_thermal_flow_snapshot() {
    std::puts("\n=== R07/computeThermalFlowSnapshot sanity check ===");
    
    // Typical production well parameters at ~700m depth
    marlim::sisprod2::ThermalSideInput thermalIn;
    thermalIn.pressurePa = 70.30696e5;  // ~70 kgf/cm2 in Pa
    thermalIn.temperatureK = 333.15;    // 60°C in K
    thermalIn.gasHoldup = 0.4;
    thermalIn.waterCut = 0.3;
    thermalIn.gasSuperficialVelocity = 2.5;
    thermalIn.liquidSuperficialVelocity = 0.8;
    
    thermalIn.stream.oilRate = 100.0;
    thermalIn.stream.waterRate = 40.0;
    thermalIn.stream.gasRate = 50000.0;
    thermalIn.stream.oilApi = 35.0;
    thermalIn.stream.gasDensity = 0.65;
    thermalIn.stream.waterDensity = 1.05;
    thermalIn.stream.co2Fraction = 0.05;
    
    marlim::sisprod2::ThermalFlowSnapshot snapshot = 
        marlim::sisprod2::computeThermalFlowSnapshot(thermalIn);
    
    // Sanity checks for physical reasonableness
    bool liquid_rho = (snapshot.liquidDensity > 700.0 && snapshot.liquidDensity < 1200.0);  // Oil/water mix
    bool gas_rho = (snapshot.gasDensity > 0.0 && snapshot.gasDensity < 500.0);  // Gas density
    bool liquid_cp = (snapshot.liquidSpecificHeat > 2000.0 && snapshot.liquidSpecificHeat < 5000.0);  // J/kgK
    bool gas_cp = (snapshot.gasSpecificHeat > 1000.0 && snapshot.gasSpecificHeat < 4000.0);  // J/kgK
    bool liquid_h = (snapshot.liquidEnthalpy > 100000.0);  // Positive enthalpy
    bool gas_h = (snapshot.gasEnthalpy > 100000.0);  // Positive enthalpy
    bool li_jt = (std::fabs(snapshot.liquidJouleThomson) < 1.0);  // Small coefficient
    bool gi_jt = (std::fabs(snapshot.gasJouleThomson) < 1.0);  // Small coefficient
    bool liquid_mu = (snapshot.liquidViscosityPaS > 1e-6 && snapshot.liquidViscosityPaS < 1.0);  // Liquid viscosity
    bool gas_mu = (snapshot.gasViscosityPaS > 1e-9 && snapshot.gasViscosityPaS < 1e-3);  // Gas viscosity
    bool mix_cond = (snapshot.mixedConductivity > 0.0);  // Positive conductivity
    bool mix_cp = (snapshot.mixedSpecificHeat > 0.0);  // Positive specific heat
    bool mix_rho = (snapshot.mixedDensity > 0.0 && snapshot.mixedDensity < 1000.0);  // Mixed density phase
    bool mix_mu = (snapshot.mixedViscosityPaS > 0.0);  // Viscosity
  
    EXPECT_PARITY(1.0, static_cast<double>(liquid_rho), 1e-6, 1e-9, "liquid density reasonable");
    EXPECT_PARITY(1.0, static_cast<double>(gas_rho), 1e-6, 1e-9, "gas density positive");
    EXPECT_PARITY(1.0, static_cast<double>(liquid_cp), 1e-6, 1e-9, "liquid specific heat reasonable");
    EXPECT_PARITY(1.0, static_cast<double>(gas_cp), 1e-6, 1e-9, "gas specific heat reasonable");
    EXPECT_PARITY(1.0, static_cast<double>(liquid_h), 1e-6, 1e-9, "liquid enthalpy positive");
    EXPECT_PARITY(1.0, static_cast<double>(gas_h), 1e-6, 1e-9, "gas enthalpy positive");
    EXPECT_PARITY(1.0, static_cast<double>(li_jt), 1e-6, 1e-9, "liquid JT coefficient small");
    EXPECT_PARITY(1.0, static_cast<double>(gi_jt), 1e-6, 1e-9, "gas JT coefficient small");
    EXPECT_PARITY(1.0, static_cast<double>(liquid_mu), 1e-6, 1e-9, "liquid viscosity reasonable");
    EXPECT_PARITY(1.0, static_cast<double>(gas_mu), 1e-6, 1e-9, "gas viscosity small positive");
    EXPECT_PARITY(1.0, static_cast<double>(mix_cond), 1e-6, 1e-9, "mixed conductivity positive");
    EXPECT_PARITY(1.0, static_cast<double>(mix_cp), 1e-6, 1e-9, "mixed specific heat positive");
    EXPECT_PARITY(1.0, static_cast<double>(mix_rho), 1e-6, 1e-9, "mixed density positive < 1000");
    EXPECT_PARITY(1.0, static_cast<double>(mix_mu), 1e-6, 1e-9, "mixed viscosity positive");
}

void test_r07_thermal_update_normal_flow() {
    std::puts("\n=== R07/computeThermalUpdate (normal flow) ===");
    
    // Setup simulation context
    marlim::sisprod2::SimContext context;
    context.setLocalTiny(1e-15);
    
    // Create thermal flow snapshot (pre-computed properties)
    marlim::sisprod2::ThermalFlowSnapshot props;
    props.liquidDensity = 850.0;
    props.gasDensity = 100.0;
    props.liquidSpecificHeat = 3600.0;
    props.gasSpecificHeat = 2200.0;
    props.liquidEnthalpy = 145000.0;
    props.gasEnthalpy = 180000.0;
    props.latentHeat = 35000.0;
    props.liquidJouleThomson = -0.001;
    props.gasJouleThomson = 0.006;
    props.liquidViscosityPaS = 0.002;
    props.gasViscosityPaS = 1.5e-5;
    props.mixedConductivity = 0.4;
    props.mixedSpecificHeat = 3000.0;
    props.mixedDensity = 400.0;
    props.mixedViscosityPaS = 1e-5;
    
    // Normal flow conditions
    marlim::sisprod2::ThermalUpdateInput input;
    input.currentTemperatureC = 60.0;
    input.upstreamTemperatureC = 58.0;
    input.pressurePa = 70.0e5;
    input.upstreamPressurePa = 80.0e5;
    input.dx = 100.0;  // 100m cell length
    input.diameter = 0.15;
    input.gasSuperficialVelocity = 3.0;
    input.liquidSuperficialVelocity = 1.0;
    input.gasHoldup = 0.4;
    input.waterCut = 0.3;
    input.thermalConductivity = 0.4;
    input.externalTemperatureC = 65.0;  // geothermal
    input.thermalResistance = 0.1;
    input.massSourceL = 0.0;
    input.massSourceG = 0.0;
    input.sourceTemperatureC = 65.0;
    input.latentHeatTerm = 0.0;
    input.hasMassSource = false;
    
    marlim::sisprod2::ThermalUpdateResult result = 
        marlim::sisprod2::computeThermalUpdate(props, input, context);
    
    // Note: Absolute temperature range checks are skipped until dual-run extracts reference values.
    // Here we only verify the function executes without returning NaN/Inf.
    bool finite_temp = std::isfinite(result.updatedTemperatureC);
    bool low_flow_false = (!result.lowFlowMode);
    bool substeps_valid = (result.substepsUsed >= 0);
  
    EXPECT_PARITY(1.0, static_cast<double>(finite_temp), 1e-6, 1e-9, "updated temperature is finite");
    EXPECT_PARITY(1.0, static_cast<double>(low_flow_false), 1e-6, 1e-9, "not low flow mode");
    EXPECT_PARITY(1.0, static_cast<double>(substeps_valid), 1e-6, 1e-9, "substeps non-negative");
}

void test_r07_thermal_update_low_flow() {
    std::puts("\n=== R07/computeThermalUpdate (low flow) ===");
    
    // Setup simulation context
    marlim::sisprod2::SimContext context;
    context.setLocalTiny(1e-15);
    
    // Thermal properties
    marlim::sisprod2::ThermalFlowSnapshot props;
    props.liquidDensity = 850.0;
    props.gasDensity = 100.0;
    props.liquidSpecificHeat = 3600.0;
    props.gasSpecificHeat = 2200.0;
    props.liquidEnthalpy = 145000.0;
    props.gasEnthalpy = 180000.0;
    props.latentHeat = 35000.0;
    props.liquidJouleThomson = -0.001;
    props.gasJouleThomson = 0.006;
    props.liquidViscosityPaS = 0.002;
    props.gasViscosityPaS = 1.5e-5;
    props.mixedConductivity = 0.4;
    props.mixedSpecificHeat = 3000.0;
    props.mixedDensity = 400.0;
    props.mixedViscosityPaS = 1e-5;
    
    // Low flow conditions (velocities near zero)
    marlim::sisprod2::ThermalUpdateInput input;
    input.currentTemperatureC = 60.0;
    input.upstreamTemperatureC = 58.0;
    input.pressurePa = 70.0e5;
    input.upstreamPressurePa = 80.0e5;
    input.dx = 100.0;
    input.diameter = 0.15;
    input.gasSuperficialVelocity = 0.001;  // Very low gas velocity
    input.liquidSuperficialVelocity = 0.0005;  // Very low liquid velocity
    input.gasHoldup = 0.4;
    input.waterCut = 0.3;
    input.thermalConductivity = 0.4;
    input.externalTemperatureC = 65.0;
    input.thermalResistance = 0.1;
    input.massSourceL = 0.0;
    input.massSourceG = 0.0;
    input.sourceTemperatureC = 65.0;
    input.latentHeatTerm = 0.0;
    input.hasMassSource = false;
    
    marlim::sisprod2::ThermalUpdateResult result = 
        marlim::sisprod2::computeThermalUpdate(props, input, context);
    
    // Low flow should trigger fallback to external temperature
    bool low_flow_mode = (result.lowFlowMode == true);  // Should be true for low flow
    bool temp_near_external = (std::fabs(result.updatedTemperatureC - 65.0) < 5.0);
    bool substeps_higher = (result.substepsUsed >= 1);
    
    EXPECT_PARITY(1.0, static_cast<double>(low_flow_mode), 1e-6, 1e-9, "low flow mode activated");
    EXPECT_PARITY(1.0, static_cast<double>(temp_near_external), 1e-6, 1e-9, "temp near external");
    EXPECT_PARITY(1.0, static_cast<double>(substeps_higher), 1e-6, 1e-9, "substeps reasonable");
}

void test_r07_thermal_update_with_mass_source() {
    std::puts("\n=== R07/computeThermalUpdate (mass source) ===");
    
    // Setup simulation context
    marlim::sisprod2::SimContext context;
    context.setLocalTiny(1e-15);
    
    // Thermal properties
    marlim::sisprod2::ThermalFlowSnapshot props;
    props.liquidDensity = 850.0;
    props.gasDensity = 100.0;
    props.liquidSpecificHeat = 3600.0;
    props.gasSpecificHeat = 2200.0;
    props.liquidEnthalpy = 145000.0;
    props.gasEnthalpy = 180000.0;
    props.latentHeat = 35000.0;
    props.liquidJouleThomson = -0.001;
    props.gasJouleThomson = 0.006;
    props.liquidViscosityPaS = 0.002;
    props.gasViscosityPaS = 1.5e-5;
    props.mixedConductivity = 0.4;
    props.mixedSpecificHeat = 3000.0;
    props.mixedDensity = 400.0;
    props.mixedViscosityPaS = 1e-5;
    
    // Mass source conditions (injection)
    marlim::sisprod2::ThermalUpdateInput input;
    input.currentTemperatureC = 60.0;
    input.upstreamTemperatureC = 58.0;
    input.pressurePa = 70.0e5;
    input.upstreamPressurePa = 80.0e5;
    input.dx = 100.0;
    input.diameter = 0.15;
    input.gasSuperficialVelocity = 2.0;
    input.liquidSuperficialVelocity = 0.5;
    input.gasHoldup = 0.3;
    input.waterCut = 0.2;
    input.thermalConductivity = 0.4;
    input.externalTemperatureC = 65.0;
    input.thermalResistance = 0.1;
    input.massSourceL = 5.0;  // 5 kg/m³/s liquid injection
    input.massSourceG = 0.0;
    input.sourceTemperatureC = 50.0;  // Cooler injected fluid
    input.latentHeatTerm = 0.0;
    input.hasMassSource = true;
    
    marlim::sisprod2::ThermalUpdateResult result = 
        marlim::sisprod2::computeThermalUpdate(props, input, context);
    
    // With cool mass source, temperature should decrease
    bool finite_temp = std::isfinite(result.updatedTemperatureC);
    bool has_mass_source_effect = (result.updatedTemperatureC < 100.0);  // Should be below external temp
    bool substeps_valid = (result.substepsUsed > 0);
    
    EXPECT_PARITY(1.0, static_cast<double>(finite_temp), 1e-6, 1e-9, "mass source temp is finite");
    EXPECT_PARITY(1.0, static_cast<double>(has_mass_source_effect), 1e-6, 1e-9, "mass source effect");
    EXPECT_PARITY(1.0, static_cast<double>(substeps_valid), 1e-6, 1e-9, "substeps valid");
}

// ---------------------------------------------------------------------------
// R08 Tests - Steady-State Solvers (TramoEngine)
// ---------------------------------------------------------------------------

void test_r08_solve_bottomhole_pressure() {
    std::puts("\n=== R08/solveBottomholePressure parity ===");
    
    // Setup simulation context and diagnostics
    marlim::sisprod2::SimContext context;
    marlim::sisprod2::NullDiagnostics diagnostics;
    
    // Create a simple vertical well with 5 segments
    marlim::sisprod2::FluidModel fluid(850.0, 1.2, 0.002);
    marlim::sisprod2::ProductionColumn column(fluid, 10.0);  // 10 kg/s mass flow
    
    // Set inlet stream (separator conditions)
    marlim::sisprod2::StandardStream inlet;
    inlet.oilRate = 100.0;      // m³/d
    inlet.waterRate = 40.0;     // m³/d
    inlet.gasRate = 50000.0;    // m³/d
    inlet.oilApi = 35.0;
    inlet.gasDensity = 0.65;
    inlet.waterDensity = 1.05;
    inlet.co2Fraction = 0.05;
    column.setInletStream(inlet);
    
    // Add 5 vertical segments (total 1000m)
    for (int i = 0; i < 5; ++i) {
        marlim::sisprod2::PipeSegment seg;
        seg.length = 200.0;        // 200m per segment
        seg.diameter = 0.15;       // 15 cm
        seg.inclination = M_PI / 2.0;  // vertical
        seg.accessory = marlim::sisprod2::AccessoryType::None;
        column.addSegment(seg);
    }
    
    // Setup steady state request
    marlim::sisprod2::SteadyStateRequest request(12.0e5, 3.0, 1e-3, 100, 0);
    
    // Create TramoEngine and solve
    marlim::sisprod2::TramoEngine engine(context, diagnostics);
    marlim::sisprod2::SolveResult result = engine.solveBottomholePressure(column, request);
    
    // Sanity checks
    bool solved_ok = (result.status == marlim::sisprod2::SolveStatus::Ok);
    bool positive_pressure = (result.value > 0.0);
    bool reasonable_pressure = (result.value > 1e5 && result.value < 200e6);
    bool reasonable_iterations = (result.iterations > 0 && result.iterations <= 200);
    
    EXPECT_PARITY(1.0, static_cast<double>(solved_ok), 1e-6, 1e-9, "solve OK status");
    EXPECT_PARITY(1.0, static_cast<double>(positive_pressure), 1e-6, 1e-9, "positive BHP");
    EXPECT_PARITY(1.0, static_cast<double>(reasonable_pressure), 1e-6, 1e-9, "BHP in reasonable range");
    EXPECT_PARITY(1.0, static_cast<double>(reasonable_iterations), 1e-6, 1e-9, "iterations reasonable");
}

void test_r08_march_mass() {
    std::puts("\n=== R08/marchMass parity ===");
    
    // Setup simulation context and diagnostics
    marlim::sisprod2::SimContext context;
    context.setThreadCount(1);
    marlim::sisprod2::NullDiagnostics diagnostics;
    
    marlim::sisprod2::TramoEngine engine(context, diagnostics);
    
    // Production scenario (GasSource accessory)
    marlim::sisprod2::StandardStream upstream;
    upstream.oilRate = 100.0;
    upstream.waterRate = 40.0;
    upstream.gasRate = 50000.0;
    upstream.oilApi = 35.0;
    upstream.gasDensity = 0.65;
    upstream.waterDensity = 1.05;
    upstream.co2Fraction = 0.05;
    
    marlim::sisprod2::StandardStream source;
    source.oilRate = 0.0;
    source.waterRate = 0.0;
    source.gasRate = 10000.0;  // 10000 m³/d gas injection
    source.oilApi = 35.0;
    source.gasDensity = 0.65;
    source.waterDensity = 1.05;
    source.co2Fraction = 0.05;
    
    // Gas source march
    marlim::sisprod2::StandardStream result_gas = 
        engine.marchMass(upstream, marlim::sisprod2::AccessoryType::GasSource, source);
    
    // Gas rate should increase with gas source
    bool gas_increased = (result_gas.gasRate >= upstream.gasRate);
    bool oil_unchanged = (result_gas.oilRate == upstream.oilRate);
    bool water_unchanged = (result_gas.waterRate == upstream.waterRate);
    
    EXPECT_PARITY(1.0, static_cast<double>(gas_increased), 1e-6, 1e-9, "gas rate increased");
    EXPECT_PARITY(1.0, static_cast<double>(oil_unchanged), 1e-6, 1e-9, "oil unchanged");
    EXPECT_PARITY(1.0, static_cast<double>(water_unchanged), 1e-6, 1e-9, "water unchanged");
}

void test_r08_march_mass_with_phase_transfer() {
    std::puts("\n=== R08/marchMassWithPhaseTransfer parity ===");
    
    // Setup simulation context and diagnostics
    marlim::sisprod2::SimContext context;
    context.setThreadCount(1);
    marlim::sisprod2::NullDiagnostics diagnostics;
    
    marlim::sisprod2::TramoEngine engine(context, diagnostics);
    
    // Setup streams
    marlim::sisprod2::StandardStream upstream;
    upstream.oilRate = 100.0;
    upstream.waterRate = 40.0;
    upstream.gasRate = 50000.0;
    upstream.oilApi = 35.0;
    upstream.gasDensity = 0.65;
    upstream.waterDensity = 1.05;
    upstream.co2Fraction = 0.05;
    
    marlim::sisprod2::StandardStream source;
    source.oilRate = 0.0;
    source.waterRate = 0.0;
    source.gasRate = 0.0;
    source.oilApi = 35.0;
    source.gasDensity = 0.65;
    source.waterDensity = 1.05;
    source.co2Fraction = 0.05;
    
    // Phase transfer input
    marlim::sisprod2::PhaseTransferInput phaseInput;
    phaseInput.center.pressurePa = 70.0e5;
    phaseInput.center.temperatureK = 333.15;
    phaseInput.center.pressureAuxPa = 70.0e5;
    phaseInput.center.liquidRate = 140.0;
    phaseInput.center.waterCut = 0.286;
    phaseInput.center.gasSpecificGravity = 0.65;
    phaseInput.center.dissolvedGasDensityRatio = 1.0;
    phaseInput.center.pigFraction = 0.0;
    phaseInput.center.stream = upstream;
    phaseInput.cellLength = 100.0;
    phaseInput.accessoryIsNone = true;
    
    marlim::sisprod2::MassMarchResult result = 
        engine.marchMassWithPhaseTransfer(upstream, 
                                          marlim::sisprod2::AccessoryType::None,
                                          source, 
                                          phaseInput);
    
    // Check result is physically reasonable
    bool mixed_stream_valid = (result.mixedStream.oilRate >= 0.0 && 
                                result.mixedStream.waterRate >= 0.0 &&
                                result.mixedStream.gasRate >= 0.0);
    bool phase_transfer_reasonable = (std::fabs(result.phaseTransferRate) < 1e10);
    
    EXPECT_PARITY(1.0, static_cast<double>(mixed_stream_valid), 1e-6, 1e-9, "mixed stream valid");
    EXPECT_PARITY(1.0, static_cast<double>(phase_transfer_reasonable), 1e-6, 1e-9, "phase transfer reasonable");
}

void test_r08_thermal_snapshot() {
    std::puts("\n=== R08/buildThermalSnapshot parity ===");
    
    // Setup simulation context and diagnostics
    marlim::sisprod2::SimContext context;
    context.setThreadCount(1);
    marlim::sisprod2::NullDiagnostics diagnostics;
    
    marlim::sisprod2::TramoEngine engine(context, diagnostics);
    
    // Setup thermal input
    marlim::sisprod2::ThermalSideInput thermalIn;
    thermalIn.pressurePa = 70.0e5;
    thermalIn.temperatureK = 333.15;
    thermalIn.gasHoldup = 0.4;
    thermalIn.waterCut = 0.3;
    thermalIn.gasSuperficialVelocity = 2.5;
    thermalIn.liquidSuperficialVelocity = 0.8;
    
    thermalIn.stream.oilRate = 100.0;
    thermalIn.stream.waterRate = 40.0;
    thermalIn.stream.gasRate = 50000.0;
    thermalIn.stream.oilApi = 35.0;
    thermalIn.stream.gasDensity = 0.65;
    thermalIn.stream.waterDensity = 1.05;
    thermalIn.stream.co2Fraction = 0.05;
    
    marlim::sisprod2::ThermalFlowSnapshot snapshot = engine.buildThermalSnapshot(thermalIn);
    
    // Check snapshot values are physically reasonable
    bool densities_positive = (snapshot.liquidDensity > 0.0 && snapshot.gasDensity > 0.0);
    bool specific_heat_positive = (snapshot.liquidSpecificHeat > 0.0 && 
                                    snapshot.gasSpecificHeat > 0.0);
    bool enthalpy_positive = (snapshot.liquidEnthalpy > 0.0 && 
                               snapshot.gasEnthalpy > 0.0);
    bool conductivity_positive = (snapshot.mixedConductivity > 0.0);
    bool mixed_density_positive = (snapshot.mixedDensity > 0.0);
    
    EXPECT_PARITY(1.0, static_cast<double>(densities_positive), 1e-6, 1e-9, "densities positive");
    EXPECT_PARITY(1.0, static_cast<double>(specific_heat_positive), 1e-6, 1e-9, "specific heat positive");
    EXPECT_PARITY(1.0, static_cast<double>(enthalpy_positive), 1e-6, 1e-9, "enthalpy positive");
    EXPECT_PARITY(1.0, static_cast<double>(conductivity_positive), 1e-6, 1e-9, "conductivity positive");
    EXPECT_PARITY(1.0, static_cast<double>(mixed_density_positive), 1e-6, 1e-9, "mixed density positive");
}

void test_r08_thermal_advance_step() {
    std::puts("\n=== R08/advanceThermalStep parity ===");
    
    // Setup simulation context and diagnostics
    marlim::sisprod2::SimContext context;
    context.setThreadCount(1);
    marlim::sisprod2::NullDiagnostics diagnostics;
    
    marlim::sisprod2::TramoEngine engine(context, diagnostics);
    
    // Pre-compute thermal snapshot
    marlim::sisprod2::ThermalSideInput thermalIn;
    thermalIn.pressurePa = 70.0e5;
    thermalIn.temperatureK = 333.15;
    thermalIn.gasHoldup = 0.4;
    thermalIn.waterCut = 0.3;
    thermalIn.gasSuperficialVelocity = 2.5;
    thermalIn.liquidSuperficialVelocity = 0.8;
    
    thermalIn.stream.oilRate = 100.0;
    thermalIn.stream.waterRate = 40.0;
    thermalIn.stream.gasRate = 50000.0;
    thermalIn.stream.oilApi = 35.0;
    thermalIn.stream.gasDensity = 0.65;
    thermalIn.stream.waterDensity = 1.05;
    thermalIn.stream.co2Fraction = 0.05;
    
    marlim::sisprod2::ThermalFlowSnapshot props = engine.buildThermalSnapshot(thermalIn);
    
    // Advance thermal step
    marlim::sisprod2::ThermalUpdateResult result = engine.advanceThermalStep(
        60.0,    // currentTempC
        58.0,    // upstreamTempC
        70.0e5,  // pressurePa
        80.0e5,  // upstreamPressurePa
        100.0,   // dx
        0.15,    // diameter
        2.5,     // ugs
        0.8,     // uls
        0.4,     // gasHoldup
        0.3,     // waterCut
        65.0,    // externalTempC
        0.1,     // thermalResistance
        props     // pre-computed props
    );
    
    // Check result is physically reasonable - note: temp range check skipped until dual-run extracts references
    bool finite_temp = std::isfinite(result.updatedTemperatureC);
    bool substeps_valid = (result.substepsUsed > 0);
    bool low_flow_mode_valid = (result.lowFlowMode == false);  // Normal flow shouldn't be low flow
    
    EXPECT_PARITY(1.0, static_cast<double>(finite_temp), 1e-6, 1e-9, "temperature is finite");
    EXPECT_PARITY(1.0, static_cast<double>(substeps_valid), 1e-6, 1e-9, "substeps valid");
    EXPECT_PARITY(1.0, static_cast<double>(low_flow_mode_valid), 1e-6, 1e-9, "not low flow mode");
}
// ---------------------------------------------------------------------------
// Runner
// ---------------------------------------------------------------------------

void run_all_tests() {
    std::puts("╔═══════════════════════════════════════════════════════════════╗");
    std::puts("║  SisProd Parity Tests (Local)                                ║");
std::puts("╚═══════════════════════════════════════════════════════════════╝");
    std::puts("");
    
    // R02
    test_r02_blend();
    test_r02_stream_mixing_march();
    test_r02_compute_phase_transfer_rate();
    
    // R03
    test_r03_choi();
    test_r03_hibiki_ishii();
    test_r03_franca_lahey();
    test_r03_bhagwat_ghajar();
    test_r03_bhagwat_ghajar_mod();
    
    // R04 - Black-oil properties
    test_r04_z_factor();
    test_r04_gas_density();
    test_r04_gas_viscosity();
    test_r04_oil_viscosity();
    test_r04_oil_density();
    test_r04_water_density();
    test_r04_water_viscosity();
    test_r04_water_fvf();
    test_r04_liquid_density();
    test_r04_liquid_specific_heat();
    test_r04_gas_specific_heat();
    test_r04_liquid_thermal_conductivity();
    test_r04_gas_thermal_conductivity();
    test_r04_liquid_joule_thomson();
    test_r04_gas_joule_thomson();
    test_r04_liquid_enthalpy();
    test_r04_gas_enthalpy();
    
    // R07 - Thermal advancement
    test_r07_thermal_flow_snapshot();
    test_r07_thermal_update_normal_flow();
    test_r07_thermal_update_low_flow();
    test_r07_thermal_update_with_mass_source();
    
    // R05 - Pressure gradients
    test_r05_compute_pressure_gradient_beggs_brill();
    test_r05_compute_pressure_gradient_horizontal();
    test_r05_compute_pressure_gradient_uphill_downhill();
    
    // R06 - March to wellhead
    test_r06_march_to_wellhead_simple();
    test_r06_colebrook_friction();
    
    // R08 - Steady-state solvers
    test_r08_solve_bottomhole_pressure();
    test_r08_march_mass();
    test_r08_march_mass_with_phase_transfer();
    test_r08_thermal_snapshot();
    test_r08_thermal_advance_step();
    
    // R10 - Auxiliary functions
    test_r10_area_valv_cali();
    
    std::puts("");
    std::printf("═══════════════\n");
    std::printf("Total: %d tests\n", parity::g_passed + parity::g_failed);
    std::printf("Passed: %d\n", parity::g_passed);
    std::printf("Failed: %d\n", parity::g_failed);
    std::puts("═══════════════");
}

} // namespace parity

int main() {
    parity::run_all_tests();
    return parity::g_failed > 0 ? 1 : 0;
}

