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
// Runner
// ---------------------------------------------------------------------------

void run_all_tests() {
    std::puts("╔═══════════════════════════════════════════════════════════════╗");
    std::puts("║  SisProd Parity Tests (Local)                                ║");
    std::puts("╚═══════════════════════════════════════════════════════════════╝");
    std::puts("");
    
    // R03
    test_r03_choi();
    test_r03_hibiki_ishii();
    test_r03_franca_lahey();
    test_r03_bhagwat_ghajar();
    
    // R04
    test_r04_z_factor();
    test_r04_gas_density();
    
    // R06/R10
    test_r06_colebrook_friction();
    test_r10_area_valv_cali();
    
    std::puts("");
    std::puts("════════════════════════════════════════════════════════════════");
    std::printf("  RESULT: %d passed, %d failed\n", g_passed, g_failed);
    std::puts("════════════════════════════════════════════════════════════════");
    
    if (g_failed > 0) {
        std::puts("  ⚠️  Parity failures detected.");
    }
}

} // namespace parity

int main() {
    parity::run_all_tests();
    return parity::g_failed > 0 ? 1 : 0;
}
