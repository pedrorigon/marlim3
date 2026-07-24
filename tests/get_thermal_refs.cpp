/*
 * get_thermal_refs.cpp - Get reference values for thermal properties
 */

#include "SisProd.h"
#include <cstdio>

int main() {
    using namespace marlim::sisprod2;
    
    std::puts("=== R04 Thermal Properties Reference Values ===\n");
    
    // Standard test conditions
    const double pres = 100.0;   // kgf/cm2
    const double temp = 60.0;    // C
    const double API = 35.0;
    const double Deng = 0.65;
    const double Denag = 1.05;
    const double BSW = 0.2;
    const double rs = 50.0;
    const double rDgD = 0.6;
    const double PCis = 673.0;   // psia
    const double TCis = 394.0;   // R
    const double yco2 = 0.05;
    const double rDgL = 0.7;
    
    // 1. Liquid specific heat
    double cp_liq = liquidSpecificHeatBlackOil(pres, temp, API, Deng, BSW, Denag, rs, rDgD);
    std::printf("liquidSpecificHeatBlackOil: %.12f J/kg/K\n", cp_liq);
    
    // 2. Gas specific heat
    double cp_gas = gasSpecificHeatBlackOil(pres, temp, Deng, PCis, TCis, yco2, rDgL);
    std::printf("gasSpecificHeatBlackOil: %.12f J/kg/K\n", cp_gas);
    
    // 3. Liquid thermal conductivity
    double cond_liq = liquidThermalConductivityBlackOil(pres, temp, API, Deng, BSW, Denag, rs, rDgD);
    std::printf("liquidThermalConductivityBlackOil: %.12f W/m/K\n", cond_liq);
    
    // 4. Gas thermal conductivity
    double cond_gas = gasThermalConductivityBlackOil(pres, temp);
    std::printf("gasThermalConductivityBlackOil: %.12f W/m/K\n", cond_gas);
    
    // 5. Liquid Joule-Thomson
    double jt_liq = liquidJouleThomsonBlackOil(pres, temp, API, Deng, BSW, Denag, rs, rDgD, 0.0);
    std::printf("liquidJouleThomsonBlackOil: %.12f K/MPa\n", jt_liq);
    
    // 6. Gas Joule-Thomson
    double rho_gas = gasDensityBlackOil(pres, temp, Deng, PCis, TCis);
    double jt_gas = gasJouleThomsonBlackOil(pres, temp, Deng, PCis, TCis, rho_gas);
    std::printf("gasJouleThomsonBlackOil: %.12f K/MPa\n", jt_gas);
    
    // 7. Liquid enthalpy
    double h_liq = liquidEnthalpyBlackOil(pres, temp, API, Deng, BSW, Denag, rs, rDgD);
    std::printf("liquidEnthalpyBlackOil: %.12f J/kg\n", h_liq);
    
    // 8. Gas enthalpy
    double h_gas = gasEnthalpyBlackOil(pres, temp, Deng, PCis, TCis, yco2, rDgL);
    std::printf("gasEnthalpyBlackOil: %.12f J/kg\n", h_gas);
    
    return 0;
}
