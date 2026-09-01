#include "SisProdThermal.h"

#include "Leitura.h"
#include "celula3.h"
#include "solver3DPoisson.h"

#include <math.h>

namespace sisprod::thermal {

double interpolateLatentHeat(const ThermalState &state, double pressure, double temperature) {
    int divisionCount = state.input.tabent.npont - 1;
    int pressureIndex = 0.;
    int temperatureIndex = 0.;
    int pressureSearchMarker;
    int temperatureSearchMarker;
    double latentHeat;
    if (pressure < state.latentHeatTable[1][0] || pressure >= state.latentHeatTable[divisionCount + 1][0] || temperature < state.latentHeatTable[0][1] || temperature >= state.latentHeatTable[0][divisionCount + 1])
        latentHeat = 0.;

    else {
        int searchLow, searchMiddle, searchHigh;
        searchLow = 1;
        searchHigh = divisionCount + 1;
        while (searchLow <= searchHigh) {
            searchMiddle = (searchLow + searchHigh) / 2;
            pressureSearchMarker = searchMiddle;
            if (searchMiddle == 1) {
                pressureIndex = searchMiddle;
                break;
            } else if (searchMiddle == divisionCount + 1 && state.latentHeatTable[searchMiddle][0] == pressure) {
                pressureIndex = searchMiddle - 1;
                break;
            }
            if (state.latentHeatTable[searchMiddle][0] > pressure && state.latentHeatTable[searchMiddle - 1][0] <= pressure) {
                pressureIndex = searchMiddle - 1;
                break;
            }
            if (state.latentHeatTable[searchMiddle][0] < pressure)
                searchLow = searchMiddle + 1;
            else
                searchHigh = searchMiddle - 1;
        }
        searchLow = 1;
        searchHigh = divisionCount + 1;
        while (searchLow <= searchHigh) {
            searchMiddle = (searchLow + searchHigh) / 2;
            temperatureSearchMarker = searchMiddle;
            if (searchMiddle == 1) {
                temperatureIndex = searchMiddle;
                break;
            } else if (searchMiddle == divisionCount + 1 && state.latentHeatTable[0][searchMiddle] == temperature) {
                temperatureIndex = searchMiddle - 1;
                break;
            }
            if (state.latentHeatTable[0][searchMiddle] > temperature && state.latentHeatTable[0][searchMiddle - 1] <= temperature) {
                temperatureIndex = searchMiddle - 1;
                break;
            }
            if (state.latentHeatTable[0][searchMiddle - 1] < temperature)
                searchLow = searchMiddle + 1;
            else
                searchHigh = searchMiddle - 1;
        }
        double pressureRatio = (state.latentHeatTable[pressureIndex][0] - pressure) / (state.latentHeatTable[pressureIndex][0] - state.latentHeatTable[pressureIndex + 1][0]);
        double temperatureRatio = (state.latentHeatTable[0][temperatureIndex] - temperature) / (state.latentHeatTable[0][temperatureIndex] - state.latentHeatTable[0][temperatureIndex + 1]);
        double latentHeatAtTemperatureIndex = (1 - pressureRatio) * (state.latentHeatTable[pressureIndex][temperatureIndex]) + pressureRatio * (state.latentHeatTable[pressureIndex + 1][temperatureIndex]);
        double latentHeatAtNextTemperatureIndex = (1 - pressureRatio) * (state.latentHeatTable[pressureIndex][temperatureIndex + 1]) + pressureRatio * (state.latentHeatTable[pressureIndex + 1][temperatureIndex + 1]);
        latentHeat = (1 - temperatureRatio) * latentHeatAtTemperatureIndex + temperatureRatio * latentHeatAtNextTemperatureIndex;
    }
    return latentHeat;
}

double computeMixtureEnthalpy(const ThermalState &state, int cellIndex) {

    double cellLength = state.cells[cellIndex].dx;
    double meanCellLength = 0.5 * (state.cells[cellIndex].dx + state.cells[cellIndex - 1].dx);
    double diameter = state.cells[cellIndex].duto.a;
    double flowArea = 0.25 * M_PI * diameter * diameter;
    double meanVoidFraction = state.cells[cellIndex].alf;
    double betmed = state.cells[cellIndex].bet;
    double previousMeanVoidFraction = state.cells[cellIndex].alfini;
    double betmed0 = state.cells[cellIndex].betini;
    double meanPressure = state.cells[cellIndex].pres;
    double previousMeanPressure = state.cells[cellIndex].presini;
    double meanTemperature = state.cells[cellIndex].temp;
    double leftSuperficialGasVelocity;
    double leftSuperficialLiquidVelocity;
    double rightSuperficialGasVelocity;
    double rightSuperficialLiquidVelocity;

    leftSuperficialGasVelocity = state.cells[cellIndex].QG / flowArea;
    leftSuperficialLiquidVelocity = state.cells[cellIndex].QL / flowArea;
    rightSuperficialGasVelocity = state.cells[cellIndex + 1].QG / flowArea;
    rightSuperficialLiquidVelocity = state.cells[cellIndex + 1].QL / flowArea;

    double betL = state.cells[cellIndex].bet;
    double betR = state.cells[cellIndex].bet;
    if (leftSuperficialGasVelocity > 0.)
        betL = state.cells[cellIndex - 1].bet;
    if (leftSuperficialGasVelocity < 0.)
        betL = state.cells[cellIndex + 1].bet;

    double meanSuperficialGasVelocity = 0.5 * (leftSuperficialGasVelocity + rightSuperficialGasVelocity);
    double meanSuperficialLiquidVelocity = 0.5 * (leftSuperficialLiquidVelocity + rightSuperficialLiquidVelocity);

    double leftPressure = state.cells[cellIndex].presaux;
    double rightPressure = state.cells[cellIndex + 1].presaux;
    if ((state.cells[cellIndex].acsr.tipo == 5 && state.cells[cellIndex].acsr.chk.AreaGarg <= (1e-3 + state.input.master1.razareaativ) * state.cells[cellIndex].duto.area) ||
        (state.cells[cellIndex].acsr.tipo == 4 && state.cells[cellIndex].acsr.bcs.freq > 0) ||
        (state.cells[cellIndex].acsr.tipo == 8 && state.cells[cellIndex].acsr.bvol.freq > 0.) ||
        (state.cells[cellIndex].acsr.tipo == 7 && fabs(state.cells[cellIndex].acsr.delp) > 0.) ||
        (state.cells[cellIndex].acsr.tipo == 17 && state.cells[cellIndex].acsr.multibcs.freq > 0)) {
        rightPressure = state.cells[cellIndex].pres + (state.cells[cellIndex].pres - state.cells[cellIndex].presaux) * 0.5;
    }

    double leftTemperature = state.cells[cellIndex].tempL;
    if (state.cells[cellIndex].VTemper < 0.)
        leftTemperature = state.cells[cellIndex].temp;
    double rightTemperature = state.cells[cellIndex].temp;
    if (state.cells[cellIndex + 1].VTemper < 0.)
        rightTemperature = state.cells[cellIndex].tempR;

    double gasDensity = state.cells[cellIndex].flui.MasEspGas(previousMeanPressure, meanTemperature);
    double liquidDensity = state.cells[cellIndex].flui.MasEspLiq(previousMeanPressure, meanTemperature);
    double rhoc = state.cells[cellIndex].fluicol.MasEspFlu(previousMeanPressure, meanTemperature);
    double gasDensityAtCellPressure = state.cells[cellIndex].flui.MasEspGas(meanPressure, meanTemperature);
    double liquidDensityAtCellPressure = state.cells[cellIndex].flui.MasEspLiq(meanPressure, meanTemperature);
    double gasEnthalpy = state.cells[cellIndex].flui.EntalpGas(previousMeanPressure, meanTemperature);
    double liquidEnthalpy = state.cells[cellIndex].flui.EntalpLiq(previousMeanPressure, meanTemperature);
    double hc = state.cells[cellIndex].fluicol.CalorLiq(previousMeanPressure, meanTemperature) * meanTemperature; // corrigir entalpia

    double leftGasDensity = state.cells[cellIndex].flui.MasEspGas(leftPressure, leftTemperature);
    double leftLiquidDensity = state.cells[cellIndex].flui.MasEspLiq(leftPressure, leftTemperature);
    double rhocL = state.cells[cellIndex].fluicol.MasEspFlu(leftPressure, leftTemperature);
    double leftGasEnthalpy = state.cells[cellIndex].flui.EntalpGas(leftPressure, leftTemperature);
    double leftLiquidEnthalpy = state.cells[cellIndex].flui.EntalpLiq(leftPressure, leftTemperature);
    double hcL = state.cells[cellIndex].fluicol.CalorLiq(leftPressure, leftTemperature) * leftTemperature; // corrigir entalpia

    double rightGasDensity = state.cells[cellIndex].flui.MasEspGas(rightPressure, rightTemperature);
    double rightLiquidDensity = state.cells[cellIndex].flui.MasEspLiq(rightPressure, rightTemperature);
    double rhocR = state.cells[cellIndex].fluicol.MasEspFlu(rightPressure, rightTemperature);
    double rightGasEnthalpy = state.cells[cellIndex].flui.EntalpGas(rightPressure, rightTemperature);
    double rightLiquidEnthalpy = state.cells[cellIndex].flui.EntalpLiq(rightPressure, rightTemperature);
    double hcR = state.cells[cellIndex].fluicol.CalorLiq(rightPressure, rightTemperature) * rightTemperature; // corrigir entalpia

    double previousMixtureInternalEnergy = gasDensity * previousMeanVoidFraction * (gasEnthalpy - previousMeanPressure * 98066.5 / gasDensity) + liquidDensity * (1 - betmed0) * (1. - previousMeanVoidFraction) * (liquidEnthalpy - previousMeanPressure * 98066.5 / liquidDensity) +
                           betmed0 * rhoc * (1. - previousMeanVoidFraction) * (hc - previousMeanPressure * 98066.5 / rhoc);
    double leftEnthalpyFlux = leftGasDensity * leftSuperficialGasVelocity * leftGasEnthalpy + leftLiquidDensity * (1. - betL) * leftSuperficialLiquidVelocity * leftLiquidEnthalpy + rhocL * betL * leftSuperficialLiquidVelocity * hcL;
    double rightEnthalpyFlux = rightGasDensity * rightSuperficialGasVelocity * rightGasEnthalpy + rightLiquidDensity * (1. - betR) * rightSuperficialLiquidVelocity * rightLiquidEnthalpy + rhocR * betR * rightSuperficialLiquidVelocity * hcR;
    double enthalpyFluxDivergence = (rightEnthalpyFlux - leftEnthalpyFlux) / cellLength;
    double hydrostaticTerm = (gasDensityAtCellPressure * meanSuperficialGasVelocity + (1 - betmed) * meanSuperficialLiquidVelocity * liquidDensityAtCellPressure + betmed * meanSuperficialLiquidVelocity) * sin(state.cells[cellIndex].duto.teta) * 9.82;

    double gasMassSourceTerm = 0.;
    double liquidMassSourceTerm = 0.;
    double fontemassC = 0.;
    double sourceTemperature = state.cells[cellIndex].temp;
    double sourceGasEnthalpy;
    double sourceLiquidEnthalpy;
    double hcF = 0.;

    if (state.cells[cellIndex].acsr.tipo == 1) {
        sourceTemperature = state.cells[cellIndex].acsr.injg.temp;
        sourceGasEnthalpy = state.cells[cellIndex].acsr.injg.FluidoPro.EntalpGas(state.cells[cellIndex].pres, sourceTemperature);
        sourceLiquidEnthalpy = 0;
        hcF = 0;
    } else if (state.cells[cellIndex].acsr.tipo == 2) {
        sourceTemperature = state.cells[cellIndex].acsr.injl.temp;
        sourceGasEnthalpy = state.cells[cellIndex].acsr.injl.FluidoPro.EntalpGas(state.cells[cellIndex].pres, sourceTemperature);
        sourceLiquidEnthalpy = state.cells[cellIndex].acsr.injl.FluidoPro.EntalpLiq(state.cells[cellIndex].pres, sourceTemperature);
        hcF = state.cells[cellIndex].acsr.injl.fluidocol.CalorLiq(state.cells[cellIndex].pres, sourceTemperature) * sourceTemperature
            /*entalpia fluido complementar a ser corrigida*/;
    } else if (state.cells[cellIndex].acsr.tipo == 3) {
        sourceTemperature = state.cells[cellIndex].acsr.ipr.Tres;
        sourceGasEnthalpy = state.cells[cellIndex].acsr.ipr.FluidoPro.EntalpGas(state.cells[cellIndex].pres, sourceTemperature);
        sourceLiquidEnthalpy = state.cells[cellIndex].acsr.ipr.FluidoPro.EntalpLiq(state.cells[cellIndex].pres, sourceTemperature);
        hcF = 0;
    } else if (state.cells[cellIndex].acsrL != 0) {
        if ((*state.cells[cellIndex].acsrL).tipo == 5) {
            if ((*state.cells[cellIndex].acsrL).chk.AreaGarg < state.input.master1.razareaativ * state.cells[cellIndex].dutoL.area && (*state.cells[cellIndex].acsrL).chk.AreaGarg > 1e-5 * state.cells[cellIndex].dutoL.area) {
                double chokeUpstreamTemperature = state.cells[cellIndex - 1].temp;
                double chokeUpstreamVoidFraction = state.cells[cellIndex - 1].alf;
                double betE = state.cells[cellIndex - 1].bet;

                double upstreamLiquidDensity = state.cells[cellIndex - 1].flui.MasEspLiq(state.cells[cellIndex - 1].pres, state.cells[cellIndex - 1].temp);
                double rholc = state.cells[cellIndex - 1].fluicol.MasEspFlu(state.cells[cellIndex - 1].pres, state.cells[cellIndex - 1].temp);
                double upstreamMixtureLiquidDensity = (1 - betE) * upstreamLiquidDensity + betE * rholc;

                double chokeDownstreamVoidFraction = state.cells[cellIndex].alf;
                double betJ = state.cells[cellIndex].bet;
                double downstreamLiquidDensity = state.cells[cellIndex].flui.MasEspLiq(state.cells[cellIndex].pres, state.cells[cellIndex].temp);
                double rholcJ = state.cells[cellIndex].fluicol.MasEspFlu(state.cells[cellIndex].pres, state.cells[cellIndex].temp);
                double downstreamMixtureLiquidDensity = (1 - betJ) * downstreamLiquidDensity + betJ * rholcJ;

                double upstreamHydrostaticHead = sin(state.cells[cellIndex - 1].duto.teta) * (0.5 * state.cells[cellIndex - 1].dx) * (upstreamMixtureLiquidDensity * (1 - chokeUpstreamVoidFraction) + chokeUpstreamVoidFraction * state.cells[cellIndex - 1].flui.MasEspGas(state.cells[cellIndex - 1].pres, state.cells[cellIndex - 1].temp)) * 9.82 / 98600.;
                double downstreamHydrostaticHead = sin(state.cells[cellIndex].duto.teta) * (0.5 * state.cells[cellIndex].dx) * (downstreamMixtureLiquidDensity * (1 - chokeDownstreamVoidFraction) + chokeDownstreamVoidFraction * state.cells[cellIndex].flui.MasEspGas(state.cells[cellIndex].pres, state.cells[cellIndex].temp)) * 9.82 / 98600.;

                double quality = chokeUpstreamVoidFraction * state.cells[cellIndex - 1].flui.MasEspGas(state.cells[cellIndex - 1].pres, state.cells[cellIndex - 1].temp) / (state.cells[cellIndex - 1].flui.MasEspGas(state.cells[cellIndex - 1].pres, state.cells[cellIndex - 1].temp) * chokeUpstreamVoidFraction + upstreamMixtureLiquidDensity * (1. - chokeUpstreamVoidFraction));

                double upstreamLiquidJouleThomson = (1. - betE) * state.cells[cellIndex - 1].flui.JTL(state.cells[cellIndex - 1].pres - upstreamHydrostaticHead, state.cells[cellIndex - 1].temp) - betE / rholcJ;
                double upstreamGasJouleThomson = state.cells[cellIndex - 1].flui.JTG(state.cells[cellIndex - 1].pres - upstreamHydrostaticHead, state.cells[cellIndex - 1].temp);
                sourceTemperature = chokeUpstreamTemperature + ((1. - quality) * upstreamLiquidJouleThomson + quality * upstreamGasJouleThomson) * (state.cells[cellIndex].pres + downstreamHydrostaticHead - state.cells[cellIndex - 1].pres - upstreamHydrostaticHead);

                sourceGasEnthalpy = state.cells[cellIndex - 1].flui.EntalpGas(state.cells[cellIndex - 1].pres, sourceTemperature);
                sourceLiquidEnthalpy = state.cells[cellIndex - 1].flui.EntalpLiq(state.cells[cellIndex - 1].pres, sourceTemperature);
                hcF = state.cells[cellIndex - 1].fluicol.CalorLiq(state.cells[cellIndex - 1].pres, sourceTemperature) * sourceTemperature;
                /*entalpia fluido complementar a ser corrigida*/

            } else {
                sourceTemperature = state.cells[cellIndex].temp;
                sourceGasEnthalpy = state.cells[cellIndex - 1].flui.EntalpGas(state.cells[cellIndex - 1].pres, state.cells[cellIndex - 1].temp);
                sourceLiquidEnthalpy = state.cells[cellIndex - 1].flui.EntalpLiq(state.cells[cellIndex - 1].pres, state.cells[cellIndex - 1].temp);
                hcF = state.cells[cellIndex - 1].fluicol.CalorLiq(state.cells[cellIndex - 1].pres, state.cells[cellIndex - 1].temp) * sourceTemperature;
                /*entalpia fluido complementar a ser corrigida*/
            }
        } else if ((*state.cells[cellIndex].acsrL).tipo == 8) {
            double pumpUpstreamVoidFraction = state.cells[cellIndex - 1].alf;
            double betM = state.cells[cellIndex - 1].bet;

            double polytropicExponent = (*state.cells[cellIndex].acsrL).bvol.npoli;
            double pumpPressureRatio = (state.cells[cellIndex].pres) / (state.cells[cellIndex - 1].pres);
            sourceTemperature = state.cells[cellIndex - 1].temp * pow(pumpPressureRatio, (polytropicExponent - 1) / polytropicExponent);

            sourceGasEnthalpy = state.cells[cellIndex - 1].flui.EntalpGas(state.cells[cellIndex - 1].pres, sourceTemperature);
            sourceLiquidEnthalpy = state.cells[cellIndex - 1].flui.EntalpLiq(state.cells[cellIndex - 1].pres, sourceTemperature);
            hcF = state.cells[cellIndex - 1].fluicol.CalorLiq(state.cells[cellIndex - 1].pres, sourceTemperature) * sourceTemperature;
            /*entalpia fluido complementar a ser corrigida*/
        } else {
            sourceGasEnthalpy = 0.;
            sourceLiquidEnthalpy = 0.;
        }
    } else {
        sourceGasEnthalpy = 0.;
        sourceLiquidEnthalpy = 0.;
    }

    liquidMassSourceTerm = 0;
    gasMassSourceTerm = 0;
    if (state.cells[cellIndex].fontemassLR > 0.)
        liquidMassSourceTerm = sourceLiquidEnthalpy * state.cells[cellIndex].fontemassLR / cellLength;
    if (state.cells[cellIndex].fontemassCR > 0.)
        liquidMassSourceTerm += hcF * state.cells[cellIndex].fontemassCR / cellLength;
    if (gasMassSourceTerm > 0.)
        gasMassSourceTerm = sourceGasEnthalpy * state.cells[cellIndex].fontemassGR / cellLength;

    return previousMixtureInternalEnergy - (enthalpyFluxDivergence + hydrostaticTerm - (gasMassSourceTerm + liquidMassSourceTerm) / flowArea) * state.cells[cellIndex].dt;
}

double interpolateMixtureEnergy(const ThermalState &state, int cellIndex, int pressureIndex, int temperatureIndex, double pressureRatio) {
    int nextPressureIndex = pressureIndex + 1;

    double lowerPressure = state.cells[cellIndex].flui.rhogF[pressureIndex][0];
    double upperPressure = state.cells[cellIndex].flui.rhogF[nextPressureIndex][0];
    double temperature = state.cells[cellIndex].flui.rhogF[0][temperatureIndex];
    double meanVoidFraction = state.cells[cellIndex].alf;
    double betmed = state.cells[cellIndex].bet;

    double lowerGasDensity = state.cells[cellIndex].flui.rhogF[pressureIndex][temperatureIndex];
    double upperGasDensity = state.cells[cellIndex].flui.rhogF[nextPressureIndex][temperatureIndex];
    double lowerGasEnthalpy = state.cells[cellIndex].flui.HgF[pressureIndex][temperatureIndex];
    double upperGasEnthalpy = state.cells[cellIndex].flui.HgF[nextPressureIndex][temperatureIndex];

    double lowerLiquidDensity = state.cells[cellIndex].flui.rholF[pressureIndex][temperatureIndex];
    double upperLiquidDensity = state.cells[cellIndex].flui.rholF[nextPressureIndex][temperatureIndex];
    double lowerLiquidEnthalpy = state.cells[cellIndex].flui.HlF[pressureIndex][temperatureIndex];
    double upperLiquidEnthalpy = state.cells[cellIndex].flui.HlF[nextPressureIndex][temperatureIndex];

    double rhocp0 = state.cells[cellIndex].fluicol.MasEspFlu(lowerPressure, temperature);
    double rhocp1 = state.cells[cellIndex].fluicol.MasEspFlu(upperPressure, temperature);
    double hlc0 = state.cells[cellIndex].fluicol.CalorLiq(lowerPressure, temperature) * temperature; // corrigir entalpia
    double hlc1 = state.cells[cellIndex].fluicol.CalorLiq(upperPressure, temperature) * temperature; // corrigir en;talpia

    double lowerMixtureEnergy = meanVoidFraction * lowerGasDensity * (lowerGasEnthalpy - lowerPressure * 98066.5 / lowerGasDensity) +
                    (1 - meanVoidFraction) * (1 - betmed) * lowerLiquidDensity * (lowerLiquidEnthalpy - lowerPressure * 98066.5 / lowerLiquidDensity) +
                    (1 - meanVoidFraction) * (betmed)*rhocp0 * (hlc0 - lowerPressure * 98066.5 / rhocp0);

    double upperMixtureEnergy = meanVoidFraction * upperGasDensity * (upperGasEnthalpy - upperPressure * 98066.5 / lowerGasDensity) +
                    (1 - meanVoidFraction) * (1 - betmed) * upperLiquidDensity * (upperLiquidEnthalpy - upperPressure * 98066.5 / upperLiquidDensity) +
                    (1 - meanVoidFraction) * (betmed)*rhocp1 * (hlc1 - upperPressure * 98066.5 / rhocp1);

    return pressureRatio * lowerMixtureEnergy + (1. - pressureRatio) * upperMixtureEnergy;
}

void updateTemperatureFromEnthalpy(const ThermalState &state, int cellIndex) {
    double pressure = state.cells[cellIndex].pres;
    double **propertyTable = state.cells[cellIndex].flui.rholF;

    int pressureIndex = 0;
    int pressureSearchMarker;
    int divisionCount = state.cells[cellIndex].flui.npontos - 1;
    if (pressure < propertyTable[1][0] || pressure >= propertyTable[divisionCount + 1][0]) {
        cout << "pressÃƒÂ£o fora dos limites de tabela";
        getchar();
    }

    int searchLow, searchMiddle, searchHigh;
    searchLow = 1;
    searchHigh = divisionCount + 1;
    while (searchLow <= searchHigh) {
        searchMiddle = (searchLow + searchHigh) / 2;
        pressureSearchMarker = searchMiddle;
        if (searchMiddle == 1) {
            pressureIndex = searchMiddle;
            break;
        } else if (searchMiddle == divisionCount + 1 && propertyTable[searchMiddle][0] == pressure) {
            pressureIndex = searchMiddle - 1;
            break;
        }
        if (propertyTable[searchMiddle][0] > pressure && propertyTable[searchMiddle - 1][0] <= pressure) {
            pressureIndex = searchMiddle - 1;
            break;
        }
        if (propertyTable[searchMiddle][0] < pressure)
            searchLow = searchMiddle + 1;
        else
            searchHigh = searchMiddle - 1;
    }

    double pressureRatio = 1. - (propertyTable[pressureIndex][0] - pressure) / (propertyTable[pressureIndex][0] - propertyTable[pressureIndex + 1][0]);

    double mixtureInternalEnergy = computeMixtureEnthalpy(state, cellIndex);
    int temperatureIndex = 1;
    double lowerEnergy = -1;
    double upperEnergy = 0;
    while (temperatureIndex < divisionCount + 1 || (mixtureInternalEnergy >= lowerEnergy && mixtureInternalEnergy <= upperEnergy) || (mixtureInternalEnergy <= lowerEnergy && mixtureInternalEnergy >= upperEnergy)) {
        lowerEnergy = interpolateMixtureEnergy(state, cellIndex, pressureIndex, temperatureIndex, pressureRatio);
        upperEnergy = interpolateMixtureEnergy(state, cellIndex, pressureIndex, temperatureIndex + 1, pressureRatio);
        temperatureIndex++;
    }
    double temperatureRatio = 1. - (lowerEnergy - mixtureInternalEnergy) / (lowerEnergy - upperEnergy);
    state.cells[cellIndex].temp = propertyTable[0][temperatureIndex] * temperatureRatio + (1. - temperatureRatio) * propertyTable[0][temperatureIndex + 1];

    if (state.cells[cellIndex].temp < -50.)
        state.cells[cellIndex].temp = -50.;
    if (state.cells[cellIndex].temp > 200.)
        state.cells[cellIndex].temp = 200.;
}

namespace {

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

TemperatureBalance prepareTemperatureBalance(const ThermalState &state,
                                             int cellIndex,
                                             int steadyStateMode) {
    double cellLength = state.cells[cellIndex].dx;
    double meanCellLength;
    if (cellIndex > 0)
        meanCellLength = 0.5 * (state.cells[cellIndex].dx + state.cells[cellIndex - 1].dx);
    else
        meanCellLength = 0.5 * state.cells[cellIndex].dx;
    double diameter = state.cells[cellIndex].duto.a;
    double flowArea = 0.25 * M_PI * diameter * diameter;
    double meanVoidFraction = state.cells[cellIndex].alf;
    double betmed = state.cells[cellIndex].bet;
    double meanSuperficialGasVelocity;
    double meanSuperficialLiquidVelocity;
    if (cellIndex > 0 && (state.cells[cellIndex - 1].acsr.tipo != 5 ||
                          state.cells[cellIndex - 1].acsr.chk.AreaGarg > (1e-3 + state.input.master1.razareaativ) * state.cells[cellIndex - 1].duto.area)) {
        if (state.cells[cellIndex].alf > (*state.globals).localtiny)
            meanSuperficialGasVelocity = state.cells[cellIndex].QG / flowArea;
        else {
            meanSuperficialGasVelocity = 0.;
        }
        if (state.cells[cellIndex].alf < 1. - (*state.globals).localtiny)
            meanSuperficialLiquidVelocity = state.cells[cellIndex].QL / flowArea;
        else {
            meanSuperficialLiquidVelocity = 0.;
        }
    } else {
        if (state.cells[cellIndex].alf > (*state.globals).localtiny)
            meanSuperficialGasVelocity = state.cells[cellIndex + 1].QG / flowArea;
        else {
            meanSuperficialGasVelocity = 0.;
        }
        if (state.cells[cellIndex].alf < 1. - (*state.globals).localtiny)
            meanSuperficialLiquidVelocity = state.cells[cellIndex + 1].QL / flowArea;
        else {
            meanSuperficialLiquidVelocity = 0.;
        }
    }
    double referenceMixtureVelocity = fabs(meanSuperficialGasVelocity) + fabs(meanSuperficialLiquidVelocity);
    double rp = state.cells[cellIndex].rpC;
    double rc = state.cells[cellIndex].rcC;
    double liquidDensity = (1. - betmed) * rp + betmed * rc;
    double gasDensity = state.cells[cellIndex].rgC;
    double liquidSpecificHeat = (1. - betmed) * state.cells[cellIndex].flui.CalorLiq(state.cells[cellIndex].presini, state.cells[cellIndex].temp) + betmed * state.cells[cellIndex].fluicol.CalorLiq(state.cells[cellIndex].presini, state.cells[cellIndex].temp);
    double liquidSpecificHeatConstantVolume = liquidSpecificHeat;
    double gasSpecificHeat = state.cells[cellIndex].flui.CalorGas(state.cells[cellIndex].presini, state.cells[cellIndex].temp);
    double gasSpecificHeatConstantVolume = state.cells[cellIndex].flui.CalorGasVolMod(state.cells[cellIndex].presini, state.cells[cellIndex].temp, state.cells[cellIndex].rgC);
    double liquidJouleThomson = (1. - betmed) * state.cells[cellIndex].flui.JTL(state.cells[cellIndex].presini, state.cells[cellIndex].temp) - betmed / rc;
    double gasJouleThomson = state.cells[cellIndex].flui.JTG(state.cells[cellIndex].presini, state.cells[cellIndex].temp, state.cells[cellIndex].rgC);
    double hydrostaticPower = (liquidDensity * meanSuperficialLiquidVelocity + gasDensity * meanSuperficialGasVelocity) * flowArea * 9.82 * sin(state.cells[cellIndex].duto.teta);

    double heatFlux = 0.;
    state.cells[cellIndex].calor.Tint = state.cells[cellIndex].temp;
    if (cellIndex > 0)
        state.cells[cellIndex].calor.dtL = state.cells[cellIndex].temp - state.cells[cellIndex - 1].tempini;
    else
        state.cells[cellIndex].calor.dtL = 0;
    state.cells[cellIndex].calor.Vint = meanSuperficialGasVelocity + meanSuperficialLiquidVelocity;
    state.cells[cellIndex].calor.dt = state.cells[cellIndex].dt;
    double liquidConductivity = (1. - betmed) * state.cells[cellIndex].flui.CondLiq(state.cells[cellIndex].presini, state.cells[cellIndex].temp) + betmed * state.cells[cellIndex].fluicol.CondLiq(state.cells[cellIndex].presini, state.cells[cellIndex].temp);
    state.cells[cellIndex].calor.kint = liquidConductivity * (1 - meanVoidFraction) + state.cells[cellIndex].flui.CondGas(state.cells[cellIndex].presini, state.cells[cellIndex].temp) * meanVoidFraction;
    state.cells[cellIndex].calor.cpint = liquidSpecificHeat * (1 - meanVoidFraction) + gasSpecificHeat * meanVoidFraction;
    state.cells[cellIndex].calor.rhoint = liquidDensity * (1 - meanVoidFraction) + gasDensity * meanVoidFraction;
    double liquidViscosity = (1. - betmed) * state.cells[cellIndex].mipC + betmed * state.cells[cellIndex].micC;
    state.cells[cellIndex].calor.viscint = liquidViscosity * (1 - meanVoidFraction) * 1.e-3 + state.cells[cellIndex].migC * meanVoidFraction * 1.e-3;
    double temperaturePerturbation = state.cells[cellIndex].temp * 0.01;
    if (fabs(state.cells[cellIndex].temp) < 1e-15)
        temperaturePerturbation = 0.1;
    double liquidDensityChange = (1. - betmed) * state.cells[cellIndex].flui.MasEspLiq(state.cells[cellIndex].presini, state.cells[cellIndex].temp + temperaturePerturbation) +
                    betmed * state.cells[cellIndex].fluicol.MasEspFlu(state.cells[cellIndex].presini, state.cells[cellIndex].temp + temperaturePerturbation) - liquidDensity;
    double gasDensityChange = state.cells[cellIndex].flui.MasEspGas(state.cells[cellIndex].presini, state.cells[cellIndex].temp + temperaturePerturbation) - gasDensity;
    state.cells[cellIndex].calor.betint = -(1 / state.cells[cellIndex].calor.rhoint) * (liquidDensityChange * (1 - meanVoidFraction) + gasDensityChange * meanVoidFraction) / (temperaturePerturbation);
    if (state.input.modoDifus3D == 0 || steadyStateMode != 0) {
        if (steadyStateMode == 0)
            heatFlux = state.cells[cellIndex].calor.transtrans();
        else
            heatFlux = state.cells[cellIndex].calor.transperm();
        if (state.productionNetworkCoupled == 1 && cellIndex >= state.primarySectionStart && cellIndex <= state.primarySectionEnd) {
            heatFlux -= state.cells[cellIndex].fluxcalAcopRedeP;
        }
    } else if (state.input.modoDifus3D == 1) {
        int coupled = -1;
        int icelAcop;
        for (int iacop = 0; iacop < state.input.nacop; iacop++) {
            icelAcop = state.input.celAcop[iacop].indCel;
            if (cellIndex == icelAcop) {
                coupled = iacop;
                break;
            }
        }
        if (coupled == -1)
            heatFlux = state.cells[cellIndex].calor.transtrans();
        else {
            int iacop1 = state.coupledCellIndices[coupled];
            heatFlux = -state.input.celAcop[coupled].FE * state.poissonSolver.dados.qTotal[iacop1] / state.cells[cellIndex].dx;
        }
    }

    state.cells[cellIndex].fluxcalmed = heatFlux;

    double timeCoefficient = (liquidDensity * (1 - meanVoidFraction) * liquidSpecificHeatConstantVolume + gasDensity * meanVoidFraction * gasSpecificHeatConstantVolume) * flowArea;
    double pressureTimeCoefficient;
    pressureTimeCoefficient = -state.cells[cellIndex].flui.CalorGasPresMod(state.cells[cellIndex].presini, state.cells[cellIndex].temp) * (gasDensity * meanVoidFraction * flowArea); //-state.cells[cellIndex].flui.CalorGasPresMod(state.cells[cellIndex].pres, state.cells[cellIndex].temp) * (gasDensity * meanVoidFraction * flowArea);
    double temperatureSpatialCoefficient = (liquidDensity * meanSuperficialLiquidVelocity * liquidSpecificHeat + gasDensity * meanSuperficialGasVelocity * gasSpecificHeat) * flowArea;
    double pressureSpatialCoefficient = (liquidDensity * meanSuperficialLiquidVelocity * liquidJouleThomson + gasDensity * meanSuperficialGasVelocity * gasJouleThomson) * flowArea;

    return TemperatureBalance{
        .cellLength = cellLength,
        .meanCellLength = meanCellLength,
        .flowArea = flowArea,
        .voidFraction = meanVoidFraction,
        .bet = betmed,
        .gasSuperficialVelocity = meanSuperficialGasVelocity,
        .liquidSuperficialVelocity = meanSuperficialLiquidVelocity,
        .referenceMixtureVelocity = referenceMixtureVelocity,
        .rp = rp,
        .rc = rc,
        .liquidDensity = liquidDensity,
        .gasDensity = gasDensity,
        .liquidHeatCapacity = liquidSpecificHeat,
        .liquidIsochoricHeatCapacity = liquidSpecificHeatConstantVolume,
        .gasHeatCapacity = gasSpecificHeat,
        .gasIsochoricHeatCapacity = gasSpecificHeatConstantVolume,
        .liquidJouleThomson = liquidJouleThomson,
        .gasJouleThomson = gasJouleThomson,
        .hydrostaticPower = hydrostaticPower,
        .heatFlux = heatFlux,
        .timeCoefficient = timeCoefficient,
        .pressureTimeCoefficient = pressureTimeCoefficient,
        .temperatureSpatialCoefficient = temperatureSpatialCoefficient,
        .pressureSpatialCoefficient = pressureSpatialCoefficient,
    };
}

double computeKineticTemperatureTerm(const ThermalState &state, int cellIndex,
                                     TemperatureBalance &balance) {
    double kineticTerm = 0;
    double upstreamMeanGasVelocity = 0;
    double upstreamMeanLiquidVelocity = 0;
    double meanGasVelocity = 0;
    double meanLiquidVelocity = 0;
    double initialMeanGasVelocity = 0;
    double initialMeanLiquidVelocity = 0;
    if (cellIndex <= state.lastCell - 1 && state.cells[cellIndex].acsr.tipo == 0 && state.cells[cellIndex + 1].acsr.tipo == 0) {

        double kineticCellLength = state.cells[cellIndex].dx;
        double upstreamDiameter = state.cells[cellIndex].duto.a;
        double upstreamFlowArea = 0.25 * M_PI * upstreamDiameter * upstreamDiameter;

        double rpcin = state.cells[cellIndex + 1].rpCi;
        double rccin = state.cells[cellIndex + 1].rcCi;
        double kineticLiquidDensity = (1. - balance.bet) * rpcin + balance.bet * rccin;
        double kineticGasDensity = state.cells[cellIndex + 1].rgCi;

        double faceVoidFraction;
        double initialVoidFraction;

        balance.gasSuperficialVelocity = state.cells[cellIndex + 1].QG / balance.flowArea;
        balance.liquidSuperficialVelocity = state.cells[cellIndex + 1].QL / balance.flowArea;

        if (state.cells[cellIndex + 1].QG > 0)
            faceVoidFraction = state.cells[cellIndex].alf;
        else
            faceVoidFraction = state.cells[cellIndex + 1].alf;

        if (state.cells[cellIndex + 1].QGini > 0)
            initialVoidFraction = state.cells[cellIndex].alfini;
        else
            initialVoidFraction = state.cells[cellIndex + 1].alfini;

        double upstreamFaceVoidFraction;
        if (state.cells[cellIndex].QG > 0)
            upstreamFaceVoidFraction = state.cells[cellIndex - 1].alf;
        else
            upstreamFaceVoidFraction = state.cells[cellIndex].alf;

        if (upstreamFaceVoidFraction > 1e-3) {
            upstreamMeanGasVelocity = state.cells[cellIndex].QG / (upstreamFlowArea);
            upstreamMeanGasVelocity /= upstreamFaceVoidFraction;
        }

        if (upstreamFaceVoidFraction < 1. - 1e-3) {
            upstreamMeanLiquidVelocity = state.cells[cellIndex].QL / (upstreamFlowArea);
            upstreamMeanLiquidVelocity /= (1. - upstreamFaceVoidFraction);
        }

        if (faceVoidFraction > 1e-3) {
            meanGasVelocity = balance.gasSuperficialVelocity;
            meanGasVelocity /= faceVoidFraction;
        }

        if (faceVoidFraction < 1. - 1e-3) {
            meanLiquidVelocity = balance.liquidSuperficialVelocity;
            meanLiquidVelocity /= (1. - faceVoidFraction);
        }

        if (initialVoidFraction > 1e-3) {

            initialMeanGasVelocity = state.cells[cellIndex + 1].QGini / (upstreamFlowArea);
            initialMeanGasVelocity /= initialVoidFraction;
        }

        if (initialVoidFraction < 1. - 1e-3) {

            initialMeanLiquidVelocity = state.cells[cellIndex + 1].QLini / (upstreamFlowArea);
            initialMeanLiquidVelocity /= (1. - initialVoidFraction);
        }

        kineticTerm = kineticLiquidDensity * (1 - faceVoidFraction) * upstreamFlowArea * (0.5 * (meanLiquidVelocity * meanLiquidVelocity - initialMeanLiquidVelocity * initialMeanLiquidVelocity)) / state.cells[cellIndex].dt +
                   kineticGasDensity * faceVoidFraction * upstreamFlowArea * (0.5 * (meanGasVelocity * meanGasVelocity - initialMeanGasVelocity * initialMeanGasVelocity)) / state.cells[cellIndex].dt +
                   ((state.cells[cellIndex + 1].MC - state.cells[cellIndex + 1].Mliqini) * meanGasVelocity * (meanGasVelocity - upstreamMeanGasVelocity) / kineticCellLength +
                    state.cells[cellIndex + 1].Mliqini * meanLiquidVelocity * (meanLiquidVelocity - upstreamMeanLiquidVelocity) / kineticCellLength);
    }
    return kineticTerm;
}

struct TemperatureSourceTerms {
    double gas;
    double liquid;
};

TemperatureSourceTerms computeTemperatureSourceTerms(const ThermalState &state,
                                                      int cellIndex,
                                                      double cellLength) {
    double gasMassSourceTerm = 0.;
    double liquidMassSourceTerm = 0.;
    double fontemassC = 0.;
    double sourceTemperature = state.cells[cellIndex].temp;
    double sourceGasSpecificHeat;
    double sourceSpecificHeatRatio = 0.;
    double sourceLiquidSpecificHeat;
    if (state.cells[cellIndex].acsr.tipo == 1) {
        sourceTemperature = state.cells[cellIndex].acsr.injg.temp;
        sourceGasSpecificHeat = state.cells[cellIndex].acsr.injg.FluidoPro.CalorGas(state.cells[cellIndex].presini, state.cells[cellIndex].temp);
        sourceSpecificHeatRatio = state.cells[cellIndex].acsr.injg.FluidoPro.ConstAdG(state.cells[cellIndex].presini, state.cells[cellIndex].temp);
        sourceLiquidSpecificHeat = 0.;
    } else if (state.cells[cellIndex].acsr.tipo == 2) {
        sourceTemperature = state.cells[cellIndex].acsr.injl.temp;
        sourceGasSpecificHeat = state.cells[cellIndex].acsr.injl.FluidoPro.CalorGas(state.cells[cellIndex].presini, state.cells[cellIndex].temp);                                                                                                                                        // state.cells[cellIndex].acsr.injl.FluidoPro.CalorGas(state.cells[cellIndex].pres, state.cells[cellIndex].temp);
        sourceSpecificHeatRatio = state.cells[cellIndex].acsr.injl.FluidoPro.CalorGas(state.cells[cellIndex].presini, state.cells[cellIndex].temp);                                                                                                                                      // state.cells[cellIndex].acsr.injl.FluidoPro.ConstAdG(state.cells[cellIndex].pres, state.cells[cellIndex].temp);
        sourceLiquidSpecificHeat = (1. - state.cells[cellIndex].acsr.injl.bet) * state.cells[cellIndex].acsr.injl.FluidoPro.CalorLiq(state.cells[cellIndex].presini, state.cells[cellIndex].temp) + state.cells[cellIndex].acsr.injl.bet * state.cells[cellIndex].acsr.injl.fluidocol.CalorLiq(state.cells[cellIndex].presini, state.cells[cellIndex].temp); //(1. - state.cells[cellIndex].acsr.injl.bet) * state.cells[cellIndex].acsr.injl.FluidoPro.CalorLiq(state.cells[cellIndex].pres, state.cells[cellIndex].temp)
    } else if (state.cells[cellIndex].acsr.tipo == 10) {
        sourceTemperature = state.cells[cellIndex].acsr.injm.temp;
        sourceGasSpecificHeat = state.cells[cellIndex].acsr.injm.FluidoPro.CalorGas(state.cells[cellIndex].presini, state.cells[cellIndex].temp);
        sourceSpecificHeatRatio = state.cells[cellIndex].acsr.injm.FluidoPro.ConstAdG(state.cells[cellIndex].presini, state.cells[cellIndex].temp);
        if ((state.cells[cellIndex].fontemassLR + state.cells[cellIndex].fontemassCR) > 0.) {
            double titbet = state.cells[cellIndex].fontemassCR / (state.cells[cellIndex].fontemassLR + state.cells[cellIndex].fontemassCR);
            sourceLiquidSpecificHeat = (1. - titbet) * state.cells[cellIndex].acsr.injl.FluidoPro.CalorLiq(state.cells[cellIndex].presini, state.cells[cellIndex].temp) +
                   titbet * state.cells[cellIndex].acsr.injl.fluidocol.CalorLiq(state.cells[cellIndex].presini, state.cells[cellIndex].temp);
        } else
            sourceLiquidSpecificHeat = 0.;
    } else if (state.cells[cellIndex].acsr.tipo == 3) {
        sourceTemperature = state.cells[cellIndex].acsr.ipr.Tres;
        sourceGasSpecificHeat = state.cells[cellIndex].acsr.ipr.FluidoPro.CalorGas(state.cells[cellIndex].presini, state.cells[cellIndex].temp);
        sourceSpecificHeatRatio = state.cells[cellIndex].acsr.ipr.FluidoPro.ConstAdG(state.cells[cellIndex].presini, state.cells[cellIndex].temp);
        sourceLiquidSpecificHeat = state.cells[cellIndex].acsr.ipr.FluidoPro.CalorLiq(state.cells[cellIndex].presini, state.cells[cellIndex].temp);
    } else if (state.cells[cellIndex].acsr.tipo == 15) {
        sourceTemperature = state.cells[cellIndex].acsr.radialPoro.tRes;
        sourceGasSpecificHeat = state.cells[cellIndex].acsr.radialPoro.flup.CalorGas(state.cells[cellIndex].presini, state.cells[cellIndex].temp);
        sourceSpecificHeatRatio = state.cells[cellIndex].acsr.radialPoro.flup.ConstAdG(state.cells[cellIndex].presini, state.cells[cellIndex].temp);
        sourceLiquidSpecificHeat = state.cells[cellIndex].acsr.radialPoro.flup.CalorLiq(state.cells[cellIndex].presini, state.cells[cellIndex].temp);
    } else if (state.cells[cellIndex].acsr.tipo == 9) {
        double betM = state.cells[cellIndex].acsr.fontechk.betISamb;
        double ambientPressure = state.cells[cellIndex].acsr.fontechk.pamb;
        double ambientTemperature = state.cells[cellIndex].acsr.fontechk.tamb;
        sourceGasSpecificHeat = state.cells[cellIndex].acsr.fontechk.fluidoPamb.CalorGas(ambientPressure, ambientTemperature);
        sourceLiquidSpecificHeat = (1. - betM) * state.cells[cellIndex].acsr.fontechk.fluidoPamb.CalorLiq(ambientPressure, ambientTemperature) + betM * state.cells[cellIndex].acsr.fontechk.fluidocol.CalorLiq(ambientPressure, ambientTemperature);
        sourceSpecificHeatRatio = state.cells[cellIndex].acsr.fontechk.fluidoPamb.ConstAdG(ambientPressure, ambientTemperature);
        sourceTemperature = state.cells[cellIndex].acsr.fontechk.tamb;
    } else if (state.cells[cellIndex].acsrL != 0 && cellIndex < state.lastCell) {
        if ((*state.cells[cellIndex].acsrL).tipo == 5 && (state.cells[cellIndex].fontemassLR + state.cells[cellIndex].fontemassCR + state.cells[cellIndex].fontemassGR) > 0.) {
            if ((*state.cells[cellIndex].acsrL).chk.AreaGarg < state.input.master1.razareaativ * state.cells[cellIndex].dutoL.area && (*state.cells[cellIndex].acsrL).chk.AreaGarg > 1e-5 * state.cells[cellIndex].dutoL.area) {
                double upstreamTemperature = state.cells[cellIndex - 1].tempini;
                double chokeUpstreamVoidFraction = state.cells[cellIndex - 1].alf;
                double betE = state.cells[cellIndex - 1].bet;

                double upstreamLiquidDensity = state.cells[cellIndex - 1].rpC;
                double rholc = state.cells[cellIndex - 1].rcC;
                double upstreamMixtureLiquidDensity = (1 - betE) * upstreamLiquidDensity + betE * rholc;

                double chokeDownstreamVoidFraction = state.cells[cellIndex].alf;
                double betJ = state.cells[cellIndex].bet;
                double downstreamLiquidDensity = state.cells[cellIndex].rpC;
                double rholcJ = state.cells[cellIndex].rcC;
                double downstreamMixtureLiquidDensity = (1 - betJ) * downstreamLiquidDensity + betJ * rholcJ;

                double upstreamHydrostaticHead = sin(state.cells[cellIndex - 1].duto.teta) * (0.5 * state.cells[cellIndex - 1].dx) * (upstreamMixtureLiquidDensity * (1 - chokeUpstreamVoidFraction) + chokeUpstreamVoidFraction * state.cells[cellIndex - 1].rgC) * 9.82 / 98600.;
                double downstreamHydrostaticHead = sin(state.cells[cellIndex].duto.teta) * (0.5 * state.cells[cellIndex].dx) * (downstreamMixtureLiquidDensity * (1 - chokeDownstreamVoidFraction) + chokeDownstreamVoidFraction * state.cells[cellIndex].rgC) * 9.82 / 98600.;
                double quality = chokeUpstreamVoidFraction * state.cells[cellIndex - 1].rgC / (state.cells[cellIndex - 1].rgC * chokeUpstreamVoidFraction + upstreamMixtureLiquidDensity * (1. - chokeUpstreamVoidFraction));
                sourceGasSpecificHeat = state.cells[cellIndex - 1].flui.CalorGas(state.cells[cellIndex - 1].presini, state.cells[cellIndex - 1].tempini);                                                                                                                              // state.cells[cellIndex - 1].flui.CalorGas(state.cells[cellIndex - 1].pres, state.cells[cellIndex - 1].temp);
                sourceLiquidSpecificHeat = (1 - state.cells[cellIndex - 1].bet) * state.cells[cellIndex - 1].flui.CalorLiq(state.cells[cellIndex - 1].presini, state.cells[cellIndex - 1].tempini) + state.cells[cellIndex - 1].bet * state.cells[cellIndex - 1].fluicol.CalorLiq(state.cells[cellIndex - 1].presini, state.cells[cellIndex - 1].tempini); //(1 - state.cells[cellIndex - 1].bet) * state.cells[cellIndex - 1].flui.CalorLiq(state.cells[cellIndex - 1].pres, state.cells[cellIndex - 1].temp)
                double upstreamLiquidJouleThomson = (1. - betE) * state.cells[cellIndex - 1].flui.JTL(state.cells[cellIndex - 1].presini - upstreamHydrostaticHead, state.cells[cellIndex - 1].tempini) - betE / rholcJ;                                                                                     //(1. - betE) * state.cells[cellIndex - 1].flui.JTL(state.cells[cellIndex - 1].pres - upstreamHydrostaticHead, state.cells[cellIndex - 1].temp)
                double upstreamGasJouleThomson = state.cells[cellIndex - 1].flui.JTG(state.cells[cellIndex - 1].presini - upstreamHydrostaticHead, state.cells[cellIndex - 1].tempini);                                                                                                                   // state.cells[cellIndex - 1].flui.JTG(state.cells[cellIndex - 1].pres - upstreamHydrostaticHead, state.cells[cellIndex - 1].temp);
                sourceTemperature = upstreamTemperature + ((1. - quality) * upstreamLiquidJouleThomson / sourceLiquidSpecificHeat + quality * upstreamGasJouleThomson / sourceGasSpecificHeat) *
                                  (state.cells[cellIndex].pres + downstreamHydrostaticHead - state.cells[cellIndex - 1].pres - upstreamHydrostaticHead) * 98066.52;

                sourceSpecificHeatRatio = state.cells[cellIndex - 1].flui.ConstAdG(state.cells[cellIndex - 1].presini, state.cells[cellIndex - 1].tempini);

            } else {
                sourceGasSpecificHeat = state.cells[cellIndex - 1].flui.CalorGas(state.cells[cellIndex - 1].presini, state.cells[cellIndex - 1].tempini);
                sourceLiquidSpecificHeat = (1 - state.cells[cellIndex - 1].bet) * state.cells[cellIndex - 1].flui.CalorLiq(state.cells[cellIndex - 1].presini, state.cells[cellIndex - 1].tempini) + state.cells[cellIndex - 1].bet * state.cells[cellIndex - 1].fluicol.CalorLiq(state.cells[cellIndex - 1].presini, state.cells[cellIndex - 1].tempini);
                sourceTemperature = state.cells[cellIndex].temp;
                sourceSpecificHeatRatio = state.cells[cellIndex - 1].flui.ConstAdG(state.cells[cellIndex - 1].presini, state.cells[cellIndex - 1].tempini);
            }
        } else if ((*state.cells[cellIndex].acsrL).tipo == 8) {
            double betM = state.cells[cellIndex - 1].bet;
            sourceGasSpecificHeat = state.cells[cellIndex - 1].flui.CalorGas(state.cells[cellIndex - 1].presini, state.cells[cellIndex - 1].tempini);
            sourceLiquidSpecificHeat = (1. - betM) * state.cells[cellIndex - 1].flui.CalorLiq(state.cells[cellIndex - 1].presini, state.cells[cellIndex - 1].tempini) + betM * state.cells[cellIndex - 1].fluicol.CalorLiq(state.cells[cellIndex - 1].presini, state.cells[cellIndex - 1].tempini);
            sourceSpecificHeatRatio = state.cells[cellIndex - 1].flui.ConstAdG(state.cells[cellIndex - 1].presini, state.cells[cellIndex - 1].tempini);
            double n = (*state.cells[cellIndex].acsrL).bvol.npoli;
            double pumpPressureRatio = (state.cells[cellIndex].pres) / (state.cells[cellIndex - 1].pres);
            sourceTemperature = state.cells[cellIndex - 1].tempini * pow(pumpPressureRatio, (n - 1) / n);
        } else {
            sourceGasSpecificHeat = 0.;
            sourceSpecificHeatRatio = 1.;
            sourceLiquidSpecificHeat = 0.;
        }
    } else if (cellIndex == state.lastCell && (state.cells[cellIndex].fontemassLR + state.cells[cellIndex].fontemassCR + state.cells[cellIndex].fontemassGR) > 0.) {
        if ((*state.globals).chaverede == 0 || state.networkEndpoint == 1 || (*state.globals).chaveRedeParalela == 1)
            sourceTemperature = state.cells[cellIndex].calor.Textern1;
        else
            sourceTemperature = state.gasSurfaceTemperature;
        sourceGasSpecificHeat = state.cells[cellIndex].flui.CalorGas(state.cells[cellIndex].presini, sourceTemperature);
        double betM = 0.;
        if (fabs(state.cells[cellIndex].fontemassCR) > 1e-15)
            betM = state.cells[cellIndex].fontemassCR /
                   (state.cells[cellIndex].fontemassCR + state.cells[cellIndex].fontemassLR);
        sourceLiquidSpecificHeat = (1. - betM) * state.cells[cellIndex].flui.CalorLiq(state.cells[cellIndex].presini, sourceTemperature) + betM * state.cells[cellIndex].fluicol.CalorLiq(state.cells[cellIndex].presini, sourceTemperature);
        sourceSpecificHeatRatio = state.cells[cellIndex].flui.ConstAdG(state.cells[cellIndex].pres, sourceTemperature);
    } else {
        sourceGasSpecificHeat = 0.;
        sourceSpecificHeatRatio = 1.;
        sourceLiquidSpecificHeat = 0.;
    }

    liquidMassSourceTerm = 0;
    if (state.cells[cellIndex].fontemassLR > 0.)
        liquidMassSourceTerm = state.cells[cellIndex].fontemassLR / cellLength;
    if (state.cells[cellIndex].fontemassCR > 0.)
        liquidMassSourceTerm += state.cells[cellIndex].fontemassCR / cellLength;
    liquidMassSourceTerm *= (sourceLiquidSpecificHeat) * (sourceTemperature - state.cells[cellIndex].temp);

    gasMassSourceTerm = state.cells[cellIndex].fontemassGR / cellLength;
    if (gasMassSourceTerm > 0.)
        gasMassSourceTerm *= (sourceGasSpecificHeat / sourceSpecificHeatRatio) * (sourceTemperature - state.cells[cellIndex].temp);
    else
        gasMassSourceTerm = 0;

    return TemperatureSourceTerms{
        .gas = gasMassSourceTerm,
        .liquid = liquidMassSourceTerm,
    };
}

TemperatureSourceTerms computeThermalMassTransferSourceTerms(
    const ThermalState &state, int cellIndex, double cellLength) {
    double gasMassSourceTerm = 0.;
    double liquidMassSourceTerm = 0.;
    double sourceTemperature = state.cells[cellIndex].temp;
    double sourceGasSpecificHeat;
    double sourceSpecificHeatRatio = 0.;
    double sourceLiquidSpecificHeat;
    if (state.cells[cellIndex].acsr.tipo == 1) {
        sourceTemperature = state.cells[cellIndex].acsr.injg.temp;
        sourceGasSpecificHeat = state.cells[cellIndex].acsr.injg.FluidoPro.CalorGas(state.cells[cellIndex].pres, state.cells[cellIndex].temp);
        sourceSpecificHeatRatio = state.cells[cellIndex].acsr.injg.FluidoPro.ConstAdG(state.cells[cellIndex].pres, state.cells[cellIndex].temp);
        sourceLiquidSpecificHeat = 0.;
    } else if (state.cells[cellIndex].acsr.tipo == 2) {
        sourceTemperature = state.cells[cellIndex].acsr.injl.temp;
        ;
        sourceGasSpecificHeat = state.cells[cellIndex].acsr.injl.FluidoPro.CalorGas(state.cells[cellIndex].pres, state.cells[cellIndex].temp);
        sourceSpecificHeatRatio = state.cells[cellIndex].acsr.injl.FluidoPro.ConstAdG(state.cells[cellIndex].pres, state.cells[cellIndex].temp);
        sourceLiquidSpecificHeat = (1. - state.cells[cellIndex].acsr.injl.bet) * state.cells[cellIndex].acsr.injl.FluidoPro.CalorLiq(state.cells[cellIndex].pres, state.cells[cellIndex].temp) + state.cells[cellIndex].acsr.injl.bet * state.cells[cellIndex].acsr.injl.fluidocol.CalorLiq(state.cells[cellIndex].pres, state.cells[cellIndex].temp);
    } else if (state.cells[cellIndex].acsr.tipo == 3) {
        sourceTemperature = state.cells[cellIndex].acsr.ipr.Tres;
        sourceGasSpecificHeat = state.cells[cellIndex].acsr.ipr.FluidoPro.CalorGas(state.cells[cellIndex].pres, state.cells[cellIndex].temp);
        sourceSpecificHeatRatio = state.cells[cellIndex].acsr.ipr.FluidoPro.ConstAdG(state.cells[cellIndex].pres, state.cells[cellIndex].temp);
        sourceLiquidSpecificHeat = state.cells[cellIndex].acsr.ipr.FluidoPro.CalorLiq(state.cells[cellIndex].pres, state.cells[cellIndex].temp);
    } else if (state.cells[cellIndex].acsrL != 0) {
        if ((*state.cells[cellIndex].acsrL).tipo == 5) {
            if ((*state.cells[cellIndex].acsrL).chk.AreaGarg < state.input.master1.razareaativ * state.cells[cellIndex].dutoL.area && (*state.cells[cellIndex].acsrL).chk.AreaGarg > 1e-5 * state.cells[cellIndex].dutoL.area) {
                double upstreamTemperature = state.cells[cellIndex - 1].temp;
                double chokeUpstreamVoidFraction = state.cells[cellIndex - 1].alf;
                double betE = state.cells[cellIndex - 1].bet;

                double upstreamLiquidDensity = state.cells[cellIndex - 1].flui.MasEspLiq(state.cells[cellIndex - 1].pres, state.cells[cellIndex - 1].temp);
                double rholc = state.cells[cellIndex - 1].fluicol.MasEspFlu(state.cells[cellIndex - 1].pres, state.cells[cellIndex - 1].temp);
                double upstreamMixtureLiquidDensity = (1 - betE) * upstreamLiquidDensity + betE * rholc;

                double chokeDownstreamVoidFraction = state.cells[cellIndex].alf;
                double betJ = state.cells[cellIndex].bet;
                double downstreamLiquidDensity = state.cells[cellIndex].flui.MasEspLiq(state.cells[cellIndex].pres, state.cells[cellIndex].temp);
                double rholcJ = state.cells[cellIndex].fluicol.MasEspFlu(state.cells[cellIndex].pres, state.cells[cellIndex].temp);
                double downstreamMixtureLiquidDensity = (1 - betJ) * downstreamLiquidDensity + betJ * rholcJ;

                double upstreamHydrostaticHead = sin(state.cells[cellIndex - 1].duto.teta) * (0.5 * state.cells[cellIndex - 1].dx) * (upstreamMixtureLiquidDensity * (1 - chokeUpstreamVoidFraction) + chokeUpstreamVoidFraction * state.cells[cellIndex - 1].flui.MasEspGas(state.cells[cellIndex - 1].pres, state.cells[cellIndex - 1].temp)) * 9.82 / 98600.;
                double downstreamHydrostaticHead = sin(state.cells[cellIndex].duto.teta) * (0.5 * state.cells[cellIndex].dx) * (downstreamMixtureLiquidDensity * (1 - chokeDownstreamVoidFraction) + chokeDownstreamVoidFraction * state.cells[cellIndex].flui.MasEspGas(state.cells[cellIndex].pres, state.cells[cellIndex].temp)) * 9.82 / 98600.;

                double quality = chokeUpstreamVoidFraction * state.cells[cellIndex - 1].flui.MasEspGas(state.cells[cellIndex - 1].pres, state.cells[cellIndex - 1].temp) / (state.cells[cellIndex - 1].flui.MasEspGas(state.cells[cellIndex - 1].pres, state.cells[cellIndex - 1].temp) * chokeUpstreamVoidFraction + upstreamMixtureLiquidDensity * (1. - chokeUpstreamVoidFraction));

                double upstreamLiquidJouleThomson = (1. - betE) * state.cells[cellIndex - 1].flui.JTL(state.cells[cellIndex - 1].pres - upstreamHydrostaticHead, state.cells[cellIndex - 1].temp) - betE / rholcJ;
                double upstreamGasJouleThomson = state.cells[cellIndex - 1].flui.JTG(state.cells[cellIndex - 1].pres - upstreamHydrostaticHead, state.cells[cellIndex - 1].temp);
                sourceTemperature = upstreamTemperature + ((1. - quality) * upstreamLiquidJouleThomson + quality * upstreamGasJouleThomson) * (state.cells[cellIndex].pres + downstreamHydrostaticHead - state.cells[cellIndex - 1].pres - upstreamHydrostaticHead);

                sourceGasSpecificHeat = state.cells[cellIndex - 1].flui.CalorGas(state.cells[cellIndex - 1].pres, state.cells[cellIndex - 1].temp);
                sourceLiquidSpecificHeat = (1 - state.cells[cellIndex - 1].bet) * state.cells[cellIndex - 1].flui.CalorLiq(state.cells[cellIndex - 1].pres, state.cells[cellIndex - 1].temp) + state.cells[cellIndex - 1].bet * state.cells[cellIndex - 1].fluicol.CalorLiq(state.cells[cellIndex - 1].pres, state.cells[cellIndex - 1].temp);
                sourceSpecificHeatRatio = state.cells[cellIndex - 1].flui.ConstAdG(state.cells[cellIndex - 1].pres, state.cells[cellIndex - 1].temp);

            } else {
                sourceGasSpecificHeat = state.cells[cellIndex - 1].flui.CalorGas(state.cells[cellIndex - 1].pres, state.cells[cellIndex - 1].temp);
                sourceLiquidSpecificHeat = (1 - state.cells[cellIndex - 1].bet) * state.cells[cellIndex - 1].flui.CalorLiq(state.cells[cellIndex - 1].pres, state.cells[cellIndex - 1].temp) + state.cells[cellIndex - 1].bet * state.cells[cellIndex - 1].fluicol.CalorLiq(state.cells[cellIndex - 1].pres, state.cells[cellIndex - 1].temp);
                sourceTemperature = state.cells[cellIndex].temp;
                sourceSpecificHeatRatio = state.cells[cellIndex - 1].flui.ConstAdG(state.cells[cellIndex - 1].pres, state.cells[cellIndex - 1].temp);
            }
        } else if ((*state.cells[cellIndex].acsrL).tipo == 8) {
            double pumpUpstreamVoidFraction = state.cells[cellIndex - 1].alf;
            double betM = state.cells[cellIndex - 1].bet;
            sourceGasSpecificHeat = state.cells[cellIndex - 1].flui.CalorGas(state.cells[cellIndex - 1].pres, state.cells[cellIndex - 1].temp);
            sourceLiquidSpecificHeat = (1. - betM) * state.cells[cellIndex - 1].flui.CalorLiq(state.cells[cellIndex - 1].pres, state.cells[cellIndex - 1].temp) + betM * state.cells[cellIndex - 1].fluicol.CalorLiq(state.cells[cellIndex - 1].pres, state.cells[cellIndex - 1].temp);
            double n = (*state.cells[cellIndex].acsrL).bvol.npoli;
            double pumpPressureRatio = (state.cells[cellIndex].pres) / (state.cells[cellIndex - 1].pres);
            sourceTemperature = state.cells[cellIndex - 1].temp * pow(pumpPressureRatio, (n - 1) / n);
        } else {
            sourceGasSpecificHeat = 0.;
            sourceSpecificHeatRatio = 1.;
            sourceLiquidSpecificHeat = 0.;
        }
    } else {
        sourceGasSpecificHeat = 0.;
        sourceSpecificHeatRatio = 1.;
        sourceLiquidSpecificHeat = 0.;
    }

    liquidMassSourceTerm = 0;
    if (state.cells[cellIndex].fontemassLR > 0.)
        liquidMassSourceTerm = state.cells[cellIndex].fontemassLR / cellLength;
    if (state.cells[cellIndex].fontemassCR > 0.)
        liquidMassSourceTerm += state.cells[cellIndex].fontemassCR / cellLength;
    liquidMassSourceTerm *= (sourceLiquidSpecificHeat) * (sourceTemperature - state.cells[cellIndex].temp);

    gasMassSourceTerm = state.cells[cellIndex].fontemassGR / cellLength;
    if (gasMassSourceTerm > 0.)
        gasMassSourceTerm *= (sourceGasSpecificHeat / sourceSpecificHeatRatio) * (sourceTemperature - state.cells[cellIndex].temp);
    else
        gasMassSourceTerm = 0;

    return TemperatureSourceTerms{
        .gas = gasMassSourceTerm,
        .liquid = liquidMassSourceTerm,
    };
}

}  // namespace

void computeTemperature(const ThermalState &state, int cellIndex, double previousTemperature, int steadyMode) {
    if (state.thermalSourceDisabled == 0) {
        TemperatureBalance balance = prepareTemperatureBalance(state, cellIndex, steadyMode);
        double pressureGradient;
        if (cellIndex < state.lastCell - 1)
            pressureGradient = 2. * (state.cells[cellIndex + 1].presaux - state.cells[cellIndex].pres) * 98066.5 / state.cells[cellIndex].dx;
        else if (cellIndex == state.lastCell - 1 && state.surfaceChoke.AreaGarg > 0.5 * state.cells[state.lastCell - 1].duto.area)
            pressureGradient = 2. * (state.cells[cellIndex + 1].presaux - state.cells[cellIndex].pres) * 98066.5 / state.cells[cellIndex].dx;
        else
            pressureGradient = 2. * (state.cells[cellIndex].pres - state.cells[cellIndex].presaux) * 98066.5 / state.cells[cellIndex].dx;
        if (cellIndex > 0 && state.cells[cellIndex - 1].acsr.tipo == 5 &&
            state.cells[cellIndex - 1].acsr.chk.AreaGarg <= (1e-3 + state.input.master1.razareaativ) * state.cells[cellIndex].duto.area)
            pressureGradient = 2. * (state.cells[cellIndex + 1].presaux - state.cells[cellIndex].pres) * 98066.5 / state.cells[cellIndex].dx;
        else if (cellIndex > 0 && state.cells[cellIndex].acsr.tipo == 5 &&
                 state.cells[cellIndex].acsr.chk.AreaGarg <= (1e-3 + state.input.master1.razareaativ) * state.cells[cellIndex].duto.area) {
            pressureGradient = (state.cells[cellIndex].pres - state.cells[cellIndex - 1].pres) * 98066.5 / balance.meanCellLength;
        } else if (cellIndex == 0)
            pressureGradient = 2. * (state.cells[cellIndex + 1].presaux - state.cells[cellIndex].pres) * 98066.5 / state.cells[cellIndex].dx;

        state.cells[cellIndex].VTemper = balance.temperatureSpatialCoefficient / balance.timeCoefficient;
        if ((cellIndex == 0 && state.cells[cellIndex].VTemper < 0.) || (((cellIndex < state.lastCell || state.cells[cellIndex].VTemper >= 0.) && cellIndex > 0) ||
                                                   (cellIndex == state.lastCell && state.surfaceChokeMassCondition == 1) || (cellIndex == state.lastCell && state.input.chkv == 1))) {
            double temperatureGradient = 0.;
            if (cellIndex > 0)
                temperatureGradient = (state.cells[cellIndex].temp - state.cells[cellIndex - 1].tempini) / balance.meanCellLength;

            if (cellIndex < state.lastCell)
                if (state.cells[cellIndex].VTemper < 0)
                    temperatureGradient = (state.cells[cellIndex + 1].tempini - state.cells[cellIndex].temp) / (0.5 * (state.cells[cellIndex + 1].dx + state.cells[cellIndex].dx));
            if (state.cells[cellIndex].acsr.tipo == 5 && state.cells[cellIndex].acsr.chk.AreaGarg <= (1e-3 + state.input.master1.razareaativ) * state.cells[cellIndex].duto.area && state.cells[cellIndex].VTemper <= 0)
                temperatureGradient = 0.;
            if (cellIndex > 0 && state.cells[cellIndex - 1].acsr.tipo == 5 &&
                state.cells[cellIndex - 1].acsr.chk.AreaGarg <= (1e-3 + state.input.master1.razareaativ) * state.cells[cellIndex - 1].duto.area && state.cells[cellIndex].VTemper >= 0)
                temperatureGradient = 0 * (state.cells[cellIndex + 1].tempini - state.cells[cellIndex].temp) / (0.5 * (state.cells[cellIndex + 1].dx + state.cells[cellIndex].dx));
            if ((cellIndex <= 1 && state.cells[cellIndex].VTemper <= 0) || (cellIndex == state.lastCell && state.cells[cellIndex].VTemper <= 0))
                temperatureGradient = 0.;
            if (state.cells[cellIndex].acsr.tipo == 8 && state.cells[cellIndex].acsr.bvol.freq > 1.) {
                pressureGradient = (state.cells[cellIndex].pres - state.cells[cellIndex - 1].pres) * 98066.5 / balance.meanCellLength;
            }

            double kineticTerm = computeKineticTemperatureTerm(state, cellIndex, balance);

            TemperatureSourceTerms sourceTerms =
                computeTemperatureSourceTerms(state, cellIndex, balance.cellLength);

            double latentHeatTerm;
            double phaseChangeMassRate = fabs(state.cells[cellIndex].FonteMudaFase);
            double phaseChangeSign = 1.;
            if (phaseChangeMassRate > 1e-25)
                phaseChangeSign = state.cells[cellIndex].FonteMudaFase / phaseChangeMassRate;
            if (state.input.limTransMass < phaseChangeMassRate)
                phaseChangeMassRate = phaseChangeSign * state.input.limTransMass;
            else
                phaseChangeMassRate *= phaseChangeSign;
            if (state.latentHeatEnabled > 0 && state.input.flashCompleto == 0) {
                latentHeatTerm = interpolateLatentHeat(state, state.cells[cellIndex].presini, state.cells[cellIndex].temp) * phaseChangeMassRate;
            } else if ((state.input.flashCompleto == 1 || state.input.flashCompleto == 2) && state.latentHeatEnabled > 0)
                latentHeatTerm = (state.cells[cellIndex].flui.EntalpGas(state.cells[cellIndex].presini, state.cells[cellIndex].temp) -
                           state.cells[cellIndex].flui.EntalpLiq(state.cells[cellIndex].presini, state.cells[cellIndex].temp)) * phaseChangeMassRate;
            else
                latentHeatTerm = 0;

            if (state.input.latente == 0)
                latentHeatTerm = 0.;
            else if (state.input.condlatente == 0 && latentHeatTerm < 0)
                latentHeatTerm = 0.;

            double interfaceVoidFraction;
            double leftInterfaceVoidFraction;
            if (balance.gasSuperficialVelocity > 0) {
                interfaceVoidFraction = state.cells[cellIndex].alf;
                leftInterfaceVoidFraction = state.cells[cellIndex].alfL;
            } else {
                interfaceVoidFraction = state.cells[cellIndex].alfR;
                leftInterfaceVoidFraction = state.cells[cellIndex].alf;
            }
            double slipVelocity;
            if (interfaceVoidFraction > (*state.globals).localtiny && interfaceVoidFraction < (1. - (*state.globals).localtiny))
                slipVelocity = balance.gasSuperficialVelocity / interfaceVoidFraction - balance.liquidSuperficialVelocity / (1. - interfaceVoidFraction);
            else if (interfaceVoidFraction > (*state.globals).localtiny)
                slipVelocity = balance.gasSuperficialVelocity;
            else
                slipVelocity = balance.liquidSuperficialVelocity;
            double interfacialWorkTerm = balance.flowArea * state.cells[cellIndex].pres * 98066.5 * slipVelocity * (interfaceVoidFraction - leftInterfaceVoidFraction) / state.cells[cellIndex].dx;

            double upstreamThermalPower = 0.;
            if (cellIndex > 0)
                upstreamThermalPower = state.cells[cellIndex - 1].potTermo;
            state.cells[cellIndex].temp = ((balance.timeCoefficient / state.cells[cellIndex].dt) * state.cells[cellIndex].temp - (-balance.pressureTimeCoefficient * (state.cells[cellIndex].pres - state.cells[cellIndex].presini) * 98066.5 / state.cells[cellIndex].dt) + state.cells[cellIndex].dTdLCor * (-balance.temperatureSpatialCoefficient * temperatureGradient + balance.pressureSpatialCoefficient * pressureGradient - kineticTerm - (balance.hydrostaticPower - 0. * interfacialWorkTerm) + (upstreamThermalPower + state.cells[cellIndex].fonteCal) / balance.meanCellLength + sourceTerms.liquid + sourceTerms.gas + balance.heatFlux - latentHeatTerm) - (balance.rc - balance.rp) * (1 - balance.voidFraction) * state.cells[cellIndex].pres * 98066.5 * balance.flowArea * (state.cells[cellIndex].bet - state.cells[cellIndex].betini) / (balance.liquidDensity * state.cells[cellIndex].dt)) / (balance.timeCoefficient / state.cells[cellIndex].dt);

            if ((cellIndex < 148 && cellIndex > 144) && (*state.globals).lixo5 > 53879) {
                int debugStop;
                debugStop = 0.;
            }

            if (fabs(state.cells[cellIndex].temp - state.cells[cellIndex].tempini) / state.cells[cellIndex].dt > 10.) {
                state.cells[cellIndex].temp = state.cells[cellIndex].tempini +
                                 (fabs(state.cells[cellIndex].temp - state.cells[cellIndex].tempini) / (state.cells[cellIndex].temp - state.cells[cellIndex].tempini)) * 10 * state.cells[cellIndex].dt;
            } else if (fabs(state.cells[cellIndex].temp - state.cells[cellIndex].tempini) / state.cells[cellIndex].dt > 1. && fabs(state.cells[cellIndex].VTemper) > 10 * balance.referenceMixtureVelocity)
                state.cells[cellIndex].temp = state.cells[cellIndex].tempini;
            double minimumTemperature = -50.;
            if (state.input.usaTabela == 1)
                minimumTemperature = state.input.tabent.tmin + 1.;
            if (state.cells[cellIndex].temp < minimumTemperature)
                state.cells[cellIndex].temp = minimumTemperature;
            if (state.cells[cellIndex].temp > 200.)
                state.cells[cellIndex].temp = 200.;
        } else if (cellIndex == state.lastCell) {
            if ((*state.globals).chaverede == 0 || state.networkEndpoint == 1 || (*state.globals).chaveRedeParalela == 1)
                state.cells[cellIndex].temp = state.cells[cellIndex].calor.Textern1;
            else
                state.cells[cellIndex].temp = state.gasSurfaceTemperature;
        }
    } else {
        state.cells[cellIndex].temp = state.cells[cellIndex].calor.Textern1;
    }
}

void computeThermalMassTransfer(const ThermalState &state, int cellIndex) {

    double cellLength = state.cells[cellIndex].dx;
    double meanCellLength = 0.5 * (state.cells[cellIndex].dx + state.cells[cellIndex - 1].dx);
    double diameter = state.cells[cellIndex].duto.a;
    double flowArea = 0.25 * M_PI * diameter * diameter;
    double meanVoidFraction = state.cells[cellIndex].alf;
    double betmed = state.cells[cellIndex].bet;
    double meanSuperficialGasVelocity;

    if (state.cells[cellIndex].alf > (*state.globals).localtiny)
        meanSuperficialGasVelocity = state.cells[cellIndex].QG / flowArea;
    else {
        meanSuperficialGasVelocity = 0.;
    }
    double meanSuperficialLiquidVelocity;
    if (state.cells[cellIndex].alf < 1. - (*state.globals).localtiny)
        meanSuperficialLiquidVelocity = state.cells[cellIndex].QL / flowArea;
    else {
        meanSuperficialLiquidVelocity = 0.;
    }
    double rp = state.cells[cellIndex].flui.MasEspLiq(state.cells[cellIndex].pres, state.cells[cellIndex].temp);
    double rc = state.cells[cellIndex].fluicol.MasEspFlu(state.cells[cellIndex].pres, state.cells[cellIndex].temp);
    double liquidDensity = (1. - betmed) * rp + betmed * rc;
    double gasDensity = state.cells[cellIndex].flui.MasEspGas(state.cells[cellIndex].pres, state.cells[cellIndex].temp);
    double liquidSpecificHeat = (1. - betmed) * state.cells[cellIndex].flui.CalorLiq(state.cells[cellIndex].pres, state.cells[cellIndex].temp) + betmed * state.cells[cellIndex].fluicol.CalorLiq(state.cells[cellIndex].pres, state.cells[cellIndex].temp);
    double liquidSpecificHeatConstantVolume = liquidSpecificHeat;
    double gasSpecificHeat = state.cells[cellIndex].flui.CalorGas(state.cells[cellIndex].pres, state.cells[cellIndex].temp);
    double gasSpecificHeatConstantVolume = state.cells[cellIndex].flui.CalorGasVolMod(state.cells[cellIndex].pres, state.cells[cellIndex].temp);
    double liquidJouleThomson = (1. - betmed) * state.cells[cellIndex].flui.JTL(state.cells[cellIndex].pres, state.cells[cellIndex].temp) - betmed / rc;
    double gasJouleThomson = state.cells[cellIndex].flui.JTG(state.cells[cellIndex].pres, state.cells[cellIndex].temp);
    double hydrostaticPower = (liquidDensity * meanSuperficialLiquidVelocity + gasDensity * meanSuperficialGasVelocity) * flowArea * 9.82 * sin(state.cells[cellIndex].duto.teta);

    state.cells[cellIndex].calor.Tint = state.cells[cellIndex].temp;
    state.cells[cellIndex].calor.Vint = meanSuperficialGasVelocity + meanSuperficialLiquidVelocity;
    state.cells[cellIndex].calor.dt = state.cells[cellIndex].dt;
    double liquidConductivity = (1. - betmed) * state.cells[cellIndex].flui.CondLiq(state.cells[cellIndex].pres, state.cells[cellIndex].temp) + betmed * state.cells[cellIndex].fluicol.CondLiq(state.cells[cellIndex].pres, state.cells[cellIndex].temp);
    state.cells[cellIndex].calor.kint = liquidConductivity * (1 - meanVoidFraction) + state.cells[cellIndex].flui.CondGas(state.cells[cellIndex].pres, state.cells[cellIndex].temp) * meanVoidFraction;
    state.cells[cellIndex].calor.cpint = liquidSpecificHeat * (1 - meanVoidFraction) + gasSpecificHeat * meanVoidFraction;
    state.cells[cellIndex].calor.rhoint = liquidDensity * (1 - meanVoidFraction) + gasDensity * meanVoidFraction;
    double liquidViscosity = (1. - betmed) * state.cells[cellIndex].flui.ViscOleo(state.cells[cellIndex].pres, state.cells[cellIndex].temp) + betmed * state.cells[cellIndex].fluicol.VisFlu(state.cells[cellIndex].pres, state.cells[cellIndex].temp);
    state.cells[cellIndex].calor.viscint = liquidViscosity * (1 - meanVoidFraction) * 1.e-3 + state.cells[cellIndex].flui.ViscGas(state.cells[cellIndex].pres, state.cells[cellIndex].temp) * meanVoidFraction * 1.e-3;
    double heatFlux = state.cells[cellIndex].calor.transtrans();

    double timeCoefficient = (liquidDensity * (1 - meanVoidFraction) * liquidSpecificHeatConstantVolume + gasDensity * meanVoidFraction * gasSpecificHeatConstantVolume) * flowArea;
    double pressureTimeCoefficient;
    pressureTimeCoefficient = -state.cells[cellIndex].flui.CalorGasPresMod(state.cells[cellIndex].pres, state.cells[cellIndex].temp) * (gasDensity * meanVoidFraction * flowArea);
    double temperatureSpatialCoefficient = (liquidDensity * meanSuperficialLiquidVelocity * liquidSpecificHeat + gasDensity * meanSuperficialGasVelocity * gasSpecificHeat) * flowArea;
    double pressureSpatialCoefficient = (liquidDensity * meanSuperficialLiquidVelocity * liquidJouleThomson + gasDensity * meanSuperficialGasVelocity * gasJouleThomson) * flowArea;
    double pressureGradient;
    if (cellIndex < state.lastCell)
        pressureGradient = 2. * (state.cells[cellIndex + 1].presaux - state.cells[cellIndex].pres) * 98066.5 / state.cells[cellIndex].dx;
    else
        pressureGradient = 2. * (state.cells[cellIndex].pres - state.cells[cellIndex].presaux) * 98066.5 / state.cells[cellIndex].dx;
    if (state.cells[cellIndex].acsr.tipo == 5 && state.cells[cellIndex].acsr.chk.AreaGarg <= (1e-3 + state.input.master1.razareaativ) * state.cells[cellIndex].duto.area)
        pressureGradient = 2. * (state.cells[cellIndex].pres - state.cells[cellIndex].presaux) * 98066.5 / state.cells[cellIndex].dx;
    state.cells[cellIndex].VTemper = temperatureSpatialCoefficient / timeCoefficient;
    double temperatureGradient = (state.cells[cellIndex].temp - state.cells[cellIndex - 1].temp) / meanCellLength;
    if (cellIndex < state.lastCell)
        if (state.cells[cellIndex].VTemper < 0)
            temperatureGradient = (state.cells[cellIndex + 1].temp - state.cells[cellIndex].temp) / (0.5 * (state.cells[cellIndex + 1].dx + state.cells[cellIndex].dx));
    if (state.cells[cellIndex].acsr.tipo == 5 && state.cells[cellIndex].acsr.chk.AreaGarg <= (1e-3 + state.input.master1.razareaativ) * state.cells[cellIndex].duto.area && state.cells[cellIndex].VTemper <= 0)
        temperatureGradient = 0.;
    if (state.cells[cellIndex - 1].acsr.tipo == 5 && state.cells[cellIndex - 1].acsr.chk.AreaGarg <= (1e-3 + state.input.master1.razareaativ) * state.cells[cellIndex - 1].duto.area && state.cells[cellIndex].VTemper >= 0)
        temperatureGradient = 0.;
    if ((cellIndex == 1 && state.cells[cellIndex].VTemper <= 0) || (cellIndex == state.lastCell && state.cells[cellIndex].VTemper <= 0))
        temperatureGradient = 0.;
    if (state.cells[cellIndex].acsr.tipo == 4 && state.cells[cellIndex].acsr.bcs.freqnova > 1.) {
        pressureGradient = 0.;
    }
    if (state.cells[cellIndex].acsr.tipo == 8 && state.cells[cellIndex].acsr.bvol.freq > 1.) {
        pressureGradient = (state.cells[cellIndex].pres - state.cells[cellIndex - 1].pres) * 98066.5 / meanCellLength;
    }
    if (state.cells[cellIndex].acsr.tipo == 17 && state.cells[cellIndex].acsr.multibcs.freqnova > 1.) {
        pressureGradient = 0.;
    }

    double kineticTerm = 0;
    double upstreamMeanGasVelocity = 0;
    double upstreamMeanLiquidVelocity = 0;
    double meanGasVelocity = 0;
    double meanLiquidVelocity = 0.;
    if (state.cells[cellIndex].acsr.tipo == 0 && state.cells[cellIndex - 1].acsr.tipo == 0 && cellIndex > 1) {

        double kineticCellLength = state.cells[cellIndex - 1].dx;
        double upstreamDiameter = state.cells[cellIndex - 1].duto.a;
        double upstreamFlowArea = 0.25 * M_PI * upstreamDiameter * upstreamDiameter;

        double faceVoidFraction;
        if (state.cells[cellIndex].QL > 0)
            faceVoidFraction = state.cells[cellIndex - 1].alf;
        else
            faceVoidFraction = state.cells[cellIndex].alf;
        double upstreamFaceVoidFraction;
        if (state.cells[cellIndex - 1].QL > 0)
            upstreamFaceVoidFraction = state.cells[cellIndex - 2].alf;
        else
            upstreamFaceVoidFraction = state.cells[cellIndex - 1].alf;

        if (upstreamFaceVoidFraction > 1e-3) {
            upstreamMeanGasVelocity = state.cells[cellIndex - 1].QG / (upstreamFlowArea);
            upstreamMeanGasVelocity /= upstreamFaceVoidFraction;
        }

        if (upstreamFaceVoidFraction < 1. - 1e-3) {
            upstreamMeanLiquidVelocity = state.cells[cellIndex - 1].QL / (upstreamFlowArea);
            upstreamMeanLiquidVelocity /= (1. - upstreamFaceVoidFraction);
        }

        if (faceVoidFraction > 1e-3) {
            meanGasVelocity = meanSuperficialGasVelocity;
            meanGasVelocity /= faceVoidFraction;
        }

        if (faceVoidFraction < 1. - 1e-3) {
            meanLiquidVelocity = meanSuperficialLiquidVelocity;
            meanLiquidVelocity /= (1. - faceVoidFraction);
        }

        kineticTerm = (state.cells[cellIndex].MC - state.cells[cellIndex - 1].Mliqini) * meanGasVelocity * (meanGasVelocity - upstreamMeanGasVelocity) / kineticCellLength + state.cells[cellIndex - 1].Mliqini * meanLiquidVelocity * (meanSuperficialLiquidVelocity - upstreamMeanLiquidVelocity) / kineticCellLength;
    }

    TemperatureSourceTerms sourceTerms =
        computeThermalMassTransferSourceTerms(state, cellIndex, cellLength);

    double latentHeatTerm;
    if (state.latentHeatEnabled > 0 && state.input.flashCompleto == 0) {
        latentHeatTerm = interpolateLatentHeat(state, state.cells[cellIndex].pres, state.cells[cellIndex].temp);
    } else if (state.input.flashCompleto == 1)
        latentHeatTerm = (state.cells[cellIndex].flui.EntalpGas(state.cells[cellIndex].pres, state.cells[cellIndex].temp) -
                   state.cells[cellIndex].flui.EntalpLiq(state.cells[cellIndex].pres, state.cells[cellIndex].temp));
    else
        latentHeatTerm = 0;

    double interfaceVoidFraction;
    double leftInterfaceVoidFraction;
    if (meanSuperficialGasVelocity > 0) {
        interfaceVoidFraction = state.cells[cellIndex].alf;
        leftInterfaceVoidFraction = state.cells[cellIndex].alfL;
    } else {
        interfaceVoidFraction = state.cells[cellIndex].alfR;
        leftInterfaceVoidFraction = state.cells[cellIndex].alf;
    }
    double slipVelocity;
    if (interfaceVoidFraction > (*state.globals).localtiny && interfaceVoidFraction < (1. - (*state.globals).localtiny))
        slipVelocity = meanSuperficialGasVelocity / interfaceVoidFraction - meanSuperficialLiquidVelocity / (1. - interfaceVoidFraction);
    else if (interfaceVoidFraction > (*state.globals).localtiny)
        slipVelocity = meanSuperficialGasVelocity;
    else
        slipVelocity = meanSuperficialLiquidVelocity;
    double interfacialWorkTerm = flowArea * state.cells[cellIndex].pres * 98066.5 * slipVelocity * (interfaceVoidFraction - leftInterfaceVoidFraction) / state.cells[cellIndex].dx;

    if (state.cells[cellIndex].temp < 0 || state.cells[cellIndex].temp > 100) {
        int debugStop;
        debugStop = 0.;
    }

    state.cells[cellIndex].FonteMudaFase = (-(timeCoefficient / state.cells[cellIndex].dt) * (state.cells[cellIndex].temp - state.cells[cellIndex].tempini) - (pressureTimeCoefficient * (state.cells[cellIndex].pres - state.cells[cellIndex].presini) * 98066.5 / state.cells[cellIndex].dt) - temperatureSpatialCoefficient * temperatureGradient + pressureSpatialCoefficient * pressureGradient - kineticTerm - (hydrostaticPower - 0. * interfacialWorkTerm) + state.cells[cellIndex - 1].potB / meanCellLength + sourceTerms.liquid + sourceTerms.gas + heatFlux - (rc - rp) * (1 - meanVoidFraction) * state.cells[cellIndex].pres * 98066.5 * flowArea * (state.cells[cellIndex].bet - state.cells[cellIndex].betini) / (liquidDensity * state.cells[cellIndex].dt)); // / (timeCoefficient / state.cells[cellIndex].dt);

    state.cells[cellIndex].FonteMudaFase /= latentHeatTerm;
}

namespace {

void initializeDistributedMassTransferInlet(
    const ThermalState &state, int cellIndex,
    double &previousLiquidDensity, double &previousOilVolumeFactor,
    double &previousSolutionGasRatio,
    double &previousSolutionGasPressureDerivative) {
    state.cells[0].transmassLini = state.cells[0].transmassL;
    state.cells[0].transmassL = 0.;
    state.sourceUpdater(cellIndex);
    state.cells[1].fontemassLLini = state.cells[1].fontemassLL;
    state.cells[1].fontemassLL = state.cells[0].fontemassLR;
    state.cells[1].fontemassCLini = state.cells[1].fontemassCL;
    state.cells[1].fontemassCL = state.cells[0].fontemassCR;
    state.cells[1].fontemassGLini = state.cells[1].fontemassGL;
    state.cells[1].fontemassGL = state.cells[0].fontemassGR;
    previousLiquidDensity =
        (1 - state.cells[cellIndex].bet) *
            state.cells[cellIndex].flui.MasEspLiq(
                state.cells[cellIndex].pres, state.cells[cellIndex].temp) +
        state.cells[cellIndex].bet *
            state.cells[cellIndex].fluicol.MasEspFlu(
                state.cells[cellIndex].pres, state.cells[cellIndex].temp);
    previousOilVolumeFactor = state.cells[cellIndex].flui.BOFunc(
        state.cells[cellIndex].pres, state.cells[cellIndex].temp);
    previousSolutionGasRatio = state.cells[cellIndex].flui.RS(
        state.cells[cellIndex].pres, state.cells[cellIndex].temp);
    if (state.input.flashCompleto != 2 || state.input.miniTabAtraso > 0) {
        double previousUpstreamOilVolumeFactor = state.cells[cellIndex].flui.BOFunc(
            state.cells[cellIndex].pres * 0.999,
            state.cells[cellIndex].temp);
        double previousUpstreamSolutionGasRatio = state.cells[cellIndex].flui.RS(
            state.cells[cellIndex].pres * 0.999,
            state.cells[cellIndex].temp);
        previousSolutionGasPressureDerivative =
            (previousSolutionGasRatio / previousOilVolumeFactor -
             previousUpstreamSolutionGasRatio / previousUpstreamOilVolumeFactor) /
            (state.cells[cellIndex].pres * 0.001);
    } else {
        ProFlu flutemp = state.cells[cellIndex].flui;
        flutemp.atualizaPropComp(
            state.cells[cellIndex].pres * 0.999,
            state.cells[cellIndex].temp, flutemp.dCalculatedBeta,
            flutemp.oCalculatedLiqComposition,
            flutemp.oCalculatedVapComposition, state.input.pocinjec);
        double previousUpstreamOilVolumeFactor = flutemp.BOFunc(
            state.cells[cellIndex].pres * 0.999,
            state.cells[cellIndex].temp);
        double previousUpstreamSolutionGasRatio = flutemp.RS(
            state.cells[cellIndex].pres * 0.999,
            state.cells[cellIndex].temp);
        previousSolutionGasPressureDerivative =
            (previousSolutionGasRatio / previousOilVolumeFactor -
             previousUpstreamSolutionGasRatio / previousUpstreamOilVolumeFactor) /
            (state.cells[cellIndex].pres * 0.001);
    }
    ProFlu flutemp = state.cells[cellIndex].flui;
}

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

DistributedMassTransferProperties prepareDistributedMassTransferProperties(
    const ThermalState &state, int cellIndex, double meanTemperature, ProFlu &flue,
    ProFlu &flud) {
    double downstreamWaterFraction;
    double upstreamWaterFraction;
    double cellOilFormationVolumeFactor = state.cells[cellIndex - 1].flui.BOFunc(
        state.cells[cellIndex - 1].pres, state.cells[cellIndex - 1].temp);
    double cellWaterFormationVolumeFactor = state.cells[cellIndex - 1].flui.BAFunc(
        state.cells[cellIndex - 1].pres, state.cells[cellIndex - 1].temp);
    double cellWaterFraction = state.cells[cellIndex - 1].flui.BSW * cellWaterFormationVolumeFactor /
                 (cellOilFormationVolumeFactor + cellWaterFormationVolumeFactor * state.cells[cellIndex - 1].flui.BSW -
                  state.cells[cellIndex - 1].flui.BSW * cellOilFormationVolumeFactor);
    if (state.cells[cellIndex].Mliqini < 0.) {
        flud = state.cells[cellIndex].flui;
        double downstreamOilFormationVolumeFactor = flud.BOFunc(state.cells[cellIndex].pres, state.cells[cellIndex].temp);
        double downstreamWaterFormationVolumeFactor = flud.BAFunc(state.cells[cellIndex].pres, state.cells[cellIndex].temp);
        downstreamWaterFraction = flud.BSW * downstreamWaterFormationVolumeFactor / (downstreamOilFormationVolumeFactor + downstreamWaterFormationVolumeFactor * flud.BSW - flud.BSW * downstreamOilFormationVolumeFactor);
    } else {
        flud = state.cells[cellIndex - 1].flui;
        double downstreamOilFormationVolumeFactor = flud.BOFunc(
            state.cells[cellIndex - 1].pres, state.cells[cellIndex - 1].temp);
        double downstreamWaterFormationVolumeFactor = flud.BAFunc(
            state.cells[cellIndex - 1].pres, state.cells[cellIndex - 1].temp);
        downstreamWaterFraction = flud.BSW * downstreamWaterFormationVolumeFactor / (downstreamOilFormationVolumeFactor + downstreamWaterFormationVolumeFactor * flud.BSW - flud.BSW * downstreamOilFormationVolumeFactor);
    }
    if (state.cells[cellIndex - 1].Mliqini < 0) {
        flue = state.cells[cellIndex - 1].flui;
        double upstreamOilFormationVolumeFactor = flue.BOFunc(
            state.cells[cellIndex - 1].pres, state.cells[cellIndex - 1].temp);
        double upstreamWaterFormationVolumeFactor = flue.BAFunc(
            state.cells[cellIndex - 1].pres, state.cells[cellIndex - 1].temp);
        upstreamWaterFraction = flue.BSW * upstreamWaterFormationVolumeFactor / (upstreamOilFormationVolumeFactor + upstreamWaterFormationVolumeFactor * flue.BSW - flue.BSW * upstreamOilFormationVolumeFactor);

    } else {
        if (cellIndex > 1) {
            flue = state.cells[cellIndex - 2].flui;
            double upstreamOilFormationVolumeFactor = flue.BOFunc(
                state.cells[cellIndex - 2].pres, state.cells[cellIndex - 2].temp);
            double upstreamWaterFormationVolumeFactor = flue.BAFunc(
                state.cells[cellIndex - 2].pres, state.cells[cellIndex - 2].temp);
            upstreamWaterFraction = flue.BSW * upstreamWaterFormationVolumeFactor /
                  (upstreamOilFormationVolumeFactor + upstreamWaterFormationVolumeFactor * flue.BSW - flue.BSW * upstreamOilFormationVolumeFactor);
        } else {
            flue = state.cells[cellIndex - 1].flui;
            double upstreamOilFormationVolumeFactor = flue.BOFunc(
                state.cells[cellIndex - 1].pres, state.cells[cellIndex - 1].temp);
            double upstreamWaterFormationVolumeFactor = flue.BAFunc(
                state.cells[cellIndex - 1].pres, state.cells[cellIndex - 1].temp);
            upstreamWaterFraction = flue.BSW * upstreamWaterFormationVolumeFactor /
                  (upstreamOilFormationVolumeFactor + upstreamWaterFormationVolumeFactor * flue.BSW - flue.BSW * upstreamOilFormationVolumeFactor);
        }
    }

    double liquidDensity;
    double gasDensity;
    double betI;

    // casoComp

    liquidDensity = flud.MasEspLiq(state.cells[cellIndex].presaux, meanTemperature);
    if (state.cells[cellIndex].Mliqini < 0)
        betI = state.cells[cellIndex].bet; // testeBeta
    else
        betI = state.cells[cellIndex].betL;

    gasDensity = flud.MasEspGas(state.cells[cellIndex].presaux, meanTemperature);

    double betL = state.cells[cellIndex - 1].betL;
    if (state.cells[cellIndex - 1].Mliqini < 0)
        betL = state.cells[cellIndex - 1].bet; // testeBeta

    if (cellIndex > 0)
        betI = state.cells[cellIndex - 1].betPigD;
    if (state.cells[cellIndex].Mliqini < 0)
        betI = state.cells[cellIndex].betPigE; // testeBeta
    if (cellIndex > 1)
        betL = state.cells[cellIndex - 2].betPigD;
    if (state.cells[cellIndex - 1].Mliqini < 0)
        betL = state.cells[cellIndex - 1].betPigE; // testebeta

    double mixtureLiquidDensity = (1 - betI) * liquidDensity +
                  betI * state.cells[cellIndex].fluicol.MasEspFlu(
                             state.cells[cellIndex].presaux, meanTemperature);
    double downstreamOilVolumeFactor;
    double downstreamSolutionGasRatio;
    double previousDownstreamOilVolumeFactor;
    double previousDownstreamSolutionGasRatio;
    if (state.input.flashCompleto != 2 || state.input.miniTabAtraso > 0) {
        downstreamOilVolumeFactor = flud.BOFunc(state.cells[cellIndex].presaux, meanTemperature);
        downstreamSolutionGasRatio = flud.RS(state.cells[cellIndex].presaux, meanTemperature);
        previousDownstreamOilVolumeFactor = flud.BOFunc(state.cells[cellIndex].presaux * 0.999, meanTemperature);
        previousDownstreamSolutionGasRatio = flud.RS(state.cells[cellIndex].presaux * 0.999, meanTemperature);
    } else {
        downstreamOilVolumeFactor = flud.BOFunc(state.cells[cellIndex].presaux, meanTemperature);
        downstreamSolutionGasRatio = flud.RS(state.cells[cellIndex].presaux, meanTemperature);
        flud.atualizaPropComp(
            state.cells[cellIndex].presaux * 0.999, meanTemperature, flud.dCalculatedBeta,
            flud.oCalculatedLiqComposition,
            flud.oCalculatedVapComposition, state.input.pocinjec);
        previousDownstreamOilVolumeFactor = flud.BOFunc(state.cells[cellIndex].presaux * 0.999, meanTemperature);
        previousDownstreamSolutionGasRatio = flud.RS(state.cells[cellIndex].presaux * 0.999, meanTemperature);
    } // casoComp
    double downstreamSolutionGasPressureDerivative =
        (downstreamSolutionGasRatio / downstreamOilVolumeFactor - previousDownstreamSolutionGasRatio / previousDownstreamOilVolumeFactor) / (state.cells[cellIndex].presaux * 0.001);
    double cellOilVolumeFactor;
    double cellSolutionGasRatio;
    double previousCellOilVolumeFactor;
    double previousCellSolutionGasRatio;
    if (state.input.flashCompleto != 2 || state.input.miniTabAtraso > 0) {
        cellOilVolumeFactor = state.cells[cellIndex - 1].flui.BOFunc(
            state.cells[cellIndex - 1].pres, state.cells[cellIndex - 1].temp);
        cellSolutionGasRatio = state.cells[cellIndex - 1].flui.RS(
            state.cells[cellIndex - 1].pres, state.cells[cellIndex - 1].temp);
        previousCellOilVolumeFactor = state.cells[cellIndex - 1].flui.BOFunc(
            state.cells[cellIndex - 1].pres * 0.999, state.cells[cellIndex - 1].temp);
        previousCellSolutionGasRatio = state.cells[cellIndex - 1].flui.RS(
            state.cells[cellIndex - 1].pres * 0.999, state.cells[cellIndex - 1].temp);
    } else {
        cellOilVolumeFactor = state.cells[cellIndex - 1].flui.BOFunc(
            state.cells[cellIndex - 1].pres, state.cells[cellIndex - 1].temp);
        cellSolutionGasRatio = state.cells[cellIndex - 1].flui.RS(
            state.cells[cellIndex - 1].pres, state.cells[cellIndex - 1].temp);
        ProFlu flutemp = state.cells[cellIndex - 1].flui;
        flutemp.atualizaPropComp(
            state.cells[cellIndex - 1].pres * 0.999,
            state.cells[cellIndex - 1].temp, flutemp.dCalculatedBeta,
            flutemp.oCalculatedLiqComposition,
            flutemp.oCalculatedVapComposition, state.input.pocinjec);
        previousCellOilVolumeFactor = flutemp.BOFunc(
            state.cells[cellIndex - 1].pres * 0.999, state.cells[cellIndex - 1].temp);
        previousCellSolutionGasRatio = flutemp.RS(
            state.cells[cellIndex - 1].pres * 0.999, state.cells[cellIndex - 1].temp);
    } // casoComp
    double cellSolutionGasPressureDerivative =
        (cellSolutionGasRatio / cellOilVolumeFactor - previousCellSolutionGasRatio / previousCellOilVolumeFactor) / (state.cells[cellIndex - 1].pres * 0.001);
    double shiftedCellOilVolumeFactor = 0.;
    double shiftedCellSolutionGasRatio = 0.;
    double cellSolutionGasTemperatureDerivative = 0;
    if (state.input.cicloAcopTerm == 1) {
        if (state.input.flashCompleto != 2 ||
            state.input.miniTabAtraso > 0) {
            shiftedCellOilVolumeFactor = state.cells[cellIndex - 1].flui.BOFunc(
                state.cells[cellIndex - 1].pres,
                state.cells[cellIndex - 1].temp * 0.999);
            shiftedCellSolutionGasRatio = state.cells[cellIndex - 1].flui.RS(
                state.cells[cellIndex - 1].pres,
                state.cells[cellIndex - 1].temp * 0.999);
        } else {
            ProFlu flutemp = state.cells[cellIndex - 1].flui;
            flutemp.atualizaPropComp(
                state.cells[cellIndex - 1].pres,
                state.cells[cellIndex - 1].temp * 0.999,
                flutemp.dCalculatedBeta,
                flutemp.oCalculatedLiqComposition,
                flutemp.oCalculatedVapComposition, state.input.pocinjec);
            shiftedCellOilVolumeFactor = flutemp.BOFunc(
                state.cells[cellIndex - 1].pres,
                state.cells[cellIndex - 1].temp * 0.999);
            shiftedCellSolutionGasRatio = flutemp.RS(
                state.cells[cellIndex - 1].pres,
                state.cells[cellIndex - 1].temp * 0.999);
        } // casoComp
        cellSolutionGasTemperatureDerivative = (cellSolutionGasRatio / cellOilVolumeFactor - shiftedCellSolutionGasRatio / shiftedCellOilVolumeFactor) /
                   (state.cells[cellIndex - 1].temp * 0.001);
    }

    return DistributedMassTransferProperties{
        .downstreamWaterFraction = downstreamWaterFraction,
        .upstreamWaterFraction = upstreamWaterFraction,
        .cellWaterFraction = cellWaterFraction,
        .liquidDensity = liquidDensity,
        .gasDensity = gasDensity,
        .downstreamComposition = betI,
        .upstreamComposition = betL,
        .mixtureLiquidDensity = mixtureLiquidDensity,
        .downstreamOilVolumeFactor = downstreamOilVolumeFactor,
        .downstreamSolutionGasRatio = downstreamSolutionGasRatio,
        .downstreamSolutionGasPressureDerivative = downstreamSolutionGasPressureDerivative,
        .cellOilVolumeFactor = cellOilVolumeFactor,
        .cellSolutionGasRatio = cellSolutionGasRatio,
        .cellSolutionGasPressureDerivative = cellSolutionGasPressureDerivative,
        .cellSolutionGasTemperatureDerivative = cellSolutionGasTemperatureDerivative,
    };
}

struct DistributedMassTransferCoefficients {
    double activeDerivative;
    double spatialCoupling;
    double flowArea;
};

DistributedMassTransferCoefficients updateDistributedMassTransferDerivatives(
    const ThermalState &state, int cellIndex, double cellWaterFraction, const ProFlu &flud,
    double cellSolutionGasPressureDerivative, double cellSolutionGasTemperatureDerivative) {
    double activeDerivative = 1.;
    double pressureLimit = 10;
    if (state.completeModel == 1)
        pressureLimit = 0;
    if (state.cells[cellIndex - 1].pres < pressureLimit ||
        state.massTransferModel != 0)
        activeDerivative = 0.;
    double spatialCoupling = 1.;
    if (cellIndex < 2 || cellIndex == state.lastCell)
        spatialCoupling = 0;

    double coefficientFlowArea = state.cells[cellIndex].dutoL.area;
    state.cells[cellIndex - 1].ativaDeri = activeDerivative;
    state.cells[cellIndex - 1].DTransDtp =
        activeDerivative * coefficientFlowArea * (1. - state.cells[cellIndex - 1].alf) *
        (1. - state.cells[cellIndex - 1].bet) * (1. - cellWaterFraction) * flud.Deng *
        1.225 * cellSolutionGasPressureDerivative * (6.29 / 35.31467);
    state.cells[cellIndex].DTransDtpL = state.cells[cellIndex - 1].DTransDtp;
    if (cellIndex == state.lastCell) {
        double cellOilVolumeFactor;
        double cellSolutionGasRatio;
        double previousCellOilVolumeFactor;
        double previousCellSolutionGasRatio;
        if (state.input.flashCompleto != 2 ||
            state.input.miniTabAtraso > 0) {
            cellOilVolumeFactor = state.cells[cellIndex].flui.BOFunc(
                state.cells[cellIndex].pres, state.cells[cellIndex].temp);
            cellSolutionGasRatio = state.cells[cellIndex].flui.RS(
                state.cells[cellIndex].pres, state.cells[cellIndex].temp);
            previousCellOilVolumeFactor = state.cells[cellIndex].flui.BOFunc(
                state.cells[cellIndex].pres * 0.999, state.cells[cellIndex].temp);
            previousCellSolutionGasRatio = state.cells[cellIndex].flui.RS(
                state.cells[cellIndex].pres * 0.999, state.cells[cellIndex].temp);
        } else {
            cellOilVolumeFactor = state.cells[cellIndex].flui.BOFunc(
                state.cells[cellIndex].pres, state.cells[cellIndex].temp);
            cellSolutionGasRatio = state.cells[cellIndex].flui.RS(
                state.cells[cellIndex].pres, state.cells[cellIndex].temp);
            ProFlu flutemp = state.cells[cellIndex].flui;
            flutemp.atualizaPropComp(
                state.cells[cellIndex].pres * 0.999, state.cells[cellIndex].temp,
                flutemp.dCalculatedBeta,
                flutemp.oCalculatedLiqComposition,
                flutemp.oCalculatedVapComposition, state.input.pocinjec);
            previousCellOilVolumeFactor = flutemp.BOFunc(
                state.cells[cellIndex].pres * 0.999, state.cells[cellIndex].temp);
            previousCellSolutionGasRatio = flutemp.RS(
                state.cells[cellIndex].pres * 0.999, state.cells[cellIndex].temp);
        } // casoComp
        double cellSolutionGasPressureDerivative =
            (cellSolutionGasRatio / cellOilVolumeFactor - previousCellSolutionGasRatio / previousCellOilVolumeFactor) / (state.cells[cellIndex].pres * 0.001);
        state.cells[cellIndex].DTransDtp =
            activeDerivative * coefficientFlowArea * (1. - state.cells[cellIndex].alf) *
            (1. - state.cells[cellIndex].bet) * (1. - cellWaterFraction) * flud.Deng *
            1.225 * cellSolutionGasPressureDerivative * (6.29 / 35.31467);
    }
    if (state.input.cicloAcopTerm == 1) {
        state.cells[cellIndex - 1].DTransDtT =
            activeDerivative * coefficientFlowArea * (1. - state.cells[cellIndex - 1].alf) *
            (1. - state.cells[cellIndex - 1].bet) * (1. - cellWaterFraction) * flud.Deng *
            1.225 * cellSolutionGasTemperatureDerivative * (6.29 / 35.31467);
        state.cells[cellIndex].DTransDtTL = state.cells[cellIndex - 1].DTransDtT;
        if (cellIndex == state.lastCell) {
            double cellOilVolumeFactor;
            double cellSolutionGasRatio;
            double previousCellOilVolumeFactor;
            double previousCellSolutionGasRatio;
            if (state.input.flashCompleto != 2 ||
                state.input.miniTabAtraso > 0) {
                cellOilVolumeFactor = state.cells[cellIndex].flui.BOFunc(
                    state.cells[cellIndex].pres, state.cells[cellIndex].temp);
                cellSolutionGasRatio = state.cells[cellIndex].flui.RS(
                    state.cells[cellIndex].pres, state.cells[cellIndex].temp);
                previousCellOilVolumeFactor = state.cells[cellIndex].flui.BOFunc(
                    state.cells[cellIndex].pres * 0.999, state.cells[cellIndex].temp);
                previousCellSolutionGasRatio = state.cells[cellIndex].flui.RS(
                    state.cells[cellIndex].pres * 0.999, state.cells[cellIndex].temp);
            } else {
                cellOilVolumeFactor = state.cells[cellIndex].flui.BOFunc(
                    state.cells[cellIndex].pres, state.cells[cellIndex].temp);
                cellSolutionGasRatio = state.cells[cellIndex].flui.RS(
                    state.cells[cellIndex].pres, state.cells[cellIndex].temp);
                ProFlu flutemp = state.cells[cellIndex].flui;
                flutemp.atualizaPropComp(
                    state.cells[cellIndex].pres * 0.999, state.cells[cellIndex].temp,
                    flutemp.dCalculatedBeta,
                    flutemp.oCalculatedLiqComposition,
                    flutemp.oCalculatedVapComposition,
                    state.input.pocinjec);
                previousCellOilVolumeFactor = flutemp.BOFunc(
                    state.cells[cellIndex].pres * 0.999, state.cells[cellIndex].temp);
                previousCellSolutionGasRatio = flutemp.RS(
                    state.cells[cellIndex].pres * 0.999, state.cells[cellIndex].temp);
            } // casoComp
            double cellSolutionGasPressureDerivative =
                (cellSolutionGasRatio / cellOilVolumeFactor - previousCellSolutionGasRatio / previousCellOilVolumeFactor) /
                (state.cells[cellIndex].pres * 0.001);
            state.cells[cellIndex].DTransDtT =
                activeDerivative * coefficientFlowArea * (1. - state.cells[cellIndex].alf) *
                (1. - state.cells[cellIndex].bet) * (1. - cellWaterFraction) * flud.Deng *
                1.225 * cellSolutionGasTemperatureDerivative * (6.29 / 35.31467);
        }
    }

    return DistributedMassTransferCoefficients{
        .activeDerivative = activeDerivative,
        .spatialCoupling = spatialCoupling,
        .flowArea = coefficientFlowArea,
    };
}

void selectDistributedMassTransferModel(
    const ThermalState &state, int cellIndex, double &meanTemperature, double &leftMeanTemperature,
    double leftAbsoluteSuperficialVelocity) {
    meanTemperature = state.cells[cellIndex - 1].temp;
    if (state.cells[cellIndex].VTemper < 0.)
        meanTemperature = state.cells[cellIndex].temp;
    leftMeanTemperature = state.cells[cellIndex - 1].tempL;
    if (state.cells[cellIndex - 1].VTemper < 0.)
        leftMeanTemperature = state.cells[cellIndex - 1].temp;

    state.cells[cellIndex - 1].TMModel = state.massTransferModel;
    if (state.massTransferModel != 3) {
        if ((((state.cells[cellIndex - 1].alf < 0.001) ||
              (state.cells[cellIndex - 1].alf > 0.999) ||
              (state.cells[cellIndex - 1].bet > 0.999 &&
               state.cells[cellIndex - 1].alf < 0.999)) &&
             leftAbsoluteSuperficialVelocity < 0.1) ||
            state.cells[cellIndex - 1].flui.RGO >= (*state.globals).RGOMax)
            state.cells[cellIndex - 1].TMModel = 3;
        else if (state.cells[cellIndex - 1].estadoPig == 1)
            state.cells[cellIndex - 1].TMModel = 3;
        else if (state.cells[cellIndex - 1].acsr.tipo == 2 ||
                 state.cells[cellIndex - 1].acsr.tipo == 3 ||
                 state.cells[cellIndex - 1].acsr.tipo == 9 ||
                 state.cells[cellIndex - 1].acsr.tipo == 15 ||
                 state.cells[cellIndex - 1].acsr.tipo == 16)
            state.cells[cellIndex - 1].TMModel = 3;
        else if (cellIndex >= 2) {
            if (state.cells[cellIndex - 2].acsr.tipo == 5 &&
                (state.cells[cellIndex - 2].acsr.chk.AreaGarg <
                 (1e-3 + state.input.master1.razareaativ) *
                     state.cells[cellIndex - 2].duto.area))
                state.cells[cellIndex - 1].TMModel = 3;
            else if (state.cells[cellIndex - 2].acsr.tipo == 4 ||
                     state.cells[cellIndex - 2].acsr.tipo == 7 ||
                     state.cells[cellIndex - 2].acsr.tipo == 17)
                state.cells[cellIndex - 1].TMModel = 0;
        }
        if (state.input.flashCompleto == 2) {
            double candidateMassFraction = state.cells[cellIndex - 1].flui.FracMass(
                state.cells[cellIndex - 1].pres, state.cells[cellIndex - 1].temp);
            if (candidateMassFraction > 1.0 - 1e-2 || candidateMassFraction < 1e-2)
                state.cells[cellIndex - 1].TMModel = 3;
        }
    }
    if (state.cells[cellIndex - 1].TMModel == 0 &&
        state.cells[cellIndex - 1].alf <= (*state.globals).CritCond)
        state.cells[cellIndex - 1].TMModel = 1;
    if (state.cells[cellIndex - 1].TMModel == 0 && cellIndex == state.lastCell)
        state.cells[cellIndex - 1].TMModel = 1;
    if (cellIndex == state.lastCell)
        state.cells[cellIndex].FonteMudaFase = 0.;
}

void applyDistributedMassTransferModel(
    const ThermalState &state, int cellIndex, double &meanTemperature, double &leftMeanTemperature,
    double leftAbsoluteSuperficialVelocity, const ProFlu &flue, const ProFlu &flud, double downstreamWaterFraction,
    double upstreamWaterFraction, double cellWaterFraction, double betI, double betL, double liquidDensity,
    double downstreamOilVolumeFactor, double downstreamSolutionGasRatio, double downstreamSolutionGasPressureDerivative, double cellOilVolumeFactor, double cellSolutionGasRatio,
    double activeDerivative, double spatialCoupling, double coefficientFlowArea, double previousMixtureLiquidDensity, double upstreamOilVolumeFactor,
    double upstreamSolutionGasRatio, double upstreamSolutionGasPressureDerivative) {
    selectDistributedMassTransferModel(state, cellIndex, meanTemperature, leftMeanTemperature, leftAbsoluteSuperficialVelocity);

    state.cells[cellIndex - 1].fontedissolv = 0.;

    state.cells[cellIndex - 1].transmassRini = state.cells[cellIndex - 1].transmassR;
    state.cells[cellIndex - 1].FonteMudaFaseini =
        state.cells[cellIndex - 1].FonteMudaFase;
    state.cells[cellIndex].DTransDxRini = state.cells[cellIndex].DTransDxR;
    state.cells[cellIndex].DTransDxLini = state.cells[cellIndex].DTransDxL;
    state.cells[cellIndex].DTransDt1ini = state.cells[cellIndex].DTransDt1;
    state.cells[cellIndex].DTransDt0ini = state.cells[cellIndex].DTransDt0;
    state.cells[cellIndex].DTransDxRpini = state.cells[cellIndex].DTransDxRp;
    state.cells[cellIndex].DTransDxLpini = state.cells[cellIndex].DTransDxLp;
    state.cells[cellIndex - 1].CoefDTLini = state.cells[cellIndex - 1].CoefDTL;
    state.cells[cellIndex - 1].coefTransBetini =
        state.cells[cellIndex - 1].coefTransBet;
    state.cells[cellIndex].transmassLini = state.cells[cellIndex].transmassL;

    state.cells[cellIndex].TMModelL = state.cells[cellIndex - 1].TMModel;
    state.cells[cellIndex - 1].FonteMudaFase = 0.;
    if (state.cells[cellIndex - 1].TMModel == 0 ||
        state.cells[cellIndex - 1].TMModel == 1) {

        state.cells[cellIndex - 1].transmassR =
            -(state.cells[cellIndex].QL * (1 - betI) * (flud.rDgD) * flud.Deng *
              1.225 * (1. - downstreamWaterFraction) * downstreamSolutionGasRatio * (6.29 / 35.31467) / downstreamOilVolumeFactor) +
            (state.cells[cellIndex - 1].QL * (1 - betL) * (flue.rDgD) *
             flue.Deng * 1.225 * (1. - upstreamWaterFraction) * upstreamSolutionGasRatio * (6.29 / 35.31467) /
             upstreamOilVolumeFactor);

        state.cells[cellIndex - 1].transmassR /= state.cells[cellIndex - 1].dx;
        state.cells[cellIndex - 1].transmassR += state.cells[cellIndex - 1].fontedissolv;
        state.cells[cellIndex].transmassL = state.cells[cellIndex - 1].transmassR;
        state.cells[cellIndex - 1].FonteMudaFase =
            state.cells[cellIndex - 1].transmassR -
            state.cells[cellIndex - 1].DTransDtp * state.cells[cellIndex - 1].d2pdt2 -
            state.cells[cellIndex - 1].DTransDtT * state.cells[cellIndex - 1].dTdtIni;
        if (state.cells[cellIndex - 1].TMModel == 1) {
            state.cells[cellIndex - 1].transmassR -=
                activeDerivative * ((1. - state.cells[cellIndex - 1].bet) *
                         (1. - state.cells[cellIndex - 1].alf) * (1. - cellWaterFraction) * coefficientFlowArea *
                         (state.cells[cellIndex - 1].flui.rDgD) *
                         state.cells[cellIndex - 1].flui.Deng * 1.225 * cellSolutionGasRatio *
                         (6.29 / 35.31467) / cellOilVolumeFactor) /
                state.cells[cellIndex - 1].dt;
            state.cells[cellIndex - 1].transmassR +=
                activeDerivative * ((1. - state.cells[cellIndex - 1].betini) *
                         (1. - state.cells[cellIndex - 1].alfini) * (1. - cellWaterFraction) * coefficientFlowArea *
                         (state.cells[cellIndex - 1].flui.rDgD) *
                         state.cells[cellIndex - 1].flui.Deng * 1.225 * cellSolutionGasRatio *
                         (6.29 / 35.31467) / cellOilVolumeFactor) /
                state.cells[cellIndex - 1].dt;

            state.cells[cellIndex].transmassL = state.cells[cellIndex - 1].transmassR;
            state.cells[cellIndex - 1].FonteMudaFase =
                state.cells[cellIndex - 1].transmassR;
        }
        if (state.cells[cellIndex - 1].TMModel == 0) {
            state.cells[cellIndex - 1].FonteMudaFase -=
                activeDerivative * ((1. - state.cells[cellIndex - 1].bet) *
                         (1. - state.cells[cellIndex - 1].alf) * (1. - cellWaterFraction) * coefficientFlowArea *
                         (state.cells[cellIndex - 1].flui.rDgD) *
                         state.cells[cellIndex - 1].flui.Deng * 1.225 * cellSolutionGasRatio *
                         (6.29 / 35.31467) / cellOilVolumeFactor) /
                state.cells[cellIndex - 1].dt;
            state.cells[cellIndex - 1].FonteMudaFase +=
                activeDerivative * ((1. - state.cells[cellIndex - 1].betini) *
                         (1. - state.cells[cellIndex - 1].alfini) * (1. - cellWaterFraction) * coefficientFlowArea *
                         (state.cells[cellIndex - 1].flui.rDgD) *
                         state.cells[cellIndex - 1].flui.Deng * 1.225 * cellSolutionGasRatio *
                         (6.29 / 35.31467) / cellOilVolumeFactor) /
                state.cells[cellIndex - 1].dt;
        }

        if (state.cells[cellIndex - 1].TMModel == 0) {
            state.cells[cellIndex].DTransDxR =
                -((1 - betI) * (flud.rDgD) * flud.Deng * 1.225 *
                  (1. - downstreamWaterFraction) * downstreamSolutionGasRatio * (6.29 / 35.31467) / downstreamOilVolumeFactor) /
                (liquidDensity * state.cells[cellIndex - 1].dx);
            state.cells[cellIndex].DtransDxLinear =
                -spatialCoupling * state.cells[cellIndex].QL *
                    ((1 - betI) * (flud.rDgD) * flud.Deng * 1.225 *
                     (1. - downstreamWaterFraction) * (6.29 / 35.31467) *
                     (downstreamSolutionGasPressureDerivative * state.cells[cellIndex].dpresaux)) /
                    (state.cells[cellIndex - 1].dx) +
                spatialCoupling * state.cells[cellIndex].QL *
                    ((1 - betI) * (flud.rDgD) * flud.Deng * 1.225 *
                     (1. - downstreamWaterFraction) * (6.29 / 35.31467) *
                     (downstreamSolutionGasPressureDerivative * state.cells[cellIndex].presaux)) /
                    (state.cells[cellIndex - 1].dx);
            state.cells[cellIndex].DTransDxRp =
                -spatialCoupling * state.cells[cellIndex].QL *
                ((1 - betI) * (flud.rDgD) * flud.Deng * 1.225 *
                 (1. - downstreamWaterFraction) * (6.29 / 35.31467) * 0.5 * downstreamSolutionGasPressureDerivative) /
                (state.cells[cellIndex - 1].dx);
            state.cells[cellIndex].DTransDxL =
                ((1 - betL) * (flue.rDgD) * flue.Deng * 1.225 *
                 (1. - upstreamWaterFraction) * upstreamSolutionGasRatio * (6.29 / 35.31467) / upstreamOilVolumeFactor) /
                (previousMixtureLiquidDensity * state.cells[cellIndex - 1].dx);
            state.cells[cellIndex].DtransDxLinear =
                state.cells[cellIndex].DtransDxLinear +
                spatialCoupling * state.cells[cellIndex - 1].QL *
                    ((1 - betL) * (flue.rDgD) * flue.Deng * 1.225 *
                     (1. - upstreamWaterFraction) * (6.29 / 35.31467) *
                     (upstreamSolutionGasPressureDerivative * state.cells[cellIndex - 1].dpresaux)) /
                    (state.cells[cellIndex - 1].dx) -
                spatialCoupling * state.cells[cellIndex - 1].QL *
                    ((1 - betL) * (flue.rDgD) * flue.Deng * 1.225 *
                     (1. - upstreamWaterFraction) * (6.29 / 35.31467) *
                     (upstreamSolutionGasPressureDerivative * state.cells[cellIndex - 1].presaux)) /
                    (state.cells[cellIndex - 1].dx);
            state.cells[cellIndex].DTransDxLp =
                spatialCoupling * state.cells[cellIndex - 1].QL *
                ((1 - betL) * (flue.rDgD) * flue.Deng * 1.225 *
                 (1. - upstreamWaterFraction) * (6.29 / 35.31467) * 0.5 * upstreamSolutionGasPressureDerivative) /
                (state.cells[cellIndex - 1].dx);
            state.cells[cellIndex].DTransDt1 =
                -activeDerivative * ((1. - cellWaterFraction) * coefficientFlowArea *
                          (state.cells[cellIndex - 1].flui.rDgD) *
                          state.cells[cellIndex - 1].flui.Deng * 1.225 * cellSolutionGasRatio *
                          (6.29 / 35.31467) / cellOilVolumeFactor);
            state.cells[cellIndex].DTransDt0 = -state.cells[cellIndex].DTransDt1;

            state.cells[cellIndex - 1].CoefDTR =
                -((1. - state.cells[cellIndex - 1].bet) * (1. - cellWaterFraction) * coefficientFlowArea *
                  (state.cells[cellIndex - 1].flui.rDgD) *
                  state.cells[cellIndex - 1].flui.Deng * 1.225 * cellSolutionGasRatio *
                  (6.29 / 35.31467) / cellOilVolumeFactor);
            state.cells[cellIndex - 1].CoefDTL = -state.cells[cellIndex - 1].CoefDTR;
            state.cells[cellIndex - 1].coefTransBet =
                ((1. - cellWaterFraction) * coefficientFlowArea * (state.cells[cellIndex - 1].flui.rDgD) *
                 state.cells[cellIndex - 1].flui.Deng * 1.225 * cellSolutionGasRatio *
                 (6.29 / 35.31467) / cellOilVolumeFactor);

            // state.cells[cellIndex-1].DTransDtp=((1. - state.cells[cellIndex - 1].bet) * (1. - state.cells[cellIndex - 1].alf) * (1. - cellWaterFraction)*coefficientFlowArea
            state.cells[cellIndex].transmassL -=
                activeDerivative * ((1. - state.cells[cellIndex - 1].bet) *
                         (1. - state.cells[cellIndex - 1].alf) * (1. - cellWaterFraction) * coefficientFlowArea *
                         (state.cells[cellIndex - 1].flui.rDgD) *
                         state.cells[cellIndex - 1].flui.Deng * 1.225 * cellSolutionGasRatio *
                         (6.29 / 35.31467) / cellOilVolumeFactor) /
                state.cells[cellIndex - 1].dt;
            state.cells[cellIndex].transmassL +=
                activeDerivative * ((1. - state.cells[cellIndex - 1].betini) *
                         (1. - state.cells[cellIndex - 1].alfini) * (1. - cellWaterFraction) * coefficientFlowArea *
                         (state.cells[cellIndex - 1].flui.rDgD) *
                         state.cells[cellIndex - 1].flui.Deng * 1.225 * cellSolutionGasRatio *
                         (6.29 / 35.31467) / cellOilVolumeFactor) /
                state.cells[cellIndex - 1].dt;

        } else {
            state.cells[cellIndex].DTransDxR = 0.;
            state.cells[cellIndex].DTransDxL = 0.;
            state.cells[cellIndex].DTransDt1 = 0.;
            state.cells[cellIndex].DTransDxRp = 0.;
            state.cells[cellIndex].DTransDxLp = 0.;
            if (state.input.desligaDeriTransMassDTemp == 1) {
                state.cells[cellIndex - 1].DTransDtT = 0;
                state.cells[cellIndex].DTransDtTL = 0.;
            }
            state.cells[cellIndex - 1].CoefDTR = 0.;
            state.cells[cellIndex - 1].CoefDTL = 0.;
            state.cells[cellIndex - 1].coefTransBet = 0.;
        }
    }
}

}  // namespace

void updateDistributedMassTransfer(const ThermalState &state) {
    // #pragma omp parallel for num_threads(state.input.nthrd)
    double previousMixtureLiquidDensity = 0.;
    double upstreamOilVolumeFactor = 0.;
    double upstreamSolutionGasRatio = 0.;
    double upstreamSolutionGasPressureDerivative = 0.;
    for (int cellIndex = 0; cellIndex <= state.lastCell; cellIndex++) {
        if (cellIndex != 0 && cellIndex != state.lastCell + 1) {

            if ((*state.globals).lixo5 >= 12264.7 && cellIndex == state.lastCell - 1) {
                int debugStop;
                debugStop = 0;
            }

            state.sourceUpdater(cellIndex);

            if (cellIndex < state.lastCell) {
                state.cells[cellIndex + 1].fontemassLLini = state.cells[cellIndex].fontemassLR;
                state.cells[cellIndex + 1].fontemassCLini = state.cells[cellIndex].fontemassCR;
                state.cells[cellIndex + 1].fontemassGLini = state.cells[cellIndex].fontemassGR;
                state.cells[cellIndex + 1].fontemassLL = state.cells[cellIndex].fontemassLR;
                state.cells[cellIndex + 1].fontemassCL = state.cells[cellIndex].fontemassCR;
                state.cells[cellIndex + 1].fontemassGL = state.cells[cellIndex].fontemassGR;
            }
            double lengthRatio = state.cells[cellIndex - 1].dx / (state.cells[cellIndex].dx + state.cells[cellIndex].dxL);
            double leftLengthRatio = state.cells[cellIndex].dx / (state.cells[cellIndex - 1].dx + state.cells[cellIndex - 1].dxL);
            double meanTemperature = lengthRatio * state.cells[cellIndex].temp + (1 - lengthRatio) * state.cells[cellIndex - 1].temp;
            double leftMeanTemperature = leftLengthRatio * state.cells[cellIndex].tempL + (1 - leftLengthRatio) * state.cells[cellIndex - 1].tempL;
            meanTemperature = state.cells[cellIndex - 1].temp;
            if (state.cells[cellIndex].VTemper < 0.)
                meanTemperature = state.cells[cellIndex].temp;
            leftMeanTemperature = state.cells[cellIndex - 1].tempL;
            if (state.cells[cellIndex - 1].VTemper < 0.)
                leftMeanTemperature = state.cells[cellIndex - 1].temp;

            double diameter = state.cells[cellIndex].duto.a;
            double flowArea = 0.25 * M_PI * diameter * diameter;
            double meanSuperficialGasVelocity = (state.cells[cellIndex].QG) / (flowArea);
            double meanSuperficialLiquidVelocity = state.cells[cellIndex].QL / (flowArea);
            double mixtureSuperficialVelocity = meanSuperficialGasVelocity + meanSuperficialLiquidVelocity;
            double leftAbsoluteSuperficialVelocity = (fabs(state.cells[cellIndex - 1].QG) + fabs(state.cells[cellIndex - 1].QL)) / state.cells[cellIndex - 1].duto.area;

            ProFlu upstreamFluid;
            ProFlu downstreamFluid;
            DistributedMassTransferProperties properties =
                prepareDistributedMassTransferProperties(
                    state, cellIndex, meanTemperature, upstreamFluid, downstreamFluid);
            double downstreamWaterFraction = properties.downstreamWaterFraction;
            double upstreamWaterFraction = properties.upstreamWaterFraction;
            double cellWaterFraction = properties.cellWaterFraction;
            double liquidDensity = properties.liquidDensity;
            double gasDensity = properties.gasDensity;
            double betI = properties.downstreamComposition;
            double betL = properties.upstreamComposition;
            double mixtureLiquidDensity = properties.mixtureLiquidDensity;
            double downstreamOilVolumeFactor = properties.downstreamOilVolumeFactor;
            double downstreamSolutionGasRatio = properties.downstreamSolutionGasRatio;
            double downstreamSolutionGasPressureDerivative =
                properties.downstreamSolutionGasPressureDerivative;
            double cellOilVolumeFactor = properties.cellOilVolumeFactor;
            double cellSolutionGasRatio = properties.cellSolutionGasRatio;
            double cellSolutionGasPressureDerivative = properties.cellSolutionGasPressureDerivative;
            double cellSolutionGasTemperatureDerivative =
                properties.cellSolutionGasTemperatureDerivative;
            DistributedMassTransferCoefficients coefficients =
                updateDistributedMassTransferDerivatives(
                    state, cellIndex, cellWaterFraction, downstreamFluid, cellSolutionGasPressureDerivative, cellSolutionGasTemperatureDerivative);
            double activeDerivative = coefficients.activeDerivative;
            double spatialCoupling = coefficients.spatialCoupling;
            double coefficientFlowArea = coefficients.flowArea;

            applyDistributedMassTransferModel(
                state, cellIndex, meanTemperature, leftMeanTemperature, leftAbsoluteSuperficialVelocity, upstreamFluid, downstreamFluid, downstreamWaterFraction, upstreamWaterFraction,
                cellWaterFraction, betI, betL, mixtureLiquidDensity, downstreamOilVolumeFactor, downstreamSolutionGasRatio, downstreamSolutionGasPressureDerivative, cellOilVolumeFactor, cellSolutionGasRatio,
                activeDerivative, spatialCoupling, coefficientFlowArea, previousMixtureLiquidDensity, upstreamOilVolumeFactor, upstreamSolutionGasRatio, upstreamSolutionGasPressureDerivative);
            if (state.cells[cellIndex - 1].TMModel == -2) {
                double hydrateTransportVelocity;
                if (state.cells[cellIndex].alfL > (*state.globals).localtiny && betI < (1. - (*state.globals).localtiny))
                    hydrateTransportVelocity = (state.cells[cellIndex].QG * gasDensity + state.cells[cellIndex].QL * (1. - betI) * liquidDensity) / (coefficientFlowArea * (state.cells[cellIndex].alfL * gasDensity + (1. - state.cells[cellIndex].alfL) * (1. - betI) * liquidDensity));
                else
                    hydrateTransportVelocity = 0.;
                double hydrateMassFraction = state.cells[cellIndex - 1].flui.FracMassHidra(state.cells[cellIndex - 1].pres, state.cells[cellIndex - 1].temp);
                double pressurePerturbationRatio = 0.999;
                double hydrateFractionPressureDerivative = (hydrateMassFraction - state.cells[cellIndex - 1].flui.FracMassHidra(state.cells[cellIndex - 1].pres * pressurePerturbationRatio, state.cells[cellIndex - 1].temp)) / ((1 - pressurePerturbationRatio) * state.cells[cellIndex - 1].pres);
                double pressureGradient;
                pressureGradient = (state.cells[cellIndex].presaux - state.cells[cellIndex - 1].presaux) / state.cells[cellIndex].dxL;
                state.cells[cellIndex].transmassL = state.cells[cellIndex - 1].transmassR = 1 * (state.cells[cellIndex].alfL * gasDensity + (1. - state.cells[cellIndex].alfL) * (1. - betI) * liquidDensity) * hydrateTransportVelocity * (pressureGradient * hydrateFractionPressureDerivative) * coefficientFlowArea;

                state.cells[cellIndex].DTransDxR = 0.;
                state.cells[cellIndex].DTransDxL = 0.;
                state.cells[cellIndex].DTransDt1 = 0.;
                state.cells[cellIndex].DTransDt0 = 0.;
                if (state.input.desligaDeriTransMassDTemp == 1) {
                    state.cells[cellIndex - 1].DTransDtT = 0;
                    state.cells[cellIndex].DTransDtTL = 0.;
                }
                state.cells[cellIndex - 1].CoefDTR = 0.;
                state.cells[cellIndex - 1].CoefDTL = 0.;
                state.cells[cellIndex - 1].coefTransBet = 0.;
            }
            if (state.cells[cellIndex - 1].TMModel == 3) {
                state.cells[cellIndex].transmassL = state.cells[cellIndex - 1].transmassR = 0.;
                state.cells[cellIndex].DTransDxR = 0.;
                state.cells[cellIndex].DTransDxL = 0.;
                state.cells[cellIndex].DTransDt1 = 0.;
                state.cells[cellIndex].DTransDt0 = 0.;
                state.cells[cellIndex].DTransDxRp = 0.;
                state.cells[cellIndex].DTransDxLp = 0.;
                state.cells[cellIndex - 1].CoefDTR = 0.;
                state.cells[cellIndex - 1].CoefDTL = 0.;
                state.cells[cellIndex - 1].coefTransBet = 0.;
                if (state.input.desligaDeriTransMassDTemp == 1) {
                    state.cells[cellIndex - 1].DTransDtT = 0;
                    state.cells[cellIndex].DTransDtTL = 0.;
                }
            }
            if (state.cells[cellIndex - 1].transmassR > 0 && (state.cells[cellIndex].alfL > (1. - (*state.globals).localtiny) || state.cells[cellIndex].betL > (1. - (*state.globals).localtiny))) {
                state.cells[cellIndex].transmassL = state.cells[cellIndex - 1].transmassR = -(*state.globals).localtiny;
                state.cells[cellIndex].DTransDxR = 0.;
                state.cells[cellIndex].DTransDxL = 0.;
                state.cells[cellIndex].DTransDt1 = 0.;
                state.cells[cellIndex].DTransDt0 = 0.;
                state.cells[cellIndex].DTransDxRp = 0.;
                state.cells[cellIndex].DTransDxLp = 0.;
                state.cells[cellIndex - 1].CoefDTR = 0.;
                state.cells[cellIndex - 1].CoefDTL = 0.;
                state.cells[cellIndex - 1].coefTransBet = 0.;
                if (state.input.desligaDeriTransMassDTemp == 1) {
                    state.cells[cellIndex - 1].DTransDtT = 0;
                    state.cells[cellIndex].DTransDtTL = 0.;
                }
            }
            if (state.cells[cellIndex - 1].transmassR < 0 && state.cells[cellIndex].alfL < (*state.globals).localtiny) {
                state.cells[cellIndex].transmassL = state.cells[cellIndex - 1].transmassR = (*state.globals).localtiny;
                state.cells[cellIndex].DTransDxR = 0.;
                state.cells[cellIndex].DTransDxL = 0.;
                state.cells[cellIndex].DTransDt1 = 0.;
                state.cells[cellIndex].DTransDt0 = 0.;
                state.cells[cellIndex].DTransDxRp = 0.;
                state.cells[cellIndex].DTransDxLp = 0.;
                state.cells[cellIndex - 1].CoefDTR = 0.;
                state.cells[cellIndex - 1].CoefDTL = 0.;
                state.cells[cellIndex - 1].coefTransBet = 0.;
                if (state.input.desligaDeriTransMassDTemp == 1) {
                    state.cells[cellIndex - 1].DTransDtT = 0;
                    state.cells[cellIndex].DTransDtTL = 0.;
                }
            }
            previousMixtureLiquidDensity = mixtureLiquidDensity;
            upstreamOilVolumeFactor = downstreamOilVolumeFactor;
            upstreamSolutionGasRatio = downstreamSolutionGasRatio;
            upstreamSolutionGasPressureDerivative = downstreamSolutionGasPressureDerivative;

        } else if (cellIndex == 0)
            initializeDistributedMassTransferInlet(
                state, cellIndex, previousMixtureLiquidDensity, upstreamOilVolumeFactor, upstreamSolutionGasRatio, upstreamSolutionGasPressureDerivative);
    }
}

namespace {

void selectAndApplyInteriorFlowRegime(
    const ThermalState &state, int cellIndex, Vcr<int> &bif, double superficialGasVelocity,
    double superficialLiquidVelocity, double leftSuperficialGasVelocity, double leftSuperficialLiquidVelocity, double rightSuperficialGasVelocity, double rightSuperficialLiquidVelocity,
    double rightGasDensity, double rightLiquidDensity, double gasDensity, double liquidDensity, double flowArea) {
    bif[cellIndex] = 1;

    if ((*state.globals).lixo5 > 29900) {
        int debugStop;
        debugStop = 0;
    }

    if (state.cells[cellIndex - 1].alfPigD <= (*state.globals).localtiny && state.cells[cellIndex].alfPigE <= (*state.globals).localtiny && state.cells[cellIndex].fontemassGL <= (*state.globals).localtiny * 1e-5 && state.cells[cellIndex].fontemassGR <= (*state.globals).localtiny * 1e-5) {
        state.cells[cellIndex].term1 = 1.;
        state.cells[cellIndex].term2 = 0.;
        bif[cellIndex] = 0;
        state.cells[cellIndex].c0 = 1;
        state.cells[cellIndex].ud = 0;
        state.cells[cellIndex].arranjo = 0;
    } else if (state.cells[cellIndex - 1].alfPigD >= (1. - (*state.globals).localtiny) && state.cells[cellIndex].alfPigE >= (1. - (*state.globals).localtiny) && ((state.cells[cellIndex].fontemassLL + state.cells[cellIndex].fontemassCL) <= (*state.globals).localtiny * 1e-5 && (state.cells[cellIndex].fontemassLR + state.cells[cellIndex].fontemassCR) <= (*state.globals).localtiny * 1e-5)) {
        state.cells[cellIndex].term1 = 0.;
        state.cells[cellIndex].term2 = 0.;
        bif[cellIndex] = 0;
        state.cells[cellIndex].c0 = 1;
        state.cells[cellIndex].ud = 0;
        state.cells[cellIndex].arranjo = 0;
    } else if (state.cells[cellIndex - 1].acsr.tipo == 5 && state.cells[cellIndex - 1].acsr.chk.AreaGarg <= (1e-3)) {
        state.cells[cellIndex].term1 = 0.;
        state.cells[cellIndex].term2 = 0.;
        bif[cellIndex] = 0;
        state.cells[cellIndex].c0 = 1;
        state.cells[cellIndex].ud = 0;
        state.cells[cellIndex].arranjo = 0;
    }

    else if (superficialGasVelocity >= 0 && state.cells[cellIndex - 1].alfPigD <= (*state.globals).localtiny && (state.cells[cellIndex].fontemassGL <= (*state.globals).localtiny * 1e-5 && state.cells[cellIndex].fontemassGR <= (*state.globals).localtiny * 1e-5)) {
        state.cells[cellIndex].term1 = 1.;
        state.cells[cellIndex].term2 = 0.;
        bif[cellIndex] = 0;
        state.cells[cellIndex].c0 = 1;
        state.cells[cellIndex].ud = 0;
        state.cells[cellIndex].arranjo = 0;
        if (fabs(superficialGasVelocity) <= 1e-15 && state.cells[cellIndex].alfPigE > (1. - (*state.globals).localtiny) && superficialLiquidVelocity < 0 && leftSuperficialLiquidVelocity < 0 && state.cells[cellIndex].duto.teta > 0) {
            state.cells[cellIndex].term1 = 0.;
            state.cells[cellIndex].term2 = 0.;
            bif[cellIndex] = 0;
        }
        if (fabs(superficialGasVelocity) <= 1e-15 && state.cells[cellIndex].alfPigE > (1. - 10 * (*state.globals).localtiny) && state.cells[cellIndex].duto.teta < 0) {
            state.cells[cellIndex].term1 = 0.;
            state.cells[cellIndex].term2 = 0.;
            bif[cellIndex] = 1;
        }
        if (fabs(superficialGasVelocity) <= 1e-15 && state.cells[cellIndex].alf >= state.cells[cellIndex + 1].alf && rightSuperficialGasVelocity < 0) {
            bif[cellIndex] = 1;
        }
    } else if (superficialGasVelocity <= 0 && state.cells[cellIndex].alfPigE <= (*state.globals).localtiny && (state.cells[cellIndex].fontemassGL <= (*state.globals).localtiny * 1e-5 && state.cells[cellIndex].fontemassGR <= (*state.globals).localtiny * 1e-5)) {
        state.cells[cellIndex].term1 = 1.;
        state.cells[cellIndex].term2 = 0.;
        bif[cellIndex] = 0;
        state.cells[cellIndex].c0 = 1;
        state.cells[cellIndex].ud = 0;
        state.cells[cellIndex].arranjo = 0;
        if (fabs(superficialGasVelocity) <= 1e-15 && state.cells[cellIndex - 1].alfPigD > (1. - (*state.globals).localtiny) && superficialLiquidVelocity > 0 && rightSuperficialLiquidVelocity > 0) {
            state.cells[cellIndex].term1 = 0.;
            state.cells[cellIndex].term2 = 0.;
            bif[cellIndex] = 0;
        }
        if (fabs(superficialGasVelocity) <= 1e-15 && state.cells[cellIndex].alfL >= state.cells[cellIndex - 1].alfL && leftSuperficialGasVelocity > 0) {
            bif[cellIndex] = 1;
        }
        if (fabs(superficialGasVelocity) <= 1e-15 && fabs(superficialLiquidVelocity) <= 1e-15 && state.cells[cellIndex].alf < state.cells[cellIndex - 1].alf && state.cells[cellIndex].duto.teta > 0) {
            bif[cellIndex] = 1;
        }
    } else if (superficialLiquidVelocity >= 0 && state.cells[cellIndex - 1].alfPigD >= 1. - 1 * (*state.globals).localtiny && ((state.cells[cellIndex].fontemassLL + state.cells[cellIndex].fontemassCL) <= (*state.globals).localtiny * 1e-5 && (state.cells[cellIndex].fontemassLR + state.cells[cellIndex].fontemassCR) <= (*state.globals).localtiny * 1e-5)) {
        state.cells[cellIndex].term1 = 0.;
        state.cells[cellIndex].term2 = 0.;
        bif[cellIndex] = 0;
        state.cells[cellIndex].c0 = 1;
        state.cells[cellIndex].ud = 0;
        state.cells[cellIndex].arranjo = 0;
        if (fabs(superficialLiquidVelocity) <= 1e-15 && state.cells[cellIndex].alfPigE < (*state.globals).localtiny && rightSuperficialLiquidVelocity < 0) { // ATENCAO!!!!!!!!!!!!!!! não teria de ser bifásico, mono-liq só seo ângulo fosse negativo, não?
            state.cells[cellIndex].term1 = 1.;
            state.cells[cellIndex].term2 = 0.;
            bif[cellIndex] = 0;
        } else if (fabs(rightSuperficialLiquidVelocity) < (*state.globals).localtiny * 1e-5) {
            if (fabs(superficialLiquidVelocity) < (*state.globals).localtiny * 1e-5 && state.cells[cellIndex].alfPigE < (*state.globals).localtiny && state.cells[cellIndex].fontemassGR >= (*state.globals).localtiny * 1e-5) { // ATENCAO!!!!!!!!!!!!!!!  sem sentido isto aqui
                state.cells[cellIndex].term1 = 1.;
                state.cells[cellIndex].term2 = 0.;
                bif[cellIndex] = 0;
            }
            if (fabs(superficialLiquidVelocity) < (*state.globals).localtiny * 1e-5 && state.cells[cellIndex].alfPigE < (*state.globals).localtiny && state.cells[cellIndex].duto.teta >= 0) { // ATENCAO!!!!!!!!!!!!!!! alteracao 11/08/24, adicionado
                bif[cellIndex] = 1;
            }
            if (fabs(superficialLiquidVelocity) < (*state.globals).localtiny * 1e-5 && state.cells[cellIndex].alfPigE < (*state.globals).localtiny && state.cells[cellIndex].duto.teta < 0) { // ATENCAO!!!!!!!!!!!!!!! alteracao 11/08/24, adicionado
                state.cells[cellIndex].term1 = 0.;
                state.cells[cellIndex].term2 = 0.;
                bif[cellIndex] = 0;
            }
        } else if ((fabs(superficialLiquidVelocity) < 1e-15 && (rightSuperficialLiquidVelocity < 0 || state.cells[cellIndex].duto.teta > 0) // ATENCAO!!!!!!!!!!!!!!! alteracao 11/08/24, estava || mudado para &&
                    && ((state.cells[cellIndex].alfPigE <= (1 - 10 * (*state.globals).localtiny + .0 * state.cells[cellIndex].alfPigER) &&
                         state.cells[cellIndex].alfPigER < 1 - 1 * (*state.globals).localtiny) ||
                        state.cells[cellIndex].alfPigE <= 0.7)))
            bif[cellIndex] = 1;

    } else if (superficialLiquidVelocity <= 0 && state.cells[cellIndex].alfPigE >= 1. - (*state.globals).localtiny && ((state.cells[cellIndex].fontemassLL + state.cells[cellIndex].fontemassCL) <= (*state.globals).localtiny * 1e-5 && (state.cells[cellIndex].fontemassLR + state.cells[cellIndex].fontemassCR) <= (*state.globals).localtiny * 1e-5)) {
        state.cells[cellIndex].term1 = 0.;
        state.cells[cellIndex].term2 = 0.;
        bif[cellIndex] = 0;
        state.cells[cellIndex].c0 = 1;
        state.cells[cellIndex].ud = 0;
        state.cells[cellIndex].arranjo = 0;
        if (fabs(superficialLiquidVelocity) <= (*state.globals).localtiny * 1e-5 && state.cells[cellIndex - 1].alfPigD < (*state.globals).localtiny && leftSuperficialLiquidVelocity > 0) {
            state.cells[cellIndex].term1 = 1.;
            state.cells[cellIndex].term2 = 0.;
            bif[cellIndex] = 0;
        }

        if (fabs(leftSuperficialLiquidVelocity) < (*state.globals).localtiny * 1e-5) {
            if (fabs(superficialLiquidVelocity) < (*state.globals).localtiny * 1e-5 && state.cells[cellIndex - 1].alfPigD < (*state.globals).localtiny && state.cells[cellIndex - 1].fontemassGR >= (*state.globals).localtiny * 1e-5) {
                state.cells[cellIndex].term1 = 1.;
                state.cells[cellIndex].term2 = 0.;
                bif[cellIndex] = 0;
            }
        } else if ((cellIndex > 1 && fabs(superficialLiquidVelocity) < (*state.globals).localtiny * 1e-5 && state.cells[cellIndex].duto.teta < 0.95 * M_PI / 2. && leftSuperficialLiquidVelocity > 0 && ((state.cells[cellIndex - 1].alfPigD <= 1.0 * state.cells[cellIndex - 2].alfPigD && state.cells[cellIndex - 2].alfPigD < 0.99) || state.cells[cellIndex - 1].alfPigD < 0.7)) && superficialGasVelocity >= 0.)
            bif[cellIndex] = 1;
        else if ((cellIndex > 1 && fabs(superficialLiquidVelocity) < (*state.globals).localtiny * 1e-5 && state.cells[cellIndex].duto.teta >= 0.95 * M_PI / 2. && leftSuperficialLiquidVelocity > 0 && ((state.cells[cellIndex - 1].alfPigD <= 1.0 * state.cells[cellIndex - 2].alfPigD && state.cells[cellIndex - 2].alfPigD < 0.99) || state.cells[cellIndex - 1].alfPigD < 0.7)) && superficialGasVelocity >= 0)
            bif[cellIndex] = 1;
        else {
            if ((fabs(superficialLiquidVelocity) < (*state.globals).localtiny * 1e-5 && state.cells[cellIndex].duto.teta < 0.95 * M_PI / 2. && leftSuperficialLiquidVelocity > 0 && ((state.cells[cellIndex - 1].alfPigD <= 1.0 * state.cells[cellIndex - 1].alfL) || state.cells[cellIndex - 1].alfPigD < 0.7)) && superficialGasVelocity >= 0)
                bif[cellIndex] = 1;
            else if ((fabs(superficialLiquidVelocity) < (*state.globals).localtiny * 1e-5 && state.cells[cellIndex].duto.teta >= 0.95 * M_PI / 2. && leftSuperficialLiquidVelocity > 0 && ((state.cells[cellIndex - 1].alfPigD <= 1.0 * state.cells[cellIndex - 1].alfL) || state.cells[cellIndex - 1].alfPigD < 0.7)) && superficialGasVelocity >= 0)
                bif[cellIndex] = 1;
            else if ((fabs(superficialLiquidVelocity) < (*state.globals).localtiny * 1e-5 && state.cells[cellIndex].duto.teta >= 0.95 * M_PI / 2. && ((state.cells[cellIndex - 1].alfPigD <= 1.0 * state.cells[cellIndex - 1].alfL) || state.cells[cellIndex - 1].alfPigD < 0.7)) && superficialGasVelocity > 0)
                bif[cellIndex] = 1;
        }
    }

    if (bif[cellIndex] == 1) {

        double c0;
        double ud;
        double meanVoidFraction;

        meanVoidFraction = state.cells[cellIndex - 1].alfPigD;
        if (superficialGasVelocity < 0)
            meanVoidFraction = state.cells[cellIndex].alfPigE;
        c0 = 1.2;
        double meanDiameter = state.cells[cellIndex].duto.a;
        if (state.cells[cellIndex].MC >= 0)
            meanDiameter = state.cells[cellIndex].dutoL.a;
        double inclinationSign = 1.;
        if (state.cells[cellIndex].duto.teta < 0.)
            inclinationSign = -1.;
        ud = inclinationSign * 0.32 * sqrt(9.82 * meanDiameter);
        if (fabs(rightGasDensity) / rightLiquidDensity > 0.9) {
            c0 = 1.;
            ud = 0.;
        }
        if (fabs(state.cells[cellIndex].QG / (0.25 * M_PI * meanDiameter * meanDiameter * meanVoidFraction)) > 100. ||
            fabs(state.cells[cellIndex].QL / (0.25 * M_PI * meanDiameter * meanDiameter * (1. - meanVoidFraction))) > 100.) {
            c0 = 1.;
            ud = 0.;
        } else
            state.closureUpdater.instantaneous(cellIndex, c0, ud);
        state.cells[cellIndex].c0 = c0;
        state.cells[cellIndex].ud = ud;
        if (cellIndex == state.lastCell) {
            double numerator = (1. - meanVoidFraction * c0);
            double denominator = 1 + c0 * meanVoidFraction * ((gasDensity / liquidDensity) - 1.);
            state.cells[cellIndex].term1 = numerator / denominator;
            state.cells[cellIndex].term2 = (-flowArea * meanVoidFraction * gasDensity * ud) / denominator;
            double previousLiquidDriftFlux = (superficialGasVelocity - meanVoidFraction * ud) / (meanVoidFraction * c0) - superficialGasVelocity;
            double liquidDriftFlux = (superficialGasVelocity + superficialLiquidVelocity) * (1. - c0 * meanVoidFraction) - meanVoidFraction * ud;
            if ((liquidDriftFlux > 0. || previousLiquidDriftFlux > 0.) && state.cells[cellIndex - 1].alfPigD > 1 - 1e-15) {
                state.cells[cellIndex].term1 = 0.;
                state.cells[cellIndex].term2 = 0.;
            }
            if ((liquidDriftFlux < 0. || previousLiquidDriftFlux < 0.) && state.cells[cellIndex].alfPigE > 1 - 1e-15) {
                state.cells[cellIndex].term1 = 0.;
                state.cells[cellIndex].term2 = 0.;
            }
            if (state.cells[cellIndex - 1].acsr.tipo == 5 && state.cells[cellIndex - 1].acsr.chk.AreaGarg <= (1e-3 + state.input.master1.razareaativ) * state.cells[cellIndex - 1].duto.area) {
                state.cells[cellIndex].term1 = 0.;
                state.cells[cellIndex].term2 = 0.;
            }
        }
    }
}

void updateInteriorFlowPartitionCell(
    const ThermalState &state, int cellIndex, Vcr<int> &bif, Vcr<int> &valv) {
    if (cellIndex < state.lastCell) {
        state.cells[cellIndex - 1].alfR = state.cells[cellIndex + 1].alfL = state.cells[cellIndex].alf;
        state.cells[cellIndex - 1].betR = state.cells[cellIndex + 1].betL = state.cells[cellIndex].bet;
    } else {
        state.cells[cellIndex - 1].alfR = state.cells[cellIndex].alf;
        state.cells[cellIndex - 1].betR = state.cells[cellIndex].bet;
    }
    valv[cellIndex] = 1;
    if (state.cells[cellIndex - 1].acsr.tipo == 5 || state.cells[cellIndex - 1].acsr.tipo == 8) {
        if ((*state.cells[cellIndex].acsrL).tipo == 5 && (*state.cells[cellIndex].acsrL).chk.AreaGarg <= (1e-3 + state.input.master1.razareaativ) * state.cells[cellIndex - 1].duto.area)
            valv[cellIndex] = 0;
        if ((*state.cells[cellIndex].acsrL).tipo == 8 && fabs((*state.cells[cellIndex].acsrL).bvol.freq) > 1)
            valv[cellIndex] = 0;
    }
    if (valv[cellIndex] == 1) {
        double lengthRatio = state.cells[cellIndex].dxL / (state.cells[cellIndex].dx + state.cells[cellIndex].dxL);
        double meanPressure;
        if (cellIndex < state.lastCell)
            meanPressure = state.cells[cellIndex].presaux;
        else
            meanPressure = state.cells[cellIndex].pres;
        double meanTemperature;
        if (cellIndex < state.lastCell)
            meanTemperature = state.cells[cellIndex].temp * lengthRatio + state.cells[cellIndex - 1].temp * (1. - lengthRatio);
        else
            meanTemperature = state.gasSurfaceTemperature;
        if (state.cells[cellIndex].VTemper < 0.) {
            if (cellIndex < state.lastCell)
                meanTemperature = state.cells[cellIndex].temp;
            else
                meanTemperature = state.gasSurfaceTemperature;
        }
        double betI = state.cells[cellIndex - 1].betPigD;
        double liquidDensity;
        if (state.cells[cellIndex].QL < 0.) { // testeBeta
            betI = state.cells[cellIndex].betPigE;
            liquidDensity = (1 - betI) * state.cells[cellIndex].rpCi + betI * state.cells[cellIndex].rcCi;
        } else {
            betI = state.cells[cellIndex - 1].betPigD;
            liquidDensity = (1 - betI) * state.cells[cellIndex - 1].rpCi + betI * state.cells[cellIndex - 1].rcCi;
            // viscl1 = (1 - betI) * state.cells[cellIndex - 1].flui.ViscOleo(meanPressure, meanTemperature)
            // tensup1 = (1 - betI) * state.cells[cellIndex - 1].flui.TensSuper(meanPressure, meanTemperature)
        }
        double gasDensity;
        double flowArea;
        double noSlipLiquidHoldup;
        if (state.cells[cellIndex].QG >= 0) {
            flowArea = state.cells[cellIndex].dutoL.area;
            noSlipLiquidHoldup = 1. - state.cells[cellIndex].alfL;
            gasDensity = state.cells[cellIndex - 1].rgCi;
        } else {
            gasDensity = state.cells[cellIndex].rgCi;
            flowArea = state.cells[cellIndex].duto.area;
            noSlipLiquidHoldup = 1. - state.cells[cellIndex].alf;
        }
        double superficialGasVelocity = state.cells[cellIndex].QG / (flowArea);
        double superficialLiquidVelocity = state.cells[cellIndex].QL / (flowArea);
        double diameter = state.cells[cellIndex].duto.a;
        if (superficialGasVelocity >= 0)
            diameter = state.cells[cellIndex - 1].duto.a;

        double mixtureDensity = noSlipLiquidHoldup * liquidDensity + (1 - noSlipLiquidHoldup) * gasDensity;
        double inclinationAngle = state.cells[cellIndex].duto.teta;
        if (cellIndex >= 2) {
            if (state.cells[cellIndex - 2].acsr.tipo == 5 && state.cells[cellIndex - 2].acsr.chk.AreaGarg <= (1e-3)) {
                if (state.cells[cellIndex].QG >= 0)
                    inclinationAngle = state.cells[cellIndex].duto.teta;
                else
                    inclinationAngle = state.cells[cellIndex].dutoR.teta;
            } else {
                if (state.cells[cellIndex].QG >= 0)
                    inclinationAngle = state.cells[cellIndex].dutoL.teta;
                else
                    inclinationAngle = state.cells[cellIndex].duto.teta;
            }
        }
        double inclinationSign = 1.;
        if (inclinationAngle < 0.)
            inclinationSign = -1.;

        double leftFlowArea = state.cells[cellIndex].dutoL.area;
        double leftLengthRatio = state.cells[cellIndex - 1].dxL / (state.cells[cellIndex - 1].dx + state.cells[cellIndex - 1].dxL);
        double betIL = 0.;
        if (cellIndex < 2)
            betIL = state.cells[cellIndex - 1].betL;
        else
            betIL = state.cells[cellIndex - 2].betPigD;
        if (state.cells[cellIndex - 1].QL < 0.)
            betIL = state.cells[cellIndex - 1].betPigE; // testeBeta
        // betIL = state.cells[cellIndex - 1].betPigE;        //duvidabeta
        double leftGasDensity = state.cells[cellIndex].rgLi;
        double leftLiquidDensity = (1 - betIL) * state.cells[cellIndex].rpLi + betIL * state.cells[cellIndex].rcLi;
        double leftSuperficialGasVelocity = (state.cells[cellIndex].ML - state.cells[cellIndex].MliqiniL) / (leftGasDensity * leftFlowArea);
        double leftSuperficialLiquidVelocity = (state.cells[cellIndex].MliqiniL) / (leftLiquidDensity * leftFlowArea);

        double rightFlowArea = state.cells[cellIndex].dutoR.area;
        double rightLengthRatio = state.cells[cellIndex].dxR / (state.cells[cellIndex].dxR + state.cells[cellIndex].dx);
        double rightMeanPressure = state.cells[cellIndex].presauxR;
        double rightMeanTemperature = state.cells[cellIndex].temp * rightLengthRatio + state.cells[cellIndex].tempR * (1. - rightLengthRatio);
        double betIR = state.cells[cellIndex].betPigD;
        if (state.cells[cellIndex].QLR < 0.) { // testeBeta
            if (cellIndex > state.lastCell - 2)
                betIR = state.cells[cellIndex].betR;
            else
                betIR = state.cells[cellIndex + 1].betPigE;
        }
        double rightGasDensity = state.cells[cellIndex].rgRi;
        double rightLiquidDensity = (1 - betIR) * state.cells[cellIndex].rpRi + betIR * state.cells[cellIndex].rcRi;
        double rightSuperficialGasVelocity = (state.cells[cellIndex].MR - state.cells[cellIndex].MliqiniR) / (rightGasDensity * rightFlowArea);
        double rightSuperficialLiquidVelocity = (state.cells[cellIndex].MliqiniR) / (rightLiquidDensity * rightFlowArea);

        selectAndApplyInteriorFlowRegime(
            state, cellIndex, bif, superficialGasVelocity, superficialLiquidVelocity, leftSuperficialGasVelocity, leftSuperficialLiquidVelocity, rightSuperficialGasVelocity, rightSuperficialLiquidVelocity,
            rightGasDensity, rightLiquidDensity, gasDensity, liquidDensity, flowArea);
    } else {
        state.cells[cellIndex].c0 = 1.;
        state.cells[cellIndex].ud = 0.;
        state.cells[cellIndex].term1 = 0.;
        state.cells[cellIndex].term2 = 0.;
        state.cells[cellIndex].term1L = state.cells[cellIndex - 1].term1;
        state.cells[cellIndex].term2L = state.cells[cellIndex - 1].term2;
        state.cells[cellIndex - 1].term1R = state.cells[cellIndex].term1;
        state.cells[cellIndex - 1].term2R = state.cells[cellIndex].term2;
    }
}

void updateOutletBoundaryFlowPartition(
    const ThermalState &state, int cellIndex, Vcr<int> &bif) {
    state.cells[state.lastCell - 1].alfR = state.cells[state.lastCell].alf;
    state.cells[state.lastCell].alfR = state.cells[state.lastCell].alf;
    state.cells[state.lastCell - 1].betR = state.cells[state.lastCell].bet;
    state.cells[state.lastCell].betR = state.cells[state.lastCell].bet;

    double lengthRatio = state.cells[cellIndex].dxL / (state.cells[cellIndex].dx + state.cells[cellIndex].dxL);
    double meanPressure = state.cells[cellIndex].presaux;
    double meanTemperature = state.cells[cellIndex].temp * lengthRatio + state.cells[cellIndex - 1].temp * (1. - lengthRatio);
    double betI = state.cells[cellIndex].betL;
    if (state.cells[cellIndex].QL < 0.)
        betI = state.cells[cellIndex].bet; // testeBeta
    // betI = state.cells[cellIndex].bet;            //duvidabeta
    double gasDensity = state.cells[cellIndex].flui.MasEspGas(meanPressure, meanTemperature);
    double liquidDensity = (1 - betI) * state.cells[cellIndex].flui.MasEspLiq(meanPressure, meanTemperature) + betI * state.cells[cellIndex].fluicol.MasEspFlu(meanPressure, meanTemperature);
    double flowArea = state.cells[cellIndex].duto.area;
    if (state.cells[cellIndex].MC >= 0)
        flowArea = state.cells[cellIndex].dutoL.area;
    double superficialGasVelocity = state.cells[cellIndex].QG / (flowArea);
    double superficialLiquidVelocity = state.cells[cellIndex].QL / (flowArea);

    double leftFlowArea = state.cells[cellIndex].dutoL.area;
    double leftLengthRatio = state.cells[cellIndex - 1].dxL / (state.cells[cellIndex - 1].dx + state.cells[cellIndex - 1].dxL);
    double leftMeanPressure = state.cells[cellIndex - 1].presaux;
    double leftMeanTemperature = state.cells[cellIndex - 1].temp * leftLengthRatio + state.cells[cellIndex - 1].tempL * (1. - leftLengthRatio);
    double betIL = state.cells[cellIndex - 1].betL;
    if (state.cells[cellIndex - 1].QL < 0.)
        betIL = state.cells[cellIndex - 1].bet; // testeBeta
    // betIL = state.cells[cellIndex - 1].bet;            //duvidabeta
    double leftGasDensity = state.cells[cellIndex].flui.MasEspGas(leftMeanPressure, leftMeanTemperature);
    double leftLiquidDensity = (1 - betIL) * state.cells[cellIndex].flui.MasEspLiq(leftMeanPressure, leftMeanTemperature) + betIL * state.cells[cellIndex].fluicol.MasEspFlu(leftMeanPressure, leftMeanTemperature);
    double leftSuperficialLiquidVelocity = (state.cells[cellIndex].MliqiniL) / (leftLiquidDensity * leftFlowArea);

    double rightFlowArea = state.cells[cellIndex].dutoR.area;
    double rightLengthRatio = state.cells[cellIndex].dxR / (state.cells[cellIndex].dxR + state.cells[cellIndex].dx);
    double rightMeanPressure = state.cells[cellIndex].pres;
    double rightMeanTemperature = state.cells[cellIndex].temp * rightLengthRatio + state.cells[cellIndex].tempR * (1. - rightLengthRatio);
    double betIR = state.cells[cellIndex].bet;
    double rightGasDensity = state.cells[cellIndex].flui.MasEspGas(rightMeanPressure, rightMeanTemperature);
    double rightLiquidDensity = (1 - betIR) * state.cells[cellIndex].flui.MasEspLiq(rightMeanPressure, rightMeanTemperature) + betIR * state.cells[cellIndex].fluicol.MasEspFlu(rightMeanPressure, rightMeanTemperature);
    double rightSuperficialLiquidVelocity = (state.cells[cellIndex].MliqiniR) / (rightLiquidDensity * rightFlowArea);

    bif[cellIndex] = 1;

    if (state.cells[cellIndex].alfL <= (*state.globals).localtiny && state.cells[cellIndex].alf <= (*state.globals).localtiny && state.cells[cellIndex].fontemassGL <= 0 && state.cells[cellIndex].fontemassGR <= 0) {
        state.cells[cellIndex].term1 = 1.;
        state.cells[cellIndex].term2 = 0.;
        bif[cellIndex] = 0;
        state.cells[cellIndex].c0 = 1. - (*state.globals).localtiny;
        state.cells[cellIndex].ud = 0.;
        state.cells[cellIndex].arranjo = 0;
    } else if (state.cells[cellIndex].alfL >= (1. - (*state.globals).localtiny) && state.cells[cellIndex].alf >= (1. - (*state.globals).localtiny) && state.cells[cellIndex].fontemassLL <= 0 && state.cells[cellIndex].fontemassLR <= 0) {
        state.cells[cellIndex].term1 = 0.;
        state.cells[cellIndex].term2 = 0.;
        bif[cellIndex] = 0;
        state.cells[cellIndex].c0 = 1. - (*state.globals).localtiny;
        state.cells[cellIndex].ud = 0.;
        state.cells[cellIndex].arranjo = 0;
    } else if (superficialGasVelocity > 0 && state.cells[cellIndex].alfL <= (*state.globals).localtiny && (state.cells[cellIndex].fontemassGL <= 0. && state.cells[cellIndex].fontemassGR <= 0.)) {
        state.cells[cellIndex].term1 = 1.;
        state.cells[cellIndex].term2 = 0.;
        bif[cellIndex] = 0;
        state.cells[cellIndex].c0 = 1. - (*state.globals).localtiny;
        state.cells[cellIndex].ud = 0.;
        state.cells[cellIndex].arranjo = 0;
        if (fabs(superficialGasVelocity) <= 1e-15 && state.cells[cellIndex].alf > (1. - (*state.globals).localtiny) && superficialLiquidVelocity < 0 && leftSuperficialLiquidVelocity < 0) {
            state.cells[cellIndex].term1 = 0.;
            state.cells[cellIndex].term2 = 0.;
            bif[cellIndex] = 0;
        }
    } else if (superficialGasVelocity < 0 && state.cells[cellIndex].alf <= (*state.globals).localtiny && (state.cells[cellIndex].fontemassGL <= 0. && state.cells[cellIndex].fontemassGR <= 0.)) {
        state.cells[cellIndex].term1 = 1.;
        state.cells[cellIndex].term2 = 0.;
        bif[cellIndex] = 0;
        state.cells[cellIndex].c0 = 1. - (*state.globals).localtiny;
        state.cells[cellIndex].ud = 0.;
        state.cells[cellIndex].arranjo = 0;
        if (fabs(superficialGasVelocity) <= 1e-15 && state.cells[cellIndex].alfL > (1. - (*state.globals).localtiny) && superficialLiquidVelocity > 0 && rightSuperficialLiquidVelocity > 0) {
            state.cells[cellIndex].term1 = 0.;
            state.cells[cellIndex].term2 = 0.;
            bif[cellIndex] = 0;
        }
    } else if (superficialLiquidVelocity >= 0 && state.cells[cellIndex - 1].alf >= 1. - (*state.globals).localtiny && ((state.cells[cellIndex].fontemassLL + state.cells[cellIndex].fontemassCL) <= (*state.globals).localtiny * 1e-5 && (state.cells[cellIndex].fontemassLR + state.cells[cellIndex].fontemassCR) <= (*state.globals).localtiny * 1e-5)) {
        state.cells[cellIndex].term1 = 0.;
        state.cells[cellIndex].term2 = 0.;
        bif[cellIndex] = 0;
        state.cells[cellIndex].c0 = 1;
        state.cells[cellIndex].ud = 0;
        state.cells[cellIndex].arranjo = 0;
        if (fabs(superficialLiquidVelocity) <= 1e-15 && state.cells[cellIndex].alf < (*state.globals).localtiny && rightSuperficialLiquidVelocity < 0) {
            state.cells[cellIndex].term1 = 1.;
            state.cells[cellIndex].term2 = 0.;
            bif[cellIndex] = 0;
        } else if (fabs(rightSuperficialLiquidVelocity) < (*state.globals).localtiny * 1e-5) {
            if (fabs(superficialLiquidVelocity) < (*state.globals).localtiny * 1e-5 && state.cells[cellIndex].alf < (*state.globals).localtiny && state.cells[cellIndex].fontemassGR >= (*state.globals).localtiny * 1e-5) {
                state.cells[cellIndex].term1 = 1.;
                state.cells[cellIndex].term2 = 0.;
                bif[cellIndex] = 0;
            }
        } else if (fabs(superficialLiquidVelocity) < 1e-15 && rightSuperficialLiquidVelocity < 0 && state.cells[cellIndex].alf <= (1 - 1 * (*state.globals).localtiny))
            bif[cellIndex] = 1;
    }

    else if (superficialLiquidVelocity <= 0 && state.cells[cellIndex].alf >= 1. - (*state.globals).localtiny && ((state.cells[cellIndex].fontemassLL + state.cells[cellIndex].fontemassCL) <= (*state.globals).localtiny * 1e-5 && (state.cells[cellIndex].fontemassLR + state.cells[cellIndex].fontemassCR) <= (*state.globals).localtiny * 1e-5)) {
        state.cells[cellIndex].term1 = 0.;
        state.cells[cellIndex].term2 = 0.;
        bif[cellIndex] = 0;
        state.cells[cellIndex].c0 = 1;
        state.cells[cellIndex].ud = 0;
        state.cells[cellIndex].arranjo = 0;
        if (fabs(superficialLiquidVelocity) <= (*state.globals).localtiny * 1e-5 && state.cells[cellIndex - 1].alf < (*state.globals).localtiny && leftSuperficialLiquidVelocity > 0) {
            state.cells[cellIndex].term1 = 1.;
            state.cells[cellIndex].term2 = 0.;
            bif[cellIndex] = 0;
        }
        if (fabs(leftSuperficialLiquidVelocity) < (*state.globals).localtiny * 1e-5) {
            if (fabs(superficialLiquidVelocity) < (*state.globals).localtiny * 1e-5 && state.cells[cellIndex - 1].alf < (*state.globals).localtiny && state.cells[cellIndex - 1].fontemassGR >= (*state.globals).localtiny * 1e-5) {
                state.cells[cellIndex].term1 = 1.;
                state.cells[cellIndex].term2 = 0.;
                bif[cellIndex] = 0;
            }
        } else if (fabs(superficialLiquidVelocity) < (*state.globals).localtiny * 1e-5 && state.cells[cellIndex].duto.teta < 0.95 * M_PI / 2. && leftSuperficialLiquidVelocity > 0 && ((state.cells[cellIndex - 1].alfPigD <= 1.0 * state.cells[cellIndex - 2].alfPigD && state.cells[cellIndex - 2].alfPigD < 0.99) || state.cells[cellIndex - 1].alfPigD < 0.7))
            bif[cellIndex] = 1;
        else if (fabs(superficialLiquidVelocity) < (*state.globals).localtiny * 1e-5 && state.cells[cellIndex].duto.teta >= 0.95 * M_PI / 2. && leftSuperficialLiquidVelocity > 0 && ((state.cells[cellIndex - 1].alfPigD <= 1.0 * state.cells[cellIndex - 2].alfPigD && state.cells[cellIndex - 2].alfPigD < 0.99) || state.cells[cellIndex - 1].alfPigD < 0.7))
            bif[cellIndex] = 1;
        else {
            if (fabs(superficialLiquidVelocity) < (*state.globals).localtiny * 1e-5 && state.cells[cellIndex].duto.teta < 0.95 * M_PI / 2. && leftSuperficialLiquidVelocity > 0 && ((state.cells[cellIndex - 1].alfPigD <= 1.0 * state.cells[cellIndex - 1].alfL) || state.cells[cellIndex - 1].alfPigD < 0.7))
                bif[cellIndex] = 1;
            else if (fabs(superficialLiquidVelocity) < (*state.globals).localtiny * 1e-5 && state.cells[cellIndex].duto.teta >= 0.95 * M_PI / 2. && leftSuperficialLiquidVelocity > 0 && ((state.cells[cellIndex - 1].alfPigD <= 1.0 * state.cells[cellIndex - 1].alfL) || state.cells[cellIndex - 1].alfPigD < 0.7))
                bif[cellIndex] = 1;
        }
    }

    if (bif[cellIndex] == 1) {
        double meanVoidFraction;
        meanVoidFraction = state.cells[cellIndex].alfL;
        double c0 = 1.2;
        double meanDiameter = state.cells[cellIndex].duto.a;
        if (state.cells[cellIndex].MC >= 0)
            meanDiameter = state.cells[cellIndex].dutoL.a;
        double inclinationSign = 1.;
        if (state.cells[cellIndex].duto.teta < 0.)
            inclinationSign = 1.;
        double ud = inclinationSign * 0.32 * sqrt(9.82 * meanDiameter);
        if (fabs(rightGasDensity) / rightLiquidDensity > 0.9) {
            c0 = 1.;
            ud = 0.;
        }
        state.closureUpdater.instantaneous(cellIndex, c0, ud);
        if (state.input.escorregamentoCelulaContorno == 0) {
            c0 = 1.;
            ud = 0.;
        }
        state.cells[cellIndex].c0 = c0;
        state.cells[cellIndex].ud = ud;
        double numerator = (1. - meanVoidFraction * c0);
        double denominator = 1. + meanVoidFraction * (gasDensity / liquidDensity) * c0 - meanVoidFraction * c0;
        state.cells[cellIndex].term1 = numerator / denominator;
        state.cells[cellIndex].term2 = (-flowArea * meanVoidFraction * gasDensity * ud) / denominator;

        // teste!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        // teste!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    }

    state.cells[cellIndex].term1L = state.cells[cellIndex - 1].term1;
    state.cells[cellIndex].term2L = state.cells[cellIndex - 1].term2;
    state.cells[cellIndex - 1].term1R = state.cells[cellIndex].term1;
    state.cells[cellIndex - 1].term2R = state.cells[cellIndex].term2;
}

void finalizeFlowPartitionTerms(
    const ThermalState &state, const Vcr<int> &bif, const Vcr<int> &valv) {
    int isShutIn = 0;
    for (int cellIndex = 1; cellIndex < state.lastCell; cellIndex++) {
        if (state.cells[cellIndex].acsr.tipo == 5 && state.cells[cellIndex].acsr.chk.AreaGarg <= 1e-15 * state.cells[cellIndex].acsr.chk.AreaTub)
            isShutIn = 1;
        else if (state.surfaceChoke.AreaGarg <= 1.e-15 * state.surfaceChoke.AreaTub)
            isShutIn = 1;
    }

    Vcr<double> c0V(state.lastCell, 0.);
    Vcr<double> udV(state.lastCell, 0.);
    for (int cellIndex = 1; cellIndex < state.lastCell; cellIndex++) {
        c0V[cellIndex] = state.cells[cellIndex].c0;
        udV[cellIndex] = state.cells[cellIndex].ud;
        if (bif[cellIndex] == 1 && valv[cellIndex] == 1 && cellIndex > 2) {
            int neighborIndex = cellIndex - 1;
            int secondNeighborIndex = cellIndex - 2;
            if (state.cells[cellIndex].QG < 0) {
                neighborIndex = cellIndex + 1;
                secondNeighborIndex = cellIndex;
            }
            if ((state.cells[neighborIndex].acsr.tipo == 0 && (state.cells[secondNeighborIndex].acsr.tipo != 5 || state.cells[secondNeighborIndex].acsr.chk.AreaGarg > (1e-3))) &&
                (state.cells[cellIndex].arranjo != state.cells[neighborIndex].arranjo && bif[neighborIndex] != 0)) {
                c0V[cellIndex] = (state.cells[cellIndex].dx * state.cells[cellIndex].c0 + state.cells[neighborIndex].dx * state.cells[neighborIndex].c0) / (state.cells[cellIndex].dx + state.cells[neighborIndex].dx);
                if (state.cells[cellIndex].duto.teta * state.cells[neighborIndex].duto.teta >= 0)
                    udV[cellIndex] = (state.cells[cellIndex].dx * state.cells[cellIndex].ud + state.cells[neighborIndex].dx * state.cells[neighborIndex].ud) / (state.cells[cellIndex].dx + state.cells[neighborIndex].dx);
            } else if (cellIndex > 2) {
                double inclinationAngle;
                double leftInclinationAngle;
                double diameter;
                double leftDiameter;
                if ((state.cells[cellIndex - 1].acsr.tipo == 0 && (state.cells[cellIndex - 2].acsr.tipo != 5 || state.cells[cellIndex - 2].acsr.chk.AreaGarg > (1e-3))) &&
                    state.cells[cellIndex].QG >= 0) {
                    inclinationAngle = state.cells[cellIndex].duto.teta;
                    leftInclinationAngle = state.cells[cellIndex - 1].duto.teta;
                    diameter = state.cells[cellIndex].duto.dia;
                    leftDiameter = state.cells[cellIndex - 1].duto.dia;
                    if (((inclinationAngle != leftInclinationAngle) && bif[neighborIndex] != 0) ||
                        (state.cells[cellIndex].arranjo != state.cells[neighborIndex].arranjo && bif[neighborIndex] != 0)) {
                        c0V[cellIndex] = (state.cells[cellIndex].dx * state.cells[cellIndex].c0 + state.cells[cellIndex - 1].dx * state.cells[cellIndex - 1].c0) / (state.cells[cellIndex].dx + state.cells[cellIndex - 1].dx);
                        if (inclinationAngle * leftInclinationAngle >= 0)
                            udV[cellIndex] = (state.cells[cellIndex].dx * state.cells[cellIndex].ud + state.cells[cellIndex - 1].dx * state.cells[cellIndex - 1].ud) / (state.cells[cellIndex].dx + state.cells[cellIndex - 1].dx);
                    }
                } else if ((state.cells[cellIndex + 1].acsr.tipo == 0 && (state.cells[cellIndex].acsr.tipo != 5 || state.cells[cellIndex].acsr.chk.AreaGarg > (1e-3))) &&
                           state.cells[cellIndex].QG < 0) {
                    inclinationAngle = state.cells[cellIndex].duto.teta;
                    leftInclinationAngle = state.cells[cellIndex + 1].duto.teta;
                    diameter = state.cells[cellIndex].duto.dia;
                    leftDiameter = state.cells[cellIndex + 1].duto.dia;
                    if (((inclinationAngle != leftInclinationAngle) && bif[neighborIndex] != 0) ||
                        (state.cells[cellIndex].arranjo != state.cells[neighborIndex].arranjo && bif[neighborIndex] != 0)) {
                        c0V[cellIndex] = (state.cells[cellIndex].dx * state.cells[cellIndex].c0 + state.cells[cellIndex + 1].dx * state.cells[cellIndex + 1].c0) / (state.cells[cellIndex].dx + state.cells[cellIndex + 1].dx);
                        if (inclinationAngle * leftInclinationAngle >= 0)
                            udV[cellIndex] = (state.cells[cellIndex].dx * state.cells[cellIndex].ud + state.cells[cellIndex + 1].dx * state.cells[cellIndex + 1].ud) / (state.cells[cellIndex].dx + state.cells[cellIndex + 1].dx);
                    }
                }
            }
        }
    }

    for (int cellIndex = 1; cellIndex < state.lastCell; cellIndex++) {
        if (bif[cellIndex] == 1 && valv[cellIndex] == 1) {
            double betI = state.cells[cellIndex - 1].betPigD;
            double liquidDensity;

            if (state.cells[cellIndex].QL < 0.) {
                betI = state.cells[cellIndex].betPigE;
                liquidDensity = (1 - betI) * state.cells[cellIndex].rpCi + betI * state.cells[cellIndex].rcCi;
            } else {
                betI = state.cells[cellIndex - 1].betPigD;
                liquidDensity = (1 - betI) * state.cells[cellIndex - 1].rpCi + betI * state.cells[cellIndex - 1].rcCi;
            }

            double gasDensity;
            double flowArea;
            double noSlipLiquidHoldup;
            if (state.cells[cellIndex].QG >= 0) {
                flowArea = state.cells[cellIndex].dutoL.area;
                noSlipLiquidHoldup = 1. - state.cells[cellIndex].alfL;
                gasDensity = state.cells[cellIndex - 1].rgCi;
            } else {
                gasDensity = state.cells[cellIndex].rgCi;
                flowArea = state.cells[cellIndex].duto.area;
                noSlipLiquidHoldup = 1. - state.cells[cellIndex].alf;
            }
            double superficialGasVelocity = state.cells[cellIndex].QG / (flowArea);
            double superficialLiquidVelocity = state.cells[cellIndex].QL / (flowArea);

            double meanVoidFraction;
            meanVoidFraction = state.cells[cellIndex - 1].alfPigD;
            if (superficialGasVelocity < 0)
                meanVoidFraction = state.cells[cellIndex].alfPigE;
            double numerator = (1. - meanVoidFraction * state.cells[cellIndex].c0);
            double denominator = 1 + state.cells[cellIndex].c0 * meanVoidFraction * ((gasDensity / liquidDensity) - 1.);
            state.cells[cellIndex].term1 = numerator / denominator;
            state.cells[cellIndex].term2 = (-flowArea * meanVoidFraction * gasDensity * state.cells[cellIndex].ud) / denominator;
            double previousLiquidDriftFlux = (superficialGasVelocity - meanVoidFraction * state.cells[cellIndex].ud) / (meanVoidFraction * state.cells[cellIndex].c0) - superficialGasVelocity;
            double liquidDriftFlux = (superficialGasVelocity + superficialLiquidVelocity) * (1. - state.cells[cellIndex].c0 * meanVoidFraction) - meanVoidFraction * state.cells[cellIndex].ud;
            double candidateLiquidMass = state.cells[cellIndex].term1 * state.cells[cellIndex].MC + state.cells[cellIndex].term2;
            double candidateGasMass = (1 - state.cells[cellIndex].term1) * state.cells[cellIndex].MC - state.cells[cellIndex].term2;
            if ((liquidDriftFlux > 0. || previousLiquidDriftFlux > 0.) && state.cells[cellIndex - 1].alfPigD > 1 - 1e-15) {
                state.cells[cellIndex].term1 = 0.;
                state.cells[cellIndex].term2 = 0.;
            }
            if ((liquidDriftFlux < 0. || previousLiquidDriftFlux < 0.) && state.cells[cellIndex].alfPigE > 1 - 1e-15) {
                state.cells[cellIndex].term1 = 0.;
                state.cells[cellIndex].term2 = 0.;
            }
            if (state.cells[cellIndex - 1].acsr.tipo == 5 && state.cells[cellIndex - 1].acsr.chk.AreaGarg <= (1e-3 + state.input.master1.razareaativ) * state.cells[cellIndex - 1].duto.area) {
                state.cells[cellIndex].term1 = 0.;
                state.cells[cellIndex].term2 = 0.;
            }
            if ((liquidDriftFlux < 0. || previousLiquidDriftFlux < 0.) && ((fabs(state.cells[cellIndex - 1].QG / (flowArea)) + fabs(state.cells[cellIndex - 1].QL / (flowArea))) < 0.1) &&
                (isShutIn == 1 && state.input.modoSegrega == 1) && (state.cells[cellIndex].duto.teta > 0 && state.cells[cellIndex - 1].duto.teta <= 0) &&
                (candidateLiquidMass > 0 && state.cells[cellIndex].term2 > 0) &&
                ((1. - state.cells[cellIndex].alfL) /**fabs(sin(state.cells[cellIndex-1].duto.teta))*/ < (1. - state.cells[cellIndex].alf) /**fabs(sin(state.cells[cellIndex].duto.teta))*/)) {
                state.cells[cellIndex].term2 = 0.;
            }
            if (state.cells[cellIndex].duto.teta > 0 && candidateGasMass < 0 && ((isShutIn == 1 && state.input.modoSegrega == 1)) && state.cells[cellIndex].alfPigE > 1 - 1e-15) {
                state.cells[cellIndex].term1 = 0.;
                state.cells[cellIndex].term2 = 0.;
            }
        }
    }
    for (int cellIndex = 1; cellIndex <= state.lastCell; cellIndex++) {
        state.cells[cellIndex].term1L = state.cells[cellIndex - 1].term1;
        state.cells[cellIndex].term2L = state.cells[cellIndex - 1].term2;
        state.cells[cellIndex - 1].term1R = state.cells[cellIndex].term1;
        state.cells[cellIndex - 1].term2R = state.cells[cellIndex].term2;
    }
}

void selectAndApplyInletBoundaryFlowRegime(
    const ThermalState &state, int cellIndex, Vcr<int> &bif, double distributionCoefficient,
    double driftVelocity, double superficialGasVelocity, double superficialLiquidVelocity, double rightSuperficialLiquidVelocity, double rightGasDensity,
    double rightLiquidDensity, double gasDensity, double liquidDensity, double flowArea) {
    bif[cellIndex] = 1;

    if ((*state.globals).lixo5 > 29900) {
        int debugStop;
        debugStop = 0;
    }

    if (state.inletVoidFraction < (*state.globals).localtiny && state.cells[cellIndex].alfPigE <= (*state.globals).localtiny && state.cells[cellIndex].fontemassGL <= (*state.globals).localtiny * 1e-5 && state.cells[cellIndex].fontemassGR <= (*state.globals).localtiny * 1e-5) {
        state.cells[cellIndex].term1 = 1.;
        state.cells[cellIndex].term2 = 0.;
        bif[cellIndex] = 0;
        state.cells[cellIndex].c0 = 1 + 0 * distributionCoefficient;
        state.cells[cellIndex].ud = 0 * driftVelocity;
        state.cells[cellIndex].arranjo = 0;
    } else if (state.inletVoidFraction >= (1. - (*state.globals).localtiny) && state.cells[cellIndex].alfPigE >= (1. - (*state.globals).localtiny) && ((state.cells[cellIndex].fontemassLL + state.cells[cellIndex].fontemassCL) <= (*state.globals).localtiny * 1e-5 && (state.cells[cellIndex].fontemassLR + state.cells[cellIndex].fontemassCR) <= (*state.globals).localtiny * 1e-5)) {
        state.cells[cellIndex].term1 = 0.;
        state.cells[cellIndex].term2 = 0.;
        bif[cellIndex] = 0;
        state.cells[cellIndex].c0 = 1 + 0 * distributionCoefficient;
        state.cells[cellIndex].ud = 0 * driftVelocity;
        state.cells[cellIndex].arranjo = 0;
    } else if (superficialGasVelocity >= 0 && state.inletVoidFraction <= (*state.globals).localtiny && (state.cells[cellIndex].fontemassGL <= (*state.globals).localtiny * 1e-5 && state.cells[cellIndex].fontemassGR <= (*state.globals).localtiny * 1e-5)) {
        state.cells[cellIndex].term1 = 1.;
        state.cells[cellIndex].term2 = 0.;
        bif[cellIndex] = 0;
        state.cells[cellIndex].c0 = 1 + 0 * distributionCoefficient;
        state.cells[cellIndex].ud = 0 * driftVelocity;
        state.cells[cellIndex].arranjo = 0;
        if (fabs(superficialGasVelocity) <= 1e-15 && state.cells[cellIndex].alfPigE > (1. - (*state.globals).localtiny) && superficialLiquidVelocity < 0 && state.cells[cellIndex].duto.teta > 0) {
            state.cells[cellIndex].term1 = 0.;
            state.cells[cellIndex].term2 = 0.;
            bif[cellIndex] = 0;
        }
        if (fabs(superficialGasVelocity) <= 1e-15 && state.cells[cellIndex].alfPigE > (1. - 1 * (*state.globals).localtiny) && state.cells[cellIndex].duto.teta < 0 && superficialLiquidVelocity > 0) {
            state.cells[cellIndex].term1 = 0.;
            state.cells[cellIndex].term2 = 0.;
            bif[cellIndex] = 1;
        }
        if (fabs(superficialGasVelocity) <= 1e-15 && state.cells[cellIndex].alf >= state.inletVoidFraction && superficialLiquidVelocity < 0) {
            bif[cellIndex] = 1;
        }
    } else if (superficialGasVelocity <= 0 && state.cells[cellIndex].alfPigE <= (*state.globals).localtiny && (state.cells[cellIndex].fontemassGL <= (*state.globals).localtiny * 1e-5 && state.cells[cellIndex].fontemassGR <= (*state.globals).localtiny * 1e-5)) {
        state.cells[cellIndex].term1 = 1.;
        state.cells[cellIndex].term2 = 0.;
        bif[cellIndex] = 0;
        state.cells[cellIndex].c0 = 1 + 0 * distributionCoefficient;
        state.cells[cellIndex].ud = 0 * driftVelocity;
        state.cells[cellIndex].arranjo = 0;
        if (fabs(superficialGasVelocity) <= 1e-15 && state.inletVoidFraction > (1. - (*state.globals).localtiny) && superficialLiquidVelocity > 0) {
            state.cells[cellIndex].term1 = 0.;
            state.cells[cellIndex].term2 = 0.;
            bif[cellIndex] = 0;
        }
        if (fabs(superficialGasVelocity) <= 1e-15 && state.inletVoidFraction > (*state.globals).localtiny && superficialLiquidVelocity > 0) {
            bif[cellIndex] = 1;
        }
    } else if (superficialLiquidVelocity >= 0 && state.inletVoidFraction >= 1. - (*state.globals).localtiny && ((state.cells[cellIndex].fontemassLL + state.cells[cellIndex].fontemassCL) <= (*state.globals).localtiny * 1e-5 && (state.cells[cellIndex].fontemassLR + state.cells[cellIndex].fontemassCR) <= (*state.globals).localtiny * 1e-5)) {
        state.cells[cellIndex].term1 = 0.;
        state.cells[cellIndex].term2 = 0.;
        bif[cellIndex] = 0;
        state.cells[cellIndex].c0 = 1 + 0 * distributionCoefficient;
        state.cells[cellIndex].ud = 0 * driftVelocity;
        state.cells[cellIndex].arranjo = 0;
        if (fabs(superficialLiquidVelocity) <= 1e-15 && state.cells[cellIndex].alfPigE < (*state.globals).localtiny && rightSuperficialLiquidVelocity < 0) {
            state.cells[cellIndex].term1 = 1.;
            state.cells[cellIndex].term2 = 0.;
            bif[cellIndex] = 0;
        } else if (fabs(rightSuperficialLiquidVelocity) < (*state.globals).localtiny * 1e-5) {
            if (fabs(superficialLiquidVelocity) < (*state.globals).localtiny * 1e-5 && state.cells[cellIndex].alfPigE < (*state.globals).localtiny && state.cells[cellIndex].fontemassGR >= (*state.globals).localtiny * 1e-5) {
                state.cells[cellIndex].term1 = 1.;
                state.cells[cellIndex].term2 = 0.;
                bif[cellIndex] = 0;
            }
        } else if (fabs(superficialLiquidVelocity) < 1e-15 && rightSuperficialLiquidVelocity < 0 && ((state.cells[cellIndex].alfPigE <= (1 - 1 * (*state.globals).localtiny + .0 * state.cells[cellIndex].alfPigER) && state.cells[cellIndex].alfPigER < 1 - 1 * (*state.globals).localtiny) || state.cells[cellIndex].alfPigE <= 0.7))
            bif[cellIndex] = 1;

    } else if (superficialLiquidVelocity <= 0 && state.cells[cellIndex].alfPigE >= 1. - (*state.globals).localtiny && ((state.cells[cellIndex].fontemassLL + state.cells[cellIndex].fontemassCL) <= (*state.globals).localtiny * 1e-5 && (state.cells[cellIndex].fontemassLR + state.cells[cellIndex].fontemassCR) <= (*state.globals).localtiny * 1e-5)) {
        state.cells[cellIndex].term1 = 0.;
        state.cells[cellIndex].term2 = 0.;
        bif[cellIndex] = 0;
        state.cells[cellIndex].c0 = 1 + 0 * distributionCoefficient;
        state.cells[cellIndex].ud = 0 * driftVelocity;
        state.cells[cellIndex].arranjo = 0;
        if (fabs(superficialLiquidVelocity) <= (*state.globals).localtiny * 1e-5 && state.inletVoidFraction < (*state.globals).localtiny && rightSuperficialLiquidVelocity < 0) {
            state.cells[cellIndex].term1 = 1.;
            state.cells[cellIndex].term2 = 0.;
            bif[cellIndex] = 0;
        }

        if (fabs(rightSuperficialLiquidVelocity) < (*state.globals).localtiny * 1e-5) {
            if (fabs(superficialLiquidVelocity) < (*state.globals).localtiny * 1e-5 && state.inletVoidFraction < (*state.globals).localtiny) {
                state.cells[cellIndex].term1 = 1.;
                state.cells[cellIndex].term2 = 0.;
                bif[cellIndex] = 0;
            }
        } else {
            if (fabs(superficialLiquidVelocity) < (*state.globals).localtiny * 1e-5 && state.cells[cellIndex].duto.teta < 0.95 * M_PI / 2. && superficialGasVelocity < 0 && (state.inletVoidFraction < 0.7))
                bif[cellIndex] = 1;
            else if (fabs(superficialLiquidVelocity) < (*state.globals).localtiny * 1e-5 && state.cells[cellIndex].duto.teta >= 0.95 * M_PI / 2. && superficialGasVelocity < 0 && (state.inletVoidFraction < 0.7))
                bif[cellIndex] = 1;
        }
    }
    if (superficialLiquidVelocity > 0 && fabs(superficialGasVelocity) <= (*state.globals).localtiny * 1e-5 && state.inletVoidFraction > (*state.globals).localtiny && state.inletVoidFraction < 1 - (*state.globals).localtiny)
        bif[cellIndex] = 1;
    if (superficialGasVelocity > 0 && fabs(superficialLiquidVelocity) <= (*state.globals).localtiny * 1e-5 && state.inletVoidFraction > (*state.globals).localtiny && state.inletVoidFraction < 1 - (*state.globals).localtiny)
        bif[cellIndex] = 1;
    if (superficialLiquidVelocity >= 0 && state.inletVoidFraction > 1 - (*state.globals).localtiny) {
        state.cells[cellIndex].term1 = 0.;
        state.cells[cellIndex].term2 = 0.;
        bif[cellIndex] = 0;
    }
    if (superficialGasVelocity >= 0 && state.inletVoidFraction < (*state.globals).localtiny) {
        state.cells[cellIndex].term1 = 1.;
        state.cells[cellIndex].term2 = 0.;
        bif[cellIndex] = 0;
    }

    if (bif[cellIndex] == 1) {

        double c0;
        double ud;
        double meanVoidFraction;

        meanVoidFraction = state.inletVoidFraction;
        if (superficialGasVelocity < 0)
            meanVoidFraction = state.cells[cellIndex].alfPigE;
        c0 = 1.2;
        double meanDiameter = state.cells[cellIndex].duto.a;
        if (state.cells[cellIndex].MC >= 0)
            meanDiameter = state.cells[cellIndex].dutoL.a;
        double inclinationSign = 1.;
        if (state.cells[cellIndex].duto.teta < 0.)
            inclinationSign = -1.;
        ud = inclinationSign * 0.32 * sqrt(9.82 * meanDiameter);
        if (fabs(rightGasDensity) / rightLiquidDensity > 0.9) {
            c0 = 1.;
            ud = 0.;
        }
        state.closureUpdater.initialization(cellIndex, c0, ud);
        state.cells[cellIndex].c0 = c0;
        state.cells[cellIndex].ud = ud;
        double numerator = (1. - meanVoidFraction * c0);
        double denominator = 1 + c0 * meanVoidFraction * ((gasDensity / liquidDensity) - 1.);
        state.cells[cellIndex].term1 = numerator / denominator;
        state.cells[cellIndex].term2 = (-flowArea * meanVoidFraction * gasDensity * ud) / denominator;
    }
}

void selectAndApplyBufferedOutletFlowRegime(
    const ThermalState &state, int cellIndex, double distributionCoefficient, double driftVelocity,
    double superficialGasVelocity, double superficialLiquidVelocity, double leftSuperficialGasVelocity, double leftSuperficialLiquidVelocity, double rightSuperficialGasVelocity,
    double rightSuperficialLiquidVelocity, double rightGasDensity, double rightLiquidDensity, double gasDensity, double liquidDensity,
    double flowArea) {
    int bif = 1;

    if ((*state.globals).lixo5 > 29900) {
        int debugStop;
        debugStop = 0;
    }

    if (state.cells[cellIndex - 1].alfPigD <= (*state.globals).localtiny && state.cells[cellIndex].alfPigE <= (*state.globals).localtiny && state.cells[cellIndex].fontemassGL <= (*state.globals).localtiny * 1e-5 && state.cells[cellIndex].fontemassGR <= (*state.globals).localtiny * 1e-5) {
        state.cells[cellIndex].term1 = 1.;
        state.cells[cellIndex].term2 = 0.;
        bif = 0;
        state.cells[cellIndex].c0 = 1 + 0 * distributionCoefficient;
        state.cells[cellIndex].ud = 0 * driftVelocity;
        state.cells[cellIndex].arranjo = 0;
    } else if (state.cells[cellIndex - 1].alfPigD >= (1. - (*state.globals).localtiny) && state.cells[cellIndex].alfPigE >= (1. - (*state.globals).localtiny) && ((state.cells[cellIndex].fontemassLL + state.cells[cellIndex].fontemassCL) <= (*state.globals).localtiny * 1e-5 && (state.cells[cellIndex].fontemassLR + state.cells[cellIndex].fontemassCR) <= (*state.globals).localtiny * 1e-5)) {
        state.cells[cellIndex].term1 = 0.;
        state.cells[cellIndex].term2 = 0.;
        bif = 0;
        state.cells[cellIndex].c0 = 1 + 0 * distributionCoefficient;
        state.cells[cellIndex].ud = 0 * driftVelocity;
        state.cells[cellIndex].arranjo = 0;
    } else if (state.cells[cellIndex - 1].acsr.tipo == 5 && state.cells[cellIndex - 1].acsr.chk.AreaGarg <= (1e-3)) {
        state.cells[cellIndex].term1 = 0.;
        state.cells[cellIndex].term2 = 0.;
        bif = 0;
        state.cells[cellIndex].c0 = 1 + 0 * distributionCoefficient;
        state.cells[cellIndex].ud = 0 * driftVelocity;
        state.cells[cellIndex].arranjo = 0;
    }

    else if (superficialGasVelocity >= 0 && state.cells[cellIndex - 1].alfPigD <= (*state.globals).localtiny && (state.cells[cellIndex].fontemassGL <= (*state.globals).localtiny * 1e-5 && state.cells[cellIndex].fontemassGR <= (*state.globals).localtiny * 1e-5)) {
        state.cells[cellIndex].term1 = 1.;
        state.cells[cellIndex].term2 = 0.;
        bif = 0;
        state.cells[cellIndex].c0 = 1 + 0 * distributionCoefficient;
        state.cells[cellIndex].ud = 0 * driftVelocity;
        state.cells[cellIndex].arranjo = 0;
        if (fabs(superficialGasVelocity) <= 1e-15 && state.cells[cellIndex].alfPigE > (1. - (*state.globals).localtiny) && superficialLiquidVelocity < 0 && leftSuperficialLiquidVelocity < 0 && state.cells[cellIndex].duto.teta > 0) {
            state.cells[cellIndex].term1 = 0.;
            state.cells[cellIndex].term2 = 0.;
            bif = 0;
        }
        if (fabs(superficialGasVelocity) <= 1e-15 && state.cells[cellIndex].alfPigE > (1. - 10 * (*state.globals).localtiny) && state.cells[cellIndex].duto.teta < 0) {
            state.cells[cellIndex].term1 = 0.;
            state.cells[cellIndex].term2 = 0.;
            bif = 1;
        }
        if (fabs(superficialGasVelocity) <= 1e-15 && state.cells[cellIndex].alf >= state.cells[cellIndex + 1].alf && rightSuperficialGasVelocity < 0) {
            bif = 1;
        }
    } else if (superficialGasVelocity <= 0 && state.cells[cellIndex].alfPigE <= (*state.globals).localtiny && (state.cells[cellIndex].fontemassGL <= (*state.globals).localtiny * 1e-5 && state.cells[cellIndex].fontemassGR <= (*state.globals).localtiny * 1e-5)) {
        state.cells[cellIndex].term1 = 1.;
        state.cells[cellIndex].term2 = 0.;
        bif = 0;
        state.cells[cellIndex].c0 = 1 + 0 * distributionCoefficient;
        state.cells[cellIndex].ud = 0 * driftVelocity;
        state.cells[cellIndex].arranjo = 0;
        if (fabs(superficialGasVelocity) <= 1e-15 && state.cells[cellIndex - 1].alfPigD > (1. - (*state.globals).localtiny) && superficialLiquidVelocity > 0 && rightSuperficialLiquidVelocity > 0) {
            state.cells[cellIndex].term1 = 0.;
            state.cells[cellIndex].term2 = 0.;
            bif = 0;
        }
        if (fabs(superficialGasVelocity) <= 1e-15 && state.cells[cellIndex].alfL >= state.cells[cellIndex - 1].alfL && leftSuperficialGasVelocity > 0) {
            bif = 1;
        }
    } else if (superficialLiquidVelocity >= 0 && state.cells[cellIndex - 1].alfPigD >= 1. - (*state.globals).localtiny && ((state.cells[cellIndex].fontemassLL + state.cells[cellIndex].fontemassCL) <= (*state.globals).localtiny * 1e-5 && (state.cells[cellIndex].fontemassLR + state.cells[cellIndex].fontemassCR) <= (*state.globals).localtiny * 1e-5)) {
        state.cells[cellIndex].term1 = 0.;
        state.cells[cellIndex].term2 = 0.;
        bif = 0;
        state.cells[cellIndex].c0 = 1 + 0 * distributionCoefficient;
        state.cells[cellIndex].ud = 0 * driftVelocity;
        state.cells[cellIndex].arranjo = 0;
        if (fabs(superficialLiquidVelocity) <= 1e-15 && state.cells[cellIndex].alfPigE < (*state.globals).localtiny && rightSuperficialLiquidVelocity < 0) {
            state.cells[cellIndex].term1 = 1.;
            state.cells[cellIndex].term2 = 0.;
            bif = 0;
        } else if (fabs(rightSuperficialLiquidVelocity) < (*state.globals).localtiny * 1e-5) {
            if (fabs(superficialLiquidVelocity) < (*state.globals).localtiny * 1e-5 && state.cells[cellIndex].alfPigE < (*state.globals).localtiny && state.cells[cellIndex].fontemassGR >= (*state.globals).localtiny * 1e-5) {
                state.cells[cellIndex].term1 = 1.;
                state.cells[cellIndex].term2 = 0.;
                bif = 0;
            }
        } else if (fabs(superficialLiquidVelocity) < 1e-15 && rightSuperficialLiquidVelocity < 0 && ((state.cells[cellIndex].alfPigE <= (1 - 1 * (*state.globals).localtiny + .0 * state.cells[cellIndex].alfPigER) && state.cells[cellIndex].alfPigER < 1 - 1 * (*state.globals).localtiny) || state.cells[cellIndex].alfPigE <= 0.7))
            bif = 1;

    } else if (superficialLiquidVelocity <= 0 && state.cells[cellIndex].alfPigE >= 1. - (*state.globals).localtiny && ((state.cells[cellIndex].fontemassLL + state.cells[cellIndex].fontemassCL) <= (*state.globals).localtiny * 1e-5 && (state.cells[cellIndex].fontemassLR + state.cells[cellIndex].fontemassCR) <= (*state.globals).localtiny * 1e-5)) {
        state.cells[cellIndex].term1 = 0.;
        state.cells[cellIndex].term2 = 0.;
        bif = 0;
        state.cells[cellIndex].c0 = 1 + 0 * distributionCoefficient;
        state.cells[cellIndex].ud = 0 * driftVelocity;
        state.cells[cellIndex].arranjo = 0;
        if (fabs(superficialLiquidVelocity) <= (*state.globals).localtiny * 1e-5 && state.cells[cellIndex - 1].alfPigD < (*state.globals).localtiny && leftSuperficialLiquidVelocity > 0) {
            state.cells[cellIndex].term1 = 1.;
            state.cells[cellIndex].term2 = 0.;
            bif = 0;
        }

        if (fabs(leftSuperficialLiquidVelocity) < (*state.globals).localtiny * 1e-5) {
            if (fabs(superficialLiquidVelocity) < (*state.globals).localtiny * 1e-5 && state.cells[cellIndex - 1].alfPigD < (*state.globals).localtiny && state.cells[cellIndex - 1].fontemassGR >= (*state.globals).localtiny * 1e-5) {
                state.cells[cellIndex].term1 = 1.;
                state.cells[cellIndex].term2 = 0.;
                bif = 0;
            }
        } else if (cellIndex > 1 && fabs(superficialLiquidVelocity) < (*state.globals).localtiny * 1e-5 && state.cells[cellIndex].duto.teta < 0.95 * M_PI / 2. && leftSuperficialLiquidVelocity > 0 && ((state.cells[cellIndex - 1].alfPigD <= 1.0 * state.cells[cellIndex - 2].alfPigD && state.cells[cellIndex - 2].alfPigD < 0.99) || state.cells[cellIndex - 1].alfPigD < 0.7))
            bif = 1;
        else if (cellIndex > 1 && fabs(superficialLiquidVelocity) < (*state.globals).localtiny * 1e-5 && state.cells[cellIndex].duto.teta >= 0.95 * M_PI / 2. && leftSuperficialLiquidVelocity > 0 && ((state.cells[cellIndex - 1].alfPigD <= 1.0 * state.cells[cellIndex - 2].alfPigD && state.cells[cellIndex - 2].alfPigD < 0.99) || state.cells[cellIndex - 1].alfPigD < 0.7))
            bif = 1;
        else {
            if (fabs(superficialLiquidVelocity) < (*state.globals).localtiny * 1e-5 && state.cells[cellIndex].duto.teta < 0.95 * M_PI / 2. && leftSuperficialLiquidVelocity > 0 && ((state.cells[cellIndex - 1].alfPigD <= 1.0 * state.cells[cellIndex - 1].alfL) || state.cells[cellIndex - 1].alfPigD < 0.7))
                bif = 1;
            else if (fabs(superficialLiquidVelocity) < (*state.globals).localtiny * 1e-5 && state.cells[cellIndex].duto.teta >= 0.95 * M_PI / 2. && leftSuperficialLiquidVelocity > 0 && ((state.cells[cellIndex - 1].alfPigD <= 1.0 * state.cells[cellIndex - 1].alfL) || state.cells[cellIndex - 1].alfPigD < 0.7))
                bif = 1;
        }
    }

    if (superficialGasVelocity < 0 && state.cells[cellIndex].alf > (*state.globals).localtiny && state.cells[cellIndex].alf < 1. - (*state.globals).localtiny)
        bif = 1;

    if (bif == 1) {

        double c0;
        double ud;
        double meanVoidFraction;

        meanVoidFraction = state.cells[cellIndex - 1].alfPigD;
        if (superficialGasVelocity < 0)
            meanVoidFraction = state.cells[cellIndex].alfPigE;
        c0 = 1.2;
        double meanDiameter = state.cells[cellIndex].duto.a;
        if (state.cells[cellIndex].MC >= 0)
            meanDiameter = state.cells[cellIndex].dutoL.a;
        double inclinationSign = 1.;
        if (state.cells[cellIndex].duto.teta < 0.)
            inclinationSign = -1.;
        ud = inclinationSign * 0.32 * sqrt(9.82 * meanDiameter);
        if (fabs(rightGasDensity) / rightLiquidDensity > 0.9) {
            c0 = 1.;
            ud = 0.;
        }
        state.closureUpdater.buffered(cellIndex, c0, ud);
        if (state.input.escorregamentoCelulaContorno == 0) {
            c0 = 1.;
            ud = 0.;
        }
        double numerator = (1. - meanVoidFraction * c0);
        double denominator = 1 + c0 * meanVoidFraction * ((gasDensity / liquidDensity) - 1.);
        state.cells[cellIndex].term1 = numerator / denominator;
        state.cells[cellIndex].term2 = (-flowArea * meanVoidFraction * gasDensity * ud) / denominator;
        if (state.cells[cellIndex - 1].acsr.tipo == 5 && state.cells[cellIndex - 1].acsr.chk.AreaGarg <= (1e-3 + state.input.master1.razareaativ) * state.cells[cellIndex - 1].duto.area) {
            state.cells[cellIndex].term1 = 0.;
            state.cells[cellIndex].term2 = 0.;
        }
    }
}

void selectAndApplyBufferedInletFlowRegime(
    const ThermalState &state, int cellIndex, double distributionCoefficient, double driftVelocity,
    double superficialGasVelocity, double superficialLiquidVelocity, double rightSuperficialLiquidVelocity, double rightGasDensity, double rightLiquidDensity,
    double gasDensity, double liquidDensity, double flowArea) {
    int bif = 1;

    if ((*state.globals).lixo5 > 29900) {
        int debugStop;
        debugStop = 0;
    }

    if (state.inletVoidFraction < (*state.globals).localtiny && state.cells[cellIndex].alfPigE <= (*state.globals).localtiny && state.cells[cellIndex].fontemassGL <= (*state.globals).localtiny * 1e-5 && state.cells[cellIndex].fontemassGR <= (*state.globals).localtiny * 1e-5) {
        state.cells[cellIndex].term1 = 1.;
        state.cells[cellIndex].term2 = 0.;
        bif = 0;
        state.cells[cellIndex].c0 = 1 + 0 * distributionCoefficient;
        state.cells[cellIndex].ud = 0 * driftVelocity;
        state.cells[cellIndex].arranjo = 0;
    } else if (state.inletVoidFraction >= (1. - (*state.globals).localtiny) && state.cells[cellIndex].alfPigE >= (1. - (*state.globals).localtiny) && ((state.cells[cellIndex].fontemassLL + state.cells[cellIndex].fontemassCL) <= (*state.globals).localtiny * 1e-5 && (state.cells[cellIndex].fontemassLR + state.cells[cellIndex].fontemassCR) <= (*state.globals).localtiny * 1e-5)) {
        state.cells[cellIndex].term1 = 0.;
        state.cells[cellIndex].term2 = 0.;
        bif = 0;
        state.cells[cellIndex].c0 = 1 + 0 * distributionCoefficient;
        state.cells[cellIndex].ud = 0 * driftVelocity;
        state.cells[cellIndex].arranjo = 0;
    } else if (superficialGasVelocity >= 0 && state.inletVoidFraction <= (*state.globals).localtiny && (state.cells[cellIndex].fontemassGL <= (*state.globals).localtiny * 1e-5 && state.cells[cellIndex].fontemassGR <= (*state.globals).localtiny * 1e-5)) {
        state.cells[cellIndex].term1 = 1.;
        state.cells[cellIndex].term2 = 0.;
        bif = 0;
        state.cells[cellIndex].c0 = 1 + 0 * distributionCoefficient;
        state.cells[cellIndex].ud = 0 * driftVelocity;
        state.cells[cellIndex].arranjo = 0;
        if (fabs(superficialGasVelocity) <= 1e-15 && state.cells[cellIndex].alfPigE > (1. - (*state.globals).localtiny) && superficialLiquidVelocity < 0 && state.cells[cellIndex].duto.teta > 0) {
            state.cells[cellIndex].term1 = 0.;
            state.cells[cellIndex].term2 = 0.;
            bif = 0;
        }
        if (fabs(superficialGasVelocity) <= 1e-15 && state.cells[cellIndex].alfPigE > (1. - 1 * (*state.globals).localtiny) && state.cells[cellIndex].duto.teta < 0 && superficialLiquidVelocity > 0) {
            state.cells[cellIndex].term1 = 0.;
            state.cells[cellIndex].term2 = 0.;
            bif = 1;
        }
        if (fabs(superficialGasVelocity) <= 1e-15 && state.cells[cellIndex].alf >= state.inletVoidFraction && superficialLiquidVelocity < 0) {
            bif = 1;
        }
    } else if (superficialGasVelocity <= 0 && state.cells[cellIndex].alfPigE <= (*state.globals).localtiny && (state.cells[cellIndex].fontemassGL <= (*state.globals).localtiny * 1e-5 && state.cells[cellIndex].fontemassGR <= (*state.globals).localtiny * 1e-5)) {
        state.cells[cellIndex].term1 = 1.;
        state.cells[cellIndex].term2 = 0.;
        bif = 0;
        state.cells[cellIndex].c0 = 1 + 0 * distributionCoefficient;
        state.cells[cellIndex].ud = 0 * driftVelocity;
        state.cells[cellIndex].arranjo = 0;
        if (fabs(superficialGasVelocity) <= 1e-15 && state.inletVoidFraction > (1. - (*state.globals).localtiny) && superficialLiquidVelocity > 0) {
            state.cells[cellIndex].term1 = 0.;
            state.cells[cellIndex].term2 = 0.;
            bif = 0;
        }
        if (fabs(superficialGasVelocity) <= 1e-15 && state.inletVoidFraction > (*state.globals).localtiny && superficialLiquidVelocity > 0) {
            bif = 1;
        }
    } else if (superficialLiquidVelocity >= 0 && state.inletVoidFraction >= 1. - (*state.globals).localtiny && ((state.cells[cellIndex].fontemassLL + state.cells[cellIndex].fontemassCL) <= (*state.globals).localtiny * 1e-5 && (state.cells[cellIndex].fontemassLR + state.cells[cellIndex].fontemassCR) <= (*state.globals).localtiny * 1e-5)) {
        state.cells[cellIndex].term1 = 0.;
        state.cells[cellIndex].term2 = 0.;
        bif = 0;
        state.cells[cellIndex].c0 = 1 + 0 * distributionCoefficient;
        state.cells[cellIndex].ud = 0 * driftVelocity;
        state.cells[cellIndex].arranjo = 0;
        if (fabs(superficialLiquidVelocity) <= 1e-15 && state.cells[cellIndex].alfPigE < (*state.globals).localtiny && rightSuperficialLiquidVelocity < 0) {
            state.cells[cellIndex].term1 = 1.;
            state.cells[cellIndex].term2 = 0.;
            bif = 0;
        } else if (fabs(rightSuperficialLiquidVelocity) < (*state.globals).localtiny * 1e-5) {
            if (fabs(superficialLiquidVelocity) < (*state.globals).localtiny * 1e-5 && state.cells[cellIndex].alfPigE < (*state.globals).localtiny && state.cells[cellIndex].fontemassGR >= (*state.globals).localtiny * 1e-5) {
                state.cells[cellIndex].term1 = 1.;
                state.cells[cellIndex].term2 = 0.;
                bif = 0;
            }
        } else if (fabs(superficialLiquidVelocity) < 1e-15 && rightSuperficialLiquidVelocity < 0 && ((state.cells[cellIndex].alfPigE <= (1 - 1 * (*state.globals).localtiny + .0 * state.cells[cellIndex].alfPigER) && state.cells[cellIndex].alfPigER < 1 - 1 * (*state.globals).localtiny) || state.cells[cellIndex].alfPigE <= 0.7))
            bif = 1;

    } else if (superficialLiquidVelocity <= 0 && state.cells[cellIndex].alfPigE >= 1. - (*state.globals).localtiny && ((state.cells[cellIndex].fontemassLL + state.cells[cellIndex].fontemassCL) <= (*state.globals).localtiny * 1e-5 && (state.cells[cellIndex].fontemassLR + state.cells[cellIndex].fontemassCR) <= (*state.globals).localtiny * 1e-5)) {
        state.cells[cellIndex].term1 = 0.;
        state.cells[cellIndex].term2 = 0.;
        bif = 0;
        state.cells[cellIndex].c0 = 1 + 0 * distributionCoefficient;
        state.cells[cellIndex].ud = 0 * driftVelocity;
        state.cells[cellIndex].arranjo = 0;
        if (fabs(superficialLiquidVelocity) <= (*state.globals).localtiny * 1e-5 && state.inletVoidFraction < (*state.globals).localtiny && rightSuperficialLiquidVelocity < 0) {
            state.cells[cellIndex].term1 = 1.;
            state.cells[cellIndex].term2 = 0.;
            bif = 0;
        }

        if (fabs(rightSuperficialLiquidVelocity) < (*state.globals).localtiny * 1e-5) {
            if (fabs(superficialLiquidVelocity) < (*state.globals).localtiny * 1e-5 && state.inletVoidFraction < (*state.globals).localtiny) {
                state.cells[cellIndex].term1 = 1.;
                state.cells[cellIndex].term2 = 0.;
                bif = 0;
            }
        } else {
            if (fabs(superficialLiquidVelocity) < (*state.globals).localtiny * 1e-5 && state.cells[cellIndex].duto.teta < 0.95 * M_PI / 2. && superficialGasVelocity < 0 && (state.inletVoidFraction < 0.7))
                bif = 1;
            else if (fabs(superficialLiquidVelocity) < (*state.globals).localtiny * 1e-5 && state.cells[cellIndex].duto.teta >= 0.95 * M_PI / 2. && superficialGasVelocity < 0 && (state.inletVoidFraction < 0.7))
                bif = 1;
        }
    }
    if (superficialLiquidVelocity > 0 && fabs(superficialGasVelocity) <= (*state.globals).localtiny * 1e-5 && state.inletVoidFraction > (*state.globals).localtiny && state.inletVoidFraction < 1 - (*state.globals).localtiny)
        bif = 1;
    if (superficialGasVelocity > 0 && fabs(superficialLiquidVelocity) <= (*state.globals).localtiny * 1e-5 && state.inletVoidFraction > (*state.globals).localtiny && state.inletVoidFraction < 1 - (*state.globals).localtiny)
        bif = 1;

    if (bif == 1) {

        double c0;
        double ud;
        double meanVoidFraction;

        meanVoidFraction = state.inletVoidFraction;
        if (superficialGasVelocity < 0)
            meanVoidFraction = state.cells[cellIndex].alfPigE;
        c0 = 1.2;
        double meanDiameter = state.cells[cellIndex].duto.a;
        if (state.cells[cellIndex].MC >= 0)
            meanDiameter = state.cells[cellIndex].dutoL.a;
        double inclinationSign = 1.;
        if (state.cells[cellIndex].duto.teta < 0.)
            inclinationSign = -1.;
        ud = inclinationSign * 0.32 * sqrt(9.82 * meanDiameter);
        if (fabs(rightGasDensity) / rightLiquidDensity > 0.9) {
            c0 = 1.;
            ud = 0.;
        }
        state.closureUpdater.bufferedInitialization(cellIndex, c0, ud);

        double numerator = (1. - meanVoidFraction * c0);
        double denominator = 1 + c0 * meanVoidFraction * ((gasDensity / liquidDensity) - 1.);
        state.cells[cellIndex].term1 = numerator / denominator;
        state.cells[cellIndex].term2 = (-flowArea * meanVoidFraction * gasDensity * ud) / denominator;
    }
}

}  // namespace

void updateFlowPartitionTerms(const ThermalState &state, int aflu) {
    // #pragma omp parallel for num_threads(numthreads)
    aflu = 0;
    Vcr<int> bifurcationFlags(state.lastCell + 1, 0);
    Vcr<int> valveFlags(state.lastCell + 1, 1);
    int lastFlowCell = state.lastCell;
    if (aflu == 1)
        lastFlowCell = state.lastCell + 1;
    for (int cellIndex = 0; cellIndex <= state.lastCell; cellIndex++) {
        state.cells[cellIndex].c0ini = state.cells[cellIndex].c0;
        state.cells[cellIndex].udini = state.cells[cellIndex].ud;
        if (cellIndex == state.lastCell && aflu == 1 && (*state.globals).lixo5 >= 1560) {
            int debugStop;
            debugStop = 0;
        }
        if (cellIndex != 0 && cellIndex != lastFlowCell) {
            updateInteriorFlowPartitionCell(state, cellIndex, bifurcationFlags, valveFlags);
        } else if (cellIndex == 0) {
            if (state.input.ConContEntrada == 0) {
                state.cells[1].alfL = state.cells[0].alf;
                state.cells[0].alfL = state.cells[0].alf;
                state.cells[1].betL = state.cells[0].bet;
                state.cells[0].betL = state.cells[0].bet;
                state.cells[0].term1 = 0.;
                state.cells[0].term2 = 0.;
                state.cells[0].term1L = 0.;
                state.cells[0].term2L = 0.;
            } else {
                if (state.inletMassFraction < 1) {
                    int debugStop;
                    debugStop = 0;
                }

                double lengthRatio = 0.5;
                double meanPressure;
                meanPressure = state.inletPressure;

                double meanTemperature;
                if (state.cells[0].QL < 0.)
                    meanTemperature = state.cells[cellIndex].temp;
                else
                    meanTemperature = state.inletTemperature;

                double gasDensity = state.cells[0].flui.MasEspGas(state.inletPressure, state.inletTemperature);
                double liquidDensity = state.cells[0].flui.MasEspLiq(state.inletPressure, state.inletTemperature);
                double rcis = state.cells[0].fluicol.MasEspFlu(state.inletPressure, state.inletTemperature);

                double liquidMixtureDensity = state.inletComposition * rcis + (1 - state.inletComposition) * liquidDensity;
                state.inletVoidFraction = (-state.inletMassFraction * liquidMixtureDensity / (state.inletMassFraction * gasDensity - gasDensity - state.inletMassFraction * liquidMixtureDensity)) / (state.cells[0].c0);

                double betI;
                double liquidViscosity;
                double surfaceTension;
                if (state.cells[cellIndex].QL < 0.) { // testeBeta
                    betI = state.cells[cellIndex].betPigE;
                    liquidDensity = (1 - betI) * state.cells[cellIndex].flui.MasEspLiq(meanPressure, meanTemperature) + betI * state.cells[cellIndex].fluicol.MasEspFlu(meanPressure, meanTemperature);
                    liquidViscosity = (1 - betI) * state.cells[cellIndex].flui.ViscOleo(meanPressure, meanTemperature) + betI * state.cells[cellIndex].fluicol.VisFlu(meanPressure, meanTemperature);
                    surfaceTension = (1 - betI) * state.cells[cellIndex].flui.TensSuper(meanPressure, meanTemperature) + betI * state.cells[cellIndex].fluicol.TensSuper(meanPressure, meanTemperature);
                } else {
                    betI = state.inletComposition;
                    liquidDensity = (1 - betI) * (*state.cells[cellIndex].fluiL).MasEspLiq(meanPressure, meanTemperature) + betI * state.cells[cellIndex].fluicol.MasEspFlu(meanPressure, meanTemperature);
                    liquidViscosity = (1 - betI) * (*state.cells[cellIndex].fluiL).ViscOleo(meanPressure, meanTemperature) + betI * state.cells[cellIndex].fluicol.VisFlu(meanPressure, meanTemperature);
                    surfaceTension = (1 - betI) * (*state.cells[cellIndex].fluiL).TensSuper(meanPressure, meanTemperature) + betI * state.cells[cellIndex].fluicol.TensSuper(meanPressure, meanTemperature);
                }

                double gasViscosity;
                double flowArea;
                double noSlipLiquidHoldup;
                if (state.cells[cellIndex].QG >= 0) {
                    flowArea = state.cells[cellIndex].duto.area;
                    gasDensity = (*state.cells[cellIndex].fluiL).MasEspGas(meanPressure, meanTemperature);
                    gasViscosity = (*state.cells[cellIndex].fluiL).ViscGas(meanPressure, meanTemperature);
                    noSlipLiquidHoldup = 1. - state.inletVoidFraction;
                } else {
                    gasDensity = state.cells[cellIndex].flui.MasEspGas(meanPressure, meanTemperature);
                    gasViscosity = state.cells[cellIndex].flui.ViscGas(meanPressure, meanTemperature);
                    flowArea = state.cells[cellIndex].duto.area;
                    noSlipLiquidHoldup = 1. - state.cells[cellIndex].alf;
                }
                double superficialGasVelocity = state.cells[cellIndex].QG / (flowArea);
                double superficialLiquidVelocity = state.cells[cellIndex].QL / (flowArea);
                double diameter = state.cells[cellIndex].duto.a;

                double mixtureDensity = noSlipLiquidHoldup * liquidDensity + (1 - noSlipLiquidHoldup) * gasDensity;
                double mixtureViscosity = (noSlipLiquidHoldup * liquidViscosity + (1 - noSlipLiquidHoldup) * gasViscosity) / pow(10., 3.);
                double inclinationAngle = state.cells[cellIndex].duto.teta;
                double inclinationSign = 1.;
                if (inclinationAngle < 0.)
                    inclinationSign = -1.;
                double c0 = 2.;
                double ud = inclinationSign * 0.0246 * cos(inclinationAngle) + 1.606 * pow(9.82 * surfaceTension * (liquidDensity - gasDensity) / (liquidDensity * liquidDensity), 0.25) * sin(inclinationAngle);

                double rightFlowArea = state.cells[cellIndex].dutoR.area;
                double rightLengthRatio = state.cells[cellIndex].dx / (state.cells[cellIndex].dxR + state.cells[cellIndex].dx);
                double rightMeanPressure = state.cells[cellIndex].presauxR;
                double rightMeanTemperature = state.cells[cellIndex].temp * rightLengthRatio + state.cells[cellIndex].tempL * (1. - rightLengthRatio);
                double betIR = state.cells[cellIndex].betPigD;
                if (state.cells[cellIndex].QLR < 0.) // testeBeta
                    betIR = state.cells[cellIndex + 1].betPigE;

                double rightGasDensity = state.cells[cellIndex].flui.MasEspGas(rightMeanPressure, rightMeanTemperature);
                double rightLiquidDensity = (1 - betIR) * state.cells[cellIndex].flui.MasEspLiq(rightMeanPressure, rightMeanTemperature) + betIR * state.cells[cellIndex].fluicol.MasEspFlu(rightMeanPressure, rightMeanTemperature);
                double rightSuperficialGasVelocity = (state.cells[cellIndex].MR - state.cells[cellIndex].MliqiniR) / (rightGasDensity * rightFlowArea);
                double rightSuperficialLiquidVelocity = (state.cells[cellIndex].MliqiniR) / (rightLiquidDensity * rightFlowArea);

                selectAndApplyInletBoundaryFlowRegime(
                    state, cellIndex, bifurcationFlags, c0, ud, superficialGasVelocity, superficialLiquidVelocity, rightSuperficialLiquidVelocity,
                    rightGasDensity, rightLiquidDensity, gasDensity, liquidDensity, flowArea);
                state.cells[1].term1L = state.cells[cellIndex].term1;
                state.cells[1].term2L = state.cells[cellIndex].term2;

                state.cells[1].alfL = state.cells[0].alf;
                state.cells[0].alfL = state.inletVoidFraction;
                state.cells[1].betL = state.cells[0].bet;
                state.cells[0].betL = state.inletComposition;
            }

        } else if (aflu == 0) {
            updateOutletBoundaryFlowPartition(state, cellIndex, bifurcationFlags);
        }
    }

    finalizeFlowPartitionTerms(state, bifurcationFlags, valveFlags);
}
void updateOutletFlowPartitionTerms(const ThermalState &state) {

    int cellIndex = state.lastCell;

    double lengthRatio = state.cells[cellIndex].dx / (state.cells[cellIndex].dx + state.cells[cellIndex].dxL);
    double meanPressure = state.cells[cellIndex].presBuf;
    double meanTemperature = state.gasSurfaceTemperature;
    meanTemperature = state.cells[cellIndex - 1].temp;
    if (state.cells[cellIndex].VTemper < 0.)
        meanTemperature = state.gasSurfaceTemperature;
    double betI = state.cells[cellIndex - 1].betPigD;
    double liquidDensity;
    double liquidViscosity;
    double surfaceTension;
    if (state.cells[cellIndex].MliqiniBuf < 0.) { // testeBeta
        betI = state.cells[cellIndex].betPigE;
        liquidDensity = (1 - betI) * state.cells[cellIndex].flui.MasEspLiq(meanPressure, meanTemperature) + betI * state.cells[cellIndex].fluicol.MasEspFlu(meanPressure, meanTemperature);
        liquidViscosity = (1 - betI) * state.cells[cellIndex].flui.ViscOleo(meanPressure, meanTemperature) + betI * state.cells[cellIndex].fluicol.VisFlu(meanPressure, meanTemperature);
        surfaceTension = (1 - betI) * state.cells[cellIndex].flui.TensSuper(meanPressure, meanTemperature) + betI * state.cells[cellIndex].fluicol.TensSuper(meanPressure, meanTemperature);
    } else {
        betI = state.cells[cellIndex - 1].betPigD;
        liquidDensity = (1 - betI) * state.cells[cellIndex - 1].flui.MasEspLiq(meanPressure, meanTemperature) + betI * state.cells[cellIndex - 1].fluicol.MasEspFlu(meanPressure, meanTemperature);
        liquidViscosity = (1 - betI) * state.cells[cellIndex - 1].flui.ViscOleo(meanPressure, meanTemperature) + betI * state.cells[cellIndex - 1].fluicol.VisFlu(meanPressure, meanTemperature);
        surfaceTension = (1 - betI) * state.cells[cellIndex - 1].flui.TensSuper(meanPressure, meanTemperature) + betI * state.cells[cellIndex - 1].fluicol.TensSuper(meanPressure, meanTemperature);
    }
    double gasDensity;
    double gasViscosity;
    double flowArea;
    double noSlipLiquidHoldup;
    if ((state.cells[cellIndex].MCBuf - state.cells[cellIndex].MliqiniBuf) >= 0) {
        flowArea = state.cells[cellIndex].dutoL.area;
        noSlipLiquidHoldup = 1. - state.cells[cellIndex].alfL;
        gasDensity = state.cells[cellIndex - 1].flui.MasEspGas(meanPressure, meanTemperature);
        gasViscosity = state.cells[cellIndex - 1].flui.ViscGas(meanPressure, meanTemperature);
    } else {
        gasDensity = state.cells[cellIndex].flui.MasEspGas(meanPressure, meanTemperature);
        gasViscosity = state.cells[cellIndex].flui.ViscGas(meanPressure, meanTemperature);
        flowArea = state.cells[cellIndex].duto.area;
        noSlipLiquidHoldup = 1. - state.cells[cellIndex].alf;
    }
    double superficialGasVelocity = (state.cells[cellIndex].MCBuf - state.cells[cellIndex].MliqiniBuf) / (gasDensity * flowArea);
    double superficialLiquidVelocity = (state.cells[cellIndex].MliqiniBuf) / (liquidDensity * flowArea);
    double diameter = state.cells[cellIndex].duto.a;
    if (superficialGasVelocity >= 0)
        diameter = state.cells[cellIndex - 1].duto.a;

    double mixtureDensity = noSlipLiquidHoldup * liquidDensity + (1 - noSlipLiquidHoldup) * gasDensity;
    double mixtureViscosity = (noSlipLiquidHoldup * liquidViscosity + (1 - noSlipLiquidHoldup) * gasViscosity) / pow(10., 3.);
    double mixtureReynolds = diameter * mixtureDensity * (fabs(superficialGasVelocity) / flowArea + fabs(superficialLiquidVelocity) / flowArea) / mixtureViscosity;
    double inclinationAngle = state.cells[cellIndex].duto.teta;
    if (cellIndex >= 2) {
        if (state.cells[cellIndex - 2].acsr.tipo == 5 && state.cells[cellIndex - 2].acsr.chk.AreaGarg <= (1e-3)) {
            if ((state.cells[cellIndex].MCBuf - state.cells[cellIndex].MliqiniBuf) >= 0)
                inclinationAngle = state.cells[cellIndex].duto.teta;
            else
                inclinationAngle = state.cells[cellIndex].dutoR.teta;
        } else {
            if ((state.cells[cellIndex].MCBuf - state.cells[cellIndex].MliqiniBuf) >= 0)
                inclinationAngle = state.cells[cellIndex].dutoL.teta;
            else
                inclinationAngle = state.cells[cellIndex].duto.teta;
        }
    }
    double inclinationSign = 1.;
    if (inclinationAngle < 0.)
        inclinationSign = -1.;
    double c0 = 2.;
    double ud = inclinationSign * 0.0246 * cos(inclinationAngle) + 1.606 * pow(9.82 * surfaceTension * (liquidDensity - gasDensity) / (liquidDensity * liquidDensity), 0.25) * sin(inclinationAngle);

    double leftFlowArea = state.cells[cellIndex].dutoL.area;
    double leftLengthRatio = state.cells[cellIndex - 1].dxL / (state.cells[cellIndex - 1].dx + state.cells[cellIndex - 1].dxL);
    double leftMeanPressure = state.cells[cellIndex - 1].presaux;
    double leftMeanTemperature = state.cells[cellIndex - 1].temp * leftLengthRatio + state.cells[cellIndex - 1].tempL * (1. - leftLengthRatio);
    double betIL = 0.;
    if (cellIndex < 2)
        betIL = state.cells[cellIndex - 1].betL;
    else
        betIL = state.cells[cellIndex - 2].betPigD;
    if (state.cells[cellIndex - 1].MliqiniBuf < 0.)
        betIL = state.cells[cellIndex - 1].betPigE; // testeBeta
    // betIL = state.cells[cellIndex - 1].betPigE;    //duvidabeta
    double leftGasDensity = state.cells[cellIndex].flui.MasEspGas(leftMeanPressure, leftMeanTemperature);
    double leftLiquidDensity = (1 - betIL) * state.cells[cellIndex].flui.MasEspLiq(leftMeanPressure, leftMeanTemperature) + betIL * state.cells[cellIndex].fluicol.MasEspFlu(leftMeanPressure, leftMeanTemperature);
    double leftSuperficialGasVelocity = (state.cells[cellIndex].MLBuf - state.cells[cellIndex].MliqiniLBuf) / (leftGasDensity * leftFlowArea);
    double leftSuperficialLiquidVelocity = (state.cells[cellIndex].MliqiniLBuf) / (leftLiquidDensity * leftFlowArea);

    double rightFlowArea = state.cells[cellIndex].dutoR.area;
    double rightLengthRatio = state.cells[cellIndex].dxR / (state.cells[cellIndex].dxR + state.cells[cellIndex].dx);
    double rightMeanPressure = state.cells[cellIndex].presRBuf;
    double rightMeanTemperature = state.cells[cellIndex].temp * rightLengthRatio + state.cells[cellIndex].tempR * (1. - rightLengthRatio);
    double betIR = state.cells[cellIndex].betPigD;
    if (state.cells[cellIndex].MliqiniRBuf < 0.) { // testeBeta
        if (cellIndex > state.lastCell - 2)
            betIR = state.cells[cellIndex].betR;
        else
            betIR = state.cells[cellIndex + 1].betPigE;
    }
    double rightGasDensity = state.cells[cellIndex].flui.MasEspGas(rightMeanPressure, rightMeanTemperature);
    double rightLiquidDensity = (1 - betIR) * state.cells[cellIndex].flui.MasEspLiq(rightMeanPressure, rightMeanTemperature) + betIR * state.cells[cellIndex].fluicol.MasEspFlu(rightMeanPressure, rightMeanTemperature);
    double rightSuperficialGasVelocity = (state.cells[cellIndex].MRBuf - state.cells[cellIndex].MliqiniRBuf) / (rightGasDensity * rightFlowArea);
    double rightSuperficialLiquidVelocity = (state.cells[cellIndex].MliqiniRBuf) / (rightLiquidDensity * rightFlowArea);

    selectAndApplyBufferedOutletFlowRegime(
        state, cellIndex, c0, ud, superficialGasVelocity, superficialLiquidVelocity, leftSuperficialGasVelocity, leftSuperficialLiquidVelocity, rightSuperficialGasVelocity, rightSuperficialLiquidVelocity,
        rightGasDensity, rightLiquidDensity, gasDensity, liquidDensity, flowArea);
    state.cells[cellIndex].term1L = state.cells[cellIndex - 1].term1;
    state.cells[cellIndex].term2L = state.cells[cellIndex - 1].term2;
    state.cells[cellIndex - 1].term1R = state.cells[cellIndex].term1;
    state.cells[cellIndex - 1].term2R = state.cells[cellIndex].term2;
}
void updateInletFlowPartitionTerms(const ThermalState &state) {

    if (state.inletMassFraction < 1) {
        int debugStop;
        debugStop = 0;
    }

    int cellIndex = 0;

    double lengthRatio = 0.5;
    double meanPressure;
    meanPressure = state.inletPressure;

    double meanTemperature;
    if (state.cells[0].MliqiniBuf < 0.)
        meanTemperature = state.cells[cellIndex].temp;
    else
        meanTemperature = state.inletTemperature;

    double betI;
    double liquidViscosity;
    double surfaceTension;
    double gasDensity = state.cells[0].flui.MasEspGas(state.inletPressure, state.inletTemperature);
    double liquidDensity = state.cells[0].flui.MasEspLiq(state.inletPressure, state.inletTemperature);
    double rcis = state.cells[0].fluicol.MasEspFlu(state.inletPressure, state.inletTemperature);

    double liquidMixtureDensity = state.inletComposition * rcis + (1 - state.inletComposition) * liquidDensity;
    state.inletVoidFraction = (-state.inletMassFraction * liquidMixtureDensity / (state.inletMassFraction * gasDensity - gasDensity - state.inletMassFraction * liquidMixtureDensity)) / (state.cells[0].c0);

    if ((state.cells[cellIndex].MCBuf - state.cells[0].MliqiniBuf) * 0 + 1 * state.cells[0].MliqiniBuf < 0.) { // duvidabeta
        betI = state.cells[cellIndex].betPigE;
        liquidDensity = (1 - betI) * state.cells[cellIndex].flui.MasEspLiq(meanPressure, meanTemperature) + betI * state.cells[cellIndex].fluicol.MasEspFlu(meanPressure, meanTemperature);
        liquidViscosity = (1 - betI) * state.cells[cellIndex].flui.ViscOleo(meanPressure, meanTemperature) + betI * state.cells[cellIndex].fluicol.VisFlu(meanPressure, meanTemperature);
        surfaceTension = (1 - betI) * state.cells[cellIndex].flui.TensSuper(meanPressure, meanTemperature) + betI * state.cells[cellIndex].fluicol.TensSuper(meanPressure, meanTemperature);
    } else {
        betI = state.inletComposition;
        liquidDensity = (1 - betI) * (*state.cells[cellIndex].fluiL).MasEspLiq(meanPressure, meanTemperature) + betI * state.cells[cellIndex].fluicol.MasEspFlu(meanPressure, meanTemperature);
        liquidViscosity = (1 - betI) * (*state.cells[cellIndex].fluiL).ViscOleo(meanPressure, meanTemperature) + betI * state.cells[cellIndex].fluicol.VisFlu(meanPressure, meanTemperature);
        surfaceTension = (1 - betI) * (*state.cells[cellIndex].fluiL).TensSuper(meanPressure, meanTemperature) + betI * state.cells[cellIndex].fluicol.TensSuper(meanPressure, meanTemperature);
    }
    double gasViscosity;
    double flowArea;
    double noSlipLiquidHoldup;
    if (state.cells[cellIndex].MCBuf - state.cells[0].MliqiniBuf >= 0) {
        flowArea = state.cells[cellIndex].dutoL.area;
        gasDensity = (*state.cells[cellIndex].fluiL).MasEspGas(meanPressure, meanTemperature);
        gasViscosity = (*state.cells[cellIndex].fluiL).ViscGas(meanPressure, meanTemperature);
        noSlipLiquidHoldup = 1. - state.inletVoidFraction;
    } else {
        gasDensity = state.cells[cellIndex].flui.MasEspGas(meanPressure, meanTemperature);
        gasViscosity = state.cells[cellIndex].flui.ViscGas(meanPressure, meanTemperature);
        flowArea = state.cells[cellIndex].duto.area;
        noSlipLiquidHoldup = 1. - state.cells[cellIndex].alf;
    }
    double superficialGasVelocity = (state.cells[cellIndex].MCBuf - state.cells[cellIndex].MliqiniBuf) / (gasDensity * flowArea);
    double superficialLiquidVelocity = state.cells[cellIndex].MliqiniBuf / (liquidDensity * flowArea);
    double diameter = state.cells[cellIndex].duto.a;
    if (superficialGasVelocity >= 0)
        diameter = state.cells[cellIndex].duto.a;

    double mixtureDensity = noSlipLiquidHoldup * liquidDensity + (1 - noSlipLiquidHoldup) * gasDensity;
    double mixtureViscosity = (noSlipLiquidHoldup * liquidViscosity + (1 - noSlipLiquidHoldup) * gasViscosity) / pow(10., 3.);
    double inclinationAngle = state.cells[cellIndex].duto.teta;
    double inclinationSign = 1.;
    if (inclinationAngle < 0.)
        inclinationSign = -1.;
    double c0 = 2.;
    double ud = inclinationSign * 0.0246 * cos(inclinationAngle) + 1.606 * pow(9.82 * surfaceTension * (liquidDensity - gasDensity) / (liquidDensity * liquidDensity), 0.25) * sin(inclinationAngle);

    double rightFlowArea = state.cells[cellIndex].dutoR.area;
    double rightLengthRatio = state.cells[cellIndex].dxR / (state.cells[cellIndex].dxR + state.cells[cellIndex].dx);
    double rightMeanPressure = state.cells[cellIndex].presRBuf * rightLengthRatio + state.cells[cellIndex].presBuf * (1. - rightLengthRatio);
    double rightMeanTemperature = state.cells[cellIndex].temp * rightLengthRatio + state.cells[cellIndex].tempR * (1. - rightLengthRatio);
    double betIR = state.cells[cellIndex].betPigD;
    if (state.cells[cellIndex].QLR < 0.) // testeBeta
        betIR = state.cells[cellIndex + 1].betPigE;

    double rightGasDensity = state.cells[cellIndex].flui.MasEspGas(rightMeanPressure, rightMeanTemperature);
    double rightLiquidDensity = (1 - betIR) * state.cells[cellIndex].flui.MasEspLiq(rightMeanPressure, rightMeanTemperature) + betIR * state.cells[cellIndex].fluicol.MasEspFlu(rightMeanPressure, rightMeanTemperature);
    double rightSuperficialGasVelocity = (state.cells[cellIndex].MRBuf - state.cells[cellIndex].MliqiniRBuf) / (rightGasDensity * rightFlowArea);
    double rightSuperficialLiquidVelocity = (state.cells[cellIndex].MliqiniRBuf) / (rightLiquidDensity * rightFlowArea);

    selectAndApplyBufferedInletFlowRegime(
        state, cellIndex, c0, ud, superficialGasVelocity, superficialLiquidVelocity, rightSuperficialLiquidVelocity, rightGasDensity, rightLiquidDensity, gasDensity, liquidDensity,
        flowArea);
    state.cells[1].term1L = state.cells[cellIndex].term1;
    state.cells[1].term2L = state.cells[cellIndex].term2;
}
void prepareNonDimensionalHeatDiffusion(const ThermalState &state, int cellIndex) {
    double diameter = state.cells[cellIndex].duto.a;
    double flowArea = 0.25 * M_PI * diameter * diameter;
    double meanVoidFraction = state.cells[cellIndex].alf;
    double betmed = state.cells[cellIndex].bet;
    double meanSuperficialGasVelocity;
    double meanSuperficialLiquidVelocity;
    if (cellIndex > 0 && (state.cells[cellIndex - 1].acsr.tipo != 5 ||
                  state.cells[cellIndex - 1].acsr.chk.AreaGarg > (1e-3 + state.input.master1.razareaativ) * state.cells[cellIndex - 1].duto.area)) {
        if (state.cells[cellIndex].alf > (*state.globals).localtiny)
            meanSuperficialGasVelocity = state.cells[cellIndex].QG / flowArea;
        else {
            meanSuperficialGasVelocity = 0.;
        }
        if (state.cells[cellIndex].alf < 1. - (*state.globals).localtiny)
            meanSuperficialLiquidVelocity = state.cells[cellIndex].QL / flowArea;
        else {
            meanSuperficialLiquidVelocity = 0.;
        }
    } else {
        if (state.cells[cellIndex].alf > (*state.globals).localtiny)
            meanSuperficialGasVelocity = state.cells[cellIndex + 1].QG / flowArea;
        else {
            meanSuperficialGasVelocity = 0.;
        }

        if (state.cells[cellIndex].alf < 1. - (*state.globals).localtiny)
            meanSuperficialLiquidVelocity = state.cells[cellIndex + 1].QL / flowArea;
        else {
            meanSuperficialLiquidVelocity = 0.;
        }
    }
    double rp = state.cells[cellIndex].rpC;
    double rc = state.cells[cellIndex].rcC;
    double liquidDensity = (1. - betmed) * rp + betmed * rc;
    double gasDensity = state.cells[cellIndex].rgC;
    double liquidSpecificHeat = (1. - betmed) * state.cells[cellIndex].flui.CalorLiq(state.cells[cellIndex].presini, state.cells[cellIndex].temp) + betmed * state.cells[cellIndex].fluicol.CalorLiq(state.cells[cellIndex].presini, state.cells[cellIndex].temp);
    double gasSpecificHeat = state.cells[cellIndex].flui.CalorGas(state.cells[cellIndex].presini, state.cells[cellIndex].temp);

    state.cells[cellIndex].calor.Tint = state.cells[cellIndex].temp;
    state.cells[cellIndex].calor.dtL = state.cells[cellIndex].temp - state.cells[cellIndex - 1].tempini;
    state.cells[cellIndex].calor.Vint = meanSuperficialGasVelocity + meanSuperficialLiquidVelocity;
    state.cells[cellIndex].calor.dt = state.cells[cellIndex].dt;
    double liquidConductivity = (1. - betmed) * state.cells[cellIndex].flui.CondLiq(state.cells[cellIndex].presini, state.cells[cellIndex].temp) + betmed * state.cells[cellIndex].fluicol.CondLiq(state.cells[cellIndex].presini, state.cells[cellIndex].temp); //(1. - betmed) * celula[cellIndex].flui.CondLiq(celula[cellIndex].pres, celula[cellIndex].temp)
    state.cells[cellIndex].calor.kint = liquidConductivity * (1 - meanVoidFraction) + state.cells[cellIndex].flui.CondGas(state.cells[cellIndex].presini, state.cells[cellIndex].temp) * meanVoidFraction;                                                 // liquidConductivity * (1 - meanVoidFraction) + celula[cellIndex].flui.CondGas(celula[cellIndex].pres, celula[cellIndex].temp) * meanVoidFraction;
    state.cells[cellIndex].calor.cpint = liquidSpecificHeat * (1 - meanVoidFraction) + gasSpecificHeat * meanVoidFraction;
    state.cells[cellIndex].calor.rhoint = liquidDensity * (1 - meanVoidFraction) + gasDensity * meanVoidFraction;
    //(1. - betmed) * celula[cellIndex].flui.ViscOleo(celula[cellIndex].pres, celula[cellIndex].temp)
    double liquidViscosity = (1. - betmed) * state.cells[cellIndex].mipC + betmed * state.cells[cellIndex].micC;
    // liquidViscosity * (1 - meanVoidFraction) * 1.e-3
    state.cells[cellIndex].calor.viscint = liquidViscosity * (1 - meanVoidFraction) * 1.e-3 + state.cells[cellIndex].migC * meanVoidFraction * 1.e-3;
    double temperaturePerturbation = state.cells[cellIndex].temp * 0.01;
    if (fabs(state.cells[cellIndex].temp) < 1e-15)
        temperaturePerturbation = 0.1;
    double liquidDensityChange = (1. - betmed) * state.cells[cellIndex].flui.MasEspLiq(state.cells[cellIndex].presini, state.cells[cellIndex].temp + temperaturePerturbation) +
                    betmed * state.cells[cellIndex].fluicol.MasEspFlu(state.cells[cellIndex].presini, state.cells[cellIndex].temp + temperaturePerturbation) - liquidDensity; //(1. - betmed) * celula[cellIndex].flui.MasEspLiq(celula[cellIndex].pres, celula[cellIndex].temp+temperaturePerturbation) +
    double gasDensityChange = state.cells[cellIndex].flui.MasEspGas(state.cells[cellIndex].presini, state.cells[cellIndex].temp + temperaturePerturbation) - gasDensity;             // celula[cellIndex].flui.MasEspGas(celula[cellIndex].pres, celula[cellIndex].temp+temperaturePerturbation)-gasDensity;
    state.cells[cellIndex].calor.betint = -(1 / state.cells[cellIndex].calor.rhoint) * (liquidDensityChange * (1 - meanVoidFraction) + gasDensityChange * meanVoidFraction) / (temperaturePerturbation);
}

void advanceTransientEnergy(const ThermalState &state, int cycle, int maximumCycle) {
    if (((*state.globals).chaverede == 0 || state.networkEndpoint == 1 || (*state.globals).chaveRedeParalela == 1) && state.input.chkv == 0) {
        if (state.cells[state.lastCell - 1].MliqiniR < 0) {
            state.cells[state.lastCell - 1].MR = state.cells[state.lastCell - 1].MR - state.cells[state.lastCell - 1].MliqiniR;
            state.cells[state.lastCell - 1].MliqiniR = 0;
            state.cells[state.lastCell - 1].term1R = 0;
            state.cells[state.lastCell - 1].term2R = 0;
            state.cells[state.lastCell].MC = state.cells[state.lastCell].MC - state.cells[state.lastCell].Mliqini;
            state.cells[state.lastCell].Mliqini = 0;
            state.cells[state.lastCell].term1 = 0;
            state.cells[state.lastCell].term2 = 0;
        }
        if (state.cells[state.lastCell - 1].QLR < 0 && (state.surfaceChokeMassCondition == 0 || state.surfaceChokeOpen == 1)) {
            state.cells[state.lastCell - 1].QLR = 0;
            state.cells[state.lastCell].QL = 0;
        }
    } else if ((state.input.chkv == 1 && state.surfaceChokeMassCondition == 0) || (state.input.chkv == 1 && state.surfaceChokeMassCondition == 1)) {
        if (state.cells[state.lastCell - 1].MliqiniR < 0) {
            state.cells[state.lastCell - 1].MR = 0.;
            state.cells[state.lastCell - 1].MliqiniR = 0;
            state.cells[state.lastCell - 1].term1R = 0;
            state.cells[state.lastCell - 1].term2R = 0;
            state.cells[state.lastCell].MC = 0.;
            state.cells[state.lastCell].Mliqini = 0;
            state.cells[state.lastCell].term1 = 0;
            state.cells[state.lastCell].term2 = 0;
        }
        if (state.cells[state.lastCell - 1].QLR < 0 && (state.surfaceChokeMassCondition == 0 || state.surfaceChokeOpen == 1)) {
            state.cells[state.lastCell - 1].QLR = 0;
            state.cells[state.lastCell].QL = 0;
        }
    }

    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // atencao!!!!!!!!!!!!!!!!!!
    // existe uma questao que parece mal resolvida na resolucao desta marcha, nao foi feito nenhum teste para
    // o caso em que a velocidade de transporte da temperatura é <0 neste caso, a temperatura na celula de indice
    // não deveria entrar no metodo calctemp, ja que não e mais o caso de ser uma celula com condicao de
    // contorno para temperatura??????????????????????????????????????????????????????????????????????/
    //"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    state.cells[0].tempini = state.cells[0].temp;
    if (state.input.ConContEntrada > 0)
        state.cells[0].temp = state.inletTemperature;
    else
        state.cells[0].temp = state.defaultInletTemperature;
    state.cells[0].dTdt = (state.cells[0].temp - state.cells[0].tempini) / state.timeStep;
    state.cells[0].dTdtL = state.cells[0].dTdt;
    state.cells[1].tempLini = state.cells[1].tempL;
    state.cells[1].tempL = state.cells[0].temp;
    for (int cellIndex = 1; cellIndex <= state.lastCell; cellIndex++) {
        state.cells[cellIndex].tempini = state.cells[cellIndex].temp;
    }
    if (state.input.modoDifus3D == 0) {
        if (state.poisson2DCellCount > 1) {
#pragma omp parallel for num_threads((*state.globals).ntrd)
            for (int poisson2DIndex = 0; poisson2DIndex < state.poisson2DCellCount; poisson2DIndex++) {
                int cellIndex = state.poisson2DCellIndices[poisson2DIndex];
                if (cellIndex <= state.lastCell) {
                    computeTemperature(state, cellIndex, state.cells[cellIndex].tempini);
                } else {
                    if ((*state.globals).chaverede == 0 || state.networkEndpoint == 1 || (*state.globals).chaveRedeParalela == 1)
                        state.cells[cellIndex].temp = state.cells[cellIndex].calor.Textern1;
                    else
                        state.cells[cellIndex].temp = state.gasSurfaceTemperature;
                }
                state.cells[cellIndex].dTdt = (state.cells[cellIndex].temp - state.cells[cellIndex].tempini) / state.timeStep;
                state.cells[cellIndex].dTdtIni = state.cells[cellIndex].dTdt;
            }
        }
#pragma omp parallel for num_threads((*state.globals).ntrd)
        for (int cellIndex = 0; cellIndex <= state.lastCell; cellIndex++) {
            if (cellIndex == 217) {
                int debugStop;
                debugStop = 0;
            }
            if (state.cells[cellIndex].calor.difus2D == 0) {
                if (cellIndex <= state.lastCell) {
                    computeTemperature(state, cellIndex, state.cells[cellIndex].tempini);
                } else {
                    if ((*state.globals).chaverede == 0 || state.networkEndpoint == 1 || (*state.globals).chaveRedeParalela == 1)
                        state.cells[cellIndex].temp = state.cells[cellIndex].calor.Textern1;
                    else
                        state.cells[cellIndex].temp = state.gasSurfaceTemperature;
                }
                state.cells[cellIndex].dTdt = (state.cells[cellIndex].temp - state.cells[cellIndex].tempini) / state.timeStep;
                state.cells[cellIndex].dTdtIni = state.cells[cellIndex].dTdt;
            }
        }
    } else if (state.input.modoDifus3D == 1) {
#pragma omp parallel for num_threads((*state.globals).ntrd)
        for (int cellIndex = 0; cellIndex <= state.lastCell; cellIndex++) {
            int coupled = -1;
            int coupledCellIndex;
            for (int couplingIndex = 0; couplingIndex < state.input.nacop; couplingIndex++) {
                coupledCellIndex = state.input.celAcop[couplingIndex].indCel;
                if (cellIndex == coupledCellIndex) {
                    coupled = couplingIndex;
                    break;
                }
            }
            if (coupled != -1) {
                prepareNonDimensionalHeatDiffusion(state, cellIndex);
                int nextCouplingIndex = state.coupledCellIndices[coupled];
                state.poissonSolver.dados.tInt[nextCouplingIndex] = state.cells[coupledCellIndex].temp;
                double cellInterfaceCoefficient = state.cells[coupledCellIndex].calor.hInt();
                state.poissonSolver.dados.hI[nextCouplingIndex] = cellInterfaceCoefficient;
            }
        }
        if (cycle < maximumCycle || maximumCycle == 0 || state.minimumCycleTimeStep != state.timeStep) {
            if (cycle == maximumCycle && maximumCycle > 0)
                state.poissonSolver.FeiticoDoTempo();
            state.poissonSolver.transientePoisson(state.timeStep);
        }
        if (state.poisson2DCellCount > 1) {
#pragma omp parallel for num_threads((*state.globals).ntrd)
            for (int poisson2DIndex = 0; poisson2DIndex < state.poisson2DCellCount; poisson2DIndex++) {
                int cellIndex = state.poisson2DCellIndices[poisson2DIndex];
                if (cellIndex <= state.lastCell) {
                    computeTemperature(state, cellIndex, state.cells[cellIndex].tempini);
                } else {
                    if ((*state.globals).chaverede == 0 || state.networkEndpoint == 1 || (*state.globals).chaveRedeParalela == 1)
                        state.cells[cellIndex].temp = state.cells[cellIndex].calor.Textern1;
                    else
                        state.cells[cellIndex].temp = state.gasSurfaceTemperature;
                }
                state.cells[cellIndex].dTdt = (state.cells[cellIndex].temp - state.cells[cellIndex].tempini) / state.timeStep;
                state.cells[cellIndex].dTdtIni = state.cells[cellIndex].dTdt;
            }
        }
#pragma omp parallel for num_threads((*state.globals).ntrd)
        for (int cellIndex = 0; cellIndex <= state.lastCell; cellIndex++) {
            if (state.cells[cellIndex].calor.difus2D == 0) {
                if (cellIndex <= state.lastCell) {
                    computeTemperature(state, cellIndex, state.cells[cellIndex].tempini);
                } else {
                    if ((*state.globals).chaverede == 0 || state.networkEndpoint == 1 || (*state.globals).chaveRedeParalela == 1)
                        state.cells[cellIndex].temp = state.cells[cellIndex].calor.Textern1;
                    else
                        state.cells[cellIndex].temp = state.gasSurfaceTemperature;
                }
                state.cells[cellIndex].dTdt = (state.cells[cellIndex].temp - state.cells[cellIndex].tempini) / state.timeStep;
                state.cells[cellIndex].dTdtIni = state.cells[cellIndex].dTdt;
            }
        }
    }
    for (int cellIndex = 1; cellIndex <= state.lastCell; cellIndex++) {
        state.cells[cellIndex].dTdtL = state.cells[cellIndex - 1].dTdt;
        if (cellIndex < state.lastCell) {
            state.cells[cellIndex + 1].tempLini = state.cells[cellIndex + 1].tempL;
            state.cells[cellIndex + 1].tempL = state.cells[cellIndex].temp;
        }
        state.cells[cellIndex - 1].tempRini = state.cells[cellIndex - 1].tempR;
        state.cells[cellIndex - 1].tempR = state.cells[cellIndex].temp;
    }
    if (cycle < maximumCycle) {
        for (int k = 0; k <= state.lastCell; k++) {
            state.cells[k].FeiticoDoTempo();
        }
    } else if (state.input.modoDifus3D == 1)
        state.poissonSolver.renova();
    if (state.completeModel == 0) {
        if (cycle < maximumCycle) {
            state.evolutionUpdater.solvePressureVelocityCoupling(cycle);
            state.evolutionUpdater.renew();
        }
    }
}

void advanceSteadyTemperature(const ThermalState &state, int cellIndex, int rungeKuttaStage) {
    double cellLength = 0.5 * (state.cells[cellIndex].dx + state.cells[cellIndex - 1].dx);
    double meanCellLength = 0.5 * (state.cells[cellIndex].dx + state.cells[cellIndex - 1].dx);
    double meanTemperatureGradientCorrection = (state.cells[cellIndex].dx * state.cells[cellIndex].dTdLCor + state.cells[cellIndex - 1].dx * state.cells[cellIndex - 1].dTdLCor) / (state.cells[cellIndex].dx + state.cells[cellIndex - 1].dx);
    double diameter = state.cells[cellIndex - 1].duto.a;
    double flowArea = 0.25 * M_PI * diameter * diameter;
    double meanVoidFraction;
    double betmed;
    double meanPressure;
    double meanTemperature;
    if (rungeKuttaStage == 0) {
        meanVoidFraction = state.cells[cellIndex - 1].alf;
        betmed = state.cells[cellIndex - 1].bet;
        meanPressure = state.cells[cellIndex - 1].pres;
        meanTemperature = state.cells[cellIndex - 1].temp;
    } else {
        meanVoidFraction = state.cells[cellIndex].alf;
        betmed = state.cells[cellIndex].bet;
        meanPressure = state.cells[cellIndex].pres;
        meanTemperature = state.cells[cellIndex].temp;
    }
    double meanSuperficialGasVelocity;
    meanSuperficialGasVelocity = state.cells[cellIndex].QG / flowArea; // velocidade superficial de gas
    double meanSuperficialLiquidVelocity;
    meanSuperficialLiquidVelocity = state.cells[cellIndex].QL / flowArea; // velocidade superficial de liquido
    double mixtureFluxSign = 1.;
    if (fabs(meanSuperficialGasVelocity + meanSuperficialLiquidVelocity) > 0.05 && state.thermalSourceDisabled == 0) { // calculo termico e feito para velocidades de mistura superiores a 0,1 m/s,
        // para velocidades inferiores se admite que a temperatura do fluidoÃƒÂ© igual ÃƒÂ  temperatura ambiente
        mixtureFluxSign = (meanSuperficialGasVelocity + meanSuperficialLiquidVelocity) / fabs(meanSuperficialGasVelocity + meanSuperficialLiquidVelocity);
        double interfaceMeanPressure = state.cells[cellIndex].presaux + 0 * state.cells[cellIndex - 1].dpB / 98066.5;
        double interfaceMeanTemperature;
        if (state.steadyIteration != 0 && state.input.AceleraConvergPerm == 0) // temperatura na interface esquerda da celula
            // caso em que se considera que ja foi
            // feita uma iteracao e ja se tem a temperatura na celula cellIndex vinda da iteracao anterior
            // isto pode dificultar a convergencia, quando se deseja a aceleracao da convergencia
            // admite-se que a temperatura na fronteira esquerda ÃƒÂ© a temperatura da celula esquerda
            interfaceMeanTemperature = (state.cells[cellIndex].dx * state.cells[cellIndex].temp + state.cells[cellIndex].dxL * state.cells[cellIndex].tempL) / (state.cells[cellIndex].dx + state.cells[cellIndex].dxL);
        else
            interfaceMeanTemperature = state.cells[cellIndex - 1].temp;
        double rp = state.cells[cellIndex - 1].flui.MasEspLiq(meanPressure, meanTemperature);    // celula[cellIndex].rpCi;
        double rc = state.cells[cellIndex - 1].fluicol.MasEspFlu(meanPressure, meanTemperature); // celula[cellIndex].rcCi;
        double liquidDensity = (1. - betmed) * rp + betmed * rc;
        double gasDensity = state.cells[cellIndex - 1].flui.MasEspGas(meanPressure, meanTemperature); // celula[cellIndex].rgCi;
        double liquidSpecificHeat = (1. - betmed) * state.cells[cellIndex - 1].flui.CalorLiq(interfaceMeanPressure, interfaceMeanTemperature) + betmed * state.cells[cellIndex - 1].fluicol.CalorLiq(interfaceMeanPressure, interfaceMeanTemperature);
        double gasSpecificHeat = state.cells[cellIndex - 1].flui.CalorGas(interfaceMeanPressure, interfaceMeanTemperature);
        // liquidJouleThomson=Joule Thomson do liquido X cp
        // gasJouleThomson=Joule Thomson do gas X cp
        /////??????????????????????????????????????????????????????????????????????????????????????????
        double liquidJouleThomson = (1. - betmed) * state.cells[cellIndex - 1].flui.JTL(interfaceMeanPressure, interfaceMeanTemperature) - betmed / rc;
        /////??????????????????????????????????????????????????????????????????????????????????????????
        if (state.input.pocinjec > 0 && state.input.condpocinj.tipoFlui == 2) {
            liquidJouleThomson = -(1 + (interfaceMeanTemperature + 273.14) * state.cells[cellIndex - 1].fluicol.DrhoDtFlu(interfaceMeanPressure, interfaceMeanTemperature) / rc) / rc;
        }
        double gasJouleThomson = state.cells[cellIndex - 1].flui.JTG(interfaceMeanPressure, interfaceMeanTemperature);
        // energia potencial:
        double hydrostaticPower = (liquidDensity * meanSuperficialLiquidVelocity + gasDensity * meanSuperficialGasVelocity) * flowArea * 9.82 * sin(state.cells[cellIndex - 1].duto.teta);

        if (cellIndex > 120) {
            int debugStop;
            debugStop = 0;
        }

        // definicao dos parametros internos na tubulacao para se obter a troca termica om o meio ambiente
        state.cells[cellIndex - 1].calor.Tint = interfaceMeanTemperature;
        state.cells[cellIndex - 1].calor.Vint = fabs(meanSuperficialGasVelocity + meanSuperficialLiquidVelocity);
        double liquidConductivity = (1. - betmed) * state.cells[cellIndex - 1].flui.CondLiq(interfaceMeanPressure, interfaceMeanTemperature) + betmed * state.cells[cellIndex - 1].fluicol.CondLiq(interfaceMeanPressure, interfaceMeanTemperature);
        state.cells[cellIndex - 1].calor.kint = liquidConductivity * (1 - meanVoidFraction) + state.cells[cellIndex - 1].flui.CondGas(interfaceMeanPressure, interfaceMeanTemperature) * meanVoidFraction;
        state.cells[cellIndex - 1].calor.cpint = liquidSpecificHeat * (1 - meanVoidFraction) + gasSpecificHeat * meanVoidFraction;
        state.cells[cellIndex - 1].calor.rhoint = liquidDensity * (1 - meanVoidFraction) + gasDensity * meanVoidFraction;
        double liquidViscosity = (1. - betmed) * state.cells[cellIndex - 1].flui.ViscOleo(interfaceMeanPressure, interfaceMeanTemperature) + betmed * state.cells[cellIndex - 1].fluicol.VisFlu(interfaceMeanPressure, interfaceMeanTemperature);
        state.cells[cellIndex - 1].calor.viscint = liquidViscosity * (1 - meanVoidFraction) * 1.e-3 + state.cells[cellIndex - 1].flui.ViscGas(interfaceMeanPressure, interfaceMeanTemperature) * meanVoidFraction * 1.e-3;

        double relaxedHeatFlux = 0; // variavel nao utilizada
        if (state.steadyIteration > 0)
            relaxedHeatFlux = state.cells[cellIndex - 1].fluxcalmed;
        double annulusResistance = 0.;
        double heatFlux;
        double gasHeatFlux;
        // interfacialWorkTerm se existe acoplamento com o anular:
        if (state.input.lingas == 1 && (cellIndex - 1 <= state.annulusTubingStart && cellIndex - 1 >= state.annulusTubingEnd)) {
            int j = state.annulusTubingStart + state.tubingAnnulusStart - (cellIndex - 1);
            // caso tenha acoplamento:

            if (state.steadyIteration == 0 && (cellIndex - 1) < state.annulusTubingStart) {
                // na primeira iteracao, considera-se a resistencia do revestimento + cimento +
                // formacao, nas outras iteracoes
                // bastara determinar a troca termica entre a coluna e o gas do anular
                // no modelo acoplado, a coluna nao tem a definicao da parede revestimento+cimento+formacao
                // para se obter esta resistencia, precisa-se recorrer ao modelo de troca termica
                // do anular, o que e feito aqui:
                // obs:isto Ã© feito ate uma celula antes de se chegar na master (cellIndex-1)<ColunaAnulaIni,
                // se esta fazendo igual ao que se faz no simulador involuta, para melhorar
                // a estimativa da temperatura na ANM. Na celula da anm, so se considera a troca termica com o gas
                // desde a primeira iteracao
                state.gasCells[j].calor.Tint = interfaceMeanTemperature;
                state.gasCells[j].calor.Vint = 100;
                state.gasCells[j].calor.kint = state.cells[cellIndex - 1].flui.CondGas(interfaceMeanPressure, interfaceMeanTemperature);
                state.gasCells[j].calor.cpint = gasSpecificHeat;
                state.gasCells[j].calor.rhoint = gasDensity;
                state.gasCells[j].calor.viscint = state.cells[cellIndex - 1].flui.ViscGas(interfaceMeanPressure, interfaceMeanTemperature) * 1.e-3;
                state.gasCells[j].fluxcal = state.gasCells[j].calor.transperm(); // troca termica no anular
                annulusResistance = state.gasCells[j].calor.resGlob;                // resistencia das paredes
                // observe que nÃ£o se esta de fato interessado no fluxo de calor, mas apenas em obter a resistencia
                // termica do conjunto de paredes a partir do revestimento em direcao aa formacao
                state.cells[cellIndex - 1].calor.Vextern1 = 100.;
                state.cells[cellIndex - 1].calor.kextern1 = state.gasCells[j].calor.kint;
                state.cells[cellIndex - 1].calor.cpextern1 = state.gasCells[j].calor.cpint;
                state.cells[cellIndex - 1].calor.rhoextern1 = state.gasCells[j].calor.rhoint;
                state.cells[cellIndex - 1].calor.viscextern1 = state.gasCells[j].calor.viscint;
            } else { // apos a primeira iteracao, considera-se apenas a troca termica entre a coluna e o gas do anular
                // passando pela parede da coluna, claro
                annulusResistance = 0.;
                state.cells[cellIndex - 1].calor.Vextern1 = state.gasCells[j].VGasR / state.gasCells[j].u1L;
                state.cells[cellIndex - 1].calor.kextern1 = state.gasCells[j].calor.kint;
                state.cells[cellIndex - 1].calor.cpextern1 = state.gasCells[j].calor.cpint;
                state.cells[cellIndex - 1].calor.rhoextern1 = state.gasCells[j].calor.rhoint;
                state.cells[cellIndex - 1].calor.viscextern1 = state.gasCells[j].calor.viscint;
            }
            if (state.steadyIteration == 0)
                state.cells[cellIndex - 1].calor.Textern1 = state.gasCells[j].calor.Textern1; // na primeira iteracao, como se
            // usa toda a resistencia termica do poco, a temperatura externa utilizada Ã© a geotermica
            else
                state.cells[cellIndex - 1].calor.Textern1 = state.gasCells[j].temp; // nas iteracoes seguintes, a temperatura ambiente
            // e a temperatura do gas
        }
        state.cells[cellIndex - 1].fluxcalmed = 0;
        if (state.productionNetworkHeatCoupled == 1 && (cellIndex - 1) >= state.primaryNetworkSectionEnd && (cellIndex - 1) <= state.primaryNetworkSectionStart) {
            heatFlux = mixtureFluxSign * state.cells[cellIndex - 1].calor.transperm(state.cells[cellIndex - 1].resAcopRedeP);
        } else
            heatFlux = mixtureFluxSign * state.cells[cellIndex - 1].calor.transperm(annulusResistance);
        state.cells[cellIndex - 1].fluxcalmed = heatFlux; // fluxo de calor na coluna

        double temperatureSpatialCoefficient = (liquidDensity * meanSuperficialLiquidVelocity * liquidSpecificHeat + gasDensity * meanSuperficialGasVelocity * gasSpecificHeat) * flowArea; // termo que multiplica
        // a derivada Dt/Dx
        double cappedMeanSuperficialLiquidVelocity = meanSuperficialLiquidVelocity;
        if ((*state.globals).blackOilTemp == 1 && fabs(meanSuperficialLiquidVelocity) > 5)
            cappedMeanSuperficialLiquidVelocity = 5 * meanSuperficialLiquidVelocity / meanSuperficialLiquidVelocity;
        double cappedMeanSuperficialGasVelocity = meanSuperficialGasVelocity;
        if ((*state.globals).blackOilTemp == 1 && fabs(meanSuperficialGasVelocity) > 5)
            cappedMeanSuperficialGasVelocity = 5 * meanSuperficialGasVelocity / fabs(meanSuperficialGasVelocity);
        double pressureSpatialCoefficient = 1 * (liquidDensity * cappedMeanSuperficialLiquidVelocity * liquidJouleThomson + gasDensity * cappedMeanSuperficialGasVelocity * gasJouleThomson) * flowArea; // termo que multiplica
        // a derivada Dp/Dx
        double pressureGradient;
        if ((state.cells[cellIndex - 1].acsr.tipo != 4 || state.cells[cellIndex - 1].acsr.bcs.freq < 1) && state.cells[cellIndex - 1].acsr.tipo != 7)
            // caso nao tenha BCS ou incremento de pressao  utiliza-se a pressao na fronteira esquerda
            //  e a pressao no centro de celula para o calculo de Dp/Dx
            pressureGradient = 2. * (state.cells[cellIndex].presaux - state.cells[cellIndex - 1].pres) * 98066.5 / state.cells[cellIndex - 1].dx;
        else {
            // caso tenha BCS ou incremento de pressao  utiliza-se a pressao da celulaa esquerda
            //  e a pressao no centro de celula para o calculoi de Dp/Dx
            pressureGradient = 2. * (state.cells[cellIndex].presaux - state.cells[cellIndex - 1].pres) * 98066.5 / state.cells[cellIndex - 1].dx;
        }
        state.cells[cellIndex].VTemper = meanSuperficialLiquidVelocity; // esta velocidade so e util no caso transiente, Ã© armazenada aqui
        // apenas para se ter um valor quando a simulacao transiente se iniciar
        double temperatureGradient = (-state.cells[cellIndex - 1].temp) / meanCellLength;

        double kineticTerm = 0;
        double upstreamMeanGasVelocity = 0;
        double upstreamMeanLiquidVelocity = 0;
        double meanGasVelocity = 0;
        double meanLiquidVelocity = 0;
        // termo de energia cinetica:
        if (state.cells[cellIndex].acsr.tipo == 0 && state.cells[cellIndex - 1].acsr.tipo == 0 && cellIndex > 2) {
            double kineticCellLength = state.cells[cellIndex - 1].dx;
            double upstreamDiameter = state.cells[cellIndex - 1].duto.a;
            double upstreamFlowArea = 0.25 * M_PI * upstreamDiameter * upstreamDiameter;

            if (state.cells[cellIndex - 1].alf > 1e-3)
                meanGasVelocity = meanSuperficialGasVelocity / state.cells[cellIndex - 1].alf;
            if (state.cells[cellIndex - 1].alf < (1. - 1e-3))
                meanLiquidVelocity = meanSuperficialLiquidVelocity / (1. - state.cells[cellIndex - 1].alf);

            if (state.cells[cellIndex - 2].alf > 1e-3) {
                upstreamMeanGasVelocity = state.cells[cellIndex - 1].QG / (upstreamFlowArea);
                upstreamMeanGasVelocity /= state.cells[cellIndex - 2].alf;
            }
            if (state.cells[cellIndex - 2].alf < (1. - 1e-3)) {
                upstreamMeanLiquidVelocity = state.cells[cellIndex - 1].QL / (upstreamFlowArea);
                upstreamMeanLiquidVelocity /= (1. - state.cells[cellIndex - 2].alf);
            }

            if (state.input.nCompTotalUnidadesP / kineticCellLength < 1e6)
                kineticTerm = (state.cells[cellIndex].MC - state.cells[cellIndex].Mliqini) * meanGasVelocity * (meanGasVelocity - upstreamMeanGasVelocity) / kineticCellLength + state.cells[cellIndex].Mliqini * meanLiquidVelocity * (meanLiquidVelocity - upstreamMeanLiquidVelocity) / kineticCellLength;
            else
                kineticTerm = 0;
        }

        double gasMassSourceTerm = 0.;
        double liquidMassSourceTerm = 0.;
        double fontemassC = 0.;
        double sourceTemperature = state.cells[cellIndex - 1].temp;
        double sourceGasSpecificHeat;
        double sourceSpecificHeatRatio = 0.;
        double sourceLiquidSpecificHeat;

        // calculo da energia adicionada no sistema devido a fontes de massa
        if (state.cells[cellIndex - 1].acsr.tipo == 1) { // caso fonte de gas
            sourceTemperature = state.cells[cellIndex - 1].acsr.injg.temp;
            sourceGasSpecificHeat = state.cells[cellIndex - 1].acsr.injg.FluidoPro.CalorGas(meanPressure, sourceTemperature);
            sourceSpecificHeatRatio = state.cells[cellIndex - 1].acsr.injg.FluidoPro.ConstAdG(meanPressure, sourceTemperature);
            sourceLiquidSpecificHeat = 0.;
        } else if (state.cells[cellIndex - 1].acsr.tipo == 2) { // caso fonte de liquido
            sourceTemperature = state.cells[cellIndex - 1].acsr.injl.temp;
            sourceGasSpecificHeat = 0.;
            sourceSpecificHeatRatio = 1.;
            sourceLiquidSpecificHeat = (1. - state.cells[cellIndex - 1].acsr.injl.bet) * state.cells[cellIndex - 1].acsr.injl.FluidoPro.CalorLiq(meanPressure, meanTemperature) + state.cells[cellIndex - 1].acsr.injl.bet * state.cells[cellIndex - 1].acsr.injl.fluidocol.CalorLiq(meanPressure, meanTemperature);
        } else if (state.cells[cellIndex - 1].acsr.tipo == 3) { // caso IPR
            sourceTemperature = state.cells[cellIndex - 1].acsr.ipr.Tres;
            sourceGasSpecificHeat = state.cells[cellIndex - 1].acsr.ipr.FluidoPro.CalorGas(meanPressure, meanTemperature);
            sourceSpecificHeatRatio = state.cells[cellIndex - 1].acsr.ipr.FluidoPro.ConstAdG(meanPressure, meanTemperature);
            sourceLiquidSpecificHeat = state.cells[cellIndex - 1].acsr.ipr.FluidoPro.CalorLiq(meanPressure, meanTemperature);
        } else if (state.cells[cellIndex - 1].acsr.tipo == 9 && state.cells[cellIndex - 1].acsr.fontechk.abertura > 1e-6 &&
                   (state.cells[cellIndex - 1].fontemassCR + state.cells[cellIndex - 1].fontemassGR + state.cells[cellIndex - 1].fontemassLR) > 1e-9) {
            // caso vazamento
            sourceTemperature = state.cells[cellIndex - 1].acsr.fontechk.tamb;
            sourceGasSpecificHeat = state.cells[cellIndex - 1].acsr.fontechk.fluidoPamb.CalorGas(meanPressure, meanTemperature);
            sourceSpecificHeatRatio = state.cells[cellIndex - 1].acsr.fontechk.fluidoPamb.ConstAdG(meanPressure, meanTemperature);
            sourceLiquidSpecificHeat = (1. - state.cells[cellIndex - 1].acsr.fontechk.betISamb) *
                       state.cells[cellIndex - 1].acsr.fontechk.fluidoPamb.CalorLiq(meanPressure, meanTemperature) +
                   state.cells[cellIndex - 1].acsr.fontechk.betISamb * state.cells[cellIndex - 1].acsr.fontechk.fluidocol.CalorLiq(meanPressure, meanTemperature);
        } else if (state.cells[cellIndex - 1].acsr.tipo == 15) { // caso IPR
            sourceTemperature = state.cells[cellIndex - 1].acsr.radialPoro.tRes;
            sourceGasSpecificHeat = state.cells[cellIndex - 1].acsr.radialPoro.flup.CalorGas(meanPressure, meanTemperature);
            sourceSpecificHeatRatio = state.cells[cellIndex - 1].acsr.radialPoro.flup.ConstAdG(meanPressure, meanTemperature);
            sourceLiquidSpecificHeat = state.cells[cellIndex - 1].acsr.radialPoro.flup.CalorLiq(meanPressure, meanTemperature);
        } else if (state.cells[cellIndex - 1].acsr.tipo == 16) { // caso IPR
            sourceTemperature = state.cells[cellIndex - 1].acsr.poroso2D.dados.tRes;
            sourceGasSpecificHeat = state.cells[cellIndex - 1].acsr.poroso2D.dados.flup.CalorGas(meanPressure, meanTemperature);
            sourceSpecificHeatRatio = state.cells[cellIndex - 1].acsr.poroso2D.dados.flup.ConstAdG(meanPressure, meanTemperature);
            sourceLiquidSpecificHeat = state.cells[cellIndex - 1].acsr.poroso2D.dados.flup.CalorLiq(meanPressure, meanTemperature);
        } else if (state.cells[cellIndex].acsrL != 0) {

            sourceGasSpecificHeat = 0.;
            sourceSpecificHeatRatio = 1.;
            sourceLiquidSpecificHeat = 0.;
        } else {
            sourceGasSpecificHeat = 0.;
            sourceSpecificHeatRatio = 1.;
            sourceLiquidSpecificHeat = 0.;
        }

        liquidMassSourceTerm = 0;
        if (state.cells[cellIndex - 1].fontemassLR > 0.)
            liquidMassSourceTerm = state.cells[cellIndex - 1].fontemassLR / cellLength;
        if (state.cells[cellIndex - 1].fontemassCR > 0.)
            liquidMassSourceTerm += state.cells[cellIndex - 1].fontemassCR / cellLength;
        liquidMassSourceTerm *= (sourceLiquidSpecificHeat) * (sourceTemperature - state.cells[cellIndex - 1].temp);

        gasMassSourceTerm = state.cells[cellIndex - 1].fontemassGR / cellLength;
        if (gasMassSourceTerm > 0.)
            gasMassSourceTerm *= (sourceGasSpecificHeat / sourceSpecificHeatRatio) * (sourceTemperature - state.cells[cellIndex - 1].temp);
        else
            gasMassSourceTerm = 0;

        // efeito do calor latentHeatTerm, quando este for solicitado
        double latentHeatTerm = 0.;
        if (isnan(state.cells[cellIndex - 1].FonteMudaFase))
            state.cells[cellIndex - 1].FonteMudaFase = 0.;
        double phaseChangeMassRate = fabs(state.cells[cellIndex - 1].FonteMudaFase);
        double phaseChangeSign = 1.;
        if (phaseChangeMassRate > 1e-25)
            phaseChangeSign = state.cells[cellIndex - 1].FonteMudaFase / phaseChangeMassRate;
        if (state.input.limTransMass < phaseChangeMassRate)
            phaseChangeMassRate = phaseChangeSign * state.input.limTransMass;
        else
            phaseChangeMassRate *= phaseChangeSign;
        if (state.cells[cellIndex].flui.dVaporMassFraction < (1 - 1e-15) && state.cells[cellIndex].flui.dVaporMassFraction > (1e-15)) {
            if (state.cells[cellIndex - 1].acsr.tipo == 1 || state.cells[cellIndex - 1].acsr.tipo == 2 || state.cells[cellIndex - 1].acsr.tipo == 3 || state.cells[cellIndex - 1].acsr.tipo == 15 || state.cells[cellIndex - 1].acsr.tipo == 16) {
                state.cells[cellIndex - 1].FonteMudaFase = 0.;
                phaseChangeMassRate = 0.;
            }
            if (state.latentHeatEnabled > 0 && state.steadyIteration != 0 && state.input.flashCompleto == 0) {
                latentHeatTerm = interpolateLatentHeat(state, meanPressure, meanTemperature) * phaseChangeMassRate;
            } else if ((state.input.flashCompleto == 1 || state.input.flashCompleto == 2) && state.latentHeatEnabled > 0 && state.steadyIteration != 0) {
                latentHeatTerm = (state.cells[cellIndex].flui.EntalpGas(meanPressure, meanTemperature) -
                           state.cells[cellIndex].flui.EntalpLiq(meanPressure, meanTemperature)) *
                          phaseChangeMassRate;
            } else
                latentHeatTerm = 0;
        }

        double interfaceVoidFraction;
        double leftInterfaceVoidFraction;
        if (meanSuperficialGasVelocity > 0) {
            interfaceVoidFraction = state.cells[cellIndex - 1].alf;
            leftInterfaceVoidFraction = state.cells[cellIndex - 1].alfL;
        } else {
            interfaceVoidFraction = state.cells[cellIndex - 1].alfR;
            leftInterfaceVoidFraction = state.cells[cellIndex - 1].alf;
        }
        double slipVelocity;
        if (interfaceVoidFraction > (*state.globals).localtiny && interfaceVoidFraction < (1. - (*state.globals).localtiny))
            slipVelocity = meanSuperficialGasVelocity / interfaceVoidFraction - meanSuperficialLiquidVelocity / (1. - interfaceVoidFraction);
        else if (interfaceVoidFraction > (*state.globals).localtiny)
            slipVelocity = meanSuperficialGasVelocity;
        else
            slipVelocity = meanSuperficialLiquidVelocity;
        double interfacialWorkTerm = flowArea * state.cells[cellIndex].pres * 98600 * slipVelocity * (interfaceVoidFraction - leftInterfaceVoidFraction) / state.cells[cellIndex].dx;

        if (fabs(temperatureSpatialCoefficient) > (*state.globals).localtiny) {
            // Energy terms for boundary work, potential and kinetic energy,
            // mass sources, latent heat, and shaft work.
            if (state.input.latente == 0)
                latentHeatTerm = 0.;
            else if (state.input.condlatente == 0 && latentHeatTerm < 0)
                latentHeatTerm = 0.;
            double sourceTemperatureGradient = meanTemperatureGradientCorrection * (pressureSpatialCoefficient * pressureGradient - kineticTerm - (hydrostaticPower) + (liquidMassSourceTerm + gasMassSourceTerm) - latentHeatTerm + (state.cells[cellIndex - 1].potTermo + state.cells[cellIndex - 1].fonteCal) / meanCellLength) / temperatureSpatialCoefficient;
            // portion of energy related to heat exchange
            double heatFluxTemperatureGradient = meanTemperatureGradientCorrection * (heatFlux) / temperatureSpatialCoefficient;

            // Check whether heat transfer is too fast for the explicit temperature update.
            // If so, use additional substeps to prevent thermal instability.
            int subStepCount;
            double subStepLength;
            double stabilityLength = fabs(temperatureSpatialCoefficient) / meanTemperatureGradientCorrection;
            if (meanCellLength / (state.cells[cellIndex - 1].calor.resGlob + annulusResistance) < (stabilityLength + 0. * 1000.)) { // Thermal resistance is not low.
                subStepCount = 1;
                subStepLength = meanCellLength;
            } else {
                // For low thermal resistance, determine the number of temperature update steps.
                subStepCount = (meanCellLength / (state.cells[cellIndex - 1].calor.resGlob + annulusResistance)) / (stabilityLength + 0 * 1000) + 1;
                subStepLength = meanCellLength / subStepCount; // cell divided into npassos
            }
            double subStepTemperature = state.cells[cellIndex - 1].temp;

            subStepTemperature = subStepLength * (-(-state.cells[cellIndex - 1].temp) / subStepLength + sourceTemperatureGradient + heatFluxTemperatureGradient); // first step
            for (int j = 1; j < subStepCount; j++) {                                                  // next steps
                // Keep all energy terms except heat flow, which is recalculated at each step.
                state.cells[cellIndex - 1].calor.Tint = subStepTemperature;
                state.cells[cellIndex - 1].calor.Vint = meanSuperficialGasVelocity + meanSuperficialLiquidVelocity;
                liquidConductivity = (1. - betmed) * state.cells[cellIndex - 1].flui.CondLiq(meanPressure, subStepTemperature) + betmed * state.cells[cellIndex - 1].fluicol.CondLiq(meanPressure, subStepTemperature);
                state.cells[cellIndex - 1].calor.kint = liquidConductivity * (1 - meanVoidFraction) + state.cells[cellIndex - 1].flui.CondGas(meanPressure, subStepTemperature) * meanVoidFraction;
                state.cells[cellIndex - 1].calor.cpint = liquidSpecificHeat * (1 - meanVoidFraction) + gasSpecificHeat * meanVoidFraction;
                state.cells[cellIndex - 1].calor.rhoint = liquidDensity * (1 - meanVoidFraction) + gasDensity * meanVoidFraction;
                liquidViscosity = (1. - betmed) * state.cells[cellIndex - 1].flui.ViscOleo(meanPressure, subStepTemperature) + betmed * state.cells[cellIndex - 1].fluicol.VisFlu(meanPressure, subStepTemperature);
                state.cells[cellIndex - 1].calor.viscint = liquidViscosity * (1 - meanVoidFraction) * 1.e-3 + state.cells[cellIndex - 1].flui.ViscGas(meanPressure, subStepTemperature) * meanVoidFraction * 1.e-3;
                if (state.steadyIteration != 0 && state.input.lingas == 1 && (cellIndex - 1 <= state.annulusTubingStart && cellIndex - 1 >= state.annulusTubingEnd)) {
                    int k = state.annulusTubingStart + state.tubingAnnulusStart - (cellIndex - 1);
                    double externalTemperatureStep = (state.gasCells[k - 1].temp - state.gasCells[k].temp) / subStepCount;
                    state.cells[cellIndex - 1].calor.Textern1 = state.gasCells[k].temp + (j - 1) * externalTemperatureStep;
                }
                heatFlux = mixtureFluxSign * state.cells[cellIndex - 1].calor.transperm(annulusResistance);
                state.cells[cellIndex - 1].fluxcalmed += heatFlux;
                heatFluxTemperatureGradient = meanTemperatureGradientCorrection * (heatFlux) / temperatureSpatialCoefficient;
                subStepTemperature = subStepLength * (-(-subStepTemperature) / subStepLength + sourceTemperatureGradient + heatFluxTemperatureGradient);
            }

            state.cells[cellIndex].temp = subStepTemperature;
            state.cells[cellIndex - 1].fluxcalmed /= subStepCount;
            state.cells[cellIndex - 1].fluxcalmed = 1. * state.cells[cellIndex - 1].fluxcalmed;

            if (state.cells[cellIndex].temp < -50.)
                state.cells[cellIndex].temp = -50.;
            if (state.cells[cellIndex].temp > 200.)
                state.cells[cellIndex].temp = 200.;
        } else
            state.cells[cellIndex].temp = state.cells[cellIndex].calor.Textern1;
    } else { // case where the mixing speed is too low
        state.cells[cellIndex].temp = state.cells[cellIndex - 1].calor.Textern1;
        if (state.input.lingas == 1 && (cellIndex - 1 <= state.annulusTubingStart && cellIndex - 1 >= state.annulusTubingEnd)) {
            int j = state.annulusTubingStart + state.tubingAnnulusStart - (cellIndex - 1);
            state.cells[cellIndex - 1].calor.Vextern1 = 100.;
            state.cells[cellIndex - 1].calor.kextern1 = state.gasCells[j].calor.kint;
            state.cells[cellIndex - 1].calor.cpextern1 = state.gasCells[j].calor.cpint;
            state.cells[cellIndex - 1].calor.rhoextern1 = state.gasCells[j].calor.rhoint;
            state.cells[cellIndex - 1].calor.viscextern1 = state.gasCells[j].calor.viscint;
            double heatFlux = state.cells[cellIndex - 1].calor.transperm(0);
        }
    }
}

void advanceReverseSteadyTemperature(const ThermalState &state, int cellIndex, int rungeKuttaStage) {
    double cellLength = 0.5 * (state.cells[cellIndex].dx + state.cells[cellIndex + 1].dx);
    double meanCellLength = 0.5 * (state.cells[cellIndex].dx + state.cells[cellIndex + 1].dx);
    double meanTemperatureGradientCorrection = (state.cells[cellIndex].dx * state.cells[cellIndex].dTdLCor + state.cells[cellIndex + 1].dx * state.cells[cellIndex + 1].dTdLCor) / (state.cells[cellIndex].dx + state.cells[cellIndex + 1].dx);
    double diameter = state.cells[cellIndex + 1].duto.a;
    double flowArea = 0.25 * M_PI * diameter * diameter;
    double meanVoidFraction;
    double betmed;
    double meanPressure;
    double meanTemperature;
    if (rungeKuttaStage == 0) {
        meanVoidFraction = state.cells[cellIndex + 1].alf;
        betmed = state.cells[cellIndex + 1].bet;
        meanPressure = state.cells[cellIndex + 1].pres;
        meanTemperature = state.cells[cellIndex + 1].temp;
    } else {
        meanVoidFraction = state.cells[cellIndex].alf;
        betmed = state.cells[cellIndex].bet;
        meanPressure = state.cells[cellIndex].pres;
        meanTemperature = state.cells[cellIndex].temp;
    }
    double meanSuperficialGasVelocity;
    meanSuperficialGasVelocity = fabs(state.cells[cellIndex + 1].QG) / flowArea; // Superficial gas velocity
    double meanSuperficialLiquidVelocity;
    meanSuperficialLiquidVelocity = fabs(state.cells[cellIndex + 1].QL) / flowArea; // Superficial liquid velocity
    double mixtureFluxSign = 1.;
    if (fabs(meanSuperficialGasVelocity + meanSuperficialLiquidVelocity) > state.slowHeatTransferThreshold) {
        // Perform the thermal calculation for mixture velocities above 0.1 m/s.
        // Otherwise, assume the fluid temperature equals the ambient temperature.
        mixtureFluxSign = (meanSuperficialGasVelocity + meanSuperficialLiquidVelocity) / fabs(meanSuperficialGasVelocity + meanSuperficialLiquidVelocity);
        double interfaceMeanPressure = state.cells[cellIndex + 1].presaux - state.cells[cellIndex].dpB / 98066.5;
        double interfaceMeanTemperature;
        if (state.steadyIteration != 0 && state.input.AceleraConvergPerm == 0)
            // Set the left interface temperature from the previous iteration.
            // Use the adjacent left cell value to accelerate convergence.
            interfaceMeanTemperature = (state.cells[cellIndex].dx * state.cells[cellIndex].temp + state.cells[cellIndex].dxR * state.cells[cellIndex].tempR) / (state.cells[cellIndex].dx + state.cells[cellIndex].dxR);
        else
            interfaceMeanTemperature = state.cells[cellIndex + 1].temp;
        double rp = state.cells[cellIndex + 1].flui.MasEspLiq(meanPressure, meanTemperature);
        double rc = state.cells[cellIndex + 1].fluicol.MasEspFlu(meanPressure, meanTemperature);
        double liquidDensity = (1. - betmed) * rp + betmed * rc;
        double gasDensity = state.cells[cellIndex + 1].flui.MasEspGas(meanPressure, meanTemperature);
        double liquidSpecificHeat = (1. - betmed) * state.cells[cellIndex + 1].flui.CalorLiq(interfaceMeanPressure, interfaceMeanTemperature) + betmed * state.cells[cellIndex + 1].fluicol.CalorLiq(interfaceMeanPressure, interfaceMeanTemperature);
        double gasSpecificHeat = state.cells[cellIndex + 1].flui.CalorGas(interfaceMeanPressure, interfaceMeanTemperature);
        // Liquid Joule-Thomson coefficient multiplied by cp.
        // Gas Joule-Thomson coefficient multiplied by cp.p
        double liquidJouleThomson = (1. - betmed) * state.cells[cellIndex + 1].flui.JTL(interfaceMeanPressure, interfaceMeanTemperature) - betmed / rc;
        if (state.input.pocinjec > 0 && state.input.condpocinj.tipoFlui == 2) {
            liquidJouleThomson = -(1 + (interfaceMeanTemperature + 273.14) * state.cells[cellIndex + 1].fluicol.DrhoDtFlu(interfaceMeanPressure, interfaceMeanTemperature) / rc) / rc;
        }
        double gasJouleThomson = state.cells[cellIndex + 1].flui.JTG(interfaceMeanPressure, interfaceMeanTemperature);
        // potential energy:
        double hydrostaticPower = -(liquidDensity * meanSuperficialLiquidVelocity + gasDensity * meanSuperficialGasVelocity) * flowArea * 9.82 * sin(state.cells[cellIndex + 1].duto.teta);

        if (cellIndex > 120) {
            int debugStop;
            debugStop = 0;
        }

        // Definition of internal parameters in the piping to achieve heat exchange with the environment
        state.cells[cellIndex + 1].calor.Tint = interfaceMeanTemperature;
        state.cells[cellIndex + 1].calor.Vint = fabs(meanSuperficialGasVelocity + meanSuperficialLiquidVelocity);
        double liquidConductivity = (1. - betmed) * state.cells[cellIndex + 1].flui.CondLiq(interfaceMeanPressure, interfaceMeanTemperature) + betmed * state.cells[cellIndex + 1].fluicol.CondLiq(interfaceMeanPressure, interfaceMeanTemperature);
        state.cells[cellIndex + 1].calor.kint = liquidConductivity * (1 - meanVoidFraction) + state.cells[cellIndex + 1].flui.CondGas(interfaceMeanPressure, interfaceMeanTemperature) * meanVoidFraction;
        state.cells[cellIndex + 1].calor.cpint = liquidSpecificHeat * (1 - meanVoidFraction) + gasSpecificHeat * meanVoidFraction;
        state.cells[cellIndex + 1].calor.rhoint = liquidDensity * (1 - meanVoidFraction) + gasDensity * meanVoidFraction;
        double liquidViscosity = (1. - betmed) * state.cells[cellIndex + 1].flui.ViscOleo(interfaceMeanPressure, interfaceMeanTemperature) + betmed * state.cells[cellIndex + 1].fluicol.VisFlu(interfaceMeanPressure, interfaceMeanTemperature);
        state.cells[cellIndex + 1].calor.viscint = liquidViscosity * (1 - meanVoidFraction) * 1.e-3 + state.cells[cellIndex + 1].flui.ViscGas(interfaceMeanPressure, interfaceMeanTemperature) * meanVoidFraction * 1.e-3;

        double relaxedHeatFlux = 0;
        if (state.steadyIteration > 0)
            relaxedHeatFlux = state.cells[cellIndex + 1].fluxcalmed;
        double annulusResistance = 0.;
        double heatFlux;
        double gasHeatFlux;
        // checks if coupling with the annular space exists:
        if (state.input.lingas == 1 && (cellIndex + 1 <= state.annulusTubingStart && cellIndex + 1 >= state.annulusTubingEnd)) {
            int j = state.annulusTubingStart + state.tubingAnnulusStart - (cellIndex + 1);
            // in case there is coupling:

            if (state.steadyIteration == 0 && (cellIndex + 1) > state.annulusTubingEnd) {
                // On the first iteration, include casing, cement, and formation resistance.
                // Later iterations consider only heat transfer between the tubing and annular gas.
                // For the coupled model, obtain the external resistance from the annulus heat-transfer model.
                // Apply this up to the cell before the master valve to improve the ANM temperature estimate.
                // At the ANM cell, consider only heat transfer with the gas from the first iteration.
                state.gasCells[j].calor.Tint = interfaceMeanTemperature;
                state.gasCells[j].calor.Vint = 100;
                state.gasCells[j].calor.kint = state.cells[cellIndex + 1].flui.CondGas(interfaceMeanPressure, interfaceMeanTemperature);
                state.gasCells[j].calor.cpint = gasSpecificHeat;
                state.gasCells[j].calor.rhoint = gasDensity;
                state.gasCells[j].calor.viscint = state.cells[cellIndex + 1].flui.ViscGas(interfaceMeanPressure, interfaceMeanTemperature) * 1.e-3;
                state.gasCells[j].fluxcal = state.gasCells[j].calor.transperm(); // troca termica no anular
                annulusResistance = state.gasCells[j].calor.resGlob;                // resistencia das paredes
                // Only the thermal resistance from the casing to the formation is needed, not the actual heat flow.
                state.cells[cellIndex + 1].calor.Vextern1 = 100.;
                state.cells[cellIndex + 1].calor.kextern1 = state.gasCells[j].calor.kint;
                state.cells[cellIndex + 1].calor.cpextern1 = state.gasCells[j].calor.cpint;
                state.cells[cellIndex + 1].calor.rhoextern1 = state.gasCells[j].calor.rhoint;
                state.cells[cellIndex + 1].calor.viscextern1 = state.gasCells[j].calor.viscint;
            } else { // After the first iteration, consider only heat transfer between the tubing and annular gas.
                // Through the tubing wall.
                annulusResistance = 0.;
                state.cells[cellIndex + 1].calor.Vextern1 = state.gasCells[j].VGasR / state.gasCells[j].u1L;
                state.cells[cellIndex + 1].calor.kextern1 = state.gasCells[j].calor.kint;
                state.cells[cellIndex + 1].calor.cpextern1 = state.gasCells[j].calor.cpint;
                state.cells[cellIndex + 1].calor.rhoextern1 = state.gasCells[j].calor.rhoint;
                state.cells[cellIndex + 1].calor.viscextern1 = state.gasCells[j].calor.viscint;
            }
            if (state.steadyIteration == 0)
                state.cells[cellIndex + 1].calor.Textern1 = state.gasCells[j].calor.Textern1; // na primeira iteracao, como se
            // usa toda a resistencia termica do poco, a temperatura externa utilizada Ã© a geotermica
            else
                state.cells[cellIndex + 1].calor.Textern1 = state.gasCells[j].temp; // nas iteracoes seguintes, a temperatura ambiente
            // e a temperatura do gas
        }
        state.cells[cellIndex + 1].fluxcalmed = 0;
        heatFlux = mixtureFluxSign * state.cells[cellIndex + 1].calor.transperm(annulusResistance);
        state.cells[cellIndex + 1].fluxcalmed = heatFlux; // fluxo de calor na coluna

        double temperatureSpatialCoefficient = (liquidDensity * meanSuperficialLiquidVelocity * liquidSpecificHeat + gasDensity * meanSuperficialGasVelocity * gasSpecificHeat) * flowArea; // termo que multiplica
        // a derivada Dt/Dx
        double pressureSpatialCoefficient = 1. * (liquidDensity * meanSuperficialLiquidVelocity * liquidJouleThomson + gasDensity * meanSuperficialGasVelocity * gasJouleThomson) * flowArea; // termo que multiplica
        // a derivada Dp/Dx
        double pressureGradient;
        if ((state.cells[cellIndex + 1].acsr.tipo != 4 || state.cells[cellIndex + 1].acsr.bcs.freq < 1) && state.cells[cellIndex + 1].acsr.tipo != 7)
            // caso nao tenha BCS ou incremento de pressao  utiliza-se a pressao na fronteira esquerda
            //  e a pressao no centro de celula para o calculo de Dp/Dx
            pressureGradient = 2. * (state.cells[cellIndex + 1].presaux - state.cells[cellIndex + 1].pres) * 98066.5 / state.cells[cellIndex + 1].dx;
        else {
            // caso tenha BCS ou incremento de pressao  utiliza-se a pressao da celulaa esquerda
            //  e a pressao no centro de celula para o calculoi de Dp/Dx
            pressureGradient = (interfaceMeanPressure - state.cells[cellIndex + 1].pres) * 98600. / cellLength;
        }
        state.cells[cellIndex].VTemper = meanSuperficialLiquidVelocity; // esta velocidade so e util no caso transiente, Ã© armazenada aqui
        // apenas para se ter um valor quando a simulacao transiente se iniciar
        double temperatureGradient = (-state.cells[cellIndex + 1].temp) / meanCellLength;

        double kineticTerm = 0;
        double upstreamMeanGasVelocity = 0;
        double upstreamMeanLiquidVelocity = 0;
        double meanGasVelocity = 0;
        double meanLiquidVelocity = 0;
        // termo de energia cinetica:
        if (state.cells[cellIndex].acsr.tipo == 0 && state.cells[cellIndex + 1].acsr.tipo == 0 && cellIndex < state.lastCell) {
            double kineticCellLength = state.cells[cellIndex + 1].dx;
            double upstreamDiameter = state.cells[cellIndex + 1].duto.a;
            double upstreamFlowArea = 0.25 * M_PI * upstreamDiameter * upstreamDiameter;

            if (state.cells[cellIndex + 1].alf > 1e-3)
                meanGasVelocity = meanSuperficialGasVelocity / state.cells[cellIndex + 1].alf;
            if (state.cells[cellIndex + 1].alf < (1. - 1e-3))
                meanLiquidVelocity = meanSuperficialLiquidVelocity / (1. - state.cells[cellIndex + 1].alf);

            if (state.cells[cellIndex].alf > 1e-3) {
                upstreamMeanGasVelocity = fabs(state.cells[cellIndex].QG) / (upstreamFlowArea);
                upstreamMeanGasVelocity /= state.cells[cellIndex].alf;
            }
            if (state.cells[cellIndex].alf < (1. - 1e-3)) {
                upstreamMeanLiquidVelocity = fabs(state.cells[cellIndex].QL) / (upstreamFlowArea);
                upstreamMeanLiquidVelocity /= (1. - state.cells[cellIndex].alf);
            }

            kineticTerm = -fabs(state.cells[cellIndex].MC - state.cells[cellIndex].Mliqini) * meanGasVelocity * (meanGasVelocity - upstreamMeanGasVelocity) / kineticCellLength -
                       fabs(state.cells[cellIndex].Mliqini) * meanLiquidVelocity * (meanLiquidVelocity - upstreamMeanLiquidVelocity) / kineticCellLength;
        }

        double gasMassSourceTerm = 0.;
        double liquidMassSourceTerm = 0.;
        double sourceTemperature = state.cells[cellIndex + 1].temp;
        double sourceGasSpecificHeat;
        double sourceSpecificHeatRatio = 0.;
        double sourceLiquidSpecificHeat;

        // calculo da energia adicionada no sistema devido a fontes de massa
        if (state.cells[cellIndex + 1].acsr.tipo == 1) { // caso fonte de gas
            sourceTemperature = state.cells[cellIndex + 1].acsr.injg.temp;
            sourceGasSpecificHeat = state.cells[cellIndex + 1].acsr.injg.FluidoPro.CalorGas(meanPressure, sourceTemperature);
            sourceSpecificHeatRatio = state.cells[cellIndex + 1].acsr.injg.FluidoPro.ConstAdG(meanPressure, sourceTemperature);
            sourceLiquidSpecificHeat = 0.;
        } else if (state.cells[cellIndex + 1].acsr.tipo == 2) { // caso fonte de liquido
            sourceTemperature = state.cells[cellIndex + 1].acsr.injl.temp;
            sourceGasSpecificHeat = 0.;
            sourceSpecificHeatRatio = 1.;
            sourceLiquidSpecificHeat = (1. - state.cells[cellIndex + 1].acsr.injl.bet) * state.cells[cellIndex + 1].acsr.injl.FluidoPro.CalorLiq(meanPressure, meanTemperature) + state.cells[cellIndex + 1].acsr.injl.bet * state.cells[cellIndex + 1].acsr.injl.fluidocol.CalorLiq(meanPressure, meanTemperature);
        } else if (state.cells[cellIndex + 1].acsr.tipo == 3) { // caso IPR
            sourceTemperature = state.cells[cellIndex + 1].acsr.ipr.Tres;
            sourceGasSpecificHeat = state.cells[cellIndex + 1].acsr.ipr.FluidoPro.CalorGas(meanPressure, meanTemperature);
            sourceSpecificHeatRatio = state.cells[cellIndex + 1].acsr.ipr.FluidoPro.ConstAdG(meanPressure, meanTemperature);
            sourceLiquidSpecificHeat = state.cells[cellIndex + 1].acsr.ipr.FluidoPro.CalorLiq(meanPressure, meanTemperature);
        } else if (state.cells[cellIndex + 1].acsr.tipo == 15) { // caso IPR
            sourceTemperature = state.cells[cellIndex + 1].acsr.radialPoro.tRes;
            sourceGasSpecificHeat = state.cells[cellIndex + 1].acsr.radialPoro.flup.CalorGas(meanPressure, meanTemperature);
            sourceSpecificHeatRatio = state.cells[cellIndex + 1].acsr.radialPoro.flup.ConstAdG(meanPressure, meanTemperature);
            sourceLiquidSpecificHeat = state.cells[cellIndex + 1].acsr.radialPoro.flup.CalorLiq(meanPressure, meanTemperature);
        } else if (state.cells[cellIndex + 1].acsr.tipo == 16) { // caso IPR
            sourceTemperature = state.cells[cellIndex + 1].acsr.poroso2D.dados.tRes;
            sourceGasSpecificHeat = state.cells[cellIndex + 1].acsr.poroso2D.dados.flup.CalorGas(meanPressure, meanTemperature);
            sourceSpecificHeatRatio = state.cells[cellIndex + 1].acsr.poroso2D.dados.flup.ConstAdG(meanPressure, meanTemperature);
            sourceLiquidSpecificHeat = state.cells[cellIndex + 1].acsr.poroso2D.dados.flup.CalorLiq(meanPressure, meanTemperature);
        } else if (state.cells[cellIndex + 1].acsr.tipo == 9 && state.cells[cellIndex + 1].acsr.fontechk.abertura > 1e-6 &&
                   (state.cells[cellIndex + 1].fontemassCR + state.cells[cellIndex + 1].fontemassGR + state.cells[cellIndex + 1].fontemassLR) > 1e-9) {
            // caso vazamento
            sourceTemperature = state.cells[cellIndex + 1].acsr.fontechk.tamb;
            sourceGasSpecificHeat = state.cells[cellIndex + 1].acsr.fontechk.fluidoPamb.CalorGas(meanPressure, meanTemperature);
            sourceSpecificHeatRatio = state.cells[cellIndex + 1].acsr.fontechk.fluidoPamb.ConstAdG(meanPressure, meanTemperature);
            sourceLiquidSpecificHeat = (1. - state.cells[cellIndex + 1].acsr.fontechk.betISamb) *
                       state.cells[cellIndex + 1].acsr.fontechk.fluidoPamb.CalorLiq(meanPressure, meanTemperature) +
                   state.cells[cellIndex + 1].acsr.fontechk.betISamb * state.cells[cellIndex + 1].acsr.fontechk.fluidocol.CalorLiq(meanPressure, meanTemperature);
        } else if (state.cells[cellIndex + 1].acsr.tipo != 0) {

            sourceGasSpecificHeat = 0.;
            sourceSpecificHeatRatio = 1.;
            sourceLiquidSpecificHeat = 0.;
        } else {
            sourceGasSpecificHeat = 0.;
            sourceSpecificHeatRatio = 1.;
            sourceLiquidSpecificHeat = 0.;
        }

        liquidMassSourceTerm = 0;
        if (state.cells[cellIndex + 1].fontemassLR > 0.)
            liquidMassSourceTerm = state.cells[cellIndex + 1].fontemassLR / cellLength;
        if (state.cells[cellIndex + 1].fontemassCR > 0.)
            liquidMassSourceTerm += state.cells[cellIndex + 1].fontemassCR / cellLength;
        liquidMassSourceTerm *= (sourceLiquidSpecificHeat) * (sourceTemperature - state.cells[cellIndex + 1].temp);

        gasMassSourceTerm = state.cells[cellIndex + 1].fontemassGR / cellLength;
        if (gasMassSourceTerm > 0.)
            gasMassSourceTerm *= (sourceGasSpecificHeat / sourceSpecificHeatRatio) * (sourceTemperature - state.cells[cellIndex + 1].temp);
        else
            gasMassSourceTerm = 0;

        // efeito do calor latentHeatTerm, quando este for solicitado
        double latentHeatTerm = 0.;
        if (state.cells[cellIndex].flui.dVaporMassFraction < (1 - 1e-15) && state.cells[cellIndex].flui.dVaporMassFraction > (1e-15)) {
            if (state.cells[cellIndex + 1].acsr.tipo == 1 || state.cells[cellIndex + 1].acsr.tipo == 2 || state.cells[cellIndex + 1].acsr.tipo == 3 || state.cells[cellIndex + 1].acsr.tipo == 15 || state.cells[cellIndex + 1].acsr.tipo == 16)
                state.cells[cellIndex + 1].FonteMudaFase =
                    0.;
            if (state.latentHeatEnabled > 0 && state.steadyIteration != 0 && state.input.flashCompleto == 0) {
                latentHeatTerm = -interpolateLatentHeat(state, meanPressure, meanTemperature) * state.cells[cellIndex + 1].FonteMudaFase;
            } else if ((state.input.flashCompleto == 1 || state.input.flashCompleto == 2) && state.latentHeatEnabled > 0 && state.steadyIteration != 0) {
                latentHeatTerm = -(state.cells[cellIndex].flui.EntalpGas(meanPressure, meanTemperature) -
                            state.cells[cellIndex].flui.EntalpLiq(meanPressure, meanTemperature)) *
                          state.cells[cellIndex + 1].FonteMudaFase;
            } else
                latentHeatTerm = 0;
        }

        //////trecho sem utilidade////////////////////////////////////////////////////////
        double interfaceVoidFraction;
        double leftInterfaceVoidFraction;
        if (meanSuperficialGasVelocity > 0) {
            interfaceVoidFraction = state.cells[cellIndex + 1].alf;
            leftInterfaceVoidFraction = state.cells[cellIndex + 1].alfL;
        } else {
            interfaceVoidFraction = state.cells[cellIndex + 1].alfR;
            leftInterfaceVoidFraction = state.cells[cellIndex + 1].alf;
        }
        double slipVelocity;
        if (interfaceVoidFraction > (*state.globals).localtiny && interfaceVoidFraction < (1. - (*state.globals).localtiny))
            slipVelocity = meanSuperficialGasVelocity / interfaceVoidFraction - meanSuperficialLiquidVelocity / (1. - interfaceVoidFraction);
        else if (interfaceVoidFraction > (*state.globals).localtiny)
            slipVelocity = meanSuperficialGasVelocity;
        else
            slipVelocity = meanSuperficialLiquidVelocity;
        double interfacialWorkTerm = flowArea * state.cells[cellIndex + 1].pres * 98600 * slipVelocity * (interfaceVoidFraction - leftInterfaceVoidFraction) / state.cells[cellIndex + 1].dx;

        if (fabs(temperatureSpatialCoefficient) > (*state.globals).localtiny) {
            // parecela de energia relacionada ao ytrabalho de fronteira, energia potencial,
            // energia cinetica, fontes de massa, calor latentHeatTerm e trabalho de eixo:
            double sourceTemperatureGradient = meanTemperatureGradientCorrection * (pressureSpatialCoefficient * pressureGradient - kineticTerm - (hydrostaticPower) + (liquidMassSourceTerm + gasMassSourceTerm) - latentHeatTerm - state.cells[cellIndex + 1].potBT / meanCellLength) / temperatureSpatialCoefficient;
            // parcela de energia relacionada aa troca termica
            double heatFluxTemperatureGradient = meanTemperatureGradientCorrection * (heatFlux) / temperatureSpatialCoefficient;

            // e feita uma avaliacao se a troca termica esta se dando de maneira muito rapida
            // como avanco da temperatura e explicita, isto pode levar a instabilidade no calculo termico
            // se for verificado que a troca termica esta ocorrendo de maneira rapida, o avanco e
            // feito em um numero maior de passos de uma celula para outra
            int subStepCount;
            double subStepLength;
            double stabilityLength = fabs(temperatureSpatialCoefficient) / meanTemperatureGradientCorrection;
            if (meanCellLength / (state.cells[cellIndex + 1].calor.resGlob + annulusResistance) < (stabilityLength + 0. * 1000.)) {
                subStepCount = 1;
                subStepLength = meanCellLength;
            } else { // resistencia termica e pequena, se determinara em quantos passos se dara
                // o avanco de temperatrura
                subStepCount = (meanCellLength / (state.cells[cellIndex + 1].calor.resGlob + annulusResistance)) / (stabilityLength + 0 * 1000) + 1;
                subStepLength = meanCellLength / subStepCount; // celula divida em npassos
            }
            double subStepTemperature = state.cells[cellIndex + 1].temp;

            subStepTemperature = subStepLength * (-(-state.cells[cellIndex + 1].temp) / subStepLength + sourceTemperatureGradient + heatFluxTemperatureGradient); // primeiro avanco
            for (int j = 1; j < subStepCount; j++) {                                                  // avancos seguintes
                // as pacelas de energia sao mantidas, com excessao do fluxo de calor que e
                // reavalkiado a cada passo:
                state.cells[cellIndex + 1].calor.Tint = subStepTemperature;
                state.cells[cellIndex + 1].calor.Vint = meanSuperficialGasVelocity + meanSuperficialLiquidVelocity;
                liquidConductivity = (1. - betmed) * state.cells[cellIndex + 1].flui.CondLiq(meanPressure, subStepTemperature) + betmed * state.cells[cellIndex + 1].fluicol.CondLiq(meanPressure, subStepTemperature);
                state.cells[cellIndex + 1].calor.kint = liquidConductivity * (1 - meanVoidFraction) + state.cells[cellIndex + 1].flui.CondGas(meanPressure, subStepTemperature) * meanVoidFraction;
                state.cells[cellIndex + 1].calor.cpint = liquidSpecificHeat * (1 - meanVoidFraction) + gasSpecificHeat * meanVoidFraction;
                state.cells[cellIndex + 1].calor.rhoint = liquidDensity * (1 - meanVoidFraction) + gasDensity * meanVoidFraction;
                liquidViscosity = (1. - betmed) * state.cells[cellIndex + 1].flui.ViscOleo(meanPressure, subStepTemperature) + betmed * state.cells[cellIndex + 1].fluicol.VisFlu(meanPressure, subStepTemperature);
                state.cells[cellIndex + 1].calor.viscint = liquidViscosity * (1 - meanVoidFraction) * 1.e-3 + state.cells[cellIndex + 1].flui.ViscGas(meanPressure, subStepTemperature) * meanVoidFraction * 1.e-3;
                if (state.steadyIteration != 0 && state.input.lingas == 1 && (cellIndex + 1 <= state.annulusTubingStart && cellIndex + 1 >= state.annulusTubingEnd)) {
                    int k = state.annulusTubingStart + state.tubingAnnulusStart - (cellIndex + 1);
                    double externalTemperatureStep = (state.gasCells[k].temp - state.gasCells[k - 1].temp) / subStepCount;
                    state.cells[cellIndex + 1].calor.Textern1 = state.gasCells[k].temp + (j - 1) * externalTemperatureStep;
                }
                heatFlux = mixtureFluxSign * state.cells[cellIndex + 1].calor.transperm(annulusResistance);
                state.cells[cellIndex + 1].fluxcalmed += heatFlux;
                heatFluxTemperatureGradient = meanTemperatureGradientCorrection * (heatFlux) / temperatureSpatialCoefficient;
                subStepTemperature = subStepLength * (-(-subStepTemperature) / subStepLength + sourceTemperatureGradient + heatFluxTemperatureGradient);
            }

            state.cells[cellIndex].temp = subStepTemperature;
            state.cells[cellIndex + 1].fluxcalmed /= subStepCount;
            state.cells[cellIndex + 1].fluxcalmed = 1. * state.cells[cellIndex + 1].fluxcalmed;

            if (state.cells[cellIndex].temp < -50.)
                state.cells[cellIndex].temp = -50.;
            if (state.cells[cellIndex].temp > 200.)
                state.cells[cellIndex].temp = 200.;
        } else
            state.cells[cellIndex].temp = state.cells[cellIndex].calor.Textern1;
    } else { // caso em que a velocidade da mistura e muito baixa
        state.cells[cellIndex].temp = state.cells[cellIndex + 1].calor.Textern1;
        if (state.input.lingas == 1 && (cellIndex + 1 <= state.annulusTubingStart && cellIndex + 1 >= state.annulusTubingEnd)) {
            int j = state.annulusTubingStart + state.tubingAnnulusStart - (cellIndex + 1);
            state.cells[cellIndex + 1].calor.Vextern1 = 100.;
            state.cells[cellIndex + 1].calor.kextern1 = state.gasCells[j].calor.kint;
            state.cells[cellIndex + 1].calor.cpextern1 = state.gasCells[j].calor.cpint;
            state.cells[cellIndex + 1].calor.rhoextern1 = state.gasCells[j].calor.rhoint;
            state.cells[cellIndex + 1].calor.viscextern1 = state.gasCells[j].calor.viscint;
            double heatFlux = state.cells[cellIndex + 1].calor.transperm(0);
        }
    }
}

void computeGasTemperature(const ThermalState &state, int cellIndex, double previousTemperature, int steadyMode) {

    if (state.thermalSourceDisabled == 0) {
        double cellLength = state.gasCells[cellIndex].dx0;
        double meanCellLength = 0.5 * (state.gasCells[cellIndex].dx0 + state.gasCells[cellIndex - 1].dx0);
        double flowArea = state.gasCells[cellIndex].duto.area;
        double meanSuperficialGasVelocity;
        if (cellIndex < state.gasCellCount)
            meanSuperficialGasVelocity = state.gasCells[cellIndex].VGasR / state.gasCells[cellIndex].u1L;
        else {
            meanSuperficialGasVelocity = state.gasCells[cellIndex].VGasL / state.gasCells[cellIndex].u1L;
        }
        double gasDensity = state.gasCells[cellIndex].rg;
        double gasSpecificHeat = state.gasCells[cellIndex].flui.CalorGas(state.gasCells[cellIndex].presini, state.gasCells[cellIndex].tempini);
        double gasSpecificHeatConstantVolume = state.gasCells[cellIndex].flui.CalorGasVolMod(state.gasCells[cellIndex].presini, state.gasCells[cellIndex].tempini);
        double gasJouleThomson = state.gasCells[cellIndex].flui.JTG(state.gasCells[cellIndex].presini, state.gasCells[cellIndex].tempini);
        double hydrostaticPower = (gasDensity * meanSuperficialGasVelocity) * flowArea * 9.82 * sin(state.gasCells[cellIndex].duto.teta);

        state.gasCells[cellIndex].calor.Tint = state.gasCells[cellIndex].tempini;
        state.gasCells[cellIndex].calor.dtL = state.gasCells[cellIndex].tempini - state.gasCells[cellIndex - 1].tempini;
        state.gasCells[cellIndex].calor.Vint = meanSuperficialGasVelocity;
        state.gasCells[cellIndex].calor.dt = state.gasCells[cellIndex].dt;
        state.gasCells[cellIndex].calor.kint = state.gasCells[cellIndex].flui.CondGas(state.gasCells[cellIndex].presini, state.gasCells[cellIndex].tempini);
        state.gasCells[cellIndex].calor.cpint = gasSpecificHeat;
        state.gasCells[cellIndex].calor.rhoint = gasDensity;
        state.gasCells[cellIndex].calor.viscint = state.gasCells[cellIndex].flui.ViscGas(state.gasCells[cellIndex].presini, state.gasCells[cellIndex].tempini) * 1.e-3;
        double temperaturePerturbation = state.gasCells[cellIndex].temp * 0.01;
        if (fabs(state.gasCells[cellIndex].temp) < 1e-15)
            temperaturePerturbation = 0.1;
        double gasDensityChange = state.gasCells[cellIndex].flui.MasEspGas(state.gasCells[cellIndex].presini, state.gasCells[cellIndex].tempini + temperaturePerturbation) - gasDensity;
        state.gasCells[cellIndex].calor.betint = -(1 / state.gasCells[cellIndex].calor.rhoint) * gasDensityChange / (temperaturePerturbation);
        if (steadyMode == 0)
            state.gasCells[cellIndex].fluxcal = state.gasCells[cellIndex].calor.transtrans();
        else
            state.gasCells[cellIndex].fluxcal = state.gasCells[cellIndex].calor.transperm();
        if (cellIndex >= state.tubingAnnulusStart && cellIndex <= state.tubingAnnulusEnd && state.networkCoupled == 1) {
            int kconecte = cellIndex - state.tubingAnnulusStart;
            int iconecte = state.annulusTubingStart - kconecte;
            state.gasCells[cellIndex].fluxcal -= state.cells[iconecte].calor.fluxFim;
        }

        double lengthRatio;
        if (cellIndex < state.gasCellCount)
            lengthRatio = cellLength / (cellLength + state.gasCells[cellIndex + 1].dx0);
        else
            lengthRatio = cellLength / (cellLength + state.gasCells[cellIndex - 1].dx0);
        double timeCoefficient = gasDensity * gasSpecificHeatConstantVolume * flowArea;
        double pressureTimeCoefficient = state.gasCells[cellIndex].flui.CalorGasPresMod(state.gasCells[cellIndex].presini, state.gasCells[cellIndex].tempini, state.gasCells[cellIndex].rg) *
                               (gasDensity * flowArea);

        double temperatureSpatialCoefficient = gasDensity * meanSuperficialGasVelocity * gasSpecificHeat * flowArea;
        double pressureSpatialCoefficient = gasDensity * meanSuperficialGasVelocity * gasJouleThomson * flowArea;
        double pressureGradient;
        if (cellIndex < state.gasCellCount)
            pressureGradient = 2. * (((1 - lengthRatio) * state.gasCells[cellIndex + 1].presini + lengthRatio * state.gasCells[cellIndex].presini) - state.gasCells[cellIndex].presini) * 98066.5 / cellLength;
        else
            pressureGradient = 2. * (state.gasCells[cellIndex].presini - ((1 - lengthRatio) * state.gasCells[cellIndex - 1].presini + lengthRatio * state.gasCells[cellIndex].presini)) * 98066.5 / cellLength;
        double temperatureGradient = (state.gasCells[cellIndex].tempini - state.gasCells[cellIndex - 1].tempini) / meanCellLength;
        if (cellIndex < state.gasCellCount)
            if (meanSuperficialGasVelocity < 0)
                temperatureGradient = (state.gasCells[cellIndex + 1].tempini - state.gasCells[cellIndex].tempini) / meanCellLength;
        if ((cellIndex == 1 && meanSuperficialGasVelocity <= 0) || (cellIndex == state.gasCellCount && meanSuperficialGasVelocity <= 0))
            temperatureGradient = 0.;

        double kineticTerm;
        double mixtureFluxDifference = 0.;
        double mixtureDensity = gasDensity;
        double upstreamMeanSuperficialGasVelocity = state.gasCells[cellIndex].VGasL / state.gasCells[cellIndex - 1].u1L;
        mixtureFluxDifference = (meanSuperficialGasVelocity - upstreamMeanSuperficialGasVelocity) / cellLength;

        kineticTerm = mixtureDensity * flowArea * meanSuperficialGasVelocity * meanSuperficialGasVelocity * mixtureFluxDifference;

        double gasMassSourceTerm = 0.;
        double liquidMassSourceTerm = 0.;

        double pressureWorkFactor = 1.;
        if ((*state.globals).lixo5 < 1000.)
            pressureWorkFactor = 1.;

        state.gasCells[cellIndex].temp = ((timeCoefficient / state.gasCells[cellIndex].dt) * state.gasCells[cellIndex].temp - (pressureWorkFactor) * (pressureTimeCoefficient * (state.gasCells[cellIndex].pres - state.gasCells[cellIndex].presini) * 98066.5 / state.gasCells[cellIndex].dt) + state.gasCells[cellIndex].dTdLCor * (-temperatureSpatialCoefficient * temperatureGradient + pressureSpatialCoefficient * pressureGradient - kineticTerm - hydrostaticPower + liquidMassSourceTerm + gasMassSourceTerm + state.gasCells[cellIndex].fluxcal)) / (timeCoefficient / state.gasCells[cellIndex].dt);

        if (state.gasCells[cellIndex].temp < -50.)
            state.gasCells[cellIndex].temp = -50.;
        if (state.gasCells[cellIndex].temp > 200.)
            state.gasCells[cellIndex].temp = 200.;

        if (cellIndex > 0)
            state.gasCells[cellIndex - 1].tempR = state.gasCells[cellIndex].temp;
        if (cellIndex < state.gasCellCount)
            state.gasCells[cellIndex + 1].tempL = state.gasCells[cellIndex].temp;
    } else {
        state.gasCells[cellIndex].temp = state.gasCells[cellIndex].calor.Textern1;

        if (cellIndex > 0)
            state.gasCells[cellIndex - 1].tempR = state.gasCells[cellIndex].temp;
        if (cellIndex < state.gasCellCount)
            state.gasCells[cellIndex + 1].tempL = state.gasCells[cellIndex].temp;
    }
}

void computeDischargeTemperature(const ThermalState &state, int cellIndex) {
    double leftHalfLength = 0.5 * state.gasCells[cellIndex].dxL;
    double rightHalfLength = 0.5 * state.gasCells[cellIndex].dx0;
    double rightGasFraction = 0.;
    if (state.gasCells[cellIndex].razInter <= 0.5)
        rightGasFraction = 2 * state.gasCells[cellIndex].razInter;
    double leftGasFraction = 0.;
    if (state.gasCells[cellIndex - 1].razInter >= 0.5)
        leftGasFraction = 2 * (state.gasCells[cellIndex - 1].razInter - 0.5);
    double leftGasLength = leftHalfLength * leftGasFraction;
    double rightGasLength = rightHalfLength * rightGasFraction;
    double leftLiquidLength = leftHalfLength - leftGasLength;
    double rightLiquidLength = rightHalfLength - rightGasLength;
    double totalLength = leftLiquidLength + rightLiquidLength + leftGasLength + rightGasLength;
    double diameter = state.gasCells[cellIndex - 1].duto.a;
    double flowArea = 0.25 * M_PI * diameter * diameter;
    double pressure;
    double temperature;

    pressure = state.gasCells[cellIndex - 1].pres;
    temperature = state.gasCells[cellIndex - 1].temp;
    double mixtureDensity = ((leftLiquidLength + rightLiquidLength) * state.gasCells[cellIndex].MasEspFlu(pressure, temperature) + (leftGasLength + rightGasLength) * state.gasCells[cellIndex].flui.MasEspGas(pressure, temperature)) / totalLength;

    double interfaceVelocity = state.gasCells[cellIndex].VGasL / (state.gasCells[cellIndex].MasEspFlu(pressure, temperature) * state.gasCells[cellIndex - 1].duto.area);
    if (state.gasCells[cellIndex].razInter > (*state.globals).localtiny)
        interfaceVelocity = state.gasCells[cellIndex].VGasL / (state.gasCells[cellIndex].flui.MasEspGas(pressure, temperature) * state.gasCells[cellIndex - 1].duto.area);
    double liquidSpecificHeat = ((leftLiquidLength + rightLiquidLength) * state.gasCells[cellIndex].CalorLiq(pressure, temperature) + (leftGasLength + rightGasLength) * state.gasCells[cellIndex].flui.CalorGas(pressure, temperature)) / totalLength;

    state.gasCells[cellIndex - 1].calor.Tint = temperature;
    state.gasCells[cellIndex - 1].calor.Vint = interfaceVelocity;
    double liquidConductivity = ((leftLiquidLength + rightLiquidLength) * state.gasCells[cellIndex - 1].CondLiq(pressure, temperature) + (leftGasLength + rightGasLength) * state.gasCells[cellIndex].flui.CondGas(pressure, temperature)) / totalLength;
    state.gasCells[cellIndex - 1].calor.kint = liquidConductivity;
    state.gasCells[cellIndex - 1].calor.cpint = liquidSpecificHeat;
    state.cells[cellIndex - 1].calor.rhoint = mixtureDensity;
    double liquidViscosity = (((leftLiquidLength + rightLiquidLength) * state.gasCells[cellIndex].VisFlu(pressure, temperature) + (leftGasLength + rightGasLength) * state.gasCells[cellIndex].flui.ViscGas(pressure, temperature)) * 1e-3) / totalLength;
    state.cells[cellIndex - 1].calor.viscint = liquidViscosity;

    double heatFlux = state.cells[cellIndex - 1].calor.transtrans();
    if ((cellIndex - 1) >= state.tubingAnnulusStart && (cellIndex - 1) <= state.tubingAnnulusEnd && state.networkCoupled == 1) {
        int kconecte = (cellIndex - 1) - state.tubingAnnulusStart;
        int iconecte = state.annulusTubingStart - kconecte;
        heatFlux -= state.cells[iconecte].calor.fluxFim;
    }

    state.gasCells[cellIndex - 1].tempR = state.gasCells[cellIndex].temp;
    state.gasCells[cellIndex].tempL = state.gasCells[cellIndex - 1].temp;
}

double computeGasLiftDischargeTemperature(const ThermalState &state, int valveIndex) {
    int stepCount = floor((state.gasLiftChokes[valveIndex].presEstag - state.gasLiftChokes[valveIndex].presGarg) / 50) + 1;
    double pressureStep = -(state.gasLiftChokes[valveIndex].presEstag - state.gasLiftChokes[valveIndex].presGarg) / stepCount;
    double stageInletPressure = state.gasLiftChokes[valveIndex].presEstag;
    double stageOutletPressure = state.gasLiftChokes[valveIndex].presEstag + pressureStep;
    double stageInletTemperature = state.gasLiftChokes[valveIndex].tempEstag;
    double stageOutletTemperature;
    for (int cellIndex = 0; cellIndex < stepCount; cellIndex++) {
        double gasSpecificHeat = state.gasLiftChokes[valveIndex].flui.CalorGas(stageInletPressure, stageInletTemperature);
        double compressibilityTemperatureDerivative = state.gasLiftChokes[valveIndex].flui.DZDT(stageInletPressure, stageInletTemperature);
        double inletTemperatureKelvin = stageInletTemperature + 273.23;
        stageOutletTemperature = 1.0 / (1.0 / inletTemperatureKelvin - ((286.998 / state.gasLiftChokes[valveIndex].flui.Deng) * compressibilityTemperatureDerivative / gasSpecificHeat) * log((stageOutletPressure) / (stageInletPressure)));
        stageInletPressure = stageOutletPressure;
        stageOutletPressure = stageOutletPressure - pressureStep;
        stageInletTemperature = stageOutletTemperature - 273.23;
    }
    if (stepCount == 0) {
        if ((*state.globals).lixo5 > 1000) {
            int debugStop;
            debugStop = 0;
        }
        double gasSpecificHeat = state.gasLiftChokes[valveIndex].flui.CalorGas(stageInletPressure, stageInletTemperature);
        double gasJouleThomson = state.gasLiftChokes[valveIndex].flui.JTG(stageInletPressure, stageInletTemperature) / gasSpecificHeat;
        stageInletTemperature -= gasJouleThomson * (state.gasLiftChokes[valveIndex].presEstag - state.gasLiftChokes[valveIndex].presGarg) * 98066.52;
    }
    return stageInletTemperature;
}

void updateProductionTemperaturePeriphery(const ThermalState &state, int cellIndex) {
    if (cellIndex > 0)
        state.cells[cellIndex - 1].tempR = state.cells[cellIndex].temp;
    if (cellIndex < state.lastCell)
        state.cells[cellIndex + 1].tempL = state.cells[cellIndex].temp;
    state.cells[cellIndex].tempini = state.cells[cellIndex].temp;
}

void computeOutletTemperature(const ThermalState &state) {

    if (state.input.chokep.abertura[0] <= 0.6 && state.input.chokep.abertura[0] > (*state.globals).localtiny && state.outletPressure < state.gasSurfacePressure) {
        double inletMassFlow = state.cells[state.lastCell - 1].MR;
        double gasMassFlow = state.cells[state.lastCell - 1].MR - state.cells[state.lastCell - 1].MliqiniR;
        double liquidDensity = state.cells[state.lastCell].flui.MasEspLiq(state.cells[state.lastCell].pres, state.cells[state.lastCell].temp);
        double rholc = state.cells[state.lastCell].fluicol.MasEspFlu(state.cells[state.lastCell].pres, state.cells[state.lastCell].temp);
        double betEF = state.cells[state.lastCell].bet;
        double quality = fabs(gasMassFlow / inletMassFlow);

        double liquidJouleThomson = (1. - betEF) * state.cells[state.lastCell].flui.JTL(state.cells[state.lastCell].pres, state.cells[state.lastCell].temp) - betEF / rholc; // alteraacao2
        double gasJouleThomson = state.cells[state.lastCell].flui.JTG(state.cells[state.lastCell].pres, state.cells[state.lastCell].temp);
        state.surfaceTemperature = state.cells[state.lastCell].temp + ((1. - quality) * liquidJouleThomson + quality * gasJouleThomson) * (state.cells[state.lastCell].pres - state.cells[state.lastCell].pres); //????????
                                                                                                                  //???????????????????????????????celula[ncel].pres - celula[ncel].pres????????????????????????????????????????????

    } else
        state.surfaceTemperature = state.cells[state.lastCell - 1].temp;
}

}  // namespace sisprod::thermal
