#include "SisProdThermal.h"

#include "SisProdConstants.h"

#include "Leitura.h"
#include "celula3.h"
#include "solver3DPoisson.h"

#include <math.h>

namespace sisprod::thermal {

namespace {

/// Bisects a monotonically increasing table axis, returning the lower end of the
/// interval bracketing `value`. `axis(i)` yields the i-th coordinate.
///
/// `descendProbeOffset` is 0 for the pressure axes and -1 for the temperature
/// axis; they pick different intervals at grid points, so it is not optional.
/// The midpoint and the == guards must stay as written for the same reason.
template <typename Axis>
[[nodiscard]] int bracketingIndex(double value, int divisionCount, Axis &&axis,
                                  int descendProbeOffset) {
    int index = 0;
    int searchLow = 1;
    int searchHigh = divisionCount + 1;
    int searchMiddle;
    while (searchLow <= searchHigh) {
        searchMiddle = (searchLow + searchHigh) / 2;
        if (searchMiddle == 1) {
            index = searchMiddle;
            break;
        } else if (searchMiddle == divisionCount + 1 && axis(searchMiddle) == value) {
            index = searchMiddle - 1;
            break;
        }
        if (axis(searchMiddle) > value && axis(searchMiddle - 1) <= value) {
            index = searchMiddle - 1;
            break;
        }
        if (axis(searchMiddle + descendProbeOffset) < value)
            searchLow = searchMiddle + 1;
        else
            searchHigh = searchMiddle - 1;
    }
    return index;
}

}  // namespace


double interpolateLatentHeat(const ThermalState &state, double pressure, double temperature) {
    int divisionCount = state.input.tabent.npont - 1;
    int pressureIndex = 0.;
    int temperatureIndex = 0.;
    double latentHeat;
    if (pressure < state.latentHeatTable[1][0] || pressure >= state.latentHeatTable[divisionCount + 1][0] || temperature < state.latentHeatTable[0][1] || temperature >= state.latentHeatTable[0][divisionCount + 1])
        latentHeat = 0.;

    else {
        pressureIndex = bracketingIndex(
            pressure, divisionCount,
            [&](int index) { return state.latentHeatTable[index][0]; }, 0);
        temperatureIndex = bracketingIndex(
            temperature, divisionCount,
            [&](int index) { return state.latentHeatTable[0][index]; }, -1);
        double pressureRatio = (state.latentHeatTable[pressureIndex][0] - pressure) / (state.latentHeatTable[pressureIndex][0] - state.latentHeatTable[pressureIndex + 1][0]);
        double temperatureRatio = (state.latentHeatTable[0][temperatureIndex] - temperature) / (state.latentHeatTable[0][temperatureIndex] - state.latentHeatTable[0][temperatureIndex + 1]);
        double latentHeatAtTemperatureIndex = (1 - pressureRatio) * (state.latentHeatTable[pressureIndex][temperatureIndex]) + pressureRatio * (state.latentHeatTable[pressureIndex + 1][temperatureIndex]);
        double latentHeatAtNextTemperatureIndex = (1 - pressureRatio) * (state.latentHeatTable[pressureIndex][temperatureIndex + 1]) + pressureRatio * (state.latentHeatTable[pressureIndex + 1][temperatureIndex + 1]);
        latentHeat = (1 - temperatureRatio) * latentHeatAtTemperatureIndex + temperatureRatio * latentHeatAtNextTemperatureIndex;
    }
    return latentHeat;
}

namespace {

[[nodiscard]] SourceEnthalpy sourceEnthalpyOf(const ThermalState &state, int cellIndex) {
    Cel &cell = state.cells[cellIndex];
    Cel &leftCell = state.cells[cellIndex - 1];
    double sourceTemperature = cell.temp;
    double sourceGasEnthalpy;
    double sourceLiquidEnthalpy;
    double hcF = 0.;

    if (cell.acsr.tipo == kAccessoryGasInjection) {
        sourceTemperature = cell.acsr.injg.temp;
        sourceGasEnthalpy = cell.acsr.injg.FluidoPro.EntalpGas(cell.pres, sourceTemperature);
        sourceLiquidEnthalpy = 0;
        hcF = 0;
    } else if (cell.acsr.tipo == kAccessoryLiquidInjection) {
        sourceTemperature = cell.acsr.injl.temp;
        sourceGasEnthalpy = cell.acsr.injl.FluidoPro.EntalpGas(cell.pres, sourceTemperature);
        sourceLiquidEnthalpy = cell.acsr.injl.FluidoPro.EntalpLiq(cell.pres, sourceTemperature);
        hcF = cell.acsr.injl.fluidocol.CalorLiq(cell.pres, sourceTemperature) * sourceTemperature
            /*complementary-fluid enthalpy, to be corrected*/;
    } else if (cell.acsr.tipo == kAccessoryInflowPerformance) {
        sourceTemperature = cell.acsr.ipr.Tres;
        sourceGasEnthalpy = cell.acsr.ipr.FluidoPro.EntalpGas(cell.pres, sourceTemperature);
        sourceLiquidEnthalpy = cell.acsr.ipr.FluidoPro.EntalpLiq(cell.pres, sourceTemperature);
        hcF = 0;
    } else if (cell.acsrL != 0) {
        if ((*cell.acsrL).tipo == kAccessoryChoke) {
            if ((*cell.acsrL).chk.AreaGarg < state.input.master1.razareaativ * cell.dutoL.area && (*cell.acsrL).chk.AreaGarg > 1e-5 * cell.dutoL.area) {
                double chokeUpstreamTemperature = leftCell.temp;
                double chokeUpstreamVoidFraction = leftCell.alf;
                double betE = leftCell.bet;

                double upstreamLiquidDensity = leftCell.flui.MasEspLiq(leftCell.pres, leftCell.temp);
                double rholc = leftCell.fluicol.MasEspFlu(leftCell.pres, leftCell.temp);
                double upstreamMixtureLiquidDensity = (1 - betE) * upstreamLiquidDensity + betE * rholc;

                double chokeDownstreamVoidFraction = cell.alf;
                double betJ = cell.bet;
                double downstreamLiquidDensity = cell.flui.MasEspLiq(cell.pres, cell.temp);
                double rholcJ = cell.fluicol.MasEspFlu(cell.pres, cell.temp);
                double downstreamMixtureLiquidDensity = (1 - betJ) * downstreamLiquidDensity + betJ * rholcJ;

                const double upstreamGasDensity = leftCell.flui.MasEspGas(leftCell.pres, leftCell.temp);
                double upstreamHydrostaticHead = sin(leftCell.duto.teta) * (0.5 * leftCell.dx) * (upstreamMixtureLiquidDensity * (1 - chokeUpstreamVoidFraction) + chokeUpstreamVoidFraction * upstreamGasDensity) * kGravity / kPascalPerKgfPerCm2Coarse;
                double downstreamHydrostaticHead = sin(cell.duto.teta) * (0.5 * cell.dx) * (downstreamMixtureLiquidDensity * (1 - chokeDownstreamVoidFraction) + chokeDownstreamVoidFraction * cell.flui.MasEspGas(cell.pres, cell.temp)) * kGravity / kPascalPerKgfPerCm2Coarse;

                double quality = chokeUpstreamVoidFraction * upstreamGasDensity / (upstreamGasDensity * chokeUpstreamVoidFraction + upstreamMixtureLiquidDensity * (1. - chokeUpstreamVoidFraction));

                double upstreamLiquidJouleThomson = (1. - betE) * leftCell.flui.JTL(leftCell.pres - upstreamHydrostaticHead, leftCell.temp) - betE / rholcJ;
                double upstreamGasJouleThomson = leftCell.flui.JTG(leftCell.pres - upstreamHydrostaticHead, leftCell.temp);
                sourceTemperature = chokeUpstreamTemperature + ((1. - quality) * upstreamLiquidJouleThomson + quality * upstreamGasJouleThomson) * (cell.pres + downstreamHydrostaticHead - leftCell.pres - upstreamHydrostaticHead);

                sourceGasEnthalpy = leftCell.flui.EntalpGas(leftCell.pres, sourceTemperature);
                sourceLiquidEnthalpy = leftCell.flui.EntalpLiq(leftCell.pres, sourceTemperature);
                hcF = leftCell.fluicol.CalorLiq(leftCell.pres, sourceTemperature) * sourceTemperature;
                /*complementary-fluid enthalpy, to be corrected*/

            } else {
                sourceTemperature = cell.temp;
                sourceGasEnthalpy = leftCell.flui.EntalpGas(leftCell.pres, leftCell.temp);
                sourceLiquidEnthalpy = leftCell.flui.EntalpLiq(leftCell.pres, leftCell.temp);
                hcF = leftCell.fluicol.CalorLiq(leftCell.pres, leftCell.temp) * sourceTemperature;
                /*complementary-fluid enthalpy, to be corrected*/
            }
        } else if ((*cell.acsrL).tipo == kAccessoryVolumetricPump) {

            double polytropicExponent = (*cell.acsrL).bvol.npoli;
            double pumpPressureRatio = (cell.pres) / (leftCell.pres);
            sourceTemperature = leftCell.temp * pow(pumpPressureRatio, (polytropicExponent - 1) / polytropicExponent);

            sourceGasEnthalpy = leftCell.flui.EntalpGas(leftCell.pres, sourceTemperature);
            sourceLiquidEnthalpy = leftCell.flui.EntalpLiq(leftCell.pres, sourceTemperature);
            hcF = leftCell.fluicol.CalorLiq(leftCell.pres, sourceTemperature) * sourceTemperature;
            /*complementary-fluid enthalpy, to be corrected*/
        } else {
            sourceGasEnthalpy = 0.;
            sourceLiquidEnthalpy = 0.;
        }
    } else {
        sourceGasEnthalpy = 0.;
        sourceLiquidEnthalpy = 0.;
    }
    return SourceEnthalpy{.temperature = sourceTemperature,
                          .gasEnthalpy = sourceGasEnthalpy,
                          .liquidEnthalpy = sourceLiquidEnthalpy,
                          .hcF = hcF};
}

/// Whether the accessory forces a pressure discontinuity, so the downstream
/// pressure must be extrapolated rather than read.
[[nodiscard]] bool accessoryImposesPressureJump(const ThermalState &state,
                                                int cellIndex) {
    Cel &cell = state.cells[cellIndex];
    return (cell.acsr.tipo == kAccessoryChoke &&
            cell.acsr.chk.AreaGarg <=
                (1e-3 + state.input.master1.razareaativ) *
                    cell.duto.area) ||
           (cell.acsr.tipo == kAccessoryPump &&
            cell.acsr.bcs.freq > 0) ||
           (cell.acsr.tipo == kAccessoryVolumetricPump &&
            cell.acsr.bvol.freq > 0.) ||
           (cell.acsr.tipo == 7 &&
            fabs(cell.acsr.delp) > 0.) ||
           (cell.acsr.tipo == kAccessoryMultiPump &&
            cell.acsr.multibcs.freq > 0);
}

}  // namespace

double computeMixtureEnthalpy(const ThermalState &state, int cellIndex) {
    Cel &cell = state.cells[cellIndex];
    Cel &leftCell = state.cells[cellIndex - 1];
    Cel &rightCell = state.cells[cellIndex + 1];

    double cellLength = cell.dx;
    double diameter = cell.duto.a;
    double flowArea = 0.25 * M_PI * diameter * diameter;
    double betmed = cell.bet;
    double previousMeanVoidFraction = cell.alfini;
    double betmed0 = cell.betini;
    double meanPressure = cell.pres;
    double previousMeanPressure = cell.presini;
    double meanTemperature = cell.temp;
    double leftSuperficialGasVelocity;
    double leftSuperficialLiquidVelocity;
    double rightSuperficialGasVelocity;
    double rightSuperficialLiquidVelocity;

    leftSuperficialGasVelocity = cell.QG / flowArea;
    leftSuperficialLiquidVelocity = cell.QL / flowArea;
    rightSuperficialGasVelocity = rightCell.QG / flowArea;
    rightSuperficialLiquidVelocity = rightCell.QL / flowArea;

    double betL = cell.bet;
    double betR = cell.bet;
    if (leftSuperficialGasVelocity > 0.)
        betL = leftCell.bet;
    if (leftSuperficialGasVelocity < 0.)
        betL = rightCell.bet;

    double meanSuperficialGasVelocity = 0.5 * (leftSuperficialGasVelocity + rightSuperficialGasVelocity);
    double meanSuperficialLiquidVelocity = 0.5 * (leftSuperficialLiquidVelocity + rightSuperficialLiquidVelocity);

    double leftPressure = cell.presaux;
    double rightPressure = rightCell.presaux;
    if (accessoryImposesPressureJump(state, cellIndex)) {
        rightPressure = cell.pres + (cell.pres - cell.presaux) * 0.5;
    }

    double leftTemperature = cell.tempL;
    if (cell.VTemper < 0.)
        leftTemperature = cell.temp;
    double rightTemperature = cell.temp;
    if (rightCell.VTemper < 0.)
        rightTemperature = cell.tempR;

    double gasDensity = cell.flui.MasEspGas(previousMeanPressure, meanTemperature);
    double liquidDensity = cell.flui.MasEspLiq(previousMeanPressure, meanTemperature);
    double rhoc = cell.fluicol.MasEspFlu(previousMeanPressure, meanTemperature);
    double gasDensityAtCellPressure = cell.flui.MasEspGas(meanPressure, meanTemperature);
    double liquidDensityAtCellPressure = cell.flui.MasEspLiq(meanPressure, meanTemperature);
    double gasEnthalpy = cell.flui.EntalpGas(previousMeanPressure, meanTemperature);
    double liquidEnthalpy = cell.flui.EntalpLiq(previousMeanPressure, meanTemperature);
    double hc = cell.fluicol.CalorLiq(previousMeanPressure, meanTemperature) * meanTemperature; // enthalpy to be corrected

    double leftGasDensity = cell.flui.MasEspGas(leftPressure, leftTemperature);
    double leftLiquidDensity = cell.flui.MasEspLiq(leftPressure, leftTemperature);
    double rhocL = cell.fluicol.MasEspFlu(leftPressure, leftTemperature);
    double leftGasEnthalpy = cell.flui.EntalpGas(leftPressure, leftTemperature);
    double leftLiquidEnthalpy = cell.flui.EntalpLiq(leftPressure, leftTemperature);
    double hcL = cell.fluicol.CalorLiq(leftPressure, leftTemperature) * leftTemperature; // enthalpy to be corrected

    double rightGasDensity = cell.flui.MasEspGas(rightPressure, rightTemperature);
    double rightLiquidDensity = cell.flui.MasEspLiq(rightPressure, rightTemperature);
    double rhocR = cell.fluicol.MasEspFlu(rightPressure, rightTemperature);
    double rightGasEnthalpy = cell.flui.EntalpGas(rightPressure, rightTemperature);
    double rightLiquidEnthalpy = cell.flui.EntalpLiq(rightPressure, rightTemperature);
    double hcR = cell.fluicol.CalorLiq(rightPressure, rightTemperature) * rightTemperature; // enthalpy to be corrected

    double previousMixtureInternalEnergy = gasDensity * previousMeanVoidFraction * (gasEnthalpy - previousMeanPressure * kPascalPerKgfPerCm2 / gasDensity) + liquidDensity * (1 - betmed0) * (1. - previousMeanVoidFraction) * (liquidEnthalpy - previousMeanPressure * kPascalPerKgfPerCm2 / liquidDensity) +
                           betmed0 * rhoc * (1. - previousMeanVoidFraction) * (hc - previousMeanPressure * kPascalPerKgfPerCm2 / rhoc);
    double leftEnthalpyFlux = leftGasDensity * leftSuperficialGasVelocity * leftGasEnthalpy + leftLiquidDensity * (1. - betL) * leftSuperficialLiquidVelocity * leftLiquidEnthalpy + rhocL * betL * leftSuperficialLiquidVelocity * hcL;
    double rightEnthalpyFlux = rightGasDensity * rightSuperficialGasVelocity * rightGasEnthalpy + rightLiquidDensity * (1. - betR) * rightSuperficialLiquidVelocity * rightLiquidEnthalpy + rhocR * betR * rightSuperficialLiquidVelocity * hcR;
    double enthalpyFluxDivergence = (rightEnthalpyFlux - leftEnthalpyFlux) / cellLength;
    double hydrostaticTerm = (gasDensityAtCellPressure * meanSuperficialGasVelocity + (1 - betmed) * meanSuperficialLiquidVelocity * liquidDensityAtCellPressure + betmed * meanSuperficialLiquidVelocity) * sin(cell.duto.teta) * kGravity;

    double gasMassSourceTerm = 0.;
    double liquidMassSourceTerm = 0.;
    const SourceEnthalpy source = sourceEnthalpyOf(state, cellIndex);
    // source.temperature is unused here: this function never reads it.
    const double sourceGasEnthalpy = source.gasEnthalpy;
    const double sourceLiquidEnthalpy = source.liquidEnthalpy;
    const double hcF = source.hcF;

    liquidMassSourceTerm = 0;
    gasMassSourceTerm = 0;
    if (cell.fontemassLR > 0.)
        liquidMassSourceTerm = sourceLiquidEnthalpy * cell.fontemassLR / cellLength;
    if (cell.fontemassCR > 0.)
        liquidMassSourceTerm += hcF * cell.fontemassCR / cellLength;
    if (gasMassSourceTerm > 0.)
        gasMassSourceTerm = sourceGasEnthalpy * cell.fontemassGR / cellLength;

    return previousMixtureInternalEnergy - (enthalpyFluxDivergence + hydrostaticTerm - (gasMassSourceTerm + liquidMassSourceTerm) / flowArea) * cell.dt;
}

double interpolateMixtureEnergy(const ThermalState &state, int cellIndex, int pressureIndex, int temperatureIndex, double pressureRatio) {
    Cel &cell = state.cells[cellIndex];
    int nextPressureIndex = pressureIndex + 1;

    double lowerPressure = cell.flui.rhogF[pressureIndex][0];
    double upperPressure = cell.flui.rhogF[nextPressureIndex][0];
    double temperature = cell.flui.rhogF[0][temperatureIndex];
    double meanVoidFraction = cell.alf;
    double betmed = cell.bet;

    double lowerGasDensity = cell.flui.rhogF[pressureIndex][temperatureIndex];
    double upperGasDensity = cell.flui.rhogF[nextPressureIndex][temperatureIndex];
    double lowerGasEnthalpy = cell.flui.HgF[pressureIndex][temperatureIndex];
    double upperGasEnthalpy = cell.flui.HgF[nextPressureIndex][temperatureIndex];

    double lowerLiquidDensity = cell.flui.rholF[pressureIndex][temperatureIndex];
    double upperLiquidDensity = cell.flui.rholF[nextPressureIndex][temperatureIndex];
    double lowerLiquidEnthalpy = cell.flui.HlF[pressureIndex][temperatureIndex];
    double upperLiquidEnthalpy = cell.flui.HlF[nextPressureIndex][temperatureIndex];

    double rhocp0 = cell.fluicol.MasEspFlu(lowerPressure, temperature);
    double rhocp1 = cell.fluicol.MasEspFlu(upperPressure, temperature);
    double hlc0 = cell.fluicol.CalorLiq(lowerPressure, temperature) * temperature; // enthalpy to be corrected
    double hlc1 = cell.fluicol.CalorLiq(upperPressure, temperature) * temperature; // enthalpy to be corrected

    double lowerMixtureEnergy = meanVoidFraction * lowerGasDensity * (lowerGasEnthalpy - lowerPressure * kPascalPerKgfPerCm2 / lowerGasDensity) +
                    (1 - meanVoidFraction) * (1 - betmed) * lowerLiquidDensity * (lowerLiquidEnthalpy - lowerPressure * kPascalPerKgfPerCm2 / lowerLiquidDensity) +
                    (1 - meanVoidFraction) * (betmed)*rhocp0 * (hlc0 - lowerPressure * kPascalPerKgfPerCm2 / rhocp0);

    double upperMixtureEnergy = meanVoidFraction * upperGasDensity * (upperGasEnthalpy - upperPressure * kPascalPerKgfPerCm2 / lowerGasDensity) +
                    (1 - meanVoidFraction) * (1 - betmed) * upperLiquidDensity * (upperLiquidEnthalpy - upperPressure * kPascalPerKgfPerCm2 / upperLiquidDensity) +
                    (1 - meanVoidFraction) * (betmed)*rhocp1 * (hlc1 - upperPressure * kPascalPerKgfPerCm2 / rhocp1);

    return pressureRatio * lowerMixtureEnergy + (1. - pressureRatio) * upperMixtureEnergy;
}

void updateTemperatureFromEnthalpy(const ThermalState &state, int cellIndex) {
    Cel &cell = state.cells[cellIndex];
    double pressure = cell.pres;
    double **propertyTable = cell.flui.rholF;

    int pressureIndex = 0;
    int divisionCount = cell.flui.npontos - 1;
    if (pressure < propertyTable[1][0] || pressure >= propertyTable[divisionCount + 1][0]) {
        cout << "pressure outside the table bounds";
    }

    pressureIndex = bracketingIndex(
        pressure, divisionCount,
        [&](int index) { return propertyTable[index][0]; }, 0);

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
    cell.temp = propertyTable[0][temperatureIndex] * temperatureRatio + (1. - temperatureRatio) * propertyTable[0][temperatureIndex + 1];

    if (cell.temp < kMinimumTemperatureCelsius)
        cell.temp = kMinimumTemperatureCelsius;
    if (cell.temp > kMaximumTemperatureCelsius)
        cell.temp = kMaximumTemperatureCelsius;
}

namespace {

[[nodiscard]] CellFlowBasis cellFlowBasisOf(const ThermalState &state, int cellIndex) {
    varGlob1D &globals = *state.globals;
    Cel &cell = state.cells[cellIndex];
    Cel &leftCell = state.cells[cellIndex - 1];
    Cel &rightCell = state.cells[cellIndex + 1];
    const bool cellHasGas = cell.alf > globals.localtiny;
    const bool cellHasLiquid = cell.alf < 1. - globals.localtiny;
    double diameter = cell.duto.a;
    double flowArea = 0.25 * M_PI * diameter * diameter;
    double meanVoidFraction = cell.alf;
    double betmed = cell.bet;
    double meanSuperficialGasVelocity;
    double meanSuperficialLiquidVelocity;
    if (cellIndex > 0 && (leftCell.acsr.tipo != kAccessoryChoke ||
                          leftCell.acsr.chk.AreaGarg > (1e-3 + state.input.master1.razareaativ) * leftCell.duto.area)) {
        if (cellHasGas)
            meanSuperficialGasVelocity = cell.QG / flowArea;
        else {
            meanSuperficialGasVelocity = 0.;
        }
        if (cellHasLiquid)
            meanSuperficialLiquidVelocity = cell.QL / flowArea;
        else {
            meanSuperficialLiquidVelocity = 0.;
        }
    } else {
        if (cellHasGas)
            meanSuperficialGasVelocity = rightCell.QG / flowArea;
        else {
            meanSuperficialGasVelocity = 0.;
        }
        if (cellHasLiquid)
            meanSuperficialLiquidVelocity = rightCell.QL / flowArea;
        else {
            meanSuperficialLiquidVelocity = 0.;
        }
    }
    return CellFlowBasis{.flowArea = flowArea,
                         .voidFraction = meanVoidFraction,
                         .betmed = betmed,
                         .gasSuperficialVelocity = meanSuperficialGasVelocity,
                         .liquidSuperficialVelocity = meanSuperficialLiquidVelocity};
}

TemperatureBalance prepareTemperatureBalance(const ThermalState &state,
                                             int cellIndex,
                                             int steadyStateMode) {
    Cel &cell = state.cells[cellIndex];
    Cel &leftCell = state.cells[cellIndex - 1];
    double cellLength = cell.dx;
    double meanCellLength;
    if (cellIndex > 0)
        meanCellLength = 0.5 * (cell.dx + leftCell.dx);
    else
        meanCellLength = 0.5 * cell.dx;
    const CellFlowBasis basis = cellFlowBasisOf(state, cellIndex);
    const double flowArea = basis.flowArea;
    const double meanVoidFraction = basis.voidFraction;
    const double betmed = basis.betmed;
    const double meanSuperficialGasVelocity = basis.gasSuperficialVelocity;
    const double meanSuperficialLiquidVelocity = basis.liquidSuperficialVelocity;
    double referenceMixtureVelocity = fabs(meanSuperficialGasVelocity) + fabs(meanSuperficialLiquidVelocity);
    double producedLiquidDensity = cell.rpC;
    double completionFluidDensity = cell.rcC;
    double liquidDensity = (1. - betmed) * producedLiquidDensity + betmed * completionFluidDensity;
    double gasDensity = cell.rgC;
    double liquidSpecificHeat = (1. - betmed) * cell.flui.CalorLiq(cell.presini, cell.temp) + betmed * cell.fluicol.CalorLiq(cell.presini, cell.temp);
    double liquidSpecificHeatConstantVolume = liquidSpecificHeat;
    double gasSpecificHeat = cell.flui.CalorGas(cell.presini, cell.temp);
    double gasSpecificHeatConstantVolume = cell.flui.CalorGasVolMod(cell.presini, cell.temp, cell.rgC);
    double liquidJouleThomson = (1. - betmed) * cell.flui.JTL(cell.presini, cell.temp) - betmed / completionFluidDensity;
    double gasJouleThomson = cell.flui.JTG(cell.presini, cell.temp, cell.rgC);
    double hydrostaticPower = (liquidDensity * meanSuperficialLiquidVelocity + gasDensity * meanSuperficialGasVelocity) * flowArea * kGravity * sin(cell.duto.teta);

    double heatFlux = 0.;
    cell.calor.Tint = cell.temp;
    if (cellIndex > 0)
        cell.calor.dtL = cell.temp - leftCell.tempini;
    else
        cell.calor.dtL = 0;
    cell.calor.Vint = meanSuperficialGasVelocity + meanSuperficialLiquidVelocity;
    cell.calor.dt = cell.dt;
    double liquidConductivity = (1. - betmed) * cell.flui.CondLiq(cell.presini, cell.temp) + betmed * cell.fluicol.CondLiq(cell.presini, cell.temp);
    cell.calor.kint = liquidConductivity * (1 - meanVoidFraction) + cell.flui.CondGas(cell.presini, cell.temp) * meanVoidFraction;
    cell.calor.cpint = liquidSpecificHeat * (1 - meanVoidFraction) + gasSpecificHeat * meanVoidFraction;
    cell.calor.rhoint = liquidDensity * (1 - meanVoidFraction) + gasDensity * meanVoidFraction;
    double liquidViscosity = (1. - betmed) * cell.mipC + betmed * cell.micC;
    cell.calor.viscint = liquidViscosity * (1 - meanVoidFraction) * kPascalSecondPerCentipoise + cell.migC * meanVoidFraction * kPascalSecondPerCentipoise;
    double temperaturePerturbation = cell.temp * 0.01;
    if (fabs(cell.temp) < 1e-15)
        temperaturePerturbation = 0.1;
    double liquidDensityChange = (1. - betmed) * cell.flui.MasEspLiq(cell.presini, cell.temp + temperaturePerturbation) +
                    betmed * cell.fluicol.MasEspFlu(cell.presini, cell.temp + temperaturePerturbation) - liquidDensity;
    double gasDensityChange = cell.flui.MasEspGas(cell.presini, cell.temp + temperaturePerturbation) - gasDensity;
    cell.calor.betint = -(1 / cell.calor.rhoint) * (liquidDensityChange * (1 - meanVoidFraction) + gasDensityChange * meanVoidFraction) / (temperaturePerturbation);
    if (state.input.modoDifus3D == 0 || steadyStateMode != 0) {
        if (steadyStateMode == 0)
            heatFlux = cell.calor.transtrans();
        else
            heatFlux = cell.calor.transperm();
        if (state.productionNetworkCoupled == 1 && cellIndex >= state.primarySectionStart && cellIndex <= state.primarySectionEnd) {
            heatFlux -= cell.fluxcalAcopRedeP;
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
            heatFlux = cell.calor.transtrans();
        else {
            int iacop1 = state.coupledCellIndices[coupled];
            heatFlux = -state.input.celAcop[coupled].FE * state.poissonSolver3D.dados.qTotal[iacop1] / cell.dx;
        }
    }

    cell.fluxcalmed = heatFlux;

    double timeCoefficient = (liquidDensity * (1 - meanVoidFraction) * liquidSpecificHeatConstantVolume + gasDensity * meanVoidFraction * gasSpecificHeatConstantVolume) * flowArea;
    double pressureTimeCoefficient;
    pressureTimeCoefficient = -cell.flui.CalorGasPresMod(cell.presini, cell.temp) * (gasDensity * meanVoidFraction * flowArea); //-cell.flui.CalorGasPresMod(cell.pres, cell.temp) * (gasDensity * meanVoidFraction * flowArea);
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
        .rp = producedLiquidDensity,
        .rc = completionFluidDensity,
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
    Cel &cell = state.cells[cellIndex];
    Cel &leftCell = state.cells[cellIndex - 1];
    Cel &rightCell = state.cells[cellIndex + 1];
    double kineticTerm = 0;
    double leftFaceMeanGasVelocity = 0;
    double leftFaceMeanLiquidVelocity = 0;
    double rightFaceMeanGasVelocity = 0;
    double rightFaceMeanLiquidVelocity = 0;
    double previousRightFaceMeanGasVelocity = 0;
    double previousRightFaceMeanLiquidVelocity = 0;
    if (cellIndex <= state.lastCell - 1 && cell.acsr.tipo == kAccessoryNone && rightCell.acsr.tipo == kAccessoryNone) {

        double kineticCellLength = cell.dx;
        double cellDiameter = cell.duto.a;
        double cellFlowArea = 0.25 * M_PI * cellDiameter * cellDiameter;

        double rpcin = rightCell.rpCi;
        double rccin = rightCell.rcCi;
        double kineticLiquidDensity = (1. - balance.bet) * rpcin + balance.bet * rccin;
        double kineticGasDensity = rightCell.rgCi;

        double rightFaceVoidFraction;
        double previousRightFaceVoidFraction;

        balance.gasSuperficialVelocity = rightCell.QG / balance.flowArea;
        balance.liquidSuperficialVelocity = rightCell.QL / balance.flowArea;

        if (rightCell.QG > 0)
            rightFaceVoidFraction = cell.alf;
        else
            rightFaceVoidFraction = rightCell.alf;

        if (rightCell.QGini > 0)
            previousRightFaceVoidFraction = cell.alfini;
        else
            previousRightFaceVoidFraction = rightCell.alfini;

        double leftFaceVoidFraction;
        if (cell.QG > 0)
            leftFaceVoidFraction = leftCell.alf;
        else
            leftFaceVoidFraction = cell.alf;

        if (leftFaceVoidFraction > 1e-3) {
            leftFaceMeanGasVelocity = cell.QG / (cellFlowArea);
            leftFaceMeanGasVelocity /= leftFaceVoidFraction;
        }

        if (leftFaceVoidFraction < 1. - 1e-3) {
            leftFaceMeanLiquidVelocity = cell.QL / (cellFlowArea);
            leftFaceMeanLiquidVelocity /= (1. - leftFaceVoidFraction);
        }

        if (rightFaceVoidFraction > 1e-3) {
            rightFaceMeanGasVelocity = balance.gasSuperficialVelocity;
            rightFaceMeanGasVelocity /= rightFaceVoidFraction;
        }

        if (rightFaceVoidFraction < 1. - 1e-3) {
            rightFaceMeanLiquidVelocity = balance.liquidSuperficialVelocity;
            rightFaceMeanLiquidVelocity /= (1. - rightFaceVoidFraction);
        }

        if (previousRightFaceVoidFraction > 1e-3) {

            previousRightFaceMeanGasVelocity = rightCell.QGini / (cellFlowArea);
            previousRightFaceMeanGasVelocity /= previousRightFaceVoidFraction;
        }

        if (previousRightFaceVoidFraction < 1. - 1e-3) {

            previousRightFaceMeanLiquidVelocity = rightCell.QLini / (cellFlowArea);
            previousRightFaceMeanLiquidVelocity /= (1. - previousRightFaceVoidFraction);
        }

        kineticTerm = kineticLiquidDensity * (1 - rightFaceVoidFraction) * cellFlowArea * (0.5 * (rightFaceMeanLiquidVelocity * rightFaceMeanLiquidVelocity - previousRightFaceMeanLiquidVelocity * previousRightFaceMeanLiquidVelocity)) / cell.dt +
                   kineticGasDensity * rightFaceVoidFraction * cellFlowArea * (0.5 * (rightFaceMeanGasVelocity * rightFaceMeanGasVelocity - previousRightFaceMeanGasVelocity * previousRightFaceMeanGasVelocity)) / cell.dt +
                   ((rightCell.MC - rightCell.Mliqini) * rightFaceMeanGasVelocity * (rightFaceMeanGasVelocity - leftFaceMeanGasVelocity) / kineticCellLength +
                    rightCell.Mliqini * rightFaceMeanLiquidVelocity * (rightFaceMeanLiquidVelocity - leftFaceMeanLiquidVelocity) / kineticCellLength);
    }
    return kineticTerm;
}

/// Resets the source specific heats to the inert defaults: no gas or liquid
/// heat capacity, unit adiabatic ratio.
void clearSourceSpecificHeats(double &sourceGasSpecificHeat,
                              double &sourceSpecificHeatRatio,
                              double &sourceLiquidSpecificHeat) {
    sourceGasSpecificHeat = 0.;
    sourceSpecificHeatRatio = 1.;
    sourceLiquidSpecificHeat = 0.;
}

TemperatureSourceTerms computeTemperatureSourceTerms(const ThermalState &state,
                                                      int cellIndex,
                                                      double cellLength) {
    Cel &cell = state.cells[cellIndex];
    Cel &leftCell = state.cells[cellIndex - 1];
    double gasMassSourceTerm = 0.;
    double liquidMassSourceTerm = 0.;
    double sourceTemperature = cell.temp;
    double sourceGasSpecificHeat;
    double sourceSpecificHeatRatio = 0.;
    double sourceLiquidSpecificHeat;
    if (cell.acsr.tipo == kAccessoryGasInjection) {
        sourceTemperature = cell.acsr.injg.temp;
        sourceGasSpecificHeat = cell.acsr.injg.FluidoPro.CalorGas(cell.presini, cell.temp);
        sourceSpecificHeatRatio = cell.acsr.injg.FluidoPro.ConstAdG(cell.presini, cell.temp);
        sourceLiquidSpecificHeat = 0.;
    } else if (cell.acsr.tipo == kAccessoryLiquidInjection) {
        sourceTemperature = cell.acsr.injl.temp;
        sourceGasSpecificHeat = cell.acsr.injl.FluidoPro.CalorGas(cell.presini, cell.temp);                                                                                                                                        // cell.acsr.injl.FluidoPro.CalorGas(cell.pres, cell.temp);
        sourceSpecificHeatRatio = cell.acsr.injl.FluidoPro.CalorGas(cell.presini, cell.temp);                                                                                                                                      // cell.acsr.injl.FluidoPro.ConstAdG(cell.pres, cell.temp);
        sourceLiquidSpecificHeat = (1. - cell.acsr.injl.bet) * cell.acsr.injl.FluidoPro.CalorLiq(cell.presini, cell.temp) + cell.acsr.injl.bet * cell.acsr.injl.fluidocol.CalorLiq(cell.presini, cell.temp); //(1. - cell.acsr.injl.bet) * cell.acsr.injl.FluidoPro.CalorLiq(cell.pres, cell.temp)
    } else if (cell.acsr.tipo == kAccessoryMultipleSource) {
        sourceTemperature = cell.acsr.injm.temp;
        sourceGasSpecificHeat = cell.acsr.injm.FluidoPro.CalorGas(cell.presini, cell.temp);
        sourceSpecificHeatRatio = cell.acsr.injm.FluidoPro.ConstAdG(cell.presini, cell.temp);
        if ((cell.fontemassLR + cell.fontemassCR) > 0.) {
            double titbet = cell.fontemassCR / (cell.fontemassLR + cell.fontemassCR);
            sourceLiquidSpecificHeat = (1. - titbet) * cell.acsr.injl.FluidoPro.CalorLiq(cell.presini, cell.temp) +
                   titbet * cell.acsr.injl.fluidocol.CalorLiq(cell.presini, cell.temp);
        } else
            sourceLiquidSpecificHeat = 0.;
    } else if (cell.acsr.tipo == kAccessoryInflowPerformance) {
        sourceTemperature = cell.acsr.ipr.Tres;
        sourceGasSpecificHeat = cell.acsr.ipr.FluidoPro.CalorGas(cell.presini, cell.temp);
        sourceSpecificHeatRatio = cell.acsr.ipr.FluidoPro.ConstAdG(cell.presini, cell.temp);
        sourceLiquidSpecificHeat = cell.acsr.ipr.FluidoPro.CalorLiq(cell.presini, cell.temp);
    } else if (cell.acsr.tipo == kAccessoryRadialPorous) {
        sourceTemperature = cell.acsr.radialPoro.tRes;
        sourceGasSpecificHeat = cell.acsr.radialPoro.flup.CalorGas(cell.presini, cell.temp);
        sourceSpecificHeatRatio = cell.acsr.radialPoro.flup.ConstAdG(cell.presini, cell.temp);
        sourceLiquidSpecificHeat = cell.acsr.radialPoro.flup.CalorLiq(cell.presini, cell.temp);
    } else if (cell.acsr.tipo == kAccessoryLeak) {
        double betM = cell.acsr.fontechk.betISamb;
        double ambientPressure = cell.acsr.fontechk.pamb;
        double ambientTemperature = cell.acsr.fontechk.tamb;
        sourceGasSpecificHeat = cell.acsr.fontechk.fluidoPamb.CalorGas(ambientPressure, ambientTemperature);
        sourceLiquidSpecificHeat = (1. - betM) * cell.acsr.fontechk.fluidoPamb.CalorLiq(ambientPressure, ambientTemperature) + betM * cell.acsr.fontechk.fluidocol.CalorLiq(ambientPressure, ambientTemperature);
        sourceSpecificHeatRatio = cell.acsr.fontechk.fluidoPamb.ConstAdG(ambientPressure, ambientTemperature);
        sourceTemperature = cell.acsr.fontechk.tamb;
    } else if (cell.acsrL != 0 && cellIndex < state.lastCell) {
        if ((*cell.acsrL).tipo == kAccessoryChoke && (cell.fontemassLR + cell.fontemassCR + cell.fontemassGR) > 0.) {
            if ((*cell.acsrL).chk.AreaGarg < state.input.master1.razareaativ * cell.dutoL.area && (*cell.acsrL).chk.AreaGarg > 1e-5 * cell.dutoL.area) {
                double upstreamTemperature = leftCell.tempini;
                double chokeUpstreamVoidFraction = leftCell.alf;
                double betE = leftCell.bet;

                double upstreamLiquidDensity = leftCell.rpC;
                double rholc = leftCell.rcC;
                double upstreamMixtureLiquidDensity = (1 - betE) * upstreamLiquidDensity + betE * rholc;

                double chokeDownstreamVoidFraction = cell.alf;
                double betJ = cell.bet;
                double downstreamLiquidDensity = cell.rpC;
                double rholcJ = cell.rcC;
                double downstreamMixtureLiquidDensity = (1 - betJ) * downstreamLiquidDensity + betJ * rholcJ;

                double upstreamHydrostaticHead = sin(leftCell.duto.teta) * (0.5 * leftCell.dx) * (upstreamMixtureLiquidDensity * (1 - chokeUpstreamVoidFraction) + chokeUpstreamVoidFraction * leftCell.rgC) * kGravity / kPascalPerKgfPerCm2Coarse;
                double downstreamHydrostaticHead = sin(cell.duto.teta) * (0.5 * cell.dx) * (downstreamMixtureLiquidDensity * (1 - chokeDownstreamVoidFraction) + chokeDownstreamVoidFraction * cell.rgC) * kGravity / kPascalPerKgfPerCm2Coarse;
                double quality = chokeUpstreamVoidFraction * leftCell.rgC / (leftCell.rgC * chokeUpstreamVoidFraction + upstreamMixtureLiquidDensity * (1. - chokeUpstreamVoidFraction));
                sourceGasSpecificHeat = leftCell.flui.CalorGas(leftCell.presini, leftCell.tempini);                                                                                                                              // leftCell.flui.CalorGas(leftCell.pres, leftCell.temp);
                sourceLiquidSpecificHeat = (1 - leftCell.bet) * leftCell.flui.CalorLiq(leftCell.presini, leftCell.tempini) + leftCell.bet * leftCell.fluicol.CalorLiq(leftCell.presini, leftCell.tempini); //(1 - leftCell.bet) * leftCell.flui.CalorLiq(leftCell.pres, leftCell.temp)
                double upstreamLiquidJouleThomson = (1. - betE) * leftCell.flui.JTL(leftCell.presini - upstreamHydrostaticHead, leftCell.tempini) - betE / rholcJ;                                                                                     //(1. - betE) * leftCell.flui.JTL(leftCell.pres - upstreamHydrostaticHead, leftCell.temp)
                double upstreamGasJouleThomson = leftCell.flui.JTG(leftCell.presini - upstreamHydrostaticHead, leftCell.tempini);                                                                                                                   // leftCell.flui.JTG(leftCell.pres - upstreamHydrostaticHead, leftCell.temp);
                sourceTemperature = upstreamTemperature + ((1. - quality) * upstreamLiquidJouleThomson / sourceLiquidSpecificHeat + quality * upstreamGasJouleThomson / sourceGasSpecificHeat) *
                                  (cell.pres + downstreamHydrostaticHead - leftCell.pres - upstreamHydrostaticHead) * kPascalPerKgfPerCm2Variant;

                sourceSpecificHeatRatio = leftCell.flui.ConstAdG(leftCell.presini, leftCell.tempini);

            } else {
                sourceGasSpecificHeat = leftCell.flui.CalorGas(leftCell.presini, leftCell.tempini);
                sourceLiquidSpecificHeat = (1 - leftCell.bet) * leftCell.flui.CalorLiq(leftCell.presini, leftCell.tempini) + leftCell.bet * leftCell.fluicol.CalorLiq(leftCell.presini, leftCell.tempini);
                sourceTemperature = cell.temp;
                sourceSpecificHeatRatio = leftCell.flui.ConstAdG(leftCell.presini, leftCell.tempini);
            }
        } else if ((*cell.acsrL).tipo == kAccessoryVolumetricPump) {
            double betM = leftCell.bet;
            sourceGasSpecificHeat = leftCell.flui.CalorGas(leftCell.presini, leftCell.tempini);
            sourceLiquidSpecificHeat = (1. - betM) * leftCell.flui.CalorLiq(leftCell.presini, leftCell.tempini) + betM * leftCell.fluicol.CalorLiq(leftCell.presini, leftCell.tempini);
            sourceSpecificHeatRatio = leftCell.flui.ConstAdG(leftCell.presini, leftCell.tempini);
            double polytropicExponent = (*cell.acsrL).bvol.npoli;
            double pumpPressureRatio = (cell.pres) / (leftCell.pres);
            sourceTemperature = leftCell.tempini * pow(pumpPressureRatio, (polytropicExponent - 1) / polytropicExponent);
        } else {
            clearSourceSpecificHeats(sourceGasSpecificHeat, sourceSpecificHeatRatio,
                                     sourceLiquidSpecificHeat);
        }
    } else if (cellIndex == state.lastCell && (cell.fontemassLR + cell.fontemassCR + cell.fontemassGR) > 0.) {
        if ((*state.globals).chaverede == 0 || state.endNode == 1 || (*state.globals).chaveRedeParalela == 1)
            sourceTemperature = cell.calor.Textern1;
        else
            sourceTemperature = state.gasSurfaceTemperature;
        sourceGasSpecificHeat = cell.flui.CalorGas(cell.presini, sourceTemperature);
        double betM = 0.;
        if (fabs(cell.fontemassCR) > 1e-15)
            betM = cell.fontemassCR /
                   (cell.fontemassCR + cell.fontemassLR);
        sourceLiquidSpecificHeat = (1. - betM) * cell.flui.CalorLiq(cell.presini, sourceTemperature) + betM * cell.fluicol.CalorLiq(cell.presini, sourceTemperature);
        sourceSpecificHeatRatio = cell.flui.ConstAdG(cell.pres, sourceTemperature);
    } else {
        clearSourceSpecificHeats(sourceGasSpecificHeat, sourceSpecificHeatRatio,
                                 sourceLiquidSpecificHeat);
    }

    liquidMassSourceTerm = 0;
    if (cell.fontemassLR > 0.)
        liquidMassSourceTerm = cell.fontemassLR / cellLength;
    if (cell.fontemassCR > 0.)
        liquidMassSourceTerm += cell.fontemassCR / cellLength;
    liquidMassSourceTerm *= (sourceLiquidSpecificHeat) * (sourceTemperature - cell.temp);

    gasMassSourceTerm = cell.fontemassGR / cellLength;
    if (gasMassSourceTerm > 0.)
        gasMassSourceTerm *= (sourceGasSpecificHeat / sourceSpecificHeatRatio) * (sourceTemperature - cell.temp);
    else
        gasMassSourceTerm = 0;

    return TemperatureSourceTerms{
        .gas = gasMassSourceTerm,
        .liquid = liquidMassSourceTerm,
    };
}

TemperatureSourceTerms computeThermalMassTransferSourceTerms(
    const ThermalState &state, int cellIndex, double cellLength) {
    Cel &cell = state.cells[cellIndex];
    Cel &leftCell = state.cells[cellIndex - 1];
    double gasMassSourceTerm = 0.;
    double liquidMassSourceTerm = 0.;
    double sourceTemperature = cell.temp;
    double sourceGasSpecificHeat;
    double sourceSpecificHeatRatio = 0.;
    double sourceLiquidSpecificHeat;
    if (cell.acsr.tipo == kAccessoryGasInjection) {
        sourceTemperature = cell.acsr.injg.temp;
        sourceGasSpecificHeat = cell.acsr.injg.FluidoPro.CalorGas(cell.pres, cell.temp);
        sourceSpecificHeatRatio = cell.acsr.injg.FluidoPro.ConstAdG(cell.pres, cell.temp);
        sourceLiquidSpecificHeat = 0.;
    } else if (cell.acsr.tipo == kAccessoryLiquidInjection) {
        sourceTemperature = cell.acsr.injl.temp;
        ;
        sourceGasSpecificHeat = cell.acsr.injl.FluidoPro.CalorGas(cell.pres, cell.temp);
        sourceSpecificHeatRatio = cell.acsr.injl.FluidoPro.ConstAdG(cell.pres, cell.temp);
        sourceLiquidSpecificHeat = (1. - cell.acsr.injl.bet) * cell.acsr.injl.FluidoPro.CalorLiq(cell.pres, cell.temp) + cell.acsr.injl.bet * cell.acsr.injl.fluidocol.CalorLiq(cell.pres, cell.temp);
    } else if (cell.acsr.tipo == kAccessoryInflowPerformance) {
        sourceTemperature = cell.acsr.ipr.Tres;
        sourceGasSpecificHeat = cell.acsr.ipr.FluidoPro.CalorGas(cell.pres, cell.temp);
        sourceSpecificHeatRatio = cell.acsr.ipr.FluidoPro.ConstAdG(cell.pres, cell.temp);
        sourceLiquidSpecificHeat = cell.acsr.ipr.FluidoPro.CalorLiq(cell.pres, cell.temp);
    } else if (cell.acsrL != 0) {
        if ((*cell.acsrL).tipo == kAccessoryChoke) {
            if ((*cell.acsrL).chk.AreaGarg < state.input.master1.razareaativ * cell.dutoL.area && (*cell.acsrL).chk.AreaGarg > 1e-5 * cell.dutoL.area) {
                double upstreamTemperature = leftCell.temp;
                double chokeUpstreamVoidFraction = leftCell.alf;
                double betE = leftCell.bet;

                double upstreamLiquidDensity = leftCell.flui.MasEspLiq(leftCell.pres, leftCell.temp);
                double rholc = leftCell.fluicol.MasEspFlu(leftCell.pres, leftCell.temp);
                double upstreamMixtureLiquidDensity = (1 - betE) * upstreamLiquidDensity + betE * rholc;

                double chokeDownstreamVoidFraction = cell.alf;
                double betJ = cell.bet;
                double downstreamLiquidDensity = cell.flui.MasEspLiq(cell.pres, cell.temp);
                double rholcJ = cell.fluicol.MasEspFlu(cell.pres, cell.temp);
                double downstreamMixtureLiquidDensity = (1 - betJ) * downstreamLiquidDensity + betJ * rholcJ;

                double upstreamHydrostaticHead = sin(leftCell.duto.teta) * (0.5 * leftCell.dx) * (upstreamMixtureLiquidDensity * (1 - chokeUpstreamVoidFraction) + chokeUpstreamVoidFraction * leftCell.flui.MasEspGas(leftCell.pres, leftCell.temp)) * kGravity / kPascalPerKgfPerCm2Coarse;
                double downstreamHydrostaticHead = sin(cell.duto.teta) * (0.5 * cell.dx) * (downstreamMixtureLiquidDensity * (1 - chokeDownstreamVoidFraction) + chokeDownstreamVoidFraction * cell.flui.MasEspGas(cell.pres, cell.temp)) * kGravity / kPascalPerKgfPerCm2Coarse;

                double quality = chokeUpstreamVoidFraction * leftCell.flui.MasEspGas(leftCell.pres, leftCell.temp) / (leftCell.flui.MasEspGas(leftCell.pres, leftCell.temp) * chokeUpstreamVoidFraction + upstreamMixtureLiquidDensity * (1. - chokeUpstreamVoidFraction));

                double upstreamLiquidJouleThomson = (1. - betE) * leftCell.flui.JTL(leftCell.pres - upstreamHydrostaticHead, leftCell.temp) - betE / rholcJ;
                double upstreamGasJouleThomson = leftCell.flui.JTG(leftCell.pres - upstreamHydrostaticHead, leftCell.temp);
                sourceTemperature = upstreamTemperature + ((1. - quality) * upstreamLiquidJouleThomson + quality * upstreamGasJouleThomson) * (cell.pres + downstreamHydrostaticHead - leftCell.pres - upstreamHydrostaticHead);

                sourceGasSpecificHeat = leftCell.flui.CalorGas(leftCell.pres, leftCell.temp);
                sourceLiquidSpecificHeat = (1 - leftCell.bet) * leftCell.flui.CalorLiq(leftCell.pres, leftCell.temp) + leftCell.bet * leftCell.fluicol.CalorLiq(leftCell.pres, leftCell.temp);
                sourceSpecificHeatRatio = leftCell.flui.ConstAdG(leftCell.pres, leftCell.temp);

            } else {
                sourceGasSpecificHeat = leftCell.flui.CalorGas(leftCell.pres, leftCell.temp);
                sourceLiquidSpecificHeat = (1 - leftCell.bet) * leftCell.flui.CalorLiq(leftCell.pres, leftCell.temp) + leftCell.bet * leftCell.fluicol.CalorLiq(leftCell.pres, leftCell.temp);
                sourceTemperature = cell.temp;
                sourceSpecificHeatRatio = leftCell.flui.ConstAdG(leftCell.pres, leftCell.temp);
            }
        } else if ((*cell.acsrL).tipo == kAccessoryVolumetricPump) {
            [[maybe_unused]] double pumpUpstreamVoidFraction = leftCell.alf;
            double betM = leftCell.bet;
            sourceGasSpecificHeat = leftCell.flui.CalorGas(leftCell.pres, leftCell.temp);
            sourceLiquidSpecificHeat = (1. - betM) * leftCell.flui.CalorLiq(leftCell.pres, leftCell.temp) + betM * leftCell.fluicol.CalorLiq(leftCell.pres, leftCell.temp);
            double polytropicExponent = (*cell.acsrL).bvol.npoli;
            double pumpPressureRatio = (cell.pres) / (leftCell.pres);
            sourceTemperature = leftCell.temp * pow(pumpPressureRatio, (polytropicExponent - 1) / polytropicExponent);
        } else {
            clearSourceSpecificHeats(sourceGasSpecificHeat, sourceSpecificHeatRatio,
                                     sourceLiquidSpecificHeat);
        }
    } else {
        clearSourceSpecificHeats(sourceGasSpecificHeat, sourceSpecificHeatRatio,
                                 sourceLiquidSpecificHeat);
    }

    liquidMassSourceTerm = 0;
    if (cell.fontemassLR > 0.)
        liquidMassSourceTerm = cell.fontemassLR / cellLength;
    if (cell.fontemassCR > 0.)
        liquidMassSourceTerm += cell.fontemassCR / cellLength;
    liquidMassSourceTerm *= (sourceLiquidSpecificHeat) * (sourceTemperature - cell.temp);

    gasMassSourceTerm = cell.fontemassGR / cellLength;
    if (gasMassSourceTerm > 0.)
        gasMassSourceTerm *= (sourceGasSpecificHeat / sourceSpecificHeatRatio) * (sourceTemperature - cell.temp);
    else
        gasMassSourceTerm = 0;

    return TemperatureSourceTerms{
        .gas = gasMassSourceTerm,
        .liquid = liquidMassSourceTerm,
    };
}

}  // namespace

void computeTemperature(const ThermalState &state, int cellIndex, double previousTemperature, int steadyMode) {
    varGlob1D &globals = *state.globals;
    Cel &cell = state.cells[cellIndex];
    Cel &leftCell = state.cells[cellIndex - 1];
    Cel &rightCell = state.cells[cellIndex + 1];
    if (state.thermalSourceDisabled == 0) {
        TemperatureBalance balance = prepareTemperatureBalance(state, cellIndex, steadyMode);
        double pressureGradient;
        if (cellIndex < state.lastCell - 1)
            pressureGradient = 2. * (rightCell.presaux - cell.pres) * kPascalPerKgfPerCm2 / cell.dx;
        else if (cellIndex == state.lastCell - 1 && state.surfaceChoke.AreaGarg > 0.5 * state.cells[state.lastCell - 1].duto.area)
            pressureGradient = 2. * (rightCell.presaux - cell.pres) * kPascalPerKgfPerCm2 / cell.dx;
        else
            pressureGradient = 2. * (cell.pres - cell.presaux) * kPascalPerKgfPerCm2 / cell.dx;
        if (cellIndex > 0 && leftCell.acsr.tipo == kAccessoryChoke &&
            leftCell.acsr.chk.AreaGarg <= (1e-3 + state.input.master1.razareaativ) * cell.duto.area)
            pressureGradient = 2. * (rightCell.presaux - cell.pres) * kPascalPerKgfPerCm2 / cell.dx;
        else if (cellIndex > 0 && cell.acsr.tipo == kAccessoryChoke &&
                 cell.acsr.chk.AreaGarg <= (1e-3 + state.input.master1.razareaativ) * cell.duto.area) {
            pressureGradient = (cell.pres - leftCell.pres) * kPascalPerKgfPerCm2 / balance.meanCellLength;
        } else if (cellIndex == 0)
            pressureGradient = 2. * (rightCell.presaux - cell.pres) * kPascalPerKgfPerCm2 / cell.dx;

        cell.VTemper = balance.temperatureSpatialCoefficient / balance.timeCoefficient;
        if ((cellIndex == 0 && cell.VTemper < 0.) || (((cellIndex < state.lastCell || cell.VTemper >= 0.) && cellIndex > 0) ||
                                                   (cellIndex == state.lastCell && state.surfaceChokeMassFlag == 1) || (cellIndex == state.lastCell && state.input.chkv == 1))) {
            double temperatureGradient = 0.;
            if (cellIndex > 0)
                temperatureGradient = (cell.temp - leftCell.tempini) / balance.meanCellLength;

            if (cellIndex < state.lastCell)
                if (cell.VTemper < 0)
                    temperatureGradient = (rightCell.tempini - cell.temp) / (0.5 * (rightCell.dx + cell.dx));
            if (cell.acsr.tipo == kAccessoryChoke && cell.acsr.chk.AreaGarg <= (1e-3 + state.input.master1.razareaativ) * cell.duto.area && cell.VTemper <= 0)
                temperatureGradient = 0.;
            if (cellIndex > 0 && leftCell.acsr.tipo == kAccessoryChoke &&
                leftCell.acsr.chk.AreaGarg <= (1e-3 + state.input.master1.razareaativ) * leftCell.duto.area && cell.VTemper >= 0)
                temperatureGradient = 0 * (rightCell.tempini - cell.temp) / (0.5 * (rightCell.dx + cell.dx));
            if ((cellIndex <= 1 && cell.VTemper <= 0) || (cellIndex == state.lastCell && cell.VTemper <= 0))
                temperatureGradient = 0.;
            if (cell.acsr.tipo == kAccessoryVolumetricPump && cell.acsr.bvol.freq > 1.) {
                pressureGradient = (cell.pres - leftCell.pres) * kPascalPerKgfPerCm2 / balance.meanCellLength;
            }

            double kineticTerm = computeKineticTemperatureTerm(state, cellIndex, balance);

            TemperatureSourceTerms sourceTerms =
                computeTemperatureSourceTerms(state, cellIndex, balance.cellLength);

            double latentHeatTerm;
            double phaseChangeMassRate = fabs(cell.FonteMudaFase);
            double phaseChangeSign = 1.;
            if (phaseChangeMassRate > kPhaseChangeFloor)
                phaseChangeSign = cell.FonteMudaFase / phaseChangeMassRate;
            if (state.input.limTransMass < phaseChangeMassRate)
                phaseChangeMassRate = phaseChangeSign * state.input.limTransMass;
            else
                phaseChangeMassRate *= phaseChangeSign;
            if (state.latentHeatEnabled > 0 && state.input.flashCompleto == 0) {
                latentHeatTerm = interpolateLatentHeat(state, cell.presini, cell.temp) * phaseChangeMassRate;
            } else if ((state.input.flashCompleto == 1 || state.input.flashCompleto == 2) && state.latentHeatEnabled > 0)
                latentHeatTerm = (cell.flui.EntalpGas(cell.presini, cell.temp) -
                           cell.flui.EntalpLiq(cell.presini, cell.temp)) * phaseChangeMassRate;
            else
                latentHeatTerm = 0;

            if (state.input.latente == 0)
                latentHeatTerm = 0.;
            else if (state.input.condlatente == 0 && latentHeatTerm < 0)
                latentHeatTerm = 0.;

            double interfaceVoidFraction;
            double leftInterfaceVoidFraction;
            if (balance.gasSuperficialVelocity > 0) {
                interfaceVoidFraction = cell.alf;
                leftInterfaceVoidFraction = cell.alfL;
            } else {
                interfaceVoidFraction = cell.alfR;
                leftInterfaceVoidFraction = cell.alf;
            }
            double slipVelocity;
            if (interfaceVoidFraction > globals.localtiny && interfaceVoidFraction < (1. - globals.localtiny))
                slipVelocity = balance.gasSuperficialVelocity / interfaceVoidFraction - balance.liquidSuperficialVelocity / (1. - interfaceVoidFraction);
            else if (interfaceVoidFraction > globals.localtiny)
                slipVelocity = balance.gasSuperficialVelocity;
            else
                slipVelocity = balance.liquidSuperficialVelocity;
            double interfacialWorkTerm = balance.flowArea * cell.pres * kPascalPerKgfPerCm2 * slipVelocity * (interfaceVoidFraction - leftInterfaceVoidFraction) / cell.dx;

            double leftCellThermalPower = 0.;
            if (cellIndex > 0)
                leftCellThermalPower = leftCell.potTermo;
            cell.temp = ((balance.timeCoefficient / cell.dt) * cell.temp - (-balance.pressureTimeCoefficient * (cell.pres - cell.presini) * kPascalPerKgfPerCm2 / cell.dt) + cell.dTdLCor * (-balance.temperatureSpatialCoefficient * temperatureGradient + balance.pressureSpatialCoefficient * pressureGradient - kineticTerm - (balance.hydrostaticPower - 0. * interfacialWorkTerm) + (leftCellThermalPower + cell.fonteCal) / balance.meanCellLength + sourceTerms.liquid + sourceTerms.gas + balance.heatFlux - latentHeatTerm) - (balance.rc - balance.rp) * (1 - balance.voidFraction) * cell.pres * kPascalPerKgfPerCm2 * balance.flowArea * (cell.bet - cell.betini) / (balance.liquidDensity * cell.dt)) / (balance.timeCoefficient / cell.dt);


            if (fabs(cell.temp - cell.tempini) / cell.dt > 10.) {
                cell.temp = cell.tempini +
                                 (fabs(cell.temp - cell.tempini) / (cell.temp - cell.tempini)) * 10 * cell.dt;
            } else if (fabs(cell.temp - cell.tempini) / cell.dt > 1. && fabs(cell.VTemper) > 10 * balance.referenceMixtureVelocity)
                cell.temp = cell.tempini;
            double minimumTemperature = kMinimumTemperatureCelsius;
            if (state.input.usaTabela == 1)
                minimumTemperature = state.input.tabent.tmin + 1.;
            if (cell.temp < minimumTemperature)
                cell.temp = minimumTemperature;
            if (cell.temp > kMaximumTemperatureCelsius)
                cell.temp = kMaximumTemperatureCelsius;
        } else if (cellIndex == state.lastCell) {
            if (globals.chaverede == 0 || state.endNode == 1 || globals.chaveRedeParalela == 1)
                cell.temp = cell.calor.Textern1;
            else
                cell.temp = state.gasSurfaceTemperature;
        }
    } else {
        cell.temp = cell.calor.Textern1;
    }
}

void computeThermalMassTransfer(const ThermalState &state, int cellIndex) {
    varGlob1D &globals = *state.globals;
    Cel &cell = state.cells[cellIndex];
    Cel &leftCell = state.cells[cellIndex - 1];
    Cel &rightCell = state.cells[cellIndex + 1];
    const bool hasRightCell = cellIndex < state.lastCell;

    double cellLength = cell.dx;
    double meanCellLength = 0.5 * (cell.dx + leftCell.dx);
    double diameter = cell.duto.a;
    double flowArea = 0.25 * M_PI * diameter * diameter;
    double meanVoidFraction = cell.alf;
    double betmed = cell.bet;
    double meanSuperficialGasVelocity;

    if (cell.alf > globals.localtiny)
        meanSuperficialGasVelocity = cell.QG / flowArea;
    else {
        meanSuperficialGasVelocity = 0.;
    }
    double meanSuperficialLiquidVelocity;
    if (cell.alf < 1. - globals.localtiny)
        meanSuperficialLiquidVelocity = cell.QL / flowArea;
    else {
        meanSuperficialLiquidVelocity = 0.;
    }
    double producedLiquidDensity = cell.flui.MasEspLiq(cell.pres, cell.temp);
    double completionFluidDensity = cell.fluicol.MasEspFlu(cell.pres, cell.temp);
    double liquidDensity = (1. - betmed) * producedLiquidDensity + betmed * completionFluidDensity;
    double gasDensity = cell.flui.MasEspGas(cell.pres, cell.temp);
    double liquidSpecificHeat = (1. - betmed) * cell.flui.CalorLiq(cell.pres, cell.temp) + betmed * cell.fluicol.CalorLiq(cell.pres, cell.temp);
    double liquidSpecificHeatConstantVolume = liquidSpecificHeat;
    double gasSpecificHeat = cell.flui.CalorGas(cell.pres, cell.temp);
    double gasSpecificHeatConstantVolume = cell.flui.CalorGasVolMod(cell.pres, cell.temp);
    double liquidJouleThomson = (1. - betmed) * cell.flui.JTL(cell.pres, cell.temp) - betmed / completionFluidDensity;
    double gasJouleThomson = cell.flui.JTG(cell.pres, cell.temp);
    double hydrostaticPower = (liquidDensity * meanSuperficialLiquidVelocity + gasDensity * meanSuperficialGasVelocity) * flowArea * kGravity * sin(cell.duto.teta);

    cell.calor.Tint = cell.temp;
    cell.calor.Vint = meanSuperficialGasVelocity + meanSuperficialLiquidVelocity;
    cell.calor.dt = cell.dt;
    double liquidConductivity = (1. - betmed) * cell.flui.CondLiq(cell.pres, cell.temp) + betmed * cell.fluicol.CondLiq(cell.pres, cell.temp);
    cell.calor.kint = liquidConductivity * (1 - meanVoidFraction) + cell.flui.CondGas(cell.pres, cell.temp) * meanVoidFraction;
    cell.calor.cpint = liquidSpecificHeat * (1 - meanVoidFraction) + gasSpecificHeat * meanVoidFraction;
    cell.calor.rhoint = liquidDensity * (1 - meanVoidFraction) + gasDensity * meanVoidFraction;
    double liquidViscosity = (1. - betmed) * cell.flui.ViscOleo(cell.pres, cell.temp) + betmed * cell.fluicol.VisFlu(cell.pres, cell.temp);
    cell.calor.viscint = liquidViscosity * (1 - meanVoidFraction) * kPascalSecondPerCentipoise + cell.flui.ViscGas(cell.pres, cell.temp) * meanVoidFraction * kPascalSecondPerCentipoise;
    double heatFlux = cell.calor.transtrans();

    double timeCoefficient = (liquidDensity * (1 - meanVoidFraction) * liquidSpecificHeatConstantVolume + gasDensity * meanVoidFraction * gasSpecificHeatConstantVolume) * flowArea;
    double pressureTimeCoefficient;
    pressureTimeCoefficient = -cell.flui.CalorGasPresMod(cell.pres, cell.temp) * (gasDensity * meanVoidFraction * flowArea);
    double temperatureSpatialCoefficient = (liquidDensity * meanSuperficialLiquidVelocity * liquidSpecificHeat + gasDensity * meanSuperficialGasVelocity * gasSpecificHeat) * flowArea;
    double pressureSpatialCoefficient = (liquidDensity * meanSuperficialLiquidVelocity * liquidJouleThomson + gasDensity * meanSuperficialGasVelocity * gasJouleThomson) * flowArea;
    double pressureGradient;
    if (hasRightCell)
        pressureGradient = 2. * (rightCell.presaux - cell.pres) * kPascalPerKgfPerCm2 / cell.dx;
    else
        pressureGradient = 2. * (cell.pres - cell.presaux) * kPascalPerKgfPerCm2 / cell.dx;
    if (cell.acsr.tipo == kAccessoryChoke && cell.acsr.chk.AreaGarg <= (1e-3 + state.input.master1.razareaativ) * cell.duto.area)
        pressureGradient = 2. * (cell.pres - cell.presaux) * kPascalPerKgfPerCm2 / cell.dx;
    cell.VTemper = temperatureSpatialCoefficient / timeCoefficient;
    double temperatureGradient = (cell.temp - leftCell.temp) / meanCellLength;
    if (hasRightCell)
        if (cell.VTemper < 0)
            temperatureGradient = (rightCell.temp - cell.temp) / (0.5 * (rightCell.dx + cell.dx));
    if (cell.acsr.tipo == kAccessoryChoke && cell.acsr.chk.AreaGarg <= (1e-3 + state.input.master1.razareaativ) * cell.duto.area && cell.VTemper <= 0)
        temperatureGradient = 0.;
    if (leftCell.acsr.tipo == kAccessoryChoke && leftCell.acsr.chk.AreaGarg <= (1e-3 + state.input.master1.razareaativ) * leftCell.duto.area && cell.VTemper >= 0)
        temperatureGradient = 0.;
    if ((cellIndex == 1 && cell.VTemper <= 0) || (cellIndex == state.lastCell && cell.VTemper <= 0))
        temperatureGradient = 0.;
    if (cell.acsr.tipo == kAccessoryPump && cell.acsr.bcs.freqnova > 1.) {
        pressureGradient = 0.;
    }
    if (cell.acsr.tipo == kAccessoryVolumetricPump && cell.acsr.bvol.freq > 1.) {
        pressureGradient = (cell.pres - leftCell.pres) * kPascalPerKgfPerCm2 / meanCellLength;
    }
    if (cell.acsr.tipo == kAccessoryMultiPump && cell.acsr.multibcs.freqnova > 1.) {
        pressureGradient = 0.;
    }

    double kineticTerm = 0;
    double leftCellLeftFaceMeanGasVelocity = 0;
    double leftCellLeftFaceMeanLiquidVelocity = 0;
    double leftFaceMeanGasVelocity = 0;
    double leftFaceMeanLiquidVelocity = 0.;
    if (cell.acsr.tipo == kAccessoryNone && leftCell.acsr.tipo == kAccessoryNone && cellIndex > 1) {

        double kineticCellLength = leftCell.dx;
        double leftCellDiameter = leftCell.duto.a;
        double leftCellFlowArea = 0.25 * M_PI * leftCellDiameter * leftCellDiameter;

        double leftFaceVoidFraction;
        if (cell.QL > 0)
            leftFaceVoidFraction = leftCell.alf;
        else
            leftFaceVoidFraction = cell.alf;
        double leftCellLeftFaceVoidFraction;
        if (leftCell.QL > 0)
            leftCellLeftFaceVoidFraction = state.cells[cellIndex - 2].alf;
        else
            leftCellLeftFaceVoidFraction = leftCell.alf;

        if (leftCellLeftFaceVoidFraction > 1e-3) {
            leftCellLeftFaceMeanGasVelocity = leftCell.QG / (leftCellFlowArea);
            leftCellLeftFaceMeanGasVelocity /= leftCellLeftFaceVoidFraction;
        }

        if (leftCellLeftFaceVoidFraction < 1. - 1e-3) {
            leftCellLeftFaceMeanLiquidVelocity = leftCell.QL / (leftCellFlowArea);
            leftCellLeftFaceMeanLiquidVelocity /= (1. - leftCellLeftFaceVoidFraction);
        }

        if (leftFaceVoidFraction > 1e-3) {
            leftFaceMeanGasVelocity = meanSuperficialGasVelocity;
            leftFaceMeanGasVelocity /= leftFaceVoidFraction;
        }

        if (leftFaceVoidFraction < 1. - 1e-3) {
            leftFaceMeanLiquidVelocity = meanSuperficialLiquidVelocity;
            leftFaceMeanLiquidVelocity /= (1. - leftFaceVoidFraction);
        }

        kineticTerm = (cell.MC - leftCell.Mliqini) * leftFaceMeanGasVelocity * (leftFaceMeanGasVelocity - leftCellLeftFaceMeanGasVelocity) / kineticCellLength + leftCell.Mliqini * leftFaceMeanLiquidVelocity * (meanSuperficialLiquidVelocity - leftCellLeftFaceMeanLiquidVelocity) / kineticCellLength;
    }

    TemperatureSourceTerms sourceTerms =
        computeThermalMassTransferSourceTerms(state, cellIndex, cellLength);

    double latentHeatTerm;
    if (state.latentHeatEnabled > 0 && state.input.flashCompleto == 0) {
        latentHeatTerm = interpolateLatentHeat(state, cell.pres, cell.temp);
    } else if (state.input.flashCompleto == 1)
        latentHeatTerm = (cell.flui.EntalpGas(cell.pres, cell.temp) -
                   cell.flui.EntalpLiq(cell.pres, cell.temp));
    else
        latentHeatTerm = 0;

    double interfaceVoidFraction;
    double leftInterfaceVoidFraction;
    if (meanSuperficialGasVelocity > 0) {
        interfaceVoidFraction = cell.alf;
        leftInterfaceVoidFraction = cell.alfL;
    } else {
        interfaceVoidFraction = cell.alfR;
        leftInterfaceVoidFraction = cell.alf;
    }
    double slipVelocity;
    if (interfaceVoidFraction > globals.localtiny && interfaceVoidFraction < (1. - globals.localtiny))
        slipVelocity = meanSuperficialGasVelocity / interfaceVoidFraction - meanSuperficialLiquidVelocity / (1. - interfaceVoidFraction);
    else if (interfaceVoidFraction > globals.localtiny)
        slipVelocity = meanSuperficialGasVelocity;
    else
        slipVelocity = meanSuperficialLiquidVelocity;
    double interfacialWorkTerm = flowArea * cell.pres * kPascalPerKgfPerCm2 * slipVelocity * (interfaceVoidFraction - leftInterfaceVoidFraction) / cell.dx;


    cell.FonteMudaFase = (-(timeCoefficient / cell.dt) * (cell.temp - cell.tempini) - (pressureTimeCoefficient * (cell.pres - cell.presini) * kPascalPerKgfPerCm2 / cell.dt) - temperatureSpatialCoefficient * temperatureGradient + pressureSpatialCoefficient * pressureGradient - kineticTerm - (hydrostaticPower - 0. * interfacialWorkTerm) + leftCell.potB / meanCellLength + sourceTerms.liquid + sourceTerms.gas + heatFlux - (completionFluidDensity - producedLiquidDensity) * (1 - meanVoidFraction) * cell.pres * kPascalPerKgfPerCm2 * flowArea * (cell.bet - cell.betini) / (liquidDensity * cell.dt)); // / (timeCoefficient / cell.dt);

    cell.FonteMudaFase /= latentHeatTerm;
}

namespace {

InletMassTransferSeed initializeDistributedMassTransferInlet(
    const ThermalState &state, int cellIndex) {
    double inletLiquidDensity;
    double inletOilVolumeFactor;
    double inletSolutionGasRatio;
    double inletSolutionGasPressureDerivative;
    Cel &cell = state.cells[cellIndex];
    state.cells[0].transmassLini = state.cells[0].transmassL;
    state.cells[0].transmassL = 0.;
    state.sourceUpdater(cellIndex);
    state.cells[1].fontemassLLini = state.cells[1].fontemassLL;
    state.cells[1].fontemassLL = state.cells[0].fontemassLR;
    state.cells[1].fontemassCLini = state.cells[1].fontemassCL;
    state.cells[1].fontemassCL = state.cells[0].fontemassCR;
    state.cells[1].fontemassGLini = state.cells[1].fontemassGL;
    state.cells[1].fontemassGL = state.cells[0].fontemassGR;
    inletLiquidDensity =
        (1 - cell.bet) *
            cell.flui.MasEspLiq(
                cell.pres, cell.temp) +
        cell.bet *
            cell.fluicol.MasEspFlu(
                cell.pres, cell.temp);
    inletOilVolumeFactor = cell.flui.BOFunc(
        cell.pres, cell.temp);
    inletSolutionGasRatio = cell.flui.RS(
        cell.pres, cell.temp);
    if (state.input.flashCompleto != 2 || state.input.miniTabAtraso > 0) {
        double inletOilVolumeFactorAtPerturbedPressure = cell.flui.BOFunc(
            cell.pres * kDerivativePerturbationFactor,
            cell.temp);
        double inletSolutionGasRatioAtPerturbedPressure = cell.flui.RS(
            cell.pres * kDerivativePerturbationFactor,
            cell.temp);
        inletSolutionGasPressureDerivative =
            (inletSolutionGasRatio / inletOilVolumeFactor -
             inletSolutionGasRatioAtPerturbedPressure / inletOilVolumeFactorAtPerturbedPressure) /
            (cell.pres * 0.001);
    } else {
        ProFlu flutemp = cell.flui;
        flutemp.atualizaPropComp(
            cell.pres * kDerivativePerturbationFactor,
            cell.temp, flutemp.dCalculatedBeta,
            flutemp.oCalculatedLiqComposition,
            flutemp.oCalculatedVapComposition, state.input.pocinjec);
        double inletOilVolumeFactorAtPerturbedPressure = flutemp.BOFunc(
            cell.pres * kDerivativePerturbationFactor,
            cell.temp);
        double inletSolutionGasRatioAtPerturbedPressure = flutemp.RS(
            cell.pres * kDerivativePerturbationFactor,
            cell.temp);
        inletSolutionGasPressureDerivative =
            (inletSolutionGasRatio / inletOilVolumeFactor -
             inletSolutionGasRatioAtPerturbedPressure / inletOilVolumeFactorAtPerturbedPressure) /
            (cell.pres * 0.001);
    }
    ProFlu flutemp = cell.flui;
    return InletMassTransferSeed{
        .liquidDensity = inletLiquidDensity,
        .oilVolumeFactor = inletOilVolumeFactor,
        .solutionGasRatio = inletSolutionGasRatio,
        .solutionGasPressureDerivative = inletSolutionGasPressureDerivative};
}

DistributedMassTransferProperties prepareDistributedMassTransferProperties(
    const ThermalState &state, int cellIndex, double meanTemperature, ProFlu &leftCellLeftFaceFluid,
    ProFlu &leftFaceFluid) {
    Cel &cell = state.cells[cellIndex];
    Cel &leftCell = state.cells[cellIndex - 1];
    const bool useTabulatedPvt =
        state.input.flashCompleto != 2 || state.input.miniTabAtraso > 0;
    const bool cellHasNegativeLiquidMassFlowRate = cell.Mliqini < 0;
    const bool leftCellHasNegativeLiquidMassFlowRate = leftCell.Mliqini < 0;
    double leftFaceWaterCut;
    double leftCellLeftFaceWaterCut;
    double leftCellOilVolumeFactorForWaterCut = leftCell.flui.BOFunc(
        leftCell.pres, leftCell.temp);
    double leftCellWaterVolumeFactorForWaterCut = leftCell.flui.BAFunc(
        leftCell.pres, leftCell.temp);
    double leftCellWaterCut = leftCell.flui.BSW * leftCellWaterVolumeFactorForWaterCut /
                 (leftCellOilVolumeFactorForWaterCut + leftCellWaterVolumeFactorForWaterCut * leftCell.flui.BSW -
                  leftCell.flui.BSW * leftCellOilVolumeFactorForWaterCut);
    if (cell.Mliqini < 0.) {
        leftFaceFluid = cell.flui;
        double leftFaceDonorOilVolumeFactor = leftFaceFluid.BOFunc(cell.pres, cell.temp);
        double leftFaceDonorWaterVolumeFactor = leftFaceFluid.BAFunc(cell.pres, cell.temp);
        leftFaceWaterCut = leftFaceFluid.BSW * leftFaceDonorWaterVolumeFactor / (leftFaceDonorOilVolumeFactor + leftFaceDonorWaterVolumeFactor * leftFaceFluid.BSW - leftFaceFluid.BSW * leftFaceDonorOilVolumeFactor);
    } else {
        leftFaceFluid = leftCell.flui;
        double leftFaceDonorOilVolumeFactor = leftFaceFluid.BOFunc(
            leftCell.pres, leftCell.temp);
        double leftFaceDonorWaterVolumeFactor = leftFaceFluid.BAFunc(
            leftCell.pres, leftCell.temp);
        leftFaceWaterCut = leftFaceFluid.BSW * leftFaceDonorWaterVolumeFactor / (leftFaceDonorOilVolumeFactor + leftFaceDonorWaterVolumeFactor * leftFaceFluid.BSW - leftFaceFluid.BSW * leftFaceDonorOilVolumeFactor);
    }
    if (leftCellHasNegativeLiquidMassFlowRate) {
        leftCellLeftFaceFluid = leftCell.flui;
        double leftCellLeftFaceDonorOilVolumeFactor = leftCellLeftFaceFluid.BOFunc(
            leftCell.pres, leftCell.temp);
        double leftCellLeftFaceDonorWaterVolumeFactor = leftCellLeftFaceFluid.BAFunc(
            leftCell.pres, leftCell.temp);
        leftCellLeftFaceWaterCut = leftCellLeftFaceFluid.BSW * leftCellLeftFaceDonorWaterVolumeFactor / (leftCellLeftFaceDonorOilVolumeFactor + leftCellLeftFaceDonorWaterVolumeFactor * leftCellLeftFaceFluid.BSW - leftCellLeftFaceFluid.BSW * leftCellLeftFaceDonorOilVolumeFactor);

    } else {
        if (cellIndex > 1) {
            leftCellLeftFaceFluid = state.cells[cellIndex - 2].flui;
            double leftCellLeftFaceDonorOilVolumeFactor = leftCellLeftFaceFluid.BOFunc(
                state.cells[cellIndex - 2].pres, state.cells[cellIndex - 2].temp);
            double leftCellLeftFaceDonorWaterVolumeFactor = leftCellLeftFaceFluid.BAFunc(
                state.cells[cellIndex - 2].pres, state.cells[cellIndex - 2].temp);
            leftCellLeftFaceWaterCut = leftCellLeftFaceFluid.BSW * leftCellLeftFaceDonorWaterVolumeFactor /
                  (leftCellLeftFaceDonorOilVolumeFactor + leftCellLeftFaceDonorWaterVolumeFactor * leftCellLeftFaceFluid.BSW - leftCellLeftFaceFluid.BSW * leftCellLeftFaceDonorOilVolumeFactor);
        } else {
            leftCellLeftFaceFluid = leftCell.flui;
            double leftCellLeftFaceDonorOilVolumeFactor = leftCellLeftFaceFluid.BOFunc(
                leftCell.pres, leftCell.temp);
            double leftCellLeftFaceDonorWaterVolumeFactor = leftCellLeftFaceFluid.BAFunc(
                leftCell.pres, leftCell.temp);
            leftCellLeftFaceWaterCut = leftCellLeftFaceFluid.BSW * leftCellLeftFaceDonorWaterVolumeFactor /
                  (leftCellLeftFaceDonorOilVolumeFactor + leftCellLeftFaceDonorWaterVolumeFactor * leftCellLeftFaceFluid.BSW - leftCellLeftFaceFluid.BSW * leftCellLeftFaceDonorOilVolumeFactor);
        }
    }

    double liquidDensity;
    double gasDensity;
    double betI;

    // compositional case

    liquidDensity = leftFaceFluid.MasEspLiq(cell.presaux, meanTemperature);
    if (cellHasNegativeLiquidMassFlowRate)
        betI = cell.bet; // beta test
    else
        betI = cell.betL;

    gasDensity = leftFaceFluid.MasEspGas(cell.presaux, meanTemperature);

    double betL = leftCell.betL;
    if (leftCellHasNegativeLiquidMassFlowRate)
        betL = leftCell.bet; // beta test

    if (cellIndex > 0)
        betI = leftCell.betPigD;
    if (cellHasNegativeLiquidMassFlowRate)
        betI = cell.betPigE; // beta test
    if (cellIndex > 1)
        betL = state.cells[cellIndex - 2].betPigD;
    if (leftCellHasNegativeLiquidMassFlowRate)
        betL = leftCell.betPigE; // beta test

    double mixtureLiquidDensity = (1 - betI) * liquidDensity +
                  betI * cell.fluicol.MasEspFlu(
                             cell.presaux, meanTemperature);
    double leftFaceOilVolumeFactor;
    double leftFaceSolutionGasRatio;
    double leftFaceOilVolumeFactorAtPerturbedPressure;
    double leftFaceSolutionGasRatioAtPerturbedPressure;
    if (useTabulatedPvt) {
        leftFaceOilVolumeFactor = leftFaceFluid.BOFunc(cell.presaux, meanTemperature);
        leftFaceSolutionGasRatio = leftFaceFluid.RS(cell.presaux, meanTemperature);
        leftFaceOilVolumeFactorAtPerturbedPressure = leftFaceFluid.BOFunc(cell.presaux * kDerivativePerturbationFactor, meanTemperature);
        leftFaceSolutionGasRatioAtPerturbedPressure = leftFaceFluid.RS(cell.presaux * kDerivativePerturbationFactor, meanTemperature);
    } else {
        leftFaceOilVolumeFactor = leftFaceFluid.BOFunc(cell.presaux, meanTemperature);
        leftFaceSolutionGasRatio = leftFaceFluid.RS(cell.presaux, meanTemperature);
        leftFaceFluid.atualizaPropComp(
            cell.presaux * kDerivativePerturbationFactor, meanTemperature, leftFaceFluid.dCalculatedBeta,
            leftFaceFluid.oCalculatedLiqComposition,
            leftFaceFluid.oCalculatedVapComposition, state.input.pocinjec);
        leftFaceOilVolumeFactorAtPerturbedPressure = leftFaceFluid.BOFunc(cell.presaux * kDerivativePerturbationFactor, meanTemperature);
        leftFaceSolutionGasRatioAtPerturbedPressure = leftFaceFluid.RS(cell.presaux * kDerivativePerturbationFactor, meanTemperature);
    } // compositional case
    double leftFaceSolutionGasPressureDerivative =
        (leftFaceSolutionGasRatio / leftFaceOilVolumeFactor - leftFaceSolutionGasRatioAtPerturbedPressure / leftFaceOilVolumeFactorAtPerturbedPressure) / (cell.presaux * 0.001);
    double leftCellOilVolumeFactor;
    double leftCellSolutionGasRatio;
    double leftCellOilVolumeFactorAtPerturbedPressure;
    double leftCellSolutionGasRatioAtPerturbedPressure;
    if (useTabulatedPvt) {
        leftCellOilVolumeFactor = leftCell.flui.BOFunc(
            leftCell.pres, leftCell.temp);
        leftCellSolutionGasRatio = leftCell.flui.RS(
            leftCell.pres, leftCell.temp);
        leftCellOilVolumeFactorAtPerturbedPressure = leftCell.flui.BOFunc(
            leftCell.pres * kDerivativePerturbationFactor, leftCell.temp);
        leftCellSolutionGasRatioAtPerturbedPressure = leftCell.flui.RS(
            leftCell.pres * kDerivativePerturbationFactor, leftCell.temp);
    } else {
        leftCellOilVolumeFactor = leftCell.flui.BOFunc(
            leftCell.pres, leftCell.temp);
        leftCellSolutionGasRatio = leftCell.flui.RS(
            leftCell.pres, leftCell.temp);
        ProFlu flutemp = leftCell.flui;
        flutemp.atualizaPropComp(
            leftCell.pres * kDerivativePerturbationFactor,
            leftCell.temp, flutemp.dCalculatedBeta,
            flutemp.oCalculatedLiqComposition,
            flutemp.oCalculatedVapComposition, state.input.pocinjec);
        leftCellOilVolumeFactorAtPerturbedPressure = flutemp.BOFunc(
            leftCell.pres * kDerivativePerturbationFactor, leftCell.temp);
        leftCellSolutionGasRatioAtPerturbedPressure = flutemp.RS(
            leftCell.pres * kDerivativePerturbationFactor, leftCell.temp);
    } // compositional case
    double leftCellSolutionGasPressureDerivative =
        (leftCellSolutionGasRatio / leftCellOilVolumeFactor - leftCellSolutionGasRatioAtPerturbedPressure / leftCellOilVolumeFactorAtPerturbedPressure) / (leftCell.pres * 0.001);
    double leftCellOilVolumeFactorAtPerturbedTemperature = 0.;
    double leftCellSolutionGasRatioAtPerturbedTemperature = 0.;
    double leftCellSolutionGasTemperatureDerivative = 0;
    if (state.input.cicloAcopTerm == 1) {
        if (useTabulatedPvt) {
            leftCellOilVolumeFactorAtPerturbedTemperature = leftCell.flui.BOFunc(
                leftCell.pres,
                leftCell.temp * kDerivativePerturbationFactor);
            leftCellSolutionGasRatioAtPerturbedTemperature = leftCell.flui.RS(
                leftCell.pres,
                leftCell.temp * kDerivativePerturbationFactor);
        } else {
            ProFlu flutemp = leftCell.flui;
            flutemp.atualizaPropComp(
                leftCell.pres,
                leftCell.temp * kDerivativePerturbationFactor,
                flutemp.dCalculatedBeta,
                flutemp.oCalculatedLiqComposition,
                flutemp.oCalculatedVapComposition, state.input.pocinjec);
            leftCellOilVolumeFactorAtPerturbedTemperature = flutemp.BOFunc(
                leftCell.pres,
                leftCell.temp * kDerivativePerturbationFactor);
            leftCellSolutionGasRatioAtPerturbedTemperature = flutemp.RS(
                leftCell.pres,
                leftCell.temp * kDerivativePerturbationFactor);
        } // compositional case
        leftCellSolutionGasTemperatureDerivative = (leftCellSolutionGasRatio / leftCellOilVolumeFactor - leftCellSolutionGasRatioAtPerturbedTemperature / leftCellOilVolumeFactorAtPerturbedTemperature) /
                   (leftCell.temp * 0.001);
    }

    return DistributedMassTransferProperties{
        .leftFaceWaterCut = leftFaceWaterCut,
        .leftCellLeftFaceWaterCut = leftCellLeftFaceWaterCut,
        .leftCellWaterCut = leftCellWaterCut,
        .liquidDensity = liquidDensity,
        .gasDensity = gasDensity,
        .leftFaceComposition = betI,
        .leftCellLeftFaceComposition = betL,
        .mixtureLiquidDensity = mixtureLiquidDensity,
        .leftFaceOilVolumeFactor = leftFaceOilVolumeFactor,
        .leftFaceSolutionGasRatio = leftFaceSolutionGasRatio,
        .leftFaceSolutionGasPressureDerivative = leftFaceSolutionGasPressureDerivative,
        .leftCellOilVolumeFactor = leftCellOilVolumeFactor,
        .leftCellSolutionGasRatio = leftCellSolutionGasRatio,
        .leftCellSolutionGasPressureDerivative = leftCellSolutionGasPressureDerivative,
        .leftCellSolutionGasTemperatureDerivative = leftCellSolutionGasTemperatureDerivative,
    };
}

/// The solution-gas pressure derivative recomputed on the LAST cell. Not pure:
/// BOFunc and RS run against the cell's own fluid object.
///
/// Known defect: the DTransDtT caller discards this and scales by the
/// temperature derivative instead.
[[nodiscard]] double lastCellSolutionGasPressureDerivative(const ThermalState &state,
                                                           int cellIndex) {
    Cel &cell = state.cells[cellIndex];
    double cellOilVolumeFactor;
    double cellSolutionGasRatio;
    double cellOilVolumeFactorAtPerturbedPressure;
    double cellSolutionGasRatioAtPerturbedPressure;
    if (state.input.flashCompleto != 2 ||
        state.input.miniTabAtraso > 0) {
        cellOilVolumeFactor = cell.flui.BOFunc(
            cell.pres, cell.temp);
        cellSolutionGasRatio = cell.flui.RS(
            cell.pres, cell.temp);
        cellOilVolumeFactorAtPerturbedPressure = cell.flui.BOFunc(
            cell.pres * kDerivativePerturbationFactor, cell.temp);
        cellSolutionGasRatioAtPerturbedPressure = cell.flui.RS(
            cell.pres * kDerivativePerturbationFactor, cell.temp);
    } else {
        cellOilVolumeFactor = cell.flui.BOFunc(
            cell.pres, cell.temp);
        cellSolutionGasRatio = cell.flui.RS(
            cell.pres, cell.temp);
        ProFlu flutemp = cell.flui;
        flutemp.atualizaPropComp(
            cell.pres * kDerivativePerturbationFactor, cell.temp,
            flutemp.dCalculatedBeta,
            flutemp.oCalculatedLiqComposition,
            flutemp.oCalculatedVapComposition, state.input.pocinjec);
        cellOilVolumeFactorAtPerturbedPressure = flutemp.BOFunc(
            cell.pres * kDerivativePerturbationFactor, cell.temp);
        cellSolutionGasRatioAtPerturbedPressure = flutemp.RS(
            cell.pres * kDerivativePerturbationFactor, cell.temp);
    } // compositional case
    double cellSolutionGasPressureDerivative =
        (cellSolutionGasRatio / cellOilVolumeFactor - cellSolutionGasRatioAtPerturbedPressure / cellOilVolumeFactorAtPerturbedPressure) / (cell.pres * 0.001);
    return cellSolutionGasPressureDerivative;
}

DistributedMassTransferCoefficients updateDistributedMassTransferDerivatives(
    const ThermalState &state, int cellIndex, double leftCellWaterCut, const ProFlu &leftFaceFluid,
    double leftCellSolutionGasPressureDerivative, double leftCellSolutionGasTemperatureDerivative) {
    Cel &cell = state.cells[cellIndex];
    Cel &leftCell = state.cells[cellIndex - 1];
    const bool isLastCell = cellIndex == state.lastCell;
    double activeDerivative = 1.;
    double pressureLimit = 10;
    if (state.fullModel == 1)
        pressureLimit = 0;
    if (leftCell.pres < pressureLimit ||
        state.massTransferModel != 0)
        activeDerivative = 0.;
    double spatialCoupling = 1.;
    if (cellIndex < 2 || isLastCell)
        spatialCoupling = 0;

    double coefficientFlowArea = cell.dutoL.area;
    leftCell.ativaDeri = activeDerivative;
    leftCell.DTransDtp =
        activeDerivative * coefficientFlowArea * (1. - leftCell.alf) *
        (1. - leftCell.bet) * (1. - leftCellWaterCut) * leftFaceFluid.Deng *
        kAirDensityAtStandardConditions * leftCellSolutionGasPressureDerivative * (kBarrelPerCubicMetre / kCubicFootPerCubicMetre);
    cell.DTransDtpL = leftCell.DTransDtp;
    if (isLastCell) {
        const double leftCellSolutionGasPressureDerivative =
            lastCellSolutionGasPressureDerivative(state, cellIndex);
        cell.DTransDtp =
            activeDerivative * coefficientFlowArea * (1. - cell.alf) *
            (1. - cell.bet) * (1. - leftCellWaterCut) * leftFaceFluid.Deng *
            kAirDensityAtStandardConditions * leftCellSolutionGasPressureDerivative * (kBarrelPerCubicMetre / kCubicFootPerCubicMetre);
    }
    if (state.input.cicloAcopTerm == 1) {
        leftCell.DTransDtT =
            activeDerivative * coefficientFlowArea * (1. - leftCell.alf) *
            (1. - leftCell.bet) * (1. - leftCellWaterCut) * leftFaceFluid.Deng *
            kAirDensityAtStandardConditions * leftCellSolutionGasTemperatureDerivative * (kBarrelPerCubicMetre / kCubicFootPerCubicMetre);
        cell.DTransDtTL = leftCell.DTransDtT;
        if (isLastCell) {
            // Not used: the assignment below reaches for the temperature derivative.
            [[maybe_unused]] const double recomputedPressureDerivative =
                lastCellSolutionGasPressureDerivative(state, cellIndex);
            cell.DTransDtT =
                activeDerivative * coefficientFlowArea * (1. - cell.alf) *
                (1. - cell.bet) * (1. - leftCellWaterCut) * leftFaceFluid.Deng *
                kAirDensityAtStandardConditions * leftCellSolutionGasTemperatureDerivative * (kBarrelPerCubicMetre / kCubicFootPerCubicMetre);
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
    Cel &cell = state.cells[cellIndex];
    Cel &leftCell = state.cells[cellIndex - 1];
    meanTemperature = leftCell.temp;
    if (cell.VTemper < 0.)
        meanTemperature = cell.temp;
    leftMeanTemperature = leftCell.tempL;
    if (leftCell.VTemper < 0.)
        leftMeanTemperature = leftCell.temp;

    leftCell.TMModel = state.massTransferModel;
    if (state.massTransferModel != 3) {
        if ((((leftCell.alf < 0.001) ||
              (leftCell.alf > kDerivativePerturbationFactor) ||
              (leftCell.bet > kDerivativePerturbationFactor &&
               leftCell.alf < kDerivativePerturbationFactor)) &&
             leftAbsoluteSuperficialVelocity < 0.1) ||
            leftCell.flui.RGO >= (*state.globals).RGOMax)
            leftCell.TMModel = 3;
        else if (leftCell.estadoPig == 1)
            leftCell.TMModel = 3;
        else if (leftCell.acsr.tipo == kAccessoryLiquidInjection ||
                 leftCell.acsr.tipo == kAccessoryInflowPerformance ||
                 leftCell.acsr.tipo == kAccessoryLeak ||
                 leftCell.acsr.tipo == kAccessoryRadialPorous ||
                 leftCell.acsr.tipo == kAccessoryPorous2D)
            leftCell.TMModel = 3;
        else if (cellIndex >= 2) {
            if (state.cells[cellIndex - 2].acsr.tipo == kAccessoryChoke &&
                (state.cells[cellIndex - 2].acsr.chk.AreaGarg <
                 (1e-3 + state.input.master1.razareaativ) *
                     state.cells[cellIndex - 2].duto.area))
                leftCell.TMModel = 3;
            else if (state.cells[cellIndex - 2].acsr.tipo == kAccessoryPump ||
                     state.cells[cellIndex - 2].acsr.tipo == 7 ||
                     state.cells[cellIndex - 2].acsr.tipo == kAccessoryMultiPump)
                leftCell.TMModel = 0;
        }
        if (state.input.flashCompleto == 2) {
            double candidateMassFraction = leftCell.flui.FracMass(
                leftCell.pres, leftCell.temp);
            if (candidateMassFraction > 1.0 - 1e-2 || candidateMassFraction < 1e-2)
                leftCell.TMModel = 3;
        }
    }
    if (leftCell.TMModel == 0 &&
        leftCell.alf <= (*state.globals).CritCond)
        leftCell.TMModel = 1;
    if (leftCell.TMModel == 0 && cellIndex == state.lastCell)
        leftCell.TMModel = 1;
    if (cellIndex == state.lastCell)
        cell.FonteMudaFase = 0.;
}

/// Copies the mass-transfer state of the face into its *ini counterparts.
void snapshotMassTransferState(const ThermalState &state, int cellIndex) {
    Cel &cell = state.cells[cellIndex];
    Cel &leftCell = state.cells[cellIndex - 1];
    leftCell.transmassRini = leftCell.transmassR;
    leftCell.FonteMudaFaseini =
        leftCell.FonteMudaFase;
    cell.DTransDxRini = cell.DTransDxR;
    cell.DTransDxLini = cell.DTransDxL;
    cell.DTransDt1ini = cell.DTransDt1;
    cell.DTransDt0ini = cell.DTransDt0;
    cell.DTransDxRpini = cell.DTransDxRp;
    cell.DTransDxLpini = cell.DTransDxLp;
    leftCell.CoefDTLini = leftCell.CoefDTL;
    leftCell.coefTransBetini =
        leftCell.coefTransBet;
    cell.transmassLini = cell.transmassL;
}

void applyDistributedMassTransferModel(
    const ThermalState &state, int cellIndex, double &meanTemperature, double &leftMeanTemperature,
    double leftAbsoluteSuperficialVelocity, const ProFlu &leftCellLeftFaceFluid, const ProFlu &leftFaceFluid,
    const DistributedMassTransferProperties &properties,
    double activeDerivative, double spatialCoupling, double coefficientFlowArea, double leftCellLeftFaceMixtureLiquidDensity, double leftCellLeftFaceOilVolumeFactor,
    double leftCellLeftFaceSolutionGasRatio, double leftCellLeftFaceSolutionGasPressureDerivative) {
    const double leftFaceWaterCut = properties.leftFaceWaterCut;
    const double leftCellLeftFaceWaterCut = properties.leftCellLeftFaceWaterCut;
    const double leftCellWaterCut = properties.leftCellWaterCut;
    const double betI = properties.leftFaceComposition;
    const double betL = properties.leftCellLeftFaceComposition;
    const double liquidDensity = properties.mixtureLiquidDensity;
    const double leftFaceOilVolumeFactor = properties.leftFaceOilVolumeFactor;
    const double leftFaceSolutionGasRatio = properties.leftFaceSolutionGasRatio;
    const double leftFaceSolutionGasPressureDerivative =
        properties.leftFaceSolutionGasPressureDerivative;
    const double leftCellOilVolumeFactor = properties.leftCellOilVolumeFactor;
    const double leftCellSolutionGasRatio = properties.leftCellSolutionGasRatio;
    Cel &cell = state.cells[cellIndex];
    Cel &leftCell = state.cells[cellIndex - 1];
    selectDistributedMassTransferModel(state, cellIndex, meanTemperature, leftMeanTemperature, leftAbsoluteSuperficialVelocity);

    leftCell.fontedissolv = 0.;

    snapshotMassTransferState(state, cellIndex);

    cell.TMModelL = leftCell.TMModel;
    leftCell.FonteMudaFase = 0.;
    if (leftCell.TMModel == 0 ||
        leftCell.TMModel == 1) {

        leftCell.transmassR =
            -(cell.QL * (1 - betI) * (leftFaceFluid.rDgD) * leftFaceFluid.Deng *
              kAirDensityAtStandardConditions * (1. - leftFaceWaterCut) * leftFaceSolutionGasRatio * (kBarrelPerCubicMetre / kCubicFootPerCubicMetre) / leftFaceOilVolumeFactor) +
            (leftCell.QL * (1 - betL) * (leftCellLeftFaceFluid.rDgD) *
             leftCellLeftFaceFluid.Deng * kAirDensityAtStandardConditions * (1. - leftCellLeftFaceWaterCut) * leftCellLeftFaceSolutionGasRatio * (kBarrelPerCubicMetre / kCubicFootPerCubicMetre) /
             leftCellLeftFaceOilVolumeFactor);

        leftCell.transmassR /= leftCell.dx;
        leftCell.transmassR += leftCell.fontedissolv;
        cell.transmassL = leftCell.transmassR;
        leftCell.FonteMudaFase =
            leftCell.transmassR -
            leftCell.DTransDtp * leftCell.d2pdt2 -
            leftCell.DTransDtT * leftCell.dTdtIni;
        if (leftCell.TMModel == 1) {
            leftCell.transmassR -=
                activeDerivative * ((1. - leftCell.bet) *
                         (1. - leftCell.alf) * (1. - leftCellWaterCut) * coefficientFlowArea *
                         (leftCell.flui.rDgD) *
                         leftCell.flui.Deng * kAirDensityAtStandardConditions * leftCellSolutionGasRatio *
                         (kBarrelPerCubicMetre / kCubicFootPerCubicMetre) / leftCellOilVolumeFactor) /
                leftCell.dt;
            leftCell.transmassR +=
                activeDerivative * ((1. - leftCell.betini) *
                         (1. - leftCell.alfini) * (1. - leftCellWaterCut) * coefficientFlowArea *
                         (leftCell.flui.rDgD) *
                         leftCell.flui.Deng * kAirDensityAtStandardConditions * leftCellSolutionGasRatio *
                         (kBarrelPerCubicMetre / kCubicFootPerCubicMetre) / leftCellOilVolumeFactor) /
                leftCell.dt;

            cell.transmassL = leftCell.transmassR;
            leftCell.FonteMudaFase =
                leftCell.transmassR;
        }
        if (leftCell.TMModel == 0) {
            leftCell.FonteMudaFase -=
                activeDerivative * ((1. - leftCell.bet) *
                         (1. - leftCell.alf) * (1. - leftCellWaterCut) * coefficientFlowArea *
                         (leftCell.flui.rDgD) *
                         leftCell.flui.Deng * kAirDensityAtStandardConditions * leftCellSolutionGasRatio *
                         (kBarrelPerCubicMetre / kCubicFootPerCubicMetre) / leftCellOilVolumeFactor) /
                leftCell.dt;
            leftCell.FonteMudaFase +=
                activeDerivative * ((1. - leftCell.betini) *
                         (1. - leftCell.alfini) * (1. - leftCellWaterCut) * coefficientFlowArea *
                         (leftCell.flui.rDgD) *
                         leftCell.flui.Deng * kAirDensityAtStandardConditions * leftCellSolutionGasRatio *
                         (kBarrelPerCubicMetre / kCubicFootPerCubicMetre) / leftCellOilVolumeFactor) /
                leftCell.dt;
        }

        if (leftCell.TMModel == 0) {
            cell.DTransDxR =
                -((1 - betI) * (leftFaceFluid.rDgD) * leftFaceFluid.Deng * kAirDensityAtStandardConditions *
                  (1. - leftFaceWaterCut) * leftFaceSolutionGasRatio * (kBarrelPerCubicMetre / kCubicFootPerCubicMetre) / leftFaceOilVolumeFactor) /
                (liquidDensity * leftCell.dx);
            cell.DtransDxLinear =
                -spatialCoupling * cell.QL *
                    ((1 - betI) * (leftFaceFluid.rDgD) * leftFaceFluid.Deng * kAirDensityAtStandardConditions *
                     (1. - leftFaceWaterCut) * (kBarrelPerCubicMetre / kCubicFootPerCubicMetre) *
                     (leftFaceSolutionGasPressureDerivative * cell.dpresaux)) /
                    (leftCell.dx) +
                spatialCoupling * cell.QL *
                    ((1 - betI) * (leftFaceFluid.rDgD) * leftFaceFluid.Deng * kAirDensityAtStandardConditions *
                     (1. - leftFaceWaterCut) * (kBarrelPerCubicMetre / kCubicFootPerCubicMetre) *
                     (leftFaceSolutionGasPressureDerivative * cell.presaux)) /
                    (leftCell.dx);
            cell.DTransDxRp =
                -spatialCoupling * cell.QL *
                ((1 - betI) * (leftFaceFluid.rDgD) * leftFaceFluid.Deng * kAirDensityAtStandardConditions *
                 (1. - leftFaceWaterCut) * (kBarrelPerCubicMetre / kCubicFootPerCubicMetre) * 0.5 * leftFaceSolutionGasPressureDerivative) /
                (leftCell.dx);
            cell.DTransDxL =
                ((1 - betL) * (leftCellLeftFaceFluid.rDgD) * leftCellLeftFaceFluid.Deng * kAirDensityAtStandardConditions *
                 (1. - leftCellLeftFaceWaterCut) * leftCellLeftFaceSolutionGasRatio * (kBarrelPerCubicMetre / kCubicFootPerCubicMetre) / leftCellLeftFaceOilVolumeFactor) /
                (leftCellLeftFaceMixtureLiquidDensity * leftCell.dx);
            cell.DtransDxLinear =
                cell.DtransDxLinear +
                spatialCoupling * leftCell.QL *
                    ((1 - betL) * (leftCellLeftFaceFluid.rDgD) * leftCellLeftFaceFluid.Deng * kAirDensityAtStandardConditions *
                     (1. - leftCellLeftFaceWaterCut) * (kBarrelPerCubicMetre / kCubicFootPerCubicMetre) *
                     (leftCellLeftFaceSolutionGasPressureDerivative * leftCell.dpresaux)) /
                    (leftCell.dx) -
                spatialCoupling * leftCell.QL *
                    ((1 - betL) * (leftCellLeftFaceFluid.rDgD) * leftCellLeftFaceFluid.Deng * kAirDensityAtStandardConditions *
                     (1. - leftCellLeftFaceWaterCut) * (kBarrelPerCubicMetre / kCubicFootPerCubicMetre) *
                     (leftCellLeftFaceSolutionGasPressureDerivative * leftCell.presaux)) /
                    (leftCell.dx);
            cell.DTransDxLp =
                spatialCoupling * leftCell.QL *
                ((1 - betL) * (leftCellLeftFaceFluid.rDgD) * leftCellLeftFaceFluid.Deng * kAirDensityAtStandardConditions *
                 (1. - leftCellLeftFaceWaterCut) * (kBarrelPerCubicMetre / kCubicFootPerCubicMetre) * 0.5 * leftCellLeftFaceSolutionGasPressureDerivative) /
                (leftCell.dx);
            cell.DTransDt1 =
                -activeDerivative * ((1. - leftCellWaterCut) * coefficientFlowArea *
                          (leftCell.flui.rDgD) *
                          leftCell.flui.Deng * kAirDensityAtStandardConditions * leftCellSolutionGasRatio *
                          (kBarrelPerCubicMetre / kCubicFootPerCubicMetre) / leftCellOilVolumeFactor);
            cell.DTransDt0 = -cell.DTransDt1;

            leftCell.CoefDTR =
                -((1. - leftCell.bet) * (1. - leftCellWaterCut) * coefficientFlowArea *
                  (leftCell.flui.rDgD) *
                  leftCell.flui.Deng * kAirDensityAtStandardConditions * leftCellSolutionGasRatio *
                  (kBarrelPerCubicMetre / kCubicFootPerCubicMetre) / leftCellOilVolumeFactor);
            leftCell.CoefDTL = -leftCell.CoefDTR;
            leftCell.coefTransBet =
                ((1. - leftCellWaterCut) * coefficientFlowArea * (leftCell.flui.rDgD) *
                 leftCell.flui.Deng * kAirDensityAtStandardConditions * leftCellSolutionGasRatio *
                 (kBarrelPerCubicMetre / kCubicFootPerCubicMetre) / leftCellOilVolumeFactor);

            cell.transmassL -=
                activeDerivative * ((1. - leftCell.bet) *
                         (1. - leftCell.alf) * (1. - leftCellWaterCut) * coefficientFlowArea *
                         (leftCell.flui.rDgD) *
                         leftCell.flui.Deng * kAirDensityAtStandardConditions * leftCellSolutionGasRatio *
                         (kBarrelPerCubicMetre / kCubicFootPerCubicMetre) / leftCellOilVolumeFactor) /
                leftCell.dt;
            cell.transmassL +=
                activeDerivative * ((1. - leftCell.betini) *
                         (1. - leftCell.alfini) * (1. - leftCellWaterCut) * coefficientFlowArea *
                         (leftCell.flui.rDgD) *
                         leftCell.flui.Deng * kAirDensityAtStandardConditions * leftCellSolutionGasRatio *
                         (kBarrelPerCubicMetre / kCubicFootPerCubicMetre) / leftCellOilVolumeFactor) /
                leftCell.dt;

        } else {
            cell.DTransDxR = 0.;
            cell.DTransDxL = 0.;
            cell.DTransDt1 = 0.;
            cell.DTransDxRp = 0.;
            cell.DTransDxLp = 0.;
            if (state.input.desligaDeriTransMassDTemp == 1) {
                leftCell.DTransDtT = 0;
                cell.DTransDtTL = 0.;
            }
            leftCell.CoefDTR = 0.;
            leftCell.CoefDTL = 0.;
            leftCell.coefTransBet = 0.;
        }
    }
}


/// Zeroes every mass-transfer derivative on the face between cellIndex - 1 and
/// cellIndex. The transmass value itself is left to the caller.
void clearMassTransferDerivatives(const ThermalState &state, int cellIndex) {
    Cel &cell = state.cells[cellIndex];
    Cel &leftCell = state.cells[cellIndex - 1];
    cell.DTransDxR = 0.;
    cell.DTransDxL = 0.;
    cell.DTransDt1 = 0.;
    cell.DTransDt0 = 0.;
    cell.DTransDxRp = 0.;
    cell.DTransDxLp = 0.;
    leftCell.CoefDTR = 0.;
    leftCell.CoefDTL = 0.;
    leftCell.coefTransBet = 0.;
    if (state.input.desligaDeriTransMassDTemp == 1) {
    leftCell.DTransDtT = 0;
    cell.DTransDtTL = 0.;
    }
}

}  // namespace

void updateDistributedMassTransfer(const ThermalState &state) {
    // #pragma omp parallel for num_threads(state.input.nthrd)
    double leftCellLeftFaceMixtureLiquidDensity = 0.;
    double leftCellLeftFaceOilVolumeFactor = 0.;
    double leftCellLeftFaceSolutionGasRatio = 0.;
    double leftCellLeftFaceSolutionGasPressureDerivative = 0.;
    for (int cellIndex = 0; cellIndex <= state.lastCell; cellIndex++) {
        Cel &cell = state.cells[cellIndex];
        Cel &leftCell = state.cells[cellIndex - 1];
        Cel &rightCell = state.cells[cellIndex + 1];
        if (cellIndex != 0 && cellIndex != state.lastCell + 1) {


            state.sourceUpdater(cellIndex);

            if (cellIndex < state.lastCell) {
                rightCell.fontemassLLini = cell.fontemassLR;
                rightCell.fontemassCLini = cell.fontemassCR;
                rightCell.fontemassGLini = cell.fontemassGR;
                rightCell.fontemassLL = cell.fontemassLR;
                rightCell.fontemassCL = cell.fontemassCR;
                rightCell.fontemassGL = cell.fontemassGR;
            }
            double lengthRatio = leftCell.dx / (cell.dx + cell.dxL);
            double leftLengthRatio = cell.dx / (leftCell.dx + leftCell.dxL);
            double meanTemperature = lengthRatio * cell.temp + (1 - lengthRatio) * leftCell.temp;
            double leftMeanTemperature = leftLengthRatio * cell.tempL + (1 - leftLengthRatio) * leftCell.tempL;
            meanTemperature = leftCell.temp;
            if (cell.VTemper < 0.)
                meanTemperature = cell.temp;
            leftMeanTemperature = leftCell.tempL;
            if (leftCell.VTemper < 0.)
                leftMeanTemperature = leftCell.temp;

            double diameter = cell.duto.a;
            double flowArea = 0.25 * M_PI * diameter * diameter;
            double meanSuperficialGasVelocity = (cell.QG) / (flowArea);
            double meanSuperficialLiquidVelocity = cell.QL / (flowArea);
            [[maybe_unused]] double mixtureSuperficialVelocity = meanSuperficialGasVelocity + meanSuperficialLiquidVelocity;
            double leftAbsoluteSuperficialVelocity = (fabs(leftCell.QG) + fabs(leftCell.QL)) / leftCell.duto.area;

            ProFlu leftCellLeftFaceFluid;
            ProFlu leftFaceFluid;
            DistributedMassTransferProperties properties =
                prepareDistributedMassTransferProperties(
                    state, cellIndex, meanTemperature, leftCellLeftFaceFluid, leftFaceFluid);
            double leftCellWaterCut = properties.leftCellWaterCut;
            double liquidDensity = properties.liquidDensity;
            double gasDensity = properties.gasDensity;
            double betI = properties.leftFaceComposition;
            double mixtureLiquidDensity = properties.mixtureLiquidDensity;
            double leftFaceOilVolumeFactor = properties.leftFaceOilVolumeFactor;
            double leftFaceSolutionGasRatio = properties.leftFaceSolutionGasRatio;
            double leftFaceSolutionGasPressureDerivative =
                properties.leftFaceSolutionGasPressureDerivative;
            [[maybe_unused]] double leftCellSolutionGasRatio = properties.leftCellSolutionGasRatio;
            double leftCellSolutionGasPressureDerivative = properties.leftCellSolutionGasPressureDerivative;
            double leftCellSolutionGasTemperatureDerivative =
                properties.leftCellSolutionGasTemperatureDerivative;
            DistributedMassTransferCoefficients coefficients =
                updateDistributedMassTransferDerivatives(
                    state, cellIndex, leftCellWaterCut, leftFaceFluid, leftCellSolutionGasPressureDerivative, leftCellSolutionGasTemperatureDerivative);
            double activeDerivative = coefficients.activeDerivative;
            double spatialCoupling = coefficients.spatialCoupling;
            double coefficientFlowArea = coefficients.flowArea;

            applyDistributedMassTransferModel(
                state, cellIndex, meanTemperature, leftMeanTemperature, leftAbsoluteSuperficialVelocity, leftCellLeftFaceFluid, leftFaceFluid, properties,
                activeDerivative, spatialCoupling, coefficientFlowArea, leftCellLeftFaceMixtureLiquidDensity, leftCellLeftFaceOilVolumeFactor, leftCellLeftFaceSolutionGasRatio, leftCellLeftFaceSolutionGasPressureDerivative);
            if (leftCell.TMModel == -2) {
                double hydrateTransportVelocity;
                if (cell.alfL > (*state.globals).localtiny && betI < (1. - (*state.globals).localtiny))
                    hydrateTransportVelocity = (cell.QG * gasDensity + cell.QL * (1. - betI) * liquidDensity) / (coefficientFlowArea * (cell.alfL * gasDensity + (1. - cell.alfL) * (1. - betI) * liquidDensity));
                else
                    hydrateTransportVelocity = 0.;
                double hydrateMassFraction = leftCell.flui.FracMassHidra(leftCell.pres, leftCell.temp);
                double pressurePerturbationRatio = kDerivativePerturbationFactor;
                double hydrateFractionPressureDerivative = (hydrateMassFraction - leftCell.flui.FracMassHidra(leftCell.pres * pressurePerturbationRatio, leftCell.temp)) / ((1 - pressurePerturbationRatio) * leftCell.pres);
                double pressureGradient;
                pressureGradient = (cell.presaux - leftCell.presaux) / cell.dxL;
                cell.transmassL = leftCell.transmassR = 1 * (cell.alfL * gasDensity + (1. - cell.alfL) * (1. - betI) * liquidDensity) * hydrateTransportVelocity * (pressureGradient * hydrateFractionPressureDerivative) * coefficientFlowArea;

                cell.DTransDxR = 0.;
                cell.DTransDxL = 0.;
                cell.DTransDt1 = 0.;
                cell.DTransDt0 = 0.;
                if (state.input.desligaDeriTransMassDTemp == 1) {
                    leftCell.DTransDtT = 0;
                    cell.DTransDtTL = 0.;
                }
                leftCell.CoefDTR = 0.;
                leftCell.CoefDTL = 0.;
                leftCell.coefTransBet = 0.;
            }
            if (leftCell.TMModel == 3) {
                cell.transmassL = leftCell.transmassR = 0.;
                clearMassTransferDerivatives(state, cellIndex);
            }
            if (leftCell.transmassR > 0 && (cell.alfL > (1. - (*state.globals).localtiny) || cell.betL > (1. - (*state.globals).localtiny))) {
                cell.transmassL = leftCell.transmassR = -(*state.globals).localtiny;
                clearMassTransferDerivatives(state, cellIndex);
            }
            if (leftCell.transmassR < 0 && cell.alfL < (*state.globals).localtiny) {
                cell.transmassL = leftCell.transmassR = (*state.globals).localtiny;
                clearMassTransferDerivatives(state, cellIndex);
            }
            leftCellLeftFaceMixtureLiquidDensity = mixtureLiquidDensity;
            leftCellLeftFaceOilVolumeFactor = leftFaceOilVolumeFactor;
            leftCellLeftFaceSolutionGasRatio = leftFaceSolutionGasRatio;
            leftCellLeftFaceSolutionGasPressureDerivative = leftFaceSolutionGasPressureDerivative;

        } else if (cellIndex == 0) {
            const InletMassTransferSeed seed =
                initializeDistributedMassTransferInlet(state, cellIndex);
            leftCellLeftFaceMixtureLiquidDensity = seed.liquidDensity;
            leftCellLeftFaceOilVolumeFactor = seed.oilVolumeFactor;
            leftCellLeftFaceSolutionGasRatio = seed.solutionGasRatio;
            leftCellLeftFaceSolutionGasPressureDerivative = seed.solutionGasPressureDerivative;
        }
    }
}

namespace {

/// meanVoidFraction is in-out: counter-current gas re-seeds it from alfPigE.
[[nodiscard]] SlugClosure applySlugDriftClosure(
    const ThermalState &state, int cellIndex, double superficialGasVelocity,
    double rightGasDensity, double rightLiquidDensity,
    double &meanVoidFraction) {
    Cel &cell = state.cells[cellIndex];
    if (superficialGasVelocity < 0)
        meanVoidFraction = cell.alfPigE;
    double c0 = 1.2;
    double meanDiameter = cell.duto.a;
    if (cell.MC >= 0)
        meanDiameter = cell.dutoL.a;
    double inclinationSign = 1.;
    if (cell.duto.teta < 0.)
        inclinationSign = -1.;
    double ud = inclinationSign * 0.32 * sqrt(kGravity * meanDiameter);
    if (fabs(rightGasDensity) / rightLiquidDensity > 0.9) {
        c0 = 1.;
        ud = 0.;
    }
    return SlugClosure{.c0 = c0, .ud = ud, .meanDiameter = meanDiameter};
}

/// Which face the upstream properties come from, decided by the sign of the gas
/// flow rate: the left duct and cell when QG >= 0, this cell's duct otherwise.
[[nodiscard]] UpstreamFaceBasis upstreamFaceBasisOf(const ThermalState &state, int cellIndex) {
    Cel &cell = state.cells[cellIndex];
    Cel &leftCell = state.cells[cellIndex - 1];
    double gasDensity;
    double flowArea;
    double noSlipLiquidHoldup;
    if (cell.QG >= 0) {
    flowArea = cell.dutoL.area;
    noSlipLiquidHoldup = 1. - cell.alfL;
    gasDensity = leftCell.rgCi;
    } else {
    gasDensity = cell.rgCi;
    flowArea = cell.duto.area;
    noSlipLiquidHoldup = 1. - cell.alf;
    }
    return UpstreamFaceBasis{.gasDensity = gasDensity,
                             .flowArea = flowArea,
                             .noSlipLiquidHoldup = noSlipLiquidHoldup};
}

/// Sets the flow-partition terms and the branch flag. `term1` is 1. when the
/// liquid partition is full, 0. when empty.
void setFlowPartitionTerms(const ThermalState &state, int cellIndex, double term1,
                           int &branchFlag, int flagValue) {
    state.cells[cellIndex].term1 = term1;
    state.cells[cellIndex].term2 = 0.;
    branchFlag = flagValue;
}

/// Clears the drift-flux closure of a cell: unit distribution, no drift
/// velocity, no flow-pattern label.
void clearDriftClosure(const ThermalState &state, int cellIndex) {
    Cel &cell = state.cells[cellIndex];
    cell.c0 = 1;
    cell.ud = 0;
    cell.arranjo = 0;
}

void selectAndApplyInteriorFlowRegime(
    const ThermalState &state, int cellIndex, Vcr<int> &bif, double superficialGasVelocity,
    double superficialLiquidVelocity, double leftSuperficialGasVelocity, double leftSuperficialLiquidVelocity, double rightSuperficialGasVelocity, double rightSuperficialLiquidVelocity,
    double rightGasDensity, double rightLiquidDensity, double gasDensity, double liquidDensity, double flowArea) {
    varGlob1D &globals = *state.globals;
    Cel &cell = state.cells[cellIndex];
    Cel &leftCell = state.cells[cellIndex - 1];
    Cel &rightCell = state.cells[cellIndex + 1];
    const bool noGasSourceOnEitherFace =
        cell.fontemassGL <= globals.localtiny * 1e-5 && cell.fontemassGR <= globals.localtiny * 1e-5;
    const bool noLiquidSourceOnEitherFace =
        (cell.fontemassLL + cell.fontemassCL) <= globals.localtiny * 1e-5 && (cell.fontemassLR + cell.fontemassCR) <= globals.localtiny * 1e-5;
    bif[cellIndex] = 1;


    if (leftCell.alfPigD <= globals.localtiny && cell.alfPigE <= globals.localtiny && noGasSourceOnEitherFace) {
        setFlowPartitionTerms(state, cellIndex, 1., bif[cellIndex], 0);
        clearDriftClosure(state, cellIndex);
    } else if (leftCell.alfPigD >= (1. - globals.localtiny) && cell.alfPigE >= (1. - globals.localtiny) && (noLiquidSourceOnEitherFace)) {
        setFlowPartitionTerms(state, cellIndex, 0., bif[cellIndex], 0);
        clearDriftClosure(state, cellIndex);
    } else if (leftCell.acsr.tipo == kAccessoryChoke && leftCell.acsr.chk.AreaGarg <= (1e-3)) {
        setFlowPartitionTerms(state, cellIndex, 0., bif[cellIndex], 0);
        clearDriftClosure(state, cellIndex);
    }

    else if (superficialGasVelocity >= 0 && leftCell.alfPigD <= globals.localtiny && (noGasSourceOnEitherFace)) {
        setFlowPartitionTerms(state, cellIndex, 1., bif[cellIndex], 0);
        clearDriftClosure(state, cellIndex);
        if (fabs(superficialGasVelocity) <= 1e-15 && cell.alfPigE > (1. - globals.localtiny) && superficialLiquidVelocity < 0 && leftSuperficialLiquidVelocity < 0 && cell.duto.teta > 0) {
            setFlowPartitionTerms(state, cellIndex, 0., bif[cellIndex], 0);
        }
        if (fabs(superficialGasVelocity) <= 1e-15 && cell.alfPigE > (1. - 10 * globals.localtiny) && cell.duto.teta < 0) {
            setFlowPartitionTerms(state, cellIndex, 0., bif[cellIndex], 1);
        }
        if (fabs(superficialGasVelocity) <= 1e-15 && cell.alf >= rightCell.alf && rightSuperficialGasVelocity < 0) {
            bif[cellIndex] = 1;
        }
    } else if (superficialGasVelocity <= 0 && cell.alfPigE <= globals.localtiny && (noGasSourceOnEitherFace)) {
        setFlowPartitionTerms(state, cellIndex, 1., bif[cellIndex], 0);
        clearDriftClosure(state, cellIndex);
        if (fabs(superficialGasVelocity) <= 1e-15 && leftCell.alfPigD > (1. - globals.localtiny) && superficialLiquidVelocity > 0 && rightSuperficialLiquidVelocity > 0) {
            setFlowPartitionTerms(state, cellIndex, 0., bif[cellIndex], 0);
        }
        if (fabs(superficialGasVelocity) <= 1e-15 && cell.alfL >= leftCell.alfL && leftSuperficialGasVelocity > 0) {
            bif[cellIndex] = 1;
        }
        if (fabs(superficialGasVelocity) <= 1e-15 && fabs(superficialLiquidVelocity) <= 1e-15 && cell.alf < leftCell.alf && cell.duto.teta > 0) {
            bif[cellIndex] = 1;
        }
    } else if (superficialLiquidVelocity >= 0 && leftCell.alfPigD >= 1. - 1 * globals.localtiny && (noLiquidSourceOnEitherFace)) {
        setFlowPartitionTerms(state, cellIndex, 0., bif[cellIndex], 0);
        clearDriftClosure(state, cellIndex);
        if (fabs(superficialLiquidVelocity) <= 1e-15 && cell.alfPigE < globals.localtiny && rightSuperficialLiquidVelocity < 0) { // ATTENTION!!!!!!!!!!!!!!! shouldn't this be two-phase, single-phase liquid only if the angle were negative?
            setFlowPartitionTerms(state, cellIndex, 1., bif[cellIndex], 0);
        } else if (fabs(rightSuperficialLiquidVelocity) < globals.localtiny * 1e-5) {
            if (fabs(superficialLiquidVelocity) < globals.localtiny * 1e-5 && cell.alfPigE < globals.localtiny && cell.fontemassGR >= globals.localtiny * 1e-5) { // ATTENTION!!!!!!!!!!!!!!!  this makes no sense here
                setFlowPartitionTerms(state, cellIndex, 1., bif[cellIndex], 0);
            }
            if (fabs(superficialLiquidVelocity) < globals.localtiny * 1e-5 && cell.alfPigE < globals.localtiny && cell.duto.teta >= 0) { // ATTENTION!!!!!!!!!!!!!!! change of 11/08/24, added
                bif[cellIndex] = 1;
            }
            if (fabs(superficialLiquidVelocity) < globals.localtiny * 1e-5 && cell.alfPigE < globals.localtiny && cell.duto.teta < 0) { // ATTENTION!!!!!!!!!!!!!!! change of 11/08/24, added
                setFlowPartitionTerms(state, cellIndex, 0., bif[cellIndex], 0);
            }
        } else if ((fabs(superficialLiquidVelocity) < 1e-15 && (rightSuperficialLiquidVelocity < 0 || cell.duto.teta > 0) // ATTENTION!!!!!!!!!!!!!!! change of 11/08/24: it was ||, changed to &&
                    && ((cell.alfPigE <= (1 - 10 * globals.localtiny + .0 * cell.alfPigER) &&
                         cell.alfPigER < 1 - 1 * globals.localtiny) ||
                        cell.alfPigE <= 0.7)))
            bif[cellIndex] = 1;

    } else if (superficialLiquidVelocity <= 0 && cell.alfPigE >= 1. - globals.localtiny && (noLiquidSourceOnEitherFace)) {
        setFlowPartitionTerms(state, cellIndex, 0., bif[cellIndex], 0);
        clearDriftClosure(state, cellIndex);
        if (fabs(superficialLiquidVelocity) <= globals.localtiny * 1e-5 && leftCell.alfPigD < globals.localtiny && leftSuperficialLiquidVelocity > 0) {
            setFlowPartitionTerms(state, cellIndex, 1., bif[cellIndex], 0);
        }

        if (fabs(leftSuperficialLiquidVelocity) < globals.localtiny * 1e-5) {
            if (fabs(superficialLiquidVelocity) < globals.localtiny * 1e-5 && leftCell.alfPigD < globals.localtiny && leftCell.fontemassGR >= globals.localtiny * 1e-5) {
                setFlowPartitionTerms(state, cellIndex, 1., bif[cellIndex], 0);
            }
        } else if ((cellIndex > 1 && fabs(superficialLiquidVelocity) < globals.localtiny * 1e-5 && cell.duto.teta < 0.95 * M_PI / 2. && leftSuperficialLiquidVelocity > 0 && ((leftCell.alfPigD <= 1.0 * state.cells[cellIndex - 2].alfPigD && state.cells[cellIndex - 2].alfPigD < 0.99) || leftCell.alfPigD < 0.7)) && superficialGasVelocity >= 0.)
            bif[cellIndex] = 1;
        else if ((cellIndex > 1 && fabs(superficialLiquidVelocity) < globals.localtiny * 1e-5 && cell.duto.teta >= 0.95 * M_PI / 2. && leftSuperficialLiquidVelocity > 0 && ((leftCell.alfPigD <= 1.0 * state.cells[cellIndex - 2].alfPigD && state.cells[cellIndex - 2].alfPigD < 0.99) || leftCell.alfPigD < 0.7)) && superficialGasVelocity >= 0)
            bif[cellIndex] = 1;
        else {
            if ((fabs(superficialLiquidVelocity) < globals.localtiny * 1e-5 && cell.duto.teta < 0.95 * M_PI / 2. && leftSuperficialLiquidVelocity > 0 && ((leftCell.alfPigD <= 1.0 * leftCell.alfL) || leftCell.alfPigD < 0.7)) && superficialGasVelocity >= 0)
                bif[cellIndex] = 1;
            else if ((fabs(superficialLiquidVelocity) < globals.localtiny * 1e-5 && cell.duto.teta >= 0.95 * M_PI / 2. && leftSuperficialLiquidVelocity > 0 && ((leftCell.alfPigD <= 1.0 * leftCell.alfL) || leftCell.alfPigD < 0.7)) && superficialGasVelocity >= 0)
                bif[cellIndex] = 1;
            else if ((fabs(superficialLiquidVelocity) < globals.localtiny * 1e-5 && cell.duto.teta >= 0.95 * M_PI / 2. && ((leftCell.alfPigD <= 1.0 * leftCell.alfL) || leftCell.alfPigD < 0.7)) && superficialGasVelocity > 0)
                bif[cellIndex] = 1;
        }
    }

    if (bif[cellIndex] == 1) {

        double c0;
        double ud;
        double meanVoidFraction;

        meanVoidFraction = leftCell.alfPigD;
        const SlugClosure slug = applySlugDriftClosure(
            state, cellIndex, superficialGasVelocity, rightGasDensity,
            rightLiquidDensity, meanVoidFraction);
        c0 = slug.c0;
        ud = slug.ud;
        if (fabs(cell.QG / (0.25 * M_PI * slug.meanDiameter * slug.meanDiameter * meanVoidFraction)) > 100. ||
            fabs(cell.QL / (0.25 * M_PI * slug.meanDiameter * slug.meanDiameter * (1. - meanVoidFraction))) > 100.) {
            c0 = 1.;
            ud = 0.;
        } else
            state.closureUpdater.instantaneous(cellIndex, c0, ud);
        cell.c0 = c0;
        cell.ud = ud;
        if (cellIndex == state.lastCell) {
            double numerator = (1. - meanVoidFraction * c0);
            double denominator = 1 + c0 * meanVoidFraction * ((gasDensity / liquidDensity) - 1.);
            cell.term1 = numerator / denominator;
            cell.term2 = (-flowArea * meanVoidFraction * gasDensity * ud) / denominator;
            double previousLiquidDriftFlux = (superficialGasVelocity - meanVoidFraction * ud) / (meanVoidFraction * c0) - superficialGasVelocity;
            double liquidDriftFlux = (superficialGasVelocity + superficialLiquidVelocity) * (1. - c0 * meanVoidFraction) - meanVoidFraction * ud;
            if ((liquidDriftFlux > 0. || previousLiquidDriftFlux > 0.) && leftCell.alfPigD > 1 - 1e-15) {
                cell.term1 = 0.;
                cell.term2 = 0.;
            }
            if ((liquidDriftFlux < 0. || previousLiquidDriftFlux < 0.) && cell.alfPigE > 1 - 1e-15) {
                cell.term1 = 0.;
                cell.term2 = 0.;
            }
            if (leftCell.acsr.tipo == kAccessoryChoke && leftCell.acsr.chk.AreaGarg <= (1e-3 + state.input.master1.razareaativ) * leftCell.duto.area) {
                cell.term1 = 0.;
                cell.term2 = 0.;
            }
        }
    }
}

void updateInteriorFlowPartitionCell(
    const ThermalState &state, int cellIndex, Vcr<int> &bif, Vcr<int> &valv) {
    Cel &cell = state.cells[cellIndex];
    Cel &leftCell = state.cells[cellIndex - 1];
    Cel &rightCell = state.cells[cellIndex + 1];
    const bool hasRightCell = cellIndex < state.lastCell;
    if (hasRightCell) {
        leftCell.alfR = rightCell.alfL = cell.alf;
        leftCell.betR = rightCell.betL = cell.bet;
    } else {
        leftCell.alfR = cell.alf;
        leftCell.betR = cell.bet;
    }
    valv[cellIndex] = 1;
    if (leftCell.acsr.tipo == kAccessoryChoke || leftCell.acsr.tipo == kAccessoryVolumetricPump) {
        if ((*cell.acsrL).tipo == kAccessoryChoke && (*cell.acsrL).chk.AreaGarg <= (1e-3 + state.input.master1.razareaativ) * leftCell.duto.area)
            valv[cellIndex] = 0;
        if ((*cell.acsrL).tipo == kAccessoryVolumetricPump && fabs((*cell.acsrL).bvol.freq) > 1)
            valv[cellIndex] = 0;
    }
    if (valv[cellIndex] == 1) {
        double betI = leftCell.betPigD;
        double liquidDensity;
        if (cell.QL < 0.) { // beta test
            betI = cell.betPigE;
            liquidDensity = (1 - betI) * cell.rpCi + betI * cell.rcCi;
        } else {
            betI = leftCell.betPigD;
            liquidDensity = (1 - betI) * leftCell.rpCi + betI * leftCell.rcCi;
        }
        const UpstreamFaceBasis upstream = upstreamFaceBasisOf(state, cellIndex);
        const double gasDensity = upstream.gasDensity;
        const double flowArea = upstream.flowArea;
        double superficialGasVelocity = cell.QG / (flowArea);
        double superficialLiquidVelocity = cell.QL / (flowArea);


        double leftFlowArea = cell.dutoL.area;
        double betIL = 0.;
        if (cellIndex < 2)
            betIL = leftCell.betL;
        else
            betIL = state.cells[cellIndex - 2].betPigD;
        if (leftCell.QL < 0.)
            betIL = leftCell.betPigE; // beta test
        // betIL = leftCell.betPigE;        // beta doubt
        double leftGasDensity = cell.rgLi;
        double leftLiquidDensity = (1 - betIL) * cell.rpLi + betIL * cell.rcLi;
        double leftSuperficialGasVelocity = (cell.ML - cell.MliqiniL) / (leftGasDensity * leftFlowArea);
        double leftSuperficialLiquidVelocity = (cell.MliqiniL) / (leftLiquidDensity * leftFlowArea);

        double rightFlowArea = cell.dutoR.area;
        double betIR = cell.betPigD;
        if (cell.QLR < 0.) { // beta test
            if (cellIndex > state.lastCell - 2)
                betIR = cell.betR;
            else
                betIR = rightCell.betPigE;
        }
        double rightGasDensity = cell.rgRi;
        double rightLiquidDensity = (1 - betIR) * cell.rpRi + betIR * cell.rcRi;
        double rightSuperficialGasVelocity = (cell.MR - cell.MliqiniR) / (rightGasDensity * rightFlowArea);
        double rightSuperficialLiquidVelocity = (cell.MliqiniR) / (rightLiquidDensity * rightFlowArea);

        selectAndApplyInteriorFlowRegime(
            state, cellIndex, bif, superficialGasVelocity, superficialLiquidVelocity, leftSuperficialGasVelocity, leftSuperficialLiquidVelocity, rightSuperficialGasVelocity, rightSuperficialLiquidVelocity,
            rightGasDensity, rightLiquidDensity, gasDensity, liquidDensity, flowArea);
    } else {
        cell.c0 = 1.;
        cell.ud = 0.;
        cell.term1 = 0.;
        cell.term2 = 0.;
        cell.term1L = leftCell.term1;
        cell.term2L = leftCell.term2;
        leftCell.term1R = cell.term1;
        leftCell.term2R = cell.term2;
    }
}

void updateOutletBoundaryFlowPartition(
    const ThermalState &state, int cellIndex, Vcr<int> &bif) {
    varGlob1D &globals = *state.globals;
    Cel &cell = state.cells[cellIndex];
    Cel &leftCell = state.cells[cellIndex - 1];
    const bool noLiquidSourceOnEitherFace =
        (cell.fontemassLL + cell.fontemassCL) <= globals.localtiny * 1e-5 && (cell.fontemassLR + cell.fontemassCR) <= globals.localtiny * 1e-5;
    state.cells[state.lastCell - 1].alfR = state.cells[state.lastCell].alf;
    state.cells[state.lastCell].alfR = state.cells[state.lastCell].alf;
    state.cells[state.lastCell - 1].betR = state.cells[state.lastCell].bet;
    state.cells[state.lastCell].betR = state.cells[state.lastCell].bet;

    double lengthRatio = cell.dxL / (cell.dx + cell.dxL);
    double meanPressure = cell.presaux;
    double meanTemperature = cell.temp * lengthRatio + leftCell.temp * (1. - lengthRatio);
    double betI = cell.betL;
    if (cell.QL < 0.)
        betI = cell.bet; // beta test
    // betI = cell.bet;            // beta doubt
    double gasDensity = cell.flui.MasEspGas(meanPressure, meanTemperature);
    double liquidDensity = (1 - betI) * cell.flui.MasEspLiq(meanPressure, meanTemperature) + betI * cell.fluicol.MasEspFlu(meanPressure, meanTemperature);
    double flowArea = cell.duto.area;
    if (cell.MC >= 0)
        flowArea = cell.dutoL.area;
    double superficialGasVelocity = cell.QG / (flowArea);
    double superficialLiquidVelocity = cell.QL / (flowArea);

    double leftFlowArea = cell.dutoL.area;
    double leftLengthRatio = leftCell.dxL / (leftCell.dx + leftCell.dxL);
    double leftMeanPressure = leftCell.presaux;
    double leftMeanTemperature = leftCell.temp * leftLengthRatio + leftCell.tempL * (1. - leftLengthRatio);
    double betIL = leftCell.betL;
    if (leftCell.QL < 0.)
        betIL = leftCell.bet; // beta test
    // betIL = leftCell.bet;            // beta doubt
    double leftLiquidDensity = (1 - betIL) * cell.flui.MasEspLiq(leftMeanPressure, leftMeanTemperature) + betIL * cell.fluicol.MasEspFlu(leftMeanPressure, leftMeanTemperature);
    double leftSuperficialLiquidVelocity = (cell.MliqiniL) / (leftLiquidDensity * leftFlowArea);

    double rightFlowArea = cell.dutoR.area;
    double rightLengthRatio = cell.dxR / (cell.dxR + cell.dx);
    double rightMeanPressure = cell.pres;
    double rightMeanTemperature = cell.temp * rightLengthRatio + cell.tempR * (1. - rightLengthRatio);
    double betIR = cell.bet;
    double rightGasDensity = cell.flui.MasEspGas(rightMeanPressure, rightMeanTemperature);
    double rightLiquidDensity = (1 - betIR) * cell.flui.MasEspLiq(rightMeanPressure, rightMeanTemperature) + betIR * cell.fluicol.MasEspFlu(rightMeanPressure, rightMeanTemperature);
    double rightSuperficialLiquidVelocity = (cell.MliqiniR) / (rightLiquidDensity * rightFlowArea);

    bif[cellIndex] = 1;

    if (cell.alfL <= globals.localtiny && cell.alf <= globals.localtiny && cell.fontemassGL <= 0 && cell.fontemassGR <= 0) {
        setFlowPartitionTerms(state, cellIndex, 1., bif[cellIndex], 0);
        cell.c0 = 1. - globals.localtiny;
        cell.ud = 0.;
        cell.arranjo = 0;
    } else if (cell.alfL >= (1. - globals.localtiny) && cell.alf >= (1. - globals.localtiny) && cell.fontemassLL <= 0 && cell.fontemassLR <= 0) {
        setFlowPartitionTerms(state, cellIndex, 0., bif[cellIndex], 0);
        cell.c0 = 1. - globals.localtiny;
        cell.ud = 0.;
        cell.arranjo = 0;
    } else if (superficialGasVelocity > 0 && cell.alfL <= globals.localtiny && (cell.fontemassGL <= 0. && cell.fontemassGR <= 0.)) {
        setFlowPartitionTerms(state, cellIndex, 1., bif[cellIndex], 0);
        cell.c0 = 1. - globals.localtiny;
        cell.ud = 0.;
        cell.arranjo = 0;
        if (fabs(superficialGasVelocity) <= 1e-15 && cell.alf > (1. - globals.localtiny) && superficialLiquidVelocity < 0 && leftSuperficialLiquidVelocity < 0) {
            setFlowPartitionTerms(state, cellIndex, 0., bif[cellIndex], 0);
        }
    } else if (superficialGasVelocity < 0 && cell.alf <= globals.localtiny && (cell.fontemassGL <= 0. && cell.fontemassGR <= 0.)) {
        setFlowPartitionTerms(state, cellIndex, 1., bif[cellIndex], 0);
        cell.c0 = 1. - globals.localtiny;
        cell.ud = 0.;
        cell.arranjo = 0;
        if (fabs(superficialGasVelocity) <= 1e-15 && cell.alfL > (1. - globals.localtiny) && superficialLiquidVelocity > 0 && rightSuperficialLiquidVelocity > 0) {
            setFlowPartitionTerms(state, cellIndex, 0., bif[cellIndex], 0);
        }
    } else if (superficialLiquidVelocity >= 0 && leftCell.alf >= 1. - globals.localtiny && (noLiquidSourceOnEitherFace)) {
        setFlowPartitionTerms(state, cellIndex, 0., bif[cellIndex], 0);
        clearDriftClosure(state, cellIndex);
        if (fabs(superficialLiquidVelocity) <= 1e-15 && cell.alf < globals.localtiny && rightSuperficialLiquidVelocity < 0) {
            setFlowPartitionTerms(state, cellIndex, 1., bif[cellIndex], 0);
        } else if (fabs(rightSuperficialLiquidVelocity) < globals.localtiny * 1e-5) {
            if (fabs(superficialLiquidVelocity) < globals.localtiny * 1e-5 && cell.alf < globals.localtiny && cell.fontemassGR >= globals.localtiny * 1e-5) {
                setFlowPartitionTerms(state, cellIndex, 1., bif[cellIndex], 0);
            }
        } else if (fabs(superficialLiquidVelocity) < 1e-15 && rightSuperficialLiquidVelocity < 0 && cell.alf <= (1 - 1 * globals.localtiny))
            bif[cellIndex] = 1;
    }

    else if (superficialLiquidVelocity <= 0 && cell.alf >= 1. - globals.localtiny && (noLiquidSourceOnEitherFace)) {
        setFlowPartitionTerms(state, cellIndex, 0., bif[cellIndex], 0);
        clearDriftClosure(state, cellIndex);
        if (fabs(superficialLiquidVelocity) <= globals.localtiny * 1e-5 && leftCell.alf < globals.localtiny && leftSuperficialLiquidVelocity > 0) {
            setFlowPartitionTerms(state, cellIndex, 1., bif[cellIndex], 0);
        }
        if (fabs(leftSuperficialLiquidVelocity) < globals.localtiny * 1e-5) {
            if (fabs(superficialLiquidVelocity) < globals.localtiny * 1e-5 && leftCell.alf < globals.localtiny && leftCell.fontemassGR >= globals.localtiny * 1e-5) {
                setFlowPartitionTerms(state, cellIndex, 1., bif[cellIndex], 0);
            }
        } else if (fabs(superficialLiquidVelocity) < globals.localtiny * 1e-5 && cell.duto.teta < 0.95 * M_PI / 2. && leftSuperficialLiquidVelocity > 0 && ((leftCell.alfPigD <= 1.0 * state.cells[cellIndex - 2].alfPigD && state.cells[cellIndex - 2].alfPigD < 0.99) || leftCell.alfPigD < 0.7))
            bif[cellIndex] = 1;
        else if (fabs(superficialLiquidVelocity) < globals.localtiny * 1e-5 && cell.duto.teta >= 0.95 * M_PI / 2. && leftSuperficialLiquidVelocity > 0 && ((leftCell.alfPigD <= 1.0 * state.cells[cellIndex - 2].alfPigD && state.cells[cellIndex - 2].alfPigD < 0.99) || leftCell.alfPigD < 0.7))
            bif[cellIndex] = 1;
        else {
            if (fabs(superficialLiquidVelocity) < globals.localtiny * 1e-5 && cell.duto.teta < 0.95 * M_PI / 2. && leftSuperficialLiquidVelocity > 0 && ((leftCell.alfPigD <= 1.0 * leftCell.alfL) || leftCell.alfPigD < 0.7))
                bif[cellIndex] = 1;
            else if (fabs(superficialLiquidVelocity) < globals.localtiny * 1e-5 && cell.duto.teta >= 0.95 * M_PI / 2. && leftSuperficialLiquidVelocity > 0 && ((leftCell.alfPigD <= 1.0 * leftCell.alfL) || leftCell.alfPigD < 0.7))
                bif[cellIndex] = 1;
        }
    }

    if (bif[cellIndex] == 1) {
        double meanVoidFraction;
        meanVoidFraction = cell.alfL;
        double c0 = 1.2;
        double meanDiameter = cell.duto.a;
        if (cell.MC >= 0)
            meanDiameter = cell.dutoL.a;
        double inclinationSign = 1.;
        if (cell.duto.teta < 0.)
            inclinationSign = 1.;
        double ud = inclinationSign * 0.32 * sqrt(kGravity * meanDiameter);
        if (fabs(rightGasDensity) / rightLiquidDensity > 0.9) {
            c0 = 1.;
            ud = 0.;
        }
        state.closureUpdater.instantaneous(cellIndex, c0, ud);
        if (state.input.escorregamentoCelulaContorno == 0) {
            c0 = 1.;
            ud = 0.;
        }
        cell.c0 = c0;
        cell.ud = ud;
        double numerator = (1. - meanVoidFraction * c0);
        double denominator = 1. + meanVoidFraction * (gasDensity / liquidDensity) * c0 - meanVoidFraction * c0;
        cell.term1 = numerator / denominator;
        cell.term2 = (-flowArea * meanVoidFraction * gasDensity * ud) / denominator;

        // test!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        // test!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    }

    cell.term1L = leftCell.term1;
    cell.term2L = leftCell.term2;
    leftCell.term1R = cell.term1;
    leftCell.term2R = cell.term2;
}

/// Blends the drift closure of a cell with a neighbour when the two disagree on
/// inclination or on flow pattern. `blendIndex` is the cell whose c0V/udV entry
/// the caller is filling.
void blendNeighbourDriftClosure(const ThermalState &state, int cellIndex,
                                int neighbourIndex, int blendIndex,
                                const Vcr<int> &bif,
                                Vcr<double> &c0V, Vcr<double> &udV) {
    double inclinationAngle;
    double leftInclinationAngle;
    inclinationAngle = state.cells[cellIndex].duto.teta;
    leftInclinationAngle = state.cells[neighbourIndex].duto.teta;
    if (((inclinationAngle != leftInclinationAngle) && bif[blendIndex] != 0) ||
        (state.cells[cellIndex].arranjo != state.cells[blendIndex].arranjo && bif[blendIndex] != 0)) {
        c0V[cellIndex] = (state.cells[cellIndex].dx * state.cells[cellIndex].c0 + state.cells[neighbourIndex].dx * state.cells[neighbourIndex].c0) / (state.cells[cellIndex].dx + state.cells[neighbourIndex].dx);
        if (inclinationAngle * leftInclinationAngle >= 0)
            udV[cellIndex] = (state.cells[cellIndex].dx * state.cells[cellIndex].ud + state.cells[neighbourIndex].dx * state.cells[neighbourIndex].ud) / (state.cells[cellIndex].dx + state.cells[neighbourIndex].dx);
    }
}

void finalizeFlowPartitionTerms(
    const ThermalState &state, const Vcr<int> &bif, const Vcr<int> &valv) {
    int isShutIn = 0;
    for (int cellIndex = 1; cellIndex < state.lastCell; cellIndex++) {
        Cel &cell = state.cells[cellIndex];
        if (cell.acsr.tipo == kAccessoryChoke && cell.acsr.chk.AreaGarg <= 1e-15 * cell.acsr.chk.AreaTub)
            isShutIn = 1;
        else if (state.surfaceChoke.AreaGarg <= 1.e-15 * state.surfaceChoke.AreaTub)
            isShutIn = 1;
    }

    Vcr<double> c0V(state.lastCell, 0.);
    Vcr<double> udV(state.lastCell, 0.);
    for (int cellIndex = 1; cellIndex < state.lastCell; cellIndex++) {
        Cel &cell = state.cells[cellIndex];
        c0V[cellIndex] = cell.c0;
        udV[cellIndex] = cell.ud;
        if (bif[cellIndex] == 1 && valv[cellIndex] == 1 && cellIndex > 2) {
            int neighborIndex = cellIndex - 1;
            int secondNeighborIndex = cellIndex - 2;
            if (cell.QG < 0) {
                neighborIndex = cellIndex + 1;
                secondNeighborIndex = cellIndex;
            }
            if ((state.cells[neighborIndex].acsr.tipo == kAccessoryNone && (state.cells[secondNeighborIndex].acsr.tipo != kAccessoryChoke || state.cells[secondNeighborIndex].acsr.chk.AreaGarg > (1e-3))) &&
                (cell.arranjo != state.cells[neighborIndex].arranjo && bif[neighborIndex] != 0)) {
                c0V[cellIndex] = (cell.dx * cell.c0 + state.cells[neighborIndex].dx * state.cells[neighborIndex].c0) / (cell.dx + state.cells[neighborIndex].dx);
                if (cell.duto.teta * state.cells[neighborIndex].duto.teta >= 0)
                    udV[cellIndex] = (cell.dx * cell.ud + state.cells[neighborIndex].dx * state.cells[neighborIndex].ud) / (cell.dx + state.cells[neighborIndex].dx);
            } else if (cellIndex > 2) {
                if ((state.cells[cellIndex - 1].acsr.tipo == kAccessoryNone && (state.cells[cellIndex - 2].acsr.tipo != kAccessoryChoke || state.cells[cellIndex - 2].acsr.chk.AreaGarg > (1e-3))) &&
                    cell.QG >= 0) {
                blendNeighbourDriftClosure(state, cellIndex, cellIndex - 1,
                                           neighborIndex, bif, c0V, udV);
                } else if ((state.cells[cellIndex + 1].acsr.tipo == kAccessoryNone && (cell.acsr.tipo != kAccessoryChoke || cell.acsr.chk.AreaGarg > (1e-3))) &&
                           cell.QG < 0) {
                blendNeighbourDriftClosure(state, cellIndex, cellIndex + 1,
                                           neighborIndex, bif, c0V, udV);
                }
            }
        }
    }

    for (int cellIndex = 1; cellIndex < state.lastCell; cellIndex++) {
        Cel &cell = state.cells[cellIndex];
        Cel &leftCell = state.cells[cellIndex - 1];
        if (bif[cellIndex] == 1 && valv[cellIndex] == 1) {
            double betI = leftCell.betPigD;
            double liquidDensity;

            if (cell.QL < 0.) {
                betI = cell.betPigE;
                liquidDensity = (1 - betI) * cell.rpCi + betI * cell.rcCi;
            } else {
                betI = leftCell.betPigD;
                liquidDensity = (1 - betI) * leftCell.rpCi + betI * leftCell.rcCi;
            }

            const UpstreamFaceBasis upstream = upstreamFaceBasisOf(state, cellIndex);
            const double gasDensity = upstream.gasDensity;
            const double flowArea = upstream.flowArea;
            double superficialGasVelocity = cell.QG / (flowArea);
            double superficialLiquidVelocity = cell.QL / (flowArea);

            double meanVoidFraction;
            meanVoidFraction = leftCell.alfPigD;
            if (superficialGasVelocity < 0)
                meanVoidFraction = cell.alfPigE;
            double numerator = (1. - meanVoidFraction * cell.c0);
            double denominator = 1 + cell.c0 * meanVoidFraction * ((gasDensity / liquidDensity) - 1.);
            cell.term1 = numerator / denominator;
            cell.term2 = (-flowArea * meanVoidFraction * gasDensity * cell.ud) / denominator;
            double previousLiquidDriftFlux = (superficialGasVelocity - meanVoidFraction * cell.ud) / (meanVoidFraction * cell.c0) - superficialGasVelocity;
            double liquidDriftFlux = (superficialGasVelocity + superficialLiquidVelocity) * (1. - cell.c0 * meanVoidFraction) - meanVoidFraction * cell.ud;
            double candidateLiquidMass = cell.term1 * cell.MC + cell.term2;
            double candidateGasMass = (1 - cell.term1) * cell.MC - cell.term2;
            if ((liquidDriftFlux > 0. || previousLiquidDriftFlux > 0.) && leftCell.alfPigD > 1 - 1e-15) {
                cell.term1 = 0.;
                cell.term2 = 0.;
            }
            if ((liquidDriftFlux < 0. || previousLiquidDriftFlux < 0.) && cell.alfPigE > 1 - 1e-15) {
                cell.term1 = 0.;
                cell.term2 = 0.;
            }
            if (leftCell.acsr.tipo == kAccessoryChoke && leftCell.acsr.chk.AreaGarg <= (1e-3 + state.input.master1.razareaativ) * leftCell.duto.area) {
                cell.term1 = 0.;
                cell.term2 = 0.;
            }
            if ((liquidDriftFlux < 0. || previousLiquidDriftFlux < 0.) && ((fabs(leftCell.QG / (flowArea)) + fabs(leftCell.QL / (flowArea))) < 0.1) &&
                (isShutIn == 1 && state.input.modoSegrega == 1) && (cell.duto.teta > 0 && leftCell.duto.teta <= 0) &&
                (candidateLiquidMass > 0 && cell.term2 > 0) &&
                ((1. - cell.alfL) /**fabs(sin(state.cells[cellIndex-1].duto.teta))*/ < (1. - cell.alf) /**fabs(sin(cell.duto.teta))*/)) {
                cell.term2 = 0.;
            }
            if (cell.duto.teta > 0 && candidateGasMass < 0 && ((isShutIn == 1 && state.input.modoSegrega == 1)) && cell.alfPigE > 1 - 1e-15) {
                cell.term1 = 0.;
                cell.term2 = 0.;
            }
        }
    }
    for (int cellIndex = 1; cellIndex <= state.lastCell; cellIndex++) {
        Cel &cell = state.cells[cellIndex];
        Cel &leftCell = state.cells[cellIndex - 1];
        cell.term1L = leftCell.term1;
        cell.term2L = leftCell.term2;
        leftCell.term1R = cell.term1;
        leftCell.term2R = cell.term2;
    }
}

/// Collapses a cell onto the no-slip closure. `term1` is 1. for single-phase
/// liquid, 0. for single-phase gas.
///
/// Do not fold the two products to 1 and 0: 0 * x is NaN for non-finite x.
void applyNoSlipClosure(const ThermalState &state, int cellIndex, double term1,
                        int &branchFlag, double distributionCoefficient,
                        double driftVelocity) {
    Cel &cell = state.cells[cellIndex];
    cell.term1 = term1;
    cell.term2 = 0.;
    branchFlag = 0;
    cell.c0 = 1 + 0 * distributionCoefficient;
    cell.ud = 0 * driftVelocity;
    cell.arranjo = 0;
}

/// The two inlet regime selectors differ in three places only; each policy
/// supplies those three and shares everything else.
struct InletBoundaryRegime {
    static void applyInletOverrides(const ThermalState &state, int cellIndex,
                                    int &branchFlag, double superficialGasVelocity,
                                    double superficialLiquidVelocity) {
        varGlob1D &globals = *state.globals;
        if (superficialLiquidVelocity >= 0 && state.inletVoidFraction > 1 - globals.localtiny) {
            setFlowPartitionTerms(state, cellIndex, 0., branchFlag, 0);
        }
        if (superficialGasVelocity >= 0 && state.inletVoidFraction < globals.localtiny) {
            setFlowPartitionTerms(state, cellIndex, 1., branchFlag, 0);
        }
    }

    static void updateClosure(const ThermalState &state, int cellIndex, double &c0,
                              double &ud) {
        state.closureUpdater.initialization(cellIndex, c0, ud);
        state.cells[cellIndex].c0 = c0;
        state.cells[cellIndex].ud = ud;
    }
};

struct BufferedInletRegime {
    static void applyInletOverrides(const ThermalState &, int, int &, double, double) {}

    static void updateClosure(const ThermalState &state, int cellIndex, double &c0,
                              double &ud) {
        state.closureUpdater.bufferedInitialization(cellIndex, c0, ud);
    }
};

/// Chooses the flow regime at the inlet boundary and applies its closure.
/// `branchFlag` is the caller's bif: an array element for the instantaneous
/// pass, a local for the buffered one.
template <typename Regime>
void selectAndApplyInletRegime(
    const ThermalState &state, int cellIndex, int &branchFlag, double distributionCoefficient,
    double driftVelocity, double superficialGasVelocity, double superficialLiquidVelocity, double rightSuperficialLiquidVelocity, double rightGasDensity,
    double rightLiquidDensity, double gasDensity, double liquidDensity, double flowArea) {
    varGlob1D &globals = *state.globals;
    Cel &cell = state.cells[cellIndex];
    const bool noGasSourceOnEitherFace =
        cell.fontemassGL <= globals.localtiny * 1e-5 && cell.fontemassGR <= globals.localtiny * 1e-5;
    const bool noLiquidSourceOnEitherFace =
        (cell.fontemassLL + cell.fontemassCL) <= globals.localtiny * 1e-5 && (cell.fontemassLR + cell.fontemassCR) <= globals.localtiny * 1e-5;
    branchFlag = 1;


    if (state.inletVoidFraction < globals.localtiny && cell.alfPigE <= globals.localtiny && noGasSourceOnEitherFace) {
        applyNoSlipClosure(state, cellIndex, 1., branchFlag, distributionCoefficient,
                           driftVelocity);
    } else if (state.inletVoidFraction >= (1. - globals.localtiny) && cell.alfPigE >= (1. - globals.localtiny) && (noLiquidSourceOnEitherFace)) {
        applyNoSlipClosure(state, cellIndex, 0., branchFlag, distributionCoefficient,
                           driftVelocity);
    } else if (superficialGasVelocity >= 0 && state.inletVoidFraction <= globals.localtiny && (noGasSourceOnEitherFace)) {
        applyNoSlipClosure(state, cellIndex, 1., branchFlag, distributionCoefficient,
                           driftVelocity);
        if (fabs(superficialGasVelocity) <= 1e-15 && cell.alfPigE > (1. - globals.localtiny) && superficialLiquidVelocity < 0 && cell.duto.teta > 0) {
            setFlowPartitionTerms(state, cellIndex, 0., branchFlag, 0);
        }
        if (fabs(superficialGasVelocity) <= 1e-15 && cell.alfPigE > (1. - 1 * globals.localtiny) && cell.duto.teta < 0 && superficialLiquidVelocity > 0) {
            setFlowPartitionTerms(state, cellIndex, 0., branchFlag, 1);
        }
        if (fabs(superficialGasVelocity) <= 1e-15 && cell.alf >= state.inletVoidFraction && superficialLiquidVelocity < 0) {
            branchFlag = 1;
        }
    } else if (superficialGasVelocity <= 0 && cell.alfPigE <= globals.localtiny && (noGasSourceOnEitherFace)) {
        applyNoSlipClosure(state, cellIndex, 1., branchFlag, distributionCoefficient,
                           driftVelocity);
        if (fabs(superficialGasVelocity) <= 1e-15 && state.inletVoidFraction > (1. - globals.localtiny) && superficialLiquidVelocity > 0) {
            setFlowPartitionTerms(state, cellIndex, 0., branchFlag, 0);
        }
        if (fabs(superficialGasVelocity) <= 1e-15 && state.inletVoidFraction > globals.localtiny && superficialLiquidVelocity > 0) {
            branchFlag = 1;
        }
    } else if (superficialLiquidVelocity >= 0 && state.inletVoidFraction >= 1. - globals.localtiny && (noLiquidSourceOnEitherFace)) {
        applyNoSlipClosure(state, cellIndex, 0., branchFlag, distributionCoefficient,
                           driftVelocity);
        if (fabs(superficialLiquidVelocity) <= 1e-15 && cell.alfPigE < globals.localtiny && rightSuperficialLiquidVelocity < 0) {
            setFlowPartitionTerms(state, cellIndex, 1., branchFlag, 0);
        } else if (fabs(rightSuperficialLiquidVelocity) < globals.localtiny * 1e-5) {
            if (fabs(superficialLiquidVelocity) < globals.localtiny * 1e-5 && cell.alfPigE < globals.localtiny && cell.fontemassGR >= globals.localtiny * 1e-5) {
                setFlowPartitionTerms(state, cellIndex, 1., branchFlag, 0);
            }
        } else if (fabs(superficialLiquidVelocity) < 1e-15 && rightSuperficialLiquidVelocity < 0 && ((cell.alfPigE <= (1 - 1 * globals.localtiny + .0 * cell.alfPigER) && cell.alfPigER < 1 - 1 * globals.localtiny) || cell.alfPigE <= 0.7))
            branchFlag = 1;

    } else if (superficialLiquidVelocity <= 0 && cell.alfPigE >= 1. - globals.localtiny && (noLiquidSourceOnEitherFace)) {
        applyNoSlipClosure(state, cellIndex, 0., branchFlag, distributionCoefficient,
                           driftVelocity);
        if (fabs(superficialLiquidVelocity) <= globals.localtiny * 1e-5 && state.inletVoidFraction < globals.localtiny && rightSuperficialLiquidVelocity < 0) {
            setFlowPartitionTerms(state, cellIndex, 1., branchFlag, 0);
        }

        if (fabs(rightSuperficialLiquidVelocity) < globals.localtiny * 1e-5) {
            if (fabs(superficialLiquidVelocity) < globals.localtiny * 1e-5 && state.inletVoidFraction < globals.localtiny) {
                setFlowPartitionTerms(state, cellIndex, 1., branchFlag, 0);
            }
        } else {
            if (fabs(superficialLiquidVelocity) < globals.localtiny * 1e-5 && cell.duto.teta < 0.95 * M_PI / 2. && superficialGasVelocity < 0 && (state.inletVoidFraction < 0.7))
                branchFlag = 1;
            else if (fabs(superficialLiquidVelocity) < globals.localtiny * 1e-5 && cell.duto.teta >= 0.95 * M_PI / 2. && superficialGasVelocity < 0 && (state.inletVoidFraction < 0.7))
                branchFlag = 1;
        }
    }
    if (superficialLiquidVelocity > 0 && fabs(superficialGasVelocity) <= globals.localtiny * 1e-5 && state.inletVoidFraction > globals.localtiny && state.inletVoidFraction < 1 - globals.localtiny)
        branchFlag = 1;
    if (superficialGasVelocity > 0 && fabs(superficialLiquidVelocity) <= globals.localtiny * 1e-5 && state.inletVoidFraction > globals.localtiny && state.inletVoidFraction < 1 - globals.localtiny)
        branchFlag = 1;
    Regime::applyInletOverrides(state, cellIndex, branchFlag,
                                superficialGasVelocity, superficialLiquidVelocity);

    if (branchFlag == 1) {

        double c0;
        double ud;
        double meanVoidFraction;

        meanVoidFraction = state.inletVoidFraction;
        const SlugClosure slug = applySlugDriftClosure(
            state, cellIndex, superficialGasVelocity, rightGasDensity,
            rightLiquidDensity, meanVoidFraction);
        c0 = slug.c0;
        ud = slug.ud;
        Regime::updateClosure(state, cellIndex, c0, ud);
        double numerator = (1. - meanVoidFraction * c0);
        double denominator = 1 + c0 * meanVoidFraction * ((gasDensity / liquidDensity) - 1.);
        cell.term1 = numerator / denominator;
        cell.term2 = (-flowArea * meanVoidFraction * gasDensity * ud) / denominator;
    }
}

void selectAndApplyInletBoundaryFlowRegime(
    const ThermalState &state, int cellIndex, Vcr<int> &bif, double distributionCoefficient,
    double driftVelocity, double superficialGasVelocity, double superficialLiquidVelocity, double rightSuperficialLiquidVelocity, double rightGasDensity,
    double rightLiquidDensity, double gasDensity, double liquidDensity, double flowArea) {
    selectAndApplyInletRegime<InletBoundaryRegime>(
        state, cellIndex, bif[cellIndex], distributionCoefficient, driftVelocity,
        superficialGasVelocity, superficialLiquidVelocity, rightSuperficialLiquidVelocity,
        rightGasDensity, rightLiquidDensity, gasDensity, liquidDensity, flowArea);
}


void selectAndApplyBufferedOutletFlowRegime(
    const ThermalState &state, int cellIndex, double distributionCoefficient, double driftVelocity,
    double superficialGasVelocity, double superficialLiquidVelocity, double leftSuperficialGasVelocity, double leftSuperficialLiquidVelocity, double rightSuperficialGasVelocity,
    double rightSuperficialLiquidVelocity, double rightGasDensity, double rightLiquidDensity, double gasDensity, double liquidDensity,
    double flowArea) {
    varGlob1D &globals = *state.globals;
    Cel &cell = state.cells[cellIndex];
    Cel &leftCell = state.cells[cellIndex - 1];
    Cel &rightCell = state.cells[cellIndex + 1];
    const bool noGasSourceOnEitherFace =
        cell.fontemassGL <= globals.localtiny * 1e-5 && cell.fontemassGR <= globals.localtiny * 1e-5;
    const bool noLiquidSourceOnEitherFace =
        (cell.fontemassLL + cell.fontemassCL) <= globals.localtiny * 1e-5 && (cell.fontemassLR + cell.fontemassCR) <= globals.localtiny * 1e-5;
    int branchFlag = 1;


    if (leftCell.alfPigD <= globals.localtiny && cell.alfPigE <= globals.localtiny && noGasSourceOnEitherFace) {
        applyNoSlipClosure(state, cellIndex, 1., branchFlag, distributionCoefficient,
                           driftVelocity);
    } else if (leftCell.alfPigD >= (1. - globals.localtiny) && cell.alfPigE >= (1. - globals.localtiny) && (noLiquidSourceOnEitherFace)) {
        applyNoSlipClosure(state, cellIndex, 0., branchFlag, distributionCoefficient,
                           driftVelocity);
    } else if (leftCell.acsr.tipo == kAccessoryChoke && leftCell.acsr.chk.AreaGarg <= (1e-3)) {
        applyNoSlipClosure(state, cellIndex, 0., branchFlag, distributionCoefficient,
                           driftVelocity);
    }

    else if (superficialGasVelocity >= 0 && leftCell.alfPigD <= globals.localtiny && (noGasSourceOnEitherFace)) {
        applyNoSlipClosure(state, cellIndex, 1., branchFlag, distributionCoefficient,
                           driftVelocity);
        if (fabs(superficialGasVelocity) <= 1e-15 && cell.alfPigE > (1. - globals.localtiny) && superficialLiquidVelocity < 0 && leftSuperficialLiquidVelocity < 0 && cell.duto.teta > 0) {
            setFlowPartitionTerms(state, cellIndex, 0., branchFlag, 0);
        }
        if (fabs(superficialGasVelocity) <= 1e-15 && cell.alfPigE > (1. - 10 * globals.localtiny) && cell.duto.teta < 0) {
            setFlowPartitionTerms(state, cellIndex, 0., branchFlag, 1);
        }
        if (fabs(superficialGasVelocity) <= 1e-15 && cell.alf >= rightCell.alf && rightSuperficialGasVelocity < 0) {
            branchFlag = 1;
        }
    } else if (superficialGasVelocity <= 0 && cell.alfPigE <= globals.localtiny && (noGasSourceOnEitherFace)) {
        applyNoSlipClosure(state, cellIndex, 1., branchFlag, distributionCoefficient,
                           driftVelocity);
        if (fabs(superficialGasVelocity) <= 1e-15 && leftCell.alfPigD > (1. - globals.localtiny) && superficialLiquidVelocity > 0 && rightSuperficialLiquidVelocity > 0) {
            setFlowPartitionTerms(state, cellIndex, 0., branchFlag, 0);
        }
        if (fabs(superficialGasVelocity) <= 1e-15 && cell.alfL >= leftCell.alfL && leftSuperficialGasVelocity > 0) {
            branchFlag = 1;
        }
    } else if (superficialLiquidVelocity >= 0 && leftCell.alfPigD >= 1. - globals.localtiny && (noLiquidSourceOnEitherFace)) {
        applyNoSlipClosure(state, cellIndex, 0., branchFlag, distributionCoefficient,
                           driftVelocity);
        if (fabs(superficialLiquidVelocity) <= 1e-15 && cell.alfPigE < globals.localtiny && rightSuperficialLiquidVelocity < 0) {
            setFlowPartitionTerms(state, cellIndex, 1., branchFlag, 0);
        } else if (fabs(rightSuperficialLiquidVelocity) < globals.localtiny * 1e-5) {
            if (fabs(superficialLiquidVelocity) < globals.localtiny * 1e-5 && cell.alfPigE < globals.localtiny && cell.fontemassGR >= globals.localtiny * 1e-5) {
                setFlowPartitionTerms(state, cellIndex, 1., branchFlag, 0);
            }
        } else if (fabs(superficialLiquidVelocity) < 1e-15 && rightSuperficialLiquidVelocity < 0 && ((cell.alfPigE <= (1 - 1 * globals.localtiny + .0 * cell.alfPigER) && cell.alfPigER < 1 - 1 * globals.localtiny) || cell.alfPigE <= 0.7))
            branchFlag = 1;

    } else if (superficialLiquidVelocity <= 0 && cell.alfPigE >= 1. - globals.localtiny && (noLiquidSourceOnEitherFace)) {
        applyNoSlipClosure(state, cellIndex, 0., branchFlag, distributionCoefficient,
                           driftVelocity);
        if (fabs(superficialLiquidVelocity) <= globals.localtiny * 1e-5 && leftCell.alfPigD < globals.localtiny && leftSuperficialLiquidVelocity > 0) {
            setFlowPartitionTerms(state, cellIndex, 1., branchFlag, 0);
        }

        if (fabs(leftSuperficialLiquidVelocity) < globals.localtiny * 1e-5) {
            if (fabs(superficialLiquidVelocity) < globals.localtiny * 1e-5 && leftCell.alfPigD < globals.localtiny && leftCell.fontemassGR >= globals.localtiny * 1e-5) {
                setFlowPartitionTerms(state, cellIndex, 1., branchFlag, 0);
            }
        } else if (cellIndex > 1 && fabs(superficialLiquidVelocity) < globals.localtiny * 1e-5 && cell.duto.teta < 0.95 * M_PI / 2. && leftSuperficialLiquidVelocity > 0 && ((leftCell.alfPigD <= 1.0 * state.cells[cellIndex - 2].alfPigD && state.cells[cellIndex - 2].alfPigD < 0.99) || leftCell.alfPigD < 0.7))
            branchFlag = 1;
        else if (cellIndex > 1 && fabs(superficialLiquidVelocity) < globals.localtiny * 1e-5 && cell.duto.teta >= 0.95 * M_PI / 2. && leftSuperficialLiquidVelocity > 0 && ((leftCell.alfPigD <= 1.0 * state.cells[cellIndex - 2].alfPigD && state.cells[cellIndex - 2].alfPigD < 0.99) || leftCell.alfPigD < 0.7))
            branchFlag = 1;
        else {
            if (fabs(superficialLiquidVelocity) < globals.localtiny * 1e-5 && cell.duto.teta < 0.95 * M_PI / 2. && leftSuperficialLiquidVelocity > 0 && ((leftCell.alfPigD <= 1.0 * leftCell.alfL) || leftCell.alfPigD < 0.7))
                branchFlag = 1;
            else if (fabs(superficialLiquidVelocity) < globals.localtiny * 1e-5 && cell.duto.teta >= 0.95 * M_PI / 2. && leftSuperficialLiquidVelocity > 0 && ((leftCell.alfPigD <= 1.0 * leftCell.alfL) || leftCell.alfPigD < 0.7))
                branchFlag = 1;
        }
    }

    if (superficialGasVelocity < 0 && cell.alf > globals.localtiny && cell.alf < 1. - globals.localtiny)
        branchFlag = 1;

    if (branchFlag == 1) {

        double c0;
        double ud;
        double meanVoidFraction;

        meanVoidFraction = leftCell.alfPigD;
        const SlugClosure slug = applySlugDriftClosure(
            state, cellIndex, superficialGasVelocity, rightGasDensity,
            rightLiquidDensity, meanVoidFraction);
        c0 = slug.c0;
        ud = slug.ud;
        state.closureUpdater.buffered(cellIndex, c0, ud);
        if (state.input.escorregamentoCelulaContorno == 0) {
            c0 = 1.;
            ud = 0.;
        }
        double numerator = (1. - meanVoidFraction * c0);
        double denominator = 1 + c0 * meanVoidFraction * ((gasDensity / liquidDensity) - 1.);
        cell.term1 = numerator / denominator;
        cell.term2 = (-flowArea * meanVoidFraction * gasDensity * ud) / denominator;
        if (leftCell.acsr.tipo == kAccessoryChoke && leftCell.acsr.chk.AreaGarg <= (1e-3 + state.input.master1.razareaativ) * leftCell.duto.area) {
            cell.term1 = 0.;
            cell.term2 = 0.;
        }
    }
}

void selectAndApplyBufferedInletFlowRegime(
    const ThermalState &state, int cellIndex, double distributionCoefficient, double driftVelocity,
    double superficialGasVelocity, double superficialLiquidVelocity, double rightSuperficialLiquidVelocity, double rightGasDensity, double rightLiquidDensity,
    double gasDensity, double liquidDensity, double flowArea) {
    int branchFlag = 1;
    selectAndApplyInletRegime<BufferedInletRegime>(
        state, cellIndex, branchFlag, distributionCoefficient, driftVelocity,
        superficialGasVelocity, superficialLiquidVelocity, rightSuperficialLiquidVelocity,
        rightGasDensity, rightLiquidDensity, gasDensity, liquidDensity, flowArea);
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
        Cel &cell = state.cells[cellIndex];
        cell.c0ini = cell.c0;
        cell.udini = cell.ud;
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

                double meanPressure;
                meanPressure = state.inletPressure;

                double meanTemperature;
                if (state.cells[0].QL < 0.)
                    meanTemperature = cell.temp;
                else
                    meanTemperature = state.inletTemperature;

                double gasDensity = state.cells[0].flui.MasEspGas(state.inletPressure, state.inletTemperature);
                double liquidDensity = state.cells[0].flui.MasEspLiq(state.inletPressure, state.inletTemperature);
                double rcis = state.cells[0].fluicol.MasEspFlu(state.inletPressure, state.inletTemperature);

                double liquidMixtureDensity = state.inletCompletionFraction * rcis + (1 - state.inletCompletionFraction) * liquidDensity;
                state.inletVoidFraction = (-state.inletQuality * liquidMixtureDensity / (state.inletQuality * gasDensity - gasDensity - state.inletQuality * liquidMixtureDensity)) / (state.cells[0].c0);

                double betI;
                double surfaceTension;
                if (cell.QL < 0.) { // beta test
                    betI = cell.betPigE;
                    liquidDensity = (1 - betI) * cell.flui.MasEspLiq(meanPressure, meanTemperature) + betI * cell.fluicol.MasEspFlu(meanPressure, meanTemperature);
                    surfaceTension = (1 - betI) * cell.flui.TensSuper(meanPressure, meanTemperature) + betI * cell.fluicol.TensSuper(meanPressure, meanTemperature);
                } else {
                    betI = state.inletCompletionFraction;
                    liquidDensity = (1 - betI) * (*cell.fluiL).MasEspLiq(meanPressure, meanTemperature) + betI * cell.fluicol.MasEspFlu(meanPressure, meanTemperature);
                    surfaceTension = (1 - betI) * (*cell.fluiL).TensSuper(meanPressure, meanTemperature) + betI * cell.fluicol.TensSuper(meanPressure, meanTemperature);
                }

                double flowArea;
                if (cell.QG >= 0) {
                    flowArea = cell.duto.area;
                    gasDensity = (*cell.fluiL).MasEspGas(meanPressure, meanTemperature);
                } else {
                    gasDensity = cell.flui.MasEspGas(meanPressure, meanTemperature);
                    flowArea = cell.duto.area;
                }
                double superficialGasVelocity = cell.QG / (flowArea);
                double superficialLiquidVelocity = cell.QL / (flowArea);

                double inclinationAngle = cell.duto.teta;
                double inclinationSign = 1.;
                if (inclinationAngle < 0.)
                    inclinationSign = -1.;
                double c0 = 2.;
                double ud = inclinationSign * 0.0246 * cos(inclinationAngle) + 1.606 * pow(kGravity * surfaceTension * (liquidDensity - gasDensity) / (liquidDensity * liquidDensity), 0.25) * sin(inclinationAngle);

                double rightFlowArea = cell.dutoR.area;
                double rightLengthRatio = cell.dx / (cell.dxR + cell.dx);
                double rightMeanPressure = cell.presauxR;
                double rightMeanTemperature = cell.temp * rightLengthRatio + cell.tempL * (1. - rightLengthRatio);
                double betIR = cell.betPigD;
                if (cell.QLR < 0.) // beta test
                    betIR = state.cells[cellIndex + 1].betPigE;

                double rightGasDensity = cell.flui.MasEspGas(rightMeanPressure, rightMeanTemperature);
                double rightLiquidDensity = (1 - betIR) * cell.flui.MasEspLiq(rightMeanPressure, rightMeanTemperature) + betIR * cell.fluicol.MasEspFlu(rightMeanPressure, rightMeanTemperature);
                double rightSuperficialLiquidVelocity = (cell.MliqiniR) / (rightLiquidDensity * rightFlowArea);

                selectAndApplyInletBoundaryFlowRegime(
                    state, cellIndex, bifurcationFlags, c0, ud, superficialGasVelocity, superficialLiquidVelocity, rightSuperficialLiquidVelocity,
                    rightGasDensity, rightLiquidDensity, gasDensity, liquidDensity, flowArea);
                state.cells[1].term1L = cell.term1;
                state.cells[1].term2L = cell.term2;

                state.cells[1].alfL = state.cells[0].alf;
                state.cells[0].alfL = state.inletVoidFraction;
                state.cells[1].betL = state.cells[0].bet;
                state.cells[0].betL = state.inletCompletionFraction;
            }

        } else if (aflu == 0) {
            updateOutletBoundaryFlowPartition(state, cellIndex, bifurcationFlags);
        }
    }

    finalizeFlowPartitionTerms(state, bifurcationFlags, valveFlags);
}
void updateOutletFlowPartitionTerms(const ThermalState &state) {

    int cellIndex = state.lastCell;
    Cel &cell = state.cells[cellIndex];
    Cel &leftCell = state.cells[cellIndex - 1];
    const bool bufferedGasMassIsPositive = (cell.MCBuf - cell.MliqiniBuf) >= 0;

    double meanPressure = cell.presBuf;
    double meanTemperature = state.gasSurfaceTemperature;
    meanTemperature = leftCell.temp;
    if (cell.VTemper < 0.)
        meanTemperature = state.gasSurfaceTemperature;
    double betI = leftCell.betPigD;
    double liquidDensity;
    double surfaceTension;
    if (cell.MliqiniBuf < 0.) { // beta test
        betI = cell.betPigE;
        liquidDensity = (1 - betI) * cell.flui.MasEspLiq(meanPressure, meanTemperature) + betI * cell.fluicol.MasEspFlu(meanPressure, meanTemperature);
        surfaceTension = (1 - betI) * cell.flui.TensSuper(meanPressure, meanTemperature) + betI * cell.fluicol.TensSuper(meanPressure, meanTemperature);
    } else {
        betI = leftCell.betPigD;
        liquidDensity = (1 - betI) * leftCell.flui.MasEspLiq(meanPressure, meanTemperature) + betI * leftCell.fluicol.MasEspFlu(meanPressure, meanTemperature);
        surfaceTension = (1 - betI) * leftCell.flui.TensSuper(meanPressure, meanTemperature) + betI * leftCell.fluicol.TensSuper(meanPressure, meanTemperature);
    }
    double gasDensity;
    double flowArea;
    if (bufferedGasMassIsPositive) {
        flowArea = cell.dutoL.area;
        gasDensity = leftCell.flui.MasEspGas(meanPressure, meanTemperature);
    } else {
        gasDensity = cell.flui.MasEspGas(meanPressure, meanTemperature);
        flowArea = cell.duto.area;
    }
    double superficialGasVelocity = (cell.MCBuf - cell.MliqiniBuf) / (gasDensity * flowArea);
    double superficialLiquidVelocity = (cell.MliqiniBuf) / (liquidDensity * flowArea);

    double inclinationAngle = cell.duto.teta;
    if (cellIndex >= 2) {
        if (state.cells[cellIndex - 2].acsr.tipo == kAccessoryChoke && state.cells[cellIndex - 2].acsr.chk.AreaGarg <= (1e-3)) {
            if (bufferedGasMassIsPositive)
                inclinationAngle = cell.duto.teta;
            else
                inclinationAngle = cell.dutoR.teta;
        } else {
            if (bufferedGasMassIsPositive)
                inclinationAngle = cell.dutoL.teta;
            else
                inclinationAngle = cell.duto.teta;
        }
    }
    double inclinationSign = 1.;
    if (inclinationAngle < 0.)
        inclinationSign = -1.;
    double c0 = 2.;
    double ud = inclinationSign * 0.0246 * cos(inclinationAngle) + 1.606 * pow(kGravity * surfaceTension * (liquidDensity - gasDensity) / (liquidDensity * liquidDensity), 0.25) * sin(inclinationAngle);

    double leftFlowArea = cell.dutoL.area;
    double leftLengthRatio = leftCell.dxL / (leftCell.dx + leftCell.dxL);
    double leftMeanPressure = leftCell.presaux;
    double leftMeanTemperature = leftCell.temp * leftLengthRatio + leftCell.tempL * (1. - leftLengthRatio);
    double betIL = 0.;
    if (cellIndex < 2)
        betIL = leftCell.betL;
    else
        betIL = state.cells[cellIndex - 2].betPigD;
    if (leftCell.MliqiniBuf < 0.)
        betIL = leftCell.betPigE; // beta test
    // betIL = leftCell.betPigE;    // beta doubt
    double leftGasDensity = cell.flui.MasEspGas(leftMeanPressure, leftMeanTemperature);
    double leftLiquidDensity = (1 - betIL) * cell.flui.MasEspLiq(leftMeanPressure, leftMeanTemperature) + betIL * cell.fluicol.MasEspFlu(leftMeanPressure, leftMeanTemperature);
    double leftSuperficialGasVelocity = (cell.MLBuf - cell.MliqiniLBuf) / (leftGasDensity * leftFlowArea);
    double leftSuperficialLiquidVelocity = (cell.MliqiniLBuf) / (leftLiquidDensity * leftFlowArea);

    double rightFlowArea = cell.dutoR.area;
    double rightLengthRatio = cell.dxR / (cell.dxR + cell.dx);
    double rightMeanPressure = cell.presRBuf;
    double rightMeanTemperature = cell.temp * rightLengthRatio + cell.tempR * (1. - rightLengthRatio);
    double betIR = cell.betPigD;
    if (cell.MliqiniRBuf < 0.) { // beta test
        if (cellIndex > state.lastCell - 2)
            betIR = cell.betR;
        else
            betIR = state.cells[cellIndex + 1].betPigE;
    }
    double rightGasDensity = cell.flui.MasEspGas(rightMeanPressure, rightMeanTemperature);
    double rightLiquidDensity = (1 - betIR) * cell.flui.MasEspLiq(rightMeanPressure, rightMeanTemperature) + betIR * cell.fluicol.MasEspFlu(rightMeanPressure, rightMeanTemperature);
    double rightSuperficialGasVelocity = (cell.MRBuf - cell.MliqiniRBuf) / (rightGasDensity * rightFlowArea);
    double rightSuperficialLiquidVelocity = (cell.MliqiniRBuf) / (rightLiquidDensity * rightFlowArea);

    selectAndApplyBufferedOutletFlowRegime(
        state, cellIndex, c0, ud, superficialGasVelocity, superficialLiquidVelocity, leftSuperficialGasVelocity, leftSuperficialLiquidVelocity, rightSuperficialGasVelocity, rightSuperficialLiquidVelocity,
        rightGasDensity, rightLiquidDensity, gasDensity, liquidDensity, flowArea);
    cell.term1L = leftCell.term1;
    cell.term2L = leftCell.term2;
    leftCell.term1R = cell.term1;
    leftCell.term2R = cell.term2;
}
void updateInletFlowPartitionTerms(const ThermalState &state) {


    int cellIndex = 0;
    Cel &cell = state.cells[cellIndex];

    double meanPressure;
    meanPressure = state.inletPressure;

    double meanTemperature;
    if (state.cells[0].MliqiniBuf < 0.)
        meanTemperature = cell.temp;
    else
        meanTemperature = state.inletTemperature;

    double betI;
    double surfaceTension;
    double gasDensity = state.cells[0].flui.MasEspGas(state.inletPressure, state.inletTemperature);
    double liquidDensity = state.cells[0].flui.MasEspLiq(state.inletPressure, state.inletTemperature);
    double rcis = state.cells[0].fluicol.MasEspFlu(state.inletPressure, state.inletTemperature);

    double liquidMixtureDensity = state.inletCompletionFraction * rcis + (1 - state.inletCompletionFraction) * liquidDensity;
    state.inletVoidFraction = (-state.inletQuality * liquidMixtureDensity / (state.inletQuality * gasDensity - gasDensity - state.inletQuality * liquidMixtureDensity)) / (state.cells[0].c0);

    if ((cell.MCBuf - state.cells[0].MliqiniBuf) * 0 + 1 * state.cells[0].MliqiniBuf < 0.) { // beta doubt
        betI = cell.betPigE;
        liquidDensity = (1 - betI) * cell.flui.MasEspLiq(meanPressure, meanTemperature) + betI * cell.fluicol.MasEspFlu(meanPressure, meanTemperature);
        surfaceTension = (1 - betI) * cell.flui.TensSuper(meanPressure, meanTemperature) + betI * cell.fluicol.TensSuper(meanPressure, meanTemperature);
    } else {
        betI = state.inletCompletionFraction;
        liquidDensity = (1 - betI) * (*cell.fluiL).MasEspLiq(meanPressure, meanTemperature) + betI * cell.fluicol.MasEspFlu(meanPressure, meanTemperature);
        surfaceTension = (1 - betI) * (*cell.fluiL).TensSuper(meanPressure, meanTemperature) + betI * cell.fluicol.TensSuper(meanPressure, meanTemperature);
    }
    double flowArea;
    if (cell.MCBuf - state.cells[0].MliqiniBuf >= 0) {
        flowArea = cell.dutoL.area;
        gasDensity = (*cell.fluiL).MasEspGas(meanPressure, meanTemperature);
    } else {
        gasDensity = cell.flui.MasEspGas(meanPressure, meanTemperature);
        flowArea = cell.duto.area;
    }
    double superficialGasVelocity = (cell.MCBuf - cell.MliqiniBuf) / (gasDensity * flowArea);
    double superficialLiquidVelocity = cell.MliqiniBuf / (liquidDensity * flowArea);

    double inclinationAngle = cell.duto.teta;
    double inclinationSign = 1.;
    if (inclinationAngle < 0.)
        inclinationSign = -1.;
    double c0 = 2.;
    double ud = inclinationSign * 0.0246 * cos(inclinationAngle) + 1.606 * pow(kGravity * surfaceTension * (liquidDensity - gasDensity) / (liquidDensity * liquidDensity), 0.25) * sin(inclinationAngle);

    double rightFlowArea = cell.dutoR.area;
    double rightLengthRatio = cell.dxR / (cell.dxR + cell.dx);
    double rightMeanPressure = cell.presRBuf * rightLengthRatio + cell.presBuf * (1. - rightLengthRatio);
    double rightMeanTemperature = cell.temp * rightLengthRatio + cell.tempR * (1. - rightLengthRatio);
    double betIR = cell.betPigD;
    if (cell.QLR < 0.) // beta test
        betIR = state.cells[cellIndex + 1].betPigE;

    double rightGasDensity = cell.flui.MasEspGas(rightMeanPressure, rightMeanTemperature);
    double rightLiquidDensity = (1 - betIR) * cell.flui.MasEspLiq(rightMeanPressure, rightMeanTemperature) + betIR * cell.fluicol.MasEspFlu(rightMeanPressure, rightMeanTemperature);
    double rightSuperficialLiquidVelocity = (cell.MliqiniRBuf) / (rightLiquidDensity * rightFlowArea);

    selectAndApplyBufferedInletFlowRegime(
        state, cellIndex, c0, ud, superficialGasVelocity, superficialLiquidVelocity, rightSuperficialLiquidVelocity, rightGasDensity, rightLiquidDensity, gasDensity, liquidDensity,
        flowArea);
    state.cells[1].term1L = cell.term1;
    state.cells[1].term2L = cell.term2;
}
void prepareNonDimensionalHeatDiffusion(const ThermalState &state, int cellIndex) {
    Cel &cell = state.cells[cellIndex];
    Cel &leftCell = state.cells[cellIndex - 1];
    const CellFlowBasis basis = cellFlowBasisOf(state, cellIndex);
    const double meanVoidFraction = basis.voidFraction;
    const double betmed = basis.betmed;
    const double meanSuperficialGasVelocity = basis.gasSuperficialVelocity;
    const double meanSuperficialLiquidVelocity = basis.liquidSuperficialVelocity;
    double rp = cell.rpC;
    double rc = cell.rcC;
    double liquidDensity = (1. - betmed) * rp + betmed * rc;
    double gasDensity = cell.rgC;
    double liquidSpecificHeat = (1. - betmed) * cell.flui.CalorLiq(cell.presini, cell.temp) + betmed * cell.fluicol.CalorLiq(cell.presini, cell.temp);
    double gasSpecificHeat = cell.flui.CalorGas(cell.presini, cell.temp);

    cell.calor.Tint = cell.temp;
    cell.calor.dtL = cell.temp - leftCell.tempini;
    cell.calor.Vint = meanSuperficialGasVelocity + meanSuperficialLiquidVelocity;
    cell.calor.dt = cell.dt;
    double liquidConductivity = (1. - betmed) * cell.flui.CondLiq(cell.presini, cell.temp) + betmed * cell.fluicol.CondLiq(cell.presini, cell.temp); //(1. - betmed) * celula[cellIndex].flui.CondLiq(celula[cellIndex].pres, celula[cellIndex].temp)
    cell.calor.kint = liquidConductivity * (1 - meanVoidFraction) + cell.flui.CondGas(cell.presini, cell.temp) * meanVoidFraction;                                                 // liquidConductivity * (1 - meanVoidFraction) + celula[cellIndex].flui.CondGas(celula[cellIndex].pres, celula[cellIndex].temp) * meanVoidFraction;
    cell.calor.cpint = liquidSpecificHeat * (1 - meanVoidFraction) + gasSpecificHeat * meanVoidFraction;
    cell.calor.rhoint = liquidDensity * (1 - meanVoidFraction) + gasDensity * meanVoidFraction;
    //(1. - betmed) * celula[cellIndex].flui.ViscOleo(celula[cellIndex].pres, celula[cellIndex].temp)
    double liquidViscosity = (1. - betmed) * cell.mipC + betmed * cell.micC;
    // liquidViscosity * (1 - meanVoidFraction) * 1.e-3
    cell.calor.viscint = liquidViscosity * (1 - meanVoidFraction) * kPascalSecondPerCentipoise + cell.migC * meanVoidFraction * kPascalSecondPerCentipoise;
    double temperaturePerturbation = cell.temp * 0.01;
    if (fabs(cell.temp) < 1e-15)
        temperaturePerturbation = 0.1;
    double liquidDensityChange = (1. - betmed) * cell.flui.MasEspLiq(cell.presini, cell.temp + temperaturePerturbation) +
                    betmed * cell.fluicol.MasEspFlu(cell.presini, cell.temp + temperaturePerturbation) - liquidDensity; //(1. - betmed) * celula[cellIndex].flui.MasEspLiq(celula[cellIndex].pres, celula[cellIndex].temp+temperaturePerturbation) +
    double gasDensityChange = cell.flui.MasEspGas(cell.presini, cell.temp + temperaturePerturbation) - gasDensity;             // celula[cellIndex].flui.MasEspGas(celula[cellIndex].pres, celula[cellIndex].temp+temperaturePerturbation)-gasDensity;
    cell.calor.betint = -(1 / cell.calor.rhoint) * (liquidDensityChange * (1 - meanVoidFraction) + gasDensityChange * meanVoidFraction) / (temperaturePerturbation);
}

namespace {

/// Refreshes one cell's temperature and its rate of change. Beyond lastCell the
/// temperature comes from the external boundary or the gas surface.
void refreshCellTemperatureAndRate(const ThermalState &state, int cellIndex) {
    Cel &cell = state.cells[cellIndex];
    if (cellIndex <= state.lastCell) {
        computeTemperature(state, cellIndex, cell.tempini);
    } else {
        if ((*state.globals).chaverede == 0 || state.endNode == 1 || (*state.globals).chaveRedeParalela == 1)
            cell.temp = cell.calor.Textern1;
        else
            cell.temp = state.gasSurfaceTemperature;
    }
    cell.dTdt = (cell.temp - cell.tempini) / state.timeStep;
    cell.dTdtIni = cell.dTdt;
}

/// Refreshes the temperature field over the 2D-Poisson cells and then over the
/// column. Both loops run on the global thread count.
void refreshTemperatureField(const ThermalState &state) {
if (state.poisson2DCellCount > 1) {
#pragma omp parallel for num_threads((*state.globals).ntrd)
    for (int poisson2DIndex = 0; poisson2DIndex < state.poisson2DCellCount; poisson2DIndex++) {
        int cellIndex = state.poisson2DCellIndices[poisson2DIndex];
        refreshCellTemperatureAndRate(state, cellIndex);
    }
}
#pragma omp parallel for num_threads((*state.globals).ntrd)
for (int cellIndex = 0; cellIndex <= state.lastCell; cellIndex++) {
    if (state.cells[cellIndex].calor.difus2D == 0) {
        refreshCellTemperatureAndRate(state, cellIndex);
    }
}
}

}  // namespace

void advanceTransientEnergy(const ThermalState &state, int cycle, int maximumCycle) {
    if (((*state.globals).chaverede == 0 || state.endNode == 1 || (*state.globals).chaveRedeParalela == 1) && state.input.chkv == 0) {
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
        if (state.cells[state.lastCell - 1].QLR < 0 && (state.surfaceChokeMassFlag == 0 || state.surfaceChokeOpen == 1)) {
            state.cells[state.lastCell - 1].QLR = 0;
            state.cells[state.lastCell].QL = 0;
        }
    } else if ((state.input.chkv == 1 && state.surfaceChokeMassFlag == 0) || (state.input.chkv == 1 && state.surfaceChokeMassFlag == 1)) {
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
        if (state.cells[state.lastCell - 1].QLR < 0 && (state.surfaceChokeMassFlag == 0 || state.surfaceChokeOpen == 1)) {
            state.cells[state.lastCell - 1].QLR = 0;
            state.cells[state.lastCell].QL = 0;
        }
    }

    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // attention!!!!!!!!!!!!!!!!!!
    // there is an issue that looks unresolved in this march: nothing was tested for
    // the case where the temperature transport velocity is <0; then the temperature in the cell with index
    // should not enter calctemp, since it is no longer a cell with a boundary
    // condition for temperature??????????????????????????????????????????????????????????????????????/
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
        refreshTemperatureField(state);
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
                state.poissonSolver3D.dados.tInt[nextCouplingIndex] = state.cells[coupledCellIndex].temp;
                double cellInterfaceCoefficient = state.cells[coupledCellIndex].calor.hInt();
                state.poissonSolver3D.dados.hI[nextCouplingIndex] = cellInterfaceCoefficient;
            }
        }
        if (cycle < maximumCycle || maximumCycle == 0 || state.minimumCycleTimeStep != state.timeStep) {
            if (cycle == maximumCycle && maximumCycle > 0)
                state.poissonSolver3D.FeiticoDoTempo();
            state.poissonSolver3D.transientePoisson(state.timeStep);
        }
        refreshTemperatureField(state);
    }
    for (int cellIndex = 1; cellIndex <= state.lastCell; cellIndex++) {
        Cel &cell = state.cells[cellIndex];
        Cel &leftCell = state.cells[cellIndex - 1];
        Cel &rightCell = state.cells[cellIndex + 1];
        cell.dTdtL = leftCell.dTdt;
        if (cellIndex < state.lastCell) {
            rightCell.tempLini = rightCell.tempL;
            rightCell.tempL = cell.temp;
        }
        leftCell.tempRini = leftCell.tempR;
        leftCell.tempR = cell.temp;
    }
    if (cycle < maximumCycle) {
        for (int k = 0; k <= state.lastCell; k++) {
            state.cells[k].FeiticoDoTempo();
        }
    } else if (state.input.modoDifus3D == 1)
        state.poissonSolver3D.renova();
    if (state.fullModel == 0) {
        if (cycle < maximumCycle) {
            state.evolutionUpdater.solvePressureVelocityCoupling(cycle);
            state.evolutionUpdater.renew();
        }
    }
}

/// Energy added to the cell by its mass sources, for the steady march.
namespace {
// Helpers of the steady march; internal linkage keeps them unexported.

TemperatureSourceTerms computeSteadySourceTerms(const ThermalState &state,
                                                int cellIndex,
                                                double meanPressure,
                                                double meanTemperature,
                                                double cellLength) {
    Cel &cell = state.cells[cellIndex];
    Cel &leftCell = state.cells[cellIndex - 1];
    double gasMassSourceTerm = 0.;
    double liquidMassSourceTerm = 0.;
    [[maybe_unused]] double completionMassSourceTerm = 0.;
    double sourceTemperature = leftCell.temp;
    double sourceGasSpecificHeat;
    double sourceSpecificHeatRatio = 0.;
    double sourceLiquidSpecificHeat;

    // energy added to the system by mass sources
    if (leftCell.acsr.tipo == kAccessoryGasInjection) { // gas source case
        sourceTemperature = leftCell.acsr.injg.temp;
        sourceGasSpecificHeat = leftCell.acsr.injg.FluidoPro.CalorGas(meanPressure, sourceTemperature);
        sourceSpecificHeatRatio = leftCell.acsr.injg.FluidoPro.ConstAdG(meanPressure, sourceTemperature);
        sourceLiquidSpecificHeat = 0.;
    } else if (leftCell.acsr.tipo == kAccessoryLiquidInjection) { // liquid source case
        sourceTemperature = leftCell.acsr.injl.temp;
        sourceGasSpecificHeat = 0.;
        sourceSpecificHeatRatio = 1.;
        sourceLiquidSpecificHeat = (1. - leftCell.acsr.injl.bet) * leftCell.acsr.injl.FluidoPro.CalorLiq(meanPressure, meanTemperature) + leftCell.acsr.injl.bet * leftCell.acsr.injl.fluidocol.CalorLiq(meanPressure, meanTemperature);
    } else if (leftCell.acsr.tipo == kAccessoryInflowPerformance) { // IPR case
        sourceTemperature = leftCell.acsr.ipr.Tres;
        sourceGasSpecificHeat = leftCell.acsr.ipr.FluidoPro.CalorGas(meanPressure, meanTemperature);
        sourceSpecificHeatRatio = leftCell.acsr.ipr.FluidoPro.ConstAdG(meanPressure, meanTemperature);
        sourceLiquidSpecificHeat = leftCell.acsr.ipr.FluidoPro.CalorLiq(meanPressure, meanTemperature);
    } else if (leftCell.acsr.tipo == kAccessoryLeak && leftCell.acsr.fontechk.abertura > 1e-6 &&
               (leftCell.fontemassCR + leftCell.fontemassGR + leftCell.fontemassLR) > 1e-9) {
        // leak case
        sourceTemperature = leftCell.acsr.fontechk.tamb;
        sourceGasSpecificHeat = leftCell.acsr.fontechk.fluidoPamb.CalorGas(meanPressure, meanTemperature);
        sourceSpecificHeatRatio = leftCell.acsr.fontechk.fluidoPamb.ConstAdG(meanPressure, meanTemperature);
        sourceLiquidSpecificHeat = (1. - leftCell.acsr.fontechk.betISamb) *
                   leftCell.acsr.fontechk.fluidoPamb.CalorLiq(meanPressure, meanTemperature) +
               leftCell.acsr.fontechk.betISamb * leftCell.acsr.fontechk.fluidocol.CalorLiq(meanPressure, meanTemperature);
    } else if (leftCell.acsr.tipo == kAccessoryRadialPorous) { // IPR case
        sourceTemperature = leftCell.acsr.radialPoro.tRes;
        sourceGasSpecificHeat = leftCell.acsr.radialPoro.flup.CalorGas(meanPressure, meanTemperature);
        sourceSpecificHeatRatio = leftCell.acsr.radialPoro.flup.ConstAdG(meanPressure, meanTemperature);
        sourceLiquidSpecificHeat = leftCell.acsr.radialPoro.flup.CalorLiq(meanPressure, meanTemperature);
    } else if (leftCell.acsr.tipo == kAccessoryPorous2D) { // IPR case
        sourceTemperature = leftCell.acsr.poroso2D.dados.tRes;
        sourceGasSpecificHeat = leftCell.acsr.poroso2D.dados.flup.CalorGas(meanPressure, meanTemperature);
        sourceSpecificHeatRatio = leftCell.acsr.poroso2D.dados.flup.ConstAdG(meanPressure, meanTemperature);
        sourceLiquidSpecificHeat = leftCell.acsr.poroso2D.dados.flup.CalorLiq(meanPressure, meanTemperature);
    } else if (cell.acsrL != 0) {

        clearSourceSpecificHeats(sourceGasSpecificHeat, sourceSpecificHeatRatio,
                                 sourceLiquidSpecificHeat);
    } else {
        clearSourceSpecificHeats(sourceGasSpecificHeat, sourceSpecificHeatRatio,
                                 sourceLiquidSpecificHeat);
    }

    liquidMassSourceTerm = 0;
    if (leftCell.fontemassLR > 0.)
        liquidMassSourceTerm = leftCell.fontemassLR / cellLength;
    if (leftCell.fontemassCR > 0.)
        liquidMassSourceTerm += leftCell.fontemassCR / cellLength;
    liquidMassSourceTerm *= (sourceLiquidSpecificHeat) * (sourceTemperature - leftCell.temp);

    gasMassSourceTerm = leftCell.fontemassGR / cellLength;
    if (gasMassSourceTerm > 0.)
        gasMassSourceTerm *= (sourceGasSpecificHeat / sourceSpecificHeatRatio) * (sourceTemperature - leftCell.temp);
    else
        gasMassSourceTerm = 0;
    return TemperatureSourceTerms{.gas = gasMassSourceTerm,
                                  .liquid = liquidMassSourceTerm};
}

/// Annulus coupling for the steady march: fills the cell's external heat
/// transfer state from the gas-line cell facing it, and returns the wall
/// resistance.
double applySteadyAnnulusCoupling(const ThermalState &state, int cellIndex,
                                  double interfaceMeanPressure,
                                  double interfaceMeanTemperature,
                                  double gasSpecificHeat,
                                  double gasDensity) {
    Cel &leftCell = state.cells[cellIndex - 1];
    double annulusResistance = 0.;
    // checks whether there is coupling with the annulus:
    if (state.input.lingas == 1 && (cellIndex - 1 <= state.annulusTubingStart && cellIndex - 1 >= state.annulusTubingEnd)) {
        int j = state.annulusTubingStart + state.tubingAnnulusStart - (cellIndex - 1);
        // with coupling:

        if (state.steadyIteration == 0 && (cellIndex - 1) < state.annulusTubingStart) {
            // on the first iteration, the resistance of casing + cement +
            // formation is considered; on the other iterations
            // it is enough to determine the heat exchange between the tubing and the annulus gas
            // in the coupled model, the tubing has no casing+cement+formation wall definition;
            // to obtain this resistance, the heat-exchange model of the
            // annulus has to be used, which is done here:
            // Note: done up to one cell before reaching the master, (cellIndex-1)<ColunaAnulaIni,
            // the same as the involuta simulator does, to improve
            // the temperature estimate at the ANM. In the ANM cell only the heat exchange with the gas is considered,
            // from the first iteration on
            state.gasCells[j].calor.Tint = interfaceMeanTemperature;
            state.gasCells[j].calor.Vint = 100;
            state.gasCells[j].calor.kint = leftCell.flui.CondGas(interfaceMeanPressure, interfaceMeanTemperature);
            state.gasCells[j].calor.cpint = gasSpecificHeat;
            state.gasCells[j].calor.rhoint = gasDensity;
            state.gasCells[j].calor.viscint = leftCell.flui.ViscGas(interfaceMeanPressure, interfaceMeanTemperature) * kPascalSecondPerCentipoise;
            state.gasCells[j].fluxcal = state.gasCells[j].calor.transperm(); // heat exchange in the annulus
            annulusResistance = state.gasCells[j].calor.resGlob;                // wall resistance
            // The heat flux itself is not of interest here, only the resistance
            // of the walls, from the casing towards the formation.
            leftCell.calor.Vextern1 = 100.;
            leftCell.calor.kextern1 = state.gasCells[j].calor.kint;
            leftCell.calor.cpextern1 = state.gasCells[j].calor.cpint;
            leftCell.calor.rhoextern1 = state.gasCells[j].calor.rhoint;
            leftCell.calor.viscextern1 = state.gasCells[j].calor.viscint;
        } else { // after the first iteration, only the heat exchange between the tubing and the annulus gas is considered,
            // through the tubing wall, of course
            annulusResistance = 0.;
            leftCell.calor.Vextern1 = state.gasCells[j].VGasR / state.gasCells[j].u1L;
            leftCell.calor.kextern1 = state.gasCells[j].calor.kint;
            leftCell.calor.cpextern1 = state.gasCells[j].calor.cpint;
            leftCell.calor.rhoextern1 = state.gasCells[j].calor.rhoint;
            leftCell.calor.viscextern1 = state.gasCells[j].calor.viscint;
        }
        if (state.steadyIteration == 0)
            leftCell.calor.Textern1 = state.gasCells[j].calor.Textern1; // on the first iteration, as the calculation
        // uses the whole well resistance, so the external temperature is the geothermal one
        else
            leftCell.calor.Textern1 = state.gasCells[j].temp; // on the following iterations, the ambient temperature
        // is the gas temperature
    }
    return annulusResistance;
}

/// Kinetic-energy term of the steady march.
double computeSteadyKineticTerm(const ThermalState &state, int cellIndex,
                                double meanSuperficialGasVelocity,
                                double meanSuperficialLiquidVelocity) {
    Cel &cell = state.cells[cellIndex];
    Cel &leftCell = state.cells[cellIndex - 1];
    double kineticTerm = 0;
    double upstreamMeanGasVelocity = 0;
    double upstreamMeanLiquidVelocity = 0;
    double meanGasVelocity = 0;
    double meanLiquidVelocity = 0;
    // kinetic energy term:
    if (cell.acsr.tipo == kAccessoryNone && leftCell.acsr.tipo == kAccessoryNone && cellIndex > 2) {
        double kineticCellLength = leftCell.dx;
        double upstreamDiameter = leftCell.duto.a;
        double upstreamFlowArea = 0.25 * M_PI * upstreamDiameter * upstreamDiameter;

        if (leftCell.alf > 1e-3)
            meanGasVelocity = meanSuperficialGasVelocity / leftCell.alf;
        if (leftCell.alf < (1. - 1e-3))
            meanLiquidVelocity = meanSuperficialLiquidVelocity / (1. - leftCell.alf);

        if (state.cells[cellIndex - 2].alf > 1e-3) {
            upstreamMeanGasVelocity = leftCell.QG / (upstreamFlowArea);
            upstreamMeanGasVelocity /= state.cells[cellIndex - 2].alf;
        }
        if (state.cells[cellIndex - 2].alf < (1. - 1e-3)) {
            upstreamMeanLiquidVelocity = leftCell.QL / (upstreamFlowArea);
            upstreamMeanLiquidVelocity /= (1. - state.cells[cellIndex - 2].alf);
        }

        if (state.input.nCompTotalUnidadesP / kineticCellLength < 1e6)
            kineticTerm = (cell.MC - cell.Mliqini) * meanGasVelocity * (meanGasVelocity - upstreamMeanGasVelocity) / kineticCellLength + cell.Mliqini * meanLiquidVelocity * (meanLiquidVelocity - upstreamMeanLiquidVelocity) / kineticCellLength;
        else
            kineticTerm = 0;
    }
    return kineticTerm;
}

/// Latent-heat term of the steady march. interfacialWorkTerm is computed here
/// and never read.
double computeSteadyLatentHeatTerm(const ThermalState &state, int cellIndex,
                                   double flowArea, double meanPressure,
                                   double meanTemperature,
                                   double meanSuperficialGasVelocity,
                                   double meanSuperficialLiquidVelocity) {
    varGlob1D &globals = *state.globals;
    Cel &cell = state.cells[cellIndex];
    Cel &leftCell = state.cells[cellIndex - 1];
    double latentHeatTerm = 0.;
    if (isnan(leftCell.FonteMudaFase))
        leftCell.FonteMudaFase = 0.;
    double phaseChangeMassRate = fabs(leftCell.FonteMudaFase);
    double phaseChangeSign = 1.;
    if (phaseChangeMassRate > kPhaseChangeFloor)
        phaseChangeSign = leftCell.FonteMudaFase / phaseChangeMassRate;
    if (state.input.limTransMass < phaseChangeMassRate)
        phaseChangeMassRate = phaseChangeSign * state.input.limTransMass;
    else
        phaseChangeMassRate *= phaseChangeSign;
    if (cell.flui.dVaporMassFraction < (1 - 1e-15) && cell.flui.dVaporMassFraction > (1e-15)) {
        if (leftCell.acsr.tipo == kAccessoryGasInjection || leftCell.acsr.tipo == kAccessoryLiquidInjection || leftCell.acsr.tipo == kAccessoryInflowPerformance || leftCell.acsr.tipo == kAccessoryRadialPorous || leftCell.acsr.tipo == kAccessoryPorous2D) {
            leftCell.FonteMudaFase = 0.;
            phaseChangeMassRate = 0.;
        }
        if (state.latentHeatEnabled > 0 && state.steadyIteration != 0 && state.input.flashCompleto == 0) {
            latentHeatTerm = interpolateLatentHeat(state, meanPressure, meanTemperature) * phaseChangeMassRate;
        } else if ((state.input.flashCompleto == 1 || state.input.flashCompleto == 2) && state.latentHeatEnabled > 0 && state.steadyIteration != 0) {
            latentHeatTerm = (cell.flui.EntalpGas(meanPressure, meanTemperature) -
                       cell.flui.EntalpLiq(meanPressure, meanTemperature)) *
                      phaseChangeMassRate;
        } else
            latentHeatTerm = 0;
    }

    double interfaceVoidFraction;
    double leftInterfaceVoidFraction;
    if (meanSuperficialGasVelocity > 0) {
        interfaceVoidFraction = leftCell.alf;
        leftInterfaceVoidFraction = leftCell.alfL;
    } else {
        interfaceVoidFraction = leftCell.alfR;
        leftInterfaceVoidFraction = leftCell.alf;
    }
    double slipVelocity;
    if (interfaceVoidFraction > globals.localtiny && interfaceVoidFraction < (1. - globals.localtiny))
        slipVelocity = meanSuperficialGasVelocity / interfaceVoidFraction - meanSuperficialLiquidVelocity / (1. - interfaceVoidFraction);
    else if (interfaceVoidFraction > globals.localtiny)
        slipVelocity = meanSuperficialGasVelocity;
    else
        slipVelocity = meanSuperficialLiquidVelocity;
    [[maybe_unused]] double interfacialWorkTerm = flowArea * cell.pres * kPascalPerKgfPerCm2Coarse * slipVelocity * (interfaceVoidFraction - leftInterfaceVoidFraction) / cell.dx;
    return latentHeatTerm;
}

}  // namespace

void advanceSteadyTemperature(const ThermalState &state, int cellIndex, int rungeKuttaStage) {
    Cel &cell = state.cells[cellIndex];
    varGlob1D &globals = *state.globals;
    Cel &leftCell = state.cells[cellIndex - 1];
    double cellLength = 0.5 * (cell.dx + leftCell.dx);
    double meanCellLength = 0.5 * (cell.dx + leftCell.dx);
    double meanTemperatureGradientCorrection = (cell.dx * cell.dTdLCor + leftCell.dx * leftCell.dTdLCor) / (cell.dx + leftCell.dx);
    double diameter = leftCell.duto.a;
    double flowArea = 0.25 * M_PI * diameter * diameter;
    double meanVoidFraction;
    double betmed;
    double meanPressure;
    double meanTemperature;
    if (rungeKuttaStage == 0) {
        meanVoidFraction = leftCell.alf;
        betmed = leftCell.bet;
        meanPressure = leftCell.pres;
        meanTemperature = leftCell.temp;
    } else {
        meanVoidFraction = cell.alf;
        betmed = cell.bet;
        meanPressure = cell.pres;
        meanTemperature = cell.temp;
    }
    double meanSuperficialGasVelocity;
    meanSuperficialGasVelocity = cell.QG / flowArea; // gas superficial velocity
    double meanSuperficialLiquidVelocity;
    meanSuperficialLiquidVelocity = cell.QL / flowArea; // liquid superficial velocity
    double mixtureFluxSign = 1.;
    if (fabs(meanSuperficialGasVelocity + meanSuperficialLiquidVelocity) > 0.05 && state.thermalSourceDisabled == 0) { // the thermal calculation is done for mixture velocities above 0.1 m/s,
        // at lower velocities the fluid temperature is taken as the ambient temperature
        mixtureFluxSign = (meanSuperficialGasVelocity + meanSuperficialLiquidVelocity) / fabs(meanSuperficialGasVelocity + meanSuperficialLiquidVelocity);
        double interfaceMeanPressure = cell.presaux + 0 * leftCell.dpB / kPascalPerKgfPerCm2;
        double interfaceMeanTemperature;
        if (state.steadyIteration != 0 && state.input.AceleraConvergPerm == 0) // temperature at the cell's left interface
            // case in which one iteration is considered to have been
            // done already, with the temperature of cell cellIndex known from the previous iteration;
            // this can hinder convergence when convergence acceleration is wanted
            // the temperature at the left boundary is taken as that of the left cell
            interfaceMeanTemperature = (cell.dx * cell.temp + cell.dxL * cell.tempL) / (cell.dx + cell.dxL);
        else
            interfaceMeanTemperature = leftCell.temp;
        double rp = leftCell.flui.MasEspLiq(meanPressure, meanTemperature);    // celula[cellIndex].rpCi;
        double rc = leftCell.fluicol.MasEspFlu(meanPressure, meanTemperature); // celula[cellIndex].rcCi;
        double liquidDensity = (1. - betmed) * rp + betmed * rc;
        double gasDensity = leftCell.flui.MasEspGas(meanPressure, meanTemperature); // celula[cellIndex].rgCi;
        double liquidSpecificHeat = (1. - betmed) * leftCell.flui.CalorLiq(interfaceMeanPressure, interfaceMeanTemperature) + betmed * leftCell.fluicol.CalorLiq(interfaceMeanPressure, interfaceMeanTemperature);
        double gasSpecificHeat = leftCell.flui.CalorGas(interfaceMeanPressure, interfaceMeanTemperature);
        // liquidJouleThomson = liquid Joule-Thomson coefficient x cp
        // gasJouleThomson = gas Joule-Thomson coefficient x cp
        /////??????????????????????????????????????????????????????????????????????????????????????????
        double liquidJouleThomson = (1. - betmed) * leftCell.flui.JTL(interfaceMeanPressure, interfaceMeanTemperature) - betmed / rc;
        /////??????????????????????????????????????????????????????????????????????????????????????????
        if (state.input.pocinjec > 0 && state.input.condpocinj.tipoFlui == 2) {
            liquidJouleThomson = -(1 + (interfaceMeanTemperature + 273.14) * leftCell.fluicol.DrhoDtFlu(interfaceMeanPressure, interfaceMeanTemperature) / rc) / rc;
        }
        double gasJouleThomson = leftCell.flui.JTG(interfaceMeanPressure, interfaceMeanTemperature);
        // potential energy:
        double hydrostaticPower = (liquidDensity * meanSuperficialLiquidVelocity + gasDensity * meanSuperficialGasVelocity) * flowArea * kGravity * sin(leftCell.duto.teta);


        // internal tubing parameters for the heat exchange with the surroundings
        leftCell.calor.Tint = interfaceMeanTemperature;
        leftCell.calor.Vint = fabs(meanSuperficialGasVelocity + meanSuperficialLiquidVelocity);
        double liquidConductivity = (1. - betmed) * leftCell.flui.CondLiq(interfaceMeanPressure, interfaceMeanTemperature) + betmed * leftCell.fluicol.CondLiq(interfaceMeanPressure, interfaceMeanTemperature);
        leftCell.calor.kint = liquidConductivity * (1 - meanVoidFraction) + leftCell.flui.CondGas(interfaceMeanPressure, interfaceMeanTemperature) * meanVoidFraction;
        leftCell.calor.cpint = liquidSpecificHeat * (1 - meanVoidFraction) + gasSpecificHeat * meanVoidFraction;
        leftCell.calor.rhoint = liquidDensity * (1 - meanVoidFraction) + gasDensity * meanVoidFraction;
        double liquidViscosity = (1. - betmed) * leftCell.flui.ViscOleo(interfaceMeanPressure, interfaceMeanTemperature) + betmed * leftCell.fluicol.VisFlu(interfaceMeanPressure, interfaceMeanTemperature);
        leftCell.calor.viscint = liquidViscosity * (1 - meanVoidFraction) * kPascalSecondPerCentipoise + leftCell.flui.ViscGas(interfaceMeanPressure, interfaceMeanTemperature) * meanVoidFraction * kPascalSecondPerCentipoise;

        [[maybe_unused]] double relaxedHeatFlux = 0; // unused variable
        if (state.steadyIteration > 0)
            relaxedHeatFlux = leftCell.fluxcalmed;
        double annulusResistance = 0.;
        double heatFlux;
        [[maybe_unused]] double gasHeatFlux;
        annulusResistance = applySteadyAnnulusCoupling(
            state, cellIndex, interfaceMeanPressure, interfaceMeanTemperature,
            gasSpecificHeat, gasDensity);
        leftCell.fluxcalmed = 0;
        if (state.productionNetworkHeatCoupled == 1 && (cellIndex - 1) >= state.primaryNetworkSectionEnd && (cellIndex - 1) <= state.primaryNetworkSectionStart) {
            heatFlux = mixtureFluxSign * leftCell.calor.transperm(leftCell.resAcopRedeP);
        } else
            heatFlux = mixtureFluxSign * leftCell.calor.transperm(annulusResistance);
        leftCell.fluxcalmed = heatFlux; // heat flux in the tubing

        double temperatureSpatialCoefficient = (liquidDensity * meanSuperficialLiquidVelocity * liquidSpecificHeat + gasDensity * meanSuperficialGasVelocity * gasSpecificHeat) * flowArea; // multiplying term
        // the derivative Dt/Dx
        double cappedMeanSuperficialLiquidVelocity = meanSuperficialLiquidVelocity;
        if (globals.blackOilTemp == 1 && fabs(meanSuperficialLiquidVelocity) > 5)
            cappedMeanSuperficialLiquidVelocity = 5 * meanSuperficialLiquidVelocity / fabs(meanSuperficialLiquidVelocity);
        double cappedMeanSuperficialGasVelocity = meanSuperficialGasVelocity;
        if (globals.blackOilTemp == 1 && fabs(meanSuperficialGasVelocity) > 5)
            cappedMeanSuperficialGasVelocity = 5 * meanSuperficialGasVelocity / fabs(meanSuperficialGasVelocity);
        double pressureSpatialCoefficient = 1 * (liquidDensity * cappedMeanSuperficialLiquidVelocity * liquidJouleThomson + gasDensity * cappedMeanSuperficialGasVelocity * gasJouleThomson) * flowArea; // multiplying term
        // the derivative Dp/Dx
        double pressureGradient;
        if ((leftCell.acsr.tipo != kAccessoryPump || leftCell.acsr.bcs.freq < 1) && leftCell.acsr.tipo != 7)
            // without an ESP or a pressure increment, the pressure at the left boundary
            // and the pressure at the cell centre are used to compute Dp/Dx
            pressureGradient = 2. * (cell.presaux - leftCell.pres) * kPascalPerKgfPerCm2 / leftCell.dx;
        else {
            // with an ESP or a pressure increment, the pressure of the cell to the left
            // and the pressure at the cell centre are used to compute Dp/Dx
            pressureGradient = 2. * (cell.presaux - leftCell.pres) * kPascalPerKgfPerCm2 / leftCell.dx;
        }
        cell.VTemper = meanSuperficialLiquidVelocity; // this velocity is only useful in the transient case
        // only so that there is a value when the transient simulation starts
        [[maybe_unused]] double temperatureGradient = (-leftCell.temp) / meanCellLength;

        double kineticTerm = computeSteadyKineticTerm(
            state, cellIndex, meanSuperficialGasVelocity,
            meanSuperficialLiquidVelocity);

        TemperatureSourceTerms steadySources = computeSteadySourceTerms(
            state, cellIndex, meanPressure, meanTemperature, cellLength);
        double gasMassSourceTerm = steadySources.gas;
        double liquidMassSourceTerm = steadySources.liquid;

        // latent heat effect, when requested
        double latentHeatTerm = computeSteadyLatentHeatTerm(
            state, cellIndex, flowArea, meanPressure, meanTemperature,
            meanSuperficialGasVelocity, meanSuperficialLiquidVelocity);

        if (fabs(temperatureSpatialCoefficient) > globals.localtiny) {
            // Energy terms for boundary work, potential and kinetic energy,
            // mass sources, latent heat, and shaft work.
            if (state.input.latente == 0)
                latentHeatTerm = 0.;
            else if (state.input.condlatente == 0 && latentHeatTerm < 0)
                latentHeatTerm = 0.;
            double sourceTemperatureGradient = meanTemperatureGradientCorrection * (pressureSpatialCoefficient * pressureGradient - kineticTerm - (hydrostaticPower) + (liquidMassSourceTerm + gasMassSourceTerm) - latentHeatTerm + (leftCell.potTermo + leftCell.fonteCal) / meanCellLength) / temperatureSpatialCoefficient;
            // portion of energy related to heat exchange
            double heatFluxTemperatureGradient = meanTemperatureGradientCorrection * (heatFlux) / temperatureSpatialCoefficient;

            // Check whether heat transfer is too fast for the explicit temperature update.
            // If so, use additional substeps to prevent thermal instability.
            int subStepCount;
            double subStepLength;
            double stabilityLength = fabs(temperatureSpatialCoefficient) / meanTemperatureGradientCorrection;
            if (meanCellLength / (leftCell.calor.resGlob + annulusResistance) < (stabilityLength + 0. * 1000.)) { // Thermal resistance is not low.
                subStepCount = 1;
                subStepLength = meanCellLength;
            } else {
                // For low thermal resistance, determine the number of temperature update steps.
                subStepCount = (meanCellLength / (leftCell.calor.resGlob + annulusResistance)) / (stabilityLength + 0 * 1000) + 1;
                subStepLength = meanCellLength / subStepCount; // cell divided into npassos
            }
            double subStepTemperature = leftCell.temp;

            subStepTemperature = subStepLength * (-(-leftCell.temp) / subStepLength + sourceTemperatureGradient + heatFluxTemperatureGradient); // first step
            for (int j = 1; j < subStepCount; j++) {                                                  // next steps
                // Keep all energy terms except heat flow, which is recalculated at each step.
                leftCell.calor.Tint = subStepTemperature;
                leftCell.calor.Vint = meanSuperficialGasVelocity + meanSuperficialLiquidVelocity;
                liquidConductivity = (1. - betmed) * leftCell.flui.CondLiq(meanPressure, subStepTemperature) + betmed * leftCell.fluicol.CondLiq(meanPressure, subStepTemperature);
                leftCell.calor.kint = liquidConductivity * (1 - meanVoidFraction) + leftCell.flui.CondGas(meanPressure, subStepTemperature) * meanVoidFraction;
                leftCell.calor.cpint = liquidSpecificHeat * (1 - meanVoidFraction) + gasSpecificHeat * meanVoidFraction;
                leftCell.calor.rhoint = liquidDensity * (1 - meanVoidFraction) + gasDensity * meanVoidFraction;
                liquidViscosity = (1. - betmed) * leftCell.flui.ViscOleo(meanPressure, subStepTemperature) + betmed * leftCell.fluicol.VisFlu(meanPressure, subStepTemperature);
                leftCell.calor.viscint = liquidViscosity * (1 - meanVoidFraction) * kPascalSecondPerCentipoise + leftCell.flui.ViscGas(meanPressure, subStepTemperature) * meanVoidFraction * kPascalSecondPerCentipoise;
                if (state.steadyIteration != 0 && state.input.lingas == 1 && (cellIndex - 1 <= state.annulusTubingStart && cellIndex - 1 >= state.annulusTubingEnd)) {
                    int k = state.annulusTubingStart + state.tubingAnnulusStart - (cellIndex - 1);
                    double externalTemperatureStep = (state.gasCells[k - 1].temp - state.gasCells[k].temp) / subStepCount;
                    leftCell.calor.Textern1 = state.gasCells[k].temp + (j - 1) * externalTemperatureStep;
                }
                heatFlux = mixtureFluxSign * leftCell.calor.transperm(annulusResistance);
                leftCell.fluxcalmed += heatFlux;
                heatFluxTemperatureGradient = meanTemperatureGradientCorrection * (heatFlux) / temperatureSpatialCoefficient;
                subStepTemperature = subStepLength * (-(-subStepTemperature) / subStepLength + sourceTemperatureGradient + heatFluxTemperatureGradient);
            }

            cell.temp = subStepTemperature;
            leftCell.fluxcalmed /= subStepCount;
            leftCell.fluxcalmed = 1. * leftCell.fluxcalmed;

            if (cell.temp < kMinimumTemperatureCelsius)
                cell.temp = kMinimumTemperatureCelsius;
            if (cell.temp > kMaximumTemperatureCelsius)
                cell.temp = kMaximumTemperatureCelsius;
        } else
            cell.temp = cell.calor.Textern1;
    } else { // case where the mixing speed is too low
        cell.temp = leftCell.calor.Textern1;
        if (state.input.lingas == 1 && (cellIndex - 1 <= state.annulusTubingStart && cellIndex - 1 >= state.annulusTubingEnd)) {
            int j = state.annulusTubingStart + state.tubingAnnulusStart - (cellIndex - 1);
            leftCell.calor.Vextern1 = 100.;
            leftCell.calor.kextern1 = state.gasCells[j].calor.kint;
            leftCell.calor.cpextern1 = state.gasCells[j].calor.cpint;
            leftCell.calor.rhoextern1 = state.gasCells[j].calor.rhoint;
            leftCell.calor.viscextern1 = state.gasCells[j].calor.viscint;
            [[maybe_unused]] double heatFlux = leftCell.calor.transperm(0);
        }
    }
}

/// Energy added to the cell by its mass sources, for the REVERSE march.
namespace {
// Helpers of the steady march; internal linkage keeps them unexported.

TemperatureSourceTerms computeReverseSteadySourceTerms(
    const ThermalState &state, int cellIndex, double meanPressure,
    double meanTemperature, double cellLength) {
    Cel &rightCell = state.cells[cellIndex + 1];
    double gasMassSourceTerm = 0.;
    double liquidMassSourceTerm = 0.;
    double sourceTemperature = rightCell.temp;
    double sourceGasSpecificHeat;
    double sourceSpecificHeatRatio = 0.;
    double sourceLiquidSpecificHeat;

    // energy added to the system by mass sources
    if (rightCell.acsr.tipo == kAccessoryGasInjection) { // gas source case
        sourceTemperature = rightCell.acsr.injg.temp;
        sourceGasSpecificHeat = rightCell.acsr.injg.FluidoPro.CalorGas(meanPressure, sourceTemperature);
        sourceSpecificHeatRatio = rightCell.acsr.injg.FluidoPro.ConstAdG(meanPressure, sourceTemperature);
        sourceLiquidSpecificHeat = 0.;
    } else if (rightCell.acsr.tipo == kAccessoryLiquidInjection) { // liquid source case
        sourceTemperature = rightCell.acsr.injl.temp;
        sourceGasSpecificHeat = 0.;
        sourceSpecificHeatRatio = 1.;
        sourceLiquidSpecificHeat = (1. - rightCell.acsr.injl.bet) * rightCell.acsr.injl.FluidoPro.CalorLiq(meanPressure, meanTemperature) + rightCell.acsr.injl.bet * rightCell.acsr.injl.fluidocol.CalorLiq(meanPressure, meanTemperature);
    } else if (rightCell.acsr.tipo == kAccessoryInflowPerformance) { // IPR case
        sourceTemperature = rightCell.acsr.ipr.Tres;
        sourceGasSpecificHeat = rightCell.acsr.ipr.FluidoPro.CalorGas(meanPressure, meanTemperature);
        sourceSpecificHeatRatio = rightCell.acsr.ipr.FluidoPro.ConstAdG(meanPressure, meanTemperature);
        sourceLiquidSpecificHeat = rightCell.acsr.ipr.FluidoPro.CalorLiq(meanPressure, meanTemperature);
    } else if (rightCell.acsr.tipo == kAccessoryRadialPorous) { // IPR case
        sourceTemperature = rightCell.acsr.radialPoro.tRes;
        sourceGasSpecificHeat = rightCell.acsr.radialPoro.flup.CalorGas(meanPressure, meanTemperature);
        sourceSpecificHeatRatio = rightCell.acsr.radialPoro.flup.ConstAdG(meanPressure, meanTemperature);
        sourceLiquidSpecificHeat = rightCell.acsr.radialPoro.flup.CalorLiq(meanPressure, meanTemperature);
    } else if (rightCell.acsr.tipo == kAccessoryPorous2D) { // IPR case
        sourceTemperature = rightCell.acsr.poroso2D.dados.tRes;
        sourceGasSpecificHeat = rightCell.acsr.poroso2D.dados.flup.CalorGas(meanPressure, meanTemperature);
        sourceSpecificHeatRatio = rightCell.acsr.poroso2D.dados.flup.ConstAdG(meanPressure, meanTemperature);
        sourceLiquidSpecificHeat = rightCell.acsr.poroso2D.dados.flup.CalorLiq(meanPressure, meanTemperature);
    } else if (rightCell.acsr.tipo == kAccessoryLeak && rightCell.acsr.fontechk.abertura > 1e-6 &&
               (rightCell.fontemassCR + rightCell.fontemassGR + rightCell.fontemassLR) > 1e-9) {
        // leak case
        sourceTemperature = rightCell.acsr.fontechk.tamb;
        sourceGasSpecificHeat = rightCell.acsr.fontechk.fluidoPamb.CalorGas(meanPressure, meanTemperature);
        sourceSpecificHeatRatio = rightCell.acsr.fontechk.fluidoPamb.ConstAdG(meanPressure, meanTemperature);
        sourceLiquidSpecificHeat = (1. - rightCell.acsr.fontechk.betISamb) *
                   rightCell.acsr.fontechk.fluidoPamb.CalorLiq(meanPressure, meanTemperature) +
               rightCell.acsr.fontechk.betISamb * rightCell.acsr.fontechk.fluidocol.CalorLiq(meanPressure, meanTemperature);
    } else if (rightCell.acsr.tipo != kAccessoryNone) {

        clearSourceSpecificHeats(sourceGasSpecificHeat, sourceSpecificHeatRatio,
                                 sourceLiquidSpecificHeat);
    } else {
        clearSourceSpecificHeats(sourceGasSpecificHeat, sourceSpecificHeatRatio,
                                 sourceLiquidSpecificHeat);
    }

    liquidMassSourceTerm = 0;
    if (rightCell.fontemassLR > 0.)
        liquidMassSourceTerm = rightCell.fontemassLR / cellLength;
    if (rightCell.fontemassCR > 0.)
        liquidMassSourceTerm += rightCell.fontemassCR / cellLength;
    liquidMassSourceTerm *= (sourceLiquidSpecificHeat) * (sourceTemperature - rightCell.temp);

    gasMassSourceTerm = rightCell.fontemassGR / cellLength;
    if (gasMassSourceTerm > 0.)
        gasMassSourceTerm *= (sourceGasSpecificHeat / sourceSpecificHeatRatio) * (sourceTemperature - rightCell.temp);
    else
        gasMassSourceTerm = 0;
    return TemperatureSourceTerms{.gas = gasMassSourceTerm,
                                  .liquid = liquidMassSourceTerm};
}

/// Kinetic-energy term of the REVERSE march.
double computeReverseSteadyKineticTerm(
    const ThermalState &state, int cellIndex,
    double meanSuperficialGasVelocity,
    double meanSuperficialLiquidVelocity) {
    Cel &cell = state.cells[cellIndex];
    Cel &rightCell = state.cells[cellIndex + 1];
    double kineticTerm = 0;
    double upstreamMeanGasVelocity = 0;
    double upstreamMeanLiquidVelocity = 0;
    double meanGasVelocity = 0;
    double meanLiquidVelocity = 0;
    // kinetic energy term:
    if (cell.acsr.tipo == kAccessoryNone && rightCell.acsr.tipo == kAccessoryNone && cellIndex < state.lastCell) {
        double kineticCellLength = rightCell.dx;
        double upstreamDiameter = rightCell.duto.a;
        double upstreamFlowArea = 0.25 * M_PI * upstreamDiameter * upstreamDiameter;

        if (rightCell.alf > 1e-3)
            meanGasVelocity = meanSuperficialGasVelocity / rightCell.alf;
        if (rightCell.alf < (1. - 1e-3))
            meanLiquidVelocity = meanSuperficialLiquidVelocity / (1. - rightCell.alf);

        if (cell.alf > 1e-3) {
            upstreamMeanGasVelocity = fabs(cell.QG) / (upstreamFlowArea);
            upstreamMeanGasVelocity /= cell.alf;
        }
        if (cell.alf < (1. - 1e-3)) {
            upstreamMeanLiquidVelocity = fabs(cell.QL) / (upstreamFlowArea);
            upstreamMeanLiquidVelocity /= (1. - cell.alf);
        }

        kineticTerm = -fabs(cell.MC - cell.Mliqini) * meanGasVelocity * (meanGasVelocity - upstreamMeanGasVelocity) / kineticCellLength -
                   fabs(cell.Mliqini) * meanLiquidVelocity * (meanLiquidVelocity - upstreamMeanLiquidVelocity) / kineticCellLength;
    }
    return kineticTerm;
}

/// Latent-heat term of the REVERSE march. Unlike the direct march, it does not
/// sanitise NaN or clamp by limTransMass.
double computeReverseSteadyLatentHeatTerm(
    const ThermalState &state, int cellIndex, double flowArea,
    double meanPressure, double meanTemperature,
    double meanSuperficialGasVelocity,
    double meanSuperficialLiquidVelocity) {
    varGlob1D &globals = *state.globals;
    Cel &cell = state.cells[cellIndex];
    Cel &rightCell = state.cells[cellIndex + 1];
    double latentHeatTerm = 0.;
    if (cell.flui.dVaporMassFraction < (1 - 1e-15) && cell.flui.dVaporMassFraction > (1e-15)) {
        if (rightCell.acsr.tipo == kAccessoryGasInjection || rightCell.acsr.tipo == kAccessoryLiquidInjection || rightCell.acsr.tipo == kAccessoryInflowPerformance || rightCell.acsr.tipo == kAccessoryRadialPorous || rightCell.acsr.tipo == kAccessoryPorous2D)
            rightCell.FonteMudaFase =
                0.;
        if (state.latentHeatEnabled > 0 && state.steadyIteration != 0 && state.input.flashCompleto == 0) {
            latentHeatTerm = -interpolateLatentHeat(state, meanPressure, meanTemperature) * rightCell.FonteMudaFase;
        } else if ((state.input.flashCompleto == 1 || state.input.flashCompleto == 2) && state.latentHeatEnabled > 0 && state.steadyIteration != 0) {
            latentHeatTerm = -(cell.flui.EntalpGas(meanPressure, meanTemperature) -
                        cell.flui.EntalpLiq(meanPressure, meanTemperature)) *
                      rightCell.FonteMudaFase;
        } else
            latentHeatTerm = 0;
    }

    ////// useless section////////////////////////////////////////////////////////
    double interfaceVoidFraction;
    double leftInterfaceVoidFraction;
    if (meanSuperficialGasVelocity > 0) {
        interfaceVoidFraction = rightCell.alf;
        leftInterfaceVoidFraction = rightCell.alfL;
    } else {
        interfaceVoidFraction = rightCell.alfR;
        leftInterfaceVoidFraction = rightCell.alf;
    }
    double slipVelocity;
    if (interfaceVoidFraction > globals.localtiny && interfaceVoidFraction < (1. - globals.localtiny))
        slipVelocity = meanSuperficialGasVelocity / interfaceVoidFraction - meanSuperficialLiquidVelocity / (1. - interfaceVoidFraction);
    else if (interfaceVoidFraction > globals.localtiny)
        slipVelocity = meanSuperficialGasVelocity;
    else
        slipVelocity = meanSuperficialLiquidVelocity;
    [[maybe_unused]] double interfacialWorkTerm = flowArea * rightCell.pres * kPascalPerKgfPerCm2Coarse * slipVelocity * (interfaceVoidFraction - leftInterfaceVoidFraction) / rightCell.dx;
    return latentHeatTerm;
}

/// Annulus coupling for the REVERSE march: faces cellIndex + 1 and never
/// consults the production network resistance.
double applyReverseSteadyAnnulusCoupling(
    const ThermalState &state, int cellIndex,
    double interfaceMeanPressure, double interfaceMeanTemperature,
    double gasSpecificHeat, double gasDensity) {
    Cel &rightCell = state.cells[cellIndex + 1];
    double annulusResistance = 0.;
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
            state.gasCells[j].calor.kint = rightCell.flui.CondGas(interfaceMeanPressure, interfaceMeanTemperature);
            state.gasCells[j].calor.cpint = gasSpecificHeat;
            state.gasCells[j].calor.rhoint = gasDensity;
            state.gasCells[j].calor.viscint = rightCell.flui.ViscGas(interfaceMeanPressure, interfaceMeanTemperature) * kPascalSecondPerCentipoise;
            state.gasCells[j].fluxcal = state.gasCells[j].calor.transperm(); // heat exchange in the annulus
            annulusResistance = state.gasCells[j].calor.resGlob;                // wall resistance
            // Only the thermal resistance from the casing to the formation is needed, not the actual heat flow.
            rightCell.calor.Vextern1 = 100.;
            rightCell.calor.kextern1 = state.gasCells[j].calor.kint;
            rightCell.calor.cpextern1 = state.gasCells[j].calor.cpint;
            rightCell.calor.rhoextern1 = state.gasCells[j].calor.rhoint;
            rightCell.calor.viscextern1 = state.gasCells[j].calor.viscint;
        } else { // After the first iteration, consider only heat transfer between the tubing and annular gas.
            // Through the tubing wall.
            annulusResistance = 0.;
            rightCell.calor.Vextern1 = state.gasCells[j].VGasR / state.gasCells[j].u1L;
            rightCell.calor.kextern1 = state.gasCells[j].calor.kint;
            rightCell.calor.cpextern1 = state.gasCells[j].calor.cpint;
            rightCell.calor.rhoextern1 = state.gasCells[j].calor.rhoint;
            rightCell.calor.viscextern1 = state.gasCells[j].calor.viscint;
        }
        if (state.steadyIteration == 0)
            rightCell.calor.Textern1 = state.gasCells[j].calor.Textern1; // on the first iteration, as the calculation
        // uses the whole well resistance, so the external temperature is the geothermal one
        else
            rightCell.calor.Textern1 = state.gasCells[j].temp; // on the following iterations, the ambient temperature
        // is the gas temperature
    }
    return annulusResistance;
}

}  // namespace

void advanceReverseSteadyTemperature(const ThermalState &state, int cellIndex, int rungeKuttaStage) {
    Cel &cell = state.cells[cellIndex];
    varGlob1D &globals = *state.globals;
    Cel &rightCell = state.cells[cellIndex + 1];
    double cellLength = 0.5 * (cell.dx + rightCell.dx);
    double meanCellLength = 0.5 * (cell.dx + rightCell.dx);
    double meanTemperatureGradientCorrection = (cell.dx * cell.dTdLCor + rightCell.dx * rightCell.dTdLCor) / (cell.dx + rightCell.dx);
    double diameter = rightCell.duto.a;
    double flowArea = 0.25 * M_PI * diameter * diameter;
    double meanVoidFraction;
    double betmed;
    double meanPressure;
    double meanTemperature;
    if (rungeKuttaStage == 0) {
        meanVoidFraction = rightCell.alf;
        betmed = rightCell.bet;
        meanPressure = rightCell.pres;
        meanTemperature = rightCell.temp;
    } else {
        meanVoidFraction = cell.alf;
        betmed = cell.bet;
        meanPressure = cell.pres;
        meanTemperature = cell.temp;
    }
    double meanSuperficialGasVelocity;
    meanSuperficialGasVelocity = fabs(rightCell.QG) / flowArea; // Superficial gas velocity
    double meanSuperficialLiquidVelocity;
    meanSuperficialLiquidVelocity = fabs(rightCell.QL) / flowArea; // Superficial liquid velocity
    double mixtureFluxSign = 1.;
    if (fabs(meanSuperficialGasVelocity + meanSuperficialLiquidVelocity) > state.slowHeatTransferThreshold) {
        // Perform the thermal calculation for mixture velocities above 0.1 m/s.
        // Otherwise, assume the fluid temperature equals the ambient temperature.
        mixtureFluxSign = (meanSuperficialGasVelocity + meanSuperficialLiquidVelocity) / fabs(meanSuperficialGasVelocity + meanSuperficialLiquidVelocity);
        double interfaceMeanPressure = rightCell.presaux - cell.dpB / kPascalPerKgfPerCm2;
        double interfaceMeanTemperature;
        if (state.steadyIteration != 0 && state.input.AceleraConvergPerm == 0)
            // Set the left interface temperature from the previous iteration.
            // Use the adjacent left cell value to accelerate convergence.
            interfaceMeanTemperature = (cell.dx * cell.temp + cell.dxR * cell.tempR) / (cell.dx + cell.dxR);
        else
            interfaceMeanTemperature = rightCell.temp;
        double rp = rightCell.flui.MasEspLiq(meanPressure, meanTemperature);
        double rc = rightCell.fluicol.MasEspFlu(meanPressure, meanTemperature);
        double liquidDensity = (1. - betmed) * rp + betmed * rc;
        double gasDensity = rightCell.flui.MasEspGas(meanPressure, meanTemperature);
        double liquidSpecificHeat = (1. - betmed) * rightCell.flui.CalorLiq(interfaceMeanPressure, interfaceMeanTemperature) + betmed * rightCell.fluicol.CalorLiq(interfaceMeanPressure, interfaceMeanTemperature);
        double gasSpecificHeat = rightCell.flui.CalorGas(interfaceMeanPressure, interfaceMeanTemperature);
        // Liquid Joule-Thomson coefficient multiplied by cp.
        // Gas Joule-Thomson coefficient multiplied by cp.p
        double liquidJouleThomson = (1. - betmed) * rightCell.flui.JTL(interfaceMeanPressure, interfaceMeanTemperature) - betmed / rc;
        if (state.input.pocinjec > 0 && state.input.condpocinj.tipoFlui == 2) {
            liquidJouleThomson = -(1 + (interfaceMeanTemperature + 273.14) * rightCell.fluicol.DrhoDtFlu(interfaceMeanPressure, interfaceMeanTemperature) / rc) / rc;
        }
        double gasJouleThomson = rightCell.flui.JTG(interfaceMeanPressure, interfaceMeanTemperature);
        // potential energy:
        double hydrostaticPower = -(liquidDensity * meanSuperficialLiquidVelocity + gasDensity * meanSuperficialGasVelocity) * flowArea * kGravity * sin(rightCell.duto.teta);


        // Definition of internal parameters in the piping to achieve heat exchange with the environment
        rightCell.calor.Tint = interfaceMeanTemperature;
        rightCell.calor.Vint = fabs(meanSuperficialGasVelocity + meanSuperficialLiquidVelocity);
        double liquidConductivity = (1. - betmed) * rightCell.flui.CondLiq(interfaceMeanPressure, interfaceMeanTemperature) + betmed * rightCell.fluicol.CondLiq(interfaceMeanPressure, interfaceMeanTemperature);
        rightCell.calor.kint = liquidConductivity * (1 - meanVoidFraction) + rightCell.flui.CondGas(interfaceMeanPressure, interfaceMeanTemperature) * meanVoidFraction;
        rightCell.calor.cpint = liquidSpecificHeat * (1 - meanVoidFraction) + gasSpecificHeat * meanVoidFraction;
        rightCell.calor.rhoint = liquidDensity * (1 - meanVoidFraction) + gasDensity * meanVoidFraction;
        double liquidViscosity = (1. - betmed) * rightCell.flui.ViscOleo(interfaceMeanPressure, interfaceMeanTemperature) + betmed * rightCell.fluicol.VisFlu(interfaceMeanPressure, interfaceMeanTemperature);
        rightCell.calor.viscint = liquidViscosity * (1 - meanVoidFraction) * kPascalSecondPerCentipoise + rightCell.flui.ViscGas(interfaceMeanPressure, interfaceMeanTemperature) * meanVoidFraction * kPascalSecondPerCentipoise;

        [[maybe_unused]] double relaxedHeatFlux = 0;
        if (state.steadyIteration > 0)
            relaxedHeatFlux = rightCell.fluxcalmed;
        double annulusResistance = 0.;
        double heatFlux;
        [[maybe_unused]] double gasHeatFlux;
        annulusResistance = applyReverseSteadyAnnulusCoupling(
            state, cellIndex, interfaceMeanPressure, interfaceMeanTemperature,
            gasSpecificHeat, gasDensity);
        rightCell.fluxcalmed = 0;
        heatFlux = mixtureFluxSign * rightCell.calor.transperm(annulusResistance);
        rightCell.fluxcalmed = heatFlux; // heat flux in the tubing

        double temperatureSpatialCoefficient = (liquidDensity * meanSuperficialLiquidVelocity * liquidSpecificHeat + gasDensity * meanSuperficialGasVelocity * gasSpecificHeat) * flowArea; // multiplying term
        // the derivative Dt/Dx
        double pressureSpatialCoefficient = 1. * (liquidDensity * meanSuperficialLiquidVelocity * liquidJouleThomson + gasDensity * meanSuperficialGasVelocity * gasJouleThomson) * flowArea; // multiplying term
        // the derivative Dp/Dx
        double pressureGradient;
        if ((rightCell.acsr.tipo != kAccessoryPump || rightCell.acsr.bcs.freq < 1) && rightCell.acsr.tipo != 7)
            // without an ESP or a pressure increment, the pressure at the left boundary
            // and the pressure at the cell centre are used to compute Dp/Dx
            pressureGradient = 2. * (rightCell.presaux - rightCell.pres) * kPascalPerKgfPerCm2 / rightCell.dx;
        else {
            // with an ESP or a pressure increment, the pressure of the cell to the left
            // and the pressure at the cell centre are used to compute Dp/Dx
            pressureGradient = (interfaceMeanPressure - rightCell.pres) * kPascalPerKgfPerCm2Coarse / cellLength;
        }
        cell.VTemper = meanSuperficialLiquidVelocity; // this velocity is only useful in the transient case
        // only so that there is a value when the transient simulation starts
        [[maybe_unused]] double temperatureGradient = (-rightCell.temp) / meanCellLength;

        double kineticTerm = computeReverseSteadyKineticTerm(
            state, cellIndex, meanSuperficialGasVelocity,
            meanSuperficialLiquidVelocity);

        TemperatureSourceTerms reverseSources = computeReverseSteadySourceTerms(
            state, cellIndex, meanPressure, meanTemperature, cellLength);
        double gasMassSourceTerm = reverseSources.gas;
        double liquidMassSourceTerm = reverseSources.liquid;

        // latent heat effect, when requested
        double latentHeatTerm = computeReverseSteadyLatentHeatTerm(
            state, cellIndex, flowArea, meanPressure, meanTemperature,
            meanSuperficialGasVelocity, meanSuperficialLiquidVelocity);

        if (fabs(temperatureSpatialCoefficient) > globals.localtiny) {
            // energy term from boundary work, potential energy,
            // kinetic energy, mass sources, latent heat and shaft work:
            double sourceTemperatureGradient = meanTemperatureGradientCorrection * (pressureSpatialCoefficient * pressureGradient - kineticTerm - (hydrostaticPower) + (liquidMassSourceTerm + gasMassSourceTerm) - latentHeatTerm - rightCell.potBT / meanCellLength) / temperatureSpatialCoefficient;
            // energy term from the heat exchange
            double heatFluxTemperatureGradient = meanTemperatureGradientCorrection * (heatFlux) / temperatureSpatialCoefficient;

            // checks whether the heat exchange is happening too fast:
            // as the temperature advance is explicit, this can make the thermal calculation unstable;
            // if the heat exchange is found to be fast, the advance is
            // made in a larger number of steps from one cell to the next
            int subStepCount;
            double subStepLength;
            double stabilityLength = fabs(temperatureSpatialCoefficient) / meanTemperatureGradientCorrection;
            if (meanCellLength / (rightCell.calor.resGlob + annulusResistance) < (stabilityLength + 0. * 1000.)) {
                subStepCount = 1;
                subStepLength = meanCellLength;
            } else { // small thermal resistance: determines in how many steps
                // the temperature advance is made
                subStepCount = (meanCellLength / (rightCell.calor.resGlob + annulusResistance)) / (stabilityLength + 0 * 1000) + 1;
                subStepLength = meanCellLength / subStepCount; // cell divided into npassos
            }
            double subStepTemperature = rightCell.temp;

            subStepTemperature = subStepLength * (-(-rightCell.temp) / subStepLength + sourceTemperatureGradient + heatFluxTemperatureGradient); // first step
            for (int j = 1; j < subStepCount; j++) {                                                  // following steps
                // the energy terms are kept, except the heat flux, which is
                // re-evaluated at each step:
                rightCell.calor.Tint = subStepTemperature;
                rightCell.calor.Vint = meanSuperficialGasVelocity + meanSuperficialLiquidVelocity;
                liquidConductivity = (1. - betmed) * rightCell.flui.CondLiq(meanPressure, subStepTemperature) + betmed * rightCell.fluicol.CondLiq(meanPressure, subStepTemperature);
                rightCell.calor.kint = liquidConductivity * (1 - meanVoidFraction) + rightCell.flui.CondGas(meanPressure, subStepTemperature) * meanVoidFraction;
                rightCell.calor.cpint = liquidSpecificHeat * (1 - meanVoidFraction) + gasSpecificHeat * meanVoidFraction;
                rightCell.calor.rhoint = liquidDensity * (1 - meanVoidFraction) + gasDensity * meanVoidFraction;
                liquidViscosity = (1. - betmed) * rightCell.flui.ViscOleo(meanPressure, subStepTemperature) + betmed * rightCell.fluicol.VisFlu(meanPressure, subStepTemperature);
                rightCell.calor.viscint = liquidViscosity * (1 - meanVoidFraction) * kPascalSecondPerCentipoise + rightCell.flui.ViscGas(meanPressure, subStepTemperature) * meanVoidFraction * kPascalSecondPerCentipoise;
                if (state.steadyIteration != 0 && state.input.lingas == 1 && (cellIndex + 1 <= state.annulusTubingStart && cellIndex + 1 >= state.annulusTubingEnd)) {
                    int k = state.annulusTubingStart + state.tubingAnnulusStart - (cellIndex + 1);
                    double externalTemperatureStep = (state.gasCells[k].temp - state.gasCells[k - 1].temp) / subStepCount;
                    rightCell.calor.Textern1 = state.gasCells[k].temp + (j - 1) * externalTemperatureStep;
                }
                heatFlux = mixtureFluxSign * rightCell.calor.transperm(annulusResistance);
                rightCell.fluxcalmed += heatFlux;
                heatFluxTemperatureGradient = meanTemperatureGradientCorrection * (heatFlux) / temperatureSpatialCoefficient;
                subStepTemperature = subStepLength * (-(-subStepTemperature) / subStepLength + sourceTemperatureGradient + heatFluxTemperatureGradient);
            }

            cell.temp = subStepTemperature;
            rightCell.fluxcalmed /= subStepCount;
            rightCell.fluxcalmed = 1. * rightCell.fluxcalmed;

            if (cell.temp < kMinimumTemperatureCelsius)
                cell.temp = kMinimumTemperatureCelsius;
            if (cell.temp > kMaximumTemperatureCelsius)
                cell.temp = kMaximumTemperatureCelsius;
        } else
            cell.temp = cell.calor.Textern1;
    } else { // case where the mixture velocity is very low
        cell.temp = rightCell.calor.Textern1;
        if (state.input.lingas == 1 && (cellIndex + 1 <= state.annulusTubingStart && cellIndex + 1 >= state.annulusTubingEnd)) {
            int j = state.annulusTubingStart + state.tubingAnnulusStart - (cellIndex + 1);
            rightCell.calor.Vextern1 = 100.;
            rightCell.calor.kextern1 = state.gasCells[j].calor.kint;
            rightCell.calor.cpextern1 = state.gasCells[j].calor.cpint;
            rightCell.calor.rhoextern1 = state.gasCells[j].calor.rhoint;
            rightCell.calor.viscextern1 = state.gasCells[j].calor.viscint;
            [[maybe_unused]] double heatFlux = rightCell.calor.transperm(0);
        }
    }
}

void computeGasTemperature(const ThermalState &state, int cellIndex, double previousTemperature, int steadyMode) {
    const bool hasRightGasCell = cellIndex < state.gasCellCount;

    if (state.thermalSourceDisabled == 0) {
        double cellLength = state.gasCells[cellIndex].dx0;
        double meanCellLength = 0.5 * (state.gasCells[cellIndex].dx0 + state.gasCells[cellIndex - 1].dx0);
        double flowArea = state.gasCells[cellIndex].duto.area;
        double meanSuperficialGasVelocity;
        if (hasRightGasCell)
            meanSuperficialGasVelocity = state.gasCells[cellIndex].VGasR / state.gasCells[cellIndex].u1L;
        else {
            meanSuperficialGasVelocity = state.gasCells[cellIndex].VGasL / state.gasCells[cellIndex].u1L;
        }
        double gasDensity = state.gasCells[cellIndex].rg;
        double gasSpecificHeat = state.gasCells[cellIndex].flui.CalorGas(state.gasCells[cellIndex].presini, state.gasCells[cellIndex].tempini);
        double gasSpecificHeatConstantVolume = state.gasCells[cellIndex].flui.CalorGasVolMod(state.gasCells[cellIndex].presini, state.gasCells[cellIndex].tempini);
        double gasJouleThomson = state.gasCells[cellIndex].flui.JTG(state.gasCells[cellIndex].presini, state.gasCells[cellIndex].tempini);
        double hydrostaticPower = (gasDensity * meanSuperficialGasVelocity) * flowArea * kGravity * sin(state.gasCells[cellIndex].duto.teta);

        state.gasCells[cellIndex].calor.Tint = state.gasCells[cellIndex].tempini;
        state.gasCells[cellIndex].calor.dtL = state.gasCells[cellIndex].tempini - state.gasCells[cellIndex - 1].tempini;
        state.gasCells[cellIndex].calor.Vint = meanSuperficialGasVelocity;
        state.gasCells[cellIndex].calor.dt = state.gasCells[cellIndex].dt;
        state.gasCells[cellIndex].calor.kint = state.gasCells[cellIndex].flui.CondGas(state.gasCells[cellIndex].presini, state.gasCells[cellIndex].tempini);
        state.gasCells[cellIndex].calor.cpint = gasSpecificHeat;
        state.gasCells[cellIndex].calor.rhoint = gasDensity;
        state.gasCells[cellIndex].calor.viscint = state.gasCells[cellIndex].flui.ViscGas(state.gasCells[cellIndex].presini, state.gasCells[cellIndex].tempini) * kPascalSecondPerCentipoise;
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
        if (hasRightGasCell)
            lengthRatio = cellLength / (cellLength + state.gasCells[cellIndex + 1].dx0);
        else
            lengthRatio = cellLength / (cellLength + state.gasCells[cellIndex - 1].dx0);
        double timeCoefficient = gasDensity * gasSpecificHeatConstantVolume * flowArea;
        double pressureTimeCoefficient = state.gasCells[cellIndex].flui.CalorGasPresMod(state.gasCells[cellIndex].presini, state.gasCells[cellIndex].tempini, state.gasCells[cellIndex].rg) *
                               (gasDensity * flowArea);

        double temperatureSpatialCoefficient = gasDensity * meanSuperficialGasVelocity * gasSpecificHeat * flowArea;
        double pressureSpatialCoefficient = gasDensity * meanSuperficialGasVelocity * gasJouleThomson * flowArea;
        double pressureGradient;
        if (hasRightGasCell)
            pressureGradient = 2. * (((1 - lengthRatio) * state.gasCells[cellIndex + 1].presini + lengthRatio * state.gasCells[cellIndex].presini) - state.gasCells[cellIndex].presini) * kPascalPerKgfPerCm2 / cellLength;
        else
            pressureGradient = 2. * (state.gasCells[cellIndex].presini - ((1 - lengthRatio) * state.gasCells[cellIndex - 1].presini + lengthRatio * state.gasCells[cellIndex].presini)) * kPascalPerKgfPerCm2 / cellLength;
        double temperatureGradient = (state.gasCells[cellIndex].tempini - state.gasCells[cellIndex - 1].tempini) / meanCellLength;
        if (hasRightGasCell)
            if (meanSuperficialGasVelocity < 0)
                temperatureGradient = (state.gasCells[cellIndex + 1].tempini - state.gasCells[cellIndex].tempini) / meanCellLength;
        if ((cellIndex == 1 && meanSuperficialGasVelocity <= 0) || (cellIndex == state.gasCellCount && meanSuperficialGasVelocity <= 0))
            temperatureGradient = 0.;

        double kineticTerm;
        double mixtureFluxDifference = 0.;
        double mixtureDensity = gasDensity;
        double leftFaceMeanSuperficialGasVelocity = state.gasCells[cellIndex].VGasL / state.gasCells[cellIndex - 1].u1L;
        mixtureFluxDifference = (meanSuperficialGasVelocity - leftFaceMeanSuperficialGasVelocity) / cellLength;

        kineticTerm = mixtureDensity * flowArea * meanSuperficialGasVelocity * meanSuperficialGasVelocity * mixtureFluxDifference;

        double gasMassSourceTerm = 0.;
        double liquidMassSourceTerm = 0.;

        double pressureWorkFactor = 1.;
        if ((*state.globals).lixo5 < 1000.)
            pressureWorkFactor = 1.;

        state.gasCells[cellIndex].temp = ((timeCoefficient / state.gasCells[cellIndex].dt) * state.gasCells[cellIndex].temp - (pressureWorkFactor) * (pressureTimeCoefficient * (state.gasCells[cellIndex].pres - state.gasCells[cellIndex].presini) * kPascalPerKgfPerCm2 / state.gasCells[cellIndex].dt) + state.gasCells[cellIndex].dTdLCor * (-temperatureSpatialCoefficient * temperatureGradient + pressureSpatialCoefficient * pressureGradient - kineticTerm - hydrostaticPower + liquidMassSourceTerm + gasMassSourceTerm + state.gasCells[cellIndex].fluxcal)) / (timeCoefficient / state.gasCells[cellIndex].dt);

        if (state.gasCells[cellIndex].temp < kMinimumTemperatureCelsius)
            state.gasCells[cellIndex].temp = kMinimumTemperatureCelsius;
        if (state.gasCells[cellIndex].temp > kMaximumTemperatureCelsius)
            state.gasCells[cellIndex].temp = kMaximumTemperatureCelsius;

        if (cellIndex > 0)
            state.gasCells[cellIndex - 1].tempR = state.gasCells[cellIndex].temp;
        if (hasRightGasCell)
            state.gasCells[cellIndex + 1].tempL = state.gasCells[cellIndex].temp;
    } else {
        state.gasCells[cellIndex].temp = state.gasCells[cellIndex].calor.Textern1;

        if (cellIndex > 0)
            state.gasCells[cellIndex - 1].tempR = state.gasCells[cellIndex].temp;
        if (hasRightGasCell)
            state.gasCells[cellIndex + 1].tempL = state.gasCells[cellIndex].temp;
    }
}

void computeDischargeTemperature(const ThermalState &state, int cellIndex) {
    Cel &leftCell = state.cells[cellIndex - 1];
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
    double pressure;
    double temperature;

    pressure = state.gasCells[cellIndex - 1].pres;
    temperature = state.gasCells[cellIndex - 1].temp;
    const double completionFluidDensity = state.gasCells[cellIndex].MasEspFlu(pressure, temperature);
    double mixtureDensity = ((leftLiquidLength + rightLiquidLength) * completionFluidDensity + (leftGasLength + rightGasLength) * state.gasCells[cellIndex].flui.MasEspGas(pressure, temperature)) / totalLength;

    double interfaceVelocity = state.gasCells[cellIndex].VGasL / (completionFluidDensity * state.gasCells[cellIndex - 1].duto.area);
    if (state.gasCells[cellIndex].razInter > (*state.globals).localtiny)
        interfaceVelocity = state.gasCells[cellIndex].VGasL / (state.gasCells[cellIndex].flui.MasEspGas(pressure, temperature) * state.gasCells[cellIndex - 1].duto.area);
    double liquidSpecificHeat = ((leftLiquidLength + rightLiquidLength) * state.gasCells[cellIndex].CalorLiq(pressure, temperature) + (leftGasLength + rightGasLength) * state.gasCells[cellIndex].flui.CalorGas(pressure, temperature)) / totalLength;

    state.gasCells[cellIndex - 1].calor.Tint = temperature;
    state.gasCells[cellIndex - 1].calor.Vint = interfaceVelocity;
    double liquidConductivity = ((leftLiquidLength + rightLiquidLength) * state.gasCells[cellIndex - 1].CondLiq(pressure, temperature) + (leftGasLength + rightGasLength) * state.gasCells[cellIndex].flui.CondGas(pressure, temperature)) / totalLength;
    state.gasCells[cellIndex - 1].calor.kint = liquidConductivity;
    state.gasCells[cellIndex - 1].calor.cpint = liquidSpecificHeat;
    leftCell.calor.rhoint = mixtureDensity;
    double liquidViscosity = (((leftLiquidLength + rightLiquidLength) * state.gasCells[cellIndex].VisFlu(pressure, temperature) + (leftGasLength + rightGasLength) * state.gasCells[cellIndex].flui.ViscGas(pressure, temperature)) * 1e-3) / totalLength;
    leftCell.calor.viscint = liquidViscosity;

    double heatFlux = leftCell.calor.transtrans();
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
        double gasSpecificHeat = state.gasLiftChokes[valveIndex].flui.CalorGas(stageInletPressure, stageInletTemperature);
        double gasJouleThomson = state.gasLiftChokes[valveIndex].flui.JTG(stageInletPressure, stageInletTemperature) / gasSpecificHeat;
        stageInletTemperature -= gasJouleThomson * (state.gasLiftChokes[valveIndex].presEstag - state.gasLiftChokes[valveIndex].presGarg) * kPascalPerKgfPerCm2Variant;
    }
    return stageInletTemperature;
}

void updateProductionTemperaturePeriphery(const ThermalState &state, int cellIndex) {
    Cel &cell = state.cells[cellIndex];
    Cel &leftCell = state.cells[cellIndex - 1];
    Cel &rightCell = state.cells[cellIndex + 1];
    if (cellIndex > 0)
        leftCell.tempR = cell.temp;
    if (cellIndex < state.lastCell)
        rightCell.tempL = cell.temp;
    cell.tempini = cell.temp;
}

void computeOutletTemperature(const ThermalState &state) {

    if (state.input.chokep.abertura[0] <= 0.6 && state.input.chokep.abertura[0] > (*state.globals).localtiny && state.outletPressure < state.gasSurfacePressure) {
        double inletMassFlowRate = state.cells[state.lastCell - 1].MR;
        double gasMassFlowRate = state.cells[state.lastCell - 1].MR - state.cells[state.lastCell - 1].MliqiniR;
        double rholc = state.cells[state.lastCell].fluicol.MasEspFlu(state.cells[state.lastCell].pres, state.cells[state.lastCell].temp);
        double betEF = state.cells[state.lastCell].bet;
        double quality = fabs(gasMassFlowRate / inletMassFlowRate);

        double liquidJouleThomson = (1. - betEF) * state.cells[state.lastCell].flui.JTL(state.cells[state.lastCell].pres, state.cells[state.lastCell].temp) - betEF / rholc; // change 2
        double gasJouleThomson = state.cells[state.lastCell].flui.JTG(state.cells[state.lastCell].pres, state.cells[state.lastCell].temp);
        state.surfaceTemperature = state.cells[state.lastCell].temp + ((1. - quality) * liquidJouleThomson + quality * gasJouleThomson) * (state.cells[state.lastCell].pres - state.cells[state.lastCell].pres); //????????
                                                                                                                  //???????????????????????????????celula[ncel].pres - celula[ncel].pres????????????????????????????????????????????

    } else
        state.surfaceTemperature = state.cells[state.lastCell - 1].temp;
}

}  // namespace sisprod::thermal
