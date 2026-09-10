#include "SisProdGasLift.h"

#include "Leitura.h"
#include "Matriz.h"
#include "Vetor.h"
#include "celulaGas.h"
#include "celula3.h"
#include "chokegas.h"

#include <math.h>

namespace sisprod::gaslift {

void computeGasUnloadingHydrostatics(const GasLiftState &state) {
    state.gasCells[0].massfonteCH = 0;
    double pmed;
    double tmed;

    if (state.input.gasinj.tipoCC == 0 && state.input.controDesc == 0)
        pmed = state.input.gasinj.presinj[0];
    else {
        pmed = 10.;
        if (state.input.controDesc == 1) {
            pmed = state.input.presIniDescG;
            state.initialGasPressure = state.input.presIniDescG;
        }
    }

    state.gasCells[0].presL = pmed;
    state.gasCells[0].pres = pmed;
    state.gasCells[0].presini = pmed;
    state.gasCells[1].presL = pmed;
    tmed = state.gasCells[0].calor.Textern1;
    state.gasCells[0].tempL = tmed;
    state.gasCells[0].temp = tmed;
    state.gasCells[1].tempL = tmed;
    double rho0 = state.gasCells[0].flui.MasEspGas(pmed, tmed);
    double rho1;
    state.gasCells[0].u1L = state.gasCells[0].duto.area * rho0;
    state.gasCells[0].u1LL = state.gasCells[0].u1L;
    state.gasCells[1].u1LL = state.gasCells[0].u1L;
    state.gasCells[0].VGasL = 0;
    state.gasCells[0].VGasR = 0;
    state.gasCells[1].VGasL = 0;
    state.gasCells[0].massfonteCH = 0.;
    for (int gasCellIndex = 1; gasCellIndex <= state.gasCellCount; gasCellIndex++) {
        double halfUpstreamLength = 0.5 * state.gasCells[gasCellIndex].dxL;
        double halfLocalLength = 0.5 * state.gasCells[gasCellIndex].dx0;
        pmed -= rho0 * 9.81 * halfUpstreamLength * sin(state.gasCells[gasCellIndex - 1].duto.teta) / 98066.52;
        tmed = state.gasCells[gasCellIndex].calor.Textern1;
        if (gasCellIndex < state.interfaceCell)
            rho1 = state.gasCells[gasCellIndex].flui.MasEspGas(pmed, tmed);
        else
            rho1 = state.gasCells[gasCellIndex].MasEspFlu(pmed, tmed);
        pmed -= rho1 * 9.81 * halfLocalLength * sin(state.gasCells[gasCellIndex].duto.teta) / 98066.52;
        rho0 = rho1;

        state.gasCells[gasCellIndex].pres = pmed;
        state.gasCells[gasCellIndex].presini = pmed;
        state.gasCells[gasCellIndex - 1].presR = pmed;
        state.gasCells[gasCellIndex].temp = tmed;
        state.gasCells[gasCellIndex - 1].tempR = tmed;
        state.gasCells[gasCellIndex].u1L = state.gasCells[gasCellIndex].duto.area * rho0;
        state.gasCells[gasCellIndex - 1].u1R = state.gasCells[gasCellIndex].u1L;
        state.gasCells[gasCellIndex].VGasR = 0;
        state.gasCells[gasCellIndex - 1].VGasRR = 0;
        state.gasCells[gasCellIndex].massfonteCH = 0.;
        state.gasCells[gasCellIndex - 1].u1R = state.gasCells[gasCellIndex].u1L;
        if (gasCellIndex < state.gasCellCount) {
            state.gasCells[gasCellIndex + 1].presL = pmed;
            state.gasCells[gasCellIndex + 1].tempL = tmed;
            state.gasCells[gasCellIndex + 1].u1LL = state.gasCells[gasCellIndex].u1L;
            state.gasCells[gasCellIndex + 1].VGasL = 0;
        }
    }
}

void updateGasLine(const GasLiftState &state) {
    for (int gasCellIndex = 0; gasCellIndex <= state.gasCellCount; gasCellIndex++) {
        if (gasCellIndex != 0 && gasCellIndex != state.gasCellCount) {
            state.gasCells[gasCellIndex].pres = state.gasFreeTerms[3 * gasCellIndex];
            state.gasCells[gasCellIndex].presL = state.gasFreeTerms[3 * gasCellIndex - 3];
            state.gasCells[gasCellIndex].presR = state.gasFreeTerms[3 * gasCellIndex + 3];
            state.gasCells[gasCellIndex].VGasR = state.gasFreeTerms[3 * gasCellIndex + 1];
            state.gasCells[gasCellIndex].VGasL = state.gasFreeTerms[3 * gasCellIndex - 2];
            state.gasCells[gasCellIndex].VGasRR = state.gasFreeTerms[3 * gasCellIndex + 4];
        } else if (gasCellIndex == 0) {
            state.gasCells[gasCellIndex].pres = state.gasFreeTerms[3 * gasCellIndex];
            state.gasCells[gasCellIndex].presL = state.gasCells[gasCellIndex].pres;
            state.gasCells[gasCellIndex].presR = state.gasFreeTerms[3 * gasCellIndex + 3];
            state.gasCells[gasCellIndex].VGasR = state.gasFreeTerms[3 * gasCellIndex + 1];
            state.gasCells[gasCellIndex].VGasRR = state.gasFreeTerms[3 * gasCellIndex + 4];
        } else {
            state.gasCells[gasCellIndex].pres = state.gasFreeTerms[3 * gasCellIndex];
            state.gasCells[gasCellIndex].presL = state.gasFreeTerms[3 * gasCellIndex - 3];
            state.gasCells[gasCellIndex].VGasR = state.gasFreeTerms[3 * gasCellIndex + 1];
            state.gasCells[gasCellIndex].VGasL = state.gasFreeTerms[3 * gasCellIndex - 2];
            state.gasCells[gasCellIndex].presR = state.gasCells[gasCellIndex].pres;
        }
    }
}

void updateBufferedGasLine(const GasLiftState &state) {
    for (int gasCellIndex = 0; gasCellIndex <= state.gasCellCount; gasCellIndex++) {
        if (gasCellIndex != 0 && gasCellIndex != state.gasCellCount) {
            state.gasCells[gasCellIndex].VGasRBuf = state.gasFreeTerms[3 * gasCellIndex + 1];
        } else if (gasCellIndex == 0) {
            state.gasCells[gasCellIndex].VGasRBuf = state.gasFreeTerms[3 * gasCellIndex + 1];
        } else {
            state.gasCells[gasCellIndex].VGasRBuf = state.gasFreeTerms[3 * gasCellIndex + 1];
        }
    }
}

double calibratedValveArea(double calibrationPressure, double calibrationTemperature, double valveOpeningPressure, double tubingPressure,
                           double externalDiameter, double throatArea, double valveRatio, double bottomHoleTemperatureFahrenheit) {
    // Imperial throughout: pressures in psi, temperature in Fahrenheit. The
    // valve opening pressure is the casing pressure; the ratio is of areas.

    double bellowsPressureAt80F = calibrationPressure * (1 - valveRatio);
    bellowsPressureAt80F = (bellowsPressureAt80F + 14.6959488) * (80 + 460.67) / (calibrationTemperature * 1.8 + 491.67) - 14.6959488;
    double bellowsPressure = bellowsPressureAt80F * (1 + 0.00215 * (bottomHoleTemperatureFahrenheit - 80));
    double openingCriterion = valveOpeningPressure * (1 - valveRatio) + tubingPressure * valveRatio;
    double openingFraction = 0.;
    if (openingCriterion > bellowsPressure)
        openingFraction = 1.;

    double bellowsArea = throatArea / valveRatio;
    double bellowsSpringRate;
    if (externalDiameter * 100. / 2.54 > 1.1)
        bellowsSpringRate = 500.0 * bellowsArea;
    else
        bellowsSpringRate = 1950.0 * bellowsArea;
    double stemTravel = ((valveOpeningPressure - bellowsPressure) * bellowsArea - (valveOpeningPressure - tubingPressure) * throatArea) / bellowsSpringRate;

    // IF THE VALVE IS CLOSED, QG = 0
    if (stemTravel <= 0.0)
        openingFraction = 0.;

    double throatDiameter = sqrt(throatArea * 4. / M_PI);
    double throatRadius = throatDiameter / 2.0;
    double bellowsRadius = sqrt(bellowsArea / M_PI);
    double openingArea = M_PI * throatRadius * stemTravel * (stemTravel + 2.0 * sqrt(bellowsRadius * bellowsRadius - throatRadius * throatRadius));
    openingArea = openingArea / sqrt((stemTravel + sqrt(bellowsRadius * bellowsRadius - throatRadius * throatRadius)) * (stemTravel + sqrt(bellowsRadius * bellowsRadius - throatRadius * throatRadius)) + throatRadius * throatRadius);
    if (openingArea > throatArea)
        openingArea = throatArea;
    openingFraction = openingArea / throatArea;

    return openingFraction;
}

double unloadingPressureCorrection(const GasLiftState &state, double maximumFlowRate, int valveIndex, double factor, int sign) {
    int gasValveCell = state.gasValveCellIndices[valveIndex];
    state.gasLiftChokes[valveIndex].presEstag = state.gasCells[gasValveCell].pres;
    state.gasLiftChokes[valveIndex].tempEstag = state.gasCells[gasValveCell].temp;
    double rho0 = state.gasCells[gasValveCell].MasEspFlu(state.gasLiftChokes[valveIndex].presEstag, state.gasLiftChokes[valveIndex].tempEstag);
    double massica = (factor * state.input.vazDescControl - maximumFlowRate) * rho0;
    double precorr = 0.;
    precorr = pow(massica / state.gasLiftChokes[valveIndex].areagarg, 2.) / (2. * rho0 * 98066.52);
    return sign * precorr;
}

double computeUnloadingValvePressure(const GasLiftState &state, double throatFlowRate, int valveIndex) {

    double velmax = 0;
    double laz1 = 0.1;
    double laz2 = 0.4;

    double pmed;
    pmed = state.gasSurfacePressure;
    state.input.presMaxDesc = 100000.;
    for (int cellIndex = state.lastCell; cellIndex >= 0; cellIndex--) {

        double tmed = state.cells[cellIndex].temp;
        double flowArea = state.cells[cellIndex].duto.area;
        double perimeter = state.cells[cellIndex].duto.peri;
        double cellLength = state.cells[cellIndex].dx;
        double rhoG = state.cells[cellIndex].flui.MasEspGas(pmed, tmed);
        double rhoP = state.cells[cellIndex].flui.MasEspLiq(pmed, tmed);
        double rhoC = state.cells[cellIndex].fluicol.MasEspFlu(pmed, tmed);
        double viscG = state.cells[cellIndex].flui.ViscGas(pmed, tmed);
        double viscP = state.cells[cellIndex].flui.ViscOleo(pmed, tmed);
        double viscC = state.cells[cellIndex].fluicol.VisFlu(pmed, tmed);
        double voidFraction = state.cells[cellIndex].alf;
        double composition = state.cells[cellIndex].bet;
        double rholiq = composition * rhoC + (1. - composition) * rhoP;
        double viscliq = composition * viscC + (1. - composition) * viscP;
        double rhomix = voidFraction * rhoG + (1. - voidFraction) * rholiq;
        double viscmix = voidFraction * viscG + (1. - voidFraction) * viscliq;
        double vel1 = state.cells[cellIndex].QL / (flowArea * rholiq) + state.cells[cellIndex].QG / (flowArea * rhoG);
        double reynolds;
        if (fabs(vel1) > 1e-15) {
            if (state.cells[cellIndex].duto.revest == 0)
                reynolds = state.cells[cellIndex].Rey(state.cells[cellIndex].duto.a, vel1, rhomix, viscmix);
            else {
                double dhid = 4 * flowArea / perimeter;
                reynolds = state.cells[cellIndex].Rey(dhid, vel1, rhomix, viscmix);
            }
        }
        double frictionFactor;
        if (fabs(vel1) > 1e-15)
            frictionFactor = state.cells[cellIndex].fric(reynolds, state.cells[cellIndex].duto.rug / state.cells[cellIndex].duto.a);
        else
            frictionFactor = 0.;
        double tens1 = frictionFactor * rhomix * vel1 * fabs(vel1) / 2.;
        pmed -= (-9.82 * rhomix * sin(state.cells[cellIndex].duto.teta) - tens1 * perimeter / flowArea) * cellLength / 98066.5;
        if (state.cells[cellIndex].acsr.tipo == 3) {
            double auxpresmax = state.cells[cellIndex].acsr.ipr.Pres - (pmed - state.gasSurfacePressure);
            if (auxpresmax < state.input.presMaxDesc)
                state.input.presMaxDesc = auxpresmax;
        }
    }
    if (state.gasSurfacePressure >= state.input.presMaxDesc * 0.9999999 && throatFlowRate > (1 - laz1) * state.input.vazDescControl && state.gasCells[0].VGasR > 0.) {
        double precorr = unloadingPressureCorrection(state, throatFlowRate, valveIndex, (1 - laz1), -1);
        if (fabs(precorr) > 0.01 * state.initialGasPressure * state.cells[0].dt)
            precorr = (fabs(precorr) / precorr) * 0.01 * state.initialGasPressure * state.cells[0].dt;
        state.initialGasPressure += precorr;
        if (state.initialGasPressure < state.input.presMinDescG) {
            state.initialGasPressure = state.input.presMinDescG;
        }

    } else if (throatFlowRate <= (1 - laz2) * state.input.vazDescControl && state.gasSurfacePressure <= state.input.presMinDesc * 1.0000001) {
        double precorr = unloadingPressureCorrection(state, throatFlowRate, valveIndex, 1 - laz2, 1);
        if (fabs(precorr) > 0.01 * state.initialGasPressure * state.cells[0].dt)
            precorr = (fabs(precorr) / precorr) * 0.01 * state.initialGasPressure * state.cells[0].dt;
        state.initialGasPressure += precorr;
        if (state.initialGasPressure > state.input.presMaxDescG) {
            state.initialGasPressure = state.input.presMaxDescG;
        }
    }

    return velmax;
}

void solveUnloading(const GasLiftState &state) {
    for (int gasCellIndex = state.interfaceCell; gasCellIndex <= state.gasCellCount; gasCellIndex++) {
        double halfUpstreamLength = 0.5 * state.gasCells[gasCellIndex].dxL;
        double halfLocalLength = 0.5 * state.gasCells[gasCellIndex].dx0;
        double RgasR = 1.;
        if (state.gasCells[gasCellIndex].razInter <= 0.5)
            RgasR = 2 * state.gasCells[gasCellIndex].razInter;
        double RgasL = 0.;
        if (state.gasCells[gasCellIndex - 1].razInter >= 0.5)
            RgasL = 2 * (state.gasCells[gasCellIndex - 1].razInter - 0.5);
        double LGasL = halfUpstreamLength * RgasL;
        double LGasR = halfLocalLength * RgasR;
        double LLiqL = halfUpstreamLength - LGasL;
        double LLiqR = halfLocalLength - LGasR;
        double temp = state.gasCells[gasCellIndex - 1].temp;
        double pres = state.gasCells[gasCellIndex - 1].pres;
        double rhoL = state.gasCells[gasCellIndex].MasEspFlu(pres, temp);
        double viscL = state.gasCells[gasCellIndex].VisFlu(pres, temp);
        double rhoG = state.gasCells[gasCellIndex].flui.MasEspGas(pres, temp);
        double viscG = state.gasCells[gasCellIndex].flui.ViscGas(pres, temp);
        double vel1 = state.gasCells[gasCellIndex].VGasL / (rhoL * state.gasCells[gasCellIndex - 1].duto.area);
        if (state.gasCells[gasCellIndex].razInter > (*state.globals).localtiny)
            vel1 = state.gasCells[gasCellIndex].VGasL / (rhoG * state.gasCells[gasCellIndex - 1].duto.area);
        double vel2 = state.gasCells[gasCellIndex].VGasL / (rhoL * state.gasCells[gasCellIndex].duto.area);
        if (state.gasCells[gasCellIndex].razInter > (*state.globals).localtiny)
            vel2 = state.gasCells[gasCellIndex].VGasL / (rhoG * state.gasCells[gasCellIndex].duto.area);
        double re1G;
        double re1L;
        double re2G;
        double re2L;
        if (state.gasCells[gasCellIndex - 1].duto.revest == 0)
            re1L = state.gasCells[gasCellIndex - 1].Rey(state.gasCells[gasCellIndex - 1].duto.a, vel1, rhoL, viscL);
        else {
            double dhid = 4 * state.gasCells[gasCellIndex - 1].duto.area / state.gasCells[gasCellIndex - 1].duto.peri;
            re1L = state.gasCells[gasCellIndex - 1].Rey(dhid, vel1, rhoL, viscL);
        }
        if (state.gasCells[gasCellIndex].duto.revest == 0)
            re2L = state.gasCells[gasCellIndex].Rey(state.gasCells[gasCellIndex].duto.a, vel2, rhoL, viscL);
        else {
            double dhid = 4 * state.gasCells[gasCellIndex].duto.area / state.gasCells[gasCellIndex].duto.peri;
            re2L = state.gasCells[gasCellIndex].Rey(dhid, vel2, rhoL, viscL);
        }
        if (state.gasCells[gasCellIndex - 1].duto.revest == 0)
            re1G = state.gasCells[gasCellIndex - 1].Rey(state.gasCells[gasCellIndex - 1].duto.a, vel1, rhoG, viscG);
        else {
            double dhid = 4 * state.gasCells[gasCellIndex - 1].duto.area / state.gasCells[gasCellIndex - 1].duto.peri;
            re1G = state.gasCells[gasCellIndex - 1].Rey(dhid, vel1, rhoG, viscG);
        }
        if (state.gasCells[gasCellIndex].duto.revest == 0)
            re2G = state.gasCells[gasCellIndex].Rey(state.gasCells[gasCellIndex].duto.a, vel2, rhoG, viscG);
        else {
            double dhid = 4 * state.gasCells[gasCellIndex].duto.area / state.gasCells[gasCellIndex].duto.peri;
            re2G = state.gasCells[gasCellIndex].Rey(dhid, vel2, rhoG, viscG);
        }
        double upstreamLiquidFriction = state.gasCells[gasCellIndex - 1].fric(re1L, state.gasCells[gasCellIndex - 1].duto.rug / state.gasCells[gasCellIndex - 1].duto.a) * LLiqL;
        double localLiquidFriction = state.gasCells[gasCellIndex].fric(re2L, state.gasCells[gasCellIndex].duto.rug / state.gasCells[gasCellIndex].duto.a) * LLiqR;
        double hidro1L = 1 * (9.82 * sin(state.gasCells[gasCellIndex - 1].duto.teta) * rhoL) * LLiqL;
        double hidro2L = 1 * (9.82 * sin(state.gasCells[gasCellIndex].duto.teta) * rhoL) * LLiqR;
        double upstreamGasFriction = state.gasCells[gasCellIndex - 1].fric(re1G, state.gasCells[gasCellIndex - 1].duto.rug / state.gasCells[gasCellIndex - 1].duto.a) * LGasL;
        double localGasFriction = state.gasCells[gasCellIndex].fric(re2G, state.gasCells[gasCellIndex].duto.rug / state.gasCells[gasCellIndex].duto.a) * LGasR;
        double hidro1G = 1 * (9.82 * sin(state.gasCells[gasCellIndex - 1].duto.teta) * rhoG) * LGasL;
        double hidro2G = 1 * (9.82 * sin(state.gasCells[gasCellIndex].duto.teta) * rhoG) * LGasR;
        state.gasCells[gasCellIndex].pres = state.gasCells[gasCellIndex - 1].pres + (-0.5 * (upstreamLiquidFriction * rhoL + upstreamGasFriction * rhoG) * vel1 * fabs(vel1) * state.gasCells[gasCellIndex - 1].duto.peri / state.gasCells[gasCellIndex - 1].duto.area - 0.5 * (localLiquidFriction * rhoL + localGasFriction * rhoG) * vel2 * fabs(vel2) * state.gasCells[gasCellIndex].duto.peri / state.gasCells[gasCellIndex].duto.area - hidro1L - hidro2L - hidro1G - hidro2G) / 98066.52;
        state.gasCells[gasCellIndex].presL = state.gasCells[gasCellIndex - 1].pres;
        state.gasCells[gasCellIndex - 1].presR = state.gasCells[gasCellIndex].pres;

        state.temperatureUpdater.dischargeTemperature(gasCellIndex);

        state.gasCells[gasCellIndex].u1L = ((1. - state.gasCells[gasCellIndex].razInter) * state.gasCells[gasCellIndex].MasEspFlu(state.gasCells[gasCellIndex].pres, state.gasCells[gasCellIndex].temp) + state.gasCells[gasCellIndex].razInter * state.gasCells[gasCellIndex].flui.MasEspGas(pres, temp)) * state.gasCells[gasCellIndex].duto.area;
        state.gasCells[gasCellIndex - 1].u1R = state.gasCells[gasCellIndex].u1L;
        state.gasCells[gasCellIndex].u1LL = state.gasCells[gasCellIndex - 1].u1L;
    }

    double Qtotal = 0.;
    for (int gasCellIndex = state.interfaceCell; gasCellIndex <= state.gasCellCount; gasCellIndex++) {
        double temp = state.gasCells[gasCellIndex].temp;
        double pres = state.gasCells[gasCellIndex].pres;
        double rhoL = state.gasCells[gasCellIndex].MasEspFlu(pres, temp);
        double rhoG = state.gasCells[gasCellIndex].flui.MasEspGas(pres, temp);
        double qfonte = state.gasCells[gasCellIndex].massfonteCH / rhoL;
        if (state.gasCells[gasCellIndex].razInter > 0.5)
            qfonte = state.gasCells[gasCellIndex].massfonteCH / rhoG;
        Qtotal += qfonte;
    }
    double temp = state.gasCells[state.interfaceCell].temp;
    double pres = state.gasCells[state.interfaceCell].pres;
    double rhoG = state.gasCells[state.interfaceCell].flui.MasEspGas(pres, temp);
    state.gasCells[state.interfaceCell].VGasL = Qtotal * rhoG;
    state.gasCells[state.interfaceCell - 1].VGasR = state.gasCells[state.interfaceCell].VGasL;
    state.gasCells[state.interfaceCell - 2].VGasRR = state.gasCells[state.interfaceCell].VGasL;
    for (int gasCellIndex = state.interfaceCell; gasCellIndex <= state.gasCellCount; gasCellIndex++) {
        double temp = state.gasCells[gasCellIndex].temp;
        double pres = state.gasCells[gasCellIndex].pres;
        double rhoL = state.gasCells[gasCellIndex].MasEspFlu(pres, temp);
        double rhoG = state.gasCells[gasCellIndex].flui.MasEspGas(pres, temp);
        double qfonte = state.gasCells[gasCellIndex].massfonteCH / rhoL;
        if (state.gasCells[gasCellIndex].razInter > 0.5)
            qfonte = state.gasCells[gasCellIndex].massfonteCH / rhoG;
        Qtotal -= qfonte;
        if (gasCellIndex < state.gasCellCount) {
            state.gasCells[gasCellIndex + 1].VGasL = (Qtotal)*rhoL;
            state.gasCells[gasCellIndex].VGasR = state.gasCells[gasCellIndex + 1].VGasL;
        } else
            state.gasCells[gasCellIndex].VGasR = 0.;
        state.gasCells[gasCellIndex - 1].VGasRR = state.gasCells[gasCellIndex].VGasL;
    }

    state.interfaceVelocity = state.gasCells[state.interfaceCell + 1].VGasL / (state.gasCells[state.interfaceCell + 1].MasEspFlu(state.gasCells[state.interfaceCell + 1].pres, state.gasCells[state.interfaceCell + 1].temp) * state.gasCells[state.interfaceCell].duto.area);
}

void advanceInterface(const GasLiftState &state) {

    state.gasCells[state.interfaceCell].razInter = (state.gasCells[state.interfaceCell].razInterIni * state.gasCells[state.interfaceCell].dx0 + state.interfaceVelocity * state.timeStep) / state.gasCells[state.interfaceCell].dx0;
    if (state.interfaceCell == (state.gasCellCount - 1) && state.gasCells[state.interfaceCell].razInter >= 0.99) {

        double pmed = state.gasCells[state.interfaceCell].pres;
        double tmed = state.gasCells[state.interfaceCell].temp;
        double interfaceFlowArea = state.gasCells[state.interfaceCell].duto.area;
        double rho1 = state.gasCells[state.interfaceCell].flui.MasEspGas(pmed, tmed);
        double rhoL = state.gasCells[state.interfaceCell].MasEspFlu(pmed, tmed);
        state.gasCells[state.interfaceCell].VGasR = state.gasCells[state.interfaceCell].VGasR * rho1 / rhoL;
        state.gasCells[state.interfaceCell].u1L = interfaceFlowArea * rho1;
        state.gasCells[state.interfaceCell - 1].u1R = state.gasCells[state.interfaceCell].u1L;
        state.gasCells[state.interfaceCell - 1].VGasRR = state.gasCells[state.interfaceCell].VGasR;
        state.gasCells[state.interfaceCell + 1].u1LL = state.gasCells[state.interfaceCell].u1L;
        state.gasCells[state.interfaceCell + 1].VGasL = state.gasCells[state.interfaceCell].VGasR;

        state.gasCells[state.interfaceCell].razInter = 1.0;
        state.gasCells[state.interfaceCell + 1].razInter = 1.0;

        state.interfaceCell++;

        state.gasCells[state.interfaceCell].presini = state.gasCells[state.interfaceCell].pres;
        pmed = state.gasCells[state.interfaceCell].pres;
        tmed = state.gasCells[state.interfaceCell].temp;
        interfaceFlowArea = state.gasCells[state.interfaceCell].duto.area;
        rho1 = state.gasCells[state.interfaceCell].flui.MasEspGas(pmed, tmed);
        rhoL = state.gasCells[state.interfaceCell].MasEspFlu(pmed, tmed);
        state.gasCells[state.interfaceCell].VGasR = 0 * state.gasCells[state.interfaceCell].VGasR * rho1 / rhoL;
        state.gasCells[state.interfaceCell].u1L = interfaceFlowArea * rho1;
        state.gasCells[state.interfaceCell].u1R = interfaceFlowArea * rho1;
        state.gasCells[state.interfaceCell - 1].u1R = state.gasCells[state.interfaceCell].u1L;
        state.gasCells[state.interfaceCell - 1].VGasRR = state.gasCells[state.interfaceCell].VGasR;

        state.interfaceCell = 1e7;
        state.input.descarga = 0;
    }

    if (state.input.descarga == 1) {
        if (((state.gasCells[state.interfaceCell].razInter <= (*state.globals).localtiny) && (state.gasCells[state.interfaceCell].razInter >= -(*state.globals).localtiny)))
            state.gasCells[state.interfaceCell].razInter =
                0;
        else if (state.gasCells[state.interfaceCell].razInter < -(*state.globals).localtiny) {
            double dtaux;
            dtaux = -state.gasCells[state.interfaceCell].razInterIni * state.gasCells[state.interfaceCell].dx0 / state.interfaceVelocity;
            if (dtaux > (*state.globals).localtiny) {
                state.interfaceTimeStep = dtaux;
                state.gasCells[state.interfaceCell].razInter = fabs(0.);
            } else
                state.gasCells[state.interfaceCell].razInter = fabs(0.);
        } else if ((state.gasCells[state.interfaceCell].razInter >= (1. - (*state.globals).localtiny) && state.gasCells[state.interfaceCell].razInter <= (1. + (*state.globals).localtiny))) {
            state.gasCells[state.interfaceCell].razInter = 1.;
        } else if (state.gasCells[state.interfaceCell].razInter > (1. + (*state.globals).localtiny)) {
            double dtaux;
            dtaux = (1. - state.gasCells[state.interfaceCell].razInterIni) * state.gasCells[state.interfaceCell].dx0 / state.interfaceVelocity;
            if (dtaux > (*state.globals).localtiny) {
                state.interfaceTimeStep = dtaux;
                state.gasCells[state.interfaceCell].razInter = 1.;
            } else
                state.gasCells[state.interfaceCell].razInter = 1.;
        }

        if (fabs(state.gasCells[state.interfaceCell].razInter) < (*state.globals).localtiny && state.interfaceVelocity < -fabs((*state.globals).localtiny)) {
            state.gasCells[state.interfaceCell].razInter = 0.;
            if (state.interfaceCell > 0) {
                state.interfaceCell--;
                state.gasCells[state.interfaceCell].razInter = 1.;
            }
        } else if (fabs(state.gasCells[state.interfaceCell].razInter - 1.) < (*state.globals).localtiny && state.interfaceVelocity > (*state.globals).localtiny) {
            state.gasCells[state.interfaceCell].razInter = 1.;
            if (state.interfaceCell < state.gasCellCount) {
                state.interfaceCell++;
                state.gasCells[state.interfaceCell].razInter = 0.;
            }
        }
    }
}

void updateTransientGasValves(const GasLiftState &state) {
    int nvalv = state.input.nvalvgas;
    int posicM2 = state.input.master2.posic;
    for (int valveIndex = 0; valveIndex < nvalv; valveIndex++) {
        if (state.gasValveCellIndices[valveIndex] < state.interfaceCell) {
            state.gasLiftChokes[valveIndex].presEstag = state.gasCells[state.gasValveCellIndices[valveIndex]].pres;
            state.gasLiftChokes[valveIndex].presGarg = (state.cells[state.productionValveCellIndices[valveIndex]].pres - state.gasLiftChokes[valveIndex].presEstag * state.gasLiftChokes[valveIndex].frec) / (1. - state.gasLiftChokes[valveIndex].frec);
            state.gasLiftChokes[valveIndex].tempEstag = state.gasCells[state.gasValveCellIndices[valveIndex]].temp;
            if (state.cells[state.productionValveCellIndices[valveIndex]].pres < state.gasCells[state.gasValveCellIndices[valveIndex]].pres)
                state.gasCells[state.gasValveCellIndices[valveIndex]].massfonteCH =
                    state.gasLiftChokes[valveIndex].massica();
            else
                state.gasCells[state.gasValveCellIndices[valveIndex]].massfonteCH = 0.;
            if (state.gasLiftChokes[valveIndex].tipo == 1) {
                double abre = calibratedValveArea(state.gasLiftChokes[valveIndex].pcalib * 14.223595, state.gasLiftChokes[valveIndex].tcalib, (state.gasLiftChokes[valveIndex].presEstag - 1.033211) * 14.223595,
                                           (state.gasLiftChokes[valveIndex].presGarg - 1.033211) * 14.223595, state.gasLiftChokes[valveIndex].dextern, state.gasLiftChokes[valveIndex].areagarg, state.gasLiftChokes[valveIndex].areagarg / state.gasLiftChokes[valveIndex].areafole,
                                           1.8 * state.gasLiftChokes[valveIndex].tempEstag + 32);
                state.gasCells[state.gasValveCellIndices[valveIndex]].massfonteCH *= abre;
            }
        }
    }
    for (int gasCellIndex = state.interfaceCell; gasCellIndex <= state.gasCellCount; gasCellIndex++) {

        state.gasCells[gasCellIndex].massfonteCH = 0.;
        if (state.gasCells[gasCellIndex].razInter < 0.5) {
            for (int candidateValveIndex = 0; candidateValveIndex < nvalv; candidateValveIndex++) {
                if (gasCellIndex == state.gasValveCellIndices[candidateValveIndex]) {
                    state.gasLiftChokes[candidateValveIndex].presEstag = state.gasCells[state.gasValveCellIndices[candidateValveIndex]].pres;
                    state.gasLiftChokes[candidateValveIndex].presGarg = (state.cells[state.productionValveCellIndices[candidateValveIndex]].pres - state.gasLiftChokes[candidateValveIndex].presEstag * state.gasLiftChokes[candidateValveIndex].frecliq) / (1. - state.gasLiftChokes[candidateValveIndex].frecliq);
                    state.gasLiftChokes[candidateValveIndex].tempEstag = state.gasCells[state.gasValveCellIndices[candidateValveIndex]].temp;
                    if (state.cells[state.productionValveCellIndices[candidateValveIndex]].pres < state.gasCells[state.gasValveCellIndices[candidateValveIndex]].pres)
                        state.gasCells[state.gasValveCellIndices[candidateValveIndex]].massfonteCH =
                            state.gasLiftChokes[candidateValveIndex].massica(1, state.input.salinDescarga);
                    else
                        state.gasCells[state.gasValveCellIndices[candidateValveIndex]].massfonteCH = 0.;
                    if (state.gasLiftChokes[candidateValveIndex].tipo == 1) {
                        double abre = calibratedValveArea(state.gasLiftChokes[candidateValveIndex].pcalib * 14.223595, state.gasLiftChokes[candidateValveIndex].tcalib,
                                                   (state.gasLiftChokes[candidateValveIndex].presEstag - 1.033211) * 14.223595, (state.gasLiftChokes[candidateValveIndex].presGarg - 1.033211) * 14.223595, state.gasLiftChokes[candidateValveIndex].dextern,
                                                   state.gasLiftChokes[candidateValveIndex].areagarg, state.gasLiftChokes[candidateValveIndex].areagarg / state.gasLiftChokes[candidateValveIndex].areafole, 1.8 * state.gasLiftChokes[candidateValveIndex].tempEstag + 32);
                        state.gasCells[state.gasValveCellIndices[candidateValveIndex]].massfonteCH *= abre;
                    }
                }
            }
        } else {
            for (int candidateValveIndex = 0; candidateValveIndex < nvalv; candidateValveIndex++) {
                if (gasCellIndex == state.gasValveCellIndices[candidateValveIndex]) {
                    state.gasLiftChokes[candidateValveIndex].presEstag = state.gasCells[state.gasValveCellIndices[candidateValveIndex]].pres;
                    state.gasLiftChokes[candidateValveIndex].presGarg = (state.cells[state.productionValveCellIndices[candidateValveIndex]].pres - state.gasLiftChokes[candidateValveIndex].presEstag * state.gasLiftChokes[candidateValveIndex].frec) / (1. - state.gasLiftChokes[candidateValveIndex].frec);
                    state.gasLiftChokes[candidateValveIndex].tempEstag = state.gasCells[state.gasValveCellIndices[candidateValveIndex]].temp;
                    if (state.cells[state.productionValveCellIndices[candidateValveIndex]].pres < state.gasCells[state.gasValveCellIndices[candidateValveIndex]].pres)
                        state.gasCells[state.gasValveCellIndices[candidateValveIndex]].massfonteCH =
                            state.gasLiftChokes[candidateValveIndex].massica();
                    else
                        state.gasCells[state.gasValveCellIndices[candidateValveIndex]].massfonteCH = 0.;
                    if (state.gasLiftChokes[candidateValveIndex].tipo == 1) {
                        double abre = calibratedValveArea(state.gasLiftChokes[candidateValveIndex].pcalib * 14.223595, state.gasLiftChokes[candidateValveIndex].tcalib,
                                                   (state.gasLiftChokes[candidateValveIndex].presEstag - 1.033211) * 14.223595, (state.gasLiftChokes[candidateValveIndex].presGarg - 1.033211) * 14.223595, state.gasLiftChokes[candidateValveIndex].dextern,
                                                   state.gasLiftChokes[candidateValveIndex].areagarg, state.gasLiftChokes[candidateValveIndex].areagarg / state.gasLiftChokes[candidateValveIndex].areafole, 1.8 * state.gasLiftChokes[candidateValveIndex].tempEstag + 32);
                        state.gasCells[state.gasValveCellIndices[candidateValveIndex]].massfonteCH *= abre;
                    }
                }
            }
        }
    }

    double areamenor = state.gasCells[posicM2].duto.area;
    if (areamenor > state.gasCells[posicM2].dutoR.area)
        areamenor = state.gasCells[posicM2].dutoR.area;
    if (state.gasCells[posicM2].chkcell.areagarg <= 0.01 * areamenor &&
        state.gasCells[posicM2].chkcell.areagarg >= 1e-5 * areamenor) {
        state.gasCells[posicM2].chkcell.tempEstag = state.gasCells[posicM2].temp;
        double rhoM = state.gasCells[posicM2].flui.MasEspGas(state.gasCells[posicM2].pres, state.gasCells[posicM2].temp);
        double rhoJ = state.gasCells[posicM2 + 1].flui.MasEspGas(state.gasCells[posicM2 + 1].pres, state.gasCells[posicM2 + 1].temp);
        double hidroM = -0.5 * rhoM * state.gasCells[posicM2].dx0 * sin(state.gasCells[posicM2].duto.teta) / 98066.52;
        double hidroJ = 0.5 * rhoJ * state.gasCells[posicM2 + 1].dx0 * sin(state.gasCells[posicM2 + 1].duto.teta) / 98066.52;
        state.gasCells[posicM2].chkcell.presEstag = state.gasCells[posicM2].pres + hidroM;
        state.gasCells[posicM2].chkcell.presGarg = state.gasCells[posicM2 + 1].pres + hidroJ;
        if (state.gasCells[posicM2 + 1].pres < state.gasCells[posicM2].pres)
            state.gasCells[posicM2].fonteM2 = -state.gasCells[posicM2].chkcell.massica();
        else {
            state.gasCells[posicM2].chkcell.presEstag = state.gasCells[posicM2 + 1].pres;
            state.gasCells[posicM2].chkcell.presGarg = state.gasCells[posicM2].pres;
            state.gasCells[posicM2].chkcell.tempEstag = state.gasCells[posicM2 + 1].temp;
            state.gasCells[posicM2].fonteM2 = state.gasCells[posicM2].chkcell.massica();
        }
        state.gasCells[posicM2 + 1].fonteM2 = -state.gasCells[posicM2].fonteM2;
    } else {
        state.gasCells[posicM2].fonteM2 = 0.;

        state.gasCells[posicM2 + 1].fonteM2 = 0.;
    }
}

double searchUnloadingInjectionPressure(const GasLiftState &state) {
    double vazmax = 0;
    double vazmaxinst = 0;

    double laz1 = 0.1;
    double laz2 = 0.4;
    if (state.input.descarga == 1) {
        double presDescini = state.cells[state.lastCell].pres;
        if ((*state.globals).lixo5 > state.input.tempoLatenciaDesc) {
            int nvalv = state.input.nvalvgas;
            int ivalv = 0;
            for (int valveIndex = 0; valveIndex < nvalv; valveIndex++) {
                int gasValveCell = state.gasValveCellIndices[valveIndex];
                int productionValveCell = state.productionValveCellIndices[valveIndex];
                double pmed = state.cells[productionValveCell].pres;
                if (gasValveCell > state.interfaceCell) {
                    double rho1 = state.gasCells[gasValveCell].MasEspFlu(state.gasCells[gasValveCell].pres, state.gasCells[gasValveCell].temp);
                    state.gasLiftChokes[valveIndex].presEstag = state.gasCells[gasValveCell].pres;
                    state.gasLiftChokes[valveIndex].presGarg = (pmed - state.gasLiftChokes[valveIndex].presEstag * state.gasLiftChokes[valveIndex].frec) /
                                           (1. - state.gasLiftChokes[valveIndex].frec);
                    state.gasLiftChokes[valveIndex].tempEstag = state.gasCells[gasValveCell].temp;
                    double massica;
                    if (pmed < state.gasCells[gasValveCell].pres)
                        massica = state.gasLiftChokes[valveIndex].massica(1, state.input.salinDescarga);
                    else
                        massica = 0.;
                    if (state.gasLiftChokes[valveIndex].tipo == 1) {

                        double abre = calibratedValveArea(state.gasLiftChokes[valveIndex].pcalib * 14.223595, state.gasLiftChokes[valveIndex].tcalib, (state.gasLiftChokes[valveIndex].presEstag - 1.033211) * 14.223595,
                                                   (state.gasLiftChokes[valveIndex].presGarg - 1.033211) * 14.223595, state.gasLiftChokes[valveIndex].dextern, state.gasLiftChokes[valveIndex].areagarg, state.gasLiftChokes[valveIndex].areagarg / state.gasLiftChokes[valveIndex].areafole,
                                                   1.8 * state.gasLiftChokes[valveIndex].tempEstag + 32);
                        massica *= abre;
                    }
                    double vazmaxAux = massica / (rho1);
                    if (vazmaxAux > vazmax) {
                        vazmax = vazmaxAux;
                        ivalv = valveIndex;
                    }
                }
            }
            if (state.unloadingTimeSteps.size() > state.maximumContinuousUnloadingCount || state.meanUnloadingTemperature > state.continuousMeanUnloadingTemperature) {
                state.meanUnloadingTemperature -= state.unloadingTimeSteps[0];
                state.unloadingTimeSteps.erase(state.unloadingTimeSteps.begin());
                state.meanUnloadingFlowRate -= state.maximumMeanUnloadingFlowRates[0];
                state.maximumMeanUnloadingFlowRates.erase(state.maximumMeanUnloadingFlowRates.begin());
            }

            state.maximumMeanUnloadingFlowRates.push_back(vazmax * state.timeStep);
            state.meanUnloadingFlowRate += vazmax * state.timeStep;
            state.unloadingTimeSteps.push_back(state.timeStep);
            state.meanUnloadingTemperature += state.timeStep;

            double pondera = 0.5;
            vazmaxinst = vazmax;
            vazmax = pondera * vazmax + (1. - pondera) * state.meanUnloadingFlowRate / state.meanUnloadingTemperature;
            computeUnloadingValvePressure(state, vazmax, ivalv);
            if (state.gasValveCellIndices[ivalv] >= state.interfaceCell) {
                if (state.cells[state.lastCell - 1].MC > 0.0) {
                    if (vazmax > (1. - laz1) * state.input.vazDescControl) {

                        double precorr = 0.;
                        precorr = unloadingPressureCorrection(state, vazmax, ivalv, 1. - laz1, -1);
                        if (fabs(precorr) > 0.01 * state.gasSurfacePressure * state.cells[0].dt)
                            precorr = (fabs(precorr) / precorr) * 0.01 * state.gasSurfacePressure * state.cells[0].dt;
                        state.gasSurfacePressure -= precorr;
                        if (state.gasSurfacePressure > state.input.presMaxDesc)
                            state.gasSurfacePressure = state.input.presMaxDesc;
                    } else if (vazmax <= (1 - laz2) * state.input.vazDescControl) {

                        double precorr = 0.;
                        precorr = unloadingPressureCorrection(state, vazmax, ivalv, 1. - laz2, 1);
                        if (fabs(precorr) > 0.01 * state.gasSurfacePressure * state.cells[0].dt)
                            precorr = (fabs(precorr) / precorr) * 0.01 * state.gasSurfacePressure * state.cells[0].dt;
                        state.gasSurfacePressure -= precorr;
                        if (state.gasSurfacePressure < state.input.presMinDesc)
                            state.gasSurfacePressure = state.input.presMinDesc;
                    }
                } else {
                    state.gasSurfacePressure *= (1 - 0.01 * state.timeStep);
                    if (state.gasSurfacePressure < state.input.presMinDesc)
                        state.gasSurfacePressure = state.input.presMinDesc;
                }
            }
        } else {
            if ((*state.globals).lixo5 < 0.5 * state.input.tempoLatenciaDesc)
                state.gasSurfacePressure = presDescini - (presDescini - state.input.presMinDesc) * state.cells[0].dt / (0.5 * state.input.tempoLatenciaDesc - (*state.globals).lixo5);
            else
                state.gasSurfacePressure = state.input.presMinDesc;
            state.initialGasPressure = state.input.presIniDescG;
        }
    } else {
        state.gasSurfacePressure *= 0.95;
        if (state.gasSurfacePressure < state.input.presMinDesc)
            state.gasSurfacePressure = state.input.presMinDesc;
        if (state.initialGasPressure > state.input.presIniDesc) {
            state.initialGasPressure *= 0.95;
            if (state.initialGasPressure < state.input.presIniDesc)
                state.initialGasPressure = state.input.presIniDesc;
        } else if (state.initialGasPressure < state.input.presIniDesc) {
            state.initialGasPressure *= 1.05;
            if (state.initialGasPressure > state.input.presIniDesc)
                state.initialGasPressure = state.input.presIniDesc;
        }
    }
    return vazmaxinst;
}

void advanceGasSubStep(const GasLiftState &state) {

    for (int gasCellIndex = 0; gasCellIndex <= state.gasCellCount; gasCellIndex++)
        state.gasCells[gasCellIndex].DeVoltaParaoFuturo();
    state.interfaceTimeStep = state.timeStep;
    if (state.input.descarga == 1) {
        state.initialInterfaceCell = state.interfaceCell;
        state.initialInterfaceTimeStep = state.interfaceTimeStep;
        state.initialInterfaceVelocity = state.interfaceVelocity;
        advanceInterface(state);
    }
    if (state.interfaceTimeStep < state.timeStep) {
        state.timeStep = state.interfaceTimeStep;
        for (int cellIndex = 0; cellIndex <= state.lastCell; cellIndex++) {
            state.cells[cellIndex].dt = state.timeStep;
            state.cells[cellIndex].dt2 = state.timeStep;
            state.cells[cellIndex].dtPig = state.timeStep;
        }
    }

    for (int gasCellIndex = 0; gasCellIndex <= state.gasCellCount; gasCellIndex++) {
        state.gasCells[gasCellIndex].dt = state.timeStep;
    }

    int cellLimit;
    if (state.interfaceCell <= state.gasCellCount)
        cellLimit = state.interfaceCell - 1;
    else
        cellLimit = state.gasCellCount + 1;
    for (int gasCellIndex = 0; gasCellIndex < state.gasCellCount + 1; gasCellIndex++)
        state.gasCells[gasCellIndex].dTdt = 0;

    int ciclomax = state.input.cicloAcopTerm;
    for (int ciclo = 0; ciclo <= ciclomax; ciclo++) {
        double abertoChk = 1.;
        if (state.gasCells[0].tipoCC == 0) {
            abertoChk = state.injectionChoke.areagarg / state.gasCells[0].duto.area;
            if (abertoChk < 0.2) {
                state.injectionChoke.presEstag = state.initialGasPressure;
                state.injectionChoke.tempEstag = state.initialGasTemperature;
                state.injectionChoke.presGarg = state.gasCells[0].pres;
                state.gasCells[0].massfonteCH = state.injectionChoke.massica();
            }
        }
#pragma omp parallel for num_threads((*state.globals).ntrd)
        for (int gasCellIndex = 0; gasCellIndex <= state.gasCellCount; gasCellIndex++) {

            state.gasCells[gasCellIndex].GeraLocal(state.gasCellCount, state.initialGasPressure, state.initialGasTemperature, abertoChk);
            for (int bandIndex = 0; bandIndex < 9; bandIndex++) {
                state.gasSystemMatrix[3 * gasCellIndex][bandIndex - 3] = state.gasCells[gasCellIndex].local[0][bandIndex];
                state.gasSystemMatrix[3 * gasCellIndex + 1][bandIndex - 4] = state.gasCells[gasCellIndex].local[1][bandIndex];
                state.gasSystemMatrix[3 * gasCellIndex + 2][bandIndex - 5] = state.gasCells[gasCellIndex].local[2][bandIndex];
            }
            state.gasFreeTerms[3 * gasCellIndex] = state.gasCells[gasCellIndex].TL[0];
            state.gasFreeTerms[3 * gasCellIndex + 1] = state.gasCells[gasCellIndex].TL[1];
            state.gasFreeTerms[3 * gasCellIndex + 2] = state.gasCells[gasCellIndex].TL[2];
        }
        state.gasSystemMatrix.GaussElimPP(state.gasFreeTerms);
        updateGasLine(state);
        state.gasCells[0].temp = state.initialGasTemperature;
        state.gasCells[0].dTdt = (state.gasCells[0].temp - state.gasCells[0].tempini) / state.timeStep;
        if (state.gasCells[state.gasCellCount].VGasR < 0.)
            state.gasCells[state.gasCellCount].temp = 20.;
#pragma omp parallel for num_threads((*state.globals).ntrd)
        for (int gasCellIndex = 1; gasCellIndex < cellLimit; gasCellIndex++) {
            state.temperatureUpdater.gasTemperature(gasCellIndex, state.gasCells[gasCellIndex - 1].tempini);
            state.gasCells[gasCellIndex].dTdt = (state.gasCells[gasCellIndex].temp - state.gasCells[gasCellIndex].tempini) / state.timeStep;
        }
        if (ciclo < ciclomax)
            for (int renewedCellIndex = 0; renewedCellIndex <= state.gasCellCount; renewedCellIndex++)
                state.gasCells[renewedCellIndex].FeiticoDoTempo();
    }

#pragma omp parallel for num_threads((*state.globals).ntrd)
    for (int gasCellIndex = 0; gasCellIndex < cellLimit; gasCellIndex++)
        state.gasCells[gasCellIndex].rg = state.gasCells[gasCellIndex].flui.MasEspGas(state.gasCells[gasCellIndex].pres, state.gasCells[gasCellIndex].temp);

    for (int gasCellIndex = 0; gasCellIndex < cellLimit; gasCellIndex++) {
        if (gasCellIndex > 0)
            state.gasCells[gasCellIndex - 1].rgR = state.gasCells[gasCellIndex].rg;
        state.gasCells[gasCellIndex].u1L = state.gasCells[gasCellIndex].rg * state.gasCells[gasCellIndex].duto.area;
        if (gasCellIndex == 0)
            state.gasCells[gasCellIndex].u1LL = state.gasCells[gasCellIndex].u1L;
        else
            state.gasCells[gasCellIndex].u1LL = state.gasCells[gasCellIndex - 1].u1L;
        if (gasCellIndex > 0)
            state.gasCells[gasCellIndex - 1].u1R = state.gasCells[gasCellIndex].u1L;
        if (gasCellIndex == state.gasCellCount) {
            state.gasCells[gasCellIndex].u1R = state.gasCells[gasCellIndex].u1L;
            state.gasCells[state.gasCellCount].rgR = state.gasCells[state.gasCellCount].rg;
        }
    }

    for (int gasCellIndex = 0; gasCellIndex <= state.gasCellCount; gasCellIndex++)
        state.gasCells[gasCellIndex].presini = state.gasCells[gasCellIndex].pres;

    if (state.input.descarga == 1)
        solveUnloading(state);
}

void advanceBufferedGasSubStep(const GasLiftState &state) {

    for (int gasCellIndex = 0; gasCellIndex < state.gasCellCount + 1; gasCellIndex++)
        state.gasCells[gasCellIndex].dTdt = 0;
    double abertoChk = 1.;
    if (state.gasCells[0].tipoCC == 0) {
        abertoChk = state.injectionChoke.areagarg / state.gasCells[0].duto.area;
        if (abertoChk < 0.2) {
            state.injectionChoke.presEstag = state.initialGasPressure;
            state.injectionChoke.tempEstag = state.initialGasTemperature;
            state.injectionChoke.presGarg = state.gasCells[0].pres;
            state.gasCells[0].massfonteCH = state.injectionChoke.massica();
        }
    }
#pragma omp parallel for num_threads((*state.globals).ntrd)
    for (int gasCellIndex = 0; gasCellIndex <= state.gasCellCount; gasCellIndex++) {

        state.gasCells[gasCellIndex].GeraLocal(state.gasCellCount, state.initialGasPressure, state.initialGasTemperature, abertoChk);
        for (int bandIndex = 0; bandIndex < 9; bandIndex++) {
            state.gasSystemMatrix[3 * gasCellIndex][bandIndex - 3] = state.gasCells[gasCellIndex].local[0][bandIndex];
            state.gasSystemMatrix[3 * gasCellIndex + 1][bandIndex - 4] = state.gasCells[gasCellIndex].local[1][bandIndex];
            state.gasSystemMatrix[3 * gasCellIndex + 2][bandIndex - 5] = state.gasCells[gasCellIndex].local[2][bandIndex];
        }
        state.gasFreeTerms[3 * gasCellIndex] = state.gasCells[gasCellIndex].TL[0];
        state.gasFreeTerms[3 * gasCellIndex + 1] = state.gasCells[gasCellIndex].TL[1];
        state.gasFreeTerms[3 * gasCellIndex + 2] = state.gasCells[gasCellIndex].TL[2];
    }
    state.gasSystemMatrix.GaussElimPP(state.gasFreeTerms);
    updateBufferedGasLine(state);
}

void connectTubing(const GasLiftState &state) {
    for (int annulusCellIndex = state.annulusTubingStart; annulusCellIndex >= state.annulusTubingEnd; annulusCellIndex--) {
        int tubingCellIndex = state.annulusTubingStart + state.tubingAnnulusStart - annulusCellIndex;
        state.cells[annulusCellIndex].calor.Textern2 = state.gasCells[tubingCellIndex].calor.Tcamada[0][0];
        state.cells[annulusCellIndex].calor.betext = state.gasCells[tubingCellIndex].calor.betint;
        int icam = state.cells[annulusCellIndex].calor.geom.ncamadas - 1;
        int idisc = state.cells[annulusCellIndex].calor.ncamada[icam] - 1;
        state.gasCells[tubingCellIndex].calor.Tint2 = state.cells[annulusCellIndex].calor.Tcamada[icam][idisc];
        state.cells[annulusCellIndex].calor.colunaDia = state.gasCells[tubingCellIndex].duto.dia;
        state.cells[annulusCellIndex].calor.geom.b = state.gasCells[tubingCellIndex].calor.geom.a;

        state.cells[annulusCellIndex].calor.Textern1 = state.gasCells[tubingCellIndex].temp;
        state.cells[annulusCellIndex].calor.Vextern1 = state.gasCells[tubingCellIndex].VGasR / state.gasCells[tubingCellIndex].u1L;
        if (tubingCellIndex < state.interfaceCell) {
            double gasSpecificHeat = state.gasCells[tubingCellIndex].flui.CalorGas(state.gasCells[tubingCellIndex].pres, state.gasCells[tubingCellIndex].temp);
            double rhog = state.gasCells[tubingCellIndex].flui.MasEspGas(state.gasCells[tubingCellIndex].pres, state.gasCells[tubingCellIndex].temp);
            state.cells[annulusCellIndex].calor.kextern1 = state.gasCells[tubingCellIndex].flui.CondGas(state.gasCells[tubingCellIndex].pres, state.gasCells[tubingCellIndex].temp);
            state.cells[annulusCellIndex].calor.cpextern1 = gasSpecificHeat;
            state.cells[annulusCellIndex].calor.rhoextern1 = rhog;
            state.cells[annulusCellIndex].calor.viscextern1 = state.gasCells[tubingCellIndex].flui.ViscGas(state.gasCells[tubingCellIndex].pres, state.gasCells[tubingCellIndex].temp) * 1.e-3;
        } else {
            double gasSpecificHeat = state.gasCells[tubingCellIndex].CalorLiq(state.gasCells[tubingCellIndex].pres, state.gasCells[tubingCellIndex].temp);
            double rhog = state.gasCells[tubingCellIndex].MasEspFlu(state.gasCells[tubingCellIndex].pres, state.gasCells[tubingCellIndex].temp);
            state.cells[annulusCellIndex].calor.kextern1 = state.gasCells[tubingCellIndex].CondLiq(state.gasCells[tubingCellIndex].pres, state.gasCells[tubingCellIndex].temp);
            state.cells[annulusCellIndex].calor.cpextern1 = gasSpecificHeat;
            state.cells[annulusCellIndex].calor.rhoextern1 = rhog;
            state.cells[annulusCellIndex].calor.viscextern1 = state.gasCells[tubingCellIndex].VisFlu(state.gasCells[tubingCellIndex].pres, state.gasCells[tubingCellIndex].temp) * 1.e-3;
        }
    }
}

void solveGasLine(const GasLiftState &state) {
    if (state.input.lingas > 0) {
        double gasTemperature;
        // dt,celula,ColunaAnulaIni,ColunaAnulaFim,
        updateTransientGasValves(state);
        for (int valveIndex = 0; valveIndex < state.input.nvalvgas; valveIndex++) {
            int posGLP = state.productionValveCellIndices[valveIndex];
            int posGLG = state.gasValveCellIndices[valveIndex];
            if (posGLG < state.interfaceCell || (posGLG == state.interfaceCell && state.gasCells[posGLG].razInter > 0.5)) {
                state.cells[posGLP].acsr.injg.QGas = state.gasCells[posGLG].massfonteCH * 86400. / state.gasCells[posGLG].flui.MasEspGas(1., 15.);
                state.cells[posGLP].acsr.injg.tipoflu = 0;
            } else {
                state.cells[posGLP].acsr.injg.QGas = state.gasCells[posGLG].massfonteCH;
                state.cells[posGLP].acsr.injg.tipoflu = 1;
            }
            if (posGLG < state.interfaceCell || (posGLG == state.interfaceCell && state.gasCells[posGLG].razInter > 0.5)) {
                if (state.gasLiftChokes[valveIndex].presEstag > state.gasLiftChokes[valveIndex].presGarg) {
                    gasTemperature = state.temperatureUpdater.gasLiftDischargeTemperature(valveIndex);
                } else
                    gasTemperature = state.cells[posGLP].temp;
                state.cells[posGLP].acsr.injg.temp = gasTemperature;
                if (state.cells[posGLP].acsr.injg.temp < -50)
                    state.cells[posGLP].acsr.injg.temp = -50;
            } else
                state.cells[posGLP].acsr.injg.temp = state.gasLiftChokes[valveIndex].tempEstag;
            state.gasCells[posGLG].pEstag = state.gasLiftChokes[valveIndex].presEstag;
            state.gasCells[posGLG].tEstag = state.gasLiftChokes[valveIndex].tempEstag;
            state.gasCells[posGLG].pGarg = state.gasLiftChokes[valveIndex].presGarg;
            state.gasCells[posGLG].tGarg = state.gasLiftChokes[valveIndex].tempGarg;
            state.gasCells[posGLG].qGarg = state.gasLiftChokes[valveIndex].qGarg;
            state.gasCells[posGLG].areaGarg = state.gasLiftChokes[valveIndex].areagarg;
        }
        advanceGasSubStep(state);
    }
}

void connectTubingSteady(const GasLiftState &state) {
    for (int annulusCellIndex = state.annulusTubingStart; annulusCellIndex >= state.annulusTubingEnd; annulusCellIndex--) {
        int tubingCellIndex = state.annulusTubingStart + state.tubingAnnulusStart - annulusCellIndex;

        state.cells[annulusCellIndex].calor.Textern2 = state.gasCells[tubingCellIndex].calor.Textern1;
        state.cells[annulusCellIndex].calor.colunaDia = state.gasCells[tubingCellIndex].duto.dia;

        double gasSpecificHeat = state.gasCells[tubingCellIndex].flui.CalorGas(state.gasCells[tubingCellIndex].pres, state.gasCells[tubingCellIndex].temp);
        double rhog = state.gasCells[tubingCellIndex].flui.MasEspGas(state.gasCells[tubingCellIndex].pres, state.gasCells[tubingCellIndex].temp);
        state.cells[annulusCellIndex].calor.Textern1 = state.gasCells[tubingCellIndex].calor.Textern1;
        state.cells[annulusCellIndex].calor.Vextern1 = state.gasCells[tubingCellIndex].VGasR / state.gasCells[tubingCellIndex].u1L;
        state.cells[annulusCellIndex].calor.kextern1 = state.gasCells[tubingCellIndex].flui.CondGas(state.gasCells[tubingCellIndex].pres, state.gasCells[tubingCellIndex].temp);
        state.cells[annulusCellIndex].calor.cpextern1 = gasSpecificHeat;
        state.cells[annulusCellIndex].calor.rhoextern1 = rhog;
        state.cells[annulusCellIndex].calor.viscextern1 = state.gasCells[tubingCellIndex].flui.ViscGas(state.gasCells[tubingCellIndex].pres, state.gasCells[tubingCellIndex].temp) * 1.e-3;
    }
}

void initializeTubingConnectionSteady(const GasLiftState &state) {
    for (int annulusCellIndex = state.annulusTubingStart; annulusCellIndex >= state.annulusTubingEnd; annulusCellIndex--) {
        int tubingCellIndex = state.annulusTubingStart + state.tubingAnnulusStart - annulusCellIndex;
        state.cells[annulusCellIndex].calor.Textern2 = state.gasCells[tubingCellIndex].calor.Textern2;
        state.cells[annulusCellIndex].calor.colunaDia = state.gasCells[tubingCellIndex].duto.dia;

        double gasSpecificHeat = state.gasCells[tubingCellIndex].flui.CalorGas(state.cells[annulusCellIndex].pres, state.input.celg[tubingCellIndex].textern);
        double rhog = state.gasCells[tubingCellIndex].flui.MasEspGas(state.cells[annulusCellIndex].pres, state.input.celg[tubingCellIndex].textern);
        state.cells[annulusCellIndex].calor.Textern1 = state.gasCells[tubingCellIndex].calor.Textern1;
        state.cells[annulusCellIndex].calor.Vextern1 = 1.;
        state.cells[annulusCellIndex].calor.kextern1 = state.gasCells[tubingCellIndex].flui.CondGas(state.cells[annulusCellIndex].pres, state.input.celg[tubingCellIndex].textern);
        state.cells[annulusCellIndex].calor.cpextern1 = gasSpecificHeat;
        state.cells[annulusCellIndex].calor.rhoextern1 = rhog;
        state.cells[annulusCellIndex].calor.viscextern1 = 0.16 * 1.e-3;
    }
}

double steadyGasPressureDrop(const GasLiftState &state, int cellIndex) {
    double dx = 0.5 * state.gasCells[cellIndex].dx0;
    double diameter = state.gasCells[cellIndex].duto.a;
    double area = 0.25 * M_PI * diameter * diameter;
    double perimeter = state.gasCells[cellIndex].duto.peri;
    double rhog = state.gasCells[cellIndex].flui.MasEspGas(state.gasCells[cellIndex].pres, state.gasCells[cellIndex].temp);

    double VGasmedR;
    VGasmedR = state.gasCells[cellIndex - 1].VGasR;

    double vel1 = VGasmedR / state.gasCells[cellIndex].u1L;

    double visc = state.gasCells[cellIndex].flui.ViscGas(state.gasCells[cellIndex].pres, state.gasCells[cellIndex].temp);

    double reynolds;
    if (state.gasCells[cellIndex].duto.revest == 0)
        reynolds = state.gasCells[cellIndex].Rey(state.gasCells[cellIndex].duto.a, vel1, rhog, visc);
    else {
        double dhid = 4 * area / perimeter;
        reynolds = state.gasCells[cellIndex].Rey(dhid, vel1, rhog, visc);
    }
    double frictionFactor = state.gasCells[cellIndex].fric(reynolds, state.gasCells[cellIndex].duto.rug / diameter);
    double gradfric = state.gasCells[cellIndex].dPdLFric * (0.5 * frictionFactor * rhog * (fabs(vel1) * vel1) * perimeter * dx / area);
    double gradhidro = state.gasCells[cellIndex].dPdLHidro * (9.82 * sin(state.gasCells[cellIndex].duto.teta) * rhog * dx);

    double difpres = (gradfric + gradhidro) / 98066.5;
    double pmed = state.gasCells[cellIndex].pres + difpres;

    double tmed;
    tmed = (state.gasCells[cellIndex].dx0 * state.gasCells[cellIndex].temp + state.gasCells[cellIndex].dxL * state.gasCells[cellIndex].tempL) / (state.gasCells[cellIndex].dx0 + state.gasCells[cellIndex].dxL);

    dx = 0.5 * state.gasCells[cellIndex].dxL;
    diameter = state.gasCells[cellIndex].dutoL.a;
    area = 0.25 * M_PI * diameter * diameter;
    perimeter = state.gasCells[cellIndex].dutoL.peri;
    rhog = state.gasCells[cellIndex].flui.MasEspGas(pmed, tmed);

    VGasmedR = state.gasCells[cellIndex - 1].VGasR;
    vel1 = VGasmedR / (rhog * area);

    if (cellIndex > 0)
        visc = state.gasCells[cellIndex - 1].flui.ViscGas(pmed, tmed);
    else
        visc = state.gasCells[cellIndex].flui.ViscGas(pmed, tmed);

    if (state.gasCells[cellIndex].dutoL.revest == 0) {
        if (cellIndex > 0)
            reynolds = state.gasCells[cellIndex].Rey(state.gasCells[cellIndex].dutoL.a, vel1, rhog, visc);
        else
            reynolds = state.gasCells[cellIndex - 1].Rey(state.gasCells[cellIndex].dutoL.a, vel1, rhog, visc);
    } else {
        double dhid = 4 * area / perimeter;
        if (cellIndex > 0)
            reynolds = state.gasCells[cellIndex - 1].Rey(dhid, vel1, rhog, visc);
        else
            reynolds = state.gasCells[cellIndex].Rey(dhid, vel1, rhog, visc);
    }
    double dpFLoc;
    double dpHLoc;
    if (cellIndex > 0) {
        frictionFactor = state.gasCells[cellIndex - 1].fric(reynolds, state.gasCells[cellIndex].dutoL.rug / diameter);
        dpFLoc = state.gasCells[cellIndex - 1].dPdLFric;
        dpHLoc = state.gasCells[cellIndex - 1].dPdLHidro;
    } else {
        frictionFactor = state.gasCells[cellIndex].fric(reynolds, state.gasCells[cellIndex].dutoL.rug / diameter);
        dpFLoc = state.gasCells[cellIndex].dPdLFric;
        dpHLoc = state.gasCells[cellIndex].dPdLHidro;
    }
    gradfric = dpFLoc * (0.5 * frictionFactor * rhog * (fabs(vel1) * vel1) * perimeter * dx / area);
    gradhidro = dpHLoc * (9.82 * sin(state.gasCells[cellIndex].dutoL.teta) * rhog * dx);
    difpres += (gradfric + gradhidro) / 98066.5;
    return difpres;
}

double steadyInjectionPressureDrop(const GasLiftState &state, int cellIndex) {
    double dx = 0.5 * state.cells[cellIndex].dx;
    double diameter = state.cells[cellIndex].duto.a;
    double area = 0.25 * M_PI * diameter * diameter;
    double perimeter = state.cells[cellIndex].duto.peri;
    double tmed;
    tmed = (state.cells[cellIndex].dx * state.cells[cellIndex].temp + state.cells[cellIndex].dxL * state.cells[cellIndex - 1].temp) / (state.cells[cellIndex].dx + state.cells[cellIndex].dxL);
    double completionFluidDensity = state.cells[cellIndex].fluicol.MasEspFlu(state.cells[cellIndex].presaux, tmed);

    double vel1 = state.cells[cellIndex - 1].QL / (area);

    double visc = state.cells[cellIndex].fluicol.VisFlu(state.cells[cellIndex].presaux, tmed);

    double reynolds;
    if (state.cells[cellIndex].duto.revest == 0)
        reynolds = state.cells[cellIndex].Rey(state.cells[cellIndex].duto.a, vel1, completionFluidDensity, visc);
    else {
        double dhid = 4 * area / perimeter;
        reynolds = state.cells[cellIndex].Rey(dhid, vel1, completionFluidDensity, visc);
    }
    double frictionFactor = state.cells[cellIndex].fric(reynolds, state.cells[cellIndex].duto.rug / diameter);
    double gradfric = state.cells[cellIndex].dPdLFric * 0.5 * frictionFactor * completionFluidDensity * (fabs(vel1) * vel1) * perimeter * dx / area;
    double gradhidro = state.cells[cellIndex].dPdLHidro * 9.82 * sin(state.cells[cellIndex].duto.teta) * completionFluidDensity * dx;

    double difpres = (gradfric + gradhidro) / 98066.5;
    double pmed = state.cells[cellIndex - 1].pres;

    dx = 0.5 * state.cells[cellIndex].dxL;
    diameter = state.cells[cellIndex].dutoL.a;
    area = 0.25 * M_PI * diameter * diameter;
    perimeter = state.cells[cellIndex].dutoL.peri;
    tmed = state.cells[cellIndex - 1].temp;
    completionFluidDensity = state.cells[cellIndex].fluicol.MasEspFlu(pmed, tmed);

    if (cellIndex > 1)
        vel1 = state.cells[cellIndex - 2].QL / (area);

    visc = state.cells[cellIndex].fluicol.VisFlu(pmed, tmed);

    if (state.cells[cellIndex].dutoL.revest == 0) {
        if (cellIndex > 0)
            reynolds = state.cells[cellIndex].Rey(state.cells[cellIndex].dutoL.a, vel1, completionFluidDensity, visc);
        else
            reynolds = state.cells[cellIndex - 1].Rey(state.cells[cellIndex].dutoL.a, vel1, completionFluidDensity, visc);
    } else {
        double dhid = 4 * area / perimeter;
        if (cellIndex > 0)
            reynolds = state.cells[cellIndex - 1].Rey(dhid, vel1, completionFluidDensity, visc);
        else
            reynolds = state.cells[cellIndex].Rey(dhid, vel1, completionFluidDensity, visc);
    }
    frictionFactor = state.cells[cellIndex - 1].fric(reynolds, state.cells[cellIndex].dutoL.rug / diameter);
    gradfric = state.cells[cellIndex - 1].dPdLFric * 0.5 * frictionFactor * completionFluidDensity * (fabs(vel1) * vel1) * perimeter * dx / area;
    gradhidro = state.cells[cellIndex - 1].dPdLHidro * 9.82 * sin(state.cells[cellIndex].dutoL.teta) * completionFluidDensity * dx;
    difpres += (gradfric + gradhidro) / 98066.5;
    return difpres;
}

void updateSteadyGasPressure(const GasLiftState &state, int cellIndex) {

    double dx = 0.5 * state.gasCells[cellIndex].dxL;
    double diameter = state.gasCells[cellIndex].dutoL.a;
    double area = 0.25 * M_PI * diameter * diameter;
    double perimeter = state.gasCells[cellIndex].dutoL.peri;
    double rhog = state.gasCells[cellIndex - 1].flui.MasEspGas(state.gasCells[cellIndex - 1].pres, state.gasCells[cellIndex - 1].temp);

    double VGasmedL;
    VGasmedL = state.gasCells[cellIndex - 1].VGasR;
    double vel1 = VGasmedL / state.gasCells[cellIndex - 1].u1L;

    double visc = state.gasCells[cellIndex - 1].flui.ViscGas(state.gasCells[cellIndex - 1].pres, state.gasCells[cellIndex - 1].temp);

    double reynolds;
    if (state.gasCells[cellIndex].dutoL.revest == 0)
        reynolds = state.gasCells[cellIndex - 1].Rey(state.gasCells[cellIndex].dutoL.a, vel1, rhog, visc);
    else {
        double dhid = 4 * area / perimeter;
        reynolds = state.gasCells[cellIndex - 1].Rey(dhid, vel1, rhog, visc);
    }
    double frictionFactor = state.gasCells[cellIndex - 1].fric(reynolds, state.gasCells[cellIndex].dutoL.rug / diameter);
    double gradfric = state.gasCells[cellIndex - 1].dPdLFric * 0.5 * frictionFactor * rhog * (fabs(vel1) * vel1) * perimeter * dx / area;
    double gradhidro = state.gasCells[cellIndex - 1].dPdLHidro * 9.82 * sin(state.gasCells[cellIndex].dutoL.teta) * rhog * dx;
    state.gasCells[cellIndex - 1].termoFric = gradfric / dx;
    state.gasCells[cellIndex - 1].termoHidro = gradhidro / dx;

    double pmed = state.gasCells[cellIndex - 1].pres - (gradfric + gradhidro) / 98066.5;

    double tmed;
    if (state.steadyIteration != 0)
        tmed = (state.gasCells[cellIndex].dxL * state.gasCells[cellIndex - 1].temp + state.gasCells[cellIndex].dx0 * state.gasCells[cellIndex].temp) /
               (state.gasCells[cellIndex].dxL + state.gasCells[cellIndex].dx0);
    else
        tmed = state.gasCells[cellIndex - 1].temp;
    tmed = state.gasCells[cellIndex - 1].temp;

    dx = 0.5 * state.gasCells[cellIndex].dx0;
    diameter = state.gasCells[cellIndex].duto.a;
    area = 0.25 * M_PI * diameter * diameter;
    perimeter = state.gasCells[cellIndex].duto.peri;
    rhog = state.gasCells[cellIndex].flui.MasEspGas(pmed, tmed);
    double VarArea = 1 / pow(state.gasCells[cellIndex].dutoL.area, 2.) - 1 / pow(state.gasCells[cellIndex].duto.area, 2.);
    double dpDina = 0.5 * VGasmedL * VGasmedL * VarArea / rhog;

    VGasmedL = state.gasCells[cellIndex - 1].VGasR;
    vel1 = VGasmedL / (rhog * area);

    visc = state.gasCells[cellIndex].flui.ViscGas(pmed, tmed);

    if (state.gasCells[cellIndex].duto.revest == 0)
        reynolds = state.gasCells[cellIndex].Rey(state.gasCells[cellIndex].duto.a, vel1, rhog, visc);
    else {
        double dhid = 4 * area / perimeter;
        reynolds = state.gasCells[cellIndex].Rey(dhid, vel1, rhog, visc);
    }
    frictionFactor = state.gasCells[cellIndex].fric(reynolds, state.gasCells[cellIndex].duto.rug / diameter);
    gradfric = state.gasCells[cellIndex].dPdLFric * 0.5 * frictionFactor * rhog * (fabs(vel1) * vel1) * perimeter * dx / area;
    gradhidro = state.gasCells[cellIndex].dPdLHidro * 9.82 * sin(state.gasCells[cellIndex].duto.teta) * rhog * dx;

    state.gasCells[cellIndex].pres = pmed - (gradfric + gradhidro - dpDina) / 98066.5;
    state.gasCells[cellIndex - 1].presR = state.gasCells[cellIndex].pres;
    if (cellIndex < state.gasCellCount)
        state.gasCells[cellIndex + 1].presL = state.gasCells[cellIndex].pres;
    state.gasCells[cellIndex].presini = state.gasCells[cellIndex].pres;
}

void computeSteadyGasFlowRate(const GasLiftState &state, int cellIndex) {

    int nvalv = state.input.nvalvgas;
    int match = 0;
    for (int valveIndex = 0; valveIndex < nvalv; valveIndex++) {
        if (state.gasValveCellIndices[valveIndex] == cellIndex) {
            match = 1;
            state.gasLiftChokes[valveIndex].presEstag = state.gasCells[state.gasValveCellIndices[valveIndex]].pres;
            state.gasLiftChokes[valveIndex].presGarg = (state.cells[state.productionValveCellIndices[valveIndex]].pres - state.gasLiftChokes[valveIndex].presEstag * state.gasLiftChokes[valveIndex].frec) /
                                   (1. - state.gasLiftChokes[valveIndex].frec);
            state.gasLiftChokes[valveIndex].tempEstag = state.gasCells[state.gasValveCellIndices[valveIndex]].temp;
            if (state.cells[state.productionValveCellIndices[valveIndex]].pres < state.gasCells[state.gasValveCellIndices[valveIndex]].pres)
                state.gasCells[state.gasValveCellIndices[valveIndex]].massfonteCH = state.gasLiftChokes[valveIndex].massica();
            else
                state.gasCells[state.gasValveCellIndices[valveIndex]].massfonteCH = 0.;
            if (state.gasLiftChokes[valveIndex].tipo == 1) {
                double abre = calibratedValveArea(state.gasLiftChokes[valveIndex].pcalib * 14.223595, state.gasLiftChokes[valveIndex].tcalib, (state.gasLiftChokes[valveIndex].presEstag - 1.033211) * 14.223595,
                                           (state.gasLiftChokes[valveIndex].presGarg - 1.033211) * 14.223595, state.gasLiftChokes[valveIndex].dextern, state.gasLiftChokes[valveIndex].areagarg,
                                           state.gasLiftChokes[valveIndex].areagarg / state.gasLiftChokes[valveIndex].areafole,
                                           1.8 * state.gasLiftChokes[valveIndex].tempEstag + 32);
                state.gasCells[state.gasValveCellIndices[valveIndex]].massfonteCH *= abre;
            }

            if (cellIndex > 0)
                state.gasCells[cellIndex].VGasR = state.gasCells[cellIndex - 1].VGasR - state.gasCells[cellIndex].massfonteCH;
            else
                state.gasCells[cellIndex].VGasR = state.gasCells[cellIndex].massfonteCH;

            int posGLP = state.productionValveCellIndices[valveIndex];
            int posGLG = state.gasValveCellIndices[valveIndex];
            double gasTemperature;
            state.cells[posGLP].acsr.injg.QGas = state.gasCells[posGLG].massfonteCH * 86400. /
                                            (1.225 * state.gasCells[posGLG].flui.Deng); // celulaG[posGLG].flui.MasEspGas(1.,15.);
            if (state.gasLiftChokes[valveIndex].presEstag > state.gasLiftChokes[valveIndex].presGarg) {
                gasTemperature = state.temperatureUpdater.gasLiftDischargeTemperature(valveIndex);
            } else
                gasTemperature = state.cells[posGLP].temp;
            state.cells[posGLP].acsr.injg.temp = gasTemperature;
            if (state.cells[posGLP].acsr.injg.temp < -50)
                state.cells[posGLP].acsr.injg.temp = -50;
            state.cells[posGLP].fontemassGR = state.gasCells[posGLG].massfonteCH;

            if (state.gasCells[cellIndex].VGasR <= 0. && cellIndex < state.gasCellCount) {
                state.gasCells[cellIndex].VGasR = 0.;
                for (int laterValveIndex = valveIndex + 1; laterValveIndex < nvalv; laterValveIndex++) {
                    state.gasCells[state.gasValveCellIndices[laterValveIndex]].massfonteCH = 0;
                    state.cells[state.productionValveCellIndices[laterValveIndex]].fontemassGR = 0.;
                }
            }

            state.gasCells[cellIndex - 1].VGasRR = state.gasCells[cellIndex].VGasR;
            if (cellIndex < state.gasCellCount)
                state.gasCells[cellIndex + 1].VGasL = state.gasCells[cellIndex].VGasR;
        }
    }
    if (match == 0) {
        if (cellIndex > 0)
            state.gasCells[cellIndex].VGasR = state.gasCells[cellIndex - 1].VGasR - 0 * state.gasCells[cellIndex].massfonteCH;
        else
            state.gasCells[cellIndex].VGasR = state.gasCells[cellIndex - 1].VGasR + 0 * state.gasCells[cellIndex].massfonteCH;
        state.gasCells[cellIndex - 1].VGasRR = state.gasCells[cellIndex].VGasR;
        if (cellIndex < state.gasCellCount)
            state.gasCells[cellIndex + 1].VGasL = state.gasCells[cellIndex].VGasR;
    }
}

void initializeSteadyValveGasFlowRate(const GasLiftState &state, int cellIndex) {
    int nvalv = state.input.nvalvgas;
    double vazvalv;
    double presRev = state.initialGasPressure;
    double tempRev = state.gasCells[0].calor.Textern1;
    for (int valveIndex = 0; valveIndex < nvalv; valveIndex++) {
        int posGLP = state.productionValveCellIndices[valveIndex];
        int posGLG = state.gasValveCellIndices[valveIndex];
        if (cellIndex == posGLP) {
            if (state.gasCells[0].tipoCC == 1) {
                vazvalv = (state.input.gasinj.vazgas[0] * state.gasCells[0].flui.MasEspGas(1., 15.) / 86400.) / nvalv;
                state.gasCells[state.gasValveCellIndices[valveIndex]].massfonteCH = vazvalv;
                state.cells[posGLP].acsr.injg.QGas = state.gasCells[posGLG].massfonteCH * 86400. /
                                                state.gasCells[posGLG].flui.MasEspGas(1., 15.);
                state.cells[posGLP].acsr.injg.temp = tempRev;
                state.cells[posGLP].fontemassGR = state.gasCells[posGLG].massfonteCH;
            } else {
                presRev = state.gasCells[posGLG].pres;
                state.gasLiftChokes[valveIndex].presEstag = presRev;
                state.gasLiftChokes[valveIndex].presGarg = state.cells[state.productionValveCellIndices[valveIndex]].pres;
                state.gasLiftChokes[valveIndex].tempEstag = tempRev;
                if (state.cells[state.productionValveCellIndices[valveIndex]].pres < presRev)
                    state.gasCells[state.gasValveCellIndices[valveIndex]].massfonteCH = state.gasLiftChokes[valveIndex].massica();
                else
                    state.gasCells[state.gasValveCellIndices[valveIndex]].massfonteCH = 0.;
                state.cells[posGLP].acsr.injg.QGas = state.gasCells[posGLG].massfonteCH * 86400. /
                                                state.gasCells[posGLG].flui.MasEspGas(1., 15.);
                state.cells[posGLP].acsr.injg.temp = tempRev;
                state.cells[posGLP].fontemassGR = state.gasCells[posGLG].massfonteCH;
            }
        }
    }
}

void updateSteadyGasTemperature(const GasLiftState &state, int cellIndex) {
    if (state.thermalSourceDisabled == 0) {
        if (((cellIndex) < state.tubingAnnulusStart || (cellIndex) > state.tubingAnnulusEnd)) {
            double dx = state.gasCells[cellIndex].dx0;
            double dxmed = 0.5 * (state.gasCells[cellIndex].dx0 + state.gasCells[cellIndex - 1].dx0);
            double dTdLMed = (state.gasCells[cellIndex].dx0 * state.gasCells[cellIndex].dTdLCor + state.gasCells[cellIndex - 1].dx0 * state.gasCells[cellIndex - 1].dTdLCor) / (state.gasCells[cellIndex].dx0 + state.gasCells[cellIndex - 1].dx0);
            double area = state.gasCells[cellIndex].duto.area;
            double ugsmed;
            double rhog = state.gasCells[cellIndex - 1].flui.MasEspGas(state.gasCells[cellIndex - 1].pres, state.gasCells[cellIndex - 1].temp);
            ugsmed = state.gasCells[cellIndex - 1].VGasR / state.gasCells[cellIndex - 1].u1L;
            double gasSpecificHeat = state.gasCells[cellIndex - 1].flui.CalorGas(state.gasCells[cellIndex - 1].pres, state.gasCells[cellIndex - 1].temp);
            double gasJouleThomson = state.gasCells[cellIndex - 1].flui.JTG(state.gasCells[cellIndex - 1].pres, state.gasCells[cellIndex - 1].temp);
            double hidro = (rhog * ugsmed) * area * 9.82 * sin(state.gasCells[cellIndex].duto.teta);

            state.gasCells[cellIndex - 1].calor.Tint = state.gasCells[cellIndex - 1].temp;
            state.gasCells[cellIndex - 1].calor.Vint = ugsmed;
            state.gasCells[cellIndex - 1].calor.kint = state.gasCells[cellIndex - 1].flui.CondGas(state.gasCells[cellIndex - 1].pres, state.gasCells[cellIndex - 1].temp);
            state.gasCells[cellIndex - 1].calor.cpint = gasSpecificHeat;
            state.gasCells[cellIndex - 1].calor.rhoint = rhog;
            state.gasCells[cellIndex - 1].calor.viscint = state.gasCells[cellIndex - 1].flui.ViscGas(state.gasCells[cellIndex - 1].pres, state.gasCells[cellIndex - 1].temp) * 1.e-3;
            state.gasCells[cellIndex - 1].fluxcal = state.gasCells[cellIndex - 1].calor.transperm();
            double flucal = state.gasCells[cellIndex - 1].fluxcal;
            double flucalCol = 0.;

            double razdx;
            if (cellIndex > 0)
                razdx = dx / (dx + state.gasCells[cellIndex - 1].dx0);
            else
                razdx = 0.5;
            double coefdxT = rhog * ugsmed * gasSpecificHeat * area;
            double coefdxP = rhog * ugsmed * gasJouleThomson * area;
            double dpdx;
            dpdx = 2. * (state.gasCells[cellIndex].pres - ((1 - razdx) * state.gasCells[cellIndex - 1].pres + razdx * state.gasCells[cellIndex].pres)) * 98600. / dx;

            double cinetico;
            double deljmix = 0.;
            double rhomix = rhog;
            double ugsmed0 = state.gasCells[cellIndex - 1].VGasL / state.gasCells[cellIndex - 1].u1L;
            deljmix = (ugsmed - ugsmed0) / dx;

            if (state.input.nCompTotalUnidadesG / dx < 1e6)
                cinetico = rhomix * area * ugsmed * ugsmed * deljmix;
            else
                cinetico = 0.;

            double fontemassG = 0.;
            double fontemassL = 0.;


            if (ugsmed > 1e-3 && (state.gasCells[cellIndex].duto.a / state.gasCells[cellIndex - 1].duto.a > 0.5 &&
                                  state.gasCells[cellIndex - 1].duto.a / state.gasCells[cellIndex].duto.a > 0.5)) {
                double parcenerg1 = dTdLMed * (coefdxP * dpdx - cinetico - hidro + fontemassL + fontemassG) / coefdxT;
                double parcenerg2 = dTdLMed * (flucal + flucalCol) / coefdxT;
                int npasso;
                double dxpasso;
                if (dxmed / state.gasCells[cellIndex - 1].calor.resGlob < 1000.) {
                    npasso = 0.;
                    dxpasso = dxmed;
                } else {
                    npasso = 2 * (dxmed / state.gasCells[cellIndex - 1].calor.resGlob) / 1000 + 1;
                    dxpasso = dxmed / npasso;
                }
                double temppasso = state.gasCells[cellIndex - 1].temp;

                temppasso = dxpasso * (-(-state.gasCells[cellIndex - 1].temp) / dxpasso + parcenerg1 + parcenerg2);
                for (int stepIndex = 1; stepIndex < npasso; stepIndex++) {
                    state.gasCells[cellIndex - 1].calor.Tint = temppasso;
                    state.gasCells[cellIndex - 1].calor.Vint = ugsmed;
                    state.gasCells[cellIndex - 1].calor.kint = state.gasCells[cellIndex - 1].flui.CondGas(state.gasCells[cellIndex - 1].pres, temppasso);
                    state.gasCells[cellIndex - 1].calor.cpint = gasSpecificHeat;
                    state.gasCells[cellIndex - 1].calor.rhoint = rhog;
                    state.gasCells[cellIndex - 1].calor.viscint = state.gasCells[cellIndex - 1].flui.ViscGas(state.gasCells[cellIndex - 1].pres, temppasso) * 1.e-3;
                    state.gasCells[cellIndex - 1].fluxcal = state.gasCells[cellIndex - 1].calor.transperm();
                    flucal = state.gasCells[cellIndex - 1].fluxcal;
                    parcenerg2 = dTdLMed * (flucal + flucalCol) / coefdxT;
                    temppasso = dxpasso * (-(-temppasso) / dxpasso + parcenerg1 + parcenerg2);
                }

                state.gasCells[cellIndex].temp = temppasso;
            } else {
                state.gasCells[cellIndex].temp = state.gasCells[cellIndex].calor.Textern1;
            }

            if (state.gasCells[cellIndex].temp < -50.)
                state.gasCells[cellIndex].temp = -50.;
            if (state.gasCells[cellIndex].temp > 200.)
                state.gasCells[cellIndex].temp = 200.;

            if (cellIndex > 0)
                state.gasCells[cellIndex - 1].tempR = state.gasCells[cellIndex].temp;
            if (cellIndex < state.gasCellCount)
                state.gasCells[cellIndex + 1].tempL = state.gasCells[cellIndex].temp;
        } else {
            if (state.networkCoupled == 1) {
                int kconecte = (cellIndex)-state.tubingAnnulusStart;
                int iconecte = state.annulusTubingStart - kconecte;
                int iconecte2 = iconecte;
                if (iconecte2 == state.lastCell)
                    iconecte2 -= 1;

                state.gasCells[cellIndex].temp = state.cells[iconecte2].calor.resGlob * state.cells[iconecte].calor.fluxFim + state.cells[iconecte].temp;

                if (state.input.correcaoContracorPerm == 1) {

                    double dxmed = 0.5 * (state.gasCells[cellIndex].dx0 + state.gasCells[cellIndex - 1].dx0);
                    double gasSpecificHeat = state.gasCells[cellIndex - 1].flui.CalorGas(state.gasCells[cellIndex - 1].pres, state.gasCells[cellIndex - 1].temp);
                    double invCap = 1. / (state.cells[iconecte].calor.cpint * state.cells[iconecte].MC) + 1. / (gasSpecificHeat * state.gasCells[cellIndex - 1].VGasR);
                    double templog = -exp(-(1. / state.cells[iconecte2].calor.resGlob) * (1 + 0 * state.cells[iconecte].duto.peri) * dxmed * invCap) *
                                         (state.cells[iconecte2 + 1].temp - state.gasCells[cellIndex - 1].temp) +
                                     state.cells[iconecte2].temp;

                    if (templog < state.gasCells[cellIndex].temp)
                        state.gasCells[cellIndex].temp = templog;
                }
            }

            if (state.gasCells[cellIndex].temp < -50.)
                state.gasCells[cellIndex].temp = -50.;
            if (state.gasCells[cellIndex].temp > 200.)
                state.gasCells[cellIndex].temp = 200.;

            if (cellIndex > 0)
                state.gasCells[cellIndex - 1].tempR = state.gasCells[cellIndex].temp;
            if (cellIndex < state.gasCellCount)
                state.gasCells[cellIndex + 1].tempL = state.gasCells[cellIndex].temp;

            double ugsmed;
            ugsmed = state.gasCells[cellIndex - 1].VGasR / state.gasCells[cellIndex - 1].u1L;
            double rhog = state.gasCells[cellIndex - 1].flui.MasEspGas(state.gasCells[cellIndex - 1].pres, state.gasCells[cellIndex - 1].temp);
            double gasSpecificHeat = state.gasCells[cellIndex - 1].flui.CalorGas(state.gasCells[cellIndex - 1].pres, state.gasCells[cellIndex - 1].temp);

            state.gasCells[cellIndex].calor.Tint = state.gasCells[cellIndex].temp;
            state.gasCells[cellIndex].calor.Vint = ugsmed;
            state.gasCells[cellIndex].calor.kint = state.gasCells[cellIndex - 1].flui.CondGas(state.gasCells[cellIndex - 1].pres, state.gasCells[cellIndex - 1].temp);
            state.gasCells[cellIndex].calor.cpint = gasSpecificHeat;
            state.gasCells[cellIndex].calor.rhoint = rhog;
            state.gasCells[cellIndex].calor.viscint = state.gasCells[cellIndex - 1].flui.ViscGas(state.gasCells[cellIndex - 1].pres, state.gasCells[cellIndex - 1].temp) * 1.e-3;
            state.gasCells[cellIndex].fluxcal = state.gasCells[cellIndex].calor.transperm();
        }
    } else {
        state.gasCells[cellIndex].temp = state.gasCells[cellIndex].calor.Textern1;
    }
}

}  // namespace sisprod::gaslift
