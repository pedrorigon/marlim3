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
    for (int i = 1; i <= state.gasCellCount; i++) {
        double A0 = state.gasCells[i - 1].duto.area;
        double dx0 = 0.5 * state.gasCells[i].dxL;
        double A1 = state.gasCells[i].duto.area;
        double dx1 = 0.5 * state.gasCells[i].dx0;
        pmed -= rho0 * 9.81 * dx0 * sin(state.gasCells[i - 1].duto.teta) / 98066.52;
        tmed = state.gasCells[i].calor.Textern1;
        if (i < state.interfaceCell)
            rho1 = state.gasCells[i].flui.MasEspGas(pmed, tmed);
        else
            rho1 = state.gasCells[i].MasEspFlu(pmed, tmed);
        pmed -= rho1 * 9.81 * dx1 * sin(state.gasCells[i].duto.teta) / 98066.52;
        rho0 = rho1;

        state.gasCells[i].pres = pmed;
        state.gasCells[i].presini = pmed;
        state.gasCells[i - 1].presR = pmed;
        state.gasCells[i].temp = tmed;
        state.gasCells[i - 1].tempR = tmed;
        state.gasCells[i].u1L = state.gasCells[i].duto.area * rho0;
        state.gasCells[i - 1].u1R = state.gasCells[i].u1L;
        state.gasCells[i].VGasR = 0;
        state.gasCells[i - 1].VGasRR = 0;
        state.gasCells[i].massfonteCH = 0.;
        state.gasCells[i - 1].u1R = state.gasCells[i].u1L;
        if (i < state.gasCellCount) {
            state.gasCells[i + 1].presL = pmed;
            state.gasCells[i + 1].tempL = tmed;
            state.gasCells[i + 1].u1LL = state.gasCells[i].u1L;
            state.gasCells[i + 1].VGasL = 0;
        }
    }
}

void updateGasLine(const GasLiftState &state) {
    for (int i = 0; i <= state.gasCellCount; i++) {
        if (i != 0 && i != state.gasCellCount) {
            state.gasCells[i].pres = state.gasFreeTerms[3 * i];
            state.gasCells[i].presL = state.gasFreeTerms[3 * i - 3];
            state.gasCells[i].presR = state.gasFreeTerms[3 * i + 3];
            state.gasCells[i].VGasR = state.gasFreeTerms[3 * i + 1];
            state.gasCells[i].VGasL = state.gasFreeTerms[3 * i - 2];
            state.gasCells[i].VGasRR = state.gasFreeTerms[3 * i + 4];
        } else if (i == 0) {
            state.gasCells[i].pres = state.gasFreeTerms[3 * i];
            state.gasCells[i].presL = state.gasCells[i].pres;
            state.gasCells[i].presR = state.gasFreeTerms[3 * i + 3];
            state.gasCells[i].VGasR = state.gasFreeTerms[3 * i + 1];
            state.gasCells[i].VGasRR = state.gasFreeTerms[3 * i + 4];
            double auxpres = state.gasCells[i].presR;
        } else {
            state.gasCells[i].pres = state.gasFreeTerms[3 * i];
            state.gasCells[i].presL = state.gasFreeTerms[3 * i - 3];
            state.gasCells[i].VGasR = state.gasFreeTerms[3 * i + 1];
            state.gasCells[i].VGasL = state.gasFreeTerms[3 * i - 2];
            state.gasCells[i].presR = state.gasCells[i].pres;
        }
    }
}

void updateBufferedGasLine(const GasLiftState &state) {
    for (int i = 0; i <= state.gasCellCount; i++) {
        if (i != 0 && i != state.gasCellCount) {
            state.gasCells[i].VGasRBuf = state.gasFreeTerms[3 * i + 1];
        } else if (i == 0) {
            state.gasCells[i].VGasRBuf = state.gasFreeTerms[3 * i + 1];
        } else {
            state.gasCells[i].VGasRBuf = state.gasFreeTerms[3 * i + 1];
        }
    }
}

double calibratedValveArea(double PCal, double TCal, double PVO, double PT,
                           double dextern, double areagarg, double Rvalv, double Temp) {
    // PCal: Calibration pressure.
    // PVO: Casing pressure.
    // PT: Tubing pressure.
    // Rvalv: Area ratio.
    // Temp: Bottom-hole temperature in Fahrenheit.

    double PB80 = PCal * (1 - Rvalv);
    PB80 = (PB80 + 14.6959488) * (80 + 460.67) / (TCal * 1.8 + 491.67) - 14.6959488;
    double PBT = PB80 * (1 + 0.00215 * (Temp - 80));
    double compara = PVO * (1 - Rvalv) + PT * Rvalv;
    double abertura = 0.;
    if (compara > PBT)
        abertura = 1.;

    double areafol = areagarg / Rvalv;
    double BSR;
    if (dextern * 100. / 2.54 > 1.1)
        BSR = 500.0 * areafol;
    else
        BSR = 1950.0 * areafol;
    double XMVS = ((PVO - PBT) * areafol - (PVO - PT) * areagarg) / BSR;

    // IF THE VALVE IS CLOSED, QG = 0
    if (XMVS <= 0.0)
        abertura = 0.;

    double DP = sqrt(areagarg * 4. / M_PI);
    double RP = DP / 2.0;
    double RB = sqrt(areafol / M_PI);
    double APE = M_PI * RP * XMVS * (XMVS + 2.0 * sqrt(RB * RB - RP * RP));
    APE = APE / sqrt((XMVS + sqrt(RB * RB - RP * RP)) * (XMVS + sqrt(RB * RB - RP * RP)) + RP * RP);
    if (APE > areagarg)
        APE = areagarg;
    abertura = APE / areagarg;

    return abertura;
}

double unloadingPressureCorrection(const GasLiftState &state, double vazmax, int ivalv, double fator, int sinal) {
    int iG = state.gasValveCellIndices[ivalv];
    int iP = state.productionValveCellIndices[ivalv];
    double pmed = state.cells[iP].pres;
    state.gasLiftChokes[ivalv].presEstag = state.gasCells[iG].pres;
    state.gasLiftChokes[ivalv].tempEstag = state.gasCells[iG].temp;
    double rho0 = state.gasCells[iG].MasEspFlu(state.gasLiftChokes[ivalv].presEstag, state.gasLiftChokes[ivalv].tempEstag);
    double massica = (fator * state.input.vazDescControl - vazmax) * rho0;
    double precorr = 0.;
    precorr = pow(massica / state.gasLiftChokes[ivalv].areagarg, 2.) / (2. * rho0 * 98066.52);
    return sinal * precorr;
}

double computeUnloadingValvePressure(const GasLiftState &state, double vazGarg, int ivalv) {

    double velmax = 0;
    double laz1 = 0.1;
    double laz2 = 0.4;

    double pmed;
    double tmed;
    double massica;
    pmed = state.gasSurfacePressure;
    state.input.presMaxDesc = 100000.;
    for (int i = state.lastCell; i >= 0; i--) {

        double tmed = state.cells[i].temp;
        double A1 = state.cells[i].duto.area;
        double S1 = state.cells[i].duto.peri;
        double dx1 = state.cells[i].dx;
        double rhoG = state.cells[i].flui.MasEspGas(pmed, tmed);
        double rhoP = state.cells[i].flui.MasEspLiq(pmed, tmed);
        double rhoC = state.cells[i].fluicol.MasEspFlu(pmed, tmed);
        double viscG = state.cells[i].flui.ViscGas(pmed, tmed);
        double viscP = state.cells[i].flui.ViscOleo(pmed, tmed);
        double viscC = state.cells[i].fluicol.VisFlu(pmed, tmed);
        double alf = state.cells[i].alf;
        double bet = state.cells[i].bet;
        double rholiq = bet * rhoC + (1. - bet) * rhoP;
        double viscliq = bet * viscC + (1. - bet) * viscP;
        double rhomix = alf * rhoG + (1. - alf) * rholiq;
        double viscmix = alf * viscG + (1. - alf) * viscliq;
        double vel1 = state.cells[i].QL / (A1 * rholiq) + state.cells[i].QG / (A1 * rhoG);
        double re1;
        if (fabs(vel1) > 1e-15) {
            if (state.cells[i].duto.revest == 0)
                re1 = state.cells[i].Rey(state.cells[i].duto.a, vel1, rhomix, viscmix);
            else {
                double dhid = 4 * A1 / S1;
                re1 = state.cells[i].Rey(dhid, vel1, rhomix, viscmix);
            }
        }
        double f1;
        if (fabs(vel1) > 1e-15)
            f1 = state.cells[i].fric(re1, state.cells[i].duto.rug / state.cells[i].duto.a);
        else
            f1 = 0.;
        double tens1 = f1 * rhomix * vel1 * fabs(vel1) / 2.;
        pmed -= (-9.82 * rhomix * sin(state.cells[i].duto.teta) - tens1 * S1 / A1) * dx1 / 98066.5;
        if (state.cells[i].acsr.tipo == 3) {
            double auxpresmax = state.cells[i].acsr.ipr.Pres - (pmed - state.gasSurfacePressure);
            if (auxpresmax < state.input.presMaxDesc)
                state.input.presMaxDesc = auxpresmax;
        }
    }
    if (state.gasSurfacePressure >= state.input.presMaxDesc * 0.9999999 && vazGarg > (1 - laz1) * state.input.vazDescControl && state.gasCells[0].VGasR > 0.) {
        double precorr = unloadingPressureCorrection(state, vazGarg, ivalv, (1 - laz1), -1);
        if (fabs(precorr) > 0.01 * state.initialGasPressure * state.cells[0].dt)
            precorr = (fabs(precorr) / precorr) * 0.01 * state.initialGasPressure * state.cells[0].dt;
        state.initialGasPressure += precorr;
        if (state.initialGasPressure < state.input.presMinDescG) {
            state.initialGasPressure = state.input.presMinDescG;
        }

    } else if (vazGarg <= (1 - laz2) * state.input.vazDescControl && state.gasSurfacePressure <= state.input.presMinDesc * 1.0000001) {
        double precorr = unloadingPressureCorrection(state, vazGarg, ivalv, 1 - laz2, 1);
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
    int nvalv = state.input.nvalvgas;
    for (int i = state.interfaceCell; i <= state.gasCellCount; i++) {
        double dx0 = 0.5 * state.gasCells[i].dxL;
        double dx1 = 0.5 * state.gasCells[i].dx0;
        double RgasR = 1.;
        if (state.gasCells[i].razInter <= 0.5)
            RgasR = 2 * state.gasCells[i].razInter;
        double RgasL = 0.;
        if (state.gasCells[i - 1].razInter >= 0.5)
            RgasL = 2 * (state.gasCells[i - 1].razInter - 0.5);
        double LGasL = dx0 * RgasL;
        double LGasR = dx1 * RgasR;
        double LLiqL = dx0 - LGasL;
        double LLiqR = dx1 - LGasR;
        double temp = state.gasCells[i - 1].temp;
        double pres = state.gasCells[i - 1].pres;
        double rhoL = state.gasCells[i].MasEspFlu(pres, temp);
        double viscL = state.gasCells[i].VisFlu(pres, temp);
        double rhoG = state.gasCells[i].flui.MasEspGas(pres, temp);
        double viscG = state.gasCells[i].flui.ViscGas(pres, temp);
        double vel1 = state.gasCells[i].VGasL / (rhoL * state.gasCells[i - 1].duto.area);
        if (state.gasCells[i].razInter > (*state.globals).localtiny)
            vel1 = state.gasCells[i].VGasL / (rhoG * state.gasCells[i - 1].duto.area);
        double vel2 = state.gasCells[i].VGasL / (rhoL * state.gasCells[i].duto.area);
        if (state.gasCells[i].razInter > (*state.globals).localtiny)
            vel2 = state.gasCells[i].VGasL / (rhoG * state.gasCells[i].duto.area);
        double re1G;
        double re1L;
        double re2G;
        double re2L;
        if (state.gasCells[i - 1].duto.revest == 0)
            re1L = state.gasCells[i - 1].Rey(state.gasCells[i - 1].duto.a, vel1, rhoL, viscL);
        else {
            double dhid = 4 * state.gasCells[i - 1].duto.area / state.gasCells[i - 1].duto.peri;
            re1L = state.gasCells[i - 1].Rey(dhid, vel1, rhoL, viscL);
        }
        if (state.gasCells[i].duto.revest == 0)
            re2L = state.gasCells[i].Rey(state.gasCells[i].duto.a, vel2, rhoL, viscL);
        else {
            double dhid = 4 * state.gasCells[i].duto.area / state.gasCells[i].duto.peri;
            re2L = state.gasCells[i].Rey(dhid, vel2, rhoL, viscL);
        }
        if (state.gasCells[i - 1].duto.revest == 0)
            re1G = state.gasCells[i - 1].Rey(state.gasCells[i - 1].duto.a, vel1, rhoG, viscG);
        else {
            double dhid = 4 * state.gasCells[i - 1].duto.area / state.gasCells[i - 1].duto.peri;
            re1G = state.gasCells[i - 1].Rey(dhid, vel1, rhoG, viscG);
        }
        if (state.gasCells[i].duto.revest == 0)
            re2G = state.gasCells[i].Rey(state.gasCells[i].duto.a, vel2, rhoG, viscG);
        else {
            double dhid = 4 * state.gasCells[i].duto.area / state.gasCells[i].duto.peri;
            re2G = state.gasCells[i].Rey(dhid, vel2, rhoG, viscG);
        }
        double f1L = state.gasCells[i - 1].fric(re1L, state.gasCells[i - 1].duto.rug / state.gasCells[i - 1].duto.a) * LLiqL;
        double f2L = state.gasCells[i].fric(re2L, state.gasCells[i].duto.rug / state.gasCells[i].duto.a) * LLiqR;
        double hidro1L = 1 * (9.82 * sin(state.gasCells[i - 1].duto.teta) * rhoL) * LLiqL;
        double hidro2L = 1 * (9.82 * sin(state.gasCells[i].duto.teta) * rhoL) * LLiqR;
        double f1G = state.gasCells[i - 1].fric(re1G, state.gasCells[i - 1].duto.rug / state.gasCells[i - 1].duto.a) * LGasL;
        double f2G = state.gasCells[i].fric(re2G, state.gasCells[i].duto.rug / state.gasCells[i].duto.a) * LGasR;
        double hidro1G = 1 * (9.82 * sin(state.gasCells[i - 1].duto.teta) * rhoG) * LGasL;
        double hidro2G = 1 * (9.82 * sin(state.gasCells[i].duto.teta) * rhoG) * LGasR;
        state.gasCells[i].pres = state.gasCells[i - 1].pres + (-0.5 * (f1L * rhoL + f1G * rhoG) * vel1 * fabs(vel1) * state.gasCells[i - 1].duto.peri / state.gasCells[i - 1].duto.area - 0.5 * (f2L * rhoL + f2G * rhoG) * vel2 * fabs(vel2) * state.gasCells[i].duto.peri / state.gasCells[i].duto.area - hidro1L - hidro2L - hidro1G - hidro2G) / 98066.52;
        state.gasCells[i].presL = state.gasCells[i - 1].pres;
        state.gasCells[i - 1].presR = state.gasCells[i].pres;

        state.temperatureUpdater.dischargeTemperature(i);

        state.gasCells[i].u1L = ((1. - state.gasCells[i].razInter) * state.gasCells[i].MasEspFlu(state.gasCells[i].pres, state.gasCells[i].temp) + state.gasCells[i].razInter * state.gasCells[i].flui.MasEspGas(pres, temp)) * state.gasCells[i].duto.area;
        state.gasCells[i - 1].u1R = state.gasCells[i].u1L;
        state.gasCells[i].u1LL = state.gasCells[i - 1].u1L;
    }

    double Qtotal = 0.;
    for (int i = state.interfaceCell; i <= state.gasCellCount; i++) {
        double temp = state.gasCells[i].temp;
        double pres = state.gasCells[i].pres;
        double rhoL = state.gasCells[i].MasEspFlu(pres, temp);
        double rhoG = state.gasCells[i].flui.MasEspGas(pres, temp);
        double qfonte = state.gasCells[i].massfonteCH / rhoL;
        if (state.gasCells[i].razInter > 0.5)
            qfonte = state.gasCells[i].massfonteCH / rhoG;
        Qtotal += qfonte;
    }
    double temp = state.gasCells[state.interfaceCell].temp;
    double pres = state.gasCells[state.interfaceCell].pres;
    double rhoL = state.gasCells[state.interfaceCell].MasEspFlu(pres, temp);
    double rhoG = state.gasCells[state.interfaceCell].flui.MasEspGas(pres, temp);
    state.gasCells[state.interfaceCell].VGasL = Qtotal * rhoG;
    state.gasCells[state.interfaceCell - 1].VGasR = state.gasCells[state.interfaceCell].VGasL;
    state.gasCells[state.interfaceCell - 2].VGasRR = state.gasCells[state.interfaceCell].VGasL;
    for (int i = state.interfaceCell; i <= state.gasCellCount; i++) {
        double temp = state.gasCells[i].temp;
        double pres = state.gasCells[i].pres;
        double rhoL = state.gasCells[i].MasEspFlu(pres, temp);
        double rhoG = state.gasCells[i].flui.MasEspGas(pres, temp);
        double qfonte = state.gasCells[i].massfonteCH / rhoL;
        if (state.gasCells[i].razInter > 0.5)
            qfonte = state.gasCells[i].massfonteCH / rhoG;
        Qtotal -= qfonte;
        if (i < state.gasCellCount) {
            state.gasCells[i + 1].VGasL = (Qtotal)*rhoL;
            state.gasCells[i].VGasR = state.gasCells[i + 1].VGasL;
        } else
            state.gasCells[i].VGasR = 0.;
        state.gasCells[i - 1].VGasRR = state.gasCells[i].VGasL;
    }

    state.interfaceVelocity = state.gasCells[state.interfaceCell + 1].VGasL / (state.gasCells[state.interfaceCell + 1].MasEspFlu(state.gasCells[state.interfaceCell + 1].pres, state.gasCells[state.interfaceCell + 1].temp) * state.gasCells[state.interfaceCell].duto.area);
}

void advanceInterface(const GasLiftState &state) {

    state.gasCells[state.interfaceCell].razInter = (state.gasCells[state.interfaceCell].razInterIni * state.gasCells[state.interfaceCell].dx0 + state.interfaceVelocity * state.timeStep) / state.gasCells[state.interfaceCell].dx0;
    if (state.interfaceCell == (state.gasCellCount - 1) && state.gasCells[state.interfaceCell].razInter >= 0.99) {

        double pmed = state.gasCells[state.interfaceCell].pres;
        double tmed = state.gasCells[state.interfaceCell].temp;
        double A1 = state.gasCells[state.interfaceCell].duto.area;
        double rho1 = state.gasCells[state.interfaceCell].flui.MasEspGas(pmed, tmed);
        double rhoL = state.gasCells[state.interfaceCell].MasEspFlu(pmed, tmed);
        state.gasCells[state.interfaceCell].VGasR = state.gasCells[state.interfaceCell].VGasR * rho1 / rhoL;
        state.gasCells[state.interfaceCell].u1L = A1 * rho1;
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
        A1 = state.gasCells[state.interfaceCell].duto.area;
        rho1 = state.gasCells[state.interfaceCell].flui.MasEspGas(pmed, tmed);
        rhoL = state.gasCells[state.interfaceCell].MasEspFlu(pmed, tmed);
        state.gasCells[state.interfaceCell].VGasR = 0 * state.gasCells[state.interfaceCell].VGasR * rho1 / rhoL;
        state.gasCells[state.interfaceCell].u1L = A1 * rho1;
        state.gasCells[state.interfaceCell].u1R = A1 * rho1;
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

}  // namespace sisprod::gaslift
