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

}  // namespace sisprod::gaslift
