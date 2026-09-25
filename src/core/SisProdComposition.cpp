#include "SisProdComposition.h"

#include "Leitura.h"
#include "celula3.h"
#include "variaveisGlobais1D.h"

#include <math.h>

namespace sisprod::composition {

void storePreviousFractionsAndMovePigs(const CompositionState &state) {

    state.cells[0].alfini = state.cells[0].alf;
    state.cells[0].betini = state.cells[0].bet;
    state.cells[0].alfPigDini = state.cells[0].alfPigD;
    state.cells[0].betPigDini = state.cells[0].betPigD;
    state.cells[0].alfPigEini = state.cells[0].alfPigE;
    state.cells[0].betPigEini = state.cells[0].betPigE;
    state.cells[1].alfLini = state.cells[0].alfini;
    state.cells[1].betLini = state.cells[0].betini;

    if (state.input.ConContEntrada == 0) {
        state.cells[0].alfLini = state.cells[0].alfini;
        state.cells[0].betLini = state.cells[0].betini;
    } else {
        state.cells[0].alfLini = state.inletQuality;
        state.cells[0].betLini = state.inletCompletionFraction;
    }

    state.previousMovingPigCount = state.movingPigCount;
    state.movingPigCount = 0;

    for (int i = 0; i <= state.lastCell; i++) {
        state.cells[i].estadoPigini = state.cells[i].estadoPig;
        state.cells[i].alfini = state.cells[i].alf;
        if (i > 0)
            state.cells[i - 1].alfRini = state.cells[i].alfini;
        if (i < state.lastCell)
            state.cells[i + 1].alfLini = state.cells[i].alfini;
        state.cells[i].betini = state.cells[i].bet;
        if (i > 0)
            state.cells[i - 1].betRini = state.cells[i].betini;
        if (i < state.lastCell)
            state.cells[i + 1].betLini = state.cells[i].betini;

        state.cells[i].alfPigDini = state.cells[i].alfPigD;
        state.cells[i].betPigDini = state.cells[i].betPigD;
        state.cells[i].alfPigEini = state.cells[i].alfPigE;
        state.cells[i].betPigEini = state.cells[i].betPigE;

        state.cells[i].alfPigERini = state.cells[i].alfPigER;
        state.cells[i].betIini = state.cells[i].betI;
        state.cells[i].betRIini = state.cells[i].betRI;
        state.cells[i].betLIini = state.cells[i].betLI;

        if (i > 0)
            state.cells[i - 1].alfPigER = state.cells[i].alfPigE;
        if (i == state.lastCell)
            state.cells[i].alfPigER = state.cells[i].alf;
        state.cells[i].DelPig = 0.;
        state.cells[i].RazAreaPig = 0.;
        state.cells[i].cdpig = 1.;
        if (state.cells[i].estadoPig == 1) {
            state.movingPigCount++;
            int ipig = state.cells[i].indpig;
            state.cells[i].DelPig = state.input.pig[ipig].delpres;
            state.cells[i].RazAreaPig = state.input.pig[ipig].razarea;
            state.cells[i].cdpig = state.input.pig[ipig].cdPig;
            double AC = state.cells[i].duto.area;
            double jL = (state.cells[i].QL + state.cells[i].QG) / AC;
            double jR;
            if (i < state.lastCell)
                jR = (state.cells[i + 1].QL + state.cells[i + 1].QG) / AC;
            else
                jR = (state.cells[i].QL + state.cells[i].QG) / AC;
            state.cells[i].velPigini = state.cells[i].velPig;
            state.cells[i].velPig = jL * state.cells[i].razPig + jR * (1. - state.cells[i].razPig) - state.cells[i].VazaPig / AC;
            for (int j = 0; j < state.scheduledPigCount; j++) {
                if (state.pigReceiverCells[j] == i) {
                    state.cells[i].velPig = 0.;
                    state.cells[i].estadoPig = 0;
                    state.cells[i].razPig = 0.;
                    state.cells[i].razPigini = 0.;
                    state.cells[i].alfPigDini = state.cells[i].alf;
                    state.cells[i].betPigDini = state.cells[i].bet;
                    state.cells[i].alfPigEini = state.cells[i].alf;
                    state.cells[i].betPigEini = state.cells[i].bet;
                    state.cells[i].alfPigD = state.cells[i].alf;
                    state.cells[i].betPigD = state.cells[i].bet;
                    state.cells[i].alfPigE = state.cells[i].alf;
                    state.cells[i].betPigE = state.cells[i].bet;
                    if (i > 0)
                        state.cells[i - 1].alfPigER = state.cells[i].alfPigE;
                    if (i == state.lastCell)
                        state.cells[i].alfPigER = state.cells[i].alf;

                    state.movingPigCount--;
                }
            }
        }
        state.cells[i].razPigini = state.cells[i].razPig;
    }
    state.cells[state.lastCell].alfRini = state.cells[state.lastCell].alfini;
    state.cells[state.lastCell].betRini = state.cells[state.lastCell].betini;
}

void cacheCellAndFaceDensities(const CompositionState &state) {

    state.cells[0].rpC = state.cells[0].flui.MasEspLiq(state.cells[0].pres, state.cells[0].temp);
    state.cells[0].rgC = state.cells[0].flui.MasEspGas(state.cells[0].pres, state.cells[0].temp /*,1*/);
    state.cells[0].rcC = state.cells[0].fluicol.MasEspFlu(state.cells[0].pres, state.cells[0].temp);
    state.cells[0].rpL = state.cells[0].rpC;
    state.cells[0].rgL = state.cells[0].rgC;
    state.cells[0].rcL = state.cells[0].rcC;

    state.cells[0].rpCi = state.cells[0].rpC;
    state.cells[0].rgCi = state.cells[0].rgC;
    state.cells[0].rcCi = state.cells[0].rcC;
    state.cells[0].rpLi = state.cells[0].rpCi;
    state.cells[0].rgLi = state.cells[0].rgCi;
    state.cells[0].rcLi = state.cells[0].rcCi;

    state.cells[0].mipC = state.cells[0].flui.ViscOleo(state.cells[0].pres, state.cells[0].temp);
    state.cells[0].migC = state.cells[0].flui.ViscGas(state.cells[0].pres, state.cells[0].temp /*,1*/);
    state.cells[0].micC = state.cells[0].fluicol.VisFlu(state.cells[0].pres, state.cells[0].temp);

#pragma omp parallel for num_threads((*state.globals).ntrd)
    for (int i = 1; i <= state.lastCell; i++) {
        double p;
        double t;
        p = state.cells[i].pres;
        t = state.cells[i].temp;
        state.cells[i].rpC = state.cells[i].flui.MasEspLiq(p, t);
        state.cells[i].rgC = state.cells[i].flui.MasEspGas(p, t);
        state.cells[i].rcC = state.cells[i].fluicol.MasEspFlu(p, t);

        state.cells[i].mipC = state.cells[i].flui.ViscOleo(p, t);
        state.cells[i].migC = state.cells[i].flui.ViscGas(p, t);
        state.cells[i].micC = state.cells[i].fluicol.VisFlu(p, t);

        double tmed = state.cells[i - 1].temp;
        if (state.cells[i].VTemper < 0.)
            tmed = state.cells[i].temp;
        ProFlu flu;
        if (state.cells[i].QL < 0.)
            flu = state.cells[i].flui;
        else
            flu = state.cells[i - 1].flui;
        state.cells[i].rpCi = flu.MasEspLiq(state.cells[i].presaux, tmed);
        state.cells[i].rgCi = flu.MasEspGas(state.cells[i].presaux, tmed);
        state.cells[i].rcCi = state.cells[i].fluicol.MasEspFlu(state.cells[i].presaux, tmed);
    }
    for (int i = 1; i <= state.lastCell; i++) {

        state.cells[i - 1].rpR = state.cells[i].rpC;
        state.cells[i - 1].rgR = state.cells[i].rgC;
        state.cells[i - 1].rcR = state.cells[i].rcC;
        state.cells[i].rpL = state.cells[i - 1].rpC;
        state.cells[i].rgL = state.cells[i - 1].rgC;
        state.cells[i].rcL = state.cells[i - 1].rcC;

        state.cells[i - 1].mipR = state.cells[i].mipC;
        state.cells[i - 1].migR = state.cells[i].migC;
        state.cells[i - 1].micR = state.cells[i].micC;

        state.cells[i - 1].rpRi = state.cells[i].rpCi;
        state.cells[i - 1].rgRi = state.cells[i].rgCi;
        state.cells[i - 1].rcRi = state.cells[i].rcCi;
        state.cells[i].rpLi = state.cells[i - 1].rpCi;
        state.cells[i].rgLi = state.cells[i - 1].rgCi;
        state.cells[i].rcLi = state.cells[i - 1].rcCi;
    }
    state.cells[state.lastCell].rpR = state.cells[state.lastCell].rpC;
    state.cells[state.lastCell].rgR = state.cells[state.lastCell].rgC;
    state.cells[state.lastCell].rcR = state.cells[state.lastCell].rcC;

    state.cells[state.lastCell].mipR = state.cells[state.lastCell].mipC;
    state.cells[state.lastCell].migR = state.cells[state.lastCell].migC;
    state.cells[state.lastCell].micR = state.cells[state.lastCell].micC;
}

void evaluateWaxDeposition(const CompositionState &state) {
    for (int i = 0; i <= state.lastCell; i++) {
        if (i < state.lastCell)
            state.cells[i].WaxDeposition(state.input.detalParafina, state.lastCell);
        else {
            if (state.surfaceChokeMassFlag == 0 || state.cells[state.lastCell].Mliqini > 0) {
                double cpDep = state.cells[i - 1].duto.cp[0];
                double kDep = state.cells[i - 1].duto.cond[0];
                double rhoDep = state.cells[i - 1].duto.rhoC[0];
                state.cells[i].deltaPar = state.cells[i - 1].deltaPar;
                if (state.cells[i].parafinado == 0 && state.cells[i].deltaPar > 0.) {
                    state.cells[i].duto.atualizaCamada(state.cells[i].deltaPar, state.input.detalParafina.rug, cpDep, kDep, rhoDep);
                    state.cells[i].calor.atualiza(state.cells[i].duto, 1);
                    state.cells[i].parafinado = 1;
                } else if (state.cells[i].deltaPar > 0.) {
                    state.cells[i].duto.atualizaCamada2(state.cells[i].deltaPar, cpDep, kDep, rhoDep);
                    state.cells[i].calor.atualiza2(state.cells[i].duto);
                }
            } else
                state.cells[i].WaxDeposition(state.input.detalParafina, state.lastCell);
        }
    }
}

}  // namespace sisprod::composition
