#include "SisProdSteadyState.h"

#include "GradientCorrelations.h"
#include "Leitura.h"
#include "SisProdConstants.h"
#include "celula3.h"
#include "celulaGas.h"
#include "chokegas.h"
#include "variaveisGlobais1D.h"

#include <math.h>

namespace sisprod::steady {

void correctGasSpecificGravity(const SteadyStateState &state, int i) {
    if (state.input.corrDeng == 0) { // nesta chave, para o caso black oil, se faz a disntincao
        // entre a densidade do gas dissolvido e do gas livre
        // rDgD= razao entre a densidade do gas dissolvido e o gas nas condicoes standard
        // rDgL= razao entre a densidade do gas livre e o gas nas condicoes standard
        // arq.corrDeng==0 implica em nÃƒÂ£o haver distincao entre as densidades
        state.cells[i].flui.rDgD = 1.;
        state.cells[i].flui.rDgL = 1.;
        state.cells[i].flui.PCis = state.cells[i].flui.PC;
        state.cells[i].flui.TCis = state.cells[i].flui.TC;
        if (state.cells[i].acsr.tipo == 1 && state.cells[i].acsr.injg.seco == 0) {
            state.cells[i].acsr.injg.FluidoPro.rDgD = 1.;
            state.cells[i].acsr.injg.FluidoPro.rDgL = 1.;
        }
        if (state.cells[i].acsr.tipo == 2) {
            state.cells[i].acsr.injl.FluidoPro.rDgD = 1.;
            state.cells[i].acsr.injl.FluidoPro.rDgL = 1.;
        } else if (state.cells[i].acsr.tipo == 3) {
            state.cells[i].acsr.ipr.FluidoPro.rDgD = 1.;
            state.cells[i].acsr.ipr.FluidoPro.rDgL = 1.;
        } else if (state.cells[i].acsr.tipo == 10) {
            state.cells[i].acsr.injm.FluidoPro.rDgD = 1.;
            state.cells[i].acsr.injm.FluidoPro.rDgL = 1.;
        }
        if (state.cells[i].acsr.tipo == 15) {
            state.cells[i].acsr.radialPoro.flup.rDgD = 1.;
            state.cells[i].acsr.radialPoro.flup.rDgL = 1.;
            for (int porousCell = 0; porousCell < state.cells[i].acsr.radialPoro.ncel; porousCell++) {
                state.cells[i].acsr.radialPoro.celula[porousCell].flup.rDgD = 1.;
                state.cells[i].acsr.radialPoro.celula[porousCell].flup.rDgL = 1.;
            }
        }
        if (state.cells[i].acsr.tipo == 16) {
            state.cells[i].acsr.poroso2D.dados.flup.rDgD = 1.;
            state.cells[i].acsr.poroso2D.dados.flup.rDgL = 1.;
            for (int porousCell = 0; porousCell < state.cells[i].acsr.poroso2D.dados.transfer.ncel; porousCell++) {
                state.cells[i].acsr.poroso2D.dados.transfer.celula[porousCell].flup.rDgD = 1.;
                state.cells[i].acsr.poroso2D.dados.transfer.celula[porousCell].flup.rDgL = 1.;
            }
            for (int porousCell = 0; porousCell < state.cells[i].acsr.poroso2D.malha.nele; porousCell++) {
                state.cells[i].acsr.poroso2D.malha.mlh2d[porousCell].flup.rDgD = 1.;
                state.cells[i].acsr.poroso2D.malha.mlh2d[porousCell].flup.rDgL = 1.;
            }
        }
    } else { // caso em que arq.corrDeng==1 que implica em haver distincao entre as densidades
        int k = i - 1;
        if (state.steadyIteration != 0)
            k = i;
        state.cells[i].flui.razDegD(state.cells[k].pres, state.cells[k].temp);
        state.cells[i].flui.rzDegL(state.cells[k].pres, state.cells[k].temp);
        state.cells[i].flui.PcTcIS(); // como a densidade de gas livre mudou, suas pressoes
        // criticas devem ser recalculadas
        // caso tenha alguma fonte, se calcula as densidades in situ das fontes
        if (state.cells[i].acsr.tipo == 1 && state.cells[i].acsr.injg.seco == 0) { // fonte de liquido
            state.cells[i].acsr.injg.FluidoPro.razDegD(state.cells[k].pres, state.cells[k].temp);
            state.cells[i].acsr.injg.FluidoPro.rzDegL(state.cells[k].pres, state.cells[k].temp);
        }
        if (state.cells[i].acsr.tipo == 2) { // fonte de liquido
            state.cells[i].acsr.injl.FluidoPro.razDegD(state.cells[k].pres, state.cells[k].temp);
            state.cells[i].acsr.injl.FluidoPro.rzDegL(state.cells[k].pres, state.cells[k].temp);
        } else if (state.cells[i].acsr.tipo == 3) { // ipr
            state.cells[i].acsr.ipr.FluidoPro.razDegD(state.cells[k].pres, state.cells[k].temp);
            state.cells[i].acsr.ipr.FluidoPro.rzDegL(state.cells[k].pres, state.cells[k].temp);
            ;
        } else if (state.cells[i].acsr.tipo == 10) { // fonte de massa generica
            state.cells[i].acsr.injm.FluidoPro.razDegD(state.cells[k].pres, state.cells[k].temp);
            state.cells[i].acsr.injm.FluidoPro.rzDegL(state.cells[k].pres, state.cells[k].temp);
        }
        if (state.cells[i].acsr.tipo == 15) {
            state.cells[i].acsr.radialPoro.flup.razDegD(state.cells[i].pres, state.cells[i].temp);
            state.cells[i].acsr.radialPoro.flup.rzDegL(state.cells[i].pres, state.cells[i].temp);
            for (int porousCell = 0; porousCell < state.cells[i].acsr.radialPoro.ncel; porousCell++) {
                double pres = state.cells[i].acsr.radialPoro.celula[porousCell].Pcamada;
                double temp = state.cells[i].acsr.radialPoro.tRes;
                state.cells[i].acsr.radialPoro.celula[porousCell].flup.razDegD(pres, temp);
                state.cells[i].acsr.radialPoro.celula[porousCell].flup.rzDegL(pres, temp);
            }
        }
        if (state.cells[i].acsr.tipo == 16) {
            state.cells[i].acsr.poroso2D.dados.flup.razDegD(state.cells[i].pres, state.cells[i].temp);
            state.cells[i].acsr.poroso2D.dados.flup.rzDegL(state.cells[i].pres, state.cells[i].temp);
            for (int porousCell = 0; porousCell < state.cells[i].acsr.poroso2D.dados.transfer.ncel; porousCell++) {
                double pres = state.cells[i].acsr.poroso2D.dados.transfer.celula[porousCell].Pcamada;
                double temp = state.cells[i].acsr.poroso2D.dados.transfer.tRes;
                state.cells[i].acsr.poroso2D.dados.transfer.celula[porousCell].flup.razDegD(pres, temp);
                state.cells[i].acsr.poroso2D.dados.transfer.celula[porousCell].flup.rzDegL(pres, temp);
            }
            for (int porousCell = 0; porousCell < state.cells[i].acsr.poroso2D.malha.nele; porousCell++) {
                double pres = state.cells[i].acsr.poroso2D.malha.mlh2d[porousCell].cel2D.presC;
                double temp = state.cells[i].acsr.poroso2D.malha.mlh2d[porousCell].tRes;
                state.cells[i].acsr.poroso2D.malha.mlh2d[porousCell].flup.razDegD(pres, temp);
                state.cells[i].acsr.poroso2D.malha.mlh2d[porousCell].flup.rzDegL(pres, temp);
            }
        }
    }
}

namespace {

/// Resolves the phase fractions of cell i for the two reverse mass marches.
///
/// Four regimes -- no flow, liquid only, gas only, two-phase -- and the drift
/// closure the two-phase arm needs. advanceReverseSteadyMass and
/// advanceReverseCompositionalSteadyMass carried these 182 lines twice, with no
/// difference at all between the copies.
///
/// T084 section 5.1 recommended sharing this and T086 moved both variants
/// without doing it, so it is done here rather than left as a duplicate the
/// stage itself created.
///
/// What comes AFTER this block is not shared and must not be: that is H15 of
/// evidencia/renovamass-diff.md, where Rev zeroes the completion-fluid
/// residence time and CompRev propagates it from the cell upstream.
void applyReverseSteadyPhaseFractions(const SteadyStateState &state, int i, double rhol, double rhog) {
    if (fabs(state.cells[i].QG + state.cells[i].QL) < (*state.globals).localtiny) {
        if (state.input.tipoFluido == 1) {
            state.cells[i].alf = 0.;
        } else {
            state.cells[i].alf = 1.;
        }
        state.cells[i].alfini = state.cells[i].alf;
        state.cells[i - 1].alfR = state.cells[i].alf;
        state.cells[i - 1].alfRini = state.cells[i].alf;
        if (i < state.lastCell)
            state.cells[i + 1].alfL = state.cells[i].alf;
        if (i < state.lastCell)
            state.cells[i + 1].alfLini = state.cells[i].alf;
        state.cells[i].alfPigD = state.cells[i].alf;
        state.cells[i].alfPigDini = state.cells[i].alf;
        state.cells[i].alfPigE = state.cells[i].alf;
        state.cells[i].alfPigEini = state.cells[i].alf;
        if (i < state.lastCell) {
            state.cells[i + 1].betL = state.cells[i].bet;
            state.cells[i + 1].betLini = state.cells[i].bet;
            state.cells[i + 1].betL = state.cells[i].bet;
            state.cells[i + 1].betLini = state.cells[i].bet;
            state.cells[i + 1].betLI = state.cells[i].bet;
        }
        state.cells[i].betini = state.cells[i].bet;
        state.cells[i - 1].betR = state.cells[i].bet;
        state.cells[i - 1].betRini = state.cells[i].bet;
        state.cells[i].betPigD = state.cells[i].bet;
        state.cells[i].betPigDini = state.cells[i].bet;
        state.cells[i].betPigE = state.cells[i].bet;
        state.cells[i].betPigEini = state.cells[i].bet;
        state.cells[i].betI = state.cells[i].bet;
        state.cells[i - 1].betRI = state.cells[i].bet;
    } else if (fabs(state.cells[i].QG) < (*state.globals).localtiny) { // caso so exista liquido:
        state.cells[i].alf = 0.;
        state.cells[i].alfini = state.cells[i].alf;
        state.cells[i - 1].alfR = state.cells[i].alf;
        state.cells[i - 1].alfRini = state.cells[i].alf;
        if (i < state.lastCell)
            state.cells[i + 1].alfL = state.cells[i].alf;
        if (i < state.lastCell)
            state.cells[i + 1].alfLini = state.cells[i].alf;
        state.cells[i].alfPigD = state.cells[i].alf;
        state.cells[i].alfPigDini = state.cells[i].alf;
        state.cells[i].alfPigE = state.cells[i].alf;
        state.cells[i].alfPigEini = state.cells[i].alf;
        if (i < state.lastCell) {
            state.cells[i + 1].betL = state.cells[i].bet;
            state.cells[i + 1].betLini = state.cells[i].bet;
            state.cells[i + 1].betL = state.cells[i].bet;
            state.cells[i + 1].betLini = state.cells[i].bet;
            state.cells[i + 1].betLI = state.cells[i].bet;
        }
        state.cells[i].betini = state.cells[i].bet;
        state.cells[i - 1].betR = state.cells[i].bet;
        state.cells[i - 1].betRini = state.cells[i].bet;
        state.cells[i].betPigD = state.cells[i].bet;
        state.cells[i].betPigDini = state.cells[i].bet;
        state.cells[i].betPigE = state.cells[i].bet;
        state.cells[i].betPigEini = state.cells[i].bet;
        state.cells[i].betI = state.cells[i].bet;
        state.cells[i - 1].betRI = state.cells[i].bet;
    } else if (fabs(state.cells[i].QL) < (*state.globals).localtiny * 1e-6) { // caso so exista gas:
        state.cells[i].alf = 1.;
        state.cells[i].alfini = state.cells[i].alf;
        state.cells[i - 1].alfR = state.cells[i].alf;
        state.cells[i - 1].alfRini = state.cells[i].alf;
        if (i < state.lastCell)
            state.cells[i + 1].alfL = state.cells[i].alf;
        if (i < state.lastCell)
            state.cells[i + 1].alfLini = state.cells[i].alf;
        state.cells[i].alfPigD = state.cells[i].alf;
        state.cells[i].alfPigDini = state.cells[i].alf;
        state.cells[i].alfPigE = state.cells[i].alf;
        state.cells[i].alfPigEini = state.cells[i].alf;
        state.cells[i].c0 = 1.;
        state.cells[i].ud = 0.;
        if (i < state.lastCell) {
            state.cells[i + 1].betL = state.cells[i].bet;
            state.cells[i + 1].betLini = state.cells[i].bet;
            state.cells[i + 1].betL = state.cells[i].bet;
            state.cells[i + 1].betLini = state.cells[i].bet;
            state.cells[i + 1].betLI = state.cells[i].bet;
        }
        state.cells[i].betini = state.cells[i].bet;
        state.cells[i - 1].betR = state.cells[i].bet;
        state.cells[i - 1].betRini = state.cells[i].bet;
        state.cells[i].betPigD = state.cells[i].bet;
        state.cells[i].betPigDini = state.cells[i].bet;
        state.cells[i].betPigE = state.cells[i].bet;
        state.cells[i].betPigEini = state.cells[i].bet;
        state.cells[i].betI = state.cells[i].bet;
        state.cells[i - 1].betRI = state.cells[i].bet;
    } else { // caso bifasico
        double c0 = 1.;
        double ud = 0.;
        if (fabs(state.cells[i].QL) > (*state.globals).localtiny * 1e-6) {
            if (state.steadyIteration == 0) { // primeira estimativa, primeira iteracao
                // utiliza-se a fracao de vazio sem escorregamento, pois a propria correlacao para se obter a
                // fracao de vazio depende do valor da fracao de vazio
                if ((fabs(state.cells[i].QG) + fabs(state.cells[i].QL)) > (*state.globals).localtiny) {
                    state.cells[i].alf = fabs(state.cells[i].QG) /
                                    (fabs(state.cells[i].QG) + fabs(state.cells[i].QL));
                    state.updaters.steadyDriftClosure(i, c0, ud);
                } else
                    state.cells[i].alf = 0.;
                state.cells[i].alfini = state.cells[i].alf;
                state.cells[i - 1].alfR = state.cells[i].alf;
                state.cells[i - 1].alfRini = state.cells[i].alf;
                if (i < state.lastCell)
                    state.cells[i + 1].alfL = state.cells[i].alf;
                if (i < state.lastCell)
                    state.cells[i + 1].alfLini = state.cells[i].alf;
                state.cells[i].alfPigD = state.cells[i].alf;
                state.cells[i].alfPigDini = state.cells[i].alf;
                state.cells[i].alfPigE = state.cells[i].alf;
                state.cells[i].alfPigEini = state.cells[i].alf;
            }
            if (fabs(rhog) / rhol > 0.9) {
                c0 = 1.;
                ud = 0.;
            }

            // para o caso permanente, a fracao de vazio Ã© obtida a partir das relacoes de escorregamento
            // portanto, e neste ponto que se obtem Co e Ud:
            else if (fabs(state.cells[i].QG) > (*state.globals).localtiny && fabs(state.cells[i].QL) > (*state.globals).localtiny * 1e-6)
                state.updaters.steadyDriftClosure(i, c0, ud);
            state.cells[i].c0 = c0;
            state.cells[i].ud = ud;
            double area = state.cells[i].duto.area;
            if (fabs(state.cells[i].QG + state.cells[i].QL) > (*state.globals).localtiny) {
                // alfa com escorregamento:
                state.cells[i].alf = state.cells[i].QG / (c0 * (state.cells[i].QG + state.cells[i].QL) + ud * area);
                double alfHomo = state.cells[i].QG / (state.cells[i].QG + state.cells[i].QL);
                if (state.cells[i].alf > 1. - 1e-15 || state.cells[i].alf < 1e-15)
                    state.cells[i].alf = alfHomo;

            } else
                state.cells[i].alf = 0.;
            if (state.cells[i].alf > (1 - (*state.globals).localtiny))
                state.cells[i].alf = state.cells[i].QG / (state.cells[i].QG + state.cells[i].QL);
        } else {
            c0 = 1.;
            ud = 0.;
            state.cells[i].alf = 1.;
        }
        if (state.cells[i].alf < 0.)
            state.cells[i].alf = 0.;
        else if (state.cells[i].alf > 1.)
            state.cells[i].alf = 1.;
        // atualizacoes dos valores das fracoes volumetricas da celula i armazendadas em
        // outras celulas, e inclusiove armazendo os valores para "tempo anterior", que nao
        // sao relevantes para o problema permanente mas importantes se o resultado permanente
        // der partida na solucao transiente:
        state.cells[i].alfini = state.cells[i].alf;
        state.cells[i - 1].alfR = state.cells[i].alf;
        state.cells[i - 1].alfRini = state.cells[i].alf;
        if (i < state.lastCell)
            state.cells[i + 1].alfL = state.cells[i].alf;
        if (i < state.lastCell)
            state.cells[i + 1].alfLini = state.cells[i].alf;
        state.cells[i].alfPigD = state.cells[i].alf;
        state.cells[i].alfPigDini = state.cells[i].alf;
        state.cells[i].alfPigE = state.cells[i].alf;
        state.cells[i].alfPigEini = state.cells[i].alf;
        if (i < state.lastCell) {
            state.cells[i + 1].betL = state.cells[i].bet;
            state.cells[i + 1].betLini = state.cells[i].bet;
            state.cells[i + 1].betL = state.cells[i].bet;
            state.cells[i + 1].betLini = state.cells[i].bet;
            state.cells[i + 1].betLI = state.cells[i].bet;
        }
        state.cells[i].betini = state.cells[i].bet;
        state.cells[i - 1].betR = state.cells[i].bet;
        state.cells[i - 1].betRini = state.cells[i].bet;
        state.cells[i].betPigD = state.cells[i].bet;
        state.cells[i].betPigDini = state.cells[i].bet;
        state.cells[i].betPigE = state.cells[i].bet;
        state.cells[i].betPigEini = state.cells[i].bet;
        state.cells[i].betI = state.cells[i].bet;
        state.cells[i - 1].betRI = state.cells[i].bet;
    }
}

}  // namespace

void advanceReverseSteadyMass(const SteadyStateState &state, int i) {
    int mudaRGO = 1;
    if (state.input.flashCompleto == 1)
        mudaRGO = 1;
    ProFlu fluF;
    double residenceTimeSource = 0.;
    if (i == 1) {
        state.updaters.updateSource(i - 1); // metodo que verifica se existe uma fonte na celula e calcula o valor das
        // vazoes massicas de liquido produzido (oleo+agua), gas e liquido complementar
        // relacao entre fonte a esquerda e a direita de uma celula:
        state.cells[i].fontemassLL = state.cells[i - 1].fontemassLR;
        state.cells[i].fontemassCL = state.cells[i - 1].fontemassCR;
        state.cells[i].fontemassGL = state.cells[i - 1].fontemassGR;

        state.cells[i].MC = state.cells[i - 1].fontemassLR + state.cells[i - 1].fontemassCR + state.cells[i - 1].fontemassGR;

        state.cells[i - 1].MR = state.cells[i].MC;
        if (i < state.lastCell)
            state.cells[i + 1].ML = state.cells[i].MC;
        state.cells[i - 1].MRini = state.cells[i - 1].MR;
        state.cells[i].MComp = state.cells[i - 1].MComp + state.cells[i - 1].fontemassCR;

        if (state.cells[i - 1].acsr.tipo == 1) {
            fluF = state.cells[i - 1].acsr.injg.FluidoPro;
            residenceTimeSource = state.cells[i - 1].acsr.injg.fluidocol.TR;
        } else if (state.cells[i - 1].acsr.tipo == 2) {
            fluF = state.cells[i - 1].acsr.injl.FluidoPro;
            residenceTimeSource = state.cells[i - 1].acsr.injl.fluidocol.TR;
        } else if (state.cells[i - 1].acsr.tipo == 3) {
            fluF = state.cells[i - 1].acsr.ipr.FluidoPro;
            residenceTimeSource = 0.;
        } else if (state.cells[i - 1].acsr.tipo == 15) {
            fluF = state.cells[i - 1].acsr.radialPoro.flup;
            residenceTimeSource = 0.;
        } else if (state.cells[i - 1].acsr.tipo == 16) {
            fluF = state.cells[i - 1].acsr.poroso2D.dados.flup;
            residenceTimeSource = 0.;
        } else if (state.cells[i - 1].acsr.tipo == 10) {
            fluF = state.cells[i - 1].acsr.injm.FluidoPro;
            residenceTimeSource = state.cells[i - 1].acsr.injm.fluidocol.TR;
        } else if (state.cells[i - 1].acsr.tipo == 9 && state.cells[i - 1].acsr.fontechk.abertura > 1e-6) {
            if (state.cells[i - 1].acsr.fontechk.presT > state.cells[i - 1].acsr.fontechk.pamb) {
                fluF = state.cells[i - 1].acsr.fontechk.fluidoP;
            } else {
                fluF = state.cells[i - 1].acsr.fontechk.fluidoPamb;
            }
        }
        state.cells[i].flui = fluF;
        state.cells[i].flui.RenovaFluido();
        correctGasSpecificGravity(state, i);
    } else {
        state.cells[i - 1].fontemassLR = 0.;
        state.cells[i - 1].fontemassCR = 0.;
        state.cells[i - 1].fontemassGR = 0.;

        state.cells[i].fontemassLL = state.cells[i - 1].fontemassLR;
        state.cells[i].fontemassCL = state.cells[i - 1].fontemassCR;
        state.cells[i].fontemassGL = state.cells[i - 1].fontemassGR;

        state.cells[i].MC = state.cells[i - 1].MC;

        state.cells[i - 1].MR = state.cells[i].MC;
        if (i < state.lastCell)
            state.cells[i + 1].ML = state.cells[i].MC;
        state.cells[i - 1].MRini = state.cells[i - 1].MR;

        state.cells[i].MComp = state.cells[i - 1].MComp;

        state.cells[i].flui = state.cells[i - 1].flui;
        state.cells[i].flui.RenovaFluido();
        correctGasSpecificGravity(state, i);
    }
    if (i == 0) {
        if (state.cells[0].acsr.tipo == 15)
            state.cells[0].flui.BSW = state.cells[0].acsr.radialPoro.BSW;
        else if (state.cells[0].acsr.tipo == 16)
            state.cells[0].flui.BSW = state.cells[0].acsr.poroso2D.dados.transfer.BSW;
    }

    double tmed;
    // temperatura na fronteira entre a celula i-1 e a celula i, da segunda iteracao em diante, pode-se
    // usar a temperatura da celula i, pois ja a tem calculada, mas isto pode ser um complicador de
    // convergencia, mais seguro manter o criterio em todas as iteracoes, neste caso, a temperatura na
    // fronteira e admitida = a temperatura da celula i-1
    if (state.steadyIteration != 0 && state.input.AceleraConvergPerm == 0)
        tmed = (state.cells[i].dx * state.cells[i].temp + state.cells[i].dxR * state.cells[i].tempR) / (state.cells[i].dx + state.cells[i].dxR);
    else
        tmed = state.cells[i - 1].temp;

    double pmed = state.cells[i].presaux;

    state.cells[i].rpCi = state.cells[i].flui.MasEspLiq(pmed, tmed);
    state.cells[i].rcCi = state.cells[i].fluicol.MasEspFlu(pmed, tmed);
    state.cells[i].rgCi = state.cells[i].flui.MasEspGas(pmed, tmed);

    double MasCarb = state.cells[i].MC - state.cells[i].MComp;
    double MasGas = MasCarb * state.cells[i].flui.FracMassHidra(state.cells[i].presaux, tmed);
    double MasLiqProd = MasCarb - MasGas;
    double oilVolumeFactorInSitu;
    double waterVolumeFactorInSitu;
    if (state.cells[i].flui.RGO < 1e6)
        oilVolumeFactorInSitu = state.cells[i].flui.BOFunc(state.cells[i].presaux, tmed);
    else
        oilVolumeFactorInSitu = 1.;
    waterVolumeFactorInSitu = state.cells[i].flui.BAFunc(state.cells[i].presaux, tmed);
    double waterCutInSitu = state.cells[i].flui.BSW * waterVolumeFactorInSitu / (oilVolumeFactorInSitu + waterVolumeFactorInSitu * state.cells[i].flui.BSW - state.cells[i].flui.BSW * oilVolumeFactorInSitu);
    double oilFlowRate = MasLiqProd * (1 - state.cells[i - 1].FW) / state.cells[i].rpCi;
    double waterFlowRate;
    if (waterCutInSitu < (1 - (*state.globals).localtiny))
        waterFlowRate = waterCutInSitu * oilFlowRate / (1 - waterCutInSitu);
    else
        waterFlowRate = MasLiqProd / state.cells[i].flui.MasEspAgua(pmed, tmed);
    double completionFlowRate;
    completionFlowRate = state.cells[i].MComp / state.cells[i].rcCi;
    ////////////////////////esperar//////////////////////////////////////
    state.cells[i].bet = completionFlowRate / (oilFlowRate + waterFlowRate + completionFlowRate);

    state.cells[i].Mliqini = MasLiqProd;
    state.cells[i - 1].MliqiniR = state.cells[i].Mliqini;
    if (i < state.lastCell)
        state.cells[i + 1].MliqiniL = state.cells[i].Mliqini;

    double rhol = (1 - state.cells[i].bet) * state.cells[i].rpCi + state.cells[i].bet * state.cells[i].rcCi;
    double rhog = state.cells[i].rgCi;
    // vazoes volumetricas:
    state.cells[i].QL = state.cells[i].Mliqini / rhol;
    if (i < state.lastCell)
        state.cells[i + 1].QLL = state.cells[i].QL;
    state.cells[i - 1].QLR = state.cells[i].QL;
    state.cells[i].QG = (state.cells[i].MC - state.cells[i].Mliqini) / rhog;

    if (i < state.lastCell)
        state.cells[i + 1].QLL = state.cells[i].QL;
    state.cells[i - 1].QLR = state.cells[i].QL;

    // definicao de temperaturas maximas e minimas para uma eventual reavaliacao do metodo
    // ASTM quando ocorre mistura de fluidos e se deseja atualizar o modelo
    // de viscosidade de oleo morto se este for o ASTM que trabalha com um par de temperatura
    double temperatureLow;
    if (state.input.flashCompleto == 1 && (state.input.tabent.tmin - 0) > (*state.globals).localtiny)
        temperatureLow = state.input.tabent.tmin + 0.1;
    else
        temperatureLow = 0.;
    double temperatureHigh;
    if (state.input.flashCompleto == 1 && (70 - state.input.tabent.tmax) < (*state.globals).localtiny)
        temperatureHigh = state.input.tabent.tmax - 0.1;
    else
        temperatureHigh = 70.;
    double oilVolumeFactor;
    double waterVolumeFactor;
    double solutionGasRatio;
    // calcula os valores de RS, Bo e Ba na celula anterior a i-esima celula
    if (state.cells[i - 1].flui.RGO < 1e7) {
        solutionGasRatio = state.cells[i - 1].flui.RS(state.cells[i - 1].pres, state.cells[i - 1].temp);
        oilVolumeFactor = state.cells[i - 1].flui.BOFunc(state.cells[i - 1].pres, state.cells[i - 1].temp, solutionGasRatio);
        waterVolumeFactor = state.cells[i - 1].flui.BAFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        solutionGasRatio = solutionGasRatio * 6.29 / 35.31467;
    } else {
        oilVolumeFactor = 1;
        solutionGasRatio = 0;
        waterVolumeFactor = 0.;
    }
    // BSW in-situ da celula anterior, na marcha, a i-esima celula
    state.cells[i - 1].FW = state.cells[i - 1].flui.BSW * waterVolumeFactor / (oilVolumeFactor + waterVolumeFactor * state.cells[i - 1].flui.BSW - state.cells[i - 1].flui.BSW * oilVolumeFactor);
    state.cells[i - 1].FWini = state.cells[i - 1].FW;
    // Definicao das fracoes volumetricas:
    applyReverseSteadyPhaseFractions(state, i, rhol, rhog);

    if (fabs(state.cells[i].QL) > 1e-15 && fabs(state.cells[i].MComp) > 1e-15) {
        double hol0 = 1. - state.cells[i - 1].alf;
        double hol1 = 1. - state.cells[i].alf;
        state.cells[i].fluicol.TR = (state.cells[i - 1].fluicol.TR * state.cells[i - 1].MComp + state.cells[i - 1].fontemassCR * residenceTimeSource) / state.cells[i].MComp;
        state.cells[i].fluicol.TR = state.cells[i].fluicol.TR + (0.5 * state.cells[i - 1].duto.area * state.cells[i - 1].dx * hol0 +
                                                       0.5 * state.cells[i].duto.area * state.cells[i].dx * hol1) /
                                                          state.cells[i].QL;
    } else
        state.cells[i].fluicol.TR = 0.;
}

namespace {

/// Refreshes the compositional flash of cell i, seeded from the nearest cell
/// upstream or downstream that has a usable calculated beta.
///
/// Carried by advanceReverseCompositionalSteadyMass and
/// advanceCompositionalSteadyMass. The two copies differed by four lines, and
/// those four were a debug anchor -- an empty `int para; para = 0;` behind a
/// test on flui.iIER -- of the same class stage 6 removed from the gas-lift
/// module. It has no effect and is not carried in here.
void refreshCompositionalFlashFromNeighbour(const SteadyStateState &state, int i, double pmed, double tmed) {
    if (((state.steadyIteration == 0 && i <= 1 && state.searchOrigin == 0) && ((*state.globals).chaverede == 0 || (*state.globals).iterRede == 1)) && state.input.tabelaDinamica == 0) {
        state.cells[i].flui.atualizaPropComp(pmed, tmed, -1, NULL, NULL, state.input.pocinjec);
    } else if (state.input.tabelaDinamica == 0) {
        if ((state.steadyIteration == 0 && state.searchOrigin == 0) && ((*state.globals).chaverede == 0 || (*state.globals).iterRede == 1)) {
            int veriI = i - 1;
            while ((state.cells[veriI].flui.dCalculatedBeta > 1 - (0.0 + 1e-15) || state.cells[veriI].flui.dCalculatedBeta < (0.0 + 1e-15)) && (veriI > 0))
                veriI--;
            if ((veriI == 0) && (state.cells[veriI].flui.dCalculatedBeta > 1 - (0.0 + 1e-15) || state.cells[veriI].flui.dCalculatedBeta < (0.0 + 1e-15)))
                state.cells[i].flui.atualizaPropComp(pmed, tmed, -1, NULL, NULL, state.input.pocinjec);
            else
                state.cells[i].flui.atualizaPropComp(pmed, tmed, state.cells[veriI].flui.dCalculatedBeta, state.cells[veriI].flui.oCalculatedLiqComposition,
                                                state.cells[veriI].flui.oCalculatedVapComposition, state.input.pocinjec);
        } else {
            int veriI;
            int veriIinf = i;
            while ((state.cells[veriIinf].flui.dCalculatedBeta > 1 - 1e-15 ||
                    state.cells[veriIinf].flui.dCalculatedBeta < 1e-15) &&
                   veriIinf > 0)
                veriIinf--;
            int veriIsup = i;
            while ((state.cells[veriIsup].flui.dCalculatedBeta > 1 - 1e-15 ||
                    state.cells[veriIsup].flui.dCalculatedBeta < 1e-15) &&
                   veriIsup < state.lastCell)
                veriIsup++;
            if ((state.cells[veriIinf].flui.dCalculatedBeta > 1 - 1e-15 ||
                 state.cells[veriIinf].flui.dCalculatedBeta < 1e-15) &&
                (state.cells[veriIsup].flui.dCalculatedBeta < 1. &&
                 state.cells[veriIsup].flui.dCalculatedBeta > 0.))
                veriI = veriIsup;
            else if ((state.cells[veriIsup].flui.dCalculatedBeta > 1 - 1e-15 || state.cells[veriIsup].flui.dCalculatedBeta < 1e-15) &&
                     (state.cells[veriIinf].flui.dCalculatedBeta < 1. && state.cells[veriIinf].flui.dCalculatedBeta > 0.))
                veriI = veriIinf;
            else if (fabs(i - veriIsup) < fabs(i - veriIinf))
                veriI = veriIsup;
            else
                veriI = veriIinf;
            if ((state.cells[veriI].flui.dCalculatedBeta > 1 - (0.0 + 1e-15) || state.cells[veriI].flui.dCalculatedBeta < (0.0 + 1e-15)))
                state.cells[i].flui.atualizaPropComp(pmed, tmed, -1, NULL, NULL, state.input.pocinjec);
            else
                state.cells[i].flui.atualizaPropComp(pmed, tmed, state.cells[veriI].flui.dCalculatedBeta, state.cells[veriI].flui.oCalculatedLiqComposition,
                                                state.cells[veriI].flui.oCalculatedVapComposition, state.input.pocinjec);
        }
    }
}

}  // namespace

void advanceReverseCompositionalSteadyMass(const SteadyStateState &state, int i) {
    int mudaRGO = 1;
    if (state.input.flashCompleto == 1)
        mudaRGO = 1;
    ProFlu fluF;
    double residenceTimeSource = 0.;
    if (i == 1) {
        state.updaters.updateSource(i - 1); // metodo que verifica se existe uma fonte na celula e calcula o valor das
        // vazoes massicas de liquido produzido (oleo+agua), gas e liquido complementar
        // relacao entre fonte a esquerda e a direita de uma celula:
        state.cells[i].fontemassLL = state.cells[i - 1].fontemassLR;
        state.cells[i].fontemassCL = state.cells[i - 1].fontemassCR;
        state.cells[i].fontemassGL = state.cells[i - 1].fontemassGR;

        state.cells[i].MC = state.cells[i - 1].fontemassLR + state.cells[i - 1].fontemassCR + state.cells[i - 1].fontemassGR;

        state.cells[i - 1].MR = state.cells[i].MC;
        if (i < state.lastCell)
            state.cells[i + 1].ML = state.cells[i].MC;
        state.cells[i - 1].MRini = state.cells[i - 1].MR;
        state.cells[i].MComp = state.cells[i - 1].MComp + state.cells[i - 1].fontemassCR;

        if (state.cells[i - 1].acsr.tipo == 1) {
            state.cells[i - 1].acsr.injg.FluidoPro.atualizaPropComp(state.cells[i].pres, state.cells[i].temp, -1, NULL, NULL, state.input.pocinjec);
            fluF = state.cells[i - 1].acsr.injg.FluidoPro;
            residenceTimeSource = state.cells[i - 1].acsr.injg.fluidocol.TR;
        } else if (state.cells[i - 1].acsr.tipo == 2) {
            state.cells[i - 1].acsr.injl.FluidoPro.atualizaPropComp(state.cells[i].pres, state.cells[i].temp, -1, NULL, NULL, state.input.pocinjec);
            fluF = state.cells[i - 1].acsr.injl.FluidoPro;
            residenceTimeSource = state.cells[i - 1].acsr.injl.fluidocol.TR;
        } else if (state.cells[i - 1].acsr.tipo == 3) {
            state.cells[i - 1].acsr.ipr.FluidoPro.atualizaPropComp(state.cells[i].pres, state.cells[i].temp, -1, NULL, NULL, state.input.pocinjec);
            fluF = state.cells[i - 1].acsr.ipr.FluidoPro;
            residenceTimeSource = 0.;
        } else if (state.cells[i - 1].acsr.tipo == 10) {
            state.cells[i - 1].acsr.injm.FluidoPro.atualizaPropComp(state.cells[i].pres, state.cells[i].temp, -1, NULL, NULL, state.input.pocinjec);
            fluF = state.cells[i - 1].acsr.injm.FluidoPro;
            residenceTimeSource = state.cells[i - 1].acsr.injm.fluidocol.TR;
        } else if (state.cells[i - 1].acsr.tipo == 9 && state.cells[i - 1].acsr.fontechk.abertura > 1e-6) {
            state.cells[i - 1].acsr.fontechk.fluidoP.atualizaPropComp(state.cells[i].pres, state.cells[i].temp, -1, NULL, NULL, state.input.pocinjec);
            state.cells[i - 1].acsr.fontechk.fluidoPamb.atualizaPropComp(state.cells[i].pres, state.cells[i].temp, -1, NULL, NULL, state.input.pocinjec);
            if (state.cells[i - 1].acsr.fontechk.presT > state.cells[i - 1].acsr.fontechk.pamb) {
                fluF = state.cells[i - 1].acsr.fontechk.fluidoP;
            } else {
                fluF = state.cells[i - 1].acsr.fontechk.fluidoPamb;
            }
        } else if (state.cells[i - 1].acsr.tipo == 15) {
            state.cells[i - 1].acsr.radialPoro.flup.atualizaPropComp(state.cells[i].pres, state.cells[i].temp, -1, NULL, NULL, state.input.pocinjec);
            fluF = state.cells[i - 1].acsr.radialPoro.flup;
            residenceTimeSource = 0.;
        } else if (state.cells[i - 1].acsr.tipo == 16) {
            state.cells[i - 1].acsr.poroso2D.dados.flup.atualizaPropComp(state.cells[i].pres, state.cells[i].temp, -1, NULL, NULL, state.input.pocinjec);
            fluF = state.cells[i - 1].acsr.poroso2D.dados.flup;
            residenceTimeSource = 0.;
        }
        state.cells[i].flui = fluF;

    } else {
        state.cells[i - 1].fontemassLR = 0.;
        state.cells[i - 1].fontemassCR = 0.;
        state.cells[i - 1].fontemassGR = 0.;

        state.cells[i].fontemassLL = state.cells[i - 1].fontemassLR;
        state.cells[i].fontemassCL = state.cells[i - 1].fontemassCR;
        state.cells[i].fontemassGL = state.cells[i - 1].fontemassGR;

        state.cells[i].MC = state.cells[i - 1].MC;

        state.cells[i - 1].MR = state.cells[i].MC;
        if (i < state.lastCell)
            state.cells[i + 1].ML = state.cells[i].MC;
        state.cells[i - 1].MRini = state.cells[i - 1].MR;

        state.cells[i].MComp = state.cells[i - 1].MComp;

        state.cells[i].flui = state.cells[i - 1].flui;
    }
    if (i == 0) {
        if (state.cells[0].acsr.tipo == 15)
            state.cells[0].flui.BSW = state.cells[0].acsr.radialPoro.BSW;
        else if (state.cells[0].acsr.tipo == 16)
            state.cells[0].flui.BSW = state.cells[0].acsr.poroso2D.dados.transfer.BSW;
    }

    double tmed;
    // temperatura na fronteira entre a celula i-1 e a celula i, da segunda iteracao em diante, pode-se
    // usar a temperatura da celula i, pois ja a tem calculada, mas isto pode ser um complicador de
    // convergencia, mais seguro manter o criterio em todas as iteracoes, neste caso, a temperatura na
    // fronteira e admitida = a temperatura da celula i-1
    if (state.steadyIteration != 0 && state.input.AceleraConvergPerm == 0)
        tmed = (state.cells[i].dx * state.cells[i].temp + state.cells[i].dxR * state.cells[i].tempR) / (state.cells[i].dx + state.cells[i].dxR);
    else
        tmed = state.cells[i - 1].temp;

    double pmed = state.cells[i].presaux;

    refreshCompositionalFlashFromNeighbour(state, i, pmed, tmed);

    double titulo = state.cells[i].flui.dVaporMassFraction;

    state.cells[i].rpCi = state.cells[i].flui.MasEspLiq(pmed, tmed);
    state.cells[i].rcCi = state.cells[i].fluicol.MasEspFlu(pmed, tmed);
    state.cells[i].rgCi = state.cells[i].flui.MasEspGas(pmed, tmed);

    double MasCarb = state.cells[i].MC - state.cells[i].MComp;
    double MasGas = MasCarb * state.cells[i].flui.FracMassHidra(state.cells[i].presaux, tmed);
    double MasLiqProd = MasCarb - MasGas;
    double oilVolumeFactorInSitu;
    double waterVolumeFactorInSitu;
    if (titulo < 1 - 1e-15)
        oilVolumeFactorInSitu = state.cells[i].flui.BOFunc(state.cells[i].presaux, tmed);
    else
        oilVolumeFactorInSitu = 1.;
    waterVolumeFactorInSitu = state.cells[i].flui.BAFunc(state.cells[i].presaux, tmed);
    double waterCutInSitu = state.cells[i].flui.BSW * waterVolumeFactorInSitu / (oilVolumeFactorInSitu + waterVolumeFactorInSitu * state.cells[i].flui.BSW - state.cells[i].flui.BSW * oilVolumeFactorInSitu);
    double oilFlowRate = MasLiqProd * (1 - state.cells[i - 1].FW) / state.cells[i].rpCi;
    double waterFlowRate;
    if (waterCutInSitu < (1 - (*state.globals).localtiny))
        waterFlowRate = waterCutInSitu * oilFlowRate / (1 - waterCutInSitu);
    else
        waterFlowRate = MasLiqProd / state.cells[i].flui.MasEspAgua(pmed, tmed);
    double completionFlowRate;
    completionFlowRate = state.cells[i].MComp / state.cells[i].rcCi;
    ////////////////////////esperar//////////////////////////////////////
    state.cells[i].bet = completionFlowRate / (oilFlowRate + waterFlowRate + completionFlowRate);

    state.cells[i].Mliqini = MasLiqProd;
    state.cells[i - 1].MliqiniR = state.cells[i].Mliqini;
    if (i < state.lastCell)
        state.cells[i + 1].MliqiniL = state.cells[i].Mliqini;

    double rhol = (1 - state.cells[i].bet) * state.cells[i].rpCi + state.cells[i].bet * state.cells[i].rcCi;
    double rhog = state.cells[i].rgCi;
    // vazoes volumetricas:
    state.cells[i].QL = state.cells[i].Mliqini / rhol;
    if (i < state.lastCell)
        state.cells[i + 1].QLL = state.cells[i].QL;
    state.cells[i - 1].QLR = state.cells[i].QL;
    state.cells[i].QG = (state.cells[i].MC - state.cells[i].Mliqini) / rhog;

    if (i < state.lastCell)
        state.cells[i + 1].QLL = state.cells[i].QL;
    state.cells[i - 1].QLR = state.cells[i].QL;

    // definicao de temperaturas maximas e minimas para uma eventual reavaliacao do metodo
    // ASTM quando ocorre mistura de fluidos e se deseja atualizar o modelo
    // de viscosidade de oleo morto se este for o ASTM que trabalha com um par de temperatura
    double temperatureLow;
    if (state.input.flashCompleto == 1 && (state.input.tabent.tmin - 0) > (*state.globals).localtiny)
        temperatureLow = state.input.tabent.tmin + 0.1;
    else
        temperatureLow = 0.;
    double temperatureHigh;
    if (state.input.flashCompleto == 1 && (70 - state.input.tabent.tmax) < (*state.globals).localtiny)
        temperatureHigh = state.input.tabent.tmax - 0.1;
    else
        temperatureHigh = 70.;
    double oilVolumeFactor;
    double waterVolumeFactor;
    double solutionGasRatio;
    // calcula os valores de RS, Bo e Ba na celula anterior a i-esima celula
    if (titulo < 1 - 1e-15) {
        solutionGasRatio = state.cells[i - 1].flui.RS(state.cells[i - 1].pres, state.cells[i - 1].temp);
        oilVolumeFactor = state.cells[i - 1].flui.BOFunc(state.cells[i - 1].pres, state.cells[i - 1].temp, solutionGasRatio);
        waterVolumeFactor = state.cells[i - 1].flui.BAFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        solutionGasRatio = solutionGasRatio * 6.29 / 35.31467;
    } else {
        oilVolumeFactor = 1;
        solutionGasRatio = 0;
        waterVolumeFactor = 0.;
    }
    // BSW in-situ da celula anterior, na marcha, a i-esima celula
    state.cells[i - 1].FW = state.cells[i - 1].flui.BSW * waterVolumeFactor / (oilVolumeFactor + waterVolumeFactor * state.cells[i - 1].flui.BSW - state.cells[i - 1].flui.BSW * oilVolumeFactor);
    state.cells[i - 1].FWini = state.cells[i - 1].FW;
    // Definicao das fracoes volumetricas:
    applyReverseSteadyPhaseFractions(state, i, rhol, rhog);
    if (fabs(state.cells[i].QL) > 1e-15 && fabs(state.cells[i].MComp) > 1e-15) {
        double dxmed = 0.5 * (state.cells[i - 1].dx + state.cells[i].dx);
        double hol0 = 1. - state.cells[i - 1].alf;
        double hol1 = 1. - state.cells[i].alf;
        state.cells[i].fluicol.TR = (state.cells[i - 1].fluicol.TR * state.cells[i - 1].MComp + state.cells[i - 1].fontemassCR * residenceTimeSource) / state.cells[i].MComp;
        state.cells[i].fluicol.TR = state.cells[i].fluicol.TR + (0.5 * state.cells[i - 1].duto.area * state.cells[i - 1].dx * hol0 +
                                                       0.5 * state.cells[i].duto.area * state.cells[i].dx * hol1) /
                                                          state.cells[i].QL;
    } else
        state.cells[i].fluicol.TR = state.cells[i - 1].fluicol.TR;
}

namespace {

/// Two-phase branch of the compositional phase-regime chain.
///
/// Split out of applyCompositionalPhaseRegime only so that neither is over 200
/// lines. It is one arm of a four-way chain and has no other caller.
void applyCompositionalTwoPhaseRegime(const SteadyStateState &state, int i, double rhol, double rhog) {
    double c0 = 1.;
    double ud = 0.;
    if (fabs(state.cells[i].QL) > (*state.globals).localtiny * 1e-5) {
        if (state.steadyIteration == 0) { // primeira estimativa, primeira iteracao
            // utiliza-se a fracao de vazio sem escorregamento, pois a propria correlacao para se obter a
            // fracao de vazio depende do valor da fracao de vazio
            if ((fabs(state.cells[i].QG) + fabs(state.cells[i].QL)) > (*state.globals).localtiny) {
                if (state.convergenceMonitor > 0.01)
                    state.cells[i].alf = fabs(state.cells[i].QG) /
                                    (fabs(state.cells[i].QG) + fabs(state.cells[i].QL));
                state.updaters.steadyDriftClosure(i, c0, ud);
            } else
                state.cells[i].alf = 0.;
            state.cells[i].alfini = state.cells[i].alf;
            state.cells[i - 1].alfR = state.cells[i].alf;
            state.cells[i - 1].alfRini = state.cells[i].alf;
            if (i < state.lastCell)
                state.cells[i + 1].alfL = state.cells[i].alf;
            if (i < state.lastCell)
                state.cells[i + 1].alfLini = state.cells[i].alf;
            state.cells[i].alfPigD = state.cells[i].alf;
            state.cells[i].alfPigDini = state.cells[i].alf;
            state.cells[i].alfPigE = state.cells[i].alf;
            state.cells[i].alfPigEini = state.cells[i].alf;
        }
        if (state.input.tipoModeloDrift == 1) {
            if (fabs(rhog) / rhol > 0.9) {
                c0 = 1.;
                ud = 0.;
            }
            // para o caso permanente, a fracao de vazio Ã© obtida a partir das relacoes de escorregamento
            // portanto, e neste ponto que se obtem Co e Ud:
            else if (fabs(state.cells[i].QG) > (*state.globals).localtiny * 1e-5 && fabs(state.cells[i].QL) > (*state.globals).localtiny * 1e-5)
                state.updaters.steadyDriftClosure(i, c0, ud);
            state.cells[i].c0 = c0;
            state.cells[i].ud = ud;
            double area = state.cells[i].duto.area;
            if (fabs(state.cells[i].QG + state.cells[i].QL) > (*state.globals).localtiny * 1e-5) {
                // alfa com escorregamento:
                state.cells[i].alf = state.cells[i].QG / (c0 * (state.cells[i].QG + state.cells[i].QL) + ud * area);
                double alfHomo = state.cells[i].QG / (state.cells[i].QG + state.cells[i].QL);
                if (state.cells[i].alf > 1. - 1e-15 || state.cells[i].alf < 1e-15)
                    state.cells[i].alf = alfHomo;
            } else
                state.cells[i].alf = 0.;
            if (state.cells[i].alf > (1 - (*state.globals).localtiny) && fabs(state.cells[i].QG + state.cells[i].QL) > (*state.globals).localtiny * 1e-5)
                state.cells[i].alf = state.cells[i].QG / (state.cells[i].QG + state.cells[i].QL);
            else if (state.cells[i].alf > (1 - (*state.globals).localtiny))
                state.cells[i].alf = 1.;
        } else {
            c0 = 1.;
            ud = 0.;
            state.cells[i].alf = 1.;
        }
    } else {
        double holdup;
        double frictionGrad;
        double gravityGrad;
        double totalGrad;
        double reynolds;
        unsigned char flowType;
        char *errorMsg;
        unsigned char errorFlag;
        executarCorrelacao(state.cells, i, 0, state.input.AceleraConvergPerm,
                           state.cells[i - 1].correlacaoMR2,
                           holdup, frictionGrad, gravityGrad, totalGrad,
                           reynolds, flowType);
        if(flowType==1 || flowType==2)state.cells[i].arranjo=0;
        if(flowType==3)state.cells[i].arranjo=1;
        if(flowType==4)state.cells[i].arranjo=2;
        if(flowType==6)state.cells[i].arranjo=-1;
        if(flowType==5)state.cells[i].arranjo=-2;
        state.cells[i].alf = 1. - holdup;
    }
    if (state.cells[i].alf < 0.)
        state.cells[i].alf = 0.;
    else if (state.cells[i].alf > 1.)
        state.cells[i].alf = 1.;
    // atualizacoes dos valores das fracoes volumetricas da celula i armazendadas em
    // outras celulas, e inclusiove armazendo os valores para "tempo anterior", que nao
    // sao relevantes para o problema permanente mas importantes se o resultado permanente
    // der partida na solucao transiente:
    state.cells[i].alfini = state.cells[i].alf;
    state.cells[i - 1].alfR = state.cells[i].alf;
    state.cells[i - 1].alfRini = state.cells[i].alf;
    if (i < state.lastCell)
        state.cells[i + 1].alfL = state.cells[i].alf;
    if (i < state.lastCell)
        state.cells[i + 1].alfLini = state.cells[i].alf;
    state.cells[i].alfPigD = state.cells[i].alf;
    state.cells[i].alfPigDini = state.cells[i].alf;
    state.cells[i].alfPigE = state.cells[i].alf;
    state.cells[i].alfPigEini = state.cells[i].alf;
    if (i < state.lastCell) {
        state.cells[i + 1].betL = state.cells[i].bet;
        state.cells[i + 1].betLini = state.cells[i].bet;
        state.cells[i + 1].betL = state.cells[i].bet;
        state.cells[i + 1].betLini = state.cells[i].bet;
        state.cells[i + 1].betLI = state.cells[i].bet;
    }
    state.cells[i].betini = state.cells[i].bet;
    state.cells[i - 1].betR = state.cells[i].bet;
    state.cells[i - 1].betRini = state.cells[i].bet;
    state.cells[i].betPigD = state.cells[i].bet;
    state.cells[i].betPigDini = state.cells[i].bet;
    state.cells[i].betPigE = state.cells[i].bet;
    state.cells[i].betPigEini = state.cells[i].bet;
    state.cells[i].betI = state.cells[i].bet;
    state.cells[i - 1].betRI = state.cells[i].bet;

}
}  // namespace

namespace {

/// Resolves the phase fractions of cell i for the compositional mass march.
///
/// Same four regimes as applyReverseSteadyPhaseFractions and NOT the same
/// function: this one compares against localtiny * 1e-5 where the reverse
/// variants compare against localtiny. Five orders of magnitude, with nothing
/// in the code saying why. Sharing the two would mean choosing one of the
/// thresholds, and that is a decision, not a move.
void applyCompositionalPhaseRegime(const SteadyStateState &state, int i, double rhol, double rhog) {
    if (fabs(state.cells[i].QG + state.cells[i].QL) < (*state.globals).localtiny * 1e-5) {
        if (state.input.tipoFluido == 1) {
            state.cells[i].alf = 0.;
        } else {
            state.cells[i].alf = 1.;
        }
        state.cells[i].alfini = state.cells[i].alf;
        state.cells[i - 1].alfR = state.cells[i].alf;
        state.cells[i - 1].alfRini = state.cells[i].alf;
        if (i < state.lastCell)
            state.cells[i + 1].alfL = state.cells[i].alf;
        if (i < state.lastCell)
            state.cells[i + 1].alfLini = state.cells[i].alf;
        state.cells[i].alfPigD = state.cells[i].alf;
        state.cells[i].alfPigDini = state.cells[i].alf;
        state.cells[i].alfPigE = state.cells[i].alf;
        state.cells[i].alfPigEini = state.cells[i].alf;
        if (i < state.lastCell) {
            state.cells[i + 1].betL = state.cells[i].bet;
            state.cells[i + 1].betLini = state.cells[i].bet;
            state.cells[i + 1].betL = state.cells[i].bet;
            state.cells[i + 1].betLini = state.cells[i].bet;
            state.cells[i + 1].betLI = state.cells[i].bet;
        }
        state.cells[i].betini = state.cells[i].bet;
        state.cells[i - 1].betR = state.cells[i].bet;
        state.cells[i - 1].betRini = state.cells[i].bet;
        state.cells[i].betPigD = state.cells[i].bet;
        state.cells[i].betPigDini = state.cells[i].bet;
        state.cells[i].betPigE = state.cells[i].bet;
        state.cells[i].betPigEini = state.cells[i].bet;
        state.cells[i].betI = state.cells[i].bet;
        state.cells[i - 1].betRI = state.cells[i].bet;
    } else if (fabs(state.cells[i].QG) < (*state.globals).localtiny * 1e-5) { // caso so exista liquido:
        state.cells[i].alf = 0.;
        state.cells[i].alfini = state.cells[i].alf;
        state.cells[i - 1].alfR = state.cells[i].alf;
        state.cells[i - 1].alfRini = state.cells[i].alf;
        if (i < state.lastCell)
            state.cells[i + 1].alfL = state.cells[i].alf;
        if (i < state.lastCell)
            state.cells[i + 1].alfLini = state.cells[i].alf;
        state.cells[i].alfPigD = state.cells[i].alf;
        state.cells[i].alfPigDini = state.cells[i].alf;
        state.cells[i].alfPigE = state.cells[i].alf;
        state.cells[i].alfPigEini = state.cells[i].alf;
        if (i < state.lastCell) {
            state.cells[i + 1].betL = state.cells[i].bet;
            state.cells[i + 1].betLini = state.cells[i].bet;
            state.cells[i + 1].betL = state.cells[i].bet;
            state.cells[i + 1].betLini = state.cells[i].bet;
            state.cells[i + 1].betLI = state.cells[i].bet;
        }
        state.cells[i].betini = state.cells[i].bet;
        state.cells[i - 1].betR = state.cells[i].bet;
        state.cells[i - 1].betRini = state.cells[i].bet;
        state.cells[i].betPigD = state.cells[i].bet;
        state.cells[i].betPigDini = state.cells[i].bet;
        state.cells[i].betPigE = state.cells[i].bet;
        state.cells[i].betPigEini = state.cells[i].bet;
        state.cells[i].betI = state.cells[i].bet;
        state.cells[i - 1].betRI = state.cells[i].bet;
    } else if (fabs(state.cells[i].QL) < (*state.globals).localtiny * 1e-5) { // caso so exista gas:
        state.cells[i].alf = 1.;
        state.cells[i].alfini = state.cells[i].alf;
        state.cells[i - 1].alfR = state.cells[i].alf;
        state.cells[i - 1].alfRini = state.cells[i].alf;
        if (i < state.lastCell)
            state.cells[i + 1].alfL = state.cells[i].alf;
        if (i < state.lastCell)
            state.cells[i + 1].alfLini = state.cells[i].alf;
        state.cells[i].alfPigD = state.cells[i].alf;
        state.cells[i].alfPigDini = state.cells[i].alf;
        state.cells[i].alfPigE = state.cells[i].alf;
        state.cells[i].alfPigEini = state.cells[i].alf;
        state.cells[i].c0 = 1.;
        state.cells[i].ud = 0.;
        if (i < state.lastCell) {
            state.cells[i + 1].betL = state.cells[i].bet;
            state.cells[i + 1].betLini = state.cells[i].bet;
            state.cells[i + 1].betL = state.cells[i].bet;
            state.cells[i + 1].betLini = state.cells[i].bet;
            state.cells[i + 1].betLI = state.cells[i].bet;
        }
        state.cells[i].betini = state.cells[i].bet;
        state.cells[i - 1].betR = state.cells[i].bet;
        state.cells[i - 1].betRini = state.cells[i].bet;
        state.cells[i].betPigD = state.cells[i].bet;
        state.cells[i].betPigDini = state.cells[i].bet;
        state.cells[i].betPigE = state.cells[i].bet;
        state.cells[i].betPigEini = state.cells[i].bet;
        state.cells[i].betI = state.cells[i].bet;
        state.cells[i - 1].betRI = state.cells[i].bet;
    } else { // caso bifasico
        applyCompositionalTwoPhaseRegime(state, i, rhol, rhog);
    }
}
}  // namespace

namespace {

/// Reads the black-oil reference properties of the accessory feeding cell i - 1.
///
/// Seven arms on acsr.tipo, each writing the same seven outputs from a
/// different source object. The wide out-parameter list is the shape T085 gave
/// the helpers it carved out of the base variant; it is kept rather than
/// improved so the two families still read alike.
void readCompositionalSourceProperties(const SteadyStateState &state, int i, double &titF, ProFlu &fluF,
                                       double &oilVolumeFactorSource, double &waterVolumeFactorSource, double &waterCutSource, double &rhoOF, double &rhoWF,
                                       double &residenceTimeSource) {
    if (state.cells[i - 1].acsr.tipo == 1) {
        if (i > 0 && (state.cells[i - 1].flui.dCalculatedBeta > 0. && state.cells[i - 1].flui.dCalculatedBeta < 1.))
            state.cells[i - 1].acsr.injg.FluidoPro.atualizaPropComp(state.cells[i - 1].pres, state.cells[i - 1].temp, state.cells[i - 1].flui.dCalculatedBeta,
                                                               state.cells[i - 1].flui.oCalculatedLiqComposition, state.cells[i - 1].flui.oCalculatedVapComposition, state.input.pocinjec);
        else
            state.cells[i - 1].acsr.injg.FluidoPro.atualizaPropComp(state.cells[i - 1].pres, state.cells[i - 1].temp, -1, NULL, NULL, state.input.pocinjec);
        fluF = state.cells[i - 1].acsr.injg.FluidoPro;
        waterCutSource = 0.;
        titF = 1.;
        residenceTimeSource = state.cells[i - 1].acsr.injg.fluidocol.TR;
    } else if (state.cells[i - 1].acsr.tipo == 2) {
        if (i > 0 && (state.cells[i - 1].flui.dCalculatedBeta > 0. && state.cells[i - 1].flui.dCalculatedBeta < 1.))
            state.cells[i - 1].acsr.injl.FluidoPro.atualizaPropComp(state.cells[i - 1].pres, state.cells[i - 1].temp, state.cells[i - 1].flui.dCalculatedBeta,
                                                               state.cells[i - 1].flui.oCalculatedLiqComposition, state.cells[i - 1].flui.oCalculatedVapComposition, state.input.pocinjec);
        else
            state.cells[i - 1].acsr.injl.FluidoPro.atualizaPropComp(state.cells[i - 1].pres, state.cells[i - 1].temp, -1, NULL, NULL, state.input.pocinjec);
        fluF = state.cells[i - 1].acsr.injl.FluidoPro;
        oilVolumeFactorSource = fluF.BOFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        waterVolumeFactorSource = fluF.BAFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        waterCutSource = fluF.BSW * waterVolumeFactorSource / (oilVolumeFactorSource + waterVolumeFactorSource * fluF.BSW - fluF.BSW * oilVolumeFactorSource);
        rhoOF = fluF.MasEspoleo(state.cells[i - 1].pres, state.cells[i - 1].temp);
        rhoWF = fluF.MasEspAgua(state.cells[i - 1].pres, state.cells[i - 1].temp);
        titF = (1 - waterCutSource) * rhoOF / ((1 - waterCutSource) * rhoOF + waterCutSource * rhoWF);
        residenceTimeSource = state.cells[i - 1].acsr.injl.fluidocol.TR;
    } else if (state.cells[i - 1].acsr.tipo == 3) {
        if (i > 0 && (state.cells[i - 1].flui.dCalculatedBeta > 0. && state.cells[i - 1].flui.dCalculatedBeta < 1.))
            state.cells[i - 1].acsr.ipr.FluidoPro.atualizaPropComp(state.cells[i - 1].pres, state.cells[i - 1].temp, state.cells[i - 1].flui.dCalculatedBeta,
                                                              state.cells[i - 1].flui.oCalculatedLiqComposition, state.cells[i - 1].flui.oCalculatedVapComposition, state.input.pocinjec);
        else
            state.cells[i - 1].acsr.ipr.FluidoPro.atualizaPropComp(state.cells[i - 1].pres, state.cells[i - 1].temp, -1, NULL, NULL, state.input.pocinjec);
        fluF = state.cells[i - 1].acsr.ipr.FluidoPro;
        oilVolumeFactorSource = fluF.BOFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        waterVolumeFactorSource = fluF.BAFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        waterCutSource = fluF.BSW * waterVolumeFactorSource / (oilVolumeFactorSource + waterVolumeFactorSource * fluF.BSW - fluF.BSW * oilVolumeFactorSource);
        rhoOF = fluF.MasEspoleo(state.cells[i - 1].pres, state.cells[i - 1].temp);
        rhoWF = fluF.MasEspAgua(state.cells[i - 1].pres, state.cells[i - 1].temp);
        titF = (1 - waterCutSource) * rhoOF / ((1 - waterCutSource) * rhoOF + waterCutSource * rhoWF);
        residenceTimeSource = 0.;
    } else if (state.cells[i - 1].acsr.tipo == 10) {
        if (i > 0 && (state.cells[i - 1].flui.dCalculatedBeta > 0. && state.cells[i - 1].flui.dCalculatedBeta < 1.))
            state.cells[i - 1].acsr.injm.FluidoPro.atualizaPropComp(state.cells[i - 1].pres, state.cells[i - 1].temp, state.cells[i - 1].flui.dCalculatedBeta,
                                                               state.cells[i - 1].flui.oCalculatedLiqComposition, state.cells[i - 1].flui.oCalculatedVapComposition, state.input.pocinjec);
        else
            state.cells[i - 1].acsr.injm.FluidoPro.atualizaPropComp(state.cells[i - 1].pres, state.cells[i - 1].temp, -1, NULL, NULL, state.input.pocinjec);
        fluF = state.cells[i - 1].acsr.injm.FluidoPro;
        oilVolumeFactorSource = fluF.BOFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        waterVolumeFactorSource = fluF.BAFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        waterCutSource = fluF.BSW * waterVolumeFactorSource / (oilVolumeFactorSource + waterVolumeFactorSource * fluF.BSW - fluF.BSW * oilVolumeFactorSource);
        rhoOF = fluF.MasEspoleo(state.cells[i - 1].pres, state.cells[i - 1].temp);
        rhoWF = fluF.MasEspAgua(state.cells[i - 1].pres, state.cells[i - 1].temp);
        titF = (1 - waterCutSource) * rhoOF / ((1 - waterCutSource) * rhoOF + waterCutSource * rhoWF);
        residenceTimeSource = state.cells[i - 1].acsr.injm.fluidocol.TR;
    } else if (state.cells[i - 1].acsr.tipo == 9 && state.cells[i - 1].acsr.fontechk.abertura > 1e-6) {

        if (i > 0 && (state.cells[i - 1].flui.dCalculatedBeta > 0. && state.cells[i - 1].flui.dCalculatedBeta < 1.))
            state.cells[i - 1].acsr.fontechk.fluidoP.atualizaPropComp(state.cells[i - 1].pres, state.cells[i - 1].temp, state.cells[i - 1].flui.dCalculatedBeta,
                                                                 state.cells[i - 1].flui.oCalculatedLiqComposition, state.cells[i - 1].flui.oCalculatedVapComposition, state.input.pocinjec);
        else
            state.cells[i - 1].acsr.fontechk.fluidoP.atualizaPropComp(state.cells[i - 1].pres, state.cells[i - 1].temp, -1, NULL, NULL, state.input.pocinjec);
        if (i > 0 && (state.cells[i - 1].flui.dCalculatedBeta > 0. && state.cells[i - 1].flui.dCalculatedBeta < 1.))
            state.cells[i - 1].acsr.fontechk.fluidoPamb.atualizaPropComp(state.cells[i - 1].pres, state.cells[i - 1].temp, state.cells[i - 1].flui.dCalculatedBeta,
                                                                    state.cells[i - 1].flui.oCalculatedLiqComposition, state.cells[i - 1].flui.oCalculatedVapComposition, state.input.pocinjec);
        else
            state.cells[i - 1].acsr.fontechk.fluidoPamb.atualizaPropComp(state.cells[i - 1].pres, state.cells[i - 1].temp, -1, NULL, NULL, state.input.pocinjec);
        if (state.cells[i - 1].acsr.fontechk.presT > state.cells[i - 1].acsr.fontechk.pamb) {
            fluF = state.cells[i - 1].acsr.fontechk.fluidoP;
        } else {
            fluF = state.cells[i - 1].acsr.fontechk.fluidoPamb;
        }
        oilVolumeFactorSource = fluF.BOFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        waterVolumeFactorSource = fluF.BAFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        waterCutSource = fluF.BSW * waterVolumeFactorSource / (oilVolumeFactorSource + waterVolumeFactorSource * fluF.BSW - fluF.BSW * oilVolumeFactorSource);
        rhoOF = fluF.MasEspoleo(state.cells[i - 1].pres, state.cells[i - 1].temp);
        rhoWF = fluF.MasEspAgua(state.cells[i - 1].pres, state.cells[i - 1].temp);
        titF = (1 - waterCutSource) * rhoOF / ((1 - waterCutSource) * rhoOF + waterCutSource * rhoWF);
    } else if (state.cells[i - 1].acsr.tipo == 15) {
        double tRes = state.cells[i - 1].acsr.radialPoro.tRes;
        if (i > 0 && (state.cells[i - 1].flui.dCalculatedBeta > 0. && state.cells[i - 1].flui.dCalculatedBeta < 1.))
            state.cells[i - 1].acsr.radialPoro.flup.atualizaPropComp(state.cells[i - 1].pres, tRes,
                                                                state.cells[i - 1].flui.dCalculatedBeta,
                                                                state.cells[i - 1].flui.oCalculatedLiqComposition, state.cells[i - 1].flui.oCalculatedVapComposition, state.input.pocinjec);
        else
            state.cells[i - 1].acsr.radialPoro.flup.atualizaPropComp(state.cells[i - 1].pres, tRes, -1, NULL, NULL, state.input.pocinjec);
        fluF = state.cells[i - 1].acsr.radialPoro.flup;
        oilVolumeFactorSource = fluF.BOFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        waterVolumeFactorSource = fluF.BAFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        waterCutSource = fluF.BSW * waterVolumeFactorSource / (oilVolumeFactorSource + waterVolumeFactorSource * fluF.BSW - fluF.BSW * oilVolumeFactorSource);
        rhoOF = fluF.MasEspoleo(state.cells[i - 1].pres, state.cells[i - 1].temp);
        rhoWF = fluF.MasEspAgua(state.cells[i - 1].pres, state.cells[i - 1].temp);
        titF = (1 - waterCutSource) * rhoOF / ((1 - waterCutSource) * rhoOF + waterCutSource * rhoWF);
        residenceTimeSource = 0.;
    } else if (state.cells[i - 1].acsr.tipo == 16) {
        double tRes = state.cells[i - 1].acsr.poroso2D.dados.transfer.tRes;
        if (i > 0 && (state.cells[i - 1].flui.dCalculatedBeta > 0. && state.cells[i - 1].flui.dCalculatedBeta < 1.))
            state.cells[i - 1].acsr.poroso2D.dados.flup.atualizaPropComp(state.cells[i - 1].pres, tRes,
                                                                    state.cells[i - 1].flui.dCalculatedBeta,
                                                                    state.cells[i - 1].flui.oCalculatedLiqComposition, state.cells[i - 1].flui.oCalculatedVapComposition, state.input.pocinjec);
        else
            state.cells[i - 1].acsr.poroso2D.dados.flup.atualizaPropComp(state.cells[i - 1].pres, tRes, -1, NULL, NULL, state.input.pocinjec);
        fluF = state.cells[i - 1].acsr.poroso2D.dados.flup;
        oilVolumeFactorSource = fluF.BOFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        waterVolumeFactorSource = fluF.BAFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        waterCutSource = fluF.BSW * waterVolumeFactorSource / (oilVolumeFactorSource + waterVolumeFactorSource * fluF.BSW - fluF.BSW * oilVolumeFactorSource);
        rhoOF = fluF.MasEspoleo(state.cells[i - 1].pres, state.cells[i - 1].temp);
        rhoWF = fluF.MasEspAgua(state.cells[i - 1].pres, state.cells[i - 1].temp);
        titF = (1 - waterCutSource) * rhoOF / ((1 - waterCutSource) * rhoOF + waterCutSource * rhoWF);
        residenceTimeSource = 0.;
    }
}

/// Cell i-1 holds a source: the stream it passes to cell i mixes its own flow
/// with the source's, weighted by the standard oil and water rates before and
/// from the source, and the complementary-liquid fraction is re-evaluated.
/// Cut from advanceCompositionalSteadyMass (SC-004).
void mixUpstreamSourceIntoCell(const SteadyStateState &state, int i, double tmed, double &titV, double solutionGasRatio, double waterVolumeFactor, double oilVolumeFactor, double temperatureHigh, double temperatureLow, double vazMasGasL, double &vazMasLiqL, double fonteMasGasL, double &fonteMasLiqL, double mComp, double &mHidro, ProFlu &fluF, double titF, int mudaRGO) {
    double waterCutCarried = state.cells[i - 1].FW;
    double rhoOV = state.cells[i - 1].flui.MasEspoleo(state.cells[i - 1].pres, state.cells[i - 1].temp);
    double rhoWV = state.cells[i - 1].flui.MasEspAgua(state.cells[i - 1].pres, state.cells[i - 1].temp);
    titV = (1 - waterCutCarried) * rhoOV / ((1 - waterCutCarried) * rhoOV + waterCutCarried * rhoWV);
    vazMasLiqL *= titV;
    fonteMasLiqL *= titF;
    if (state.input.tabelaDinamica == 0) {
        double pesoMolV = 0;
        double pesoMolF = 0;
        for (int j = 0; j < state.cells[i].flui.npseudo; j++) {
            pesoMolV += state.cells[i - 1].flui.masMol[j] * state.cells[i - 1].flui.fracMol[j];
            pesoMolF += fluF.masMol[j] * fluF.fracMol[j];
        }
        double vazMolV = (vazMasLiqL + vazMasGasL) / pesoMolV;
        double vazMolF = (fonteMasLiqL + fonteMasGasL) / pesoMolF;
        double razMolV = 1.;
        if (fabs(vazMolV + vazMolF) > 1e-15)
            razMolV = vazMolV / (vazMolV + vazMolF);
        if (vazMolF > 0.) {
            for (int j = 0; j < state.cells[i].flui.npseudo; j++) {
                state.cells[i].flui.fracMol[j] = razMolV * state.cells[i - 1].flui.fracMol[j] +
                                            (1. - razMolV) * fluF.fracMol[j];
            }
        } else {
            for (int j = 0; j < state.cells[i].flui.npseudo; j++)
                state.cells[i].flui.fracMol[j] = state.cells[i - 1].flui.fracMol[j];
        }
        state.cells[i].flui.Pmol = 0.;
        for (int j = 0; j < state.cells[i].flui.npseudo; j++)
            state.cells[i].flui.Pmol += state.cells[i].flui.fracMol[j] * state.cells[i].flui.masMol[j];

        state.cells[i].flui.atualizaPropCompStandard();
    }

    double apiGravity = state.cells[i - 1].flui.API;
    double rholis = state.cells[i - 1].flui.MasEspLiq(state.cells[i - 1].pres, state.cells[i - 1].temp);
    double rholisF;
    double boinjl;
    double bainjl;
    double fwinjl;
    if (fabs(fonteMasLiqL + fonteMasGasL) > 1e-15) {
        rholisF = fluF.MasEspLiq(state.cells[i - 1].pres, state.cells[i - 1].temp);
        boinjl = fluF.BOFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        bainjl = fluF.BAFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        fwinjl = fluF.BSW * bainjl / (boinjl + bainjl * fluF.BSW - fluF.BSW * boinjl);
    } else {
        rholisF = rholis;
        boinjl = oilVolumeFactor;
        bainjl = waterVolumeFactor;
        fwinjl = state.cells[i - 1].FW;
    }
    // vazao de oleo sytandard antes da fonte:
    double qostd1;
    if (state.cells[i - 1].flui.dStockTankVaporMassFraction < 1. - 1e-15)
        qostd1 = state.cells[i - 1].Mliqini * (1. - state.cells[i - 1].FW) * (1. - state.cells[i - 1].bet) / (oilVolumeFactor * rholis);
    else
        qostd1 = 0.;

    // vazao de oleo standard da fonte:
    double qostd2;
    if (fluF.dStockTankVaporMassFraction < 1. - 1e-15)
        qostd2 = state.cells[i - 1].fontemassLR * (1. - fwinjl) / (boinjl * rholisF);
    else
        qostd2 = 0.;

    double liquidHoldup = 1 - state.cells[i - 1].alf;
    double completionFraction = state.cells[i - 1].bet;
    double bsw = state.cells[i - 1].FW;
    double rhog = state.cells[i - 1].flui.MasEspGas(state.cells[i - 1].pres, state.cells[i - 1].temp);
    double rhogST = state.cells[i - 1].flui.Deng * 1.225;
    // calculo do volume de leve na celula i, seguindo o equacionamento mostrado em relatorio,
    // nÃ£o Ã© relevante para o permanente, mas deve ser calculado,
    // pois Ã© utilizado como entrada no transiente
    state.cells[i - 1].VolLeveST = (((1 - liquidHoldup) * rhog / rhogST) + liquidHoldup * (1 - completionFraction) * (1. - bsw) * solutionGasRatio / oilVolumeFactor);
    if (state.cells[i - 1].VolLeveST < 1e-15)
        state.cells[i - 1].VolLeveST = 0.;
    if (state.input.trackRGO == -1) {
        // esta chave nÃ£o Ã© utilizada, Ã© mantida aqui como reserva, atualmente este
        // calculo nÃ£o esta funcionando, lembrando que arq.trackRGO assume apenas 2 valores, 0 ou 1
        if (liquidHoldup > (*state.globals).localtiny && completionFraction < (1. - (*state.globals).localtiny) && bsw < (1. - (*state.globals).localtiny) && mudaRGO == 1)
            state.cells[i].flui.RGO = state.cells[i - 1].VolLeveST * oilVolumeFactor / (liquidHoldup * (1 - completionFraction) * (1 - bsw));
        if (state.cells[i].flui.RGO > (*state.globals).RGOMax && mudaRGO == 1 && state.cells[i].flui.RGO < 1e6)
            state.cells[i].flui.RGO = (*state.globals).RGOMax;
    }

    double waterFlowRate1;
    if ((1. - state.cells[i - 1].flui.BSW) > 0)
        waterFlowRate1 = qostd1 * state.cells[i - 1].flui.BSW / (1. - state.cells[i - 1].flui.BSW);
    else
        waterFlowRate1 = state.cells[i - 1].Mliqini * state.cells[i - 1].FW * (1. - state.cells[i - 1].bet) / (1000. * state.cells[i - 1].flui.Denag);
    double waterFlowRate2;
    if ((1. - fluF.BSW) > 0)
        waterFlowRate2 = qostd2 * fluF.BSW / (1. - fluF.BSW);
    else
        waterFlowRate2 = state.cells[i - 1].fontemassLR * fwinjl / (1000. * fluF.Denag);

    if (fabs(waterFlowRate1 + waterFlowRate2 + qostd1 + qostd2) > 1e-15 && fabs(waterFlowRate2 + qostd2) > 1e-15 && state.cells[i - 1].fontemassLR > 1e-15)
        state.cells[i].flui.BSW = (waterFlowRate1 + waterFlowRate2) / (waterFlowRate1 + waterFlowRate2 + qostd1 + qostd2);
    else
        state.cells[i].flui.BSW = state.cells[i - 1].flui.BSW;

    if (fabs(waterFlowRate1 + waterFlowRate2) > 1e-15 && state.cells[i - 1].fontemassLR > 1e-15)
        state.cells[i].flui.Denag = (waterFlowRate1 * state.cells[i - 1].flui.Denag +
                                waterFlowRate2 * fluF.Denag) /
                               (waterFlowRate1 + waterFlowRate2);
    else
        state.cells[i].flui.Denag = state.cells[i - 1].flui.Denag;

    if (fabs(qostd1 + qostd2) > 1e-15 && fabs(qostd2) > 1e-15 && state.cells[i - 1].fontemassLR > 1e-15) { // vazao de liquido >0
        state.cells[i].flui.TempL = temperatureLow;
        state.cells[i].flui.TempH = temperatureHigh;
        state.cells[i].flui.LVisL = (qostd1 * state.cells[i - 1].flui.VisOM(temperatureLow) + qostd2 * fluF.VisOM(temperatureLow)) / (qostd1 + qostd2);
        state.cells[i].flui.LVisH = (qostd1 * state.cells[i - 1].flui.VisOM(temperatureHigh) + qostd2 * fluF.VisOM(temperatureHigh)) / (qostd1 + qostd2);
    } else {
        state.cells[i].flui.TempL = state.cells[i - 1].flui.TempL;
        state.cells[i].flui.TempH = state.cells[i - 1].flui.TempH;
        state.cells[i].flui.LVisL = state.cells[i - 1].flui.LVisL;
        state.cells[i].flui.LVisH = state.cells[i - 1].flui.LVisH;
    }

    if (state.input.tipoFluido == 0) {
        // reavaliacao da fracao volumetrica do liquido complementar
        // observar que a fonte de liquido pode ter uma fracao de liquido complementar
        // distinta da onservada a esquerda da celula i-1
        double oilVolumeFactorNeighbour = state.cells[i].flui.BOFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        double waterVolumeFactorNeighbour = state.cells[i].flui.BAFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        double waterCutNeighbour = state.cells[i].flui.BSW * waterVolumeFactorNeighbour /
                     (oilVolumeFactorNeighbour + waterVolumeFactorNeighbour * state.cells[i].flui.BSW - state.cells[i].flui.BSW * oilVolumeFactorNeighbour);
        double oilVolumeFactorInSitu = state.cells[i].flui.BOFunc(state.cells[i].presaux, tmed);
        double waterVolumeFactorInSitu = state.cells[i].flui.BAFunc(state.cells[i].presaux, tmed);
        double waterCutInSitu = state.cells[i].flui.BSW * waterVolumeFactorInSitu /
                     (oilVolumeFactorInSitu + waterVolumeFactorInSitu * state.cells[i].flui.BSW - state.cells[i].flui.BSW * oilVolumeFactorInSitu);
        double qlpF = state.cells[i - 1].fontemassLR / fluF.MasEspLiq(state.cells[i - 1].pres, state.cells[i - 1].temp);
        double qlcF = state.cells[i - 1].fontemassCR / state.cells[i].fluicol.MasEspFlu(state.cells[i - 1].pres, state.cells[i - 1].temp);
        double betN = (state.cells[i - 1].QL * state.cells[i - 1].bet + qlcF) /
                      (state.cells[i - 1].QL + qlpF + qlcF);
        double oilFlowRate = (state.cells[i - 1].QL + qlpF + qlcF) * (1 - waterCutNeighbour) * (1 - betN) * oilVolumeFactorInSitu / oilVolumeFactor;
        double waterFlowRate;
        if (waterCutInSitu < (1 - (*state.globals).localtiny * 1e-5))
            waterFlowRate = waterCutInSitu * oilFlowRate / (1 - waterCutInSitu);
        else
            waterFlowRate = state.cells[i - 1].QL * (1 - state.cells[i - 1].bet) + qlpF;
        double completionFlowRate = (state.cells[i - 1].QL * (state.cells[i - 1].bet) + qlcF) * state.cells[i - 1].fluicol.MasEspFlu(state.cells[i - 1].pres, state.cells[i - 1].temp) / state.cells[i].fluicol.MasEspFlu(state.cells[i].presaux, tmed);
        if (fabs(oilFlowRate + waterFlowRate + completionFlowRate) > 1e-15)
            state.cells[i].bet = fabs(completionFlowRate) / (fabs(oilFlowRate) + fabs(waterFlowRate) + fabs(completionFlowRate));
        else
            state.cells[i].bet = state.cells[i - 1].bet;
    } else {
        mHidro = state.cells[i].MC - mComp;
    }
}

/// Cell i-1 holds no source: the separator gas-oil ratio, BSW, API, gas density
/// and the rest pass unchanged into cell i, and the volume of light components is
/// refreshed, because the transient reads it.
/// Cut from advanceCompositionalSteadyMass (SC-004).
void carryUpstreamCompositionIntoCell(const SteadyStateState &state, int i, double tmed, double &titV, double solutionGasRatio, double oilVolumeFactor, double mComp, double &mHidro) {
    // neste caso, variaveis como RGO de separador, BSW, API, densidade de gas e outras nÃ£o muda, sao iguais
    // aos valores da celula i-1
    double liquidHoldup = 1 - state.cells[i - 1].alf;
    double completionFraction = state.cells[i - 1].bet;
    double bsw = state.cells[i - 1].FW;
    double rhog = state.cells[i - 1].flui.MasEspGas(state.cells[i - 1].pres, state.cells[i - 1].temp);
    double rhogST = state.cells[i - 1].flui.Deng * 1.225;
    // o valor de volume de leve Ã© atualizado neste ponto,
    // seguindo o equacionamento mostrado em relatorio, nÃ£o Ã© relevante para o permanente, mas
    // deve ser calculado, pois Ã© utilizado como entrada no transiente
    state.cells[i - 1].VolLeveST = (((1 - liquidHoldup) * rhog / rhogST) + liquidHoldup * (1 - completionFraction) * (1. - bsw) * solutionGasRatio / oilVolumeFactor);
    if (state.cells[i - 1].VolLeveST < 1e-15)
        state.cells[i - 1].VolLeveST = 0.;

    state.cells[i].flui.BSW = state.cells[i - 1].flui.BSW;

    double waterCutCarried = state.cells[i - 1].FW;
    double rhoOV = state.cells[i - 1].flui.MasEspoleo(state.cells[i - 1].pres, state.cells[i - 1].temp);
    double rhoWV = state.cells[i - 1].flui.MasEspAgua(state.cells[i - 1].pres, state.cells[i - 1].temp);
    titV = (1 - waterCutCarried) * rhoOV / ((1 - waterCutCarried) * rhoOV + waterCutCarried * rhoWV);

    for (int j = 0; j < state.cells[i].flui.npseudo; j++)
        state.cells[i].flui.fracMol[j] = state.cells[i - 1].flui.fracMol[j];
    if (i > 1) {
        state.cells[i].flui.iCalculatedStockTankThermodynamicCondition = state.cells[i - 1].flui.iCalculatedStockTankThermodynamicCondition;
        state.cells[i].flui.dStockTankVaporMassFraction = state.cells[i - 1].flui.dStockTankVaporMassFraction;
        state.cells[i].flui.dStockTankLiquidDensity = state.cells[i - 1].flui.dStockTankLiquidDensity;
        state.cells[i].flui.dStockTankVaporDensity = state.cells[i - 1].flui.dStockTankVaporDensity;

        if (state.cells[i].flui.dStockTankLiquidDensity > 0.01) {
            state.cells[i].flui.API = 141.5 / (state.cells[i].flui.dStockTankLiquidDensity / 1000.) - 131.5;
        } else
            state.cells[i].flui.API = 50;
        state.cells[i].flui.Deng = state.cells[i].flui.dStockTankVaporDensity / 1.225;
        state.cells[i].flui.RGO = state.cells[i - 1].flui.RGO;
        state.cells[i].flui.IRGO = state.cells[i - 1].flui.IRGO;
    } else {
        state.cells[i].flui.atualizaPropCompStandard();
    }

    if (state.input.tipoFluido == 0) { // reavaliacao da fracao volumetrica do liquido complementar
        // mesmo que nÃ£o tenha fonte, ela pode mudar, devido ao encolhimento do liquido produzido
        double oilVolumeFactorInSitu;
        double waterVolumeFactorInSitu;
        if (state.cells[i].flui.RGO < 1e6)
            oilVolumeFactorInSitu = state.cells[i].flui.BOFunc(state.cells[i].presaux, tmed);
        else
            oilVolumeFactorInSitu = 1.;
        waterVolumeFactorInSitu = state.cells[i].flui.BAFunc(state.cells[i].presaux, tmed);
        double waterCutInSitu = state.cells[i].flui.BSW * waterVolumeFactorInSitu / (oilVolumeFactorInSitu + waterVolumeFactorInSitu * state.cells[i].flui.BSW - state.cells[i].flui.BSW * oilVolumeFactorInSitu);
        double oilFlowRate = state.cells[i - 1].QL * (1 - state.cells[i - 1].FW) * (1 - state.cells[i - 1].bet) * oilVolumeFactorInSitu / oilVolumeFactor;
        double waterFlowRate;
        if (waterCutInSitu < (1 - (*state.globals).localtiny * 1e-5))
            waterFlowRate = waterCutInSitu * oilFlowRate / (1 - waterCutInSitu);
        else
            waterFlowRate = state.cells[i - 1].QL * (1 - state.cells[i - 1].bet);
        double completionFlowRate;
        completionFlowRate = state.cells[i - 1].QL * (state.cells[i - 1].bet) * state.cells[i - 1].fluicol.MasEspFlu(state.cells[i - 1].pres, state.cells[i - 1].temp) / state.cells[i].fluicol.MasEspFlu(state.cells[i].presaux, tmed);

        ////////////////////////esperar//////////////////////////////////////
        if (fabs(oilFlowRate + waterFlowRate + completionFlowRate) > 1e-15)
            state.cells[i].bet = fabs(completionFlowRate) / (fabs(oilFlowRate) + fabs(waterFlowRate) + fabs(completionFlowRate));
        else
            state.cells[i].bet = state.cells[i - 1].bet;
    } else {
        mHidro = state.cells[i].MC - mComp;
    }
}

}  // namespace

void advanceCompositionalSteadyMass(const SteadyStateState &state, int i) {
    int mudaRGO = 1;
    double titF = 0.;
    ProFlu fluF;

    double oilVolumeFactorSource = 1.;
    double waterVolumeFactorSource = 1.;
    double waterCutSource = 1.;
    double rhoOF = 900.;
    double rhoWF = 1000.;

    double mHidro = 0.;
    double mComp;
    double residenceTimeSource = 0.;

    readCompositionalSourceProperties(state, i, titF, fluF, oilVolumeFactorSource, waterVolumeFactorSource, waterCutSource, rhoOF, rhoWF, residenceTimeSource);

    state.updaters.updateSource(i - 1); // metodo que verifica se existe uma fonte na celula e calcula o valor das
    // vazoes massicas de liquido produzido (oleo+agua), gas e liquido complementar
    // relacao entre fonte a esquerda e a direita de uma celula:
    state.cells[i].fontemassLL = state.cells[i - 1].fontemassLR;
    state.cells[i].fontemassCL = state.cells[i - 1].fontemassCR;
    state.cells[i].fontemassGL = state.cells[i - 1].fontemassGR;
    mComp = state.cells[i].MComp = state.cells[i - 1].MComp + state.cells[i - 1].fontemassCR;
    state.cells[i].MC = state.cells[i - 1].MC + state.cells[i - 1].fontemassLR + state.cells[i - 1].fontemassCR + state.cells[i - 1].fontemassGR;
    if (i == 0) {
        if (state.cells[0].acsr.tipo == 15)
            state.cells[0].flui.BSW = state.cells[0].acsr.radialPoro.BSW;
        else if (state.cells[0].acsr.tipo == 16)
            state.cells[0].flui.BSW = state.cells[0].acsr.poroso2D.dados.transfer.BSW;
    }

    double fonteMasLiqL = state.cells[i].fontemassLL;
    double fonteMasGasL = state.cells[i].fontemassGL;
    double vazMasLiqL = state.cells[i].MliqiniL;
    double vazMasGasL = state.cells[i].ML - state.cells[i].MliqiniL;

    // definicao de temperaturas maximas e minimas para uma eventual reavaliacao do metodo
    // ASTM quando ocorre mistura de fluidos e se deseja atualizar o modelo
    // de viscosidade de oleo morto se este for o ASTM que trabalha com um par de temperatura
    double temperatureLow = 0;
    double temperatureHigh = 70.;

    double oilVolumeFactor;
    double waterVolumeFactor;
    double solutionGasRatio;
    double titV = 0.;
    // calcula os valores de RS, Bo e Ba na celula anterior a i-esima celula
    if (state.cells[i - 1].flui.RGO < 1e7) {
        solutionGasRatio = state.cells[i - 1].flui.RS(state.cells[i - 1].pres, state.cells[i - 1].temp);
        oilVolumeFactor = state.cells[i - 1].flui.BOFunc(state.cells[i - 1].pres, state.cells[i - 1].temp, solutionGasRatio);
        waterVolumeFactor = state.cells[i - 1].flui.BAFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        solutionGasRatio = solutionGasRatio * 6.29 / 35.31467;
    } else {
        oilVolumeFactor = 1;
        solutionGasRatio = 0;
        waterVolumeFactor = 0.;
    }
    // BSW in-situ da celula anterior, na marcha, a i-esima celula
    state.cells[i - 1].FW = state.cells[i - 1].flui.BSW * waterVolumeFactor / (oilVolumeFactor + waterVolumeFactor * state.cells[i - 1].flui.BSW - state.cells[i - 1].flui.BSW * oilVolumeFactor);
    state.cells[i - 1].FWini = state.cells[i - 1].FW;
    double tmed;
    // temperatura na fronteira entre a celula i-1 e a celula i, da segunda iteracao em diante, pode-se
    // usar a temperatura da celula i, pois jÃ¡ a tem calculada, mas isto pode ser um complicador de
    // convergencia, mais seguro manter o criterio em todas as iteracoes, neste caso, a temperatura na
    // fronteira e admitida = a temperatura da celula i-1
    if (state.steadyIteration != 0 && state.input.AceleraConvergPerm == 0)
        tmed = (state.cells[i].dx * state.cells[i].temp + state.cells[i].dxL * state.cells[i].tempL) / (state.cells[i].dx + state.cells[i].dxL);
    else
        tmed = state.cells[i - 1].temp;
    // primeiro teste: nÃ£o hÃ¡ fontes na celula i-1:
    if (state.cells[i - 1].acsr.tipo != 1 && state.cells[i - 1].acsr.tipo != 2 && state.cells[i - 1].acsr.tipo != 3 && state.cells[i - 1].acsr.tipo != 10 && (state.cells[i - 1].acsr.tipo != 9 || (state.cells[i - 1].acsr.tipo == 9 && state.cells[i - 1].acsr.fontechk.abertura <= 1e-6)) &&
        state.cells[i - 1].acsr.tipo != 15 && state.cells[i - 1].acsr.tipo != 16) {
        carryUpstreamCompositionIntoCell(state, i, tmed, titV, solutionGasRatio, oilVolumeFactor, mComp, mHidro);
    } else {
        mixUpstreamSourceIntoCell(state, i, tmed, titV, solutionGasRatio, waterVolumeFactor, oilVolumeFactor, temperatureHigh, temperatureLow, vazMasGasL, vazMasLiqL, fonteMasGasL, fonteMasLiqL, mComp, mHidro, fluF, titF, mudaRGO);
    }
    // atualizaÃ§Ã£o da vazao massica da mistura:
    state.cells[i].MC = state.cells[i - 1].MC + state.cells[i - 1].fontemassCR + state.cells[i - 1].fontemassLR + state.cells[i - 1].fontemassGR;
    state.cells[i - 1].MR = state.cells[i].MC;
    if (i < state.lastCell)
        state.cells[i + 1].ML = state.cells[i].MC;
    state.cells[i - 1].MRini = state.cells[i - 1].MR;

    double pmed = state.cells[i].presaux + state.cells[i - 1].dpB / 98066.5;
    double rhog;
    double rhol;
    // calculo da vazao massica de liquido e das vazoes volumetricas:

    refreshCompositionalFlashFromNeighbour(state, i, pmed, tmed);
    double titulo = state.cells[i].flui.dVaporMassFraction;

    double betI;
    double waterCutInSitu;
    ProFlu fluI;
    ProFluCol fluCI;
    if (state.cells[i].MC >= 0 && i > 0) {
        betI = state.cells[i].betL;
        fluI = state.cells[i - 1].flui;
        fluCI = state.cells[i - 1].fluicol;
    } else {
        betI = state.cells[i].bet;
        fluI = state.cells[i].flui;
        fluCI = state.cells[i].fluicol;
    }
    if (state.input.tipoFluido == 0) {
        solutionGasRatio = state.cells[i].flui.RS(pmed, tmed);
        oilVolumeFactor = state.cells[i].flui.BOFunc(pmed, tmed, solutionGasRatio);
        solutionGasRatio = solutionGasRatio * 6.29 / 35.31467;
        waterVolumeFactor = state.cells[i].flui.BAFunc(pmed, tmed);
        double rhogstd = state.cells[i].flui.Deng * 1.225;
        double rhololeostd = 1000. * 141.5 / (131.5 + state.cells[i].flui.API);
        double rhoa = state.cells[i].flui.Denag * 1000.;

        if (titulo < 1. - 1e-15)
            state.cells[i].FW = state.cells[i].flui.BSW * waterVolumeFactor / (oilVolumeFactor + waterVolumeFactor * state.cells[i].flui.BSW - state.cells[i].flui.BSW * oilVolumeFactor);
        else
            state.cells[i].FW = 1.;
        double oilDensity = state.cells[i].flui.MasEspoleo(pmed, tmed);
        double waterDensity = state.cells[i].flui.MasEspAgua(pmed, tmed);
        double completionDensity = state.cells[i].fluicol.MasEspFlu(pmed, tmed);
        double denom = completionDensity * state.cells[i].bet + (1 - state.cells[i].bet) * state.cells[i].FW * waterDensity;
        if (titulo < 1. - 1e-15) {
            denom += (1 - state.cells[i].bet) * (1. - state.cells[i].FW) * oilDensity / (1. - titulo);
            state.cells[i].QL = state.cells[i].MC / denom;

            state.cells[i].Mliqini = state.cells[i].QL * (state.cells[i].bet * completionDensity + (1 - state.cells[i].bet) * (state.cells[i].FW * waterDensity +
                                                                                            (1. - state.cells[i].FW) * oilDensity));
        } else if (state.cells[i].flui.BSW > 1.e-15 || state.cells[i - 1].betI > 1.e-15) {
            state.cells[i].Mliqini = state.cells[i - 1].QL * state.cells[i - 1].rpCi * (1. - state.cells[i - 1].betI) * (1. - titV) +
                                state.cells[i - 1].fontemassLR * (1. - titF) +
                                state.cells[i - 1].QL * state.cells[i - 1].betI * state.cells[i - 1].rcCi + state.cells[i - 1].fontemassCR;
            state.cells[i].QL = state.cells[i].Mliqini / denom;
        } else {
            state.cells[i].Mliqini = 0.;
            state.cells[i].QL = 0.;
        }
    } else {
        double oilDensity = state.cells[i].flui.MasEspoleo(pmed, tmed);
        double completionDensity = state.cells[i].fluicol.MasEspFlu(pmed, tmed);
        double oilFlowRate = mHidro * (1. - titulo) / oilDensity;
        double completionFlowRate = mComp / completionDensity;
        if (fabs(oilFlowRate + completionFlowRate) > 1e-15)
            state.cells[i].bet = completionFlowRate / (oilFlowRate + completionFlowRate);
        else
            state.cells[i].bet = state.cells[i - 1].bet;
        double completionFraction = state.cells[i].bet;
        double denom;
        if (titulo < 1. - 1e-15) {
            denom = completionDensity * completionFraction + (1. - completionFraction) * oilDensity;
            state.cells[i].QL = (mHidro * (1. - titulo) + mComp) / denom;
            state.cells[i].Mliqini = state.cells[i].QL * (state.cells[i].bet * completionDensity + (1 - state.cells[i].bet) * oilDensity);
        } else if (completionFraction > (1 - 1e-15)) {
            state.cells[i].Mliqini = state.cells[i - 1].QL * state.cells[i - 1].betI * state.cells[i - 1].rcCi + state.cells[i - 1].fontemassCR;
            state.cells[i].QL = state.cells[i].Mliqini / (completionDensity * completionFraction + (1. - completionFraction) * oilDensity);
        } else {
            state.cells[i].Mliqini = 0.;
            state.cells[i].QL = 0.;
        }
    }

    state.cells[i - 1].MliqiniR = state.cells[i].Mliqini;
    if (i < state.lastCell)
        state.cells[i + 1].MliqiniL = state.cells[i].Mliqini;
    rhog = state.cells[i].rgCi = state.cells[i].flui.MasEspGas(pmed, tmed);

    // calculo das massas especificas na interface a esquerda da celula
    state.cells[i].rpCi = state.cells[i].flui.MasEspLiq(pmed, tmed);
    state.cells[i].rcCi = state.cells[i].fluicol.MasEspFlu(pmed, tmed);
    rhol = (1 - state.cells[i].bet) * state.cells[i].rpCi + state.cells[i].bet * state.cells[i].rcCi;
    // rhol = (1 - celula[i].bet) * celula[i].flui.MasEspLiq(pmed, tmed)
    // vazoes volumetricas:
    if (i < state.lastCell)
        state.cells[i + 1].QLL = state.cells[i].QL;
    state.cells[i - 1].QLR = state.cells[i].QL;
    state.cells[i].QG = (state.cells[i].MC - state.cells[i].Mliqini) / rhog;
    // Definicao das fracoes volumetricas:
    applyCompositionalPhaseRegime(state, i, rhol, rhog);
    if (fabs(state.cells[i].QL) > 1e-15 && fabs(state.cells[i].MComp) > 1e-15) {
        double dxmed = 0.5 * (state.cells[i - 1].dx + state.cells[i].dx);
        double hol0 = 1. - state.cells[i - 1].alf;
        double hol1 = 1. - state.cells[i].alf;
        state.cells[i].fluicol.TR = (state.cells[i - 1].fluicol.TR * state.cells[i - 1].MComp + state.cells[i - 1].fontemassCR * residenceTimeSource) / state.cells[i].MComp;
        state.cells[i].fluicol.TR = state.cells[i].fluicol.TR + (0.5 * state.cells[i - 1].duto.area * state.cells[i - 1].dx * hol0 +
                                                       0.5 * state.cells[i].duto.area * state.cells[i].dx * hol1) /
                                                          state.cells[i].QL;
    } else
        state.cells[i].fluicol.TR = state.cells[i - 1].fluicol.TR;
}

// ------------------------------------------------- mass march helpers ----
//
// The fourteen bodies T085 carved out of RenovaMassPerm. They were private to
// SProd and they stay private here: nothing outside advanceSteadyMass calls
// them, so they get internal linkage rather than a line in the header.
//
// Nine are the arms of the accessory dispatch, one per kind attached to the
// upstream cell. Five close the march once the sources are known, chosen by
// which phases are actually flowing.
namespace {

void applySteadyMassWithoutSource(const SteadyStateState &state, int i, int mudaRGO, double oilVolumeFactor, double solutionGasRatio, double tmed, double &oilVolumeFactorInSitu, double &waterVolumeFactorInSitu, double &waterCutInSitu) {
    // neste caso, variaveis como RGO de separador, BSW, API, densidade de gas e outras nÃƒÂ£o muda, sao iguais
    // aos valores da celula i-1
    double liquidHoldup = 1 - state.cells[i - 1].alf;
    double completionFraction = state.cells[i - 1].bet;
    double bsw = state.cells[i - 1].FW;
    double rhog = state.cells[i - 1].flui.MasEspGas(state.cells[i - 1].pres, state.cells[i - 1].temp);
    double rhogST = state.cells[i - 1].flui.Deng * 1.225;
    // o valor de volume de leve ÃƒÂ© atualizado neste ponto,
    // seguindo o equacionamento mostrado em relatorio, nÃƒÂ£o ÃƒÂ© relevante para o permanente, mas
    // deve ser calculado, pois ÃƒÂ© utilizado como entrada no transiente
    state.cells[i - 1].VolLeveST = (((1 - liquidHoldup) * rhog / rhogST) + liquidHoldup * (1 - completionFraction) * (1. - bsw) * solutionGasRatio / oilVolumeFactor);
    if (state.cells[i - 1].VolLeveST < 1e-15)
        state.cells[i - 1].VolLeveST = 0.;

    if (state.input.trackRGO == -1) {
        // esta chave nÃƒÂ£o ÃƒÂ© utilizada, ÃƒÂ© mantida aqui como reserva, atualmente este
        // calculo nÃƒÂ£o esta funcionando, lembrando que arq.trackRGO assume apenas 2 valores, 0 ou 1
        if (liquidHoldup > (*state.globals).localtiny && completionFraction < (1. - (*state.globals).localtiny) && bsw < (1. - (*state.globals).localtiny))
            state.cells[i].flui.RGO = state.cells[i - 1].VolLeveST * oilVolumeFactor / (liquidHoldup * (1 - completionFraction) * (1 - bsw));
        if (state.cells[i].flui.RGO > (*state.globals).RGOMax && mudaRGO == 1 && state.cells[i].flui.RGO < 1e7)
            state.cells[i].flui.RGO = (*state.globals).RGOMax;
    } else if (mudaRGO == 1)
        state.cells[i].flui.RGO = state.cells[i - 1].flui.RGO;

    state.cells[i].flui.BSW = state.cells[i - 1].flui.BSW;

    if (state.input.flashCompleto == 0) { // nesta chave se faz o carregamento na celula i
        // de variaveis importantes para o modelo black oil
        state.cells[i].flui.Deng = state.cells[i - 1].flui.Deng;
        state.cells[i].flui.yco2 = state.cells[i - 1].flui.yco2;
        if (state.productionFluidCount > 1 || (*state.globals).chaverede == 1) {
            state.cells[i].flui.API = state.cells[i - 1].flui.API;
            state.cells[i].flui.Denag = state.cells[i - 1].flui.Denag;
            state.cells[i].flui.TempL = state.cells[i - 1].flui.TempL;
            state.cells[i].flui.TempH = state.cells[i - 1].flui.TempH;
            state.cells[i].flui.LVisL = state.cells[i - 1].flui.LVisL;
            state.cells[i].flui.LVisH = state.cells[i - 1].flui.LVisH;
        }
        state.cells[i].flui.RenovaFluido();
        correctGasSpecificGravity(state, i);
    }
    if (state.cells[i - 1].bet > (*state.globals).localtiny * 1e-6) { // reavaliacao da fracao volumetrica do liquido complementar
        // mesmo que nÃƒÂ£o tenha fonte, ela pode mudar, devido ao encolhimento do liquido produzido
        if (state.cells[i].flui.RGO < 1e6)
            oilVolumeFactorInSitu = state.cells[i].flui.BOFunc(state.cells[i].presaux, tmed);
        else
            oilVolumeFactorInSitu = 1.;
        waterVolumeFactorInSitu = state.cells[i].flui.BAFunc(state.cells[i].presaux, tmed);
        waterCutInSitu = state.cells[i].flui.BSW * waterVolumeFactorInSitu / (oilVolumeFactorInSitu + waterVolumeFactorInSitu * state.cells[i].flui.BSW - state.cells[i].flui.BSW * oilVolumeFactorInSitu);
        double oilFlowRate = state.cells[i - 1].QL * (1 - state.cells[i - 1].FW) * (1 - state.cells[i - 1].bet) * oilVolumeFactorInSitu / oilVolumeFactor;
        double waterFlowRate;
        if (waterCutInSitu < (1 - (*state.globals).localtiny))
            waterFlowRate = waterCutInSitu * oilFlowRate / (1 - waterCutInSitu);
        else
            waterFlowRate = state.cells[i - 1].QL * (1 - state.cells[i - 1].bet);
        double completionFlowRate;
        if (state.cells[i].flui.RGO < 1e7)
            completionFlowRate = state.cells[i - 1].QL * (state.cells[i - 1].bet) * state.cells[i - 1].fluicol.MasEspFlu(state.cells[i - 1].pres, state.cells[i - 1].temp) / state.cells[i].fluicol.MasEspFlu(state.cells[i].presaux, tmed);
        else
            completionFlowRate = 0.;
        ////////////////////////esperar//////////////////////////////////////
        if ((fabs(oilFlowRate) + fabs(waterFlowRate) + fabs(completionFlowRate)) > 0)
            state.cells[i].bet = fabs(completionFlowRate) / (fabs(oilFlowRate) + fabs(waterFlowRate) + fabs(completionFlowRate));
        else
            state.cells[i].bet = 0.;
    } else
        state.cells[i].bet = state.cells[i - 1].bet;
}

void applySteadyMassDryGasInjection(const SteadyStateState &state, int i, int mudaRGO, double temperatureLow, double temperatureHigh, double oilVolumeFactor, double waterVolumeFactor, double solutionGasRatio, double tmed, double &residenceTimeSource, double &oilVolumeFactorInSitu, double &waterVolumeFactorInSitu, double &waterCutInSitu) {
    double apiGravity = state.cells[i - 1].flui.API;
    double rhololeo = 1000. * 141.5 / (131.5 + apiGravity);
    double rholis = state.cells[i - 1].flui.MasEspLiq(state.cells[i - 1].pres, state.cells[i - 1].temp);
    waterCutInSitu = 0.;
    // vazao de oleo standard
    double qostd;
    if (state.cells[i - 1].flui.dStockTankVaporMassFraction < 1. - 1e-15)
        qostd = state.cells[i - 1].Mliqini * (1. - state.cells[i - 1].FW) * (1. - state.cells[i - 1].bet) / (oilVolumeFactor * rholis);
    else
        qostd = 0.;
    // calculo da nova vazao de gas standard com a soma da fonte de gas:
    double qgstd = qostd * state.cells[i - 1].flui.RGO + state.cells[i - 1].acsr.injg.QGas / 86400;
    double deng;
    double yco2;
    // balanco que define a densidade de gas e a fracao de CO2 devido  aa fonte de gas
    if (fabs(qgstd) > (*state.globals).localtiny && fabs(state.cells[i - 1].acsr.injg.QGas) > 1e-15) {
        deng = (qostd * state.cells[i - 1].flui.RGO * state.cells[i - 1].flui.Deng + state.cells[i - 1].acsr.injg.QGas * state.cells[i - 1].acsr.injg.FluidoPro.Deng / 86400) / qgstd;
        yco2 = (qostd * state.cells[i - 1].flui.RGO * state.cells[i - 1].flui.yco2 + state.cells[i - 1].acsr.injg.QGas * state.cells[i - 1].acsr.injg.FluidoPro.yco2 / 86400) / qgstd;
    } else {
        deng = state.cells[i - 1].flui.Deng;
        yco2 = state.cells[i - 1].flui.yco2;
    }

    double liquidHoldup = 1 - state.cells[i - 1].alf;
    double completionFraction = state.cells[i - 1].bet;
    double bsw = state.cells[i - 1].FW;
    double rhog = state.cells[i - 1].flui.MasEspGas(state.cells[i - 1].pres, state.cells[i - 1].temp);
    double rhogST = state.cells[i - 1].flui.Deng * 1.225;
    // calculo do volume de leve na celula i, seguindo o equacionamento mostrado em relatorio,
    // nÃƒÂ£o ÃƒÂ© relevante para o permanente, mas deve ser calculado,
    // pois ÃƒÂ© utilizado como entrada no transiente
    state.cells[i - 1].VolLeveST = (((1 - liquidHoldup) * rhog / rhogST) + liquidHoldup * (1 - completionFraction) * (1. - bsw) * solutionGasRatio / oilVolumeFactor);
    if (state.cells[i - 1].VolLeveST < 1e-15)
        state.cells[i - 1].VolLeveST = 0.;
    if (state.input.trackRGO == -1) {
        // esta chave nÃƒÂ£o ÃƒÂ© utilizada, ÃƒÂ© mantida aqui como reserva, atualmente este
        // calculo nÃƒÂ£o esta funcionando, lembrando que arq.trackRGO assume apenas 2 valores, 0 ou 1
        if (liquidHoldup > (*state.globals).localtiny && completionFraction < (1. - (*state.globals).localtiny) && bsw < (1. - (*state.globals).localtiny) && mudaRGO == 1)
            state.cells[i].flui.RGO = state.cells[i - 1].VolLeveST * oilVolumeFactor / (liquidHoldup * (1 - completionFraction) * (1 - bsw));
        if (state.cells[i].flui.RGO > (*state.globals).RGOMax && mudaRGO == 1 && state.cells[i].flui.RGO < 1e7)
            state.cells[i].flui.RGO = (*state.globals).RGOMax;
    }
    double rgo;
    // calculo de novo RGO de separador - sem escorregamento - a partir das vazÃƒÂµes satndard de gas e oleo
    if (qostd > (*state.globals).localtiny && mudaRGO == 1)
        rgo = qgstd / qostd;
    else
        rgo = state.cells[i - 1].flui.RGO;
    state.cells[i].flui.RGO = rgo;

    state.cells[i].flui.BSW = state.cells[i - 1].flui.BSW; // bsw nÃƒÂ£o muda devido a uma fonte de gas

    if (state.input.flashCompleto == 0) { // nesta chave se faz o carregamento na celula i
        // de variaveis importantes para o modelo black oil
        state.cells[i].flui.Deng = deng;
        state.cells[i].flui.yco2 = yco2;
        if (state.productionFluidCount > 1 || (*state.globals).chaverede == 1) {
            state.cells[i].flui.API = state.cells[i - 1].flui.API;
            state.cells[i].flui.TempL = state.cells[i - 1].flui.TempL;
            state.cells[i].flui.TempH = state.cells[i - 1].flui.TempH;
            state.cells[i].flui.LVisL = state.cells[i - 1].flui.LVisL;
            state.cells[i].flui.LVisH = state.cells[i - 1].flui.LVisH;
        }
        state.cells[i].flui.RenovaFluido();
        correctGasSpecificGravity(state, i);
    }

    if (state.cells[i - 1].bet > (*state.globals).localtiny * 1e-6 && state.cells[i - 1].alf < 1. - (*state.globals).localtiny * 1e-6 && fabs(state.cells[i - 1].acsr.injg.QGas) > 1e-15) {
        // reavaliacao da fracao volumetrica do liquido complementar
        // mesmo que nÃƒÂ£o tenha fonte de liquido, ela pode mudar,
        // devido ao encolhimento do liquido produzido
        oilVolumeFactorInSitu = state.cells[i].flui.BOFunc(state.cells[i].presaux, tmed);
        waterVolumeFactorInSitu = state.cells[i].flui.BAFunc(state.cells[i].presaux, tmed);
        waterCutInSitu = state.cells[i].flui.BSW * waterVolumeFactorInSitu / (oilVolumeFactorInSitu + waterVolumeFactorInSitu * state.cells[i].flui.BSW - state.cells[i].flui.BSW * oilVolumeFactorInSitu);
        double oilFlowRate = state.cells[i - 1].QL * (1 - state.cells[i - 1].FW) * (1 - state.cells[i - 1].bet) * oilVolumeFactorInSitu / oilVolumeFactor;
        double waterFlowRate;
        if (waterCutInSitu < (1 - (*state.globals).localtiny))
            waterFlowRate = waterCutInSitu * oilFlowRate / (1 - waterCutInSitu);
        else
            waterFlowRate = state.cells[i - 1].QL * (1 - state.cells[i - 1].bet);
        double completionFlowRate = state.cells[i - 1].QL * (state.cells[i - 1].bet) * state.cells[i - 1].fluicol.MasEspFlu(state.cells[i - 1].pres, state.cells[i - 1].temp) / state.cells[i].fluicol.MasEspFlu(state.cells[i].presaux, tmed);
        state.cells[i].bet = fabs(completionFlowRate) / (fabs(oilFlowRate) + fabs(waterFlowRate) + fabs(completionFlowRate));
    } else
        state.cells[i].bet = state.cells[i - 1].bet;
}

void applySteadyMassWetGasInjection(const SteadyStateState &state, int i, int mudaRGO, double temperatureLow, double temperatureHigh, double oilVolumeFactor, double waterVolumeFactor, double solutionGasRatio, double tmed, double &residenceTimeSource, double &oilVolumeFactorInSitu, double &waterVolumeFactorInSitu, double &waterCutInSitu) {
    double apiGravity = state.cells[i - 1].flui.API;
    double rhololeo = 1000. * 141.5 / (131.5 + apiGravity);
    double rholis = state.cells[i - 1].flui.MasEspLiq(state.cells[i - 1].pres, state.cells[i - 1].temp);
    double rholisF;
    double boinjl;
    double bainjl;
    double fwinjl;
    residenceTimeSource = state.cells[i - 1].acsr.injg.fluidocol.TR;
    if (fabs(state.cells[i - 1].acsr.injg.QGas) > 1e-15) {
        rholisF = state.cells[i - 1].acsr.injg.FluidoPro.MasEspLiq(state.cells[i - 1].pres, state.cells[i - 1].temp);
        boinjl = state.cells[i - 1].acsr.injg.FluidoPro.BOFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        bainjl = state.cells[i - 1].acsr.injg.FluidoPro.BAFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        fwinjl = state.cells[i - 1].acsr.injg.FluidoPro.BSW * bainjl / (boinjl + bainjl * state.cells[i - 1].acsr.injg.FluidoPro.BSW - state.cells[i - 1].acsr.injg.FluidoPro.BSW * boinjl);
    } else {
        rholisF = rholis;
        boinjl = oilVolumeFactor;
        bainjl = waterVolumeFactor;
        fwinjl = state.cells[i - 1].FW;
    }
    // vazao de oleo standard
    double qostd = 0.;
    if (state.cells[i - 1].flui.dStockTankVaporMassFraction < 1. - 1e-15)
        qostd = state.cells[i - 1].Mliqini * (1. - state.cells[i - 1].FW) *
                (1. - state.cells[i - 1].bet) / (oilVolumeFactor * rholis);
    // vazao de oleo standard da fonte:
    double qostd2 = 0.;
    if (state.cells[i - 1].acsr.injg.FluidoPro.dStockTankVaporMassFraction < 1. - 1e-15)
        qostd2 = state.cells[i - 1].fontemassLR * (1. - fwinjl) / (boinjl * rholisF);
    // calculo da nova vazao de gas standard com a soma da fonte de gas:
    double qgstd = qostd * state.cells[i - 1].flui.RGO + state.cells[i - 1].acsr.injg.QGas / 86400;
    double deng;
    double yco2;
    // balanco que define a densidade de gas e a fracao de CO2 devido  aa fonte de gas
    if (fabs(qgstd) > (*state.globals).localtiny * 1e-10 && fabs(state.cells[i - 1].acsr.injg.QGas) > 1e-15) {
        deng = (qostd * state.cells[i - 1].flui.RGO * state.cells[i - 1].flui.Deng + state.cells[i - 1].acsr.injg.QGas * state.cells[i - 1].acsr.injg.FluidoPro.Deng / 86400) / qgstd;
        yco2 = (qostd * state.cells[i - 1].flui.RGO * state.cells[i - 1].flui.yco2 + state.cells[i - 1].acsr.injg.QGas * state.cells[i - 1].acsr.injg.FluidoPro.yco2 / 86400) / qgstd;
    } else {
        deng = state.cells[i - 1].flui.Deng;
        yco2 = state.cells[i - 1].flui.yco2;
    }

    double liquidHoldup = 1 - state.cells[i - 1].alf;
    double completionFraction = state.cells[i - 1].bet;
    double bsw = state.cells[i - 1].FW;
    double rhog = state.cells[i - 1].flui.MasEspGas(state.cells[i - 1].pres, state.cells[i - 1].temp);
    double rhogST = state.cells[i - 1].flui.Deng * 1.225;
    // calculo do volume de leve na celula i, seguindo o equacionamento mostrado em relatorio,
    // nÃƒÂ£o ÃƒÂ© relevante para o permanente, mas deve ser calculado,
    // pois ÃƒÂ© utilizado como entrada no transiente
    state.cells[i - 1].VolLeveST = (((1 - liquidHoldup) * rhog / rhogST) + liquidHoldup * (1 - completionFraction) * (1. - bsw) * solutionGasRatio / oilVolumeFactor);
    if (state.cells[i - 1].VolLeveST < 1e-15)
        state.cells[i - 1].VolLeveST = 0.;
    if (state.input.trackRGO == -1) {
        // esta chave nÃƒÂ£o ÃƒÂ© utilizada, ÃƒÂ© mantida aqui como reserva, atualmente este
        // calculo nÃƒÂ£o esta funcionando, lembrando que arq.trackRGO assume apenas 2 valores, 0 ou 1
        if (liquidHoldup > (*state.globals).localtiny && completionFraction < (1. - (*state.globals).localtiny) && bsw < (1. - (*state.globals).localtiny) && mudaRGO == 1)
            state.cells[i].flui.RGO = state.cells[i - 1].VolLeveST * oilVolumeFactor / (liquidHoldup * (1 - completionFraction) * (1 - bsw));
        if (state.cells[i].flui.RGO > (*state.globals).RGOMax && mudaRGO == 1 && state.cells[i].flui.RGO < 1e7)
            state.cells[i].flui.RGO = (*state.globals).RGOMax;
    }

    // calculo do novo RGO do separador modeificado pela fonte de liquido
    if (fabs(qostd + qostd2) > (*state.globals).localtiny && mudaRGO == 1 && fabs(state.cells[i - 1].acsr.injg.QGas) > 1e-15)
        state.cells[i].flui.RGO = qgstd / (qostd + qostd2);
    else if (mudaRGO == 1 && fabs(state.cells[i - 1].acsr.injg.QGas) > 1e-15)
        state.cells[i].flui.RGO = (*state.globals).RGOMax;
    else
        state.cells[i].flui.RGO = state.cells[i - 1].flui.RGO;

    if (state.productionFluidCount > 1 || (*state.globals).chaverede == 1) {
        // o modelo seja ASTM. Observar que isto sÃƒÂ³ faz sentido se se tiver mais de um fluido
        // de producao cadastrado no JSON, ou a simulacao se insere em um sistema de redes
        // com varios tramos alimentando outros tramos com fluidos com propriedades diferentes
        // entre cada tramo
        if (fabs(qostd + qostd2) > 1e-15 && state.input.flashCompleto == 0 && fabs(state.cells[i - 1].acsr.injg.QGas) > 1e-15) {
            double denmixSTD = (141.5 / (131.5 + state.cells[i - 1].flui.API)) * qostd + (141.5 / (131.5 + state.cells[i].flui.API)) * qostd2;
            denmixSTD /= (qostd + qostd2);
            state.cells[i].flui.API = 141.5 / denmixSTD - 131.5;
        } else if (state.input.flashCompleto == 0)
            state.cells[i].flui.API = state.cells[i - 1].flui.API;
        double waterFlowRate1;
        if ((1. - state.cells[i - 1].flui.BSW) > 0)
            waterFlowRate1 = qostd * state.cells[i - 1].flui.BSW / (1. - state.cells[i - 1].flui.BSW);
        else
            waterFlowRate1 = state.cells[i - 1].Mliqini * state.cells[i - 1].FW * (1. - state.cells[i - 1].bet) / (1000. * state.cells[i - 1].flui.Denag);
        double waterFlowRate2;
        if ((1. - state.cells[i - 1].acsr.injg.FluidoPro.BSW) > 0)
            waterFlowRate2 = qostd2 * state.cells[i - 1].acsr.injg.FluidoPro.BSW / (1. - state.cells[i - 1].acsr.injg.FluidoPro.BSW);
        else
            waterFlowRate2 = state.cells[i - 1].fontemassLR * fwinjl / (1000. * state.cells[i - 1].acsr.injg.FluidoPro.Denag);

        if (fabs(waterFlowRate1 + waterFlowRate2 + qostd + qostd2) > 1e-15 && fabs(state.cells[i - 1].acsr.injg.QGas) > 1e-15)
            state.cells[i].flui.BSW = (waterFlowRate1 + waterFlowRate2) / (waterFlowRate1 + waterFlowRate2 + qostd + qostd2);
        else
            state.cells[i].flui.BSW = state.cells[i - 1].flui.BSW;
        if (fabs(waterFlowRate1 + waterFlowRate2) > 1e-15)
            state.cells[i].flui.Denag = (waterFlowRate1 * state.cells[i - 1].flui.Denag +
                                    waterFlowRate2 * state.cells[i - 1].acsr.injg.FluidoPro.Denag) /
                                   (waterFlowRate1 + waterFlowRate2);
        else
            state.cells[i].flui.Denag = state.cells[i - 1].flui.Denag;
        // reavaliacao dos pares necessarios para o modelo de viscosidade de oleo morto
        // ASTM. Isto ÃƒÂ© feito tendo sempre os mesmos valores de temperatura maxima e
        // e minima para a construcao dos pares
        if (fabs(qostd + qostd2) > 1e-15 && state.input.flashCompleto == 0 && fabs(state.cells[i - 1].acsr.injg.QGas) > 1e-15) { // vazao de liquido >0
            state.cells[i].flui.TempL = temperatureLow;
            state.cells[i].flui.TempH = temperatureHigh;
            state.cells[i].flui.LVisL = (qostd * state.cells[i - 1].flui.VisOM(temperatureLow) + qostd2 * state.cells[i - 1].acsr.injg.FluidoPro.VisOM(temperatureLow)) / (qostd + qostd2);
            state.cells[i].flui.LVisH = (qostd * state.cells[i - 1].flui.VisOM(temperatureHigh) + qostd2 * state.cells[i - 1].acsr.injg.FluidoPro.VisOM(temperatureHigh)) / (qostd + qostd2);
        } else if (state.input.flashCompleto == 0 || state.cells[i - 1].acsr.injg.QGas <= 0.) { // se,m vazao de liquido
            state.cells[i].flui.TempL = state.cells[i - 1].flui.TempL;
            state.cells[i].flui.TempH = state.cells[i - 1].flui.TempH;
            state.cells[i].flui.LVisL = state.cells[i - 1].flui.LVisL;
            state.cells[i].flui.LVisH = state.cells[i - 1].flui.LVisH;
        }
    }

    if (state.input.flashCompleto == 0) { // nesta chave se faz o carregamento na celula i
        // de variaveis importantes para o modelo black oil
        state.cells[i].flui.Deng = deng;
        state.cells[i].flui.yco2 = yco2;
        if (state.productionFluidCount > 1 || (*state.globals).chaverede == 1) {
            state.cells[i].flui.API = state.cells[i - 1].flui.API;
            state.cells[i].flui.TempL = state.cells[i - 1].flui.TempL;
            state.cells[i].flui.TempH = state.cells[i - 1].flui.TempH;
            state.cells[i].flui.LVisL = state.cells[i - 1].flui.LVisL;
            state.cells[i].flui.LVisH = state.cells[i - 1].flui.LVisH;
        }
        state.cells[i].flui.RenovaFluido();
        correctGasSpecificGravity(state, i);
    }

    if (state.cells[i - 1].bet > (*state.globals).localtiny && state.cells[i - 1].alf < 1. - (*state.globals).localtiny * 1e-10 && fabs(state.cells[i - 1].acsr.injg.QGas) > 1e-15) {
        // reavaliacao da fracao volumetrica do liquido complementar
        double oilVolumeFactorNeighbour = state.cells[i].flui.BOFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        double waterVolumeFactorNeighbour = state.cells[i].flui.BAFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        double waterCutNeighbour = state.cells[i].flui.BSW * waterVolumeFactorNeighbour / (oilVolumeFactorNeighbour + waterVolumeFactorNeighbour * state.cells[i].flui.BSW - state.cells[i].flui.BSW * oilVolumeFactorNeighbour);
        oilVolumeFactorInSitu = state.cells[i].flui.BOFunc(state.cells[i].presaux, tmed);
        waterVolumeFactorInSitu = state.cells[i].flui.BAFunc(state.cells[i].presaux, tmed);
        waterCutInSitu = state.cells[i].flui.BSW * waterVolumeFactorInSitu / (oilVolumeFactorInSitu + waterVolumeFactorInSitu * state.cells[i].flui.BSW - state.cells[i].flui.BSW * oilVolumeFactorInSitu);
        double qlpF = state.cells[i - 1].fontemassLR / state.cells[i - 1].acsr.injg.FluidoPro.MasEspLiq(state.cells[i - 1].pres, state.cells[i - 1].temp);
        double qlcF = state.cells[i - 1].fontemassCR / state.cells[i - 1].acsr.injg.fluidocol.MasEspFlu(state.cells[i - 1].pres, state.cells[i - 1].temp);

        double betN = (state.cells[i - 1].QL * state.cells[i - 1].bet + qlcF) / (state.cells[i - 1].QL + qlpF + qlcF);
        double oilFlowRate = (state.cells[i - 1].QL + qlpF + qlcF) * (1 - waterCutNeighbour) * (1 - betN) * oilVolumeFactorInSitu / oilVolumeFactor;

        double waterFlowRate;
        if (waterCutInSitu < (1 - (*state.globals).localtiny))
            waterFlowRate = waterCutInSitu * oilFlowRate / (1 - waterCutInSitu);
        else
            waterFlowRate = state.cells[i - 1].QL * (1 - state.cells[i - 1].bet) + qlpF;
        //* celula[i - 1].fluicol.MasEspFlu(celula[i - 1].pres, celula[i - 1].temp)

        double completionFlowRate = (state.cells[i - 1].QL * (state.cells[i - 1].bet) + qlcF) * state.cells[i - 1].fluicol.MasEspFlu(state.cells[i - 1].pres, state.cells[i - 1].temp) / state.cells[i].fluicol.MasEspFlu(state.cells[i].presaux, tmed);

        state.cells[i].bet = fabs(completionFlowRate) / (fabs(oilFlowRate) + fabs(waterFlowRate) + fabs(completionFlowRate));
    } else
        state.cells[i].bet = state.cells[i - 1].bet;
}

void applySteadyMassLiquidInjection(const SteadyStateState &state, int i, int mudaRGO, double temperatureLow, double temperatureHigh, double oilVolumeFactor, double waterVolumeFactor, double solutionGasRatio, double tmed, double &residenceTimeSource, double &oilVolumeFactorInSitu, double &waterVolumeFactorInSitu, double &waterCutInSitu) {
    double apiGravity = state.cells[i - 1].flui.API;
    double rhololeo = 1000. * 141.5 / (131.5 + apiGravity);
    double rhololeoF = 1000. * 141.5 / (131.5 + state.cells[i - 1].acsr.injl.FluidoPro.API);
    double rholis = state.cells[i - 1].flui.MasEspLiq(state.cells[i - 1].pres, state.cells[i - 1].temp);
    double rholisF;
    double boinjl;
    double bainjl;
    double fwinjl;
    residenceTimeSource = state.cells[i - 1].acsr.injl.fluidocol.TR;
    if (fabs(state.cells[i - 1].acsr.injl.QLiq) > 1e-15) {
        rholisF = state.cells[i - 1].acsr.injl.FluidoPro.MasEspLiq(state.cells[i - 1].pres, state.cells[i - 1].temp);
        boinjl = state.cells[i - 1].acsr.injl.FluidoPro.BOFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        bainjl = state.cells[i - 1].acsr.injl.FluidoPro.BAFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        fwinjl = state.cells[i - 1].acsr.injl.FluidoPro.BSW * bainjl / (boinjl + bainjl * state.cells[i - 1].acsr.injl.FluidoPro.BSW - state.cells[i - 1].acsr.injl.FluidoPro.BSW * boinjl);
    } else {
        rholisF = rholis;
        boinjl = oilVolumeFactor;
        bainjl = waterVolumeFactor;
        fwinjl = state.cells[i - 1].FW;
    }
    // vazao de oleo sytandard antes da fonte:
    double qostd1 = 0.;
    if (state.cells[i - 1].flui.dStockTankVaporMassFraction < 1. - 1e-15)
        qostd1 = state.cells[i - 1].Mliqini * (1. - state.cells[i - 1].FW) * (1. - state.cells[i - 1].bet) / (oilVolumeFactor * rholis);
    // vazao de oleo standard da fonte:
    double qostd2 = 0.;
    if (state.cells[i - 1].acsr.injl.FluidoPro.dStockTankVaporMassFraction < 1. - 1e-15)
        qostd2 = state.cells[i - 1].fontemassLR * (1. - fwinjl) / (boinjl * rholisF);
    // vazao total de gas nas condicoes standard fruto da soma do gas transportado
    // da fronteira a esquerda
    // da celula e do gas associado da fonte de liquido
    double qgstd = qostd1 * state.cells[i - 1].flui.RGO + qostd2 * state.cells[i - 1].acsr.injl.FluidoPro.RGO;
    double liquidHoldup = 1 - state.cells[i - 1].alf;
    double completionFraction = state.cells[i - 1].bet;
    double bsw = state.cells[i - 1].FW;
    double rhog = state.cells[i - 1].flui.MasEspGas(state.cells[i - 1].pres, state.cells[i - 1].temp);
    double rhogST = state.cells[i - 1].flui.Deng * 1.225;
    // calculo do volume de leve na celula i, seguindo o equacionamento mostrado em relatorio,
    // nÃƒÂ£o ÃƒÂ© relevante para o permanente, mas deve ser calculado,
    // pois ÃƒÂ© utilizado como entrada no transiente
    state.cells[i - 1].VolLeveST = (((1 - liquidHoldup) * rhog / rhogST) + liquidHoldup * (1 - completionFraction) * (1. - bsw) * solutionGasRatio / oilVolumeFactor);
    if (state.cells[i - 1].VolLeveST < 1e-15)
        state.cells[i - 1].VolLeveST = 0.;

    if (fabs(state.cells[i - 1].acsr.injl.QLiq) > 1e-15) {
        if (state.input.trackRGO == -1) {
            // esta chave nÃƒÂ£o ÃƒÂ© utilizada, ÃƒÂ© mantida aqui como reserva, atualmente este
            // calculo nÃƒÂ£o esta funcionando, lembrando que arq.trackRGO assume apenas 2 valores, 0 ou 1
            if (liquidHoldup > (*state.globals).localtiny && completionFraction < (1. - (*state.globals).localtiny) && bsw < (1. - (*state.globals).localtiny) && mudaRGO == 1)
                state.cells[i].flui.RGO = state.cells[i - 1].VolLeveST * oilVolumeFactor / (liquidHoldup * (1 - completionFraction) * (1 - bsw));
            if (state.cells[i].flui.RGO > (*state.globals).RGOMax && mudaRGO == 1)
                state.cells[i].flui.RGO = (*state.globals).RGOMax;
        } else {
            // calculo do novo RGO do separador modeificado pela fonte de liquido
            if (fabs(qostd1 + qostd2) > (*state.globals).localtiny && mudaRGO == 1)
                state.cells[i].flui.RGO = qgstd / (qostd1 + qostd2);
            else if (mudaRGO == 1)
                state.cells[i].flui.RGO = (*state.globals).RGOMax;
        }
    } else
        state.cells[i].flui.RGO = state.cells[i - 1].flui.RGO;
    // calculo da nova densidade de gas nas condicoes standard modificada pela fonte de liquido
    if (fabs(qgstd) > (*state.globals).localtiny && state.input.flashCompleto == 0 && fabs(state.cells[i - 1].acsr.injl.QLiq) > 1e-15 &&
        fabs(qgstd) > 1e-15)
        state.cells[i].flui.Deng = (qostd1 * state.cells[i - 1].flui.RGO * state.cells[i - 1].flui.Deng + qostd2 * state.cells[i - 1].acsr.injl.FluidoPro.RGO * state.cells[i - 1].acsr.injl.FluidoPro.Deng) / qgstd;
    else if (fabs(state.cells[i - 1].acsr.injl.QLiq) <= 1e-15)
        state.cells[i].flui.Deng = state.cells[i - 1].flui.Deng;
    // calculo da nova fracao de co2 modificad pela fonte de liquido
    if (fabs(qgstd) > (*state.globals).localtiny && state.input.flashCompleto == 0 && fabs(state.cells[i - 1].acsr.injl.QLiq) > 1e-15 &&
        fabs(qgstd) > 1e-15)
        state.cells[i].flui.yco2 = (qostd1 * state.cells[i - 1].flui.RGO * state.cells[i - 1].flui.yco2 + qostd2 * state.cells[i - 1].acsr.injl.FluidoPro.RGO * state.cells[i - 1].acsr.injl.FluidoPro.yco2) / qgstd;
    else if (fabs(state.cells[i - 1].acsr.injl.QLiq) <= 1e-15)
        state.cells[i].flui.yco2 = state.cells[i - 1].flui.yco2;
    if (state.productionFluidCount > 1 || (*state.globals).chaverede == 1) {
        // o modelo seja ASTM. Observar que isto sÃƒÂ³ faz sentido se se tiver mais de um fluido
        // de producao cadastrado no JSON, ou a simulacao se insere em um sistema de redes
        // com varios tramos alimentando outros tramos com fluidos com propriedades diferentes
        // entre cada tramo
        if (fabs(qostd1 + qostd2) > 1e-15 && state.input.flashCompleto == 0 &&
            fabs(state.cells[i - 1].acsr.injl.QLiq) > 1e-15) {
            double denmixSTD = (141.5 / (131.5 + state.cells[i - 1].flui.API)) * qostd1 +
                               (141.5 / (131.5 + state.cells[i - 1].acsr.injl.FluidoPro.API)) * qostd2;
            denmixSTD /= (qostd1 + qostd2);
            state.cells[i].flui.API = 141.5 / denmixSTD - 131.5;
        } else if (state.input.flashCompleto == 0 ||
                   fabs(state.cells[i - 1].acsr.injl.QLiq) <= 1e-15)
            state.cells[i].flui.API = state.cells[i - 1].flui.API;
        double waterFlowRate1;
        if ((1. - state.cells[i - 1].flui.BSW) > 0)
            waterFlowRate1 = qostd1 * state.cells[i - 1].flui.BSW / (1. - state.cells[i - 1].flui.BSW);
        else
            waterFlowRate1 = state.cells[i - 1].Mliqini * state.cells[i - 1].FW * (1. - state.cells[i - 1].bet) / (1000. * state.cells[i - 1].flui.Denag);
        double waterFlowRate2;
        if ((1. - state.cells[i - 1].acsr.injl.FluidoPro.BSW) > 0)
            waterFlowRate2 = qostd2 * state.cells[i - 1].acsr.injl.FluidoPro.BSW / (1. - state.cells[i - 1].acsr.injl.FluidoPro.BSW);
        else
            waterFlowRate2 = state.cells[i - 1].fontemassLR * fwinjl / (1000. * state.cells[i - 1].acsr.injl.FluidoPro.Denag);

        if (fabs(waterFlowRate1 + waterFlowRate2 + qostd1 + qostd2) > 1e-15 &&
            fabs(state.cells[i - 1].acsr.injl.QLiq) > 1e-15)
            state.cells[i].flui.BSW = (waterFlowRate1 + waterFlowRate2) / (waterFlowRate1 + waterFlowRate2 + qostd1 + qostd2);
        else
            state.cells[i].flui.BSW = state.cells[i - 1].flui.BSW;
        if (fabs(waterFlowRate1 + waterFlowRate2) > 1e-15)
            state.cells[i].flui.Denag = (waterFlowRate1 * state.cells[i - 1].flui.Denag +
                                    waterFlowRate2 * state.cells[i - 1].acsr.injl.FluidoPro.Denag) /
                                   (waterFlowRate1 + waterFlowRate2);
        else
            state.cells[i].flui.Denag = state.cells[i - 1].flui.Denag;
        // reavaliacao dos pares necessarios para o modelo de viscosidade de oleo morto
        // ASTM. Isto ÃƒÂ© feito tendo sempre os mesmos valores de temperatura maxima e
        // e minima para a construcao dos pares
        if (fabs(qostd1 + qostd2) > 1e-15 && state.input.flashCompleto == 0 &&
            fabs(state.cells[i - 1].acsr.injl.QLiq) > 1e-15) { // vazao de liquido >0
            state.cells[i].flui.TempL = temperatureLow;
            state.cells[i].flui.TempH = temperatureHigh;
            state.cells[i].flui.LVisL = (qostd1 * state.cells[i - 1].flui.VisOM(temperatureLow) + qostd2 * state.cells[i - 1].acsr.injl.FluidoPro.VisOM(temperatureLow)) / (qostd1 + qostd2);
            state.cells[i].flui.LVisH = (qostd1 * state.cells[i - 1].flui.VisOM(temperatureHigh) + qostd2 * state.cells[i - 1].acsr.injl.FluidoPro.VisOM(temperatureHigh)) / (qostd1 + qostd2);
        } else if (state.input.flashCompleto == 0) { // se,m vazao de liquido
            state.cells[i].flui.TempL = state.cells[i - 1].flui.TempL;
            state.cells[i].flui.TempH = state.cells[i - 1].flui.TempH;
            state.cells[i].flui.LVisL = state.cells[i - 1].flui.LVisL;
            state.cells[i].flui.LVisH = state.cells[i - 1].flui.LVisH;
        }
    }
    if (state.input.flashCompleto == 0) {
        state.cells[i].flui.RenovaFluido();
        correctGasSpecificGravity(state, i);
    }
    if ((state.cells[i - 1].bet > (*state.globals).localtiny * 1e-6 || fabs(state.cells[i - 1].fontemassCR) > (*state.globals).localtiny * 1e-6) &&
        fabs(state.cells[i - 1].acsr.injl.QLiq) > 1e-15) {
        // reavaliacao da fracao volumetrica do liquido complementar
        // observar que a fonte de liquido pode ter uma fracao de liquido complementar
        // distinta da onservada a esquerda da celula i-1
        double oilVolumeFactorNeighbour = state.cells[i].flui.BOFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        double waterVolumeFactorNeighbour = state.cells[i].flui.BAFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        double waterCutNeighbour = state.cells[i].flui.BSW * waterVolumeFactorNeighbour / (oilVolumeFactorNeighbour + waterVolumeFactorNeighbour * state.cells[i].flui.BSW - state.cells[i].flui.BSW * oilVolumeFactorNeighbour);
        oilVolumeFactorInSitu = state.cells[i].flui.BOFunc(state.cells[i].presaux, tmed);
        waterVolumeFactorInSitu = state.cells[i].flui.BAFunc(state.cells[i].presaux, tmed);
        waterCutInSitu = state.cells[i].flui.BSW * waterVolumeFactorInSitu / (oilVolumeFactorInSitu + waterVolumeFactorInSitu * state.cells[i].flui.BSW - state.cells[i].flui.BSW * oilVolumeFactorInSitu);
        double qlpF = state.cells[i - 1].fontemassLR / state.cells[i - 1].acsr.injl.FluidoPro.MasEspLiq(state.cells[i - 1].pres, state.cells[i - 1].temp);
        double qlcF = state.cells[i - 1].fontemassCR / state.cells[i - 1].acsr.injl.fluidocol.MasEspFlu(state.cells[i - 1].pres, state.cells[i - 1].temp);
        double betN = (state.cells[i - 1].QL * state.cells[i - 1].bet + qlcF) / (state.cells[i - 1].QL + qlpF + qlcF);
        double oilFlowRate = (state.cells[i - 1].QL + qlpF + qlcF) * (1 - waterCutNeighbour) * (1 - betN) * oilVolumeFactorInSitu / oilVolumeFactor;
        double waterFlowRate;
        if (waterCutInSitu < (1 - (*state.globals).localtiny))
            waterFlowRate = waterCutInSitu * oilFlowRate / (1 - waterCutInSitu);
        else
            waterFlowRate = state.cells[i - 1].QL * (1 - state.cells[i - 1].bet) + qlpF;
        double completionFlowRate = (state.cells[i - 1].QL * (state.cells[i - 1].bet) + qlcF) * state.cells[i - 1].fluicol.MasEspFlu(state.cells[i - 1].pres, state.cells[i - 1].temp) / state.cells[i].fluicol.MasEspFlu(state.cells[i].presaux, tmed);
        state.cells[i].bet = fabs(completionFlowRate) / (fabs(oilFlowRate) + fabs(waterFlowRate) + fabs(completionFlowRate));
    } else
        state.cells[i].bet = state.cells[i - 1].bet;
}

void applySteadyMassInflowPerformance(const SteadyStateState &state, int i, int mudaRGO, double temperatureLow, double temperatureHigh, double oilVolumeFactor, double waterVolumeFactor, double solutionGasRatio, double tmed, double &residenceTimeSource, double &oilVolumeFactorInSitu, double &waterVolumeFactorInSitu, double &waterCutInSitu) {
    double apiGravity = state.cells[i - 1].flui.API;
    double rhololeo = 1000. * 141.5 / (131.5 + apiGravity);
    double rhololeoF = 1000. * 141.5 / (131.5 + state.cells[i - 1].acsr.ipr.FluidoPro.API);
    double rholis = state.cells[i - 1].flui.MasEspLiq(state.cells[i - 1].pres, state.cells[i - 1].temp);
    double rholisF;
    double boipr;
    double baipr;
    double fwipr;
    residenceTimeSource = 0.;
    if (fabs(state.cells[i - 1].fontemassLR) > 1e-15) {
        rholisF = state.cells[i - 1].acsr.ipr.FluidoPro.MasEspLiq(state.cells[i - 1].pres, state.cells[i - 1].temp);
        boipr = state.cells[i - 1].acsr.ipr.FluidoPro.BOFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        baipr = state.cells[i - 1].acsr.ipr.FluidoPro.BAFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        fwipr = state.cells[i - 1].acsr.ipr.FluidoPro.BSW * baipr / (boipr + baipr * state.cells[i - 1].acsr.ipr.FluidoPro.BSW - state.cells[i - 1].acsr.ipr.FluidoPro.BSW * boipr);
    } else {
        rholisF = rholis;
        boipr = oilVolumeFactor;
        baipr = waterVolumeFactor;
        fwipr = state.cells[i - 1].FW;
    }
    // vazao de oleo standard antes da fonte:
    double qostd1 = 0.;
    if (state.cells[i - 1].flui.dStockTankVaporMassFraction < 1. - 1e-15)
        qostd1 = state.cells[i - 1].Mliqini * (1. - state.cells[i - 1].FW) * (1. - state.cells[i - 1].bet) / (oilVolumeFactor * rholis);
    // vazao de oleo standard da fonte:
    double qostd2 = 0.;
    if (state.cells[i - 1].acsr.ipr.FluidoPro.dStockTankVaporMassFraction < 1. - 1e-15)
        qostd2 = state.cells[i - 1].fontemassLR * (1. - fwipr) / (boipr * rholisF);
    // vazao total de gas nas condicoes standard fruto da soma do gas transportado
    // da fronteira a esquerda
    // da celula e do gas associado da IPR
    double qgstd = qostd1 * state.cells[i - 1].flui.RGO + qostd2 * state.cells[i - 1].acsr.ipr.FluidoPro.RGO;
    double liquidHoldup = 1 - state.cells[i - 1].alf;
    double completionFraction = state.cells[i - 1].bet;
    double bsw = state.cells[i - 1].FW;
    double rhog = state.cells[i - 1].flui.MasEspGas(state.cells[i - 1].pres, state.cells[i - 1].temp);
    double rhogST = state.cells[i - 1].flui.Deng * 1.225;
    // calculo do volume de leve na celula i, seguindo o equacionamento mostrado em relatorio,
    // nÃƒÂ£o ÃƒÂ© relevante para o permanente, mas deve ser calculado,
    // pois ÃƒÂ© utilizado como entrada no transiente
    state.cells[i - 1].VolLeveST = (((1 - liquidHoldup) * rhog / rhogST) + liquidHoldup * (1 - completionFraction) * (1. - bsw) * solutionGasRatio / oilVolumeFactor);
    if (state.cells[i - 1].VolLeveST < 1e-15)
        state.cells[i - 1].VolLeveST = 0.;
    if (state.input.trackRGO == -1) {
        // esta chave nÃƒÂ£o ÃƒÂ© utilizada, ÃƒÂ© mantida aqui como reserva, atualmente este
        // calculo nÃƒÂ£o esta funcionando, lembrando que arq.trackRGO assume apenas 2 valores, 0 ou 1
        if ((state.cells[i - 1].fontemassLR) > 1e-15 && (state.cells[i - 1].fontemassGR) > 1e-15) {
            if (liquidHoldup > (*state.globals).localtiny && completionFraction < (1. - (*state.globals).localtiny) && bsw < (1. - (*state.globals).localtiny) &&
                mudaRGO == 1 && fabs(state.cells[i - 1].fontemassLR) > 1e-15)
                state.cells[i].flui.RGO = state.cells[i - 1].VolLeveST * oilVolumeFactor / (liquidHoldup * (1 - completionFraction) * (1 - bsw));
            if (state.cells[i].flui.RGO > (*state.globals).RGOMax && mudaRGO == 1)
                state.cells[i].flui.RGO = (*state.globals).RGOMax;
        } else
            state.cells[i].flui.RGO = state.cells[i - 1].flui.RGO;

    } else { // calculo do novo RGO do separador modificado pela IPR
        if ((state.cells[i - 1].fontemassLR) > 1e-15 && (state.cells[i - 1].fontemassGR) > 1e-15) {
            if (fabs(qostd1 + qostd2) > (*state.globals).localtiny * 1e-10 && mudaRGO == 1)
                state.cells[i].flui.RGO = qgstd / (qostd1 + qostd2);
            else if (mudaRGO == 1)
                state.cells[i].flui.RGO = (*state.globals).RGOMax;
        } else
            state.cells[i].flui.RGO = state.cells[i - 1].flui.RGO;
    }
    // calculo da nova densidade de gas nas condicoes standard modificada pela IPR
    if (state.input.flashCompleto == 0 && (state.cells[i - 1].fontemassLR) > 1e-15 && (state.cells[i - 1].fontemassGR) > 1e-15 &&
        fabs(qgstd) > 1e-15)
        state.cells[i].flui.Deng = (qostd1 * state.cells[i - 1].flui.RGO * state.cells[i - 1].flui.Deng + qostd2 * state.cells[i - 1].acsr.ipr.FluidoPro.RGO * state.cells[i - 1].acsr.ipr.FluidoPro.Deng) / qgstd;
    else
        state.cells[i].flui.Deng = state.cells[i - 1].flui.Deng;
    // calculo da nova fracao de co2 modificada pela fonte de liquido
    if (state.input.flashCompleto == 0 && (state.cells[i - 1].fontemassLR) > 1e-15 && (state.cells[i - 1].fontemassGR) > 1e-15 && fabs(qgstd) > 1e-15)
        state.cells[i].flui.yco2 = (qostd1 * state.cells[i - 1].flui.RGO * state.cells[i - 1].flui.yco2 + qostd2 * state.cells[i - 1].acsr.ipr.FluidoPro.RGO * state.cells[i - 1].acsr.ipr.FluidoPro.yco2) / qgstd;
    else
        state.cells[i].flui.yco2 = state.cells[i - 1].flui.yco2;
    if (state.productionFluidCount > 1 || (*state.globals).chaverede == 1) {
        // o modelo seja ASTM. Observar que isto sÃƒÂ³ faz sentido se se tiver mais de um fluido
        // de producao cadastrado no JSON, ou a simulacao se insere em um sistema de redes
        // com varios tramos alimentando outros tramos com fluidos com propriedades diferentes
        // entre cada tramo
        if (fabs(qostd1 + qostd2) > 1e-15 && state.input.flashCompleto == 0 &&
            (state.cells[i - 1].fontemassLR) > 1e-15) {
            double denmixSTD = (141.5 / (131.5 + state.cells[i - 1].flui.API)) * qostd1 +
                               (141.5 / (131.5 + state.cells[i - 1].acsr.ipr.FluidoPro.API)) * qostd2;
            denmixSTD /= (qostd1 + qostd2);
            state.cells[i].flui.API = 141.5 / denmixSTD - 131.5;
        } else if (state.input.flashCompleto == 0 || (state.cells[i - 1].fontemassLR) <= 1e-15)
            state.cells[i].flui.API = state.cells[i - 1].flui.API;
        double waterFlowRate1;
        if ((1. - state.cells[i - 1].flui.BSW) > 0)
            waterFlowRate1 = qostd1 * state.cells[i - 1].flui.BSW / (1. - state.cells[i - 1].flui.BSW);
        else
            waterFlowRate1 = state.cells[i - 1].Mliqini * state.cells[i - 1].FW * (1. - state.cells[i - 1].bet) / (1000. * state.cells[i - 1].flui.Denag);
        double waterFlowRate2;
        if ((1. - state.cells[i - 1].acsr.ipr.FluidoPro.BSW) > 0)
            waterFlowRate2 = qostd2 * state.cells[i - 1].acsr.ipr.FluidoPro.BSW / (1. - state.cells[i - 1].acsr.ipr.FluidoPro.BSW);
        else
            waterFlowRate2 = state.cells[i - 1].fontemassLR * fwipr / (1000. * state.cells[i - 1].acsr.ipr.FluidoPro.Denag);

        if (fabs(waterFlowRate1 + waterFlowRate2 + qostd1 + qostd2) > 1e-15 &&
            (state.cells[i - 1].fontemassLR) > 1e-15)
            state.cells[i].flui.BSW = ((waterFlowRate1 + qostd1) * state.cells[i - 1].flui.BSW + (waterFlowRate2 + qostd2) * state.cells[i - 1].acsr.ipr.FluidoPro.BSW) / (waterFlowRate1 + waterFlowRate2 + qostd1 + qostd2);
        else
            state.cells[i].flui.BSW = state.cells[i - 1].flui.BSW;
        if ((waterFlowRate1 + waterFlowRate2) > 1e-15)
            state.cells[i].flui.Denag = (waterFlowRate1 * state.cells[i - 1].flui.Denag +
                                    waterFlowRate2 * state.cells[i - 1].acsr.ipr.FluidoPro.Denag) /
                                   (waterFlowRate1 + waterFlowRate2);
        else
            state.cells[i].flui.Denag = state.cells[i - 1].flui.Denag;
        if ((qostd1 + qostd2) > 1e-15 && state.input.flashCompleto == 0 && (state.cells[i - 1].fontemassLR) > 1e-15) {
            // reavaliacao dos pares necessarios para o modelo de viscosidade de oleo morto
            // ASTM. Isto ÃƒÂ© feito tendo sempre os mesmos valores de temperatura maxima e
            // e minima para a construcao dos pares
            state.cells[i].flui.TempL = temperatureLow;
            state.cells[i].flui.TempH = temperatureHigh;
            state.cells[i].flui.LVisL = (qostd1 * state.cells[i - 1].flui.VisOM(temperatureLow) + qostd2 * state.cells[i - 1].acsr.ipr.FluidoPro.VisOM(temperatureLow)) / (qostd1 + qostd2);
            state.cells[i].flui.LVisH = (qostd1 * state.cells[i - 1].flui.VisOM(temperatureHigh) + qostd2 * state.cells[i - 1].acsr.ipr.FluidoPro.VisOM(temperatureHigh)) / (qostd1 + qostd2);
        } else if (state.input.flashCompleto == 0 || fabs(state.cells[i - 1].fontemassLR) <= 1e-15) { // sem vazao de liquido
            state.cells[i].flui.TempL = state.cells[i - 1].flui.TempL;
            state.cells[i].flui.TempH = state.cells[i - 1].flui.TempH;
            state.cells[i].flui.LVisL = state.cells[i - 1].flui.LVisL;
            state.cells[i].flui.LVisH = state.cells[i - 1].flui.LVisH;
        }
    }
    if (state.input.flashCompleto == 0) {
        state.cells[i].flui.RenovaFluido();
        correctGasSpecificGravity(state, i);
    }
    if (state.cells[i - 1].bet > (*state.globals).localtiny * 1e-6 && fabs(state.cells[i - 1].fontemassLR) > 1e-15) {
        // reavaliacao da fracao volumetrica do liquido complementar
        // observar que a IPR pode ter uma fracao de liquido complementar
        // distinta da onservada a esquerda da celula i-1
        double oilVolumeFactorNeighbour = state.cells[i].flui.BOFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        double waterVolumeFactorNeighbour = state.cells[i].flui.BAFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        double waterCutNeighbour = state.cells[i].flui.BSW * waterVolumeFactorNeighbour / (oilVolumeFactorNeighbour + waterVolumeFactorNeighbour * state.cells[i].flui.BSW - state.cells[i].flui.BSW * oilVolumeFactorNeighbour);
        oilVolumeFactorInSitu = state.cells[i].flui.BOFunc(state.cells[i].presaux, tmed);
        waterVolumeFactorInSitu = state.cells[i].flui.BAFunc(state.cells[i].presaux, tmed);
        waterCutInSitu = state.cells[i].flui.BSW * waterVolumeFactorInSitu / (oilVolumeFactorInSitu + waterVolumeFactorInSitu * state.cells[i].flui.BSW - state.cells[i].flui.BSW * oilVolumeFactorInSitu);
        double qlpF = state.cells[i - 1].fontemassLR / state.cells[i - 1].acsr.ipr.FluidoPro.MasEspLiq(state.cells[i - 1].pres, state.cells[i - 1].temp);
        double qlcF = 0.;
        double betN = (state.cells[i - 1].QL * state.cells[i - 1].bet + qlcF) / (state.cells[i - 1].QL + qlpF + qlcF);
        double oilFlowRate = (state.cells[i - 1].QL + qlpF + qlcF) * (1 - waterCutNeighbour) * (1 - betN) * oilVolumeFactorInSitu / oilVolumeFactor;
        double waterFlowRate;
        if (waterCutInSitu < (1 - (*state.globals).localtiny))
            waterFlowRate = waterCutInSitu * oilFlowRate / (1 - waterCutInSitu);
        else
            waterFlowRate = state.cells[i - 1].QL * (1 - state.cells[i - 1].bet) + qlpF;
        double completionFlowRate = (state.cells[i - 1].QL * (state.cells[i - 1].bet) + qlcF) * state.cells[i - 1].fluicol.MasEspFlu(state.cells[i - 1].pres, state.cells[i - 1].temp) / state.cells[i].fluicol.MasEspFlu(state.cells[i].presaux, tmed);
        state.cells[i].bet = fabs(completionFlowRate) / (fabs(oilFlowRate) + fabs(waterFlowRate) + fabs(completionFlowRate));
    } else
        state.cells[i].bet = state.cells[i - 1].bet;
}

void applySteadyMassMultipleSource(const SteadyStateState &state, int i, int mudaRGO, double temperatureLow, double temperatureHigh, double oilVolumeFactor, double waterVolumeFactor, double solutionGasRatio, double tmed, double &residenceTimeSource, double &oilVolumeFactorInSitu, double &waterVolumeFactorInSitu, double &waterCutInSitu) {
    double apiGravity = state.cells[i - 1].flui.API;
    double rhololeo = 1000. * 141.5 / (131.5 + apiGravity);
    double rhololeoF = 1000. * 141.5 / (131.5 + state.cells[i - 1].acsr.injm.FluidoPro.API);
    double rholis = state.cells[i - 1].flui.MasEspLiq(state.cells[i - 1].pres, state.cells[i - 1].temp);
    double rholisF;
    double boinjl;
    double bainjl;
    double fwinjl;
    residenceTimeSource = state.cells[i - 1].acsr.injm.fluidocol.TR;
    if (fabs(state.cells[i - 1].acsr.injm.MassG + state.cells[i - 1].acsr.injm.MassP) > 1e-15) {
        rholisF = state.cells[i - 1].acsr.injm.FluidoPro.MasEspLiq(state.cells[i - 1].pres, state.cells[i - 1].temp);
        boinjl = state.cells[i - 1].acsr.injm.FluidoPro.BOFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        bainjl = state.cells[i - 1].acsr.injm.FluidoPro.BAFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        fwinjl = state.cells[i - 1].acsr.injm.FluidoPro.BSW * bainjl / (boinjl + bainjl * state.cells[i - 1].acsr.injm.FluidoPro.BSW - state.cells[i - 1].acsr.injm.FluidoPro.BSW * boinjl);
    } else {
        rholisF = rholis;
        boinjl = oilVolumeFactor;
        bainjl = waterVolumeFactor;
        fwinjl = state.cells[i - 1].FW;
    }
    // vazao de oleo sytandard antes da fonte:
    double qostd1 = 0.;
    if (state.cells[i - 1].flui.dStockTankVaporMassFraction < 1. - 1e-15)
        qostd1 = state.cells[i - 1].Mliqini * (1. - state.cells[i - 1].FW) * (1. - state.cells[i - 1].bet) / (oilVolumeFactor * rholis);
    // vazao de oleo standard da fonte:
    double qostd2 = 0.;
    if (state.cells[i - 1].acsr.injm.FluidoPro.dStockTankVaporMassFraction < 1. - 1e-15)
        qostd2 = state.cells[i - 1].fontemassLR * (1. - fwinjl) / (boinjl * rholisF);
    // vazao total de gas nas condicoes standard fruto da soma do gas transportado
    // da fronteira a esquerda
    // da celula e do gas associado da fonte de liquido
    double qgstd = qostd1 * state.cells[i - 1].flui.RGO + qostd2 * state.cells[i - 1].acsr.injm.FluidoPro.RGO;
    double liquidHoldup = 1 - state.cells[i - 1].alf;
    double completionFraction = state.cells[i - 1].bet;
    double bsw = state.cells[i - 1].FW;
    double rhog = state.cells[i - 1].flui.MasEspGas(state.cells[i - 1].pres, state.cells[i - 1].temp);
    double rhogST = state.cells[i - 1].flui.Deng * 1.225;
    // calculo do volume de leve na celula i, seguindo o equacionamento mostrado em relatorio,
    // nÃƒÂ£o ÃƒÂ© relevante para o permanente, mas deve ser calculado,
    // pois ÃƒÂ© utilizado como entrada no transiente
    state.cells[i - 1].VolLeveST = (((1 - liquidHoldup) * rhog / rhogST) + liquidHoldup * (1 - completionFraction) * (1. - bsw) * solutionGasRatio / oilVolumeFactor);
    if (state.cells[i - 1].VolLeveST < 1e-15)
        state.cells[i - 1].VolLeveST = 0.;

    if (state.input.trackRGO == -1) {
        // esta chave nÃƒÂ£o ÃƒÂ© utilizada, ÃƒÂ© mantida aqui como reserva, atualmente este
        // calculo nÃƒÂ£o esta funcionando, lembrando que arq.trackRGO assume apenas 2 valores, 0 ou 1
        if ((state.cells[i - 1].acsr.injm.MassG + state.cells[i - 1].acsr.injm.MassP) > 1e-15) {
            if (liquidHoldup > (*state.globals).localtiny && completionFraction < (1. - (*state.globals).localtiny) && bsw < (1. - (*state.globals).localtiny) && mudaRGO == 1)
                state.cells[i].flui.RGO = state.cells[i - 1].VolLeveST * oilVolumeFactor / (liquidHoldup * (1 - completionFraction) * (1 - bsw));
            if (state.cells[i].flui.RGO > (*state.globals).RGOMax && mudaRGO == 1)
                state.cells[i].flui.RGO = (*state.globals).RGOMax;
        } else
            state.cells[i].flui.RGO = state.cells[i - 1].flui.RGO;
    } else {
        // calculo do novo RGO do separador modeificado pela fonte de liquido
        if ((state.cells[i - 1].acsr.injm.MassG + state.cells[i - 1].acsr.injm.MassP) > 1e-15) {
            if (fabs(qostd1 + qostd2) > (*state.globals).localtiny && mudaRGO == 1)
                state.cells[i].flui.RGO = qgstd / (qostd1 + qostd2);
            else if (mudaRGO == 1)
                state.cells[i].flui.RGO = (*state.globals).RGOMax;
        } else
            state.cells[i].flui.RGO = state.cells[i - 1].flui.RGO;
    }
    // calculo da nova densidade de gas nas condicoes standard modificada pela fonte de liquido
    if (fabs(qgstd) > (*state.globals).localtiny * 1e-10 && state.input.flashCompleto == 0 &&
        (state.cells[i - 1].acsr.injm.MassG + state.cells[i - 1].acsr.injm.MassP) > 1e-15)
        state.cells[i].flui.Deng = (qostd1 * state.cells[i - 1].flui.RGO * state.cells[i - 1].flui.Deng + qostd2 * state.cells[i - 1].acsr.injm.FluidoPro.RGO * state.cells[i - 1].acsr.injm.FluidoPro.Deng) / qgstd;
    else
        state.cells[i].flui.RGO = state.cells[i - 1].flui.RGO;
    // calculo da nova fracao de co2 modificad pela fonte de liquido
    if (fabs(qgstd) > (*state.globals).localtiny * 1e-10 && state.input.flashCompleto == 0 &&
        (state.cells[i - 1].acsr.injm.MassG + state.cells[i - 1].acsr.injm.MassP) > 1e-15)
        state.cells[i].flui.yco2 = (qostd1 * state.cells[i - 1].flui.RGO * state.cells[i - 1].flui.yco2 + qostd2 * state.cells[i - 1].acsr.injm.FluidoPro.RGO * state.cells[i - 1].acsr.injm.FluidoPro.yco2) / qgstd;
    else
        state.cells[i].flui.yco2 = state.cells[i - 1].flui.yco2;
    if (state.productionFluidCount > 1 || (*state.globals).chaverede == 1) {
        // o modelo seja ASTM. Observar que isto sÃƒÂ³ faz sentido se se tiver mais de um fluido
        // de producao cadastrado no JSON, ou a simulacao se insere em um sistema de redes
        // com varios tramos alimentando outros tramos com fluidos com propriedades diferentes
        // entre cada tramo
        if (fabs(qostd1 + qostd2) > 1e-15 && state.input.flashCompleto == 0 &&
            fabs(state.cells[i - 1].acsr.injm.MassG + state.cells[i - 1].acsr.injm.MassP) > 1e-15) {
            double denmixSTD = (141.5 / (131.5 + state.cells[i - 1].flui.API)) * qostd1 +
                               (141.5 / (131.5 + state.cells[i - 1].acsr.injm.FluidoPro.API)) * qostd2;
            denmixSTD /= (qostd1 + qostd2);
            state.cells[i].flui.API = 141.5 / denmixSTD - 131.5;
        } else if (state.input.flashCompleto == 0 ||
                   fabs(state.cells[i - 1].acsr.injm.MassG + state.cells[i - 1].acsr.injm.MassP) <= 1e-15)
            state.cells[i].flui.API = state.cells[i - 1].flui.API;
        double waterFlowRate1;
        if ((1. - state.cells[i - 1].flui.BSW) > 0)
            waterFlowRate1 = qostd1 * state.cells[i - 1].flui.BSW / (1. - state.cells[i - 1].flui.BSW);
        else
            waterFlowRate1 = state.cells[i - 1].Mliqini * state.cells[i - 1].FW * (1. - state.cells[i - 1].bet) / (1000. * state.cells[i - 1].flui.Denag);
        double waterFlowRate2;
        if ((1. - state.cells[i - 1].acsr.injm.FluidoPro.BSW) > 0)
            waterFlowRate2 = qostd2 * state.cells[i - 1].acsr.injm.FluidoPro.BSW / (1. - state.cells[i - 1].acsr.injm.FluidoPro.BSW);
        else
            waterFlowRate2 = state.cells[i - 1].fontemassLR * fwinjl / (1000. * state.cells[i - 1].acsr.injm.FluidoPro.Denag);

        if (fabs(waterFlowRate1 + waterFlowRate2 + qostd1 + qostd2) > 1e-15 &&
            (state.cells[i - 1].acsr.injm.MassG + state.cells[i - 1].acsr.injm.MassP) > 1e-15)
            state.cells[i].flui.BSW = (waterFlowRate1 + waterFlowRate2) / (waterFlowRate1 + waterFlowRate2 + qostd1 + qostd2);
        else
            state.cells[i].flui.BSW = state.cells[i - 1].flui.BSW;
        if (fabs(waterFlowRate1 + waterFlowRate2) > 1e-15 &&
            (state.cells[i - 1].acsr.injm.MassG + state.cells[i - 1].acsr.injm.MassP) > 1e-15)
            state.cells[i].flui.Denag = (waterFlowRate1 * state.cells[i - 1].flui.Denag +
                                    waterFlowRate2 * state.cells[i - 1].acsr.injm.FluidoPro.Denag) /
                                   (waterFlowRate1 + waterFlowRate2);
        else
            state.cells[i].flui.Denag = state.cells[i - 1].flui.Denag;
        // reavaliacao dos pares necessarios para o modelo de viscosidade de oleo morto
        // ASTM. Isto ÃƒÂ© feito tendo sempre os mesmos valores de temperatura maxima e
        // e minima para a construcao dos pares
        if (fabs(qostd1 + qostd2) > 1e-15 && state.input.flashCompleto == 0 &&
            (state.cells[i - 1].acsr.injm.MassG + state.cells[i - 1].acsr.injm.MassP) > 1e-15) { // vazao de liquido >0
            state.cells[i].flui.TempL = temperatureLow;
            state.cells[i].flui.TempH = temperatureHigh;
            state.cells[i].flui.LVisL = (qostd1 * state.cells[i - 1].flui.VisOM(temperatureLow) + qostd2 * state.cells[i - 1].acsr.injm.FluidoPro.VisOM(temperatureLow)) / (qostd1 + qostd2);
            state.cells[i].flui.LVisH = (qostd1 * state.cells[i - 1].flui.VisOM(temperatureHigh) + qostd2 * state.cells[i - 1].acsr.injm.FluidoPro.VisOM(temperatureHigh)) / (qostd1 + qostd2);
        } else if (state.input.flashCompleto == 0 ||
                   fabs(state.cells[i - 1].acsr.injm.MassG + state.cells[i - 1].acsr.injm.MassP) <= 1e-15) { // se,m vazao de liquido
            state.cells[i].flui.TempL = state.cells[i - 1].flui.TempL;
            state.cells[i].flui.TempH = state.cells[i - 1].flui.TempH;
            state.cells[i].flui.LVisL = state.cells[i - 1].flui.LVisL;
            state.cells[i].flui.LVisH = state.cells[i - 1].flui.LVisH;
        }
    }
    if (state.input.flashCompleto == 0) {
        state.cells[i].flui.RenovaFluido();
        correctGasSpecificGravity(state, i);
    }
    if ((state.cells[i - 1].bet > (*state.globals).localtiny * 1e-6 || fabs(state.cells[i - 1].fontemassCR) > (*state.globals).localtiny) &&
        fabs(state.cells[i - 1].acsr.injm.MassG + state.cells[i - 1].acsr.injm.MassP + state.cells[i - 1].acsr.injm.MassC) > 1e-15) {
        // reavaliacao da fracao volumetrica do liquido complementar
        // observar que a fonte de liquido pode ter uma fracao de liquido complementar
        // distinta da onservada a esquerda da celula i-1
        double oilVolumeFactorNeighbour = state.cells[i].flui.BOFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        double waterVolumeFactorNeighbour = state.cells[i].flui.BAFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        double waterCutNeighbour = state.cells[i].flui.BSW * waterVolumeFactorNeighbour / (oilVolumeFactorNeighbour + waterVolumeFactorNeighbour * state.cells[i].flui.BSW - state.cells[i].flui.BSW * oilVolumeFactorNeighbour);
        oilVolumeFactorInSitu = state.cells[i].flui.BOFunc(state.cells[i].presaux, tmed);
        waterVolumeFactorInSitu = state.cells[i].flui.BAFunc(state.cells[i].presaux, tmed);
        waterCutInSitu = state.cells[i].flui.BSW * waterVolumeFactorInSitu / (oilVolumeFactorInSitu + waterVolumeFactorInSitu * state.cells[i].flui.BSW - state.cells[i].flui.BSW * oilVolumeFactorInSitu);
        double qlpF = state.cells[i - 1].fontemassLR / state.cells[i - 1].acsr.injm.FluidoPro.MasEspLiq(state.cells[i - 1].pres, state.cells[i - 1].temp);
        double qlcF = state.cells[i - 1].fontemassCR / state.cells[i - 1].acsr.injm.fluidocol.MasEspFlu(state.cells[i - 1].pres, state.cells[i - 1].temp);
        double betN = (state.cells[i - 1].QL * state.cells[i - 1].bet + qlcF) / (state.cells[i - 1].QL + qlpF + qlcF);
        double oilFlowRate = (state.cells[i - 1].QL + qlpF + qlcF) * (1 - waterCutNeighbour) * (1 - betN) * oilVolumeFactorInSitu / oilVolumeFactor;
        double waterFlowRate;
        if (waterCutInSitu < (1 - (*state.globals).localtiny))
            waterFlowRate = waterCutInSitu * oilFlowRate / (1 - waterCutInSitu);
        else
            waterFlowRate = state.cells[i - 1].QL * (1 - state.cells[i - 1].bet) + qlpF;
        double completionFlowRate = (state.cells[i - 1].QL * (state.cells[i - 1].bet) + qlcF) * state.cells[i - 1].fluicol.MasEspFlu(state.cells[i - 1].pres, state.cells[i - 1].temp) / state.cells[i].fluicol.MasEspFlu(state.cells[i].presaux, tmed);
        state.cells[i].bet = fabs(completionFlowRate) / (fabs(oilFlowRate) + fabs(waterFlowRate) + fabs(completionFlowRate));
    } else
        state.cells[i].bet = state.cells[i - 1].bet;
}

void applySteadyMassLeakSource(const SteadyStateState &state, int i, int mudaRGO, double temperatureLow, double temperatureHigh, double oilVolumeFactor, double waterVolumeFactor, double solutionGasRatio, double tmed, double &residenceTimeSource, double &oilVolumeFactorInSitu, double &waterVolumeFactorInSitu, double &waterCutInSitu) {
    ProFlu fluF;
    // define qual o fluido envolvido na fonte, se a pressÃƒÂ£o ambiente for menor do que a
    // pressao da tubulacao, fluido da tubulacao, senao, fluido definido como
    // fluido ambiente
    if (state.cells[i - 1].acsr.fontechk.presT > state.cells[i - 1].acsr.fontechk.pamb) {
        fluF = state.cells[i - 1].acsr.fontechk.fluidoP;
    } else {
        fluF = state.cells[i - 1].acsr.fontechk.fluidoPamb;
    }
    double apiGravity = state.cells[i - 1].flui.API;
    double rhololeo = 1000. * 141.5 / (131.5 + apiGravity);
    double rhololeoF = 1000. * 141.5 / (131.5 + fluF.API);
    double rholis = state.cells[i - 1].flui.MasEspLiq(state.cells[i - 1].pres, state.cells[i - 1].temp);
    double rholisF = fluF.MasEspLiq(state.cells[i - 1].pres, state.cells[i - 1].temp);
    double boinjl = fluF.BOFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
    double bainjl = fluF.BAFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
    double fwinjl = fluF.BSW * bainjl / (boinjl + bainjl * fluF.BSW - fluF.BSW * boinjl);
    // vazao de oleo standard antes da fonte:
    double qostd1 = 0.;
    if (state.cells[i - 1].flui.dStockTankVaporMassFraction < 1. - 1e-15)
        qostd1 = state.cells[i - 1].Mliqini * (1. - state.cells[i - 1].FW) * (1. - state.cells[i - 1].bet) / (oilVolumeFactor * rholis);
    // vazao de oleo standard da fonte:
    double qostd2 = 0.;
    if (fluF.dStockTankVaporMassFraction < 1. - 1e-15)
        qostd2 = state.cells[i - 1].fontemassLR * (1. - fwinjl) / (boinjl * rholisF);
    // vazao total de gas nas condicoes standard fruto da soma do gas transportado
    // da fronteira a esquerda
    // da celula e do gas associado ao vazamento
    double qgstd;
    if (fabs(qostd2) > 0.)
        qgstd = qostd1 * state.cells[i - 1].flui.RGO + qostd2 * fluF.RGO;
    else {
        qgstd = qostd1 * state.cells[i - 1].flui.RGO + state.cells[i - 1].fontemassGR / (fluF.Deng * 1.225);
    }
    double liquidHoldup = 1 - state.cells[i - 1].alf;
    double completionFraction = state.cells[i - 1].bet;
    double bsw = state.cells[i - 1].FW;
    double rhog = state.cells[i - 1].flui.MasEspGas(state.cells[i - 1].pres, state.cells[i - 1].temp);
    double rhogST = state.cells[i - 1].flui.Deng * 1.225;
    // calculo do volume de leve na celula i, seguindo o equacionamento mostrado em relatorio,
    // nÃƒÂ£o ÃƒÂ© relevante para o permanente, mas deve ser calculado,
    // pois ÃƒÂ© utilizado como entrada no transiente
    state.cells[i - 1].VolLeveST = (((1 - liquidHoldup) * rhog / rhogST) + liquidHoldup * (1 - completionFraction) * (1. - bsw) * solutionGasRatio / oilVolumeFactor);
    if (state.cells[i - 1].VolLeveST < 1e-15)
        state.cells[i - 1].VolLeveST = 0.;

    if (state.input.trackRGO == -1) {
        // esta chave nÃƒÂ£o ÃƒÂ© utilizada, ÃƒÂ© mantida aqui como reserva, atualmente este
        // calculo nÃƒÂ£o esta funcionando, lembrando que arq.trackRGO assume apenas 2 valores, 0 ou 1
        if (liquidHoldup > (*state.globals).localtiny && completionFraction < (1. - (*state.globals).localtiny * 1e-6) && bsw < (1. - (*state.globals).localtiny) && mudaRGO == 1)
            state.cells[i].flui.RGO = state.cells[i - 1].VolLeveST * oilVolumeFactor / (liquidHoldup * (1 - completionFraction) * (1 - bsw));
        if (state.cells[i].flui.RGO > (*state.globals).RGOMax && mudaRGO == 1)
            state.cells[i].flui.RGO = (*state.globals).RGOMax;
    } else { // calculo do novo RGO do separador modificado pela IPR
        if (fabs(qostd1 + qostd2) > (*state.globals).localtiny * 1e-6 && mudaRGO == 1)
            state.cells[i].flui.RGO = qgstd / (qostd1 + qostd2);
        else if (mudaRGO == 1)
            state.cells[i].flui.RGO = (*state.globals).RGOMax;
    }
    // calculo da nova densidade de gas nas condicoes standard modificada pelo vazamento
    if (fabs(qgstd) > (*state.globals).localtiny && state.input.flashCompleto == 0) {
        if (fabs(qostd2) > 0)
            state.cells[i].flui.Deng = (qostd1 * state.cells[i - 1].flui.RGO * state.cells[i - 1].flui.Deng + qostd2 * fluF.RGO * fluF.Deng) / qgstd;
        else
            state.cells[i].flui.Deng = (qostd1 * state.cells[i - 1].flui.RGO * state.cells[i - 1].flui.Deng + state.cells[i - 1].fontemassGR / (fluF.Deng * 1.225)) / qgstd;
    }
    // calculo da nova fracao de co2 modificada pelo vazamento
    if (fabs(qgstd) > (*state.globals).localtiny && state.input.flashCompleto == 0) {
        if (fabs(qostd2) > 0)
            state.cells[i].flui.yco2 = (qostd1 * state.cells[i - 1].flui.RGO * state.cells[i - 1].flui.yco2 + qostd2 * fluF.RGO * fluF.yco2) / qgstd;
        else
            state.cells[i].flui.yco2 = (qostd1 * state.cells[i - 1].flui.RGO * state.cells[i - 1].flui.yco2 + fluF.yco2 * state.cells[i - 1].fontemassGR / (fluF.Deng * 1.225)) / qgstd;
    }
    if (state.productionFluidCount > 1 || (*state.globals).chaverede == 1) {
        // o modelo seja ASTM. Observar que isto sÃƒÂ³ faz sentido se se tiver mais de um fluido
        // de producao cadastrado no JSON, ou a simulacao se insere em um sistema de redes
        // com varios tramos alimentando outros tramos com fluidos com propriedades diferentes
        // entre cada tramo
        if (fabs(qostd1 + qostd2) > 1e-15 && state.input.flashCompleto == 0) {
            double denmixSTD = (141.5 / (131.5 + state.cells[i - 1].flui.API)) * qostd1 +
                               (141.5 / (131.5 + fluF.API)) * qostd2;
            denmixSTD /= (qostd1 + qostd2);
            state.cells[i].flui.API = 141.5 / denmixSTD - 131.5;
        } else if (state.input.flashCompleto == 0)
            state.cells[i].flui.API = state.cells[i - 1].flui.API;
        double waterFlowRate1;
        if ((1. - state.cells[i - 1].flui.BSW) > 0)
            waterFlowRate1 = qostd1 * state.cells[i - 1].flui.BSW / (1. - state.cells[i - 1].flui.BSW);
        else
            waterFlowRate1 = state.cells[i - 1].Mliqini * state.cells[i - 1].FW * (1. - state.cells[i - 1].bet) / (1000. * state.cells[i - 1].flui.Denag);
        double waterFlowRate2;
        if ((1. - fluF.BSW) > 0)
            waterFlowRate2 = qostd2 * fluF.BSW / (1. - fluF.BSW);
        else
            waterFlowRate2 = state.cells[i - 1].fontemassLR * fwinjl / (1000. * fluF.Denag);

        if (fabs(waterFlowRate1 + waterFlowRate2 + qostd1 + qostd2) > 1e-15)
            state.cells[i].flui.BSW = (waterFlowRate1 + waterFlowRate2) / (waterFlowRate1 + waterFlowRate2 + qostd1 + qostd2);
        else
            state.cells[i].flui.BSW = state.cells[i - 1].flui.BSW;

        if (fabs(qostd1 + qostd2) > 1e-15 && state.input.flashCompleto == 0) {
            // reavaliacao dos pares necessarios para o modelo de viscosidade de oleo morto
            // ASTM. Isto ÃƒÂ© feito tendo sempre os mesmos valores de temperatura maxima e
            // e minima para a construcao dos pares
            state.cells[i].flui.TempL = temperatureLow;
            state.cells[i].flui.TempH = temperatureHigh;
            state.cells[i].flui.LVisL = (qostd1 * state.cells[i - 1].flui.VisOM(temperatureLow) + qostd2 * fluF.VisOM(temperatureLow)) / (qostd1 + qostd2);
            state.cells[i].flui.LVisH = (qostd1 * state.cells[i - 1].flui.VisOM(temperatureHigh) + qostd2 * fluF.VisOM(temperatureHigh)) / (qostd1 + qostd2);
        } else if (state.input.flashCompleto == 0) {
            state.cells[i].flui.TempL = state.cells[i - 1].flui.TempL;
            state.cells[i].flui.TempH = state.cells[i - 1].flui.TempH;
            state.cells[i].flui.LVisL = state.cells[i - 1].flui.LVisL;
            state.cells[i].flui.LVisH = state.cells[i - 1].flui.LVisH;
        }
    }
    if (state.input.flashCompleto == 0) {
        state.cells[i].flui.RenovaFluido();
        correctGasSpecificGravity(state, i);
    }
    if (state.cells[i - 1].bet > (*state.globals).localtiny * 1e-6 || fabs(state.cells[i - 1].fontemassCR) > (*state.globals).localtiny * 1e-6) {
        // reavaliacao da fracao volumetrica do liquido complementar
        // observar a fonte de vazamento pode ter uma fracao de liquido complementar
        // distinta da observada a esquerda da celula i-1
        double oilVolumeFactorNeighbour = state.cells[i].flui.BOFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        double waterVolumeFactorNeighbour = state.cells[i].flui.BAFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        double waterCutNeighbour = state.cells[i].flui.BSW * waterVolumeFactorNeighbour / (oilVolumeFactorNeighbour + waterVolumeFactorNeighbour * state.cells[i].flui.BSW - state.cells[i].flui.BSW * oilVolumeFactorNeighbour);
        oilVolumeFactorInSitu = state.cells[i].flui.BOFunc(state.cells[i].presaux, tmed);
        waterVolumeFactorInSitu = state.cells[i].flui.BAFunc(state.cells[i].presaux, tmed);
        waterCutInSitu = state.cells[i].flui.BSW * waterVolumeFactorInSitu / (oilVolumeFactorInSitu + waterVolumeFactorInSitu * state.cells[i].flui.BSW - state.cells[i].flui.BSW * oilVolumeFactorInSitu);
        double qlpF = state.cells[i - 1].fontemassLR / fluF.MasEspLiq(state.cells[i - 1].pres, state.cells[i - 1].temp);
        double qlcF = state.cells[i - 1].fontemassCR / state.cells[i - 1].acsr.fontechk.fluidocol.MasEspFlu(state.cells[i - 1].pres, state.cells[i - 1].temp);
        double betN = (state.cells[i - 1].QL * state.cells[i - 1].bet + qlcF) / (state.cells[i - 1].QL + qlpF + qlcF);
        double oilFlowRate = (state.cells[i - 1].QL + qlpF + qlcF) * (1 - waterCutNeighbour) * (1 - betN) * oilVolumeFactorInSitu / oilVolumeFactor;
        double waterFlowRate;
        if (waterCutInSitu < (1 - (*state.globals).localtiny))
            waterFlowRate = waterCutInSitu * oilFlowRate / (1 - waterCutInSitu);
        else
            waterFlowRate = state.cells[i - 1].QL * (1 - state.cells[i - 1].bet) + qlpF;
        double completionFlowRate = (state.cells[i - 1].QL * (state.cells[i - 1].bet) + qlcF) * state.cells[i - 1].fluicol.MasEspFlu(state.cells[i - 1].pres, state.cells[i - 1].temp) / state.cells[i].fluicol.MasEspFlu(state.cells[i].presaux, tmed);
        state.cells[i].bet = fabs(completionFlowRate) / (fabs(oilFlowRate) + fabs(waterFlowRate) + fabs(completionFlowRate));
    } else
        state.cells[i].bet = state.cells[i - 1].bet;
}

void applySteadyMassRadialPorous(const SteadyStateState &state, int i, int mudaRGO, double temperatureLow, double temperatureHigh, double oilVolumeFactor, double waterVolumeFactor, double solutionGasRatio, double tmed, double &residenceTimeSource, double &oilVolumeFactorInSitu, double &waterVolumeFactorInSitu, double &waterCutInSitu) {
    double apiGravity = state.cells[i - 1].flui.API;
    double rhololeo = 1000. * 141.5 / (131.5 + apiGravity);
    double rhololeoF = 1000. * 141.5 / (131.5 + state.cells[i - 1].acsr.radialPoro.flup.API);
    double rholis = state.cells[i - 1].flui.MasEspLiq(state.cells[i - 1].pres, state.cells[i - 1].temp);
    double rholisF;
    double boipr;
    double baipr;
    double fwipr;
    residenceTimeSource = 0.;
    if (fabs(state.cells[i - 1].fontemassLR) > 1e-15) {
        rholisF = state.cells[i - 1].acsr.radialPoro.flup.MasEspLiq(state.cells[i - 1].pres, state.cells[i - 1].temp);
        boipr = state.cells[i - 1].acsr.radialPoro.flup.BOFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        baipr = state.cells[i - 1].acsr.radialPoro.flup.BAFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        fwipr = state.cells[i - 1].acsr.radialPoro.BSW * baipr / (boipr + baipr * state.cells[i - 1].acsr.radialPoro.BSW - state.cells[i - 1].acsr.radialPoro.BSW * boipr);
    } else {
        rholisF = rholis;
        boipr = oilVolumeFactor;
        baipr = waterVolumeFactor;
        fwipr = state.cells[i - 1].FW;
    }
    // vazao de oleo standard antes da fonte:
    double qostd1 = 0.;
    if (state.cells[i - 1].flui.dStockTankVaporMassFraction < 1. - 1e-15)
        qostd1 = state.cells[i - 1].Mliqini * (1. - state.cells[i - 1].FW) * (1. - state.cells[i - 1].bet) / (oilVolumeFactor * rholis);
    // vazao de oleo standard da fonte:
    double qostd2 = 0.;
    if (state.cells[i - 1].acsr.radialPoro.flup.dStockTankVaporMassFraction < 1. - 1e-15)
        qostd2 = state.cells[i - 1].fontemassLR * (1. - fwipr) / (boipr * rholisF);
    // vazao total de gas nas condicoes standard fruto da soma do gas transportado
    // da fronteira a esquerda
    // da celula e do gas associado da IPR
    double qgstd = qostd1 * state.cells[i - 1].flui.RGO + qostd2 * state.cells[i - 1].acsr.radialPoro.flup.RGO;
    double liquidHoldup = 1 - state.cells[i - 1].alf;
    double completionFraction = state.cells[i - 1].bet;
    double bsw = state.cells[i - 1].FW;
    double rhog = state.cells[i - 1].flui.MasEspGas(state.cells[i - 1].pres, state.cells[i - 1].temp);
    double rhogST = state.cells[i - 1].flui.Deng * 1.225;
    // calculo do volume de leve na celula i, seguindo o equacionamento mostrado em relatorio,
    // nÃƒÂ£o ÃƒÂ© relevante para o permanente, mas deve ser calculado,
    // pois ÃƒÂ© utilizado como entrada no transiente
    state.cells[i - 1].VolLeveST = (((1 - liquidHoldup) * rhog / rhogST) + liquidHoldup * (1 - completionFraction) * (1. - bsw) * solutionGasRatio / oilVolumeFactor);
    if (state.cells[i - 1].VolLeveST < 1e-15)
        state.cells[i - 1].VolLeveST = 0.;
    if (state.input.trackRGO == -1) {
        // esta chave nÃƒÂ£o ÃƒÂ© utilizada, ÃƒÂ© mantida aqui como reserva, atualmente este
        // calculo nÃƒÂ£o esta funcionando, lembrando que arq.trackRGO assume apenas 2 valores, 0 ou 1
        if (fabs(state.cells[i - 1].fontemassLR) > 1e-15) {
            if (liquidHoldup > (*state.globals).localtiny && completionFraction < (1. - (*state.globals).localtiny) && bsw < (1. - (*state.globals).localtiny) &&
                mudaRGO == 1 && fabs(state.cells[i - 1].fontemassLR) > 1e-15)
                state.cells[i].flui.RGO = state.cells[i - 1].VolLeveST * oilVolumeFactor / (liquidHoldup * (1 - completionFraction) * (1 - bsw));
            if (state.cells[i].flui.RGO > (*state.globals).RGOMax && mudaRGO == 1)
                state.cells[i].flui.RGO = (*state.globals).RGOMax;
        } else
            state.cells[i].flui.RGO = state.cells[i - 1].flui.RGO;

    } else { // calculo do novo RGO do separador modificado pela IPR
        if ((state.cells[i - 1].fontemassLR) > 1e-15 && (state.cells[i - 1].fontemassGR) > 1e-15) {
            if (fabs(qostd1 + qostd2) > (*state.globals).localtiny * 1e-10 && mudaRGO == 1)
                state.cells[i].flui.RGO = qgstd / (qostd1 + qostd2);
            else if (mudaRGO == 1)
                state.cells[i].flui.RGO = (*state.globals).RGOMax;
        } else
            state.cells[i].flui.RGO = state.cells[i - 1].flui.RGO;
    }
    // calculo da nova densidade de gas nas condicoes standard modificada pela IPR
    if (state.input.flashCompleto == 0 && (state.cells[i - 1].fontemassLR) > 1e-15 && (state.cells[i - 1].fontemassGR) > 1e-15 &&
        fabs(qgstd) > 1e-15)
        state.cells[i].flui.Deng = (qostd1 * state.cells[i - 1].flui.RGO * state.cells[i - 1].flui.Deng + qostd2 * state.cells[i - 1].acsr.radialPoro.flup.RGO * state.cells[i - 1].acsr.radialPoro.flup.Deng) / qgstd;
    else
        state.cells[i].flui.Deng = state.cells[i - 1].flui.Deng;
    // calculo da nova fracao de co2 modificada pela fonte de liquido
    if (state.input.flashCompleto == 0 && (state.cells[i - 1].fontemassLR) > 1e-15 && (state.cells[i - 1].fontemassGR) > 1e-15 && fabs(qgstd) > 1e-15)
        state.cells[i].flui.yco2 = (qostd1 * state.cells[i - 1].flui.RGO * state.cells[i - 1].flui.yco2 + qostd2 * state.cells[i - 1].acsr.radialPoro.flup.RGO * state.cells[i - 1].acsr.radialPoro.flup.yco2) / qgstd;
    else
        state.cells[i].flui.yco2 = state.cells[i - 1].flui.yco2;
    if (state.productionFluidCount > 1 || (*state.globals).chaverede == 1) {
        // o modelo seja ASTM. Observar que isto sÃƒÂ³ faz sentido se se tiver mais de um fluido
        // de producao cadastrado no JSON, ou a simulacao se insere em um sistema de redes
        // com varios tramos alimentando outros tramos com fluidos com propriedades diferentes
        // entre cada tramo
        if (fabs(qostd1 + qostd2) > 1e-15 && state.input.flashCompleto == 0 &&
            (state.cells[i - 1].fontemassLR) > 1e-15) {
            double denmixSTD = (141.5 / (131.5 + state.cells[i - 1].flui.API)) * qostd1 +
                               (141.5 / (131.5 + state.cells[i - 1].acsr.radialPoro.flup.API)) * qostd2;
            denmixSTD /= (qostd1 + qostd2);
            state.cells[i].flui.API = 141.5 / denmixSTD - 131.5;
        } else if (state.input.flashCompleto == 0 || fabs(state.cells[i - 1].fontemassLR) <= 1e-15)
            state.cells[i].flui.API = state.cells[i - 1].flui.API;
        double waterFlowRate1;
        if ((1. - state.cells[i - 1].flui.BSW) > 0)
            waterFlowRate1 = qostd1 * state.cells[i - 1].flui.BSW / (1. - state.cells[i - 1].flui.BSW);
        else
            waterFlowRate1 = state.cells[i - 1].Mliqini * state.cells[i - 1].FW * (1. - state.cells[i - 1].bet) / (1000. * state.cells[i - 1].flui.Denag);
        double waterFlowRate2;
        if ((1. - state.cells[i - 1].acsr.radialPoro.BSW) > 0)
            waterFlowRate2 = qostd2 * state.cells[i - 1].acsr.radialPoro.BSW / (1. - state.cells[i - 1].acsr.radialPoro.BSW);
        else
            waterFlowRate2 = state.cells[i - 1].fontemassLR * fwipr / (1000. * state.cells[i - 1].acsr.radialPoro.flup.Denag);

        if (fabs(waterFlowRate1 + waterFlowRate2 + qostd1 + qostd2) > 1e-15 &&
            fabs(state.cells[i - 1].fontemassLR) > 1e-15)
            state.cells[i].flui.BSW = ((waterFlowRate1 + qostd1) * state.cells[i - 1].flui.BSW + (waterFlowRate2 + qostd2) * state.cells[i - 1].acsr.radialPoro.BSW) / (waterFlowRate1 + waterFlowRate2 + qostd1 + qostd2);
        else
            state.cells[i].flui.BSW = state.cells[i - 1].flui.BSW;
        if ((waterFlowRate1 + waterFlowRate2) > 1e-15 && (state.cells[i - 1].fontemassLR) > 1e-15)
            state.cells[i].flui.Denag = (waterFlowRate1 * state.cells[i - 1].flui.Denag +
                                    waterFlowRate2 * state.cells[i - 1].acsr.radialPoro.flup.Denag) /
                                   (waterFlowRate1 + waterFlowRate2);
        else
            state.cells[i].flui.Denag = state.cells[i - 1].flui.Denag;
        if (fabs(qostd1 + qostd2) > 1e-15 && state.input.flashCompleto == 0 && (state.cells[i - 1].fontemassLR) > 1e-15) {
            // reavaliacao dos pares necessarios para o modelo de viscosidade de oleo morto
            // ASTM. Isto ÃƒÂ© feito tendo sempre os mesmos valores de temperatura maxima e
            // e minima para a construcao dos pares
            state.cells[i].flui.TempL = temperatureLow;
            state.cells[i].flui.TempH = temperatureHigh;
            state.cells[i].flui.LVisL = (qostd1 * state.cells[i - 1].flui.VisOM(temperatureLow) + qostd2 * state.cells[i - 1].acsr.radialPoro.flup.VisOM(temperatureLow)) / (qostd1 + qostd2);
            state.cells[i].flui.LVisH = (qostd1 * state.cells[i - 1].flui.VisOM(temperatureHigh) + qostd2 * state.cells[i - 1].acsr.radialPoro.flup.VisOM(temperatureHigh)) / (qostd1 + qostd2);
        } else if (state.input.flashCompleto == 0 || fabs(state.cells[i - 1].fontemassLR) <= 1e-15) { // sem vazao de liquido
            state.cells[i].flui.TempL = state.cells[i - 1].flui.TempL;
            state.cells[i].flui.TempH = state.cells[i - 1].flui.TempH;
            state.cells[i].flui.LVisL = state.cells[i - 1].flui.LVisL;
            state.cells[i].flui.LVisH = state.cells[i - 1].flui.LVisH;
        }
    }
    if (state.input.flashCompleto == 0) {
        state.cells[i].flui.RenovaFluido();
        correctGasSpecificGravity(state, i);
    }
    if (state.cells[i - 1].bet > (*state.globals).localtiny * 1e-6 && fabs(state.cells[i - 1].fontemassLR) > 1e-15) {
        // reavaliacao da fracao volumetrica do liquido complementar
        // observar que a IPR pode ter uma fracao de liquido complementar
        // distinta da onservada a esquerda da celula i-1
        double oilVolumeFactorNeighbour = state.cells[i].flui.BOFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        double waterVolumeFactorNeighbour = state.cells[i].flui.BAFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        double waterCutNeighbour = state.cells[i].flui.BSW * waterVolumeFactorNeighbour / (oilVolumeFactorNeighbour + waterVolumeFactorNeighbour * state.cells[i].flui.BSW - state.cells[i].flui.BSW * oilVolumeFactorNeighbour);
        oilVolumeFactorInSitu = state.cells[i].flui.BOFunc(state.cells[i].presaux, tmed);
        waterVolumeFactorInSitu = state.cells[i].flui.BAFunc(state.cells[i].presaux, tmed);
        waterCutInSitu = state.cells[i].flui.BSW * waterVolumeFactorInSitu / (oilVolumeFactorInSitu + waterVolumeFactorInSitu * state.cells[i].flui.BSW - state.cells[i].flui.BSW * oilVolumeFactorInSitu);
        double qlpF = state.cells[i - 1].fontemassLR / state.cells[i - 1].acsr.radialPoro.flup.MasEspLiq(state.cells[i - 1].pres, state.cells[i - 1].temp);
        double qlcF = 0.;
        double betN = (state.cells[i - 1].QL * state.cells[i - 1].bet + qlcF) / (state.cells[i - 1].QL + qlpF + qlcF);
        double oilFlowRate = (state.cells[i - 1].QL + qlpF + qlcF) * (1 - waterCutNeighbour) * (1 - betN) * oilVolumeFactorInSitu / oilVolumeFactor;
        double waterFlowRate;
        if (waterCutInSitu < (1 - (*state.globals).localtiny))
            waterFlowRate = waterCutInSitu * oilFlowRate / (1 - waterCutInSitu);
        else
            waterFlowRate = state.cells[i - 1].QL * (1 - state.cells[i - 1].bet) + qlpF;
        double completionFlowRate = (state.cells[i - 1].QL * (state.cells[i - 1].bet) + qlcF) * state.cells[i - 1].fluicol.MasEspFlu(state.cells[i - 1].pres, state.cells[i - 1].temp) / state.cells[i].fluicol.MasEspFlu(state.cells[i].presaux, tmed);
        state.cells[i].bet = fabs(completionFlowRate) / (fabs(oilFlowRate) + fabs(waterFlowRate) + fabs(completionFlowRate));
    } else
        state.cells[i].bet = state.cells[i - 1].bet;
}

void applySteadyMassPorous2D(const SteadyStateState &state, int i, int mudaRGO, double temperatureLow, double temperatureHigh, double oilVolumeFactor, double waterVolumeFactor, double solutionGasRatio, double tmed, double &residenceTimeSource, double &oilVolumeFactorInSitu, double &waterVolumeFactorInSitu, double &waterCutInSitu) {
    double apiGravity = state.cells[i - 1].flui.API;
    double rhololeo = 1000. * 141.5 / (131.5 + apiGravity);
    double rhololeoF = 1000. * 141.5 / (131.5 + state.cells[i - 1].acsr.poroso2D.dados.flup.API);
    double rholis = state.cells[i - 1].flui.MasEspLiq(state.cells[i - 1].pres, state.cells[i - 1].temp);
    double rholisF;
    double boipr;
    double baipr;
    double fwipr;
    residenceTimeSource = 0.;
    if (fabs(state.cells[i - 1].fontemassLR) > 1e-15) {
        rholisF = state.cells[i - 1].acsr.poroso2D.dados.flup.MasEspLiq(state.cells[i - 1].pres, state.cells[i - 1].temp);
        boipr = state.cells[i - 1].acsr.poroso2D.dados.flup.BOFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        baipr = state.cells[i - 1].acsr.poroso2D.dados.flup.BAFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        fwipr = state.cells[i - 1].acsr.poroso2D.dados.transfer.BSW * baipr / (boipr + baipr * state.cells[i - 1].acsr.poroso2D.dados.transfer.BSW - state.cells[i - 1].acsr.poroso2D.dados.transfer.BSW * boipr);
    } else {
        rholisF = rholis;
        boipr = oilVolumeFactor;
        baipr = waterVolumeFactor;
        fwipr = state.cells[i - 1].FW;
    }
    // vazao de oleo standard antes da fonte:
    double qostd1 = 0.;
    if (state.cells[i - 1].flui.dStockTankVaporMassFraction < 1. - 1e-15)
        qostd1 = state.cells[i - 1].Mliqini * (1. - state.cells[i - 1].FW) * (1. - state.cells[i - 1].bet) / (oilVolumeFactor * rholis);
    // vazao de oleo standard da fonte:
    double qostd2 = 0.;
    if (state.cells[i - 1].acsr.poroso2D.dados.flup.dStockTankVaporMassFraction < 1. - 1e-15)
        qostd2 = state.cells[i - 1].fontemassLR * (1. - fwipr) / (boipr * rholisF);
    // vazao total de gas nas condicoes standard fruto da soma do gas transportado
    // da fronteira a esquerda
    // da celula e do gas associado da IPR
    double qgstd = qostd1 * state.cells[i - 1].flui.RGO + qostd2 * state.cells[i - 1].acsr.poroso2D.dados.flup.RGO;
    double liquidHoldup = 1 - state.cells[i - 1].alf;
    double completionFraction = state.cells[i - 1].bet;
    double bsw = state.cells[i - 1].FW;
    double rhog = state.cells[i - 1].flui.MasEspGas(state.cells[i - 1].pres, state.cells[i - 1].temp);
    double rhogST = state.cells[i - 1].flui.Deng * 1.225;
    // calculo do volume de leve na celula i, seguindo o equacionamento mostrado em relatorio,
    // nÃƒÂ£o ÃƒÂ© relevante para o permanente, mas deve ser calculado,
    // pois ÃƒÂ© utilizado como entrada no transiente
    state.cells[i - 1].VolLeveST = (((1 - liquidHoldup) * rhog / rhogST) + liquidHoldup * (1 - completionFraction) * (1. - bsw) * solutionGasRatio / oilVolumeFactor);
    if (state.cells[i - 1].VolLeveST < 1e-15)
        state.cells[i - 1].VolLeveST = 0.;
    if (state.input.trackRGO == -1) {
        // esta chave nÃƒÂ£o ÃƒÂ© utilizada, ÃƒÂ© mantida aqui como reserva, atualmente este
        // calculo nÃƒÂ£o esta funcionando, lembrando que arq.trackRGO assume apenas 2 valores, 0 ou 1
        if (fabs(state.cells[i - 1].fontemassLR) > 1e-15) {
            if (liquidHoldup > (*state.globals).localtiny && completionFraction < (1. - (*state.globals).localtiny) && bsw < (1. - (*state.globals).localtiny) &&
                mudaRGO == 1 && fabs(state.cells[i - 1].fontemassLR) > 1e-15)
                state.cells[i].flui.RGO = state.cells[i - 1].VolLeveST * oilVolumeFactor / (liquidHoldup * (1 - completionFraction) * (1 - bsw));
            if (state.cells[i].flui.RGO > (*state.globals).RGOMax && mudaRGO == 1)
                state.cells[i].flui.RGO = (*state.globals).RGOMax;
        } else
            state.cells[i].flui.RGO = state.cells[i - 1].flui.RGO;

    } else { // calculo do novo RGO do separador modificado pela IPR
        if ((state.cells[i - 1].fontemassLR) > 1e-15 && (state.cells[i - 1].fontemassGR) > 1e-15) {
            if (fabs(qostd1 + qostd2) > (*state.globals).localtiny * 1e-10 && mudaRGO == 1)
                state.cells[i].flui.RGO = qgstd / (qostd1 + qostd2);
            else if (mudaRGO == 1)
                state.cells[i].flui.RGO = (*state.globals).RGOMax;
        } else
            state.cells[i].flui.RGO = state.cells[i - 1].flui.RGO;
    }
    // calculo da nova densidade de gas nas condicoes standard modificada pela IPR
    if (state.input.flashCompleto == 0 && (state.cells[i - 1].fontemassLR) > 1e-15 && (state.cells[i - 1].fontemassGR) > 1e-15 &&
        fabs(qgstd) > 1e-15)
        state.cells[i].flui.Deng = (qostd1 * state.cells[i - 1].flui.RGO * state.cells[i - 1].flui.Deng + qostd2 * state.cells[i - 1].acsr.poroso2D.dados.flup.RGO * state.cells[i - 1].acsr.poroso2D.dados.flup.Deng) / qgstd;
    else
        state.cells[i].flui.Deng = state.cells[i - 1].flui.Deng;
    // calculo da nova fracao de co2 modificada pela fonte de liquido
    if (state.input.flashCompleto == 0 && (state.cells[i - 1].fontemassLR) > 1e-15 && (state.cells[i - 1].fontemassGR) > 1e-15 && fabs(qgstd) > 1e-15)
        state.cells[i].flui.yco2 = (qostd1 * state.cells[i - 1].flui.RGO * state.cells[i - 1].flui.yco2 + qostd2 * state.cells[i - 1].acsr.poroso2D.dados.flup.RGO * state.cells[i - 1].acsr.poroso2D.dados.flup.yco2) / qgstd;
    else
        state.cells[i].flui.yco2 = state.cells[i - 1].flui.yco2;
    if (state.productionFluidCount > 1 || (*state.globals).chaverede == 1) {
        // o modelo seja ASTM. Observar que isto sÃƒÂ³ faz sentido se se tiver mais de um fluido
        // de producao cadastrado no JSON, ou a simulacao se insere em um sistema de redes
        // com varios tramos alimentando outros tramos com fluidos com propriedades diferentes
        // entre cada tramo
        if (fabs(qostd1 + qostd2) > 1e-15 && state.input.flashCompleto == 0 &&
            (state.cells[i - 1].fontemassLR) > 1e-15) {
            double denmixSTD = (141.5 / (131.5 + state.cells[i - 1].flui.API)) * qostd1 +
                               (141.5 / (131.5 + state.cells[i - 1].acsr.poroso2D.dados.flup.API)) * qostd2;
            denmixSTD /= (qostd1 + qostd2);
            state.cells[i].flui.API = 141.5 / denmixSTD - 131.5;
        } else if (state.input.flashCompleto == 0 || (state.cells[i - 1].fontemassLR) <= 1e-15)
            state.cells[i].flui.API = state.cells[i - 1].flui.API;
        double waterFlowRate1;
        if ((1. - state.cells[i - 1].flui.BSW) > 0)
            waterFlowRate1 = qostd1 * state.cells[i - 1].flui.BSW / (1. - state.cells[i - 1].flui.BSW);
        else
            waterFlowRate1 = state.cells[i - 1].Mliqini * state.cells[i - 1].FW * (1. - state.cells[i - 1].bet) / (1000. * state.cells[i - 1].flui.Denag);
        double waterFlowRate2;
        if ((1. - state.cells[i - 1].acsr.poroso2D.dados.transfer.BSW) > 0)
            waterFlowRate2 = qostd2 * state.cells[i - 1].acsr.poroso2D.dados.transfer.BSW / (1. - state.cells[i - 1].acsr.poroso2D.dados.transfer.BSW);
        else
            waterFlowRate2 = state.cells[i - 1].fontemassLR * fwipr / (1000. * state.cells[i - 1].acsr.poroso2D.dados.flup.Denag);

        if (fabs(waterFlowRate1 + waterFlowRate2 + qostd1 + qostd2) > 1e-15 &&
            (state.cells[i - 1].fontemassLR) > 1e-15)
            state.cells[i].flui.BSW = ((waterFlowRate1 + qostd1) * state.cells[i - 1].flui.BSW + (waterFlowRate2 + qostd2) * state.cells[i - 1].acsr.poroso2D.dados.transfer.BSW) / (waterFlowRate1 + waterFlowRate2 + qostd1 + qostd2);
        else
            state.cells[i].flui.BSW = state.cells[i - 1].flui.BSW;
        if (fabs(waterFlowRate1 + waterFlowRate2) > 1e-15 &&
            (state.cells[i - 1].fontemassLR) > 1e-15)
            state.cells[i].flui.Denag = (waterFlowRate1 * state.cells[i - 1].flui.Denag +
                                    waterFlowRate2 * state.cells[i - 1].acsr.poroso2D.dados.flup.Denag) /
                                   (waterFlowRate1 + waterFlowRate2);
        else
            state.cells[i].flui.Denag = state.cells[i - 1].flui.Denag;
        if (fabs(qostd1 + qostd2) > 1e-15 && state.input.flashCompleto == 0 && (state.cells[i - 1].fontemassLR) > 1e-15) {
            // reavaliacao dos pares necessarios para o modelo de viscosidade de oleo morto
            // ASTM. Isto ÃƒÂ© feito tendo sempre os mesmos valores de temperatura maxima e
            // e minima para a construcao dos pares
            state.cells[i].flui.TempL = temperatureLow;
            state.cells[i].flui.TempH = temperatureHigh;
            state.cells[i].flui.LVisL = (qostd1 * state.cells[i - 1].flui.VisOM(temperatureLow) + qostd2 * state.cells[i - 1].acsr.poroso2D.dados.flup.VisOM(temperatureLow)) / (qostd1 + qostd2);
            state.cells[i].flui.LVisH = (qostd1 * state.cells[i - 1].flui.VisOM(temperatureHigh) + qostd2 * state.cells[i - 1].acsr.poroso2D.dados.flup.VisOM(temperatureHigh)) / (qostd1 + qostd2);
        } else if (state.input.flashCompleto == 0 || fabs(state.cells[i - 1].fontemassLR) <= 1e-15) { // sem vazao de liquido
            state.cells[i].flui.TempL = state.cells[i - 1].flui.TempL;
            state.cells[i].flui.TempH = state.cells[i - 1].flui.TempH;
            state.cells[i].flui.LVisL = state.cells[i - 1].flui.LVisL;
            state.cells[i].flui.LVisH = state.cells[i - 1].flui.LVisH;
        }
    }
    if (state.input.flashCompleto == 0) {
        state.cells[i].flui.RenovaFluido();
        correctGasSpecificGravity(state, i);
    }
    if (state.cells[i - 1].bet > (*state.globals).localtiny * 1e-6 && fabs(state.cells[i - 1].fontemassLR) > 1e-15) {
        // reavaliacao da fracao volumetrica do liquido complementar
        // observar que a IPR pode ter uma fracao de liquido complementar
        // distinta da onservada a esquerda da celula i-1
        double oilVolumeFactorNeighbour = state.cells[i].flui.BOFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        double waterVolumeFactorNeighbour = state.cells[i].flui.BAFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        double waterCutNeighbour = state.cells[i].flui.BSW * waterVolumeFactorNeighbour / (oilVolumeFactorNeighbour + waterVolumeFactorNeighbour * state.cells[i].flui.BSW - state.cells[i].flui.BSW * oilVolumeFactorNeighbour);
        oilVolumeFactorInSitu = state.cells[i].flui.BOFunc(state.cells[i].presaux, tmed);
        waterVolumeFactorInSitu = state.cells[i].flui.BAFunc(state.cells[i].presaux, tmed);
        waterCutInSitu = state.cells[i].flui.BSW * waterVolumeFactorInSitu / (oilVolumeFactorInSitu + waterVolumeFactorInSitu * state.cells[i].flui.BSW - state.cells[i].flui.BSW * oilVolumeFactorInSitu);
        double qlpF = state.cells[i - 1].fontemassLR / state.cells[i - 1].acsr.poroso2D.dados.flup.MasEspLiq(state.cells[i - 1].pres, state.cells[i - 1].temp);
        double qlcF = 0.;
        double betN = (state.cells[i - 1].QL * state.cells[i - 1].bet + qlcF) / (state.cells[i - 1].QL + qlpF + qlcF);
        double oilFlowRate = (state.cells[i - 1].QL + qlpF + qlcF) * (1 - waterCutNeighbour) * (1 - betN) * oilVolumeFactorInSitu / oilVolumeFactor;
        double waterFlowRate;
        if (waterCutInSitu < (1 - (*state.globals).localtiny))
            waterFlowRate = waterCutInSitu * oilFlowRate / (1 - waterCutInSitu);
        else
            waterFlowRate = state.cells[i - 1].QL * (1 - state.cells[i - 1].bet) + qlpF;
        double completionFlowRate = (state.cells[i - 1].QL * (state.cells[i - 1].bet) + qlcF) * state.cells[i - 1].fluicol.MasEspFlu(state.cells[i - 1].pres, state.cells[i - 1].temp) / state.cells[i].fluicol.MasEspFlu(state.cells[i].presaux, tmed);
        state.cells[i].bet = fabs(completionFlowRate) / (fabs(oilFlowRate) + fabs(waterFlowRate) + fabs(completionFlowRate));
    } else
        state.cells[i].bet = state.cells[i - 1].bet;
}

void finalizeSteadyMassWithLiquid(const SteadyStateState &state, int i, double waterCutInSitu, double tmed, double &oilVolumeFactor, double &solutionGasRatio, double &pmed,
                                            double &rhog, double &rhol, double &oilFlowRate, double &masoleo) {
    if (state.input.tipoFluido == 0) {
        solutionGasRatio = state.cells[i].flui.RS(pmed, tmed);
        oilVolumeFactor = state.cells[i].flui.BOFunc(pmed, tmed, solutionGasRatio);
        solutionGasRatio = solutionGasRatio * 6.29 / 35.31467;
        double rhogstd = state.cells[i].flui.Deng * 1.225;
        double rhololeostd = 1000. * 141.5 / (131.5 + state.cells[i].flui.API);
        double rhoa = state.cells[i].flui.Denag * 1000.;

        double denom = state.cells[i].flui.RGO * rhogstd + rhololeostd;
        if (fabs(1 - state.cells[i].flui.BSW) > (*state.globals).localtiny) {
            denom += rhoa * state.cells[i].flui.BSW / (1 - state.cells[i].flui.BSW);
            if (fabs(1 - state.cells[i].bet) > (*state.globals).localtiny * 1e-6 && fabs(state.cells[i].bet) > (*state.globals).localtiny * 1e-10) {
                denom += (state.cells[i].fluicol.MasEspFlu(pmed, tmed) * state.cells[i].bet / (1 - state.cells[i].bet)) * oilVolumeFactor * (1. + waterCutInSitu / (1 - waterCutInSitu));
            }
        }
        if (fabs(1 - state.cells[i].bet) > (*state.globals).localtiny * 1e-6 && fabs(1 - state.cells[i].flui.BSW) > (*state.globals).localtiny)
            oilFlowRate = state.cells[i].MC / denom;

        if (fabs(1 - state.cells[i].bet) > (*state.globals).localtiny * 1e-6 && fabs(1 - state.cells[i].flui.BSW) > (*state.globals).localtiny) {
            if ((*state.globals).tipoFluidoRedeGlob == 0)
                state.cells[i].Mliqini = state.cells[i].MC - (state.cells[i].flui.RGO - solutionGasRatio) * rhogstd * oilFlowRate;
            else {
                masoleo = oilFlowRate * oilVolumeFactor * state.cells[i].flui.MasEspoleo(pmed, tmed, solutionGasRatio) + oilFlowRate * rhoa * state.cells[i].flui.BSW / (1 - state.cells[i].flui.BSW);
                state.cells[i].Mliqini = masoleo + state.cells[i].MComp;
            }
        } else
            state.cells[i].Mliqini = state.cells[i].MC - (state.cells[i - 1].MC - state.cells[i - 1].Mliqini) - state.cells[i - 1].fontemassGR;
    } else {
        double quality = state.cells[i].flui.FracMass(pmed, tmed);
        state.cells[i].Mliqini = state.cells[i].MC * (1. - quality);
    }
    state.cells[i - 1].MliqiniR = state.cells[i].Mliqini;
    if (i < state.lastCell)
        state.cells[i + 1].MliqiniL = state.cells[i].Mliqini;
    pmed = state.cells[i].presaux + state.cells[i - 1].dpB / 98066.5;
    rhog = state.cells[i].rgCi = state.cells[i].flui.MasEspGas(pmed, tmed);

    // calculo das massas especificas na interface a esquerda da celula
    state.cells[i].rpCi = state.cells[i].flui.MasEspLiq(pmed, tmed);
    state.cells[i].rcCi = state.cells[i].fluicol.MasEspFlu(pmed, tmed);
    rhol = (1 - state.cells[i].bet) * state.cells[i].rpCi + state.cells[i].bet * state.cells[i].rcCi;
    // rhol = (1 - celula[i].bet) * celula[i].flui.MasEspLiq(pmed, tmed)
    // vazoes volumetricas:
    state.cells[i].QL = state.cells[i].Mliqini / rhol;
    if (i < state.lastCell)
        state.cells[i + 1].QLL = state.cells[i].QL;
    state.cells[i - 1].QLR = state.cells[i].QL;
    state.cells[i].QG = (state.cells[i].MC - state.cells[i].Mliqini) / rhog;
}

void finalizeSteadyMassNoFlow(const SteadyStateState &state, int i) {
    if (state.input.tipoFluido == 1) {
        state.cells[i].alf = 0.;
    } else {
        state.cells[i].alf = 1.;
    }
    state.cells[i].alfini = state.cells[i].alf;
    state.cells[i - 1].alfR = state.cells[i].alf;
    state.cells[i - 1].alfRini = state.cells[i].alf;
    if (i < state.lastCell)
        state.cells[i + 1].alfL = state.cells[i].alf;
    if (i < state.lastCell)
        state.cells[i + 1].alfLini = state.cells[i].alf;
    state.cells[i].alfPigD = state.cells[i].alf;
    state.cells[i].alfPigDini = state.cells[i].alf;
    state.cells[i].alfPigE = state.cells[i].alf;
    state.cells[i].alfPigEini = state.cells[i].alf;
    if (i < state.lastCell) {
        state.cells[i + 1].betL = state.cells[i].bet;
        state.cells[i + 1].betLini = state.cells[i].bet;
        state.cells[i + 1].betL = state.cells[i].bet;
        state.cells[i + 1].betLini = state.cells[i].bet;
        state.cells[i + 1].betLI = state.cells[i].bet;
    }
    state.cells[i].betini = state.cells[i].bet;
    state.cells[i - 1].betR = state.cells[i].bet;
    state.cells[i - 1].betRini = state.cells[i].bet;
    state.cells[i].betPigD = state.cells[i].bet;
    state.cells[i].betPigDini = state.cells[i].bet;
    state.cells[i].betPigE = state.cells[i].bet;
    state.cells[i].betPigEini = state.cells[i].bet;
    state.cells[i].betI = state.cells[i].bet;
    state.cells[i - 1].betRI = state.cells[i].bet;
}

void finalizeSteadyMassLiquidOnly(const SteadyStateState &state, int i) {
    state.cells[i].alf = 0.;
    state.cells[i].alfini = state.cells[i].alf;
    state.cells[i - 1].alfR = state.cells[i].alf;
    state.cells[i - 1].alfRini = state.cells[i].alf;
    if (i < state.lastCell)
        state.cells[i + 1].alfL = state.cells[i].alf;
    if (i < state.lastCell)
        state.cells[i + 1].alfLini = state.cells[i].alf;
    state.cells[i].alfPigD = state.cells[i].alf;
    state.cells[i].alfPigDini = state.cells[i].alf;
    state.cells[i].alfPigE = state.cells[i].alf;
    state.cells[i].alfPigEini = state.cells[i].alf;
    if (i < state.lastCell) {
        state.cells[i + 1].betL = state.cells[i].bet;
        state.cells[i + 1].betLini = state.cells[i].bet;
        state.cells[i + 1].betL = state.cells[i].bet;
        state.cells[i + 1].betLini = state.cells[i].bet;
        state.cells[i + 1].betLI = state.cells[i].bet;
    }
    state.cells[i].betini = state.cells[i].bet;
    state.cells[i - 1].betR = state.cells[i].bet;
    state.cells[i - 1].betRini = state.cells[i].bet;
    state.cells[i].betPigD = state.cells[i].bet;
    state.cells[i].betPigDini = state.cells[i].bet;
    state.cells[i].betPigE = state.cells[i].bet;
    state.cells[i].betPigEini = state.cells[i].bet;
    state.cells[i].betI = state.cells[i].bet;
    state.cells[i - 1].betRI = state.cells[i].bet;
}

void finalizeSteadyMassGasOnly(const SteadyStateState &state, int i) {
    state.cells[i].alf = 1.;
    state.cells[i].alfini = state.cells[i].alf;
    state.cells[i - 1].alfR = state.cells[i].alf;
    state.cells[i - 1].alfRini = state.cells[i].alf;
    if (i < state.lastCell)
        state.cells[i + 1].alfL = state.cells[i].alf;
    if (i < state.lastCell)
        state.cells[i + 1].alfLini = state.cells[i].alf;
    state.cells[i].alfPigD = state.cells[i].alf;
    state.cells[i].alfPigDini = state.cells[i].alf;
    state.cells[i].alfPigE = state.cells[i].alf;
    state.cells[i].alfPigEini = state.cells[i].alf;
    state.cells[i].c0 = 1.;
    state.cells[i].ud = 0.;
    if (i < state.lastCell) {
        state.cells[i + 1].betL = state.cells[i].bet;
        state.cells[i + 1].betLini = state.cells[i].bet;
        state.cells[i + 1].betL = state.cells[i].bet;
        state.cells[i + 1].betLini = state.cells[i].bet;
        state.cells[i + 1].betLI = state.cells[i].bet;
    }
    state.cells[i].betini = state.cells[i].bet;
    state.cells[i - 1].betR = state.cells[i].bet;
    state.cells[i - 1].betRini = state.cells[i].bet;
    state.cells[i].betPigD = state.cells[i].bet;
    state.cells[i].betPigDini = state.cells[i].bet;
    state.cells[i].betPigE = state.cells[i].bet;
    state.cells[i].betPigEini = state.cells[i].bet;
    state.cells[i].betI = state.cells[i].bet;
    state.cells[i - 1].betRI = state.cells[i].bet;
}

void finalizeSteadyMassTwoPhase(const SteadyStateState &state, int i, double rhog, double rhol) {
    double c0 = 1.;
    double ud = 0.;
    if (fabs(state.cells[i].QL) > (*state.globals).localtiny * 1e-6) {
        if (state.input.tipoModeloDrift == 1) {
            if (state.steadyIteration == 0) { // primeira estimativa, primeira iteracao
                // utiliza-se a fracao de vazio sem escorregamento, pois a propria correlacao para se obter a
                // fracao de vazio depende do valor da fracao de vazio
                if ((fabs(state.cells[i].QG) + fabs(state.cells[i].QL)) > (*state.globals).localtiny) {
                    if (state.convergenceMonitor > 0.01)
                        state.cells[i].alf = fabs(state.cells[i].QG) /
                                        (fabs(state.cells[i].QG) + fabs(state.cells[i].QL));
                    state.updaters.steadyDriftClosure(i, c0, ud);
                } else
                    state.cells[i].alf = 0.;
                state.cells[i].alfini = state.cells[i].alf;
                state.cells[i - 1].alfR = state.cells[i].alf;
                state.cells[i - 1].alfRini = state.cells[i].alf;
                if (i < state.lastCell)
                    state.cells[i + 1].alfL = state.cells[i].alf;
                if (i < state.lastCell)
                    state.cells[i + 1].alfLini = state.cells[i].alf;
                state.cells[i].alfPigD = state.cells[i].alf;
                state.cells[i].alfPigDini = state.cells[i].alf;
                state.cells[i].alfPigE = state.cells[i].alf;
                state.cells[i].alfPigEini = state.cells[i].alf;
            }
            if (fabs(rhog) / rhol > 0.9) {
                c0 = 1.;
                ud = 0.;
            }

            // para o caso permanente, a fracao de vazio e obtida a partir das relacoes de escorregamento
            // portanto, e neste ponto que se obtem Co e Ud:
            else if (fabs(state.cells[i].QG) > (*state.globals).localtiny && fabs(state.cells[i].QL) > (*state.globals).localtiny * 1e-6)
                state.updaters.steadyDriftClosure(i, c0, ud);
            state.cells[i].c0 = c0;
            state.cells[i].ud = ud;
            double area = state.cells[i].duto.area;
            if (fabs(state.cells[i].QG + state.cells[i].QL) > (*state.globals).localtiny) {
                // alfa com escorregamento:
                state.cells[i].alf = state.cells[i].QG / (c0 * (state.cells[i].QG + state.cells[i].QL) + ud * area);
                double alfHomo = state.cells[i].QG / (state.cells[i].QG + state.cells[i].QL);
                if (state.cells[i].alf > 1. - 1e-15 || state.cells[i].alf < 1e-15)
                    state.cells[i].alf = alfHomo;

            } else
                state.cells[i].alf = 0.;
            if (state.cells[i].alf > (1 - (*state.globals).localtiny) && fabs(state.cells[i].QG + state.cells[i].QL) > (*state.globals).localtiny)
                state.cells[i].alf = fabs(state.cells[i].QG) / fabs(state.cells[i].QG + state.cells[i].QL);
            else if (state.cells[i].alf > (1 - (*state.globals).localtiny))
                state.cells[i].alf = 1.;
        } else {
            double holdup;
            double frictionGrad;
            double gravityGrad;
            double totalGrad;
            double reynolds;
            unsigned char flowType;
            char *errorMsg;
            unsigned char errorFlag;
            executarCorrelacao(state.cells, i, 0, state.input.AceleraConvergPerm,
                               state.cells[i - 1].correlacaoMR2,
                               holdup, frictionGrad, gravityGrad, totalGrad,
                               reynolds, flowType);
           if(state.cells[i - 1].correlacaoMR2==16){
            if(flowType=='1' || flowType=='2')state.cells[i].arranjo=0;
            if(flowType=='3')state.cells[i].arranjo=1;
            if(flowType=='4')state.cells[i].arranjo=2;
            if(flowType=='6')state.cells[i].arranjo=-1;
            if(flowType=='5')state.cells[i].arranjo=-2;
           }
           else{
        	   if(flowType=='1')state.cells[i].arranjo=1;
        	   if(flowType=='2')state.cells[i].arranjo=2;
               if(flowType=='3')state.cells[i].arranjo=3;
               if(flowType=='4')state.cells[i].arranjo=4;
               if(flowType=='6')state.cells[i].arranjo=5;
               if(flowType=='5')state.cells[i].arranjo=6;
           }
            state.cells[i].alf = 1. - holdup;
        }
    } else {
        c0 = 1.;
        ud = 0.;
        state.cells[i].alf = 1.;
    }
    if (state.cells[i].alf < 0.)
        state.cells[i].alf = 0.;
    else if (state.cells[i].alf > 1.)
        state.cells[i].alf = 1.;
    // atualizacoes dos valores das fracoes volumetricas da celula i armazendadas em
    // outras celulas, e inclusiove armazendo os valores para "tempo anterior", que nao
    // sao relevantes para o problema permanente mas importantes se o resultado permanente
    // der partida na solucao transiente:
    state.cells[i].alfini = state.cells[i].alf;
    state.cells[i - 1].alfR = state.cells[i].alf;
    state.cells[i - 1].alfRini = state.cells[i].alf;
    if (i < state.lastCell)
        state.cells[i + 1].alfL = state.cells[i].alf;
    if (i < state.lastCell)
        state.cells[i + 1].alfLini = state.cells[i].alf;
    state.cells[i].alfPigD = state.cells[i].alf;
    state.cells[i].alfPigDini = state.cells[i].alf;
    state.cells[i].alfPigE = state.cells[i].alf;
    state.cells[i].alfPigEini = state.cells[i].alf;
    if (i < state.lastCell) {
        state.cells[i + 1].betL = state.cells[i].bet;
        state.cells[i + 1].betLini = state.cells[i].bet;
        state.cells[i + 1].betL = state.cells[i].bet;
        state.cells[i + 1].betLini = state.cells[i].bet;
        state.cells[i + 1].betLI = state.cells[i].bet;
    }
    state.cells[i].betini = state.cells[i].bet;
    state.cells[i - 1].betR = state.cells[i].bet;
    state.cells[i - 1].betRini = state.cells[i].bet;
    state.cells[i].betPigD = state.cells[i].bet;
    state.cells[i].betPigDini = state.cells[i].bet;
    state.cells[i].betPigE = state.cells[i].bet;
    state.cells[i].betPigEini = state.cells[i].bet;
    state.cells[i].betI = state.cells[i].bet;
    state.cells[i - 1].betRI = state.cells[i].bet;
}

}  // namespace

void advanceSteadyMass(const SteadyStateState &state, int i) {
    int mudaRGO = 1;
    if (state.input.flashCompleto == 1)
        mudaRGO = 1;

    if (state.cells[i - 1].acsr.tipo == 1 && state.cells[i - 1].acsr.injg.QGas < 0) {
        state.cells[i - 1].acsr.injg.FluidoPro = state.cells[i - 1].flui;
    } else if (state.cells[i - 1].acsr.tipo == 2 && state.cells[i - 1].acsr.injl.QLiq < 0) {
        state.cells[i - 1].acsr.injl.FluidoPro = state.cells[i - 1].flui;
    } else if (state.cells[i - 1].acsr.tipo == 10 && state.cells[i - 1].acsr.injm.MassC + state.cells[i - 1].acsr.injm.MassG + state.cells[i - 1].acsr.injm.MassP < 0) {
        state.cells[i - 1].acsr.injm.FluidoPro = state.cells[i - 1].flui;
    }

    state.updaters.updateSource(i - 1); // metodo que verifica se existe uma fonte na celula e calcula o valor das
    // vazoes massicas de liquido produzido (oleo+agua), gas e liquido complementar
    // relacao entre fonte a esquerda e a direita de uma celula:
    state.cells[i].fontemassLL = state.cells[i - 1].fontemassLR;
    state.cells[i].fontemassCL = state.cells[i - 1].fontemassCR;
    state.cells[i].fontemassGL = state.cells[i - 1].fontemassGR;
    state.cells[i].MComp = state.cells[i - 1].MComp + state.cells[i - 1].fontemassCR;

    double residenceTimeSource = 0.;
    double sinalQ = 1.;
    if (i > 1) {
        if (fabs(state.cells[i - 1].QL + state.cells[i - 1].QG) > 1e-15)
            sinalQ = (state.cells[i - 1].QL + state.cells[i - 1].QG) / fabs(state.cells[i - 1].QL + state.cells[i - 1].QG);
    } else {
        if (fabs(state.cells[i - 1].fontemassLR + state.cells[i - 1].fontemassCR + state.cells[i - 1].fontemassGR) > 1e-15)
            sinalQ = (state.cells[i - 1].fontemassLR + state.cells[i - 1].fontemassCR + state.cells[i - 1].fontemassGR) /
                     fabs(state.cells[i - 1].fontemassLR + state.cells[i - 1].fontemassCR + state.cells[i - 1].fontemassGR);
    }
    double oilVolumeFactorInSitu;
    double waterVolumeFactorInSitu;
    double waterCutInSitu = 0;
    // definicao de temperaturas maximas e minimas para uma eventual reavaliacao do metodo
    // ASTM quando ocorre mistura de fluidos e se deseja atualizar o modelo
    // de viscosidade de oleo morto se este for o ASTM que trabalha com um par de temperatura
    double temperatureLow;
    if (state.input.flashCompleto == 1 && (state.input.tabent.tmin - 0) > (*state.globals).localtiny)
        temperatureLow = state.input.tabent.tmin + 0.1;
    else
        temperatureLow = 0.;
    double temperatureHigh;
    if (state.input.flashCompleto == 1 && (70 - state.input.tabent.tmax) < (*state.globals).localtiny)
        temperatureHigh = state.input.tabent.tmax - 0.1;
    else
        temperatureHigh = 70.;
    double oilVolumeFactor;
    double waterVolumeFactor;
    double solutionGasRatio;
    // calcula os valores de RS, Bo e Ba na celula anterior a i-esima celula
    if (state.cells[i - 1].flui.RGO < 1e7) {
        solutionGasRatio = state.cells[i - 1].flui.RS(state.cells[i - 1].pres, state.cells[i - 1].temp);
        oilVolumeFactor = state.cells[i - 1].flui.BOFunc(state.cells[i - 1].pres, state.cells[i - 1].temp, solutionGasRatio);
        waterVolumeFactor = state.cells[i - 1].flui.BAFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        solutionGasRatio = solutionGasRatio * 6.29 / 35.31467;
    } else {
        oilVolumeFactor = 1;
        solutionGasRatio = 0;
        waterVolumeFactor = 0.;
    }
    // BSW in-situ da celula anterior, na marcha, a i-esima celula
    state.cells[i - 1].FW = state.cells[i - 1].flui.BSW * waterVolumeFactor / (oilVolumeFactor + waterVolumeFactor * state.cells[i - 1].flui.BSW - state.cells[i - 1].flui.BSW * oilVolumeFactor);
    state.cells[i - 1].FWini = state.cells[i - 1].FW;
    double tmed;
    // temperatura na fronteira entre a celula i-1 e a celula i, da segunda iteracao em diante, pode-se
    // usar a temperatura da celula i, pois ja a tem calculada, mas isto pode ser um complicador de
    // convergencia, mais seguro manter o criterio em todas as iteracoes, neste caso, a temperatura na
    // fronteira e admitida = a temperatura da celula i-1
    if (state.steadyIteration != 0 && state.input.AceleraConvergPerm == 0)
        tmed = (state.cells[i].dx * state.cells[i].temp + state.cells[i].dxL * state.cells[i].tempL) / (state.cells[i].dx + state.cells[i].dxL);
    else
        tmed = state.cells[i - 1].temp;
    // primeiro teste: nÃƒÂ£o ha fontes na celula i-1:
    if (state.cells[i - 1].acsr.tipo != 1 && state.cells[i - 1].acsr.tipo != 2 && state.cells[i - 1].acsr.tipo != 3 && state.cells[i - 1].acsr.tipo != 10 && (state.cells[i - 1].acsr.tipo != 9 || (state.cells[i - 1].acsr.tipo == 9 && state.cells[i - 1].acsr.fontechk.abertura <= 1e-6)) && state.cells[i - 1].acsr.tipo != 15 && state.cells[i - 1].acsr.tipo != 16) {
        applySteadyMassWithoutSource(state, i, mudaRGO, oilVolumeFactor, solutionGasRatio, tmed, oilVolumeFactorInSitu, waterVolumeFactorInSitu, waterCutInSitu);
    }
    // caso em que se tem uma fonte de gas na celula i-1, o que mudara a RGO e a densidade de gas em i
    else if (state.cells[i - 1].acsr.tipo == 1 && state.cells[i - 1].acsr.injg.seco == 1) {
        applySteadyMassDryGasInjection(state, i, mudaRGO, temperatureLow, temperatureHigh, oilVolumeFactor, waterVolumeFactor, solutionGasRatio, tmed, residenceTimeSource, oilVolumeFactorInSitu, waterVolumeFactorInSitu, waterCutInSitu);
    } else if (state.cells[i - 1].acsr.tipo == 1 && state.cells[i - 1].acsr.injg.seco == 0) {
        applySteadyMassWetGasInjection(state, i, mudaRGO, temperatureLow, temperatureHigh, oilVolumeFactor, waterVolumeFactor, solutionGasRatio, tmed, residenceTimeSource, oilVolumeFactorInSitu, waterVolumeFactorInSitu, waterCutInSitu);
    }
    // caso de fonte de liquido na celula i-1:
    else if (state.cells[i - 1].acsr.tipo == 2) {
        applySteadyMassLiquidInjection(state, i, mudaRGO, temperatureLow, temperatureHigh, oilVolumeFactor, waterVolumeFactor, solutionGasRatio, tmed, residenceTimeSource, oilVolumeFactorInSitu, waterVolumeFactorInSitu, waterCutInSitu);
    }
    // caso de IPR na celula i-1:
    else if (state.cells[i - 1].acsr.tipo == 3) {
        applySteadyMassInflowPerformance(state, i, mudaRGO, temperatureLow, temperatureHigh, oilVolumeFactor, waterVolumeFactor, solutionGasRatio, tmed, residenceTimeSource, oilVolumeFactorInSitu, waterVolumeFactorInSitu, waterCutInSitu);
    }
    // caso de fonte de massa na celula i-1:
    else if (state.cells[i - 1].acsr.tipo == 10) {
        applySteadyMassMultipleSource(state, i, mudaRGO, temperatureLow, temperatureHigh, oilVolumeFactor, waterVolumeFactor, solutionGasRatio, tmed, residenceTimeSource, oilVolumeFactorInSitu, waterVolumeFactorInSitu, waterCutInSitu);
    }
    // caso especial, fonte de vazamento:
    else if (state.cells[i - 1].acsr.tipo == 9 && state.cells[i - 1].acsr.fontechk.abertura > 1e-6) {
        applySteadyMassLeakSource(state, i, mudaRGO, temperatureLow, temperatureHigh, oilVolumeFactor, waterVolumeFactor, solutionGasRatio, tmed, residenceTimeSource, oilVolumeFactorInSitu, waterVolumeFactorInSitu, waterCutInSitu);
    } else if (state.cells[i - 1].acsr.tipo == 15) {
        applySteadyMassRadialPorous(state, i, mudaRGO, temperatureLow, temperatureHigh, oilVolumeFactor, waterVolumeFactor, solutionGasRatio, tmed, residenceTimeSource, oilVolumeFactorInSitu, waterVolumeFactorInSitu, waterCutInSitu);
    } else if (state.cells[i - 1].acsr.tipo == 16) {
        applySteadyMassPorous2D(state, i, mudaRGO, temperatureLow, temperatureHigh, oilVolumeFactor, waterVolumeFactor, solutionGasRatio, tmed, residenceTimeSource, oilVolumeFactorInSitu, waterVolumeFactorInSitu, waterCutInSitu);
    }

    // atualizaÃ§Ã£o da vazao massica da mistura:
    state.cells[i].MC = state.cells[i - 1].MC + state.cells[i - 1].fontemassCR + state.cells[i - 1].fontemassLR + state.cells[i - 1].fontemassGR;
    state.cells[i - 1].MR = state.cells[i].MC;
    if (i < state.lastCell)
        state.cells[i + 1].ML = state.cells[i].MC;
    state.cells[i - 1].MRini = state.cells[i - 1].MR;
    double pmed = state.cells[i].presaux;
    double rhog;
    double rhol = 1000.;
    // calculo da vazao massica de liquido e das vazoes volumetricas:
    double oilFlowRate = 0.;
    double masoleo;
    if (state.cells[i].flui.RGO < 1e7) { // Caso exista liquido
        finalizeSteadyMassWithLiquid(state, i, waterCutInSitu, tmed, oilVolumeFactor, solutionGasRatio, pmed, rhog, rhol, oilFlowRate, masoleo);
    } else {
        state.cells[i].Mliqini = 0.;
        state.cells[i - 1].MliqiniR = state.cells[i].Mliqini;
        if (i < state.lastCell)
            state.cells[i + 1].MliqiniL = state.cells[i].Mliqini;
        state.cells[i].QL = 0.;
        if (i < state.lastCell)
            state.cells[i + 1].QLL = state.cells[i].QL;
        state.cells[i - 1].QLR = state.cells[i].QL;
        pmed = state.cells[i].presaux + state.cells[i - 1].dpB / 98066.5;
        state.cells[i].rpCi = state.cells[i].flui.MasEspLiq(pmed, tmed);
        state.cells[i].rcCi = state.cells[i].fluicol.MasEspFlu(pmed, tmed);
        rhog = state.cells[i].rgCi = state.cells[i].flui.MasEspGas(pmed, tmed);
        state.cells[i].QG = (state.cells[i].MC) / rhog;
    }
    // Definicao das fracoes volumetricas:
    if (fabs(state.cells[i].QG + state.cells[i].QL) < (*state.globals).localtiny) {
        finalizeSteadyMassNoFlow(state, i);
    } else if (fabs(state.cells[i].QG) < (*state.globals).localtiny) { // caso so exista liquido:
        finalizeSteadyMassLiquidOnly(state, i);
    } else if (fabs(state.cells[i].QL) < (*state.globals).localtiny * 1e-6) { // caso so exista gas:
        finalizeSteadyMassGasOnly(state, i);
    } else { // caso bifasico
        finalizeSteadyMassTwoPhase(state, i, rhog, rhol);
    }

    if (fabs(state.cells[i].QL) > 1e-15 && fabs(state.cells[i].MComp) > 1e-15) {
        double dxmed = 0.5 * (state.cells[i - 1].dx + state.cells[i].dx);
        double hol0 = 1. - state.cells[i - 1].alf;
        double hol1 = 1. - state.cells[i].alf;
        state.cells[i].fluicol.TR = (state.cells[i - 1].fluicol.TR * state.cells[i - 1].MComp + state.cells[i - 1].fontemassCR * residenceTimeSource) / state.cells[i].MComp;
        state.cells[i].fluicol.TR = state.cells[i].fluicol.TR + (0.5 * state.cells[i - 1].duto.area * state.cells[i - 1].dx * hol0 +
                                                       0.5 * state.cells[i].duto.area * state.cells[i].dx * hol1) /
                                                          state.cells[i].QL;
    } else
        state.cells[i].fluicol.TR = state.cells[i - 1].fluicol.TR;
}

double areaChangePressureDrop(const SteadyStateState &state, int i, double rhomix, double reynolds, double jmix) {
    double dpArea = 0.;
    if ((state.cells[i].duto.area != state.cells[i].dutoR.area && reynolds > 2400) && (state.cells[i].mudaArea == 1 && (fabs(jmix) >= 0.1))) {
        double areaMenor;
        double areaMaior;
        double bernou;
        double perdLoc;
        bernou = 0.5 * (state.cells[i].MR * state.cells[i].MR / (state.cells[i].duto.area * state.cells[i].duto.area * rhomix)) *
                 (1. - pow(state.cells[i].duto.area / state.cells[i].dutoR.area, 2.));
        if (state.cells[i].duto.area > state.cells[i].dutoR.area) {
            areaMenor = state.cells[i].dutoR.area;
            areaMaior = state.cells[i].duto.area;
            if (state.cells[i].MR > 0.) {
                perdLoc = 0.42 * 0.5 * (state.cells[i].MR * state.cells[i].MR / (state.cells[i].duto.area * state.cells[i].duto.area * rhomix)) *
                          (1. - pow(state.cells[i].duto.area / state.cells[i].dutoR.area, 2.));
            } else {
                double expanse = (1. - pow(areaMenor / areaMaior, 2.));
                perdLoc = 0.5 * (state.cells[i].MR * state.cells[i].MR / (areaMenor * areaMenor * rhomix)) * expanse * expanse;
            }
        } else {
            areaMenor = state.cells[i].duto.area;
            areaMaior = state.cells[i].dutoR.area;
            if (state.cells[i].MR < 0.) {
                perdLoc = 0.42 * 0.5 * (state.cells[i].MR * state.cells[i].MR / (state.cells[i].duto.area * state.cells[i].duto.area * rhomix)) *
                          (1. - pow(state.cells[i].duto.area / state.cells[i].dutoR.area, 2.));
            } else {
                double expanse = (1. - pow(areaMenor / areaMaior, 2.));
                perdLoc = -0.5 * (state.cells[i].MR * state.cells[i].MR / (areaMenor * areaMenor * rhomix)) * expanse * expanse;
            }
        }
        dpArea = bernou + perdLoc;
    }
    return dpArea / 98066.5;
}

double steadyPressureAtLastCell(const SteadyStateState &state) {

    double dx = 0.5 * state.cells[state.lastCell].dx;
    double diameter = state.cells[state.lastCell].duto.a;
    double area = 0.25 * M_PI * diameter * diameter;
    double perimeter = state.cells[state.lastCell].duto.peri;
    double alfmed = state.cells[state.lastCell].alf;
    double betmed = state.cells[state.lastCell].bet;
    double rhog = state.cells[state.lastCell].flui.MasEspGas(state.cells[state.lastCell].pres, state.cells[state.lastCell].temp);
    double rhol = (1. - betmed) * state.cells[state.lastCell].flui.MasEspLiq(state.cells[state.lastCell].pres, state.cells[state.lastCell].temp) + betmed * state.cells[state.lastCell].fluicol.MasEspFlu(state.cells[state.lastCell].pres, state.cells[state.lastCell].temp);
    double ugsmed = (state.cells[state.lastCell].MC - state.cells[state.lastCell].Mliqini) / (area * rhog);
    double ulsmed = state.cells[state.lastCell].Mliqini / (area * rhol);

    double j = ugsmed + ulsmed;
    double sinalJ = 1.;
    if (fabs(j) > 1e-15)
        sinalJ = j / fabs(j);

    double rhomix = alfmed * rhog + (1 - alfmed) * rhol;
    double viscmix = alfmed * state.cells[state.lastCell].flui.ViscGas(state.cells[state.lastCell].pres, state.cells[state.lastCell].temp) + (1 - alfmed) * ((1. - betmed) * state.cells[state.lastCell].flui.ViscOleo(state.cells[state.lastCell].pres, state.cells[state.lastCell].temp) + betmed * state.cells[state.lastCell].fluicol.VisFlu(state.cells[state.lastCell].pres, state.cells[state.lastCell].temp));

    double mixtureReynolds;
    if (state.cells[state.lastCell].duto.revest == 0)
        mixtureReynolds = state.cells[state.lastCell].Rey(state.cells[state.lastCell].duto.a, j, rhomix, viscmix);
    else {
        double dhid = 4 * area / perimeter;
        mixtureReynolds = state.cells[state.lastCell].Rey(dhid, j, rhomix, viscmix);
    }
    double frictionFactor = state.cells[state.lastCell].fric(mixtureReynolds, state.cells[state.lastCell].duto.rug / diameter);
    double gradfric = state.cells[state.lastCell].dPdLFric * (0.5 * frictionFactor * rhomix * (fabs(j) * j) * perimeter * dx / area);
    double gradhidro = state.cells[state.lastCell].dPdLHidro * (9.82 * sin(state.cells[state.lastCell].duto.teta) * rhomix * dx);
    return -(1. * gradfric + gradhidro) / 98066.5;
}

void advanceUpstreamSteadyPressure(const SteadyStateState &state, int i, int rungeKuttaStage) {

    if (state.input.tipoModeloDrift == 1) {
        double dx = 0.5 * state.cells[i].dxL;
        double diameter = state.cells[i].dutoL.a;
        double area = 0.25 * M_PI * diameter * diameter;
        double perimeter = state.cells[i].dutoL.peri;
        if (rungeKuttaStage == 0) {
            double alfmed = state.cells[i - 1].alf;
            double betmed = state.cells[i - 1].bet;
            double rhog = state.cells[i - 1].flui.MasEspGas(state.cells[i - 1].pres, state.cells[i - 1].temp);
            double rhol = (1. - betmed) * state.cells[i - 1].flui.MasEspLiq(state.cells[i - 1].pres, state.cells[i - 1].temp) + betmed * state.cells[i - 1].fluicol.MasEspFlu(state.cells[i - 1].pres, state.cells[i - 1].temp);
            double ugsmed = (state.cells[i - 1].MC - state.cells[i - 1].Mliqini) / (area * rhog);
            double ulsmed = state.cells[i - 1].Mliqini / (area * rhol);
            double j = ugsmed + ulsmed;
            double sinalJ = 1.;
            if (fabs(j) > 1e-15)
                sinalJ = j / fabs(j);

            double rhomix = alfmed * rhog + (1 - alfmed) * rhol;
            double visl = ((1. - betmed) * state.cells[i - 1].flui.ViscOleo(state.cells[i - 1].pres, state.cells[i - 1].temp) + betmed * state.cells[i - 1].fluicol.VisFlu(state.cells[i - 1].pres, state.cells[i - 1].temp));
            double viscmix = alfmed * state.cells[i - 1].flui.ViscGas(state.cells[i - 1].pres, state.cells[i - 1].temp) + (1 - alfmed) * visl;

            double mixtureReynolds;
            if (state.cells[i].dutoL.revest == 0)
                mixtureReynolds = state.cells[i - 1].Rey(state.cells[i].dutoL.a, j, rhomix, viscmix);
            else {
                double dhid = 4 * area / perimeter;
                mixtureReynolds = state.cells[i - 1].Rey(dhid, j, rhomix, viscmix);
            }
            double frictionFactor = state.cells[i - 1].fric(mixtureReynolds, state.cells[i].dutoL.rug / diameter);
            if (i > 1 && state.cells[i - 1].fluicol.tipoF == 2) {
                frictionFactor *= (1 - state.cells[i - 1].dR);
            }
            double gradfric = state.cells[i - 1].dPdLFric * (0.5 * frictionFactor * rhomix * (fabs(j) * j) * perimeter * dx / area);
            double gradhidro = state.cells[i - 1].dPdLHidro * (9.82 * sin(state.cells[i].dutoL.teta) * rhomix * dx);
            state.cells[i].presaux = state.cells[i - 1].pres - (gradfric + gradhidro) / 98066.5;

            state.cells[i].termoHidro = sinalJ * gradhidro / dx;
            state.cells[i].termoFric = gradfric / dx;
        } else {
            double alfmed = state.cells[i].alf;
            double betmed = state.cells[i].bet;
            double razdx = state.cells[i].dx / (state.cells[i].dx + state.cells[i].dxL);
            double tmed = razdx * state.cells[i].temp + (1. - razdx) * state.cells[i - 1].temp;
            double pmed = state.cells[i].presaux;
            double rhog = state.cells[i - 1].flui.MasEspGas(pmed, tmed);
            double rhol = (1. - betmed) * state.cells[i - 1].flui.MasEspLiq(pmed, tmed) + betmed * state.cells[i - 1].fluicol.MasEspFlu(pmed, tmed);
            double ugsmed = (state.cells[i].QG) / (area);
            double ulsmed = state.cells[i].QL / (area);
            double j = ugsmed + ulsmed;
            double sinalJ = 1.;
            if (fabs(j) > 1e-15)
                sinalJ = j / fabs(j);

            double rhomix = alfmed * rhog + (1 - alfmed) * rhol;
            double viscmix = alfmed * state.cells[i - 1].flui.ViscGas(pmed, tmed) + (1 - alfmed) * ((1. - betmed) * state.cells[i - 1].flui.ViscOleo(pmed, tmed) + betmed * state.cells[i - 1].fluicol.VisFlu(pmed, tmed));

            double mixtureReynolds;
            if (state.cells[i].dutoL.revest == 0)
                mixtureReynolds = state.cells[i - 1].Rey(state.cells[i].dutoL.a, j, rhomix, viscmix);
            else {
                double dhid = 4 * area / perimeter;
                mixtureReynolds = state.cells[i - 1].Rey(dhid, j, rhomix, viscmix);
            }
            double frictionFactor = state.cells[i - 1].fric(mixtureReynolds, state.cells[i].dutoL.rug / diameter);
            double gradfric = state.cells[i - 1].dPdLFric * (0.5 * frictionFactor * rhomix * (fabs(j) * j) * perimeter * dx / area);
            double gradhidro = state.cells[i - 1].dPdLHidro * (9.82 * sin(state.cells[i].dutoL.teta) * rhomix * dx);
            state.cells[i].presaux = state.cells[i - 1].pres - (gradfric + gradhidro) / 98066.5;

            state.cells[i].termoHidro = sinalJ * gradhidro / dx;
            state.cells[i].termoFric = gradfric / dx;
        }
    } else {
        double holdup;
        double frictionGrad;
        double gravityGrad;
        double totalGrad;
        double reynolds;
        unsigned char flowType;
        char *errorMsg;
        unsigned char errorFlag;

        double dx = 0.5 * state.cells[i].dxL;
        double diameter = state.cells[i].dutoL.a;
        double area = 0.25 * M_PI * diameter * diameter;
        double perimeter = state.cells[i].dutoL.peri;
        executarCorrelacao(state.cells, i, 1, state.input.AceleraConvergPerm,
                           state.cells[i - 1].correlacaoMR2,
                           holdup, frictionGrad, gravityGrad, totalGrad,
                           reynolds, flowType);

        double gradfric = state.cells[i - 1].dPdLFric * frictionGrad * 22620.6 * dx;
        double gradhidro = state.cells[i - 1].dPdLHidro * gravityGrad * 22620.6 * dx;
        state.cells[i].presaux = state.cells[i - 1].pres - (gradfric + gradhidro) / 98066.5;

        state.cells[i].termoHidro = gradhidro / dx;
        state.cells[i].termoFric = gradfric / dx;
    }
}

void advanceDownstreamSteadyPressure(const SteadyStateState &state, int i, int rungeKuttaStage) {

    double dx;
    double diameter;
    double area;
    double perimeter;
    double alfmed;
    double betmed;
    double rhog;
    double rhol;
    double ugsmed;
    double ulsmed;
    double j;

    double rhomix;
    double viscmix;
    double mixtureReynolds;
    double frictionFactor;
    double gradfric;
    double gradhidro;

    double tmed;

    double dpArea = 0.;

    dx = 0.5 * state.cells[i].dx;
    diameter = state.cells[i].duto.a;
    area = 0.25 * M_PI * diameter * diameter;
    perimeter = state.cells[i].duto.peri;

    if (state.input.tipoModeloDrift == 1) {
        if (rungeKuttaStage == 0) {
            alfmed = state.cells[i].alf;
            betmed = state.cells[i].bet;
            double razdx = state.cells[i].dx / (state.cells[i].dx + state.cells[i].dxL);
            if (state.input.AceleraConvergPerm == 0)
                tmed = razdx * state.cells[i].temp + (1. - razdx) * state.cells[i - 1].temp;
            else
                tmed = state.cells[i - 1].temp;
            double pmed = state.cells[i].presaux + state.cells[i - 1].dpB / 98066.5;
            rhog = state.cells[i].rgCi; // celula[i].flui.MasEspGas(pmed, tmed);
            rhol = (1. - betmed) * state.cells[i].rpCi + betmed * state.cells[i].rcCi;
            ugsmed = (state.cells[i].QG) / (area);
            ulsmed = state.cells[i].QL / (area);
            j = ugsmed + ulsmed;
            double sinalJ = 1.;
            if (fabs(j) > 1e-15)
                sinalJ = j / fabs(j);

            rhomix = alfmed * rhog + (1 - alfmed) * rhol;
            double visl = ((1. - betmed) * state.cells[i].flui.ViscOleo(pmed, tmed) + betmed * state.cells[i].fluicol.VisFlu(pmed, tmed));
            viscmix = alfmed * state.cells[i].flui.ViscGas(pmed, tmed) + (1 - alfmed) * visl;

            if (state.cells[i].duto.revest == 0)
                mixtureReynolds = state.cells[i].Rey(state.cells[i].duto.a, j, rhomix, viscmix);
            else {
                double dhid = 4 * area / perimeter;
                mixtureReynolds = state.cells[i].Rey(dhid, j, rhomix, viscmix);
            }
            frictionFactor = state.cells[i].fric(mixtureReynolds, state.cells[i].duto.rug / diameter);
            if (state.cells[i].fluicol.tipoF == 2) {
                double ulmed = 0.;
                double liquidReynolds;
                if (fabs(state.cells[i].MComp) > 1e-15) {
                    if (alfmed < 1 - 1e-15)
                        ulmed = ulsmed / (1 - alfmed);
                    else
                        ulmed = 0.;
                    if (state.cells[i].duto.revest == 0)
                        liquidReynolds = state.cells[i].Rey(state.cells[i].duto.a, ulmed, rhol, visl);
                    else {
                        double dhid = 4 * area / perimeter;
                        liquidReynolds = state.cells[i].Rey(dhid, ulmed, rhol, visl);
                    }
                    state.cells[i].dR = state.cells[i].fluicol.calcDR(liquidReynolds);
                } else
                    state.cells[i].dR = 0.;
                frictionFactor *= (1 - state.cells[i].dR);
            }
            gradfric = state.cells[i].dPdLFric * (0.5 * frictionFactor * rhomix * (fabs(j) * j) * perimeter * dx / area);
            gradhidro = state.cells[i].dPdLHidro * (9.82 * sin(state.cells[i].duto.teta) * rhomix * dx);
            if (state.cells[i].mudaArea == 1)
                dpArea = areaChangePressureDrop(state, i - 1, rhomix, mixtureReynolds, fabs(j));

            state.cells[i].pres = pmed - (gradfric + gradhidro) / 98066.5 + dpArea;
        } else {
            alfmed = state.cells[i].alf;
            betmed = state.cells[i].bet;
            double razdx = state.cells[i].dx / (state.cells[i].dx + state.cells[i].dxL);
            tmed = state.cells[i].temp;
            double pmed = state.cells[i].pres;
            rhog = state.cells[i].flui.MasEspGas(pmed, tmed);
            rhol = (1. - betmed) * state.cells[i].flui.MasEspLiq(pmed, tmed) + betmed * state.cells[i].fluicol.MasEspFlu(pmed, tmed);
            ugsmed = (state.cells[i].MC - state.cells[i].Mliqini) / (area * rhog);
            ulsmed = state.cells[i].Mliqini / (area * rhol);
            j = ugsmed + ulsmed;
            double sinalJ = 1.;
            if (fabs(j) > 1e-15)
                sinalJ = j / fabs(j);

            rhomix = alfmed * rhog + (1 - alfmed) * rhol;
            viscmix = alfmed * state.cells[i].flui.ViscGas(pmed, tmed) + (1 - alfmed) * ((1. - betmed) * state.cells[i].flui.ViscOleo(pmed, tmed) + betmed * state.cells[i].fluicol.VisFlu(pmed, tmed));

            if (state.cells[i].duto.revest == 0)
                mixtureReynolds = state.cells[i].Rey(state.cells[i].duto.a, j, rhomix, viscmix);
            else {
                double dhid = 4 * area / perimeter;
                mixtureReynolds = state.cells[i].Rey(dhid, j, rhomix, viscmix);
            }
            frictionFactor = state.cells[i].fric(mixtureReynolds, state.cells[i].duto.rug / diameter);
            gradfric = state.cells[i].dPdLFric * (0.5 * frictionFactor * rhomix * (fabs(j) * j) * perimeter * dx / area);
            gradhidro = state.cells[i].dPdLHidro * (9.82 * sin(state.cells[i].duto.teta) * rhomix * dx);
            if (state.cells[i].mudaArea == 1)
                dpArea = areaChangePressureDrop(state, i - 1, rhomix, mixtureReynolds, fabs(j));

            state.cells[i].pres = state.cells[i].presaux + state.cells[i - 1].dpB / 98066.5 - (gradfric + gradhidro) / 98066.5 + dpArea;
        }
    } else {
        double holdup;
        double frictionGrad;
        double gravityGrad;
        double totalGrad;
        double reynolds;
        unsigned char flowType;
        char *errorMsg;
        unsigned char errorFlag;
        executarCorrelacao(state.cells, i, 2, state.input.AceleraConvergPerm,
                           state.cells[i - 1].correlacaoMR2,
                           holdup, frictionGrad, gravityGrad, totalGrad,
                           reynolds, flowType);

        double gradfric = state.cells[i - 1].dPdLFric * frictionGrad * 22620.6 * dx;
        double gradhidro = state.cells[i - 1].dPdLHidro * gravityGrad * 22620.6 * dx;
        state.cells[i].pres = state.cells[i].presaux + state.cells[i - 1].dpB / 98066.5 - (gradfric + gradhidro) / 98066.5;

        state.cells[i].termoHidro = gradhidro / dx;
        state.cells[i].termoFric = gradfric / dx;
    }
}

void advanceSteadyMassTransfer(const SteadyStateState &state, int i) {

    double waterCutLocal;
    double waterCutLeftCell;

    double razdx = state.cells[i].dxR / (state.cells[i].dx + state.cells[i].dxR);
    double razdxL = state.cells[i].dx / (state.cells[i].dx + state.cells[i].dxL);
    double tmed;
    tmed = razdx * state.cells[i + 1].temp + (1 - razdx) * state.cells[i].temp;
    double tmed0;
    if (i > 1) {
        tmed0 = razdxL * state.cells[i].temp + (1 - razdxL) * state.cells[i - 1].temp;
    } else
        tmed0 = state.cells[i].temp;

    double oilVolumeFactorLocal = state.cells[i].flui.BOFunc(state.cells[i].pres, state.cells[i].temp);
    double waterVolumeFactorLocal = state.cells[i].flui.BAFunc(state.cells[i].pres, state.cells[i].temp);
    waterCutLocal = state.cells[i].flui.BSW * waterVolumeFactorLocal / (oilVolumeFactorLocal + waterVolumeFactorLocal * state.cells[i].flui.BSW - state.cells[i].flui.BSW * oilVolumeFactorLocal);
    double solutionGasRatioRightFace = state.cells[i].flui.RS(state.cells[i + 1].presaux, tmed);
    double oilVolumeFactorRightFace = state.cells[i].flui.BOFunc(state.cells[i + 1].presaux, tmed);
    double oilVolumeFactorLeftCell;
    double waterVolumeFactorLeftCell;
    double solutionGasRatioLeftFace;
    double oilVolumeFactorLeftFace;
    double waterVolumeFactorLeftFace;
    double dengD = state.cells[i].flui.Deng;
    double dengE;
    double dissolvedGasGravityRatio = state.cells[i].flui.rDgD;
    double freeGasGravityRatio;
    if (i > 0) {
        oilVolumeFactorLeftCell = state.cells[i - 1].flui.BOFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        waterVolumeFactorLeftCell = state.cells[i - 1].flui.BAFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        waterCutLeftCell = state.cells[i - 1].flui.BSW * waterVolumeFactorLeftCell / (oilVolumeFactorLeftCell + waterVolumeFactorLeftCell * state.cells[i - 1].flui.BSW - state.cells[i - 1].flui.BSW * oilVolumeFactorLeftCell);
        solutionGasRatioLeftFace = state.cells[i - 1].flui.RS(state.cells[i].presaux, tmed0);
        oilVolumeFactorLeftFace = state.cells[i - 1].flui.BOFunc(state.cells[i].presaux, tmed0);
        waterVolumeFactorLeftFace = state.cells[i - 1].flui.BAFunc(state.cells[i].presaux, tmed0);
        dengE = state.cells[i - 1].flui.Deng;
        freeGasGravityRatio = state.cells[i - 1].flui.rDgD;
    } else {
        oilVolumeFactorLeftCell = state.cells[i].flui.BOFunc(state.cells[i].pres, state.cells[i].temp);
        waterVolumeFactorLeftCell = state.cells[i].flui.BAFunc(state.cells[i].pres, state.cells[i].temp);

        waterCutLeftCell = state.cells[i].flui.BSW * waterVolumeFactorLeftCell / (oilVolumeFactorLeftCell + waterVolumeFactorLeftCell * state.cells[i].flui.BSW - state.cells[i].flui.BSW * oilVolumeFactorLeftCell);
        solutionGasRatioLeftFace = state.cells[i].flui.RS(state.cells[i].presaux, tmed0);
        oilVolumeFactorLeftFace = state.cells[i].flui.BOFunc(state.cells[i].presaux, tmed0);
        waterVolumeFactorLeftFace = state.cells[i].flui.BAFunc(state.cells[i].presaux, tmed0);
        dengE = state.cells[i].flui.Deng;
        freeGasGravityRatio = state.cells[i].flui.rDgD;
    }
    double betI;
    double betL;
    if (i < state.lastCell) {
        betI = state.cells[i].betPigD;
        if (state.cells[i + 1].Mliqini < 0)
            betI = state.cells[i + 1].betPigE;
    } else
        betI = state.cells[i].betPigD;
    if (i > 0) {
        betL = state.cells[i - 1].betPigD;
        if (state.cells[i].Mliqini < 0)
            betL = state.cells[i].betPigE;
    } else
        betL = state.cells[i].betPigE;

    if (state.cells[i].acsr.tipo == 0)
        state.cells[i].transmassR = (-(state.cells[i + 1].QL * (1. - betI) * dissolvedGasGravityRatio * dengD * 1.225 * (1. - waterCutLocal) * solutionGasRatioRightFace * (6.29 / 35.31467) / oilVolumeFactorRightFace) + (state.cells[i].QL * (1. - betL) * freeGasGravityRatio * dengE * 1.225 * (1. - waterCutLeftCell) * solutionGasRatioLeftFace * (6.29 / 35.31467) / oilVolumeFactorLeftFace));
    else
        state.cells[i].transmassR = 0.;

    state.cells[i].transmassR /= state.cells[i].dx;
    state.cells[i + 1].transmassL = state.cells[i].transmassR;
    state.cells[i].FonteMudaFase = state.cells[i].transmassR;
    if (i == state.lastCell - 1) {
        state.cells[i + 1].transmassR = state.cells[i].transmassR;
        state.cells[i + 1].FonteMudaFase = state.cells[i].transmassR;
    }
}

void advanceSteadyGasMassTransfer(const SteadyStateState &state, int i) {

    if (state.cells[i].acsr.tipo == 0 && fabs(state.cells[i + 1].flui.dVaporMassFraction - state.cells[i].flui.dVaporMassFraction) < 0.2) {
        state.cells[i].transmassR = ((state.cells[i + 1].MC - state.cells[i + 1].Mliqini) -
                                (state.cells[i].MC - state.cells[i].Mliqini));
    } else
        state.cells[i].transmassR = 0.;

    state.cells[i].transmassR /= state.cells[i].dx;
    state.cells[i + 1].transmassL = state.cells[i].transmassR;
    state.cells[i].FonteMudaFase = state.cells[i].transmassR;
    if (i == state.lastCell - 1) {
        state.cells[i + 1].transmassR = state.cells[i].transmassR;
        state.cells[i + 1].FonteMudaFase = state.cells[i].transmassR;
    }
}

void refreshProperties(const SteadyStateState &state) {
    for (int i = 0; i <= state.lastCell; i++) {
        state.cells[i].rpC = state.cells[i].rpCi =
            state.cells[i].flui.MasEspLiq(state.cells[i].pres, state.cells[i].temp);
        state.cells[i].rgC = state.cells[i].rgCi =
            state.cells[i].flui.MasEspGas(state.cells[i].pres, state.cells[i].temp);
        state.cells[i].rcC = state.cells[i].rcCi =
            state.cells[i].fluicol.MasEspFlu(state.cells[i].pres, state.cells[i].temp);
        state.cells[i].rpL = state.cells[i].rpLi =
            state.cells[i].flui.MasEspLiq(state.cells[i].presL, state.cells[i].tempL);
        state.cells[i].rgL = state.cells[i].rgLi =
            state.cells[i].flui.MasEspGas(state.cells[i].presL, state.cells[i].tempL);
        state.cells[i].rcL = state.cells[i].rcLi =
            state.cells[i].fluicol.MasEspFlu(state.cells[i].presL, state.cells[i].tempL);

        state.cells[i].mipC = state.cells[i].flui.ViscOleo(state.cells[i].pres, state.cells[i].temp);
        state.cells[i].migC = state.cells[i].flui.ViscGas(state.cells[i].pres, state.cells[i].temp);
        state.cells[i].micC = state.cells[i].fluicol.VisFlu(state.cells[i].pres, state.cells[i].temp);

        if (i > 0) {
            state.cells[i - 1].rpR = state.cells[i - 1].rpRi = state.cells[i].rpC;
            state.cells[i - 1].rgR = state.cells[i - 1].rgRi = state.cells[i].rgC;
            state.cells[i - 1].rcR = state.cells[i - 1].rcRi = state.cells[i].rcC;

            state.cells[i - 1].mipR = state.cells[i].mipC;
            state.cells[i - 1].migR = state.cells[i].migC;
            state.cells[i - 1].micR = state.cells[i].micC;
        }
        if (i == state.lastCell) {
            state.cells[i].rpR = state.cells[i].rpRi = state.cells[i].rpC;
            state.cells[i].rgR = state.cells[i].rgRi = state.cells[i].rgC;
            state.cells[i].rcR = state.cells[i].rcRi = state.cells[i].rcC;

            state.cells[i].mipR = state.cells[i].mipC;
            state.cells[i].migR = state.cells[i].migC;
            state.cells[i].micR = state.cells[i].micC;
        }
    }
}

void refreshSteadyThermalVelocities(const SteadyStateState &state) {
    for (int i = 0; i <= state.lastCell; i++) {
        double diameter = state.cells[i].duto.a;
        double area = 0.25 * M_PI * diameter * diameter;
        double alfmed = state.cells[i].alf;
        double betmed = state.cells[i].bet;
        double rhol = (1. - betmed) * state.cells[i].rpC + betmed * state.cells[i].rcC;
        double rhog = state.cells[i].rgC;
        double liquidSpecificHeat = (1. - betmed) * state.cells[i].flui.CalorLiq(state.cells[i].presini, state.cells[i].temp) + betmed * state.cells[i].fluicol.CalorLiq(state.cells[i].presini, state.cells[i].temp);
        double liquidSpecificHeatVolume = liquidSpecificHeat;
        double gasSpecificHeat = state.cells[i].flui.CalorGas(state.cells[i].presini, state.cells[i].temp);
        double gasSpecificHeatVolume = state.cells[i].flui.CalorGasVolMod(state.cells[i].presini, state.cells[i].temp, state.cells[i].rgC);
        double coefTempo = (rhol * (1 - alfmed) * liquidSpecificHeatVolume + rhog * alfmed * gasSpecificHeatVolume) * area;
        double ugsmed = state.cells[i].QG / area;
        double ulsmed = state.cells[i].QL / area;
        double coefdxT = (rhol * ulsmed * liquidSpecificHeat + rhog * ugsmed * gasSpecificHeat) * area;
        state.cells[i].VTemper = coefdxT / coefTempo;
    }
}

void computePseudoTransientTimeStep(const SteadyStateState &state) {
    state.timeStep = 100.;
    for (int i = 0; i <= state.lastCell; i++) {
        double jmix = 0.;
        double dtaux;
        if (fabs(state.cells[i].VTemper) > jmix)
            jmix = fabs(state.cells[i].VTemper);
        if (fabs(jmix) > 1e-5)
            dtaux = 1.0 * state.cells[i].dx / jmix;
        else
            dtaux = state.timeStep;
        if (dtaux < state.timeStep)
            state.timeStep = dtaux;
    }
    if (state.input.lingas == 1) {
        for (int i = 0; i <= state.gasCellCount; i++) {
            double dtaux;
            double jmix = fabs(state.gasCells[i].VGasR / state.gasCells[i].u1L);
            if (fabs(jmix) > 1e-5)
                dtaux = state.gasCells[i].dx0 / jmix;
            else
                dtaux = state.timeStep;
            if (dtaux < state.timeStep)
                state.timeStep = dtaux;
        }
    }
    state.timeStep = 0.8 * state.timeStep;
    for (int i = 0; i <= state.lastCell; i++)
        state.cells[i].dt = state.timeStep;
    if (state.input.lingas == 1) {
        for (int i = 0; i <= state.gasCellCount; i++)
            state.gasCells[i].dt = state.timeStep;
    }
}

void refreshUpstreamProductionPeriphery(const SteadyStateState &state, int i) {
    if (state.cells[i].presaux < 0.5)
        state.cells[i].presaux = 0.5;
    state.cells[i - 1].presauxR = state.cells[i].presaux;
    if (i < state.lastCell)
        state.cells[i + 1].presauxL = state.cells[i].presaux;

    state.cells[i - 1].dpB = 0.;
    state.cells[i - 1].potB = 0.;
    state.cells[i - 1].potBT = 0.;
    double tmed;
    if (state.steadyIteration != 0 && state.input.AceleraConvergPerm == 0)
        tmed = (state.cells[i].dx * state.cells[i].temp + state.cells[i].dxL * state.cells[i].tempL) / (state.cells[i].dx + state.cells[i].dxL);
    else
        tmed = state.cells[i - 1].temp;
    double sinalQ = 1.;
    if (fabs(state.cells[i - 1].QL + state.cells[i - 1].QG) > 1e-15)
        sinalQ = fabs(state.cells[i - 1].QL + state.cells[i - 1].QG) / (state.cells[i - 1].QL + state.cells[i - 1].QG);
    if (state.cells[i - 1].acsr.tipo == 4 && state.cells[i - 1].acsr.bcs.freqnova > 1.) {
        double vazmix = fabs(state.cells[i - 1].QL + state.cells[i - 1].QG);
        double alf0 = state.cells[i - 1].alf;
        double bet0 = state.cells[i - 1].bet;
        double rhomis = (alf0 * state.cells[i - 1].flui.MasEspGas(state.cells[i].presaux, tmed) + (1 - alf0) * ((1 - bet0) * state.cells[i - 1].flui.MasEspLiq(state.cells[i].presaux, tmed) + bet0 * state.cells[i - 1].fluicol.MasEspFlu(state.cells[i].presaux, tmed)));
        double vismis = alf0 * state.cells[i - 1].flui.ViscGas(state.cells[i].presaux, tmed) + (1 - alf0) * ((1. - bet0) * state.cells[i - 1].flui.ViscOleo(state.cells[i].presaux, tmed) + bet0 * state.cells[i].fluicol.VisFlu(state.cells[i].presaux, tmed));
        vazmix *= (86400 / 0.1589876);
        state.cells[i - 1].acsr.bcs.NovaVis(vismis, rhomis, vazmix);
        state.cells[i - 1].dpB = sinalQ * 0.3048 * state.cells[i - 1].acsr.bcs.Hvis * rhomis * 9.82;
        state.cells[i - 1].potB = state.cells[i - 1].acsr.bcs.Pvis * 745.7;
        state.cells[i - 1].potTermo = (1. - state.cells[i - 1].acsr.bcs.Evis / 100.) * state.cells[i - 1].potB;
        if (state.cells[i - 1].acsr.bcs.eficM > 0.)
            state.cells[i - 1].potBT = (1. + 100. * (1. - state.cells[i - 1].acsr.bcs.eficM / 100.) / state.cells[i - 1].acsr.bcs.eficM) * state.cells[i - 1].potB;
        else
            state.cells[i - 1].potBT = 0.;
        state.cells[i - 1].potTermo += state.cells[i - 1].potBT * (1. - state.cells[i - 1].acsr.bcs.eficM / 100.) * state.cells[i - 1].acsr.bcs.fracTermMotorEfic;

    } else if (state.cells[i - 1].acsr.tipo == 7) {
        state.cells[i - 1].dpB = sinalQ * state.cells[i - 1].acsr.delp * 98066.5;
        double bet0 = state.cells[i - 1].bet;
        double rhog = state.cells[i - 1].flui.MasEspGas(state.cells[i].presaux, tmed);
        double rhol = ((1 - bet0) * state.cells[i - 1].flui.MasEspLiq(state.cells[i].presaux, tmed) + bet0 * state.cells[i - 1].fluicol.MasEspFlu(state.cells[i].presaux, tmed));
        double qgMon = (state.cells[i - 1].MC - state.cells[i - 1].Mliqini) / rhog;
        double qlmon = state.cells[i - 1].Mliqini / rhol;
        double npoli = 1.;
        if (state.cells[i - 1].acsr.tipoCompGas == 0) {
            npoli = state.cells[i - 1].flui.ConstAdG(state.cells[i - 1].pres, state.cells[i - 1].temp);
            if ((npoli - 1.05) < 1e-2)
                npoli = 1.05;
        } else if (state.cells[i - 1].acsr.tipoCompGas == 1)
            npoli = state.cells[i - 1].acsr.fatPoli;
        double Wcomp;
        double Wbomb;
        Wbomb = (100. / state.cells[i - 1].acsr.eficLiq) * state.cells[i - 1].dpB * qlmon;
        double wcompiso = -state.cells[i].presaux * 98066.5 * qgMon * log(state.cells[i].presaux / (state.cells[i].presaux + state.cells[i - 1].acsr.delp));
        if (state.cells[i - 1].acsr.tipoCompGas != 2)
            Wcomp = -(state.cells[i].presaux * 98066.5 * qgMon / (1. - npoli)) *
                    (pow(1 + sinalQ * state.cells[i - 1].acsr.delp / state.cells[i].presaux, (npoli - 1.) / npoli) - 1.);
        else
            Wcomp = wcompiso;
        Wcomp *= (100. / state.cells[i - 1].acsr.eficGas);
        state.cells[i - 1].potB = sinalQ * (Wcomp + Wbomb);
        state.cells[i - 1].potTermo = sinalQ * ((1. - state.cells[i - 1].acsr.eficLiq / 100.) * Wbomb + (1. - state.cells[i - 1].acsr.eficGas / 100.) * Wcomp);
        state.cells[i - 1].potBT = state.cells[i - 1].potB;
    } else if (state.cells[i - 1].acsr.tipo == 17 && state.cells[i - 1].acsr.multibcs.freqnova > 1.) {
        double alf0 = state.cells[i - 1].alf;
        double bet0 = state.cells[i - 1].bet;
        state.cells[i - 1].acsr.multibcs.flui = state.cells[i - 1].flui;
        state.cells[i - 1].acsr.multibcs.fluicol = state.cells[i - 1].fluicol;
        state.cells[i - 1].acsr.multibcs.marchaMultiBcs(state.cells[i - 1].QG, state.cells[i - 1].QL, state.cells[i].presaux, tmed, alf0, bet0);
        state.cells[i - 1].dpB = state.cells[i - 1].acsr.multibcs.dpB * 98066.52;
        state.cells[i - 1].potB = state.cells[i - 1].acsr.multibcs.potBT;
        state.cells[i - 1].potBT = state.cells[i - 1].acsr.multibcs.potBT;
        state.cells[i - 1].potTermo = state.cells[i - 1].acsr.multibcs.potTermo;
        state.cells[i - 1].potTermo += state.cells[i - 1].potBT * (1. - state.cells[i - 1].acsr.multibcs.eficM / 100.) * state.cells[i - 1].acsr.multibcs.fracTermMotorEfic;
    }
}

void refreshDownstreamProductionPeriphery(const SteadyStateState &state, int i) {
    if (state.cells[i].pres < 0.1)
        state.cells[i].pres = 0.1;
    state.cells[i - 1].presR = state.cells[i].pres;
    if (i < state.lastCell) {
        state.cells[i + 1].presL = state.cells[i].pres;
        state.cells[i + 1].presLini = state.cells[i + 1].presL;
    }
    state.cells[i].presini = state.cells[i].pres;
}

void gasLineHydrostatic(const SteadyStateState &state) {
    double pchute = state.input.gasinj.presinj[0];
    state.gasCells[0].pres = pchute;
    double taux;
    taux = state.input.celg[0].textern;
    state.cells[0].temp = taux;
    double rhog = state.gasCells[0].flui.MasEspGas(state.gasCells[0].pres, taux);
    for (int i = 0; i < state.gasCellCount; i++) {
        taux = state.input.celg[i].textern;
        rhog = state.gasCells[i].flui.MasEspGas(pchute, taux);
        double dxmed = 0.5 * (state.gasCells[i].dx0 + state.gasCells[i + 1].dx0);
        pchute -= ((rhog * 9.81 * sinl(state.gasCells[i].duto.teta) * dxmed) / 98066.5);
        state.gasCells[i + 1].pres = pchute;
        state.gasCells[i + 1].temp = taux;
    }
}

double marchGasSteady(const SteadyStateState &state, double chutemass) {
    int valveCount = state.input.nvalvgas;
    double massGas = 0.;
    double erro = 10.;
    double erro1 = 10.;
    if (chutemass < 0) {
        for (int j = 0; j < valveCount; j++)
            massGas += state.gasCells[state.gasValveCellIndices[j]].massfonteCH; // caso nÃƒÂ£o se tenha um chute da vazao massica
        // na entrada da linha de injecao, faz-se o somatÃƒÂ³rio dos valores de vazao em cada VGL
    } else
        massGas = chutemass;
    int itera = 0;
    double relaxa = 0.5; // relaxacao para a iteracao da vazao total de injecao estimada a cada iteracao
    double presteste = -10;
    double presteste0 = -10;
    int injectionFlowRateIsLow = 0;
    while (((erro > 0.00001 && injectionFlowRateIsLow == 0) || erro1 > 0.00001) && itera < 600) { // iteracao de convergencia
        presteste0 = presteste;
        if (itera > 20)
            relaxa = 0.1;
        else if (itera > 50)
            relaxa = 0.05;
        else if (itera > 100)
            relaxa = 0.01;

        state.gasCells[0].presL = state.initialGasPressure;
        state.gasCells[0].pres = state.initialGasPressure;
        state.gasCells[0].presini = state.initialGasPressure;
        state.gasCells[1].presL = state.initialGasPressure;
        state.gasCells[0].tempL = state.initialGasTemperature;
        state.gasCells[0].temp = state.initialGasTemperature;
        state.gasCells[1].tempL = state.initialGasTemperature;
        state.gasCells[0].rg = state.gasCells[0].flui.MasEspGas(state.initialGasPressure, state.initialGasTemperature);
        state.gasCells[0].u1L = state.gasCells[0].duto.area * state.gasCells[0].rg;
        state.gasCells[0].u1LL = state.gasCells[0].u1L;
        state.gasCells[1].u1LL = state.gasCells[0].u1L;
        state.gasCells[0].VGasL = 0;
        state.gasCells[0].massfonteCH = massGas;
        state.gasCells[0].VGasR = massGas;
        state.gasCells[1].VGasL = massGas;
        for (int i = 1; i <= state.gasCellCount; i++) { // marcha
            state.updaters.updateSteadyGasPressure(i);            // avanco do valor de pressao de uma celula para outra, no centro da celula
            state.updaters.updateSteadyGasTemperature(i);            // avanco do valor de temperatura de uma celula para outra, centro da celula
            if (isnan(state.gasCells[i].temp))
                NumError("Temperatrura na linha de servico com valor NaN");
            state.updaters.computeSteadyGasFlowRate(i); // verifica se no centro desta celula tem uma VGL, calcula a vazao da VGL
            // retira este valor da vazÃƒÂ£o total na linha
            state.gasCells[i].rg = state.gasCells[i].flui.MasEspGas(state.gasCells[i].pres, state.gasCells[i].temp);
            state.gasCells[i - 1].rgR = state.gasCells[i].rg;
            state.gasCells[i].u1L = state.gasCells[i].duto.area * state.gasCells[i].rg;
            state.gasCells[i - 1].u1R = state.gasCells[i].u1L;
            if (i < state.gasCellCount)
                state.gasCells[i + 1].u1LL = state.gasCells[i].u1L;
        }
        state.gasCells[state.gasCellCount].rgR = state.gasCells[state.gasCellCount].rg;
        double massGas2 = massGas; // guarda o valor antigo de injecao de gas
        massGas = 0;
        for (int j = 0; j < valveCount; j++)
            massGas += state.gasCells[state.gasValveCellIndices[j]].massfonteCH; // atualiza a vazao de injecao a partir
        // das vazoes calculadas nas VGLs
        itera++;
        massGas = (relaxa * massGas + (1. - relaxa) * massGas2); // relaxacao da estimativa de vazao de injecao
        if (0.05 * state.gasCells[0].duto.area * state.gasCells[0].rg > massGas)
            injectionFlowRateIsLow = 1;
        else
            injectionFlowRateIsLow = 0;
        if (fabs(massGas) > 1e-15)
            erro = fabs(massGas - massGas2) / fabs(massGas); // erro na estimativa de vazao
        // de uma iteracao para a outra
        else if (fabs(massGas2) > 1e-15)
            erro = fabs(massGas - massGas2) / fabs(massGas2); // erro = erro/2.;
        else
            erro = 0.;
        presteste = state.gasCells[state.gasCellCount].pres;
        erro1 = fabs(presteste - presteste0) / presteste0; // erro no valor da pressao na ultima celula
    }
    if (itera >= 600) {
        if ((*state.globals).chaverede == 0 && state.input.AP == 0) {
            NumError("Busca de valores iniciais para calculo de zero de funcao em marchaGasPerm1 atingiu maximo de iteracoes");
            return 0.0;
        } else {
            if ((*state.globals).iterRede > 0)
                return -1.1e10;
            else
                return 1.1e10;
        }
    } else
        return erro1;
}

namespace {

/// Which fluid receives the dry-gas aware (six-argument) compositional flash in
/// the gas-source arm of the steady production marches.
///
/// It exists because marchProductionSteady and marchProductionSteadySecondary
/// disagreed about it and nothing in the code said why. Naming the disagreement
/// is not the same as resolving it: see H4 in evidencia/marchaprod-diff.md.
enum class DryGasFlashTarget {
    /// marchProductionSteady, marchReverseProductionSteady: the flag goes to
    /// the injected gas.
    sourceFluid,
    /// marchProductionSteadySecondary: the flag goes to the cell fluid.
    cellFluid,
};

void seedFirstCellVoidFraction(const SteadyStateState &state, double pchute, double &alfini, double &betini,
                                      DryGasFlashTarget dryGasFlashTarget) {
    // estimativa da fracao de vazio na primeira celula do sistema
    if (state.cells[0].acsr.tipo == 0) { // sem nenhuma fonte
        state.cells[0].temp = state.input.celp[0].textern;
        alfini = 1.;
        betini = 0.;
        if (state.input.flashCompleto == 2) {
            if (state.input.tabelaDinamica == 0)
                state.cells[0].flui.atualizaPropComp(pchute, state.cells[0].temp);
        }
    } else if (state.cells[0].acsr.tipo == 1) { // fonte de gas
        state.cells[0].temp = state.cells[0].acsr.injg.temp;
        if (state.input.flashCompleto == 2) {
            // The only place the three marches disagree. See H4 in
            // evidencia/marchaprod-diff.md: marchaProdPerm1 and
            // marchaProdPerm1Rev hand the dry-gas flag to the SOURCE fluid,
            // marchaProdPerm2 to the CELL fluid. Both forms are kept, on
            // purpose, because at most one of them can be right and this
            // refactoring is not the place to decide which.
            if (dryGasFlashTarget == DryGasFlashTarget::sourceFluid) {
                if (state.input.tabelaDinamica == 0)
                    state.cells[0].flui.atualizaPropComp(pchute, state.cells[0].temp);
                state.cells[0].acsr.injg.FluidoPro.atualizaPropComp(pchute, state.cells[0].temp, -1, NULL, NULL, state.cells[0].acsr.injg.seco);
            } else {
                if (state.input.tabelaDinamica == 0)
                    state.cells[0].flui.atualizaPropComp(pchute, state.cells[0].temp, -1, NULL, NULL, state.cells[0].acsr.injg.seco);
                state.cells[0].acsr.injg.FluidoPro.atualizaPropComp(pchute, state.cells[0].temp);
            }
        }
        if (state.cells[0].acsr.injg.seco == 1) {
            alfini = 1.;
            betini = 0.;
        } else {
            double masgas = state.cells[0].acsr.injg.VMas(pchute, state.cells[0].temp);
            double quality;
            if (state.input.flashCompleto != 2)
                quality = state.cells[0].acsr.injg.FluidoPro.FracMassHidra(1., 20.);
            else
                quality = state.cells[0].acsr.injg.FluidoPro.dStockTankVaporMassFraction;
            double masT = masgas / quality;
            quality = state.cells[0].acsr.injg.FluidoPro.FracMassHidra(pchute, state.cells[0].temp);
            double qgas = masT * quality /
                          state.cells[0].acsr.injg.FluidoPro.MasEspGas(pchute, state.cells[0].temp);
            double qliq = masT * (1. - quality) /
                          state.cells[0].acsr.injg.FluidoPro.MasEspLiq(pchute, state.cells[0].temp);
            double qcomp = state.cells[0].acsr.injg.razCompGas *
                           state.cells[0].acsr.injg.QGas * state.cells[0].acsr.injg.fluidocol.MasEspFlu(1., 20.) /
                           state.cells[0].acsr.injg.fluidocol.MasEspFlu(pchute, state.cells[0].temp);
            qcomp /= 86400.;
            alfini = qgas / (qliq + qcomp + qgas);
            if ((fabs(qcomp) + fabs(qliq)) > 1e-15)
                betini = fabs(qcomp) / (fabs(qcomp) + fabs(qliq));
            else
                betini = 0.;
        }
    } else if (state.cells[0].acsr.tipo == 2) { // fonte de liquido, faz-se uma estimativa a partir
        // da vazÃƒÂ£o volumÃƒÂ©trica das fases, fracao de vazio = fracao de vazio sem escorregamento
        state.cells[0].temp = state.cells[0].acsr.injl.temp;
        if (state.input.flashCompleto == 2) {
            if (state.input.tabelaDinamica == 0)
                state.cells[0].flui.atualizaPropComp(pchute, state.cells[0].temp);
            state.cells[0].acsr.injl.FluidoPro.atualizaPropComp(pchute, state.cells[0].temp);
        }
        double qgas = state.cells[0].acsr.injl.QLiq * (1 - state.cells[0].acsr.injl.bet) *
                      (1. - state.cells[0].acsr.injl.FluidoPro.BSW) *
                      (state.cells[0].acsr.injl.FluidoPro.RGO -
                       state.cells[0].acsr.injl.FluidoPro.rDgD * state.cells[0].acsr.injl.FluidoPro.RS(pchute, state.cells[0].temp) * 6.29 / 35.31467) *
                      state.cells[0].acsr.injl.FluidoPro.Deng * 1.225 / state.cells[0].acsr.injl.FluidoPro.MasEspGas(pchute, state.cells[0].temp);
        double qliq = state.cells[0].acsr.injl.QLiq * (1 - state.cells[0].acsr.injl.bet) *
                          (1. - state.cells[0].acsr.injl.FluidoPro.BSW) * state.cells[0].acsr.injl.FluidoPro.BOFunc(pchute, state.cells[0].temp) +
                      state.cells[0].acsr.injl.QLiq * (1 - state.cells[0].acsr.injl.bet) *
                          state.cells[0].acsr.injl.FluidoPro.BSW * state.cells[0].acsr.injl.FluidoPro.BAFunc(pchute, state.cells[0].temp) +
                      state.cells[0].acsr.injl.QLiq * state.cells[0].acsr.injl.bet;
        alfini = qgas / (qliq + qgas);
        betini = state.cells[0].acsr.injl.bet;
    } else if (state.cells[0].acsr.tipo == 3) { // IPR no inicio da tubulacao
        // da mesma maneira que no caso de fonte de liquido, fracao de vazio= sem escorregamento
        state.cells[0].temp = state.cells[0].acsr.ipr.Tres;
        if (state.input.flashCompleto == 2) {
            if (state.input.tabelaDinamica == 0)
                state.cells[0].flui.atualizaPropComp(pchute, state.cells[0].temp);
            state.cells[0].acsr.ipr.FluidoPro.atualizaPropComp(pchute, state.cells[0].temp);
        }
        double qgas = state.cells[0].acsr.ipr.MasG(pchute, state.cells[0].temp) /
                      state.cells[0].acsr.ipr.FluidoPro.MasEspGas(pchute, state.cells[0].temp);
        double qliq = state.cells[0].acsr.ipr.MasL(pchute, state.cells[0].temp) /
                      state.cells[0].acsr.ipr.FluidoPro.MasEspLiq(pchute, state.cells[0].temp);
        alfini = qgas / (qliq + qgas);
        betini = 0.;
    } else if (state.cells[0].acsr.tipo == 15) { // IPR no inicio da tubulacao
        // da mesma maneira que no caso de fonte de liquido, fracao de vazio= sem escorregamento
        state.cells[0].temp = state.cells[0].acsr.radialPoro.tRes;
        state.cells[0].pres = pchute;
        if (state.input.flashCompleto == 2) {
            if (state.input.tabelaDinamica == 0)
                state.cells[0].flui.atualizaPropComp(pchute, state.cells[0].temp);
            state.cells[0].acsr.radialPoro.flup.atualizaPropComp(pchute, state.cells[0].temp);
        }
        state.updaters.updateSource(0);
        double qgas = state.cells[0].acsr.radialPoro.fluxIniG /
                      state.cells[0].acsr.radialPoro.flup.MasEspGas(pchute, state.cells[0].temp);
        double qliq = (state.cells[0].acsr.radialPoro.fluxIni + state.cells[0].acsr.radialPoro.fluxIniA) /
                      state.cells[0].acsr.radialPoro.flup.MasEspLiq(pchute, state.cells[0].temp);
        alfini = qgas / (qliq + qgas);
        betini = 0.;
    } else if (state.cells[0].acsr.tipo == 16) { // IPR no inicio da tubulacao
        // da mesma maneira que no caso de fonte de liquido, fracao de vazio= sem escorregamento
        state.cells[0].temp = state.cells[0].acsr.poroso2D.dados.tRes;
        state.cells[0].pres = pchute;
        if (state.input.flashCompleto == 2) {
            if (state.input.tabelaDinamica == 0)
                state.cells[0].flui.atualizaPropComp(pchute, state.cells[0].temp);
            state.cells[0].acsr.poroso2D.dados.flup.atualizaPropComp(pchute, state.cells[0].temp);
        }
        state.updaters.updateSource(0);
        double qgas = state.cells[0].acsr.poroso2D.dados.transfer.fluxIniG /
                      state.cells[0].acsr.poroso2D.dados.flup.MasEspGas(pchute, state.cells[0].temp);
        double qliq = (state.cells[0].acsr.poroso2D.dados.transfer.fluxIni + state.cells[0].acsr.poroso2D.dados.transfer.fluxIniA) /
                      state.cells[0].acsr.poroso2D.dados.flup.MasEspLiq(pchute, state.cells[0].temp);
        alfini = qgas / (qliq + qgas);
        betini = 0.;
    } else if (state.cells[0].acsr.tipo == 10) { // fonte de massa no inicio da tubulacao
        // da mesma maneira que no caso de fonte de liquido, fracao de vazio= sem escorregamento
        state.cells[0].temp = state.cells[0].acsr.injm.temp;
        if (state.input.flashCompleto == 2) {
            if (state.input.tabelaDinamica == 0)
                state.cells[0].flui.atualizaPropComp(pchute, state.cells[0].temp);
            state.cells[0].acsr.injm.FluidoPro.atualizaPropComp(pchute, state.cells[0].temp);
        }
        if (state.cells[0].acsr.injm.condTermo == 0) {
            state.cells[0].pres = pchute;
            state.updaters.updateSource(0);
        }
        double qgas = state.cells[0].acsr.injm.MassG /
                      state.cells[0].acsr.injm.FluidoPro.MasEspGas(pchute, state.cells[0].temp);
        double qliq = state.cells[0].acsr.injm.MassP /
                          state.cells[0].acsr.injm.FluidoPro.MasEspLiq(pchute, state.cells[0].temp) +
                      state.cells[0].acsr.injm.MassC /
                          state.cells[0].acsr.injm.fluidocol.MasEspFlu(pchute, state.cells[0].temp);
        alfini = qgas / (qliq + qgas);
        betini = 0.;
    } else { // se nenhuma das opcoes, fracao de vazio=1
        state.cells[0].temp = state.input.celp[0].textern;
        if (state.input.flashCompleto == 2) {
            if (state.input.tabelaDinamica == 0)
                state.cells[0].flui.atualizaPropComp(pchute, state.cells[0].temp);
        }
        alfini = 1.;
        betini = 0.;
    }
    if (state.cells[0].acsr.tipo == 15) {
        state.cells[0].flui.BSW = state.cells[0].acsr.radialPoro.BSW;
    } else if (state.cells[0].acsr.tipo == 16) {
        state.cells[0].flui.BSW = state.cells[0].acsr.poroso2D.dados.transfer.BSW;
    }
}

/// Solves the gas line at steady state under its inlet condition. With the
/// injection pressure given, marches the line when the injection choke is not
/// throttling, and searches its pressure (tertiary) when it is. With the flow
/// rate given, searches its pressure (secondary).
void solveGasLineSteady(const SteadyStateState &state, InjectionPressureCondition) {
    // marcha para o caso, pressao de injecao
    if (state.input.chokes.abertura[0] >= 0.2) { // choke de injecao inativo
        for (int iter = 0; iter < 1; iter++) {
            marchGasSteady(state);
        }
    } else
        state.updaters.searchGasPressureSteadyTertiary(); // choke de injecao ativo
}

void solveGasLineSteady(const SteadyStateState &state, InjectionFlowRateCondition) {
    state.updaters.searchGasPressureSteadySecondary(); // marcha na linha de gas para o caso de vazao de injecao
}

/// Before the first march, a pressure condition lets the gas-line pressure at
/// the valves be estimated by hydrostatics, and the valves' injection computed
/// from it. Under a flow-rate condition there is nothing to estimate: the valves
/// share the given rate.
void estimateValveGasPressure(const SteadyStateState &state, InjectionPressureCondition) {
    gasLineHydrostatic(state);
}

void estimateValveGasPressure(const SteadyStateState &, InjectionFlowRateCondition) {}

void marchGasLineAndCoupleAnnulus(const SteadyStateState &state, double pchute) {
    if (pchute > 0 && state.input.lingas > 0 && state.input.nvalvgas > 0) {
        withGasInletCondition(state.gasCells[0].tipoCC,
                              [&](auto condition) { solveGasLineSteady(state, condition); });
    }

    if (state.input.acopColAnulPermForte > 0 && state.input.lingas > 0 && state.thermalSourceDisabled == 0) {
        refreshProperties(state);
        refreshSteadyThermalVelocities(state);
        computePseudoTransientTimeStep(state);
        state.updaters.connectTubing();
        for (int kontaPseudo = 0; kontaPseudo < state.input.acopColAnulPermForte; kontaPseudo++) {
            for (int iterm = 0; iterm <= state.gasCellCount; iterm++)
                state.gasCells[iterm].tempini = state.gasCells[iterm].temp;
            for (int iterm = 1; iterm <= state.gasCellCount; iterm++) {
                state.updaters.computeGasTemperature(iterm, state.gasCells[iterm - 1].tempini, 1);
            }
            state.cells[0].tempini = state.cells[0].temp;
            state.cells[1].tempLini = state.cells[1].tempL;
            state.cells[1].tempL = state.cells[0].temp;
            for (int iterm = 1; iterm <= state.lastCell; iterm++) {
                state.cells[iterm].tempini = state.cells[iterm].temp;
            }
            for (int iterm = 1; iterm <= state.lastCell; iterm++) {
                state.updaters.computeTemperature(iterm, state.cells[iterm].tempini, 1);
            }
            refreshProperties(state);
            computePseudoTransientTimeStep(state);
            state.updaters.connectTubing();
        }
    }
}

bool advanceProductionCells(const SteadyStateState &state, double pchute, int &i, double &abortValue) {
    while (i <= state.lastCell && state.cells[i - 1].pres >= 0.1 && fabs(pchute - state.cells[0].pres) < (*state.globals).localtiny) {

        advanceUpstreamSteadyPressure(state, i, 0); // avanco da marcha para obter a pressao na fronteira esquerda
        // da celula i
        // teste para ver se ocorreu algum problema:
        if (state.input.usaTabela == 1 && (state.input.tabent.pmax - state.cells[i].presaux) < (*state.globals).localtiny)
            {
                abortValue = 1e10;
                return true;
            }
        if (state.cells[i].presaux <= 0.1 ||
            (state.input.usaTabela == 1 && (state.input.tabent.pmin - state.cells[i].presaux) > (*state.globals).localtiny)) {
            {
                abortValue = -1e10;
                return true;
            }
        }
        refreshUpstreamProductionPeriphery(state, i); // atualizacao da pressao da fronteira esquerda,
        // caso exista alguma BCS ou incremento de pressao
        if (i == 312) {
            int para;
            para = 0;
        }
        if (state.input.flashCompleto != 2)
            advanceSteadyMass(state, i); // verifica se existe alguma fonte na celula anterior, com isto, atualiza
        // as vazoes massica na fronteira a esquerda, alÃ©m das propriedades dos fluidos,
        // densidade do gas, RGO, API, BSW, beta
        else
            advanceCompositionalSteadyMass(state, i);

        if (state.input.acopColAnulPermForte == 0 || state.input.lingas == 0 || state.convergenceMonitor > 0.3)
            state.updaters.advanceSteadyTemperature(i, 0); // faz o avanco da temperatura, da celula i-1 para a celula i
        // verifica se teve algum problema nos limites de temperatura
        // caso se esteja trabalhando com tabela PVTSim
        if (state.input.usaTabela == 1 && (state.cells[i].temp - state.input.tabent.tmin) < (*state.globals).localtiny)
            state.cells[i].temp = state.input.tabent.tmin;
        state.updaters.updateProductionTemperaturePeriphery(i); // mera atualizacao de atributos de temperatura a esquerda e a direita
        advanceDownstreamSteadyPressure(state, i, 0); // evolui a pressao  fronteira a esquerda da celula i para o
        // seu centro de celula
        // verifica se ocorreu algum problema nesta evolucao de de pressao no centro da
        // celula
        if (state.input.usaTabela == 1 && (state.input.tabent.pmax - state.cells[i].pres) < (*state.globals).localtiny) {
            {
                abortValue = 1e10;
                return true;
            }
        }
        if (state.cells[i].pres <= 0.1 ||
            (state.input.usaTabela == 1 && (state.input.tabent.pmin - state.cells[i].pres) > (*state.globals).localtiny)) {
            {
                abortValue = -1e10;
                return true;
            }
        }
        refreshDownstreamProductionPeriphery(state, i); // mera atualizacao de atributos que guardam valores de pressao
        // das celulas a esquerda e a direita
        for (int j = 0; j < state.input.nvalvgas; j++) { // reavaliacao da vazao da valvula de gas lift, quando
            // a celula tem uma.
            // P.S. parece uma acao desnecessÃ¡ria e talvez atÃ© um complicador
            // densecessario, em vavliacao
            if (state.productionValveCellIndices[j] == i) {
                int k = state.gasValveCellIndices[j];
                state.updaters.computeSteadyGasFlowRate(k);
            }
        }
        if (state.input.tipoFluido == 0)
            advanceSteadyMassTransfer(state, i - 1);
        else
            advanceSteadyGasMassTransfer(state, i - 1); // caso seja uma tabela PVTSim, calcula-se a
        // taxa de transferÃªncia de massa entre as fases para o uso no calculo de
        // calor latente da equacao de energia
        if (state.input.ordperm > 1) { // correcao de segunda ordem
            double D0presaux = state.cells[i].presaux - state.cells[i - 1].pres;
            double D0pres = state.cells[i].pres - state.cells[i].presaux;
            double D0temp = state.cells[i].temp - state.cells[i - 1].temp;
            advanceUpstreamSteadyPressure(state, i, 1);
            refreshUpstreamProductionPeriphery(state, i);
            advanceSteadyMass(state, i);
            state.updaters.advanceSteadyTemperature(i, 1);
            if (isnan(state.cells[i].temp)) {
                if (state.input.transiente == 0 && state.input.AP == 0)
                    // neste caso, se finaliza a simulacao, nao tem um transiente
                    // a ser feito a seguir e nem se estÃ¡ em uma rede
                    NumError(
                        "Temperatrura na linha de producao com valor NaN em marchaProdPerm1");
                else {
                    // apresenta apenas um aviso
                    cout << "#################PERMANENTE FALHOU EM SUA CONVERGENCIA##############################" << endl;
                    // se for em uma iteracao de rede, apos a primeira iteracao
                    if ((*state.globals).iterRede > 0)
                        {
                            abortValue = -1.1e10;
                            return true;
                        }
                    // se logo apos tem uma simulacao transiente ou se esta na primeira iteracao de rede
                    else
                        {
                            abortValue = 1.1e10;
                            return true;
                        }
                }
            }
            advanceDownstreamSteadyPressure(state, i, 1);
            state.cells[i].pres = 0.5 * (state.cells[i].presaux + D0pres + state.cells[i].pres);
            state.cells[i].presaux = 0.5 * (state.cells[i - 1].pres + D0presaux + state.cells[i].presaux);
            state.cells[i].temp = 0.5 * (state.cells[i - 1].temp + D0temp + state.cells[i].temp);
            refreshDownstreamProductionPeriphery(state, i);
            refreshUpstreamProductionPeriphery(state, i);
            state.updaters.updateProductionTemperaturePeriphery(i);
            advanceSteadyMass(state, i);
            if (state.input.tipoFluido == 0)
                advanceSteadyMassTransfer(state, i - 1);
            else
                advanceSteadyGasMassTransfer(state, i - 1);
        }

        // apÃ³s se atingir a pressao no centro da celula i, primeira iteracao de marcha
        // verifica-se se existe uma VGL em i e faz-se uma estimativa inicial da Vazao de
        // GL (caso exista linha de gas). Observar que isto sÃ³ Ã© feito para a iteracao zero.
        if (state.input.lingas > 0 && state.input.nvalvgas > 0 && state.steadyIteration == 0 && state.convergenceMonitor > 0.1)
            state.updaters.initializeSteadyValveGasFlowRate(i);
        i++;

        if (isnan(state.cells[i - 1].pres) || isnan(state.cells[i - 1].temp) || isnan(state.cells[i - 1].alf)) {
            {
                abortValue = 1e10;
                return true;
            }
        }

        // teste para verificar se a pressao do centro de celula ficou acima
        // da pressao estatica de uma eventual IPR
        if (state.cells[i - 1].pres <= 0.1 ||
            (state.input.usaTabela == 1 && (state.input.tabent.pmin - state.cells[i - 1].pres) > (*state.globals).localtiny)) {
            {
                abortValue = -1e10;
                return true;
            }
        } else if ((state.cells[i - 1].acsr.tipo == 3 &&
                    ((state.cells[i - 1].acsr.ipr.Pres - state.cells[i - 1].pres) < (*state.globals).localtiny) && i == 1)) {
            {
                abortValue = 1e10;
                return true;
            }
        }
    }
    return false;
}

bool advanceReverseProductionCells(const SteadyStateState &state, double pchute, int &i, double &abortValue) {
    while (i <= state.lastCell && state.cells[i - 1].pres >= 0.1 && fabs(pchute - state.cells[0].pres) < (*state.globals).localtiny) {

        advanceUpstreamSteadyPressure(state, i, 0); // avanco da marcha para obter a pressao na fronteira esquerda
        // da celula i
        // teste para ver se ocorreu algum problema:
        if (state.input.usaTabela == 1 && (state.input.tabent.pmax - state.cells[i].presaux) < (*state.globals).localtiny)
            {
                abortValue = 1e10;
                return true;
            }
        if (state.cells[i].presaux <= 0.1 ||
            (state.input.usaTabela == 1 && (state.input.tabent.pmin - state.cells[i].presaux) > (*state.globals).localtiny)) {
            {
                abortValue = -1e10;
                return true;
            }
        }
        refreshUpstreamProductionPeriphery(state, i); // atualizacao da pressao da fronteira esquerda,
        // caso exista alguma BCS ou incremento de pressao
        if (state.input.flashCompleto != 2)
            advanceSteadyMass(state, i); // verifica se existe alguma fonte na celula anterior, com isto, atualiza
        // as vazoes massica na fronteira a esquerda, alÃ©m das propriedades dos fluidos,
        // densidade do gas, RGO, API, BSW, beta
        else
            advanceCompositionalSteadyMass(state, i);
        // RenovaTempPerm(i, 0);//faz o avanco da temperatura, da celula i-1 para a celula i
        // verifica se teve algum problema nos limites de temperatura
        // caso se esteja trabalhando com tabela PVTSim
        if (isnan(state.cells[i].temp))
            NumError("Temperatrura na linha de producao com valor NaN");
        if (state.input.usaTabela == 1 && (state.cells[i].temp - state.input.tabent.tmin) < (*state.globals).localtiny)
            state.cells[i].temp = state.input.tabent.tmin;
        state.updaters.updateProductionTemperaturePeriphery(i); // mera atualizacao de atributos de temperatura a esquerda e a direita
        advanceDownstreamSteadyPressure(state, i, 0); // evolui a pressao  fronteira a esquerda da celula i para o
        // seu centro de celula
        // verifica se ocorreu algum problema nesta evolucao de de pressao no centro da
        // celula
        if (state.input.usaTabela == 1 && (state.input.tabent.pmax - state.cells[i].pres) < (*state.globals).localtiny) {
            {
                abortValue = 1e10;
                return true;
            }
        }
        if (state.cells[i].pres <= 0.1 ||
            (state.input.usaTabela == 1 && (state.input.tabent.pmin - state.cells[i].pres) > (*state.globals).localtiny)) {
            {
                abortValue = -1e10;
                return true;
            }
        }
        refreshDownstreamProductionPeriphery(state, i); // mera atualizacao de atributos que guardam valores de pressao
        // das celulas a esquerda e a direita
        if (state.input.tipoFluido == 0)
            advanceSteadyMassTransfer(state, i - 1);
        else
            advanceSteadyGasMassTransfer(state, i - 1); // caso seja uma tabela PVTSim, calcula-se a
        // taxa de transferÃªncia de massa entre as fases para o uso no calculo de
        // calor latente da equacao de energia
        if (state.input.ordperm > 1) { // correcao de segunda ordem
            double D0presaux = state.cells[i].presaux - state.cells[i - 1].pres;
            double D0pres = state.cells[i].pres - state.cells[i].presaux;
            double D0temp = state.cells[i].temp - state.cells[i - 1].temp;
            advanceUpstreamSteadyPressure(state, i, 1);
            refreshUpstreamProductionPeriphery(state, i);
            advanceSteadyMass(state, i);
            advanceDownstreamSteadyPressure(state, i, 1);
            state.cells[i].pres = 0.5 * (state.cells[i].presaux + D0pres + state.cells[i].pres);
            state.cells[i].presaux = 0.5 * (state.cells[i - 1].pres + D0presaux + state.cells[i].presaux);
            state.cells[i].temp = 0.5 * (state.cells[i - 1].temp + D0temp + state.cells[i].temp);
            refreshDownstreamProductionPeriphery(state, i);
            refreshUpstreamProductionPeriphery(state, i);
            state.updaters.updateProductionTemperaturePeriphery(i);
            advanceSteadyMass(state, i);
            if (state.input.tipoFluido == 0)
                advanceSteadyMassTransfer(state, i - 1);
            else
                advanceSteadyGasMassTransfer(state, i - 1);
        }

        i++;

        if (state.cells[i - 1].pres <= 0.1 ||
            (state.input.usaTabela == 1 && (state.input.tabent.pmin - state.cells[i - 1].pres) > (*state.globals).localtiny)) {
            {
                abortValue = -1e10;
                return true;
            }
        }
    }
    return false;
}

bool advanceProductionCellsSecondary(const SteadyStateState &state, double pchute, int &i, double &abortValue) {
    while (i <= state.lastCell && state.cells[i - 1].pres >= 0.1 && fabs(pchute - state.cells[0].pres) < (*state.globals).localtiny) {

        advanceUpstreamSteadyPressure(state, i, 0); // avanco da marcha para obter a pressao na fronteira esquerda
        // da celula i
        // teste para ver se ocorreu algum problema:
        if (state.input.usaTabela == 1 && (state.input.tabent.pmax - state.cells[i].presaux) < (*state.globals).localtiny)
            {
                abortValue = 1e10;
                return true;
            }
        refreshUpstreamProductionPeriphery(state, i); // atualizacao da pressao da fronteira esquerda,
        // caso exista alguma BCS ou incremento de pressao
        if (state.input.flashCompleto != 2)
            advanceSteadyMass(state, i); // verifica se existe alguma fonte na celula anterior, com isto, atualiza
        // as vazoes massica na fronteira a esquerda, alÃ©m das propriedades dos fluidos,
        // densidade do gas, RGO, API, BSW, beta
        else
            advanceCompositionalSteadyMass(state, i);

        if (state.input.acopColAnulPermForte == 0 || state.input.lingas == 0 || state.convergenceMonitor > 0.3)
            state.updaters.advanceSteadyTemperature(i, 0); // faz o avanco da temperatura, da celula i-1 para a celula i
        // verifica se teve algum problema nos limites de temperatura
        // caso se esteja trabalhando com tabela PVTSim
        if (isnan(state.cells[i].temp)) {
            if (state.input.transiente == 0 && state.input.AP == 0)
                // neste caso, se finaliza a simulacao, nao tem um transiente
                // a ser feito a seguir e nem se estÃ¡ em uma rede
                NumError(
                    "Temperatrura na linha de producao com valor NaN em marchaProdPerm2");
            else {
                // apresenta apenas um aviso
                cout << "#################PERMANENTE FALHOU EM SUA CONVERGENCIA##############################" << endl;
                // se for em uma iteracao de rede, apos a primeira iteracao
                if ((*state.globals).iterRede > 0)
                    {
                        abortValue = -1.1e10;
                        return true;
                    }
                // se logo apos tem uma simulacao transiente ou se esta na primeira iteracao de rede
                else
                    {
                        abortValue = 1.1e10;
                        return true;
                    }
            }
        }
        if (state.input.usaTabela == 1 && (state.cells[i].temp - state.input.tabent.tmin) < (*state.globals).localtiny)
            state.cells[i].temp = state.input.tabent.tmin;
        state.updaters.updateProductionTemperaturePeriphery(i); // mera atualizacao de atributos de temperatura a esquerda e a direita
        advanceDownstreamSteadyPressure(state, i, 0); // evolui a pressao  fronteira a esquerda da celula i para o
        // seu centro de celula
        // verifica se ocorreu algum problema nesta evolucao de de pressao no centro da
        // celula
        if (state.input.usaTabela == 1 && (state.input.tabent.pmax - state.cells[i].pres) < (*state.globals).localtiny)
            {
                abortValue = 1e10;
                return true;
            }
        refreshDownstreamProductionPeriphery(state, i); // mera atualizacao de atributos que guardam valores de pressao
        // das celulas a esquerda e a direita
        for (int j = 0; j < state.input.nvalvgas; j++) { // reavaliacao da vazao da valvula de gas lift, quando
            // a celula tem uma.
            // P.S. parece uma acao desnecessÃ¡ria e talvez atÃ© um complicador
            // densecessario, em vavliacao
            if (state.productionValveCellIndices[j] == i) {
                int k = state.gasValveCellIndices[j];
                state.updaters.computeSteadyGasFlowRate(k);
            }
        }
        if (state.input.tipoFluido == 0)
            advanceSteadyMassTransfer(state, i - 1);
        else
            advanceSteadyGasMassTransfer(state, i - 1); // caso seja uma tabela PVTSim, calcula-se a
        // taxa de transferÃªncia de massa entre as fases para o uso no calculo de
        // calor latente da equacao de energia
        if (state.input.ordperm > 1) { // correcao de segunda ordem
            double D0presaux = state.cells[i].presaux - state.cells[i - 1].pres;
            double D0pres = state.cells[i].pres - state.cells[i].presaux;
            double D0temp = state.cells[i].temp - state.cells[i - 1].temp;
            advanceUpstreamSteadyPressure(state, i, 1);
            refreshUpstreamProductionPeriphery(state, i);
            advanceSteadyMass(state, i);
            state.updaters.advanceSteadyTemperature(i, 1);
            advanceDownstreamSteadyPressure(state, i, 1);
            state.cells[i].pres = 0.5 * (state.cells[i].presaux + D0pres + state.cells[i].pres);
            state.cells[i].presaux = 0.5 * (state.cells[i - 1].pres + D0presaux + state.cells[i].presaux);
            state.cells[i].temp = 0.5 * (state.cells[i - 1].temp + D0temp + state.cells[i].temp);
            refreshDownstreamProductionPeriphery(state, i);
            refreshUpstreamProductionPeriphery(state, i);
            state.updaters.updateProductionTemperaturePeriphery(i);
            advanceSteadyMass(state, i);
            if (state.input.tipoFluido == 0)
                advanceSteadyMassTransfer(state, i - 1);
            else
                advanceSteadyGasMassTransfer(state, i - 1);
        }
        // apÃ³s se atingir a pressao no centro da celula i, primeira iteracao de marcha
        // verifica-se se existe uma VGL em i e faz-se uma estimativa inicial da Vazao de
        // GL (caso exista linha de gas). Observar que isto sÃ³ Ã© feito para a iteracao zero.
        if (state.input.lingas > 0 && state.input.nvalvgas > 0 && state.steadyIteration == 0 && state.convergenceMonitor > 0.1)
            state.updaters.initializeSteadyValveGasFlowRate(i);
        i++;
        // teste para verificar se a pressao do centro de celula ficou acima
        // da pressao estatica de uma eventual IPR ou ficou baixa demais
        if (state.cells[i - 1].pres <= 0.1 ||
            (state.input.usaTabela == 1 && (state.input.tabent.pmin - state.cells[i - 1].pres) > (*state.globals).localtiny))
            {
                abortValue = -1e10;
                return true;
            }
        else if (state.cells[i - 1].acsr.tipo == 3 && ((state.cells[i - 1].acsr.ipr.Pres - state.cells[i - 1].pres) < 1e-15 && i == 1))
            {
                abortValue = 1e10;
                return true;
            }
        else if (state.cells[i - 1].acsr.tipo == 15 && ((state.cells[i - 1].acsr.radialPoro.pRes[0] - state.cells[i - 1].pres) < 1e-15 && i == 1))
            {
                abortValue = 1e10;
                return true;
            }
        else if (state.cells[i - 1].acsr.tipo == 16 && ((state.cells[i - 1].acsr.poroso2D.dados.pRes - state.cells[i - 1].pres) < 1e-15 && i == 1))
            {
                abortValue = 1e10;
                return true;
            }
        else if (state.input.usaTabela == 1) {
            if ((state.input.tabent.pmax - state.cells[i - 1].pres) < (*state.globals).localtiny)
                {
                    abortValue = 1e10;
                    return true;
                }
        }
    }
    return false;
}

double surfaceChokeMassFlowRate(const SteadyStateState &state) {
    double maxSup = 0.;
    if (state.annulusDrift != 0 && state.cells[state.lastCell].pres > state.gasSurfacePressure) { // se for o anel de GL, a vazao no final deve ser zero
        double tESup = state.cells[state.lastCell].temp;
        double alfSup = state.cells[state.lastCell].alf;
        double betSup = state.cells[state.lastCell].bet;

        double masentrada = state.cells[state.lastCell - 1].MR;
        double massgas = state.cells[state.lastCell - 1].MR - state.cells[state.lastCell - 1].MliqiniR;

        double quality;
        state.outletPressure = state.cells[state.lastCell].pres;
        double rholp = state.cells[state.lastCell].flui.MasEspLiq(state.cells[state.lastCell].pres, state.cells[state.lastCell].temp);
        double rholc = state.cells[state.lastCell].fluicol.MasEspFlu(state.cells[state.lastCell].pres, state.cells[state.lastCell].temp);
        quality = fabs(massgas / masentrada);

        double masChk;

        double ypres = state.gasSurfacePressure / state.outletPressure;
        masChk = state.surfaceChoke.vazmassSachd(ypres, state.outletPressure, tESup, alfSup, betSup, quality, state.cells[state.lastCell - 1].flui,
                                       state.cells[state.lastCell - 1].fluicol);
        maxSup = state.surfaceChoke.vazmaxSachd(state.outletPressure, tESup, alfSup, betSup, quality, state.cells[state.lastCell - 1].flui, state.cells[state.lastCell - 1].fluicol);
        if (fabs(ypres) > fabs(state.surfaceChoke.razpres))
            maxSup = masChk;
        // maxSup Ã© a vazao total passando pelo choke

        if (state.surfaceChoke.AreaGarg > (1e-3) * state.cells[state.lastCell - 1].duto.area && ypres < 1.) {
            double cplM = (1. - betSup) * state.cells[state.lastCell].flui.CalorLiq(state.outletPressure, tESup) -
                          betSup * state.cells[state.lastCell].fluicol.CalorLiq(state.outletPressure, tESup);
            double jtlM = (1. - betSup) * state.cells[state.lastCell].flui.JTL(state.outletPressure, tESup) - betSup / rholc;
            double gasSpecificHeat = state.cells[state.lastCell].flui.CalorGas(state.outletPressure, tESup);
            double jtgM = state.cells[state.lastCell].flui.JTG(state.outletPressure, tESup);
            state.input.valTempChokeJus = tESup + ((1. - quality) * jtlM / cplM + quality * jtgM / gasSpecificHeat) * (state.gasSurfacePressure - state.outletPressure) * 98066.52;
        }

    } else {
        maxSup = 0.;
        state.input.valTempChokeJus = state.cells[state.lastCell].temp;
    }

    return maxSup;
}
}  // namespace


double marchProductionSteady(const SteadyStateState &state, double pchute) {

    int guessNeedsCorrection = 1;
    double alfini = 0.;
    double betini = 0.;

    seedFirstCellVoidFraction(state, pchute, alfini, betini, DryGasFlashTarget::sourceFluid);
    if (fabs(alfini) < 1e-6)
        alfini = 0.;
    if (fabs(betini) < 1e-6)
        betini = 0.;

    // esta marcha e feita para quando se tem alguma fonte no inicio da tubulacao,
    // portanto, admite-se que o duto esta fechado e coloca-se uma fonte no centro da
    // primeira celula. As vazoes na fronteira esquerda da celula sÃ£o portanto = 0
    state.cells[0].tempL = state.cells[0].temp;
    state.cells[1].tempL = state.cells[0].temp;
    state.cells[0].tempini = state.cells[0].temp;

    state.cells[0].ML = 0.;
    state.cells[0].MC = 0.;
    state.cells[1].ML = 0.;
    state.cells[0].MliqiniL = 0.;
    state.cells[0].Mliqini = 0.;
    state.cells[1].MliqiniL = 0.;
    state.cells[0].MComp = 0.;
    state.cells[0].QLL = 0.;
    state.cells[0].QL = 0;
    state.cells[1].QLL = 0.;
    state.cells[0].QG = 0.;
    double masfim = 1.;
    double masfim0 = -10;
    double presteste = 1.;
    double presteste0 = -10;
    state.steadyIteration = 0;

    int limIter = 2; // limite de iteracoes quando a opcao de aceleracao da convergencia esta desligado. desaconselhavel
    // desligar esta opcao, os ganhos sao pouco e a convergencia se torna instavel, principalmente
    // quando o tramo faz perte de um sistema de redes
    if (state.input.AceleraConvergPerm == 1) { // opcao aceleracao de convergencia ligada
        limIter = 1;                   // em geral faz-se apenas duas iteracoes de marcha para um determinado chute
        // Here the original went on to test for a pressure condition on the gas
        // line and set limIter to 1 again: the same value, so the test selected
        // nothing, and SC-012 removed it (evidencia/anomalias.md, A9-02). The
        // note that went with it, kept below, describes an extra march that
        // only searchProductionBottomHolePressureTertiary makes.
        //
        // no caso de se ter
        // uma condicao de contorno na injecao de gas = pressao, observou-se que o acoplamento dinamico
        // entre a linha de gas e de producao e mais difoicil, para se conseguir um sistema
        // melhor acoplado, deve-se fazer uma marcha iterativa a mais
    }
    // arq.CriterioConvergPerm Ã© um criterio de convergencia da marcha, so faz sentido
    // quando a aceleracao de convergencia esta desligada. Observe que esta convergencia nÃ£o
    // e de fato a conevregencia do problema, e apenas um repeticao de marcha para um determinado
    // chute de pressao ou de vazao no inicio da tubulacao. O que a convergencia busca de fato e
    // determinar qual a pressao ou vazao de fundo que satisfaz as condicoes de contorno no fim
    // da tubulacao, esta busca e feita nos metodos de busca.
    while ((fabs(masfim - masfim0) / fabs(masfim) > state.input.CriterioConvergPerm ||
            fabs(presteste - presteste0) / fabs(presteste) > state.input.CriterioConvergPerm) &&
           state.steadyIteration < limIter) {

        masfim0 = masfim;
        presteste0 = presteste;
        if (state.steadyIteration == 0 && state.input.lingas > 0 && state.input.nvalvgas > 0 && state.networkCoupled == 1)
            state.updaters.initializeTubingConnectionSteady(); // antes de iniciar a primeira iteracao de marcha,
        // faz-se uma estimativa inicial de como se da o acopamento termico entre a coluna e o anular,
        // caso se tenha linha de gas
        else if (state.input.lingas > 0 && state.input.nvalvgas > 0 && state.networkCoupled == 1)
            state.updaters.connectTubingSteady(); // acoplamento termico feito com valores de pressao e temperatura
        // obtidas na primeira iteracao de marcha
        int i;
        guessNeedsCorrection = 1;
        while (guessNeedsCorrection == 1) { // opcao antiga, ja nao tem mais efeito
            // efetivamente, este while sempre so e feito uma vez, quando a marcha consegue ir ate
            // a ultima celula sem problemas, caso ocorra algum problema, a marcha e finalizada e
            // sai do metodo retornando ou 1e10 ou -1e10

            // inicializando as pressoes e fracoes volumetricas das celulas iniciais, centro de celula
            // e fronteira de celula
            state.cells[0].presauxL = pchute;
            state.cells[0].presLini = pchute;
            state.cells[0].presL = pchute;
            state.cells[0].pres = pchute;
            state.cells[1].presL = pchute;
            state.cells[0].presini = pchute;
            state.cells[1].presLini = pchute;
            state.cells[0].presaux = pchute;
            state.cells[1].presauxL = pchute;

            state.cells[0].alf = alfini;
            state.cells[0].alfini = alfini;
            state.cells[0].bet = betini;
            state.cells[0].betini = betini;
            state.cells[1].alfL = state.cells[0].alf;
            state.cells[1].alfLini = state.cells[0].alf;
            state.cells[0].alfPigD = state.cells[0].alf;
            state.cells[0].alfPigDini = state.cells[0].alf;
            state.cells[0].alfPigE = state.cells[0].alf;
            state.cells[0].alfPigEini = state.cells[0].alf;
            state.cells[1].betL = state.cells[0].bet;
            state.cells[1].betLini = state.cells[0].bet;
            state.cells[0].betPigD = state.cells[0].bet;
            state.cells[0].betPigDini = state.cells[0].bet;
            state.cells[0].betPigE = state.cells[0].bet;
            state.cells[0].betPigEini = state.cells[0].bet;
            state.cells[0].betI = state.cells[0].bet;
            state.cells[1].betLI = state.cells[0].bet;
            // verifica se ja existe algum problema no inicio da marcha
            if (state.cells[0].pres <= 0.1 ||
                (state.input.usaTabela == 1 && (state.input.tabent.pmin - state.cells[0].pres) > (*state.globals).localtiny))
                return -1e10;
            else if ((state.cells[0].acsr.tipo == 3 &&
                      (state.cells[0].acsr.ipr.Pres - state.cells[0].pres) < (*state.globals).localtiny))
                return 1e10;
            // IniciaVazValvGasPerm e um metodo que faz uma estimativa da vazao na valvula de GL
            // quando ainda nao foi feita a marcha na linha de gas. Neste caso, ele recebe o
            // indice da celula de producao e verifica se nesta celula existe uma VGL, se existir,
            // caso a condicao na linha de gas seja vazao injetada, divide a vazao injetada pelo numero de
            // valvulas e indica este valor para a VGL relacionada a celula de producao
            // caso a condicao seja pressao de injecao, faz-se uma estimativa da pressao na linha de gas
            // na posicao da VGL por hidrotatica e com isto se calcula a vazao de injecao da VGL
            if (state.input.lingas > 0 && state.input.nvalvgas > 0 && state.steadyIteration == 0 && state.convergenceMonitor > 0.1) {
                withGasInletCondition(state.gasCells[0].tipoCC,
                                      [&](auto condition) { estimateValveGasPressure(state, condition); });
                state.updaters.initializeSteadyValveGasFlowRate(0);
            }
            i = 1;
            // inicio da marcha propriamente dita
            double abortValue;
            if (advanceProductionCells(state, pchute, i, abortValue))
                return abortValue;
            if (i == state.lastCell + 1)
                guessNeedsCorrection = 0; // fim da marcha
        }
        // apÃ³s o fim da marcha da linha de produÃ§Ã£o, Ã© feita a marcha da linha de gas
        // caso exista
        marchGasLineAndCoupleAnnulus(state, pchute);

        masfim = state.cells[state.lastCell - 1].MC; // guarda valor de vazao para se calcular o erro, quando a
        // opcao de acelerador de convergencia esta desligado
        presteste = state.cells[state.lastCell].pres; // guarda valor de pressao para se calcular o erro, quando a
        // opcao de acelerador de convergencia esta desligado
        state.steadyIteration++; // atualiza a ieteracao da marcha
        if (state.steadyIteration > 200 && state.input.AP == 0)
            NumError("ConvergÃƒÂªncia em marchaProdPerm1 atingiu maximo de iteracoes");
        else if (state.steadyIteration > 200)
            return 1.1e10;
    }

    double corrigePresF = 0.;
    if (((*state.globals).chaverede == 0 || state.endNode == 1 || (*state.globals).chaveRedeParalela == 1) && state.input.corrigeContSep == 1)
        corrigePresF = steadyPressureAtLastCell(state);

    state.baseConvergenceMonitor = state.gasSurfacePressure;
    return state.gasSurfacePressure - (state.cells[state.lastCell].pres + corrigePresF); // caso a marcha tenha conseguido ir atÃ© a Ãºltima celula,
    // retorna a diferenca entre a pressao a montante do choke e a pressao da ultima celula
    // calculada pela marcha
}

double marchReverseProductionSteady(const SteadyStateState &state, double pchute) {

    int guessNeedsCorrection = 1;
    double alfini = 0.;
    double betini = 0.;
    state.slowHeatTransferThreshold = 0.1;

    seedFirstCellVoidFraction(state, pchute, alfini, betini, DryGasFlashTarget::sourceFluid);
    // esta marcha e feita para quando se tem alguma fonte no inicio da tubulacao,
    // portanto, admite-se que o duto esta fechado e coloca-se uma fonte no centro da
    // primeira celula. As vazoes na fronteira esquerda da celula sÃ£o portanto = 0
    state.cells[0].tempL = state.cells[0].temp;
    state.cells[1].tempL = state.cells[0].temp;
    state.cells[0].tempini = state.cells[0].temp;

    state.cells[0].ML = 0.;
    state.cells[0].MC = 0.;
    state.cells[1].ML = 0.;
    state.cells[0].MliqiniL = 0.;
    state.cells[0].Mliqini = 0.;
    state.cells[1].MliqiniL = 0.;
    state.cells[0].MComp = 0.;
    state.cells[0].QLL = 0.;
    state.cells[0].QL = 0;
    state.cells[1].QLL = 0.;
    state.cells[0].QG = 0.;
    double masfim = 1.;
    double masfim0 = -10;
    double presteste = 1.;
    double presteste0 = -10;
    double tempteste = state.cells[0].temp;
    double tempteste0 = -1000;
    state.steadyIteration = 0;

    double pchute0 = pchute;
    int limIter = 2; // limite de iteracoes quando a opcao de aceleracao da convergencia esta desligado. desaconselhavel
    // desligar esta opcao, os ganhos sao pouco e a convergencia se torna instavel, principalmente
    // quando o tramo faz perte de um sistema de redes
    if (state.input.AceleraConvergPerm == 1) { // opcao aceleracao de convergencia ligada
        limIter = 1;                   // em geral faz-se apenas duas iteracoes de marcha para um determinado chute
        // Here the original went on to test for a pressure condition on the gas
        // line and set limIter to 1 again: the same value, so the test selected
        // nothing, and SC-012 removed it (evidencia/anomalias.md, A9-02). The
        // note that went with it, kept below, describes an extra march that
        // only searchProductionBottomHolePressureTertiary makes.
        //
        // no caso de se ter
        // uma condicao de contorno na injecao de gas = pressao, observou-se que o acoplamento dinamico
        // entre a linha de gas e de producao e mais difoicil, para se conseguir um sistema
        // melhor acoplado, deve-se fazer uma marcha iterativa a mais
    }
    // arq.CriterioConvergPerm Ã© um criterio de convergencia da marcha, so faz sentido
    // quando a aceleracao de convergencia esta desligada. Observe que esta convergencia nÃ£o
    // e de fato a conevregencia do problema, e apenas um repeticao de marcha para um determinado
    // chute de pressao ou de vazao no inicio da tubulacao. O que a convergencia busca de fato e
    // determinar qual a pressao ou vazao de fundo que satisfaz as condicoes de contorno no fim
    // da tubulacao, esta busca e feita nos metodos de busca.
    while ((fabs(masfim - masfim0) / fabs(masfim) > state.input.CriterioConvergPerm ||
            fabs(presteste - presteste0) / fabs(presteste) > state.input.CriterioConvergPerm) &&
           (state.steadyIteration < limIter || fabs(tempteste - tempteste0) / ((tempteste) + 273) > 0.001)) {

        masfim0 = masfim;
        presteste0 = presteste;
        tempteste0 = tempteste;
        if (state.steadyIteration == 0 && state.input.lingas > 0 && state.input.nvalvgas > 0 && state.networkCoupled == 1)
            state.updaters.initializeTubingConnectionSteady(); // antes de iniciar a primeira iteracao de marcha,
        // faz-se uma estimativa inicial de como se da o acopamento termico entre a coluna e o anular,
        // caso se tenha linha de gas
        else if (state.input.lingas > 0 && state.input.nvalvgas > 0 && state.networkCoupled == 1)
            state.updaters.connectTubingSteady(); // acoplamento termico feito com valores de pressao e temperatura
        // obtidas na primeira iteracao de marcha
        int i;
        guessNeedsCorrection = 1;
        while (guessNeedsCorrection == 1) { // opcao antiga, ja nao tem mais efeito
            // efetivamente, este while sempre so e feito uma vez, quando a marcha consegue ir ate
            // a ultima celula sem problemas, caso ocorra algum problema, a marcha e finalizada e
            // sai do metodo retornando ou 1e10 ou -1e10

            // inicializando as pressoes e fracoes volumetricas das celulas iniciais, centro de celula
            // e fronteira de celula
            state.cells[0].presauxL = pchute;
            state.cells[0].presLini = pchute;
            state.cells[0].presL = pchute;
            state.cells[0].pres = pchute;
            state.cells[1].presL = pchute;
            state.cells[0].presini = pchute;
            state.cells[1].presLini = pchute;
            state.cells[0].presaux = pchute;
            state.cells[1].presauxL = pchute;

            state.cells[0].alf = alfini;
            state.cells[0].alfini = alfini;
            state.cells[0].bet = betini;
            state.cells[0].betini = betini;
            state.cells[1].alfL = state.cells[0].alf;
            state.cells[1].alfLini = state.cells[0].alf;
            state.cells[0].alfPigD = state.cells[0].alf;
            state.cells[0].alfPigDini = state.cells[0].alf;
            state.cells[0].alfPigE = state.cells[0].alf;
            state.cells[0].alfPigEini = state.cells[0].alf;
            state.cells[1].betL = state.cells[0].bet;
            state.cells[1].betLini = state.cells[0].bet;
            state.cells[0].betPigD = state.cells[0].bet;
            state.cells[0].betPigDini = state.cells[0].bet;
            state.cells[0].betPigE = state.cells[0].bet;
            state.cells[0].betPigEini = state.cells[0].bet;
            state.cells[0].betI = state.cells[0].bet;
            state.cells[1].betLI = state.cells[0].bet;
            // verifica se ja existe algum problema no inicio da marcha
            if (state.cells[0].pres <= 0.1 ||
                (state.input.usaTabela == 1 && (state.input.tabent.pmin - state.cells[0].pres) > (*state.globals).localtiny))
                return -1e10;

            if (state.input.lingas > 0 && state.input.nvalvgas > 0 && state.steadyIteration == 0)
                state.updaters.initializeSteadyValveGasFlowRate(0);
            i = 1;
            // inicio da marcha propriamente dita
            double abortValue;
            if (advanceReverseProductionCells(state, pchute, i, abortValue))
                return abortValue;
            if (i == state.lastCell + 1)
                guessNeedsCorrection = 0; // fim da marcha
        }

        masfim = state.cells[state.lastCell - 1].MC; // guarda valor de vazao para se calcular o erro, quando a
        // opcao de acelerador de convergencia esta desligado
        presteste = state.cells[state.lastCell].pres; // guarda valor de pressao para se calcular o erro, quando a
        // opcao de acelerador de convergencia esta desligado
        state.steadyIteration++; // atualiza a ieteracao da marcha
        if (state.steadyIteration > 200 && state.input.AP == 0)
            NumError("ConvergÃƒÂªncia em marchaProdPerm1 atingiu maximo de iteracoes");
        else if (state.steadyIteration > 200)
            return 1.1e10;
        if (state.input.tipoFluido != 10000) {
            state.cells[state.lastCell].temp = state.casingTemperature;
            if (state.steadyIteration < 100) {
                int lento = 0;
                double media = 0.;
                double desvio = 0;
                for (int ktemp = state.lastCell - 1; ktemp >= 1; ktemp--) {

                    double area = state.cells[ktemp].duto.area;
                    double ugsmed;
                    ugsmed = fabs(state.cells[ktemp].QG) / area; // velocidade superficial de gas
                    double ulsmed;
                    ulsmed = fabs(state.cells[ktemp].QL) / area; // velocidade superficial de liquido
                    if (fabs(ugsmed + ulsmed) <= state.slowHeatTransferThreshold)
                        lento += 1;
                    media += fabs(ugsmed + ulsmed);
                }
                media /= (state.lastCell - 1);
                desvio = (media - 0.1);
                if (lento > 0 && lento < state.lastCell) {
                    int para;
                    para = 0;
                    state.slowHeatTransferThreshold = 100.;
                }
                for (int ktemp = state.lastCell - 1; ktemp >= 0; ktemp--) {
                    state.updaters.advanceReverseSteadyTemperature(ktemp, 0); // faz o avanco da temperatura, da celula i-1 para a celula i
                    // verifica se teve algum problema nos limites de temperatura
                    // caso se esteja trabalhando com tabela PVTSim
                    if (state.input.usaTabela == 1 && (state.cells[ktemp].temp - state.input.tabent.tmin) < (*state.globals).localtiny)
                        state.cells[ktemp].temp = state.input.tabent.tmin;
                    state.updaters.updateProductionTemperaturePeriphery(ktemp); // mera atualizacao de atributos de temperatura a esquerda e a direita
                }
            }
            tempteste = state.cells[0].temp; // 0.5*(celula[0].temp+tempteste);
        } else {
            refreshProperties(state);
            refreshSteadyThermalVelocities(state);
            computePseudoTransientTimeStep(state);
            for (int kontaPseudo = 0; kontaPseudo < 20; kontaPseudo++) {
                state.cells[0].tempini = state.cells[0].temp;
                state.cells[1].tempLini = state.cells[1].tempL;
                state.cells[1].tempL = state.cells[0].temp;
                state.cells[state.lastCell].temp = state.casingTemperature;
                for (int i = 0; i < state.lastCell; i++) {
                    state.cells[i].tempini = state.cells[i].temp;
                }
                for (int itemp = 1; itemp <= state.lastCell; itemp++) {
                    state.updaters.computeTemperature(itemp, state.cells[itemp].tempini, 1);
                }
                refreshProperties(state);
                computePseudoTransientTimeStep(state);
            }
        }
    }

    double corrigePresF = 0.;
    if (((*state.globals).chaverede == 0 || state.endNode == 1 || (*state.globals).chaveRedeParalela == 1) && state.input.corrigeContSep == 1)
        corrigePresF = steadyPressureAtLastCell(state);

    state.baseConvergenceMonitor = state.gasSurfacePressure;
    return state.gasSurfacePressure - (state.cells[state.lastCell].pres + corrigePresF); // caso a marcha tenha conseguido ir atÃ© a Ãºltima celula,
    // retorna a diferenca entre a pressao a montante do choke e a pressao da ultima celula
    // calculada pela marcha
}

double marchProductionSteadySecondary(const SteadyStateState &state, double pchute) {

    int guessNeedsCorrection = 1;

    double alfini = 0.;
    double betini = 0.;

    seedFirstCellVoidFraction(state, pchute, alfini, betini, DryGasFlashTarget::cellFluid);

    // esta marcha e feita para quando se tem alguma fonte no inicio da tubulacao,
    // portanto, admite-se que o duto esta fechado e coloca-se uma fonte no centro da
    // primeira celula. As vazoes na fronteira esquerda da celula sÃ£o portanto = 0
    state.cells[0].tempL = state.cells[0].temp;
    state.cells[1].tempL = state.cells[0].temp;
    state.cells[0].tempini = state.cells[0].temp;

    state.cells[0].ML = 0.;
    state.cells[0].MC = 0.;
    state.cells[1].ML = 0.;
    state.cells[0].MliqiniL = 0.;
    state.cells[0].Mliqini = 0.;
    state.cells[1].MliqiniL = 0.;
    state.cells[0].MComp = 0.;
    state.cells[0].QLL = 0.;
    state.cells[0].QL = 0;
    state.cells[1].QLL = 0.;
    state.cells[0].QG = 0.;
    double masfim = 1.;
    double masfim0 = -10;
    double presteste = 1.;
    double presteste0 = -10;
    state.steadyIteration = 0;

    int limIter = 2; // limite de iteracoes quando a opcao de aceleracao da convergencia esta desligado. desaconselhavel
    // desligar esta opcao, os ganhos sao pouco e a convergencia se torna instavel, principalmente
    // quando o tramo faz perte de um sistema de redes
    if (state.input.AceleraConvergPerm == 1) { // opcao aceleracao de convergencia ligada
        limIter = 1;                   // em geral faz-se apenas duas iteracoes de marcha para um determinado chute
        // Here the original went on to test for a pressure condition on the gas
        // line and set limIter to 1 again: the same value, so the test selected
        // nothing, and SC-012 removed it (evidencia/anomalias.md, A9-02). The
        // note that went with it, kept below, describes an extra march that
        // only searchProductionBottomHolePressureTertiary makes.
        //
        // no caso de se ter
        // uma condicao de contorno na injecao de gas = pressao, observou-se que o acoplamento dinamico
        // entre a linha de gas e de producao e mais difoicil, para se conseguir um sistema
        // melhor acoplado, deve-se fazer uma marcha iterativa a mais
    }
    // arq.CriterioConvergPerm Ã© um criterio de convergencia da marcha, so faz sentido
    // quando a aceleracao de convergencia esta desligada. Observe que esta convergencia nÃ£o
    // e de fato a conevregencia do problema, e apenas um repeticao de marcha para um determinado
    // chute de pressao ou de vazao no inicio da tubulacao. O que a convergencia busca de fato e
    // determinar qual a pressao ou vazao de fundo que satisfaz as condicoes de contorno no fim
    // da tubulacao, esta busca e feita nos metodos de busca.
    while ((fabs(masfim - masfim0) / masfim > state.input.CriterioConvergPerm ||
            fabs(presteste - presteste0) / fabs(presteste) > state.input.CriterioConvergPerm) &&
           (state.steadyIteration < limIter)) {
        masfim0 = masfim;
        presteste0 = presteste;
        if (state.steadyIteration == 0 && state.input.lingas > 0 && state.input.nvalvgas > 0 && state.networkCoupled == 1)
            state.updaters.initializeTubingConnectionSteady(); // antes de iniciar a primeira iteracao de marcha,
        // faz-se uma estimativa inicial de como se da o acopamento termico entre a coluna e o anular,
        // caso se tenha linha de gas
        else if (state.input.lingas > 0 && state.input.nvalvgas > 0 && state.networkCoupled == 1)
            state.updaters.connectTubingSteady(); // acoplamento termico feito com valores de pressao e temperatura
        // obtidas na primeira iteracao de marcha
        int i;
        guessNeedsCorrection = 1;
        while (guessNeedsCorrection == 1) { // opcao antiga, ja nao tem mais efeito
            // efetivamente, este while sempre so e feito uma vez, quando a marcha consegue ir ate
            // a ultima celula sem problemas, caso ocorra algum problema, a marcha e finalizada e
            // sai do metodo retornando ou 1e10 ou -1e10

            // inicializando as pressoes e fracoes volumetricas das celulas iniciais, centro de celula
            // e fronteira de celula
            state.cells[0].presauxL = pchute;
            state.cells[0].presLini = pchute;
            state.cells[0].presL = pchute;
            state.cells[0].pres = pchute;
            state.cells[1].presL = pchute;
            state.cells[0].presini = pchute;
            state.cells[1].presLini = pchute;
            state.cells[0].presaux = pchute;
            state.cells[1].presauxL = pchute;

            state.cells[0].alf = alfini;
            state.cells[0].alfini = alfini;
            state.cells[0].bet = betini;
            state.cells[0].betini = betini;
            state.cells[1].alfL = state.cells[0].alf;
            state.cells[1].alfLini = state.cells[0].alf;
            state.cells[0].alfPigD = state.cells[0].alf;
            state.cells[0].alfPigDini = state.cells[0].alf;
            state.cells[0].alfPigE = state.cells[0].alf;
            state.cells[0].alfPigEini = state.cells[0].alf;
            state.cells[1].betL = state.cells[0].bet;
            state.cells[1].betLini = state.cells[0].bet;
            state.cells[0].betPigD = state.cells[0].bet;
            state.cells[0].betPigDini = state.cells[0].bet;
            state.cells[0].betPigE = state.cells[0].bet;
            state.cells[0].betPigEini = state.cells[0].bet;
            state.cells[0].betI = state.cells[0].bet;
            state.cells[1].betLI = state.cells[0].bet;

            // verifica se ja existe algum problema no inicio da marcha
            if (state.cells[0].pres <= 0.1 ||
                (state.input.usaTabela == 1 && (state.input.tabent.pmin - state.cells[0].presaux) > (*state.globals).localtiny))
                return -1e10;
            else if ((state.cells[0].acsr.tipo == 3 && (state.cells[0].acsr.ipr.Pres - state.cells[0].pres) < 1e-15)) {
                return 1e10;
            } else if ((state.cells[0].acsr.tipo == 15 && (state.cells[0].acsr.radialPoro.pRes[0] - state.cells[0].pres) < 1e-15)) {
                return 1e10;
            } else if ((state.cells[0].acsr.tipo == 16 && (state.cells[0].acsr.poroso2D.dados.pRes - state.cells[0].pres) < 1e-15)) {
                return 1e10;
            }
            // IniciaVazValvGasPerm e um metodo que faz uma estimativa da vazao na valvula de GL
            // quando ainda nao foi feita a marcha na linha de gas. Neste caso, ele recebe o
            // indice da celula de producao e verifica se nesta celula existe uma VGL, se existir,
            // caso a condicao na linha de gas seja vazao injetada, divide a vazao injetada pelo numero de
            // valvulas e indica este valor para a VGL relacionada a celula de producao
            // caso a condicao seja pressao de injecao, faz-se uma estimativa da pressao na linha de gas
            // na posicao da VGL por hidrotatica e com isto se calcula a vazao de injecao da VGL
            if (state.input.lingas > 0 && state.input.nvalvgas > 0 && state.steadyIteration == 0 && state.convergenceMonitor > 0.1) {
                withGasInletCondition(state.gasCells[0].tipoCC,
                                      [&](auto condition) { estimateValveGasPressure(state, condition); });
                state.updaters.initializeSteadyValveGasFlowRate(0);
            }
            i = 1;
            // inicio da marcha propriamente dita
            double abortValue;
            if (advanceProductionCellsSecondary(state, pchute, i, abortValue))
                return abortValue;
            if (i == state.lastCell + 1)
                guessNeedsCorrection = 0; // fim da marcha
        }
        // apÃ³s o fim da marcha da linha de produÃ§Ã£o, Ã© feita a marcha da linha de gas
        // caso exista
        marchGasLineAndCoupleAnnulus(state, pchute);

        masfim = state.cells[state.lastCell - 1].MC; // guarda valor de vazao para se calcular o erro, quando a
        // opcao de acelerador de convergencia esta desligado
        presteste = state.cells[state.lastCell].pres; // guarda valor de pressao para se calcular o erro, quando a
        // opcao de acelerador de convergencia esta desligado
        state.steadyIteration++; // atualiza a ieteracao da marcha
        if (state.steadyIteration > 200 && state.input.AP == 0)
            NumError("Convergencia em marchaProdPerm2 atingiu maximo de iteracoes");
        else if (state.steadyIteration > 200)
            return 1e10;
    }

    // nesta marcha, a condicao no fim da tubulacao nÃ£o e a pressao a montante do choke,
    // mas a pressao a jusante do choke. Neste caso, a condicao que se deseja convergir
    // e a vazao que passa pelo choke, definida a aprtir da diferenca entre a pressao
    // na ultima celula e a pressao a jusante do choke. Primeiro, portanto, deve-se
    // calcular a vazao que passa pelo choke e compara-la com a vazao massica total
    // na ultima celula
    double maxSup = surfaceChokeMassFlowRate(state);
    if (fabs(masfim) > 1e-15)
        state.baseConvergenceMonitor = fabs(masfim);
    else if (fabs(maxSup) > 1e-15)
        state.baseConvergenceMonitor = fabs(maxSup);
    else
        state.baseConvergenceMonitor = 1.;
    return (masfim - maxSup); // diferenca entre a vazao total na fronteira a esquerda da
    // penultima celula e a vazao que passa pelo choke. Se o chute de pressao for alto
    // maxSup>masfim, retorna valor negativo,
    // se for uma estimativa baixa de pressao de fundo, maxSup<masfim retirna valor positivo
}

namespace {

void seedFirstCellFromFlowRateGuess(const SteadyStateState &state, double mchute, double &alfini, double &betini) {
    if (state.cells[0].acsr.tipo == 1) {
        state.cells[0].temp = state.cells[0].acsr.injg.temp;
        if (state.input.flashCompleto == 2) {
            if (state.input.tabelaDinamica == 0)
                state.cells[0].flui.atualizaPropComp(state.cells[0].pres, state.cells[0].temp);
            state.cells[0].acsr.injg.FluidoPro.atualizaPropComp(state.cells[0].pres, state.cells[0].temp, -1, NULL, NULL, state.cells[0].acsr.injg.seco);
        }
        state.cells[0].acsr.injg.QGas = mchute;
    } else if (state.cells[0].acsr.tipo == 2) {
        state.cells[0].temp = state.cells[0].acsr.injl.temp;
        if (state.input.flashCompleto == 2) {
            if (state.input.tabelaDinamica == 0)
                state.cells[0].flui.atualizaPropComp(state.cells[0].pres, state.cells[0].temp);
            state.cells[0].acsr.injl.FluidoPro.atualizaPropComp(state.cells[0].pres, state.cells[0].temp);
        }
        state.cells[0].acsr.injl.QLiq = mchute;
    }

    if (state.cells[0].acsr.tipo == 1) {
        if (state.cells[0].acsr.injg.seco == 1) {
            alfini = 1.;
            betini = 0.;
        } else {
            double masgas = state.cells[0].acsr.injg.VMas(state.cells[0].pres, state.cells[0].temp);
            double quality;
            if (state.input.flashCompleto != 2)
                quality = state.cells[0].acsr.injg.FluidoPro.FracMassHidra(1., 20.);
            else
                quality = state.cells[0].acsr.injg.FluidoPro.dStockTankVaporMassFraction;
            double masT = masgas / quality;
            quality = state.cells[0].acsr.injg.FluidoPro.FracMassHidra(state.cells[0].pres, state.cells[0].temp);
            double qgas = masT * quality /
                          state.cells[0].acsr.injg.FluidoPro.MasEspGas(state.cells[0].pres, state.cells[0].temp);
            double qliq = masT * (1. - quality) /
                          state.cells[0].acsr.injg.FluidoPro.MasEspLiq(state.cells[0].pres, state.cells[0].temp);
            double qcomp = state.cells[0].acsr.injg.razCompGas *
                           state.cells[0].acsr.injg.QGas * state.cells[0].acsr.injg.fluidocol.MasEspFlu(1., 20.) /
                           state.cells[0].acsr.injg.fluidocol.MasEspFlu(state.cells[0].pres, state.cells[0].temp);
            qcomp /= 86400.;
            alfini = qgas / (qliq + qcomp + qgas);
            if ((fabs(qcomp) + fabs(qliq)) > 1e-15)
                betini = fabs(qcomp) / (fabs(qcomp) + fabs(qliq));
            else
                betini = 0.;
        }
    } else if (state.cells[0].acsr.tipo == 2) {
        double qgas = state.cells[0].acsr.injl.QLiq * (1 - state.cells[0].acsr.injl.bet) *
                      (1. - state.cells[0].acsr.injl.FluidoPro.BSW) *
                      (state.cells[0].acsr.injl.FluidoPro.RGO -
                       state.cells[0].acsr.injl.FluidoPro.rDgD * state.cells[0].acsr.injl.FluidoPro.RS(state.cells[0].pres, state.cells[0].temp) * 6.29 / 35.31467) *
                      state.cells[0].acsr.injl.FluidoPro.Deng * 1.225 / state.cells[0].acsr.injl.FluidoPro.MasEspGas(state.cells[0].pres, state.cells[0].temp);
        double qliq = state.cells[0].acsr.injl.QLiq * (1 - state.cells[0].acsr.injl.bet) *
                          (1. - state.cells[0].acsr.injl.FluidoPro.BSW) * state.cells[0].acsr.injl.FluidoPro.BOFunc(state.cells[0].pres, state.cells[0].temp) +
                      state.cells[0].acsr.injl.QLiq * (1 - state.cells[0].acsr.injl.bet) *
                          state.cells[0].acsr.injl.FluidoPro.BSW * state.cells[0].acsr.injl.FluidoPro.BAFunc(state.cells[0].pres, state.cells[0].temp) +
                      state.cells[0].acsr.injl.QLiq * state.cells[0].acsr.injl.bet;
        alfini = qgas / (qliq + qgas);
        betini = state.cells[0].acsr.injl.bet;
    }
}

void advanceProductionCellsPressureToPressureSecondary(const SteadyStateState &state, int &i) {
    while (i <= state.lastCell && state.cells[i - 1].pres >= 1.) {

        advanceUpstreamSteadyPressure(state, i, 0);
        refreshUpstreamProductionPeriphery(state, i);
        if (state.input.flashCompleto != 2)
            advanceSteadyMass(state, i);
        else
            advanceCompositionalSteadyMass(state, i);
        state.updaters.advanceSteadyTemperature(i, 0);
        if (isnan(state.cells[i].temp))
            NumError("Temperatrura na linha de producao com valor NaN");
        if (state.input.usaTabela == 1 && (state.cells[i].temp - state.input.tabent.tmin) < (*state.globals).localtiny)
            state.cells[i].temp = state.input.tabent.tmin;
        state.updaters.updateProductionTemperaturePeriphery(i);
        advanceDownstreamSteadyPressure(state, i, 0);
        refreshDownstreamProductionPeriphery(state, i);
        if (state.input.tipoFluido == 0)
            advanceSteadyMassTransfer(state, i - 1);
        else
            advanceSteadyGasMassTransfer(state, i - 1);
        if (state.input.ordperm > 1) {
            double D0presaux = state.cells[i].presaux - state.cells[i - 1].pres;
            double D0pres = state.cells[i].pres - state.cells[i].presaux;
            double D0temp = state.cells[i].temp - state.cells[i - 1].temp;
            advanceUpstreamSteadyPressure(state, i, 1);
            refreshUpstreamProductionPeriphery(state, i);
            advanceSteadyMass(state, i);
            state.updaters.advanceSteadyTemperature(i, 1);
            advanceDownstreamSteadyPressure(state, i, 1);
            state.cells[i].pres = 0.5 * (state.cells[i].presaux + D0pres + state.cells[i].pres);
            state.cells[i].presaux = 0.5 * (state.cells[i - 1].pres + D0presaux + state.cells[i].presaux);
            state.cells[i].temp = 0.5 * (state.cells[i - 1].temp + D0temp + state.cells[i].temp);
            refreshDownstreamProductionPeriphery(state, i);
            refreshUpstreamProductionPeriphery(state, i);
            state.updaters.updateProductionTemperaturePeriphery(i);
            advanceSteadyMass(state, i);
            if (state.input.tipoFluido == 0)
                advanceSteadyMassTransfer(state, i - 1);
            else
                advanceSteadyGasMassTransfer(state, i - 1);
        }

        i++;
    }
}
}  // namespace


double marchProductionPressureToPressure(const SteadyStateState &state, double mchute) {

    int guessNeedsCorrection = 1;
    double alfini = 0.;
    double betini = 0.;

    seedFirstCellFromFlowRateGuess(state, mchute, alfini, betini);

    state.cells[0].tempL = state.cells[0].temp;
    state.cells[1].tempL = state.cells[0].temp;
    state.cells[0].tempini = state.cells[0].temp;

    state.cells[0].ML = 0.;
    state.cells[0].MC = 0.;
    state.cells[1].ML = 0.;
    state.cells[0].MliqiniL = 0.;
    state.cells[0].Mliqini = 0.;
    state.cells[1].MliqiniL = 0.;
    state.cells[0].MComp = 0.;
    state.cells[0].QLL = 0.;
    state.cells[0].QL = 0;
    state.cells[1].QLL = 0.;
    state.cells[0].QG = 0.;

    int i;

    state.cells[0].presauxL = state.cells[0].pres;
    state.cells[0].presLini = state.cells[0].pres;
    state.cells[0].presL = state.cells[0].pres;
    state.cells[1].presL = state.cells[0].pres;
    state.cells[0].presini = state.cells[0].pres;
    state.cells[1].presLini = state.cells[0].pres;
    state.cells[0].presaux = state.cells[0].pres;
    state.cells[1].presauxL = state.cells[0].pres;

    state.cells[0].alf = alfini;
    state.cells[0].alfini = alfini;
    state.cells[0].bet = betini;
    state.cells[0].betini = betini;
    state.cells[1].alfL = state.cells[0].alf;
    state.cells[1].alfLini = state.cells[0].alf;
    state.cells[0].alfPigD = state.cells[0].alf;
    state.cells[0].alfPigDini = state.cells[0].alf;
    state.cells[0].alfPigE = state.cells[0].alf;
    state.cells[0].alfPigEini = state.cells[0].alf;
    state.cells[1].betL = state.cells[0].bet;
    state.cells[1].betLini = state.cells[0].bet;
    state.cells[0].betPigD = state.cells[0].bet;
    state.cells[0].betPigDini = state.cells[0].bet;
    state.cells[0].betPigE = state.cells[0].bet;
    state.cells[0].betPigEini = state.cells[0].bet;
    state.cells[0].betI = state.cells[0].bet;
    state.cells[1].betLI = state.cells[0].bet;

    state.steadyIteration = 0;
    while (state.steadyIteration < 3) {
        i = 1;
        while (i <= state.lastCell && state.cells[i - 1].pres >= 1.) {

            advanceUpstreamSteadyPressure(state, i, 0);
            refreshUpstreamProductionPeriphery(state, i);
            if (state.input.flashCompleto != 2)
                advanceSteadyMass(state, i);
            else
                advanceCompositionalSteadyMass(state, i);
            state.updaters.advanceSteadyTemperature(i, 0);
            if (i > 1160) {
                int para;
                para = 0;
            }
            if (state.input.usaTabela == 1 && (state.cells[i].temp - state.input.tabent.tmin) < (*state.globals).localtiny)
                state.cells[i].temp = state.input.tabent.tmin;
            state.updaters.updateProductionTemperaturePeriphery(i);
            advanceDownstreamSteadyPressure(state, i, 0);
            refreshDownstreamProductionPeriphery(state, i);
            if (state.input.tipoFluido == 0)
                advanceSteadyMassTransfer(state, i - 1);
            else
                advanceSteadyGasMassTransfer(state, i - 1);
            if (state.input.ordperm > 1) {
                double D0presaux = state.cells[i].presaux - state.cells[i - 1].pres;
                double D0pres = state.cells[i].pres - state.cells[i].presaux;
                double D0temp = state.cells[i].temp - state.cells[i - 1].temp;
                advanceUpstreamSteadyPressure(state, i, 1);
                refreshUpstreamProductionPeriphery(state, i);
                advanceSteadyMass(state, i);
                state.updaters.advanceSteadyTemperature(i, 1);
                if (isnan(state.cells[i].temp))
                    NumError("Temperatrura na linha de producao com valor NaN");
                advanceDownstreamSteadyPressure(state, i, 1);
                state.cells[i].pres = 0.5 * (state.cells[i].presaux + D0pres + state.cells[i].pres);
                state.cells[i].presaux = 0.5 * (state.cells[i - 1].pres + D0presaux + state.cells[i].presaux);
                state.cells[i].temp = 0.5 * (state.cells[i - 1].temp + D0temp + state.cells[i].temp);
                refreshDownstreamProductionPeriphery(state, i);
                refreshUpstreamProductionPeriphery(state, i);
                state.updaters.updateProductionTemperaturePeriphery(i);
                advanceSteadyMass(state, i);
                if (state.input.tipoFluido == 0)
                    advanceSteadyMassTransfer(state, i - 1);
                else
                    advanceSteadyGasMassTransfer(state, i - 1);
            }
            i++;
            if (state.cells[i - 1].pres <= 1)
                return -1e10;
        }

        state.steadyIteration++;
    }

    double corrigePresF = 0.;
    if ((*state.globals).chaverede == 0 || state.endNode == 1 || (*state.globals).chaveRedeParalela == 1)
        corrigePresF = steadyPressureAtLastCell(state);

    state.baseConvergenceMonitor = state.gasSurfacePressure;
    return state.cells[state.lastCell].pres + corrigePresF - state.gasSurfacePressure;
}

double marchReverseProductionPressureToPressure(const SteadyStateState &state, double mchute) {

    double alfini = 0.;
    double betini = 0.;
    state.slowHeatTransferThreshold = 0.05;

    seedFirstCellFromFlowRateGuess(state, mchute, alfini, betini);

    state.cells[0].tempL = state.cells[0].temp;
    state.cells[1].tempL = state.cells[0].temp;
    state.cells[0].tempini = state.cells[0].temp;

    state.cells[0].ML = 0.;
    state.cells[0].MC = 0.;
    state.cells[1].ML = 0.;
    state.cells[0].MliqiniL = 0.;
    state.cells[0].Mliqini = 0.;
    state.cells[1].MliqiniL = 0.;
    state.cells[0].MComp = 0.;
    state.cells[0].QLL = 0.;
    state.cells[0].QL = 0;
    state.cells[1].QLL = 0.;
    state.cells[0].QG = 0.;

    int i;

    state.cells[0].presauxL = state.cells[0].pres;
    state.cells[0].presLini = state.cells[0].pres;
    state.cells[0].presL = state.cells[0].pres;
    state.cells[1].presL = state.cells[0].pres;
    state.cells[0].presini = state.cells[0].pres;
    state.cells[1].presLini = state.cells[0].pres;
    state.cells[0].presaux = state.cells[0].pres;
    state.cells[1].presauxL = state.cells[0].pres;

    state.cells[0].alf = alfini;
    state.cells[0].alfini = alfini;
    state.cells[0].bet = betini;
    state.cells[0].betini = betini;
    state.cells[1].alfL = state.cells[0].alf;
    state.cells[1].alfLini = state.cells[0].alf;
    state.cells[0].alfPigD = state.cells[0].alf;
    state.cells[0].alfPigDini = state.cells[0].alf;
    state.cells[0].alfPigE = state.cells[0].alf;
    state.cells[0].alfPigEini = state.cells[0].alf;
    state.cells[1].betL = state.cells[0].bet;
    state.cells[1].betLini = state.cells[0].bet;
    state.cells[0].betPigD = state.cells[0].bet;
    state.cells[0].betPigDini = state.cells[0].bet;
    state.cells[0].betPigE = state.cells[0].bet;
    state.cells[0].betPigEini = state.cells[0].bet;
    state.cells[0].betI = state.cells[0].bet;
    state.cells[1].betLI = state.cells[0].bet;

    state.steadyIteration = 0;
    double tempteste = state.cells[0].temp;
    double tempteste0 = -1000;
    while (state.steadyIteration < 3 || fabs(tempteste - tempteste0) / ((tempteste) + 273) > 0.001) {
        i = 1;
        tempteste0 = tempteste;
        while (i <= state.lastCell && state.cells[i - 1].pres >= 1.) {

            advanceUpstreamSteadyPressure(state, i, 0);
            refreshUpstreamProductionPeriphery(state, i);
            if (state.input.flashCompleto != 2)
                advanceSteadyMass(state, i);
            else
                advanceCompositionalSteadyMass(state, i);
            if (isnan(state.cells[i].temp))
                NumError("Temperatrura na linha de producao com valor NaN");
            if (state.input.usaTabela == 1 && (state.cells[i].temp - state.input.tabent.tmin) < (*state.globals).localtiny)
                state.cells[i].temp = state.input.tabent.tmin;
            state.updaters.updateProductionTemperaturePeriphery(i);
            advanceDownstreamSteadyPressure(state, i, 0);
            refreshDownstreamProductionPeriphery(state, i);
            if (state.input.tipoFluido == 0)
                advanceSteadyMassTransfer(state, i - 1);
            else
                advanceSteadyGasMassTransfer(state, i - 1);
            if (state.input.ordperm > 1) {
                double D0presaux = state.cells[i].presaux - state.cells[i - 1].pres;
                double D0pres = state.cells[i].pres - state.cells[i].presaux;
                double D0temp = state.cells[i].temp - state.cells[i - 1].temp;
                advanceUpstreamSteadyPressure(state, i, 1);
                refreshUpstreamProductionPeriphery(state, i);
                advanceSteadyMass(state, i);
                advanceDownstreamSteadyPressure(state, i, 1);
                state.cells[i].pres = 0.5 * (state.cells[i].presaux + D0pres + state.cells[i].pres);
                state.cells[i].presaux = 0.5 * (state.cells[i - 1].pres + D0presaux + state.cells[i].presaux);
                state.cells[i].temp = 0.5 * (state.cells[i - 1].temp + D0temp + state.cells[i].temp);
                refreshDownstreamProductionPeriphery(state, i);
                refreshUpstreamProductionPeriphery(state, i);
                state.updaters.updateProductionTemperaturePeriphery(i);
                advanceSteadyMass(state, i);
                if (state.input.tipoFluido == 0)
                    advanceSteadyMassTransfer(state, i - 1);
                else
                    advanceSteadyGasMassTransfer(state, i - 1);
            }
            i++;
            if (state.cells[i - 1].pres <= 1)
                return -1e10;
        }

        state.steadyIteration++;
        if (state.steadyIteration > 200 && state.input.AP == 0)
            NumError("ConvergÃƒÂªncia em marchaProdPerm1 atingiu maximo de iteracoes");
        else if (state.steadyIteration > 200)
            return 1.1e10;
        if (state.input.tipoFluido != 10000) {
            state.cells[state.lastCell].temp = state.casingTemperature;
            if (state.steadyIteration < 100) {
                int lento = 0;
                state.slowHeatTransferThreshold = 0;
                for (int ktemp = state.lastCell - 1; ktemp >= 0; ktemp--) {

                    double area = state.cells[ktemp].duto.area;
                    double ugsmed;
                    ugsmed = fabs(state.cells[ktemp].QG) / area; // velocidade superficial de gas
                    double ulsmed;
                    ulsmed = fabs(state.cells[ktemp].QL) / area; // velocidade superficial de liquido
                    if (fabs(ugsmed + ulsmed) <= 0.1)
                        lento += 1;
                }
                if (state.steadyIteration > 10 && state.slowHeatTransferThreshold < 0.01) {
                    state.slowHeatTransferThreshold -= 0.01;
                    state.steadyIteration = 0;
                }
                for (int ktemp = state.lastCell - 1; ktemp >= 0; ktemp--) {
                    state.updaters.advanceReverseSteadyTemperature(ktemp, 0);
                    if (state.input.usaTabela == 1 && (state.cells[ktemp].temp - state.input.tabent.tmin) < (*state.globals).localtiny)
                        state.cells[ktemp].temp = state.input.tabent.tmin;
                    state.updaters.updateProductionTemperaturePeriphery(ktemp); // mera atualizacao de atributos de temperatura a esquerda e a direita
                }
            }
            tempteste = state.cells[0].temp;
        } else {
            refreshProperties(state);
            refreshSteadyThermalVelocities(state);
            computePseudoTransientTimeStep(state);
            for (int kontaPseudo = 0; kontaPseudo < 20; kontaPseudo++) {
                state.cells[0].tempini = state.cells[0].temp;
                state.cells[1].tempLini = state.cells[1].tempL;
                state.cells[1].tempL = state.cells[0].temp;
                state.cells[state.lastCell].temp = state.casingTemperature;
                for (int i = 0; i < state.lastCell; i++) {
                    state.cells[i].tempini = state.cells[i].temp;
                }
                for (int itemp = 1; itemp <= state.lastCell; itemp++) {
                    state.updaters.computeTemperature(itemp, state.cells[itemp].tempini, 1);
                }
                refreshProperties(state);
                computePseudoTransientTimeStep(state);
            }
        }
    }
    double corrigePresF = 0.;
    if ((*state.globals).chaverede == 0 || state.endNode == 1 || (*state.globals).chaveRedeParalela == 1)
        corrigePresF = steadyPressureAtLastCell(state);
    state.baseConvergenceMonitor = state.gasSurfacePressure;
    return state.cells[state.lastCell].pres + corrigePresF - state.gasSurfacePressure;
}

double marchProductionPressureToPressureSecondary(const SteadyStateState &state, double mchute) {

    int guessNeedsCorrection = 1;
    double alfini = 0.;
    double betini = 0.;

    if (state.cells[0].acsr.tipo == 1) {
        state.cells[0].temp = state.cells[0].acsr.injg.temp;
        if (state.input.flashCompleto == 2) {
            if (state.input.tabelaDinamica == 0)
                state.cells[0].flui.atualizaPropComp(state.cells[0].pres, state.cells[0].temp);
            state.cells[0].acsr.injg.FluidoPro.atualizaPropComp(state.cells[0].pres, state.cells[0].temp, -1, NULL, NULL, state.cells[0].acsr.injg.seco);
        }
        state.cells[0].acsr.injg.QGas = mchute;
    } else if (state.cells[0].acsr.tipo == 2) {
        state.cells[0].temp = state.cells[0].acsr.injl.temp;
        if (state.input.flashCompleto == 2) {
            if (state.input.tabelaDinamica == 0)
                state.cells[0].flui.atualizaPropComp(state.cells[0].pres, state.cells[0].temp);
            state.cells[0].acsr.injl.FluidoPro.atualizaPropComp(state.cells[0].pres, state.cells[0].temp);
        }
        state.cells[0].acsr.injl.QLiq = mchute;
    }

    if (state.cells[0].acsr.tipo == 1) {
        if (state.cells[0].acsr.injg.seco == 1) {
            alfini = 1.;
            betini = 0.;
        } else {
            double masgas = state.cells[0].acsr.injg.VMas(state.cells[0].pres, state.cells[0].temp);
            double quality;
            double masT;
            if (state.input.ConContEntrada != 1) {
                if (state.input.flashCompleto != 2)
                    quality = state.cells[0].acsr.injg.FluidoPro.FracMassHidra(1., 20.);
                else
                    quality = state.cells[0].acsr.injg.FluidoPro.dStockTankVaporMassFraction;
                masT = masgas / quality;
                quality = state.cells[0].acsr.injg.FluidoPro.FracMassHidra(state.cells[0].pres, state.cells[0].temp);
            } else {
                quality = state.inletQuality;
                masT = masgas / quality;
            }
            double qgas = masT * quality /
                          state.cells[0].acsr.injg.FluidoPro.MasEspGas(state.cells[0].pres, state.cells[0].temp);
            double qliq = masT * (1. - quality) /
                          state.cells[0].acsr.injg.FluidoPro.MasEspLiq(state.cells[0].pres, state.cells[0].temp);
            double qcomp = state.cells[0].acsr.injg.razCompGas *
                           state.cells[0].acsr.injg.QGas * state.cells[0].acsr.injg.fluidocol.MasEspFlu(1., 20.) /
                           state.cells[0].acsr.injg.fluidocol.MasEspFlu(state.cells[0].pres, state.cells[0].temp);
            qcomp /= 86400.;
            if (fabs(qliq + qgas) > 1e-15)
                alfini = qgas / (qliq + qcomp + qgas);
            else
                alfini = state.inletQuality;
            if ((fabs(qcomp) + fabs(qliq)) > 1e-15)
                betini = fabs(qcomp) / (fabs(qcomp) + fabs(qliq));
            else
                betini = 0.;
        }
    } else if (state.cells[0].acsr.tipo == 2) {
        double qgas = state.cells[0].acsr.injl.QLiq * (1 - state.cells[0].acsr.injl.bet) *
                      (1. - state.cells[0].acsr.injl.FluidoPro.BSW) *
                      (state.cells[0].acsr.injl.FluidoPro.RGO -
                       state.cells[0].acsr.injl.FluidoPro.rDgD * state.cells[0].acsr.injl.FluidoPro.RS(state.cells[0].pres, state.cells[0].temp) * 6.29 / 35.31467) *
                      state.cells[0].acsr.injl.FluidoPro.Deng * 1.225 / state.cells[0].acsr.injl.FluidoPro.MasEspGas(state.cells[0].pres, state.cells[0].temp);
        double qliq = state.cells[0].acsr.injl.QLiq * (1 - state.cells[0].acsr.injl.bet) *
                          (1. - state.cells[0].acsr.injl.FluidoPro.BSW) * state.cells[0].acsr.injl.FluidoPro.BOFunc(state.cells[0].pres, state.cells[0].temp) +
                      state.cells[0].acsr.injl.QLiq * (1 - state.cells[0].acsr.injl.bet) *
                          state.cells[0].acsr.injl.FluidoPro.BSW * state.cells[0].acsr.injl.FluidoPro.BAFunc(state.cells[0].pres, state.cells[0].temp) +
                      state.cells[0].acsr.injl.QLiq * state.cells[0].acsr.injl.bet;
        if (fabs(qliq + qgas) > 1e-15)
            alfini = fabs(qgas / (qliq + qgas));
        else
            alfini = state.inletQuality;
        betini = state.cells[0].acsr.injl.bet;
    }

    state.cells[0].tempL = state.cells[0].temp;
    state.cells[1].tempL = state.cells[0].temp;
    state.cells[0].tempini = state.cells[0].temp;

    state.cells[0].ML = 0.;
    state.cells[0].MC = 0.;
    state.cells[1].ML = 0.;
    state.cells[0].MliqiniL = 0.;
    state.cells[0].Mliqini = 0.;
    state.cells[1].MliqiniL = 0.;
    state.cells[0].MComp = 0.;
    state.cells[0].QLL = 0.;
    state.cells[0].QL = 0;
    state.cells[1].QLL = 0.;
    state.cells[0].QG = 0.;

    int i;

    state.cells[0].presauxL = state.cells[0].pres;
    state.cells[0].presLini = state.cells[0].pres;
    state.cells[0].presL = state.cells[0].pres;
    state.cells[1].presL = state.cells[0].pres;
    state.cells[0].presini = state.cells[0].pres;
    state.cells[1].presLini = state.cells[0].pres;
    state.cells[0].presaux = state.cells[0].pres;
    state.cells[1].presauxL = state.cells[0].pres;

    state.cells[0].alf = alfini;
    state.cells[0].alfini = alfini;
    state.cells[0].bet = betini;
    state.cells[0].betini = betini;
    state.cells[1].alfL = state.cells[0].alf;
    state.cells[1].alfLini = state.cells[0].alf;
    state.cells[0].alfPigD = state.cells[0].alf;
    state.cells[0].alfPigDini = state.cells[0].alf;
    state.cells[0].alfPigE = state.cells[0].alf;
    state.cells[0].alfPigEini = state.cells[0].alf;
    state.cells[1].betL = state.cells[0].bet;
    state.cells[1].betLini = state.cells[0].bet;
    state.cells[0].betPigD = state.cells[0].bet;
    state.cells[0].betPigDini = state.cells[0].bet;
    state.cells[0].betPigE = state.cells[0].bet;
    state.cells[0].betPigEini = state.cells[0].bet;
    state.cells[0].betI = state.cells[0].bet;
    state.cells[1].betLI = state.cells[0].bet;

    state.steadyIteration = 0;
    while (state.steadyIteration < 3) {
        i = 1;
        advanceProductionCellsPressureToPressureSecondary(state, i);
        state.steadyIteration++;
    }

    double maxSup;
    if (state.cells[state.lastCell].pres >= state.gasSurfacePressure) {
        double tESup = state.cells[state.lastCell].temp;
        double alfSup = state.cells[state.lastCell].alf;
        double betSup = state.cells[state.lastCell].bet;

        double masentrada = state.cells[state.lastCell - 1].MR;
        double massgas = state.cells[state.lastCell - 1].MR - state.cells[state.lastCell - 1].MliqiniR;
        maxSup = 0.;

        double quality;
        state.outletPressure = state.cells[state.lastCell].pres;
        double rholp = state.cells[state.lastCell].flui.MasEspLiq(state.cells[state.lastCell].pres, state.cells[state.lastCell].temp);
        double rholc = state.cells[state.lastCell].fluicol.MasEspFlu(state.cells[state.lastCell].pres, state.cells[state.lastCell].temp);
        quality = fabs(massgas / masentrada);

        double masChk;

        double ypres = state.gasSurfacePressure / state.outletPressure;
        masChk = state.surfaceChoke.vazmassSachd(ypres, state.outletPressure, tESup, alfSup, betSup, quality, state.cells[state.lastCell - 1].flui,
                                       state.cells[state.lastCell - 1].fluicol);
        maxSup = state.surfaceChoke.vazmaxSachd(state.outletPressure, tESup, alfSup, betSup, quality, state.cells[state.lastCell - 1].flui, state.cells[state.lastCell - 1].fluicol);
        if (fabs(ypres) > fabs(state.surfaceChoke.razpres))
            maxSup = masChk;

        if (state.surfaceChoke.AreaGarg > (1e-3) * state.cells[state.lastCell - 1].duto.area && ypres < 1.) {
            double cplM = (1. - betSup) * state.cells[state.lastCell].flui.CalorLiq(state.outletPressure, tESup) -
                          betSup * state.cells[state.lastCell].fluicol.CalorLiq(state.outletPressure, tESup);
            double jtlM = (1. - betSup) * state.cells[state.lastCell].flui.JTL(state.outletPressure, tESup) - betSup / rholc;
            double gasSpecificHeat = state.cells[state.lastCell].flui.CalorGas(state.outletPressure, tESup);
            double jtgM = state.cells[state.lastCell].flui.JTG(state.outletPressure, tESup);
            state.input.valTempChokeJus = tESup + ((1. - quality) * jtlM / cplM + quality * jtgM / gasSpecificHeat) * (state.gasSurfacePressure - state.outletPressure) * 98066.52;
        }
    } else {
        maxSup = 0.;
        state.input.valTempChokeJus = state.cells[state.lastCell].temp;
    }

    if (fabs(maxSup) > 1e-15)
        state.baseConvergenceMonitor = fabs(maxSup);
    else
        state.baseConvergenceMonitor = 1.;
    return maxSup - state.cells[state.lastCell - 1].MR;
}

double marchProductionPressureToPressureTertiary(const SteadyStateState &state, double mchute) {

    double alfini = 0.;
    double betini = 0.;

    seedFirstCellFromFlowRateGuess(state, mchute, alfini, betini);

    state.cells[0].tempL = state.cells[0].temp;
    state.cells[1].tempL = state.cells[0].temp;
    state.cells[0].tempini = state.cells[0].temp;

    state.cells[0].ML = 0.;
    state.cells[0].MC = 0.;
    state.cells[1].ML = 0.;
    state.cells[0].MliqiniL = 0.;
    state.cells[0].Mliqini = 0.;
    state.cells[1].MliqiniL = 0.;
    state.cells[0].MComp = 0.;
    state.cells[0].QLL = 0.;
    state.cells[0].QL = 0;
    state.cells[1].QLL = 0.;
    state.cells[0].QG = 0.;

    int i;

    state.cells[0].presauxL = state.cells[0].pres;
    state.cells[0].presLini = state.cells[0].pres;
    state.cells[0].presL = state.cells[0].pres;
    state.cells[1].presL = state.cells[0].pres;
    state.cells[0].presini = state.cells[0].pres;
    state.cells[1].presLini = state.cells[0].pres;
    state.cells[0].presaux = state.cells[0].pres;
    state.cells[1].presauxL = state.cells[0].pres;

    state.cells[0].alf = alfini;
    state.cells[0].alfini = alfini;
    state.cells[0].bet = betini;
    state.cells[0].betini = betini;
    state.cells[1].alfL = state.cells[0].alf;
    state.cells[1].alfLini = state.cells[0].alf;
    state.cells[0].alfPigD = state.cells[0].alf;
    state.cells[0].alfPigDini = state.cells[0].alf;
    state.cells[0].alfPigE = state.cells[0].alf;
    state.cells[0].alfPigEini = state.cells[0].alf;
    state.cells[1].betL = state.cells[0].bet;
    state.cells[1].betLini = state.cells[0].bet;
    state.cells[0].betPigD = state.cells[0].bet;
    state.cells[0].betPigDini = state.cells[0].bet;
    state.cells[0].betPigE = state.cells[0].bet;
    state.cells[0].betPigEini = state.cells[0].bet;
    state.cells[0].betI = state.cells[0].bet;
    state.cells[1].betLI = state.cells[0].bet;

    state.steadyIteration = 0;
    while (state.steadyIteration < 3) {
        i = 1;
        while (i <= state.lastCell && state.cells[i - 1].pres >= 1.) {

            advanceUpstreamSteadyPressure(state, i, 0);
            refreshUpstreamProductionPeriphery(state, i);
            if (state.input.flashCompleto != 2)
                advanceSteadyMass(state, i);
            else
                advanceCompositionalSteadyMass(state, i);
            state.updaters.advanceSteadyTemperature(i, 0);
            if (isnan(state.cells[i].temp))
                NumError("Temperatrura na linha de producao com valor NaN");
            if (state.input.usaTabela == 1 && (state.cells[i].temp - state.input.tabent.tmin) < (*state.globals).localtiny)
                state.cells[i].temp = state.input.tabent.tmin;
            state.updaters.updateProductionTemperaturePeriphery(i);
            advanceDownstreamSteadyPressure(state, i, 0);
            refreshDownstreamProductionPeriphery(state, i);
            if (state.input.tipoFluido == 0)
                advanceSteadyMassTransfer(state, i - 1);
            else
                advanceSteadyGasMassTransfer(state, i - 1);
            if (state.input.ordperm > 1) {
                double D0presaux = state.cells[i].presaux - state.cells[i - 1].pres;
                double D0pres = state.cells[i].pres - state.cells[i].presaux;
                double D0temp = state.cells[i].temp - state.cells[i - 1].temp;
                advanceUpstreamSteadyPressure(state, i, 1);
                refreshUpstreamProductionPeriphery(state, i);
                advanceSteadyMass(state, i);
                state.updaters.advanceSteadyTemperature(i, 1);
                advanceDownstreamSteadyPressure(state, i, 1);
                state.cells[i].pres = 0.5 * (state.cells[i].presaux + D0pres + state.cells[i].pres);
                state.cells[i].presaux = 0.5 * (state.cells[i - 1].pres + D0presaux + state.cells[i].presaux);
                state.cells[i].temp = 0.5 * (state.cells[i - 1].temp + D0temp + state.cells[i].temp);
                refreshDownstreamProductionPeriphery(state, i);
                refreshUpstreamProductionPeriphery(state, i);
                state.updaters.updateProductionTemperaturePeriphery(i);
                advanceSteadyMass(state, i);
                if (state.input.tipoFluido == 0)
                    advanceSteadyMassTransfer(state, i - 1);
                else
                    advanceSteadyGasMassTransfer(state, i - 1);
            }

            i++;
            if (state.cells[i - 1].pres <= state.gasSurfacePressure)
                return -1e10;
        }
        state.steadyIteration++;
    }

    state.baseConvergenceMonitor = 1.;
    return 0. - state.cells[state.lastCell - 1].MR;
}

double marchGasSteadySecondary(const SteadyStateState &state, double pchute, double chutemass) {
    state.gasCells[0].presL = pchute;
    state.gasCells[0].pres = pchute;
    state.gasCells[0].presini = pchute;
    state.gasCells[1].presL = pchute;
    state.gasCells[0].tempL = state.initialGasTemperature;
    state.gasCells[0].temp = state.initialGasTemperature;
    state.gasCells[1].tempL = state.initialGasTemperature;
    state.gasCells[0].rg = state.gasCells[0].flui.MasEspGas(pchute, state.initialGasTemperature);
    state.gasCells[0].u1L = state.gasCells[0].duto.area * state.gasCells[0].rg;
    state.gasCells[0].u1LL = state.gasCells[0].u1L;
    state.gasCells[1].u1LL = state.gasCells[0].u1L;
    state.gasCells[0].VGasL = 0.;
    if (chutemass < 0) // se nenhum valor de chutemass for colocado na lista de parÃ£metro,
                       // usa o valor dado no json para a injecao de gas
        state.gasCells[0].massfonteCH = state.input.gasinj.vazgas[0] * state.gasCells[0].flui.MasEspGas(1., 15.6) / 86400.;
    else
        state.gasCells[0].massfonteCH = chutemass * state.gasCells[0].flui.MasEspGas(1., 15.6) / 86400.;
    state.gasCells[0].VGasR = state.gasCells[0].massfonteCH;
    state.gasCells[1].VGasL = state.gasCells[0].massfonteCH;

    for (int i = 1; i <= state.gasCellCount; i++) { // marcha na linha de servico
        state.updaters.updateSteadyGasPressure(i);            // avanco do valor de pressao de uma celula para outra, no centro da celula
        state.updaters.updateSteadyGasTemperature(i);            // avanco do valor de temperatura de uma celula para outra, centro da celula
        if (isnan(state.gasCells[i].temp))
            NumError("Temperatrura na linha de servico com valor NaN");
        state.updaters.computeSteadyGasFlowRate(i); // verifica se no centro desta celula tem uma VGL, calcula a vazao da VGL
        // retira este valor da vazÃƒÂ£o total na linha
        state.gasCells[i].rg = state.gasCells[i].flui.MasEspGas(state.gasCells[i].pres, state.gasCells[i].temp);
        state.gasCells[i - 1].rgR = state.gasCells[i].rg;
        state.gasCells[i].u1L = state.gasCells[i].duto.area * state.gasCells[i].rg;
        state.gasCells[i - 1].u1R = state.gasCells[i].u1L;
        if (i < state.gasCellCount)
            state.gasCells[i + 1].u1LL = state.gasCells[i].u1L;
    }
    state.gasCells[state.gasCellCount].rgR = state.gasCells[state.gasCellCount].rg;

    int valveCount = state.input.nvalvgas;
    double mastot = 0.;
    for (int j = 0; j < valveCount; j++)
        mastot += state.gasCells[state.gasValveCellIndices[j]].massfonteCH;
    return mastot - state.gasCells[0].massfonteCH; // diferenca entre a soma das vazoes nas VGL
    // e a vazao de injecao na linha
}

double marchGasSteadyTertiary(const SteadyStateState &state, double pchute) {
    state.gasCells[0].presL = pchute; // pressao a jusante do choque de injecao
    state.gasCells[0].pres = pchute;
    state.gasCells[0].presini = pchute;
    state.gasCells[1].presL = pchute;
    state.gasCells[0].tempL = state.initialGasTemperature;
    state.gasCells[0].temp = state.initialGasTemperature;
    state.gasCells[1].tempL = state.initialGasTemperature;
    state.gasCells[0].rg = state.gasCells[0].flui.MasEspGas(pchute, state.initialGasTemperature);
    state.gasCells[0].u1L = state.gasCells[0].duto.area * state.gasCells[0].rg;
    state.gasCells[0].u1LL = state.gasCells[0].u1L;
    state.gasCells[1].u1LL = state.gasCells[0].u1L;
    state.gasCells[0].VGasL = 0.;
    state.injectionChoke.presGarg = pchute;
    double chutemass = state.injectionChoke.massica(); // vazao de injecao na linha obtido a partir da vazao
    // massica do choke de injecao
    state.gasCells[0].massfonteCH = chutemass;
    state.gasCells[0].VGasR = state.gasCells[0].massfonteCH;
    state.gasCells[1].VGasL = state.gasCells[0].massfonteCH;

    for (int i = 1; i <= state.gasCellCount; i++) { // marcha na linha de gas
        state.updaters.updateSteadyGasPressure(i);            // avanco do valor de pressao de uma celula para outra, no centro da celula
        state.updaters.updateSteadyGasTemperature(i);            // avanco do valor de temperatura de uma celula para outra, centro da celula
        if (isnan(state.gasCells[i].temp))
            NumError("Temperatura na linha de servico com valor NaN");
        state.updaters.computeSteadyGasFlowRate(i); // verifica se no centro desta celula tem uma VGL, calcula a vazao da VGL
        // retira este valor da vazÃƒÂ£o total na linha
        state.gasCells[i].rg = state.gasCells[i].flui.MasEspGas(state.gasCells[i].pres, state.gasCells[i].temp);
        state.gasCells[i - 1].rgR = state.gasCells[i].rg;
        state.gasCells[i].u1L = state.gasCells[i].duto.area * state.gasCells[i].rg;
        state.gasCells[i - 1].u1R = state.gasCells[i].u1L;
        if (i < state.gasCellCount)
            state.gasCells[i + 1].u1LL = state.gasCells[i].u1L;
    }
    state.gasCells[state.gasCellCount].rgR = state.gasCells[state.gasCellCount].rg;

    int valveCount = state.input.nvalvgas;
    double mastot = 0.;
    for (int j = 0; j < valveCount; j++)
        mastot += state.gasCells[state.gasValveCellIndices[j]].massfonteCH;
    return mastot - chutemass; // diferenca entre a soma das vazoes nas VGL
    // e a vazao de injecao na linha
}

double marchInjectionSteady(const SteadyStateState &state, double chute) {

    int guessNeedsCorrection = 1;

    if (state.input.flashCompleto < 1)
        state.cells[0].temp = state.cells[0].acsr.injl.temp;
    else
        state.cells[0].temp = state.cells[0].acsr.injg.temp;
    double alfini = 0.;
    double betini = 1.;
    double delp = 0.;
    if (state.input.condpocinj.CC == 1 || state.input.condpocinj.CC == 2 || state.input.condpocinj.CC == 3) {
        if ((state.surfaceChoke.AreaGarg / state.surfaceChoke.AreaTub) >= 0.6) {
            if (state.input.flashCompleto < 1) {
                state.cells[0].acsr.injl.QLiq = chute;
            } else {
                state.cells[0].acsr.injg.QGas = chute;
                if (state.cells[0].acsr.injg.seco == 1) {
                    alfini = 1.;
                    betini = 0.;
                } else {
                    if (state.input.flashCompleto == 2) {
                        if (state.input.tabelaDinamica == 0)
                            state.cells[0].flui.atualizaPropComp(state.gasSurfacePressure, state.cells[0].acsr.injg.temp, -1, NULL, NULL, state.input.pocinjec);
                        state.cells[0].acsr.injg.FluidoPro.atualizaPropComp(state.gasSurfacePressure, state.cells[0].acsr.injg.temp, -1, NULL, NULL, state.input.pocinjec);
                    }
                    double masgas = state.cells[0].acsr.injg.VMas(state.gasSurfacePressure, state.cells[0].acsr.injg.temp);
                    double quality;
                    if (state.input.flashCompleto != 2)
                        quality = state.cells[0].acsr.injg.FluidoPro.FracMassHidra(1., 20.);
                    else
                        quality = state.cells[0].acsr.injg.FluidoPro.dStockTankVaporMassFraction;
                    double masT = masgas / quality;
                    quality = state.cells[0].acsr.injg.FluidoPro.FracMassHidra(state.gasSurfacePressure, state.cells[0].acsr.injg.temp);
                    double qgas = masT * quality /
                                  state.cells[0].acsr.injg.FluidoPro.MasEspGas(state.gasSurfacePressure, state.cells[0].acsr.injg.temp);
                    double qliq = masT * (1. - quality) /
                                  state.cells[0].acsr.injg.FluidoPro.MasEspLiq(state.gasSurfacePressure, state.cells[0].acsr.injg.temp);
                    alfini = qgas / (qliq + qgas);
                    betini = 0.;
                }
            }
        } else {
            if (state.input.flashCompleto < 1) {
                state.cells[0].acsr.injl.QLiq = chute;
                delp = (1 / (state.cells[0].fluicol.MasEspFlu(state.gasSurfacePressure, state.cells[0].acsr.injl.temp) * state.surfaceChoke.cdchk * 2.)) * pow((chute * state.cells[0].fluicol.MasEspFlu(1.01, 15.) / 86400) / (state.surfaceChoke.AreaGarg), 2.) / 98066.5;
            } else {
                state.cells[0].acsr.injg.QGas = chute;
                if (state.input.flashCompleto == 2) {
                    if (state.input.tabelaDinamica == 0)
                        state.cells[0].flui.atualizaPropComp(state.gasSurfacePressure, state.cells[0].acsr.injg.temp, -1, NULL, NULL, state.input.pocinjec);
                    state.cells[0].acsr.injg.FluidoPro.atualizaPropComp(state.gasSurfacePressure, state.cells[0].acsr.injg.temp, -1, NULL, NULL, state.input.pocinjec);
                }
                double rhogstd = state.cells[0].flui.Deng * 1.225;
                delp = (1 / (state.cells[0].flui.MasEspGas(state.gasSurfacePressure, state.cells[0].acsr.injg.temp) * state.surfaceChoke.cdchk * 2.)) * pow((chute * rhogstd / 86400) / (state.surfaceChoke.AreaGarg), 2.) / 98066.5;
                if (state.cells[0].acsr.injg.seco == 1) {
                    alfini = 1.;
                    betini = 0.;
                } else {
                    double masgas = state.cells[0].acsr.injg.VMas(state.gasSurfacePressure, state.cells[0].acsr.injg.temp);
                    double quality;
                    if (state.input.flashCompleto < 2)
                        quality = state.cells[0].acsr.injg.FluidoPro.FracMassHidra(1., 20.);
                    else
                        quality = state.cells[0].acsr.injg.FluidoPro.dStockTankVaporMassFraction;
                    double masT = masgas / quality;
                    quality = state.cells[0].acsr.injg.FluidoPro.FracMassHidra(state.gasSurfacePressure, state.cells[0].acsr.injg.temp);
                    double qgas = masT * quality /
                                  state.cells[0].acsr.injg.FluidoPro.MasEspGas(state.gasSurfacePressure, state.cells[0].acsr.injg.temp);
                    double qliq = masT * (1. - quality) /
                                  state.cells[0].acsr.injg.FluidoPro.MasEspLiq(state.gasSurfacePressure, state.cells[0].acsr.injg.temp);
                    alfini = qgas / (qliq + qgas);
                    betini = 0.;
                }
            }
        }
    }
    state.cells[0].tempL = state.cells[0].temp;
    state.cells[1].tempL = state.cells[0].temp;
    state.cells[0].tempini = state.cells[0].temp;

    state.cells[0].ML = 0.;
    state.cells[0].MC = 0.;
    state.cells[1].ML = 0.;
    state.cells[0].MliqiniL = 0.;
    state.cells[0].Mliqini = 0.;
    state.cells[1].MliqiniL = 0.;
    state.cells[0].QLL = 0.;
    state.cells[0].QL = 0;
    state.cells[1].QLL = 0.;
    state.cells[0].QG = 0.;
    if (state.input.condpocinj.CC == 0 || state.input.condpocinj.CC == 5) {
        state.gasSurfacePressure = chute;
    }
    state.steadyIteration = 0;

    double masfim = 0.;

    int i;
    guessNeedsCorrection = 1;
    while (guessNeedsCorrection == 1) {
        state.cells[0].presauxL = state.gasSurfacePressure - delp;
        state.cells[0].presLini = state.gasSurfacePressure - delp;
        state.cells[0].presL = state.gasSurfacePressure - delp;
        state.cells[0].pres = state.gasSurfacePressure - delp;
        state.cells[1].presL = state.gasSurfacePressure - delp;
        state.cells[0].presini = state.gasSurfacePressure - delp;
        state.cells[1].presLini = state.gasSurfacePressure - delp;
        state.cells[0].presaux = state.gasSurfacePressure - delp;
        state.cells[1].presauxL = state.gasSurfacePressure - delp;
        state.cells[0].alf = alfini;
        state.cells[0].alfini = alfini;
        state.cells[0].bet = betini;
        state.cells[0].betini = betini;
        state.cells[1].alfL = state.cells[0].alf;
        state.cells[1].alfLini = state.cells[0].alf;
        state.cells[0].alfPigD = state.cells[0].alf;
        state.cells[0].alfPigDini = state.cells[0].alf;
        state.cells[0].alfPigE = state.cells[0].alf;
        state.cells[0].alfPigEini = state.cells[0].alf;
        state.cells[1].betL = state.cells[0].bet;
        state.cells[1].betLini = state.cells[0].bet;
        state.cells[0].betPigD = state.cells[0].bet;
        state.cells[0].betPigDini = state.cells[0].bet;
        state.cells[0].betPigE = state.cells[0].bet;
        state.cells[0].betPigEini = state.cells[0].bet;
        state.cells[0].betI = state.cells[0].bet;
        state.cells[1].betLI = state.cells[0].bet;
        i = 1;
        while (i <= state.lastCell && state.cells[i - 1].pres >= 1. && masfim >= -(*state.globals).localtiny) {

            advanceUpstreamSteadyPressure(state, i, 0);
            refreshUpstreamProductionPeriphery(state, i);
            if (state.input.flashCompleto != 2)
                advanceSteadyMass(state, i);
            else
                advanceCompositionalSteadyMass(state, i);
            state.updaters.advanceSteadyTemperature(i, 0);
            state.updaters.updateProductionTemperaturePeriphery(i);
            advanceDownstreamSteadyPressure(state, i, 0);
            refreshDownstreamProductionPeriphery(state, i);
            advanceSteadyMassTransfer(state, i - 1);
            if (state.input.ordperm > 1) { // correcao de segunda ordem
                double D0presaux = state.cells[i].presaux - state.cells[i - 1].pres;
                double D0pres = state.cells[i].pres - state.cells[i].presaux;
                double D0temp = state.cells[i].temp - state.cells[i - 1].temp;
                advanceUpstreamSteadyPressure(state, i, 1);
                refreshUpstreamProductionPeriphery(state, i);
                advanceSteadyMass(state, i);
                state.updaters.advanceSteadyTemperature(i, 1);
                advanceDownstreamSteadyPressure(state, i, 1);
                state.cells[i].pres = 0.5 * (state.cells[i].presaux + D0pres + state.cells[i].pres);
                state.cells[i].presaux = 0.5 * (state.cells[i - 1].pres + D0presaux + state.cells[i].presaux);
                state.cells[i].temp = 0.5 * (state.cells[i - 1].temp + D0temp + state.cells[i].temp);
                refreshDownstreamProductionPeriphery(state, i);
                refreshUpstreamProductionPeriphery(state, i);
                state.updaters.updateProductionTemperaturePeriphery(i);
                advanceSteadyMass(state, i);
                advanceSteadyMassTransfer(state, i - 1);
            }
            masfim += (state.cells[i - 1].fontemassCR + state.cells[i - 1].fontemassLR + state.cells[i - 1].fontemassGR);
            i++;
            if (state.cells[i - 1].pres <= 1 || (state.cells[i - 1].acsr.tipo == 3 && (state.cells[i - 1].acsr.ipr.Pres - state.cells[i - 1].pres) > -(*state.globals).localtiny))
                return -1e10;
            else if (i < (state.lastCell + 1) && masfim < -(*state.globals).localtiny)
                return 1e10;
        }
        if (i == state.lastCell + 1)
            guessNeedsCorrection = 0;
    }
    state.updaters.updateSource(state.lastCell);
    masfim += (state.cells[state.lastCell].fontemassCR + state.cells[state.lastCell].fontemassLR + state.cells[state.lastCell].fontemassGR);

    if (state.input.condpocinj.CC != 3 && state.input.condpocinj.CC != 5)
        return masfim;
    else
        return state.input.condpocinj.presfundo - state.cells[state.lastCell].pres;
}

double reverseHydrostatic(const SteadyStateState &state, double liquidHoldup, double liquidFlowRate, double gasFlowRate) {
    double pchute = state.gasSurfacePressure;
    double taux;
    state.cells[state.lastCell].pres = pchute;
    double j = 0.;
    double rmis = 0.;
    double frictionFactor = 0.;
    for (int i = state.lastCell; i > 0; i--) {
        if (liquidFlowRate > 0. || gasFlowRate > 0.) {
            taux = state.input.celp[0].textern;
            double completionFraction = 0.;
            double visC = state.cells[i].fluicol.VisFlu(pchute, taux);
            double visP = state.cells[i].flui.ViscOleo(pchute, taux);
            double visG = state.cells[i].flui.ViscGas(pchute, taux);
            double visMis = (1 - completionFraction) * visP + completionFraction * visC;
            double completionDensityAtGuess = state.cells[i].fluicol.MasEspFlu(pchute, taux);
            double liquidDensityAtGuess = state.cells[i].flui.MasEspLiq(pchute, taux);
            double gasDensityAtGuess = state.cells[i].flui.MasEspGas(pchute, taux);
            rmis = (1 - completionFraction) * liquidDensityAtGuess + completionFraction * completionDensityAtGuess;
            double rlpA = state.cells[i].flui.MasEspLiq(1., 15.);
            double rlcA = state.cells[i].fluicol.MasEspFlu(1.001, 15.);
            double massicC = rlcA * liquidFlowRate * completionFraction;
            double massic = rlpA * liquidFlowRate * (1. - completionFraction);
            double Rhogs = state.cells[i].flui.Deng * 1.225;
            double Rhols = (1000 * 141.5 / (131.5 + state.cells[i].flui.API)) * (1 - state.cells[i].flui.BSW) + 1000. * state.cells[i].flui.Denag * state.cells[i].flui.BSW;
            double multiplicador = (Rhols + state.cells[i].flui.RGO * Rhogs * (1 - state.cells[i].flui.BSW));
            massic = 1 * liquidFlowRate * (1. - completionFraction) * multiplicador;
            double fracmasshidra = state.cells[i].flui.FracMassHidra(pchute, taux);
            double massicP = (1. - fracmasshidra) * massic;
            double massicG = fracmasshidra * massic + gasFlowRate * Rhogs;
            j = (massicP / liquidDensityAtGuess + massicC / completionDensityAtGuess + massicG / gasDensityAtGuess) / state.cells[i].duto.area;
            double alfmis = (massicG / gasDensityAtGuess) / (massicP / liquidDensityAtGuess + massicC / completionDensityAtGuess + massicG / gasDensityAtGuess);
            rmis = (1 - alfmis) * rmis + alfmis * gasDensityAtGuess;
            visMis = (1 - alfmis) * visMis + alfmis * visG;
            double mixtureReynolds;
            if (state.cells[i].duto.revest == 0)
                mixtureReynolds = state.cells[0].Rey(state.cells[i].duto.a, j, rmis, visMis);
            else {
                double dhid = 4 * state.cells[i].duto.area / state.cells[0].duto.peri;
                mixtureReynolds = state.cells[i].Rey(dhid, j, rmis, visMis);
            }
            frictionFactor = state.cells[i].fric(mixtureReynolds, state.cells[0].duto.rug / state.cells[0].duto.a);
        }
        double perdafric = (frictionFactor * rmis * j * fabs(j) / 2.) * state.cells[i].duto.peri / state.cells[i].duto.area;
        taux = state.input.celp[i].textern;
        if (i == 500) {
            int para;
            para = 0;
        }
        double rhol = state.cells[i].flui.MasEspLiq(pchute, taux);
        double rhog = state.cells[i].flui.MasEspGas(pchute, taux);
        double alfa = 1. - liquidHoldup;
        double rhomix = (1. - alfa) * rhol + alfa * rhog;
        double dxmed = 0.5 * (state.cells[i].dx + state.cells[i - 1].dx);
        pchute += ((rhomix * 9.81 * sin(state.cells[i].duto.teta) * dxmed + perdafric * dxmed) / 98066.5);
        if (state.cells[i - 1].acsr.tipo == 7)
            pchute -= state.cells[i - 1].acsr.delp;
        if (state.cells[i - 1].acsr.tipo == 3 && (state.cells[i - 1].acsr.ipr.Pres - pchute) < (*state.globals).localtiny)
            pchute = 0.99 * state.cells[i - 1].acsr.ipr.Pres;
        if (state.cells[i - 1].acsr.tipo == 15 && (state.cells[i - 1].acsr.radialPoro.pRes[0] - pchute) < (*state.globals).localtiny)
            pchute = 0.99 * state.cells[i - 1].acsr.radialPoro.pRes[0];
        if (state.cells[i - 1].acsr.tipo == 16 && (state.cells[i - 1].acsr.poroso2D.dados.pRes - pchute) < (*state.globals).localtiny)
            pchute = 0.99 * state.cells[i - 1].acsr.poroso2D.dados.pRes;

        state.cells[i - 1].dpB = 0.;
        if (state.cells[i - 1].acsr.tipo == 4 && state.cells[i - 1].acsr.bcs.freqnova > 1. && liquidFlowRate >= 0.) {
            double vazmix = j * state.cells[i - 1].dutoL.area;
            double rhomis = state.cells[i - 1].flui.MasEspLiq(pchute, taux);
            double vismis = state.cells[i - 1].flui.ViscOleo(pchute, taux);
            vazmix *= (86400 / 0.1589876);
            state.cells[i - 1].acsr.bcs.NovaVis(vismis, rhomis, vazmix);
            state.cells[i - 1].dpB = 0.3048 * state.cells[i - 1].acsr.bcs.Hvis * rhomis * 9.82;
        }
        if (state.cells[i - 1].acsr.tipo == 17 && state.cells[i - 1].acsr.multibcs.freqnova > 1. && liquidFlowRate >= 0.) {
            double alf0 = state.cells[i - 1].alf;
            double bet0 = state.cells[i - 1].bet;
            state.cells[i - 1].acsr.multibcs.flui = state.cells[i - 1].flui;
            state.cells[i - 1].acsr.multibcs.fluicol = state.cells[i - 1].fluicol;
            state.cells[i - 1].acsr.multibcs.marchaMultiBcs(state.cells[i - 1].QG, state.cells[i - 1].QL,
                                                       pchute, taux, alf0, bet0);
            state.cells[i - 1].dpB = state.cells[i - 1].acsr.multibcs.dpB * 98066.52;
        }
        pchute -= state.cells[i - 1].dpB / 98066.5;
        state.cells[i - 1].pres = pchute;
    }
    return pchute;
}

double reverseInjectionHydrostatic(const SteadyStateState &state, double liquidHoldup, double liquidFlowRate) {
    double pchute = 0.;
    if (state.input.condpocinj.presfundo > 1e-5)
        pchute = state.input.condpocinj.presfundo;
    else if (state.cells[state.lastCell].acsr.tipo == 3)
        pchute = state.cells[state.lastCell].acsr.ipr.Pres;
    else
        NumError("Sem pressao no fim do tramo e sem IPR-metodo hidroreversoInj");
    double taux;
    double rmis = 0.;
    double j = 0.;
    double frictionFactor = 0.;

    if (liquidFlowRate > 0.) {
        taux = state.input.celp[0].textern;
        state.cells[state.lastCell].pres = pchute;
        double visC = state.cells[0].acsr.injl.fluidocol.VisFlu(pchute, taux);
        double visMis = visC;
        double completionDensityAtGuess = state.cells[0].acsr.injl.fluidocol.MasEspFlu(pchute, taux);
        rmis = completionDensityAtGuess;
        double rlcA = state.cells[0].acsr.injl.fluidocol.MasEspFlu(1.001, 15.);
        double massicC = rlcA * liquidFlowRate;
        j = (massicC / completionDensityAtGuess) / state.cells[0].duto.area;
        double mixtureReynolds;
        if (state.cells[0].duto.revest == 0)
            mixtureReynolds = state.cells[0].Rey(state.cells[0].duto.a, j, rmis, visMis);
        else {
            double dhid = 4 * state.cells[0].duto.area / state.cells[0].duto.peri;
            mixtureReynolds = state.cells[0].Rey(dhid, j, rmis, visMis);
        }
        frictionFactor = state.cells[0].fric(mixtureReynolds, state.cells[0].duto.rug / state.cells[0].duto.a);
    }
    double perdafric = (frictionFactor * rmis * j * fabs(j) / 2.) * state.cells[0].duto.peri / state.cells[0].duto.area;
    for (int i = state.lastCell; i > 0; i--) {
        taux = state.input.celp[i].textern;
        double rhol = state.cells[i].fluicol.MasEspFlu(pchute, taux);
        double rhog = 0.;
        double alfa = 1. - liquidHoldup;
        double rhomix = (1. - alfa) * rhol + alfa * rhog;
        double dxmed = 0.5 * (state.cells[i].dx + state.cells[i - 1].dx);
        pchute += ((rhomix * 9.81 * sin(state.cells[i].duto.teta) * dxmed + perdafric * dxmed) / 98066.5);
        if (state.cells[i - 1].acsr.tipo == 7)
            pchute -= state.cells[i - 1].acsr.delp;
        if (state.cells[i - 1].acsr.tipo == 3 && (state.cells[i - 1].acsr.ipr.Pres - pchute) < (*state.globals).localtiny)
            pchute = 0.99 * state.cells[i - 1].acsr.ipr.Pres;
        state.cells[i - 1].pres = pchute;
    }
    return pchute;
}

double secondaryBranchHydrostatic(const SteadyStateState &state, double titulo) {
    double pchute;
    if (state.input.ConContEntrada == 1)
        pchute = state.input.CCPres.pres[0];
    else
        pchute = state.cells[0].pres;
    double taux;
    double liquidHoldup;
    double completionFraction;
    if (state.input.ConContEntrada == 1)
        completionFraction = state.input.CCPres.bet[0];
    else
        completionFraction = state.cells[0].bet;
    state.cells[0].bet = completionFraction;
    taux = state.input.celp[0].textern;
    double titRef = state.cells[0].flui.FracMassHidra(pchute, taux);
    if (state.input.ConContEntrada == 1)
        state.cells[0].pres = state.input.CCPres.pres[0];
    state.cells[0].temp = taux;
    double rhol = state.cells[0].flui.MasEspLiq(state.cells[0].pres, taux);
    double rhoc = state.cells[0].fluicol.MasEspFlu(state.cells[0].pres, taux);
    double rhog = state.cells[0].flui.MasEspGas(state.cells[0].pres, taux);
    double titEntra;
    if (state.input.ConContEntrada == 1)
        titEntra = state.input.CCPres.tit[0];
    else
        titEntra = titulo;
    state.cells[0].alf = titEntra * (rhol * (1. - completionFraction) + rhoc * completionFraction) / (rhog * (1. - titEntra) + titEntra * (rhol * (1. - completionFraction) + rhoc * completionFraction));
    for (int i = 0; i < state.lastCell; i++) {
        taux = state.input.celp[i].textern;
        rhol = state.cells[i].flui.MasEspLiq(pchute, taux);
        rhoc = state.cells[i].fluicol.MasEspFlu(pchute, taux);
        rhog = state.cells[i].flui.MasEspGas(pchute, taux);
        double tit0 = state.cells[i].flui.FracMassHidra(pchute, taux);
        double delTit = tit0 - titRef;
        double novoTit = titulo + delTit * (titulo / titRef);
        double alfa = novoTit * (rhol * (1. - completionFraction) + rhoc * completionFraction) / (rhog * (1. - novoTit) + novoTit * (rhol * (1. - completionFraction) + rhoc * completionFraction));
        state.cells[i + 1].alf = alfa;
        state.cells[i + 1].bet = completionFraction;

        double rhomix = (1. - alfa) * ((1. - completionFraction) * rhol + completionFraction * rhoc) + alfa * rhog;
        double dxmed = 0.5 * (state.cells[i].dx + state.cells[i + 1].dx);
        pchute -= ((rhomix * 9.81 * sin(state.cells[i].duto.teta) * dxmed) / 98066.5);
        state.cells[i + 1].pres = pchute;
        state.cells[i + 1].temp = taux;
    }
    return pchute;
}

}  // namespace sisprod::steady
