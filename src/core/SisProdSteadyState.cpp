#include "SisProdSteadyState.h"

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
            for (int iRP = 0; iRP < state.cells[i].acsr.radialPoro.ncel; iRP++) {
                state.cells[i].acsr.radialPoro.celula[iRP].flup.rDgD = 1.;
                state.cells[i].acsr.radialPoro.celula[iRP].flup.rDgL = 1.;
            }
        }
        if (state.cells[i].acsr.tipo == 16) {
            state.cells[i].acsr.poroso2D.dados.flup.rDgD = 1.;
            state.cells[i].acsr.poroso2D.dados.flup.rDgL = 1.;
            for (int iRP = 0; iRP < state.cells[i].acsr.poroso2D.dados.transfer.ncel; iRP++) {
                state.cells[i].acsr.poroso2D.dados.transfer.celula[iRP].flup.rDgD = 1.;
                state.cells[i].acsr.poroso2D.dados.transfer.celula[iRP].flup.rDgL = 1.;
            }
            for (int iRP = 0; iRP < state.cells[i].acsr.poroso2D.malha.nele; iRP++) {
                state.cells[i].acsr.poroso2D.malha.mlh2d[iRP].flup.rDgD = 1.;
                state.cells[i].acsr.poroso2D.malha.mlh2d[iRP].flup.rDgL = 1.;
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
            for (int iRP = 0; iRP < state.cells[i].acsr.radialPoro.ncel; iRP++) {
                double pres = state.cells[i].acsr.radialPoro.celula[iRP].Pcamada;
                double temp = state.cells[i].acsr.radialPoro.tRes;
                state.cells[i].acsr.radialPoro.celula[iRP].flup.razDegD(pres, temp);
                state.cells[i].acsr.radialPoro.celula[iRP].flup.rzDegL(pres, temp);
            }
        }
        if (state.cells[i].acsr.tipo == 16) {
            state.cells[i].acsr.poroso2D.dados.flup.razDegD(state.cells[i].pres, state.cells[i].temp);
            state.cells[i].acsr.poroso2D.dados.flup.rzDegL(state.cells[i].pres, state.cells[i].temp);
            for (int iRP = 0; iRP < state.cells[i].acsr.poroso2D.dados.transfer.ncel; iRP++) {
                double pres = state.cells[i].acsr.poroso2D.dados.transfer.celula[iRP].Pcamada;
                double temp = state.cells[i].acsr.poroso2D.dados.transfer.tRes;
                state.cells[i].acsr.poroso2D.dados.transfer.celula[iRP].flup.razDegD(pres, temp);
                state.cells[i].acsr.poroso2D.dados.transfer.celula[iRP].flup.rzDegL(pres, temp);
            }
            for (int iRP = 0; iRP < state.cells[i].acsr.poroso2D.malha.nele; iRP++) {
                double pres = state.cells[i].acsr.poroso2D.malha.mlh2d[iRP].cel2D.presC;
                double temp = state.cells[i].acsr.poroso2D.malha.mlh2d[iRP].tRes;
                state.cells[i].acsr.poroso2D.malha.mlh2d[iRP].flup.razDegD(pres, temp);
                state.cells[i].acsr.poroso2D.malha.mlh2d[iRP].flup.rzDegL(pres, temp);
            }
        }
    }
}

void advanceReverseSteadyMass(const SteadyStateState &state, int i) {
    int mudaRGO = 1;
    if (state.input.flashCompleto == 1)
        mudaRGO = 1;
    ProFlu fluF;
    double trF = 0.;
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
            trF = state.cells[i - 1].acsr.injg.fluidocol.TR;
        } else if (state.cells[i - 1].acsr.tipo == 2) {
            fluF = state.cells[i - 1].acsr.injl.FluidoPro;
            trF = state.cells[i - 1].acsr.injl.fluidocol.TR;
        } else if (state.cells[i - 1].acsr.tipo == 3) {
            fluF = state.cells[i - 1].acsr.ipr.FluidoPro;
            trF = 0.;
        } else if (state.cells[i - 1].acsr.tipo == 15) {
            fluF = state.cells[i - 1].acsr.radialPoro.flup;
            trF = 0.;
        } else if (state.cells[i - 1].acsr.tipo == 16) {
            fluF = state.cells[i - 1].acsr.poroso2D.dados.flup;
            trF = 0.;
        } else if (state.cells[i - 1].acsr.tipo == 10) {
            fluF = state.cells[i - 1].acsr.injm.FluidoPro;
            trF = state.cells[i - 1].acsr.injm.fluidocol.TR;
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
    double boI;
    double baI;
    if (state.cells[i].flui.RGO < 1e6)
        boI = state.cells[i].flui.BOFunc(state.cells[i].presaux, tmed);
    else
        boI = 1.;
    baI = state.cells[i].flui.BAFunc(state.cells[i].presaux, tmed);
    double fwI = state.cells[i].flui.BSW * baI / (boI + baI * state.cells[i].flui.BSW - state.cells[i].flui.BSW * boI);
    double qo = MasLiqProd * (1 - state.cells[i - 1].FW) / state.cells[i].rpCi;
    double qw;
    if (fwI < (1 - (*state.globals).localtiny))
        qw = fwI * qo / (1 - fwI);
    else
        qw = MasLiqProd / state.cells[i].flui.MasEspAgua(pmed, tmed);
    double qc;
    qc = state.cells[i].MComp / state.cells[i].rcCi;
    ////////////////////////esperar//////////////////////////////////////
    state.cells[i].bet = qc / (qo + qw + qc);

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
    double tL;
    if (state.input.flashCompleto == 1 && (state.input.tabent.tmin - 0) > (*state.globals).localtiny)
        tL = state.input.tabent.tmin + 0.1;
    else
        tL = 0.;
    double tH;
    if (state.input.flashCompleto == 1 && (70 - state.input.tabent.tmax) < (*state.globals).localtiny)
        tH = state.input.tabent.tmax - 0.1;
    else
        tH = 70.;
    double bo;
    double ba;
    double rs;
    // calcula os valores de RS, Bo e Ba na celula anterior a i-esima celula
    if (state.cells[i - 1].flui.RGO < 1e7) {
        rs = state.cells[i - 1].flui.RS(state.cells[i - 1].pres, state.cells[i - 1].temp);
        bo = state.cells[i - 1].flui.BOFunc(state.cells[i - 1].pres, state.cells[i - 1].temp, rs);
        ba = state.cells[i - 1].flui.BAFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        rs = rs * 6.29 / 35.31467;
    } else {
        bo = 1;
        rs = 0;
        ba = 0.;
    }
    // BSW in-situ da celula anterior, na marcha, a i-esima celula
    state.cells[i - 1].FW = state.cells[i - 1].flui.BSW * ba / (bo + ba * state.cells[i - 1].flui.BSW - state.cells[i - 1].flui.BSW * bo);
    state.cells[i - 1].FWini = state.cells[i - 1].FW;
    // Definicao das fracoes volumetricas:
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

    if (fabs(state.cells[i].QL) > 1e-15 && fabs(state.cells[i].MComp) > 1e-15) {
        double hol0 = 1. - state.cells[i - 1].alf;
        double hol1 = 1. - state.cells[i].alf;
        state.cells[i].fluicol.TR = (state.cells[i - 1].fluicol.TR * state.cells[i - 1].MComp + state.cells[i - 1].fontemassCR * trF) / state.cells[i].MComp;
        state.cells[i].fluicol.TR = state.cells[i].fluicol.TR + (0.5 * state.cells[i - 1].duto.area * state.cells[i - 1].dx * hol0 +
                                                       0.5 * state.cells[i].duto.area * state.cells[i].dx * hol1) /
                                                          state.cells[i].QL;
    } else
        state.cells[i].fluicol.TR = 0.;
}

}  // namespace sisprod::steady
