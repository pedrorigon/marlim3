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

void advanceReverseCompositionalSteadyMass(const SteadyStateState &state, int i) {
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
            state.cells[i - 1].acsr.injg.FluidoPro.atualizaPropComp(state.cells[i].pres, state.cells[i].temp, -1, NULL, NULL, state.input.pocinjec);
            fluF = state.cells[i - 1].acsr.injg.FluidoPro;
            trF = state.cells[i - 1].acsr.injg.fluidocol.TR;
        } else if (state.cells[i - 1].acsr.tipo == 2) {
            state.cells[i - 1].acsr.injl.FluidoPro.atualizaPropComp(state.cells[i].pres, state.cells[i].temp, -1, NULL, NULL, state.input.pocinjec);
            fluF = state.cells[i - 1].acsr.injl.FluidoPro;
            trF = state.cells[i - 1].acsr.injl.fluidocol.TR;
        } else if (state.cells[i - 1].acsr.tipo == 3) {
            state.cells[i - 1].acsr.ipr.FluidoPro.atualizaPropComp(state.cells[i].pres, state.cells[i].temp, -1, NULL, NULL, state.input.pocinjec);
            fluF = state.cells[i - 1].acsr.ipr.FluidoPro;
            trF = 0.;
        } else if (state.cells[i - 1].acsr.tipo == 10) {
            state.cells[i - 1].acsr.injm.FluidoPro.atualizaPropComp(state.cells[i].pres, state.cells[i].temp, -1, NULL, NULL, state.input.pocinjec);
            fluF = state.cells[i - 1].acsr.injm.FluidoPro;
            trF = state.cells[i - 1].acsr.injm.fluidocol.TR;
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
            trF = 0.;
        } else if (state.cells[i - 1].acsr.tipo == 16) {
            state.cells[i - 1].acsr.poroso2D.dados.flup.atualizaPropComp(state.cells[i].pres, state.cells[i].temp, -1, NULL, NULL, state.input.pocinjec);
            fluF = state.cells[i - 1].acsr.poroso2D.dados.flup;
            trF = 0.;
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

    double titulo = state.cells[i].flui.dVaporMassFraction;

    state.cells[i].rpCi = state.cells[i].flui.MasEspLiq(pmed, tmed);
    state.cells[i].rcCi = state.cells[i].fluicol.MasEspFlu(pmed, tmed);
    state.cells[i].rgCi = state.cells[i].flui.MasEspGas(pmed, tmed);

    double MasCarb = state.cells[i].MC - state.cells[i].MComp;
    double MasGas = MasCarb * state.cells[i].flui.FracMassHidra(state.cells[i].presaux, tmed);
    double MasLiqProd = MasCarb - MasGas;
    double boI;
    double baI;
    if (titulo < 1 - 1e-15)
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
    if (titulo < 1 - 1e-15) {
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
        double dxmed = 0.5 * (state.cells[i - 1].dx + state.cells[i].dx);
        double hol0 = 1. - state.cells[i - 1].alf;
        double hol1 = 1. - state.cells[i].alf;
        state.cells[i].fluicol.TR = (state.cells[i - 1].fluicol.TR * state.cells[i - 1].MComp + state.cells[i - 1].fontemassCR * trF) / state.cells[i].MComp;
        state.cells[i].fluicol.TR = state.cells[i].fluicol.TR + (0.5 * state.cells[i - 1].duto.area * state.cells[i - 1].dx * hol0 +
                                                       0.5 * state.cells[i].duto.area * state.cells[i].dx * hol1) /
                                                          state.cells[i].QL;
    } else
        state.cells[i].fluicol.TR = state.cells[i - 1].fluicol.TR;
}

void advanceCompositionalSteadyMass(const SteadyStateState &state, int i) {
    int mudaRGO = 1;
    double titF = 0.;
    ProFlu fluF;

    double boF = 1.;
    double baF = 1.;
    double fwF = 1.;
    double rhoOF = 900.;
    double rhoWF = 1000.;

    double mHidro = 0.;
    double mComp;
    double trF = 0.;

    if (state.cells[i - 1].acsr.tipo == 1) {
        if (i > 0 && (state.cells[i - 1].flui.dCalculatedBeta > 0. && state.cells[i - 1].flui.dCalculatedBeta < 1.))
            state.cells[i - 1].acsr.injg.FluidoPro.atualizaPropComp(state.cells[i - 1].pres, state.cells[i - 1].temp, state.cells[i - 1].flui.dCalculatedBeta,
                                                               state.cells[i - 1].flui.oCalculatedLiqComposition, state.cells[i - 1].flui.oCalculatedVapComposition, state.input.pocinjec);
        else
            state.cells[i - 1].acsr.injg.FluidoPro.atualizaPropComp(state.cells[i - 1].pres, state.cells[i - 1].temp, -1, NULL, NULL, state.input.pocinjec);
        fluF = state.cells[i - 1].acsr.injg.FluidoPro;
        fwF = 0.;
        titF = 1.;
        trF = state.cells[i - 1].acsr.injg.fluidocol.TR;
    } else if (state.cells[i - 1].acsr.tipo == 2) {
        if (i > 0 && (state.cells[i - 1].flui.dCalculatedBeta > 0. && state.cells[i - 1].flui.dCalculatedBeta < 1.))
            state.cells[i - 1].acsr.injl.FluidoPro.atualizaPropComp(state.cells[i - 1].pres, state.cells[i - 1].temp, state.cells[i - 1].flui.dCalculatedBeta,
                                                               state.cells[i - 1].flui.oCalculatedLiqComposition, state.cells[i - 1].flui.oCalculatedVapComposition, state.input.pocinjec);
        else
            state.cells[i - 1].acsr.injl.FluidoPro.atualizaPropComp(state.cells[i - 1].pres, state.cells[i - 1].temp, -1, NULL, NULL, state.input.pocinjec);
        fluF = state.cells[i - 1].acsr.injl.FluidoPro;
        boF = fluF.BOFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        baF = fluF.BAFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        fwF = fluF.BSW * baF / (boF + baF * fluF.BSW - fluF.BSW * boF);
        rhoOF = fluF.MasEspoleo(state.cells[i - 1].pres, state.cells[i - 1].temp);
        rhoWF = fluF.MasEspAgua(state.cells[i - 1].pres, state.cells[i - 1].temp);
        titF = (1 - fwF) * rhoOF / ((1 - fwF) * rhoOF + fwF * rhoWF);
        trF = state.cells[i - 1].acsr.injl.fluidocol.TR;
    } else if (state.cells[i - 1].acsr.tipo == 3) {
        if (i > 0 && (state.cells[i - 1].flui.dCalculatedBeta > 0. && state.cells[i - 1].flui.dCalculatedBeta < 1.))
            state.cells[i - 1].acsr.ipr.FluidoPro.atualizaPropComp(state.cells[i - 1].pres, state.cells[i - 1].temp, state.cells[i - 1].flui.dCalculatedBeta,
                                                              state.cells[i - 1].flui.oCalculatedLiqComposition, state.cells[i - 1].flui.oCalculatedVapComposition, state.input.pocinjec);
        else
            state.cells[i - 1].acsr.ipr.FluidoPro.atualizaPropComp(state.cells[i - 1].pres, state.cells[i - 1].temp, -1, NULL, NULL, state.input.pocinjec);
        fluF = state.cells[i - 1].acsr.ipr.FluidoPro;
        boF = fluF.BOFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        baF = fluF.BAFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        fwF = fluF.BSW * baF / (boF + baF * fluF.BSW - fluF.BSW * boF);
        rhoOF = fluF.MasEspoleo(state.cells[i - 1].pres, state.cells[i - 1].temp);
        rhoWF = fluF.MasEspAgua(state.cells[i - 1].pres, state.cells[i - 1].temp);
        titF = (1 - fwF) * rhoOF / ((1 - fwF) * rhoOF + fwF * rhoWF);
        trF = 0.;
    } else if (state.cells[i - 1].acsr.tipo == 10) {
        if (i > 0 && (state.cells[i - 1].flui.dCalculatedBeta > 0. && state.cells[i - 1].flui.dCalculatedBeta < 1.))
            state.cells[i - 1].acsr.injm.FluidoPro.atualizaPropComp(state.cells[i - 1].pres, state.cells[i - 1].temp, state.cells[i - 1].flui.dCalculatedBeta,
                                                               state.cells[i - 1].flui.oCalculatedLiqComposition, state.cells[i - 1].flui.oCalculatedVapComposition, state.input.pocinjec);
        else
            state.cells[i - 1].acsr.injm.FluidoPro.atualizaPropComp(state.cells[i - 1].pres, state.cells[i - 1].temp, -1, NULL, NULL, state.input.pocinjec);
        fluF = state.cells[i - 1].acsr.injm.FluidoPro;
        boF = fluF.BOFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        baF = fluF.BAFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        fwF = fluF.BSW * baF / (boF + baF * fluF.BSW - fluF.BSW * boF);
        rhoOF = fluF.MasEspoleo(state.cells[i - 1].pres, state.cells[i - 1].temp);
        rhoWF = fluF.MasEspAgua(state.cells[i - 1].pres, state.cells[i - 1].temp);
        titF = (1 - fwF) * rhoOF / ((1 - fwF) * rhoOF + fwF * rhoWF);
        trF = state.cells[i - 1].acsr.injm.fluidocol.TR;
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
        boF = fluF.BOFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        baF = fluF.BAFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        fwF = fluF.BSW * baF / (boF + baF * fluF.BSW - fluF.BSW * boF);
        rhoOF = fluF.MasEspoleo(state.cells[i - 1].pres, state.cells[i - 1].temp);
        rhoWF = fluF.MasEspAgua(state.cells[i - 1].pres, state.cells[i - 1].temp);
        titF = (1 - fwF) * rhoOF / ((1 - fwF) * rhoOF + fwF * rhoWF);
    } else if (state.cells[i - 1].acsr.tipo == 15) {
        double tRes = state.cells[i - 1].acsr.radialPoro.tRes;
        if (i > 0 && (state.cells[i - 1].flui.dCalculatedBeta > 0. && state.cells[i - 1].flui.dCalculatedBeta < 1.))
            state.cells[i - 1].acsr.radialPoro.flup.atualizaPropComp(state.cells[i - 1].pres, tRes,
                                                                state.cells[i - 1].flui.dCalculatedBeta,
                                                                state.cells[i - 1].flui.oCalculatedLiqComposition, state.cells[i - 1].flui.oCalculatedVapComposition, state.input.pocinjec);
        else
            state.cells[i - 1].acsr.radialPoro.flup.atualizaPropComp(state.cells[i - 1].pres, tRes, -1, NULL, NULL, state.input.pocinjec);
        fluF = state.cells[i - 1].acsr.radialPoro.flup;
        boF = fluF.BOFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        baF = fluF.BAFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        fwF = fluF.BSW * baF / (boF + baF * fluF.BSW - fluF.BSW * boF);
        rhoOF = fluF.MasEspoleo(state.cells[i - 1].pres, state.cells[i - 1].temp);
        rhoWF = fluF.MasEspAgua(state.cells[i - 1].pres, state.cells[i - 1].temp);
        titF = (1 - fwF) * rhoOF / ((1 - fwF) * rhoOF + fwF * rhoWF);
        trF = 0.;
    } else if (state.cells[i - 1].acsr.tipo == 16) {
        double tRes = state.cells[i - 1].acsr.poroso2D.dados.transfer.tRes;
        if (i > 0 && (state.cells[i - 1].flui.dCalculatedBeta > 0. && state.cells[i - 1].flui.dCalculatedBeta < 1.))
            state.cells[i - 1].acsr.poroso2D.dados.flup.atualizaPropComp(state.cells[i - 1].pres, tRes,
                                                                    state.cells[i - 1].flui.dCalculatedBeta,
                                                                    state.cells[i - 1].flui.oCalculatedLiqComposition, state.cells[i - 1].flui.oCalculatedVapComposition, state.input.pocinjec);
        else
            state.cells[i - 1].acsr.poroso2D.dados.flup.atualizaPropComp(state.cells[i - 1].pres, tRes, -1, NULL, NULL, state.input.pocinjec);
        fluF = state.cells[i - 1].acsr.poroso2D.dados.flup;
        boF = fluF.BOFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        baF = fluF.BAFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        fwF = fluF.BSW * baF / (boF + baF * fluF.BSW - fluF.BSW * boF);
        rhoOF = fluF.MasEspoleo(state.cells[i - 1].pres, state.cells[i - 1].temp);
        rhoWF = fluF.MasEspAgua(state.cells[i - 1].pres, state.cells[i - 1].temp);
        titF = (1 - fwF) * rhoOF / ((1 - fwF) * rhoOF + fwF * rhoWF);
        trF = 0.;
    }

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
    double tL = 0;
    double tH = 70.;

    double bo;
    double ba;
    double rs;
    double titV = 0.;
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
        // neste caso, variaveis como RGO de separador, BSW, API, densidade de gas e outras nÃ£o muda, sao iguais
        // aos valores da celula i-1
        double hol = 1 - state.cells[i - 1].alf;
        double bet = state.cells[i - 1].bet;
        double bsw = state.cells[i - 1].FW;
        double rhog = state.cells[i - 1].flui.MasEspGas(state.cells[i - 1].pres, state.cells[i - 1].temp);
        double rhogST = state.cells[i - 1].flui.Deng * 1.225;
        // o valor de volume de leve Ã© atualizado neste ponto,
        // seguindo o equacionamento mostrado em relatorio, nÃ£o Ã© relevante para o permanente, mas
        // deve ser calculado, pois Ã© utilizado como entrada no transiente
        state.cells[i - 1].VolLeveST = (((1 - hol) * rhog / rhogST) + hol * (1 - bet) * (1. - bsw) * rs / bo);
        if (state.cells[i - 1].VolLeveST < 1e-15)
            state.cells[i - 1].VolLeveST = 0.;

        state.cells[i].flui.BSW = state.cells[i - 1].flui.BSW;

        double fwV = state.cells[i - 1].FW;
        double rhoOV = state.cells[i - 1].flui.MasEspoleo(state.cells[i - 1].pres, state.cells[i - 1].temp);
        double rhoWV = state.cells[i - 1].flui.MasEspAgua(state.cells[i - 1].pres, state.cells[i - 1].temp);
        titV = (1 - fwV) * rhoOV / ((1 - fwV) * rhoOV + fwV * rhoWV);

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
            double boI;
            double baI;
            if (state.cells[i].flui.RGO < 1e6)
                boI = state.cells[i].flui.BOFunc(state.cells[i].presaux, tmed);
            else
                boI = 1.;
            baI = state.cells[i].flui.BAFunc(state.cells[i].presaux, tmed);
            double fwI = state.cells[i].flui.BSW * baI / (boI + baI * state.cells[i].flui.BSW - state.cells[i].flui.BSW * boI);
            double qo = state.cells[i - 1].QL * (1 - state.cells[i - 1].FW) * (1 - state.cells[i - 1].bet) * boI / bo;
            double qw;
            if (fwI < (1 - (*state.globals).localtiny * 1e-5))
                qw = fwI * qo / (1 - fwI);
            else
                qw = state.cells[i - 1].QL * (1 - state.cells[i - 1].bet);
            double qc;
            qc = state.cells[i - 1].QL * (state.cells[i - 1].bet) * state.cells[i - 1].fluicol.MasEspFlu(state.cells[i - 1].pres, state.cells[i - 1].temp) / state.cells[i].fluicol.MasEspFlu(state.cells[i].presaux, tmed);

            ////////////////////////esperar//////////////////////////////////////
            if (fabs(qo + qw + qc) > 1e-15)
                state.cells[i].bet = fabs(qc) / (fabs(qo) + fabs(qw) + fabs(qc));
            else
                state.cells[i].bet = state.cells[i - 1].bet;
        } else {
            mHidro = state.cells[i].MC - mComp;
        }
    } else {
        double fwV = state.cells[i - 1].FW;
        double rhoOV = state.cells[i - 1].flui.MasEspoleo(state.cells[i - 1].pres, state.cells[i - 1].temp);
        double rhoWV = state.cells[i - 1].flui.MasEspAgua(state.cells[i - 1].pres, state.cells[i - 1].temp);
        titV = (1 - fwV) * rhoOV / ((1 - fwV) * rhoOV + fwV * rhoWV);
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

        double api = state.cells[i - 1].flui.API;
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
            boinjl = bo;
            bainjl = ba;
            fwinjl = state.cells[i - 1].FW;
        }
        // vazao de oleo sytandard antes da fonte:
        double qostd1;
        if (state.cells[i - 1].flui.dStockTankVaporMassFraction < 1. - 1e-15)
            qostd1 = state.cells[i - 1].Mliqini * (1. - state.cells[i - 1].FW) * (1. - state.cells[i - 1].bet) / (bo * rholis);
        else
            qostd1 = 0.;

        // vazao de oleo standard da fonte:
        double qostd2;
        if (fluF.dStockTankVaporMassFraction < 1. - 1e-15)
            qostd2 = state.cells[i - 1].fontemassLR * (1. - fwinjl) / (boinjl * rholisF);
        else
            qostd2 = 0.;

        double hol = 1 - state.cells[i - 1].alf;
        double bet = state.cells[i - 1].bet;
        double bsw = state.cells[i - 1].FW;
        double rhog = state.cells[i - 1].flui.MasEspGas(state.cells[i - 1].pres, state.cells[i - 1].temp);
        double rhogST = state.cells[i - 1].flui.Deng * 1.225;
        // calculo do volume de leve na celula i, seguindo o equacionamento mostrado em relatorio,
        // nÃ£o Ã© relevante para o permanente, mas deve ser calculado,
        // pois Ã© utilizado como entrada no transiente
        state.cells[i - 1].VolLeveST = (((1 - hol) * rhog / rhogST) + hol * (1 - bet) * (1. - bsw) * rs / bo);
        if (state.cells[i - 1].VolLeveST < 1e-15)
            state.cells[i - 1].VolLeveST = 0.;
        if (state.input.trackRGO == -1) {
            // esta chave nÃ£o Ã© utilizada, Ã© mantida aqui como reserva, atualmente este
            // calculo nÃ£o esta funcionando, lembrando que arq.trackRGO assume apenas 2 valores, 0 ou 1
            if (hol > (*state.globals).localtiny && bet < (1. - (*state.globals).localtiny) && bsw < (1. - (*state.globals).localtiny) && mudaRGO == 1)
                state.cells[i].flui.RGO = state.cells[i - 1].VolLeveST * bo / (hol * (1 - bet) * (1 - bsw));
            if (state.cells[i].flui.RGO > (*state.globals).RGOMax && mudaRGO == 1 && state.cells[i].flui.RGO < 1e6)
                state.cells[i].flui.RGO = (*state.globals).RGOMax;
        }

        double qw1;
        if ((1. - state.cells[i - 1].flui.BSW) > 0)
            qw1 = qostd1 * state.cells[i - 1].flui.BSW / (1. - state.cells[i - 1].flui.BSW);
        else
            qw1 = state.cells[i - 1].Mliqini * state.cells[i - 1].FW * (1. - state.cells[i - 1].bet) / (1000. * state.cells[i - 1].flui.Denag);
        double qw2;
        if ((1. - fluF.BSW) > 0)
            qw2 = qostd2 * fluF.BSW / (1. - fluF.BSW);
        else
            qw2 = state.cells[i - 1].fontemassLR * fwinjl / (1000. * fluF.Denag);

        if (fabs(qw1 + qw2 + qostd1 + qostd2) > 1e-15 && fabs(qw2 + qostd2) > 1e-15 && state.cells[i - 1].fontemassLR > 1e-15)
            state.cells[i].flui.BSW = (qw1 + qw2) / (qw1 + qw2 + qostd1 + qostd2);
        else
            state.cells[i].flui.BSW = state.cells[i - 1].flui.BSW;

        if (fabs(qw1 + qw2) > 1e-15 && state.cells[i - 1].fontemassLR > 1e-15)
            state.cells[i].flui.Denag = (qw1 * state.cells[i - 1].flui.Denag +
                                    qw2 * fluF.Denag) /
                                   (qw1 + qw2);
        else
            state.cells[i].flui.Denag = state.cells[i - 1].flui.Denag;

        if (fabs(qostd1 + qostd2) > 1e-15 && fabs(qostd2) > 1e-15 && state.cells[i - 1].fontemassLR > 1e-15) { // vazao de liquido >0
            state.cells[i].flui.TempL = tL;
            state.cells[i].flui.TempH = tH;
            state.cells[i].flui.LVisL = (qostd1 * state.cells[i - 1].flui.VisOM(tL) + qostd2 * fluF.VisOM(tL)) / (qostd1 + qostd2);
            state.cells[i].flui.LVisH = (qostd1 * state.cells[i - 1].flui.VisOM(tH) + qostd2 * fluF.VisOM(tH)) / (qostd1 + qostd2);
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
            double boN = state.cells[i].flui.BOFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
            double baN = state.cells[i].flui.BAFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
            double fwN = state.cells[i].flui.BSW * baN /
                         (boN + baN * state.cells[i].flui.BSW - state.cells[i].flui.BSW * boN);
            double boI = state.cells[i].flui.BOFunc(state.cells[i].presaux, tmed);
            double baI = state.cells[i].flui.BAFunc(state.cells[i].presaux, tmed);
            double fwI = state.cells[i].flui.BSW * baI /
                         (boI + baI * state.cells[i].flui.BSW - state.cells[i].flui.BSW * boI);
            double qlpF = state.cells[i - 1].fontemassLR / fluF.MasEspLiq(state.cells[i - 1].pres, state.cells[i - 1].temp);
            double qlcF = state.cells[i - 1].fontemassCR / state.cells[i].fluicol.MasEspFlu(state.cells[i - 1].pres, state.cells[i - 1].temp);
            double betN = (state.cells[i - 1].QL * state.cells[i - 1].bet + qlcF) /
                          (state.cells[i - 1].QL + qlpF + qlcF);
            double qo = (state.cells[i - 1].QL + qlpF + qlcF) * (1 - fwN) * (1 - betN) * boI / bo;
            double qw;
            if (fwI < (1 - (*state.globals).localtiny * 1e-5))
                qw = fwI * qo / (1 - fwI);
            else
                qw = state.cells[i - 1].QL * (1 - state.cells[i - 1].bet) + qlpF;
            double qc = (state.cells[i - 1].QL * (state.cells[i - 1].bet) + qlcF) * state.cells[i - 1].fluicol.MasEspFlu(state.cells[i - 1].pres, state.cells[i - 1].temp) / state.cells[i].fluicol.MasEspFlu(state.cells[i].presaux, tmed);
            if (fabs(qo + qw + qc) > 1e-15)
                state.cells[i].bet = fabs(qc) / (fabs(qo) + fabs(qw) + fabs(qc));
            else
                state.cells[i].bet = state.cells[i - 1].bet;
        } else {
            mHidro = state.cells[i].MC - mComp;
        }
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
        if (state.cells[i].flui.iIER != 0) {
            int para;
            para = 0;
        }
    }
    double titulo = state.cells[i].flui.dVaporMassFraction;

    double betI;
    double fwI;
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
        rs = state.cells[i].flui.RS(pmed, tmed);
        bo = state.cells[i].flui.BOFunc(pmed, tmed, rs);
        rs = rs * 6.29 / 35.31467;
        ba = state.cells[i].flui.BAFunc(pmed, tmed);
        double rhogstd = state.cells[i].flui.Deng * 1.225;
        double rhololeostd = 1000. * 141.5 / (131.5 + state.cells[i].flui.API);
        double rhoa = state.cells[i].flui.Denag * 1000.;

        if (titulo < 1. - 1e-15)
            state.cells[i].FW = state.cells[i].flui.BSW * ba / (bo + ba * state.cells[i].flui.BSW - state.cells[i].flui.BSW * bo);
        else
            state.cells[i].FW = 1.;
        double ro = state.cells[i].flui.MasEspoleo(pmed, tmed);
        double ra = state.cells[i].flui.MasEspAgua(pmed, tmed);
        double rc = state.cells[i].fluicol.MasEspFlu(pmed, tmed);
        double denom = rc * state.cells[i].bet + (1 - state.cells[i].bet) * state.cells[i].FW * ra;
        if (titulo < 1. - 1e-15) {
            denom += (1 - state.cells[i].bet) * (1. - state.cells[i].FW) * ro / (1. - titulo);
            state.cells[i].QL = state.cells[i].MC / denom;

            state.cells[i].Mliqini = state.cells[i].QL * (state.cells[i].bet * rc + (1 - state.cells[i].bet) * (state.cells[i].FW * ra +
                                                                                            (1. - state.cells[i].FW) * ro));
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
        double ro = state.cells[i].flui.MasEspoleo(pmed, tmed);
        double rc = state.cells[i].fluicol.MasEspFlu(pmed, tmed);
        double qo = mHidro * (1. - titulo) / ro;
        double qc = mComp / rc;
        if (fabs(qo + qc) > 1e-15)
            state.cells[i].bet = qc / (qo + qc);
        else
            state.cells[i].bet = state.cells[i - 1].bet;
        double bet = state.cells[i].bet;
        double denom;
        if (titulo < 1. - 1e-15) {
            denom = rc * bet + (1. - bet) * ro;
            state.cells[i].QL = (mHidro * (1. - titulo) + mComp) / denom;
            state.cells[i].Mliqini = state.cells[i].QL * (state.cells[i].bet * rc + (1 - state.cells[i].bet) * ro);
        } else if (bet > (1 - 1e-15)) {
            state.cells[i].Mliqini = state.cells[i - 1].QL * state.cells[i - 1].betI * state.cells[i - 1].rcCi + state.cells[i - 1].fontemassCR;
            state.cells[i].QL = state.cells[i].Mliqini / (rc * bet + (1. - bet) * ro);
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
    if (fabs(state.cells[i].QL) > 1e-15 && fabs(state.cells[i].MComp) > 1e-15) {
        double dxmed = 0.5 * (state.cells[i - 1].dx + state.cells[i].dx);
        double hol0 = 1. - state.cells[i - 1].alf;
        double hol1 = 1. - state.cells[i].alf;
        state.cells[i].fluicol.TR = (state.cells[i - 1].fluicol.TR * state.cells[i - 1].MComp + state.cells[i - 1].fontemassCR * trF) / state.cells[i].MComp;
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

void applySteadyMassWithoutSource(const SteadyStateState &state, int i, int mudaRGO, double bo, double rs, double tmed, double &boI, double &baI, double &fwI) {
    // neste caso, variaveis como RGO de separador, BSW, API, densidade de gas e outras nÃƒÂ£o muda, sao iguais
    // aos valores da celula i-1
    double hol = 1 - state.cells[i - 1].alf;
    double bet = state.cells[i - 1].bet;
    double bsw = state.cells[i - 1].FW;
    double rhog = state.cells[i - 1].flui.MasEspGas(state.cells[i - 1].pres, state.cells[i - 1].temp);
    double rhogST = state.cells[i - 1].flui.Deng * 1.225;
    // o valor de volume de leve ÃƒÂ© atualizado neste ponto,
    // seguindo o equacionamento mostrado em relatorio, nÃƒÂ£o ÃƒÂ© relevante para o permanente, mas
    // deve ser calculado, pois ÃƒÂ© utilizado como entrada no transiente
    state.cells[i - 1].VolLeveST = (((1 - hol) * rhog / rhogST) + hol * (1 - bet) * (1. - bsw) * rs / bo);
    if (state.cells[i - 1].VolLeveST < 1e-15)
        state.cells[i - 1].VolLeveST = 0.;

    if (state.input.trackRGO == -1) {
        // esta chave nÃƒÂ£o ÃƒÂ© utilizada, ÃƒÂ© mantida aqui como reserva, atualmente este
        // calculo nÃƒÂ£o esta funcionando, lembrando que arq.trackRGO assume apenas 2 valores, 0 ou 1
        if (hol > (*state.globals).localtiny && bet < (1. - (*state.globals).localtiny) && bsw < (1. - (*state.globals).localtiny))
            state.cells[i].flui.RGO = state.cells[i - 1].VolLeveST * bo / (hol * (1 - bet) * (1 - bsw));
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
            boI = state.cells[i].flui.BOFunc(state.cells[i].presaux, tmed);
        else
            boI = 1.;
        baI = state.cells[i].flui.BAFunc(state.cells[i].presaux, tmed);
        fwI = state.cells[i].flui.BSW * baI / (boI + baI * state.cells[i].flui.BSW - state.cells[i].flui.BSW * boI);
        double qo = state.cells[i - 1].QL * (1 - state.cells[i - 1].FW) * (1 - state.cells[i - 1].bet) * boI / bo;
        double qw;
        if (fwI < (1 - (*state.globals).localtiny))
            qw = fwI * qo / (1 - fwI);
        else
            qw = state.cells[i - 1].QL * (1 - state.cells[i - 1].bet);
        double qc;
        if (state.cells[i].flui.RGO < 1e7)
            qc = state.cells[i - 1].QL * (state.cells[i - 1].bet) * state.cells[i - 1].fluicol.MasEspFlu(state.cells[i - 1].pres, state.cells[i - 1].temp) / state.cells[i].fluicol.MasEspFlu(state.cells[i].presaux, tmed);
        else
            qc = 0.;
        ////////////////////////esperar//////////////////////////////////////
        if ((fabs(qo) + fabs(qw) + fabs(qc)) > 0)
            state.cells[i].bet = fabs(qc) / (fabs(qo) + fabs(qw) + fabs(qc));
        else
            state.cells[i].bet = 0.;
    } else
        state.cells[i].bet = state.cells[i - 1].bet;
}

void applySteadyMassDryGasInjection(const SteadyStateState &state, int i, int mudaRGO, double tL, double tH, double bo, double ba, double rs, double tmed, double &trF, double &boI, double &baI, double &fwI) {
    double api = state.cells[i - 1].flui.API;
    double rhololeo = 1000. * 141.5 / (131.5 + api);
    double rholis = state.cells[i - 1].flui.MasEspLiq(state.cells[i - 1].pres, state.cells[i - 1].temp);
    fwI = 0.;
    // vazao de oleo standard
    double qostd;
    if (state.cells[i - 1].flui.dStockTankVaporMassFraction < 1. - 1e-15)
        qostd = state.cells[i - 1].Mliqini * (1. - state.cells[i - 1].FW) * (1. - state.cells[i - 1].bet) / (bo * rholis);
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

    double hol = 1 - state.cells[i - 1].alf;
    double bet = state.cells[i - 1].bet;
    double bsw = state.cells[i - 1].FW;
    double rhog = state.cells[i - 1].flui.MasEspGas(state.cells[i - 1].pres, state.cells[i - 1].temp);
    double rhogST = state.cells[i - 1].flui.Deng * 1.225;
    // calculo do volume de leve na celula i, seguindo o equacionamento mostrado em relatorio,
    // nÃƒÂ£o ÃƒÂ© relevante para o permanente, mas deve ser calculado,
    // pois ÃƒÂ© utilizado como entrada no transiente
    state.cells[i - 1].VolLeveST = (((1 - hol) * rhog / rhogST) + hol * (1 - bet) * (1. - bsw) * rs / bo);
    if (state.cells[i - 1].VolLeveST < 1e-15)
        state.cells[i - 1].VolLeveST = 0.;
    if (state.input.trackRGO == -1) {
        // esta chave nÃƒÂ£o ÃƒÂ© utilizada, ÃƒÂ© mantida aqui como reserva, atualmente este
        // calculo nÃƒÂ£o esta funcionando, lembrando que arq.trackRGO assume apenas 2 valores, 0 ou 1
        if (hol > (*state.globals).localtiny && bet < (1. - (*state.globals).localtiny) && bsw < (1. - (*state.globals).localtiny) && mudaRGO == 1)
            state.cells[i].flui.RGO = state.cells[i - 1].VolLeveST * bo / (hol * (1 - bet) * (1 - bsw));
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
        boI = state.cells[i].flui.BOFunc(state.cells[i].presaux, tmed);
        baI = state.cells[i].flui.BAFunc(state.cells[i].presaux, tmed);
        fwI = state.cells[i].flui.BSW * baI / (boI + baI * state.cells[i].flui.BSW - state.cells[i].flui.BSW * boI);
        double qo = state.cells[i - 1].QL * (1 - state.cells[i - 1].FW) * (1 - state.cells[i - 1].bet) * boI / bo;
        double qw;
        if (fwI < (1 - (*state.globals).localtiny))
            qw = fwI * qo / (1 - fwI);
        else
            qw = state.cells[i - 1].QL * (1 - state.cells[i - 1].bet);
        double qc = state.cells[i - 1].QL * (state.cells[i - 1].bet) * state.cells[i - 1].fluicol.MasEspFlu(state.cells[i - 1].pres, state.cells[i - 1].temp) / state.cells[i].fluicol.MasEspFlu(state.cells[i].presaux, tmed);
        state.cells[i].bet = fabs(qc) / (fabs(qo) + fabs(qw) + fabs(qc));
    } else
        state.cells[i].bet = state.cells[i - 1].bet;
}

void applySteadyMassWetGasInjection(const SteadyStateState &state, int i, int mudaRGO, double tL, double tH, double bo, double ba, double rs, double tmed, double &trF, double &boI, double &baI, double &fwI) {
    double api = state.cells[i - 1].flui.API;
    double rhololeo = 1000. * 141.5 / (131.5 + api);
    double rholis = state.cells[i - 1].flui.MasEspLiq(state.cells[i - 1].pres, state.cells[i - 1].temp);
    double rholisF;
    double boinjl;
    double bainjl;
    double fwinjl;
    trF = state.cells[i - 1].acsr.injg.fluidocol.TR;
    if (fabs(state.cells[i - 1].acsr.injg.QGas) > 1e-15) {
        rholisF = state.cells[i - 1].acsr.injg.FluidoPro.MasEspLiq(state.cells[i - 1].pres, state.cells[i - 1].temp);
        boinjl = state.cells[i - 1].acsr.injg.FluidoPro.BOFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        bainjl = state.cells[i - 1].acsr.injg.FluidoPro.BAFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        fwinjl = state.cells[i - 1].acsr.injg.FluidoPro.BSW * bainjl / (boinjl + bainjl * state.cells[i - 1].acsr.injg.FluidoPro.BSW - state.cells[i - 1].acsr.injg.FluidoPro.BSW * boinjl);
    } else {
        rholisF = rholis;
        boinjl = bo;
        bainjl = ba;
        fwinjl = state.cells[i - 1].FW;
    }
    // vazao de oleo standard
    double qostd = 0.;
    if (state.cells[i - 1].flui.dStockTankVaporMassFraction < 1. - 1e-15)
        qostd = state.cells[i - 1].Mliqini * (1. - state.cells[i - 1].FW) *
                (1. - state.cells[i - 1].bet) / (bo * rholis);
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

    double hol = 1 - state.cells[i - 1].alf;
    double bet = state.cells[i - 1].bet;
    double bsw = state.cells[i - 1].FW;
    double rhog = state.cells[i - 1].flui.MasEspGas(state.cells[i - 1].pres, state.cells[i - 1].temp);
    double rhogST = state.cells[i - 1].flui.Deng * 1.225;
    // calculo do volume de leve na celula i, seguindo o equacionamento mostrado em relatorio,
    // nÃƒÂ£o ÃƒÂ© relevante para o permanente, mas deve ser calculado,
    // pois ÃƒÂ© utilizado como entrada no transiente
    state.cells[i - 1].VolLeveST = (((1 - hol) * rhog / rhogST) + hol * (1 - bet) * (1. - bsw) * rs / bo);
    if (state.cells[i - 1].VolLeveST < 1e-15)
        state.cells[i - 1].VolLeveST = 0.;
    if (state.input.trackRGO == -1) {
        // esta chave nÃƒÂ£o ÃƒÂ© utilizada, ÃƒÂ© mantida aqui como reserva, atualmente este
        // calculo nÃƒÂ£o esta funcionando, lembrando que arq.trackRGO assume apenas 2 valores, 0 ou 1
        if (hol > (*state.globals).localtiny && bet < (1. - (*state.globals).localtiny) && bsw < (1. - (*state.globals).localtiny) && mudaRGO == 1)
            state.cells[i].flui.RGO = state.cells[i - 1].VolLeveST * bo / (hol * (1 - bet) * (1 - bsw));
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
        double qw1;
        if ((1. - state.cells[i - 1].flui.BSW) > 0)
            qw1 = qostd * state.cells[i - 1].flui.BSW / (1. - state.cells[i - 1].flui.BSW);
        else
            qw1 = state.cells[i - 1].Mliqini * state.cells[i - 1].FW * (1. - state.cells[i - 1].bet) / (1000. * state.cells[i - 1].flui.Denag);
        double qw2;
        if ((1. - state.cells[i - 1].acsr.injg.FluidoPro.BSW) > 0)
            qw2 = qostd2 * state.cells[i - 1].acsr.injg.FluidoPro.BSW / (1. - state.cells[i - 1].acsr.injg.FluidoPro.BSW);
        else
            qw2 = state.cells[i - 1].fontemassLR * fwinjl / (1000. * state.cells[i - 1].acsr.injg.FluidoPro.Denag);

        if (fabs(qw1 + qw2 + qostd + qostd2) > 1e-15 && fabs(state.cells[i - 1].acsr.injg.QGas) > 1e-15)
            state.cells[i].flui.BSW = (qw1 + qw2) / (qw1 + qw2 + qostd + qostd2);
        else
            state.cells[i].flui.BSW = state.cells[i - 1].flui.BSW;
        if (fabs(qw1 + qw2) > 1e-15)
            state.cells[i].flui.Denag = (qw1 * state.cells[i - 1].flui.Denag +
                                    qw2 * state.cells[i - 1].acsr.injg.FluidoPro.Denag) /
                                   (qw1 + qw2);
        else
            state.cells[i].flui.Denag = state.cells[i - 1].flui.Denag;
        // reavaliacao dos pares necessarios para o modelo de viscosidade de oleo morto
        // ASTM. Isto ÃƒÂ© feito tendo sempre os mesmos valores de temperatura maxima e
        // e minima para a construcao dos pares
        if (fabs(qostd + qostd2) > 1e-15 && state.input.flashCompleto == 0 && fabs(state.cells[i - 1].acsr.injg.QGas) > 1e-15) { // vazao de liquido >0
            state.cells[i].flui.TempL = tL;
            state.cells[i].flui.TempH = tH;
            state.cells[i].flui.LVisL = (qostd * state.cells[i - 1].flui.VisOM(tL) + qostd2 * state.cells[i - 1].acsr.injg.FluidoPro.VisOM(tL)) / (qostd + qostd2);
            state.cells[i].flui.LVisH = (qostd * state.cells[i - 1].flui.VisOM(tH) + qostd2 * state.cells[i - 1].acsr.injg.FluidoPro.VisOM(tH)) / (qostd + qostd2);
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
        double boN = state.cells[i].flui.BOFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        double baN = state.cells[i].flui.BAFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        double fwN = state.cells[i].flui.BSW * baN / (boN + baN * state.cells[i].flui.BSW - state.cells[i].flui.BSW * boN);
        boI = state.cells[i].flui.BOFunc(state.cells[i].presaux, tmed);
        baI = state.cells[i].flui.BAFunc(state.cells[i].presaux, tmed);
        fwI = state.cells[i].flui.BSW * baI / (boI + baI * state.cells[i].flui.BSW - state.cells[i].flui.BSW * boI);
        double qlpF = state.cells[i - 1].fontemassLR / state.cells[i - 1].acsr.injg.FluidoPro.MasEspLiq(state.cells[i - 1].pres, state.cells[i - 1].temp);
        double qlcF = state.cells[i - 1].fontemassCR / state.cells[i - 1].acsr.injg.fluidocol.MasEspFlu(state.cells[i - 1].pres, state.cells[i - 1].temp);

        double betN = (state.cells[i - 1].QL * state.cells[i - 1].bet + qlcF) / (state.cells[i - 1].QL + qlpF + qlcF);
        double qo = (state.cells[i - 1].QL + qlpF + qlcF) * (1 - fwN) * (1 - betN) * boI / bo;

        double qw;
        if (fwI < (1 - (*state.globals).localtiny))
            qw = fwI * qo / (1 - fwI);
        else
            qw = state.cells[i - 1].QL * (1 - state.cells[i - 1].bet) + qlpF;
        //* celula[i - 1].fluicol.MasEspFlu(celula[i - 1].pres, celula[i - 1].temp)

        double qc = (state.cells[i - 1].QL * (state.cells[i - 1].bet) + qlcF) * state.cells[i - 1].fluicol.MasEspFlu(state.cells[i - 1].pres, state.cells[i - 1].temp) / state.cells[i].fluicol.MasEspFlu(state.cells[i].presaux, tmed);

        state.cells[i].bet = fabs(qc) / (fabs(qo) + fabs(qw) + fabs(qc));
    } else
        state.cells[i].bet = state.cells[i - 1].bet;
}

void applySteadyMassLiquidInjection(const SteadyStateState &state, int i, int mudaRGO, double tL, double tH, double bo, double ba, double rs, double tmed, double &trF, double &boI, double &baI, double &fwI) {
    double api = state.cells[i - 1].flui.API;
    double rhololeo = 1000. * 141.5 / (131.5 + api);
    double rhololeoF = 1000. * 141.5 / (131.5 + state.cells[i - 1].acsr.injl.FluidoPro.API);
    double rholis = state.cells[i - 1].flui.MasEspLiq(state.cells[i - 1].pres, state.cells[i - 1].temp);
    double rholisF;
    double boinjl;
    double bainjl;
    double fwinjl;
    trF = state.cells[i - 1].acsr.injl.fluidocol.TR;
    if (fabs(state.cells[i - 1].acsr.injl.QLiq) > 1e-15) {
        rholisF = state.cells[i - 1].acsr.injl.FluidoPro.MasEspLiq(state.cells[i - 1].pres, state.cells[i - 1].temp);
        boinjl = state.cells[i - 1].acsr.injl.FluidoPro.BOFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        bainjl = state.cells[i - 1].acsr.injl.FluidoPro.BAFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        fwinjl = state.cells[i - 1].acsr.injl.FluidoPro.BSW * bainjl / (boinjl + bainjl * state.cells[i - 1].acsr.injl.FluidoPro.BSW - state.cells[i - 1].acsr.injl.FluidoPro.BSW * boinjl);
    } else {
        rholisF = rholis;
        boinjl = bo;
        bainjl = ba;
        fwinjl = state.cells[i - 1].FW;
    }
    // vazao de oleo sytandard antes da fonte:
    double qostd1 = 0.;
    if (state.cells[i - 1].flui.dStockTankVaporMassFraction < 1. - 1e-15)
        qostd1 = state.cells[i - 1].Mliqini * (1. - state.cells[i - 1].FW) * (1. - state.cells[i - 1].bet) / (bo * rholis);
    // vazao de oleo standard da fonte:
    double qostd2 = 0.;
    if (state.cells[i - 1].acsr.injl.FluidoPro.dStockTankVaporMassFraction < 1. - 1e-15)
        qostd2 = state.cells[i - 1].fontemassLR * (1. - fwinjl) / (boinjl * rholisF);
    // vazao total de gas nas condicoes standard fruto da soma do gas transportado
    // da fronteira a esquerda
    // da celula e do gas associado da fonte de liquido
    double qgstd = qostd1 * state.cells[i - 1].flui.RGO + qostd2 * state.cells[i - 1].acsr.injl.FluidoPro.RGO;
    double hol = 1 - state.cells[i - 1].alf;
    double bet = state.cells[i - 1].bet;
    double bsw = state.cells[i - 1].FW;
    double rhog = state.cells[i - 1].flui.MasEspGas(state.cells[i - 1].pres, state.cells[i - 1].temp);
    double rhogST = state.cells[i - 1].flui.Deng * 1.225;
    // calculo do volume de leve na celula i, seguindo o equacionamento mostrado em relatorio,
    // nÃƒÂ£o ÃƒÂ© relevante para o permanente, mas deve ser calculado,
    // pois ÃƒÂ© utilizado como entrada no transiente
    state.cells[i - 1].VolLeveST = (((1 - hol) * rhog / rhogST) + hol * (1 - bet) * (1. - bsw) * rs / bo);
    if (state.cells[i - 1].VolLeveST < 1e-15)
        state.cells[i - 1].VolLeveST = 0.;

    if (fabs(state.cells[i - 1].acsr.injl.QLiq) > 1e-15) {
        if (state.input.trackRGO == -1) {
            // esta chave nÃƒÂ£o ÃƒÂ© utilizada, ÃƒÂ© mantida aqui como reserva, atualmente este
            // calculo nÃƒÂ£o esta funcionando, lembrando que arq.trackRGO assume apenas 2 valores, 0 ou 1
            if (hol > (*state.globals).localtiny && bet < (1. - (*state.globals).localtiny) && bsw < (1. - (*state.globals).localtiny) && mudaRGO == 1)
                state.cells[i].flui.RGO = state.cells[i - 1].VolLeveST * bo / (hol * (1 - bet) * (1 - bsw));
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
        double qw1;
        if ((1. - state.cells[i - 1].flui.BSW) > 0)
            qw1 = qostd1 * state.cells[i - 1].flui.BSW / (1. - state.cells[i - 1].flui.BSW);
        else
            qw1 = state.cells[i - 1].Mliqini * state.cells[i - 1].FW * (1. - state.cells[i - 1].bet) / (1000. * state.cells[i - 1].flui.Denag);
        double qw2;
        if ((1. - state.cells[i - 1].acsr.injl.FluidoPro.BSW) > 0)
            qw2 = qostd2 * state.cells[i - 1].acsr.injl.FluidoPro.BSW / (1. - state.cells[i - 1].acsr.injl.FluidoPro.BSW);
        else
            qw2 = state.cells[i - 1].fontemassLR * fwinjl / (1000. * state.cells[i - 1].acsr.injl.FluidoPro.Denag);

        if (fabs(qw1 + qw2 + qostd1 + qostd2) > 1e-15 &&
            fabs(state.cells[i - 1].acsr.injl.QLiq) > 1e-15)
            state.cells[i].flui.BSW = (qw1 + qw2) / (qw1 + qw2 + qostd1 + qostd2);
        else
            state.cells[i].flui.BSW = state.cells[i - 1].flui.BSW;
        if (fabs(qw1 + qw2) > 1e-15)
            state.cells[i].flui.Denag = (qw1 * state.cells[i - 1].flui.Denag +
                                    qw2 * state.cells[i - 1].acsr.injl.FluidoPro.Denag) /
                                   (qw1 + qw2);
        else
            state.cells[i].flui.Denag = state.cells[i - 1].flui.Denag;
        // reavaliacao dos pares necessarios para o modelo de viscosidade de oleo morto
        // ASTM. Isto ÃƒÂ© feito tendo sempre os mesmos valores de temperatura maxima e
        // e minima para a construcao dos pares
        if (fabs(qostd1 + qostd2) > 1e-15 && state.input.flashCompleto == 0 &&
            fabs(state.cells[i - 1].acsr.injl.QLiq) > 1e-15) { // vazao de liquido >0
            state.cells[i].flui.TempL = tL;
            state.cells[i].flui.TempH = tH;
            state.cells[i].flui.LVisL = (qostd1 * state.cells[i - 1].flui.VisOM(tL) + qostd2 * state.cells[i - 1].acsr.injl.FluidoPro.VisOM(tL)) / (qostd1 + qostd2);
            state.cells[i].flui.LVisH = (qostd1 * state.cells[i - 1].flui.VisOM(tH) + qostd2 * state.cells[i - 1].acsr.injl.FluidoPro.VisOM(tH)) / (qostd1 + qostd2);
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
        double boN = state.cells[i].flui.BOFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        double baN = state.cells[i].flui.BAFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        double fwN = state.cells[i].flui.BSW * baN / (boN + baN * state.cells[i].flui.BSW - state.cells[i].flui.BSW * boN);
        boI = state.cells[i].flui.BOFunc(state.cells[i].presaux, tmed);
        baI = state.cells[i].flui.BAFunc(state.cells[i].presaux, tmed);
        fwI = state.cells[i].flui.BSW * baI / (boI + baI * state.cells[i].flui.BSW - state.cells[i].flui.BSW * boI);
        double qlpF = state.cells[i - 1].fontemassLR / state.cells[i - 1].acsr.injl.FluidoPro.MasEspLiq(state.cells[i - 1].pres, state.cells[i - 1].temp);
        double qlcF = state.cells[i - 1].fontemassCR / state.cells[i - 1].acsr.injl.fluidocol.MasEspFlu(state.cells[i - 1].pres, state.cells[i - 1].temp);
        double betN = (state.cells[i - 1].QL * state.cells[i - 1].bet + qlcF) / (state.cells[i - 1].QL + qlpF + qlcF);
        double qo = (state.cells[i - 1].QL + qlpF + qlcF) * (1 - fwN) * (1 - betN) * boI / bo;
        double qw;
        if (fwI < (1 - (*state.globals).localtiny))
            qw = fwI * qo / (1 - fwI);
        else
            qw = state.cells[i - 1].QL * (1 - state.cells[i - 1].bet) + qlpF;
        double qc = (state.cells[i - 1].QL * (state.cells[i - 1].bet) + qlcF) * state.cells[i - 1].fluicol.MasEspFlu(state.cells[i - 1].pres, state.cells[i - 1].temp) / state.cells[i].fluicol.MasEspFlu(state.cells[i].presaux, tmed);
        state.cells[i].bet = fabs(qc) / (fabs(qo) + fabs(qw) + fabs(qc));
    } else
        state.cells[i].bet = state.cells[i - 1].bet;
}

void applySteadyMassInflowPerformance(const SteadyStateState &state, int i, int mudaRGO, double tL, double tH, double bo, double ba, double rs, double tmed, double &trF, double &boI, double &baI, double &fwI) {
    double api = state.cells[i - 1].flui.API;
    double rhololeo = 1000. * 141.5 / (131.5 + api);
    double rhololeoF = 1000. * 141.5 / (131.5 + state.cells[i - 1].acsr.ipr.FluidoPro.API);
    double rholis = state.cells[i - 1].flui.MasEspLiq(state.cells[i - 1].pres, state.cells[i - 1].temp);
    double rholisF;
    double boipr;
    double baipr;
    double fwipr;
    trF = 0.;
    if (fabs(state.cells[i - 1].fontemassLR) > 1e-15) {
        rholisF = state.cells[i - 1].acsr.ipr.FluidoPro.MasEspLiq(state.cells[i - 1].pres, state.cells[i - 1].temp);
        boipr = state.cells[i - 1].acsr.ipr.FluidoPro.BOFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        baipr = state.cells[i - 1].acsr.ipr.FluidoPro.BAFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        fwipr = state.cells[i - 1].acsr.ipr.FluidoPro.BSW * baipr / (boipr + baipr * state.cells[i - 1].acsr.ipr.FluidoPro.BSW - state.cells[i - 1].acsr.ipr.FluidoPro.BSW * boipr);
    } else {
        rholisF = rholis;
        boipr = bo;
        baipr = ba;
        fwipr = state.cells[i - 1].FW;
    }
    // vazao de oleo standard antes da fonte:
    double qostd1 = 0.;
    if (state.cells[i - 1].flui.dStockTankVaporMassFraction < 1. - 1e-15)
        qostd1 = state.cells[i - 1].Mliqini * (1. - state.cells[i - 1].FW) * (1. - state.cells[i - 1].bet) / (bo * rholis);
    // vazao de oleo standard da fonte:
    double qostd2 = 0.;
    if (state.cells[i - 1].acsr.ipr.FluidoPro.dStockTankVaporMassFraction < 1. - 1e-15)
        qostd2 = state.cells[i - 1].fontemassLR * (1. - fwipr) / (boipr * rholisF);
    // vazao total de gas nas condicoes standard fruto da soma do gas transportado
    // da fronteira a esquerda
    // da celula e do gas associado da IPR
    double qgstd = qostd1 * state.cells[i - 1].flui.RGO + qostd2 * state.cells[i - 1].acsr.ipr.FluidoPro.RGO;
    double hol = 1 - state.cells[i - 1].alf;
    double bet = state.cells[i - 1].bet;
    double bsw = state.cells[i - 1].FW;
    double rhog = state.cells[i - 1].flui.MasEspGas(state.cells[i - 1].pres, state.cells[i - 1].temp);
    double rhogST = state.cells[i - 1].flui.Deng * 1.225;
    // calculo do volume de leve na celula i, seguindo o equacionamento mostrado em relatorio,
    // nÃƒÂ£o ÃƒÂ© relevante para o permanente, mas deve ser calculado,
    // pois ÃƒÂ© utilizado como entrada no transiente
    state.cells[i - 1].VolLeveST = (((1 - hol) * rhog / rhogST) + hol * (1 - bet) * (1. - bsw) * rs / bo);
    if (state.cells[i - 1].VolLeveST < 1e-15)
        state.cells[i - 1].VolLeveST = 0.;
    if (state.input.trackRGO == -1) {
        // esta chave nÃƒÂ£o ÃƒÂ© utilizada, ÃƒÂ© mantida aqui como reserva, atualmente este
        // calculo nÃƒÂ£o esta funcionando, lembrando que arq.trackRGO assume apenas 2 valores, 0 ou 1
        if ((state.cells[i - 1].fontemassLR) > 1e-15 && (state.cells[i - 1].fontemassGR) > 1e-15) {
            if (hol > (*state.globals).localtiny && bet < (1. - (*state.globals).localtiny) && bsw < (1. - (*state.globals).localtiny) &&
                mudaRGO == 1 && fabs(state.cells[i - 1].fontemassLR) > 1e-15)
                state.cells[i].flui.RGO = state.cells[i - 1].VolLeveST * bo / (hol * (1 - bet) * (1 - bsw));
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
        double qw1;
        if ((1. - state.cells[i - 1].flui.BSW) > 0)
            qw1 = qostd1 * state.cells[i - 1].flui.BSW / (1. - state.cells[i - 1].flui.BSW);
        else
            qw1 = state.cells[i - 1].Mliqini * state.cells[i - 1].FW * (1. - state.cells[i - 1].bet) / (1000. * state.cells[i - 1].flui.Denag);
        double qw2;
        if ((1. - state.cells[i - 1].acsr.ipr.FluidoPro.BSW) > 0)
            qw2 = qostd2 * state.cells[i - 1].acsr.ipr.FluidoPro.BSW / (1. - state.cells[i - 1].acsr.ipr.FluidoPro.BSW);
        else
            qw2 = state.cells[i - 1].fontemassLR * fwipr / (1000. * state.cells[i - 1].acsr.ipr.FluidoPro.Denag);

        if (fabs(qw1 + qw2 + qostd1 + qostd2) > 1e-15 &&
            (state.cells[i - 1].fontemassLR) > 1e-15)
            state.cells[i].flui.BSW = ((qw1 + qostd1) * state.cells[i - 1].flui.BSW + (qw2 + qostd2) * state.cells[i - 1].acsr.ipr.FluidoPro.BSW) / (qw1 + qw2 + qostd1 + qostd2);
        else
            state.cells[i].flui.BSW = state.cells[i - 1].flui.BSW;
        if ((qw1 + qw2) > 1e-15)
            state.cells[i].flui.Denag = (qw1 * state.cells[i - 1].flui.Denag +
                                    qw2 * state.cells[i - 1].acsr.ipr.FluidoPro.Denag) /
                                   (qw1 + qw2);
        else
            state.cells[i].flui.Denag = state.cells[i - 1].flui.Denag;
        if ((qostd1 + qostd2) > 1e-15 && state.input.flashCompleto == 0 && (state.cells[i - 1].fontemassLR) > 1e-15) {
            // reavaliacao dos pares necessarios para o modelo de viscosidade de oleo morto
            // ASTM. Isto ÃƒÂ© feito tendo sempre os mesmos valores de temperatura maxima e
            // e minima para a construcao dos pares
            state.cells[i].flui.TempL = tL;
            state.cells[i].flui.TempH = tH;
            state.cells[i].flui.LVisL = (qostd1 * state.cells[i - 1].flui.VisOM(tL) + qostd2 * state.cells[i - 1].acsr.ipr.FluidoPro.VisOM(tL)) / (qostd1 + qostd2);
            state.cells[i].flui.LVisH = (qostd1 * state.cells[i - 1].flui.VisOM(tH) + qostd2 * state.cells[i - 1].acsr.ipr.FluidoPro.VisOM(tH)) / (qostd1 + qostd2);
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
        double boN = state.cells[i].flui.BOFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        double baN = state.cells[i].flui.BAFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        double fwN = state.cells[i].flui.BSW * baN / (boN + baN * state.cells[i].flui.BSW - state.cells[i].flui.BSW * boN);
        boI = state.cells[i].flui.BOFunc(state.cells[i].presaux, tmed);
        baI = state.cells[i].flui.BAFunc(state.cells[i].presaux, tmed);
        fwI = state.cells[i].flui.BSW * baI / (boI + baI * state.cells[i].flui.BSW - state.cells[i].flui.BSW * boI);
        double qlpF = state.cells[i - 1].fontemassLR / state.cells[i - 1].acsr.ipr.FluidoPro.MasEspLiq(state.cells[i - 1].pres, state.cells[i - 1].temp);
        double qlcF = 0.;
        double betN = (state.cells[i - 1].QL * state.cells[i - 1].bet + qlcF) / (state.cells[i - 1].QL + qlpF + qlcF);
        double qo = (state.cells[i - 1].QL + qlpF + qlcF) * (1 - fwN) * (1 - betN) * boI / bo;
        double qw;
        if (fwI < (1 - (*state.globals).localtiny))
            qw = fwI * qo / (1 - fwI);
        else
            qw = state.cells[i - 1].QL * (1 - state.cells[i - 1].bet) + qlpF;
        double qc = (state.cells[i - 1].QL * (state.cells[i - 1].bet) + qlcF) * state.cells[i - 1].fluicol.MasEspFlu(state.cells[i - 1].pres, state.cells[i - 1].temp) / state.cells[i].fluicol.MasEspFlu(state.cells[i].presaux, tmed);
        state.cells[i].bet = fabs(qc) / (fabs(qo) + fabs(qw) + fabs(qc));
    } else
        state.cells[i].bet = state.cells[i - 1].bet;
}

void applySteadyMassMultipleSource(const SteadyStateState &state, int i, int mudaRGO, double tL, double tH, double bo, double ba, double rs, double tmed, double &trF, double &boI, double &baI, double &fwI) {
    double api = state.cells[i - 1].flui.API;
    double rhololeo = 1000. * 141.5 / (131.5 + api);
    double rhololeoF = 1000. * 141.5 / (131.5 + state.cells[i - 1].acsr.injm.FluidoPro.API);
    double rholis = state.cells[i - 1].flui.MasEspLiq(state.cells[i - 1].pres, state.cells[i - 1].temp);
    double rholisF;
    double boinjl;
    double bainjl;
    double fwinjl;
    trF = state.cells[i - 1].acsr.injm.fluidocol.TR;
    if (fabs(state.cells[i - 1].acsr.injm.MassG + state.cells[i - 1].acsr.injm.MassP) > 1e-15) {
        rholisF = state.cells[i - 1].acsr.injm.FluidoPro.MasEspLiq(state.cells[i - 1].pres, state.cells[i - 1].temp);
        boinjl = state.cells[i - 1].acsr.injm.FluidoPro.BOFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        bainjl = state.cells[i - 1].acsr.injm.FluidoPro.BAFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        fwinjl = state.cells[i - 1].acsr.injm.FluidoPro.BSW * bainjl / (boinjl + bainjl * state.cells[i - 1].acsr.injm.FluidoPro.BSW - state.cells[i - 1].acsr.injm.FluidoPro.BSW * boinjl);
    } else {
        rholisF = rholis;
        boinjl = bo;
        bainjl = ba;
        fwinjl = state.cells[i - 1].FW;
    }
    // vazao de oleo sytandard antes da fonte:
    double qostd1 = 0.;
    if (state.cells[i - 1].flui.dStockTankVaporMassFraction < 1. - 1e-15)
        qostd1 = state.cells[i - 1].Mliqini * (1. - state.cells[i - 1].FW) * (1. - state.cells[i - 1].bet) / (bo * rholis);
    // vazao de oleo standard da fonte:
    double qostd2 = 0.;
    if (state.cells[i - 1].acsr.injm.FluidoPro.dStockTankVaporMassFraction < 1. - 1e-15)
        qostd2 = state.cells[i - 1].fontemassLR * (1. - fwinjl) / (boinjl * rholisF);
    // vazao total de gas nas condicoes standard fruto da soma do gas transportado
    // da fronteira a esquerda
    // da celula e do gas associado da fonte de liquido
    double qgstd = qostd1 * state.cells[i - 1].flui.RGO + qostd2 * state.cells[i - 1].acsr.injm.FluidoPro.RGO;
    double hol = 1 - state.cells[i - 1].alf;
    double bet = state.cells[i - 1].bet;
    double bsw = state.cells[i - 1].FW;
    double rhog = state.cells[i - 1].flui.MasEspGas(state.cells[i - 1].pres, state.cells[i - 1].temp);
    double rhogST = state.cells[i - 1].flui.Deng * 1.225;
    // calculo do volume de leve na celula i, seguindo o equacionamento mostrado em relatorio,
    // nÃƒÂ£o ÃƒÂ© relevante para o permanente, mas deve ser calculado,
    // pois ÃƒÂ© utilizado como entrada no transiente
    state.cells[i - 1].VolLeveST = (((1 - hol) * rhog / rhogST) + hol * (1 - bet) * (1. - bsw) * rs / bo);
    if (state.cells[i - 1].VolLeveST < 1e-15)
        state.cells[i - 1].VolLeveST = 0.;

    if (state.input.trackRGO == -1) {
        // esta chave nÃƒÂ£o ÃƒÂ© utilizada, ÃƒÂ© mantida aqui como reserva, atualmente este
        // calculo nÃƒÂ£o esta funcionando, lembrando que arq.trackRGO assume apenas 2 valores, 0 ou 1
        if ((state.cells[i - 1].acsr.injm.MassG + state.cells[i - 1].acsr.injm.MassP) > 1e-15) {
            if (hol > (*state.globals).localtiny && bet < (1. - (*state.globals).localtiny) && bsw < (1. - (*state.globals).localtiny) && mudaRGO == 1)
                state.cells[i].flui.RGO = state.cells[i - 1].VolLeveST * bo / (hol * (1 - bet) * (1 - bsw));
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
        double qw1;
        if ((1. - state.cells[i - 1].flui.BSW) > 0)
            qw1 = qostd1 * state.cells[i - 1].flui.BSW / (1. - state.cells[i - 1].flui.BSW);
        else
            qw1 = state.cells[i - 1].Mliqini * state.cells[i - 1].FW * (1. - state.cells[i - 1].bet) / (1000. * state.cells[i - 1].flui.Denag);
        double qw2;
        if ((1. - state.cells[i - 1].acsr.injm.FluidoPro.BSW) > 0)
            qw2 = qostd2 * state.cells[i - 1].acsr.injm.FluidoPro.BSW / (1. - state.cells[i - 1].acsr.injm.FluidoPro.BSW);
        else
            qw2 = state.cells[i - 1].fontemassLR * fwinjl / (1000. * state.cells[i - 1].acsr.injm.FluidoPro.Denag);

        if (fabs(qw1 + qw2 + qostd1 + qostd2) > 1e-15 &&
            (state.cells[i - 1].acsr.injm.MassG + state.cells[i - 1].acsr.injm.MassP) > 1e-15)
            state.cells[i].flui.BSW = (qw1 + qw2) / (qw1 + qw2 + qostd1 + qostd2);
        else
            state.cells[i].flui.BSW = state.cells[i - 1].flui.BSW;
        if (fabs(qw1 + qw2) > 1e-15 &&
            (state.cells[i - 1].acsr.injm.MassG + state.cells[i - 1].acsr.injm.MassP) > 1e-15)
            state.cells[i].flui.Denag = (qw1 * state.cells[i - 1].flui.Denag +
                                    qw2 * state.cells[i - 1].acsr.injm.FluidoPro.Denag) /
                                   (qw1 + qw2);
        else
            state.cells[i].flui.Denag = state.cells[i - 1].flui.Denag;
        // reavaliacao dos pares necessarios para o modelo de viscosidade de oleo morto
        // ASTM. Isto ÃƒÂ© feito tendo sempre os mesmos valores de temperatura maxima e
        // e minima para a construcao dos pares
        if (fabs(qostd1 + qostd2) > 1e-15 && state.input.flashCompleto == 0 &&
            (state.cells[i - 1].acsr.injm.MassG + state.cells[i - 1].acsr.injm.MassP) > 1e-15) { // vazao de liquido >0
            state.cells[i].flui.TempL = tL;
            state.cells[i].flui.TempH = tH;
            state.cells[i].flui.LVisL = (qostd1 * state.cells[i - 1].flui.VisOM(tL) + qostd2 * state.cells[i - 1].acsr.injm.FluidoPro.VisOM(tL)) / (qostd1 + qostd2);
            state.cells[i].flui.LVisH = (qostd1 * state.cells[i - 1].flui.VisOM(tH) + qostd2 * state.cells[i - 1].acsr.injm.FluidoPro.VisOM(tH)) / (qostd1 + qostd2);
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
        double boN = state.cells[i].flui.BOFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        double baN = state.cells[i].flui.BAFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        double fwN = state.cells[i].flui.BSW * baN / (boN + baN * state.cells[i].flui.BSW - state.cells[i].flui.BSW * boN);
        boI = state.cells[i].flui.BOFunc(state.cells[i].presaux, tmed);
        baI = state.cells[i].flui.BAFunc(state.cells[i].presaux, tmed);
        fwI = state.cells[i].flui.BSW * baI / (boI + baI * state.cells[i].flui.BSW - state.cells[i].flui.BSW * boI);
        double qlpF = state.cells[i - 1].fontemassLR / state.cells[i - 1].acsr.injm.FluidoPro.MasEspLiq(state.cells[i - 1].pres, state.cells[i - 1].temp);
        double qlcF = state.cells[i - 1].fontemassCR / state.cells[i - 1].acsr.injm.fluidocol.MasEspFlu(state.cells[i - 1].pres, state.cells[i - 1].temp);
        double betN = (state.cells[i - 1].QL * state.cells[i - 1].bet + qlcF) / (state.cells[i - 1].QL + qlpF + qlcF);
        double qo = (state.cells[i - 1].QL + qlpF + qlcF) * (1 - fwN) * (1 - betN) * boI / bo;
        double qw;
        if (fwI < (1 - (*state.globals).localtiny))
            qw = fwI * qo / (1 - fwI);
        else
            qw = state.cells[i - 1].QL * (1 - state.cells[i - 1].bet) + qlpF;
        double qc = (state.cells[i - 1].QL * (state.cells[i - 1].bet) + qlcF) * state.cells[i - 1].fluicol.MasEspFlu(state.cells[i - 1].pres, state.cells[i - 1].temp) / state.cells[i].fluicol.MasEspFlu(state.cells[i].presaux, tmed);
        state.cells[i].bet = fabs(qc) / (fabs(qo) + fabs(qw) + fabs(qc));
    } else
        state.cells[i].bet = state.cells[i - 1].bet;
}

void applySteadyMassLeakSource(const SteadyStateState &state, int i, int mudaRGO, double tL, double tH, double bo, double ba, double rs, double tmed, double &trF, double &boI, double &baI, double &fwI) {
    ProFlu fluF;
    // define qual o fluido envolvido na fonte, se a pressÃƒÂ£o ambiente for menor do que a
    // pressao da tubulacao, fluido da tubulacao, senao, fluido definido como
    // fluido ambiente
    if (state.cells[i - 1].acsr.fontechk.presT > state.cells[i - 1].acsr.fontechk.pamb) {
        fluF = state.cells[i - 1].acsr.fontechk.fluidoP;
    } else {
        fluF = state.cells[i - 1].acsr.fontechk.fluidoPamb;
    }
    double api = state.cells[i - 1].flui.API;
    double rhololeo = 1000. * 141.5 / (131.5 + api);
    double rhololeoF = 1000. * 141.5 / (131.5 + fluF.API);
    double rholis = state.cells[i - 1].flui.MasEspLiq(state.cells[i - 1].pres, state.cells[i - 1].temp);
    double rholisF = fluF.MasEspLiq(state.cells[i - 1].pres, state.cells[i - 1].temp);
    double boinjl = fluF.BOFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
    double bainjl = fluF.BAFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
    double fwinjl = fluF.BSW * bainjl / (boinjl + bainjl * fluF.BSW - fluF.BSW * boinjl);
    // vazao de oleo standard antes da fonte:
    double qostd1 = 0.;
    if (state.cells[i - 1].flui.dStockTankVaporMassFraction < 1. - 1e-15)
        qostd1 = state.cells[i - 1].Mliqini * (1. - state.cells[i - 1].FW) * (1. - state.cells[i - 1].bet) / (bo * rholis);
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
    double hol = 1 - state.cells[i - 1].alf;
    double bet = state.cells[i - 1].bet;
    double bsw = state.cells[i - 1].FW;
    double rhog = state.cells[i - 1].flui.MasEspGas(state.cells[i - 1].pres, state.cells[i - 1].temp);
    double rhogST = state.cells[i - 1].flui.Deng * 1.225;
    // calculo do volume de leve na celula i, seguindo o equacionamento mostrado em relatorio,
    // nÃƒÂ£o ÃƒÂ© relevante para o permanente, mas deve ser calculado,
    // pois ÃƒÂ© utilizado como entrada no transiente
    state.cells[i - 1].VolLeveST = (((1 - hol) * rhog / rhogST) + hol * (1 - bet) * (1. - bsw) * rs / bo);
    if (state.cells[i - 1].VolLeveST < 1e-15)
        state.cells[i - 1].VolLeveST = 0.;

    if (state.input.trackRGO == -1) {
        // esta chave nÃƒÂ£o ÃƒÂ© utilizada, ÃƒÂ© mantida aqui como reserva, atualmente este
        // calculo nÃƒÂ£o esta funcionando, lembrando que arq.trackRGO assume apenas 2 valores, 0 ou 1
        if (hol > (*state.globals).localtiny && bet < (1. - (*state.globals).localtiny * 1e-6) && bsw < (1. - (*state.globals).localtiny) && mudaRGO == 1)
            state.cells[i].flui.RGO = state.cells[i - 1].VolLeveST * bo / (hol * (1 - bet) * (1 - bsw));
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
        double qw1;
        if ((1. - state.cells[i - 1].flui.BSW) > 0)
            qw1 = qostd1 * state.cells[i - 1].flui.BSW / (1. - state.cells[i - 1].flui.BSW);
        else
            qw1 = state.cells[i - 1].Mliqini * state.cells[i - 1].FW * (1. - state.cells[i - 1].bet) / (1000. * state.cells[i - 1].flui.Denag);
        double qw2;
        if ((1. - fluF.BSW) > 0)
            qw2 = qostd2 * fluF.BSW / (1. - fluF.BSW);
        else
            qw2 = state.cells[i - 1].fontemassLR * fwinjl / (1000. * fluF.Denag);

        if (fabs(qw1 + qw2 + qostd1 + qostd2) > 1e-15)
            state.cells[i].flui.BSW = (qw1 + qw2) / (qw1 + qw2 + qostd1 + qostd2);
        else
            state.cells[i].flui.BSW = state.cells[i - 1].flui.BSW;

        if (fabs(qostd1 + qostd2) > 1e-15 && state.input.flashCompleto == 0) {
            // reavaliacao dos pares necessarios para o modelo de viscosidade de oleo morto
            // ASTM. Isto ÃƒÂ© feito tendo sempre os mesmos valores de temperatura maxima e
            // e minima para a construcao dos pares
            state.cells[i].flui.TempL = tL;
            state.cells[i].flui.TempH = tH;
            state.cells[i].flui.LVisL = (qostd1 * state.cells[i - 1].flui.VisOM(tL) + qostd2 * fluF.VisOM(tL)) / (qostd1 + qostd2);
            state.cells[i].flui.LVisH = (qostd1 * state.cells[i - 1].flui.VisOM(tH) + qostd2 * fluF.VisOM(tH)) / (qostd1 + qostd2);
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
        double boN = state.cells[i].flui.BOFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        double baN = state.cells[i].flui.BAFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        double fwN = state.cells[i].flui.BSW * baN / (boN + baN * state.cells[i].flui.BSW - state.cells[i].flui.BSW * boN);
        boI = state.cells[i].flui.BOFunc(state.cells[i].presaux, tmed);
        baI = state.cells[i].flui.BAFunc(state.cells[i].presaux, tmed);
        fwI = state.cells[i].flui.BSW * baI / (boI + baI * state.cells[i].flui.BSW - state.cells[i].flui.BSW * boI);
        double qlpF = state.cells[i - 1].fontemassLR / fluF.MasEspLiq(state.cells[i - 1].pres, state.cells[i - 1].temp);
        double qlcF = state.cells[i - 1].fontemassCR / state.cells[i - 1].acsr.fontechk.fluidocol.MasEspFlu(state.cells[i - 1].pres, state.cells[i - 1].temp);
        double betN = (state.cells[i - 1].QL * state.cells[i - 1].bet + qlcF) / (state.cells[i - 1].QL + qlpF + qlcF);
        double qo = (state.cells[i - 1].QL + qlpF + qlcF) * (1 - fwN) * (1 - betN) * boI / bo;
        double qw;
        if (fwI < (1 - (*state.globals).localtiny))
            qw = fwI * qo / (1 - fwI);
        else
            qw = state.cells[i - 1].QL * (1 - state.cells[i - 1].bet) + qlpF;
        double qc = (state.cells[i - 1].QL * (state.cells[i - 1].bet) + qlcF) * state.cells[i - 1].fluicol.MasEspFlu(state.cells[i - 1].pres, state.cells[i - 1].temp) / state.cells[i].fluicol.MasEspFlu(state.cells[i].presaux, tmed);
        state.cells[i].bet = fabs(qc) / (fabs(qo) + fabs(qw) + fabs(qc));
    } else
        state.cells[i].bet = state.cells[i - 1].bet;
}

void applySteadyMassRadialPorous(const SteadyStateState &state, int i, int mudaRGO, double tL, double tH, double bo, double ba, double rs, double tmed, double &trF, double &boI, double &baI, double &fwI) {
    double api = state.cells[i - 1].flui.API;
    double rhololeo = 1000. * 141.5 / (131.5 + api);
    double rhololeoF = 1000. * 141.5 / (131.5 + state.cells[i - 1].acsr.radialPoro.flup.API);
    double rholis = state.cells[i - 1].flui.MasEspLiq(state.cells[i - 1].pres, state.cells[i - 1].temp);
    double rholisF;
    double boipr;
    double baipr;
    double fwipr;
    trF = 0.;
    if (fabs(state.cells[i - 1].fontemassLR) > 1e-15) {
        rholisF = state.cells[i - 1].acsr.radialPoro.flup.MasEspLiq(state.cells[i - 1].pres, state.cells[i - 1].temp);
        boipr = state.cells[i - 1].acsr.radialPoro.flup.BOFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        baipr = state.cells[i - 1].acsr.radialPoro.flup.BAFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        fwipr = state.cells[i - 1].acsr.radialPoro.BSW * baipr / (boipr + baipr * state.cells[i - 1].acsr.radialPoro.BSW - state.cells[i - 1].acsr.radialPoro.BSW * boipr);
    } else {
        rholisF = rholis;
        boipr = bo;
        baipr = ba;
        fwipr = state.cells[i - 1].FW;
    }
    // vazao de oleo standard antes da fonte:
    double qostd1 = 0.;
    if (state.cells[i - 1].flui.dStockTankVaporMassFraction < 1. - 1e-15)
        qostd1 = state.cells[i - 1].Mliqini * (1. - state.cells[i - 1].FW) * (1. - state.cells[i - 1].bet) / (bo * rholis);
    // vazao de oleo standard da fonte:
    double qostd2 = 0.;
    if (state.cells[i - 1].acsr.radialPoro.flup.dStockTankVaporMassFraction < 1. - 1e-15)
        qostd2 = state.cells[i - 1].fontemassLR * (1. - fwipr) / (boipr * rholisF);
    // vazao total de gas nas condicoes standard fruto da soma do gas transportado
    // da fronteira a esquerda
    // da celula e do gas associado da IPR
    double qgstd = qostd1 * state.cells[i - 1].flui.RGO + qostd2 * state.cells[i - 1].acsr.radialPoro.flup.RGO;
    double hol = 1 - state.cells[i - 1].alf;
    double bet = state.cells[i - 1].bet;
    double bsw = state.cells[i - 1].FW;
    double rhog = state.cells[i - 1].flui.MasEspGas(state.cells[i - 1].pres, state.cells[i - 1].temp);
    double rhogST = state.cells[i - 1].flui.Deng * 1.225;
    // calculo do volume de leve na celula i, seguindo o equacionamento mostrado em relatorio,
    // nÃƒÂ£o ÃƒÂ© relevante para o permanente, mas deve ser calculado,
    // pois ÃƒÂ© utilizado como entrada no transiente
    state.cells[i - 1].VolLeveST = (((1 - hol) * rhog / rhogST) + hol * (1 - bet) * (1. - bsw) * rs / bo);
    if (state.cells[i - 1].VolLeveST < 1e-15)
        state.cells[i - 1].VolLeveST = 0.;
    if (state.input.trackRGO == -1) {
        // esta chave nÃƒÂ£o ÃƒÂ© utilizada, ÃƒÂ© mantida aqui como reserva, atualmente este
        // calculo nÃƒÂ£o esta funcionando, lembrando que arq.trackRGO assume apenas 2 valores, 0 ou 1
        if (fabs(state.cells[i - 1].fontemassLR) > 1e-15) {
            if (hol > (*state.globals).localtiny && bet < (1. - (*state.globals).localtiny) && bsw < (1. - (*state.globals).localtiny) &&
                mudaRGO == 1 && fabs(state.cells[i - 1].fontemassLR) > 1e-15)
                state.cells[i].flui.RGO = state.cells[i - 1].VolLeveST * bo / (hol * (1 - bet) * (1 - bsw));
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
        double qw1;
        if ((1. - state.cells[i - 1].flui.BSW) > 0)
            qw1 = qostd1 * state.cells[i - 1].flui.BSW / (1. - state.cells[i - 1].flui.BSW);
        else
            qw1 = state.cells[i - 1].Mliqini * state.cells[i - 1].FW * (1. - state.cells[i - 1].bet) / (1000. * state.cells[i - 1].flui.Denag);
        double qw2;
        if ((1. - state.cells[i - 1].acsr.radialPoro.BSW) > 0)
            qw2 = qostd2 * state.cells[i - 1].acsr.radialPoro.BSW / (1. - state.cells[i - 1].acsr.radialPoro.BSW);
        else
            qw2 = state.cells[i - 1].fontemassLR * fwipr / (1000. * state.cells[i - 1].acsr.radialPoro.flup.Denag);

        if (fabs(qw1 + qw2 + qostd1 + qostd2) > 1e-15 &&
            fabs(state.cells[i - 1].fontemassLR) > 1e-15)
            state.cells[i].flui.BSW = ((qw1 + qostd1) * state.cells[i - 1].flui.BSW + (qw2 + qostd2) * state.cells[i - 1].acsr.radialPoro.BSW) / (qw1 + qw2 + qostd1 + qostd2);
        else
            state.cells[i].flui.BSW = state.cells[i - 1].flui.BSW;
        if ((qw1 + qw2) > 1e-15 && (state.cells[i - 1].fontemassLR) > 1e-15)
            state.cells[i].flui.Denag = (qw1 * state.cells[i - 1].flui.Denag +
                                    qw2 * state.cells[i - 1].acsr.radialPoro.flup.Denag) /
                                   (qw1 + qw2);
        else
            state.cells[i].flui.Denag = state.cells[i - 1].flui.Denag;
        if (fabs(qostd1 + qostd2) > 1e-15 && state.input.flashCompleto == 0 && (state.cells[i - 1].fontemassLR) > 1e-15) {
            // reavaliacao dos pares necessarios para o modelo de viscosidade de oleo morto
            // ASTM. Isto ÃƒÂ© feito tendo sempre os mesmos valores de temperatura maxima e
            // e minima para a construcao dos pares
            state.cells[i].flui.TempL = tL;
            state.cells[i].flui.TempH = tH;
            state.cells[i].flui.LVisL = (qostd1 * state.cells[i - 1].flui.VisOM(tL) + qostd2 * state.cells[i - 1].acsr.radialPoro.flup.VisOM(tL)) / (qostd1 + qostd2);
            state.cells[i].flui.LVisH = (qostd1 * state.cells[i - 1].flui.VisOM(tH) + qostd2 * state.cells[i - 1].acsr.radialPoro.flup.VisOM(tH)) / (qostd1 + qostd2);
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
        double boN = state.cells[i].flui.BOFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        double baN = state.cells[i].flui.BAFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        double fwN = state.cells[i].flui.BSW * baN / (boN + baN * state.cells[i].flui.BSW - state.cells[i].flui.BSW * boN);
        boI = state.cells[i].flui.BOFunc(state.cells[i].presaux, tmed);
        baI = state.cells[i].flui.BAFunc(state.cells[i].presaux, tmed);
        fwI = state.cells[i].flui.BSW * baI / (boI + baI * state.cells[i].flui.BSW - state.cells[i].flui.BSW * boI);
        double qlpF = state.cells[i - 1].fontemassLR / state.cells[i - 1].acsr.radialPoro.flup.MasEspLiq(state.cells[i - 1].pres, state.cells[i - 1].temp);
        double qlcF = 0.;
        double betN = (state.cells[i - 1].QL * state.cells[i - 1].bet + qlcF) / (state.cells[i - 1].QL + qlpF + qlcF);
        double qo = (state.cells[i - 1].QL + qlpF + qlcF) * (1 - fwN) * (1 - betN) * boI / bo;
        double qw;
        if (fwI < (1 - (*state.globals).localtiny))
            qw = fwI * qo / (1 - fwI);
        else
            qw = state.cells[i - 1].QL * (1 - state.cells[i - 1].bet) + qlpF;
        double qc = (state.cells[i - 1].QL * (state.cells[i - 1].bet) + qlcF) * state.cells[i - 1].fluicol.MasEspFlu(state.cells[i - 1].pres, state.cells[i - 1].temp) / state.cells[i].fluicol.MasEspFlu(state.cells[i].presaux, tmed);
        state.cells[i].bet = fabs(qc) / (fabs(qo) + fabs(qw) + fabs(qc));
    } else
        state.cells[i].bet = state.cells[i - 1].bet;
}

void applySteadyMassPorous2D(const SteadyStateState &state, int i, int mudaRGO, double tL, double tH, double bo, double ba, double rs, double tmed, double &trF, double &boI, double &baI, double &fwI) {
    double api = state.cells[i - 1].flui.API;
    double rhololeo = 1000. * 141.5 / (131.5 + api);
    double rhololeoF = 1000. * 141.5 / (131.5 + state.cells[i - 1].acsr.poroso2D.dados.flup.API);
    double rholis = state.cells[i - 1].flui.MasEspLiq(state.cells[i - 1].pres, state.cells[i - 1].temp);
    double rholisF;
    double boipr;
    double baipr;
    double fwipr;
    trF = 0.;
    if (fabs(state.cells[i - 1].fontemassLR) > 1e-15) {
        rholisF = state.cells[i - 1].acsr.poroso2D.dados.flup.MasEspLiq(state.cells[i - 1].pres, state.cells[i - 1].temp);
        boipr = state.cells[i - 1].acsr.poroso2D.dados.flup.BOFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        baipr = state.cells[i - 1].acsr.poroso2D.dados.flup.BAFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        fwipr = state.cells[i - 1].acsr.poroso2D.dados.transfer.BSW * baipr / (boipr + baipr * state.cells[i - 1].acsr.poroso2D.dados.transfer.BSW - state.cells[i - 1].acsr.poroso2D.dados.transfer.BSW * boipr);
    } else {
        rholisF = rholis;
        boipr = bo;
        baipr = ba;
        fwipr = state.cells[i - 1].FW;
    }
    // vazao de oleo standard antes da fonte:
    double qostd1 = 0.;
    if (state.cells[i - 1].flui.dStockTankVaporMassFraction < 1. - 1e-15)
        qostd1 = state.cells[i - 1].Mliqini * (1. - state.cells[i - 1].FW) * (1. - state.cells[i - 1].bet) / (bo * rholis);
    // vazao de oleo standard da fonte:
    double qostd2 = 0.;
    if (state.cells[i - 1].acsr.poroso2D.dados.flup.dStockTankVaporMassFraction < 1. - 1e-15)
        qostd2 = state.cells[i - 1].fontemassLR * (1. - fwipr) / (boipr * rholisF);
    // vazao total de gas nas condicoes standard fruto da soma do gas transportado
    // da fronteira a esquerda
    // da celula e do gas associado da IPR
    double qgstd = qostd1 * state.cells[i - 1].flui.RGO + qostd2 * state.cells[i - 1].acsr.poroso2D.dados.flup.RGO;
    double hol = 1 - state.cells[i - 1].alf;
    double bet = state.cells[i - 1].bet;
    double bsw = state.cells[i - 1].FW;
    double rhog = state.cells[i - 1].flui.MasEspGas(state.cells[i - 1].pres, state.cells[i - 1].temp);
    double rhogST = state.cells[i - 1].flui.Deng * 1.225;
    // calculo do volume de leve na celula i, seguindo o equacionamento mostrado em relatorio,
    // nÃƒÂ£o ÃƒÂ© relevante para o permanente, mas deve ser calculado,
    // pois ÃƒÂ© utilizado como entrada no transiente
    state.cells[i - 1].VolLeveST = (((1 - hol) * rhog / rhogST) + hol * (1 - bet) * (1. - bsw) * rs / bo);
    if (state.cells[i - 1].VolLeveST < 1e-15)
        state.cells[i - 1].VolLeveST = 0.;
    if (state.input.trackRGO == -1) {
        // esta chave nÃƒÂ£o ÃƒÂ© utilizada, ÃƒÂ© mantida aqui como reserva, atualmente este
        // calculo nÃƒÂ£o esta funcionando, lembrando que arq.trackRGO assume apenas 2 valores, 0 ou 1
        if (fabs(state.cells[i - 1].fontemassLR) > 1e-15) {
            if (hol > (*state.globals).localtiny && bet < (1. - (*state.globals).localtiny) && bsw < (1. - (*state.globals).localtiny) &&
                mudaRGO == 1 && fabs(state.cells[i - 1].fontemassLR) > 1e-15)
                state.cells[i].flui.RGO = state.cells[i - 1].VolLeveST * bo / (hol * (1 - bet) * (1 - bsw));
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
        double qw1;
        if ((1. - state.cells[i - 1].flui.BSW) > 0)
            qw1 = qostd1 * state.cells[i - 1].flui.BSW / (1. - state.cells[i - 1].flui.BSW);
        else
            qw1 = state.cells[i - 1].Mliqini * state.cells[i - 1].FW * (1. - state.cells[i - 1].bet) / (1000. * state.cells[i - 1].flui.Denag);
        double qw2;
        if ((1. - state.cells[i - 1].acsr.poroso2D.dados.transfer.BSW) > 0)
            qw2 = qostd2 * state.cells[i - 1].acsr.poroso2D.dados.transfer.BSW / (1. - state.cells[i - 1].acsr.poroso2D.dados.transfer.BSW);
        else
            qw2 = state.cells[i - 1].fontemassLR * fwipr / (1000. * state.cells[i - 1].acsr.poroso2D.dados.flup.Denag);

        if (fabs(qw1 + qw2 + qostd1 + qostd2) > 1e-15 &&
            (state.cells[i - 1].fontemassLR) > 1e-15)
            state.cells[i].flui.BSW = ((qw1 + qostd1) * state.cells[i - 1].flui.BSW + (qw2 + qostd2) * state.cells[i - 1].acsr.poroso2D.dados.transfer.BSW) / (qw1 + qw2 + qostd1 + qostd2);
        else
            state.cells[i].flui.BSW = state.cells[i - 1].flui.BSW;
        if (fabs(qw1 + qw2) > 1e-15 &&
            (state.cells[i - 1].fontemassLR) > 1e-15)
            state.cells[i].flui.Denag = (qw1 * state.cells[i - 1].flui.Denag +
                                    qw2 * state.cells[i - 1].acsr.poroso2D.dados.flup.Denag) /
                                   (qw1 + qw2);
        else
            state.cells[i].flui.Denag = state.cells[i - 1].flui.Denag;
        if (fabs(qostd1 + qostd2) > 1e-15 && state.input.flashCompleto == 0 && (state.cells[i - 1].fontemassLR) > 1e-15) {
            // reavaliacao dos pares necessarios para o modelo de viscosidade de oleo morto
            // ASTM. Isto ÃƒÂ© feito tendo sempre os mesmos valores de temperatura maxima e
            // e minima para a construcao dos pares
            state.cells[i].flui.TempL = tL;
            state.cells[i].flui.TempH = tH;
            state.cells[i].flui.LVisL = (qostd1 * state.cells[i - 1].flui.VisOM(tL) + qostd2 * state.cells[i - 1].acsr.poroso2D.dados.flup.VisOM(tL)) / (qostd1 + qostd2);
            state.cells[i].flui.LVisH = (qostd1 * state.cells[i - 1].flui.VisOM(tH) + qostd2 * state.cells[i - 1].acsr.poroso2D.dados.flup.VisOM(tH)) / (qostd1 + qostd2);
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
        double boN = state.cells[i].flui.BOFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        double baN = state.cells[i].flui.BAFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        double fwN = state.cells[i].flui.BSW * baN / (boN + baN * state.cells[i].flui.BSW - state.cells[i].flui.BSW * boN);
        boI = state.cells[i].flui.BOFunc(state.cells[i].presaux, tmed);
        baI = state.cells[i].flui.BAFunc(state.cells[i].presaux, tmed);
        fwI = state.cells[i].flui.BSW * baI / (boI + baI * state.cells[i].flui.BSW - state.cells[i].flui.BSW * boI);
        double qlpF = state.cells[i - 1].fontemassLR / state.cells[i - 1].acsr.poroso2D.dados.flup.MasEspLiq(state.cells[i - 1].pres, state.cells[i - 1].temp);
        double qlcF = 0.;
        double betN = (state.cells[i - 1].QL * state.cells[i - 1].bet + qlcF) / (state.cells[i - 1].QL + qlpF + qlcF);
        double qo = (state.cells[i - 1].QL + qlpF + qlcF) * (1 - fwN) * (1 - betN) * boI / bo;
        double qw;
        if (fwI < (1 - (*state.globals).localtiny))
            qw = fwI * qo / (1 - fwI);
        else
            qw = state.cells[i - 1].QL * (1 - state.cells[i - 1].bet) + qlpF;
        double qc = (state.cells[i - 1].QL * (state.cells[i - 1].bet) + qlcF) * state.cells[i - 1].fluicol.MasEspFlu(state.cells[i - 1].pres, state.cells[i - 1].temp) / state.cells[i].fluicol.MasEspFlu(state.cells[i].presaux, tmed);
        state.cells[i].bet = fabs(qc) / (fabs(qo) + fabs(qw) + fabs(qc));
    } else
        state.cells[i].bet = state.cells[i - 1].bet;
}

void finalizeSteadyMassWithLiquid(const SteadyStateState &state, int i, double fwI, double tmed, double &bo, double &rs, double &pmed,
                                            double &rhog, double &rhol, double &qo, double &masoleo) {
    if (state.input.tipoFluido == 0) {
        rs = state.cells[i].flui.RS(pmed, tmed);
        bo = state.cells[i].flui.BOFunc(pmed, tmed, rs);
        rs = rs * 6.29 / 35.31467;
        double rhogstd = state.cells[i].flui.Deng * 1.225;
        double rhololeostd = 1000. * 141.5 / (131.5 + state.cells[i].flui.API);
        double rhoa = state.cells[i].flui.Denag * 1000.;

        double denom = state.cells[i].flui.RGO * rhogstd + rhololeostd;
        if (fabs(1 - state.cells[i].flui.BSW) > (*state.globals).localtiny) {
            denom += rhoa * state.cells[i].flui.BSW / (1 - state.cells[i].flui.BSW);
            if (fabs(1 - state.cells[i].bet) > (*state.globals).localtiny * 1e-6 && fabs(state.cells[i].bet) > (*state.globals).localtiny * 1e-10) {
                denom += (state.cells[i].fluicol.MasEspFlu(pmed, tmed) * state.cells[i].bet / (1 - state.cells[i].bet)) * bo * (1. + fwI / (1 - fwI));
            }
        }
        if (fabs(1 - state.cells[i].bet) > (*state.globals).localtiny * 1e-6 && fabs(1 - state.cells[i].flui.BSW) > (*state.globals).localtiny)
            qo = state.cells[i].MC / denom;

        if (fabs(1 - state.cells[i].bet) > (*state.globals).localtiny * 1e-6 && fabs(1 - state.cells[i].flui.BSW) > (*state.globals).localtiny) {
            if ((*state.globals).tipoFluidoRedeGlob == 0)
                state.cells[i].Mliqini = state.cells[i].MC - (state.cells[i].flui.RGO - rs) * rhogstd * qo;
            else {
                masoleo = qo * bo * state.cells[i].flui.MasEspoleo(pmed, tmed, rs) + qo * rhoa * state.cells[i].flui.BSW / (1 - state.cells[i].flui.BSW);
                state.cells[i].Mliqini = masoleo + state.cells[i].MComp;
            }
        } else
            state.cells[i].Mliqini = state.cells[i].MC - (state.cells[i - 1].MC - state.cells[i - 1].Mliqini) - state.cells[i - 1].fontemassGR;
    } else {
        double tit = state.cells[i].flui.FracMass(pmed, tmed);
        state.cells[i].Mliqini = state.cells[i].MC * (1. - tit);
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

    double trF = 0.;
    double sinalQ = 1.;
    if (i > 1) {
        if (fabs(state.cells[i - 1].QL + state.cells[i - 1].QG) > 1e-15)
            sinalQ = (state.cells[i - 1].QL + state.cells[i - 1].QG) / fabs(state.cells[i - 1].QL + state.cells[i - 1].QG);
    } else {
        if (fabs(state.cells[i - 1].fontemassLR + state.cells[i - 1].fontemassCR + state.cells[i - 1].fontemassGR) > 1e-15)
            sinalQ = (state.cells[i - 1].fontemassLR + state.cells[i - 1].fontemassCR + state.cells[i - 1].fontemassGR) /
                     fabs(state.cells[i - 1].fontemassLR + state.cells[i - 1].fontemassCR + state.cells[i - 1].fontemassGR);
    }
    double boI;
    double baI;
    double fwI = 0;
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
        applySteadyMassWithoutSource(state, i, mudaRGO, bo, rs, tmed, boI, baI, fwI);
    }
    // caso em que se tem uma fonte de gas na celula i-1, o que mudara a RGO e a densidade de gas em i
    else if (state.cells[i - 1].acsr.tipo == 1 && state.cells[i - 1].acsr.injg.seco == 1) {
        applySteadyMassDryGasInjection(state, i, mudaRGO, tL, tH, bo, ba, rs, tmed, trF, boI, baI, fwI);
    } else if (state.cells[i - 1].acsr.tipo == 1 && state.cells[i - 1].acsr.injg.seco == 0) {
        applySteadyMassWetGasInjection(state, i, mudaRGO, tL, tH, bo, ba, rs, tmed, trF, boI, baI, fwI);
    }
    // caso de fonte de liquido na celula i-1:
    else if (state.cells[i - 1].acsr.tipo == 2) {
        applySteadyMassLiquidInjection(state, i, mudaRGO, tL, tH, bo, ba, rs, tmed, trF, boI, baI, fwI);
    }
    // caso de IPR na celula i-1:
    else if (state.cells[i - 1].acsr.tipo == 3) {
        applySteadyMassInflowPerformance(state, i, mudaRGO, tL, tH, bo, ba, rs, tmed, trF, boI, baI, fwI);
    }
    // caso de fonte de massa na celula i-1:
    else if (state.cells[i - 1].acsr.tipo == 10) {
        applySteadyMassMultipleSource(state, i, mudaRGO, tL, tH, bo, ba, rs, tmed, trF, boI, baI, fwI);
    }
    // caso especial, fonte de vazamento:
    else if (state.cells[i - 1].acsr.tipo == 9 && state.cells[i - 1].acsr.fontechk.abertura > 1e-6) {
        applySteadyMassLeakSource(state, i, mudaRGO, tL, tH, bo, ba, rs, tmed, trF, boI, baI, fwI);
    } else if (state.cells[i - 1].acsr.tipo == 15) {
        applySteadyMassRadialPorous(state, i, mudaRGO, tL, tH, bo, ba, rs, tmed, trF, boI, baI, fwI);
    } else if (state.cells[i - 1].acsr.tipo == 16) {
        applySteadyMassPorous2D(state, i, mudaRGO, tL, tH, bo, ba, rs, tmed, trF, boI, baI, fwI);
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
    double qo = 0.;
    double masoleo;
    if (state.cells[i].flui.RGO < 1e7) { // Caso exista liquido
        finalizeSteadyMassWithLiquid(state, i, fwI, tmed, bo, rs, pmed, rhog, rhol, qo, masoleo);
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
        state.cells[i].fluicol.TR = (state.cells[i - 1].fluicol.TR * state.cells[i - 1].MComp + state.cells[i - 1].fontemassCR * trF) / state.cells[i].MComp;
        state.cells[i].fluicol.TR = state.cells[i].fluicol.TR + (0.5 * state.cells[i - 1].duto.area * state.cells[i - 1].dx * hol0 +
                                                       0.5 * state.cells[i].duto.area * state.cells[i].dx * hol1) /
                                                          state.cells[i].QL;
    } else
        state.cells[i].fluicol.TR = state.cells[i - 1].fluicol.TR;
}

double areaChangePressureDrop(const SteadyStateState &state, int i, double rhomix, double rey, double jmix) {
    double dpArea = 0.;
    if ((state.cells[i].duto.area != state.cells[i].dutoR.area && rey > 2400) && (state.cells[i].mudaArea == 1 && (fabs(jmix) >= 0.1))) {
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
    double dia = state.cells[state.lastCell].duto.a;
    double area = 0.25 * M_PI * dia * dia;
    double si = state.cells[state.lastCell].duto.peri;
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

    double re1;
    if (state.cells[state.lastCell].duto.revest == 0)
        re1 = state.cells[state.lastCell].Rey(state.cells[state.lastCell].duto.a, j, rhomix, viscmix);
    else {
        double dhid = 4 * area / si;
        re1 = state.cells[state.lastCell].Rey(dhid, j, rhomix, viscmix);
    }
    double f1 = state.cells[state.lastCell].fric(re1, state.cells[state.lastCell].duto.rug / dia);
    double gradfric = state.cells[state.lastCell].dPdLFric * (0.5 * f1 * rhomix * (fabs(j) * j) * si * dx / area);
    double gradhidro = state.cells[state.lastCell].dPdLHidro * (9.82 * sin(state.cells[state.lastCell].duto.teta) * rhomix * dx);
    return -(1. * gradfric + gradhidro) / 98066.5;
}

void advanceUpstreamSteadyPressure(const SteadyStateState &state, int i, int RK) {

    if (state.input.tipoModeloDrift == 1) {
        double dx = 0.5 * state.cells[i].dxL;
        double dia = state.cells[i].dutoL.a;
        double area = 0.25 * M_PI * dia * dia;
        double si = state.cells[i].dutoL.peri;
        if (RK == 0) {
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

            double re1;
            if (state.cells[i].dutoL.revest == 0)
                re1 = state.cells[i - 1].Rey(state.cells[i].dutoL.a, j, rhomix, viscmix);
            else {
                double dhid = 4 * area / si;
                re1 = state.cells[i - 1].Rey(dhid, j, rhomix, viscmix);
            }
            double f1 = state.cells[i - 1].fric(re1, state.cells[i].dutoL.rug / dia);
            if (i > 1 && state.cells[i - 1].fluicol.tipoF == 2) {
                f1 *= (1 - state.cells[i - 1].dR);
            }
            double gradfric = state.cells[i - 1].dPdLFric * (0.5 * f1 * rhomix * (fabs(j) * j) * si * dx / area);
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

            double re1;
            if (state.cells[i].dutoL.revest == 0)
                re1 = state.cells[i - 1].Rey(state.cells[i].dutoL.a, j, rhomix, viscmix);
            else {
                double dhid = 4 * area / si;
                re1 = state.cells[i - 1].Rey(dhid, j, rhomix, viscmix);
            }
            double f1 = state.cells[i - 1].fric(re1, state.cells[i].dutoL.rug / dia);
            double gradfric = state.cells[i - 1].dPdLFric * (0.5 * f1 * rhomix * (fabs(j) * j) * si * dx / area);
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
        double dia = state.cells[i].dutoL.a;
        double area = 0.25 * M_PI * dia * dia;
        double si = state.cells[i].dutoL.peri;
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

void advanceDownstreamSteadyPressure(const SteadyStateState &state, int i, int RK) {

    double dx;
    double dia;
    double area;
    double si;
    double alfmed;
    double betmed;
    double rhog;
    double rhol;
    double ugsmed;
    double ulsmed;
    double j;

    double rhomix;
    double viscmix;
    double re1;
    double f1;
    double gradfric;
    double gradhidro;

    double tmed;

    double dpArea = 0.;

    dx = 0.5 * state.cells[i].dx;
    dia = state.cells[i].duto.a;
    area = 0.25 * M_PI * dia * dia;
    si = state.cells[i].duto.peri;

    if (state.input.tipoModeloDrift == 1) {
        if (RK == 0) {
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
                re1 = state.cells[i].Rey(state.cells[i].duto.a, j, rhomix, viscmix);
            else {
                double dhid = 4 * area / si;
                re1 = state.cells[i].Rey(dhid, j, rhomix, viscmix);
            }
            f1 = state.cells[i].fric(re1, state.cells[i].duto.rug / dia);
            if (state.cells[i].fluicol.tipoF == 2) {
                double ulmed = 0.;
                double reL;
                if (fabs(state.cells[i].MComp) > 1e-15) {
                    if (alfmed < 1 - 1e-15)
                        ulmed = ulsmed / (1 - alfmed);
                    else
                        ulmed = 0.;
                    if (state.cells[i].duto.revest == 0)
                        reL = state.cells[i].Rey(state.cells[i].duto.a, ulmed, rhol, visl);
                    else {
                        double dhid = 4 * area / si;
                        reL = state.cells[i].Rey(dhid, ulmed, rhol, visl);
                    }
                    state.cells[i].dR = state.cells[i].fluicol.calcDR(reL);
                } else
                    state.cells[i].dR = 0.;
                f1 *= (1 - state.cells[i].dR);
            }
            gradfric = state.cells[i].dPdLFric * (0.5 * f1 * rhomix * (fabs(j) * j) * si * dx / area);
            gradhidro = state.cells[i].dPdLHidro * (9.82 * sin(state.cells[i].duto.teta) * rhomix * dx);
            if (state.cells[i].mudaArea == 1)
                dpArea = areaChangePressureDrop(state, i - 1, rhomix, re1, fabs(j));

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
                re1 = state.cells[i].Rey(state.cells[i].duto.a, j, rhomix, viscmix);
            else {
                double dhid = 4 * area / si;
                re1 = state.cells[i].Rey(dhid, j, rhomix, viscmix);
            }
            f1 = state.cells[i].fric(re1, state.cells[i].duto.rug / dia);
            gradfric = state.cells[i].dPdLFric * (0.5 * f1 * rhomix * (fabs(j) * j) * si * dx / area);
            gradhidro = state.cells[i].dPdLHidro * (9.82 * sin(state.cells[i].duto.teta) * rhomix * dx);
            if (state.cells[i].mudaArea == 1)
                dpArea = areaChangePressureDrop(state, i - 1, rhomix, re1, fabs(j));

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

    double fwd;
    double fwe;

    double razdx = state.cells[i].dxR / (state.cells[i].dx + state.cells[i].dxR);
    double razdxL = state.cells[i].dx / (state.cells[i].dx + state.cells[i].dxL);
    double tmed;
    tmed = razdx * state.cells[i + 1].temp + (1 - razdx) * state.cells[i].temp;
    double tmed0;
    if (i > 1) {
        tmed0 = razdxL * state.cells[i].temp + (1 - razdxL) * state.cells[i - 1].temp;
    } else
        tmed0 = state.cells[i].temp;

    double bo1 = state.cells[i].flui.BOFunc(state.cells[i].pres, state.cells[i].temp);
    double ba1 = state.cells[i].flui.BAFunc(state.cells[i].pres, state.cells[i].temp);
    fwd = state.cells[i].flui.BSW * ba1 / (bo1 + ba1 * state.cells[i].flui.BSW - state.cells[i].flui.BSW * bo1);
    double rsR = state.cells[i].flui.RS(state.cells[i + 1].presaux, tmed);
    double boR = state.cells[i].flui.BOFunc(state.cells[i + 1].presaux, tmed);
    double bo0;
    double ba0;
    double rsL;
    double boL;
    double baL;
    double dengD = state.cells[i].flui.Deng;
    double dengE;
    double rD = state.cells[i].flui.rDgD;
    double rE;
    if (i > 0) {
        bo0 = state.cells[i - 1].flui.BOFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        ba0 = state.cells[i - 1].flui.BAFunc(state.cells[i - 1].pres, state.cells[i - 1].temp);
        fwe = state.cells[i - 1].flui.BSW * ba0 / (bo0 + ba0 * state.cells[i - 1].flui.BSW - state.cells[i - 1].flui.BSW * bo0);
        rsL = state.cells[i - 1].flui.RS(state.cells[i].presaux, tmed0);
        boL = state.cells[i - 1].flui.BOFunc(state.cells[i].presaux, tmed0);
        baL = state.cells[i - 1].flui.BAFunc(state.cells[i].presaux, tmed0);
        dengE = state.cells[i - 1].flui.Deng;
        rE = state.cells[i - 1].flui.rDgD;
    } else {
        bo0 = state.cells[i].flui.BOFunc(state.cells[i].pres, state.cells[i].temp);
        ba0 = state.cells[i].flui.BAFunc(state.cells[i].pres, state.cells[i].temp);

        fwe = state.cells[i].flui.BSW * ba0 / (bo0 + ba0 * state.cells[i].flui.BSW - state.cells[i].flui.BSW * bo0);
        rsL = state.cells[i].flui.RS(state.cells[i].presaux, tmed0);
        boL = state.cells[i].flui.BOFunc(state.cells[i].presaux, tmed0);
        baL = state.cells[i].flui.BAFunc(state.cells[i].presaux, tmed0);
        dengE = state.cells[i].flui.Deng;
        rE = state.cells[i].flui.rDgD;
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
        state.cells[i].transmassR = (-(state.cells[i + 1].QL * (1. - betI) * rD * dengD * 1.225 * (1. - fwd) * rsR * (6.29 / 35.31467) / boR) + (state.cells[i].QL * (1. - betL) * rE * dengE * 1.225 * (1. - fwe) * rsL * (6.29 / 35.31467) / boL));
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

}  // namespace sisprod::steady
