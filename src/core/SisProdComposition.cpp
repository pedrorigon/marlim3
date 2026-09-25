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

void transportBlackOilProperties(const CompositionState &state, ProFlu fluiRev) {
    Vcr<double> rgo(state.lastCell);
    Vcr<double> dg(state.lastCell);
    Vcr<double> yco2(state.lastCell);
    Vcr<double> API(state.lastCell);
    Vcr<double> BSW(state.lastCell);
    Vcr<double> denag(state.lastCell);
    Vcr<double> VISCL(state.lastCell);
    Vcr<double> VISCH(state.lastCell);
    double dt = state.cells[1].dt;
    double tL;
    if (state.input.flashCompleto == 1 && (state.input.tabent.tmin - 0) > (*state.globals).localtiny)
        tL = state.input.tabent.tmin + 0.1;
    else
        tL = 0.;
    double tH;
    if (state.input.flashCompleto == 1 && (70 - state.input.tabent.tmax) > (*state.globals).localtiny)
        tH = state.input.tabent.tmax - 0.1;
    else
        tH = 70.;
    int imin = 1;
    if ((*state.globals).chaverede != 0 && (state.cells[0].acsr.tipo == 10 || state.input.ConContEntrada > 0))
        imin = 0;
    if (state.cells[0].acsr.tipo == 15 && (state.cells[0].fontemassCR + state.cells[0].fontemassGR + state.cells[0].fontemassLR) > 0.)
        state.cells[0].flui.BSW = state.cells[0].acsr.radialPoro.BSW;
    else if (state.cells[0].acsr.tipo == 16 && (state.cells[0].fontemassCR + state.cells[0].fontemassGR + state.cells[0].fontemassLR) > 0.)
        state.cells[0].flui.BSW = state.cells[0].acsr.poroso2D.dados.transfer.BSW;
#pragma omp parallel for num_threads((*state.globals).ntrd)
    for (int i = imin; i < state.lastCell; i++) {
        double MultOe;
        double MultOd;
        double at = state.cells[i].duto.area;
        double dx = state.cells[i].dx;
        double ti0;
        if (state.cells[i].VTemper < 0.)
            ti0 = state.cells[i].temp;
        else {
            if (i > 0)
                ti0 = state.cells[i - 1].temp;
            else if (state.input.ConContEntrada == 1)
                ti0 = state.inletTemperature;
            else
                ti0 = state.cells[i].temp;
        }
        double ti1 = state.cells[i].temp;
        if (state.cells[i + 1].VTemper < 0.)
            ti1 = state.cells[i + 1].temp;

        double rgo0;
        double betI0;
        double bo0;
        double ba0;
        double bsw0;
        double rs0;
        double dg0O;
        double yco20O;
        double API0;
        double BSW0;
        double denag0;
        double viscL0;
        double viscH0;
        double razdgd0;
        double razdgl0;
        if (i > 0 || state.input.ConContEntrada == 0) {
            if (i > 0 && state.cells[i].QG >= 0.)
                betI0 = state.cells[i - 1].betPigD;
            else
                betI0 = state.cells[i].betPigE;
        } else {
            if (state.cells[i].QG >= 0.)
                betI0 = state.inletCompletionFraction;
            else
                betI0 = state.cells[i].betPigE;
        }
        if ((i > 0 || state.input.ConContEntrada == 1) && state.cells[i].QL >= 0.) {
            double pe;
            double te;
            if (i > 0) {
                pe = state.cells[i - 1].pres;
                te = state.cells[i - 1].temp;
            } else {
                pe = state.inletPressure;
                te = state.inletTemperature;
            }
            rgo0 = (*state.cells[i].fluiL).RGO;
            if (state.input.ConContEntrada == 0)
                betI0 = state.cells[i - 1].betPigD; // testeBeta
            else
                betI0 = state.inletCompletionFraction; // testeBeta
            rs0 = (*state.cells[i].fluiL).RS(pe, te);
            bo0 = (*state.cells[i].fluiL).BOFunc(pe, te, rs0);
            ba0 = (*state.cells[i].fluiL).BAFunc(pe, te);
            bsw0 = (*state.cells[i].fluiL).BSW * ba0 / (bo0 + ba0 * (*state.cells[i].fluiL).BSW - (*state.cells[i].fluiL).BSW * bo0);
            rs0 = rs0 * 6.29 / 35.31467;
            dg0O = (*state.cells[i].fluiL).Deng;
            razdgd0 = 1 / (*state.cells[i].fluiL).rDgD;
            razdgl0 = 1 / (*state.cells[i].fluiL).rDgL;
            yco20O = (*state.cells[i].fluiL).yco2;
            API0 = (*state.cells[i].fluiL).API;
            BSW0 = (*state.cells[i].fluiL).BSW;
            denag0 = (*state.cells[i].fluiL).Denag;
            viscL0 = 0 * 30 + 1 * (*state.cells[i].fluiL).VisOM(tL);
            viscH0 = 0 * 20 + 1 * (*state.cells[i].fluiL).VisOM(tH);
        } else {
            betI0 = state.cells[i].betPigE; // testebeta
            rgo0 = state.cells[i].flui.RGO;
            rs0 = state.cells[i].flui.RS(state.cells[i].pres, state.cells[i].temp);
            bo0 = state.cells[i].flui.BOFunc(state.cells[i].pres, state.cells[i].temp, rs0);
            ba0 = state.cells[i].flui.BAFunc(state.cells[i].pres, state.cells[i].temp);
            bsw0 = state.cells[i].flui.BSW * ba0 / (bo0 + ba0 * state.cells[i].flui.BSW - state.cells[i].flui.BSW * bo0);
            rs0 = rs0 * 6.29 / 35.31467;
            dg0O = state.cells[i].flui.Deng;
            razdgd0 = 1 / state.cells[i].flui.rDgD;
            razdgl0 = 1 / state.cells[i].flui.rDgL;
            yco20O = state.cells[i].flui.yco2;
            API0 = state.cells[i].flui.API;
            BSW0 = state.cells[i].flui.BSW;
            denag0 = state.cells[i].flui.Denag;
            viscL0 = 0 * 30 + 1 * state.cells[i].flui.VisOM(tL);
            viscH0 = 0 * 20 + 1 * state.cells[i].flui.VisOM(tH);
        }
        if (bo0 < 1e-15)
            bo0 = 1e-15;

        double rgo1 = state.cells[i].flui.RGO;
        double betI1 = state.cells[i].betPigD;
        double rs1 = state.cells[i].flui.RS(state.cells[i].pres, state.cells[i].temp);
        double bo1 = state.cells[i].flui.BOFunc(state.cells[i].pres, state.cells[i].temp, rs1);
        double ba1 = state.cells[i].flui.BAFunc(state.cells[i].pres, state.cells[i].temp);
        double bsw1 = state.cells[i].flui.BSW * ba1 / (bo1 + ba1 * state.cells[i].flui.BSW - state.cells[i].flui.BSW * bo1);
        rs1 = rs1 * 6.29 / 35.31467;
        double dg1O = state.cells[i].flui.Deng;
        double yco21O = state.cells[i].flui.yco2;
        double API1 = state.cells[i].flui.API;
        double BSW1 = state.cells[i].flui.BSW;
        double denag1 = state.cells[i].flui.Denag;
        double viscL1 = 0 * 30 + 1 * state.cells[i].flui.VisOM(tL);
        double viscH1 = 0 * 20 + 1 * state.cells[i].flui.VisOM(tH);
        double razdgd1 = 1 / state.cells[i].flui.rDgD;
        double razdgl1 = 1 / state.cells[i].flui.rDgL;

        // betI1 = celula[i + 1].betPigE;    //duvidabeta
        if (state.cells[i + 1].QL < 0.) {
            betI1 = state.cells[i + 1].betPigE; // testeBeta
            rgo1 = state.cells[i + 1].flui.RGO;
            rs1 = state.cells[i + 1].flui.RS(state.cells[i + 1].pres, state.cells[i + 1].temp);
            bo1 = state.cells[i + 1].flui.BOFunc(state.cells[i + 1].pres, state.cells[i + 1].temp, rs1);
            ba1 = state.cells[i + 1].flui.BAFunc(state.cells[i + 1].pres, state.cells[i + 1].temp);
            bsw1 = state.cells[i + 1].flui.BSW * ba1 / (bo1 + ba1 * state.cells[i + 1].flui.BSW - state.cells[i + 1].flui.BSW * bo1);
            rs1 = rs1 * 6.29 / 35.31467;
            dg1O = state.cells[i + 1].flui.Deng;
            razdgd1 = 1 / state.cells[i + 1].flui.rDgD;
            razdgl1 = 1 / state.cells[i + 1].flui.rDgL;
            yco21O = state.cells[i + 1].flui.yco2;
            API1 = state.cells[i + 1].flui.API;
            BSW1 = state.cells[i + 1].flui.BSW;
            denag1 = state.cells[i + 1].flui.Denag;
            viscL1 = 0 * 30 + 1 * state.cells[i + 1].flui.VisOM(tL);
            viscH1 = 0 * 20 + 1 * state.cells[i + 1].flui.VisOM(tH);
        }
        if (bo1 < 1e-15)
            bo1 = 1e-15;

        double rhog0;
        double rhogST0;
        double dg0G;
        double yco20G;
        if ((i > 0 || state.input.ConContEntrada == 1) && state.cells[i].QG > 0) {
            double pe;
            double te;
            if (i > 0) {
                pe = state.cells[i - 1].pres;
                te = state.cells[i - 1].temp;
            } else {
                pe = state.inletPressure;
                te = state.inletTemperature;
            }
            rhog0 = state.cells[i].rgL;
            rhogST0 = (*state.cells[i].fluiL).Deng * 1.225;
            dg0G = (*state.cells[i].fluiL).Deng;
            yco20G = (*state.cells[i].fluiL).yco2;
        } else {
            rhog0 = state.cells[i].rgC;
            rhogST0 = state.cells[i].flui.Deng * 1.225;
            dg0G = state.cells[i].flui.Deng;
            yco20G = state.cells[i].flui.yco2;
        }

        double rhog1 = state.cells[i].rgC;
        double rhogST1 = state.cells[i].flui.Deng * 1.225;
        double dg1G = state.cells[i].flui.Deng;
        double yco21G = state.cells[i].flui.yco2;
        if (state.cells[i + 1].QG <= 0.) {
            rhog1 = state.cells[i].rgR;
            rhogST1 = state.cells[i + 1].flui.Deng * 1.225;
            dg1G = state.cells[i + 1].flui.Deng;
            yco21G = state.cells[i + 1].flui.yco2;
        }

        if (i == 237) {
            int para;
            para = 0;
        }

        double hol = 1. - state.cells[i].alf;
        double bet = state.cells[i].bet;
        double rholST = (1 - state.cells[i].flui.BSW) * (1000 * 141.5 / (131.5 + state.cells[i].flui.API)) + state.cells[i].flui.BSW * 1000 * state.cells[i].flui.Denag;
        double rhog = state.cells[i].rgC;
        double rhogST = state.cells[i].flui.Deng * 1.225;
        double rs = state.cells[i].flui.RS(state.cells[i].pres, state.cells[i].temp);
        double razdgd = 1 / state.cells[i].flui.rDgD;
        double razdgl = 1 / state.cells[i].flui.rDgL;
        double bo = state.cells[i].flui.BOFunc(state.cells[i].pres, state.cells[i].temp, rs);
        double ba = state.cells[i].flui.BAFunc(state.cells[i].pres, state.cells[i].temp);
        double bsw = state.cells[i].flui.BSW * ba / (bo + ba * state.cells[i].flui.BSW - state.cells[i].flui.BSW * bo);
        rs = rs * 6.29 / 35.31467;

        double dgini = state.cells[i].flui.Deng;
        double yco2ini = state.cells[i].flui.yco2;
        double rgoini = state.cells[i].flui.RGO;
        double fonteO = state.cells[i].fontemassLR;
        double fonteG = state.cells[i].fontemassGR;
        double fonteP = state.cells[i].fontemassLR;
        double fonteA = state.cells[i].fontemassLR;
        double APIini = state.cells[i].flui.API;
        double BSWini = state.cells[i].flui.BSW;
        double denagini = state.cells[i].flui.Denag;
        double APIF = APIini;
        double BSWF = BSWini;
        double denagF = denagini;
        double dgFO = dgini;
        double dgFG = dgini;
        double yco2FO = yco2ini;
        double yco2FG = yco2ini;
        double rgoFO = rgoini;
        double viscLini = 0 * 30 + 1 * state.cells[i].flui.VisOM(tL);
        double viscHini = 0 * 20 + 1 * state.cells[i].flui.VisOM(tH);
        double viscLF = viscLini;
        double viscHF = viscHini;
        double rholSTF = rholST;
        double rhogSTF = rhogST;
        double razdgdF = 1.;
        double razdglF = 1.;
        double titFonte = 0.;
        ProFlu fluF;
        if (state.cells[i].acsr.tipo == 1 && state.cells[i].acsr.injg.seco == 1) {
            if (state.cells[i].acsr.injg.QGas > 0.)
                fluF = state.cells[i].acsr.injg.FluidoPro;
            else
                fluF = state.cells[i].flui;
            dgFG = fluF.Deng;
            yco2FG = fluF.yco2;
            rhogSTF = fluF.Deng * 1.225;
        } else if (state.cells[i].acsr.tipo == 1 && state.cells[i].acsr.injg.seco == 0) {
            if (state.cells[i].acsr.injg.QGas > 0.)
                fluF = state.cells[i].acsr.injg.FluidoPro;
            else
                fluF = state.cells[i].flui;

            titFonte = fluF.dStockTankVaporMassFraction;
            dgFO = fluF.Deng;
            yco2FO = fluF.yco2;
            rgoFO = fluF.RGO;

            double rsF = fluF.RS(state.cells[i].pres, state.cells[i].temp) * (6.29 / 35.31467);
            razdgdF = 1 / fluF.rDgD;
            razdglF = 1 / fluF.rDgL;
            if (state.cells[i].acsr.injg.FluidoPro.BSW < 1 - (*state.globals).localtiny)
                rholSTF = (1 - fluF.BSW) * (1000 * 141.5 / (131.5 + fluF.API)) + fluF.BSW * 1000 * fluF.Denag + fluF.Deng * 1.225 * rgoFO * (1. - fluF.BSW);
            else
                rholSTF = fluF.BSW * 1000 * fluF.Denag;
            fonteO = (fonteO + fonteG) * razdgdF * (fluF.RS(state.cells[i].pres, state.cells[i].temp) * (6.29 / 35.31467) * (1. - fluF.BSW) / rholSTF);

            if (state.input.nfluP > 1 || (*state.globals).chaverede != 0) {
                double bswaux = fluF.BSW;
                double contrabsw = 1. - bswaux;
                if (contrabsw < (*state.globals).localtiny)
                    contrabsw = 0.9 * (*state.globals).localtiny;
                double rhoPSTF;
                if (contrabsw > (*state.globals).localtiny)
                    rhoPSTF = (1000 * 141.5 / (131.5 + fluF.API)) + (bswaux / contrabsw) * 1000 * fluF.Denag + fluF.Deng * 1.225 * rsF;
                else
                    rhoPSTF = 1000 * fluF.Denag;
                if (contrabsw > (*state.globals).localtiny) {
                    fonteP *= (1. / rhoPSTF);
                    fonteA *= ((bswaux / contrabsw) / rhoPSTF);
                } else {
                    fonteP = 0.;
                    fonteA *= (1 / (fluF.BSW * 1000 * fluF.Denag));
                }
                APIF = fluF.API;
                BSWF = fluF.BSW;
                denagF = state.cells[i].acsr.injg.FluidoPro.Denag;
                viscLF = 0 * 30 + 1 * fluF.VisOM(tL);
                viscHF = 0 * 20 + 1 * fluF.VisOM(tH);
            }
        } else if (state.cells[i].acsr.tipo == 2) {
            if (state.cells[i].acsr.injl.QLiq > 0.)
                fluF = state.cells[i].acsr.injl.FluidoPro;
            else
                fluF = state.cells[i].flui;

            titFonte = fluF.dStockTankVaporMassFraction;
            dgFO = fluF.Deng;
            yco2FO = fluF.yco2;
            rgoFO = fluF.RGO;

            double rsF = fluF.RS(state.cells[i].pres, state.cells[i].temp) * (6.29 / 35.31467);
            razdgdF = 1 / fluF.rDgD;
            razdglF = 1 / fluF.rDgL;
            if (fluF.BSW < 1 - (*state.globals).localtiny)
                rholSTF = (1 - fluF.BSW) * (1000 * 141.5 / (131.5 + fluF.API)) + fluF.BSW * 1000 * fluF.Denag + fluF.Deng * 1.225 * rgoFO * (1. - fluF.BSW);
            else
                rholSTF = fluF.BSW * 1000 * fluF.Denag;
            fonteO = (fonteO + fonteG) * razdgdF * (fluF.RS(state.cells[i].pres, state.cells[i].temp) * (6.29 / 35.31467) * (1. - fluF.BSW) / rholSTF);

            if (state.input.nfluP > 1 || (*state.globals).chaverede != 0) {
                double bswaux = fluF.BSW;
                double contrabsw = 1. - bswaux;
                if (contrabsw < (*state.globals).localtiny)
                    contrabsw = 0.9 * (*state.globals).localtiny;
                double rhoPSTF;
                if (contrabsw > (*state.globals).localtiny)
                    rhoPSTF = (1000 * 141.5 / (131.5 + fluF.API)) + (bswaux / contrabsw) * 1000 * fluF.Denag + fluF.Deng * 1.225 * rsF;
                else
                    rhoPSTF = 1000 * fluF.Denag;
                if (contrabsw > (*state.globals).localtiny) {
                    fonteP *= (1. / rhoPSTF);
                    fonteA *= ((bswaux / contrabsw) / rhoPSTF);
                } else {
                    fonteP = 0.;
                    fonteA *= (1 / (fluF.BSW * 1000 * fluF.Denag));
                }
                APIF = fluF.API;
                BSWF = fluF.BSW;
                denagF = fluF.Denag;
                viscLF = 0 * 30 + 1 * fluF.VisOM(tL);
                viscHF = 0 * 20 + 1 * fluF.VisOM(tH);
            }
        } else if (state.cells[i].acsr.tipo == 10) {
            if ((state.cells[i].acsr.injm.MassC + state.cells[i].acsr.injm.MassG + state.cells[i].acsr.injm.MassP) > 0.)
                fluF = state.cells[i].acsr.injm.FluidoPro;
            else
                fluF = state.cells[i].flui;

            titFonte = fluF.dStockTankVaporMassFraction;
            dgFO = fluF.Deng;
            yco2FO = fluF.yco2;
            rgoFO = fluF.RGO;

            double rsF = fluF.RS(state.cells[i].pres, state.cells[i].temp) * (6.29 / 35.31467);
            razdgdF = 1 / fluF.rDgD;
            razdglF = 1 / fluF.rDgL;
            rholSTF = (1 - fluF.BSW) * (1000 * 141.5 / (131.5 + fluF.API)) +
                      fluF.BSW * 1000 *
                          fluF.Denag +
                      fluF.Deng * 1.225 * rgoFO * (1. - fluF.BSW);
            fonteO = (fonteO + fonteG) * razdgdF * (fluF.RS(state.cells[i].pres, state.cells[i].temp) * (6.29 / 35.31467) * (1. - fluF.BSW) / rholSTF);

            if (state.input.nfluP > 1 || (*state.globals).chaverede != 0) {
                double bswaux = fluF.BSW;
                double contrabsw = 1. - bswaux;
                if (contrabsw < (*state.globals).localtiny)
                    contrabsw = 0.9 * (*state.globals).localtiny;
                double rhoPSTF;
                if (contrabsw > (*state.globals).localtiny)
                    rhoPSTF = (1000 * 141.5 / (131.5 + fluF.API)) +
                              (bswaux / contrabsw) * 1000 *
                                  fluF.Denag +
                              fluF.Deng * 1.225 * rsF;
                else
                    rhoPSTF = 1000 * fluF.Denag;
                if (contrabsw > (*state.globals).localtiny) {
                    fonteP *= (1. / rhoPSTF);
                    fonteA *= ((bswaux / contrabsw) / rhoPSTF);
                } else {
                    fonteP = 0.;
                    fonteA *= (1 / (fluF.BSW * 1000 * fluF.Denag));
                }
                APIF = fluF.API;
                BSWF = fluF.BSW;
                denagF = fluF.Denag;
                viscLF = 0 * 30 + 1 * fluF.VisOM(tL);
                viscHF = 0 * 20 + 1 * fluF.VisOM(tH);
            }
        } else if (state.cells[i].acsr.tipo == 3) {
            if ((state.cells[i].acsr.ipr.Pres) > state.cells[i].pres)
                fluF = state.cells[i].acsr.ipr.FluidoPro;
            else
                fluF = state.cells[i].flui;

            titFonte = fluF.dStockTankVaporMassFraction;
            dgFO = fluF.Deng;
            yco2FO = fluF.yco2;
            rgoFO = fluF.RGO;

            rholSTF = (1 - fluF.BSW) * (1000 * 141.5 / (131.5 + fluF.API)) + fluF.BSW * 1000 * fluF.Denag + fluF.Deng * 1.225 * rgoFO * (1. - fluF.BSW);
            double rsF = fluF.RS(state.cells[i].pres, state.cells[i].temp) * (6.29 / 35.31467);
            razdgdF = 1 / fluF.rDgD;
            razdglF = 1 / fluF.rDgL;
            fonteO = (fonteO + fonteG) * razdgdF * (fluF.RS(state.cells[i].pres, state.cells[i].temp) * (6.29 / 35.31467) * (1. - fluF.BSW) / rholSTF);
            if (state.input.nfluP > 1 || (*state.globals).chaverede != 0) {
                double bswaux = fluF.BSW;
                double contrabsw = 1. - bswaux;
                if (contrabsw <= (*state.globals).localtiny)
                    contrabsw = 0.9 * (*state.globals).localtiny;
                double rhoPSTF;
                if (contrabsw > 0)
                    rhoPSTF = (1000 * 141.5 / (131.5 + fluF.API)) + (bswaux / contrabsw) * 1000 * fluF.Denag + fluF.Deng * 1.225 * rsF;
                else
                    rhoPSTF = 1000 * fluF.Denag;
                if (contrabsw > (*state.globals).localtiny) {
                    fonteP *= (1. / rhoPSTF);
                    fonteA *= ((bswaux / contrabsw) / rhoPSTF);
                } else {
                    fonteP = 0.;
                    fonteA *= (1 / (fluF.BSW * 1000 * fluF.Denag));
                }
                APIF = fluF.API;
                BSWF = fluF.BSW;
                denagF = fluF.Denag;
                viscLF = 0 * 30 + 1 * fluF.VisOM(tL);
                viscHF = 0 * 20 + 1 * fluF.VisOM(tH);
            }
        } else if (state.cells[i].acsr.tipo == 9) {
            ProFlu fluF;
            if (state.cells[i].acsr.fontechk.presT > state.cells[i].acsr.fontechk.pamb) {
                fluF = state.cells[i].acsr.fontechk.fluidoP;
            } else {
                fluF = state.cells[i].acsr.fontechk.fluidoPamb;
            }
            dgFO = fluF.Deng;
            yco2FO = fluF.yco2;
            rgoFO = fluF.RGO;
            titFonte = fluF.dStockTankVaporMassFraction;

            if (fluF.BSW < 1 - (*state.globals).localtiny)
                rholSTF = (1 - fluF.BSW) * (1000 * 141.5 / (131.5 + fluF.API)) + fluF.BSW * 1000 * fluF.Denag + fluF.Deng * 1.225 * rgoFO * (1. - fluF.BSW);
            else
                rholSTF = fluF.BSW * 1000 * fluF.Denag;

            double rsF = fluF.RS(state.cells[i].pres, state.cells[i].temp) * (6.29 / 35.31467);
            razdgdF = 1 / fluF.rDgD;
            razdglF = 1 / fluF.rDgL;
            if (state.cells[i].acsr.fontechk.ambGas != 1 || (fonteO + fonteG) < 0.)
                fonteO = (fonteO + fonteG) * razdgdF * (fluF.RS(state.cells[i].pres, state.cells[i].temp) * (6.29 / 35.31467) * (1. - fluF.BSW) / rholSTF);
            else
                fonteO = 0.;
            if (state.input.nfluP > 1 || (*state.globals).chaverede != 0) {
                double bswaux = fluF.BSW;
                double contrabsw = 1. - bswaux;
                if (contrabsw <= (*state.globals).localtiny)
                    contrabsw = 0.9 * (*state.globals).localtiny;
                double rhoPSTF;
                if (contrabsw > (*state.globals).localtiny)
                    rhoPSTF = (1000 * 141.5 / (131.5 + fluF.API)) + (bswaux / contrabsw) * 1000 * fluF.Denag + fluF.Deng * 1.225 * rsF;
                else
                    rhoPSTF = 1000 * fluF.Denag;
                if (contrabsw > (*state.globals).localtiny) {
                    fonteP *= (1. / rhoPSTF);
                    fonteA *= ((bswaux / contrabsw) / rhoPSTF);
                } else {
                    fonteP = 0.;
                    fonteA *= (1 / (fluF.BSW * 1000 * fluF.Denag));
                }
                APIF = fluF.API;
                BSWF = fluF.BSW;
                denagF = fluF.Denag;
                viscLF = 0 * 30 + 1 * fluF.VisOM(tL);
                viscHF = 0 * 20 + 1 * fluF.VisOM(tH);
            }
        } else if (state.cells[i].acsr.tipo == 15) {
            if ((state.cells[i].fontemassLR + state.cells[i].fontemassGR) > 1e-15)
                fluF = state.cells[i].acsr.radialPoro.flup;
            else
                fluF = state.cells[i].flui;

            titFonte = fluF.dStockTankVaporMassFraction;
            dgFO = fluF.Deng;
            yco2FO = fluF.yco2;
            rgoFO = fluF.RGO;

            rholSTF = (1 - state.cells[i].acsr.radialPoro.BSW) * (1000 * 141.5 / (131.5 + fluF.API)) + state.cells[i].acsr.radialPoro.BSW * 1000 * fluF.Denag + fluF.Deng * 1.225 * rgoFO * (1. - state.cells[i].acsr.radialPoro.BSW);
            double rsF = fluF.RS(state.cells[i].pres, state.cells[i].temp) * (6.29 / 35.31467);
            razdgdF = 1 / fluF.rDgD;
            razdglF = 1 / fluF.rDgL;
            fonteO = (fonteO + fonteG) * razdgdF * (fluF.RS(state.cells[i].pres, state.cells[i].temp) * (6.29 / 35.31467) * (1. - state.cells[i].acsr.radialPoro.BSW) / rholSTF);
            if (state.input.nfluP > 1 || (*state.globals).chaverede != 0) {
                double bswaux = state.cells[i].acsr.radialPoro.BSW;
                double contrabsw = 1. - bswaux;
                if (contrabsw <= (*state.globals).localtiny)
                    contrabsw = 0.9 * (*state.globals).localtiny;
                double rhoPSTF;
                if (contrabsw > 0)
                    rhoPSTF = (1000 * 141.5 / (131.5 + fluF.API)) + (bswaux / contrabsw) * 1000 * fluF.Denag + fluF.Deng * 1.225 * rsF;
                else
                    rhoPSTF = 1000 * fluF.Denag;
                if (contrabsw > (*state.globals).localtiny) {
                    fonteP *= (1. / rhoPSTF);
                    fonteA *= ((bswaux / contrabsw) / rhoPSTF);
                } else {
                    fonteP = 0.;
                    fonteA *= (1 / (state.cells[i].acsr.radialPoro.BSW * 1000 * fluF.Denag));
                }
                APIF = fluF.API;
                BSWF = state.cells[i].acsr.radialPoro.BSW;
                denagF = fluF.Denag;
                viscLF = 0 * 30 + 1 * fluF.VisOM(tL);
                viscHF = 0 * 20 + 1 * fluF.VisOM(tH);
            }
        } else if (state.cells[i].acsr.tipo == 16) {
            if ((state.cells[i].fontemassLR + state.cells[i].fontemassGR) > 1e-15)
                fluF = state.cells[i].acsr.radialPoro.flup;
            else
                fluF = state.cells[i].flui;

            titFonte = fluF.dStockTankVaporMassFraction;
            dgFO = fluF.Deng;
            yco2FO = fluF.yco2;
            rgoFO = fluF.RGO;

            rholSTF = (1 - state.cells[i].acsr.poroso2D.dados.transfer.BSW) * (1000 * 141.5 / (131.5 + fluF.API)) + state.cells[i].acsr.poroso2D.dados.transfer.BSW * 1000 * fluF.Denag + fluF.Deng * 1.225 * rgoFO * (1. - state.cells[i].acsr.poroso2D.dados.transfer.BSW);
            double rsF = fluF.RS(state.cells[i].pres, state.cells[i].temp) * (6.29 / 35.31467);
            razdgdF = 1 / fluF.rDgD;
            razdglF = 1 / fluF.rDgL;
            fonteO = (fonteO + fonteG) * razdgdF * (fluF.RS(state.cells[i].pres, state.cells[i].temp) * (6.29 / 35.31467) * (1. - state.cells[i].acsr.poroso2D.dados.transfer.BSW) / rholSTF);
            if (state.input.nfluP > 1 || (*state.globals).chaverede != 0) {
                double bswaux = state.cells[i].acsr.poroso2D.dados.transfer.BSW;
                double contrabsw = 1. - bswaux;
                if (contrabsw <= (*state.globals).localtiny)
                    contrabsw = 0.9 * (*state.globals).localtiny;
                double rhoPSTF;
                if (contrabsw > 0)
                    rhoPSTF = (1000 * 141.5 / (131.5 + fluF.API)) + (bswaux / contrabsw) * 1000 * fluF.Denag + fluF.Deng * 1.225 * rsF;
                else
                    rhoPSTF = 1000 * fluF.Denag;
                if (contrabsw > (*state.globals).localtiny) {
                    fonteP *= (1. / rhoPSTF);
                    fonteA *= ((bswaux / contrabsw) / rhoPSTF);
                } else {
                    fonteP = 0.;
                    fonteA *= (1 / (state.cells[i].acsr.poroso2D.dados.transfer.BSW * 1000 * fluF.Denag));
                }
                APIF = fluF.API;
                BSWF = state.cells[i].acsr.poroso2D.dados.transfer.BSW;
                denagF = fluF.Denag;
                viscLF = 0 * 30 + 1 * fluF.VisOM(tL);
                viscHF = 0 * 20 + 1 * fluF.VisOM(tH);
            }
        } else if ((fabs(fonteO) > (*state.globals).localtiny && state.cells[i].acsr.tipo != 2 && state.cells[i].acsr.tipo != 3 &&
                    state.cells[i].acsr.tipo != 9 && state.cells[i].acsr.tipo != 15 && state.cells[i].acsr.tipo != 16) ||
                   (fabs(fonteG) > (*state.globals).localtiny && state.cells[i].acsr.tipo != 1 && state.cells[i].acsr.tipo != 2 && state.cells[i].acsr.tipo != 3 && state.cells[i].acsr.tipo != 9 && state.cells[i].acsr.tipo != 15 && state.cells[i].acsr.tipo != 16)) {
            if (state.cells[i].acsr.tipo == 5 || state.cells[i].acsr.tipo == 8) {
                dgFO = state.cells[i].flui.Deng;
                yco2FO = state.cells[i].flui.yco2;
                rgoFO = state.cells[i].flui.RGO;
                titFonte = state.cells[i].flui.dStockTankVaporMassFraction;
                double rholiq = (1 - state.cells[i].flui.BSW) * (1000 * 141.5 / (131.5 + state.cells[i].flui.API)) + state.cells[i].flui.BSW * 1000 * state.cells[i].flui.Denag;
                double rhogas = state.cells[i].flui.Deng * 1.225;
                rholSTF = rholiq + rhogas * rgoFO * (1. - state.cells[i].flui.BSW);
                rhogSTF = rhogas;
                double rsF = state.cells[i].flui.RS(state.cells[i].pres, state.cells[i].temp) * (6.29 / 35.31467);
                razdgdF = 1 / state.cells[i].flui.rDgD;
                razdglF = 1 / state.cells[i].flui.rDgL;
                fonteO = (fonteO + fonteG) * (razdgdF * rsF * (1. - state.cells[i].flui.BSW) / rholSTF);
                if (state.input.nfluP > 1 || (*state.globals).chaverede != 0) {
                    double rhoPSTF = (1 - state.cells[i].flui.BSW) * (1000 * 141.5 / (131.5 + state.cells[i].flui.API)) + state.cells[i].flui.BSW * 1000 * state.cells[i].flui.Denag;
                    fonteP *= ((1 - state.cells[i].flui.BSW) / rhoPSTF);
                    fonteA *= (state.cells[i].flui.BSW / rhoPSTF);
                    APIF = state.cells[i].flui.API;
                    BSWF = state.cells[i].flui.BSW;
                    viscLF = 0 * 30 + 1 * state.cells[i].flui.VisOM(tL);
                    viscHF = 0 * 20 + 1 * state.cells[i].flui.VisOM(tH);
                }
            } else if ((*state.cells[i].acsrL).tipo == 5 || (*state.cells[i].acsrL).tipo == 8) {
                double rholiq;
                double rhogas;

                if (i > 0) {
                    dgFO = state.cells[i - 1].flui.Deng;
                    yco2FO = state.cells[i - 1].flui.yco2;
                    rgoFO = state.cells[i - 1].flui.RGO;
                    rholiq = (1 - state.cells[i - 1].flui.BSW) * (1000 * 141.5 / (131.5 + state.cells[i - 1].flui.API)) +
                             state.cells[i - 1].flui.BSW * 1000 * state.cells[i - 1].flui.Denag;
                    rhogas = state.cells[i - 1].flui.Deng * 1.225;
                    titFonte = state.cells[i - 1].flui.dStockTankVaporMassFraction;
                    rholSTF = rholiq + rhogas * rgoFO * (1. - state.cells[i - 1].flui.BSW);
                    rhogSTF = rhogas;
                    double rsF = state.cells[i - 1].flui.RS(state.cells[i].pres, state.cells[i].temp) * (6.29 / 35.31467);
                    razdgdF = 1 / state.cells[i - 1].flui.rDgD;
                    razdglF = 1 / state.cells[i - 1].flui.rDgL;
                    fonteO = (fonteO + fonteG) * (razdgdF * rsF *
                                                  (1. - state.cells[i - 1].flui.BSW) / rholSTF);
                    if (state.input.nfluP > 1 || (*state.globals).chaverede != 0) {
                        double rhoPSTF = (1 - state.cells[i - 1].flui.BSW) * (1000 * 141.5 / (131.5 + state.cells[i - 1].flui.API)) +
                                         state.cells[i - 1].flui.BSW * 1000 *
                                             state.cells[i - 1].flui.Denag;
                        fonteP *= ((1 - state.cells[i - 1].flui.BSW) / rhoPSTF);
                        fonteA *= (state.cells[i - 1].flui.BSW / rhoPSTF);
                        APIF = state.cells[i - 1].flui.API;
                        BSWF = state.cells[i - 1].flui.BSW;
                        denagF = state.cells[i - 1].flui.Denag;
                        viscLF = 0 * 30 + 1 * state.cells[i - 1].flui.VisOM(tL);
                        viscHF = 0 * 20 + 1 * state.cells[i - 1].flui.VisOM(tH);
                    }
                } else {
                    dgFO = state.cells[i].flui.Deng;
                    yco2FO = state.cells[i].flui.yco2;
                    rgoFO = state.cells[i].flui.RGO;
                    rholiq = (1 - state.cells[i].flui.BSW) * (1000 * 141.5 / (131.5 + state.cells[i].flui.API)) +
                             state.cells[i].flui.BSW * 1000 * state.cells[i].flui.Denag;
                    rhogas = state.cells[i].flui.Deng * 1.225;
                    titFonte = state.cells[i].flui.dStockTankVaporMassFraction;
                    rholSTF = rholiq + rhogas * rgoFO * (1. - state.cells[i].flui.BSW);
                    rhogSTF = rhogas;
                    double rsF = state.cells[i].flui.RS(state.cells[i].pres, state.cells[i].temp) * (6.29 / 35.31467);
                    razdgdF = 1 / state.cells[i].flui.rDgD;
                    razdglF = 1 / state.cells[i].flui.rDgL;
                    fonteO = (fonteO + fonteG) * (razdgdF * state.cells[i].flui.RS(state.cells[i].pres, state.cells[i].temp) * (6.29 / 35.31467) *
                                                  (1. - state.cells[i].flui.BSW) / rholSTF);
                    if (state.input.nfluP > 1 || (*state.globals).chaverede != 0) {
                        double rhoPSTF = (1 - state.cells[i].flui.BSW) * (1000 * 141.5 / (131.5 + state.cells[i].flui.API)) +
                                         state.cells[i].flui.BSW * 1000 *
                                             state.cells[i].flui.Denag;
                        fonteP *= ((1 - state.cells[i].flui.BSW) / rhoPSTF);
                        fonteA *= (state.cells[i].flui.BSW / rhoPSTF);
                        APIF = state.cells[i].flui.API;
                        BSWF = state.cells[i].flui.BSW;
                        denagF = state.cells[i].flui.Denag;
                        viscLF = 0 * 30 + 1 * state.cells[i].flui.VisOM(tL);
                        viscHF = 0 * 20 + 1 * state.cells[i].flui.VisOM(tH);
                    }
                }
            }
        }
        fonteG *= (razdglF / (rhogSTF));

        if (titFonte > 1. - 1e-15) {
            fonteO = 0.;
            fonteP = 0.;
            fonteA = 0.;
        }

        MultOe = 0.;
        if (state.cells[i].flui.dStockTankVaporMassFraction < 1. - 1e-15)
            MultOe = state.cells[i].QL * (1 - betI0) * (1 - bsw0) * razdgd0 * rs0 / bo0;
        MultOd = 0.;
        if (state.cells[i + 1].flui.dStockTankVaporMassFraction < 1. - 1e-15)
            MultOd = state.cells[i + 1].QL * (1 - betI1) * (1 - bsw1) * razdgd1 * rs1 / bo1;
        double MultGe = (state.cells[i].MC - state.cells[i].Mliqini) * razdgl0 / (rhogST0);
        double MultGd = (state.cells[i + 1].MC - state.cells[i + 1].Mliqini) * razdgl1 / (rhogST1);
        double volleveFim = (((1 - hol) * rhog * razdgl / (rhogST)) + hol * (1 - bet) * (1. - bsw) * rs * razdgd / (bo));
        if (volleveFim < 1e-15)
            volleveFim = 0.;
        double residuo = (volleveFim - state.cells[i].VolLeveST) * at / dt + (MultOd - MultOe) / dx + (MultGd - MultGe) / dx - (fonteO / dx + fonteG / dx);
        double volpesFim = hol * (1 - bet) * (1 - bsw) / bo;
        double volaguaFim = hol * (1 - bet) * bsw; // nao deveria ser dividido por Bo???????????
        double MultPe = 0.;
        double MultPd = 0.;
        double residuoP = 0.;
        double MultAe = 0.;
        double MultAd = 0.;
        double residuoA;
        if (state.input.nfluP > 1 || (*state.globals).chaverede != 0) {
            MultPe = 0.;
            if (state.cells[i].flui.dStockTankVaporMassFraction < 1. - 1e-15)
                MultPe = state.cells[i].QL * (1 - betI0) * (1 - bsw0) / bo0;
            MultPd = 0.;
            if (state.cells[i + 1].flui.dStockTankVaporMassFraction < 1. - 1e-15)
                MultPd = state.cells[i + 1].QL * (1 - betI1) * (1 - bsw1) / bo1;
            residuoP = (volpesFim - state.cells[i].VolPesaST) * at / dt + (MultPd - MultPe) / dx - fonteP / dx;
            MultAe = 0.;
            if (state.cells[i].flui.dStockTankVaporMassFraction < 1. - 1e-15)
                MultAe = state.cells[i].QL * (1 - betI0) * bsw0 / bo0;
            MultAd = 0.;
            if (state.cells[i + 1].flui.dStockTankVaporMassFraction < 1. - 1e-15)
                MultAd = state.cells[i + 1].QL * (1 - betI1) * bsw1 / bo1;
            residuoA = (volaguaFim - state.cells[i].VolAguaST) * at / dt + (MultAd - MultAe) / dx - fonteA / dx;
        }
        rgo[i] = (*state.globals).RGOMax;
        if (hol > (*state.globals).localtiny && bet < (1. - (*state.globals).localtiny) && bsw < (1. - (*state.globals).localtiny)) {
            rgo[i] = (volleveFim - residuo * dt / at) * bo / (hol * (1 - bet) * (1 - bsw));
            if (rgo[i] > (*state.globals).RGOMax)
                rgo[i] = (*state.globals).RGOMax;
        } else if (bet >= (1. - (*state.globals).localtiny) || bsw >= (1. - (*state.globals).localtiny))
            rgo[i] = 0.;
        else
            rgo[i] = (*state.globals).RGOMax;

        if (volleveFim > 1e-5 && state.input.flashCompleto == 0 && ((fonteG >= 0 || fonteO > 0) || ((MultGd < 0 || MultGe > 0) || (MultOd < 0 || MultOe > 0)))) {
            dg[i] = (dt * (dgFO * fonteO / dx + dgFG * fonteG / dx + 1. * dgini * residuo - (dg1O * MultOd - dg0O * MultOe) / dx - (dg1G * MultGd - dg0G * MultGe) / dx) + dgini * state.cells[i].VolLeveST * at) /
                    (volleveFim * at - 0. * residuo * dt);
            yco2[i] = (dt * (yco2FO * fonteO / dx + yco2FG * fonteG / dx + 1. * yco2ini * residuo - (yco21O * MultOd - yco20O * MultOe) / dx - (yco21G * MultGd - yco20G * MultGe) / dx) + yco2ini * state.cells[i].VolLeveST * at) /
                      (volleveFim * at - 0. * residuo * dt);
            if (yco2[i] < 0.)
                yco2[i] = 0.;
            else if (yco2[i] > 1.)
                yco2[i] = 1.;
        } else {
            dg[i] = dgini;
            yco2[i] = yco2ini;
        }
        if ((state.input.nfluP > 1 && state.input.flashCompleto == 0) || (*state.globals).chaverede != 0) {
            if (volpesFim > 1e-3 && (fonteP > 0 || (MultPd < 0 || MultPe > 0))) {
                double denmixSTDF = 141.5 / (131.5 + APIF);
                double denmixSTDini = 141.5 / (131.5 + APIini);
                double denmixSTD1 = 141.5 / (131.5 + API1);
                double denmixSTD0 = 141.5 / (131.5 + API0);
                API[i] = (dt * (denmixSTDF * fonteP / dx + 1. * denmixSTDini * residuoP - (denmixSTD1 * MultPd - denmixSTD0 * MultPe) / dx) + denmixSTDini * state.cells[i].VolPesaST * at) / (volpesFim * at - 0. * residuoP * dt);
                API[i] = 141.5 / API[i] - 131.5;
                VISCL[i] = (dt * (viscLF * fonteP / dx + 1. * viscLini * residuoP - (viscL1 * MultPd - viscL0 * MultPe) / dx) + viscLini * state.cells[i].VolPesaST * at) / (volpesFim * at - 0. * residuoP * dt);
                VISCH[i] = (dt * (viscHF * fonteP / dx + 1. * viscHini * residuoP - (viscH1 * MultPd - viscH0 * MultPe) / dx) + viscHini * state.cells[i].VolPesaST * at) / (volpesFim * at - 0. * residuoP * dt);
            } else {
                API[i] = APIini;
                VISCL[i] = viscLini;
                VISCH[i] = viscHini;
            }
        }
        if (state.input.nfluP > 1 || (*state.globals).chaverede != 0) {
            if ((volaguaFim + volpesFim) > 1e-3 && ((fonteP > 0 || fonteA > 0) ||
                                                    ((MultPd < 0 || MultPe > 0) || (MultAd < 0 || MultAe > 0)))) {
                BSW[i] = (dt * (BSWF * (fonteA + fonteP) / dx + 1. * BSWini * (residuoA + residuoP) - (BSW1 * (MultAd + MultPd) - BSW0 * (MultAe + MultPe)) / dx) + BSWini * (state.cells[i].VolAguaST + state.cells[i].VolPesaST) * at) /
                         ((volaguaFim + volpesFim) * at - 0. * (residuoA + residuoP) * dt);
                denag[i] = (dt * (denagF * (fonteA) / dx + 1. * denagini * (residuoA) - (denag1 * (MultAd)-denag0 * (MultAe)) / dx) + denagini * (state.cells[i].VolAguaST) * at) /
                           ((volaguaFim)*at - 0. * (residuoA)*dt);
                if (BSW[i] < 0.)
                    BSW[i] = 0.;
                else if (BSW[i] > 1.)
                    BSW[i] = 1.;
                if (denag[i] < 1.)
                    denag[i] = 1.;
            } else
                BSW[i] = BSWini;
            denag[i] = denagini;
        }
        state.cells[i].VolLeveST = volleveFim;
        state.cells[i].VolPesaST = volpesFim;
        state.cells[i].VolAguaST = volaguaFim;
    }
    for (int i = imin; i <= state.lastCell - 1; i++) {
        if (state.trackGasOilRatio > 0 && state.cells[i].flui.corrSat != 4)
            state.cells[i].flui.RGO = rgo[i];
        if (state.trackGasGravity > 0 && state.input.flashCompleto == 0) {
            state.cells[i].flui.Deng = dg[i];
            if ((state.cells[i].flui.Deng > 5 || state.cells[i].flui.Deng < 0) && i > 0)
                state.cells[i].flui.Deng = state.cells[i - 1].flui.Deng;
            state.cells[i].flui.yco2 = yco2[i];
        }
        if (state.input.nfluP > 1 || (*state.globals).chaverede != 0) {
            state.cells[i].flui.BSW = BSW[i];
            state.cells[i].flui.Denag = denag[i];
            if (state.input.flashCompleto == 0) {
                state.cells[i].flui.API = API[i];
                state.cells[i].flui.LVisL = VISCL[i];
                state.cells[i].flui.LVisH = VISCH[i];
                state.cells[i].flui.TempL = tL;
                state.cells[i].flui.TempH = tH;
            }
        }
        if (state.input.flashCompleto == 0)
            state.cells[i].flui.RenovaFluido();
        state.updaters.correctGasSpecificGravity(i);
    }

    if ((*state.globals).chaverede == 0 || state.endNode == 1 || state.cells[state.lastCell].Mliqini > -(*state.globals).localtiny) {
        if (state.trackGasOilRatio > 0 && state.cells[state.lastCell].flui.corrSat != 4)
            state.cells[state.lastCell].flui.RGO = rgo[state.lastCell - 1];
        if (state.trackGasGravity > 0 && state.input.flashCompleto == 0) {
            state.cells[state.lastCell].flui.Deng = dg[state.lastCell - 1];
            state.cells[state.lastCell].flui.yco2 = yco2[state.lastCell - 1];
        }
        if (state.input.nfluP > 1 || (*state.globals).chaverede != 0) {
            state.cells[state.lastCell].flui.BSW = BSW[state.lastCell - 1];
            state.cells[state.lastCell].flui.Denag = denag[state.lastCell - 1];
            if (state.input.flashCompleto == 0) {
                state.cells[state.lastCell].flui.API = API[state.lastCell - 1];
                state.cells[state.lastCell].flui.LVisL = VISCL[state.lastCell - 1];
                state.cells[state.lastCell].flui.LVisH = VISCH[state.lastCell - 1];
                state.cells[state.lastCell].flui.TempL = tL;
                state.cells[state.lastCell].flui.TempH = tH;
            }
        }
    } else {
        state.cells[state.lastCell].flui.RGO = fluiRev.RGO;
        state.cells[state.lastCell].flui.Deng = fluiRev.Deng;
        state.cells[state.lastCell].flui.yco2 = fluiRev.yco2;
        state.cells[state.lastCell].flui.API = fluiRev.API;
        state.cells[state.lastCell].flui.BSW = fluiRev.BSW;
        state.cells[state.lastCell].flui.Denag = fluiRev.Denag;
        state.cells[state.lastCell].flui.LVisL = fluiRev.LVisL;
        state.cells[state.lastCell].flui.LVisH = fluiRev.LVisH;
        state.cells[state.lastCell].flui.TempL = fluiRev.TempL;
        state.cells[state.lastCell].flui.TempH = fluiRev.TempH;
    }
    if (state.input.flashCompleto == 0)
        state.cells[state.lastCell].flui.RenovaFluido();
    if (state.input.corrDeng == 0) {
        state.cells[state.lastCell].flui.rDgD = 1.;
        state.cells[state.lastCell].flui.rDgL = 1.;
        state.cells[state.lastCell].flui.PCis = state.cells[state.lastCell].flui.PC;
        state.cells[state.lastCell].flui.TCis = state.cells[state.lastCell].flui.TC;
    } else {
        state.cells[state.lastCell].flui.razDegD(state.cells[state.lastCell].pres, state.cells[state.lastCell].temp);
        state.cells[state.lastCell].flui.rzDegL(state.cells[state.lastCell].pres, state.cells[state.lastCell].temp);
        state.cells[state.lastCell].flui.PcTcIS();
    }
}

void transportOverallMolarFractions(const CompositionState &state, ProFlu fluiRev) {
    Vcr<double> BSW(state.lastCell);
    Vcr<double> VISCL(state.lastCell);
    Vcr<double> VISCH(state.lastCell);
    int ncomp = state.input.npseudo;
    ProFlu *fluC;
    fluC = new ProFlu[state.lastCell];
    Vcr<double> fracMol0(ncomp);
    Vcr<double> fracMol1(ncomp);
    Vcr<double> fracMolF(ncomp);
    double pesoMol0;
    double pesoMol1;
    double pesoMolF;
    double pesoMolC;
    double dt = state.cells[1].dt;
    double tL = 0.;
    double tH = 70.;
    if (state.cells[0].acsr.tipo == 15 && (state.cells[0].fontemassCR + state.cells[0].fontemassGR + state.cells[0].fontemassLR) > 0.)
        state.cells[0].flui.BSW = state.cells[0].acsr.radialPoro.BSW;
    else if (state.cells[0].acsr.tipo == 16 && (state.cells[0].fontemassCR + state.cells[0].fontemassGR + state.cells[0].fontemassLR) > 0.)
        state.cells[0].flui.BSW = state.cells[0].acsr.poroso2D.dados.transfer.BSW;
    int imin = 1;
    if ((*state.globals).chaverede != 0 && (state.cells[0].acsr.tipo == 10 || state.input.ConContEntrada > 0))
        imin = 0;
    for (int i = imin; i < state.lastCell; i++) {
        double at = state.cells[i].duto.area;
        double dx = state.cells[i].dx;
        double ti0;
        if (state.cells[i].VTemper < 0.)
            ti0 = state.cells[i].temp;
        else {
            if (i > 0)
                ti0 = state.cells[i - 1].temp;
            else if (state.input.ConContEntrada == 1)
                ti0 = state.inletTemperature;
            else
                ti0 = state.cells[i].temp;
        }
        double ti1 = state.cells[i].temp;
        if (state.cells[i + 1].VTemper < 0.)
            ti1 = state.cells[i + 1].temp;

        double titF = 0.;
        ProFlu fluF;

        double boF = 1.;
        double baF = 1.;
        double fwF = 1.;
        double rhoOF = 900.;
        double rhoWF = 1000.;

        double betI0;
        double bo0;
        double ba0;
        double bsw0;
        double BSW0;
        double viscL0;
        double viscH0;
        if (i > 0 || state.input.ConContEntrada == 0) {
            if (i > 0 && state.cells[i].QG >= 0.)
                betI0 = state.cells[i - 1].betPigD;
            else
                betI0 = state.cells[i].betPigE;
        } else {
            if (state.cells[i].QG >= 0.)
                betI0 = state.inletCompletionFraction;
            else
                betI0 = state.cells[i].betPigE;
        }
        if ((i > 0 || state.input.ConContEntrada == 1) && state.cells[i].QL >= 0.) {
            double pe;
            double te;
            if (i > 0) {
                pe = state.cells[i - 1].pres;
                te = state.cells[i - 1].temp;
            } else {
                pe = state.inletPressure;
                te = state.inletTemperature;
            }
            if (state.input.ConContEntrada == 0)
                betI0 = state.cells[i - 1].betPigD; // testeBeta
            else
                betI0 = state.inletCompletionFraction; // testeBeta
            double rs0 = (*state.cells[i].fluiL).RS(pe, te);
            bo0 = (*state.cells[i].fluiL).BOFunc(pe, te, rs0);
            ba0 = (*state.cells[i].fluiL).BAFunc(pe, te);
            bsw0 = (*state.cells[i].fluiL).BSW * ba0 / (bo0 + ba0 * (*state.cells[i].fluiL).BSW - (*state.cells[i].fluiL).BSW * bo0);
            BSW0 = (*state.cells[i].fluiL).BSW;
            viscL0 = 0 * 30 + 1 * (*state.cells[i].fluiL).VisOM(tL);
            viscH0 = 0 * 20 + 1 * (*state.cells[i].fluiL).VisOM(tH);
        } else {
            betI0 = state.cells[i].betPigE; // testebeta
            double rs0 = state.cells[i].flui.RS(state.cells[i].pres, state.cells[i].temp);
            bo0 = state.cells[i].flui.BOFunc(state.cells[i].pres, state.cells[i].temp, rs0);
            ba0 = state.cells[i].flui.BAFunc(state.cells[i].pres, state.cells[i].temp);
            bsw0 = state.cells[i].flui.BSW * ba0 / (bo0 + ba0 * state.cells[i].flui.BSW - state.cells[i].flui.BSW * bo0);
            BSW0 = state.cells[i].flui.BSW;
            viscL0 = 0 * 30 + 1 * state.cells[i].flui.VisOM(tL);
            viscH0 = 0 * 20 + 1 * state.cells[i].flui.VisOM(tH);
        }
        if (bo0 < 1e-15)
            bo0 = 1e-15;

        double betI1 = state.cells[i].betPigD;
        double rs1 = state.cells[i].flui.RS(state.cells[i].pres, state.cells[i].temp);
        double bo1 = state.cells[i].flui.BOFunc(state.cells[i].pres, state.cells[i].temp, rs1);
        double ba1 = state.cells[i].flui.BAFunc(state.cells[i].pres, state.cells[i].temp);
        double bsw1 = state.cells[i].flui.BSW * ba1 / (bo1 + ba1 * state.cells[i].flui.BSW - state.cells[i].flui.BSW * bo1);
        double BSW1 = state.cells[i].flui.BSW;
        double viscL1 = 0 * 30 + 1 * state.cells[i].flui.VisOM(tL);
        double viscH1 = 0 * 20 + 1 * state.cells[i].flui.VisOM(tH);

        // betI1 = celula[i + 1].betPigE;    //duvidabeta
        if (state.cells[i + 1].QL < 0.) {
            betI1 = state.cells[i + 1].betPigE; // testeBet
            rs1 = state.cells[i + 1].flui.RS(state.cells[i + 1].pres, state.cells[i + 1].temp);
            bo1 = state.cells[i + 1].flui.BOFunc(state.cells[i + 1].pres, state.cells[i + 1].temp, rs1);
            ba1 = state.cells[i + 1].flui.BAFunc(state.cells[i + 1].pres, state.cells[i + 1].temp);
            bsw1 = state.cells[i + 1].flui.BSW * ba1 / (bo1 + ba1 * state.cells[i + 1].flui.BSW - state.cells[i + 1].flui.BSW * bo1);
            BSW1 = state.cells[i + 1].flui.BSW;
            viscL1 = 0 * 30 + 1 * state.cells[i + 1].flui.VisOM(tL);
            viscH1 = 0 * 20 + 1 * state.cells[i + 1].flui.VisOM(tH);
        }
        if (bo1 < 1e-15)
            bo1 = 1e-15;
        double titV0;
        double titV1;
        double vazMasLiq0 = state.cells[i].MliqiniL;
        double vazMasLiq1 = state.cells[i + 1].MliqiniL;
        double vazMasGas0 = state.cells[i].MC - state.cells[i].MliqiniL;
        double vazMasGas1 = state.cells[i + 1].MC - state.cells[i + 1].MliqiniL;
        pesoMolC = 0.;
        for (int kfrac = 0; kfrac < ncomp; kfrac++) {
            pesoMolC += state.cells[i].flui.masMol[kfrac] * state.cells[i].flui.fracMol[kfrac];
        }
        if ((i > 0 || state.input.ConContEntrada == 1) && state.cells[i].MC >= 0.) {
            double pe;
            double te;
            if (i > 0) {
                pe = state.cells[i - 1].pres;
                te = state.cells[i - 1].temp;
            } else {
                pe = state.inletPressure;
                te = state.inletTemperature;
            }
            double fwV = bsw0;
            double rhoOV = (*state.cells[i].fluiL).MasEspoleo(pe, te);
            double rhoWV = (*state.cells[i].fluiL).MasEspAgua(pe, te);
            titV0 = (1 - fwV) * rhoOV / ((1 - fwV) * rhoOV + fwV * rhoWV);
            vazMasLiq0 -= betI0 * state.cells[i].QL * state.cells[i].fluicol.MasEspFlu(pe, te);
            vazMasLiq0 *= titV0;
            pesoMol0 = 0;
            for (int kfrac = 0; kfrac < ncomp; kfrac++) {
                fracMol0[kfrac] = (*state.cells[i].fluiL).fracMol[kfrac];
                pesoMol0 += (*state.cells[i].fluiL).masMol[kfrac] * (*state.cells[i].fluiL).fracMol[kfrac];
            }
        } else {
            double fwV = bsw0;
            double rhoOV = state.cells[i].flui.MasEspoleo(state.cells[i].pres, state.cells[i].temp);
            double rhoWV = state.cells[i].flui.MasEspAgua(state.cells[i].pres, state.cells[i].temp);
            titV0 = (1 - fwV) * rhoOV / ((1 - fwV) * rhoOV + fwV * rhoWV);
            vazMasLiq0 -= betI0 * state.cells[i].QL * state.cells[i].fluicol.MasEspFlu(state.cells[i].pres, state.cells[i].temp);
            vazMasLiq0 *= titV0;
            pesoMol0 = 0;
            for (int kfrac = 0; kfrac < ncomp; kfrac++) {
                fracMol0[kfrac] = state.cells[i].flui.fracMol[kfrac];
                pesoMol0 += state.cells[i].flui.masMol[kfrac] * state.cells[i].flui.fracMol[kfrac];
            }
        }
        if (state.cells[i + 1].MC >= 0.) {
            double fwV = bsw1;
            double rhoOV = state.cells[i].flui.MasEspoleo(state.cells[i].pres, state.cells[i].temp);
            double rhoWV = state.cells[i].flui.MasEspAgua(state.cells[i].pres, state.cells[i].temp);
            titV1 = (1 - fwV) * rhoOV / ((1 - fwV) * rhoOV + fwV * rhoWV);
            vazMasLiq1 -= betI1 * state.cells[i].QL * state.cells[i].fluicol.MasEspFlu(state.cells[i].pres, state.cells[i].temp);
            vazMasLiq1 *= titV1;
            pesoMol1 = 0;
            for (int kfrac = 0; kfrac < ncomp; kfrac++) {
                fracMol1[kfrac] = state.cells[i].flui.fracMol[kfrac];
                pesoMol1 += state.cells[i].flui.masMol[kfrac] * state.cells[i].flui.fracMol[kfrac];
            }
        } else {
            double fwV = bsw1;
            double rhoOV = state.cells[i + 1].flui.MasEspoleo(state.cells[i + 1].pres, state.cells[i + 1].temp);
            double rhoWV = state.cells[i + 1].flui.MasEspAgua(state.cells[i + 1].pres, state.cells[i + 1].temp);
            titV1 = (1 - fwV) * rhoOV / ((1 - fwV) * rhoOV + fwV * rhoWV);
            vazMasLiq1 -= betI1 * state.cells[i + 1].QL * state.cells[i].fluicol.MasEspFlu(state.cells[i + 1].pres, state.cells[i + 1].temp);
            vazMasLiq1 *= titV1;
            pesoMol1 = 0;
            for (int kfrac = 0; kfrac < ncomp; kfrac++) {
                fracMol1[kfrac] = state.cells[i + 1].flui.fracMol[kfrac];
                pesoMol1 += state.cells[i + 1].flui.masMol[kfrac] * state.cells[i + 1].flui.fracMol[kfrac];
            }
        }
        double fonteO = state.cells[i].fontemassLR;
        double fonteG = state.cells[i].fontemassGR;
        double fonteP = state.cells[i].fontemassLR;
        double fonteA = state.cells[i].fontemassLR;

        if (state.cells[i].acsr.tipo == 1) {
            state.cells[i].acsr.injg.FluidoPro.atualizaPropComp(state.cells[i].pres, state.cells[i].temp, state.cells[i].flui.dCalculatedBeta,
                                                           state.cells[i].flui.oCalculatedLiqComposition, state.cells[i].flui.oCalculatedVapComposition, state.input.pocinjec);
            if (state.cells[i].acsr.injg.QGas > 0.)
                fluF = state.cells[i].acsr.injg.FluidoPro;
            else
                fluF = state.cells[i].flui;
            fwF = 0.;
            titF = 1.;
        } else if (state.cells[i].acsr.tipo == 2) {
            state.cells[i].acsr.injl.FluidoPro.atualizaPropComp(state.cells[i].pres, state.cells[i].temp, state.cells[i].flui.dCalculatedBeta,
                                                           state.cells[i].flui.oCalculatedLiqComposition, state.cells[i].flui.oCalculatedVapComposition, state.input.pocinjec);
            if (state.cells[i].acsr.injl.QLiq > 0.)
                fluF = state.cells[i].acsr.injl.FluidoPro;
            else
                fluF = state.cells[i].flui;
            boF = fluF.BOFunc(state.cells[i].pres, state.cells[i].temp);
            baF = fluF.BAFunc(state.cells[i].pres, state.cells[i].temp);
            fwF = fluF.BSW * baF / (boF + baF * fluF.BSW - fluF.BSW * boF);
            rhoOF = fluF.MasEspoleo(state.cells[i].pres, state.cells[i].temp);
            rhoWF = fluF.MasEspAgua(state.cells[i].pres, state.cells[i].temp);
            titF = (1 - fwF) * rhoOF / ((1 - fwF) * rhoOF + fwF * rhoWF);
        } else if (state.cells[i].acsr.tipo == 3) {
            state.cells[i].acsr.ipr.FluidoPro.atualizaPropComp(state.cells[i].pres, state.cells[i].temp, state.cells[i].flui.dCalculatedBeta,
                                                          state.cells[i].flui.oCalculatedLiqComposition, state.cells[i].flui.oCalculatedVapComposition, state.input.pocinjec);
            if ((state.cells[i].acsr.ipr.Pres) > state.cells[i].pres)
                fluF = state.cells[i].acsr.ipr.FluidoPro;
            else
                fluF = state.cells[i].flui;
            boF = fluF.BOFunc(state.cells[i].pres, state.cells[i].temp);
            baF = fluF.BAFunc(state.cells[i].pres, state.cells[i].temp);
            fwF = fluF.BSW * baF / (boF + baF * fluF.BSW - fluF.BSW * boF);
            rhoOF = fluF.MasEspoleo(state.cells[i].pres, state.cells[i].temp);
            rhoWF = fluF.MasEspAgua(state.cells[i].pres, state.cells[i].temp);
            titF = (1 - fwF) * rhoOF / ((1 - fwF) * rhoOF + fwF * rhoWF);
        } else if (state.cells[i].acsr.tipo == 10) {
            state.cells[i].acsr.injm.FluidoPro.atualizaPropComp(state.cells[i].pres, state.cells[i].temp, state.cells[i].flui.dCalculatedBeta,
                                                           state.cells[i].flui.oCalculatedLiqComposition, state.cells[i].flui.oCalculatedVapComposition, state.input.pocinjec);
            if ((state.cells[i].acsr.injm.MassC + state.cells[i].acsr.injm.MassG + state.cells[i].acsr.injm.MassP) > 0.)
                fluF = state.cells[i].acsr.injm.FluidoPro;
            else
                fluF = state.cells[i].flui;
            boF = fluF.BOFunc(state.cells[i].pres, state.cells[i].temp);
            baF = fluF.BAFunc(state.cells[i].pres, state.cells[i].temp);
            fwF = fluF.BSW * baF / (boF + baF * fluF.BSW - fluF.BSW * boF);
            rhoOF = fluF.MasEspoleo(state.cells[i].pres, state.cells[i].temp);
            rhoWF = fluF.MasEspAgua(state.cells[i].pres, state.cells[i].temp);
            titF = (1 - fwF) * rhoOF / ((1 - fwF) * rhoOF + fwF * rhoWF);
        } else if (state.cells[i].acsr.tipo == 9 && state.cells[i].acsr.fontechk.abertura > 1e-6) {
            state.cells[i].acsr.fontechk.fluidoP.atualizaPropComp(state.cells[i].pres, state.cells[i].temp, state.cells[i].flui.dCalculatedBeta,
                                                             state.cells[i].flui.oCalculatedLiqComposition, state.cells[i].flui.oCalculatedVapComposition, state.input.pocinjec);
            state.cells[i].acsr.fontechk.fluidoPamb.atualizaPropComp(state.cells[i].pres, state.cells[i].temp, state.cells[i].flui.dCalculatedBeta,
                                                                state.cells[i].flui.oCalculatedLiqComposition, state.cells[i].flui.oCalculatedVapComposition, state.input.pocinjec);
            if (state.cells[i].acsr.fontechk.presT > state.cells[i].acsr.fontechk.pamb) {
                fluF = state.cells[i].acsr.fontechk.fluidoP;
            } else {
                fluF = state.cells[i].acsr.fontechk.fluidoPamb;
            }
            boF = fluF.BOFunc(state.cells[i].pres, state.cells[i].temp);
            baF = fluF.BAFunc(state.cells[i].pres, state.cells[i].temp);
            fwF = fluF.BSW * baF / (boF + baF * fluF.BSW - fluF.BSW * boF);
            rhoOF = fluF.MasEspoleo(state.cells[i].pres, state.cells[i].temp);
            rhoWF = fluF.MasEspAgua(state.cells[i].pres, state.cells[i].temp);
            titF = (1 - fwF) * rhoOF / ((1 - fwF) * rhoOF + fwF * rhoWF);
        } else if (state.cells[i].acsr.tipo == 15) {
            double tRes = state.cells[i].acsr.radialPoro.tRes;
            state.cells[i].acsr.radialPoro.flup.atualizaPropComp(state.cells[i].pres, tRes, state.cells[i].flui.dCalculatedBeta,
                                                            state.cells[i].flui.oCalculatedLiqComposition, state.cells[i].flui.oCalculatedVapComposition, state.input.pocinjec);
            if ((state.cells[i].fontemassLR + state.cells[i].fontemassGR) > 1e-15)
                fluF = state.cells[i].acsr.radialPoro.flup;
            else
                fluF = state.cells[i].flui;
            boF = fluF.BOFunc(state.cells[i].pres, state.cells[i].temp);
            baF = fluF.BAFunc(state.cells[i].pres, state.cells[i].temp);
            fwF = fluF.BSW * baF / (boF + baF * fluF.BSW - fluF.BSW * boF);
            rhoOF = fluF.MasEspoleo(state.cells[i].pres, state.cells[i].temp);
            rhoWF = fluF.MasEspAgua(state.cells[i].pres, state.cells[i].temp);
            titF = (1 - fwF) * rhoOF / ((1 - fwF) * rhoOF + fwF * rhoWF);
        } else if (state.cells[i].acsr.tipo == 16) {
            double tRes = state.cells[i].acsr.poroso2D.dados.tRes;
            state.cells[i].acsr.poroso2D.dados.flup.atualizaPropComp(state.cells[i].pres, tRes, state.cells[i].flui.dCalculatedBeta,
                                                                state.cells[i].flui.oCalculatedLiqComposition, state.cells[i].flui.oCalculatedVapComposition, state.input.pocinjec);
            if ((state.cells[i].fontemassLR + state.cells[i].fontemassGR) > 1e-15)
                fluF = state.cells[i].acsr.poroso2D.dados.flup;
            else
                fluF = state.cells[i].flui;
            boF = fluF.BOFunc(state.cells[i].pres, state.cells[i].temp);
            baF = fluF.BAFunc(state.cells[i].pres, state.cells[i].temp);
            fwF = fluF.BSW * baF / (boF + baF * fluF.BSW - fluF.BSW * boF);
            rhoOF = fluF.MasEspoleo(state.cells[i].pres, state.cells[i].temp);
            rhoWF = fluF.MasEspAgua(state.cells[i].pres, state.cells[i].temp);
            titF = (1 - fwF) * rhoOF / ((1 - fwF) * rhoOF + fwF * rhoWF);
        }

        pesoMolF = 0;
        for (int j = 0; j < ncomp; j++)
            pesoMolF += fluF.masMol[j] * fluF.fracMol[j];
        fonteO *= titF;
        double tempMol;
        fluC[i] = state.cells[i].flui;
        for (int corrige = 0; corrige < 2; corrige++) {
            tempMol = (fluC[i].MasEspoleo(state.cells[i].pres, state.cells[i].temp) * (1 - state.cells[i].alf) +
                       fluC[i].MasEspGas(state.cells[i].pres, state.cells[i].temp) *
                           state.cells[i].alf) *
                      state.cells[i].duto.area / pesoMolC;
            for (int kfrac = 0; kfrac < ncomp; kfrac++) {
                fluC[i].fracMol[kfrac] = (state.cells[i].nMolIni * state.cells[i].flui.fracMol[kfrac] +
                                          ((fonteO + fonteG) * fracMolF[kfrac] / pesoMolF -
                                           (vazMasLiq1 + vazMasGas1) * fracMol1[kfrac] / pesoMol1 -
                                           (vazMasLiq0 + vazMasGas0) * fracMol0[kfrac] / pesoMol0) *
                                              dt / dx) /
                                         tempMol;
            }
            pesoMolC = 0;
            for (int kfrac = 0; kfrac < ncomp; kfrac++)
                pesoMolC += fluC[i].masMol[kfrac] * state.cells[i].flui.fracMol[kfrac];
            fluC[i].atualizaPropCompStandard();
            fluC[i].atualizaPropComp(state.cells[i].pres, state.cells[i].temp, fluC[i].dCalculatedBeta,
                                     fluC[i].oCalculatedLiqComposition,
                                     fluC[i].oCalculatedVapComposition, state.input.pocinjec);
        }
        state.cells[i].nMol = (fluC[i].MasEspoleo(state.cells[i].pres, state.cells[i].temp) * (1 - state.cells[i].alf) +
                          fluC[i].MasEspGas(state.cells[i].pres, state.cells[i].temp) *
                              state.cells[i].alf) *
                         state.cells[i].duto.area / pesoMolC;
        fluC[i].Pmol = pesoMolC;
        double hol = 1. - state.cells[i].alf;
        double bet = state.cells[i].bet;

        double rs = state.cells[i].flui.RS(state.cells[i].pres, state.cells[i].temp);
        double bo = state.cells[i].flui.BOFunc(state.cells[i].pres, state.cells[i].temp, rs);
        double ba = state.cells[i].flui.BAFunc(state.cells[i].pres, state.cells[i].temp);
        double bsw = state.cells[i].flui.BSW * ba / (bo + ba * state.cells[i].flui.BSW - state.cells[i].flui.BSW * bo);

        double BSWini = state.cells[i].flui.BSW;
        double BSWF = 0.;
        double viscLini = 0 * 30 + 1 * state.cells[i].flui.VisOM(tL);
        double viscHini = 0 * 20 + 1 * state.cells[i].flui.VisOM(tH);
        double viscLF = viscLini;
        double viscHF = viscHini;

        if (state.cells[i].acsr.tipo == 2) {
            double rsF = state.cells[i].acsr.injl.FluidoPro.RS(state.cells[i].pres, state.cells[i].temp) * (6.29 / 35.31467);
            if (state.input.nfluP > 1 || (*state.globals).chaverede != 0) {
                double bswaux = state.cells[i].acsr.injl.FluidoPro.BSW;
                double contrabsw = 1. - bswaux;
                if (contrabsw < (*state.globals).localtiny)
                    contrabsw = 0.9 * (*state.globals).localtiny;
                double rhoPSTF;
                if (contrabsw > (*state.globals).localtiny)
                    rhoPSTF = (1000 * 141.5 / (131.5 + state.cells[i].acsr.injl.FluidoPro.API)) + (bswaux / contrabsw) * 1000 * state.cells[i].acsr.injl.FluidoPro.Denag + state.cells[i].acsr.injl.FluidoPro.Deng * 1.225 * rsF;
                else
                    rhoPSTF = 1000 * state.cells[i].acsr.injl.FluidoPro.Denag;
                if (contrabsw > (*state.globals).localtiny) {
                    fonteP *= (1. / rhoPSTF);
                    fonteA *= ((bswaux / contrabsw) / rhoPSTF);
                } else {
                    fonteP = 0.;
                    fonteA *= (1 / (state.cells[i].acsr.injl.FluidoPro.BSW * 1000 * state.cells[i].acsr.injl.FluidoPro.Denag));
                }
                BSWF = state.cells[i].acsr.injl.FluidoPro.BSW;
                viscLF = 0 * 30 + 1 * state.cells[i].acsr.injl.FluidoPro.VisOM(tL);
                viscHF = 0 * 20 + 1 * state.cells[i].acsr.injl.FluidoPro.VisOM(tH);
            }
        } else if (state.cells[i].acsr.tipo == 10) {
            double rsF = state.cells[i].acsr.injm.FluidoPro.RS(state.cells[i].pres, state.cells[i].temp) * (6.29 / 35.31467);
            if (state.input.nfluP > 1 || (*state.globals).chaverede != 0) {
                double bswaux = state.cells[i].acsr.injm.FluidoPro.BSW;
                double contrabsw = 1. - bswaux;
                if (contrabsw < (*state.globals).localtiny)
                    contrabsw = 0.9 * (*state.globals).localtiny;
                double rhoPSTF;
                if (contrabsw > (*state.globals).localtiny)
                    rhoPSTF = (1000 * 141.5 / (131.5 + state.cells[i].acsr.injm.FluidoPro.API)) +
                              (bswaux / contrabsw) * 1000 *
                                  state.cells[i].acsr.injm.FluidoPro.Denag +
                              state.cells[i].acsr.injm.FluidoPro.Deng * 1.225 * rsF;
                else
                    rhoPSTF = 1000 * state.cells[i].acsr.injm.FluidoPro.Denag;
                if (contrabsw > (*state.globals).localtiny) {
                    fonteP *= (1. / rhoPSTF);
                    fonteA *= ((bswaux / contrabsw) / rhoPSTF);
                } else {
                    fonteP = 0.;
                    fonteA *= (1 / (state.cells[i].acsr.injm.FluidoPro.BSW * 1000 * state.cells[i].acsr.injm.FluidoPro.Denag));
                }
                BSWF = state.cells[i].acsr.injm.FluidoPro.BSW;
                viscLF = 0 * 30 + 1 * state.cells[i].acsr.injm.FluidoPro.VisOM(tL);
                viscHF = 0 * 20 + 1 * state.cells[i].acsr.injm.FluidoPro.VisOM(tH);
            }
        } else if (state.cells[i].acsr.tipo == 3) {
            double rsF = state.cells[i].acsr.ipr.FluidoPro.RS(state.cells[i].pres, state.cells[i].temp) * (6.29 / 35.31467);
            if (state.input.nfluP > 1 || (*state.globals).chaverede != 0) {
                double bswaux = state.cells[i].acsr.ipr.FluidoPro.BSW;
                double contrabsw = 1. - bswaux;
                if (contrabsw <= (*state.globals).localtiny)
                    contrabsw = 0.9 * (*state.globals).localtiny;
                double rhoPSTF;
                if (contrabsw > 0)
                    rhoPSTF = (1000 * 141.5 / (131.5 + state.cells[i].acsr.ipr.FluidoPro.API)) + (bswaux / contrabsw) * 1000 * state.cells[i].acsr.ipr.FluidoPro.Denag + state.cells[i].acsr.ipr.FluidoPro.Deng * 1.225 * rsF;
                else
                    rhoPSTF = 1000 * state.cells[i].acsr.ipr.FluidoPro.Denag;
                if (contrabsw > (*state.globals).localtiny) {
                    fonteP *= (1. / rhoPSTF);
                    fonteA *= ((bswaux / contrabsw) / rhoPSTF);
                } else {
                    fonteP = 0.;
                    fonteA *= (1 / (state.cells[i].acsr.ipr.FluidoPro.BSW * 1000 * state.cells[i].acsr.ipr.FluidoPro.Denag));
                }
                BSWF = state.cells[i].acsr.ipr.FluidoPro.BSW;
                viscLF = 0 * 30 + 1 * state.cells[i].acsr.ipr.FluidoPro.VisOM(tL);
                viscHF = 0 * 20 + 1 * state.cells[i].acsr.ipr.FluidoPro.VisOM(tH);
            }
        } else if (state.cells[i].acsr.tipo == 9) {

            double rsF = fluF.RS(state.cells[i].pres, state.cells[i].temp) * (6.29 / 35.31467);
            if (state.input.nfluP > 1 || (*state.globals).chaverede != 0) {
                double bswaux = fluF.BSW;
                double contrabsw = 1. - bswaux;
                if (contrabsw <= (*state.globals).localtiny)
                    contrabsw = 0.9 * (*state.globals).localtiny;
                double rhoPSTF;
                if (contrabsw > (*state.globals).localtiny)
                    rhoPSTF = (1000 * 141.5 / (131.5 + fluF.API)) + (bswaux / contrabsw) * 1000 * fluF.Denag + fluF.Deng * 1.225 * rsF;
                else
                    rhoPSTF = 1000 * fluF.Denag;
                if (contrabsw > (*state.globals).localtiny) {
                    fonteP *= (1. / rhoPSTF);
                    fonteA *= ((bswaux / contrabsw) / rhoPSTF);
                } else {
                    fonteP = 0.;
                    fonteA *= (1 / (fluF.BSW * 1000 * fluF.Denag));
                }
                BSWF = fluF.BSW;
                viscLF = 0 * 30 + 1 * fluF.VisOM(tL);
                viscHF = 0 * 20 + 1 * fluF.VisOM(tH);
            }
        } else if (state.cells[i].acsr.tipo == 15) {
            double rsF = state.cells[i].acsr.radialPoro.flup.RS(state.cells[i].pres, state.cells[i].temp) * (6.29 / 35.31467);
            if (state.input.nfluP > 1 || (*state.globals).chaverede != 0) {
                double bswaux = state.cells[i].acsr.radialPoro.BSW;
                double contrabsw = 1. - bswaux;
                if (contrabsw <= (*state.globals).localtiny)
                    contrabsw = 0.9 * (*state.globals).localtiny;
                double rhoPSTF;
                if (contrabsw > 0)
                    rhoPSTF = (1000 * 141.5 / (131.5 + state.cells[i].acsr.radialPoro.flup.API)) + (bswaux / contrabsw) * 1000 * state.cells[i].acsr.radialPoro.flup.Denag + state.cells[i].acsr.radialPoro.flup.Deng * 1.225 * rsF;
                else
                    rhoPSTF = 1000 * state.cells[i].acsr.radialPoro.flup.Denag;
                if (contrabsw > (*state.globals).localtiny) {
                    fonteP *= (1. / rhoPSTF);
                    fonteA *= ((bswaux / contrabsw) / rhoPSTF);
                } else {
                    fonteP = 0.;
                    fonteA *= (1 / (state.cells[i].acsr.radialPoro.BSW * 1000 * state.cells[i].acsr.radialPoro.flup.Denag));
                }
                BSWF = state.cells[i].acsr.radialPoro.BSW;
                viscLF = 0 * 30 + 1 * state.cells[i].acsr.radialPoro.flup.VisOM(tL);
                viscHF = 0 * 20 + 1 * state.cells[i].acsr.radialPoro.flup.VisOM(tH);
            }
        } else if (state.cells[i].acsr.tipo == 16) {
            double rsF = state.cells[i].acsr.poroso2D.dados.flup.RS(state.cells[i].pres, state.cells[i].temp) * (6.29 / 35.31467);
            if (state.input.nfluP > 1 || (*state.globals).chaverede != 0) {
                double bswaux = state.cells[i].acsr.poroso2D.dados.transfer.BSW;
                double contrabsw = 1. - bswaux;
                if (contrabsw <= (*state.globals).localtiny)
                    contrabsw = 0.9 * (*state.globals).localtiny;
                double rhoPSTF;
                if (contrabsw > 0)
                    rhoPSTF = (1000 * 141.5 / (131.5 + state.cells[i].acsr.poroso2D.dados.flup.API)) + (bswaux / contrabsw) * 1000 * state.cells[i].acsr.poroso2D.dados.flup.Denag + state.cells[i].acsr.poroso2D.dados.flup.Deng * 1.225 * rsF;
                else
                    rhoPSTF = 1000 * state.cells[i].acsr.poroso2D.dados.flup.Denag;
                if (contrabsw > (*state.globals).localtiny) {
                    fonteP *= (1. / rhoPSTF);
                    fonteA *= ((bswaux / contrabsw) / rhoPSTF);
                } else {
                    fonteP = 0.;
                    fonteA *= (1 / (state.cells[i].acsr.poroso2D.dados.transfer.BSW * 1000 * state.cells[i].acsr.poroso2D.dados.flup.Denag));
                }
                BSWF = state.cells[i].acsr.poroso2D.dados.transfer.BSW;
                viscLF = 0 * 30 + 1 * state.cells[i].acsr.poroso2D.dados.flup.VisOM(tL);
                viscHF = 0 * 20 + 1 * state.cells[i].acsr.poroso2D.dados.flup.VisOM(tH);
            }
        } else if ((fabs(fonteO) > (*state.globals).localtiny && state.cells[i].acsr.tipo != 2 && state.cells[i].acsr.tipo != 3 &&
                    state.cells[i].acsr.tipo != 9 && state.cells[i].acsr.tipo != 15 && state.cells[i].acsr.tipo != 16) ||
                   (fabs(fonteG) > (*state.globals).localtiny && state.cells[i].acsr.tipo != 1 && state.cells[i].acsr.tipo != 2 && state.cells[i].acsr.tipo != 3 && state.cells[i].acsr.tipo != 9 && state.cells[i].acsr.tipo != 15 && state.cells[i].acsr.tipo != 16)) {
            if (state.cells[i].acsr.tipo == 5 || state.cells[i].acsr.tipo == 8) {
                double rsF = state.cells[i].flui.RS(state.cells[i].pres, state.cells[i].temp) * (6.29 / 35.31467);
                if (state.input.nfluP > 1 || (*state.globals).chaverede != 0) {
                    double rhoPSTF = (1 - state.cells[i].flui.BSW) * (1000 * 141.5 / (131.5 + state.cells[i].flui.API)) + state.cells[i].flui.BSW * 1000 * state.cells[i].flui.Denag;
                    fonteP *= ((1 - state.cells[i].flui.BSW) / rhoPSTF);
                    fonteA *= (state.cells[i].flui.BSW / rhoPSTF);
                    BSWF = state.cells[i].flui.BSW;
                    viscLF = 0 * 30 + 1 * state.cells[i].flui.VisOM(tL);
                    viscHF = 0 * 20 + 1 * state.cells[i].flui.VisOM(tH);
                }
            } else if ((*state.cells[i].acsrL).tipo == 5 || (*state.cells[i].acsrL).tipo == 8) {
                if (i > 0) {
                    double rsF = state.cells[i - 1].flui.RS(state.cells[i].pres, state.cells[i].temp) * (6.29 / 35.31467);
                    if (state.input.nfluP > 1 || (*state.globals).chaverede != 0) {
                        double rhoPSTF = (1 - state.cells[i - 1].flui.BSW) * (1000 * 141.5 / (131.5 + state.cells[i - 1].flui.API)) +
                                         state.cells[i - 1].flui.BSW * 1000 *
                                             state.cells[i - 1].flui.Denag;
                        fonteP *= ((1 - state.cells[i - 1].flui.BSW) / rhoPSTF);
                        fonteA *= (state.cells[i - 1].flui.BSW / rhoPSTF);
                        BSWF = state.cells[i - 1].flui.BSW;
                        viscLF = 0 * 30 + 1 * state.cells[i - 1].flui.VisOM(tL);
                        viscHF = 0 * 20 + 1 * state.cells[i - 1].flui.VisOM(tH);
                    }
                } else {
                    double rsF = state.cells[i].flui.RS(state.cells[i].pres, state.cells[i].temp) * (6.29 / 35.31467);
                    if (state.input.nfluP > 1 || (*state.globals).chaverede != 0) {
                        double rhoPSTF = (1 - state.cells[i].flui.BSW) * (1000 * 141.5 / (131.5 + state.cells[i].flui.API)) +
                                         state.cells[i].flui.BSW * 1000 *
                                             state.cells[i].flui.Denag;
                        fonteP *= ((1 - state.cells[i].flui.BSW) / rhoPSTF);
                        fonteA *= (state.cells[i].flui.BSW / rhoPSTF);
                        BSWF = state.cells[i].flui.BSW;
                        viscLF = 0 * 30 + 1 * state.cells[i].flui.VisOM(tL);
                        viscHF = 0 * 20 + 1 * state.cells[i].flui.VisOM(tH);
                    }
                }
            }
        }
        if (titF > 1. - 1e-15) {
            fonteO = 0.;
            fonteP = 0.;
            fonteA = 0.;
        }

        double volpesFim = hol * (1 - bet) * (1 - bsw) / bo;
        double volaguaFim = hol * (1 - bet) * bsw; // nao deveria ser dividido por Bo???????????
        double MultPe;
        double MultPd = 0.;
        double residuoP;
        double MultAe = 0.;
        double MultAd;
        double residuoA;
        if (state.input.nfluP > 1 || (*state.globals).chaverede != 0) {
            MultPe = 0.;
            if (state.cells[i].flui.dStockTankVaporMassFraction < 1. - 1e-15)
                MultPe = state.cells[i].QL * (1 - betI0) * (1 - bsw0) / bo0;
            MultPd = 0.;
            if (state.cells[i + 1].flui.dStockTankVaporMassFraction < 1. - 1e-15)
                MultPd = state.cells[i + 1].QL * (1 - betI1) * (1 - bsw1) / bo1;
            residuoP = (volpesFim - state.cells[i].VolPesaST) * at / dt + (MultPd - MultPe) / dx - fonteP / dx;
            MultAe = 0.;
            if (state.cells[i].flui.dStockTankVaporMassFraction < 1. - 1e-15)
                MultAe = state.cells[i].QL * (1 - betI0) * bsw0 / bo0;
            MultAd = 0.;
            if (state.cells[i + 1].flui.dStockTankVaporMassFraction < 1. - 1e-15)
                MultAd = state.cells[i + 1].QL * (1 - betI1) * bsw1 / bo1;
            residuoA = (volaguaFim - state.cells[i].VolAguaST) * at / dt + (MultAd - MultAe) / dx - fonteA / dx;
        }
        if ((state.input.nfluP > 1 && (state.input.flashCompleto == 0 || state.cells[i].flui.viscBlackOil == 1)) || (*state.globals).chaverede != 0) {
            if (volpesFim > 1e-3 && (fonteP > 0 || (MultPd < 0 || MultPe > 0))) {
                VISCL[i] = (dt * (viscLF * fonteP / dx + viscLini * residuoP - (viscL1 * MultPd - viscL0 * MultPe) / dx) + viscLini * state.cells[i].VolPesaST * at) / (volpesFim * at);
                VISCH[i] = (dt * (viscHF * fonteP / dx + viscHini * residuoP - (viscH1 * MultPd - viscH0 * MultPe) / dx) + viscHini * state.cells[i].VolPesaST * at) / (volpesFim * at);
            } else {
                VISCL[i] = viscLini;
                VISCH[i] = viscHini;
            }
        }
        if (state.input.nfluP > 1 || (*state.globals).chaverede != 0) {
            if ((volaguaFim + volpesFim) > 1e-3 && ((fonteP > 0 || fonteA > 0) ||
                                                    ((MultPd < 0 || MultPe > 0) || (MultAd < 0 || MultAe > 0))))
                BSW[i] = (dt * (BSWF * (fonteA + fonteP) / dx + BSWini * (residuoA + residuoP) - (BSW1 * (MultAd + MultPd) - BSW0 * (MultAe + MultPe)) / dx) + BSWini * (state.cells[i].VolAguaST + state.cells[i].VolPesaST) * at) / ((volaguaFim + volpesFim) * at);
            else
                BSW[i] = BSWini;
        }
        state.cells[i].VolPesaST = volpesFim;
        state.cells[i].VolAguaST = volaguaFim;
    }
    for (int i = imin; i <= state.lastCell - 1; i++) {
        if (state.input.nfluP > 1 || (*state.globals).chaverede != 0) {
            state.cells[i].flui = fluC[i];
            state.cells[i].flui.BSW = BSW[i];
            if (state.input.flashCompleto == 0 || state.cells[i].flui.viscBlackOil == 1) {
                state.cells[i].flui.LVisL = VISCL[i];
                state.cells[i].flui.LVisH = VISCH[i];
                state.cells[i].flui.TempL = tL;
                state.cells[i].flui.TempH = tH;
            }
        }
        if (state.input.flashCompleto == 0)
            state.cells[i].flui.RenovaFluido();
        state.updaters.correctGasSpecificGravity(i);
    }

    if ((*state.globals).chaverede == 0 || state.endNode == 1 || state.cells[state.lastCell].Mliqini > (*state.globals).localtiny) {
        if (state.input.nfluP > 1 || (*state.globals).chaverede != 0) {
            state.cells[state.lastCell].flui = state.cells[state.lastCell - 1].flui;
            state.cells[state.lastCell].flui.BSW = BSW[state.lastCell - 1];
            if (state.input.flashCompleto == 0) {
                state.cells[state.lastCell].flui.LVisL = VISCL[state.lastCell - 1];
                state.cells[state.lastCell].flui.LVisH = VISCH[state.lastCell - 1];
                state.cells[state.lastCell].flui.TempL = tL;
                state.cells[state.lastCell].flui.TempH = tH;
            }
        }
    } else {
        state.cells[state.lastCell].flui = fluiRev;
        state.cells[state.lastCell].flui.BSW = fluiRev.BSW;
        state.cells[state.lastCell].flui.LVisL = fluiRev.LVisL;
        state.cells[state.lastCell].flui.LVisH = fluiRev.LVisH;
        state.cells[state.lastCell].flui.TempL = fluiRev.TempL;
        state.cells[state.lastCell].flui.TempH = fluiRev.TempH;
    }
    delete[] fluC;
}

void transportPhaseMolarFractions(const CompositionState &state, ProFlu fluiRev) {
    Vcr<double> BSW(state.lastCell);
    Vcr<double> denag(state.lastCell);
    Vcr<double> VISCL(state.lastCell);
    Vcr<double> VISCH(state.lastCell);
    int ncomp = state.input.npseudo;
    ProFlu *fluC;
    fluC = new ProFlu[state.lastCell];
    Vcr<double> fracMol0(ncomp);
    Vcr<double> fracMol1(ncomp);
    Vcr<double> fracMol0O(ncomp);
    Vcr<double> fracMol1O(ncomp);
    Vcr<double> fracMol0G(ncomp);
    Vcr<double> fracMol1G(ncomp);
    Vcr<double> fracMolF(ncomp);
    if (state.cells[0].acsr.tipo == 15 && (state.cells[0].fontemassCR + state.cells[0].fontemassGR + state.cells[0].fontemassLR) > 0.)
        state.cells[0].flui.BSW = state.cells[0].acsr.radialPoro.BSW;
    else if (state.cells[0].acsr.tipo == 16 && (state.cells[0].fontemassCR + state.cells[0].fontemassGR + state.cells[0].fontemassLR) > 0.)
        state.cells[0].flui.BSW = state.cells[0].acsr.poroso2D.dados.transfer.BSW;
    double dt = state.cells[1].dt;
    double tL = 0.;
    double tH = 70.;
    if (state.compositionalRefreshCounter == state.input.miniTabAtraso  && state.input.miniTabAtraso > 0)
        (*state.globals).modoTransiente = 0;
    int imin = 1;
    if ((*state.globals).chaverede != 0 && (state.cells[0].acsr.tipo == 10 || state.input.ConContEntrada > 0))
        imin = 0;

    auto normalizarFracoes = [](vector<double>& fracMolFase, double* fracMolOriginal, int npseudo) {
        for (int kfrac = 0; kfrac < npseudo; kfrac++) {
        	fracMolFase[kfrac] = fracMolOriginal[kfrac];
        }
    	// Encontrar o menor valor
        double menorFracFase = fracMolFase[0];
        for (int kfrac = 1; kfrac < npseudo; kfrac++) {
            if (menorFracFase > fracMolFase[kfrac]) {
                menorFracFase = fracMolFase[kfrac];
            }
        }
        menorFracFase=0.;

        // Subtrair o menor valor
        for (int kfrac = 0; kfrac < npseudo; kfrac++) {
        	fracMolFase[kfrac] -= menorFracFase;
        }

        // Calcular soma total
        double fracTotFase = 0.;
        for (int kfrac = 0; kfrac < npseudo; kfrac++) {
            fracTotFase += fracMolFase[kfrac];
        }

        // Normalizar
        if (fracTotFase > 0) {
        	for (int kfrac = 0; kfrac < npseudo; kfrac++) {
        		fracMolFase[kfrac] /= fracTotFase;
            }
        }
    };

    for (int i = imin; i < state.lastCell; i++) {

        double pesoMol0;
        double pesoMol1;
        double pesoMol0O;
        double pesoMol1O;
        double pesoMol0G;
        double pesoMol1G;
        double pesoMolF;
        double pesoMolC;
        double at = state.cells[i].duto.area;
        double dx = state.cells[i].dx;
        double ti0;
        if (state.cells[i].VTemper < 0.)
            ti0 = state.cells[i].temp;
        else {
            if (i > 0)
                ti0 = state.cells[i - 1].temp;
            else if (state.input.ConContEntrada == 1)
                ti0 = state.inletTemperature;
            else
                ti0 = state.cells[i].temp;
        }
        double ti1 = state.cells[i].temp;
        if (state.cells[i + 1].VTemper < 0.)
            ti1 = state.cells[i + 1].temp;

        double titF = 0.;
        ProFlu fluF;

        double boF = 1.;
        double baF = 1.;
        double fwF = 1.;
        double rhoOF = 900.;
        double rhoWF = 1000.;

        double betI0;
        double bo0;
        double ba0;
        double bsw0;
        double BSW0;
        double denag0;
        double viscL0;
        double viscH0;
        if (i > 0 || state.input.ConContEntrada == 0) {
            if (i > 0 && state.cells[i].QG >= 0.)
                betI0 = state.cells[i - 1].betPigD;
            else
                betI0 = state.cells[i].betPigE;
        } else {
            if (state.cells[i].QG >= 0.)
                betI0 = state.inletCompletionFraction;
            else
                betI0 = state.cells[i].betPigE;
        }
        if ((i > 0 || state.input.ConContEntrada == 1) && state.cells[i].QL >= 0.) {
            double pe;
            double te;
            if (i > 0) {
                pe = state.cells[i - 1].pres;
                te = state.cells[i - 1].temp;
            } else {
                pe = state.inletPressure;
                te = state.inletTemperature;
            }
            if (state.input.ConContEntrada == 0)
                betI0 = state.cells[i - 1].betPigD; // testeBeta
            else
                betI0 = state.inletCompletionFraction; // testeBeta
            double rs0 = (*state.cells[i].fluiL).RS(pe, te);
            bo0 = (*state.cells[i].fluiL).BOFunc(pe, te, rs0);
            ba0 = (*state.cells[i].fluiL).BAFunc(pe, te);
            bsw0 = (*state.cells[i].fluiL).BSW * ba0 / (bo0 + ba0 * (*state.cells[i].fluiL).BSW - (*state.cells[i].fluiL).BSW * bo0);
            BSW0 = (*state.cells[i].fluiL).BSW;
            denag0 = (*state.cells[i].fluiL).Denag;
            viscL0 = 0 * 30 + 1 * (*state.cells[i].fluiL).VisOM(tL);
            viscH0 = 0 * 20 + 1 * (*state.cells[i].fluiL).VisOM(tH);
        } else {
            betI0 = state.cells[i].betPigE; // testebeta
            double rs0 = state.cells[i].flui.RS(state.cells[i].pres, state.cells[i].temp);
            bo0 = state.cells[i].flui.BOFunc(state.cells[i].pres, state.cells[i].temp, rs0);
            ba0 = state.cells[i].flui.BAFunc(state.cells[i].pres, state.cells[i].temp);
            bsw0 = state.cells[i].flui.BSW * ba0 / (bo0 + ba0 * state.cells[i].flui.BSW - state.cells[i].flui.BSW * bo0);
            BSW0 = state.cells[i].flui.BSW;
            denag0 = state.cells[i].flui.Denag;
            viscL0 = 0 * 30 + 1 * state.cells[i].flui.VisOM(tL);
            viscH0 = 0 * 20 + 1 * state.cells[i].flui.VisOM(tH);
        }
        if (bo0 < 1e-15)
            bo0 = 1e-15;

        double betI1 = state.cells[i].betPigD;
        double rs1 = state.cells[i].flui.RS(state.cells[i].pres, state.cells[i].temp);
        double bo1 = state.cells[i].flui.BOFunc(state.cells[i].pres, state.cells[i].temp, rs1);
        double ba1 = state.cells[i].flui.BAFunc(state.cells[i].pres, state.cells[i].temp);
        double bsw1 = state.cells[i].flui.BSW * ba1 / (bo1 + ba1 * state.cells[i].flui.BSW - state.cells[i].flui.BSW * bo1);
        double BSW1 = state.cells[i].flui.BSW;
        double denag1 = state.cells[i].flui.Denag;
        double viscL1 = 0 * 30 + 1 * state.cells[i].flui.VisOM(tL);
        double viscH1 = 0 * 20 + 1 * state.cells[i].flui.VisOM(tH);

        // betI1 = celula[i + 1].betPigE;    //duvidabeta
        if (state.cells[i + 1].QL < 0.) {
            betI1 = state.cells[i + 1].betPigE; // testeBet
            rs1 = state.cells[i + 1].flui.RS(state.cells[i + 1].pres, state.cells[i + 1].temp);
            bo1 = state.cells[i + 1].flui.BOFunc(state.cells[i + 1].pres, state.cells[i + 1].temp, rs1);
            ba1 = state.cells[i + 1].flui.BAFunc(state.cells[i + 1].pres, state.cells[i + 1].temp);
            bsw1 = state.cells[i + 1].flui.BSW * ba1 / (bo1 + ba1 * state.cells[i + 1].flui.BSW - state.cells[i + 1].flui.BSW * bo1);
            BSW1 = state.cells[i + 1].flui.BSW;
            denag1 = state.cells[i + 1].flui.Denag;
            viscL1 = 0 * 30 + 1 * state.cells[i + 1].flui.VisOM(tL);
            viscH1 = 0 * 20 + 1 * state.cells[i + 1].flui.VisOM(tH);
        }
        if (bo1 < 1e-15)
            bo1 = 1e-15;
        double titV0;
        double titV1;
        double vazMasLiq0 = state.cells[i].Mliqini;
        double vazMasLiq1 = state.cells[i + 1].Mliqini;
        double vazMasGas0 = state.cells[i].MC - state.cells[i].Mliqini;
        double vazMasGas1 = state.cells[i + 1].MC - state.cells[i + 1].Mliqini;
        pesoMolC = 0.;
        for (int kfrac = 0; kfrac < ncomp; kfrac++) {
            pesoMolC += state.cells[i].flui.masMol[kfrac] * state.cells[i].flui.fracMol[kfrac];
        }

        /*double menorFracFase = 0.;
        for (int kfrac = 0; kfrac < ncomp; kfrac++) {
            if (menorFracFase > fluC[i].fracMol[kfrac]) {
                menorFracFase = fluC[i].fracMol[kfrac];
            }
        }
        double fracTotFase = 0.;
        for (int kfrac = 0; kfrac < ncomp; kfrac++)
            fluC[i].fracMol[kfrac] -= menorFracFase;
        for (int kfrac = 0; kfrac < ncomp; kfrac++)
            fracTotFase += fluC[i].fracMol[kfrac];
        for (int kfrac = 0; kfrac < ncomp; kfrac++)
            fluC[i].fracMol[kfrac] /= fracTotFase;*/

        if ((i > 0 || state.input.ConContEntrada == 1) && state.cells[i].Mliqini >= 0.) {
            double pe;
            double te;
            if (i > 0) {
                pe = state.cells[i - 1].pres;
                te = state.cells[i - 1].temp;
            } else {
                pe = state.inletPressure;
                te = state.inletTemperature;
            }
            double fwV = bsw0;
            double rhoOV = (*state.cells[i].fluiL).MasEspoleo(pe, te);
            double rhoWV = (*state.cells[i].fluiL).MasEspAgua(pe, te);
            double titLocal=(*state.cells[i].fluiL).FracMass(pe, te);
            titV0 = (1 - fwV) * rhoOV / ((1 - fwV) * rhoOV + fwV * rhoWV);
            vazMasLiq0 -= betI0 * state.cells[i].QL * state.cells[i].fluicol.MasEspFlu(pe, te);
            vazMasLiq0 *= titV0;
            pesoMol0O = 0;
            vector<double> fracMolFase(ncomp);
            normalizarFracoes(fracMolFase, (*state.cells[i].fluiL).oCalculatedLiqComposition, ncomp);
            for (int kfrac = 0; kfrac < ncomp; kfrac++) {
                fracMol0O[kfrac] = fracMolFase[kfrac];
                pesoMol0O += (*state.cells[i].fluiL).masMol[kfrac] * fracMolFase[kfrac];
            }
            if(pesoMol0O<1e-3 || titLocal>1.-1e-3 || state.cells[i].alfL>1.-1e-3){
                pesoMol0O = 0;
                for (int kfrac = 0; kfrac < ncomp; kfrac++) {
                    fracMol0O[kfrac] = (*state.cells[i].fluiL).fracMol[kfrac];
                    pesoMol0O += (*state.cells[i].fluiL).masMol[kfrac] * (*state.cells[i].fluiL).fracMol[kfrac];
                }
            }
        } else {
            double fwV = bsw0;
            double rhoOV = state.cells[i].flui.MasEspoleo(state.cells[i].pres, state.cells[i].temp);
            double rhoWV = state.cells[i].flui.MasEspAgua(state.cells[i].pres, state.cells[i].temp);
            double titLocal=state.cells[i].flui.FracMass(state.cells[i].pres, state.cells[i].temp);
            titV0 = (1 - fwV) * rhoOV / ((1 - fwV) * rhoOV + fwV * rhoWV);
            vazMasLiq0 -= betI0 * state.cells[i].QL * state.cells[i].fluicol.MasEspFlu(state.cells[i].pres, state.cells[i].temp);
            vazMasLiq0 *= titV0;
            pesoMol0O = 0;
            vector<double> fracMolFase(ncomp);
            normalizarFracoes(fracMolFase, state.cells[i].flui.oCalculatedLiqComposition, ncomp);
            for (int kfrac = 0; kfrac < ncomp; kfrac++) {
                fracMol0O[kfrac] = fracMolFase[kfrac];
                pesoMol0O += state.cells[i].flui.masMol[kfrac] * fracMolFase[kfrac];
            }
            if(pesoMol0O<1e-3 || titLocal>1.-1e-3 || state.cells[i].alf>1.-1e-3){
                pesoMol0O = 0;
                for (int kfrac = 0; kfrac < ncomp; kfrac++) {
                    fracMol0O[kfrac] = state.cells[i].flui.fracMol[kfrac];
                    pesoMol0O += state.cells[i].flui.masMol[kfrac] * state.cells[i].flui.fracMol[kfrac];
                }
            }
        }
        if (state.cells[i + 1].Mliqini >= 0.) {
            double fwV = bsw1;
            double rhoOV = state.cells[i].flui.MasEspoleo(state.cells[i].pres, state.cells[i].temp);
            double rhoWV = state.cells[i].flui.MasEspAgua(state.cells[i].pres, state.cells[i].temp);
            double titLocal=state.cells[i].flui.FracMass(state.cells[i].pres, state.cells[i].temp);
            titV1 = (1 - fwV) * rhoOV / ((1 - fwV) * rhoOV + fwV * rhoWV);
            vazMasLiq1 -= betI1 * state.cells[i].QL * state.cells[i].fluicol.MasEspFlu(state.cells[i].pres, state.cells[i].temp);
            vazMasLiq1 *= titV1;
            pesoMol1O = 0;
            vector<double> fracMolFase(ncomp);
            normalizarFracoes(fracMolFase, state.cells[i].flui.oCalculatedLiqComposition, ncomp);
            for (int kfrac = 0; kfrac < ncomp; kfrac++) {
                fracMol1O[kfrac] = fracMolFase[kfrac];
                pesoMol1O += state.cells[i].flui.masMol[kfrac] * fracMolFase[kfrac];
            }
            if(pesoMol1O<1e-3 || titLocal>1.-1e-3 || state.cells[i].alf>1-1e-3){
                pesoMol1O = 0;
                for (int kfrac = 0; kfrac < ncomp; kfrac++) {
                    fracMol1O[kfrac] = state.cells[i].flui.fracMol[kfrac];
                    pesoMol1O += state.cells[i].flui.masMol[kfrac] * state.cells[i].flui.fracMol[kfrac];
                }
            }
        } else {
            double fwV = bsw1;
            double rhoOV = state.cells[i + 1].flui.MasEspoleo(state.cells[i + 1].pres, state.cells[i + 1].temp);
            double rhoWV = state.cells[i + 1].flui.MasEspAgua(state.cells[i + 1].pres, state.cells[i + 1].temp);
            double titLocal=state.cells[i+1].flui.FracMass(state.cells[i+1].pres, state.cells[i+1].temp);
            titV1 = (1 - fwV) * rhoOV / ((1 - fwV) * rhoOV + fwV * rhoWV);
            vazMasLiq1 -= betI1 * state.cells[i + 1].QL * state.cells[i].fluicol.MasEspFlu(state.cells[i + 1].pres, state.cells[i + 1].temp);
            vazMasLiq1 *= titV1;
            pesoMol1O = 0;
            vector<double> fracMolFase(ncomp);
            normalizarFracoes(fracMolFase, state.cells[i + 1].flui.oCalculatedLiqComposition, ncomp);
            for (int kfrac = 0; kfrac < ncomp; kfrac++) {
                fracMol1O[kfrac] = fracMolFase[kfrac];
                pesoMol1O += state.cells[i + 1].flui.masMol[kfrac] * fracMolFase[kfrac];
            }
            if(pesoMol1O<1e-3 || titLocal>1.-1e-3 || state.cells[i+1].alf>1-1e-3){
                pesoMol1O = 0;
                for (int kfrac = 0; kfrac < ncomp; kfrac++) {
                    fracMol1O[kfrac] = state.cells[i + 1].flui.fracMol[kfrac];
                    pesoMol1O += state.cells[i + 1].flui.masMol[kfrac] * state.cells[i + 1].flui.fracMol[kfrac];
                }
            }
        }

        if ((i > 0 || state.input.ConContEntrada == 1) && (state.cells[i].MC - state.cells[i].Mliqini) >= 0.) {
            double pe;
            double te;
            if (i > 0) {
                pe = state.cells[i - 1].pres;
                te = state.cells[i - 1].temp;
            } else {
                pe = state.inletPressure;
                te = state.inletTemperature;
            }
        	pesoMol0G = 0;
            double titLocal=(*state.cells[i].fluiL).FracMass(pe, te);
            vector<double> fracMolFase(ncomp);
            normalizarFracoes(fracMolFase, (*state.cells[i].fluiL).oCalculatedVapComposition, ncomp);
            for (int kfrac = 0; kfrac < ncomp; kfrac++) {
                fracMol0G[kfrac] = fracMolFase[kfrac];
                pesoMol0G += (*state.cells[i].fluiL).masMol[kfrac] * fracMolFase[kfrac];
            }
            if(pesoMol0G<1e-3 || titLocal<1e-3 || state.cells[i].alfL<1e-3){
                pesoMol0G = 0;
                for (int kfrac = 0; kfrac < ncomp; kfrac++) {
                    fracMol0G[kfrac] = (*state.cells[i].fluiL).fracMol[kfrac];
                    pesoMol0G += (*state.cells[i].fluiL).masMol[kfrac] * (*state.cells[i].fluiL).fracMol[kfrac];
                }
            }
        } else {
            pesoMol0G = 0;
            double titLocal=state.cells[i].flui.FracMass(state.cells[i].pres, state.cells[i].temp);
            vector<double> fracMolFase(ncomp);
            normalizarFracoes(fracMolFase, state.cells[i].flui.oCalculatedVapComposition, ncomp);
            for (int kfrac = 0; kfrac < ncomp; kfrac++) {
                fracMol0G[kfrac] = fracMolFase[kfrac];
                pesoMol0G += state.cells[i].flui.masMol[kfrac] * fracMolFase[kfrac];
            }
            if(pesoMol0G<1e-3 || titLocal<1e-3 || state.cells[i].alf<1e-3){
                pesoMol0G = 0;
                for (int kfrac = 0; kfrac < ncomp; kfrac++) {
                    fracMol0G[kfrac] = state.cells[i].flui.fracMol[kfrac];
                    pesoMol0G += state.cells[i].flui.masMol[kfrac] * state.cells[i].flui.fracMol[kfrac];
                }
            }
        }
        if ((state.cells[i + 1].MC - state.cells[i + 1].Mliqini) >= 0.) {
            pesoMol1G = 0;
            double titLocal=state.cells[i].flui.FracMass(state.cells[i].pres, state.cells[i].temp);
            vector<double> fracMolFase(ncomp);
            normalizarFracoes(fracMolFase, state.cells[i].flui.oCalculatedVapComposition, ncomp);
            for (int kfrac = 0; kfrac < ncomp; kfrac++) {
                fracMol1G[kfrac] = fracMolFase[kfrac];
                pesoMol1G += state.cells[i].flui.masMol[kfrac] * fracMolFase[kfrac];
            }
            if(pesoMol1G<1e-3 || titLocal<1e-3 || state.cells[i].alf<1e-3){
                pesoMol1G = 0;
                for (int kfrac = 0; kfrac < ncomp; kfrac++) {
                    fracMol1G[kfrac] = state.cells[i].flui.fracMol[kfrac];
                    pesoMol1G += state.cells[i].flui.masMol[kfrac] * state.cells[i].flui.fracMol[kfrac];
                }
            }
        } else {
            pesoMol1G = 0;
            double titLocal=state.cells[i+1].flui.FracMass(state.cells[i+1].pres, state.cells[i+1].temp);
            vector<double> fracMolFase(ncomp);
            normalizarFracoes(fracMolFase, state.cells[i + 1].flui.oCalculatedVapComposition, ncomp);
            for (int kfrac = 0; kfrac < ncomp; kfrac++) {
                fracMol1G[kfrac] = fracMolFase[kfrac];
                pesoMol1G += state.cells[i + 1].flui.masMol[kfrac] * fracMolFase[kfrac];
            }
            if(pesoMol1G<1e-3 || titLocal<1e-3 || state.cells[i+1].alf<1e-3){
                pesoMol1G = 0;
                for (int kfrac = 0; kfrac < ncomp; kfrac++) {
                    fracMol1G[kfrac] = state.cells[i + 1].flui.fracMol[kfrac];
                    pesoMol1G += state.cells[i + 1].flui.masMol[kfrac] * state.cells[i + 1].flui.fracMol[kfrac];
                }
            }
        }

        if ((i > 0 || state.input.ConContEntrada == 1) && (vazMasLiq0 + vazMasGas0) >= 0) {
            pesoMol0 = 0;
            for (int kfrac = 0; kfrac < ncomp; kfrac++) {
                fracMol0[kfrac] = (*state.cells[i].fluiL).fracMol[kfrac];
                pesoMol0 += (*state.cells[i].fluiL).masMol[kfrac] * (*state.cells[i].fluiL).fracMol[kfrac];
            }
        } else {
            pesoMol0 = 0;
            for (int kfrac = 0; kfrac < ncomp; kfrac++) {
                fracMol0[kfrac] = state.cells[i].flui.fracMol[kfrac];
                pesoMol0 += state.cells[i].flui.masMol[kfrac] * state.cells[i].flui.fracMol[kfrac];
            }
        }
        if ((vazMasLiq1 + vazMasGas1) >= 0) {
            pesoMol1 = 0;
            for (int kfrac = 0; kfrac < ncomp; kfrac++) {
                fracMol1[kfrac] = state.cells[i].flui.fracMol[kfrac];
                pesoMol1 += state.cells[i].flui.masMol[kfrac] * state.cells[i].flui.fracMol[kfrac];
            }
        } else {
            pesoMol1 = 0;
            for (int kfrac = 0; kfrac < ncomp; kfrac++) {
                fracMol1[kfrac] = state.cells[i + 1].flui.fracMol[kfrac];
                pesoMol1 += state.cells[i + 1].flui.masMol[kfrac] * state.cells[i + 1].flui.fracMol[kfrac];
            }
        }

        double fonteO = state.cells[i].fontemassLR;
        double fonteG = state.cells[i].fontemassGR;
        double fonteP = state.cells[i].fontemassLR;
        double fonteA = state.cells[i].fontemassLR;

        if (state.cells[i].acsr.tipo == 1) {
            if (state.cells[i].acsr.injg.FluidoPro.dCalculatedBeta < 0. || state.cells[i].acsr.injg.FluidoPro.dCalculatedBeta > 1.)
                state.cells[i].acsr.injg.FluidoPro.atualizaPropComp(state.cells[i].pres, state.cells[i].temp, -1, NULL, NULL, state.input.pocinjec);
            else
                state.cells[i].acsr.injg.FluidoPro.atualizaPropComp(state.cells[i].pres, state.cells[i].temp,
                                                               state.cells[i].acsr.injg.FluidoPro.dCalculatedBeta, state.cells[i].acsr.injg.FluidoPro.oCalculatedLiqComposition,
                                                               state.cells[i].acsr.injg.FluidoPro.oCalculatedVapComposition, state.input.pocinjec);
            if (state.cells[i].acsr.injg.QGas > 0.)
                fluF = state.cells[i].acsr.injg.FluidoPro;
            else
                fluF = state.cells[i].flui;
            fwF = 0.;
            titF = 1.;
        } else if (state.cells[i].acsr.tipo == 2) {
            if (state.cells[i].acsr.injl.FluidoPro.dCalculatedBeta < 0. || state.cells[i].acsr.injl.FluidoPro.dCalculatedBeta > 1.)
                state.cells[i].acsr.injl.FluidoPro.atualizaPropComp(state.cells[i].pres, state.cells[i].temp, -1, NULL, NULL, state.input.pocinjec);
            else
                state.cells[i].acsr.injl.FluidoPro.atualizaPropComp(state.cells[i].pres, state.cells[i].temp,
                                                               state.cells[i].acsr.injl.FluidoPro.dCalculatedBeta, state.cells[i].acsr.injl.FluidoPro.oCalculatedLiqComposition,
                                                               state.cells[i].acsr.injl.FluidoPro.oCalculatedVapComposition, state.input.pocinjec);
            if (state.cells[i].acsr.injl.QLiq > 0.)
                fluF = state.cells[i].acsr.injl.FluidoPro;
            else
                fluF = state.cells[i].flui;
            boF = fluF.BOFunc(state.cells[i].pres, state.cells[i].temp);
            baF = fluF.BAFunc(state.cells[i].pres, state.cells[i].temp);
            fwF = fluF.BSW * baF / (boF + baF * fluF.BSW - fluF.BSW * boF);
            rhoOF = fluF.MasEspoleo(state.cells[i].pres, state.cells[i].temp);
            rhoWF = fluF.MasEspAgua(state.cells[i].pres, state.cells[i].temp);
            titF = (1 - fwF) * rhoOF / ((1 - fwF) * rhoOF + fwF * rhoWF);
        } else if (state.cells[i].acsr.tipo == 3) {
            if (state.cells[i].acsr.ipr.FluidoPro.dCalculatedBeta < 0. || state.cells[i].acsr.ipr.FluidoPro.dCalculatedBeta > 1.)
                state.cells[i].acsr.ipr.FluidoPro.atualizaPropComp(state.cells[i].pres, state.cells[i].temp, -1, NULL, NULL, state.input.pocinjec);
            else
                state.cells[i].acsr.ipr.FluidoPro.atualizaPropComp(state.cells[i].pres, state.cells[i].temp,
                                                              state.cells[i].acsr.ipr.FluidoPro.dCalculatedBeta, state.cells[i].acsr.ipr.FluidoPro.oCalculatedLiqComposition,
                                                              state.cells[i].acsr.ipr.FluidoPro.oCalculatedVapComposition, state.input.pocinjec);
            if ((state.cells[i].acsr.ipr.Pres) > state.cells[i].pres)
                fluF = state.cells[i].acsr.ipr.FluidoPro;
            else
                fluF = state.cells[i].flui;
            boF = fluF.BOFunc(state.cells[i].pres, state.cells[i].temp);
            baF = fluF.BAFunc(state.cells[i].pres, state.cells[i].temp);
            fwF = fluF.BSW * baF / (boF + baF * fluF.BSW - fluF.BSW * boF);
            rhoOF = fluF.MasEspoleo(state.cells[i].pres, state.cells[i].temp);
            rhoWF = fluF.MasEspAgua(state.cells[i].pres, state.cells[i].temp);
            titF = (1 - fwF) * rhoOF / ((1 - fwF) * rhoOF + fwF * rhoWF);
        } else if (state.cells[i].acsr.tipo == 10) {
            if (state.cells[i].acsr.injm.FluidoPro.dCalculatedBeta < 0. || state.cells[i].acsr.injm.FluidoPro.dCalculatedBeta > 1.)
                state.cells[i].acsr.injm.FluidoPro.atualizaPropComp(state.cells[i].pres, state.cells[i].temp, -1, NULL, NULL, state.input.pocinjec);
            else
                state.cells[i].acsr.injm.FluidoPro.atualizaPropComp(state.cells[i].pres, state.cells[i].temp,
                                                               state.cells[i].acsr.injm.FluidoPro.dCalculatedBeta, state.cells[i].acsr.injm.FluidoPro.oCalculatedLiqComposition,
                                                               state.cells[i].acsr.injm.FluidoPro.oCalculatedVapComposition, state.input.pocinjec);
            if ((state.cells[i].acsr.injm.MassC + state.cells[i].acsr.injm.MassG + state.cells[i].acsr.injm.MassP) > 0.)
                fluF = state.cells[i].acsr.injm.FluidoPro;
            else
                fluF = state.cells[i].flui;
            boF = fluF.BOFunc(state.cells[i].pres, state.cells[i].temp);
            baF = fluF.BAFunc(state.cells[i].pres, state.cells[i].temp);
            fwF = fluF.BSW * baF / (boF + baF * fluF.BSW - fluF.BSW * boF);
            rhoOF = fluF.MasEspoleo(state.cells[i].pres, state.cells[i].temp);
            rhoWF = fluF.MasEspAgua(state.cells[i].pres, state.cells[i].temp);
            titF = (1 - fwF) * rhoOF / ((1 - fwF) * rhoOF + fwF * rhoWF);
        } else if (state.cells[i].acsr.tipo == 9 && state.cells[i].acsr.fontechk.abertura > 1e-6) {
            if (state.cells[i].acsr.fontechk.fluidoP.dCalculatedBeta < 0. || state.cells[i].acsr.fontechk.fluidoP.dCalculatedBeta > 1.)
                state.cells[i].acsr.fontechk.fluidoP.atualizaPropComp(state.cells[i].pres, state.cells[i].temp, -1, NULL, NULL, state.input.pocinjec);
            else
                state.cells[i].acsr.fontechk.fluidoP.atualizaPropComp(state.cells[i].pres, state.cells[i].temp,
                                                                 state.cells[i].acsr.fontechk.fluidoP.dCalculatedBeta, state.cells[i].acsr.fontechk.fluidoP.oCalculatedLiqComposition,
                                                                 state.cells[i].acsr.fontechk.fluidoP.oCalculatedVapComposition, state.input.pocinjec);
            if (state.cells[i].acsr.fontechk.fluidoPamb.dCalculatedBeta < 0. || state.cells[i].acsr.fontechk.fluidoPamb.dCalculatedBeta > 1.)
                state.cells[i].acsr.fontechk.fluidoPamb.atualizaPropComp(state.cells[i].pres, state.cells[i].temp, -1, NULL, NULL, state.input.pocinjec);
            else
                state.cells[i].acsr.fontechk.fluidoPamb.atualizaPropComp(state.cells[i].pres, state.cells[i].temp,
                                                                    state.cells[i].acsr.fontechk.fluidoPamb.dCalculatedBeta, state.cells[i].acsr.fontechk.fluidoPamb.oCalculatedLiqComposition,
                                                                    state.cells[i].acsr.fontechk.fluidoPamb.oCalculatedVapComposition, state.input.pocinjec);
            if (state.cells[i].acsr.fontechk.presT > state.cells[i].acsr.fontechk.pamb) {
                fluF = state.cells[i].acsr.fontechk.fluidoP;
            } else {
                fluF = state.cells[i].acsr.fontechk.fluidoPamb;
            }
            boF = fluF.BOFunc(state.cells[i].pres, state.cells[i].temp);
            baF = fluF.BAFunc(state.cells[i].pres, state.cells[i].temp);
            fwF = fluF.BSW * baF / (boF + baF * fluF.BSW - fluF.BSW * boF);
            rhoOF = fluF.MasEspoleo(state.cells[i].pres, state.cells[i].temp);
            rhoWF = fluF.MasEspAgua(state.cells[i].pres, state.cells[i].temp);
            titF = (1 - fwF) * rhoOF / ((1 - fwF) * rhoOF + fwF * rhoWF);
        } else if (state.cells[i].acsr.tipo == 15) {
            double tRes = state.cells[i].acsr.radialPoro.tRes;
            if (state.cells[i].acsr.radialPoro.flup.dCalculatedBeta < 0. || state.cells[i].acsr.radialPoro.flup.dCalculatedBeta > 1.)
                state.cells[i].acsr.radialPoro.flup.atualizaPropComp(state.cells[i].pres, tRes, -1, NULL, NULL, state.input.pocinjec);
            else
                state.cells[i].acsr.radialPoro.flup.atualizaPropComp(state.cells[i].pres, tRes,
                                                                state.cells[i].acsr.radialPoro.flup.dCalculatedBeta, state.cells[i].acsr.radialPoro.flup.oCalculatedLiqComposition,
                                                                state.cells[i].acsr.radialPoro.flup.oCalculatedVapComposition, state.input.pocinjec);
            if ((state.cells[i].fontemassLR + state.cells[i].fontemassGR) > 1e-15)
                fluF = state.cells[i].acsr.radialPoro.flup;
            else
                fluF = state.cells[i].flui;
            boF = fluF.BOFunc(state.cells[i].pres, state.cells[i].temp);
            baF = fluF.BAFunc(state.cells[i].pres, state.cells[i].temp);
            fwF = fluF.BSW * baF / (boF + baF * fluF.BSW - fluF.BSW * boF);
            rhoOF = fluF.MasEspoleo(state.cells[i].pres, state.cells[i].temp);
            rhoWF = fluF.MasEspAgua(state.cells[i].pres, state.cells[i].temp);
            titF = (1 - fwF) * rhoOF / ((1 - fwF) * rhoOF + fwF * rhoWF);
        } else if (state.cells[i].acsr.tipo == 16) {
            double tRes = state.cells[i].acsr.poroso2D.dados.tRes;
            if (state.cells[i].acsr.poroso2D.dados.flup.dCalculatedBeta < 0. || state.cells[i].acsr.poroso2D.dados.flup.dCalculatedBeta > 1.)
                state.cells[i].acsr.poroso2D.dados.flup.atualizaPropComp(state.cells[i].pres, tRes, -1, NULL, NULL, state.input.pocinjec);
            else
                state.cells[i].acsr.poroso2D.dados.flup.atualizaPropComp(state.cells[i].pres, tRes,
                                                                    state.cells[i].acsr.poroso2D.dados.flup.dCalculatedBeta, state.cells[i].acsr.poroso2D.dados.flup.oCalculatedLiqComposition,
                                                                    state.cells[i].acsr.poroso2D.dados.flup.oCalculatedVapComposition, state.input.pocinjec);
            if ((state.cells[i].fontemassLR + state.cells[i].fontemassGR) > 1e-15)
                fluF = state.cells[i].acsr.poroso2D.dados.flup;
            else
                fluF = state.cells[i].flui;
            boF = fluF.BOFunc(state.cells[i].pres, state.cells[i].temp);
            baF = fluF.BAFunc(state.cells[i].pres, state.cells[i].temp);
            fwF = fluF.BSW * baF / (boF + baF * fluF.BSW - fluF.BSW * boF);
            rhoOF = fluF.MasEspoleo(state.cells[i].pres, state.cells[i].temp);
            rhoWF = fluF.MasEspAgua(state.cells[i].pres, state.cells[i].temp);
            titF = (1 - fwF) * rhoOF / ((1 - fwF) * rhoOF + fwF * rhoWF);
        } else {
            fluF = state.cells[i].flui;
            titF = 0.;
        }

        pesoMolF = 0;
        for (int j = 0; j < ncomp; j++) {
            pesoMolF += fluF.masMol[j] * fluF.fracMol[j];
        }
        fonteO *= titF;
        double tempMol;
        fluC[i] = state.cells[i].flui;
        double betIV;
        double rsV;
        double boV;
        double baV;
        double bswV;
        double rhoOVol;
        double rhoWVol;
        double titVol;
        int limCorrige = 1;
        if (state.compositionalRefreshCounter == state.input.miniTabAtraso  || state.input.miniTabAtraso == 0)
            limCorrige = 2;
        for (int corrige = 0; corrige < limCorrige; corrige++) {

            betIV = state.cells[i].bet;
            rsV = fluC[i].RS(state.cells[i].pres, state.cells[i].temp);
            boV = fluC[i].BOFunc(state.cells[i].pres, state.cells[i].temp, rsV);
            baV = fluC[i].BAFunc(state.cells[i].pres, state.cells[i].temp);
            bswV = fluC[i].BSW * baV / (boV + baV * state.cells[i].flui.BSW - state.cells[i].flui.BSW * boV);
            rhoOVol = fluC[i].MasEspoleo(state.cells[i].pres, state.cells[i].temp);
            rhoWVol = fluC[i].MasEspAgua(state.cells[i].pres, state.cells[i].temp);
            titVol = (1 - bswV) * rhoOVol / ((1 - bswV) * rhoOVol + bswV * rhoWVol);

            tempMol = (fluC[i].MasEspLiq(state.cells[i].pres, state.cells[i].temp) * (1. - state.cells[i].alf) * (1. - betIV) * titVol +
                       fluC[i].MasEspGas(state.cells[i].pres, state.cells[i].temp) * state.cells[i].alf) *
                      state.cells[i].duto.area * state.cells[i].dx / pesoMolC;

            double menorFrac = 0.;
            for (int kfrac = 0; kfrac < ncomp; kfrac++) {

                /*double vazMol1O;
                double vazMol0O;
                double vazMol1G;
                double vazMol0G;
                vazMol1O=vazMasLiq1*fracMol1O[kfrac]/pesoMol1O;
                vazMol0O=vazMasLiq0*fracMol0O[kfrac]/pesoMol0O;
                vazMol1G=vazMasGas1*fracMol1G[kfrac]/pesoMol1G;
                vazMol0G=vazMasGas0*fracMol0G[kfrac]/pesoMol0G;
            	fluC[i].fracMol[kfrac] = (celula[i].nMolIni * celula[i].flui.fracMol[kfrac] +
                                          ((fonteO + fonteG) * fluF.fracMol[kfrac] / pesoMolF -
                                           (vazMol1O + vazMol1G) +
                                           (vazMol0O + vazMol0G)) *
                                              dt) /
                                         tempMol;*/


                fluC[i].fracMol[kfrac] = (state.cells[i].nMolIni * state.cells[i].flui.fracMol[kfrac] +
                                          ((fonteO + fonteG) * fluF.fracMol[kfrac] / pesoMolF -
                                           (vazMasLiq1 + vazMasGas1) * fracMol1[kfrac] / pesoMol1 +
                                           (vazMasLiq0 + vazMasGas0) * fracMol0[kfrac] / pesoMol0) *
                                              dt) /
                                         tempMol;
                if (menorFrac > fluC[i].fracMol[kfrac]) {
                    menorFrac = fluC[i].fracMol[kfrac];
                }
            }
            double fracTot = 0.;
            if(menorFrac<0.){
            	for (int kfrac = 0; kfrac < ncomp; kfrac++)
            		fluC[i].fracMol[kfrac] -= menorFrac;
            }
            for (int kfrac = 0; kfrac < ncomp; kfrac++)
                fracTot += fluC[i].fracMol[kfrac];
            for (int kfrac = 0; kfrac < ncomp; kfrac++)
                fluC[i].fracMol[kfrac] /= fracTot;
            pesoMolC = 0;
            for (int kfrac = 0; kfrac < ncomp; kfrac++)
                pesoMolC += fluC[i].masMol[kfrac] * state.cells[i].flui.fracMol[kfrac];
            if (state.compositionalRefreshCounter == state.input.miniTabAtraso || state.input.miniTabAtraso == 0) {
                fluC[i].atualizaPropCompStandard();
                if (fluC[i].dCalculatedBeta < 0. || fluC[i].dCalculatedBeta > 1.)
                    fluC[i].atualizaPropComp(state.cells[i].pres, state.cells[i].temp, -1, NULL, NULL, state.input.pocinjec);
                else
                    fluC[i].atualizaPropComp(state.cells[i].pres, state.cells[i].temp, fluC[i].dCalculatedBeta,
                                             fluC[i].oCalculatedLiqComposition,
                                             fluC[i].oCalculatedVapComposition, state.input.pocinjec);
                if (fluC[i].iIER != 0) {
                    int para;
                    para = 0;
                }
            }
        }
        betIV = state.cells[i].bet;
        rsV = state.cells[i].flui.RS(state.cells[i].pres, state.cells[i].temp);
        boV = state.cells[i].flui.BOFunc(state.cells[i].pres, state.cells[i].temp, rsV);
        baV = state.cells[i].flui.BAFunc(state.cells[i].pres, state.cells[i].temp);
        bswV = state.cells[i].flui.BSW * baV / (boV + baV * state.cells[i].flui.BSW - state.cells[i].flui.BSW * boV);
        rhoOVol = state.cells[i].flui.MasEspoleo(state.cells[i].pres, state.cells[i].temp);
        rhoWVol = state.cells[i].flui.MasEspAgua(state.cells[i].pres, state.cells[i].temp);
        titVol = (1 - bswV) * rhoOVol / ((1 - bswV) * rhoOVol + bswV * rhoWVol);
        state.cells[i].nMol = (fluC[i].MasEspLiq(state.cells[i].pres, state.cells[i].temp) * (1. - state.cells[i].alf) * (1. - betIV) * titVol +
                          fluC[i].MasEspGas(state.cells[i].pres, state.cells[i].temp) *
                              state.cells[i].alf) *
                         state.cells[i].duto.area * state.cells[i].dx / pesoMolC;
        fluC[i].Pmol = pesoMolC;
        double hol = 1. - state.cells[i].alf;
        double bet = state.cells[i].bet;

        double rs = state.cells[i].flui.RS(state.cells[i].pres, state.cells[i].temp);
        double bo = state.cells[i].flui.BOFunc(state.cells[i].pres, state.cells[i].temp, rs);
        double ba = state.cells[i].flui.BAFunc(state.cells[i].pres, state.cells[i].temp);
        double bsw = state.cells[i].flui.BSW * ba / (bo + ba * state.cells[i].flui.BSW - state.cells[i].flui.BSW * bo);

        double BSWini = state.cells[i].flui.BSW;
        double BSWF = state.cells[i].flui.BSW;
        double denagini = state.cells[i].flui.Denag;
        double denagF = state.cells[i].flui.Denag;
        double viscLini = 0 * 30 + 1 * state.cells[i].flui.VisOM(tL);
        double viscHini = 0 * 20 + 1 * state.cells[i].flui.VisOM(tH);
        double viscLF = viscLini;
        double viscHF = viscHini;

        if (state.cells[i].acsr.tipo == 2) {
            double rsF = state.cells[i].acsr.injl.FluidoPro.RS(state.cells[i].pres, state.cells[i].temp) * (6.29 / 35.31467);
            if (state.input.nfluP > 0 || (*state.globals).chaverede != 0) {
                double bswaux = state.cells[i].acsr.injl.FluidoPro.BSW;
                double contrabsw = 1. - bswaux;
                if (contrabsw < (*state.globals).localtiny)
                    contrabsw = 0.9 * (*state.globals).localtiny;
                double rhoPSTF;
                if (contrabsw > (*state.globals).localtiny)
                    rhoPSTF = (1000 * 141.5 / (131.5 + state.cells[i].acsr.injl.FluidoPro.API)) + (bswaux / contrabsw) * 1000 * state.cells[i].acsr.injl.FluidoPro.Denag + state.cells[i].acsr.injl.FluidoPro.Deng * 1.225 * rsF;
                else
                    rhoPSTF = 1000 * state.cells[i].acsr.injl.FluidoPro.Denag;
                if (contrabsw > (*state.globals).localtiny) {
                    fonteP *= (1. / rhoPSTF);
                    fonteA *= ((bswaux / contrabsw) / rhoPSTF);
                } else {
                    fonteP = 0.;
                    fonteA *= (1 / (state.cells[i].acsr.injl.FluidoPro.BSW * 1000 * state.cells[i].acsr.injl.FluidoPro.Denag));
                }
                BSWF = state.cells[i].acsr.injl.FluidoPro.BSW;
                denagF = state.cells[i].acsr.injl.FluidoPro.Denag;
                viscLF = 0 * 30 + 1 * state.cells[i].acsr.injl.FluidoPro.VisOM(tL);
                viscHF = 0 * 20 + 1 * state.cells[i].acsr.injl.FluidoPro.VisOM(tH);
            }
        } else if (state.cells[i].acsr.tipo == 10) {
            double rsF = state.cells[i].acsr.injm.FluidoPro.RS(state.cells[i].pres, state.cells[i].temp) * (6.29 / 35.31467);
            if (state.input.nfluP > 0 || (*state.globals).chaverede != 0) {
                double bswaux = state.cells[i].acsr.injm.FluidoPro.BSW;
                double contrabsw = 1. - bswaux;
                if (contrabsw < (*state.globals).localtiny)
                    contrabsw = 0.9 * (*state.globals).localtiny;
                double rhoPSTF;
                if (contrabsw > (*state.globals).localtiny)
                    rhoPSTF = (1000 * 141.5 / (131.5 + state.cells[i].acsr.injm.FluidoPro.API)) +
                              (bswaux / contrabsw) * 1000 *
                                  state.cells[i].acsr.injm.FluidoPro.Denag +
                              state.cells[i].acsr.injm.FluidoPro.Deng * 1.225 * rsF;
                else
                    rhoPSTF = 1000 * state.cells[i].acsr.injm.FluidoPro.Denag;
                if (contrabsw > (*state.globals).localtiny) {
                    fonteP *= (1. / rhoPSTF);
                    fonteA *= ((bswaux / contrabsw) / rhoPSTF);
                } else {
                    fonteP = 0.;
                    fonteA *= (1 / (state.cells[i].acsr.injm.FluidoPro.BSW * 1000 * state.cells[i].acsr.injm.FluidoPro.Denag));
                }
                BSWF = state.cells[i].acsr.injm.FluidoPro.BSW;
                denagF = state.cells[i].acsr.injl.FluidoPro.Denag;
                viscLF = 0 * 30 + 1 * state.cells[i].acsr.injm.FluidoPro.VisOM(tL);
                viscHF = 0 * 20 + 1 * state.cells[i].acsr.injm.FluidoPro.VisOM(tH);
            }
        } else if (state.cells[i].acsr.tipo == 3) {
            double rsF = state.cells[i].acsr.ipr.FluidoPro.RS(state.cells[i].pres, state.cells[i].temp) * (6.29 / 35.31467);
            if (state.input.nfluP > 0 || (*state.globals).chaverede != 0) {
                double bswaux = state.cells[i].acsr.ipr.FluidoPro.BSW;
                double contrabsw = 1. - bswaux;
                if (contrabsw <= (*state.globals).localtiny)
                    contrabsw = 0.9 * (*state.globals).localtiny;
                double rhoPSTF;
                if (contrabsw > 0)
                    rhoPSTF = (1000 * 141.5 / (131.5 + state.cells[i].acsr.ipr.FluidoPro.API)) + (bswaux / contrabsw) * 1000 * state.cells[i].acsr.ipr.FluidoPro.Denag + state.cells[i].acsr.ipr.FluidoPro.Deng * 1.225 * rsF;
                else
                    rhoPSTF = 1000 * state.cells[i].acsr.ipr.FluidoPro.Denag;
                if (contrabsw > (*state.globals).localtiny) {
                    fonteP *= (1. / rhoPSTF);
                    fonteA *= ((bswaux / contrabsw) / rhoPSTF);
                } else {
                    fonteP = 0.;
                    fonteA *= (1 / (state.cells[i].acsr.ipr.FluidoPro.BSW * 1000 * state.cells[i].acsr.ipr.FluidoPro.Denag));
                }
                BSWF = state.cells[i].acsr.ipr.FluidoPro.BSW;
                denagF = state.cells[i].acsr.ipr.FluidoPro.Denag;
                viscLF = 0 * 30 + 1 * state.cells[i].acsr.ipr.FluidoPro.VisOM(tL);
                viscHF = 0 * 20 + 1 * state.cells[i].acsr.ipr.FluidoPro.VisOM(tH);
            }
        } else if (state.cells[i].acsr.tipo == 9) {

            double rsF = fluF.RS(state.cells[i].pres, state.cells[i].temp) * (6.29 / 35.31467);
            if (state.input.nfluP > 0 || (*state.globals).chaverede != 0) {
                double bswaux = fluF.BSW;
                double contrabsw = 1. - bswaux;
                if (contrabsw <= (*state.globals).localtiny)
                    contrabsw = 0.9 * (*state.globals).localtiny;
                double rhoPSTF;
                if (contrabsw > (*state.globals).localtiny)
                    rhoPSTF = (1000 * 141.5 / (131.5 + fluF.API)) + (bswaux / contrabsw) * 1000 * fluF.Denag + fluF.Deng * 1.225 * rsF;
                else
                    rhoPSTF = 1000 * fluF.Denag;
                if (contrabsw > (*state.globals).localtiny) {
                    fonteP *= (1. / rhoPSTF);
                    fonteA *= ((bswaux / contrabsw) / rhoPSTF);
                } else {
                    fonteP = 0.;
                    fonteA *= (1 / (fluF.BSW * 1000 * fluF.Denag));
                }
                BSWF = fluF.BSW;
                denagF = fluF.Denag;
                viscLF = 0 * 30 + 1 * fluF.VisOM(tL);
                viscHF = 0 * 20 + 1 * fluF.VisOM(tH);
            }
        } else if (state.cells[i].acsr.tipo == 15) {
            double rsF = state.cells[i].acsr.radialPoro.flup.RS(state.cells[i].pres, state.cells[i].temp) * (6.29 / 35.31467);
            if (state.input.nfluP > 0 || (*state.globals).chaverede != 0) {
                double bswaux = state.cells[i].acsr.radialPoro.BSW;
                double contrabsw = 1. - bswaux;
                if (contrabsw <= (*state.globals).localtiny)
                    contrabsw = 0.9 * (*state.globals).localtiny;
                double rhoPSTF;
                if (contrabsw > 0)
                    rhoPSTF = (1000 * 141.5 / (131.5 + state.cells[i].acsr.radialPoro.flup.API)) + (bswaux / contrabsw) * 1000 * state.cells[i].acsr.radialPoro.flup.Denag + state.cells[i].acsr.radialPoro.flup.Deng * 1.225 * rsF;
                else
                    rhoPSTF = 1000 * state.cells[i].acsr.radialPoro.flup.Denag;
                if (contrabsw > (*state.globals).localtiny) {
                    fonteP *= (1. / rhoPSTF);
                    fonteA *= ((bswaux / contrabsw) / rhoPSTF);
                } else {
                    fonteP = 0.;
                    fonteA *= (1 / (state.cells[i].acsr.radialPoro.BSW * 1000 * state.cells[i].acsr.radialPoro.flup.Denag));
                }
                BSWF = state.cells[i].acsr.radialPoro.BSW;
                denagF = state.cells[i].acsr.radialPoro.flup.Denag;
                viscLF = 0 * 30 + 1 * state.cells[i].acsr.radialPoro.flup.VisOM(tL);
                viscHF = 0 * 20 + 1 * state.cells[i].acsr.radialPoro.flup.VisOM(tH);
            }
        } else if (state.cells[i].acsr.tipo == 16) {
            double rsF = state.cells[i].acsr.poroso2D.dados.flup.RS(state.cells[i].pres, state.cells[i].temp) * (6.29 / 35.31467);
            if (state.input.nfluP > 0 || (*state.globals).chaverede != 0) {
                double bswaux = state.cells[i].acsr.poroso2D.dados.transfer.BSW;
                double contrabsw = 1. - bswaux;
                if (contrabsw <= (*state.globals).localtiny)
                    contrabsw = 0.9 * (*state.globals).localtiny;
                double rhoPSTF;
                if (contrabsw > 0)
                    rhoPSTF = (1000 * 141.5 / (131.5 + state.cells[i].acsr.poroso2D.dados.flup.API)) + (bswaux / contrabsw) * 1000 * state.cells[i].acsr.poroso2D.dados.flup.Denag + state.cells[i].acsr.poroso2D.dados.flup.Deng * 1.225 * rsF;
                else
                    rhoPSTF = 1000 * state.cells[i].acsr.poroso2D.dados.flup.Denag;
                if (contrabsw > (*state.globals).localtiny) {
                    fonteP *= (1. / rhoPSTF);
                    fonteA *= ((bswaux / contrabsw) / rhoPSTF);
                } else {
                    fonteP = 0.;
                    fonteA *= (1 / (state.cells[i].acsr.poroso2D.dados.transfer.BSW * 1000 * state.cells[i].acsr.poroso2D.dados.flup.Denag));
                }
                BSWF = state.cells[i].acsr.poroso2D.dados.transfer.BSW;
                denagF = state.cells[i].acsr.poroso2D.dados.flup.Denag;
                viscLF = 0 * 30 + 1 * state.cells[i].acsr.poroso2D.dados.flup.VisOM(tL);
                viscHF = 0 * 20 + 1 * state.cells[i].acsr.poroso2D.dados.flup.VisOM(tH);
            }
        } else if ((fabs(fonteO) > (*state.globals).localtiny && state.cells[i].acsr.tipo != 2 && state.cells[i].acsr.tipo != 3 &&
                    state.cells[i].acsr.tipo != 9 && state.cells[i].acsr.tipo != 15 && state.cells[i].acsr.tipo != 16) ||
                   (fabs(fonteG) > (*state.globals).localtiny && state.cells[i].acsr.tipo != 1 && state.cells[i].acsr.tipo != 2 && state.cells[i].acsr.tipo != 3 && state.cells[i].acsr.tipo != 9 && state.cells[i].acsr.tipo != 15 && state.cells[i].acsr.tipo != 16)) {
            if (state.cells[i].acsr.tipo == 5 || state.cells[i].acsr.tipo == 8) {
                double rsF = state.cells[i].flui.RS(state.cells[i].pres, state.cells[i].temp) * (6.29 / 35.31467);
                if (state.input.nfluP > 0 || (*state.globals).chaverede != 0) {
                    double rhoPSTF = (1 - state.cells[i].flui.BSW) * (1000 * 141.5 / (131.5 + state.cells[i].flui.API)) + state.cells[i].flui.BSW * 1000 * state.cells[i].flui.Denag;
                    fonteP *= ((1 - state.cells[i].flui.BSW) / rhoPSTF);
                    fonteA *= (state.cells[i].flui.BSW / rhoPSTF);
                    BSWF = state.cells[i].flui.BSW;
                    denagF = fluF.Denag;
                    viscLF = 0 * 30 + 1 * state.cells[i].flui.VisOM(tL);
                    viscHF = 0 * 20 + 1 * state.cells[i].flui.VisOM(tH);
                }
            } else if ((*state.cells[i].acsrL).tipo == 5 || (*state.cells[i].acsrL).tipo == 8) {
                if (i > 0) {
                    double rsF = state.cells[i - 1].flui.RS(state.cells[i].pres, state.cells[i].temp) * (6.29 / 35.31467);
                    if (state.input.nfluP > 0 || (*state.globals).chaverede != 0) {
                        double rhoPSTF = (1 - state.cells[i - 1].flui.BSW) * (1000 * 141.5 / (131.5 + state.cells[i - 1].flui.API)) +
                                         state.cells[i - 1].flui.BSW * 1000 *
                                             state.cells[i - 1].flui.Denag;
                        fonteP *= ((1 - state.cells[i - 1].flui.BSW) / rhoPSTF);
                        fonteA *= (state.cells[i - 1].flui.BSW / rhoPSTF);
                        BSWF = state.cells[i - 1].flui.BSW;
                        denagF = state.cells[i - 1].flui.Denag;
                        viscLF = 0 * 30 + 1 * state.cells[i - 1].flui.VisOM(tL);
                        viscHF = 0 * 20 + 1 * state.cells[i - 1].flui.VisOM(tH);
                    }
                } else {
                    double rsF = state.cells[i].flui.RS(state.cells[i].pres, state.cells[i].temp) * (6.29 / 35.31467);
                    if (state.input.nfluP > 0 || (*state.globals).chaverede != 0) {
                        double rhoPSTF = (1 - state.cells[i].flui.BSW) * (1000 * 141.5 / (131.5 + state.cells[i].flui.API)) +
                                         state.cells[i].flui.BSW * 1000 *
                                             state.cells[i].flui.Denag;
                        fonteP *= ((1 - state.cells[i].flui.BSW) / rhoPSTF);
                        fonteA *= (state.cells[i].flui.BSW / rhoPSTF);
                        BSWF = state.cells[i].flui.BSW;
                        denagF = state.cells[i - 1].flui.Denag;
                        viscLF = 0 * 30 + 1 * state.cells[i].flui.VisOM(tL);
                        viscHF = 0 * 20 + 1 * state.cells[i].flui.VisOM(tH);
                    }
                }
            }
        }
        if (titF > 1. - 1e-15) {
            fonteO = 0.;
            fonteP = 0.;
            fonteA = 0.;
        }

        double volpesFim = hol * (1 - bet) * (1 - bsw) / bo;
        double volaguaFim = hol * (1 - bet) * bsw / bo; // nao deveria ser dividido por Bo???????????
        double MultPe;
        double MultPd = 0.;
        double residuoP = 0.;
        double MultAe = 0.;
        double MultAd;
        double residuoA = 0.;
        if ((*state.globals).lixo5 < 1e-15 && state.input.snaps != 1) {
            state.cells[i].VolPesaST = volpesFim;
            state.cells[i].VolAguaST = volaguaFim;
        }
        if (state.input.nfluP > 1 || (*state.globals).chaverede != 0) {
            MultPe = 0.;
            if (state.cells[i].flui.dStockTankVaporMassFraction < 1. - 1e-15)
                MultPe = state.cells[i].QL * (1 - betI0) * (1 - bsw0) / bo0;
            MultPd = 0.;
            if (state.cells[i + 1].flui.dStockTankVaporMassFraction < 1. - 1e-15)
                MultPd = state.cells[i + 1].QL * (1 - betI1) * (1 - bsw1) / bo1;
            residuoP = (volpesFim - state.cells[i].VolPesaST) * at / dt + (MultPd - MultPe) / dx - fonteP / dx;
            MultAe = 0.;
            if (state.cells[i].flui.dStockTankVaporMassFraction < 1. - 1e-15)
                MultAe = state.cells[i].QL * (1 - betI0) * bsw0 / bo0;
            MultAd = 0.;
            if (state.cells[i + 1].flui.dStockTankVaporMassFraction < 1. - 1e-15)
                MultAd = state.cells[i + 1].QL * (1 - betI1) * bsw1 / bo1;
            residuoA = (volaguaFim - state.cells[i].VolAguaST) * at / dt + (MultAd - MultAe) / dx - fonteA / dx;
        }
        if ((state.input.nfluP > 1) || (*state.globals).chaverede != 0) {
            if (((volpesFim > 1e-3) && (fonteP > 0 || (MultPd < 0 || MultPe > 0))) && state.input.nfluP > 0) {
                VISCL[i] = (dt * (viscLF * fonteP / dx + 1. * viscLini * residuoP - (viscL1 * MultPd - viscL0 * MultPe) / dx) + viscLini * state.cells[i].VolPesaST * at) / (volpesFim * at - 0. * residuoP * dt);
                VISCH[i] = (dt * (viscHF * fonteP / dx + 1. * viscHini * residuoP - (viscH1 * MultPd - viscH0 * MultPe) / dx) + viscHini * state.cells[i].VolPesaST * at) / (volpesFim * at - 0. * residuoP * dt);
            } else {
                VISCL[i] = viscLini;
                VISCH[i] = viscHini;
            }
        }
        if (state.input.nfluP > 1 || (*state.globals).chaverede != 0) {
            if ((
                    (volaguaFim + volpesFim) > 1e-3) &&
                ((fonteP > 0 || fonteA > 0) ||
                 ((MultPd < 0 || MultPe > 0) || (MultAd < 0 || MultAe > 0)))) {
                BSW[i] = (dt * (BSWF * (fonteA + fonteP) / dx + 1. * BSWini * (residuoA + residuoP) - (BSW1 * (MultAd + MultPd) - BSW0 * (MultAe + MultPe)) / dx) + BSWini * (state.cells[i].VolAguaST + state.cells[i].VolPesaST) * at) / ((volaguaFim + volpesFim) * at -
                                                                                                                                                                                                                                  0. * dt * (residuoA + residuoP));
                denag[i] = (dt * (denagF * (fonteA) / dx + 1. * denagini * (residuoA) - (denag1 * (MultAd)-denag0 * (MultAe)) / dx) + denagini * (state.cells[i].VolAguaST) * at) / ((volaguaFim)*at -
                                                                                                                                                                                0. * dt * (residuoA));
                if (BSW[i] < 0.)
                    BSW[i] = 0.;
                else if (BSW[i] > 1.)
                    BSW[i] = 1.;
                else if (isnan(BSW[i]))
                    BSW[i] = BSWini;
                if (denag[i] < 1.)
                    denag[i] = 1.;
                else if (isnan(denag[i]))
                    denag[i] = denagini;
            } else
                BSW[i] = BSWini;
            denag[i] = denagini;
        }
        state.cells[i].VolPesaST = volpesFim - 0. * residuoP * dt / at;
        state.cells[i].VolAguaST = volaguaFim - 0. * residuoA * dt / at;
    }
    for (int i = imin; i <= state.lastCell - 1; i++) {
        if (state.input.nfluP > 0 || (*state.globals).chaverede != 0) {
            state.cells[i].flui = fluC[i];
            if (state.input.nfluP > 1 || (*state.globals).chaverede != 0) {
                state.cells[i].flui.BSW = BSW[i];
                state.cells[i].flui.Denag = denag[i];
                state.cells[i].flui.LVisL = VISCL[i];
                state.cells[i].flui.LVisH = VISCH[i];
                state.cells[i].flui.TempL = tL;
                state.cells[i].flui.TempH = tH;
            }
        }
    }

    if ((*state.globals).chaverede == 0 || state.endNode == 1 || state.cells[state.lastCell].Mliqini > -(*state.globals).localtiny) {
        if (state.input.nfluP > 0 || (*state.globals).chaverede != 0) {
            state.cells[state.lastCell].flui = state.cells[state.lastCell - 1].flui;
            if (state.input.nfluP > 1 || (*state.globals).chaverede != 0) {
                state.cells[state.lastCell].flui.BSW = BSW[state.lastCell - 1];
                state.cells[state.lastCell].flui.Denag = denag[state.lastCell - 1];
                state.cells[state.lastCell].flui.LVisL = VISCL[state.lastCell - 1];
                state.cells[state.lastCell].flui.LVisH = VISCH[state.lastCell - 1];
                state.cells[state.lastCell].flui.TempL = tL;
                state.cells[state.lastCell].flui.TempH = tH;
            }
        }
    } else {
        state.cells[state.lastCell].flui = fluiRev;
        state.cells[state.lastCell].flui.BSW = fluiRev.BSW;
        state.cells[state.lastCell].flui.Denag = fluiRev.Denag;
        state.cells[state.lastCell].flui.LVisL = fluiRev.LVisL;
        state.cells[state.lastCell].flui.LVisH = fluiRev.LVisH;
        state.cells[state.lastCell].flui.TempL = fluiRev.TempL;
        state.cells[state.lastCell].flui.TempH = fluiRev.TempH;
    }
    delete[] fluC;
    if(state.input.miniTabAtraso > 0)(*state.globals).modoTransiente = 1;
}

}  // namespace sisprod::composition
