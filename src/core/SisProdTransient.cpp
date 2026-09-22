#include "SisProdTransient.h"

#include "Acidentes2.h"
#include "Leitura.h"
#include "Matriz.h"
#include "Vetor.h"
#include "celula3.h"
#include "celulaGas.h"
#include "variaveisGlobais1D.h"

#include <math.h>

namespace sisprod::transient {

void updateInteriorCell(const TransientStepState &state, int i, int expli) {
    state.cells[i].presini = state.cells[i].pres;
    if (expli == 0)
        state.cells[i].pres = state.productionSolution[2 * i + 1];
    if (isnan(state.cells[i].pres))
        NumError("Pressao na linha com valor NaN");
    state.cells[i].d2pdt2 = state.cells[i].dpdt;
    state.cells[i].dpdt = 0. * (state.cells[i].pres - state.cells[i].presini) / state.cells[i].dt;
    state.cells[i].MCini = state.cells[i].MC;
    if (expli == 0)
        state.cells[i].MC = state.productionSolution[2 * i];
    if (isnan(state.cells[i].MC))
        NumError("Vazao massica da mistura na linha com valor NaN");
    state.cells[i + 1].presLini = state.cells[i + 1].presL;
    state.cells[i - 1].presRini = state.cells[i - 1].presR;
    state.cells[i + 1].presL = state.cells[i - 1].presR = state.cells[i].pres;
    state.cells[i + 1].MLini = state.cells[i + 1].ML;
    state.cells[i - 1].MRini = state.cells[i - 1].MR;
    state.cells[i + 1].ML = state.cells[i - 1].MR = state.cells[i].MC;

    state.cells[i].Mliqini0 = state.cells[i].Mliqini;
    state.cells[i].Mliqini = state.cells[i].term1 * state.cells[i].MC + state.cells[i].term2;
    state.cells[i + 1].MliqiniL0 = state.cells[i + 1].MliqiniL;
    state.cells[i - 1].MliqiniR0 = state.cells[i - 1].MliqiniR;
    state.cells[i + 1].MliqiniL = state.cells[i - 1].MliqiniR = state.cells[i].Mliqini;

    double dx = 0.5 * state.cells[i].dx;
    double dia = state.cells[i].duto.a;
    double area = 0.25 * M_PI * dia * dia;
    double si = state.cells[i].duto.peri;
    double alfmed = state.cells[i].alf;
    double rhog = state.cells[i].flui.MasEspGas(state.cells[i].pres, state.cells[i].temp);
    double rhol = (1 - state.cells[i].bet) * state.cells[i].flui.MasEspLiq(state.cells[i].pres, state.cells[i].temp) + state.cells[i].bet * state.cells[i].fluicol.MasEspFlu(state.cells[i].pres, state.cells[i].temp);
    double ugsmed = (state.cells[i].MC - state.cells[i].Mliqini) / (area * rhog);
    double ulsmed = state.cells[i].Mliqini / (area * rhol);
    double j = ugsmed + ulsmed;

    double rhomix = alfmed * rhog + (1 - alfmed) * rhol;
    double viscmix = alfmed * state.cells[i].flui.ViscGas(state.cells[i].pres, state.cells[i].temp) + (1 - alfmed) * ((1 - state.cells[i].bet) * state.cells[i].flui.ViscOleo(state.cells[i].pres, state.cells[i].temp) + state.cells[i].bet * state.cells[i].fluicol.VisFlu(state.cells[i].pres, state.cells[i].temp));

    double re1;
    if (state.cells[i].duto.revest == 0)
        re1 = state.cells[i].Rey(state.cells[i].duto.a, j, rhomix, viscmix);
    else {
        double dhid = 4 * area / si;
        re1 = state.cells[i].Rey(dhid, j, rhomix, viscmix);
    }
    double f1 = state.cells[i].fric(re1, state.cells[i].duto.rug / dia);
    double medpres = 0;
    if (state.input.MedSimpPresFront == 0) {
        if (state.cells[i].presaux <= 10)
            medpres = 1;
        else
            medpres = 0.;
    } else
        medpres = 1;
    double gradfric = (1 - medpres) * 0.5 * f1 * rhomix * (fabs(j) * j) * si * dx / area;
    double gradhidro = (1 - medpres) * 9.82 * sin(state.cells[i].duto.teta) * rhomix * dx;
    state.cells[i].presauxini = state.cells[i].presaux;
    state.cells[i].presaux = state.cells[i].pres + (gradfric + gradhidro - state.cells[i - 1].dpB) / 98066.5;
    state.cells[i].dpresaux = 0.5 * (gradfric + gradhidro - state.cells[i - 1].dpB) / 98066.5;
    dx = 0.5 * state.cells[i].dxL;
    dia = state.cells[i - 1].duto.a;
    area = 0.25 * M_PI * dia * dia;
    si = state.cells[i - 1].duto.peri;
    alfmed = state.cells[i - 1].alf;
    rhog = state.cells[i - 1].flui.MasEspGas(state.cells[i - 1].pres, state.cells[i - 1].temp);
    rhol = (1 - state.cells[i - 1].bet) * state.cells[i - 1].flui.MasEspLiq(state.cells[i - 1].pres, state.cells[i - 1].temp) + state.cells[i - 1].bet * state.cells[i - 1].fluicol.MasEspFlu(state.cells[i - 1].pres, state.cells[i - 1].temp);
    ugsmed = (state.cells[i - 1].MC - state.cells[i - 1].Mliqini) / (area * rhog);
    ulsmed = state.cells[i - 1].Mliqini / (area * rhol);
    j = ugsmed + ulsmed;

    rhomix = alfmed * rhog + (1 - alfmed) * rhol;
    viscmix = alfmed * state.cells[i - 1].flui.ViscGas(state.cells[i - 1].pres, state.cells[i - 1].temp) + (1 - alfmed) * ((1 - state.cells[i - 1].bet) * state.cells[i - 1].flui.ViscOleo(state.cells[i - 1].pres, state.cells[i - 1].temp) + state.cells[i - 1].bet * state.cells[i - 1].fluicol.VisFlu(state.cells[i - 1].pres, state.cells[i - 1].temp));

    if (state.cells[i - 1].duto.revest == 0)
        re1 = state.cells[i - 1].Rey(state.cells[i - 1].duto.a, j, rhomix, viscmix);
    else {
        double dhid = 4 * area / si;
        re1 = state.cells[i - 1].Rey(dhid, j, rhomix, viscmix);
    }
    f1 = state.cells[i - 1].fric(re1, state.cells[i - 1].duto.rug / dia);

    gradfric = (1 - medpres) * 0.5 * f1 * rhomix * (fabs(j) * j) * si * dx / area;
    gradhidro = (1 - medpres) * 9.82 * sin(state.cells[i - 1].duto.teta) * rhomix * dx;

    if (state.cells[i - 1].acsr.tipo != 5 || state.cells[i - 1].acsr.chk.AreaGarg > state.cells[i - 1].acsr.chk.AreaTub * 0.5)
        state.cells[i].presaux = 0.5 * (state.cells[i].presaux) +
                            0.5 * (state.cells[i - 1].pres - (gradfric + gradhidro) / 98066.5);
    state.cells[i].dpresaux -= 0.5 * (gradfric + gradhidro) / 98066.5;
    state.cells[i - 1].presauxRini = state.cells[i - 1].presauxR;
    state.cells[i - 1].presauxR = state.cells[i].presaux;
    if (i < state.lastCell) {
        state.cells[i + 1].presauxLini = state.cells[i + 1].presauxL;
        state.cells[i + 1].presauxL = state.cells[i].presaux;
    }

    double tmed = state.cells[i - 1].temp;
    if (state.cells[i].VTemper < 0.)
        tmed = state.cells[i].temp;

    ProFlu flud;
    if (state.cells[i].Mliqini < 0.)
        flud = state.cells[i].flui;
    else
        flud = state.cells[i - 1].flui;

    double betI;
    if (((state.cells[i].MC - state.cells[i].Mliqini) * 0 + 1 * state.cells[i].Mliqini) < 0)
        betI = state.cells[i].bet; // duvidabeta
    else
        betI = state.cells[i].betL;

    double rl = flud.MasEspLiq(state.cells[i].presaux, tmed);
    rhol = (1 - betI) * rl + betI * state.cells[i].fluicol.MasEspFlu(state.cells[i].presaux, tmed);

    double rg = flud.MasEspGas(state.cells[i].presaux, tmed);
    double vLiqTest = 1 + 0 * fabs(state.cells[i].Mliqini / (rhol * area));
    double vGasTest = 1 + 0 * fabs((state.cells[i].MC - state.cells[i].Mliqini) / (rg * area));

    state.cells[i].QLini = state.cells[i].QL;
    if (vLiqTest > 1e-3)
        state.cells[i].QL = state.cells[i].Mliqini / rhol;
    else {
        state.cells[i].QL = 0.;
        state.cells[i].MC = (state.cells[i].MC - state.cells[i].Mliqini);
        state.cells[i].Mliqini = 0;
        state.cells[i + 1].ML = state.cells[i - 1].MR = state.cells[i].MC;
        state.cells[i + 1].MliqiniL = state.cells[i - 1].MliqiniR = state.cells[i].Mliqini;
    }
    state.cells[i].QGini = state.cells[i].QG;
    if (vGasTest > 1e-3)
        state.cells[i].QG = (state.cells[i].MC - state.cells[i].Mliqini) / rg;
    else {
        state.cells[i].QG = 0;
        state.cells[i].MC = state.cells[i].Mliqini;
        state.cells[i + 1].ML = state.cells[i - 1].MR = state.cells[i].MC;
    }
    state.cells[i - 1].QLRini = state.cells[i - 1].QLR;
    state.cells[i - 1].QLR = state.cells[i].QL;
    if (i < state.lastCell) {
        state.cells[i + 1].QLLini = state.cells[i + 1].QLL;
        state.cells[i + 1].QLL = state.cells[i].QL;
    }

}

void updateFirstCell(const TransientStepState &state, int i, int expli) {
    state.cells[0].presini = state.cells[0].pres;
    if (expli == 0)
        state.cells[0].pres = state.productionSolution[1];
    state.cells[i].d2pdt2 = state.cells[i].dpdt;
    state.cells[i].dpdt = 0 * (state.cells[i].pres - state.cells[i].presini) / state.cells[i].dt;
    state.cells[0].presauxini = state.cells[0].presaux;
    state.cells[0].presaux = state.cells[0].pres;
    state.cells[0].dpresaux = 0.;
    state.cells[0].MCini = state.cells[0].MC;
    if (expli == 0)
        state.cells[0].MC = state.productionSolution[0];
    state.cells[1].presLini = state.cells[1].presL;
    state.cells[1].presL = state.cells[0].pres;
    state.cells[1].MLini = state.cells[1].ML;
    state.cells[1].ML = state.cells[0].MC;
    state.cells[0].Mliqini0 = state.cells[0].Mliqini;
    state.cells[0].Mliqini = state.cells[i].term1 * state.cells[i].MC + state.cells[i].term2;
    state.cells[1].MliqiniL0 = state.cells[1].MliqiniL;
    state.cells[1].MliqiniL = state.cells[0].Mliqini;

    state.cells[0].QLini = state.cells[0].QL;
    state.cells[0].QL = 0.;
    state.cells[1].QLLini = state.cells[1].QLL;
    state.cells[1].QLL = state.cells[0].QL;
    state.cells[0].QGini = state.cells[0].QG;
    state.cells[0].QG = 0.;

    if (state.input.ConContEntrada > 0) {
        double rhogC = state.cells[i].flui.MasEspGas(state.inletPressure, state.inletTemperature);
        double rhopC = state.cells[i].flui.MasEspLiq(state.inletPressure, state.inletTemperature);
        double rhocC = state.cells[i].fluicol.MasEspFlu(state.inletPressure, state.inletTemperature);
        double rholC = rhopC * (1 - state.inletCompletionFraction) + rhocC * state.inletCompletionFraction;
        state.cells[0].QL = state.cells[0].Mliqini / rholC;
        state.cells[1].QLL = state.cells[0].QL;
        state.cells[0].QG = (state.cells[0].MC - state.cells[0].Mliqini) / rhogC;
    }
}

void updateLastCell(const TransientStepState &state, int i, int expli) {
    state.cells[i].presini = state.cells[i].pres;
    if (expli == 0)
        state.cells[state.lastCell].pres = state.productionSolution[2 * state.lastCell + 1];
    state.cells[i].d2pdt2 = state.cells[i].dpdt;
    state.cells[i].dpdt = 0 * (state.cells[i].pres - state.cells[i].presini) / state.cells[i].dt;
    state.cells[i].d2pdt2 = (state.cells[i].dpdt - state.cells[i].d2pdt2) / state.cells[i].dt;
    state.cells[state.lastCell].MCini = state.cells[state.lastCell].MC;
    if (expli == 0)
        state.cells[state.lastCell].MC = state.productionSolution[2 * state.lastCell];
    // teste!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // teste!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    state.cells[state.lastCell - 1].presRini = state.cells[state.lastCell - 1].presR;
    state.cells[state.lastCell - 1].presR = state.cells[state.lastCell].pres;
    state.cells[state.lastCell - 1].MRini = state.cells[state.lastCell - 1].MR;
    state.cells[state.lastCell - 1].MR = state.cells[state.lastCell].MC;
    state.cells[state.lastCell].MRini = state.cells[state.lastCell].MR;
    state.cells[state.lastCell].MR = state.cells[state.lastCell].MC;

    state.cells[state.lastCell].Mliqini0 = state.cells[state.lastCell].Mliqini;
    state.cells[state.lastCell].Mliqini = state.cells[state.lastCell].term1 * state.cells[state.lastCell].MC + state.cells[state.lastCell].term2;
    state.cells[i - 1].MliqiniR0 = state.cells[i - 1].MliqiniR;
    state.cells[state.lastCell - 1].MliqiniR = state.cells[state.lastCell].Mliqini;
    state.cells[state.lastCell].MliqiniR0 = state.cells[state.lastCell].MliqiniR;
    state.cells[state.lastCell].MliqiniR = state.cells[state.lastCell].Mliqini;

    double dx = 0.5 * state.cells[i].dx;
    double dia = state.cells[i].duto.a;
    double area = 0.25 * M_PI * dia * dia;
    double si = state.cells[i].duto.peri;
    double alfmed = state.cells[i].alf;
    double rhog = state.cells[i].flui.MasEspGas(state.cells[i].pres, state.cells[i].temp);
    double rhol = (1 - state.cells[i].bet) * state.cells[i].flui.MasEspLiq(state.cells[i].pres, state.cells[i].temp) + state.cells[i].bet * state.cells[i].fluicol.MasEspFlu(state.cells[i].pres, state.cells[i].temp);
    double ugsmed = (state.cells[i].QG) / (area);
    double ulsmed = state.cells[i].QL / (area);
    double j = ugsmed + ulsmed;
    double ABSjL = (fabs(state.cells[i - 1].QG) + fabs(state.cells[i - 1].QL)) / state.cells[i - 1].duto.area;

    double rhomix = alfmed * rhog + (1 - alfmed) * rhol;
    double viscmix = alfmed * state.cells[i].flui.ViscGas(state.cells[i].pres, state.cells[i].temp) + (1 - alfmed) * ((1 - state.cells[i].bet) * state.cells[i].flui.ViscOleo(state.cells[i].pres, state.cells[i].temp) + state.cells[i].bet * state.cells[i].fluicol.VisFlu(state.cells[i].pres, state.cells[i].temp));

    double re1;
    if (state.cells[i].duto.revest == 0)
        re1 = state.cells[i].Rey(state.cells[i].duto.a, j, rhomix, viscmix);
    else {
        double dhid = 4 * area / si;
        re1 = state.cells[i].Rey(dhid, j, rhomix, viscmix);
    }
    double f1 = state.cells[i].fric(re1, state.cells[i].duto.rug / dia);
    double medpres = 0;
    double gradfric = (1 - medpres) * 0.5 * f1 * rhomix * (fabs(j) * j) * si * dx / area;
    double gradhidro = (1 - medpres) * 9.82 * sin(state.cells[i].duto.teta) * rhomix * dx;

    state.cells[i].presauxini = state.cells[i].presaux;
    state.cells[i].presaux = state.cells[i].pres + (gradfric + gradhidro - state.cells[i - 1].dpB) / 98066.5;
    state.cells[i].dpresaux = 0.5 * (gradfric + gradhidro - state.cells[i - 1].dpB) / 98066.5;
    dx = 0.5 * state.cells[i].dxL;
    dia = state.cells[i - 1].duto.a;
    area = 0.25 * M_PI * dia * dia;
    si = state.cells[i - 1].duto.peri;
    alfmed = state.cells[i - 1].alf;
    rhog = state.cells[i - 1].flui.MasEspGas(state.cells[i - 1].pres, state.cells[i - 1].temp);
    rhol = (1 - state.cells[i - 1].bet) * state.cells[i - 1].flui.MasEspLiq(state.cells[i - 1].pres, state.cells[i - 1].temp) + state.cells[i - 1].bet * state.cells[i - 1].fluicol.MasEspFlu(state.cells[i - 1].pres, state.cells[i - 1].temp);
    ugsmed = (state.cells[i - 1].MC - state.cells[i - 1].Mliqini) / (area * rhog);
    ulsmed = state.cells[i - 1].Mliqini / (area * rhol);
    j = ugsmed + ulsmed;

    rhomix = alfmed * rhog + (1 - alfmed) * rhol;
    viscmix = alfmed * state.cells[i - 1].flui.ViscGas(state.cells[i - 1].pres, state.cells[i - 1].temp) + (1 - alfmed) * ((1 - state.cells[i - 1].bet) * state.cells[i - 1].flui.ViscOleo(state.cells[i - 1].pres, state.cells[i - 1].temp) + state.cells[i - 1].bet * state.cells[i - 1].fluicol.VisFlu(state.cells[i - 1].pres, state.cells[i - 1].temp));

    if (state.cells[i - 1].duto.revest == 0)
        re1 = state.cells[i - 1].Rey(state.cells[i - 1].duto.a, j, rhomix, viscmix);
    else {
        double dhid = 4 * area / si;
        re1 = state.cells[i - 1].Rey(dhid, j, rhomix, viscmix);
    }
    f1 = state.cells[i - 1].fric(re1, state.cells[i - 1].duto.rug / dia);
    gradfric = (1 - medpres) * 0.5 * f1 * rhomix * (fabs(j) * j) * si * dx / area;
    gradhidro = (1 - medpres) * 9.82 * sin(state.cells[i - 1].duto.teta) * rhomix * dx;

    state.cells[i].presaux = 0.5 * (state.cells[i].presaux) +
                        0.5 * (state.cells[i - 1].pres - (gradfric + gradhidro) / 98066.5);
    state.cells[i].dpresaux -= 0.5 * (gradfric + gradhidro) / 98066.5;
    state.cells[i - 1].presauxRini = state.cells[i - 1].presauxR;
    state.cells[i - 1].presauxR = state.cells[i].presaux;
    if (i < state.lastCell) {
        state.cells[i + 1].presauxLini = state.cells[i + 1].presauxL;
        state.cells[i + 1].presauxL = state.cells[i].presaux;
    }

    double tmed = state.cells[i - 1].temp;
    if (state.cells[i].VTemper < 0.)
        tmed = state.cells[i].temp;

    ProFlu flud;
    if (state.cells[i].Mliqini < 0.)
        flud = state.cells[i].flui;
    else
        flud = state.cells[i - 1].flui;

    double betI;
    if (((state.cells[i].MC - state.cells[i].Mliqini) * 0.99 + 0.01 * state.cells[i].Mliqini) < 0)
        betI = state.cells[i].bet; // duvidabeta
    else
        betI = state.cells[i].betL;

    double rl = flud.MasEspLiq(state.cells[i].presaux, tmed);
    rhol = (1 - betI) * rl + betI * state.cells[i].fluicol.MasEspFlu(state.cells[i].presaux, tmed);

    double rg = flud.MasEspGas(state.cells[i].presaux, tmed);

    state.cells[i].QL = state.cells[i].Mliqini / rhol;
    state.cells[i].QG = (state.cells[i].MC - state.cells[i].Mliqini) / rg;
    state.cells[i - 1].QLR = state.cells[i].QL;
}

void updateCells(const TransientStepState &state, int expli) {
    for (int i = 0; i <= state.lastCell; i++) {
        if (i != 0 && i != state.lastCell) {
            updateInteriorCell(state, i, expli);
        } else if (i == 0) {
            updateFirstCell(state, i, expli);
        } else {
            updateLastCell(state, i, expli);
        }
    }
}

void updateFlowRates(const TransientStepState &state) {
    for (int i = 0; i <= state.lastCell; i++) {
        if (i != 0 && i != state.lastCell) {
            state.cells[i].Mliqini = state.cells[i].term1 * state.cells[i].MC + state.cells[i].term2;
            state.cells[i + 1].MliqiniL = state.cells[i - 1].MliqiniR = state.cells[i].Mliqini;

            double tmed = state.cells[i - 1].temp;
            if (state.cells[i].VTemper < 0.)
                tmed = state.cells[i].temp;
            ProFlu flud;
            if (state.cells[i].Mliqini < 0.)
                flud = state.cells[i].flui;
            else
                flud = state.cells[i - 1].flui;

            double betI;
            if (((state.cells[i].MC - state.cells[i].Mliqini) * 0 + 1 * state.cells[i].Mliqini) < 0)
                betI = state.cells[i].bet; // duvidabeta
            else
                betI = state.cells[i].betL;

            double rl = flud.MasEspLiq(state.cells[i].presaux, tmed);
            double rhol = (1 - betI) * rl + betI * state.cells[i].fluicol.MasEspFlu(state.cells[i].presaux, tmed);

            double rg = flud.MasEspGas(state.cells[i].presaux, tmed);
            state.cells[i].QL = state.cells[i].Mliqini / rhol;
            state.cells[i].QG = (state.cells[i].MC - state.cells[i].Mliqini) / rg;
            state.cells[i - 1].QLR = state.cells[i].QL;
            if (i < state.lastCell) {
                state.cells[i + 1].QLL = state.cells[i].QL;
            }

        } else if (i == 0) {

            state.cells[0].Mliqini = state.cells[i].term1 * state.cells[i].MC + state.cells[i].term2;
            state.cells[1].MliqiniL = state.cells[0].Mliqini;

            if (state.input.ConContEntrada > 0) {
                double rhogC = state.cells[i].flui.MasEspGas(state.inletPressure, state.inletTemperature);
                double rhopC = state.cells[i].flui.MasEspLiq(state.inletPressure, state.inletTemperature);
                double rhocC = state.cells[i].fluicol.MasEspFlu(state.inletPressure, state.inletTemperature);
                double rholC = rhopC * (1 - state.inletCompletionFraction) + rhocC * state.inletCompletionFraction;
                state.cells[0].QL = state.cells[0].Mliqini / rholC;
                state.cells[1].QLL = state.cells[0].QL;
                state.cells[0].QG = (state.cells[0].MC - state.cells[0].Mliqini) / rhogC;
            }
        } else {

            state.cells[state.lastCell].Mliqini = state.cells[state.lastCell].term1 * state.cells[state.lastCell].MC + state.cells[state.lastCell].term2;
            state.cells[state.lastCell - 1].MliqiniR = state.cells[state.lastCell].Mliqini;
            state.cells[state.lastCell].MliqiniR = state.cells[state.lastCell].Mliqini;

            double tmed = state.cells[i - 1].temp;
            if (state.cells[i].VTemper < 0.)
                tmed = state.cells[i].temp;

            ProFlu flud;
            if (state.cells[i].Mliqini < 0.)
                flud = state.cells[i].flui;
            else
                flud = state.cells[i - 1].flui;

            double betI;
            if (((state.cells[i].MC - state.cells[i].Mliqini) * 0.99 + 0.01 * state.cells[i].Mliqini) < 0)
                betI = state.cells[i].bet; // duvidabeta
            else
                betI = state.cells[i].betL;

            double rl = flud.MasEspLiq(state.cells[i].presaux, tmed);
            double rhol = (1 - betI) * rl + betI * state.cells[i].fluicol.MasEspFlu(state.cells[i].presaux, tmed);

            double rg = flud.MasEspGas(state.cells[i].presaux, tmed);

            state.cells[i].QL = state.cells[i].Mliqini / rhol;
            state.cells[i].QG = (state.cells[i].MC - state.cells[i].Mliqini) / rg;
            state.cells[i - 1].QLR = state.cells[i].QL;
        }
    }
}

void updateBufferFromSolution(const TransientStepState &state) {

    state.cells[0].presBuf = state.productionSolution[1];
    state.cells[0].presauxBuf = state.cells[0].presBuf;
    state.cells[0].MCBuf = state.productionSolution[0];
    state.cells[1].presLiniBuf = state.cells[1].presLBuf;
    state.cells[1].presLBuf = state.cells[0].presBuf;
    state.cells[1].MLBuf = state.cells[0].MCBuf;
    state.cells[0].MliqiniBuf = state.cells[0].term1 * state.cells[0].MCBuf + state.cells[0].term2;
    state.cells[1].MliqiniLBuf = state.cells[0].MliqiniBuf;

    int i = 1;
    state.cells[i].presBuf = state.productionSolution[2 * i + 1];
    state.cells[i].MCBuf = state.productionSolution[2 * i];
    state.cells[i + 1].presLiniBuf = state.cells[i + 1].presLBuf;
    state.cells[i + 1].presLBuf = state.cells[i - 1].presRBuf = state.cells[i].presBuf;
    state.cells[i - 1].MRiniBuf = state.cells[i - 1].MRBuf;
    state.cells[i + 1].MLBuf = state.cells[i - 1].MRBuf = state.cells[i].MCBuf;

    state.cells[i].MliqiniBuf = state.cells[i].term1 * state.cells[i].MCBuf + state.cells[i].term2;
    state.cells[i + 1].MliqiniLBuf = state.cells[i - 1].MliqiniRBuf = state.cells[i].MliqiniBuf;

    state.cells[state.lastCell].presBuf = state.productionSolution[2 * state.lastCell + 1];
    state.cells[state.lastCell].MCBuf = state.productionSolution[2 * state.lastCell];
    state.cells[state.lastCell - 1].presRBuf = state.cells[state.lastCell].presBuf;
    state.cells[state.lastCell - 1].MRBuf = state.cells[state.lastCell].MCBuf;
    state.cells[state.lastCell].MRBuf = state.cells[state.lastCell].MCBuf;

    state.cells[state.lastCell].MliqiniBuf = state.cells[state.lastCell].term1 * state.cells[state.lastCell].MCBuf + state.cells[state.lastCell].term2;
    state.cells[state.lastCell - 1].MliqiniRBuf = state.cells[state.lastCell].MliqiniBuf;
    state.cells[state.lastCell].MliqiniRBuf = state.cells[state.lastCell].MliqiniBuf;
}

void updateBufferFromCells(const TransientStepState &state) {

    state.cells[0].presBuf = state.cells[0].pres;
    state.cells[0].presauxBuf = state.cells[0].presaux;
    state.cells[0].MCBuf = state.cells[0].MC;
    state.cells[1].presLiniBuf = state.cells[1].presLini;
    state.cells[1].presLBuf = state.cells[1].presL;
    state.cells[1].MLBuf = state.cells[1].ML;
    state.cells[0].MliqiniBuf = state.cells[0].Mliqini;
    state.cells[1].MliqiniLBuf = state.cells[1].MliqiniL;

    int i = 1;
    state.cells[i].presBuf = state.cells[i].pres;
    state.cells[i].MCBuf = state.cells[i].MC;
    state.cells[i + 1].presLiniBuf = state.cells[i + 1].presLini;
    state.cells[i + 1].presLBuf = state.cells[i - 1].presRBuf = state.cells[i].pres;
    state.cells[i - 1].MRiniBuf = state.cells[i - 1].MR;
    state.cells[i + 1].MLBuf = state.cells[i - 1].MRBuf = state.cells[i].MC;

    state.cells[i].MliqiniBuf = state.cells[i].Mliqini;
    state.cells[i + 1].MliqiniLBuf = state.cells[i - 1].MliqiniRBuf = state.cells[i].Mliqini;

    state.cells[state.lastCell].presBuf = state.cells[state.lastCell].pres;
    state.cells[state.lastCell].MCBuf = state.cells[state.lastCell].MC;
    state.cells[state.lastCell - 1].presRBuf = state.cells[state.lastCell].pres;
    state.cells[state.lastCell - 1].MRBuf = state.cells[state.lastCell].MC;

    state.cells[state.lastCell].MliqiniBuf = state.cells[state.lastCell].Mliqini;
    state.cells[state.lastCell - 1].MliqiniRBuf = state.cells[state.lastCell].Mliqini;

    int fim = state.lastCell - 1;
    state.bufferedCompletionMassSource = state.cells[fim + 1].fontemassCR;

    state.bufferedLiquidMassSource = state.cells[fim + 1].fontemassLR;

    state.bufferedGasMassSource = state.cells[fim + 1].fontemassGR;
}

bool surfaceChokeIsOpen(const TransientStepState &state) {
    return state.surfaceChoke.AreaGarg > (1e-3) * state.cells[state.lastCell - 1].duto.area;
}

bool surfaceChokeIsShut(const TransientStepState &state) {
    return state.surfaceChoke.AreaGarg < (1e-3) * state.cells[state.lastCell - 1].duto.area;
}

void applyOutletPressureCondition(const TransientStepState &state, double titRev, double alfRev, double betRev) {

    double tESup = state.cells[state.lastCell].temp;
    double alfSup = state.cells[state.lastCell].alf;
    double betSup = state.cells[state.lastCell].bet;

    double masentrada = state.cells[state.lastCell - 1].MR;
    double massgas = state.cells[state.lastCell - 1].MR - state.cells[state.lastCell - 1].MliqiniR;
    double maxSup = 0.;
    double chokemas = 0;

    double rholp = state.cells[state.lastCell].rpC;
    double rholc = state.cells[state.lastCell].rcC;
    double rholmix = (1 - betSup) * rholp + betSup * rholc;
    double romix = alfSup * state.cells[state.lastCell].rgC + (1 - alfSup) * rholmix;

    double tit;
    if ((massgas >= 0 && state.cells[state.lastCell - 1].MliqiniR <= 0) || (massgas < 0 && state.cells[state.lastCell - 1].MliqiniR == 0))
        tit = 1.;
    else if (massgas <= 0 && state.cells[state.lastCell - 1].MliqiniR > 0)
        tit = 0.;
    else if (masentrada < 0) {
        tit = 1.;
    } else if (fabs(masentrada) < 1e-15)
        tit = alfSup * state.cells[state.lastCell].flui.MasEspGas(state.cells[state.lastCell].pres, state.cells[state.lastCell].temp) / romix;
    else
        tit = fabs(massgas / masentrada);
    if (tit > 1)
        tit = 1;
    if (state.finalPressure < state.gasSurfacePressure) {
        tit = 1.;
    }

    if (tit == 0 && alfSup > 0.05)
        tit = alfSup * state.cells[state.lastCell].flui.MasEspGas(state.cells[state.lastCell].pres, state.cells[state.lastCell].temp) / romix;

    romix = tit * (1. / state.cells[state.lastCell].rgC) + (1 - tit) * (1. / rholmix);
    romix = 1 / romix;

    double masChk;

    double sinal = 1.;
    double pmon = state.finalPressure;

    double ypres = state.gasSurfacePressure / state.finalPressure;
    if (surfaceChokeIsOpen(state) && ypres < 1.) {
        double cplM = (1. - betSup) * state.cells[state.lastCell].flui.CalorLiq(state.finalPressure, tESup) -
                      betSup * state.cells[state.lastCell].fluicol.CalorLiq(state.finalPressure, tESup);
        double jtlM = (1. - betSup) * state.cells[state.lastCell].flui.JTL(state.finalPressure, tESup) - betSup / rholc;
        double cpg = state.cells[state.lastCell].flui.CalorGas(state.finalPressure, tESup);
        double jtgM = state.cells[state.lastCell].flui.JTG(state.finalPressure, tESup);
        state.input.valTempChokeJus = tESup + ((1. - tit) * jtlM / cplM + tit * jtgM / cpg) * (state.gasSurfacePressure - state.finalPressure) * 98066.52;
    }
    if (ypres > 1.) {
        if (state.input.chkv == 0)
            sinal = -1.;
        else
            sinal = 0.;
        tit = 1.;
        pmon = state.gasSurfacePressure;
        ypres = 1. / ypres;
    }

    masChk = state.surfaceChoke.vazmassSachd(ypres, pmon, tESup, alfSup, betSup, tit, state.cells[state.lastCell - 1].flui,
                                   state.cells[state.lastCell - 1].fluicol);
    maxSup = state.surfaceChoke.vazmaxSachd(pmon, tESup, alfSup, betSup, tit, state.cells[state.lastCell - 1].flui, state.cells[state.lastCell - 1].fluicol);

    int fluxcri = 1;
    if (tit <= 0.01 || fabs(ypres) > fabs(state.surfaceChoke.razpres)) {
        fluxcri = 0;
        maxSup = masChk;
    }
    if (surfaceChokeIsShut(state))
        maxSup = 0.;

    double sinal2 = 1.;
    pmon = state.finalPressure * 1.0001;
    ypres = state.gasSurfacePressure / pmon;
    if (ypres > 1.) {
        sinal2 = -1.;
        tit = 1;
        pmon = state.gasSurfacePressure;
        ypres = 1. / ypres;
        if (state.input.chkv == 0)
            sinal2 = -1.;
        else
            sinal2 = 0.;
    }

    double masChk2 = state.surfaceChoke.vazmassSachd(ypres, pmon, tESup, alfSup, betSup, tit, state.cells[state.lastCell - 1].flui,
                                           state.cells[state.lastCell - 1].fluicol);
    double maxSup2 = state.surfaceChoke.vazmaxSachd(pmon, tESup, alfSup, betSup, tit, state.cells[state.lastCell - 1].flui,
                                          state.cells[state.lastCell - 1].fluicol);

    fluxcri = 1;
    if (tit <= 0.01 || fabs(ypres) > fabs(state.surfaceChoke.razpres)) {
        fluxcri = 0;
        maxSup2 = masChk2;
    }
    if (surfaceChokeIsShut(state))
        maxSup2 = 0.;
    double dmaxsup = (sinal2 * maxSup2 - sinal * maxSup) / (state.finalPressure * 0.0001);

    double masliq;
    double masgas;
    int abertoini = state.open;
    double delp;
    if (surfaceChokeIsOpen(state))
        delp = (0.5 / 98066.5) * (1 / romix) * (1 / (state.surfaceChoke.AreaGarg * state.surfaceChoke.AreaGarg * state.surfaceChoke.cdchk * state.surfaceChoke.cdchk)) * masentrada * masentrada;
    else
        delp = 0;

    int masChkSup0 = state.surfaceChokeMassFlag;
    state.chokeModeChanged = 0;

    double difdelp = state.finalPressure - state.gasSurfacePressure;
    if (((tit < 1e-7 && surfaceChokeIsOpen(state)) ||
         (tit < 0.01 && surfaceChokeIsOpen(state) &&
          fabs(difdelp) / delp < 1.2 && fabs(difdelp) / delp > 0.8 &&
          ((fabs(maxSup) > 0 && fabs((masentrada - maxSup) / maxSup) < 0.2) ||
           (fabs(masentrada) > 0 && fabs((masentrada - maxSup) / masentrada) < 0.2))))) {
        state.open = 1;
        double sens = 1.;
        if (masentrada < 0. && state.cells[state.lastCell - 1].MliqiniR <= 0)
            sens = 0.;
        state.finalPressure = state.gasSurfacePressure + sens * delp;
        state.surfaceChokeMassFlag = 0;
        if (state.surfaceChokeMassFlag != masChkSup0)
            state.chokeModeChanged = 1;

    } else {
        if (surfaceChokeIsOpen(state) &&
            ((((*state.globals).lixo5 - 2 * state.input.dtmax > state.movingMeanTemperature || state.input.perm == 2) &&
              fabs(state.movingMeanPressure - state.gasSurfacePressure) / state.movingMeanPressure < 0.05 && fabs(state.movingMeanFlux) < 0.5) ||
             (((*state.globals).lixo5 - 2 * state.input.dtmax > state.movingMeanTemperature || state.input.perm == 2) &&
              fabs(state.movingMeanPressure - state.gasSurfacePressure) < (0.05 * state.gasSurfacePressure) && fabs(state.movingMeanFlux) < 5. && delp < 0.01 * state.gasSurfacePressure) ||
             (((*state.globals).lixo5 - 2 * state.input.dtmax > state.movingMeanTemperature || state.input.perm == 2) && (state.gasSurfacePressure - state.movingMeanPressure) / state.movingMeanPressure > 0.001 && state.input.chkv == 0) || (fabs(delp) < 0.1 && (state.finalPressure - state.gasSurfacePressure) / state.finalPressure < 0.05 && state.input.chkv == 0))) {
            state.open = 1;
            if (abertoini != state.open)
                state.openTime = 1;
        } else {
            state.open = 0;
            if (state.openTime > 60)
                state.openTime = 0;
        }
        if ((((tit > -0.01 && state.cells[state.lastCell].alf > -0.01) || surfaceChokeIsShut(state)) && state.surfaceChoke.AreaGarg < 0.6 * state.cells[state.lastCell - 1].duto.area && (state.open == 0 && (state.openTime == 0 || state.openTime > 60)))) {

            state.open = 0;
            state.openTime = 0;
            masliq = sinal * maxSup * (1. - tit);
            masgas = sinal * maxSup * tit;

            state.cells[state.lastCell].DmasschokeG = -1 * (1. - tit) * dmaxsup;
            state.cells[state.lastCell].DmasschokeL = -1 * tit * ((1 - betSup) * rholp / rholmix) * dmaxsup;
            state.cells[state.lastCell].DmasschokeC = -1 * tit * (betSup * rholc / rholmix) * dmaxsup;

            state.surfaceChokeMassFlag = 1;
            if (state.surfaceChokeMassFlag != masChkSup0)
                state.chokeModeChanged = 1;
            state.cells[state.lastCell].fontemassLR = -masliq * (1 - betSup) * rholp / rholmix;
            state.cells[state.lastCell].fontemassCR = -masliq * betSup * rholc / rholmix;
            state.cells[state.lastCell].fontemassGR = -masgas;

        } else {
            if (state.openTime != 0)
                state.openTime++;
            if (state.surfaceChoke.AreaGarg >= 0.601 * state.cells[state.lastCell - 1].duto.area)
                state.finalPressure = state.gasSurfacePressure;
            else {

                state.open = 1;
                state.surfaceChokeMassFlag = 0;
                if (state.surfaceChokeMassFlag != masChkSup0)
                    state.chokeModeChanged = 1;
                state.finalPressure = state.gasSurfacePressure;
            }
            state.open = 1;
            state.surfaceChokeMassFlag = 0;
            if (state.surfaceChokeMassFlag != masChkSup0)
                state.chokeModeChanged = 1;
        }
    }

    if (state.surfaceChokeMassFlag == 0 && (*state.globals).chaverede == 1) {
        double betloc;
        if ((state.cells[state.lastCell - 1].MR - state.cells[state.lastCell - 1].MliqiniR) * 0 + 1 * state.cells[state.lastCell - 1].MliqiniR > 0.)
            betloc = state.cells[state.lastCell - 1].bet; // testeBeta//duvidabeta
        else
            betloc = betRev;

        if (state.input.chkv == 0 || state.cells[state.lastCell - 1].MR > 0.)
            sinal = 1.;
        else
            sinal = 0.;

        state.cells[state.lastCell].fontemassCR = -sinal * state.cells[state.lastCell - 1].QLR * (betloc)*state.cells[state.lastCell - 1].rcC;
        state.cells[state.lastCell].fontemassLR = -sinal * (state.cells[state.lastCell - 1].MliqiniR + state.cells[state.lastCell].fontemassCR);
        state.cells[state.lastCell].fontemassGR = -sinal * (state.cells[state.lastCell - 1].MR - state.cells[state.lastCell - 1].MliqiniR);
    }
}

void applyOutletBufferCondition(const TransientStepState &state, double titRev, double alfRev, double betRev) {

    double tESup = state.cells[state.lastCell].temp;
    double alfSup = state.cells[state.lastCell].alf;
    double betSup = state.cells[state.lastCell].bet;

    double masentrada = state.cells[state.lastCell - 1].MRBuf;
    double massgas = state.cells[state.lastCell - 1].MRBuf - state.cells[state.lastCell - 1].MliqiniRBuf;
    double maxSup = 0.;

    double rholp = state.cells[state.lastCell].flui.MasEspLiq(state.cells[state.lastCell].presBuf, state.cells[state.lastCell].temp);
    double rholc = state.cells[state.lastCell].fluicol.MasEspFlu(state.cells[state.lastCell].presBuf, state.cells[state.lastCell].temp);
    double rholmix = (1 - betSup) * rholp + betSup * rholc;
    double romix = alfSup * state.cells[state.lastCell].flui.MasEspGas(state.cells[state.lastCell].presBuf, state.cells[state.lastCell].temp) + (1 - alfSup) * rholmix;

    double tit;
    if (massgas > 0 && state.cells[state.lastCell - 1].MliqiniRBuf < 0)
        tit = 1.;
    else if (massgas <= 0 && state.cells[state.lastCell - 1].MliqiniRBuf >= 0)
        tit = 0.;
    else if (masentrada < 0) {
        if ((*state.globals).chaverede == 0 || state.endNode == 1)
            tit = 1.;
        else {
            tit = titRev;
            alfSup = alfRev;
            betSup = betRev;
        }
    } else if (fabs(masentrada) < 1e-15)
        tit = alfSup * state.cells[state.lastCell].flui.MasEspGas(state.cells[state.lastCell].presBuf, state.cells[state.lastCell].temp) / romix;
    else
        tit = fabs(massgas / masentrada);
    if (tit > 1)
        tit = 1;
    if (state.cells[state.lastCell].presBuf < state.gasSurfacePressure) {
        if ((*state.globals).chaverede == 0 || state.endNode == 1)
            tit = 1.;
        else {
            tit = titRev;
            alfSup = alfRev;
            betSup = betRev;
        }
    }

    if (tit == 0 && alfSup > 0.05)
        tit = alfSup * state.cells[state.lastCell].flui.MasEspGas(state.cells[state.lastCell].presBuf, state.cells[state.lastCell].temp) / romix;

    romix = tit * (1. / state.cells[state.lastCell].rgC) + (1 - tit) * (1. / rholmix);
    romix = 1 / romix;

    double masChk;

    double sinal = 1.;
    double pmon = state.cells[state.lastCell].presBuf;

    double ypres = state.gasSurfacePressure / state.cells[state.lastCell].presBuf;
    if (ypres > 1.) {
        if (state.input.chkv == 0)
            sinal = -1.;
        else
            sinal = 0.;
        if ((*state.globals).chaverede == 0 || state.endNode == 1)
            tit = 1.;
        else {
            tit = titRev;
            alfSup = alfRev;
            betSup = betRev;
        }
        pmon = state.gasSurfacePressure;
        ypres = 1. / ypres;
    }

    masChk = state.surfaceChoke.vazmassSachd(ypres, pmon, tESup, alfSup, betSup, tit, state.cells[state.lastCell - 1].flui,
                                   state.cells[state.lastCell - 1].fluicol);
    maxSup = state.surfaceChoke.vazmaxSachd(pmon, tESup, alfSup, betSup, tit, state.cells[state.lastCell - 1].flui,
                                  state.cells[state.lastCell - 1].fluicol);

    int fluxcri = 1;
    if (tit <= 0.01 || fabs(ypres) > fabs(state.surfaceChoke.razpres)) {
        fluxcri = 0;
        maxSup = masChk;
    }
    if (surfaceChokeIsShut(state))
        maxSup = 0.;

    double masliq;
    double masgas;
    int abertoini = state.open;

    double delp;
    if (surfaceChokeIsOpen(state))
        delp = (0.5 / 98066.5) * (1 / romix) *
               (1 / (state.surfaceChoke.AreaGarg * state.surfaceChoke.AreaGarg * state.surfaceChoke.cdchk * state.surfaceChoke.cdchk)) * masentrada * masentrada;
    else
        delp = 0.;

    double difdelp = fabs(fabs(state.cells[state.lastCell].presBuf - state.gasSurfacePressure) - delp);
    if (((tit < 1e-7 && surfaceChokeIsOpen(state)) ||
         (tit < 0.01 && surfaceChokeIsOpen(state) &&
          difdelp / delp < 0.2 &&
          ((fabs(maxSup) > 0 && fabs((masentrada - maxSup) / maxSup) < 0.2) ||
           (fabs(masentrada) > 0 && fabs((masentrada - maxSup) / masentrada) < 0.2))))) {
        state.open = 1;
        if (state.input.chkv == 0 || state.cells[state.lastCell - 1].MRBuf > 0.)
            sinal = 1.;
        else
            sinal = 0.;
        double betloc;
        if ((state.cells[state.lastCell - 1].MRBuf - state.cells[state.lastCell - 1].MliqiniRBuf) * 0 + 1 * state.cells[state.lastCell - 1].MliqiniRBuf > 0.)
            betloc = state.cells[state.lastCell - 1].bet; // testeBeta//duvidabeta
        else
            betloc = betRev;
        double rhomistBuf = betloc *
                                state.cells[state.lastCell - 1].fluicol.MasEspFlu(state.cells[state.lastCell].presBuf, state.cells[state.lastCell].temp) +
                            (1. - betloc) * state.cells[state.lastCell - 1].flui.MasEspLiq(state.cells[state.lastCell].presBuf, state.cells[state.lastCell].temp);
        double QLbuf = state.cells[state.lastCell - 1].MliqiniRBuf / rhomistBuf;

        state.bufferedCompletionMassSource = -sinal * QLbuf * (betloc)*state.cells[state.lastCell - 1].rcC;
        state.bufferedLiquidMassSource = -sinal * (state.cells[state.lastCell - 1].MliqiniRBuf + state.bufferedCompletionMassSource);
        state.bufferedGasMassSource = -sinal * (state.cells[state.lastCell - 1].MRBuf - state.cells[state.lastCell - 1].MliqiniRBuf);
    } else {
        if (surfaceChokeIsOpen(state) &&
            ((((*state.globals).lixo5 - 2 * state.input.dtmax > state.movingMeanTemperature || state.input.perm == 2) && fabs(state.movingMeanPressure - state.gasSurfacePressure) / state.movingMeanPressure < 0.05 && fabs(state.movingMeanFlux) < 0.5) || (((*state.globals).lixo5 - 2 * state.input.dtmax > state.movingMeanTemperature || state.input.perm == 2) && fabs(state.movingMeanPressure - state.gasSurfacePressure) < (0.05 * state.gasSurfacePressure) && fabs(state.movingMeanFlux) < 5. && delp < 0.01 * state.gasSurfacePressure) || (((*state.globals).lixo5 - 2 * state.input.dtmax > state.movingMeanTemperature || state.input.perm == 2) && (state.gasSurfacePressure - state.movingMeanPressure) > 0.01 && state.input.chkv == 0) || (fabs(delp) < 0.1 && (state.cells[state.lastCell].presBuf - state.gasSurfacePressure) / state.cells[state.lastCell].presBuf < 0.05 && state.input.chkv == 0))) {
        } else {
            abertoini = state.open;
        }
        if (((tit > -0.01 && state.cells[state.lastCell].alf > -0.01) || surfaceChokeIsShut(state)) &&
            state.surfaceChoke.AreaGarg < 0.6 * state.cells[state.lastCell - 1].duto.area && (state.open == 0 && (state.openTime == 0 || state.openTime > 60))) {

            masliq = sinal * maxSup * (1. - tit);
            masgas = sinal * maxSup * tit;
            state.bufferedLiquidMassSource = -masliq * (1 - betSup) * rholp / rholmix;
            state.bufferedCompletionMassSource = -masliq * betSup * rholc / rholmix;
            state.bufferedGasMassSource = -masgas;
        }
        if (state.surfaceChokeMassFlag == 1 && (*state.globals).chaverede == 1) {
            state.cells[state.lastCell].fontemassCR = state.bufferedCompletionMassSource;
            state.cells[state.lastCell].fontemassLR = state.bufferedLiquidMassSource;
            state.cells[state.lastCell].fontemassGR = state.bufferedGasMassSource;

            state.cells[state.lastCell].DmasschokeG = 0.;
            state.cells[state.lastCell].DmasschokeL = 0.;
            state.cells[state.lastCell].DmasschokeC = 0.;
        } else if (state.surfaceChokeMassFlag == 0 && (*state.globals).chaverede == 1) {
            if (state.input.chkv == 0 || state.cells[state.lastCell - 1].MRBuf > 0.)
                sinal = 1.;
            else
                sinal = 0.;
            double betloc;
            if ((state.cells[state.lastCell - 1].MRBuf - state.cells[state.lastCell - 1].MliqiniRBuf) * 0.0 + 1.0 * state.cells[state.lastCell - 1].MliqiniRBuf > 0.)
                betloc = state.cells[state.lastCell - 1].bet; // testeBeta//duvidabeta
            else
                betloc = betRev;
            double rhomistBuf = betloc *
                                    state.cells[state.lastCell - 1].fluicol.MasEspFlu(state.cells[state.lastCell].presBuf, state.cells[state.lastCell].temp) +
                                (1. - betloc) * state.cells[state.lastCell - 1].flui.MasEspLiq(state.cells[state.lastCell].presBuf, state.cells[state.lastCell].temp);
            double QLbuf = state.cells[state.lastCell - 1].MliqiniRBuf / rhomistBuf;

            state.bufferedCompletionMassSource = -sinal * QLbuf * (betloc)*state.cells[state.lastCell - 1].rcC;
            state.bufferedLiquidMassSource = -sinal * (state.cells[state.lastCell - 1].MliqiniRBuf + state.bufferedCompletionMassSource);
            state.bufferedGasMassSource = -sinal * (state.cells[state.lastCell - 1].MRBuf - state.cells[state.lastCell - 1].MliqiniRBuf);
        }
    }
}

void computeImplicitTimeStep(const TransientStepState &state) {
    int multChoke = 1.;
    if (state.cells[state.lastCell - 1].alf < 0.9)
        multChoke = 1;
    double mgas = state.cells[state.lastCell].MC - state.cells[state.lastCell].Mliqini;
    double mgas0 = state.cells[state.lastCell].MCini - state.cells[state.lastCell].Mliqini0;
    if (state.input.RelaxaDTChoke == 0 &&
        (((*state.globals).lixo5 > 1e-15 && state.surfaceChokeMassFlag == 1 && state.surfaceChoke.AreaGarg < 0.6 * state.cells[state.lastCell - 1].duto.area &&
          state.surfaceChoke.AreaGarg > 1e-3 * state.cells[state.lastCell - 1].duto.area &&
          ((state.cells[state.lastCell].Mliqini >= 0 && state.cells[state.lastCell].Mliqini0 < 0 && mgas > 0) || (state.cells[state.lastCell].Mliqini >= 0 && mgas < 0 && mgas0 > 0))) ||
         state.slugCount < multChoke * 200)) {

        if (state.slugCount > multChoke * 200) {
            state.slugCount = 0;
        }
        if (((state.cells[state.lastCell].Mliqini >= 0 && state.cells[state.lastCell].Mliqini0 < 0 && mgas > 0) ||
             (state.cells[state.lastCell].Mliqini >= 0 && mgas < 0 && mgas0 > 0)) &&
            state.slugCount > multChoke * 100)
            state.slugCount = multChoke * 100;
        double progres = 1.;
        if (state.slugCount > multChoke * 100) {
            progres = state.slugCount - multChoke * 100;
        }
        state.timeStep *= (progres / (multChoke * 100.));
        state.slugCount++;
    }
    if (state.surfaceChoke.AreaGarg >= 0.6 * state.cells[state.lastCell - 1].duto.area)
        state.open = 1;
    if (fabs(state.finalPressure - state.gasSurfacePressure) / state.finalPressure < 0.05 && state.open == 0 && state.timeStep > 1. && state.surfaceChoke.AreaGarg < 0.6 * state.cells[state.lastCell - 1].duto.area &&
        state.surfaceChoke.AreaGarg > 1e-15 * state.cells[state.lastCell - 1].duto.area)
        state.timeStep =
            1.;
    for (int i = 0; i < state.input.eventoabre; i++) {
        if ((*state.globals).lixo5 > state.input.Tevento[i] - state.input.dtmax && (*state.globals).lixo5 < state.input.Tevento[i] + 30) {
            if (state.timeStep > state.smallestCellLength / 100.) {
                state.timeStep = state.smallestCellLength / 100.;
            }
        }
    }
    state.multiplier = 0.8;
    if ((((*state.globals).lixo5 - 2 * state.input.dtmax) > state.movingMeanTemperature && (fabs(state.movingMeanPressure - state.finalPressure) / state.finalPressure > 0.4)) || (state.finalPressure < state.gasSurfacePressure && state.surfaceChokeMassFlag == 1 && state.surfaceChoke.AreaGarg > (1e-3) * state.cells[state.lastCell - 1].duto.area) || (state.cells[state.lastCell].alf <= 0.1 && state.surfaceChoke.AreaGarg < 0.6 * state.cells[state.lastCell - 1].duto.area && state.surfaceChoke.AreaGarg > 1e-3 * state.cells[state.lastCell - 1].duto.area && state.open == 0)) {
        double denominador = 10.;
        double progres = 0.95;

        if ((state.finalPressure < state.gasSurfacePressure && state.surfaceChokeMassFlag == 1)) {
            denominador = 100.;
            progres = 0.5;
            state.multiplier = pow(progres, 10. * fabs(state.gasSurfacePressure - state.finalPressure) / state.finalPressure) * state.multiplier;
        } else if ((state.cells[state.lastCell].alf <= 0.1 && state.surfaceChoke.AreaGarg < 0.6 * state.cells[state.lastCell - 1].duto.area)) {
            denominador = 1000.;
            progres = 0.5;
            state.multiplier = (state.cells[state.lastCell].alf + 1e-5) * state.multiplier;
        } else {
            state.multiplier = pow(progres, 10 * fabs(state.movingMeanPressure - state.finalPressure) / state.finalPressure) * state.multiplier;
        }
        if (state.input.dtmax * state.multiplier < state.smallestCellLength / denominador)
            state.multiplier = (state.smallestCellLength / denominador) / state.input.dtmax;
        if (state.multiplier > 0.8)
            state.multiplier = 0.8;
    } else {
        if (state.cells[state.lastCell].alf >= 0.11 || state.surfaceChoke.AreaGarg >= 0.6 * state.cells[state.lastCell - 1].duto.area ||
            state.surfaceChoke.AreaGarg <= 1e-3 * state.cells[state.lastCell - 1].duto.area || state.open == 1) {
            state.multiplier = state.multiplier / 0.95;
            if (state.multiplier > 0.8)
                state.multiplier = 0.8;
        }
    }

    for (int i = 0; i <= state.lastCell; i++) {
        double jmix = 0.;
        double A1 = state.cells[i].duto.area;
        double dtaux;
        // celula[i].Mliqini
        double alfteste = state.cells[i].alf;
        if (i > 0 && (state.cells[i].MC - state.cells[i].Mliqini) > 0)
            alfteste = state.cells[i - 1].alf;
        if (alfteste > 1e-9)
            jmix += fabs(
                (state.cells[i].MC - state.cells[i].Mliqini) / (A1 * state.cells[i].flui.MasEspGas(state.cells[i].pres, state.cells[i].temp) * (0 + 1 * alfteste)));
        if (fabs(state.cells[i].VTemper) > jmix)
            jmix = fabs(state.cells[i].VTemper);

        double jmixL = 0.;
        if (alfteste < 1 - 1e-5)
            jmixL = fabs(
                state.cells[i].Mliqini / (A1 * ((1. - state.cells[i].bet) * state.cells[i].rpC + state.cells[i].bet * state.cells[i].rcC) * (1. - 0 * alfteste)));
        if (fabs(jmixL) > jmix)
            jmix = fabs(jmixL);
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

    if (state.timeStepRestricted == 1 && state.input.desligaPenalizaDT == 0) {
        if (state.timeChanged < 2)
            state.timeStep /= 10.;
        state.timeChanged++;
        if (state.timeChanged > 10)
            state.timeChanged = 0;
    }

    state.timeStep = state.multiplier * state.timeStep;
    if (state.timeStep > state.input.dtmax)
        state.timeStep = state.input.dtmax;

    if (state.input.evento.size() > state.eventIndex) {
        if ((*state.globals).lixo5 < state.input.evento[state.eventIndex] && ((*state.globals).lixo5 + state.timeStep) > (state.input.evento[state.eventIndex] + 0.1)) {
            state.timeStep = state.input.evento[state.eventIndex] - (*state.globals).lixo5;
            state.eventIndex++;
        } else if ((*state.globals).lixo5 < state.input.evento[state.eventIndex] && ((*state.globals).lixo5 + state.timeStep) >= (state.input.evento[state.eventIndex]))
            state.eventIndex++;
    }

    state.masterState = 1;
    for (int i = 0; i < state.masterCloseCount; i++) {
        if ((*state.globals).lixo5 > state.masterCloseSchedule[i]) {
            for (int j = 0; j < state.masterOpenCount; j++) {
                if ((*state.globals).lixo5 < state.masterOpenSchedule[j]) {
                    state.masterState = 0;
                    state.masterCounter++;
                    if (state.masterCounter > 19)
                        state.masterCounter = 20;
                    if ((*state.globals).lixo5 >= state.masterOpenSchedule[j] - state.timeStep)
                        state.masterCounter = 0;
                    break;
                }
            }
        }
    }

    for (int i = 0; i <= state.lastCell; i++) {
        state.cells[i].dt = state.timeStep;
        state.cells[i].dt1 = state.timeStep;
        state.cells[i].dt2 = state.timeStep;
        state.cells[i].dtPig = state.timeStep;
    }

    for (int i = 0; i <= state.lastCell; i++) {
        if (state.cells[i].acsr.tipo == 15) {
            state.cells[i].acsr.radialPoro.dt = state.timeStep;
        } else if (state.cells[i].acsr.tipo == 16) {
            state.cells[i].acsr.poroso2D.dt = state.timeStep;
        }
    }
    // dtInter=dt;//alteracao2
}

void computeExplicitTimeStep(const TransientStepState &state) {

    for (int i = 0; i <= state.lastCell; i++) {
        double velAux = state.cells[i].termAdSomVel();
        double som = state.cells[i].somVel();
        double velpropag1 = velAux + som;
        double velpropag2 = fabs(velAux - som);
        double velMax = velpropag1;
        if (velpropag2 > velpropag1)
            velMax = velpropag2;
        double dtaux;
        double alfteste = state.cells[i].alf;
        dtaux = state.cells[i].dx / velMax;
        if (dtaux < state.timeStep)
            state.timeStep = dtaux;
    }
    for (int i = 0; i <= state.lastCell; i++) {
        state.cells[i].dt = state.timeStep;
        state.cells[i].dt1 = state.timeStep;
        state.cells[i].dt2 = state.timeStep;
        state.cells[i].dtPig = state.timeStep;
    }
}

void computeTimeStep(const TransientStepState &state, int vexpli) {

    if ((*state.globals).lixo5 < (*state.globals).localtiny && (*state.globals).chaverede == 0) {
        state.input.imprimeProfile(state.cells, state.productionFreeTerms, (*state.globals).lixo5, state.branchIndex);
        if (state.input.lingas > 0 && state.input.nvalvgas > 0)
            state.input.imprimeProfileG(state.gasCells, state.gasFreeTerms, (*state.globals).lixo5, state.branchIndex);
        state.input.imprimeProfileTrans(state.cells, state.productionCrossSectionCount, (*state.globals).lixo5, state.branchIndex);
    }
    state.timeStep = state.input.dtmax;

    int parada = 0;
    for (int i = 1; i < state.lastCell; i++) {
        if (state.cells[i].acsr.tipo == 5 && state.cells[i].acsr.chk.AreaGarg <= 1e-15 * state.cells[i].acsr.chk.AreaTub)
            parada = 1;
        else if (state.surfaceChoke.AreaGarg <= 1.e-15 * state.surfaceChoke.AreaTub)
            parada = 1;
    }
    if (parada == 1) {
        Vcr<int> oscila(state.lastCell, 0);
        for (int i = 0; i <= state.lastCell; i++) {
            if (i > 0 && i < state.lastCell) {
                if (fabs(state.cells[i - 1].alf - state.cells[i + 1].alf) < fabs(state.cells[i].alf - state.cells[i - 1].alf)) {
                    double area = state.cells[i].duto.area;
                    double vLiqTest = fabs(state.cells[i].QL / (area));
                    double vGasTest = fabs(state.cells[i].QG / (area));
                    if ((vLiqTest + vGasTest) > 0.1)
                        oscila[i] = 1;
                }
            }
        }
        int alarmOscila = 0;
        for (int i = 1; i <= state.lastCell - 4; i++) {
            int kontaOsc = 0;
            int multOsc = 0;
            while (kontaOsc < 3) {
                multOsc += oscila[i + kontaOsc];
                kontaOsc++;
            }
            if (multOsc == 3)
                alarmOscila = 1;
        }

        if (alarmOscila == 1 && state.input.desligaPenalizaDT == 0)
            state.timeStep /= 10.;
    }

    for (int i = 0; i <= state.lastCell; i++) {
        if (state.cells[i].acsr.tipo == 15) {
            state.cells[i].acsr.radialPoro.defineDT(0);
            if (state.cells[i].acsr.radialPoro.dt < state.timeStep)
                state.timeStep = state.cells[i].acsr.radialPoro.dt;
        }
        if (state.cells[i].acsr.tipo == 16) {
            state.cells[i].acsr.poroso2D.defineDT(0);
            if (state.cells[i].acsr.poroso2D.dt < state.timeStep)
                state.timeStep = state.cells[i].acsr.poroso2D.dt;
        }
    }

    if (vexpli == 1) {
        computeExplicitTimeStep(state);
    } else {

        computeImplicitTimeStep(state);
    }
}

void valveOpeningLow(const TransientStepState &state) {
    int celpos = state.input.master1.posic;
    state.masterRatio0[0] = state.cells[celpos].acsr.chk.AreaGarg / state.cells[celpos].duto.area;
    for (int i = 1; i <= state.input.nvalv; i++) {
        celpos = state.input.valv[i - 1].posicP;
        state.masterRatio0[i] = state.cells[celpos].acsr.chk.AreaGarg / state.cells[celpos].duto.area;
    }
}

void valveOpeningHigh(const TransientStepState &state) {
    int celpos = state.input.master1.posic;
    state.masterRatio1[0] = state.cells[celpos].acsr.chk.AreaGarg / state.cells[celpos].duto.area;
    for (int i = 1; i <= state.input.nvalv; i++) {
        celpos = state.input.valv[i - 1].posicP;
        state.masterRatio1[i] = state.cells[celpos].acsr.chk.AreaGarg / state.cells[celpos].duto.area;
    }
}

void dampMaximumTimeStep(const TransientStepState &state) {

    state.cflTimeSteps.push_back(state.auxiliaryCflTimeStep);
    state.simulationTimeSteps.push_back(state.finalAuxiliaryTimeStep);
    state.totalCflTimeStep += state.auxiliaryCflTimeStep;
    state.totalSimulationTimeStep += state.finalAuxiliaryTimeStep;
    state.timeStepRestrictionCount++;
    if (state.stepIndex > 10) {
        state.totalCflTimeStep -= state.cflTimeSteps.front();
        state.cflTimeSteps.erase(state.cflTimeSteps.begin());
        state.totalSimulationTimeStep -= state.simulationTimeSteps.front();
        state.simulationTimeSteps.erase(state.simulationTimeSteps.begin());
        if (state.timeStepRestrictionCount > 10 && state.timeStepRestricted == 1) {
            state.timeStepRestricted = 0;
            state.timeStepRestrictionCount = 0;
        }
    }
    state.meanCflTimeStep = state.totalCflTimeStep / 10;
    state.meanSimulationTimeStep = state.totalSimulationTimeStep / 10.;
    if ((state.meanSimulationTimeStep < state.meanCflTimeStep / 2 && state.timeStepRestrictionCount > 10)) {
        state.timeStepRestricted = 1;
        state.timeStepRestrictionCount = 0;
    }
}

void evaluatePressureRateOfChange(const TransientStepState &state, double razMast, double razMast0, int vexpli) {
    double dpdtRef = 2 * state.input.taxaDespre;
    double modDpDt = 0.;
    double modDTDt = 0.;
    double dpdtMax = 0.;
    double dTdtMax = 0.;
    if (state.fullModel == 1) {
        state.fullModel = 0;
        int i = 0;
        while (i < state.lastCell) {
            int nBloco = 10;
            if (state.lastCell - i < nBloco)
                nBloco = state.lastCell - i;
            for (int j = 0; j < nBloco; j++) {
                modDpDt += fabs(state.cells[i + j].pres - state.cells[i + j].presini) / (state.timeStep);
                modDTDt += fabs(state.cells[i + j].temp - state.cells[i + j].tempini) / (state.timeStep);
            }
            i += nBloco;
            modDpDt /= nBloco;
            modDTDt /= nBloco;
            if (modDpDt > dpdtMax)
                dpdtMax = modDpDt;
            if (modDTDt > dTdtMax)
                dTdtMax = modDTDt;
            if (modDpDt > dpdtRef)
                state.fullModel = 1;
            modDpDt = 0.;
            modDTDt = 0.;
        }
        state.maximumPressureRates.push_back(dpdtMax);
        state.maximumTimeStepRates.push_back(dTdtMax);
        if (state.maximumPressureRates.size() > 10)
            state.maximumPressureRates.erase(state.maximumPressureRates.begin());
        if (state.maximumTimeStepRates.size() > 10)
            state.maximumTimeStepRates.erase(state.maximumTimeStepRates.begin());
        int nvec = state.maximumPressureRates.size();
        int nvecT = state.maximumTimeStepRates.size();
        state.meanMaximumPressureChange = 0.;
        state.meanMaximumTimeStep = 0.;
        for (int i = 0; i < nvec; i++)
            state.meanMaximumPressureChange += state.maximumPressureRates[i];
        for (int i = 0; i < nvecT; i++)
            state.meanMaximumTimeStep += state.maximumTimeStepRates[i];
        state.meanMaximumPressureChange /= nvec;
        state.meanMaximumTimeStep /= nvecT;
        if (state.meanMaximumPressureChange > dpdtRef)
            state.fullModel = 1;
        if (state.surfaceChoke.AreaGarg / state.cells[state.lastCell - 1].duto.area < 1e-3 &&
            (state.meanMaximumPressureChange > state.input.taxaDespre / 10. || state.meanMaximumTimeStep > 0.001))
            state.fullModel = 1;
        int linAberta = 1; // caso varias valvulas
        for (int i = 0; i <= state.input.nvalv; i++)
            if (state.masterRatio1[i] <= 1e-3)
                linAberta = 0; // caso varias valvulas
        if ((linAberta == 1 && state.surfaceChoke.AreaGarg / state.cells[state.lastCell - 1].duto.area > 1e-3))
            state.fullModel = 1; // caso varias valvulas
    }
    for (int i = 0; i <= state.input.nvalv; i++)
        if (state.masterRatio1[i] != state.masterRatio0[i])
            state.fullModel = 0; // caso varias valvulas
    if ((state.meanMaximumPressureChange > 10 || state.meanMaximumTimeStep > 1) && state.fullModel == 1 && vexpli == 0) {
        state.fullModel = 0;
    }
}

void restrictTimeStepByValve(const TransientStepState &state) {

    if (state.input.ConContEntrada == 0 && (*state.globals).chaveRedeParalela == 0) {
        int celpos;
        double dtaux = state.cells[0].dt;
        if (state.restart == -1) {
            for (int i = 1; i <= state.lastCell; i++)
                if (dtaux > state.cells[i].dt)
                    dtaux = state.cells[i].dt;
        }
        double dtvec[state.input.nvalv + 1];
        for (int i = 0; i <= state.input.nvalv; i++)
            dtvec[i] = state.cells[0].dt;
        for (int i = 0; i <= state.input.nvalv; i++) {
            if (i == 0)
                celpos = state.input.master1.posic;
            else
                celpos = state.input.valv[i - 1].posicP;
            if (state.masterRatio1[i] < state.masterRatio0[i] && (state.masterRatio1[i] <= 1.1 * state.input.master1.razareaativ && state.masterRatio1[i] >= 1e-3 * state.input.master1.razareaativ)) {
                if (state.cells[celpos].alf < 0.05) {
                    if (state.cells[celpos].alf < 0.01)
                        state.desperationMoment += 1.5;
                    else
                        state.desperationMoment = 1;
                    if (state.desperationMoment < 1.)
                        state.desperationMoment = 1.;
                    state.cells[celpos].fontemassGL += state.desperationMoment * 10000 * state.cells[celpos].flui.Deng * 1.225 / 86400;
                    state.cells[celpos - 1].fontemassGR = state.cells[celpos].fontemassGL;
                }
                if (state.cells[celpos].alf < 0.5)
                    state.masterCriticalRatio[i] = 0.01;
            } else
                state.desperationMoment = 0.;
            if (state.masterRatio1[i] < state.masterRatio0[i] && (state.masterRatio1[i] <= 1.1 * state.input.master1.razareaativ && state.masterRatio1[i] >= state.masterCriticalRatio[i] * state.input.master1.razareaativ)) {
                double raz = 20.;
                if (state.cells[celpos].alf < 0.5)
                    raz = 40.;
                dtvec[i] = dtaux;
                if (dtvec[i] > 1)
                    dtvec[i] = 1.;
                if (state.cells[celpos].alf < 0.5)
                    dtvec[i] = 0.1 + 0.9 * (state.cells[celpos].alf) / 0.5;
                dtvec[i] /= raz;
                state.restart = -1;
            }
        }
        dtaux = dtvec[0];
        for (int i = 1; i <= state.input.nvalv; i++)
            if (dtvec[i] < dtaux)
                dtaux = dtvec[i];
        state.cells[0].dt = dtaux;
    }
}

void restartFractionEvolutionInitial(const TransientStepState &state) {
    for (int i = 0; i <= state.lastCell; i++) {
        if (state.cells[i].dt < state.timeStep)
            state.timeStep = state.cells[i].dt;
        if (state.cells[i].dt1 < state.timeStep)
            state.timeStep = state.cells[i].dt1;
        if (state.cells[i].dt2 < state.timeStep)
            state.timeStep = state.cells[i].dt2;
        if (state.cells[i].dtPig < state.timeStep)
            state.timeStep = state.cells[i].dtPig;
    }
}

void restartFractionEvolutionSub(const TransientStepState &state) {
    for (int i = 0; i <= state.lastCell; i++) {
        if (state.cells[i].pres > -10.) {
            state.cells[i].alf = state.cells[i].alfini;
            state.cells[i].alfPigE = state.cells[i].alfPigEini;
            state.cells[i].alfPigD = state.cells[i].alfPigDini;
        }
    }
}

void restartFractionEvolution(const TransientStepState &state) {
    for (int i = 0; i <= state.lastCell; i++) {
        state.cells[i].dt = state.timeStep;
        state.cells[i].dt1 = state.timeStep;
        state.cells[i].dt2 = state.timeStep;
        state.cells[i].dtPig = state.timeStep;
    }
    if (state.input.lingas > 0) {
        for (int i = 0; i <= state.gasCellCount; i++)
            state.gasCells[i].FeiticoDoTempo();
        if (state.input.descarga == 1) {
            state.interfaceCell = state.initialInterfaceCell;
            state.interfaceTimeStep = state.initialInterfaceTimeStep;
            state.interfaceVelocity = state.initialInterfaceVelocity;
        }
        state.updaters.advanceGasSubStep();
    }
    restartFractionEvolutionSub(state);
    for (int i = 0; i <= state.lastCell; i++) {
        state.cells[i].bet = state.cells[i].betini;
        state.cells[i].razPig = state.cells[i].razPigini;
        state.cells[i].betPigE = state.cells[i].betPigEini;
        state.cells[i].betPigD = state.cells[i].betPigDini;
    }
}

void evolveFractions(const TransientStepState &state, double alfrev, double betrev, int ciclo) {

#pragma omp parallel for num_threads((*state.globals).ntrd)
    for (int i = 0; i <= state.lastCell; i++) {

        if (i < state.lastCell) {
            state.cells[i].avancalf(state.restart, state.lastCell);
        }
        if (i == state.lastCell) {

            if (((*state.globals).chaverede == 0 || state.endNode == 1 || (*state.globals).chaveRedeParalela == 1)) {
                if (state.surfaceChokeMassFlag == 0 || state.cells[state.lastCell].Mliqini > 0)
                    state.cells[i].alf = state.cells[i - 1].alf;
                else
                    state.cells[i].avancalf(state.restart, state.lastCell);
            } else {
                if (state.cells[state.lastCell].Mliqini > 0)
                    state.cells[i].alf = state.cells[i - 1].alf;
                else if (state.surfaceChokeMassFlag == 0 && state.input.chkv == 0)
                    state.cells[i].alf = alfrev;
                else
                    state.cells[i].avancalf(state.restart, state.lastCell);
            }
        }
    }
#pragma omp parallel for num_threads((*state.globals).ntrd)
    for (int i = 0; i <= state.lastCell; i++) {
        if (i < state.lastCell) {
            state.cells[i].avancbet(state.restart, state.lastCell);
        }
        if (i == state.lastCell) {
            if (((*state.globals).chaverede == 0 || state.endNode == 1 || (*state.globals).chaveRedeParalela == 1)) {
                if (state.surfaceChokeMassFlag == 0 || state.cells[state.lastCell].Mliqini > 0)
                    state.cells[i].bet = state.cells[i - 1].bet;
                else
                    state.cells[i].avancbet(state.restart, state.lastCell);
            } else {
                if (state.surfaceChokeMassFlag == 1)
                    state.cells[i].avancbet(state.restart, state.lastCell);
                else if (state.input.chkv == 1)
                    state.cells[i].bet = state.cells[i - 1].bet;
            }
        }
    }
#pragma omp parallel for num_threads((*state.globals).ntrd)
    for (int i = 0; i <= state.lastCell; i++) {
        if (i < state.lastCell && state.cells[i].estadoPig == 1) {
            state.cells[i].avancPig(state.restart);
            state.cells[i].avancalfPig();
            state.cells[i].avancbetPig();
        } else {
            state.cells[i].alfPigE = state.cells[i].alf;
            state.cells[i].betPigE = state.cells[i].bet;
            state.cells[i].alfPigD = state.cells[i].alf;
            state.cells[i].betPigD = state.cells[i].bet;
        }
    }
    for (int i = 0; i <= state.lastCell; i++) {
        if (state.cells[i].correrGlobHol == 1) {
            if (state.cells[i].reiniciaAlf < 0 || state.cells[i].reiniciaBet < 0 || state.cells[i].reiniciaPig < 0)
                state.restart = -1;
            state.cells[i].reiniciaAlf = 0;
            state.cells[i].reiniciaBet = 0;
            state.cells[i].reiniciaPig = 0;
        }
    }
}

void updatePig(const TransientStepState &state) {
    for (int i = 1; i <= state.lastCell; i++) {
        state.cells[i].velPigini = state.cells[i].velPig;
        state.cells[i].estadoPigini = state.cells[i].estadoPig;
        state.cells[i].indpigini = state.cells[i].indpig;
    }
    for (int i = 1; i <= state.lastCell; i++) {
        if (i < state.lastCell && state.cells[i].estadoPig == 1) {
            if (state.cells[i].velPig >= 0) {
                if (state.cells[i].razPig >= 1. - (*state.globals).localtiny) {
                    state.cells[i].estadoPig = 0;
                    state.cells[i + 1].estadoPig = 1;

                    state.cells[i].razPig = 0.;
                    state.cells[i + 1].razPig = 0.;
                    state.cells[i + 1].indpig = state.cells[i].indpig;
                    state.cells[i].indpig = -1;
                    state.cells[i + 1].velPig = state.cells[i].velPig;
                    state.cells[i + 1].alfPigE = state.cells[i].alf;
                    state.cells[i + 1].betPigE = state.cells[i].bet;
                }
            } else {
                if (state.cells[i].razPig <= (*state.globals).localtiny) {
                    state.cells[i].estadoPig = 0;
                    state.cells[i - 1].estadoPig = 1;
                    state.cells[i].razPig = 0.;
                    state.cells[i - 1].razPig = 1.;
                    state.cells[i - 1].indpig = state.cells[i].indpig;
                    state.cells[i].indpig = -1;
                    state.cells[i - 1].velPig = state.cells[i].velPig;
                    state.cells[i - 1].alfPigD = state.cells[i].alf;
                    state.cells[i - 1].betPigD = state.cells[i].bet;
                }
            }
        }
    }
    if (state.input.ConContEntrada == 0) {
        state.cells[0].betI = state.cells[0].bet;
        state.cells[0].betLI = state.cells[0].bet;
    } else {
        state.cells[0].betI = state.inletCompletionFraction;
        if ((state.cells[0].MC - state.cells[0].Mliqini) * 0 + state.cells[0].Mliqini < 0.)
            state.cells[0].betI = state.cells[0].betPigE; // testeBeta
        state.cells[0].betLI = state.cells[0].betI;
    }
    for (int i = 1; i <= state.lastCell; i++) {
        double betLI;
        double betI;
        double betRI;
        state.cells[i].betI = state.cells[i].betPigE;
        if (state.cells[i].QL > 0.)
            state.cells[i].betI = state.cells[i - 1].betPigD; // testeBeta
        state.cells[i - 1].betRI = state.cells[i].betI;
        if (i < state.lastCell)
            state.cells[i + 1].betLI = state.cells[i].betI;
    }
}

void solvePressureVolumeCoupling(const TransientStepState &state, int vexpli, int ciclo) {
#pragma omp parallel for num_threads((*state.globals).ntrd)
    for (int i = 0; i <= state.lastCell; i++) {
        state.cells[i].GeraLocal(state.finalPressure, state.surfaceChokeMassFlag, state.lastCell, state.input.master1.razareaativ, state.inletPressure, state.inletTemperature, state.inletQuality, state.inletCompletionFraction, ciclo,
                            state.fullModel, state.endNode, state.input.corrigeContSep, state.surfaceChoke.AreaGarg, vexpli);
        for (int j = 0; j < 6; j++) {
            state.productionMatrix[2 * i][j - 3] = state.cells[i].local[0][j];
            state.productionMatrix[2 * i + 1][j - 3] = state.cells[i].local[1][j];
            state.productionSolution[2 * i] = state.cells[i].TL[0];
            state.productionSolution[2 * i + 1] = state.cells[i].TL[1];
        }
    }

    state.productionMatrix.GaussElimPP(state.productionSolution);
}

void refreshFluidMiniTable(const TransientStepState &state) {
    //if(arq.miniTabAtraso>0)
    	state.updaters.generateFluidMiniTable();
    double betIV;
    double rsV;
    double boV;
    double baV;
    double bswV;
    double rhoOVol;
    double rhoWVol;
    double titVol;
    double rhoGVol;
    double ZGVol;
    double DZDPGVol;
    double DZDTGVol;
    for (int i = 0; i < state.lastCell; i++) {

        double pres = state.cells[i].pres;
        double temp = state.cells[i].temp;

        betIV = state.cells[i].bet;
        rsV = state.cells[i].flui.RS(state.cells[i].pres, state.cells[i].temp);
        boV = state.cells[i].flui.BOFunc(state.cells[i].pres, state.cells[i].temp, rsV);
        baV = state.cells[i].flui.BAFunc(state.cells[i].pres, state.cells[i].temp);
        bswV = state.cells[i].flui.BSW * baV / (boV + baV * state.cells[i].flui.BSW - state.cells[i].flui.BSW * boV);
        rhoOVol = state.cells[i].flui.MasEspoleo(state.cells[i].pres, state.cells[i].temp);
        rhoWVol = state.cells[i].flui.MasEspAgua(state.cells[i].pres, state.cells[i].temp);
        titVol = (1 - bswV) * rhoOVol / ((1 - bswV) * rhoOVol + bswV * rhoWVol);

        rhoGVol = state.cells[i].flui.MasEspGas(state.cells[i].pres, state.cells[i].temp);
        ZGVol = state.cells[i].flui.Zdran(state.cells[i].pres, state.cells[i].temp);
        DZDPGVol = state.cells[i].flui.FracMassHidra(state.cells[i].pres, state.cells[i].temp);
        DZDTGVol = state.cells[i].flui.PB(state.cells[i].pres, state.cells[i].temp);
        state.cells[i].nMol = (state.cells[i].flui.MasEspLiq(pres, temp) * (1. - state.cells[i].alf) * (1. - betIV) * titVol +
                          state.cells[i].rgC * state.cells[i].alf) *
                         state.cells[i].duto.area * state.cells[i].dx / state.cells[i].flui.Pmol;
        state.cells[i].nMolIni = state.cells[i].nMol;
    }
}

void refreshInletCondition(const TransientStepState &state) {
    if (state.input.ConContEntrada == 1) {
        if (state.input.tipoFluido == 0 && state.input.flashCompleto == 2) {
            double rgST = state.cells[0].flui.Deng * 1.225;
            double roST = 141.5 * 1000. / (131.5 + state.cells[0].flui.API);
            double rg = state.cells[0].flui.MasEspGas(state.inletPressure, state.inletTemperature);
            double rl = state.cells[0].flui.MasEspLiq(state.inletPressure, state.inletTemperature);
            double titH = state.cells[0].flui.FracMassHidra(state.inletPressure, state.inletTemperature);
            double rcST = state.cells[0].fluicol.MasEspFlu(1.01, 20.);
            double rc = state.cells[0].fluicol.MasEspFlu(state.inletPressure, state.inletTemperature);
            double rlMix = state.inletCompletionFraction * rc + (1. - state.inletCompletionFraction) * rl;
            double val1 = ((1. - state.inletCompletionFraction) * rl * titH / (1. - titH));
            state.inletQuality = val1 / (rlMix + val1);
        } else if (state.input.tipoFluido == 1) {
            double rgST = state.cells[0].flui.Deng * 1.225;
            double roST = 141.5 * 1000. / (131.5 + state.cells[0].flui.API);
            double rg = state.cells[0].flui.MasEspGas(state.inletPressure, state.inletTemperature);
            double rl = state.cells[0].flui.MasEspoleo(state.inletPressure, state.inletTemperature);
            double tit = state.cells[0].flui.FracMass(state.inletPressure, state.inletTemperature);
            double rcST = state.cells[0].fluicol.MasEspFlu(1.01, 20.);
            double rc = state.cells[0].fluicol.MasEspFlu(state.inletPressure, state.inletTemperature);
            double val1 = (rcST / rc) * (rg / rgST) * state.input.CCPres.bet[0] / tit;
            double val2 = (rg / rl) * (1 - tit) / tit;
            double titT = rg / (((1. - tit) / tit) * (rg / rl) + rg + val1);
            state.inletQuality = titT;
            state.inletCompletionFraction = val1 / (val2 + val1);
        }
    }
}

}  // namespace sisprod::transient
