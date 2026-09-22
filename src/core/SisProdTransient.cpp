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

}  // namespace sisprod::transient
