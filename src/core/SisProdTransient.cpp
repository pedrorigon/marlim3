#include "SisProdTransient.h"

#include "Acidentes2.h"
#include "Leitura.h"
#include "SisProdConstants.h"
#include "Matriz.h"
#include "Vetor.h"
#include "solver3DPoisson.h"
#include "celula3.h"
#include "celulaGas.h"
#include "variaveisGlobais1D.h"

#include <math.h>

// For the TimeStepPolicy concept.
#include <concepts>

// The run's start date, which the progress report prints. Globals defined
// elsewhere and declared extern in SisProd.h -- a header this module does not
// include on purpose, so the four declarations it needs are repeated here, and
// only those four.
extern int diaIni;
extern int horaIni;
extern int minutoIni;
extern int segundoIni;

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
    double diameter = state.cells[i].duto.a;
    double area = 0.25 * M_PI * diameter * diameter;
    double perimeter = state.cells[i].duto.peri;
    double alfmed = state.cells[i].alf;
    double rhog = state.cells[i].flui.MasEspGas(state.cells[i].pres, state.cells[i].temp);
    double rhol = (1 - state.cells[i].bet) * state.cells[i].flui.MasEspLiq(state.cells[i].pres, state.cells[i].temp) + state.cells[i].bet * state.cells[i].fluicol.MasEspFlu(state.cells[i].pres, state.cells[i].temp);
    double ugsmed = (state.cells[i].MC - state.cells[i].Mliqini) / (area * rhog);
    double ulsmed = state.cells[i].Mliqini / (area * rhol);
    double j = ugsmed + ulsmed;

    double rhomix = alfmed * rhog + (1 - alfmed) * rhol;
    double viscmix = alfmed * state.cells[i].flui.ViscGas(state.cells[i].pres, state.cells[i].temp) + (1 - alfmed) * ((1 - state.cells[i].bet) * state.cells[i].flui.ViscOleo(state.cells[i].pres, state.cells[i].temp) + state.cells[i].bet * state.cells[i].fluicol.VisFlu(state.cells[i].pres, state.cells[i].temp));

    double reynolds;
    if (state.cells[i].duto.revest == 0)
        reynolds = state.cells[i].Rey(state.cells[i].duto.a, j, rhomix, viscmix);
    else {
        double dhid = 4 * area / perimeter;
        reynolds = state.cells[i].Rey(dhid, j, rhomix, viscmix);
    }
    double frictionFactor = state.cells[i].fric(reynolds, state.cells[i].duto.rug / diameter);
    double medpres = 0;
    if (state.input.MedSimpPresFront == 0) {
        if (state.cells[i].presaux <= 10)
            medpres = 1;
        else
            medpres = 0.;
    } else
        medpres = 1;
    double gradfric = (1 - medpres) * 0.5 * frictionFactor * rhomix * (fabs(j) * j) * perimeter * dx / area;
    double gradhidro = (1 - medpres) * kGravity * sin(state.cells[i].duto.teta) * rhomix * dx;
    state.cells[i].presauxini = state.cells[i].presaux;
    state.cells[i].presaux = state.cells[i].pres + (gradfric + gradhidro - state.cells[i - 1].dpB) / kPascalPerKgfPerCm2;
    state.cells[i].dpresaux = 0.5 * (gradfric + gradhidro - state.cells[i - 1].dpB) / kPascalPerKgfPerCm2;
    dx = 0.5 * state.cells[i].dxL;
    diameter = state.cells[i - 1].duto.a;
    area = 0.25 * M_PI * diameter * diameter;
    perimeter = state.cells[i - 1].duto.peri;
    alfmed = state.cells[i - 1].alf;
    rhog = state.cells[i - 1].flui.MasEspGas(state.cells[i - 1].pres, state.cells[i - 1].temp);
    rhol = (1 - state.cells[i - 1].bet) * state.cells[i - 1].flui.MasEspLiq(state.cells[i - 1].pres, state.cells[i - 1].temp) + state.cells[i - 1].bet * state.cells[i - 1].fluicol.MasEspFlu(state.cells[i - 1].pres, state.cells[i - 1].temp);
    ugsmed = (state.cells[i - 1].MC - state.cells[i - 1].Mliqini) / (area * rhog);
    ulsmed = state.cells[i - 1].Mliqini / (area * rhol);
    j = ugsmed + ulsmed;

    rhomix = alfmed * rhog + (1 - alfmed) * rhol;
    viscmix = alfmed * state.cells[i - 1].flui.ViscGas(state.cells[i - 1].pres, state.cells[i - 1].temp) + (1 - alfmed) * ((1 - state.cells[i - 1].bet) * state.cells[i - 1].flui.ViscOleo(state.cells[i - 1].pres, state.cells[i - 1].temp) + state.cells[i - 1].bet * state.cells[i - 1].fluicol.VisFlu(state.cells[i - 1].pres, state.cells[i - 1].temp));

    if (state.cells[i - 1].duto.revest == 0)
        reynolds = state.cells[i - 1].Rey(state.cells[i - 1].duto.a, j, rhomix, viscmix);
    else {
        double dhid = 4 * area / perimeter;
        reynolds = state.cells[i - 1].Rey(dhid, j, rhomix, viscmix);
    }
    frictionFactor = state.cells[i - 1].fric(reynolds, state.cells[i - 1].duto.rug / diameter);

    gradfric = (1 - medpres) * 0.5 * frictionFactor * rhomix * (fabs(j) * j) * perimeter * dx / area;
    gradhidro = (1 - medpres) * kGravity * sin(state.cells[i - 1].duto.teta) * rhomix * dx;

    if (state.cells[i - 1].acsr.tipo != kAccessoryChoke || state.cells[i - 1].acsr.chk.AreaGarg > state.cells[i - 1].acsr.chk.AreaTub * 0.5)
        state.cells[i].presaux = 0.5 * (state.cells[i].presaux) +
                            0.5 * (state.cells[i - 1].pres - (gradfric + gradhidro) / kPascalPerKgfPerCm2);
    state.cells[i].dpresaux -= 0.5 * (gradfric + gradhidro) / kPascalPerKgfPerCm2;
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
        betI = state.cells[i].bet; // beta doubt
    else
        betI = state.cells[i].betL;

    double liquidDensity = flud.MasEspLiq(state.cells[i].presaux, tmed);
    rhol = (1 - betI) * liquidDensity + betI * state.cells[i].fluicol.MasEspFlu(state.cells[i].presaux, tmed);

    double gasDensity = flud.MasEspGas(state.cells[i].presaux, tmed);
    double vLiqTest = 1 + 0 * fabs(state.cells[i].Mliqini / (rhol * area));
    double vGasTest = 1 + 0 * fabs((state.cells[i].MC - state.cells[i].Mliqini) / (gasDensity * area));

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
        state.cells[i].QG = (state.cells[i].MC - state.cells[i].Mliqini) / gasDensity;
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
    // test!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // test!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
    double diameter = state.cells[i].duto.a;
    double area = 0.25 * M_PI * diameter * diameter;
    double perimeter = state.cells[i].duto.peri;
    double alfmed = state.cells[i].alf;
    double rhog = state.cells[i].flui.MasEspGas(state.cells[i].pres, state.cells[i].temp);
    double rhol = (1 - state.cells[i].bet) * state.cells[i].flui.MasEspLiq(state.cells[i].pres, state.cells[i].temp) + state.cells[i].bet * state.cells[i].fluicol.MasEspFlu(state.cells[i].pres, state.cells[i].temp);
    double ugsmed = (state.cells[i].QG) / (area);
    double ulsmed = state.cells[i].QL / (area);
    double j = ugsmed + ulsmed;
    double ABSjL = (fabs(state.cells[i - 1].QG) + fabs(state.cells[i - 1].QL)) / state.cells[i - 1].duto.area;

    double rhomix = alfmed * rhog + (1 - alfmed) * rhol;
    double viscmix = alfmed * state.cells[i].flui.ViscGas(state.cells[i].pres, state.cells[i].temp) + (1 - alfmed) * ((1 - state.cells[i].bet) * state.cells[i].flui.ViscOleo(state.cells[i].pres, state.cells[i].temp) + state.cells[i].bet * state.cells[i].fluicol.VisFlu(state.cells[i].pres, state.cells[i].temp));

    double reynolds;
    if (state.cells[i].duto.revest == 0)
        reynolds = state.cells[i].Rey(state.cells[i].duto.a, j, rhomix, viscmix);
    else {
        double dhid = 4 * area / perimeter;
        reynolds = state.cells[i].Rey(dhid, j, rhomix, viscmix);
    }
    double frictionFactor = state.cells[i].fric(reynolds, state.cells[i].duto.rug / diameter);
    double medpres = 0;
    double gradfric = (1 - medpres) * 0.5 * frictionFactor * rhomix * (fabs(j) * j) * perimeter * dx / area;
    double gradhidro = (1 - medpres) * kGravity * sin(state.cells[i].duto.teta) * rhomix * dx;

    state.cells[i].presauxini = state.cells[i].presaux;
    state.cells[i].presaux = state.cells[i].pres + (gradfric + gradhidro - state.cells[i - 1].dpB) / kPascalPerKgfPerCm2;
    state.cells[i].dpresaux = 0.5 * (gradfric + gradhidro - state.cells[i - 1].dpB) / kPascalPerKgfPerCm2;
    dx = 0.5 * state.cells[i].dxL;
    diameter = state.cells[i - 1].duto.a;
    area = 0.25 * M_PI * diameter * diameter;
    perimeter = state.cells[i - 1].duto.peri;
    alfmed = state.cells[i - 1].alf;
    rhog = state.cells[i - 1].flui.MasEspGas(state.cells[i - 1].pres, state.cells[i - 1].temp);
    rhol = (1 - state.cells[i - 1].bet) * state.cells[i - 1].flui.MasEspLiq(state.cells[i - 1].pres, state.cells[i - 1].temp) + state.cells[i - 1].bet * state.cells[i - 1].fluicol.MasEspFlu(state.cells[i - 1].pres, state.cells[i - 1].temp);
    ugsmed = (state.cells[i - 1].MC - state.cells[i - 1].Mliqini) / (area * rhog);
    ulsmed = state.cells[i - 1].Mliqini / (area * rhol);
    j = ugsmed + ulsmed;

    rhomix = alfmed * rhog + (1 - alfmed) * rhol;
    viscmix = alfmed * state.cells[i - 1].flui.ViscGas(state.cells[i - 1].pres, state.cells[i - 1].temp) + (1 - alfmed) * ((1 - state.cells[i - 1].bet) * state.cells[i - 1].flui.ViscOleo(state.cells[i - 1].pres, state.cells[i - 1].temp) + state.cells[i - 1].bet * state.cells[i - 1].fluicol.VisFlu(state.cells[i - 1].pres, state.cells[i - 1].temp));

    if (state.cells[i - 1].duto.revest == 0)
        reynolds = state.cells[i - 1].Rey(state.cells[i - 1].duto.a, j, rhomix, viscmix);
    else {
        double dhid = 4 * area / perimeter;
        reynolds = state.cells[i - 1].Rey(dhid, j, rhomix, viscmix);
    }
    frictionFactor = state.cells[i - 1].fric(reynolds, state.cells[i - 1].duto.rug / diameter);
    gradfric = (1 - medpres) * 0.5 * frictionFactor * rhomix * (fabs(j) * j) * perimeter * dx / area;
    gradhidro = (1 - medpres) * kGravity * sin(state.cells[i - 1].duto.teta) * rhomix * dx;

    state.cells[i].presaux = 0.5 * (state.cells[i].presaux) +
                        0.5 * (state.cells[i - 1].pres - (gradfric + gradhidro) / kPascalPerKgfPerCm2);
    state.cells[i].dpresaux -= 0.5 * (gradfric + gradhidro) / kPascalPerKgfPerCm2;
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
        betI = state.cells[i].bet; // beta doubt
    else
        betI = state.cells[i].betL;

    double liquidDensity = flud.MasEspLiq(state.cells[i].presaux, tmed);
    rhol = (1 - betI) * liquidDensity + betI * state.cells[i].fluicol.MasEspFlu(state.cells[i].presaux, tmed);

    double gasDensity = flud.MasEspGas(state.cells[i].presaux, tmed);

    state.cells[i].QL = state.cells[i].Mliqini / rhol;
    state.cells[i].QG = (state.cells[i].MC - state.cells[i].Mliqini) / gasDensity;
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
                betI = state.cells[i].bet; // beta doubt
            else
                betI = state.cells[i].betL;

            double liquidDensity = flud.MasEspLiq(state.cells[i].presaux, tmed);
            double rhol = (1 - betI) * liquidDensity + betI * state.cells[i].fluicol.MasEspFlu(state.cells[i].presaux, tmed);

            double gasDensity = flud.MasEspGas(state.cells[i].presaux, tmed);
            state.cells[i].QL = state.cells[i].Mliqini / rhol;
            state.cells[i].QG = (state.cells[i].MC - state.cells[i].Mliqini) / gasDensity;
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
                betI = state.cells[i].bet; // beta doubt
            else
                betI = state.cells[i].betL;

            double liquidDensity = flud.MasEspLiq(state.cells[i].presaux, tmed);
            double rhol = (1 - betI) * liquidDensity + betI * state.cells[i].fluicol.MasEspFlu(state.cells[i].presaux, tmed);

            double gasDensity = flud.MasEspGas(state.cells[i].presaux, tmed);

            state.cells[i].QL = state.cells[i].Mliqini / rhol;
            state.cells[i].QG = (state.cells[i].MC - state.cells[i].Mliqini) / gasDensity;
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

    int lastInteriorCell = state.lastCell - 1;
    state.bufferedCompletionMassSource = state.cells[lastInteriorCell + 1].fontemassCR;

    state.bufferedLiquidMassSource = state.cells[lastInteriorCell + 1].fontemassLR;

    state.bufferedGasMassSource = state.cells[lastInteriorCell + 1].fontemassGR;
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

    double rholp = state.cells[state.lastCell].rpC;
    double rholc = state.cells[state.lastCell].rcC;
    double rholmix = (1 - betSup) * rholp + betSup * rholc;
    double romix = alfSup * state.cells[state.lastCell].rgC + (1 - alfSup) * rholmix;

    double quality;
    if ((massgas >= 0 && state.cells[state.lastCell - 1].MliqiniR <= 0) || (massgas < 0 && state.cells[state.lastCell - 1].MliqiniR == 0))
        quality = 1.;
    else if (massgas <= 0 && state.cells[state.lastCell - 1].MliqiniR > 0)
        quality = 0.;
    else if (masentrada < 0) {
        quality = 1.;
    } else if (fabs(masentrada) < 1e-15)
        quality = alfSup * state.cells[state.lastCell].flui.MasEspGas(state.cells[state.lastCell].pres, state.cells[state.lastCell].temp) / romix;
    else
        quality = fabs(massgas / masentrada);
    if (quality > 1)
        quality = 1;
    if (state.outletPressure < state.gasSurfacePressure) {
        quality = 1.;
    }

    if (quality == 0 && alfSup > 0.05)
        quality = alfSup * state.cells[state.lastCell].flui.MasEspGas(state.cells[state.lastCell].pres, state.cells[state.lastCell].temp) / romix;

    romix = quality * (1. / state.cells[state.lastCell].rgC) + (1 - quality) * (1. / rholmix);
    romix = 1 / romix;

    double masChk;

    double sinal = 1.;
    double pmon = state.outletPressure;

    double ypres = state.gasSurfacePressure / state.outletPressure;
    if (surfaceChokeIsOpen(state) && ypres < 1.) {
        double cplM = (1. - betSup) * state.cells[state.lastCell].flui.CalorLiq(state.outletPressure, tESup) -
                      betSup * state.cells[state.lastCell].fluicol.CalorLiq(state.outletPressure, tESup);
        double jtlM = (1. - betSup) * state.cells[state.lastCell].flui.JTL(state.outletPressure, tESup) - betSup / rholc;
        double gasSpecificHeat = state.cells[state.lastCell].flui.CalorGas(state.outletPressure, tESup);
        double jtgM = state.cells[state.lastCell].flui.JTG(state.outletPressure, tESup);
        state.input.valTempChokeJus = tESup + ((1. - quality) * jtlM / cplM + quality * jtgM / gasSpecificHeat) * (state.gasSurfacePressure - state.outletPressure) * kPascalPerKgfPerCm2Variant;
    }
    if (ypres > 1.) {
        if (state.input.chkv == 0)
            sinal = -1.;
        else
            sinal = 0.;
        quality = 1.;
        pmon = state.gasSurfacePressure;
        ypres = 1. / ypres;
    }

    masChk = state.surfaceChoke.vazmassSachd(ypres, pmon, tESup, alfSup, betSup, quality, state.cells[state.lastCell - 1].flui,
                                   state.cells[state.lastCell - 1].fluicol);
    maxSup = state.surfaceChoke.vazmaxSachd(pmon, tESup, alfSup, betSup, quality, state.cells[state.lastCell - 1].flui, state.cells[state.lastCell - 1].fluicol);

    int fluxcri = 1;
    if (quality <= 0.01 || fabs(ypres) > fabs(state.surfaceChoke.razpres)) {
        fluxcri = 0;
        maxSup = masChk;
    }
    if (surfaceChokeIsShut(state))
        maxSup = 0.;

    double sinal2 = 1.;
    pmon = state.outletPressure * 1.0001;
    ypres = state.gasSurfacePressure / pmon;
    if (ypres > 1.) {
        sinal2 = -1.;
        quality = 1;
        pmon = state.gasSurfacePressure;
        ypres = 1. / ypres;
        if (state.input.chkv == 0)
            sinal2 = -1.;
        else
            sinal2 = 0.;
    }

    double masChk2 = state.surfaceChoke.vazmassSachd(ypres, pmon, tESup, alfSup, betSup, quality, state.cells[state.lastCell - 1].flui,
                                           state.cells[state.lastCell - 1].fluicol);
    double maxSup2 = state.surfaceChoke.vazmaxSachd(pmon, tESup, alfSup, betSup, quality, state.cells[state.lastCell - 1].flui,
                                          state.cells[state.lastCell - 1].fluicol);

    fluxcri = 1;
    if (quality <= 0.01 || fabs(ypres) > fabs(state.surfaceChoke.razpres)) {
        fluxcri = 0;
        maxSup2 = masChk2;
    }
    if (surfaceChokeIsShut(state))
        maxSup2 = 0.;
    double dmaxsup = (sinal2 * maxSup2 - sinal * maxSup) / (state.outletPressure * 0.0001);

    double masliq;
    double masgas;
    int abertoini = state.surfaceChokeOpen;
    double delp;
    if (surfaceChokeIsOpen(state))
        delp = (0.5 / kPascalPerKgfPerCm2) * (1 / romix) * (1 / (state.surfaceChoke.AreaGarg * state.surfaceChoke.AreaGarg * state.surfaceChoke.cdchk * state.surfaceChoke.cdchk)) * masentrada * masentrada;
    else
        delp = 0;

    int masChkSup0 = state.surfaceChokeMassFlag;
    state.chokeModeChanged = 0;

    double difdelp = state.outletPressure - state.gasSurfacePressure;
    if (((quality < 1e-7 && surfaceChokeIsOpen(state)) ||
         (quality < 0.01 && surfaceChokeIsOpen(state) &&
          fabs(difdelp) / delp < 1.2 && fabs(difdelp) / delp > 0.8 &&
          ((fabs(maxSup) > 0 && fabs((masentrada - maxSup) / maxSup) < 0.2) ||
           (fabs(masentrada) > 0 && fabs((masentrada - maxSup) / masentrada) < 0.2))))) {
        state.surfaceChokeOpen = 1;
        double sens = 1.;
        if (masentrada < 0. && state.cells[state.lastCell - 1].MliqiniR <= 0)
            sens = 0.;
        state.outletPressure = state.gasSurfacePressure + sens * delp;
        state.surfaceChokeMassFlag = 0;
        if (state.surfaceChokeMassFlag != masChkSup0)
            state.chokeModeChanged = 1;

    } else {
        if (surfaceChokeIsOpen(state) &&
            ((((*state.globals).lixo5 - 2 * state.input.dtmax > state.movingMeanTemperature || state.input.perm == 2) &&
              fabs(state.movingMeanPressure - state.gasSurfacePressure) / state.movingMeanPressure < 0.05 && fabs(state.movingMeanFlux) < 0.5) ||
             (((*state.globals).lixo5 - 2 * state.input.dtmax > state.movingMeanTemperature || state.input.perm == 2) &&
              fabs(state.movingMeanPressure - state.gasSurfacePressure) < (0.05 * state.gasSurfacePressure) && fabs(state.movingMeanFlux) < 5. && delp < 0.01 * state.gasSurfacePressure) ||
             (((*state.globals).lixo5 - 2 * state.input.dtmax > state.movingMeanTemperature || state.input.perm == 2) && (state.gasSurfacePressure - state.movingMeanPressure) / state.movingMeanPressure > 0.001 && state.input.chkv == 0) || (fabs(delp) < 0.1 && (state.outletPressure - state.gasSurfacePressure) / state.outletPressure < 0.05 && state.input.chkv == 0))) {
            state.surfaceChokeOpen = 1;
            if (abertoini != state.surfaceChokeOpen)
                state.openTime = 1;
        } else {
            state.surfaceChokeOpen = 0;
            if (state.openTime > 60)
                state.openTime = 0;
        }
        if ((((quality > -0.01 && state.cells[state.lastCell].alf > -0.01) || surfaceChokeIsShut(state)) && state.surfaceChoke.AreaGarg < 0.6 * state.cells[state.lastCell - 1].duto.area && (state.surfaceChokeOpen == 0 && (state.openTime == 0 || state.openTime > 60)))) {

            state.surfaceChokeOpen = 0;
            state.openTime = 0;
            masliq = sinal * maxSup * (1. - quality);
            masgas = sinal * maxSup * quality;

            state.cells[state.lastCell].DmasschokeG = -1 * (1. - quality) * dmaxsup;
            state.cells[state.lastCell].DmasschokeL = -1 * quality * ((1 - betSup) * rholp / rholmix) * dmaxsup;
            state.cells[state.lastCell].DmasschokeC = -1 * quality * (betSup * rholc / rholmix) * dmaxsup;

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
                state.outletPressure = state.gasSurfacePressure;
            else {

                state.surfaceChokeOpen = 1;
                state.surfaceChokeMassFlag = 0;
                if (state.surfaceChokeMassFlag != masChkSup0)
                    state.chokeModeChanged = 1;
                state.outletPressure = state.gasSurfacePressure;
            }
            state.surfaceChokeOpen = 1;
            state.surfaceChokeMassFlag = 0;
            if (state.surfaceChokeMassFlag != masChkSup0)
                state.chokeModeChanged = 1;
        }
    }

    if (state.surfaceChokeMassFlag == 0 && (*state.globals).chaverede == 1) {
        double betloc;
        if ((state.cells[state.lastCell - 1].MR - state.cells[state.lastCell - 1].MliqiniR) * 0 + 1 * state.cells[state.lastCell - 1].MliqiniR > 0.)
            betloc = state.cells[state.lastCell - 1].bet; // beta test // beta doubt
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

    double quality;
    if (massgas > 0 && state.cells[state.lastCell - 1].MliqiniRBuf < 0)
        quality = 1.;
    else if (massgas <= 0 && state.cells[state.lastCell - 1].MliqiniRBuf >= 0)
        quality = 0.;
    else if (masentrada < 0) {
        if ((*state.globals).chaverede == 0 || state.endNode == 1)
            quality = 1.;
        else {
            quality = titRev;
            alfSup = alfRev;
            betSup = betRev;
        }
    } else if (fabs(masentrada) < 1e-15)
        quality = alfSup * state.cells[state.lastCell].flui.MasEspGas(state.cells[state.lastCell].presBuf, state.cells[state.lastCell].temp) / romix;
    else
        quality = fabs(massgas / masentrada);
    if (quality > 1)
        quality = 1;
    if (state.cells[state.lastCell].presBuf < state.gasSurfacePressure) {
        if ((*state.globals).chaverede == 0 || state.endNode == 1)
            quality = 1.;
        else {
            quality = titRev;
            alfSup = alfRev;
            betSup = betRev;
        }
    }

    if (quality == 0 && alfSup > 0.05)
        quality = alfSup * state.cells[state.lastCell].flui.MasEspGas(state.cells[state.lastCell].presBuf, state.cells[state.lastCell].temp) / romix;

    romix = quality * (1. / state.cells[state.lastCell].rgC) + (1 - quality) * (1. / rholmix);
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
            quality = 1.;
        else {
            quality = titRev;
            alfSup = alfRev;
            betSup = betRev;
        }
        pmon = state.gasSurfacePressure;
        ypres = 1. / ypres;
    }

    masChk = state.surfaceChoke.vazmassSachd(ypres, pmon, tESup, alfSup, betSup, quality, state.cells[state.lastCell - 1].flui,
                                   state.cells[state.lastCell - 1].fluicol);
    maxSup = state.surfaceChoke.vazmaxSachd(pmon, tESup, alfSup, betSup, quality, state.cells[state.lastCell - 1].flui,
                                  state.cells[state.lastCell - 1].fluicol);

    int fluxcri = 1;
    if (quality <= 0.01 || fabs(ypres) > fabs(state.surfaceChoke.razpres)) {
        fluxcri = 0;
        maxSup = masChk;
    }
    if (surfaceChokeIsShut(state))
        maxSup = 0.;

    double masliq;
    double masgas;
    int abertoini = state.surfaceChokeOpen;

    double delp;
    if (surfaceChokeIsOpen(state))
        delp = (0.5 / kPascalPerKgfPerCm2) * (1 / romix) *
               (1 / (state.surfaceChoke.AreaGarg * state.surfaceChoke.AreaGarg * state.surfaceChoke.cdchk * state.surfaceChoke.cdchk)) * masentrada * masentrada;
    else
        delp = 0.;

    double difdelp = fabs(fabs(state.cells[state.lastCell].presBuf - state.gasSurfacePressure) - delp);
    if (((quality < 1e-7 && surfaceChokeIsOpen(state)) ||
         (quality < 0.01 && surfaceChokeIsOpen(state) &&
          difdelp / delp < 0.2 &&
          ((fabs(maxSup) > 0 && fabs((masentrada - maxSup) / maxSup) < 0.2) ||
           (fabs(masentrada) > 0 && fabs((masentrada - maxSup) / masentrada) < 0.2))))) {
        state.surfaceChokeOpen = 1;
        if (state.input.chkv == 0 || state.cells[state.lastCell - 1].MRBuf > 0.)
            sinal = 1.;
        else
            sinal = 0.;
        double betloc;
        if ((state.cells[state.lastCell - 1].MRBuf - state.cells[state.lastCell - 1].MliqiniRBuf) * 0 + 1 * state.cells[state.lastCell - 1].MliqiniRBuf > 0.)
            betloc = state.cells[state.lastCell - 1].bet; // beta test // beta doubt
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
            abertoini = state.surfaceChokeOpen;
        }
        if (((quality > -0.01 && state.cells[state.lastCell].alf > -0.01) || surfaceChokeIsShut(state)) &&
            state.surfaceChoke.AreaGarg < 0.6 * state.cells[state.lastCell - 1].duto.area && (state.surfaceChokeOpen == 0 && (state.openTime == 0 || state.openTime > 60))) {

            masliq = sinal * maxSup * (1. - quality);
            masgas = sinal * maxSup * quality;
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
                betloc = state.cells[state.lastCell - 1].bet; // beta test // beta doubt
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
        state.surfaceChokeOpen = 1;
    if (fabs(state.outletPressure - state.gasSurfacePressure) / state.outletPressure < 0.05 && state.surfaceChokeOpen == 0 && state.timeStep > 1. && state.surfaceChoke.AreaGarg < 0.6 * state.cells[state.lastCell - 1].duto.area &&
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
    if ((((*state.globals).lixo5 - 2 * state.input.dtmax) > state.movingMeanTemperature && (fabs(state.movingMeanPressure - state.outletPressure) / state.outletPressure > 0.4)) || (state.outletPressure < state.gasSurfacePressure && state.surfaceChokeMassFlag == 1 && state.surfaceChoke.AreaGarg > (1e-3) * state.cells[state.lastCell - 1].duto.area) || (state.cells[state.lastCell].alf <= 0.1 && state.surfaceChoke.AreaGarg < 0.6 * state.cells[state.lastCell - 1].duto.area && state.surfaceChoke.AreaGarg > 1e-3 * state.cells[state.lastCell - 1].duto.area && state.surfaceChokeOpen == 0)) {
        double denominador = 10.;
        double progres = 0.95;

        if ((state.outletPressure < state.gasSurfacePressure && state.surfaceChokeMassFlag == 1)) {
            denominador = 100.;
            progres = 0.5;
            state.multiplier = pow(progres, 10. * fabs(state.gasSurfacePressure - state.outletPressure) / state.outletPressure) * state.multiplier;
        } else if ((state.cells[state.lastCell].alf <= 0.1 && state.surfaceChoke.AreaGarg < 0.6 * state.cells[state.lastCell - 1].duto.area)) {
            denominador = 1000.;
            progres = 0.5;
            state.multiplier = (state.cells[state.lastCell].alf + 1e-5) * state.multiplier;
        } else {
            state.multiplier = pow(progres, 10 * fabs(state.movingMeanPressure - state.outletPressure) / state.outletPressure) * state.multiplier;
        }
        if (state.input.dtmax * state.multiplier < state.smallestCellLength / denominador)
            state.multiplier = (state.smallestCellLength / denominador) / state.input.dtmax;
        if (state.multiplier > 0.8)
            state.multiplier = 0.8;
    } else {
        if (state.cells[state.lastCell].alf >= 0.11 || state.surfaceChoke.AreaGarg >= 0.6 * state.cells[state.lastCell - 1].duto.area ||
            state.surfaceChoke.AreaGarg <= 1e-3 * state.cells[state.lastCell - 1].duto.area || state.surfaceChokeOpen == 1) {
            state.multiplier = state.multiplier / 0.95;
            if (state.multiplier > 0.8)
                state.multiplier = 0.8;
        }
    }

    for (int i = 0; i <= state.lastCell; i++) {
        double jmix = 0.;
        double pipeArea = state.cells[i].duto.area;
        double dtaux;
        // celula[i].Mliqini
        double alfteste = state.cells[i].alf;
        if (i > 0 && (state.cells[i].MC - state.cells[i].Mliqini) > 0)
            alfteste = state.cells[i - 1].alf;
        if (alfteste > 1e-9)
            jmix += fabs(
                (state.cells[i].MC - state.cells[i].Mliqini) / (pipeArea * state.cells[i].flui.MasEspGas(state.cells[i].pres, state.cells[i].temp) * (0 + 1 * alfteste)));
        if (fabs(state.cells[i].VTemper) > jmix)
            jmix = fabs(state.cells[i].VTemper);

        double jmixL = 0.;
        if (alfteste < 1 - 1e-5)
            jmixL = fabs(
                state.cells[i].Mliqini / (pipeArea * ((1. - state.cells[i].bet) * state.cells[i].rpC + state.cells[i].bet * state.cells[i].rcC) * (1. - 0 * alfteste)));
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
        if (state.cells[i].acsr.tipo == kAccessoryRadialPorous) {
            state.cells[i].acsr.radialPoro.dt = state.timeStep;
        } else if (state.cells[i].acsr.tipo == kAccessoryPorous2D) {
            state.cells[i].acsr.poroso2D.dt = state.timeStep;
        }
    }
    // dtInter=dt; // change 2
}

void computeExplicitTimeStep(const TransientStepState &state) {

    for (int i = 0; i <= state.lastCell; i++) {
        double velAux = state.cells[i].termAdSomVel();
        double speedOfSound = state.cells[i].somVel();
        double velpropag1 = velAux + speedOfSound;
        double velpropag2 = fabs(velAux - speedOfSound);
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
        if (state.cells[i].acsr.tipo == kAccessoryChoke && state.cells[i].acsr.chk.AreaGarg <= 1e-15 * state.cells[i].acsr.chk.AreaTub)
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
        if (state.cells[i].acsr.tipo == kAccessoryRadialPorous) {
            state.cells[i].acsr.radialPoro.defineDT(0);
            if (state.cells[i].acsr.radialPoro.dt < state.timeStep)
                state.timeStep = state.cells[i].acsr.radialPoro.dt;
        }
        if (state.cells[i].acsr.tipo == kAccessoryPorous2D) {
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
        int linAberta = 1; // several-valve case
        for (int i = 0; i <= state.input.nvalv; i++)
            if (state.masterRatio1[i] <= 1e-3)
                linAberta = 0; // several-valve case
        if ((linAberta == 1 && state.surfaceChoke.AreaGarg / state.cells[state.lastCell - 1].duto.area > 1e-3))
            state.fullModel = 1; // several-valve case
    }
    for (int i = 0; i <= state.input.nvalv; i++)
        if (state.masterRatio1[i] != state.masterRatio0[i])
            state.fullModel = 0; // several-valve case
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
                    state.cells[celpos].fontemassGL += state.desperationMoment * 10000 * state.cells[celpos].flui.Deng * kAirDensityAtStandardConditions / kSecondsPerDay;
                    state.cells[celpos - 1].fontemassGR = state.cells[celpos].fontemassGL;
                }
                if (state.cells[celpos].alf < 0.5)
                    state.masterCriticalRatio[i] = 0.01;
            } else
                state.desperationMoment = 0.;
            if (state.masterRatio1[i] < state.masterRatio0[i] && (state.masterRatio1[i] <= 1.1 * state.input.master1.razareaativ && state.masterRatio1[i] >= state.masterCriticalRatio[i] * state.input.master1.razareaativ)) {
                double timeStepDivisor = 20.;
                if (state.cells[celpos].alf < 0.5)
                    timeStepDivisor = 40.;
                dtvec[i] = dtaux;
                if (dtvec[i] > 1)
                    dtvec[i] = 1.;
                if (state.cells[celpos].alf < 0.5)
                    dtvec[i] = 0.1 + 0.9 * (state.cells[celpos].alf) / 0.5;
                dtvec[i] /= timeStepDivisor;
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
            state.cells[0].betI = state.cells[0].betPigE; // beta test
        state.cells[0].betLI = state.cells[0].betI;
    }
    for (int i = 1; i <= state.lastCell; i++) {
        state.cells[i].betI = state.cells[i].betPigE;
        if (state.cells[i].QL > 0.)
            state.cells[i].betI = state.cells[i - 1].betPigD; // beta test
        state.cells[i - 1].betRI = state.cells[i].betI;
        if (i < state.lastCell)
            state.cells[i + 1].betLI = state.cells[i].betI;
    }
}

void solvePressureVolumeCoupling(const TransientStepState &state, int vexpli, int ciclo) {
#pragma omp parallel for num_threads((*state.globals).ntrd)
    for (int i = 0; i <= state.lastCell; i++) {
        state.cells[i].GeraLocal(state.outletPressure, state.surfaceChokeMassFlag, state.lastCell, state.input.master1.razareaativ, state.inletPressure, state.inletTemperature, state.inletQuality, state.inletCompletionFraction, ciclo,
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
    	generateFluidMiniTables(state);
    double betIV;
    double solutionGasRatioInSitu;
    double oilVolumeFactorInSitu;
    double waterVolumeFactorInSitu;
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
        solutionGasRatioInSitu = state.cells[i].flui.RS(state.cells[i].pres, state.cells[i].temp);
        oilVolumeFactorInSitu = state.cells[i].flui.BOFunc(state.cells[i].pres, state.cells[i].temp, solutionGasRatioInSitu);
        waterVolumeFactorInSitu = state.cells[i].flui.BAFunc(state.cells[i].pres, state.cells[i].temp);
        bswV = state.cells[i].flui.BSW * waterVolumeFactorInSitu / (oilVolumeFactorInSitu + waterVolumeFactorInSitu * state.cells[i].flui.BSW - state.cells[i].flui.BSW * oilVolumeFactorInSitu);
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
            double rgST = state.cells[0].flui.Deng * kAirDensityAtStandardConditions;
            double roST = 141.5 * 1000. / (131.5 + state.cells[0].flui.API);
            double gasDensity = state.cells[0].flui.MasEspGas(state.inletPressure, state.inletTemperature);
            double liquidDensity = state.cells[0].flui.MasEspLiq(state.inletPressure, state.inletTemperature);
            double titH = state.cells[0].flui.FracMassHidra(state.inletPressure, state.inletTemperature);
            double rcST = state.cells[0].fluicol.MasEspFlu(1.01, 20.);
            double completionDensity = state.cells[0].fluicol.MasEspFlu(state.inletPressure, state.inletTemperature);
            double rlMix = state.inletCompletionFraction * completionDensity + (1. - state.inletCompletionFraction) * liquidDensity;
            double val1 = ((1. - state.inletCompletionFraction) * liquidDensity * titH / (1. - titH));
            state.inletQuality = val1 / (rlMix + val1);
        } else if (state.input.tipoFluido == 1) {
            double rgST = state.cells[0].flui.Deng * kAirDensityAtStandardConditions;
            double roST = 141.5 * 1000. / (131.5 + state.cells[0].flui.API);
            double gasDensity = state.cells[0].flui.MasEspGas(state.inletPressure, state.inletTemperature);
            double liquidDensity = state.cells[0].flui.MasEspoleo(state.inletPressure, state.inletTemperature);
            double quality = state.cells[0].flui.FracMass(state.inletPressure, state.inletTemperature);
            double rcST = state.cells[0].fluicol.MasEspFlu(1.01, 20.);
            double completionDensity = state.cells[0].fluicol.MasEspFlu(state.inletPressure, state.inletTemperature);
            double val1 = (rcST / completionDensity) * (gasDensity / rgST) * state.input.CCPres.bet[0] / quality;
            double val2 = (gasDensity / liquidDensity) * (1 - quality) / quality;
            double titT = gasDensity / (((1. - quality) / quality) * (gasDensity / liquidDensity) + gasDensity + val1);
            state.inletQuality = titT;
            state.inletCompletionFraction = val1 / (val2 + val1);
        }
    }
}

namespace {

// ---------------------------------------------- time-step policy registry --
//
// The restrictions on the time step, as a list of policies instead of calls
// written into the step. A policy names the point of the step where it acts
// (its hook), the condition under which it acts, and the action. The step asks
// the registry for the policies of a hook, and they run in the order they are
// listed. Adding a policy is writing it and listing it in TimeStepPolicies:
// neither the step nor any other policy changes.
//
// Resolved at compile time: `if constexpr` drops the policies of every other
// hook, so a hook costs what its own conditions cost and nothing is called
// through a pointer.
//
// The registry covers the step this module runs, a single line (chaverede ==
// 0). A network run is sequenced by Num4Main.cpp, which calls the same three
// restrictions through SProd's surface, so a policy listed here does not
// reach a network run.
enum class TimeStepHook {
    CouplingIterationStart, // start of a pressure-volume coupling iteration
    AfterPigUpdate,         // after the pig moves, in the same iteration
    AfterValveOpenings,     // after the step re-reads the valve openings
};

/// What a policy may read. couplingIteration is kontaAcop inside the coupling
/// loop, and kOutsideCouplingLoop at a hook outside it.
struct TimeStepPolicyContext {
    const TransientStepState &step;
    int couplingIteration;
    int explicitScheme;
};

constexpr int kOutsideCouplingLoop = -1;

/// What the registry requires of a policy. A listed type that lacks one of the
/// three fails here, in one line, instead of inside the fold below.
template <typename Policy>
concept TimeStepPolicy = requires(const TimeStepPolicyContext &context) {
    { Policy::hook } -> std::convertible_to<TimeStepHook>;
    { Policy::applies(context) } -> std::same_as<bool>;
    { Policy::apply(context) } -> std::same_as<void>;
};

/// With valve time-step control on, the valves' travel restricts the time
/// step, once per step, at the first coupling iteration.
struct RestrictTimeStepByValvePolicy {
    static constexpr TimeStepHook hook = TimeStepHook::CouplingIterationStart;
    static bool applies(const TimeStepPolicyContext &context) {
        return context.couplingIteration == 0 && context.step.input.controleDTvalv == 1;
    }
    static void apply(const TimeStepPolicyContext &context) {
        restrictTimeStepByValve(context.step);
    }
};

/// Damps the maximum time step at the coupling iteration numbered by the
/// full-model flag: the only iteration of the reduced model, the second of the
/// full one.
struct DampMaximumTimeStepPolicy {
    static constexpr TimeStepHook hook = TimeStepHook::AfterPigUpdate;
    static bool applies(const TimeStepPolicyContext &context) {
        return context.couplingIteration == 1 * context.step.fullModel;
    }
    static void apply(const TimeStepPolicyContext &context) {
        dampMaximumTimeStep(context.step);
    }
};

/// With the full model on, the rates of change of pressure and temperature
/// decide whether it stays on, and with it how many coupling iterations the
/// step takes.
struct PressureRateOfChangePolicy {
    static constexpr TimeStepHook hook = TimeStepHook::AfterValveOpenings;
    static bool applies(const TimeStepPolicyContext &context) {
        return context.step.fullModel == 1;
    }
    static void apply(const TimeStepPolicyContext &context) {
        evaluatePressureRateOfChange(context.step, 0, 0, context.explicitScheme);
    }
};

/// The policies of a hook, applied in the order they are listed.
template <TimeStepPolicy... Policies>
struct TimeStepPolicyRegistry {
    template <TimeStepHook Hook>
    static void apply(const TimeStepPolicyContext &context) {
        (applyIfAt<Policies, Hook>(context), ...);
    }

  private:
    template <typename Policy, TimeStepHook Hook>
    static void applyIfAt(const TimeStepPolicyContext &context) {
        if constexpr (Policy::hook == Hook) {
            if (Policy::applies(context))
                Policy::apply(context);
        }
    }
};

/// Every time-step policy of the step. Registering one is listing it here.
using TimeStepPolicies = TimeStepPolicyRegistry<
    RestrictTimeStepByValvePolicy,
    DampMaximumTimeStepPolicy,
    PressureRateOfChangePolicy>;

void advanceCouplingIteration(const TransientSolveState &state, int kontaAcop, int celpos, int vExpli, int ciclomax, double titRev, double alfRev, double betRev) {
    if (state.step.fullModel == 0) {
        for (int i = 0; i <= state.step.lastCell; i++)
            state.step.cells[i].m2d = 0.;
    } else {
        for (int i = 0; i <= state.step.lastCell; i++) {
            double area = state.step.cells[i].duto.area;
            double vLiqTest = fabs(state.step.cells[i].QL / (area));
            double vGasTest = fabs(state.step.cells[i].QG / (area));
            double razDp = 0.1;
            double razDT = 1;
            if (i < celpos && state.step.cells[celpos].acsr.chk.AreaGarg < 1e-15 * state.step.cells[celpos].acsr.chk.AreaTub) {
                razDT = 1;
            } else if (i == celpos + 1 && state.step.cells[celpos].acsr.chk.AreaGarg < 1e-15 * state.step.cells[celpos].acsr.chk.AreaTub) {
                razDT = 1;
            }
            if ((fabs(state.step.cells[i].dpdtIni) / state.step.cells[i].pres < razDp) && fabs(state.step.cells[i].dTdtIni) < razDT) {
                if (state.massTransferModel == 0)
                    state.step.cells[i].m2d = 1.;
                else
                    state.step.cells[i].m2d = 0.;
                state.step.cells[i].mudaDT = 1.;
            } else {
                state.step.cells[i].m2d = 0.;
                state.step.cells[i].mudaDT = 0.;
            }
        }
    }
    if (state.step.input.estabCol == 1) {
        for (int i = 0; i <= celpos; i++) {
            state.step.cells[i].m2d = 0.;
            state.step.cells[i].mudaDT = 0.;
            state.step.cells[i].estabCol = 1;
        }
    }
    evolveFractions(state.step, alfRev, betRev, kontaAcop);
    for (int i = 0; i <= state.step.lastCell; i++) {
        if (state.step.cells[i].acsr.tipo == kAccessoryRadialPorous) {
            state.step.cells[i].acsr.radialPoro.avancoSW(state.step.timeStep);
            if (state.step.cells[i].acsr.radialPoro.reinicia == -1) {
                if (state.step.restart > -1)
                    state.step.restart = -1;
                // celula[i].acsr.radialPoro.reavaliaDT(Ndt)
            }
        } else if (state.step.cells[i].acsr.tipo == kAccessoryPorous2D) {
            state.step.cells[i].acsr.poroso2D.avancoSW(state.step.timeStep);
            if (state.step.cells[i].acsr.poroso2D.reinicia == -1) {
                if (state.step.restart > -1)
                    state.step.restart = -1;
                // celula[i].acsr.radialPoro.reavaliaDT(Ndt)
            }
        }
    }

    if (state.step.input.correcaoMassaEspLiq == 1) {
        for (int i = 0; i < state.step.lastCell; i++)
            state.step.cells[i + 1].mudaDTL = state.step.cells[i].mudaDT;
    }

    // master-only case
    // master-only case
    TimeStepPolicies::apply<TimeStepHook::CouplingIterationStart>(
        {state.step, kontaAcop, vExpli}); // several-valve case
    if (state.step.restart == -1) {
        restartFractionEvolutionInitial(state.step);
        for (int i = 0; i <= state.step.lastCell; i++) {
            if (state.step.cells[i].acsr.tipo == kAccessoryRadialPorous) {
                state.step.cells[i].acsr.radialPoro.reavaliaDT(state.step.timeStep);
            } else if (state.step.cells[i].acsr.tipo == kAccessoryPorous2D) {
                state.step.cells[i].acsr.poroso2D.reavaliaDT(state.step.timeStep);
            }
        }
        for (int i = 0; i <= state.step.lastCell; i++) {
            if (state.step.cells[i].acsr.tipo == kAccessoryRadialPorous) {
                state.step.cells[i].acsr.radialPoro.reiniciaEvoluiSW(state.step.timeStep);
            }
            if (state.step.cells[i].acsr.tipo == kAccessoryPorous2D) {
                state.step.cells[i].acsr.poroso2D.reiniciaEvoluiSW(state.step.timeStep);
            }
        }
        state.step.finalAuxiliaryTimeStep = state.step.timeStep;
        restartFractionEvolution(state.step);
        evolveFractions(state.step, alfRev, betRev, kontaAcop);
        state.step.restart = 0;
        for (int i = 0; i <= state.step.lastCell; i++) {
            if (state.step.cells[i].acsr.tipo == kAccessoryRadialPorous) {
                state.step.cells[i].acsr.radialPoro.avancoSWcorrec();
            } else if (state.step.cells[i].acsr.tipo == kAccessoryPorous2D) {
                state.step.cells[i].acsr.poroso2D.avancoSWcorrec();
            }
        }
    }
    updatePig(state.step);

    if (kontaAcop == 0)
        state.minimumCycleTimeStep = state.step.timeStep;

    TimeStepPolicies::apply<TimeStepHook::AfterPigUpdate>({state.step, kontaAcop, vExpli});

    double gasMassSource = 0.;
    double liquidMassSource = 0.;
    double completionMassSource = 0.;
    if (state.step.fullModel == 1) {
        completionMassSource = state.step.cells[state.step.lastCell].fontemassCR;
        liquidMassSource = state.step.cells[state.step.lastCell].fontemassLR;
        gasMassSource = state.step.cells[state.step.lastCell].fontemassGR;
    }

    applyOutletPressureCondition(state.step, titRev, alfRev, betRev);
    state.updaters.updateThermal();

    if (state.step.cells[state.step.lastCell].alf < 0.05 && state.step.surfaceChokeMassFlag == 1)
        state.step.cells[state.step.lastCell].alf = 0.05;
    // several-valve case
    for (int j = 0; j <= state.step.input.nvalv; j++) {
        int celposAux;
        if (j > 0)
            celposAux = state.step.input.valv[j - 1].posicP;
        else
            celposAux = celpos;
        if (state.step.cells[celposAux].alf < 0.05 && state.step.masterRatio1[j] <= state.step.input.master1.razareaativ)
            state.step.cells[celposAux].alf = 0.05;
    }
    // several-valve case
    solvePressureVolumeCoupling(state.step, vExpli);

    if (kontaAcop < 1 * state.step.fullModel) {
        for (int i = 0; i <= state.step.lastCell; i++) {
            state.step.cells[i].dpdt = 1 * (state.step.productionSolution[2 * i + 1] - state.step.cells[i].pres) / state.step.cells[i].dt;
            state.step.cells[i].dpdtIni = state.step.cells[i].dpdt;
        }
    }
    if (kontaAcop == 1 * state.step.fullModel || state.step.input.cicloAcopTerm == 1) {
        updateCells(state.step);
    }
    if (state.step.input.cicloAcopTerm == 1 && state.step.fullModel == 1) {
        if (kontaAcop < 1 * state.step.fullModel)
            for (int i = 0; i <= state.step.lastCell; i++)
                state.step.cells[i].dpdt = state.step.cells[i].d2pdt2;
        state.updaters.marchTransientEnergy(kontaAcop, ciclomax);
    }
    if (kontaAcop != 1 * state.step.fullModel) {
        for (int i = 0; i <= state.step.lastCell; i++) {
            state.step.cells[i].FeiticoDoTempo2();
            if (state.step.cells[i].acsr.tipo == kAccessoryRadialPorous) {
                state.step.cells[i].acsr.radialPoro.FeiticoDoTempoSW();
            } else if (state.step.cells[i].acsr.tipo == kAccessoryPorous2D) {
                state.step.cells[i].acsr.poroso2D.FeiticoDoTempoSW();
            }
        }

        state.step.cells[state.step.lastCell].fontemassCR = completionMassSource;
        state.step.cells[state.step.lastCell].fontemassLR = liquidMassSource;
        state.step.cells[state.step.lastCell].fontemassGR = gasMassSource;

        state.step.surfaceChokeOpen = state.step.initiallyOpen;
        state.step.openTime = state.initialOpenTime;
    }
}

void writeProfiles(const TransientSolveState &state, int nrede) {
    if (state.step.input.nperfisp > 0) {
        if (((*state.step.globals).lixo5 > (*state.step.globals).localtiny && (*state.step.globals).lixo5 <= state.step.input.profp.tempo[state.productionProfileTimeCounter] && (*state.step.globals).lixo5 + state.step.timeStep >= state.step.input.profp.tempo[state.productionProfileTimeCounter])) {
            state.step.input.imprimeProfile(state.step.cells, state.step.productionFreeTerms, (*state.step.globals).lixo5, state.step.branchIndex, nrede);
            state.step.input.profp.tempo[state.productionProfileTimeCounter] = (*state.step.globals).lixo5;
            state.productionProfileTimeCounter++;
            if (state.productionProfileTimeCounter >= state.step.input.profp.n)
                state.productionProfileTimeCounter--;
        }
    }
    if (state.step.input.nperfisg > 0 && state.step.input.lingas > 0) {
        if (((*state.step.globals).lixo5 > (*state.step.globals).localtiny && (*state.step.globals).lixo5 <= state.step.input.profg.tempo[state.gasProfileTimeCounter] && (*state.step.globals).lixo5 + state.step.timeStep >= state.step.input.profg.tempo[state.gasProfileTimeCounter])) {
            state.step.input.imprimeProfileG(state.step.gasCells, state.step.gasFreeTerms, (*state.step.globals).lixo5, state.step.branchIndex, nrede);
            state.step.input.profg.tempo[state.gasProfileTimeCounter] = (*state.step.globals).lixo5;
            state.gasProfileTimeCounter++;
            if (state.gasProfileTimeCounter >= state.step.input.profg.n)
                state.gasProfileTimeCounter--;
        }
    }
    if (state.step.input.nperfistransp > 0) {
        if (((*state.step.globals).lixo5 > (*state.step.globals).localtiny && (*state.step.globals).lixo5 <= state.step.input.proftransp.tempo[state.productionCrossSectionProfileTimeCounter] && (*state.step.globals).lixo5 + state.step.timeStep >= state.step.input.proftransp.tempo[state.productionCrossSectionProfileTimeCounter])) {
            state.step.input.imprimeProfileTrans(state.step.cells, state.step.productionCrossSectionCount, (*state.step.globals).lixo5, state.step.branchIndex, nrede);
            state.step.input.proftransp.tempo[state.productionCrossSectionProfileTimeCounter] = (*state.step.globals).lixo5;
            state.productionCrossSectionProfileTimeCounter++;
            if (state.productionCrossSectionProfileTimeCounter >= state.step.input.proftransp.n)
                state.productionCrossSectionProfileTimeCounter--;
        }
    }
    if (state.step.input.nperfistransg > 0 && state.step.input.lingas > 0) {
        if (((*state.step.globals).lixo5 > (*state.step.globals).localtiny && (*state.step.globals).lixo5 <= state.step.input.proftransg.tempo[state.gasCrossSectionProfileTimeCounter] && (*state.step.globals).lixo5 + state.step.timeStep >= state.step.input.proftransg.tempo[state.gasCrossSectionProfileTimeCounter])) {
            state.step.input.imprimeProfileTransG(state.step.gasCells, state.gasCrossSectionCellCounts, (*state.step.globals).lixo5, state.step.branchIndex, nrede);
            state.step.input.proftransg.tempo[state.gasCrossSectionProfileTimeCounter] = (*state.step.globals).lixo5;
            state.gasCrossSectionProfileTimeCounter++;
            if (state.gasCrossSectionProfileTimeCounter >= state.step.input.proftransg.n)
                state.gasCrossSectionProfileTimeCounter--;
        }
    }
    if (state.step.input.nCelUnit>0) {
    	for(int iCelU=0;iCelU<state.step.input.nCelUnit;iCelU++){
    		if (((*state.step.globals).lixo5 > (*state.step.globals).localtiny && (*state.step.globals).lixo5 <= state.step.input.celUnit[iCelU].tempo[state.unitCellTimeCounters[iCelU]] &&
    				(*state.step.globals).lixo5 + state.step.timeStep >= state.step.input.celUnit[iCelU].tempo[state.unitCellTimeCounters[iCelU]])) {
    			state.step.input.relatorioCelulaUnitaria(state.step.cells,state.step.input.celUnit[iCelU].posicP, state.step.branchIndex,nrede);
    			state.step.input.celUnit[iCelU].tempo[state.unitCellTimeCounters[iCelU]] = (*state.step.globals).lixo5;
    			state.unitCellTimeCounters[iCelU]++;
    			if (state.unitCellTimeCounters[iCelU] >= state.step.input.celUnit[iCelU].parserie)
    				state.unitCellTimeCounters[iCelU]--;
    		}
    	}
    }
}

void writeTrends(const TransientSolveState &state, int ordemImpT, double velmaxdesc, int nrede) {
    if (state.step.input.ntendp > 0) {
        for (int i = 0; i < state.step.input.ntendp; i++) {
            if (state.productionTrendResetTimers[i] == 0) {
                if ((*state.step.globals).lixo5 < 1e-15)
                    state.updaters.writeProductionTrendHeader(i, nrede);
                if ((*state.step.globals).lixo5 > 1e-15)
                    state.step.input.imprimeTrend(state.step.cells, state.productionTrendMatrix[i], (*state.step.globals).lixo5, i, state.productionTrendCounts[i]);
                state.productionTrendCounts[i]++;
            }
            if (ordemImpT == 1) {
                if ((*state.step.globals).lixo5 >= 800.8000000000000445) {
                    int para;
                    para == 1;
                }
                state.updaters.writeProductionTrendRows(i, nrede);
                state.productionTrendBufferedCounts[i] = state.productionTrendCounts[i];
            }
            state.productionTrendResetTimers[i] += state.step.timeStep;
            if (state.productionTrendResetTimers[i] > state.step.input.trendp[i].dt)
                state.productionTrendResetTimers[i] = 0;
        }
    }
    if (state.step.input.ntendg > 0 && state.step.input.lingas > 0) {
        for (int i = 0; i < state.step.input.ntendg; i++) {
            if (state.gasTrendResetTimers[i] == 0 || (*state.step.globals).lixo5 < 1e-15) {
                if ((*state.step.globals).lixo5 < 1e-15)
                    state.updaters.writeGasTrendHeader(i, nrede);
                state.step.input.imprimeTrendG(state.step.gasCells, state.gasTrendMatrix[i], (*state.step.globals).lixo5, i, state.gasTrendCounts[i], velmaxdesc);
                state.gasTrendCounts[i]++;
            }
            if (ordemImpT == 1) {
                state.updaters.writeGasTrendRows(i, nrede);
                state.gasTrendBufferedCounts[i] = state.gasTrendCounts[i];
            }
            state.gasTrendResetTimers[i] += state.step.timeStep;
            if (state.gasTrendResetTimers[i] > state.step.input.trendg[i].dt)
                state.gasTrendResetTimers[i] = 0;
        }
    }
    if (state.step.input.ntendtransp > 0) {
        for (int i = 0; i < state.step.input.ntendtransp; i++) {
            if (state.productionCrossSectionTrendResetTimers[i] == 0 || (*state.step.globals).lixo5 < 1e-15) {
                if ((*state.step.globals).lixo5 < 1e-15)
                    state.updaters.writeProductionCrossSectionTrendHeader(i);
                state.productionCrossSectionTrendMatrix[i][state.productionCrossSectionTrendCounts[i]][0] = (*state.step.globals).lixo5;
                int poscel = state.step.input.trendtransp[i].posic;
                int poscam = state.step.input.trendtransp[i].camada - 1;
                int posdiscre = state.step.input.trendtransp[i].discre - 1;
                state.productionCrossSectionTrendMatrix[i][state.productionCrossSectionTrendCounts[i]][1] = state.step.cells[poscel].calor.Tcamada[poscam][posdiscre];
                state.productionCrossSectionTrendCounts[i]++;
            }
            if (ordemImpT == 1) {
                state.updaters.writeProductionCrossSectionTrendRows(i);
                state.productionCrossSectionTrendBufferedCounts[i] = state.productionCrossSectionTrendCounts[i];
            }
            state.productionCrossSectionTrendResetTimers[i] += state.step.timeStep;
            if (state.productionCrossSectionTrendResetTimers[i] > state.step.input.trendtransp[i].dt)
                state.productionCrossSectionTrendResetTimers[i] = 0;
        }
    }
    if (state.step.input.ntendtransg > 0 && state.step.input.lingas > 0) {
        for (int i = 0; i < state.step.input.ntendtransg; i++) {
            if (state.gasCrossSectionTrendResetTimers[i] == 0 || (*state.step.globals).lixo5 < 1e-15) {
                if ((*state.step.globals).lixo5 < 1e-15)
                    state.updaters.writeGasCrossSectionTrendHeader(i);
                state.gasCrossSectionTrendMatrix[i][state.gasCrossSectionTrendCounts[i]][0] = (*state.step.globals).lixo5;
                int poscel = state.step.input.trendtransg[i].posic;
                int poscam = state.step.input.trendtransg[i].camada - 1;
                int posdiscre = state.step.input.trendtransg[i].discre - 1;
                state.gasCrossSectionTrendMatrix[i][state.gasCrossSectionTrendCounts[i]][1] = state.step.gasCells[poscel].calor.Tcamada[poscam][posdiscre];
                state.gasCrossSectionTrendCounts[i]++;
            }
            if (ordemImpT == 1) {
                state.updaters.writeGasCrossSectionTrendRows(i);
                state.gasCrossSectionTrendBufferedCounts[i] = state.gasCrossSectionTrendCounts[i];
            }
            state.gasCrossSectionTrendResetTimers[i] += state.step.timeStep;
            if (state.gasCrossSectionTrendResetTimers[i] > state.step.input.trendtransg[i].dt)
                state.gasCrossSectionTrendResetTimers[i] = 0;
        }
    }
}

void writeScreenOutput(const TransientSolveState &state, const chrono::steady_clock::time_point &begin, const chrono::steady_clock::time_point &end) {
    if (state.step.input.saidaTela == 1) {
        cout << state.step.stepIndex << "  " << (*state.step.globals).lixo5 << " " << state.step.timeStep;
        for (int i = 0; i < state.step.input.ntela; i++) {
            int posic = state.step.input.tela[i].posic;
            if (state.step.input.tela[i].col == 1) {
                switch (state.step.input.tela[i].var) {
                case 1:
                    cout << " " << state.step.fullModel;
                    break;
                case 2:
                    cout << " " << state.step.cells[posic].temp;
                    break;
                case 3:
                    cout << " " << state.step.cells[posic].alf;
                    break;
                case 4:
                    cout << " " << state.step.cells[posic].bet;
                    break;
                case 5:
                    cout << " " << (state.step.cells[posic].QG / state.step.cells[posic].duto.area);
                    break;
                case 6:
                    cout << " " << (state.step.cells[posic].QL / state.step.cells[posic].duto.area);
                    break;
                }
            } else {
                switch (state.step.input.tela[i].var) {
                case 1:
                    cout << " " << state.step.gasCells[posic].pres;
                    break;
                case 2:
                    cout << " " << state.step.gasCells[posic].temp;
                    break;
                case 3:
                    cout << " " << (state.step.gasCells[posic].VGasR / state.step.gasCells[posic].duto.area);
                    break;
                }
            }
        }
        cout << " " << chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
        cout << endl;
    }
}

void writeEventLog(const TransientSolveState &state, int maxEvento) {
    if (state.logCounter < maxEvento) {
        while (fabs(state.step.input.logevento[state.logCounter].instante - (*state.step.globals).lixo5) < state.step.timeStep) {
            // current date/time based on current system
            time_t now = time(0);
            tm *ltm = localtime(&now); /////////// Taken from https://www.tutorialspoint.com/cplusplus/cpp_date_time.htm
            ostringstream saidaT;
            if (state.step.branchIndex < 0) {
                saidaT << state.logBuffer;
            } else {
                saidaT << "Tramo" << state.step.branchIndex << "-" << state.logBuffer;
            }
            string tmp = saidaT.str();
            ofstream escreveIni(tmp.c_str(), ios_base::app);
            escreveIni << "************************************************************************************************"
                       << endl;
            escreveIni << "Evento Externo = ";
            escreveIni << state.step.input.logevento[state.logCounter].instante << " ; ";
            escreveIni << state.step.input.logevento[state.logCounter].duracao << " ; ";
            escreveIni << state.step.input.logevento[state.logCounter].estIni << " ; ";
            escreveIni << state.step.input.logevento[state.logCounter].estFim << " ; ";
            escreveIni << state.step.input.logevento[state.logCounter].descricao << " ; ";
            escreveIni << "datahora = ";
            escreveIni << ltm->tm_mday << "/";
            escreveIni << 1 + ltm->tm_mon << "/";
            escreveIni << 1900 + ltm->tm_year << " ";
            escreveIni << 0 + ltm->tm_hour << ":";
            escreveIni << 0 + ltm->tm_min << ":";
            escreveIni << 0 + ltm->tm_sec;
            escreveIni << endl;
            state.logCounter++;

            escreveIni.close();
        }
    }
}

void writeProgressReport(const TransientSolveState &state, int MaxKontaImpres) {
    if ((fabs((*state.step.globals).lixo5 * (100. / 5.) / state.step.input.tfinal - round((*state.step.globals).lixo5 * (100. / 5.) / state.step.input.tfinal)) < 0.5 * state.step.timeStep * (100 / 5.) / state.step.input.tfinal) || state.printCounter > MaxKontaImpres || ((*state.step.globals).lixo5 + state.step.timeStep >= state.step.input.tfinal)) {
        if (state.step.input.saidaTela == 0)
            cout << (*state.step.globals).lixo5 * (100.) / state.step.input.tfinal << " % da simulacao alcancado" << endl;
        state.printCounter = 0;
        ostringstream saidaT;
        if (state.step.branchIndex < 0) {
            saidaT << state.logBuffer;
        } else {
            saidaT << "Tramo" << state.step.branchIndex << "-" << state.logBuffer;
        }
        string tmp = saidaT.str();
        ofstream escreveIni(tmp.c_str(), ios_base::app);
        escreveIni << "************************************************************************************************"
                   << endl;
        escreveIni << "Percentual alcancado = " << (*state.step.globals).lixo5 * (100.) / state.step.input.tfinal << " % da simulacao alcancado" << endl;
        escreveIni << "| Passo de Tempo = " << state.step.stepIndex << "| Tempo (s) = " << (*state.step.globals).lixo5 << "| Incremento de Tempo (s) = " << state.step.timeStep
                   << " |" << " Incremento de Tempo Medio CFL (s) = "
                   << state.step.meanCflTimeStep << "| Incremento de Tempo Medio Simulado (s) = " << state.step.meanSimulationTimeStep
                   << " |" << endl;
        for (int i = 0; i < state.step.input.ntela; i++) {
            int posic = state.step.input.tela[i].posic;
            if (state.step.input.tela[i].col == 1) {
                switch (state.step.input.tela[i].var) {
                case 1:
                    escreveIni << " Pressao na Linha de Producao (kgf/cm2), Celula " << posic << " = " << state.step.cells[posic].pres
                               << endl;
                    break;
                case 2:
                    escreveIni << " Temperatura na Linha de Producao (C), Celula " << posic << " = " << state.step.cells[posic].temp
                               << endl;
                    break;
                case 3:
                    escreveIni << " Fracao de Vazio na Linha de Producao (-), Celula " << posic << " = " << state.step.cells[posic].alf
                               << endl;
                    break;
                case 4:
                    escreveIni << " Fracao Beta na Linha de Producao (-), Celula " << posic << " = " << state.step.cells[posic].bet
                               << endl;
                    break;
                case 5:
                    escreveIni << " Velocidade Superficial de Gas na Linha de Producao (m/s), Celula " << posic << " = "
                               << (state.step.cells[posic].QG / state.step.cells[posic].duto.area) << endl;
                    break;
                case 6:
                    escreveIni << " Velocidade Superficial de Liquido na Linha de Producao (m/s), Celula " << posic << " = "
                               << (state.step.cells[posic].QL / state.step.cells[posic].duto.area) << endl;
                    break;
                }
            } else {
                switch (state.step.input.tela[i].var) {
                case 1:
                    escreveIni << " Pressao na Linha de Servico (kgf/cm2), Celula " << posic << " = " << state.step.gasCells[posic].pres
                               << endl;
                    break;
                case 2:
                    escreveIni << " Temperatura na Linha de Servico (C), Ceula " << posic << " = " << state.step.gasCells[posic].temp
                               << endl;
                    break;
                case 3:
                    escreveIni << " Velocidade de Gas na Linha de Servico (m/s), Celula " << posic << " = "
                               << (state.step.gasCells[posic].VGasR / state.step.gasCells[posic].duto.area) << endl;
                    break;
                }
            }
        }
        if (fabs((*state.step.globals).lixo5 - state.step.input.tfinal) <= state.step.timeStep) {
            time_t now = time(0);
            tm *ltm = localtime(&now);
            int diaFim = (ltm->tm_mday);
            int horaFim;
            if (diaFim == diaIni)
                horaFim = ltm->tm_hour;
            else
                horaFim = ltm->tm_hour + 24;
            horaFim *= 3600;
            int minutoFim = 60 * ltm->tm_min;
            int segundoFim = ltm->tm_sec;
            int totalFim = horaFim + minutoFim + segundoFim;
            int totalIni = horaIni * 3600 + minutoIni * 60 + segundoIni;
            escreveIni << "     DURACAO    " << totalFim - totalIni << " segundos " << endl;
            escreveIni << "     Versao    " << versao << endl;
            if (state.step.input.saidaClassica == 1) {
                srand(time(NULL));
                int frase = rand() % 16;
                escreveIni << "*******************************************************************************" << endl;
                escreveIni << "                                  UFA!!!!!!!!                                  " << endl;
                escreveIni << state.closingTitles[frase] << endl;
                escreveIni << state.closingSubtitles[frase] << endl;
                escreveIni << "*******************************************************************************" << endl;
            } else
                escreveIni << "                                 FIM                                  " << endl;
        }
        time_t now = time(0);
        tm *ltm = localtime(&now); /////////// Taken from https://www.tutorialspoint.com/cplusplus/cpp_date_time.htm
        escreveIni << "datahora = ";
        escreveIni << ltm->tm_mday << "/";
        escreveIni << 1 + ltm->tm_mon << "/";
        escreveIni << 1900 + ltm->tm_year << " ";
        escreveIni << 0 + ltm->tm_hour << ":";
        escreveIni << 0 + ltm->tm_min << ":";
        escreveIni << 0 + ltm->tm_sec;
        escreveIni << endl;

        escreveIni.close();
    }
}
}  // namespace


void solveTransientStep(const TransientSolveState &state, double titRev, double alfRev, double betRev, int nrede, ProFlu fluiRev) {
    chrono::steady_clock::time_point begin, end;
    begin = chrono::steady_clock::now();
    double velmaxdesc = 0;

    if ((*state.step.globals).chaverede == 0) {

        state.updaters.solveHydrateEnvelopes();

        if (state.step.input.flashCompleto == 2 && (*state.step.globals).lixo5 < 1e-15 && state.step.input.miniTabAtraso>0) {
            refreshFluidMiniTable(state.step);
        }
        if ((*state.step.globals).lixo5 >= 0) {
            int para;
            para = 0;
           // arq.imprimeProfile(celula, flut, (*vg1dSP).lixo5, indTramo, nrede);
        }

        if ((*state.step.globals).lixo5 < 1e-15) {
        	for(int iCelU=0;iCelU<state.step.input.nCelUnit;iCelU++)state.unitCellTimeCounters[iCelU]=1;
            for (int i = 0; i < state.step.input.ntendp; i++) {
                state.step.input.imprimeTrend(state.step.cells, state.productionTrendMatrix[i], (*state.step.globals).lixo5, i, state.productionTrendCounts[i]);
            }
            state.updaters.updateTemperatures();
        }
        int ciclomax = state.step.input.cicloAcopTerm;

        int vExpli = 0;
        state.step.fullModel = state.step.input.correcaoMassaEspLiq;
        state.step.input.atualizaSonico((*state.step.globals).lixo5, vExpli);
        computeTimeStep(state.step, vExpli);
        state.step.auxiliaryCflTimeStep = state.step.timeStep;
        state.step.finalAuxiliaryTimeStep = state.step.timeStep;

        state.step.restart = 0;
        int celpos = state.step.input.master1.posic;
        // razMast0=celula[celpos].acsr.chk.AreaGarg/celula[celpos].duto.area;//master-only case
        valveOpeningLow(state.step); // several-valve case

        if (state.step.input.controDesc == 1)
            velmaxdesc = state.updaters.searchUnloadingInjectionPressure();
        state.updaters.solveGasLine();
        state.initialGasSurfacePressure = state.step.gasSurfacePressure;
        state.step.input.atualiza(state.startNode, state.step.endNode, state.annulusDrift, state.step.surfaceChoke, state.injectionChoke, state.step.cells, state.step.gasCells, state.step.gasSurfacePressure,
                     state.defaultInletTemperature, state.initialGasPressure, state.initialGasTemperature,
                     state.step.inletPressure, state.step.inletTemperature, state.step.inletQuality, state.step.inletCompletionFraction, (*state.step.globals).lixo5, state.step.timeStep);
        refreshInletCondition(state.step);

        for (int i = 0; i <= state.step.input.nvalv; i++)
            state.step.masterCriticalRatio[i] = 0.5; // several-valve case
        valveOpeningHigh(state.step);            // several-valve case
        // razMast=celula[celpos].acsr.chk.AreaGarg/celula[celpos].duto.area;//master-only case
        for (int i = 0; i <= state.step.input.nvalv; i++)
            if (state.step.masterRatio1[i] != state.step.masterRatio0[i])
                state.step.fullModel = 0; // several-valve case
        TimeStepPolicies::apply<TimeStepHook::AfterValveOpenings>(
            {state.step, kOutsideCouplingLoop, vExpli}); // several-valve case
        if (state.step.fullModel == 0)
            state.step.input.cicloAcopTerm = 0;
        else
            state.step.input.cicloAcopTerm = 1;
        ciclomax = state.step.input.cicloAcopTerm;

        state.step.initiallyOpen = state.step.surfaceChokeOpen;
        state.initialOpenTime = state.step.openTime;
        for (int kontaAcop = 0; kontaAcop <= 1 * state.step.fullModel; kontaAcop++) {
            advanceCouplingIteration(state, kontaAcop, celpos, vExpli, ciclomax, titRev, alfRev, betRev);
        }

        if (state.step.fullModel == 0 || state.step.input.cicloAcopTerm == 0) {
            for (int ciclo = 0; ciclo <= ciclomax; ciclo++) {
                state.updaters.marchTransientEnergy(ciclo, ciclomax);
            }
        }
    }
    if (state.poissonSolver3D.itera > 7) {
        state.poissonSolver3D.penalizaDt = 20;
    }

    for (int i = 1; i <= state.step.lastCell; i++) {
        state.step.cells[i].dTdt = 0.;
        state.step.cells[i].dTdtL = 0.;
        if (state.step.fullModel == 0 || state.step.cells[i].estabCol == 1) {
            state.step.cells[i].dTdtIni = 0.;
            state.step.cells[i].d2pdt2 = 0.;
        }
    }

    if (state.step.input.modoParafina == 1)
        state.updaters.evaluateParaffin();

    storePreviousSources(state.step);
    state.updaters.updateTemperatures();

    state.step.outletPressure = state.step.cells[state.step.lastCell].pres;
    if (state.step.input.lingas > 0 && state.networkCoupled == 1)
        state.updaters.connectTubing();

    if (state.step.input.flashCompleto == 2) {
        for (int i = 0; i < state.step.lastCell; i++) {
            state.step.cells[i].nMolIni = state.step.cells[i].nMol;
        }
    }
    double totbet = 0.;
    for (int i = 0; i < state.step.lastCell; i++)
        totbet += fabs(state.step.cells[i].bet);
    totbet /= state.step.lastCell;
    state.updaters.updateInitialFractions();
    if ((state.trackGasOilRatio > 0 || state.trackGasGravity > 0) && state.step.input.flashCompleto != 2)
        state.updaters.updateGasOilRatioAndCo2(fluiRev);
    if (state.step.input.flashCompleto == 2) {
        state.updaters.updateMolarFractions(fluiRev);
        state.compositionalRefreshCounter++;
        if (state.compositionalRefreshCounter == state.step.input.miniTabAtraso + 1 && state.step.input.miniTabAtraso > 0) {
            generateFluidMiniTables(state.step);
            state.compositionalRefreshCounter = 0;
        }
    } // compositional case
    state.updaters.updateDensities();

    state.temperatureHistory.push_back(state.step.timeStep);
    state.pressureHistory.push_back(state.step.outletPressure * state.step.timeStep);
    double jtemporario = state.step.cells[state.step.lastCell - 1].Mliqini / (state.step.cells[state.step.lastCell - 1].duto.area * ((1. - state.step.cells[state.step.lastCell - 1].bet) * state.step.cells[state.step.lastCell - 1].rpC + state.step.cells[state.step.lastCell - 1].bet * state.step.cells[state.step.lastCell - 1].rcC));
    jtemporario += (state.step.cells[state.step.lastCell - 1].MC - state.step.cells[state.step.lastCell - 1].Mliqini) / (state.step.cells[state.step.lastCell - 1].duto.area * state.step.cells[state.step.lastCell - 1].flui.MasEspGas(state.step.cells[state.step.lastCell - 1].pres, state.step.cells[state.step.lastCell - 1].temp));
    state.fluxHistory.push_back(jtemporario * state.step.timeStep);
    state.voidFractionHistory.push_back(state.step.cells[state.step.lastCell - 1].alf * state.step.timeStep);
    state.movingMeanCounter += state.step.timeStep;
    state.totalPressure += state.step.outletPressure * state.step.timeStep;
    state.totalFlux += jtemporario * state.step.timeStep;
    state.totalVoidFraction += state.step.cells[state.step.lastCell - 1].alf * state.step.timeStep;
    if ((*state.step.globals).lixo5 > state.step.movingMeanTemperature) {
        state.totalPressure -= state.pressureHistory.front();
        state.pressureHistory.erase(state.pressureHistory.begin());
        state.totalFlux -= state.fluxHistory.front();
        state.fluxHistory.erase(state.fluxHistory.begin());
        state.totalVoidFraction -= state.voidFractionHistory.front();
        state.voidFractionHistory.erase(state.voidFractionHistory.begin());
        state.movingMeanCounter -= state.temperatureHistory.front();
        state.temperatureHistory.erase(state.temperatureHistory.begin());
        state.step.movingMeanPressure = state.totalPressure / state.movingMeanCounter;
        state.step.movingMeanFlux = state.totalFlux / state.movingMeanCounter;
        state.movingMeanVoidFraction = state.totalVoidFraction / state.movingMeanCounter;
    }
    // buried pipe
    for (int j = 0; j <= state.step.lastCell; j++) {
        if (state.step.cells[j].calor.difus2D == 1) {
            state.step.cells[j].calor.poisson2D.finalizaPassoTransiente(state.step.timeStep, state.step.branchIndex);
        }
    }

    int MaxKontaImpres = 1000;
    int ordemImpT = 0;
    if ((fabs(state.step.input.logevento[state.logCounter].instante - (*state.step.globals).lixo5) < state.step.timeStep) ||
        (fabs((*state.step.globals).lixo5 * (100. / 5.) / state.step.input.tfinal - round((*state.step.globals).lixo5 * (100. / 5.) / state.step.input.tfinal)) < 0.5 * state.step.timeStep * (100 / 5.) / state.step.input.tfinal) || state.printCounter > MaxKontaImpres || ((*state.step.globals).lixo5 + state.step.timeStep >= state.step.input.tfinal) ||
        (*state.step.globals).lixo5 < 1e-15) {
        state.printPassCount++;
        ordemImpT = 1;
    }

    if ((*state.step.globals).chaverede != 0)
        (*state.step.globals).lixo5 = (*state.step.globals).lixo5R;
    writeProfiles(state, nrede);
    writeTrends(state, ordemImpT, velmaxdesc, nrede);
    // (*vg1dSP).lixo5 += dt; // change 7
    end = chrono::steady_clock::now();
    (*state.step.globals).contador = state.step.stepIndex;
    writeScreenOutput(state, begin, end);

    int maxEvento = state.step.input.logevento.size();
    writeEventLog(state, maxEvento);
    writeProgressReport(state, MaxKontaImpres);

    state.step.stepIndex++;
    state.printCounter++;
    if ((*state.step.globals).chaverede == 0)
        (*state.step.globals).lixo5 += state.step.timeStep;
}

// The per-cell loops run every time step -- the fluid mini-tables and the
// mass sources at the previous time level -- and the helpers only they call.
namespace {

/// Evaluates the fluid at the two minimum-pressure corners of the mini-table,
/// (pmin, tmin) and (pmin, tmax): refreshes fluC's composition at each corner and
/// writes every tabulated property of flui's mini-table there, plus the
/// bubble-point pressure at each temperature.
void fillMiniTableCornersAtMinPressure(const TransientStepState &state, ProFlu &fluC, ProFlu &flui) {
    if (fluC.dCalculatedBeta > 0. && fluC.dCalculatedBeta < 1.)
        fluC.atualizaPropComp(flui.miniTabDin.pmin, flui.miniTabDin.tmin,
                              fluC.dCalculatedBeta, fluC.oCalculatedLiqComposition,
                              fluC.oCalculatedVapComposition, state.input.pocinjec);
    else
        fluC.atualizaPropComp(flui.miniTabDin.pmin, flui.miniTabDin.tmin, -1, NULL, NULL, state.input.pocinjec);
    flui.miniTabDin.rholF[0][0] =
        fluC.MasEspoleo(flui.miniTabDin.pmin, flui.miniTabDin.tmin);
    flui.miniTabDin.rhogF[0][0] =
        fluC.MasEspGas(flui.miniTabDin.pmin, flui.miniTabDin.tmin);
    flui.miniTabDin.DrhogDpF[0][0] =
        fluC.drhodp(flui.miniTabDin.pmin, flui.miniTabDin.tmin);
    flui.miniTabDin.DrholDtF[0][0] =
        fluC.DrholDT(flui.miniTabDin.pmin, flui.miniTabDin.tmin);
    flui.miniTabDin.DrhogDtF[0][0] =
        fluC.drhodt(flui.miniTabDin.pmin, flui.miniTabDin.tmin);
    flui.miniTabDin.valBO[0][0] =
        fluC.BOFunc(flui.miniTabDin.pmin, flui.miniTabDin.tmin);
    flui.miniTabDin.valZ[0][0] =
        fluC.Zdran(flui.miniTabDin.pmin, flui.miniTabDin.tmin);
    flui.miniTabDin.valdZdT[0][0] =
        fluC.DZDT(flui.miniTabDin.pmin, flui.miniTabDin.tmin);
    flui.miniTabDin.valdZdP[0][0] =
        fluC.DZDP(flui.miniTabDin.pmin, flui.miniTabDin.tmin);
    flui.miniTabDin.tit[0][0] =
        fluC.FracMass(flui.miniTabDin.pmin, flui.miniTabDin.tmin);
    flui.miniTabDin.rs[0][0] =
        fluC.RS(flui.miniTabDin.pmin, flui.miniTabDin.tmin);
    flui.miniTabDin.cplF[0][0] =
        fluC.CalorLiq(flui.miniTabDin.pmin, flui.miniTabDin.tmin);
    flui.miniTabDin.cpgF[0][0] =
        fluC.CalorGas(flui.miniTabDin.pmin, flui.miniTabDin.tmin);
    flui.miniTabDin.HlF[0][0] =
        fluC.EntalpLiq(flui.miniTabDin.pmin, flui.miniTabDin.tmin);
    flui.miniTabDin.HgF[0][0] =
        fluC.EntalpGas(flui.miniTabDin.pmin, flui.miniTabDin.tmin);
    flui.miniTabDin.PBF[0] =
        fluC.PB(flui.miniTabDin.pmin, flui.miniTabDin.tmin);
    if (fluC.dCalculatedBeta > 0. && fluC.dCalculatedBeta < 1.)
        fluC.atualizaPropComp(flui.miniTabDin.pmin, flui.miniTabDin.tmax,
                              fluC.dCalculatedBeta, fluC.oCalculatedLiqComposition,
                              fluC.oCalculatedVapComposition, state.input.pocinjec);
    else
        fluC.atualizaPropComp(flui.miniTabDin.pmin, flui.miniTabDin.tmax, -1, NULL, NULL, state.input.pocinjec);
    flui.miniTabDin.rholF[0][1] =
        fluC.MasEspoleo(flui.miniTabDin.pmin, flui.miniTabDin.tmax);
    flui.miniTabDin.rhogF[0][1] =
        fluC.MasEspGas(flui.miniTabDin.pmin, flui.miniTabDin.tmax);
    flui.miniTabDin.DrhogDpF[0][1] =
        fluC.drhodp(flui.miniTabDin.pmin, flui.miniTabDin.tmax);
    flui.miniTabDin.DrholDtF[0][1] =
        fluC.DrholDT(flui.miniTabDin.pmin, flui.miniTabDin.tmax);
    flui.miniTabDin.DrhogDtF[0][1] =
        fluC.drhodt(flui.miniTabDin.pmin, flui.miniTabDin.tmax);
    flui.miniTabDin.valBO[0][1] =
        fluC.BOFunc(flui.miniTabDin.pmin, flui.miniTabDin.tmax);
    flui.miniTabDin.valZ[0][1] =
        fluC.Zdran(flui.miniTabDin.pmin, flui.miniTabDin.tmax);
    flui.miniTabDin.valdZdT[0][1] =
        fluC.DZDT(flui.miniTabDin.pmin, flui.miniTabDin.tmax);
    flui.miniTabDin.valdZdP[0][1] =
        fluC.DZDP(flui.miniTabDin.pmin, flui.miniTabDin.tmax);
    flui.miniTabDin.tit[0][1] =
        fluC.FracMass(flui.miniTabDin.pmin, flui.miniTabDin.tmax);
    flui.miniTabDin.rs[0][1] =
        fluC.RS(flui.miniTabDin.pmin, flui.miniTabDin.tmax);
    flui.miniTabDin.cplF[0][1] =
        fluC.CalorLiq(flui.miniTabDin.pmin, flui.miniTabDin.tmax);
    flui.miniTabDin.cpgF[0][1] =
        fluC.CalorGas(flui.miniTabDin.pmin, flui.miniTabDin.tmax);
    flui.miniTabDin.HlF[0][1] =
        fluC.EntalpLiq(flui.miniTabDin.pmin, flui.miniTabDin.tmax);
    flui.miniTabDin.HgF[0][1] =
        fluC.EntalpGas(flui.miniTabDin.pmin, flui.miniTabDin.tmax);
    flui.miniTabDin.PBF[1] =
        fluC.PB(flui.miniTabDin.pmin, flui.miniTabDin.tmax);
}

/// Evaluates the fluid at the two maximum-pressure corners of the mini-table,
/// (pmax, tmin) and (pmax, tmax): refreshes fluC's composition at each corner and
/// writes every tabulated property of flui's mini-table there.
void fillMiniTableCornersAtMaxPressure(const TransientStepState &state, ProFlu &fluC, ProFlu &flui) {
    if (fluC.dCalculatedBeta > 0. && fluC.dCalculatedBeta < 1.)
        fluC.atualizaPropComp(flui.miniTabDin.pmax, flui.miniTabDin.tmin,
                              fluC.dCalculatedBeta, fluC.oCalculatedLiqComposition,
                              fluC.oCalculatedVapComposition, state.input.pocinjec);
    else
        fluC.atualizaPropComp(flui.miniTabDin.pmax, flui.miniTabDin.tmin, -1, NULL, NULL, state.input.pocinjec);
    flui.miniTabDin.rholF[1][0] =
        fluC.MasEspoleo(flui.miniTabDin.pmax, flui.miniTabDin.tmin);
    flui.miniTabDin.rhogF[1][0] =
        fluC.MasEspGas(flui.miniTabDin.pmax, flui.miniTabDin.tmin);
    flui.miniTabDin.DrhogDpF[1][0] =
        fluC.drhodp(flui.miniTabDin.pmax, flui.miniTabDin.tmin);
    flui.miniTabDin.DrholDtF[1][0] =
        fluC.DrholDT(flui.miniTabDin.pmax, flui.miniTabDin.tmin);
    flui.miniTabDin.DrhogDtF[1][0] =
        fluC.drhodt(flui.miniTabDin.pmax, flui.miniTabDin.tmin);
    flui.miniTabDin.valBO[1][0] =
        fluC.BOFunc(flui.miniTabDin.pmax, flui.miniTabDin.tmin);
    flui.miniTabDin.valZ[1][0] =
        fluC.Zdran(flui.miniTabDin.pmax, flui.miniTabDin.tmin);
    flui.miniTabDin.valdZdT[1][0] =
        fluC.DZDT(flui.miniTabDin.pmax, flui.miniTabDin.tmin);
    flui.miniTabDin.valdZdP[1][0] =
        fluC.DZDP(flui.miniTabDin.pmax, flui.miniTabDin.tmin);
    flui.miniTabDin.tit[1][0] =
        fluC.FracMass(flui.miniTabDin.pmax, flui.miniTabDin.tmin);
    flui.miniTabDin.rs[1][0] =
        fluC.RS(flui.miniTabDin.pmax, flui.miniTabDin.tmin);
    flui.miniTabDin.cplF[1][0] =
        fluC.CalorLiq(flui.miniTabDin.pmax, flui.miniTabDin.tmin);
    flui.miniTabDin.cpgF[1][0] =
        fluC.CalorGas(flui.miniTabDin.pmax, flui.miniTabDin.tmin);
    flui.miniTabDin.HlF[1][0] =
        fluC.EntalpLiq(flui.miniTabDin.pmax, flui.miniTabDin.tmin);
    flui.miniTabDin.HgF[1][0] =
        fluC.EntalpGas(flui.miniTabDin.pmax, flui.miniTabDin.tmin);
    if (fluC.dCalculatedBeta > 0. && fluC.dCalculatedBeta < 1.)
        fluC.atualizaPropComp(flui.miniTabDin.pmax, flui.miniTabDin.tmax,
                              fluC.dCalculatedBeta, fluC.oCalculatedLiqComposition,
                              fluC.oCalculatedVapComposition, state.input.pocinjec);
    else
        fluC.atualizaPropComp(flui.miniTabDin.pmax, flui.miniTabDin.tmax, -1, NULL, NULL, state.input.pocinjec);
    flui.miniTabDin.rholF[1][1] =
        fluC.MasEspoleo(flui.miniTabDin.pmax, flui.miniTabDin.tmax);
    flui.miniTabDin.rhogF[1][1] =
        fluC.MasEspGas(flui.miniTabDin.pmax, flui.miniTabDin.tmax);
    flui.miniTabDin.DrhogDpF[1][1] =
        fluC.drhodp(flui.miniTabDin.pmax, flui.miniTabDin.tmax);
    flui.miniTabDin.DrholDtF[1][1] =
        fluC.DrholDT(flui.miniTabDin.pmax, flui.miniTabDin.tmax);
    flui.miniTabDin.DrhogDtF[1][1] =
        fluC.drhodt(flui.miniTabDin.pmax, flui.miniTabDin.tmax);
    flui.miniTabDin.valBO[1][1] =
        fluC.BOFunc(flui.miniTabDin.pmax, flui.miniTabDin.tmax);
    flui.miniTabDin.valZ[1][1] =
        fluC.Zdran(flui.miniTabDin.pmax, flui.miniTabDin.tmax);
    flui.miniTabDin.valdZdT[1][1] =
        fluC.DZDT(flui.miniTabDin.pmax, flui.miniTabDin.tmax);
    flui.miniTabDin.valdZdP[1][1] =
        fluC.DZDP(flui.miniTabDin.pmax, flui.miniTabDin.tmax);
    flui.miniTabDin.tit[1][1] =
        fluC.FracMass(flui.miniTabDin.pmax, flui.miniTabDin.tmax);
    flui.miniTabDin.rs[1][1] =
        fluC.RS(flui.miniTabDin.pmax, flui.miniTabDin.tmax);
    flui.miniTabDin.cplF[1][1] =
        fluC.CalorLiq(flui.miniTabDin.pmax, flui.miniTabDin.tmax);
    flui.miniTabDin.cpgF[1][1] =
        fluC.CalorGas(flui.miniTabDin.pmax, flui.miniTabDin.tmax);
    flui.miniTabDin.HlF[1][1] =
        fluC.EntalpLiq(flui.miniTabDin.pmax, flui.miniTabDin.tmax);
    flui.miniTabDin.HgF[1][1] =
        fluC.EntalpGas(flui.miniTabDin.pmax, flui.miniTabDin.tmax);
}

}  // namespace

/// Builds a fluid's dynamic mini-table: evaluates it at the four corners
/// (pmin/pmax x tmin/tmax) and reorders the quality corners.
void fillFluidMiniTable(const TransientStepState &state, ProFlu &flui) {
    ProFlu fluC;
    fluC = flui;
    fluC.atualizaPropCompStandard();
    fillMiniTableCornersAtMinPressure(state, fluC, flui);
    fillMiniTableCornersAtMaxPressure(state, fluC, flui);


    std::pair<double, int> titVec[4];

    for(int j=0;j<2;j++){
    	for(int k=0;k<2;k++){
    		titVec[2*j+k]={flui.miniTabDin.tit[j][k],2*j+k};
    	}
    }
    std::sort(titVec, titVec + 4);
    if(titVec[0].first<1e-3){
    	int busca=1;
    	while(busca<4 && titVec[busca].first<1e-3)busca++;
    	if(busca<4){
    		int jtroca;
    		int ktroca;
    		if(titVec[busca].second==0){
    			jtroca=0;
    			ktroca=0;
    		}
    		else if(titVec[busca].second==1){
    			jtroca=0;
    			ktroca=1;
    		}
    		else if(titVec[busca].second==2){
    			jtroca=1;
    			ktroca=0;
    		}
    		else if(titVec[busca].second==3){
    			jtroca=1;
    			ktroca=1;
    		}
    	    for(int j=0;j<2;j++){
    	    	for(int k=0;k<2;k++){
    	    		if(flui.miniTabDin.tit[j][k]<1e-3){
    	    			flui.miniTabDin.rhogF[j][k]=flui.miniTabDin.rhogF[jtroca][ktroca];
    	    			flui.miniTabDin.DrhogDpF[j][k]=flui.miniTabDin.DrhogDpF[jtroca][ktroca];
    	    			flui.miniTabDin.DrhogDtF[j][k]=flui.miniTabDin.DrhogDtF[jtroca][ktroca];
    	    			flui.miniTabDin.valZ[j][k]=flui.miniTabDin.valZ[jtroca][ktroca];
    	    			flui.miniTabDin.valdZdT[j][k]=flui.miniTabDin.valdZdT[jtroca][ktroca];
    	    			flui.miniTabDin.valdZdP[j][k]=flui.miniTabDin.valdZdP[jtroca][ktroca];
    	    			flui.miniTabDin.cpgF[j][k]=flui.miniTabDin.cpgF[jtroca][ktroca];
    	    			flui.miniTabDin.HgF[j][k]=flui.miniTabDin.HgF[jtroca][ktroca];
    	    		}
    	    	}
    	    }
    	}
    }
    if(titVec[3].first>1.-1e-3){
    	int busca=2;
    	while(busca>=0 && titVec[busca].first>1.-1e-3)busca--;
    	if(busca>=0){
    		int jtroca;
    		int ktroca;
    		if(titVec[busca].second==0){
    			jtroca=0;
    			ktroca=0;
    		}
    		else if(titVec[busca].second==1){
    			jtroca=0;
    			ktroca=1;
    		}
    		else if(titVec[busca].second==2){
    			jtroca=1;
    			ktroca=0;
    		}
    		else if(titVec[busca].second==3){
    			jtroca=1;
    			ktroca=1;
    		}
    	    for(int j=0;j<2;j++){
    	    	for(int k=0;k<2;k++){
    	    		if(flui.miniTabDin.tit[j][k]>1.-1e-3){
    	    			flui.miniTabDin.rholF[j][k]=flui.miniTabDin.rholF[jtroca][ktroca];
    	    			flui.miniTabDin.valBO[j][k]=flui.miniTabDin.valBO[jtroca][ktroca];
    	    			flui.miniTabDin.DrholDtF[j][k]=flui.miniTabDin.DrholDtF[jtroca][ktroca];
    	    			flui.miniTabDin.rs[j][k]=flui.miniTabDin.rs[jtroca][ktroca];
    	    			flui.miniTabDin.cplF[j][k]=flui.miniTabDin.cplF[jtroca][ktroca];
    	    			flui.miniTabDin.HlF[j][k]=flui.miniTabDin.HlF[jtroca][ktroca];
    	    		}
    	    	}
    	    }
    	}
    }
}

/// Recentres every cell's mini-table, and the mini-tables of the accessory
/// fluids, on the cell's current pressure and temperature.
void generateFluidMiniTables(const TransientStepState &state) {
    (*state.globals).modoTransiente = 0;
#pragma omp parallel for num_threads((*state.globals).ntrd)
    for (int i = 0; i <= state.lastCell; i++) {
        double delp;
        double delt;
        delp = 0.5 * state.cells[i].pres;
        if (delp > state.input.miniTabDp)
            delp = state.input.miniTabDp;
        if (delp < 5) {
            state.cells[i].flui.miniTabDin.pmax = state.cells[i].pres + 5.;
            state.cells[i].flui.miniTabDin.pmin = state.cells[i].pres - delp;
            if (state.cells[i].flui.miniTabDin.pmin < 0.9)
                state.cells[i].flui.miniTabDin.pmin = 0.9;
        } else {
            state.cells[i].flui.miniTabDin.pmax = state.cells[i].pres + delp;
            state.cells[i].flui.miniTabDin.pmin = state.cells[i].pres - delp;
        }
        delt = state.input.miniTabDt;
        state.cells[i].flui.miniTabDin.tmax = state.cells[i].temp + delt;
        state.cells[i].flui.miniTabDin.tmin = state.cells[i].temp - delt;
        if(state.input.miniTabAtraso > 0)fillFluidMiniTable(state, state.cells[i].flui);
        if (state.cells[i].acsr.tipo == kAccessoryGasInjection) {
            state.cells[i].acsr.injg.FluidoPro.miniTabDin.pmax = state.cells[i].flui.miniTabDin.pmax;
            state.cells[i].acsr.injg.FluidoPro.miniTabDin.pmin = state.cells[i].flui.miniTabDin.pmin;
            state.cells[i].acsr.injg.FluidoPro.miniTabDin.tmax = state.cells[i].flui.miniTabDin.tmax;
            state.cells[i].acsr.injg.FluidoPro.miniTabDin.tmin = state.cells[i].flui.miniTabDin.tmin;
            state.cells[i].acsr.injg.FluidoPro.atualizaPropCompStandard();
            if (state.cells[i].acsr.injg.FluidoPro.dCalculatedBeta > 0. && state.cells[i].acsr.injg.FluidoPro.dCalculatedBeta < 1.)
                state.cells[i].acsr.injg.FluidoPro.atualizaPropComp(state.cells[i].pres, state.cells[i].temp,
                                                               state.cells[i].acsr.injg.FluidoPro.dCalculatedBeta, state.cells[i].acsr.injg.FluidoPro.oCalculatedLiqComposition,
                                                               state.cells[i].acsr.injg.FluidoPro.oCalculatedVapComposition, state.input.pocinjec);
            else
                state.cells[i].acsr.injg.FluidoPro.atualizaPropComp(state.cells[i].pres, state.cells[i].temp, -1, NULL, NULL, state.input.pocinjec);
            if(state.input.miniTabAtraso > 0)fillFluidMiniTable(state, state.cells[i].acsr.injg.FluidoPro);
        } else if (state.cells[i].acsr.tipo == kAccessoryLiquidInjection) {
            state.cells[i].acsr.injl.FluidoPro.miniTabDin.pmax = state.cells[i].flui.miniTabDin.pmax;
            state.cells[i].acsr.injl.FluidoPro.miniTabDin.pmin = state.cells[i].flui.miniTabDin.pmin;
            state.cells[i].acsr.injl.FluidoPro.miniTabDin.tmax = state.cells[i].flui.miniTabDin.tmax;
            state.cells[i].acsr.injl.FluidoPro.miniTabDin.tmin = state.cells[i].flui.miniTabDin.tmin;
            state.cells[i].acsr.injl.FluidoPro.atualizaPropCompStandard();
            if (state.cells[i].acsr.injl.FluidoPro.dCalculatedBeta > 0. && state.cells[i].acsr.injl.FluidoPro.dCalculatedBeta < 1.)
                state.cells[i].acsr.injl.FluidoPro.atualizaPropComp(state.cells[i].pres, state.cells[i].temp,
                                                               state.cells[i].acsr.injl.FluidoPro.dCalculatedBeta, state.cells[i].acsr.injl.FluidoPro.oCalculatedLiqComposition,
                                                               state.cells[i].acsr.injl.FluidoPro.oCalculatedVapComposition, state.input.pocinjec);
            else
                state.cells[i].acsr.injl.FluidoPro.atualizaPropComp(state.cells[i].pres, state.cells[i].temp, -1, NULL, NULL, state.input.pocinjec);
            if(state.input.miniTabAtraso > 0)fillFluidMiniTable(state, state.cells[i].acsr.injl.FluidoPro);
        } else if (state.cells[i].acsr.tipo == kAccessoryInflowPerformance) {
            state.cells[i].acsr.ipr.FluidoPro.miniTabDin.pmax = state.cells[i].flui.miniTabDin.pmax;
            state.cells[i].acsr.ipr.FluidoPro.miniTabDin.pmin = state.cells[i].flui.miniTabDin.pmin;
            state.cells[i].acsr.ipr.FluidoPro.miniTabDin.tmax = state.cells[i].flui.miniTabDin.tmax;
            state.cells[i].acsr.ipr.FluidoPro.miniTabDin.tmin = state.cells[i].flui.miniTabDin.tmin;
            state.cells[i].acsr.ipr.FluidoPro.atualizaPropCompStandard();
            if (state.cells[i].acsr.ipr.FluidoPro.dCalculatedBeta > 0. && state.cells[i].acsr.ipr.FluidoPro.dCalculatedBeta < 1.)
                state.cells[i].acsr.ipr.FluidoPro.atualizaPropComp(state.cells[i].pres, state.cells[i].temp,
                                                              state.cells[i].acsr.ipr.FluidoPro.dCalculatedBeta, state.cells[i].acsr.ipr.FluidoPro.oCalculatedLiqComposition,
                                                              state.cells[i].acsr.ipr.FluidoPro.oCalculatedVapComposition, state.input.pocinjec);
            else
                state.cells[i].acsr.ipr.FluidoPro.atualizaPropComp(state.cells[i].pres, state.cells[i].temp, -1, NULL, NULL, state.input.pocinjec);
            if(state.input.miniTabAtraso > 0)fillFluidMiniTable(state, state.cells[i].acsr.ipr.FluidoPro);
        } else if (state.cells[i].acsr.tipo == kAccessoryRadialPorous) {
        	if(state.input.miniTabAtraso > 0)state.cells[i].acsr.radialPoro.geraMiniTabFlu();
        } else if (state.cells[i].acsr.tipo == kAccessoryPorous2D) {
        	if(state.input.miniTabAtraso > 0)state.cells[i].acsr.poroso2D.geraMiniTabFlu();
        } else if (state.cells[i].acsr.tipo == kAccessoryLeak) {
            state.cells[i].acsr.fontechk.fluidoP.miniTabDin.pmax = state.cells[i].flui.miniTabDin.pmax;
            state.cells[i].acsr.fontechk.fluidoP.miniTabDin.pmin = state.cells[i].flui.miniTabDin.pmin;
            state.cells[i].acsr.fontechk.fluidoP.miniTabDin.tmax = state.cells[i].flui.miniTabDin.tmax;
            state.cells[i].acsr.fontechk.fluidoP.miniTabDin.tmin = state.cells[i].flui.miniTabDin.tmin;
            state.cells[i].acsr.fontechk.fluidoP.atualizaPropCompStandard();
            if (state.cells[i].acsr.fontechk.fluidoP.dCalculatedBeta > 0. && state.cells[i].acsr.fontechk.fluidoP.dCalculatedBeta < 1.)
                state.cells[i].acsr.fontechk.fluidoP.atualizaPropComp(state.cells[i].pres, state.cells[i].temp,
                                                                 state.cells[i].acsr.fontechk.fluidoP.dCalculatedBeta, state.cells[i].acsr.fontechk.fluidoP.oCalculatedLiqComposition,
                                                                 state.cells[i].acsr.fontechk.fluidoP.oCalculatedVapComposition, state.input.pocinjec);
            else
                state.cells[i].acsr.fontechk.fluidoP.atualizaPropComp(state.cells[i].pres, state.cells[i].temp, -1, NULL, NULL, state.input.pocinjec);
            if(state.input.miniTabAtraso > 0)fillFluidMiniTable(state, state.cells[i].acsr.fontechk.fluidoP);
        } else if (state.cells[i].acsr.tipo == kAccessoryMultipleSource) {
            state.cells[i].acsr.injm.FluidoPro.miniTabDin.pmax = state.cells[i].flui.miniTabDin.pmax;
            state.cells[i].acsr.injm.FluidoPro.miniTabDin.pmin = state.cells[i].flui.miniTabDin.pmin;
            state.cells[i].acsr.injm.FluidoPro.miniTabDin.tmax = state.cells[i].flui.miniTabDin.tmax;
            state.cells[i].acsr.injm.FluidoPro.miniTabDin.tmin = state.cells[i].flui.miniTabDin.tmin;
            state.cells[i].acsr.injm.FluidoPro.atualizaPropCompStandard();
            if (state.cells[i].acsr.injm.FluidoPro.dCalculatedBeta > 0. && state.cells[i].acsr.injm.FluidoPro.dCalculatedBeta < 1.)
                state.cells[i].acsr.injm.FluidoPro.atualizaPropComp(state.cells[i].pres, state.cells[i].temp,
                                                               state.cells[i].acsr.injm.FluidoPro.dCalculatedBeta, state.cells[i].acsr.injm.FluidoPro.oCalculatedLiqComposition,
                                                               state.cells[i].acsr.injm.FluidoPro.oCalculatedVapComposition, state.input.pocinjec);
            else
                state.cells[i].acsr.injm.FluidoPro.atualizaPropComp(state.cells[i].pres, state.cells[i].temp, -1, NULL, NULL, state.input.pocinjec);
            if(state.input.miniTabAtraso > 0)fillFluidMiniTable(state, state.cells[i].acsr.injm.FluidoPro);
        }
    }
    (*state.globals).modoTransiente = 1;
}

/// Copies each cell's mass sources to the previous-time-level fields (...ini).
void storePreviousSources(const TransientStepState &state) {
    for (int i = 0; i <= state.lastCell; i++) {
        state.cells[i].fontemassLRini = state.cells[i].fontemassLR;
        state.cells[i].fontemassCRini = state.cells[i].fontemassCR;
        state.cells[i].fontemassGRini = state.cells[i].fontemassGR;
    }
}

}  // namespace sisprod::transient
