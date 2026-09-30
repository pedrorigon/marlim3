#include "SisProdComposition.h"

#include "Leitura.h"
#include "SisProdConstants.h"
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
            double flowArea = state.cells[i].duto.area;
            double superficialMixtureVelocityLeft = (state.cells[i].QL + state.cells[i].QG) / flowArea;
            double superficialMixtureVelocityRight;
            if (i < state.lastCell)
                superficialMixtureVelocityRight = (state.cells[i + 1].QL + state.cells[i + 1].QG) / flowArea;
            else
                superficialMixtureVelocityRight = (state.cells[i].QL + state.cells[i].QG) / flowArea;
            state.cells[i].velPigini = state.cells[i].velPig;
            state.cells[i].velPig = superficialMixtureVelocityLeft * state.cells[i].razPig + superficialMixtureVelocityRight * (1. - state.cells[i].razPig) - state.cells[i].VazaPig / flowArea;
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
        double pressure;
        double temperature;
        pressure = state.cells[i].pres;
        temperature = state.cells[i].temp;
        state.cells[i].rpC = state.cells[i].flui.MasEspLiq(pressure, temperature);
        state.cells[i].rgC = state.cells[i].flui.MasEspGas(pressure, temperature);
        state.cells[i].rcC = state.cells[i].fluicol.MasEspFlu(pressure, temperature);

        state.cells[i].mipC = state.cells[i].flui.ViscOleo(pressure, temperature);
        state.cells[i].migC = state.cells[i].flui.ViscGas(pressure, temperature);
        state.cells[i].micC = state.cells[i].fluicol.VisFlu(pressure, temperature);

        double tmed = state.cells[i - 1].temp;
        if (state.cells[i].VTemper < 0.)
            tmed = state.cells[i].temp;
        ProFlu upwindFluid;
        if (state.cells[i].QL < 0.)
            upwindFluid = state.cells[i].flui;
        else
            upwindFluid = state.cells[i - 1].flui;
        state.cells[i].rpCi = upwindFluid.MasEspLiq(state.cells[i].presaux, tmed);
        state.cells[i].rgCi = upwindFluid.MasEspGas(state.cells[i].presaux, tmed);
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

namespace {

/// What upwinding gives at one face of cell i: the liquid's black-oil properties
/// and the gas's density and composition, each taken from the side its phase
/// flows from.
struct BlackOilFace {
    double rgo;
    double betI;
    double oilVolumeFactor;
    double waterVolumeFactor;
    double bsw;
    double solutionGasRatio;
    double dgO;
    double yco2O;
    double API;
    double BSW;
    double denag;
    double viscL;
    double viscH;
    double razdgd;
    double razdgl;
    /// Written for both faces and read for neither. As a local, the right
    /// face's (rhog1) drew GCC's "set but not used"; a member draws nothing.
    double rhog;
    double rhogST;
    double dgG;
    double yco2G;
};

/// Cell i itself: its liquid holdup, completion fraction and in-situ black-oil
/// properties, and the fluid properties it holds before this step (the *ini
/// members).
struct BlackOilCell {
    double liquidHoldup;
    double completionFraction;
    double rholST;
    double rhog;
    double rhogST;
    double solutionGasRatio;
    double razdgd;
    double razdgl;
    double oilVolumeFactor;
    double waterVolumeFactor;
    double bsw;
    double dgini;
    double yco2ini;
    double rgoini;
    double APIini;
    double BSWini;
    double denagini;
    double viscLini;
    double viscHini;
};

/// The source of cell i, by accessory: the mass rates it brings (dissolvedGas,
/// freeGas, deadOil, water; they were dissolvedGasSource and so on) and the
/// properties of the fluid it brings them in (they were APIF, dgFO, dgFG and so
/// on, and hold cell i's own values when there is no source).
struct BlackOilSource {
    double dissolvedGas;
    double freeGas;
    double deadOil;
    double water;
    double API;
    double BSW;
    double denag;
    double dgO;
    double dgG;
    double yco2O;
    double yco2G;
    double rgo;
    double viscL;
    double viscH;
    double rholST;
    double rhogST;
    double razdgd;
    double razdgl;
    double stockTankQuality;
};

/// The balances of cell i: the face fluxes (Mult*), the volumes at the end of
/// the step (vol*Fim) and the residuals (residuo*).
struct BlackOilBalance {
    double MultOe;
    double MultOd;
    double MultGe;
    double MultGd;
    double volleveFim;
    double residuo;
    double volpesFim;
    double volaguaFim;
    double MultPe;
    double MultPd;
    double residuoP;
    double MultAe;
    double MultAd;
    double residuoA;
};

/// The balances of cell i: its gas-oil ratio, gas density and CO2 fraction, its API,
/// and, with more than one production fluid or in a network, its BSW, water
/// density and dead-oil viscosities, each transported from the faces and the
/// source.
void transportBlackOilBalances(const CompositionState &state, const BlackOilFace &left, const BlackOilFace &right, const BlackOilCell &cell, const BlackOilSource &source, BlackOilBalance &balance, double dx, double flowArea, int i, Vcr<double> &rgo, Vcr<double> &dg, Vcr<double> &yco2, Vcr<double> &API, Vcr<double> &BSW, Vcr<double> &denag, Vcr<double> &VISCL, Vcr<double> &VISCH, double dt) {
    if (state.input.nfluP > 1 || (*state.globals).chaverede != 0) {
        balance.MultPe = 0.;
        if (state.cells[i].flui.dStockTankVaporMassFraction < 1. - 1e-15)
            balance.MultPe = state.cells[i].QL * (1 - left.betI) * (1 - left.bsw) / left.oilVolumeFactor;
        balance.MultPd = 0.;
        if (state.cells[i + 1].flui.dStockTankVaporMassFraction < 1. - 1e-15)
            balance.MultPd = state.cells[i + 1].QL * (1 - right.betI) * (1 - right.bsw) / right.oilVolumeFactor;
        balance.residuoP = (balance.volpesFim - state.cells[i].VolPesaST) * flowArea / dt + (balance.MultPd - balance.MultPe) / dx - source.deadOil / dx;
        balance.MultAe = 0.;
        if (state.cells[i].flui.dStockTankVaporMassFraction < 1. - 1e-15)
            balance.MultAe = state.cells[i].QL * (1 - left.betI) * left.bsw / left.oilVolumeFactor;
        balance.MultAd = 0.;
        if (state.cells[i + 1].flui.dStockTankVaporMassFraction < 1. - 1e-15)
            balance.MultAd = state.cells[i + 1].QL * (1 - right.betI) * right.bsw / right.oilVolumeFactor;
        balance.residuoA = (balance.volaguaFim - state.cells[i].VolAguaST) * flowArea / dt + (balance.MultAd - balance.MultAe) / dx - source.water / dx;
    }
    rgo[i] = (*state.globals).RGOMax;
    if (cell.liquidHoldup > (*state.globals).localtiny && cell.completionFraction < (1. - (*state.globals).localtiny) && cell.bsw < (1. - (*state.globals).localtiny)) {
        rgo[i] = (balance.volleveFim - balance.residuo * dt / flowArea) * cell.oilVolumeFactor / (cell.liquidHoldup * (1 - cell.completionFraction) * (1 - cell.bsw));
        if (rgo[i] > (*state.globals).RGOMax)
            rgo[i] = (*state.globals).RGOMax;
    } else if (cell.completionFraction >= (1. - (*state.globals).localtiny) || cell.bsw >= (1. - (*state.globals).localtiny))
        rgo[i] = 0.;
    else
        rgo[i] = (*state.globals).RGOMax;

    if (balance.volleveFim > 1e-5 && state.input.flashCompleto == 0 && ((source.freeGas >= 0 || source.dissolvedGas > 0) || ((balance.MultGd < 0 || balance.MultGe > 0) || (balance.MultOd < 0 || balance.MultOe > 0)))) {
        dg[i] = (dt * (source.dgO * source.dissolvedGas / dx + source.dgG * source.freeGas / dx + 1. * cell.dgini * balance.residuo - (right.dgO * balance.MultOd - left.dgO * balance.MultOe) / dx - (right.dgG * balance.MultGd - left.dgG * balance.MultGe) / dx) + cell.dgini * state.cells[i].VolLeveST * flowArea) /
                (balance.volleveFim * flowArea - 0. * balance.residuo * dt);
        yco2[i] = (dt * (source.yco2O * source.dissolvedGas / dx + source.yco2G * source.freeGas / dx + 1. * cell.yco2ini * balance.residuo - (right.yco2O * balance.MultOd - left.yco2O * balance.MultOe) / dx - (right.yco2G * balance.MultGd - left.yco2G * balance.MultGe) / dx) + cell.yco2ini * state.cells[i].VolLeveST * flowArea) /
                  (balance.volleveFim * flowArea - 0. * balance.residuo * dt);
        if (yco2[i] < 0.)
            yco2[i] = 0.;
        else if (yco2[i] > 1.)
            yco2[i] = 1.;
    } else {
        dg[i] = cell.dgini;
        yco2[i] = cell.yco2ini;
    }
    if ((state.input.nfluP > 1 && state.input.flashCompleto == 0) || (*state.globals).chaverede != 0) {
        if (balance.volpesFim > 1e-3 && (source.deadOil > 0 || (balance.MultPd < 0 || balance.MultPe > 0))) {
            double denmixSTDF = 141.5 / (131.5 + source.API);
            double denmixSTDini = 141.5 / (131.5 + cell.APIini);
            double denmixSTD1 = 141.5 / (131.5 + right.API);
            double denmixSTD0 = 141.5 / (131.5 + left.API);
            API[i] = (dt * (denmixSTDF * source.deadOil / dx + 1. * denmixSTDini * balance.residuoP - (denmixSTD1 * balance.MultPd - denmixSTD0 * balance.MultPe) / dx) + denmixSTDini * state.cells[i].VolPesaST * flowArea) / (balance.volpesFim * flowArea - 0. * balance.residuoP * dt);
            API[i] = 141.5 / API[i] - 131.5;
            VISCL[i] = (dt * (source.viscL * source.deadOil / dx + 1. * cell.viscLini * balance.residuoP - (right.viscL * balance.MultPd - left.viscL * balance.MultPe) / dx) + cell.viscLini * state.cells[i].VolPesaST * flowArea) / (balance.volpesFim * flowArea - 0. * balance.residuoP * dt);
            VISCH[i] = (dt * (source.viscH * source.deadOil / dx + 1. * cell.viscHini * balance.residuoP - (right.viscH * balance.MultPd - left.viscH * balance.MultPe) / dx) + cell.viscHini * state.cells[i].VolPesaST * flowArea) / (balance.volpesFim * flowArea - 0. * balance.residuoP * dt);
        } else {
            API[i] = cell.APIini;
            VISCL[i] = cell.viscLini;
            VISCH[i] = cell.viscHini;
        }
    }
    if (state.input.nfluP > 1 || (*state.globals).chaverede != 0) {
        if ((balance.volaguaFim + balance.volpesFim) > 1e-3 && ((source.deadOil > 0 || source.water > 0) ||
                                                ((balance.MultPd < 0 || balance.MultPe > 0) || (balance.MultAd < 0 || balance.MultAe > 0)))) {
            BSW[i] = (dt * (source.BSW * (source.water + source.deadOil) / dx + 1. * cell.BSWini * (balance.residuoA + balance.residuoP) - (right.BSW * (balance.MultAd + balance.MultPd) - left.BSW * (balance.MultAe + balance.MultPe)) / dx) + cell.BSWini * (state.cells[i].VolAguaST + state.cells[i].VolPesaST) * flowArea) /
                     ((balance.volaguaFim + balance.volpesFim) * flowArea - 0. * (balance.residuoA + balance.residuoP) * dt);
            denag[i] = (dt * (source.denag * (source.water) / dx + 1. * cell.denagini * (balance.residuoA) - (right.denag * (balance.MultAd)-left.denag * (balance.MultAe)) / dx) + cell.denagini * (state.cells[i].VolAguaST) * flowArea) /
                       ((balance.volaguaFim)*flowArea - 0. * (balance.residuoA)*dt);
            if (BSW[i] < 0.)
                BSW[i] = 0.;
            else if (BSW[i] > 1.)
                BSW[i] = 1.;
            if (denag[i] < 1.)
                denag[i] = 1.;
        } else
            BSW[i] = cell.BSWini;
        denag[i] = cell.denagini;
    }
}

/// No accessory, but mass coming in: the source is the cell's own fluid.
void blackOilSourceCellFluid(const CompositionState &state, BlackOilSource &source, int i, double temperatureHigh, double temperatureLow) {
    if (state.cells[i].acsr.tipo == kAccessoryChoke || state.cells[i].acsr.tipo == kAccessoryVolumetricPump) {
        source.dgO = state.cells[i].flui.Deng;
        source.yco2O = state.cells[i].flui.yco2;
        source.rgo = state.cells[i].flui.RGO;
        source.stockTankQuality = state.cells[i].flui.dStockTankVaporMassFraction;
        double rholiq = (1 - state.cells[i].flui.BSW) * (1000 * 141.5 / (131.5 + state.cells[i].flui.API)) + state.cells[i].flui.BSW * 1000 * state.cells[i].flui.Denag;
        double rhogas = state.cells[i].flui.Deng * kAirDensityAtStandardConditions;
        source.rholST = rholiq + rhogas * source.rgo * (1. - state.cells[i].flui.BSW);
        source.rhogST = rhogas;
        double solutionGasRatioSource = state.cells[i].flui.RS(state.cells[i].pres, state.cells[i].temp) * (kBarrelPerCubicMetre / kCubicFootPerCubicMetre);
        source.razdgd = 1 / state.cells[i].flui.rDgD;
        source.razdgl = 1 / state.cells[i].flui.rDgL;
        source.dissolvedGas = (source.dissolvedGas + source.freeGas) * (source.razdgd * solutionGasRatioSource * (1. - state.cells[i].flui.BSW) / source.rholST);
        if (state.input.nfluP > 1 || (*state.globals).chaverede != 0) {
            double rhoPSTF = (1 - state.cells[i].flui.BSW) * (1000 * 141.5 / (131.5 + state.cells[i].flui.API)) + state.cells[i].flui.BSW * 1000 * state.cells[i].flui.Denag;
            source.deadOil *= ((1 - state.cells[i].flui.BSW) / rhoPSTF);
            source.water *= (state.cells[i].flui.BSW / rhoPSTF);
            source.API = state.cells[i].flui.API;
            source.BSW = state.cells[i].flui.BSW;
            source.viscL = 0 * 30 + 1 * state.cells[i].flui.VisOM(temperatureLow);
            source.viscH = 0 * 20 + 1 * state.cells[i].flui.VisOM(temperatureHigh);
        }
    } else if ((*state.cells[i].acsrL).tipo == kAccessoryChoke || (*state.cells[i].acsrL).tipo == kAccessoryVolumetricPump) {
        double rholiq;
        double rhogas;

        if (i > 0) {
            source.dgO = state.cells[i - 1].flui.Deng;
            source.yco2O = state.cells[i - 1].flui.yco2;
            source.rgo = state.cells[i - 1].flui.RGO;
            rholiq = (1 - state.cells[i - 1].flui.BSW) * (1000 * 141.5 / (131.5 + state.cells[i - 1].flui.API)) +
                     state.cells[i - 1].flui.BSW * 1000 * state.cells[i - 1].flui.Denag;
            rhogas = state.cells[i - 1].flui.Deng * kAirDensityAtStandardConditions;
            source.stockTankQuality = state.cells[i - 1].flui.dStockTankVaporMassFraction;
            source.rholST = rholiq + rhogas * source.rgo * (1. - state.cells[i - 1].flui.BSW);
            source.rhogST = rhogas;
            double solutionGasRatioSource = state.cells[i - 1].flui.RS(state.cells[i].pres, state.cells[i].temp) * (kBarrelPerCubicMetre / kCubicFootPerCubicMetre);
            source.razdgd = 1 / state.cells[i - 1].flui.rDgD;
            source.razdgl = 1 / state.cells[i - 1].flui.rDgL;
            source.dissolvedGas = (source.dissolvedGas + source.freeGas) * (source.razdgd * solutionGasRatioSource *
                                          (1. - state.cells[i - 1].flui.BSW) / source.rholST);
            if (state.input.nfluP > 1 || (*state.globals).chaverede != 0) {
                double rhoPSTF = (1 - state.cells[i - 1].flui.BSW) * (1000 * 141.5 / (131.5 + state.cells[i - 1].flui.API)) +
                                 state.cells[i - 1].flui.BSW * 1000 *
                                     state.cells[i - 1].flui.Denag;
                source.deadOil *= ((1 - state.cells[i - 1].flui.BSW) / rhoPSTF);
                source.water *= (state.cells[i - 1].flui.BSW / rhoPSTF);
                source.API = state.cells[i - 1].flui.API;
                source.BSW = state.cells[i - 1].flui.BSW;
                source.denag = state.cells[i - 1].flui.Denag;
                source.viscL = 0 * 30 + 1 * state.cells[i - 1].flui.VisOM(temperatureLow);
                source.viscH = 0 * 20 + 1 * state.cells[i - 1].flui.VisOM(temperatureHigh);
            }
        } else {
            source.dgO = state.cells[i].flui.Deng;
            source.yco2O = state.cells[i].flui.yco2;
            source.rgo = state.cells[i].flui.RGO;
            rholiq = (1 - state.cells[i].flui.BSW) * (1000 * 141.5 / (131.5 + state.cells[i].flui.API)) +
                     state.cells[i].flui.BSW * 1000 * state.cells[i].flui.Denag;
            rhogas = state.cells[i].flui.Deng * kAirDensityAtStandardConditions;
            source.stockTankQuality = state.cells[i].flui.dStockTankVaporMassFraction;
            source.rholST = rholiq + rhogas * source.rgo * (1. - state.cells[i].flui.BSW);
            source.rhogST = rhogas;
            double solutionGasRatioSource = state.cells[i].flui.RS(state.cells[i].pres, state.cells[i].temp) * (kBarrelPerCubicMetre / kCubicFootPerCubicMetre);
            source.razdgd = 1 / state.cells[i].flui.rDgD;
            source.razdgl = 1 / state.cells[i].flui.rDgL;
            source.dissolvedGas = (source.dissolvedGas + source.freeGas) * (source.razdgd * state.cells[i].flui.RS(state.cells[i].pres, state.cells[i].temp) * (kBarrelPerCubicMetre / kCubicFootPerCubicMetre) *
                                          (1. - state.cells[i].flui.BSW) / source.rholST);
            if (state.input.nfluP > 1 || (*state.globals).chaverede != 0) {
                double rhoPSTF = (1 - state.cells[i].flui.BSW) * (1000 * 141.5 / (131.5 + state.cells[i].flui.API)) +
                                 state.cells[i].flui.BSW * 1000 *
                                     state.cells[i].flui.Denag;
                source.deadOil *= ((1 - state.cells[i].flui.BSW) / rhoPSTF);
                source.water *= (state.cells[i].flui.BSW / rhoPSTF);
                source.API = state.cells[i].flui.API;
                source.BSW = state.cells[i].flui.BSW;
                source.denag = state.cells[i].flui.Denag;
                source.viscL = 0 * 30 + 1 * state.cells[i].flui.VisOM(temperatureLow);
                source.viscH = 0 * 20 + 1 * state.cells[i].flui.VisOM(temperatureHigh);
            }
        }
    }
}

/// A 2D porous-medium source (accessory 16): its transfer fluid.
void blackOilSourcePorous2D(const CompositionState &state, BlackOilSource &source, ProFlu &fluF, int i, double temperatureHigh, double temperatureLow) {
    if ((state.cells[i].fontemassLR + state.cells[i].fontemassGR) > 1e-15)
        fluF = state.cells[i].acsr.radialPoro.flup;
    else
        fluF = state.cells[i].flui;

    source.stockTankQuality = fluF.dStockTankVaporMassFraction;
    source.dgO = fluF.Deng;
    source.yco2O = fluF.yco2;
    source.rgo = fluF.RGO;

    source.rholST = (1 - state.cells[i].acsr.poroso2D.dados.transfer.BSW) * (1000 * 141.5 / (131.5 + fluF.API)) + state.cells[i].acsr.poroso2D.dados.transfer.BSW * 1000 * fluF.Denag + fluF.Deng * kAirDensityAtStandardConditions * source.rgo * (1. - state.cells[i].acsr.poroso2D.dados.transfer.BSW);
    double solutionGasRatioSource = fluF.RS(state.cells[i].pres, state.cells[i].temp) * (kBarrelPerCubicMetre / kCubicFootPerCubicMetre);
    source.razdgd = 1 / fluF.rDgD;
    source.razdgl = 1 / fluF.rDgL;
    source.dissolvedGas = (source.dissolvedGas + source.freeGas) * source.razdgd * (fluF.RS(state.cells[i].pres, state.cells[i].temp) * (kBarrelPerCubicMetre / kCubicFootPerCubicMetre) * (1. - state.cells[i].acsr.poroso2D.dados.transfer.BSW) / source.rholST);
    if (state.input.nfluP > 1 || (*state.globals).chaverede != 0) {
        double bswaux = state.cells[i].acsr.poroso2D.dados.transfer.BSW;
        double contrabsw = 1. - bswaux;
        if (contrabsw <= (*state.globals).localtiny)
            contrabsw = 0.9 * (*state.globals).localtiny;
        double rhoPSTF;
        if (contrabsw > 0)
            rhoPSTF = (1000 * 141.5 / (131.5 + fluF.API)) + (bswaux / contrabsw) * 1000 * fluF.Denag + fluF.Deng * kAirDensityAtStandardConditions * solutionGasRatioSource;
        else
            rhoPSTF = 1000 * fluF.Denag;
        if (contrabsw > (*state.globals).localtiny) {
            source.deadOil *= (1. / rhoPSTF);
            source.water *= ((bswaux / contrabsw) / rhoPSTF);
        } else {
            source.deadOil = 0.;
            source.water *= (1 / (state.cells[i].acsr.poroso2D.dados.transfer.BSW * 1000 * fluF.Denag));
        }
        source.API = fluF.API;
        source.BSW = state.cells[i].acsr.poroso2D.dados.transfer.BSW;
        source.denag = fluF.Denag;
        source.viscL = 0 * 30 + 1 * fluF.VisOM(temperatureLow);
        source.viscH = 0 * 20 + 1 * fluF.VisOM(temperatureHigh);
    }
}

/// A radial porous-medium source (accessory 15): its fluid.
void blackOilSourceRadialPorous(const CompositionState &state, BlackOilSource &source, ProFlu &fluF, int i, double temperatureHigh, double temperatureLow) {
    if ((state.cells[i].fontemassLR + state.cells[i].fontemassGR) > 1e-15)
        fluF = state.cells[i].acsr.radialPoro.flup;
    else
        fluF = state.cells[i].flui;

    source.stockTankQuality = fluF.dStockTankVaporMassFraction;
    source.dgO = fluF.Deng;
    source.yco2O = fluF.yco2;
    source.rgo = fluF.RGO;

    source.rholST = (1 - state.cells[i].acsr.radialPoro.BSW) * (1000 * 141.5 / (131.5 + fluF.API)) + state.cells[i].acsr.radialPoro.BSW * 1000 * fluF.Denag + fluF.Deng * kAirDensityAtStandardConditions * source.rgo * (1. - state.cells[i].acsr.radialPoro.BSW);
    double solutionGasRatioSource = fluF.RS(state.cells[i].pres, state.cells[i].temp) * (kBarrelPerCubicMetre / kCubicFootPerCubicMetre);
    source.razdgd = 1 / fluF.rDgD;
    source.razdgl = 1 / fluF.rDgL;
    source.dissolvedGas = (source.dissolvedGas + source.freeGas) * source.razdgd * (fluF.RS(state.cells[i].pres, state.cells[i].temp) * (kBarrelPerCubicMetre / kCubicFootPerCubicMetre) * (1. - state.cells[i].acsr.radialPoro.BSW) / source.rholST);
    if (state.input.nfluP > 1 || (*state.globals).chaverede != 0) {
        double bswaux = state.cells[i].acsr.radialPoro.BSW;
        double contrabsw = 1. - bswaux;
        if (contrabsw <= (*state.globals).localtiny)
            contrabsw = 0.9 * (*state.globals).localtiny;
        double rhoPSTF;
        if (contrabsw > 0)
            rhoPSTF = (1000 * 141.5 / (131.5 + fluF.API)) + (bswaux / contrabsw) * 1000 * fluF.Denag + fluF.Deng * kAirDensityAtStandardConditions * solutionGasRatioSource;
        else
            rhoPSTF = 1000 * fluF.Denag;
        if (contrabsw > (*state.globals).localtiny) {
            source.deadOil *= (1. / rhoPSTF);
            source.water *= ((bswaux / contrabsw) / rhoPSTF);
        } else {
            source.deadOil = 0.;
            source.water *= (1 / (state.cells[i].acsr.radialPoro.BSW * 1000 * fluF.Denag));
        }
        source.API = fluF.API;
        source.BSW = state.cells[i].acsr.radialPoro.BSW;
        source.denag = fluF.Denag;
        source.viscL = 0 * 30 + 1 * fluF.VisOM(temperatureLow);
        source.viscH = 0 * 20 + 1 * fluF.VisOM(temperatureHigh);
    }
}

/// A leak to or from the annulus (accessory 9): the fluid of the side at the
/// higher pressure. The block declares its own fluF, shadowing the step's.
void blackOilSourceLeak(const CompositionState &state, BlackOilSource &source, int i, double temperatureHigh, double temperatureLow) {
    ProFlu fluF;
    if (state.cells[i].acsr.fontechk.presT > state.cells[i].acsr.fontechk.pamb) {
        fluF = state.cells[i].acsr.fontechk.fluidoP;
    } else {
        fluF = state.cells[i].acsr.fontechk.fluidoPamb;
    }
    source.dgO = fluF.Deng;
    source.yco2O = fluF.yco2;
    source.rgo = fluF.RGO;
    source.stockTankQuality = fluF.dStockTankVaporMassFraction;

    if (fluF.BSW < 1 - (*state.globals).localtiny)
        source.rholST = (1 - fluF.BSW) * (1000 * 141.5 / (131.5 + fluF.API)) + fluF.BSW * 1000 * fluF.Denag + fluF.Deng * kAirDensityAtStandardConditions * source.rgo * (1. - fluF.BSW);
    else
        source.rholST = fluF.BSW * 1000 * fluF.Denag;

    double solutionGasRatioSource = fluF.RS(state.cells[i].pres, state.cells[i].temp) * (kBarrelPerCubicMetre / kCubicFootPerCubicMetre);
    source.razdgd = 1 / fluF.rDgD;
    source.razdgl = 1 / fluF.rDgL;
    if (state.cells[i].acsr.fontechk.ambGas != 1 || (source.dissolvedGas + source.freeGas) < 0.)
        source.dissolvedGas = (source.dissolvedGas + source.freeGas) * source.razdgd * (fluF.RS(state.cells[i].pres, state.cells[i].temp) * (kBarrelPerCubicMetre / kCubicFootPerCubicMetre) * (1. - fluF.BSW) / source.rholST);
    else
        source.dissolvedGas = 0.;
    if (state.input.nfluP > 1 || (*state.globals).chaverede != 0) {
        double bswaux = fluF.BSW;
        double contrabsw = 1. - bswaux;
        if (contrabsw <= (*state.globals).localtiny)
            contrabsw = 0.9 * (*state.globals).localtiny;
        double rhoPSTF;
        if (contrabsw > (*state.globals).localtiny)
            rhoPSTF = (1000 * 141.5 / (131.5 + fluF.API)) + (bswaux / contrabsw) * 1000 * fluF.Denag + fluF.Deng * kAirDensityAtStandardConditions * solutionGasRatioSource;
        else
            rhoPSTF = 1000 * fluF.Denag;
        if (contrabsw > (*state.globals).localtiny) {
            source.deadOil *= (1. / rhoPSTF);
            source.water *= ((bswaux / contrabsw) / rhoPSTF);
        } else {
            source.deadOil = 0.;
            source.water *= (1 / (fluF.BSW * 1000 * fluF.Denag));
        }
        source.API = fluF.API;
        source.BSW = fluF.BSW;
        source.denag = fluF.Denag;
        source.viscL = 0 * 30 + 1 * fluF.VisOM(temperatureLow);
        source.viscH = 0 * 20 + 1 * fluF.VisOM(temperatureHigh);
    }
}

/// A reservoir inflow (accessory 3): the IPR's fluid.
void blackOilSourceInflowPerformance(const CompositionState &state, BlackOilSource &source, ProFlu &fluF, int i, double temperatureHigh, double temperatureLow) {
    if ((state.cells[i].acsr.ipr.Pres) > state.cells[i].pres)
        fluF = state.cells[i].acsr.ipr.FluidoPro;
    else
        fluF = state.cells[i].flui;

    source.stockTankQuality = fluF.dStockTankVaporMassFraction;
    source.dgO = fluF.Deng;
    source.yco2O = fluF.yco2;
    source.rgo = fluF.RGO;

    source.rholST = (1 - fluF.BSW) * (1000 * 141.5 / (131.5 + fluF.API)) + fluF.BSW * 1000 * fluF.Denag + fluF.Deng * kAirDensityAtStandardConditions * source.rgo * (1. - fluF.BSW);
    double solutionGasRatioSource = fluF.RS(state.cells[i].pres, state.cells[i].temp) * (kBarrelPerCubicMetre / kCubicFootPerCubicMetre);
    source.razdgd = 1 / fluF.rDgD;
    source.razdgl = 1 / fluF.rDgL;
    source.dissolvedGas = (source.dissolvedGas + source.freeGas) * source.razdgd * (fluF.RS(state.cells[i].pres, state.cells[i].temp) * (kBarrelPerCubicMetre / kCubicFootPerCubicMetre) * (1. - fluF.BSW) / source.rholST);
    if (state.input.nfluP > 1 || (*state.globals).chaverede != 0) {
        double bswaux = fluF.BSW;
        double contrabsw = 1. - bswaux;
        if (contrabsw <= (*state.globals).localtiny)
            contrabsw = 0.9 * (*state.globals).localtiny;
        double rhoPSTF;
        if (contrabsw > 0)
            rhoPSTF = (1000 * 141.5 / (131.5 + fluF.API)) + (bswaux / contrabsw) * 1000 * fluF.Denag + fluF.Deng * kAirDensityAtStandardConditions * solutionGasRatioSource;
        else
            rhoPSTF = 1000 * fluF.Denag;
        if (contrabsw > (*state.globals).localtiny) {
            source.deadOil *= (1. / rhoPSTF);
            source.water *= ((bswaux / contrabsw) / rhoPSTF);
        } else {
            source.deadOil = 0.;
            source.water *= (1 / (fluF.BSW * 1000 * fluF.Denag));
        }
        source.API = fluF.API;
        source.BSW = fluF.BSW;
        source.denag = fluF.Denag;
        source.viscL = 0 * 30 + 1 * fluF.VisOM(temperatureLow);
        source.viscH = 0 * 20 + 1 * fluF.VisOM(temperatureHigh);
    }
}

/// A multiple source (accessory 10): its fluid.
void blackOilSourceMultipleSource(const CompositionState &state, BlackOilSource &source, ProFlu &fluF, int i, double temperatureHigh, double temperatureLow) {
    if ((state.cells[i].acsr.injm.MassC + state.cells[i].acsr.injm.MassG + state.cells[i].acsr.injm.MassP) > 0.)
        fluF = state.cells[i].acsr.injm.FluidoPro;
    else
        fluF = state.cells[i].flui;

    source.stockTankQuality = fluF.dStockTankVaporMassFraction;
    source.dgO = fluF.Deng;
    source.yco2O = fluF.yco2;
    source.rgo = fluF.RGO;

    double solutionGasRatioSource = fluF.RS(state.cells[i].pres, state.cells[i].temp) * (kBarrelPerCubicMetre / kCubicFootPerCubicMetre);
    source.razdgd = 1 / fluF.rDgD;
    source.razdgl = 1 / fluF.rDgL;
    source.rholST = (1 - fluF.BSW) * (1000 * 141.5 / (131.5 + fluF.API)) +
              fluF.BSW * 1000 *
                  fluF.Denag +
              fluF.Deng * kAirDensityAtStandardConditions * source.rgo * (1. - fluF.BSW);
    source.dissolvedGas = (source.dissolvedGas + source.freeGas) * source.razdgd * (fluF.RS(state.cells[i].pres, state.cells[i].temp) * (kBarrelPerCubicMetre / kCubicFootPerCubicMetre) * (1. - fluF.BSW) / source.rholST);

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
                      fluF.Deng * kAirDensityAtStandardConditions * solutionGasRatioSource;
        else
            rhoPSTF = 1000 * fluF.Denag;
        if (contrabsw > (*state.globals).localtiny) {
            source.deadOil *= (1. / rhoPSTF);
            source.water *= ((bswaux / contrabsw) / rhoPSTF);
        } else {
            source.deadOil = 0.;
            source.water *= (1 / (fluF.BSW * 1000 * fluF.Denag));
        }
        source.API = fluF.API;
        source.BSW = fluF.BSW;
        source.denag = fluF.Denag;
        source.viscL = 0 * 30 + 1 * fluF.VisOM(temperatureLow);
        source.viscH = 0 * 20 + 1 * fluF.VisOM(temperatureHigh);
    }
}

/// A liquid injection (accessory 2): the injected fluid.
void blackOilSourceLiquidInjection(const CompositionState &state, BlackOilSource &source, ProFlu &fluF, int i, double temperatureHigh, double temperatureLow) {
    if (state.cells[i].acsr.injl.QLiq > 0.)
        fluF = state.cells[i].acsr.injl.FluidoPro;
    else
        fluF = state.cells[i].flui;

    source.stockTankQuality = fluF.dStockTankVaporMassFraction;
    source.dgO = fluF.Deng;
    source.yco2O = fluF.yco2;
    source.rgo = fluF.RGO;

    double solutionGasRatioSource = fluF.RS(state.cells[i].pres, state.cells[i].temp) * (kBarrelPerCubicMetre / kCubicFootPerCubicMetre);
    source.razdgd = 1 / fluF.rDgD;
    source.razdgl = 1 / fluF.rDgL;
    if (fluF.BSW < 1 - (*state.globals).localtiny)
        source.rholST = (1 - fluF.BSW) * (1000 * 141.5 / (131.5 + fluF.API)) + fluF.BSW * 1000 * fluF.Denag + fluF.Deng * kAirDensityAtStandardConditions * source.rgo * (1. - fluF.BSW);
    else
        source.rholST = fluF.BSW * 1000 * fluF.Denag;
    source.dissolvedGas = (source.dissolvedGas + source.freeGas) * source.razdgd * (fluF.RS(state.cells[i].pres, state.cells[i].temp) * (kBarrelPerCubicMetre / kCubicFootPerCubicMetre) * (1. - fluF.BSW) / source.rholST);

    if (state.input.nfluP > 1 || (*state.globals).chaverede != 0) {
        double bswaux = fluF.BSW;
        double contrabsw = 1. - bswaux;
        if (contrabsw < (*state.globals).localtiny)
            contrabsw = 0.9 * (*state.globals).localtiny;
        double rhoPSTF;
        if (contrabsw > (*state.globals).localtiny)
            rhoPSTF = (1000 * 141.5 / (131.5 + fluF.API)) + (bswaux / contrabsw) * 1000 * fluF.Denag + fluF.Deng * kAirDensityAtStandardConditions * solutionGasRatioSource;
        else
            rhoPSTF = 1000 * fluF.Denag;
        if (contrabsw > (*state.globals).localtiny) {
            source.deadOil *= (1. / rhoPSTF);
            source.water *= ((bswaux / contrabsw) / rhoPSTF);
        } else {
            source.deadOil = 0.;
            source.water *= (1 / (fluF.BSW * 1000 * fluF.Denag));
        }
        source.API = fluF.API;
        source.BSW = fluF.BSW;
        source.denag = fluF.Denag;
        source.viscL = 0 * 30 + 1 * fluF.VisOM(temperatureLow);
        source.viscH = 0 * 20 + 1 * fluF.VisOM(temperatureHigh);
    }
}

/// A gas injection that carries liquid (accessory 1, not dry): the injected
/// fluid, or the cell's when the injection is idle.
void blackOilSourceWetGasInjection(const CompositionState &state, BlackOilSource &source, ProFlu &fluF, int i, double temperatureHigh, double temperatureLow) {
    if (state.cells[i].acsr.injg.QGas > 0.)
        fluF = state.cells[i].acsr.injg.FluidoPro;
    else
        fluF = state.cells[i].flui;

    source.stockTankQuality = fluF.dStockTankVaporMassFraction;
    source.dgO = fluF.Deng;
    source.yco2O = fluF.yco2;
    source.rgo = fluF.RGO;

    double solutionGasRatioSource = fluF.RS(state.cells[i].pres, state.cells[i].temp) * (kBarrelPerCubicMetre / kCubicFootPerCubicMetre);
    source.razdgd = 1 / fluF.rDgD;
    source.razdgl = 1 / fluF.rDgL;
    if (state.cells[i].acsr.injg.FluidoPro.BSW < 1 - (*state.globals).localtiny)
        source.rholST = (1 - fluF.BSW) * (1000 * 141.5 / (131.5 + fluF.API)) + fluF.BSW * 1000 * fluF.Denag + fluF.Deng * kAirDensityAtStandardConditions * source.rgo * (1. - fluF.BSW);
    else
        source.rholST = fluF.BSW * 1000 * fluF.Denag;
    source.dissolvedGas = (source.dissolvedGas + source.freeGas) * source.razdgd * (fluF.RS(state.cells[i].pres, state.cells[i].temp) * (kBarrelPerCubicMetre / kCubicFootPerCubicMetre) * (1. - fluF.BSW) / source.rholST);

    if (state.input.nfluP > 1 || (*state.globals).chaverede != 0) {
        double bswaux = fluF.BSW;
        double contrabsw = 1. - bswaux;
        if (contrabsw < (*state.globals).localtiny)
            contrabsw = 0.9 * (*state.globals).localtiny;
        double rhoPSTF;
        if (contrabsw > (*state.globals).localtiny)
            rhoPSTF = (1000 * 141.5 / (131.5 + fluF.API)) + (bswaux / contrabsw) * 1000 * fluF.Denag + fluF.Deng * kAirDensityAtStandardConditions * solutionGasRatioSource;
        else
            rhoPSTF = 1000 * fluF.Denag;
        if (contrabsw > (*state.globals).localtiny) {
            source.deadOil *= (1. / rhoPSTF);
            source.water *= ((bswaux / contrabsw) / rhoPSTF);
        } else {
            source.deadOil = 0.;
            source.water *= (1 / (fluF.BSW * 1000 * fluF.Denag));
        }
        source.API = fluF.API;
        source.BSW = fluF.BSW;
        source.denag = state.cells[i].acsr.injg.FluidoPro.Denag;
        source.viscL = 0 * 30 + 1 * fluF.VisOM(temperatureLow);
        source.viscH = 0 * 20 + 1 * fluF.VisOM(temperatureHigh);
    }
}

/// The gas at cell i's left face, from the side it flows from: the upstream cell
/// or the inlet when it flows in, cell i itself otherwise.
void upwindLeftFaceGasProperties(const CompositionState &state, BlackOilFace &left, int i) {
    if ((i > 0 || state.input.ConContEntrada == 1) && state.cells[i].QG > 0) {
        double upstreamPressure;
        double upstreamTemperature;
        if (i > 0) {
            upstreamPressure = state.cells[i - 1].pres;
            upstreamTemperature = state.cells[i - 1].temp;
        } else {
            upstreamPressure = state.inletPressure;
            upstreamTemperature = state.inletTemperature;
        }
        left.rhog = state.cells[i].rgL;
        left.rhogST = (*state.cells[i].fluiL).Deng * kAirDensityAtStandardConditions;
        left.dgG = (*state.cells[i].fluiL).Deng;
        left.yco2G = (*state.cells[i].fluiL).yco2;
    } else {
        left.rhog = state.cells[i].rgC;
        left.rhogST = state.cells[i].flui.Deng * kAirDensityAtStandardConditions;
        left.dgG = state.cells[i].flui.Deng;
        left.yco2G = state.cells[i].flui.yco2;
    }
}

/// The liquid at cell i's right face when it flows back from cell i+1.
void upwindRightFaceLiquidProperties(const CompositionState &state, BlackOilFace &right, int i, double temperatureHigh, double temperatureLow) {
    if (state.cells[i + 1].QL < 0.) {
        right.betI = state.cells[i + 1].betPigE; // testeBeta
        right.rgo = state.cells[i + 1].flui.RGO;
        right.solutionGasRatio = state.cells[i + 1].flui.RS(state.cells[i + 1].pres, state.cells[i + 1].temp);
        right.oilVolumeFactor = state.cells[i + 1].flui.BOFunc(state.cells[i + 1].pres, state.cells[i + 1].temp, right.solutionGasRatio);
        right.waterVolumeFactor = state.cells[i + 1].flui.BAFunc(state.cells[i + 1].pres, state.cells[i + 1].temp);
        right.bsw = state.cells[i + 1].flui.BSW * right.waterVolumeFactor / (right.oilVolumeFactor + right.waterVolumeFactor * state.cells[i + 1].flui.BSW - state.cells[i + 1].flui.BSW * right.oilVolumeFactor);
        right.solutionGasRatio = right.solutionGasRatio * kBarrelPerCubicMetre / kCubicFootPerCubicMetre;
        right.dgO = state.cells[i + 1].flui.Deng;
        right.razdgd = 1 / state.cells[i + 1].flui.rDgD;
        right.razdgl = 1 / state.cells[i + 1].flui.rDgL;
        right.yco2O = state.cells[i + 1].flui.yco2;
        right.API = state.cells[i + 1].flui.API;
        right.BSW = state.cells[i + 1].flui.BSW;
        right.denag = state.cells[i + 1].flui.Denag;
        right.viscL = 0 * 30 + 1 * state.cells[i + 1].flui.VisOM(temperatureLow);
        right.viscH = 0 * 20 + 1 * state.cells[i + 1].flui.VisOM(temperatureHigh);
    }
}

/// The liquid's black-oil properties at cell i's left face, from the side it flows
/// from: the upstream cell or the inlet when it flows in, cell i itself when it
/// flows back.
void upwindLeftFaceBlackOilLiquid(const CompositionState &state, BlackOilFace &left, int i, double temperatureHigh, double temperatureLow) {
    if ((i > 0 || state.input.ConContEntrada == 1) && state.cells[i].QL >= 0.) {
        double upstreamPressure;
        double upstreamTemperature;
        if (i > 0) {
            upstreamPressure = state.cells[i - 1].pres;
            upstreamTemperature = state.cells[i - 1].temp;
        } else {
            upstreamPressure = state.inletPressure;
            upstreamTemperature = state.inletTemperature;
        }
        left.rgo = (*state.cells[i].fluiL).RGO;
        if (state.input.ConContEntrada == 0)
            left.betI = state.cells[i - 1].betPigD; // testeBeta
        else
            left.betI = state.inletCompletionFraction; // testeBeta
        left.solutionGasRatio = (*state.cells[i].fluiL).RS(upstreamPressure, upstreamTemperature);
        left.oilVolumeFactor = (*state.cells[i].fluiL).BOFunc(upstreamPressure, upstreamTemperature, left.solutionGasRatio);
        left.waterVolumeFactor = (*state.cells[i].fluiL).BAFunc(upstreamPressure, upstreamTemperature);
        left.bsw = (*state.cells[i].fluiL).BSW * left.waterVolumeFactor / (left.oilVolumeFactor + left.waterVolumeFactor * (*state.cells[i].fluiL).BSW - (*state.cells[i].fluiL).BSW * left.oilVolumeFactor);
        left.solutionGasRatio = left.solutionGasRatio * kBarrelPerCubicMetre / kCubicFootPerCubicMetre;
        left.dgO = (*state.cells[i].fluiL).Deng;
        left.razdgd = 1 / (*state.cells[i].fluiL).rDgD;
        left.razdgl = 1 / (*state.cells[i].fluiL).rDgL;
        left.yco2O = (*state.cells[i].fluiL).yco2;
        left.API = (*state.cells[i].fluiL).API;
        left.BSW = (*state.cells[i].fluiL).BSW;
        left.denag = (*state.cells[i].fluiL).Denag;
        left.viscL = 0 * 30 + 1 * (*state.cells[i].fluiL).VisOM(temperatureLow);
        left.viscH = 0 * 20 + 1 * (*state.cells[i].fluiL).VisOM(temperatureHigh);
    } else {
        left.betI = state.cells[i].betPigE; // testebeta
        left.rgo = state.cells[i].flui.RGO;
        left.solutionGasRatio = state.cells[i].flui.RS(state.cells[i].pres, state.cells[i].temp);
        left.oilVolumeFactor = state.cells[i].flui.BOFunc(state.cells[i].pres, state.cells[i].temp, left.solutionGasRatio);
        left.waterVolumeFactor = state.cells[i].flui.BAFunc(state.cells[i].pres, state.cells[i].temp);
        left.bsw = state.cells[i].flui.BSW * left.waterVolumeFactor / (left.oilVolumeFactor + left.waterVolumeFactor * state.cells[i].flui.BSW - state.cells[i].flui.BSW * left.oilVolumeFactor);
        left.solutionGasRatio = left.solutionGasRatio * kBarrelPerCubicMetre / kCubicFootPerCubicMetre;
        left.dgO = state.cells[i].flui.Deng;
        left.razdgd = 1 / state.cells[i].flui.rDgD;
        left.razdgl = 1 / state.cells[i].flui.rDgL;
        left.yco2O = state.cells[i].flui.yco2;
        left.API = state.cells[i].flui.API;
        left.BSW = state.cells[i].flui.BSW;
        left.denag = state.cells[i].flui.Denag;
        left.viscL = 0 * 30 + 1 * state.cells[i].flui.VisOM(temperatureLow);
        left.viscH = 0 * 20 + 1 * state.cells[i].flui.VisOM(temperatureHigh);
    }
}

/// One cell's step of the black-oil property transport: the face properties, the
/// source by accessory, and the transported gas-oil ratio, gas density, CO2
/// fraction, API, BSW, water density and dead-oil viscosities.
void transportCellBlackOilProperties(const CompositionState &state, int i, Vcr<double> &rgo, Vcr<double> &dg, Vcr<double> &yco2, Vcr<double> &API, Vcr<double> &BSW, Vcr<double> &denag, Vcr<double> &VISCL, Vcr<double> &VISCH, double temperatureHigh, double temperatureLow, double dt) {
    BlackOilBalance balance;
    double flowArea = state.cells[i].duto.area;
    double dx = state.cells[i].dx;
    double temperatureLeft;
    if (state.cells[i].VTemper < 0.)
        temperatureLeft = state.cells[i].temp;
    else {
        if (i > 0)
            temperatureLeft = state.cells[i - 1].temp;
        else if (state.input.ConContEntrada == 1)
            temperatureLeft = state.inletTemperature;
        else
            temperatureLeft = state.cells[i].temp;
    }
    double temperatureRight = state.cells[i].temp;
    if (state.cells[i + 1].VTemper < 0.)
        temperatureRight = state.cells[i + 1].temp;

    BlackOilFace left;
    if (i > 0 || state.input.ConContEntrada == 0) {
        if (i > 0 && state.cells[i].QG >= 0.)
            left.betI = state.cells[i - 1].betPigD;
        else
            left.betI = state.cells[i].betPigE;
    } else {
        if (state.cells[i].QG >= 0.)
            left.betI = state.inletCompletionFraction;
        else
            left.betI = state.cells[i].betPigE;
    }
    upwindLeftFaceBlackOilLiquid(state, left, i, temperatureHigh, temperatureLow);
    if (left.oilVolumeFactor < 1e-15)
        left.oilVolumeFactor = 1e-15;

    BlackOilFace right;
    right.rgo = state.cells[i].flui.RGO;
    right.betI = state.cells[i].betPigD;
    right.solutionGasRatio = state.cells[i].flui.RS(state.cells[i].pres, state.cells[i].temp);
    right.oilVolumeFactor = state.cells[i].flui.BOFunc(state.cells[i].pres, state.cells[i].temp, right.solutionGasRatio);
    right.waterVolumeFactor = state.cells[i].flui.BAFunc(state.cells[i].pres, state.cells[i].temp);
    right.bsw = state.cells[i].flui.BSW * right.waterVolumeFactor / (right.oilVolumeFactor + right.waterVolumeFactor * state.cells[i].flui.BSW - state.cells[i].flui.BSW * right.oilVolumeFactor);
    right.solutionGasRatio = right.solutionGasRatio * kBarrelPerCubicMetre / kCubicFootPerCubicMetre;
    right.dgO = state.cells[i].flui.Deng;
    right.yco2O = state.cells[i].flui.yco2;
    right.API = state.cells[i].flui.API;
    right.BSW = state.cells[i].flui.BSW;
    right.denag = state.cells[i].flui.Denag;
    right.viscL = 0 * 30 + 1 * state.cells[i].flui.VisOM(temperatureLow);
    right.viscH = 0 * 20 + 1 * state.cells[i].flui.VisOM(temperatureHigh);
    right.razdgd = 1 / state.cells[i].flui.rDgD;
    right.razdgl = 1 / state.cells[i].flui.rDgL;

    // betI1 = celula[i + 1].betPigE;    //duvidabeta
    upwindRightFaceLiquidProperties(state, right, i, temperatureHigh, temperatureLow);
    if (right.oilVolumeFactor < 1e-15)
        right.oilVolumeFactor = 1e-15;

    upwindLeftFaceGasProperties(state, left, i);

    right.rhog = state.cells[i].rgC;
    right.rhogST = state.cells[i].flui.Deng * kAirDensityAtStandardConditions;
    right.dgG = state.cells[i].flui.Deng;
    right.yco2G = state.cells[i].flui.yco2;
    if (state.cells[i + 1].QG <= 0.) {
        right.rhog = state.cells[i].rgR;
        right.rhogST = state.cells[i + 1].flui.Deng * kAirDensityAtStandardConditions;
        right.dgG = state.cells[i + 1].flui.Deng;
        right.yco2G = state.cells[i + 1].flui.yco2;
    }

    if (i == 237) {
        int para;
        para = 0;
    }

    BlackOilCell cell;
    cell.liquidHoldup = 1. - state.cells[i].alf;
    cell.completionFraction = state.cells[i].bet;
    cell.rholST = (1 - state.cells[i].flui.BSW) * (1000 * 141.5 / (131.5 + state.cells[i].flui.API)) + state.cells[i].flui.BSW * 1000 * state.cells[i].flui.Denag;
    cell.rhog = state.cells[i].rgC;
    cell.rhogST = state.cells[i].flui.Deng * kAirDensityAtStandardConditions;
    cell.solutionGasRatio = state.cells[i].flui.RS(state.cells[i].pres, state.cells[i].temp);
    cell.razdgd = 1 / state.cells[i].flui.rDgD;
    cell.razdgl = 1 / state.cells[i].flui.rDgL;
    cell.oilVolumeFactor = state.cells[i].flui.BOFunc(state.cells[i].pres, state.cells[i].temp, cell.solutionGasRatio);
    cell.waterVolumeFactor = state.cells[i].flui.BAFunc(state.cells[i].pres, state.cells[i].temp);
    cell.bsw = state.cells[i].flui.BSW * cell.waterVolumeFactor / (cell.oilVolumeFactor + cell.waterVolumeFactor * state.cells[i].flui.BSW - state.cells[i].flui.BSW * cell.oilVolumeFactor);
    cell.solutionGasRatio = cell.solutionGasRatio * kBarrelPerCubicMetre / kCubicFootPerCubicMetre;

    cell.dgini = state.cells[i].flui.Deng;
    cell.yco2ini = state.cells[i].flui.yco2;
    cell.rgoini = state.cells[i].flui.RGO;
    BlackOilSource source;
    source.dissolvedGas = state.cells[i].fontemassLR;
    source.freeGas = state.cells[i].fontemassGR;
    source.deadOil = state.cells[i].fontemassLR;
    source.water = state.cells[i].fontemassLR;
    cell.APIini = state.cells[i].flui.API;
    cell.BSWini = state.cells[i].flui.BSW;
    cell.denagini = state.cells[i].flui.Denag;
    source.API = cell.APIini;
    source.BSW = cell.BSWini;
    source.denag = cell.denagini;
    source.dgO = cell.dgini;
    source.dgG = cell.dgini;
    source.yco2O = cell.yco2ini;
    source.yco2G = cell.yco2ini;
    source.rgo = cell.rgoini;
    cell.viscLini = 0 * 30 + 1 * state.cells[i].flui.VisOM(temperatureLow);
    cell.viscHini = 0 * 20 + 1 * state.cells[i].flui.VisOM(temperatureHigh);
    source.viscL = cell.viscLini;
    source.viscH = cell.viscHini;
    source.rholST = cell.rholST;
    source.rhogST = cell.rhogST;
    source.razdgd = 1.;
    source.razdgl = 1.;
    source.stockTankQuality = 0.;
    ProFlu fluF;
    if (state.cells[i].acsr.tipo == kAccessoryGasInjection && state.cells[i].acsr.injg.seco == 1) {
        if (state.cells[i].acsr.injg.QGas > 0.)
            fluF = state.cells[i].acsr.injg.FluidoPro;
        else
            fluF = state.cells[i].flui;
        source.dgG = fluF.Deng;
        source.yco2G = fluF.yco2;
        source.rhogST = fluF.Deng * kAirDensityAtStandardConditions;
    } else if (state.cells[i].acsr.tipo == kAccessoryGasInjection && state.cells[i].acsr.injg.seco == 0) {
        blackOilSourceWetGasInjection(state, source, fluF, i, temperatureHigh, temperatureLow);
    } else if (state.cells[i].acsr.tipo == kAccessoryLiquidInjection) {
        blackOilSourceLiquidInjection(state, source, fluF, i, temperatureHigh, temperatureLow);
    } else if (state.cells[i].acsr.tipo == kAccessoryMultipleSource) {
        blackOilSourceMultipleSource(state, source, fluF, i, temperatureHigh, temperatureLow);
    } else if (state.cells[i].acsr.tipo == kAccessoryInflowPerformance) {
        blackOilSourceInflowPerformance(state, source, fluF, i, temperatureHigh, temperatureLow);
    } else if (state.cells[i].acsr.tipo == kAccessoryLeak) {
        blackOilSourceLeak(state, source, i, temperatureHigh, temperatureLow);
    } else if (state.cells[i].acsr.tipo == kAccessoryRadialPorous) {
        blackOilSourceRadialPorous(state, source, fluF, i, temperatureHigh, temperatureLow);
    } else if (state.cells[i].acsr.tipo == kAccessoryPorous2D) {
        blackOilSourcePorous2D(state, source, fluF, i, temperatureHigh, temperatureLow);
    } else if ((fabs(source.dissolvedGas) > (*state.globals).localtiny && state.cells[i].acsr.tipo != kAccessoryLiquidInjection && state.cells[i].acsr.tipo != kAccessoryInflowPerformance &&
                state.cells[i].acsr.tipo != kAccessoryLeak && state.cells[i].acsr.tipo != kAccessoryRadialPorous && state.cells[i].acsr.tipo != kAccessoryPorous2D) ||
               (fabs(source.freeGas) > (*state.globals).localtiny && state.cells[i].acsr.tipo != kAccessoryGasInjection && state.cells[i].acsr.tipo != kAccessoryLiquidInjection && state.cells[i].acsr.tipo != kAccessoryInflowPerformance && state.cells[i].acsr.tipo != kAccessoryLeak && state.cells[i].acsr.tipo != kAccessoryRadialPorous && state.cells[i].acsr.tipo != kAccessoryPorous2D)) {
        blackOilSourceCellFluid(state, source, i, temperatureHigh, temperatureLow);
    }
    source.freeGas *= (source.razdgl / (source.rhogST));

    if (source.stockTankQuality > 1. - 1e-15) {
        source.dissolvedGas = 0.;
        source.deadOil = 0.;
        source.water = 0.;
    }

    balance.MultOe = 0.;
    if (state.cells[i].flui.dStockTankVaporMassFraction < 1. - 1e-15)
        balance.MultOe = state.cells[i].QL * (1 - left.betI) * (1 - left.bsw) * left.razdgd * left.solutionGasRatio / left.oilVolumeFactor;
    balance.MultOd = 0.;
    if (state.cells[i + 1].flui.dStockTankVaporMassFraction < 1. - 1e-15)
        balance.MultOd = state.cells[i + 1].QL * (1 - right.betI) * (1 - right.bsw) * right.razdgd * right.solutionGasRatio / right.oilVolumeFactor;
    balance.MultGe = (state.cells[i].MC - state.cells[i].Mliqini) * left.razdgl / (left.rhogST);
    balance.MultGd = (state.cells[i + 1].MC - state.cells[i + 1].Mliqini) * right.razdgl / (right.rhogST);
    balance.volleveFim = (((1 - cell.liquidHoldup) * cell.rhog * cell.razdgl / (cell.rhogST)) + cell.liquidHoldup * (1 - cell.completionFraction) * (1. - cell.bsw) * cell.solutionGasRatio * cell.razdgd / (cell.oilVolumeFactor));
    if (balance.volleveFim < 1e-15)
        balance.volleveFim = 0.;
    balance.residuo = (balance.volleveFim - state.cells[i].VolLeveST) * flowArea / dt + (balance.MultOd - balance.MultOe) / dx + (balance.MultGd - balance.MultGe) / dx - (source.dissolvedGas / dx + source.freeGas / dx);
    balance.volpesFim = cell.liquidHoldup * (1 - cell.completionFraction) * (1 - cell.bsw) / cell.oilVolumeFactor;
    balance.volaguaFim = cell.liquidHoldup * (1 - cell.completionFraction) * cell.bsw; // nao deveria ser dividido por Bo???????????
    balance.MultPe = 0.;
    balance.MultPd = 0.;
    balance.residuoP = 0.;
    balance.MultAe = 0.;
    balance.MultAd = 0.;
    transportBlackOilBalances(state, left, right, cell, source, balance, dx, flowArea, i, rgo, dg, yco2, API, BSW, denag, VISCL, VISCH, dt);
    state.cells[i].VolLeveST = balance.volleveFim;
    state.cells[i].VolPesaST = balance.volpesFim;
    state.cells[i].VolAguaST = balance.volaguaFim;
}

}  // namespace

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
    double temperatureLow;
    if (state.input.flashCompleto == 1 && (state.input.tabent.tmin - 0) > (*state.globals).localtiny)
        temperatureLow = state.input.tabent.tmin + 0.1;
    else
        temperatureLow = 0.;
    double temperatureHigh;
    if (state.input.flashCompleto == 1 && (70 - state.input.tabent.tmax) > (*state.globals).localtiny)
        temperatureHigh = state.input.tabent.tmax - 0.1;
    else
        temperatureHigh = 70.;
    int imin = 1;
    if ((*state.globals).chaverede != 0 && (state.cells[0].acsr.tipo == kAccessoryMultipleSource || state.input.ConContEntrada > 0))
        imin = 0;
    if (state.cells[0].acsr.tipo == kAccessoryRadialPorous && (state.cells[0].fontemassCR + state.cells[0].fontemassGR + state.cells[0].fontemassLR) > 0.)
        state.cells[0].flui.BSW = state.cells[0].acsr.radialPoro.BSW;
    else if (state.cells[0].acsr.tipo == kAccessoryPorous2D && (state.cells[0].fontemassCR + state.cells[0].fontemassGR + state.cells[0].fontemassLR) > 0.)
        state.cells[0].flui.BSW = state.cells[0].acsr.poroso2D.dados.transfer.BSW;
#pragma omp parallel for num_threads((*state.globals).ntrd)
    for (int i = imin; i < state.lastCell; i++) {
        transportCellBlackOilProperties(state, i, rgo, dg, yco2, API, BSW, denag, VISCL, VISCH, temperatureHigh, temperatureLow, dt);
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
                state.cells[i].flui.TempL = temperatureLow;
                state.cells[i].flui.TempH = temperatureHigh;
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
                state.cells[state.lastCell].flui.TempL = temperatureLow;
                state.cells[state.lastCell].flui.TempH = temperatureHigh;
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

namespace {

/// What upwinding gives at one face of cell i in the phase molar-fraction
/// transport: the liquid's properties, the phase mass flow rates and quality, and
/// the molar weights of the mixture and of each phase.
struct PhaseFace {
    double betI;
    double oilVolumeFactor;
    double waterVolumeFactor;
    double bsw;
    double BSW;
    double denag;
    double viscL;
    double viscH;
    double solutionGasRatio;
    double vazMasLiq;
    double vazMasGas;
    double titV;
    double pesoMol;
    double pesoMolO;
    double pesoMolG;
};

/// Cell i itself: its liquid holdup, completion fraction and in-situ properties,
/// the properties it holds before this step (the *ini members), and its molar
/// weight.
struct PhaseCell {
    double liquidHoldup;
    double completionFraction;
    double solutionGasRatio;
    double oilVolumeFactor;
    double waterVolumeFactor;
    double bsw;
    double BSWini;
    double denagini;
    double viscLini;
    double viscHini;
    double pesoMol;
};

/// Cell i's properties once its molar fractions are solved for. The
/// *Transported locals lost that suffix; betIV, bswV, rhoOVol, rhoWVol and
/// titVol keep their names.
struct PhaseTransported {
    double betIV;
    double solutionGasRatio;
    double oilVolumeFactor;
    double waterVolumeFactor;
    double bswV;
    double rhoOVol;
    double rhoWVol;
    double titVol;
    double tempMol;
};

/// The source of cell i, by accessory: the mass rates it brings and the
/// properties of the fluid it brings them in (they were BSWF, rhoOF,
/// waterCutSource, titF, pesoMolF and so on).
struct PhaseSource {
    double dissolvedGas;
    double freeGas;
    double deadOil;
    double water;
    double BSW;
    double denag;
    double viscL;
    double viscH;
    double rhoW;
    double rhoO;
    double waterCut;
    double waterVolumeFactor;
    double oilVolumeFactor;
    double tit;
    double pesoMol;
};

/// The balances of cell i: the face fluxes (Mult*), the volumes at the end of
/// the step (vol*Fim) and the residuals (residuo*).
struct PhaseBalance {
    double volpesFim;
    double volaguaFim;
    double MultPe;
    double MultPd;
    double residuoP;
    double MultAe;
    double MultAd;
    double residuoA;
};

/// Copies a phase's molar fractions into fracMolFase and renormalises them to
/// sum one when the sum is positive. The smallest fraction is found and zeroed
/// before it is subtracted, so the subtraction takes nothing away.
void normalizePhaseFractions(vector<double>& fracMolFase, double* fracMolOriginal, int npseudo) {
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
}


/// The stock-tank volume balances of cell i, and from them the transported BSW,
/// water density and dead-oil viscosities, when there is more than one production
/// fluid or the line is part of a network.
void transportPhaseVolumesAndViscosities(const CompositionState &state, const PhaseFace &left, const PhaseFace &right, const PhaseCell &cell, const PhaseSource &source, PhaseBalance &balance, double dx, double flowArea, int i, Vcr<double> &BSW, Vcr<double> &denag, Vcr<double> &VISCL, Vcr<double> &VISCH, double dt) {
    if ((*state.globals).lixo5 < 1e-15 && state.input.snaps != 1) {
        state.cells[i].VolPesaST = balance.volpesFim;
        state.cells[i].VolAguaST = balance.volaguaFim;
    }
    if (state.input.nfluP > 1 || (*state.globals).chaverede != 0) {
        balance.MultPe = 0.;
        if (state.cells[i].flui.dStockTankVaporMassFraction < 1. - 1e-15)
            balance.MultPe = state.cells[i].QL * (1 - left.betI) * (1 - left.bsw) / left.oilVolumeFactor;
        balance.MultPd = 0.;
        if (state.cells[i + 1].flui.dStockTankVaporMassFraction < 1. - 1e-15)
            balance.MultPd = state.cells[i + 1].QL * (1 - right.betI) * (1 - right.bsw) / right.oilVolumeFactor;
        balance.residuoP = (balance.volpesFim - state.cells[i].VolPesaST) * flowArea / dt + (balance.MultPd - balance.MultPe) / dx - source.deadOil / dx;
        balance.MultAe = 0.;
        if (state.cells[i].flui.dStockTankVaporMassFraction < 1. - 1e-15)
            balance.MultAe = state.cells[i].QL * (1 - left.betI) * left.bsw / left.oilVolumeFactor;
        balance.MultAd = 0.;
        if (state.cells[i + 1].flui.dStockTankVaporMassFraction < 1. - 1e-15)
            balance.MultAd = state.cells[i + 1].QL * (1 - right.betI) * right.bsw / right.oilVolumeFactor;
        balance.residuoA = (balance.volaguaFim - state.cells[i].VolAguaST) * flowArea / dt + (balance.MultAd - balance.MultAe) / dx - source.water / dx;
    }
    if ((state.input.nfluP > 1) || (*state.globals).chaverede != 0) {
        if (((balance.volpesFim > 1e-3) && (source.deadOil > 0 || (balance.MultPd < 0 || balance.MultPe > 0))) && state.input.nfluP > 0) {
            VISCL[i] = (dt * (source.viscL * source.deadOil / dx + 1. * cell.viscLini * balance.residuoP - (right.viscL * balance.MultPd - left.viscL * balance.MultPe) / dx) + cell.viscLini * state.cells[i].VolPesaST * flowArea) / (balance.volpesFim * flowArea - 0. * balance.residuoP * dt);
            VISCH[i] = (dt * (source.viscH * source.deadOil / dx + 1. * cell.viscHini * balance.residuoP - (right.viscH * balance.MultPd - left.viscH * balance.MultPe) / dx) + cell.viscHini * state.cells[i].VolPesaST * flowArea) / (balance.volpesFim * flowArea - 0. * balance.residuoP * dt);
        } else {
            VISCL[i] = cell.viscLini;
            VISCH[i] = cell.viscHini;
        }
    }
    if (state.input.nfluP > 1 || (*state.globals).chaverede != 0) {
        if ((
                (balance.volaguaFim + balance.volpesFim) > 1e-3) &&
            ((source.deadOil > 0 || source.water > 0) ||
             ((balance.MultPd < 0 || balance.MultPe > 0) || (balance.MultAd < 0 || balance.MultAe > 0)))) {
            BSW[i] = (dt * (source.BSW * (source.water + source.deadOil) / dx + 1. * cell.BSWini * (balance.residuoA + balance.residuoP) - (right.BSW * (balance.MultAd + balance.MultPd) - left.BSW * (balance.MultAe + balance.MultPe)) / dx) + cell.BSWini * (state.cells[i].VolAguaST + state.cells[i].VolPesaST) * flowArea) / ((balance.volaguaFim + balance.volpesFim) * flowArea -
                                                                                                                                                                                                                              0. * dt * (balance.residuoA + balance.residuoP));
            denag[i] = (dt * (source.denag * (source.water) / dx + 1. * cell.denagini * (balance.residuoA) - (right.denag * (balance.MultAd)-left.denag * (balance.MultAe)) / dx) + cell.denagini * (state.cells[i].VolAguaST) * flowArea) / ((balance.volaguaFim)*flowArea -
                                                                                                                                                                            0. * dt * (balance.residuoA));
            if (BSW[i] < 0.)
                BSW[i] = 0.;
            else if (BSW[i] > 1.)
                BSW[i] = 1.;
            else if (isnan(BSW[i]))
                BSW[i] = cell.BSWini;
            if (denag[i] < 1.)
                denag[i] = 1.;
            else if (isnan(denag[i]))
                denag[i] = cell.denagini;
        } else
            BSW[i] = cell.BSWini;
        denag[i] = cell.denagini;
    }
}

/// For the accessory in cell i, turns the source's mass rates into stock-tank volume
/// rates of dead oil and water, and takes the source fluid's BSW, water density and
/// dead-oil viscosities; with no accessory and dissolved gas coming in, the same
/// from the cell's own fluid.
void phaseSourceStandardRatesByAccessory(const CompositionState &state, PhaseSource &source, ProFlu &fluF, int i, double temperatureHigh, double temperatureLow) {
    if (state.cells[i].acsr.tipo == kAccessoryLiquidInjection) {
        double solutionGasRatioSource = state.cells[i].acsr.injl.FluidoPro.RS(state.cells[i].pres, state.cells[i].temp) * (kBarrelPerCubicMetre / kCubicFootPerCubicMetre);
        if (state.input.nfluP > 0 || (*state.globals).chaverede != 0) {
            double bswaux = state.cells[i].acsr.injl.FluidoPro.BSW;
            double contrabsw = 1. - bswaux;
            if (contrabsw < (*state.globals).localtiny)
                contrabsw = 0.9 * (*state.globals).localtiny;
            double rhoPSTF;
            if (contrabsw > (*state.globals).localtiny)
                rhoPSTF = (1000 * 141.5 / (131.5 + state.cells[i].acsr.injl.FluidoPro.API)) + (bswaux / contrabsw) * 1000 * state.cells[i].acsr.injl.FluidoPro.Denag + state.cells[i].acsr.injl.FluidoPro.Deng * kAirDensityAtStandardConditions * solutionGasRatioSource;
            else
                rhoPSTF = 1000 * state.cells[i].acsr.injl.FluidoPro.Denag;
            if (contrabsw > (*state.globals).localtiny) {
                source.deadOil *= (1. / rhoPSTF);
                source.water *= ((bswaux / contrabsw) / rhoPSTF);
            } else {
                source.deadOil = 0.;
                source.water *= (1 / (state.cells[i].acsr.injl.FluidoPro.BSW * 1000 * state.cells[i].acsr.injl.FluidoPro.Denag));
            }
            source.BSW = state.cells[i].acsr.injl.FluidoPro.BSW;
            source.denag = state.cells[i].acsr.injl.FluidoPro.Denag;
            source.viscL = 0 * 30 + 1 * state.cells[i].acsr.injl.FluidoPro.VisOM(temperatureLow);
            source.viscH = 0 * 20 + 1 * state.cells[i].acsr.injl.FluidoPro.VisOM(temperatureHigh);
        }
    } else if (state.cells[i].acsr.tipo == kAccessoryMultipleSource) {
        double solutionGasRatioSource = state.cells[i].acsr.injm.FluidoPro.RS(state.cells[i].pres, state.cells[i].temp) * (kBarrelPerCubicMetre / kCubicFootPerCubicMetre);
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
                          state.cells[i].acsr.injm.FluidoPro.Deng * kAirDensityAtStandardConditions * solutionGasRatioSource;
            else
                rhoPSTF = 1000 * state.cells[i].acsr.injm.FluidoPro.Denag;
            if (contrabsw > (*state.globals).localtiny) {
                source.deadOil *= (1. / rhoPSTF);
                source.water *= ((bswaux / contrabsw) / rhoPSTF);
            } else {
                source.deadOil = 0.;
                source.water *= (1 / (state.cells[i].acsr.injm.FluidoPro.BSW * 1000 * state.cells[i].acsr.injm.FluidoPro.Denag));
            }
            source.BSW = state.cells[i].acsr.injm.FluidoPro.BSW;
            source.denag = state.cells[i].acsr.injl.FluidoPro.Denag;
            source.viscL = 0 * 30 + 1 * state.cells[i].acsr.injm.FluidoPro.VisOM(temperatureLow);
            source.viscH = 0 * 20 + 1 * state.cells[i].acsr.injm.FluidoPro.VisOM(temperatureHigh);
        }
    } else if (state.cells[i].acsr.tipo == kAccessoryInflowPerformance) {
        double solutionGasRatioSource = state.cells[i].acsr.ipr.FluidoPro.RS(state.cells[i].pres, state.cells[i].temp) * (kBarrelPerCubicMetre / kCubicFootPerCubicMetre);
        if (state.input.nfluP > 0 || (*state.globals).chaverede != 0) {
            double bswaux = state.cells[i].acsr.ipr.FluidoPro.BSW;
            double contrabsw = 1. - bswaux;
            if (contrabsw <= (*state.globals).localtiny)
                contrabsw = 0.9 * (*state.globals).localtiny;
            double rhoPSTF;
            if (contrabsw > 0)
                rhoPSTF = (1000 * 141.5 / (131.5 + state.cells[i].acsr.ipr.FluidoPro.API)) + (bswaux / contrabsw) * 1000 * state.cells[i].acsr.ipr.FluidoPro.Denag + state.cells[i].acsr.ipr.FluidoPro.Deng * kAirDensityAtStandardConditions * solutionGasRatioSource;
            else
                rhoPSTF = 1000 * state.cells[i].acsr.ipr.FluidoPro.Denag;
            if (contrabsw > (*state.globals).localtiny) {
                source.deadOil *= (1. / rhoPSTF);
                source.water *= ((bswaux / contrabsw) / rhoPSTF);
            } else {
                source.deadOil = 0.;
                source.water *= (1 / (state.cells[i].acsr.ipr.FluidoPro.BSW * 1000 * state.cells[i].acsr.ipr.FluidoPro.Denag));
            }
            source.BSW = state.cells[i].acsr.ipr.FluidoPro.BSW;
            source.denag = state.cells[i].acsr.ipr.FluidoPro.Denag;
            source.viscL = 0 * 30 + 1 * state.cells[i].acsr.ipr.FluidoPro.VisOM(temperatureLow);
            source.viscH = 0 * 20 + 1 * state.cells[i].acsr.ipr.FluidoPro.VisOM(temperatureHigh);
        }
    } else if (state.cells[i].acsr.tipo == kAccessoryLeak) {

        double solutionGasRatioSource = fluF.RS(state.cells[i].pres, state.cells[i].temp) * (kBarrelPerCubicMetre / kCubicFootPerCubicMetre);
        if (state.input.nfluP > 0 || (*state.globals).chaverede != 0) {
            double bswaux = fluF.BSW;
            double contrabsw = 1. - bswaux;
            if (contrabsw <= (*state.globals).localtiny)
                contrabsw = 0.9 * (*state.globals).localtiny;
            double rhoPSTF;
            if (contrabsw > (*state.globals).localtiny)
                rhoPSTF = (1000 * 141.5 / (131.5 + fluF.API)) + (bswaux / contrabsw) * 1000 * fluF.Denag + fluF.Deng * kAirDensityAtStandardConditions * solutionGasRatioSource;
            else
                rhoPSTF = 1000 * fluF.Denag;
            if (contrabsw > (*state.globals).localtiny) {
                source.deadOil *= (1. / rhoPSTF);
                source.water *= ((bswaux / contrabsw) / rhoPSTF);
            } else {
                source.deadOil = 0.;
                source.water *= (1 / (fluF.BSW * 1000 * fluF.Denag));
            }
            source.BSW = fluF.BSW;
            source.denag = fluF.Denag;
            source.viscL = 0 * 30 + 1 * fluF.VisOM(temperatureLow);
            source.viscH = 0 * 20 + 1 * fluF.VisOM(temperatureHigh);
        }
    } else if (state.cells[i].acsr.tipo == kAccessoryRadialPorous) {
        double solutionGasRatioSource = state.cells[i].acsr.radialPoro.flup.RS(state.cells[i].pres, state.cells[i].temp) * (kBarrelPerCubicMetre / kCubicFootPerCubicMetre);
        if (state.input.nfluP > 0 || (*state.globals).chaverede != 0) {
            double bswaux = state.cells[i].acsr.radialPoro.BSW;
            double contrabsw = 1. - bswaux;
            if (contrabsw <= (*state.globals).localtiny)
                contrabsw = 0.9 * (*state.globals).localtiny;
            double rhoPSTF;
            if (contrabsw > 0)
                rhoPSTF = (1000 * 141.5 / (131.5 + state.cells[i].acsr.radialPoro.flup.API)) + (bswaux / contrabsw) * 1000 * state.cells[i].acsr.radialPoro.flup.Denag + state.cells[i].acsr.radialPoro.flup.Deng * kAirDensityAtStandardConditions * solutionGasRatioSource;
            else
                rhoPSTF = 1000 * state.cells[i].acsr.radialPoro.flup.Denag;
            if (contrabsw > (*state.globals).localtiny) {
                source.deadOil *= (1. / rhoPSTF);
                source.water *= ((bswaux / contrabsw) / rhoPSTF);
            } else {
                source.deadOil = 0.;
                source.water *= (1 / (state.cells[i].acsr.radialPoro.BSW * 1000 * state.cells[i].acsr.radialPoro.flup.Denag));
            }
            source.BSW = state.cells[i].acsr.radialPoro.BSW;
            source.denag = state.cells[i].acsr.radialPoro.flup.Denag;
            source.viscL = 0 * 30 + 1 * state.cells[i].acsr.radialPoro.flup.VisOM(temperatureLow);
            source.viscH = 0 * 20 + 1 * state.cells[i].acsr.radialPoro.flup.VisOM(temperatureHigh);
        }
    } else if (state.cells[i].acsr.tipo == kAccessoryPorous2D) {
        double solutionGasRatioSource = state.cells[i].acsr.poroso2D.dados.flup.RS(state.cells[i].pres, state.cells[i].temp) * (kBarrelPerCubicMetre / kCubicFootPerCubicMetre);
        if (state.input.nfluP > 0 || (*state.globals).chaverede != 0) {
            double bswaux = state.cells[i].acsr.poroso2D.dados.transfer.BSW;
            double contrabsw = 1. - bswaux;
            if (contrabsw <= (*state.globals).localtiny)
                contrabsw = 0.9 * (*state.globals).localtiny;
            double rhoPSTF;
            if (contrabsw > 0)
                rhoPSTF = (1000 * 141.5 / (131.5 + state.cells[i].acsr.poroso2D.dados.flup.API)) + (bswaux / contrabsw) * 1000 * state.cells[i].acsr.poroso2D.dados.flup.Denag + state.cells[i].acsr.poroso2D.dados.flup.Deng * kAirDensityAtStandardConditions * solutionGasRatioSource;
            else
                rhoPSTF = 1000 * state.cells[i].acsr.poroso2D.dados.flup.Denag;
            if (contrabsw > (*state.globals).localtiny) {
                source.deadOil *= (1. / rhoPSTF);
                source.water *= ((bswaux / contrabsw) / rhoPSTF);
            } else {
                source.deadOil = 0.;
                source.water *= (1 / (state.cells[i].acsr.poroso2D.dados.transfer.BSW * 1000 * state.cells[i].acsr.poroso2D.dados.flup.Denag));
            }
            source.BSW = state.cells[i].acsr.poroso2D.dados.transfer.BSW;
            source.denag = state.cells[i].acsr.poroso2D.dados.flup.Denag;
            source.viscL = 0 * 30 + 1 * state.cells[i].acsr.poroso2D.dados.flup.VisOM(temperatureLow);
            source.viscH = 0 * 20 + 1 * state.cells[i].acsr.poroso2D.dados.flup.VisOM(temperatureHigh);
        }
    } else if ((fabs(source.dissolvedGas) > (*state.globals).localtiny && state.cells[i].acsr.tipo != kAccessoryLiquidInjection && state.cells[i].acsr.tipo != kAccessoryInflowPerformance &&
                state.cells[i].acsr.tipo != kAccessoryLeak && state.cells[i].acsr.tipo != kAccessoryRadialPorous && state.cells[i].acsr.tipo != kAccessoryPorous2D) ||
               (fabs(source.freeGas) > (*state.globals).localtiny && state.cells[i].acsr.tipo != kAccessoryGasInjection && state.cells[i].acsr.tipo != kAccessoryLiquidInjection && state.cells[i].acsr.tipo != kAccessoryInflowPerformance && state.cells[i].acsr.tipo != kAccessoryLeak && state.cells[i].acsr.tipo != kAccessoryRadialPorous && state.cells[i].acsr.tipo != kAccessoryPorous2D)) {
        if (state.cells[i].acsr.tipo == kAccessoryChoke || state.cells[i].acsr.tipo == kAccessoryVolumetricPump) {
            double solutionGasRatioSource = state.cells[i].flui.RS(state.cells[i].pres, state.cells[i].temp) * (kBarrelPerCubicMetre / kCubicFootPerCubicMetre);
            if (state.input.nfluP > 0 || (*state.globals).chaverede != 0) {
                double rhoPSTF = (1 - state.cells[i].flui.BSW) * (1000 * 141.5 / (131.5 + state.cells[i].flui.API)) + state.cells[i].flui.BSW * 1000 * state.cells[i].flui.Denag;
                source.deadOil *= ((1 - state.cells[i].flui.BSW) / rhoPSTF);
                source.water *= (state.cells[i].flui.BSW / rhoPSTF);
                source.BSW = state.cells[i].flui.BSW;
                source.denag = fluF.Denag;
                source.viscL = 0 * 30 + 1 * state.cells[i].flui.VisOM(temperatureLow);
                source.viscH = 0 * 20 + 1 * state.cells[i].flui.VisOM(temperatureHigh);
            }
        } else if ((*state.cells[i].acsrL).tipo == kAccessoryChoke || (*state.cells[i].acsrL).tipo == kAccessoryVolumetricPump) {
            if (i > 0) {
                double solutionGasRatioSource = state.cells[i - 1].flui.RS(state.cells[i].pres, state.cells[i].temp) * (kBarrelPerCubicMetre / kCubicFootPerCubicMetre);
                if (state.input.nfluP > 0 || (*state.globals).chaverede != 0) {
                    double rhoPSTF = (1 - state.cells[i - 1].flui.BSW) * (1000 * 141.5 / (131.5 + state.cells[i - 1].flui.API)) +
                                     state.cells[i - 1].flui.BSW * 1000 *
                                         state.cells[i - 1].flui.Denag;
                    source.deadOil *= ((1 - state.cells[i - 1].flui.BSW) / rhoPSTF);
                    source.water *= (state.cells[i - 1].flui.BSW / rhoPSTF);
                    source.BSW = state.cells[i - 1].flui.BSW;
                    source.denag = state.cells[i - 1].flui.Denag;
                    source.viscL = 0 * 30 + 1 * state.cells[i - 1].flui.VisOM(temperatureLow);
                    source.viscH = 0 * 20 + 1 * state.cells[i - 1].flui.VisOM(temperatureHigh);
                }
            } else {
                double solutionGasRatioSource = state.cells[i].flui.RS(state.cells[i].pres, state.cells[i].temp) * (kBarrelPerCubicMetre / kCubicFootPerCubicMetre);
                if (state.input.nfluP > 0 || (*state.globals).chaverede != 0) {
                    double rhoPSTF = (1 - state.cells[i].flui.BSW) * (1000 * 141.5 / (131.5 + state.cells[i].flui.API)) +
                                     state.cells[i].flui.BSW * 1000 *
                                         state.cells[i].flui.Denag;
                    source.deadOil *= ((1 - state.cells[i].flui.BSW) / rhoPSTF);
                    source.water *= (state.cells[i].flui.BSW / rhoPSTF);
                    source.BSW = state.cells[i].flui.BSW;
                    source.denag = state.cells[i - 1].flui.Denag;
                    source.viscL = 0 * 30 + 1 * state.cells[i].flui.VisOM(temperatureLow);
                    source.viscH = 0 * 20 + 1 * state.cells[i].flui.VisOM(temperatureHigh);
                }
            }
        }
    }
}

/// The phase molar fractions of cell i from the liquid and vapour mole balances,
/// corrected limCorrige times, and the cell's fluid refreshed from them.
void solveCellPhaseMolarFractions(const CompositionState &state, const PhaseFace &left, const PhaseFace &right, PhaseCell &cell, PhaseTransported &transported, const PhaseSource &source, int limCorrige, ProFlu &fluF, int i, Vcr<double> &fracMol0, Vcr<double> &fracMol1, double dt, ProFlu *fluC, int ncomp) {
    for (int corrige = 0; corrige < limCorrige; corrige++) {

        transported.betIV = state.cells[i].bet;
        transported.solutionGasRatio = fluC[i].RS(state.cells[i].pres, state.cells[i].temp);
        transported.oilVolumeFactor = fluC[i].BOFunc(state.cells[i].pres, state.cells[i].temp, transported.solutionGasRatio);
        transported.waterVolumeFactor = fluC[i].BAFunc(state.cells[i].pres, state.cells[i].temp);
        transported.bswV = fluC[i].BSW * transported.waterVolumeFactor / (transported.oilVolumeFactor + transported.waterVolumeFactor * state.cells[i].flui.BSW - state.cells[i].flui.BSW * transported.oilVolumeFactor);
        transported.rhoOVol = fluC[i].MasEspoleo(state.cells[i].pres, state.cells[i].temp);
        transported.rhoWVol = fluC[i].MasEspAgua(state.cells[i].pres, state.cells[i].temp);
        transported.titVol = (1 - transported.bswV) * transported.rhoOVol / ((1 - transported.bswV) * transported.rhoOVol + transported.bswV * transported.rhoWVol);

        transported.tempMol = (fluC[i].MasEspLiq(state.cells[i].pres, state.cells[i].temp) * (1. - state.cells[i].alf) * (1. - transported.betIV) * transported.titVol +
                   fluC[i].MasEspGas(state.cells[i].pres, state.cells[i].temp) * state.cells[i].alf) *
                  state.cells[i].duto.area * state.cells[i].dx / cell.pesoMol;

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
                                      ((dissolvedGasSource + freeGasSource) * fluF.fracMol[kfrac] / pesoMolF -
                                       (vazMol1O + vazMol1G) +
                                       (vazMol0O + vazMol0G)) *
                                          dt) /
                                     tempMol;*/


            fluC[i].fracMol[kfrac] = (state.cells[i].nMolIni * state.cells[i].flui.fracMol[kfrac] +
                                      ((source.dissolvedGas + source.freeGas) * fluF.fracMol[kfrac] / source.pesoMol -
                                       (right.vazMasLiq + right.vazMasGas) * fracMol1[kfrac] / right.pesoMol +
                                       (left.vazMasLiq + left.vazMasGas) * fracMol0[kfrac] / left.pesoMol) *
                                          dt) /
                                     transported.tempMol;
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
        cell.pesoMol = 0;
        for (int kfrac = 0; kfrac < ncomp; kfrac++)
            cell.pesoMol += fluC[i].masMol[kfrac] * state.cells[i].flui.fracMol[kfrac];
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
}

/// The fluid a source in cell i brings, with its phase compositions: the
/// accessory's own fluid, flashed at the cell's pressure and temperature, or the
/// cell's when the source is idle, with its water cut, volume factors, densities
/// and quality.
void readPhaseSourceFluid(const CompositionState &state, PhaseSource &source, ProFlu &fluF, int i) {
    if (state.cells[i].acsr.tipo == kAccessoryGasInjection) {
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
        source.waterCut = 0.;
        source.tit = 1.;
    } else if (state.cells[i].acsr.tipo == kAccessoryLiquidInjection) {
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
        source.oilVolumeFactor = fluF.BOFunc(state.cells[i].pres, state.cells[i].temp);
        source.waterVolumeFactor = fluF.BAFunc(state.cells[i].pres, state.cells[i].temp);
        source.waterCut = fluF.BSW * source.waterVolumeFactor / (source.oilVolumeFactor + source.waterVolumeFactor * fluF.BSW - fluF.BSW * source.oilVolumeFactor);
        source.rhoO = fluF.MasEspoleo(state.cells[i].pres, state.cells[i].temp);
        source.rhoW = fluF.MasEspAgua(state.cells[i].pres, state.cells[i].temp);
        source.tit = (1 - source.waterCut) * source.rhoO / ((1 - source.waterCut) * source.rhoO + source.waterCut * source.rhoW);
    } else if (state.cells[i].acsr.tipo == kAccessoryInflowPerformance) {
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
        source.oilVolumeFactor = fluF.BOFunc(state.cells[i].pres, state.cells[i].temp);
        source.waterVolumeFactor = fluF.BAFunc(state.cells[i].pres, state.cells[i].temp);
        source.waterCut = fluF.BSW * source.waterVolumeFactor / (source.oilVolumeFactor + source.waterVolumeFactor * fluF.BSW - fluF.BSW * source.oilVolumeFactor);
        source.rhoO = fluF.MasEspoleo(state.cells[i].pres, state.cells[i].temp);
        source.rhoW = fluF.MasEspAgua(state.cells[i].pres, state.cells[i].temp);
        source.tit = (1 - source.waterCut) * source.rhoO / ((1 - source.waterCut) * source.rhoO + source.waterCut * source.rhoW);
    } else if (state.cells[i].acsr.tipo == kAccessoryMultipleSource) {
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
        source.oilVolumeFactor = fluF.BOFunc(state.cells[i].pres, state.cells[i].temp);
        source.waterVolumeFactor = fluF.BAFunc(state.cells[i].pres, state.cells[i].temp);
        source.waterCut = fluF.BSW * source.waterVolumeFactor / (source.oilVolumeFactor + source.waterVolumeFactor * fluF.BSW - fluF.BSW * source.oilVolumeFactor);
        source.rhoO = fluF.MasEspoleo(state.cells[i].pres, state.cells[i].temp);
        source.rhoW = fluF.MasEspAgua(state.cells[i].pres, state.cells[i].temp);
        source.tit = (1 - source.waterCut) * source.rhoO / ((1 - source.waterCut) * source.rhoO + source.waterCut * source.rhoW);
    } else if (state.cells[i].acsr.tipo == kAccessoryLeak && state.cells[i].acsr.fontechk.abertura > 1e-6) {
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
        source.oilVolumeFactor = fluF.BOFunc(state.cells[i].pres, state.cells[i].temp);
        source.waterVolumeFactor = fluF.BAFunc(state.cells[i].pres, state.cells[i].temp);
        source.waterCut = fluF.BSW * source.waterVolumeFactor / (source.oilVolumeFactor + source.waterVolumeFactor * fluF.BSW - fluF.BSW * source.oilVolumeFactor);
        source.rhoO = fluF.MasEspoleo(state.cells[i].pres, state.cells[i].temp);
        source.rhoW = fluF.MasEspAgua(state.cells[i].pres, state.cells[i].temp);
        source.tit = (1 - source.waterCut) * source.rhoO / ((1 - source.waterCut) * source.rhoO + source.waterCut * source.rhoW);
    } else if (state.cells[i].acsr.tipo == kAccessoryRadialPorous) {
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
        source.oilVolumeFactor = fluF.BOFunc(state.cells[i].pres, state.cells[i].temp);
        source.waterVolumeFactor = fluF.BAFunc(state.cells[i].pres, state.cells[i].temp);
        source.waterCut = fluF.BSW * source.waterVolumeFactor / (source.oilVolumeFactor + source.waterVolumeFactor * fluF.BSW - fluF.BSW * source.oilVolumeFactor);
        source.rhoO = fluF.MasEspoleo(state.cells[i].pres, state.cells[i].temp);
        source.rhoW = fluF.MasEspAgua(state.cells[i].pres, state.cells[i].temp);
        source.tit = (1 - source.waterCut) * source.rhoO / ((1 - source.waterCut) * source.rhoO + source.waterCut * source.rhoW);
    } else if (state.cells[i].acsr.tipo == kAccessoryPorous2D) {
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
        source.oilVolumeFactor = fluF.BOFunc(state.cells[i].pres, state.cells[i].temp);
        source.waterVolumeFactor = fluF.BAFunc(state.cells[i].pres, state.cells[i].temp);
        source.waterCut = fluF.BSW * source.waterVolumeFactor / (source.oilVolumeFactor + source.waterVolumeFactor * fluF.BSW - fluF.BSW * source.oilVolumeFactor);
        source.rhoO = fluF.MasEspoleo(state.cells[i].pres, state.cells[i].temp);
        source.rhoW = fluF.MasEspAgua(state.cells[i].pres, state.cells[i].temp);
        source.tit = (1 - source.waterCut) * source.rhoO / ((1 - source.waterCut) * source.rhoO + source.waterCut * source.rhoW);
    } else {
        fluF = state.cells[i].flui;
        source.tit = 0.;
    }
}

/// The vapour molar fractions and masses at cell i's two faces, each normalised
/// from the side the gas comes from, and the mixture's from the side the mixture
/// comes from.
void upwindVapourAndMixtureFractions(const CompositionState &state, PhaseFace &left, PhaseFace &right, int i, Vcr<double> &fracMol0, Vcr<double> &fracMol1, Vcr<double> &fracMol0G, Vcr<double> &fracMol1G, int ncomp) {
    if ((i > 0 || state.input.ConContEntrada == 1) && (state.cells[i].MC - state.cells[i].Mliqini) >= 0.) {
        double upstreamPressure;
        double upstreamTemperature;
        if (i > 0) {
            upstreamPressure = state.cells[i - 1].pres;
            upstreamTemperature = state.cells[i - 1].temp;
        } else {
            upstreamPressure = state.inletPressure;
            upstreamTemperature = state.inletTemperature;
        }
    	left.pesoMolG = 0;
        double titLocal=(*state.cells[i].fluiL).FracMass(upstreamPressure, upstreamTemperature);
        vector<double> fracMolFase(ncomp);
        normalizePhaseFractions(fracMolFase, (*state.cells[i].fluiL).oCalculatedVapComposition, ncomp);
        for (int kfrac = 0; kfrac < ncomp; kfrac++) {
            fracMol0G[kfrac] = fracMolFase[kfrac];
            left.pesoMolG += (*state.cells[i].fluiL).masMol[kfrac] * fracMolFase[kfrac];
        }
        if(left.pesoMolG<1e-3 || titLocal<1e-3 || state.cells[i].alfL<1e-3){
            left.pesoMolG = 0;
            for (int kfrac = 0; kfrac < ncomp; kfrac++) {
                fracMol0G[kfrac] = (*state.cells[i].fluiL).fracMol[kfrac];
                left.pesoMolG += (*state.cells[i].fluiL).masMol[kfrac] * (*state.cells[i].fluiL).fracMol[kfrac];
            }
        }
    } else {
        left.pesoMolG = 0;
        double titLocal=state.cells[i].flui.FracMass(state.cells[i].pres, state.cells[i].temp);
        vector<double> fracMolFase(ncomp);
        normalizePhaseFractions(fracMolFase, state.cells[i].flui.oCalculatedVapComposition, ncomp);
        for (int kfrac = 0; kfrac < ncomp; kfrac++) {
            fracMol0G[kfrac] = fracMolFase[kfrac];
            left.pesoMolG += state.cells[i].flui.masMol[kfrac] * fracMolFase[kfrac];
        }
        if(left.pesoMolG<1e-3 || titLocal<1e-3 || state.cells[i].alf<1e-3){
            left.pesoMolG = 0;
            for (int kfrac = 0; kfrac < ncomp; kfrac++) {
                fracMol0G[kfrac] = state.cells[i].flui.fracMol[kfrac];
                left.pesoMolG += state.cells[i].flui.masMol[kfrac] * state.cells[i].flui.fracMol[kfrac];
            }
        }
    }
    if ((state.cells[i + 1].MC - state.cells[i + 1].Mliqini) >= 0.) {
        right.pesoMolG = 0;
        double titLocal=state.cells[i].flui.FracMass(state.cells[i].pres, state.cells[i].temp);
        vector<double> fracMolFase(ncomp);
        normalizePhaseFractions(fracMolFase, state.cells[i].flui.oCalculatedVapComposition, ncomp);
        for (int kfrac = 0; kfrac < ncomp; kfrac++) {
            fracMol1G[kfrac] = fracMolFase[kfrac];
            right.pesoMolG += state.cells[i].flui.masMol[kfrac] * fracMolFase[kfrac];
        }
        if(right.pesoMolG<1e-3 || titLocal<1e-3 || state.cells[i].alf<1e-3){
            right.pesoMolG = 0;
            for (int kfrac = 0; kfrac < ncomp; kfrac++) {
                fracMol1G[kfrac] = state.cells[i].flui.fracMol[kfrac];
                right.pesoMolG += state.cells[i].flui.masMol[kfrac] * state.cells[i].flui.fracMol[kfrac];
            }
        }
    } else {
        right.pesoMolG = 0;
        double titLocal=state.cells[i+1].flui.FracMass(state.cells[i+1].pres, state.cells[i+1].temp);
        vector<double> fracMolFase(ncomp);
        normalizePhaseFractions(fracMolFase, state.cells[i + 1].flui.oCalculatedVapComposition, ncomp);
        for (int kfrac = 0; kfrac < ncomp; kfrac++) {
            fracMol1G[kfrac] = fracMolFase[kfrac];
            right.pesoMolG += state.cells[i + 1].flui.masMol[kfrac] * fracMolFase[kfrac];
        }
        if(right.pesoMolG<1e-3 || titLocal<1e-3 || state.cells[i+1].alf<1e-3){
            right.pesoMolG = 0;
            for (int kfrac = 0; kfrac < ncomp; kfrac++) {
                fracMol1G[kfrac] = state.cells[i + 1].flui.fracMol[kfrac];
                right.pesoMolG += state.cells[i + 1].flui.masMol[kfrac] * state.cells[i + 1].flui.fracMol[kfrac];
            }
        }
    }

    if ((i > 0 || state.input.ConContEntrada == 1) && (left.vazMasLiq + left.vazMasGas) >= 0) {
        left.pesoMol = 0;
        for (int kfrac = 0; kfrac < ncomp; kfrac++) {
            fracMol0[kfrac] = (*state.cells[i].fluiL).fracMol[kfrac];
            left.pesoMol += (*state.cells[i].fluiL).masMol[kfrac] * (*state.cells[i].fluiL).fracMol[kfrac];
        }
    } else {
        left.pesoMol = 0;
        for (int kfrac = 0; kfrac < ncomp; kfrac++) {
            fracMol0[kfrac] = state.cells[i].flui.fracMol[kfrac];
            left.pesoMol += state.cells[i].flui.masMol[kfrac] * state.cells[i].flui.fracMol[kfrac];
        }
    }
    if ((right.vazMasLiq + right.vazMasGas) >= 0) {
        right.pesoMol = 0;
        for (int kfrac = 0; kfrac < ncomp; kfrac++) {
            fracMol1[kfrac] = state.cells[i].flui.fracMol[kfrac];
            right.pesoMol += state.cells[i].flui.masMol[kfrac] * state.cells[i].flui.fracMol[kfrac];
        }
    } else {
        right.pesoMol = 0;
        for (int kfrac = 0; kfrac < ncomp; kfrac++) {
            fracMol1[kfrac] = state.cells[i + 1].flui.fracMol[kfrac];
            right.pesoMol += state.cells[i + 1].flui.masMol[kfrac] * state.cells[i + 1].flui.fracMol[kfrac];
        }
    }
}

/// The liquid molar fractions and masses at cell i's two faces, each normalised
/// from the side the liquid comes from, with the liquid's mass rates.
void upwindLiquidFractions(const CompositionState &state, PhaseFace &left, PhaseFace &right, int i, Vcr<double> &fracMol0O, Vcr<double> &fracMol1O, int ncomp) {
    if ((i > 0 || state.input.ConContEntrada == 1) && state.cells[i].Mliqini >= 0.) {
        double upstreamPressure;
        double upstreamTemperature;
        if (i > 0) {
            upstreamPressure = state.cells[i - 1].pres;
            upstreamTemperature = state.cells[i - 1].temp;
        } else {
            upstreamPressure = state.inletPressure;
            upstreamTemperature = state.inletTemperature;
        }
        double waterCutCarried = left.bsw;
        double rhoOV = (*state.cells[i].fluiL).MasEspoleo(upstreamPressure, upstreamTemperature);
        double rhoWV = (*state.cells[i].fluiL).MasEspAgua(upstreamPressure, upstreamTemperature);
        double titLocal=(*state.cells[i].fluiL).FracMass(upstreamPressure, upstreamTemperature);
        left.titV = (1 - waterCutCarried) * rhoOV / ((1 - waterCutCarried) * rhoOV + waterCutCarried * rhoWV);
        left.vazMasLiq -= left.betI * state.cells[i].QL * state.cells[i].fluicol.MasEspFlu(upstreamPressure, upstreamTemperature);
        left.vazMasLiq *= left.titV;
        left.pesoMolO = 0;
        vector<double> fracMolFase(ncomp);
        normalizePhaseFractions(fracMolFase, (*state.cells[i].fluiL).oCalculatedLiqComposition, ncomp);
        for (int kfrac = 0; kfrac < ncomp; kfrac++) {
            fracMol0O[kfrac] = fracMolFase[kfrac];
            left.pesoMolO += (*state.cells[i].fluiL).masMol[kfrac] * fracMolFase[kfrac];
        }
        if(left.pesoMolO<1e-3 || titLocal>1.-1e-3 || state.cells[i].alfL>1.-1e-3){
            left.pesoMolO = 0;
            for (int kfrac = 0; kfrac < ncomp; kfrac++) {
                fracMol0O[kfrac] = (*state.cells[i].fluiL).fracMol[kfrac];
                left.pesoMolO += (*state.cells[i].fluiL).masMol[kfrac] * (*state.cells[i].fluiL).fracMol[kfrac];
            }
        }
    } else {
        double waterCutCarried = left.bsw;
        double rhoOV = state.cells[i].flui.MasEspoleo(state.cells[i].pres, state.cells[i].temp);
        double rhoWV = state.cells[i].flui.MasEspAgua(state.cells[i].pres, state.cells[i].temp);
        double titLocal=state.cells[i].flui.FracMass(state.cells[i].pres, state.cells[i].temp);
        left.titV = (1 - waterCutCarried) * rhoOV / ((1 - waterCutCarried) * rhoOV + waterCutCarried * rhoWV);
        left.vazMasLiq -= left.betI * state.cells[i].QL * state.cells[i].fluicol.MasEspFlu(state.cells[i].pres, state.cells[i].temp);
        left.vazMasLiq *= left.titV;
        left.pesoMolO = 0;
        vector<double> fracMolFase(ncomp);
        normalizePhaseFractions(fracMolFase, state.cells[i].flui.oCalculatedLiqComposition, ncomp);
        for (int kfrac = 0; kfrac < ncomp; kfrac++) {
            fracMol0O[kfrac] = fracMolFase[kfrac];
            left.pesoMolO += state.cells[i].flui.masMol[kfrac] * fracMolFase[kfrac];
        }
        if(left.pesoMolO<1e-3 || titLocal>1.-1e-3 || state.cells[i].alf>1.-1e-3){
            left.pesoMolO = 0;
            for (int kfrac = 0; kfrac < ncomp; kfrac++) {
                fracMol0O[kfrac] = state.cells[i].flui.fracMol[kfrac];
                left.pesoMolO += state.cells[i].flui.masMol[kfrac] * state.cells[i].flui.fracMol[kfrac];
            }
        }
    }
    if (state.cells[i + 1].Mliqini >= 0.) {
        double waterCutCarried = right.bsw;
        double rhoOV = state.cells[i].flui.MasEspoleo(state.cells[i].pres, state.cells[i].temp);
        double rhoWV = state.cells[i].flui.MasEspAgua(state.cells[i].pres, state.cells[i].temp);
        double titLocal=state.cells[i].flui.FracMass(state.cells[i].pres, state.cells[i].temp);
        right.titV = (1 - waterCutCarried) * rhoOV / ((1 - waterCutCarried) * rhoOV + waterCutCarried * rhoWV);
        right.vazMasLiq -= right.betI * state.cells[i].QL * state.cells[i].fluicol.MasEspFlu(state.cells[i].pres, state.cells[i].temp);
        right.vazMasLiq *= right.titV;
        right.pesoMolO = 0;
        vector<double> fracMolFase(ncomp);
        normalizePhaseFractions(fracMolFase, state.cells[i].flui.oCalculatedLiqComposition, ncomp);
        for (int kfrac = 0; kfrac < ncomp; kfrac++) {
            fracMol1O[kfrac] = fracMolFase[kfrac];
            right.pesoMolO += state.cells[i].flui.masMol[kfrac] * fracMolFase[kfrac];
        }
        if(right.pesoMolO<1e-3 || titLocal>1.-1e-3 || state.cells[i].alf>1-1e-3){
            right.pesoMolO = 0;
            for (int kfrac = 0; kfrac < ncomp; kfrac++) {
                fracMol1O[kfrac] = state.cells[i].flui.fracMol[kfrac];
                right.pesoMolO += state.cells[i].flui.masMol[kfrac] * state.cells[i].flui.fracMol[kfrac];
            }
        }
    } else {
        double waterCutCarried = right.bsw;
        double rhoOV = state.cells[i + 1].flui.MasEspoleo(state.cells[i + 1].pres, state.cells[i + 1].temp);
        double rhoWV = state.cells[i + 1].flui.MasEspAgua(state.cells[i + 1].pres, state.cells[i + 1].temp);
        double titLocal=state.cells[i+1].flui.FracMass(state.cells[i+1].pres, state.cells[i+1].temp);
        right.titV = (1 - waterCutCarried) * rhoOV / ((1 - waterCutCarried) * rhoOV + waterCutCarried * rhoWV);
        right.vazMasLiq -= right.betI * state.cells[i + 1].QL * state.cells[i].fluicol.MasEspFlu(state.cells[i + 1].pres, state.cells[i + 1].temp);
        right.vazMasLiq *= right.titV;
        right.pesoMolO = 0;
        vector<double> fracMolFase(ncomp);
        normalizePhaseFractions(fracMolFase, state.cells[i + 1].flui.oCalculatedLiqComposition, ncomp);
        for (int kfrac = 0; kfrac < ncomp; kfrac++) {
            fracMol1O[kfrac] = fracMolFase[kfrac];
            right.pesoMolO += state.cells[i + 1].flui.masMol[kfrac] * fracMolFase[kfrac];
        }
        if(right.pesoMolO<1e-3 || titLocal>1.-1e-3 || state.cells[i+1].alf>1-1e-3){
            right.pesoMolO = 0;
            for (int kfrac = 0; kfrac < ncomp; kfrac++) {
                fracMol1O[kfrac] = state.cells[i + 1].flui.fracMol[kfrac];
                right.pesoMolO += state.cells[i + 1].flui.masMol[kfrac] * state.cells[i + 1].flui.fracMol[kfrac];
            }
        }
    }
}

/// The liquid's properties at cell i's left face, water density included, from the
/// side the liquid comes from: the upstream cell or the inlet when it flows in,
/// cell i itself when it flows back.
void upwindLeftFacePhaseLiquidProperties(const CompositionState &state, PhaseFace &left, int i, double temperatureHigh, double temperatureLow) {
    if ((i > 0 || state.input.ConContEntrada == 1) && state.cells[i].QL >= 0.) {
        double upstreamPressure;
        double upstreamTemperature;
        if (i > 0) {
            upstreamPressure = state.cells[i - 1].pres;
            upstreamTemperature = state.cells[i - 1].temp;
        } else {
            upstreamPressure = state.inletPressure;
            upstreamTemperature = state.inletTemperature;
        }
        if (state.input.ConContEntrada == 0)
            left.betI = state.cells[i - 1].betPigD; // testeBeta
        else
            left.betI = state.inletCompletionFraction; // testeBeta
        double solutionGasRatioLeft = (*state.cells[i].fluiL).RS(upstreamPressure, upstreamTemperature);
        left.oilVolumeFactor = (*state.cells[i].fluiL).BOFunc(upstreamPressure, upstreamTemperature, solutionGasRatioLeft);
        left.waterVolumeFactor = (*state.cells[i].fluiL).BAFunc(upstreamPressure, upstreamTemperature);
        left.bsw = (*state.cells[i].fluiL).BSW * left.waterVolumeFactor / (left.oilVolumeFactor + left.waterVolumeFactor * (*state.cells[i].fluiL).BSW - (*state.cells[i].fluiL).BSW * left.oilVolumeFactor);
        left.BSW = (*state.cells[i].fluiL).BSW;
        left.denag = (*state.cells[i].fluiL).Denag;
        left.viscL = 0 * 30 + 1 * (*state.cells[i].fluiL).VisOM(temperatureLow);
        left.viscH = 0 * 20 + 1 * (*state.cells[i].fluiL).VisOM(temperatureHigh);
    } else {
        left.betI = state.cells[i].betPigE; // testebeta
        double solutionGasRatioLeft = state.cells[i].flui.RS(state.cells[i].pres, state.cells[i].temp);
        left.oilVolumeFactor = state.cells[i].flui.BOFunc(state.cells[i].pres, state.cells[i].temp, solutionGasRatioLeft);
        left.waterVolumeFactor = state.cells[i].flui.BAFunc(state.cells[i].pres, state.cells[i].temp);
        left.bsw = state.cells[i].flui.BSW * left.waterVolumeFactor / (left.oilVolumeFactor + left.waterVolumeFactor * state.cells[i].flui.BSW - state.cells[i].flui.BSW * left.oilVolumeFactor);
        left.BSW = state.cells[i].flui.BSW;
        left.denag = state.cells[i].flui.Denag;
        left.viscL = 0 * 30 + 1 * state.cells[i].flui.VisOM(temperatureLow);
        left.viscH = 0 * 20 + 1 * state.cells[i].flui.VisOM(temperatureHigh);
    }
}

/// One cell's step of the phase molar-fraction transport: the face properties,
/// the upwinded phase fractions, the source, the cell's balance and the transported
/// BSW, water density and viscosities.
void transportCellPhaseMolarFractions(const CompositionState &state, int i, Vcr<double> &BSW, Vcr<double> &denag, Vcr<double> &VISCL, Vcr<double> &VISCH, Vcr<double> &fracMol0, Vcr<double> &fracMol1, Vcr<double> &fracMol0O, Vcr<double> &fracMol1O, Vcr<double> &fracMol0G, Vcr<double> &fracMol1G, double temperatureHigh, double temperatureLow, double dt, ProFlu *fluC, int ncomp) {

    PhaseFace left;
    PhaseFace right;
    PhaseSource source;
    PhaseCell cell;
    double flowArea = state.cells[i].duto.area;
    double dx = state.cells[i].dx;
    double temperatureLeft;
    if (state.cells[i].VTemper < 0.)
        temperatureLeft = state.cells[i].temp;
    else {
        if (i > 0)
            temperatureLeft = state.cells[i - 1].temp;
        else if (state.input.ConContEntrada == 1)
            temperatureLeft = state.inletTemperature;
        else
            temperatureLeft = state.cells[i].temp;
    }
    double temperatureRight = state.cells[i].temp;
    if (state.cells[i + 1].VTemper < 0.)
        temperatureRight = state.cells[i + 1].temp;

    source.tit = 0.;
    ProFlu fluF;

    source.oilVolumeFactor = 1.;
    source.waterVolumeFactor = 1.;
    source.waterCut = 1.;
    source.rhoO = 900.;
    source.rhoW = 1000.;

    if (i > 0 || state.input.ConContEntrada == 0) {
        if (i > 0 && state.cells[i].QG >= 0.)
            left.betI = state.cells[i - 1].betPigD;
        else
            left.betI = state.cells[i].betPigE;
    } else {
        if (state.cells[i].QG >= 0.)
            left.betI = state.inletCompletionFraction;
        else
            left.betI = state.cells[i].betPigE;
    }
    upwindLeftFacePhaseLiquidProperties(state, left, i, temperatureHigh, temperatureLow);
    if (left.oilVolumeFactor < 1e-15)
        left.oilVolumeFactor = 1e-15;

    right.betI = state.cells[i].betPigD;
    right.solutionGasRatio = state.cells[i].flui.RS(state.cells[i].pres, state.cells[i].temp);
    right.oilVolumeFactor = state.cells[i].flui.BOFunc(state.cells[i].pres, state.cells[i].temp, right.solutionGasRatio);
    right.waterVolumeFactor = state.cells[i].flui.BAFunc(state.cells[i].pres, state.cells[i].temp);
    right.bsw = state.cells[i].flui.BSW * right.waterVolumeFactor / (right.oilVolumeFactor + right.waterVolumeFactor * state.cells[i].flui.BSW - state.cells[i].flui.BSW * right.oilVolumeFactor);
    right.BSW = state.cells[i].flui.BSW;
    right.denag = state.cells[i].flui.Denag;
    right.viscL = 0 * 30 + 1 * state.cells[i].flui.VisOM(temperatureLow);
    right.viscH = 0 * 20 + 1 * state.cells[i].flui.VisOM(temperatureHigh);

    // betI1 = celula[i + 1].betPigE;    //duvidabeta
    if (state.cells[i + 1].QL < 0.) {
        right.betI = state.cells[i + 1].betPigE; // testeBet
        right.solutionGasRatio = state.cells[i + 1].flui.RS(state.cells[i + 1].pres, state.cells[i + 1].temp);
        right.oilVolumeFactor = state.cells[i + 1].flui.BOFunc(state.cells[i + 1].pres, state.cells[i + 1].temp, right.solutionGasRatio);
        right.waterVolumeFactor = state.cells[i + 1].flui.BAFunc(state.cells[i + 1].pres, state.cells[i + 1].temp);
        right.bsw = state.cells[i + 1].flui.BSW * right.waterVolumeFactor / (right.oilVolumeFactor + right.waterVolumeFactor * state.cells[i + 1].flui.BSW - state.cells[i + 1].flui.BSW * right.oilVolumeFactor);
        right.BSW = state.cells[i + 1].flui.BSW;
        right.denag = state.cells[i + 1].flui.Denag;
        right.viscL = 0 * 30 + 1 * state.cells[i + 1].flui.VisOM(temperatureLow);
        right.viscH = 0 * 20 + 1 * state.cells[i + 1].flui.VisOM(temperatureHigh);
    }
    if (right.oilVolumeFactor < 1e-15)
        right.oilVolumeFactor = 1e-15;
    left.vazMasLiq = state.cells[i].Mliqini;
    right.vazMasLiq = state.cells[i + 1].Mliqini;
    left.vazMasGas = state.cells[i].MC - state.cells[i].Mliqini;
    right.vazMasGas = state.cells[i + 1].MC - state.cells[i + 1].Mliqini;
    cell.pesoMol = 0.;
    for (int kfrac = 0; kfrac < ncomp; kfrac++) {
        cell.pesoMol += state.cells[i].flui.masMol[kfrac] * state.cells[i].flui.fracMol[kfrac];
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

    upwindLiquidFractions(state, left, right, i, fracMol0O, fracMol1O, ncomp);

    upwindVapourAndMixtureFractions(state, left, right, i, fracMol0, fracMol1, fracMol0G, fracMol1G, ncomp);

    source.dissolvedGas = state.cells[i].fontemassLR;
    source.freeGas = state.cells[i].fontemassGR;
    source.deadOil = state.cells[i].fontemassLR;
    source.water = state.cells[i].fontemassLR;

    readPhaseSourceFluid(state, source, fluF, i);

    source.pesoMol = 0;
    for (int j = 0; j < ncomp; j++) {
        source.pesoMol += fluF.masMol[j] * fluF.fracMol[j];
    }
    source.dissolvedGas *= source.tit;
    PhaseTransported transported;
    fluC[i] = state.cells[i].flui;
    int limCorrige = 1;
    if (state.compositionalRefreshCounter == state.input.miniTabAtraso  || state.input.miniTabAtraso == 0)
        limCorrige = 2;
    solveCellPhaseMolarFractions(state, left, right, cell, transported, source, limCorrige, fluF, i, fracMol0, fracMol1, dt, fluC, ncomp);
    transported.betIV = state.cells[i].bet;
    transported.solutionGasRatio = state.cells[i].flui.RS(state.cells[i].pres, state.cells[i].temp);
    transported.oilVolumeFactor = state.cells[i].flui.BOFunc(state.cells[i].pres, state.cells[i].temp, transported.solutionGasRatio);
    transported.waterVolumeFactor = state.cells[i].flui.BAFunc(state.cells[i].pres, state.cells[i].temp);
    transported.bswV = state.cells[i].flui.BSW * transported.waterVolumeFactor / (transported.oilVolumeFactor + transported.waterVolumeFactor * state.cells[i].flui.BSW - state.cells[i].flui.BSW * transported.oilVolumeFactor);
    transported.rhoOVol = state.cells[i].flui.MasEspoleo(state.cells[i].pres, state.cells[i].temp);
    transported.rhoWVol = state.cells[i].flui.MasEspAgua(state.cells[i].pres, state.cells[i].temp);
    transported.titVol = (1 - transported.bswV) * transported.rhoOVol / ((1 - transported.bswV) * transported.rhoOVol + transported.bswV * transported.rhoWVol);
    state.cells[i].nMol = (fluC[i].MasEspLiq(state.cells[i].pres, state.cells[i].temp) * (1. - state.cells[i].alf) * (1. - transported.betIV) * transported.titVol +
                      fluC[i].MasEspGas(state.cells[i].pres, state.cells[i].temp) *
                          state.cells[i].alf) *
                     state.cells[i].duto.area * state.cells[i].dx / cell.pesoMol;
    fluC[i].Pmol = cell.pesoMol;
    cell.liquidHoldup = 1. - state.cells[i].alf;
    cell.completionFraction = state.cells[i].bet;

    cell.solutionGasRatio = state.cells[i].flui.RS(state.cells[i].pres, state.cells[i].temp);
    cell.oilVolumeFactor = state.cells[i].flui.BOFunc(state.cells[i].pres, state.cells[i].temp, cell.solutionGasRatio);
    cell.waterVolumeFactor = state.cells[i].flui.BAFunc(state.cells[i].pres, state.cells[i].temp);
    cell.bsw = state.cells[i].flui.BSW * cell.waterVolumeFactor / (cell.oilVolumeFactor + cell.waterVolumeFactor * state.cells[i].flui.BSW - state.cells[i].flui.BSW * cell.oilVolumeFactor);

    cell.BSWini = state.cells[i].flui.BSW;
    source.BSW = state.cells[i].flui.BSW;
    cell.denagini = state.cells[i].flui.Denag;
    source.denag = state.cells[i].flui.Denag;
    cell.viscLini = 0 * 30 + 1 * state.cells[i].flui.VisOM(temperatureLow);
    cell.viscHini = 0 * 20 + 1 * state.cells[i].flui.VisOM(temperatureHigh);
    source.viscL = cell.viscLini;
    source.viscH = cell.viscHini;

    phaseSourceStandardRatesByAccessory(state, source, fluF, i, temperatureHigh, temperatureLow);
    if (source.tit > 1. - 1e-15) {
        source.dissolvedGas = 0.;
        source.deadOil = 0.;
        source.water = 0.;
    }

    PhaseBalance balance;
    balance.volpesFim = cell.liquidHoldup * (1 - cell.completionFraction) * (1 - cell.bsw) / cell.oilVolumeFactor;
    balance.volaguaFim = cell.liquidHoldup * (1 - cell.completionFraction) * cell.bsw / cell.oilVolumeFactor; // nao deveria ser dividido por Bo???????????
    balance.MultPd = 0.;
    balance.residuoP = 0.;
    balance.MultAe = 0.;
    balance.residuoA = 0.;
    transportPhaseVolumesAndViscosities(state, left, right, cell, source, balance, dx, flowArea, i, BSW, denag, VISCL, VISCH, dt);
    state.cells[i].VolPesaST = balance.volpesFim - 0. * balance.residuoP * dt / flowArea;
    state.cells[i].VolAguaST = balance.volaguaFim - 0. * balance.residuoA * dt / flowArea;
}

}  // namespace

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
    if (state.cells[0].acsr.tipo == kAccessoryRadialPorous && (state.cells[0].fontemassCR + state.cells[0].fontemassGR + state.cells[0].fontemassLR) > 0.)
        state.cells[0].flui.BSW = state.cells[0].acsr.radialPoro.BSW;
    else if (state.cells[0].acsr.tipo == kAccessoryPorous2D && (state.cells[0].fontemassCR + state.cells[0].fontemassGR + state.cells[0].fontemassLR) > 0.)
        state.cells[0].flui.BSW = state.cells[0].acsr.poroso2D.dados.transfer.BSW;
    double dt = state.cells[1].dt;
    double temperatureLow = 0.;
    double temperatureHigh = 70.;
    if (state.compositionalRefreshCounter == state.input.miniTabAtraso  && state.input.miniTabAtraso > 0)
        (*state.globals).modoTransiente = 0;
    int imin = 1;
    if ((*state.globals).chaverede != 0 && (state.cells[0].acsr.tipo == kAccessoryMultipleSource || state.input.ConContEntrada > 0))
        imin = 0;

    for (int i = imin; i < state.lastCell; i++) {
transportCellPhaseMolarFractions(state, i, BSW, denag, VISCL, VISCH, fracMol0, fracMol1, fracMol0O, fracMol1O, fracMol0G, fracMol1G, temperatureHigh, temperatureLow, dt, fluC, ncomp);
    }
    for (int i = imin; i <= state.lastCell - 1; i++) {
        if (state.input.nfluP > 0 || (*state.globals).chaverede != 0) {
            state.cells[i].flui = fluC[i];
            if (state.input.nfluP > 1 || (*state.globals).chaverede != 0) {
                state.cells[i].flui.BSW = BSW[i];
                state.cells[i].flui.Denag = denag[i];
                state.cells[i].flui.LVisL = VISCL[i];
                state.cells[i].flui.LVisH = VISCH[i];
                state.cells[i].flui.TempL = temperatureLow;
                state.cells[i].flui.TempH = temperatureHigh;
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
                state.cells[state.lastCell].flui.TempL = temperatureLow;
                state.cells[state.lastCell].flui.TempH = temperatureHigh;
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
