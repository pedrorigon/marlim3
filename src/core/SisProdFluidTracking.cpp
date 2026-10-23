#include "SisProdFluidTracking.h"
#include "SisProdFluidTrackingInternal.h"

#include "Leitura.h"
#include "SisProdConstants.h"
#include "celula3.h"
#include "variaveisGlobais1D.h"
#include "UnitConversions.h"

#include <math.h>

namespace sisprod::fluidtracking {

namespace {

void transportBlackOilBalances(const FluidTrackingState &state, const BlackOilFace &left, const BlackOilFace &right, const BlackOilCell &cell, const BlackOilSource &source, BlackOilBalance &balance, double dx, double flowArea, int i, Vcr<double> &rgo, Vcr<double> &dg, Vcr<double> &yco2, Vcr<double> &API, Vcr<double> &BSW, Vcr<double> &denag, Vcr<double> &VISCL, Vcr<double> &VISCH, double dt) {
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
    if (cell.liquidHoldup > (*state.globals).localtiny && cell.complementaryFraction < (1. - (*state.globals).localtiny) && cell.bsw < (1. - (*state.globals).localtiny)) {
        rgo[i] = (balance.volleveFim - balance.residuo * dt / flowArea) * cell.oilVolumeFactor / (cell.liquidHoldup * (1 - cell.complementaryFraction) * (1 - cell.bsw));
        if (rgo[i] > (*state.globals).RGOMax)
            rgo[i] = (*state.globals).RGOMax;
    } else if (cell.complementaryFraction >= (1. - (*state.globals).localtiny) || cell.bsw >= (1. - (*state.globals).localtiny))
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

void blackOilSourceCellFluid(const FluidTrackingState &state, BlackOilSource &source, int i, double temperatureHigh, double temperatureLow) {
    if (state.cells[i].acsr.tipo == kAccessoryChoke || state.cells[i].acsr.tipo == kAccessoryVolumetricPump) {
        source.dgO = state.cells[i].flui.Deng;
        source.yco2O = state.cells[i].flui.yco2;
        source.rgo = state.cells[i].flui.RGO;
        source.stockTankQuality = state.cells[i].flui.dStockTankVaporMassFraction;
        double rholiq = (1 - state.cells[i].flui.BSW) * (1000 * 141.5 / (131.5 + state.cells[i].flui.API)) + state.cells[i].flui.BSW * 1000 * state.cells[i].flui.Denag;
        double rhogas = state.cells[i].flui.Deng * kAirDensityAtStandardConditions;
        source.rholST = rholiq + rhogas * source.rgo * (1. - state.cells[i].flui.BSW);
        source.rhogST = rhogas;
        double solutionGasRatioSource = state.cells[i].flui.RS(state.cells[i].pres, state.cells[i].temp) * (units::kBarrelPerCubicMetre / units::kCubicFootPerCubicMetre);
        source.razdgd = 1 / state.cells[i].flui.rDgD;
        source.razdgl = 1 / state.cells[i].flui.rDgL;
        source.dissolvedGas = (source.dissolvedGas + source.freeGas) * (source.razdgd * solutionGasRatioSource * (1. - state.cells[i].flui.BSW) / source.rholST);
        if (state.input.nfluP > 1 || (*state.globals).chaverede != 0) {
            double rhoPSTF = (1 - state.cells[i].flui.BSW) * (1000 * 141.5 / (131.5 + state.cells[i].flui.API)) + state.cells[i].flui.BSW * 1000 * state.cells[i].flui.Denag;
            source.deadOil *= ((1 - state.cells[i].flui.BSW) / rhoPSTF);
            source.water *= (state.cells[i].flui.BSW / rhoPSTF);
            source.API = state.cells[i].flui.API;
            source.BSW = state.cells[i].flui.BSW;
            source.viscL = state.cells[i].flui.VisOM(temperatureLow);
            source.viscH = state.cells[i].flui.VisOM(temperatureHigh);
        }
    } else if (state.cells[i].acsrL != 0 && ((*state.cells[i].acsrL).tipo == kAccessoryChoke || (*state.cells[i].acsrL).tipo == kAccessoryVolumetricPump)) {
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
            double solutionGasRatioSource = state.cells[i - 1].flui.RS(state.cells[i].pres, state.cells[i].temp) * (units::kBarrelPerCubicMetre / units::kCubicFootPerCubicMetre);
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
                source.viscL = state.cells[i - 1].flui.VisOM(temperatureLow);
                source.viscH = state.cells[i - 1].flui.VisOM(temperatureHigh);
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
            (void)(state.cells[i].flui.RS(state.cells[i].pres, state.cells[i].temp) * (units::kBarrelPerCubicMetre / units::kCubicFootPerCubicMetre));
            source.razdgd = 1 / state.cells[i].flui.rDgD;
            source.razdgl = 1 / state.cells[i].flui.rDgL;
            source.dissolvedGas = (source.dissolvedGas + source.freeGas) * (source.razdgd * state.cells[i].flui.RS(state.cells[i].pres, state.cells[i].temp) * (units::kBarrelPerCubicMetre / units::kCubicFootPerCubicMetre) *
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
                source.viscL = state.cells[i].flui.VisOM(temperatureLow);
                source.viscH = state.cells[i].flui.VisOM(temperatureHigh);
            }
        }
    }
}

void blackOilSourcePorous2D(const FluidTrackingState &state, BlackOilSource &source, ProFlu &fluF, int i, double temperatureHigh, double temperatureLow) {
    if ((state.cells[i].fontemassLR + state.cells[i].fontemassGR) > 1e-15)
        fluF = state.cells[i].acsr.radialPoro.flup;
    else
        fluF = state.cells[i].flui;

    source.stockTankQuality = fluF.dStockTankVaporMassFraction;
    source.dgO = fluF.Deng;
    source.yco2O = fluF.yco2;
    source.rgo = fluF.RGO;

    source.rholST = (1 - state.cells[i].acsr.poroso2D.dados.transfer.BSW) * (1000 * 141.5 / (131.5 + fluF.API)) + state.cells[i].acsr.poroso2D.dados.transfer.BSW * 1000 * fluF.Denag + fluF.Deng * kAirDensityAtStandardConditions * source.rgo * (1. - state.cells[i].acsr.poroso2D.dados.transfer.BSW);
    double solutionGasRatioSource = fluF.RS(state.cells[i].pres, state.cells[i].temp) * (units::kBarrelPerCubicMetre / units::kCubicFootPerCubicMetre);
    source.razdgd = 1 / fluF.rDgD;
    source.razdgl = 1 / fluF.rDgL;
    source.dissolvedGas = (source.dissolvedGas + source.freeGas) * source.razdgd * (fluF.RS(state.cells[i].pres, state.cells[i].temp) * (units::kBarrelPerCubicMetre / units::kCubicFootPerCubicMetre) * (1. - state.cells[i].acsr.poroso2D.dados.transfer.BSW) / source.rholST);
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
        source.viscL = fluF.VisOM(temperatureLow);
        source.viscH = fluF.VisOM(temperatureHigh);
    }
}

void blackOilSourceRadialPorous(const FluidTrackingState &state, BlackOilSource &source, ProFlu &fluF, int i, double temperatureHigh, double temperatureLow) {
    if ((state.cells[i].fontemassLR + state.cells[i].fontemassGR) > 1e-15)
        fluF = state.cells[i].acsr.radialPoro.flup;
    else
        fluF = state.cells[i].flui;

    source.stockTankQuality = fluF.dStockTankVaporMassFraction;
    source.dgO = fluF.Deng;
    source.yco2O = fluF.yco2;
    source.rgo = fluF.RGO;

    source.rholST = (1 - state.cells[i].acsr.radialPoro.BSW) * (1000 * 141.5 / (131.5 + fluF.API)) + state.cells[i].acsr.radialPoro.BSW * 1000 * fluF.Denag + fluF.Deng * kAirDensityAtStandardConditions * source.rgo * (1. - state.cells[i].acsr.radialPoro.BSW);
    double solutionGasRatioSource = fluF.RS(state.cells[i].pres, state.cells[i].temp) * (units::kBarrelPerCubicMetre / units::kCubicFootPerCubicMetre);
    source.razdgd = 1 / fluF.rDgD;
    source.razdgl = 1 / fluF.rDgL;
    source.dissolvedGas = (source.dissolvedGas + source.freeGas) * source.razdgd * (fluF.RS(state.cells[i].pres, state.cells[i].temp) * (units::kBarrelPerCubicMetre / units::kCubicFootPerCubicMetre) * (1. - state.cells[i].acsr.radialPoro.BSW) / source.rholST);
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
        source.viscL = fluF.VisOM(temperatureLow);
        source.viscH = fluF.VisOM(temperatureHigh);
    }
}

void blackOilSourceLeak(const FluidTrackingState &state, BlackOilSource &source, int i, double temperatureHigh, double temperatureLow) {
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

    double solutionGasRatioSource = fluF.RS(state.cells[i].pres, state.cells[i].temp) * (units::kBarrelPerCubicMetre / units::kCubicFootPerCubicMetre);
    source.razdgd = 1 / fluF.rDgD;
    source.razdgl = 1 / fluF.rDgL;
    if (state.cells[i].acsr.fontechk.ambGas != 1 || (source.dissolvedGas + source.freeGas) < 0.)
        source.dissolvedGas = (source.dissolvedGas + source.freeGas) * source.razdgd * (fluF.RS(state.cells[i].pres, state.cells[i].temp) * (units::kBarrelPerCubicMetre / units::kCubicFootPerCubicMetre) * (1. - fluF.BSW) / source.rholST);
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
        source.viscL = fluF.VisOM(temperatureLow);
        source.viscH = fluF.VisOM(temperatureHigh);
    }
}

void blackOilSourceInflowPerformance(const FluidTrackingState &state, BlackOilSource &source, ProFlu &fluF, int i, double temperatureHigh, double temperatureLow) {
    if ((state.cells[i].acsr.ipr.Pres) > state.cells[i].pres)
        fluF = state.cells[i].acsr.ipr.FluidoPro;
    else
        fluF = state.cells[i].flui;

    source.stockTankQuality = fluF.dStockTankVaporMassFraction;
    source.dgO = fluF.Deng;
    source.yco2O = fluF.yco2;
    source.rgo = fluF.RGO;

    source.rholST = (1 - fluF.BSW) * (1000 * 141.5 / (131.5 + fluF.API)) + fluF.BSW * 1000 * fluF.Denag + fluF.Deng * kAirDensityAtStandardConditions * source.rgo * (1. - fluF.BSW);
    double solutionGasRatioSource = fluF.RS(state.cells[i].pres, state.cells[i].temp) * (units::kBarrelPerCubicMetre / units::kCubicFootPerCubicMetre);
    source.razdgd = 1 / fluF.rDgD;
    source.razdgl = 1 / fluF.rDgL;
    source.dissolvedGas = (source.dissolvedGas + source.freeGas) * source.razdgd * (fluF.RS(state.cells[i].pres, state.cells[i].temp) * (units::kBarrelPerCubicMetre / units::kCubicFootPerCubicMetre) * (1. - fluF.BSW) / source.rholST);
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
        source.viscL = fluF.VisOM(temperatureLow);
        source.viscH = fluF.VisOM(temperatureHigh);
    }
}

void blackOilSourceMultipleSource(const FluidTrackingState &state, BlackOilSource &source, ProFlu &fluF, int i, double temperatureHigh, double temperatureLow) {
    if ((state.cells[i].acsr.injm.MassC + state.cells[i].acsr.injm.MassG + state.cells[i].acsr.injm.MassP) > 0.)
        fluF = state.cells[i].acsr.injm.FluidoPro;
    else
        fluF = state.cells[i].flui;

    source.stockTankQuality = fluF.dStockTankVaporMassFraction;
    source.dgO = fluF.Deng;
    source.yco2O = fluF.yco2;
    source.rgo = fluF.RGO;

    double solutionGasRatioSource = fluF.RS(state.cells[i].pres, state.cells[i].temp) * (units::kBarrelPerCubicMetre / units::kCubicFootPerCubicMetre);
    source.razdgd = 1 / fluF.rDgD;
    source.razdgl = 1 / fluF.rDgL;
    source.rholST = (1 - fluF.BSW) * (1000 * 141.5 / (131.5 + fluF.API)) +
              fluF.BSW * 1000 *
                  fluF.Denag +
              fluF.Deng * kAirDensityAtStandardConditions * source.rgo * (1. - fluF.BSW);
    source.dissolvedGas = (source.dissolvedGas + source.freeGas) * source.razdgd * (fluF.RS(state.cells[i].pres, state.cells[i].temp) * (units::kBarrelPerCubicMetre / units::kCubicFootPerCubicMetre) * (1. - fluF.BSW) / source.rholST);

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
        source.viscL = fluF.VisOM(temperatureLow);
        source.viscH = fluF.VisOM(temperatureHigh);
    }
}

void blackOilSourceLiquidInjection(const FluidTrackingState &state, BlackOilSource &source, ProFlu &fluF, int i, double temperatureHigh, double temperatureLow) {
    if (state.cells[i].acsr.injl.QLiq > 0.)
        fluF = state.cells[i].acsr.injl.FluidoPro;
    else
        fluF = state.cells[i].flui;

    source.stockTankQuality = fluF.dStockTankVaporMassFraction;
    source.dgO = fluF.Deng;
    source.yco2O = fluF.yco2;
    source.rgo = fluF.RGO;

    double solutionGasRatioSource = fluF.RS(state.cells[i].pres, state.cells[i].temp) * (units::kBarrelPerCubicMetre / units::kCubicFootPerCubicMetre);
    source.razdgd = 1 / fluF.rDgD;
    source.razdgl = 1 / fluF.rDgL;
    if (fluF.BSW < 1 - (*state.globals).localtiny)
        source.rholST = (1 - fluF.BSW) * (1000 * 141.5 / (131.5 + fluF.API)) + fluF.BSW * 1000 * fluF.Denag + fluF.Deng * kAirDensityAtStandardConditions * source.rgo * (1. - fluF.BSW);
    else
        source.rholST = fluF.BSW * 1000 * fluF.Denag;
    source.dissolvedGas = (source.dissolvedGas + source.freeGas) * source.razdgd * (fluF.RS(state.cells[i].pres, state.cells[i].temp) * (units::kBarrelPerCubicMetre / units::kCubicFootPerCubicMetre) * (1. - fluF.BSW) / source.rholST);

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
        source.viscL = fluF.VisOM(temperatureLow);
        source.viscH = fluF.VisOM(temperatureHigh);
    }
}

void blackOilSourceWetGasInjection(const FluidTrackingState &state, BlackOilSource &source, ProFlu &fluF, int i, double temperatureHigh, double temperatureLow) {
    if (state.cells[i].acsr.injg.QGas > 0.)
        fluF = state.cells[i].acsr.injg.FluidoPro;
    else
        fluF = state.cells[i].flui;

    source.stockTankQuality = fluF.dStockTankVaporMassFraction;
    source.dgO = fluF.Deng;
    source.yco2O = fluF.yco2;
    source.rgo = fluF.RGO;

    double solutionGasRatioSource = fluF.RS(state.cells[i].pres, state.cells[i].temp) * (units::kBarrelPerCubicMetre / units::kCubicFootPerCubicMetre);
    source.razdgd = 1 / fluF.rDgD;
    source.razdgl = 1 / fluF.rDgL;
    if (state.cells[i].acsr.injg.FluidoPro.BSW < 1 - (*state.globals).localtiny)
        source.rholST = (1 - fluF.BSW) * (1000 * 141.5 / (131.5 + fluF.API)) + fluF.BSW * 1000 * fluF.Denag + fluF.Deng * kAirDensityAtStandardConditions * source.rgo * (1. - fluF.BSW);
    else
        source.rholST = fluF.BSW * 1000 * fluF.Denag;
    source.dissolvedGas = (source.dissolvedGas + source.freeGas) * source.razdgd * (fluF.RS(state.cells[i].pres, state.cells[i].temp) * (units::kBarrelPerCubicMetre / units::kCubicFootPerCubicMetre) * (1. - fluF.BSW) / source.rholST);

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
        source.viscL = fluF.VisOM(temperatureLow);
        source.viscH = fluF.VisOM(temperatureHigh);
    }
}

void upwindLeftFaceGasProperties(const FluidTrackingState &state, BlackOilFace &left, int i) {
    if ((i > 0 || state.input.ConContEntrada == 1) && state.cells[i].QG > 0) {
        left.rhogST = (*state.cells[i].fluiL).Deng * kAirDensityAtStandardConditions;
        left.dgG = (*state.cells[i].fluiL).Deng;
        left.yco2G = (*state.cells[i].fluiL).yco2;
    } else {
        left.rhogST = state.cells[i].flui.Deng * kAirDensityAtStandardConditions;
        left.dgG = state.cells[i].flui.Deng;
        left.yco2G = state.cells[i].flui.yco2;
    }
}

void upwindRightFaceLiquidProperties(const FluidTrackingState &state, BlackOilFace &right, int i, double temperatureHigh, double temperatureLow) {
    if (state.cells[i + 1].QL < 0.) {
        right.betI = state.cells[i + 1].betPigE;
        right.solutionGasRatio = state.cells[i + 1].flui.RS(state.cells[i + 1].pres, state.cells[i + 1].temp);
        right.oilVolumeFactor = state.cells[i + 1].flui.BOFunc(state.cells[i + 1].pres, state.cells[i + 1].temp, right.solutionGasRatio);
        right.waterVolumeFactor = state.cells[i + 1].flui.BAFunc(state.cells[i + 1].pres, state.cells[i + 1].temp);
        right.bsw = state.cells[i + 1].flui.BSW * right.waterVolumeFactor / (right.oilVolumeFactor + right.waterVolumeFactor * state.cells[i + 1].flui.BSW - state.cells[i + 1].flui.BSW * right.oilVolumeFactor);
        right.solutionGasRatio = right.solutionGasRatio * units::kBarrelPerCubicMetre / units::kCubicFootPerCubicMetre;
        right.dgO = state.cells[i + 1].flui.Deng;
        right.razdgd = 1 / state.cells[i + 1].flui.rDgD;
        right.razdgl = 1 / state.cells[i + 1].flui.rDgL;
        right.yco2O = state.cells[i + 1].flui.yco2;
        right.API = state.cells[i + 1].flui.API;
        right.BSW = state.cells[i + 1].flui.BSW;
        right.denag = state.cells[i + 1].flui.Denag;
        right.viscL = state.cells[i + 1].flui.VisOM(temperatureLow);
        right.viscH = state.cells[i + 1].flui.VisOM(temperatureHigh);
    }
}

void upwindLeftFaceBlackOilLiquid(const FluidTrackingState &state, BlackOilFace &left, int i, double temperatureHigh, double temperatureLow) {
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
            left.betI = state.cells[i - 1].betPigD;
        else
            left.betI = state.inletComplementaryFraction;
        left.solutionGasRatio = (*state.cells[i].fluiL).RS(upstreamPressure, upstreamTemperature);
        left.oilVolumeFactor = (*state.cells[i].fluiL).BOFunc(upstreamPressure, upstreamTemperature, left.solutionGasRatio);
        left.waterVolumeFactor = (*state.cells[i].fluiL).BAFunc(upstreamPressure, upstreamTemperature);
        left.bsw = (*state.cells[i].fluiL).BSW * left.waterVolumeFactor / (left.oilVolumeFactor + left.waterVolumeFactor * (*state.cells[i].fluiL).BSW - (*state.cells[i].fluiL).BSW * left.oilVolumeFactor);
        left.solutionGasRatio = left.solutionGasRatio * units::kBarrelPerCubicMetre / units::kCubicFootPerCubicMetre;
        left.dgO = (*state.cells[i].fluiL).Deng;
        left.razdgd = 1 / (*state.cells[i].fluiL).rDgD;
        left.razdgl = 1 / (*state.cells[i].fluiL).rDgL;
        left.yco2O = (*state.cells[i].fluiL).yco2;
        left.API = (*state.cells[i].fluiL).API;
        left.BSW = (*state.cells[i].fluiL).BSW;
        left.denag = (*state.cells[i].fluiL).Denag;
        left.viscL = (*state.cells[i].fluiL).VisOM(temperatureLow);
        left.viscH = (*state.cells[i].fluiL).VisOM(temperatureHigh);
    } else {
        left.betI = state.cells[i].betPigE;
        left.solutionGasRatio = state.cells[i].flui.RS(state.cells[i].pres, state.cells[i].temp);
        left.oilVolumeFactor = state.cells[i].flui.BOFunc(state.cells[i].pres, state.cells[i].temp, left.solutionGasRatio);
        left.waterVolumeFactor = state.cells[i].flui.BAFunc(state.cells[i].pres, state.cells[i].temp);
        left.bsw = state.cells[i].flui.BSW * left.waterVolumeFactor / (left.oilVolumeFactor + left.waterVolumeFactor * state.cells[i].flui.BSW - state.cells[i].flui.BSW * left.oilVolumeFactor);
        left.solutionGasRatio = left.solutionGasRatio * units::kBarrelPerCubicMetre / units::kCubicFootPerCubicMetre;
        left.dgO = state.cells[i].flui.Deng;
        left.razdgd = 1 / state.cells[i].flui.rDgD;
        left.razdgl = 1 / state.cells[i].flui.rDgL;
        left.yco2O = state.cells[i].flui.yco2;
        left.API = state.cells[i].flui.API;
        left.BSW = state.cells[i].flui.BSW;
        left.denag = state.cells[i].flui.Denag;
        left.viscL = state.cells[i].flui.VisOM(temperatureLow);
        left.viscH = state.cells[i].flui.VisOM(temperatureHigh);
    }
}

void transportCellBlackOilProperties(const FluidTrackingState &state, int i, Vcr<double> &rgo, Vcr<double> &dg, Vcr<double> &yco2, Vcr<double> &API, Vcr<double> &BSW, Vcr<double> &denag, Vcr<double> &VISCL, Vcr<double> &VISCH, double temperatureHigh, double temperatureLow, double dt) {
    BlackOilBalance balance;
    double flowArea = state.cells[i].duto.area;
    double dx = state.cells[i].dx;
    BlackOilFace left;
    upwindLeftFaceBlackOilLiquid(state, left, i, temperatureHigh, temperatureLow);
    if (left.oilVolumeFactor < 1e-15)
        left.oilVolumeFactor = 1e-15;

    BlackOilFace right;
    right.betI = state.cells[i].betPigD;
    right.solutionGasRatio = state.cells[i].flui.RS(state.cells[i].pres, state.cells[i].temp);
    right.oilVolumeFactor = state.cells[i].flui.BOFunc(state.cells[i].pres, state.cells[i].temp, right.solutionGasRatio);
    right.waterVolumeFactor = state.cells[i].flui.BAFunc(state.cells[i].pres, state.cells[i].temp);
    right.bsw = state.cells[i].flui.BSW * right.waterVolumeFactor / (right.oilVolumeFactor + right.waterVolumeFactor * state.cells[i].flui.BSW - state.cells[i].flui.BSW * right.oilVolumeFactor);
    right.solutionGasRatio = right.solutionGasRatio * units::kBarrelPerCubicMetre / units::kCubicFootPerCubicMetre;
    right.dgO = state.cells[i].flui.Deng;
    right.yco2O = state.cells[i].flui.yco2;
    right.API = state.cells[i].flui.API;
    right.BSW = state.cells[i].flui.BSW;
    right.denag = state.cells[i].flui.Denag;
    right.viscL = state.cells[i].flui.VisOM(temperatureLow);
    right.viscH = state.cells[i].flui.VisOM(temperatureHigh);
    right.razdgd = 1 / state.cells[i].flui.rDgD;
    right.razdgl = 1 / state.cells[i].flui.rDgL;

    upwindRightFaceLiquidProperties(state, right, i, temperatureHigh, temperatureLow);
    if (right.oilVolumeFactor < 1e-15)
        right.oilVolumeFactor = 1e-15;

    upwindLeftFaceGasProperties(state, left, i);

    right.rhogST = state.cells[i].flui.Deng * kAirDensityAtStandardConditions;
    right.dgG = state.cells[i].flui.Deng;
    right.yco2G = state.cells[i].flui.yco2;
    if (state.cells[i + 1].QG <= 0.) {
        right.rhogST = state.cells[i + 1].flui.Deng * kAirDensityAtStandardConditions;
        right.dgG = state.cells[i + 1].flui.Deng;
        right.yco2G = state.cells[i + 1].flui.yco2;
    }

    BlackOilCell cell;
    cell.liquidHoldup = 1. - state.cells[i].alf;
    cell.complementaryFraction = state.cells[i].bet;
    cell.rholST = (1 - state.cells[i].flui.BSW) * (1000 * 141.5 / (131.5 + state.cells[i].flui.API)) + state.cells[i].flui.BSW * 1000 * state.cells[i].flui.Denag;
    cell.rhog = state.cells[i].rgC;
    cell.rhogST = state.cells[i].flui.Deng * kAirDensityAtStandardConditions;
    cell.solutionGasRatio = state.cells[i].flui.RS(state.cells[i].pres, state.cells[i].temp);
    cell.razdgd = 1 / state.cells[i].flui.rDgD;
    cell.razdgl = 1 / state.cells[i].flui.rDgL;
    cell.oilVolumeFactor = state.cells[i].flui.BOFunc(state.cells[i].pres, state.cells[i].temp, cell.solutionGasRatio);
    cell.waterVolumeFactor = state.cells[i].flui.BAFunc(state.cells[i].pres, state.cells[i].temp);
    cell.bsw = state.cells[i].flui.BSW * cell.waterVolumeFactor / (cell.oilVolumeFactor + cell.waterVolumeFactor * state.cells[i].flui.BSW - state.cells[i].flui.BSW * cell.oilVolumeFactor);
    cell.solutionGasRatio = cell.solutionGasRatio * units::kBarrelPerCubicMetre / units::kCubicFootPerCubicMetre;

    cell.dgini = state.cells[i].flui.Deng;
    cell.yco2ini = state.cells[i].flui.yco2;
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
    cell.viscLini = state.cells[i].flui.VisOM(temperatureLow);
    cell.viscHini = state.cells[i].flui.VisOM(temperatureHigh);
    source.viscL = cell.viscLini;
    source.viscH = cell.viscHini;
    source.rholST = cell.rholST;
    source.rhogST = cell.rhogST;
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
    balance.volleveFim = (((1 - cell.liquidHoldup) * cell.rhog * cell.razdgl / (cell.rhogST)) + cell.liquidHoldup * (1 - cell.complementaryFraction) * (1. - cell.bsw) * cell.solutionGasRatio * cell.razdgd / (cell.oilVolumeFactor));
    if (balance.volleveFim < 1e-15)
        balance.volleveFim = 0.;
    balance.residuo = (balance.volleveFim - state.cells[i].VolLeveST) * flowArea / dt + (balance.MultOd - balance.MultOe) / dx + (balance.MultGd - balance.MultGe) / dx - (source.dissolvedGas / dx + source.freeGas / dx);
    balance.volpesFim = cell.liquidHoldup * (1 - cell.complementaryFraction) * (1 - cell.bsw) / cell.oilVolumeFactor;
    balance.volaguaFim = cell.liquidHoldup * (1 - cell.complementaryFraction) * cell.bsw;
    transportBlackOilBalances(state, left, right, cell, source, balance, dx, flowArea, i, rgo, dg, yco2, API, BSW, denag, VISCL, VISCH, dt);
    state.cells[i].VolLeveST = balance.volleveFim;
    state.cells[i].VolPesaST = balance.volpesFim;
    state.cells[i].VolAguaST = balance.volaguaFim;
}

}  // namespace

void transportBlackOilProperties(const FluidTrackingState &state, ProFlu fluiRev) {
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

void normalizePhaseFractions(std::vector<double>& fracMolFase, double* fracMolOriginal, int npseudo) {
    for (int kfrac = 0; kfrac < npseudo; kfrac++) {
    	fracMolFase[kfrac] = fracMolOriginal[kfrac];
    }
	// Find the smallest value
    double menorFracFase = fracMolFase[0];
    for (int kfrac = 1; kfrac < npseudo; kfrac++) {
        if (menorFracFase > fracMolFase[kfrac]) {
            menorFracFase = fracMolFase[kfrac];
        }
    }
    menorFracFase=0.;

    // Subtract the smallest value
    for (int kfrac = 0; kfrac < npseudo; kfrac++) {
    	fracMolFase[kfrac] -= menorFracFase;
    }

    // Compute the total sum
    double fracTotFase = 0.;
    for (int kfrac = 0; kfrac < npseudo; kfrac++) {
        fracTotFase += fracMolFase[kfrac];
    }

    // Normalise
    if (fracTotFase > 0) {
    	for (int kfrac = 0; kfrac < npseudo; kfrac++) {
    		fracMolFase[kfrac] /= fracTotFase;
        }
    }
}

void transportPhaseVolumesAndViscosities(const FluidTrackingState &state, const PhaseFace &left, const PhaseFace &right, const PhaseCell &cell, const PhaseSource &source, PhaseBalance &balance, double dx, double flowArea, int i, Vcr<double> &BSW, Vcr<double> &denag, Vcr<double> &VISCL, Vcr<double> &VISCH, double dt) {
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

void phaseSourceStandardRatesByAccessory(const FluidTrackingState &state, PhaseSource &source, ProFlu &fluF, int i, double temperatureHigh, double temperatureLow) {
    if (state.cells[i].acsr.tipo == kAccessoryLiquidInjection) {
        double solutionGasRatioSource = state.cells[i].acsr.injl.FluidoPro.RS(state.cells[i].pres, state.cells[i].temp) * (units::kBarrelPerCubicMetre / units::kCubicFootPerCubicMetre);
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
            source.viscL = state.cells[i].acsr.injl.FluidoPro.VisOM(temperatureLow);
            source.viscH = state.cells[i].acsr.injl.FluidoPro.VisOM(temperatureHigh);
        }
    } else if (state.cells[i].acsr.tipo == kAccessoryMultipleSource) {
        double solutionGasRatioSource = state.cells[i].acsr.injm.FluidoPro.RS(state.cells[i].pres, state.cells[i].temp) * (units::kBarrelPerCubicMetre / units::kCubicFootPerCubicMetre);
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
            source.denag = state.cells[i].acsr.injm.FluidoPro.Denag;
            source.viscL = state.cells[i].acsr.injm.FluidoPro.VisOM(temperatureLow);
            source.viscH = state.cells[i].acsr.injm.FluidoPro.VisOM(temperatureHigh);
        }
    } else if (state.cells[i].acsr.tipo == kAccessoryInflowPerformance) {
        double solutionGasRatioSource = state.cells[i].acsr.ipr.FluidoPro.RS(state.cells[i].pres, state.cells[i].temp) * (units::kBarrelPerCubicMetre / units::kCubicFootPerCubicMetre);
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
            source.viscL = state.cells[i].acsr.ipr.FluidoPro.VisOM(temperatureLow);
            source.viscH = state.cells[i].acsr.ipr.FluidoPro.VisOM(temperatureHigh);
        }
    } else if (state.cells[i].acsr.tipo == kAccessoryLeak) {

        double solutionGasRatioSource = fluF.RS(state.cells[i].pres, state.cells[i].temp) * (units::kBarrelPerCubicMetre / units::kCubicFootPerCubicMetre);
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
            source.viscL = fluF.VisOM(temperatureLow);
            source.viscH = fluF.VisOM(temperatureHigh);
        }
    } else if (state.cells[i].acsr.tipo == kAccessoryRadialPorous) {
        double solutionGasRatioSource = state.cells[i].acsr.radialPoro.flup.RS(state.cells[i].pres, state.cells[i].temp) * (units::kBarrelPerCubicMetre / units::kCubicFootPerCubicMetre);
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
            source.viscL = state.cells[i].acsr.radialPoro.flup.VisOM(temperatureLow);
            source.viscH = state.cells[i].acsr.radialPoro.flup.VisOM(temperatureHigh);
        }
    } else if (state.cells[i].acsr.tipo == kAccessoryPorous2D) {
        double solutionGasRatioSource = state.cells[i].acsr.poroso2D.dados.flup.RS(state.cells[i].pres, state.cells[i].temp) * (units::kBarrelPerCubicMetre / units::kCubicFootPerCubicMetre);
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
            source.viscL = state.cells[i].acsr.poroso2D.dados.flup.VisOM(temperatureLow);
            source.viscH = state.cells[i].acsr.poroso2D.dados.flup.VisOM(temperatureHigh);
        }
    } else if ((fabs(source.dissolvedGas) > (*state.globals).localtiny && state.cells[i].acsr.tipo != kAccessoryLiquidInjection && state.cells[i].acsr.tipo != kAccessoryInflowPerformance &&
                state.cells[i].acsr.tipo != kAccessoryLeak && state.cells[i].acsr.tipo != kAccessoryRadialPorous && state.cells[i].acsr.tipo != kAccessoryPorous2D) ||
               (fabs(source.freeGas) > (*state.globals).localtiny && state.cells[i].acsr.tipo != kAccessoryGasInjection && state.cells[i].acsr.tipo != kAccessoryLiquidInjection && state.cells[i].acsr.tipo != kAccessoryInflowPerformance && state.cells[i].acsr.tipo != kAccessoryLeak && state.cells[i].acsr.tipo != kAccessoryRadialPorous && state.cells[i].acsr.tipo != kAccessoryPorous2D)) {
        if (state.cells[i].acsr.tipo == kAccessoryChoke || state.cells[i].acsr.tipo == kAccessoryVolumetricPump) {
            (void)(state.cells[i].flui.RS(state.cells[i].pres, state.cells[i].temp) * (units::kBarrelPerCubicMetre / units::kCubicFootPerCubicMetre));
            if (state.input.nfluP > 0 || (*state.globals).chaverede != 0) {
                double rhoPSTF = (1 - state.cells[i].flui.BSW) * (1000 * 141.5 / (131.5 + state.cells[i].flui.API)) + state.cells[i].flui.BSW * 1000 * state.cells[i].flui.Denag;
                source.deadOil *= ((1 - state.cells[i].flui.BSW) / rhoPSTF);
                source.water *= (state.cells[i].flui.BSW / rhoPSTF);
                source.BSW = state.cells[i].flui.BSW;
                source.denag = fluF.Denag;
                source.viscL = state.cells[i].flui.VisOM(temperatureLow);
                source.viscH = state.cells[i].flui.VisOM(temperatureHigh);
            }
        } else if (state.cells[i].acsrL != 0 && ((*state.cells[i].acsrL).tipo == kAccessoryChoke || (*state.cells[i].acsrL).tipo == kAccessoryVolumetricPump)) {
            if (i > 0) {
                (void)(state.cells[i - 1].flui.RS(state.cells[i].pres, state.cells[i].temp) * (units::kBarrelPerCubicMetre / units::kCubicFootPerCubicMetre));
                if (state.input.nfluP > 0 || (*state.globals).chaverede != 0) {
                    double rhoPSTF = (1 - state.cells[i - 1].flui.BSW) * (1000 * 141.5 / (131.5 + state.cells[i - 1].flui.API)) +
                                     state.cells[i - 1].flui.BSW * 1000 *
                                         state.cells[i - 1].flui.Denag;
                    source.deadOil *= ((1 - state.cells[i - 1].flui.BSW) / rhoPSTF);
                    source.water *= (state.cells[i - 1].flui.BSW / rhoPSTF);
                    source.BSW = state.cells[i - 1].flui.BSW;
                    source.denag = state.cells[i - 1].flui.Denag;
                    source.viscL = state.cells[i - 1].flui.VisOM(temperatureLow);
                    source.viscH = state.cells[i - 1].flui.VisOM(temperatureHigh);
                }
            } else {
                (void)(state.cells[i].flui.RS(state.cells[i].pres, state.cells[i].temp) * (units::kBarrelPerCubicMetre / units::kCubicFootPerCubicMetre));
                if (state.input.nfluP > 0 || (*state.globals).chaverede != 0) {
                    double rhoPSTF = (1 - state.cells[i].flui.BSW) * (1000 * 141.5 / (131.5 + state.cells[i].flui.API)) +
                                     state.cells[i].flui.BSW * 1000 *
                                         state.cells[i].flui.Denag;
                    source.deadOil *= ((1 - state.cells[i].flui.BSW) / rhoPSTF);
                    source.water *= (state.cells[i].flui.BSW / rhoPSTF);
                    source.BSW = state.cells[i].flui.BSW;
                    source.denag = state.cells[i].flui.Denag;
                    source.viscL = state.cells[i].flui.VisOM(temperatureLow);
                    source.viscH = state.cells[i].flui.VisOM(temperatureHigh);
                }
            }
        }
    }
}

void solveCellPhaseMolarFractions(const FluidTrackingState &state, const PhaseFace &left, const PhaseFace &right, PhaseCell &cell, PhaseTransported &transported, const PhaseSource &source, int limCorrige, ProFlu &fluF, int i, Vcr<double> &fracMol0, Vcr<double> &fracMol1, double dt, ProFlu *fluC, int ncomp) {
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
        }
    }
}

void readPhaseSourceFluid(const FluidTrackingState &state, PhaseSource &source, ProFlu &fluF, int i) {
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

void upwindVapourAndMixtureFractions(const FluidTrackingState &state, PhaseFace &left, PhaseFace &right, int i, Vcr<double> &fracMol0, Vcr<double> &fracMol1, Vcr<double> &fracMol0G, Vcr<double> &fracMol1G, int ncomp) {
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

void upwindLiquidFractions(const FluidTrackingState &state, PhaseFace &left, PhaseFace &right, int i, Vcr<double> &fracMol0O, Vcr<double> &fracMol1O, int ncomp) {
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

void upwindLeftFacePhaseLiquidProperties(const FluidTrackingState &state, PhaseFace &left, int i, double temperatureHigh, double temperatureLow) {
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
            left.betI = state.cells[i - 1].betPigD;
        else
            left.betI = state.inletComplementaryFraction;
        double solutionGasRatioLeft = (*state.cells[i].fluiL).RS(upstreamPressure, upstreamTemperature);
        left.oilVolumeFactor = (*state.cells[i].fluiL).BOFunc(upstreamPressure, upstreamTemperature, solutionGasRatioLeft);
        left.waterVolumeFactor = (*state.cells[i].fluiL).BAFunc(upstreamPressure, upstreamTemperature);
        left.bsw = (*state.cells[i].fluiL).BSW * left.waterVolumeFactor / (left.oilVolumeFactor + left.waterVolumeFactor * (*state.cells[i].fluiL).BSW - (*state.cells[i].fluiL).BSW * left.oilVolumeFactor);
        left.BSW = (*state.cells[i].fluiL).BSW;
        left.denag = (*state.cells[i].fluiL).Denag;
        left.viscL = (*state.cells[i].fluiL).VisOM(temperatureLow);
        left.viscH = (*state.cells[i].fluiL).VisOM(temperatureHigh);
    } else {
        left.betI = state.cells[i].betPigE;
        double solutionGasRatioLeft = state.cells[i].flui.RS(state.cells[i].pres, state.cells[i].temp);
        left.oilVolumeFactor = state.cells[i].flui.BOFunc(state.cells[i].pres, state.cells[i].temp, solutionGasRatioLeft);
        left.waterVolumeFactor = state.cells[i].flui.BAFunc(state.cells[i].pres, state.cells[i].temp);
        left.bsw = state.cells[i].flui.BSW * left.waterVolumeFactor / (left.oilVolumeFactor + left.waterVolumeFactor * state.cells[i].flui.BSW - state.cells[i].flui.BSW * left.oilVolumeFactor);
        left.BSW = state.cells[i].flui.BSW;
        left.denag = state.cells[i].flui.Denag;
        left.viscL = state.cells[i].flui.VisOM(temperatureLow);
        left.viscH = state.cells[i].flui.VisOM(temperatureHigh);
    }
}

void transportCellPhaseMolarFractions(const FluidTrackingState &state, int i, Vcr<double> &BSW, Vcr<double> &denag, Vcr<double> &VISCL, Vcr<double> &VISCH, Vcr<double> &fracMol0, Vcr<double> &fracMol1, Vcr<double> &fracMol0O, Vcr<double> &fracMol1O, Vcr<double> &fracMol0G, Vcr<double> &fracMol1G, double temperatureHigh, double temperatureLow, double dt, ProFlu *fluC, int ncomp) {

    PhaseFace left;
    PhaseFace right;
    PhaseSource source;
    PhaseCell cell;
    double flowArea = state.cells[i].duto.area;
    double dx = state.cells[i].dx;
    ProFlu fluF;

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
    right.viscL = state.cells[i].flui.VisOM(temperatureLow);
    right.viscH = state.cells[i].flui.VisOM(temperatureHigh);

    if (state.cells[i + 1].QL < 0.) {
        right.betI = state.cells[i + 1].betPigE;
        right.solutionGasRatio = state.cells[i + 1].flui.RS(state.cells[i + 1].pres, state.cells[i + 1].temp);
        right.oilVolumeFactor = state.cells[i + 1].flui.BOFunc(state.cells[i + 1].pres, state.cells[i + 1].temp, right.solutionGasRatio);
        right.waterVolumeFactor = state.cells[i + 1].flui.BAFunc(state.cells[i + 1].pres, state.cells[i + 1].temp);
        right.bsw = state.cells[i + 1].flui.BSW * right.waterVolumeFactor / (right.oilVolumeFactor + right.waterVolumeFactor * state.cells[i + 1].flui.BSW - state.cells[i + 1].flui.BSW * right.oilVolumeFactor);
        right.BSW = state.cells[i + 1].flui.BSW;
        right.denag = state.cells[i + 1].flui.Denag;
        right.viscL = state.cells[i + 1].flui.VisOM(temperatureLow);
        right.viscH = state.cells[i + 1].flui.VisOM(temperatureHigh);
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
    cell.complementaryFraction = state.cells[i].bet;

    cell.solutionGasRatio = state.cells[i].flui.RS(state.cells[i].pres, state.cells[i].temp);
    cell.oilVolumeFactor = state.cells[i].flui.BOFunc(state.cells[i].pres, state.cells[i].temp, cell.solutionGasRatio);
    cell.waterVolumeFactor = state.cells[i].flui.BAFunc(state.cells[i].pres, state.cells[i].temp);
    cell.bsw = state.cells[i].flui.BSW * cell.waterVolumeFactor / (cell.oilVolumeFactor + cell.waterVolumeFactor * state.cells[i].flui.BSW - state.cells[i].flui.BSW * cell.oilVolumeFactor);

    cell.BSWini = state.cells[i].flui.BSW;
    source.BSW = state.cells[i].flui.BSW;
    cell.denagini = state.cells[i].flui.Denag;
    source.denag = state.cells[i].flui.Denag;
    cell.viscLini = state.cells[i].flui.VisOM(temperatureLow);
    cell.viscHini = state.cells[i].flui.VisOM(temperatureHigh);
    source.viscL = cell.viscLini;
    source.viscH = cell.viscHini;

    phaseSourceStandardRatesByAccessory(state, source, fluF, i, temperatureHigh, temperatureLow);
    if (source.tit > 1. - 1e-15) {
        source.dissolvedGas = 0.;
        source.deadOil = 0.;
        source.water = 0.;
    }

    PhaseBalance balance;
    balance.volpesFim = cell.liquidHoldup * (1 - cell.complementaryFraction) * (1 - cell.bsw) / cell.oilVolumeFactor;
    balance.volaguaFim = cell.liquidHoldup * (1 - cell.complementaryFraction) * cell.bsw / cell.oilVolumeFactor;
    balance.residuoP = 0.;
    balance.residuoA = 0.;
    transportPhaseVolumesAndViscosities(state, left, right, cell, source, balance, dx, flowArea, i, BSW, denag, VISCL, VISCH, dt);
    state.cells[i].VolPesaST = balance.volpesFim - 0. * balance.residuoP * dt / flowArea;
    state.cells[i].VolAguaST = balance.volaguaFim - 0. * balance.residuoA * dt / flowArea;
}

}  // namespace

void transportPhaseMolarFractions(const FluidTrackingState &state, ProFlu fluiRev) {
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
    }
    delete[] fluC;
    if(state.input.miniTabAtraso > 0)(*state.globals).modoTransiente = 1;
}

}  // namespace sisprod::fluidtracking
