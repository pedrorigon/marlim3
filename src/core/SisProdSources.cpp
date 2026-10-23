#include "SisProdSources.h"
#include "SisProdSourcesInternal.h"

#include "Leitura.h"
#include "SisProdConstants.h"
#include "celula3.h"
#include "variaveisGlobais1D.h"
#include "UnitConversions.h"

#include <algorithm>
#include <math.h>

namespace sisprod::sources {

using enum sisprod::AccessoryKind;
using sisprod::kAirDensityAtStandardConditions;
using sisprod::kGravity;

namespace {

void consumeHydrateFormationMass(const SourceState &state, double &gas_consumido_Mg, double &agua_consumida_Mw, int cellIndex) {
    const bool model2 = state.input.calculaEnvelope == 1 && state.input.tipoHmodel == 2 && (*state.globals).lixo5 > 0.01;
    const bool model3 = state.input.calculaEnvelope == 1 && state.input.tipoHmodel == 3 && (*state.globals).lixo5 > 0.01;
    if (model2 || model3) {

        agua_consumida_Mw = state.cells[cellIndex].agua_consumida_massa_step;

        gas_consumido_Mg = state.cells[cellIndex].gas_consumido_massa_step;

        // Update the BSW
        double A_cross = state.cells[cellIndex].duto.area;
        double Lcel    = state.cells[cellIndex].dx;
        double Vlivre  = std::max(A_cross * Lcel - (model2 ? state.cells[cellIndex].V_h : state.cells[cellIndex].V_h_total), 1e-12);

        double frac_agua = std::max((1-state.cells[cellIndex].alfR)*(1-state.cells[cellIndex].betR)*state.cells[cellIndex].FW, 1e-12);
        double frac_oleo = std::max((1-state.cells[cellIndex].alfR)*(1-state.cells[cellIndex].betR)*(1-state.cells[cellIndex].FW), 1e-12);

        double Vagua = frac_agua * Vlivre;
        double Voil  = frac_oleo * Vlivre;

        double rho_w = std::max(state.cells[cellIndex].flui.MasEspAgua(state.cells[cellIndex].pres, state.cells[cellIndex].temp), 1e-12);
        double Vagua_new = Vagua - agua_consumida_Mw / rho_w;
        if (Vagua_new < 0.0) Vagua_new = 0.0;

        double den = Voil + Vagua_new;
        if (den > 1e-12) {
            state.cells[cellIndex].flui.BSW = Vagua_new / den;
        }

    }
}

double inSituWaterFraction(const Cel &cell, Cel &fluidCell) {
    double bo;
    double ba;
    if (cell.flui.RGO < 1e7) {
        double rs = fluidCell.flui.RS(cell.pres, cell.temp);
        bo = fluidCell.flui.BOFunc(cell.pres, cell.temp, rs);
        ba = fluidCell.flui.BAFunc(cell.pres, cell.temp);
    } else {
        bo = 1;
        ba = 0.;
    }
    return fluidCell.flui.BSW * ba / (bo + ba * fluidCell.flui.BSW - fluidCell.flui.BSW * bo);
}

void refreshChokeMultipleAndPorousSources(const SourceState &state, int cellIndex) {
    if (state.cells[cellIndex].acsr.tipo == kAccessoryLeak) {
        state.cells[cellIndex].acsr.fontechk.fluidoP = state.cells[cellIndex].flui;
        state.cells[cellIndex].acsr.fontechk.presT = state.cells[cellIndex].pres;
        state.cells[cellIndex].acsr.fontechk.tempT = state.cells[cellIndex].temp;
        double pres = state.cells[cellIndex].pres;
        double temp = state.cells[cellIndex].temp;
        double alf = state.cells[cellIndex].alf;
        double bet = state.cells[cellIndex].bet;
        double rhog = state.cells[cellIndex].flui.MasEspGas(pres, temp);
        double rhoP = state.cells[cellIndex].flui.MasEspLiq(pres, temp);
        double rhoC = state.cells[cellIndex].fluicol.MasEspFlu(pres, temp);
        state.cells[cellIndex].acsr.fontechk.titT = alf * rhog / (alf * rhog + (1 - alf) * (bet * rhoC + (1 - bet) * rhoP));
        state.cells[cellIndex].acsr.fontechk.betIST = bet;
        int iconex = -1;
        if (!((*state.globals).chaveRedeParalela == 0 || (*state.globals).iterRede > 0 || state.parallelSecondaryBranch == 1 || state.parallelSecondaryBoundaryCondition == 1)) {
            int nfonte = state.parallelSourceCells.size();
            for (int ifonte = 0; ifonte < nfonte; ifonte++) {
                if (state.parallelSourceCells[ifonte] == cellIndex) {
                    iconex = ifonte;
                    break;
                }
            }
        }
        if (iconex < 0) {
            state.cells[cellIndex].acsr.fontechk.VMas();
            state.cells[cellIndex].fontemassLR += state.cells[cellIndex].acsr.fontechk.masP;
            state.cells[cellIndex].fontemassCR += state.cells[cellIndex].acsr.fontechk.masC;
            state.cells[cellIndex].fontemassGR += state.cells[cellIndex].acsr.fontechk.masG;
        } else {
            state.cells[cellIndex].fontemassLR += state.parallelSourceProductionLiquid[iconex];
            state.cells[cellIndex].fontemassCR += state.parallelSourceComplementaryLiquid[iconex];
            state.cells[cellIndex].fontemassGR += state.parallelSourceGas[iconex];
        }
    }
    if (state.cells[cellIndex].acsr.tipo == kAccessoryMultipleSource) {
        if (state.cells[cellIndex].acsr.injm.condTermo == 1) {
            state.cells[cellIndex].fontemassLR = state.cells[cellIndex].acsr.injm.MassP;
            state.cells[cellIndex].fontemassCR = state.cells[cellIndex].acsr.injm.MassC;
            state.cells[cellIndex].fontemassGR = state.cells[cellIndex].acsr.injm.MassG;
        } else {
            state.cells[cellIndex].fontemassCR = state.cells[cellIndex].acsr.injm.MassC;
            double masT = state.cells[cellIndex].acsr.injm.MassP + state.cells[cellIndex].acsr.injm.MassG;
            double tit = state.cells[cellIndex].acsr.injm.FluidoPro.FracMassHidra(state.cells[cellIndex].pres, state.cells[cellIndex].temp);
            state.cells[cellIndex].acsr.injm.MassG = state.cells[cellIndex].fontemassGR = masT * tit;
            state.cells[cellIndex].acsr.injm.MassP = state.cells[cellIndex].fontemassLR = masT - state.cells[cellIndex].fontemassGR;
        }
    }
    if (state.cells[cellIndex].acsr.tipo == kAccessoryRadialPorous) {
        state.cells[cellIndex].acsr.radialPoro.pW.val[0] = state.cells[cellIndex].pres;
        // Use the next array cell's fluid at this cell's pressure and temperature.
        double vfw = inSituWaterFraction(state.cells[cellIndex], state.cells[cellIndex + 1]);
        state.cells[cellIndex].acsr.radialPoro.sWPoc = vfw * (1. - state.cells[cellIndex].acsr.radialPoro.satAconat) + state.cells[cellIndex].acsr.radialPoro.satAconat;
        state.cells[cellIndex].acsr.radialPoro.Pint = state.cells[cellIndex].pres;
        if (state.steadyMode == 0)
            state.cells[cellIndex].acsr.radialPoro.avancoPressao();
        else
            state.cells[cellIndex].acsr.radialPoro.pseudoTrans();
        state.cells[cellIndex].fontemassLR = state.cells[cellIndex].acsr.radialPoro.fluxIni + state.cells[cellIndex].acsr.radialPoro.fluxIniA;
        state.cells[cellIndex].fontemassCR = 0.;
        state.cells[cellIndex].fontemassGR = state.cells[cellIndex].acsr.radialPoro.fluxIniG;
    }
    if (state.cells[cellIndex].acsr.tipo == kAccessoryPorous2D) {
        state.cells[cellIndex].acsr.poroso2D.dados.pW.val[0] = state.cells[cellIndex].pres;
        // Use the next array cell's fluid at this cell's pressure and temperature.
        double vfw = inSituWaterFraction(state.cells[cellIndex], state.cells[cellIndex + 1]);
        state.cells[cellIndex].acsr.poroso2D.sWPoc = vfw * (1. - state.cells[cellIndex].acsr.poroso2D.dados.satAconat) + state.cells[cellIndex].acsr.poroso2D.dados.satAconat;
        state.cells[cellIndex].acsr.poroso2D.dados.transfer.sWPoc = state.cells[cellIndex].acsr.poroso2D.sWPoc;
        state.cells[cellIndex].acsr.poroso2D.dados.pInt = state.cells[cellIndex].pres;
        state.cells[cellIndex].acsr.poroso2D.dados.transfer.Pint = state.cells[cellIndex].pres;
        if (state.steadyMode == 0)
            state.cells[cellIndex].acsr.poroso2D.avancoPressao();
        else
            state.cells[cellIndex].acsr.poroso2D.pseudoTransientePoroso();
        state.cells[cellIndex].fontemassLR = state.cells[cellIndex].acsr.poroso2D.dados.transfer.fluxIni + state.cells[cellIndex].acsr.poroso2D.dados.transfer.fluxIniA;
        state.cells[cellIndex].fontemassCR = 0.;
        state.cells[cellIndex].fontemassGR = state.cells[cellIndex].acsr.poroso2D.dados.transfer.fluxIniG;
    }
}

}  // namespace

void addMasterValveFlow(const SourceState &state, int cellIndex) {
    Cel &cell = state.cells[cellIndex];
    if (cell.acsr.chk.AreaGarg < (1e-3 + state.input.master1.razareaativ) * cell.duto.area && cell.acsr.chk.AreaGarg > 1e-3 * cell.duto.area) {
        Cel &nextCell = state.cells[cellIndex + 1];
        const double tE = cell.temp;
        const double alfE = cell.alf;
        const double betE = cell.bet;
        const double sense = 1.;

        double rholp = cell.flui.MasEspLiq(cell.pres, cell.temp);
        double rholc = cell.fluicol.MasEspFlu(cell.pres, cell.temp);
        double rholmix = (1 - betE) * rholp + betE * rholc;

        double alfJ = nextCell.alf;
        double betJ = nextCell.bet;
        double rholpJ = nextCell.flui.MasEspLiq(nextCell.pres, nextCell.temp);
        double rholcJ = nextCell.fluicol.MasEspFlu(nextCell.pres, nextCell.temp);
        double rholmixJ = (1 - betJ) * rholpJ + betJ * rholcJ;

        double hidroM = sin(cell.duto.teta) * (0.5 * cell.dx) * (rholmix * (1 - alfE) + alfE * cell.flui.MasEspGas(cell.pres, cell.temp)) * kGravity / units::kPascalPerKgfPerCm2;
        double hidroJ = sin(nextCell.duto.teta) * (0.5 * nextCell.dx) * (rholmixJ * (1 - alfJ) + alfJ * nextCell.flui.MasEspGas(nextCell.pres, nextCell.temp)) * kGravity / units::kPascalPerKgfPerCm2;

        double masentrada = cell.MC;
        double massgas = cell.MC - cell.Mliqini;
        double tit;
        if (fabs(masentrada) < 1e-9 && fabs(massgas) < 1e-9) {
            tit = alfE * cell.flui.MasEspGas(cell.pres, cell.temp) / (cell.flui.MasEspGas(cell.pres, cell.temp) * alfE + rholmix * (1. - alfE));
        } else {
            if ((massgas >= 0 && cell.Mliqini <= 0) || (massgas < 0 && cell.Mliqini == 0))
                tit = 1.;
            else if (massgas <= 0 && cell.Mliqini > 0)
                tit = 0.;
            else if (masentrada < 0)
                tit = 1.;
            else
                tit = fabs(massgas / masentrada);
            if (tit > 1)
                tit = 1;
        }

        const double ypres = (nextCell.pres + hidroJ) / (cell.pres - hidroM);
        const double masChk = cell.acsr.chk.vazmassSachd(ypres, cell.pres - hidroM, tE, alfE,
                                                      betE, tit, cell.flui, cell.fluicol);
        double maxSup = cell.acsr.chk.vazmaxSachd(cell.pres - hidroM, tE, alfE,
                                               betE, tit, cell.flui, cell.fluicol);

        if (fabs(ypres) > fabs(cell.acsr.chk.razpres))
            maxSup = masChk;

        if (cell.acsr.chk.AreaGarg < (1e-3) * cell.duto.area || (cell.pres < nextCell.pres))
            maxSup = 0.;

        const double masliq = sense * maxSup * (1. - tit);
        const double masgas = sense * maxSup * tit;

        cell.fontemassLR += (-1 * masliq * (1 - betE) * rholp / rholmix);
        cell.fontemassCR += (-1 * masliq * betE * rholc / rholmix);
        cell.fontemassGR += (-1 * masgas);

        nextCell.fontemassLR += (masliq * (1 - betE) * rholp / rholmix);
        nextCell.fontemassCR += (masliq * betE * rholc / rholmix);
        nextCell.fontemassGR += (masgas);
    }
}

void renewSourceTerms(const SourceState &state, int cellIndex) {

    double pr = state.cells[cellIndex].pres;
    double tr = state.cells[cellIndex].temp;
    if (cellIndex == 0) {
        state.cells[cellIndex].fontemassLR = 0.;
        state.cells[cellIndex].fontemassCR = 0.;
        state.cells[cellIndex].fontemassGR = 0.;
    } else {
        if (state.cells[cellIndex - 1].acsr.tipo != kAccessoryChoke && state.cells[cellIndex - 1].acsr.tipo != kAccessoryVolumetricPump && (cellIndex < state.lastCell || state.cells[state.lastCell].acsr.tipo == kAccessoryInflowPerformance)) {
            state.cells[cellIndex].fontemassLR = 0.;
            state.cells[cellIndex].fontemassCR = 0.;
            state.cells[cellIndex].fontemassGR = 0.;
        }
    }

    double agua_consumida_Mw = 0.;
    double gas_consumido_Mg = 0.;

    consumeHydrateFormationMass(state, gas_consumido_Mg, agua_consumida_Mw, cellIndex);
    if (state.cells[cellIndex].acsr.tipo == kAccessoryGasInjection) {
        if (state.cells[cellIndex].acsr.injg.tipoflu == 0) {
            double masgas = state.cells[cellIndex].acsr.injg.VMas(pr, tr);
            if (fabs(masgas) < (*state.globals).localtiny)
                masgas = 0.;
            if (state.cells[cellIndex].acsr.injg.seco == 1) {
                state.cells[cellIndex].fontemassGR += masgas;
                state.cells[cellIndex].fontemassLR = 0.;
                state.cells[cellIndex].fontemassCR = 0.;
            } else {
                double tit;
                if (state.input.flashCompleto != 2)
                    tit = state.cells[cellIndex].acsr.injg.FluidoPro.FracMassHidra(1., 20.);
                else
                    tit = state.cells[cellIndex].acsr.injg.FluidoPro.dStockTankVaporMassFraction;
                double masT = masgas / tit;
                tit = state.cells[cellIndex].acsr.injg.FluidoPro.FracMassHidra(pr, tr);
                state.cells[cellIndex].fontemassGR += masT * tit;
                state.cells[cellIndex].fontemassLR += masT * (1. - tit);
                double rcomp = state.cells[cellIndex].acsr.injg.fluidocol.standardDensity();
                state.cells[cellIndex].fontemassCR += rcomp * state.cells[cellIndex].acsr.injg.razCompGas * state.cells[cellIndex].acsr.injg.QGas / units::kSecondsPerDay;
            }
        } else {
            state.cells[cellIndex].fontemassGR += 0.;
            state.cells[cellIndex].fontemassLR = 0.;
            state.cells[cellIndex].fontemassCR = state.cells[cellIndex].acsr.injg.QGas;
        }
    }
    if (state.cells[cellIndex].acsr.tipo == kAccessoryLiquidInjection) {
        double rlcA = state.cells[cellIndex].acsr.injl.fluidocol.standardDensity();
        state.cells[cellIndex].fontemassCR += rlcA * state.cells[cellIndex].acsr.injl.QLiq * state.cells[cellIndex].acsr.injl.bet / units::kSecondsPerDay;
        double massic = state.cells[cellIndex].acsr.injl.QLiq * (1. - state.cells[cellIndex].acsr.injl.bet) / units::kSecondsPerDay;
        double Rhogs = state.cells[cellIndex].acsr.injl.FluidoPro.Deng * kAirDensityAtStandardConditions;
        double Rhols = (1000 * 141.5 / (131.5 + state.cells[cellIndex].acsr.injl.FluidoPro.API)) * (1 - state.cells[cellIndex].acsr.injl.FluidoPro.BSW) + 1000. * state.cells[cellIndex].acsr.injl.FluidoPro.Denag * state.cells[cellIndex].acsr.injl.FluidoPro.BSW;
        double multiplicador = (Rhols + state.cells[cellIndex].acsr.injl.FluidoPro.RGO * Rhogs * (1 - state.cells[cellIndex].acsr.injl.FluidoPro.BSW));
        massic *= multiplicador;
        double fracmasshidra = state.cells[cellIndex].acsr.injl.FluidoPro.FracMassHidra(pr, tr);
        state.cells[cellIndex].fontemassLR += (1. - fracmasshidra) * massic;
        state.cells[cellIndex].fontemassGR += fracmasshidra * massic;
    }
    if (state.cells[cellIndex].acsr.tipo == kAccessoryInflowPerformance) {
        if (pr < state.cells[cellIndex].acsr.ipr.Pres) {
            state.cells[cellIndex].fontemassLR += state.cells[cellIndex].acsr.ipr.MasL(pr, tr);
            state.cells[cellIndex].fontemassCR = 0.;
            state.cells[cellIndex].fontemassGR += state.cells[cellIndex].acsr.ipr.MasG(pr, tr);
        } else {
            double tit;
            tit = state.cells[cellIndex].alf * state.cells[cellIndex].flui.MasEspGas(pr, tr) /
                  (state.cells[cellIndex].alf * state.cells[cellIndex].flui.MasEspGas(pr, tr) +
                   (1. - state.cells[cellIndex].alf) * (1. - state.cells[cellIndex].bet) * state.cells[cellIndex].flui.MasEspLiq(pr, tr) +
                   (1. - state.cells[cellIndex].alf) * state.cells[cellIndex].bet * state.cells[cellIndex].fluicol.MasEspFlu(pr, tr));
            state.cells[cellIndex].fontemassLR += (1. - state.cells[cellIndex].alf) * (1. - state.cells[cellIndex].bet) * state.cells[cellIndex].acsr.ipr.VMas(pr, tr) * state.cells[cellIndex].flui.MasEspLiq(pr, tr);
            state.cells[cellIndex].fontemassCR += (1. - state.cells[cellIndex].alf) * state.cells[cellIndex].bet * state.cells[cellIndex].acsr.ipr.VMas(pr, tr) * state.cells[cellIndex].fluicol.MasEspFlu(pr, tr);
            state.cells[cellIndex].fontemassGR += 1 * state.cells[cellIndex].alf * state.cells[cellIndex].acsr.ipr.VMas(pr, tr) * state.cells[cellIndex].flui.MasEspGas(pr, tr);
            state.cells[cellIndex].acsr.ipr.deriP *= (1. - state.cells[cellIndex].alf) * (1. - state.cells[cellIndex].bet) * state.cells[cellIndex].flui.MasEspLiq(pr, tr);
            state.cells[cellIndex].acsr.ipr.deriC *= (1. - state.cells[cellIndex].alf) * state.cells[cellIndex].bet * state.cells[cellIndex].fluicol.MasEspFlu(pr, tr);
            state.cells[cellIndex].acsr.ipr.deriG *= 1 * state.cells[cellIndex].alf * state.cells[cellIndex].flui.MasEspGas(pr, tr);
        }
    }
    if (state.cells[cellIndex].acsr.tipo == kAccessoryChoke) {
        state.cells[cellIndex + 1].fontemassLR = 0.;
        state.cells[cellIndex + 1].fontemassCR = 0.;
        state.cells[cellIndex + 1].fontemassGR = 0.;
        if (state.steadyMode == 0)
            addMasterValveFlow(state, cellIndex);
    }
    if (state.cells[cellIndex].acsr.tipo == kAccessoryVolumetricPump) {
        state.cells[cellIndex + 1].fontemassLR = 0.;
        state.cells[cellIndex + 1].fontemassCR = 0.;
        state.cells[cellIndex + 1].fontemassGR = 0.;
        state.cells[cellIndex].acsr.bvol.fluido = state.cells[cellIndex].flui;
        state.cells[cellIndex].acsr.bvol.fluicol = state.cells[cellIndex].fluicol;
        if (fabs(state.cells[cellIndex].acsr.bvol.freq) > 1) {
            double alfM = state.cells[cellIndex].alf;
            double betM = state.cells[cellIndex].bet;
            double presM = state.cells[cellIndex].pres;
            double tempM = state.cells[cellIndex].temp;
            double presM1 = state.cells[cellIndex].presR;
            double tempM1 = state.cells[cellIndex].tempR;
            state.cells[cellIndex].acsr.bvol.vazmass(presM, tempM, presM1, tempM1, betM, alfM);
            state.cells[cellIndex].fontemassLR -= state.cells[cellIndex].acsr.bvol.MLiqP;
            state.cells[cellIndex].fontemassCR -= state.cells[cellIndex].acsr.bvol.MLiqC;
            state.cells[cellIndex].fontemassGR -= state.cells[cellIndex].acsr.bvol.MGas;
            state.cells[cellIndex + 1].fontemassLR += state.cells[cellIndex].acsr.bvol.MLiqP;
            state.cells[cellIndex + 1].fontemassCR += state.cells[cellIndex].acsr.bvol.MLiqC;
            state.cells[cellIndex + 1].fontemassGR += state.cells[cellIndex].acsr.bvol.MGas;
        }
    }
    refreshChokeMultipleAndPorousSources(state, cellIndex);

    const bool hydrateModel = state.input.calculaEnvelope == 1 && (state.input.tipoHmodel == 2 || state.input.tipoHmodel == 3);
    if (hydrateModel && state.cells[cellIndex].flui.BSW > 1e-12 && (*state.globals).lixo5 > 0.01) {
        state.cells[cellIndex].fontemassLR -= (agua_consumida_Mw / (*state.globals).lixo5);
        state.cells[cellIndex].fontemassGR -= (gas_consumido_Mg / (*state.globals).lixo5);
    }
}

}  // namespace sisprod::sources
