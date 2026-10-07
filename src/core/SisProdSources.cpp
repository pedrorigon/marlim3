#include "SisProdSources.h"

#include "Leitura.h"
#include "SisProdConstants.h"
#include "celula3.h"
#include "variaveisGlobais1D.h"

#include <algorithm>
#include <math.h>

namespace sisprod::sources {

using enum sisprod::AccessoryKind;
using sisprod::kAirDensityAtStandardConditions;
using sisprod::kBarrelPerCubicMetre;
using sisprod::kCubicFootPerCubicMetre;
using sisprod::kGravity;
using sisprod::kPascalPerKgfPerCm2Coarse;
using sisprod::kSecondsPerDay;

namespace {

/// With the hydrate envelope on (models 2 and 3) and past the first 0.01 s, takes
/// the water and gas that hydrate formation consumed in the cell during the step,
/// hands them back through the two out-parameters, and lowers the cell's BSW for
/// the free water that is gone.
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

        double BSW_old = state.cells[cellIndex].flui.BSW;
        double den = Voil + Vagua_new;
        if (den > 1e-12) {
            state.cells[cellIndex].flui.BSW = Vagua_new / den;
        } else {
            state.cells[cellIndex].flui.BSW = BSW_old;
        }

    } // hydrate change
}

/// In-situ water fraction of the fluid of fluidCell, at the pressure and
/// temperature of cell.
double inSituWaterFraction(const Cel &cell, Cel &fluidCell) {
    double rs;
    double bo;
    double ba;
    if (cell.flui.RGO < 1e7) {
        rs = fluidCell.flui.RS(cell.pres, cell.temp);
        bo = fluidCell.flui.BOFunc(cell.pres, cell.temp, rs);
        ba = fluidCell.flui.BAFunc(cell.pres, cell.temp);
        rs = rs * kBarrelPerCubicMetre / kCubicFootPerCubicMetre;
    } else {
        bo = 1;
        rs = 0;
        ba = 0.;
    }
    return fluidCell.flui.BSW * ba / (bo + ba * fluidCell.flui.BSW - fluidCell.flui.BSW * bo);
}

/// Adds to the cell the mass its source delivers this step when the source is a
/// choke source (type 9 -- on the first iteration of a parallel network, on its
/// primary side, the flow recorded for that connection instead), a multiple
/// source (10) or a radial or 2D porous medium (15, 16).
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
        // in-situ BSW of the previous cell (in the march, cell i)
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
        // in-situ BSW of the previous cell (in the march, cell i)
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
    if (state.cells[cellIndex].acsr.chk.AreaGarg < (1e-3 + state.input.master1.razareaativ) * state.cells[cellIndex].duto.area && state.cells[cellIndex].acsr.chk.AreaGarg > 1e-3 * state.cells[cellIndex].duto.area) {
        double tE = state.cells[cellIndex].temp;
        double alfE = state.cells[cellIndex].alf;
        double betE = state.cells[cellIndex].bet;
        double sense = 1.;

        double maxSup = 0.;

        double rholp = state.cells[cellIndex].flui.MasEspLiq(state.cells[cellIndex].pres, state.cells[cellIndex].temp);
        double rholc = state.cells[cellIndex].fluicol.MasEspFlu(state.cells[cellIndex].pres, state.cells[cellIndex].temp);
        double rholmix = (1 - betE) * rholp + betE * rholc;

        double alfJ = state.cells[cellIndex + 1].alf;
        double betJ = state.cells[cellIndex + 1].bet;
        double rholpJ = state.cells[cellIndex + 1].flui.MasEspLiq(state.cells[cellIndex + 1].pres, state.cells[cellIndex + 1].temp);
        double rholcJ = state.cells[cellIndex + 1].fluicol.MasEspFlu(state.cells[cellIndex + 1].pres, state.cells[cellIndex + 1].temp);
        double rholmixJ = (1 - betJ) * rholpJ + betJ * rholcJ;

        double hidroM = sin(state.cells[cellIndex].duto.teta) * (0.5 * state.cells[cellIndex].dx) * (rholmix * (1 - alfE) + alfE * state.cells[cellIndex].flui.MasEspGas(state.cells[cellIndex].pres, state.cells[cellIndex].temp)) * kGravity / kPascalPerKgfPerCm2Coarse;
        double hidroJ = sin(state.cells[cellIndex + 1].duto.teta) * (0.5 * state.cells[cellIndex + 1].dx) * (rholmixJ * (1 - alfJ) + alfJ * state.cells[cellIndex + 1].flui.MasEspGas(state.cells[cellIndex + 1].pres, state.cells[cellIndex + 1].temp)) * kGravity / kPascalPerKgfPerCm2Coarse;

        double masentrada = state.cells[cellIndex].MC;
        double massgas = state.cells[cellIndex].MC - state.cells[cellIndex].Mliqini;
        double tit;
        if (fabs(masentrada) < 1e-9 && fabs(massgas) < 1e-9) {
            tit = alfE * state.cells[cellIndex].flui.MasEspGas(state.cells[cellIndex].pres, state.cells[cellIndex].temp) / (state.cells[cellIndex].flui.MasEspGas(state.cells[cellIndex].pres, state.cells[cellIndex].temp) * alfE + rholmix * (1. - alfE));
        } else {
            if ((massgas >= 0 && state.cells[cellIndex].Mliqini <= 0) || (massgas < 0 && state.cells[cellIndex].Mliqini == 0))
                tit = 1.;
            else if (massgas <= 0 && state.cells[cellIndex].Mliqini > 0)
                tit = 0.;
            else if (masentrada < 0)
                tit = 1.;
            else
                tit = fabs(massgas / masentrada);
            if (tit > 1)
                tit = 1;
        }

        double masChk;

        double ypres = (state.cells[cellIndex + 1].pres + hidroJ) / (state.cells[cellIndex].pres - hidroM);
        int check = 1;
        if (ypres < 1. || check == 1) {
            masChk = state.cells[cellIndex].acsr.chk.vazmassSachd(ypres, state.cells[cellIndex].pres - hidroM, tE, alfE,
                                                       betE, tit, state.cells[cellIndex].flui, state.cells[cellIndex].fluicol);
            maxSup = state.cells[cellIndex].acsr.chk.vazmaxSachd(state.cells[cellIndex].pres - hidroM, tE, alfE,
                                                      betE, tit, state.cells[cellIndex].flui, state.cells[cellIndex].fluicol);
        } else {
            tE = state.cells[cellIndex + 1].temp;
            alfE = state.cells[cellIndex + 1].alf;
            betE = state.cells[cellIndex + 1].bet;
            ypres = 1. / ypres;
            sense = -1;
            masChk = state.cells[cellIndex].acsr.chk.vazmassSachd(ypres, state.cells[cellIndex + 1].pres + hidroJ, tE, alfE, betE, tit, state.cells[cellIndex + 1].flui,
                                                       state.cells[cellIndex + 1].fluicol);
            maxSup = state.cells[cellIndex].acsr.chk.vazmaxSachd(state.cells[cellIndex + 1].pres + hidroJ, tE, alfE, betE, tit, state.cells[cellIndex + 1].flui,
                                                      state.cells[cellIndex + 1].fluicol);
        }

        if (fabs(ypres) > fabs(state.cells[cellIndex].acsr.chk.razpres))
            maxSup = masChk;

        if (state.cells[cellIndex].acsr.chk.AreaGarg < (1e-3) * state.cells[cellIndex].duto.area || (state.cells[cellIndex].pres < state.cells[cellIndex + 1].pres && check == 1))
            maxSup =
                0.;

        double masliq;
        double masgas;
        masliq = sense * maxSup * (1. - tit);
        masgas = sense * maxSup * tit;

        state.cells[cellIndex].fontemassLR += (-1 * masliq * (1 - betE) * rholp / rholmix);
        state.cells[cellIndex].fontemassCR += (-1 * masliq * betE * rholc / rholmix);
        state.cells[cellIndex].fontemassGR += (-1 * masgas);

        state.cells[cellIndex + 1].fontemassLR += (masliq * (1 - betE) * rholp / rholmix);
        state.cells[cellIndex + 1].fontemassCR += (masliq * betE * rholc / rholmix);
        state.cells[cellIndex + 1].fontemassGR += (masgas);
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
                double rcomp = state.cells[cellIndex].acsr.injg.fluidocol.MasEspFlu(1., 20.);
                state.cells[cellIndex].fontemassCR += rcomp * state.cells[cellIndex].acsr.injg.razCompGas * state.cells[cellIndex].acsr.injg.QGas / kSecondsPerDay;
            }
        } else {
            state.cells[cellIndex].fontemassGR += 0.;
            state.cells[cellIndex].fontemassLR = 0.;
            state.cells[cellIndex].fontemassCR = state.cells[cellIndex].acsr.injg.QGas;
        }
    }
    if (state.cells[cellIndex].acsr.tipo == kAccessoryLiquidInjection) {
        double rlcA = state.cells[cellIndex].acsr.injl.fluidocol.MasEspFlu(1.001, 15.);
        state.cells[cellIndex].fontemassCR += rlcA * state.cells[cellIndex].acsr.injl.QLiq * state.cells[cellIndex].acsr.injl.bet / kSecondsPerDay;
        double massic = state.cells[cellIndex].acsr.injl.QLiq * (1. - state.cells[cellIndex].acsr.injl.bet) / kSecondsPerDay;
        double Rhogs = state.cells[cellIndex].acsr.injl.FluidoPro.Deng * kAirDensityAtStandardConditions; // cel[cellIndex].acsr.injl.FluidoPro.MasEspGas(1, 15);
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

    const bool hydrateModel = state.input.tipoHmodel == 2 || state.input.tipoHmodel == 3;
    if (state.input.calculaEnvelope == 1 && hydrateModel && state.cells[cellIndex].flui.BSW > 1e-12 && (*state.globals).lixo5 > 0.01) {
        state.cells[cellIndex].fontemassLR -= (agua_consumida_Mw / (*state.globals).lixo5);
        state.cells[cellIndex].fontemassGR -= (gas_consumido_Mg / (*state.globals).lixo5);
    }
}

}  // namespace sisprod::sources
