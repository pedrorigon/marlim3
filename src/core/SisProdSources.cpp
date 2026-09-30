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
/// the water and gas that hydrate formation consumed in cell ind during the step,
/// hands them back through the two out-parameters, and lowers the cell's BSW for
/// the free water that is gone.
void consumeHydrateFormationMass(const SourceState &state, double &gas_consumido_Mg, double &agua_consumida_Mw, int ind) {
    if (state.input.calculaEnvelope == 1 && state.input.tipoHmodel == 2 && (*state.globals).lixo5 > 0.01) {

        agua_consumida_Mw = state.cells[ind].agua_consumida_massa_step;

        gas_consumido_Mg = state.cells[ind].gas_consumido_massa_step;

    // Update the BSW
    double A_cross = state.cells[ind].duto.area;
    double Lcel    = state.cells[ind].dx;
    double Vlivre  = std::max(A_cross * Lcel - state.cells[ind].V_h, 1e-12);

    double frac_agua = std::max((1-state.cells[ind].alfR)*(1-state.cells[ind].betR)*state.cells[ind].FW, 1e-12);
    double frac_oleo = std::max((1-state.cells[ind].alfR)*(1-state.cells[ind].betR)*(1-state.cells[ind].FW), 1e-12);

    double Vagua = frac_agua * Vlivre;
    double Voil  = frac_oleo * Vlivre;

    double rho_w = std::max(state.cells[ind].flui.MasEspAgua(state.cells[ind].pres, state.cells[ind].temp), 1e-12);
    double Vagua_new = Vagua - agua_consumida_Mw / rho_w;
    if (Vagua_new < 0.0) Vagua_new = 0.0;

    double BSW_old = state.cells[ind].flui.BSW;
    double den = Voil + Vagua_new;
    if (den > 1e-12) {
    state.cells[ind].flui.BSW = Vagua_new / den;
    } else {
    state.cells[ind].flui.BSW = BSW_old;
    }
    //state.cells[ind].FW=state.cells[ind].flui.BSW;
    if (ind==3) cout << " t [s]: " << (*state.globals).lixo5 << " BSW: " << BSW_old << " FW: " << state.cells[ind].FW << " frac_agua: " << frac_agua << " BSW atualizada apos acoplamento " << state.cells[ind].flui.BSW << endl;
    //if (ind==3) system("pause");

    } // hydrate change

    if (state.input.calculaEnvelope==1 && state.input.tipoHmodel==3 && (*state.globals).lixo5>0.01) { // hydrate change

    agua_consumida_Mw  = state.cells[ind].agua_consumida_massa_step;

    gas_consumido_Mg   = state.cells[ind].gas_consumido_massa_step;

    // Update the BSW
    double A_cross = state.cells[ind].duto.area;
    double Lcel    = state.cells[ind].dx;
    double Vlivre  = std::max(A_cross * Lcel - state.cells[ind].V_h_total, 1e-12);

    double frac_agua = std::max((1-state.cells[ind].alfR)*(1-state.cells[ind].betR)*state.cells[ind].FW, 1e-12);
    double frac_oleo = std::max((1-state.cells[ind].alfR)*(1-state.cells[ind].betR)*(1-state.cells[ind].FW), 1e-12);

    double Vagua = frac_agua * Vlivre;
    double Voil  = frac_oleo * Vlivre;

    double rho_w = std::max(state.cells[ind].flui.MasEspAgua(state.cells[ind].pres, state.cells[ind].temp), 1e-12);
    double Vagua_new = Vagua - agua_consumida_Mw / rho_w;
    if (Vagua_new < 0.0) Vagua_new = 0.0;

    double BSW_old = state.cells[ind].flui.BSW;
    double den = Voil + Vagua_new;
    if (den > 1e-12) {
    state.cells[ind].flui.BSW = Vagua_new / den;
    } else {
    state.cells[ind].flui.BSW = BSW_old;
    }

    //if (ind==3) cout << " t [s]: " << (*state.globals).lixo5 << " BSW: " << BSW_old << " FW: " << state.cells[ind].FW << " frac_agua: " << frac_agua << " BSW atualizada apos acoplamento " << state.cells[ind].flui.BSW << endl;

    } // hydrate change
}

/// Adds to cell ind the mass its source delivers this step when the source is a
/// choke source (type 9 -- on the first iteration of a parallel network, on its
/// primary side, the flow recorded for that connection instead), a multiple
/// source (10) or a radial or 2D porous medium (15, 16).
void refreshChokeMultipleAndPorousSources(const SourceState &state, int ind) {
    if (state.cells[ind].acsr.tipo == kAccessoryLeak) {
        state.cells[ind].acsr.fontechk.fluidoP = state.cells[ind].flui;
        state.cells[ind].acsr.fontechk.presT = state.cells[ind].pres;
        state.cells[ind].acsr.fontechk.tempT = state.cells[ind].temp;
        double pres = state.cells[ind].pres;
        double temp = state.cells[ind].temp;
        double alf = state.cells[ind].alf;
        double bet = state.cells[ind].bet;
        double rhog = state.cells[ind].flui.MasEspGas(pres, temp);
        double rhoP = state.cells[ind].flui.MasEspLiq(pres, temp);
        double rhoC = state.cells[ind].fluicol.MasEspFlu(pres, temp);
        state.cells[ind].acsr.fontechk.titT = alf * rhog / (alf * rhog + (1 - alf) * (bet * rhoC + (1 - bet) * rhoP));
        state.cells[ind].acsr.fontechk.betIST = bet;
        if ((*state.globals).chaveRedeParalela == 0 || (*state.globals).iterRede > 0 || state.parallelSecondaryBranch == 1 || state.parallelSecondaryBoundaryCondition == 1) {
            state.cells[ind].acsr.fontechk.VMas();
            state.cells[ind].fontemassLR += state.cells[ind].acsr.fontechk.masP;
            state.cells[ind].fontemassCR += state.cells[ind].acsr.fontechk.masC;
            state.cells[ind].fontemassGR += state.cells[ind].acsr.fontechk.masG;
        } else {
            int nfonte = state.parallelSourceCells.size();
            int match = 0;
            int iconex;
            for (int ifonte = 0; ifonte < nfonte; ifonte++) {
                if (state.parallelSourceCells[ifonte] == ind) {
                    match = 1;
                    iconex = ifonte;
                    break;
                }
            }
            if (match == 0) {
                state.cells[ind].acsr.fontechk.VMas();
                state.cells[ind].fontemassLR += state.cells[ind].acsr.fontechk.masP;
                state.cells[ind].fontemassCR += state.cells[ind].acsr.fontechk.masC;
                state.cells[ind].fontemassGR += state.cells[ind].acsr.fontechk.masG;
            } else {
                state.cells[ind].fontemassLR += state.parallelSourceProductionLiquid[iconex];
                state.cells[ind].fontemassCR += state.parallelSourceComplementaryLiquid[iconex];
                state.cells[ind].fontemassGR += state.parallelSourceGas[iconex];
            }
        }
    }
    if (state.cells[ind].acsr.tipo == kAccessoryMultipleSource) {
        if (state.cells[ind].acsr.injm.condTermo == 1) {
            state.cells[ind].fontemassLR = state.cells[ind].acsr.injm.MassP;
            state.cells[ind].fontemassCR = state.cells[ind].acsr.injm.MassC;
            state.cells[ind].fontemassGR = state.cells[ind].acsr.injm.MassG;
        } else {
            state.cells[ind].fontemassCR = state.cells[ind].acsr.injm.MassC;
            double masT = state.cells[ind].acsr.injm.MassP + state.cells[ind].acsr.injm.MassG;
            double tit = state.cells[ind].acsr.injm.FluidoPro.FracMassHidra(state.cells[ind].pres, state.cells[ind].temp);
            state.cells[ind].acsr.injm.MassG = state.cells[ind].fontemassGR = masT * tit;
            state.cells[ind].acsr.injm.MassP = state.cells[ind].fontemassLR = masT - state.cells[ind].fontemassGR;
        }
    }
    if (state.cells[ind].acsr.tipo == kAccessoryRadialPorous) {
        state.cells[ind].acsr.radialPoro.pW.val[0] = state.cells[ind].pres;
        double rs;
        double bo;
        double ba;
        if (state.cells[ind].flui.RGO < 1e7) {
            rs = state.cells[ind + 1].flui.RS(state.cells[ind].pres, state.cells[ind].temp);
            bo = state.cells[ind + 1].flui.BOFunc(state.cells[ind].pres, state.cells[ind].temp, rs);
            ba = state.cells[ind + 1].flui.BAFunc(state.cells[ind].pres, state.cells[ind].temp);
            rs = rs * kBarrelPerCubicMetre / kCubicFootPerCubicMetre;
        } else {
            bo = 1;
            rs = 0;
            ba = 0.;
        }
        // in-situ BSW of the previous cell (in the march, cell i)
        double vfw = state.cells[ind + 1].flui.BSW * ba / (bo + ba * state.cells[ind + 1].flui.BSW - state.cells[ind + 1].flui.BSW * bo);
        state.cells[ind].acsr.radialPoro.sWPoc = vfw * (1. - state.cells[ind].acsr.radialPoro.satAconat) + state.cells[ind].acsr.radialPoro.satAconat;
        state.cells[ind].acsr.radialPoro.Pint = state.cells[ind].pres;
        if (state.steadyMode == 0) {
            state.cells[ind].acsr.radialPoro.avancoPressao();
            state.cells[ind].fontemassLR = state.cells[ind].acsr.radialPoro.fluxIni + state.cells[ind].acsr.radialPoro.fluxIniA;
            state.cells[ind].fontemassCR = 0.;
            state.cells[ind].fontemassGR = state.cells[ind].acsr.radialPoro.fluxIniG;
        } else {
            state.cells[ind].acsr.radialPoro.pseudoTrans();
            state.cells[ind].fontemassLR = state.cells[ind].acsr.radialPoro.fluxIni + state.cells[ind].acsr.radialPoro.fluxIniA;
            state.cells[ind].fontemassCR = 0.;
            state.cells[ind].fontemassGR = state.cells[ind].acsr.radialPoro.fluxIniG;
        }
    }
    if (state.cells[ind].acsr.tipo == kAccessoryPorous2D) {
        state.cells[ind].acsr.poroso2D.dados.pW.val[0] = state.cells[ind].pres;
        double rs;
        double bo;
        double ba;
        if (state.cells[ind].flui.RGO < 1e7) {
            rs = state.cells[ind + 1].flui.RS(state.cells[ind].pres, state.cells[ind].temp);
            bo = state.cells[ind + 1].flui.BOFunc(state.cells[ind].pres, state.cells[ind].temp, rs);
            ba = state.cells[ind + 1].flui.BAFunc(state.cells[ind].pres, state.cells[ind].temp);
            rs = rs * kBarrelPerCubicMetre / kCubicFootPerCubicMetre;
        } else {
            bo = 1;
            rs = 0;
            ba = 0.;
        }
        // in-situ BSW of the previous cell (in the march, cell i)
        double vfw = state.cells[ind + 1].flui.BSW * ba / (bo + ba * state.cells[ind + 1].flui.BSW - state.cells[ind + 1].flui.BSW * bo);
        state.cells[ind].acsr.poroso2D.sWPoc = vfw * (1. - state.cells[ind].acsr.poroso2D.dados.satAconat) + state.cells[ind].acsr.poroso2D.dados.satAconat;
        state.cells[ind].acsr.poroso2D.dados.transfer.sWPoc = state.cells[ind].acsr.poroso2D.sWPoc;
        state.cells[ind].acsr.poroso2D.dados.pInt = state.cells[ind].pres;
        state.cells[ind].acsr.poroso2D.dados.transfer.Pint = state.cells[ind].pres;
        if (state.steadyMode == 0) {
            state.cells[ind].acsr.poroso2D.avancoPressao();
            state.cells[ind].fontemassLR = state.cells[ind].acsr.poroso2D.dados.transfer.fluxIni + state.cells[ind].acsr.poroso2D.dados.transfer.fluxIniA;
            state.cells[ind].fontemassCR = 0.;
            state.cells[ind].fontemassGR = state.cells[ind].acsr.poroso2D.dados.transfer.fluxIniG;
        } else {
            state.cells[ind].acsr.poroso2D.pseudoTransientePoroso();
            state.cells[ind].fontemassLR = state.cells[ind].acsr.poroso2D.dados.transfer.fluxIni + state.cells[ind].acsr.poroso2D.dados.transfer.fluxIniA;
            state.cells[ind].fontemassCR = 0.;
            state.cells[ind].fontemassGR = state.cells[ind].acsr.poroso2D.dados.transfer.fluxIniG;
        }
    }
}

}  // namespace

void addMasterValveFlow(const SourceState &state, int ind) {
    if (state.cells[ind].acsr.chk.AreaGarg < (1e-3 + state.input.master1.razareaativ) * state.cells[ind].duto.area && state.cells[ind].acsr.chk.AreaGarg > 1e-3 * state.cells[ind].duto.area) {
        double tE = state.cells[ind].temp;
        double alfE = state.cells[ind].alf;
        double betE = state.cells[ind].bet;
        double sense = 1.;

        double maxSup = 0.;

        double rholp = state.cells[ind].flui.MasEspLiq(state.cells[ind].pres, state.cells[ind].temp);
        double rholc = state.cells[ind].fluicol.MasEspFlu(state.cells[ind].pres, state.cells[ind].temp);
        double rholmix = (1 - betE) * rholp + betE * rholc;

        double alfJ = state.cells[ind + 1].alf;
        double betJ = state.cells[ind + 1].bet;
        double rholpJ = state.cells[ind + 1].flui.MasEspLiq(state.cells[ind + 1].pres, state.cells[ind + 1].temp);
        double rholcJ = state.cells[ind + 1].fluicol.MasEspFlu(state.cells[ind + 1].pres, state.cells[ind + 1].temp);
        double rholmixJ = (1 - betJ) * rholpJ + betJ * rholcJ;

        double hidroM = sin(state.cells[ind].duto.teta) * (0.5 * state.cells[ind].dx) * (rholmix * (1 - alfE) + alfE * state.cells[ind].flui.MasEspGas(state.cells[ind].pres, state.cells[ind].temp)) * kGravity / kPascalPerKgfPerCm2Coarse;
        double hidroJ = sin(state.cells[ind + 1].duto.teta) * (0.5 * state.cells[ind + 1].dx) * (rholmixJ * (1 - alfJ) + alfJ * state.cells[ind + 1].flui.MasEspGas(state.cells[ind + 1].pres, state.cells[ind + 1].temp)) * kGravity / kPascalPerKgfPerCm2Coarse;

        double masentrada = state.cells[ind].MC;
        double massgas = state.cells[ind].MC - state.cells[ind].Mliqini;
        double tit;
        if (fabs(masentrada < 1e-9) && fabs(massgas < 1e-9)) {
            tit = alfE * state.cells[ind].flui.MasEspGas(state.cells[ind].pres, state.cells[ind].temp) / (state.cells[ind].flui.MasEspGas(state.cells[ind].pres, state.cells[ind].temp) * alfE + rholmix * (1. - alfE));
        } else {
            if ((massgas >= 0 && state.cells[ind].Mliqini <= 0) || (massgas < 0 && state.cells[ind].Mliqini == 0))
                tit = 1.;
            else if (massgas <= 0 && state.cells[ind].Mliqini > 0)
                tit = 0.;
            else if (masentrada < 0)
                tit = 1.;
            else
                tit = fabs(massgas / masentrada);
            if (tit > 1)
                tit = 1;
        }

        double masChk;

        double ypres = (state.cells[ind + 1].pres + hidroJ) / (state.cells[ind].pres - hidroM);
        int check = 1;
        if (ypres < 1. || check == 1) {
            if ((*state.globals).lixo5 > 7059) {
                int para;
                para = 1;
            }
            masChk = state.cells[ind].acsr.chk.vazmassSachd(ypres, state.cells[ind].pres - hidroM, tE, alfE,
                                                       betE, tit, state.cells[ind].flui, state.cells[ind].fluicol);
            maxSup = state.cells[ind].acsr.chk.vazmaxSachd(state.cells[ind].pres - hidroM, tE, alfE,
                                                      betE, tit, state.cells[ind].flui, state.cells[ind].fluicol);
        } else {
            tE = state.cells[ind + 1].temp;
            alfE = state.cells[ind + 1].alf;
            betE = state.cells[ind + 1].bet;
            ypres = 1. / ypres;
            sense = -1;
            if ((*state.globals).lixo5 > 7059) {
                int para;
                para = 1;
            }
            masChk = state.cells[ind].acsr.chk.vazmassSachd(ypres, state.cells[ind + 1].pres + hidroJ, tE, alfE, betE, tit, state.cells[ind + 1].flui,
                                                       state.cells[ind + 1].fluicol);
            maxSup = state.cells[ind].acsr.chk.vazmaxSachd(state.cells[ind + 1].pres + hidroJ, tE, alfE, betE, tit, state.cells[ind + 1].flui,
                                                      state.cells[ind + 1].fluicol);
        }

        if (fabs(ypres) > fabs(state.cells[ind].acsr.chk.razpres))
            maxSup = masChk;

        if (state.cells[ind].acsr.chk.AreaGarg < (1e-3) * state.cells[ind].duto.area || (state.cells[ind].pres < state.cells[ind + 1].pres && check == 1))
            maxSup =
                0.;

        double masliq;
        double masgas;
        masliq = sense * maxSup * (1. - tit);
        masgas = sense * maxSup * tit;

        state.cells[ind].fontemassLR += (-1 * masliq * (1 - betE) * rholp / rholmix);
        state.cells[ind].fontemassCR += (-1 * masliq * betE * rholc / rholmix);
        state.cells[ind].fontemassGR += (-1 * masgas);

        state.cells[ind + 1].fontemassLR += (masliq * (1 - betE) * rholp / rholmix);
        state.cells[ind + 1].fontemassCR += (masliq * betE * rholc / rholmix);
        state.cells[ind + 1].fontemassGR += (masgas);
    }
}

void renewSourceTerms(const SourceState &state, int ind) {

    double pr = state.cells[ind].pres;
    double tr = state.cells[ind].temp;
    if (ind == 0) {
        state.cells[ind].fontemassLR = 0.;
        state.cells[ind].fontemassCR = 0.;
        state.cells[ind].fontemassGR = 0.;
    } else {
        if (state.cells[ind - 1].acsr.tipo != kAccessoryChoke && state.cells[ind - 1].acsr.tipo != kAccessoryVolumetricPump && (ind < state.lastCell || state.cells[state.lastCell].acsr.tipo == kAccessoryInflowPerformance)) {
            state.cells[ind].fontemassLR = 0.;
            state.cells[ind].fontemassCR = 0.;
            state.cells[ind].fontemassGR = 0.;
        }
    }

    double agua_consumida_Mw = 0.;
    double gas_consumido_Mg = 0.;

    consumeHydrateFormationMass(state, gas_consumido_Mg, agua_consumida_Mw, ind);
    if (state.cells[ind].acsr.tipo == kAccessoryGasInjection) {
        if (state.cells[ind].acsr.injg.tipoflu == 0) {
            double masgas = state.cells[ind].acsr.injg.VMas(pr, tr);
            if (fabs(masgas) < (*state.globals).localtiny)
                masgas = 0.;
            if (state.cells[ind].acsr.injg.seco == 1) {
                state.cells[ind].fontemassGR += masgas;
                state.cells[ind].fontemassLR = 0.;
                state.cells[ind].fontemassCR = 0.;
            } else {
                double tit;
                if (state.input.flashCompleto != 2)
                    tit = state.cells[ind].acsr.injg.FluidoPro.FracMassHidra(1., 20.);
                else
                    tit = state.cells[ind].acsr.injg.FluidoPro.dStockTankVaporMassFraction;
                double masT = masgas / tit;
                tit = state.cells[ind].acsr.injg.FluidoPro.FracMassHidra(pr, tr);
                state.cells[ind].fontemassGR += masT * tit;
                state.cells[ind].fontemassLR += masT * (1. - tit);
                double rcomp = state.cells[ind].acsr.injg.fluidocol.MasEspFlu(1., 20.);
                state.cells[ind].fontemassCR += rcomp * state.cells[ind].acsr.injg.razCompGas * state.cells[ind].acsr.injg.QGas / kSecondsPerDay;
            }
        } else {
            state.cells[ind].fontemassGR += 0.;
            state.cells[ind].fontemassLR = 0.;
            state.cells[ind].fontemassCR = state.cells[ind].acsr.injg.QGas;
        }
    }
    if (state.cells[ind].acsr.tipo == kAccessoryLiquidInjection) {
        double rlcA = state.cells[ind].acsr.injl.fluidocol.MasEspFlu(1.001, 15.);
        state.cells[ind].fontemassCR += rlcA * state.cells[ind].acsr.injl.QLiq * state.cells[ind].acsr.injl.bet / kSecondsPerDay;
        double massic = state.cells[ind].acsr.injl.QLiq * (1. - state.cells[ind].acsr.injl.bet) / kSecondsPerDay;
        double Rhogs = state.cells[ind].acsr.injl.FluidoPro.Deng * kAirDensityAtStandardConditions; // cel[ind].acsr.injl.FluidoPro.MasEspGas(1, 15);
        double Rhols = (1000 * 141.5 / (131.5 + state.cells[ind].acsr.injl.FluidoPro.API)) * (1 - state.cells[ind].acsr.injl.FluidoPro.BSW) + 1000. * state.cells[ind].acsr.injl.FluidoPro.Denag * state.cells[ind].acsr.injl.FluidoPro.BSW;
        double multiplicador = (Rhols + state.cells[ind].acsr.injl.FluidoPro.RGO * Rhogs * (1 - state.cells[ind].acsr.injl.FluidoPro.BSW));
        massic *= multiplicador;
        double fracmasshidra = state.cells[ind].acsr.injl.FluidoPro.FracMassHidra(pr, tr);
        state.cells[ind].fontemassLR += (1. - fracmasshidra) * massic;
        state.cells[ind].fontemassGR += fracmasshidra * massic;
    }
    if (state.cells[ind].acsr.tipo == kAccessoryInflowPerformance) {
        if (pr < state.cells[ind].acsr.ipr.Pres) {
            state.cells[ind].fontemassLR += state.cells[ind].acsr.ipr.MasL(pr, tr);
            state.cells[ind].fontemassCR = 0.;
            state.cells[ind].fontemassGR += state.cells[ind].acsr.ipr.MasG(pr, tr);
        } else {
            double tit;
            tit = state.cells[ind].alf * state.cells[ind].flui.MasEspGas(pr, tr) /
                  (state.cells[ind].alf * state.cells[ind].flui.MasEspGas(pr, tr) +
                   (1. - state.cells[ind].alf) * (1. - state.cells[ind].bet) * state.cells[ind].flui.MasEspLiq(pr, tr) +
                   (1. - state.cells[ind].alf) * state.cells[ind].bet * state.cells[ind].fluicol.MasEspFlu(pr, tr));
            state.cells[ind].fontemassLR += (1. - state.cells[ind].alf) * (1. - state.cells[ind].bet) * state.cells[ind].acsr.ipr.VMas(pr, tr) * state.cells[ind].flui.MasEspLiq(pr, tr);
            state.cells[ind].fontemassCR += (1. - state.cells[ind].alf) * state.cells[ind].bet * state.cells[ind].acsr.ipr.VMas(pr, tr) * state.cells[ind].fluicol.MasEspFlu(pr, tr);
            state.cells[ind].fontemassGR += 1 * state.cells[ind].alf * state.cells[ind].acsr.ipr.VMas(pr, tr) * state.cells[ind].flui.MasEspGas(pr, tr);
            state.cells[ind].acsr.ipr.deriP *= (1. - state.cells[ind].alf) * (1. - state.cells[ind].bet) * state.cells[ind].flui.MasEspLiq(pr, tr);
            state.cells[ind].acsr.ipr.deriC *= (1. - state.cells[ind].alf) * state.cells[ind].bet * state.cells[ind].fluicol.MasEspFlu(pr, tr);
            state.cells[ind].acsr.ipr.deriG *= 1 * state.cells[ind].alf * state.cells[ind].flui.MasEspGas(pr, tr);
        }
    }
    if (state.cells[ind].acsr.tipo == kAccessoryChoke) {
        state.cells[ind + 1].fontemassLR = 0.;
        state.cells[ind + 1].fontemassCR = 0.;
        state.cells[ind + 1].fontemassGR = 0.;
        if (state.steadyMode == 0)
            addMasterValveFlow(state, ind);
    }
    if (state.cells[ind].acsr.tipo == kAccessoryVolumetricPump) {
        state.cells[ind + 1].fontemassLR = 0.;
        state.cells[ind + 1].fontemassCR = 0.;
        state.cells[ind + 1].fontemassGR = 0.;
        state.cells[ind].acsr.bvol.fluido = state.cells[ind].flui;
        state.cells[ind].acsr.bvol.fluicol = state.cells[ind].fluicol;
        if (fabs(state.cells[ind].acsr.bvol.freq) > 1) {
            double alfM = state.cells[ind].alf;
            double betM = state.cells[ind].bet;
            double presM = state.cells[ind].pres;
            double tempM = state.cells[ind].temp;
            double presM1 = state.cells[ind].presR;
            double tempM1 = state.cells[ind].tempR;
            state.cells[ind].acsr.bvol.vazmass(presM, tempM, presM1, tempM1, betM, alfM);
            state.cells[ind].fontemassLR -= state.cells[ind].acsr.bvol.MLiqP;
            state.cells[ind].fontemassCR -= state.cells[ind].acsr.bvol.MLiqC;
            state.cells[ind].fontemassGR -= state.cells[ind].acsr.bvol.MGas;
            state.cells[ind + 1].fontemassLR += state.cells[ind].acsr.bvol.MLiqP;
            state.cells[ind + 1].fontemassCR += state.cells[ind].acsr.bvol.MLiqC;
            state.cells[ind + 1].fontemassGR += state.cells[ind].acsr.bvol.MGas;
        }
    }
    refreshChokeMultipleAndPorousSources(state, ind);

    if (state.input.calculaEnvelope == 1 && state.input.tipoHmodel == 2 && state.cells[ind].flui.BSW > 1e-12 && (*state.globals).lixo5 > 0.01) {
        state.cells[ind].fontemassLR -= (agua_consumida_Mw / (*state.globals).lixo5);
        state.cells[ind].fontemassGR -= (gas_consumido_Mg / (*state.globals).lixo5);
    }
  	// hydrate change 3
  	if (state.input.calculaEnvelope==1 && state.input.tipoHmodel==3 && state.cells[ind].flui.BSW>1e-14 && (*state.globals).lixo5>0.01) { // hydrate change
  		state.cells[ind].fontemassLR -= (agua_consumida_Mw / (*state.globals).lixo5);
  		state.cells[ind].fontemassGR -= (gas_consumido_Mg / (*state.globals).lixo5);
  	}
}

}  // namespace sisprod::sources
