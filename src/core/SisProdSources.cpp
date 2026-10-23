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

void refreshLeakSource(const SourceState &state, int cellIndex) {
    Cel &cell = state.cells[cellIndex];
    auto &source = cell.acsr.fontechk;
    source.fluidoP = cell.flui;
    source.presT = cell.pres;
    source.tempT = cell.temp;
    double pres = cell.pres;
    double temp = cell.temp;
    double alf = cell.alf;
    double bet = cell.bet;
    double rhog = cell.flui.MasEspGas(pres, temp);
    double rhoP = cell.flui.MasEspLiq(pres, temp);
    double rhoC = cell.fluicol.MasEspFlu(pres, temp);
    source.titT = alf * rhog / (alf * rhog + (1 - alf) * (bet * rhoC + (1 - bet) * rhoP));
    source.betIST = bet;
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
        source.VMas();
        cell.fontemassLR += source.masP;
        cell.fontemassCR += source.masC;
        cell.fontemassGR += source.masG;
    } else {
        cell.fontemassLR += state.parallelSourceProductionLiquid[iconex];
        cell.fontemassCR += state.parallelSourceComplementaryLiquid[iconex];
        cell.fontemassGR += state.parallelSourceGas[iconex];
    }
}

void advanceRadialPorousSource(const SourceState &state, int cellIndex) {
    Cel &cell = state.cells[cellIndex];
    auto &source = cell.acsr.radialPoro;
    source.pW.val[0] = cell.pres;
    // Use the next array cell's fluid at this cell's pressure and temperature.
    double vfw = inSituWaterFraction(cell, state.cells[cellIndex + 1]);
    source.sWPoc = vfw * (1. - source.satAconat) + source.satAconat;
    source.Pint = cell.pres;
    if (state.steadyMode == 0)
        source.avancoPressao();
    else
        source.pseudoTrans();
    cell.fontemassLR = source.fluxIni + source.fluxIniA;
    cell.fontemassCR = 0.;
    cell.fontemassGR = source.fluxIniG;
}

void advancePorous2DSource(const SourceState &state, int cellIndex) {
    Cel &cell = state.cells[cellIndex];
    auto &source = cell.acsr.poroso2D;
    source.dados.pW.val[0] = cell.pres;
    // Use the next array cell's fluid at this cell's pressure and temperature.
    double vfw = inSituWaterFraction(cell, state.cells[cellIndex + 1]);
    source.sWPoc = vfw * (1. - source.dados.satAconat) + source.dados.satAconat;
    source.dados.transfer.sWPoc = source.sWPoc;
    source.dados.pInt = cell.pres;
    source.dados.transfer.Pint = cell.pres;
    if (state.steadyMode == 0)
        source.avancoPressao();
    else
        source.pseudoTransientePoroso();
    cell.fontemassLR = source.dados.transfer.fluxIni + source.dados.transfer.fluxIniA;
    cell.fontemassCR = 0.;
    cell.fontemassGR = source.dados.transfer.fluxIniG;
}

void refreshChokeMultipleAndPorousSources(const SourceState &state, int cellIndex) {
    if (state.cells[cellIndex].acsr.tipo == kAccessoryLeak) {
        refreshLeakSource(state, cellIndex);
    }
    if (state.cells[cellIndex].acsr.tipo == kAccessoryMultipleSource) {
        Cel &cell = state.cells[cellIndex];
        if (cell.acsr.injm.condTermo == 1) {
            cell.fontemassLR = cell.acsr.injm.MassP;
            cell.fontemassCR = cell.acsr.injm.MassC;
            cell.fontemassGR = cell.acsr.injm.MassG;
        } else {
            cell.fontemassCR = cell.acsr.injm.MassC;
            double masT = cell.acsr.injm.MassP + cell.acsr.injm.MassG;
            double tit = cell.acsr.injm.FluidoPro.FracMassHidra(cell.pres, cell.temp);
            cell.acsr.injm.MassG = cell.fontemassGR = masT * tit;
            cell.acsr.injm.MassP = cell.fontemassLR = masT - cell.fontemassGR;
        }
    }
    if (state.cells[cellIndex].acsr.tipo == kAccessoryRadialPorous) {
        advanceRadialPorousSource(state, cellIndex);
    }
    if (state.cells[cellIndex].acsr.tipo == kAccessoryPorous2D) {
        advancePorous2DSource(state, cellIndex);
    }
}

void applyGasInjectionSource(const SourceState &state, int cellIndex, double pr, double tr) {
    Cel &cell = state.cells[cellIndex];
    auto &source = cell.acsr.injg;
    if (source.tipoflu == 0) {
        double masgas = source.VMas(pr, tr);
        if (fabs(masgas) < (*state.globals).localtiny)
            masgas = 0.;
        if (source.seco == 1) {
            cell.fontemassGR += masgas;
            cell.fontemassLR = 0.;
            cell.fontemassCR = 0.;
        } else {
            double tit;
            if (state.input.flashCompleto != 2)
                tit = source.FluidoPro.FracMassHidra(1., 20.);
            else
                tit = source.FluidoPro.dStockTankVaporMassFraction;
            double masT = masgas / tit;
            tit = source.FluidoPro.FracMassHidra(pr, tr);
            cell.fontemassGR += masT * tit;
            cell.fontemassLR += masT * (1. - tit);
            double rcomp = source.fluidocol.standardDensity();
            cell.fontemassCR += rcomp * source.razCompGas * source.QGas / units::kSecondsPerDay;
        }
    } else {
        cell.fontemassGR += 0.;
        cell.fontemassLR = 0.;
        cell.fontemassCR = source.QGas;
    }
}

void applyLiquidInjectionSource(const SourceState &state, int cellIndex, double pr, double tr) {
    Cel &cell = state.cells[cellIndex];
    auto &source = cell.acsr.injl;
    double rlcA = source.fluidocol.standardDensity();
    cell.fontemassCR += rlcA * source.QLiq * source.bet / units::kSecondsPerDay;
    double massic = source.QLiq * (1. - source.bet) / units::kSecondsPerDay;
    double Rhogs = source.FluidoPro.Deng * kAirDensityAtStandardConditions;
    double Rhols = (1000 * 141.5 / (131.5 + source.FluidoPro.API)) * (1 - source.FluidoPro.BSW) + 1000. * source.FluidoPro.Denag * source.FluidoPro.BSW;
    double multiplicador = (Rhols + source.FluidoPro.RGO * Rhogs * (1 - source.FluidoPro.BSW));
    massic *= multiplicador;
    double fracmasshidra = source.FluidoPro.FracMassHidra(pr, tr);
    cell.fontemassLR += (1. - fracmasshidra) * massic;
    cell.fontemassGR += fracmasshidra * massic;
}

void applyInflowPerformanceSource(const SourceState &state, int cellIndex, double pr, double tr) {
    Cel &cell = state.cells[cellIndex];
    auto &source = cell.acsr.ipr;
    if (pr < source.Pres) {
        cell.fontemassLR += source.MasL(pr, tr);
        cell.fontemassCR = 0.;
        cell.fontemassGR += source.MasG(pr, tr);
    } else {
        double tit;
        tit = cell.alf * cell.flui.MasEspGas(pr, tr) /
              (cell.alf * cell.flui.MasEspGas(pr, tr) +
               (1. - cell.alf) * (1. - cell.bet) * cell.flui.MasEspLiq(pr, tr) +
               (1. - cell.alf) * cell.bet * cell.fluicol.MasEspFlu(pr, tr));
        cell.fontemassLR += (1. - cell.alf) * (1. - cell.bet) * source.VMas(pr, tr) * cell.flui.MasEspLiq(pr, tr);
        cell.fontemassCR += (1. - cell.alf) * cell.bet * source.VMas(pr, tr) * cell.fluicol.MasEspFlu(pr, tr);
        cell.fontemassGR += 1 * cell.alf * source.VMas(pr, tr) * cell.flui.MasEspGas(pr, tr);
        source.deriP *= (1. - cell.alf) * (1. - cell.bet) * cell.flui.MasEspLiq(pr, tr);
        source.deriC *= (1. - cell.alf) * cell.bet * cell.fluicol.MasEspFlu(pr, tr);
        source.deriG *= 1 * cell.alf * cell.flui.MasEspGas(pr, tr);
    }
}

void applyVolumetricPumpSource(const SourceState &state, int cellIndex) {
    Cel &cell = state.cells[cellIndex];
    auto &source = cell.acsr.bvol;
    state.cells[cellIndex + 1].fontemassLR = 0.;
    state.cells[cellIndex + 1].fontemassCR = 0.;
    state.cells[cellIndex + 1].fontemassGR = 0.;
    source.fluido = cell.flui;
    source.fluicol = cell.fluicol;
    if (fabs(source.freq) > 1) {
        double alfM = cell.alf;
        double betM = cell.bet;
        double presM = cell.pres;
        double tempM = cell.temp;
        double presM1 = cell.presR;
        double tempM1 = cell.tempR;
        source.vazmass(presM, tempM, presM1, tempM1, betM, alfM);
        cell.fontemassLR -= source.MLiqP;
        cell.fontemassCR -= source.MLiqC;
        cell.fontemassGR -= source.MGas;
        state.cells[cellIndex + 1].fontemassLR += source.MLiqP;
        state.cells[cellIndex + 1].fontemassCR += source.MLiqC;
        state.cells[cellIndex + 1].fontemassGR += source.MGas;
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
    if (cellIndex == 0 ||
        (state.cells[cellIndex - 1].acsr.tipo != kAccessoryChoke &&
         state.cells[cellIndex - 1].acsr.tipo != kAccessoryVolumetricPump &&
         (cellIndex < state.lastCell || state.cells[state.lastCell].acsr.tipo == kAccessoryInflowPerformance))) {
        state.cells[cellIndex].fontemassLR = 0.;
        state.cells[cellIndex].fontemassCR = 0.;
        state.cells[cellIndex].fontemassGR = 0.;
    }

    double agua_consumida_Mw = 0.;
    double gas_consumido_Mg = 0.;

    consumeHydrateFormationMass(state, gas_consumido_Mg, agua_consumida_Mw, cellIndex);
    if (state.cells[cellIndex].acsr.tipo == kAccessoryGasInjection) {
        applyGasInjectionSource(state, cellIndex, pr, tr);
    }
    if (state.cells[cellIndex].acsr.tipo == kAccessoryLiquidInjection) {
        applyLiquidInjectionSource(state, cellIndex, pr, tr);
    }
    if (state.cells[cellIndex].acsr.tipo == kAccessoryInflowPerformance) {
        applyInflowPerformanceSource(state, cellIndex, pr, tr);
    }
    if (state.cells[cellIndex].acsr.tipo == kAccessoryChoke) {
        state.cells[cellIndex + 1].fontemassLR = 0.;
        state.cells[cellIndex + 1].fontemassCR = 0.;
        state.cells[cellIndex + 1].fontemassGR = 0.;
        if (state.steadyMode == 0)
            addMasterValveFlow(state, cellIndex);
    }
    if (state.cells[cellIndex].acsr.tipo == kAccessoryVolumetricPump) {
        applyVolumetricPumpSource(state, cellIndex);
    }
    refreshChokeMultipleAndPorousSources(state, cellIndex);

    const bool hydrateModel = state.input.calculaEnvelope == 1 && (state.input.tipoHmodel == 2 || state.input.tipoHmodel == 3);
    if (hydrateModel && state.cells[cellIndex].flui.BSW > 1e-12 && (*state.globals).lixo5 > 0.01) {
        state.cells[cellIndex].fontemassLR -= (agua_consumida_Mw / (*state.globals).lixo5);
        state.cells[cellIndex].fontemassGR -= (gas_consumido_Mg / (*state.globals).lixo5);
    }
}

}  // namespace sisprod::sources
