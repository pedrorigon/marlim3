#include "SisProdThermal.h"

#include "Leitura.h"
#include "celula3.h"

#include <math.h>

namespace sisprod::thermal {

double interpolateLatentHeat(const ThermalState &state, double pres, double temp) {
    int ndiv = state.input.tabent.npont - 1;
    int ipres = 0.;
    int itemp = 0.;
    int ipmarcador;
    int itmarcador;
    double latt;
    if (pres < state.latentHeatTable[1][0] || pres >= state.latentHeatTable[ndiv + 1][0] || temp < state.latentHeatTable[0][1] || temp >= state.latentHeatTable[0][ndiv + 1])
        latt = 0.;

    else {
        int e, m, d;
        e = 1;
        d = ndiv + 1;
        while (e <= d) {
            m = (e + d) / 2;
            ipmarcador = m;
            if (m == 1) {
                ipres = m;
                break;
            } else if (m == ndiv + 1 && state.latentHeatTable[m][0] == pres) {
                ipres = m - 1;
                break;
            }
            if (state.latentHeatTable[m][0] > pres && state.latentHeatTable[m - 1][0] <= pres) {
                ipres = m - 1;
                break;
            }
            if (state.latentHeatTable[m][0] < pres)
                e = m + 1;
            else
                d = m - 1;
        }
        e = 1;
        d = ndiv + 1;
        while (e <= d) {
            m = (e + d) / 2;
            itmarcador = m;
            if (m == 1) {
                itemp = m;
                break;
            } else if (m == ndiv + 1 && state.latentHeatTable[0][m] == temp) {
                itemp = m - 1;
                break;
            }
            if (state.latentHeatTable[0][m] > temp && state.latentHeatTable[0][m - 1] <= temp) {
                itemp = m - 1;
                break;
            }
            if (state.latentHeatTable[0][m - 1] < temp)
                e = m + 1;
            else
                d = m - 1;
        }
        double razpres = (state.latentHeatTable[ipres][0] - pres) / (state.latentHeatTable[ipres][0] - state.latentHeatTable[ipres + 1][0]);
        double raztemp = (state.latentHeatTable[0][itemp] - temp) / (state.latentHeatTable[0][itemp] - state.latentHeatTable[0][itemp + 1]);
        double latp1 = (1 - razpres) * (state.latentHeatTable[ipres][itemp]) + razpres * (state.latentHeatTable[ipres + 1][itemp]);
        double latp2 = (1 - razpres) * (state.latentHeatTable[ipres][itemp + 1]) + razpres * (state.latentHeatTable[ipres + 1][itemp + 1]);
        latt = (1 - raztemp) * latp1 + raztemp * latp2;
    }
    return latt;
}

double computeMixtureEnthalpy(const ThermalState &state, int i) {

    double dx = state.cells[i].dx;
    double dxmed = 0.5 * (state.cells[i].dx + state.cells[i - 1].dx);
    double dia = state.cells[i].duto.a;
    double area = 0.25 * M_PI * dia * dia;
    double alfmed = state.cells[i].alf;
    double betmed = state.cells[i].bet;
    double alfmed0 = state.cells[i].alfini;
    double betmed0 = state.cells[i].betini;
    double pmed = state.cells[i].pres;
    double pmed0 = state.cells[i].presini;
    double tmed = state.cells[i].temp;
    double ugsL;
    double ulsL;
    double ugsR;
    double ulsR;

    ugsL = state.cells[i].QG / area;
    ulsL = state.cells[i].QL / area;
    ugsR = state.cells[i + 1].QG / area;
    ulsR = state.cells[i + 1].QL / area;

    double betL = state.cells[i].bet;
    double betR = state.cells[i].bet;
    if (ugsL > 0.)
        betL = state.cells[i - 1].bet;
    if (ugsL < 0.)
        betL = state.cells[i + 1].bet;

    double ugsmed = 0.5 * (ugsL + ugsR);
    double ulsmed = 0.5 * (ulsL + ulsR);

    double presL = state.cells[i].presaux;
    double presR = state.cells[i + 1].presaux;
    if ((state.cells[i].acsr.tipo == 5 && state.cells[i].acsr.chk.AreaGarg <= (1e-3 + state.input.master1.razareaativ) * state.cells[i].duto.area) ||
        (state.cells[i].acsr.tipo == 4 && state.cells[i].acsr.bcs.freq > 0) ||
        (state.cells[i].acsr.tipo == 8 && state.cells[i].acsr.bvol.freq > 0.) ||
        (state.cells[i].acsr.tipo == 7 && fabs(state.cells[i].acsr.delp) > 0.) ||
        (state.cells[i].acsr.tipo == 17 && state.cells[i].acsr.multibcs.freq > 0)) {
        presR = state.cells[i].pres + (state.cells[i].pres - state.cells[i].presaux) * 0.5;
    }

    double tmedL = state.cells[i].tempL;
    if (state.cells[i].VTemper < 0.)
        tmedL = state.cells[i].temp;
    double tmedR = state.cells[i].temp;
    if (state.cells[i + 1].VTemper < 0.)
        tmedR = state.cells[i].tempR;

    double rhog = state.cells[i].flui.MasEspGas(pmed0, tmed);
    double rhop = state.cells[i].flui.MasEspLiq(pmed0, tmed);
    double rhoc = state.cells[i].fluicol.MasEspFlu(pmed0, tmed);
    double rhog1 = state.cells[i].flui.MasEspGas(pmed, tmed);
    double rhop1 = state.cells[i].flui.MasEspLiq(pmed, tmed);
    double hg = state.cells[i].flui.EntalpGas(pmed0, tmed);
    double hp = state.cells[i].flui.EntalpLiq(pmed0, tmed);
    double hc = state.cells[i].fluicol.CalorLiq(pmed0, tmed) * tmed; // corrigir entalpia

    double rhogL = state.cells[i].flui.MasEspGas(presL, tmedL);
    double rhopL = state.cells[i].flui.MasEspLiq(presL, tmedL);
    double rhocL = state.cells[i].fluicol.MasEspFlu(presL, tmedL);
    double hgL = state.cells[i].flui.EntalpGas(presL, tmedL);
    double hpL = state.cells[i].flui.EntalpLiq(presL, tmedL);
    double hcL = state.cells[i].fluicol.CalorLiq(presL, tmedL) * tmedL; // corrigir entalpia

    double rhogR = state.cells[i].flui.MasEspGas(presR, tmedR);
    double rhopR = state.cells[i].flui.MasEspLiq(presR, tmedR);
    double rhocR = state.cells[i].fluicol.MasEspFlu(presR, tmedR);
    double hgR = state.cells[i].flui.EntalpGas(presR, tmedR);
    double hpR = state.cells[i].flui.EntalpLiq(presR, tmedR);
    double hcR = state.cells[i].fluicol.CalorLiq(presR, tmedR) * tmedR; // corrigir entalpia

    double energintmixT0 = rhog * alfmed0 * (hg - pmed0 * 98066.5 / rhog) + rhop * (1 - betmed0) * (1. - alfmed0) * (hp - pmed0 * 98066.5 / rhop) +
                           betmed0 * rhoc * (1. - alfmed0) * (hc - pmed0 * 98066.5 / rhoc);
    double fluxHL = rhogL * ugsL * hgL + rhopL * (1. - betL) * ulsL * hpL + rhocL * betL * ulsL * hcL;
    double fluxHR = rhogR * ugsR * hgR + rhopR * (1. - betR) * ulsR * hpR + rhocR * betR * ulsR * hcR;
    double delFlux = (fluxHR - fluxHL) / dx;
    double hidro = (rhog1 * ugsmed + (1 - betmed) * ulsmed * rhop1 + betmed * ulsmed) * sin(state.cells[i].duto.teta) * 9.82;

    double fontemassG = 0.;
    double fontemassL = 0.;
    double fontemassC = 0.;
    double tfonte = state.cells[i].temp;
    double hgF;
    double hlF;
    double hcF = 0.;

    if (state.cells[i].acsr.tipo == 1) {
        tfonte = state.cells[i].acsr.injg.temp;
        hgF = state.cells[i].acsr.injg.FluidoPro.EntalpGas(state.cells[i].pres, tfonte);
        hlF = 0;
        hcF = 0;
    } else if (state.cells[i].acsr.tipo == 2) {
        tfonte = state.cells[i].acsr.injl.temp;
        hgF = state.cells[i].acsr.injl.FluidoPro.EntalpGas(state.cells[i].pres, tfonte);
        hlF = state.cells[i].acsr.injl.FluidoPro.EntalpLiq(state.cells[i].pres, tfonte);
        hcF = state.cells[i].acsr.injl.fluidocol.CalorLiq(state.cells[i].pres, tfonte) * tfonte
            /*entalpia fluido complementar a ser corrigida*/;
    } else if (state.cells[i].acsr.tipo == 3) {
        tfonte = state.cells[i].acsr.ipr.Tres;
        hgF = state.cells[i].acsr.ipr.FluidoPro.EntalpGas(state.cells[i].pres, tfonte);
        hlF = state.cells[i].acsr.ipr.FluidoPro.EntalpLiq(state.cells[i].pres, tfonte);
        hcF = 0;
    } else if (state.cells[i].acsrL != 0) {
        if ((*state.cells[i].acsrL).tipo == 5) {
            if ((*state.cells[i].acsrL).chk.AreaGarg < state.input.master1.razareaativ * state.cells[i].dutoL.area && (*state.cells[i].acsrL).chk.AreaGarg > 1e-5 * state.cells[i].dutoL.area) {
                double tE = state.cells[i - 1].temp;
                double alfE = state.cells[i - 1].alf;
                double betE = state.cells[i - 1].bet;

                double rholp = state.cells[i - 1].flui.MasEspLiq(state.cells[i - 1].pres, state.cells[i - 1].temp);
                double rholc = state.cells[i - 1].fluicol.MasEspFlu(state.cells[i - 1].pres, state.cells[i - 1].temp);
                double rholmix = (1 - betE) * rholp + betE * rholc;

                double alfJ = state.cells[i].alf;
                double betJ = state.cells[i].bet;
                double rholpJ = state.cells[i].flui.MasEspLiq(state.cells[i].pres, state.cells[i].temp);
                double rholcJ = state.cells[i].fluicol.MasEspFlu(state.cells[i].pres, state.cells[i].temp);
                double rholmixJ = (1 - betJ) * rholpJ + betJ * rholcJ;

                double hidroM = sin(state.cells[i - 1].duto.teta) * (0.5 * state.cells[i - 1].dx) * (rholmix * (1 - alfE) + alfE * state.cells[i - 1].flui.MasEspGas(state.cells[i - 1].pres, state.cells[i - 1].temp)) * 9.82 / 98600.;
                double hidroJ = sin(state.cells[i].duto.teta) * (0.5 * state.cells[i].dx) * (rholmixJ * (1 - alfJ) + alfJ * state.cells[i].flui.MasEspGas(state.cells[i].pres, state.cells[i].temp)) * 9.82 / 98600.;

                double tit = alfE * state.cells[i - 1].flui.MasEspGas(state.cells[i - 1].pres, state.cells[i - 1].temp) / (state.cells[i - 1].flui.MasEspGas(state.cells[i - 1].pres, state.cells[i - 1].temp) * alfE + rholmix * (1. - alfE));

                double jtlM = (1. - betE) * state.cells[i - 1].flui.JTL(state.cells[i - 1].pres - hidroM, state.cells[i - 1].temp) - betE / rholcJ;
                double jtgM = state.cells[i - 1].flui.JTG(state.cells[i - 1].pres - hidroM, state.cells[i - 1].temp);
                tfonte = tE + ((1. - tit) * jtlM + tit * jtgM) * (state.cells[i].pres + hidroJ - state.cells[i - 1].pres - hidroM);

                hgF = state.cells[i - 1].flui.EntalpGas(state.cells[i - 1].pres, tfonte);
                hlF = state.cells[i - 1].flui.EntalpLiq(state.cells[i - 1].pres, tfonte);
                hcF = state.cells[i - 1].fluicol.CalorLiq(state.cells[i - 1].pres, tfonte) * tfonte;
                /*entalpia fluido complementar a ser corrigida*/

            } else {
                tfonte = state.cells[i].temp;
                hgF = state.cells[i - 1].flui.EntalpGas(state.cells[i - 1].pres, state.cells[i - 1].temp);
                hlF = state.cells[i - 1].flui.EntalpLiq(state.cells[i - 1].pres, state.cells[i - 1].temp);
                hcF = state.cells[i - 1].fluicol.CalorLiq(state.cells[i - 1].pres, state.cells[i - 1].temp) * tfonte;
                /*entalpia fluido complementar a ser corrigida*/
            }
        } else if ((*state.cells[i].acsrL).tipo == 8) {
            double alfM = state.cells[i - 1].alf;
            double betM = state.cells[i - 1].bet;

            double n = (*state.cells[i].acsrL).bvol.npoli;
            double ypres = (state.cells[i].pres) / (state.cells[i - 1].pres);
            tfonte = state.cells[i - 1].temp * pow(ypres, (n - 1) / n);

            hgF = state.cells[i - 1].flui.EntalpGas(state.cells[i - 1].pres, tfonte);
            hlF = state.cells[i - 1].flui.EntalpLiq(state.cells[i - 1].pres, tfonte);
            hcF = state.cells[i - 1].fluicol.CalorLiq(state.cells[i - 1].pres, tfonte) * tfonte;
            /*entalpia fluido complementar a ser corrigida*/
        } else {
            hgF = 0.;
            hlF = 0.;
        }
    } else {
        hgF = 0.;
        hlF = 0.;
    }

    fontemassL = 0;
    fontemassG = 0;
    if (state.cells[i].fontemassLR > 0.)
        fontemassL = hlF * state.cells[i].fontemassLR / dx;
    if (state.cells[i].fontemassCR > 0.)
        fontemassL += hcF * state.cells[i].fontemassCR / dx;
    if (fontemassG > 0.)
        fontemassG = hgF * state.cells[i].fontemassGR / dx;

    return energintmixT0 - (delFlux + hidro - (fontemassG + fontemassL) / area) * state.cells[i].dt;
}

double interpolateMixtureEnergy(const ThermalState &state, int i, int jp0, int jt, double razp) {
    int jp1 = jp0 + 1;

    double pres0 = state.cells[i].flui.rhogF[jp0][0];
    double pres1 = state.cells[i].flui.rhogF[jp1][0];
    double temp = state.cells[i].flui.rhogF[0][jt];
    double alfmed = state.cells[i].alf;
    double betmed = state.cells[i].bet;

    double rhogp0 = state.cells[i].flui.rhogF[jp0][jt];
    double rhogp1 = state.cells[i].flui.rhogF[jp1][jt];
    double hgp0 = state.cells[i].flui.HgF[jp0][jt];
    double hgp1 = state.cells[i].flui.HgF[jp1][jt];

    double rholp0 = state.cells[i].flui.rholF[jp0][jt];
    double rholp1 = state.cells[i].flui.rholF[jp1][jt];
    double hlp0 = state.cells[i].flui.HlF[jp0][jt];
    double hlp1 = state.cells[i].flui.HlF[jp1][jt];

    double rhocp0 = state.cells[i].fluicol.MasEspFlu(pres0, temp);
    double rhocp1 = state.cells[i].fluicol.MasEspFlu(pres1, temp);
    double hlc0 = state.cells[i].fluicol.CalorLiq(pres0, temp) * temp; // corrigir entalpia
    double hlc1 = state.cells[i].fluicol.CalorLiq(pres1, temp) * temp; // corrigir en;talpia

    double energ0 = alfmed * rhogp0 * (hgp0 - pres0 * 98066.5 / rhogp0) +
                    (1 - alfmed) * (1 - betmed) * rholp0 * (hlp0 - pres0 * 98066.5 / rholp0) +
                    (1 - alfmed) * (betmed)*rhocp0 * (hlc0 - pres0 * 98066.5 / rhocp0);

    double energ1 = alfmed * rhogp1 * (hgp1 - pres1 * 98066.5 / rhogp0) +
                    (1 - alfmed) * (1 - betmed) * rholp1 * (hlp1 - pres1 * 98066.5 / rholp1) +
                    (1 - alfmed) * (betmed)*rhocp1 * (hlc1 - pres1 * 98066.5 / rhocp1);

    return razp * energ0 + (1. - razp) * energ1;
}
}  // namespace sisprod::thermal
