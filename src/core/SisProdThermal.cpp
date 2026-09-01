#include "SisProdThermal.h"

#include "Leitura.h"
#include "celula3.h"
#include "solver3DPoisson.h"

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

void updateTemperatureFromEnthalpy(const ThermalState &state, int i) {
    double pres = state.cells[i].pres;
    double **Var = state.cells[i].flui.rholF;

    int ipres = 0;
    int ipmarcador;
    int ndiv = state.cells[i].flui.npontos - 1;
    if (pres < Var[1][0] || pres >= Var[ndiv + 1][0]) {
        cout << "pressÃƒÂ£o fora dos limites de tabela";
        getchar();
    }

    int e, m, d;
    e = 1;
    d = ndiv + 1;
    while (e <= d) {
        m = (e + d) / 2;
        ipmarcador = m;
        if (m == 1) {
            ipres = m;
            break;
        } else if (m == ndiv + 1 && Var[m][0] == pres) {
            ipres = m - 1;
            break;
        }
        if (Var[m][0] > pres && Var[m - 1][0] <= pres) {
            ipres = m - 1;
            break;
        }
        if (Var[m][0] < pres)
            e = m + 1;
        else
            d = m - 1;
    }

    double razpres = 1. - (Var[ipres][0] - pres) / (Var[ipres][0] - Var[ipres + 1][0]);

    double energint = computeMixtureEnthalpy(state, i);
    int j = 1;
    double val1 = -1;
    double val2 = 0;
    while (j < ndiv + 1 || (energint >= val1 && energint <= val2) || (energint <= val1 && energint >= val2)) {
        val1 = interpolateMixtureEnergy(state, i, ipres, j, razpres);
        val2 = interpolateMixtureEnergy(state, i, ipres, j + 1, razpres);
        j++;
    }
    double raztemp = 1. - (val1 - energint) / (val1 - val2);
    state.cells[i].temp = Var[0][j] * raztemp + (1. - raztemp) * Var[0][j + 1];

    if (state.cells[i].temp < -50.)
        state.cells[i].temp = -50.;
    if (state.cells[i].temp > 200.)
        state.cells[i].temp = 200.;
}

namespace {

struct TemperatureBalance {
    double cellLength;
    double meanCellLength;
    double flowArea;
    double voidFraction;
    double bet;
    double gasSuperficialVelocity;
    double liquidSuperficialVelocity;
    double referenceMixtureVelocity;
    double rp;
    double rc;
    double liquidDensity;
    double gasDensity;
    double liquidHeatCapacity;
    double liquidIsochoricHeatCapacity;
    double gasHeatCapacity;
    double gasIsochoricHeatCapacity;
    double liquidJouleThomson;
    double gasJouleThomson;
    double hydrostaticPower;
    double heatFlux;
    double timeCoefficient;
    double pressureTimeCoefficient;
    double temperatureSpatialCoefficient;
    double pressureSpatialCoefficient;
};

TemperatureBalance prepareTemperatureBalance(const ThermalState &state,
                                             int cellIndex,
                                             int steadyStateMode) {
    double dx = state.cells[cellIndex].dx;
    double dxmed;
    if (cellIndex > 0)
        dxmed = 0.5 * (state.cells[cellIndex].dx + state.cells[cellIndex - 1].dx);
    else
        dxmed = 0.5 * state.cells[cellIndex].dx;
    double dia = state.cells[cellIndex].duto.a;
    double area = 0.25 * M_PI * dia * dia;
    double alfmed = state.cells[cellIndex].alf;
    double betmed = state.cells[cellIndex].bet;
    double ugsmed;
    double ulsmed;
    if (cellIndex > 0 && (state.cells[cellIndex - 1].acsr.tipo != 5 ||
                          state.cells[cellIndex - 1].acsr.chk.AreaGarg > (1e-3 + state.input.master1.razareaativ) * state.cells[cellIndex - 1].duto.area)) {
        if (state.cells[cellIndex].alf > (*state.globals).localtiny)
            ugsmed = state.cells[cellIndex].QG / area;
        else {
            ugsmed = 0.;
        }
        if (state.cells[cellIndex].alf < 1. - (*state.globals).localtiny)
            ulsmed = state.cells[cellIndex].QL / area;
        else {
            ulsmed = 0.;
        }
    } else {
        if (state.cells[cellIndex].alf > (*state.globals).localtiny)
            ugsmed = state.cells[cellIndex + 1].QG / area;
        else {
            ugsmed = 0.;
        }
        if (state.cells[cellIndex].alf < 1. - (*state.globals).localtiny)
            ulsmed = state.cells[cellIndex + 1].QL / area;
        else {
            ulsmed = 0.;
        }
    }
    double JmixRef = fabs(ugsmed) + fabs(ulsmed);
    double rp = state.cells[cellIndex].rpC;
    double rc = state.cells[cellIndex].rcC;
    double rhol = (1. - betmed) * rp + betmed * rc;
    double rhog = state.cells[cellIndex].rgC;
    double cpl = (1. - betmed) * state.cells[cellIndex].flui.CalorLiq(state.cells[cellIndex].presini, state.cells[cellIndex].temp) + betmed * state.cells[cellIndex].fluicol.CalorLiq(state.cells[cellIndex].presini, state.cells[cellIndex].temp);
    double cvl = cpl;
    double cpg = state.cells[cellIndex].flui.CalorGas(state.cells[cellIndex].presini, state.cells[cellIndex].temp);
    double cvg = state.cells[cellIndex].flui.CalorGasVolMod(state.cells[cellIndex].presini, state.cells[cellIndex].temp, state.cells[cellIndex].rgC);
    double jtl = (1. - betmed) * state.cells[cellIndex].flui.JTL(state.cells[cellIndex].presini, state.cells[cellIndex].temp) - betmed / rc;
    double jtg = state.cells[cellIndex].flui.JTG(state.cells[cellIndex].presini, state.cells[cellIndex].temp, state.cells[cellIndex].rgC);
    double hidro = (rhol * ulsmed + rhog * ugsmed) * area * 9.82 * sin(state.cells[cellIndex].duto.teta);

    double fluxcal = 0.;
    state.cells[cellIndex].calor.Tint = state.cells[cellIndex].temp;
    if (cellIndex > 0)
        state.cells[cellIndex].calor.dtL = state.cells[cellIndex].temp - state.cells[cellIndex - 1].tempini;
    else
        state.cells[cellIndex].calor.dtL = 0;
    state.cells[cellIndex].calor.Vint = ugsmed + ulsmed;
    state.cells[cellIndex].calor.dt = state.cells[cellIndex].dt;
    double condliq = (1. - betmed) * state.cells[cellIndex].flui.CondLiq(state.cells[cellIndex].presini, state.cells[cellIndex].temp) + betmed * state.cells[cellIndex].fluicol.CondLiq(state.cells[cellIndex].presini, state.cells[cellIndex].temp);
    state.cells[cellIndex].calor.kint = condliq * (1 - alfmed) + state.cells[cellIndex].flui.CondGas(state.cells[cellIndex].presini, state.cells[cellIndex].temp) * alfmed;
    state.cells[cellIndex].calor.cpint = cpl * (1 - alfmed) + cpg * alfmed;
    state.cells[cellIndex].calor.rhoint = rhol * (1 - alfmed) + rhog * alfmed;
    double viscliq = (1. - betmed) * state.cells[cellIndex].mipC + betmed * state.cells[cellIndex].micC;
    state.cells[cellIndex].calor.viscint = viscliq * (1 - alfmed) * 1.e-3 + state.cells[cellIndex].migC * alfmed * 1.e-3;
    double dtemp = state.cells[cellIndex].temp * 0.01;
    if (fabs(state.cells[cellIndex].temp) < 1e-15)
        dtemp = 0.1;
    double rholdT = (1. - betmed) * state.cells[cellIndex].flui.MasEspLiq(state.cells[cellIndex].presini, state.cells[cellIndex].temp + dtemp) +
                    betmed * state.cells[cellIndex].fluicol.MasEspFlu(state.cells[cellIndex].presini, state.cells[cellIndex].temp + dtemp) - rhol;
    double rhogdT = state.cells[cellIndex].flui.MasEspGas(state.cells[cellIndex].presini, state.cells[cellIndex].temp + dtemp) - rhog;
    state.cells[cellIndex].calor.betint = -(1 / state.cells[cellIndex].calor.rhoint) * (rholdT * (1 - alfmed) + rhogdT * alfmed) / (dtemp);
    if (state.input.modoDifus3D == 0 || steadyStateMode != 0) {
        if (steadyStateMode == 0)
            fluxcal = state.cells[cellIndex].calor.transtrans();
        else
            fluxcal = state.cells[cellIndex].calor.transperm();
        if (state.productionNetworkCoupled == 1 && cellIndex >= state.primarySectionStart && cellIndex <= state.primarySectionEnd) {
            fluxcal -= state.cells[cellIndex].fluxcalAcopRedeP;
        }
    } else if (state.input.modoDifus3D == 1) {
        int acoplado = -1;
        int icelAcop;
        for (int iacop = 0; iacop < state.input.nacop; iacop++) {
            icelAcop = state.input.celAcop[iacop].indCel;
            if (cellIndex == icelAcop) {
                acoplado = iacop;
                break;
            }
        }
        if (acoplado == -1)
            fluxcal = state.cells[cellIndex].calor.transtrans();
        else {
            int iacop1 = state.coupledCellIndices[acoplado];
            fluxcal = -state.input.celAcop[acoplado].FE * state.poissonSolver.dados.qTotal[iacop1] / state.cells[cellIndex].dx;
        }
    }

    state.cells[cellIndex].fluxcalmed = fluxcal;

    double coefTempo = (rhol * (1 - alfmed) * cvl + rhog * alfmed * cvg) * area;
    double coefPresTempo;
    coefPresTempo = -state.cells[cellIndex].flui.CalorGasPresMod(state.cells[cellIndex].presini, state.cells[cellIndex].temp) * (rhog * alfmed * area); //-state.cells[cellIndex].flui.CalorGasPresMod(state.cells[cellIndex].pres, state.cells[cellIndex].temp) * (rhog * alfmed * area);
    double coefdxT = (rhol * ulsmed * cpl + rhog * ugsmed * cpg) * area;
    double coefdxP = (rhol * ulsmed * jtl + rhog * ugsmed * jtg) * area;

    return TemperatureBalance{
        .cellLength = dx,
        .meanCellLength = dxmed,
        .flowArea = area,
        .voidFraction = alfmed,
        .bet = betmed,
        .gasSuperficialVelocity = ugsmed,
        .liquidSuperficialVelocity = ulsmed,
        .referenceMixtureVelocity = JmixRef,
        .rp = rp,
        .rc = rc,
        .liquidDensity = rhol,
        .gasDensity = rhog,
        .liquidHeatCapacity = cpl,
        .liquidIsochoricHeatCapacity = cvl,
        .gasHeatCapacity = cpg,
        .gasIsochoricHeatCapacity = cvg,
        .liquidJouleThomson = jtl,
        .gasJouleThomson = jtg,
        .hydrostaticPower = hidro,
        .heatFlux = fluxcal,
        .timeCoefficient = coefTempo,
        .pressureTimeCoefficient = coefPresTempo,
        .temperatureSpatialCoefficient = coefdxT,
        .pressureSpatialCoefficient = coefdxP,
    };
}

double computeKineticTemperatureTerm(const ThermalState &state, int cellIndex,
                                     TemperatureBalance &balance) {
    double cinetico = 0;
    double ugmed0 = 0;
    double ulmed0 = 0;
    double ugmed = 0;
    double ulmed = 0;
    double ugmedini = 0;
    double ulmedini = 0;
    if (cellIndex <= state.lastCell - 1 && state.cells[cellIndex].acsr.tipo == 0 && state.cells[cellIndex + 1].acsr.tipo == 0) {

        double dxCin = state.cells[cellIndex].dx;
        double diaaux = state.cells[cellIndex].duto.a;
        double areaaux = 0.25 * M_PI * diaaux * diaaux;

        double rpcin = state.cells[cellIndex + 1].rpCi;
        double rccin = state.cells[cellIndex + 1].rcCi;
        double rholcin = (1. - balance.bet) * rpcin + balance.bet * rccin;
        double rhogcin = state.cells[cellIndex + 1].rgCi;

        double alf;
        double alfini;

        balance.gasSuperficialVelocity = state.cells[cellIndex + 1].QG / balance.flowArea;
        balance.liquidSuperficialVelocity = state.cells[cellIndex + 1].QL / balance.flowArea;

        if (state.cells[cellIndex + 1].QG > 0)
            alf = state.cells[cellIndex].alf;
        else
            alf = state.cells[cellIndex + 1].alf;

        if (state.cells[cellIndex + 1].QGini > 0)
            alfini = state.cells[cellIndex].alfini;
        else
            alfini = state.cells[cellIndex + 1].alfini;

        double alf0;
        if (state.cells[cellIndex].QG > 0)
            alf0 = state.cells[cellIndex - 1].alf;
        else
            alf0 = state.cells[cellIndex].alf;

        if (alf0 > 1e-3) {
            ugmed0 = state.cells[cellIndex].QG / (areaaux);
            ugmed0 /= alf0;
        }

        if (alf0 < 1. - 1e-3) {
            ulmed0 = state.cells[cellIndex].QL / (areaaux);
            ulmed0 /= (1. - alf0);
        }

        if (alf > 1e-3) {
            ugmed = balance.gasSuperficialVelocity;
            ugmed /= alf;
        }

        if (alf < 1. - 1e-3) {
            ulmed = balance.liquidSuperficialVelocity;
            ulmed /= (1. - alf);
        }

        if (alfini > 1e-3) {

            ugmedini = state.cells[cellIndex + 1].QGini / (areaaux);
            ugmedini /= alfini;
        }

        if (alfini < 1. - 1e-3) {

            ulmedini = state.cells[cellIndex + 1].QLini / (areaaux);
            ulmedini /= (1. - alfini);
        }

        cinetico = rholcin * (1 - alf) * areaaux * (0.5 * (ulmed * ulmed - ulmedini * ulmedini)) / state.cells[cellIndex].dt +
                   rhogcin * alf * areaaux * (0.5 * (ugmed * ugmed - ugmedini * ugmedini)) / state.cells[cellIndex].dt +
                   ((state.cells[cellIndex + 1].MC - state.cells[cellIndex + 1].Mliqini) * ugmed * (ugmed - ugmed0) / dxCin +
                    state.cells[cellIndex + 1].Mliqini * ulmed * (ulmed - ulmed0) / dxCin);
    }
    return cinetico;
}

struct TemperatureSourceTerms {
    double gas;
    double liquid;
};

TemperatureSourceTerms computeTemperatureSourceTerms(const ThermalState &state,
                                                      int cellIndex,
                                                      double cellLength) {
    double fontemassG = 0.;
    double fontemassL = 0.;
    double fontemassC = 0.;
    double tfonte = state.cells[cellIndex].temp;
    double cpgF;
    double razcpF = 0.;
    double cplF;
    if (state.cells[cellIndex].acsr.tipo == 1) {
        tfonte = state.cells[cellIndex].acsr.injg.temp;
        cpgF = state.cells[cellIndex].acsr.injg.FluidoPro.CalorGas(state.cells[cellIndex].presini, state.cells[cellIndex].temp);
        razcpF = state.cells[cellIndex].acsr.injg.FluidoPro.ConstAdG(state.cells[cellIndex].presini, state.cells[cellIndex].temp);
        cplF = 0.;
    } else if (state.cells[cellIndex].acsr.tipo == 2) {
        tfonte = state.cells[cellIndex].acsr.injl.temp;
        cpgF = state.cells[cellIndex].acsr.injl.FluidoPro.CalorGas(state.cells[cellIndex].presini, state.cells[cellIndex].temp);                                                                                                                                        // state.cells[cellIndex].acsr.injl.FluidoPro.CalorGas(state.cells[cellIndex].pres, state.cells[cellIndex].temp);
        razcpF = state.cells[cellIndex].acsr.injl.FluidoPro.CalorGas(state.cells[cellIndex].presini, state.cells[cellIndex].temp);                                                                                                                                      // state.cells[cellIndex].acsr.injl.FluidoPro.ConstAdG(state.cells[cellIndex].pres, state.cells[cellIndex].temp);
        cplF = (1. - state.cells[cellIndex].acsr.injl.bet) * state.cells[cellIndex].acsr.injl.FluidoPro.CalorLiq(state.cells[cellIndex].presini, state.cells[cellIndex].temp) + state.cells[cellIndex].acsr.injl.bet * state.cells[cellIndex].acsr.injl.fluidocol.CalorLiq(state.cells[cellIndex].presini, state.cells[cellIndex].temp); //(1. - state.cells[cellIndex].acsr.injl.bet) * state.cells[cellIndex].acsr.injl.FluidoPro.CalorLiq(state.cells[cellIndex].pres, state.cells[cellIndex].temp)
    } else if (state.cells[cellIndex].acsr.tipo == 10) {
        tfonte = state.cells[cellIndex].acsr.injm.temp;
        cpgF = state.cells[cellIndex].acsr.injm.FluidoPro.CalorGas(state.cells[cellIndex].presini, state.cells[cellIndex].temp);
        razcpF = state.cells[cellIndex].acsr.injm.FluidoPro.ConstAdG(state.cells[cellIndex].presini, state.cells[cellIndex].temp);
        if ((state.cells[cellIndex].fontemassLR + state.cells[cellIndex].fontemassCR) > 0.) {
            double titbet = state.cells[cellIndex].fontemassCR / (state.cells[cellIndex].fontemassLR + state.cells[cellIndex].fontemassCR);
            cplF = (1. - titbet) * state.cells[cellIndex].acsr.injl.FluidoPro.CalorLiq(state.cells[cellIndex].presini, state.cells[cellIndex].temp) +
                   titbet * state.cells[cellIndex].acsr.injl.fluidocol.CalorLiq(state.cells[cellIndex].presini, state.cells[cellIndex].temp);
        } else
            cplF = 0.;
    } else if (state.cells[cellIndex].acsr.tipo == 3) {
        tfonte = state.cells[cellIndex].acsr.ipr.Tres;
        cpgF = state.cells[cellIndex].acsr.ipr.FluidoPro.CalorGas(state.cells[cellIndex].presini, state.cells[cellIndex].temp);
        razcpF = state.cells[cellIndex].acsr.ipr.FluidoPro.ConstAdG(state.cells[cellIndex].presini, state.cells[cellIndex].temp);
        cplF = state.cells[cellIndex].acsr.ipr.FluidoPro.CalorLiq(state.cells[cellIndex].presini, state.cells[cellIndex].temp);
    } else if (state.cells[cellIndex].acsr.tipo == 15) {
        tfonte = state.cells[cellIndex].acsr.radialPoro.tRes;
        cpgF = state.cells[cellIndex].acsr.radialPoro.flup.CalorGas(state.cells[cellIndex].presini, state.cells[cellIndex].temp);
        razcpF = state.cells[cellIndex].acsr.radialPoro.flup.ConstAdG(state.cells[cellIndex].presini, state.cells[cellIndex].temp);
        cplF = state.cells[cellIndex].acsr.radialPoro.flup.CalorLiq(state.cells[cellIndex].presini, state.cells[cellIndex].temp);
    } else if (state.cells[cellIndex].acsr.tipo == 9) {
        double betM = state.cells[cellIndex].acsr.fontechk.betISamb;
        double pamb = state.cells[cellIndex].acsr.fontechk.pamb;
        double tamb = state.cells[cellIndex].acsr.fontechk.tamb;
        cpgF = state.cells[cellIndex].acsr.fontechk.fluidoPamb.CalorGas(pamb, tamb);
        cplF = (1. - betM) * state.cells[cellIndex].acsr.fontechk.fluidoPamb.CalorLiq(pamb, tamb) + betM * state.cells[cellIndex].acsr.fontechk.fluidocol.CalorLiq(pamb, tamb);
        razcpF = state.cells[cellIndex].acsr.fontechk.fluidoPamb.ConstAdG(pamb, tamb);
        tfonte = state.cells[cellIndex].acsr.fontechk.tamb;
    } else if (state.cells[cellIndex].acsrL != 0 && cellIndex < state.lastCell) {
        if ((*state.cells[cellIndex].acsrL).tipo == 5 && (state.cells[cellIndex].fontemassLR + state.cells[cellIndex].fontemassCR + state.cells[cellIndex].fontemassGR) > 0.) {
            if ((*state.cells[cellIndex].acsrL).chk.AreaGarg < state.input.master1.razareaativ * state.cells[cellIndex].dutoL.area && (*state.cells[cellIndex].acsrL).chk.AreaGarg > 1e-5 * state.cells[cellIndex].dutoL.area) {
                double tE = state.cells[cellIndex - 1].tempini;
                double alfE = state.cells[cellIndex - 1].alf;
                double betE = state.cells[cellIndex - 1].bet;

                double rholp = state.cells[cellIndex - 1].rpC;
                double rholc = state.cells[cellIndex - 1].rcC;
                double rholmix = (1 - betE) * rholp + betE * rholc;

                double alfJ = state.cells[cellIndex].alf;
                double betJ = state.cells[cellIndex].bet;
                double rholpJ = state.cells[cellIndex].rpC;
                double rholcJ = state.cells[cellIndex].rcC;
                double rholmixJ = (1 - betJ) * rholpJ + betJ * rholcJ;

                double hidroM = sin(state.cells[cellIndex - 1].duto.teta) * (0.5 * state.cells[cellIndex - 1].dx) * (rholmix * (1 - alfE) + alfE * state.cells[cellIndex - 1].rgC) * 9.82 / 98600.;
                double hidroJ = sin(state.cells[cellIndex].duto.teta) * (0.5 * state.cells[cellIndex].dx) * (rholmixJ * (1 - alfJ) + alfJ * state.cells[cellIndex].rgC) * 9.82 / 98600.;
                double tit = alfE * state.cells[cellIndex - 1].rgC / (state.cells[cellIndex - 1].rgC * alfE + rholmix * (1. - alfE));
                cpgF = state.cells[cellIndex - 1].flui.CalorGas(state.cells[cellIndex - 1].presini, state.cells[cellIndex - 1].tempini);                                                                                                                              // state.cells[cellIndex - 1].flui.CalorGas(state.cells[cellIndex - 1].pres, state.cells[cellIndex - 1].temp);
                cplF = (1 - state.cells[cellIndex - 1].bet) * state.cells[cellIndex - 1].flui.CalorLiq(state.cells[cellIndex - 1].presini, state.cells[cellIndex - 1].tempini) + state.cells[cellIndex - 1].bet * state.cells[cellIndex - 1].fluicol.CalorLiq(state.cells[cellIndex - 1].presini, state.cells[cellIndex - 1].tempini); //(1 - state.cells[cellIndex - 1].bet) * state.cells[cellIndex - 1].flui.CalorLiq(state.cells[cellIndex - 1].pres, state.cells[cellIndex - 1].temp)
                double jtlM = (1. - betE) * state.cells[cellIndex - 1].flui.JTL(state.cells[cellIndex - 1].presini - hidroM, state.cells[cellIndex - 1].tempini) - betE / rholcJ;                                                                                     //(1. - betE) * state.cells[cellIndex - 1].flui.JTL(state.cells[cellIndex - 1].pres - hidroM, state.cells[cellIndex - 1].temp)
                double jtgM = state.cells[cellIndex - 1].flui.JTG(state.cells[cellIndex - 1].presini - hidroM, state.cells[cellIndex - 1].tempini);                                                                                                                   // state.cells[cellIndex - 1].flui.JTG(state.cells[cellIndex - 1].pres - hidroM, state.cells[cellIndex - 1].temp);
                tfonte = tE + ((1. - tit) * jtlM / cplF + tit * jtgM / cpgF) *
                                  (state.cells[cellIndex].pres + hidroJ - state.cells[cellIndex - 1].pres - hidroM) * 98066.52;

                razcpF = state.cells[cellIndex - 1].flui.ConstAdG(state.cells[cellIndex - 1].presini, state.cells[cellIndex - 1].tempini);

            } else {
                cpgF = state.cells[cellIndex - 1].flui.CalorGas(state.cells[cellIndex - 1].presini, state.cells[cellIndex - 1].tempini);
                cplF = (1 - state.cells[cellIndex - 1].bet) * state.cells[cellIndex - 1].flui.CalorLiq(state.cells[cellIndex - 1].presini, state.cells[cellIndex - 1].tempini) + state.cells[cellIndex - 1].bet * state.cells[cellIndex - 1].fluicol.CalorLiq(state.cells[cellIndex - 1].presini, state.cells[cellIndex - 1].tempini);
                tfonte = state.cells[cellIndex].temp;
                razcpF = state.cells[cellIndex - 1].flui.ConstAdG(state.cells[cellIndex - 1].presini, state.cells[cellIndex - 1].tempini);
            }
        } else if ((*state.cells[cellIndex].acsrL).tipo == 8) {
            double betM = state.cells[cellIndex - 1].bet;
            cpgF = state.cells[cellIndex - 1].flui.CalorGas(state.cells[cellIndex - 1].presini, state.cells[cellIndex - 1].tempini);
            cplF = (1. - betM) * state.cells[cellIndex - 1].flui.CalorLiq(state.cells[cellIndex - 1].presini, state.cells[cellIndex - 1].tempini) + betM * state.cells[cellIndex - 1].fluicol.CalorLiq(state.cells[cellIndex - 1].presini, state.cells[cellIndex - 1].tempini);
            razcpF = state.cells[cellIndex - 1].flui.ConstAdG(state.cells[cellIndex - 1].presini, state.cells[cellIndex - 1].tempini);
            double n = (*state.cells[cellIndex].acsrL).bvol.npoli;
            double ypres = (state.cells[cellIndex].pres) / (state.cells[cellIndex - 1].pres);
            tfonte = state.cells[cellIndex - 1].tempini * pow(ypres, (n - 1) / n);
        } else {
            cpgF = 0.;
            razcpF = 1.;
            cplF = 0.;
        }
    } else if (cellIndex == state.lastCell && (state.cells[cellIndex].fontemassLR + state.cells[cellIndex].fontemassCR + state.cells[cellIndex].fontemassGR) > 0.) {
        if ((*state.globals).chaverede == 0 || state.networkEndpoint == 1 || (*state.globals).chaveRedeParalela == 1)
            tfonte = state.cells[cellIndex].calor.Textern1;
        else
            tfonte = state.gasSurfaceTemperature;
        cpgF = state.cells[cellIndex].flui.CalorGas(state.cells[cellIndex].presini, tfonte);
        double betM = 0.;
        if (fabs(state.cells[cellIndex].fontemassCR) > 1e-15)
            betM = state.cells[cellIndex].fontemassCR /
                   (state.cells[cellIndex].fontemassCR + state.cells[cellIndex].fontemassLR);
        cplF = (1. - betM) * state.cells[cellIndex].flui.CalorLiq(state.cells[cellIndex].presini, tfonte) + betM * state.cells[cellIndex].fluicol.CalorLiq(state.cells[cellIndex].presini, tfonte);
        razcpF = state.cells[cellIndex].flui.ConstAdG(state.cells[cellIndex].pres, tfonte);
    } else {
        cpgF = 0.;
        razcpF = 1.;
        cplF = 0.;
    }

    fontemassL = 0;
    if (state.cells[cellIndex].fontemassLR > 0.)
        fontemassL = state.cells[cellIndex].fontemassLR / cellLength;
    if (state.cells[cellIndex].fontemassCR > 0.)
        fontemassL += state.cells[cellIndex].fontemassCR / cellLength;
    fontemassL *= (cplF) * (tfonte - state.cells[cellIndex].temp);

    fontemassG = state.cells[cellIndex].fontemassGR / cellLength;
    if (fontemassG > 0.)
        fontemassG *= (cpgF / razcpF) * (tfonte - state.cells[cellIndex].temp);
    else
        fontemassG = 0;

    return TemperatureSourceTerms{
        .gas = fontemassG,
        .liquid = fontemassL,
    };
}

TemperatureSourceTerms computeThermalMassTransferSourceTerms(
    const ThermalState &state, int i, double dx) {
    double fontemassG = 0.;
    double fontemassL = 0.;
    double tfonte = state.cells[i].temp;
    double cpgF;
    double razcpF = 0.;
    double cplF;
    if (state.cells[i].acsr.tipo == 1) {
        tfonte = state.cells[i].acsr.injg.temp;
        cpgF = state.cells[i].acsr.injg.FluidoPro.CalorGas(state.cells[i].pres, state.cells[i].temp);
        razcpF = state.cells[i].acsr.injg.FluidoPro.ConstAdG(state.cells[i].pres, state.cells[i].temp);
        cplF = 0.;
    } else if (state.cells[i].acsr.tipo == 2) {
        tfonte = state.cells[i].acsr.injl.temp;
        ;
        cpgF = state.cells[i].acsr.injl.FluidoPro.CalorGas(state.cells[i].pres, state.cells[i].temp);
        razcpF = state.cells[i].acsr.injl.FluidoPro.ConstAdG(state.cells[i].pres, state.cells[i].temp);
        cplF = (1. - state.cells[i].acsr.injl.bet) * state.cells[i].acsr.injl.FluidoPro.CalorLiq(state.cells[i].pres, state.cells[i].temp) + state.cells[i].acsr.injl.bet * state.cells[i].acsr.injl.fluidocol.CalorLiq(state.cells[i].pres, state.cells[i].temp);
    } else if (state.cells[i].acsr.tipo == 3) {
        tfonte = state.cells[i].acsr.ipr.Tres;
        cpgF = state.cells[i].acsr.ipr.FluidoPro.CalorGas(state.cells[i].pres, state.cells[i].temp);
        razcpF = state.cells[i].acsr.ipr.FluidoPro.ConstAdG(state.cells[i].pres, state.cells[i].temp);
        cplF = state.cells[i].acsr.ipr.FluidoPro.CalorLiq(state.cells[i].pres, state.cells[i].temp);
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

                cpgF = state.cells[i - 1].flui.CalorGas(state.cells[i - 1].pres, state.cells[i - 1].temp);
                cplF = (1 - state.cells[i - 1].bet) * state.cells[i - 1].flui.CalorLiq(state.cells[i - 1].pres, state.cells[i - 1].temp) + state.cells[i - 1].bet * state.cells[i - 1].fluicol.CalorLiq(state.cells[i - 1].pres, state.cells[i - 1].temp);
                razcpF = state.cells[i - 1].flui.ConstAdG(state.cells[i - 1].pres, state.cells[i - 1].temp);

            } else {
                cpgF = state.cells[i - 1].flui.CalorGas(state.cells[i - 1].pres, state.cells[i - 1].temp);
                cplF = (1 - state.cells[i - 1].bet) * state.cells[i - 1].flui.CalorLiq(state.cells[i - 1].pres, state.cells[i - 1].temp) + state.cells[i - 1].bet * state.cells[i - 1].fluicol.CalorLiq(state.cells[i - 1].pres, state.cells[i - 1].temp);
                tfonte = state.cells[i].temp;
                razcpF = state.cells[i - 1].flui.ConstAdG(state.cells[i - 1].pres, state.cells[i - 1].temp);
            }
        } else if ((*state.cells[i].acsrL).tipo == 8) {
            double alfM = state.cells[i - 1].alf;
            double betM = state.cells[i - 1].bet;
            cpgF = state.cells[i - 1].flui.CalorGas(state.cells[i - 1].pres, state.cells[i - 1].temp);
            cplF = (1. - betM) * state.cells[i - 1].flui.CalorLiq(state.cells[i - 1].pres, state.cells[i - 1].temp) + betM * state.cells[i - 1].fluicol.CalorLiq(state.cells[i - 1].pres, state.cells[i - 1].temp);
            double n = (*state.cells[i].acsrL).bvol.npoli;
            double ypres = (state.cells[i].pres) / (state.cells[i - 1].pres);
            tfonte = state.cells[i - 1].temp * pow(ypres, (n - 1) / n);
        } else {
            cpgF = 0.;
            razcpF = 1.;
            cplF = 0.;
        }
    } else {
        cpgF = 0.;
        razcpF = 1.;
        cplF = 0.;
    }

    fontemassL = 0;
    if (state.cells[i].fontemassLR > 0.)
        fontemassL = state.cells[i].fontemassLR / dx;
    if (state.cells[i].fontemassCR > 0.)
        fontemassL += state.cells[i].fontemassCR / dx;
    fontemassL *= (cplF) * (tfonte - state.cells[i].temp);

    fontemassG = state.cells[i].fontemassGR / dx;
    if (fontemassG > 0.)
        fontemassG *= (cpgF / razcpF) * (tfonte - state.cells[i].temp);
    else
        fontemassG = 0;

    return TemperatureSourceTerms{
        .gas = fontemassG,
        .liquid = fontemassL,
    };
}

}  // namespace

void computeTemperature(const ThermalState &state, int i, double tempantiga, int modoPerm) {
    if (state.thermalSourceDisabled == 0) {
        TemperatureBalance balance = prepareTemperatureBalance(state, i, modoPerm);
        double dpdx;
        if (i < state.lastCell - 1)
            dpdx = 2. * (state.cells[i + 1].presaux - state.cells[i].pres) * 98066.5 / state.cells[i].dx;
        else if (i == state.lastCell - 1 && state.surfaceChoke.AreaGarg > 0.5 * state.cells[state.lastCell - 1].duto.area)
            dpdx = 2. * (state.cells[i + 1].presaux - state.cells[i].pres) * 98066.5 / state.cells[i].dx;
        else
            dpdx = 2. * (state.cells[i].pres - state.cells[i].presaux) * 98066.5 / state.cells[i].dx;
        if (i > 0 && state.cells[i - 1].acsr.tipo == 5 &&
            state.cells[i - 1].acsr.chk.AreaGarg <= (1e-3 + state.input.master1.razareaativ) * state.cells[i].duto.area)
            dpdx = 2. * (state.cells[i + 1].presaux - state.cells[i].pres) * 98066.5 / state.cells[i].dx;
        else if (i > 0 && state.cells[i].acsr.tipo == 5 &&
                 state.cells[i].acsr.chk.AreaGarg <= (1e-3 + state.input.master1.razareaativ) * state.cells[i].duto.area) {
            dpdx = (state.cells[i].pres - state.cells[i - 1].pres) * 98066.5 / balance.meanCellLength;
        } else if (i == 0)
            dpdx = 2. * (state.cells[i + 1].presaux - state.cells[i].pres) * 98066.5 / state.cells[i].dx;

        state.cells[i].VTemper = balance.temperatureSpatialCoefficient / balance.timeCoefficient;
        if ((i == 0 && state.cells[i].VTemper < 0.) || (((i < state.lastCell || state.cells[i].VTemper >= 0.) && i > 0) ||
                                                   (i == state.lastCell && state.surfaceChokeMassCondition == 1) || (i == state.lastCell && state.input.chkv == 1))) {
            double dtdx = 0.;
            if (i > 0)
                dtdx = (state.cells[i].temp - state.cells[i - 1].tempini) / balance.meanCellLength;

            if (i < state.lastCell)
                if (state.cells[i].VTemper < 0)
                    dtdx = (state.cells[i + 1].tempini - state.cells[i].temp) / (0.5 * (state.cells[i + 1].dx + state.cells[i].dx));
            if (state.cells[i].acsr.tipo == 5 && state.cells[i].acsr.chk.AreaGarg <= (1e-3 + state.input.master1.razareaativ) * state.cells[i].duto.area && state.cells[i].VTemper <= 0)
                dtdx = 0.;
            if (i > 0 && state.cells[i - 1].acsr.tipo == 5 &&
                state.cells[i - 1].acsr.chk.AreaGarg <= (1e-3 + state.input.master1.razareaativ) * state.cells[i - 1].duto.area && state.cells[i].VTemper >= 0)
                dtdx = 0 * (state.cells[i + 1].tempini - state.cells[i].temp) / (0.5 * (state.cells[i + 1].dx + state.cells[i].dx));
            if ((i <= 1 && state.cells[i].VTemper <= 0) || (i == state.lastCell && state.cells[i].VTemper <= 0))
                dtdx = 0.;
            if (state.cells[i].acsr.tipo == 8 && state.cells[i].acsr.bvol.freq > 1.) {
                dpdx = (state.cells[i].pres - state.cells[i - 1].pres) * 98066.5 / balance.meanCellLength;
            }

            double cinetico = computeKineticTemperatureTerm(state, i, balance);

            TemperatureSourceTerms sourceTerms =
                computeTemperatureSourceTerms(state, i, balance.cellLength);

            double latente;
            double valTransMass = fabs(state.cells[i].FonteMudaFase);
            double sigTransMass = 1.;
            if (valTransMass > 1e-25)
                sigTransMass = state.cells[i].FonteMudaFase / valTransMass;
            if (state.input.limTransMass < valTransMass)
                valTransMass = sigTransMass * state.input.limTransMass;
            else
                valTransMass *= sigTransMass;
            if (state.latentHeatEnabled > 0 && state.input.flashCompleto == 0) {
                latente = interpolateLatentHeat(state, state.cells[i].presini, state.cells[i].temp) * valTransMass;
            } else if ((state.input.flashCompleto == 1 || state.input.flashCompleto == 2) && state.latentHeatEnabled > 0)
                latente = (state.cells[i].flui.EntalpGas(state.cells[i].presini, state.cells[i].temp) -
                           state.cells[i].flui.EntalpLiq(state.cells[i].presini, state.cells[i].temp)) * valTransMass;
            else
                latente = 0;

            if (state.input.latente == 0)
                latente = 0.;
            else if (state.input.condlatente == 0 && latente < 0)
                latente = 0.;

            double alfinter;
            double alfinterL;
            if (balance.gasSuperficialVelocity > 0) {
                alfinter = state.cells[i].alf;
                alfinterL = state.cells[i].alfL;
            } else {
                alfinter = state.cells[i].alfR;
                alfinterL = state.cells[i].alf;
            }
            double delvel;
            if (alfinter > (*state.globals).localtiny && alfinter < (1. - (*state.globals).localtiny))
                delvel = balance.gasSuperficialVelocity / alfinter - balance.liquidSuperficialVelocity / (1. - alfinter);
            else if (alfinter > (*state.globals).localtiny)
                delvel = balance.gasSuperficialVelocity;
            else
                delvel = balance.liquidSuperficialVelocity;
            double verifica = balance.flowArea * state.cells[i].pres * 98066.5 * delvel * (alfinter - alfinterL) / state.cells[i].dx;

            double vPot = 0.;
            if (i > 0)
                vPot = state.cells[i - 1].potTermo;
            state.cells[i].temp = ((balance.timeCoefficient / state.cells[i].dt) * state.cells[i].temp - (-balance.pressureTimeCoefficient * (state.cells[i].pres - state.cells[i].presini) * 98066.5 / state.cells[i].dt) + state.cells[i].dTdLCor * (-balance.temperatureSpatialCoefficient * dtdx + balance.pressureSpatialCoefficient * dpdx - cinetico - (balance.hydrostaticPower - 0. * verifica) + (vPot + state.cells[i].fonteCal) / balance.meanCellLength + sourceTerms.liquid + sourceTerms.gas + balance.heatFlux - latente) - (balance.rc - balance.rp) * (1 - balance.voidFraction) * state.cells[i].pres * 98066.5 * balance.flowArea * (state.cells[i].bet - state.cells[i].betini) / (balance.liquidDensity * state.cells[i].dt)) / (balance.timeCoefficient / state.cells[i].dt);

            if ((i < 148 && i > 144) && (*state.globals).lixo5 > 53879) {
                int para;
                para = 0.;
            }

            if (fabs(state.cells[i].temp - state.cells[i].tempini) / state.cells[i].dt > 10.) {
                state.cells[i].temp = state.cells[i].tempini +
                                 (fabs(state.cells[i].temp - state.cells[i].tempini) / (state.cells[i].temp - state.cells[i].tempini)) * 10 * state.cells[i].dt;
            } else if (fabs(state.cells[i].temp - state.cells[i].tempini) / state.cells[i].dt > 1. && fabs(state.cells[i].VTemper) > 10 * balance.referenceMixtureVelocity)
                state.cells[i].temp = state.cells[i].tempini;
            double tminimo = -50.;
            if (state.input.usaTabela == 1)
                tminimo = state.input.tabent.tmin + 1.;
            if (state.cells[i].temp < tminimo)
                state.cells[i].temp = tminimo;
            if (state.cells[i].temp > 200.)
                state.cells[i].temp = 200.;
        } else if (i == state.lastCell) {
            if ((*state.globals).chaverede == 0 || state.networkEndpoint == 1 || (*state.globals).chaveRedeParalela == 1)
                state.cells[i].temp = state.cells[i].calor.Textern1;
            else
                state.cells[i].temp = state.gasSurfaceTemperature;
        }
    } else {
        state.cells[i].temp = state.cells[i].calor.Textern1;
    }
}

void computeThermalMassTransfer(const ThermalState &state, int i) {

    double dx = state.cells[i].dx;
    double dxmed = 0.5 * (state.cells[i].dx + state.cells[i - 1].dx);
    double dia = state.cells[i].duto.a;
    double area = 0.25 * M_PI * dia * dia;
    double alfmed = state.cells[i].alf;
    double betmed = state.cells[i].bet;
    double ugsmed;

    if (state.cells[i].alf > (*state.globals).localtiny)
        ugsmed = state.cells[i].QG / area;
    else {
        ugsmed = 0.;
    }
    double ulsmed;
    if (state.cells[i].alf < 1. - (*state.globals).localtiny)
        ulsmed = state.cells[i].QL / area;
    else {
        ulsmed = 0.;
    }
    double rp = state.cells[i].flui.MasEspLiq(state.cells[i].pres, state.cells[i].temp);
    double rc = state.cells[i].fluicol.MasEspFlu(state.cells[i].pres, state.cells[i].temp);
    double rhol = (1. - betmed) * rp + betmed * rc;
    double rhog = state.cells[i].flui.MasEspGas(state.cells[i].pres, state.cells[i].temp);
    double cpl = (1. - betmed) * state.cells[i].flui.CalorLiq(state.cells[i].pres, state.cells[i].temp) + betmed * state.cells[i].fluicol.CalorLiq(state.cells[i].pres, state.cells[i].temp);
    double cvl = cpl;
    double cpg = state.cells[i].flui.CalorGas(state.cells[i].pres, state.cells[i].temp);
    double cvg = state.cells[i].flui.CalorGasVolMod(state.cells[i].pres, state.cells[i].temp);
    double jtl = (1. - betmed) * state.cells[i].flui.JTL(state.cells[i].pres, state.cells[i].temp) - betmed / rc;
    double jtg = state.cells[i].flui.JTG(state.cells[i].pres, state.cells[i].temp);
    double hidro = (rhol * ulsmed + rhog * ugsmed) * area * 9.82 * sin(state.cells[i].duto.teta);

    state.cells[i].calor.Tint = state.cells[i].temp;
    state.cells[i].calor.Vint = ugsmed + ulsmed;
    state.cells[i].calor.dt = state.cells[i].dt;
    double condliq = (1. - betmed) * state.cells[i].flui.CondLiq(state.cells[i].pres, state.cells[i].temp) + betmed * state.cells[i].fluicol.CondLiq(state.cells[i].pres, state.cells[i].temp);
    state.cells[i].calor.kint = condliq * (1 - alfmed) + state.cells[i].flui.CondGas(state.cells[i].pres, state.cells[i].temp) * alfmed;
    state.cells[i].calor.cpint = cpl * (1 - alfmed) + cpg * alfmed;
    state.cells[i].calor.rhoint = rhol * (1 - alfmed) + rhog * alfmed;
    double viscliq = (1. - betmed) * state.cells[i].flui.ViscOleo(state.cells[i].pres, state.cells[i].temp) + betmed * state.cells[i].fluicol.VisFlu(state.cells[i].pres, state.cells[i].temp);
    state.cells[i].calor.viscint = viscliq * (1 - alfmed) * 1.e-3 + state.cells[i].flui.ViscGas(state.cells[i].pres, state.cells[i].temp) * alfmed * 1.e-3;
    double fluxcal = state.cells[i].calor.transtrans();

    double coefTempo = (rhol * (1 - alfmed) * cvl + rhog * alfmed * cvg) * area;
    double coefPresTempo;
    coefPresTempo = -state.cells[i].flui.CalorGasPresMod(state.cells[i].pres, state.cells[i].temp) * (rhog * alfmed * area);
    double coefdxT = (rhol * ulsmed * cpl + rhog * ugsmed * cpg) * area;
    double coefdxP = (rhol * ulsmed * jtl + rhog * ugsmed * jtg) * area;
    double dpdx;
    if (i < state.lastCell)
        dpdx = 2. * (state.cells[i + 1].presaux - state.cells[i].pres) * 98066.5 / state.cells[i].dx;
    else
        dpdx = 2. * (state.cells[i].pres - state.cells[i].presaux) * 98066.5 / state.cells[i].dx;
    if (state.cells[i].acsr.tipo == 5 && state.cells[i].acsr.chk.AreaGarg <= (1e-3 + state.input.master1.razareaativ) * state.cells[i].duto.area)
        dpdx = 2. * (state.cells[i].pres - state.cells[i].presaux) * 98066.5 / state.cells[i].dx;
    state.cells[i].VTemper = coefdxT / coefTempo;
    double dtdx = (state.cells[i].temp - state.cells[i - 1].temp) / dxmed;
    if (i < state.lastCell)
        if (state.cells[i].VTemper < 0)
            dtdx = (state.cells[i + 1].temp - state.cells[i].temp) / (0.5 * (state.cells[i + 1].dx + state.cells[i].dx));
    if (state.cells[i].acsr.tipo == 5 && state.cells[i].acsr.chk.AreaGarg <= (1e-3 + state.input.master1.razareaativ) * state.cells[i].duto.area && state.cells[i].VTemper <= 0)
        dtdx = 0.;
    if (state.cells[i - 1].acsr.tipo == 5 && state.cells[i - 1].acsr.chk.AreaGarg <= (1e-3 + state.input.master1.razareaativ) * state.cells[i - 1].duto.area && state.cells[i].VTemper >= 0)
        dtdx = 0.;
    if ((i == 1 && state.cells[i].VTemper <= 0) || (i == state.lastCell && state.cells[i].VTemper <= 0))
        dtdx = 0.;
    if (state.cells[i].acsr.tipo == 4 && state.cells[i].acsr.bcs.freqnova > 1.) {
        dpdx = 0.;
    }
    if (state.cells[i].acsr.tipo == 8 && state.cells[i].acsr.bvol.freq > 1.) {
        dpdx = (state.cells[i].pres - state.cells[i - 1].pres) * 98066.5 / dxmed;
    }
    if (state.cells[i].acsr.tipo == 17 && state.cells[i].acsr.multibcs.freqnova > 1.) {
        dpdx = 0.;
    }

    double cinetico = 0;
    double ugmed0 = 0;
    double ulmed0 = 0;
    double ugmed = 0;
    double ulmed = 0.;
    if (state.cells[i].acsr.tipo == 0 && state.cells[i - 1].acsr.tipo == 0 && i > 1) {

        double dxCin = state.cells[i - 1].dx;
        double diaaux = state.cells[i - 1].duto.a;
        double areaaux = 0.25 * M_PI * diaaux * diaaux;

        double alf;
        if (state.cells[i].QL > 0)
            alf = state.cells[i - 1].alf;
        else
            alf = state.cells[i].alf;
        double alf0;
        if (state.cells[i - 1].QL > 0)
            alf0 = state.cells[i - 2].alf;
        else
            alf0 = state.cells[i - 1].alf;

        if (alf0 > 1e-3) {
            ugmed0 = state.cells[i - 1].QG / (areaaux);
            ugmed0 /= alf0;
        }

        if (alf0 < 1. - 1e-3) {
            ulmed0 = state.cells[i - 1].QL / (areaaux);
            ulmed0 /= (1. - alf0);
        }

        if (alf > 1e-3) {
            ugmed = ugsmed;
            ugmed /= alf;
        }

        if (alf < 1. - 1e-3) {
            ulmed = ulsmed;
            ulmed /= (1. - alf);
        }

        cinetico = (state.cells[i].MC - state.cells[i - 1].Mliqini) * ugmed * (ugmed - ugmed0) / dxCin + state.cells[i - 1].Mliqini * ulmed * (ulsmed - ulmed0) / dxCin;
    }

    TemperatureSourceTerms sourceTerms =
        computeThermalMassTransferSourceTerms(state, i, dx);

    double latente;
    if (state.latentHeatEnabled > 0 && state.input.flashCompleto == 0) {
        latente = interpolateLatentHeat(state, state.cells[i].pres, state.cells[i].temp);
    } else if (state.input.flashCompleto == 1)
        latente = (state.cells[i].flui.EntalpGas(state.cells[i].pres, state.cells[i].temp) -
                   state.cells[i].flui.EntalpLiq(state.cells[i].pres, state.cells[i].temp));
    else
        latente = 0;

    double alfinter;
    double alfinterL;
    if (ugsmed > 0) {
        alfinter = state.cells[i].alf;
        alfinterL = state.cells[i].alfL;
    } else {
        alfinter = state.cells[i].alfR;
        alfinterL = state.cells[i].alf;
    }
    double delvel;
    if (alfinter > (*state.globals).localtiny && alfinter < (1. - (*state.globals).localtiny))
        delvel = ugsmed / alfinter - ulsmed / (1. - alfinter);
    else if (alfinter > (*state.globals).localtiny)
        delvel = ugsmed;
    else
        delvel = ulsmed;
    double verifica = area * state.cells[i].pres * 98066.5 * delvel * (alfinter - alfinterL) / state.cells[i].dx;

    if (state.cells[i].temp < 0 || state.cells[i].temp > 100) {
        int para;
        para = 0.;
    }

    state.cells[i].FonteMudaFase = (-(coefTempo / state.cells[i].dt) * (state.cells[i].temp - state.cells[i].tempini) - (coefPresTempo * (state.cells[i].pres - state.cells[i].presini) * 98066.5 / state.cells[i].dt) - coefdxT * dtdx + coefdxP * dpdx - cinetico - (hidro - 0. * verifica) + state.cells[i - 1].potB / dxmed + sourceTerms.liquid + sourceTerms.gas + fluxcal - (rc - rp) * (1 - alfmed) * state.cells[i].pres * 98066.5 * area * (state.cells[i].bet - state.cells[i].betini) / (rhol * state.cells[i].dt)); // / (coefTempo / state.cells[i].dt);

    state.cells[i].FonteMudaFase /= latente;
}

namespace {

void initializeDistributedMassTransferInlet(
    const ThermalState &state, int cellIndex,
    double &previousLiquidDensity, double &previousOilVolumeFactor,
    double &previousSolutionGasRatio,
    double &previousSolutionGasPressureDerivative) {
    state.cells[0].transmassLini = state.cells[0].transmassL;
    state.cells[0].transmassL = 0.;
    state.sourceUpdater(cellIndex);
    state.cells[1].fontemassLLini = state.cells[1].fontemassLL;
    state.cells[1].fontemassLL = state.cells[0].fontemassLR;
    state.cells[1].fontemassCLini = state.cells[1].fontemassCL;
    state.cells[1].fontemassCL = state.cells[0].fontemassCR;
    state.cells[1].fontemassGLini = state.cells[1].fontemassGL;
    state.cells[1].fontemassGL = state.cells[0].fontemassGR;
    previousLiquidDensity =
        (1 - state.cells[cellIndex].bet) *
            state.cells[cellIndex].flui.MasEspLiq(
                state.cells[cellIndex].pres, state.cells[cellIndex].temp) +
        state.cells[cellIndex].bet *
            state.cells[cellIndex].fluicol.MasEspFlu(
                state.cells[cellIndex].pres, state.cells[cellIndex].temp);
    previousOilVolumeFactor = state.cells[cellIndex].flui.BOFunc(
        state.cells[cellIndex].pres, state.cells[cellIndex].temp);
    previousSolutionGasRatio = state.cells[cellIndex].flui.RS(
        state.cells[cellIndex].pres, state.cells[cellIndex].temp);
    if (state.input.flashCompleto != 2 || state.input.miniTabAtraso > 0) {
        double boL0 = state.cells[cellIndex].flui.BOFunc(
            state.cells[cellIndex].pres * 0.999,
            state.cells[cellIndex].temp);
        double rsL0 = state.cells[cellIndex].flui.RS(
            state.cells[cellIndex].pres * 0.999,
            state.cells[cellIndex].temp);
        previousSolutionGasPressureDerivative =
            (previousSolutionGasRatio / previousOilVolumeFactor -
             rsL0 / boL0) /
            (state.cells[cellIndex].pres * 0.001);
    } else {
        ProFlu flutemp = state.cells[cellIndex].flui;
        flutemp.atualizaPropComp(
            state.cells[cellIndex].pres * 0.999,
            state.cells[cellIndex].temp, flutemp.dCalculatedBeta,
            flutemp.oCalculatedLiqComposition,
            flutemp.oCalculatedVapComposition, state.input.pocinjec);
        double boL0 = flutemp.BOFunc(
            state.cells[cellIndex].pres * 0.999,
            state.cells[cellIndex].temp);
        double rsL0 = flutemp.RS(
            state.cells[cellIndex].pres * 0.999,
            state.cells[cellIndex].temp);
        previousSolutionGasPressureDerivative =
            (previousSolutionGasRatio / previousOilVolumeFactor -
             rsL0 / boL0) /
            (state.cells[cellIndex].pres * 0.001);
    }
    ProFlu flutemp = state.cells[cellIndex].flui;
}

struct DistributedMassTransferProperties {
    double downstreamWaterFraction;
    double upstreamWaterFraction;
    double cellWaterFraction;
    double liquidDensity;
    double gasDensity;
    double downstreamComposition;
    double upstreamComposition;
    double mixtureLiquidDensity;
    double downstreamOilVolumeFactor;
    double downstreamSolutionGasRatio;
    double downstreamSolutionGasPressureDerivative;
    double cellOilVolumeFactor;
    double cellSolutionGasRatio;
    double cellSolutionGasPressureDerivative;
    double cellSolutionGasTemperatureDerivative;
};

DistributedMassTransferProperties prepareDistributedMassTransferProperties(
    const ThermalState &state, int i, double tmed, ProFlu &flue,
    ProFlu &flud) {
    double fwd;
    double fwe;
    double boC = state.cells[i - 1].flui.BOFunc(
        state.cells[i - 1].pres, state.cells[i - 1].temp);
    double baC = state.cells[i - 1].flui.BAFunc(
        state.cells[i - 1].pres, state.cells[i - 1].temp);
    double fwC = state.cells[i - 1].flui.BSW * baC /
                 (boC + baC * state.cells[i - 1].flui.BSW -
                  state.cells[i - 1].flui.BSW * boC);
    if (state.cells[i].Mliqini < 0.) {
        flud = state.cells[i].flui;
        double bo1 = flud.BOFunc(state.cells[i].pres, state.cells[i].temp);
        double ba1 = flud.BAFunc(state.cells[i].pres, state.cells[i].temp);
        fwd = flud.BSW * ba1 / (bo1 + ba1 * flud.BSW - flud.BSW * bo1);
    } else {
        flud = state.cells[i - 1].flui;
        double bo1 = flud.BOFunc(
            state.cells[i - 1].pres, state.cells[i - 1].temp);
        double ba1 = flud.BAFunc(
            state.cells[i - 1].pres, state.cells[i - 1].temp);
        fwd = flud.BSW * ba1 / (bo1 + ba1 * flud.BSW - flud.BSW * bo1);
    }
    if (state.cells[i - 1].Mliqini < 0) {
        flue = state.cells[i - 1].flui;
        double bo0 = flue.BOFunc(
            state.cells[i - 1].pres, state.cells[i - 1].temp);
        double ba0 = flue.BAFunc(
            state.cells[i - 1].pres, state.cells[i - 1].temp);
        fwe = flue.BSW * ba0 / (bo0 + ba0 * flue.BSW - flue.BSW * bo0);

    } else {
        if (i > 1) {
            flue = state.cells[i - 2].flui;
            double bo0 = flue.BOFunc(
                state.cells[i - 2].pres, state.cells[i - 2].temp);
            double ba0 = flue.BAFunc(
                state.cells[i - 2].pres, state.cells[i - 2].temp);
            fwe = flue.BSW * ba0 /
                  (bo0 + ba0 * flue.BSW - flue.BSW * bo0);
        } else {
            flue = state.cells[i - 1].flui;
            double bo0 = flue.BOFunc(
                state.cells[i - 1].pres, state.cells[i - 1].temp);
            double ba0 = flue.BAFunc(
                state.cells[i - 1].pres, state.cells[i - 1].temp);
            fwe = flue.BSW * ba0 /
                  (bo0 + ba0 * flue.BSW - flue.BSW * bo0);
        }
    }

    double rl;
    double rg;
    double betI;

    // casoComp

    rl = flud.MasEspLiq(state.cells[i].presaux, tmed);
    if (state.cells[i].Mliqini < 0)
        betI = state.cells[i].bet; // testeBeta
    else
        betI = state.cells[i].betL;

    rg = flud.MasEspGas(state.cells[i].presaux, tmed);

    double betL = state.cells[i - 1].betL;
    if (state.cells[i - 1].Mliqini < 0)
        betL = state.cells[i - 1].bet; // testeBeta

    if (i > 0)
        betI = state.cells[i - 1].betPigD;
    if (state.cells[i].Mliqini < 0)
        betI = state.cells[i].betPigE; // testeBeta
    if (i > 1)
        betL = state.cells[i - 2].betPigD;
    if (state.cells[i - 1].Mliqini < 0)
        betL = state.cells[i - 1].betPigE; // testebeta

    double rhol = (1 - betI) * rl +
                  betI * state.cells[i].fluicol.MasEspFlu(
                             state.cells[i].presaux, tmed);
    double boR;
    double rsR;
    double boR0;
    double rsR0;
    if (state.input.flashCompleto != 2 || state.input.miniTabAtraso > 0) {
        boR = flud.BOFunc(state.cells[i].presaux, tmed);
        rsR = flud.RS(state.cells[i].presaux, tmed);
        boR0 = flud.BOFunc(state.cells[i].presaux * 0.999, tmed);
        rsR0 = flud.RS(state.cells[i].presaux * 0.999, tmed);
    } else {
        boR = flud.BOFunc(state.cells[i].presaux, tmed);
        rsR = flud.RS(state.cells[i].presaux, tmed);
        flud.atualizaPropComp(
            state.cells[i].presaux * 0.999, tmed, flud.dCalculatedBeta,
            flud.oCalculatedLiqComposition,
            flud.oCalculatedVapComposition, state.input.pocinjec);
        boR0 = flud.BOFunc(state.cells[i].presaux * 0.999, tmed);
        rsR0 = flud.RS(state.cells[i].presaux * 0.999, tmed);
    } // casoComp
    double DRsBoR =
        (rsR / boR - rsR0 / boR0) / (state.cells[i].presaux * 0.001);
    double boM;
    double rsM;
    double boM0;
    double rsM0;
    if (state.input.flashCompleto != 2 || state.input.miniTabAtraso > 0) {
        boM = state.cells[i - 1].flui.BOFunc(
            state.cells[i - 1].pres, state.cells[i - 1].temp);
        rsM = state.cells[i - 1].flui.RS(
            state.cells[i - 1].pres, state.cells[i - 1].temp);
        boM0 = state.cells[i - 1].flui.BOFunc(
            state.cells[i - 1].pres * 0.999, state.cells[i - 1].temp);
        rsM0 = state.cells[i - 1].flui.RS(
            state.cells[i - 1].pres * 0.999, state.cells[i - 1].temp);
    } else {
        boM = state.cells[i - 1].flui.BOFunc(
            state.cells[i - 1].pres, state.cells[i - 1].temp);
        rsM = state.cells[i - 1].flui.RS(
            state.cells[i - 1].pres, state.cells[i - 1].temp);
        ProFlu flutemp = state.cells[i - 1].flui;
        flutemp.atualizaPropComp(
            state.cells[i - 1].pres * 0.999,
            state.cells[i - 1].temp, flutemp.dCalculatedBeta,
            flutemp.oCalculatedLiqComposition,
            flutemp.oCalculatedVapComposition, state.input.pocinjec);
        boM0 = flutemp.BOFunc(
            state.cells[i - 1].pres * 0.999, state.cells[i - 1].temp);
        rsM0 = flutemp.RS(
            state.cells[i - 1].pres * 0.999, state.cells[i - 1].temp);
    } // casoComp
    double DRsBoM =
        (rsM / boM - rsM0 / boM0) / (state.cells[i - 1].pres * 0.001);
    double boM0T = 0.;
    double rsM0T = 0.;
    double DRsBoMT = 0;
    if (state.input.cicloAcopTerm == 1) {
        if (state.input.flashCompleto != 2 ||
            state.input.miniTabAtraso > 0) {
            boM0T = state.cells[i - 1].flui.BOFunc(
                state.cells[i - 1].pres,
                state.cells[i - 1].temp * 0.999);
            rsM0T = state.cells[i - 1].flui.RS(
                state.cells[i - 1].pres,
                state.cells[i - 1].temp * 0.999);
        } else {
            ProFlu flutemp = state.cells[i - 1].flui;
            flutemp.atualizaPropComp(
                state.cells[i - 1].pres,
                state.cells[i - 1].temp * 0.999,
                flutemp.dCalculatedBeta,
                flutemp.oCalculatedLiqComposition,
                flutemp.oCalculatedVapComposition, state.input.pocinjec);
            boM0T = flutemp.BOFunc(
                state.cells[i - 1].pres,
                state.cells[i - 1].temp * 0.999);
            rsM0T = flutemp.RS(
                state.cells[i - 1].pres,
                state.cells[i - 1].temp * 0.999);
        } // casoComp
        DRsBoMT = (rsM / boM - rsM0T / boM0T) /
                   (state.cells[i - 1].temp * 0.001);
    }

    return DistributedMassTransferProperties{
        .downstreamWaterFraction = fwd,
        .upstreamWaterFraction = fwe,
        .cellWaterFraction = fwC,
        .liquidDensity = rl,
        .gasDensity = rg,
        .downstreamComposition = betI,
        .upstreamComposition = betL,
        .mixtureLiquidDensity = rhol,
        .downstreamOilVolumeFactor = boR,
        .downstreamSolutionGasRatio = rsR,
        .downstreamSolutionGasPressureDerivative = DRsBoR,
        .cellOilVolumeFactor = boM,
        .cellSolutionGasRatio = rsM,
        .cellSolutionGasPressureDerivative = DRsBoM,
        .cellSolutionGasTemperatureDerivative = DRsBoMT,
    };
}

struct DistributedMassTransferCoefficients {
    double activeDerivative;
    double spatialCoupling;
    double flowArea;
};

DistributedMassTransferCoefficients updateDistributedMassTransferDerivatives(
    const ThermalState &state, int i, double fwC, const ProFlu &flud,
    double DRsBoM, double DRsBoMT) {
    double ativa = 1.;
    double limipres = 10;
    if (state.completeModel == 1)
        limipres = 0;
    if (state.cells[i - 1].pres < limipres ||
        state.massTransferModel != 0)
        ativa = 0.;
    double acop = 1.;
    if (i < 2 || i == state.lastCell)
        acop = 0;

    double A1 = state.cells[i].dutoL.area;
    state.cells[i - 1].ativaDeri = ativa;
    state.cells[i - 1].DTransDtp =
        ativa * A1 * (1. - state.cells[i - 1].alf) *
        (1. - state.cells[i - 1].bet) * (1. - fwC) * flud.Deng *
        1.225 * DRsBoM * (6.29 / 35.31467);
    state.cells[i].DTransDtpL = state.cells[i - 1].DTransDtp;
    if (i == state.lastCell) {
        double boM;
        double rsM;
        double boM0;
        double rsM0;
        if (state.input.flashCompleto != 2 ||
            state.input.miniTabAtraso > 0) {
            boM = state.cells[i].flui.BOFunc(
                state.cells[i].pres, state.cells[i].temp);
            rsM = state.cells[i].flui.RS(
                state.cells[i].pres, state.cells[i].temp);
            boM0 = state.cells[i].flui.BOFunc(
                state.cells[i].pres * 0.999, state.cells[i].temp);
            rsM0 = state.cells[i].flui.RS(
                state.cells[i].pres * 0.999, state.cells[i].temp);
        } else {
            boM = state.cells[i].flui.BOFunc(
                state.cells[i].pres, state.cells[i].temp);
            rsM = state.cells[i].flui.RS(
                state.cells[i].pres, state.cells[i].temp);
            ProFlu flutemp = state.cells[i].flui;
            flutemp.atualizaPropComp(
                state.cells[i].pres * 0.999, state.cells[i].temp,
                flutemp.dCalculatedBeta,
                flutemp.oCalculatedLiqComposition,
                flutemp.oCalculatedVapComposition, state.input.pocinjec);
            boM0 = flutemp.BOFunc(
                state.cells[i].pres * 0.999, state.cells[i].temp);
            rsM0 = flutemp.RS(
                state.cells[i].pres * 0.999, state.cells[i].temp);
        } // casoComp
        double DRsBoM =
            (rsM / boM - rsM0 / boM0) / (state.cells[i].pres * 0.001);
        state.cells[i].DTransDtp =
            ativa * A1 * (1. - state.cells[i].alf) *
            (1. - state.cells[i].bet) * (1. - fwC) * flud.Deng *
            1.225 * DRsBoM * (6.29 / 35.31467);
    }
    if (state.input.cicloAcopTerm == 1) {
        state.cells[i - 1].DTransDtT =
            ativa * A1 * (1. - state.cells[i - 1].alf) *
            (1. - state.cells[i - 1].bet) * (1. - fwC) * flud.Deng *
            1.225 * DRsBoMT * (6.29 / 35.31467);
        state.cells[i].DTransDtTL = state.cells[i - 1].DTransDtT;
        if (i == state.lastCell) {
            double boM;
            double rsM;
            double boM0;
            double rsM0;
            if (state.input.flashCompleto != 2 ||
                state.input.miniTabAtraso > 0) {
                boM = state.cells[i].flui.BOFunc(
                    state.cells[i].pres, state.cells[i].temp);
                rsM = state.cells[i].flui.RS(
                    state.cells[i].pres, state.cells[i].temp);
                boM0 = state.cells[i].flui.BOFunc(
                    state.cells[i].pres * 0.999, state.cells[i].temp);
                rsM0 = state.cells[i].flui.RS(
                    state.cells[i].pres * 0.999, state.cells[i].temp);
            } else {
                boM = state.cells[i].flui.BOFunc(
                    state.cells[i].pres, state.cells[i].temp);
                rsM = state.cells[i].flui.RS(
                    state.cells[i].pres, state.cells[i].temp);
                ProFlu flutemp = state.cells[i].flui;
                flutemp.atualizaPropComp(
                    state.cells[i].pres * 0.999, state.cells[i].temp,
                    flutemp.dCalculatedBeta,
                    flutemp.oCalculatedLiqComposition,
                    flutemp.oCalculatedVapComposition,
                    state.input.pocinjec);
                boM0 = flutemp.BOFunc(
                    state.cells[i].pres * 0.999, state.cells[i].temp);
                rsM0 = flutemp.RS(
                    state.cells[i].pres * 0.999, state.cells[i].temp);
            } // casoComp
            double DRsBoM =
                (rsM / boM - rsM0 / boM0) /
                (state.cells[i].pres * 0.001);
            state.cells[i].DTransDtT =
                ativa * A1 * (1. - state.cells[i].alf) *
                (1. - state.cells[i].bet) * (1. - fwC) * flud.Deng *
                1.225 * DRsBoMT * (6.29 / 35.31467);
        }
    }

    return DistributedMassTransferCoefficients{
        .activeDerivative = ativa,
        .spatialCoupling = acop,
        .flowArea = A1,
    };
}

void selectDistributedMassTransferModel(
    const ThermalState &state, int i, double &tmed, double &tmedL,
    double ABSjL) {
    tmed = state.cells[i - 1].temp;
    if (state.cells[i].VTemper < 0.)
        tmed = state.cells[i].temp;
    tmedL = state.cells[i - 1].tempL;
    if (state.cells[i - 1].VTemper < 0.)
        tmedL = state.cells[i - 1].temp;

    state.cells[i - 1].TMModel = state.massTransferModel;
    if (state.massTransferModel != 3) {
        if ((((state.cells[i - 1].alf < 0.001) ||
              (state.cells[i - 1].alf > 0.999) ||
              (state.cells[i - 1].bet > 0.999 &&
               state.cells[i - 1].alf < 0.999)) &&
             ABSjL < 0.1) ||
            state.cells[i - 1].flui.RGO >= (*state.globals).RGOMax)
            state.cells[i - 1].TMModel = 3;
        else if (state.cells[i - 1].estadoPig == 1)
            state.cells[i - 1].TMModel = 3;
        else if (state.cells[i - 1].acsr.tipo == 2 ||
                 state.cells[i - 1].acsr.tipo == 3 ||
                 state.cells[i - 1].acsr.tipo == 9 ||
                 state.cells[i - 1].acsr.tipo == 15 ||
                 state.cells[i - 1].acsr.tipo == 16)
            state.cells[i - 1].TMModel = 3;
        else if (i >= 2) {
            if (state.cells[i - 2].acsr.tipo == 5 &&
                (state.cells[i - 2].acsr.chk.AreaGarg <
                 (1e-3 + state.input.master1.razareaativ) *
                     state.cells[i - 2].duto.area))
                state.cells[i - 1].TMModel = 3;
            else if (state.cells[i - 2].acsr.tipo == 4 ||
                     state.cells[i - 2].acsr.tipo == 7 ||
                     state.cells[i - 2].acsr.tipo == 17)
                state.cells[i - 1].TMModel = 0;
        }
        if (state.input.flashCompleto == 2) {
            double titTeste = state.cells[i - 1].flui.FracMass(
                state.cells[i - 1].pres, state.cells[i - 1].temp);
            if (titTeste > 1.0 - 1e-2 || titTeste < 1e-2)
                state.cells[i - 1].TMModel = 3;
        }
    }
    if (state.cells[i - 1].TMModel == 0 &&
        state.cells[i - 1].alf <= (*state.globals).CritCond)
        state.cells[i - 1].TMModel = 1;
    if (state.cells[i - 1].TMModel == 0 && i == state.lastCell)
        state.cells[i - 1].TMModel = 1;
    if (i == state.lastCell)
        state.cells[i].FonteMudaFase = 0.;
}

void applyDistributedMassTransferModel(
    const ThermalState &state, int i, double &tmed, double &tmedL,
    double ABSjL, const ProFlu &flue, const ProFlu &flud, double fwd,
    double fwe, double fwC, double betI, double betL, double rhol,
    double boR, double rsR, double DRsBoR, double boM, double rsM,
    double ativa, double acop, double A1, double rhol0, double boL,
    double rsL, double DRsBoL) {
    selectDistributedMassTransferModel(state, i, tmed, tmedL, ABSjL);

    state.cells[i - 1].fontedissolv = 0.;

    state.cells[i - 1].transmassRini = state.cells[i - 1].transmassR;
    state.cells[i - 1].FonteMudaFaseini =
        state.cells[i - 1].FonteMudaFase;
    state.cells[i].DTransDxRini = state.cells[i].DTransDxR;
    state.cells[i].DTransDxLini = state.cells[i].DTransDxL;
    state.cells[i].DTransDt1ini = state.cells[i].DTransDt1;
    state.cells[i].DTransDt0ini = state.cells[i].DTransDt0;
    state.cells[i].DTransDxRpini = state.cells[i].DTransDxRp;
    state.cells[i].DTransDxLpini = state.cells[i].DTransDxLp;
    state.cells[i - 1].CoefDTLini = state.cells[i - 1].CoefDTL;
    state.cells[i - 1].coefTransBetini =
        state.cells[i - 1].coefTransBet;
    state.cells[i].transmassLini = state.cells[i].transmassL;

    state.cells[i].TMModelL = state.cells[i - 1].TMModel;
    state.cells[i - 1].FonteMudaFase = 0.;
    if (state.cells[i - 1].TMModel == 0 ||
        state.cells[i - 1].TMModel == 1) {

        state.cells[i - 1].transmassR =
            -(state.cells[i].QL * (1 - betI) * (flud.rDgD) * flud.Deng *
              1.225 * (1. - fwd) * rsR * (6.29 / 35.31467) / boR) +
            (state.cells[i - 1].QL * (1 - betL) * (flue.rDgD) *
             flue.Deng * 1.225 * (1. - fwe) * rsL * (6.29 / 35.31467) /
             boL);

        state.cells[i - 1].transmassR /= state.cells[i - 1].dx;
        state.cells[i - 1].transmassR += state.cells[i - 1].fontedissolv;
        state.cells[i].transmassL = state.cells[i - 1].transmassR;
        state.cells[i - 1].FonteMudaFase =
            state.cells[i - 1].transmassR -
            state.cells[i - 1].DTransDtp * state.cells[i - 1].d2pdt2 -
            state.cells[i - 1].DTransDtT * state.cells[i - 1].dTdtIni;
        if (state.cells[i - 1].TMModel == 1) {
            state.cells[i - 1].transmassR -=
                ativa * ((1. - state.cells[i - 1].bet) *
                         (1. - state.cells[i - 1].alf) * (1. - fwC) * A1 *
                         (state.cells[i - 1].flui.rDgD) *
                         state.cells[i - 1].flui.Deng * 1.225 * rsM *
                         (6.29 / 35.31467) / boM) /
                state.cells[i - 1].dt;
            state.cells[i - 1].transmassR +=
                ativa * ((1. - state.cells[i - 1].betini) *
                         (1. - state.cells[i - 1].alfini) * (1. - fwC) * A1 *
                         (state.cells[i - 1].flui.rDgD) *
                         state.cells[i - 1].flui.Deng * 1.225 * rsM *
                         (6.29 / 35.31467) / boM) /
                state.cells[i - 1].dt;

            state.cells[i].transmassL = state.cells[i - 1].transmassR;
            state.cells[i - 1].FonteMudaFase =
                state.cells[i - 1].transmassR;
        }
        if (state.cells[i - 1].TMModel == 0) {
            state.cells[i - 1].FonteMudaFase -=
                ativa * ((1. - state.cells[i - 1].bet) *
                         (1. - state.cells[i - 1].alf) * (1. - fwC) * A1 *
                         (state.cells[i - 1].flui.rDgD) *
                         state.cells[i - 1].flui.Deng * 1.225 * rsM *
                         (6.29 / 35.31467) / boM) /
                state.cells[i - 1].dt;
            state.cells[i - 1].FonteMudaFase +=
                ativa * ((1. - state.cells[i - 1].betini) *
                         (1. - state.cells[i - 1].alfini) * (1. - fwC) * A1 *
                         (state.cells[i - 1].flui.rDgD) *
                         state.cells[i - 1].flui.Deng * 1.225 * rsM *
                         (6.29 / 35.31467) / boM) /
                state.cells[i - 1].dt;
        }

        if (state.cells[i - 1].TMModel == 0) {
            state.cells[i].DTransDxR =
                -((1 - betI) * (flud.rDgD) * flud.Deng * 1.225 *
                  (1. - fwd) * rsR * (6.29 / 35.31467) / boR) /
                (rhol * state.cells[i - 1].dx);
            state.cells[i].DtransDxLinear =
                -acop * state.cells[i].QL *
                    ((1 - betI) * (flud.rDgD) * flud.Deng * 1.225 *
                     (1. - fwd) * (6.29 / 35.31467) *
                     (DRsBoR * state.cells[i].dpresaux)) /
                    (state.cells[i - 1].dx) +
                acop * state.cells[i].QL *
                    ((1 - betI) * (flud.rDgD) * flud.Deng * 1.225 *
                     (1. - fwd) * (6.29 / 35.31467) *
                     (DRsBoR * state.cells[i].presaux)) /
                    (state.cells[i - 1].dx);
            state.cells[i].DTransDxRp =
                -acop * state.cells[i].QL *
                ((1 - betI) * (flud.rDgD) * flud.Deng * 1.225 *
                 (1. - fwd) * (6.29 / 35.31467) * 0.5 * DRsBoR) /
                (state.cells[i - 1].dx);
            state.cells[i].DTransDxL =
                ((1 - betL) * (flue.rDgD) * flue.Deng * 1.225 *
                 (1. - fwe) * rsL * (6.29 / 35.31467) / boL) /
                (rhol0 * state.cells[i - 1].dx);
            state.cells[i].DtransDxLinear =
                state.cells[i].DtransDxLinear +
                acop * state.cells[i - 1].QL *
                    ((1 - betL) * (flue.rDgD) * flue.Deng * 1.225 *
                     (1. - fwe) * (6.29 / 35.31467) *
                     (DRsBoL * state.cells[i - 1].dpresaux)) /
                    (state.cells[i - 1].dx) -
                acop * state.cells[i - 1].QL *
                    ((1 - betL) * (flue.rDgD) * flue.Deng * 1.225 *
                     (1. - fwe) * (6.29 / 35.31467) *
                     (DRsBoL * state.cells[i - 1].presaux)) /
                    (state.cells[i - 1].dx);
            state.cells[i].DTransDxLp =
                acop * state.cells[i - 1].QL *
                ((1 - betL) * (flue.rDgD) * flue.Deng * 1.225 *
                 (1. - fwe) * (6.29 / 35.31467) * 0.5 * DRsBoL) /
                (state.cells[i - 1].dx);
            state.cells[i].DTransDt1 =
                -ativa * ((1. - fwC) * A1 *
                          (state.cells[i - 1].flui.rDgD) *
                          state.cells[i - 1].flui.Deng * 1.225 * rsM *
                          (6.29 / 35.31467) / boM);
            state.cells[i].DTransDt0 = -state.cells[i].DTransDt1;

            state.cells[i - 1].CoefDTR =
                -((1. - state.cells[i - 1].bet) * (1. - fwC) * A1 *
                  (state.cells[i - 1].flui.rDgD) *
                  state.cells[i - 1].flui.Deng * 1.225 * rsM *
                  (6.29 / 35.31467) / boM);
            state.cells[i - 1].CoefDTL = -state.cells[i - 1].CoefDTR;
            state.cells[i - 1].coefTransBet =
                ((1. - fwC) * A1 * (state.cells[i - 1].flui.rDgD) *
                 state.cells[i - 1].flui.Deng * 1.225 * rsM *
                 (6.29 / 35.31467) / boM);

            // state.cells[i-1].DTransDtp=((1. - state.cells[i - 1].bet) * (1. - state.cells[i - 1].alf) * (1. - fwC)*A1
            state.cells[i].transmassL -=
                ativa * ((1. - state.cells[i - 1].bet) *
                         (1. - state.cells[i - 1].alf) * (1. - fwC) * A1 *
                         (state.cells[i - 1].flui.rDgD) *
                         state.cells[i - 1].flui.Deng * 1.225 * rsM *
                         (6.29 / 35.31467) / boM) /
                state.cells[i - 1].dt;
            state.cells[i].transmassL +=
                ativa * ((1. - state.cells[i - 1].betini) *
                         (1. - state.cells[i - 1].alfini) * (1. - fwC) * A1 *
                         (state.cells[i - 1].flui.rDgD) *
                         state.cells[i - 1].flui.Deng * 1.225 * rsM *
                         (6.29 / 35.31467) / boM) /
                state.cells[i - 1].dt;

        } else {
            state.cells[i].DTransDxR = 0.;
            state.cells[i].DTransDxL = 0.;
            state.cells[i].DTransDt1 = 0.;
            state.cells[i].DTransDxRp = 0.;
            state.cells[i].DTransDxLp = 0.;
            if (state.input.desligaDeriTransMassDTemp == 1) {
                state.cells[i - 1].DTransDtT = 0;
                state.cells[i].DTransDtTL = 0.;
            }
            state.cells[i - 1].CoefDTR = 0.;
            state.cells[i - 1].CoefDTL = 0.;
            state.cells[i - 1].coefTransBet = 0.;
        }
    }
}

}  // namespace

void updateDistributedMassTransfer(const ThermalState &state) {
    // #pragma omp parallel for num_threads(state.input.nthrd)
    double rhol0 = 0.;
    double boL = 0.;
    double rsL = 0.;
    double DRsBoL = 0.;
    for (int i = 0; i <= state.lastCell; i++) {
        if (i != 0 && i != state.lastCell + 1) {

            if ((*state.globals).lixo5 >= 12264.7 && i == state.lastCell - 1) {
                int para;
                para = 0;
            }

            state.sourceUpdater(i);

            if (i < state.lastCell) {
                state.cells[i + 1].fontemassLLini = state.cells[i].fontemassLR;
                state.cells[i + 1].fontemassCLini = state.cells[i].fontemassCR;
                state.cells[i + 1].fontemassGLini = state.cells[i].fontemassGR;
                state.cells[i + 1].fontemassLL = state.cells[i].fontemassLR;
                state.cells[i + 1].fontemassCL = state.cells[i].fontemassCR;
                state.cells[i + 1].fontemassGL = state.cells[i].fontemassGR;
            }
            double razdx = state.cells[i - 1].dx / (state.cells[i].dx + state.cells[i].dxL);
            double razdxL = state.cells[i].dx / (state.cells[i - 1].dx + state.cells[i - 1].dxL);
            double tmed = razdx * state.cells[i].temp + (1 - razdx) * state.cells[i - 1].temp;
            double tmedL = razdxL * state.cells[i].tempL + (1 - razdxL) * state.cells[i - 1].tempL;
            tmed = state.cells[i - 1].temp;
            if (state.cells[i].VTemper < 0.)
                tmed = state.cells[i].temp;
            tmedL = state.cells[i - 1].tempL;
            if (state.cells[i - 1].VTemper < 0.)
                tmedL = state.cells[i - 1].temp;

            double dia = state.cells[i].duto.a;
            double area = 0.25 * M_PI * dia * dia;
            double ugsmed = (state.cells[i].QG) / (area);
            double ulsmed = state.cells[i].QL / (area);
            double j = ugsmed + ulsmed;
            double ABSjL = (fabs(state.cells[i - 1].QG) + fabs(state.cells[i - 1].QL)) / state.cells[i - 1].duto.area;

            ProFlu flue;
            ProFlu flud;
            DistributedMassTransferProperties properties =
                prepareDistributedMassTransferProperties(
                    state, i, tmed, flue, flud);
            double fwd = properties.downstreamWaterFraction;
            double fwe = properties.upstreamWaterFraction;
            double fwC = properties.cellWaterFraction;
            double rl = properties.liquidDensity;
            double rg = properties.gasDensity;
            double betI = properties.downstreamComposition;
            double betL = properties.upstreamComposition;
            double rhol = properties.mixtureLiquidDensity;
            double boR = properties.downstreamOilVolumeFactor;
            double rsR = properties.downstreamSolutionGasRatio;
            double DRsBoR =
                properties.downstreamSolutionGasPressureDerivative;
            double boM = properties.cellOilVolumeFactor;
            double rsM = properties.cellSolutionGasRatio;
            double DRsBoM = properties.cellSolutionGasPressureDerivative;
            double DRsBoMT =
                properties.cellSolutionGasTemperatureDerivative;
            DistributedMassTransferCoefficients coefficients =
                updateDistributedMassTransferDerivatives(
                    state, i, fwC, flud, DRsBoM, DRsBoMT);
            double ativa = coefficients.activeDerivative;
            double acop = coefficients.spatialCoupling;
            double A1 = coefficients.flowArea;

            applyDistributedMassTransferModel(
                state, i, tmed, tmedL, ABSjL, flue, flud, fwd, fwe,
                fwC, betI, betL, rhol, boR, rsR, DRsBoR, boM, rsM,
                ativa, acop, A1, rhol0, boL, rsL, DRsBoL);
            if (state.cells[i - 1].TMModel == -2) {
                double veltit;
                if (state.cells[i].alfL > (*state.globals).localtiny && betI < (1. - (*state.globals).localtiny))
                    veltit = (state.cells[i].QG * rg + state.cells[i].QL * (1. - betI) * rl) / (A1 * (state.cells[i].alfL * rg + (1. - state.cells[i].alfL) * (1. - betI) * rl));
                else
                    veltit = 0.;
                double tit = state.cells[i - 1].flui.FracMassHidra(state.cells[i - 1].pres, state.cells[i - 1].temp);
                double raz = 0.999;
                double dtit = (tit - state.cells[i - 1].flui.FracMassHidra(state.cells[i - 1].pres * raz, state.cells[i - 1].temp)) / ((1 - raz) * state.cells[i - 1].pres);
                double dpres;
                dpres = (state.cells[i].presaux - state.cells[i - 1].presaux) / state.cells[i].dxL;
                state.cells[i].transmassL = state.cells[i - 1].transmassR = 1 * (state.cells[i].alfL * rg + (1. - state.cells[i].alfL) * (1. - betI) * rl) * veltit * (dpres * dtit) * A1;

                state.cells[i].DTransDxR = 0.;
                state.cells[i].DTransDxL = 0.;
                state.cells[i].DTransDt1 = 0.;
                state.cells[i].DTransDt0 = 0.;
                if (state.input.desligaDeriTransMassDTemp == 1) {
                    state.cells[i - 1].DTransDtT = 0;
                    state.cells[i].DTransDtTL = 0.;
                }
                state.cells[i - 1].CoefDTR = 0.;
                state.cells[i - 1].CoefDTL = 0.;
                state.cells[i - 1].coefTransBet = 0.;
            }
            if (state.cells[i - 1].TMModel == 3) {
                state.cells[i].transmassL = state.cells[i - 1].transmassR = 0.;
                state.cells[i].DTransDxR = 0.;
                state.cells[i].DTransDxL = 0.;
                state.cells[i].DTransDt1 = 0.;
                state.cells[i].DTransDt0 = 0.;
                state.cells[i].DTransDxRp = 0.;
                state.cells[i].DTransDxLp = 0.;
                state.cells[i - 1].CoefDTR = 0.;
                state.cells[i - 1].CoefDTL = 0.;
                state.cells[i - 1].coefTransBet = 0.;
                if (state.input.desligaDeriTransMassDTemp == 1) {
                    state.cells[i - 1].DTransDtT = 0;
                    state.cells[i].DTransDtTL = 0.;
                }
            }
            if (state.cells[i - 1].transmassR > 0 && (state.cells[i].alfL > (1. - (*state.globals).localtiny) || state.cells[i].betL > (1. - (*state.globals).localtiny))) {
                state.cells[i].transmassL = state.cells[i - 1].transmassR = -(*state.globals).localtiny;
                state.cells[i].DTransDxR = 0.;
                state.cells[i].DTransDxL = 0.;
                state.cells[i].DTransDt1 = 0.;
                state.cells[i].DTransDt0 = 0.;
                state.cells[i].DTransDxRp = 0.;
                state.cells[i].DTransDxLp = 0.;
                state.cells[i - 1].CoefDTR = 0.;
                state.cells[i - 1].CoefDTL = 0.;
                state.cells[i - 1].coefTransBet = 0.;
                if (state.input.desligaDeriTransMassDTemp == 1) {
                    state.cells[i - 1].DTransDtT = 0;
                    state.cells[i].DTransDtTL = 0.;
                }
            }
            if (state.cells[i - 1].transmassR < 0 && state.cells[i].alfL < (*state.globals).localtiny) {
                state.cells[i].transmassL = state.cells[i - 1].transmassR = (*state.globals).localtiny;
                state.cells[i].DTransDxR = 0.;
                state.cells[i].DTransDxL = 0.;
                state.cells[i].DTransDt1 = 0.;
                state.cells[i].DTransDt0 = 0.;
                state.cells[i].DTransDxRp = 0.;
                state.cells[i].DTransDxLp = 0.;
                state.cells[i - 1].CoefDTR = 0.;
                state.cells[i - 1].CoefDTL = 0.;
                state.cells[i - 1].coefTransBet = 0.;
                if (state.input.desligaDeriTransMassDTemp == 1) {
                    state.cells[i - 1].DTransDtT = 0;
                    state.cells[i].DTransDtTL = 0.;
                }
            }
            rhol0 = rhol;
            boL = boR;
            rsL = rsR;
            DRsBoL = DRsBoR;

        } else if (i == 0)
            initializeDistributedMassTransferInlet(
                state, i, rhol0, boL, rsL, DRsBoL);
    }
}

namespace {

void selectAndApplyInteriorFlowRegime(
    const ThermalState &state, int i, Vcr<int> &bif, double ugs,
    double uls, double ugs0, double uls0, double ugs1, double uls1,
    double rgR, double rlR, double rg, double rl, double amed) {
    bif[i] = 1;

    if ((*state.globals).lixo5 > 29900) {
        int para;
        para = 0;
    }

    if (state.cells[i - 1].alfPigD <= (*state.globals).localtiny && state.cells[i].alfPigE <= (*state.globals).localtiny && state.cells[i].fontemassGL <= (*state.globals).localtiny * 1e-5 && state.cells[i].fontemassGR <= (*state.globals).localtiny * 1e-5) {
        state.cells[i].term1 = 1.;
        state.cells[i].term2 = 0.;
        bif[i] = 0;
        state.cells[i].c0 = 1;
        state.cells[i].ud = 0;
        state.cells[i].arranjo = 0;
    } else if (state.cells[i - 1].alfPigD >= (1. - (*state.globals).localtiny) && state.cells[i].alfPigE >= (1. - (*state.globals).localtiny) && ((state.cells[i].fontemassLL + state.cells[i].fontemassCL) <= (*state.globals).localtiny * 1e-5 && (state.cells[i].fontemassLR + state.cells[i].fontemassCR) <= (*state.globals).localtiny * 1e-5)) {
        state.cells[i].term1 = 0.;
        state.cells[i].term2 = 0.;
        bif[i] = 0;
        state.cells[i].c0 = 1;
        state.cells[i].ud = 0;
        state.cells[i].arranjo = 0;
    } else if (state.cells[i - 1].acsr.tipo == 5 && state.cells[i - 1].acsr.chk.AreaGarg <= (1e-3)) {
        state.cells[i].term1 = 0.;
        state.cells[i].term2 = 0.;
        bif[i] = 0;
        state.cells[i].c0 = 1;
        state.cells[i].ud = 0;
        state.cells[i].arranjo = 0;
    }

    else if (ugs >= 0 && state.cells[i - 1].alfPigD <= (*state.globals).localtiny && (state.cells[i].fontemassGL <= (*state.globals).localtiny * 1e-5 && state.cells[i].fontemassGR <= (*state.globals).localtiny * 1e-5)) {
        state.cells[i].term1 = 1.;
        state.cells[i].term2 = 0.;
        bif[i] = 0;
        state.cells[i].c0 = 1;
        state.cells[i].ud = 0;
        state.cells[i].arranjo = 0;
        if (fabs(ugs) <= 1e-15 && state.cells[i].alfPigE > (1. - (*state.globals).localtiny) && uls < 0 && uls0 < 0 && state.cells[i].duto.teta > 0) {
            state.cells[i].term1 = 0.;
            state.cells[i].term2 = 0.;
            bif[i] = 0;
        }
        if (fabs(ugs) <= 1e-15 && state.cells[i].alfPigE > (1. - 10 * (*state.globals).localtiny) && state.cells[i].duto.teta < 0) {
            state.cells[i].term1 = 0.;
            state.cells[i].term2 = 0.;
            bif[i] = 1;
        }
        if (fabs(ugs) <= 1e-15 && state.cells[i].alf >= state.cells[i + 1].alf && ugs1 < 0) {
            bif[i] = 1;
        }
    } else if (ugs <= 0 && state.cells[i].alfPigE <= (*state.globals).localtiny && (state.cells[i].fontemassGL <= (*state.globals).localtiny * 1e-5 && state.cells[i].fontemassGR <= (*state.globals).localtiny * 1e-5)) {
        state.cells[i].term1 = 1.;
        state.cells[i].term2 = 0.;
        bif[i] = 0;
        state.cells[i].c0 = 1;
        state.cells[i].ud = 0;
        state.cells[i].arranjo = 0;
        if (fabs(ugs) <= 1e-15 && state.cells[i - 1].alfPigD > (1. - (*state.globals).localtiny) && uls > 0 && uls1 > 0) {
            state.cells[i].term1 = 0.;
            state.cells[i].term2 = 0.;
            bif[i] = 0;
        }
        if (fabs(ugs) <= 1e-15 && state.cells[i].alfL >= state.cells[i - 1].alfL && ugs0 > 0) {
            bif[i] = 1;
        }
        if (fabs(ugs) <= 1e-15 && fabs(uls) <= 1e-15 && state.cells[i].alf < state.cells[i - 1].alf && state.cells[i].duto.teta > 0) {
            bif[i] = 1;
        }
    } else if (uls >= 0 && state.cells[i - 1].alfPigD >= 1. - 1 * (*state.globals).localtiny && ((state.cells[i].fontemassLL + state.cells[i].fontemassCL) <= (*state.globals).localtiny * 1e-5 && (state.cells[i].fontemassLR + state.cells[i].fontemassCR) <= (*state.globals).localtiny * 1e-5)) {
        state.cells[i].term1 = 0.;
        state.cells[i].term2 = 0.;
        bif[i] = 0;
        state.cells[i].c0 = 1;
        state.cells[i].ud = 0;
        state.cells[i].arranjo = 0;
        if (fabs(uls) <= 1e-15 && state.cells[i].alfPigE < (*state.globals).localtiny && uls1 < 0) { // ATENCAO!!!!!!!!!!!!!!! não teria de ser bifásico, mono-liq só seo ângulo fosse negativo, não?
            state.cells[i].term1 = 1.;
            state.cells[i].term2 = 0.;
            bif[i] = 0;
        } else if (fabs(uls1) < (*state.globals).localtiny * 1e-5) {
            if (fabs(uls) < (*state.globals).localtiny * 1e-5 && state.cells[i].alfPigE < (*state.globals).localtiny && state.cells[i].fontemassGR >= (*state.globals).localtiny * 1e-5) { // ATENCAO!!!!!!!!!!!!!!!  sem sentido isto aqui
                state.cells[i].term1 = 1.;
                state.cells[i].term2 = 0.;
                bif[i] = 0;
            }
            if (fabs(uls) < (*state.globals).localtiny * 1e-5 && state.cells[i].alfPigE < (*state.globals).localtiny && state.cells[i].duto.teta >= 0) { // ATENCAO!!!!!!!!!!!!!!! alteracao 11/08/24, adicionado
                bif[i] = 1;
            }
            if (fabs(uls) < (*state.globals).localtiny * 1e-5 && state.cells[i].alfPigE < (*state.globals).localtiny && state.cells[i].duto.teta < 0) { // ATENCAO!!!!!!!!!!!!!!! alteracao 11/08/24, adicionado
                state.cells[i].term1 = 0.;
                state.cells[i].term2 = 0.;
                bif[i] = 0;
            }
        } else if ((fabs(uls) < 1e-15 && (uls1 < 0 || state.cells[i].duto.teta > 0) // ATENCAO!!!!!!!!!!!!!!! alteracao 11/08/24, estava || mudado para &&
                    && ((state.cells[i].alfPigE <= (1 - 10 * (*state.globals).localtiny + .0 * state.cells[i].alfPigER) &&
                         state.cells[i].alfPigER < 1 - 1 * (*state.globals).localtiny) ||
                        state.cells[i].alfPigE <= 0.7)))
            bif[i] = 1;

    } else if (uls <= 0 && state.cells[i].alfPigE >= 1. - (*state.globals).localtiny && ((state.cells[i].fontemassLL + state.cells[i].fontemassCL) <= (*state.globals).localtiny * 1e-5 && (state.cells[i].fontemassLR + state.cells[i].fontemassCR) <= (*state.globals).localtiny * 1e-5)) {
        state.cells[i].term1 = 0.;
        state.cells[i].term2 = 0.;
        bif[i] = 0;
        state.cells[i].c0 = 1;
        state.cells[i].ud = 0;
        state.cells[i].arranjo = 0;
        if (fabs(uls) <= (*state.globals).localtiny * 1e-5 && state.cells[i - 1].alfPigD < (*state.globals).localtiny && uls0 > 0) {
            state.cells[i].term1 = 1.;
            state.cells[i].term2 = 0.;
            bif[i] = 0;
        }

        if (fabs(uls0) < (*state.globals).localtiny * 1e-5) {
            if (fabs(uls) < (*state.globals).localtiny * 1e-5 && state.cells[i - 1].alfPigD < (*state.globals).localtiny && state.cells[i - 1].fontemassGR >= (*state.globals).localtiny * 1e-5) {
                state.cells[i].term1 = 1.;
                state.cells[i].term2 = 0.;
                bif[i] = 0;
            }
        } else if ((i > 1 && fabs(uls) < (*state.globals).localtiny * 1e-5 && state.cells[i].duto.teta < 0.95 * M_PI / 2. && uls0 > 0 && ((state.cells[i - 1].alfPigD <= 1.0 * state.cells[i - 2].alfPigD && state.cells[i - 2].alfPigD < 0.99) || state.cells[i - 1].alfPigD < 0.7)) && ugs >= 0.)
            bif[i] = 1;
        else if ((i > 1 && fabs(uls) < (*state.globals).localtiny * 1e-5 && state.cells[i].duto.teta >= 0.95 * M_PI / 2. && uls0 > 0 && ((state.cells[i - 1].alfPigD <= 1.0 * state.cells[i - 2].alfPigD && state.cells[i - 2].alfPigD < 0.99) || state.cells[i - 1].alfPigD < 0.7)) && ugs >= 0)
            bif[i] = 1;
        else {
            if ((fabs(uls) < (*state.globals).localtiny * 1e-5 && state.cells[i].duto.teta < 0.95 * M_PI / 2. && uls0 > 0 && ((state.cells[i - 1].alfPigD <= 1.0 * state.cells[i - 1].alfL) || state.cells[i - 1].alfPigD < 0.7)) && ugs >= 0)
                bif[i] = 1;
            else if ((fabs(uls) < (*state.globals).localtiny * 1e-5 && state.cells[i].duto.teta >= 0.95 * M_PI / 2. && uls0 > 0 && ((state.cells[i - 1].alfPigD <= 1.0 * state.cells[i - 1].alfL) || state.cells[i - 1].alfPigD < 0.7)) && ugs >= 0)
                bif[i] = 1;
            else if ((fabs(uls) < (*state.globals).localtiny * 1e-5 && state.cells[i].duto.teta >= 0.95 * M_PI / 2. && ((state.cells[i - 1].alfPigD <= 1.0 * state.cells[i - 1].alfL) || state.cells[i - 1].alfPigD < 0.7)) && ugs > 0)
                bif[i] = 1;
        }
    }

    if (bif[i] == 1) {

        double c0;
        double ud;
        double alfmed;

        alfmed = state.cells[i - 1].alfPigD;
        if (ugs < 0)
            alfmed = state.cells[i].alfPigE;
        c0 = 1.2;
        double dmed = state.cells[i].duto.a;
        if (state.cells[i].MC >= 0)
            dmed = state.cells[i].dutoL.a;
        double sinal = 1.;
        if (state.cells[i].duto.teta < 0.)
            sinal = -1.;
        ud = sinal * 0.32 * sqrt(9.82 * dmed);
        if (fabs(rgR) / rlR > 0.9) {
            c0 = 1.;
            ud = 0.;
        }
        if (fabs(state.cells[i].QG / (0.25 * M_PI * dmed * dmed * alfmed)) > 100. ||
            fabs(state.cells[i].QL / (0.25 * M_PI * dmed * dmed * (1. - alfmed))) > 100.) {
            c0 = 1.;
            ud = 0.;
        } else
            state.closureUpdater.instantaneous(i, c0, ud);
        state.cells[i].c0 = c0;
        state.cells[i].ud = ud;
        if (i == state.lastCell) {
            double num = (1. - alfmed * c0);
            double den = 1 + c0 * alfmed * ((rg / rl) - 1.);
            state.cells[i].term1 = num / den;
            state.cells[i].term2 = (-amed * alfmed * rg * ud) / den;
            double jlTeste0 = (ugs - alfmed * ud) / (alfmed * c0) - ugs;
            double jlTeste = (ugs + uls) * (1. - c0 * alfmed) - alfmed * ud;
            if ((jlTeste > 0. || jlTeste0 > 0.) && state.cells[i - 1].alfPigD > 1 - 1e-15) {
                state.cells[i].term1 = 0.;
                state.cells[i].term2 = 0.;
            }
            if ((jlTeste < 0. || jlTeste0 < 0.) && state.cells[i].alfPigE > 1 - 1e-15) {
                state.cells[i].term1 = 0.;
                state.cells[i].term2 = 0.;
            }
            if (state.cells[i - 1].acsr.tipo == 5 && state.cells[i - 1].acsr.chk.AreaGarg <= (1e-3 + state.input.master1.razareaativ) * state.cells[i - 1].duto.area) {
                state.cells[i].term1 = 0.;
                state.cells[i].term2 = 0.;
            }
        }
    }
}

void updateInteriorFlowPartitionCell(
    const ThermalState &state, int i, Vcr<int> &bif, Vcr<int> &valv) {
    if (i < state.lastCell) {
        state.cells[i - 1].alfR = state.cells[i + 1].alfL = state.cells[i].alf;
        state.cells[i - 1].betR = state.cells[i + 1].betL = state.cells[i].bet;
    } else {
        state.cells[i - 1].alfR = state.cells[i].alf;
        state.cells[i - 1].betR = state.cells[i].bet;
    }
    valv[i] = 1;
    if (state.cells[i - 1].acsr.tipo == 5 || state.cells[i - 1].acsr.tipo == 8) {
        if ((*state.cells[i].acsrL).tipo == 5 && (*state.cells[i].acsrL).chk.AreaGarg <= (1e-3 + state.input.master1.razareaativ) * state.cells[i - 1].duto.area)
            valv[i] = 0;
        if ((*state.cells[i].acsrL).tipo == 8 && fabs((*state.cells[i].acsrL).bvol.freq) > 1)
            valv[i] = 0;
    }
    if (valv[i] == 1) {
        double razdx = state.cells[i].dxL / (state.cells[i].dx + state.cells[i].dxL);
        double pmed;
        if (i < state.lastCell)
            pmed = state.cells[i].presaux;
        else
            pmed = state.cells[i].pres;
        double tmed;
        if (i < state.lastCell)
            tmed = state.cells[i].temp * razdx + state.cells[i - 1].temp * (1. - razdx);
        else
            tmed = state.gasSurfaceTemperature;
        if (state.cells[i].VTemper < 0.) {
            if (i < state.lastCell)
                tmed = state.cells[i].temp;
            else
                tmed = state.gasSurfaceTemperature;
        }
        double betI = state.cells[i - 1].betPigD;
        double rl;
        if (state.cells[i].QL < 0.) { // testeBeta
            betI = state.cells[i].betPigE;
            rl = (1 - betI) * state.cells[i].rpCi + betI * state.cells[i].rcCi;
        } else {
            betI = state.cells[i - 1].betPigD;
            rl = (1 - betI) * state.cells[i - 1].rpCi + betI * state.cells[i - 1].rcCi;
            // viscl1 = (1 - betI) * state.cells[i - 1].flui.ViscOleo(pmed, tmed)
            // tensup1 = (1 - betI) * state.cells[i - 1].flui.TensSuper(pmed, tmed)
        }
        double rg;
        double amed;
        double hns;
        if (state.cells[i].QG >= 0) {
            amed = state.cells[i].dutoL.area;
            hns = 1. - state.cells[i].alfL;
            rg = state.cells[i - 1].rgCi;
        } else {
            rg = state.cells[i].rgCi;
            amed = state.cells[i].duto.area;
            hns = 1. - state.cells[i].alf;
        }
        double ugs = state.cells[i].QG / (amed);
        double uls = state.cells[i].QL / (amed);
        double dia1 = state.cells[i].duto.a;
        if (ugs >= 0)
            dia1 = state.cells[i - 1].duto.a;

        double rmed = hns * rl + (1 - hns) * rg;
        double ang = state.cells[i].duto.teta;
        if (i >= 2) {
            if (state.cells[i - 2].acsr.tipo == 5 && state.cells[i - 2].acsr.chk.AreaGarg <= (1e-3)) {
                if (state.cells[i].QG >= 0)
                    ang = state.cells[i].duto.teta;
                else
                    ang = state.cells[i].dutoR.teta;
            } else {
                if (state.cells[i].QG >= 0)
                    ang = state.cells[i].dutoL.teta;
                else
                    ang = state.cells[i].duto.teta;
            }
        }
        double sinal = 1.;
        if (ang < 0.)
            sinal = -1.;

        double amedL = state.cells[i].dutoL.area;
        double razdxL = state.cells[i - 1].dxL / (state.cells[i - 1].dx + state.cells[i - 1].dxL);
        double betIL = 0.;
        if (i < 2)
            betIL = state.cells[i - 1].betL;
        else
            betIL = state.cells[i - 2].betPigD;
        if (state.cells[i - 1].QL < 0.)
            betIL = state.cells[i - 1].betPigE; // testeBeta
        // betIL = state.cells[i - 1].betPigE;        //duvidabeta
        double rgL = state.cells[i].rgLi;
        double rlL = (1 - betIL) * state.cells[i].rpLi + betIL * state.cells[i].rcLi;
        double ugs0 = (state.cells[i].ML - state.cells[i].MliqiniL) / (rgL * amedL);
        double uls0 = (state.cells[i].MliqiniL) / (rlL * amedL);

        double amedR = state.cells[i].dutoR.area;
        double razdxR = state.cells[i].dxR / (state.cells[i].dxR + state.cells[i].dx);
        double pmedR = state.cells[i].presauxR;
        double tmedR = state.cells[i].temp * razdxR + state.cells[i].tempR * (1. - razdxR);
        double betIR = state.cells[i].betPigD;
        if (state.cells[i].QLR < 0.) { // testeBeta
            if (i > state.lastCell - 2)
                betIR = state.cells[i].betR;
            else
                betIR = state.cells[i + 1].betPigE;
        }
        double rgR = state.cells[i].rgRi;
        double rlR = (1 - betIR) * state.cells[i].rpRi + betIR * state.cells[i].rcRi;
        double ugs1 = (state.cells[i].MR - state.cells[i].MliqiniR) / (rgR * amedR);
        double uls1 = (state.cells[i].MliqiniR) / (rlR * amedR);

        selectAndApplyInteriorFlowRegime(
            state, i, bif, ugs, uls, ugs0, uls0, ugs1, uls1,
            rgR, rlR, rg, rl, amed);
    } else {
        state.cells[i].c0 = 1.;
        state.cells[i].ud = 0.;
        state.cells[i].term1 = 0.;
        state.cells[i].term2 = 0.;
        state.cells[i].term1L = state.cells[i - 1].term1;
        state.cells[i].term2L = state.cells[i - 1].term2;
        state.cells[i - 1].term1R = state.cells[i].term1;
        state.cells[i - 1].term2R = state.cells[i].term2;
    }
}

void updateOutletBoundaryFlowPartition(
    const ThermalState &state, int i, Vcr<int> &bif) {
    state.cells[state.lastCell - 1].alfR = state.cells[state.lastCell].alf;
    state.cells[state.lastCell].alfR = state.cells[state.lastCell].alf;
    state.cells[state.lastCell - 1].betR = state.cells[state.lastCell].bet;
    state.cells[state.lastCell].betR = state.cells[state.lastCell].bet;

    double razdx = state.cells[i].dxL / (state.cells[i].dx + state.cells[i].dxL);
    double pmed = state.cells[i].presaux;
    double tmed = state.cells[i].temp * razdx + state.cells[i - 1].temp * (1. - razdx);
    double betI = state.cells[i].betL;
    if (state.cells[i].QL < 0.)
        betI = state.cells[i].bet; // testeBeta
    // betI = state.cells[i].bet;            //duvidabeta
    double rg = state.cells[i].flui.MasEspGas(pmed, tmed);
    double rl = (1 - betI) * state.cells[i].flui.MasEspLiq(pmed, tmed) + betI * state.cells[i].fluicol.MasEspFlu(pmed, tmed);
    double amed = state.cells[i].duto.area;
    if (state.cells[i].MC >= 0)
        amed = state.cells[i].dutoL.area;
    double ugs = state.cells[i].QG / (amed);
    double uls = state.cells[i].QL / (amed);

    double amedL = state.cells[i].dutoL.area;
    double razdxL = state.cells[i - 1].dxL / (state.cells[i - 1].dx + state.cells[i - 1].dxL);
    double pmedL = state.cells[i - 1].presaux;
    double tmedL = state.cells[i - 1].temp * razdxL + state.cells[i - 1].tempL * (1. - razdxL);
    double betIL = state.cells[i - 1].betL;
    if (state.cells[i - 1].QL < 0.)
        betIL = state.cells[i - 1].bet; // testeBeta
    // betIL = state.cells[i - 1].bet;            //duvidabeta
    double rgL = state.cells[i].flui.MasEspGas(pmedL, tmedL);
    double rlL = (1 - betIL) * state.cells[i].flui.MasEspLiq(pmedL, tmedL) + betIL * state.cells[i].fluicol.MasEspFlu(pmedL, tmedL);
    double uls0 = (state.cells[i].MliqiniL) / (rlL * amedL);

    double amedR = state.cells[i].dutoR.area;
    double razdxR = state.cells[i].dxR / (state.cells[i].dxR + state.cells[i].dx);
    double pmedR = state.cells[i].pres;
    double tmedR = state.cells[i].temp * razdxR + state.cells[i].tempR * (1. - razdxR);
    double betIR = state.cells[i].bet;
    double rgR = state.cells[i].flui.MasEspGas(pmedR, tmedR);
    double rlR = (1 - betIR) * state.cells[i].flui.MasEspLiq(pmedR, tmedR) + betIR * state.cells[i].fluicol.MasEspFlu(pmedR, tmedR);
    double uls1 = (state.cells[i].MliqiniR) / (rlR * amedR);

    bif[i] = 1;

    if (state.cells[i].alfL <= (*state.globals).localtiny && state.cells[i].alf <= (*state.globals).localtiny && state.cells[i].fontemassGL <= 0 && state.cells[i].fontemassGR <= 0) {
        state.cells[i].term1 = 1.;
        state.cells[i].term2 = 0.;
        bif[i] = 0;
        state.cells[i].c0 = 1. - (*state.globals).localtiny;
        state.cells[i].ud = 0.;
        state.cells[i].arranjo = 0;
    } else if (state.cells[i].alfL >= (1. - (*state.globals).localtiny) && state.cells[i].alf >= (1. - (*state.globals).localtiny) && state.cells[i].fontemassLL <= 0 && state.cells[i].fontemassLR <= 0) {
        state.cells[i].term1 = 0.;
        state.cells[i].term2 = 0.;
        bif[i] = 0;
        state.cells[i].c0 = 1. - (*state.globals).localtiny;
        state.cells[i].ud = 0.;
        state.cells[i].arranjo = 0;
    } else if (ugs > 0 && state.cells[i].alfL <= (*state.globals).localtiny && (state.cells[i].fontemassGL <= 0. && state.cells[i].fontemassGR <= 0.)) {
        state.cells[i].term1 = 1.;
        state.cells[i].term2 = 0.;
        bif[i] = 0;
        state.cells[i].c0 = 1. - (*state.globals).localtiny;
        state.cells[i].ud = 0.;
        state.cells[i].arranjo = 0;
        if (fabs(ugs) <= 1e-15 && state.cells[i].alf > (1. - (*state.globals).localtiny) && uls < 0 && uls0 < 0) {
            state.cells[i].term1 = 0.;
            state.cells[i].term2 = 0.;
            bif[i] = 0;
        }
    } else if (ugs < 0 && state.cells[i].alf <= (*state.globals).localtiny && (state.cells[i].fontemassGL <= 0. && state.cells[i].fontemassGR <= 0.)) {
        state.cells[i].term1 = 1.;
        state.cells[i].term2 = 0.;
        bif[i] = 0;
        state.cells[i].c0 = 1. - (*state.globals).localtiny;
        state.cells[i].ud = 0.;
        state.cells[i].arranjo = 0;
        if (fabs(ugs) <= 1e-15 && state.cells[i].alfL > (1. - (*state.globals).localtiny) && uls > 0 && uls1 > 0) {
            state.cells[i].term1 = 0.;
            state.cells[i].term2 = 0.;
            bif[i] = 0;
        }
    } else if (uls >= 0 && state.cells[i - 1].alf >= 1. - (*state.globals).localtiny && ((state.cells[i].fontemassLL + state.cells[i].fontemassCL) <= (*state.globals).localtiny * 1e-5 && (state.cells[i].fontemassLR + state.cells[i].fontemassCR) <= (*state.globals).localtiny * 1e-5)) {
        state.cells[i].term1 = 0.;
        state.cells[i].term2 = 0.;
        bif[i] = 0;
        state.cells[i].c0 = 1;
        state.cells[i].ud = 0;
        state.cells[i].arranjo = 0;
        if (fabs(uls) <= 1e-15 && state.cells[i].alf < (*state.globals).localtiny && uls1 < 0) {
            state.cells[i].term1 = 1.;
            state.cells[i].term2 = 0.;
            bif[i] = 0;
        } else if (fabs(uls1) < (*state.globals).localtiny * 1e-5) {
            if (fabs(uls) < (*state.globals).localtiny * 1e-5 && state.cells[i].alf < (*state.globals).localtiny && state.cells[i].fontemassGR >= (*state.globals).localtiny * 1e-5) {
                state.cells[i].term1 = 1.;
                state.cells[i].term2 = 0.;
                bif[i] = 0;
            }
        } else if (fabs(uls) < 1e-15 && uls1 < 0 && state.cells[i].alf <= (1 - 1 * (*state.globals).localtiny))
            bif[i] = 1;
    }

    else if (uls <= 0 && state.cells[i].alf >= 1. - (*state.globals).localtiny && ((state.cells[i].fontemassLL + state.cells[i].fontemassCL) <= (*state.globals).localtiny * 1e-5 && (state.cells[i].fontemassLR + state.cells[i].fontemassCR) <= (*state.globals).localtiny * 1e-5)) {
        state.cells[i].term1 = 0.;
        state.cells[i].term2 = 0.;
        bif[i] = 0;
        state.cells[i].c0 = 1;
        state.cells[i].ud = 0;
        state.cells[i].arranjo = 0;
        if (fabs(uls) <= (*state.globals).localtiny * 1e-5 && state.cells[i - 1].alf < (*state.globals).localtiny && uls0 > 0) {
            state.cells[i].term1 = 1.;
            state.cells[i].term2 = 0.;
            bif[i] = 0;
        }
        if (fabs(uls0) < (*state.globals).localtiny * 1e-5) {
            if (fabs(uls) < (*state.globals).localtiny * 1e-5 && state.cells[i - 1].alf < (*state.globals).localtiny && state.cells[i - 1].fontemassGR >= (*state.globals).localtiny * 1e-5) {
                state.cells[i].term1 = 1.;
                state.cells[i].term2 = 0.;
                bif[i] = 0;
            }
        } else if (fabs(uls) < (*state.globals).localtiny * 1e-5 && state.cells[i].duto.teta < 0.95 * M_PI / 2. && uls0 > 0 && ((state.cells[i - 1].alfPigD <= 1.0 * state.cells[i - 2].alfPigD && state.cells[i - 2].alfPigD < 0.99) || state.cells[i - 1].alfPigD < 0.7))
            bif[i] = 1;
        else if (fabs(uls) < (*state.globals).localtiny * 1e-5 && state.cells[i].duto.teta >= 0.95 * M_PI / 2. && uls0 > 0 && ((state.cells[i - 1].alfPigD <= 1.0 * state.cells[i - 2].alfPigD && state.cells[i - 2].alfPigD < 0.99) || state.cells[i - 1].alfPigD < 0.7))
            bif[i] = 1;
        else {
            if (fabs(uls) < (*state.globals).localtiny * 1e-5 && state.cells[i].duto.teta < 0.95 * M_PI / 2. && uls0 > 0 && ((state.cells[i - 1].alfPigD <= 1.0 * state.cells[i - 1].alfL) || state.cells[i - 1].alfPigD < 0.7))
                bif[i] = 1;
            else if (fabs(uls) < (*state.globals).localtiny * 1e-5 && state.cells[i].duto.teta >= 0.95 * M_PI / 2. && uls0 > 0 && ((state.cells[i - 1].alfPigD <= 1.0 * state.cells[i - 1].alfL) || state.cells[i - 1].alfPigD < 0.7))
                bif[i] = 1;
        }
    }

    if (bif[i] == 1) {
        double alfmed;
        alfmed = state.cells[i].alfL;
        double c0 = 1.2;
        double dmed = state.cells[i].duto.a;
        if (state.cells[i].MC >= 0)
            dmed = state.cells[i].dutoL.a;
        double sinal = 1.;
        if (state.cells[i].duto.teta < 0.)
            sinal = 1.;
        double ud = sinal * 0.32 * sqrt(9.82 * dmed);
        if (fabs(rgR) / rlR > 0.9) {
            c0 = 1.;
            ud = 0.;
        }
        state.closureUpdater.instantaneous(i, c0, ud);
        if (state.input.escorregamentoCelulaContorno == 0) {
            c0 = 1.;
            ud = 0.;
        }
        state.cells[i].c0 = c0;
        state.cells[i].ud = ud;
        double num = (1. - alfmed * c0);
        double den = 1. + alfmed * (rg / rl) * c0 - alfmed * c0;
        state.cells[i].term1 = num / den;
        state.cells[i].term2 = (-amed * alfmed * rg * ud) / den;

        // teste!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        // teste!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    }

    state.cells[i].term1L = state.cells[i - 1].term1;
    state.cells[i].term2L = state.cells[i - 1].term2;
    state.cells[i - 1].term1R = state.cells[i].term1;
    state.cells[i - 1].term2R = state.cells[i].term2;
}

void finalizeFlowPartitionTerms(
    const ThermalState &state, const Vcr<int> &bif, const Vcr<int> &valv) {
    int parada = 0;
    for (int i = 1; i < state.lastCell; i++) {
        if (state.cells[i].acsr.tipo == 5 && state.cells[i].acsr.chk.AreaGarg <= 1e-15 * state.cells[i].acsr.chk.AreaTub)
            parada = 1;
        else if (state.surfaceChoke.AreaGarg <= 1.e-15 * state.surfaceChoke.AreaTub)
            parada = 1;
    }

    Vcr<double> c0V(state.lastCell, 0.);
    Vcr<double> udV(state.lastCell, 0.);
    for (int i = 1; i < state.lastCell; i++) {
        c0V[i] = state.cells[i].c0;
        udV[i] = state.cells[i].ud;
        if (bif[i] == 1 && valv[i] == 1 && i > 2) {
            int iViz = i - 1;
            int iViz2 = i - 2;
            if (state.cells[i].QG < 0) {
                iViz = i + 1;
                iViz2 = i;
            }
            if ((state.cells[iViz].acsr.tipo == 0 && (state.cells[iViz2].acsr.tipo != 5 || state.cells[iViz2].acsr.chk.AreaGarg > (1e-3))) &&
                (state.cells[i].arranjo != state.cells[iViz].arranjo && bif[iViz] != 0)) {
                c0V[i] = (state.cells[i].dx * state.cells[i].c0 + state.cells[iViz].dx * state.cells[iViz].c0) / (state.cells[i].dx + state.cells[iViz].dx);
                if (state.cells[i].duto.teta * state.cells[iViz].duto.teta >= 0)
                    udV[i] = (state.cells[i].dx * state.cells[i].ud + state.cells[iViz].dx * state.cells[iViz].ud) / (state.cells[i].dx + state.cells[iViz].dx);
            } else if (i > 2) {
                double ang;
                double angL;
                double dia;
                double diaL;
                if ((state.cells[i - 1].acsr.tipo == 0 && (state.cells[i - 2].acsr.tipo != 5 || state.cells[i - 2].acsr.chk.AreaGarg > (1e-3))) &&
                    state.cells[i].QG >= 0) {
                    ang = state.cells[i].duto.teta;
                    angL = state.cells[i - 1].duto.teta;
                    dia = state.cells[i].duto.dia;
                    diaL = state.cells[i - 1].duto.dia;
                    if (((ang != angL) && bif[iViz] != 0) ||
                        (state.cells[i].arranjo != state.cells[iViz].arranjo && bif[iViz] != 0)) {
                        c0V[i] = (state.cells[i].dx * state.cells[i].c0 + state.cells[i - 1].dx * state.cells[i - 1].c0) / (state.cells[i].dx + state.cells[i - 1].dx);
                        if (ang * angL >= 0)
                            udV[i] = (state.cells[i].dx * state.cells[i].ud + state.cells[i - 1].dx * state.cells[i - 1].ud) / (state.cells[i].dx + state.cells[i - 1].dx);
                    }
                } else if ((state.cells[i + 1].acsr.tipo == 0 && (state.cells[i].acsr.tipo != 5 || state.cells[i].acsr.chk.AreaGarg > (1e-3))) &&
                           state.cells[i].QG < 0) {
                    ang = state.cells[i].duto.teta;
                    angL = state.cells[i + 1].duto.teta;
                    dia = state.cells[i].duto.dia;
                    diaL = state.cells[i + 1].duto.dia;
                    if (((ang != angL) && bif[iViz] != 0) ||
                        (state.cells[i].arranjo != state.cells[iViz].arranjo && bif[iViz] != 0)) {
                        c0V[i] = (state.cells[i].dx * state.cells[i].c0 + state.cells[i + 1].dx * state.cells[i + 1].c0) / (state.cells[i].dx + state.cells[i + 1].dx);
                        if (ang * angL >= 0)
                            udV[i] = (state.cells[i].dx * state.cells[i].ud + state.cells[i + 1].dx * state.cells[i + 1].ud) / (state.cells[i].dx + state.cells[i + 1].dx);
                    }
                }
            }
        }
    }

    for (int i = 1; i < state.lastCell; i++) {
        if (bif[i] == 1 && valv[i] == 1) {
            double betI = state.cells[i - 1].betPigD;
            double rl;

            if (state.cells[i].QL < 0.) {
                betI = state.cells[i].betPigE;
                rl = (1 - betI) * state.cells[i].rpCi + betI * state.cells[i].rcCi;
            } else {
                betI = state.cells[i - 1].betPigD;
                rl = (1 - betI) * state.cells[i - 1].rpCi + betI * state.cells[i - 1].rcCi;
            }

            double rg;
            double amed;
            double hns;
            if (state.cells[i].QG >= 0) {
                amed = state.cells[i].dutoL.area;
                hns = 1. - state.cells[i].alfL;
                rg = state.cells[i - 1].rgCi;
            } else {
                rg = state.cells[i].rgCi;
                amed = state.cells[i].duto.area;
                hns = 1. - state.cells[i].alf;
            }
            double ugs = state.cells[i].QG / (amed);
            double uls = state.cells[i].QL / (amed);

            double alfmed;
            alfmed = state.cells[i - 1].alfPigD;
            if (ugs < 0)
                alfmed = state.cells[i].alfPigE;
            double num = (1. - alfmed * state.cells[i].c0);
            double den = 1 + state.cells[i].c0 * alfmed * ((rg / rl) - 1.);
            state.cells[i].term1 = num / den;
            state.cells[i].term2 = (-amed * alfmed * rg * state.cells[i].ud) / den;
            double jlTeste0 = (ugs - alfmed * state.cells[i].ud) / (alfmed * state.cells[i].c0) - ugs;
            double jlTeste = (ugs + uls) * (1. - state.cells[i].c0 * alfmed) - alfmed * state.cells[i].ud;
            double MLTeste = state.cells[i].term1 * state.cells[i].MC + state.cells[i].term2;
            double MGTeste = (1 - state.cells[i].term1) * state.cells[i].MC - state.cells[i].term2;
            if ((jlTeste > 0. || jlTeste0 > 0.) && state.cells[i - 1].alfPigD > 1 - 1e-15) {
                state.cells[i].term1 = 0.;
                state.cells[i].term2 = 0.;
            }
            if ((jlTeste < 0. || jlTeste0 < 0.) && state.cells[i].alfPigE > 1 - 1e-15) {
                state.cells[i].term1 = 0.;
                state.cells[i].term2 = 0.;
            }
            if (state.cells[i - 1].acsr.tipo == 5 && state.cells[i - 1].acsr.chk.AreaGarg <= (1e-3 + state.input.master1.razareaativ) * state.cells[i - 1].duto.area) {
                state.cells[i].term1 = 0.;
                state.cells[i].term2 = 0.;
            }
            if ((jlTeste < 0. || jlTeste0 < 0.) && ((fabs(state.cells[i - 1].QG / (amed)) + fabs(state.cells[i - 1].QL / (amed))) < 0.1) &&
                (parada == 1 && state.input.modoSegrega == 1) && (state.cells[i].duto.teta > 0 && state.cells[i - 1].duto.teta <= 0) &&
                (MLTeste > 0 && state.cells[i].term2 > 0) &&
                ((1. - state.cells[i].alfL) /**fabs(sin(state.cells[i-1].duto.teta))*/ < (1. - state.cells[i].alf) /**fabs(sin(state.cells[i].duto.teta))*/)) {
                state.cells[i].term2 = 0.;
            }
            if (state.cells[i].duto.teta > 0 && MGTeste < 0 && ((parada == 1 && state.input.modoSegrega == 1)) && state.cells[i].alfPigE > 1 - 1e-15) {
                state.cells[i].term1 = 0.;
                state.cells[i].term2 = 0.;
            }
        }
    }
    for (int i = 1; i <= state.lastCell; i++) {
        state.cells[i].term1L = state.cells[i - 1].term1;
        state.cells[i].term2L = state.cells[i - 1].term2;
        state.cells[i - 1].term1R = state.cells[i].term1;
        state.cells[i - 1].term2R = state.cells[i].term2;
    }
}

void selectAndApplyInletBoundaryFlowRegime(
    const ThermalState &state, int i, Vcr<int> &bif, double xc0,
    double xud, double ugs, double uls, double uls1, double rgR,
    double rlR, double rg, double rl, double amed) {
    bif[i] = 1;

    if ((*state.globals).lixo5 > 29900) {
        int para;
        para = 0;
    }

    if (state.inletVoidFraction < (*state.globals).localtiny && state.cells[i].alfPigE <= (*state.globals).localtiny && state.cells[i].fontemassGL <= (*state.globals).localtiny * 1e-5 && state.cells[i].fontemassGR <= (*state.globals).localtiny * 1e-5) {
        state.cells[i].term1 = 1.;
        state.cells[i].term2 = 0.;
        bif[i] = 0;
        state.cells[i].c0 = 1 + 0 * xc0;
        state.cells[i].ud = 0 * xud;
        state.cells[i].arranjo = 0;
    } else if (state.inletVoidFraction >= (1. - (*state.globals).localtiny) && state.cells[i].alfPigE >= (1. - (*state.globals).localtiny) && ((state.cells[i].fontemassLL + state.cells[i].fontemassCL) <= (*state.globals).localtiny * 1e-5 && (state.cells[i].fontemassLR + state.cells[i].fontemassCR) <= (*state.globals).localtiny * 1e-5)) {
        state.cells[i].term1 = 0.;
        state.cells[i].term2 = 0.;
        bif[i] = 0;
        state.cells[i].c0 = 1 + 0 * xc0;
        state.cells[i].ud = 0 * xud;
        state.cells[i].arranjo = 0;
    } else if (ugs >= 0 && state.inletVoidFraction <= (*state.globals).localtiny && (state.cells[i].fontemassGL <= (*state.globals).localtiny * 1e-5 && state.cells[i].fontemassGR <= (*state.globals).localtiny * 1e-5)) {
        state.cells[i].term1 = 1.;
        state.cells[i].term2 = 0.;
        bif[i] = 0;
        state.cells[i].c0 = 1 + 0 * xc0;
        state.cells[i].ud = 0 * xud;
        state.cells[i].arranjo = 0;
        if (fabs(ugs) <= 1e-15 && state.cells[i].alfPigE > (1. - (*state.globals).localtiny) && uls < 0 && state.cells[i].duto.teta > 0) {
            state.cells[i].term1 = 0.;
            state.cells[i].term2 = 0.;
            bif[i] = 0;
        }
        if (fabs(ugs) <= 1e-15 && state.cells[i].alfPigE > (1. - 1 * (*state.globals).localtiny) && state.cells[i].duto.teta < 0 && uls > 0) {
            state.cells[i].term1 = 0.;
            state.cells[i].term2 = 0.;
            bif[i] = 1;
        }
        if (fabs(ugs) <= 1e-15 && state.cells[i].alf >= state.inletVoidFraction && uls < 0) {
            bif[i] = 1;
        }
    } else if (ugs <= 0 && state.cells[i].alfPigE <= (*state.globals).localtiny && (state.cells[i].fontemassGL <= (*state.globals).localtiny * 1e-5 && state.cells[i].fontemassGR <= (*state.globals).localtiny * 1e-5)) {
        state.cells[i].term1 = 1.;
        state.cells[i].term2 = 0.;
        bif[i] = 0;
        state.cells[i].c0 = 1 + 0 * xc0;
        state.cells[i].ud = 0 * xud;
        state.cells[i].arranjo = 0;
        if (fabs(ugs) <= 1e-15 && state.inletVoidFraction > (1. - (*state.globals).localtiny) && uls > 0) {
            state.cells[i].term1 = 0.;
            state.cells[i].term2 = 0.;
            bif[i] = 0;
        }
        if (fabs(ugs) <= 1e-15 && state.inletVoidFraction > (*state.globals).localtiny && uls > 0) {
            bif[i] = 1;
        }
    } else if (uls >= 0 && state.inletVoidFraction >= 1. - (*state.globals).localtiny && ((state.cells[i].fontemassLL + state.cells[i].fontemassCL) <= (*state.globals).localtiny * 1e-5 && (state.cells[i].fontemassLR + state.cells[i].fontemassCR) <= (*state.globals).localtiny * 1e-5)) {
        state.cells[i].term1 = 0.;
        state.cells[i].term2 = 0.;
        bif[i] = 0;
        state.cells[i].c0 = 1 + 0 * xc0;
        state.cells[i].ud = 0 * xud;
        state.cells[i].arranjo = 0;
        if (fabs(uls) <= 1e-15 && state.cells[i].alfPigE < (*state.globals).localtiny && uls1 < 0) {
            state.cells[i].term1 = 1.;
            state.cells[i].term2 = 0.;
            bif[i] = 0;
        } else if (fabs(uls1) < (*state.globals).localtiny * 1e-5) {
            if (fabs(uls) < (*state.globals).localtiny * 1e-5 && state.cells[i].alfPigE < (*state.globals).localtiny && state.cells[i].fontemassGR >= (*state.globals).localtiny * 1e-5) {
                state.cells[i].term1 = 1.;
                state.cells[i].term2 = 0.;
                bif[i] = 0;
            }
        } else if (fabs(uls) < 1e-15 && uls1 < 0 && ((state.cells[i].alfPigE <= (1 - 1 * (*state.globals).localtiny + .0 * state.cells[i].alfPigER) && state.cells[i].alfPigER < 1 - 1 * (*state.globals).localtiny) || state.cells[i].alfPigE <= 0.7))
            bif[i] = 1;

    } else if (uls <= 0 && state.cells[i].alfPigE >= 1. - (*state.globals).localtiny && ((state.cells[i].fontemassLL + state.cells[i].fontemassCL) <= (*state.globals).localtiny * 1e-5 && (state.cells[i].fontemassLR + state.cells[i].fontemassCR) <= (*state.globals).localtiny * 1e-5)) {
        state.cells[i].term1 = 0.;
        state.cells[i].term2 = 0.;
        bif[i] = 0;
        state.cells[i].c0 = 1 + 0 * xc0;
        state.cells[i].ud = 0 * xud;
        state.cells[i].arranjo = 0;
        if (fabs(uls) <= (*state.globals).localtiny * 1e-5 && state.inletVoidFraction < (*state.globals).localtiny && uls1 < 0) {
            state.cells[i].term1 = 1.;
            state.cells[i].term2 = 0.;
            bif[i] = 0;
        }

        if (fabs(uls1) < (*state.globals).localtiny * 1e-5) {
            if (fabs(uls) < (*state.globals).localtiny * 1e-5 && state.inletVoidFraction < (*state.globals).localtiny) {
                state.cells[i].term1 = 1.;
                state.cells[i].term2 = 0.;
                bif[i] = 0;
            }
        } else {
            if (fabs(uls) < (*state.globals).localtiny * 1e-5 && state.cells[i].duto.teta < 0.95 * M_PI / 2. && ugs < 0 && (state.inletVoidFraction < 0.7))
                bif[i] = 1;
            else if (fabs(uls) < (*state.globals).localtiny * 1e-5 && state.cells[i].duto.teta >= 0.95 * M_PI / 2. && ugs < 0 && (state.inletVoidFraction < 0.7))
                bif[i] = 1;
        }
    }
    if (uls > 0 && fabs(ugs) <= (*state.globals).localtiny * 1e-5 && state.inletVoidFraction > (*state.globals).localtiny && state.inletVoidFraction < 1 - (*state.globals).localtiny)
        bif[i] = 1;
    if (ugs > 0 && fabs(uls) <= (*state.globals).localtiny * 1e-5 && state.inletVoidFraction > (*state.globals).localtiny && state.inletVoidFraction < 1 - (*state.globals).localtiny)
        bif[i] = 1;
    if (uls >= 0 && state.inletVoidFraction > 1 - (*state.globals).localtiny) {
        state.cells[i].term1 = 0.;
        state.cells[i].term2 = 0.;
        bif[i] = 0;
    }
    if (ugs >= 0 && state.inletVoidFraction < (*state.globals).localtiny) {
        state.cells[i].term1 = 1.;
        state.cells[i].term2 = 0.;
        bif[i] = 0;
    }

    if (bif[i] == 1) {

        double c0;
        double ud;
        double alfmed;

        alfmed = state.inletVoidFraction;
        if (ugs < 0)
            alfmed = state.cells[i].alfPigE;
        c0 = 1.2;
        double dmed = state.cells[i].duto.a;
        if (state.cells[i].MC >= 0)
            dmed = state.cells[i].dutoL.a;
        double sinal = 1.;
        if (state.cells[i].duto.teta < 0.)
            sinal = -1.;
        ud = sinal * 0.32 * sqrt(9.82 * dmed);
        if (fabs(rgR) / rlR > 0.9) {
            c0 = 1.;
            ud = 0.;
        }
        state.closureUpdater.initialization(i, c0, ud);
        state.cells[i].c0 = c0;
        state.cells[i].ud = ud;
        double num = (1. - alfmed * c0);
        double den = 1 + c0 * alfmed * ((rg / rl) - 1.);
        state.cells[i].term1 = num / den;
        state.cells[i].term2 = (-amed * alfmed * rg * ud) / den;
    }
}

void selectAndApplyBufferedOutletFlowRegime(
    const ThermalState &state, int i, double xc0, double xud,
    double ugs, double uls, double ugs0, double uls0, double ugs1,
    double uls1, double rgR, double rlR, double rg, double rl,
    double amed) {
    int bif = 1;

    if ((*state.globals).lixo5 > 29900) {
        int para;
        para = 0;
    }

    if (state.cells[i - 1].alfPigD <= (*state.globals).localtiny && state.cells[i].alfPigE <= (*state.globals).localtiny && state.cells[i].fontemassGL <= (*state.globals).localtiny * 1e-5 && state.cells[i].fontemassGR <= (*state.globals).localtiny * 1e-5) {
        state.cells[i].term1 = 1.;
        state.cells[i].term2 = 0.;
        bif = 0;
        state.cells[i].c0 = 1 + 0 * xc0;
        state.cells[i].ud = 0 * xud;
        state.cells[i].arranjo = 0;
    } else if (state.cells[i - 1].alfPigD >= (1. - (*state.globals).localtiny) && state.cells[i].alfPigE >= (1. - (*state.globals).localtiny) && ((state.cells[i].fontemassLL + state.cells[i].fontemassCL) <= (*state.globals).localtiny * 1e-5 && (state.cells[i].fontemassLR + state.cells[i].fontemassCR) <= (*state.globals).localtiny * 1e-5)) {
        state.cells[i].term1 = 0.;
        state.cells[i].term2 = 0.;
        bif = 0;
        state.cells[i].c0 = 1 + 0 * xc0;
        state.cells[i].ud = 0 * xud;
        state.cells[i].arranjo = 0;
    } else if (state.cells[i - 1].acsr.tipo == 5 && state.cells[i - 1].acsr.chk.AreaGarg <= (1e-3)) {
        state.cells[i].term1 = 0.;
        state.cells[i].term2 = 0.;
        bif = 0;
        state.cells[i].c0 = 1 + 0 * xc0;
        state.cells[i].ud = 0 * xud;
        state.cells[i].arranjo = 0;
    }

    else if (ugs >= 0 && state.cells[i - 1].alfPigD <= (*state.globals).localtiny && (state.cells[i].fontemassGL <= (*state.globals).localtiny * 1e-5 && state.cells[i].fontemassGR <= (*state.globals).localtiny * 1e-5)) {
        state.cells[i].term1 = 1.;
        state.cells[i].term2 = 0.;
        bif = 0;
        state.cells[i].c0 = 1 + 0 * xc0;
        state.cells[i].ud = 0 * xud;
        state.cells[i].arranjo = 0;
        if (fabs(ugs) <= 1e-15 && state.cells[i].alfPigE > (1. - (*state.globals).localtiny) && uls < 0 && uls0 < 0 && state.cells[i].duto.teta > 0) {
            state.cells[i].term1 = 0.;
            state.cells[i].term2 = 0.;
            bif = 0;
        }
        if (fabs(ugs) <= 1e-15 && state.cells[i].alfPigE > (1. - 10 * (*state.globals).localtiny) && state.cells[i].duto.teta < 0) {
            state.cells[i].term1 = 0.;
            state.cells[i].term2 = 0.;
            bif = 1;
        }
        if (fabs(ugs) <= 1e-15 && state.cells[i].alf >= state.cells[i + 1].alf && ugs1 < 0) {
            bif = 1;
        }
    } else if (ugs <= 0 && state.cells[i].alfPigE <= (*state.globals).localtiny && (state.cells[i].fontemassGL <= (*state.globals).localtiny * 1e-5 && state.cells[i].fontemassGR <= (*state.globals).localtiny * 1e-5)) {
        state.cells[i].term1 = 1.;
        state.cells[i].term2 = 0.;
        bif = 0;
        state.cells[i].c0 = 1 + 0 * xc0;
        state.cells[i].ud = 0 * xud;
        state.cells[i].arranjo = 0;
        if (fabs(ugs) <= 1e-15 && state.cells[i - 1].alfPigD > (1. - (*state.globals).localtiny) && uls > 0 && uls1 > 0) {
            state.cells[i].term1 = 0.;
            state.cells[i].term2 = 0.;
            bif = 0;
        }
        if (fabs(ugs) <= 1e-15 && state.cells[i].alfL >= state.cells[i - 1].alfL && ugs0 > 0) {
            bif = 1;
        }
    } else if (uls >= 0 && state.cells[i - 1].alfPigD >= 1. - (*state.globals).localtiny && ((state.cells[i].fontemassLL + state.cells[i].fontemassCL) <= (*state.globals).localtiny * 1e-5 && (state.cells[i].fontemassLR + state.cells[i].fontemassCR) <= (*state.globals).localtiny * 1e-5)) {
        state.cells[i].term1 = 0.;
        state.cells[i].term2 = 0.;
        bif = 0;
        state.cells[i].c0 = 1 + 0 * xc0;
        state.cells[i].ud = 0 * xud;
        state.cells[i].arranjo = 0;
        if (fabs(uls) <= 1e-15 && state.cells[i].alfPigE < (*state.globals).localtiny && uls1 < 0) {
            state.cells[i].term1 = 1.;
            state.cells[i].term2 = 0.;
            bif = 0;
        } else if (fabs(uls1) < (*state.globals).localtiny * 1e-5) {
            if (fabs(uls) < (*state.globals).localtiny * 1e-5 && state.cells[i].alfPigE < (*state.globals).localtiny && state.cells[i].fontemassGR >= (*state.globals).localtiny * 1e-5) {
                state.cells[i].term1 = 1.;
                state.cells[i].term2 = 0.;
                bif = 0;
            }
        } else if (fabs(uls) < 1e-15 && uls1 < 0 && ((state.cells[i].alfPigE <= (1 - 1 * (*state.globals).localtiny + .0 * state.cells[i].alfPigER) && state.cells[i].alfPigER < 1 - 1 * (*state.globals).localtiny) || state.cells[i].alfPigE <= 0.7))
            bif = 1;

    } else if (uls <= 0 && state.cells[i].alfPigE >= 1. - (*state.globals).localtiny && ((state.cells[i].fontemassLL + state.cells[i].fontemassCL) <= (*state.globals).localtiny * 1e-5 && (state.cells[i].fontemassLR + state.cells[i].fontemassCR) <= (*state.globals).localtiny * 1e-5)) {
        state.cells[i].term1 = 0.;
        state.cells[i].term2 = 0.;
        bif = 0;
        state.cells[i].c0 = 1 + 0 * xc0;
        state.cells[i].ud = 0 * xud;
        state.cells[i].arranjo = 0;
        if (fabs(uls) <= (*state.globals).localtiny * 1e-5 && state.cells[i - 1].alfPigD < (*state.globals).localtiny && uls0 > 0) {
            state.cells[i].term1 = 1.;
            state.cells[i].term2 = 0.;
            bif = 0;
        }

        if (fabs(uls0) < (*state.globals).localtiny * 1e-5) {
            if (fabs(uls) < (*state.globals).localtiny * 1e-5 && state.cells[i - 1].alfPigD < (*state.globals).localtiny && state.cells[i - 1].fontemassGR >= (*state.globals).localtiny * 1e-5) {
                state.cells[i].term1 = 1.;
                state.cells[i].term2 = 0.;
                bif = 0;
            }
        } else if (i > 1 && fabs(uls) < (*state.globals).localtiny * 1e-5 && state.cells[i].duto.teta < 0.95 * M_PI / 2. && uls0 > 0 && ((state.cells[i - 1].alfPigD <= 1.0 * state.cells[i - 2].alfPigD && state.cells[i - 2].alfPigD < 0.99) || state.cells[i - 1].alfPigD < 0.7))
            bif = 1;
        else if (i > 1 && fabs(uls) < (*state.globals).localtiny * 1e-5 && state.cells[i].duto.teta >= 0.95 * M_PI / 2. && uls0 > 0 && ((state.cells[i - 1].alfPigD <= 1.0 * state.cells[i - 2].alfPigD && state.cells[i - 2].alfPigD < 0.99) || state.cells[i - 1].alfPigD < 0.7))
            bif = 1;
        else {
            if (fabs(uls) < (*state.globals).localtiny * 1e-5 && state.cells[i].duto.teta < 0.95 * M_PI / 2. && uls0 > 0 && ((state.cells[i - 1].alfPigD <= 1.0 * state.cells[i - 1].alfL) || state.cells[i - 1].alfPigD < 0.7))
                bif = 1;
            else if (fabs(uls) < (*state.globals).localtiny * 1e-5 && state.cells[i].duto.teta >= 0.95 * M_PI / 2. && uls0 > 0 && ((state.cells[i - 1].alfPigD <= 1.0 * state.cells[i - 1].alfL) || state.cells[i - 1].alfPigD < 0.7))
                bif = 1;
        }
    }

    if (ugs < 0 && state.cells[i].alf > (*state.globals).localtiny && state.cells[i].alf < 1. - (*state.globals).localtiny)
        bif = 1;

    if (bif == 1) {

        double c0;
        double ud;
        double alfmed;

        alfmed = state.cells[i - 1].alfPigD;
        if (ugs < 0)
            alfmed = state.cells[i].alfPigE;
        c0 = 1.2;
        double dmed = state.cells[i].duto.a;
        if (state.cells[i].MC >= 0)
            dmed = state.cells[i].dutoL.a;
        double sinal = 1.;
        if (state.cells[i].duto.teta < 0.)
            sinal = -1.;
        ud = sinal * 0.32 * sqrt(9.82 * dmed);
        if (fabs(rgR) / rlR > 0.9) {
            c0 = 1.;
            ud = 0.;
        }
        state.closureUpdater.buffered(i, c0, ud);
        if (state.input.escorregamentoCelulaContorno == 0) {
            c0 = 1.;
            ud = 0.;
        }
        double num = (1. - alfmed * c0);
        double den = 1 + c0 * alfmed * ((rg / rl) - 1.);
        state.cells[i].term1 = num / den;
        state.cells[i].term2 = (-amed * alfmed * rg * ud) / den;
        if (state.cells[i - 1].acsr.tipo == 5 && state.cells[i - 1].acsr.chk.AreaGarg <= (1e-3 + state.input.master1.razareaativ) * state.cells[i - 1].duto.area) {
            state.cells[i].term1 = 0.;
            state.cells[i].term2 = 0.;
        }
    }
}

void selectAndApplyBufferedInletFlowRegime(
    const ThermalState &state, int i, double xc0, double xud,
    double ugs, double uls, double uls1, double rgR, double rlR,
    double rg, double rl, double amed) {
    int bif = 1;

    if ((*state.globals).lixo5 > 29900) {
        int para;
        para = 0;
    }

    if (state.inletVoidFraction < (*state.globals).localtiny && state.cells[i].alfPigE <= (*state.globals).localtiny && state.cells[i].fontemassGL <= (*state.globals).localtiny * 1e-5 && state.cells[i].fontemassGR <= (*state.globals).localtiny * 1e-5) {
        state.cells[i].term1 = 1.;
        state.cells[i].term2 = 0.;
        bif = 0;
        state.cells[i].c0 = 1 + 0 * xc0;
        state.cells[i].ud = 0 * xud;
        state.cells[i].arranjo = 0;
    } else if (state.inletVoidFraction >= (1. - (*state.globals).localtiny) && state.cells[i].alfPigE >= (1. - (*state.globals).localtiny) && ((state.cells[i].fontemassLL + state.cells[i].fontemassCL) <= (*state.globals).localtiny * 1e-5 && (state.cells[i].fontemassLR + state.cells[i].fontemassCR) <= (*state.globals).localtiny * 1e-5)) {
        state.cells[i].term1 = 0.;
        state.cells[i].term2 = 0.;
        bif = 0;
        state.cells[i].c0 = 1 + 0 * xc0;
        state.cells[i].ud = 0 * xud;
        state.cells[i].arranjo = 0;
    } else if (ugs >= 0 && state.inletVoidFraction <= (*state.globals).localtiny && (state.cells[i].fontemassGL <= (*state.globals).localtiny * 1e-5 && state.cells[i].fontemassGR <= (*state.globals).localtiny * 1e-5)) {
        state.cells[i].term1 = 1.;
        state.cells[i].term2 = 0.;
        bif = 0;
        state.cells[i].c0 = 1 + 0 * xc0;
        state.cells[i].ud = 0 * xud;
        state.cells[i].arranjo = 0;
        if (fabs(ugs) <= 1e-15 && state.cells[i].alfPigE > (1. - (*state.globals).localtiny) && uls < 0 && state.cells[i].duto.teta > 0) {
            state.cells[i].term1 = 0.;
            state.cells[i].term2 = 0.;
            bif = 0;
        }
        if (fabs(ugs) <= 1e-15 && state.cells[i].alfPigE > (1. - 1 * (*state.globals).localtiny) && state.cells[i].duto.teta < 0 && uls > 0) {
            state.cells[i].term1 = 0.;
            state.cells[i].term2 = 0.;
            bif = 1;
        }
        if (fabs(ugs) <= 1e-15 && state.cells[i].alf >= state.inletVoidFraction && uls < 0) {
            bif = 1;
        }
    } else if (ugs <= 0 && state.cells[i].alfPigE <= (*state.globals).localtiny && (state.cells[i].fontemassGL <= (*state.globals).localtiny * 1e-5 && state.cells[i].fontemassGR <= (*state.globals).localtiny * 1e-5)) {
        state.cells[i].term1 = 1.;
        state.cells[i].term2 = 0.;
        bif = 0;
        state.cells[i].c0 = 1 + 0 * xc0;
        state.cells[i].ud = 0 * xud;
        state.cells[i].arranjo = 0;
        if (fabs(ugs) <= 1e-15 && state.inletVoidFraction > (1. - (*state.globals).localtiny) && uls > 0) {
            state.cells[i].term1 = 0.;
            state.cells[i].term2 = 0.;
            bif = 0;
        }
        if (fabs(ugs) <= 1e-15 && state.inletVoidFraction > (*state.globals).localtiny && uls > 0) {
            bif = 1;
        }
    } else if (uls >= 0 && state.inletVoidFraction >= 1. - (*state.globals).localtiny && ((state.cells[i].fontemassLL + state.cells[i].fontemassCL) <= (*state.globals).localtiny * 1e-5 && (state.cells[i].fontemassLR + state.cells[i].fontemassCR) <= (*state.globals).localtiny * 1e-5)) {
        state.cells[i].term1 = 0.;
        state.cells[i].term2 = 0.;
        bif = 0;
        state.cells[i].c0 = 1 + 0 * xc0;
        state.cells[i].ud = 0 * xud;
        state.cells[i].arranjo = 0;
        if (fabs(uls) <= 1e-15 && state.cells[i].alfPigE < (*state.globals).localtiny && uls1 < 0) {
            state.cells[i].term1 = 1.;
            state.cells[i].term2 = 0.;
            bif = 0;
        } else if (fabs(uls1) < (*state.globals).localtiny * 1e-5) {
            if (fabs(uls) < (*state.globals).localtiny * 1e-5 && state.cells[i].alfPigE < (*state.globals).localtiny && state.cells[i].fontemassGR >= (*state.globals).localtiny * 1e-5) {
                state.cells[i].term1 = 1.;
                state.cells[i].term2 = 0.;
                bif = 0;
            }
        } else if (fabs(uls) < 1e-15 && uls1 < 0 && ((state.cells[i].alfPigE <= (1 - 1 * (*state.globals).localtiny + .0 * state.cells[i].alfPigER) && state.cells[i].alfPigER < 1 - 1 * (*state.globals).localtiny) || state.cells[i].alfPigE <= 0.7))
            bif = 1;

    } else if (uls <= 0 && state.cells[i].alfPigE >= 1. - (*state.globals).localtiny && ((state.cells[i].fontemassLL + state.cells[i].fontemassCL) <= (*state.globals).localtiny * 1e-5 && (state.cells[i].fontemassLR + state.cells[i].fontemassCR) <= (*state.globals).localtiny * 1e-5)) {
        state.cells[i].term1 = 0.;
        state.cells[i].term2 = 0.;
        bif = 0;
        state.cells[i].c0 = 1 + 0 * xc0;
        state.cells[i].ud = 0 * xud;
        state.cells[i].arranjo = 0;
        if (fabs(uls) <= (*state.globals).localtiny * 1e-5 && state.inletVoidFraction < (*state.globals).localtiny && uls1 < 0) {
            state.cells[i].term1 = 1.;
            state.cells[i].term2 = 0.;
            bif = 0;
        }

        if (fabs(uls1) < (*state.globals).localtiny * 1e-5) {
            if (fabs(uls) < (*state.globals).localtiny * 1e-5 && state.inletVoidFraction < (*state.globals).localtiny) {
                state.cells[i].term1 = 1.;
                state.cells[i].term2 = 0.;
                bif = 0;
            }
        } else {
            if (fabs(uls) < (*state.globals).localtiny * 1e-5 && state.cells[i].duto.teta < 0.95 * M_PI / 2. && ugs < 0 && (state.inletVoidFraction < 0.7))
                bif = 1;
            else if (fabs(uls) < (*state.globals).localtiny * 1e-5 && state.cells[i].duto.teta >= 0.95 * M_PI / 2. && ugs < 0 && (state.inletVoidFraction < 0.7))
                bif = 1;
        }
    }
    if (uls > 0 && fabs(ugs) <= (*state.globals).localtiny * 1e-5 && state.inletVoidFraction > (*state.globals).localtiny && state.inletVoidFraction < 1 - (*state.globals).localtiny)
        bif = 1;
    if (ugs > 0 && fabs(uls) <= (*state.globals).localtiny * 1e-5 && state.inletVoidFraction > (*state.globals).localtiny && state.inletVoidFraction < 1 - (*state.globals).localtiny)
        bif = 1;

    if (bif == 1) {

        double c0;
        double ud;
        double alfmed;

        alfmed = state.inletVoidFraction;
        if (ugs < 0)
            alfmed = state.cells[i].alfPigE;
        c0 = 1.2;
        double dmed = state.cells[i].duto.a;
        if (state.cells[i].MC >= 0)
            dmed = state.cells[i].dutoL.a;
        double sinal = 1.;
        if (state.cells[i].duto.teta < 0.)
            sinal = -1.;
        ud = sinal * 0.32 * sqrt(9.82 * dmed);
        if (fabs(rgR) / rlR > 0.9) {
            c0 = 1.;
            ud = 0.;
        }
        state.closureUpdater.bufferedInitialization(i, c0, ud);

        double num = (1. - alfmed * c0);
        double den = 1 + c0 * alfmed * ((rg / rl) - 1.);
        state.cells[i].term1 = num / den;
        state.cells[i].term2 = (-amed * alfmed * rg * ud) / den;
    }
}

}  // namespace

void updateFlowPartitionTerms(const ThermalState &state, int aflu) {
    // #pragma omp parallel for num_threads(numthreads)
    aflu = 0;
    Vcr<int> bif(state.lastCell + 1, 0);
    Vcr<int> valv(state.lastCell + 1, 1);
    int imax = state.lastCell;
    if (aflu == 1)
        imax = state.lastCell + 1;
    for (int i = 0; i <= state.lastCell; i++) {
        state.cells[i].c0ini = state.cells[i].c0;
        state.cells[i].udini = state.cells[i].ud;
        if (i == state.lastCell && aflu == 1 && (*state.globals).lixo5 >= 1560) {
            int para;
            para = 0;
        }
        if (i != 0 && i != imax) {
            updateInteriorFlowPartitionCell(state, i, bif, valv);
        } else if (i == 0) {
            if (state.input.ConContEntrada == 0) {
                state.cells[1].alfL = state.cells[0].alf;
                state.cells[0].alfL = state.cells[0].alf;
                state.cells[1].betL = state.cells[0].bet;
                state.cells[0].betL = state.cells[0].bet;
                state.cells[0].term1 = 0.;
                state.cells[0].term2 = 0.;
                state.cells[0].term1L = 0.;
                state.cells[0].term2L = 0.;
            } else {
                if (state.inletMassFraction < 1) {
                    int para;
                    para = 0;
                }

                double razdx = 0.5;
                double pmed;
                pmed = state.inletPressure;

                double tmed;
                if (state.cells[0].QL < 0.)
                    tmed = state.cells[i].temp;
                else
                    tmed = state.inletTemperature;

                double rg = state.cells[0].flui.MasEspGas(state.inletPressure, state.inletTemperature);
                double rl = state.cells[0].flui.MasEspLiq(state.inletPressure, state.inletTemperature);
                double rcis = state.cells[0].fluicol.MasEspFlu(state.inletPressure, state.inletTemperature);

                double rlmist = state.inletComposition * rcis + (1 - state.inletComposition) * rl;
                state.inletVoidFraction = (-state.inletMassFraction * rlmist / (state.inletMassFraction * rg - rg - state.inletMassFraction * rlmist)) / (state.cells[0].c0);

                double betI;
                double viscl1;
                double tensup1;
                if (state.cells[i].QL < 0.) { // testeBeta
                    betI = state.cells[i].betPigE;
                    rl = (1 - betI) * state.cells[i].flui.MasEspLiq(pmed, tmed) + betI * state.cells[i].fluicol.MasEspFlu(pmed, tmed);
                    viscl1 = (1 - betI) * state.cells[i].flui.ViscOleo(pmed, tmed) + betI * state.cells[i].fluicol.VisFlu(pmed, tmed);
                    tensup1 = (1 - betI) * state.cells[i].flui.TensSuper(pmed, tmed) + betI * state.cells[i].fluicol.TensSuper(pmed, tmed);
                } else {
                    betI = state.inletComposition;
                    rl = (1 - betI) * (*state.cells[i].fluiL).MasEspLiq(pmed, tmed) + betI * state.cells[i].fluicol.MasEspFlu(pmed, tmed);
                    viscl1 = (1 - betI) * (*state.cells[i].fluiL).ViscOleo(pmed, tmed) + betI * state.cells[i].fluicol.VisFlu(pmed, tmed);
                    tensup1 = (1 - betI) * (*state.cells[i].fluiL).TensSuper(pmed, tmed) + betI * state.cells[i].fluicol.TensSuper(pmed, tmed);
                }

                double viscg1;
                double amed;
                double hns;
                if (state.cells[i].QG >= 0) {
                    amed = state.cells[i].duto.area;
                    rg = (*state.cells[i].fluiL).MasEspGas(pmed, tmed);
                    viscg1 = (*state.cells[i].fluiL).ViscGas(pmed, tmed);
                    hns = 1. - state.inletVoidFraction;
                } else {
                    rg = state.cells[i].flui.MasEspGas(pmed, tmed);
                    viscg1 = state.cells[i].flui.ViscGas(pmed, tmed);
                    amed = state.cells[i].duto.area;
                    hns = 1. - state.cells[i].alf;
                }
                double ugs = state.cells[i].QG / (amed);
                double uls = state.cells[i].QL / (amed);
                double dia1 = state.cells[i].duto.a;

                double rmed = hns * rl + (1 - hns) * rg;
                double visc = (hns * viscl1 + (1 - hns) * viscg1) / pow(10., 3.);
                double ang = state.cells[i].duto.teta;
                double sinal = 1.;
                if (ang < 0.)
                    sinal = -1.;
                double xc0 = 2.;
                double xud = sinal * 0.0246 * cos(ang) + 1.606 * pow(9.82 * tensup1 * (rl - rg) / (rl * rl), 0.25) * sin(ang);

                double amedR = state.cells[i].dutoR.area;
                double razdxR = state.cells[i].dx / (state.cells[i].dxR + state.cells[i].dx);
                double pmedR = state.cells[i].presauxR;
                double tmedR = state.cells[i].temp * razdxR + state.cells[i].tempL * (1. - razdxR);
                double betIR = state.cells[i].betPigD;
                if (state.cells[i].QLR < 0.) // testeBeta
                    betIR = state.cells[i + 1].betPigE;

                double rgR = state.cells[i].flui.MasEspGas(pmedR, tmedR);
                double rlR = (1 - betIR) * state.cells[i].flui.MasEspLiq(pmedR, tmedR) + betIR * state.cells[i].fluicol.MasEspFlu(pmedR, tmedR);
                double ugs1 = (state.cells[i].MR - state.cells[i].MliqiniR) / (rgR * amedR);
                double uls1 = (state.cells[i].MliqiniR) / (rlR * amedR);

                selectAndApplyInletBoundaryFlowRegime(
                    state, i, bif, xc0, xud, ugs, uls, uls1,
                    rgR, rlR, rg, rl, amed);
                state.cells[1].term1L = state.cells[i].term1;
                state.cells[1].term2L = state.cells[i].term2;

                state.cells[1].alfL = state.cells[0].alf;
                state.cells[0].alfL = state.inletVoidFraction;
                state.cells[1].betL = state.cells[0].bet;
                state.cells[0].betL = state.inletComposition;
            }

        } else if (aflu == 0) {
            updateOutletBoundaryFlowPartition(state, i, bif);
        }
    }

    finalizeFlowPartitionTerms(state, bif, valv);
}
void updateOutletFlowPartitionTerms(const ThermalState &state) {

    int i = state.lastCell;

    double razdx = state.cells[i].dx / (state.cells[i].dx + state.cells[i].dxL);
    double pmed = state.cells[i].presBuf;
    double tmed = state.gasSurfaceTemperature;
    tmed = state.cells[i - 1].temp;
    if (state.cells[i].VTemper < 0.)
        tmed = state.gasSurfaceTemperature;
    double betI = state.cells[i - 1].betPigD;
    double rl;
    double viscl1;
    double tensup1;
    if (state.cells[i].MliqiniBuf < 0.) { // testeBeta
        betI = state.cells[i].betPigE;
        rl = (1 - betI) * state.cells[i].flui.MasEspLiq(pmed, tmed) + betI * state.cells[i].fluicol.MasEspFlu(pmed, tmed);
        viscl1 = (1 - betI) * state.cells[i].flui.ViscOleo(pmed, tmed) + betI * state.cells[i].fluicol.VisFlu(pmed, tmed);
        tensup1 = (1 - betI) * state.cells[i].flui.TensSuper(pmed, tmed) + betI * state.cells[i].fluicol.TensSuper(pmed, tmed);
    } else {
        betI = state.cells[i - 1].betPigD;
        rl = (1 - betI) * state.cells[i - 1].flui.MasEspLiq(pmed, tmed) + betI * state.cells[i - 1].fluicol.MasEspFlu(pmed, tmed);
        viscl1 = (1 - betI) * state.cells[i - 1].flui.ViscOleo(pmed, tmed) + betI * state.cells[i - 1].fluicol.VisFlu(pmed, tmed);
        tensup1 = (1 - betI) * state.cells[i - 1].flui.TensSuper(pmed, tmed) + betI * state.cells[i - 1].fluicol.TensSuper(pmed, tmed);
    }
    double rg;
    double viscg1;
    double amed;
    double hns;
    if ((state.cells[i].MCBuf - state.cells[i].MliqiniBuf) >= 0) {
        amed = state.cells[i].dutoL.area;
        hns = 1. - state.cells[i].alfL;
        rg = state.cells[i - 1].flui.MasEspGas(pmed, tmed);
        viscg1 = state.cells[i - 1].flui.ViscGas(pmed, tmed);
    } else {
        rg = state.cells[i].flui.MasEspGas(pmed, tmed);
        viscg1 = state.cells[i].flui.ViscGas(pmed, tmed);
        amed = state.cells[i].duto.area;
        hns = 1. - state.cells[i].alf;
    }
    double ugs = (state.cells[i].MCBuf - state.cells[i].MliqiniBuf) / (rg * amed);
    double uls = (state.cells[i].MliqiniBuf) / (rl * amed);
    double dia1 = state.cells[i].duto.a;
    if (ugs >= 0)
        dia1 = state.cells[i - 1].duto.a;

    double rmed = hns * rl + (1 - hns) * rg;
    double visc = (hns * viscl1 + (1 - hns) * viscg1) / pow(10., 3.);
    double nrey = dia1 * rmed * (fabs(ugs) / amed + fabs(uls) / amed) / visc;
    double ang = state.cells[i].duto.teta;
    if (i >= 2) {
        if (state.cells[i - 2].acsr.tipo == 5 && state.cells[i - 2].acsr.chk.AreaGarg <= (1e-3)) {
            if ((state.cells[i].MCBuf - state.cells[i].MliqiniBuf) >= 0)
                ang = state.cells[i].duto.teta;
            else
                ang = state.cells[i].dutoR.teta;
        } else {
            if ((state.cells[i].MCBuf - state.cells[i].MliqiniBuf) >= 0)
                ang = state.cells[i].dutoL.teta;
            else
                ang = state.cells[i].duto.teta;
        }
    }
    double sinal = 1.;
    if (ang < 0.)
        sinal = -1.;
    double xc0 = 2.;
    double xud = sinal * 0.0246 * cos(ang) + 1.606 * pow(9.82 * tensup1 * (rl - rg) / (rl * rl), 0.25) * sin(ang);

    double amedL = state.cells[i].dutoL.area;
    double razdxL = state.cells[i - 1].dxL / (state.cells[i - 1].dx + state.cells[i - 1].dxL);
    double pmedL = state.cells[i - 1].presaux;
    double tmedL = state.cells[i - 1].temp * razdxL + state.cells[i - 1].tempL * (1. - razdxL);
    double betIL = 0.;
    if (i < 2)
        betIL = state.cells[i - 1].betL;
    else
        betIL = state.cells[i - 2].betPigD;
    if (state.cells[i - 1].MliqiniBuf < 0.)
        betIL = state.cells[i - 1].betPigE; // testeBeta
    // betIL = state.cells[i - 1].betPigE;    //duvidabeta
    double rgL = state.cells[i].flui.MasEspGas(pmedL, tmedL);
    double rlL = (1 - betIL) * state.cells[i].flui.MasEspLiq(pmedL, tmedL) + betIL * state.cells[i].fluicol.MasEspFlu(pmedL, tmedL);
    double ugs0 = (state.cells[i].MLBuf - state.cells[i].MliqiniLBuf) / (rgL * amedL);
    double uls0 = (state.cells[i].MliqiniLBuf) / (rlL * amedL);

    double amedR = state.cells[i].dutoR.area;
    double razdxR = state.cells[i].dxR / (state.cells[i].dxR + state.cells[i].dx);
    double pmedR = state.cells[i].presRBuf;
    double tmedR = state.cells[i].temp * razdxR + state.cells[i].tempR * (1. - razdxR);
    double betIR = state.cells[i].betPigD;
    if (state.cells[i].MliqiniRBuf < 0.) { // testeBeta
        if (i > state.lastCell - 2)
            betIR = state.cells[i].betR;
        else
            betIR = state.cells[i + 1].betPigE;
    }
    double rgR = state.cells[i].flui.MasEspGas(pmedR, tmedR);
    double rlR = (1 - betIR) * state.cells[i].flui.MasEspLiq(pmedR, tmedR) + betIR * state.cells[i].fluicol.MasEspFlu(pmedR, tmedR);
    double ugs1 = (state.cells[i].MRBuf - state.cells[i].MliqiniRBuf) / (rgR * amedR);
    double uls1 = (state.cells[i].MliqiniRBuf) / (rlR * amedR);

    selectAndApplyBufferedOutletFlowRegime(
        state, i, xc0, xud, ugs, uls, ugs0, uls0, ugs1, uls1,
        rgR, rlR, rg, rl, amed);
    state.cells[i].term1L = state.cells[i - 1].term1;
    state.cells[i].term2L = state.cells[i - 1].term2;
    state.cells[i - 1].term1R = state.cells[i].term1;
    state.cells[i - 1].term2R = state.cells[i].term2;
}
void updateInletFlowPartitionTerms(const ThermalState &state) {

    if (state.inletMassFraction < 1) {
        int para;
        para = 0;
    }

    int i = 0;

    double razdx = 0.5;
    double pmed;
    pmed = state.inletPressure;

    double tmed;
    if (state.cells[0].MliqiniBuf < 0.)
        tmed = state.cells[i].temp;
    else
        tmed = state.inletTemperature;

    double betI;
    double viscl1;
    double tensup1;
    double rg = state.cells[0].flui.MasEspGas(state.inletPressure, state.inletTemperature);
    double rl = state.cells[0].flui.MasEspLiq(state.inletPressure, state.inletTemperature);
    double rcis = state.cells[0].fluicol.MasEspFlu(state.inletPressure, state.inletTemperature);

    double rlmist = state.inletComposition * rcis + (1 - state.inletComposition) * rl;
    state.inletVoidFraction = (-state.inletMassFraction * rlmist / (state.inletMassFraction * rg - rg - state.inletMassFraction * rlmist)) / (state.cells[0].c0);

    if ((state.cells[i].MCBuf - state.cells[0].MliqiniBuf) * 0 + 1 * state.cells[0].MliqiniBuf < 0.) { // duvidabeta
        betI = state.cells[i].betPigE;
        rl = (1 - betI) * state.cells[i].flui.MasEspLiq(pmed, tmed) + betI * state.cells[i].fluicol.MasEspFlu(pmed, tmed);
        viscl1 = (1 - betI) * state.cells[i].flui.ViscOleo(pmed, tmed) + betI * state.cells[i].fluicol.VisFlu(pmed, tmed);
        tensup1 = (1 - betI) * state.cells[i].flui.TensSuper(pmed, tmed) + betI * state.cells[i].fluicol.TensSuper(pmed, tmed);
    } else {
        betI = state.inletComposition;
        rl = (1 - betI) * (*state.cells[i].fluiL).MasEspLiq(pmed, tmed) + betI * state.cells[i].fluicol.MasEspFlu(pmed, tmed);
        viscl1 = (1 - betI) * (*state.cells[i].fluiL).ViscOleo(pmed, tmed) + betI * state.cells[i].fluicol.VisFlu(pmed, tmed);
        tensup1 = (1 - betI) * (*state.cells[i].fluiL).TensSuper(pmed, tmed) + betI * state.cells[i].fluicol.TensSuper(pmed, tmed);
    }
    double viscg1;
    double amed;
    double hns;
    if (state.cells[i].MCBuf - state.cells[0].MliqiniBuf >= 0) {
        amed = state.cells[i].dutoL.area;
        rg = (*state.cells[i].fluiL).MasEspGas(pmed, tmed);
        viscg1 = (*state.cells[i].fluiL).ViscGas(pmed, tmed);
        hns = 1. - state.inletVoidFraction;
    } else {
        rg = state.cells[i].flui.MasEspGas(pmed, tmed);
        viscg1 = state.cells[i].flui.ViscGas(pmed, tmed);
        amed = state.cells[i].duto.area;
        hns = 1. - state.cells[i].alf;
    }
    double ugs = (state.cells[i].MCBuf - state.cells[i].MliqiniBuf) / (rg * amed);
    double uls = state.cells[i].MliqiniBuf / (rl * amed);
    double dia1 = state.cells[i].duto.a;
    if (ugs >= 0)
        dia1 = state.cells[i].duto.a;

    double rmed = hns * rl + (1 - hns) * rg;
    double visc = (hns * viscl1 + (1 - hns) * viscg1) / pow(10., 3.);
    double ang = state.cells[i].duto.teta;
    double sinal = 1.;
    if (ang < 0.)
        sinal = -1.;
    double xc0 = 2.;
    double xud = sinal * 0.0246 * cos(ang) + 1.606 * pow(9.82 * tensup1 * (rl - rg) / (rl * rl), 0.25) * sin(ang);

    double amedR = state.cells[i].dutoR.area;
    double razdxR = state.cells[i].dxR / (state.cells[i].dxR + state.cells[i].dx);
    double pmedR = state.cells[i].presRBuf * razdxR + state.cells[i].presBuf * (1. - razdxR);
    double tmedR = state.cells[i].temp * razdxR + state.cells[i].tempR * (1. - razdxR);
    double betIR = state.cells[i].betPigD;
    if (state.cells[i].QLR < 0.) // testeBeta
        betIR = state.cells[i + 1].betPigE;

    double rgR = state.cells[i].flui.MasEspGas(pmedR, tmedR);
    double rlR = (1 - betIR) * state.cells[i].flui.MasEspLiq(pmedR, tmedR) + betIR * state.cells[i].fluicol.MasEspFlu(pmedR, tmedR);
    double ugs1 = (state.cells[i].MRBuf - state.cells[i].MliqiniRBuf) / (rgR * amedR);
    double uls1 = (state.cells[i].MliqiniRBuf) / (rlR * amedR);

    selectAndApplyBufferedInletFlowRegime(
        state, i, xc0, xud, ugs, uls, uls1, rgR, rlR, rg, rl,
        amed);
    state.cells[1].term1L = state.cells[i].term1;
    state.cells[1].term2L = state.cells[i].term2;
}
void prepareNonDimensionalHeatDiffusion(const ThermalState &state, int i) {
    double dia = state.cells[i].duto.a;
    double area = 0.25 * M_PI * dia * dia;
    double alfmed = state.cells[i].alf;
    double betmed = state.cells[i].bet;
    double ugsmed;
    double ulsmed;
    if (i > 0 && (state.cells[i - 1].acsr.tipo != 5 ||
                  state.cells[i - 1].acsr.chk.AreaGarg > (1e-3 + state.input.master1.razareaativ) * state.cells[i - 1].duto.area)) {
        if (state.cells[i].alf > (*state.globals).localtiny)
            ugsmed = state.cells[i].QG / area;
        else {
            ugsmed = 0.;
        }
        if (state.cells[i].alf < 1. - (*state.globals).localtiny)
            ulsmed = state.cells[i].QL / area;
        else {
            ulsmed = 0.;
        }
    } else {
        if (state.cells[i].alf > (*state.globals).localtiny)
            ugsmed = state.cells[i + 1].QG / area;
        else {
            ugsmed = 0.;
        }

        if (state.cells[i].alf < 1. - (*state.globals).localtiny)
            ulsmed = state.cells[i + 1].QL / area;
        else {
            ulsmed = 0.;
        }
    }
    double rp = state.cells[i].rpC;
    double rc = state.cells[i].rcC;
    double rhol = (1. - betmed) * rp + betmed * rc;
    double rhog = state.cells[i].rgC;
    double cpl = (1. - betmed) * state.cells[i].flui.CalorLiq(state.cells[i].presini, state.cells[i].temp) + betmed * state.cells[i].fluicol.CalorLiq(state.cells[i].presini, state.cells[i].temp);
    double cpg = state.cells[i].flui.CalorGas(state.cells[i].presini, state.cells[i].temp);

    state.cells[i].calor.Tint = state.cells[i].temp;
    state.cells[i].calor.dtL = state.cells[i].temp - state.cells[i - 1].tempini;
    state.cells[i].calor.Vint = ugsmed + ulsmed;
    state.cells[i].calor.dt = state.cells[i].dt;
    double condliq = (1. - betmed) * state.cells[i].flui.CondLiq(state.cells[i].presini, state.cells[i].temp) + betmed * state.cells[i].fluicol.CondLiq(state.cells[i].presini, state.cells[i].temp); //(1. - betmed) * celula[i].flui.CondLiq(celula[i].pres, celula[i].temp)
    state.cells[i].calor.kint = condliq * (1 - alfmed) + state.cells[i].flui.CondGas(state.cells[i].presini, state.cells[i].temp) * alfmed;                                                 // condliq * (1 - alfmed) + celula[i].flui.CondGas(celula[i].pres, celula[i].temp) * alfmed;
    state.cells[i].calor.cpint = cpl * (1 - alfmed) + cpg * alfmed;
    state.cells[i].calor.rhoint = rhol * (1 - alfmed) + rhog * alfmed;
    //(1. - betmed) * celula[i].flui.ViscOleo(celula[i].pres, celula[i].temp)
    double viscliq = (1. - betmed) * state.cells[i].mipC + betmed * state.cells[i].micC;
    // viscliq * (1 - alfmed) * 1.e-3
    state.cells[i].calor.viscint = viscliq * (1 - alfmed) * 1.e-3 + state.cells[i].migC * alfmed * 1.e-3;
    double dtemp = state.cells[i].temp * 0.01;
    if (fabs(state.cells[i].temp) < 1e-15)
        dtemp = 0.1;
    double rholdT = (1. - betmed) * state.cells[i].flui.MasEspLiq(state.cells[i].presini, state.cells[i].temp + dtemp) +
                    betmed * state.cells[i].fluicol.MasEspFlu(state.cells[i].presini, state.cells[i].temp + dtemp) - rhol; //(1. - betmed) * celula[i].flui.MasEspLiq(celula[i].pres, celula[i].temp+dtemp) +
    double rhogdT = state.cells[i].flui.MasEspGas(state.cells[i].presini, state.cells[i].temp + dtemp) - rhog;             // celula[i].flui.MasEspGas(celula[i].pres, celula[i].temp+dtemp)-rhog;
    state.cells[i].calor.betint = -(1 / state.cells[i].calor.rhoint) * (rholdT * (1 - alfmed) + rhogdT * alfmed) / (dtemp);
}

void advanceTransientEnergy(const ThermalState &state, int ciclo, int ciclomax) {
    if (((*state.globals).chaverede == 0 || state.networkEndpoint == 1 || (*state.globals).chaveRedeParalela == 1) && state.input.chkv == 0) {
        if (state.cells[state.lastCell - 1].MliqiniR < 0) {
            state.cells[state.lastCell - 1].MR = state.cells[state.lastCell - 1].MR - state.cells[state.lastCell - 1].MliqiniR;
            state.cells[state.lastCell - 1].MliqiniR = 0;
            state.cells[state.lastCell - 1].term1R = 0;
            state.cells[state.lastCell - 1].term2R = 0;
            state.cells[state.lastCell].MC = state.cells[state.lastCell].MC - state.cells[state.lastCell].Mliqini;
            state.cells[state.lastCell].Mliqini = 0;
            state.cells[state.lastCell].term1 = 0;
            state.cells[state.lastCell].term2 = 0;
        }
        if (state.cells[state.lastCell - 1].QLR < 0 && (state.surfaceChokeMassCondition == 0 || state.surfaceChokeOpen == 1)) {
            state.cells[state.lastCell - 1].QLR = 0;
            state.cells[state.lastCell].QL = 0;
        }
    } else if ((state.input.chkv == 1 && state.surfaceChokeMassCondition == 0) || (state.input.chkv == 1 && state.surfaceChokeMassCondition == 1)) {
        if (state.cells[state.lastCell - 1].MliqiniR < 0) {
            state.cells[state.lastCell - 1].MR = 0.;
            state.cells[state.lastCell - 1].MliqiniR = 0;
            state.cells[state.lastCell - 1].term1R = 0;
            state.cells[state.lastCell - 1].term2R = 0;
            state.cells[state.lastCell].MC = 0.;
            state.cells[state.lastCell].Mliqini = 0;
            state.cells[state.lastCell].term1 = 0;
            state.cells[state.lastCell].term2 = 0;
        }
        if (state.cells[state.lastCell - 1].QLR < 0 && (state.surfaceChokeMassCondition == 0 || state.surfaceChokeOpen == 1)) {
            state.cells[state.lastCell - 1].QLR = 0;
            state.cells[state.lastCell].QL = 0;
        }
    }

    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // atencao!!!!!!!!!!!!!!!!!!
    // existe uma questao que parece mal resolvida na resolucao desta marcha, nao foi feito nenhum teste para
    // o caso em que a velocidade de transporte da temperatura é <0 neste caso, a temperatura na celula de indice
    // não deveria entrar no metodo calctemp, ja que não e mais o caso de ser uma celula com condicao de
    // contorno para temperatura??????????????????????????????????????????????????????????????????????/
    //"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    state.cells[0].tempini = state.cells[0].temp;
    if (state.input.ConContEntrada > 0)
        state.cells[0].temp = state.inletTemperature;
    else
        state.cells[0].temp = state.defaultInletTemperature;
    state.cells[0].dTdt = (state.cells[0].temp - state.cells[0].tempini) / state.timeStep;
    state.cells[0].dTdtL = state.cells[0].dTdt;
    state.cells[1].tempLini = state.cells[1].tempL;
    state.cells[1].tempL = state.cells[0].temp;
    for (int i = 1; i <= state.lastCell; i++) {
        state.cells[i].tempini = state.cells[i].temp;
    }
    if (state.input.modoDifus3D == 0) {
        if (state.poisson2DCellCount > 1) {
#pragma omp parallel for num_threads((*state.globals).ntrd)
            for (int iP2D = 0; iP2D < state.poisson2DCellCount; iP2D++) {
                int i = state.poisson2DCellIndices[iP2D];
                if (i <= state.lastCell) {
                    computeTemperature(state, i, state.cells[i].tempini);
                } else {
                    if ((*state.globals).chaverede == 0 || state.networkEndpoint == 1 || (*state.globals).chaveRedeParalela == 1)
                        state.cells[i].temp = state.cells[i].calor.Textern1;
                    else
                        state.cells[i].temp = state.gasSurfaceTemperature;
                }
                state.cells[i].dTdt = (state.cells[i].temp - state.cells[i].tempini) / state.timeStep;
                state.cells[i].dTdtIni = state.cells[i].dTdt;
            }
        }
#pragma omp parallel for num_threads((*state.globals).ntrd)
        for (int i = 0; i <= state.lastCell; i++) {
            if (i == 217) {
                int para;
                para = 0;
            }
            if (state.cells[i].calor.difus2D == 0) {
                if (i <= state.lastCell) {
                    computeTemperature(state, i, state.cells[i].tempini);
                } else {
                    if ((*state.globals).chaverede == 0 || state.networkEndpoint == 1 || (*state.globals).chaveRedeParalela == 1)
                        state.cells[i].temp = state.cells[i].calor.Textern1;
                    else
                        state.cells[i].temp = state.gasSurfaceTemperature;
                }
                state.cells[i].dTdt = (state.cells[i].temp - state.cells[i].tempini) / state.timeStep;
                state.cells[i].dTdtIni = state.cells[i].dTdt;
            }
        }
    } else if (state.input.modoDifus3D == 1) {
#pragma omp parallel for num_threads((*state.globals).ntrd)
        for (int i = 0; i <= state.lastCell; i++) {
            int acoplado = -1;
            int icelAcop;
            for (int iacop = 0; iacop < state.input.nacop; iacop++) {
                icelAcop = state.input.celAcop[iacop].indCel;
                if (i == icelAcop) {
                    acoplado = iacop;
                    break;
                }
            }
            if (acoplado != -1) {
                prepareNonDimensionalHeatDiffusion(state, i);
                int iacop1 = state.coupledCellIndices[acoplado];
                state.poissonSolver.dados.tInt[iacop1] = state.cells[icelAcop].temp;
                double hiCel = state.cells[icelAcop].calor.hInt();
                state.poissonSolver.dados.hI[iacop1] = hiCel;
            }
        }
        if (ciclo < ciclomax || ciclomax == 0 || state.minimumCycleTimeStep != state.timeStep) {
            if (ciclo == ciclomax && ciclomax > 0)
                state.poissonSolver.FeiticoDoTempo();
            state.poissonSolver.transientePoisson(state.timeStep);
        }
        if (state.poisson2DCellCount > 1) {
#pragma omp parallel for num_threads((*state.globals).ntrd)
            for (int iP2D = 0; iP2D < state.poisson2DCellCount; iP2D++) {
                int i = state.poisson2DCellIndices[iP2D];
                if (i <= state.lastCell) {
                    computeTemperature(state, i, state.cells[i].tempini);
                } else {
                    if ((*state.globals).chaverede == 0 || state.networkEndpoint == 1 || (*state.globals).chaveRedeParalela == 1)
                        state.cells[i].temp = state.cells[i].calor.Textern1;
                    else
                        state.cells[i].temp = state.gasSurfaceTemperature;
                }
                state.cells[i].dTdt = (state.cells[i].temp - state.cells[i].tempini) / state.timeStep;
                state.cells[i].dTdtIni = state.cells[i].dTdt;
            }
        }
#pragma omp parallel for num_threads((*state.globals).ntrd)
        for (int i = 0; i <= state.lastCell; i++) {
            if (state.cells[i].calor.difus2D == 0) {
                if (i <= state.lastCell) {
                    computeTemperature(state, i, state.cells[i].tempini);
                } else {
                    if ((*state.globals).chaverede == 0 || state.networkEndpoint == 1 || (*state.globals).chaveRedeParalela == 1)
                        state.cells[i].temp = state.cells[i].calor.Textern1;
                    else
                        state.cells[i].temp = state.gasSurfaceTemperature;
                }
                state.cells[i].dTdt = (state.cells[i].temp - state.cells[i].tempini) / state.timeStep;
                state.cells[i].dTdtIni = state.cells[i].dTdt;
            }
        }
    }
    for (int i = 1; i <= state.lastCell; i++) {
        state.cells[i].dTdtL = state.cells[i - 1].dTdt;
        if (i < state.lastCell) {
            state.cells[i + 1].tempLini = state.cells[i + 1].tempL;
            state.cells[i + 1].tempL = state.cells[i].temp;
        }
        state.cells[i - 1].tempRini = state.cells[i - 1].tempR;
        state.cells[i - 1].tempR = state.cells[i].temp;
    }
    if (ciclo < ciclomax) {
        for (int k = 0; k <= state.lastCell; k++) {
            state.cells[k].FeiticoDoTempo();
        }
    } else if (state.input.modoDifus3D == 1)
        state.poissonSolver.renova();
    if (state.completeModel == 0) {
        if (ciclo < ciclomax) {
            state.evolutionUpdater.solvePressureVelocityCoupling(ciclo);
            state.evolutionUpdater.renew();
        }
    }
}

void advanceSteadyTemperature(const ThermalState &state, int i, int RK) {
    double dx = 0.5 * (state.cells[i].dx + state.cells[i - 1].dx);
    double dxmed = 0.5 * (state.cells[i].dx + state.cells[i - 1].dx);
    double dTdLmed = (state.cells[i].dx * state.cells[i].dTdLCor + state.cells[i - 1].dx * state.cells[i - 1].dTdLCor) / (state.cells[i].dx + state.cells[i - 1].dx);
    double dia = state.cells[i - 1].duto.a;
    double area = 0.25 * M_PI * dia * dia;
    double alfmed;
    double betmed;
    double pmed;
    double tmed;
    if (RK == 0) {
        alfmed = state.cells[i - 1].alf;
        betmed = state.cells[i - 1].bet;
        pmed = state.cells[i - 1].pres;
        tmed = state.cells[i - 1].temp;
    } else {
        alfmed = state.cells[i].alf;
        betmed = state.cells[i].bet;
        pmed = state.cells[i].pres;
        tmed = state.cells[i].temp;
    }
    double ugsmed;
    ugsmed = state.cells[i].QG / area; // velocidade superficial de gas
    double ulsmed;
    ulsmed = state.cells[i].QL / area; // velocidade superficial de liquido
    double sinalJ = 1.;
    if (fabs(ugsmed + ulsmed) > 0.05 && state.thermalSourceDisabled == 0) { // calculo termico e feito para velocidades de mistura superiores a 0,1 m/s,
        // para velocidades inferiores se admite que a temperatura do fluidoÃƒÂ© igual ÃƒÂ  temperatura ambiente
        sinalJ = (ugsmed + ulsmed) / fabs(ugsmed + ulsmed);
        double pmedi = state.cells[i].presaux + 0 * state.cells[i - 1].dpB / 98066.5;
        double tmedi;
        if (state.steadyIteration != 0 && state.input.AceleraConvergPerm == 0) // temperatura na interface esquerda da celula
            // caso em que se considera que ja foi
            // feita uma iteracao e ja se tem a temperatura na celula i vinda da iteracao anterior
            // isto pode dificultar a convergencia, quando se deseja a aceleracao da convergencia
            // admite-se que a temperatura na fronteira esquerda ÃƒÂ© a temperatura da celula esquerda
            tmedi = (state.cells[i].dx * state.cells[i].temp + state.cells[i].dxL * state.cells[i].tempL) / (state.cells[i].dx + state.cells[i].dxL);
        else
            tmedi = state.cells[i - 1].temp;
        double rp = state.cells[i - 1].flui.MasEspLiq(pmed, tmed);    // celula[i].rpCi;
        double rc = state.cells[i - 1].fluicol.MasEspFlu(pmed, tmed); // celula[i].rcCi;
        double rhol = (1. - betmed) * rp + betmed * rc;
        double rhog = state.cells[i - 1].flui.MasEspGas(pmed, tmed); // celula[i].rgCi;
        double cpl = (1. - betmed) * state.cells[i - 1].flui.CalorLiq(pmedi, tmedi) + betmed * state.cells[i - 1].fluicol.CalorLiq(pmedi, tmedi);
        double cpg = state.cells[i - 1].flui.CalorGas(pmedi, tmedi);
        // jtl=Joule Thomson do liquido X cp
        // jtg=Joule Thomson do gas X cp
        /////??????????????????????????????????????????????????????????????????????????????????????????
        double jtl = (1. - betmed) * state.cells[i - 1].flui.JTL(pmedi, tmedi) - betmed / rc;
        /////??????????????????????????????????????????????????????????????????????????????????????????
        if (state.input.pocinjec > 0 && state.input.condpocinj.tipoFlui == 2) {
            jtl = -(1 + (tmedi + 273.14) * state.cells[i - 1].fluicol.DrhoDtFlu(pmedi, tmedi) / rc) / rc;
        }
        double jtg = state.cells[i - 1].flui.JTG(pmedi, tmedi);
        // energia potencial:
        double hidro = (rhol * ulsmed + rhog * ugsmed) * area * 9.82 * sin(state.cells[i - 1].duto.teta);

        if (i > 120) {
            int para;
            para = 0;
        }

        // definicao dos parametros internos na tubulacao para se obter a troca termica om o meio ambiente
        state.cells[i - 1].calor.Tint = tmedi;
        state.cells[i - 1].calor.Vint = fabs(ugsmed + ulsmed);
        double condliq = (1. - betmed) * state.cells[i - 1].flui.CondLiq(pmedi, tmedi) + betmed * state.cells[i - 1].fluicol.CondLiq(pmedi, tmedi);
        state.cells[i - 1].calor.kint = condliq * (1 - alfmed) + state.cells[i - 1].flui.CondGas(pmedi, tmedi) * alfmed;
        state.cells[i - 1].calor.cpint = cpl * (1 - alfmed) + cpg * alfmed;
        state.cells[i - 1].calor.rhoint = rhol * (1 - alfmed) + rhog * alfmed;
        double viscliq = (1. - betmed) * state.cells[i - 1].flui.ViscOleo(pmedi, tmedi) + betmed * state.cells[i - 1].fluicol.VisFlu(pmedi, tmedi);
        state.cells[i - 1].calor.viscint = viscliq * (1 - alfmed) * 1.e-3 + state.cells[i - 1].flui.ViscGas(pmedi, tmedi) * alfmed * 1.e-3;

        double fluxrelax = 0; // variavel nao utilizada
        if (state.steadyIteration > 0)
            fluxrelax = state.cells[i - 1].fluxcalmed;
        double resanul = 0.;
        double fluxcal;
        double fluxcalG;
        // verifica se existe acoplamento com o anular:
        if (state.input.lingas == 1 && (i - 1 <= state.annulusTubingStart && i - 1 >= state.annulusTubingEnd)) {
            int j = state.annulusTubingStart + state.tubingAnnulusStart - (i - 1);
            // caso tenha acoplamento:

            if (state.steadyIteration == 0 && (i - 1) < state.annulusTubingStart) {
                // na primeira iteracao, considera-se a resistencia do revestimento + cimento +
                // formacao, nas outras iteracoes
                // bastara determinar a troca termica entre a coluna e o gas do anular
                // no modelo acoplado, a coluna nao tem a definicao da parede revestimento+cimento+formacao
                // para se obter esta resistencia, precisa-se recorrer ao modelo de troca termica
                // do anular, o que e feito aqui:
                // obs:isto Ã© feito ate uma celula antes de se chegar na master (i-1)<ColunaAnulaIni,
                // se esta fazendo igual ao que se faz no simulador involuta, para melhorar
                // a estimativa da temperatura na ANM. Na celula da anm, so se considera a troca termica com o gas
                // desde a primeira iteracao
                state.gasCells[j].calor.Tint = tmedi;
                state.gasCells[j].calor.Vint = 100;
                state.gasCells[j].calor.kint = state.cells[i - 1].flui.CondGas(pmedi, tmedi);
                state.gasCells[j].calor.cpint = cpg;
                state.gasCells[j].calor.rhoint = rhog;
                state.gasCells[j].calor.viscint = state.cells[i - 1].flui.ViscGas(pmedi, tmedi) * 1.e-3;
                state.gasCells[j].fluxcal = state.gasCells[j].calor.transperm(); // troca termica no anular
                resanul = state.gasCells[j].calor.resGlob;                // resistencia das paredes
                // observe que nÃ£o se esta de fato interessado no fluxo de calor, mas apenas em obter a resistencia
                // termica do conjunto de paredes a partir do revestimento em direcao aa formacao
                state.cells[i - 1].calor.Vextern1 = 100.;
                state.cells[i - 1].calor.kextern1 = state.gasCells[j].calor.kint;
                state.cells[i - 1].calor.cpextern1 = state.gasCells[j].calor.cpint;
                state.cells[i - 1].calor.rhoextern1 = state.gasCells[j].calor.rhoint;
                state.cells[i - 1].calor.viscextern1 = state.gasCells[j].calor.viscint;
            } else { // apos a primeira iteracao, considera-se apenas a troca termica entre a coluna e o gas do anular
                // passando pela parede da coluna, claro
                resanul = 0.;
                state.cells[i - 1].calor.Vextern1 = state.gasCells[j].VGasR / state.gasCells[j].u1L;
                state.cells[i - 1].calor.kextern1 = state.gasCells[j].calor.kint;
                state.cells[i - 1].calor.cpextern1 = state.gasCells[j].calor.cpint;
                state.cells[i - 1].calor.rhoextern1 = state.gasCells[j].calor.rhoint;
                state.cells[i - 1].calor.viscextern1 = state.gasCells[j].calor.viscint;
            }
            if (state.steadyIteration == 0)
                state.cells[i - 1].calor.Textern1 = state.gasCells[j].calor.Textern1; // na primeira iteracao, como se
            // usa toda a resistencia termica do poco, a temperatura externa utilizada Ã© a geotermica
            else
                state.cells[i - 1].calor.Textern1 = state.gasCells[j].temp; // nas iteracoes seguintes, a temperatura ambiente
            // e a temperatura do gas
        }
        state.cells[i - 1].fluxcalmed = 0;
        if (state.productionNetworkHeatCoupled == 1 && (i - 1) >= state.primaryNetworkSectionEnd && (i - 1) <= state.primaryNetworkSectionStart) {
            fluxcal = sinalJ * state.cells[i - 1].calor.transperm(state.cells[i - 1].resAcopRedeP);
        } else
            fluxcal = sinalJ * state.cells[i - 1].calor.transperm(resanul);
        state.cells[i - 1].fluxcalmed = fluxcal; // fluxo de calor na coluna

        double coefdxT = (rhol * ulsmed * cpl + rhog * ugsmed * cpg) * area; // termo que multiplica
        // a derivada Dt/Dx
        double ulsmedTemp = ulsmed;
        if ((*state.globals).blackOilTemp == 1 && fabs(ulsmed) > 5)
            ulsmedTemp = 5 * ulsmed / ulsmed;
        double ugsmedTemp = ugsmed;
        if ((*state.globals).blackOilTemp == 1 && fabs(ugsmed) > 5)
            ugsmedTemp = 5 * ugsmed / fabs(ugsmed);
        double coefdxP = 1 * (rhol * ulsmedTemp * jtl + rhog * ugsmedTemp * jtg) * area; // termo que multiplica
        // a derivada Dp/Dx
        double dpdx;
        if ((state.cells[i - 1].acsr.tipo != 4 || state.cells[i - 1].acsr.bcs.freq < 1) && state.cells[i - 1].acsr.tipo != 7)
            // caso nao tenha BCS ou incremento de pressao  utiliza-se a pressao na fronteira esquerda
            //  e a pressao no centro de celula para o calculo de Dp/Dx
            dpdx = 2. * (state.cells[i].presaux - state.cells[i - 1].pres) * 98066.5 / state.cells[i - 1].dx;
        else {
            // caso tenha BCS ou incremento de pressao  utiliza-se a pressao da celulaa esquerda
            //  e a pressao no centro de celula para o calculoi de Dp/Dx
            dpdx = 2. * (state.cells[i].presaux - state.cells[i - 1].pres) * 98066.5 / state.cells[i - 1].dx;
        }
        state.cells[i].VTemper = ulsmed; // esta velocidade so e util no caso transiente, Ã© armazenada aqui
        // apenas para se ter um valor quando a simulacao transiente se iniciar
        double dtdx = (-state.cells[i - 1].temp) / dxmed;

        double cinetico = 0;
        double ugmed0 = 0;
        double ulmed0 = 0;
        double ugmed = 0;
        double ulmed = 0;
        // termo de energia cinetica:
        if (state.cells[i].acsr.tipo == 0 && state.cells[i - 1].acsr.tipo == 0 && i > 2) {
            double dxCin = state.cells[i - 1].dx;
            double diaaux = state.cells[i - 1].duto.a;
            double areaaux = 0.25 * M_PI * diaaux * diaaux;

            if (state.cells[i - 1].alf > 1e-3)
                ugmed = ugsmed / state.cells[i - 1].alf;
            if (state.cells[i - 1].alf < (1. - 1e-3))
                ulmed = ulsmed / (1. - state.cells[i - 1].alf);

            if (state.cells[i - 2].alf > 1e-3) {
                ugmed0 = state.cells[i - 1].QG / (areaaux);
                ugmed0 /= state.cells[i - 2].alf;
            }
            if (state.cells[i - 2].alf < (1. - 1e-3)) {
                ulmed0 = state.cells[i - 1].QL / (areaaux);
                ulmed0 /= (1. - state.cells[i - 2].alf);
            }

            if (state.input.nCompTotalUnidadesP / dxCin < 1e6)
                cinetico = (state.cells[i].MC - state.cells[i].Mliqini) * ugmed * (ugmed - ugmed0) / dxCin + state.cells[i].Mliqini * ulmed * (ulmed - ulmed0) / dxCin;
            else
                cinetico = 0;
        }

        double fontemassG = 0.;
        double fontemassL = 0.;
        double fontemassC = 0.;
        double tfonte = state.cells[i - 1].temp;
        double cpgF;
        double razcpF = 0.;
        double cplF;

        // calculo da energia adicionada no sistema devido a fontes de massa
        if (state.cells[i - 1].acsr.tipo == 1) { // caso fonte de gas
            tfonte = state.cells[i - 1].acsr.injg.temp;
            cpgF = state.cells[i - 1].acsr.injg.FluidoPro.CalorGas(pmed, tfonte);
            razcpF = state.cells[i - 1].acsr.injg.FluidoPro.ConstAdG(pmed, tfonte);
            cplF = 0.;
        } else if (state.cells[i - 1].acsr.tipo == 2) { // caso fonte de liquido
            tfonte = state.cells[i - 1].acsr.injl.temp;
            cpgF = 0.;
            razcpF = 1.;
            cplF = (1. - state.cells[i - 1].acsr.injl.bet) * state.cells[i - 1].acsr.injl.FluidoPro.CalorLiq(pmed, tmed) + state.cells[i - 1].acsr.injl.bet * state.cells[i - 1].acsr.injl.fluidocol.CalorLiq(pmed, tmed);
        } else if (state.cells[i - 1].acsr.tipo == 3) { // caso IPR
            tfonte = state.cells[i - 1].acsr.ipr.Tres;
            cpgF = state.cells[i - 1].acsr.ipr.FluidoPro.CalorGas(pmed, tmed);
            razcpF = state.cells[i - 1].acsr.ipr.FluidoPro.ConstAdG(pmed, tmed);
            cplF = state.cells[i - 1].acsr.ipr.FluidoPro.CalorLiq(pmed, tmed);
        } else if (state.cells[i - 1].acsr.tipo == 9 && state.cells[i - 1].acsr.fontechk.abertura > 1e-6 &&
                   (state.cells[i - 1].fontemassCR + state.cells[i - 1].fontemassGR + state.cells[i - 1].fontemassLR) > 1e-9) {
            // caso vazamento
            tfonte = state.cells[i - 1].acsr.fontechk.tamb;
            cpgF = state.cells[i - 1].acsr.fontechk.fluidoPamb.CalorGas(pmed, tmed);
            razcpF = state.cells[i - 1].acsr.fontechk.fluidoPamb.ConstAdG(pmed, tmed);
            cplF = (1. - state.cells[i - 1].acsr.fontechk.betISamb) *
                       state.cells[i - 1].acsr.fontechk.fluidoPamb.CalorLiq(pmed, tmed) +
                   state.cells[i - 1].acsr.fontechk.betISamb * state.cells[i - 1].acsr.fontechk.fluidocol.CalorLiq(pmed, tmed);
        } else if (state.cells[i - 1].acsr.tipo == 15) { // caso IPR
            tfonte = state.cells[i - 1].acsr.radialPoro.tRes;
            cpgF = state.cells[i - 1].acsr.radialPoro.flup.CalorGas(pmed, tmed);
            razcpF = state.cells[i - 1].acsr.radialPoro.flup.ConstAdG(pmed, tmed);
            cplF = state.cells[i - 1].acsr.radialPoro.flup.CalorLiq(pmed, tmed);
        } else if (state.cells[i - 1].acsr.tipo == 16) { // caso IPR
            tfonte = state.cells[i - 1].acsr.poroso2D.dados.tRes;
            cpgF = state.cells[i - 1].acsr.poroso2D.dados.flup.CalorGas(pmed, tmed);
            razcpF = state.cells[i - 1].acsr.poroso2D.dados.flup.ConstAdG(pmed, tmed);
            cplF = state.cells[i - 1].acsr.poroso2D.dados.flup.CalorLiq(pmed, tmed);
        } else if (state.cells[i].acsrL != 0) {

            cpgF = 0.;
            razcpF = 1.;
            cplF = 0.;
        } else {
            cpgF = 0.;
            razcpF = 1.;
            cplF = 0.;
        }

        fontemassL = 0;
        if (state.cells[i - 1].fontemassLR > 0.)
            fontemassL = state.cells[i - 1].fontemassLR / dx;
        if (state.cells[i - 1].fontemassCR > 0.)
            fontemassL += state.cells[i - 1].fontemassCR / dx;
        fontemassL *= (cplF) * (tfonte - state.cells[i - 1].temp);

        fontemassG = state.cells[i - 1].fontemassGR / dx;
        if (fontemassG > 0.)
            fontemassG *= (cpgF / razcpF) * (tfonte - state.cells[i - 1].temp);
        else
            fontemassG = 0;

        // efeito do calor latente, quando este for solicitado
        double latente = 0.;
        if (isnan(state.cells[i - 1].FonteMudaFase))
            state.cells[i - 1].FonteMudaFase = 0.;
        double valTransMass = fabs(state.cells[i - 1].FonteMudaFase);
        double sigTransMass = 1.;
        if (valTransMass > 1e-25)
            sigTransMass = state.cells[i - 1].FonteMudaFase / valTransMass;
        if (state.input.limTransMass < valTransMass)
            valTransMass = sigTransMass * state.input.limTransMass;
        else
            valTransMass *= sigTransMass;
        if (state.cells[i].flui.dVaporMassFraction < (1 - 1e-15) && state.cells[i].flui.dVaporMassFraction > (1e-15)) {
            if (state.cells[i - 1].acsr.tipo == 1 || state.cells[i - 1].acsr.tipo == 2 || state.cells[i - 1].acsr.tipo == 3 || state.cells[i - 1].acsr.tipo == 15 || state.cells[i - 1].acsr.tipo == 16) {
                state.cells[i - 1].FonteMudaFase = 0.;
                valTransMass = 0.;
            }
            if (state.latentHeatEnabled > 0 && state.steadyIteration != 0 && state.input.flashCompleto == 0) {
                latente = interpolateLatentHeat(state, pmed, tmed) * valTransMass;
            } else if ((state.input.flashCompleto == 1 || state.input.flashCompleto == 2) && state.latentHeatEnabled > 0 && state.steadyIteration != 0) {
                latente = (state.cells[i].flui.EntalpGas(pmed, tmed) -
                           state.cells[i].flui.EntalpLiq(pmed, tmed)) *
                          valTransMass;
            } else
                latente = 0;
        }

        double alfinter;
        double alfinterL;
        if (ugsmed > 0) {
            alfinter = state.cells[i - 1].alf;
            alfinterL = state.cells[i - 1].alfL;
        } else {
            alfinter = state.cells[i - 1].alfR;
            alfinterL = state.cells[i - 1].alf;
        }
        double delvel;
        if (alfinter > (*state.globals).localtiny && alfinter < (1. - (*state.globals).localtiny))
            delvel = ugsmed / alfinter - ulsmed / (1. - alfinter);
        else if (alfinter > (*state.globals).localtiny)
            delvel = ugsmed;
        else
            delvel = ulsmed;
        double verifica = area * state.cells[i].pres * 98600 * delvel * (alfinter - alfinterL) / state.cells[i].dx;

        if (fabs(coefdxT) > (*state.globals).localtiny) {
            // Energy terms for boundary work, potential and kinetic energy,
            // mass sources, latent heat, and shaft work.
            if (state.input.latente == 0)
                latente = 0.;
            else if (state.input.condlatente == 0 && latente < 0)
                latente = 0.;
            double parcenerg1 = dTdLmed * (coefdxP * dpdx - cinetico - (hidro) + (fontemassL + fontemassG) - latente + (state.cells[i - 1].potTermo + state.cells[i - 1].fonteCal) / dxmed) / coefdxT;
            // portion of energy related to heat exchange
            double parcenerg2 = dTdLmed * (fluxcal) / coefdxT;

            // Check whether heat transfer is too fast for the explicit temperature update.
            // If so, use additional substeps to prevent thermal instability.
            int npasso;
            double dxpasso;
            double limiteEstab = fabs(coefdxT) / dTdLmed;
            if (dxmed / (state.cells[i - 1].calor.resGlob + resanul) < (limiteEstab + 0. * 1000.)) { // Thermal resistance is not low.
                npasso = 1;
                dxpasso = dxmed;
            } else {
                // For low thermal resistance, determine the number of temperature update steps.
                npasso = (dxmed / (state.cells[i - 1].calor.resGlob + resanul)) / (limiteEstab + 0 * 1000) + 1;
                dxpasso = dxmed / npasso; // cell divided into npassos
            }
            double temppasso = state.cells[i - 1].temp;

            temppasso = dxpasso * (-(-state.cells[i - 1].temp) / dxpasso + parcenerg1 + parcenerg2); // first step
            for (int j = 1; j < npasso; j++) {                                                  // next steps
                // Keep all energy terms except heat flow, which is recalculated at each step.
                state.cells[i - 1].calor.Tint = temppasso;
                state.cells[i - 1].calor.Vint = ugsmed + ulsmed;
                condliq = (1. - betmed) * state.cells[i - 1].flui.CondLiq(pmed, temppasso) + betmed * state.cells[i - 1].fluicol.CondLiq(pmed, temppasso);
                state.cells[i - 1].calor.kint = condliq * (1 - alfmed) + state.cells[i - 1].flui.CondGas(pmed, temppasso) * alfmed;
                state.cells[i - 1].calor.cpint = cpl * (1 - alfmed) + cpg * alfmed;
                state.cells[i - 1].calor.rhoint = rhol * (1 - alfmed) + rhog * alfmed;
                viscliq = (1. - betmed) * state.cells[i - 1].flui.ViscOleo(pmed, temppasso) + betmed * state.cells[i - 1].fluicol.VisFlu(pmed, temppasso);
                state.cells[i - 1].calor.viscint = viscliq * (1 - alfmed) * 1.e-3 + state.cells[i - 1].flui.ViscGas(pmed, temppasso) * alfmed * 1.e-3;
                if (state.steadyIteration != 0 && state.input.lingas == 1 && (i - 1 <= state.annulusTubingStart && i - 1 >= state.annulusTubingEnd)) {
                    int k = state.annulusTubingStart + state.tubingAnnulusStart - (i - 1);
                    double dtext = (state.gasCells[k - 1].temp - state.gasCells[k].temp) / npasso;
                    state.cells[i - 1].calor.Textern1 = state.gasCells[k].temp + (j - 1) * dtext;
                }
                fluxcal = sinalJ * state.cells[i - 1].calor.transperm(resanul);
                state.cells[i - 1].fluxcalmed += fluxcal;
                parcenerg2 = dTdLmed * (fluxcal) / coefdxT;
                temppasso = dxpasso * (-(-temppasso) / dxpasso + parcenerg1 + parcenerg2);
            }

            state.cells[i].temp = temppasso;
            state.cells[i - 1].fluxcalmed /= npasso;
            state.cells[i - 1].fluxcalmed = 1. * state.cells[i - 1].fluxcalmed;

            if (state.cells[i].temp < -50.)
                state.cells[i].temp = -50.;
            if (state.cells[i].temp > 200.)
                state.cells[i].temp = 200.;
        } else
            state.cells[i].temp = state.cells[i].calor.Textern1;
    } else { // case where the mixing speed is too low
        state.cells[i].temp = state.cells[i - 1].calor.Textern1;
        if (state.input.lingas == 1 && (i - 1 <= state.annulusTubingStart && i - 1 >= state.annulusTubingEnd)) {
            int j = state.annulusTubingStart + state.tubingAnnulusStart - (i - 1);
            state.cells[i - 1].calor.Vextern1 = 100.;
            state.cells[i - 1].calor.kextern1 = state.gasCells[j].calor.kint;
            state.cells[i - 1].calor.cpextern1 = state.gasCells[j].calor.cpint;
            state.cells[i - 1].calor.rhoextern1 = state.gasCells[j].calor.rhoint;
            state.cells[i - 1].calor.viscextern1 = state.gasCells[j].calor.viscint;
            double fluxcal = state.cells[i - 1].calor.transperm(0);
        }
    }
}

void advanceReverseSteadyTemperature(const ThermalState &state, int i, int RK) {
    double dx = 0.5 * (state.cells[i].dx + state.cells[i + 1].dx);
    double dxmed = 0.5 * (state.cells[i].dx + state.cells[i + 1].dx);
    double dTdLmed = (state.cells[i].dx * state.cells[i].dTdLCor + state.cells[i + 1].dx * state.cells[i + 1].dTdLCor) / (state.cells[i].dx + state.cells[i + 1].dx);
    double dia = state.cells[i + 1].duto.a;
    double area = 0.25 * M_PI * dia * dia;
    double alfmed;
    double betmed;
    double pmed;
    double tmed;
    if (RK == 0) {
        alfmed = state.cells[i + 1].alf;
        betmed = state.cells[i + 1].bet;
        pmed = state.cells[i + 1].pres;
        tmed = state.cells[i + 1].temp;
    } else {
        alfmed = state.cells[i].alf;
        betmed = state.cells[i].bet;
        pmed = state.cells[i].pres;
        tmed = state.cells[i].temp;
    }
    double ugsmed;
    ugsmed = fabs(state.cells[i + 1].QG) / area; // Superficial gas velocity
    double ulsmed;
    ulsmed = fabs(state.cells[i + 1].QL) / area; // Superficial liquid velocity
    double sinalJ = 1.;
    if (fabs(ugsmed + ulsmed) > state.slowHeatTransferThreshold) {
        // Perform the thermal calculation for mixture velocities above 0.1 m/s.
        // Otherwise, assume the fluid temperature equals the ambient temperature.
        sinalJ = (ugsmed + ulsmed) / fabs(ugsmed + ulsmed);
        double pmedi = state.cells[i + 1].presaux - state.cells[i].dpB / 98066.5;
        double tmedi;
        if (state.steadyIteration != 0 && state.input.AceleraConvergPerm == 0)
            // Set the left interface temperature from the previous iteration.
            // Use the adjacent left cell value to accelerate convergence.
            tmedi = (state.cells[i].dx * state.cells[i].temp + state.cells[i].dxR * state.cells[i].tempR) / (state.cells[i].dx + state.cells[i].dxR);
        else
            tmedi = state.cells[i + 1].temp;
        double rp = state.cells[i + 1].flui.MasEspLiq(pmed, tmed);
        double rc = state.cells[i + 1].fluicol.MasEspFlu(pmed, tmed);
        double rhol = (1. - betmed) * rp + betmed * rc;
        double rhog = state.cells[i + 1].flui.MasEspGas(pmed, tmed);
        double cpl = (1. - betmed) * state.cells[i + 1].flui.CalorLiq(pmedi, tmedi) + betmed * state.cells[i + 1].fluicol.CalorLiq(pmedi, tmedi);
        double cpg = state.cells[i + 1].flui.CalorGas(pmedi, tmedi);
        // Liquid Joule-Thomson coefficient multiplied by cp.
        // Gas Joule-Thomson coefficient multiplied by cp.p
        double jtl = (1. - betmed) * state.cells[i + 1].flui.JTL(pmedi, tmedi) - betmed / rc;
        if (state.input.pocinjec > 0 && state.input.condpocinj.tipoFlui == 2) {
            jtl = -(1 + (tmedi + 273.14) * state.cells[i + 1].fluicol.DrhoDtFlu(pmedi, tmedi) / rc) / rc;
        }
        double jtg = state.cells[i + 1].flui.JTG(pmedi, tmedi);
        // potential energy:
        double hidro = -(rhol * ulsmed + rhog * ugsmed) * area * 9.82 * sin(state.cells[i + 1].duto.teta);

        if (i > 120) {
            int para;
            para = 0;
        }

        // Definition of internal parameters in the piping to achieve heat exchange with the environment
        state.cells[i + 1].calor.Tint = tmedi;
        state.cells[i + 1].calor.Vint = fabs(ugsmed + ulsmed);
        double condliq = (1. - betmed) * state.cells[i + 1].flui.CondLiq(pmedi, tmedi) + betmed * state.cells[i + 1].fluicol.CondLiq(pmedi, tmedi);
        state.cells[i + 1].calor.kint = condliq * (1 - alfmed) + state.cells[i + 1].flui.CondGas(pmedi, tmedi) * alfmed;
        state.cells[i + 1].calor.cpint = cpl * (1 - alfmed) + cpg * alfmed;
        state.cells[i + 1].calor.rhoint = rhol * (1 - alfmed) + rhog * alfmed;
        double viscliq = (1. - betmed) * state.cells[i + 1].flui.ViscOleo(pmedi, tmedi) + betmed * state.cells[i + 1].fluicol.VisFlu(pmedi, tmedi);
        state.cells[i + 1].calor.viscint = viscliq * (1 - alfmed) * 1.e-3 + state.cells[i + 1].flui.ViscGas(pmedi, tmedi) * alfmed * 1.e-3;

        double fluxrelax = 0;
        if (state.steadyIteration > 0)
            fluxrelax = state.cells[i + 1].fluxcalmed;
        double resanul = 0.;
        double fluxcal;
        double fluxcalG;
        // checks if coupling with the annular space exists:
        if (state.input.lingas == 1 && (i + 1 <= state.annulusTubingStart && i + 1 >= state.annulusTubingEnd)) {
            int j = state.annulusTubingStart + state.tubingAnnulusStart - (i + 1);
            // in case there is coupling:

            if (state.steadyIteration == 0 && (i + 1) > state.annulusTubingEnd) {
                // On the first iteration, include casing, cement, and formation resistance.
                // Later iterations consider only heat transfer between the tubing and annular gas.
                // For the coupled model, obtain the external resistance from the annulus heat-transfer model.
                // Apply this up to the cell before the master valve to improve the ANM temperature estimate.
                // At the ANM cell, consider only heat transfer with the gas from the first iteration.
                state.gasCells[j].calor.Tint = tmedi;
                state.gasCells[j].calor.Vint = 100;
                state.gasCells[j].calor.kint = state.cells[i + 1].flui.CondGas(pmedi, tmedi);
                state.gasCells[j].calor.cpint = cpg;
                state.gasCells[j].calor.rhoint = rhog;
                state.gasCells[j].calor.viscint = state.cells[i + 1].flui.ViscGas(pmedi, tmedi) * 1.e-3;
                state.gasCells[j].fluxcal = state.gasCells[j].calor.transperm(); // troca termica no anular
                resanul = state.gasCells[j].calor.resGlob;                // resistencia das paredes
                // Only the thermal resistance from the casing to the formation is needed, not the actual heat flow.
                state.cells[i + 1].calor.Vextern1 = 100.;
                state.cells[i + 1].calor.kextern1 = state.gasCells[j].calor.kint;
                state.cells[i + 1].calor.cpextern1 = state.gasCells[j].calor.cpint;
                state.cells[i + 1].calor.rhoextern1 = state.gasCells[j].calor.rhoint;
                state.cells[i + 1].calor.viscextern1 = state.gasCells[j].calor.viscint;
            } else { // After the first iteration, consider only heat transfer between the tubing and annular gas.
                // Through the tubing wall.
                resanul = 0.;
                state.cells[i + 1].calor.Vextern1 = state.gasCells[j].VGasR / state.gasCells[j].u1L;
                state.cells[i + 1].calor.kextern1 = state.gasCells[j].calor.kint;
                state.cells[i + 1].calor.cpextern1 = state.gasCells[j].calor.cpint;
                state.cells[i + 1].calor.rhoextern1 = state.gasCells[j].calor.rhoint;
                state.cells[i + 1].calor.viscextern1 = state.gasCells[j].calor.viscint;
            }
            if (state.steadyIteration == 0)
                state.cells[i + 1].calor.Textern1 = state.gasCells[j].calor.Textern1; // na primeira iteracao, como se
            // usa toda a resistencia termica do poco, a temperatura externa utilizada Ã© a geotermica
            else
                state.cells[i + 1].calor.Textern1 = state.gasCells[j].temp; // nas iteracoes seguintes, a temperatura ambiente
            // e a temperatura do gas
        }
        state.cells[i + 1].fluxcalmed = 0;
        fluxcal = sinalJ * state.cells[i + 1].calor.transperm(resanul);
        state.cells[i + 1].fluxcalmed = fluxcal; // fluxo de calor na coluna

        double coefdxT = (rhol * ulsmed * cpl + rhog * ugsmed * cpg) * area; // termo que multiplica
        // a derivada Dt/Dx
        double coefdxP = 1. * (rhol * ulsmed * jtl + rhog * ugsmed * jtg) * area; // termo que multiplica
        // a derivada Dp/Dx
        double dpdx;
        if ((state.cells[i + 1].acsr.tipo != 4 || state.cells[i + 1].acsr.bcs.freq < 1) && state.cells[i + 1].acsr.tipo != 7)
            // caso nao tenha BCS ou incremento de pressao  utiliza-se a pressao na fronteira esquerda
            //  e a pressao no centro de celula para o calculo de Dp/Dx
            dpdx = 2. * (state.cells[i + 1].presaux - state.cells[i + 1].pres) * 98066.5 / state.cells[i + 1].dx;
        else {
            // caso tenha BCS ou incremento de pressao  utiliza-se a pressao da celulaa esquerda
            //  e a pressao no centro de celula para o calculoi de Dp/Dx
            dpdx = (pmedi - state.cells[i + 1].pres) * 98600. / dx;
        }
        state.cells[i].VTemper = ulsmed; // esta velocidade so e util no caso transiente, Ã© armazenada aqui
        // apenas para se ter um valor quando a simulacao transiente se iniciar
        double dtdx = (-state.cells[i + 1].temp) / dxmed;

        double cinetico = 0;
        double ugmed0 = 0;
        double ulmed0 = 0;
        double ugmed = 0;
        double ulmed = 0;
        // termo de energia cinetica:
        if (state.cells[i].acsr.tipo == 0 && state.cells[i + 1].acsr.tipo == 0 && i < state.lastCell) {
            double dxCin = state.cells[i + 1].dx;
            double diaaux = state.cells[i + 1].duto.a;
            double areaaux = 0.25 * M_PI * diaaux * diaaux;

            if (state.cells[i + 1].alf > 1e-3)
                ugmed = ugsmed / state.cells[i + 1].alf;
            if (state.cells[i + 1].alf < (1. - 1e-3))
                ulmed = ulsmed / (1. - state.cells[i + 1].alf);

            if (state.cells[i].alf > 1e-3) {
                ugmed0 = fabs(state.cells[i].QG) / (areaaux);
                ugmed0 /= state.cells[i].alf;
            }
            if (state.cells[i].alf < (1. - 1e-3)) {
                ulmed0 = fabs(state.cells[i].QL) / (areaaux);
                ulmed0 /= (1. - state.cells[i].alf);
            }

            cinetico = -fabs(state.cells[i].MC - state.cells[i].Mliqini) * ugmed * (ugmed - ugmed0) / dxCin -
                       fabs(state.cells[i].Mliqini) * ulmed * (ulmed - ulmed0) / dxCin;
        }

        double fontemassG = 0.;
        double fontemassL = 0.;
        double tfonte = state.cells[i + 1].temp;
        double cpgF;
        double razcpF = 0.;
        double cplF;

        // calculo da energia adicionada no sistema devido a fontes de massa
        if (state.cells[i + 1].acsr.tipo == 1) { // caso fonte de gas
            tfonte = state.cells[i + 1].acsr.injg.temp;
            cpgF = state.cells[i + 1].acsr.injg.FluidoPro.CalorGas(pmed, tfonte);
            razcpF = state.cells[i + 1].acsr.injg.FluidoPro.ConstAdG(pmed, tfonte);
            cplF = 0.;
        } else if (state.cells[i + 1].acsr.tipo == 2) { // caso fonte de liquido
            tfonte = state.cells[i + 1].acsr.injl.temp;
            cpgF = 0.;
            razcpF = 1.;
            cplF = (1. - state.cells[i + 1].acsr.injl.bet) * state.cells[i + 1].acsr.injl.FluidoPro.CalorLiq(pmed, tmed) + state.cells[i + 1].acsr.injl.bet * state.cells[i + 1].acsr.injl.fluidocol.CalorLiq(pmed, tmed);
        } else if (state.cells[i + 1].acsr.tipo == 3) { // caso IPR
            tfonte = state.cells[i + 1].acsr.ipr.Tres;
            cpgF = state.cells[i + 1].acsr.ipr.FluidoPro.CalorGas(pmed, tmed);
            razcpF = state.cells[i + 1].acsr.ipr.FluidoPro.ConstAdG(pmed, tmed);
            cplF = state.cells[i + 1].acsr.ipr.FluidoPro.CalorLiq(pmed, tmed);
        } else if (state.cells[i + 1].acsr.tipo == 15) { // caso IPR
            tfonte = state.cells[i + 1].acsr.radialPoro.tRes;
            cpgF = state.cells[i + 1].acsr.radialPoro.flup.CalorGas(pmed, tmed);
            razcpF = state.cells[i + 1].acsr.radialPoro.flup.ConstAdG(pmed, tmed);
            cplF = state.cells[i + 1].acsr.radialPoro.flup.CalorLiq(pmed, tmed);
        } else if (state.cells[i + 1].acsr.tipo == 16) { // caso IPR
            tfonte = state.cells[i + 1].acsr.poroso2D.dados.tRes;
            cpgF = state.cells[i + 1].acsr.poroso2D.dados.flup.CalorGas(pmed, tmed);
            razcpF = state.cells[i + 1].acsr.poroso2D.dados.flup.ConstAdG(pmed, tmed);
            cplF = state.cells[i + 1].acsr.poroso2D.dados.flup.CalorLiq(pmed, tmed);
        } else if (state.cells[i + 1].acsr.tipo == 9 && state.cells[i + 1].acsr.fontechk.abertura > 1e-6 &&
                   (state.cells[i + 1].fontemassCR + state.cells[i + 1].fontemassGR + state.cells[i + 1].fontemassLR) > 1e-9) {
            // caso vazamento
            tfonte = state.cells[i + 1].acsr.fontechk.tamb;
            cpgF = state.cells[i + 1].acsr.fontechk.fluidoPamb.CalorGas(pmed, tmed);
            razcpF = state.cells[i + 1].acsr.fontechk.fluidoPamb.ConstAdG(pmed, tmed);
            cplF = (1. - state.cells[i + 1].acsr.fontechk.betISamb) *
                       state.cells[i + 1].acsr.fontechk.fluidoPamb.CalorLiq(pmed, tmed) +
                   state.cells[i + 1].acsr.fontechk.betISamb * state.cells[i + 1].acsr.fontechk.fluidocol.CalorLiq(pmed, tmed);
        } else if (state.cells[i + 1].acsr.tipo != 0) {

            cpgF = 0.;
            razcpF = 1.;
            cplF = 0.;
        } else {
            cpgF = 0.;
            razcpF = 1.;
            cplF = 0.;
        }

        fontemassL = 0;
        if (state.cells[i + 1].fontemassLR > 0.)
            fontemassL = state.cells[i + 1].fontemassLR / dx;
        if (state.cells[i + 1].fontemassCR > 0.)
            fontemassL += state.cells[i + 1].fontemassCR / dx;
        fontemassL *= (cplF) * (tfonte - state.cells[i + 1].temp);

        fontemassG = state.cells[i + 1].fontemassGR / dx;
        if (fontemassG > 0.)
            fontemassG *= (cpgF / razcpF) * (tfonte - state.cells[i + 1].temp);
        else
            fontemassG = 0;

        // efeito do calor latente, quando este for solicitado
        double latente = 0.;
        if (state.cells[i].flui.dVaporMassFraction < (1 - 1e-15) && state.cells[i].flui.dVaporMassFraction > (1e-15)) {
            if (state.cells[i + 1].acsr.tipo == 1 || state.cells[i + 1].acsr.tipo == 2 || state.cells[i + 1].acsr.tipo == 3 || state.cells[i + 1].acsr.tipo == 15 || state.cells[i + 1].acsr.tipo == 16)
                state.cells[i + 1].FonteMudaFase =
                    0.;
            if (state.latentHeatEnabled > 0 && state.steadyIteration != 0 && state.input.flashCompleto == 0) {
                latente = -interpolateLatentHeat(state, pmed, tmed) * state.cells[i + 1].FonteMudaFase;
            } else if ((state.input.flashCompleto == 1 || state.input.flashCompleto == 2) && state.latentHeatEnabled > 0 && state.steadyIteration != 0) {
                latente = -(state.cells[i].flui.EntalpGas(pmed, tmed) -
                            state.cells[i].flui.EntalpLiq(pmed, tmed)) *
                          state.cells[i + 1].FonteMudaFase;
            } else
                latente = 0;
        }

        //////trecho sem utilidade////////////////////////////////////////////////////////
        double alfinter;
        double alfinterL;
        if (ugsmed > 0) {
            alfinter = state.cells[i + 1].alf;
            alfinterL = state.cells[i + 1].alfL;
        } else {
            alfinter = state.cells[i + 1].alfR;
            alfinterL = state.cells[i + 1].alf;
        }
        double delvel;
        if (alfinter > (*state.globals).localtiny && alfinter < (1. - (*state.globals).localtiny))
            delvel = ugsmed / alfinter - ulsmed / (1. - alfinter);
        else if (alfinter > (*state.globals).localtiny)
            delvel = ugsmed;
        else
            delvel = ulsmed;
        double verifica = area * state.cells[i + 1].pres * 98600 * delvel * (alfinter - alfinterL) / state.cells[i + 1].dx;

        if (fabs(coefdxT) > (*state.globals).localtiny) {
            // parecela de energia relacionada ao ytrabalho de fronteira, energia potencial,
            // energia cinetica, fontes de massa, calor latente e trabalho de eixo:
            double parcenerg1 = dTdLmed * (coefdxP * dpdx - cinetico - (hidro) + (fontemassL + fontemassG) - latente - state.cells[i + 1].potBT / dxmed) / coefdxT;
            // parcela de energia relacionada aa troca termica
            double parcenerg2 = dTdLmed * (fluxcal) / coefdxT;

            // e feita uma avaliacao se a troca termica esta se dando de maneira muito rapida
            // como avanco da temperatura e explicita, isto pode levar a instabilidade no calculo termico
            // se for verificado que a troca termica esta ocorrendo de maneira rapida, o avanco e
            // feito em um numero maior de passos de uma celula para outra
            int npasso;
            double dxpasso;
            double limiteEstab = fabs(coefdxT) / dTdLmed;
            if (dxmed / (state.cells[i + 1].calor.resGlob + resanul) < (limiteEstab + 0. * 1000.)) {
                npasso = 1;
                dxpasso = dxmed;
            } else { // resistencia termica e pequena, se determinara em quantos passos se dara
                // o avanco de temperatrura
                npasso = (dxmed / (state.cells[i + 1].calor.resGlob + resanul)) / (limiteEstab + 0 * 1000) + 1;
                dxpasso = dxmed / npasso; // celula divida em npassos
            }
            double temppasso = state.cells[i + 1].temp;

            temppasso = dxpasso * (-(-state.cells[i + 1].temp) / dxpasso + parcenerg1 + parcenerg2); // primeiro avanco
            for (int j = 1; j < npasso; j++) {                                                  // avancos seguintes
                // as pacelas de energia sao mantidas, com excessao do fluxo de calor que e
                // reavalkiado a cada passo:
                state.cells[i + 1].calor.Tint = temppasso;
                state.cells[i + 1].calor.Vint = ugsmed + ulsmed;
                condliq = (1. - betmed) * state.cells[i + 1].flui.CondLiq(pmed, temppasso) + betmed * state.cells[i + 1].fluicol.CondLiq(pmed, temppasso);
                state.cells[i + 1].calor.kint = condliq * (1 - alfmed) + state.cells[i + 1].flui.CondGas(pmed, temppasso) * alfmed;
                state.cells[i + 1].calor.cpint = cpl * (1 - alfmed) + cpg * alfmed;
                state.cells[i + 1].calor.rhoint = rhol * (1 - alfmed) + rhog * alfmed;
                viscliq = (1. - betmed) * state.cells[i + 1].flui.ViscOleo(pmed, temppasso) + betmed * state.cells[i + 1].fluicol.VisFlu(pmed, temppasso);
                state.cells[i + 1].calor.viscint = viscliq * (1 - alfmed) * 1.e-3 + state.cells[i + 1].flui.ViscGas(pmed, temppasso) * alfmed * 1.e-3;
                if (state.steadyIteration != 0 && state.input.lingas == 1 && (i + 1 <= state.annulusTubingStart && i + 1 >= state.annulusTubingEnd)) {
                    int k = state.annulusTubingStart + state.tubingAnnulusStart - (i + 1);
                    double dtext = (state.gasCells[k].temp - state.gasCells[k - 1].temp) / npasso;
                    state.cells[i + 1].calor.Textern1 = state.gasCells[k].temp + (j - 1) * dtext;
                }
                fluxcal = sinalJ * state.cells[i + 1].calor.transperm(resanul);
                state.cells[i + 1].fluxcalmed += fluxcal;
                parcenerg2 = dTdLmed * (fluxcal) / coefdxT;
                temppasso = dxpasso * (-(-temppasso) / dxpasso + parcenerg1 + parcenerg2);
            }

            state.cells[i].temp = temppasso;
            state.cells[i + 1].fluxcalmed /= npasso;
            state.cells[i + 1].fluxcalmed = 1. * state.cells[i + 1].fluxcalmed;

            if (state.cells[i].temp < -50.)
                state.cells[i].temp = -50.;
            if (state.cells[i].temp > 200.)
                state.cells[i].temp = 200.;
        } else
            state.cells[i].temp = state.cells[i].calor.Textern1;
    } else { // caso em que a velocidade da mistura e muito baixa
        state.cells[i].temp = state.cells[i + 1].calor.Textern1;
        if (state.input.lingas == 1 && (i + 1 <= state.annulusTubingStart && i + 1 >= state.annulusTubingEnd)) {
            int j = state.annulusTubingStart + state.tubingAnnulusStart - (i + 1);
            state.cells[i + 1].calor.Vextern1 = 100.;
            state.cells[i + 1].calor.kextern1 = state.gasCells[j].calor.kint;
            state.cells[i + 1].calor.cpextern1 = state.gasCells[j].calor.cpint;
            state.cells[i + 1].calor.rhoextern1 = state.gasCells[j].calor.rhoint;
            state.cells[i + 1].calor.viscextern1 = state.gasCells[j].calor.viscint;
            double fluxcal = state.cells[i + 1].calor.transperm(0);
        }
    }
}

void computeGasTemperature(const ThermalState &state, int i, double tempantiga, int modoPerm) {

    if (state.thermalSourceDisabled == 0) {
        double dx = state.gasCells[i].dx0;
        double dxmed = 0.5 * (state.gasCells[i].dx0 + state.gasCells[i - 1].dx0);
        double area = state.gasCells[i].duto.area;
        double ugsmed;
        if (i < state.gasCellCount)
            ugsmed = state.gasCells[i].VGasR / state.gasCells[i].u1L;
        else {
            ugsmed = state.gasCells[i].VGasL / state.gasCells[i].u1L;
        }
        double rhog = state.gasCells[i].rg;
        double cpg = state.gasCells[i].flui.CalorGas(state.gasCells[i].presini, state.gasCells[i].tempini);
        double cvg = state.gasCells[i].flui.CalorGasVolMod(state.gasCells[i].presini, state.gasCells[i].tempini);
        double jtg = state.gasCells[i].flui.JTG(state.gasCells[i].presini, state.gasCells[i].tempini);
        double hidro = (rhog * ugsmed) * area * 9.82 * sin(state.gasCells[i].duto.teta);

        state.gasCells[i].calor.Tint = state.gasCells[i].tempini;
        state.gasCells[i].calor.dtL = state.gasCells[i].tempini - state.gasCells[i - 1].tempini;
        state.gasCells[i].calor.Vint = ugsmed;
        state.gasCells[i].calor.dt = state.gasCells[i].dt;
        state.gasCells[i].calor.kint = state.gasCells[i].flui.CondGas(state.gasCells[i].presini, state.gasCells[i].tempini);
        state.gasCells[i].calor.cpint = cpg;
        state.gasCells[i].calor.rhoint = rhog;
        state.gasCells[i].calor.viscint = state.gasCells[i].flui.ViscGas(state.gasCells[i].presini, state.gasCells[i].tempini) * 1.e-3;
        double dtemp = state.gasCells[i].temp * 0.01;
        if (fabs(state.gasCells[i].temp) < 1e-15)
            dtemp = 0.1;
        double rhogdT = state.gasCells[i].flui.MasEspGas(state.gasCells[i].presini, state.gasCells[i].tempini + dtemp) - rhog;
        state.gasCells[i].calor.betint = -(1 / state.gasCells[i].calor.rhoint) * rhogdT / (dtemp);
        if (modoPerm == 0)
            state.gasCells[i].fluxcal = state.gasCells[i].calor.transtrans();
        else
            state.gasCells[i].fluxcal = state.gasCells[i].calor.transperm();
        if (i >= state.tubingAnnulusStart && i <= state.tubingAnnulusEnd && state.networkCoupled == 1) {
            int kconecte = i - state.tubingAnnulusStart;
            int iconecte = state.annulusTubingStart - kconecte;
            state.gasCells[i].fluxcal -= state.cells[iconecte].calor.fluxFim;
        }

        double razdx;
        if (i < state.gasCellCount)
            razdx = dx / (dx + state.gasCells[i + 1].dx0);
        else
            razdx = dx / (dx + state.gasCells[i - 1].dx0);
        double coefTempo = rhog * cvg * area;
        double coefPresTempo = state.gasCells[i].flui.CalorGasPresMod(state.gasCells[i].presini, state.gasCells[i].tempini, state.gasCells[i].rg) *
                               (rhog * area);

        double coefdxT = rhog * ugsmed * cpg * area;
        double coefdxP = rhog * ugsmed * jtg * area;
        double dpdx;
        if (i < state.gasCellCount)
            dpdx = 2. * (((1 - razdx) * state.gasCells[i + 1].presini + razdx * state.gasCells[i].presini) - state.gasCells[i].presini) * 98066.5 / dx;
        else
            dpdx = 2. * (state.gasCells[i].presini - ((1 - razdx) * state.gasCells[i - 1].presini + razdx * state.gasCells[i].presini)) * 98066.5 / dx;
        double dtdx = (state.gasCells[i].tempini - state.gasCells[i - 1].tempini) / dxmed;
        if (i < state.gasCellCount)
            if (ugsmed < 0)
                dtdx = (state.gasCells[i + 1].tempini - state.gasCells[i].tempini) / dxmed;
        if ((i == 1 && ugsmed <= 0) || (i == state.gasCellCount && ugsmed <= 0))
            dtdx = 0.;

        double cinetico;
        double deljmix = 0.;
        double rhomix = rhog;
        double ugsmed0 = state.gasCells[i].VGasL / state.gasCells[i - 1].u1L;
        deljmix = (ugsmed - ugsmed0) / dx;

        cinetico = rhomix * area * ugsmed * ugsmed * deljmix;

        double fontemassG = 0.;
        double fontemassL = 0.;

        double fator = 1.;
        if ((*state.globals).lixo5 < 1000.)
            fator = 1.;

        state.gasCells[i].temp = ((coefTempo / state.gasCells[i].dt) * state.gasCells[i].temp - (fator) * (coefPresTempo * (state.gasCells[i].pres - state.gasCells[i].presini) * 98066.5 / state.gasCells[i].dt) + state.gasCells[i].dTdLCor * (-coefdxT * dtdx + coefdxP * dpdx - cinetico - hidro + fontemassL + fontemassG + state.gasCells[i].fluxcal)) / (coefTempo / state.gasCells[i].dt);

        if (state.gasCells[i].temp < -50.)
            state.gasCells[i].temp = -50.;
        if (state.gasCells[i].temp > 200.)
            state.gasCells[i].temp = 200.;

        if (i > 0)
            state.gasCells[i - 1].tempR = state.gasCells[i].temp;
        if (i < state.gasCellCount)
            state.gasCells[i + 1].tempL = state.gasCells[i].temp;
    } else {
        state.gasCells[i].temp = state.gasCells[i].calor.Textern1;

        if (i > 0)
            state.gasCells[i - 1].tempR = state.gasCells[i].temp;
        if (i < state.gasCellCount)
            state.gasCells[i + 1].tempL = state.gasCells[i].temp;
    }
}

void computeDischargeTemperature(const ThermalState &state, int i) {
    double dx0 = 0.5 * state.gasCells[i].dxL;
    double dx1 = 0.5 * state.gasCells[i].dx0;
    double RgasR = 0.;
    if (state.gasCells[i].razInter <= 0.5)
        RgasR = 2 * state.gasCells[i].razInter;
    double RgasL = 0.;
    if (state.gasCells[i - 1].razInter >= 0.5)
        RgasL = 2 * (state.gasCells[i - 1].razInter - 0.5);
    double LGasL = dx0 * RgasL;
    double LGasR = dx1 * RgasR;
    double LLiqL = dx0 - LGasL;
    double LLiqR = dx1 - LGasR;
    double LTotal = LLiqL + LLiqR + LGasL + LGasR;
    double dia = state.gasCells[i - 1].duto.a;
    double area = 0.25 * M_PI * dia * dia;
    double pres;
    double temp;

    pres = state.gasCells[i - 1].pres;
    temp = state.gasCells[i - 1].temp;
    double rho = ((LLiqL + LLiqR) * state.gasCells[i].MasEspFlu(pres, temp) + (LGasL + LGasR) * state.gasCells[i].flui.MasEspGas(pres, temp)) / LTotal;

    double vel1 = state.gasCells[i].VGasL / (state.gasCells[i].MasEspFlu(pres, temp) * state.gasCells[i - 1].duto.area);
    if (state.gasCells[i].razInter > (*state.globals).localtiny)
        vel1 = state.gasCells[i].VGasL / (state.gasCells[i].flui.MasEspGas(pres, temp) * state.gasCells[i - 1].duto.area);
    double cpl = ((LLiqL + LLiqR) * state.gasCells[i].CalorLiq(pres, temp) + (LGasL + LGasR) * state.gasCells[i].flui.CalorGas(pres, temp)) / LTotal;

    state.gasCells[i - 1].calor.Tint = temp;
    state.gasCells[i - 1].calor.Vint = vel1;
    double condliq = ((LLiqL + LLiqR) * state.gasCells[i - 1].CondLiq(pres, temp) + (LGasL + LGasR) * state.gasCells[i].flui.CondGas(pres, temp)) / LTotal;
    state.gasCells[i - 1].calor.kint = condliq;
    state.gasCells[i - 1].calor.cpint = cpl;
    state.cells[i - 1].calor.rhoint = rho;
    double viscliq = (((LLiqL + LLiqR) * state.gasCells[i].VisFlu(pres, temp) + (LGasL + LGasR) * state.gasCells[i].flui.ViscGas(pres, temp)) * 1e-3) / LTotal;
    state.cells[i - 1].calor.viscint = viscliq;

    double fluxcal = state.cells[i - 1].calor.transtrans();
    if ((i - 1) >= state.tubingAnnulusStart && (i - 1) <= state.tubingAnnulusEnd && state.networkCoupled == 1) {
        int kconecte = (i - 1) - state.tubingAnnulusStart;
        int iconecte = state.annulusTubingStart - kconecte;
        fluxcal -= state.cells[iconecte].calor.fluxFim;
    }

    state.gasCells[i - 1].tempR = state.gasCells[i].temp;
    state.gasCells[i].tempL = state.gasCells[i - 1].temp;
}

double computeGasLiftDischargeTemperature(const ThermalState &state, int igl) {
    int passo = floor((state.gasLiftChokes[igl].presEstag - state.gasLiftChokes[igl].presGarg) / 50) + 1;
    double deltaP = -(state.gasLiftChokes[igl].presEstag - state.gasLiftChokes[igl].presGarg) / passo;
    double PD0 = state.gasLiftChokes[igl].presEstag;
    double PD = state.gasLiftChokes[igl].presEstag + deltaP;
    double TD0 = state.gasLiftChokes[igl].tempEstag;
    double T1;
    for (int i = 0; i < passo; i++) {
        double cpg = state.gasLiftChokes[igl].flui.CalorGas(PD0, TD0);
        double DZDT = state.gasLiftChokes[igl].flui.DZDT(PD0, TD0);
        double tK = TD0 + 273.23;
        T1 = 1.0 / (1.0 / tK - ((286.998 / state.gasLiftChokes[igl].flui.Deng) * DZDT / cpg) * log((PD) / (PD0)));
        PD0 = PD;
        PD = PD - deltaP;
        TD0 = T1 - 273.23;
    }
    if (passo == 0) {
        if ((*state.globals).lixo5 > 1000) {
            int para;
            para = 0;
        }
        double cpg = state.gasLiftChokes[igl].flui.CalorGas(PD0, TD0);
        double jtg = state.gasLiftChokes[igl].flui.JTG(PD0, TD0) / cpg;
        TD0 -= jtg * (state.gasLiftChokes[igl].presEstag - state.gasLiftChokes[igl].presGarg) * 98066.52;
    }
    return TD0;
}

void updateProductionTemperaturePeriphery(const ThermalState &state, int i) {
    if (i > 0)
        state.cells[i - 1].tempR = state.cells[i].temp;
    if (i < state.lastCell)
        state.cells[i + 1].tempL = state.cells[i].temp;
    state.cells[i].tempini = state.cells[i].temp;
}

void computeOutletTemperature(const ThermalState &state) {

    if (state.input.chokep.abertura[0] <= 0.6 && state.input.chokep.abertura[0] > (*state.globals).localtiny && state.outletPressure < state.gasSurfacePressure) {
        double masentrada = state.cells[state.lastCell - 1].MR;
        double massgas = state.cells[state.lastCell - 1].MR - state.cells[state.lastCell - 1].MliqiniR;
        double rholp = state.cells[state.lastCell].flui.MasEspLiq(state.cells[state.lastCell].pres, state.cells[state.lastCell].temp);
        double rholc = state.cells[state.lastCell].fluicol.MasEspFlu(state.cells[state.lastCell].pres, state.cells[state.lastCell].temp);
        double betEF = state.cells[state.lastCell].bet;
        double tit = fabs(massgas / masentrada);

        double jtlM = (1. - betEF) * state.cells[state.lastCell].flui.JTL(state.cells[state.lastCell].pres, state.cells[state.lastCell].temp) - betEF / rholc; // alteraacao2
        double jtgM = state.cells[state.lastCell].flui.JTG(state.cells[state.lastCell].pres, state.cells[state.lastCell].temp);
        state.surfaceTemperature = state.cells[state.lastCell].temp + ((1. - tit) * jtlM + tit * jtgM) * (state.cells[state.lastCell].pres - state.cells[state.lastCell].pres); //????????
                                                                                                                  //???????????????????????????????celula[ncel].pres - celula[ncel].pres????????????????????????????????????????????

    } else
        state.surfaceTemperature = state.cells[state.lastCell - 1].temp;
}

}  // namespace sisprod::thermal
