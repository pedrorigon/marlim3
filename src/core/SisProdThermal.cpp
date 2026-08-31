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

}  // namespace sisprod::thermal
