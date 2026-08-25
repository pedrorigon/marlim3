/*
 * SisProdTrendOutput.cpp
 *
 * Trend-file output for the production and service lines, moved out of
 * SisProd.cpp. See SisProdTrendOutput.h for why the writers take their state
 * through TrendState instead of reading it from SProd.
 *
 * The bodies below were moved token for token. The five pre-existing anomalies
 * they carry -- the service cross-section writer reading the production
 * counters, the TENDTRANSG file named TENDTRANSP on the branch path, the
 * asymmetric round() in the TENDP name, the cross-section copy that ignores the
 * base offset, and the three stray endl -- are preserved deliberately, not
 * overlooked. They are catalogued in
 * specs/001-refatoracao-sisprod/evidencia/trend-diff.md.
 */
#include "SisProdTrendOutput.h"

#include "Leitura.h"
#include "Matriz.h"
#include "OutputI18n.h"
#include "variaveisGlobais1D.h"

#include <fstream>
#include <math.h>
#include <sstream>
#include <string>

namespace trendoutput {

void writeProductionTrendHeader(const TrendState &state, int i, int nrede) {
    if (state.inputData.ntendp > 0) {
        const auto t = [&state](const char *pt, const char *en) {
            return output_i18n::tr(state.inputData.idiomaSaida, pt, en);
        };
        ostringstream saidaT;
        if (state.branchIndex < 0 && state.inputData.AP == 0) {
            saidaT << pathPrefixoArqSaida << "TENDP" << "-" << round(state.inputData.trendp[i].comp) << ".dat";
        } else if (state.branchIndex < 0 && state.inputData.AP == 1) {
            saidaT << pathPrefixoArqSaida << "TENDP-AP-" << round(state.inputData.trendp[i].comp) << ".dat";
        } else {
            saidaT << pathPrefixoArqSaida << "Tramo" << state.branchIndex << "-R-" << nrede << "-" << "TENDP" << "-" << state.inputData.trendp[i].comp << ".dat";
        }
        string tmp = saidaT.str();
        ofstream escreveTrend(tmp.c_str(), ios_base::out);
        if (state.printPassCount == 1) {
            double comprimento = 0;
            int posicn = state.inputData.trendp[i].posic;
            for (int k = 0; k <= posicn; k++)
                comprimento += state.inputData.celp[k].dx;
            escreveTrend << t("# Comprimento a partir do Fundo de Poco (m) = ", "# Length from Bottomhole (m) = ") << comprimento << endl;
            escreveTrend << t("# Rotulo = ", "# Label = ") << state.inputData.trendp[i].rotulo << endl;
            escreveTrend << t("# Indice da Celula = ", "# Cell index = ") << state.inputData.trendp[i].posic << endl;
            if (state.branchIndex < 0 && state.inputData.AP == 1)
                escreveTrend << t(" Sequencia AP ;", " SA sequence ;");
            escreveTrend << t(" Tempo (s) ;", " Time (s) ;");
            if (state.inputData.trendp[i].pres == 1)
                escreveTrend << t(" Pressao (kgf/cm2) ;", " Pressure (kgf/cm2) ;");
            if (state.inputData.trendp[i].temp == 1)
                escreveTrend << t(" Temperatura (C) ;", " Temperature (C) ;");
            if (state.inputData.trendp[i].hol == 1)
                escreveTrend << t(" Holdup de liquido (-) ;", " Liquid holdup (-) ;");
            if (state.inputData.trendp[i].FVH == 1)
                escreveTrend << t(" Fracao Volumetrica Hidrato (-) ;", " Hydrate volumetric fraction (-) ;"); // solver de Hidratos - chris
            if (state.inputData.trendp[i].bet == 1)
                escreveTrend << t(" Fracao vol. de liquido complementar (-) ;", " Complementary liquid vol. fraction (-) ;");
            if (state.inputData.trendp[i].ugs == 1)
                escreveTrend << t(" Velocidade superficial do gas (m/s) ;", " Gas superficial velocity (m/s) ;");
            if (state.inputData.trendp[i].uls == 1)
                escreveTrend << t(" Velocidade superficial do liquido (m/s) ;", " Liquid superficial velocity (m/s) ;");
            if (state.inputData.trendp[i].ug == 1)
                escreveTrend << t(" Velocidade do gas (m/s) ;", " Gas velocity (m/s) ;");
            if (state.inputData.trendp[i].ul == 1)
                escreveTrend << t(" Velocidade do liquido (m/s) ;", " Liquid velocity (m/s) ;");
            if (state.inputData.trendp[i].arra == 1)
                escreveTrend << t(" Indicador de arranjo de fases (-) ;", " Phase pattern indicator (-) ;");
            if (state.inputData.trendp[i].viscl == 1)
                escreveTrend << t(" Viscosidade do Liquido (cP) ;", " Liquid viscosity (cP) ;");
            if (state.inputData.trendp[i].viscg == 1)
                escreveTrend << t(" Viscosidade do Gas (cP) ;", " Gas viscosity (cP) ;");
            if (state.inputData.trendp[i].rhog == 1)
                escreveTrend << t(" Massa Especifica do Gas (kg/m3) ;", " Gas density (kg/m3) ;");
            if (state.inputData.trendp[i].rhol == 1)
                escreveTrend << t(" Massa Especifica do Liquido (kg/m3) ;", " Liquid density (kg/m3) ;");
            if (state.inputData.trendp[i].rhoMix == 1)
                escreveTrend << t(" Massa Especifica da Mistura (kg/m3) ;", " Mixture density (kg/m3) ;");
            if (state.inputData.trendp[i].masg == 1)
                escreveTrend << t(" Vazao Massica do Gas (kg/s) ;", " Gas mass flow rate (kg/s) ;");
            if (state.inputData.trendp[i].masl == 1)
                escreveTrend << t(" Vazao Massica do Liquido (kg/s) ;", " Liquid mass flow rate (kg/s) ;");
            if (state.inputData.trendp[i].c0 == 1)
                escreveTrend << t(" Coeficiente de distribuição: C0 (-) ;", " Distribution coefficient: C0 (-) ;");
            if (state.inputData.trendp[i].ud == 1)
                escreveTrend << t(" Velocidade de escorregamento: Ud (m/s) ;", " Slip velocity: Ud (m/s) ;");
            if (state.inputData.trendp[i].RGO == 1)
                escreveTrend << " RGO (Sm3/Sm3) ;";
            if (state.inputData.trendp[i].deng == 1)
                escreveTrend << t(" Densidade do Gas (-) ;", " Gas specific gravity (-) ;");
            if (state.inputData.trendp[i].yco2 == 1)
                escreveTrend << t(" Fracao Molar de CO2 (-) ;", " CO2 molar fraction (-) ;");
            if (state.inputData.trendp[i].calor == 1)
                escreveTrend << t(" Fluxo de calor entre escoamento e parede (W/m) ;", " Heat flow between flow and wall (W/m) ;");
            if (state.inputData.trendp[i].masstrans == 1)
                escreveTrend << t(" Transferencia de Massa entre Fases (kg / [s m]) ;", " Interphase mass transfer (kg / [s m]) ;");
            if (state.inputData.trendp[i].qlst == 1)
                escreveTrend << t(" Vazao volumetrica standard de oleo morto (Sm3/d) ;", " Standard dead oil volumetric flow rate (Sm3/d) ;");
            if (state.inputData.trendp[i].qlwst == 1)
                escreveTrend << t(" Vazao volumetrica standard de oleo morto + agua (Sm3/d) ;", " Standard dead oil + water volumetric flow rate (Sm3/d) ;");
            if (state.inputData.trendp[i].qlstTot == 1)
                escreveTrend << t(" Vazao volumetrica standard de oleo morto + agua + liquido complementar (Sm3/d) ;", " Standard dead oil + water + complementary liquid volumetric flow rate (Sm3/d) ;");
            if (state.inputData.trendp[i].qgst == 1)
                escreveTrend << t(" Vazao volumetrica standard de gas livre + dissolvido (Sm3/d) ;", " Standard free + dissolved gas volumetric flow rate (Sm3/d) ;");
            if (state.inputData.trendp[i].api == 1)
                escreveTrend << t(" Grau API (-) ;", " API gravity (-) ;");
            if (state.inputData.trendp[i].bsw == 1)
                escreveTrend << " BSW (-) ;";
            if (state.inputData.trendp[i].hidro == 1)
                escreveTrend << t(" Termo Hidrostatico (Pa/m) ;", " Hydrostatic term (Pa/m) ;");
            if (state.inputData.trendp[i].fric == 1)
                escreveTrend << t(" Termo Friccao (Pa/m) ;", " Friction term (Pa/m) ;");
            if (state.inputData.trendp[i].dengD == 1)
                escreveTrend << t(" Densidade Gas Dissolvido In Situ (-) ;", " In-situ dissolved gas specific gravity (-) ;");
            if (state.inputData.trendp[i].dengL == 1)
                escreveTrend << t(" Densidade Gas Livre In Situ (-) ;", " In-situ free gas specific gravity (-) ;");
            if (state.inputData.trendp[i].mlFonte == 1)
                escreveTrend << t(" Fonte massica - Liq. (Hidrocarb+Agua) (kg/s);", " Mass source - Liq. (Hydrocarbon+Water) (kg/s);");
            if (state.inputData.trendp[i].mgFonte == 1)
                escreveTrend << t(" Fonte massica - Gas (kg/s);", " Mass source - Gas (kg/s);");
            if (state.inputData.trendp[i].mcFonte == 1)
                escreveTrend << t(" Fonte massica - Liq. Complementar (kg/s);", " Mass source - Complementary liquid (kg/s);");
            if (state.inputData.trendp[i].dpB == 1)
                escreveTrend << t(" Incremento de pressao de Bombeio (kgf/cm2);", " Pump pressure increment (kgf/cm2);");
            if (state.inputData.trendp[i].potB == 1)
                escreveTrend << t(" Potencia de Bombeio (kW);", " Pump power (kW);");
            if (state.inputData.trendp[i].tempChokeJus == 1)
                escreveTrend << t(" Temperatura a Jusante do Choke de Superficie (C);", " Surface choke downstream temperature (C);");
            if (state.inputData.trendp[i].reyi == 1)
                escreveTrend << t(" Reynolds interno da mistura (-) ;", " Internal mixture Reynolds (-) ;");
            if (state.inputData.trendp[i].reye == 1)
                escreveTrend << t(" Reynolds externo (-) ;", " External Reynolds (-) ;");
            if (state.inputData.trendp[i].Fr == 1)
                escreveTrend << t(" Froud (-) ;", " Froude (-) ;");
            if (state.inputData.trendp[i].grashi == 1)
                escreveTrend << t(" Grashof interno da mistura (-) ;", " Internal mixture Grashof (-) ;");
            if (state.inputData.trendp[i].grashe == 1)
                escreveTrend << t(" Grashof externo (-) ;", " External Grashof (-) ;");
            if (state.inputData.trendp[i].nusi == 1)
                escreveTrend << t(" Nusselt interno da mistura (-) ;", " Internal mixture Nusselt (-) ;");
            if (state.inputData.trendp[i].nuse == 1)
                escreveTrend << t(" Nusselt externo (-) ;", " External Nusselt (-) ;");
            if (state.inputData.trendp[i].hi == 1)
                escreveTrend << t(" Coeficiente de pelicula interno da mistura (W/(m2.K)) ;", " Internal mixture film coefficient (W/(m2.K)) ;");
            if (state.inputData.trendp[i].he == 1)
                escreveTrend << t(" Coeficiente de pelicula externo (W/(m2.K)) ;", " External film coefficient (W/(m2.K)) ;");
            if (state.inputData.trendp[i].pri == 1)
                escreveTrend << t(" Prandtl interno da mistura (-) ;", " Internal mixture Prandtl (-) ;");
            if (state.inputData.trendp[i].pre == 1)
                escreveTrend << t(" Prandtl externo (-) ;", " External Prandtl (-) ;");
            if (state.inputData.trendp[i].Rs == 1)
                escreveTrend << t(" Razao de Solubilidade (-) ;", " Solubility ratio (-) ;");
            if (state.inputData.trendp[i].Bo == 1)
                escreveTrend << t(" Fator Volume de Formacao (-) ;", " Formation volume factor (-) ;");
            if (state.inputData.trendp[i].volMonM1PT == 1)
                escreveTrend << t(" Volume de liquido a montante da Master1, a PT, m3 ;", " Liquid volume upstream of Master1, at PT, m3 ;");
            if (state.inputData.trendp[i].volJusM1PT == 1)
                escreveTrend << t(" Volume de liquido a jusante da Master1, a PT, m3 ;", " Liquid volume downstream of Master1, at PT, m3 ;");
            if (state.inputData.trendp[i].volMonM1ST == 1)
                escreveTrend << t(" Volume de liquido a montante da Master1, standard, m3 ;", " Liquid volume upstream of Master1, standard, m3 ;");
            if (state.inputData.trendp[i].volJusM1ST == 1)
                escreveTrend << t(" Volume de liquido a jusante da Master1, standard, m3 ;", " Liquid volume downstream of Master1, standard, m3 ;");
            if (state.inputData.trendp[i].inventarioGas == 1)
                escreveTrend << t(" Inventario de Gas em toda a tubulação, standard, m3 ;", " Gas inventory in whole tubing, standard, m3 ;");
            if (state.inputData.trendp[i].inventarioLiq == 1)
                escreveTrend << t(" Inventario de Liquido em toda a tubulação, standard, m3 ;", " Liquid inventory in whole tubing, standard, m3 ;");
            if (state.inputData.trendp[i].diamInt == 1)
                escreveTrend << t(" Diametro Interno da tubulacao, m ;", " Tubing inner diameter, m ;");
            if (state.inputData.trendp[i].TempParede == 1)
                escreveTrend << t(" Temperatura Interna da Parede, C ;", " Internal wall temperature, C ;");
            if (state.inputData.trendp[i].subResfria == 1)
                escreveTrend << t(" Subresfriamento, C ;", " Subcooling, C ;");
            if (state.inputData.trendp[i].dadosParafina == 1) {
                escreveTrend << t(" TIAC (C) C;", " TIAC (C) C;");
                escreveTrend << t(" Cp Parafina (J/[kg C]) C;", " Paraffin Cp (J/[kg C]) C;");
                escreveTrend << t(" Condutividade Termica Parafina (W / [m K]) C;", " Paraffin thermal conductivity (W / [m K]) C;");
                escreveTrend << t(" Massa Especifica Parafina (kg/m3) C;", " Paraffin density (kg/m3) C;");
                escreveTrend << t(" Massa molar do Liquido Parafina (kg/mol) C;", " Paraffin liquid molar mass (kg/mol) C;");
                escreveTrend << t(" Difusividade Massica Parafina (m2/s) C;", " Paraffin mass diffusivity (m2/s) C;");
                escreveTrend << t(" Fluxo Massico de Parafina Total (kg/(m2-s)) C;", " Total paraffin mass flux (kg/(m2-s)) C;");
                escreveTrend << t(" Fluxo Massico de Parafina por Difusao (kg/(m2-s)) C;", " Paraffin diffusive mass flux (kg/(m2-s)) C;");
                escreveTrend << t(" Vazao Massica de Parafina por Difusao (kg/(s)) C;", " Paraffin mass rate (kg/(s)) C;");
                escreveTrend << t(" Gradiente de concentracao de parafina (1/m) C;", " Paraffin concentration gradient (1/m) C;");
                escreveTrend << t(" Condutividade do deposito (W/(m-K)) C;", " Deposit conductivity (W/(m-K)) C;");
                escreveTrend << t(" Temperatura da Interface do deposito (C) C;", " Deposit interface temperature (C) C;");
            }
            if (state.inputData.trendp[i].autoVal == 1) {
                for (int konta1 = 0; konta1 < 3; konta1++) {
                    escreveTrend << t(" Celeridade, familia de onda ", " Celerity, wave family ") << konta1 << t(" m/s ;", " m/s ;");
                }
            }
            if (state.inputData.trendp[i].autoVel == 1) {
                for (int konta1 = 0; konta1 < 3; konta1++) {
                    for (int konta2 = 0; konta2 < 3; konta2++)
                        escreveTrend << t(" Componente do autovetor, condicao adiabatica, familia de onda = ", " Eigenvector component, adiabatic condition, wave family = ") << konta1 << t("termo = ", " term = ") << konta2 << " ;";
                }
            }
            if (state.inputData.trendp[i].flutuacao == 1) {
                for (int konta1 = 0; konta1 < 3; konta1++) {
                    escreveTrend << t(" Componente de flutuacao da familia de onda ", " Fluctuation component of wave family ") << konta1 << " ;";
                }
            }
            escreveTrend << endl;
        }
        escreveTrend.close();
        // caso nao seja simulacao POCO_INJETOR
        if (state.inputData.tipoSimulacao != tipoSimulacao_t::poco_injetor) {
            arqRelatorioPerfis << tmp.c_str() << endl;
            arqRelatorioPerfis.flush();
        }
    }
}

void writeProductionTrendRows(const TrendState &state, int i, int nrede) {
    if (state.inputData.ntendp > 0) {
        int size = state.productionCount[i] - state.productionCountBase[i] + 1;
        FullMtx<double> saidatrend(size, state.inputData.nvartrendp[i] + 1);
        for (int k = 0; k < size; k++)
            for (int j = 0; j <= state.inputData.nvartrendp[i]; j++)
                saidatrend[k][j] = state.productionBuffer[i][k + state.productionCountBase[i]][j];
        ostringstream saidaT;
        if (state.branchIndex < 0 && state.inputData.AP == 0) {
            saidaT << pathPrefixoArqSaida << "TENDP" << "-" << state.inputData.trendp[i].comp << ".dat";
        } else if (state.branchIndex < 0 && state.inputData.AP == 1) {
            saidaT << pathPrefixoArqSaida << "TENDP-AP-" << round(state.inputData.trendp[i].comp) << ".dat";
        } else {
            saidaT << pathPrefixoArqSaida << "Tramo" << state.branchIndex << "-R-" << nrede << "-" << "TENDP" << "-" << state.inputData.trendp[i].comp << ".dat";
        }
        string tmp = saidaT.str();
        ofstream escreveTrend(tmp.c_str(), ios_base::app);
        int nc = saidatrend.col();
        int nl = saidatrend.lin();
        if (state.branchIndex < 0 && state.inputData.AP == 1)
            nc++;
        for (int i = 0; i < nl; i++) {
            if (saidatrend[i][0] <= -9999)
                break;
            for (int j = 0; j < nc; j++) {
                escreveTrend.width(20);
                escreveTrend.precision(19);
                if (state.branchIndex < 0 && state.inputData.AP == 1) {
                    if (j == 0)
                        escreveTrend << (*state.globals).sequenciaAP << " ; ";
                    else
                        escreveTrend << saidatrend[i][j - 1] << " ; ";
                } else
                    escreveTrend << saidatrend[i][j] << " ; ";
            }
            escreveTrend << endl;
        }
        escreveTrend.close();
        // caso nao seja simulacao POCO_INJETOR
        if (state.inputData.tipoSimulacao != tipoSimulacao_t::poco_injetor) {
            arqRelatorioPerfis << tmp.c_str() << endl;
            arqRelatorioPerfis.flush();
        }
    }
}

void writeServiceTrendHeader(const TrendState &state, int i, int nrede) {
    if (state.inputData.ntendg > 0 && state.inputData.lingas > 0) {
        const auto t = [&state](const char *pt, const char *en) {
            return output_i18n::tr(state.inputData.idiomaSaida, pt, en);
        };
        ostringstream saidaT;
        if (state.branchIndex < 0 && state.inputData.AP == 0)
            saidaT << pathPrefixoArqSaida << "TENDG" << "-" << state.inputData.trendg[i].comp << ".dat";
        else if (state.branchIndex < 0 && state.inputData.AP == 1) {
            saidaT << pathPrefixoArqSaida << "TENDG-AP-" << state.inputData.trendg[i].comp << ".dat";
        } else
            saidaT << pathPrefixoArqSaida << "Tramo" << state.branchIndex << "-R-" << nrede << "-" << "TENDG" << "-" << state.inputData.trendg[i].comp << ".dat";
        string tmp = saidaT.str();
        ofstream escreveTrend(tmp.c_str(), ios_base::out);
        if (state.printPassCount == 1) {
            double comprimento = 0;
            int posicn = state.inputData.trendg[i].posic;
            for (int k = 0; k <= posicn; k++)
                comprimento += state.inputData.celg[k].dx;
            escreveTrend << t("# Comprimento a partir da Plataforma (m) = ", "# Length from Platform (m) = ") << comprimento << endl;
            escreveTrend << t("# Rotulo = ", "# Label = ") << state.inputData.trendg[i].rotulo << endl;
            escreveTrend << t("# Indice da Celula = ", "# Cell index = ") << state.inputData.trendg[i].posic << endl;
            if (state.branchIndex < 0 && state.inputData.AP == 1)
                escreveTrend << t(" Sequencia AP ;", " SA sequence ;");
            escreveTrend << t(" Tempo (s) ;", " Time (s) ;");
            if (state.inputData.trendg[i].pres == 1)
                escreveTrend << t(" Pressao (kgf/cm2) ;", " Pressure (kgf/cm2) ;");
            if (state.inputData.trendg[i].temp == 1)
                escreveTrend << t(" Temperatura (C) ;", " Temperature (C) ;");
            if (state.inputData.trendg[i].ugs == 1)
                escreveTrend << t(" Velocidade superficial do gas (m/s) ;", " Gas superficial velocity (m/s) ;");
            if (state.inputData.trendg[i].ug == 1)
                escreveTrend << t(" Velocidade do gas (m/s) ;", " Gas velocity (m/s) ;");
            if (state.inputData.trendg[i].tens == 1)
                escreveTrend << t(" Tensao Cisalhante (N/m2) ;", " Shear stress (N/m2) ;");
            if (state.inputData.trendg[i].viscg == 1)
                escreveTrend << t(" Viscosidade do Gas (cP) ;", " Gas viscosity (cP) ;");
            if (state.inputData.trendg[i].rhog == 1)
                escreveTrend << t(" Massa Especifica do Gas (kg/m3) ;", " Gas density (kg/m3) ;");
            if (state.inputData.trendg[i].masg == 1)
                escreveTrend << t(" Vazao Massica do Gas (kg/s) ;", " Gas mass flow rate (kg/s) ;");
            if (state.inputData.trendg[i].hidro == 1)
                escreveTrend << t(" Termo Hidrostatico (Pa/m) ;", " Hydrostatic term (Pa/m) ;");
            if (state.inputData.trendg[i].FVHG == 1)
                escreveTrend << t(" FVHG (-) ;", " FVHG (-) ;");
            if (state.inputData.trendg[i].fric == 1)
                escreveTrend << t(" Termo Friccao (Pa/m) ;", " Friction term (Pa/m) ;");
            if (state.inputData.trendg[i].calor == 1)
                escreveTrend << t(" Fluxo de calor entre escoamento e parede (W/m) ;", " Heat flow between flow and wall (W/m) ;");
            if (state.inputData.trendg[i].qgst == 1)
                escreveTrend << t(" Vazao volumetrica standard de Gas (Sm3/d) ;", " Standard gas volumetric flow rate (Sm3/d) ;");
            if (state.inputData.trendg[i].pEstagVGL == 1)
                escreveTrend << t(" Pressao de Estagnacao VGL (kgf/cm²) ;", " VGL stagnation pressure (kgf/cm2) ;");
            if (state.inputData.trendg[i].tEstagVGL == 1)
                escreveTrend << t(" Temperatura de Estagnacao VGL (C) ;", " VGL stagnation temperature (C) ;");
            if (state.inputData.trendg[i].pGargVGL == 1)
                escreveTrend << t(" Pressao na Garganta VGL (kgf/cm²) ;", " VGL throat pressure (kgf/cm2) ;");
            if (state.inputData.trendg[i].tGargVGL == 1)
                escreveTrend << t(" Temperatura na Garganta VGL (C) ;", " VGL throat temperature (C) ;");
            if (state.inputData.trendg[i].velgarg == 1)
                escreveTrend << t(" Velocidade na VGL (m/s) ;", " Velocity in VGL (m/s) ;");
            if (state.inputData.trendg[i].qVGL == 1)
                escreveTrend << t(" Vazao volumetrica na VGL (mÂ³/s) ;", " Volumetric flow rate in VGL (m3/s) ;");
            if (state.inputData.trendg[i].reyi == 1)
                escreveTrend << t(" Reynolds interno (-) ;", " Internal Reynolds (-) ;");
            if (state.inputData.trendg[i].reye == 1)
                escreveTrend << t(" Reynolds externo (-) ;", " External Reynolds (-) ;");
            if (state.inputData.trendg[i].grashi == 1)
                escreveTrend << t(" Grashof interno (-) ;", " Internal Grashof (-) ;");
            if (state.inputData.trendg[i].grashe == 1)
                escreveTrend << t(" Grashof externo (-) ;", " External Grashof (-) ;");
            if (state.inputData.trendg[i].nusi == 1)
                escreveTrend << t(" Nusselt interno (-) ;", " Internal Nusselt (-) ;");
            if (state.inputData.trendg[i].nuse == 1)
                escreveTrend << t(" Nusselt externo (-) ;", " External Nusselt (-) ;");
            if (state.inputData.trendg[i].hi == 1)
                escreveTrend << t(" Coeficiente de pelicula interno (W/[m2 K]) ;", " Internal film coefficient (W/[m2 K]) ;");
            if (state.inputData.trendg[i].he == 1)
                escreveTrend << t(" Coeficiente de pelicula externo (W/[m2 K]) ;", " External film coefficient (W/[m2 K]) ;");
            if (state.inputData.trendg[i].pri == 1)
                escreveTrend << t(" Prandtl interno (-) ;", " Internal Prandtl (-) ;");
            if (state.inputData.trendg[i].pre == 1)
                escreveTrend << t(" Prandtl externo (-) ;", " External Prandtl (-) ;");
            if (state.inputData.trendg[i].diamInt == 1)
                escreveTrend << t(" Diametro Interno da tubulacao, m ;", " Tubing inner diameter, m ;");
            if (state.inputData.trendg[i].TempParede == 1)
                escreveTrend << t(" Temperatura Interna da Parede, C ;", " Internal wall temperature, C ;");
            if (state.inputData.trendg[i].subResfria == 1)
                escreveTrend << t(" Subresfriamento, C ;", " Subcooling, C ;");
            escreveTrend << endl;
        }
        escreveTrend << endl;
        escreveTrend.close();
        // caso nao seja simulacao POCO_INJETOR
        if (state.inputData.tipoSimulacao != tipoSimulacao_t::poco_injetor) {
            arqRelatorioPerfis << tmp.c_str() << endl;
            arqRelatorioPerfis.flush();
        }
    }
}

void writeServiceTrendRows(const TrendState &state, int i, int nrede) {
    if (state.inputData.ntendg > 0 && state.inputData.lingas > 0) {
        int size = state.serviceCount[i] - state.serviceCountBase[i] + 1;
        FullMtx<double> saidatrend(size, state.inputData.nvartrendg[i] + 1);
        for (int k = 0; k < size; k++)
            for (int j = 0; j <= state.inputData.nvartrendg[i]; j++)
                saidatrend[k][j] = state.serviceBuffer[i][k + state.serviceCountBase[i]][j];
        ostringstream saidaT;
        if (state.branchIndex < 0 && state.inputData.AP == 0)
            saidaT << pathPrefixoArqSaida << "TENDG" << "-" << state.inputData.trendg[i].comp << ".dat";
        else if (state.branchIndex < 0 && state.inputData.AP == 1) {
            saidaT << pathPrefixoArqSaida << "TENDG-AP-" << state.inputData.trendg[i].comp << ".dat";
        } else
            saidaT << pathPrefixoArqSaida << "Tramo" << state.branchIndex << "-R-" << nrede << "-" << "TENDG" << "-" << state.inputData.trendg[i].comp << ".dat";
        string tmp = saidaT.str();
        ofstream escreveTrend(tmp.c_str(), ios_base::app);
        int nc = saidatrend.col();
        int nl = saidatrend.lin();
        if (state.branchIndex < 0 && state.inputData.AP == 1)
            nc++;
        for (int i = 0; i < nl; i++) {
            if (saidatrend[i][0] <= -9999)
                break;
            for (int j = 0; j < nc; j++) {
                escreveTrend.width(20);
                escreveTrend.precision(19);
                if (state.branchIndex < 0 && state.inputData.AP == 1) {
                    if (j == 0)
                        escreveTrend << (*state.globals).sequenciaAP << " ; ";
                    else
                        escreveTrend << saidatrend[i][j - 1] << " ; ";
                } else
                    escreveTrend << saidatrend[i][j] << " ; ";
            }
            escreveTrend << endl;
        }
        escreveTrend.close();
        // caso nao seja simulacao POCO_INJETOR
        if (state.inputData.tipoSimulacao != tipoSimulacao_t::poco_injetor) {
            arqRelatorioPerfis << tmp.c_str() << endl;
            arqRelatorioPerfis.flush();
        }
    }
}

void writeProductionCrossSectionTrendHeader(const TrendState &state, int i) {

    if (state.inputData.ntendtransp > 0) {
        const auto t = [&state](const char *pt, const char *en) {
            return output_i18n::tr(state.inputData.idiomaSaida, pt, en);
        };
        int poscel = state.inputData.trendtransp[i].comp;
        int poscam = state.inputData.trendtransp[i].camada - 1;
        int posdiscre = state.inputData.trendtransp[i].discre - 1;
        ostringstream saidaT;
        if (state.branchIndex < 0)
            saidaT << pathPrefixoArqSaida << "TENDTRANSP" << "-" << poscel << "-" << poscam << "-" << posdiscre << ".dat";
        else
            saidaT << pathPrefixoArqSaida << "Tramo" << state.branchIndex << "-" << "TENDTRANSP" << "-" << poscel << "-" << poscam << "-" << posdiscre << ".dat";
        string tmp = saidaT.str();
        ofstream escreveTrend(tmp.c_str(), ios_base::out);
        if (state.printPassCount == 1) {
            escreveTrend << t("# Rotulo = ", "# Label = ") << state.inputData.trendtransp[i].rotulo << endl;
            escreveTrend << t(" Tempo (s) ; ", " Time (s) ; ");
            escreveTrend << t(" Temperatura (C) ;", " Temperature (C) ;");
            escreveTrend << endl;
        }
        escreveTrend.close();
        // caso nao seja simulacao POCO_INJETOR
        if (state.inputData.tipoSimulacao != tipoSimulacao_t::poco_injetor) {
            arqRelatorioPerfis << tmp.c_str() << endl;
            arqRelatorioPerfis.flush();
        }
    }
}

void writeProductionCrossSectionTrendRows(const TrendState &state, int i) {
    if (state.inputData.ntendtransp > 0) {
        int size = state.crossSectionCount[i] - state.crossSectionCountBase[i] + 1;
        FullMtx<double> saidatrend(size, 2);
        for (int k = 0; k < size; k++)
            for (int j = 0; j < 2; j++)
                saidatrend[k][j] = state.productionCrossSectionBuffer[i][k][j];
        int poscel = state.inputData.trendtransp[i].comp;
        int poscam = state.inputData.trendtransp[i].camada - 1;
        int posdiscre = state.inputData.trendtransp[i].discre - 1;
        ostringstream saidaT;
        if (state.branchIndex < 0)
            saidaT << pathPrefixoArqSaida << "TENDTRANSP" << "-" << poscel << "-" << poscam << "-" << posdiscre << ".dat";
        else
            saidaT << pathPrefixoArqSaida << "Tramo" << state.branchIndex << "-" << "TENDTRANSP" << "-" << poscel << "-" << poscam << "-" << posdiscre << ".dat";
        string tmp = saidaT.str();
        ofstream escreveTrend(tmp.c_str(), ios_base::app);
        int nc = saidatrend.col();
        int nl = saidatrend.lin();
        for (int i = 0; i < nl; i++) {
            if (saidatrend[i][0] <= -9999)
                break;
            for (int j = 0; j < nc; j++) {
                escreveTrend.width(20);
                escreveTrend.precision(19);
                escreveTrend << saidatrend[i][j] << " ; ";
            }
            escreveTrend << endl;
        }
        escreveTrend << endl;
        escreveTrend.close();
        // caso nao seja simulacao POCO_INJETOR
        if (state.inputData.tipoSimulacao != tipoSimulacao_t::poco_injetor) {
            arqRelatorioPerfis << tmp.c_str() << endl;
            arqRelatorioPerfis.flush();
        }
    }
}

void writeServiceCrossSectionTrendHeader(const TrendState &state, int i) {
    if (state.inputData.ntendtransg > 0 && state.inputData.lingas > 0) {
        const auto t = [&state](const char *pt, const char *en) {
            return output_i18n::tr(state.inputData.idiomaSaida, pt, en);
        };
        int poscel = state.inputData.trendtransg[i].comp;
        int poscam = state.inputData.trendtransg[i].camada - 1;
        int posdiscre = state.inputData.trendtransg[i].discre - 1;
        ostringstream saidaT;
        if (state.branchIndex < 0)
            saidaT << pathPrefixoArqSaida << "TENDTRANSG" << "-" << poscel << "-" << poscam << "-" << posdiscre << ".dat";
        else
            saidaT << pathPrefixoArqSaida << "Tramo" << state.branchIndex << "-" << "TENDTRANSP" << "-" << poscel << "-" << poscam << "-" << posdiscre << ".dat";
        string tmp = saidaT.str();
        ofstream escreveTrend(tmp.c_str(), ios_base::out);
        if (state.printPassCount == 1) {
            escreveTrend << t("# Rotulo = ", "# Label = ") << state.inputData.trendtransg[i].rotulo << endl;
            escreveTrend << t(" Tempo (s) ; ", " Time (s) ; ");
            escreveTrend << t(" Temperatura (C) ;", " Temperature (C) ;");
            escreveTrend << endl;
        }
        escreveTrend.close();
        // caso nao seja simulacao POCO_INJETOR
        if (state.inputData.tipoSimulacao != tipoSimulacao_t::poco_injetor) {
            arqRelatorioPerfis << tmp.c_str() << endl;
            arqRelatorioPerfis.flush();
        }
    }
}

void writeServiceCrossSectionTrendRows(const TrendState &state, int i) {
    if (state.inputData.ntendtransg > 0 && state.inputData.lingas > 0) {
        int size = state.crossSectionCount[i] - state.crossSectionCountBase[i] + 1;
        FullMtx<double> saidatrend(size, 2);
        for (int k = 0; k < size; k++)
            for (int j = 0; j < 2; j++)
                saidatrend[k][j] = state.serviceCrossSectionBuffer[i][k][j];
        int poscel = state.inputData.trendtransg[i].comp;
        int poscam = state.inputData.trendtransg[i].camada - 1;
        int posdiscre = state.inputData.trendtransg[i].discre - 1;
        ostringstream saidaT;
        if (state.branchIndex < 0)
            saidaT << pathPrefixoArqSaida << "TENDTRANSG" << "-" << poscel << "-" << poscam << "-" << posdiscre << ".dat";
        else
            saidaT << pathPrefixoArqSaida << "Tramo" << state.branchIndex << "-" << "TENDTRANSP" << "-" << poscel << "-" << poscam << "-" << posdiscre << ".dat";
        string tmp = saidaT.str();
        ofstream escreveTrend(tmp.c_str(), ios_base::app);
        int nc = saidatrend.col();
        int nl = saidatrend.lin();
        for (int i = 0; i < nl; i++) {
            if (saidatrend[i][0] <= -9999)
                break;
            for (int j = 0; j < nc; j++) {
                escreveTrend.width(20);
                escreveTrend.precision(19);
                escreveTrend << saidatrend[i][j] << " ; ";
            }
            escreveTrend << endl;
        }
        escreveTrend << endl;

        escreveTrend.close();
        // caso nao seja simulacao POCO_INJETOR
        if (state.inputData.tipoSimulacao != tipoSimulacao_t::poco_injetor) {
            arqRelatorioPerfis << tmp.c_str() << endl;
            arqRelatorioPerfis.flush();
        }
    }
}
} // namespace trendoutput
