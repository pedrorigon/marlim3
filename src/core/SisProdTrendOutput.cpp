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
namespace {

/// Records a produced file in the profile report.
///
/// Identical in all eight writers, including the comment, which is why it is
/// the one piece shared without a hook.
void reportProducedFile(const TrendState &state, const string &fileName) {
    // caso nao seja simulacao POCO_INJETOR
    if (state.inputData.tipoSimulacao != tipoSimulacao_t::poco_injetor) {
        arqRelatorioPerfis << fileName.c_str() << endl;
        arqRelatorioPerfis.flush();
    }
}

/// The skeleton the four header writers share.
///
/// buildName and writeCaptions are the hooks. They are template parameters, not
/// function pointers or virtuals: the compiler inlines them, so the skeleton
/// costs nothing over the four copies it replaces (FR-022).
///
/// blankLineBeforeClose is a value, not a mode flag. Only the service header
/// emits that line, and it emits it OUTSIDE the print-pass guard -- so a service
/// trend file starts with a blank line even on the passes that write no header.
/// That asymmetry is in the baseline and is preserved here (A4-05).
template <typename BuildName, typename WriteCaptions>
void writeTrendHeaderFile(const TrendState &state, BuildName buildName,
                          WriteCaptions writeCaptions, bool blankLineBeforeClose) {
    string fileName = buildName();
    ofstream trendFile(fileName.c_str(), ios_base::out);
    if (state.printPassCount == 1) {
        writeCaptions(trendFile);
        trendFile << endl;
    }
    if (blankLineBeforeClose)
        trendFile << endl;
    trendFile.close();
    reportProducedFile(state, fileName);
}

/// The skeleton the four row writers share.
///
/// windowSize, columnCount, series and sourceBase are resolved by the caller,
/// which is what keeps the four data sources -- and the two anomalies they
/// carry -- visible at the call site instead of buried in a condition here.
///
/// The AP sequence column exists only for the production and service lines; the
/// two cross-section writers pass supportsApSequence = false, which makes the
/// guard permanently false exactly as in the baseline, where the condition is
/// absent. The condition itself is evaluated once rather than per column: it
/// reads only branchIndex and AP, and neither changes inside the loop.
template <typename BuildName>
void writeTrendRowsFile(const TrendState &state, int windowSize, int columnCount,
                        double **series, int sourceBase, BuildName buildName,
                        bool supportsApSequence, bool blankLineBeforeClose) {
    FullMtx<double> window(windowSize, columnCount);
    for (int sourceRow = 0; sourceRow < windowSize; sourceRow++)
        for (int sourceColumn = 0; sourceColumn < columnCount; sourceColumn++)
            window[sourceRow][sourceColumn] = series[sourceRow + sourceBase][sourceColumn];
    string fileName = buildName();
    ofstream trendFile(fileName.c_str(), ios_base::app);
    int columnsToWrite = window.col();
    int rowsToWrite = window.lin();
    const bool apSequenceColumn =
        supportsApSequence && state.branchIndex < 0 && state.inputData.AP == 1;
    if (apSequenceColumn)
        columnsToWrite++;
    for (int rowIndex = 0; rowIndex < rowsToWrite; rowIndex++) {
        if (window[rowIndex][0] <= -9999)
            break;
        for (int columnIndex = 0; columnIndex < columnsToWrite; columnIndex++) {
            trendFile.width(20);
            trendFile.precision(19);
            if (apSequenceColumn) {
                if (columnIndex == 0)
                    trendFile << (*state.globals).sequenciaAP << " ; ";
                else
                    trendFile << window[rowIndex][columnIndex - 1] << " ; ";
            } else
                trendFile << window[rowIndex][columnIndex] << " ; ";
        }
        trendFile << endl;
    }
    if (blankLineBeforeClose)
        trendFile << endl;
    trendFile.close();
    reportProducedFile(state, fileName);
}

} // namespace

/// Writes the column captions of a production-line trend file.
void writeProductionTrendHeader(const TrendState &state, int trendIndex, int networkIndex) {
    if (state.inputData.ntendp > 0) {
        const auto translate = [&state](const char *ptBrText, const char *enText) {
            return output_i18n::tr(state.inputData.idiomaSaida, ptBrText, enText);
        };
        writeTrendHeaderFile(
            state,
            [&] {
                ostringstream fileNameStream;
                if (state.branchIndex < 0 && state.inputData.AP == 0) {
                    fileNameStream << pathPrefixoArqSaida << "TENDP" << "-" << round(state.inputData.trendp[trendIndex].comp) << ".dat";
                } else if (state.branchIndex < 0 && state.inputData.AP == 1) {
                    fileNameStream << pathPrefixoArqSaida << "TENDP-AP-" << round(state.inputData.trendp[trendIndex].comp) << ".dat";
                } else {
                    fileNameStream << pathPrefixoArqSaida << "Tramo" << state.branchIndex << "-R-" << networkIndex << "-" << "TENDP" << "-" << state.inputData.trendp[trendIndex].comp << ".dat";
                }
                return fileNameStream.str();
            },
            [&](ofstream &trendFile) {
                    double lengthFromOrigin = 0;
                    int lastCellIndex = state.inputData.trendp[trendIndex].posic;
                    for (int cell = 0; cell <= lastCellIndex; cell++)
                        lengthFromOrigin += state.inputData.celp[cell].dx;
                    trendFile << translate("# Comprimento a partir do Fundo de Poco (m) = ", "# Length from Bottomhole (m) = ") << lengthFromOrigin << endl;
                    trendFile << translate("# Rotulo = ", "# Label = ") << state.inputData.trendp[trendIndex].rotulo << endl;
                    trendFile << translate("# Indice da Celula = ", "# Cell index = ") << state.inputData.trendp[trendIndex].posic << endl;
                    if (state.branchIndex < 0 && state.inputData.AP == 1)
                        trendFile << translate(" Sequencia AP ;", " SA sequence ;");
                    trendFile << translate(" Tempo (s) ;", " Time (s) ;");
                    if (state.inputData.trendp[trendIndex].pres == 1)
                        trendFile << translate(" Pressao (kgf/cm2) ;", " Pressure (kgf/cm2) ;");
                    if (state.inputData.trendp[trendIndex].temp == 1)
                        trendFile << translate(" Temperatura (C) ;", " Temperature (C) ;");
                    if (state.inputData.trendp[trendIndex].hol == 1)
                        trendFile << translate(" Holdup de liquido (-) ;", " Liquid holdup (-) ;");
                    if (state.inputData.trendp[trendIndex].FVH == 1)
                        trendFile << translate(" Fracao Volumetrica Hidrato (-) ;", " Hydrate volumetric fraction (-) ;"); // solver de Hidratos - chris
                    if (state.inputData.trendp[trendIndex].bet == 1)
                        trendFile << translate(" Fracao vol. de liquido complementar (-) ;", " Complementary liquid vol. fraction (-) ;");
                    if (state.inputData.trendp[trendIndex].ugs == 1)
                        trendFile << translate(" Velocidade superficial do gas (m/s) ;", " Gas superficial velocity (m/s) ;");
                    if (state.inputData.trendp[trendIndex].uls == 1)
                        trendFile << translate(" Velocidade superficial do liquido (m/s) ;", " Liquid superficial velocity (m/s) ;");
                    if (state.inputData.trendp[trendIndex].ug == 1)
                        trendFile << translate(" Velocidade do gas (m/s) ;", " Gas velocity (m/s) ;");
                    if (state.inputData.trendp[trendIndex].ul == 1)
                        trendFile << translate(" Velocidade do liquido (m/s) ;", " Liquid velocity (m/s) ;");
                    if (state.inputData.trendp[trendIndex].arra == 1)
                        trendFile << translate(" Indicador de arranjo de fases (-) ;", " Phase pattern indicator (-) ;");
                    if (state.inputData.trendp[trendIndex].viscl == 1)
                        trendFile << translate(" Viscosidade do Liquido (cP) ;", " Liquid viscosity (cP) ;");
                    if (state.inputData.trendp[trendIndex].viscg == 1)
                        trendFile << translate(" Viscosidade do Gas (cP) ;", " Gas viscosity (cP) ;");
                    if (state.inputData.trendp[trendIndex].rhog == 1)
                        trendFile << translate(" Massa Especifica do Gas (kg/m3) ;", " Gas density (kg/m3) ;");
                    if (state.inputData.trendp[trendIndex].rhol == 1)
                        trendFile << translate(" Massa Especifica do Liquido (kg/m3) ;", " Liquid density (kg/m3) ;");
                    if (state.inputData.trendp[trendIndex].rhoMix == 1)
                        trendFile << translate(" Massa Especifica da Mistura (kg/m3) ;", " Mixture density (kg/m3) ;");
                    if (state.inputData.trendp[trendIndex].masg == 1)
                        trendFile << translate(" Vazao Massica do Gas (kg/s) ;", " Gas mass flow rate (kg/s) ;");
                    if (state.inputData.trendp[trendIndex].masl == 1)
                        trendFile << translate(" Vazao Massica do Liquido (kg/s) ;", " Liquid mass flow rate (kg/s) ;");
                    if (state.inputData.trendp[trendIndex].c0 == 1)
                        trendFile << translate(" Coeficiente de distribuição: C0 (-) ;", " Distribution coefficient: C0 (-) ;");
                    if (state.inputData.trendp[trendIndex].ud == 1)
                        trendFile << translate(" Velocidade de escorregamento: Ud (m/s) ;", " Slip velocity: Ud (m/s) ;");
                    if (state.inputData.trendp[trendIndex].RGO == 1)
                        trendFile << " RGO (Sm3/Sm3) ;";
                    if (state.inputData.trendp[trendIndex].deng == 1)
                        trendFile << translate(" Densidade do Gas (-) ;", " Gas specific gravity (-) ;");
                    if (state.inputData.trendp[trendIndex].yco2 == 1)
                        trendFile << translate(" Fracao Molar de CO2 (-) ;", " CO2 molar fraction (-) ;");
                    if (state.inputData.trendp[trendIndex].calor == 1)
                        trendFile << translate(" Fluxo de calor entre escoamento e parede (W/m) ;", " Heat flow between flow and wall (W/m) ;");
                    if (state.inputData.trendp[trendIndex].masstrans == 1)
                        trendFile << translate(" Transferencia de Massa entre Fases (kg / [s m]) ;", " Interphase mass transfer (kg / [s m]) ;");
                    if (state.inputData.trendp[trendIndex].qlst == 1)
                        trendFile << translate(" Vazao volumetrica standard de oleo morto (Sm3/d) ;", " Standard dead oil volumetric flow rate (Sm3/d) ;");
                    if (state.inputData.trendp[trendIndex].qlwst == 1)
                        trendFile << translate(" Vazao volumetrica standard de oleo morto + agua (Sm3/d) ;", " Standard dead oil + water volumetric flow rate (Sm3/d) ;");
                    if (state.inputData.trendp[trendIndex].qlstTot == 1)
                        trendFile << translate(" Vazao volumetrica standard de oleo morto + agua + liquido complementar (Sm3/d) ;", " Standard dead oil + water + complementary liquid volumetric flow rate (Sm3/d) ;");
                    if (state.inputData.trendp[trendIndex].qgst == 1)
                        trendFile << translate(" Vazao volumetrica standard de gas livre + dissolvido (Sm3/d) ;", " Standard free + dissolved gas volumetric flow rate (Sm3/d) ;");
                    if (state.inputData.trendp[trendIndex].api == 1)
                        trendFile << translate(" Grau API (-) ;", " API gravity (-) ;");
                    if (state.inputData.trendp[trendIndex].bsw == 1)
                        trendFile << " BSW (-) ;";
                    if (state.inputData.trendp[trendIndex].hidro == 1)
                        trendFile << translate(" Termo Hidrostatico (Pa/m) ;", " Hydrostatic term (Pa/m) ;");
                    if (state.inputData.trendp[trendIndex].fric == 1)
                        trendFile << translate(" Termo Friccao (Pa/m) ;", " Friction term (Pa/m) ;");
                    if (state.inputData.trendp[trendIndex].dengD == 1)
                        trendFile << translate(" Densidade Gas Dissolvido In Situ (-) ;", " In-situ dissolved gas specific gravity (-) ;");
                    if (state.inputData.trendp[trendIndex].dengL == 1)
                        trendFile << translate(" Densidade Gas Livre In Situ (-) ;", " In-situ free gas specific gravity (-) ;");
                    if (state.inputData.trendp[trendIndex].mlFonte == 1)
                        trendFile << translate(" Fonte massica - Liq. (Hidrocarb+Agua) (kg/s);", " Mass source - Liq. (Hydrocarbon+Water) (kg/s);");
                    if (state.inputData.trendp[trendIndex].mgFonte == 1)
                        trendFile << translate(" Fonte massica - Gas (kg/s);", " Mass source - Gas (kg/s);");
                    if (state.inputData.trendp[trendIndex].mcFonte == 1)
                        trendFile << translate(" Fonte massica - Liq. Complementar (kg/s);", " Mass source - Complementary liquid (kg/s);");
                    if (state.inputData.trendp[trendIndex].dpB == 1)
                        trendFile << translate(" Incremento de pressao de Bombeio (kgf/cm2);", " Pump pressure increment (kgf/cm2);");
                    if (state.inputData.trendp[trendIndex].potB == 1)
                        trendFile << translate(" Potencia de Bombeio (kW);", " Pump power (kW);");
                    if (state.inputData.trendp[trendIndex].tempChokeJus == 1)
                        trendFile << translate(" Temperatura a Jusante do Choke de Superficie (C);", " Surface choke downstream temperature (C);");
                    if (state.inputData.trendp[trendIndex].reyi == 1)
                        trendFile << translate(" Reynolds interno da mistura (-) ;", " Internal mixture Reynolds (-) ;");
                    if (state.inputData.trendp[trendIndex].reye == 1)
                        trendFile << translate(" Reynolds externo (-) ;", " External Reynolds (-) ;");
                    if (state.inputData.trendp[trendIndex].Fr == 1)
                        trendFile << translate(" Froud (-) ;", " Froude (-) ;");
                    if (state.inputData.trendp[trendIndex].grashi == 1)
                        trendFile << translate(" Grashof interno da mistura (-) ;", " Internal mixture Grashof (-) ;");
                    if (state.inputData.trendp[trendIndex].grashe == 1)
                        trendFile << translate(" Grashof externo (-) ;", " External Grashof (-) ;");
                    if (state.inputData.trendp[trendIndex].nusi == 1)
                        trendFile << translate(" Nusselt interno da mistura (-) ;", " Internal mixture Nusselt (-) ;");
                    if (state.inputData.trendp[trendIndex].nuse == 1)
                        trendFile << translate(" Nusselt externo (-) ;", " External Nusselt (-) ;");
                    if (state.inputData.trendp[trendIndex].hi == 1)
                        trendFile << translate(" Coeficiente de pelicula interno da mistura (W/(m2.K)) ;", " Internal mixture film coefficient (W/(m2.K)) ;");
                    if (state.inputData.trendp[trendIndex].he == 1)
                        trendFile << translate(" Coeficiente de pelicula externo (W/(m2.K)) ;", " External film coefficient (W/(m2.K)) ;");
                    if (state.inputData.trendp[trendIndex].pri == 1)
                        trendFile << translate(" Prandtl interno da mistura (-) ;", " Internal mixture Prandtl (-) ;");
                    if (state.inputData.trendp[trendIndex].pre == 1)
                        trendFile << translate(" Prandtl externo (-) ;", " External Prandtl (-) ;");
                    if (state.inputData.trendp[trendIndex].Rs == 1)
                        trendFile << translate(" Razao de Solubilidade (-) ;", " Solubility ratio (-) ;");
                    if (state.inputData.trendp[trendIndex].Bo == 1)
                        trendFile << translate(" Fator Volume de Formacao (-) ;", " Formation volume factor (-) ;");
                    if (state.inputData.trendp[trendIndex].volMonM1PT == 1)
                        trendFile << translate(" Volume de liquido a montante da Master1, a PT, m3 ;", " Liquid volume upstream of Master1, at PT, m3 ;");
                    if (state.inputData.trendp[trendIndex].volJusM1PT == 1)
                        trendFile << translate(" Volume de liquido a jusante da Master1, a PT, m3 ;", " Liquid volume downstream of Master1, at PT, m3 ;");
                    if (state.inputData.trendp[trendIndex].volMonM1ST == 1)
                        trendFile << translate(" Volume de liquido a montante da Master1, standard, m3 ;", " Liquid volume upstream of Master1, standard, m3 ;");
                    if (state.inputData.trendp[trendIndex].volJusM1ST == 1)
                        trendFile << translate(" Volume de liquido a jusante da Master1, standard, m3 ;", " Liquid volume downstream of Master1, standard, m3 ;");
                    if (state.inputData.trendp[trendIndex].inventarioGas == 1)
                        trendFile << translate(" Inventario de Gas em toda a tubulação, standard, m3 ;", " Gas inventory in whole tubing, standard, m3 ;");
                    if (state.inputData.trendp[trendIndex].inventarioLiq == 1)
                        trendFile << translate(" Inventario de Liquido em toda a tubulação, standard, m3 ;", " Liquid inventory in whole tubing, standard, m3 ;");
                    if (state.inputData.trendp[trendIndex].diamInt == 1)
                        trendFile << translate(" Diametro Interno da tubulacao, m ;", " Tubing inner diameter, m ;");
                    if (state.inputData.trendp[trendIndex].TempParede == 1)
                        trendFile << translate(" Temperatura Interna da Parede, C ;", " Internal wall temperature, C ;");
                    if (state.inputData.trendp[trendIndex].subResfria == 1)
                        trendFile << translate(" Subresfriamento, C ;", " Subcooling, C ;");
                    if (state.inputData.trendp[trendIndex].dadosParafina == 1) {
                        trendFile << translate(" TIAC (C) C;", " TIAC (C) C;");
                        trendFile << translate(" Cp Parafina (J/[kg C]) C;", " Paraffin Cp (J/[kg C]) C;");
                        trendFile << translate(" Condutividade Termica Parafina (W / [m K]) C;", " Paraffin thermal conductivity (W / [m K]) C;");
                        trendFile << translate(" Massa Especifica Parafina (kg/m3) C;", " Paraffin density (kg/m3) C;");
                        trendFile << translate(" Massa molar do Liquido Parafina (kg/mol) C;", " Paraffin liquid molar mass (kg/mol) C;");
                        trendFile << translate(" Difusividade Massica Parafina (m2/s) C;", " Paraffin mass diffusivity (m2/s) C;");
                        trendFile << translate(" Fluxo Massico de Parafina Total (kg/(m2-s)) C;", " Total paraffin mass flux (kg/(m2-s)) C;");
                        trendFile << translate(" Fluxo Massico de Parafina por Difusao (kg/(m2-s)) C;", " Paraffin diffusive mass flux (kg/(m2-s)) C;");
                        trendFile << translate(" Vazao Massica de Parafina por Difusao (kg/(s)) C;", " Paraffin mass rate (kg/(s)) C;");
                        trendFile << translate(" Gradiente de concentracao de parafina (1/m) C;", " Paraffin concentration gradient (1/m) C;");
                        trendFile << translate(" Condutividade do deposito (W/(m-K)) C;", " Deposit conductivity (W/(m-K)) C;");
                        trendFile << translate(" Temperatura da Interface do deposito (C) C;", " Deposit interface temperature (C) C;");
                    }
                    if (state.inputData.trendp[trendIndex].autoVal == 1) {
                        for (int waveFamily = 0; waveFamily < 3; waveFamily++) {
                            trendFile << translate(" Celeridade, familia de onda ", " Celerity, wave family ") << waveFamily << translate(" m/s ;", " m/s ;");
                        }
                    }
                    if (state.inputData.trendp[trendIndex].autoVel == 1) {
                        for (int waveFamily = 0; waveFamily < 3; waveFamily++) {
                            for (int eigenvectorTerm = 0; eigenvectorTerm < 3; eigenvectorTerm++)
                                trendFile << translate(" Componente do autovetor, condicao adiabatica, familia de onda = ", " Eigenvector component, adiabatic condition, wave family = ") << waveFamily << translate("termo = ", " term = ") << eigenvectorTerm << " ;";
                        }
                    }
                    if (state.inputData.trendp[trendIndex].flutuacao == 1) {
                        for (int waveFamily = 0; waveFamily < 3; waveFamily++) {
                            trendFile << translate(" Componente de flutuacao da familia de onda ", " Fluctuation component of wave family ") << waveFamily << " ;";
                        }
                    }
            },
            /*blankLineBeforeClose=*/false);
    }
}

/// Appends the buffered rows of a production-line trend file.
void writeProductionTrendRows(const TrendState &state, int trendIndex, int networkIndex) {
    if (state.inputData.ntendp > 0) {
        writeTrendRowsFile(
            state,
            state.productionCount[trendIndex] - state.productionCountBase[trendIndex] + 1,
            state.inputData.nvartrendp[trendIndex] + 1,
            state.productionBuffer[trendIndex],
            state.productionCountBase[trendIndex],
            [&] {
                ostringstream fileNameStream;
                if (state.branchIndex < 0 && state.inputData.AP == 0) {
                    fileNameStream << pathPrefixoArqSaida << "TENDP" << "-" << state.inputData.trendp[trendIndex].comp << ".dat";
                } else if (state.branchIndex < 0 && state.inputData.AP == 1) {
                    fileNameStream << pathPrefixoArqSaida << "TENDP-AP-" << round(state.inputData.trendp[trendIndex].comp) << ".dat";
                } else {
                    fileNameStream << pathPrefixoArqSaida << "Tramo" << state.branchIndex << "-R-" << networkIndex << "-" << "TENDP" << "-" << state.inputData.trendp[trendIndex].comp << ".dat";
                }
                return fileNameStream.str();
            },
            /*supportsApSequence=*/true,
            /*blankLineBeforeClose=*/false);
    }
}

/// Writes the column captions of a service-line trend file.
void writeServiceTrendHeader(const TrendState &state, int trendIndex, int networkIndex) {
    if (state.inputData.ntendg > 0 && state.inputData.lingas > 0) {
        const auto translate = [&state](const char *ptBrText, const char *enText) {
            return output_i18n::tr(state.inputData.idiomaSaida, ptBrText, enText);
        };
        writeTrendHeaderFile(
            state,
            [&] {
                ostringstream fileNameStream;
                if (state.branchIndex < 0 && state.inputData.AP == 0)
                    fileNameStream << pathPrefixoArqSaida << "TENDG" << "-" << state.inputData.trendg[trendIndex].comp << ".dat";
                else if (state.branchIndex < 0 && state.inputData.AP == 1) {
                    fileNameStream << pathPrefixoArqSaida << "TENDG-AP-" << state.inputData.trendg[trendIndex].comp << ".dat";
                } else
                    fileNameStream << pathPrefixoArqSaida << "Tramo" << state.branchIndex << "-R-" << networkIndex << "-" << "TENDG" << "-" << state.inputData.trendg[trendIndex].comp << ".dat";
                return fileNameStream.str();
            },
            [&](ofstream &trendFile) {
                    double lengthFromOrigin = 0;
                    int lastCellIndex = state.inputData.trendg[trendIndex].posic;
                    for (int cell = 0; cell <= lastCellIndex; cell++)
                        lengthFromOrigin += state.inputData.celg[cell].dx;
                    trendFile << translate("# Comprimento a partir da Plataforma (m) = ", "# Length from Platform (m) = ") << lengthFromOrigin << endl;
                    trendFile << translate("# Rotulo = ", "# Label = ") << state.inputData.trendg[trendIndex].rotulo << endl;
                    trendFile << translate("# Indice da Celula = ", "# Cell index = ") << state.inputData.trendg[trendIndex].posic << endl;
                    if (state.branchIndex < 0 && state.inputData.AP == 1)
                        trendFile << translate(" Sequencia AP ;", " SA sequence ;");
                    trendFile << translate(" Tempo (s) ;", " Time (s) ;");
                    if (state.inputData.trendg[trendIndex].pres == 1)
                        trendFile << translate(" Pressao (kgf/cm2) ;", " Pressure (kgf/cm2) ;");
                    if (state.inputData.trendg[trendIndex].temp == 1)
                        trendFile << translate(" Temperatura (C) ;", " Temperature (C) ;");
                    if (state.inputData.trendg[trendIndex].ugs == 1)
                        trendFile << translate(" Velocidade superficial do gas (m/s) ;", " Gas superficial velocity (m/s) ;");
                    if (state.inputData.trendg[trendIndex].ug == 1)
                        trendFile << translate(" Velocidade do gas (m/s) ;", " Gas velocity (m/s) ;");
                    if (state.inputData.trendg[trendIndex].tens == 1)
                        trendFile << translate(" Tensao Cisalhante (N/m2) ;", " Shear stress (N/m2) ;");
                    if (state.inputData.trendg[trendIndex].viscg == 1)
                        trendFile << translate(" Viscosidade do Gas (cP) ;", " Gas viscosity (cP) ;");
                    if (state.inputData.trendg[trendIndex].rhog == 1)
                        trendFile << translate(" Massa Especifica do Gas (kg/m3) ;", " Gas density (kg/m3) ;");
                    if (state.inputData.trendg[trendIndex].masg == 1)
                        trendFile << translate(" Vazao Massica do Gas (kg/s) ;", " Gas mass flow rate (kg/s) ;");
                    if (state.inputData.trendg[trendIndex].hidro == 1)
                        trendFile << translate(" Termo Hidrostatico (Pa/m) ;", " Hydrostatic term (Pa/m) ;");
                    if (state.inputData.trendg[trendIndex].FVHG == 1)
                        trendFile << translate(" FVHG (-) ;", " FVHG (-) ;");
                    if (state.inputData.trendg[trendIndex].fric == 1)
                        trendFile << translate(" Termo Friccao (Pa/m) ;", " Friction term (Pa/m) ;");
                    if (state.inputData.trendg[trendIndex].calor == 1)
                        trendFile << translate(" Fluxo de calor entre escoamento e parede (W/m) ;", " Heat flow between flow and wall (W/m) ;");
                    if (state.inputData.trendg[trendIndex].qgst == 1)
                        trendFile << translate(" Vazao volumetrica standard de Gas (Sm3/d) ;", " Standard gas volumetric flow rate (Sm3/d) ;");
                    if (state.inputData.trendg[trendIndex].pEstagVGL == 1)
                        trendFile << translate(" Pressao de Estagnacao VGL (kgf/cm²) ;", " VGL stagnation pressure (kgf/cm2) ;");
                    if (state.inputData.trendg[trendIndex].tEstagVGL == 1)
                        trendFile << translate(" Temperatura de Estagnacao VGL (C) ;", " VGL stagnation temperature (C) ;");
                    if (state.inputData.trendg[trendIndex].pGargVGL == 1)
                        trendFile << translate(" Pressao na Garganta VGL (kgf/cm²) ;", " VGL throat pressure (kgf/cm2) ;");
                    if (state.inputData.trendg[trendIndex].tGargVGL == 1)
                        trendFile << translate(" Temperatura na Garganta VGL (C) ;", " VGL throat temperature (C) ;");
                    if (state.inputData.trendg[trendIndex].velgarg == 1)
                        trendFile << translate(" Velocidade na VGL (m/s) ;", " Velocity in VGL (m/s) ;");
                    if (state.inputData.trendg[trendIndex].qVGL == 1)
                        trendFile << translate(" Vazao volumetrica na VGL (mÂ³/s) ;", " Volumetric flow rate in VGL (m3/s) ;");
                    if (state.inputData.trendg[trendIndex].reyi == 1)
                        trendFile << translate(" Reynolds interno (-) ;", " Internal Reynolds (-) ;");
                    if (state.inputData.trendg[trendIndex].reye == 1)
                        trendFile << translate(" Reynolds externo (-) ;", " External Reynolds (-) ;");
                    if (state.inputData.trendg[trendIndex].grashi == 1)
                        trendFile << translate(" Grashof interno (-) ;", " Internal Grashof (-) ;");
                    if (state.inputData.trendg[trendIndex].grashe == 1)
                        trendFile << translate(" Grashof externo (-) ;", " External Grashof (-) ;");
                    if (state.inputData.trendg[trendIndex].nusi == 1)
                        trendFile << translate(" Nusselt interno (-) ;", " Internal Nusselt (-) ;");
                    if (state.inputData.trendg[trendIndex].nuse == 1)
                        trendFile << translate(" Nusselt externo (-) ;", " External Nusselt (-) ;");
                    if (state.inputData.trendg[trendIndex].hi == 1)
                        trendFile << translate(" Coeficiente de pelicula interno (W/[m2 K]) ;", " Internal film coefficient (W/[m2 K]) ;");
                    if (state.inputData.trendg[trendIndex].he == 1)
                        trendFile << translate(" Coeficiente de pelicula externo (W/[m2 K]) ;", " External film coefficient (W/[m2 K]) ;");
                    if (state.inputData.trendg[trendIndex].pri == 1)
                        trendFile << translate(" Prandtl interno (-) ;", " Internal Prandtl (-) ;");
                    if (state.inputData.trendg[trendIndex].pre == 1)
                        trendFile << translate(" Prandtl externo (-) ;", " External Prandtl (-) ;");
                    if (state.inputData.trendg[trendIndex].diamInt == 1)
                        trendFile << translate(" Diametro Interno da tubulacao, m ;", " Tubing inner diameter, m ;");
                    if (state.inputData.trendg[trendIndex].TempParede == 1)
                        trendFile << translate(" Temperatura Interna da Parede, C ;", " Internal wall temperature, C ;");
                    if (state.inputData.trendg[trendIndex].subResfria == 1)
                        trendFile << translate(" Subresfriamento, C ;", " Subcooling, C ;");
            },
            /*blankLineBeforeClose=*/true);
    }
}

/// Appends the buffered rows of a service-line trend file.
void writeServiceTrendRows(const TrendState &state, int trendIndex, int networkIndex) {
    if (state.inputData.ntendg > 0 && state.inputData.lingas > 0) {
        writeTrendRowsFile(
            state,
            state.serviceCount[trendIndex] - state.serviceCountBase[trendIndex] + 1,
            state.inputData.nvartrendg[trendIndex] + 1,
            state.serviceBuffer[trendIndex],
            state.serviceCountBase[trendIndex],
            [&] {
                ostringstream fileNameStream;
                if (state.branchIndex < 0 && state.inputData.AP == 0)
                    fileNameStream << pathPrefixoArqSaida << "TENDG" << "-" << state.inputData.trendg[trendIndex].comp << ".dat";
                else if (state.branchIndex < 0 && state.inputData.AP == 1) {
                    fileNameStream << pathPrefixoArqSaida << "TENDG-AP-" << state.inputData.trendg[trendIndex].comp << ".dat";
                } else
                    fileNameStream << pathPrefixoArqSaida << "Tramo" << state.branchIndex << "-R-" << networkIndex << "-" << "TENDG" << "-" << state.inputData.trendg[trendIndex].comp << ".dat";
                return fileNameStream.str();
            },
            /*supportsApSequence=*/true,
            /*blankLineBeforeClose=*/false);
    }
}

/// Writes the captions of a production cross-section temperature trend.
void writeProductionCrossSectionTrendHeader(const TrendState &state, int trendIndex) {
    if (state.inputData.ntendtransp > 0) {
        const auto translate = [&state](const char *ptBrText, const char *enText) {
            return output_i18n::tr(state.inputData.idiomaSaida, ptBrText, enText);
        };
        writeTrendHeaderFile(
            state,
            [&] {
                int cellPosition = state.inputData.trendtransp[trendIndex].comp;
                int layerIndex = state.inputData.trendtransp[trendIndex].camada - 1;
                int discretizationIndex = state.inputData.trendtransp[trendIndex].discre - 1;
                ostringstream fileNameStream;
                if (state.branchIndex < 0)
                    fileNameStream << pathPrefixoArqSaida << "TENDTRANSP" << "-" << cellPosition << "-" << layerIndex << "-" << discretizationIndex << ".dat";
                else
                    fileNameStream << pathPrefixoArqSaida << "Tramo" << state.branchIndex << "-" << "TENDTRANSP" << "-" << cellPosition << "-" << layerIndex << "-" << discretizationIndex << ".dat";
                return fileNameStream.str();
            },
            [&](ofstream &trendFile) {
                    trendFile << translate("# Rotulo = ", "# Label = ") << state.inputData.trendtransp[trendIndex].rotulo << endl;
                    trendFile << translate(" Tempo (s) ; ", " Time (s) ; ");
                    trendFile << translate(" Temperatura (C) ;", " Temperature (C) ;");
            },
            /*blankLineBeforeClose=*/false);
    }
}

/// Appends the buffered rows of a production cross-section trend.
///
/// sourceBase is 0, not crossSectionCountBase[i]: the baseline sizes the window
/// from the counters but always copies from the start of the buffer (A4-04).
/// Passing the zero explicitly is what keeps that visible.
void writeProductionCrossSectionTrendRows(const TrendState &state, int trendIndex) {
    if (state.inputData.ntendtransp > 0) {
        writeTrendRowsFile(
            state,
            state.crossSectionCount[trendIndex] - state.crossSectionCountBase[trendIndex] + 1,
            2,
            state.productionCrossSectionBuffer[trendIndex],
            0,
            [&] {
                int cellPosition = state.inputData.trendtransp[trendIndex].comp;
                int layerIndex = state.inputData.trendtransp[trendIndex].camada - 1;
                int discretizationIndex = state.inputData.trendtransp[trendIndex].discre - 1;
                ostringstream fileNameStream;
                if (state.branchIndex < 0)
                    fileNameStream << pathPrefixoArqSaida << "TENDTRANSP" << "-" << cellPosition << "-" << layerIndex << "-" << discretizationIndex << ".dat";
                else
                    fileNameStream << pathPrefixoArqSaida << "Tramo" << state.branchIndex << "-" << "TENDTRANSP" << "-" << cellPosition << "-" << layerIndex << "-" << discretizationIndex << ".dat";
                return fileNameStream.str();
            },
            /*supportsApSequence=*/false,
            /*blankLineBeforeClose=*/true);
    }
}

/// Writes the captions of a service cross-section temperature trend.
void writeServiceCrossSectionTrendHeader(const TrendState &state, int trendIndex) {
    if (state.inputData.ntendtransg > 0 && state.inputData.lingas > 0) {
        const auto translate = [&state](const char *ptBrText, const char *enText) {
            return output_i18n::tr(state.inputData.idiomaSaida, ptBrText, enText);
        };
        writeTrendHeaderFile(
            state,
            [&] {
                int cellPosition = state.inputData.trendtransg[trendIndex].comp;
                int layerIndex = state.inputData.trendtransg[trendIndex].camada - 1;
                int discretizationIndex = state.inputData.trendtransg[trendIndex].discre - 1;
                ostringstream fileNameStream;
                if (state.branchIndex < 0)
                    fileNameStream << pathPrefixoArqSaida << "TENDTRANSG" << "-" << cellPosition << "-" << layerIndex << "-" << discretizationIndex << ".dat";
                else
                    fileNameStream << pathPrefixoArqSaida << "Tramo" << state.branchIndex << "-" << "TENDTRANSP" << "-" << cellPosition << "-" << layerIndex << "-" << discretizationIndex << ".dat";
                return fileNameStream.str();
            },
            [&](ofstream &trendFile) {
                    trendFile << translate("# Rotulo = ", "# Label = ") << state.inputData.trendtransg[trendIndex].rotulo << endl;
                    trendFile << translate(" Tempo (s) ; ", " Time (s) ; ");
                    trendFile << translate(" Temperatura (C) ;", " Temperature (C) ;");
            },
            /*blankLineBeforeClose=*/false);
    }
}

/// Appends the buffered rows of a service cross-section trend.
///
/// Reads the PRODUCTION cross-section counters, and copies from the start of the
/// buffer. Both are baseline behaviour (A4-01, A4-04), preserved deliberately.
void writeServiceCrossSectionTrendRows(const TrendState &state, int trendIndex) {
    if (state.inputData.ntendtransg > 0 && state.inputData.lingas > 0) {
        writeTrendRowsFile(
            state,
            state.crossSectionCount[trendIndex] - state.crossSectionCountBase[trendIndex] + 1,
            2,
            state.serviceCrossSectionBuffer[trendIndex],
            0,
            [&] {
                int cellPosition = state.inputData.trendtransg[trendIndex].comp;
                int layerIndex = state.inputData.trendtransg[trendIndex].camada - 1;
                int discretizationIndex = state.inputData.trendtransg[trendIndex].discre - 1;
                ostringstream fileNameStream;
                if (state.branchIndex < 0)
                    fileNameStream << pathPrefixoArqSaida << "TENDTRANSG" << "-" << cellPosition << "-" << layerIndex << "-" << discretizationIndex << ".dat";
                else
                    fileNameStream << pathPrefixoArqSaida << "Tramo" << state.branchIndex << "-" << "TENDTRANSP" << "-" << cellPosition << "-" << layerIndex << "-" << discretizationIndex << ".dat";
                return fileNameStream.str();
            },
            /*supportsApSequence=*/false,
            /*blankLineBeforeClose=*/true);
    }
}

} // namespace trendoutput
