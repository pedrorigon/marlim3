/*
 * SisProdTrendOutput.cpp
 *
 * Trend-file output for the production and service lines, extracted from
 * SisProd.cpp. See SisProdTrendOutput.h for why the writers take their state
 * through TrendState instead of reading it from SProd.
 *
 * Layout, top to bottom: language selection, the caption tables, the file-name
 * builders, and the two skeletons the eight entry points share. Nothing above a
 * section is allowed to depend on anything below it.
 *
 * The caption tables are data, not code. Sixty-two production captions and
 * thirty-two service captions used to be sixty-two and thirty-two consecutive
 * if statements; expressing them as a table of (flag, pt-BR text, en text)
 * makes adding a column an edit to one line of data, and removes the only place
 * in this module where a long conditional chain was doing bookkeeping rather
 * than deciding anything.
 */
#include "SisProdTrendOutput.h"

#include "Leitura.h"
#include "OutputI18n.h"
#include "variaveisGlobais1D.h"

#include <cstddef>
#include <fstream>
#include <sstream>
#include <string>

namespace trendoutput {
namespace {

// ---------------------------------------------------------------- language --

/// Picks the caption spelling for the configured output language.
///
/// A callable rather than a free function so the caption blocks read as
/// translate(pt, en) instead of threading the language code through every one
/// of the ninety-odd call sites.
struct CaptionTranslator {
    int languageCode;

    const char *operator()(const char *ptBrText, const char *enText) const {
        return output_i18n::tr(languageCode, ptBrText, enText);
    }
};

// ---------------------------------------------------------------- captions --

/// A column written only when its flag is set in the trend descriptor.
template <typename Trend>
struct OptionalColumn {
    int Trend::*enabled;
    const char *ptBrText;
    const char *enText;
};

/// A column belonging to a group that is written or skipped as a whole.
struct FixedColumn {
    const char *ptBrText;
    const char *enText;
};

/// Order is the file format: these are emitted left to right, and the row
/// writers rely on the data columns arriving in exactly this sequence.
constexpr OptionalColumn<detTRENDP> kProductionColumns[] = {
    {&detTRENDP::pres, " Pressao (kgf/cm2) ;", " Pressure (kgf/cm2) ;"},
    {&detTRENDP::temp, " Temperatura (C) ;", " Temperature (C) ;"},
    {&detTRENDP::hol, " Holdup de liquido (-) ;", " Liquid holdup (-) ;"},
    {&detTRENDP::FVH, " Fracao Volumetrica Hidrato (-) ;", " Hydrate volumetric fraction (-) ;"},
    {&detTRENDP::bet, " Fracao vol. de liquido complementar (-) ;", " Complementary liquid vol. fraction (-) ;"},
    {&detTRENDP::ugs, " Velocidade superficial do gas (m/s) ;", " Gas superficial velocity (m/s) ;"},
    {&detTRENDP::uls, " Velocidade superficial do liquido (m/s) ;", " Liquid superficial velocity (m/s) ;"},
    {&detTRENDP::ug, " Velocidade do gas (m/s) ;", " Gas velocity (m/s) ;"},
    {&detTRENDP::ul, " Velocidade do liquido (m/s) ;", " Liquid velocity (m/s) ;"},
    {&detTRENDP::arra, " Indicador de arranjo de fases (-) ;", " Phase pattern indicator (-) ;"},
    {&detTRENDP::viscl, " Viscosidade do Liquido (cP) ;", " Liquid viscosity (cP) ;"},
    {&detTRENDP::viscg, " Viscosidade do Gas (cP) ;", " Gas viscosity (cP) ;"},
    {&detTRENDP::rhog, " Massa Especifica do Gas (kg/m3) ;", " Gas density (kg/m3) ;"},
    {&detTRENDP::rhol, " Massa Especifica do Liquido (kg/m3) ;", " Liquid density (kg/m3) ;"},
    {&detTRENDP::rhoMix, " Massa Especifica da Mistura (kg/m3) ;", " Mixture density (kg/m3) ;"},
    {&detTRENDP::masg, " Vazao Massica do Gas (kg/s) ;", " Gas mass flow rate (kg/s) ;"},
    {&detTRENDP::masl, " Vazao Massica do Liquido (kg/s) ;", " Liquid mass flow rate (kg/s) ;"},
    {&detTRENDP::c0, " Coeficiente de distribuição: C0 (-) ;", " Distribution coefficient: C0 (-) ;"},
    {&detTRENDP::ud, " Velocidade de escorregamento: Ud (m/s) ;", " Slip velocity: Ud (m/s) ;"},
    {&detTRENDP::RGO, " RGO (Sm3/Sm3) ;", " RGO (Sm3/Sm3) ;"},
    {&detTRENDP::deng, " Densidade do Gas (-) ;", " Gas specific gravity (-) ;"},
    {&detTRENDP::yco2, " Fracao Molar de CO2 (-) ;", " CO2 molar fraction (-) ;"},
    {&detTRENDP::calor, " Fluxo de calor entre escoamento e parede (W/m) ;", " Heat flow between flow and wall (W/m) ;"},
    {&detTRENDP::masstrans, " Transferencia de Massa entre Fases (kg / [s m]) ;", " Interphase mass transfer (kg / [s m]) ;"},
    {&detTRENDP::qlst, " Vazao volumetrica standard de oleo morto (Sm3/d) ;", " Standard dead oil volumetric flow rate (Sm3/d) ;"},
    {&detTRENDP::qlwst, " Vazao volumetrica standard de oleo morto + agua (Sm3/d) ;", " Standard dead oil + water volumetric flow rate (Sm3/d) ;"},
    {&detTRENDP::qlstTot, " Vazao volumetrica standard de oleo morto + agua + liquido complementar (Sm3/d) ;", " Standard dead oil + water + complementary liquid volumetric flow rate (Sm3/d) ;"},
    {&detTRENDP::qgst, " Vazao volumetrica standard de gas livre + dissolvido (Sm3/d) ;", " Standard free + dissolved gas volumetric flow rate (Sm3/d) ;"},
    {&detTRENDP::api, " Grau API (-) ;", " API gravity (-) ;"},
    {&detTRENDP::bsw, " BSW (-) ;", " BSW (-) ;"},
    {&detTRENDP::hidro, " Termo Hidrostatico (Pa/m) ;", " Hydrostatic term (Pa/m) ;"},
    {&detTRENDP::fric, " Termo Friccao (Pa/m) ;", " Friction term (Pa/m) ;"},
    {&detTRENDP::dengD, " Densidade Gas Dissolvido In Situ (-) ;", " In-situ dissolved gas specific gravity (-) ;"},
    {&detTRENDP::dengL, " Densidade Gas Livre In Situ (-) ;", " In-situ free gas specific gravity (-) ;"},
    {&detTRENDP::mlFonte, " Fonte massica - Liq. (Hidrocarb+Agua) (kg/s);", " Mass source - Liq. (Hydrocarbon+Water) (kg/s);"},
    {&detTRENDP::mgFonte, " Fonte massica - Gas (kg/s);", " Mass source - Gas (kg/s);"},
    {&detTRENDP::mcFonte, " Fonte massica - Liq. Complementar (kg/s);", " Mass source - Complementary liquid (kg/s);"},
    {&detTRENDP::dpB, " Incremento de pressao de Bombeio (kgf/cm2);", " Pump pressure increment (kgf/cm2);"},
    {&detTRENDP::potB, " Potencia de Bombeio (kW);", " Pump power (kW);"},
    {&detTRENDP::tempChokeJus, " Temperatura a Jusante do Choke de Superficie (C);", " Surface choke downstream temperature (C);"},
    {&detTRENDP::reyi, " Reynolds interno da mistura (-) ;", " Internal mixture Reynolds (-) ;"},
    {&detTRENDP::reye, " Reynolds externo (-) ;", " External Reynolds (-) ;"},
    {&detTRENDP::Fr, " Froud (-) ;", " Froude (-) ;"},
    {&detTRENDP::grashi, " Grashof interno da mistura (-) ;", " Internal mixture Grashof (-) ;"},
    {&detTRENDP::grashe, " Grashof externo (-) ;", " External Grashof (-) ;"},
    {&detTRENDP::nusi, " Nusselt interno da mistura (-) ;", " Internal mixture Nusselt (-) ;"},
    {&detTRENDP::nuse, " Nusselt externo (-) ;", " External Nusselt (-) ;"},
    {&detTRENDP::hi, " Coeficiente de pelicula interno da mistura (W/(m2.K)) ;", " Internal mixture film coefficient (W/(m2.K)) ;"},
    {&detTRENDP::he, " Coeficiente de pelicula externo (W/(m2.K)) ;", " External film coefficient (W/(m2.K)) ;"},
    {&detTRENDP::pri, " Prandtl interno da mistura (-) ;", " Internal mixture Prandtl (-) ;"},
    {&detTRENDP::pre, " Prandtl externo (-) ;", " External Prandtl (-) ;"},
    {&detTRENDP::Rs, " Razao de Solubilidade (-) ;", " Solubility ratio (-) ;"},
    {&detTRENDP::Bo, " Fator Volume de Formacao (-) ;", " Formation volume factor (-) ;"},
    {&detTRENDP::volMonM1PT, " Volume de liquido a montante da Master1, a PT, m3 ;", " Liquid volume upstream of Master1, at PT, m3 ;"},
    {&detTRENDP::volJusM1PT, " Volume de liquido a jusante da Master1, a PT, m3 ;", " Liquid volume downstream of Master1, at PT, m3 ;"},
    {&detTRENDP::volMonM1ST, " Volume de liquido a montante da Master1, standard, m3 ;", " Liquid volume upstream of Master1, standard, m3 ;"},
    {&detTRENDP::volJusM1ST, " Volume de liquido a jusante da Master1, standard, m3 ;", " Liquid volume downstream of Master1, standard, m3 ;"},
    {&detTRENDP::inventarioGas, " Inventario de Gas em toda a tubulação, standard, m3 ;", " Gas inventory in whole tubing, standard, m3 ;"},
    {&detTRENDP::inventarioLiq, " Inventario de Liquido em toda a tubulação, standard, m3 ;", " Liquid inventory in whole tubing, standard, m3 ;"},
    {&detTRENDP::diamInt, " Diametro Interno da tubulacao, m ;", " Tubing inner diameter, m ;"},
    {&detTRENDP::TempParede, " Temperatura Interna da Parede, C ;", " Internal wall temperature, C ;"},
    {&detTRENDP::subResfria, " Subresfriamento, C ;", " Subcooling, C ;"},
};

constexpr OptionalColumn<detTRENDG> kServiceColumns[] = {
    {&detTRENDG::pres, " Pressao (kgf/cm2) ;", " Pressure (kgf/cm2) ;"},
    {&detTRENDG::temp, " Temperatura (C) ;", " Temperature (C) ;"},
    {&detTRENDG::ugs, " Velocidade superficial do gas (m/s) ;", " Gas superficial velocity (m/s) ;"},
    {&detTRENDG::ug, " Velocidade do gas (m/s) ;", " Gas velocity (m/s) ;"},
    {&detTRENDG::tens, " Tensao Cisalhante (N/m2) ;", " Shear stress (N/m2) ;"},
    {&detTRENDG::viscg, " Viscosidade do Gas (cP) ;", " Gas viscosity (cP) ;"},
    {&detTRENDG::rhog, " Massa Especifica do Gas (kg/m3) ;", " Gas density (kg/m3) ;"},
    {&detTRENDG::masg, " Vazao Massica do Gas (kg/s) ;", " Gas mass flow rate (kg/s) ;"},
    {&detTRENDG::hidro, " Termo Hidrostatico (Pa/m) ;", " Hydrostatic term (Pa/m) ;"},
    {&detTRENDG::FVHG, " FVHG (-) ;", " FVHG (-) ;"},
    {&detTRENDG::fric, " Termo Friccao (Pa/m) ;", " Friction term (Pa/m) ;"},
    {&detTRENDG::calor, " Fluxo de calor entre escoamento e parede (W/m) ;", " Heat flow between flow and wall (W/m) ;"},
    {&detTRENDG::qgst, " Vazao volumetrica standard de Gas (Sm3/d) ;", " Standard gas volumetric flow rate (Sm3/d) ;"},
    {&detTRENDG::pEstagVGL, " Pressao de Estagnacao VGL (kgf/cm²) ;", " VGL stagnation pressure (kgf/cm2) ;"},
    {&detTRENDG::tEstagVGL, " Temperatura de Estagnacao VGL (C) ;", " VGL stagnation temperature (C) ;"},
    {&detTRENDG::pGargVGL, " Pressao na Garganta VGL (kgf/cm²) ;", " VGL throat pressure (kgf/cm2) ;"},
    {&detTRENDG::tGargVGL, " Temperatura na Garganta VGL (C) ;", " VGL throat temperature (C) ;"},
    {&detTRENDG::velgarg, " Velocidade na VGL (m/s) ;", " Velocity in VGL (m/s) ;"},
    {&detTRENDG::qVGL, " Vazao volumetrica na VGL (mÂ³/s) ;", " Volumetric flow rate in VGL (m3/s) ;"},
    {&detTRENDG::reyi, " Reynolds interno (-) ;", " Internal Reynolds (-) ;"},
    {&detTRENDG::reye, " Reynolds externo (-) ;", " External Reynolds (-) ;"},
    {&detTRENDG::grashi, " Grashof interno (-) ;", " Internal Grashof (-) ;"},
    {&detTRENDG::grashe, " Grashof externo (-) ;", " External Grashof (-) ;"},
    {&detTRENDG::nusi, " Nusselt interno (-) ;", " Internal Nusselt (-) ;"},
    {&detTRENDG::nuse, " Nusselt externo (-) ;", " External Nusselt (-) ;"},
    {&detTRENDG::hi, " Coeficiente de pelicula interno (W/[m2 K]) ;", " Internal film coefficient (W/[m2 K]) ;"},
    {&detTRENDG::he, " Coeficiente de pelicula externo (W/[m2 K]) ;", " External film coefficient (W/[m2 K]) ;"},
    {&detTRENDG::pri, " Prandtl interno (-) ;", " Internal Prandtl (-) ;"},
    {&detTRENDG::pre, " Prandtl externo (-) ;", " External Prandtl (-) ;"},
    {&detTRENDG::diamInt, " Diametro Interno da tubulacao, m ;", " Tubing inner diameter, m ;"},
    {&detTRENDG::TempParede, " Temperatura Interna da Parede, C ;", " Internal wall temperature, C ;"},
    {&detTRENDG::subResfria, " Subresfriamento, C ;", " Subcooling, C ;"},
};

/// Written as one group when dadosParafina is set.
constexpr FixedColumn kParaffinColumns[] = {
    {" TIAC (C) C;", " TIAC (C) C;"},
    {" Cp Parafina (J/[kg C]) C;", " Paraffin Cp (J/[kg C]) C;"},
    {" Condutividade Termica Parafina (W / [m K]) C;", " Paraffin thermal conductivity (W / [m K]) C;"},
    {" Massa Especifica Parafina (kg/m3) C;", " Paraffin density (kg/m3) C;"},
    {" Massa molar do Liquido Parafina (kg/mol) C;", " Paraffin liquid molar mass (kg/mol) C;"},
    {" Difusividade Massica Parafina (m2/s) C;", " Paraffin mass diffusivity (m2/s) C;"},
    {" Fluxo Massico de Parafina Total (kg/(m2-s)) C;", " Total paraffin mass flux (kg/(m2-s)) C;"},
    {" Fluxo Massico de Parafina por Difusao (kg/(m2-s)) C;", " Paraffin diffusive mass flux (kg/(m2-s)) C;"},
    {" Vazao Massica de Parafina por Difusao (kg/(s)) C;", " Paraffin mass rate (kg/(s)) C;"},
    {" Gradiente de concentracao de parafina (1/m) C;", " Paraffin concentration gradient (1/m) C;"},
    {" Condutividade do deposito (W/(m-K)) C;", " Deposit conductivity (W/(m-K)) C;"},
    {" Temperatura da Interface do deposito (C) C;", " Deposit interface temperature (C) C;"},
};

template <typename Trend, std::size_t Count>
void writeOptionalColumns(ostream &trendFile, const CaptionTranslator &translate,
                          const Trend &trend,
                          const OptionalColumn<Trend> (&columns)[Count]) {
    for (const OptionalColumn<Trend> &column : columns)
        if (trend.*column.enabled == 1)
            trendFile << translate(column.ptBrText, column.enText);
}

template <std::size_t Count>
void writeFixedColumns(ostream &trendFile, const CaptionTranslator &translate,
                       const FixedColumn (&columns)[Count]) {
    for (const FixedColumn &column : columns)
        trendFile << translate(column.ptBrText, column.enText);
}

/// The three wave families the characteristic analysis reports.
constexpr int kWaveFamilyCount = 3;

void writeWaveCelerityColumns(ostream &trendFile, const CaptionTranslator &translate) {
    for (int waveFamily = 0; waveFamily < kWaveFamilyCount; waveFamily++) {
        trendFile << translate(" Celeridade, familia de onda ", " Celerity, wave family ") << waveFamily << translate(" m/s ;", " m/s ;");
    }
}

void writeEigenvectorColumns(ostream &trendFile, const CaptionTranslator &translate) {
    for (int waveFamily = 0; waveFamily < kWaveFamilyCount; waveFamily++) {
        for (int eigenvectorTerm = 0; eigenvectorTerm < kWaveFamilyCount; eigenvectorTerm++)
            trendFile << translate(" Componente do autovetor, condicao adiabatica, familia de onda = ", " Eigenvector component, adiabatic condition, wave family = ") << waveFamily << translate("termo = ", " term = ") << eigenvectorTerm << " ;";
    }
}

void writeFluctuationColumns(ostream &trendFile, const CaptionTranslator &translate) {
    for (int waveFamily = 0; waveFamily < kWaveFamilyCount; waveFamily++) {
        trendFile << translate(" Componente de flutuacao da familia de onda ", " Fluctuation component of wave family ") << waveFamily << " ;";
    }
}

/// Distance from the line origin to the trend cell, accumulated cell by cell.
///
/// Templated on the cell type because the production and service meshes are
/// different structs that happen to share a dx. The accumulation order is the
/// baseline's and must stay that way: reordering a floating-point sum changes
/// its result.
template <typename Cell>
double lengthFromOrigin(const Cell *cells, int lastCellIndex) {
    double length = 0;
    for (int cell = 0; cell <= lastCellIndex; cell++)
        length += cells[cell].dx;
    return length;
}

/// The metadata block both line-trend headers open with.
///
/// The two differ only in which mesh measures the distance and what that
/// distance is measured from -- the wellbore for production, the platform for
/// service.
template <typename Trend, typename Cell>
void writeLineTrendMetadata(ostream &trendFile, const CaptionTranslator &translate,
                            const TrendState &state, const Trend &trend,
                            const Cell *cells, const char *originPtBrText,
                            const char *originEnText) {
    const double length = lengthFromOrigin(cells, trend.posic);
    trendFile << translate(originPtBrText, originEnText) << length << endl;
    trendFile << translate("# Rotulo = ", "# Label = ") << trend.rotulo << endl;
    trendFile << translate("# Indice da Celula = ", "# Cell index = ") << trend.posic << endl;
    if (state.branchIndex < 0 && state.inputData.AP == 1)
        trendFile << translate(" Sequencia AP ;", " SA sequence ;");
    trendFile << translate(" Tempo (s) ;", " Time (s) ;");
}

/// Every cross-section trend carries the same two columns, whichever line it
/// belongs to.
void writeCrossSectionCaptions(ostream &trendFile, const CaptionTranslator &translate,
                               const detTRENDTrans &trend) {
    trendFile << translate("# Rotulo = ", "# Label = ") << trend.rotulo << endl;
    trendFile << translate(" Tempo (s) ; ", " Time (s) ; ");
    trendFile << translate(" Temperatura (C) ;", " Temperature (C) ;");
}

// ------------------------------------------------------------- file names --

/// Builds the name of a production- or service-line trend file.
///
/// One builder for all four line writers. It could only become one once the
/// captions file and the data file agreed on how to spell the position: the
/// header used to round it, and round() on an int yields a double, which
/// ostream prints in exponent form at or above 1e6 -- sending the two halves of
/// the same trend to two different files.
string lineTrendFileName(const TrendState &state, const char *prefix,
                         int position, int networkIndex) {
    ostringstream fileNameStream;
    fileNameStream << pathPrefixoArqSaida;
    if (state.branchIndex < 0)
        fileNameStream << prefix << (state.inputData.AP == 1 ? "-AP-" : "-") << position;
    else
        fileNameStream << "Tramo" << state.branchIndex << "-R-" << networkIndex
                       << "-" << prefix << "-" << position;
    fileNameStream << ".dat";
    return fileNameStream.str();
}

/// Where a cross-section trend is sampled: the cell, the wall layer, and the
/// discretisation point inside that layer.
struct CrossSectionPosition {
    int cell;
    int layer;
    int discretization;
};

CrossSectionPosition crossSectionPositionOf(const detTRENDTrans &trend) {
    return {trend.comp, trend.camada - 1, trend.discre - 1};
}

/// Builds the name of a cross-section temperature trend file.
///
/// Also one builder for both writers, and also only after a fix: the service
/// writer spelled its branch-path file TENDTRANSP, so inside a network it
/// overwrote the production trend of the same branch.
string crossSectionTrendFileName(const TrendState &state, const char *prefix,
                                 const CrossSectionPosition &position) {
    ostringstream fileNameStream;
    fileNameStream << pathPrefixoArqSaida;
    if (state.branchIndex >= 0)
        fileNameStream << "Tramo" << state.branchIndex << "-";
    fileNameStream << prefix << "-" << position.cell << "-" << position.layer
                   << "-" << position.discretization << ".dat";
    return fileNameStream.str();
}

// --------------------------------------------------------------- reporting --

/// Records a produced file in the profile report.
void reportProducedFile(const TrendState &state, const string &fileName) {
    // caso nao seja simulacao POCO_INJETOR
    if (state.inputData.tipoSimulacao != tipoSimulacao_t::poco_injetor) {
        arqRelatorioPerfis << fileName.c_str() << endl;
        arqRelatorioPerfis.flush();
    }
}

// --------------------------------------------------------------- skeletons --

/// Opens a trend file, lets the caller fill it, then closes and reports it.
///
/// Every trend file this module produces goes through here, which is the point:
/// reporting the file to the profile index is not something a writer can forget
/// to do, because it is not a writer's job any more.
template <typename WriteContent>
void produceTrendFile(const TrendState &state, const string &fileName,
                      ios_base::openmode mode, bool blankLineBeforeClose,
                      WriteContent writeContent) {
    ofstream trendFile(fileName.c_str(), mode);
    writeContent(trendFile);
    if (blankLineBeforeClose)
        trendFile << endl;
    trendFile.close();
    reportProducedFile(state, fileName);
}

/// The block of samples one row writer emits.
///
/// Rows are read straight out of the trend buffer. An earlier form copied the
/// window into a FullMtx first, allocating one vector per row on every call --
/// some four thousand calls per run -- to hand the same values to the same loop
/// in the same order.
struct SampleWindow {
    double **samples;
    int firstRow;
    int rowCount;
    int columnCount;
};

/// Rows recorded since the last write, for one series of a trend group.
int pendingRowCount(const TrendSeries &series, int trendIndex) {
    return series.count[trendIndex] - series.countBase[trendIndex] + 1;
}

SampleWindow windowOf(const TrendSeries &series, int trendIndex, int columnCount) {
    return {series.samples[trendIndex], series.countBase[trendIndex],
            pendingRowCount(series, trendIndex), columnCount};
}

/// Rows at or below this value in column zero end the file.
constexpr double kEndOfSamplesMarker = -9999;

/// Field width and precision of every value written to a trend file.
constexpr int kValueWidth = 20;
constexpr int kValuePrecision = 19;

/// The skeleton the four header writers share.
///
/// writeCaptions is the hook, and it is a template parameter rather than a
/// function pointer or a virtual, so the compiler inlines it and no indirect
/// call appears in the generated code (FR-022).
///
/// blankLineBeforeClose is a value, not a mode flag. Only the service header
/// asks for it, and it lands OUTSIDE the print-pass guard, so a service trend
/// file starts with a blank line even on the passes that write no captions.
/// That asymmetry comes from the baseline and is preserved.
template <typename WriteCaptions>
void writeTrendHeaderFile(const TrendState &state, const string &fileName,
                          WriteCaptions writeCaptions, bool blankLineBeforeClose) {
    produceTrendFile(state, fileName, ios_base::out, blankLineBeforeClose,
                     [&](ofstream &trendFile) {
                         if (state.printPassCount == 1) {
                             writeCaptions(trendFile);
                             trendFile << endl;
                         }
                     });
}

/// The skeleton the four row writers share.
///
/// The AP sequence column exists only for the line trends. It is emitted once
/// ahead of the data columns rather than tested inside the column loop, which
/// is one conditional fewer per value and a truer description: it is a prefix,
/// not a special case of a data column.
void writeTrendRowsFile(const TrendState &state, const string &fileName,
                        const SampleWindow &window, bool supportsApSequence,
                        bool blankLineBeforeClose) {
    produceTrendFile(
        state, fileName, ios_base::app, blankLineBeforeClose,
        [&](ofstream &trendFile) {
            const bool apSequenceColumn =
                supportsApSequence && state.branchIndex < 0 && state.inputData.AP == 1;

            trendFile.precision(kValuePrecision);
            for (int rowIndex = 0; rowIndex < window.rowCount; rowIndex++) {
                const double *row = window.samples[window.firstRow + rowIndex];
                if (row[0] <= kEndOfSamplesMarker)
                    break;
                if (apSequenceColumn) {
                    trendFile.width(kValueWidth);
                    trendFile << (*state.globals).sequenciaAP << " ; ";
                }
                for (int columnIndex = 0; columnIndex < window.columnCount; columnIndex++) {
                    trendFile.width(kValueWidth);
                    trendFile << row[columnIndex] << " ; ";
                }
                trendFile << endl;
            }
        });
}

/// The two columns every cross-section trend carries: time and temperature.
constexpr int kCrossSectionColumnCount = 2;

} // namespace

// ------------------------------------------------------------ entry points --

void writeProductionTrendHeader(const TrendState &state, int trendIndex, int networkIndex) {
    if (state.inputData.ntendp > 0) {
        const CaptionTranslator translate{state.inputData.idiomaSaida};
        const detTRENDP &trend = state.inputData.trendp[trendIndex];
        writeTrendHeaderFile(
            state, lineTrendFileName(state, "TENDP", trend.comp, networkIndex),
            [&](ofstream &trendFile) {
                writeLineTrendMetadata(trendFile, translate, state, trend,
                                       state.inputData.celp,
                                       "# Comprimento a partir do Fundo de Poco (m) = ",
                                       "# Length from Bottomhole (m) = ");
                writeOptionalColumns(trendFile, translate, trend, kProductionColumns);
                if (trend.dadosParafina == 1)
                    writeFixedColumns(trendFile, translate, kParaffinColumns);
                if (trend.autoVal == 1)
                    writeWaveCelerityColumns(trendFile, translate);
                if (trend.autoVel == 1)
                    writeEigenvectorColumns(trendFile, translate);
                if (trend.flutuacao == 1)
                    writeFluctuationColumns(trendFile, translate);
            },
            /*blankLineBeforeClose=*/false);
    }
}

void writeProductionTrendRows(const TrendState &state, int trendIndex, int networkIndex) {
    if (state.inputData.ntendp > 0) {
        const detTRENDP &trend = state.inputData.trendp[trendIndex];
        writeTrendRowsFile(
            state, lineTrendFileName(state, "TENDP", trend.comp, networkIndex),
            windowOf(state.production, trendIndex, state.inputData.nvartrendp[trendIndex] + 1),
            /*supportsApSequence=*/true, /*blankLineBeforeClose=*/false);
    }
}

void writeServiceTrendHeader(const TrendState &state, int trendIndex, int networkIndex) {
    if (state.inputData.ntendg > 0 && state.inputData.lingas > 0) {
        const CaptionTranslator translate{state.inputData.idiomaSaida};
        const detTRENDG &trend = state.inputData.trendg[trendIndex];
        writeTrendHeaderFile(
            state, lineTrendFileName(state, "TENDG", trend.comp, networkIndex),
            [&](ofstream &trendFile) {
                writeLineTrendMetadata(trendFile, translate, state, trend,
                                       state.inputData.celg,
                                       "# Comprimento a partir da Plataforma (m) = ",
                                       "# Length from Platform (m) = ");
                writeOptionalColumns(trendFile, translate, trend, kServiceColumns);
            },
            /*blankLineBeforeClose=*/true);
    }
}

void writeServiceTrendRows(const TrendState &state, int trendIndex, int networkIndex) {
    if (state.inputData.ntendg > 0 && state.inputData.lingas > 0) {
        const detTRENDG &trend = state.inputData.trendg[trendIndex];
        writeTrendRowsFile(
            state, lineTrendFileName(state, "TENDG", trend.comp, networkIndex),
            windowOf(state.service, trendIndex, state.inputData.nvartrendg[trendIndex] + 1),
            /*supportsApSequence=*/true, /*blankLineBeforeClose=*/false);
    }
}

void writeProductionCrossSectionTrendHeader(const TrendState &state, int trendIndex) {
    if (state.inputData.ntendtransp > 0) {
        const CaptionTranslator translate{state.inputData.idiomaSaida};
        const detTRENDTrans &trend = state.inputData.trendtransp[trendIndex];
        writeTrendHeaderFile(
            state,
            crossSectionTrendFileName(state, "TENDTRANSP", crossSectionPositionOf(trend)),
            [&](ofstream &trendFile) { writeCrossSectionCaptions(trendFile, translate, trend); },
            /*blankLineBeforeClose=*/false);
    }
}

void writeProductionCrossSectionTrendRows(const TrendState &state, int trendIndex) {
    if (state.inputData.ntendtransp > 0) {
        const detTRENDTrans &trend = state.inputData.trendtransp[trendIndex];
        writeTrendRowsFile(
            state,
            crossSectionTrendFileName(state, "TENDTRANSP", crossSectionPositionOf(trend)),
            windowOf(state.productionCrossSection, trendIndex, kCrossSectionColumnCount),
            /*supportsApSequence=*/false, /*blankLineBeforeClose=*/true);
    }
}

void writeServiceCrossSectionTrendHeader(const TrendState &state, int trendIndex) {
    if (state.inputData.ntendtransg > 0 && state.inputData.lingas > 0) {
        const CaptionTranslator translate{state.inputData.idiomaSaida};
        const detTRENDTrans &trend = state.inputData.trendtransg[trendIndex];
        writeTrendHeaderFile(
            state,
            crossSectionTrendFileName(state, "TENDTRANSG", crossSectionPositionOf(trend)),
            [&](ofstream &trendFile) { writeCrossSectionCaptions(trendFile, translate, trend); },
            /*blankLineBeforeClose=*/false);
    }
}

void writeServiceCrossSectionTrendRows(const TrendState &state, int trendIndex) {
    if (state.inputData.ntendtransg > 0 && state.inputData.lingas > 0) {
        const detTRENDTrans &trend = state.inputData.trendtransg[trendIndex];
        writeTrendRowsFile(
            state,
            crossSectionTrendFileName(state, "TENDTRANSG", crossSectionPositionOf(trend)),
            windowOf(state.serviceCrossSection, trendIndex, kCrossSectionColumnCount),
            /*supportsApSequence=*/false, /*blankLineBeforeClose=*/true);
    }
}

} // namespace trendoutput
