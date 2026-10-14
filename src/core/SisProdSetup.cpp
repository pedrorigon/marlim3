/// How a production system is assembled from its input: the cells and their
/// sources, the gas-lift line, the tables and the output bookkeeping, in the
/// order montasistema runs them. These are SProd's own methods; SisProd.cpp
/// holds its lifecycle.
#include "SisProd.h"
#include "OutputI18n.h"
#include "SisProdConstants.h"
#include "UnitConversions.h"
#include <math.h>
#include <numbers>

using enum sisprod::AccessoryKind;

namespace {

/// Calls apply on each fluid the cell's source carries: that of a liquid,
/// inflow-performance or multiple source, or, for a porous reservoir, its own
/// and those of its cells. A gas injection's fluid is not visited. An apply that
/// also takes a pressure and a temperature gets the cell's for the source's own
/// fluid and each reservoir cell's for that cell's.
template <typename Apply>
void forEachSourceFluid(Cel &cell, Apply &&apply) {
    auto visit = [&](ProFlu &fluid, double pressure, double temperature) {
        if constexpr (std::invocable<Apply &, ProFlu &>)
            apply(fluid);
        else
            apply(fluid, pressure, temperature);
    };
    acessorio &source = cell.acsr;
    if (source.tipo == kAccessoryLiquidInjection) {
        visit(source.injl.FluidoPro, cell.pres, cell.temp);
    } else if (source.tipo == kAccessoryInflowPerformance) {
        visit(source.ipr.FluidoPro, cell.pres, cell.temp);
    } else if (source.tipo == kAccessoryMultipleSource) {
        visit(source.injm.FluidoPro, cell.pres, cell.temp);
    } else if (source.tipo == kAccessoryRadialPorous) {
        visit(source.radialPoro.flup, cell.pres, cell.temp);
        for (int k = 0; k < source.radialPoro.ncel; k++)
            visit(source.radialPoro.celula[k].flup, source.radialPoro.celula[k].Pcamada, source.radialPoro.tRes);
    } else if (source.tipo == kAccessoryPorous2D) {
        visit(source.poroso2D.dados.flup, cell.pres, cell.temp);
        for (int k = 0; k < source.poroso2D.dados.transfer.ncel; k++)
            visit(source.poroso2D.dados.transfer.celula[k].flup, source.poroso2D.dados.transfer.celula[k].Pcamada,
                  source.poroso2D.dados.transfer.tRes);
        for (int k = 0; k < source.poroso2D.malha.nele; k++)
            visit(source.poroso2D.malha.mlh2d[k].flup, source.poroso2D.malha.mlh2d[k].cel2D.presC,
                  source.poroso2D.malha.mlh2d[k].tRes);
    }
}

}  // namespace

/// Points every cell fluid, and every source fluid it carries, at the bubble-point
/// tables read from the PVTSim file, with their size, and switches them to saturation
/// model 4.
void SProd::assignPvtSimBubbleTablesToCells() {
    auto assignTables = [&](ProFlu &fluid) {
        fluid.PBPVTSim = tables.pvtSimBubblePressures.data();
        fluid.TBPVTSim = tables.pvtSimBubbleTemperatures.data();
        fluid.npontosB = static_cast<int>(tables.pvtSimBubbleTemperatures.size());
        fluid.corrSat = 4;
        fluid.awaitingSaturationTable = false;
    };
    for (int i = 0; i <= ncel; i++) {
        assignTables(celula[i].flui);
        forEachSourceFluid(celula[i], assignTables);
    }
}

namespace {

/// Reads one row of a PVTSim table file: skips words up to `name`, then converts every number on the rest of that
/// line. A missing row, or one with fewer than two numbers, ends the run naming the file.
template <typename Convert>
vector<double> readPvtSimRow(ifstream &file, string &word, const char *name, const string &fileName, Convert convert) {
    while (word != name && file >> word) {
    }
    if (word != name)
        logger.log_write_logs_and_exit(LOGGER_FALHA, LOG_ERR_PARSE_BUSINESS_RULE_VALIDATION,
                                       string("PVTSim file without the row ") + name, "pvtsimArq", fileName);
    string line;
    getline(file, line);
    vector<double> values;
    for (char *token = strtok(line.data(), " ,()=\\"); token != nullptr; token = strtok(nullptr, " ,()=\\")) {
        char *end = nullptr;
        double value = strtod(token, &end);
        if (end != token)
            values.push_back(convert(value));
    }
    if (values.size() < 2)
        logger.log_write_logs_and_exit(LOGGER_FALHA, LOG_ERR_PARSE_BUSINESS_RULE_VALIDATION,
                                       string("PVTSim row with fewer than two values: ") + name, "pvtsimArq", fileName);
    return values;
}

/// True when every value is above the one before it.
bool strictlyIncreasing(const vector<double> &axis) {
    for (size_t k = 1; k < axis.size(); k++)
        if (!(axis[k] > axis[k - 1]))
            return false;
    return true;
}

/// Writes a table to the file at `path`, replacing it, after a `# header` line naming the quantity when one is given.
void writeTable(const string &path, const FullMtx<double> &table, const string &header = "") {
    ofstream file(path.c_str(), ios_base::out);
    if (!header.empty())
        file << "# " << header << "\n";
    file << table;
}

}  // namespace

/// Reads the bubble-point curve from the PVTSim file, points the cell fluids at it
/// and writes perfilBolha; with tabRSPB on, also reads the solution gas-oil ratio
/// table and writes perfilRSLivia.
void SProd::loadPvtSimSaturationTables() {
    string impfile;
    impfile = arq.pvtsimarq;
    string dadosMR = impfile;
    ifstream lendoPVTSim(dadosMR.c_str(), ios_base::in);
    string chave;
    char line[4000];
    lendoPVTSim.get(line, 4000);
    lendoPVTSim >> chave;
    // The file's own grid: its rows set the table sizes, apart from the input table's resolution.
    vector<double> PresPVTSim = readPvtSimRow(lendoPVTSim, chave, "PRESSURE", dadosMR,
                                              [](double value) { return value / units::kPascalPerKgfPerCm2; });
    tables.pvtSimBubblePressures = readPvtSimRow(lendoPVTSim, chave, "BUBBLEPRESSURES", dadosMR,
                                    [](double value) { return value / units::kPascalPerPsi; });
    tables.pvtSimBubbleTemperatures =
        readPvtSimRow(lendoPVTSim, chave, "BUBBLETEMPERATURES", dadosMR, [](double value) { return value; });
    if (tables.pvtSimBubblePressures.size() != tables.pvtSimBubbleTemperatures.size())
        logger.log_write_logs_and_exit(LOGGER_FALHA, LOG_ERR_PARSE_BUSINESS_RULE_VALIDATION,
                                       "PVTSim bubble-point curve with different numbers of pressures and temperatures",
                                       "pvtsimArq", dadosMR);
    if (!strictlyIncreasing(PresPVTSim) || !strictlyIncreasing(tables.pvtSimBubbleTemperatures))
        logger.log_write_logs_and_exit(LOGGER_FALHA, LOG_ERR_PARSE_BUSINESS_RULE_VALIDATION,
                                       "PVTSim pressure or bubble-temperature axis not strictly increasing", "pvtsimArq",
                                       dadosMR);
    const int pressureCount = static_cast<int>(PresPVTSim.size());
    const int bubbleCount = static_cast<int>(tables.pvtSimBubbleTemperatures.size());
    assignPvtSimBubbleTablesToCells();

    FullMtx<double> BolhaTemp(bubbleCount + 1, 2);
    for (int i = 0; i < bubbleCount; i++) {
        BolhaTemp[i][0] = tables.pvtSimBubbleTemperatures[i];
        BolhaTemp[i][1] = tables.pvtSimBubblePressures[i];
    }
    writeTable(pathPrefixoArqSaida + "perfilBolha", BolhaTemp);

    if (arq.tabRSPB == 1) {
        tables.liviaSolutionGasRatioTable.allocate(pressureCount + 1, bubbleCount + 1);
        for (int j = 1; j <= bubbleCount; j++)
            tables.liviaSolutionGasRatioTable[0][j] = tables.pvtSimBubbleTemperatures[j - 1];
        for (int i = 1; i <= pressureCount; i++) {
            tables.liviaSolutionGasRatioTable[i][0] = PresPVTSim[i - 1];
            for (int j = 1; j <= bubbleCount; j++) {
                if (tables.pvtSimBubbleTemperatures[j - 1] > -10) {
                    double TFa = 1.8 * tables.pvtSimBubbleTemperatures[j - 1] + 32;

                    constexpr double A0 = 6542.69213 * 1e-11;
                    constexpr double A1 = 2.60464618;
                    constexpr double A2 = 1.21334544;
                    constexpr double A3 = 0.16464125;
                    constexpr double A4 = 1.19382967;
                    constexpr double B0 = 87 * 1e-8;
                    constexpr double B1 = 4.77390016;
                    constexpr double B2 = 1.77267703;
                    constexpr double B3 = -0.73072038;
                    constexpr double B4 = 1.58093857;
                    constexpr double C0 = 3226.09 * 1e-6;
                    constexpr double C1 = 0.09281075;
                    constexpr double C2 = 0.13633665;
                    constexpr double C3 = 0.09634381;
                    constexpr double C4 = 0.53238728;
                    constexpr double D0 = 1.00544053;
                    constexpr double D1 = -0.00134177;
                    constexpr double D2 = 0.51839397;

                    double yco2 = celula[0].flui.yco2;
                    double Deng = celula[0].flui.Deng;
                    double API = celula[0].flui.API;
                    double multCor = (D0 + D1 * yco2 * pow(TFa, D2));
                    double pbtemp = tables.pvtSimBubblePressures[j - 1];

                    double pcor = tables.liviaSolutionGasRatioTable[i][0] * units::kPsiPerKgfPerCm2 * multCor;
                    double pr = pcor / pbtemp;

                    double a1 = A0 * pow(Deng, A1) * pow(API, A2) * pow(TFa, A3) * pow(pbtemp, A4);
                    double a2 = B0 * pow(Deng, B1) * pow(API, B2) * pow(TFa, B3) * pow(pbtemp, B4);
                    double a3 = C0 * pow(Deng, C1) * pow(API, C2) * pow(TFa, C3) * pow(pbtemp, C4);

                    double Rsr = a1 * pow(pr, a2) + (1. - a1) * pow(pr, a3);
                    double rstemp = Rsr;
                    if (rstemp < 0.)
                        rstemp = 0.;
                    else if (tables.liviaSolutionGasRatioTable[i][0] * units::kPsiPerKgfPerCm2 > pbtemp)
                        rstemp = 1.;
                    else if (rstemp > 1.)
                        rstemp = 1.;
                    else
                        rstemp = Rsr;
                    tables.liviaSolutionGasRatioTable[i][j] = rstemp;
                } else
                    tables.liviaSolutionGasRatioTable[i][j] = 0;
            }
        }
        FullMtx<double> RSTemp(pressureCount + 1, bubbleCount + 1);
        for (int j = 1; j <= bubbleCount; j++)
            RSTemp[0][j] = tables.liviaSolutionGasRatioTable[0][j];
        for (int i = 1; i <= pressureCount; i++) {
            RSTemp[i][0] = tables.liviaSolutionGasRatioTable[i][0];
            for (int j = 1; j <= bubbleCount; j++)
                RSTemp[i][j] = tables.liviaSolutionGasRatioTable[i][j];
        }
        // The table holds the reduced ratio Rs/Rsb, which has no unit to convert.
        writeTable(pathPrefixoArqSaida + "perfilRSLivia", RSTemp,
                   output_i18n::tr(arq.idiomaSaida,
                                   "Razao de solubilidade reduzida Rs/Rsb (-); primeira coluna: pressao (kgf/cm2); primeira linha: temperatura (C)",
                                   "Reduced solution gas-oil ratio Rs/Rsb (-); first column: pressure (kgf/cm2); first row: temperature (C)"));

        auto pointAtRatioTable = [&](ProFlu &fluid) {
            fluid.TabRSLivia = tables.liviaSolutionGasRatioTable.data();
            fluid.rsTablePressureCount = pressureCount;
            fluid.rsTableTemperatureCount = bubbleCount;
            fluid.tabRSPB = 1;
        };
        for (int i = 0; i <= ncel; i++) {
            pointAtRatioTable(celula[i].flui);
            forEachSourceFluid(celula[i], pointAtRatioTable);
        }
    }
}

/// Builds the bubble-point curve and the solution gas-oil ratio table from the
/// fluid correlations over the input table's pressure-temperature grid, writes
/// perfilBolha and perfilRSLivia, and points the cell fluids at both, with the
/// bubble curve's size.
void SProd::generateSaturationTablesFromCorrelations() {
    int ndiv = arq.tabent.npont - 1;
    double pteste = arq.tabent.pmin;
    double tteste = arq.tabent.tmin;
    double dpteste = (arq.tabent.pmax - pteste) / ndiv;
    double dtteste = (arq.tabent.tmax - tteste) / ndiv;
    FullMtx<double> RSTemp(ndiv + 2, ndiv + 2);
    FullMtx<double> PBTemp(ndiv + 2, 2);
    tables.pvtSimBubbleTemperatures = vector<double>(ndiv + 1);
    tables.pvtSimBubblePressures = vector<double>(ndiv + 1);
    tables.liviaSolutionGasRatioTable.allocate(ndiv + 2, ndiv + 2);
    double ttestepb = tteste;
    for (int i = 0; i <= ndiv; i++) {
        tables.pvtSimBubbleTemperatures[i] = ttestepb;
        tables.pvtSimBubblePressures[i] = arq.flup[0].PB(pteste, ttestepb);
        PBTemp[i][0] = tables.pvtSimBubbleTemperatures[i];
        PBTemp[i][1] = tables.pvtSimBubblePressures[i];
        ttestepb += dtteste;
    }

    for (int i = 1; i <= ndiv + 1; i++) {
        tables.liviaSolutionGasRatioTable[i][0] = pteste;
        RSTemp[i][0] = tables.liviaSolutionGasRatioTable[i][0];
        ttestepb = tteste;
        for (int j = 1; j <= ndiv + 1; j++) {
            tables.liviaSolutionGasRatioTable[0][j] = ttestepb;
            tables.liviaSolutionGasRatioTable[i][j] = arq.flup[0].RS(pteste, ttestepb);
            RSTemp[0][j] = tables.liviaSolutionGasRatioTable[0][j];
            RSTemp[i][j] =
                tables.liviaSolutionGasRatioTable[i][j] * units::kBarrelPerCubicMetre / units::kCubicFootPerCubicMetre;
            ttestepb += dtteste;
        }
        pteste += dpteste;
    }

    writeTable(pathPrefixoArqSaida + "perfilBolha", PBTemp);
    writeTable(pathPrefixoArqSaida + "perfilRSLivia", RSTemp,
               output_i18n::tr(arq.idiomaSaida,
                               "Razao de solubilidade Rs (m3/m3); primeira coluna: pressao (kgf/cm2); primeira linha: temperatura (C)",
                               "Solution gas-oil ratio Rs (m3/m3); first column: pressure (kgf/cm2); first row: temperature (C)"));

    auto pointAtTables = [&](ProFlu &fluid) {
        fluid.PBPVTSim = tables.pvtSimBubblePressures.data();
        fluid.TBPVTSim = tables.pvtSimBubbleTemperatures.data();
        fluid.npontosB = static_cast<int>(tables.pvtSimBubbleTemperatures.size());
        fluid.TabRSLivia = tables.liviaSolutionGasRatioTable.data();
        fluid.rsTablePressureCount = ndiv + 1;
        fluid.rsTableTemperatureCount = ndiv + 1;
        fluid.tabRSPB = 1;
    };
    for (int i = 0; i <= ncel; i++) {
        pointAtTables(celula[i].flui);
        forEachSourceFluid(celula[i], pointAtTables);
    }

}

/// Copies the run configuration into the members, checks that an injection well
/// has the IPR its boundary condition needs, builds the Cp and JTL tables, hands
/// the fluid constants to every fluid, generates the pipe and the production
/// cells, and places every source -- including the extra gas sources at
/// extraSourceLengths, whose cell indices it records in extraSourceCells.
void SProd::buildProductionCells(double *extraSourceLengths, int *extraSourceCells, int extraSourceCount) {
    indTramo = -1;
    transient.printCounter = 0;
    transient.activeChokeDelay = 0.;
    transient.previousActiveChokeDelay = 0.;
    modoPerm = 0;
    steadyIteration = 0;
    masChkSup = 0;
    dt = arq.dtmax;
    gasLift.interfaceTimeStep = arq.dtmax;
    transient.finalTime = arq.tfinal;
    transient.massTransferModel = 0;
    transient.movingPigCount = 0;
    transient.previousMovingPigCount = transient.movingPigCount;
    transient.rollbackRequested = 0;
    networkCoupling.previousSurfaceChokeActive = 0;
    networkCoupling.reverseLiquidFraction = 0.;
    transient.trackGasGravity = arq.trackDeng;
    steadySearch.productionFluidCount = arq.nfluP;
    CalcLat = arq.latente;
    if (arq.flashCompleto == 1)
        arq.trackRGO = 1;
    transient.trackGasOilRatio = arq.trackRGO;
    const int injPoc = arq.pocinjec + arq.condpocinj.tipoFlui;
    arq.fluc.injPoc = injPoc;
    if (arq.flashCompleto == 1)
        arq.fluc.injPoc = 0;
    (*vg1dSP).localtiny = arq.mono;
    (*vg1dSP).CritCond = arq.critcond;
    networkCoupling.reverseQuality = -1;
    if (injPoc >= 1 && arq.nipr == 0 && arq.condpocinj.CC != 3 && arq.condpocinj.CC != 4 && arq.condpocinj.CC != 5)
        NumError(
            "O simulador esta no modo Injecao de Agua em uma condicao de contorno que pede uma IPR e nao foi incluido nenhuma IPR no sistema");
    else if (injPoc >= 1 && arq.condpocinj.CC != 3 && arq.condpocinj.CC != 4 && arq.condpocinj.CC != 5) {
        int testaIPR = 0;
        for (int i = 0; i < arq.nipr; i++)
            if (arq.IPRS[i].indcel == arq.ncelp - 1)
                testaIPR = 1;
        if (testaIPR == 0)
            NumError(
                "O simulador esta no modo Injecao de Agua em uma condicao de contorno que pede uma IPR e a ultima celula nao contem uma IPR, e necessario neste modo se ter uma IPR na ultima celula");
    }

    if (arq.modelcp > 0)
        arq.geraTabCp();
    if (arq.modelJTL == 1)
        arq.geraTabDrholDt();
    tables.gasSpecificHeatTable = arq.cpg;
    tables.liquidSpecificHeatTable = arq.cpl;
    tables.liquidDensityTemperatureDerivative = arq.drholdT;

    tables.gasCompressibilityTable = arq.zdranP;
    tables.gasCompressibilityPressureDerivative = arq.dzdpP;
    tables.gasCompressibilityTemperatureDerivative = arq.dzdtP;
    steadySearch.reverseNetworkFluid = arq.flup[0];
    if (arq.flashCompleto == 2)
        steadySearch.reverseNetworkFluid.atualizaPropCompStandard();
    for (int i = 0; i < arq.nfluP; i++) {
        arq.flup[i].npontos = arq.tabent.npont;
        arq.flup[i].cpg = tables.gasSpecificHeatTable;
        arq.flup[i].cpl = tables.liquidSpecificHeatTable;
        arq.flup[i].drholdT = tables.liquidDensityTemperatureDerivative;
        arq.flup[i].nfluP = steadySearch.productionFluidCount;
    }
    arq.flug.npontos = arq.tabent.npont;
    arq.flug.cpg = tables.gasSpecificHeatTable;
    arq.flug.cpl = tables.liquidSpecificHeatTable;
    arq.flug.drholdT = tables.liquidDensityTemperatureDerivative;
    arq.flug.nfluP = steadySearch.productionFluidCount;

    arq.fluc.npontos = arq.tabent.npont;
    arq.fluc.RhoInj = arq.RhoInj;
    arq.fluc.ViscInj = arq.ViscInj;
    arq.fluc.CondInj = arq.CondInj;
    arq.fluc.CpInj = arq.CpInj;
    arq.fluc.DrhoDtInj = arq.DrhoDtInj;

    if (arq.modelcp > 0 || arq.modelJTL == 1 || arq.latente > 0 || (arq.pocinjec == 1 && (arq.condpocinj.tipoFlui == 2 || arq.condpocinj.tipoFlui == 3)))
        transient.massTransferModel = arq.transmass;
    arq.geraduto();
    ncel = arq.ncelp - 1;
    temperatura = arq.celp[0].textern;
    steadySearch.reverseFlowTemperature = arq.tempReves;
    productionCells = vector<Cel>(ncel + 1);
    celula = productionCells.data();
    arq.geracelp(celula);
    if (arq.nipr > 0)
        arq.geraipr(celula);
    if (arq.nvalvgas > 0 && arq.lingas > 0)
        arq.gerafgasVGL(celula);
    if (arq.ninjgas > 0)
        arq.gerafgasFonte(celula);
    for (int i = arq.ninjgas; i < arq.ninjgas + extraSourceCount; i++) {
        int iposp = arq.buscaIndiceMeioP(extraSourceLengths[i - arq.ninjgas]);
        extraSourceCells[i - arq.ninjgas] = iposp;
        InjGas injgasMRT(0, 0, arq.flug);
        celula[iposp].acsr.tipo = kAccessoryGasInjection;
        celula[iposp].acsr.injg = injgasMRT;
    }
    if (arq.ninjliq > 0)
        arq.gerafliqFonte(celula);
    if (arq.ninjmass > 0)
        arq.gerafmassFonte(celula);
    if (arq.nfuro > 0)
        arq.geraFuro(celula);
    if (arq.nPoroRad > 0)
        arq.gerafPoroRadFonte(celula);
    if (arq.nPoro2D > 0)
        arq.gerafPoro2DFonte(celula);
}

namespace {

/// Gives `cell` a multiple, a liquid and a gas source, all with zero flow,
/// carrying the fluids of `inlet` (and, for the multiple source, its
/// temperature).
void attachEmptySources(Cel &cell, Cel &inlet) {
    cell.acsr.injm = InjMult(0, 0, 0, inlet.temp, inlet.flui, inlet.fluicol);
    cell.acsr.injl = InjLiq(0, 0, 0, inlet.flui, inlet.fluicol);
    cell.acsr.injg = InjGas(0, 0, inlet.flui, inlet.fluicol);
}

}  // namespace

/// Gives the inlet the sources its boundary condition needs (and the second cell,
/// under a blockage), then places the accessories -- pumps, volumetric pumps,
/// pressure-drop requirements, heat sources, the master valve and the other
/// valves -- and sets up the outlet pressure, the surface and injection chokes and
/// the pigs.
void SProd::configureInletSourcesAndAccessories(int extraSourceCount) {
    if (celula[0].acsr.tipo != kAccessoryNone && arq.ConContEntrada == 1) {
        // RN-302: The inlet has both a pressure boundary condition and a mass source. Report a warning.
        logger.log(LOGGER_AVISO, LOG_ERR_PARSE_BUSINESS_RULE_VALIDATION,
                   "Foi escolhida uma condicao de contorno de pressao no inicio da tubulacao e foi colocada uma fonte de massa neste inicio, isto pode causar inconsistencias",
                   "", "");
    }
    if ((arq.perm == 1 || (*vg1dSP).chaverede > 0) && celula[0].acsr.tipo == kAccessoryNone) {
        if ((arq.ConContEntrada == 0 && extraSourceCount == 0 && (*vg1dSP).chaverede == 0)) {
            NumError("O simulador pede calculo de permanente com condicao de vazao na entrada, mas sem nenhuma fonte na entrada");
        }
        if ((arq.ConContEntrada == 0 && extraSourceCount == 0 && (*vg1dSP).chaverede > 0) ||
            arq.CCPres.tit[0] < (1 - (*vg1dSP).localtiny)) {
            if (arq.tipoFluido == 1)
                celula[0].acsr.tipo = kAccessoryGasInjection;
            else if (arq.tipoFluido == 0)
                celula[0].acsr.tipo = kAccessoryLiquidInjection;
            else
                celula[0].acsr.tipo = kAccessoryMultipleSource;
            attachEmptySources(celula[0], celula[0]);
        } else if (extraSourceCount == 0 || arq.ConContEntrada == 1) {
            attachEmptySources(celula[0], celula[0]);
            celula[0].acsr.injg.seco = 0;
            if (arq.tipoFluido == 1)
                celula[0].acsr.tipo = kAccessoryGasInjection;
            else if (arq.tipoFluido == 0)
                celula[0].acsr.tipo = kAccessoryLiquidInjection;
            else
                celula[0].acsr.tipo = kAccessoryMultipleSource;
        }
    }
    if (bloq == 1 && celula[1].acsr.tipo == kAccessoryNone) {
        if ((*vg1dSP).fluidoRede == 1)
            celula[1].acsr.tipo = kAccessoryLiquidInjection;
        else if ((*vg1dSP).fluidoRede == 0)
            celula[1].acsr.tipo = kAccessoryGasInjection;
        else
            celula[1].acsr.tipo = kAccessoryMultipleSource;
        attachEmptySources(celula[1], celula[0]);
    } else if (bloq == 1) {
        logger.log(LOGGER_AVISO, LOG_ERR_PARSE_BUSINESS_RULE_VALIDATION,
                   "Foi escolhida uma condicao de bloqueio que afeta este tramo, para que a selecao de fontes seja"
                   "feita corretamente, e necessario que a segunda celula do tramo nao tenha fontes catalogadas",
                   "", "");
    }
    if (arq.nbcs > 0)
        arq.gerafBCS(celula);
    if (arq.nmultibcs > 0)
        arq.gerafmultiBCS(celula);
    if (arq.nbvol > 0)
        arq.gerafBVOL(celula);
    if (arq.ndpreq > 0)
        arq.geraDPReq(celula);
    if (arq.ncalor > 0)
        arq.geraFonteCalor(celula);
    if (arq.master1.posic > 0.)
        arq.geraMaster1(celula);
    else {
        int verifica = 0;
        double Lverifica = 0.5 * celula[verifica].dx;
        while (celula[verifica].acsr.tipo != kAccessoryNone) {
            verifica++;
            Lverifica += 0.5 * (celula[verifica].dx + celula[verifica - 1].dx);
        }
        arq.master1.posic = verifica;
        arq.master1.comp = Lverifica;
        arq.geraMaster1(celula);
    }
    if (arq.nvalv > 0)
        arq.geraValv(celula);
    arq.gerapresfim(presfim, gasLift.separatorPressure);
    transient.previousSeparatorPressure = gasLift.separatorPressure;
    chokeSup = choke(1., 1.);
    gasLift.injectionChoke = ChokeGas();
    arq.gerachokesup(chokeSup);
    if (arq.npig > 0) {
        transient.pigReceiverCells = vector<int>(arq.npig);
        for (int i = 0; i < arq.npig; i++)
            transient.pigReceiverCells[i] = arq.pig[i].receb;
    }
}

/// Builds the gas-lift line when there is one: its cells, the second master valve,
/// the injection choke, one gas-lift valve choke per valve with its position on
/// both lines, and each gas cell's share of the annulus above the discharge cell.
/// Without a gas line, only zeroes the gas cell count.
void SProd::buildGasLiftLine() {
    if (arq.lingas == 0)
        gasLift.cellCount = 0;
    else if (arq.lingas > 0) {
        gasLift.cellCount = arq.ncelg;
        gasLift.cellCount--;
        gasLift.interfaceCell = arq.celdescarga;
        gasLift.gasCells = vector<CelG>(gasLift.cellCount + 1);
        gasLift.cells = gasLift.gasCells.data();
        arq.geracelg(gasLift.cells);
        arq.geraMaster2(gasLift.cells);
        arq.gerachokeinj(gasLift.injectionChoke);
        if (arq.nvalvgas > 0) {
            gasLift.valveChokes = vector<ChokeGas>(arq.nvalvgas);
            gasLift.productionValveCells = vector<int>(arq.nvalvgas);
            gasLift.gasValveCells = vector<int>(arq.nvalvgas);
        }
        for (int i = 0; i < arq.nvalvgas; i++) {
            double diaG = arq.valvgl[i].diagarg;
            double presEstag = gasLift.cells[arq.valvgl[i].posicG].pres;
            double tempEstag = gasLift.cells[arq.valvgl[i].posicG].temp;
            double presGarg = celula[arq.valvgl[i].posicP].pres;
            gasLift.productionValveCells[i] = arq.valvgl[i].posicP;
            gasLift.gasValveCells[i] = arq.valvgl[i].posicG;
            gasLift.cells[gasLift.gasValveCells[i]].vgl = 1;
            if (arq.valvgl[i].tipo == 2) {
                arq.valvgl[i].frec = 0.;
                arq.valvgl[i].cd = 1.;
            }
            double throatArea = std::numbers::pi * diaG * diaG / 4.;
            gasLift.valveChokes[i] = ChokeGas(arq.flug, throatArea,
                                   arq.valvgl[i].diaexter, arq.valvgl[i].cd, presEstag, presGarg,
                                   tempEstag, arq.valvgl[i].frec, arq.valvgl[i].tipo,
                                   throatArea / arq.valvgl[i].razarea,
                                   arq.valvgl[i].pcali, arq.valvgl[i].tcali, arq.valvgl[i].cdLiq, arq.valvgl[i].frecLiq);
            gasLift.cells[gasLift.gasValveCells[i]].pEstag = gasLift.valveChokes[i].presEstag;
            gasLift.cells[gasLift.gasValveCells[i]].tEstag = gasLift.valveChokes[i].tempEstag;
            gasLift.cells[gasLift.gasValveCells[i]].pGarg = gasLift.valveChokes[i].presGarg;
            gasLift.cells[gasLift.gasValveCells[i]].tGarg = gasLift.valveChokes[i].tempGarg;
            gasLift.cells[gasLift.gasValveCells[i]].qGarg = gasLift.valveChokes[i].qGarg;
            gasLift.cells[gasLift.gasValveCells[i]].areaGarg = gasLift.valveChokes[i].areagarg;
        }
        for (int i = 0; i <= gasLift.cellCount; i++) {
            gasLift.cells[i].celInter = &gasLift.interfaceCell;
            gasLift.cells[i].razInter = 1.;
            gasLift.cells[i].razInterIni = 1.;
            if (i >= gasLift.interfaceCell) {
                gasLift.cells[i].razInter = 0.;
                gasLift.cells[i].razInterIni = 0.;
            }
        }
    }
}

/// Rejects or warns about source and boundary-condition combinations the run
/// cannot honour (RN-300, RN-301, gas-lift discharge without an IPR, no outlet
/// pressure, accessories in the last two cells), builds the gas-lift discharge
/// hydrostatics and the event log, applies the initial state of a production
/// well, and records the surface temperature and mass flow it starts from.
void SProd::validateSetupAndApplyInitialState() {
    if (celula[0].acsr.tipo != kAccessoryGasInjection && celula[0].acsr.tipo != kAccessoryLiquidInjection && celula[0].acsr.tipo != kAccessoryInflowPerformance && celula[0].acsr.tipo != kAccessoryMultipleSource && celula[0].acsr.tipo != kAccessoryRadialPorous && celula[0].acsr.tipo != kAccessoryPorous2D) {
        if (arq.perm == 0 && arq.ConContEntrada == 0) {
            // RN-300: No source is defined in the first production system cell. Report a warning.
            logger.log(LOGGER_AVISO, LOG_ERR_PARSE_BUSINESS_RULE_VALIDATION,
                       "Nao existe nenhuma fonte na primeira celula do sistema de producao, esta nao e a condicao mais adequada para a simulacao transiente, problemas podem ocorrer nesta situacao, embora nao seja impeditivo para o modo",
                       "", "");
        } else {
            // RN-301: No source is defined in the first production system cell
            // under steady-state conditions. Report a failure.
            NumError("Nao existe nenhuma fonte na primeira celula do sistema de producao, no modo permanente esta e uma condicao necessaria");
        }
    }

    if (arq.descarga == 1 && celula[0].acsr.tipo != kAccessoryInflowPerformance) {
        NumError(
            "Nao existe nenhuma IPR na primeira celula do sistema de producao, no modo descarga de gas lift esta e uma condicao necessaria");
    }
    if (arq.descarga == 1) {
        HidroDescargaG();
        HidroDescargaP();
    }

    arq.geraevento(noinicial, noextremo);
    if (arq.saidaClassica == 1 && (arq.transiente == 1 || (*vg1dSP).chaveredeT == 1)) {
        cout << "------------------------------------------   Eventos no Tramo -------------------------------------------" << endl;
        int contaImp = 0;
        for (int i = 0; i < arq.logevento.size(); i++)
            if (arq.logevento[i].descricao != "Gravando Perfil Linha de Producao" && arq.logevento[i].descricao != "Gravando Perfil Linha de Gas") {
                cout << contaImp << "  momento = " << arq.logevento[i].instante << " segundos " << " tipo de evento = " << arq.logevento[i].descricao << endl;
                contaImp++;
            }
    }

    if (arq.pocinjec == 0 && arq.perm == 1 && ((*vg1dSP).chaverede == 0 || noextremo == 1) && arq.ConContEntrada != 2 && arq.psep.pres[0] < -1e6) {
        NumError("sem pressao a jusante como condição de contorno para um tramo onde esta condicao é necessaria ");
    }

    if ((celula[ncel].acsr.tipo != kAccessoryNone || celula[ncel - 1].acsr.tipo != kAccessoryNone) && arq.transiente == 1) {
        NumError("As duas celulas finais não devem ter acessorios, estas duas celulas devem ser preservadas, principalmente para simulacoes transientes ");
    }

    (*vg1dSP).lixo5 = 0.;
    (*vg1dSP).contador = 0;
    if (arq.pocinjec == 0) {
        arq.atualiza(noinicial, noextremo, derivaAnel, chokeSup, gasLift.injectionChoke, celula, gasLift.cells,
                     gasLift.separatorPressure, temperatura, gasLift.injectionPressure, gasLift.injectionTemperature,
                     presE, tempE, titE, betaE, (*vg1dSP).lixo5, dt);
        transient.previousSeparatorPressure = gasLift.separatorPressure;
        if (chokeSup.AreaGarg >= 0.6 * celula[ncel - 1].duto.area) {
            transient.surfaceChokeOpen = 1;
            transient.previousSurfaceChokeOpen = 1;
        } else {
            transient.surfaceChokeOpen = 0;
            transient.previousSurfaceChokeOpen = 1;
        }
    }

    networkCoupling.chokeOutletTemperature = celula[ncel].temp;

    tGSup = celula[ncel].calor.Textern1;
}

/// With the dynamic property table on, splits the pipe into table segments that
/// end at every source cell. Then gives each cell its left and right inclination:
/// its own, or, where its pipe is horizontal, the nearest inclined neighbour's on
/// that side.
void SProd::buildDynamicTablesAndInclinations() {
    if (arq.tabelaDinamica == 1) {
        tables.dynamicTables.push_back(tabelaDinamica());
        for (int i = 1; i < ncel; i++) {
            if ((celula[i].acsr.tipo == kAccessoryGasInjection) ||
                (celula[i].acsr.tipo == kAccessoryLiquidInjection) ||
                celula[i].acsr.tipo == kAccessoryInflowPerformance ||
                (celula[i].acsr.tipo == kAccessoryMultipleSource) || celula[i].acsr.tipo == kAccessoryLeak || celula[i].acsr.tipo == kAccessoryRadialPorous || celula[i].acsr.tipo == kAccessoryPorous2D) {
                tabelaDinamica temp;
                temp.celIni = i + 1;
                tables.dynamicTables.back().celFim = i;
                tables.dynamicTables.push_back(std::move(temp));
            }
        }
        tables.dynamicTables.back().celFim = ncel;
    }

    celula[0].angEsq = celula[0].duto.teta;
    for (int i = 1; i <= ncel; i++) {
        if (fabs(celula[i].duto.teta) < 1e-10)
            celula[i].angEsq = celula[i - 1].angEsq;
        else
            celula[i].angEsq = celula[i].duto.teta;
    }
    celula[ncel].angEsq = celula[ncel].duto.teta;
    for (int i = ncel - 1; i >= 0; i--) {
        if (fabs(celula[i].duto.teta) < 1e-10)
            celula[i].angDir = celula[i + 1].angDir;
        else
            celula[i].angDir = celula[i].duto.teta;
    }
}

/// Sets the latent-heat switch from the input and, when it is on and flashCompleto
/// is 0, reads the latent-heat table from the PVTSim file into HLat and writes
/// perfilLatente.
void SProd::configureLatentHeat() {
    CalcLat = arq.latente;
    if (arq.latente > 0) {
        if (arq.flashCompleto == 0) {
            string impfile;
            impfile = arq.pvtsimarq;
            string dadosMR = impfile;
            ifstream lendoPVTSim(dadosMR.c_str(), ios_base::in);
            string chave;
            char *tenta;
            double testatok;
            int ndiv = arq.tabent.npont - 1;
            Vcr<double> presPVTSim(ndiv + 1, 0.);
            Vcr<double> tempPVTSim(ndiv + 1, 0.);

            FullMtx<double> HLatTemp(ndiv + 2, ndiv + 2);
            char line[4000];
            lendoPVTSim.get(line, 4000);
            tenta = strtok(line, " ,()=");
            while (strcmp(tenta, "PHASE") != 0) {
                tenta = strtok(NULL, " ,()=");
            }
            tenta = strtok(NULL, " ,()=");
            int lacoleitura = 12;
            if (strcmp(tenta, "THREE") == 0)
                lacoleitura = 18;

            vector<double> pressuresRead = readPvtSimRow(lendoPVTSim, chave, "PRESSURE", dadosMR,
                                                         [](double value) { return value / units::kPascalPerKgfPerCm2; });
            vector<double> temperaturesRead = readPvtSimRow(lendoPVTSim, chave, "TEMPERATURE", dadosMR,
                                                            [](double value) { return value; });
            if (static_cast<int>(pressuresRead.size()) <= ndiv || static_cast<int>(temperaturesRead.size()) <= ndiv)
                logger.log_write_logs_and_exit(LOGGER_FALHA, LOG_ERR_PARSE_BUSINESS_RULE_VALIDATION,
                                               "PVTSim file with fewer table points than the input table", "pvtsimArq",
                                               dadosMR);
            for (int k = 0; k <= ndiv; k++) {
                presPVTSim[k] = pressuresRead[k];
                tempPVTSim[k] = temperaturesRead[k];
            }

            for (int i = 1; i <= ndiv + 1; i++) {
                HLatTemp[i][0] = presPVTSim[i - 1];
                for (int j = 1; j <= ndiv + 1; j++) {
                    HLatTemp[0][j] = tempPVTSim[j - 1];
                    while (chave != "POINT")
                        lendoPVTSim >> chave;
                    lendoPVTSim.get(line, 4000);
                    tenta = strtok(line, " ,()=");
                    for (int kontaPVT = 0; kontaPVT < lacoleitura; kontaPVT++) {
                        tenta = strtok(NULL, " ,");
                        testatok = atof(tenta);
                    }
                    tenta = strtok(NULL, " ,");
                    testatok = atof(tenta);
                    HLatTemp[i][j] = testatok;
                    tenta = strtok(NULL, " ,");
                    testatok = atof(tenta);
                    HLatTemp[i][j] -= testatok;
                    if (i == ndiv + 1 && j == ndiv + 1)
                        break;
                    while (chave != "PVTTABLE")
                        lendoPVTSim >> chave;
                }
            }
            lendoPVTSim.close();
            tables.latentHeatTable.allocate(ndiv + 2, ndiv + 2);
            for (int i = 0; i < ndiv + 2; i++) {
                for (int j = 0; j < ndiv + 2; j++)
                    tables.latentHeatTable[i][j] = HLatTemp[i][j];
            }
            string tmp = pathPrefixoArqSaida + "perfilLatente.dat";
            writeTable(tmp, HLatTemp);
            // unless this is an injection-well (POCO_INJETOR) simulation
            if (arq.tipoSimulacao != tipoSimulacao_t::poco_injetor) {
                arqRelatorioPerfis << tmp.c_str() << endl;
                arqRelatorioPerfis.flush();
            }
        }
    }
}

/// Sets the gas-density correction factors of every fluid in every cell -- the
/// cell's, its source's and its reservoir cells' -- to 1 when the correction is
/// off, or evaluates them at local pressure and temperature when it is on; then
/// makes the first cell's fluid the fluid of its source.
void SProd::applyDensityCorrectionsAndInletFluid() {
    for (int i = 0; i <= ncel; i++) {
        if (arq.corrDeng == 0) {
            celula[i].flui.rDgD = 1.;
            celula[i].flui.rDgL = 1.;
            celula[i].flui.PCis = celula[i].flui.PC;
            celula[i].flui.TCis = celula[i].flui.TC;
            if (celula[i].acsr.tipo == kAccessoryLiquidInjection) {
                celula[i].acsr.injl.FluidoPro.rDgD = 1.;
                celula[i].acsr.injl.FluidoPro.rDgL = 1.;
            } else if (celula[i].acsr.tipo == kAccessoryInflowPerformance) {
                celula[i].acsr.ipr.FluidoPro.rDgD = 1.;
                celula[i].acsr.ipr.FluidoPro.rDgL = 1.;
            } else if (celula[i].acsr.tipo == kAccessoryMultipleSource) {
                celula[i].acsr.injm.FluidoPro.rDgD = 1.;
                celula[i].acsr.injm.FluidoPro.rDgL = 1.;
            }
            if (celula[i].acsr.tipo == kAccessoryRadialPorous) {
                celula[i].acsr.radialPoro.flup.rDgD = 1.;
                celula[i].acsr.radialPoro.flup.rDgL = 1.;
                for (int iRP = 0; iRP < celula[i].acsr.radialPoro.ncel; iRP++) {
                    celula[i].acsr.radialPoro.celula[iRP].flup.rDgD = 1.;
                    celula[i].acsr.radialPoro.celula[iRP].flup.rDgL = 1.;
                }
            }
            if (celula[i].acsr.tipo == kAccessoryPorous2D) {
                celula[i].acsr.poroso2D.dados.flup.rDgD = 1.;
                celula[i].acsr.poroso2D.dados.flup.rDgL = 1.;
                for (int iRP = 0; iRP < celula[i].acsr.poroso2D.dados.transfer.ncel; iRP++) {
                    celula[i].acsr.poroso2D.dados.transfer.celula[iRP].flup.rDgD = 1.;
                    celula[i].acsr.poroso2D.dados.transfer.celula[iRP].flup.rDgL = 1.;
                }
                for (int iRP = 0; iRP < celula[i].acsr.poroso2D.malha.nele; iRP++) {
                    celula[i].acsr.poroso2D.malha.mlh2d[iRP].flup.rDgD = 1.;
                    celula[i].acsr.poroso2D.malha.mlh2d[iRP].flup.rDgL = 1.;
                }
            }
        } else {
            celula[i].flui.razDegD(celula[i].pres, celula[i].temp);
            celula[i].flui.rzDegL(celula[i].pres, celula[i].temp);
            celula[i].flui.PcTcIS();
            forEachSourceFluid(celula[i], [](ProFlu &fluid, double pressure, double temperature) {
                fluid.razDegD(pressure, temperature);
                fluid.rzDegL(pressure, temperature);
            });
        }
    }

    if (celula[0].acsr.tipo == kAccessoryGasInjection) {
        celula[0].flui = celula[0].acsr.injg.FluidoPro;
    } else if (celula[0].acsr.tipo == kAccessoryLiquidInjection) {
        celula[0].flui = celula[0].acsr.injl.FluidoPro;
    } else if (celula[0].acsr.tipo == kAccessoryInflowPerformance) {
        celula[0].flui = celula[0].acsr.ipr.FluidoPro;
    } else if (celula[0].acsr.tipo == kAccessoryMultipleSource) {
        celula[0].flui = celula[0].acsr.injm.FluidoPro;
    } else if (celula[0].acsr.tipo == kAccessoryRadialPorous) {
        celula[0].flui = celula[0].acsr.radialPoro.flup;
    } else if (celula[0].acsr.tipo == kAccessoryPorous2D) {
        celula[0].flui = celula[0].acsr.poroso2D.dados.flup;
    }
}

/// Lists the cells with two-dimensional heat diffusion, copies the master-valve
/// opening and closing times, maps the transient profile positions to global
/// thermal-node indices, and -- unless the network is only provisional --
/// allocates the trend matrices of the production line, the gas line and the
/// transient trends, sized to the longest simulated time and filled with the
/// -10000 sentinel.
void SProd::allocateEventProfileAndTrendArrays() {
    for (int i = 0; i < ncel; i++) {
        if (celula[i].calor.difus2D == 1) {
            transient.poisson2DCells.push_back(i);
        }
    }

    transient.masterValveOpenTimes = vector<double>(arq.eventoabre);
    transient.masterValveCloseTimes = vector<double>(arq.eventofecha);
    for (int i = 0; i < arq.eventoabre; i++)
        transient.masterValveOpenTimes[i] = arq.Tevento[i];
    for (int i = 0; i < arq.eventofecha; i++)
        transient.masterValveCloseTimes[i] = arq.Teventof[i];

    int ntempGas = 0;
    if (arq.lingas > 0)
        ntempGas = gasLift.cellCount;

    if (arq.nperfistransp > 0) {
        trends.productionWallProfileCells = vector<int>(arq.nperfistransp);
        for (int i = 0; i < arq.nperfistransp; i++) {
            int posiccel = arq.proftransp.posic[i];
            trends.productionWallProfileCells[i] = celula[posiccel].calor.nglobal;
        }
    }
    if (arq.nperfistransg > 0 && arq.lingas > 0) {
        trends.gasWallProfileCells = vector<int>(arq.nperfistransg);
        for (int i = 0; i < arq.nperfistransg; i++) {
            int posiccel = arq.proftransg.posic[i];
            trends.gasWallProfileCells[i] = gasLift.cells[posiccel].calor.nglobal;
        }
    }

    if ((*vg1dSP).chaverede == 0)
        (*vg1dSP).TmaxR = arq.tfinal;
    if (arq.ntendp > 0) {
        trends.productionTrendSet.allocate(
            arq.ntendp,
            [&](int i) -> int { return 1 + 1 + ceil((*vg1dSP).TmaxR / arq.trendp[i].dt); }, // round(arq.tfinal / arq.trendp[i].dt);
            [&](int i) { return arq.nvartrendp[i] + 2; }, [&](int i) { return arq.nvartrendp[i] + 1; },
            trends.productionTrendMatrices, trends.productionTrendResetTimes, trends.productionTrendCounts,
            trends.productionTrendBufferedCounts);
    }
    if (arq.ntendg > 0) {
        trends.gasTrendSet.allocate(
            arq.ntendg, [&](int i) -> int { return 1 + 1 + ceil((*vg1dSP).TmaxR / arq.trendg[i].dt); },
            [&](int i) { return arq.nvartrendg[i] + 2; }, [&](int i) { return arq.nvartrendg[i] + 1; },
            trends.gasTrendMatrices, trends.gasTrendResetTimes, trends.gasTrendCounts, trends.gasTrendBufferedCounts);
    }
    if (arq.ntendtransp > 0) {
        trends.productionWallTrendSet.allocate(
            arq.ntendtransp, [&](int i) -> int { return 1 + 1 + ceil((*vg1dSP).TmaxR / arq.trendtransp[i].dt); },
            [](int) { return 2; }, [](int) { return 2; }, trends.productionWallTrendMatrices,
            trends.productionWallTrendResetTimes, trends.productionWallTrendCounts,
            trends.productionWallTrendBufferedCounts);
    }
    if (arq.ntendtransg > 0) {
        trends.gasWallTrendSet.allocate(
            arq.ntendtransg, [&](int i) -> int { return 1 + 1 + ceil((*vg1dSP).TmaxR / arq.trendtransg[i].dt); },
            [](int) { return 2; }, [](int) { return 2; }, trends.gasWallTrendMatrices, trends.gasWallTrendResetTimes,
            trends.gasWallTrendCounts, trends.gasWallTrendBufferedCounts);
    }
}

/// Resets the column-annulus and network coupling flags, locates the
/// column-annulus coupling range on the gas line, sets the steady and transient
/// profile counters and their first output times, opens the event log with the
/// events known at start, finds the smallest cell length, and zeroes the
/// moving-average and running-total state the transient loop starts from.
void SProd::resetCouplingAndOutputState() {
    gasLift.annulusThermalCoupling = 0;

    if (arq.lingas > 0) {
        gasLift.annulusCouplingStart = arq.anulcoluini();
        gasLift.annulusCouplingEnd = arq.anulcolufim();
        gasLift.tubingCouplingStart = arq.coluanulini();
        gasLift.tubingCouplingEnd = arq.coluanulfim();
        if (gasLift.tubingCouplingEnd == -1) {
            gasLift.tubingCouplingEnd = 0;
            gasLift.annulusCouplingEnd--;
        }
        if (gasLift.annulusCouplingStart >= 0 && gasLift.annulusCouplingEnd >= 0 && gasLift.tubingCouplingStart >= 0 &&
            gasLift.tubingCouplingEnd >= 0)
            gasLift.annulusThermalCoupling = 1;
    }

    networkCoupling.primaryThermalCoupling = 0;
    networkCoupling.secondaryThermalCoupling = 0;
    networkCoupling.secondaryCouplingStart = 0;
    networkCoupling.secondaryCouplingEnd = 0;
    networkCoupling.primaryCouplingStart = 0;
    networkCoupling.primaryCouplingEnd = 0;

    if(arq.AP==0){
    	if (arq.nperfisp > 0) {
    		arq.imprimeProfile(celula, flut, (*vg1dSP).lixo5, indTramo);
    	}
    	if (arq.nperfisg > 0 && arq.lingas > 0) {
    		arq.imprimeProfileG(gasLift.cells, flutG, (*vg1dSP).lixo5, indTramo);
    	}
    	if (arq.nperfistransp > 0) {
    		arq.imprimeProfileTrans(celula, trends.productionWallProfileCells.data(), (*vg1dSP).lixo5, indTramo);
    	}
    	if (arq.nperfistransg > 0 && arq.lingas > 0) {
    		arq.imprimeProfileTransG(gasLift.cells, trends.gasWallProfileCells.data(), (*vg1dSP).lixo5, indTramo);
    	}
    }
    trends.productionProfileTimeIndex = 0;
    trends.gasProfileTimeIndex = 0;
    if (arq.nperfisp > 0) {
        if (arq.profp.tempo[0] <= 0 + (*vg1dSP).localtiny)
            trends.productionProfileTimeIndex++;
    }
    if (arq.nperfisg > 0 && arq.lingas > 0) {
        if (arq.profg.tempo[0] <= 0 + (*vg1dSP).localtiny)
            trends.gasProfileTimeIndex++;
    }
    trends.productionWallProfileTimeIndex = 0;
    trends.gasWallProfileTimeIndex = 0;
    if (arq.nperfistransp > 0) {
        if (arq.proftransp.tempo[0] <= 0 + (*vg1dSP).localtiny)
            trends.productionWallProfileTimeIndex++;
    }
    if (arq.nperfistransg > 0 && arq.lingas > 0) {
        if (arq.proftransg.tempo[0] <= 0 + (*vg1dSP).localtiny)
            trends.gasWallProfileTimeIndex++;
    }

    ostringstream saidaLog;
    saidaLog << pathPrefixoArqSaida << "LogEvento" << ".dat";
    tmpLog = saidaLog.str();
    // if it's not a simulation POCO_INJETOR
    if (arq.tipoSimulacao != tipoSimulacao_t::poco_injetor) {
        transient.eventLogCount = 0;
        int nevent = arq.logevento.size();
        while (transient.eventLogCount < nevent &&
               fabs(arq.logevento[transient.eventLogCount].instante - (*vg1dSP).lixo5) < dt) {
            time_t now = time(0);
            tm *ltm = localtime(&now); // Taken from https://www.tutorialspoint.com/cplusplus/cpp_date_time.htm
            ofstream escreveIni(tmpLog.c_str(), ios_base::app);
            escreveIni << "Evento Externo = ";
            escreveIni << arq.logevento[transient.eventLogCount].instante << " ; ";
            escreveIni << arq.logevento[transient.eventLogCount].duracao << " ; ";
            escreveIni << arq.logevento[transient.eventLogCount].descricao << " ; ";
            escreveIni << "datahora = ";
            escreveIni << ltm->tm_mday << "/";
            escreveIni << 1 + ltm->tm_mon << "/";
            escreveIni << 1900 + ltm->tm_year << " ";
            escreveIni << 0 + ltm->tm_hour << ":";
            escreveIni << 0 + ltm->tm_min << ":";
            escreveIni << 0 + ltm->tm_sec;
            escreveIni << endl;
            transient.eventLogCount++;
            escreveIni.close();
        }
    }

    transient.smallestCellLength = 1e10;
    for (int i = 0; i <= ncel; i++) {
        if (celula[i].dx < transient.smallestCellLength)
            transient.smallestCellLength = celula[i].dx;
        if (i > 0)
            celula[i].razdxTM = celula[i - 1].dx / (celula[i - 1].dx + celula[i].dx);
        if (i > 1)
            celula[i].razdxTM0 = celula[i - 2].dx / (celula[i - 1].dx + celula[i - 2].dx);
    }

    if (arq.nCelUnit>0) {
    	for(int iCelU=0;iCelU<arq.nCelUnit;iCelU++){
    		trends.singleCellOutputTimeIndices.push_back(0);
    	}
    }

    transient.stepCount = 0;
    transient.eventIndex = 1;
    transient.cflSafetyFactor = 0.8;
    transient.movingMeanPressure = 0.;
    transient.movingMeanFlux = 0.;
    transient.movingMeanStartTime = 60.;
    transient.movingMeanWindow = 0.;
    transient.accumulatedPressure = 0.;
    transient.accumulatedFlux = 0.;
    transient.accumulatedVoidFraction = 0.;
}

void SProd::montasistema(double *extraSourceLengths, int *extraSourceCells, int extraSourceCount) {
    try {
        buildProductionCells(extraSourceLengths, extraSourceCells, extraSourceCount);
        configureInletSourcesAndAccessories(extraSourceCount);
        buildGasLiftLine();
        validateSetupAndApplyInitialState();
        buildDynamicTablesAndInclinations();
        configureLatentHeat();
        if (celula[0].flui.awaitingSaturationTable) {
            loadPvtSimSaturationTables();
        } else if (arq.tabRSPB == 1) {
            generateSaturationTablesFromCorrelations();
        }
        applyDensityCorrectionsAndInletFluid();
        allocateEventProfileAndTrendArrays();
        resetCouplingAndOutputState();
    } catch (exception &excInt) {
        cout << "EXCECAO INESPERADA: " << excInt.what() << endl;
        // report the failure
        logger.log(LOGGER_FALHA, LOG_ERR_UNEXPECTED_EXCEPTION, "", "", excInt.what());
        // write the log file
        logger.writeOutputLog();
        // end the application
        exit(EXIT_SUCCESS);
    }
}
