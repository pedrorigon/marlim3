/*
 * JSONKeyTranslator.cpp
 *
 * Implements the English→Portuguese key translation for Marlim3 JSON input
 * files. The translation is a flat map: every occurrence of an English key
 * at any nesting level is renamed to the corresponding Portuguese key before
 * the typed JSON_entrada schema processes the document.
 *
 * Key design decisions
 * --------------------
 * - Flat map: the same English name translates to the same Portuguese name
 *   regardless of nesting level. This is safe because Portuguese sibling keys
 *   within a single object are always distinct.
 * - Keys that are already language-neutral (acronyms such as "api", "bsw",
 *   "bcs", "pig", "ipr") are kept identical in both languages.
 * - Keys added here but absent from the document are silently ignored.
 * - Unknown English keys (not in the map) are left as-is; the schema will
 *   ignore them just as it would ignore unknown Portuguese keys.
 */

#include "JSONKeyTranslator.h"
#include "rapidjson/document.h"

#include <algorithm>
#include <string>
#include <unordered_map>
#include <vector>
#include <utility>

using namespace rapidjson;

// ---------------------------------------------------------------------------
// Static English → Portuguese key map
// ---------------------------------------------------------------------------
// Format: { "englishKey", "portugueseKey" }
// Keys that are identical in both languages are omitted.
// ---------------------------------------------------------------------------
static const std::unordered_map<std::string, std::string> EN_TO_PT = {

    // -----------------------------------------------------------------------
    // Root-level section names
    // -----------------------------------------------------------------------
    { "system",                         "sistema" },
    { "jsonVersion",                    "versaoJSON" },
    { "initialConfig",                  "configuracaoInicial" },
    { "timeSettings",                   "tempo" },
    { "pvtTable",                       "tabela" },
    { "wax",                            "parafina" },
    { "gasFluid",                       "fluidoGas" },
    { "productionFluids",               "fluidosProducao" },
    { "complementaryFluid",             "fluidoComplementar" },
    { "valve",                          "valvula" },
    { "liquidSource",                   "fonteLiquido" },
    { "massSource",                     "fonteMassa" },
    { "gasSource",                      "fonteGas" },
    { "porousRadialSource",             "fontePoroRadial" },
    { "porous2DSource",                 "fontePoro2D" },
    { "gasLiftSource",                  "fonteGasLift" },
    { "crossSection",                   "secaoTransversal" },
    { "productionDucts",                "dutosProducao" },
    { "serviceDucts",                   "dutosServico" },
    { "hydrate",                        "hidrato" },
    { "chokeSource",                    "fonteChoke" },
    { "multiBcs",                       "multibcs" },
    { "volumetricPump",                 "bombaVolumetrica" },
    { "pressureDrop",                   "deltaPressao" },
    { "heatSource",                     "fonteCalor" },
    { "masterValve",                    "master1" },
    { "masterValve2",                   "master2" },
    { "pressureSource",                 "fontePressao" },
    { "productionTrend",                "tendP" },
    { "transientProductionTrend",       "tendTransP" },
    { "serviceTrend",                   "tendS" },
    { "transientServiceTrend",          "tendTransS" },
    { "screenConfig",                   "tela" },
    { "gasInjection",                   "gasInj" },
    { "productionProfile",              "perfilProducao" },
    { "serviceProfile",                 "perfilServico" },
    { "transientProductionProfiles",    "perfisTransP" },
    { "transientServiceProfiles",       "perfisTransS" },
    { "separator",                      "separador" },
    { "correction",                     "correcao" },
    { "surfaceChoke",                   "chokeSup" },
    { "injectionChoke",                 "chokeInj" },
    { "injectionWellBC",                "CondicaoContPocInjec" },
    { "severeSlugging",                 "intermitenciaSevera" },

    // -----------------------------------------------------------------------
    // Shared keys (appear in multiple objects at different nesting levels)
    // -----------------------------------------------------------------------
    { "active",                         "ativo" },
    // "id" is identical – omitted
    { "time",                           "tempo" },
    { "temperature",                    "temperatura" },
    { "pressure",                       "pressao" },
    { "measuredLength",                 "comprimentoMedido" },
    { "opening",                        "abertura" },
    { "prodFluidId",                    "indiFluidoPro" },
    { "outerDiameter",                  "diametroExterno" },
    { "roughness",                      "rugosidade" },
    { "conductivity",                   "condutividade" },
    { "specificHeat",                   "calorEspecifico" },
    { "density",                        "massaEspecifica" },
    { "annular",                        "anular" },
    { "discretization",                 "discretizacao" },
    { "flowRate",                       "vazao" },
    { "frequency",                      "frequencia" },
    { "power",                          "potencia" },
    { "efficiency",                     "eficiencia" },
    { "stages",                         "nestag" },
    { "manufacturerStages",             "nestagFab" },
    { "motorEfficiency",                "EficienciaMotor" },
    { "minFrequency",                   "FrequenciaMinima" },
    { "referenceFreq",                  "freqref" },
    { "hiCorrection",                   "correcHI" },
    { "motorHeatFraction",              "fracTermMotorEfic" },
    { "curve",                          "curva" },

    // -----------------------------------------------------------------------
    // configuracaoInicial
    // -----------------------------------------------------------------------
    { "geometryOrigin",                 "origemGeometria" },
    { "classicOutput",                  "saidaClassica" },
    { "geometryFollowsFlow",            "sentidoGeometriaSegueEscoamento" },
    { "gasLine",                        "linhaGas" },
    { "screenPrint",                    "saidaTela" },
    { "thermalEquilibrium",             "equilibrioTermico" },
    { "latentHeat",                     "latente" },
    { "latentHeatCond",                 "condlatente" },
    { "pvtFile",                        "pvtsimArq" },
    { "flashTableFluidModel",           "modeloFluidoTabelaFlash" },
    { "compositionalFluidModel",        "modeloFluidoComposicional" },
    { "dynamicTableModel",              "modeloTabelaDinamica" },
    { "cpModel",                        "modeloCp" },
    { "jtlModel",                       "modeloJTL" },
    { "pressureTable",                  "tabP" },
    { "sensitivityAnalysis",            "AS" },
    { "parallelizeSA",                  "paralelizaAS" },
    { "trackGOR",                       "trackRgo" },
    { "trackGasDensity",                "trackDensidadeGas" },
    { "freeGasDensityCorrection",       "correcaoDenGasLivreBlackOil" },
    { "rspbTable",                      "tabelaRSPB" },
    { "fluidProperties",                "propFluido" },
    { "initialFluidId",                 "iniFluidoP" },
    { "gasTable",                       "tabG" },
    { "steadyStateSlip",                "escorregamentoPermanente" },
    { "transientSlip",                  "escorregamentoTransiente" },
    { "flowRegimeMap",                  "mapaArranjo" },
    { "initialCondition",               "condicaoInicial" },
    { "steadyStateOrder",               "ordemperm" },
    { "snapshotFile",                   "SnapShotArq" },
    { "hisep",                          "HISEP" },
    { "fluidSalinity",                  "SalinidadeFluido" },
    { "gasLineInterfaceLength",         "comprimentoMedidoInterfaceLinhaGas" },
    { "prodLineInterfaceLength",        "comprimentoMedidoInterfaceLinhaProd" },
    { "dischargeControl",               "controleDescarga" },
    { "dischargeParameters",            "parametrosDescarga" },
    { "transient",                      "transiente" },
    { "massTransfer",                   "transferenciaMassa" },
    { "checkValve",                     "CheckValve" },
    { "advanced",                       "Avancado" },
    { "pressureCondition",              "condicaoPressao" },
    { "flowPressureCondition",          "condicaoVazPres" },
    { "slipCorrelations",               "correlacoesEscorregamento" },
    { "correlationsByRegime",           "correlacoesPorArranjo" },
    { "formation",                      "Formacao" },
    { "fluidType",                      "tipoFluido" },
    { "reverseTemp",                    "tempReves" },
    { "reverseGasRatio",                "razCompGasReves" },
    { "steadyGuess",                    "chutePerm" },
    { "xyMode",                         "modoXY" },
    { "xProdStart",                     "xProdInicio" },
    { "yProdStart",                     "yProdInicio" },
    { "xServiceStart",                  "xServInicio" },
    { "yServiceStart",                  "yServInicio" },
    { "waxMode",                        "modoParafina" },
    { "driftModel",                     "tipoModeloDrift" },
    { "diffusion3dMode",                "modoDifus3D" },
    { "diffusion3dThreads",             "threadP3D" },
    { "diffusion3dJson",                "modoDifus3DJson" },

    // dischargeParameters sub-keys
    { "maxDischargeFlow",               "vazaoLimiteDescarga" },
    { "maxDischargePressure",           "pressaoLimiteDescarga" },
    { "minDischargePressure",           "pressaoMinimaDescarga" },
    { "workGasChargePressure",          "pressaoTrabalhoDescargaGas" },
    { "maxGasChargePressure",           "pressaoLimiteDescargaGas" },
    { "minGasChargePressure",           "pressaoMinimaDescargaGas" },
    { "initialGasChargePressure",       "pressaoInicialDescargaGas" },
    { "dischargeTemperature",           "temperaturaDescarga" },
    { "latencyTime",                    "tempoLatencia" },

    // advanced (Avancado) sub-keys
    { "monophasicCriterion",            "CriterioMonofasico" },
    { "condensationCriterion",          "CriterioCondensacao" },
    { "minTimestepCriterion",           "CriterioDTMin" },
    { "falseCordSearchCriterion",       "CriterioBuscaFalsaCorda" },
    { "neglectRate",                    "taxaDespre" },
    { "simplePressureFrontier",         "MedSimpPresFront" },
    { "liquidJTSimple",                 "JTLiquidoSimple" },
    { "massTransferLimit",              "limTransMass" },
    { "relaxChokeTimestep",             "RelaxaDTChoke" },
    { "disablePenalizeTimestep",        "desligaPenalizaDT" },
    { "valveTimestepControl",           "controleDTvalv" },
    { "steadyConvergenceCriterion",     "CriterioConvergPerm" },
    { "accelerateSteadyConvergence",    "AceleraConvergPerm" },
    { "slipBoundaryCell",               "escorregamentoCelulaContorno" },
    { "counterflowCorrectionSteady",    "correcaoContracorPerm" },
    { "columnStabilization",            "estabCol" },
    { "compModelCorrectionTime",        "TcorrecaoModComp" },
    { "compModelCorrectionFlag",        "correcaoModComp" },
    { "disableMassTransferTempDeriv",   "desligaDeriTransMassDTemp" },
    { "correctSepCondition",            "corrigeContSep" },
    { "strongAnnularColCoupling",       "acopColAnulPermForte" },
    { "areaChange",                     "mudaArea" },
    { "threads",                        "nthrd" },
    { "matrixThreads",                  "nthrdMatriz" },
    { "dynTableMinDelay",               "miniTabDinAtraso" },
    { "dynTableMinDp",                  "miniTabDinDp" },
    { "dynTableMinDt",                  "miniTabDinDt" },
    { "sonicTimes",                     "Tsonico" },
    { "sonicFlags",                     "sonico" },

    // pressureCondition sub-keys
    { "fluidId",                        "indFluido" },
    { "voidFraction",                   "titulo" },
    { "waterCut",                       "razaoBeta" },

    // flowPressureCondition sub-keys
    { "massFlowRate",                   "VazMass" },

    // correlationsByRegime sub-keys
    { "stratified",                     "estratificado" },
    { "slugBubble",                     "bolhaGolfada" },
    { "annularChurn",                   "anularChurn" },

    // formation sub-keys
    { "properties",                     "Propriedades" },
    { "productionTime",                 "TempoProducao" },

    // -----------------------------------------------------------------------
    // timeSettings (tempo)
    // -----------------------------------------------------------------------
    { "finalTime",                      "tempoFinal" },
    { "times",                          "tempos" },
    { "maxTimestep",                    "dtmax" },
    { "segregationTime",                "tempoSegrega" },
    { "segregation",                    "segrega" },
    { "saveSnapshot",                   "gravaMomento" },

    // -----------------------------------------------------------------------
    // pvtTable (tabela)
    // -----------------------------------------------------------------------
    { "numPoints",                      "nPontos" },
    { "maxPressure",                    "pressaoMaxima" },
    { "minPressure",                    "pressaoMinima" },
    { "maxTemperature",                 "temperaturaMaxima" },
    { "minTemperature",                 "temperaturaMinima" },

    // -----------------------------------------------------------------------
    // productionFluids Item
    // -----------------------------------------------------------------------
    { "gor",                            "rgo" },
    { "gasDensity",                     "densidadeGas" },
    { "waterDensity",                   "densidadeAgua" },
    { "emulsionType",                   "tipoEmul" },
    { "emulsionCoefA",                  "coefAModeloExp" },
    { "emulsionCoefB",                  "coefBModeloExp" },
    { "phi100",                         "PHI100" },
    { "bswCut",                         "bswCorte" },
    { "bswVec",                         "BSWVec" },
    { "co2Fraction",                    "fracCO2" },
    { "criticalCorrelation",            "correlacaoCritica" },
    { "rsPbModel",                      "modeloRsPb" },
    { "deadOilModel",                   "modeloOleoMorto" },
    { "deadOilTemps",                   "tempOleoMorto" },
    { "deadOilViscs",                   "viscOleoMorto" },
    { "liveOilModel",                   "modeloOleoVivo" },
    { "undersaturatedOilModel",         "modeloOleoSubSaturado" },
    { "blackOilViscModel",              "modeloViscBlackOil" },
    { "blackOilWaterModel",             "modeloAguaBlackOil" },
    { "userMolarFraction",              "fracMolarUsuario" },
    { "molarFraction",                  "fracMolar" },
    { "userGORComp",                    "RGOCompUsuario" },

    // -----------------------------------------------------------------------
    // gasFluid (fluidoGas)
    // -----------------------------------------------------------------------
    { "useFlashTable",                  "usaTabelaFlash" },

    // -----------------------------------------------------------------------
    // complementaryFluid (fluidoComplementar)
    // -----------------------------------------------------------------------
    { "pressureCompressibility",        "compP" },
    { "thermalCompressibility",         "compT" },
    { "surfaceTension",                 "tensup" },
    { "salinity",                       "salinidade" },
    { "complementaryFluidType",         "tipoF" },

    // -----------------------------------------------------------------------
    // material Item
    // -----------------------------------------------------------------------
    { "type",                           "tipo" },

    // -----------------------------------------------------------------------
    // crossSection (secaoTransversal) Item and layers
    // -----------------------------------------------------------------------
    { "innerDiameter",                  "diametroInterno" },
    { "layers",                         "camadas" },
    { "layerMeasurementType",           "tipoMedicaoCamada" },
    { "diameter",                       "diametro" },
    { "thickness",                      "espessura" },
    { "materialId",                     "idMaterial" },

    // -----------------------------------------------------------------------
    // productionDucts / serviceDucts Item
    // -----------------------------------------------------------------------
    { "mr2Correlation",                 "correlacaoMR2" },
    { "angle",                          "angulo" },
    { "crossSectionId",                 "idCorte" },
    { "formationId",                    "idFormacao" },
    { "externalEnvironment",            "ambienteExterno" },
    { "convectionDirection",            "direcaoConveccao" },
    { "thermalCoupling",                "acoplamentoTermico" },
    { "parallelNetworkThermalCoupling", "acoplamentoTermicoRedeParalela" },
    { "grouping",                       "agrupamento" },
    { "cellDx",                         "dxCelula" },
    { "initialConditions",              "condicoesIniciais" },
    { "initialConditionsAndAmbient",    "condicoesIniciaisEAmbiente" },
    { "hydroGradient",                  "dPdLHidro" },
    { "fricGradient",                   "dPdLFric" },
    { "tempGradient",                   "dTdL" },
    { "diffusion2d",                    "difusTerm2D" },
    { "diffusion2dJson",                "difusTerm2DJSON" },
    { "diffusion3d",                    "difusTerm3D" },
    { "diffusion3dFE",                  "difusTerm3DFE" },
    { "diffusion3dCoupling",            "difusTerm3DAcop" },
    { "xCoord",                         "xCoor" },
    { "yCoord",                         "yCoor" },
    { "cellsXY",                        "nCelulas_XY" },

    // discretizacao items
    { "numCells",                       "nCelulas" },
    { "length",                         "comprimento" },

    // initialConditions sub-keys
    { "measuredPositions",              "compInter" },
    // "temp" is kept as-is (it's already English-style)
    { "holdup",                         "holdup" },  // same
    { "waterCutIC",                     "bet" },
    { "superficialLiquidVel",           "uls" },
    { "superficialGasVel",              "ugs" },
    { "externalTemp",                   "tempExterna" },
    { "externalVel",                    "velExterna" },
    { "externalConductivity",           "kExterna" },
    { "externalSpecificHeat",           "calorEspecificoExterno" },
    { "externalDensity",                "rhoExterno" },
    { "externalVisc",                   "viscExterna" },
    { "gasMassFlowRate",                "vazaoMassicaGas" },

    // -----------------------------------------------------------------------
    // valve (valvula) Item
    // -----------------------------------------------------------------------
    { "cvCurve",                        "curvaCV" },
    { "dynamicCurve",                   "curvaDinamic" },

    // -----------------------------------------------------------------------
    // gasSource Item
    // -----------------------------------------------------------------------
    { "gasFlowRate",                    "vazaoGas" },
    { "complementaryFluidFlowRate",     "vazaoFluidoComplementar" },
    { "dry",                            "seco" },

    // -----------------------------------------------------------------------
    // liquidSource Item
    // -----------------------------------------------------------------------
    { "liquidFlowRate",                 "vazaoLiquido" },

    // -----------------------------------------------------------------------
    // massSource Item
    // -----------------------------------------------------------------------
    { "termType",                       "tipoTermo" },
    { "totalMassFlowRate",              "vazaoMassT" },
    { "complementaryMassFlowRate",      "vazaoMassC" },
    // gasMassFlowRate → vazaoMassG (already mapped above to vazaoMassicaGas
    // for service ducts; massSource uses vazaoMassG). Use separate key:
    { "massMassFlowGas",                "vazaoMassG" },

    // -----------------------------------------------------------------------
    // gasLiftSource (fonteGasLift) Item
    // -----------------------------------------------------------------------
    { "annularColumn",                  "colunaEanular" },
    { "prodMeasuredLength",             "comprimentoMedidoProducao" },
    { "serviceMeasuredLength",          "comprimentoMedidoServico" },
    { "valveType",                      "tipoValvula" },
    { "orificeDiameter",                "diametroOrificio" },
    { "vglDischCoef",                   "cdvgl" },
    { "recoveryFreq",                   "frecupera" },
    { "liquidDischCoef",                "cdvLiq" },
    { "liquidRecoveryFreq",             "frecuperaLiq" },
    { "areaRatio",                      "razaoArea" },
    { "calibrationPressure",            "pressaoCalibracao" },
    { "calibrationTemperature",         "temperaturaCalibracao" },

    // -----------------------------------------------------------------------
    // porousRadialSource / porous2DSource Item
    // -----------------------------------------------------------------------
    { "file",                           "arquivo" },

    // -----------------------------------------------------------------------
    // IPR Item
    // -----------------------------------------------------------------------
    { "iprType",                        "tipoIPR" },
    { "staticPressure",                 "pressaoEstatica" },
    { "staticPressureTime",             "tempoPressaoEstatica" },
    { "qMaxTime",                       "tempoqMax" },
    { "temperatures",                   "temperaturas" },
    { "temperaturesTime",               "tempoTemperaturas" },
    { "ipTime",                         "tempoip" },
    { "iiTime",                         "tempoii" },

    // -----------------------------------------------------------------------
    // chokeSource (fonteChoke)
    // -----------------------------------------------------------------------
    { "dischargeCoefficient",           "coeficienteDescarga" },
    { "model",                          "modelo" },

    // -----------------------------------------------------------------------
    // pressureDrop (deltaPressao) Item
    // -----------------------------------------------------------------------
    // Note: the section name "pressureDrop" → "deltaPressao" is already
    // mapped above. The inner field also named "deltaPressao" is written
    // "pressureDropValue" in English to avoid ambiguity.
    { "pressureDropValue",              "deltaPressao" },
    { "gasCompType",                    "tipoCompGas" },
    { "polyCorrFactor",                 "fatPoli" },
    { "liquidEfficiency",               "eficLiq" },
    { "gasEfficiency",                  "eficGas" },

    // -----------------------------------------------------------------------
    // heatSource (fonteCalor) Item
    // -----------------------------------------------------------------------
    { "heat",                           "calor" },

    // -----------------------------------------------------------------------
    // masterValve (master1)
    // -----------------------------------------------------------------------
    { "activeAreaRatio",                "razaoAreaAtiva" },

    // -----------------------------------------------------------------------
    // pig Item
    // -----------------------------------------------------------------------
    { "launcher",                       "lancador" },
    { "receiver",                       "recebedor" },
    { "areaGap",                        "folgaArea" },
    { "pigDischCoef",                   "cdPig" },

    // -----------------------------------------------------------------------
    // hydrate (hidrato)
    // -----------------------------------------------------------------------
    { "internalCalc",                   "calculoInterno" },
    { "hydrateModel",                   "modeloHidrato" },
    { "hydrateFluidProps",              "PropFluHidrato" },
    { "TurnerModel",                    "ModeloTurner" },
    { "inhibitor",                      "inibidor" },
    { "loadedWaterFraction",            "fracFWcarregada" },
    { "stoichCoef",                     "coefEsteq" },
    { "hydrateStructure",               "estruturaHidratos" },
    { "hydrateModelType",               "tipoHmodel" },

    // -----------------------------------------------------------------------
    // wax (parafina)
    // -----------------------------------------------------------------------
    { "waxFile",                        "arquivoWax" },
    { "userPorosity",                   "usuarioPorosidade" },
    { "porosity",                       "porosidade" },
    { "userC2C3",                       "usuarioC2C3" },
    { "userDiffusion",                  "usuarioDifus" },
    { "changeFluidVisc",                "alteraViscFlu" },
    { "diffusivity",                    "difus" },
    { "viscMultiplier",                 "multVis" },
    { "waxMultD",                       "DmultipWax" },
    { "waxMultE",                       "EmultipWax" },
    { "waxMultF",                       "FmultipWax" },

    // -----------------------------------------------------------------------
    // volumetricPump (bombaVolumetrica)
    // -----------------------------------------------------------------------
    { "capacity",                       "capacidade" },
    { "polyFactor",                     "fatorpoli" },

    // -----------------------------------------------------------------------
    // BCS / multiBcs shared
    // -----------------------------------------------------------------------
    { "thermalEquilBcs",                "equilTerm" },

    // -----------------------------------------------------------------------
    // injectionWellBC (CondicaoContPocInjec)
    // -----------------------------------------------------------------------
    { "pvtsimFile",                     "arquivoPvtsim" },
    { "boundaryCondition",              "condContorno" },
    { "bcType",                         "tipoCC" },
    { "injectionPressure",              "presInjec" },
    { "injectionPressures",             "pressaoInjecao" },
    { "initialFlowGuess",               "chuteVazaoInjecao" },
    { "openingType",                    "TipoAbertura" },

    // -----------------------------------------------------------------------
    // complementaryFluid extra
    // -----------------------------------------------------------------------
    { "gasAmbient",                     "ambienteGas" },

    // -----------------------------------------------------------------------
    // Screen output config (tela) — item identifiers
    // -----------------------------------------------------------------------
    { "tubeColumn",                     "coluna" },
    { "cellIndex",                      "celula" },
    { "layerIndex",                     "camada" },
    { "outputVariable",                 "variavel" },
    { "label",                          "rotulo" },
    { "ambientHoldup",                  "titAmb" },
    { "chokeTime",                      "tempoChk" },

    // -----------------------------------------------------------------------
    // Screen output config — boolean output flags
    // (users set these to true/false in the tela section)
    // -----------------------------------------------------------------------
    { "flowPattern",                    "arra" },
    { "frictionPressure",               "fric" },
    { "hydrostaticPressure",            "hidro" },
    { "cpGas",                          "cpgas" },
    { "cpLiquid",                       "cpliq" },
    { "massTransferRate",               "masstrans" },
    { "sourceGasMassFlow",              "mgFonte" },
    { "sourceLiqMassFlow",              "mlFonte" },
    { "sourceCompMassFlow",             "mcFonte" },
    { "waveVelocity",                   "autoVal" },
    { "eigenVelocity",                  "autoVel" },
    { "fluctuation",                    "flutuacao" },
    { "ambientTemperature",             "temperaturaAmbiente" },
    { "yCO2",                           "yco2" },
    { "frictionReducer",                "RedutorAtrito" },
    { "formationVolumeFactor",          "Bo" },
    { "froudeNumber",                   "Froud" },
    { "externalGrashof",                "GrashExterno" },
    { "internalGrashof",                "GrashInterno" },
    { "externalHeatCoef",               "Hext" },
    { "internalHeatCoef",               "Hint" },
    { "externalNusselt",                "NusselExterno" },
    { "internalNusselt",                "NusselInterno" },
    { "externalPrandtl",                "PrandtlExterno" },
    { "gasPrandtl",                     "PrandtlGas" },
    { "internalPrandtl",                "PrandtlInterno" },
    { "liquidPrandtl",                  "PrandtlLiquido" },
    { "stdGasFlowRate",                 "QGstd" },
    { "stdLiqWaterFlowRate",            "QLWstd" },
    { "stdLiqFlowRate",                 "QLstd" },
    { "stdTotalLiqFlowRate",            "QLstdTotal" },
    { "gasOilRatioOut",                 "RGO" },
    { "gasOilRatioStd",                 "RS" },
    { "gasSolubility",                  "Rs" },
    { "externalReynolds",               "ReyExterno" },
    { "internalReynolds",               "ReyInterno" },
    { "residenceTime",                  "TResi" },
    { "term1",                          "Term1" },
    { "term2",                          "Term2" },
    { "maxThroatVelocity",              "VelocidadeMaximaGarganta" },
    { "gasStdDensity",                  "deng" },
    { "downstreamGasDensity",           "dengD" },
    { "upstreamGasDensity",             "dengL" },
    { "pumpDeltaP",                     "deltaPBomba" },
    { "pumpHead",                       "head" },
    { "pumpPower",                      "potenciaBomba" },
    { "vglStagePressure",               "presEstagVGL" },
    { "vglThroatPressure",              "presGargVGL" },
    { "bottomholePressure",             "presFundo" },
    { "pseudoLiquid",                   "pseudoLiquido" },
    { "mixtureDensity",                 "rhoMix" },
    { "gasInSituDensity",               "rhog" },
    { "liquidInSituDensity",            "rhol" },
    { "chokeDownstreamTemp",            "tempChokeJusante" },
    { "vglStageTemp",                   "tempEstagVGL" },
    { "vglThroatTemp",                  "tempGargVGL" },
    { "injectionTemp",                  "tempInj" },
    { "shearStress",                    "tensaoCisalhamento" },
    { "driftVelocity",                  "ud" },
    { "gasVelocity",                    "ug" },
    { "liquidVelocity",                 "ul" },
    { "liquidFlowRateSc",               "vazLiq" },
    { "liquidMassFlowRate",             "vazaoMassicaLiquido" },
    { "vglFlowRate",                    "vazaoVGL" },
    { "gasViscosity",                   "viscosidadeGas" },
    { "liquidViscosity",                "viscosidadeLiquido" },
    { "m1ProdDownstreamVol",            "volJusM1PT" },
    { "m1ServDownstreamVol",            "volJusM1ST" },
    { "m1ProdUpstreamVol",              "volMonM1PT" },
    { "m1ServUpstreamVol",              "volMonM1ST" },

    // -----------------------------------------------------------------------
    // Production fluid — extra emulsion
    // -----------------------------------------------------------------------
    { "emulsionVec",                    "emulVec" },

    // -----------------------------------------------------------------------
    // IPR / sources — fluid index (different from indiFluidoPro)
    // -----------------------------------------------------------------------
    { "prodFluidIndex",                 "indFluidoPro" },

    // -----------------------------------------------------------------------
    // Severe slugging (intermitenciaSevera)
    // -----------------------------------------------------------------------
    { "accumStartLength",               "inicioTrechoAcumula" },
    { "accumEndLength",                 "fimTrechoAcumula" },
    { "columnEndLength",                "fimTrechoColuna" },
    { "voidFractionPenetration",        "fracaoVazioPenetracao" },
    { "criterion",                      "criterio" },

    // -----------------------------------------------------------------------
    // Sensitivity analysis — variable descriptor
    // -----------------------------------------------------------------------
    { "variable",                       "variavel" },

}; // EN_TO_PT


// ---------------------------------------------------------------------------
// Value translations: for specific English keys, translate string values
// to their Portuguese equivalents before the key itself is renamed.
// Matching is case-insensitive (values are upper-cased before lookup).
// ---------------------------------------------------------------------------
static const std::unordered_map<
    std::string,
    std::unordered_map<std::string, std::string>> VALUE_TRANSLATIONS = {
    { "layerMeasurementType", {
        { "THICKNESS", "ESPESSURA" },
        { "DIAMETER",  "DIAMETRO"  },
    }},
};

// ---------------------------------------------------------------------------
// Internal implementation
// ---------------------------------------------------------------------------

static void renameObjectKeys(Value& obj, Document::AllocatorType& alloc) {
    // Collect keys that need renaming (avoid modifying map while iterating)
    std::vector<std::pair<std::string, std::string>> renames;
    renames.reserve(obj.MemberCount());

    for (Value::MemberIterator it = obj.MemberBegin();
         it != obj.MemberEnd(); ++it)
    {
        const std::string key(it->name.GetString());

        // Translate string values for specific keys (e.g. layerMeasurementType)
        auto vit = VALUE_TRANSLATIONS.find(key);
        if (vit != VALUE_TRANSLATIONS.end() && it->value.IsString()) {
            std::string val(it->value.GetString());
            std::transform(val.begin(), val.end(), val.begin(), ::toupper);
            auto vfound = vit->second.find(val);
            if (vfound != vit->second.end()) {
                it->value.SetString(vfound->second.c_str(),
                                    static_cast<SizeType>(vfound->second.size()),
                                    alloc);
            }
        }

        auto found = EN_TO_PT.find(key);
        if (found != EN_TO_PT.end() && found->second != key) {
            renames.emplace_back(key, found->second);
        }
    }

    // Apply renames: extract value, remove old member, add with new key name
    for (const auto& rename : renames) {
        const std::string& oldKey = rename.first;
        const std::string& newKey = rename.second;

        Value::MemberIterator it = obj.FindMember(oldKey.c_str());
        if (it == obj.MemberEnd()) {
            continue; // already removed (shouldn't happen, but be safe)
        }

        Value val(kNullType);
        val.Swap(it->value);          // take the value before removal
        obj.RemoveMember(it);         // remove old key

        Value nameStr(newKey.c_str(),
                      static_cast<SizeType>(newKey.size()), alloc);
        obj.AddMember(nameStr, val, alloc);
    }
}

static void translateValueRecursive(Value& v,
                                    Document::AllocatorType& alloc) {
    if (v.IsObject()) {
        // Recurse into children first (depth-first)
        for (Value::MemberIterator it = v.MemberBegin();
             it != v.MemberEnd(); ++it)
        {
            translateValueRecursive(it->value, alloc);
        }
        // Then rename this object's own keys
        renameObjectKeys(v, alloc);
    } else if (v.IsArray()) {
        for (Value::ValueIterator it = v.Begin(); it != v.End(); ++it) {
            translateValueRecursive(*it, alloc);
        }
    }
    // Scalars (string, number, bool, null) have no keys – nothing to do
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

namespace JSONKeyTranslator {

void translateEnToPt(Value& node, Document::AllocatorType& allocator) {
    translateValueRecursive(node, allocator);
}

} // namespace JSONKeyTranslator
