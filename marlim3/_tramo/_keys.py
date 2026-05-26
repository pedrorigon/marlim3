"""_keys.py — Python-side Portuguese → English JSON key translator.

Used by Branch.from_json() to convert old Portuguese-keyed JSON files to
the English key scheme before populating Branch attributes.

The C++ engine handles the reverse direction (EN→PT) at simulation time
when ``"language": "en"`` is present in the JSON.
"""

# ---------------------------------------------------------------------------
# PT_TO_EN: flat map from Portuguese JSON key → English JSON key.
# Language-neutral keys (ipr, pig, api, bsw, rho, ip, ii, etc.) are
# intentionally omitted — they pass through unchanged.
# ---------------------------------------------------------------------------
PT_TO_EN = {
    # Root-level / top-level object keys
    # Legacy alias kept for backward compatibility with existing inputs.
    "sistema":              "system",
    "versaoJson":           "jsonVersion",
    "versaoJSON":           "jsonVersion",
    "configuracaoInicial":  "initialConfig",
    "tabela":               "pvtTable",
    "correcao":             "correction",
    "fluidosProducao":      "productionFluids",
    "fluidoComplementar":   "complementaryFluid",
    "fluidoGas":            "gasFluid",
    "secaoTransversal":     "crossSection",
    "dutosProducao":        "productionDucts",
    "dutosServico":         "serviceDucts",
    "hidrato":              "hydrate",
    "fonteChoke":           "chokeSource",
    "multibcs":             "multiBcs",
    "bombaVolumetrica":     "volumetricPump",
    "deltaPressao":         "pressureDrop",
    "fonteCalor":           "heatSource",
    "master1":              "masterValve",
    "master2":              "masterValve2",
    "fontePressao":         "pressureSource",
    "tendP":                "productionTrend",
    "tendTransP":           "transientProductionTrend",
    "tendS":                "serviceTrend",
    "tendTransS":           "transientServiceTrend",
    "tela":                 "screenConfig",
    "gasInj":               "gasInjection",
    "perfilProducao":       "productionProfile",
    "perfilServico":        "serviceProfile",
    "perfisTransP":         "transientProductionProfiles",
    "perfisTransS":         "transientServiceProfiles",
    "separador":            "separator",
    "chokeSup":             "surfaceChoke",
    "chokeInj":             "injectionChoke",
    "CondicaoContPocInjec": "injectionWellBC",
    "intermitenciaSevera":  "severeSlugging",
    "parafina":             "wax",
    "valvula":              "valve",
    "fonteLiquido":         "liquidSource",
    "fonteMassa":           "massSource",
    "fonteGas":             "gasSource",
    "fontePoroRadial":      "porousRadialSource",
    "fontePoro2D":          "porous2DSource",
    "fonteGasLift":         "gasLiftSource",
    "bcs":                  "esp",

    # Shared keys (appear at multiple nesting levels)
    "ativo":                "active",
    # "tempo" as inner dict key maps to "time"; at root level it is handled
    # separately via ROOT_PT_TO_EN (mapped to "timeSettings").
    "tempo":                "time",
    "temperatura":          "temperature",
    "pressao":              "pressure",
    "holdup":               "holdup",
    "comprimentoMedido":    "measuredLength",
    "abertura":             "opening",
    "indiFluidoPro":        "prodFluidId",
    "diametroExterno":      "outerDiameter",
    "rugosidade":           "roughness",
    "condutividade":        "conductivity",
    "calorEspecifico":      "specificHeat",
    "massaEspecifica":      "density",
    "anular":               "annular",
    "discretizacao":        "discretization",
    "vazao":                "flowRate",
    "frequencia":           "frequency",
    "potencia":             "power",
    "eficiencia":           "efficiency",
    "nestag":               "stages",
    "nestagFab":            "manufacturerStages",
    "EficienciaMotor":      "motorEfficiency",
    "FrequenciaMinima":     "minFrequency",
    "freqref":              "referenceFreq",
    "correcHI":             "hiCorrection",
    "fracTermMotorEfic":    "motorHeatFraction",
    "curva":                "curve",

    # configuracaoInicial
    "origemGeometria":      "geometryOrigin",
    "saidaClassica":        "classicOutput",
    "sentidoGeometriaSegueEscoamento": "geometryFollowsFlow",
    "linhaGas":             "gasLine",
    "saidaTela":            "screenPrint",
    "equilibrioTermico":    "thermalEquilibrium",
    "latente":              "latentHeat",
    "condlatente":          "latentHeatCond",
    "pvtsimArq":            "pvtFile",
    "modeloFluidoTabelaFlash": "flashTableFluidModel",
    "modeloFluidoComposicional": "compositionalFluidModel",
    "modeloTabelaDinamica": "dynamicTableModel",
    "modeloCp":             "cpModel",
    "modeloJTL":            "jtlModel",
    "tabP":                 "pressureTable",
    "AS":                   "sensitivityAnalysis",
    "paralelizaAS":         "parallelizeSA",
    "trackRgo":             "trackGOR",
    "trackDensidadeGas":    "trackGasDensity",
    "correcaoDenGasLivreBlackOil": "freeGasDensityCorrection",
    "tabelaRSPB":           "rspbTable",
    "propFluido":           "fluidProperties",
    "iniFluidoP":           "initialFluidId",
    "tabG":                 "gasTable",
    "escorregamentoPermanente": "steadyStateSlip",
    "escorregamentoTransiente": "transientSlip",
    "mapaArranjo":          "flowRegimeMap",
    "condicaoInicial":      "initialCondition",
    "ordemperm":            "steadyStateOrder",
    "SnapShotArq":          "snapshotFile",
    "HISEP":                "hisep",
    "SalinidadeFluido":     "fluidSalinity",
    "comprimentoMedidoInterfaceLinhaGas": "gasLineInterfaceLength",
    "comprimentoMedidoInterfaceLinhaProd": "prodLineInterfaceLength",
    "controleDescarga":     "dischargeControl",
    "parametrosDescarga":   "dischargeParameters",
    "transiente":           "transient",
    "transferenciaMassa":   "massTransfer",
    "CheckValve":           "checkValve",
    "Avancado":             "advanced",
    "condicaoPressao":      "pressureCondition",
    "condicaoVazPres":      "flowPressureCondition",
    "correlacoesEscorregamento": "slipCorrelations",
    "correlacoesPorArranjo": "correlationsByRegime",
    "Formacao":             "formation",
    "tipoFluido":           "fluidType",
    "tempReves":            "reverseTemp",
    "razCompGasReves":      "reverseGasRatio",
    "chutePerm":            "steadyGuess",
    "modoXY":               "xyMode",
    "xProdInicio":          "xProdStart",
    "yProdInicio":          "yProdStart",
    "xServInicio":          "xServiceStart",
    "yServInicio":          "yServiceStart",
    "modoParafina":         "waxMode",
    "tipoModeloDrift":      "driftModel",
    "modoDifus3D":          "diffusion3dMode",
    "threadP3D":            "diffusion3dThreads",
    "modoDifus3DJson":      "diffusion3dJson",

    # dischargeParameters sub-keys
    "vazaoLimiteDescarga":  "maxDischargeFlow",
    "pressaoLimiteDescarga": "maxDischargePressure",
    "pressaoMinimaDescarga": "minDischargePressure",
    "pressaoTrabalhoDescargaGas": "workGasChargePressure",
    "pressaoLimiteDescargaGas": "maxGasChargePressure",
    "pressaoMinimaDescargaGas": "minGasChargePressure",
    "pressaoInicialDescargaGas": "initialGasChargePressure",
    "temperaturaDescarga":  "dischargeTemperature",
    "tempoLatencia":        "latencyTime",

    # advanced sub-keys
    "CriterioMonofasico":   "monophasicCriterion",
    "CriterioCondensacao":  "condensationCriterion",
    "CriterioDTMin":        "minTimestepCriterion",
    "CriterioBuscaFalsaCorda": "falseCordSearchCriterion",
    "taxaDespre":           "neglectRate",
    "MedSimpPresFront":     "simplePressureFrontier",
    "JTLiquidoSimple":      "liquidJTSimple",
    "limTransMass":         "massTransferLimit",
    "RelaxaDTChoke":        "relaxChokeTimestep",
    "desligaPenalizaDT":    "disablePenalizeTimestep",
    "controleDTvalv":       "valveTimestepControl",
    "CriterioConvergPerm":  "steadyConvergenceCriterion",
    "AceleraConvergPerm":   "accelerateSteadyConvergence",
    "escorregamentoCelulaContorno": "slipBoundaryCell",
    "correcaoContracorPerm": "counterflowCorrectionSteady",
    "estabCol":             "columnStabilization",
    "TcorrecaoModComp":     "compModelCorrectionTime",
    "correcaoModComp":      "compModelCorrectionFlag",
    "desligaDeriTransMassDTemp": "disableMassTransferTempDeriv",
    "corrigeContSep":       "correctSepCondition",
    "acopColAnulPermForte": "strongAnnularColCoupling",
    "mudaArea":             "areaChange",
    "nthrd":                "threads",
    "nthrdMatriz":          "matrixThreads",
    "miniTabDinAtraso":     "dynTableMinDelay",
    "miniTabDinDp":         "dynTableMinDp",
    "miniTabDinDt":         "dynTableMinDt",
    "Tsonico":              "sonicTimes",
    "sonico":               "sonicFlags",

    # pressureCondition / flowPressureCondition sub-keys
    "indFluido":            "fluidId",
    "titulo":               "voidFraction",
    "razaoBeta":            "waterCut",
    "VazMass":              "massFlowRate",

    # correlationsByRegime sub-keys
    "estratificado":        "stratified",
    "bolhaGolfada":         "slugBubble",
    "anularChurn":          "annularChurn",

    # formation sub-keys
    "Propriedades":         "properties",
    "TempoProducao":        "productionTime",

    # pvtTable (tabela) sub-keys
    "nPontos":              "numPoints",
    "pressaoMaxima":        "maxPressure",
    "pressaoMinima":        "minPressure",
    "temperaturaMaxima":    "maxTemperature",
    "temperaturaMinima":    "minTemperature",

    # timeSettings (tempo) sub-keys — NOTE: root-level "tempo" maps to
    # "timeSettings" via ROOT_PT_TO_EN; these are the inner sub-keys.
    "tempoFinal":           "finalTime",
    "tempos":               "times",
    "dtmax":                "maxTimestep",
    "tempoSegrega":         "segregationTime",
    "segrega":              "segregation",
    "gravaMomento":         "saveSnapshot",

    # productionFluids item
    "rgo":                  "gor",
    "densidadeGas":         "gasDensity",
    "densidadeAgua":        "waterDensity",
    "tipoEmul":             "emulsionType",
    "coefAModeloExp":       "emulsionCoefA",
    "coefBModeloExp":       "emulsionCoefB",
    "PHI100":               "phi100",
    "bswCorte":             "bswCut",
    "BSWVec":               "bswVec",
    "fracCO2":              "co2Fraction",
    "correlacaoCritica":    "criticalCorrelation",
    "modeloRsPb":           "rsPbModel",
    "modeloOleoMorto":      "deadOilModel",
    "tempOleoMorto":        "deadOilTemps",
    "viscOleoMorto":        "deadOilViscs",
    "modeloOleoVivo":       "liveOilModel",
    "modeloOleoSubSaturado": "undersaturatedOilModel",
    "modeloViscBlackOil":   "blackOilViscModel",
    "modeloAguaBlackOil":   "blackOilWaterModel",
    "fracMolarUsuario":     "userMolarFraction",
    "fracMolar":            "molarFraction",
    "RGOCompUsuario":       "userGORComp",
    "emulVec":              "emulsionVec",

    # gasFluid sub-keys
    "usaTabelaFlash":       "useFlashTable",

    # complementaryFluid sub-keys
    "compP":                "pressureCompressibility",
    "compT":                "thermalCompressibility",
    "tensup":               "surfaceTension",
    "salinidade":           "salinity",
    "tipoF":                "complementaryFluidType",
    "ambienteGas":          "gasAmbient",

    # material item
    "tipo":                 "type",

    # crossSection item and layers
    "diametroInterno":      "innerDiameter",
    "camadas":              "layers",
    "tipoMedicaoCamada":    "layerMeasurementType",
    "diametro":             "diameter",
    "espessura":            "thickness",
    "idMaterial":           "materialId",

    # productionDucts / serviceDucts item
    "correlacaoMR2":        "mr2Correlation",
    "angulo":               "angle",
    "idCorte":              "crossSectionId",
    "idFormacao":           "formationId",
    "ambienteExterno":      "externalEnvironment",
    "direcaoConveccao":     "convectionDirection",
    "acoplamentoTermico":   "thermalCoupling",
    "acoplamentoTermicoRedeParalela": "parallelNetworkThermalCoupling",
    "agrupamento":          "grouping",
    "dxCelula":             "cellDx",
    "condicoesIniciais":    "initialConditions",
    "condicoesIniciaisEAmbiente": "initialConditionsAndAmbient",
    "dPdLHidro":            "hydroGradient",
    "dPdLFric":             "fricGradient",
    "dTdL":                 "tempGradient",
    "difusTerm2D":          "diffusion2d",
    "difusTerm2DJSON":      "diffusion2dJson",
    "difusTerm3D":          "diffusion3d",
    "difusTerm3DFE":        "diffusion3dFE",
    "difusTerm3DAcop":      "diffusion3dCoupling",
    "xCoor":                "xCoord",
    "yCoor":                "yCoord",
    "nCelulas_XY":          "cellsXY",
    "nCelulas":             "numCells",
    "comprimento":          "length",

    # initialConditions sub-keys
    "compInter":            "measuredPositions",
    "bet":                  "waterCutIC",
    "uls":                  "superficialLiquidVel",
    "ugs":                  "superficialGasVel",
    "tempExterna":          "externalTemp",
    "velExterna":           "externalVel",
    "kExterna":             "externalConductivity",
    "calorEspecificoExterno": "externalSpecificHeat",
    "rhoExterno":           "externalDensity",
    "viscExterna":          "externalVisc",
    "vazaoMassicaGas":      "gasMassFlowRate",

    # valve item
    "curvaCV":              "cvCurve",
    "curvaDinamic":         "dynamicCurve",

    # gasSource item
    "vazaoGas":             "gasFlowRate",
    "vazaoFluidoComplementar": "complementaryFluidFlowRate",
    "seco":                 "dry",

    # liquidSource item
    "vazaoLiquido":         "liquidFlowRate",

    # massSource item
    "tipoTermo":            "termType",
    "vazaoMassT":           "totalMassFlowRate",
    "vazaoMassC":           "complementaryMassFlowRate",
    "vazaoMassG":           "massMassFlowGas",

    # gasLiftSource item
    "colunaEanular":        "annularColumn",
    "comprimentoMedidoProducao": "prodMeasuredLength",
    "comprimentoMedidoServico": "serviceMeasuredLength",
    "tipoValvula":          "valveType",
    "diametroOrificio":     "orificeDiameter",
    "cdvgl":                "vglDischCoef",
    "frecupera":            "recoveryFreq",
    "cdvLiq":               "liquidDischCoef",
    "frecuperaLiq":         "liquidRecoveryFreq",
    "razaoArea":            "areaRatio",
    "pressaoCalibracao":    "calibrationPressure",
    "temperaturaCalibracao": "calibrationTemperature",

    # porousRadialSource / porous2DSource
    "arquivo":              "file",

    # IPR item
    "tipoIPR":              "iprType",
    "pressaoEstatica":      "staticPressure",
    "tempoPressaoEstatica": "staticPressureTime",
    "tempoqMax":            "qMaxTime",
    "temperaturas":         "temperatures",
    "tempoTemperaturas":    "temperaturesTime",
    "tempoip":              "ipTime",
    "tempoii":              "iiTime",

    # chokeSource
    "coeficienteDescarga":  "dischargeCoefficient",
    "modelo":               "model",

    # pressureDrop item
    "tipoCompGas":          "gasCompType",
    "fatPoli":              "polyCorrFactor",
    "eficLiq":              "liquidEfficiency",
    "eficGas":              "gasEfficiency",

    # heatSource
    "calor":                "heat",

    # masterValve
    "razaoAreaAtiva":       "activeAreaRatio",

    # pig item
    "lancador":             "launcher",
    "recebedor":            "receiver",
    "folgaArea":            "areaGap",
    "cdPig":                "pigDischCoef",

    # hydrate
    "calculoInterno":       "internalCalc",
    "modeloHidrato":        "hydrateModel",
    "PropFluHidrato":       "hydrateFluidProps",
    "ModeloTurner":         "TurnerModel",
    "inibidor":             "inhibitor",
    "fracFWcarregada":      "loadedWaterFraction",
    "coefEsteq":            "stoichCoef",
    "estruturaHidratos":    "hydrateStructure",
    "tipoHmodel":           "hydrateModelType",

    # wax (parafina)
    "arquivoWax":           "waxFile",
    "usuarioPorosidade":    "userPorosity",
    "porosidade":           "porosity",
    "usuarioC2C3":          "userC2C3",
    "usuarioDifus":         "userDiffusion",
    "alteraViscFlu":        "changeFluidVisc",
    "difus":                "diffusivity",
    "multVis":              "viscMultiplier",
    "DmultipWax":           "waxMultD",
    "EmultipWax":           "waxMultE",
    "FmultipWax":           "waxMultF",

    # volumetricPump
    "capacidade":           "capacity",
    "fatorpoli":            "polyFactor",

    # ESP
    "equilTerm":            "thermalEquilBcs",

    # injectionWellBC (CondicaoContPocInjec)
    "arquivoPvtsim":        "pvtsimFile",
    "condContorno":         "boundaryCondition",
    "tipoCC":               "bcType",
    "presInjec":            "injectionPressure",
    "pressaoInjecao":       "injectionPressures",
    "chuteVazaoInjecao":    "initialFlowGuess",
    "TipoAbertura":         "openingType",

    # screenConfig (tela) item identifiers
    "coluna":               "tubeColumn",
    "celula":               "cellIndex",
    "camada":               "layerIndex",
    "variavel":             "outputVariable",
    "rotulo":               "label",
    "titAmb":               "ambientHoldup",
    "tempoChk":             "chokeTime",

    # Screen output boolean flags
    "arra":                 "flowPattern",
    "fric":                 "frictionPressure",
    "hidro":                "hydrostaticPressure",
    "cpgas":                "cpGas",
    "cpliq":                "cpLiquid",
    "masstrans":            "massTransferRate",
    "mgFonte":              "sourceGasMassFlow",
    "mlFonte":              "sourceLiqMassFlow",
    "mcFonte":              "sourceCompMassFlow",
    "autoVal":              "waveVelocity",
    "autoVel":              "eigenVelocity",
    "flutuacao":            "fluctuation",
    "temperaturaAmbiente":  "ambientTemperature",
    "yco2":                 "yCO2",
    "RedutorAtrito":        "frictionReducer",
    "Bo":                   "formationVolumeFactor",
    "Froud":                "froudeNumber",
    "GrashExterno":         "externalGrashof",
    "GrashInterno":         "internalGrashof",
    "Hext":                 "externalHeatCoef",
    "Hint":                 "internalHeatCoef",
    "NusselExterno":        "externalNusselt",
    "NusselInterno":        "internalNusselt",
    "PrandtlExterno":       "externalPrandtl",
    "PrandtlGas":           "gasPrandtl",
    "PrandtlInterno":       "internalPrandtl",
    "PrandtlLiquido":       "liquidPrandtl",
    "QGstd":                "stdGasFlowRate",
    "QLWstd":               "stdLiqWaterFlowRate",
    "QLstd":                "stdLiqFlowRate",
    "QLstdTotal":           "stdTotalLiqFlowRate",
    "RGO":                  "gasOilRatioOut",
    "RS":                   "gasOilRatioStd",
    "Rs":                   "gasSolubility",
    "ReyExterno":           "externalReynolds",
    "ReyInterno":           "internalReynolds",
    "TResi":                "residenceTime",
    "Term1":                "term1",
    "Term2":                "term2",
    "VelocidadeMaximaGarganta": "maxThroatVelocity",
    "deng":                 "gasStdDensity",
    "dengD":                "downstreamGasDensity",
    "dengL":                "upstreamGasDensity",
    "deltaPBomba":          "pumpDeltaP",
    "head":                 "pumpHead",
    "potenciaBomba":        "pumpPower",
    "presEstagVGL":         "vglStagePressure",
    "presGargVGL":          "vglThroatPressure",
    "presFundo":            "bottomholePressure",
    "pseudoLiquido":        "pseudoLiquid",
    "rhoMix":               "mixtureDensity",
    "rhog":                 "gasInSituDensity",
    "rhol":                 "liquidInSituDensity",
    "tempChokeJusante":     "chokeDownstreamTemp",
    "tempEstagVGL":         "vglStageTemp",
    "tempGargVGL":          "vglThroatTemp",
    "tempInj":              "injectionTemp",
    "tensaoCisalhamento":   "shearStress",
    "ud":                   "driftVelocity",
    "ug":                   "gasVelocity",
    "ul":                   "liquidVelocity",
    "vazLiq":               "liquidFlowRateSc",
    "vazaoMassicaLiquido":  "liquidMassFlowRate",
    "vazaoVGL":             "vglFlowRate",
    "viscosidadeGas":       "gasViscosity",
    "viscosidadeLiquido":   "liquidViscosity",
    "volJusM1PT":           "m1ProdDownstreamVol",
    "volJusM1ST":           "m1ServDownstreamVol",
    "volMonM1PT":           "m1ProdUpstreamVol",
    "volMonM1ST":           "m1ServUpstreamVol",

    # IPR extra
    "indFluidoPro":         "prodFluidIndex",

    # severeSlugging sub-keys
    "inicioTrechoAcumula":  "accumStartLength",
    "fimTrechoAcumula":     "accumEndLength",
    "fimTrechoColuna":      "columnEndLength",
    "fracaoVazioPenetracao": "voidFractionPenetration",
    "criterio":             "criterion",
}

# ---------------------------------------------------------------------------
# ROOT_PT_TO_EN: overrides for the root-level dict.
# The key difference is "tempo" → "timeSettings" (not "time") at root level.
# ---------------------------------------------------------------------------
ROOT_PT_TO_EN = {**PT_TO_EN, "tempo": "timeSettings"}

# ---------------------------------------------------------------------------
# Value translations for specific Portuguese keys.
# When loading a Portuguese JSON, certain string VALUES also need translation.
# ---------------------------------------------------------------------------
_VALUE_TRANSLATIONS = {
    "tipoMedicaoCamada": {
        "ESPESSURA": "THICKNESS",
        "DIAMETRO":  "DIAMETER",
    },
}


def translate(data, _root=False):
    """Recursively translate Portuguese JSON keys (and some values) to English.

    Args:
        data: A JSON-like structure (dict, list, or scalar).
        _root: Internal flag — set True only for the top-level call so that
               ``"tempo"`` maps to ``"timeSettings"`` instead of ``"time"``.

    Returns:
        A new structure with all Portuguese keys replaced by English equivalents.
    """
    if isinstance(data, dict):
        key_map = ROOT_PT_TO_EN if _root else PT_TO_EN
        result = {}
        for k, v in data.items():
            en_key = key_map.get(k, k)
            # Translate string values for specific keys
            if k in _VALUE_TRANSLATIONS and isinstance(v, str):
                v = _VALUE_TRANSLATIONS[k].get(v.upper(), v)
            result[en_key] = translate(v, _root=False)
        return result
    elif isinstance(data, list):
        return [translate(item, _root=False) for item in data]
    else:
        return data
