"""
_branch.py — English-API mirror of :class:`Tramo`.

``Branch`` stores all configuration under English attribute names.  When
serialised to JSON (via :meth:`to_json` or :meth:`simular`) the key
``"language": "en"`` is written at the JSON root, which instructs the Marlim3
C++ engine to translate every English JSON key to its Portuguese equivalent
before schema validation and simulation.

Inner-object keys (e.g. keys inside a duct dict or a fluid dict) must also
use English names; they are translated by the same mechanism.

Inherited methods (:meth:`simular`, :meth:`to_json`, :meth:`plotar_geometria`,
:meth:`from_mr2`, …) work unchanged because Portuguese attribute names are
transparently redirected to their English counterparts via ``__getattr__``
and ``__setattr__``.
"""

import json as _json
from ._tramo import Tramo


# ---------------------------------------------------------------------------
# Reverse attribute-name mapping: Portuguese Tramo name -> English Branch name.
# Only attributes that differ between the two classes are listed here.
# Language-neutral names (bcs, ipr, pig, material, etc.) are omitted.
# ---------------------------------------------------------------------------
_PT_TO_EN: dict = {
    "sistema":              "system",
    "versaoJson":           "jsonVersion",
    "configuracaoInicial":  "initialConfig",
    "tabela":               "pvtTable",
    "tempo":                "timeSettings",
    "correcao":             "correction",
    "fluidosProducao":      "productionFluids",
    "fluidoComplementar":   "complementaryFluid",
    "fluidoGas":            "gasFluid",
    "secaoTransversal":     "crossSection",
    "dutosProducao":        "productionDucts",
    "dutosServico":         "serviceDucts",
    "separador":            "separator",
    "gasInj":               "gasInjection",
    "CondicaoContPocInjec": "injectionWellBC",
    "fonteMassa":           "massSource",
    "fonteGas":             "gasSource",
    "fonteLiquido":         "liquidSource",
    "fontePressao":         "pressureSource",
    "fontePoroRadial":      "porousRadialSource",
    "fontePoro2D":          "porous2DSource",
    "valvula":              "valve",
    "fonteGasLift":         "gasLiftSource",
    "chokeSup":             "surfaceChoke",
    "chokeInj":             "injectionChoke",
    "master1":              "masterValve",
    "master2":              "masterValve2",
    "bombaVolumetrica":     "volumetricPump",
    "deltaPressao":         "pressureDrop",
    "parafina":             "wax",
    "perfilProducao":       "productionProfile",
    "perfilServico":        "serviceProfile",
    "tendP":                "productionTrend",
    "tendS":                "serviceTrend",
    "perfisTransP":         "transientProductionProfiles",
    "perfisTransS":         "transientServiceProfiles",
    "tendTransP":           "transientProductionTrend",
    "tendTransS":           "transientServiceTrend",
    "tela":                 "screenConfig",
}


class Branch(Tramo):
    """English-API mirror of :class:`Tramo`.

    Usage::

        import marlim3

        b = marlim3.Branch()
        b.productionFluids = [{"id": 0, "api": 32, "gor": 100,
                                "gasDensity": 0.7, "bsw": 0.0}]
        b.material         = [{"id": 0, "type": 0, "conductivity": 58,
                                "specificHeat": 480, "rho": 7850}]
        # … configure cross-sections, ducts, sources, etc. …
        b.simular()
    """

    # ------------------------------------------------------------------
    # Transparent Portuguese <-> English attribute bridging
    # ------------------------------------------------------------------

    def __getattr__(self, name: str):
        """Redirect Portuguese attribute reads to their English counterpart.

        Only called when normal lookup (__dict__ + class chain) fails,
        i.e. only for Portuguese names that do not exist as real attributes.
        """
        en_name = _PT_TO_EN.get(name)
        if en_name is not None:
            return object.__getattribute__(self, en_name)
        raise AttributeError(f"'Branch' object has no attribute '{name}'")

    def __setattr__(self, name: str, value):
        """Redirect Portuguese attribute writes to their English counterpart."""
        en_name = _PT_TO_EN.get(name)
        object.__setattr__(self, en_name if en_name is not None else name, value)

    # ------------------------------------------------------------------
    # Constructor
    # ------------------------------------------------------------------

    def __init__(self,
                 system='MULTIFASICO',
                 jsonVersion='1.3.9',
                 initialConfig=None,
                 pvtTable=None,
                 timeSettings=None,
                 correction=None,
                 productionFluids=None,
                 complementaryFluid=None,
                 gasFluid=None,
                 material=None,
                 crossSection=None,
                 productionDucts=None,
                 serviceDucts=None,
                 separator=None,
                 gasInjection=None,
                 injectionWellBC=None,
                 ipr=None,
                 massSource=None,
                 gasSource=None,
                 liquidSource=None,
                 pressureSource=None,
                 porousRadialSource=None,
                 porous2DSource=None,
                 valve=None,
                 gasLiftSource=None,
                 surfaceChoke=None,
                 injectionChoke=None,
                 masterValve=None,
                 masterValve2=None,
                 bcs=None,
                 volumetricPump=None,
                 pressureDrop=None,
                 pig=None,
                 wax=None,
                 productionProfile=None,
                 serviceProfile=None,
                 productionTrend=None,
                 serviceTrend=None,
                 transientProductionProfiles=None,
                 transientServiceProfiles=None,
                 transientProductionTrend=None,
                 transientServiceTrend=None,
                 screenConfig=None,
                 name=None):

        # "language": "en" triggers C++ key translation on every JSON read.
        # Must be set BEFORE json_entrada_keys is captured.
        self.language = "en"

        self.system          = system
        self.jsonVersion     = jsonVersion

        # general config
        self.initialConfig   = initialConfig   if initialConfig   is not None else {}
        self.pvtTable        = pvtTable        if pvtTable        is not None else {}
        self.timeSettings    = timeSettings
        self.correction      = correction

        # fluids
        self.productionFluids    = productionFluids    if productionFluids    is not None else []
        self.complementaryFluid  = complementaryFluid
        self.gasFluid            = gasFluid

        # geometry
        self.material        = material        if material        is not None else []
        self.crossSection    = crossSection    if crossSection    is not None else []
        self.productionDucts = productionDucts if productionDucts is not None else []
        self.serviceDucts    = serviceDucts    if serviceDucts    is not None else []

        # boundary conditions / sinks+sources
        self.separator       = separator
        self.gasInjection    = gasInjection
        self.injectionWellBC = injectionWellBC
        self.ipr             = ipr             if ipr             is not None else []
        self.massSource      = massSource      if massSource      is not None else []
        self.gasSource       = gasSource       if gasSource       is not None else []
        self.liquidSource    = liquidSource    if liquidSource    is not None else []
        self.pressureSource  = pressureSource  if pressureSource  is not None else []
        self.porousRadialSource = porousRadialSource if porousRadialSource is not None else []
        self.porous2DSource  = porous2DSource  if porous2DSource  is not None else []

        # valves / chokes
        self.valve           = valve           if valve           is not None else []
        self.gasLiftSource   = gasLiftSource   if gasLiftSource   is not None else []
        self.surfaceChoke    = surfaceChoke
        self.injectionChoke  = injectionChoke
        self.masterValve     = masterValve
        self.masterValve2    = masterValve2

        # equipment
        self.bcs             = bcs             if bcs             is not None else []
        self.volumetricPump  = volumetricPump  if volumetricPump  is not None else []
        self.pressureDrop    = pressureDrop    if pressureDrop    is not None else []
        self.pig             = pig             if pig             is not None else []
        self.wax             = wax             if wax             is not None else {}

        # output profiles / trends
        self.productionProfile          = productionProfile          if productionProfile          is not None else {}
        self.serviceProfile             = serviceProfile             if serviceProfile             is not None else {}
        self.productionTrend            = productionTrend            if productionTrend            is not None else []
        self.serviceTrend               = serviceTrend               if serviceTrend               is not None else []
        self.transientProductionProfiles = transientProductionProfiles if transientProductionProfiles is not None else {}
        self.transientServiceProfiles   = transientServiceProfiles   if transientServiceProfiles   is not None else {}
        self.transientProductionTrend   = transientProductionTrend   if transientProductionTrend   is not None else []
        self.transientServiceTrend      = transientServiceTrend      if transientServiceTrend      is not None else []
        self.screenConfig               = screenConfig               if screenConfig               is not None else []

        # mark which attributes are JSON input keys (everything set above)
        self.json_entrada_keys = set(self.__dict__.keys())

        # non-JSON attributes (set after the snapshot)
        self.resultados  = {}
        self.nome_tramo  = name

    # ------------------------------------------------------------------
    # from_json: handles both English-keyed and Portuguese-keyed JSON.
    # ------------------------------------------------------------------

    def from_json(self, json_input, is_string=False):
        """Load configuration from a JSON file or dict.

        Accepts:
        - English-keyed JSON (``"language": "en"``) — as written by
          :meth:`to_json` when ``language`` is ``"en"``.
        - Portuguese-keyed JSON (legacy Marlim3 / :class:`Tramo` output).
        """
        if is_string:
            data = json_input
        else:
            with open(json_input, 'r') as fh:
                data = _json.load(fh)
            if not hasattr(self, 'label'):
                label = json_input[:-5] if json_input.endswith('.json') else json_input
                self.label = label

        lang = data.get('language', 'pt-br').lower() if isinstance(data, dict) else 'pt-br'
        en = lang == 'en'

        def _get(key_en, key_pt, default=None):
            """Try English key first when reading English JSON, otherwise Portuguese."""
            if en:
                return data.get(key_en, data.get(key_pt, default))
            return data.get(key_pt, data.get(key_en, default))

        self.system          = _get('system',           'sistema',              'MULTIFASICO')
        self.jsonVersion     = _get('jsonVersion',      'versaoJson')
        self.initialConfig   = _get('initialConfig',    'configuracaoInicial')  or {}
        self.pvtTable        = _get('pvtTable',         'tabela')               or {}
        self.timeSettings    = _get('timeSettings',     'tempo')
        self.correction      = _get('correction',       'correcao')

        self.productionFluids   = _get('productionFluids',   'fluidosProducao')   or []
        self.complementaryFluid = _get('complementaryFluid', 'fluidoComplementar')
        self.gasFluid           = _get('gasFluid',           'fluidoGas')

        self.material        = _get('material',         'material')             or []
        self.crossSection    = _get('crossSection',     'secaoTransversal')     or []
        self.productionDucts = _get('productionDucts',  'dutosProducao')        or []
        self.serviceDucts    = _get('serviceDucts',     'dutosServico')         or []

        self.separator       = _get('separator',        'separador')
        self.gasInjection    = _get('gasInjection',     'gasInj')
        self.injectionWellBC = _get('injectionWellBC',  'CondicaoContPocInjec')

        self.ipr             = _get('ipr',              'ipr')                  or []
        self.massSource      = _get('massSource',       'fonteMassa')           or []
        self.gasSource       = _get('gasSource',        'fonteGas')             or []
        self.liquidSource    = _get('liquidSource',     'fonteLiquido')         or []
        self.pressureSource  = _get('pressureSource',   'fontePressao')         or []
        self.porousRadialSource = _get('porousRadialSource', 'fontePoroRadial') or []
        self.porous2DSource  = _get('porous2DSource',   'fontePoro2D')          or []

        self.valve           = _get('valve',            'valvula')              or []
        self.gasLiftSource   = _get('gasLiftSource',    'fonteGasLift')         or []
        self.surfaceChoke    = _get('surfaceChoke',     'chokeSup')
        self.injectionChoke  = _get('injectionChoke',   'chokeInj')
        self.masterValve     = _get('masterValve',      'master1')
        self.masterValve2    = _get('masterValve2',     'master2')

        self.bcs             = _get('bcs',              'bcs')                  or []
        self.volumetricPump  = _get('volumetricPump',   'bombaVolumetrica')     or []
        self.pressureDrop    = _get('pressureDrop',     'deltaPressao')         or []
        self.pig             = _get('pig',              'pig')                  or []
        self.wax             = _get('wax',              'parafina')             or {}

        self.productionProfile           = _get('productionProfile',           'perfilProducao')
        self.serviceProfile              = _get('serviceProfile',              'perfilServico')
        self.productionTrend             = _get('productionTrend',             'tendP')
        self.serviceTrend                = _get('serviceTrend',                'tendS')
        self.transientProductionProfiles = _get('transientProductionProfiles', 'perfisTransP')
        self.transientServiceProfiles    = _get('transientServiceProfiles',    'perfisTransS')
        self.transientProductionTrend    = _get('transientProductionTrend',    'tendTransP')
        self.transientServiceTrend       = _get('transientServiceTrend',       'tendTransS')
        self.screenConfig                = _get('screenConfig',                'tela')
