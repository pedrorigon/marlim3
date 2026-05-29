"""_keys.py — Bilingual (Portuguese ↔ English) JSON key translator.

Used by Branch.from_json() to convert Portuguese-keyed JSON files to
the English key scheme, and by Branch.to_json() to optionally output
in Portuguese.

The C++ engine handles EN→PT translation at simulation time when
``"language": "en"`` is present in the JSON.

Both this module and the C++ JSONKeyTranslator load from the same
source of truth: ``marlim3/translations.json``.
"""
import json
from pathlib import Path

_TRANSLATIONS_FILE = Path(__file__).resolve().parent.parent / "translations.json"


def _load_translations():
    """Load translations.json and build bidirectional maps."""
    with open(_TRANSLATIONS_FILE, encoding="utf-8") as f:
        data = json.load(f)

    en_to_pt = data["keys"]

    # --- PT → EN (for loading Portuguese JSON) ---
    pt_to_en = {}
    for en, pt in en_to_pt.items():
        pt_to_en[pt] = en

    # Legacy aliases: alternative Portuguese spellings accepted in old files.
    pt_to_en["versaoJson"] = "jsonVersion"  # alternate casing of versaoJSON

    # Identity overrides: abbreviated keys that are the same in both languages.
    _IDENTITY_KEYS = [
        "cdPig", "compP", "compT", "correlacoesPorArranjo",
        "difusTerm3DFE", "fracTermMotorEfic", "mapaArranjo",
        "tabG", "tabP", "taxaDespre", "xCoor", "yCoor",
    ]
    for k in _IDENTITY_KEYS:
        pt_to_en[k] = k

    # --- Value translations (both directions) ---
    # PT→EN values (for loading)
    value_translations_pt_en = {}
    for en_key, mapping in data.get("values", {}).items():
        pt_key = en_to_pt.get(en_key, en_key)
        value_translations_pt_en[pt_key] = {v: k for k, v in mapping.items()}

    # EN→PT values (for saving)
    value_translations_en_pt = {}
    for en_key, mapping in data.get("values", {}).items():
        value_translations_en_pt[en_key] = dict(mapping)

    return pt_to_en, en_to_pt, value_translations_pt_en, value_translations_en_pt


# ---------------------------------------------------------------------------
# Module-level maps (loaded once on import)
# ---------------------------------------------------------------------------
PT_TO_EN, EN_TO_PT, _VALUE_TRANSLATIONS, _VALUE_TRANSLATIONS_EN_PT = (
    _load_translations()
)


def translate(data, _root=False):
    """Recursively translate Portuguese JSON keys (and some values) to English.

    Args:
        data: A JSON-like structure (dict, list, or scalar).
        _root: Internal flag — kept for API compatibility.

    Returns:
        A new structure with all Portuguese keys replaced by English equivalents.
    """
    if isinstance(data, dict):
        result = {}
        for k, v in data.items():
            en_key = PT_TO_EN.get(k, k)
            if k in _VALUE_TRANSLATIONS and isinstance(v, str):
                v = _VALUE_TRANSLATIONS[k].get(v.upper(), v)
            result[en_key] = translate(v, _root=False)
        return result
    elif isinstance(data, list):
        return [translate(item, _root=False) for item in data]
    else:
        return data


def translate_en_to_pt(data):
    """Recursively translate English JSON keys (and enum values) to Portuguese.

    Args:
        data: A JSON-like structure (dict, list, or scalar).

    Returns:
        A new structure with English keys replaced by Portuguese equivalents.
    """
    if isinstance(data, dict):
        result = {}
        for k, v in data.items():
            pt_key = EN_TO_PT.get(k, k)
            if k in _VALUE_TRANSLATIONS_EN_PT and isinstance(v, str):
                v = _VALUE_TRANSLATIONS_EN_PT[k].get(v, v)
            result[pt_key] = translate_en_to_pt(v)
        return result
    elif isinstance(data, list):
        return [translate_en_to_pt(item) for item in data]
    else:
        return data
