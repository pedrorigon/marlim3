"""_keys.py — Python-side Portuguese → English JSON key translator.

Used by Branch.from_json() to convert old Portuguese-keyed JSON files to
the English key scheme before populating Branch attributes.

The C++ engine handles the reverse direction (EN→PT) at simulation time
when ``"language": "en"`` is present in the JSON.

Both this module and the C++ JSONKeyTranslator load from the same
source of truth: ``marlim3/translations.json``.
"""
import json
from pathlib import Path

_TRANSLATIONS_FILE = Path(__file__).resolve().parent.parent / "translations.json"

def _load_translations():
    """Load translations.json and build the PT→EN maps."""
    with open(_TRANSLATIONS_FILE, encoding="utf-8") as f:
        data = json.load(f)

    # Invert EN→PT to get PT→EN.  For duplicate PT values (multiple EN keys
    # mapping to the same PT key), keep the LAST one encountered — this gives
    # inner/nested keys priority, matching the original hand-written dict order.
    en_to_pt = data["keys"]
    pt_to_en = {}
    for en, pt in en_to_pt.items():
        pt_to_en[pt] = en

    # Legacy aliases: alternative Portuguese spellings accepted in old files.
    pt_to_en["versaoJson"] = "jsonVersion"  # alternate casing of versaoJSON

    # Identity overrides: abbreviated keys that are the same in both languages.
    # The C++ translator maps longer EN names to these, but the Python side
    # keeps them unchanged for backward compatibility with Branch attributes.
    _IDENTITY_KEYS = [
        "cdPig", "compP", "compT", "correlacoesPorArranjo",
        "difusTerm3DFE", "fracTermMotorEfic", "mapaArranjo",
        "tabG", "tabP", "taxaDespre", "xCoor", "yCoor",
    ]
    for k in _IDENTITY_KEYS:
        pt_to_en[k] = k

    # Invert value translations (C++ stores EN→PT values; Python needs PT→EN).
    value_translations = {}
    for en_key, mapping in data.get("values", {}).items():
        # Find the corresponding PT key name.
        pt_key = en_to_pt.get(en_key, en_key)
        value_translations[pt_key] = {v: k for k, v in mapping.items()}

    return pt_to_en, value_translations


# ---------------------------------------------------------------------------
# Module-level maps (loaded once on import)
# ---------------------------------------------------------------------------
PT_TO_EN, _VALUE_TRANSLATIONS = _load_translations()


def translate(data, _root=False):
    """Recursively translate Portuguese JSON keys (and some values) to English.

    Args:
        data: A JSON-like structure (dict, list, or scalar).
        _root: Internal flag — kept for API compatibility but no longer
               needed (no root-level ambiguities remain).

    Returns:
        A new structure with all Portuguese keys replaced by English equivalents.
    """
    if isinstance(data, dict):
        result = {}
        for k, v in data.items():
            en_key = PT_TO_EN.get(k, k)
            # Translate string values for specific keys
            if k in _VALUE_TRANSLATIONS and isinstance(v, str):
                v = _VALUE_TRANSLATIONS[k].get(v.upper(), v)
            result[en_key] = translate(v, _root=False)
        return result
    elif isinstance(data, list):
        return [translate(item, _root=False) for item in data]
    else:
        return data
