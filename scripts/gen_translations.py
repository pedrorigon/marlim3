"""Extract EN→PT translations from JSONKeyTranslator.cpp into translations.json.

Run this ONCE to create the initial marlim3/translations.json.
After that, edit translations.json directly (single source of truth).
"""
import re
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CPP_FILE = ROOT / "src" / "core" / "JSONKeyTranslator.cpp"
OUT_FILE = ROOT / "marlim3" / "translations.json"

content = CPP_FILE.read_text(encoding="utf-8")

# --- Extract key translations ---
pattern = r'\{\s*"([^"]+)"\s*,\s*"([^"]+)"\s*\}'
start = content.find("EN_TO_PT = {")
end = content.find("}; // EN_TO_PT")
block = content[start:end]
pairs = re.findall(pattern, block)

keys_map = {}
for en, pt in pairs:
    if en in keys_map:
        print(f"  WARNING: duplicate EN key '{en}' -> '{keys_map[en]}' vs '{pt}'")
    keys_map[en] = pt

# --- Extract value translations ---
start2 = content.find("VALUE_TRANSLATIONS = {")
end2 = content.find("};", start2)
val_block = content[start2 : end2 + 2]

outer_pattern = r'\{\s*"([^"]+)"\s*,\s*\{(.*?)\}\s*\}'
val_matches = re.findall(outer_pattern, val_block, re.DOTALL)
values_map = {}
for key, inner in val_matches:
    inner_pairs = re.findall(r'\{\s*"([^"]+)"\s*,\s*"([^"]+)"\s*\}', inner)
    values_map[key] = {k: v for k, v in inner_pairs}

# --- Write JSON ---
result = {"keys": keys_map, "values": values_map}
OUT_FILE.write_text(
    json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
)
print(f"Extracted {len(keys_map)} key translations, {len(values_map)} value translations")
print(f"Written to {OUT_FILE}")
