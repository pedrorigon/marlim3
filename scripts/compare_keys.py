"""Compare existing _keys.py PT_TO_EN with inverted translations.json."""
import json, re
from pathlib import Path

# Load translations.json (EN->PT)
t = json.loads(Path("marlim3/translations.json").read_text(encoding="utf-8"))
en_to_pt = t["keys"]

# Invert to PT->EN (first EN key wins for duplicate PT values)
pt_to_en_from_json = {}
for en, pt in en_to_pt.items():
    if pt not in pt_to_en_from_json:
        pt_to_en_from_json[pt] = en

# Parse existing Python PT_TO_EN using regex (safer than ast.literal_eval)
content = Path("marlim3/_tramo/_keys.py").read_text(encoding="utf-8")
# Find entries: "pt_key": "en_key"
# Only look between 'PT_TO_EN = {' and the closing '}'
start = content.find("PT_TO_EN = {")
# Find the end by looking for '\n}' at start of line
end = content.find("\n}", start)
block = content[start:end]
pairs = re.findall(r'"([^"]+)"\s*:\s*"([^"]+)"', block)
existing = dict(pairs)

print(f"Existing Python PT_TO_EN: {len(existing)} entries")
print(f"Inverted JSON PT_TO_EN:   {len(pt_to_en_from_json)} entries")

# Find entries in Python not derivable from JSON
only_in_py = {k: v for k, v in existing.items() if k not in pt_to_en_from_json}
print(f"\nEntries only in Python (legacy aliases): {len(only_in_py)}")
for k, v in sorted(only_in_py.items()):
    print(f'  "{k}": "{v}"')

# Find entries in JSON not in Python
only_in_json = {k: v for k, v in pt_to_en_from_json.items() if k not in existing}
print(f"\nEntries only in JSON (not in Python): {len(only_in_json)}")
for k, v in sorted(only_in_json.items()):
    print(f'  "{k}": "{v}"')

# Find entries where both exist but values differ
conflicts = {k: (existing[k], pt_to_en_from_json[k])
             for k in existing if k in pt_to_en_from_json and existing[k] != pt_to_en_from_json[k]}
print(f"\nConflicts (same PT key, different EN values): {len(conflicts)}")
for k, (py_v, json_v) in sorted(conflicts.items()):
    print(f'  "{k}": Python="{py_v}" vs JSON="{json_v}"')
