# JSON Key Translation System

Marlim3 uses a bilingual JSON interface: users can write input files with **English or Portuguese keys**, and the engine internally operates with **Portuguese keys**. A single source of truth — `marlim3/translations.json` — drives both directions of translation.

## Architecture Overview

```
┌─────────────────────────────────────────────────────────────────┐
│                   marlim3/translations.json                      │
│                  (single source of truth)                        │
│                   EN key → PT key mapping                        │
└────────────────────────┬────────────────────────────────────────┘
                         │
            ┌────────────┼────────────────┐
            ▼                             ▼
┌───────────────────────┐    ┌────────────────────────────┐
│   C++ (compile-time)  │    │   Python (runtime)         │
│                       │    │                            │
│ CMake reads .json →   │    │ _keys.py loads .json →     │
│ embeds in             │    │ inverts to PT→EN map       │
│ translations_data.h → │    │                            │
│ JSONKeyTranslator.cpp │    │ Branch.from_json() uses    │
│ translates EN→PT at   │    │ PT→EN map to load old      │
│ simulation startup    │    │ Portuguese JSON files       │
└───────────────────────┘    └────────────────────────────┘
```

### C++ side (EN → PT)

When the input JSON contains `"language": "en"` at the root level, the C++ reader (`Leitura.cpp`) invokes `JSONKeyTranslator::translateEnToPt()` to recursively rename all English keys to their Portuguese equivalents before populating the internal data structures. The translation map is embedded at compile time via CMake's `configure_file()`.

### Python side (PT → EN)

Python API works exclusively with English keys. The Python API (`marlim3._tramo._keys`) inverts the translations map to create a PT→EN dictionary. This allows `Branch.from_json()` to load Portuguese-keyed JSON files and present all attributes with English names.

## File Roles

| File | Role |
|------|------|
| `marlim3/translations.json` | Single source of truth for all EN↔PT key mappings |
| `src/include/translations_data.h` | Auto-generated C++ header (do not edit manually) |
| `src/include/translations_data.h.in` | CMake template for the header |
| `src/core/JSONKeyTranslator.cpp` | C++ translation logic (EN→PT at runtime) |
| `marlim3/_tramo/_keys.py` | Python translation logic (PT→EN at import time) |

## How `translations.json` Is Structured

```json
{
  "keys": {
    "englishKeyName": "nomeChavePortugues",
    ...
  },
  "values": {
    "keyName": {
      "ENGLISH_VALUE": "VALOR_PORTUGUES",
      ...
    }
  }
}
```

- **`keys`**: flat EN→PT map for JSON key names. The same English name always maps to the same Portuguese name regardless of nesting depth.
- **`values`**: optional EN→PT map for string *values* of specific keys (e.g., `"layerMeasurementType"` values `"THICKNESS"` ↔ `"ESPESSURA"`).

### Special cases

- **Identity mappings** (EN == PT): some short technical terms are the same in both languages (e.g., `"gasInj": "gasInj"`, `"holdup": "holdup"`). The C++ translator skips renaming when source equals target.
- **No ambiguous PT keys**: each Portuguese key maps to exactly one English key, ensuring lossless round-trip translation.

## Adding a New Field

When you add a new field to the Marlim3 JSON input schema:

1. **Choose the English key name** — use camelCase, be descriptive but concise.

2. **Determine the Portuguese key name** — this is what the C++ engine expects internally. Check existing `JSON_entrada_*.cpp` code for the expected key.

3. **Add the entry to `marlim3/translations.json`**:
   ```json
   "myNewField": "meuNovoCampo"
   ```

4. **Re-run CMake configure** to regenerate `translations_data.h`:
   ```bash
   cmake --preset <your-preset>
   ```
   This embeds the updated JSON into the C++ header. You do **not** need to rebuild the full project — just re-configure.

5. **Verify** — the Python side picks up changes automatically (it reads `translations.json` at import time). No code changes needed in `_keys.py` unless the new field has special semantics.

6. **Update documentation**:
   - Add the key to `docs/schema_branch.json` if it's a branch-level field.
   - Update relevant tutorial notebooks if the field is user-facing.

### If the field has translatable *values*

Add an entry under `"values"`:

```json
"values": {
  "myNewField": {
    "ENGLISH_VALUE": "VALOR_PORTUGUES"
  }
}
```

## Modifying an Existing Translation (Renaming a Key)

Renaming an English key requires coordinated changes across the project:

1. **Update `marlim3/translations.json`** — change the EN key name (left side). The PT value (right side) must remain unchanged (the C++ engine depends on it):
   ```json
   // Before
   "oldEnglishName": "nomePortugues"
   // After
   "newEnglishName": "nomePortugues"
   ```

2. **Re-run CMake configure** to regenerate `translations_data.h`.

3. **Update Python `Branch` class** (`marlim3/_tramo/_branch.py`):
   - Rename the attribute in `__init__()`, `from_json()`, and any methods that reference it (e.g., `display_table()`).
   - The `to_json()` method writes attributes directly as JSON keys, so the attribute name *is* the output key.

4. **Update English demo files** (`demos/*.json` with `"language": "en"`):
   - Replace all occurrences of the old key with the new one.

5. **Update documentation**:
   - `docs/schema_branch.json` — rename the property key and update descriptions.
   - `docs/tutorials/*.ipynb` — update code cells and markdown references.

6. **Run tests** to verify nothing breaks:
   ```bash
   pytest tests/
   ```

### Important constraints

- **Never change the Portuguese value** (right side) unless the C++ engine code is also updated. The C++ `JSON_entrada` classes parse by Portuguese key name.
- **No backward-compatibility layer** — the Python API uses attribute names directly as JSON keys. Old files with renamed keys will fail unless loaded through the PT→EN translator (which handles PT-keyed files automatically).
- **The translation map is flat and unambiguous** — each EN key maps to exactly one PT key and vice-versa. Avoid introducing duplicate PT values (multiple EN keys mapping to the same PT key), as this breaks the Python PT→EN inversion.
