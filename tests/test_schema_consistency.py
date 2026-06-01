"""test_schema_consistency.py — Validate consistency between C++ engine,
schema_tramo.json, schema_branch.json, and translations.json.

Checks (all warnings, no errors):
1. Fields read by C++ engine (JSON_entrada.h) not present in schema_tramo.json.
2. Fields in schema_tramo.json not read by C++ engine.
3. Fields in schema_tramo.json without translation in translations.json.
4. Mismatches in schema_branch.json relative to schema_tramo.json + translations.json.
5. Fields without descriptions in both schemas.
"""
import json
import re
import warnings
from pathlib import Path

import pytest

# ============================================================================
# Paths
# ============================================================================
ROOT = Path(__file__).resolve().parent.parent
SCHEMA_TRAMO = ROOT / "docs" / "schema_tramo.json"
SCHEMA_BRANCH = ROOT / "docs" / "schema_branch.json"
TRANSLATIONS = ROOT / "marlim3" / "translations.json"
JSON_ENTRADA_H = ROOT / "src" / "include" / "JSON_entrada.h"


# ============================================================================
# Helpers
# ============================================================================

def load_json(path: Path) -> dict:
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def extract_cpp_fields(header_path: Path) -> set:
    """Parse JSON_entrada.h to extract all Portuguese field names read by the
    C++ engine.

    Strategy: extract field names from class accessor method declarations.
    Every JSON field the engine reads is exposed as a method in a JSONObject
    subclass, e.g.:
        JSON_entrada_configuracaoInicial_saidaClassica& saidaClassica();
        JSON_entrada_configuracaoInicial& configuracaoInicial();
        JSONString& inibidor();

    This captures both primitive fields (#define-based) and object/array fields
    (class-based), and correctly handles field names with underscores like
    nCelulas_XY.
    """
    content = header_path.read_text(encoding="utf-8")

    fields = set()

    # Pattern 1: Class accessor methods — lines like:
    #   JSON_entrada_xxx& fieldName();
    #   JSONString& fieldName();
    # These appear inside class bodies and declare the field accessors.
    # We want the method name (the actual JSON field name).
    accessor_pattern = re.compile(
        r'^\s+'                           # leading whitespace (inside class)
        r'(?:JSON_entrada_\S+|JSONString|JSONBoolean|JSONInteger|JSONNumber)'
        r'&\s+'                           # return type reference
        r'(\w+)'                          # method name = field name
        r'\s*\(\s*\)\s*;',               # ();
        re.MULTILINE
    )
    for m in accessor_pattern.finditer(content):
        fields.add(m.group(1))

    # Pattern 2: Root-level class (JSON_entrada) accessor methods that use
    # the same pattern but may not have JSON_entrada_ prefix in return type
    # Already covered above.

    # Pattern 3: Also extract from the root JSON_entrada class which has methods
    # like: JSON_entrada_sistema& sistema();
    # Already covered by pattern 1.

    return fields


def extract_schema_fields_recursive(schema_obj: dict, path: str = "") -> dict:
    """Recursively extract all property names from a JSON Schema object.
    Returns a dict mapping field_name -> list of paths where it appears.
    """
    result = {}
    if not isinstance(schema_obj, dict):
        return result

    properties = schema_obj.get("properties", {})
    for key, value in properties.items():
        full_path = f"{path}.{key}" if path else key
        result.setdefault(key, []).append(full_path)
        # Recurse into nested objects
        if isinstance(value, dict):
            if value.get("type") == "object":
                nested = extract_schema_fields_recursive(value, full_path)
                for k, v in nested.items():
                    result.setdefault(k, []).extend(v)
            elif value.get("type") == "array":
                items = value.get("items", {})
                if isinstance(items, dict):
                    nested = extract_schema_fields_recursive(items, full_path + "[]")
                    for k, v in nested.items():
                        result.setdefault(k, []).extend(v)
            # Also handle objects without explicit "type": "object" but with "properties"
            elif "properties" in value:
                nested = extract_schema_fields_recursive(value, full_path)
                for k, v in nested.items():
                    result.setdefault(k, []).extend(v)
    return result


def extract_schema_field_set(schema: dict) -> set:
    """Get the flat set of all field names in a schema."""
    return set(extract_schema_fields_recursive(schema).keys())


def extract_fields_without_description(schema_obj: dict, path: str = "") -> list:
    """Recursively find fields that lack a 'description' key or have empty description."""
    missing = []
    if not isinstance(schema_obj, dict):
        return missing

    properties = schema_obj.get("properties", {})
    for key, value in properties.items():
        full_path = f"{path}.{key}" if path else key
        if isinstance(value, dict) and (
            "description" not in value or value.get("description", "").strip() == ""
        ):
            missing.append(full_path)
        # Recurse
        if isinstance(value, dict):
            if value.get("type") == "object":
                missing.extend(extract_fields_without_description(value, full_path))
            elif value.get("type") == "array":
                items = value.get("items", {})
                if isinstance(items, dict):
                    missing.extend(
                        extract_fields_without_description(items, full_path + "[]")
                    )
            elif "properties" in value:
                missing.extend(extract_fields_without_description(value, full_path))
    return missing


# ============================================================================
# Fixtures
# ============================================================================

@pytest.fixture(scope="module")
def schema_tramo():
    return load_json(SCHEMA_TRAMO)


@pytest.fixture(scope="module")
def schema_branch():
    return load_json(SCHEMA_BRANCH)


@pytest.fixture(scope="module")
def translations():
    return load_json(TRANSLATIONS)


@pytest.fixture(scope="module")
def cpp_fields():
    return extract_cpp_fields(JSON_ENTRADA_H)


@pytest.fixture(scope="module")
def schema_tramo_fields(schema_tramo):
    return extract_schema_field_set(schema_tramo)


@pytest.fixture(scope="module")
def schema_branch_fields(schema_branch):
    return extract_schema_field_set(schema_branch)


# ============================================================================
# Tests
# ============================================================================

class TestCppVsSchema:
    """Check consistency between C++ engine fields and schema_tramo.json."""

    def test_cpp_fields_missing_from_schema(self, cpp_fields, schema_tramo_fields):
        """Warn about fields read by C++ engine but not in schema_tramo.json."""
        # Known exceptions: internal fields that are not user-facing
        known_exceptions = {
            "versaoJson",  # version metadata, not a user field
        }
        missing = cpp_fields - schema_tramo_fields - known_exceptions
        if missing:
            warnings.warn(
                f"Fields read by C++ engine but NOT in schema_tramo.json "
                f"({len(missing)}):\n  " + "\n  ".join(sorted(missing)),
                stacklevel=1,
            )

    def test_schema_fields_not_read_by_cpp(self, cpp_fields, schema_tramo_fields):
        """Warn about fields in schema_tramo.json not read by C++ engine."""
        # Known exceptions: metadata or documentation-only fields
        known_exceptions = {
            "$schema",
        }
        not_read = schema_tramo_fields - cpp_fields - known_exceptions
        if not_read:
            warnings.warn(
                f"Fields in schema_tramo.json NOT read by C++ engine "
                f"({len(not_read)}):\n  " + "\n  ".join(sorted(not_read)),
                stacklevel=1,
            )


class TestTranslationsCoverage:
    """Check that translations.json covers all schema_tramo.json fields."""

    def test_schema_fields_untranslated(self, schema_tramo_fields, translations):
        """Warn about schema_tramo.json fields not found in translations.json."""
        # translations["keys"] maps EN -> PT; we need the PT values to check coverage
        pt_values = set(translations["keys"].values())
        # Also include EN keys themselves (some fields are identity-mapped)
        en_keys = set(translations["keys"].keys())

        # A schema_tramo field is covered if it appears as a PT value in translations
        untranslated = set()
        for field in schema_tramo_fields:
            # Field is in PT, check if any translation maps to it
            if field not in pt_values and field not in en_keys:
                untranslated.add(field)

        if untranslated:
            warnings.warn(
                f"Fields in schema_tramo.json WITHOUT translation in "
                f"translations.json ({len(untranslated)}):\n  "
                + "\n  ".join(sorted(untranslated)),
                stacklevel=1,
            )


class TestSchemaBranchConsistency:
    """Check schema_branch.json consistency with schema_tramo.json + translations."""

    def test_branch_fields_match_translations(
        self, schema_branch_fields, schema_tramo_fields, translations
    ):
        """Warn about schema_branch.json fields that don't correspond to any
        English key in translations.json or that map to a PT key not in
        schema_tramo.json.
        """
        en_to_pt = translations["keys"]
        en_keys = set(en_to_pt.keys())

        # Fields in schema_branch that are not recognized EN keys
        unknown_en_fields = schema_branch_fields - en_keys
        # Filter out JSON-Schema meta keys
        meta_keys = {"$schema"}
        unknown_en_fields -= meta_keys

        if unknown_en_fields:
            warnings.warn(
                f"Fields in schema_branch.json NOT recognized as English keys "
                f"in translations.json ({len(unknown_en_fields)}):\n  "
                + "\n  ".join(sorted(unknown_en_fields)),
                stacklevel=1,
            )

    def test_branch_translated_fields_in_schema_tramo(
        self, schema_branch_fields, schema_tramo_fields, translations
    ):
        """Warn about EN fields in schema_branch whose PT equivalent is missing
        from schema_tramo.json.
        """
        en_to_pt = translations["keys"]
        missing_pt = []
        for en_field in schema_branch_fields:
            if en_field in en_to_pt:
                pt_field = en_to_pt[en_field]
                if pt_field not in schema_tramo_fields:
                    missing_pt.append(f"{en_field} -> {pt_field}")

        if missing_pt:
            warnings.warn(
                f"Fields in schema_branch.json whose PT translation is NOT in "
                f"schema_tramo.json ({len(missing_pt)}):\n  "
                + "\n  ".join(sorted(missing_pt)),
                stacklevel=1,
            )

    def test_schema_tramo_fields_in_branch(
        self, schema_branch_fields, schema_tramo_fields, translations
    ):
        """Warn about schema_tramo.json fields that have a translation but are
        missing from schema_branch.json.
        """
        pt_to_en = {v: k for k, v in translations["keys"].items()}

        missing_from_branch = []
        for pt_field in schema_tramo_fields:
            if pt_field in pt_to_en:
                en_field = pt_to_en[pt_field]
                if en_field not in schema_branch_fields:
                    missing_from_branch.append(f"{pt_field} -> {en_field}")

        if missing_from_branch:
            warnings.warn(
                f"Fields in schema_tramo.json (with translation) NOT present in "
                f"schema_branch.json ({len(missing_from_branch)}):\n  "
                + "\n  ".join(sorted(missing_from_branch)),
                stacklevel=1,
            )


class TestDescriptions:
    """Check that all schema fields have descriptions."""

    def test_schema_tramo_descriptions(self, schema_tramo):
        """Warn about fields without description in schema_tramo.json."""
        missing = extract_fields_without_description(schema_tramo)
        if missing:
            warnings.warn(
                f"Fields WITHOUT description in schema_tramo.json "
                f"({len(missing)}):\n  " + "\n  ".join(sorted(missing)),
                stacklevel=1,
            )

    def test_schema_branch_descriptions(self, schema_branch):
        """Warn about fields without description in schema_branch.json."""
        missing = extract_fields_without_description(schema_branch)
        if missing:
            warnings.warn(
                f"Fields WITHOUT description in schema_branch.json "
                f"({len(missing)}):\n  " + "\n  ".join(sorted(missing)),
                stacklevel=1,
            )
