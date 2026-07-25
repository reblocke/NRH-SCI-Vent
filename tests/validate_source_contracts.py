#!/usr/bin/env python3
"""Validate the public NRH-003 contracts without reading restricted data."""

from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

SCHEMA_FIELDS = [
    "source_field",
    "required",
    "storage_type",
    "semantic_type",
    "allowed_values",
    "missing_allowed",
    "sensitivity",
    "used_by",
    "notes",
]
RULE_FIELDS = [
    "rule_id",
    "source_field",
    "rule_type",
    "parameters",
    "strict_action",
    "development_action",
    "disclosure_output",
    "notes",
]

STORAGE_TYPES = {"numeric", "string"}
SEMANTIC_TYPES = {
    "identifier",
    "integer",
    "date_mdy",
    "categorical",
    "numeric_with_exceptions",
    "free_text",
    "empty_placeholder",
}
SENSITIVITY_CLASSES = {
    "direct_identifier",
    "quasi_identifier",
    "sensitive_clinical",
    "restricted_free_text",
    "operational",
}
RULE_TYPES = {
    "unique_nonmissing",
    "parseable_date",
    "integer",
    "nonnegative",
    "all_missing",
    "numeric_or_allowed_literal",
}


def read_contract(name: str, expected_fields: list[str]) -> list[dict[str, str]]:
    with (ROOT / "validation" / name).open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        assert reader.fieldnames == expected_fields, (name, reader.fieldnames)
        rows = list(reader)
    assert rows, f"{name} must not be empty"
    assert all(None not in row for row in rows), f"{name} has an over-wide row"
    return rows


def assert_unique_nonblank(rows: list[dict[str, str]], key: str) -> None:
    values = [row[key].strip() for row in rows]
    assert all(values), f"{key} must be nonblank"
    assert len(values) == len(set(values)), f"{key} must be unique"


def validate_schema(rows: list[dict[str, str]]) -> None:
    assert len(rows) == 36
    assert_unique_nonblank(rows, "source_field")
    assert all(row["required"] == "true" for row in rows)
    assert all(row["missing_allowed"] in {"true", "false"} for row in rows)
    assert all(row["storage_type"] in STORAGE_TYPES for row in rows)
    assert all(row["semantic_type"] in SEMANTIC_TYPES for row in rows)
    assert all(row["sensitivity"] in SENSITIVITY_CLASSES for row in rows)

    for row in rows:
        numeric_semantics = {"identifier", "integer", "empty_placeholder"}
        expected_storage = (
            "numeric" if row["semantic_type"] in numeric_semantics else "string"
        )
        assert row["storage_type"] == expected_storage
        allowed = row["allowed_values"]
        if row["semantic_type"] == "categorical":
            assert allowed
        if row["semantic_type"] == "free_text":
            assert not allowed
        if allowed:
            tokens = allowed.split(" | ")
            assert all(token and "|" not in token for token in tokens)
            assert len(tokens) == len(set(tokens))


def validate_rules(
    rows: list[dict[str, str]], schema_rows: list[dict[str, str]]
) -> None:
    assert_unique_nonblank(rows, "rule_id")
    fields = {row["source_field"] for row in schema_rows}
    assert all(row["source_field"] in fields for row in rows)
    assert all(row["rule_type"] in RULE_TYPES for row in rows)
    assert all(row["strict_action"] == "fail" for row in rows)
    assert all(row["development_action"] in {"fail", "warn"} for row in rows)
    assert all(row["disclosure_output"] == "aggregate_count_only" for row in rows)
    for row in rows:
        expected_parameters = {
            "parseable_date": "mask=MD20Y",
            "numeric_or_allowed_literal": "nonnegative=true",
        }.get(row["rule_type"], "")
        assert row["parameters"] == expected_parameters


def assert_public_boundary(paths: list[Path]) -> None:
    prohibited = (
        "/" + "Users/",
        "\\" + "Users\\",
        "source_" + "sha256",
        "source_" + "row_count",
        "aggregate_" + "missing_count",
        "raw_" + "value",
        "PAT_MRN_ID=",
    )
    text = "\n".join(path.read_text(encoding="utf-8") for path in paths)
    assert not any(token in text for token in prohibited)


def main() -> None:
    schema = read_contract("source_schema.csv", SCHEMA_FIELDS)
    rules = read_contract("validation_rules.csv", RULE_FIELDS)
    validate_schema(schema)
    validate_rules(rules, schema)

    # Mutation checks prove the structural validator rejects common bad edits.
    duplicate_schema = [dict(row) for row in schema]
    duplicate_schema[1]["source_field"] = duplicate_schema[0]["source_field"]
    try:
        validate_schema(duplicate_schema)
    except AssertionError:
        pass
    else:
        raise AssertionError("duplicate source fields must fail validation")

    invalid_rules = [dict(row) for row in rules]
    invalid_rules[0]["rule_type"] = "unapproved_method"
    try:
        validate_rules(invalid_rules, schema)
    except AssertionError:
        pass
    else:
        raise AssertionError("invalid rule types must fail validation")

    assert_public_boundary(
        [
            ROOT / "validation" / "data_sources.csv",
            ROOT / "validation" / "source_schema.csv",
            ROOT / "validation" / "validation_rules.csv",
        ]
    )
    print("Public NRH-003 source contracts validated.")


if __name__ == "__main__":
    main()
