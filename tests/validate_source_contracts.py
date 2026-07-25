#!/usr/bin/env python3
"""Validate the public NRH-004 source and value-mapping contracts."""

from __future__ import annotations

import csv
import re
import unicodedata
from collections import Counter, defaultdict
from pathlib import Path
from typing import Callable


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
DATA_SOURCE_FIELDS = [
    "source_id",
    "file_name",
    "format",
    "delimiter",
    "encoding",
    "line_endings",
    "header_row",
    "expected_columns",
    "classification",
    "contract_version",
    "approval_status",
    "approval_date",
    "approver_role",
    "secure_review_artifact_id",
    "notes",
]
MAPPING_FIELDS = [
    "mapping_id",
    "source_field",
    "input_kind",
    "normalized_literal",
    "target_variable",
    "action",
    "target_code",
    "target_label",
    "decision_id",
    "unlisted_action",
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
DECISION_IDS = {f"NRH004-D{number:03d}" for number in range(1, 9)}
SOURCE_DECISION_IDS = {f"NRH004-D{number:03d}" for number in range(1, 8)}
MAPPING_DECISION_IDS = {f"NRH004-D{number:03d}" for number in range(2, 8)}
TARGET_VARIABLES = {
    "sex",
    "death",
    "rib_fx",
    "ptx",
    "chest_tube",
    "init_injury_level",
    "wean_pre_trach",
    "wean_pre_trans",
    "pneumonia_prior",
    "level",
    "asia_class",
    "partial_wean_at_admit",
    "wean_during_day",
    "wean_24hr",
    "decannulate",
    "rehab_to_icu",
    "resp_icu_transfer",
    "pna_at_rehab",
    "outside_hospital",
}


def read_contract(name: str, expected_fields: list[str]) -> list[dict[str, str]]:
    with (ROOT / "validation" / name).open(
        encoding="utf-8", newline=""
    ) as handle:
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


def normalize_literal(value: str) -> str:
    return unicodedata.normalize("NFC", value).strip().upper()


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
            assert all(token == normalize_literal(token) for token in tokens)


def validate_data_sources(rows: list[dict[str, str]]) -> None:
    assert len(rows) == 1
    row = rows[0]
    assert row["source_id"] == "working_nrh_sci"
    assert row["file_name"] == "Working NRH SCI.csv"
    assert row["format"] == "csv"
    assert row["delimiter"] == "comma"
    assert row["encoding"] == "UTF-8 with BOM"
    assert row["line_endings"] == "CRLF"
    assert row["header_row"] == "1"
    assert row["expected_columns"] == "36"
    assert row["classification"] == "restricted clinical data"
    assert row["contract_version"] == "2"
    assert row["approval_status"] == "approved"
    assert row["approval_date"] == "2026-07-25"
    assert row["approver_role"] == "scientific owner and data steward"
    assert (
        row["secure_review_artifact_id"]
        == "NRH004-37fdea65-5a0b-414c-8a5b-d9356b044e61"
    )


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


MappingSignature = tuple[str, str, str, str, str, str, str, str]


def expected_mapping_signatures() -> set[MappingSignature]:
    expected: set[MappingSignature] = set()

    def literal(
        source: str,
        target: str,
        value: str,
        code: int,
        label: str,
        decision: str,
    ) -> None:
        expected.add(
            (
                source,
                target,
                "literal",
                value,
                "assign",
                str(code),
                label,
                decision,
            )
        )

    def missing_literal(
        source: str, target: str, value: str, decision: str
    ) -> None:
        expected.add(
            (
                source,
                target,
                "literal",
                value,
                "set_missing",
                "",
                "",
                decision,
            )
        )

    def blank(source: str, target: str, decision: str) -> None:
        expected.add(
            (
                source,
                target,
                "blank",
                "",
                "set_missing",
                "",
                "",
                decision,
            )
        )

    def binary(
        source: str,
        target: str,
        negative_literals: tuple[str, ...],
        affirmative_literals: tuple[str, ...],
        *,
        missing_literals: tuple[str, ...] = (),
        decision: str = "NRH004-D006",
    ) -> None:
        for value in negative_literals:
            literal(source, target, value, 0, "N", decision)
        for value in affirmative_literals:
            literal(source, target, value, 1, "Y", decision)
        for value in missing_literals:
            missing_literal(source, target, value, decision)
        blank(source, target, decision)

    literal("sex", "sex", "F", 0, "Female", "NRH004-D002")
    literal("sex", "sex", "M", 1, "Male", "NRH004-D002")
    literal("sex", "sex", "MALE", 1, "Male", "NRH004-D002")
    blank("sex", "sex", "NRH004-D002")

    binary(
        "chesttubeyorn",
        "chest_tube",
        ("N",),
        (
            "Y",
            "Y. AFTER INTUBATION",
            "Y. AFTER INTUBATION. PLEURAL EFFUSION",
        ),
        decision="NRH004-D003",
    )

    for level_number in range(1, 9):
        level = f"C{level_number}"
        literal(
            "injurylevelreportedpriortorehab",
            "init_injury_level",
            level,
            level_number,
            level,
            "NRH004-D004",
        )
    blank(
        "injurylevelreportedpriortorehab",
        "init_injury_level",
        "NRH004-D004",
    )

    grades = {"A": 1, "B": 2, "C": 3, "D": 4}
    for level_number in range(1, 9):
        level = f"C{level_number}"
        grade_literals = [
            (f"{level} AIS {grade}", grade, code)
            for grade, code in grades.items()
        ]
        grade_literals.extend(
            [
                (f"{level} AIA A", "A", 1),
                (f"{level} CENTRAL CORD", "D", 4),
            ]
        )
        for value, grade, grade_code in grade_literals:
            literal(
                "asiaclassificationatrehab",
                "level",
                value,
                level_number,
                level,
                "NRH004-D004",
            )
            literal(
                "asiaclassificationatrehab",
                "asia_class",
                value,
                grade_code,
                f"AIS {grade}",
                "NRH004-D005",
            )
    blank("asiaclassificationatrehab", "level", "NRH004-D004")
    blank("asiaclassificationatrehab", "asia_class", "NRH004-D005")

    binary("deceased", "death", ("N",), ("Y",))
    binary("ribfracturesyorn", "rib_fx", ("N",), ("Y",))
    binary("pneumothoraxyorn", "ptx", ("N",), ("Y",))
    binary(
        "attempttoweanpriortotrachyorn",
        "wean_pre_trach",
        ("N", "NO"),
        ("Y", "Y - EXTUBATED ON 8/20  BUT PT STRIDORED AND C.A."),
        missing_literals=("UNKNOWN",),
    )
    binary(
        "attempttoweanoffventpriortotrans",
        "wean_pre_trans",
        ("N",),
        ("Y", "Y (ABLE TO WEANE TO VENT AT NIGHT)"),
        missing_literals=("UNKNOWN", "UNKOWN"),
    )
    binary(
        "didpatientdeveloppneumoniapriort",
        "pneumonia_prior",
        ("N",),
        ("Y", "YES"),
    )
    binary("osh", "outside_hospital", ("N",), ("Y",))
    binary(
        "weanoffventforfull24hoursyorn",
        "wean_24hr",
        ("N", "N - BIPAP AT NIGHT"),
        ("Y",),
    )
    binary(
        "didpatientdecanulateyorn",
        "decannulate",
        ("N", "N."),
        ("Y",),
    )
    binary(
        "didtheytransfertoicufromrehabyor",
        "rehab_to_icu",
        ("N",),
        ("Y", "Y  - SICU"),
    )
    binary(
        "iftheytransferredtoicuwascausere",
        "resp_icu_transfer",
        ("N", "N (DIAPHRAM PACER PLACEMENT)"),
        ("Y",),
    )
    binary(
        "didtheydeveloppnumoniaatrehabyor",
        "pna_at_rehab",
        ("N",),
        ("Y",),
    )

    binary(
        "weanoffventduringthedayandcontin",
        "partial_wean_at_admit",
        ("N", "Y", "Y."),
        ("Y (PRIOR TO REHAB)", "Y PRIOR TO ADMISSION"),
        missing_literals=("UNKNOWN", "UNKOWN"),
        decision="NRH004-D007",
    )
    binary(
        "weanoffventduringthedayandcontin",
        "wean_during_day",
        ("N",),
        ("Y", "Y.", "Y (PRIOR TO REHAB)", "Y PRIOR TO ADMISSION"),
        missing_literals=("UNKNOWN", "UNKOWN"),
        decision="NRH004-D006",
    )
    return expected


def mapping_signature(row: dict[str, str]) -> MappingSignature:
    return (
        row["source_field"],
        row["target_variable"],
        row["input_kind"],
        row["normalized_literal"],
        row["action"],
        row["target_code"],
        row["target_label"],
        row["decision_id"],
    )


def validate_mappings(
    rows: list[dict[str, str]], schema_rows: list[dict[str, str]]
) -> None:
    assert_unique_nonblank(rows, "mapping_id")
    schema = {row["source_field"]: row for row in schema_rows}
    assert all(row["source_field"] in schema for row in rows)
    assert all(row["target_variable"] in TARGET_VARIABLES for row in rows)
    assert all(row["input_kind"] in {"literal", "blank"} for row in rows)
    assert all(row["action"] in {"assign", "set_missing"} for row in rows)
    assert all(row["decision_id"] in MAPPING_DECISION_IDS for row in rows)
    assert all(row["unlisted_action"] == "reject" for row in rows)

    keys = [
        (
            row["source_field"],
            row["target_variable"],
            row["input_kind"],
            row["normalized_literal"],
        )
        for row in rows
    ]
    assert len(keys) == len(set(keys)), "mapping lookup keys must be unique"

    pairs = {(row["source_field"], row["target_variable"]) for row in rows}
    blank_pairs = Counter(
        (row["source_field"], row["target_variable"])
        for row in rows
        if row["input_kind"] == "blank"
    )
    assert set(blank_pairs) == pairs
    assert all(count == 1 for count in blank_pairs.values())

    for row in rows:
        if row["input_kind"] == "blank":
            assert row["normalized_literal"] == ""
            assert row["action"] == "set_missing"
        else:
            assert row["normalized_literal"]
            assert row["normalized_literal"] == normalize_literal(
                row["normalized_literal"]
            )
        if row["action"] == "assign":
            assert re.fullmatch(r"-?[0-9]+", row["target_code"])
            assert row["target_label"]
        else:
            assert row["target_code"] == ""
            assert row["target_label"] == ""

    actual = {mapping_signature(row) for row in rows}
    expected = expected_mapping_signatures()
    assert actual == expected, (
        f"mapping coverage differs: missing={len(expected - actual)} "
        f"unexpected={len(actual - expected)}"
    )

    mapped_literals: dict[str, set[str]] = defaultdict(set)
    for row in rows:
        if row["input_kind"] == "literal":
            mapped_literals[row["source_field"]].add(row["normalized_literal"])
    for source_field, literals in mapped_literals.items():
        schema_literals = set(
            schema[source_field]["allowed_values"].split(" | ")
        )
        assert literals == schema_literals, (
            f"{source_field} schema and mapping literal domains differ"
        )

    compound = [
        row
        for row in rows
        if row["source_field"] == "asiaclassificationatrehab"
        and row["input_kind"] == "literal"
    ]
    assert len(compound) == 96
    assert Counter(row["target_variable"] for row in compound) == {
        "level": 48,
        "asia_class": 48,
    }


def top_level_section(text: str, name: str) -> str:
    match = re.search(
        rf"(?ms)^{re.escape(name)}:\n(.*?)(?=^[A-Za-z_][A-Za-z0-9_]*:\n|\Z)",
        text,
    )
    assert match, f"PROJECT.yml lacks {name}"
    return match.group(1)


def validate_project() -> None:
    text = (ROOT / "PROJECT.yml").read_text(encoding="utf-8")
    latest_release = top_level_section(text, "latest_release")
    assert "  version: v0.2.0\n" in latest_release
    assert (
        "  commit: 45a50136bf6b0279b2a1cf34a51c559109596c0a\n"
        in latest_release
    )

    baseline = top_level_section(text, "baseline")
    assert (
        "  authoritative_commit: "
        "6ba1caf49ac9e260723d850b2a9f189f83255258\n"
        in baseline
    )

    analysis_files = top_level_section(text, "analysis_files")
    frozen_hashes = {
        "d220ae6ab659760d7def73010ee985a14e6e2a65ad00c9c9ab309d81a0001a42",
        "7e57f1bf30b8029a67d5a24bf601df4ffab7e0a839298c28d64b53179280e866",
        "73e6b20e0f4e337286ea0b5787644c7e1c40b9452b8a5d5778201fe6a70b4480",
    }
    assert all(value in analysis_files for value in frozen_hashes)

    source_contract = top_level_section(text, "source_contract")
    assert "  version: 2\n" in source_contract
    assert "  approval_history:\n" in source_contract
    assert (
        "NRH003-4fe498dd-ac82-47fd-8677-5a28815a953f"
        in source_contract
    )
    assert (
        "NRH004-37fdea65-5a0b-414c-8a5b-d9356b044e61"
        in source_contract
    )
    assert all(decision in source_contract for decision in SOURCE_DECISION_IDS)
    assert "NRH004-D008" not in source_contract

    mappings = top_level_section(text, "value_mapping_contract")
    assert "  version: 1\n" in mappings
    assert "  path: validation/value_mappings.csv\n" in mappings
    assert "    - assign\n" in mappings
    assert "    - set_missing\n" in mappings
    assert "  unlisted_action: reject\n" in mappings
    assert "    decision_record: DECISIONS.md\n" in mappings
    assert all(decision in mappings for decision in MAPPING_DECISION_IDS)
    assert not (DECISION_IDS - MAPPING_DECISION_IDS) & set(
        re.findall(r"NRH004-D[0-9]{3}", mappings)
    )

    constants = top_level_section(text, "scientific_constants")
    assert "  nrh_admin_censor_date_iso:\n" in constants
    assert '    value: "2023-09-30"\n' in constants
    assert "    decision_id: NRH004-D008\n" in constants
    assert "      - scientific owner\n" in constants
    assert "      - data steward\n" in constants
    assert (
        "    secure_artifact_id: "
        "NRH004-37fdea65-5a0b-414c-8a5b-d9356b044e61\n"
        in constants
    )
    assert (
        "    change_policy: Any change requires a new signed decision "
        "and full baseline review.\n"
        in constants
    )


def validate_preprocessing(
    text: str, mapping_rows: list[dict[str, str]]
) -> None:
    assert not re.search(r"(?mi)^[ \t]*encode[ \t]+", text)

    call_blocks = re.findall(
        r"(?ms)^nrh_apply_value_mapping\b.*?(?=^[A-Za-z]|\Z)",
        text,
    )
    actual_pairs: list[tuple[str, str]] = []
    for block in call_blocks:
        source = re.search(r'sourcefield\("([^"]+)"\)', block)
        target = re.search(r"target\(([A-Za-z0-9_]+)\)", block)
        assert source and target, "malformed preprocessing mapping call"
        actual_pairs.append((source.group(1), target.group(1)))

    expected_pairs = {
        (row["source_field"], row["target_variable"])
        for row in mapping_rows
    }
    assert set(actual_pairs) == expected_pairs
    assert len(actual_pairs) == len(expected_pairs)

    required_derivations = (
        "replace comp_vs_part = 1 if inlist(asia_class, 1, 2)",
        "replace comp_vs_part = 2 if inlist(asia_class, 3, 4)",
        "replace high_vs_low = 1 if inrange(level, 1, 4)",
        "replace high_vs_low = 2 if inrange(level, 5, 8)",
        "replace init_high_vs_low = 1 if inrange(init_injury_level, 1, 4)",
        "replace init_high_vs_low = 2 if inrange(init_injury_level, 5, 8)",
        "replace level_and_completeness = 0 if high_vs_low == 1 & comp_vs_part == 1",
        "replace level_and_completeness = 1 if high_vs_low == 1 & comp_vs_part == 2",
        "replace level_and_completeness = 2 if high_vs_low == 2 & comp_vs_part == 1",
        "replace level_and_completeness = 3 if high_vs_low == 2 & comp_vs_part == 2",
        "drop `nrh_discharge_normalized'",
        'local nrh_admin_censor_date_iso "2023-09-30"',
        'daily("`nrh_admin_censor_date_iso\'", "YMD")',
    )
    assert all(derivation in text for derivation in required_derivations)


def assert_public_boundary(paths: list[Path]) -> None:
    prohibited_tokens = (
        "/" + "Users/",
        "\\" + "Users\\",
        "NRH " + "Secure Baselines",
        "source_" + "sha256",
        "source_" + "row_count",
        "aggregate_" + "missing_count",
        "raw_" + "source_literal",
        "raw_" + "value",
        "PAT_" + "MRN_ID=",
    )
    text = "\n".join(path.read_text(encoding="utf-8") for path in paths)
    assert not any(token in text for token in prohibited_tokens)
    assert not re.search(r"(?i)(?:^|[\\s=(])/(?:Users|home|Volumes)/", text)
    windows_user_root = r"(?i)[A-Z]:\\\\" + "Users" + r"\\\\"
    assert not re.search(windows_user_root, text)


def must_fail(call: Callable[[], None], message: str) -> None:
    try:
        call()
    except AssertionError:
        return
    raise AssertionError(message)


def main() -> None:
    schema = read_contract("source_schema.csv", SCHEMA_FIELDS)
    rules = read_contract("validation_rules.csv", RULE_FIELDS)
    data_sources = read_contract("data_sources.csv", DATA_SOURCE_FIELDS)
    mappings = read_contract("value_mappings.csv", MAPPING_FIELDS)

    validate_schema(schema)
    validate_rules(rules, schema)
    validate_data_sources(data_sources)
    validate_mappings(mappings, schema)
    validate_project()
    preprocessing = (ROOT / "NRH SCI Cohort Preprocessing.do").read_text(
        encoding="utf-8"
    )
    validate_preprocessing(preprocessing, mappings)

    duplicate_schema = [dict(row) for row in schema]
    duplicate_schema[1]["source_field"] = duplicate_schema[0]["source_field"]
    must_fail(
        lambda: validate_schema(duplicate_schema),
        "duplicate source fields must fail validation",
    )

    invalid_rules = [dict(row) for row in rules]
    invalid_rules[0]["rule_type"] = "unapproved_method"
    must_fail(
        lambda: validate_rules(invalid_rules, schema),
        "invalid rule types must fail validation",
    )

    duplicate_mapping = [dict(row) for row in mappings]
    duplicate_mapping[1]["mapping_id"] = duplicate_mapping[0]["mapping_id"]
    must_fail(
        lambda: validate_mappings(duplicate_mapping, schema),
        "duplicate mapping IDs must fail validation",
    )

    missing_blank = [
        dict(row)
        for row in mappings
        if row["mapping_id"] != "NRH004-SEX-BLANK"
    ]
    must_fail(
        lambda: validate_mappings(missing_blank, schema),
        "a missing source-target blank row must fail validation",
    )

    invalid_action = [dict(row) for row in mappings]
    invalid_action[0]["action"] = "coerce"
    must_fail(
        lambda: validate_mappings(invalid_action, schema),
        "an unapproved action must fail validation",
    )

    invalid_decision = [dict(row) for row in mappings]
    invalid_decision[0]["decision_id"] = "NRH004-D999"
    must_fail(
        lambda: validate_mappings(invalid_decision, schema),
        "an unapproved decision ID must fail validation",
    )

    invalid_unlisted_action = [dict(row) for row in mappings]
    invalid_unlisted_action[0]["unlisted_action"] = "set_missing"
    must_fail(
        lambda: validate_mappings(invalid_unlisted_action, schema),
        "an unlisted action other than reject must fail validation",
    )

    unnormalized = [dict(row) for row in mappings]
    unnormalized[0]["normalized_literal"] = " f "
    must_fail(
        lambda: validate_mappings(unnormalized, schema),
        "an unnormalized literal must fail validation",
    )

    changed_code = [dict(row) for row in mappings]
    changed_code[0]["target_code"] = "1"
    must_fail(
        lambda: validate_mappings(changed_code, schema),
        "an unapproved target code must fail validation",
    )

    widened_unknown = [dict(row) for row in mappings]
    extra = dict(widened_unknown[0])
    extra["mapping_id"] = "MUTATION-UNAPPROVED-UNKNOWN"
    extra["source_field"] = "chesttubeyorn"
    extra["target_variable"] = "chest_tube"
    extra["normalized_literal"] = "UNKNOWN"
    extra["action"] = "set_missing"
    extra["target_code"] = ""
    extra["target_label"] = ""
    extra["decision_id"] = "NRH004-D003"
    widened_unknown.append(extra)
    must_fail(
        lambda: validate_mappings(widened_unknown, schema),
        "an undocumented field-specific unknown must fail validation",
    )

    narrowed_schema = [dict(row) for row in schema]
    initial_row = next(
        row
        for row in narrowed_schema
        if row["source_field"] == "injurylevelreportedpriortorehab"
    )
    initial_row["allowed_values"] = initial_row["allowed_values"].replace(
        " | C8", ""
    )
    must_fail(
        lambda: validate_mappings(mappings, narrowed_schema),
        "schema and mapping literal drift must fail validation",
    )

    missing_mapping_call = preprocessing.replace(
        'sourcefield("sex") target(sex)',
        'sourcefield("sex") target(sex_removed)',
        1,
    )
    must_fail(
        lambda: validate_preprocessing(missing_mapping_call, mappings),
        "preprocessing mapping-call drift must fail validation",
    )

    changed_grouping = preprocessing.replace(
        "replace high_vs_low = 2 if inrange(level, 5, 8)",
        "replace high_vs_low = 2 if inrange(level, 5, 7)",
        1,
    )
    must_fail(
        lambda: validate_preprocessing(changed_grouping, mappings),
        "preprocessing C8 grouping drift must fail validation",
    )

    reintroduced_encode = preprocessing + "\nencode sex, generate(male)\n"
    must_fail(
        lambda: validate_preprocessing(reintroduced_encode, mappings),
        "order-dependent encode must fail validation",
    )

    assert_public_boundary(
        [
            ROOT / ".gitignore",
            ROOT / "CHANGELOG.md",
            ROOT / "DECISIONS.md",
            ROOT / "NRH SCI Cohort Preprocessing.do",
            ROOT / "PROJECT.yml",
            ROOT / "README.md",
            ROOT / "code" / "lib" / "data_contracts.do",
            ROOT / "code" / "lib" / "value_mappings.do",
            ROOT / "data_dictionary.csv",
            ROOT / "data_dictionary.md",
            ROOT / "llms.txt",
            ROOT / "run_all.do",
            ROOT / "tests" / "test_source_contract.do",
            ROOT / "tests" / "test_value_mappings.do",
            ROOT / "tests" / "validate_source_contracts.py",
            ROOT / "validation" / "README.md",
            ROOT / "validation" / "data_sources.csv",
            ROOT / "validation" / "source_schema.csv",
            ROOT / "validation" / "validation_rules.csv",
            ROOT / "validation" / "value_mappings.csv",
        ]
    )
    print("Public NRH-004 source and value-mapping contracts validated.")


if __name__ == "__main__":
    main()
