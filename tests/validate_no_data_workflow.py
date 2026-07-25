#!/usr/bin/env python3
"""Validate the NRH-005 no-data public-verification boundary."""

from __future__ import annotations

import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

ALLOWED_TRACKED_CSV = {
    "data_dictionary.csv",
    "validation/data_sources.csv",
    "validation/expected_results_schema.csv",
    "validation/sample_flow_schema.csv",
    "validation/source_schema.csv",
    "validation/validation_rules.csv",
    "validation/value_mappings.csv",
    "vendor/stata/manifest.csv",
}

PROHIBITED_TRACKED_SUFFIXES = {
    ".dta",
    ".ster",
    ".log",
    ".smcl",
    ".gph",
    ".tif",
    ".tiff",
    ".png",
    ".xlsx",
}


def read(path: str) -> str:
    return (ROOT / path).read_text(encoding="utf-8")


def tracked_paths() -> list[str]:
    completed = subprocess.run(
        ["git", "ls-files", "--cached", "--others", "--exclude-standard", "-z"],
        cwd=ROOT,
        check=True,
        capture_output=True,
    )
    return [
        value.decode("utf-8")
        for value in completed.stdout.split(b"\0")
        if value
    ]


def validate_tracked_boundary(paths: list[str]) -> None:
    assert all(not path.startswith("data/synthetic/") for path in paths)
    assert all(Path(path).name != "Working NRH SCI.csv" for path in paths)
    assert all(
        Path(path).suffix.lower() not in PROHIBITED_TRACKED_SUFFIXES
        for path in paths
    )
    tracked_csv = {path for path in paths if Path(path).suffix.lower() == ".csv"}
    assert tracked_csv == ALLOWED_TRACKED_CSV, tracked_csv ^ ALLOWED_TRACKED_CSV


def validate_runner() -> None:
    runner = read("run_all.do")
    assert 'local profile "verify"' in runner
    assert 'inlist("`profile\'", "verify", "full", "release")' in runner
    assert "The public smoke profile is retired" in runner
    assert 'local nrh_run_scope "no_data"' in runner
    assert 'data_accessed,`nrh_data_accessed\'' in runner
    assert 'public_contract_status,`nrh_public_contract_status\'' in runner
    assert 'mapping_unit_status,`nrh_mapping_unit_status\'' in runner
    assert 'preflight_status,`nrh_preflight_status\'' in runner
    assert runner.count('if "`profile\'" != "verify" & `nrh_overall_rc\' == 0') == 7

    verify_adapter = read("scripts/run_verify.do")
    assert 'do run_all.do verify ///' in verify_adapter
    assert '"" "`output_root\'" "`run_id\'"' in verify_adapter

    verify_config = read("config/verify.do")
    assert 'local nrh_default_output_root "Verification"' in verify_config


def validate_test_boundary() -> None:
    public_contract_test = read("tests/test_public_contracts.do")
    prohibited_record_builders = (
        "set obs",
        "generate long pat_mrn_id",
        "export delimited",
        "nrh_test_write_fixture",
    )
    assert not any(
        token in public_contract_test for token in prohibited_record_builders
    )
    assert "without source data" in public_contract_test

    mapping_test = read("tests/test_value_mappings.do")
    assert "patient- or cohort-shaped records" in mapping_test
    assert "generate long pat_mrn_id" not in mapping_test


def validate_documentation() -> None:
    synchronized_paths = (
        "README.md",
        "AGENTS.md",
        "CONTRIBUTING.md",
        "llms.txt",
        "PROJECT.yml",
        "CHANGELOG.md",
        "validation/README.md",
    )
    text = "\n".join(read(path) for path in synchronized_paths)
    assert "NRH-005 will add" not in text
    assert "data/synthetic/Working NRH SCI.csv" not in text
    assert "scripts/run_verify.sh" in text
    assert "no-data" in text.lower()


def validate_ci() -> None:
    workflow = read(".github/workflows/public-validation.yml")
    assert "permissions:\n  contents: read\n" in workflow
    assert "tests/validate_source_contracts.py" in workflow
    assert "tests/validate_no_data_workflow.py" in workflow
    assert "tests/test_run_verify.sh" in workflow
    assert "stata-mp" not in workflow.lower()
    assert "statabe" not in workflow.lower()


def main() -> None:
    paths = tracked_paths()
    validate_tracked_boundary(paths)
    validate_runner()
    validate_test_boundary()
    validate_documentation()
    validate_ci()
    print("NRH-005 no-data public-verification boundary validated.")


if __name__ == "__main__":
    main()
