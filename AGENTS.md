# AGENTS

## Project purpose
This repository contains Stata code to reproduce the analyses and figures for ventilator-weaning outcomes in adults with cervical spinal cord injury during inpatient rehabilitation at the University of Utah Craig H. Neilsen Rehabilitation Hospital.

## Data-access constraints
- The source data are restricted retrospective EHR data governed by University of Utah IRB #00153003.
- Do not add PHI, row-level clinical data, cleaned analysis datasets, generated logs containing identifiers, or exported intermediate datasets to this repository.
- Public documentation should link to the final article DOI and repository release metadata rather than copying publisher-formatted article text.

## How to orient quickly
- Start with `README.md` for citation, data-access, workflow, and paper-to-code mapping.
- Use `llms.txt` for the shortest machine-readable project summary.
- Use `CITATION.cff` for citation metadata.
- Use `data_dictionary.md` or `data_dictionary.csv` for expected private inputs, key variables, derived outcomes, and generated output artifacts.
- Use `PROJECT.yml` and `validation/` for the authorized-baseline metadata and public value-free baseline and source contracts.
- Use `vendor/stata/manifest.csv` and `vendor/stata/LICENSES.md` for exact
  community-package files, hashes, sources, and third-party licenses.
- `code/00_preflight.do` is the required offline dependency gate.
- `run_all.do` is the canonical orchestration entry point.
- The main scripts are:
  - `NRH SCI Cohort Preprocessing.do`
  - `NRH SCI Cohort Paper Analysis.do`
  - `NRH SCI Cohort Supplemental Sensitivity Analyses.do`

## Reproduction workflow
Run canonical no-data public verification from the repository root:

```bash
scripts/run_verify.sh
# Windows PowerShell:
.\scripts\run_verify.ps1
```

The launchers call `run_all.do verify` and translate Stata's return code into the process status. Verify runs dependency preflight from `vendor/stata/plus/`, checks public contract structure, and executes atomic approved-literal mapping tests. It does not accept a data directory, open a source CSV, or invoke preprocessing or analyses. Its manifest must record `run_scope=no_data`, `data_accessed=false`, and every data-consuming stage as `not_run`. The public smoke profile is retired because no synthetic cohort is authorized; do not add a fixture or reinterpret verify as an end-to-end analysis test without a new explicit governance decision. Set `STATA_BIN` for a different executable and `STATA_BATCH_FLAG` for a platform-specific batch flag.

The canonical argument order is `profile`, `data_dir`, `output_root`, and optional `run_id`. Use only generated or opaque non-sensitive run IDs. `full` and `release` require an explicit approved restricted-data directory containing `Working NRH SCI.csv`; neither may infer or substitute another input. `release` also refuses to overwrite existing generated `.dta` files:

```bash
stata-mp -b do run_all.do full "/approved/data" "/approved/output"
stata-mp -b do run_all.do release "/approved/data" "/approved/output" "path-safe-run-id"
```

On macOS, if Stata is not on `PATH`, configure the launcher:

```bash
STATA_BIN="/Applications/Stata/StataBE.app/Contents/MacOS/StataBE" \
STATA_BATCH_FLAG=-e scripts/run_verify.sh
```

Direct execution remains available for legacy compatibility:

```bash
stata-mp -b do "NRH SCI Cohort Preprocessing.do"
stata-mp -b do "NRH SCI Cohort Paper Analysis.do"
stata-mp -b do "NRH SCI Cohort Supplemental Sensitivity Analyses.do"
```

## Expected outputs
Restricted canonical runs write generated artifacts to one unique `Results and Figures/<run_id>/` directory. Verify writes value-free evidence under `Verification/<run_id>/`. `run_all.log` and `run_manifest.csv` are value-free orchestration evidence; detailed child logs are under `Logs/`. Direct two-argument legacy calls retain `Results and Figures/<date>/`. These are outputs, not source files, and should not be committed unless a future release explicitly calls for public derived artifacts.

## Documentation surfaces to keep synchronized
- `README.md`: human-readable overview and run instructions.
- `llms.txt`: compact machine-readable summary for indexing.
- `CITATION.cff`: software and final article citation metadata.
- `data_dictionary.md` and `data_dictionary.csv`: expected inputs, derived variables, and output artifacts.
- `PROJECT.yml`, `vendor/stata/manifest.csv`, and `vendor/stata/LICENSES.md`:
  development contracts, dependency hashes, provenance, and licenses.

## Verification before publishing changes
- Run `git diff --check`.
- Validate `CITATION.cff` as YAML after citation edits.
- Search for hard-coded local absolute paths before publishing.
- Confirm `Data/`, `.dta`, Stata logs, `.gph`, and generated output folders remain ignored.
- Exercise invalid profiles, missing restricted inputs, wrong working directories, unsafe and duplicate run IDs, and unwritable output roots without using restricted data.
- Confirm `verify` rejects every data-directory argument and never invokes a data-consuming stage.
- Confirm `full` and `release` reject the retired fixture location and that failed runs never substitute another input.
- Confirm every failed run admitted past directory and writability setup has a value-free top-level log and manifest with a nonzero overall status.
- Run `tests/test_public_contracts.do` and confirm both Stata contract libraries load and every public CSV contract passes structural and cross-reference checks without source data.
- Run `tests/test_value_mappings.do` only as atomic public-literal code tests; do not combine fields into patient- or cohort-shaped fabricated records.
- Do not claim public execution coverage for row-level source-content mutations. Strict content enforcement remains active in authorized `full` and `release` runs, and any new public row-shaped test requires a governance decision.
- Confirm source-contract logs contain only rule IDs, sanitized field names, status, and aggregate invalid counts, with no values, identifiers, row numbers, or source paths.
- Confirm dependency preflight runs before the source-file check, resolves only
  from `vendor/stata/plus/`, and detects missing or modified files.
- Recompute every `vendor/stata/manifest.csv` SHA-256 independently and reject
  runtime `ssc install`, `net install`, or equivalent network installation in
  executable do-files.
- Do not update, replace, or add a vendored dependency without exact upstream
  provenance, redistribution review, a new manifest hash, and focused Stata
  verification in an isolated environment.
- If Stata verification is run, inspect and then remove transient root-level `.log` files before committing.
