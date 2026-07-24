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
- Use `PROJECT.yml` and `validation/` for the authorized-baseline metadata and public value-free contracts.
- `run_all.do` is the canonical orchestration entry point.
- The main scripts are:
  - `NRH SCI Cohort Preprocessing.do`
  - `NRH SCI Cohort Paper Analysis.do`
  - `NRH SCI Cohort Supplemental Sensitivity Analyses.do`

## Reproduction workflow
Run the canonical public smoke command from the repository root:

```bash
scripts/run_smoke.sh
# Windows PowerShell:
.\scripts\run_smoke.ps1
```

The launchers call `run_all.do smoke` and translate Stata's return code into the process status. Smoke accepts only `data/synthetic/Working NRH SCI.csv`. NRH-005 will add the approved fixture, so the command currently must fail clearly without falling back to restricted data. Set `STATA_BIN` for a different executable and `STATA_BATCH_FLAG` for a platform-specific batch flag.

The canonical argument order is `profile`, `data_dir`, `output_root`, and optional `run_id`. Use only generated or opaque non-sensitive run IDs. `full` and `release` require an explicit approved restricted-data directory containing `Working NRH SCI.csv`; neither may use or fall back to the synthetic fixture. `release` also refuses to overwrite existing generated `.dta` files:

```bash
stata-mp -b do run_all.do full "/approved/data" "/approved/output"
stata-mp -b do run_all.do release "/approved/data" "/approved/output" "path-safe-run-id"
```

On macOS, if Stata is not on `PATH`, configure the launcher:

```bash
STATA_BIN="/Applications/Stata/StataBE.app/Contents/MacOS/StataBE" \
STATA_BATCH_FLAG=-e scripts/run_smoke.sh
```

Direct execution remains available for legacy compatibility:

```bash
stata-mp -b do "NRH SCI Cohort Preprocessing.do"
stata-mp -b do "NRH SCI Cohort Paper Analysis.do"
stata-mp -b do "NRH SCI Cohort Supplemental Sensitivity Analyses.do"
```

## Expected outputs
Canonical runs write generated artifacts to one unique `Results and Figures/<run_id>/` directory. `run_all.log` and `run_manifest.csv` are value-free orchestration evidence; detailed child logs are under `Logs/`. Direct two-argument legacy calls retain `Results and Figures/<date>/`. These are outputs, not source files, and should not be committed unless a future release explicitly calls for public derived artifacts.

## Documentation surfaces to keep synchronized
- `README.md`: human-readable overview and run instructions.
- `llms.txt`: compact machine-readable summary for indexing.
- `CITATION.cff`: software and final article citation metadata.
- `data_dictionary.md` and `data_dictionary.csv`: expected inputs, derived variables, and output artifacts.

## Verification before publishing changes
- Run `git diff --check`.
- Validate `CITATION.cff` as YAML after citation edits.
- Search for hard-coded local absolute paths before publishing.
- Confirm `Data/`, `.dta`, Stata logs, `.gph`, and generated output folders remain ignored.
- Exercise invalid profiles, missing restricted inputs, wrong working directories, unsafe and duplicate run IDs, and unwritable output roots without using restricted data.
- Confirm `full` and `release` reject the public synthetic path and that failed runs never substitute another input.
- Confirm every failed run admitted past directory and writability setup has a value-free top-level log and manifest with a nonzero overall status.
- If Stata verification is run, inspect and then remove transient root-level `.log` files before committing.
