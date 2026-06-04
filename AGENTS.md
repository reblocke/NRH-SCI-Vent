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
- The main scripts are:
  - `NRH SCI Cohort Preprocessing.do`
  - `NRH SCI Cohort Paper Analysis.do`
  - `NRH SCI Cohort Supplemental Sensitivity Analyses.do`

## Reproduction workflow
Run from the repository root with Stata 18 after the restricted cleaned dataset is available:

```bash
stata-mp -b do "NRH SCI Cohort Preprocessing.do"
stata-mp -b do "NRH SCI Cohort Paper Analysis.do"
stata-mp -b do "NRH SCI Cohort Supplemental Sensitivity Analyses.do"
```

The default private data directory is `Data/`. The expected source input is `Data/Working NRH SCI.csv`; preprocessing writes `Data/nrh-sci-raw.dta` and `Data/nrh-sci-cleaned.dta`. Optional script arguments are `data_dir` and `output_root`.

On macOS, if `stata-mp` is not on `PATH`, use the installed Stata app binary, for example:

```bash
/Applications/Stata/StataBE.app/Contents/MacOS/StataBE -b do "NRH SCI Cohort Supplemental Sensitivity Analyses.do"
```

## Expected outputs
Generated tables, figures, and logs are written under `Results and Figures/<date>/`. These are outputs, not source files, and should not be committed unless a future release explicitly calls for public derived artifacts.

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
- If Stata verification is run, inspect and then remove transient root-level `.log` files before committing.
