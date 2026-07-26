# NRH-SCI-Vent

[![DOI](https://img.shields.io/badge/DOI-10.1177%2F19433654261428406-blue)](https://doi.org/10.1177/19433654261428406)
[![PubMed](https://img.shields.io/badge/PMID-41823210-green)](https://pubmed.ncbi.nlm.nih.gov/41823210/)
[![Release](https://img.shields.io/github/v/release/reblocke/NRH-SCI-Vent?label=code%20release)](https://github.com/reblocke/NRH-SCI-Vent/releases/tag/v0.2.0)
[![Public validation](https://github.com/reblocke/NRH-SCI-Vent/actions/workflows/public-validation.yml/badge.svg?branch=main)](https://github.com/reblocke/NRH-SCI-Vent/actions/workflows/public-validation.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](./LICENSE)

> Code to reproduce the analyses and figures for ventilator-weaning outcomes in adults with **cervical spinal cord injury (CSCI)** during inpatient rehabilitation at the University of Utah Craig H. Neilsen Rehabilitation Hospital.

> **Development status:** `main` is the unreleased `0.3.0-dev` line.
> [`v0.2.0`](https://github.com/reblocke/NRH-SCI-Vent/releases/tag/v0.2.0)
> remains the latest citable code release; use a tagged release when exact
> reproducibility is required.

**Links & IDs**
- Repository: https://github.com/reblocke/NRH-SCI-Vent
- Paper and author-response code release: [`v0.2.0`](https://github.com/reblocke/NRH-SCI-Vent/releases/tag/v0.2.0)
- Original manuscript snapshot release: [`v0.1.0`](https://github.com/reblocke/NRH-SCI-Vent/releases/tag/v0.1.0)
- Latest analysis commit used for the final paper: `b437708ea58aba0f4da0afe73e54f64de423b7d2` (2025‑09‑29)
- Final article: Fenger C, Locke BW, Barker J, et al. **Post-Acute Ventilator Weaning and Discharge Outcomes in Individuals With Cervical Spinal Cord Injury.** *Respiratory Care*. 2026;71(4):351-357. https://doi.org/10.1177/19433654261428406
- PubMed: PMID [`41823210`](https://pubmed.ncbi.nlm.nih.gov/41823210/)
- Article license: open access under Creative Commons Attribution 4.0 on the SAGE article page.
- Author response: Fenger C, Locke B, Brown J. **Letter: Author Response to Letter to the Editor from Luo et al.** *Respiratory Care*. OnlineFirst. First published May 11, 2026. https://doi.org/10.1177/19433654261450544
- Related abstract (CHEST 2023): *Just keep trying: Prior attempts at weaning do not determine eventual liberation from tracheostomy and mechanical ventilation in high-level spinal cord patients.*
- Statistical software: **Stata 18** (StataCorp, College Station, TX).

## Cite this work
Please cite the final **Respiratory Care** article for the primary paper and cite this repository when referencing the analysis code. GitHub release [`v0.2.0`](https://github.com/reblocke/NRH-SCI-Vent/releases/tag/v0.2.0) contains the final paper code plus the supplemental analyses used for the author response. The earlier frozen manuscript snapshot remains available as [`v0.1.0`](https://github.com/reblocke/NRH-SCI-Vent/releases/tag/v0.1.0). See [`CITATION.cff`](./CITATION.cff) for citation metadata you can paste into your reference manager.

## Authors and affiliations
Article authors: Casey Fenger; Brian W. Locke ([ORCID](https://orcid.org/0000-0002-3588-5238)); James Barker; William Tang; Polly Creveling; Stormy Foster-Palmer; Alexandra Flis; Kevin Park; Jeffrey Rosenbluth; Jeanette P. Brown ([ORCID](https://orcid.org/0009-0009-8407-4034)). Primary affiliations include the University of Utah Craig H. Neilsen Rehabilitation Hospital, University of Utah Department of Physical Medicine and Rehabilitation, Spencer Fox Eccles School of Medicine, University of Utah Division of Pulmonary and Critical Care Medicine, and Intermountain Medical Center.

## Quick start

> **Requirements**: Stata 18 with graph export support for PNG/TIFF; ability
> to read/write to your working directory; and a local SHA-256 utility
> (`shasum` on macOS, `sha256sum` on Linux, or `certutil` on Windows).

1) **Clone** this repository and move into it:
```sh
git clone https://github.com/reblocke/NRH-SCI-Vent
cd NRH-SCI-Vent
```

2) **Choose a profile.** `run_all.do` is the canonical entry point and accepts, in order, `profile`, `data_dir`, `output_root`, and an optional `run_id`.

- `verify` is the default. It checks the controlled Stata dependencies, public contract structure, and atomic approved-literal mappings without accepting a data directory or running preprocessing or analyses.
- `full` requires an explicit approved restricted-data directory containing `Working NRH SCI.csv` and never infers or substitutes another input.
- `release` has the same restricted-input and filename requirement and also refuses to overwrite existing generated `.dta` files.
- `smoke` is retired because no public synthetic cohort is authorized. A future reconsideration requires an explicit governance decision.

Run explicitly no-data public verification from the repository root:

```sh
stata-mp -b do run_all.do verify
# macOS/Linux launcher
scripts/run_verify.sh
# Windows PowerShell launcher
.\scripts\run_verify.ps1
```

For an authorized restricted run, pass the input and output locations explicitly:

```sh
stata-mp -b do run_all.do full "/approved/data" "/approved/output"
stata-mp -b do run_all.do release "/approved/data" "/approved/output" "2026-07-24T120000_release"
```

The runner must start from the repository root. Every admitted invocation creates one unique, non-overwritable `<output_root>/<run_id>/` directory. The `verify` profile writes under `Verification/` by default and records `run_scope=no_data`, `data_accessed=false`, passing public-contract and mapping-unit stages, and every data-consuming stage as `not_run`. It never opens a source CSV or invokes preprocessing, paper analysis, or supplemental analysis. Full and release runs fail fast if dependency preflight, input preflight, strict source-contract validation, approved value mapping, or any analysis stage fails. For those restricted profiles, source-contract version 2 runs before preprocessing, and preprocessing executes mapping-contract version 1 from `validation/value_mappings.csv`.

Once directory and writability setup succeeds, the run contains a value-free `run_all.log` and `run_manifest.csv`; private source-validation diagnostics contain stable rule IDs, and private mapping diagnostics contain sanitized target-variable names, status, and safe aggregate counts—never input values, identifiers, row numbers, or source paths. Detailed child logs and generated outputs appear only for stages that execute. Use a generated run ID or another opaque, non-sensitive identifier; never encode a patient, source extract, or restricted location in `run_id`. Early invocation errors that occur before a safe run directory can be created—such as an invalid profile, wrong working directory, unsafe run ID, duplicate run ID, or unwritable output root—return a Stata error but do not create a manifest. The verify launchers translate Stata's return code into the process status, including on macOS Stata builds that otherwise return process status zero.

The POSIX launcher defaults to `-e` on macOS and `-b` on other Unix-like systems; the PowerShell launcher defaults to `/e` for Windows Stata. Set `STATA_BIN` for a different executable and `STATA_BATCH_FLAG` to override the platform default. For example, with the macOS app binary:

```sh
STATA_BIN="/Applications/Stata/StataBE.app/Contents/MacOS/StataBE" \
STATA_BATCH_FLAG=-e scripts/run_verify.sh
```

The legacy scripts remain directly executable for compatibility. Their first two optional arguments are `data_dir` and `output_root`; without an orchestration `run_id`, they retain the historical `Results and Figures/<date>/` layout:

```sh
stata-mp -b do "NRH SCI Cohort Preprocessing.do"
stata-mp -b do "NRH SCI Cohort Paper Analysis.do"
stata-mp -b do "NRH SCI Cohort Supplemental Sensitivity Analyses.do"
```

### Locked Stata dependencies

The project carries the exact seven-file runtime closure used by the authorized
baseline under `vendor/stata/plus/`. `code/00_preflight.do` puts that directory
first on Stata's ado-path, clears any cached user-written programs, verifies
every file against `vendor/stata/manifest.csv` with SHA-256, and confirms that
the required commands and graph schemes resolve to the controlled copies.
Canonical and standalone analysis runs perform this check before opening an
analysis dataset.

No analysis script installs packages or requires network access. Full and
release runs are therefore dependency-offline once the repository is cloned.
If the SHA-256 utility is unavailable, a file is missing or modified, or a
command resolves outside `vendor/stata/plus/`, preflight exits nonzero before
the source-data check. The dependency manifest hash and preflight status are
recorded in `run_manifest.csv`; detailed value-free evidence is written to
`Logs/dependency_preflight.log`.

The locked packages and citation sources are:

- Mark Chatfield's [`table1_mc`](https://ideas.repec.org/c/boc/bocode/s458351.html), version 3.3.
- Ben Jann's [`coefplot`](https://github.com/benjann/coefplot/tree/fe6c883881fbda70e506e6ae89a3921d7f220926), version 1.8.6 at the locked commit.
- Cameron Campbell's [`stackedcount`](https://ideas.repec.org/c/boc/bocode/s458825.html), version 1.0.
- Asjad Naqvi's [`schemepack`](https://github.com/asjadnaqvi/stata-schemepack/tree/962ca13accdc8a87aa3cfcda522ba81143c0a31f), version 1.4 at the locked commit.
- Maarten L. Buis's [`oparallel`](https://ideas.repec.org/c/boc/bocode/s457720.html), version 1.0.9, including its `gologit3_lf2` evaluator.

See [`vendor/stata/LICENSES.md`](./vendor/stata/LICENSES.md) for file-level
provenance and third-party license notices.

### Expected outputs
The analysis script exports (filenames may be updated as wording evolves):
- `Fig 2 - stacked_states.tiff` — stacked daily counts of weaning & discharge status.
- `Fig 3 - Ordinal Regressions.tiff` — ORs for ventilator independence and discharge.
- `Supp Figure - CIFs for milestones.tiff` — cumulative incidence for day‑wean, liberation from IMV, decannulation.
- `Figure 4 - KMs for death.tiff` — KM curves by discharge location and weaning milestone.
- Supplemental tables in `xlsx` under `Results and Figures/<run_id>/`.
- Orchestration evidence at `Results and Figures/<run_id>/run_all.log` and `run_manifest.csv`.
- Dependency-preflight and detailed script logs under `Results and Figures/<run_id>/Logs/`.
- No-data verification evidence under `Verification/<run_id>/`; it contains no analysis outputs.

(Figure 1 was generated separately)

### Author-response supplemental analyses
`NRH SCI Cohort Supplemental Sensitivity Analyses.do` documents the supplemental analyses and figures generated for the published author response. It recreates the paper-matched analytic cohort from the restricted cleaned dataset, writes response-focused figures and a compact correspondence table, and records the run log.

Run it from Stata 18 after the cleaned analysis dataset is available:
```sh
stata-mp -b do "NRH SCI Cohort Supplemental Sensitivity Analyses.do"
```

When run directly, expected response outputs are written under `Results and Figures/<date>/`; orchestration writes them under `<output_root>/<run_id>/`:
- `Supplemental Figure - Median Days to Milestones by Exact Injury Level.tiff`
- `Supplemental Figure - Milestone Rates by Finer Injury Groups.tiff`
- `Supplemental Figure - Observed Discharge Disposition by Age, Split by Decannulation.tiff`
- `Supplemental Figure - Discharge by Age Among Decannulated Patients.tiff`
- `Correspondence Table - Adjusted Home Discharge Probabilities by Milestone and Age.csv`
- `Logs/supplemental_sensitivity_analyses.log`

The response analyses require the same restricted dataset as the paper analyses. Source data and generated outputs are not tracked in this public repository.

## Data access & ethics
- **Population**: Adults with CSCI admitted 2015–2022 to rehabilitation on continuous IMV via tracheostomy.
- **IRB**: University of Utah IRB **#00153003**; retrospective, exempt category.
- **Human‑subject protections**: Source EHR data contain PHI/PII and are **not** shared in this repository. Researchers wishing to reproduce results must obtain appropriate IRB approval and data use permissions at their institution.
- **Data dictionary**: Public documentation of expected inputs, derived variables, and output artifacts is available in [`data_dictionary.md`](./data_dictionary.md) and [`data_dictionary.csv`](./data_dictionary.csv). The exact ordered 36-field source contract is in [`validation/source_schema.csv`](./validation/source_schema.csv), with additional checks in [`validation/validation_rules.csv`](./validation/validation_rules.csv).
- **Approved mappings**: [`DECISIONS.md`](./DECISIONS.md) records the eight NRH-004 decisions, and [`validation/value_mappings.csv`](./validation/value_mappings.csv) is the executable mapping allowlist. Public approval metadata contain the approval roles, date, and opaque secure-artifact ID; these files contain no row-level data, restricted counts, private paths, or real-data outputs.

## Repository layout
```
├── README.md                             # Human-readable project overview and run instructions
├── llms.txt                              # Machine-readable repository summary for LLM indexing
├── PROJECT.yml                           # Development, release, and authorized-baseline metadata
├── DECISIONS.md                          # Public scientific and data-governance decision record
├── CITATION.cff                          # Software and preferred article citation metadata
├── AGENTS.md                             # Repository-specific guidance for coding agents
├── run_all.do                            # Canonical verify/full/release orchestration entry point
├── code/00_preflight.do                  # Offline dependency resolution and SHA-256 gate
├── code/lib/data_contracts.do            # Strict source-schema and domain validator
├── code/lib/value_mappings.do            # Approved deterministic categorical mapping helper
├── config/                               # Profile defaults and overrides
├── scripts/                              # No-data verify launchers and retired smoke adapters
├── vendor/stata/                         # Locked runtime ado files, manifest, and licenses
├── validation/                           # Public value-free baseline, source, and mapping contracts
├── tests/test_public_contracts.do        # Data-free Stata contract-structure checks
├── tests/test_value_mappings.do          # Atomic accepted/missing/rejected mapping tests
├── tests/validate_no_data_workflow.py    # Public no-data boundary and workflow checks
├── data_dictionary.md / data_dictionary.csv
│                                          # Public variable, file, and output documentation
├── NRH SCI Cohort Preprocessing.do      # Prepares analysis dataset(s) and helper variables
├── NRH SCI Cohort Paper Analysis.do     # Produces tables/figures and model outputs
├── NRH SCI Cohort Supplemental Sensitivity Analyses.do
│                                          # Produces author-response supplemental outputs
├── Data/                                # Ignored private inputs and generated .dta files
├── Results and Figures/                 # Ignored generated outputs by run ID or legacy date
├── Verification/                        # Ignored value-free no-data verification evidence
└── LICENSE                              # Code license (MIT)
```

## Workflow
1. `run_all.do` — runs explicitly no-data public verification or, for `full` and `release`, validates paths and unique run IDs, runs dependency and input preflight, enforces strict source-contract version 2, invokes each legacy stage, checks cleaned-data and output presence, and finalizes value-free run evidence.
2. `code/00_preflight.do` — prepends the locked ado tree, verifies every manifest SHA-256, and confirms controlled command and scheme resolution without network access.
3. `NRH SCI Cohort Preprocessing.do` — independently re-enforces the strict source contract for standalone safety, reads the complete `Working NRH SCI.csv` without positional truncation, applies mapping-contract version 1 through `code/lib/value_mappings.do`, derives composites directly from stable numeric codes with missing-value guards, and writes `nrh-sci-raw.dta` and `nrh-sci-cleaned.dta` beside the private input.
4. `NRH SCI Cohort Paper Analysis.do` — reads `nrh-sci-cleaned.dta`, fits proportional‑odds and Fine–Gray models, generates figures/tables, and exports publication graphics (TIFF).
5. `NRH SCI Cohort Supplemental Sensitivity Analyses.do` — reads `nrh-sci-cleaned.dta`, generates the supplemental figures, correspondence table, and model log used to support the published author response.

## Machine-readable metadata
- [`llms.txt`](./llms.txt) summarizes the repository purpose, article identifiers, run order, data restrictions, and agent cautions.
- [`CITATION.cff`](./CITATION.cff) provides software citation metadata and the preferred citation for the final paper.
- [`data_dictionary.csv`](./data_dictionary.csv) is the machine-usable companion to the Markdown data dictionary.
- [`DECISIONS.md`](./DECISIONS.md) and [`validation/value_mappings.csv`](./validation/value_mappings.csv) provide the public NRH-004 decision record and executable normalized-literal mapping contract.

## Paper ↔ code mapping
Release [`v0.1.0`](https://github.com/reblocke/NRH-SCI-Vent/releases/tag/v0.1.0) is the original manuscript snapshot. Release [`v0.2.0`](https://github.com/reblocke/NRH-SCI-Vent/releases/tag/v0.2.0) contains the final paper code and author-response supplemental-analysis code.

| Paper item | Script | Command/section | Output |
|---|---|---|---|
| Fig 2: Time-course of ventilator-weaning milestones | `NRH SCI Cohort Paper Analysis.do` | `stackedcount state day, ...` | `Fig 2 - stacked_states.tiff` |
| Fig 3: Predictors of milestones & discharge | `NRH SCI Cohort Paper Analysis.do` | `ologit ...` + export | `Fig 3 - Ordinal Regressions.tiff` |
| Supp CIFs: day-wean / liberation / decannulation | `NRH SCI Cohort Paper Analysis.do` | `stcurve, cif ...` + export | `Supp Figure - CIFs for milestones.tiff` |
| Fig 4: Mortality Kaplan–Meier | `NRH SCI Cohort Paper Analysis.do` | `sts graph ...` + export | `Figure 4 - KMs for death.tiff` |
| Flow diagram (Fig 1) | external diagram tool | (not generated by code) | (see manuscript) |
| Author response: milestone timing by exact injury level | `NRH SCI Cohort Supplemental Sensitivity Analyses.do` | exact-level milestone timing figure block | `Supplemental Figure - Median Days to Milestones by Exact Injury Level.tiff` |
| Author response: milestone rates by finer injury groups | `NRH SCI Cohort Supplemental Sensitivity Analyses.do` | grouped injury-level milestone figure block | `Supplemental Figure - Milestone Rates by Finer Injury Groups.tiff` |
| Author response: age, decannulation, and discharge | `NRH SCI Cohort Supplemental Sensitivity Analyses.do` | observed age/discharge figure blocks | `Supplemental Figure - Observed Discharge Disposition by Age, Split by Decannulation.tiff`; `Supplemental Figure - Discharge by Age Among Decannulated Patients.tiff` |
| Author response: adjusted home-discharge probabilities | `NRH SCI Cohort Supplemental Sensitivity Analyses.do` | milestone-adjusted discharge model and margins block | `Correspondence Table - Adjusted Home Discharge Probabilities by Milestone and Age.csv` |

## Definitions used in this project
- **Weaning**: withdrawal of ventilator support over time.
- **Liberated from IMV**: no longer requiring invasive mechanical ventilation (some subjects may remain tracheostomized).
- **Decannulated**: tracheostomy removed.
- **Desirability convention**: For regression results, OR/SHR **> 1** indicate more independence / less resource‑intensive disposition / higher cumulative probability of achieving the milestone.

## Environment
- Tested with **Stata 18**. The workflow uses built-in Stata survival, regression, and graphics commands plus the exact project-controlled dependencies recorded in `vendor/stata/manifest.csv`.
- Dependency verification is offline and fail-closed. It uses the operating system's local SHA-256 utility but never downloads or updates a package at runtime.
- The launchers and orchestration paths are cross-platform. The current supplemental TIFF fallback still uses macOS command-line tools (`sips`, `grep`, and `touch`), so full end-to-end portability is not yet claimed.
- Hardware: standard laptop/desktop is sufficient (no GPU required). Runs complete in minutes on a typical workstation.

## Funding & acknowledgements
This work was supported by **NIH** Ruth L. Kirschstein NRSA **5T32HL105321** (B.W.L), **American Thoracic Society** ASPIRE Fellowship (B.W.L), and the **Intermountain Foundation** (B.W.L). We thank the respiratory therapy, PM&R, and pulmonary teams at the Craig H. Neilsen Rehabilitation Hospital.

**Conflicts**: B.W.L. reports equity in Mountain Biometrics, Inc. (unrelated to this study). J.P.B. reports consulting fees from Breas Medical and Baxter. Others report no conflicts.

## License
NRH-SCI-Vent project code is released under the **MIT License** (see
`LICENSE`). Third-party files under `vendor/stata/plus/` retain their
GPL-3.0-only or MIT licenses as documented in
[`vendor/stata/LICENSES.md`](./vendor/stata/LICENSES.md). The final SAGE article
is open access under **CC BY 4.0** on the article page. Restricted source data,
private generated datasets, local logs, and generated output folders are not
included in this repository.

## Contributing, conduct & security
We welcome fixes to documentation, clarity of variable definitions, and portability improvements. Please **do not** submit PHI/PII or any data files via issues or PRs. See:
- [`CONTRIBUTING.md`](./CONTRIBUTING.md)
- [`CODE_OF_CONDUCT.md`](./CODE_OF_CONDUCT.md)
- [`SECURITY.md`](./SECURITY.md)

## Maintainers / contact
- Maintainer: Brian W. Locke (GitHub: `@reblocke`)
- Corresponding author: Jeanette P. Brown — jeanette.brown@hsc.utah.edu
