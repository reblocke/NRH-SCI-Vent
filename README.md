# NRH-SCI-Vent

[![DOI](https://img.shields.io/badge/DOI-10.1177%2F19433654261428406-blue)](https://doi.org/10.1177/19433654261428406)
[![PubMed](https://img.shields.io/badge/PMID-41823210-green)](https://pubmed.ncbi.nlm.nih.gov/41823210/)
[![Release](https://img.shields.io/github/v/release/reblocke/NRH-SCI-Vent?label=code%20release)](https://github.com/reblocke/NRH-SCI-Vent/releases/tag/v0.2.0)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](./LICENSE)

> Code to reproduce the analyses and figures for ventilator-weaning outcomes in adults with **cervical spinal cord injury (CSCI)** during inpatient rehabilitation at the University of Utah Craig H. Neilsen Rehabilitation Hospital.

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

## Quick start (reproduce the main results)

> **Requirements**: Stata 18 with graph export support for PNG/TIFF; ability to read/write to your working directory.

1) **Clone** this repository and move into it:
```sh
git clone https://github.com/reblocke/NRH-SCI-Vent
cd NRH-SCI-Vent
```

2) **Prepare data** (restricted). The analysis uses retrospective EHR data containing PHI/PII and cannot be distributed publicly. Access was under University of Utah IRB #00153003. If you have appropriate approvals, place the restricted source CSV at `Data/Working NRH SCI.csv`. The `Data/` directory is ignored by git. The preprocessing script writes private generated Stata datasets to `Data/nrh-sci-raw.dta` and `Data/nrh-sci-cleaned.dta`.

> If you do **not** have access to the restricted data, you can still review all code and outputs. We recommend creating or substituting a synthetic/de‑identified dataset with the same variable names and types to execute the pipeline end‑to‑end for demonstration purposes.

3) **Run the pipeline** from Stata (GUI or batch). In batch on macOS/Linux/Windows:
```sh
stata-mp -b do "NRH SCI Cohort Preprocessing.do"
stata-mp -b do "NRH SCI Cohort Paper Analysis.do"
stata-mp -b do "NRH SCI Cohort Supplemental Sensitivity Analyses.do"
# or use 'stata-se' / 'stata' depending on your license
```
On completion, figures and tables will be written under `Results and Figures/<date>/`.

Optional script arguments are `data_dir` and `output_root`. For example:
```sh
stata-mp -b do "NRH SCI Cohort Paper Analysis.do" Data "Results and Figures"
```

### Expected outputs
The analysis script exports (filenames may be updated as wording evolves):
- `Fig 2 - stacked_states.tiff` — stacked daily counts of weaning & discharge status.
- `Fig 3 - Ordinal Regressions.tiff` — ORs for ventilator independence and discharge.
- `Supp Figure - CIFs for milestones.tiff` — cumulative incidence for day‑wean, liberation from IMV, decannulation.
- `Figure 4 - KMs for death.tiff` — KM curves by discharge location and weaning milestone.
- Supplemental tables in `xlsx` under `Results and Figures/<date>/`.
- Logs under `Results and Figures/<date>/Logs/`.

(Figure 1 was generated separately)

### Author-response supplemental analyses
`NRH SCI Cohort Supplemental Sensitivity Analyses.do` documents the supplemental analyses and figures generated for the published author response. It recreates the paper-matched analytic cohort from the restricted cleaned dataset, writes response-focused figures and a compact correspondence table, and records the run log.

Run it from Stata 18 after the cleaned analysis dataset is available:
```sh
stata-mp -b do "NRH SCI Cohort Supplemental Sensitivity Analyses.do"
```

Expected response outputs are written under `Results and Figures/<date>/`:
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
- **Data dictionary**: Public documentation of expected inputs, derived variables, and output artifacts is available in [`data_dictionary.md`](./data_dictionary.md) and [`data_dictionary.csv`](./data_dictionary.csv). These files do not contain row-level data.

## Repository layout
```
├── README.md                             # Human-readable project overview and run instructions
├── llms.txt                              # Machine-readable repository summary for LLM indexing
├── CITATION.cff                          # Software and preferred article citation metadata
├── AGENTS.md                             # Repository-specific guidance for coding agents
├── data_dictionary.md / data_dictionary.csv
│                                          # Public variable, file, and output documentation
├── NRH SCI Cohort Preprocessing.do      # Prepares analysis dataset(s) and helper variables
├── NRH SCI Cohort Paper Analysis.do     # Produces tables/figures and model outputs
├── NRH SCI Cohort Supplemental Sensitivity Analyses.do
│                                          # Produces author-response supplemental outputs
├── Data/                                # Ignored private inputs and generated .dta files
├── Results and Figures/                 # Ignored generated outputs by date
└── LICENSE                              # Code license (MIT)
```

## Workflow
1. `NRH SCI Cohort Preprocessing.do` — reads `Data/Working NRH SCI.csv`, constructs analysis variables, and writes `Data/nrh-sci-raw.dta` and `Data/nrh-sci-cleaned.dta`.
2. `NRH SCI Cohort Paper Analysis.do` — reads `Data/nrh-sci-cleaned.dta`, fits proportional‑odds and Fine–Gray models, generates figures/tables, and exports publication graphics (TIFF).
3. `NRH SCI Cohort Supplemental Sensitivity Analyses.do` — reads `Data/nrh-sci-cleaned.dta`, generates the supplemental figures, correspondence table, and model log used to support the published author response.

## Machine-readable metadata
- [`llms.txt`](./llms.txt) summarizes the repository purpose, article identifiers, run order, data restrictions, and agent cautions.
- [`CITATION.cff`](./CITATION.cff) provides software citation metadata and the preferred citation for the final paper.
- [`data_dictionary.csv`](./data_dictionary.csv) is the machine-usable companion to the Markdown data dictionary.

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
- Tested with **Stata 18**. The workflow uses built-in Stata survival, regression, and graphics commands, plus existing user-written commands used by the legacy analysis workflow, including `table1_mc`, `coefplot`, and `stackedcount`.
- Hardware: standard laptop/desktop is sufficient (no GPU required). Runs complete in minutes on a typical workstation.

## Funding & acknowledgements
This work was supported by **NIH** Ruth L. Kirschstein NRSA **5T32HL105321** (B.W.L), **American Thoracic Society** ASPIRE Fellowship (B.W.L), and the **Intermountain Foundation** (B.W.L). We thank the respiratory therapy, PM&R, and pulmonary teams at the Craig H. Neilsen Rehabilitation Hospital.

**Conflicts**: B.W.L. reports equity in Mountain Biometrics, Inc. (unrelated to this study). J.P.B. reports consulting fees from Breas Medical and Baxter. Others report no conflicts.

## License
Code is released under the **MIT License** (see `LICENSE`). The final SAGE article is open access under **CC BY 4.0** on the article page. Restricted source data, private generated datasets, local logs, and generated output folders are not included in this repository.

## Contributing, conduct & security
We welcome fixes to documentation, clarity of variable definitions, and portability improvements. Please **do not** submit PHI/PII or any data files via issues or PRs. See:
- [`CONTRIBUTING.md`](./CONTRIBUTING.md)
- [`CODE_OF_CONDUCT.md`](./CODE_OF_CONDUCT.md)
- [`SECURITY.md`](./SECURITY.md)

## Maintainers / contact
- Maintainer: Brian W. Locke (GitHub: `@reblocke`)
- Corresponding author: Jeanette P. Brown — jeanette.brown@hsc.utah.edu
