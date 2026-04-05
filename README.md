# NRH-SCI-Vent

> Code to reproduce the analyses and figures for ventilator-weaning outcomes in adults with **cervical spinal cord injury (CSCI)** during inpatient rehabilitation at the University of Utah Craig H. Neilsen Rehabilitation Hospital.

**Links & IDs**
- Repository: https://github.com/reblocke/NRH-SCI-Vent
- Manuscript snapshot release: [`v0.1.0`](https://github.com/reblocke/NRH-SCI-Vent/releases/tag/v0.1.0)
- Latest analysis commit used for the manuscript: `b437708ea58aba0f4da0afe73e54f64de423b7d2` (2025‑09‑29)
- Manuscript: *in review* (no DOI yet). Working title: **Predictors of Ventilator Weaning and Discharge Outcomes in Cervical Spinal Cord Injury Subjects: A Retrospective Analysis**.
- Related abstract (CHEST 2023): *Just keep trying: Prior attempts at weaning do not determine eventual liberation from tracheostomy and mechanical ventilation in high-level spinal cord patients.*
- Statistical software: **Stata 18** (StataCorp, College Station, TX).

## Cite this work
Please cite the **software** until the article is accepted and has a DOI. The frozen manuscript snapshot is GitHub release [`v0.1.0`](https://github.com/reblocke/NRH-SCI-Vent/releases/tag/v0.1.0), which points to the exact analysis commit documented above. See [`CITATION.cff`](./CITATION.cff) for citation metadata you can paste into your reference manager. When the paper is published and an archive DOI has been minted for this release, prefer citing the article and that exact archived code release DOI.

## Quick start (reproduce the main results)

> **Requirements**: Stata 18 with graph export support for PNG/TIFF; ability to read/write to your working directory.

1) **Clone** this repository and move into it:
```sh
git clone https://github.com/reblocke/NRH-SCI-Vent
cd NRH-SCI-Vent
```

2) **Prepare data** (restricted). The analysis uses retrospective EHR data containing PHI/PII and cannot be distributed publicly. Access was under University of Utah IRB #00153003. If you have appropriate approvals and the analysis dataset, place the cleaned analysis dataset(s) in the location expected by the preprocessing script (see comments at the top of `NRH SCI Cohort Preprocessing.do`).

> If you do **not** have access to the restricted data, you can still review all code and outputs. We recommend creating or substituting a synthetic/de‑identified dataset with the same variable names and types to execute the pipeline end‑to‑end for demonstration purposes.

3) **Run the pipeline** from Stata (GUI or batch). In batch on macOS/Linux/Windows:
```sh
stata-mp -b do "NRH SCI Cohort Preprocessing.do"
stata-mp -b do "NRH SCI Cohort Paper Analysis.do"
# or use 'stata-se' / 'stata' depending on your license
```
On completion, figures and tables will be written under `Results and Figures/<date>/`.

### Expected outputs
The analysis script exports (filenames may be updated as wording evolves):
- `Fig 2 - stacked_states.tiff` — stacked daily counts of weaning & discharge status.
- `Fig 3 - Ordinal Regressions.tiff` — ORs for ventilator independence and discharge.
- `Supp Figure - CIFs for milestones.tiff` — cumulative incidence for day‑wean, liberation from IMV, decannulation.
- `Figure 4 - KMs for death.tiff` — KM curves by discharge location and weaning milestone.
- Supplemental tables in `xlsx` under `Results and Figures/<date>/`.

(Figure 1 was generated separately)

## Data access & ethics
- **Population**: Adults with CSCI admitted 2015–2022 to rehabilitation on continuous IMV via tracheostomy.
- **IRB**: University of Utah IRB **#00153003**; retrospective, exempt category.
- **Human‑subject protections**: Source EHR data contain PHI/PII and are **not** shared in this repository. Researchers wishing to reproduce results must obtain appropriate IRB approval and data use permissions at their institution.

## Repository layout
```
├── NRH SCI Cohort Preprocessing.do      # Prepares analysis dataset(s) and helper variables
├── NRH SCI Cohort Paper Analysis.do     # Produces tables/figures and model outputs
├── Results and Figures/                 # Created on first run; contains outputs by date
└── LICENSE                              # Code license (MIT)
```

## Workflow
1. `NRH SCI Cohort Preprocessing.do` — reads the source analysis dataset(s), constructs analysis variables, and writes intermediate files as needed.
2. `NRH SCI Cohort Paper Analysis.do` — fits proportional‑odds and Fine–Gray models, generates figures/tables, and exports publication graphics (TIFF).

## Paper ↔ code mapping
| Paper item | Script | Command/section | Output |
|---|---|---|---|
| Fig 2: Time-course of ventilator-weaning milestones | `NRH SCI Cohort Paper Analysis.do` | `stackedcount state day, ...` | `Fig 2 - stacked_states.tiff` |
| Fig 3: Predictors of milestones & discharge | `NRH SCI Cohort Paper Analysis.do` | `ologit ...` + export | `Fig 3 - Ordinal Regressions.tiff` |
| Supp CIFs: day-wean / liberation / decannulation | `NRH SCI Cohort Paper Analysis.do` | `stcurve, cif ...` + export | `Supp Figure - CIFs for milestones.tiff` |
| Fig 4: Mortality Kaplan–Meier | `NRH SCI Cohort Paper Analysis.do` | `sts graph ...` + export | `Figure 4 - KMs for death.tiff` |
| Flow diagram (Fig 1) | external diagram tool | (not generated by code) | (see manuscript) |

## Definitions used in this project
- **Weaning**: withdrawal of ventilator support over time.
- **Liberated from IMV**: no longer requiring invasive mechanical ventilation (some subjects may remain tracheostomized).
- **Decannulated**: tracheostomy removed.
- **Desirability convention**: For regression results, OR/SHR **> 1** indicate more independence / less resource‑intensive disposition / higher cumulative probability of achieving the milestone.

## Environment
- Tested with **Stata 18**. Analyses rely on base Stata commands and `stcurve`, `ologit`, and graphics export.
- Hardware: standard laptop/desktop is sufficient (no GPU required). Runs complete in minutes on a typical workstation.

## Funding & acknowledgements
This work was supported by **NIH** Ruth L. Kirschstein NRSA **5T32HL105321** (B.W.L), **American Thoracic Society** ASPIRE Fellowship (B.W.L), and the **Intermountain Foundation** (B.W.L). We thank the respiratory therapy, PM&R, and pulmonary teams at the Craig H. Neilsen Rehabilitation Hospital.

**Conflicts**: B.W.L. reports equity in Mountain Biometrics, Inc. (unrelated to this study). J.P.B. reports consulting fees from Breas Medical and Baxter. Others report no conflicts.

## License
Code is released under the **MIT License** (see `LICENSE`). Generated figures/tables are © the authors and journal policies may apply once published.

## Contributing, conduct & security
We welcome fixes to documentation, clarity of variable definitions, and portability improvements. Please **do not** submit PHI/PII or any data files via issues or PRs. See:
- [`CONTRIBUTING.md`](./CONTRIBUTING.md)
- [`CODE_OF_CONDUCT.md`](./CODE_OF_CONDUCT.md)
- [`SECURITY.md`](./SECURITY.md)

## Maintainers / contact
- Maintainer: Brian W. Locke (GitHub: `@reblocke`)
- Corresponding author: Jeanette P. Brown — jeanette.brown@hsc.utah.edu
