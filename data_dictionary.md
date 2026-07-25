# NRH-SCI-Vent Data Dictionary

This dictionary documents the expected local inputs, generated private datasets, key analysis variables, derived variables, and public output artifacts used by the NRH-SCI-Vent Stata workflow. It is derived from the public do-files and does not contain row-level data.

## Data Access

The study data are restricted retrospective EHR data governed by University of Utah IRB #00153003. Place private inputs in the ignored `Data/` directory. Do not commit source CSVs, generated `.dta` files, logs containing identifiers, or row-level exports.

## Files

| File | Role | Public? | Notes |
|---|---|---|---|
| `Data/Working NRH SCI.csv` | Restricted source CSV for preprocessing | No | Expected input to `NRH SCI Cohort Preprocessing.do`; all 36 ordered columns must conform to `validation/source_schema.csv`. |
| `Data/nrh-sci-raw.dta` | Generated raw Stata intermediate | No | Created from the source CSV by preprocessing. |
| `Data/nrh-sci-cleaned.dta` | Generated cleaned analysis dataset | No | Consumed by the paper and supplemental analysis scripts. |
| `Results and Figures/<run_id>/` | Canonical generated tables, figures, logs, and run evidence | No | Created by `run_all.do`; direct legacy calls retain `Results and Figures/<date>/`. |

## Key Variables

The CSV companion file has one row per documented file, source field, derived analysis variable, or generated output. The most important analytic constructs are:

| Variable | Definition | Values / units |
|---|---|---|
| `age` | Age at injury or admission as used in the analysis dataset. | Years; source-derived. |
| `age_decade` | Age scaled by decade for regression. | `age / 10`. |
| `level` | Rehab-assessed cervical injury level. | Encoded C1-C7 in the legacy code. |
| `high_vs_low` | Cervical level group used in paper models. | High = C4 or above; Low = C5 or below. |
| `comp_vs_part` | Motor completeness group. | Complete motor = AIS A/B; partial motor = AIS C/D. |
| `wean_during_day` | Daytime ventilator weaning achieved. | 0 = No; 1 = Yes. |
| `wean_24hr` | Liberation from invasive mechanical ventilation achieved. | 0 = No; 1 = Yes. |
| `decannulate` | Tracheostomy decannulation achieved. | 0 = No; 1 = Yes. |
| `weaning_outcome` | Ordinal respiratory outcome at rehab discharge. | 1 = fully ventilator dependent; 2 = daytime wean; 3 = liberated from IMV; 4 = decannulated. |
| `discharge_to` | Ordinal discharge destination. | 1 = LTAC; 2 = SNF; 3 = home with home health; 4 = home. |
| `time_to_censor_death` | Follow-up time from rehabilitation admission to death or censoring. | Days. |

## Review Flags

- The canonical ordered source-field, storage-type, missingness, sensitivity, and sanitized literal domains are defined in `validation/source_schema.csv`; `validation/validation_rules.csv` adds the identifier, date, nonnegative, placeholder, and mixed numeric/literal rules.
- Preprocessing imports the complete source CSV only after strict validation. Missing, extra, reordered, mistyped, or unexpected categorical fields fail rather than being silently truncated.
- Public categorical domains freeze normalized lexical inputs only. They do not approve the unresolved clinical mappings flagged here for NRH-004.
- Do not infer additional clinical definitions from private data without updating the public documentation.
