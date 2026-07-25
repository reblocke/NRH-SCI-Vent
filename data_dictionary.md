# NRH-SCI-Vent Data Dictionary

This dictionary documents the expected local inputs, generated private datasets, key analysis variables, derived variables, and public output artifacts used by the NRH-SCI-Vent Stata workflow. It is derived from the public do-files and does not contain row-level data.

## Data Access

The study data are restricted retrospective EHR data governed by University of Utah IRB #00153003. Place private inputs in the ignored `Data/` directory. Do not commit source CSVs, generated `.dta` files, logs containing identifiers, or row-level exports.

## Files

| File | Role | Public? | Notes |
|---|---|---|---|
| `Data/Working NRH SCI.csv` | Restricted source CSV for preprocessing | No | Expected input to `NRH SCI Cohort Preprocessing.do`; all 36 ordered columns must conform to source-contract version 2 in `validation/source_schema.csv`. |
| `Data/nrh-sci-raw.dta` | Generated raw Stata intermediate | No | Created from the source CSV by preprocessing. |
| `Data/nrh-sci-cleaned.dta` | Generated cleaned analysis dataset | No | Consumed by the paper and supplemental analysis scripts. |
| `Results and Figures/<run_id>/` | Canonical generated tables, figures, logs, and run evidence | No | Created by `run_all.do`; direct legacy calls retain `Results and Figures/<date>/`. |
| `Verification/<run_id>/` | Value-free no-data verification log and manifest | No | Created by `run_all.do verify`; records `data_accessed=false` and all data-consuming stages as `not_run`. |

## Key Variables

The CSV companion file has one row per documented file, source field, derived analysis variable, or generated output. The most important analytic constructs are:

| Variable | Definition | Values / units |
|---|---|---|
| `age` | Age at inpatient rehabilitation admission. | Years; source-derived. |
| `age_decade` | Age scaled by decade for regression. | `age / 10`. |
| `sex` | Male sex indicator with female as the reference category. | 0 = Female; 1 = Male. |
| `level` | Rehabilitation-assessed cervical injury level. | C1–C8 coded 1–8. |
| `asia_class` | AIS classification at rehabilitation admission. | 1 = AIS A; 2 = AIS B; 3 = AIS C; 4 = AIS D. |
| `high_vs_low` | Cervical level group used in paper models. | High = C1–C4; Low = C5–C8. |
| `comp_vs_part` | Motor completeness group. | Complete motor = AIS A/B; partial motor = AIS C/D. |
| `partial_wean_at_admit` | Partial ventilator wean before inpatient rehabilitation admission. | 1 only for explicit prior-to-rehabilitation or prior-to-admission phrases; 0 for known non-prior states. |
| `wean_during_day` | Daytime ventilator weaning achieved. | 0 = No; 1 = Yes. |
| `wean_24hr` | Liberation from invasive mechanical ventilation achieved. | 0 = No; 1 = Yes. |
| `decannulate` | Tracheostomy decannulation achieved. | 0 = No; 1 = Yes. |
| `weaning_outcome` | Ordinal respiratory outcome at rehab discharge. | 1 = fully ventilator dependent; 2 = daytime wean; 3 = liberated from IMV; 4 = decannulated. |
| `discharge_to` | Ordinal discharge destination. | 1 = LTAC; 2 = SNF; 3 = home with home health; 4 = home. |
| `time_to_censor_death` | Follow-up time from rehabilitation admission to death or administrative censoring. | Days; non-deaths are censored at the date derived from `nrh_admin_censor_date_iso` = 2023-09-30. |

## Review Flags

- The canonical ordered source-field, storage-type, missingness, sensitivity, and sanitized literal domains are defined by source-contract version 2 in `validation/source_schema.csv`; `validation/validation_rules.csv` adds the identifier, date, nonnegative, placeholder, and mixed numeric/literal rules.
- Preprocessing imports the complete source CSV only after strict validation. Missing, extra, reordered, mistyped, or unexpected categorical fields fail rather than being silently truncated.
- The approved executable mapping contract is `validation/value_mappings.csv`, version 1. It uses normalized allowlisted literals, explicit missing actions, stable target codes and labels, and reject-unlisted policy.
- `DECISIONS.md` records the eight NRH-004 decisions approved on 2026-07-25 by the scientific owner and data steward. Public approval metadata uses only the roles, date, and opaque secure-artifact ID.
- The administrative censor date remains controlled by NRH004-D008. A future change requires a new signed decision and full baseline review.
- Do not infer additional clinical definitions or accepted literals from private data without a new reviewed public-contract update.
