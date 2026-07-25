# Changelog

All notable public documentation and release metadata changes for this repository are recorded here.

## Unreleased

- Recorded the eight NRH-004 clinical-mapping and administrative-censoring decisions approved by the scientific owner and data steward on 2026-07-25, using public role/date/opaque-artifact metadata only.
- Advanced the strict source contract to version 2 for C8 and the prospective C1–C8 by AIS A–D domain.
- Added mapping-contract version 1 with deterministic source-to-code mappings, explicit missing actions, documented narrative exceptions, and reject-unlisted behavior.
- Replaced order-dependent categorical encoding in preprocessing with approved explicit mappings and direct numeric derivations while retaining stable target codes and labels.
- Named the scientific-owner-controlled 2023-09-30 administrative censoring constant and added structural and atomic mapping validation.
- Updated the public data dictionaries and workflow documentation. No restricted data, populated aggregates, or real-data outputs are included.
- Preserved the paper and supplemental analysis scripts, dependency closure, `CITATION.cff`, release tags, and latest-release metadata. No release or tag is created by these unreleased changes.
- Re-scoped NRH-005 from a deferred synthetic-cohort smoke test to explicitly no-data public verification because no fabricated public cohort is authorized.
- Added the `verify` profile, cross-platform return-code launchers, value-free no-data manifests, structural Stata contract checks, launcher regression tests, and data-free GitHub validation.
- Retired the `smoke` profile with a clear nonzero migration message. Full and release remain the only data-consuming profiles.

## v0.2.0 - paper and author-response code

Release date: 2026-05-30

- Updated repository documentation for the final article: Fenger C, Locke BW, Barker J, et al. *Post-Acute Ventilator Weaning and Discharge Outcomes in Individuals With Cervical Spinal Cord Injury*. Respiratory Care. 2026;71(4):351-357. https://doi.org/10.1177/19433654261428406
- Added documentation for the published author response: Fenger C, Locke B, Brown J. *Letter: Author Response to Letter to the Editor from Luo et al*. Respiratory Care. OnlineFirst. First published May 11, 2026. https://doi.org/10.1177/19433654261450544
- Documented `NRH SCI Cohort Supplemental Sensitivity Analyses.do`, including its Stata 18 requirement, restricted-data requirement, expected output directory, response figures, correspondence table, and run log.
- Updated `CITATION.cff` so the preferred citation points to the final Respiratory Care article and the software version points to GitHub release `v0.2.0`.
- No generated outputs or restricted source data are included in the repository.

## v0.1.0 - manuscript snapshot

Release date: 2026-04-04

- Frozen manuscript-associated analysis snapshot at commit `b437708ea58aba0f4da0afe73e54f64de423b7d2`.
- Preserved the paper analysis and preprocessing scripts used for the submitted manuscript snapshot.
