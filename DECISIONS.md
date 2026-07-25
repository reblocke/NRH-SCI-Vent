# NRH-SCI-Vent Decision Record

This public record documents approved scientific and data-governance decisions
that are implemented in version-controlled contracts. It contains no patient
data, cohort or model counts, estimates, private source paths, private artifact
hashes, or restricted run output.

## NRH-004 approval

| Field | Public metadata |
|---|---|
| Status | Approved |
| Approval date | 2026-07-25 |
| Approval roles | Scientific owner; data steward |
| Secure artifact ID | `NRH004-37fdea65-5a0b-414c-8a5b-d9356b044e61` |
| Mapping contract | `validation/value_mappings.csv`, version 1 |
| Source contract | `validation/source_schema.csv`, version 2 |

The approval applies to the following eight decisions. The executable mapping
rows use these stable decision IDs.

### NRH004-D001 — Age timepoint

`age` is age at inpatient rehabilitation admission.

### NRH004-D002 — Male sex indicator

The cleaned variable remains named `sex` and is a male indicator: `F` maps to
`0`, `M` and `MALE` map to `1`, and female is the reference category. Blank
input is missing; unlisted input is rejected. Public labels and documentation
describe this variable as Male sex.

### NRH004-D003 — Chest-tube indicator

`N` maps to `0`. `Y`, `Y. AFTER INTUBATION`, and
`Y. AFTER INTUBATION. PLEURAL EFFUSION` map to `1`. Blank and documented
unknown input are missing; unlisted input is rejected.

### NRH004-D004 — Cervical injury level

The accepted level domain is C1 through C8, coded `1` through `8`. C1–C4 form
the High group and C5–C8 form the Low group. The supplemental C7–C8 grouping
is unchanged.

### NRH004-D005 — AIS domain and motor grouping

The prospective compound domain is C1–C8 crossed with AIS A–D. A level-specific
`AIA A` alias normalizes to AIS A, and `CENTRAL CORD` normalizes to AIS D.
AIS A/B form the Complete Motor group; AIS C/D form the Partial Motor group.

### NRH004-D006 — Remaining binary fields

Each remaining binary source field has an explicit field-specific negative-to-0
and affirmative-to-1 mapping. Blank input is missing; documented `UNKNOWN` or
`UNKOWN` literals are also missing only where explicitly allowlisted for that
field. Unlisted or conflicting input is rejected. Documented narrative
exceptions remain explicit allowlisted rows rather than broad catch-all rules.

### NRH004-D007 — Partial weaning at rehabilitation admission

`partial_wean_at_admit` equals `1` only for explicit prior-to-rehabilitation or
prior-to-admission phrases. Known non-prior states map to `0`; blank and
documented unknown input are missing; unlisted input is rejected. In this
decision, admission means inpatient rehabilitation admission.

### NRH004-D008 — Administrative censor date

The scientific-owner-controlled administrative censor date remains
2023-09-30. Its public constant is
`nrh_admin_censor_date_iso`; preprocessing derives the corresponding
numeric daily date. A future change requires a new signed decision and a full
baseline review.

## Approval boundary

The public contracts expose approved normalized literals, missing actions,
target codes, labels, decision IDs, and reject-unlisted policy without exposing
restricted observations or aggregates. Approval of NRH-004 does not itself
authorize a release, tag, or merge when secure comparison identifies a change
in a baseline value or output meaning. Such a change requires the separate
scientific-impact gate.
