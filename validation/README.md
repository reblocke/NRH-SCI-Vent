# Public validation contracts

This directory contains the public, value-free contracts for the approved
NRH-000 behavioral baseline, the approved source schema, and the approved
NRH-004 clinical mappings. The private baseline has been populated, validated,
and approved by the scientific owner. Source-contract version 1 was approved by
the data steward on 2026-07-24. Source-contract version 2 and mapping-contract
version 1 were approved by the scientific owner and data steward on 2026-07-25.
Public approval metadata expose only roles, dates, and opaque secure-artifact
IDs; restricted values and paths remain outside Git.

## Public and private artifacts

| Artifact | Location | Git status | Contents |
|---|---|---|---|
| Project and approval metadata | `PROJECT.yml` | Public | Git anchors, contract versions, script hashes, and approval state |
| Expected-results contract | `validation/expected_results_schema.csv` | Public | Stable result IDs and comparison rules; every `expected_value` is blank |
| Sample-flow contract | `validation/sample_flow_schema.csv` | Public | Stable flow IDs and comparison rules; every expected count is blank |
| Source metadata | `validation/data_sources.csv` | Public | Value-free format, classification, version, and approval metadata |
| Ordered source schema | `validation/source_schema.csv` | Public | Sanitized field/type/domain contract without source values or counts |
| Source validation rules | `validation/validation_rules.csv` | Public | Aggregate-only rule definitions and strict/development actions |
| Clinical decision record | `DECISIONS.md` | Public | Eight approved decisions plus role/date/opaque-artifact metadata |
| Value-mapping contract | `validation/value_mappings.csv` | Public | Normalized accepted literals, missing actions, stable target codes and labels, decision IDs, and reject-unlisted policy |
| Baseline audit | `docs/BASELINE_AUDIT.md` | Public | Scope, limitations, and value-free signoff metadata |
| Approved baseline bundle | Institutionally approved storage outside Git | Private | Populated contracts, run record, logs, outputs, and permitted hashes |
| Source-contract review packet | Institutionally approved storage outside Git | Private | Raw-header/category review, aggregate counts, source checksum, validation log, and integrity manifest |
| NRH-004 decision packet | Institutionally approved storage outside Git | Private | Decision evidence, impact review, approval record, and integrity manifest |

`validation/private/` is ignored as defense in depth. It is not a substitute
for approved secure storage.

## Contract rules

The private `approved_expected_results.csv` is a copy of the public
expected-results contract with `expected_value` populated. The private
`approved_sample_flow.csv` is a copy of the public sample-flow contract with
the expected before, excluded, and after counts populated. Public templates
must never contain scientific values. Private copies must preserve every
public ID, comparison method, tolerance, evidence label, requirement flag, and
disclosure class unchanged.

Allowed comparison methods are:

- `exact`: integer counts, codes, labels, statuses, return codes, and stable
  SHA-256 values. Both tolerances must be zero.
- `numeric_tolerance`: full-precision numeric results pass when
  `abs(actual - expected) <= max(absolute_tolerance, relative_tolerance *
  abs(expected))`. Defaults are `1e-8` absolute and `1e-6` relative.
- `informational`: private provenance such as current TIFF or XLSX hashes that
  is recorded but does not determine pass or fail. These rows are not required.

For exact status rows, use `present` for an expected output. For `oparallel`,
use `success` or `controlled_failure` and record the matching integer return
code separately.

If only a rounded value can be recovered, set `recorded_precision` to the final
recorded unit, set relative tolerance to zero, and set absolute tolerance to
half that unit in a reviewed update to the public value-free schema. A private
bundle may not silently change precision, tolerances, methods, or other
contract fields. Do not label a rounded value as full precision.

Every sample-flow step must preserve the identity
`expected_before_n = expected_excluded_n + expected_after_n` when all three
fields apply. Model-flow `expected_after_n` values must equal the matching
`model_n` result. All result IDs and flow IDs must be unique and stable across
later refactors.

Populate all three count fields for every required private flow row. For
`import`, `save`, `load`, and `retain`, use `expected_excluded_n = 0` and set
before equal to after. For `exclude`, `subset`, and `model_complete_case`, use
the actual parent count as before, the removed or listwise-excluded count as
excluded, and the resulting cohort or `e(sample)` count as after.

The disclosure classes used in NRH-000 are:

- `public-metadata`: public provenance already recorded elsewhere in the
  repository; its expected-value cell remains blank in the template.
- `restricted-aggregate`: a value that remains only in the approved private
  bundle unless a later disclosure review explicitly authorizes publication.

## Secure baseline procedure

1. Use a separate clean checkout or worktree of
   `6ba1caf49ac9e260723d850b2a9f189f83255258`. Confirm that the three do-file
   hashes match `PROJECT.yml` and that the working tree is clean.
2. Work only in an approved environment. Keep the source data read-only and
   create both a fresh writable secure data-staging directory and a new, empty,
   run-specific output directory outside the Git repository. Place an approved
   copy or read-only link to the source CSV in staging; point `data_dir` there
   so the audited preprocessing script cannot overwrite historical derived
   files. Do not point it at the canonical research archive.
3. Record a private run identifier; start and end timestamps; Stata edition,
   version, update/build, and operating system; `adopath`; and the resolved
   locations and hashes for `table1_mc`, `coefplot`, `stackedcount`, graph
   schemes, and `oparallel`. Record any runtime installation or network use as
   a known limitation rather than changing the audited scripts.
4. Record an opaque source-extract identifier. Record the input SHA-256 only
   when institutional policy permits. Never copy the source path into public
   metadata.
5. Run the three audited scripts unchanged in separate clean Stata processes,
   preserving their documented order. Supply the approved data directory and
   the new private output root as arguments. Inspect each batch log and record
   its exit status before continuing. The private run record must retain the
   exact commit and clean/dirty status, command and arguments, per-script start
   and end times, exit status, and all known warnings.
6. After the unchanged workflow succeeds, use a reviewed private capture-only
   harness against the same staged cleaned dataset to collect results that the
   legacy scripts do not persist at full precision. The harness must reproduce
   the audited exclusions and model/test/margins commands exactly, extract
   aggregate `e()` and `r()` results immediately after each command, and emit
   no row-level content. Store its source, SHA-256, exact invocation, and review
   record in the private bundle. It supplements the audited run and may not
   replace it; any discrepancy with the workflow logs or retained outputs is a
   blocker.
7. Populate private copies named `approved_expected_results.csv` and
   `approved_sample_flow.csv`. Record an evidence pointer for every value, but
   do not place row-level excerpts or source paths in either file.
8. Store the populated contracts, `baseline_run_record.yml`, capture harness,
   logs, staging evidence, and outputs together under one opaque secure-artifact
   ID. Generated binary hashes are informational because the audited workflow
   is not yet deterministic. Retain or destroy identifier-bearing staged
   intermediates according to the approved data-governance policy.
9. Validate both private files against the public contracts. Move the public
   status to `pending_owner_approval` only after all required values and
   provenance fields are complete and internally consistent.
10. The scientific owner reviews the secure evidence. On approval, change the
   public status to `approved` and fill `approved_at`, `approver_role`, and the
   opaque `secure_artifact_id`. Do not publish the private values or location.

## Stata result capture conventions

- Model coefficients are stored on their native estimation scales:
  `log_odds`, `log_hazard`, or `log_subhazard`. Capture `_b[term]`,
  `_se[term]`, and confidence limits on that same scale. Do not paste the
  exponentiated OR, HR, or subhazard ratio displayed by commands using `or`.
- For a 95% interval on the stored scale, use
  `b +/- invnormal(.975) * SE`; use the matching two-sided model p-value.
- Ordered-logit cut points are captured separately and are not interpreted as
  clinical effects.
- In the Fine-Gray code, `high_vs_low` is entered as a numeric variable rather
  than `i.high_vs_low`. With the audited coding, its one-unit coefficient is
  the Low-versus-High contrast. Do not silently reinterpret it as a factor
  parameter during capture.
- Capture log-rank results from the secure `sts test` output or returned
  results. The graph's hard-coded annotation text is not baseline evidence.
- Capture supplemental margins from the full-precision margins dataset on the
  probability scale. The retained correspondence CSV is rounded to one decimal
  percentage point and must be marked precision-limited if used instead.
- Record `oparallel` status and return code even when the audited script treats
  a diagnostic failure as nonfatal. Record the reviewed conclusion as exactly
  `no_violation_detected`, `violation_detected`, or `not_computable`. Do not
  guess unreported diagnostic return values.
- The output-presence rows cover the final artifacts named in the NRH-000
  contract and use the exact value `present`. Intermediate graphs and logs are
  retained in the private bundle but are deferred to the complete output
  manifest in a later ticket.

## NRH-003 source contract

The public source contract is split across three value-free files:

| Contract | Purpose |
|---|---|
| `validation/data_sources.csv` | Public format, classification, version, and approval metadata for the restricted source |
| `validation/source_schema.csv` | Exact ordered 36-field names, storage and semantic types, missingness, sensitivity classes, and sanitized lexical domains |
| `validation/validation_rules.csv` | Stable aggregate-only checks for the identifier, dates, nonnegative timings, empty placeholder, and mixed numeric/literal fields |

`code/lib/data_contracts.do` imports the complete CSV and fails on missing,
extra, reordered, or mistyped columns. In strict mode it also fails on
required-value, duplicate-identifier, date, integer, nonnegative, placeholder,
mixed-value, and categorical-domain violations. The data-consuming `full` and
`release` profiles and standalone preprocessing use strict mode. The no-data
`verify` profile does not invoke the source validator. Development mode is
available only by directly invoking the validator: it may downgrade content
drift to a warning, but structural drift remains fatal.

The detailed source-contract log is private even though it contains only
aggregate evidence. Each row is limited to a rule ID, sanitized source-field
name, status, and exact invalid count. It must never include an input value,
identifier, row number, source path, checksum, or row count. The public
top-level orchestration log and manifest record only the contract version and
stage return code.

The data steward approved source-contract version 1 through the private packet
referenced by the opaque NRH-003 secure-artifact ID in `PROJECT.yml` and
`validation/data_sources.csv` on 2026-07-24. That approval covered the ordered
schema, lexical domains, sensitivity classes, missingness rules, and
aggregate-only disclosure format; it did not itself endorse clinical mappings.
NRH-004 separately approved the source-domain extensions and clinical mappings
described below.

## NRH-004 clinical mappings

`validation/value_mappings.csv`, version 1, is the executable public allowlist
for the approved categorical mappings. `code/lib/value_mappings.do` applies
Unicode NFC normalization, surrounding-whitespace trimming, and uppercasing
before matching. Each contract row identifies the source field, normalized
literal or blank input, target variable, action, code, label, decision ID, and
reject-unlisted policy. Preprocessing derives composites directly from the
mapped numeric codes with explicit missing-value guards.

The contract implements:

- age at inpatient rehabilitation admission;
- Male sex coded 0 = Female and 1 = Male, with female as the reference;
- the exact chest-tube negative, affirmative, missing, and reject actions;
- C1–C8 numeric coding with C1–C4 High and C5–C8 Low;
- the prospective C1–C8 by AIS A–D domain, including the approved AIA A and
  CENTRAL CORD normalizations;
- field-specific binary mappings with documented narrative exceptions;
- a prior-to-rehabilitation definition for `partial_wean_at_admit`; and
- the scientific-owner-controlled 2023-09-30 administrative censor date,
  exposed as `nrh_admin_censor_date_iso`.

Unexpected, conflicting, or otherwise unlisted categorical input fails rather
than being guessed. Approved blank input becomes missing; documented unknown
literals become missing only where explicitly allowlisted for the field. The
mapping helper may emit only safe aggregate diagnostics. Those diagnostics
remain private and may contain sanitized target-variable names, status, and
aggregate counts, but never identifiers, row numbers, source paths, or input
literals.

The scientific owner and data steward approved the eight decisions on
2026-07-25. Public metadata use the opaque secure-artifact ID
`NRH004-37fdea65-5a0b-414c-8a5b-d9356b044e61`; the private packet location,
contents, counts, and integrity hashes are not public.

## NRH-005 no-data public verification

No public synthetic cohort is authorized. NRH-005 therefore replaces the
deferred end-to-end smoke concept with the explicitly no-data `verify` profile.
Verify checks the controlled Stata dependency closure, loads both executable
contract libraries, validates the public contract files structurally, and
executes isolated approved-literal mapping unit cases. It does not accept a
data directory, open `Working NRH SCI.csv`, invoke preprocessing, or run either
analysis script.

`tests/test_public_contracts.do` operates only on public contract metadata and
does not construct patient- or cohort-shaped records. The atomic cases in
`tests/test_value_mappings.do` exercise individual public mapping literals
without combining clinical fields. Public execution coverage for row-level
source-content mutations is intentionally not claimed. Strict content
enforcement remains active when an authorized `full` or `release` run opens the
restricted source.

Successful verify evidence records `run_scope=no_data`, `data_accessed=false`,
the public-contract and mapping-unit stages as passed, and every data-consuming
stage as `not_run`. Evidence is written to a unique ignored
`Verification/<run_id>/` directory. The data-free GitHub workflow runs Python
contract and mutation validation plus launcher status tests; it does not claim
to run Stata or any analysis.

The scientific owner and data steward decided on 2026-07-25 not to authorize a
fabricated public cohort. Any future patient- or cohort-shaped fixture requires
a new explicit governance decision before implementation.

## Approval boundary

NRH-000 baseline approval, NRH-003 source-contract approval, and NRH-004
clinical-mapping approval are distinct. NRH-004 approves only the eight
decisions in `DECISIONS.md` and their versioned implementation contracts. It
does not authorize a release, tag, or merge through a change in a baseline
value or output meaning. Any such difference remains subject to the separate
value-free scientific-impact gate and release-line confirmation.
