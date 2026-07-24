# Authorized baseline validation

This directory contains the public, value-free contracts for NRH-000. It does
not contain an approved scientific baseline. The baseline remains
unapproved while it moves from `pending_authorized_run` to
`pending_owner_approval`; only the scientific owner's review can set it to
`approved`.

## Public and private artifacts

| Artifact | Location | Git status | Contents |
|---|---|---|---|
| Project and approval metadata | `PROJECT.yml` | Public | Git anchors, contract versions, script hashes, and approval state |
| Expected-results contract | `validation/expected_results_schema.csv` | Public | Stable result IDs and comparison rules; every `expected_value` is blank |
| Sample-flow contract | `validation/sample_flow_schema.csv` | Public | Stable flow IDs and comparison rules; every expected count is blank |
| Baseline audit | `docs/BASELINE_AUDIT.md` | Public | Scope, limitations, and value-free signoff metadata |
| Approved baseline bundle | Institutionally approved storage outside Git | Private | Populated contracts, run record, logs, outputs, and permitted hashes |

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

## Approval boundary

Creating these contracts does not approve the source schema, clinical coding,
model specification, censoring, disclosure safety, or current dependencies.
Any result-changing difference or unresolved scientific decision remains a
blocker for later tickets. NRH-000 is complete only when the private baseline
has been populated, validated, and approved.
