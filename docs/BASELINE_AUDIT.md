# NRH-000 authorized baseline audit

## Status

**Baseline status:** `pending_authorized_run`

The public contracts are present, but no approved scientific baseline is
claimed. The NRH-000 pull request must remain draft until the human gate in
`validation/README.md` is complete.

| Public signoff field | Value |
|---|---|
| Approved at | Pending |
| Approver role | Pending |
| Secure artifact ID | Pending |

No scientific count, estimate, source path, log excerpt, or generated
real-data output belongs in this document.

## Frozen Git anchors

| Role | Version or ref | Commit |
|---|---|---|
| Authoritative development baseline | `origin/main` | `6ba1caf49ac9e260723d850b2a9f189f83255258` |
| Latest frozen release | `v0.2.0` | `45a50136bf6b0279b2a1cf34a51c559109596c0a` |

The development baseline is identified as `0.3.0-dev`. `CITATION.cff`
continues to describe the latest frozen release, `v0.2.0`; NRH-000 does not
create or modify a release.

## Audited analysis files

| File | SHA-256 at the authoritative commit |
|---|---|
| `NRH SCI Cohort Preprocessing.do` | `d220ae6ab659760d7def73010ee985a14e6e2a65ad00c9c9ab309d81a0001a42` |
| `NRH SCI Cohort Paper Analysis.do` | `7e57f1bf30b8029a67d5a24bf601df4ffab7e0a839298c28d64b53179280e866` |
| `NRH SCI Cohort Supplemental Sensitivity Analyses.do` | `73e6b20e0f4e337286ea0b5787644c7e1c40b9452b8a5d5778201fe6a70b4480` |

These files are unchanged by NRH-000. Later equivalence checks compare against
the authorized private baseline rather than treating documentation or graph
annotations as scientific evidence.

## Baseline scope

The private baseline contract covers:

- preprocessing and analytic sample flow;
- analytic and clinically relevant category counts;
- model-specific complete-case counts;
- the Fine-Gray, ordered-logit, and Cox models executed by the paper script;
- the ordered-logit models and adjusted margins executed by the supplemental
  response script;
- competing-event, failure, censoring, log-rank, joint-test, and diagnostic
  results; and
- expected output presence plus non-gating binary hashes for provenance.

Individual patient observations, plotted patient-level points, dates, free
text, source paths, and suppression-sensitive values are excluded from the
public contract.

## Known limitations frozen for comparison

The baseline records current behavior; it does not endorse or resolve these
limitations:

- no public synthetic end-to-end execution path or automated Stata tests;
- unpinned user-written ado dependencies and runtime installation behavior;
- positional source import without an exact public schema contract;
- unresolved clinical and categorical mappings documented as needing review;
- creation of an identifier-bearing raw intermediate before identifier removal;
- unseeded random jitter in supplemental figures;
- macOS/POSIX-specific supplemental image post-processing;
- hard-coded Kaplan-Meier p-value annotation text;
- date-based output directories that can overwrite same-day results; and
- incomplete output determinism and claim-level traceability.
- a sentinel baseline that does not yet structure every table cell or plotted
  aggregate as an independent claim-regression result.

These limitations are work for later tickets. They must not be silently fixed
inside NRH-000 because doing so would invalidate the pre-refactor baseline.

## Completion gate

NRH-000 becomes complete only after an authorized analyst runs the unchanged
audited code in the approved environment, populates and validates the private
contracts, and the scientific owner approves the secure bundle. Approval is
represented publicly only by the status, approval date, approver role, and an
opaque secure-artifact ID in `PROJECT.yml` and this document.
