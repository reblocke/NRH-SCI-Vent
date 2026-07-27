# NRH-000 authorized baseline audit

## Status

**Baseline status:** `approved`

An authorized secure run has populated and validated the private contracts,
and the scientific owner has approved the secure bundle as the authorized
behavioral comparison baseline. NRH-000's human gate is complete.

| Public signoff field | Value |
|---|---|
| Approved at | 2026-07-24 |
| Approver role | scientific owner |
| Secure artifact ID | `NRH000-b3941de5-7816-495e-961d-6665031606d5` |

No scientific count, estimate, source path, log excerpt, or generated
real-data output belongs in this document.

## Frozen Git anchors

| Role | Version or ref | Commit |
|---|---|---|
| Ref observed at audit time | `origin/main` | `6ba1caf49ac9e260723d850b2a9f189f83255258` |
| Immutable NRH-000 baseline anchor | commit | `6ba1caf49ac9e260723d850b2a9f189f83255258` |
| Release observed at audit time | `v0.2.0` | `45a50136bf6b0279b2a1cf34a51c559109596c0a` |

The immutable NRH-000 comparison anchor is the recorded commit, not the moving
`origin/main` ref. At audit time, the development line was identified as
`0.3.0-dev`, and `CITATION.cff` described the then-latest frozen release,
`v0.2.0`. NRH-000 did not create or modify a release.

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

At frozen baseline commit
`6ba1caf49ac9e260723d850b2a9f189f83255258`, the audited workflow had the
following limitations. This historical list does not describe current `main`
and does not endorse the underlying behavior:

- no public synthetic end-to-end execution path and no automated Stata tests
  at that commit;
- unpinned user-written ado dependencies and runtime installation behavior;
- positional source import without an exact public schema contract;
- unresolved clinical and categorical mappings documented as needing review;
- creation of an identifier-bearing raw intermediate before identifier removal;
- unseeded random jitter in supplemental figures;
- macOS/POSIX-specific supplemental image post-processing;
- hard-coded Kaplan-Meier p-value annotation text;
- date-based output directories that can overwrite same-day results;
- incomplete output determinism and claim-level traceability; and
- a sentinel baseline that does not yet structure every table cell or plotted
  aggregate as an independent claim-regression result.

Later reviewed changes have resolved some of these historical limitations while
preserving the immutable comparison anchor. Remaining limitations must not be
silently changed inside NRH-000 because doing so would invalidate the
pre-refactor baseline.

## Completion gate

The authorized analyst ran the unchanged audited code in the approved
environment and populated and validated the private contracts. The scientific
owner reviewed and approved the secure bundle on 2026-07-24, completing
NRH-000. Approval is represented publicly only by the status, approval date,
approver role, and opaque secure-artifact ID in `PROJECT.yml` and this
document.
