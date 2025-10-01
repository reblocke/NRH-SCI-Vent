# Contributing to NRH-SCI-Vent

Thanks for your interest in improving this research code! We welcome contributions that increase clarity and reproducibility. Because this work involves restricted clinical data, please read this carefully before opening a pull request.

## Ground rules

- **Never include data** in issues or PRs. Do not upload PHI/PII, raw EHR files, screenshots, or any derived tables that could identify individuals.
- Keep discussions and review **professional and respectful** (see our [Code of Conduct](CODE_OF_CONDUCT.md)).
- Prefer **small, focused PRs** with a clear problem statement and rationale.
- For changes that might affect figures or numbers in the manuscript, describe expected effects and attach regenerated outputs in a separate artifact file (not as tracked data).

## What contributions are helpful?

- Documentation: clarifying comments, variable definitions, or updating README tables.
- Portability: making paths and globals configurable; enabling non-interactive batch runs.
- Reproducibility: pinning Stata version-dependent behaviors; improving figure export steps.
- Lightweight tests: smoke tests that can run on a synthetic dataset to ensure scripts execute without error.

## Development setup

1. Install **Stata 18** and ensure `stata`, `stata-se`, or `stata-mp` is on your PATH.
2. Clone the repository and create a feature branch:
   ```sh
   git clone https://github.com/reblocke/NRH-SCI-Vent
   cd NRH-SCI-Vent
   git checkout -b feature/<short-name>
   ```
3. Run the scripts in **batch** to ensure they complete without interaction:
   ```sh
   stata-mp -b do "NRH SCI Cohort Preprocessing.do"
   stata-mp -b do "NRH SCI Cohort Paper Analysis.do"
   ```

## Style and structure

- Keep Stata code **modular**, with clear section headers and comments for each figure/table.
- Use **descriptive variable names** and document categorical encodings (e.g., `high_vs_low`, `comp_vs_part`).
- Export graphics to **TIFF** (for journals) and **PNG** where helpful for web previews, using the existing folder convention: `Results and Figures/<date>/`.

## Submitting a pull request

1. Open an issue describing the change and link any relevant manuscript section.
2. Submit your PR from your feature branch. Include:
   - A short summary of the change and **why** it’s needed.
   - Any impacts on outputs (attach images/logs to the PR conversation).
   - Confirmation that no PHI/PII or real patient-level data are included.
3. A maintainer will review for clarity and scope. We may request small edits before merge.

## Reporting problems or security concerns

Please **do not** include restricted data in a public issue. If you believe you’ve found a security/privacy risk or a problem that could cause leakage of PHI/PII, see our [Security Policy](SECURITY.md) for private reporting.
