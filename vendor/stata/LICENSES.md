# Vendored Stata dependency notices

The repository's root MIT license applies to the NRH-SCI-Vent project code. It
does not replace the licenses of the third-party files under
`vendor/stata/plus/`.

The vendored files are the smallest runtime closure used by the authorized
NRH-000 baseline. They are unmodified copies whose SHA-256 hashes are recorded
in `manifest.csv`.

| Package | Runtime files | Version/source anchor | License |
|---|---|---|---|
| `table1_mc` by Mark Chatfield | `t/table1_mc.ado` | SSC/RePEc S458351; approved copy identifies version 3.3 dated 2022-05-05 | GPL-3.0-only |
| `coefplot` by Ben Jann | `c/coefplot.ado` | Git commit `fe6c883881fbda70e506e6ae89a3921d7f220926`; version 1.8.6 | MIT |
| `stackedcount` by Cameron Campbell | `s/stackedcount.ado` | SSC/RePEc S458825; help version 1.0 dated 2020-07-28 | GPL-3.0-only |
| `schemepack` by Asjad Naqvi | `s/scheme-white_w3d.scheme`; `s/scheme-white_tableau.scheme` | Git commit `962ca13accdc8a87aa3cfcda522ba81143c0a31f`; version 1.4 | MIT |
| `oparallel` by Maarten L. Buis | `o/oparallel.ado`; `g/gologit3_lf2.ado` | SSC/RePEc S457720; version 1.0.9 revised 2019-06-13 | GPL-3.0-only |

License texts are provided in:

- `licenses/GPL-3.0-only.txt`
- `licenses/MIT-coefplot.txt`
- `licenses/MIT-schemepack.txt`

Authoritative source and citation pages:

- <https://ideas.repec.org/c/boc/bocode/s458351.html>
- <https://github.com/benjann/coefplot/tree/fe6c883881fbda70e506e6ae89a3921d7f220926>
- <https://ideas.repec.org/c/boc/bocode/s458825.html>
- <https://github.com/asjadnaqvi/stata-schemepack/tree/962ca13accdc8a87aa3cfcda522ba81143c0a31f>
- <https://ideas.repec.org/c/boc/bocode/s457720.html>

Optional documentation, examples, and unused export helpers are not vendored.
In particular, `oparallel_ex.ado` is not part of the runtime closure and is not
executed by this project.
