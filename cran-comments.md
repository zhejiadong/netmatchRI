This is a resubmission.

## Response to the previous CRAN review

The previous review asked us to declare where the optional, non-CRAN `gurobi`
package can be obtained. The package now:

- states in `DESCRIPTION` that eligible academic users can request a Gurobi
  license, links to the optimizer download page, and links to the official R
  interface installation instructions;
- gives the corresponding optional-backend setup steps in `README.md`;
- keeps the CRAN package `highs` as the required open-source default, with all
  Gurobi use conditional.

## Revision summary

The package provides network-constrained matching and randomization-based
inference for observational network data. It exposes critical eta at rho = 1
and the critical ratio through `critical_sensitivity()`. The design-based
covariance bound is evaluated by a finite sum rather than a quantile-grid
approximation. The rank statistic, ATT/ATC/ATE matching-weight options, and
matching/caliper behavior are unchanged.

## Test environment and completed checks

Windows 11 x64; R 4.5.0.

Exact artifact: `netmatchRI_0.1.0.tar.gz` (69,887 bytes).
SHA-256: `fdd7c9123496ad1c7612014dc5b50d0690d69f8dfd1478a461ed7dfda7e639f1`.

- 1,117 source test assertions across 72 test blocks passed; no failures,
  errors, warnings, or skips.
- Ordinary `R CMD check` on the exact source archive completed with 0 errors,
  0 warnings, and 0 notes.
- `R CMD check --as-cran` on the exact source archive completed with 0 errors,
  0 warnings, and 1 NOTE. Its installed test run had 1,112 passing assertions
  and 2 expected skips because optional Gurobi was unavailable.
- The exact archive built the PDF and HTML manuals, rebuilt the vignette, and
  ran all examples successfully.
- The exact archive installed in an isolated library; its installed citation,
  package metadata, and default HiGHS workflow were verified. The smoke test
  matched 60 units in 13 sets and reproduced critical eta 0.3702144390796152
  at rho = 1 and critical ratio 0.3468691726339535.
- All package URLs passed `urlchecker::url_check()`.
- GitHub Actions passed at commit
  `3debdda308ad8a3734f489473ab11cad5d90d679` on Windows R-release, macOS
  R-release, Ubuntu R-devel, Ubuntu R-release, and Ubuntu R-oldrel-1:
  <https://github.com/zhejiadong/netmatchRI/actions/runs/35813166652>.

## NOTE explanation

The sole NOTE reports a new submission and optional `gurobi` in `Enhances`.
The Description supplies the official academic-license, optimizer-download, and
R-interface installation links. Gurobi is not required to install, load, test,
or use the default HiGHS workflow.

## External checks / submission state

The exact 69,887-byte archive passed win-builder R-release and R-devel with
0 errors, 0 warnings, and the same expected NOTE described above. Installation,
examples, tests, vignettes, and PDF/HTML manuals all passed on both versions:

- R-release: <https://win-builder.r-project.org/Wr4092MLmkrN/00check.log>
- R-devel: <https://win-builder.r-project.org/3772CIG7aYDM/00check.log>

The package was uploaded through the CRAN submission form and the maintainer
confirmation step was completed on 2026-09-23. The submission is awaiting
CRAN's automated pretest/manual review. No release tag has been made.

This file is excluded from the source tarball; it records evidence for the
submission form without changing the checked payload.
