## Revision summary

The package provides network-constrained matching and randomization-based
inference for observational network data. It exposes critical eta at rho = 1
and the critical ratio through `critical_sensitivity()`. The design-based
covariance bound is evaluated by a finite sum rather than a quantile-grid
approximation. The rank statistic, ATT/ATC/ATE matching-weight options, and
matching/caliper behavior are unchanged.

## Test environment and completed checks

Windows 11 x64; R 4.5.0.

Exact artifact: `netmatchRI_0.1.0.tar.gz` (70,359 bytes).
SHA-256: `349f2e1294fdad5d2825da585ed9314d42df91879450244f862c2309fd87c048`.

- 1,117 source test assertions across 72 test blocks passed; no test failures,
  errors, warnings, or skips.
- `R CMD check --as-cran` on the exact source archive completed with 0 errors,
  0 warnings, and 1 NOTE.
- The exact archive built the PDF and HTML manuals, rebuilt the vignette, and
  ran all examples and tests successfully.
- The exact archive installed in an isolated library; its installed citation,
  package metadata, and default HiGHS workflow were verified.

## NOTE explanation

The sole NOTE reports a new submission and optional `gurobi` in `Enhances`.
The Description supplies the official academic-license, optimizer-download, and
R-interface installation links. Gurobi is not required to install, load, test,
or use the default HiGHS workflow.

## External checks / submission state

No CRAN upload, Git push, release tag, or current-artifact remote-platform
check has been made. The public GitHub repository should be aligned with this
release candidate before submission because the package Description links to it.

This file is excluded from the source tarball; it records evidence for the
submission form without changing the checked payload.
