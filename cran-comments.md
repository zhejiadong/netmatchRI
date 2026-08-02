## Resubmission

This is a resubmission in response to CRAN's request to declare where the
optional `gurobi` package can be obtained.

- The `Description` field now provides links for eligible academic users to
  request a free Gurobi license, to download Gurobi Optimizer, and to install
  the non-CRAN `gurobi` R package distributed with the optimizer.
- The README now provides step-by-step optional-backend installation guidance.
- Gurobi remains entirely optional. The package uses the open-source CRAN
  package `highs` as its default solver, and Gurobi is not required to install,
  load, test, or use the default functionality of `netmatchRI`.

## Test environments

- local Windows 11 x64, R 4.5.0
- GitHub Actions, Windows, R-release
- GitHub Actions, macOS, R-release
- GitHub Actions, Ubuntu, R-devel
- GitHub Actions, Ubuntu, R-release
- GitHub Actions, Ubuntu, R-oldrel-1
- Win-builder, R-release
- Win-builder, R-devel

## R CMD check results

0 errors | 0 warnings | 1 note

- This is a new submission.
- `gurobi` is listed in `Enhances` but is not in a mainstream R repository. The
  `Description` field now states where users can request a license, download
  Gurobi Optimizer, and install its bundled R package. The optional interface is
  used only after `requireNamespace("gurobi", quietly = TRUE)` and a license
  probe succeed. The package defaults to the open-source CRAN package `highs`;
  Gurobi is not required to install, load, test, or use `netmatchRI`.

## Additional checks

- Package URLs were checked successfully.
- The source tarball installs and the default HiGHS workflow passes a clean-library smoke test.
- Win-builder checks passed on R-release and R-devel with the same expected NOTE. PDF and HTML manuals were generated successfully.
