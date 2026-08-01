## Test environments

- local Windows 11 x64, R 4.5.0

## R CMD check results

0 errors | 0 warnings | 1 note

- This is a new submission.
- `gurobi` is listed in `Enhances` but is not in a mainstream R repository. It is an optional interface to the separately distributed commercial Gurobi Optimizer and is used only after `requireNamespace("gurobi", quietly = TRUE)` and a license probe succeed. The package defaults to the open-source CRAN package `highs`; Gurobi is not required to install, load, test, or use `netmatchRI`.
