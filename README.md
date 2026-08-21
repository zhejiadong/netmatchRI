# netmatchRI

`netmatchRI` provides matching and randomization-based inference for observational
network data. It supports dual-penalty matching, comparison designs, balance and
network diagnostics, inference under residual network dependence, and
sensitivity analysis.

## Installation

```r
install.packages("remotes")
remotes::install_github("zhejiadong/netmatchRI")
library(netmatchRI)
```

## Basic workflow

```r
sim <- simulate_netmatch_example(seed = 20260821, n = 32)

m <- netmatch(
  data = sim$data,
  treat = "Z",
  covariates = c("X1", "X2", "X3"),
  network = sim$net_dist,
  method = "dual",
  kappa = 2
)

summary(m)
diag <- diagnose_match(m)
diag$covariate_balance
diag$network_summary

ri <- RI_Sensitivity(m, "Y", eta = 0.03, rho = 0.10)
ri$result

grid <- sensitivity_grid(
  m, "Y",
  eta = c(0, 0.03),
  rho = seq(0, 1, by = 0.25)
)
grid$grid
plot_sensitivity(grid)
```

The default dual-design solver is the open-source HiGHS backend. Set
`solver = "auto"` to try Gurobi, HiGHS, and GLPK in that order, or select a
backend explicitly. Gurobi requires a separate optimizer, license, and R package;
see the [official installation guide](https://docs.gurobi.com/projects/optimizer/en/current/reference/r/setup.html).

## Main functions

- `netmatch()` constructs dual-penalty, single-penalty, or covariate-only matched designs.
- `matched_data()`, `summary()`, and `diagnose_match()` inspect a matched design.
- `RI_Naive()`, `RI_Sensitivity()`, and `RI_Design()` perform randomization-based inference.
- `sensitivity_grid()` and `critical_sensitivity()` evaluate sensitivity to residual network dependence.
- `plot_sensitivity()` plots sensitivity grids and critical curves.
- `simulate_netmatch_example()` generates example network data.

## Citation and support

Use `citation("netmatchRI")` for the installed citation. To report a problem or
request a feature, open an issue at <https://github.com/zhejiadong/netmatchRI/issues>.
