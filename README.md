# netmatchRI

`netmatchRI` provides network-constrained matching and randomization-based
inference for observational network data. Its three matching methods are
dual-penalty, single-penalty, and covariate-only matching.

## Installation

Install the source archive supplied with the release:

```r
install.packages("netmatchRI_0.1.0.tar.gz", repos = NULL, type = "source")
library(netmatchRI)
```

## Basic workflow

The input data should contain a binary treatment $Z$, baseline covariates $X$,
an outcome $Y$, and an adjacency or network-distance matrix.

```r
sim <- simulate_netmatch_example(seed = 90141, n = 60, beta_z = 2)

m <- netmatch(
  data = sim$data,
  treat = "Z",
  covariates = c("X1", "X2", "X3"),
  network = sim$net_dist,
  method = "dual",
  kappa = 2,
  estimand = "ATT"
)

summary(m)
diag <- diagnose_match(m)
diag$covariate_balance
diag$network_summary

m_covariate <- netmatch(
  sim$data, "Z", c("X1", "X2", "X3"), sim$net_dist,
  method = "covariate"
)
m_single <- netmatch(
  sim$data, "Z", c("X1", "X2", "X3"), sim$net_dist,
  method = "single", kappa = 2
)
designs <- list(
  "Dual-penalty" = m,
  "Covariate-only" = m_covariate,
  "Single-penalty" = m_single
)
plot_covariate_balance(designs)
plot_covariate_similarity(designs)

ri <- RI_adjusted(m, "Y", eta = 0.03, rho = 0.10)
ri$result

sens <- sensitivity_grid(
  m, "Y",
  eta = seq(0, 0.5, by = 0.1),
  rho = seq(0, 0.3, by = 0.05)
)
sens$grid
plot_sensitivity(sens)

critical <- critical_sensitivity(m, "Y")
critical
critical$summary
plot_sensitivity(critical)
```

`critical$summary` reports the critical eta at `rho = 1` and the critical
ratio. See `?critical_sensitivity` for their definitions.

`estimand = "ATT"`, `"ATC"`, or `"ATE"` selects unit matching weights for
matched-data summaries and balance diagnostics. It does not change the
sharp-null RI test, whose matched-set weights are chosen by `weight_type`.
The defaults are `estimand = "ATT"` and `weight_type = "ns"`.

The default dual-penalty matching solver is the open-source HiGHS backend. Set
`solver = "auto"` to try HiGHS, Gurobi, and GLPK in that order, or select a
backend explicitly. Gurobi requires a separate optimizer, license, and R package;
see the [official installation guide](https://docs.gurobi.com/projects/optimizer/en/current/reference/r/setup.html).

## Main functions

- `netmatch()` constructs dual-penalty, single-penalty, or covariate-only matched designs.
- `matched_data()` and `summary()` inspect a matched design; `diagnose_match()` checks covariate balance and network distance diagnostics.
- `plot_covariate_balance()` and `plot_covariate_similarity()` compare any selected matching designs, optionally including the unmatched sample.
- `RI_unadjusted()`, `RI_adjusted()`, and `RI_design()` perform randomization-based inference.
- `sensitivity_grid()` and `critical_sensitivity()` evaluate sensitivity to residual network dependence; the critical result includes critical eta at `rho = 1` and the critical ratio.
- `plot_sensitivity()` plots sensitivity grids and critical curves.
- `simulate_netmatch_example()` generates example network data.

## Citation and support

Use `citation("netmatchRI")` for the installed citation. To report a problem or
request a feature, open an issue at <https://github.com/zhejiadong/netmatchRI/issues>.
