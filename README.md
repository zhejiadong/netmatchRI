# netmatchRI

`netmatchRI` provides network-constrained matching and randomization inference
for observational studies where units are connected by a graph. The main design
is the dual-penalty matched design: close treated-control pairs are excluded,
and close same-arm units are also prevented from appearing in the same matched
set.

The exported API is focused on reusable public workflows: simulation,
network-constrained matching, diagnostics, randomization inference, and
sensitivity analysis.

## What The Package Does

- `simulate_netmatch_example()` creates a 300-unit strong-dependence example
  data set with an outcome, treatment, covariates, adjacency matrix, graph
  distances, and covariance matrix.
- `netmatch()` builds a matched design. The core method is
  `method = "dual"`. `method = "covariate"` and `method = "single"` are
  comparison designs.
- `diagnose_match()` checks covariate balance and within-set network distances.
- `RI_naive()`, `RI_decay()`, and `RI_design()` run randomization inference.
- `netmatch_sensitivity()` evaluates p-values over an eta-rho sensitivity grid.
- `critical_sensitivity()` computes the critical eta curve for a dual-penalty
  matched design.
- `plot_sensitivity()` plots either p-value sensitivity curves or the critical
  eta curve.

## Installation

Install from the package source directory:

```r
install.packages("remotes")
remotes::install_local("path/to/netmatchRI")
library(netmatchRI)
```

If you are working inside this repository, build and install from the parent
directory:

```r
setwd("/Users/zhejia/Brown Dropbox/Zhejia Dong/24 NetworkMatching/MatchingNetDepedence/MatchingNet")
system("R CMD build netmatchRI")
install.packages("netmatchRI_0.0.1.tar.gz", repos = NULL, type = "source")
library(netmatchRI)
```

The dual-penalty method is a mixed-integer program. Use Gurobi when available:
it is the recommended solver for the 300-unit example and repeated analyses.
Without Gurobi, `solver = "auto"` falls back to the open-source GLPK backend
through `Rglpk` and `slam`. GLPK is useful as a backup for moderate problems,
but it can be slower for difficult or infeasible dual-penalty MIPs.

## Inspect The Package

R packages do not open a window or menu after installation. Inspect the package
from the R console:

```r
library(netmatchRI)

ls("package:netmatchRI")
help(package = "netmatchRI")
?netmatchRI
?netmatch
```

If RStudio reports that the package help index is missing, restart R and
reinstall the built source tarball:

```r
.rs.restartR()
install.packages("netmatchRI_0.0.1.tar.gz", repos = NULL, type = "source")
library(netmatchRI)

file.exists(system.file("html", "00Index.html", package = "netmatchRI"))
help(package = "netmatchRI")
```

## Dual-Penalty Workflow

Start with the built-in example:

```r
sim <- simulate_netmatch_example()

dat <- sim$data
net_dist <- sim$net_dist

table(dat$Z)
summary(lm(Y ~ Z, data = dat))$coefficients
```

If you want to use your own data, `dat` should be a data frame with one row per
unit, a binary treatment column such as `Z`, covariate columns such as `X1`,
`X2`, `X3`, and an outcome column such as `Y` for inference.

Build the proposed dual-penalty design:

```r
m_dual <- netmatch(
  data = dat,
  treat = "Z",
  covariates = c("X1", "X2", "X3"),
  network = net_dist,
  method = "dual",
  kappa = 2,
  solver = "gurobi"
)

m_dual
summary(m_dual)
```

If Gurobi is not available, try:

```r
m_dual <- netmatch(
  data = dat,
  treat = "Z",
  covariates = c("X1", "X2", "X3"),
  network = net_dist,
  method = "dual",
  kappa = 2,
  solver = "glpk",
  timelimit = 120
)
```

Check diagnostics:

```r
diag <- diagnose_match(m_dual)
diag$covariate_balance
diag$network_distance
diag$within_distance_table

plot(m_dual)
```

Run randomization inference:

```r
ri_naive <- RI_naive(m_dual, outcome = "Y")
ri_design <- RI_design(m_dual, outcome = "Y")
ri_decay <- RI_decay(
  m_dual,
  outcome = "Y",
  eta = 0.03,
  rho = 0.10
)

rbind(ri_naive$result, ri_decay$result, ri_design$result)
```

Evaluate p-value sensitivity:

```r
sens <- netmatch_sensitivity(
  match = m_dual,
  outcome = "Y",
  eta = c(0, 0.01, 0.02, 0.03, 0.05, 0.10),
  rho = seq(0, 1, by = 0.05)
)

head(sens$grid)
plot_sensitivity(sens, type = "pvalue")
```

Compute and plot the critical sensitivity curve:

```r
crit <- critical_sensitivity(
  match = m_dual,
  outcome = "Y",
  rho = seq(0, 1, by = 0.01),
  alpha = 0.05
)

head(crit$curve)
crit$interpretation

plot_sensitivity(crit, type = "critical")
plot_sensitivity(crit, type = "critical", critical_ylim = c(0, 1))
```

## Comparison Designs

The package also includes two comparison methods:

```r
m_cov <- netmatch(
  dat, "Z", c("X1", "X2", "X3"), net_dist,
  method = "covariate",
  kappa = 2
)

m_single <- netmatch(
  dat, "Z", c("X1", "X2", "X3"), net_dist,
  method = "single",
  kappa = 2
)
```

Use these to compare against the proposed dual-penalty design, not as the main
workflow.

## Function Reference

| Function | Main Returned Values |
| --- | --- |
| `simulate_netmatch_example()` | `data`, `Adj`, `net_dist`, `V`, `seed`, `dep_index`, `beta_z` |
| `netmatch()` | matched `data`, `match_table`, `network_distance`, `method`, `kappa`, `solver` |
| `summary()` | matched-set size summaries |
| `plot()` | within-set network-distance plot |
| `diagnose_match()` | `covariate_balance`, `network_distance`, `within_distance_table`, aliases `covariate_smd`, `average_within_distance` |
| `RI_naive()`, `RI_decay()`, `RI_design()` | `result`, `detail`, `covariance`, `set_distance` |
| `netmatch_sensitivity()` | p-value `grid`, source `match`, `outcome`, `kappa` |
| `critical_sensitivity()` | critical `curve`, `naive`, `interpretation` |
| `plot_sensitivity()` | a `ggplot` object |

## Development Check

From a terminal:

```sh
R CMD build netmatchRI
R CMD check netmatchRI_0.0.1.tar.gz --no-manual
```

A successful check ends with `Status: OK`.
