# netmatchRI

`netmatchRI` implements the dual-penalty matching framework for observational
studies with network dependence. The proposed design excludes close
treated-control pairs and prevents close same-arm units from appearing in the
same matched set. The package also provides randomization inference and
sensitivity analysis for the resulting matched design.

Covariate-only and single-penalty designs are included as comparison methods.

## Installation

```r
install.packages("remotes")
remotes::install_local("path/to/netmatchRI")
library(netmatchRI)
```

Installing `netmatchRI` also installs the core R dependencies used by the
package. Gurobi is optional. When Gurobi is unavailable, `solver = "auto"`
falls back to the open-source GLPK backend.

## Example

```r
sim <- simulate_netmatch_example()

m_dual <- netmatch(
  data = sim$data,
  treat = "Z",
  covariates = c("X1", "X2", "X3"),
  network = sim$net_dist,
  method = "dual",
  kappa = 2,
  solver = "auto"
)

diagnose_match(m_dual)

ri_naive <- RI_naive(m_dual, "Y")
ri_decay <- RI_decay(m_dual, "Y", eta = 0.03, rho = 0.10)
ri_design <- RI_design(m_dual, "Y")

sens <- netmatch_sensitivity(
  m_dual,
  "Y",
  eta = c(0, 0.01, 0.03),
  rho = seq(0, 1, by = 0.1)
)

crit <- critical_sensitivity(m_dual, "Y")

plot_sensitivity(sens, type = "pvalue")
plot_sensitivity(crit, type = "critical")
```

If you use your own data, provide one row per unit, a binary treatment column
such as `Z`, covariate columns such as `X1`, `X2`, `X3`, and an outcome column
such as `Y` for inference.

## Main Functions

- `netmatch()` builds the dual-penalty matched design via MIP.
- `diagnose_match()` summarizes covariate balance and within-set network
  distances.
- `RI_naive()`, `RI_decay()`, and `RI_design()` run randomization inference.
- `netmatch_sensitivity()` evaluates p-values over an `(eta, rho)` grid.
- `critical_sensitivity()` computes the critical sensitivity curve.
- `plot_sensitivity()` plots p-value sensitivity curves and critical curves.
- `simulate_netmatch_example()` provides a built-in 300-unit example.
