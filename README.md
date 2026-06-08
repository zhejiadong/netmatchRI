# netmatchRI

`netmatchRI` implements network-penalty matching method and randomization
inference for observational studies when multiple variables exihibt network dependence. Its main design is
the dual-penalty matching formulation, which excludes units close in the network to be placed in the same matched set (both within and across treatment arms).
It also provide the matched-based randomization inference methods to adjust for residual network dependence across matched sets. It provides sensitivity anlaysis for evaluating how robusteness the inferential conclusions sensitive to the strength of network dependece and the critical curves for quantifing the minimal network dependence required for rendering the observed significance as a result of spurious association.
The package also includes covariate-only and single-penalty comparison designs,
diagnostic summaries of matching design.

## Installation

Install the package directly from GitHub:

```r
install.packages("remotes")
remotes::install_github("zhejiadong/netmatchRI")
library(netmatchRI)
```

The package installs its required R dependencies automatically.

## Solver recommendation

We recommend using Gurobi whenever it is available, especially for larger or
more demanding matching problems. In `netmatchRI`, the simplest way to do that
is to set `solver = "auto"`: the package will use Gurobi when it is installed
and licensed, and otherwise fall back to the open-source GLPK backend.

If you want to use Gurobi explicitly, first install the Gurobi Optimizer and
activate a valid Gurobi license on your machine, then install the `{gurobi}` R
package following the official Gurobi instructions for your platform. After
that, you can call `netmatch(..., solver = "gurobi")`. If Gurobi is not
available, `solver = "auto"` remains the recommended default for portability.

## Basic workflow

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

diagnostics <- diagnose_match(m_dual)
diagnostics$covariate_balance
diagnostics$network_distance

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

For your own study, supply one row per unit, a binary treatment indicator such
as `Z`, observed covariates such as `X1`, `X2`, `X3`, an outcome column such as
`Y`, and a network distance or adjacency representation compatible with
`netmatch()`.

## Main functions

- `netmatch()` builds a matched design under the selected matching method.
- `diagnose_match()` summarizes covariate balance and within-set network
  distances.
- `RI_naive()`, `RI_decay()`, and `RI_design()` run randomization inference for
  matched designs.
- `netmatch_sensitivity()` evaluates p-values over an `(eta, rho)` grid.
- `critical_sensitivity()` computes the critical sensitivity curve.
- `plot_sensitivity()` plots sensitivity results.
- `simulate_netmatch_example()` provides a 300-unit example dataset.

## Additional documentation

The package vignette `vignettes/netmatchRI.Rmd` gives a longer walkthrough of
simulation, matching, inference, and sensitivity analysis.

## Citation

If you use this repository or the accompanying methods in your work, please
cite:

Dong, Z., & Lee, Y. "Design and Analysis for Valid Causal Inference with
Network-Dependent Data." In process.
