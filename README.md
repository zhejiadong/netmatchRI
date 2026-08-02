# netmatchRI

`netmatchRI` implements dual-penalty matching and randomization-based inference
for observational network data with network dependence and interference. The
dual-penalty matching method incorporates network proximity across treatment
arms and among units assigned to the same matched sets. The package provides
randomization-based inference that accounts for residual network dependence
across matched sets, sensitivity analysis of inferential conclusions across
dependence parameters, and critical curves that quantify the minimum network
dependence required to render observed significance a spurious association.
It also includes covariate-only and single-penalty comparison designs and
diagnostic summaries of the matching design.

## Installation

Install the package directly from GitHub:

```r
install.packages("remotes")
remotes::install_github("zhejiadong/netmatchRI")
library(netmatchRI)
```

The package installs its required R dependencies automatically.

## Solver recommendation

The default `solver = "highs"` uses the open-source HiGHS backend. We selected
HiGHS as the preferred open-source default based on package-local benchmarks
for the current dual-penalty formulation at N = 300--500; this is not a claim
that HiGHS is universally the fastest solver. HiGHS is a required dependency
and is installed automatically with `netmatchRI`.

For explicit licensed-solver preference, `solver = "auto"` tries Gurobi first,
then HiGHS, then the open-source GLPK fallback. You can request any backend with
`solver = "gurobi"`, `solver = "highs"`, or `solver = "glpk"`. Gurobi requires
its optimizer, license, and R package; GLPK uses `{Rglpk}` and `{slam}`.

The defaults are `timelimit = 90` seconds and `mipgap = 0.01`. In package-local
N = 500 runs, HiGHS may stop at the time limit with a feasible incumbent rather
than a proven optimum. `netmatch()` accepts such a result only after independent
bounds, integrality, constraint, and objective validation, and warns with the
actual status and available relative gap.

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
  solver = "highs"
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
  eta = seq(0, 0.03, by = 0.01),
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
- `RI_naive()`, `RI_decay()`, and `RI_design()` run randomization-based inference
  for matched designs, including sensitivity analysis and design-based approaches.
- `netmatch_sensitivity()` evaluates p-values over an `(eta, rho)` grid.
- `critical_sensitivity()` computes the critical sensitivity curve.
- `plot_sensitivity()` plots sensitivity results.
- `simulate_netmatch_example()` provides a 300-unit example dataset.

## Additional documentation

The package vignette `vignettes/netmatchRI.Rmd` gives a longer walkthrough of
simulation, dual-penalty matching, randomization-based inference, and
sensitivity analysis for observational network data.

## Citation

If you use this repository or the accompanying methods in your work, please
cite:

Dong Z, Lee Y (2026). "Design and Analysis for Valid Causal Inference with
Network-Dependent Data." Manuscript in preparation.
