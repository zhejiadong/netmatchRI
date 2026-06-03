# netmatchRI

`netmatchRI` provides a small user-facing workflow for network-constrained
matching and randomization inference.

The package is designed for two uses:

- applied analysis of observational data with network dependence;
- reproduction of the simulation studies in *Design and Analysis for Valid
  Causal Inference with Network-Dependent Data*.

## Install

```r
# after the GitHub repository is published
remotes::install_github("<user>/netmatchRI")
```

## Basic workflow

```r
library(netmatchRI)

A <- matrix(0, 8, 8)
A[cbind(1:7, 2:8)] <- 1
A[cbind(2:8, 1:7)] <- 1

dat <- data.frame(
  Z = c(1, 1, 1, 1, 0, 0, 0, 0),
  X1 = c(0, 1, 2, 3, 0.1, 1.1, 2.1, 3.1),
  X2 = c(1, 1, 2, 2, 1.2, 1.1, 2.2, 2.1),
  Y  = c(3, 4, 5, 6, 2, 3, 4, 5)
)

m <- netmatch(
  data = dat,
  treat = "Z",
  covariates = c("X1", "X2"),
  network = A,
  method = "dual",
  kappa = 2
)

netmatch_test(m, outcome = "Y", method = "decay", eta = 0.03, rho = 0.10)
netmatch_sensitivity(m, outcome = "Y")
diagnose_match(m)
```

`kappa` is expressed in direct network-distance units. For example, `kappa = 2`
means units at graph distance 1 or 2 are considered too close for the relevant
network restriction.

## Reproduce the paper simulations

Use `fast = TRUE` for a quick smoke test.

```r
smoke <- reproduce_paper("table1", fast = TRUE)
paper_table1(smoke$table1)
```

Full simulation runs are intentionally not bundled as package data. They can be
regenerated:

```r
res <- reproduce_paper("all", n_rep = 500, output_dir = "paper-output")
```

The full run can take a long time because matching is repeated over many
replications, dependence levels, and treatment-effect settings.

## Solver note

The current portable backend is deterministic and dependency-light, which keeps
examples and smoke tests easy to run. The package interface reserves `solver =
"gurobi"` and `solver = "glpk"` for exact MIP backends; Gurobi remains the
preferred solver for reproducing the manuscript-scale optimized designs when
available.
