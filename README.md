# netmatchRI

`netmatchRI` provides a user-friendly workflow for network-constrained
matching and randomization inference with observational network data. The
package follows the notation in the manuscript *Design and Analysis for Valid
Causal Inference with Network-Dependent Data*.

The package is designed for two tasks:

- applied analysis of network-dependent observational data;
- reproduction of the paper simulation studies and real-data illustrations.

## Main Workflow

The core workflow has four steps:

1. Start with a data frame containing a binary treatment `Z`, outcome `Y`,
   covariates `X`, and a network adjacency matrix `A`.
2. Build matched sets with `netmatch()`.
3. Test the sharp null with `netmatch_test()`.
4. Inspect robustness to residual across-set dependence with
   `netmatch_sensitivity()`.

The network threshold is written in direct graph-distance units. If `d_ij` is
the shortest-path distance between units `i` and `j`, then `kappa = 2` means
that pairs with `d_ij <= 2` are treated as network-proximate. There is no
inverse-distance `1 / d_ij` threshold in the user-facing API.

## Install Locally Before Publication

From R, install the local package from this folder:

```r
install.packages("remotes")
remotes::install_local("C:/Users/zheji/Brown Dropbox/Zhejia Dong/24 NetworkMatching/MatchingNetDepedence/MatchingNet/netmatchRI")
library(netmatchRI)
```

From PowerShell, run package checks before pushing or releasing:

```powershell
cd "C:\Users\zheji\Brown Dropbox\Zhejia Dong\24 NetworkMatching\MatchingNetDepedence\MatchingNet\netmatchRI"

& "C:\Program Files\R\R-4.5.0\bin\Rscript.exe" -e "testthat::test_local('.')"
& "C:\Program Files\R\R-4.5.0\bin\R.exe" CMD build .
& "C:\Program Files\R\R-4.5.0\bin\R.exe" CMD check netmatchRI_0.0.1.tar.gz --no-manual
```

A clean development check should end with:

```text
Status: OK
```

## Example: A 300-Unit TRIP-Like Network

This example creates a 300-unit network with four 75-unit communities, close to
the scale used in the simulation study. It is also similar in spirit to the
TRIP application because the data include treatment `Z`, outcome `Y`, and
covariates named like the TRIP analysis:
`edu`, `employment`, `ACMDT`, `baseRisk`, and `HIV`.

```r
library(netmatchRI)

set.seed(24)
n <- 300
n_block <- 4
block_size <- n / n_block
block <- rep(seq_len(n_block), each = block_size)

# Adjacency matrix A for an undirected stochastic-block-style network.
# Within-community ties are more likely than between-community ties.
p_in <- 0.19
p_out <- 0.003
A <- matrix(0, n, n)
for (i in seq_len(n - 1)) {
  for (j in (i + 1):n) {
    p_ij <- if (block[i] == block[j]) p_in else p_out
    A[i, j] <- A[j, i] <- rbinom(1, 1, p_ij)
  }
}

# Add a few deterministic bridges so the example network is connected enough
# for shortest-path distances d_ij to be useful.
A[75, 76] <- A[76, 75] <- 1
A[150, 151] <- A[151, 150] <- 1
A[225, 226] <- A[226, 225] <- 1

edu <- sample(1:4, n, replace = TRUE)
employment <- sample(0:1, n, replace = TRUE)
ACMDT <- sample(0:1, n, replace = TRUE)
baseRisk <- 0.35 * block + rnorm(n)
HIV <- rbinom(n, 1, plogis(-1.2 + 0.25 * block + 0.25 * baseRisk))

# A toy binary treatment with about 30% prevalence, as in the simulation study.
lin_z <- 0.15 * edu + 0.25 * employment + 0.35 * HIV +
  0.20 * baseRisk + 0.20 * block + rnorm(n)
Z <- as.integer(lin_z >= stats::quantile(lin_z, 0.70))

# A toy outcome under a nonzero treatment effect.
Y <- 0.30 * Z + 0.20 * edu + 0.40 * HIV + 0.10 * baseRisk +
  0.15 * block + rnorm(n)

trip_like <- data.frame(
  Z = Z,
  Y = Y,
  edu = edu,
  employment = employment,
  ACMDT = ACMDT,
  baseRisk = baseRisk,
  HIV = HIV
)

table(trip_like$Z)
```

Run the dual-penalty matched design. In manuscript notation, this uses
`kappa = 2`, so network-proximate units with `d_ij <= 2` are not allowed inside
the relevant matched-set relations.

```r
m_dual <- netmatch(
  data = trip_like,
  treat = "Z",
  covariates = c("edu", "employment", "ACMDT", "baseRisk", "HIV"),
  network = A,
  method = "dual",
  kappa = 2
)

m_dual
summary(m_dual)
```

Run randomization inference. The statistic `T` is the weighted sum of
within-set Mann-Whitney statistics. `method = "decay"` uses the model-assisted
variance bound with sensitivity parameters `eta` and `rho`, and truncates
across-set dependence beyond `d0`.

```r
fit_decay <- netmatch_test(
  match = m_dual,
  outcome = "Y",
  method = "decay",
  eta = 0.03,
  rho = 0.10,
  d0 = 2
)

fit_decay
```

Compare with naive randomization inference, which sets all across-set
covariances to zero:

```r
fit_naive <- netmatch_test(m_dual, outcome = "Y", method = "naive")
fit_naive
```

Evaluate an `(eta, rho)` grid:

```r
sens <- netmatch_sensitivity(
  match = m_dual,
  outcome = "Y",
  eta = c(0, 0.03),
  rho = c(0, 0.10),
  d0 = 2
)

sens$grid
```

Check matching diagnostics:

```r
diag <- diagnose_match(m_dual)
diag$covariate_smd
diag$average_within_distance
diag$within_distance_table
```

You can also compare the three matched designs:

```r
m_cov <- netmatch(trip_like, "Z",
                  c("edu", "employment", "ACMDT", "baseRisk", "HIV"),
                  A, method = "covariate", kappa = 2)

m_single <- netmatch(trip_like, "Z",
                     c("edu", "employment", "ACMDT", "baseRisk", "HIV"),
                     A, method = "single", kappa = 2)

m_dual <- netmatch(trip_like, "Z",
                   c("edu", "employment", "ACMDT", "baseRisk", "HIV"),
                   A, method = "dual", kappa = 2)

rbind(
  covariate = netmatch_test(m_cov, "Y", method = "decay",
                            eta = 0.03, rho = 0.10, d0 = 2)$result,
  single = netmatch_test(m_single, "Y", method = "decay",
                         eta = 0.03, rho = 0.10, d0 = 2)$result,
  dual = netmatch_test(m_dual, "Y", method = "decay",
                       eta = 0.03, rho = 0.10, d0 = 2)$result
)
```

## Reproduce Paper Simulation Smoke Tests

Use `fast = TRUE` to verify that the reproduction pipeline works in your local
environment:

```r
smoke <- reproduce_paper(
  task = "table1",
  fast = TRUE,
  output_dir = "paper-smoke-output"
)

paper_table1(smoke$table1)
```

Run all paper reproduction tasks in smoke-test mode:

```r
smoke_all <- reproduce_paper(
  task = "all",
  fast = TRUE,
  output_dir = "paper-smoke-output"
)
```

Full manuscript-scale reproduction uses `n_rep = 500` and can take a long time
because matching is repeated across replications, dependence levels, and
treatment-effect settings:

```r
res <- reproduce_paper(
  task = "all",
  n_rep = 500,
  output_dir = "paper-output"
)
```

Large simulation outputs are intentionally not bundled in the package.

## TRIP Real-Data Application

The package includes a wrapper for the local TRIP real-data file. The raw TRIP
data are not bundled in the package or pushed to GitHub.

The wrapper expects a local `.RData` file containing:

- `ds_use`: a data frame with `EGO_ID`, `Z`, `Y`, `edu`, `employment`,
  `ACMDT`, `baseRisk`, and `HIV`;
- `Adj`: a square adjacency matrix aligned with the TRIP units.

Run:

```r
trip <- trip_application(
  data_file = "../TRIP/data_clear.RData",
  output_dir = "trip-application-output"
)

trip$tests
trip$diagnostics
```

Or from PowerShell:

```powershell
Rscript inst/reproduce/trip-application.R "../TRIP/data_clear.RData" "trip-application-output"
```

## Install From the Private GitHub Repository

After you have access to the private GitHub repository:

```r
remotes::install_github("zhejiadong/netmatchRI")
library(netmatchRI)
```

The current repository is private. This keeps the package publishable and easy
to share with selected collaborators before a public release. GitHub Actions
are included for `R CMD check`.

## Solver Note

The current package backend is deterministic and dependency-light so examples,
tests, and smoke runs work without a licensed solver. The public interface
already reserves `solver = "gurobi"` and `solver = "glpk"` for exact MIP
backends. Gurobi remains the intended backend for manuscript-scale optimized
matching once the exact solver path is wired into the package internals.
