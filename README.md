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

## Example: A TRIP-Like Network

This example creates a small network with three loosely connected communities,
similar in spirit to the TRIP application. The data include treatment `Z`,
outcome `Y`, and covariates named like the TRIP analysis:
`edu`, `employment`, `ACMDT`, `baseRisk`, and `HIV`.

```r
library(netmatchRI)

set.seed(24)
n <- 24

# Adjacency matrix A for an undirected network.
A <- matrix(0, n, n)
for (b in 0:2) {
  ids <- (1:8) + 8 * b
  for (i in seq_along(ids)) {
    j <- ifelse(i == 8, 1, i + 1)
    A[ids[i], ids[j]] <- 1
    A[ids[j], ids[i]] <- 1
  }
}

# A few bridges between communities.
A[4, 12] <- A[12, 4] <- 1
A[12, 20] <- A[20, 12] <- 1
A[8, 16] <- A[16, 8] <- 1

community <- rep(c(0, 1, 2), each = 8)
edu <- sample(1:4, n, replace = TRUE)
employment <- sample(0:1, n, replace = TRUE)
ACMDT <- sample(0:1, n, replace = TRUE)
baseRisk <- 0.4 * community + rnorm(n)
HIV <- rbinom(n, 1, plogis(-1 + 0.4 * community))

# A toy binary treatment and outcome.
Z <- as.integer(seq_len(n) %% 2 == 0)
Y <- 0.3 * Z + 0.2 * edu + 0.4 * HIV + 0.1 * baseRisk + rnorm(n)

trip_like <- data.frame(
  Z = Z,
  Y = Y,
  edu = edu,
  employment = employment,
  ACMDT = ACMDT,
  baseRisk = baseRisk,
  HIV = HIV
)
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
  eta = c(0, 0.03, 0.06),
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
