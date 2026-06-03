\name{netmatch}
\alias{netmatch}
\alias{summary.netmatch}
\alias{plot.netmatch}
\title{Build a Dual-Penalty Matched Design}
\usage{
netmatch(data, treat, covariates, network,
  method = c("dual", "single", "covariate"), kappa = 2,
  solver = c("auto", "gurobi", "glpk"), min_controls = 0.01,
  max_controls = 100, caliper = 8, timelimit = 50, mipgap = 0.01,
  threads = 1)

\method{summary}{netmatch}(object, ...)

\method{plot}{netmatch}(x, ...)
}
\arguments{
\item{data}{A data frame with one row per unit. Users supplying their own
data should include at least a binary treatment column such as \code{Z},
covariate columns such as \code{X1}, \code{X2}, \code{X3}, and an outcome
column such as \code{Y} when randomization inference will be run later.}
\item{treat}{Name of the binary treatment column coded 0/1.}
\item{covariates}{Character vector of covariate column names.}
\item{network}{Square adjacency or network-distance matrix.}
\item{method}{Matching method: \code{"dual"}, \code{"single"}, or \code{"covariate"}.}
\item{kappa}{Network-distance threshold. Pairs with distance less than or equal to \code{kappa} are disallowed where the method applies network restrictions.}
\item{solver}{Solver preference for \code{method = "dual"}. \code{"auto"} tries
Gurobi first, then the open-source GLPK backend through Rglpk. \code{"gurobi"}
uses Gurobi only. \code{"glpk"} uses Rglpk with a sparse triplet constraint
matrix. \code{"covariate"} and \code{"single"} use
\code{optmatch::fullmatch()}.}
\item{min_controls}{Minimum controls per treated unit in each matched set.}
\item{max_controls}{Maximum controls per treated unit in each matched set.}
\item{caliper}{Mahalanobis-distance caliper for feasible treated-control
edges.}
\item{timelimit}{Solver time limit in seconds.}
\item{mipgap}{Gurobi MIP gap.}
\item{threads}{Number of Gurobi threads.}
\item{object}{A \code{netmatch} object.}
\item{x}{A \code{netmatch} object.}
\item{...}{Unused.}
}
\value{
For \code{netmatch()}, a \code{netmatch} object, which is a list with:
\itemize{
\item \code{data}: matched data frame with an added \code{subclass} column.
\item \code{treat}: treatment column name.
\item \code{covariates}: covariates used for matching.
\item \code{method}: matching method used.
\item \code{method_label}: human-readable method label.
\item \code{kappa}: graph-distance cutoff used by the design.
\item \code{solver}: backend used for matching.
\item \code{network_distance}: graph-distance matrix used by the analysis.
\item \code{match_table}: unit-level table with subclass and treatment.
\item \code{solver_result}: raw solver or \pkg{optmatch} result.
}

For \code{summary.netmatch()}, a data frame with one row per matched set and
columns for subclass, total set size, treated count, and control count.

For \code{plot.netmatch()}, the diagnostics list from
\code{\link{diagnose_match}()} is returned invisibly after drawing the
within-set network-distance bar plot.
}
\description{
Builds the dual-penalty matched design via a mixed-integer program. The
covariate-only and single-penalty comparison designs use \pkg{optmatch}.
}
\examples{
\donttest{
if (requireNamespace("igraph", quietly = TRUE) &&
    requireNamespace("mvtnorm", quietly = TRUE) &&
    (requireNamespace("gurobi", quietly = TRUE) ||
     (requireNamespace("Rglpk", quietly = TRUE) &&
      requireNamespace("slam", quietly = TRUE)))) {
  sim <- simulate_netmatch_example()

  m <- netmatch(
    data = sim$data,
    treat = "Z",
    covariates = c("X1", "X2", "X3"),
    network = sim$net_dist,
    method = "dual",
    kappa = 2,
    solver = "auto"
  )
  m
  summary(m)
}
}
}
