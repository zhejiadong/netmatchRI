\name{netmatch}
\alias{netmatch}
\alias{summary.netmatch}
\alias{plot.netmatch}
\title{Build a Dual-Penalty Matched Design}
\usage{
netmatch(data, treat, covariates, network,
  method = c("dual", "single", "covariate"), kappa = 2,
  solver = c("highs", "auto", "gurobi", "glpk"), min_controls = 0.01,
  max_controls = 100, caliper = 8, timelimit = 90, mipgap = 0.01,
  threads = 1, estimand = c("ATT", "ATC", "ATE"),
  network_type = c("auto", "adjacency", "distance"),
  include_solver = FALSE)

\method{summary}{netmatch}(object, ...)

\method{plot}{netmatch}(x, ...)
}
\arguments{
\item{data}{A data frame with one row per unit.}
\item{treat}{Name of the binary treatment column coded 0/1.}
\item{covariates}{Character vector of covariate column names. Factors are allowed.}
\item{network}{Numeric square adjacency or network-distance matrix.}
\item{method}{Matching method: \code{"dual"}, \code{"single"}, or \code{"covariate"}.}
\item{kappa}{Network-distance threshold. Pairs with distance less than or equal to \code{kappa} are disallowed where the method applies network restrictions.}
\item{solver}{Solver preference for \code{method = "dual"}. \code{"highs"} is the default open-source backend. \code{"auto"} tries Gurobi, HiGHS, then GLPK. Comparison methods use \code{optmatch::fullmatch()}.}
\item{min_controls}{Minimum controls per treated unit in each matched set.}
\item{max_controls}{Maximum controls per treated unit in each matched set.}
\item{caliper}{Mahalanobis-distance caliper for feasible treated-control edges.}
\item{timelimit}{Solver time limit in seconds.}
\item{mipgap}{Relative MIP gap for Gurobi and HiGHS.}
\item{threads}{Number of Gurobi or HiGHS threads.}
\item{estimand}{Target estimand used for unit matching weights: \code{"ATT"}, \code{"ATC"}, or \code{"ATE"}. Nonzero weights are normalized within each treatment group to have mean 1.}
\item{network_type}{Interpretation of \code{network}: \code{"auto"}, \code{"adjacency"}, or \code{"distance"}. An explicit choice overrides automatic detection. Infinite values are allowed only for distance matrices.}
\item{include_solver}{If \code{TRUE}, retain the raw backend result in \code{solver_result}. The default retains only stable \code{solver_info}.}
\item{object}{A \code{netmatch} object.}
\item{x}{A \code{netmatch} object.}
\item{...}{Unused.}
}
\value{
For \code{netmatch()}, a \code{netmatch} object containing matched \code{data}; full-length \code{subclass}, \code{weights}, and \code{matched} vectors; the original \code{call} and \code{original_data}; \code{estimand}; resolved \code{network_type}; \code{network_distance}; and stable \code{solver_info}. The raw \code{solver_result} is present only when \code{include_solver = TRUE}.

For \code{summary.netmatch()}, a \code{netmatch_summary} object with sample counts, effective sample sizes, matched-set sizes, covariate balance, network summary, and solver information.

For \code{plot.netmatch()}, the diagnostics list from \code{\link{diagnose_match}()} is returned invisibly after drawing the within-set network-distance bar plot.
}
\description{
Builds the dual-penalty matched design via a mixed-integer program. The covariate-only and single-penalty comparison designs use \pkg{optmatch}.
}
\details{
\code{solver}, \code{caliper}, \code{timelimit}, \code{mipgap}, and \code{threads} are used only by \code{method = "dual"}. \code{kappa} is ignored by \code{method = "covariate"}. \code{method = "single"} uses \code{kappa} but not the mixed-integer solver controls.
}
\examples{
\dontrun{
sim <- simulate_netmatch_example()
m <- netmatch(
  data = sim$data,
  treat = "Z",
  covariates = c("X1", "X2", "X3"),
  network = sim$net_dist,
  method = "dual",
  kappa = 2,
  solver = "highs"
)
m
summary(m)
}
}
