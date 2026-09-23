\name{netmatch}
\alias{netmatch}
\title{Build a Network-Constrained Matched Design}
\usage{
netmatch(data, treat, covariates, network,
  method = c("dual", "single", "covariate"), kappa = 2,
  solver = c("highs", "auto", "gurobi", "glpk"), min_controls = 0.01,
  max_controls = 100, caliper = Inf, timelimit = 90, mipgap = 0.01,
  threads = 1, estimand = c("ATT", "ATC", "ATE"),
  network_type = c("auto", "adjacency", "distance"),
  include_solver = FALSE)

}
\arguments{
\item{data}{Data frame with one row per unit. The names \code{subclass} and \code{weights} are reserved.}
\item{treat}{Name of the binary treatment column.}
\item{covariates}{Names of the numeric or factor covariates used for matching.}
\item{network}{Symmetric adjacency or network-distance matrix aligned with \code{data}.}
\item{method}{Matching method: dual-penalty, single-penalty, or covariate-only. Default: dual-penalty (\code{"dual"}).}
\item{kappa}{Network distance threshold used by designs with network restrictions. Default: \code{2}.}
\item{solver}{Solver for dual-penalty matching. \code{"auto"} tries HiGHS, then Gurobi, then GLPK. Default: \code{"highs"}.}
\item{min_controls}{Minimum control-to-treated ratio. Default: \code{0.01}.}
\item{max_controls}{Maximum control-to-treated ratio; use \code{Inf} for no upper limit. Default: \code{100}.}
\item{caliper}{Mahalanobis-distance caliper for treat-control pairs, like propensity score caliper; use \code{Inf} for no caliper. Default: \code{Inf}.}
\item{timelimit}{Solver time limit in seconds. Default: \code{90}.}
\item{mipgap}{Target relative MIP gap. Default: \code{0.01}.}
\item{threads}{Requested solver threads. Default: \code{1}.}
\item{estimand}{Unit matching weights used in summaries and balance diagnostics: \code{"ATT"}, \code{"ATC"}, or \code{"ATE"}. Default: \code{"ATT"}. This choice does not change the sharp-null RI test.}
\item{network_type}{Format of \code{network}: \code{"auto"}, \code{"adjacency"}, or \code{"distance"}. Default: \code{"auto"}.}
\item{include_solver}{If \code{TRUE}, retain the complete solver result. Default: \code{FALSE}.}
}
\value{
\code{netmatch()} returns a \code{netmatch} object containing the matched data,
full-sample subclasses and unit weights, design settings, network distances, and
solver summary. The complete solver result is included only when requested.

\code{summary()} returns sample, balance, network, and solver summaries.
\code{plot()} draws the within-matched-set matching diagnostic.
}
\description{
Provides three network-constrained matching methods for observational network
data: dual-penalty, single-penalty, and covariate-only matching. Dual-penalty
matching minimizes total covariate dissimilarity while restricting
network-close pairs both across treatment arms and within matched sets.
Single-penalty matching restricts network-close pairs across treatment arms;
covariate-only matching has no network-distance restriction.
}
\details{
The dual-penalty matching is formulated and solved as a mixed-integer program. Covariate-only and single-penalty comparison
designs use \code{optmatch::fullmatch()}.
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
  solver = "highs",
  estimand = "ATT"
)
m
summary(m)
}
}
