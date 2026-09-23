\name{RI_adjusted}
\alias{RI_adjusted}
\title{Randomization-Based Inference Adjusting for Covariance}
\usage{
RI_adjusted(match, outcome, eta = 0.03, rho = 0.10, kappa = NULL,
  weight_type = c("ns", "ntc"))
}
\arguments{
\item{match}{A \code{netmatch} object.}
\item{outcome}{Name of the numeric outcome column.}
\item{eta}{Relative magnitude of the design-based covariance bound in \code{[0, 1]}. Default: \code{0.03}.}
\item{rho}{Decay parameter in \code{[0, 1]}. Default: \code{0.10}.}
\item{kappa}{Network distance threshold. \code{NULL} uses \code{match$kappa}.}
\item{weight_type}{Matched-set weighting: \code{"ns"} or \code{"ntc"}. Default: \code{"ns"}.}
}
\value{
A \code{netmatch_test} object containing the test result, matched-set
details, covariance and distance matrices, analysis settings, and original
matched design.
}
\description{
Tests the sharp null of no causal effect while adjusting for covariance between
matched sets.
}
\details{
The test uses a weighted sum of Mann-Whitney U statistics from matched sets and
a normal approximation. The design-based covariance bound between matched sets
is scaled by \code{eta} and decays with network distance according to
\code{rho}. The null moments assume no within-matched-set ties.
}
\examples{
\dontrun{
sim <- simulate_netmatch_example()
m <- netmatch(sim$data, "Z", c("X1", "X2", "X3"), sim$net_dist,
              method = "dual", kappa = 2, solver = "auto")
RI_adjusted(m, "Y", eta = 0.03, rho = 0.10)
}
}