\name{RI_design}
\alias{RI_design}
\title{Randomization-Based Inference with Design-Based Covariance Bounds}
\usage{
RI_design(match, outcome, kappa = NULL, weight_type = c("ns", "ntc"))
}
\arguments{
\item{match}{A \code{netmatch} object.}
\item{outcome}{Name of the numeric outcome column.}
\item{kappa}{Network distance threshold. \code{NULL} uses \code{match$kappa}.}
\item{weight_type}{Matched-set weighting: \code{"ns"} or \code{"ntc"}. Default: \code{"ns"}.}
}
\value{
A \code{netmatch_test} object containing the test result, matched-set
details, covariance and distance matrices, analysis settings, and original
matched design.
}
\description{
Tests the sharp null of no causal effect with design-based covariance bounds
for covariance between matched sets.
}
\details{
The test uses a weighted sum of Mann-Whitney U statistics from matched sets and
a normal approximation. It uses the design-based covariance bounds without
scaling or decay. The null moments assume no within-matched-set ties.
}
\examples{
\dontrun{
sim <- simulate_netmatch_example()
m <- netmatch(sim$data, "Z", c("X1", "X2", "X3"), sim$net_dist,
              method = "dual", kappa = 2, solver = "auto")
RI_design(m, "Y")
}
}