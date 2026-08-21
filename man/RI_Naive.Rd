\name{RI_Naive}
\alias{RI_Naive}
\alias{RI_Sensitivity}
\alias{RI_Design}
\title{Randomization-Based Inference for a Matched Design}
\usage{
RI_Naive(match, outcome, kappa = NULL, weight_type = c("ns", "ntc"))

RI_Sensitivity(match, outcome, eta = 0.03, rho = 0.10, kappa = NULL,
  weight_type = c("ns", "ntc"))

RI_Design(match, outcome, kappa = NULL, weight_type = c("ns", "ntc"))
}
\arguments{
\item{match}{A \code{netmatch} object.}
\item{outcome}{Name of a numeric outcome column with finite, non-missing values for all matched units.}
\item{eta}{One finite number from 0 to 1. Eta scales residual covariance between matched sets.}
\item{rho}{One finite number from 0 to 1. Rho controls how residual covariance changes with network distance.}
\item{kappa}{Analysis cutoff. Defaults to the matching cutoff stored in \code{match$kappa}.}
\item{weight_type}{Matched-set weighting. \code{"ns"} uses \eqn{1/(n_s+1)}, where \eqn{n_s} is the set size. \code{"ntc"} uses \eqn{1/(n_{ts}n_{cs})}, where \eqn{n_{ts}} and \eqn{n_{cs}} are the treated and control counts.}
}
\value{
A randomization-inference result object. Its \code{result} is a one-row data frame with the method, sensitivity values, cutoff, statistic, null expectation, variance, z-score, and p-value. The object also stores the matched-set details, covariance matrix, set distances, outcome, weighting choice, complete analysis options, and original match.
}
\description{
Runs randomization-based inference for a matched design. \code{RI_Naive()} uses only within-set variances. \code{RI_Sensitivity()} allows residual dependence between sets. \code{RI_Design()} uses the design covariance bound within the analysis cutoff.
}
\details{
For two matched sets at distance \eqn{d}, the sensitivity analysis multiplies the covariance bound by \eqn{\eta\rho^{d-1}}. Eta sets the overall size of residual dependence, and rho controls how it changes with distance.

The set statistic counts a treated outcome only when it is strictly greater than a control outcome. Tied outcomes add zero to that count. P-values use the current two-sided normal approximation.
}
\examples{
\dontrun{
sim <- simulate_netmatch_example()
m <- netmatch(sim$data, "Z", c("X1", "X2", "X3"), sim$net_dist,
              method = "dual", kappa = 2, solver = "auto")
RI_Naive(m, "Y")
RI_Sensitivity(m, "Y", eta = 0.03, rho = 0.10)
RI_Design(m, "Y")
}
}
