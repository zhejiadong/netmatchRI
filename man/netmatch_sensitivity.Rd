\name{netmatch_sensitivity}
\alias{netmatch_sensitivity}
\title{Sensitivity Grid for Network Dependence}
\usage{
netmatch_sensitivity(match, outcome, eta = seq(0, 0.10, by = 0.01),
  rho = seq(0, 0.50, by = 0.05), kappa = NULL,
  weight_type = c("ns", "ntc"))
}
\arguments{
\item{match}{A \code{netmatch} object.}
\item{outcome}{Name of the outcome column.}
\item{eta}{Finite numeric values from 0 to 1. Eta scales residual covariance between matched sets.}
\item{rho}{Finite numeric values from 0 to 1. Rho controls how residual covariance changes with network distance.}
\item{kappa}{Analysis cutoff. Defaults to the matching cutoff stored in \code{match$kappa}.}
\item{weight_type}{Matched-set weighting. \code{"ns"} uses \eqn{1/(n_s+1)}. \code{"ntc"} uses \eqn{1/(n_{ts}n_{cs})}.}
}
\value{
A \code{netmatch_sensitivity} object with a p-value \code{grid}, the original \code{match}, and all analysis-defining options: outcome, eta values, rho values, cutoff, and weighting choice.
}
\description{
Evaluates the sensitivity analysis over an eta-rho grid. For matched sets at distance \eqn{d}, the covariance bound is multiplied by \eqn{\eta\rho^{d-1}}.
}
\details{
The p-values use a two-sided normal approximation. The set statistic counts a treated outcome only when it is strictly greater than a control outcome; ties add zero.
}
\examples{
\dontrun{
sim <- simulate_netmatch_example()
m <- netmatch(sim$data, "Z", c("X1", "X2", "X3"), sim$net_dist,
              method = "dual", kappa = 2, solver = "auto")
sens <- netmatch_sensitivity(m, "Y", eta = seq(0, 0.03, by = 0.01),
                             rho = seq(0, 0.50, by = 0.25))
sens$grid
}
}
