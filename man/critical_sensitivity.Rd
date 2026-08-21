\name{critical_sensitivity}
\alias{critical_sensitivity}
\title{Critical Sensitivity Curve for Dual-Penalty Matching}
\usage{
critical_sensitivity(match, outcome, rho = seq(0, 1, by = 0.01),
  alpha = 0.05, kappa = NULL, weight_type = c("ns", "ntc"))
}
\arguments{
\item{match}{A \code{netmatch} object, ideally from \code{method = "dual"}.}
\item{outcome}{Name of a numeric outcome column with finite, non-missing values for all matched units.}
\item{rho}{Finite numeric values from 0 to 1. Rho controls how residual covariance changes with network distance.}
\item{alpha}{Test level.}
\item{kappa}{Analysis cutoff. Defaults to \code{match$kappa}.}
\item{weight_type}{Matched-set weighting. \code{"ns"} uses \eqn{1/(n_s+1)}. \code{"ntc"} uses \eqn{1/(n_{ts}n_{cs})}.}
}
\value{
A \code{netmatch_critical_sensitivity} object. The \code{curve} shows the critical eta for each rho and whether that eta is between 0 and 1. The object also retains the analysis settings and supporting results.
}
\description{
Computes the critical eta at each rho where the two-sided normal-approximation p-value equals \code{alpha}.
}
\details{
For matched sets at distance \eqn{d}, residual covariance is multiplied by \eqn{\eta\rho^{d-1}}. Eta values below the critical curve retain rejection at level \code{alpha}; values above it allow enough residual dependence to explain away the result.

Critical eta values above 1 are kept in \code{eta_critical} and marked \code{FALSE} in \code{eta_in_range}. If the naive analysis is already not significant, there is no positive critical eta and \code{eta_critical} is \code{NA}.
}
\examples{
\dontrun{
sim <- simulate_netmatch_example()
m <- netmatch(sim$data, "Z", c("X1", "X2", "X3"), sim$net_dist,
              method = "dual", kappa = 2, solver = "auto")
crit <- critical_sensitivity(m, "Y", rho = seq(0, 1, by = 0.25))
crit$curve
crit$interpretation
}
}
