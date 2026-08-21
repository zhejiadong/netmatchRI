\name{sensitivity_grid}
\alias{sensitivity_grid}
\title{Sensitivity Grid for Network Dependence}
\usage{
sensitivity_grid(match, outcome, eta = seq(0, 0.10, by = 0.01),
  rho = seq(0, 0.50, by = 0.05), kappa = NULL,
  weight_type = c("ns", "ntc"))
}
\arguments{
\item{match}{A \code{netmatch} object.}
\item{outcome}{Name of a numeric outcome column with finite, non-missing values for all matched units.}
\item{eta}{Finite numeric values from 0 to 1.}
\item{rho}{Finite numeric values from 0 to 1.}
\item{kappa}{Analysis cutoff. Defaults to the matching cutoff stored in \code{match$kappa}.}
\item{weight_type}{Matched-set weighting. \code{"ns"} uses \code{1 / (n_s + 1)}; \code{"ntc"} uses \code{1 / (n_ts * n_cs)}.}
}
\value{
A \code{netmatch_sensitivity} object with a p-value \code{grid}, source
\code{match}, and the complete analysis options.
}
\description{
Evaluates sensitivity analysis results over a grid of eta and rho values. Eta
scales residual covariance between matched sets. Rho controls how that
covariance changes with distance: for sets at distance \code{d}, the multiplier
is \code{eta * rho^(d - 1)}.
}
\details{
P-values use a two-sided normal approximation. The set statistic counts a
treated outcome only when it is strictly greater than a control outcome; tied
outcomes add zero.
}
\examples{
\dontrun{
sim <- simulate_netmatch_example()
m <- netmatch(sim$data, "Z", c("X1", "X2", "X3"), sim$net_dist,
              method = "dual", kappa = 2, solver = "auto")
sens <- sensitivity_grid(m, "Y", eta = seq(0, 0.03, by = 0.01),
                         rho = seq(0, 0.50, by = 0.25))
sens$grid
}
}
