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
\item{eta}{Numeric vector of eta values.}
\item{rho}{Numeric vector of rho values.}
\item{kappa}{Analysis cutoff. Defaults to the matching cutoff stored in
\code{match$kappa}.}
\item{weight_type}{Weighting scheme for matched-set U-statistics.}
}
\value{
A \code{netmatch_sensitivity} object, which is a list with:
\itemize{
\item \code{grid}: data frame with one row per eta-rho pair and the resulting
test statistic, variance, z-score, and p-value.
\item \code{match}: original \code{netmatch} object.
\item \code{outcome}: outcome column name.
\item \code{kappa}: analysis cutoff.
}
}
\description{
Evaluates model-assisted randomization inference over an eta-rho grid using the
matching cutoff \code{kappa} for covariance truncation.
}
\examples{
\donttest{
sim <- simulate_netmatch_example()
m <- netmatch(sim$data, "Z", c("X1", "X2", "X3"), sim$net_dist,
              method = "dual", kappa = 2, solver = "auto")
sens <- netmatch_sensitivity(
  m,
  "Y",
  eta = c(0, 0.01, 0.03),
  rho = c(0, 0.25, 0.50)
)
sens$grid
}
}
