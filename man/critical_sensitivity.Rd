\name{critical_sensitivity}
\alias{critical_sensitivity}
\title{Critical Sensitivity Curve for Dual-Penalty Matching}
\usage{
critical_sensitivity(match, outcome, rho = seq(0, 1, by = 0.01),
  alpha = 0.05, kappa = NULL, weight_type = c("ns", "ntc"))
}
\arguments{
\item{match}{A \code{netmatch} object, ideally from \code{method = "dual"}.}
\item{outcome}{Name of the outcome column.}
\item{rho}{Numeric vector of rho values.}
\item{alpha}{Test level.}
\item{kappa}{Analysis cutoff. Defaults to \code{match$kappa}.}
\item{weight_type}{Weighting scheme for matched-set U-statistics.}
}
\value{
A \code{netmatch_critical_sensitivity} object, which is a list with:
\itemize{
\item \code{curve}: data frame with \code{rho} and \code{eta_critical}.
\item \code{alpha}: test level.
\item \code{kappa}: analysis cutoff.
\item \code{statistic}: observed weighted statistic.
\item \code{expectation}: null expectation.
\item \code{diagonal_variance}: independence variance component.
\item \code{naive}: naive randomization-inference result.
\item \code{interpretation}: text summary of the critical curve.
\item \code{match}: original \code{netmatch} object.
\item \code{outcome}: outcome column name.
}
}
\description{
Computes the critical value of eta as a function of rho that solves
\code{p(eta, rho) = alpha} for the sensitivity analysis.
}
\details{
The curve is intended as the robustness display for the dual-penalty design.
For each rho value, eta values below the critical curve
retain rejection at level \code{alpha}; values above the curve allow enough
residual network dependence to explain away the result.

If the independence/naive variance is already not significant at
\code{alpha}, there is no positive critical eta and \code{eta_critical} is
\code{NA}.
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
