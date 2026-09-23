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
\item{outcome}{Name of the numeric outcome column.}
\item{eta}{Sensitivity-scale values in \code{[0, 1]}. Default: \code{seq(0, 0.10, by = 0.01)}.}
\item{rho}{Sensitivity-decay values in \code{[0, 1]}. Default: \code{seq(0, 0.50, by = 0.05)}.}
\item{kappa}{Network distance threshold. \code{NULL} uses \code{match$kappa}.}
\item{weight_type}{Matched-set weighting: \code{"ns"} or \code{"ntc"}. Default: \code{"ns"}.}
}
\value{
A \code{netmatch_sensitivity} object containing the p-value grid, original
matched design, and analysis settings.
}
\description{
Evaluates how the p-value changes over a grid of eta and
rho values. Eta scales the design-based covariance bound between matched sets. Rho controls
how that bound decays with network distance: for sets at distance \code{d}, the
multiplier is \code{eta * rho^(d - 1)}.
}
\details{
The set statistic counts a treated outcome only when it is strictly greater than a control outcome.
The null moments assume no tied outcomes within a matched set.
}
\examples{
\dontrun{
sim <- simulate_netmatch_example()
m <- netmatch(sim$data, "Z", c("X1", "X2", "X3"), sim$net_dist,
              method = "dual", kappa = 2, solver = "auto")
sens <- sensitivity_grid(m, "Y", eta = seq(0, 0.5, by = 0.1),
                         rho = seq(0, 0.3, by = 0.05))
sens$grid
}
}
