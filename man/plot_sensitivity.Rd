\name{plot_sensitivity}
\alias{plot_sensitivity}
\title{Plot Sensitivity Results}
\usage{
plot_sensitivity(x, type = NULL, alpha = 0.05, unadjusted = NULL,
  critical_ylim = NULL)

}
\arguments{
\item{x}{A sensitivity-grid or critical-sensitivity object.}
\item{type}{\code{"pvalue"} or \code{"critical"}. The default follows the class of \code{x}.}
\item{alpha}{Reference or test level. Default: \code{0.05}.}
\item{unadjusted}{Optional unadjusted p-value or test result for the p-value plot.}
\item{critical_ylim}{Optional y-axis limits for the critical plot. Default: \code{c(0, 1)}.}
}
\value{
A \code{ggplot} object.
}
\description{
Plots a p-value grid or a critical sensitivity curve. Stored analysis choices, including the weighting method, are reused for conversions and unadjusted reference lines.
}
\details{
A p-value plot shows how the p-value changes with \code{rho} for
each selected value of \code{eta}. A critical plot shows the boundary value of
\code{eta}; values below the curve retain rejection at the selected \code{alpha},
and values above it do not. If no finite boundary exists on the requested grid,
the function returns an error with the reason rather than drawing a curve.
Algebraic crossings outside the displayed range remain in the plot data but are
not critical contour points within the parameter range.
}
\examples{
\dontrun{
sim <- simulate_netmatch_example()
m <- netmatch(sim$data, "Z", c("X1", "X2", "X3"), sim$net_dist,
              method = "dual", kappa = 2, solver = "auto")
sens <- sensitivity_grid(m, "Y", eta = seq(0, 0.5, by = 0.1), rho = seq(0, 0.3, by = 0.05))
plot_sensitivity(sens)

crit <- critical_sensitivity(m, "Y")
plot_sensitivity(crit)
}
}
