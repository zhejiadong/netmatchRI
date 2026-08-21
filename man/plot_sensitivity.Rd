\name{plot_sensitivity}
\alias{plot_sensitivity}
\alias{plot.netmatch_sensitivity}
\alias{plot.netmatch_critical_sensitivity}
\title{Plot Sensitivity Results}
\usage{
plot_sensitivity(x, type = NULL, alpha = 0.05, naive = NULL,
  critical_ylim = NULL, ...)
}
\arguments{
\item{x}{A \code{netmatch_sensitivity} or \code{netmatch_critical_sensitivity} object.}
\item{type}{Plot type: \code{"pvalue"} or \code{"critical"}. By default, sensitivity grid objects use \code{"pvalue"} and critical objects use \code{"critical"}.}
\item{alpha}{Test level for the reference line.}
\item{naive}{Optional naive p-value or randomization-inference result object.}
\item{critical_ylim}{Optional coordinate limits for a critical plot. Defaults to \code{c(0, 1)}. Values outside the visible range remain unchanged in the plot data.}
\item{...}{Additional arguments passed to \code{critical_sensitivity()} when converting a grid object to a critical plot.}
}
\value{
A \code{ggplot} object. A p-value plot shows p-values against rho with one curve per eta. A critical plot shows critical eta against rho.
}
\description{
Plots a p-value grid or a critical sensitivity curve. Stored analysis choices, including the weighting method, are reused for conversions and naive reference lines.
}
\examples{
\dontrun{
sim <- simulate_netmatch_example()
m <- netmatch(sim$data, "Z", c("X1", "X2", "X3"), sim$net_dist,
              method = "dual", kappa = 2, solver = "auto")
sens <- netmatch_sensitivity(m, "Y", eta = c(0, 0.03), rho = seq(0, 1, by = 0.5))
plot_sensitivity(sens)

crit <- critical_sensitivity(m, "Y")
plot_sensitivity(crit)
}
}
