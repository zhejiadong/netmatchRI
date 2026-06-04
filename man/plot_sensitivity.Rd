\name{plot_sensitivity}
\alias{plot_sensitivity}
\alias{plot.netmatch_sensitivity}
\alias{plot.netmatch_critical_sensitivity}
\title{Plot Sensitivity Results}
\usage{
plot_sensitivity(x, type = c("critical", "pvalue"), alpha = 0.05,
  naive = NULL, critical_ylim = NULL, ...)
}
\arguments{
\item{x}{A \code{netmatch_sensitivity} or
\code{netmatch_critical_sensitivity} object.}
\item{type}{Plot type: \code{"critical"} or \code{"pvalue"}.}
\item{alpha}{Test level for the reference line.}
\item{naive}{Optional naive p-value or randomization-inference result object.}
\item{critical_ylim}{Optional y-axis limits for \code{type = "critical"}.
Defaults to an adaptive range that starts at zero. Use \code{c(0, 1)} to force
the full sensitivity-parameter range.}
\item{...}{Additional arguments passed to \code{critical_sensitivity()} when
\code{x} is a grid object and \code{type = "critical"}.}
}
\value{
A \code{ggplot} object. For \code{type = "pvalue"}, the plot shows p-values
against rho with one curve per eta. For \code{type = "critical"}, the plot
shows the critical eta boundary against rho.
}
\description{
Plots either the p-value grid from \code{netmatch_sensitivity()} or the
critical eta curve from \code{critical_sensitivity()}.
}
\examples{
\donttest{
sim <- simulate_netmatch_example()
m <- netmatch(sim$data, "Z", c("X1", "X2", "X3"), sim$net_dist,
              method = "dual", kappa = 2, solver = "auto")
sens <- netmatch_sensitivity(m, "Y", eta = c(0, 0.03), rho = c(0, 0.5, 1))
plot_sensitivity(sens, type = "pvalue")

crit <- critical_sensitivity(m, "Y")
plot_sensitivity(crit, type = "critical")
}
}
