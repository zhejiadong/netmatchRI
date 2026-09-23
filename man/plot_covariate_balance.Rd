\name{plot_covariate_balance}
\alias{plot_covariate_balance}
\title{Plot Covariate Balance Across Matching Designs}
\usage{
plot_covariate_balance(designs, include_before = TRUE)
}
\arguments{
\item{designs}{Named list of comparable \code{netmatch} objects.}
\item{include_before}{Include the unmatched sample. Default: \code{TRUE}.}
}
\value{
A \code{ggplot} object.
}
\description{
Compares absolute standardized mean differences before matching and across selected matching designs.
}
\examples{
\dontrun{
sim <- simulate_netmatch_example(seed = 123, n = 32)
m <- netmatch(sim$data, "Z", c("X1", "X2", "X3"), sim$net_dist,
              method = "dual", kappa = 2)
plot_covariate_balance(list("Dual-penalty" = m))
}
}
