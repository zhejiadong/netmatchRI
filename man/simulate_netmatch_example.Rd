\name{simulate_netmatch_example}
\alias{simulate_netmatch_example}
\title{Simulate the 300-Unit Network-Matching Example}
\usage{
simulate_netmatch_example(seed = 90141, beta_z = 0, n = 300)
}
\arguments{
\item{seed}{Random seed for the covariates, treatment, and outcome.}
\item{beta_z}{Treatment effect in the outcome model.}
\item{n}{Number of units. The built-in example uses 300.}
}
\value{
A list with:
\itemize{
\item \code{data}: data frame with \code{Y}, \code{Z}, \code{X1},
\code{X2}, and \code{X3}.
\item \code{Adj}: 300 by 300 adjacency matrix for the fixed network.
\item \code{net_dist}: graph-distance matrix for the fixed network.
\item \code{V}: strong-dependence covariance matrix used by the data
generator.
\item \code{seed}: simulation seed.
\item \code{dep_index}: dependence label, fixed at 3 for strong dependence.
\item \code{beta_z}: treatment effect used in the outcome model.
}
}
\description{
Generates one data set from the built-in strong-dependence simulation design:
a fixed 300-unit stochastic block network, network-dependent covariates,
treatment, and outcome.
}
\details{
This helper is intended for package evaluation and examples, especially when
searching for seeds that illustrate the p-value sensitivity plot and the
critical sensitivity curve.
}
\examples{
sim <- simulate_netmatch_example()
dat <- sim$data
net_dist <- sim$net_dist

dim(dat)
table(dat$Z)
summary(lm(Y ~ Z, data = dat))$coefficients
}
