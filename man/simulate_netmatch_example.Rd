\name{simulate_netmatch_example}
\alias{simulate_netmatch_example}
\title{Generate Example Network Data}
\usage{
simulate_netmatch_example(seed = 90141, beta_z = 0, n = 300,
  alpha1 = 0.9, alpha2 = 0.1, pin = 0.19, pout = 0.003)
}
\arguments{
\item{seed}{Simulation seed. Default: \code{90141}.}
\item{beta_z}{Treatment effect in the outcome model. Default: \code{0}.}
\item{n}{Number of units, divisible by four. Default: \code{300}.}
\item{alpha1}{Network-related covariance coefficient in \code{[0, 1]}. Default: \code{0.9}.}
\item{alpha2}{Unit-specific covariance coefficient in \code{[0, 1]}. Default: \code{0.1}.}
\item{pin}{Within-block edge probability. Default: \code{0.19}.}
\item{pout}{Between-block edge probability. Default: \code{0.003}.}
}
\value{
A list containing the simulated \code{data}, adjacency matrix \code{Adj},
network-distance matrix \code{net_dist}, covariance matrix \code{V}, and the
simulation settings.
}
\description{
Generates a stochastic block network and a data frame with three covariates, a
binary treatment, and an outcome. Covariates, treatment, and outcome may be
correlated across connected units.
}
\examples{
sim <- simulate_netmatch_example(seed = 123, n = 32)
dim(sim$data)
table(sim$data$Z)
}
