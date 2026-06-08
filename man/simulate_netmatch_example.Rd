\name{simulate_netmatch_example}
\alias{simulate_netmatch_example}
\title{Simulate the 300-Unit Network-Matching Example}
\usage{
simulate_netmatch_example(seed = 90141, beta_z = 0, n = 300,
  alpha1 = 0.9, alpha2 = 0.1, pin = 0.19, pout = 0.003)
}
\arguments{
\item{seed}{Random seed for the simulation.}
\item{beta_z}{Treatment effect in the outcome model.}
\item{n}{Number of units. The built-in example uses 300.}
\item{alpha1}{Dependence strength in \code{[0, 1]}; larger values induce stronger network dependence.}
\item{alpha2}{Individual variation in \code{[0, 1]} independent of the network.}
\item{pin}{Within-block edge probability for the stochastic block network.}
\item{pout}{Between-block edge probability for the stochastic block network.}
}
\value{
A list with:
\itemize{
\item \code{data}: data frame with \code{Y}, \code{Z}, \code{X1},
\code{X2}, and \code{X3}.
\item \code{Adj}: \code{n} by \code{n} adjacency matrix for the fixed network.
\item \code{net_dist}: graph-distance matrix for the fixed network.
\item \code{V}: variance-covariance matrix used by the data generator.
\item \code{seed}: simulation seed.
\item \code{beta_z}: treatment effect used in the outcome model.
\item \code{alpha1}: dependence strength used in the covariance construction.
\item \code{alpha2}: individual variation used in the covariance construction.
\item \code{pin}: within-block edge probability used to generate the network.
\item \code{pout}: between-block edge probability used to generate the network.
}
}
\description{
Generates one data set from the built-in simulation design: a fixed stochastic
block network, network-dependent covariates, treatment, and outcome. The default
values reproduce the package illustration.
}
\details{
\code{alpha1} controls the strength of network dependence, while \code{alpha2}
controls individual variation independent of the network. Both parameters are
restricted to \code{[0, 1]}. For illustration, the package examples use the
default values.
}
\examples{
sim <- simulate_netmatch_example()
dat <- sim$data
net_dist <- sim$net_dist

dim(dat)
table(dat$Z)
summary(lm(Y ~ Z, data = dat))$coefficients
}
