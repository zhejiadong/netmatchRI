\name{simulate_netmatch_example}
\alias{simulate_netmatch_example}
\title{Generate Example Network Data}
\usage{
simulate_netmatch_example(seed = 90141, beta_z = 0, n = 300,
  alpha1 = 0.9, alpha2 = 0.1, pin = 0.19, pout = 0.003)
}
\arguments{
\item{seed}{Random seed for the simulation.}
\item{beta_z}{Treatment effect in the outcome model.}
\item{n}{Number of units. Must be divisible by four.}
\item{alpha1}{Network-related contribution to the covariance structure, from 0 to 1.}
\item{alpha2}{Unit-specific contribution to the covariance structure, from 0 to 1.}
\item{pin}{Within-block edge probability for the stochastic block network.}
\item{pout}{Between-block edge probability for the stochastic block network.}
}
\value{
A list with:
\itemize{
\item \code{data}: data frame with \code{Y}, \code{Z}, \code{X1},
\code{X2}, and \code{X3}.
\item \code{Adj}: \code{n} by \code{n} adjacency matrix.
\item \code{net_dist}: graph-distance matrix.
\item \code{V}: covariance matrix used to generate the data.
\item \code{seed}: simulation seed.
\item \code{beta_z}: treatment effect used in the outcome model.
\item \code{alpha1}: network-related covariance contribution.
\item \code{alpha2}: unit-specific covariance contribution.
\item \code{pin}: within-block edge probability.
\item \code{pout}: between-block edge probability.
}
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
