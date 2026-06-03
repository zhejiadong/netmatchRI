\name{netmatch}
\alias{netmatch}
\title{Build a Network-Constrained Matched Design}
\usage{
netmatch(data, treat, covariates, network,
  method = c("dual", "single", "covariate"), kappa = 2,
  solver = c("auto", "gurobi", "glpk"), min_controls = 1,
  max_controls = Inf)
}
\arguments{
\item{data}{Data frame containing treatment and covariates.}
\item{treat}{Name of the binary treatment column coded 0/1.}
\item{covariates}{Character vector of covariate column names.}
\item{network}{Square adjacency or network-distance matrix.}
\item{method}{Matching method: \code{"dual"}, \code{"single"}, or \code{"covariate"}.}
\item{kappa}{Network-distance threshold. Pairs with distance less than or equal to \code{kappa} are disallowed where the method applies network restrictions.}
\item{solver}{Solver preference.}
\item{min_controls}{Minimum controls per treated set.}
\item{max_controls}{Maximum controls per treated set.}
}
\value{A \code{netmatch} object.}
\description{
Creates covariate-only, single-penalty, or dual-penalty matched designs using
direct graph-distance thresholds.
}
\examples{
A <- matrix(0, 8, 8)
A[cbind(1:7, 2:8)] <- 1
A[cbind(2:8, 1:7)] <- 1
dat <- data.frame(
  Z = c(1, 1, 1, 1, 0, 0, 0, 0),
  X1 = c(0, 1, 2, 3, 0.1, 1.1, 2.1, 3.1),
  X2 = c(1, 1, 2, 2, 1.2, 1.1, 2.2, 2.1)
)
m <- netmatch(dat, "Z", c("X1", "X2"), A, method = "dual", kappa = 2)
m
}
