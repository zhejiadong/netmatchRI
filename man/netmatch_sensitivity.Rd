\name{netmatch_sensitivity}
\alias{netmatch_sensitivity}
\title{Sensitivity Grid for Network Dependence}
\usage{
netmatch_sensitivity(match, outcome, eta = seq(0, 0.10, by = 0.01),
  rho = seq(0, 0.50, by = 0.05), d0 = 2)
}
\arguments{
\item{match}{A \code{netmatch} object.}
\item{outcome}{Name of the outcome column.}
\item{eta}{Numeric vector of eta values.}
\item{rho}{Numeric vector of rho values.}
\item{d0}{Maximum set distance with nonzero across-set covariance.}
}
\value{A \code{netmatch_sensitivity} object.}
\description{
Evaluates model-assisted randomization inference over an eta-rho grid.
}
