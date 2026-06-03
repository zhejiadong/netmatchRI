\name{netmatch_test}
\alias{netmatch_test}
\title{Test a Network-Matched Design}
\usage{
netmatch_test(match, outcome, method = c("decay", "naive", "design"),
  eta = 0.03, rho = 0.10, d0 = 2)
}
\arguments{
\item{match}{A \code{netmatch} object.}
\item{outcome}{Name of the outcome column.}
\item{method}{Variance method.}
\item{eta}{Decay magnitude for model-assisted bounds.}
\item{rho}{Decay rate for model-assisted bounds.}
\item{d0}{Maximum set distance with nonzero across-set covariance.}
}
\value{A \code{netmatch_test} object with tidy test results.}
\description{
Runs normal-approximation randomization inference for a matched design.
}
