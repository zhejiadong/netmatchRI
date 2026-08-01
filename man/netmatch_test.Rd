\name{netmatch_test}
\alias{RI_naive}
\alias{RI_decay}
\alias{RI_design}
\alias{netmatch_test}
\title{Randomization-Based Inference Given a Matched Design}
\usage{
netmatch_test(match, outcome, method = c("decay", "naive", "design"),
  eta = 0.03, rho = 0.10, kappa = NULL,
  weight_type = c("ns", "ntc"))

RI_naive(match, outcome, kappa = NULL, weight_type = c("ns", "ntc"))

RI_decay(match, outcome, eta = 0.03, rho = 0.10, kappa = NULL,
  weight_type = c("ns", "ntc"))

RI_design(match, outcome, kappa = NULL, weight_type = c("ns", "ntc"))
}
\arguments{
\item{match}{A \code{netmatch} object.}
\item{outcome}{Name of the outcome column.}
\item{method}{Variance method used internally by \code{netmatch_test()}.}
\item{eta}{Magnitude parameter for sensitivity analysis.}
\item{rho}{Decay-rate parameter for sensitivity analysis.}
\item{kappa}{Analysis cutoff. Defaults to the matching cutoff stored in
\code{match$kappa}.}
\item{weight_type}{Weighting scheme for matched-set U-statistics.}
}
\value{
A \code{netmatch_test} object, which is a list with:
\itemize{
\item \code{result}: one-row data frame with method, eta, rho, kappa,
test statistic, null expectation, variance, z-score, and p-value.
\item \code{detail}: matched-set U-statistic details.
\item \code{covariance}: covariance matrix used for the test statistic.
\item \code{set_distance}: graph-distance matrix between matched sets.
\item \code{kappa}: analysis cutoff.
\item \code{match}: original \code{netmatch} object.
}
}
\description{
Runs normal-approximation randomization-based inference for a matched design.
\code{RI_decay()} supports sensitivity analysis, while \code{RI_design()}
provides the design-based analysis and keeps the distance truncation at
\code{kappa}: across-set
covariance is \eqn{M_{sk}} when \eqn{d_{sk} <= kappa} and zero otherwise. For
\code{RI_decay()}, the same bound is multiplied by \eqn{\eta \rho^{d_{sk}-1}}.
}
\examples{
\dontrun{
sim <- simulate_netmatch_example()
m <- netmatch(sim$data, "Z", c("X1", "X2", "X3"), sim$net_dist,
              method = "dual", kappa = 2, solver = "auto")
RI_naive(m, "Y")
RI_decay(m, "Y", eta = 0.03, rho = 0.10)
RI_design(m, "Y")
}
}
