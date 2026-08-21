\name{diagnose_match}
\alias{diagnose_match}
\title{Diagnose a Matched Design}
\usage{
diagnose_match(match)
}
\arguments{
\item{match}{A \code{netmatch} object.}
}
\value{
A list with exactly three fields:
\itemize{
\item \code{covariate_balance}: columns \code{covariate}, \code{before_abs_smd}, and \code{after_abs_smd} for each model-matrix covariate term.
\item \code{network_summary}: one-row summary with total, finite, and disconnected within-set pair counts and finite-distance minimum, mean, and maximum.
\item \code{within_distance_table}: counts and proportions for each within-set distance, including \code{Inf} for disconnected pairs when present.
}
Distance summaries are \code{NA} when there are no finite within-set distances.
}
\description{
Computes covariate balance and within-set network-distance diagnostics.
}
\details{
Factor covariates are expanded into model-matrix indicator columns. For every expanded term, the before and after absolute standardized mean differences use one fixed pooled original-sample standard deviation,
\deqn{\sqrt{\{(n_1-1)s_1^2 + (n_0-1)s_0^2\}/(n_1+n_0-2)}.}
The before difference uses unweighted original-sample treatment-group means. The after difference uses the unit matching weights among matched units.
}
\examples{
\dontrun{
sim <- simulate_netmatch_example()
m <- netmatch(sim$data, "Z", c("X1", "X2", "X3"), sim$net_dist,
              method = "dual", kappa = 2, solver = "auto")
diag <- diagnose_match(m)
diag$covariate_balance
diag$network_summary
diag$within_distance_table
}
}
