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
A list with:
\itemize{
\item \code{covariate_balance}: one row per covariate with the matched
absolute standardized mean difference.
\item \code{network_distance}: one-row summary with the number of within-set
unit pairs and the minimum, mean, and maximum within-set graph distance.
\item \code{within_distance_table}: counts and proportions for each observed
within-set graph distance.
\item \code{covariate_smd}: backward-compatible alias of
\code{covariate_balance}.
\item \code{average_within_distance}: backward-compatible alias of the mean
within-set graph distance.
}
}
\description{
Computes compact covariate balance and within-set network-distance diagnostics.
}
\examples{
\donttest{
sim <- simulate_netmatch_example()
m <- netmatch(sim$data, "Z", c("X1", "X2", "X3"), sim$net_dist,
              method = "dual", kappa = 2, solver = "auto")
diag <- diagnose_match(m)
diag$covariate_balance
diag$network_distance
diag$within_distance_table
}
}
