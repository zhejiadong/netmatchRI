\name{plot_covariate_similarity}
\alias{plot_covariate_similarity}
\title{Plot Covariate Similarity Across Matching Designs}
\usage{
plot_covariate_similarity(designs, include_before = TRUE)
}
\arguments{
\item{designs}{Named list of comparable \code{netmatch} objects.}
\item{include_before}{Include the unmatched sample. Default: \code{TRUE}.}
}
\value{
A \code{ggplot} object.
}
\description{
Compares average treated-control Mahalanobis distance before matching and within
matched sets across selected matching designs. Smaller values indicate greater
covariate similarity.
}
\examples{
\dontrun{
sim <- simulate_netmatch_example(seed = 123, n = 32)
m <- netmatch(sim$data, "Z", c("X1", "X2", "X3"), sim$net_dist,
              method = "dual", kappa = 2)
plot_covariate_similarity(list("Dual-penalty" = m))
}
}
