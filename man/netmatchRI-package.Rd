\docType{package}
\name{netmatchRI-package}
\alias{netmatchRI}
\title{netmatchRI: Network-Constrained Matching and Randomization-Based Inference}
\description{
\pkg{netmatchRI} provides network-constrained matching and randomization-based
inference for observational network data with network dependence and
interference. The three matching methods are dual-penalty, single-penalty, and
covariate-only matching.
}
\details{
Use \code{\link{netmatch}()} to construct a dual-penalty, single-penalty, or
covariate-only matched design. Inspect the result with \code{summary()},
\code{\link{matched_data}()}, and \code{\link{diagnose_match}()}, which checks
covariate balance and network distance. Compare selected
matching designs with \code{\link{plot_covariate_balance}()} and
\code{\link{plot_covariate_similarity}()}. Use \code{\link{RI_unadjusted}()},
\code{\link{RI_adjusted}()}, or \code{\link{RI_design}()} for inference, and use
\code{\link{sensitivity_grid}()}, \code{\link{critical_sensitivity}()}, and
\code{\link{plot_sensitivity}()} to examine sensitivity to residual network
dependence.

Dual-penalty matching uses the open-source HiGHS solver by default. The
\code{solver} argument can select HiGHS, Gurobi, or GLPK;
\code{solver = "auto"} tries them in that order of availability.

\code{\link{simulate_netmatch_example}()} generates the example network data
used in the documentation.
}
