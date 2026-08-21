\docType{package}
\name{netmatchRI-package}
\alias{netmatchRI}
\alias{netmatchRI-package}
\title{netmatchRI: Dual-Penalty Matching and Randomization-Based Inference}
\description{
\pkg{netmatchRI} provides matching and randomization-based inference for
observational network data with network dependence and interference.
}
\details{
Use \code{\link{netmatch}()} to construct a dual-penalty, single-penalty, or
covariate-only matched design. Inspect the result with \code{summary()},
\code{\link{matched_data}()}, and \code{\link{diagnose_match}()}. Use
\code{\link{RI_Naive}()}, \code{\link{RI_Sensitivity}()}, or
\code{\link{RI_Design}()} for inference, and use
\code{\link{sensitivity_grid}()}, \code{\link{critical_sensitivity}()}, and
\code{\link{plot_sensitivity}()} to examine sensitivity to residual network
dependence.

Dual-penalty matching uses the open-source HiGHS solver by default. The
\code{solver} argument can select HiGHS, Gurobi, or GLPK;
\code{solver = "auto"} tries them in that order of availability.

\code{\link{simulate_netmatch_example}()} generates example data for learning
the workflow.
}
