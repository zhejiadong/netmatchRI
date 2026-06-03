\docType{package}
\name{netmatchRI-package}
\alias{netmatchRI}
\alias{netmatchRI-package}
\title{netmatchRI: Network Matching and Randomization Inference}
\description{
\pkg{netmatchRI} provides a compact workflow for matched designs and
randomization inference when observations may be related through a network.
}
\details{
The main workflow is:
\itemize{
\item \code{\link{netmatch}()} builds the dual-penalty matched design from a
data frame, binary treatment, covariates, and an adjacency or distance matrix.
\item \code{\link{diagnose_match}()} summarizes covariate balance and
within-set network distance in a matched design.
\item \code{\link{RI_naive}()}, \code{\link{RI_decay}()}, and
\code{\link{RI_design}()} run randomization inference for a matched design.
\item \code{\link{netmatch_sensitivity}()} evaluates the
randomization-inference result over a grid of residual network-dependence
parameters.
\item \code{\link{critical_sensitivity}()} computes the critical eta curve
for the dual-penalty design.
\item \code{\link{plot_sensitivity}()} plots p-value grids or the critical
curve.
\item \code{\link{simulate_netmatch_example}()} generates one 300-unit
strong-dependence simulation data set for examples and workflow checks.
}

The dual-penalty matching backend uses Gurobi when available and can fall back
to open-source GLPK through \code{Rglpk}. After installation, use
\code{library(netmatchRI)},
\code{ls("package:netmatchRI")}, \code{help(package = "netmatchRI")},
\code{?netmatchRI}, and \code{?netmatch} to inspect the package.
}
