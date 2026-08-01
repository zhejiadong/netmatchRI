#' netmatchRI: Network Matching and Randomization Inference
#'
#' `netmatchRI` provides a compact workflow for matched designs and
#' randomization inference when observations may be related through a network.
#'
#' The main workflow is:
#'
#' * `netmatch()` builds the dual-penalty matched design from a data frame,
#'   binary treatment, covariates, and an adjacency or distance matrix.
#' * `diagnose_match()` summarizes covariate balance and within-set network
#'   distance in a matched design.
#' * `RI_naive()`, `RI_decay()`, and `RI_design()` run randomization inference
#'   for a matched design.
#' * `netmatch_sensitivity()` evaluates the randomization-inference result over
#'   a grid of residual network-dependence parameters.
#' * `critical_sensitivity()` computes the critical eta curve for the
#'   dual-penalty design.
#' * `plot_sensitivity()` plots p-value grids or the critical curve.
#' * `simulate_netmatch_example()` generates one 300-unit example
#'   dataset for examples and workflow checks.
#'
#' The dual-penalty matching backend defaults to open-source HiGHS. This choice
#' is based on package-local benchmarks at N = 300--500 and is not a claim that
#' HiGHS is universally fastest. `solver = "auto"` tries Gurobi, HiGHS, then
#' GLPK; either licensed Gurobi or open-source GLPK can also be requested
#' explicitly. At N = 500, HiGHS may return a validated time-limited incumbent,
#' in which case `netmatch()` issues a warning with the reported gap. After
#' installation, use `library(netmatchRI)`,
#' `ls("package:netmatchRI")`, `help(package = "netmatchRI")`,
#' `?netmatchRI`, and `?netmatch` to inspect the package.
#'
"_PACKAGE"

utils::globalVariables(c("eta_label", "eta_plot", "p_value", "rho"))
