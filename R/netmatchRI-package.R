#' netmatchRI: Network-Constrained Matching and Randomization-Based Inference
#'
#' `netmatchRI` provides network-constrained matching and randomization-based
#' inference for observational network data with network dependence and
#' interference. The three matching methods are dual-penalty, single-penalty,
#' and covariate-only matching.
#'
#' Use `netmatch()` to construct a dual-penalty, single-penalty, or
#' covariate-only matched design. Inspect the result with `summary()`,
#' `matched_data()`, and `diagnose_match()`, which checks covariate balance and
#' network distance. Compare selected designs with
#' `plot_covariate_balance()` and `plot_covariate_similarity()`. Use `RI_unadjusted()`,
#' `RI_adjusted()`, or `RI_design()` for inference, and use
#' `sensitivity_grid()`, `critical_sensitivity()`, and `plot_sensitivity()` to
#' examine sensitivity to network dependence across matched sets.
#'
#' Dual-penalty matching uses the open-source HiGHS solver by default. The
#' `solver` argument can select HiGHS, Gurobi, or GLPK; `solver = "auto"` tries
#' them in that order of availability.
#'
#' `simulate_netmatch_example()` generates the example network data used in the
#' documentation.
#'
"_PACKAGE"

utils::globalVariables(c(
  "abs_smd", "average_mahalanobis", "covariate", "design", "eta_critical",
  "eta_label", "label", "p_value", "rho"
))
