#' netmatchRI: Dual-Penalty Matching and Randomization-Based Inference
#'
#' `netmatchRI` provides matching and randomization-based inference for
#' observational network data with network dependence and interference.
#'
#' Use `netmatch()` to construct a dual-penalty, single-penalty, or
#' covariate-only matched design. Inspect the result with `summary()`,
#' `matched_data()`, and `diagnose_match()`. Use `RI_Naive()`,
#' `RI_Sensitivity()`, or `RI_Design()` for inference, and use
#' `sensitivity_grid()`, `critical_sensitivity()`, and `plot_sensitivity()` to
#' examine sensitivity to residual network dependence.
#'
#' Dual-penalty matching uses the open-source HiGHS solver by default. The
#' `solver` argument can select HiGHS, Gurobi, or GLPK; `solver = "auto"` tries
#' them in that order of availability.
#'
#' `simulate_netmatch_example()` generates example data for learning the
#' workflow.
#'
"_PACKAGE"

utils::globalVariables(c("eta_critical", "eta_label", "p_value", "rho"))
