#' Test a Network-Matched Design
#'
#' Runs normal-approximation randomization inference for a `netmatch` object.
#' The returned data frame contains the observed weighted Mann-Whitney statistic,
#' null expectation, variance, z-score, and two-sided p-value.
#'
#' @param match A `netmatch` object.
#' @param outcome Name of the outcome column.
#' @param method Variance method: `"decay"`, `"naive"`, or `"design"`.
#' @param eta Decay magnitude for model-assisted bounds.
#' @param rho Decay rate for model-assisted bounds.
#' @param d0 Maximum set distance with nonzero across-set covariance.
#' @return A `netmatch_test` object.
#' @export
netmatch_test <- function(match,
                          outcome,
                          method = c("decay", "naive", "design"),
                          eta = 0.03,
                          rho = 0.10,
                          d0 = 2) {
  if (!inherits(match, "netmatch")) {
    stop("`match` must be a netmatch object.", call. = FALSE)
  }
  method <- match.arg(method)
  if (!outcome %in% names(match$data)) stop("`outcome` column not found.", call. = FALSE)

  obs <- .observed_stats(match$data, outcome, match$treat, "subclass")
  unit_ids <- as.integer(rownames(match$data))
  set_dist <- .set_distance_matrix(match$network_distance, as.integer(match$data$subclass), unit_ids)
  Sigma <- .covariance_matrix(obs$detail, set_dist, method, eta, rho, d0)
  w <- obs$detail$weight
  variance <- as.numeric(t(w) %*% Sigma %*% w)
  pv <- .normal_pvalue(obs$statistic, obs$expectation, variance)

  result <- data.frame(
    method = method,
    eta = if (method == "decay") eta else NA_real_,
    rho = if (method == "decay") rho else NA_real_,
    d0 = d0,
    statistic = obs$statistic,
    expectation = obs$expectation,
    variance = variance,
    z_score = pv$z,
    p_value = pv$p,
    stringsAsFactors = FALSE
  )

  out <- list(
    result = result,
    detail = obs$detail,
    covariance = Sigma,
    set_distance = set_dist,
    match = match
  )
  class(out) <- "netmatch_test"
  out
}

#' @export
print.netmatch_test <- function(x, ...) {
  print(x$result, row.names = FALSE)
  invisible(x)
}

#' Sensitivity Grid for Model-Assisted Network Dependence
#'
#' Evaluates `netmatch_test(method = "decay")` over a grid of eta and rho
#' values.
#'
#' @param match A `netmatch` object.
#' @param outcome Name of the outcome column.
#' @param eta Numeric vector of eta values.
#' @param rho Numeric vector of rho values.
#' @param d0 Maximum set distance with nonzero across-set covariance.
#' @return A `netmatch_sensitivity` object.
#' @export
netmatch_sensitivity <- function(match,
                                 outcome,
                                 eta = seq(0, 0.10, by = 0.01),
                                 rho = seq(0, 0.50, by = 0.05),
                                 d0 = 2) {
  grid <- expand.grid(eta = eta, rho = rho, KEEP.OUT.ATTRS = FALSE)
  rows <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    fit <- netmatch_test(
      match = match,
      outcome = outcome,
      method = "decay",
      eta = grid$eta[i],
      rho = grid$rho[i],
      d0 = d0
    )
    rows[[i]] <- fit$result
  }
  out <- list(grid = do.call(rbind, rows), match = match, outcome = outcome)
  class(out) <- "netmatch_sensitivity"
  out
}

#' @export
print.netmatch_sensitivity <- function(x, ...) {
  print(utils::head(x$grid), row.names = FALSE)
  cat("...\n")
  invisible(x)
}
