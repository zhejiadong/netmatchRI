#' Diagnose a Matched Design
#'
#' Checks covariate balance and network distance within matched sets for a
#' `netmatch` object. Each level of a factor covariate is reported as an
#' indicator. For each reported covariate, both standardized mean differences use
#' the same pooled original-sample standard deviation,
#' `sqrt(((n1 - 1) * var1 + (n0 - 1) * var0) / (n1 + n0 - 2))`. The before
#' difference uses unweighted original-sample means; the after difference uses
#' the unit matching weights among matched units.
#'
#' @param match A \code{netmatch} object.
#' @return A list containing covariate-balance results, a network-distance
#'   summary, and the within-matched-set distance table.
#' @examples
#' \dontrun{
#' sim <- simulate_netmatch_example()
#' m <- netmatch(sim$data, "Z", c("X1", "X2", "X3"), sim$net_dist,
#'               method = "dual", kappa = 2, solver = "auto")
#' diag <- diagnose_match(m)
#' diag$covariate_balance
#' diag$network_summary
#' diag$within_distance_table
#' }
#' @export
diagnose_match <- function(match) {
  if (!inherits(match, "netmatch")) {
    stop("`match` must be a netmatch object.", call. = FALSE)
  }
  covariate_balance <- .covariate_balance(match)
  within_d <- .within_distances(match$subclass, match$network_distance)
  finite_d <- within_d[is.finite(within_d)]

  if (length(within_d)) {
    values <- sort(unique(within_d))
    counts <- vapply(values, function(value) sum(within_d == value), integer(1))
    dist_tab <- data.frame(
      distance = values,
      count = counts,
      proportion = counts / length(within_d)
    )
  } else {
    dist_tab <- data.frame(
      distance = numeric(0), count = integer(0), proportion = numeric(0)
    )
  }

  network_summary <- data.frame(
    n_pairs = length(within_d),
    n_finite_pairs = length(finite_d),
    n_disconnected_pairs = sum(is.infinite(within_d)),
    min_distance = if (length(finite_d)) min(finite_d) else NA_real_,
    mean_distance = if (length(finite_d)) mean(finite_d) else NA_real_,
    max_distance = if (length(finite_d)) max(finite_d) else NA_real_
  )

  list(
    covariate_balance = covariate_balance,
    network_summary = network_summary,
    within_distance_table = dist_tab
  )
}

.covariate_balance <- function(match) {
  X <- .balance_covariate_matrix(match$original_data, match$covariates)
  z <- match$original_data[[match$treat]]
  w <- match$weights
  matched <- match$matched & w > 0

  rows <- lapply(seq_len(ncol(X)), function(j) {
    x <- X[, j]
    n1 <- sum(z == 1)
    n0 <- sum(z == 0)
    pooled_sd <- sqrt(
      ((n1 - 1) * stats::var(x[z == 1]) + (n0 - 1) * stats::var(x[z == 0])) /
        (n1 + n0 - 2)
    )
    before <- mean(x[z == 1]) - mean(x[z == 0])
    after <- .weighted_mean(x[matched & z == 1], w[matched & z == 1]) -
      .weighted_mean(x[matched & z == 0], w[matched & z == 0])
    if (!is.finite(pooled_sd) || pooled_sd == 0) {
      before_smd <- after_smd <- NA_real_
    } else {
      before_smd <- abs(before / pooled_sd)
      after_smd <- abs(after / pooled_sd)
    }
    data.frame(
      covariate = colnames(X)[j],
      before_abs_smd = before_smd,
      after_abs_smd = after_smd,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

.weighted_mean <- function(x, w) {
  if (!length(x) || !length(w) || sum(w) <= 0) return(NA_real_)
  sum(x * w) / sum(w)
}

.within_distances <- function(subclass, dist) {
  ids <- which(!is.na(subclass))
  sets <- split(ids, subclass[ids], drop = TRUE)
  unlist(lapply(sets, function(set_ids) {
    if (length(set_ids) < 2) return(numeric(0))
    pairs <- utils::combn(set_ids, 2)
    dist[cbind(pairs[1, ], pairs[2, ])]
  }), use.names = FALSE)
}
