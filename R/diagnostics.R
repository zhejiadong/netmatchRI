#' Diagnose a Matched Design
#'
#' Computes compact covariate and network diagnostics for a `netmatch` object.
#'
#' @param match A `netmatch` object.
#' @return A list with `covariate_balance`, `network_distance`,
#'   `within_distance_table`, and backward-compatible aliases
#'   `covariate_smd` and `average_within_distance`.
#' @examples
#' \dontrun{
#' sim <- simulate_netmatch_example()
#' m <- netmatch(sim$data, "Z", c("X1", "X2", "X3"), sim$net_dist,
#'               method = "dual", kappa = 2, solver = "auto")
#' diag <- diagnose_match(m)
#' diag$covariate_balance
#' diag$network_distance
#' diag$within_distance_table
#' }
#' @export
diagnose_match <- function(match) {
  if (!inherits(match, "netmatch")) {
    stop("`match` must be a netmatch object.", call. = FALSE)
  }
  df <- match$data
  z <- df[[match$treat]]
  cov_smd <- data.frame(
    covariate = match$covariates,
    abs_smd = vapply(match$covariates, function(x) .matched_smd(df, x, match$treat), numeric(1)),
    stringsAsFactors = FALSE
  )

  within_d <- .within_distances(df, match$network_distance)
  tab <- table(factor(within_d, levels = sort(unique(within_d))))
  dist_tab <- data.frame(
    distance = as.numeric(names(tab)),
    count = as.integer(tab),
    proportion = as.numeric(tab) / sum(tab)
  )

  network_summary <- data.frame(
    n_pairs = length(within_d),
    min_distance = min(within_d, na.rm = TRUE),
    mean_distance = mean(within_d, na.rm = TRUE),
    max_distance = max(within_d, na.rm = TRUE)
  )

  list(
    covariate_balance = cov_smd,
    network_distance = network_summary,
    within_distance_table = dist_tab,
    covariate_smd = cov_smd,
    average_within_distance = network_summary$mean_distance
  )
}

.matched_smd <- function(df, x, treat) {
  z <- df[[treat]]
  s <- df$subclass
  sd0 <- sqrt((stats::var(df[[x]][z == 1]) + stats::var(df[[x]][z == 0])) / 2)
  if (!is.finite(sd0) || sd0 == 0) return(NA_real_)
  sets <- split(seq_len(nrow(df)), s)
  diffs <- vapply(sets, function(ids) {
    zz <- z[ids]
    if (!any(zz == 1) || !any(zz == 0)) return(NA_real_)
    mean(df[[x]][ids][zz == 1]) - mean(df[[x]][ids][zz == 0])
  }, numeric(1))
  abs(mean(diffs, na.rm = TRUE) / sd0)
}

.within_distances <- function(df, dist) {
  idx <- as.integer(rownames(df))
  sets <- split(idx, df$subclass)
  out <- unlist(lapply(sets, function(ids) {
    if (length(ids) < 2) return(numeric(0))
    pr <- utils::combn(ids, 2)
    dist[cbind(pr[1, ], pr[2, ])]
  }), use.names = FALSE)
  out[is.finite(out)]
}
