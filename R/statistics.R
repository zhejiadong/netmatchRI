.mw_u <- function(y_t, y_c) {
  sum(outer(y_t, y_c, ">"))
}

.observed_stats <- function(data,
                            outcome,
                            treat = "Z",
                            subclass = "subclass",
                            weight_type = c("ns", "ntc")) {
  weight_type <- match.arg(weight_type)
  sets <- sort(unique(data[[subclass]]))
  rows <- lapply(sets, function(s) {
    sub <- data[data[[subclass]] == s, , drop = FALSE]
    z <- sub[[treat]]
    y <- sub[[outcome]]
    nt <- sum(z == 1)
    nc <- sum(z == 0)
    n <- nt + nc
    u <- if (nt > 0 && nc > 0) .mw_u(y[z == 1], y[z == 0]) else 0
    data.frame(
      subclass = s,
      n = n,
      nt = nt,
      nc = nc,
      U = u,
      mu = nt * nc / 2,
      var = nt * nc * (n + 1) / 12,
      weight = switch(weight_type,
        ns = 1 / (n + 1),
        ntc = if (nt > 0 && nc > 0) 1 / (nt * nc) else 0
      )
    )
  })
  stats_df <- do.call(rbind, rows)
  list(
    detail = stats_df,
    statistic = sum(stats_df$weight * stats_df$U),
    expectation = sum(stats_df$weight * stats_df$mu)
  )
}

.design_cov_bound <- function(nt_s, nc_s, nt_k, nc_k) {
  grid <- seq(0.0001, 0.9999, length.out = 10000)
  q_s <- stats::qwilcox(grid, m = nt_s, n = nc_s)
  q_k <- stats::qwilcox(grid, m = nt_k, n = nc_k)
  mean(q_s * q_k) - mean(q_s) * mean(q_k)
}

.set_distance_matrix <- function(unit_dist, subclass_vec, unit_ids = seq_along(subclass_vec)) {
  sets <- sort(unique(subclass_vec))
  S <- length(sets)
  out <- matrix(0, S, S, dimnames = list(sets, sets))
  ids <- lapply(sets, function(s) unit_ids[which(subclass_vec == s)])
  if (S < 2) return(out)
  for (i in seq_len(S - 1)) {
    for (j in (i + 1):S) {
      dij <- min(unit_dist[ids[[i]], ids[[j]], drop = FALSE], na.rm = TRUE)
      if (!is.finite(dij)) dij <- Inf
      out[i, j] <- out[j, i] <- dij
    }
  }
  out
}

.variance_components <- function(stats_df, set_dist, kappa) {
  S <- nrow(stats_df)
  bound <- matrix(0, S, S)
  diag(bound) <- stats_df$var
  if (S < 2) return(list(bound = bound, set_dist = set_dist, kappa = kappa))
  for (i in seq_len(S - 1)) {
    for (j in (i + 1):S) {
      val <- .design_cov_bound(stats_df$nt[i], stats_df$nc[i], stats_df$nt[j], stats_df$nc[j])
      if (!is.finite(set_dist[i, j]) || set_dist[i, j] > kappa) val <- 0
      bound[i, j] <- bound[j, i] <- val
    }
  }
  list(bound = bound, set_dist = set_dist, kappa = kappa)
}

.covariance_matrix <- function(stats_df, set_dist, method, eta = 1, rho = 1, kappa = Inf) {
  components <- .variance_components(stats_df, set_dist, kappa)
  .covariance_from_components(components, method, eta, rho)
}

.covariance_from_components <- function(components, method, eta = 1, rho = 1) {
  Sigma <- components$bound
  S <- nrow(Sigma)
  diag_vals <- diag(Sigma)
  if (S < 2 || method == "naive") {
    Sigma[,] <- 0
    diag(Sigma) <- diag_vals
    return(Sigma)
  }
  if (method == "design") return(Sigma)
  set_dist <- components$set_dist
  for (i in seq_len(S - 1)) {
    for (j in (i + 1):S) {
      if (Sigma[i, j] == 0) next
      Sigma[i, j] <- Sigma[j, i] <- Sigma[i, j] * eta * rho^(set_dist[i, j] - 1)
    }
  }
  diag(Sigma) <- diag_vals
  Sigma
}

.normal_pvalue <- function(statistic, expectation, variance) {
  if (!is.finite(variance) || variance <= 0) {
    return(list(z = NA_real_, p = NA_real_))
  }
  z <- (statistic - expectation) / sqrt(variance)
  list(z = z, p = 2 * (1 - stats::pnorm(abs(z))))
}
