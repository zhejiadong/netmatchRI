.validate_outcome <- function(data, outcome, subclass = "subclass") {
  if (!outcome %in% names(data)) {
    stop("`outcome` column not found.", call. = FALSE)
  }
  y <- data[[outcome]]
  if (!is.numeric(y) || anyNA(y) || any(!is.finite(y))) {
    stop("`outcome` must be numeric, finite and non-missing.", call. = FALSE)
  }
  invisible(TRUE)
}

.mw_u <- function(y_t, y_c) {
  sum(outer(y_t, y_c, ">"))
}

.observed_stats <- function(data,
                            outcome,
                            treat = "Z",
                            subclass = "subclass",
                            weight_type = c("ns", "ntc")) {
  .validate_outcome(data, outcome, subclass)
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

.empty_set_distance <- function(sets) {
  sets <- sort(unique(sets))
  matrix(0, length(sets), length(sets), dimnames = list(sets, sets))
}

.unadjusted_covariance_matrix <- function(stats_df) {
  Sigma <- matrix(0, nrow(stats_df), nrow(stats_df), dimnames = list(stats_df$subclass, stats_df$subclass))
  diag(Sigma) <- stats_df$var
  Sigma
}

.design_cov_bound <- function(nt_s, nc_s, nt_k, nc_k) {
  counts <- c(nt_s, nc_s, nt_k, nc_k)
  if (length(counts) != 4L || any(!is.finite(counts)) ||
      any(counts < 0 | counts != floor(counts))) {
    stop("Matched-set counts must be finite non-negative integers.", call. = FALSE)
  }
  # Use doubles for support products, including integer-valued input counts.
  counts <- as.double(counts)
  a <- sort(counts[1:2]); b <- sort(counts[3:4])
  if (a[1] == 0 || b[1] == 0) return(0)
  if (identical(a, b)) return(prod(a) * (sum(a) + 1) / 12)

  # The finite-support sum of rectangular differences of min(F_s, F_k)
  # equals the integral of centered quantile products on their common CDF
  # intervals. Full matching has uniform U marginals; integer endpoints
  # avoid numerical CDF evaluations and midpoint rounding at jumps.
  if (a[1] == 1 && b[1] == 1) {
    ns <- sum(a); nk <- sum(b)
    ends <- sort(unique(c((0:ns) * nk, (0:nk) * ns)))
    left <- ends[-length(ends)]
    value <- sum((diff(ends) / (ns * nk)) *
                   (floor(left / nk) - (ns - 1) / 2) *
                   (floor(left / ns) - (nk - 1) / 2))
  } else {
    # Wilcoxon symmetry lets us double the integral on [0, 1/2]. This
    # avoids CDFs rounded to one in the upper tail. Evaluate immediately
    # to the right of each left endpoint, not at a rounded midpoint.
    us <- prod(a); uk <- prod(b)
    fs <- stats::pwilcox(0:floor(us / 2), a[1], a[2])
    fk <- stats::pwilcox(0:floor(uk / 2), b[1], b[2])
    fs <- fs[fs < 0.5]; fk <- fk[fk < 0.5]
    ends <- sort(unique(c(0, fs, fk, 0.5)))
    left <- ends[-length(ends)]
    qs <- findInterval(left, fs)
    qk <- findInterval(left, fk)
    value <- 2 * sum(diff(ends) * (qs - us / 2) * (qk - uk / 2))
  }
  max(0, value)
}

.set_distance_matrix <- function(unit_dist, subclass_vec, unit_ids = seq_along(subclass_vec)) {
  if (!is.matrix(unit_dist) || !is.numeric(unit_dist) ||
      nrow(unit_dist) != ncol(unit_dist)) {
    stop("`network_distance` must be a square numeric distance matrix for network inference.", call. = FALSE)
  }
  if (length(unit_ids) != length(subclass_vec) || anyNA(unit_ids) ||
      anyDuplicated(unit_ids) || any(unit_ids < 1 | unit_ids > nrow(unit_dist)) ||
      any(unit_ids != floor(unit_ids))) {
    stop("Matched unit row indices must uniquely index `network_distance`.", call. = FALSE)
  }
  # Validate matched-unit distances only on network inference paths, never RI_unadjusted.
  .validate_network(unit_dist[unit_ids, unit_ids, drop = FALSE], "distance")
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
  cov_cache <- new.env(parent = emptyenv())
  for (i in seq_len(S - 1)) {
    for (j in (i + 1):S) {
      if (!is.finite(set_dist[i, j]) || set_dist[i, j] > kappa) next
      key_i <- paste(stats_df$nt[i], stats_df$nc[i], sep = ":")
      key_j <- paste(stats_df$nt[j], stats_df$nc[j], sep = ":")
      key <- paste(sort(c(key_i, key_j)), collapse = "|")
      if (exists(key, envir = cov_cache, inherits = FALSE)) {
        val <- get(key, envir = cov_cache, inherits = FALSE)
      } else {
        val <- .design_cov_bound(stats_df$nt[i], stats_df$nc[i], stats_df$nt[j], stats_df$nc[j])
        assign(key, val, envir = cov_cache)
      }
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
  if (S < 2 || method == "unadjusted") {
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

.weighted_variance_from_components <- function(components, weights, method, eta = 1, rho = 1) {
  bound <- components$bound
  diag_var <- sum((weights^2) * diag(bound))
  S <- nrow(bound)
  if (S < 2 || method == "unadjusted") return(diag_var)
  pair_idx <- which(upper.tri(bound) & bound != 0, arr.ind = TRUE)
  if (!nrow(pair_idx)) return(diag_var)
  off_diag <- 2 * weights[pair_idx[, 1]] * weights[pair_idx[, 2]] * bound[pair_idx]
  if (method == "adjusted") {
    set_dist <- components$set_dist[pair_idx]
    off_diag <- off_diag * eta * rho^(set_dist - 1)
  }
  diag_var + sum(off_diag)
}

.normal_pvalue <- function(statistic, expectation, variance) {
  if (!is.finite(variance) || variance <= 0) {
    return(list(z = NA_real_, p = NA_real_))
  }
  z <- (statistic - expectation) / sqrt(variance)
  list(z = z, p = 2 * stats::pnorm(abs(z), lower.tail = FALSE))
}
