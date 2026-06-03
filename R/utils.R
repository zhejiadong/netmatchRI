.check_binary <- function(x, name = "treat") {
  if (!all(stats::na.omit(x) %in% c(0, 1))) {
    stop(sprintf("`%s` must be coded 0/1.", name), call. = FALSE)
  }
  invisible(TRUE)
}

.as_network_distance <- function(network) {
  if (is.null(network)) stop("`network` is required.", call. = FALSE)
  network <- as.matrix(network)
  if (nrow(network) != ncol(network)) {
    stop("`network` must be a square adjacency or distance matrix.", call. = FALSE)
  }
  diag(network) <- 0

  finite_vals <- network[is.finite(network) & upper.tri(network)]
  finite_vals <- finite_vals[finite_vals > 0]
  is_adjacency <- length(finite_vals) > 0 &&
    all(finite_vals %in% c(1)) &&
    all(network %in% c(0, 1), na.rm = TRUE)

  if (!is_adjacency) {
    return(network)
  }

  .floyd_warshall(network)
}

.floyd_warshall <- function(adj) {
  n <- nrow(adj)
  d <- matrix(Inf, n, n)
  d[adj > 0] <- 1
  diag(d) <- 0
  for (k in seq_len(n)) {
    d <- pmin(d, outer(d[, k], d[k, ], "+"))
  }
  d
}

.safe_inverse <- function(S) {
  out <- tryCatch(solve(S), error = function(e) NULL)
  if (!is.null(out)) return(out)
  sv <- svd(S)
  keep <- sv$d > sqrt(.Machine$double.eps) * max(sv$d)
  if (!any(keep)) stop("Covariance matrix is singular.", call. = FALSE)
  sv$v[, keep, drop = FALSE] %*%
    diag(1 / sv$d[keep], nrow = sum(keep)) %*%
    t(sv$u[, keep, drop = FALSE])
}

.mahalanobis_matrix <- function(data, treat, covariates) {
  X <- stats::model.matrix(
    stats::reformulate(covariates),
    data = data
  )
  X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  X <- scale(X)
  S_inv <- .safe_inverse(stats::cov(X))
  Q <- X %*% S_inv %*% t(X)
  d2 <- outer(diag(Q), diag(Q), "+") - 2 * Q
  d2[d2 < 0] <- 0
  D <- sqrt(d2)
  rownames(D) <- colnames(D) <- seq_len(nrow(data))
  ids_t <- which(treat == 1)
  ids_c <- which(treat == 0)
  D[ids_t, ids_c, drop = FALSE]
}

.method_label <- function(method) {
  switch(method,
    covariate = "Covariates only",
    single = "Single penalty",
    dual = "Dual penalty",
    method
  )
}
