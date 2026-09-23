.check_binary <- function(x, name = "treat") {
  if (anyNA(x)) {
    stop(sprintf("`%s` must not contain missing values.", name), call. = FALSE)
  }
  if (!is.numeric(x) && !is.integer(x) && !is.logical(x)) {
    stop(sprintf("`%s` must be coded 0/1.", name), call. = FALSE)
  }
  if (!all(x %in% c(0, 1))) {
    stop(sprintf("`%s` must be coded 0/1.", name), call. = FALSE)
  }
  invisible(TRUE)
}

.validate_covariates <- function(data, covariates) {
  missing_cov <- setdiff(covariates, names(data))
  if (length(missing_cov)) {
    stop("Missing covariates: ", paste(missing_cov, collapse = ", "), call. = FALSE)
  }
  has_missing <- vapply(data[covariates], anyNA, logical(1))
  if (any(has_missing)) {
    stop(
      "Matching covariates must not contain missing values: ",
      paste(covariates[has_missing], collapse = ", "), ".",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.validate_reserved_columns <- function(data) {
  reserved <- intersect(c("subclass", "weights"), names(data))
  if (length(reserved)) {
    stop(
      "`data` must not contain reserved columns: ",
      paste(reserved, collapse = ", "), ". Rename them before matching.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.validate_network <- function(network,
                              network_type = c("auto", "adjacency", "distance"),
                              data_rownames = NULL) {
  network_type <- match.arg(network_type)
  if (is.null(network)) stop("`network` is required.", call. = FALSE)
  if (!is.matrix(network) || !is.numeric(network)) {
    stop("`network` must be a numeric matrix.", call. = FALSE)
  }
  if (nrow(network) != ncol(network)) {
    stop("`network` must be square.", call. = FALSE)
  }
  if (anyNA(network) || any(is.nan(network))) {
    stop("`network` must not contain missing or NaN values.", call. = FALSE)
  }

  rn <- rownames(network)
  cn <- colnames(network)
  if (!is.null(rn) || !is.null(cn)) {
    if (is.null(rn) || is.null(cn)) {
      stop("A named `network` must have both row and column names.", call. = FALSE)
    }
    if (anyDuplicated(rn) || anyDuplicated(cn) || !setequal(rn, cn)) {
      stop("`network` row and column names must identify the same unique units.", call. = FALSE)
    }
    if (!is.null(data_rownames)) {
      if (anyDuplicated(data_rownames) || !setequal(rn, data_rownames)) {
        stop("`network` dimnames must match `data` row names.", call. = FALSE)
      }
      network <- network[data_rownames, data_rownames, drop = FALSE]
    }
  }

  if (!isTRUE(all.equal(network, t(network), tolerance = sqrt(.Machine$double.eps),
                        check.attributes = FALSE))) {
    stop("`network` must be symmetric.", call. = FALSE)
  }
  if (any(!is.finite(diag(network))) || any(diag(network) != 0)) {
    stop("`network` must have a zero diagonal.", call. = FALSE)
  }
  if (any(network < 0, na.rm = TRUE)) {
    stop("`network` entries must be nonnegative.", call. = FALSE)
  }

  resolved_type <- network_type
  if (network_type == "auto") {
    resolved_type <- if (all(is.finite(network)) && all(network %in% c(0, 1))) {
      "adjacency"
    } else {
      "distance"
    }
  }
  if (resolved_type == "adjacency") {
    if (any(!is.finite(network))) {
      stop("`network` cannot contain Inf when `network_type = \"adjacency\"`.", call. = FALSE)
    }
    if (!all(network %in% c(0, 1))) {
      stop("An adjacency `network` must contain only 0 and 1.", call. = FALSE)
    }
  }
  if (resolved_type == "distance" && any(network[row(network) != col(network)] < 1)) {
    stop("A distance `network` must have off-diagonal distances >= 1 (or Inf for disconnected units).", call. = FALSE)
  }
  list(network = network, network_type = resolved_type)
}

.as_network_distance <- function(network,
                                 network_type = c("auto", "adjacency", "distance"),
                                 data_rownames = NULL) {
  checked <- .validate_network(network, network_type, data_rownames)
  network <- checked$network
  if (checked$network_type == "distance") {
    attr(network, "network_type") <- checked$network_type
    return(network)
  }
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package `igraph` is required to convert adjacency matrices to network distances.", call. = FALSE)
  }
  g <- igraph::graph_from_adjacency_matrix(network, mode = "undirected", diag = FALSE)
  out <- igraph::distances(g)
  attr(out, "network_type") <- checked$network_type
  out
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

.ginv <- function(S) {
  if (requireNamespace("MASS", quietly = TRUE)) return(MASS::ginv(S))
  .safe_inverse(S)
}

.covariate_matrix <- function(data, covariates) {
  X <- stats::model.matrix(stats::reformulate(covariates), data = data)
  X[, colnames(X) != "(Intercept)", drop = FALSE]
}

.balance_covariate_matrix <- function(data, covariates) {
  pieces <- lapply(covariates, function(name) {
    x <- data[[name]]
    if (is.factor(x) || is.character(x)) {
      x <- factor(x, levels = if (is.factor(x)) levels(x) else unique(x))
      out <- vapply(levels(x), function(level) as.numeric(x == level), numeric(length(x)))
      if (is.null(dim(out))) out <- matrix(out, ncol = 1L)
      colnames(out) <- paste0(name, "=", levels(x))
      return(out)
    }
    out <- matrix(as.numeric(x), ncol = 1L)
    colnames(out) <- name
    out
  })
  out <- do.call(cbind, pieces)
  colnames(out) <- make.unique(colnames(out))
  out
}

.mahalanobis_matrix <- function(data, treat, covariates, cov_type = c("pooled", "overall")) {
  cov_type <- match.arg(cov_type)
  X <- .covariate_matrix(data, covariates)
  keep <- apply(X, 2, stats::sd) > 0
  if (!any(keep)) {
    stop("At least one matching covariate must have nonzero variance.", call. = FALSE)
  }
  X <- scale(X[, keep, drop = FALSE])
  if (cov_type == "pooled") {
    if (sum(treat == 1) < 2 || sum(treat == 0) < 2) {
      stop("Pooled Mahalanobis distance requires at least two treated and two control units.", call. = FALSE)
    }
    S1 <- stats::cov(X[treat == 1, , drop = FALSE])
    S0 <- stats::cov(X[treat == 0, , drop = FALSE])
    S <- ((sum(treat == 1) - 1) * S1 + (sum(treat == 0) - 1) * S0) /
      (length(treat) - 2)
  } else {
    S <- stats::cov(X)
  }
  S_inv <- .ginv(S)
  Q <- X %*% S_inv %*% t(X)
  d2 <- outer(diag(Q), diag(Q), "+") - 2 * Q
  d2[d2 < 0] <- 0
  D <- sqrt(d2)
  rownames(D) <- colnames(D) <- seq_len(nrow(data))
  ids_t <- which(treat == 1)
  ids_c <- which(treat == 0)
  D[ids_t, ids_c, drop = FALSE]
}
