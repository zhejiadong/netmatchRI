#' Build a Dual-Penalty Matched Design
#'
#' `netmatch()` builds the proposed dual-penalty matched design via a mixed-
#' integer program. The covariate-only and single-penalty comparison designs
#' use `optmatch::fullmatch()`. The cutoff is a direct graph-distance
#' threshold: `kappa = 2` means network-distance pairs less than or equal to 2
#' are disallowed where the selected design applies network restrictions.
#'
#' @param data A data frame with one row per unit. Users supplying their own
#'   data should include at least a binary treatment column such as `Z`,
#'   covariate columns such as `X1`, `X2`, `X3`, and an outcome column such as
#'   `Y` when randomization inference will be run later.
#' @param treat Name of the binary treatment column coded 0/1.
#' @param covariates Character vector of covariate column names.
#' @param network Square adjacency or network-distance matrix.
#' @param method Matching method: `"dual"`, `"single"`, or `"covariate"`.
#' @param kappa Network-distance threshold for disallowed close pairs.
#' @param solver Solver preference for `method = "dual"`. `"auto"` tries
#'   Gurobi first, then the open-source GLPK backend through `Rglpk`.
#'   `"gurobi"` uses Gurobi only. `"glpk"` uses `Rglpk` with a sparse triplet
#'   constraint matrix. The `"covariate"` and `"single"` methods use
#'   `optmatch::fullmatch()`.
#' @param min_controls Minimum controls per treated unit in each matched set.
#' @param max_controls Maximum controls per treated unit in each matched set.
#' @param caliper Mahalanobis-distance caliper for feasible treated-control
#'   edges.
#' @param timelimit Solver time limit in seconds.
#' @param mipgap Gurobi MIP gap.
#' @param threads Number of Gurobi threads.
#' @return A `netmatch` object, which is a list with matched `data`, treatment
#'   and covariate names, `method`, `method_label`, `kappa`, solver backend,
#'   the graph-distance matrix in `network_distance`, a unit-level
#'   `match_table`, and the raw `solver_result`.
#' @examples
#' \dontrun{
#' sim <- simulate_netmatch_example()
#' m <- netmatch(
#'   data = sim$data,
#'   treat = "Z",
#'   covariates = c("X1", "X2", "X3"),
#'   network = sim$net_dist,
#'   method = "dual",
#'   kappa = 2,
#'   solver = "auto"
#' )
#' m
#' summary(m)
#' }
#' @export
netmatch <- function(data,
                     treat,
                     covariates,
                     network,
                     method = c("dual", "single", "covariate"),
                     kappa = 2,
                     solver = c("auto", "gurobi", "glpk"),
                     min_controls = 0.01,
                     max_controls = 100,
                     caliper = 8,
                     timelimit = 50,
                     mipgap = 0.01,
                     threads = 1) {
  method <- match.arg(method)
  solver <- match.arg(solver)
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  if (!treat %in% names(data)) stop("`treat` column not found.", call. = FALSE)
  missing_cov <- setdiff(covariates, names(data))
  if (length(missing_cov)) {
    stop("Missing covariates: ", paste(missing_cov, collapse = ", "), call. = FALSE)
  }
  z <- data[[treat]]
  .check_binary(z, treat)
  unit_dist <- .as_network_distance(network)
  if (nrow(unit_dist) != nrow(data)) {
    stop("`network` size must match `nrow(data)`.", call. = FALSE)
  }
  if (!is.numeric(kappa) || length(kappa) != 1 || !is.finite(kappa) || kappa < 0) {
    stop("`kappa` must be one non-negative graph-distance threshold.", call. = FALSE)
  }

  D <- .mahalanobis_matrix(data, z, covariates, cov_type = "pooled")
  if (method %in% c("covariate", "single")) {
    matched <- .optmatch_match(
      data = data,
      z = z,
      D = D,
      unit_dist = unit_dist,
      method = method,
      kappa = kappa,
      min_controls = min_controls,
      max_controls = max_controls
    )
    solver_used <- "optmatch"
  } else {
    solver_used <- .resolve_solver(solver)
    matched <- .solve_network_match(
      data = data,
      z = z,
      D = D,
      unit_dist = unit_dist,
      method = method,
      kappa = kappa,
      solver = solver_used,
      min_controls = min_controls,
      max_controls = max_controls,
      caliper = caliper,
      timelimit = timelimit,
      mipgap = mipgap,
      threads = threads
    )
  }

  out <- list(
    data = matched$data,
    treat = treat,
    covariates = covariates,
    method = method,
    method_label = .method_label(method),
    kappa = kappa,
    solver = solver_used,
    network_distance = unit_dist,
    match_table = matched$table,
    solver_result = matched$solver_result
  )
  class(out) <- "netmatch"
  out
}

.resolve_solver <- function(solver) {
  if (solver == "gurobi") {
    if (.gurobi_available()) return("gurobi")
    stop("Gurobi is not available or its license test failed.", call. = FALSE)
  }
  if (solver == "glpk") {
    if (.glpk_available()) return("glpk")
    stop("Rglpk and slam are required for `solver = \"glpk\"`.", call. = FALSE)
  }
  if (solver == "auto") {
    if (.gurobi_available()) return("gurobi")
    if (.glpk_available()) return("glpk")
    stop("Dual matching requires Gurobi or Rglpk.", call. = FALSE)
  }
  stop("Unknown solver.", call. = FALSE)
}

.gurobi_available <- function() {
  if (!requireNamespace("gurobi", quietly = TRUE)) return(FALSE)
  tryCatch({
    model <- list(
      A = Matrix::Matrix(c(1, 1), nrow = 1, sparse = TRUE),
      obj = c(1, 2),
      sense = ">=",
      rhs = 1,
      vtype = c("B", "B"),
      modelsense = "min"
    )
    res <- gurobi::gurobi(model, params = list(OutputFlag = 0))
    identical(res$status, "OPTIMAL")
  }, error = function(e) FALSE)
}

.glpk_available <- function() {
  requireNamespace("Rglpk", quietly = TRUE) &&
    requireNamespace("slam", quietly = TRUE)
}

.optmatch_match <- function(data,
                            z,
                            D,
                            unit_dist,
                            method,
                            kappa,
                            min_controls,
                            max_controls) {
  if (!requireNamespace("optmatch", quietly = TRUE)) {
    stop("`optmatch` is required for covariate-only and single-penalty matching.", call. = FALSE)
  }
  treat_ids <- which(z == 1)
  ctrl_ids <- which(z == 0)
  if (!length(treat_ids) || !length(ctrl_ids)) {
    stop("Need at least one treated and one control unit.", call. = FALSE)
  }
  D <- as.matrix(D)
  rownames(D) <- as.character(treat_ids)
  colnames(D) <- as.character(ctrl_ids)
  if (method == "single") {
    tc_dist <- unit_dist[treat_ids, ctrl_ids, drop = FALSE]
    D[tc_dist <= kappa] <- Inf
  }
  fm <- optmatch::fullmatch(
    D,
    data = data,
    min.controls = min_controls,
    max.controls = max_controls
  )
  subclass <- rep(NA_integer_, nrow(data))
  keep <- !is.na(fm)
  subclass[keep] <- as.integer(factor(fm[keep]))
  if (!any(keep)) stop("No feasible matched sets were found.", call. = FALSE)
  matched <- data[keep, , drop = FALSE]
  matched$subclass <- factor(subclass[keep])
  rownames(matched) <- which(keep)
  tab <- data.frame(
    subclass = as.integer(matched$subclass),
    unit = as.integer(rownames(matched)),
    treat = z[as.integer(rownames(matched))]
  )
  list(data = matched, table = tab, solver_result = fm)
}

.solve_network_match <- function(data,
                                 z,
                                 D,
                                 unit_dist,
                                 method,
                                 kappa,
                                 solver,
                                 min_controls,
                                 max_controls,
                                 caliper,
                                 timelimit,
                                 mipgap,
                                 threads) {
  last_error <- NULL
  res <- tryCatch(
    .gurobi_match_once(
      data = data,
      z = z,
      D = D,
      unit_dist = unit_dist,
      method = method,
      kappa = kappa,
      solver = solver,
      min_controls = min_controls,
      max_controls = max_controls,
      caliper = caliper,
      timelimit = timelimit,
      mipgap = mipgap,
      threads = threads
    ),
    error = function(e) {
      last_error <<- conditionMessage(e)
      NULL
    }
  )
  if (!is.null(res)) return(res)
  if (identical(solver, "glpk") &&
      grepl("GLPK|usable matched design within the current search limit", last_error, fixed = FALSE)) {
    stop(last_error, call. = FALSE)
  }
  stop("No feasible matched design was found. Last solver message: ", last_error, call. = FALSE)
}

.gurobi_match_once <- function(data,
                               z,
                               D,
                               unit_dist,
                               method,
                               kappa,
                               solver,
                               min_controls,
                               max_controls,
                               caliper,
                               timelimit,
                               mipgap,
                               threads) {
  treat_ids <- which(z == 1)
  ctrl_ids <- which(z == 0)
  if (!length(treat_ids) || !length(ctrl_ids)) {
    stop("Need at least one treated and one control unit.", call. = FALSE)
  }

  D <- as.matrix(D)
  rownames(D) <- as.character(treat_ids)
  colnames(D) <- as.character(ctrl_ids)
  tc_dist <- unit_dist[treat_ids, ctrl_ids, drop = FALSE]
  keep <- is.finite(D) & D < caliper
  if (method %in% c("single", "dual")) {
    keep <- keep & tc_dist > kappa
  }
  if (!any(keep)) stop("No feasible treated-control edges under caliper and kappa.", call. = FALSE)
  if (any(rowSums(keep) == 0)) {
    stop("Infeasible: at least one treated unit has no feasible control.", call. = FALSE)
  }
  if (min_controls != 0 && any(colSums(keep) == 0)) {
    stop("Infeasible: at least one control unit has no feasible treated unit.", call. = FALSE)
  }

  which_keep <- which(keep, arr.ind = TRUE)
  E <- data.frame(
    e = seq_len(nrow(which_keep)),
    t = which_keep[, 1],
    c = which_keep[, 2]
  )
  E$cost <- D[cbind(E$t, E$c)]
  eid <- matrix(0L, nrow = length(treat_ids), ncol = length(ctrl_ids))
  eid[cbind(E$t, E$c)] <- E$e

  minC_per_T <- max(1, min_controls)
  maxC_per_T <- if (is.finite(max_controls)) max_controls else length(ctrl_ids)
  if (min_controls == 0) {
    minT_per_C <- 0
    maxT_per_C <- length(treat_ids)
  } else {
    minT_per_C <- 1
    maxT_per_C <- floor(1 / min_controls)
    if (!is.finite(maxT_per_C)) maxT_per_C <- length(treat_ids)
    maxT_per_C <- max(1, maxT_per_C)
  }

  cc_pairs_list <- vector("list", length(treat_ids))
  tt_pairs_list <- vector("list", length(ctrl_ids))
  if (method == "dual") {
    ctrl_dist <- unit_dist[ctrl_ids, ctrl_ids, drop = FALSE]
    treat_dist <- unit_dist[treat_ids, treat_ids, drop = FALSE]
    for (t in seq_along(treat_ids)) {
      cc_pairs_list[[t]] <- .close_pairs(keep[t, ], ctrl_dist, kappa)
    }
    for (c in seq_along(ctrl_ids)) {
      tt_pairs_list[[c]] <- .close_pairs(keep[, c], treat_dist, kappa)
    }
  }

  model <- .build_gurobi_model(
    E = E,
    eid = eid,
    nT = length(treat_ids),
    nC = length(ctrl_ids),
    minC_per_T = minC_per_T,
    maxC_per_T = maxC_per_T,
    minT_per_C = minT_per_C,
    maxT_per_C = maxT_per_C,
    cc_pairs_list = cc_pairs_list,
    tt_pairs_list = tt_pairs_list
  )
  if (solver == "glpk") {
    res <- .solve_roi_glpk(model, timelimit = timelimit)
  } else {
    res <- gurobi::gurobi(
      model,
      params = list(
        TimeLimit = timelimit,
        MIPGap = mipgap,
        Threads = as.integer(max(1, threads)),
        OutputFlag = 0
      )
    )
  }
  if (solver == "glpk" && res$status %in% c("INFEASIBLE", "INF_OR_UNBD")) {
    stop(
      "GLPK could not produce a matched design within the current search limit. ",
      "Increase `timelimit` or try `solver = \"gurobi\"`. ",
      "Last GLPK status: `", res$status, "`.",
      call. = FALSE
    )
  }
  if (res$status %in% c("INFEASIBLE", "INF_OR_UNBD")) {
    stop("MIP model is infeasible.", call. = FALSE)
  }
  if (solver == "glpk" && (!res$status %in% c("OPTIMAL", "TIME_LIMIT") || is.null(res$x))) {
    stop(
      "GLPK did not return a usable matched design within the current search limit. ",
      "Increase `timelimit` or try `solver = \"gurobi\"`. ",
      "Last GLPK status: `", res$status, "`.",
      call. = FALSE
    )
  }
  if (!res$status %in% c("OPTIMAL", "TIME_LIMIT") || is.null(res$x)) {
    stop("MIP solver returned status `", res$status, "` without a usable solution.", call. = FALSE)
  }

  subclass <- .extract_subclasses(res, treat_ids, ctrl_ids, E, nrow(data))
  keep_units <- !is.na(subclass)
  if (!any(keep_units)) stop("No feasible matched sets were found.", call. = FALSE)
  matched <- data[keep_units, , drop = FALSE]
  matched$subclass <- factor(subclass[keep_units])
  rownames(matched) <- which(keep_units)
  tab <- data.frame(
    subclass = as.integer(matched$subclass),
    unit = as.integer(rownames(matched)),
    treat = z[as.integer(rownames(matched))]
  )
  list(data = matched, table = tab, solver_result = res)
}

.solve_roi_glpk <- function(model, timelimit) {
  if (!requireNamespace("slam", quietly = TRUE)) {
    stop("GLPK backend requires package `slam` to keep the MIP constraint matrix sparse.", call. = FALSE)
  }
  sol <- Rglpk::Rglpk_solve_LP(
    obj = model$obj,
    mat = slam::as.simple_triplet_matrix(model$A),
    dir = model$sense,
    rhs = model$rhs,
    types = model$vtype,
    max = FALSE,
    control = list(
      presolve = TRUE,
      tm_limit = as.integer(max(1, timelimit) * 1000),
      canonicalize_status = FALSE
    )
  )
  status <- as.character(sol$status)
  x <- as.numeric(sol$solution)
  if (identical(status, "5") || identical(status, "0")) status <- "OPTIMAL"
  list(status = status, x = x, raw = sol)
}

.close_pairs <- function(keep_slice, dist_mat, kappa) {
  cand <- which(keep_slice)
  if (length(cand) < 2) return(data.frame(c1 = integer(0), c2 = integer(0)))
  d <- dist_mat[cand, cand, drop = FALSE]
  diag(d) <- Inf
  close <- which(d <= kappa & upper.tri(d), arr.ind = TRUE)
  if (!nrow(close)) return(data.frame(c1 = integer(0), c2 = integer(0)))
  data.frame(c1 = cand[close[, "row"]], c2 = cand[close[, "col"]])
}

.build_gurobi_model <- function(E,
                                eid,
                                nT,
                                nC,
                                minC_per_T,
                                maxC_per_T,
                                minT_per_C,
                                maxT_per_C,
                                cc_pairs_list,
                                tt_pairs_list) {
  m <- nrow(E)
  Ai <- integer(0)
  Aj <- integer(0)
  Ax <- numeric(0)
  sense <- character(0)
  rhs <- numeric(0)

  add_row <- function(cols, coeffs, s, r) {
    row_idx <- length(rhs) + 1L
    Ai <<- c(Ai, rep(row_idx, length(cols)))
    Aj <<- c(Aj, cols)
    Ax <<- c(Ax, coeffs)
    sense <<- c(sense, s)
    rhs <<- c(rhs, r)
  }

  for (t in seq_len(nT)) {
    cols <- E$e[E$t == t]
    if (length(cols)) {
      if (is.finite(minC_per_T) && minC_per_T > 0) add_row(cols, rep(1, length(cols)), ">=", minC_per_T)
      if (is.finite(maxC_per_T)) add_row(cols, rep(1, length(cols)), "<=", maxC_per_T)
    }
  }
  for (c in seq_len(nC)) {
    cols <- E$e[E$c == c]
    if (length(cols)) {
      if (is.finite(minT_per_C) && minT_per_C > 0) add_row(cols, rep(1, length(cols)), ">=", minT_per_C)
      if (is.finite(maxT_per_C)) add_row(cols, rep(1, length(cols)), "<=", maxT_per_C)
    }
  }

  cc_pairs <- lapply(seq_len(nT), function(t) {
    cp <- cc_pairs_list[[t]]
    if (is.null(cp) || !nrow(cp)) return(NULL)
    e1 <- eid[t, cp$c1]
    e2 <- eid[t, cp$c2]
    keep <- e1 > 0L & e2 > 0L
    if (!any(keep)) return(NULL)
    cbind(e1[keep], e2[keep])
  })
  cc_pairs <- do.call(rbind, cc_pairs[!vapply(cc_pairs, is.null, logical(1))])
  if (!is.null(cc_pairs) && nrow(cc_pairs)) {
    k <- nrow(cc_pairs)
    rows <- (length(rhs) + 1L):(length(rhs) + k)
    Ai <- c(Ai, rep(rows, each = 2L))
    Aj <- c(Aj, as.vector(t(cc_pairs)))
    Ax <- c(Ax, rep(1, 2L * k))
    sense <- c(sense, rep("<=", k))
    rhs <- c(rhs, rep(1, k))
  }

  tt_pairs <- lapply(seq_len(nC), function(c) {
    tp <- tt_pairs_list[[c]]
    if (is.null(tp) || !nrow(tp)) return(NULL)
    e1 <- eid[tp$c1, c]
    e2 <- eid[tp$c2, c]
    keep <- e1 > 0L & e2 > 0L
    if (!any(keep)) return(NULL)
    cbind(e1[keep], e2[keep])
  })
  tt_pairs <- do.call(rbind, tt_pairs[!vapply(tt_pairs, is.null, logical(1))])
  if (!is.null(tt_pairs) && nrow(tt_pairs)) {
    k <- nrow(tt_pairs)
    rows <- (length(rhs) + 1L):(length(rhs) + k)
    Ai <- c(Ai, rep(rows, each = 2L))
    Aj <- c(Aj, as.vector(t(tt_pairs)))
    Ax <- c(Ax, rep(1, 2L * k))
    sense <- c(sense, rep("<=", k))
    rhs <- c(rhs, rep(1, k))
  }

  A <- if (length(Ai)) {
    Matrix::sparseMatrix(i = Ai, j = Aj, x = Ax, dims = c(length(rhs), m))
  } else {
    Matrix::Matrix(0, nrow = 0, ncol = m, sparse = TRUE)
  }
  list(A = A, obj = E$cost, sense = sense, rhs = rhs, vtype = rep("B", m), modelsense = "min")
}

.extract_subclasses <- function(res, treat_ids, ctrl_ids, E, n) {
  chosen <- which(res$x > 0.5)
  if (!length(chosen)) return(rep(NA_integer_, n))
  sel <- E[chosen, , drop = FALSE]
  parent <- seq_len(length(treat_ids) + length(ctrl_ids))
  find <- function(x) {
    while (parent[x] != x) {
      parent[x] <<- parent[parent[x]]
      x <- parent[x]
    }
    x
  }
  union <- function(a, b) {
    ra <- find(a)
    rb <- find(b)
    if (ra != rb) parent[rb] <<- ra
  }
  for (i in seq_len(nrow(sel))) {
    union(sel$t[i], length(treat_ids) + sel$c[i])
  }
  local_membership <- vapply(seq_along(parent), find, integer(1))
  remap <- match(local_membership, unique(local_membership))
  subclass <- rep(NA_integer_, n)
  subclass[treat_ids] <- remap[seq_along(treat_ids)]
  subclass[ctrl_ids] <- remap[length(treat_ids) + seq_along(ctrl_ids)]
  subclass
}

#' @export
print.netmatch <- function(x, ...) {
  cat("<netmatch>\n")
  cat("  Method: ", x$method_label, "\n", sep = "")
  cat("  Kappa: ", x$kappa, "\n", sep = "")
  cat("  Solver: ", x$solver, "\n", sep = "")
  cat("  Matched units: ", nrow(x$data), "\n", sep = "")
  cat("  Matched sets: ", length(unique(x$data$subclass)), "\n", sep = "")
  invisible(x)
}

#' Summarize a Network-Matched Design
#'
#' @param object A `netmatch` object.
#' @param ... Unused.
#' @return A data frame with one row per matched set and columns for subclass,
#'   total set size, treated count, and control count.
#' @rdname netmatch
#' @export
summary.netmatch <- function(object, ...) {
  df <- object$data
  z <- df[[object$treat]]
  set_tab <- stats::aggregate(z, list(subclass = df$subclass), function(v) {
    c(n = length(v), nt = sum(v == 1), nc = sum(v == 0))
  })
  stats <- do.call(data.frame, set_tab)
  names(stats) <- c("subclass", "n", "nt", "nc")
  stats
}

#' Plot Within-Set Network Distances
#'
#' @param x A `netmatch` object.
#' @param ... Unused.
#' @return Invisibly returns the diagnostics list from `diagnose_match()`.
#' @rdname netmatch
#' @export
plot.netmatch <- function(x, ...) {
  diag <- diagnose_match(x)
  graphics::barplot(diag$within_distance_table$proportion,
                    names.arg = diag$within_distance_table$distance,
                    xlab = "Network distance",
                    ylab = "Within-set proportion",
                    main = x$method_label)
  invisible(diag)
}
