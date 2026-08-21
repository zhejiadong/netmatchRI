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
#'   `Y` when randomization-based inference will be run later.
#' @param treat Name of the binary treatment column coded 0/1.
#' @param covariates Character vector of covariate column names.
#' @param network Square adjacency or network-distance matrix.
#' @param method Matching method: `"dual"`, `"single"`, or `"covariate"`.
#' @param kappa Network-distance threshold for disallowed close pairs.
#' @param solver Solver preference for `method = "dual"`. `"highs"` is the
#'   default open-source backend. `"auto"` tries Gurobi, HiGHS, then GLPK.
#'   `"gurobi"`, `"highs"`, and `"glpk"` request one backend explicitly. The
#'   `"covariate"` and `"single"` methods use `optmatch::fullmatch()`.
#' @param min_controls Minimum controls per treated unit in each matched set.
#' @param max_controls Maximum controls per treated unit in each matched set.
#' @param caliper Mahalanobis-distance caliper for feasible treated-control
#'   edges.
#' @param timelimit Solver time limit in seconds.
#' @param mipgap Relative MIP gap for Gurobi and HiGHS.
#' @param threads Number of Gurobi or HiGHS threads.
#' @param estimand Target estimand used to construct unit matching weights:
#'   `"ATT"`, `"ATC"`, or `"ATE"`. Nonzero weights are normalized within
#'   each treatment group to have mean 1.
#' @param network_type Interpretation of `network`: `"adjacency"`,
#'   `"distance"`, or `"auto"`. An explicit choice overrides automatic
#'   detection. Infinite values are allowed only for distance matrices.
#' @param include_solver If `TRUE`, retain the backend's raw result in
#'   `solver_result`. The default keeps only stable `solver_info`.
#' @return A `netmatch` object containing matched `data`, full-length
#'   `subclass`, `weights`, and `matched` vectors, the original call and data,
#'   the estimand and resolved network type, the unit-level network-distance
#'   matrix, and stable solver information. Raw solver output is included only
#'   when `include_solver = TRUE`.
#' @details `solver`, `caliper`, `timelimit`, `mipgap`, and `threads` are used
#'   only by `method = "dual"`; comparison methods use
#'   `optmatch::fullmatch()`. `kappa` is ignored by `method = "covariate"`.
#'   `method = "single"` uses `kappa` but not the mixed-integer solver controls.
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
#'   solver = "highs"
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
                     solver = c("highs", "auto", "gurobi", "glpk"),
                     min_controls = 0.01,
                     max_controls = 100,
                     caliper = 8,
                     timelimit = 90,
                     mipgap = 0.01,
                     threads = 1,
                     estimand = c("ATT", "ATC", "ATE"),
                     network_type = c("auto", "adjacency", "distance"),
                     include_solver = FALSE) {
  matched_call <- match.call()
  method <- match.arg(method)
  solver <- match.arg(solver)
  estimand <- match.arg(estimand)
  network_type <- match.arg(network_type)
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  if (!treat %in% names(data)) stop("`treat` column not found.", call. = FALSE)
  .validate_covariates(data, covariates)
  z <- data[[treat]]
  .check_binary(z, treat)
  unit_dist <- .as_network_distance(network, network_type, rownames(data))
  resolved_network_type <- attr(unit_dist, "network_type")
  attr(unit_dist, "network_type") <- NULL
  if (nrow(unit_dist) != nrow(data)) {
    stop("`network` size must match `nrow(data)`.", call. = FALSE)
  }
  if (!is.logical(include_solver) || length(include_solver) != 1L || is.na(include_solver)) {
    stop("`include_solver` must be TRUE or FALSE.", call. = FALSE)
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
    solver_candidates <- .resolve_solver_candidates(solver)
    matched <- .solve_network_match(
      data = data,
      z = z,
      D = D,
      unit_dist = unit_dist,
      method = method,
      kappa = kappa,
      solver = solver_candidates,
      min_controls = min_controls,
      max_controls = max_controls,
      caliper = caliper,
      timelimit = timelimit,
      mipgap = mipgap,
      threads = threads
    )
    solver_used <- matched$solver
  }

  subclass <- rep(NA_integer_, nrow(data))
  matched_ids <- as.integer(rownames(matched$data))
  subclass[matched_ids] <- as.integer(matched$data$subclass)
  subclass <- factor(subclass, levels = sort(unique(subclass[!is.na(subclass)])))
  matched_units <- !is.na(subclass)
  weights <- .matching_weights(z, subclass, estimand)
  solver_info <- .solver_info(solver_used, matched$solver_result)

  matched_data_frame <- data[matched_units, , drop = FALSE]
  matched_data_frame$subclass <- droplevels(subclass[matched_units])
  rownames(matched_data_frame) <- which(matched_units)

  out <- list(
    data = matched_data_frame,
    treat = treat,
    covariates = covariates,
    method = method,
    kappa = kappa,
    solver = solver_used,
    network_distance = unit_dist,
    subclass = subclass,
    weights = weights,
    matched = matched_units,
    call = matched_call,
    estimand = estimand,
    network_type = resolved_network_type,
    original_data = data,
    solver_info = solver_info
  )
  if (include_solver) out$solver_result <- matched$solver_result
  class(out) <- "netmatch"
  out
}

.matching_weights <- function(z, subclass, estimand) {
  weights <- numeric(length(z))
  matched <- !is.na(subclass)
  sets <- split(which(matched), subclass[matched], drop = TRUE)
  for (ids in sets) {
    nt <- sum(z[ids] == 1)
    nc <- sum(z[ids] == 0)
    if (nt == 0 || nc == 0) next
    if (estimand == "ATT") {
      weights[ids[z[ids] == 1]] <- 1
      weights[ids[z[ids] == 0]] <- nt / nc
    } else if (estimand == "ATC") {
      weights[ids[z[ids] == 1]] <- nc / nt
      weights[ids[z[ids] == 0]] <- 1
    } else {
      weights[ids[z[ids] == 1]] <- (nt + nc) / nt
      weights[ids[z[ids] == 0]] <- (nt + nc) / nc
    }
  }
  for (group in c(0, 1)) {
    ids <- which(matched & z == group & weights > 0)
    if (length(ids)) weights[ids] <- weights[ids] / mean(weights[ids])
  }
  weights
}

.scalar_or_na <- function(x, character = FALSE) {
  if (is.null(x) || !length(x)) return(if (character) NA_character_ else NA_real_)
  if (character) as.character(x[[1L]]) else as.numeric(x[[1L]])
}

.solver_info <- function(backend, result) {
  if (!is.list(result)) {
    return(list(
      backend = backend,
      status = if (identical(backend, "optmatch")) "completed" else NA_character_,
      objective = NA_real_, gap = NA_real_, runtime = NA_real_
    ))
  }
  status <- if (identical(backend, "optmatch")) "completed" else .scalar_or_na(result$status, TRUE)
  objective <- .scalar_or_na(if (!is.null(result$objective)) result$objective else result$objval)
  gap <- .scalar_or_na(if (!is.null(result$metadata$gap)) result$metadata$gap else result$mipgap)
  runtime <- .scalar_or_na(
    if (!is.null(result$metadata$run_time)) result$metadata$run_time else result$runtime
  )
  list(
    backend = backend,
    status = status,
    objective = objective,
    gap = gap,
    runtime = runtime
  )
}

#' Extract Data from a Matched Design
#'
#' Returns the original unit-level data with the full-length matched-set
#' membership and matching weights appended.
#'
#' @param object A `netmatch` object.
#' @param drop_unmatched If `TRUE`, omit unmatched units. If `FALSE`, retain
#'   them with missing `subclass` and zero `weights`.
#' @return A data frame containing the original columns plus `subclass` and
#'   `weights`.
#' @export
matched_data <- function(object, drop_unmatched = TRUE) {
  if (!inherits(object, "netmatch")) {
    stop("`object` must be a netmatch object.", call. = FALSE)
  }
  if (!is.logical(drop_unmatched) || length(drop_unmatched) != 1L || is.na(drop_unmatched)) {
    stop("`drop_unmatched` must be TRUE or FALSE.", call. = FALSE)
  }
  out <- object$original_data
  out$subclass <- object$subclass
  out$weights <- object$weights
  if (drop_unmatched) out <- out[object$matched, , drop = FALSE]
  out
}

.resolve_solver <- function(solver) {
  .resolve_solver_candidates(solver)[1L]
}

.resolve_solver_candidates <- function(solver) {
  if (solver == "gurobi") {
    if (.gurobi_available()) return("gurobi")
    stop("Gurobi is not available or its license test failed.", call. = FALSE)
  }
  if (solver == "highs") {
    if (.highs_available()) return("highs")
    stop("Package `highs` is required for `solver = \"highs\"`.", call. = FALSE)
  }
  if (solver == "glpk") {
    if (.glpk_available()) return("glpk")
    stop("Rglpk and slam are required for `solver = \"glpk\"`.", call. = FALSE)
  }
  if (solver == "auto") {
    candidates <- c(
      if (.gurobi_available()) "gurobi",
      if (.highs_available()) "highs",
      if (.glpk_available()) "glpk"
    )
    if (length(candidates)) return(candidates)
    stop("Dual matching requires Gurobi, HiGHS, or Rglpk.", call. = FALSE)
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

.highs_available <- function() {
  requireNamespace("highs", quietly = TRUE)
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
  rownames(D) <- rownames(data)[treat_ids]
  colnames(D) <- rownames(data)[ctrl_ids]
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
  list(data = matched, solver_result = fm)
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
  solve_once <- function(candidate) {
    .gurobi_match_once(
      data = data,
      z = z,
      D = D,
      unit_dist = unit_dist,
      method = method,
      kappa = kappa,
      solver = candidate,
      min_controls = min_controls,
      max_controls = max_controls,
      caliper = caliper,
      timelimit = timelimit,
      mipgap = mipgap,
      threads = threads
    )
  }

  if (length(solver) == 1L) {
    matched <- solve_once(solver)
    matched$solver <- solver
    return(matched)
  }

  errors <- character(0)
  for (candidate in solver) {
    matched <- tryCatch(
      solve_once(candidate),
      error = function(e) {
        errors[[candidate]] <<- conditionMessage(e)
        NULL
      }
    )
    if (!is.null(matched)) {
      matched$solver <- candidate
      return(matched)
    }
  }
  details <- paste(paste0(names(errors), ": ", unname(errors)), collapse = "; ")
  stop("No solver in the `auto` chain returned a usable matched design. ", details,
       call. = FALSE)
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
  } else if (solver == "highs") {
    res <- .solve_highs(
      model,
      timelimit = timelimit,
      mipgap = mipgap,
      threads = threads
    )
  } else {
    params <- list(
      TimeLimit = timelimit,
      MIPGap = mipgap,
      Threads = as.integer(max(1, threads)),
      OutputFlag = 0
    )
    res <- gurobi::gurobi(model, params = params)
  }
  if (solver == "glpk" && res$status %in% c("INFEASIBLE", "INF_OR_UNBD", "UNBOUNDED")) {
    stop(
      "GLPK could not produce a matched design within the current search limit. ",
      "Increase `timelimit` or try `solver = \"gurobi\"`. ",
      "Last GLPK status: `", res$status, "`.",
      call. = FALSE
    )
  }
  if (res$status %in% c("INFEASIBLE", "INF_OR_UNBD", "UNBOUNDED")) {
    stop("MIP model is infeasible or unbounded.", call. = FALSE)
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

  reported_objective <- if (identical(solver, "gurobi")) res$objval else NULL
  checked <- .validate_mip_solution(model, res$x, reported_objective = reported_objective)
  res$x <- checked$x
  res$objective <- checked$objective

  subclass <- .extract_subclasses(res, treat_ids, ctrl_ids, E, nrow(data))
  keep_units <- !is.na(subclass)
  if (!any(keep_units)) stop("No feasible matched sets were found.", call. = FALSE)
  matched <- data[keep_units, , drop = FALSE]
  matched$subclass <- factor(subclass[keep_units])
  rownames(matched) <- which(keep_units)
  list(data = matched, solver_result = res)
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

.highs_row_bounds <- function(sense, rhs) {
  if (length(sense) != length(rhs)) {
    stop("Model `sense` and `rhs` must have the same length.", call. = FALSE)
  }
  known <- sense %in% c(">=", "<=", "=")
  if (any(!known)) {
    stop(
      "Unknown model sense: ", paste(unique(sense[!known]), collapse = ", "),
      call. = FALSE
    )
  }
  lhs <- rep(-Inf, length(rhs))
  upper_rhs <- rep(Inf, length(rhs))
  lhs[sense == ">="] <- rhs[sense == ">="]
  upper_rhs[sense == "<="] <- rhs[sense == "<="]
  lhs[sense == "="] <- rhs[sense == "="]
  upper_rhs[sense == "="] <- rhs[sense == "="]
  list(lhs = lhs, rhs = upper_rhs)
}

.validate_mip_solution <- function(model,
                                   x,
                                   reported_objective = NULL,
                                   tolerance = 1e-6) {
  p <- length(model$obj)
  if (length(x) != p) {
    stop("Solver solution has the wrong length.", call. = FALSE)
  }
  if (!all(is.finite(x))) {
    stop("Solver solution contains non-finite values.", call. = FALSE)
  }
  if (any(x < -tolerance | x > 1 + tolerance)) {
    stop("Solver solution violates binary variable bounds.", call. = FALSE)
  }
  if (any(abs(x - round(x)) > tolerance)) {
    stop("Solver solution violates integrality.", call. = FALSE)
  }
  x <- as.numeric(round(x))
  if (length(model$sense) != length(model$rhs) || nrow(model$A) != length(model$rhs)) {
    stop("Model constraint dimensions are inconsistent.", call. = FALSE)
  }
  known <- model$sense %in% c(">=", "<=", "=")
  if (any(!known)) {
    stop(
      "Unknown model sense: ", paste(unique(model$sense[!known]), collapse = ", "),
      call. = FALSE
    )
  }
  activity <- as.numeric(model$A %*% x)
  violated <- (model$sense == ">=" & activity < model$rhs - tolerance) |
    (model$sense == "<=" & activity > model$rhs + tolerance) |
    (model$sense == "=" & abs(activity - model$rhs) > tolerance)
  if (any(violated)) {
    stop("Solver solution violates model constraints.", call. = FALSE)
  }
  objective <- sum(model$obj * x)
  if (!is.null(reported_objective)) {
    if (length(reported_objective) != 1L || !is.finite(reported_objective) ||
        abs(objective - reported_objective) > tolerance * max(1, abs(objective))) {
      stop("Solver objective does not match the recomputed objective.", call. = FALSE)
    }
  }
  list(x = as.numeric(round(x)), objective = objective, activity = activity)
}

.normalize_highs_result <- function(sol, model, run_time) {
  status_code <- as.integer(sol$status)[1L]
  status_message <- as.character(sol$status_message)[1L]
  info <- sol$info
  gap <- if (!is.null(info$mip_gap)) as.numeric(info$mip_gap)[1L] else NA_real_
  metadata <- list(
    status_code = status_code,
    status_message = status_message,
    gap = gap,
    run_time = as.numeric(run_time)[1L],
    time_limited_incumbent = FALSE
  )
  if (identical(status_code, 8L)) {
    return(list(status = "INFEASIBLE", x = NULL, raw = sol, metadata = metadata))
  }
  if (identical(status_code, 9L)) {
    return(list(status = "INF_OR_UNBD", x = NULL, raw = sol, metadata = metadata))
  }
  if (identical(status_code, 10L)) {
    return(list(status = "UNBOUNDED", x = NULL, raw = sol, metadata = metadata))
  }
  if (is.na(status_code) || !status_code %in% c(7L, 13L)) {
    return(list(status = "UNKNOWN", x = NULL, raw = sol, metadata = metadata))
  }
  if (!identical(info$primal_solution_status, "Feasible")) {
    status <- if (identical(status_code, 13L)) {
      "TIME_LIMIT_NO_INCUMBENT"
    } else {
      "OPTIMAL_NO_SOLUTION"
    }
    return(list(status = status, x = NULL, raw = sol, metadata = metadata))
  }
  checked <- .validate_mip_solution(
    model,
    sol$primal_solution,
    reported_objective = sol$objective_value
  )
  if (identical(status_code, 13L)) {
    metadata$time_limited_incumbent <- TRUE
    warning(
      "HiGHS reached its time limit with a validated feasible incumbent; status: ",
      status_message, "; actual relative MIP gap: ",
      if (is.finite(gap)) format(gap, digits = 6) else "unavailable", ".",
      call. = FALSE
    )
  }
  list(
    status = if (identical(status_code, 7L)) "OPTIMAL" else "TIME_LIMIT",
    x = checked$x,
    objective = checked$objective,
    raw = sol,
    metadata = metadata
  )
}

.solve_highs <- function(model, timelimit, mipgap, threads) {
  bounds <- .highs_row_bounds(model$sense, model$rhs)
  p <- length(model$obj)
  started <- proc.time()[["elapsed"]]
  sol <- highs::highs_solve(
    L = model$obj,
    lower = rep.int(0, p),
    upper = rep.int(1, p),
    A = model$A,
    lhs = bounds$lhs,
    rhs = bounds$rhs,
    types = rep.int("I", p),
    maximum = FALSE,
    control = list(
      output_flag = FALSE,
      presolve = "on",
      time_limit = as.numeric(timelimit),
      mip_rel_gap = as.numeric(mipgap),
      threads = as.integer(max(1, threads))
    )
  )
  run_time <- proc.time()[["elapsed"]] - started
  .normalize_highs_result(sol, model, run_time = run_time)
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
  selected_local <- sort(unique(c(sel$t, length(treat_ids) + sel$c)))
  local_membership <- vapply(selected_local, find, integer(1))
  remap <- match(local_membership, unique(local_membership))
  global_ids <- c(treat_ids, ctrl_ids)[selected_local]
  subclass <- rep(NA_integer_, n)
  subclass[global_ids] <- remap
  subclass
}

#' @export
print.netmatch <- function(x, ...) {
  n_sets <- length(unique(x$subclass[!is.na(x$subclass)]))
  cat("<netmatch> ", x$method, " matching; ", x$estimand, " estimand\n", sep = "")
  cat("  ", sum(x$matched), "/", length(x$matched), " units in ", n_sets,
      " matched sets; ", x$solver_info$backend, " backend\n", sep = "")
  invisible(x)
}

#' Summarize a Network-Matched Design
#'
#' @param object A `netmatch` object.
#' @param ... Unused.
#' @return A `netmatch_summary` object containing sample counts, effective
#'   sample sizes, matched-set sizes, covariate balance, network diagnostics,
#'   and stable solver information.
#' @rdname netmatch
#' @export
summary.netmatch <- function(object, ...) {
  z <- object$original_data[[object$treat]]
  groups <- list(overall = rep(TRUE, length(z)), treated = z == 1, control = z == 0)
  sample_counts <- do.call(rbind, lapply(names(groups), function(group) {
    ids <- groups[[group]]
    data.frame(
      group = group,
      original = sum(ids),
      matched = sum(ids & object$matched),
      unmatched = sum(ids & !object$matched),
      stringsAsFactors = FALSE
    )
  }))
  rownames(sample_counts) <- NULL

  ess <- do.call(rbind, lapply(names(groups), function(group) {
    w <- object$weights[groups[[group]]]
    data.frame(
      group = group,
      ess = if (sum(w^2) > 0) sum(w)^2 / sum(w^2) else 0,
      stringsAsFactors = FALSE
    )
  }))
  rownames(ess) <- NULL

  df <- object$data
  matched_z <- df[[object$treat]]
  set_tab <- stats::aggregate(matched_z, list(subclass = df$subclass), function(v) {
    c(n = length(v), nt = sum(v == 1), nc = sum(v == 0))
  })
  set_sizes <- do.call(data.frame, set_tab)
  names(set_sizes) <- c("subclass", "n", "nt", "nc")
  diagnostics <- diagnose_match(object)

  out <- list(
    sample_counts = sample_counts,
    ess = ess,
    set_sizes = set_sizes,
    covariate_balance = diagnostics$covariate_balance,
    network_summary = diagnostics$network_summary,
    solver_info = object$solver_info
  )
  class(out) <- "netmatch_summary"
  out
}

#' @export
print.netmatch_summary <- function(x, ...) {
  cat("Matched design summary\n")
  cat("\nSample counts\n")
  print(x$sample_counts, row.names = FALSE)
  cat("\nEffective sample size\n")
  print(x$ess, row.names = FALSE)
  cat("\nMatched-set sizes\n")
  print(x$set_sizes, row.names = FALSE)
  cat("\nCovariate balance (absolute SMD)\n")
  print(x$covariate_balance, row.names = FALSE)
  cat("\nNetwork summary\n")
  print(x$network_summary, row.names = FALSE)
  cat("\nSolver: ", x$solver_info$backend, " (", x$solver_info$status, ")\n", sep = "")
  invisible(x)
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
                    main = x$method)
  invisible(diag)
}
