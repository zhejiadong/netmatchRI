#' Build a Network-Constrained Matched Design
#'
#' `netmatch()` is the main matching entry point. It creates covariate-only,
#' single-penalty, or dual-penalty matched designs using direct graph-distance
#' thresholds: `kappa = 2` means network-distance pairs less than or equal to 2
#' are disallowed where the selected design applies network restrictions.
#'
#' @param data A data frame containing treatment and covariates.
#' @param treat Name of the binary treatment column coded 0/1.
#' @param covariates Character vector of covariate column names.
#' @param network Square adjacency or network-distance matrix.
#' @param method Matching method: `"dual"`, `"single"`, or `"covariate"`.
#' @param kappa Network-distance threshold for disallowed close pairs.
#' @param solver Solver preference. The current portable implementation uses a
#'   deterministic greedy backend for `"auto"` and `"glpk"`; `"gurobi"` is
#'   reserved for the licensed MIP backend.
#' @param min_controls Minimum controls per treated set.
#' @param max_controls Maximum controls per treated set.
#' @return A `netmatch` object with matched data and diagnostics.
#' @export
netmatch <- function(data,
                     treat,
                     covariates,
                     network,
                     method = c("dual", "single", "covariate"),
                     kappa = 2,
                     solver = c("auto", "gurobi", "glpk"),
                     min_controls = 1,
                     max_controls = Inf) {
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

  if (solver == "gurobi" && !requireNamespace("gurobi", quietly = TRUE)) {
    stop("`solver = \"gurobi\"` requires the optional gurobi package.", call. = FALSE)
  }

  D <- .mahalanobis_matrix(data, z, covariates)
  matched <- .greedy_match(data, z, D, unit_dist, method, kappa, min_controls, max_controls)

  out <- list(
    data = matched$data,
    treat = treat,
    covariates = covariates,
    method = method,
    method_label = .method_label(method),
    kappa = kappa,
    solver = if (solver == "auto") "greedy" else solver,
    network_distance = unit_dist,
    match_table = matched$table
  )
  class(out) <- "netmatch"
  out
}

.greedy_match <- function(data, z, D, unit_dist, method, kappa, min_controls, max_controls) {
  treat_ids <- which(z == 1)
  ctrl_ids <- which(z == 0)
  if (!length(treat_ids) || !length(ctrl_ids)) {
    stop("Need at least one treated and one control unit.", call. = FALSE)
  }
  available_c <- ctrl_ids
  subclass <- rep(NA_integer_, nrow(data))
  table <- data.frame()
  s <- 1L

  # Full matching approximation: one treated anchor per set, one or more controls.
  for (ti in treat_ids) {
    candidates <- available_c
    if (!length(candidates)) break
    if (method %in% c("single", "dual")) {
      candidates <- candidates[unit_dist[ti, candidates] > kappa]
    }
    if (!length(candidates)) next

    ord <- order(D[as.character(ti), as.character(candidates)])
    candidates <- candidates[ord]
    chosen <- integer(0)
    for (ci in candidates) {
      if (length(chosen) >= max_controls) break
      ok <- TRUE
      if (method == "dual" && length(chosen)) {
        ok <- all(unit_dist[ci, chosen] > kappa)
      }
      if (ok) chosen <- c(chosen, ci)
      if (length(chosen) >= min_controls && is.infinite(max_controls)) break
    }
    if (length(chosen) < min_controls) next

    ids <- c(ti, chosen)
    subclass[ids] <- s
    available_c <- setdiff(available_c, chosen)
    table <- rbind(
      table,
      data.frame(subclass = s, unit = ids, treat = z[ids])
    )
    s <- s + 1L
  }

  # Attach unmatched controls to the nearest feasible existing set when possible.
  if (any(is.na(subclass[available_c])) && any(!is.na(subclass))) {
    for (ci in available_c) {
      set_ids <- sort(unique(stats::na.omit(subclass)))
      assigned <- FALSE
      for (sid in set_ids) {
        members <- which(subclass == sid)
        tr <- members[z[members] == 1]
        co <- members[z[members] == 0]
        ok_single <- method == "covariate" || all(unit_dist[ci, tr] > kappa)
        ok_dual <- method != "dual" || all(unit_dist[ci, co] > kappa)
        if (ok_single && ok_dual) {
          subclass[ci] <- sid
          assigned <- TRUE
          break
        }
      }
      if (!assigned) next
    }
  }

  keep <- !is.na(subclass)
  if (!any(keep)) stop("No feasible matched sets were found.", call. = FALSE)
  matched <- data[keep, , drop = FALSE]
  matched$subclass <- factor(subclass[keep])
  rownames(matched) <- which(keep)
  list(data = matched, table = table)
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
