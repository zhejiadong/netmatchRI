#' Randomization-Based Inference Given a Matched Design
#'
#' `RI_Naive()`, `RI_Sensitivity()`, and `RI_Design()` run two-sided,
#' normal-approximation randomization-based inference for a `netmatch` object.
#' `RI_Sensitivity()` allows residual dependence between matched sets, while
#' `RI_Design()` provides the design-based analysis and keeps the distance
#' truncation at `kappa`: across-set covariance is set to zero whenever the
#' matched-set distance is greater than `kappa`.
#'
#' @param match A `netmatch` object.
#' @param outcome Name of the outcome column.
#' @param method Internal variance method: `"sensitivity"`, `"naive"`, or
#'   `"design"`.
#' @param eta A number from 0 to 1 that scales residual covariance between
#'   matched sets.
#' @param rho A number from 0 to 1 that controls how residual covariance changes
#'   with network distance. For sets at distance `d`, the sensitivity multiplier
#'   is `eta * rho^(d - 1)`.
#' @param kappa Analysis cutoff. Defaults to the matching cutoff stored in
#'   `match$kappa`.
#' @param weight_type Matched-set weighting. `"ns"` uses `1 / (n_s + 1)`,
#'   where `n_s` is the set size. `"ntc"` uses `1 / (n_ts * n_cs)`, where
#'   `n_ts` and `n_cs` are the treated and control counts in the set.
#' @return A randomization-based inference result object with a one-row `result` data frame,
#'   matched-set `detail`, covariance matrix, matched-set distance matrix,
#'   outcome, weighting choice, complete analysis options, analysis `kappa`,
#'   and the original `match`.
#' @details For two matched sets at distance `d`, sensitivity analysis multiplies
#'   the covariance bound by `eta * rho^(d - 1)`. The set statistic counts a
#'   treated outcome only when it is strictly greater than a control outcome;
#'   tied outcomes add zero. P-values use a two-sided normal approximation.
#' @examples
#' \dontrun{
#' sim <- simulate_netmatch_example()
#' m <- netmatch(sim$data, "Z", c("X1", "X2", "X3"), sim$net_dist,
#'               method = "dual", kappa = 2, solver = "auto")
#' RI_Naive(m, "Y")
#' RI_Sensitivity(m, "Y", eta = 0.03, rho = 0.10)
#' RI_Design(m, "Y")
#' }
#' @name RI_Naive
NULL

netmatch_test <- function(match,
                          outcome,
                          method = c("sensitivity", "naive", "design"),
                          eta = 0.03,
                          rho = 0.10,
                          kappa = NULL,
                          weight_type = c("ns", "ntc")) {
  if (!inherits(match, "netmatch")) {
    stop("`match` must be a netmatch object.", call. = FALSE)
  }
  method <- match.arg(method)
  weight_type <- match.arg(weight_type)
  if (method == "sensitivity") {
    .validate_unit_interval(eta, "eta", scalar = TRUE)
    .validate_unit_interval(rho, "rho", scalar = TRUE)
  }
  if (!outcome %in% names(match$data)) stop("`outcome` column not found.", call. = FALSE)
  kappa <- .analysis_kappa(match, kappa)

  obs <- .observed_stats(match$data, outcome, match$treat, "subclass", weight_type)
  w <- obs$detail$weight
  if (method == "naive") {
    set_dist <- .empty_set_distance(obs$detail$subclass)
    Sigma <- .naive_covariance_matrix(obs$detail)
    variance <- sum((w^2) * obs$detail$var)
  } else {
    unit_ids <- as.integer(rownames(match$data))
    set_dist <- .set_distance_matrix(match$network_distance, as.integer(match$data$subclass), unit_ids)
    Sigma <- .covariance_matrix(obs$detail, set_dist, method, eta, rho, kappa)
    variance <- as.numeric(t(w) %*% Sigma %*% w)
  }
  pv <- .normal_pvalue(obs$statistic, obs$expectation, variance)

  result <- data.frame(
    method = method,
    eta = if (method == "sensitivity") eta else NA_real_,
    rho = if (method == "sensitivity") rho else NA_real_,
    kappa = kappa,
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
    kappa = kappa,
    outcome = outcome,
    weight_type = weight_type,
    options = list(
      outcome = outcome, method = method, eta = result$eta, rho = result$rho,
      kappa = kappa, weight_type = weight_type
    ),
    match = match
  )
  class(out) <- "netmatch_test"
  out
}

#' @rdname RI_Naive
#' @export
RI_Naive <- function(match,
                     outcome,
                     kappa = NULL,
                     weight_type = c("ns", "ntc")) {
  netmatch_test(
    match = match,
    outcome = outcome,
    method = "naive",
    kappa = kappa,
    weight_type = weight_type
  )
}

#' @rdname RI_Naive
#' @export
RI_Sensitivity <- function(match,
                     outcome,
                     eta = 0.03,
                     rho = 0.10,
                     kappa = NULL,
                     weight_type = c("ns", "ntc")) {
  netmatch_test(
    match = match,
    outcome = outcome,
    method = "sensitivity",
    eta = eta,
    rho = rho,
    kappa = kappa,
    weight_type = weight_type
  )
}

#' @rdname RI_Naive
#' @export
RI_Design <- function(match,
                      outcome,
                      kappa = NULL,
                      weight_type = c("ns", "ntc")) {
  netmatch_test(
    match = match,
    outcome = outcome,
    method = "design",
    kappa = kappa,
    weight_type = weight_type
  )
}

#' @export
print.netmatch_test <- function(x, ...) {
  print(x$result, row.names = FALSE)
  invisible(x)
}

#' Sensitivity Analysis for Network Dependence
#'
#' Evaluates sensitivity analysis results over a grid of eta and rho values.
#' Eta scales residual covariance between matched sets. Rho controls how that
#' covariance changes with distance: for sets at distance `d`, the multiplier
#' is `eta * rho^(d - 1)`.
#'
#' @param match A `netmatch` object.
#' @param outcome Name of the outcome column.
#' @param eta Finite numeric values from 0 to 1.
#' @param rho Finite numeric values from 0 to 1.
#' @param kappa Analysis cutoff. Defaults to the matching cutoff stored in
#'   `match$kappa`.
#' @param weight_type Matched-set weighting. `"ns"` uses `1 / (n_s + 1)`;
#'   `"ntc"` uses `1 / (n_ts * n_cs)`.
#' @return A `netmatch_sensitivity` object with a p-value `grid`, source
#'   `match`, and the complete analysis options.
#' @details P-values use a two-sided normal approximation. The set statistic
#'   counts a treated outcome only when it is strictly greater than a control
#'   outcome; tied outcomes add zero.
#' @examples
#' \dontrun{
#' sim <- simulate_netmatch_example()
#' m <- netmatch(sim$data, "Z", c("X1", "X2", "X3"), sim$net_dist,
#'               method = "dual", kappa = 2, solver = "auto")
#' sens <- netmatch_sensitivity(
#'   m,
#'   "Y",
#'   eta = seq(0, 0.03, by = 0.01),
#'   rho = seq(0, 0.50, by = 0.25)
#' )
#' sens$grid
#' }
#' @export
netmatch_sensitivity <- function(match,
                                 outcome,
                                 eta = seq(0, 0.10, by = 0.01),
                                 rho = seq(0, 0.50, by = 0.05),
                                 kappa = NULL,
                                 weight_type = c("ns", "ntc")) {
  if (!inherits(match, "netmatch")) {
    stop("`match` must be a netmatch object.", call. = FALSE)
  }
  weight_type <- match.arg(weight_type)
  kappa <- .analysis_kappa(match, kappa)
  .validate_unit_interval(eta, "eta")
  .validate_unit_interval(rho, "rho")
  if (!outcome %in% names(match$data)) stop("`outcome` column not found.", call. = FALSE)

  obs <- .observed_stats(match$data, outcome, match$treat, "subclass", weight_type)
  unit_ids <- as.integer(rownames(match$data))
  set_dist <- .set_distance_matrix(match$network_distance, as.integer(match$data$subclass), unit_ids)
  components <- .variance_components(obs$detail, set_dist, kappa)
  w <- obs$detail$weight

  grid <- expand.grid(eta = eta, rho = rho, KEEP.OUT.ATTRS = FALSE)
  rows <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    variance <- .weighted_variance_from_components(components, w, "sensitivity", grid$eta[i], grid$rho[i])
    pv <- .normal_pvalue(obs$statistic, obs$expectation, variance)
    rows[[i]] <- data.frame(
      method = "sensitivity",
      eta = grid$eta[i],
      rho = grid$rho[i],
      kappa = kappa,
      statistic = obs$statistic,
      expectation = obs$expectation,
      variance = variance,
      z_score = pv$z,
      p_value = pv$p,
      stringsAsFactors = FALSE
    )
  }
  out <- list(
    grid = do.call(rbind, rows), match = match, outcome = outcome,
    kappa = kappa, weight_type = weight_type,
    options = list(
      outcome = outcome, eta = eta, rho = rho, kappa = kappa,
      weight_type = weight_type
    )
  )
  class(out) <- "netmatch_sensitivity"
  out
}

#' @export
print.netmatch_sensitivity <- function(x, ...) {
  print(utils::head(x$grid), row.names = FALSE)
  cat("...\n")
  invisible(x)
}

#' Critical Sensitivity Curve for Dual-Penalty Matching
#'
#' Computes the critical value of eta as a function of rho that solves
#' `p(eta, rho) = alpha` for the two-sided normal-approximation analysis.
#'
#' @param match A `netmatch` object, ideally from `method = "dual"`.
#' @param outcome Name of the outcome column.
#' @param rho Finite numeric values from 0 to 1. For matched sets at distance
#'   `d`, residual covariance is multiplied by `eta * rho^(d - 1)`.
#' @param alpha Test level.
#' @param kappa Analysis cutoff. Defaults to `match$kappa`.
#' @param weight_type Matched-set weighting. `"ns"` uses `1 / (n_s + 1)`;
#'   `"ntc"` uses `1 / (n_ts * n_cs)`.
#' @return A `netmatch_critical_sensitivity` object with the critical `curve`,
#'   test level, cutoff, observed statistic, null expectation, diagonal
#'   variance component, naive test result, interpretation text, source
#'   `match`, and the complete analysis options. The curve keeps every
#'   non-negative finite critical eta and marks whether it lies in `[0, 1]`.
#' @examples
#' \dontrun{
#' sim <- simulate_netmatch_example()
#' m <- netmatch(sim$data, "Z", c("X1", "X2", "X3"), sim$net_dist,
#'               method = "dual", kappa = 2, solver = "auto")
#' crit <- critical_sensitivity(m, "Y", rho = seq(0, 1, by = 0.25))
#' crit$curve
#' crit$interpretation
#' }
#' @export
critical_sensitivity <- function(match,
                                 outcome,
                                 rho = seq(0, 1, by = 0.01),
                                 alpha = 0.05,
                                 kappa = NULL,
                                 weight_type = c("ns", "ntc")) {
  if (!inherits(match, "netmatch")) {
    stop("`match` must be a netmatch object.", call. = FALSE)
  }
  .validate_unit_interval(rho, "rho")
  if (!identical(match$method, "dual")) {
    warning("The critical-curve interpretation is designed for dual-penalty matching.", call. = FALSE)
  }
  weight_type <- match.arg(weight_type)
  kappa <- .analysis_kappa(match, kappa)
  if (!is.numeric(alpha) || length(alpha) != 1 || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be one number between 0 and 1.", call. = FALSE)
  }

  if (!outcome %in% names(match$data)) stop("`outcome` column not found.", call. = FALSE)

  obs <- .observed_stats(match$data, outcome, match$treat, "subclass", weight_type)
  unit_ids <- as.integer(rownames(match$data))
  set_dist <- .set_distance_matrix(match$network_distance, as.integer(match$data$subclass), unit_ids)
  components <- .variance_components(obs$detail, set_dist, kappa)
  w <- obs$detail$weight
  v_diag <- sum((w^2) * obs$detail$var)
  delta <- obs$statistic - obs$expectation
  z_alpha <- stats::qnorm(1 - alpha / 2)
  numerator <- delta^2 / z_alpha^2 - v_diag

  curve <- data.frame(
    rho = rho,
    eta_critical = rep(NA_real_, length(rho)),
    eta_in_range = rep(NA, length(rho))
  )
  S <- nrow(components$bound)
  for (r in seq_along(rho)) {
    denom <- 0
    if (S >= 2) {
      for (i in seq_len(S - 1)) {
        for (j in (i + 1):S) {
          if (components$bound[i, j] == 0) next
          denom <- denom + 2 * w[i] * w[j] * rho[r]^(set_dist[i, j] - 1) * components$bound[i, j]
        }
      }
    }
    eta_star <- if (denom > 0) numerator / denom else NA_real_
    if (is.finite(eta_star) && eta_star >= 0) {
      curve$eta_critical[r] <- eta_star
      curve$eta_in_range[r] <- eta_star <= 1
    }
  }

  naive <- RI_Naive(match, outcome, kappa = kappa, weight_type = weight_type)
  out <- list(
    curve = curve,
    alpha = alpha,
    kappa = kappa,
    statistic = obs$statistic,
    expectation = obs$expectation,
    diagonal_variance = v_diag,
    naive = naive$result,
    interpretation = .critical_interpretation(numerator, alpha),
    match = match,
    outcome = outcome,
    weight_type = weight_type,
    options = list(
      outcome = outcome, rho = rho, alpha = alpha, kappa = kappa,
      weight_type = weight_type
    )
  )
  class(out) <- "netmatch_critical_sensitivity"
  out
}

#' @export
print.netmatch_critical_sensitivity <- function(x, ...) {
  print(utils::head(x$curve), row.names = FALSE)
  if (!is.null(x$interpretation)) cat(x$interpretation, "\n")
  cat("...\n")
  invisible(x)
}

#' Plot Sensitivity Results
#'
#' @param x A `netmatch_sensitivity` or `netmatch_critical_sensitivity` object.
#' @param type Plot type: `"critical"` or `"pvalue"`. By default, sensitivity
#'   grid objects use `"pvalue"` and critical objects use `"critical"`.
#' @param alpha Test level for the reference line.
#' @param naive Optional naive p-value or randomization-based inference result object.
#' @param critical_ylim Optional coordinate limits for `type = "critical"`.
#'   Defaults to `c(0, 1)`. Values outside the visible range remain in the data.
#' @param ... Additional arguments passed to `critical_sensitivity()` when
#'   `x` is a grid object and `type = "critical"`.
#' @return A `ggplot` object.
#' @examples
#' \dontrun{
#' sim <- simulate_netmatch_example()
#' m <- netmatch(sim$data, "Z", c("X1", "X2", "X3"), sim$net_dist,
#'               method = "dual", kappa = 2, solver = "auto")
#' sens <- netmatch_sensitivity(m, "Y", eta = seq(0, 0.03, by = 0.03),
#'                              rho = seq(0, 1, by = 0.5))
#' plot_sensitivity(sens)
#'
#' crit <- critical_sensitivity(m, "Y")
#' plot_sensitivity(crit)
#' }
#' @export
plot_sensitivity <- function(x,
                             type = NULL,
                             alpha = 0.05,
                             naive = NULL,
                             critical_ylim = NULL,
                             ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("`ggplot2` is required for sensitivity plots.", call. = FALSE)
  }
  if (is.null(type)) {
    type <- if (inherits(x, "netmatch_critical_sensitivity")) "critical" else "pvalue"
  }
  type <- match.arg(type, c("pvalue", "critical"))
  if (type == "critical") {
    if (inherits(x, "netmatch_sensitivity")) {
      stored_rho <- if (!is.null(x$options$rho)) x$options$rho else sort(unique(x$grid$rho))
      x <- critical_sensitivity(
        x$match, x$outcome, rho = stored_rho, alpha = alpha, kappa = x$kappa,
        weight_type = x$weight_type, ...
      )
    }
    if (!inherits(x, "netmatch_critical_sensitivity")) {
      stop("`type = \"critical\"` requires a critical sensitivity object or sensitivity grid.", call. = FALSE)
    }
    return(.plot_critical_sensitivity(x, ylim = critical_ylim))
  }
  if (!inherits(x, "netmatch_sensitivity")) {
    stop("`type = \"pvalue\"` requires a netmatch_sensitivity object.", call. = FALSE)
  }
  .plot_pvalue_sensitivity(x, alpha = alpha, naive = naive)
}

#' @export
plot.netmatch_sensitivity <- function(x, ...) {
  plot_sensitivity(x, type = "pvalue", ...)
}

#' @export
plot.netmatch_critical_sensitivity <- function(x, ...) {
  plot_sensitivity(x, type = "critical", ...)
}

.plot_pvalue_sensitivity <- function(x, alpha, naive = NULL) {
  grid <- x$grid
  grid$eta_label <- factor(sprintf("%.2f", grid$eta),
                           levels = sprintf("%.2f", sort(unique(grid$eta))))
  p <- ggplot2::ggplot(grid, ggplot2::aes(x = rho, y = p_value, colour = eta_label, group = eta_label)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_hline(yintercept = alpha, linetype = "dotted", colour = "red", linewidth = 0.5) +
    ggplot2::labs(x = expression(rho), y = "p-value", colour = expression(eta)) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
  naive_p <- .extract_naive_p(naive)
  if (is.null(naive_p) && inherits(x$match, "netmatch")) {
    naive_p <- RI_Naive(
      x$match, x$outcome, kappa = x$kappa, weight_type = x$weight_type
    )$result$p_value
  }
  if (!is.null(naive_p) && is.finite(naive_p)) {
    p <- p + ggplot2::geom_hline(yintercept = naive_p, linetype = "dashed", colour = "black", linewidth = 0.7)
  }
  p
}

.plot_critical_sensitivity <- function(x, ylim = NULL) {
  curve <- x$curve[is.finite(x$curve$eta_critical), , drop = FALSE]
  if (!nrow(curve)) {
    stop("No non-negative critical eta boundary is available to plot. Check `x$interpretation`.", call. = FALSE)
  }
  if (is.null(ylim)) {
    ylim <- c(0, 1)
  }
  ggplot2::ggplot(curve, ggplot2::aes(x = rho, y = eta_critical)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = 0, ymax = eta_critical),
                         fill = "#D9EAF7", alpha = 0.6) +
    ggplot2::geom_line(linewidth = 0.9, colour = "#2166AC") +
    ggplot2::coord_cartesian(ylim = ylim) +
    ggplot2::labs(
      x = expression(rho),
      y = expression(eta^"*"),
      title = NULL,
      subtitle = NULL
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
}

.critical_interpretation <- function(numerator, alpha) {
  if (is.finite(numerator) && numerator < 0) {
    return(paste0(
      "The curve (eta, rho) quantifies the minimum residual network dependence ",
      "needed to render the observed significance a spurious association. ",
      "At alpha = ", alpha, ", the naive bound is already not significant, so ",
      "there is no positive critical eta boundary."
    ))
  }
  "The curve (eta, rho) quantifies the minimum residual network dependence needed to render the observed significance a spurious association."
}

.extract_naive_p <- function(naive) {
  if (is.null(naive)) return(NULL)
  if (inherits(naive, "netmatch_test")) return(naive$result$p_value[1])
  if (is.data.frame(naive) && "p_value" %in% names(naive)) return(naive$p_value[1])
  if (is.numeric(naive) && length(naive) == 1) return(naive)
  NULL
}

.validate_unit_interval <- function(x, name, scalar = FALSE) {
  valid <- is.numeric(x) && length(x) > 0 && all(is.finite(x)) &&
    all(x >= 0 & x <= 1)
  if (scalar) valid <- valid && length(x) == 1
  if (!valid) {
    requirement <- if (scalar) "one finite number" else "finite values"
    stop(sprintf("`%s` must %s in [0, 1].", name,
                 if (scalar) paste("be", requirement) else paste("contain", requirement)),
         call. = FALSE)
  }
  invisible(x)
}

.analysis_kappa <- function(match, kappa = NULL) {
  if (is.null(kappa)) kappa <- match$kappa
  if (!is.numeric(kappa) || length(kappa) != 1 || !is.finite(kappa) || kappa < 0) {
    stop("`kappa` must be one non-negative graph-distance cutoff.", call. = FALSE)
  }
  kappa
}
