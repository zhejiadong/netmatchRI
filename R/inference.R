#' Randomization-Based Inference Given a Matched Design
#'
#' `RI_unadjusted()`, `RI_adjusted()`, and `RI_design()` test the sharp null of
#' no causal effect using a weighted sum of within-matched-set rank statistics and a
#' normal approximation. `RI_unadjusted()` uses only within-matched-set variances,
#' `RI_adjusted()` allows residual covariance between matched sets, and
#' `RI_design()` uses the design-based covariance bound within the network distance threshold.
#'
#' @param match A \code{netmatch} object.
#' @param outcome Name of the numeric outcome column.
#' @param eta Relative magnitude of the design-based covariance bound in \code{[0, 1]}. Default: \code{0.03}.
#' @param rho Decay parameter in \code{[0, 1]}. Default: \code{0.10}.
#' @param kappa Network distance threshold. \code{NULL} uses \code{match$kappa}.
#' @param weight_type Matched-set weighting: \code{"ns"} or \code{"ntc"}. Default: \code{"ns"}.
#' @return A \code{netmatch_test} object containing the test result,
#'   matched-set details, covariance and distance matrices, settings,
#'   and original matched design.
#' @details For two matched sets at distance `d`, sensitivity analysis multiplies
#'   the design-based covariance bound by `eta * rho^(d - 1)`. The set statistic counts a
#'   treated outcome only when it is strictly greater than a control outcome;
#'   the null moments assume no within-matched-set ties. The design-based covariance bound is evaluated by an exact
#'   finite-support sum, subject to floating-point arithmetic and numerical
#'   Wilcoxon probabilities, with uniform marginals for full matching.
#' @examples
#' \dontrun{
#' sim <- simulate_netmatch_example()
#' m <- netmatch(sim$data, "Z", c("X1", "X2", "X3"), sim$net_dist,
#'               method = "dual", kappa = 2, solver = "auto")
#' RI_unadjusted(m, "Y")
#' RI_adjusted(m, "Y", eta = 0.03, rho = 0.10)
#' RI_design(m, "Y")
#' }
#' @name RI_unadjusted
NULL

netmatch_test <- function(match,
                          outcome,
                          method = c("unadjusted", "adjusted", "design"),
                          eta = 0.03,
                          rho = 0.10,
                          kappa = NULL,
                          weight_type = c("ns", "ntc")) {
  if (!inherits(match, "netmatch")) {
    stop("`match` must be a netmatch object.", call. = FALSE)
  }
  method <- match.arg(method)
  weight_type <- match.arg(weight_type)
  if (method == "adjusted") {
    .validate_unit_interval(eta, "eta", scalar = TRUE)
    .validate_unit_interval(rho, "rho", scalar = TRUE)
  }
  if (!outcome %in% names(match$data)) stop("`outcome` column not found.", call. = FALSE)
  kappa <- .analysis_kappa(match, kappa)

  obs <- .observed_stats(match$data, outcome, match$treat, "subclass", weight_type)
  w <- obs$detail$weight
  if (method == "unadjusted") {
    set_dist <- .empty_set_distance(obs$detail$subclass)
    Sigma <- .unadjusted_covariance_matrix(obs$detail)
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
    eta = if (method == "adjusted") eta else NA_real_,
    rho = if (method == "adjusted") rho else NA_real_,
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

#' @rdname RI_unadjusted
#' @export
RI_unadjusted <- function(match,
                     outcome,
                     kappa = NULL,
                     weight_type = c("ns", "ntc")) {
  netmatch_test(
    match = match,
    outcome = outcome,
    method = "unadjusted",
    kappa = kappa,
    weight_type = weight_type
  )
}

#' @rdname RI_unadjusted
#' @export
RI_adjusted <- function(match,
                     outcome,
                     eta = 0.03,
                     rho = 0.10,
                     kappa = NULL,
                     weight_type = c("ns", "ntc")) {
  netmatch_test(
    match = match,
    outcome = outcome,
    method = "adjusted",
    eta = eta,
    rho = rho,
    kappa = kappa,
    weight_type = weight_type
  )
}

#' @rdname RI_unadjusted
#' @export
RI_design <- function(match,
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
#' Evaluates how the p-value changes over a grid of eta
#' and rho values. Eta scales the design-based covariance bound between matched sets. Rho
#' controls how that bound decays with network distance: for sets at distance
#' `d`, the multiplier is `eta * rho^(d - 1)`.
#'
#' @param match A \code{netmatch} object.
#' @param outcome Name of the numeric outcome column.
#' @param eta Sensitivity-scale values in \code{[0, 1]}. Default: \code{seq(0, 0.10, by = 0.01)}.
#' @param rho Sensitivity-decay values in \code{[0, 1]}. Default: \code{seq(0, 0.50, by = 0.05)}.
#' @param kappa Network distance threshold. \code{NULL} uses \code{match$kappa}.
#' @param weight_type Matched-set weighting: \code{"ns"} or \code{"ntc"}. Default: \code{"ns"}.
#' @return A \code{netmatch_sensitivity} object containing the p-value grid,
#'   original matched design, and analysis settings.
#' @details The set statistic
#'   counts a treated outcome only when it is strictly greater than a control
#'   outcome. The null moments assume no within-matched-set ties.
#' @examples
#' \dontrun{
#' sim <- simulate_netmatch_example()
#' m <- netmatch(sim$data, "Z", c("X1", "X2", "X3"), sim$net_dist,
#'               method = "dual", kappa = 2, solver = "auto")
#' sens <- sensitivity_grid(
#'   m,
#'   "Y",
#'   eta = seq(0, 0.5, by = 0.1),
#'   rho = seq(0, 0.3, by = 0.05)
#' )
#' sens$grid
#' }
#' @export
sensitivity_grid <- function(match,
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
    variance <- .weighted_variance_from_components(components, w, "adjusted", grid$eta[i], grid$rho[i])
    pv <- .normal_pvalue(obs$statistic, obs$expectation, variance)
    rows[[i]] <- data.frame(
      method = "adjusted",
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
#' `p(eta, rho) = alpha`.
#'
#' @param match A \code{netmatch} object. Critical-curve interpretation is intended for dual-penalty matching.
#' @param outcome Name of the numeric outcome column.
#' @param rho Sensitivity-decay values in \code{[0, 1]}. Default: \code{seq(0, 1, by = 0.01)}.
#' @param alpha Test level. Default: \code{0.05}.
#' @param kappa Network distance threshold. \code{NULL} uses \code{match$kappa}.
#' @param weight_type Matched-set weighting: \code{"ns"} or \code{"ntc"}. Default: \code{"ns"}.
#' @return A \code{netmatch_critical_sensitivity} object. \code{summary}
#'   reports the critical eta at \code{rho = 1} and the critical ratio;
#'   \code{curve} contains the critical eta over the requested rho values;
#'   \code{detail} contains the supporting counts and averages. The object also
#'   retains the analysis settings and matched design.
#' @details For matched sets at distance \eqn{d_{s,l}}, the sensitivity bound is
#'   \eqn{M^*_{sl}(\eta,\rho)=\eta\rho^{d_{s,l}-1}M_{sl}}. At a critical
#'   contour point, eta solves \code{p(eta, rho) = alpha}. Values below the curve
#'   retain rejection and values above it do not.
#'
#'   The critical ratio compares the average weighted critical sensitivity
#'   bound across relevant unordered matched-set pairs with the average weighted
#'   diagonal variance. Its numerator averages one unordered-pair contribution
#'   and therefore excludes the factor of two in the full variance expansion.
#'   The ratio is constant along the critical contour because the p-value
#'   threshold fixes the total weighted contribution.
#' @examples
#' \dontrun{
#' sim <- simulate_netmatch_example()
#' m <- netmatch(sim$data, "Z", c("X1", "X2", "X3"), sim$net_dist,
#'               method = "dual", kappa = 2, solver = "auto")
#' crit <- critical_sensitivity(m, "Y", rho = seq(0, 1, by = 0.25))
#' crit$curve
#' crit$summary
#' crit$detail
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
  z_alpha <- stats::qnorm(alpha / 2, lower.tail = FALSE)
  # log.p also handles the smallest positive alpha if alpha / 2 underflows.
  if (!is.finite(z_alpha)) {
    z_alpha <- stats::qnorm(log(alpha) - log(2), lower.tail = FALSE, log.p = TRUE)
  }
  target_variance <- delta^2 / z_alpha^2
  numerator <- target_variance - v_diag
  # Preserve exact baseline p = alpha after a pnorm/qnorm round trip.
  baseline_p <- .normal_pvalue(obs$statistic, obs$expectation, v_diag)$p
  if (is.finite(baseline_p) && baseline_p == alpha) numerator <- 0

  curve <- data.frame(
    rho = rho,
    eta_critical = rep(NA_real_, length(rho)),
    eta_in_range = rep(NA, length(rho))
  )
  S <- nrow(components$bound)
  denominators <- numeric(length(rho))
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
    denominators[r] <- denom
    eta_star <- if (is.finite(v_diag) && v_diag > 0 && is.finite(denom) && denom > 0) {
      numerator / denom
    } else NA_real_
    if (is.finite(eta_star) && eta_star >= 0) {
      curve$eta_critical[r] <- eta_star
      curve$eta_in_range[r] <- eta_star <= 1
    }
  }

  critical_summary <- .critical_summary(
    components = components,
    weights = w,
    diagonal_variance = v_diag,
    numerator = numerator
  )
  unadjusted <- RI_unadjusted(match, outcome, kappa = kappa, weight_type = weight_type)
  out <- list(
    curve = curve,
    summary = critical_summary$summary,
    detail = critical_summary$detail,
    alpha = alpha,
    kappa = kappa,
    statistic = obs$statistic,
    expectation = obs$expectation,
    diagonal_variance = v_diag,
    unadjusted = unadjusted$result,
    interpretation = .critical_interpretation(
      numerator, v_diag, critical_summary$denominator,
      critical_summary$detail$n_relevant_pairs, critical_summary$summary$eta_critical
    ),
    curve_message = .critical_interpretation(
      numerator, v_diag, denominators, critical_summary$detail$n_relevant_pairs,
      curve$eta_critical, requested = TRUE
    ),
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

.critical_summary <- function(components,
                              weights,
                              diagonal_variance,
                              numerator) {
  bound <- components$bound
  set_dist <- components$set_dist
  n_sets <- nrow(bound)
  relevant <- upper.tri(bound) & is.finite(set_dist) &
    set_dist <= components$kappa
  pair_index <- which(relevant, arr.ind = TRUE)
  n_pairs <- nrow(pair_index)

  denominator <- 0
  if (n_pairs > 0) {
    denominator <- sum(
      2 * weights[pair_index[, 1]] * weights[pair_index[, 2]] *
        bound[pair_index]
    )
  }
  eta_critical <- if (is.finite(diagonal_variance) && diagonal_variance > 0 &&
                      is.finite(denominator) && denominator > 0) {
    numerator / denominator
  } else {
    NA_real_
  }
  if (!is.finite(eta_critical) || eta_critical < 0) {
    eta_critical <- NA_real_
  }
  eta_in_range <- if (is.na(eta_critical)) NA else eta_critical <= 1

  average_diagonal <- if (n_sets > 0 && is.finite(diagonal_variance)) diagonal_variance / n_sets else NA_real_
  across_critical <- if (!isTRUE(eta_in_range)) {
    NA_real_
  } else {
    eta_critical * denominator
  }
  average_covariance_bound <- if (n_pairs > 0 && is.finite(across_critical)) {
    across_critical / (2 * n_pairs)
  } else {
    NA_real_
  }
  relative <- if (is.finite(average_covariance_bound) &&
                  is.finite(average_diagonal) && average_diagonal > 0) {
    average_covariance_bound / average_diagonal
  } else {
    NA_real_
  }

  list(
    summary = data.frame(
      rho = 1,
      eta_critical = eta_critical,
      critical_ratio = relative
    ),
    detail = data.frame(
      n_sets = n_sets,
      n_relevant_pairs = n_pairs,
      average_weighted_diagonal_variance = average_diagonal,
      average_weighted_critical_sensitivity_bound = average_covariance_bound
    ),
    denominator = denominator
  )
}

#' @export
print.netmatch_critical_sensitivity <- function(x, ...) {
  cat("Critical eta (rho = 1): ", format(x$summary$eta_critical, digits = 4), "\n", sep = "")
  ratio <- x$summary$critical_ratio
  cat("Critical ratio: ", if (is.finite(ratio)) sprintf("%.2f%%", 100 * ratio) else "NA", "\n", sep = "")
  cat("Settings: alpha = ", x$alpha, "; kappa = ", x$kappa, "\n", sep = "")
  invisible(x)
}

#' Plot Sensitivity Results
#'
#' @param x A sensitivity-grid or critical-sensitivity object.
#' @param type \code{"pvalue"} or \code{"critical"}. The default follows the class of \code{x}.
#' @param alpha Reference or test level. Default: \code{0.05}.
#' @param unadjusted Optional unadjusted p-value or test result for the p-value plot.
#' @param critical_ylim Optional y-axis limits for the critical plot. Default: \code{c(0, 1)}.
#' @return A `ggplot` object.
#' @details A p-value plot shows how the p-value changes with
#'   \code{rho} for each selected value of \code{eta}. A critical plot shows
#'   the boundary value of \code{eta}; values below the curve retain rejection
#'   at the selected \code{alpha}, and values above it do not. If no finite
#'   boundary exists on the requested grid, the function returns an error with
#'   the reason rather than drawing a curve. Algebraic crossings outside the
#'   displayed range remain in the plot data but are not critical contour
#'   points within the parameter range.
#' @examples
#' \dontrun{
#' sim <- simulate_netmatch_example()
#' m <- netmatch(sim$data, "Z", c("X1", "X2", "X3"), sim$net_dist,
#'               method = "dual", kappa = 2, solver = "auto")
#' sens <- sensitivity_grid(m, "Y", eta = seq(0, 0.5, by = 0.1),
#'                          rho = seq(0, 0.3, by = 0.05))
#' plot_sensitivity(sens)
#'
#' crit <- critical_sensitivity(m, "Y")
#' plot_sensitivity(crit)
#' }
#' @export
plot_sensitivity <- function(x,
                             type = NULL,
                             alpha = 0.05,
                             unadjusted = NULL,
                             critical_ylim = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("`ggplot2` is required for sensitivity plots.", call. = FALSE)
  }
  if (!is.numeric(alpha) || length(alpha) != 1 || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be one number between 0 and 1.", call. = FALSE)
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
        weight_type = x$weight_type
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
  .plot_pvalue_sensitivity(x, alpha = alpha, unadjusted = unadjusted)
}

#' @export
plot.netmatch_sensitivity <- function(x, ...) {
  plot_sensitivity(x, type = "pvalue", ...)
}

#' @export
plot.netmatch_critical_sensitivity <- function(x, ...) {
  plot_sensitivity(x, type = "critical", ...)
}

.plot_pvalue_sensitivity <- function(x, alpha, unadjusted = NULL) {
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
  unadjusted_p <- .extract_unadjusted_p(unadjusted)
  if (is.null(unadjusted_p) && inherits(x$match, "netmatch")) {
    unadjusted_p <- RI_unadjusted(
      x$match, x$outcome, kappa = x$kappa, weight_type = x$weight_type
    )$result$p_value
  }
  if (!is.null(unadjusted_p) && is.finite(unadjusted_p)) {
    p <- p + ggplot2::geom_hline(yintercept = unadjusted_p, linetype = "dashed", colour = "black", linewidth = 0.7)
  }
  p
}

.plot_critical_sensitivity <- function(x, ylim = NULL) {
  curve <- x$curve[is.finite(x$curve$eta_critical), , drop = FALSE]
  if (!nrow(curve)) {
    reason <- x$curve_message
    if (is.null(reason)) reason <- "No finite critical eta values are available on the requested rho grid."
    stop(paste("No critical curve is available.", reason), call. = FALSE)
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

.critical_interpretation <- function(numerator, diagonal_variance, denominator,
                                     n_pairs, eta, requested = FALSE) {
  if (!is.finite(diagonal_variance) || diagonal_variance <= 0) {
    return("The diagonal variance is zero or non-finite; the normal-approximation boundary and critical ratio are undefined.")
  }
  if (!is.finite(numerator) || any(!is.finite(denominator))) {
    return("Non-finite boundary components prevent a critical eta calculation.")
  }
  if (numerator < 0) {
    return("The unadjusted p-value is greater than alpha; there is no critical contour in [0, 1].")
  }
  if (all(denominator == 0)) {
    reason <- if (n_pairs == 0) "There are no relevant matched-set pairs." else if (requested) {
      "All denominators on the requested rho grid are zero."
    } else "The denominator at rho = 1 is zero."
    conclusion <- if (numerator == 0) {
      "Every eta yields p = alpha; there is no unique critical eta boundary."
    } else {
      "The p-value remains below alpha for every eta; there is no unique critical eta boundary."
    }
    return(paste(reason, conclusion))
  }
  if (numerator == 0) {
    return(paste("The unadjusted p-value equals alpha. Critical eta is 0 where the denominator is positive;",
                 "where it is zero, every eta yields p = alpha. The critical ratio is 0."))
  }
  if (!any(is.finite(eta) & eta <= 1)) {
    return("Critical eta exceeds 1: there is no crossing in [0, 1], and the p-value remains below alpha. The critical ratio is undefined in this domain.")
  }
  paste("On the critical contour, p = alpha; below it rejection is retained and above it rejection is not retained.",
        "The critical ratio compares the average weighted critical sensitivity bound with the average weighted diagonal variance;",
        "it is constant along the contour.")
}

.extract_unadjusted_p <- function(unadjusted) {
  if (is.null(unadjusted)) return(NULL)
  if (inherits(unadjusted, "netmatch_test")) return(unadjusted$result$p_value[1])
  if (is.data.frame(unadjusted) && "p_value" %in% names(unadjusted)) return(unadjusted$p_value[1])
  if (is.numeric(unadjusted) && length(unadjusted) == 1) return(unadjusted)
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
    stop("`kappa` must be one non-negative network distance threshold.", call. = FALSE)
  }
  kappa
}
