#' Plot Covariate Balance Across Matching Designs
#'
#' Compares absolute standardized mean differences before matching and across
#' selected matching designs.
#'
#' @param designs Named list of comparable \code{netmatch} objects.
#' @param include_before Include the unmatched sample. Default: \code{TRUE}.
#' @return A `ggplot` object.
#' @examples
#' \dontrun{
#' sim <- simulate_netmatch_example(seed = 123, n = 32)
#' m <- netmatch(sim$data, "Z", c("X1", "X2", "X3"), sim$net_dist,
#'               method = "dual", kappa = 2)
#' plot_covariate_balance(list("Dual-penalty" = m))
#' }
#' @export
plot_covariate_balance <- function(designs, include_before = TRUE) {
  designs <- .validate_plot_designs(designs, include_before)
  first_balance <- diagnose_match(designs[[1]])$covariate_balance
  parts <- list()
  if (include_before) {
    parts[["Before"]] <- data.frame(
      design = "Before",
      covariate = first_balance$covariate,
      abs_smd = first_balance$before_abs_smd,
      stringsAsFactors = FALSE
    )
  }
  for (name in names(designs)) {
    balance <- diagnose_match(designs[[name]])$covariate_balance
    parts[[name]] <- data.frame(
      design = name,
      covariate = balance$covariate,
      abs_smd = balance$after_abs_smd,
      stringsAsFactors = FALSE
    )
  }
  data <- do.call(rbind, parts)
  rownames(data) <- NULL
  design_levels <- c(if (include_before) "Before", names(designs))
  data$design <- factor(data$design, levels = design_levels)
  data$covariate <- factor(
    data$covariate,
    levels = rev(first_balance$covariate)
  )

  ggplot2::ggplot(
    data,
    ggplot2::aes(x = abs_smd, y = covariate, colour = design, shape = design)
  ) +
    ggplot2::geom_vline(
      xintercept = 0.10, linetype = "dashed", colour = "grey50",
      linewidth = 0.5
    ) +
    ggplot2::geom_point(size = 2.4) +
    ggplot2::labs(
      x = "Absolute standardized mean difference",
      y = NULL,
      colour = "Design",
      shape = "Design"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
}

#' Plot Covariate Similarity Across Matching Designs
#'
#' Compares average treated-control Mahalanobis distance before matching and
#' within matched sets across selected matching designs. Smaller
#' values indicate greater covariate similarity.
#'
#' @inheritParams plot_covariate_balance
#' @return A `ggplot` object.
#' @examples
#' \dontrun{
#' sim <- simulate_netmatch_example(seed = 123, n = 32)
#' m <- netmatch(sim$data, "Z", c("X1", "X2", "X3"), sim$net_dist,
#'               method = "dual", kappa = 2)
#' plot_covariate_similarity(list("Dual-penalty" = m))
#' }
#' @export
plot_covariate_similarity <- function(designs, include_before = TRUE) {
  designs <- .validate_plot_designs(designs, include_before)
  first <- designs[[1]]
  data <- first$original_data
  z <- data[[first$treat]]
  mahalanobis <- .mahalanobis_matrix(
    data = data,
    treat = z,
    covariates = first$covariates,
    cov_type = "pooled"
  )
  treated_ids <- which(z == 1)
  control_ids <- which(z == 0)
  rows <- list()

  if (include_before) {
    before_values <- mahalanobis[is.finite(mahalanobis)]
    rows[["Before"]] <- data.frame(
      design = "Before",
      average_mahalanobis = mean(before_values),
      n_pairs = length(before_values),
      stringsAsFactors = FALSE
    )
  }
  for (name in names(designs)) {
    design <- designs[[name]]
    subclass <- design$subclass
    matched_treated <- which(z == 1 & !is.na(subclass))
    matched_control <- which(z == 0 & !is.na(subclass))
    distance_block <- mahalanobis[
      match(matched_treated, treated_ids),
      match(matched_control, control_ids),
      drop = FALSE
    ]
    same_set <- outer(
      as.character(subclass[matched_treated]),
      as.character(subclass[matched_control]),
      FUN = "=="
    )
    values <- distance_block[same_set & is.finite(distance_block)]
    rows[[name]] <- data.frame(
      design = name,
      average_mahalanobis = if (length(values)) mean(values) else NA_real_,
      n_pairs = length(values),
      stringsAsFactors = FALSE
    )
  }
  plot_data <- do.call(rbind, rows)
  rownames(plot_data) <- NULL
  design_levels <- c(if (include_before) "Before", names(designs))
  plot_data$design <- factor(plot_data$design, levels = design_levels)
  plot_data$label <- ifelse(
    is.finite(plot_data$average_mahalanobis),
    sprintf("%.2f", plot_data$average_mahalanobis),
    "NA"
  )

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = design, y = average_mahalanobis)
  ) +
    ggplot2::geom_col(width = 0.65, fill = "#2166AC") +
    ggplot2::geom_text(ggplot2::aes(label = label), vjust = -0.35, size = 3.4) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.10))) +
    ggplot2::labs(
      x = NULL,
      y = "Average treated-control Mahalanobis distance"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 20, hjust = 1)
    )
}

.validate_plot_designs <- function(designs, include_before) {
  if (!is.list(designs) || !length(designs) ||
      any(!vapply(designs, inherits, logical(1), what = "netmatch"))) {
    stop("`designs` must be a non-empty named list of `netmatch` objects.",
         call. = FALSE)
  }
  design_names <- names(designs)
  if (is.null(design_names) || any(!nzchar(design_names)) ||
      anyDuplicated(design_names)) {
    stop("All elements of `designs` must have unique non-empty names.",
         call. = FALSE)
  }
  if (!is.logical(include_before) || length(include_before) != 1L ||
      is.na(include_before)) {
    stop("`include_before` must be TRUE or FALSE.", call. = FALSE)
  }
  if (include_before && "Before" %in% design_names) {
    stop("`Before` is reserved when `include_before = TRUE`.", call. = FALSE)
  }

  first <- designs[[1]]
  for (i in seq_along(designs)[-1]) {
    design <- designs[[i]]
    same_inputs <- identical(design$treat, first$treat) &&
      identical(design$covariates, first$covariates) &&
      identical(design$original_data, first$original_data)
    if (!same_inputs) {
      stop(
        "All `designs` must use the same original data, treatment, and covariates.",
        call. = FALSE
      )
    }
  }
  designs
}
