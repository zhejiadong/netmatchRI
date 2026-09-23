ri_api_fixture <- function(method = "covariate") {
  D <- matrix(1, 8, 8)
  diag(D) <- 0
  dat <- data.frame(
    Z = c(1, 1, 1, 1, 0, 0, 0, 0),
    X1 = c(0, 1, 2, 3, 0.1, 1.1, 2.1, 3.1),
    X2 = c(1, 1, 2, 2, 1.2, 1.1, 2.2, 2.1),
    Y = c(8, 7, 6, 5, 1, 2, 3, 4)
  )
  netmatch(dat, "Z", c("X1", "X2"), D, method = method, kappa = 2)
}

test_that("the public RI API uses only the canonical names", {
  exports <- getNamespaceExports("netmatchRI")
  expect_setequal(grep("^RI_", exports, value = TRUE),
                  c("RI_unadjusted", "RI_adjusted", "RI_design"))
  expect_false(any(c("RI_decay", "netmatch_test") %in% exports))
  expect_true("sensitivity_grid" %in% exports)
  expect_false("netmatch_sensitivity" %in% exports)
  expect_false(exists("netmatch_sensitivity", envir = asNamespace("netmatchRI"),
                      inherits = FALSE))
})

test_that("canonical RI functions label sensitivity results consistently", {
  m <- ri_api_fixture()
  expect_equal(RI_unadjusted(m, "Y")$result$method, "unadjusted")
  expect_equal(RI_adjusted(m, "Y", eta = 0.2, rho = 0.4)$result$method, "adjusted")
  expect_equal(RI_design(m, "Y")$result$method, "design")
})

test_that("inference method settings use only the canonical labels", {
  m <- ri_api_fixture()
  fits <- list(
    RI_unadjusted(m, "Y"),
    RI_adjusted(m, "Y", eta = 0.2, rho = 0.4),
    RI_design(m, "Y")
  )
  expect_identical(vapply(fits, function(x) x$options$method, character(1)),
                   c("unadjusted", "adjusted", "design"))
  grid <- sensitivity_grid(m, "Y", eta = 0.2, rho = 0.4)$grid
  expect_identical(unique(grid$method), "adjusted")
})

test_that("eta and rho scalar inputs have direct unit-interval errors", {
  m <- ri_api_fixture()
  expect_error(RI_adjusted(m, "Y", eta = c(0.1, 0.2)), "`eta` must be one finite number in [0, 1]", fixed = TRUE)
  expect_error(RI_adjusted(m, "Y", eta = Inf), "`eta` must be one finite number in [0, 1]", fixed = TRUE)
  expect_error(RI_adjusted(m, "Y", rho = -0.1), "`rho` must be one finite number in [0, 1]", fixed = TRUE)
})

test_that("eta and rho grid inputs have direct unit-interval errors", {
  m <- ri_api_fixture()
  expect_error(sensitivity_grid(m, "Y", eta = c(0, NA_real_)), "`eta` must contain finite values in [0, 1]", fixed = TRUE)
  expect_error(sensitivity_grid(m, "Y", rho = c(0, 1.1)), "`rho` must contain finite values in [0, 1]", fixed = TRUE)
  expect_error(critical_sensitivity(m, "Y", rho = numeric()), "`rho` must contain finite values in [0, 1]", fixed = TRUE)
})

test_that("sensitivity objects retain analysis-defining options", {
  m <- ri_api_fixture()
  sens <- sensitivity_grid(m, "Y", eta = c(0, 0.2), rho = c(0.1, 0.4), kappa = 1, weight_type = "ntc")
  expect_equal(sens$weight_type, "ntc")
  expect_equal(sens$options, list(outcome = "Y", eta = c(0, 0.2), rho = c(0.1, 0.4), kappa = 1, weight_type = "ntc"))

  crit <- suppressWarnings(critical_sensitivity(m, "Y", rho = c(0.1, 0.4), alpha = 0.1, kappa = 1, weight_type = "ntc"))
  expect_equal(crit$weight_type, "ntc")
  expect_equal(crit$options, list(outcome = "Y", rho = c(0.1, 0.4), alpha = 0.1, kappa = 1, weight_type = "ntc"))
})

test_that("sensitivity plot defaults agree and reuse stored weighting", {
  skip_if_not_installed("ggplot2")
  m <- ri_api_fixture()
  sens <- sensitivity_grid(m, "Y", eta = c(0, 0.2), rho = c(0.1, 0.4), weight_type = "ntc")
  p1 <- plot(sens)
  p2 <- plot_sensitivity(sens)
  expect_false("..." %in% names(formals(plot_sensitivity)))
  expect_equal(p1$data$p_value, p2$data$p_value)
  expect_true("p_value" %in% names(p2$data))

  critical_plot <- suppressWarnings(plot_sensitivity(sens, type = "critical"))
  expect_equal(sort(unique(critical_plot$data$rho)), c(0.1, 0.4))

  built <- ggplot2::ggplot_build(p2)
  unadjusted_ntc <- RI_unadjusted(m, "Y", weight_type = "ntc")$result$p_value
  expect_equal(unique(built$data[[4]]$yintercept), unadjusted_ntc)
})

test_that("critical plots retain out-of-range eta values and use coordinate limits", {
  skip_if_not_installed("ggplot2")
  crit <- structure(
    list(curve = data.frame(rho = c(0.2, 0.4), eta_critical = c(0.8, 1.4), eta_in_range = c(TRUE, FALSE))),
    class = "netmatch_critical_sensitivity"
  )
  p1 <- plot(crit)
  p2 <- plot_sensitivity(crit)
  expect_equal(p1$data$eta_critical, c(0.8, 1.4))
  expect_equal(p2$data$eta_critical, c(0.8, 1.4))
  expect_equal(p2$coordinates$limits$y, c(0, 1))

  built <- ggplot2::ggplot_build(p2)
  expect_length(built$data, 2)
})

test_that("critical summaries use the article definitions", {
  bound <- matrix(
    c(2, 0.5, 0.25,
      0.5, 3, 0,
      0.25, 0, 4),
    nrow = 3,
    byrow = TRUE
  )
  set_dist <- matrix(
    c(0, 1, 2,
      1, 0, 3,
      2, 3, 0),
    nrow = 3,
    byrow = TRUE
  )
  components <- list(bound = bound, set_dist = set_dist, kappa = 2)
  weights <- c(0.5, 0.25, 0.4)
  diagonal <- sum(weights^2 * diag(bound))
  numerator <- 0.1
  denominator <- 2 * (weights[1] * weights[2] * bound[1, 2] +
                        weights[1] * weights[3] * bound[1, 3])

  result <- netmatchRI:::.critical_summary(
    components, weights, diagonal, numerator
  )

  expect_equal(result$summary$eta_critical, numerator / denominator)
  expect_equal(result$detail$n_relevant_pairs, 2)
  expect_equal(
    result$detail$average_weighted_critical_sensitivity_bound,
    numerator / (2 * result$detail$n_relevant_pairs)
  )
  expect_equal(
    result$summary$critical_ratio,
    (numerator / (2 * result$detail$n_relevant_pairs)) /
      (diagonal / result$detail$n_sets)
  )
  expect_named(result$summary, c("rho", "eta_critical", "critical_ratio"))
})

test_that("critical summaries are undefined when unadjusted inference is non-significant", {
  components <- list(
    bound = matrix(c(1, 0.5, 0.5, 1), 2, 2),
    set_dist = matrix(c(0, 1, 1, 0), 2, 2),
    kappa = 1
  )
  result <- netmatchRI:::.critical_summary(
    components, weights = c(1, 1), diagonal_variance = 2,
    numerator = -0.1
  )
  expect_true(is.na(result$summary$eta_critical))
  expect_true(is.na(result$detail$average_weighted_critical_sensitivity_bound))
  expect_true(is.na(result$summary$critical_ratio))
})

test_that("rho = 1 summary is computed even when rho one is not requested", {
  m <- ri_api_fixture()
  m$method <- "dual"
  partial <- critical_sensitivity(m, "Y", rho = c(0.2, 0.6))
  at_one <- critical_sensitivity(m, "Y", rho = 1)
  expect_equal(partial$summary, at_one$summary)
  expect_named(
    partial$detail,
    c("n_sets", "n_relevant_pairs",
      "average_weighted_diagonal_variance",
      "average_weighted_critical_sensitivity_bound")
  )
})
