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
                  c("RI_Naive", "RI_Sensitivity", "RI_Design"))
  expect_false(any(c("RI_naive", "RI_decay", "RI_design", "netmatch_test") %in% exports))
})

test_that("canonical RI functions label sensitivity results consistently", {
  m <- ri_api_fixture()
  expect_equal(RI_Naive(m, "Y")$result$method, "naive")
  expect_equal(RI_Sensitivity(m, "Y", eta = 0.2, rho = 0.4)$result$method, "sensitivity")
  expect_equal(RI_Design(m, "Y")$result$method, "design")
})

test_that("eta and rho scalar inputs have direct unit-interval errors", {
  m <- ri_api_fixture()
  expect_error(RI_Sensitivity(m, "Y", eta = c(0.1, 0.2)), "`eta` must be one finite number in [0, 1]", fixed = TRUE)
  expect_error(RI_Sensitivity(m, "Y", eta = Inf), "`eta` must be one finite number in [0, 1]", fixed = TRUE)
  expect_error(RI_Sensitivity(m, "Y", rho = -0.1), "`rho` must be one finite number in [0, 1]", fixed = TRUE)
})

test_that("eta and rho grid inputs have direct unit-interval errors", {
  m <- ri_api_fixture()
  expect_error(netmatch_sensitivity(m, "Y", eta = c(0, NA_real_)), "`eta` must contain finite values in [0, 1]", fixed = TRUE)
  expect_error(netmatch_sensitivity(m, "Y", rho = c(0, 1.1)), "`rho` must contain finite values in [0, 1]", fixed = TRUE)
  expect_error(critical_sensitivity(m, "Y", rho = numeric()), "`rho` must contain finite values in [0, 1]", fixed = TRUE)
})

test_that("sensitivity objects retain analysis-defining options", {
  m <- ri_api_fixture()
  sens <- netmatch_sensitivity(m, "Y", eta = c(0, 0.2), rho = c(0.1, 0.4), kappa = 1, weight_type = "ntc")
  expect_equal(sens$weight_type, "ntc")
  expect_equal(sens$options, list(outcome = "Y", eta = c(0, 0.2), rho = c(0.1, 0.4), kappa = 1, weight_type = "ntc"))

  crit <- suppressWarnings(critical_sensitivity(m, "Y", rho = c(0.1, 0.4), alpha = 0.1, kappa = 1, weight_type = "ntc"))
  expect_equal(crit$weight_type, "ntc")
  expect_equal(crit$options, list(outcome = "Y", rho = c(0.1, 0.4), alpha = 0.1, kappa = 1, weight_type = "ntc"))
})

test_that("sensitivity plot defaults agree and reuse stored weighting", {
  skip_if_not_installed("ggplot2")
  m <- ri_api_fixture()
  sens <- netmatch_sensitivity(m, "Y", eta = c(0, 0.2), rho = c(0.1, 0.4), weight_type = "ntc")
  p1 <- plot(sens)
  p2 <- plot_sensitivity(sens)
  expect_false("..." %in% names(formals(plot_sensitivity)))
  expect_equal(p1$data$p_value, p2$data$p_value)
  expect_true("p_value" %in% names(p2$data))

  critical_plot <- suppressWarnings(plot_sensitivity(sens, type = "critical"))
  expect_equal(sort(unique(critical_plot$data$rho)), c(0.1, 0.4))

  built <- ggplot2::ggplot_build(p2)
  naive_ntc <- RI_Naive(m, "Y", weight_type = "ntc")$result$p_value
  expect_equal(unique(built$data[[4]]$yintercept), naive_ntc)
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
})
