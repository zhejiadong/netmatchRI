test_that("RI helpers return tidy inference output", {
  A <- matrix(0, 8, 8)
  A[cbind(1:7, 2:8)] <- 1
  A[cbind(2:8, 1:7)] <- 1
  dat <- data.frame(
    Z = c(1, 1, 1, 1, 0, 0, 0, 0),
    X1 = c(0, 1, 2, 3, 0.1, 1.1, 2.1, 3.1),
    X2 = c(1, 1, 2, 2, 1.2, 1.1, 2.2, 2.1),
    Y = c(3, 4, 5, 6, 2, 3, 4, 5)
  )
  m <- netmatch(dat, "Z", c("X1", "X2"), A, method = "covariate")
  fit <- RI_naive(m, "Y")
  expect_s3_class(fit, "netmatch_test")
  expect_true(all(c("statistic", "expectation", "variance", "p_value") %in% names(fit$result)))
  expect_true(is.finite(fit$result$variance))
})

test_that("decay with eta zero equals naive variance", {
  A <- matrix(0, 8, 8)
  A[cbind(1:7, 2:8)] <- 1
  A[cbind(2:8, 1:7)] <- 1
  dat <- data.frame(
    Z = c(1, 1, 1, 1, 0, 0, 0, 0),
    X1 = c(0, 1, 2, 3, 0.1, 1.1, 2.1, 3.1),
    X2 = c(1, 1, 2, 2, 1.2, 1.1, 2.2, 2.1),
    Y = c(3, 4, 5, 6, 2, 3, 4, 5)
  )
  m <- netmatch(dat, "Z", c("X1", "X2"), A, method = "covariate")
  naive <- RI_naive(m, "Y")
  decay <- RI_decay(m, "Y", eta = 0, rho = 0)
  expect_equal(decay$result$variance, naive$result$variance)
})

test_that("design equals decay with eta and rho equal to one", {
  A <- matrix(0, 8, 8)
  A[cbind(1:7, 2:8)] <- 1
  A[cbind(2:8, 1:7)] <- 1
  dat <- data.frame(
    Z = c(1, 1, 1, 1, 0, 0, 0, 0),
    X1 = c(0, 1, 2, 3, 0.1, 1.1, 2.1, 3.1),
    X2 = c(1, 1, 2, 2, 1.2, 1.1, 2.2, 2.1),
    Y = c(3, 4, 5, 6, 2, 3, 4, 5)
  )
  m <- netmatch(dat, "Z", c("X1", "X2"), A, method = "covariate", kappa = 2)
  design <- RI_design(m, "Y")
  decay <- RI_decay(m, "Y", eta = 1, rho = 1)
  expect_equal(decay$result$variance, design$result$variance)
})

test_that("design covariance truncates at kappa", {
  stats_df <- data.frame(
    subclass = 1:3,
    n = c(2, 2, 2),
    nt = c(1, 1, 1),
    nc = c(1, 1, 1),
    U = c(1, 1, 1),
    mu = c(0.5, 0.5, 0.5),
    var = c(0.25, 0.25, 0.25),
    weight = c(1, 1, 1)
  )
  set_dist <- matrix(c(
    0, 1, 3,
    1, 0, 4,
    3, 4, 0
  ), nrow = 3, byrow = TRUE)
  Sigma <- netmatchRI:::.covariance_matrix(stats_df, set_dist, method = "design", kappa = 2)
  expect_gt(Sigma[1, 2], 0)
  expect_equal(Sigma[1, 3], 0)
  expect_equal(Sigma[2, 3], 0)
})

test_that("critical sensitivity returns one curve and reaches alpha", {
  skip_without_gurobi()
  D <- matrix(4, 8, 8)
  diag(D) <- 0
  dat <- data.frame(
    Z = c(1, 1, 1, 1, 0, 0, 0, 0),
    X1 = c(0, 1, 2, 3, 0.1, 1.1, 2.1, 3.1),
    X2 = c(1, 1, 2, 2, 1.2, 1.1, 2.2, 2.1),
    Y = c(8, 7, 6, 5, 1, 2, 3, 4)
  )
  m <- netmatch(dat, "Z", c("X1", "X2"), D, method = "dual", kappa = 2)
  crit <- critical_sensitivity(m, "Y", rho = c(0.1, 0.5), alpha = 0.05)
  expect_s3_class(crit, "netmatch_critical_sensitivity")
  expect_equal(nrow(crit$curve), 2)
  eta_star <- crit$curve$eta_critical[is.finite(crit$curve$eta_critical) & crit$curve$eta_critical > 0][1]
  rho_star <- crit$curve$rho[is.finite(crit$curve$eta_critical) & crit$curve$eta_critical > 0][1]
  if (is.finite(eta_star) && eta_star <= 1) {
    fit <- RI_decay(m, "Y", eta = eta_star, rho = rho_star)
    expect_equal(fit$result$p_value, 0.05, tolerance = 1e-6)
  }
})

test_that("critical sensitivity does not report negative eta boundaries", {
  D <- matrix(4, 8, 8)
  diag(D) <- 0
  dat <- data.frame(
    Z = c(1, 1, 1, 1, 0, 0, 0, 0),
    X1 = c(0, 1, 2, 3, 0.1, 1.1, 2.1, 3.1),
    X2 = c(1, 1, 2, 2, 1.2, 1.1, 2.2, 2.1),
    Y = c(3, 4, 5, 6, 2, 3, 4, 5)
  )
  m <- netmatch(dat, "Z", c("X1", "X2"), D, method = "covariate", kappa = 2)
  crit <- suppressWarnings(critical_sensitivity(m, "Y", rho = c(0.1, 0.5)))
  expect_true(all(is.na(crit$curve$eta_critical) | crit$curve$eta_critical >= 0))
  expect_equal(names(crit$curve), c("rho", "eta_critical"))
})

test_that("sensitivity plots return ggplot objects", {
  testthat::skip_if_not_installed("ggplot2")
  A <- matrix(0, 8, 8)
  A[cbind(1:7, 2:8)] <- 1
  A[cbind(2:8, 1:7)] <- 1
  dat <- data.frame(
    Z = c(1, 1, 1, 1, 0, 0, 0, 0),
    X1 = c(0, 1, 2, 3, 0.1, 1.1, 2.1, 3.1),
    X2 = c(1, 1, 2, 2, 1.2, 1.1, 2.2, 2.1),
    Y = c(8, 7, 6, 5, 1, 2, 3, 4)
  )
  m <- netmatch(dat, "Z", c("X1", "X2"), A, method = "covariate", kappa = 2)
  sens <- netmatch_sensitivity(m, "Y", eta = c(0, 0.1), rho = c(0.1, 0.5))
  crit <- suppressWarnings(critical_sensitivity(m, "Y", rho = c(0.1, 0.5)))
  expect_s3_class(plot_sensitivity(sens, type = "pvalue"), "ggplot")
  expect_s3_class(plot_sensitivity(crit, type = "critical"), "ggplot")
  expect_s3_class(plot_sensitivity(crit, type = "critical", critical_ylim = c(0, 1)), "ggplot")
})
