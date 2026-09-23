test_that("adjacency matrices are converted with igraph distances", {
  A <- matrix(0, 4, 4)
  A[cbind(1:3, 2:4)] <- 1
  A[cbind(2:4, 1:3)] <- 1
  D <- netmatchRI:::.as_network_distance(A)
  expect_equal(D[1, 4], 3)
  expect_equal(D[1, 3], 2)
  expect_equal(diag(D), rep(0, 4))
})

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
  fit <- RI_unadjusted(m, "Y")
  expect_s3_class(fit, "netmatch_test")
  expect_true(all(c("statistic", "expectation", "variance", "p_value") %in% names(fit$result)))
  expect_true(is.finite(fit$result$variance))
})

test_that("RI_unadjusted does not require set-distance construction", {
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
  m$network_distance <- NULL
  fit <- RI_unadjusted(m, "Y")
  expect_s3_class(fit, "netmatch_test")
  expect_true(is.finite(fit$result$variance))
})

test_that("sensitivity with eta zero equals unadjusted variance", {
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
  unadjusted <- RI_unadjusted(m, "Y")
  sensitivity <- RI_adjusted(m, "Y", eta = 0, rho = 0)
  expect_equal(sensitivity$result$variance, unadjusted$result$variance)
})

test_that("design equals sensitivity with eta and rho equal to one", {
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
  sensitivity <- RI_adjusted(m, "Y", eta = 1, rho = 1)
  expect_equal(sensitivity$result$variance, design$result$variance)
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

test_that("sensitivity grid matches repeated RI_adjusted calls", {
  D <- matrix(4, 8, 8)
  diag(D) <- 0
  dat <- data.frame(
    Z = c(1, 1, 1, 1, 0, 0, 0, 0),
    X1 = c(0, 1, 2, 3, 0.1, 1.1, 2.1, 3.1),
    X2 = c(1, 1, 2, 2, 1.2, 1.1, 2.2, 2.1),
    Y = c(8, 7, 6, 5, 1, 2, 3, 4)
  )
  m <- netmatch(dat, "Z", c("X1", "X2"), D, method = "covariate", kappa = 2)
  sens <- sensitivity_grid(m, "Y", eta = seq(0, 0.03, by = 0.03), rho = seq(0.1, 0.5, by = 0.4))
  expected <- do.call(rbind, lapply(seq_len(nrow(sens$grid)), function(i) {
    RI_adjusted(m, "Y", eta = sens$grid$eta[i], rho = sens$grid$rho[i])$result
  }))
  rownames(expected) <- NULL
  rownames(sens$grid) <- NULL
  expect_equal(sens$grid, expected)
})

test_that("critical sensitivity returns one curve and reaches alpha", {
  D <- matrix(1, 8, 8); diag(D) <- 0
  dat <- data.frame(Z = rep(c(1, 0), 4), Y = rep(c(2, 1), 4),
                    subclass = rep(1:4, each = 2))
  m <- structure(list(data = dat, treat = "Z", method = "dual", kappa = 2,
                      network_distance = D), class = "netmatch")
  crit <- critical_sensitivity(m, "Y", rho = c(0.1, 0.5), alpha = 0.05)
  expect_s3_class(crit, "netmatch_critical_sensitivity")
  expect_equal(nrow(crit$curve), 2)
  expect_true(all(crit$curve$eta_in_range))
  expect_true(all(crit$curve$eta_critical > 0))
  for (i in seq_len(nrow(crit$curve))) {
    fit <- RI_adjusted(m, "Y", eta = crit$curve$eta_critical[i], rho = crit$curve$rho[i])
    expect_equal(fit$result$p_value, 0.05, tolerance = 1e-12)
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
  crit <- suppressWarnings(critical_sensitivity(m, "Y", rho = seq(0.1, 0.5, by = 0.4)))
  expect_true(all(is.na(crit$curve$eta_critical) | crit$curve$eta_critical >= 0))
  expect_equal(names(crit$curve), c("rho", "eta_critical", "eta_in_range"))
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
  sens <- sensitivity_grid(m, "Y", eta = seq(0, 0.1, by = 0.1), rho = seq(0.1, 0.5, by = 0.4))
  crit <- suppressWarnings(critical_sensitivity(m, "Y", rho = seq(0.1, 0.5, by = 0.4)))
  p_pvalue <- plot_sensitivity(sens, type = "pvalue")
  p_critical <- plot_sensitivity(crit, type = "critical")
  p_critical_full <- plot_sensitivity(crit, type = "critical", critical_ylim = c(0, 1))
  expect_s3_class(p_pvalue, "ggplot")
  expect_s3_class(p_critical, "ggplot")
  expect_s3_class(p_critical_full, "ggplot")
  expect_equal(levels(p_pvalue$data$eta_label), c("0.00", "0.10"))
  expect_equal(p_pvalue$labels$colour, expression(eta))
  expect_null(p_critical$labels$title)
  expect_null(p_critical$labels$subtitle)
  expect_error(plot_sensitivity(sens, alpha = 0), "between 0 and 1")
  expect_error(plot_sensitivity(sens, alpha = 1), "between 0 and 1")
})
