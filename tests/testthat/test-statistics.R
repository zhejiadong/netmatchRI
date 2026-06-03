test_that("netmatch_test returns tidy inference output", {
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
  fit <- netmatch_test(m, "Y", method = "naive")
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
  naive <- netmatch_test(m, "Y", method = "naive")
  decay <- netmatch_test(m, "Y", method = "decay", eta = 0, rho = 0, d0 = 2)
  expect_equal(decay$result$variance, naive$result$variance)
})
