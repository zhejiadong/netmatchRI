# Independent rectangle oracle: no quantile intervals or production helpers.
rectangle_covariance <- function(a, b) {
  x <- 0:prod(a); y <- 0:prod(b)
  fx <- c(0, stats::pwilcox(x, a[1], a[2]))
  fy <- c(0, stats::pwilcox(y, b[1], b[2]))
  F <- outer(fx, fy, pmin)
  mass <- F[-1, -1, drop = FALSE] - F[-nrow(F), -1, drop = FALSE] -
    F[-1, -ncol(F), drop = FALSE] + F[-nrow(F), -ncol(F), drop = FALSE]
  expect_gte(min(mass), -1e-12)
  expect_equal(sum(mass), 1, tolerance = 1e-12)
  sum(outer(x - prod(a) / 2, y - prod(b) / 2) * mass)
}

bound_for_counts <- function(a, b) {
  netmatchRI:::.design_cov_bound(a[1], a[2], b[1], b[2])
}

test_that("finite-support bounds agree with independent CDF rectangles", {
  types <- list(c(1, 1), c(1, 2), c(1, 4), c(6, 1), c(1, 24),
                c(49, 1), c(1, 74), c(1, 99), c(2, 2), c(2, 3),
                c(3, 4), c(4, 5))
  for (i in seq_along(types)) for (j in i:length(types)) {
    a <- types[[i]]; b <- types[[j]]
    value <- bound_for_counts(a, b)
    expect_true(is.finite(value) && value >= 0)
    expect_equal(value, rectangle_covariance(a, b), tolerance = 1e-11)
    expect_equal(value, bound_for_counts(b, a), tolerance = 1e-14)
    expect_identical(value, bound_for_counts(rev(a), b))
    expect_identical(value, bound_for_counts(a, rev(b)))
  }
  expect_equal(bound_for_counts(c(1, 74), c(1, 74)),
               468.6666666666667, tolerance = 1e-14)
})

test_that("identical marginals give their variance in either orientation", {
  for (n in c(2:20, 25, 50, 75, 100, 100001)) {
    expect_equal(bound_for_counts(c(1, n - 1), c(n - 1, 1)),
                 (n^2 - 1) / 12, tolerance = 1e-14)
  }
  for (a in list(c(2, 2), c(2, 3), c(3, 4), c(4, 5))) {
    # Enumerate rank allocations, independently of Wilcoxon CDF evaluation.
    u <- colSums(combn(seq_len(sum(a)), a[1])) - a[1] * (a[1] + 1) / 2
    expect_equal(bound_for_counts(a, rev(a)), mean((u - mean(u))^2),
                 tolerance = 1e-13)
  }
})

test_that("CDF endpoints and larger supports remain finite and nonnegative", {
  for (pair in list(list(c(1, 999), c(1, 1000)),
                   list(c(20, 20), c(19, 21)),
                   list(c(40, 41), c(39, 42)),
                   list(c(1, 74), c(10, 11)))) {
    a <- pair[[1]]; b <- pair[[2]]
    value <- bound_for_counts(a, b)
    expect_true(is.finite(value) && value >= 0)
    expect_equal(value, rectangle_covariance(a, b), tolerance = 1e-10)
    expect_equal(value, bound_for_counts(b, a), tolerance = 1e-14)
    limit <- sqrt(prod(a) * (sum(a) + 1) / 12 * prod(b) * (sum(b) + 1) / 12)
    expect_lte(value, limit * (1 + 1e-12))
  }
  expect_identical(bound_for_counts(c(0, 3), c(1, 4)), 0)
  for (bad in c(-1, 0.5, Inf, NA_real_)) {
    expect_error(bound_for_counts(c(1, bad), c(1, 2)), "finite non-negative integers")
  }
})

test_that("unadjusted inference bypasses the bound and distance helpers", {
  m <- structure(list(data = data.frame(Z = rep(c(1, 0), 4),
    Y = rep(c(2, 1), 4), subclass = rep(1:4, each = 2)),
    treat = "Z", method = "dual", kappa = 2, network_distance = NULL), class = "netmatch")
  local_mocked_bindings(
    .design_cov_bound = function(...) stop("bound called"),
    .set_distance_matrix = function(...) stop("distance called"),
    .package = "netmatchRI")
  for (wt in c("ns", "ntc")) {
    fit <- RI_unadjusted(m, "Y", weight_type = wt)
    w <- if (wt == "ns") 1 / 3 else 1
    expect_equal(fit$result$statistic, 4 * w)
    expect_equal(fit$result$expectation, 2 * w)
    expect_equal(fit$result$variance, w^2)
    expect_equal(fit$result$p_value, 2 * pnorm(2, lower.tail = FALSE))
    expect_equal(fit$covariance, diag(rep(0.25, 4)), ignore_attr = TRUE)
  }
})

test_that("repeated count pairs keep the existing covariance cache", {
  calls <- 0L
  original_bound <- netmatchRI:::.design_cov_bound
  local_mocked_bindings(.design_cov_bound = function(...) {
    calls <<- calls + 1L
    original_bound(...)
  }, .package = "netmatchRI")
  d <- data.frame(subclass = 1:4, nt = 1, nc = 74,
                  var = rep(468.6666666666667, 4))
  distances <- matrix(1, 4, 4); diag(distances) <- 0
  got <- netmatchRI:::.variance_components(d, distances, 2)
  expect_identical(calls, 1L)
  expect_equal(got$bound, matrix(468.6666666666667, 4, 4), tolerance = 1e-14)
})
