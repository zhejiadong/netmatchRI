# Solver-free matched designs; test_local loads these against the live source namespace.
contract_fixture <- function(full = FALSE, distance = 1, n_sets = 4L) {
  counts <- if (full) list(c(1, 1), c(1, 2), c(3, 1), c(1, 4)) else
    rep(list(c(1, 1)), n_sets)
  dat <- do.call(rbind, lapply(seq_along(counts), function(s) {
    nt <- counts[[s]][1]; nc <- counts[[s]][2]
    data.frame(Z = c(rep(1, nt), rep(0, nc)),
               Y = c(10 + seq_len(nt), seq_len(nc)), subclass = s)
  }))
  rownames(dat) <- seq_len(nrow(dat))
  D <- matrix(distance, nrow(dat), nrow(dat)); diag(D) <- 0
  structure(list(data = dat, treat = "Z", method = "dual", kappa = 2,
                 network_distance = D), class = "netmatch")
}

# Source-path identity is verified by the external development runner, not here:
# R CMD check correctly tests the installed namespace in its isolated library.

test_that("critical contract is concise with no prerelease aliases", {
  m <- contract_fixture()
  x <- critical_sensitivity(m, "Y", rho = c(0.2, 0.6))
  expect_named(x$summary, c("rho", "eta_critical", "critical_ratio"))
  expect_equal(nrow(x$summary), 1)
  expect_identical(x$summary$rho, 1)
  expect_named(x$curve, c("rho", "eta_critical", "eta_in_range"))
  expect_named(x$detail, c("n_sets", "n_relevant_pairs",
                          "average_weighted_diagonal_variance",
                          "average_weighted_critical_sensitivity_bound"))
  expect_false(any(c("no_decay", "pair_level", "pair_level_label") %in% names(x)))
  expect_identical(x$summary, critical_sensitivity(m, "Y", rho = 1)$summary)
  text <- capture.output(print(x))
  expect_match(text[1], "^Critical eta \\(rho = 1\\):")
  expect_match(text[2], "^Critical ratio: .*%$")
  expect_match(text[3], "alpha = 0.05; kappa = 2", fixed = TRUE)
  expect_false(any(grepl("n_relevant_pairs|no.decay|pair.level", text)))
  expect_length(text, 3L)
  expect_identical(text[3], "Settings: alpha = 0.05; kappa = 2")
})

test_that("unequal full sets obey contour equality, invariant ratio and endpoints for both weights", {
  m <- contract_fixture(full = TRUE)
  # Some distances > 1 exercise rho decay, while a distance-one pair stays active at rho=0.
  far <- m$data$subclass >= 3
  m$network_distance[far, !far] <- m$network_distance[!far, far] <- 2
  for (wt in c("ns", "ntc")) {
    unadjusted <- RI_unadjusted(m, "Y", weight_type = wt)
    design <- RI_design(m, "Y", weight_type = wt)
    alpha <- RI_adjusted(m, "Y", eta = 0.1, rho = 1, weight_type = wt)$result$p_value
    x <- critical_sensitivity(m, "Y", rho = c(0, 0.25, 0.5, 1), alpha = alpha, weight_type = wt)
    expect_true(all(x$curve$eta_in_range))
    expect_equal(x$summary$eta_critical, 0.1, tolerance = 1e-12)
    w <- design$detail$weight; V <- unadjusted$result$variance
    relevant <- upper.tri(design$covariance) & is.finite(design$set_distance) &
      design$set_distance <= m$kappa
    ij <- which(relevant, arr.ind = TRUE)
    ratios <- numeric(nrow(x$curve))
    for (i in seq_len(nrow(x$curve))) {
      r <- x$curve$rho[i]; eta <- x$curve$eta_critical[i]
      fit <- RI_adjusted(m, "Y", eta = eta, rho = r, weight_type = wt)
      expect_equal(fit$result$p_value, alpha, tolerance = 1e-12)
      ratios[i] <- mean(w[ij[, 1]] * w[ij[, 2]] * eta *
                          r^(design$set_distance[ij] - 1) * design$covariance[ij]) /
        (V / nrow(design$detail))
    }
    expect_equal(ratios, rep(x$summary$critical_ratio, length(ratios)), tolerance = 1e-12)
    grid <- sensitivity_grid(m, "Y", eta = c(0, 1), rho = c(0, 0.5, 1), weight_type = wt)$grid
    for (i in seq_len(nrow(grid))) {
      fit <- RI_adjusted(m, "Y", eta = grid$eta[i], rho = grid$rho[i], weight_type = wt)
      expect_equal(grid$variance[i], fit$result$variance, tolerance = 1e-12)
      expect_equal(grid$p_value[i], fit$result$p_value, tolerance = 1e-12)
    }
    expect_equal(grid$variance[grid$eta == 0], rep(V, 3), tolerance = 1e-12)
    expect_equal(grid$p_value[grid$eta == 0], rep(unadjusted$result$p_value, 3), tolerance = 1e-12)
    expect_equal(grid$variance[grid$eta == 1 & grid$rho == 1], design$result$variance)
    expect_equal(grid$p_value[grid$eta == 1 & grid$rho == 1], design$result$p_value)
    expect_equal(RI_adjusted(m, "Y", eta = 1, rho = 1, weight_type = wt)$covariance,
                 design$covariance)
  }
})

test_that("negative, zero and out-of-domain boundaries have truthful interpretations", {
  m <- contract_fixture()
  p0 <- RI_unadjusted(m, "Y")$result$p_value
  negative <- critical_sensitivity(m, "Y", alpha = p0 / 2)
  expect_true(all(is.na(negative$curve$eta_critical)))
  expect_true(is.na(negative$summary$critical_ratio))
  expect_match(negative$interpretation, "greater than alpha")
  expect_error(plot_sensitivity(negative), "greater than alpha")
  equal <- critical_sensitivity(m, "Y", alpha = p0, rho = c(0, 0.5, 1))
  expect_equal(equal$curve$eta_critical, rep(0, 3))
  expect_true(all(equal$curve$eta_in_range))
  expect_identical(equal$summary$critical_ratio, 0)
  expect_match(equal$interpretation, "equals alpha")
  expect_equal(RI_adjusted(m, "Y", eta = 0)$result$p_value, equal$alpha)
  p1 <- RI_design(m, "Y")$result$p_value
  outside <- critical_sensitivity(m, "Y", alpha = (1 + p1) / 2, rho = c(0.5, 1))
  expect_true(all(outside$curve$eta_critical > 1))
  expect_false(any(outside$curve$eta_in_range))
  expect_true(is.na(outside$summary$critical_ratio))
  expect_true(is.na(outside$detail$average_weighted_critical_sensitivity_bound))
  expect_match(outside$interpretation, "no crossing in [0, 1]", fixed = TRUE)
  p <- plot_sensitivity(outside)
  expect_equal(p$data$eta_critical, outside$curve$eta_critical)
  expect_equal(p$coordinates$limits$y, c(0, 1))
})

test_that("zero denominators distinguish no relevant pairs and requested rho grid", {
  for (case in c("kappa", "disconnected", "one-set")) {
    m <- contract_fixture(n_sets = if (case == "one-set") 1L else 4L)
    if (case == "kappa") m$kappa <- 0
    if (case == "disconnected") {
      m$network_distance[,] <- Inf; diag(m$network_distance) <- 0
    }
    p0 <- RI_unadjusted(m, "Y")$result$p_value
    x <- critical_sensitivity(m, "Y", alpha = (1 + p0) / 2)
    expect_equal(x$detail$n_relevant_pairs, 0)
    expect_true(is.na(x$summary$eta_critical))
    expect_true(is.na(x$summary$critical_ratio))
    expect_match(x$interpretation, "remains below alpha")
    expect_error(plot_sensitivity(x), "no relevant matched-set pairs")
    equality <- critical_sensitivity(m, "Y", alpha = p0)
    expect_true(is.na(equality$summary$eta_critical))
    expect_true(is.na(equality$summary$critical_ratio))
    expect_match(equality$interpretation, "Every eta yields p = alpha")
    expect_error(plot_sensitivity(equality), "no unique critical eta boundary")
    expect_equal(sensitivity_grid(m, "Y", eta = c(0, 1), rho = 1)$grid$p_value, rep(p0, 2))
  }
  m <- contract_fixture(distance = 2)
  x <- critical_sensitivity(m, "Y", rho = 0)
  expect_lte(x$summary$eta_critical, 1)
  expect_true(is.finite(x$summary$critical_ratio))
  expect_true(is.na(x$curve$eta_critical))
  expect_error(plot_sensitivity(x), "All denominators on the requested rho grid are zero")
  equality <- critical_sensitivity(m, "Y", rho = 0, alpha = RI_unadjusted(m, "Y")$result$p_value)
  expect_identical(equality$summary$critical_ratio, 0)
  expect_error(plot_sensitivity(equality), "Every eta yields p = alpha")
})

test_that("relevant zero-bound pairs remain in the critical ratio denominator", {
  components <- list(bound = matrix(c(1, 0.5, 0, 0.5, 1, 0, 0, 0, 1), 3),
                     set_dist = matrix(1, 3, 3), kappa = 1)
  diag(components$set_dist) <- 0
  x <- netmatchRI:::.critical_summary(components, rep(1, 3), 3, 0.5)
  expect_equal(x$detail$n_relevant_pairs, 3)
  expect_equal(x$summary$eta_critical, 0.5)
  expect_equal(x$summary$critical_ratio, (0.5 / (2 * 3)) / (3 / 3))
  components$bound[1, 2] <- components$bound[2, 1] <- 0
  for (numerator in c(-1, 0, 1)) {
    x <- netmatchRI:::.critical_summary(components, rep(1, 3), 3, numerator)
    expect_true(is.na(x$summary$eta_critical))
    expect_true(is.na(x$summary$critical_ratio))
  }
})

test_that("zero and non-finite diagonal variances are explicit, never undefined NaNs", {
  components <- list(bound = matrix(c(1, 0.5, 0.5, 1), 2),
                     set_dist = matrix(c(0, 1, 1, 0), 2), kappa = 1)
  for (v in c(0, Inf, NaN)) {
    x <- netmatchRI:::.critical_summary(components, c(1, 1), v, 0)
    expect_true(is.na(x$summary$eta_critical))
    expect_true(is.na(x$summary$critical_ratio))
    expect_match(netmatchRI:::.critical_interpretation(0, v, 1, 1, NA_real_), "zero or non-finite")
    expect_identical(netmatchRI:::.normal_pvalue(1, 0, v), list(z = NA_real_, p = NA_real_))
  }
  m <- contract_fixture(n_sets = 1L); m$data$Z[] <- 1
  x <- critical_sensitivity(m, "Y")
  expect_true(is.na(x$summary$eta_critical))
  expect_error(plot_sensitivity(x), "diagonal variance is zero or non-finite")
})

test_that("inference entry points retain strict comparisons without an outcome-ties guard", {
  functions <- list(RI_unadjusted, RI_adjusted, RI_design, sensitivity_grid, critical_sensitivity)
  m <- contract_fixture(full = TRUE) # repeated values occur across sets only
  for (fn in functions) expect_silent(fn(m, "Y"))
  for (ids in list(c(1, 2), c(4, 5), c(6, 7))) { # across arms, controls, treated
    tied <- m; tied$data$Y[ids] <- 99
    # Execution contract only: this does not establish validity with ties.
    for (fn in functions) expect_silent(fn(tied, "Y"))
  }
  expect_equal(netmatchRI:::.mw_u(c(1, 2), c(1, 2)), 1)
})

test_that("normal tails and tiny alpha avoid cancellation", {
  p <- netmatchRI:::.normal_pvalue(10, 0, 1)$p
  expect_gt(p, 0)
  expect_equal(p, 2 * pnorm(10, lower.tail = FALSE), tolerance = 1e-14)
  m <- contract_fixture(n_sets = 100L)
  expect_equal(RI_unadjusted(m, "Y")$result$p_value, p, tolerance = 1e-12)
  x <- critical_sensitivity(m, "Y", alpha = 1e-20, rho = 1)
  expect_lte(x$summary$eta_critical, 1)
  expect_gt(x$summary$eta_critical, 0)
  got <- RI_adjusted(m, "Y", eta = x$summary$eta_critical, rho = 1)$result$p_value
  expect_equal(got / 1e-20, 1, tolerance = 1e-12)
  expect_silent(critical_sensitivity(contract_fixture(), "Y", alpha = 5e-324, rho = 1))
})

test_that("distance domain is enforced without forcing integer distances or unadjusted network work", {
  m <- contract_fixture()
  network_functions <- list(RI_adjusted, RI_design, sensitivity_grid, critical_sensitivity)
  for (value in c(0, 0.5, -1, -Inf, NA_real_, NaN)) {
    bad <- m; bad$network_distance[1, 3] <- bad$network_distance[3, 1] <- value
    expect_error(netmatchRI:::.validate_network(bad$network_distance, "distance"),
                 "off-diagonal distances >= 1|nonnegative|missing or NaN")
    for (fn in network_functions) expect_error(fn(bad, "Y"),
                 "off-diagonal distances >= 1|nonnegative|missing or NaN")
    expect_silent(RI_unadjusted(bad, "Y"))
  }
  for (value in c(1.5, Inf)) {
    good <- m; good$network_distance[1, 3] <- good$network_distance[3, 1] <- value
    expect_silent(netmatchRI:::.validate_network(good$network_distance, "distance"))
    for (fn in network_functions) expect_silent(fn(good, "Y"))
  }
  fractional <- contract_fixture(distance = 1.5)
  expect_equal(RI_adjusted(fractional, "Y", eta = 0, rho = 0)$result$variance,
               RI_unadjusted(fractional, "Y")$result$variance)
  expect_lte(critical_sensitivity(fractional, "Y", rho = 0)$summary$eta_critical, 1)
  bad <- m; diag(bad$network_distance) <- 1
  expect_error(RI_design(bad, "Y"), "zero diagonal")
  bad <- m; bad$network_distance[1, 3] <- 2
  expect_error(RI_design(bad, "Y"), "symmetric")
  m$network_distance <- NULL
  expect_silent(RI_unadjusted(m, "Y"))
  expect_error(RI_design(m, "Y"), "square numeric distance matrix")
  A <- matrix(0, 4, 4)
  expect_silent(netmatchRI:::.validate_network(A, "adjacency"))
  expect_true(all(is.infinite(netmatchRI:::.as_network_distance(A, "adjacency")[upper.tri(A)])))
})
