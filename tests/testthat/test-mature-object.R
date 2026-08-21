make_mature_fixture <- function() {
  data <- data.frame(
    Z = c(1, 1, 1, 0, 0, 0),
    X = c(0, 1, 2, 0.1, 1.1, 2.1),
    F = factor(c("a", "b", "a", "a", "b", "b")),
    Y = 1:6,
    row.names = paste0("u", 1:6)
  )
  network <- matrix(4, 6, 6, dimnames = list(rev(rownames(data)), rev(rownames(data))))
  diag(network) <- 0
  list(data = data, network = network)
}

test_that("new netmatch arguments are appended to the positional API", {
  expect_identical(
    tail(names(formals(netmatch)), 3),
    c("estimand", "network_type", "include_solver")
  )
})

test_that("full-matching weights follow ATT, ATC, and ATE definitions", {
  z <- c(1, 1, 0, 1, 0, 0)
  subclass <- factor(c(1, 1, 1, 2, 2, 2))
  expect_equal(netmatchRI:::.matching_weights(z, subclass, "ATT"),
               c(1, 1, 2, 1, 0.5, 0.5))
  expect_equal(netmatchRI:::.matching_weights(z, subclass, "ATC"),
               c(0.5, 0.5, 1, 2, 1, 1))
  expect_equal(netmatchRI:::.matching_weights(z, subclass, "ATE"),
               c(0.75, 0.75, 1.5, 1.5, 0.75, 0.75))
})

test_that("netmatch validates and aligns mature public inputs", {
  fx <- make_mature_fixture()
  m <- netmatch(fx$data, "Z", c("X", "F"), fx$network,
                method = "covariate", estimand = "ATT", network_type = "distance")
  expect_identical(m$network_type, "distance")
  expect_identical(rownames(m$network_distance), rownames(fx$data))
  expect_identical(colnames(m$network_distance), rownames(fx$data))

  bad <- fx$data
  bad$Z[1] <- NA
  expect_error(netmatch(bad, "Z", c("X", "F"), fx$network,
                        method = "covariate", network_type = "distance"), "missing")
  bad <- fx$data
  bad$Z[1] <- 2
  expect_error(netmatch(bad, "Z", c("X", "F"), fx$network,
                        method = "covariate", network_type = "distance"), "0/1")
  bad <- fx$data
  bad$X[1] <- NA
  expect_error(netmatch(bad, "Z", c("X", "F"), fx$network,
                        method = "covariate", network_type = "distance"), "missing")

  nonsymmetric <- unname(fx$network)
  nonsymmetric[1, 2] <- 3
  expect_error(netmatch(fx$data, "Z", c("X", "F"), nonsymmetric,
                        method = "covariate", network_type = "distance"), "symmetric")
  nonzero_diag <- unname(fx$network)
  diag(nonzero_diag) <- 1
  expect_error(netmatch(fx$data, "Z", c("X", "F"), nonzero_diag,
                        method = "covariate", network_type = "distance"), "zero diagonal")
  negative <- unname(fx$network)
  negative[1, 2] <- negative[2, 1] <- -1
  expect_error(netmatch(fx$data, "Z", c("X", "F"), negative,
                        method = "covariate", network_type = "distance"), "nonnegative")
  disconnected <- unname(fx$network)
  disconnected[1, 2] <- disconnected[2, 1] <- Inf
  expect_error(netmatch(fx$data, "Z", c("X", "F"), disconnected,
                        method = "covariate", network_type = "adjacency"), "Inf")
  expect_s3_class(netmatch(fx$data, "Z", c("X", "F"), disconnected,
                           method = "covariate", network_type = "distance"), "netmatch")
})

test_that("matched objects expose full-length membership and normalized estimand weights", {
  fx <- make_mature_fixture()
  for (estimand in c("ATT", "ATC", "ATE")) {
    m <- netmatch(fx$data, "Z", c("X", "F"), fx$network,
                  method = "covariate", estimand = estimand,
                  network_type = "distance")
    expect_length(m$subclass, nrow(fx$data))
    expect_length(m$weights, nrow(fx$data))
    expect_length(m$matched, nrow(fx$data))
    expect_identical(m$matched, !is.na(m$subclass))
    expect_true(all(m$weights[!m$matched] == 0))
    expect_equal(mean(m$weights[m$matched & fx$data$Z == 1]), 1)
    expect_equal(mean(m$weights[m$matched & fx$data$Z == 0]), 1)
    expect_identical(m$estimand, estimand)
    expect_identical(m$original_data, fx$data)
    expect_true(inherits(m$call, "call"))
    expect_named(m$solver_info, c("backend", "status", "objective", "gap", "runtime"))
    expect_false("solver_result" %in% names(m))

    md <- matched_data(m)
    expect_equal(nrow(md), sum(m$matched))
    expect_true(all(c("subclass", "weights") %in% names(md)))
    expect_equal(nrow(matched_data(m, drop_unmatched = FALSE)), nrow(fx$data))
  }

  with_solver <- netmatch(fx$data, "Z", c("X", "F"), fx$network,
                          method = "covariate", network_type = "distance",
                          include_solver = TRUE)
  expect_true("solver_result" %in% names(with_solver))
})

test_that("diagnostics expand factors and use clear field names", {
  fx <- make_mature_fixture()
  m <- netmatch(fx$data, "Z", c("X", "F"), fx$network,
                method = "covariate", network_type = "distance")
  d <- diagnose_match(m)
  expect_named(d, c("covariate_balance", "network_summary", "within_distance_table"))
  expect_named(d$covariate_balance,
               c("covariate", "before_abs_smd", "after_abs_smd"))
  expect_true(all(c("F=a", "F=b") %in% d$covariate_balance$covariate))
  expect_true(all(c("n_pairs", "n_finite_pairs", "n_disconnected_pairs",
                    "min_distance", "mean_distance", "max_distance") %in%
                  names(d$network_summary)))
})

test_that("diagnostics use indicators for every unordered and ordered factor level", {
  fx <- make_mature_fixture()
  fx$data$F <- factor(c("a", "b", "c", "a", "b", "c"))
  fx$data$OF <- ordered(c("low", "mid", "high", "low", "mid", "high"),
                        levels = c("low", "mid", "high"))
  m <- netmatch(fx$data, "Z", c("X", "F", "OF"), fx$network,
                method = "covariate", network_type = "distance")
  balance_names <- diagnose_match(m)$covariate_balance$covariate
  expect_true(all(c("F=a", "F=b", "F=c") %in% balance_names))
  expect_true(all(c("OF=low", "OF=mid", "OF=high") %in% balance_names))
  expect_false(any(grepl("\\.L$|\\.Q$", balance_names)))
})

test_that("reserved matching column names are rejected clearly", {
  fx <- make_mature_fixture()
  for (reserved in c("subclass", "weights")) {
    bad <- fx$data
    bad[[reserved]] <- 1
    expect_error(
      netmatch(bad, "Z", c("X", "F"), fx$network,
               method = "covariate", network_type = "distance"),
      "reserved columns"
    )
  }
})

test_that("diagnostics retain disconnected pair counts with no finite distances", {
  fx <- make_mature_fixture()
  m <- netmatch(fx$data, "Z", c("X", "F"), fx$network,
                method = "covariate", network_type = "distance")
  m$network_distance[,] <- Inf
  diag(m$network_distance) <- 0
  d <- diagnose_match(m)
  expect_gt(d$network_summary$n_disconnected_pairs, 0)
  expect_equal(d$network_summary$n_finite_pairs, 0)
  expect_true(is.na(d$network_summary$mean_distance))
})

test_that("summary is a stable netmatch summary", {
  fx <- make_mature_fixture()
  m <- netmatch(fx$data, "Z", c("X", "F"), fx$network,
                method = "covariate", network_type = "distance")
  s <- summary(m)
  expect_s3_class(s, "netmatch_summary")
  expect_named(s, c("sample_counts", "ess", "set_sizes", "covariate_balance",
                    "network_summary", "solver_info"))
  expect_output(print(s), "Effective sample size")
  expect_output(print(m), "ATT")
})

test_that("randomization inference rejects missing and non-finite outcomes", {
  fx <- make_mature_fixture()
  m <- netmatch(fx$data, "Z", c("X", "F"), fx$network,
                method = "covariate", network_type = "distance")
  m$data$Y[1] <- NA_real_
  expect_error(RI_Naive(m, "Y"), "finite and non-missing")
  m$data$Y[1] <- Inf
  expect_error(RI_Naive(m, "Y"), "finite and non-missing")
})
