make_plot_design <- function(subclass, method) {
  data <- data.frame(
    Z = c(1, 1, 1, 1, 0, 0, 0, 0),
    X1 = c(0, 1, 2, 3, 0.2, 1.2, 2.4, 3.8),
    X2 = c(1, 2, 1, 2, 1.1, 1.8, 1.3, 2.5),
    Y = c(8, 7, 6, 5, 1, 2, 3, 4)
  )
  subclass <- factor(subclass)
  matched <- !is.na(subclass)
  matched_data_frame <- data[matched, , drop = FALSE]
  matched_data_frame$subclass <- droplevels(subclass[matched])
  rownames(matched_data_frame) <- which(matched)
  network_distance <- matrix(3, nrow(data), nrow(data))
  diag(network_distance) <- 0
  out <- list(
    data = matched_data_frame,
    treat = "Z",
    covariates = c("X1", "X2"),
    method = method,
    kappa = 2,
    solver = if (method == "dual") "highs" else "optmatch",
    network_distance = network_distance,
    subclass = subclass,
    weights = netmatchRI:::.matching_weights(data$Z, subclass, "ATT"),
    matched = matched,
    call = quote(netmatch()),
    estimand = "ATT",
    network_type = "distance",
    original_data = data,
    solver_info = list(
      backend = if (method == "dual") "highs" else "optmatch",
      status = "OPTIMAL", objective = NA_real_, gap = NA_real_,
      runtime = NA_real_
    )
  )
  class(out) <- "netmatch"
  out
}

plot_designs_fixture <- function() {
  list(
    "Dual-penalty" = make_plot_design(
      c(1, 2, 3, 4, 1, 2, 3, 4), "dual"
    ),
    "Covariate-only" = make_plot_design(
      c(1, 1, 2, 2, 1, 1, 2, 2), "covariate"
    ),
    "Single-penalty" = make_plot_design(
      c(1, 2, 2, 3, 1, 2, 2, 3), "single"
    )
  )
}

test_that("comparison plots validate a named common-design list", {
  designs <- plot_designs_fixture()
  expect_error(
    plot_covariate_balance(list()),
    "non-empty named list"
  )
  expect_error(
    plot_covariate_balance(unname(designs)),
    "unique non-empty names"
  )
  different <- designs
  different[[2]]$original_data$X1[1] <- 99
  expect_error(
    plot_covariate_similarity(different),
    "same original data"
  )
  expect_error(
    plot_covariate_similarity(designs, include_before = NA),
    "must be TRUE or FALSE"
  )
})

test_that("covariate-balance plot includes one Before level and selected designs", {
  skip_if_not_installed("ggplot2")
  designs <- plot_designs_fixture()[c("Dual-penalty", "Covariate-only")]
  plot <- plot_covariate_balance(designs, include_before = TRUE)

  expect_s3_class(plot, "ggplot")
  expect_equal(
    levels(plot$data$design),
    c("Before", "Dual-penalty", "Covariate-only")
  )
  expect_equal(nrow(plot$data), 3 * 2)
  expected <- diagnose_match(designs[["Dual-penalty"]])$covariate_balance
  expect_equal(
    plot$data$abs_smd[plot$data$design == "Before"],
    expected$before_abs_smd
  )
  expect_equal(plot$labels$x, "Absolute standardized mean difference")
  expect_equal(plot$labels$colour, "Design")
})

test_that("covariate-similarity plot uses all pairs before and within-matched-set pairs after", {
  skip_if_not_installed("ggplot2")
  designs <- plot_designs_fixture()[c("Dual-penalty", "Single-penalty")]
  plot <- plot_covariate_similarity(designs, include_before = TRUE)
  first <- designs[[1]]
  z <- first$original_data[[first$treat]]
  distance <- netmatchRI:::.mahalanobis_matrix(
    first$original_data, z, first$covariates, cov_type = "pooled"
  )
  expected_before <- mean(distance[is.finite(distance)])

  expect_s3_class(plot, "ggplot")
  expect_equal(
    levels(plot$data$design),
    c("Before", "Dual-penalty", "Single-penalty")
  )
  expect_equal(
    plot$data$average_mahalanobis[plot$data$design == "Before"],
    expected_before
  )

  dual <- designs[["Dual-penalty"]]
  treated <- which(z == 1)
  control <- which(z == 0)
  same_set <- outer(
    as.character(dual$subclass[treated]),
    as.character(dual$subclass[control]),
    FUN = "=="
  )
  expected_after <- mean(distance[same_set])
  expect_equal(
    plot$data$average_mahalanobis[plot$data$design == "Dual-penalty"],
    expected_after
  )
  expect_equal(
    plot$data$n_pairs[plot$data$design == "Dual-penalty"],
    sum(same_set)
  )
  expect_equal(
    plot$labels$y,
    "Average treated-control Mahalanobis distance"
  )
})

test_that("comparison plots can omit Before and preserve requested design order", {
  skip_if_not_installed("ggplot2")
  designs <- plot_designs_fixture()[c("Single-penalty", "Dual-penalty")]
  balance <- plot_covariate_balance(designs, include_before = FALSE)
  similarity <- plot_covariate_similarity(designs, include_before = FALSE)
  expect_equal(levels(balance$data$design), names(designs))
  expect_equal(levels(similarity$data$design), names(designs))
  expect_false("Before" %in% as.character(balance$data$design))
  expect_false("Before" %in% as.character(similarity$data$design))
})

test_that("comparison plots add no composition-package dependency", {
  code <- c(
    deparse(body(plot_covariate_balance)),
    deparse(body(plot_covariate_similarity))
  )
  expect_false(any(grepl("patchwork|cowplot|gridExtra", code)))
})
