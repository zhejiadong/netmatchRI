.make_sbm_network <- function(n = 300, n_block = 4, p_in = 0.19, p_out = 0.003, seed = 2026) {
  set.seed(seed)
  block <- rep(seq_len(n_block), each = n / n_block)
  A <- matrix(0, n, n)
  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) {
      p <- if (block[i] == block[j]) p_in else p_out
      A[i, j] <- A[j, i] <- stats::rbinom(1, 1, p)
    }
  }
  A
}

.row_standardize <- function(A) {
  rs <- rowSums(A)
  W <- matrix(0, nrow(A), ncol(A))
  idx <- rs > 0
  W[idx, ] <- A[idx, , drop = FALSE] / rs[idx]
  W
}

.make_cov_matrix <- function(W, alpha1, alpha2) {
  B <- alpha1 * W + alpha2 * diag(nrow(W))
  V <- B %*% t(B)
  d <- sqrt(diag(V))
  V / (d %o% d)
}

.mvnorm_one <- function(mu, Sigma) {
  L <- tryCatch(chol(Sigma), error = function(e) NULL)
  z <- stats::rnorm(length(mu))
  if (!is.null(L)) return(as.numeric(mu + crossprod(L, z)))
  ev <- eigen(Sigma, symmetric = TRUE)
  as.numeric(mu + ev$vectors %*% (sqrt(pmax(ev$values, 0)) * z))
}

.solve_threshold <- function(mu, V, prevalence) {
  sd_i <- sqrt(diag(V))
  f <- function(cut) mean(1 - stats::pnorm((cut - mu) / sd_i)) - prevalence
  stats::uniroot(f, c(-60, 60))$root
}

.gen_paper_data <- function(V, beta_z = 0, prevalence = 0.3, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(V)
  X1 <- .mvnorm_one(rep(0, n), V)
  X2 <- .mvnorm_one(rep(0, n), V)
  X3 <- stats::rnorm(n)
  X <- cbind(X1, X2, X3)
  gamma_x <- c(0.1, 0.1, 0.1)
  beta_x <- c(0.1, 0.1, 0.1)
  mu_z <- as.numeric(X %*% gamma_x)
  cut <- .solve_threshold(mu_z, V, prevalence)
  z_latent <- .mvnorm_one(mu_z, V)
  Z <- as.integer(z_latent >= cut)
  Y <- as.numeric(X %*% beta_x + Z * beta_z + .mvnorm_one(rep(0, n), V))
  data.frame(Y = Y, Z = Z, X1 = X1, X2 = X2, X3 = X3)
}

#' Reproduce Paper Simulation Results
#'
#' Runs paper-reproduction tasks. Use `fast = TRUE` for installation smoke tests;
#' full manuscript-scale runs use `n_rep = 500` and can take a long time because
#' matching is repeated many times.
#'
#' @param task One of `"table1"`, `"ri_grid"`, `"matching_diagnostics"`,
#'   `"power"`, or `"all"`.
#' @param n_rep Number of replications. Defaults to 500 unless `fast = TRUE`.
#' @param fast If `TRUE`, uses a small network and two replications.
#' @param output_dir Directory where CSV outputs are written.
#' @param seed Random seed.
#' @return A list with result data frames.
#' @export
reproduce_paper <- function(task = c("table1", "ri_grid", "matching_diagnostics", "power", "all"),
                            n_rep = if (fast) 2 else 500,
                            fast = FALSE,
                            output_dir = file.path(getwd(), "netmatchRI-paper-output"),
                            seed = 1) {
  task <- match.arg(task)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  n <- if (fast) 40 else 300
  n_block <- 4
  A <- .make_sbm_network(n = n, n_block = n_block, seed = 2026)
  Dnet <- .as_network_distance(A)
  W <- .row_standardize(A)
  dep <- list(
    independent = c(0, 1),
    moderate = c(0.7, 0.3),
    strong = c(0.9, 0.1)
  )
  eta_rho <- data.frame(
    dependence = names(dep),
    eta = c(0, 0.01, 0.03),
    rho = c(0, 0.05, 0.10),
    stringsAsFactors = FALSE
  )

  tasks <- if (task == "all") c("table1", "ri_grid", "matching_diagnostics", "power") else task
  out <- list()
  if ("table1" %in% tasks) {
    out$table1 <- .run_table1(n_rep, dep, eta_rho, W, Dnet, seed, output_dir)
  }
  if ("ri_grid" %in% tasks) {
    out$ri_grid <- .run_ri_grid(max(1, n_rep), dep, W, Dnet, seed, output_dir)
  }
  if ("matching_diagnostics" %in% tasks) {
    out$matching_diagnostics <- .run_matching_diagnostics(max(1, n_rep), dep, W, Dnet, seed, output_dir)
  }
  if ("power" %in% tasks) {
    out$power <- .run_power(max(1, n_rep), dep, eta_rho, W, Dnet, seed, output_dir)
  }
  out
}

.run_table1 <- function(n_rep, dep, eta_rho, W, Dnet, seed, output_dir) {
  rows <- list()
  methods <- c("covariate", "single", "dual")
  for (dname in names(dep)) {
    V <- .make_cov_matrix(W, dep[[dname]][1], dep[[dname]][2])
    er <- eta_rho[eta_rho$dependence == dname, ]
    for (r in seq_len(n_rep)) {
      df <- .gen_paper_data(V, beta_z = 0, seed = seed + r)
      for (m in methods) {
        mt <- tryCatch(netmatch(df, "Z", c("X1", "X2", "X3"), Dnet, method = m), error = function(e) NULL)
        if (is.null(mt)) next
        naive <- netmatch_test(mt, "Y", method = "naive")$result
        decay <- netmatch_test(mt, "Y", method = "decay", eta = er$eta, rho = er$rho, d0 = 2)$result
        lm_p <- tryCatch(summary(stats::lm(Y ~ Z + factor(subclass), data = mt$data))$coefficients["Z", "Pr(>|t|)"], error = function(e) NA_real_)
        rows[[length(rows) + 1]] <- data.frame(
          dependence = dname, rep = r, match_method = m,
          inference = c("LM_Subclass", "RI_Naive", "RI_Decay"),
          p_value = c(lm_p, naive$p_value, decay$p_value)
        )
      }
      rows[[length(rows) + 1]] <- data.frame(
        dependence = dname, rep = r, match_method = "All",
        inference = "LM_Crude",
        p_value = tryCatch(summary(stats::lm(Y ~ Z, data = df))$coefficients["Z", "Pr(>|t|)"], error = function(e) NA_real_)
      )
    }
  }
  res <- do.call(rbind, rows)
  utils::write.csv(res, file.path(output_dir, "table1_raw.csv"), row.names = FALSE)
  res
}

.run_ri_grid <- function(n_rep, dep, W, Dnet, seed, output_dir) {
  grid <- expand.grid(eta = seq(0, 0.05, by = 0.01), rho = seq(0, 0.20, by = 0.05))
  rows <- list()
  for (dname in names(dep)) {
    V <- .make_cov_matrix(W, dep[[dname]][1], dep[[dname]][2])
    for (r in seq_len(n_rep)) {
      df <- .gen_paper_data(V, beta_z = 0, seed = seed + r)
      mt <- tryCatch(netmatch(df, "Z", c("X1", "X2", "X3"), Dnet, method = "dual"), error = function(e) NULL)
      if (is.null(mt)) next
      sens <- netmatch_sensitivity(mt, "Y", eta = grid$eta, rho = grid$rho, d0 = 2)$grid
      sens$dependence <- dname
      sens$rep <- r
      rows[[length(rows) + 1]] <- sens
    }
  }
  res <- do.call(rbind, rows)
  utils::write.csv(res, file.path(output_dir, "ri_grid_raw.csv"), row.names = FALSE)
  res
}

.run_matching_diagnostics <- function(n_rep, dep, W, Dnet, seed, output_dir) {
  rows <- list()
  for (dname in names(dep)) {
    V <- .make_cov_matrix(W, dep[[dname]][1], dep[[dname]][2])
    for (r in seq_len(n_rep)) {
      df <- .gen_paper_data(V, beta_z = 0, seed = seed + r)
      for (m in c("covariate", "single", "dual")) {
        mt <- tryCatch(netmatch(df, "Z", c("X1", "X2", "X3"), Dnet, method = m), error = function(e) NULL)
        if (is.null(mt)) next
        dg <- diagnose_match(mt)
        rows[[length(rows) + 1]] <- data.frame(
          dependence = dname,
          rep = r,
          match_method = m,
          average_within_distance = dg$average_within_distance,
          mean_abs_smd = mean(dg$covariate_smd$abs_smd, na.rm = TRUE)
        )
      }
    }
  }
  res <- do.call(rbind, rows)
  utils::write.csv(res, file.path(output_dir, "matching_diagnostics_raw.csv"), row.names = FALSE)
  res
}

.run_power <- function(n_rep, dep, eta_rho, W, Dnet, seed, output_dir) {
  rows <- list()
  for (dname in names(dep)) {
    V <- .make_cov_matrix(W, dep[[dname]][1], dep[[dname]][2])
    er <- eta_rho[eta_rho$dependence == dname, ]
    for (bz in c(0, 0.1, 0.2, 0.3, 0.4)) {
      for (r in seq_len(n_rep)) {
        df <- .gen_paper_data(V, beta_z = bz, seed = seed + r + round(1000 * bz))
        mt <- tryCatch(netmatch(df, "Z", c("X1", "X2", "X3"), Dnet, method = "dual"), error = function(e) NULL)
        if (is.null(mt)) next
        naive <- netmatch_test(mt, "Y", method = "naive")$result
        decay <- netmatch_test(mt, "Y", method = "decay", eta = er$eta, rho = er$rho, d0 = 2)$result
        rows[[length(rows) + 1]] <- data.frame(
          dependence = dname, beta_z = bz, rep = r,
          inference = c("RI_Naive", "RI_Decay"),
          p_value = c(naive$p_value, decay$p_value)
        )
      }
    }
  }
  res <- do.call(rbind, rows)
  utils::write.csv(res, file.path(output_dir, "power_raw.csv"), row.names = FALSE)
  res
}

#' Summarize Table 1 Rejection Rates
#' @param raw Output from `reproduce_paper("table1")$table1`.
#' @param alpha Test level.
#' @export
paper_table1 <- function(raw, alpha = 0.05) {
  raw$reject <- raw$p_value < alpha
  stats::aggregate(reject ~ dependence + match_method + inference, raw, mean, na.rm = TRUE)
}

#' Plot Paper Power Results
#' @param raw Output from `reproduce_paper("power")$power`.
#' @param alpha Test level.
#' @export
paper_power_plot <- function(raw, alpha = 0.05) {
  raw$reject <- raw$p_value < alpha
  tab <- stats::aggregate(reject ~ dependence + beta_z + inference, raw, mean, na.rm = TRUE)
  names(tab)[names(tab) == "reject"] <- "rejection_rate"
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    ggplot2::ggplot(tab, ggplot2::aes_string("beta_z", "rejection_rate", color = "inference")) +
      ggplot2::geom_line() +
      ggplot2::geom_point() +
      ggplot2::facet_wrap(~ dependence) +
      ggplot2::labs(x = "Treatment effect", y = "Rejection rate")
  } else {
    plot(tab$beta_z, tab[[4]], xlab = "Treatment effect", ylab = "Rejection rate")
  }
}

#' Plot Paper RI Grid Results
#' @param raw Output from `reproduce_paper("ri_grid")$ri_grid`.
#' @param alpha Test level.
#' @export
paper_grid_plot <- function(raw, alpha = 0.05) {
  raw$reject <- raw$p_value < alpha
  tab <- stats::aggregate(reject ~ dependence + eta + rho, raw, mean, na.rm = TRUE)
  names(tab)[names(tab) == "reject"] <- "rejection_rate"
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    ggplot2::ggplot(tab, ggplot2::aes_string("rho", "rejection_rate", color = "factor(eta)")) +
      ggplot2::geom_line() +
      ggplot2::facet_wrap(~ dependence) +
      ggplot2::labs(color = "eta", y = "Rejection rate")
  } else {
    plot(tab$rho, tab[[4]], xlab = "rho", ylab = "Rejection rate")
  }
}
