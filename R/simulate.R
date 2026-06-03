#' Simulate the 300-Unit Network-Matching Example
#'
#' Generates one data set from the built-in strong-dependence simulation design:
#' a fixed 300-unit stochastic block network,
#' network-dependent covariates, treatment, and outcome. This helper is
#' intended for package evaluation and examples, especially when searching for
#' seeds that illustrate the sensitivity analysis.
#'
#' @param seed Random seed for the covariates, treatment, and outcome.
#' @param beta_z Treatment effect in the outcome model.
#' @param n Number of units. The built-in example uses 300.
#' @return A list with the simulated `data`, fixed adjacency matrix `Adj`,
#'   graph-distance matrix `net_dist`, strong-dependence covariance matrix
#'   `V`, simulation `seed`, fixed `dep_index = 3`, and `beta_z`.
#' @export
simulate_netmatch_example <- function(seed = 90141,
                                      beta_z = 0,
                                      n = 300) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package `igraph` is required for simulation.", call. = FALSE)
  }
  if (!requireNamespace("mvtnorm", quietly = TRUE)) {
    stop("Package `mvtnorm` is required for simulation.", call. = FALSE)
  }
  if (n %% 4 != 0) stop("`n` must be divisible by 4.", call. = FALSE)

  dep_index <- 3
  setup <- .strong_simulation_setup(n)
  V <- setup$V

  beta_x <- matrix(c(0.1, 0.1, 0.1), nrow = 3, ncol = 1)
  gamma_x <- matrix(c(0.1, 0.1, 0.1), nrow = 3, ncol = 1)
  alpha_vec <- c(1, 1, 0)
  dat <- .simulate_process(
    V = V,
    beta_x = beta_x,
    gamma_x = gamma_x,
    beta_z = beta_z,
    alpha_vec = alpha_vec,
    seed = seed,
    target_treat = 0.23,
    K_min = 0.2 * n,
    K_max = 0.3 * n,
    max_try = 100
  )

  list(
    data = dat,
    Adj = setup$Adj,
    net_dist = setup$net_dist,
    V = V,
    seed = seed,
    dep_index = dep_index,
    beta_z = beta_z
  )
}

.netmatch_sim_cache <- new.env(parent = emptyenv())

.strong_simulation_setup <- function(n) {
  key <- paste0("strong_", n)
  if (exists(key, envir = .netmatch_sim_cache, inherits = FALSE)) {
    return(get(key, envir = .netmatch_sim_cache, inherits = FALSE))
  }
  n_blocks <- 4
  block_sizes <- rep(n / n_blocks, n_blocks)
  p_in <- 0.19
  p_out <- 0.003
  pm <- matrix(p_out, nrow = n_blocks, ncol = n_blocks)
  diag(pm) <- p_in

  set.seed(2026)
  g <- igraph::sample_sbm(
    n,
    pref.matrix = pm,
    block.sizes = block_sizes,
    directed = FALSE
  )
  Adj <- as.matrix(igraph::as_adjacency_matrix(g))
  net_dist <- igraph::distances(g)

  kappa_strong <- 0.9
  alpha_strong <- 0.1
  W <- .make_W_row_standardized(Adj)
  V <- .scale_to_unit_diag(.generate_V_W_step2(
    W,
    kappa = kappa_strong,
    alpha = alpha_strong,
    addI = 0.001
  ))
  setup <- list(Adj = Adj, net_dist = net_dist, V = V)
  assign(key, setup, envir = .netmatch_sim_cache)
  setup
}

.make_W_row_standardized <- function(Adj) {
  Adj <- as.matrix(Adj)
  diag(Adj) <- 0
  rs <- rowSums(Adj)
  W <- matrix(0, nrow(Adj), ncol(Adj))
  idx <- which(rs > 0)
  W[idx, ] <- Adj[idx, , drop = FALSE] / rs[idx]
  W
}

.scale_to_unit_diag <- function(V) {
  d <- sqrt(diag(V))
  V / (d %o% d)
}

.generate_V_W_step2 <- function(W, kappa, alpha, addI = 0.001) {
  n <- nrow(W)
  B <- alpha * diag(n) + kappa * W
  B %*% t(B) + addI * diag(n)
}

.solve_c_global <- function(mu, V, p) {
  sd_i <- sqrt(diag(V))
  g <- function(cut) mean(1 - stats::pnorm((cut - mu) / sd_i)) - p
  stats::uniroot(g, c(-60, 60))$root
}

.simulate_process <- function(V,
                              beta_x,
                              gamma_x,
                              beta_z,
                              alpha_vec,
                              seed,
                              target_treat,
                              K_min,
                              K_max,
                              max_try) {
  n <- nrow(V)
  set.seed(seed)

  X <- matrix(0, nrow = n, ncol = length(alpha_vec))
  for (j in seq_along(alpha_vec)) {
    alpha <- alpha_vec[j]
    Sigma_X <- alpha * V + (1 - alpha) * diag(n)
    X[, j] <- as.numeric(mvtnorm::rmvnorm(1, mean = rep(0, n), sigma = Sigma_X))
  }

  mu <- as.numeric(X %*% gamma_x)
  c_thr <- .solve_c_global(mu, V, target_treat)

  n1 <- 0
  for (attempt in seq_len(max_try)) {
    Z_lin <- as.numeric(mvtnorm::rmvnorm(1, mean = mu, sigma = V))
    Z <- as.integer(Z_lin > c_thr)
    n1 <- sum(Z)
    if (n1 >= K_min && n1 <= K_max) break
  }
  if (n1 < K_min || n1 > K_max) {
    warning(sprintf(
      "Reached %d attempts but only %d treated units (min=%d, max=%d).",
      max_try, n1, K_min, K_max
    ), call. = FALSE)
  }

  Y <- as.numeric(mvtnorm::rmvnorm(
    1,
    mean = as.numeric(X %*% beta_x + Z * beta_z),
    sigma = V
  ))
  data.frame(Y = Y, Z = Z, X1 = X[, 1], X2 = X[, 2], X3 = X[, 3])
}
