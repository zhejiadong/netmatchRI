#' Simulate the 300-Unit Network-Matching Example
#'
#' Generates one data set from the built-in simulation design: a fixed
#' stochastic block network, network-dependent covariates, treatment, and
#' outcome. The default values reproduce the package illustration.
#'
#' @param seed Random seed for the simulation.
#' @param beta_z Treatment effect in the outcome model.
#' @param n Number of units. The built-in example uses 300.
#' @param alpha1 Dependence strength in `[0, 1]`; larger values induce stronger
#'   network dependence.
#' @param alpha2 Individual variation in `[0, 1]` independent of the network.
#' @param pin Within-block edge probability for the stochastic block network.
#' @param pout Between-block edge probability for the stochastic block network.
#' @return A list with the simulated `data`, `Adj` (an `n` by `n` adjacency
#'   matrix for the fixed network), graph-distance matrix `net_dist`, `V`
#'   (variance-covariance matrix), simulation `seed`, `beta_z`, `alpha1`,
#'   `alpha2`, `pin`, and `pout`.
#' @export
simulate_netmatch_example <- function(seed = 90141,
                                      beta_z = 0,
                                      n = 300,
                                      alpha1 = 0.9,
                                      alpha2 = 0.1,
                                      pin = 0.19,
                                      pout = 0.003) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package `igraph` is required for simulation.", call. = FALSE)
  }
  if (!requireNamespace("mvtnorm", quietly = TRUE)) {
    stop("Package `mvtnorm` is required for simulation.", call. = FALSE)
  }
  if (n %% 4 != 0) stop("`n` must be divisible by 4.", call. = FALSE)
  .check_unit_interval(alpha1, "alpha1")
  .check_unit_interval(alpha2, "alpha2")
  .check_unit_interval(pin, "pin")
  .check_unit_interval(pout, "pout")

  setup <- .strong_simulation_setup(n, pin = pin, pout = pout, alpha1 = alpha1, alpha2 = alpha2)
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
    beta_z = beta_z,
    alpha1 = alpha1,
    alpha2 = alpha2,
    pin = pin,
    pout = pout
  )
}

.check_unit_interval <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1 || !is.finite(x) || x < 0 || x > 1) {
    stop(sprintf("`%s` must be one number in [0, 1].", name), call. = FALSE)
  }
  invisible(TRUE)
}

.netmatch_sim_cache <- new.env(parent = emptyenv())

.strong_simulation_setup <- function(n, pin, pout, alpha1, alpha2) {
  key <- paste("strong", n, pin, pout, alpha1, alpha2, sep = "_")
  if (exists(key, envir = .netmatch_sim_cache, inherits = FALSE)) {
    return(get(key, envir = .netmatch_sim_cache, inherits = FALSE))
  }
  n_blocks <- 4
  block_sizes <- rep(n / n_blocks, n_blocks)
  pm <- matrix(pout, nrow = n_blocks, ncol = n_blocks)
  diag(pm) <- pin

  set.seed(2026)
  g <- igraph::sample_sbm(
    n,
    pref.matrix = pm,
    block.sizes = block_sizes,
    directed = FALSE
  )
  Adj <- as.matrix(igraph::as_adjacency_matrix(g))
  net_dist <- igraph::distances(g)

  W <- .make_W_row_standardized(Adj)
  V <- .scale_to_unit_diag(.generate_V_W_step2(
    W,
    alpha1 = alpha1,
    alpha2 = alpha2,
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

.generate_V_W_step2 <- function(W, alpha1, alpha2, addI = 0.001) {
  n <- nrow(W)
  B <- alpha2 * diag(n) + alpha1 * W
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
