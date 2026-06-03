#' Run the TRIP Real-Data Application
#'
#' Convenience wrapper for the TRIP real-data illustration. The function expects
#' a local `.RData` file containing `ds_use` and `Adj`, matching the current TRIP
#' workspace. It does not bundle or expose the raw TRIP data.
#'
#' @param data_file Path to `TRIP/data_clear.RData` or a compatible file.
#' @param output_dir Directory where application CSV outputs are written.
#' @param methods Matching methods to run.
#' @param kappa Network-distance threshold.
#' @param eta_grid Eta values for sensitivity analysis.
#' @param rho_grid Rho values for sensitivity analysis.
#' @param d0 Maximum set distance for model-assisted RI.
#' @return A named list with matched objects and result tables.
#' @export
trip_application <- function(data_file,
                             output_dir = file.path(getwd(), "trip-application-output"),
                             methods = c("covariate", "single", "dual"),
                             kappa = 2,
                             eta_grid = c(0.01, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40),
                             rho_grid = c(0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60),
                             d0 = 2) {
  if (!file.exists(data_file)) {
    stop("`data_file` does not exist: ", data_file, call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  env <- new.env(parent = emptyenv())
  load(data_file, envir = env)
  if (!exists("ds_use", envir = env) || !exists("Adj", envir = env)) {
    stop("`data_file` must contain objects `ds_use` and `Adj`.", call. = FALSE)
  }

  dat <- get("ds_use", envir = env)
  adj <- get("Adj", envir = env)
  needed <- c("EGO_ID", "Z", "Y", "edu", "employment", "ACMDT", "baseRisk", "HIV")
  missing <- setdiff(needed, names(dat))
  if (length(missing)) {
    stop("TRIP data is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  dat <- dat[, needed, drop = FALSE]
  dat$Z <- as.integer(dat$Z)
  dat$Y <- as.numeric(dat$Y)
  covariates <- c("edu", "employment", "ACMDT", "baseRisk", "HIV")

  adj <- .align_trip_network(adj, dat$EGO_ID)
  matches <- list()
  tests <- list()
  diagnostics <- list()
  sensitivities <- list()

  for (method in methods) {
    fit <- netmatch(
      data = dat,
      treat = "Z",
      covariates = covariates,
      network = adj,
      method = method,
      kappa = kappa
    )
    matches[[method]] <- fit

    naive <- netmatch_test(fit, "Y", method = "naive")
    decay <- netmatch_test(fit, "Y", method = "decay", eta = 0.03, rho = 0.10, d0 = d0)
    tmp <- rbind(naive$result, decay$result)
    tmp$match_method <- method
    tests[[method]] <- tmp

    dg <- diagnose_match(fit)
    diagnostics[[method]] <- data.frame(
      match_method = method,
      matched_units = nrow(fit$data),
      matched_sets = length(unique(fit$data$subclass)),
      average_within_distance = dg$average_within_distance,
      mean_abs_smd = mean(dg$covariate_smd$abs_smd, na.rm = TRUE)
    )

    sens <- netmatch_sensitivity(fit, "Y", eta = eta_grid, rho = rho_grid, d0 = d0)
    sens_grid <- sens$grid
    sens_grid$match_method <- method
    sensitivities[[method]] <- sens_grid
  }

  test_tbl <- do.call(rbind, tests)
  diag_tbl <- do.call(rbind, diagnostics)
  sens_tbl <- do.call(rbind, sensitivities)

  utils::write.csv(test_tbl, file.path(output_dir, "trip_tests.csv"), row.names = FALSE)
  utils::write.csv(diag_tbl, file.path(output_dir, "trip_diagnostics.csv"), row.names = FALSE)
  utils::write.csv(sens_tbl, file.path(output_dir, "trip_sensitivity.csv"), row.names = FALSE)

  list(
    matches = matches,
    tests = test_tbl,
    diagnostics = diag_tbl,
    sensitivity = sens_tbl,
    output_dir = output_dir
  )
}

.align_trip_network <- function(adj, ego_id) {
  adj <- as.matrix(adj)
  if (!is.null(rownames(adj)) && all(as.character(ego_id) %in% rownames(adj))) {
    adj <- adj[as.character(ego_id), as.character(ego_id), drop = FALSE]
  }
  if (nrow(adj) != length(ego_id)) {
    stop("TRIP adjacency matrix size does not match `ds_use`.", call. = FALSE)
  }
  rownames(adj) <- colnames(adj) <- seq_along(ego_id)
  adj
}
