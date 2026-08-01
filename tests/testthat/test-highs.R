highs_test_model <- function() {
  list(
    A = Matrix::Matrix(matrix(c(1, 1), nrow = 1), sparse = TRUE),
    obj = c(1, 2),
    sense = ">=",
    rhs = 1,
    vtype = rep("B", 2),
    modelsense = "min"
  )
}

highs_mock_result <- function(status = 7L,
                              x = c(1, 0),
                              objective = 1,
                              primal_status = "Feasible",
                              gap = 0,
                              message = "Optimal") {
  list(
    primal_solution = x,
    objective_value = objective,
    status = status,
    status_message = message,
    info = list(primal_solution_status = primal_status, mip_gap = gap)
  )
}

test_that("default and explicit HiGHS solve a small feasible design", {
  skip_if_not_installed("highs")
  D <- matrix(4, 4, 4)
  diag(D) <- 0
  dat <- data.frame(
    Z = c(1, 1, 0, 0),
    X1 = c(0, 10, 0.1, 10.1),
    X2 = c(0, 1, 0.2, 1.2)
  )

  default_fit <- netmatch(dat, "Z", c("X1", "X2"), D,
                          method = "dual", kappa = 2)
  explicit_fit <- netmatch(dat, "Z", c("X1", "X2"), D,
                           method = "dual", kappa = 2, solver = "highs")

  expect_equal(default_fit$solver, "highs")
  expect_equal(explicit_fit$solver, "highs")
  expect_equal(sort(as.integer(rownames(default_fit$data))), 1:4)
  expect_equal(sort(as.integer(rownames(explicit_fit$data))), 1:4)
  expect_equal(default_fit$solver_result$status, "OPTIMAL")
  expect_true(is.list(default_fit$solver_result$raw))
  expect_equal(default_fit$solver_result$metadata$status_code, 7L)
  expect_false(default_fit$solver_result$metadata$time_limited_incumbent)
})

test_that("solver resolution priority is independent of installed solvers", {
  local_mocked_bindings(
    .gurobi_available = function() TRUE,
    .highs_available = function() TRUE,
    .glpk_available = function() TRUE,
    .package = "netmatchRI"
  )
  expect_equal(netmatchRI:::.resolve_solver("auto"), "gurobi")
  expect_equal(
    netmatchRI:::.resolve_solver_candidates("auto"),
    c("gurobi", "highs", "glpk")
  )

  local_mocked_bindings(
    .gurobi_available = function() FALSE,
    .highs_available = function() TRUE,
    .glpk_available = function() TRUE,
    .package = "netmatchRI"
  )
  expect_equal(netmatchRI:::.resolve_solver("auto"), "highs")

  local_mocked_bindings(
    .gurobi_available = function() FALSE,
    .highs_available = function() FALSE,
    .glpk_available = function() TRUE,
    .package = "netmatchRI"
  )
  expect_equal(netmatchRI:::.resolve_solver("auto"), "glpk")
})

test_that("auto retries the next available backend after a runtime failure", {
  attempts <- character(0)
  fake_match_once <- function(data, ..., solver) {
    attempts <<- c(attempts, solver)
    if (solver == "highs") stop("simulated HiGHS runtime failure", call. = FALSE)
    list(data = data, solver_result = list(status = "OPTIMAL"))
  }
  local_mocked_bindings(
    .gurobi_match_once = fake_match_once,
    .package = "netmatchRI"
  )
  args <- list(
    data = data.frame(Z = c(1, 0)), z = c(1, 0), D = matrix(1),
    unit_dist = matrix(3, 2, 2), method = "dual", kappa = 2,
    min_controls = 0.01, max_controls = 100, caliper = 8,
    timelimit = 90, mipgap = 0.01, threads = 1
  )
  matched <- do.call(
    netmatchRI:::.solve_network_match,
    c(args, list(solver = c("highs", "glpk")))
  )
  expect_equal(attempts, c("highs", "glpk"))
  expect_equal(matched$solver, "glpk")

  attempts <- character(0)
  expect_error(
    do.call(netmatchRI:::.solve_network_match, c(args, list(solver = "highs"))),
    "^simulated HiGHS runtime failure$"
  )
  expect_equal(attempts, "highs")
})

test_that("missing explicit HiGHS errors and auto preserves GLPK fallback", {
  local_mocked_bindings(
    .highs_available = function() FALSE,
    .package = "netmatchRI"
  )
  expect_error(netmatchRI:::.resolve_solver("highs"), "Package `highs` is required")

  local_mocked_bindings(
    .gurobi_available = function() FALSE,
    .highs_available = function() FALSE,
    .glpk_available = function() TRUE,
    .package = "netmatchRI"
  )
  expect_equal(netmatchRI:::.resolve_solver("auto"), "glpk")
})

test_that("HiGHS row bounds convert every supported sense", {
  bounds <- netmatchRI:::.highs_row_bounds(c(">=", "<=", "="), c(1, 2, 3))
  expect_equal(bounds$lhs, c(1, -Inf, 3))
  expect_equal(bounds$rhs, c(Inf, 2, 3))
  expect_error(netmatchRI:::.highs_row_bounds(c(">=", "<"), c(1, 2)),
               "Unknown model sense")
})

test_that("HiGHS adapter passes sparse matrix, full types and bounds", {
  skip_if_not_installed("highs")
  captured <- NULL
  fake_solve <- function(Q = NULL, L, lower, upper, A = NULL, lhs = NULL,
                         rhs = NULL, types = rep.int(1L, length(L)),
                         maximum = FALSE, offset = 0, start = NULL,
                         control = list()) {
    captured <<- list(L = L, lower = lower, upper = upper, A = A, lhs = lhs,
                     rhs = rhs, types = types, maximum = maximum,
                     control = control)
    highs_mock_result()
  }
  local_mocked_bindings(highs_solve = fake_solve, .package = "highs")
  result <- netmatchRI:::.solve_highs(highs_test_model(), 12, 0.02, 3)

  expect_s4_class(captured$A, "dgCMatrix")
  expect_equal(captured$types, rep("I", 2))
  expect_equal(captured$lower, c(0, 0))
  expect_equal(captured$upper, c(1, 1))
  expect_equal(captured$lhs, 1)
  expect_equal(captured$rhs, Inf)
  expect_identical(captured$maximum, FALSE)
  expect_identical(captured$control$output_flag, FALSE)
  expect_identical(captured$control$presolve, "on")
  expect_equal(captured$control$time_limit, 12)
  expect_equal(captured$control$mip_rel_gap, 0.02)
  expect_equal(captured$control$threads, 3L)
  expect_equal(result$status, "OPTIMAL")
})

test_that("HiGHS statuses are normalized safely", {
  model <- highs_test_model()
  optimal <- netmatchRI:::.normalize_highs_result(highs_mock_result(), model, 0.5)
  infeasible <- netmatchRI:::.normalize_highs_result(
    highs_mock_result(status = 8L, x = NULL, objective = NA_real_,
                      primal_status = "None", message = "Infeasible"),
    model, 0.1
  )
  ambiguous <- netmatchRI:::.normalize_highs_result(
    highs_mock_result(status = 9L, x = NULL, objective = NA_real_,
                      primal_status = "None", message = "Unbounded or infeasible"),
    model, 0.1
  )
  unbounded <- netmatchRI:::.normalize_highs_result(
    highs_mock_result(status = 10L, x = NULL, objective = NA_real_,
                      primal_status = "None", message = "Unbounded"),
    model, 0.1
  )
  unknown <- netmatchRI:::.normalize_highs_result(
    highs_mock_result(status = 99L, message = "Future status"), model, 0.1
  )

  expect_equal(optimal$status, "OPTIMAL")
  expect_equal(optimal$objective, 1)
  expect_equal(infeasible$status, "INFEASIBLE")
  expect_null(infeasible$x)
  expect_equal(ambiguous$status, "INF_OR_UNBD")
  expect_null(ambiguous$x)
  expect_equal(unbounded$status, "UNBOUNDED")
  expect_null(unbounded$x)
  expect_equal(unknown$status, "UNKNOWN")
  expect_null(unknown$x)
})

test_that("optimal status without a feasible primal solution is rejected", {
  none <- highs_mock_result(
    status = 7L, x = c(0, 0), objective = 0,
    primal_status = "None", message = "Optimal"
  )
  result <- netmatchRI:::.normalize_highs_result(none, highs_test_model(), 0.1)
  expect_equal(result$status, "OPTIMAL_NO_SOLUTION")
  expect_null(result$x)
})

test_that("time limit requires and warns about a validated incumbent", {
  model <- highs_test_model()
  feasible <- highs_mock_result(status = 13L, primal_status = "Feasible",
                                gap = 0.25, message = "Time limit reached")
  expect_warning(
    result <- netmatchRI:::.normalize_highs_result(feasible, model, 2.5),
    "validated feasible incumbent.*actual relative MIP gap: 0.25"
  )
  expect_equal(result$status, "TIME_LIMIT")
  expect_true(result$metadata$time_limited_incumbent)
  expect_equal(result$metadata$gap, 0.25)
  expect_equal(result$metadata$run_time, 2.5)

  none <- highs_mock_result(status = 13L, x = c(0, 0), objective = 0,
                            primal_status = "None", gap = Inf,
                            message = "Time limit reached")
  no_incumbent <- netmatchRI:::.normalize_highs_result(none, model, 2.5)
  expect_equal(no_incumbent$status, "TIME_LIMIT_NO_INCUMBENT")
  expect_null(no_incumbent$x)
  expect_false(no_incumbent$metadata$time_limited_incumbent)
})

test_that("solution validator rejects malformed or invalid vectors", {
  model <- highs_test_model()
  validate <- netmatchRI:::.validate_mip_solution

  expect_error(validate(model, 1), "wrong length")
  expect_error(validate(model, c(1, NA_real_)), "non-finite")
  expect_error(validate(model, c(0.5, 0.5)), "integrality")
  expect_error(validate(model, c(1.1, 0)), "binary variable bounds")
  expect_error(validate(model, c(0, 0)), "model constraints")
  expect_error(validate(model, c(1, 0), reported_objective = 2),
               "objective does not match")

  bad_sense <- model
  bad_sense$sense <- "<"
  expect_error(validate(bad_sense, c(1, 0)), "Unknown model sense")
  expect_equal(validate(model, c(1, 0), reported_objective = 1)$objective, 1)
})

test_that("validator canonicalizes near-integer values before large row checks", {
  p <- 5000L
  model <- list(
    A = Matrix::Matrix(matrix(1, nrow = 1, ncol = p), sparse = TRUE),
    obj = rep(1, p), sense = ">=", rhs = p,
    vtype = rep("B", p), modelsense = "min"
  )
  checked <- netmatchRI:::.validate_mip_solution(model, rep(1 - 5e-7, p))
  expect_equal(checked$x, rep(1, p))
  expect_equal(checked$objective, p)
})

test_that("subclass extraction excludes units absent from selected edges", {
  E <- data.frame(e = 1:2, t = c(1L, 2L), c = c(1L, 2L), cost = c(1, 1))
  subclasses <- netmatchRI:::.extract_subclasses(
    list(x = c(1, 0)), treat_ids = c(1L, 2L), ctrl_ids = c(3L, 4L),
    E = E, n = 4L
  )
  expect_equal(subclasses[c(1, 3)], c(1L, 1L))
  expect_true(all(is.na(subclasses[c(2, 4)])))
})

test_that("HiGHS and GLPK agree on a deterministic proven optimum", {
  skip_if_not_installed("highs")
  skip_if_not_installed("Rglpk")
  skip_if_not_installed("slam")
  model <- highs_test_model()
  highs_result <- netmatchRI:::.solve_highs(model, 30, 0, 1)
  glpk_result <- netmatchRI:::.solve_roi_glpk(model, 30)
  expect_equal(highs_result$status, "OPTIMAL")
  expect_equal(glpk_result$status, "OPTIMAL")
  highs_objective <- netmatchRI:::.validate_mip_solution(model, highs_result$x)$objective
  glpk_objective <- netmatchRI:::.validate_mip_solution(model, glpk_result$x)$objective
  expect_equal(highs_objective, glpk_objective, tolerance = 1e-8)
})

test_that("available solvers agree on an original simulation smoke test", {
  sim <- simulate_netmatch_example(seed = 90141, n = 60)
  fit_with <- function(solver) {
    netmatch(
      data = sim$data, treat = "Z", covariates = c("X1", "X2", "X3"),
      network = sim$net_dist, method = "dual", kappa = 2,
      solver = solver, timelimit = 120, mipgap = 0, threads = 1
    )
  }

  fits <- list(highs = fit_with("highs"), glpk = fit_with("glpk"))
  objectives <- c(
    highs = fits$highs$solver_result$objective,
    glpk = fits$glpk$solver_result$raw$optimum
  )
  if (netmatchRI:::.gurobi_available()) {
    fits$gurobi <- fit_with("gurobi")
    objectives <- c(objectives, gurobi = fits$gurobi$solver_result$objval)
  }

  expect_true(all(vapply(fits, function(x) x$solver_result$status == "OPTIMAL", logical(1))))
  expect_equal(unname(objectives), rep(objectives[[1]], length(objectives)), tolerance = 1e-8)
  expect_true(all(vapply(fits, function(x) nrow(x$data) == nrow(sim$data), logical(1))))
  expect_true(all(vapply(fits, function(x) sum(x$solver_result$x > 0.5) == 47L, logical(1))))
})
