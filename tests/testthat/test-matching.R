test_that("single penalty uses direct network distance threshold", {
  testthat::skip_if_not_installed("optmatch")
  D <- matrix(3, 6, 6)
  diag(D) <- 0
  D[1, 4] <- D[4, 1] <- 1
  dat <- data.frame(
    Z = c(1, 1, 1, 0, 0, 0),
    X1 = c(0, 10, 20, 0.1, 10.1, 20.1),
    X2 = c(0, 0, 0, 0, 0, 0)
  )
  m <- netmatch(dat, "Z", c("X1", "X2"), D, method = "single", kappa = 2)
  expect_equal(m$solver, "optmatch")
  expect_false(any(rownames(m$data)[m$data$subclass == m$data$subclass[rownames(m$data) == "1"]] == "4"))
})

test_that("dual penalty keeps full matching constraints", {
  skip_without_gurobi()
  D <- matrix(4, 5, 5)
  diag(D) <- 0
  D[3, 4] <- D[4, 3] <- 1
  D[3, 5] <- D[5, 3] <- 1
  D[4, 5] <- D[5, 4] <- 1
  dat <- data.frame(
    Z = c(1, 1, 0, 0, 0),
    X1 = c(0, 10, 0.1, 10.1, 5),
    X2 = c(0, 1, 0.2, 1.2, 0.5)
  )
  expect_error(
    netmatch(dat, "Z", c("X1", "X2"), D, method = "dual", kappa = 2,
             solver = "gurobi", max_controls = 3),
    "infeasible or unbounded"
  )
})

test_that("dual penalty uses Gurobi and covers all units when feasible", {
  skip_without_gurobi()
  D <- matrix(4, 4, 4)
  diag(D) <- 0
  dat <- data.frame(
    Z = c(1, 1, 0, 0),
    X1 = c(0, 10, 0.1, 10.1),
    X2 = c(0, 0, 0, 0)
  )
  m <- netmatch(dat, "Z", c("X1", "X2"), D, method = "dual", kappa = 2,
                solver = "gurobi")
  expect_equal(m$solver, "gurobi")
  expect_false("method_label" %in% names(m))
  expect_false("match_table" %in% names(m))
  expect_equal(sort(as.integer(rownames(m$data))), 1:4)
})

test_that("dual penalty can use GLPK backend on a small feasible design", {
  testthat::skip_if_not_installed("Rglpk")
  testthat::skip_if_not_installed("slam")
  D <- matrix(4, 4, 4)
  diag(D) <- 0
  dat <- data.frame(
    Z = c(1, 1, 0, 0),
    X1 = c(0, 10, 0.1, 10.1),
    X2 = c(0, 1, 0.2, 1.2)
  )
  m <- netmatch(dat, "Z", c("X1", "X2"), D, method = "dual", kappa = 2, solver = "glpk")
  expect_equal(m$solver, "glpk")
  expect_equal(sort(as.integer(rownames(m$data))), 1:4)
})
