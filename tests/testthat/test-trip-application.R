test_that("trip_application runs on a compatible local RData file", {
  ds_use <- data.frame(
    EGO_ID = 1:8,
    Z = c(1, 1, 1, 1, 0, 0, 0, 0),
    Y = c(3, 4, 5, 6, 2, 3, 4, 5),
    edu = c(1, 2, 3, 4, 1, 2, 3, 4),
    employment = c(1, 1, 2, 2, 1, 1, 2, 2),
    ACMDT = c(0, 0, 1, 1, 0, 0, 1, 1),
    baseRisk = c(1, 2, 1, 2, 1, 2, 1, 2),
    HIV = c(0, 0, 0, 1, 0, 0, 0, 1)
  )
  Adj <- matrix(0, 8, 8)
  Adj[cbind(1:7, 2:8)] <- 1
  Adj[cbind(2:8, 1:7)] <- 1
  f <- tempfile(fileext = ".RData")
  save(ds_use, Adj, file = f)

  out <- trip_application(
    data_file = f,
    output_dir = tempdir(),
    methods = "covariate",
    eta_grid = 0,
    rho_grid = 0
  )
  expect_true(nrow(out$tests) > 0)
  expect_true(file.exists(file.path(tempdir(), "trip_tests.csv")))
})
