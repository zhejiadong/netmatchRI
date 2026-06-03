skip_without_gurobi <- function() {
  testthat::skip_if_not_installed("gurobi")
  ok <- tryCatch({
    model <- list(
      A = Matrix::Matrix(c(1, 1), nrow = 1, sparse = TRUE),
      obj = c(1, 2),
      sense = ">=",
      rhs = 1,
      vtype = c("B", "B"),
      modelsense = "min"
    )
    res <- gurobi::gurobi(model, params = list(OutputFlag = 0))
    identical(res$status, "OPTIMAL")
  }, error = function(e) FALSE)
  testthat::skip_if_not(ok, "Gurobi is installed but no working license/solver is available.")
}
