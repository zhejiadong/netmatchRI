test_that("paper reproduction smoke test runs", {
  out <- reproduce_paper("table1", fast = TRUE, output_dir = tempdir())
  expect_true("table1" %in% names(out))
  expect_true(nrow(out$table1) > 0)
  expect_true(file.exists(file.path(tempdir(), "table1_raw.csv")))
})
