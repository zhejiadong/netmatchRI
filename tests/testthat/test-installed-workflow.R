test_that("the installed basic workflow runs from source", {
  workflow <- system.file("examples", "basic-workflow.R", package = "netmatchRI")
  expect_true(nzchar(workflow))
  expect_true(file.exists(workflow))

  env <- new.env(parent = globalenv())
  expect_silent(sys.source(workflow, envir = env))

  expect_s3_class(env$match, "netmatch")
  expect_s3_class(env$match_summary, "netmatch_summary")
  expect_named(env$diagnostics,
               c("covariate_balance", "network_summary", "within_distance_table"))
  expect_s3_class(env$ri_naive, "netmatch_test")
  expect_s3_class(env$ri_sensitivity, "netmatch_test")
  expect_s3_class(env$ri_design, "netmatch_test")
  expect_s3_class(env$grid, "netmatch_sensitivity")
  expect_s3_class(env$critical, "netmatch_critical_sensitivity")
})
