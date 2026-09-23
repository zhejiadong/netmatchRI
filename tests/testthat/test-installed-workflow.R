test_that("the installed basic workflow runs from source", {
  workflow <- system.file("examples", "basic-workflow.R", package = "netmatchRI")
  expect_true(nzchar(workflow))
  expect_true(file.exists(workflow))

  env <- new.env(parent = globalenv())
  expect_silent(output <- capture.output(sys.source(workflow, envir = env)))
  expect_true(any(grepl("Critical eta (rho = 1)", output, fixed = TRUE)))
  expect_true(any(grepl("Critical ratio", output, fixed = TRUE)))

  expect_s3_class(env$m, "netmatch")
  expect_s3_class(env$match_summary, "netmatch_summary")
  expect_named(env$diagnostics,
               c("covariate_balance", "network_summary", "within_distance_table"))
  expect_s3_class(env$ri_unadjusted, "netmatch_test")
  expect_s3_class(env$ri_adjusted, "netmatch_test")
  expect_s3_class(env$ri_design, "netmatch_test")
  expect_s3_class(env$sens, "netmatch_sensitivity")
  expect_s3_class(env$critical, "netmatch_critical_sensitivity")
})
