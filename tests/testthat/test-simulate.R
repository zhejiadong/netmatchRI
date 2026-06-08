test_that("simulate_netmatch_example exposes dependence and network parameters", {
  sim <- simulate_netmatch_example(
    seed = 123,
    n = 40,
    alpha1 = 0.5,
    alpha2 = 0.2,
    pin = 0.3,
    pout = 0.01
  )

  expect_s3_class(sim$data, "data.frame")
  expect_equal(nrow(sim$data), 40)
  expect_equal(dim(sim$Adj), c(40L, 40L))
  expect_equal(dim(sim$V), c(40L, 40L))
  expect_equal(sim$seed, 123)
  expect_equal(sim$alpha1, 0.5)
  expect_equal(sim$alpha2, 0.2)
  expect_equal(sim$pin, 0.3)
  expect_equal(sim$pout, 0.01)
  expect_false("dep_index" %in% names(sim))
})

test_that("simulate_netmatch_example validates unit-interval parameters", {
  expect_error(simulate_netmatch_example(n = 40, alpha1 = 1.1), "alpha1")
  expect_error(simulate_netmatch_example(n = 40, alpha2 = -0.1), "alpha2")
  expect_error(simulate_netmatch_example(n = 40, pin = 2), "pin")
  expect_error(simulate_netmatch_example(n = 40, pout = NA_real_), "pout")
})
