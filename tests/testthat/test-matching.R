test_that("single penalty uses direct network distance threshold", {
  D <- matrix(3, 6, 6)
  diag(D) <- 0
  D[1, 4] <- D[4, 1] <- 1
  dat <- data.frame(
    Z = c(1, 1, 1, 0, 0, 0),
    X1 = c(0, 10, 20, 0.1, 10.1, 20.1),
    X2 = c(0, 0, 0, 0, 0, 0)
  )
  m <- netmatch(dat, "Z", c("X1", "X2"), D, method = "single", kappa = 2)
  expect_false(any(rownames(m$data)[m$data$subclass == m$data$subclass[rownames(m$data) == "1"]] == "4"))
})

test_that("dual penalty avoids close same-arm controls in one set", {
  D <- matrix(4, 5, 5)
  diag(D) <- 0
  D[2, 3] <- D[3, 2] <- 1
  dat <- data.frame(
    Z = c(1, 0, 0, 0, 0),
    X1 = c(0, 0.1, 0.2, 5, 6),
    X2 = c(0, 0, 0, 0, 0)
  )
  m <- netmatch(dat, "Z", c("X1", "X2"), D, method = "dual", kappa = 2, max_controls = 3)
  controls_in_set <- as.integer(rownames(m$data[m$data$Z == 0, ]))
  expect_false(all(c(2, 3) %in% controls_in_set))
})
