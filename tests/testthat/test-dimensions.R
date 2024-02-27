library(testthat)

# Test cases for LSGM function in lsgm.R
test_that("lsgm function handles invalid inputs appropriately", {
  # Mock inputs
  X1 <- toy_data[1:10, ]
  W1 <- matrix(runif(100), nrow = 10, ncol = 10) # Not a square matrix
  Y1 <- rnorm(10) # Not a binary vector
  burnin <- 5
  thin <- 3

  # Expect error for non-square W
  expect_error(lsgm(Y = Y1, W = W1[, 1:9], X = X1, burnin = burnin, thin = thin), "non-conformable arguments")
})
