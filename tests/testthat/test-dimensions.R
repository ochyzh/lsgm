# Test cases for LSGM function in lsgm.R
test_that("lsgm fails when dimensions are incompatible", {
  # Mock inputs
  X1 <- as.matrix(toy_data$X[1:10, ])
  W1 <- matrix(runif(20), nrow = 10, ncol = 2)
  Y1 <- as.matrix(rbinom(10, size=1, prob=.5))
  burnin <- 5
  thin <- 3

  # Expect error for non-square W
  expect_error(lsgm(Y = Y1, W = W1[, 1:2], X = X1, burnin = burnin, thin = thin), "non-conformable arguments")
})
