test_that("lsgm parses data to matrices", {
  # Mock inputs
  X1 <- toy_data$X[1:10, ]
  W1 <- as.data.frame(matrix(runif(100), nrow = 10, ncol = 10))
  Y1 <- rnorm(10)
  burnin <- 5
  thin <- 3

  # Expect error for non-square W
  expect_type(lsgm(Y = Y1, W = W1, X = X1, burnin = burnin, thin = thin), "double")
})

test_that("lsgm handles variable names", {
  # Mock inputs
  X1 <- as.matrix(toy_data$X[1:10, ])
  colnames(X1) <- "A"

  W1 <- matrix(runif(100), nrow = 10, ncol = 10)

  Y1 <- as.matrix(rnorm(10))

  burnin <- 5
  thin <- 3

  expect_type(lsgm(Y = Y1, W = W1, X = X1, burnin = burnin, thin = thin), "double")
})

test_that("lsgm works with a null X", {
  # Mock inputs
  W1 <- matrix(runif(25), nrow = 5, ncol = 5)

  Y1 <- as.matrix(rnorm(5))

  burnin <- 2
  thin <- 1

  expect_type(lsgm(Y = Y1, W = W1, X = NULL, burnin = burnin, thin = thin), "double")
})
