# Gradient function for loglik_C:
loglik_gr <- function(par, X, W, Y) {
  xbeta <- X %*% par[seq_len(length(par) - 1)]
  kappa <- exp(xbeta) / (1 + exp(xbeta)) # logit of Xb
  A_i <- log(kappa / (1 - kappa)) + par[length(par)] * W %*% (Y - kappa) # Eqn 2
  p_i <- exp(A_i) / (1 + exp(A_i))
  dl_d <- (Y / p_i - (1 - Y) / (1 - p_i)) / ((1 / p_i + 1 / (1 - p_i)))
  t(dl_d) %*% cbind(X, W %*% (Y - kappa))
}

loglik <- function(par, X, W, Y) {
  xbeta <- X %*% par[seq_len(length(par) - 1)]
  kappa <- exp(xbeta) / (1 + exp(xbeta)) # logit of Xb
  A_i <- log(kappa / (1 - kappa)) + par[length(par)] * W %*% (Y - kappa) # Eqn 2
  p_i <- exp(A_i) / (1 + exp(A_i)) # Eqn 1, also Eqn 4
  PL <- Y * log(p_i) + (1 - Y) * log(1 - p_i) # Eqn 3
  -sum(PL)
}

spatbin.genone <- function(coeffs, W, curys) {
  b0 <- coeffs[1]
  eta <- coeffs[2]
  xbeta <- b0
  kappa <- exp(xbeta) / (1 + exp(xbeta))
  A_i <- log(kappa / (1 - kappa)) + eta * W %*% (curys - kappa)
  p_i <- exp(A_i) / (1 + exp(A_i))
  y <- rbinom(n = length(curys), size = 1, prob = p_i)
  return(y)
}

spatbin.onegibbs <- function(coeffs, W, curys) {
  cnt <- 0
  n <- length(curys)
  newys <- NULL
  repeat{
    cnt <- cnt + 1
    ny <- spatbin.genone(coeffs = coeffs, W = W, curys = curys)
    curys[cnt] <- ny[cnt]
    if (cnt == n) break
  }
  newys <- curys
  return(newys)
}

spatbin.genfield <- function(coeffs, W, y0s, M) {
  curys <- y0s
  cnt <- 0
  res <- as.data.frame(y0s)
  repeat{
    cnt <- cnt + 1
    newys <- spatbin.onegibbs(coeffs = coeffs, W = W, curys = curys)
    curys <- newys
    res <- cbind(res, curys)
    if (cnt == M) break
  }

  return(res)
}

sim_est <- function(Y, m1, W, X) {
  if (any(class(Y) %in% "matrix")) {
    Y <- as.matrix(Y)
  }

  tryCatch(
    res <- optim(par = m1$par, loglik, gr = function(...) loglik_gr(...) * 10e-3, W = W, Y = Y, X = X),
    error = function(cond) {
      message(cond)
      # Choose a return value in case of error
      return(NULL)
    },
    warning = function(cond) {
      message(cond)
      # Choose a return value in case of warning
      return(NULL)
    }
  )
  return(c(res$par, res$convergence))
}
