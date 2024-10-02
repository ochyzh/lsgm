# Gradient function for loglik_C:
loglik_gr <- function(par, X, W, Y) {
  xbeta <- X %*% par[seq_len(length(par) - 1)]
  kappa <- exp(xbeta) / (1 + exp(xbeta)) # logit of Xb
  W_Y_kappa <- W %*% (Y - kappa)
  A_i <- log(kappa / (1 - kappa)) + par[length(par)] * W_Y_kappa # Eqn 2
  p_i <- exp(A_i) / (1 + exp(A_i))
  # dl_d <- (Y / p_i - (1 - Y) / (1 - p_i)) / ((1 / p_i + 1 / (1 - p_i)))
  dl_d <- Y - p_i
  t(dl_d) %*% cbind(X, W_Y_kappa)
}

loglik <- function(par, X, W, Y) {
  xbeta <- X %*% par[seq_len(length(par) - 1)]
  kappa <- exp(xbeta) / (1 + exp(xbeta)) # logit of Xb
  A_i <- log(kappa / (1 - kappa)) + par[length(par)] * W %*% (Y - kappa) # Eqn 2
  p_i <- exp(A_i) / (1 + exp(A_i)) # Eqn 1, also Eqn 4
  PL <- Y * log(p_i) + (1 - Y) * log(1 - p_i) # Eqn 3
  -sum(PL)
}

spatbin.genone <- function(coeffs, X, W, curys) {
  eta<- coeffs[[length(coeffs)]]
  xbeta<-  X %*% coeffs[seq_len(length(coeffs) - 1)]
  kappa <- exp(xbeta) / (1 + exp(xbeta))
  A_i <- log(kappa / (1 - kappa)) + eta * W %*% (curys - kappa)
  p_i <- exp(A_i) / (1 + exp(A_i))
  rbinom(n = length(curys), size = 1, prob = p_i)
}

spatbin.onegibbs <- function(coeffs,X, W, curys) {
  for (cnt in seq_along(curys)) {
    ny <- spatbin.genone(coeffs = coeffs, X=X, W = W, curys = curys)
    curys[cnt] <- ny[cnt]
  }
  curys
}

spatbin.genfield <- function(coeffs,X, W, y0s, M) {
  res <- matrix(nrow = length(y0s), ncol = M + 1)
  res[, 1] <- y0s
  
  for (cnt in seq_len(M)) {
    res[, cnt + 1] <- spatbin.onegibbs(coeffs = coeffs,X=X, W = W, curys = res[, cnt])
  }

  as.data.frame(res)
}

sim_est <- function(Y, m1, W, X) {
  if (any(class(Y) %in% "matrix")) {
    Y <- as.matrix(Y)
  }

  tryCatch(
    res <- optim(par = m1$par, loglik, gr = function(...) loglik_gr(...) * 10e-3, W = W, Y = Y, X = X),
    error = function(cond) {
      message(cond)
      return(NULL)
    },
    warning = function(cond) {
      message(cond)
      return(NULL)
    }
  )

  c(res$par, res$convergence)
}
