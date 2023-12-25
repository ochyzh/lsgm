#' Estimation of a Local Structure Graph Model (LSGM)
#'
#' This function estimates an LSGM, developed in Chyzh and Kaiser (2019)
#'
#' @param X a NxK data.frame that contains K independent variables for N observations. Must be of class data.frame.
#' @param W an NXN square matrix whose values represent connectivity between pairs of edges
#' @param Y a vector of length N that contains the values of the binary dependent variable.
#' @param burnin a scalar specifying burnin.
#' @param thin a scalar specifying thinning.
#' @param MCMC a scalar specifying the number of iterations for the Gibbs sampler.
#' @param seed random seed
#'
#' @return the output is a table of estimated coefficients, standard errors and z-statistics, *eta* is the estimated parameter on the spatial term
#' @keywords LSGM
#' @references{
#'  Chyzh, Olga V. and Mark S. Kaiser. 2019. "A Local Structure Graph Model: Modeling Formation
#' of Network Edges as a Function of Other Edges." \emph{Political Analysis}. 27 (4): 397--414. \url{https://doi.org/10.1017/pan.2019.8}
#'
#' Nieman, Mark David, Carla Martinez Machain, Olga V. Chyzh, and Sam Bell. 2021.
#' "An International Game of Risk: Troop Placement and Major Power Competition." \emph{The Journal of Politics}. \url{https://doi.org/10.1086/711716}
#' }
#'
#' @examples
#' data(W)
#' data(toy_data)
#' lsgm(Y = as.matrix(toy_data$Y), W = W, X = as.data.frame(toy_data$X))
#'
#' @importFrom stats optim rbinom sd
#' @export
lsgm <- function(Y, W, X = NULL, burnin = 10000, thin = 100, MCMC = 10000, seed = 4448) {
  set.seed(seed)
  MCMC <- MCMC + burnin
  n <- length(Y)


  if (!is.null(X)) {
    K <- ncol(X)
    if (!is.null(names(X))) {
      mynames <- names(X)
    } else {
      mynames <- paste("Var", seq_len(ncol(X)), sep = "")
    }


    # m0<-glm(formula = Y ~ as.matrix(X), family = binomial(link = "logit"))
    X <- cbind(1, X) # Add a column for constant
    X <- as.matrix(X)
  } else {
    K <- 0
    #  m0<-glm(formula = Y ~ 1, family = binomial(link = "logit"))
    X <- matrix(1, nrow = n, ncol = 1)
    mynames <- NULL
  }

  if (!any(class(Y) %in% "matrix")) {
    Y <- as.matrix(Y)
  }
  if (!any(class(W) %in% "matrix")) {
    W <- as.matrix(W)
  }

  pars <- rep(0, (K + 2))
  # pars<-c(m0$coef,0)

  m1 <- optim(par = pars, loglik, gr = function(...) loglik_gr(...) * 10e-3, W = W, Y = Y, X = X)

  y0s <- rbinom(n = n, size = 1, prob = .5)
  sims <- spatbin.genfield(coeffs = m1$par, W = W, y0s = y0s, M = MCMC)

  # Take every 10th simulated network, i.e. burnin=10, thinning=10
  sims <- sims[, seq(from = burnin + 1, to = ncol(sims), by = thin)]

  sim_est <- do.call("rbind", lapply(sims, function(sim) sim_est(sim, m1, W, X)))

  # Drop results if didn't converge (models that converged have convergence=0)
  sim_est <- sim_est[sim_est[, ncol(sim_est)] == 0, ]

  # Get sds of the estimates:
  boot_se <- apply(sim_est, 2, sd)
  mytable <- cbind("coeff" = m1$par, "se" = boot_se[-ncol(sim_est)], "z-value" = (m1$par / boot_se[-ncol(sim_est)]))
  row.names(mytable) <- c("constant", mynames, "eta")

  return(mytable)
}
