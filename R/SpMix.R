#' @importFrom fmlogcondens fmlcd
#' @importFrom mclust dmvnorm
#' @importFrom logcondens activeSetLogCon
#'
#' @title LocalFDR estimation for given z-values
#'
#' @description \code{SpMix} returns LocalFDR estimates and semiparametric
#' mixture density estimates for given multi-dimensional lists of z-values, which
#' are the probit-transformed p-values.
#' For the hypothesis testing \code{SpMix} uses a two-component semiparametric
#' mixture model to estimate the LocalFDR from the p-values. The two pillars of the
#' proposed approach are Efron's empirical null principle and log-concave density
#' estimation for the alternative distribution.
#'
#' @param z Matrix which column indicates z-values, probit-transformed p-values.
#' @param tol Stopping criteria for the EM algorithm. If maximum absolute difference
#' of current and previous gamma value is smaller than tol,
#' i.e. \eqn{max_i |\gamma_i^{(k+1)}-\gamma_i^{(k)} <tol}, for k-th step,
#' then optimization stops. (default: 5e-6)
#' @param max.iter Maximum number of iterations in the EM algorithm. (default: 30)
#' @param mono If TRUE, LocalFDR is in ascending order of z-values. (default: TRUE)
#' @param thre.z Threshold value which only z-values smaller than thre.z
#' are used to compute the log-concave estimates f_1 in M-step.
#' @param Uthre.gam Upper threshold of gamma which are used to compute stopping criteria for the EM algorithm.
#' @param Lthre.gam Lower threshold of gamma which are used to compute stopping criteria for the EM algorithm.
#'
#'
#' @return Estimates of semiparametric mixture model for f for given z-values.
#'
#'   \item{p.0}{Prior probability for null distribution}
#'   \item{mu.0 sig.0}{Parameter estimates of normal null distribution, N(mu.0, sig.0^2)}
#'   \item{f}{Probability estimates of semiparametric mixture model for each z-value point.}
#'   \item{localfdr}{LocalFDR estimates for given z-values}
#'   \item{iter}{Number of iterations of EM algorithm to compute LocalFDR.}
#'
#' @export
SpMix <- function(z, tol = 5e-6, max.iter = 30, mono = TRUE, thre.z = 0.9,
                  Uthre.gam = 0.9, Lthre.gam = 0.01, )
{
  # *****************DEFINITION OF INTERNAL FUNCTIONS ******************

  NormMix <- function(z, tol = 5e-3, max_iter = 10)
  {

    k <- 0; converged <- 0
    z <- as.matrix(z)
    m_dist <- mahalanobis(z, 0, cov(z))

    p0 <- mean(m_dist <= 1.65)
    mu0 <- rep(0, dim(z)[2])
    sig0 <- diag(1, dim(z)[2])
    f0 <- dmvnorm(z, mu0, sig0)
    mu1 <- apply(as.matrix(z[m_dist > 1.65,]), 2, mean)
    sig1 <- cov(as.matrix(z[m_dist > 1.65,]))
    f1 <- dmvnorm(z, mu1, sig1)

    while ((k < 3) | ((k < max_iter) & (!converged))) {
      k <- k + 1

      ## E-step
      gam <- p0 * f0 / (p0 * f0 + (1-p0) * f1)

      ## M-step
      new_p0 <- mean(gam)
      new_mu0 <- as.vector(t(z) %*% gam) / sum(gam)
      dev0 <- (z - new_mu0) * sqrt(gam)
      new_sig0 <- t(dev0) %*% dev0 / sum(gam)
      f0 <- dmvnorm(z, new_mu0, new_sig0)
      new_mu1 <- as.vector(t(z) %*% (1 - gam)) / sum(1 - gam)
      dev1 <- (z - new_mu1) * sqrt(1 - gam)
      new_sig1 <- t(dev1) %*% dev1 / sum(1 - gam)

      ## Update
      diff <- max(abs((new_mu0 - mu0)),
                  abs((new_sig0 - sig0)),
                  abs((new_mu1 - mu1)),
                  abs((new_sig1 - sig1)),
                  abs((new_p0 - p0)))
      converged <- (diff <= tol)
      p0 <- new_p0
      mu0 <- new_mu0
      sig0 <- new_sig0
      mu1 <- new_mu1
      sig1 <- new_sig1
    }

    return(list(p0 = p0, mu0 = mu0, sig0 = sig0, mu1 = mu1, sig1 = sig1))
  }



  NE <- function(x, X)
  {
    n <- nrow(X)
    p <- ncol(X)
    xx <- matrix(x, nrow = n, ncol = p, byrow = TRUE)
    ne.ind <- apply(1*(X >= xx), 1, prod)

    return((1:n)[ne.ind == 1])
  }


  MonotoneFDR <- function(z, fdr)
  {
    n <- nrow(z)
    MFDR <- numeric(n)
    for (i in 1:n) {
      MFDR[i] <- max(fdr[NE(z[i,], z)])
    }

    return(MFDR)
  }

  # ******************* MAIN FUNCTION *******************************

  z <- as.matrix(z)
  n <- dim(z)[1]

  ## Initial step: to fit normal mixture
  if (dim(z)[2] == 1) {
    q0 <- quantile(z, probs = .9)
    p0 <- mean(z <= q0)
    mu0 <- mean(z[z <= q0])
    sig0 <- sd(z[z <= q0])
    f0 <- dmvnorm(z, mu0, sig0)
    mu1 <- mean(z[z > q0])
    sig1 <- sd(z[z > q0])
    f1 <- dnorm(z, mu1, sig1)
  }
  else {
    Params <- NormMix(z)
    p0 <- Params$p0
    mu0 <- Params$mu0
    sig0 <- Params$sig0
    f1 <- dmvnorm(z, Params$mu1, Params$sig1)
  }
  gam <- f <- rep(0, n)

  ## EM-step
  k <- 0; converged <- 0
  while ( (k < 3)|((k < max.iter) & (!converged)) ) {
    k <- k + 1

    ## E-step
    new_gam <- p0 * f0 / (p0 * f0 + (1 - p0) * f1)

    if (mono) new_gam <- MonotoneFDR(z, new_gam)

    ## M-step
    sum.gam <- sum(new.gam)
    new.mu.0 <- as.vector(t(z)%*%new.gam)/sum.gam
    dev <- t(t(z)-new.mu.0)*sqrt(new.gam)
    new.sig.0 <- t(dev)%*%dev/sum.gam
    new.p.0 <- mean(new.gam)
    new.f.0 <- dmvnorm(z, new.mu.0, new.sig.0)
    weight <- 1 - new.gam
    new.f1.tilde <- rep(0, n)
    which.z <- (new.gam <= thre.z)
    lcd <- fmlogcondens::fmlcd(X=z[which.z,], w = weight[which.z]/sum(weight[which.z]))
    new.f1.tilde[which.z] <- exp(lcd$logMLE)

    ## Update
    which.gam <- (new.gam <= Uthre.gam)*(new.gam >= Lthre.gam)
    diff <- max(abs(gam - new.gam)[which.gam])
    converged <- (diff <= tol)
    cat("   EM iteration:", k, ", Change in mdfdr fit = ", round(diff, 5), "\n")
    p.0 <- new.p.0; mu.0 <- new.mu.0; sig.0 <- new.sig.0
    f1.tilde <- new.f1.tilde
    gam <- new.gam
    f <- new.f
  }

  res <- list(p.0 = p.0, mu.0 = mu.0, sig.0 = sig.0,
              f1.hat = f1.tilde, f = f, localfdr = gam, iter = k)

  return(res)
}

