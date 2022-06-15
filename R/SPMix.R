#' @importFrom fmlogcondens fmlcd
#' @importFrom mclust dmvnorm
#' @importFrom LogConcDEAD mlelcd
#' @importFrom graphics legend
#' @import stats
#' 
#' @title Semiparametric Mixture Density Estimation for given z-values
#'
#' @description \code{SPMix} returns localFDR estimates and semiparametric
#' mixture density estimates for given multi-dimensional lists of z-values, which
#' are the probit-transformed p-values.
#' For the hypothesis testing \code{SPMix} uses a two-component semiparametric
#' mixture model to estimate the localFDR from the z-values. The two pillars of the
#' proposed approach are Efron's empirical null principle and log-concave density
#' estimation for the alternative distribution.
#'
#' @param z Matrix which column indicates z-values, probit-transformed p-values.
#' @param tol Stopping criteria for the EM algorithm. If maximum absolute difference
#' of current and previous gamma value is smaller than tol,
#' i.e. \eqn{max_i |\gamma_i^{(k+1)}-\gamma_i^{(k)} <tol}, for k-th step,
#' then optimization stops. (default: 5e-6)
#' @param p_value If TRUE, the column of input indicates p-values, if FALSE, it indicates z-values or raw data. (default: FALSE)
#' @param alternative a character string specifying the alternative hypothesis, must be one of "two.sided", "greater" (default) or "less". You can specify just the initial letter. (default: "greater")
#' @param max_iter Maximum number of iterations in the EM algorithm. (default: 30)
#' @param mono If TRUE, localFDR is in ascending order of z-values. (default: TRUE)
#' @param thre_z Threshold value which only z-values smaller than thre.z
#' are used to compute the log-concave estimates f_1 in M-step.
#' @param Uthre_gam Upper threshold of gamma which are used to compute stopping criteria for the EM algorithm.
#' @param Lthre_gam Lower threshold of gamma which are used to compute stopping criteria for the EM algorithm.
#'
#'
#' @return Estimates of semiparametric mixture model for f for given z-values.
#'
#'   \item{p0}{Prior probability for null distribution}
#'   \item{mu0 sig0}{Parameter estimates of normal null distribution, N(mu0, sig0^2)}
#'   \item{f}{Probability estimates of semiparametric mixture model for each z-value point.}
#'   \item{f1}{Probability estimates of alternative distribution of mixture model for each z-value point.}
#'   \item{localFDR}{localFDR estimates for given z-values}
#'   \item{iter}{Number of iterations of EM algorithm to compute localFDR.}
#'
#' @export

SPMix <- function(z, tol = 5e-6, p_value = FALSE, alternative = "greater", max_iter = 30, mono = TRUE, thre_z = 0.9,
                  Uthre_gam = 0.9, Lthre_gam = 0.01 )
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
  d <- dim(z)[2]

  if (p_value) {
    if (alternative == "greater" | alternative == "g") {
      z = qnorm(1-z)
    } else {
      z = qnorm(z)
    }
  } else {
    raw_mean = mean(z)
    if (d == 1) {raw_sd = sd(z)} else {raw_cov = cov(z)}
    z = scale(z)
  }

  ## Initial step: to fit normal mixture
  if (d == 1) {
    z = z[,1]
    if (alternative == "greater" | alternative == "g") {
      q0 <- quantile(z, probs = .9)
      p0 <- mean(z <= q0)
      mu0 <- mean(z[z <= q0])
      sig0 <- sd(z[z <= q0])
      f0 <- dmvnorm(z, mu0, sig0)
      mu1 <- mean(z[z > q0])
      sig1 <- sd(z[z > q0])
      f1 <- dnorm(z, mu1, sig1)
    } else {
      q0 <- quantile(z, probs = .7)
      p0 <- mean(z >= q0)
      mu0 <- mean(z[z >= q0])
      sig0 <- sd(z[z >= q0])
      f0 <- dmvnorm(z, mu0, sig0)
      mu1 <- mean(z[z < q0])
      sig1 <- sd(z[z < q0])
      f1 <- dnorm(z, mu1, sig1)
    }
  } else {
    Params <- NormMix(z)
    p0 <- Params$p0
    mu0 <- Params$mu0
    sig0 <- Params$sig0
    f0 <- dmvnorm(z, mu0, sig0)
    f1 <- dmvnorm(z, Params$mu1, Params$sig1)
  }
  gam <- f <- rep(0, n)


  if (d == 1) {
    z <- as.numeric(z)
    ## EM-step
    k <- 0; converged <- 0
    while ( (k < 3) | ((k < max_iter) & (!converged)) ) {
      k <- k + 1

      ## E-step
      new_f <- (p0 * f0 + (1 - p0) * f1)
      new_gam <- p0 * f0 / new_f

      ## M-step

      w_gam <- new_gam/sum(new_gam, na.rm = TRUE)
      new_mu0 <- sum(w_gam*z, na.rm = TRUE)
      new_sig0 <- sqrt(sum(w_gam*(z-new_mu0)^2, na.rm = TRUE))
      new_p0 <- mean(new_gam, na.rm = TRUE)
      new_f0 <- dnorm(z, new_mu0, new_sig0)
      new_f1 <- rep(0, n)
      which_z <- (new_gam <= thre_z)
      weight <- 1 - new_gam[which_z]
      weight <- weight/sum(weight)
      new_f1[which_z] <- exp(LogConcDEAD::mlelcd(z[which_z], w = weight)$logMLE)

      ## Update
      which_gam <- (new_gam <= Uthre_gam) * (new_gam >= Lthre_gam)
      diff <- max(abs(gam - new_gam)[which_gam])
      converged <- (diff <= tol)
      cat("   EM iteration:", k, ", Change in fdr fit = ", round(diff, 5), "\n")
      p0 <- new_p0
      mu0 <- new_mu0
      sig0 <- new_sig0
      f1 <- new_f1
      f0 <- new_f0
      f <- p0 * f0 + (1 - p0) * f1
      gam <- new_gam

    }
  } else {
    ## EM-step
    k <- 0; converged <- 0
    while ( (k < 3)|((k < max_iter) & (!converged)) ) {
      k <- k + 1

      ## E-step
      new_f <- (p0 * f0 + (1 - p0) * f1)
      new_gam <- p0 * f0 / new_f

      if (mono) new_gam <- MonotoneFDR(z, new_gam)

      ## M-step
      sum_gam <- sum(new_gam)
      new_mu0 <- as.vector(t(z) %*% new_gam) / sum_gam
      dev <- t(t(z)-new_mu0) * sqrt(new_gam)
      new_sig0 <- t(dev) %*% dev / sum_gam
      new_p0 <- mean(new_gam)
      new_f0 <- dmvnorm(z, new_mu0, new_sig0)
      weight <- 1 - new_gam
      new_f1 <- rep(0, n)
      which_z <- (new_gam <= thre_z)
      lcd <- fmlogcondens::fmlcd(X=z[which_z,], w = weight[which_z] / sum(weight[which_z]))
      new_f1[which_z] <- exp(lcd$logMLE)

      ## Update
      which_gam <- (new_gam <= Uthre_gam) * (new_gam >= Lthre_gam)
      diff <- max(abs(gam - new_gam)[which_gam])
      converged <- (diff <= tol)
      cat("   EM iteration:", k, ", Change in mdfdr fit = ", round(diff, 5), "\n")
      p0 <- new_p0; mu0 <- new_mu0; sig0 <- new_sig0
      f1 <- new_f1
      f0 <- new_f0
      f <- p0 * f0 + (1 - p0) * f1
      gam <- new_gam
    }
  }

  if (p_value) {
    res <- list(p0 = p0, mu0 = mu0, sig0 = sig0,
                f = f, f1 = f1, localFDR = gam, iter = k)
  } else {
    if (d == 1){
      res <- list(p0 = p0, mu0 = mu0*raw_sd + raw_mean, sig0 = sig0*raw_sd,
                  f = f/raw_sd, f1 = f1/raw_sd, localFDR = gam, iter = k)
    } else {
      res <- list(p0 = p0, mu0 = mu0%*%raw_cov + raw_mean, sig0 = sig0%*%raw_cov,
                  f = f, f1 = f1, localFDR = gam, iter = k)
    }
  }
  return(res)
}
