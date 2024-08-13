#' @importFrom mclust dmvnorm
#' @importFrom mvtnorm pmvnorm
#' @importFrom LogConcDEAD mlelcd
#' @importFrom logcondens activeSetLogCon
#' @importFrom logcondens evaluateLogConDens
#' @importFrom graphics legend
#' @import stats
#' 
#' @title Parameter estimates of null(normal) distribution and fitted values for 
#' both alternative(nonparametric) and mixture density. 
#' 
#' @description \code{SPMix} returns localFDR estimates and semiparametric
#' mixture density estimates for given multi-dimensional lists of z-values, p-values or raw data. 
#' For the hypothesis testing \code{SPMix} uses a two-component semiparametric
#' mixture model to estimate the localFDR from the z-values. The two pillars of the
#' proposed approach are Efron's empirical null principle and log-concave density
#' estimation for the alternative distribution.
#'
#' @param z Matrix which each row indicates each data point (z-values, p-values, or raw data).
#' @param tol Stopping criteria for the EM algorithm. If maximum absolute difference
#' of current and previous gamma value is smaller than tol,
#' i.e. \eqn{max_i |\gamma_i^{(k+1)}-\gamma_i^{(k)} <tol}, for k-th step,
#' then optimization stops. (default: 5e-6)
#' @param p_value If TRUE, input data indicates p-values, if FALSE, it indicates z-values or raw data. (default: FALSE)
#' @param greater_alt A vector indicating if the alternative distribution is greater (`TRUE`) or less (`FALSE`) than the null distribution for each dimension. 
#' For multidimensional data, the vector should match the number of columns in `z`.
#' @param min_iter Minimum number of iterations in the EM algorithm. (default: 3)
#' @param max_iter Maximum number of iterations in the EM algorithm. (default: 30)
#' @param thre_z The upper threshold of gamma whose z-values are used in log-concave estimates in the M-step of the EM-type algorithm. (default: 1-1e-5)
#' @param Uthre_gam The upper threshold of gamma which are used to compute stopping criteria for the EM algorithm. (default: 0.99)
#' @param Lthre_gam The lower threshold of gamma which are used to compute stopping criteria for the EM algorithm. (default: 0.01)
#'
#' @return Estimates of semiparametric mixture model for given data.
#'
#'   \item{z}{Matrix which each row indicates each data point}
#'   \item{p0}{Prior probability for null distribution}
#'   \item{mu0 sig0}{Parameter estimates of Gaussian (null) distribution, N(mu0, sig0^2)}
#'   \item{f}{Probability estimates of semiparametric mixture model for given data.}
#'   \item{f1}{Probability estimates of log-concave (alternative) distribution of mixture model for given data.}
#'   \item{F}{Cumulative density estimates of mixture model for given data.}
#'   \item{localFDR}{localFDR estimates for given data.}
#'   \item{FDR}{FDR estimates for given data.}
#'   \item{iter}{Number of iterations of EM algorithm to compute localFDR.}
#'   \item{dim}{Dimension of the given data}
#'   \item{greater_alt}{A vector indicating if the alternative distribution is greater (`TRUE`) or less (`FALSE`) than the null distribution for each dimension.}
#'
#' @export

SPMix <- function(z, init_p0 = 0.5, tol = 5.0e-5,
                  greater_alt = NULL, min_iter = 3, max_iter = 50, 
                  thre_z = 0.5, Uthre_gam = 0.99, Lthre_gam = 0.01 )
{
  # *****************DEFINITION OF INTERNAL FUNCTIONS ******************
  
  NormMix <- function(z, p0, tol = 5e-3, max_iter = 5)
  {
    require(mvtnorm)
    
    k <- 0; converged <- 0
    z <- as.matrix(z)
    n <- nrow(z)
    d <- ncol(z)
    q <- qnorm(p0)
    z.scaled <- scale(z)
    ind0 <- ind1 <- rep(1, n)
    for (j in 1:d) {
      ind0 <- ind0*(z.scaled[,j] <= q)
      ind1 <- ind1*(z.scaled[,j] > -q)
    }
    
    z0 <- as.matrix(z[ind0 == 1,])
    z1 <- as.matrix(z[ind1 == 1,])
    
    mu0 <- colMeans(z0)
    sig0 <- cov(z0)
    f0 <- dmvnorm(z, mu0, sig0)
    
    mu1 <- colMeans(z1)
    sig1 <- cov(z1)
    f1 <- dmvnorm(z, mu1, sig1)
    
    f <- p0*f0 + (1 - p0)*f1
    ell <- mean(log(f))
    
    while ((k < 3) | ((k < max_iter) & (!converged))) {
      k <- k + 1
      
      ## E-step
      gam <- p0 * f0 / (p0 * f0 + (1 - p0) * f1)
      
      ## M-step
      new_p0 <- mean(gam)
      new_mu0 <- as.vector(t(z) %*% gam) / sum(gam)
      dev0 <- (z - new_mu0) * sqrt(gam)
      new_sig0 <- t(dev0) %*% dev0 / sum(gam)
      f0 <- dmvnorm(z, new_mu0, new_sig0)
      
      new_mu1 <- as.vector(t(z) %*% (1 - gam)) / sum(1 - gam)
      dev1 <- (z - new_mu1) * sqrt(1 - gam)
      new_sig1 <- t(dev1) %*% dev1 / sum(1 - gam)
      f1 <- dmvnorm(z, new_mu1, new_sig1)
      
      f <- p0*f0 + (1 - p0)*f1
      new_ell <- mean(log(f))
      
      ## Update
      diff <- abs(new_ell - ell)
      converged <- (diff <= tol)
      p0 <- new_p0
      mu0 <- new_mu0
      sig0 <- new_sig0
      mu1 <- new_mu1
      sig1 <- new_sig1
      ell <- new_ell
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
  
  mult.ecdf <- function(x)
  {
    x <- as.matrix(x)
    n <- nrow(x)
    d <- ncol(x)
    Fn <- rep(0, n)
    for ( i in 1:n ) {
      tmp <- matrix(x[i,], n, d, byrow = TRUE)
      Fn[i] <- mean(rowSums(x <= tmp) == d)
    }
    return(Fn)
  }
  
  # Normal Mixture Model for 1D data
  normal.mixture.1d <- function(z, p0, tol = 5e-3, max.iter = 5) {
    q0 <- quantile(z, probs = p0)
    mu0 <- mean(z[z <= q0])
    sig0 <- sd(z[z <= q0])
    f0 <- dnorm(z, mean = mu0, sd = sig0)
    
    mu1 <- mean(z[z > q0])
    sig1 <- sd(z[z > q0])
    f1 <- dnorm(z, mean = mu1, sd = sig1)
    ell <- sum(log(p0 * f0 + (1 - p0) * f1))
    
    k <- 0
    while (k < max.iter) {
      k <- k + 1
      
      # E-step
      term1 <- p0 * f0
      term2 <- term1 + (1 - p0) * f1
      gam <- term1 / term2
      
      # M-step
      p0 <- mean(gam)
      w0 <- gam / sum(gam)
      mu0 <- sum(z * w0)
      sig0 <- sqrt(sum(w0 * (z - mu0)^2))
      f0 <- dnorm(z, mean = mu0, sd = sig0)
      
      w1 <- (1 - gam) / sum(1 - gam)
      mu1 <- sum(z * w1)
      sig1 <- sqrt(sum(w1 * (z - mu1)^2))
      f1 <- dnorm(z, mean = mu1, sd = sig1)
      
      # Update log-likelihood
      new.ell <- sum(log(p0 * f0 + (1 - p0) * f1))
      diff <- abs(ell - new.ell)
      ell <- new.ell
      if (diff < tol) break
    }
    
    list(p0 = p0, mu0 = mu0, sigma0 = sig0, mu1 = mu1, sigma1 = sig1)
  }

  # ******************* MAIN FUNCTION *******************************
  
  z <- as.matrix(z)
  n <- nrow(z)
  d <- ncol(z)
  ell <- rep(NA, max_iter)

  if (p_value) {
    if (is.null(greater_alt)) {
      z = qnorm(1-z)
    } else {
      z = qnorm(z)
    }
  } else {
    raw_mean = mean(z)
    if (d == 1) {
      raw_sd = sd(z); z = scale(z)
    }
  }
  
  if(!is.null(greater_alt)) {
    for(j in 1:d) {
      if(greater_alt[j] == FALSE) z[,j] <- -z[,j]
    }
  }
  
  ## Initial step: to fit normal mixture
  if (d == 1) {
    z = z[,1]
    # Initial step: Fit normal mixture
    nmEM <- normal.mixture.1d(z, p0 = init_p0)
    p0 <- nmEM$p0
    f0 <- p0 * dnorm(z, mean = nmEM$mu0, sd = nmEM$sigma0)
    f1 <- dnorm(z, mean = nmEM$mu1, sd = nmEM$sigma1)
    f <- p0 * f0 + (1 - p0) * f1
  } else {
    ## Initial step: to fit normal mixture
    Params <- NormMix(z, p0 = init_p0)
    p0 <- Params$p0
    mu0 <- Params$mu0
    sig0 <- Params$sig0
    f0 <- dmvnorm(z, mu0, sig0)
    f1 <- dmvnorm(z, Params$mu1, Params$sig1)
    f <- p0*f0 + (1 - p0)*f1
    ell[1] <- mean(log(f), na.rm = TRUE)
  }
  if (d == 1) {
    z <- as.numeric(z)
    ## EM-step
    k <- 0; converged <- 0
    while (k < max.iter) {
      k <- k + 1
      
      # E-step
      gam <- p0 * f0 / f
      
      # M-step
      weight <- gam / sum(gam, na.rm = TRUE)
      mu0 <- sum(weight * z, na.rm = TRUE)
      sig0 <- sqrt(sum(weight * (z - mu0)^2, na.rm = TRUE))
      f0 <- p0 * dnorm(z, mean = mu0, sd = sig0)
      p0 <- mean(gam, na.rm = TRUE)
      
      # Update gam
      gam <- p0 * f0 / f
      which.z <- gam <= thre
      weight <- (1 - gam[which.z]) / sum(1 - gam[which.z], na.rm = TRUE)
      z1 <- z[which.z] + rnorm(sum(which.z), sd = 1e-5 * sd(z[which.z]))
      lcd <- activeSetLogCon(x = z1, w = weight)
      f1 <- numeric(n)
      f1[which.z] <- exp(lcd$phi)[rank(z1)]
      
      # Update f and log-likelihood
      f <- p0 * f0 + (1 - p0) * f1
      ell[k] <- mean(log(p0 * f0[!which.z]) + log((1 - p0) * f1[which.z]))
      cat(".")
      
      # Check for convergence
      if (k >= 7) {
        diff <- abs(ell[k] - ell[k - 1])
        converged <- (diff < tol)
      }
      if (converged) {
        cat("Converged!\n")
        break
      }
    }
    
    if (!converged) cat("Warning: Not converged!\n")
    # Final posterior probabilities
    gam <- p0 * f0 / f
    ell <- na.omit(ell)
    
  } else {
    ## EM-step
    k <- 1; converged <- 0
    while ( (k < 5) | ((k < max_iter) & (!converged)) ) {
      k <- k + 1
      
      ## E-step
      gam <- p0 * f0 / f
      gam <- MonotoneFDR(z, gam)
      
      ## M-step
      sum_gam <- sum(gam)
      mu0 <- as.vector(t(z) %*% gam) / sum_gam
      dev <- t(t(z) - mu0) * sqrt(gam)
      sig0 <- t(dev) %*% dev / sum_gam
      p0 <- mean(gam)
      f0 <- dmvnorm(z, mu0, sig0)
      f1 <- rep(0, n)
      which_z <- (gam <= thre)
      weight <- 1 - gam[which_z]
      weight <- weight/sum(weight)
      #lcd <- fmlogcondens::fmlcd(X = z[which_z,], w = weight[which_z] / sum(weight[which_z]))
      lcd <- mlelcd(z[which_z,], w = weight)
      f1[which_z] <- exp(lcd$logMLE)
      f <- p0*f0 + (1 - p0)*f1
      ell[k] <- mean(log(f), na.rm = TRUE)
      
      ## Update
      diff <- abs(ell[k] - ell[k - 1])
      converged <- (diff <= tol)
      cat(".")
      if (converged) {
        cat("Converged!\n")
        break
      }
    }
    if(!converged) cat("Warning: Not converged!\n")
  }

  # return results

  if (p_value) {
    res <- list(z = z, p0 = p0, mu0 = mu0, sig0 = sig0, f = f, f1 = f1,  
                log.likelihood = ell,
                localFDR = gam, posterior = cbind(gam, 1 - gam),
                iter = k, dim = d, greater_alt = greater_alt,
                converged = converged)
  } else {
    if (d == 1){
      res <- list(z = z*raw_sd + raw_mean, p0 = p0, mu0 = mu0*raw_sd + raw_mean, 
                  sig0 = sig0*raw_sd, f = f/raw_sd, f1 = f1/raw_sd, 
                  log.likelihood = ell,
                  localFDR = gam, posterior = cbind(gam, 1 - gam),
                  iter = k, dim = d, greater_alt = greater_alt,
                  converged = converged)
    } else {
      # return results
      res <- list(z = z, p0 = p0, mu0 = mu0, sig0 = sig0, f1 = f1, f = f,
                  iter = k, log.likelihood = ell, lcd = lcd, dim = d,
                  posterior = cbind(gam, 1 - gam),
                  converged = converged)
    }
  }
  
  return(res)
}
