#' @importFrom mclust dmvnorm
#' @importFrom mvtnorm pmvnorm
#' @importFrom LogConcDEAD mlelcd
#' @importFrom logcondens activeSetLogCon
#' @importFrom logcondens evaluateLogConDens
#' @importFrom graphics legend
#' @import stats
#' 
#' @title Parameter Estimates of Null Distribution and Fitted Values for Mixture and Nonparametric Density
#' 
#' @description \code{SpMix} computes local false discovery rate (localFDR) estimates and semiparametric mixture density estimates from multi-dimensional inputs, which may include z-values, p-values, or raw data. The function utilizes a two-component semiparametric mixture model to estimate localFDR from z-values, incorporating Efron's empirical null principle and log-concave density estimation for the alternative distribution.
#'
#' @param z A matrix where each row represents a data point (z-values, p-values, or raw data).
#' @param init_p0 Initial value for the prior probability \( p_0 \) in the EM algorithm (default: 0.5).
#' @param tol Convergence threshold for the EM algorithm. The optimization stops when the maximum absolute difference between current and previous gamma values is smaller than \code{tol}. Specifically, if \eqn{ \text{max}_i |\gamma_i^{(k+1)} - \gamma_i^{(k)} | < \text{tol} }, for the k-th step, the algorithm terminates (default: 5e-6).
#' @param p_value Logical indicating if the input data are p-values (TRUE) or z-values/raw data (FALSE) (default: FALSE).
#' @param greater_alt A vector specifying whether the alternative distribution is greater (`TRUE`) or less (`FALSE`) than the null distribution for each dimension. For multidimensional data, this vector should match the number of columns in \code{z}.
#' @param min_iter Minimum number of iterations for the EM algorithm (default: 3).
#' @param max_iter Maximum number of iterations for the EM algorithm (default: 30).
#' @param thre_z The upper threshold of gamma used for log-concave estimates in the M-step of the EM algorithm (default: 1-1e-5).
#' @param Uthre_gam The upper threshold of gamma used to determine the stopping criteria for the EM algorithm (default: 0.99).
#' @param Lthre_gam The lower threshold of gamma used to determine the stopping criteria for the EM algorithm (default: 0.01).
#'
#' @return A list containing estimates from the semiparametric mixture model for the given data, including:
#'   \item{z}{Matrix where each row represents a data point.}
#'   \item{p0}{Prior probability for the null distribution.}
#'   \item{mu0, sig0}{Parameter estimates of the Gaussian (null) distribution, \( N(\mu_0, \sigma_0^2) \).}
#'   \item{f}{Probability estimates from the semiparametric mixture model.}
#'   \item{f1}{Probability estimates from the log-concave (alternative) distribution of the mixture model.}
#'   \item{F}{Cumulative density estimates from the mixture model.}
#'   \item{localFDR}{Local FDR estimates for the given data.}
#'   \item{FDR}{FDR estimates for the given data.}
#'   \item{iter}{Number of iterations performed by the EM algorithm.}
#'   \item{dim}{Dimension of the input data.}
#'   \item{greater_alt}{Vector indicating if the alternative distribution is greater (`TRUE`) or less (`FALSE`) than the null distribution for each dimension.}
#'
#' @export
#' 
SpMix <- function(z, init_p0 = 0.5, tol = 5.0e-5, p_value = FALSE,
                  greater_alt = NULL, min_iter = 3, max_iter = 50, 
                  thre_z = 0.5, Uthre_gam = 0.99, Lthre_gam = 0.01,
                  thre = 0.2)
{
  # *****************DEFINITION OF INTERNAL FUNCTIONS ******************
  
  # Internal Functions
  fit_normal_mixture_1d <- function(z, p0, tol = 5e-3, max.iter = 5) {
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
  
  fit_normal_mixture_nd <- function(z, p0, tol = 5e-3, max_iter = 5)
  {
    require(mvtnorm)
    
    k <- 0; converged <- 0
    z <- as.matrix(z)
    n <- nrow(z)
    d <- ncol(z)
    q <- qnorm(p0)
    z_scaled <- scale(z)
    ind0 <- ind1 <- rep(1, n)
    for (j in 1:d) {
      ind0 <- ind0*(z_scaled[,j] <= q)
      ind1 <- ind1*(z_scaled[,j] > -q)
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
  
  select_NE <- function(x, X){
    n <- nrow(X)
    p <- ncol(X)
    xx <- matrix(x, nrow = n, ncol = p, byrow = TRUE)
    ne.ind <- apply(1*(X >= xx), 1, prod)

    return((1:n)[ne.ind == 1])
  }

  monotone_fdr <- function(z, fdr) {
    sapply(seq_len(nrow(z)), function(i) max(fdr[select_NE(z[i, ], z)]))
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
      if(!greater_alt[j]) z[,j] <- -z[,j]
    }
  }
  
  ## Initial step: to fit normal mixture
  if (d == 1) {
    z = z[,1]
    # Initial step: Fit normal mixture
    nmEM <- fit_normal_mixture_1d(z, p0 = init_p0)
    p0 <- nmEM$p0
    f0 <- p0 * dnorm(z, mean = nmEM$mu0, sd = nmEM$sigma0)
    f1 <- dnorm(z, mean = nmEM$mu1, sd = nmEM$sigma1)
    f <- p0 * f0 + (1 - p0) * f1
  } else {
    ## Initial step: to fit normal mixture
    Params <- fit_normal_mixture_nd(z, p0 = init_p0)
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
      gam <- monotone_fdr(z, gam)
      
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
  
  if(!is.null(greater_alt)) {
    for(j in 1:d) {
      if(!greater_alt[j]){
        mu0[j] <- -mu0[j]
        z[,j] <- -z[,j]
      }
    }
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
      res <- list(z = z, p0 = p0, mu0 = mu0, sig0 = sig0, f1 = f1, f = f,
                  iter = k, log.likelihood = ell, lcd = lcd, dim = d,
                  posterior = cbind(gam, 1 - gam),
                  converged = converged)
    }
  }
  return(res)
}



