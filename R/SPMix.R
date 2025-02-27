#' @useDynLib multiLocalFDR, .registration = TRUE
#' @importFrom mclust dmvnorm
#' @importFrom mvtnorm pmvnorm
#' @importFrom logcondens activeSetLogCon
#' @importFrom logcondens evaluateLogConDens
#' @importFrom graphics legend
#' @import stats
#' @importFrom Rcpp sourceCpp
#' 
#' @title Semiparametric Mixture Model for Local False Discovery Rate Estimation
#' 
#' @description 
#' The `SpMix` function estimates the local false discovery rate (localFDR) and semiparametric mixture density 
#' from multi-dimensional inputs such as z-values, p-values, or raw data. It employs a two-component 
#' semiparametric mixture model, integrating Efron's empirical null principle and log-concave density 
#' estimation for the alternative distribution.
#' 
#' @param z A numeric matrix where each row represents a data point (z-values, p-values, or raw data).
#' @param initial_p0 Initial prior probability \( p_0 \) for the null distribution in the EM algorithm (default: 0.5).
#' @param tol Convergence threshold for the EM algorithm. The algorithm terminates when 
#' the maximum absolute difference between consecutive gamma values falls below \code{tol} (default: 5e-6).
#' @param is_pvalue Logical; if TRUE, the input is assumed to be p-values and transformed accordingly (default: FALSE).
#' @param alternative A logical vector indicating whether the alternative distribution is greater (`TRUE`) 
#' or less (`FALSE`) than the null distribution in each dimension. It must match the number of columns in `z`.
#' @param min_iter Minimum number of iterations for the EM algorithm (default: 3).
#' @param max_iter Maximum number of iterations for the EM algorithm (default: 50).
#' @param thre_z Threshold for gamma values used in log-concave estimation during the M-step of the EM algorithm (default: 0.5).
#' @param Uthre_gam Upper threshold for gamma to determine stopping criteria in the EM algorithm (default: 0.99).
#' @param Lthre_gam Lower threshold for gamma to determine stopping criteria in the EM algorithm (default: 0.01).
#' 
#' @return A list containing estimates from the semiparametric mixture model, including:
#'   \item{z}{Matrix of input data points.}
#'   \item{p0}{Estimated prior probability for the null distribution.}
#'   \item{mu0, sig0}{Parameters of the Gaussian null distribution, \( N(\mu_0, \sigma_0^2) \).}
#'   \item{f}{Estimated mixture density.}
#'   \item{f1}{Estimated log-concave alternative density.}
#'   \item{localFDR}{Estimated local false discovery rates.}
#'   \item{FDR}{Estimated false discovery rates.}
#'   \item{iter}{Total number of EM iterations performed.}
#'   \item{dim}{Dimensionality of the input data.}
#'   \item{alternative}{Vector indicating whether the alternative distribution is greater (`TRUE`) or less (`FALSE`) than the null distribution.}
#' 
#' @export
#' 
SpMix <- function(z, initial_p0 = 0.5, tol = 5.0e-5, is_pvalue = FALSE,
                  alternative = NULL, min_iter = 3, max_iter = 50, 
                  thre_z = 0.5, Uthre_gam = 0.99, Lthre_gam = 0.01,
                  thre = 0.2)
{
  # *****************DEFINITION OF INTERNAL FUNCTIONS ******************
  
  # Internal Functions
  fit_normal_mixture_1d <- function(z, p0, tol = 5e-3, max_iter = 5) {
    q0 <- quantile(z, probs = p0)
    
    # Compute initial estimates
    mu0 <- mean(z[z <= q0])
    sig0 <- sd(z[z <= q0])
    mu1 <- mean(z[z > q0])
    sig1 <- sd(z[z > q0])
    
    f0 <- dnorm(z, mean = mu0, sd = sig0)
    f1 <- dnorm(z, mean = mu1, sd = sig1)
    ell <- sum(log(p0 * f0 + (1 - p0) * f1))
    
    for (k in 1:max_iter) {
      gam <- p0 * f0 / (p0 * f0 + (1 - p0) * f1)
      
      # M-step
      p0 <- mean(gam)
      mu0 <- sum(z * gam) / sum(gam)
      sig0 <- sqrt(sum(gam * (z - mu0)^2) / sum(gam))
      f0 <- dnorm(z, mean = mu0, sd = sig0)
      
      mu1 <- sum(z * (1 - gam)) / sum(1 - gam)
      sig1 <- sqrt(sum((1 - gam) * (z - mu1)^2) / sum(1 - gam))
      f1 <- dnorm(z, mean = mu1, sd = sig1)
      
      new_ell <- sum(log(p0 * f0 + (1 - p0) * f1))
      if (k >= 3 && abs(new_ell - ell) <= tol) break
      ell <- new_ell
    }
    
    list(p0 = p0, mu0 = mu0, sigma0 = sig0, mu1 = mu1, sigma1 = sig1)
  }
  
  fit_normal_mixture_nd <- function(z, p0, tol = 5e-3, max_iter = 5) {
    require(mvtnorm)
    
    z <- as.matrix(z)
    n <- nrow(z)
    d <- ncol(z)
    q <- qnorm(p0)
    
    # Identify indices for the null and alternative distributions
    ind0 <- rowSums(scale(z) <= q) == d
    ind1 <- rowSums(scale(z) > -q) == d
    
    # Compute initial estimates
    mu0 <- colMeans(z[ind0, , drop = FALSE])
    sig0 <- cov(z[ind0, , drop = FALSE])
    mu1 <- colMeans(z[ind1, , drop = FALSE])
    sig1 <- cov(z[ind1, , drop = FALSE])
    
    f0 <- dmvnorm(z, mu0, sig0)
    f1 <- dmvnorm(z, mu1, sig1)
    f <- p0 * f0 + (1 - p0) * f1
    ell <- mean(log(f))
    
    for (k in 1:max_iter) {
      gam <- p0 * f0 / f
      
      # M-step
      p0 <- mean(gam)
      mu0 <- colSums(z * gam) / sum(gam)
      dev0 <- (z - mu0) * sqrt(gam)
      sig0 <- t(dev0) %*% dev0 / sum(gam)
      f0 <- dmvnorm(z, mu0, sig0)
      
      mu1 <- colSums(z * (1 - gam)) / sum(1 - gam)
      dev1 <- (z - mu1) * sqrt(1 - gam)
      sig1 <- t(dev1) %*% dev1 / sum(1 - gam)
      f1 <- dmvnorm(z, mu1, sig1)
      
      f <- p0 * f0 + (1 - p0) * f1
      new_ell <- mean(log(f))
      
      if (k >= 3 && abs(new_ell - ell) <= tol) break
      ell <- new_ell
    }
    
    list(p0 = p0, mu0 = mu0, sig0 = sig0, mu1 = mu1, sig1 = sig1)
  }
  
  select_NE <- function(x, X) {
    return(which(rowSums(X >= matrix(x, nrow = nrow(X), ncol = ncol(X), byrow = TRUE)) == ncol(X)))
  }
  
  monotone_fdr <- function(z, fdr) {
    apply(z, 1, function(row) max(fdr[select_NE(row, z)]))
  }
  
  # ******************* MAIN FUNCTION *******************************
  
  z <- as.matrix(z)
  n <- nrow(z)
  d <- ncol(z)
  ell <- rep(NA, max_iter)
  
  if (is_pvalue) {
    z <- qnorm(ifelse(is.null(alternative), 1 - z, z))
  } else {
    if (ncol(z) == 1) z <- scale(z)
  }
  
  if (!is.null(alternative)) {
    z[, !alternative] <- -z[, !alternative]
  }
  
  ## Initial step: to fit normal mixture
  if (d == 1) {
    z = z[,1]
    # Initial step: Fit normal mixture
    nmEM <- fit_normal_mixture_1d(z, p0 = initial_p0)
    p0 <- nmEM$p0
    f0 <- p0 * dnorm(z, mean = nmEM$mu0, sd = nmEM$sigma0)
    f1 <- dnorm(z, mean = nmEM$mu1, sd = nmEM$sigma1)
    f <- p0 * f0 + (1 - p0) * f1
  } else {
    ## Initial step: to fit normal mixture
    Params <- fit_normal_mixture_nd(z, p0 = initial_p0)
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
      lcd <- fmlcd(X = z[which_z,], w = weight)
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
  
  if(!is.null(alternative)) {
    for(j in 1:d) {
      if(!alternative[j]){
        mu0[j] <- -mu0[j]
        z[,j] <- -z[,j]
      }
    }
  }

  # return results

  if (is_pvalue) {
    res <- list(z = z, p0 = p0, mu0 = mu0, sig0 = sig0, f = f, f1 = f1,  
                log.likelihood = ell,
                localFDR = gam, posterior = cbind(gam, 1 - gam),
                iter = k, dim = d, alternative = alternative,
                converged = converged)
  } else {
    if (d == 1){
      res <- list(z = z*raw_sd + raw_mean, p0 = p0, mu0 = mu0*raw_sd + raw_mean, 
                  sig0 = sig0*raw_sd, f = f/raw_sd, f1 = f1/raw_sd, 
                  log.likelihood = ell,
                  localFDR = gam, posterior = cbind(gam, 1 - gam),
                  iter = k, dim = d, alternative = alternative,
                  converged = converged)
    } else {
      res <- list(z = z, p0 = p0, mu0 = mu0, sig0 = sig0, f1 = f1, f = f,
                  iter = k, log.likelihood = ell, lcd = lcd, dim = d,
                  posterior = cbind(gam, 1 - gam), localFDR = gam,
                  converged = converged)
    }
  }
  class(res) <- "SpMix"
  return(res)
}



