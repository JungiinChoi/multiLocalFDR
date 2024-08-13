#' @importFrom logcondens activeSetLogCon
#'
#' @title LocalFDR estimation for 1-dimensional z-values
#'
#' @description \code{sp.mix.1D} returns LocalFDR estimates and semiparametric
#' mixture density estimates for given 1-dimensional lists of z-values, which
#' are the probit-transformed p-values.
#' For the hypothesis testing \code{sp.mix.1D} uses a two-component semiparametric
#' mixture model to estimate the LocalFDR from the p-values. The two pillars of
#' the proposed approach are Efron's empirical null principle and log-concave
#' density estimation for the alternative distribution.
#'
#'
#' @param z Vector which each element indicates z-values, probit-transformed p-values.
#' @param tol Stopping criteria for the EM algorithm. If maximum absolute difference
#' of current and previous gamma value is smaller than tol,
#' i.e. \eqn{max_i |\gamma_i^{(k+1)}-\gamma_i^{(k)} <tol}, for k-th step, then optimization stops. (default: 5e-6)
#' @param max.iter Maximum number of iterations in the EM algorithm. (default: 30)
#' @param doplot Boolean parameter that draws histogram and fitted lines for estimated density,
#' if TRUE. (default: TRUE)
#' @param thre.localFDR Threshold of LocalFDR which is used for calculation of
#' false positive rate (FPR) and sensitivity with the threshold set to be
#' \eqn{fdr \leq thre.localFDR}. (default: 0.2)
#' @param thre.z Threshold value which only z-values smaller than thre.z
#' are used to compute the log-concave estimates f_1 in M-step.
#' @param Uthre.gam Upper threshold of gamma which are used to compute stopping criteria for the EM algorithm.
#' @param Lthre.gam Lower threshold of gamma which are used to compute stopping criteria for the EM algorithm.
#'
#'
#' @return Estimates of semiparametric mixture model for f for given z-values.
#'
#'   \item{p0}{Prior probability for null distribution}
#'   \item{mu0 sig0}{Parameter estimates of normal null distribution, N(mu0, sig0^2)}
#'   \item{f}{Probability estimates of semiparametric mixture model for each z-value point.}
#'   \item{localfdr}{LocalFDR estimates for given z-values}
#'   \item{iter}{Number of iterations of EM algorithm to compute LocalFDR.}
#'
#' @export
# Load required package
require(logcondens)

# 1D SP Mixture Model Function
sp.mix.1D <- function(z, p0.init = 0.5, thre = 0.5, tol = 5e-5, max.iter = 50) {
  
  # Ensure z is numeric
  z <- as.numeric(z)
  n <- length(z)
  zs <- sort(z)
  ell <- numeric(max.iter)
  
  # Initial step: Fit normal mixture
  nmEM <- normal.mixture.1d(z, p0 = p0.init)
  p0 <- nmEM$p0
  f0 <- p0 * dnorm(z, mean = nmEM$mu0, sd = nmEM$sigma0)
  f1 <- dnorm(z, mean = nmEM$mu1, sd = nmEM$sigma1)
  f <- p0 * f0 + (1 - p0) * f1
  
  # EM algorithm for SP mixture
  k <- 0
  converged <- FALSE
  
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
  
  list(
    z = z, p0 = p0, mu0 = mu0, sig0 = sig0, f1 = f1, f = f,
    iter = k, log.likelihood = ell, lcd = lcd,
    posterior = cbind(gam, 1 - gam),
    converged = converged
  )
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

# Fit SP Mixture Model for 1D data
fit.sp.mix.1D <- function(z, np.left = "No", doplot = TRUE, p0.init = 0.5, thre = 0.5, tol = 5e-4, max.iter = 50) {
  
  # Adjust data based on np.left
  zx <- if (np.left == "Yes") -z else z
  
  fit <- sp.mix.1D(zx, tol = tol, max.iter = max.iter, p0.init = p0.init, thre = thre)
  fit$f1 <- evaluateLogConDens(zx, fit$lcd)[, "density"]
  if (np.left == "Yes") fit$mu0 <- -fit$mu0
  
  # Plotting
  if (doplot) {
    par(mfrow = c(1, 1))
    hist(z,
         nclass = min(round(length(z) / 20), 200),
         probability = TRUE,
         col = "gray", border = "white",
         xlab = "", main = "",
         sub = substitute(
           paste(p[0], " = ", p0, ", ",
                 mu[0], " = ", mu0, ", ",
                 sigma[0], " = ", sigma0),
           list(p0 = round(fit$p0, 3),
                mu0 = round(fit$mu0, 3),
                sigma0 = round(fit$sig0, 3))),
         cex.sub = 0.8)
    
    zs <- sort(z)
    comp0 <- fit$p0 * dnorm(zs, mean = fit$mu0, sd = fit$sig0)
    comp1 <- (1 - fit$p0) * fit$f1[order(z)]
    
    lcolor <- c(4, 2, 3)
    ltype <- c(3, 3, 1)
    lwidth <- c(3, 3, 2)
    
    lines(zs, comp0 + comp1, col = lcolor[3], lwd = lwidth[3], lty = ltype[3])
    lines(zs, comp0, col = lcolor[1], lwd = lwidth[1], lty = ltype[1])
    lines(zs, comp1, col = lcolor[2], lwd = lwidth[2], lty = ltype[2])
    legend("topright", c("Normal", "Log-concave", "Marginal"), 
           col = lcolor, lty = ltype, lwd = lwidth, cex = 0.8)
  }
  
  return(fit)
}

# Local FDR Calculation for 1D data
local.fdr.1D <- function(z, tol = 5e-4, max.iter = 50, doplot = TRUE, p0.init = 0.9, thre = 0.2) {
  
  fit <- sp.mix.1D(z, tol = tol, max.iter = max.iter, p0.init = p0.init, thre = 1 - thre)
  critical <- min(z[fit$posterior[, 1] <= thre])
  
  # Plotting
  if (doplot) {
    par(mfrow = c(1, 1))
    hist(z,
         nclass = min(round(length(z) / 20), 200),
         probability = TRUE,
         col = "gray", border = "white",
         xlab = "", main = "",
         sub = substitute(
           paste(p[0], " = ", p0, ", ",
                 mu[0], " = ", mu0, ", ",
                 sigma[0], " = ", sigma0, ", ",
                 "critical-value = ", critical),
           list(p0 = round(fit$p0, 3),
                mu0 = round(fit$mu0, 3),
                sigma0 = round(fit$sig0, 3),
                critical = round(critical, 3))),
         cex.sub = 0.8)
    
    zs <- sort(z)
    comp0 <- fit$p0 * dnorm(zs, mean = fit$mu0, sd = fit$sig0)
    comp1 <- (1 - fit$p0) * fit$f1[order(z)]
    
    lcolor <- c(4, 2)
    ltype <- c(3, 3)
    lwidth <- c(3, 3)
    
    lines(zs, comp0, col = lcolor[1], lwd = lwidth[1], lty = ltype[1])
    lines(zs, comp1, col = lcolor[2], lwd = lwidth[2], lty = ltype[2])
    points(critical, 0, bg = "yellow", col = "magenta", pch = 25, cex = 1.5)
    legend("topright", c("Empirical Null", "Non-null (log-concave)"), 
           col = lcolor, lty = ltype, lwd = lwidth, cex = 0.8)
  }
  
  fit$critical.value <- critical
  fit$local.fdr <- fit$posterior[, 1]
  
  return(fit)
}

##############################
#                            #
#  Examples: To fit mixture  #
#                            #
##############################

N <- 1000
p0 <- 0.7
N0 <- rbinom(1, size = N, prob = p0)
membership <- c(rep(0, N0), rep(1, N - N0))

# Simulated Data: Normal + Normal
mu1 <- 4; sig1 <- 0.5
z <- c(rnorm(N0), rnorm(N - N0, mean = mu1, sd = sig1))
x <- seq(from = min(z), to = max(z), by = 0.01)
fit <- fit.sp.mix.1D(z)
lines(x, p0*dnorm(x), lty = 2)
lines(x, (1 - p0)*dnorm(x, mean = mu1, sd = sig1), lty = 2)
plot(fit$log.likelihood, type = "b", pch = 19, col = "navy",
     xlab = "Iteration", ylab = "Log-likelihood", main = "",
     cex.axis = 0.8, cex.lab = 0.8, cex.main = 0.9)
membership.cl <- 1*(fit$posterior[,2] >= 0.5)
table(membership, membership.cl)

# Simulated Data: Normal + Gamma
k <- 2; s <- 4; r <- 2
z <- c(rnorm(N0), rgamma(N - N0, shape = s, rate = r) + k)
x <- seq(from = min(z), to = max(z), by = 0.01)
fit <- fit.sp.mix.1D(z)
lines(x, p0*dnorm(x), lty = 2) 
lines(x + k, (1 - p0)*dgamma(x, shape = s, rate = r), lty = 2)
plot(fit$log.likelihood, type = "b", pch = 19, col = "navy",
     xlab = "Iteration", ylab = "Log-likelihood", main = "",
     cex.axis = 0.8, cex.lab = 0.8, cex.main = 0.9)
membership.cl <- 1*(fit$posterior[,2] >= 0.5)
table(membership, membership.cl)

# Old Faithful geyser data
fit <- fit.sp.mix.1D(faithful$waiting, np.left = "Yes")
plot(fit$log.likelihood, type = "b", pch = 19, col = "navy",
     xlab = "Iteration", ylab = "Log-likelihood", main = "",
     cex.axis = 0.8, cex.lab = 0.8, cex.main = 0.9)

# Carina data
carina <- read.table(file = "./dat/carina.dat")
x <- carina[carina$V8 + carina$V9 > 0,]
x <- x[x$V6 < 3,]
vel <- x$V4 # represents the Radial velocity data of stars in the Carina galaxy
fit <- fit.sp.mix.1D(vel, np.left = "Yes")
plot(fit$log.likelihood, type = "b", pch = 19, col = "navy",
     xlab = "Iteration", ylab = "Log-likelihood", main = "",
     cex.axis = 0.8, cex.lab = 0.8, cex.main = 0.9)


################################
#                              #
#  Examples: multiple testing  #
#                              #
################################

require(locfdr)
N <- 10000
p0 <- 0.9
N0 <- rbinom(1, size = N, prob = p0)
membership <- factor(c(rep("null", N0), rep("nonnull", N - N0)), 
                     levels = c("null", "nonnull"))

# Simulated Data: Normal + Normal
mu1 <- 3.5; sig1 <- 0.5
z <- c(rnorm(N0), rnorm(N - N0, mean = mu1, sd = sig1)) 
x <- seq(from = min(z), to = max(z), by = 0.01)
fit <- local.fdr.1D(z)
membership.cl <- factor(rep("null", N), levels = c("null", "nonnull"))
membership.cl[fit$local.fdr <= 0.2] <- "nonnull"
table(membership, membership.cl)

efron <- locfdr(z)
min(z[efron$fdr <= 0.2])
membership.cl <- factor(rep("null", N), levels = c("null", "nonnull"))
membership.cl[efron$fdr <= 0.2] <- "nonnull"
table(membership, membership.cl)

par(mfcol = c(2, 2))
plot(sort(z), fit$local.fdr[order(z)], type = "l",
     xlab = "z", ylab = "local fdr", main = "SpMix",
     cex.axis = 0.8, cex.lab = 0.8)
abline(h = 0.2, col = 2)
plot(sort(z), efron$fdr[order(z)], type = "l",
     xlab = "z", ylab = "local fdr", main = "Efron",
     cex.axis = 0.8, cex.lab = 0.8)
abline(h = 0.2, col = 2)
boxplot(fit$local.fdr ~ membership, ylab = "local fdr", xlab = "", main = "SpMix")
abline(h = 0.2, col = 2)
boxplot(efron$fdr ~ membership, ylab = "local fdr", xlab = "", main = "Efron")
abline(h = 0.2, col = 2)


# Simulated Data: Normal + Gamma
k <- 2; s <- 4; r <- 2
z <- c(rnorm(N0), rgamma(N - N0, shape = s, rate = r) + k)
fit <- local.fdr.1D(z)
membership.cl <- factor(rep("null", N), levels = c("null", "nonnull"))
membership.cl[fit$local.fdr <= 0.2] <- "nonnull"
table(membership, membership.cl)

efron <- locfdr(z)
min(z[efron$fdr <= 0.2])
membership.cl <- factor(rep("null", N), levels = c("null", "nonnull"))
membership.cl[efron$fdr <= 0.2] <- "nonnull"
table(membership, membership.cl)

par(mfcol = c(2, 2))
plot(sort(z), fit$local.fdr[order(z)], type = "l",
     xlab = "z", ylab = "local fdr", main = "SpMix",
     cex.axis = 0.8, cex.lab = 0.8)
abline(h = 0.2, col = 2)
plot(sort(z), efron$fdr[order(z)], type = "l",
     xlab = "z", ylab = "local fdr", main = "Efron",
     cex.axis = 0.8, cex.lab = 0.8)
abline(h = 0.2, col = 2)
boxplot(fit$local.fdr ~ membership, ylab = "local fdr", xlab = "", main = "SpMix")
abline(h = 0.2, col = 2)
boxplot(efron$fdr ~ membership, ylab = "local fdr", xlab = "", main = "Efron")
abline(h = 0.2, col = 2)
