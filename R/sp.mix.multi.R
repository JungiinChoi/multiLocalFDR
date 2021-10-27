#' @importFrom fmlogcondens fmlcd
#' @importFrom mclust dmvnorm
#'
#'
#' @title LocalFDR estimation for multi-dimensional z-values
#'
#' @description \code{sp.mix.multi} returns LocalFDR estimates and semiparametric
#' mixture density estimates for given multi-dimensional lists of z-values, which
#' are the probit-transformed p-values.
#' For the hypothesis testing \code{sp.mix.multi} uses a two-component semiparametric
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
sp.mix.multi <- function(z, tol = 5e-6, max.iter = 30, mono = TRUE, thre.z = 0.9, Uthre.gam = 0.9, Lthre.gam = 0.01)
  # FOR MULTIVARIATE CASE ONLY
{

  z <- as.matrix(z)
  n <- dim(z)[1]

  ## Initial step: to fit normal mixture
  nmEM <- normal.mixture(z)
  p.0 <- nmEM$p.0
  mu.0 <- nmEM$mu.0
  sig.0 <- nmEM$Sigma.0
  f1.tilde <- dmvnorm(z, nmEM$mu.1, nmEM$Sigma.1)
  gam <- f <- rep(0, n)

  ## EM-step
  k <- 0; converged <- 0
  while ( (k < 3)|((k < max.iter) & (!converged)) ) {
    k <- k + 1
    ## E-step
    tmp <- p.0*dmvnorm(z, mu.0, sig.0)
    new.f <- tmp + (1-p.0)*f1.tilde
    new.gam <- tmp/new.f
    if(mono) new.gam <- MonotoneFDR(z, new.gam)

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

normal.mixture <- function(z, tol = 5e-3, max.iter = 10)
{

  k <- 0; diff <- 100
  z <- as.matrix(z)
  if ( dim(z)[2] == 1 ) m.dist <- z/sd(z) else m.dist <- mahalanobis(z, center = rep(0, dim(z)[2]), cov = cov(z))

  p.0 <- mean(m.dist <= 1.65)
  mu.0 <- rep(0, dim(z)[2])
  sig.0 <- diag(1, dim(z)[2])
  f.0 <- dmvnorm(z, mean = mu.0, sigma = sig.0)

  mu.1 <- apply(z[m.dist > 1.65,], 2, mean)
  sig.1 <- cov(z[m.dist > 1.65,])
  f.1 <- dmvnorm(z, mean = mu.1, sigma = sig.1)

  while ( (k < 3)|((k < max.iter) & (diff > tol)) ) {
    k <- k + 1

    ## E-step
    term1 <- p.0*f.0
    term2 <- term1 + (1-p.0)*f.1
    gam <- term1/term2

    ## M-step
    new.p.0 <- mean(gam)
    new.mu.0 <- as.vector(t(z)%*%gam)/sum(gam)
    dev <- (z - new.mu.0)*sqrt(gam)
    new.sig.0 <- t(dev)%*%dev/sum(gam)
    f.0 <- dmvnorm(z, mean = new.mu.0, sigma = new.sig.0)
    new.mu.1 <- as.vector(t(z)%*%(1 - gam))/sum(1 - gam)
    dev <- (z - new.mu.1)*sqrt(1 - gam)
    new.sig.1 <- t(dev)%*%dev/sum(1 - gam)
    f.1 <- dmvnorm(z, mean = new.mu.1, sigma = new.sig.1)

    ## Update
    diff <- max(abs((new.mu.0 - mu.0)),
                abs((new.sig.0 - sig.0)),
                abs((new.mu.1 - mu.1)),
                abs((new.sig.1 - sig.1)),
                abs((new.p.0 - p.0)))
    p.0 <- new.p.0
    mu.0 <- new.mu.0
    sig.0 <- new.sig.0
    mu.1 <- new.mu.1
    sig.1 <- new.sig.1
  }
  return(list(iter = k, p.0 = p.0, mu.0 = mu.0, Sigma.0 = sig.0, mu.1 = mu.1, Sigma.1 = sig.1))
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
