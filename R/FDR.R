#' @importFrom fmlogcondens fmlcd
#' @importFrom mclust dmvnorm
#' @importFrom logcondens activeSetLogCon
#'
#' @title FDR estimation for given z-values / p-values
#'
#' @description \code{FDR} returns FDR estimates for given multi-dimensional lists of z-values / p-values.
#' \code{FDR} imports \code{SpMix} for a two-component semiparametric
#' mixture model to estimate the FDR from the z-values / p-values.
#'
#' @param z Matrix which column indicates z-values, probit-transformed p-values.
#' @param tol Stopping criteria for the EM algorithm. If maximum absolute difference
#' of current and previous gamma value is smaller than tol,
#' i.e. \eqn{max_i |\gamma_i^{(k+1)}-\gamma_i^{(k)} <tol}, for k-th step,
#' then optimization stops. (default: 5e-6)
#' @param p_value If TRUE, input are p-values. If FALSE, input are z-values. (default: FALSE)
#' @param local If TRUE, \code{FDR} returns localFDR estimates for given z-values or p-values. IF FALSE, \code{FDR} returns FDR estimates. (default: FALSE)
#' @param alternative a character string specifying the alternative hypothesis, must be one of "two.sided", "greater" (default) or "less". You can specify just the initial letter. (default: "greater")
#' @param max_iter Maximum number of iterations in the EM algorithm. (default: 30)
#' @param mono If TRUE, FDR is in ascending order of z-values. (default: TRUE)
#' @param thre_z Threshold value which only z-values smaller than thre.z
#' are used to compute the log-concave estimates f_1 in M-step.
#' @param Uthre_gam Upper threshold of gamma which are used to compute stopping criteria for the EM algorithm.
#' @param Lthre_gam Lower threshold of gamma which are used to compute stopping criteria for the EM algorithm.
#'
#'
#' @return Estimates of FDR or localFDR for given z-values / p-values.
#'
#'   \item{FDR}{FDR estimates for given z-values / p-values}
#'
#' @export

FDR <- function(z, tol = 5e-6, p_value = FALSE, alternative = "greater", max_iter = 30,
                     mono = TRUE, thre_z = 0.9, Uthre_gam = 0.9, Lthre_gam = 0.01)
{
  SpMixParams <- SpMix(z, tol, p_value, alternative, max_iter, mono, thre_z,
                       Uthre_gam, Lthre_gam)

  z <- as.matrix(z)
  n <- dim(z)[1]
  d <- dim(z)[2]

  if (p_value) {
    if (alternative == "greater" | alternative == "g") {
      z = qnorm(1-z)
    } else {
      z = qnorm(z)
    }
  }

  if (d == 1) {
    F0 <- pnorm(z, SpMixParams$mu0, SpMixParams$sig0)
    F1 <- ecdf(SpMixParams$f1)
    p0 <- SpMixParams$p0
    FDR <- p0 * F0 / (p0 * F0 + (1-p0) * F1)
  }

  res_FDR <- list(F0 = F0, F1 = F1, FDR = FDR)
  return(res_FDR)
}
