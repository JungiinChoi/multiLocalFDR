#' @importFrom fmlogcondens fmlcd
#' @importFrom mclust dmvnorm
#' @importFrom logcondens activeSetLogCon
#'
#' @title localFDR estimation for given data
#'
#' @description \code{localFDR} returns localFDR estimates for given multi-dimensional lists of raw data, z-values, or p-values.
#' \code{localFDR} imports \code{SPMix} for a two-component semiparametric
#' mixture model to estimate the localFDR.
#'
#' @param z Matrix which column indicates z-values, probit-transformed p-values.
#' @param tol Stopping criteria for the EM algorithm. If maximum absolute difference
#' of current and previous gamma value is smaller than tol,
#' i.e. \eqn{max_i |\gamma_i^{(k+1)}-\gamma_i^{(k)} <tol}, for k-th step,
#' then optimization stops. (default: 5e-6)
#' @param p_value If TRUE, input are p-values. If FALSE, input are z-values. (default: FALSE)
#' @param alternative a character string specifying the alternative hypothesis, must be one of "two.sided", "greater" (default) or "less". You can specify just the initial letter. (default: "greater")
#' @param max_iter Maximum number of iterations in the EM algorithm. (default: 30)
#' @param mono If TRUE, FDR is in ascending order of z-values. (default: TRUE)
#' @param Uthre_gam Upper threshold of gamma which are used to compute stopping criteria for the EM algorithm.
#' @param Lthre_gam Lower threshold of gamma which are used to compute stopping criteria for the EM algorithm.
#'
#'
#' @return Estimates of localFDR for given z-values / p-values.
#'   \item{localFDR}{local FDR estimates for given z-values / p-values}
#'
#' @export
#'

localFDR <- function(z, tol = 5e-6, p_value = FALSE, alternative = "greater", max_iter = 30,
                     mono = TRUE, Uthre_gam = 0.99, Lthre_gam = 0.01)
{
  SPMixParams <- SPMix(z, tol, p_value, alternative, max_iter, mono, 
                       Uthre_gam, Lthre_gam)

  return(SPMixParams$localFDR)
}
