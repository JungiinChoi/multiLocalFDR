#' @importFrom mclust dmvnorm
#' @importFrom logcondens activeSetLogCon
#'
#' @title Local False Discovery Rate (localFDR) Estimation
#'
#' @description
#' `localFDR` estimates local false discovery rates using a semiparametric mixture model
#' via `SPMix()`, which combines Efron's empirical null principle with log-concave density estimation.
#'
#' @param z A numeric matrix or vector. Each row represents a data point (e.g., z-values or raw data).
#' @param tol Convergence threshold for the EM algorithm (default: 5e-6).
#' @param is_pvalue Logical; if `TRUE`, the input is treated as p-values and transformed (default: FALSE).
#' @param alternative Logical; if `TRUE`, assumes alternative distribution is greater than null (right-tailed); 
#'                    if `FALSE`, left-tailed; if `NULL`, determined automatically (default: NULL).
#' @param max_iter Maximum number of iterations for the EM algorithm (default: 30).
#' @param min_iter Minimum number of EM iterations before convergence is checked (default: 3).
#' @param Uthre_gam Upper threshold for gamma convergence (default: 0.99).
#' @param Lthre_gam Lower threshold for gamma convergence (default: 0.01).
#' @param thre Threshold for hypothesis rejection (used internally to determine the cutoff, default: 0.2).
#'
#' @return A numeric vector of estimated local FDR values for each observation.
#' 
#' @export
localFDR <- function(z, 
                     tol = 5e-6, 
                     is_pvalue = FALSE, 
                     alternative = NULL,
                     max_iter = 30, 
                     min_iter = 3,
                     Uthre_gam = 0.99, 
                     Lthre_gam = 0.01,
                     thre = 0.2) {
  
  # Run semiparametric mixture model
  fit <- SPMix(z = z,
               tol = tol,
               is_pvalue = is_pvalue,
               alternative = alternative,
               max_iter = max_iter,
               min_iter = min_iter,
               Uthre_gam = Uthre_gam,
               Lthre_gam = Lthre_gam,
               thre = thre)
  
  return(fit$localFDR)
}

