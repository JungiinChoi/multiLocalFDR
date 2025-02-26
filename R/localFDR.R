#' @importFrom fmlogcondens fmlcd
#' @importFrom mclust dmvnorm
#' @importFrom logcondens activeSetLogCon
#'
#' @title Local False Discovery Rate (localFDR) Estimation
#'
#' @description 
#' `localFDR` computes local false discovery rate (localFDR) estimates for multi-dimensional input data 
#' (raw data, z-values, or p-values). It uses the semiparametric mixture model from `SpMix()`, 
#' which integrates Efron's empirical null principle and log-concave density estimation.
#'
#' @param z A numeric matrix where each row represents a data point (z-values, p-values, or raw data).
#' @param tol Convergence threshold for the EM algorithm. The optimization stops when 
#' the maximum absolute difference between consecutive gamma values is below `tol`. (default: 5e-6)
#' @param is_pvalue Logical; if `TRUE`, the input is assumed to be p-values and transformed accordingly. (default: `FALSE`)
#' @param alternative A logical vector specifying whether the alternative distribution is greater (`TRUE`) 
#' or less (`FALSE`) than the null distribution in each dimension. Must match the number of columns in `z`. 
#' (default: `NULL`, assumes `greater` for all dimensions)
#' @param max_iter Maximum number of iterations for the EM algorithm. (default: 30)
#' @param monotone Logical; if `TRUE`, ensures localFDR values are non-decreasing along the z-values. (default: `TRUE`)
#' @param Uthre_gam Upper threshold for gamma to determine stopping criteria in the EM algorithm. (default: 0.99)
#' @param Lthre_gam Lower threshold for gamma to determine stopping criteria in the EM algorithm. (default: 0.01)
#'
#' @return A numeric vector of estimated local false discovery rates (localFDR).
#' 
#' @export
localFDR <- function(z, tol = 5e-6, is_pvalue = FALSE, alternative = NULL, 
                     max_iter = 30, monotone = TRUE, Uthre_gam = 0.99, Lthre_gam = 0.01) 
{
  # Run the semiparametric mixture model
  SpMixParams <- SpMix(z = z, initial_p0 = 0.5, tol = tol, 
                       is_pvalue = is_pvalue, alternative = alternative, 
                       min_iter = 3, max_iter = max_iter, thre_z = 0.5, 
                       Uthre_gam = Uthre_gam, Lthre_gam = Lthre_gam, thre = 0.2)
  
  # Return only localFDR values
  return(SpMixParams$localFDR)
}
