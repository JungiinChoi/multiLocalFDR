#' @title Parameter Initialization Based on a Kernel Density
#'
#' @description \code{paramFitKernelDensity} first fits a kernel density to a
#'   sample X with weight vector w. It then calculates the parameters of the
#'   piecewise linear function defined to be the upper convex hull of
#'   (X,log(y)).
#'
#' @param X Set of data points (one sample per row)
#' @param w Vector with weights for X (\code{sum(w) == 1})
#' @param cvh Matrix where each row is a set of indices of points in X
#'   describing one face of conv(X)
#' @param h Scalar parameter that governs the Gaussian kernel
#'
#' @return A list containing the description of the upper convex hull of
#'   (X,log(y)) in term of hyperplane parameters: \item{a}{A matrix where each
#'   row constitutes the normal vector of a face} \item{b}{A vector where each
#'   entry constitutes the offset of a face}
#'
#' @example R/Examples/paramFitKernelDensity

paramFitKernelDensity <- function(X, w, cvh, h = apply(X, 2, sd) * n ** (-1 / (d+4))) {
  n <- dim(X)[1]
  d <- dim(X)[2]

  # C-function that calculates the kernel density
  r <- .C("calcKernelDens", X = t(X), sampleWeights = w, yT = matrix(0, n), h = h, n = n, d = d)
  y = log(r$yT)

  # find upper convex hull of X and y
  finiteVals <- is.finite(y)
  # analytically normalize log-concave density
  r <- callCalcExactIntegralC(X, y, cvh, finiteVals, 1e-2)

  return(c(r$a, r$b))
  }
