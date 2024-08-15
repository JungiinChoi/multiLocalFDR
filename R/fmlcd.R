#' @title Estimates a Log-Concave Density
#'
#' @description \code{fmlcd} returns a MLE estimate of a log-concave
#'   density for \code{X}. After obtaining an initial parameter estimate the MLE
#'   objective with log-concavity and normalization constraint is optimized
#'   using a quasi-Newton approach for large scale optimization (BFGS-L). The
#'   logarithm of the optimal density f(x) is a piecewise-linear function. Its
#'   parametrization in terms of a set of hyperplanes is returned.
#'
#' @param X Matrix of data points (one sample per row)
#' @param w Vector of sample weights (default: uniform weights)
#' @param init String that sets the initialization approach. 'kernel' based on
#'   kernel density, 'smooth' based on smooth log-concave density, '' compares
#'   both and takes the optimal one. (default: '')
#' @param verbose Int determining the verboseness of the code; 0 = no output
#'   to 3. (default: 0)
#' @param intEps Stopping criteria for the numerical integration accuracy
#'   Optimization stops if both measurements are smaller than intEps and objEps
#'   Modification of this value is not recommended. (default: 1e-3)
#' @param objEps Stopping criteria for the change in the objective function
#'   Optimization stops if both measurements are smaller than intEps and objEps
#'   Modification of this value is not recommended. (default: 1e-7)
#' @param offset Smaller values correspond to slower hyperplane reduction.
#'   Offset should be a value smaller than 1. Modification of this value is not
#'   recommended. (default: 1e-1)
#' @param maxIter Number of iterations in the main optimization (default: 1e4)
#'
#' @return Parametrization of f(x) in terms of hyperplanes and function
#'   evaluations y = log(f(x)) \item{aOpt, bOpt}{Analytically normalized
#'   parameters of f(x).} \item{logLike}{Log-likelihood of f(x)} \item{logMLE}{Vector
#'   with values y_i = log(f(X_)) of the normalized density (\eqn{logLike =
#'   \sum(y_i)}).} \item{aOptSparse, bOptSparse}{Sparse parametrization
#'   normalized on the integration grid.}
#'
#' @example R/Examples/correctIntegral

fmlcd <- function(X, w=rep(1/nrow(X),nrow(X)), init='', verbose=0, intEps = 1e-3, objEps = 1e-7, offset = 1e-1, maxIter = 1e4) {

  gamma = 1000
  n <- dim(X)[1]
  d <- dim(X)[2]

  # Do some sanity checks on the inputs
  if (n <= d) {
    stop("There must be at least d+1 data points.")
  }

  if (length(w) != n) { stop("There must be exactly one weight for each X_i.") }

  if(sum(w <= 0)) stop("All weights must be strictly positive.")

  if (offset < 0) {
    offset <- 1e-1
    warning("Offset values smaller than zero are not allowed, set to 0.1.")
  }

  if (intEps <= 0) {
    intEps <- 1e-4
    warning("intEps must be larger than zero, set to 1e-4.")
  }

  if (objEps <= 0) {
    objEps <- 1e-7
    warning("objEps must be larger than zero, set to 1e-7.")
  }

  # taken from LogConcDEAD
  if (is.matrix(X)==FALSE) {
    if(is.numeric(X)) {
      X <- matrix(X, ncol=1)
    }
    else {
      stop("X must be a numeric vector or matrix.")
    }
  }


  ## renormalize weights to sum to one
  w <- w/sum(w)

  ## substract the mean from X
  mu = apply(X,2,mean)
  X = X - t(matrix(rep(mu,n), d, n))

  ## obtain faces of convex hull of X
  cvhParams <- calcCvxHullFaces(X)

  paramsKernel <- double()
  if (init == 'kernel') {
    if (verbose > 0) {
      print("Use a kernel density estimation for initialization.")
    }
      params <- paramFitKernelDensity(X, w, cvhParams$cvh)
  } else if (init == 'smooth') {
      if (verbose > 0) {
        print("Use a smooth log-concave density estimation for initialization.")
      }
      params <- paramFitGammaOne(X, w, cvhParams$ACVH, cvhParams$bCVH, cvhParams$cvh)
  } else {
    if (n < 2500) {
      if (verbose > 0) {
        print("Choose between a kernel density and a smooth log-concave density estimation for initialization.")
      }
      paramsKernel <- paramFitKernelDensity(X, w, cvhParams$cvh)
      params <- paramFitGammaOne(X, w, cvhParams$ACVH, cvhParams$bCVH, cvhParams$cvh)
    } else {
      if (verbose > 0) {
        print("Use a smooth log-concave density estimation for initialization.")
      }
      params <- paramFitGammaOne(X, w, cvhParams$ACVH, cvhParams$bCVH, cvhParams$cvh)
    }
  }

  # The larger vector has to be params (due to C conventions)
  if (length(paramsKernel) > length(params)) {
    paramsTmp <- params
    params <- paramsKernel
    paramsKernel <- params
   }

  # call C code that optimizes the MLE objective with initial parameters params and paramsKernel
  # optimal set of parameters choosen inside the function
  res <- callNewtonBFGSLC(X, w, params, paramsKernel, cvhParams, gamma, verbose, intEps, objEps, offset)
  result <- correctIntegral(X, mu, res$a, res$b, cvhParams$cvh);
  logLike = result$logMLE * t(w) * length(w);

  # update convex hull parameters for true X
  #r$bCVH <- r$bCVH + r$ACVH %*% t(mu)

  return(c(result,list("logLike" = logLike)))
}
