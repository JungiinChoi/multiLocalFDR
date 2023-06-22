callNewtonBFGSLC <- function (X, w, params, paramsKernel, cvhParams, gamma=1000, verbose=0, intEps = 1e-3, objEps = 1e-7, offset = 1e-1, maxIter = 1e4) {
  n <- dim(X)[1]
  d <- dim(X)[2]

  r <- .C(  "newtonBFGSLC",
              as.double(X),
              as.double(w),
              as.double(c(apply(X, 2, min), apply(X, 2, max))),
              params = as.double(params),
              as.double(paramsKernel),
              lenP = as.integer(length(params)),
              as.integer(length(paramsKernel)),
              as.integer(d),
              as.integer(n),
              as.double(cvhParams$ACVH),
              as.double(cvhParams$bCVH),
              as.integer(length(cvhParams$bCVH)),
              as.double(intEps),
              as.double(objEps),
              as.double(offset),
              as.integer(verbose),
              as.double(gamma),
              as.integer(maxIter))
  # res$lenP denotes the number of active parameters in res$params
  optParams <- r$params[1:r$lenP]

  nH <- r$lenP/(d+1) # number of hyperplanes
  a <- matrix(optParams[1:(d * nH)],nH,d)
  b <- optParams[(d*nH+1):length(optParams)]
  return(list("a" = a, "b" = b))
}

callNewtonBFGSLInitC <- function (X, w, params, ACVH, bCVH) {
  n <- dim(X)[1]
  d <- dim(X)[2]

  r <- .C("newtonBFGSLInitC",
          as.double(X),
          as.double(w),
          as.double(c(apply(X, 2, min),apply(X, 2, max))),
          params = as.double(params),
          as.integer(d),
          as.integer(length(params)),
          as.integer(n),
          as.double(ACVH),
          as.double(bCVH),
          as.integer(length(bCVH)),
          as.double(1e-3),
          as.double(1e-6),
          logLike = as.double(matrix(0, 1)))
  return(r)
}

callCalcExactIntegralC <- function(X, y, cvh, filter, eps) {

  idxCVH <- setdiff(unique(as.vector(cvh)), which(!filter))
  P <- matrix(c(X[filter, ], y[filter]), nrow = length(filter))
  Q <- matrix(c(X[idxCVH, ], rep(min(y[idxCVH]) - 1, length(idxCVH))), nrow = length(idxCVH))

  X = X[filter, ]
  y = y[filter]

  n <- dim(X)[1]
  d <- dim(X)[2]

  T <- geometry::convhulln(rbind(P, Q))
  T <- T[!(apply(T, 1, max) > n), ]

  # if there is only one simplex
  if (class(T)[1] != "matrix") {
    T <- t(as.matrix(T))
  }

  r <- .C(  "calcExactIntegralC",
          as.double(t(X)),
          y = as.double(y),
          as.integer(t(T-1)),
          as.integer(dim(T)[1]),
          as.integer(n),
          as.integer(d),
          as.double(1),
          as.double(eps),
          a = as.double(matrix(0, dim(T)[1] * d)),
          b = as.double(matrix(0, dim(T)[1])))
  r$a = t(matrix(r$a, d, length(r$a) / d))
  return(r)
}


#' @title Checks if AVX is active
#'
#' @description \code{compilationInfo} is a wrapper to a C function, which prints information
#'   about whether AVX extensions where activated during compilation. Active AVX speeds up
#'   computations significantly (by a factor of 4 to 8 roughly). More informations about how
#'   AVX can be enabled during installation can be found in documentation.

compilationInfo <- function() {
  invisible(.C ("printAVXInfo"))
}
