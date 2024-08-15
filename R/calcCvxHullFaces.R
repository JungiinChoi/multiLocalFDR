#' @title Calculates the convex hull of X
#'
#' @description \code{calcCvxHullFaces} returns the parameters of the convex
#'   hull of X: The indices of points in X that compose the faces of conv(X)
#'   (one row per face) as well as the parameters of the hyperplanes that
#'   describe these faces.
#'
#' @param X Set of data points (one sample per row)
#' @return A list consisting of parameters describing the convex hull of X:
#'   \item{cvh}{A matrix of indices where each row constitutes one face}
#'   \item{ACVH}{A matrix where each row constitutes the normal vector of a
#'   face} \item{bCVH}{A vector where each entry constitutes the offset of the
#'   hyperplane for a face}

calcCvxHullFaces <-function(X) {

  d <- dim(X)[2]

  # the 1-D case is easy
  if (d == 1) {
    A <- matrix(c(1,-1), ncol=1)
    b <- matrix(c(max(X), -min(X)), ncol=1)
    cvh <- matrix(c(which.min(X), which.max(X)), nrow = 1)
  } else {
    mu <-  apply(X, 2, mean)
    cvh <- geometry::convhulln(X);

    # init param matrix/vector
    A <- array(0, c(dim(cvh)[1], d))
    b <- array(0, c(dim(cvh)[1], 1))
    # iterate over all faces of \conv(X)
    for (i in 1:dim(cvh)[1]) {
      B <- X[cvh[i,1:d-1], ]-X[cvh[i,2:d], ]
      if (d == 2) {
        A[i, ] <- MASS::Null(B)
      } else {
        A[i, ] <- MASS::Null(t(B))
      }

      # test orientation with sample mean and reverse if necessary
      if (A[i, ] %*% (X[cvh[i,1], ] - mu) < 0) {
        A[i, ] <- -A[i, ]
      }
      b[i] <- A[i, ] %*% X[cvh[i,1], ]
    }
  }

  r <- list(
    cvh = cvh,
    ACVH = A,
    bCVH = b
  )
  return(r)
}
