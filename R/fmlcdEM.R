#' @title Estimates a Log-Concave Mixture Density
#'
#' @description \code{fmlcdEM} Utilizes the EM approach to obtain a mixture
#'   of log-concave densities. Utilizes Gaussian hierarchical clustering to
#'   initilize the posterior probabilities of class affiliation (as proposed
#'   by the package LogConcDEAD by Cule et al.).
#'
#' @param X Matrix of data points (one sample per row)
#' @param K Number of latent variables (default: 2)
#' @param posterior Matrix with posterior probabilities for class affiliation;
#'   Initialized if not provided using a Gaussian hierarchical clustering.
#' @param verbose Int determining the verboseness of the code; 0 = no output
#'   to 3. (default: 0)
#' @param maxIter Maximal number of EM iterations. (default: 50)
#'
#' @return Parametrization of the mixture density
#'   \item{params}{List of length K, where each entry contains the hyperplane for one
#'   density}
#'   \item{densEst}{Matrix where each row contains the marginal distribution p(x)}
#'   \item{tau}{Marginal distribution over the latent variable p(z)}

fmlcdEM <- function(X, K = 2, posterior, verbose=0, maxIter = 50) {
  n <- nrow(X); d <- ncol(X)

  # Initialize posterior distribution using Gaussian hierarchical clustering
  # as proposed by Cule et al. (see their package LogConcDEAD/R/EMmixlcd.R)
  if (missing(posterior)) {
    highclust <- mclust::hc(modelName="VVV", X)
    class <- c(mclust::hclass(highclust, K))
    props <- rep(0, K)
    y <- matrix(0, nrow=n, ncol=K)
    for(i in 1:K) {
      props[i] <- sum(class==i) / n
      ss <- X[class==i, ]
      y[, i] <- mvtnorm::dmvnorm(X, mean=apply(ss, 2, mean), sigma=var(ss), log=TRUE)
    }

    pif <- t(t(exp(y)) * props)
    posterior <- pif / apply(pif, 1, sum)
  }

  ## obtain faces of convex hull of X
  cvhParams <- calcCvxHullFaces(X)

  gamma <- 1000
  minLogLike <- -Inf
  logLike <- matrix(0, maxIter, 1)
  densEst <- matrix(0, n, K)
  params <- matrix(list(), 1, K)
  # optimization parameters
  intEps = 1e-4
  objEps = 1e-7
  offset = 1e-1

  # initialize both densities; after that use previous parameters for initialization
  for (j in 1:K) {
    sampleWeights <- posterior[, j] / sum(posterior[, j])
    r <- fmlcd(X, w = sampleWeights, verbose = verbose-1)
    densEst[, j] <- exp(apply(-r$a %*% t(X) - matrix(rep(r$b, n), length(r$b), n), 2, min))
    params[[j]] <- c(r$a, r$b)
  }
  tau <- apply(posterior, 2, sum) / sum(posterior)

  for (iter in 1:maxIter) {
    # evaluate the log likelihood p(X|\beta)
    logLike[iter] <- sum(log(densEst %*% t(t(tau))))
    if (verbose > 0) {
      print(sprintf("Iter %d: Log-Likelihood = %.4f", iter, logLike[iter]))
    }

    # Check for convergence
    if (iter > 5) {
      likeDelta = abs((logLike[iter - 2] - logLike[iter]) / logLike[iter])
      if (likeDelta < 1e-5) {
        break
      }
    }

    # E-Step: Update posterior probabilities p(z|X,\beta)
    for (j in 1:K) {
      posterior[ ,j] <- tau[j] * densEst[, j] / (densEst %*% t(t(tau)))
    }

    # M-Step Update density parameters and mixing property tau
    for (j in 1:K) {
      sampleWeights <- posterior[, j] / sum(posterior[, j])
      sampleWeights[sampleWeights < 1e-8 / n] <- 0
      r <- callNewtonBFGSLC(X, sampleWeights, params[[j]], matrix(0, 0, 0),
                    cvhParams, gamma, verbose - 1, intEps, objEps, offset)

      result <- correctIntegral(X, rep(0,2), r$a, r$b, cvhParams$cvh);
      params[[j]] <- c(result$a, result$b)
      densEst[ ,j] <- exp(result$logMLE)
    }
    tau <- apply(posterior, 2, sum) / sum(posterior)
  }

  paramsReturn <- matrix(list(), 1, K)
  for (j in 1:K) {
    nH <- length(params[[j]])/(d+1)
    paramsReturn[[j]]$a = matrix(params[[j]][1:(d*nH)],nH,d)
    paramsReturn[[j]]$b = params[[j]][(d*nH+1):((d+1)*nH)]
  }
  r <- list("params" = paramsReturn, "densEst" = densEst, "tau" = tau)
  return(r)
}
