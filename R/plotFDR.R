#' @importFrom fmlogcondens fmlcd
#' @importFrom mclust dmvnorm
#' @importFrom logcondens activeSetLogCon
#'
#' @title plotFDR estimation for given z-values
#'
#' @description \code{plotFDR} returns plotFDR estimates and semiparametric
#' mixture density estimates for given multi-dimensional lists of z-values, which
#' are the probit-transformed p-values.
#' For the hypothesis testing \code{plotFDR} uses a two-component semiparametric
#' mixture model to estimate the plotFDR from the p-values. The two pillars of the
#' proposed approach are Efron's empirical null principle and log-concave density
#' estimation for the alternative distribution.
#'
#' @param z Matrix which column indicates z-values, probit-transformed p-values.
#' @param p0 Prior probability for null distribution
#' @param mu0 sig0 Parameter estimates of normal null distribution, N(mu0, sig0^2)
#' @param f1 Probability estimates of alternative distribution of mixture model for each z-value point.
#' @param localFDR localFDR estimates for given z-values
#' @param alternative a character string specifying the alternative hypothesis, must be one of "two.sided", "greater" (default) or "less". You can specify just the initial letter. (default: "greater")
#' @param thre_localFDR Threshold of localFDR for null and alternative distribution (default: 0.2)
#' @param density TRUE if it's for density estimation, FALSE if it's for hypothesis testing
#'
#' @return Plot estimated semiparametric mixture density and return threshold value.
#'
#'   \item{thre}{Threshold z-value for null and alternative distribution}
#'
#' @export
plotFDR <- function(z, p0, mu0, sig0, f1, localFDR, alternative = "greater", thre_localFDR = 0.2, density = FALSE)
  # FOR MULTIVARIATE CASE ONLY
{
  which_z <- localFDR <= thre_localFDR

  if (alternative == "greater" | alternative == "g"){
    thre <- min(z[which_z])
  }
  else{
    thre <- max(z[which_z])
  }

  hist(z,
       nclass = max(round(length(z)/20), 24),
       probability = TRUE,
       col = "gray", border = "white",
       xlab = "",
       main = "",
       sub = substitute(
         paste(p[0], " = ", p0, ", ",
               mu[0], " = ", mu0, ", ",
               sigma[0], " = ", sigma0, ", ",
               "threshold = ", threshold,
               sep = ""),
         list(p0 = round(p0, 2),
              mu0 = round(mu0, 2),
              sigma0 = round(sig0, 2),
              threshold = round(thre, 2))))
  rug(z, col = "gray")
  rug(z[which_z], col = 2)
  z_sorted <- sort(z)
  lines(z_sorted, p0 * dnorm(z_sorted, mu0, sig0), col = 3, lwd = 2)
  lines(z_sorted, (1-p0) * f1[order(z)], col = 2, lwd = 2)
  points(thre, 0, bg = "yellow", col = 2, pch = 25)

  return(thre = thre)

}

