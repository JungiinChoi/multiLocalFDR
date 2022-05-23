#' @import ggplot2
#'
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
#' @param p_value If TRUE, the column of input indicates p-values, if FALSE, it indicates z-values or raw data. (default: FALSE)
#' @param alternative a character string specifying the alternative hypothesis, must be one of "two.sided", "greater" (default) or "less". You can specify just the initial letter. (default: "greater")
#' @param thre_localFDR Threshold of localFDR for null and alternative distribution (default: 0.2)
#' @param testing TRUE if it's for hypothesis testing, FALSEif it's for density estimation
#'
#' @return Plot estimated semiparametric mixture density and return threshold value.
#'
#'   \item{thre}{Threshold z-value for null and alternative distribution}
#'
#' @export
plotFDR <- function(z, p0, mu0, sig0, f1, localFDR, p_value = FALSE,
                    alternative = "greater", thre_localFDR = 0.2, testing = TRUE)
{
  which_z <- (localFDR <= thre_localFDR)

  if (p_value) {
    if (alternative == "greater" | alternative == "g") {
      z = qnorm(1-z)
    } else {
      z = qnorm(z)
    }
  }

  z <- as.matrix(z)
  n <- dim(z)[1]
  d <- dim(z)[2]

  if (d == 1){
    z = as.numeric(z)
    if (alternative == "greater" | alternative == "g"){
      thre <- min(z[which_z])
    } else{
      thre <- max(z[which_z])
    }

    legend_testing <- factor(which_z, levels = c(FALSE, TRUE), labels = c("Null", "Alternative"))
    legend_density <- factor(which_z, levels = c(FALSE, TRUE), labels = c("Normal", "Nonparametric"))

    sub_testing=substitute(
      paste(p[0], " = ", p0, ", ",
            mu[0], " = ", mu0, ", ",
            sigma[0], " = ", sigma0, ", ",
            "threshold = ", threshold,
            sep = ""),
      list(p0 = round(p0, 2),
           mu0 = round(mu0, digits = 2),
           sigma0 = round(sig0, digits = 2),
           threshold = round(thre, digits = 2)))

    sub_density=substitute(
      paste(p[0], " = ", p0, ", ",
            mu[0], " = ", mu0, ", ",
            sigma[0], " = ", sigma0,
            sep = ""),
      list(p0 = round(p0, 2),
           mu0 = round(mu0, digits = 2),
           sigma0 = round(sig0, digits = 2)))

    df = data.frame(z=z)
    zs <- sort(z)
    if (testing) {
      ggplot(df,aes(x=z)) +
        geom_histogram(aes(y = ..density..),colour = 1, fill = "white",bins=100) +
        geom_line(aes(sort(z), p0*dnorm(zs, mean = mu0, sd = sig0)),color = "#0072B2",lwd=1.1) +
        geom_line(aes(sort(z), ((1-p0)*f1[order(z)])),color = "#D55E00",lwd=1.1) +
        geom_vline(aes(xintercept=mu0), color="#0072B2",linetype="dashed") +
        geom_point(mapping = aes(x = thre, y = 0.01),size = 2,color='#F0E442',shape=25,fill="#F0E442") +
        labs(x="z-value", y = "density") +
        ggtitle(sub_testing) +
        theme(plot.title = element_text(margin = margin(b = -10))) +
        geom_rug(aes(z,color = legend_testing))+
        scale_color_manual(values = c("#0072B2", "#D55E00"), name="")+
        theme_classic()
    } else {
      ggplot(df,aes(x=z)) +
        geom_histogram(aes(y = ..density..),colour = 1, fill = "white",bins=100) +
        geom_line(aes(sort(z), p0*dnorm(zs, mean = mu0, sd = sig0)),color = "#0072B2",lwd=1.1) +
        geom_line(aes(sort(z), ((1-p0)*f1[order(z)])),color = "#D55E00",lwd=1.1) +
        geom_vline(aes(xintercept=mu0), color="#0072B2",linetype="dashed") +
        labs(x="x", y = "density") +
        ggtitle(sub_density) +
        theme(plot.title = element_text(margin = margin(b = -10))) +
        geom_rug(aes(z,color = legend_density))+
        scale_color_manual(values = c("#0072B2", "#D55E00"), name="")+
        theme_classic()
    }

  }

}

