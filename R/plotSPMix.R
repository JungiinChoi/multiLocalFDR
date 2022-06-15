#' @import ggplot2
#' @import scatterplot3d
#' @importFrom graphics legend
#' @import stats
#'
#' @title plotMixture estimation for given z-values
#'
#' @description \code{plotSPMix} returns plotSPMix estimates and semiparametric
#' mixture density estimates for given multi-dimensional lists of z-values, which
#' are the probit-transformed p-values.
#' For the hypothesis testing \code{plotSPMix} uses a two-component semiparametric
#' mixture model to estimate the plotSPMix from the p-values. The two pillars of the
#' proposed approach are Efron's empirical null principle and log-concave density
#' estimation for the alternative distribution.
#'
#' @param z Matrix which column indicates given data point, it can be raw-data, z-values or p-values.
#' @param p0 Prior probability for null distribution
#' @param mu0 Parameter estimates of normal null distribution, N(mu0, sig0^2)
#' @param sig0 Parameter estimates of normal null distribution, N(mu0, sig0^2)
#' @param f Probability estimates of mixture model for each given data point.
#' @param f1 Probability estimates of alternative distribution of mixture model for each given data point.
#' @param localFDR localFDR estimates for given z-values
#' @param p_value If TRUE, the column of input indicates p-values, if FALSE, it indicates z-values or raw data. (default: FALSE)
#' @param alternative a character string specifying the alternative hypothesis, must be one of "two.sided", "greater" (default) or "less". You can specify just the initial letter. (default: "greater")
#' @param thre_localFDR Threshold of localFDR for null and alternative distribution (default: 0.2)
#' @param testing TRUE if it's for hypothesis testing, FALSEif it's for density estimation
#' @param xlab Label for x-axis on histogram or 3D scatter plot (default: "x")
#' @param ylab Label for y-axis on 3D scatter plot (default: "x")
#' @param type Type of 2-dimensional density plot (3d/contour plot)(default: "3d")
#' @param coord_legend Coordinate of a legend for 3d scatter plot when given data is 2D. (default: c(8, -5, 0.2))
#'
#' @return Plot estimated semiparametric mixture density and return threshold value.
#'
#'   \item{thre}{Threshold z-value for null and alternative distribution}
#'
#' @export
plotSPMix <- function(z, p0, mu0, sig0, f, f1, localFDR, p_value = FALSE,
                      alternative = "greater", thre_localFDR = 0.2, testing = TRUE,
                      xlab = "x", ylab = "y", type = "3d", coord_legend = c(8, -5, 0.2))
{
  which_z <- (localFDR <= thre_localFDR)

  z <- as.matrix(z)
  n <- dim(z)[1]
  d <- dim(z)[2]

  if (p_value) {
    if (alternative == "greater" | alternative == "g") {
      z = qnorm(1-z)
    } else {
      z = qnorm(z)
    }
  }

  if (d == 1){
    z = as.numeric(z)
    if (alternative == "greater" | alternative == "g"){
      thre <- min(z[which_z])
    } else{
      thre <- max(z[which_z])
    }

    legend_testing <- factor(which_z, levels = c(FALSE, TRUE), labels = c("Nonsignificant", "Significant"))
    legend_density <- factor((localFDR <= 0.5), levels = c(FALSE, TRUE), labels = c("Normal", "Nonparametric"))

    sub_testing=substitute(
      paste("Hypothesis Testing: ", p[0], " = ", p0, ", ",
            mu[0], " = ", mu0, ", ",
            sigma[0], " = ", sigma0, ", ",
            "threshold = ", threshold,
            sep = ""),
      list(p0 = round(p0, 2),
           mu0 = round(mu0, digits = 2),
           sigma0 = round(sig0, digits = 2),
           threshold = round(thre, digits = 2)))

    sub_density=substitute(
      paste("Density Estimation: ", p[0], " = ", p0, ", ",
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
        geom_histogram(aes(y = ..density..,),colour = 1, fill = "white", bins=100) +
        geom_line(aes(sort(z), (f[order(z)])),color = "gray25", lwd=1.1) +
        geom_line(aes(sort(z), p0*dnorm(zs, mean = mu0, sd = sig0)),color = "#0072B2",lwd=0.7) +
        geom_line(aes(sort(z), ((1-p0)*f1[order(z)])),color = "#D55E00",lwd=0.7) +
        geom_vline(aes(xintercept=thre), color="#E69F00",linetype="dashed") +
        geom_point(mapping = aes(x = thre, y = 0.01),size = 2,color="#E69F00",shape=25, fill="#E69F00") +
        labs(x=xlab, y = "density") +
        ggtitle(sub_testing) +
        theme(plot.title = element_text(margin = margin(b = -10))) +
        geom_rug(aes(z,color = legend_testing)) +
        scale_color_manual(values = c("#999999", "#E69F00"), name="") +
        theme_classic()
    } else {
      ggplot(df,aes(x=z)) +
        geom_histogram(aes(y = ..density..),colour = 1, fill = "white",bins=100) +
        geom_line(aes(sort(z), (f[order(z)])),color = "gray25", lwd=1.1) +
        geom_line(aes(sort(z), p0*dnorm(zs, mean = mu0, sd = sig0)),color = "#0072B2",lwd=0.7) +
        geom_line(aes(sort(z), ((1-p0)*f1[order(z)])),color = "#D55E00",lwd=0.7) +
        labs(x=xlab, y = "density") +
        ggtitle(sub_density) +
        theme(plot.title = element_text(margin = margin(b = -10))) +
        geom_rug(aes(z,color = legend_density)) +
        scale_color_manual(values = c("#0072B2", "#D55E00"), name="") +
        theme_classic()
    }

  } else if (d == 2) {
    sub_testing=substitute(
      paste("Hypothesis Testing: ", p[0], " = ", p0, ", ",
            mu[0], " = (", mu01,",",mu02, "), ",
            "localFDR threshold = ", threshold,
            sep = " "),
      list(p0 = round(p0, 2),
           mu01 = round(mu0[1], digits = 2),
           mu02 = round(mu0[2], digits = 2),
           threshold = round(thre_localFDR, digits = 2)))

    sub_density=substitute(
      paste("Density Estimation: ", p[0], " = ", p0, ", ",
            mu[0], " = (", mu01,",",mu02, "), ",
            sep = ""),
      list(p0 = round(p0, 2),
           mu01 = round(mu0[1], digits = 2),
           mu02 = round(mu0[2], digits = 2)))

    sub_3d <- if (testing) {sub_testing} else {sub_density}

    if (type == "contour" | type == "c") {
      # contour
      ggplot(df, aes(x = z[,1], y = z[,2])) +
        geom_point(aes(color = legend_3d)) +
        scale_color_manual(values = c("#999999", "#E69F00"), name="") +
        geom_density_2d(colour = "black", alpha= 0.7) +
        theme(plot.title = element_text(margin = margin(b = -10))) +
        labs(title = sub_3d, x=xlab, y = ylab) +
        theme_classic()
    } else {
      # 3D scatterplot
      colors <- c("#999999", "#E69F00")
      colors <- colors[as.numeric(which_z)+1]
      scatterplot<- scatterplot3d(z[,1],z[,2],f, pch = 16, color=colors,
                                  xlab = xlab, ylab = ylab, zlab = "f", main = sub_3d)
      legend_testing <- factor(which_z, levels = c(FALSE, TRUE), labels = c("Nonsignificant", "Significant"))
      legend_density <- factor((localFDR <= 0.5), levels = c(FALSE, TRUE), labels = c("Normal", "Nonparametric"))
      legend_3d <- if (testing) {legend_testing} else {legend_density}
      legend(scatterplot$xyz.convert(coord_legend[1], coord_legend[2], coord_legend[3]),
             legend = levels(legend_3d), col = c("#999999", "#E69F00"), pch = 16)
    }
  }
}
