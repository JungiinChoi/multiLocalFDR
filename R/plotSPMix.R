#' @import ggplot2
#' @import scatterplot3d
#' @importFrom graphics legend
#' @import stats
#'
#' @title draws a histogram(univariate) or a scatterplot(2D or 3D) for fitted SPMix object.
#' 
#' @description \code{plotSPMix} draws a histogram(univariate) or a scatterplot(2D or 3D) 
#' for fitted SPMix object which can be customized with the ggplot2 package.
#' \code{plotSPMix()} can be used for both hypothesis testing and density estimation. 
#' If testing is TRUE (default), it gives fitted density for both null and alternative distribution. 
#'
#'
#' @param x object of class "SPMix"
#' @param thre_localFDR Threshold of localFDR determining significant data points. (default: 0.2)
#' @param testing If TRUE, it's for hypothesis testing. If FALSE, it's for density estimation (default: TRUE)
#' @param xlab Label for x-axis on histogram or 3D scatter plot (default: "x")
#' @param ylab Label for y-axis on 3D scatter plot (default: "y")
#' @param zlab Label for z-axis on 3D scatter plot (default: "z")
#' @param coord_legend Coordinate of a legend for a 3d scatter plot when given data is 2D or 3D. (default: c(8, -5, 0.2))
#'
#' @return Plot estimated semiparametric mixture density and return threshold value.
#'
#'   \item{thre}{Threshold z-value for null and alternative distribution}
#'
#' @export
plotSPMix <- function(x, thre_localFDR = 0.2, testing = TRUE,
                      xlab = "x", ylab = "y", zlab = "z", coord_legend = c(8, -5, 0.2))
{
  z <- x$z
  p0 <- x$p0
  mu0 <- x$mu0
  sig0 <- x$sig0
  f <- x$f
  f1 <- x$f1
  localFDR <- x$localFDR
  d <- x$dim
  alternative <- x$alternative
  
  which_z <- (localFDR <= thre_localFDR)
  

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
  
  } else if (d == 3) {
    colors <- c("#999999", "#E69F00")
    colors <- colors[as.numeric(which_z)+1]
    scatterplot<- scatterplot3d(z[,1],z[,2],z[,3], pch = 16, color=colors,
                                xlab = xlab, ylab = ylab, zlab = zlab)
    legend_testing <- factor(which_z, levels = c(FALSE, TRUE), labels = c("Nonsignificant", "Significant"))
    legend_density <- factor((localFDR <= 0.5), levels = c(FALSE, TRUE), labels = c("Normal", "Nonparametric"))
    legend_3d <- if (testing) {legend_testing} else {legend_density}
    legend(scatterplot$xyz.convert(coord_legend[1], coord_legend[2], coord_legend[3]),
           legend = levels(legend_3d), col = c("#999999", "#E69F00"), pch = 16)
  }
}
