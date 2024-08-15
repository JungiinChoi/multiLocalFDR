#' @import ggplot2
#' @import scatterplot3d
#' @importFrom graphics legend
#' @import stats
#'
#' @title draws a histogram(univariate) or a scatterplot(2D or 3D) for fitted SpMix object.
#' 
#' @description \code{plotSPMix} draws a histogram(univariate) or a scatterplot(2D or 3D) 
#' for fitted SpMix object which can be customized with the ggplot2 package.
#' \code{plotSpMix()} can be used for both hypothesis testing and density estimation. 
#' If testing is TRUE (default), it gives fitted density for both null and alternative distribution. 
#'
#'
#' @param x object of class "SpMix"
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
plotSpMix <- function(x, thre_localFDR = 0.2, testing = FALSE,
                      xlab = "x", ylab = "y", zlab = "z", 
                      coord_legend = c(8, -5, 0.2))
{
  z <- x$z
  p0 <- x$p0
  mu0 <- x$mu0
  sig0 <- x$sig0
  f <- x$f
  f1 <- x$f1
  d <- x$dim
  localFDR <- x$localFDR
  greater_alt <- x$greater_alt

  if (d == 1){
    z = as.numeric(z)
    which_z <- (localFDR <= thre_localFDR)
    if (is.null(greater_alt)){
      thre <- min(z[which_z])
    } else{
      thre <- max(z[which_z])
    }

    legend_testing <- factor(which_z, levels = c(FALSE, TRUE), labels = c("Null", "Non-null"))
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
    
    if (!testing) {
      # Contour Plot
      ngrid <- 50
      
      x1 <- seq(from = min(z[,1]), to = max(z[,1]), length = ngrid)
      x2 <- seq(from = min(z[,2]), to = max(z[,2]), length = ngrid)
      
      comp0 <- x$p0 * dmvnorm(as.matrix(expand.grid(x1, x2)), x$mu0, x$sig0)
      comp1 <- (1 - x$p0) * dlcd(as.matrix(expand.grid(x1, x2)), x$lcd, uselog = FALSE)
      comp0 <- matrix(comp0, ngrid, ngrid)
      comp1 <- matrix(comp1, ngrid, ngrid)
      den <- comp0 + comp1
      
      # Combine data for ggplot
      contour_data <- data.frame(grid, den)
      
      # Create the plot
      ggplot(contour_data) +
        geom_contour(aes(x = x1, y = x2, z = den)) +
        stat_contour(geom = "polygon", aes(fill=stat(level))) +
        scale_fill_distiller(palette = "Spectral", direction = -1)+
        geom_point(data = pointdf, aes(x = x, y = y), shape = 20) +
        labs(title = "Fancier Contour Plot", x = "X-axis", y = "Y-axis") +
        theme_minimal()
      
      layout(matrix(1))
      
      image(x1, x2, den, cex.axis = 0.7, xlab = xlab, ylab = ylab, 
            main = sub_3d, col = hcl.colors(50))
      points(x$z, pch = 20, col = "gray")
      contour(x1, x2, comp0, add = TRUE, col = "white")
      contour(x1, x2, comp1, add = TRUE, col = "white")
    } else{
      z <- x$z
      x$local.fdr <- x$posterior[,1]
      discovered <- (x$local.fdr <= thre_localFDR)  + 1
      library(scales)
      
      ngrid <- 50
      x1 <- seq(from = min(z[,1]), to = max(z[,1]), length = ngrid) 
      x2 <- seq(from = min(z[,2]), to = max(z[,2]), length = ngrid)
      par(mfrow = c(1, 2))
      comp0 <- x$p0*dmvnorm(as.matrix(expand.grid(x1, x2)), x$mu0, x$sig0)
      comp0 <- matrix(comp0, ngrid, ngrid)
      comp1 <- (1 - x$p0)*dlcd(as.matrix(expand.grid(x1, x2)), x$lcd, uselog = FALSE)
      comp1 <- matrix(comp1, ngrid, ngrid)
      den <- comp0 + comp1
      image(x1, x2, den, cex.axis = 0.7, 
            xlab = "", ylab = "", col = hcl.colors(50))
      cols <- c("#999999", "#E69F00")[discovered]
      points(z, pch = 20, col = alpha(cols, 0.4))
      contour(x1, x2, comp0, add = TRUE)
      contour(x1, x2, comp1, add = TRUE)
      
      plot(x$local.fdr, xlab = "index", ylab = "local fdr", 
           pch = 20, col = cols)
      abline(h = thre_localFDR, col = 2, lty = 2)
    }
  
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
