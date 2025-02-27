#' @import ggplot2
#' @import scatterplot3d
#' @importFrom graphics legend
#' @import stats
#' @importFrom plotly plot_ly
#' @importFrom plotly layout
#' @importFrom plotly add_surface
#' @importFrom scales alpha 
#'
#' @title Visualization of Fitted SpMix Model
#'
#' @description 
#' `plot.SpMix` visualizes the fitted semiparametric mixture model for objects of class "SpMix". 
#' It provides:
#' - A **histogram** with fitted densities for univariate data (1D)
#' - A **scatter plot** for 2D or 3D data, customizable with the `ggplot2` package.
#' - When `testing = TRUE` (default), it overlays null and alternative fitted densities for hypothesis testing.
#'
#' @param x An object of class "SpMix", obtained from `SpMix()`.
#' @param fdr_cutoff Threshold of local false discovery rate (localFDR) for highlighting significant data points (default: 0.2).
#' @param testing Logical; if `TRUE`, plots for hypothesis testing with null and alternative distributions. If `FALSE`, plots density estimation (default: `TRUE`).
#' @param xlab Label for the x-axis (default: "x").
#' @param ylab Label for the y-axis in 2D/3D plots (default: "y").
#' @param zlab Label for the z-axis in 3D plots (default: "z").
#' @param coord_legend Coordinates for the legend in a 3D scatter plot when data is 2D or 3D (default: `c(8, -5, 0.2)`).
#'
#' @return A visualization of the fitted `SpMix` model. Additionally, it returns:
#' 
#'   - `thre`: The threshold z-value for distinguishing between null and alternative distributions.
#'
#' @export
plot.SpMix <- function(x, fdr_cutoff = 0.2, testing = TRUE,
                      xlab = "x", ylab = "y", zlab = "z", 
                      coord_legend = c(8, -5, 0.2))
{
  if (!inherits(x, "SpMix")) stop("Input must be a SpMix object")
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
    which_z <- (localFDR <= fdr_cutoff)
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
      paste("Density Estimates: ", p[0], " = ", p0, ", ",
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
    sub_3d <- substitute(
      paste(italic(p)[0], " = ", p0, ", ",
        italic(mu)[0], " = (", mu01, ", ", mu02, ")"
      ),
      list(
        p0 = round(p0, 2),
        mu01 = round(mu0[1], 2),
        mu02 = round(mu0[2], 2)
      )
    )
    
    plots <- list()
    ngrid <- 50
    x1 <- seq(min(z[,1]), max(z[,1]), length = ngrid) 
    x2 <- seq(min(z[,2]), max(z[,2]), length = ngrid)
    grid_points <- as.matrix(expand.grid(x1, x2))
    
    comp0 <- x$p0 * dmvnorm(grid_points, x$mu0, x$sig0)
    comp0 <- matrix(comp0, ngrid, ngrid)
    comp1 <- (1 - x$p0) * dlcd(grid_points, x$lcd, uselog = FALSE)
    comp1 <- matrix(comp1, ngrid, ngrid)
    density <- den <- comp0 + comp1
    filled_contour_colormap <- hcl.colors(20, "YlOrRd", rev = TRUE)
    
    plots$contour3d <- function(){
      p <- plot_ly(x = x1, y = x2, z = ~density) %>% add_surface(
        colorscale = list(
          seq(0, 1, length.out = length(filled_contour_colormap)),
          filled_contour_colormap
        ),
        contours = list(
          z = list(
            show=TRUE,
            start = 0, end = max(density), size = max(density)/10,
            usecolormap=FALSE,
            highlightcolor="#ff0000",
            project=list(z=TRUE)
          )
        )
      ) %>% layout(
        title = list(
          text = "3D Contour Plot of Density Estimates",
          font = list(size = 18),   # Adjust font size
          x = 0.5,  # Center title horizontally (0 = left, 1 = right)
          y = 0.95, # Move title down (0 = bottom, 1 = top)
          xanchor = "center",
          yanchor = "top"
        ),
        scene = list(
          camera=list(
            eye = list(x=1.87, y=0.88, z=-0.64)
          )
        )
      )
      print(p)
    }

    if (!testing) {
      plots$density <- function(){
        filled.contour(x1, x2, den, plot.axes = {
          axis(1)
          axis(2)
          contour(x1, x2, den, add = TRUE, lwd = 1, col = "#999999")
          points(z, pch = 20, col = scales::alpha("black", 0.1))
        }, plot.title = title(main = "", xlab=xlab,  ylab=ylab)
        )
        title("Density Estimates", line = 3)
        title(sub_3d, line = 1.5)
      }
    } else{
      cols <- c("#999999", "#E69F00")[discovered]
      plots$density <- function(){
        filled.contour(x1, x2, den, plot.axes = {
          axis(1)
          axis(2)
          contour(x1, x2, comp0, add = TRUE, col = "#999999")
          contour(x1, x2, comp1, add = TRUE, col = "#E69F00")
          points(z, pch = 20, col = scales::alpha(cols, 0.2))
        }, plot.title = title(main = "", xlab=xlab,  ylab=ylab)
        )
        title("Density Estimates", line = 3)
        title(sub_3d, line = 1.5)
      }
      
      plots$localFDR <- function(){
        plot(x$localFDR, xlab = "Data Points", ylab = "local FDR",
             main = "Local FDR Estimates", 
             pch = 20, col = cols)
        abline(h = fdr_cutoff, col = "red", lty = 2)
      }
    }
    plots$density()
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
  return(plots)
}

