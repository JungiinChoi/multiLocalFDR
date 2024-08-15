SPMix <- function(z, p0.init = 0.5, doplot = TRUE, thre = 0.5, tol = 5e-5, max_iter = 50) {
  
  # Internal Functions
  NormMix <- function(z, p0, tol = 5e-3, max_iter = 5) {
    require(mvtnorm)
    
    k <- 0
    converged <- FALSE
    z <- as.matrix(z)
    n <- nrow(z)
    d <- ncol(z)
    q <- qnorm(p0)
    z.scaled <- scale(z)
    
    ind0 <- rowSums(z.scaled <= q) == d
    ind1 <- rowSums(z.scaled > -q) == d
    
    z0 <- z[ind0, , drop = FALSE]
    z1 <- z[ind1, , drop = FALSE]
    
    mu0 <- colMeans(z0)
    sig0 <- cov(z0)
    f0 <- dmvnorm(z, mu0, sig0)
    
    mu1 <- colMeans(z1)
    sig1 <- cov(z1)
    f1 <- dmvnorm(z, mu1, sig1)
    
    f <- p0 * f0 + (1 - p0) * f1
    ell <- mean(log(f))
    
    while (k < max_iter && !converged) {
      k <- k + 1
      
      # E-step
      gam <- p0 * f0 / f
      
      # M-step
      new_p0 <- mean(gam)
      new_mu0 <- colSums(gam * z) / sum(gam)
      dev0 <- (z - new_mu0) * sqrt(gam)
      new_sig0 <- t(dev0) %*% dev0 / sum(gam)
      f0 <- dmvnorm(z, new_mu0, new_sig0)
      
      new_mu1 <- colSums((1 - gam) * z) / sum(1 - gam)
      dev1 <- (z - new_mu1) * sqrt(1 - gam)
      new_sig1 <- t(dev1) %*% dev1 / sum(1 - gam)
      f1 <- dmvnorm(z, new_mu1, new_sig1)
      
      f <- new_p0 * f0 + (1 - new_p0) * f1
      new_ell <- mean(log(f))
      
      # Update and check convergence
      converged <- abs(new_ell - ell) <= tol
      p0 <- new_p0
      mu0 <- new_mu0
      sig0 <- new_sig0
      mu1 <- new_mu1
      sig1 <- new_sig1
      ell <- new_ell
    }
    
    return(list(p0 = p0, mu0 = mu0, sig0 = sig0, mu1 = mu1, sig1 = sig1))
  }
  
  NE <- function(x, X) {
    n <- nrow(X)
    p <- ncol(X)
    xx <- matrix(x, nrow = n, ncol = p, byrow = TRUE)
    ne.ind <- apply(X >= xx, 1, all)
    which(ne.ind)
  }
  
  MonotoneFDR <- function(z, fdr) {
    n <- nrow(z)
    MFDR <- numeric(n)
    for (i in 1:n) {
      MFDR[i] <- max(fdr[NE(z[i, ], z)])
    }
    MFDR
  }
  
  # Main Function
  require(mvtnorm)
  require(LogConcDEAD)
  
  z <- as.matrix(z)
  n <- nrow(z)
  d <- ncol(z)
  ell <- numeric(max_iter)
  
  # Initial step: fit normal mixture
  Params <- NormMix(z, p0 = p0.init)
  p0 <- Params$p0
  mu0 <- Params$mu0
  sig0 <- Params$sig0
  f0 <- dmvnorm(z, mu0, sig0)
  f1 <- dmvnorm(z, Params$mu1, Params$sig1)
  f <- p0 * f0 + (1 - p0) * f1
  ell[1] <- mean(log(f), na.rm = TRUE)
  
  # EM-step
  k <- 1
  converged <- FALSE
  while (k < max_iter && !converged) {
    k <- k + 1
    
    # E-step
    gam <- p0 * f0 / f
    gam <- MonotoneFDR(z, gam)
    
    # M-step
    sum_gam <- sum(gam)
    mu0 <- colSums(gam * z) / sum_gam
    dev <- (z - mu0) * sqrt(gam)
    sig0 <- t(dev) %*% dev / sum_gam
    p0 <- mean(gam)
    f0 <- dmvnorm(z, mu0, sig0)
    
    which_z <- gam <= thre
    weight <- 1 - gam[which_z]
    weight <- weight / sum(weight)
    
    lcd <- fmlcd(z[which_z, ], w = weight)
    f1 <- rep(0, n)
    f1[which_z] <- exp(lcd$logMLE)
    f <- p0 * f0 + (1 - p0) * f1
    ell[k] <- mean(log(f), na.rm = TRUE)
    
    # Update and check convergence
    converged <- abs(ell[k] - ell[k - 1]) <= tol
    cat(".")
    if (converged) {
      cat("Converged!\n")
    }
  }
  
  if (!converged) {
    cat("Warning: Not converged!\n")
  }
  
  # Return results
  list(
    z = z, p0 = p0, mu0 = mu0, sig0 = sig0, f1 = f1, f = f,
    iter = k, log.likelihood = ell, lcd = lcd, dim = d,
    posterior = cbind(gam, 1 - gam), converged = converged
  )
}

fit.sp.mix.mult <- function(z, np.left = NULL, doplot = TRUE, p0.init = 0.5, thre = 0.5, tol = 5e-4, max.iter = 50) {
  d <- ncol(z)
  zx <- z
  
  if (!is.null(np.left)) {
    for (j in 1:d) {
      if (np.left[j] == "Yes") {
        zx[, j] <- -z[, j]
      }
    }
  }
  
  fit <- SPMix(zx, tol = tol, max_iter = max.iter, p0.init = p0.init, thre = thre)
  
  if (!is.null(np.left)) {
    for (j in 1:d) {
      if (np.left[j] == "Yes") {
        fit$mu0[j] <- -fit$mu0[j]
      }
    }
  }
  
  if (doplot) {
    if (d == 2) {
      ngrid <- 50
      x1 <- seq(min(z[, 1]), max(z[, 1]), length = ngrid)
      x2 <- seq(min(z[, 2]), max(z[, 2]), length = ngrid)
      grid <- expand.grid(x1, x2)
      comp0 <- fit$p0 * dmvnorm(grid, fit$mu0, fit$sig0)
      comp0 <- matrix(comp0, ngrid, ngrid)
      comp1 <- (1 - fit$p0) * dlcd(grid, fit$lcd, uselog = FALSE)
      comp1 <- matrix(comp1, ngrid, ngrid)
      den <- comp0 + comp1
      
      par(mfrow = c(1, 1))
      image(x1, x2, den, cex.axis = 0.7, xlab = "", ylab = "", col = gray.colors(128))
      points(z, pch = 20)
      contour(x1, x2, comp0, add = TRUE)
      contour(x1, x2, comp1, add = TRUE)
    } else if (d == 3) {
      require(scatter3d)
      scatterplot3d(z, xlab = "", ylab = "", zlab = "")
    } else {
      cat("Warning: d > 3 and doplot = TRUE is ignored.")
    }
  }
  
  return(fit)
}

local.fdr.mult <- function(z, tol = 5e-4, max.iter = 50, doplot = TRUE, p0.init = 0.9, thre = 0.2) {
  fit <- SPMix(z, tol = tol, max_iter = max.iter, p0.init = p0.init, thre = 1 - thre)
  fit$local.fdr <- fit$posterior[, 1]
  fit$discovered <- fit$local.fdr <= thre
  
  if (doplot) {
    if (ncol(z) == 2) {
      ngrid <- 50
      x1 <- seq(min(z[, 1]), max(z[, 1]), length = ngrid)
      x2 <- seq(min(z[, 2]), max(z[, 2]), length = ngrid)
      grid <- expand.grid(x1, x2)
      comp0 <- fit$p0 * dmvnorm(grid, fit$mu0, fit$sig0)
      comp0 <- matrix(comp0, ngrid, ngrid)
      comp1 <- (1 - fit$p0) * dlcd(grid, fit$lcd, uselog = FALSE)
      comp1 <- matrix(comp1, ngrid, ngrid)
      den <- comp0 + comp1
      
      par(mfrow = c(1, 2))
      image(x1, x2, den, cex.axis = 0.7, xlab = "", ylab = "", col = gray.colors(128))
      points(z, pch = 20, col = 4 - 2 * fit$discovered)
      contour(x1, x2, comp0, add = TRUE)
      contour(x1, x2, comp1, add = TRUE)
      
      plot(fit$local.fdr, xlab = "index", ylab = "local fdr", pch = 20, col = 4 - 2 * fit$discovered)
      abline(h = thre, col = 2, lty = 2)
    } else if (ncol(z) == 3) {
      require(scatterplot3d)
      scatterplot3d(z, xlab = "", ylab = "", zlab = "")
    } else {
      plot(fit$local.fdr, xlab = "index", ylab = "local fdr", col = 4 - 2 * fit$discovered)
      abline(h = thre, col = 2, lty = 2)
    }
  }
  
  return(fit)
}

##########
# Examples #
##########

N <- 1000
p0 <- 0.3
N0 <- rbinom(1, size = N, prob = p0)
membership <- c(rep(0, N0), rep(1, N - N0))

library(mvtnorm)
d <- 2
mu0 <- rep(1, d)
sig0 <- diag(d)
mu1 <- c(3, 3)
sig1 <- 0.5 * sig0
z <- rbind(rmvnorm(N0, mu0, sig0), 
           rmvnorm(N - N0, mu1, sig1))
res <- fit.sp.mix.mult(z, np.left = c("No", "No"))



N <- 1000
p0 <- 0.9
N0 <- round(N*p0) #rbinom(1, size = N, prob = p0)
membership <- c(rep(0, N0), rep(1, N - N0))
d <- 2
mu0 <- rep(0, d)
sig0 <- diag(d)
mu1 <- rep(3, d)
sig1 <- 0.5 * sig0
z <- rbind(rmvnorm(N0, mu0, sig0), 
           rmvnorm(N - N0, mu1, sig1))
res <- local.fdr.mult(z)


### Simulations
library(scatterplot3d)
library(mvtnorm)
d <- 2
mu0 <- rep(1, d)
sig0 <- diag(d)
mu1 <- c(3, 3)
sig1 <- 0.5 * sig0
z <- rbind(rmvnorm(N0, mu0, sig0), 
           rmvnorm(N - N0, mu1, sig1))
res <- fit.sp.mix.mult(z, np.left = c("No", "No"))

par(mfrow=c(1,1))
z <- res$z
which_z <- res$local.fdr <= 0.5
f <- res$f
colors <- c("#999999", "#E69F00")
colors <- colors[as.numeric(which_z)+1]
scatterplot<- scatterplot3d(z[,1],z[,2],f, pch = 16, color=colors,
                            xlab = "", ylab = "", zlab = "f")
length(colors)

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

library(copula)
library(mvtnorm)
library(scatterplot3d)
p0 <- 0.8
n <- 1000
n0 <- rbinom(1, n, p0)
n1 <- n - n0
sig0 <- matrix(c(1, 0.25, 0.25, 1), ncol = 2, byrow = T)
V <- rCopula(n1, ellipCopula(family = "normal", param = 0.5, dim = 2))
z0 <- rmvnorm(n0, sigma = sig0)
z1 <- qnorm(V, mean = 2.0, sd = .5)
z_NN2d <- data.frame(rbind(z0, z1))
res <- fit.sp.mix.mult(z_NN2d, np.left = c("No", "No"))


fit <- res
thre <- 0.5
z <- fit$z

ngrid <- 50
x1 <- seq(from = min(z[,1]), to = max(z[,1]), length = ngrid)
x2 <- seq(from = min(z[,2]), to = max(z[,2]), length = ngrid)
par(mfrow = c(1, 1))
comp0 <- fit$p0*dmvnorm(as.matrix(expand.grid(x1, x2)), fit$mu0, fit$sig0)
comp0 <- matrix(comp0, ngrid, ngrid)
comp1 <- (1 - fit$p0)*(dlcd(as.matrix(expand.grid(x1, x2)), fit$lcd, uselog = FALSE))
comp1 <- matrix(comp0, ngrid, ngrid)
den <- comp0 + comp1
image(x1, x2, den, cex.axis = 0.7, 
      xlab = "", ylab = "", col = gray.colors(128))
points(z, pch = 20)
contour(x1, x2, comp0, add = TRUE)
contour(x1, x2, comp1, add = TRUE)

# Set the number of grid points
ngrid <- 50

# Generate x1 and x2 sequences
x1 <- seq(from = min(z[,1]), to = max(z[,1]), length = ngrid)
x2 <- seq(from = min(z[,2]), to = max(z[,2]), length = ngrid)

# Compute component densities
comp0 <- fit$p0 * dmvnorm(as.matrix(expand.grid(x1, x2)), fit$mu0, fit$sig0)
comp1 <- (1 - fit$p0) * dlcd(as.matrix(expand.grid(x1, x2)), fit$lcd, uselog = FALSE)

# Reshape component densities matrices
comp0 <- matrix(comp0, ngrid, ngrid)
comp1 <- matrix(comp1, ngrid, ngrid)

# Compute density
den <- comp0 + comp1

# Set the layout
layout(matrix(1))

# Plot the image
image(x1, x2, den, cex.axis = 0.7, xlab = "X1", ylab = "X2", col = hcl.colors(50))

# Add points
points(z, pch = 20, col = "gray")

# Add contour lines for comp0
contour(x1, x2, comp0, add = TRUE, col = "white")

# Add contour lines for comp1
contour(x1, x2, comp1, add = TRUE, col = "white")


res <- local.fdr.mult(z_NN2d, p0.init = 0.7)


library(ggplot2)
library(MASS)

# Assuming 'z' is your data frame with columns x and y

# Define the number of grid points
ngrid <- 50

# Generate grid points
x1 <- seq(min(z[,1]), max(z[,1]), length.out = ngrid)
x2 <- seq(min(z[,2]), max(z[,2]), length.out = ngrid)
grid <- expand.grid(x1 = x1, x2 = x2)

# Calculate densities for components
comp0 <- fit$p0 * dmvnorm(grid, fit$mu0, fit$sig0)
comp1 <- (1 - fit$p0) * dlcd(as.matrix(grid), fit$lcd, uselog = FALSE)
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

library(plotly)
library(reshape2)

pointdf <- data.frame(x = z[,1], y = z[,2])

p <- ggplot(contour_data, aes(x1, x2, z= den)) +
  stat_contour(geom = "polygon", aes(fill=stat(level))) +
  scale_fill_distiller(palette = "Spectral", direction = -1)

for (i in 1:nrow(z)){
  p <- p + geom_point(aes(z[i,1],z[i,2]))
}

p + geom_point(pointdf, aes(x,y))

x <- res
z <- x$z
thre_localFDR <- 0.2
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
