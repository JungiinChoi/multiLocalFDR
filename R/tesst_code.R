n = 1000

p0 <- 0.8
library(copula)
library(mvtnorm)
library(ggplot2)
library(scatterplot3d)

n0 <- rbinom(1, n, p0)
n1 <- n - n0
sig0 <- matrix(c(1, 0.25, 0.25, 1), ncol = 2, byrow = T)
V <- rCopula(n1, ellipCopula(family = "normal", param = 0.5, dim = 2))
z0 <- rmvnorm(n0, sigma = sig0)
z1 <- qnorm(V, mean = 2, sd = .5)
z_NN2d <- data.frame(rbind(z0, z1))


ggplot(z_NN2d, aes(x = z_NN2d[,1], y = z_NN2d[,2])) +
  geom_point(aes(color = factor(c(rep(1,n0), rep(2,n1)), 
                                labels = c("Normal(Null)", "Gaussian Copula(Alternative)"))), 
             size = 1) +
  scale_color_manual(values = c("#999999", "#E69F00"), name="")  +
  labs(title = "Random Generated Mixture density data", x = "x", y = "y") +
  theme_classic()

x <- SPMix(z_NN2d)

z <- x$z
p0 <- x$p0
mu0 <- x$mu0
sig0 <- x$sig0
f <- x$f
f1 <- x$f1
d <- x$dim

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

sub_3d <- sub_density

# Contour Plot
ngrid <- 50

x1 <- seq(from = min(z[,1]), to = max(z[,1]), length = ngrid)
x2 <- seq(from = min(z[,2]), to = max(z[,2]), length = ngrid)

comp0 <- x$p0 * dmvnorm(as.matrix(expand.grid(x1, x2)), x$mu0, x$sig0)
comp1 <- (1 - x$p0) * dlcd(as.matrix(expand.grid(x1, x2)), x$lcd, uselog = FALSE)
comp0 <- matrix(comp0, ngrid, ngrid)
comp1 <- matrix(comp1, ngrid, ngrid)
den <- comp0 + comp1

layout(matrix(1))

image(x1, x2, den, cex.axis = 0.7, xlab = "x", ylab = "y", 
      main = sub_3d, col = hcl.colors(50))

points(x$z, pch = 20, col = "gray")
contour(x1, x2, comp0, add = TRUE, col = "gray")
contour(x1, x2, comp1, add = TRUE, col = "gray")







z <- x$z
x$local.fdr <- x$posterior[,1]
thre_localFDR <- 0.2
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
      xlab = "", ylab = "", col = hcl.colors(50), main = sub_testing)
cols <- c("#999999", "#E69F00")[discovered]
points(z, pch = 20, col = alpha(cols, 0.7))
contour(x1, x2, comp0, add = TRUE, col = "gray")
contour(x1, x2, comp1, add = TRUE, col = "#E69F00")

plot(x$local.fdr, xlab = "index", ylab = "local fdr", 
     pch = 20, col = cols)
abline(h = thre_localFDR, col = 2, lty = 2)





z1_gamma <- qgamma(V, shape = 12, rate = 4)
z_NG <- data.frame(rbind(z0, z1_gamma))

ggplot(z_NG, aes(x = z_NG[,1], y = z_NG[,2])) +
  geom_point(aes(color = factor(c(rep(1,n0), rep(2,n1)), labels = c("Normal", "Gamma")))) +
  scale_color_manual(values = c("#999999", "#E69F00"), name="")  +
  labs(title = "Random Generated Normal/Gamma Mixture", x="x", y = "y") +
  theme_classic()

x <- SPMix(z_NG)



### Efron's Data

load("data/policez.Rdata")
hist(policez, breaks = 100)

load("data/Chisqz.Rdata")
hist(Chisqz, breaks = 100)

hist(hivz, breaks = 100)
?hist
