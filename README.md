
# multiLocalFDR 


## Overview

multiLocalFDR is a package for multi-dimensional local-FDR estimation using a semiparametric mixture method.
The two pillars of the proposed approach are Efron's empirical null principle and log-concave density estimation for the alternative distribution. A unique feature of our method is that it can be extended to compute the local false discovery rates by combining multiple lists of p-values.

  - `SPMix()` provides estimates parameters of null and alternative distribution of our semiparametric mixture method.
  - `localFDR()` provides estimates of local-FDR for given lists of z-values / p-values.
  - `FDR()` provides estimates of FDR for given lists of z-values / p-values.
  - `plotSPMix()` plots estimated semiparametric mixture distribution and provides threshold z-value for null and alternative distribution.

You can learn more about them in `vignette("multiLocalFDR")`. 

## Installation

You can install multiLocalFDR from GitHub.

``` r
# install.packages("devtools")
devtools::install_github("JungiinChoi/multiLocalFDR")
```

### fmlogcondens version

multiLocalFDR imports the modified version of fmlogcondens from my [GitHub](https://github.com/JungiinChoi/fmlogcondens).

If you already have the [original fmlogcondens](https://github.com/FabianRathke/fmlogcondens), multiLocalFDR will overwrite this package and give a warning. 

``` r
# install.packages("devtools")
devtools::install_github("JungiinChoi/multiLocalFDR")

#> Warning message:
#> package 'multiLocalFDR' overwrites 'fmlogcondens' to modified version.
```

## Usage

### 1-dimensional 

``` r
library(multiLocalFDR)

# Included Data: Radial velocities of globular clusters
data("galaxy", package = "multiLocalFDR")
head(galaxy)

# Estimate mixture density with greater values of normal distribution. 
# (alternative = "less")
z_galaxy <- galaxy$velocity
params <- SPMix(z_galaxy, alternative = "less")

# plot fitted mixture density
plotSPMix(z_galaxy, params$p0, params$mu0, params$sig0, params$f, params$f1, 
          params$localFDR, alternative = "greater", testing = FALSE,
                      xlab = "velocity(km/s)")
```

### multi-dimensional

``` r
library(multiLocalFDR)

# Randomly generated data from Gaussian mixture density with null probability 0.8
n <- 1000; p0 <- 0.8
n0 <- rbinom(1, n, p0)
n1 <- n - n0
sig0 <- matrix(c(1, 0.25, 0.25, 1), ncol = 2, byrow = T)
V <- rCopula(n1, ellipCopula(family = "normal", param = 0.5, dim = 2))
z0 <- rmvnorm(n0, sigma = sig0)
z1 <- qnorm(V, mean = 3.5, sd = .5)
z_NN2d <- data.frame(rbind(z0, z1))

# Plot generated data
ggplot(z_NN2d, aes(x = z_NN2d[,1], y = z_NN2d[,2])) +
  geom_point(aes(color = factor(c(rep(1,n0), rep(2,n1)), labels = c("Normal(Null)", "Gaussian Copula(Alternative)")))) +
  scale_color_manual(values = c("#999999", "#E69F00"), name="")  +
  labs(title = "Random Generated Mixture density data", x="x", y = "y") +
  theme_classic()

# Fit 2-dimensional mixture density and local-FDR values
NN2d <- SPMix(z_NN2d, tol = 1e-10)

# Plot 3-dimensional scatter plot for fitted density
plotSPMix(z_NN2d, NN2d$p0, NN2d$mu0, NN2d$sig0, NN2d$f, NN2d$f1, NN2d$localFDR, 
          xlab = "x", ylab = "y")
```

## Getting help

If you encounter a clear bug, please file an issue with a minimal
reproducible example on
[GitHub](https://github.com/JungiinChoi/multiLocalFDR/issues). For questions and
other discussion, feel free to contact me: Jungin Choi (serimtech07 at snu.ac.kr).

