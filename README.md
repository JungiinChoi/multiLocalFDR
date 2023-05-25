
# multiLocalFDR 


## Overview

multiLocalFDR is a package for multi-dimensional local-FDR estimation using a semiparametric mixture method.
The two pillars of the proposed approach are Efron's empirical null principle and log-concave density estimation for the alternative distribution. A unique feature of our method is that it can be extended to compute the local false discovery rates by combining multiple lists of p-values.

multiLocalFDR package is built around its main functions SPMix() and plotSPMix():

  - `SPMix()` fits the semiparametric mixture model and computes local-FDR, density estimates for given data. 
  - `plotSPMix()` draws a histogram(univariate) or a scatterplot(2D or 3D) for fitted SPMix object which can be customized with the ggplot2 package.
  
It also provides other useful functions for local-FDR estimation:
  - `localFDR()` provides localFDR estimates for given multi-dimensional lists of raw data.

There are three external data contained in multiLocalFDR package:
  - `pathways`: The list of p-values of significantly upregulated genes in three published studies
  - `galaxy`: Radial velocities of globular clusters of M104 and Milky Way stars.
  - `microarrays`: List of p-values of microarray data for three different experiments.

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
#> package 'multiLocalFDR' overwrites 'fmlogcondens' to modified version in https://github.com/JungiinChoi/fmlogcondens.
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
SPMix_galaxy <- SPMix(z_galaxy, alternative = "less")

# plot fitted mixture density
plotSPMix(SPMix_galaxy, testing = FALSE, xlab = "velocity(km/s)")
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

# Fit 2-dimensional mixture density and local-FDR values
NN2d <- SPMix(z_NN2d)

# Plot 3-dimensional scatter plot for fitted density
plotSPMix(NN2d)

```

## Getting help

If you encounter a clear bug, please file an issue with a minimal
reproducible example on
[GitHub](https://github.com/JungiinChoi/multiLocalFDR/issues). For questions and
other discussion, feel free to contact me: Jungin Choi (jchoi177 at jhu.edu).

