
# multiLocalFDR 


## Overview

multiLocalFDR is a package for multi-dimensional local-FDR estimation using a semiparametric mixture method.
The two pillars of the proposed approach are Efron's empirical null principle and log-concave density estimation for the alternative distribution. A unique feature of our method is that it can be extended to compute the local false discovery rates by combining multiple lists of p-values.

  - `localFDR()` provides estimates of local-fdr for given lists of p-values.
  - `SpMixParam()` provides estimates parameters of null and alternative distribution of our method.
  - `nMixParam()` provides estimates parameters of null and alternative distribution of normal mixture model. 
  - `arrangeNE()` arranges given data as an increasing order for multi-dimensional data.

You can learn more about them in
`vignette("multiLocalFDR")`. 

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

# z-values (probit transformed p-values) from data points

z<-Carina

# get the parameter estimates of null and alternative distribution
density<-SpMixParam(z)

# get the local-FDR estimates by semiparametric mixture
SpMix<-localFDR(z)

# get the local-FDR estimates by normal mixture
normalMix<-nMixParam(z)
```

### multi-dimensional

``` r
library(multiLocalFDR)

# sample data points
Sigma0 <- matrix(c(1, rho12, 
                   rho12, 1), 
                 ncol = 2, byrow = T)
                   
n0 <- rbinom(1, n, p0)
n1 <- n - n0
V <- rCopula(n1, ellipCopula(family = "normal", param = param, dim = 2))

# z-values (probit transformed p-values) from data points

z0 <- rmvnorm(n0, sigma = Sigma0)
z1 <- qnorm(V, mean = 3.5, sd = .5)
z <- rbind(z0, z1)

# get the parameter estimates of null and alternative distribution
density<-SpMixParam(z)

# get the local-FDR estimates by semiparametric mixture
SpMix<-localFDR(z)

# get the local-FDR estimates by normal mixture
normalMix<-nMixParam(z)
```

## Getting help

If you encounter a clear bug, please file an issue with a minimal
reproducible example on
[GitHub](https://github.com/JungiinChoi/multiLocalFDR/issues). For questions and
other discussion, feel free to contact me: Jungin Choi (serimtech07 at snu.ac.kr).

