
# multiLocalFDR 


## Overview

multiLocalFDR is a package for multi-dimensional local-FDR estimation using a semiparametric mixture method.
The two pillars of the proposed approach are Efron's empirical null principle and log-concave density estimation for the alternative distribution. A unique feature of our method is that it can be extended to compute the local false discovery rates by combining multiple lists of p-values.

  - `localFDR()` provides estimates of local-fdr for given lists of p-values.
  - `SPestimate()` provides estimates of null and alternative distribution of our method.
  - `normmix()` provides estimates of null and alternative distribution of normal mixture model. 
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

``` r
library(multiLocalFDR)

starwars %>% 
  filter(species == "Droid")
#> # A tibble: 6 x 14
#>   name   height  mass hair_color skin_color  eye_color birth_year sex   gender  
#>   <chr>   <int> <dbl> <chr>      <chr>       <chr>          <dbl> <chr> <chr>   
#> 1 C-3PO     167    75 <NA>       gold        yellow           112 none  masculi…
#> 2 R2-D2      96    32 <NA>       white, blue red               33 none  masculi…
#> 3 R5-D4      97    32 <NA>       white, red  red               NA none  masculi…
#> 4 IG-88     200   140 none       metal       red               15 none  masculi…
#> 5 R4-P17     96    NA none       silver, red red, blue         NA none  feminine
#> # … with 1 more row, and 5 more variables: homeworld <chr>, species <chr>,
#> #   films <list>, vehicles <list>, starships <list>

starwars %>% 
  select(name, ends_with("color"))
#> # A tibble: 87 x 4
#>   name           hair_color skin_color  eye_color
#>   <chr>          <chr>      <chr>       <chr>    
#> 1 Luke Skywalker blond      fair        blue     
#> 2 C-3PO          <NA>       gold        yellow   
#> 3 R2-D2          <NA>       white, blue red      
#> 4 Darth Vader    none       white       yellow   
#> 5 Leia Organa    brown      light       brown    
#> # … with 82 more rows

starwars %>% 
  mutate(name, bmi = mass / ((height / 100)  ^ 2)) %>%
  select(name:mass, bmi)
#> # A tibble: 87 x 4
#>   name           height  mass   bmi
#>   <chr>           <int> <dbl> <dbl>
#> 1 Luke Skywalker    172    77  26.0
#> 2 C-3PO             167    75  26.9
#> 3 R2-D2              96    32  34.7
#> 4 Darth Vader       202   136  33.3
#> 5 Leia Organa       150    49  21.8
#> # … with 82 more rows

starwars %>% 
  arrange(desc(mass))
#> # A tibble: 87 x 14
#>   name    height  mass hair_color skin_color  eye_color  birth_year sex   gender
#>   <chr>    <int> <dbl> <chr>      <chr>       <chr>           <dbl> <chr> <chr> 
#> 1 Jabba …    175  1358 <NA>       green-tan,… orange          600   herm… mascu…
#> 2 Grievo…    216   159 none       brown, whi… green, ye…       NA   male  mascu…
#> 3 IG-88      200   140 none       metal       red              15   none  mascu…
#> 4 Darth …    202   136 none       white       yellow           41.9 male  mascu…
#> 5 Tarfful    234   136 brown      brown       blue             NA   male  mascu…
#> # … with 82 more rows, and 5 more variables: homeworld <chr>, species <chr>,
#> #   films <list>, vehicles <list>, starships <list>

starwars %>%
  group_by(species) %>%
  summarise(
    n = n(),
    mass = mean(mass, na.rm = TRUE)
  ) %>%
  filter(
    n > 1,
    mass > 50
  )
#> # A tibble: 8 x 3
#>   species      n  mass
#>   <chr>    <int> <dbl>
#> 1 Droid        6  69.8
#> 2 Gungan       3  74  
#> 3 Human       35  82.8
#> 4 Kaminoan     2  88  
#> 5 Mirialan     2  53.1
#> # … with 3 more rows
```

## Getting help

If you encounter a clear bug, please file an issue with a minimal
reproducible example on
[GitHub](https://github.com/JungiinChoi/multiLocalFDR/issues). For questions and
other discussion, feel free to contact me: Jungin Choi (serimtech07 at snu.ac.kr).

