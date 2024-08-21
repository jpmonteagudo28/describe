
<!-- README.md is generated from README.Rmd. Please edit that file -->

# describe <a href="https://describe.jpmonteagudo.com"><img src="man/figures/logo.png" align="right" height="120" alt="describe website" /></a>

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/describe)](https://CRAN.R-project.org/package=describe)
<!-- badges: end -->

The goal of describe is to make it easy to compute descriptive
statistics for non-normal distributions without creating custom
functions or spending hours researching stack overflow. This package
also provides flexible and intuitive univariate and multiple data
imputation and visualizations for exploratory data analysis.

## Installation

You can install the development version of describe like so:

``` r
devtools::install_github("jpmonteagudo28/describe")
```

You can also install it from CRAN like so:

``` r
install.packages("describe")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(describe)

# Generate random data with NAs
data <- gen.mcar(50,rho = .467,sigma = c(1,1.3),n_vars = 2, na_prob = .10)

# Impute missing values using predictive mean matching
imputed_data <- pmean.match(data, robust = TRUE, verbose = TRUE)
#> Initial missing values: 8 (8.00%)
#> Remaining missing values: 0 (0.00%)

head(imputed_data,10)
#>            V1         V2
#> 1   1.7945695 -0.8437137
#> 2   0.8450450 -0.4223777
#> 3   1.5095428  2.0415957
#> 4   0.0367590  0.3120576
#> 5   0.5888666  0.6589248
#> 6   0.4108906 -0.2331285
#> 7  -2.3065881 -3.2270210
#> 8  -0.5356565  0.6062233
#> 9   0.9194851 -0.4223777
#> 10  0.6637492  1.4249718

summary(imputed_data)
#>        V1                 V2         
#>  Min.   :-2.30659   Min.   :-3.2270  
#>  1st Qu.:-0.94678   1st Qu.:-0.4224  
#>  Median : 0.05330   Median : 0.2962  
#>  Mean   :-0.05093   Mean   : 0.1442  
#>  3rd Qu.: 0.66375   3rd Qu.: 0.8352  
#>  Max.   : 2.43697   Max.   : 2.5081
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
# Some examples of bootstrapped correlation coefficients and geometric median
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
