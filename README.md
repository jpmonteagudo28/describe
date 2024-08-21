
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
#> Initial missing values: 11 (11.00%)
#> Remaining missing values: 0 (0.00%)

head(imputed_data,10)
#>            V1         V2
#> 1  -1.4127223 -0.5205037
#> 2   0.3090218  1.4865998
#> 3   0.3240344  0.6016828
#> 4  -0.9598153 -0.2362346
#> 5  -1.4410359 -1.2159774
#> 6  -0.4378229 -2.3591382
#> 7   0.3059329  1.0015804
#> 8  -1.0438890 -1.0072981
#> 9  -1.7135055 -1.3785545
#> 10 -0.4770137 -0.9220649

summary(imputed_data)
#>        V1                V2         
#>  Min.   :-2.3115   Min.   :-3.4326  
#>  1st Qu.:-1.0948   1st Qu.:-1.2160  
#>  Median :-0.1535   Median :-0.3153  
#>  Mean   :-0.2703   Mean   :-0.2040  
#>  3rd Qu.: 0.3237   3rd Qu.: 0.9156  
#>  Max.   : 2.8761   Max.   : 2.7516
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
