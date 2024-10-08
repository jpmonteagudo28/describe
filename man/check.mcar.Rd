% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_missing.R
\name{check.mcar}
\alias{check.mcar}
\title{Little's test for missing completely at random (MCAR)}
\usage{
check.mcar(data, digits = 3)
}
\arguments{
\item{data}{a data frame or matrix of at least two columns}

\item{digits}{significant figures used in decimals}
}
\value{
a numeric matrix
\item{statistic}{Chi-squared statistic for Little's test}
\item{df}{Degrees of freedom used to compute chi-square statistic}
\item{p_value}{P-value for the chi-square statistic}
\item{pattern}{Unique missing data patterns found}
}
\description{
In Little's test of MCAR the data are modeled as multidimensional, multivariate normal
with mean 'mu' and covariance matrix 'sigma'.
The test statistic is the sum of the squared standardized differences between sub-sample
means and the expected population means weighted by the variance-covariance matrix and the
number of observations. Under the null hypothesis, the test statistic follows a chi-square distribution with \eqn{\sum k_j - k}
degrees of freedom, where \eqn{k_j} is the number of complete variables for missing data
pattern \eqn{j}, and \eqn{k} is the total number of variables. A statistically
significant result provides evidence against MCAR.

When the normality assumption is not  satisfied, the test will work for quantitative
random variables but not categorical ones. Additionally, the test will not specify which
variable(s) are not MCAR, and it fails to identify collinearity among variables. Third, the
test can neither prove the MCAR assumption nor rule out the hypothesis of MNAR.
}
\examples{
set.seed(123)
data <- data.frame(x1 = stats::rnorm(100),x2 = stats::rnorm(100),y = stats::rnorm(100))
data$x1[sample(1:100, 20)] <- NA
data$x2[sample(1:100, 15)] <- NA
data$y[sample(1:100, 10)] <- NA
check.mcar(data,digits = 3)
}
\references{
Little, R. J. A. (1988). A test of Missing Completely at Random for multivariate data with missing values. Journal of the American Statistical Association, 83, 1198-1202. https://doi.org/10.2307/2290157

code is adapted from Eric Stemmler: \url{https://web.archive.org/web/20201120030409/https://stats-bayes.com/post/2020/08/14/r-function-for-little-s-test-for-data-missing-completely-at-random/}
and naniar's \code{mcar_test}.
}
