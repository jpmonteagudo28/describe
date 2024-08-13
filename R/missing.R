#' @description
#' Check percent of missing values in data frame or matrix pre and post-processing
#'
#' @param data data frame or matrix prior to processing
#' @param post_data data frame or matrix after processing
#' @verbose indicate verbose warning and error messages
#'
#' @return a character vector specifying the percent of missing values after processing
#'
#' @details
#' `check.missing` takes a data frame and compares the percentage of missing values
#' pre and post-processing within one imputation function. If missing values completely
#' replaced by imputed values, the remaining percent missing should be 0, otherwise
#' the user will see a warning indicating imputation wasn't performed or the regressors
#' used in the model contained missing values.
#'
#' @keywords internal
#'
#' @examples
#' \donttest
#' set.seed(123)
#' data <- data.frame(x1 = rnorm(100),x2 = rnorm(100),y = rnorm(100))
#'
#' Introduce missing values
#' data$x1[sample(1:100, 20)] <- NA
#' data$x2[sample(1:100, 15)] <- NA
#' data$y[sample(1:100, 10)] <- NA
#'
#' check.missing(data, data, verbose = TRUE)
#'

check.missing <- function(data, post_data = NULL, verbose = FALSE) {


  tot_initial_obs <- prod(dim(data))
  initial_na <- sum(is.na(data))
  tot_na_prcnt <- (initial_na / tot_initial_obs) * 100


  if (is.null(post_data)) {
    post_data <- data
  }


  tot_obs <- prod(dim(post_data))
  remain_na_count <- sum(is.na(post_data))
  remain_na_prcnt <- (remain_na_count / tot_obs) * 100


  if (verbose) {
    cat(sprintf("Initial missing values: %d (%.2f%%)\n", initial_na, tot_na_prcnt))
    cat(sprintf("Remaining missing values: %d (%.2f%%)\n", remain_na_count, remain_na_prcnt))
  }


  if (remain_na_count > 0) {
    message_text <- sprintf("Initial missing: %.2f%%. Remaining missing: %.2f%%.", tot_na_prcnt, remain_na_prcnt)
    if (verbose) {
      warning(sprintf("%s Imputation not performed or missing values in the regressors.", message_text))
    } else {
      warning("Imputation not performed or missing values in the regressors.")
    }
  }
}

#' @description MAR check using logistic regression
#'
#' @details We can informally test for data missing at random (MAR) by
#' creating a binary variable that represents missing data (1) and non-missing
#' data (0). We then perform a logistic regression using the binary variable as
#' our target to obtain p-values. Pairwise comparisons are performed on each unique
#' pair and p-values are obtained after applying the Dunn-Šidák correction.
#'
#' This test is not a formal test used to check MAR. Additionally, tests for missingness
#' may not be very practical and do not substitute knowledge of the field and other much
#' better joint tests dealing with missingness.
#'
#' @param data a data frame or matrix of at least two columns
#' @param digits significant figures used in decimals
#'
#' @return a square matrix with p-values across pairwise comparisons
#' @export
#'
#' @examples
#' set.seed(123)
#' data <- data.frame(x1 = rnorm(100),x2 = rnorm(100),y = rnorm(100))
#'
#' Introduce missing values
#' data$x1[sample(1:100, 20)] <- NA
#' data$x2[sample(1:100, 15)] <- NA
#' data$y[sample(1:100, 10)] <- NA
#'
#' check.MAR(data, digits = 2)

check.MAR <- function(data,
                      digits = 3){

  stopifnot(is.numeric(digits))

  if(!is.data.frame(data)&& !is.matrix(data))
    stop("Data must be a data frame or matrix of at least two columns")

  data <- na.col.add(data)
  p_values <- numeric(ncol(data))
  nvars <- ncol(data)

  result <- matrix(NA, nrow = nvars,
                   ncol = nvars,
                   dimnames = list(colnames(data),
                                   colnames(data)))

  for (i in 1:(nvars - 1)) {
    for (j in (i + 1):nvars) {

      glm_summ <- summary(glm(data$is_na ~ data[, i] + data[, j],
                              data, family = binomial))

      p_value <- glm_summ$coefficients[3, 4]


      adj_p_value <- use.sidak(p_value, nvars)

      result[i, j] <- adj_p_value
      result[j, i] <- adj_p_value
    }
  }
  result <- round(result, digits)

  return(result)
}

#' @description
#' Little's test for missing completely at random (MCAR)
#'
#' @details
#' In Little's test of MCAR the data are modeled as multidimensional, multivariate normal
#' with mean 'mu' and covariance matrix 'sigma'.
#' The test statistic is the sum of the squared standardized differences between sub-sample
#' means and the expected population means weighted by the variance-covariance matrix and the
#' number of observations. Under the null hypothesis, the test statistic follows a chi-square distribution with \eqn{\sum k_j - k}
#' degrees of freedom, where \eqn{k_j} is the number of complete variables for missing data
#' pattern \eqn{j}, and \eqn{k} is the total number of variables. A statistically
#' significant result provides evidence against MCAR.
#'
#' When the normality assumption is not  satisfied, the test will work for quantitative
#' random variables but not categorical ones. Additionally, the test will not specify which
#' variable(s) are not MCAR, and it fails to identify collinearity among variables. Third, the
#' test can neither prove the MCAR assumption nor rule out the hypothesis of MNAR.
#'
#' @param data a data frame or matrix of at least two columns
#' @param digits  significant figures used in decimals
#'
#' @return a numeric matrix
#' \item{statistic}{Chi-squared statistic for Little's test}
#' \item{df}{Degrees of freedom used to compute chi-square statistic}
#' \item{p_value}{P-value for the chi-square statistic}
#' \item{pattern}{Unique missing data patterns found}
#' @reference
#' Little, R. J. A. (1988). A test of Missing Completely at Random for multivariate data with missing values. Journal of the American Statistical Association, 83, 1198-1202. https://doi.org/10.2307/2290157
#'
#'code is adapted from Eric Stemmler: \url{https://web.archive.org/web/20201120030409/https://stats-bayes.com/post/2020/08/14/r-function-for-little-s-test-for-data-missing-completely-at-random/}
#' and naniar's `mcar_test`.
#'
#' @examples
#' set.seed(123)
#' data <- gen.syn.vars(100,rho = c(.8,.53,.23,-.13,.35,-.56),
#'                      sigma = c(1,2,1,2), n_vars = 4, na_prob = .17)
#'
#'
#'
check.MCAR <- function(data,
                       digits = 3) {

  stopifnot(is.numeric(digits))

  if (!is.data.frame(data) && !is.matrix(data))
    stop("Data must be a data frame or matrix of at least two columns")

  # Convert to matrix if it's a data frame
  if (is.data.frame(data)) {
    data <- data.matrix(data)
  }

  n_var <- ncol(data)
  varnames <- colnames(data)

  # Add pattern column to data
  na_pattern <- pattern.rank(data)
  miss_data <- cbind(data, na_pattern)

  # Maximum likelihood estimation from {norm}
  s <- norm::prelim.norm(data)
  ll <- norm::em.norm(s, showits = FALSE)
  fit <- norm::getparam.norm(s = s, theta = ll)
  grand_mean <- fit$mu
  grand_cov <- fit$sigma
  colnames(grand_cov) <- rownames(grand_cov) <- varnames

  # Split data by missing pattern
  split_data <- split(data.frame(miss_data), na_pattern)

  # Initialize vectors to store results
  n_patterns <- length(split_data)
  d2 <- numeric(n_patterns)
  kj <- numeric(n_patterns)

  # Loop over each group
  for (i in seq_len(n_patterns)) {


  group_data <- split_data[[i]][, -ncol(split_data[[i]])]

  # Calculate kj for the current pattern
  kj[i] <- sum(1 * !is.na(colSums(group_data)))

  mu <- colMeans(group_data, na.rm = TRUE) - grand_mean
  keep <- !is.na(mu)
  mu <- mu[keep]
  sigma <- grand_cov[keep, keep, drop = FALSE]

  d2[i] <- nrow(group_data) * (t(mu) %*% solve(sigma) %*% mu)
}

# Aggregate results
total_d2 <- sum(d2)
total_kj <- sum(kj)
df <- total_kj - n_var

p_value <- tryCatch(
  stats::pchisq(total_d2, df, lower.tail = FALSE),
  error = function(e) {
    warning("Error in p-value calculation: ", e$message)
    return(NA)
  }
)

  # Output the results
  result <- data.frame(
    statistic = round(total_d2, digits = digits),
    df = df,
    p.value = round(p_value, digits = digits),
    missing.patterns = max(na_pattern)
  )

  return(result)
}




