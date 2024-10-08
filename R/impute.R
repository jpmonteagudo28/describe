#' @title
#' Hot deck imputation
#'
#' @description Hot deck imputation (HDI) is a univariate imputation technique where, for each
#' respondent or recipient with a missing value, we find a donor with similar values across
#' a subset of categorical or numerical predictors and use it to fill the recipient's missing
#' observation. For this reason, HDI has been used with stratification across
#' categorical variables.
#'
#' The current implementation of HDI allows the user to choose one of four selection
#' methods: deterministic, random sampling from all possible donors, from k-nearest neighbors and
#' random sampling using weights as probabilities. The function will iteratively impute missing values across
#' all variables with missing observations using the selection method specified in the function arguments.
#'
#' @param data  a matrix or data frame containing missing values in at least one predictor
#' @param method selection method for imputing missing values based on donor similarity.
#' Can be one of:
#' \describe{
#'   \item{\code{"deterministic"}}{Select the same donor value for multiple repetitions of the CDI.}
#'   \item{\code{"rand_from_all"}}{Select a different donor value for each repetition of the CDI.}
#'   \item{\code{"rand_nearest_k"}}{Select one random donor value from a subset of k nearest neighbors for each repetition of the CDI.}
#'   \item{\code{"weighted_rand"}}{Select one random donor through a probability-weighted choice for each repetition of the CDI.}
#' }
#' @param k number of nearest neighbors to select from when using the `rand_nearest_k` method.
#' @param seed a numeric seed for reproducible results for every method except deterministic selection
#' @param na.rm indicates removal of NA values from every row in the matrix or data frame
#'
#' @return a matrix or data frame of imputed values
#' @export
#'
#' @examples
#' data <- gen.mcar(100,rho = c(.56,.23,.18),sigma = c(1,2,.5),n_vars = 3,na_prob = .18)
#' hot_data <- hotdeck.impute(data)
#' summary(hot_data)

hotdeck.impute <- function(data,
                           method = "deterministic",
                           k = NULL,
                           seed = NULL,
                           na.rm = TRUE){

        stopifnot(is.matrix(data)|| is.data.frame(data),
                  is.character(method),
                  is.logical(na.rm))

        if(na.rm){
          complete_data <- data[stats::complete.cases(data),]

          if(nrow(complete_data) == 0) {
            warning("No complete cases found.
              Using original data for imputation.")
            comp_data <- data
          }
        } else {
          complete_data <- data
          warning("Returning original data")
        }

        if(!is.null(seed)){
          set.seed(seed)
        }

        na_index <- which(is.na(data), arr.ind = TRUE)

        # Check for NA values again
        if (nrow(na_index) == 0) {
          return(data)  # No imputation needed
        }

        for(j in unique(na_index[,2])){

          col_na_index <- na_index[na_index[,2] == j, , drop = FALSE]
          complete_cols <- setdiff(1:ncol(data),j)
          na_rows <- data[col_na_index[,1],complete_cols, drop = FALSE]
          donor_rows <- complete_data[, complete_cols, drop = FALSE]
          dist <- calc.dist(na_rows,donor_rows)

          donor_index <- switch(method,
                  "deterministic" = apply(dist, 2, which.min),
                  "rand_from_all" = sample(1:nrow(complete_data),
                                          nrow(col_na_index),
                                          replace = TRUE),
                  "rand_nearest_k" = {
                        apply(dist, 2, function(d) sample(order(d)[1:min(k, length(d))], 1))
                          },
                  "weighted_rand" = {
                        apply(dist, 2, function(d) {
                          weights <- 1 / (d + 1e-8)
                          sample(1:length(d), 1,
                          prob = weights)
                        })
                       }
          )
          data[col_na_index[, 1], j] <- complete_data[donor_index, j]
        }
        return(data)
}

#' @title
#' Cold deck imputation
#'
#' @description Cold deck imputation (CDI) is a univariate imputation technique where, for each
#' respondent or recipient with a missing value, we use an external, pre-existing source to
#' find a donor with similar values across a subset of categorical or numerical predictors
#' and use it to fill the recipient's missing observation. For this reason, cold deck
#' imputation has been used with stratification across categorical variables.
#'
#' The current implementation of CDI allows the user to choose one of four selection
#' methods: deterministic, random sampling from all possible donors, from k-nearest neighbors and
#' random sampling using weights as probabilities. The function will iteratively impute missing values across
#' all variables with missing observations using the selection method specified in the function arguments.
#'
#' @details
#' CDI is a valuable method when a reliable external source of data is available and
#' can ensure more standardized imputations, particularly in large-scale studies or when historical
#' consistency is important. However, it requires careful consideration to avoid introducing bias
#' due to mismatches between the current dataset and the external source.
#'
#' @param data  a matrix or data frame containing missing values in at least one predictor
#' @param ext_data external data source of complete cases to be used as donor values
#' @param method selection method for imputing missing values based on donor similarity.
#' Can be one of:
#' \describe{
#'   \item{\code{"deterministic"}}{Select the same donor value for multiple repetitions of the CDI.}
#'   \item{\code{"rand_from_all"}}{Select a different donor value for each repetition of the CDI.}
#'   \item{\code{"rand_nearest_k"}}{Select one random donor value from a subset of k nearest neighbors for each repetition of the CDI.}
#'   \item{\code{"weighted_rand"}}{Select one random donor through a probability-weighted choice for each repetition of the CDI.}
#' }
#' @param k number of nearest neighbors to select from when using the `rand_nearest_k` method.
#' @param seed a numeric seed for reproducible results for every method except deterministic selection
#' @param na.rm indicates removal of NA values from every row in the matrix or data frame
#'
#' @return a matrix or data frame of imputed values
#' @export
#'
#' @examples
#' data <- gen.mcar(100,rho = c(.56,.23,.18),sigma = c(1,2,.5),n_vars = 3,na_prob = .18)
#' ext_data <- gen.mcar(100,rho = c(.45,.26,.21),sigma = c(1.67,2.23,.56),n_vars = 3,na_prob = 0)
#' cold_data <- coldeck.impute(data,ext_data)
#' summary(cold_data)
coldeck.impute <- function(data,
                           ext_data = NULL,
                           method = "deterministic",
                           k = NULL,
                           seed = NULL,
                           na.rm = TRUE){

  if(is.null(ext_data)){
    stop("`ext_data` must be provided")
  }

  stopifnot(is.data.frame(data) || is.matrix(data),
            is.data.frame(ext_data) || is.matrix(ext_data),
            is.character(method),
            is.logical(na.rm)
  )

  if (na.rm) {
    complete_ext_data <- ext_data[stats::complete.cases(ext_data),]

    if (nrow(complete_ext_data) == 0) {
      stop("No complete cases found in the external data source.")
    }
  } else {
    complete_ext_data <- ext_data
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  na_index <- which(is.na(data), arr.ind = TRUE)

  # Check for NA values in the target data
  if (nrow(na_index) == 0) {
    return(data)  # No imputation needed
  }

  for (j in unique(na_index[, 2])) {

    col_na_index <- na_index[na_index[, 2] == j, , drop = FALSE]
    complete_cols <- setdiff(1:ncol(data), j)
    na_rows <- data[col_na_index[, 1], complete_cols, drop = FALSE]
    donor_rows <- complete_ext_data[, complete_cols, drop = FALSE]
    dist <- calc.dist(na_rows, donor_rows)

    donor_index <- switch(method,
                          "deterministic" = apply(dist, 2, which.min),
                          "rand_from_all" = sample(1:nrow(complete_ext_data),
                                                   nrow(col_na_index),
                                                   replace = TRUE),
                          "rand_nearest_k" = {
                            apply(dist, 2, function(d) sample(order(d)[1:min(k, length(d))], 1))
                          },
                          "weighted_rand" = {
                            apply(dist, 2, function(d) {
                              weights <- 1 / (d + 1e-8)
                              sample(1:length(d), 1,
                                     prob = weights)
                            })
                          }
    )
    data[col_na_index[, 1], j] <- complete_ext_data[donor_index, j]
  }
  return(data)
}

#' @title
#' Predictive mean matching imputation (PMI)
#'
#' @description
#' Predictive mean matching (PMM)
#' is an imputation technique introduced by  Donald . Rubin in 1987. This imputation method aims
#' to maintain the natural variability of the data and avoid implausible imputations that can
#' occur with other univariate imputation methods.
#'
#'
#' @param data a numeric matrix or data frame of at least 2 columns.
#' @param family the distribution family of your observations. The family arguments defaults
#' to 'AUTO'; and it will automatically select a distribution family (gaussian, binomial, multinomial) based on the type of
#' variable (numeric or factor). The distribution family dictates the regression model used (lm,glm, multinom).
#' However, the user can change the family argument to match his response variable distribution
#' and the function will adapt to this input by using the generalized linear model or beta regression.
#' @param robust logical indicated whether to use robust estimation methods or ignore them. If set to 'TRUE',
#' the function will make use of robust linear and generalized linear models to make its prediction.
#' @param k numeric vector indicating number of nearest neighbors to extract for imputation.
#' Currently k defaults to 3 but can be changed.
#' @param seed numeric vector used for reproducible results. Used to sample the same predicted value over time.
#' @param char_to_factor transform character variable to unordered factor variable
#' @param verbose verbose error handling
#'
#' @return a matrix or data frame containing the imputed dataset.
#' @export
#'
#'
#' @details
#' How's predictive mean matching different from conditional mean imputation(CMI)?
#'
#' PMM is a combination of CMI and HDI.
#' Predictive mean matching (PMM) uses regression on observed variables to estimate
#' missing values, like CMI, however, PMM will also fill in the missing value by
#' randomly sampling observed values whose predicted values are closest
#' to the predicted values of the missing observation. This is currently done using
#' the nearest neighbor approach with a set number of neighbors (3) but can be changed
#' depending on your data.
#'
#' @examples
#' set.seed(123)
#' data <- data.frame(x1 = stats::rnorm(100),x2 = stats::rnorm(100),y = stats::rnorm(100))
#' data$x1[sample(1:100, 20)] <- NA
#' data$x2[sample(1:100, 15)] <- NA
#' data$y[sample(1:100, 10)] <- NA
#' fact_dat <- data.frame(data, c = gl(5,20))
#' matched_data <- pmean.match(fact_dat, robust = TRUE)
#' summary(matched_data)
pmean.match <- function(data,
                        family = "AUTO",
                        robust = FALSE,
                        k = 3,
                        char_to_factor = FALSE,
                        seed = NULL,
                        verbose = FALSE){

  if(!is.data.frame(data)){
    data <- as.data.frame(data)
  }

  stopifnot(is.character(family),
            is.logical(robust),
            is.numeric(k),
            is.logical(char_to_factor))

  if(char_to_factor){
    data <- char2factor(data)
    if(verbose){
    warning("Converted character variable to unordered factor variable")
    }
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # find where NA's located
  na_where <- lapply(data,is.na)

  # save copy of OG data for imputation
  dat2impute <- data

  # Find k obs_val based on dist from pred_obs to pred_na
  findk.near <- function(pred_na, pred_obs, obs_val, k = k) {

    dist <- abs(pred_obs - pred_na)
    k_nearest <- order(dist)[1:min(k, length(dist))]
    return(sample(obs_val[k_nearest], 1))
  }

  # Substitute missing values w/ obs_vals sampled from space of length k
  for(var in names(data)){

    na_index <- which(na_where[[var]])

    if(length(na_index) > 0){
      pred_val <- gen.predict(dat2impute, var, robust = robust, family = family)

      if(!is.null(pred_val)) {
        obs_index <- which(!na_where[[var]])
        obs_val <- dat2impute[[var]][obs_index]
        pred_obs <- pred_val[obs_index]

        for(i in na_index) {
          imputed_value <- findk.near(pred_val[i], pred_obs, obs_val, k)
          dat2impute[[var]][i] <- imputed_value
        }
      }
    }
  }
  #Final check for any remaining missing values
  check.missing(data, dat2impute, verbose = verbose)

  return(dat2impute)
}

#' @title
#' Conditional Mean Imputation (CMI)
#'
#' @param data a numeric matrix or data frame of at least 2 columns.
#' @param family the distribution family of your observations. The family arguments defaults
#' to 'AUTO'; and it will automatically select a distribution family (gaussian, binomial, multinomial) based on the type of
#' variable (numeric or factor). The distribution family dictates the regression model used (lm,glm, multinom).
#' However, the user can change the family argument to match his response variable distribution
#' and the function will adapt to this input by using the generalized linear model or beta regression.
#' @param robust logical indicated whether to use robust estimation methods or ignore them. If set to 'TRUE',
#' the function will make use of robust linear and generalized linear models to make its prediction.
#' @param char_to_factor transform character variable to unordered factor variable
#' @param verbose verbose error handling
#'
#' @return a matrix or data frame containing the imputed dataset.
#' @export
#'
#'
#' @description
#' This method replaces missing values with the expected value of the missing variable, given
#' other variables in the dataset. Predictions of conditional mean performed using specific
#' regression models.
#' CMI is a univariate imputation method that leverages the relationship between
#' variables in the data to make informed predictions about missing values. Although
#' this method reduces bias and aims to maintain the relationship among variables, it
#' will yield residuals with less variation than the original data.
#'
#' CMI will fail if a regressor also contains missing values, thus making imputation of
#' the target's missing values impossible. The user will want to understand and fill in
#' as much of the missing data as possible before imputing through this method. If the
#' regressors contain missing values, it is better to use multiple imputation techniques.
#'
#'
#' @examples
#' set.seed(123)
#' data <- data.frame(x1 = c(stats::rnorm(87),rep(NA,13)),
#' x2 = stats::rnorm(100),y = stats::rnorm(100))
#' cmi_data <- cmean.impute(data)
#' summary(cmi_data)

cmean.impute <- function(data,
                      family = "AUTO",
                      robust = FALSE,
                      char_to_factor = FALSE,
                      verbose = FALSE) {

  if (!is.data.frame(data)) {
    data <- as.data.frame(data)
  }

  stopifnot(is.character(family),
            is.logical(robust),
            is.logical(char_to_factor))

  if (char_to_factor) {
    data <- char2factor(data)
    warning("Converted character variable to unordered factor variable")
  }

  # find where NA's located
  na_where <- lapply(data, is.na)

  # save copy of OG data for imputation
  dat2impute <- data

  # Substitute missing values w/ obs_vals
  for (var in names(data)) {
    na_index <- which(na_where[[var]])

    if (length(na_index) > 0) {
      pred_val <-
        gen.predict(dat2impute, var, robust = robust, family = family)

      if (!is.null(pred_val)) {
        # Replace only the missing values with their corresponding predictions
        dat2impute[[var]][na_index] <- pred_val[na_index]
      }
    }
  }
  #Final check for any remaining missing values
  check.missing(data,dat2impute, verbose = verbose)

  return(dat2impute)
}

#' @title
#' Stochastic regression imputation (SRI)
#'
#' @description
#' This method corrects the lack of variability in conditional mean imputation (CMI) by adding
#' an error term to the conditional mean calculation. SRI is more effective than CMI in reducing
#' bias in the imputed values and works well with MCAR and MAR data.
#'
#' @param data a numeric matrix or data frame of at least 2 columns.
#' @param family the distribution family of your observations. The family arguments defaults
#' to 'AUTO'; and it will automatically select a distribution family (gaussian, binomial, multinomial) based on the type of
#' variable (numeric or factor). The distribution family dictates the regression model used (lm,glm, multinom).
#' However, the user can change the family argument to match his response variable distribution
#' and the function will adapt to this input by using the generalized linear model or beta regression.
#' @param tol tolerance,a numeric vector of length 1 used as multiplicative factor to standard deviation for generalized linear models.
#' As the sample size increases, the tolerance value should be decreased to represent the decreasing variability of the sample estimate.
#' @param robust logical indicated whether to use robust estimation methods or ignore them. If set to 'TRUE',
#' the function will make use of robust linear and generalized linear models to make its prediction.
#' @param char_to_factor transform character variable to unordered factor variable
#' @param verbose verbose error handling
#'
#' @return a matrix or data frame containing the imputed dataset.
#' @export
#'
#' @examples
#  set.seed(123)
#' data <- data.frame(x1 = c(stats::rnorm(87),rep(NA,13)),
#' x2 = stats::rnorm(100),y = stats::rnorm(100))
#' sri_data <- stoc.impute(data,tol = 1e-3)
#' summary(sri_data)

stoc.impute <- function(data,
                          family = "AUTO",
                          tol = NULL,
                          robust = FALSE,
                          char_to_factor = FALSE,
                          verbose = FALSE){

  if (!is.data.frame(data)) {
    data <- as.data.frame(data)
  }

  stopifnot(is.character(family),
            is.numeric(tol),
            is.logical(robust),
            is.logical(char_to_factor))

  if (char_to_factor) {
    data <- char2factor(data)
    warning("Converted character variable to unordered factor variable")
  }

  # find where NA's located
  na_where <- lapply(data, is.na)

  # save copy of OG data for imputation
  dat2impute <- data

  # Substitute missing values w/ obs_vals
  for (var in names(data)) {
    na_index <- which(na_where[[var]])

    if (length(na_index) > 0) {
      pred_val <-
        stochastic.predict(dat2impute, var, tol = tol,robust = robust, family = family)

      if (!is.null(pred_val)) {
        # Replace only the missing values with their corresponding predictions
        dat2impute[[var]][na_index] <- pred_val[na_index]
      }
    }
  }
  #Final check for any remaining missing values
  check.missing(data,dat2impute, verbose = verbose)

  return(dat2impute)
}
