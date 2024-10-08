#-----------------------------------------------#
#' @title
#' Generate a data frame with values missing completely at random
#'
#' @description
#' Generate correlated, synthetic normal variables with user-specified probability of MCAR.
#' Specify the column length, correlation coefficient, standard deviation, number of columns
#' and desired probability of missing values to obtain a data frame of correlated observations
#' with missing values.
#'
#' @param len number of rows per column
#' @param rho desired correlation coefficient of generated variables. The length of rho
#' must be equal to the product of `n_vars` and half of `n_vars` minus one.
#' @param sigma desired standard deviation for each generated variable
#' @param n_vars total number of variables to be generated. At least two variables must be provided.
#' @param na_prob desired probability of missingness in each variable set to 10% by default.
#'
#'
#' @return a data frame of at least 2 columns
#' @export
#'
#' @examples
#' syn_na <- gen.mcar(50,c(.25,.75,.044),c(1.1,.56,1.56),3,.47)
#' summary(syn_na)
#'
gen.mcar <- function(len,
                     rho,
                     sigma,
                     n_vars,
                     na_prob = 0.10) {

  stopifnot(is.numeric(len),
            is.numeric(rho),
            is.numeric(sigma),
            rho >= -1,
            rho <= 1,
            length(rho) %in% c(1, (n_vars * (n_vars - 1)) / 2),
            length(sigma) == n_vars,
            n_vars >= 2,
            na_prob >= -1,
            na_prob <= 1)

  n <- len
  A <- diag(n_vars)
  A[upper.tri(A)]<- rho
  A[lower.tri(A)] <- rho
  L <- chol(A)
  eps <- diag(sigma, nrow = n_vars)

  lambda <- L %*% eps
  x <- matrix(stats::rnorm(n * n_vars), nrow = n, ncol = n_vars)
  z <- x %*% lambda
  z <- as.data.frame(z)

  if (na_prob > 0) {
    na_mask <- matrix(stats::rbinom(n * n_vars,1,na_prob) == 1, n, n_vars)
    z[na_mask] <- NA
  }

  return(z)
}

#-----------------------------------------------#
#' @title
#' Generate a data frame with values missing at random
#'
#' @description
#' Generate correlated, synthetic normal variables with user-specified probability of MAR.
#' Specify the column length, correlation coefficient, standard deviation, number of columns
#' and desired probability of missing values to obtain a data frame of correlated observations
#' with missing values.
#'
#' @details
#' The `gen.mar` algorithm will randomly pick a causative feature/column and
#' use the lowest indices of this target column to assign the missing values to
#' the remaining columns. Note that there are additional ways of creating missing
#' values at random, for example, using conditional probabilities and set thresholds,
#' or logistic regression-based masking; however, the truncation algorithm
#' here is a straight-forward implementation of the MAR mechanism.
#'
#'
#' @param len number of rows per column
#' @param rho desired correlation coefficient of generated variables. The length of rho
#' must be equal to the product of `n_vars` and half of `n_vars` minus one.
#' @param sigma desired standard deviation for each generated variable.
#' @param n_vars total number of variables to be generated. At least two variables must be provided.
#' @param na_prob desired probability of missingness in each variable set to 10% by default.
#'
#'
#' @return a data frame of at least 2 columns
#' @export
#'
#' @examples
#' syn_na <- gen.mar(50,c(.25,.75,.044),c(1.1,.56,1.56),3,.15)
#' summary(syn_na)
gen.mar <- function(len,
                    rho,
                    sigma,
                    n_vars,
                    na_prob = 0.10){

  stopifnot(is.numeric(len),
            is.numeric(rho),
            is.numeric(sigma),
            rho >= -1,
            rho <= 1,
            length(rho) %in% c(1, (n_vars * (n_vars - 1)) / 2),
            length(sigma) == n_vars,
            n_vars >= 2,
            na_prob >= 0)

  if(na_prob == 1){
    stop("`na_prob` must be less than 1")
  }

  n <- len
  A <- diag(n_vars)
  A[upper.tri(A)] <- rho
  A[lower.tri(A)] <- rho
  L <- chol(A)
  eps <- diag(sigma, nrow = n_vars)

  lambda <- L %*% eps
  x <- matrix(stats::rnorm(n * n_vars), nrow = n, ncol = n_vars)
  z <- x %*% lambda
  z<- as.data.frame(z)

  if(na_prob > 0){

    # Randomly pick target variable
    target_var <- sample(z,1)

    #index it to be subset
    target_index <- find.index(z,target_var)

    subset_z <- z[,-target_index, drop = FALSE]

    n_rows <- nrow(z)
    n_cols <- ncol(subset_z)

    total_na_needed <- ceiling((n_rows * n_cols * na_prob)/n_cols)

    if(total_na_needed >= n_rows * n_cols){
      stop("Number of NAs required fills or exceeds the total number of cells")
    }

    sort_index <- order(as.matrix(target_var)) # can't order two columns in data frame
    lowest_index <- sort_index[1:total_na_needed]

    z[lowest_index,-target_index] <- NA

  }
  z <- as.data.frame(z)

  return(z)
}

#------------------------------------------------#
#' @title Generate a data frame with values missing
#' not at random
#'
#' @description
#' Generate a dataset containing values missing not at random.
#' Specify the column length, correlation coefficient, standard deviation, number of columns
#' and desired probability of missing values to obtain a data frame of correlated observations
#' with missing values.
#'
#' @details
#' The MNAR mechanism requires missingness to be related to events that are not measured or observed.
#' This type of missingness cannot be determined by statistical analysis and will
#' produce biased estimates. `gen.mnar` uses a self-selection mechanism to create values missing
#' not at random.Others methods, such as logistic regression models, and exponential decay models
#' may be suitable to generate missing values but the self-selection  mechanism provides a
#' straight-forward implementation of the MNAR mechanism.
#'
#' @param len number of rows per column
#' @param rho desired correlation coefficient of generated variables. The length of rho
#' must be equal to the product of `n_vars` and half of `n_vars` minus one.
#' @param sigma desired standard deviation for each generated variable.
#' @param n_vars total number of variables to be generated. At least two variables must be provided.
#' @param na_prob desired probability of missingness in each variable set to 10% by default.
#'
#'
#' @return a data frame of at least 2 columns
#' @export
#'
#' @examples
#' syn_na <- gen.mnar(50,c(.25,.75,.044),c(1.1,.56,1.56),3,.15)
#' summary(syn_na)

gen.mnar <- function(len,
                     rho,
                     sigma,
                     n_vars,
                     na_prob = .10){

  stopifnot(is.numeric(len),
            is.numeric(rho),
            is.numeric(sigma),
            rho >= -1,
            rho <= 1,
            length(rho) %in% c(1, (n_vars * (n_vars - 1)) / 2),
            length(sigma) == n_vars,
            n_vars >= 2,
            na_prob >= 0)

  if(na_prob == 1){
    stop("`na_prob` must be less than 1")
  }

  n <- len
  A <- diag(n_vars)
  A[upper.tri(A)] <- rho
  A[lower.tri(A)] <- rho
  L <- chol(A)
  eps <- diag(sigma, nrow = n_vars)

  lambda <- L %*% eps
  x <- matrix(stats::rnorm(n * n_vars), nrow = n, ncol = n_vars)
  z <- x %*% lambda
  z<- as.data.frame(z)

  if(na_prob > 0){

    # Randomly pick target variable
    target_var <- sample(z,1)

    #index it to be subset
    target_index <- find.index(z,target_var)

    subset_z <- z[,-target_index, drop = FALSE]

    n_rows <- nrow(z)
    n_cols <- ncol(subset_z)
    features <- names(subset_z)

    total_na_needed <- ceiling((n_rows * n_cols * na_prob)/n_cols)

    if(total_na_needed >= n_rows * n_cols){
      stop("Number of NAs required fills or exceeds the total number of cells")
    }

    for(i in features){

      feature_col <- subset_z[,i]

      sorted_index <- order(feature_col)
      lowest_index <- sorted_index[1:total_na_needed]

      z[lowest_index,i] <- NA
    }

  }
  return(z)
}

#-------------------------------------------------#
#' @title
#'  Transform data to missing values completely at random
#'
#' @description
#'  Transform a complete case dataset according to the MCAR mechanism.The MCAR
#'  mechanism assumes that the probability of missingness is the same for all cases,
#'  and therefore,the missing data are not related to the observed data. `transform.mcar`
#'  uses the binomial distribution to generate NA  an index of possible NA values that will
#'  replace a percentage of data points in the input.
#'
#' @param input data to transform using MCAR mechanism
#' @param na_prob desired probability of NA count defined as the probability of success
#' in a Bernoulli trial.
#'
#' @return a data frame or matrix containing NAs
#' @export
#'
#'
#' @examples
#' set.seed(123)
#' data <- data.frame(x1 = stats::rnorm(100),x2 = stats::rnorm(100),y = stats::rnorm(100))
#' mcar_data <- mcar.transform(data,na_prob = .25)
#'
#' summary(mcar_data)
#'

mcar.transform <- function(input,na_prob){

  if(!is.data.frame(input) && !is.matrix(input)){
    input <- as.matrix(input)
  }

  # Create copy of data to output
  output <- as.matrix(input)
  n_rows <- nrow(input)
  n_cols <- ncol(input)

  # Set threshold for max. number of NA permissible
  threshold <- 1

  na_mask <- matrix(stats::rbinom(n_rows * n_cols,1,na_prob) == 1, n_rows,n_cols)
  output[na_mask] <- NA

  col_na_risk <- apply(output, 2, function(col) sum(is.na(col)) >= (n_cols - threshold))
  row_na_risk <- apply(output, 1, function(row) sum(is.na(row)) >= (n_rows - threshold))
  zero_missing <- sum(is.na(output))

  if(any(col_na_risk)){
    warning("One or more columns are at risk of being entirely NA.")
  }
  if(any(row_na_risk)){
    warning("One or more rows are at risk of being entirely NA.")
  }

  if(zero_missing == 0){
    warning("No missing values were generated for the data")
  }

  output <- as.data.frame(output)
  return(output)
}


#-------------------------------------------------#
#' @title Transform  data to missing values at random.
#'
#' @description
#' Transform a complete case dataset according to the MAR mechanism.The MAR
#' mechanism assumes that the probability of missingness in a variable depends
#' on the observed data but not on the missing data itself. This function introduces
#' missing values in selected features of a dataset, with the missing values determined
#' by the values of the causative (target) feature.
#'
#' @param input data to transform using MAR mechanism
#' @param target variable to be used as causative feature
#' @param features variables to which NA values are introduced using a causative feature
#' @param na_rate proportion of missing values to be added to data
#'
#' @return a matrix or data frame containing NAs
#' @export
#'
#'
#' @examples
#' set.seed(123)
#' data <- gen.mcar(100,rho = c(.15,.25,.12,.45,.34,.54),sigma = c(1,2,1,2),n_vars = 4, na_prob = 0)
#'
#' mar_data <-mar.transform(data,"V1",c("V2","V3"), na_rate = .25)
#' summary(mar_data)

mar.transform <- function(input,
                          target,
                          features,
                          na_rate){

  if(!is.data.frame(input) && !is.matrix(input)){
    stop("Input must be a data frame or matrix")
  }

  # Create copy of data to output and subset to variables of interest
  output <- as.matrix(input)
  subset_output <- as.matrix(output[,c(target,features)])
  target_var <- subset_output[,target]


  n_rows <- nrow(output)
  n_cols <- ncol(subset_output)
  n_feat <- length(features)

  # Total number of NAs needed
  total_na_needed <- ceiling((n_rows * n_cols * na_rate)/n_feat)

  if(total_na_needed > n_rows){
    stop("Number of features is too low for specified na_rate")
  }

  if(n_feat >= n_cols){
    stop("Number of features exceeds current number of columns in data")
  }

  sort_index <- order(target_var)
  lowest_index <- sort_index[1:total_na_needed]

  # Set threshold for max. number of NA permissible computation
  threshold <- 1

  subset_output[lowest_index,features] <- NA
  output[,features] <- subset_output[,features]

  col_na_risk <- apply(output[, features, drop = FALSE], 2, function(col) sum(is.na(col)) >= (n_feat - threshold))
  row_na_risk <- apply(output, 1, function(row) sum(is.na(row)) >= (n_rows - threshold))

  if (any(col_na_risk)) {
    warning("One or more features are at risk of being entirely NA.")
  }
  if (any(row_na_risk)) {
    warning("One or more rows are at risk of being entirely NA.")
  }

  # Convert the output back to a data frame if necessary and return it
  output <- as.data.frame(output)
  return(output)
}
#--------------------------------------------------#
#' @title
#' Transform data to missing values not at random
#'
#' @description
#' The MNAR mechanism assumes that the probability and cause of missingness are unknown to us and lie
#' in the unobserved data. Under this mechanism we generate missing values by looking at each
#' variable's lowest values and replacing them with NA.
#'
#' @param input data to transform using MAR mechanism
#' @param target variable to be used as causative feature
#' @param features variables to which NA values are introduced using a causative feature
#' @param na_rate proportion of missing values to be added to data
#'
#' @return a matrix or data frame containing NAs
#' @export
#'
#' @examples
#' set.seed(123)
#' data <- gen.mcar(100,rho = c(.15,.25,.12,.45,.34,.54),sigma = c(1,2,1,2),n_vars = 4, na_prob = 0)
#'
#' mnar_data <- mnar.transform(data,"V1",c("V2","V3"), na_rate = .25)
#' summary(mnar_data)

mnar.transform <- function(input,
                           target,
                           features,
                           na_rate){

  if(!is.data.frame(input) && !is.matrix(input)){
    stop("Input must be a data frame or matrix")
  }

  # Create copy of data to output and subset to variables of interest
  output <- as.matrix(input)
  subset_output <- as.matrix(output[,features])



  n_rows <- nrow(output)
  n_cols <- ncol(output)
  n_feat <- length(features)

  # Total number of NAs needed
  total_na_needed <- ceiling((n_rows * n_cols * na_rate)/n_feat)

  if(total_na_needed > n_rows){
    stop("Number of features is too low for specified na_rate")
  }

  if(n_feat >= n_cols){
    stop("Number of features exceeds current number of columns in data")
  }

  for(i in features){

    feature_col <- subset_output[,i]

  sorted_index <- order(feature_col)
  lowest_index <- sorted_index[1:total_na_needed]

  output[lowest_index,i] <- NA
  }

  output <- as.data.frame(output)
  return(output)
}
