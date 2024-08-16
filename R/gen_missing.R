#-----------------------------------------------#
# Generate correlated, synthetic normal variables
# with user-specified probability of MCAR
# under the Bernoulli distribution.
gen.mcar <- function(len,
                     rho,
                     sigma,
                     n_vars,
                     na_prob = 0) {

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
  x <- matrix(rnorm(n * n_vars), nrow = n, ncol = n_vars)
  z <- x %*% lambda
  z <- as.data.frame(z)

  # Introduce NAs based on na_prob
  if (na_prob > 0) {
    na_mask <- matrix(rbinom(n * n_vars,1,na_prob) == 1, n, n_vars)
    z[na_mask] <- NA
  }

  return(z)
}

#-----------------------------------------------#
# Generate correlated, synthetic normal variables
# with user-specified probability of MAR
gen.mar <- function(len,
                    rho,
                    sigma,
                    n_vars,
                    na_prob){

}



#------------------------------------------------#
# Generate correlated, synthetic normal variables
# with user-specified probability of MNAR.



#-------------------------------------------------#
# Transform non-missing data to missing values
# with user-specified probability of MCAR

transform.mcar <- function(input,na_prob){

  if(!is.data.frame(input) && !is.matrix(input)){
    input <- as.matrix(input)
  }

  # Create copy of data to output
  output <- as.matrix(input)
  n_rows <- nrow(input)
  n_cols <- ncol(input)

  # Set threshold for max. number of NA permissible
  threshold <- 1

  na_mask <- matrix(rbinom(n_rows * n_cols,1,na_prob) == 1, n_rows,n_cols)
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
#' @description Transform non-missing data to missing values
#' with user-specified probability of MAR
#'
#' @details
#' Transform a complete case dataset according to the MAR mechanism.The MAR
#' mechanism assumes that the probability of missingness in a variable depends
#' on the observed data but not on the missing data itself. This function introduces
#' missing values in selected features of a dataset, with the missing values determined
#' by the values of the causative (target) feature.
#'
#'

transform.mar <- function(input,target,features, na_prob){

  if(!is.data.frame(input) && !is.matrix(input)){
    stop("Input must be a data frame or matrix")
  }

  # Create copy of data to output and subset to variables of interest
  output <- as.matrix(input)
  subset_output <- as.matrix(output[,c(target,features)])
  target_var <- subset_output[,target]


  n_rows <- nrow(subset_output)
  n_cols <- ncol(subset_output)
  n_feat <- length(features)

  # Total number of NAs needed
  total_na_needed <- (n_rows * n_cols * na_prob)/n_feat

  if(total_na_needed > n_rows){
    stop("NUmber of features is too low for specified na_prob")
  }

  if(n_feat >= n_cols){
    stop("Number of features exceeds current number of columns in data")
  }

  sort_index <- order(target_var)
  lowest_index <- sort_index[1:total_na_needed]

  # Set threshold for max. number of NA permissible
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
# Transform non-missing data to missing values
# with user-specified probability MNAR.
