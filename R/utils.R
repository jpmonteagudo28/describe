#-------------------------------------------#
# Author: JP Monteagudo
# Liberty University, Dept. of Public Health
#-------------------------------------------#
#
# Convert vector of type character to factor
char2factor <- function(data) {

  if(is.vector(data) && is.character(data)){

    data <- as.factor(data)
  } else if (is.data.frame(data) || is.matrix(data)) {

    char_vars <- which(sapply(data, is.character))

    for (var in char_vars) {
      data[[var]] <- as.factor(data[[var]])
    }
  } else {

    stop("Input must be a vector, data frame, or matrix.")
  }
  return(data)
}

#---------------------------------------#
# Create correlated, synthetic normal variables
gen.syn.vars <- function(len, rho, sigma) {
  stopifnot(is.numeric(len),
            is.numeric(rho),
            is.numeric(sigma),
            rho >= -1,
            rho <= 1)

  n <- len
  A <- matrix(c(1, rho, rho, 1), 2,2)
  L <- chol(A)
  eps <- diag(sigma)

  lambda <- t(L) %*% eps
  x <- cbind(rnorm(n),rnorm(n))
  z <- x %*% lambda

  return(z)
}

#----------------------------------------------#
# Generate nth samples for multiple imputation
gen.nth.samples <- function(data,
                            iters = 3,
                            seed = 2899,
                            na.rm = TRUE){
  stopifnot(is.data.frame(data) || is.matrix(data),
            is.numeric(iters),
            is.numeric(seed),
            is.logical(na.rm))

  if(na.rm){
    comp_data <- data[complete.cases(data), ]

    if(nrow(comp_data) == 0) {
      warning("No complete cases found.
              Using original data for imputation.")
      comp_data <- data
    }
  } else {
    comp_data <- data
  }

  make.values <- function(col,
                          iters){
    if (is.factor(col)) {
      obs_vals <- col[!is.na(col)]
      N <- length(col)
      sim_vals <- replicate(iters,
                            sample(levels(col),
                                   N, replace = TRUE,
                                   prob = table(obs_vals) / length(obs_vals))
      )
    } else {
      obs_vals <- col[!is.na(col)]
      N <- length(col)
      sim_vals <- replicate(iters,
                            sample(obs_vals,
                                   N,
                                   replace = TRUE)
      )
    }
    return(sim_vals)
  }

  set.seed(seed = seed)
  imp_data <- replicate(iters,comp_data,simplify = FALSE)

  for (j in 1:ncol(comp_data)) {
    if (any(is.na(comp_data[, j]))) {
      sim_vals <- make.values(comp_data[, j], iters)
      na_index <- which(is.na(comp_data[, j]))
      for (k in 1:iters) {
        imp_data[[k]][na_index, j] <- sim_vals[na_index, k]
      }
    }
  }
  return(imp_data)
}

#--------------------------------------#
# Calculate Euclidean distance for
# rows with NA and set of complete rows
# to use in imputation functions
calc.dist <- function(x,y){
  apply(x,1, function(w){
    sqrt(rowSums((w - y)^2, na.rm = TRUE))
   })
}

#----------------------------------------#
# Generate predicted values from data using
# different regression models
gen.predict <- function(data,
                        varname,
                        robust = FALSE,
                        family = "AUTO",
                        verbose = TRUE) {

  complete_cases <- complete.cases(data)

  if(sum(complete_cases) == 0) {
    if (verbose) {
      warning(sprintf("Imputation not performed. No complete cases for variable: %s", varname))
    }
    return(NULL)
  }

  # Check for sufficient degrees of freedom
  if (sum(complete_cases) <= ncol(data)) {
    if (verbose) {
      warning(sprintf("Stopping missing value prediction. Insufficient degrees of freedom for variable: %s", varname))
    }
    return(NULL)
  }

  # Check for one unique value
  if (length(unique(data[[varname]][complete_cases])) == 1) {
    if (verbose) {
      warning(sprintf("Stopping missing value prediction. Only one unique value for variable: %s", varname))
    }
    return(NULL)
  }

  #Check data distribution family
  if (family == "AUTO") {
    if (is.numeric(data[[varname]])) {
      model_type <- "gaussian"
    } else if (is.factor(data[[varname]])) {
      nLev <- length(levels(data[[varname]]))
      if (nLev == 2) {
        model_type <- "binomial"
      } else {
        model_type <- "multinomial"
      }
    }
  } else {
    model_type <- family
  }

  # Fit model to get cond. mean
  model <- tryCatch({
    if (model_type == "gaussian") {
      if (robust) {
        MASS::rlm(as.formula(paste(varname, "~ .")), data = data[complete_cases,])
      } else {
        stats::lm(as.formula(paste(varname, "~ .")), data = data[complete_cases,])
      }
    } else if (model_type == "binomial") {
      if (robust) {
        robustbase::glmrob(as.formula(paste(varname, "~ .")), data = data[complete_cases,],
                           family = binomial)
      } else {
        stats::glm(as.formula(paste(varname, "~ .")), data = data[complete_cases,],
                   family = binomial)
      }
    } else if(model_type == "multinomial") {
      nnet::multinom(as.formula(paste(varname, "~ .")), data = data[complete_cases,])

    } else if (model_type == "beta") {
      betareg::betareg(as.formula(paste(varname, "~ .")), data = data[complete_cases,])

    } else {
      stats::glm(as.formula(paste(varname, "~ .")), data = data[complete_cases,],
                 family = model_type)
    }
  }, error = function(e) {
    if (verbose) {
      warning(sprintf("Model fitting failed for variable: %s. Error: %s", varname, e$message))
    }
    return(NULL)
  })
  if(is.null(model)) return(NULL)

  predicted <- tryCatch({
    stats::predict(model, newdata = data,
                   type = if (model_type == "gaussian") "response" else "probs")
  }, error = function(e) {
    if (verbose) {
      warning(sprintf("Prediction failed for variable: %s. Error: %s", varname, e$message))
    }
    return(NULL)
  })
  return(predicted)
}

#-------------------------------------------#
# Check remaining missing values in matrix or
# data frame

check.missing <- function(data, post_data = NULL, verbose = TRUE) {


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


