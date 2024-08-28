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

#----------------------------------------------#
# Generate nth samples for multiple imputation
gen.nth.samples <- function(data,
                            iters = 3,
                            seed = NULL,
                            na.rm = TRUE){

  stopifnot(is.data.frame(data) || is.matrix(data),
            is.numeric(iters),
            is.logical(na.rm))

  if(is.null(seed)){
    stop("Provide a numeric seed argument for reproducible results.")
  }


  if(na.rm){
    complete_data <- data[stats::complete.cases(data), ]

    if(nrow(complete_data) == 0) {
      warning("No complete cases found.
              Using original data for sample generation.")
      comp_data <- data
    }
  } else {
    complete_data <- data
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
  imp_data <- replicate(iters,complete_data,simplify = FALSE)

  for (j in 1:ncol(complete_data)) {
    if (any(is.na(complete_data[, j]))) {
      sim_vals <- make.values(complete_data[, j], iters)
      na_index <- which(is.na(complete_data[, j]))
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
                        var_name,
                        robust = FALSE,
                        family = "AUTO",
                        verbose = TRUE) {

  complete_cases <- stats::complete.cases(data)

  if(sum(complete_cases) == 0) {
    if (verbose) {
      warning(sprintf("Imputation not performed. No complete cases for variable: %s", var_name))
    }
    return(NULL)
  }

  # Check for sufficient degrees of freedom
  if (sum(complete_cases) <= ncol(data)) {
    if (verbose) {
      warning(sprintf("Stopping missing value prediction. Insufficient degrees of freedom for variable: %s", var_name))
    }
    return(NULL)
  }

  # Check for one unique value
  if (length(unique(data[[var_name]][complete_cases])) == 1) {
    if (verbose) {
      warning(sprintf("Stopping missing value prediction. Only one unique value for variable: %s", var_name))
    }
    return(NULL)
  }

  #Check data distribution family
  if (family == "AUTO") {
    if (is.numeric(data[[var_name]])) {
      model_type <- "gaussian"
    } else if (is.factor(data[[var_name]])) {
      nLev <- length(levels(data[[var_name]]))
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
        MASS::rlm(stats::as.formula(paste(var_name, "~ .")), data = data[complete_cases,])
      } else {
        stats::lm(stats::as.formula(paste(var_name, "~ .")), data = data[complete_cases,])
      }
    } else if (model_type == "binomial") {
      if (robust) {
        robustbase::glmrob(stats::as.formula(paste(var_name, "~ .")), data = data[complete_cases,],
                           family = "binomial")
      } else {
        stats::glm(stats::as.formula(paste(var_name, "~ .")), data = data[complete_cases,],
                   family = "binomial")
      }
    } else if(model_type == "multinomial") {
      nnet::multinom(stats::as.formula(paste(var_name, "~ .")), data = data[complete_cases,])

    } else if (model_type == "beta") {
      betareg::betareg(stats::as.formula(paste(var_name, "~ .")), data = data[complete_cases,])

    } else {
      stats::glm(stats::as.formula(paste(var_name, "~ .")), data = data[complete_cases,],
                 family = model_type)
    }
  }, error = function(e) {
    if (verbose) {
      warning(sprintf("Model fitting failed for variable: %s. Error: %s", var_name, e$message))
    }
    return(NULL)
  })
  if(is.null(model)) return(NULL)

  predicted <- tryCatch({
    stats::predict(model, newdata = data,
                   type = if (model_type == "gaussian") "response" else "probs")
  }, error = function(e) {
    if (verbose) {
      warning(sprintf("Prediction failed for variable: %s. Error: %s", var_name, e$message))
    }
    return(NULL)
  })
  return(predicted)
}

#-------------------------------------------#
# Generate predicted values for stochastic
# regression imputation
stochastic.predict <- function(data,
                        var_name,
                        tol = NULL,
                        robust = FALSE,
                        family = "AUTO",
                        verbose = TRUE) {

  complete_cases <- stats::complete.cases(data)

  if(sum(complete_cases) == 0) {
    if (verbose) {
      warning(sprintf("Imputation not performed. No complete cases for variable: %s", var_name))
    }
    return(NULL)
  }

  # Check for sufficient degrees of freedom
  if (sum(complete_cases) <= ncol(data)) {
    if (verbose) {
      warning(sprintf("Stopping missing value prediction. Insufficient degrees of freedom for variable: %s", var_name))
    }
    return(NULL)
  }

  # Check for one unique value
  if (length(unique(data[[var_name]][complete_cases])) == 1) {
    if (verbose) {
      warning(sprintf("Stopping missing value prediction. Only one unique value for variable: %s", var_name))
    }
    return(NULL)
  }

  # Determine the model type based on the family argument
  if (family == "AUTO") {
    if (is.numeric(data[[var_name]])) {
      model_type <- "gaussian"
    } else if (is.factor(data[[var_name]])) {
      nLev <- length(levels(data[[var_name]]))
      model_type <- if (nLev == 2) "binomial" else "multinomial"
    }
  } else {
    model_type <- family
  }

  # Fit the model based on the determined type
  model <- tryCatch({
    if (model_type == "gaussian") {
      if (robust) {
        MASS::rlm(stats::as.formula(paste(var_name, "~ .")), data = data[complete_cases,])
      } else {
        stats::lm(stats::as.formula(paste(var_name, "~ .")), data = data[complete_cases,])
      }
    } else if (model_type == "binomial") {
      if (robust) {
        robustbase::glmrob(stats::as.formula(paste(var_name, "~ .")), data = data[complete_cases,],
                           family = "binomial")
      } else {
        stats::glm(stats::as.formula(paste(var_name, "~ .")), data = data[complete_cases,],
                   family = "binomial")
      }
    } else if (model_type == "multinomial") {
      nnet::multinom(stats::as.formula(paste(var_name, "~ .")), data = data[complete_cases,])

    } else if (model_type == "beta") {
      betareg::betareg(stats::as.formula(paste(var_name, "~ .")), data = data[complete_cases,])

    } else {
      stats::glm(stats::as.formula(paste(var_name, "~ .")), data = data[complete_cases,],
                 family = model_type)
    }
  }, error = function(e) {
    if (verbose) {
      warning(sprintf("Model fitting failed for variable: %s. Error: %s", var_name, e$message))
    }
    return(NULL)
  })

  if (is.null(model)) return(NULL)

  # Predict values for the variable or dataset
  predicted <- tryCatch({
    stats::predict(model, newdata = data,
                   type = if (model_type == "gaussian") "response" else "probs")
  }, error = function(e) {
    if (verbose) {
      warning(sprintf("Prediction failed for variable: %s. Error: %s", var_name, e$message))
    }
    return(NULL)
  })

  # residual standard deviation using sigma() for Gaussian models
  residual_sd <- if (model_type == "gaussian") {
    stats::sigma(model)
  } else {
    #  in other models, apply a tolerance to the residual standard deviation
    residuals <- stats::residuals(model)
    tol * sqrt(sum(residuals^2) / stats::df.residual(model))
  }

  missing_indices <- which(is.na(data[[var_name]]))

  # Add stochastic noise to predictions
  if (length(missing_indices) > 0) {
    predicted[missing_indices] <- predicted[missing_indices] +
      stats::rnorm(length(missing_indices), mean = 0, sd = residual_sd)
  }

  # Return the predicted values with stochastic imputation
  return(predicted)
}

#-------------------------------------------#
# Reshape data frame by creating binary column
# indicating whether rows contain missing values.

na.col.add <- function(x) {

  if(!is.data.frame(x) || !is.matrix(x)){
    x <- as.data.frame(x)
  }

  is_na <- ifelse(rowSums(is.na(x)) > 0,1,0)
   x <- cbind(x,is_na)

   return(x)
}

#----------------------------------------------#
# Use a Šidák correction to reduce family-wise error
# in model comparisons for MCAR/MAR check. The Šidák
# correction assumes that each comparison is independent.
# The correction assumes that for m variables, there will
# be m(m-1) tests.

use.sidak <- function(p_value,n_vars){

  stopifnot(is.numeric(p_value),
            is.numeric(n_vars))

  adj_sidak <- (1 - (1 - p_value)^(n_vars*(n_vars-1)))

  return(adj_sidak)
}

#-----------------------------------------------#
# Create custom pattern building and ranking
# function
pattern.rank <- function(data) {

  if(!is.data.frame(data)){
    data <- as.data.frame(data)
  }

  n_rows <- nrow(data)
  n_cols <- ncol(data)

  # Step 1: Convert missing patterns to integers using the custom
  # bitwise shift operator
  patterns <- integer(n_rows)
  for (i in 1:n_cols) {
    patterns <- patterns + (as.integer(is.na(data[[i]])) %<<% (i - 1))
  }

  # Step 2: Create a hash table for unique patterns
  pattern_hash <- new.env(hash = TRUE)

  # Step 3: Assign ranks to patterns
  ranks <- integer(n_rows)
  current_rank <- 1

  for (i in 1:n_rows) {
    pattern <- patterns[i]
    if (is.null(pattern_hash[[as.character(pattern)]])) {
      pattern_hash[[as.character(pattern)]] <- current_rank
      current_rank <- current_rank + 1
    }
    ranks[i] <- pattern_hash[[as.character(pattern)]]
  }

  return(ranks)
}

#---------------------------------------------------#
# Custom right/left bit-shift operator wrapper for
# base R bit-wise ops.

`%<<%` <- function(x,shift){
  bitwShiftL(x,shift)
}

`%>>%` <- function(x,shift){
  bitwShiftR(x,shift)
}

#----------------------------------------------------#
# Simple function to find column index of a subset of
# columns taken from a data frame. Used in `gen.mar`
# function
find.index <- function(data,cols){

  if(!is.data.frame(cols) && !is.matrix(cols)) {
    stop("The 'cols' argument should be a data frame or matrix.")
  }

  index <- which(colnames(data) %in% colnames(cols))

  return(index)
}

#-----------------------------------------------#
# Helper function for correlation output set to
# uppercase

upper_title <- function(s) {
  paste0(toupper(substring(s, 1L, 1L)), substring(s, 2L))
}

#-----------------------------------------------#
# Custom computation function for correlation
cor.estimate <- function(data,i,method = "kendall",use = "complete.obs"){

  stopifnot(is.data.frame(data) || is.matrix(data),
            is.character(use),
            is.character(method))

  method <- match.arg(method, c("pearson", "kendall", "spearman"))

  d <- data[i,]

  cor_val <- stats::cor(d[,1],d[,2],method = method, use = use)

  return(cor_val)
}
