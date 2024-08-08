#---------------------------------------#
# Convert character vector to factor
char2factor <- function(data) {

  if(!is.matrix(data) || !is.data.frame(data)){

    data <- as.factor(data)
  } else {
    char_vars <- which(sapply(data, is.character))

    for (var in char_vars) {
      data[[var]] <- as.factor(data[[var]])
    }
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
gen.predict <- function(data, varname, robust = FALSE, family = "AUTO") {

  complete_cases <- complete.cases(data)

  if(sum(complete_cases) == 0) {
    warning(paste("Imputation not performed. \nNo complete cases for variable:", varname))
    return(NULL)
  }

  # Check for sufficient degrees of freedom
  if (sum(complete_cases) <= ncol(data)) {
    warning(paste("Stopping missing value prediction. \nInsufficient degrees of freedom for variable:", varname))
    return(NULL)
  }

  # Check for one unique value
  if (length(unique(data[[varname]][complete_cases])) == 1) {
    warning(paste("Stopping missing value prediction. \nOnly one unique value for variable:", varname))
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
    warning(paste("Model fitting failed for variable:", varname, "Error:", e$message))
    return(NULL)
  })
  if(is.null(model)) return(NULL)

  predicted <- tryCatch({
    stats::predict(model, newdata = data,
                   type = if (model_type == "gaussian") "response" else "probs")
  }, error = function(e) {
    warning(paste("Prediction failed for variable:", varname, "\nError:", e$message))
    return(NULL)
  })

  return(predicted)
}


