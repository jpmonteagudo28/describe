hotdeck.impute <- function(data,
                           method = "deterministic",
                           k = NULL,
                           na.rm = TRUE){

  stopifnot(is.matrix(data)|| is.data.frame(data))

  if(na.rm){
    nna_data <- data[complete.cases(data),]

    if(nrow(nna_data) == 0) {
      warning("No complete cases found.
              Using original data for imputation.")
      nna_data <- data
    }
  } else {
    nna_data <- data
  }

 for(j in 1:ncol(data)){
   na_index <- which(is.na(data[,j]))

   for(i in na_index){

     nna_col <- setdiff(1:ncol(data),j)

     dist <- apply(nna_data[,nna_col, drop = F],1,function(row){
       sqrt(sum((row - data[i, -j])^2,
                na.rm = TRUE))
     })

     donor_index <- switch(method,
       "deterministic" = which.min(dist),
       "rand_from_all" = sample(1:nrow(nna_data),1),
       "rand_nearest_k" = {
         k_nearest <- order(dist)[1:min(k,length(dist))]
         sample(k_nearest,1)
       },
       "weight_rand" = {
         weights <- 1/(dist + 1e-8)
         sample(1:nrow(nna_data),1, prob = weights)
       }
     )
     data[i,j] <- nna_data[donor_index,j]
   }
  }
 return(data)
}

coldeck.impute <- function(x,
                           ext_data = NULL,
                           impute_method = "deterministic",
                           k = NULL,
                           na.rm = TRUE){
  if(is.null(ext_data)){
    stop("Either `ext_data` must be provided")
  }
  stopifnot(is.data.frame(data) || is.matrix(data),
            is.data.frame(ext_data) || is.matrix(ext_data),
            is.character(impute_method),
            is.logical(na.rm)
  )

}

multiple.imputation <- function(data,
                                imp_datasets = 5,
                                seed = NULL,
                                na.rm = TRUE){
  stopifnot(is.data.frame(data) || is.matrix(data),
            is.integer(imp_datasets),
            is.integer(seed),
            is.logical(na.rm))

    if(na.rm){
      nna_data <- data[complete.cases(data), ]

      if(nrow(nna_data) == 0) {
        warning("No complete cases found.
              Using original data for imputation.")
        nna_data <- data
      }
    } else {
      nna_data <- data
    }

  make.values <- function(col,
                          iters = imp_datasets){
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
  imp_data <- replicate(imp_datasets,nna_data,simplify = FALSE)

  for (j in 1:ncol(nna_data)) {
    if (any(is.na(nna_data[, j]))) {
      sim_vals <- make.values(nna_data[, j], iters)
      na_index <- which(is.na(nna_data[, j]))
      for (k in 1:iters) {
        imp_data[[k]][na_index, j] <- sim_vals[na_index, k]
      }
    }
  }
  return(imp_data)
}
