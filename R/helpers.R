hotdeck.impute <- function(data,
                           method = "deterministic",
                           k = NULL,
                           na.rm = TRUE){

  stopifnot(is.matrix(data)|| is.data.frame(data))

  if(na.rm){
    nan_data <- data[complete.cases(data),]

    if(nrow(nan_data) == 0) {
      warning("No complete cases found.
              Using original data for imputation.")
      nan_data <- data
    }
  } else {
    nan_data <- data
  }

 for(j in 1:ncol(data)){
   na_index <- which(is.na(data[,j]))

   for(i in na_index){

     nan_col <- setdiff(1:ncol(data),j)

     dist <- apply(nan_data[,nan_col, drop = F],1,function(row){
       sqrt(sum((row - data[i, -j])^2,
                na.rm = TRUE))
     })

     donor_index <- switch(method,
       "deterministic" = which.min(dist),
       "rand_from_all" = sample(1:nrow(nan_data),1),
       "rand_nearest_k" = {
         k_nearest <- order(dist)[1:min(k,length(dist))]
         sample(k_nearest,1)
       },
       "weight_rand" = {
         weights <- 1/(dist + 1e-8)
         sample(1:nrow(nan_data),1, prob = weights)
       }
     )
     data[i,j] <- nan_data[donor_index,j]
   }
  }
 return(data)
}
