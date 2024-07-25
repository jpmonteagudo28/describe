#> Calculate the geometric median for non-collinear points in 2D+ Euclidean Space
#>
#> Using the Weiszfeld algorithm: an iterative method to compute the geometric
#> median of a set of points in Euclidean space. In general, measures of central
#> tendency minimize the sum of the Euclidean distance from the center to each
#> point in the set,this is represented as `arg min`. We try to come up with an
#> estimate of this center by first using our best guess and then approximating
#> it by calculating  the distance from each point to a newly chosen point in
#> the set and updating the center if the distance is smaller than that
#> calculated using our best, initial guess.
#>
#> You can find an explanation of a similar algorithm here:
#> https://www.sciencedirect.com/science/article/pii/S0898122109004234

#>
#> #> params: @x: vector of non-collinear points
#>         @iters: number of iterations used in Weiszfeld algorithm
#>         @tol: tolerance value used to check convergence of changes between
#>               successive guesses.
#>
#> NOTE: If you're sample is a 1D set of points use `median` or `grouped.median` instead.
#> The `geom.median` function is used for samples of points in 2D+ space.

geom.median <- function(x,iters = 1000, tol = 1e-8,na.rm = FALSE){
  stopifnot(is.numeric(x),
            is.numeric(iters),
            is.numeric(tol),
            is.logical(na.rm)
  )

  if(na.rm){
    x <- x[!is.na(x)]
  }
  x <- as.matrix(x)
  center <- colMeans(x)

  for(i in 1:iters){
    dist <- sqrt(rowSums((x - center)^2))

    dist[dist == 0] <- tol

    weights <- 1/dist

    new_center <- colSums(weights * x)/sum(weights)

    if(sqrt(sum(new_center - center)^2)< tol){
      return(new_center)
    }
    center <- new_center
  }
  return(center)
}

#> Calculate the median of distances to the geometric mean for a 2D+ space
#> We're calculating the median absolute deviation for matrices and arrays
#> @params x: a numeric vector or matrix containing a set of points in n-dimensions
#> @params iters: number of iterations used in Weiszfeld algorithm
#> @params na.rm: option to leave or remove NAs
#> @params ... optional arguments to specify alternative tolerance for convergence
#>              current tolerance set to 1e-8

geom.mad <- function(x,
                     iters = 1000,
                     na.rm = FALSE,
                     ...){

  stopifnot(is.numeric(x),
            is.numeric(iters),
            is.logical(na.rm)
  )

  if(na.rm){
    x <- x[!is.na(x)]
  }

  inv_qnorm <- ((qnorm(c(.25,.75),0,1))^-1)[2]
  k = 1/inv.qnorm

  if(!is.matrix(x)){
    return(stats::mad(x))
  } else{
    n <- length(x)
    center <- geom.median(x,iters = iters,...)

    dist <- sqrt(rowSums((x - center)^2))

    mad <- k * median(dist)
    return(mad)
  }
}
