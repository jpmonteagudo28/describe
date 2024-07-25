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
 #> @params x: vector of non-collinear points
#>  @params iters: number of iterations used in Weiszfeld algorithm
#>  @params tol: tolerance value used to check convergence of changes between
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

#> Calculate the geometric mean for a set of observations.
#> The geometric mean is the nth root of the product of 'n' observations.
#> It multiplies rather than sums the observations in a set and
#> can be thought of as the exponential of the arithmetic mean of logarithms.
#> We use  this definition to create our function.
#>
#> The 'geom.mean' function will calculate the mean for all numeric columns in
#> a matrix or data frame.
#>
#> @params x: a numeric vector, matrix or data frame
#> @params ... additional parameters to be passed to 'mean()'
#> @params na.rm: indicate whether or not to remove NAs
#>
#> NOTE: Because we're using logarithms, our observations must be non-negative.
#> Additionally, by using the geometric mean and SD, we'll be measuring the
#> log-normal dispersion of a log-normal distribution,
#>
geom.mean <- function(x, ..., na.rm = FALSE) {
  if (is.data.frame(x) || is.matrix(x)) {
    if (na.rm) {
      x <- apply(x, 2, function(column) column[!is.na(column)])
    }
    geom_mean <- apply(x, 2, function(column) {
      g.mean <- exp(mean(log(column)))
    })
    return(geom_mean)
  } else {
    if (na.rm) {
      x <- x[!is.na(x)]
    }
    g.mean <- exp(mean(log(x)))
    return(g.mean)
  }
}

#> Calculate the geometric standard deviation (GSD) of the geometric mean. the GSD
#> is a dimensionless (no units as output) multiplicative factor. The GSD describes
#> the range from the mean divided by GSD to the mean multiplied by GSD. Because of
#> its multiplicative nature, the GSD cannot be added/subtracted from the mean.
#>
#> @params x: a numeric vector, matrix or data frame
#> @params ... additional parameters to be passed to 'geom.mean()'
#> @params na.rm: indicate whether or not to remove NAs

geom.sd <- function(x,...,na.rm = FALSE){
  if (is.data.frame(x) || is.matrix(x)) {
    if (na.rm) {
      x <- apply(x, 2, function(column) column[!is.na(column)])
    }
    geom_sd <- apply(x, 2, function(column) {
      n <- length(column)
      g_mean <- geom.mean(column)
      g.sd <- exp(sqrt(sum((log(column/g_mean))^2) / n))
    })
    return(geom_sd)
  } else {
    if (na.rm) {
      x <- x[!is.na(x)]
    }
    n <- length(x)
    g_mean <- geom.mean(x)
    g.sd <- exp(sqrt(sum((log(x/g_mean))^2) / n))
    return(g.sd)
  }
}

#> Calculate the geometric variance of a set of observations. The variance is also
#> a dimensionless factor.
#>
#> @params x: a numeric vector, matrix or data frame
#> @params ... additional parameters to be passed to 'geom.mean()'
#> @params na.rm: indicate whether or not to remove NAs

geom.var <- function(x,...,na.rm = FALSE){
  if (is.data.frame(x) || is.matrix(x)) {
    if (na.rm) {
      x <- apply(x, 2, function(column) column[!is.na(column)])
    }
    geom_var <- apply(x, 2, function(column) {
      n <- length(column)
      g_mean <- geom.mean(column)
  g.var <- exp(sum((log(column/g_mean))^2) / n)
  })
    return(geom_var)
  }  else {
    if (na.rm) {
      x <- x[!is.na(x)]
    }
    n <- length(x)
    g_mean <- geom.mean(x)
    g.var <- exp(sum((log(x/g_mean))^2) / n)
    return(g.var)
  }
}
#> Calculate the number of GSD by which the raw value is above/below the
#> mean. The formula here makes use of the change of base property of logarithms
#>
#> @params x: a numeric vector, matrix or data frame
#> @params ... additional parameters to be passed to 'geom.mean()'
#> @params na.rm: indicate whether or not to remove NAs

geom.zscore <- function(x,...,na.rm = FALSE){
  if (is.data.frame(x) || is.matrix(x)) {
    if (na.rm) {
      x <- apply(x, 2, function(column) column[!is.na(column)])
    }
    geom_z <- apply(x, 2, function(column) {
      g_mean <- geom.mean(column)
      g_sd <- geom.sd(column)
      g.z <- log(column/g_mean, base = g_sd)
    })
    return(geom_z)
  }  else {
    if (na.rm) {
      x <- x[!is.na(x)]
    }
    g_mean <- geom.mean(x)
    g_sd<- geom.sd(x)
    g.z <- log(x/g_mean, base = g_sd)
    return(g.z)
  }
}
