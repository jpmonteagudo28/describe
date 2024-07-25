#' Calculate the geometric median for non-collinear points in 2D+ Euclidean Space
#'
#' Using the Weiszfeld algorithm: an iterative method to compute the geometric
#' median of a set of points in Euclidean space. In general, measures of central
#' tendency minimize the sum of the Euclidean distance from the center to each
#' point in the set,this is represented as `arg min`. We try to come up with an
#' estimate of this center by first using our best guess and then approximating
#' it by calculating  the distance from each point to a newly chosen point in
#' the set and updating the center if the distance is smaller than that
#' calculated using our best, initial guess.
#'
#' You can find an explanation of a similar algorithm here:
#' https://www.sciencedirect.com/science/article/pii/S0898122109004234
#' NOTE: If you're sample is a 1D set of points use `median` or `grouped.median` instead.
#' The `geom.median` function is used for samples of points in 2D+ space.
#'
#' @param x: vector of non-collinear points
#'  @param iters: number of iterations used in Weiszfeld algorithm
#'  @param tol: tolerance value used to check convergence of changes between
#'               successive guesses.
#' @return a numeric vector of length 1L.
#' @export
#'
#' @examples
#' x <- matrix(rlnorm(60,1,1.4),ncol = 4)
#' geom.median(x,iters = 500, tol = 1e-10)
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

#' Calculate the median of distances to the geometric mean for a 2D+ space
#' We're calculating the median absolute deviation for matrices and arrays
#' @param x: a numeric vector or matrix containing a set of points in n-dimensions
#' @param iters: number of iterations used in Weiszfeld algorithm
#' @param na.rm: option to leave or remove NAs
#' @param ... optional arguments to specify alternative tolerance for convergence
#'              current tolerance set to 1e-8
#' @return a numeric vector of length 1L.
#' @export
#'
#' @examples
#' x <- matrix(rlnorm(60,1,1.4), ncol = 4)
#' geom.mad(x,iters = 500)
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

  inv_qnorm <- ((stats::qnorm(c(.25,.75),0,1))^-1)[2]
  k = 1/inv_qnorm

  if(!is.matrix(x)){
    return(stats::mad(x))
  } else{
    n <- length(x)
    center <- geom.median(x,iters = iters,...)

    dist <- sqrt(rowSums((x - center)^2))

    mad <- k * stats::median(dist)
    return(mad)
  }
}

#' Calculate the geometric mean for a set of observations.
#' The geometric mean is the nth root of the product of 'n' observations.
#' It multiplies rather than sums the observations in a set and
#' can be thought of as the exponential of the arithmetic mean of logarithms.
#' We use  this definition to create our function. #' The 'geom.mean'
#' function will calculate the mean for all numeric columns in  a matrix
#' or data frame.
#'
#' NOTE: Because we're using logarithms, our observations must be non-negative.
#' Additionally, by using the geometric mean and SD and variance, we'll be
#' measuring the log-normal dispersion of a log-normal distribution
#'
#' @param x: a numeric vector, matrix or data frame
#' @param ... additional parameters to be passed to 'mean()'
#' @param na.rm: indicate whether or not to remove NAs
#'
#' @return a numeric vector of length equal or greater than 1L.
#' @export
#'
#' @examples
#' x <- matrix(rlnorm(60,1,1.4), ncol = 4)
#' geom.mean(x)
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

#' Calculate the geometric standard deviation (GSD) of the geometric mean. the GSD
#' is a dimensionless (unitless) multiplicative factor. The GSD describes
#' the range from the mean divided by GSD to the mean multiplied by GSD. Because of
#' its multiplicative nature, the GSD cannot be added/subtracted from the mean.
#'
#' @param x: a numeric vector, matrix or data frame
#' @param ... additional parameters to be passed to 'geom.mean()'
#' @param na.rm: indicate whether or not to remove NAs
#'
#' @return a numeric vector of length equal or greater than 1L.
#' @export
#'
#' @examples
#' x <- matrix(rlnorm(60,1,1.4), ncol = 4)
#' geom.sd(x)

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

#' Calculate the geometric variance of a set of observations. The variance is also
#' a dimensionless factor.
#'
#' @param x: a numeric vector, matrix or data frame
#' @param ... additional parameters to be passed to 'geom.mean()'
#' @param na.rm: indicate whether or not to remove NAs
#'
#' @return a numeric vector of length equal or greater than 1L.
#' @export
#'
#' @examples
#' x <- matrix(rlnorm(60,1,1.4), ncol = 4)
#' geom.var(x)
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
#' Calculate the number of GSD by which the raw value is above/below the
#' mean. The formula here makes use of the change of base property of logarithms
#'
#' @param x: a numeric vector, matrix or data frame
#' @param ... additional parameters to be passed to 'geom.mean()'
#' @param na.rm: indicate whether or not to remove NAs
#'
#' @return a matrix or array of z-scores of length > 1L.
#' @export
#'
#' @examples
#' x <- matrix(rlnorm(60,1,1.4), ncol = 4)
#' geom.zscore(x)
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
