#' Calculate the harmonic mean for matrices, data frames and vectors.
#'
#' The harmonic mean is the reciprocal of the arithmetic mean of the
#' reciprocals of a set of observations. It is a more robust measure of central
#' tendency than the arithmetic mean in the presence of outliers since it tends
#' to the smallest elements of a set of observations.
#'
#' @param x an object containing a set of observations in a matrix, data frame or vector.
#' @param na.rm logical by default set to FALSE. To remove NAs set argument to TRUE
#' @param ... additional arguments passed to or from other methods
#'
#' @return a numeric vector of length equal or greater than 1L
#' @export
#'
#' @examples
#' x <- matrix(rlnorm(60,1,1.4),ncol = 4)
#' harm.mean(x)
harm.mean <- function(x, na.rm = FALSE,...) {
  if (is.data.frame(x) || is.matrix(x)) {
    if (na.rm) {
      x <- apply(x, 2, function(column) na.omit(column))
    }
    harm_mean <- apply(x, 2, function(column) {
      n <- length(column)
      inv_sum <- sum(column^-1)
      harm.mean <- n / inv_sum
      return(harm.mean)
    })
    return(harm_mean)
  } else {
    if (na.rm) {
      x <- x[!is.na(x)]
    }
    n <- length(x)
    inv.sum <- sum(x^-1)
    harm.mean <- n / inv.sum
    return(harm.mean)
  }
}

#' Calculate the weighted harmonic mean for matrices, data frames and vectors
#'
#' @param x an object containing a set of observations whose weighted harmonic mean
#' is to be computed
#' @param w a numerical vector of weights the same length as x
#' @param na.rm logical by default set to FALSE. To remove NAs set argument to TRUE
#' @param ... additional arguments passed to or from other methods
#'
#' @return a numeric vector or length equal or greater than 1L
#' @export
#'
#' @examples
#' x <- matrix(rlnorm(60,1,1.4),ncol = 4)
#' w <- matrix(rnorm(60,0,0.23), ncol = 4)
#' w.harm.mean(x,w)
w.harm.mean <- function(x,w = NULL,na.rm = FALSE,...){


  if(is.null(w)){
    if(na.rm){
      x <- x[!is.na(x)]
      return(harm.mean(x))
    }
  }
  if(length(w) != length(x)){
    stop("'x' and 'w' must be of the same length")
  }
  w <- as.double(w)
  if (na.rm) {
    i <- !is.na(x)
    w <- w[i]
    x <- x[i]
  }
  inv_sum <- sum(w/x)
  w_sum <- sum(x)
  harm_mean <- w_sum/inv_sum
  return(harm_mean)
}


quad.mean <- function(x, ..., na.rm = FALSE) {
  if (is.data.frame(x) || is.matrix(x)) {
    if (na.rm) {
      x <- apply(x, 2, function(column) column[!is.na(column)])
    }

    quad_means <- apply(x, 2, function(column) {
      n <- length(column)
      square.sum <- sum(column^2)
      sqrt.mean <- sqrt(square.sum / n)
      return(sqrt.mean)
    })

    return(quad_means)
  } else {
    if (na.rm) {
      x <- x[!is.na(x)]
    }
    n <- length(x)
    square.sum <- sum(x^2)
    sqrt.mean <- sqrt(square.sum / n)
    return(sqrt.mean)
  }
}
