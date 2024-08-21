#' @title
#' Calculate the harmonic mean
#'
#' @description
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
#' d <- rlnorm(60,1,1.4)
#' harm.mean(d)
#'
#' x <- matrix(rlnorm(60,1,1.4),ncol = 4)
#' mapply(harm.mean, as.data.frame(x))
harm.mean <- function(x, na.rm = FALSE,...) {

  if (is.data.frame(x) || is.matrix(x)) {
    stop("Provide a numeric vector instead of a matrix or data frame")
  }

  if(any(x < 0, na.rm = TRUE)){
    stop("Mean only defined for positive real numbers")
  }

    if (na.rm) {
      x <- x[!is.na(x)]
    }

    n <- length(x)
    inv.sum <- sum(x^-1)
    harm.mean <- n / inv.sum

    return(harm.mean)
}

#' @title
#' Calculate the weighted harmonic mean
#'
#' @param x an object containing a set of observations whose weighted harmonic mean
#' is to be computed
#' @param weights a numerical vector of weights the same length as x
#' @param na.rm logical by default set to FALSE. To remove NAs set argument to TRUE
#' @param ... additional arguments passed to or from other methods
#'
#' @return a numeric vector or length equal or greater than 1L
#' @export
#'
#' @examples
#' d <- rlnorm(60,1,1.4)
#' w <- rlnorm(60,2,3.23)
#' weight.harm.mean(d,w)
#'
#' x <- matrix(rlnorm(60,1,1.4), ncol = 4)
#' w <- matrix(rlnorm(60,1,1), ncol = 4)
#' mapply(weight.harm.mean, as.data.frame(x),as.data.frame(w))
weight.harm.mean <- function(x, weights = NULL, na.rm = FALSE, ...) {

  stopifnot(is.numeric(x),
            is.logical(na.rm)
            )

  if (!is.null(weights)) {
    stopifnot(is.numeric(weights))

    if (length(x) != length(weights)) {
      stop("'x' and 'weights' must be of the same length")
    }
  }

  if (na.rm) {
    valid <- !is.na(x)
    x <- x[valid]

    if (!is.null(weights)) {
      weights <- weights[valid]
    }
  }

  if (is.null(weights)) {
    harm_mean <- harm.mean(x)
  } else {
    inv_sum <- sum(weights / x)
    weight_sum <- sum(weights)
    harm_mean <- weight_sum / inv_sum
  }

  return(harm_mean)
}

#' @title Calculate the quadratic mean
#'
#' @description
#' Calculate the absolute magnitude for a numeric vector where greater values
#' are given more importance than lesser values. The arithmetic mean is typically
#' equal or greater than the arithmetic mean.
#'
#' @details
#' The quadratic mean or root mean square is scale invariant but not shift invariant, shifting
#' the input values by a fix amount will not change the point estimate.
#'
#' @param x a numeric vector of observations
#' @param ... additional arguments passed to or from other methods
#' @param na.rm logical by default set to FALSE. To remove NAs set argument to TRUE
#'
#' @return a numeric vector of length 1L
#' @export
#'
#' @examples
#' d <- rlnorm(60,1,1.4)
#' quad.mean(d)
#'
#' x <- matrix(rlnorm(60,1,1.4), ncol = 4)
#' mapply(quad.mean, as.data.frame(x))
#'
quad.mean <- function(x, ..., na.rm = FALSE) {

  if (is.data.frame(x) || is.matrix(x)) {
    stop("Provide a numeric vector instead of a matrix or data frame")
  }

  if(any(x < 0, na.rm = TRUE)){
    stop("Mean only defined for non-negative values")
  }
    if (na.rm) {
      x <- x[!is.na(x)]
    }
    n <- length(x)
    square.sum <- sum(x^2)
    sqrt.mean <- sqrt(square.sum / n)

    return(sqrt.mean)
}
