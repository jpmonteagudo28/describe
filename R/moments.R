#' @title Calculate skewness
#'
#' @description
#'
#' @param x a numeric vector of length +2L
#' @param ... additional arguments to be passed to function
#' @param na.rm option to remove or keep NAs, set to FALSE by default.
#'
#' @return a numeric vector of length 1
#' @export
#'
#' @examples
#' s <- rnorm(30,1,1.4)
#' skewness(s)
#'
#' air_data <- datasets::airquality
#' mapply(skewness,air_data)
#'
 skewness <- function(x,
                     ...,
                     na.rm = FALSE) {

  stopifnot(is.numeric(x),
            length(x) > 1L)

  if (is.data.frame(x) || is.matrix(x)) {

    stop("'x' must be a numeric vector")
  }

    if (na.rm) {
      x <- x[!is.na(x)]
    }

    n <- length(x)
    skew <- (sum((x - mean(x))^3) / (n * stats::sd(x)^3))

  return(skew)
}

#' @title Calculate kurtosis
#'
#' @description
#'
#' @param x a numeric vector of length +2L
#' @param ... additional arguments to be passed to function
#' @param na.rm option to remove or keep NAs, set to FALSE by default.
#'
#' @return a numeric vector of length 1
#' @export
#'
#' @examples
#' s <- rnorm(30,1,1.4)
#' kurtosis(s)
#'
#' air_data <- datasets::airquality
#' mapply(kurtosis,air_data)
kurtosis <- function(x,
                     ...,
                     na.rm = FALSE) {

  stopifnot(is.numeric(x))

  if (is.data.frame(x) || is.matrix(x)) {
    stop("'x' must be a numeric vector")
  }

    if (na.rm) {
      x <- x[!is.na(x)]
    }

    n <- length(x)
    kurt <- (sum((x - mean(x))^4) / (n * stats::sd(x)^4)) - 3

  return(kurt)
}


#' @title Calculate z-scores
#'
#' @description
#'
#' @param x a numeric vector of length +2L
#' @param ... additional arguments to be passed to function
#' @param na.rm option to remove or keep NAs, set to FALSE by default.
#'
#' @return a numeric vector of length 1
#' @export
#'
#' @examples
#' s <- rnorm(30,1,1.4)
#' z.scores(s)
#'
#' air_data <- datasets::airquality
#' mapply(z.scores,air_data)
#'
z.scores <- function(x,
                     ...,
                     na.rm = FALSE) {

  stopifnot(is.numeric(x))

  if (is.data.frame(x) || is.matrix(x)) {
    stop("'x' must be a numeric vector")
  }

    if (na.rm) {
      x <- x[!is.na(x)]
    }
    avg <- mean(x)
    stdev <- stats::sd(x)
    z_score <- (x - mean(x)) / stats::sd(x)

    return(z_score)
}

#' @title Calculate variance
#'
#' @description
#'
#' @param x a numeric vector of length +2L
#' @param ... additional arguments to be passed to function
#' @param na.rm option to remove or keep NAs, set to FALSE by default.
#'
#' @return a numeric vector of length 1
#' @export
#'
#' @examples
#' s <- rnorm(30,1,1.4)
#' variance(s)
#'
#' air_data <- datasets::airquality
#' mapply(variance,air_data)

variance <- function(x,
                     ...,
                     na.rm = FALSE) {

  stopifnot(is.numeric(x))

  if (is.data.frame(x) || is.matrix(x)) {
    stop("'x' must be a numeric vector")
  }
    if (na.rm) {
      x <- x[!is.na(x)]
    }

    n <- length(x)
    var <- sum((x - mean(x))^2) / (n - 1)

  return(var)
}
