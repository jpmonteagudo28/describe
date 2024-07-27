#' Calculate mean for grouped data using frequency, either percent or counts.
#'
#' The function will split each character element in x and find the midpoint of
#' its lower and upper bounds.
#'
#' @param x character vector of categories or ranges
#' @param freq numeric vector containing category counts or percentages
#' @param freq_format intended freq format. Use percent if you want freq
#'        converted to percent
#' @param na.rm remove `NA` from calculations
#'
#' @return a numeric vector of length 1L
#' @export
#'
#' @examples
#' age <- c("0-4","5-9","10-14",
#'          "15-19","20-24",
#'          "25-29","30-34",
#'          "35-39","40-44",
#'          "45-49","50-54",
#'          "55-59","60-64",
#'          "65-69","70-74",
#'          "75-79","80-84",
#'          "85-89")
#'
#' freq <- seq(5400, 300, by = -300)
#' grouped.mean(age,freq, freq_format = "percent")
 grouped.mean <- function(x,
                          freq,
                          freq_format = "count",
                          na.rm = FALSE){
  stopifnot(is.character(x),
            is.numeric(freq),
            length(x) == length(freq),
            is.character(freq_format))

  freq_format<- match.arg(freq_format,c("percent","count"))

  midpoint <- function(x){
    bounds <- strsplit(x,"-")[[1]]

    lower <- as.numeric(bounds[1])
    upper <- as.numeric(bounds[2])

    if (is.na(lower) | is.na(upper)) {
      return(NA)
    } else {
      midpoint <- (lower + upper)/2
      return(midpoint)
    }
  }

  freq <- abs(freq)
  if(freq_format == "percent"){
    freq <- (freq/sum(freq)* 100)
  }

  x_mid <- sapply(x, midpoint)
  if(na.rm){
    index <- !is.na(x_mid)
    x_mid <- x_mid[index]
    freq <- freq[index]
  }

  n <- sum(freq)
  mu <- sum(x_mid*freq)/n

  return(mu)
}

#' Calculate median of grouped data using frequency, either percent or counts.
#'
#' @param x character vector of categories or ranges
#' @param freq numeric vector containing category counts or percentages
#' @param freq_format intended freq format. Use percent if you want freq
#'        converted to percent
#' @param na.rm remove `NA` from calculations
#'
#' @return a numeric vector of length 1L
#' @export
#'
#' @examples
#' age <- c("0-4","5-9","10-14",
#'          "15-19","20-24",
#'          "25-29","30-34",
#'          "35-39","40-44",
#'          "45-49","50-54",
#'          "55-59","60-64",
#'          "65-69","70-74",
#'          "75-79","80-84",
#'          "85-89")
#'
#' freq <- seq(5400, 300, by = -300)
#' grouped.median(age,freq)

grouped.median <- function(x,
                           freq,
                           freq_format = "count",
                           na.rm = FALSE){
  stopifnot(is.character(x),
            is.numeric(freq),
            length(x) == length(freq),
            is.character(freq_format),
            is.logical(na.rm))

  freq_format<- match.arg(freq_format,c("percent","count"))

  na_ind <- is.na(x) | is.na(freq)
  if (any(na_ind)) {
    if (na.rm) {
      x <- x[!na_ind]
      freq <- freq[!na_ind]
    } else {
      return(NA)
    }
  }
  if (length(x) == 0 || length(freq) == 0) {
    return(NA)
  }


  freq <- abs(freq)
  if (freq_format == "percent") {
    freq <- freq * sum(freq) / 100
  }


  df <- data.frame(x,freq,cf = cumsum(freq))
  half <- sum(freq)/2
  med_class <- df$x[df$cf >= half][1]
  med_class_freq <- df$freq[df$x == med_class]
  index <- which(df$cf >= half)[1]

  if(index > 1){
    prev_class_cf <- df$cf[index -1]
  } else{
    prev_class_cf <- 0
  }


  l <- function(x) as.numeric(strsplit(x,"-")[[1]][1])
  h <- function(x) as.numeric(strsplit(x,"-")[[1]][2]) -
    as.numeric(strsplit(x,"-")[[1]][1])

  lower <- l(med_class)
  size <- h(med_class)

  median <- lower + ((half  - prev_class_cf)/med_class_freq)*size
  return(median)
}

