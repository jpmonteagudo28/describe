# Creates a publication-ready / formatted correlation matrix, using `Hmisc::rcorr`
## and ppcor::pcor.test in the back-end.
#
# @param df data-frame; containing numeric and/or logical columns to calculate
# correlations for
# @param type character; specifies the type of correlations to compute; gets passed
# to `Hmisc::rcorr`; options are `"pearson"` or `"spearman"`; defaults to `"pearson"`
# @param digits integer/double; number of decimals to show in the correlation matrix;
# gets passed to `formatC`; defaults to `3`
# @param decimal.mark character; which decimal.mark to use; gets passed to `formatC`;
# defaults to `.`
# @param use character; which part of the correlation matrix to display; options are
# `"all"`, `"upper"`, `"lower"`; defaults to `"all"`
# @param show_significance boolean; whether to add `*` to represent the significance
# levels for the correlations; defaults to `TRUE`
# @param replace_diagonal boolean; whether to replace the correlations on the
# diagonal; defaults to `FALSE`
# @param replacement character; what to replace the diagonal and/or upper/lower
# triangles with; defaults to `""` (empty string)

cor.matrix <- function(df,
                       type = "spearman",
                       digits = 3,
                       decimal.mark = ".",
                       use = "all",
                       show_significance = TRUE,
                       replace_diagonal = FALSE,
                       replacement = ""){

  # check arguments
  stopifnot(
    is.numeric(digits),
    digits >= 0,
    use %in% c("all", "upper", "lower"),
    is.logical(replace_diagonal),
    is.logical(show_significance),
    is.character(replacement)
  )

  # retain only numeric and boolean columns
  is_num_or_boolean = vapply(df, function(x) is.numeric(x) | is.logical(x),
                              logical(1))
  if (sum(!is_num_or_boolean) > 0) {
    cat('Dropping non-numeric/-boolean column(s):',
        paste(names(is_num_or_boolean)[!is_num_or_boolean], collapse = ', '),
        '\n\n')
  }
  df = df[is_num_or_boolean]

  # transform input data frame to matrix
  x <- as.matrix(df)

  # run correlation analysis using Hmisc package
  correlation_matrix <- Hmisc::rcorr(x, type = )
  R <- correlation_matrix[[r]] # Matrix of correlation coefficients
  p <- correlation_matrix[[P]] # Matrix of p-value

  # transform correlations to specific character format
  R_formatted = formatC(R, format = 'f', digits = digits,
                       decimal.mark = decimal.mark)

  # if negative numbers, we want to put a space before the positives to align all
  if (sum(R < 0) > 0) {
    R_formatted = ifelse(R > 0, paste0(' ', R_formatted), R_formatted)
  }

  # add significance levels if desired
  if (show_significance) {
    # define notions for significance levels; spacing is important.
    stars <- ifelse(is.na(p), "   ", ifelse(p < .05, "*  ", "   "))
    R_formatted = paste0(R_formatted, stars)
  }
  # build a new matrix with the formatted correlations and significance stars
  R_new <- matrix(R_formatted, ncol = ncol(x))
  rownames(R_new) <- colnames(x)
  colnames(R_new) <- paste(colnames(x), "", sep =" ")

  # replace undesired values
  if (use == 'upper') {
    Rnew[lower.tri(R_new, diag = replace_diagonal)] <- replacement
  } else if (use == 'lower') {
    Rnew[upper.tri(R_new, diag = replace_diagonal)] <- replacement
  } else if (replace_diagonal) {
    diag(Rnew) <- replacement
  }

  return(R_new)
}

## ----------------------------------------------------------------------------- ##
## Partial correlation function for data frames and tibbles
## ----------------------------------------------------------------------------- ##

part.cor.matrix <- function(df,z,
                     method = c("kendall","spearman","pearson"),
                     digits = 3,
                     decimal.mark = ".",
                     use = "all",
                     show_significance = TRUE,
                     replace_diagonal = FALSE,
                     replacement = ""){
  # check arguments
  stopifnot(
    is.numeric(digits),
    digits >= 0,
    use %in% c("all", "upper", "lower"),
    is.logical(replace_diagonal),
    is.logical(show_significance),
    is.character(replacement)
  )

  is_num_or_boolean = vapply(df, function(x) is.numeric(x) | is.logical(x),
                              logical(1))
  if (sum(!is_num_or_boolean) > 0) {
    cat('Dropping non-numeric/-boolean column(s):',
        paste(names(is_num_or_boolean)[!is_num_or_boolean], collapse = ', '),
        '\n\n')
  }
  df = df[is_num_or_boolean]

  if (!(z %in% colnames(df))) {
    stop("The specified variable to hold constant does
         not exist in the data frame.")
  }
  if (!is.numeric(df[[z]])) {
    stop("The specified variable to hold constant must be numeric.")
  }

  x <- df |> dplyr::select(-dplyr::all_of(z))
  z_mat <- df[[z]]


  n <- ncol(x)
  pairs <- utils::combn(n, 2, simplify = FALSE)


  pcor <- function(pair) {
    i <- pair[1]
    j <- pair[2]
    a <- x[[i]]
    b <- x[[j]]
    pcor_test <- ppcor::pcor.test(a, b, z_mat,
                                  method = match.arg(method,
                                                     c("pearson",
                                                       "kendall",
                                                       "spearman")))

    c(pcor_test[[estimate]], pcor_test[[p.value]])
  }

  results <- sapply(pairs,pcor)


  R <- matrix(1, ncol = n, nrow = n)
  p <- matrix(1, ncol = n, nrow = n)

  for (k in seq_along(pairs)) {
    i <- pairs[[k]][[1]]
    j <- pairs[[k]][[2]]
    R[i, j] <- results[1, k]
    R[j, i] <- results[1, k]
    p[i, j] <- results[2, k]
    p[j, i] <- results[2, k]
  }


  R_formatted = formatC(R, format = 'f', digits = digits,
                       decimal.mark = decimal.mark)

  if (sum(R < 0, na.rm = TRUE) > 0) {
    R_formatted = ifelse(R > 0, paste0(' ', R_formatted), R_formatted)
  }


  if (show_significance) {

    stars <- ifelse(is.na(p), "   ", ifelse(p < .05, "*  ", "   "))
    R_formatted = paste0(R_formatted, stars)
  }

  R_new <- matrix(R_formatted, ncol = ncol(x))
  rownames(Rnew) <- colnames(x)
  colnames(Rnew) <- paste(colnames(x), "", sep =" ")


  if (use == 'upper') {
    R_new[lower.tri(R_new, diag = replace_diagonal)] <- replacement
  } else if (use == 'lower') {
    R_new[upper.tri(R_new, diag = replace_diagonal)] <- replacement
  } else if (replace_diagonal) {
    diag(Rnew) <- replacement
  }

  return(R_new)
}

## ----------------------------------------------------------------------------- ##
## Bootstrap CI for correlation coefficient (kendall's tau)
## ----------------------------------------------------------------------------- ##
# Check this link for a Bootstrap refresher:
# http://users.stat.umn.edu/~helwig/notes/npboot-notes.html#need-for-nonparametric-bootstrap
# Use when sample distribution is unknown/not Gaussian
# x,y is your vector/matrix
# @N refers to the number of samples for the bootstrap procedure
# @ conf.level is the chosen confidence level (0 - 100)
# @alternative refers to less,greater or two-sided
# @ method refers to the calc. of correlation coefficient (pearson, spearman,
# kendall's)
# Fisher's transformation is used to calc. Pearson's r in case of small Ns and
# highly skewed data. Read here for info:
# https://blogs.sas.com/content/iml/2017/09/20/fishers-transformation-correlation.html





cor.boot.ci <- function(x,y,
                    N = NULL,
                    conf.level = 95,
                    alternative = "two.sided",
                    method = "kendall",
                    bootstrap = TRUE){

  if(length(x) != length(y))
    stop("'x' and 'y' must have the same length")

  if(!is.numeric(x)) stop("'x' must be a numeric vector")

  if(!is.numeric(y)) stop("'y' must be a numeric vector")


  stopifnot(
    is.numeric(N),
    is.numeric(conf.level),
    conf.level > 0 && conf.level < 100,
    is.character(alternative),
    is.logical(bootstrap)
  )

  alternative <- match.arg(alternative,c("two.sided", "less", "greater"))
  method <- match.arg(method, c("pearson", "kendall", "spearman"))

  if(bootstrap == TRUE && is.null(N))
    stop("Use method = 'kendall' or method = 'spearman' or provide the number of N to compute
         the confidence interval")

  if (!bootstrap && !is.null(N))
    warning("Samples are only used for bootstrap.If bootstrap = FALSE, N will be ignored.")

  if(bootstrap && is.null(N))
    stop("Set sample number to compute bootstrap estimate")

  if(bootstrap == TRUE && !is.null(N) && method == "kendall" || method == "spearman"){

    # Run correlation analysis
    cor.obs <- replicate(N, stats::cor(sample(x, replace = TRUE),
                                       sample(y, replace = TRUE),
                                       method = method))
    cor.val <- stats::cor(x,y, method = method)

    if (alternative == "less") {

      alpha <- conf.level / 100
      lower.ci <- -1
      upper.ci <- stats::quantile(cor.obs, probs = 1 - alpha)
    } else if (alternative == "greater") {

      alpha <- conf.level / 100
      lower.ci <- stats::quantile(cor.obs, probs = alpha)
      upper.ci <- 1
    } else if (alternative == "two.sided") {

      alpha <- (1 - conf.level / 100) / 2
      lower.ci <- stats::quantile(cor.obs, probs = alpha)
      upper.ci <- stats::quantile(cor.obs, probs = 1 - alpha)
    }
    result <- list(coefficient = cor.val,
                   lower.ci = lower.ci,
                   upper.ci = upper.ci)
  } else if( !bootstrap && method == "pearson") {

    cor.val <- stats::cor(x, y,method = method)

    # Fischer Transformation of Pearson's correlation coeff.
    z <- stats::qnorm(1 - (1 - conf.level/100) / 2)
    se <- 1 / sqrt(length(x) - 3)
    if (alternative == "less") {

      lower.ci <- -1
      upper.ci <- tanh(atanh(cor.val) + z * se)

    } else if (alternative == "greater") {

      lower.ci <- tanh(atanh(cor.val) - z * se)
      upper.ci <- 1

    } else if (alternative == "two.sided") {

      lower.ci <- tanh(atanh(cor.val) - z * se)
      upper.ci <- tanh(atanh(cor.val) + z * se)
    }
    result <- list(coefficient = cor.val,
                   lower.ci = lower.ci,
                   upper.ci = upper.ci)
  }
  else if(!bootstrap && method == "kendall" || method == "spearman")
    stop("Fisher's transformation applies to Pearson's correlation coefficient. Use bootstrap for other methods")

  dplyr::tibble(
    parameter = sprintf(uppercase_title(method)),
    coefficient = result$coefficient,
    lower.ci = result$lower.ci,
    upper.ci = result$upper.ci)
}
