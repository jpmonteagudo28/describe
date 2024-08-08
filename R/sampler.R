#> Gibbs sampler
#> if the target distribution has ind. components, random scan sampler
#> takes `log d` times as many updates to converge as the deterministic scan
#> as `iters` increases.
#>
sampler <- function(data,
                    iters,
                    burn_in = iters*.1,
                    ){

}


#' @references https://probability.ca/jeff/ftpdir/gibbsex.pdf
#' @references https://www.sciencedirect.com/science/article/pii/0304414994901341
#' @references https://vbn.aau.dk/files/257318283/ejbrm_volume15_issue1_article450.pdf
#'
