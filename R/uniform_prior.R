#' uniform_prior
#'
#' @description likelihood for correlation parameter with a uniform prior
#'
#' @param n sample size
#' @param R correlation matrix
#' @param j inner component of normal likelihood
#'
#' @return proportional likelihood
#'
#' @examples 
#' example(uniform_prior(n = 100, R = matrix(diag(c(1,1,1))), j = 10))

uniform_prior = function(n, R, j){ 
  return( - (n / 2) * log(det(R)) + ( - 0.5 * j))
}
