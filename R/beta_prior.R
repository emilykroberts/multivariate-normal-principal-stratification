#' beta_prior
#'
#' @description likelihood for correlation parameter with a Beta prior
#'
#' @param n sample size
#' @param R correlation matrix
#' @param j inner component of normal likelihood
#' @param coord1 row coordinate of correlation matrix
#' @param coord2 column coordinate of correlation matrix
#'
#' @return proportional likelihood
#'
#' @examples 
#' example(beta_prior(n = 100, R = matrix(diag(c(1, 1, 1))), j = 10))

beta_prior = function(n, R, j, coord1, coord2){
  p = 2.7; q = 6 # want p>q
  a = -0.4 # lower bound of beta
  b = 1 # upper bound of beta
  return( - (n / 2) * log(det(R)) + ( - 0.5 * j) + (q - 1) * log(b - R[coord1, coord2]) + # beta
            (p - 1) * log(R[coord1, coord2] - a)) # alpha
}
