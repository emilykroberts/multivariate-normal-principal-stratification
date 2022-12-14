#' fdeltBeta
#'
#' @description likelihood for correlation parameter with a Beta prior
#'
#' @param n sample size
#' @param R correlation matrix
#' @param j inner component of normal likelihood
#'
#' @return proportional likelihood
#'
#' @examples 
#' example(fdeltBeta(n = 100, R = matrix(diag(c(1,1,1))), j = 10))
     
fdeltBeta = function(n, R, j){
 p = 2.7; q = 6 # want p>q
 a = -0.4 # lower bound of beta
 b = 1 # upper bound of beta
 return( -(n/2) * log(det(R)) + (-0.5 * j) + (q - 1) * log(b - R[2,3]) + # beta
  (p - 1) * log(R[2,3] - a) ) # alpha
}
