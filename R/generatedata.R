#' simulate data for simulations
#'
#' @description simulate data
#'
#' @param n sample size
#' @param mu mean of outcomes
#' @param psi1 coefficient for X for S(0)
#' @param psi2 coefficient for X for S(1)
#' @param omega1 coefficient for X for T(0)
#' @param omega2 coefficient for X for T(1)
#' @param sig variance parameter
#' @param R correlation matrix
#'
#' @return dataset
#'
#' @examples 
#' example(generatedata(n = 100, mu = c(0, 2, 3, 4), 
#' psi2 = 0, psi1 = 0, omega1 = 1, omega2 = 1, sig = 1))
     
generatedata = function(n, mu, psi2, psi1, omega1, omega2, sig, R){

S = diag(c(sig, sig, sig, sig))

if(missing(R)){
  R= matrix(rep(0,4*4),4,4)
  R[1,1] = 1; R[2,2] = 1; R[3,3] = 1; R[4,4] = 1; R[1,2] = 0.21428; 
  R[1,3] = 0.7; R[1,4] = 0.15; R[2,3] = 0.15
  R[2,4] = 0.7; R[3,4] = 0.15
  } else {
  R = R
  }



for(i in 2:4){
 for(j in 1:(i-1)){
 R[i,j] = R[j,i] }}

Sig = S%*%R%*%S

#X = rbinom(n,1, prob=.5)
X = rep(0, n)

samp = mvrnorm(n, mu, Sig)
samp[,1] = samp[,1] + psi1 %*% X
samp[,2] = samp[,2] + psi2 %*% X
samp[,3] = samp[,3] + omega1 %*% X
samp[,4] = samp[,4] + omega2 %*% X

trt= c(rep(0,n/2), rep(1, n/2))

ST= samp
# ST[1:(n/2),2] = 0
# ST[1:(n/2),4] = 0
# ST[(n/2 + 1):n,1] = 0
# ST[(n/2 + 1):n,3] = 0

ST = cbind(ST, X)

return(ST)

}
