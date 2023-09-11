#' run the simulation
#'
#' importFrom("stats", "coef", "complete.cases", "cor", "lm", "median",
#'               "na.fail", "na.omit", "quantile", "runif", "sd", "var", "utils", "write.table")
#'
#' @description run mcmc simulation
#'
#' @param SIM number of iterations of mcmc to run
#' @param ST dataset
#' @param X covariate
#' @param n sample size
#' @param condindfit fit using conditional independence assumption
#' @param trt treatment allocation vector
#'
#' @return simulation results
#'
#' @examples 
#' example(run_sim(SIM = 1000, ST = STdat, X = X, n = 100))
run_sim = function(SIM, ST, X, n, condindfit, trt){
  
  fsig = function(n, s, j){ # put in own function if keeping
    return(- (n) * log(s) + (-0.5 * j)
    )
  }
  
  burnin = 0.3 * SIM
  # trt = c(rep(0, n / 2), rep(1, n / 2))
  
  {holdmu = matrix(rep(0, 4 * SIM), 4, SIM)
    holdmu1 = matrix(rep(0, 4 * SIM), 4, SIM)
    
    holdS = array(rep(0, 4 * 4 * SIM), dim = c(4, 4, SIM))
    holdR = array(rep(0.1, 4 * 4 * SIM), dim = c(4, 4, SIM))
    holdST = array(rep(0, 4 * n * SIM), dim = c(4, n, SIM))
    holdR[, , 1] = R = (diag(c(1, 1, 1, 1))); holdS[, , 1] = S = (diag(c(1, 1, 1, 1)))
    # identified parameters
    holdR[1, 3, 1] = holdR[3, 1, 1] = cor(ST[trt == 0, 1], ST[trt == 0, 3])
    holdR[2, 4, 1] = holdR[4, 2, 1] = cor(ST[trt == 1, 2], ST[trt == 1, 4])
    
    holdS[1, 1, 1] = sd(ST[trt == 0, 1])
    holdS[2, 2, 1] = sd(ST[trt == 1, 2])
    holdS[3, 3, 1] = sd(ST[trt == 0, 3])
    holdS[4, 4, 1] = sd(ST[trt == 1, 4])
    
    holdpsi1 = (rep(0, 1 * SIM));holdpsi2 = (rep(0, 1 * SIM))
    holdomega1 = (rep(0, 1 * SIM)); holdomega2 = (rep(0, 1 * SIM))
    holdalpha0 = (rep(0, 1 * SIM)); holdalpha01 = (rep(0, 1 * SIM))
    holdbeta0 = (rep(0, 1 * SIM)); holdbeta01 = (rep(0, 1 * SIM))
    
    slope = array(0, c((SIM), 1))
    int = array(0, c((SIM), 1))
    
    SIG0Smat = matrix(c(0.1, 0.1), nrow = 2, ncol = 2)
    theta0S = matrix(c(1, 0), nrow = 2)
    tauSS0 = tauSS1 = tauST0 = tauST1 = c(1)
    
    #prelim values
    #R[1, 2] = 0.3; R[1, 3] = 0.7; R[1, 4] = 0.15; R[2, 3] = 0.15
    #R[2, 4] = 0.7; R[3, 4] = 0.3
    #holdR[, , 1] = R
    
    for(i in 2:4){ for(j in 1:(i - 1)){ R[i, j] = R[j, i]}}
    XmatS = cbind(rep(1, n), X)
    sim = 2
    
    holdalpha0[1] = coef(lm(ST[trt == 0, 1] ~ X[trt == 0]))[1]
    holdalpha01[1] = coef(lm(ST[trt == 1, 2] ~ X[trt == 1]))[1]
    holdbeta0[1] = coef(lm(ST[trt == 0, 3] ~ X[trt == 0]))[1]
    holdbeta01[1] = coef(lm(ST[trt == 1, 4] ~ X[trt == 1]))[1]
    SIG0 = diag(rep(10 ^ 6, 4))
    
    ST = ST[, 1:4]
  }
  
  holdSmatrix = FALSE
    holdRmatrix_id = FALSE
    holdRmatrix_nonid = TRUE
    
  while(sim <= SIM){
    
    S = holdS[, , sim - 1]
    R = holdR[, , sim - 1]
    SIG = S %*% R %*% S
    
    mu = cbind(rep(holdalpha0[sim - 1], n), 
               rep(holdalpha01[sim - 1], n), 
               rep(holdbeta0[sim - 1], n), 
               rep(holdbeta01[sim - 1], n))
    
    for(i in 1:n){
      if(trt[i] == 0){
        ST[i, c(2, 4)] =
          c(mu[i, c(2, 4)] + (SIG[c(2, 4), c(1, 3)] %*% ginv(SIG[c(1, 3), c(1, 3)])) %*% (ST[i, c(1, 3)] - mu[i, c(1, 3)])) + 
          mvrnorm(1, c(0, 0), SIG[c(2, 4), c(2, 4)] - SIG[c(2, 4), c(1, 3)] %*% ginv(SIG[c(1, 3), c(1, 3)]) %*% SIG[c(1, 3), c(2, 4)])
      }
      if(trt[i] == 1){
        ST[i, c(1 , 3)] =
          c(mu[i, c(1 , 3)] + (SIG[c(1 , 3), c(2, 4)]  %*%  ginv(SIG[c(2, 4), c(2, 4)])) %*% (ST[i, c(2, 4)] - mu[i, c(2, 4)])) + 
          mvrnorm(1,  c(0, 0), SIG[c(1, 3), c(1, 3)] - SIG[c(1, 3), c(2, 4)] %*% ginv(SIG[c(2, 4), c(2, 4)]) %*% SIG[c(2, 4), c(1, 3)])
      }
    }
    
    #estimate coefficients
    Xmat = XmatS = cbind(rep(1, n), 0)
    
    mu = mvrnorm(1, Sigma = solve(n * solve(SIG) + solve(SIG0)), 
                 mu = solve(n * solve(SIG) + solve(SIG0)) %*% (n * solve(SIG) %*% (t(ST) %*% rep(1, n) / n)))
    
    holdalpha0[sim] = mu[1]
    holdalpha01[sim] = mu[2]
    holdbeta0[sim] = mu[3]
    holdbeta01[sim] = mu[4]
    
    tmp1 = (ST[, 1]) - holdalpha0[sim]
    
    tmp2 = (ST[, 2]) - holdalpha01[sim]

    tmp3 = (ST[, 3]) - holdbeta0[sim]

    tmp4 = (ST[, 4]) - holdbeta01[sim]
    
    #update entire sigma
    tmp = rbind(t(tmp1), t(tmp2), t(tmp3), t(tmp4))
    
    resid = ST - t(matrix((mu), 4, n, byrow = T))
    cor(resid)
    var(resid)
    
    a = b = 0.1
    
    if(F){
    s1 = MCMCpack::rinvgamma(1, shape = a  +  n / 2, scale = (sum(tmp1 ^ 2) / 2  +  b))
    s2 = MCMCpack::rinvgamma(1, shape = a  +  n / 2, scale = (sum(tmp2 ^ 2) / 2  +  b))
    s3 = MCMCpack::rinvgamma(1, shape = a  +  n / 2, scale = (sum(tmp3 ^ 2) / 2  +  b))
    s4 = MCMCpack::rinvgamma(1, shape = a  +  n / 2, scale = (sum(tmp4 ^ 2) / 2  +  b))

    S[1, 1] = s1
    S[2, 2] = s2
    S[3, 3] = s3
    S[4, 4] = s4
    
    # if(S[1,1] > (holdS[1,1,sim-1] + 2)) {print(S); next} # some settings are not converging well
    # if(S[2,2] > (holdS[2,2,sim-1] + 2)) {print(S); next} # some settings are not converging well
    # if(S[3,3] > (holdS[3,3,sim-1] + 2)) {print(S); next} # some settings are not converging well
    # if(S[4,4] > (holdS[4,4,sim-1] + 2)) {print(S); next} # some settings are not converging well
    }
    
    if(T){
      ####calc s's
      Rin = ginv(R)
      ##s1
      up = 250
      low = 100
      fs = matrix(rep(0, 201 * 4), 201, 4)
      for (k in low:up){
        s = k / 100
        Sig = S
        Sig[1, 1] = s
        fs[(k-low + 1), 1] = s
        
        j3 = sum(
          (resid[, 1] / Sig[1, 1]) * ((resid[, 1] / Sig[1, 1]) * Rin[1, 1] + (resid[, 2] / Sig[2, 2]) * Rin[2, 1] + (resid[, 3] / Sig[3, 3]) * Rin[3, 1] + (resid[, 4] / Sig[4, 4]) * Rin[4, 1]) + 
            (resid[, 2] / Sig[2, 2]) * ((resid[, 1] / Sig[1, 1]) * Rin[1, 2] + (resid[, 2] / Sig[2, 2]) * Rin[2, 2] + (resid[, 3] / Sig[3, 3]) * Rin[3, 2] + (resid[, 4] / Sig[4, 4]) * Rin[4, 2]) + 
            (resid[, 3] / Sig[3, 3]) * ((resid[, 1] / Sig[1, 1]) * Rin[1, 3] + (resid[, 2] / Sig[2, 2]) * Rin[2, 3] + (resid[, 3] / Sig[3, 3]) * Rin[3, 3] + (resid[, 4] / Sig[4, 4]) * Rin[4, 3]) + 
            (resid[, 4] / Sig[4, 4]) * ((resid[, 1] / Sig[1, 1]) * Rin[1, 4] + (resid[, 2] / Sig[2, 2]) * Rin[2, 4] + (resid[, 3] / Sig[3, 3]) * Rin[3, 4] + (resid[, 4] / Sig[4, 4]) * Rin[4, 4])
        )
        
        fs[(k-low + 1), 2] = fsig(n, s, j3)
      }
      fs = fs[!is.infinite(fs[, 2]), ] 
      fs[, 2] = fs[, 2] - median(fs[, 2])
      
      fs[, 2] = exp(fs[, 2])
      fs = fs[!is.infinite(fs[, 2]), ] 
      m = which(fs[, 2] == max(fs[, 2]))[1]
      x1 = max(1, (m-25))
      x2 = min(length(fs[, 1]), (m + 25))
      
      low = 1000 * fs[x1, 1]
      up = 1000 * fs[x2, 1]
      
      fs = array(0, c((up-low + 2), 4))
      
      for (k in low:up){
        s = k / 1000
        Sig = S
        Sig[1, 1] = s
        fs[(k-low + 1), 1] = s
        
        j3 = sum(
          (resid[, 1] / Sig[1, 1]) * ((resid[, 1] / Sig[1, 1]) * Rin[1, 1] + (resid[, 2] / Sig[2, 2]) * Rin[2, 1] + (resid[, 3] / Sig[3, 3]) * Rin[3, 1] + (resid[, 4] / Sig[4, 4]) * Rin[4, 1]) + 
            (resid[, 2] / Sig[2, 2]) * ((resid[, 1] / Sig[1, 1]) * Rin[1, 2] + (resid[, 2] / Sig[2, 2]) * Rin[2, 2] + (resid[, 3] / Sig[3, 3]) * Rin[3, 2] + (resid[, 4] / Sig[4, 4]) * Rin[4, 2]) + 
            (resid[, 3] / Sig[3, 3]) * ((resid[, 1] / Sig[1, 1]) * Rin[1, 3] + (resid[, 2] / Sig[2, 2]) * Rin[2, 3] + (resid[, 3] / Sig[3, 3]) * Rin[3, 3] + (resid[, 4] / Sig[4, 4]) * Rin[4, 3]) + 
            (resid[, 4] / Sig[4, 4]) * ((resid[, 1] / Sig[1, 1]) * Rin[1, 4] + (resid[, 2] / Sig[2, 2]) * Rin[2, 4] + (resid[, 3] / Sig[3, 3]) * Rin[3, 4] + (resid[, 4] / Sig[4, 4]) * Rin[4, 4])
        )
        
        fs[(k-low + 1), 2] = (fsig(n, s, j3))
      }
      fs = fs[c(fs[, 1] != 0), ]
      fs = fs[!is.infinite(fs[, 2]), ] 
      
      fs[, 2] = fs[, 2] - median(fs[, 2])
      
      fs[, 2] = exp(fs[, 2])
      fs = fs[!is.infinite(fs[, 2]), ] 
      for (k in 1: length(fs[, 1])){
        fs[k, 3] = fs[k, 2] / (sum(fs[, 2]))
      }
      
      for (k in 1:length(fs[, 1])){
        fs[k, 4] = sum(fs[1:k, 3])
      }
      
      u = runif(1, 0, 1)
      
      
      if (u<fs[1, 4]){
        s1 = fs[1, 1]
      }
      if (u>fs[length(fs[, 1]), 4]){
        s1 = fs[length(fs[, 1]), 1]
      }
      if (u >= fs[1, 4] & u <= fs[length(fs[, 1]), 4]){
        nr = nearest(fs[, 4], u)
        s1 = fs[nr, 1]
      }
      
      s1 = s1[1]
      
      ##s2
      up = 250
      low = 100
      fs = matrix(rep(0, 201 * 4), 201, 4)
      for (k in low:up){
        s = k / 100
        Sig = S
        Sig[1, 1] = s1
        Sig[2, 2] = s
        
        fs[(k-low + 1), 1] = s
        
        j3 = sum(
          (resid[, 1] / Sig[1, 1]) * ((resid[, 1] / Sig[1, 1]) * Rin[1, 1] + (resid[, 2] / Sig[2, 2]) * Rin[2, 1] + (resid[, 3] / Sig[3, 3]) * Rin[3, 1] + (resid[, 4] / Sig[4, 4]) * Rin[4, 1]) + 
            (resid[, 2] / Sig[2, 2]) * ((resid[, 1] / Sig[1, 1]) * Rin[1, 2] + (resid[, 2] / Sig[2, 2]) * Rin[2, 2] + (resid[, 3] / Sig[3, 3]) * Rin[3, 2] + (resid[, 4] / Sig[4, 4]) * Rin[4, 2]) + 
            (resid[, 3] / Sig[3, 3]) * ((resid[, 1] / Sig[1, 1]) * Rin[1, 3] + (resid[, 2] / Sig[2, 2]) * Rin[2, 3] + (resid[, 3] / Sig[3, 3]) * Rin[3, 3] + (resid[, 4] / Sig[4, 4]) * Rin[4, 3]) + 
            (resid[, 4] / Sig[4, 4]) * ((resid[, 1] / Sig[1, 1]) * Rin[1, 4] + (resid[, 2] / Sig[2, 2]) * Rin[2, 4] + (resid[, 3] / Sig[3, 3]) * Rin[3, 4] + (resid[, 4] / Sig[4, 4]) * Rin[4, 4])
        )
        
        fs[(k-low + 1), 2] = fsig(n, s, j3)
      }
      
      fs = fs[!is.infinite(fs[, 2]), ] 
      
      fs[, 2] = fs[, 2] - median(fs[, 2])
      fs[, 2] = exp(fs[, 2])
      fs = fs[!is.infinite(fs[, 2]), ] 
      m = which(fs[, 2] == max(fs[, 2]))[1]
      
      x1 = max(1, (m-25))
      x2 = min(length(fs[, 1]), (m + 25))
      
      low = 1000 * fs[x1, 1]
      up = 1000 * fs[x2, 1]
      fs = array(0, c((up-low + 2), 4))
      for (k in low:up){
        s = k / 1000
        Sig = S
        Sig[1, 1] = s1
        Sig[2, 2] = s
        
        fs[(k-low + 1), 1] = s
        
        j3 = sum(
          (resid[, 1] / Sig[1, 1]) * ((resid[, 1] / Sig[1, 1]) * Rin[1, 1] + (resid[, 2] / Sig[2, 2]) * Rin[2, 1] + (resid[, 3] / Sig[3, 3]) * Rin[3, 1] + (resid[, 4] / Sig[4, 4]) * Rin[4, 1]) + 
            (resid[, 2] / Sig[2, 2]) * ((resid[, 1] / Sig[1, 1]) * Rin[1, 2] + (resid[, 2] / Sig[2, 2]) * Rin[2, 2] + (resid[, 3] / Sig[3, 3]) * Rin[3, 2] + (resid[, 4] / Sig[4, 4]) * Rin[4, 2]) + 
            (resid[, 3] / Sig[3, 3]) * ((resid[, 1] / Sig[1, 1]) * Rin[1, 3] + (resid[, 2] / Sig[2, 2]) * Rin[2, 3] + (resid[, 3] / Sig[3, 3]) * Rin[3, 3] + (resid[, 4] / Sig[4, 4]) * Rin[4, 3]) + 
            (resid[, 4] / Sig[4, 4]) * ((resid[, 1] / Sig[1, 1]) * Rin[1, 4] + (resid[, 2] / Sig[2, 2]) * Rin[2, 4] + (resid[, 3] / Sig[3, 3]) * Rin[3, 4] + (resid[, 4] / Sig[4, 4]) * Rin[4, 4])
        )
      
        fs[(k-low + 1), 2] = fsig(n, s, j3)
      }
      fs = fs[c(fs[, 1] != 0), ]
      fs = fs[!is.infinite(fs[, 2]), ] 
      
      fs[, 2] = fs[, 2] - median(fs[, 2])
      fs[, 2] = exp(fs[, 2])
      fs = fs[!is.infinite(fs[, 2]), ] 
      
      for (k in 1: length(fs[, 1])){
        fs[k, 3] = fs[k, 2] / (sum(fs[, 2]))
      }
      
      for (k in 1:length(fs[, 1])){
        fs[k, 4] = sum(fs[1:k, 3])
      }
      
      u = runif(1, 0, 1)
      
      if (u<fs[1, 4]){
        s2 = fs[1, 1]
      }
      if (u>fs[length(fs[, 1]), 4]){
        s2 = fs[length(fs[, 1]), 1]
      }
      if (u >= fs[1, 4] & u <= fs[length(fs[, 1]), 4]){
        nr = nearest(fs[, 4], u)
        s2 = fs[nr, 1]
      }
      s2 = s2[1]
      
      ##s3
      up = 250
      low = 100
      fs = matrix(rep(0, 201 * 4), 201, 4)
      for (k in low:up){
        s = k / 100
        Sig = S
        Sig[1, 1] = s1
        Sig[2, 2] = s2
        Sig[3, 3] = s
        
        fs[(k-low + 1), 1] = s
        
        j3 = sum(
          (resid[, 1] / Sig[1, 1]) * ((resid[, 1] / Sig[1, 1]) * Rin[1, 1] + (resid[, 2] / Sig[2, 2]) * Rin[2, 1] + (resid[, 3] / Sig[3, 3]) * Rin[3, 1] + (resid[, 4] / Sig[4, 4]) * Rin[4, 1]) + 
            (resid[, 2] / Sig[2, 2]) * ((resid[, 1] / Sig[1, 1]) * Rin[1, 2] + (resid[, 2] / Sig[2, 2]) * Rin[2, 2] + (resid[, 3] / Sig[3, 3]) * Rin[3, 2] + (resid[, 4] / Sig[4, 4]) * Rin[4, 2]) + 
            (resid[, 3] / Sig[3, 3]) * ((resid[, 1] / Sig[1, 1]) * Rin[1, 3] + (resid[, 2] / Sig[2, 2]) * Rin[2, 3] + (resid[, 3] / Sig[3, 3]) * Rin[3, 3] + (resid[, 4] / Sig[4, 4]) * Rin[4, 3]) + 
            (resid[, 4] / Sig[4, 4]) * ((resid[, 1] / Sig[1, 1]) * Rin[1, 4] + (resid[, 2] / Sig[2, 2]) * Rin[2, 4] + (resid[, 3] / Sig[3, 3]) * Rin[3, 4] + (resid[, 4] / Sig[4, 4]) * Rin[4, 4])
        )
        
        fs[(k-low + 1), 2] = fsig(n, s, j3)
      }
      fs = fs[!is.infinite(fs[, 2]), ] 
      
      fs[, 2] = fs[, 2] - median(fs[, 2])
      fs[, 2] = exp(fs[, 2])
      fs = fs[!is.infinite(fs[, 2]), ] 
      m = which(fs[, 2] == max(fs[, 2]))[1]
      
      x1 = max(1, (m-25))
      x2 = min(length(fs[, 1]), (m + 25))
      
      low = 1000 * fs[x1, 1]
      up = 1000 * fs[x2, 1]
      
      fs = array(0, c((up-low + 2), 4))
      for (k in low:up){
        s = k / 1000
        Sig = S
        Sig[1, 1] = s1
        Sig[2, 2] = s2
        Sig[3, 3] = s
        
        fs[(k-low + 1), 1] = s
        
        j3 = sum(
          (resid[, 1] / Sig[1, 1]) * ((resid[, 1] / Sig[1, 1]) * Rin[1, 1] + (resid[, 2] / Sig[2, 2]) * Rin[2, 1] + (resid[, 3] / Sig[3, 3]) * Rin[3, 1] + (resid[, 4] / Sig[4, 4]) * Rin[4, 1]) + 
            (resid[, 2] / Sig[2, 2]) * ((resid[, 1] / Sig[1, 1]) * Rin[1, 2] + (resid[, 2] / Sig[2, 2]) * Rin[2, 2] + (resid[, 3] / Sig[3, 3]) * Rin[3, 2] + (resid[, 4] / Sig[4, 4]) * Rin[4, 2]) + 
            (resid[, 3] / Sig[3, 3]) * ((resid[, 1] / Sig[1, 1]) * Rin[1, 3] + (resid[, 2] / Sig[2, 2]) * Rin[2, 3] + (resid[, 3] / Sig[3, 3]) * Rin[3, 3] + (resid[, 4] / Sig[4, 4]) * Rin[4, 3]) + 
            (resid[, 4] / Sig[4, 4]) * ((resid[, 1] / Sig[1, 1]) * Rin[1, 4] + (resid[, 2] / Sig[2, 2]) * Rin[2, 4] + (resid[, 3] / Sig[3, 3]) * Rin[3, 4] + (resid[, 4] / Sig[4, 4]) * Rin[4, 4])
        )
        
        fs[(k-low + 1), 2] = fsig(n, s, j3)
      }
      fs = fs[c(fs[, 1] != 0), ]
      fs = fs[!is.infinite(fs[, 2]), ] 
      
      fs[, 2] = fs[, 2] - median(fs[, 2])
      fs[, 2] = exp(fs[, 2])
      fs = fs[!is.infinite(fs[, 2]), ] 
      for (k in 1: length(fs[, 1])){
        fs[k, 3] = fs[k, 2] / (sum(fs[, 2]))
      }
      
      for (k in 1:length(fs[, 1])){
        fs[k, 4] = sum(fs[1:k, 3])
      }
      
      u = runif(1, 0, 1)
      
      
      if (u<fs[1, 4]){
        s3 = fs[1, 1]
      }
      if (u>fs[length(fs[, 1]), 4]){
        s3 = fs[length(fs[, 1]), 1]
      }
      if (u >= fs[1, 4] & u <= fs[length(fs[, 1]), 4]){
        nr = nearest(fs[, 4], u)
        s3 = fs[nr, 1]
      }
      s3 = s3[1]
      
      ##s4
      up = 250
      low = 100
      fs = matrix(rep(0, 201 * 4), 201, 4)
      for (k in low:up){
        s = k / 100
        Sig = S
        Sig[1, 1] = s1
        Sig[2, 2] = s2
        Sig[3, 3] = s3
        Sig[4, 4] = s
        
        fs[(k-low + 1), 1] = s
        
        j3 = sum(
          (resid[, 1] / Sig[1, 1]) * ((resid[, 1] / Sig[1, 1]) * Rin[1, 1] + (resid[, 2] / Sig[2, 2]) * Rin[2, 1] + (resid[, 3] / Sig[3, 3]) * Rin[3, 1] + (resid[, 4] / Sig[4, 4]) * Rin[4, 1]) + 
            (resid[, 2] / Sig[2, 2]) * ((resid[, 1] / Sig[1, 1]) * Rin[1, 2] + (resid[, 2] / Sig[2, 2]) * Rin[2, 2] + (resid[, 3] / Sig[3, 3]) * Rin[3, 2] + (resid[, 4] / Sig[4, 4]) * Rin[4, 2]) + 
            (resid[, 3] / Sig[3, 3]) * ((resid[, 1] / Sig[1, 1]) * Rin[1, 3] + (resid[, 2] / Sig[2, 2]) * Rin[2, 3] + (resid[, 3] / Sig[3, 3]) * Rin[3, 3] + (resid[, 4] / Sig[4, 4]) * Rin[4, 3]) + 
            (resid[, 4] / Sig[4, 4]) * ((resid[, 1] / Sig[1, 1]) * Rin[1, 4] + (resid[, 2] / Sig[2, 2]) * Rin[2, 4] + (resid[, 3] / Sig[3, 3]) * Rin[3, 4] + (resid[, 4] / Sig[4, 4]) * Rin[4, 4])
        )
        
        fs[(k-low + 1), 2] = fsig(n, s, j3)
      }
      fs = fs[!is.infinite(fs[, 2]), ] 
      
      fs[, 2] = fs[, 2] - median(fs[, 2])
      fs[, 2] = exp(fs[, 2])
      fs = fs[!is.infinite(fs[, 2]), ] 
      
      m = which(fs[, 2] == max(fs[, 2]))[1]
      
      x1 = max(1, (m-25))
      x2 = min(length(fs[, 1]), (m + 25))
      
      low = 1000 * fs[x1, 1]
      up = 1000 * fs[x2, 1]
      
      fs = array(0, c((up-low + 2), 4))
      for (k in low:up){
        s = k / 1000
        Sig = S
        Sig[1, 1] = s1
        Sig[2, 2] = s2
        Sig[3, 3] = s3
        Sig[4, 4] = s
        
        fs[(k-low + 1), 1] = s
        
        j3 = sum(
          (resid[, 1] / Sig[1, 1]) * ((resid[, 1] / Sig[1, 1]) * Rin[1, 1] + (resid[, 2] / Sig[2, 2]) * Rin[2, 1] + (resid[, 3] / Sig[3, 3]) * Rin[3, 1] + (resid[, 4] / Sig[4, 4]) * Rin[4, 1]) + 
            (resid[, 2] / Sig[2, 2]) * ((resid[, 1] / Sig[1, 1]) * Rin[1, 2] + (resid[, 2] / Sig[2, 2]) * Rin[2, 2] + (resid[, 3] / Sig[3, 3]) * Rin[3, 2] + (resid[, 4] / Sig[4, 4]) * Rin[4, 2]) + 
            (resid[, 3] / Sig[3, 3]) * ((resid[, 1] / Sig[1, 1]) * Rin[1, 3] + (resid[, 2] / Sig[2, 2]) * Rin[2, 3] + (resid[, 3] / Sig[3, 3]) * Rin[3, 3] + (resid[, 4] / Sig[4, 4]) * Rin[4, 3]) + 
            (resid[, 4] / Sig[4, 4]) * ((resid[, 1] / Sig[1, 1]) * Rin[1, 4] + (resid[, 2] / Sig[2, 2]) * Rin[2, 4] + (resid[, 3] / Sig[3, 3]) * Rin[3, 4] + (resid[, 4] / Sig[4, 4]) * Rin[4, 4])
        )
        
        fs[(k-low + 1), 2] = fsig(n, s, j3)
      }
      fs = fs[c(fs[, 1] != 0), ]
      fs = fs[!is.infinite(fs[, 2]), ] 
      
      fs[, 2] = fs[, 2] - median(fs[, 2])
      fs[, 2] = exp(fs[, 2])
      fs = fs[!is.infinite(fs[, 2]), ] 
      for (k in 1: length(fs[, 1])){
        fs[k, 3] = fs[k, 2] / (sum(fs[, 2]))
      }
      
      for (k in 1:length(fs[, 1])){
        fs[k, 4] = sum(fs[1:k, 3])
      }
      
      u = runif(1, 0, 1)
      
      if (u<fs[1, 4]){
        s4 = fs[1, 1]
      }
      if (u>fs[length(fs[, 1]), 4]){
        s4 = fs[length(fs[, 1]), 1]
      }
      if (u >= fs[1, 4] & u <= fs[length(fs[, 1]), 4]){
        nr = nearest(fs[, 4], u)
        s4 = fs[nr, 1]
      }
      s4 = s4[1]
      
      s1 = sd(tmp1)
      s2 = sd(tmp2)
      s3 = sd(tmp3)
      s4 = sd(tmp4)
      
      S = diag(c(s1, s2, s3, s4))
      
    }
    
    if(holdSmatrix) S = holdS[, , 1]
    
    if(!holdRmatrix_nonid){
      ##r12
      a12 = (R[3, 4] ^ 2 - 1)
      b12 = 2 * R[1, 4] * R[2, 4] - 2 * R[1, 3] * R[2, 4] * R[3, 4] - 2 * R[1, 4] * R[2, 3] * R[3, 4] + 2 * R[1, 3] * R[2, 3]
      c12 = 1 - R[3, 4] ^ 2 - R[2, 3] ^ 2 - R[2, 4] ^ 2 - R[1, 3] ^ 2 - R[1, 4] ^ 2 + R[1, 3] ^ 2 * R[2, 4] ^ 2 + R[1, 4] ^ 2 * R[2, 3] ^ 2 + 2 * R[2, 3] * R[2, 4] * R[3, 4] - 2 * R[1, 3] * R[1, 4] * R[2, 3] * R[2, 4] + 2 * R[1, 3] * R[1, 4] * R[3, 4]
      
      L12 = (- b12 + sqrt((b12 ^ 2) - 4 * a12 * c12)) / (2 * a12)
      U12 = (- b12 - sqrt((b12 ^ 2) - 4 * a12 * c12)) / (2 * a12)
      
      low12 = ceiling(100 * max(0, min(L12, U12)))
      up12 = floor(100 * max(min(1, L12), min(1, U12)))
      
      d12 = up12 - low12 + 1
      
      fr12 = matrix(rep(0, d12 * 4), d12, 4)
      for (k in low12:up12){
        r = k / 100
        Rho = R
        Rho[1, 2] = r
        Rho[2, 1] = r
        fr12[(k - low12 + 1), 1] = r
        
        summand = apply(resid, 1, function(resid) resid %*% ginv(S %*% Rho %*% S) %*% resid)
        
        fr12[(k - low12 + 1), 2] = uniform_prior(n = n, R = Rho, j = sum(summand))
      }
      
      fr12 = matrix(fr12, ncol = 4)
      fr12 = fr12[!is.na(fr12[, 2]), ] 
      fr12 = matrix(fr12, ncol = 4)
      fr12 = fr12[!is.infinite(fr12[, 2]), ] 
      fr12 = matrix(fr12, ncol = 4)
      fr12 = fr12[complete.cases(fr12), ]
      fr12 = matrix(fr12, ncol = 4)
      fr12[, 2] = fr12[, 2] - median(fr12[, 2])
      fr12 = matrix(fr12, ncol = 4)
      fr12[, 2] = exp(fr12[, 2])
      fr12 = matrix(fr12, ncol = 4)
      fr12 = fr12[!is.infinite(fr12[, 2]), ] 
      fr12 = matrix(fr12, ncol = 4)
      
      m = which(fr12[, 2] == max(fr12[, 2]))[1]
      x1 = max(1, (m - 10))
      x2 = min(length(fr12[, 1]), (m + 10))
      low12 = 500 * fr12[x1, 1]
      up12 = 500 * fr12[x2, 1]
      d12 = up12 - low12 + 1
      
      fr12 = matrix(rep(0, d12 * 4), d12, 4)
      for (k in low12:up12){
        r = k / 500
        Rho = R
        Rho[1, 2] = r
        Rho[2, 1] = r
        fr12[(k - low12 + 1), 1] = r
        
        summand = apply(resid, 1, function(resid) resid %*% ginv(S %*% Rho %*% S) %*% resid)
        
        fr12[(k - low12 + 1), 2] = uniform_prior(n = n, R = Rho, j = sum(summand))
      }
      fr12 = matrix(fr12, ncol = 4)
      
      fr12 = matrix(fr12, ncol = 4)
      fr12 = fr12[!is.na(fr12[, 2]), ] 
      fr12 = matrix(fr12, ncol = 4)
      fr12 = fr12[!is.infinite(fr12[, 2]), ] 
      fr12 = matrix(fr12, ncol = 4)
      fr12 = fr12[complete.cases(fr12), ]
      fr12 = matrix(fr12, ncol = 4)
      fr12[, 2] = fr12[, 2] - median(fr12[, 2])
      fr12 = matrix(fr12, ncol = 4)
      fr12[, 2] = exp(fr12[, 2])
      fr12 = matrix(fr12, ncol = 4)
      fr12 = fr12[!is.infinite(fr12[, 2]), ] 
      fr12 = matrix(fr12, ncol = 4)
      
      for (k in 1: length(fr12[, 1])){
        fr12[k, 3] = fr12[k, 2] / (sum(fr12[, 2]))
      }
      
      for (k in 1:length(fr12[, 1])){
        fr12[k, 4] = sum(fr12[1:k, 3])
      }
      
      u = runif(1, 0, 1)
      
      if (u<fr12[1, 4]){
        r12 = fr12[1, 1]
      }
      if (u > fr12[length(fr12[, 1]), 4]){
        r12 = fr12[length(fr12[, 1]), 1]
      }
      if (u >= fr12[1, 4] & u <= fr12[length(fr12[, 1]), 4]){
        nr = nearest(fr12[, 4], u)
        r12 = fr12[nr, 1]
      }
      r12[r12>1] = 1
      r12 = r12[1]
      
      R[1, 2] = r12
      R[2, 1] = r12
    }
    
    if(!holdRmatrix_id){
      
      #r13
      a13 = (R[2, 4] ^ 2 - 1)
      b13 = 2 * R[1, 2] * R[2, 3] + 2 * R[1, 4] * R[3, 4] - 2 * R[1, 2] * R[2, 4] * R[3, 4] - 2 * R[1, 4] * R[2, 3] * R[2, 4]
      c13 = 1 - R[3, 4] ^ 2 - R[2, 3] ^ 2 - R[2, 4] ^ 2 - R[1, 2] ^ 2 - R[1, 4] ^ 2 + R[1, 2] ^ 2 * R[3, 4] ^ 2 + R[1, 4] ^ 2 * R[2, 3] ^ 2 + 2 * R[2, 3] * R[2, 4] * R[3, 4] + 2 * R[1, 2] * R[1, 4] * R[2, 4] - 2 * R[1, 2] * R[2, 3] * R[3, 4] * R[1, 4]
      
      L13 = (- b13 + sqrt((b13 ^ 2) - 4 * a13 * c13)) / (2 * a13)
      U13 = (- b13 - sqrt((b13 ^ 2) - 4 * a13 * c13)) / (2 * a13)
      
      a13 = - R[3, 4] * R[2, 4] * R[1, 4] - R[2, 4] ^ 2 * R[2, 3] * R[2, 4] ^ 2 * R[1, 4] - R[2, 4] * R[1, 2] + R[1, 2] * R[3, 4]
      b13 = - R[2, 4] * R[1, 4] * R[2, 3] * R[1, 2] * R[3, 4] - R[2, 4] ^ 2 * R[1, 4] * R[1, 2] - R[3, 4] * R[2, 4] ^ 2 * R[2, 3] + 
        R[1, 2] * R[3, 4] * R[1, 4] * R[2, 3] * R[2, 4] + R[3, 4] ^ 2 * R[2, 4] - R[1, 2] * R[3, 4] ^ 2 * R[1, 4]
      c13 = - (R[1, 4] + R[1, 4] * R[2, 3] ^ 2) * (R[2, 3] + R[2, 3] * R[1, 4] ^ 2) * 
        (R[3, 4] * R[1, 2] * R[2, 4] ^ 2 - R[1, 2] ^ 2 * R[2, 4] * R[1, 4] * R[2, 4] - R[1, 4] ^ 2 * R[2, 4] * R[2, 3] * R[1, 2] + R[1, 2] ^ 2 * R[3, 4] ^ 2 * R[2, 3] * R[1, 4])
      
      low13 = ceiling(100 * max(0, min(L13, U13)))
      up13 = floor(100 * max(L13, U13))
      d13 = up13 - low13 + 1
      
      fr13 = matrix(rep(0, d13 * 4), d13, 4)
      for (k in low13:up13){
        r = k / 100
        Rho = R
        Rho[1, 3] = r
        Rho[3, 1] = r
        fr13[(k - low13 + 1), 1] = r
        
        summand = apply(resid, 1, function(resid) resid %*% ginv(S %*% Rho %*% S) %*% resid)
        
        fr13[(k - low13 + 1), 2] = uniform_prior(n = n, R = Rho, j = sum(summand))
      }
      
      fr13 = fr13[!is.na(fr13[, 2]), ]
      fr13 = fr13[!is.infinite(fr13[, 2]), ] 
      fr13 = fr13[(fr13[, 2]!= 0), ] 
      
      fr13 = fr13[complete.cases(fr13), ]
      fr13[, 2] = fr13[, 2] - median(fr13[, 2])
      fr13[, 2] = exp(fr13[, 2])
      fr13 = fr13[!is.infinite(fr13[, 2]), ] 
      m = which(fr13[, 2] == max(fr13[, 2]))[1]
      x1 = max(1, (m - 10))
      x2 = min(length(fr13[, 1]), (m + 10))
      fr13[m, ]
      
      low13 = 500 * fr13[x1, 1]
      up13 = 500 * fr13[x2, 1]
      
      d13 = up13 - low13 + 1
      
      fr13 = matrix(rep(0, d13 * 4), d13, 4)
      for (k in low13:up13){
        r = k / 500
        Rho = R
        Rho[1, 3] = r
        Rho[3, 1] = r
        fr13[(k - low13 + 1), 1] = r
        summand = apply(resid, 1, function(resid) resid %*% ginv(S %*% Rho %*% S) %*% resid)
        
        fr13[(k - low13 + 1), 2] = uniform_prior(n = n, R = Rho, j = sum(summand))
      }
      
      fr13 = fr13[!is.na(fr13[, 2]), ]
      fr13 = fr13[!is.infinite(fr13[, 2]), ] 
      fr13 = fr13[complete.cases(fr13), ]
      fr13[, 2] = fr13[, 2] - median(fr13[, 2])
      fr13[, 2] = exp(fr13[, 2])
      fr13 = fr13[!is.infinite(fr13[, 2]), ] 
      for (k in 1: length(fr13[, 1])){
        fr13[k, 3] = fr13[k, 2] / (sum(fr13[, 2]))
      }
      
      for (k in 1:length(fr13[, 1])){
        fr13[k, 4] = sum(fr13[1:k, 3])
      }
      
      
      u = runif(1, 0, 1)
      if (u<fr13[1, 4]){
        r13 = fr13[1, 1]
      }
      if (u > fr13[length(fr13[, 1]), 4]){
        r13 = fr13[length(fr13[, 1]), 1]
      }
      if (u >= fr13[1, 4] & u <= fr13[length(fr13[, 1]), 4]){
        nr = nearest(fr13[, 4], u)
        r13 = fr13[nr, 1]
      }
      
      r13 = r13[1]
      R[1, 3] = r13
      R[3, 1] = r13
      
      #r24
      
      a24 = (R[1, 3] ^ 2 - 1)
      b24 = 2 * R[2, 3] * R[3, 4] + 2 * R[1, 2] * R[1, 4] - 2 * R[1, 2] * R[1, 3] * R[3, 4] - 2 * R[1, 3] * R[1, 4] * R[2, 3]
      c24 = 1 - R[3, 4] ^ 2 - R[2, 3] ^ 2 - R[1, 2] ^ 2 - R[1, 3] ^ 2 - R[1, 4] ^ 2 + R[1, 2] ^ 2 * R[3, 4] ^ 2 + R[1, 4] ^ 2 * R[2, 3] ^ 2 + 2 * R[1, 2] * R[1, 3] * R[2, 3] + 2 * R[1, 3] * R[1, 4] * R[3, 4] - 2 * R[1, 2] * R[2, 3] * R[3, 4] * R[1, 4]
      L24 = (- b24 + sqrt((b24 ^ 2) - 4 * a24 * c24)) / (2 * a24)
      U24 = (- b24 - sqrt((b24 ^ 2) - 4 * a24 * c24)) / (2 * a24)
      
      #low24 = ceiling(100 * max(0, R[2, 3], R[1, 4], min(L24, U24)))
      low24 = ceiling(100 * max(0, min(L24, U24)))
      up24 = floor(100 * max(min(1, L24), min(1, U24)))
      d24 = up24 - low24 + 1
      fr24 = matrix(rep(0, d24 * 4), d24, 4)
      for (k in low24:up24){
        r = k / 100
        Rho = R
        Rho[2, 4] = r
        Rho[4, 2] = r
        fr24[(k - low24 + 1), 1] = r
        summand = apply(resid, 1, function(resid) resid %*% ginv(S %*% Rho %*% S) %*% resid)
        
        fr24[(k - low24 + 1), 2] = uniform_prior(n = n, R = Rho, j = sum(summand))
      }
      fr24 = matrix(fr24, ncol = 4)
      
      fr24 = fr24[!is.na(fr24[, 2]), ] 
      fr24 = matrix(fr24, ncol = 4)
      fr24 = fr24[!is.infinite(fr24[, 2]), ]
      fr24 = matrix(fr24, ncol = 4)
      fr24 = fr24[complete.cases(fr24), ] 
      fr24 = matrix(fr24, ncol = 4)
      fr24[, 2] = fr24[, 2] - median(fr24[, 2])
      fr24 = matrix(fr24, ncol = 4)
      fr24[, 2] = exp(fr24[, 2])
      fr24 = matrix(fr24, ncol = 4)
      fr24 = fr24[!is.infinite(fr24[, 2]), ] 
      fr24 = matrix(fr24, ncol = 4)
      
      m = which(fr24[, 2] == max(fr24[, 2]))[1]
      x1 = max(1, (m - 10))
      x2 = min(length(fr24[, 1]), (m + 10))
      low24 = 500 * fr24[x1, 1]
      up24 = 500 * fr24[x2, 1]
      fr24[m, ]
      
      d24 = up24 - low24 + 1
      
      fr24 = matrix(rep(0, d24 * 4), d24, 4)
      for (k in low24:up24){
        r = k / 500
        Rho = R
        Rho[2, 4] = r
        Rho[4, 2] = r
        fr24[(k - low24 + 1), 1] = r
        summand = apply(resid, 1, function(resid) resid %*% ginv(S %*% Rho %*% S) %*% resid)
        
        fr24[(k - low24 + 1), 2] = uniform_prior(n = n, R = Rho, j = sum(summand))
      }
      fr24 = matrix(fr24, ncol = 4)
      fr24 = fr24[!is.na(fr24[, 2]), ] 
      fr24 = matrix(fr24, ncol = 4)
      fr24 = fr24[!is.infinite(fr24[, 2]), ]
      fr24 = matrix(fr24, ncol = 4)
      fr24 = fr24[complete.cases(fr24), ] 
      fr24 = matrix(fr24, ncol = 4)
      fr24[, 2] = fr24[, 2] - median(fr24[, 2])
      fr24 = matrix(fr24, ncol = 4)
      fr24[, 2] = exp(fr24[, 2])
      fr24 = matrix(fr24, ncol = 4)
      fr24 = fr24[!is.infinite(fr24[, 2]), ] 
      fr24 = matrix(fr24, ncol = 4)
      
      for (k in 1: length(fr24[, 1])){
        fr24[k, 3] = fr24[k, 2] / (sum(fr24[, 2]))
      }
      
      for (k in 1:length(fr24[, 1])){
        fr24[k, 4] = sum(fr24[1:k, 3])
      }
      
      u = runif(1, 0, 1)
      
      if (u<fr24[1, 4]){
        r24 = fr24[1, 1]
      }
      if (u > fr24[length(fr24[, 1]), 4]){
        r24 = fr24[length(fr24[, 1]), 1]
      }
      if (u >= fr24[1, 4] & u <= fr24[length(fr24[, 1]), 4]){
        nr = nearest(fr24[, 4], u)
        r24 = fr24[nr, 1]
      }
      r24 = r24[1]
      
      R[2, 4] = r24
      R[4, 2] = r24
    }
    
    if(!holdRmatrix_nonid){
      #r14
      a14 = (R[2, 3] ^ 2 - 1)
      b14 = 2 * R[1, 2] * R[2, 4] + 2 * R[1, 3] * R[3, 4] - 2 * R[1, 2] * R[2, 3] * R[3, 4] - 2 * R[1, 3] * R[2, 3] * R[2, 4]
      c14 = 1 - R[1, 2] ^ 2 - R[1, 3] ^ 2 - R[2, 3] ^ 2 - R[2, 4] ^ 2 - R[3, 4] ^ 2 + R[1, 2] ^ 2 * R[3, 4] ^ 2 + R[1, 3] ^ 2 * R[2, 4] ^ 2 + 2 * R[2, 3] * R[2, 4] * R[3, 4] + 2 * R[1, 2] * R[1, 3] * R[2, 3] - 2 * R[1, 2] * R[1, 3] * R[2, 4] * R[3, 4]
      
      L14 = (- b14 + sqrt((b14 ^ 2) - 4 * a14 * c14)) / (2 * a14)
      U14 = (- b14 - sqrt((b14 ^ 2) - 4 * a14 * c14)) / (2 * a14)
      L14[is.na(L14)] = - 1
      U14[is.na(U14)] = 1
      
      low14 = ceiling(100 * max(0, min(L14, U14)))
      up14 = floor(100 * max(min(1, L14), min(1, U14)))
      #up14 = floor(100 * min(R[1, 2], R[1, 3], R[2, 4], R[3, 4], max(L14, U14)))
      d14 = up14 - low14 + 1
      d14[is.na(d14)] = 1
      fr14 = matrix(rep(0, d14 * 4), d14, 4)
      for (k in low14:up14){
        r = k / 100
        Rho = R
        Rho[1, 4] = r
        Rho[4, 1] = r
        fr14[(k - low14 + 1), 1] = r
        summand = apply(resid, 1, function(resid) resid %*% ginv(S %*% Rho %*% S) %*% resid)
        
        fr14[(k - low14 + 1), 2] = uniform_prior(n = n, R = Rho, j = sum(summand))
      }
      fr14 = matrix(fr14, ncol = 4)
      fr14 = fr14[!is.na(fr14[, 2]), ] 
      fr14 = matrix(fr14, ncol = 4)
      fr14 = fr14[!is.infinite(fr14[, 2]), ] 
      fr14 = matrix(fr14, ncol = 4)
      fr14 = fr14[complete.cases(fr14), ] 
      fr14 = matrix(fr14, ncol = 4)
      fr14[, 2] = fr14[, 2] - median(fr14[, 2])
      fr14 = matrix(fr14, ncol = 4)
      
      fr14[, 2] = exp(fr14[, 2])
      fr14 = matrix(fr14, ncol = 4)
      
      fr14 = fr14[!is.infinite(fr14[, 2]), ] 
      fr14 = matrix(fr14, ncol = 4)
      
      m = which(fr14[, 2] == max(fr14[, 2]))[1]
      x1 = max(1, (m - 10))
      x2 = min(length(fr14[, 1]), (m + 10))
      low14 = 500 * fr14[x1, 1]
      up14 = 500 * fr14[x2, 1]
      
      d14 = up14 - low14 + 1
      d14[is.na(d14)] = 1
      fr14 = matrix(rep(0, d14 * 4), d14, 4)
      for (k in low14:up14){
        r = k / 500
        Rho = R
        Rho[1, 4] = r
        Rho[4, 1] = r
        fr14[(k - low14 + 1), 1] = r
        
        summand = apply(resid, 1, function(resid) resid %*% ginv(S %*% Rho %*% S) %*% resid)
        
        fr14[(k - low14 + 1), 2] = uniform_prior(n = n, R = Rho, j = sum(summand))
      }
      fr14 = matrix(fr14, ncol = 4)
      fr14 = fr14[!is.na(fr14[, 2]), ] 
      fr14 = matrix(fr14, ncol = 4)
      fr14 = fr14[!is.infinite(fr14[, 2]), ] 
      fr14 = matrix(fr14, ncol = 4)
      fr14 = fr14[complete.cases(fr14), ] 
      fr14 = matrix(fr14, ncol = 4)
      fr14[, 2] = fr14[, 2] - median(fr14[, 2])
      fr14 = matrix(fr14, ncol = 4)
      fr14[, 2] = exp(fr14[, 2])
      fr14 = matrix(fr14, ncol = 4)
      fr14 = fr14[!is.infinite(fr14[, 2]), ] 
      fr14 = matrix(fr14, ncol = 4)
      
      
      for (k in 1: length(fr14[, 1])){
        fr14[k, 3] = fr14[k, 2] / (sum(fr14[, 2]))
      }
      
      for (k in 1:length(fr14[, 1])){
        fr14[k, 4] = sum(fr14[1:k, 3])
      }
      
      u = runif(1, 0, 1)
      
      if (u<fr14[1, 4]){
        r14 = fr14[1, 1]
      }
      if (u > fr14[length(fr14[, 1]), 4]){
        r14 = fr14[length(fr14[, 1]), 1]
      }
      if (u >= fr14[1, 4] & u <= fr14[length(fr14[, 1]), 4]){
        nr = nearest(fr14[, 4], u)
        r14 = fr14[nr, 1]
      }
      
      r14 = r14[1]
      
      R[1, 4] = r14
      R[4, 1] = r14
      
      #r23
      a23 = (R[1, 4] ^ 2 - 1)
      b23 = 2 * R[2, 4] * R[3, 4] + 2 * R[1, 2] * R[1, 3] - 2 * R[1, 2] * R[3, 4] * R[1, 4] - 2 * R[1, 3] * R[1, 4] * R[2, 4]
      c23 = 1 - R[3, 4] ^ 2 - R[2, 4] ^ 2 - R[1, 2] ^ 2 - R[1, 3] ^ 2 - R[1, 4] ^ 2 + R[1, 2] ^ 2 * R[3, 4] ^ 2 + R[1, 3] ^ 2 * R[2, 4] ^ 2 + 2 * R[1, 2] * R[1, 4] * R[2, 4] + 2 * R[1, 3] * R[1, 4] * R[3, 4] - 2 * R[1, 2] * R[1, 3] * R[2, 4] * R[3, 4]
      
      L23 = (- b23 + sqrt((b23 ^ 2) - 4 * a23 * c23)) / (2 * a23)
      U23 = (- b23 - sqrt((b23 ^ 2) - 4 * a23 * c23)) / (2 * a23)
      
      low23 = ceiling(100 * max(0, min(L23, U23)))
      #up23 = floor(100 * min(R[1, 2], R[1, 3], R[2, 4], R[3, 4], max(L23, U23)))
      low23 = ceiling(100 * max(0, min(L23, U23)))
      up23 = floor(100 * max(min(1, L23), min(1, U23)))
      
      d23 = up23 - low23 + 1
      d23[is.na(d23)] = 1
      fr23 = matrix(rep(0, d23 * 4), d23, 4)
      
      for (k in low23:up23){
        r = k / 100
        Rho = R
        Rho[2, 3] = r
        Rho[3, 2] = r
        fr23[(k - low23 + 1), 1] = r
        
        summand = apply(resid, 1, function(resid) resid %*% ginv(S %*% Rho %*% S) %*% resid)
        
        fr23[(k - low23 + 1), 2] = uniform_prior(n = n, R = Rho, j = sum(summand))
      }
      
      fr23 = matrix(fr23, ncol = 4)
      
      fr23 = fr23[!is.na(fr23[, 2]), ] 
      fr23 = matrix(fr23, ncol = 4)
      
      fr23 = fr23[!is.infinite(fr23[, 2]), ]
      fr23 = matrix(fr23, ncol = 4)
      
      fr23 = fr23[complete.cases(fr23), ] 
      
      fr23 = matrix(fr23, ncol = 4)
      fr23[, 2] = fr23[, 2] - median(fr23[, 2])
      fr23 = matrix(fr23, ncol = 4)
      
      fr23[, 2] = exp(fr23[, 2])
      fr23 = matrix(fr23, ncol = 4)
      
      fr23 = fr23[!is.infinite(fr23[, 2]), ] 
      fr23 = matrix(fr23, ncol = 4)
      
      m = which(fr23[, 2] == max(fr23[, 2]))[1]
      x1 = max(1, (m - 10))
      x2 = min(length(fr23[, 1]), (m + 10))
      low23 = 500 * fr23[x1, 1]
      up23 = 500 * fr23[x2, 1]
      
      d23 = up23 - low23 + 1
      
      fr23 = matrix(rep(0, d23 * 4), d23, 4)
      
      for (k in low23:up23){
        r = k / 500
        Rho = R
        Rho[2, 3] = r
        Rho[3, 2] = r
        fr23[(k - low23 + 1), 1] = r
        
        summand = apply(resid, 1, function(resid) resid %*% ginv(S %*% Rho %*% S) %*% resid)
        
        fr23[(k - low23 + 1), 2] = uniform_prior(n = n, R = Rho, j = sum(summand))
      }
      fr23 = matrix(fr23, ncol = 4)
      
      fr23 = fr23[!is.na(fr23[, 2]), ] 
      fr23 = matrix(fr23, ncol = 4)
      
      fr23 = fr23[!is.infinite(fr23[, 2]), ]
      fr23 = matrix(fr23, ncol = 4)
      
      fr23 = fr23[complete.cases(fr23), ] 
      fr23 = matrix(fr23, ncol = 4)
      
      fr23[, 2] = fr23[, 2] - median(fr23[, 2])
      fr23 = matrix(fr23, ncol = 4)
      
      fr23[, 2] = exp(fr23[, 2])
      fr23 = matrix(fr23, ncol = 4)
      
      fr23 = fr23[!is.infinite(fr23[, 2]), ] 
      fr23 = matrix(fr23, ncol = 4)
      
      
      for (k in 1: length(fr23[, 1])){
        fr23[k, 3] = fr23[k, 2] / (sum(fr23[, 2]))
      }
      
      for (k in 1:length(fr23[, 1])){
        fr23[k, 4] = sum(fr23[1:k, 3])
      }
      
      u = runif(1, 0, 1)
      
      if (u<fr23[1, 4]){
        r23 = fr23[1, 1]
      }
      if (u > fr23[length(fr23[, 1]), 4]){
        r23 = fr23[length(fr23[, 1]), 1]
      }
      if (u >= fr23[1, 4] & u <= fr23[length(fr23[, 1]), 4]){
        nr = nearest(fr23[, 4], u)
        r23 = fr23[nr, 1]
      }
      
      r23 = r23[1]
      
      R[2, 3] = r23
      R[3, 2] = r23
      
      
      #r34
      a34 = (R[1, 2] ^ 2 - 1)
      b34 = 2 * R[2, 3] * R[2, 4] + 2 * R[1, 3] * R[1, 4] - 2 * R[1, 2] * R[1, 4] * R[2, 3] - 2 * R[1, 2] * R[1, 3] * R[2, 4]
      c34 = 1 - R[1, 2] ^ 2 - R[1, 3] ^ 2 - R[1, 4] ^ 2 - R[2, 3] ^ 2 - R[2, 4] ^ 2 + R[1, 3] ^ 2 * R[2, 4] ^ 2 + R[2, 3] ^ 2 * R[1, 4] ^ 2 + 2 * R[1, 2] * R[1, 3] * R[2, 3] + 2 * R[1, 4] * R[1, 2] * R[2, 4] - 2 * R[1, 3] * R[2, 4] * R[1, 4] * R[2, 3]
      
      L34 = (- b34 + sqrt((b34 ^ 2) - 4 * a34 * c34)) / (2 * a34)
      U34 = (- b34 - sqrt((b34 ^ 2) - 4 * a34 * c34)) / (2 * a34)
      
      low34 = ceiling(100 * max(0, min(L34, U34)))
      up34 = floor(100 * max(min(1, L34), min(1, U34)))
      
      d34 = up34 - low34 + 1
      
      fr34 = matrix(rep(0, d34 * 4), d34, 4)
      for (k in low34:up34){
        r = k / 100
        Rho = R
        Rho[3, 4] = r
        Rho[4, 3] = r
        fr34[(k - low34 + 1), 1] = r
        
        summand = apply(resid, 1, function(resid) resid %*% ginv(S %*% Rho %*% S) %*% resid)
        
        
        fr34[(k - low34 + 1), 2] = uniform_prior(n, Rho, j = sum(summand))
      }
      
      fr34 = fr34[!is.na(fr34[, 2]), ] 
      fr34 = fr34[!is.infinite(fr34[, 2]), ] 
      fr34 = fr34[complete.cases(fr34), ] 
      fr34[, 2] = fr34[, 2] - median(fr34[, 2])
      fr34[, 2] = exp(fr34[, 2])
      fr34 = fr34[!is.infinite(fr34[, 2]), ] 
      if (length(fr34) == 4){r34 = fr34[1]}
      if (length(fr34)!= 4){
        m = which(fr34[, 2] == max(fr34[, 2]))[1]
        x1 = max(1, (m - 10))
        x2 = min(length(fr34[, 1]), (m + 10))
        low34 = 500 * fr34[x1, 1]
        up34 = 500 * fr34[x2, 1]
        
        d34 = up34 - low34 + 1
        
        fr34 = matrix(rep(0, d34 * 4), d34, 4)
        for (k in low34:up34){
          r = k / 500
          Rho = R
          Rho[3, 4] = r
          Rho[4, 3] = r
          fr34[(k - low34 + 1), 1] = r
          
          summand = apply(resid, 1, function(resid) resid %*% ginv(S %*% Rho %*% S) %*% resid)
          
          fr34[(k - low34 + 1), 2] = uniform_prior(n = n, R = Rho, j = sum(summand))
        }
        
        fr34 = fr34[!is.na(fr34[, 2]), ] 
        fr34 = fr34[!is.infinite(fr34[, 2]), ] 
        fr34 = fr34[complete.cases(fr34), ] 
        fr34[, 2] = fr34[, 2] - median(fr34[, 2])
        fr34[, 2] = exp(fr34[, 2])
        fr34 = fr34[!is.infinite(fr34[, 2]), ] 
        for (k in 1: length(fr34[, 1])){
          fr34[k, 3] = fr34[k, 2] / (sum(fr34[, 2]))
        }
        
        for (k in 1:length(fr34[, 1])){
          fr34[k, 4] = sum(fr34[1:k, 3])
        }
        
        u = runif(1, 0, 1)
        
        if (u<fr34[1, 4]){
          r34 = fr34[1, 1]
        }
        if (u > fr34[length(fr34[, 1]), 4]){
          r34 = fr34[length(fr34[, 1]), 1]
        }
        if (u >= fr34[1, 4] & u <= fr34[length(fr34[, 1]), 4]){
          nr = nearest(fr34[, 4], u)
          r34 = fr34[nr, 1]
        }
        r34 = r34[1]
      }
      R[3, 4] = r34
      R[4, 3] = r34
    
    }
    
    if(condindfit){
      R[1, 2] = R[2, 1] = R[1, 4] / R[2, 4]
      R[3, 2] = R[2, 3] = R[1, 2] * R[1, 3]
    }
    
    if(any(eigen(R)$values < 0)) next; if(any(abs(R) > 1)) next
    
    holdmu[, sim] = c(holdalpha0[sim], holdalpha01[sim], holdbeta0[sim], holdbeta01[sim])
    holdS[, , sim] = S
    holdR[, , sim] = R
    holdST[, , sim] = t(ST)
    
    if(sim %% 10 == 0) print(sim)
    
    slope[sim] = (holdR[2, 4, sim] * holdS[2, 2, sim] * holdS[4, 4, sim] - holdR[1, 4, sim] * holdS[1, 1, sim] * holdS[4, 4, sim] - holdR[2, 3, sim] * holdS[2, 2, sim] * holdS[3, 3, sim] + 
                    holdR[1, 3, sim] * holdS[1, 1, sim] * holdS[3, 3, sim]) / (holdS[1, 1, sim] ^ 2 + holdS[2, 2, sim] ^ 2 - 2 * holdR[1, 2, sim] * holdS[1, 1, sim] * holdS[2, 2, sim])
    
    int[sim] = (holdmu[4, sim] - holdmu[3, sim]) - ((holdR[2, 4, sim] * holdS[2, 2, sim] * holdS[4, 4, sim] - holdR[1, 4, sim] * holdS[1, 1, sim] * holdS[4, 4, sim] - holdR[2, 3, sim] * holdS[2, 2, sim] * holdS[3, 3, sim] + 
                                                       holdR[1, 3, sim] * holdS[1, 1, sim] * holdS[3, 3, sim]) / (holdS[1, 1, sim] ^ 2 + holdS[2, 2, sim] ^ 2 - 2 * holdR[1, 2, sim] * holdS[1, 1, sim] * holdS[2, 2, sim])) * (holdmu[2, sim] - holdmu[1, sim])
    
    sim = sim  +  1
    
  } 
  
  params_matrix = data.frame(holdpsi1 = holdpsi1, holdpsi2 = holdpsi2, holdomega1 = holdomega1, holdomega2 = holdomega2, 
                             holdalpha0 = holdalpha0, holdalpha01 = holdalpha01, holdbeta0 = holdbeta0, holdbeta01 = holdbeta01, 
                             int = int, slope = slope, 
                             r12 = holdR[1, 2, ], r13 = holdR[1, 3, ], r14 = holdR[1, 4, ], 
                             r23 = holdR[2, 3, ], r24 = holdR[2, 4, ], r34 = holdR[3, 4, ], 
                             s1 = holdS[1, 1, ], s2 = holdS[2, 2, ], s3 = holdS[3, 3, ], s4 = holdS[4, 4, ], 
                             holdpsi1SE = sd(holdpsi1, na.rm = T), holdpsi2SE = sd(holdpsi2, na.rm = T), holdomega1SE = sd(holdomega1, na.rm = T), 
                             holdomega2SE = sd(holdomega2, na.rm = T), 
                             holdalpha0SE = sd(holdalpha0, na.rm = T), holdalpha01SE = sd(holdalpha01, na.rm = T), holdbeta0SE = sd(holdbeta0, na.rm = T), 
                             holdbeta01SE = sd(holdbeta01, na.rm = T), 
                             intSE = sd(int, na.rm = T), slopeSE = sd(slope, na.rm = T), 
                             r12SE = sd(holdR[1, 2, ], na.rm = T), r13SE = sd(holdR[1, 3, ], na.rm = T), r14SE = sd(holdR[1, 4, ], na.rm = T), 
                             r23SE = sd(holdR[2, 3, ], na.rm = T), r24SE = sd(holdR[2, 4, ], na.rm = T), r34SE = sd(holdR[3, 4, ], na.rm = T), 
                             s1SE = sd(holdS[1, 1, ], na.rm = T), s2SE = sd(holdS[2, 2, ], na.rm = T), s3SE = sd(holdS[3, 3, ], na.rm = T), 
                             s4SE = sd(holdS[4, 4, ], na.rm = T))
  
  params = list(holdS = holdS, holdR = holdR)
  
  colMeans(params_matrix[2:sim, ], na.rm = T)
  
  result = list(params = params, params_matrix = params_matrix, 
                args = list(SIM = SIM, burnin = burnin, n = n))
  
  return(result)
  
}

