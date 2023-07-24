#' run the simulation
#'
#' @description run mcmc simulation
#'
#' @param SIM number of iterations of mcmc to run
#' @param ST dataset
#' @param X covariate
#' @param trt treatment vector
#' @param condindfit boolean for conditional independence
#' @param grid interval of values for assessing valid covariance matrices
#' @param min minimum value of correlations
#' @param max maximum value of correlations
#'
#' @return simulation results
#'
#' @examples 
#' example(run_sim(SIM = 1000, ST = STdat, X = X, n = 100))
run_sim_obsdata = function(SIM, ST, X, trt, condindfit, grid, min, max){
  
  if(missing(condindfit)){
    condindfit = F
  }
  
  if(missing(min)){
    min = - 1
  }
  
  if(missing(max)){
    max = 1
  }
  
  if(missing(grid)){
    grid = 0.2
  }
  
  n = nrow(ST)
  burnin = 0.3 * SIM
  n1 = sum(trt == 1, na.rm = T)
  n0 = sum(trt == 0, na.rm = T)
  
  {holdmu= matrix(rep(0, 4 * SIM), 4, SIM)
    holdmu1 = matrix(rep(0, 4 * SIM), 4, SIM)
    holdS = array(rep(0, 4 * 4 * SIM), dim = c(4, 4, SIM))
    holdR = array(rep(0, 4 * 4 * SIM), dim = c(4, 4, SIM))
    holdR[, , 1] = R = (diag(c(1, 1, 1, 1))); holdS[, , 1] = S = (diag(c(1, 1, 1, 1)))
    # identified parameters
    holdR[1, 3, 1] = holdR[3, 1, 1] = cor(ST[trt == 0, 1], ST[trt == 0, 3])
    holdR[2, 4, 1] = holdR[4, 2, 1] = cor(ST[trt == 1, 2], ST[trt == 1, 4])
    
    holdS[1, 1, 1] = sd(ST[trt == 0, 1])
    holdS[2, 2, 1] = sd(ST[trt == 1, 2])
    holdS[3, 3, 1] = sd(ST[trt == 0, 3])
    holdS[4, 4, 1] = sd(ST[trt == 1, 4])
    
    holdpsi1 = (rep(0, 1 * SIM)); holdpsi2= (rep(0, 1 * SIM))
    holdomega1 = (rep(0, 1 * SIM)); holdomega2= (rep(0, 1 * SIM))
    holdalpha0 = (rep(0, 1 * SIM)); holdalpha01 = (rep(0, 1 * SIM))
    holdbeta0 = (rep(0, 1 * SIM)); holdbeta01 = (rep(0, 1 * SIM))
    
    slope= array(0, c((SIM), 1))
    int= array(0, c((SIM), 1))
    
    SIG0Smat = matrix(c(0.1, 0.1), nrow = 2, ncol = 2)
    theta0S = matrix(c(1, 0), nrow = 2)
    tauSS0 = tauSS1 = tauST0 = tauST1 = c(1)
    
    #prelim values
    R[1, 2] = 0.3; R[1, 3] = 0.7; R[1, 4] = 0.15; R[2, 3] = 0.15
    R[2, 4] = 0.7; R[3, 4] = 0.3
    for(i in 2:4){ for(j in 1:(i - 1)){ R[i, j] = R[j, i]}}
    Xmat= cbind(rep(1, n), X)
    sim = 2
    
    holdalpha0[1] = coef(lm(ST[trt == 0, 1] ~ X[trt == 0]))[1]
    holdalpha01[1] = coef(lm(ST[trt == 1, 2] ~ X[trt == 1]))[1]
    holdbeta0[1] = coef(lm(ST[trt == 0, 3] ~ X[trt == 0]))[1]
    holdbeta01[1] = coef(lm(ST[trt == 1, 4] ~ X[trt == 1]))[1]
    
    ST = ST[, 1:4]
  }
  
  holdSmatrix = F
  
  while(sim<=SIM){
    
    S = holdS[, , sim - 1]
    
    #estimate coefficients 
    Lambda0t = diag(c(rep(0.1, 2))); tauST0 = holdS[1, 1, sim - 1] ^ 2
    v = ginv(Lambda0t + as.numeric(tauST0) * (t(Xmat[trt == 0, ]) %*% Xmat[trt == 0, ]))
    m = v %*% (tauST0 * t(Xmat[trt == 0, ]) %*% as.matrix(ST[trt == 0, 1])) 
    betaS = c(rmvnorm(1, m, v / n))
    holdalpha0[sim] = betaS[1]
    holdpsi1[sim] = betaS[2]
    tmp1 = (ST[trt == 0, 1]) - Xmat[trt == 0, ] %*% betaS
    
    Lambda0t = diag(c(rep(0.1, 2))); tauST0 = holdS[2, 2, sim - 1] ^ 2
    v = ginv(Lambda0t + as.numeric(tauST0) * (t(Xmat[trt == 1, ]) %*% Xmat[trt == 1, ]))
    m = v %*% (tauST0 * t(Xmat[trt == 1, ]) %*% as.matrix(ST[trt == 1, 2])) 
    betaS = c(rmvnorm(1, m, v / n))
    holdalpha01[sim] = betaS[1]
    holdpsi2[sim] = betaS[2]
    tmp2 = (ST[trt == 1, 2]) - Xmat[trt == 1, ] %*% betaS
    
    Lambda0t = diag(c(rep(0.1, 2))); tauST0 = holdS[3, 3, sim - 1] ^ 2
    v = ginv(Lambda0t + as.numeric(tauST0) * (t(Xmat[trt == 0, ]) %*% Xmat[trt == 0, ]))
    m = v %*% (tauST0 * t(Xmat[trt == 0, ]) %*% as.matrix(ST[trt == 0, 3])) 
    betaT = c(rmvnorm(1, m, v / n))
    holdbeta0[sim] = betaT[1]
    holdomega1[sim] = betaT[2]
    tmp3 = (ST[trt == 0, 3]) - Xmat[trt == 0, ] %*% betaT
    
    Lambda0t = diag(c(rep(0.1, 2))); tauST0 = holdS[4, 4, sim - 1] ^ 2
    v = ginv(Lambda0t + as.numeric(tauST0) * (t(Xmat[trt == 1, ]) %*% Xmat[trt == 1, ]))
    m = v %*% (tauST0 * t(Xmat[trt == 1, ]) %*% as.matrix(ST[trt == 1, 4])) 
    betaT = c(rmvnorm(1, m, v / n))
    holdbeta01[sim] = betaT[1]
    holdomega2[sim] = betaT[2]
    tmp4 = (ST[trt == 1, 4]) - Xmat[trt == 1, ] %*% betaT
    
    #update entire sigma
    mu = cbind(rep(holdalpha0[sim], n) + holdpsi1[sim] * X, 
              rep(holdalpha01[sim], n) + holdpsi2[sim] * X, 
              rep(holdbeta0[sim], n) + holdomega1[sim] * X, 
              rep(holdbeta01[sim], n) + holdomega2[sim] * X)
    
    a = b = 0.1
    s1 = MCMCpack::rinvgamma(1, shape = a + n0 / 2, scale = (sum(tmp1 ^ 2) / 2 + b))
    s2 = MCMCpack::rinvgamma(1, shape = a + n1 / 2, scale = (sum(tmp2 ^ 2) / 2 + b))
    s3 = MCMCpack::rinvgamma(1, shape = a + n0 / 2, scale = (sum(tmp3 ^ 2) / 2 + b))
    s4 = MCMCpack::rinvgamma(1, shape = a + n1 / 2, scale = (sum(tmp4 ^ 2) / 2 + b))
    
    S[1, 1] = sqrt(s1)
    S[2, 2] = sqrt(s2)
    S[3, 3] = sqrt(s3)
    S[4, 4] = sqrt(s4)
    
    if(holdSmatrix){S = holdS[, , 1]}
    
    holdmu[, sim] = c(holdalpha0[sim], holdalpha01[sim], holdbeta0[sim], holdbeta01[sim])
    holdmu1[, sim] = c(holdalpha0[sim] + holdpsi1[sim], holdalpha01[sim] + holdpsi2[sim], holdbeta0[sim] + holdomega1[sim], holdbeta01[sim] + holdomega2[sim])
    
    holdS[, , sim] = S
    
    r24 = holdR[2, 4, sim - 1]
    
    # estimate correlation parameters r00
    resid = cbind(tmp1, tmp3)
    phat = cor(tmp1, tmp3)
    
    y = rnorm(1, 0.5 * log((1 + holdR[1, 3, sim - 1]) / (1 - holdR[1, 3, sim - 1])), sd = sqrt(1 / (n - 3)))
    
    R2 = holdR[, , sim - 1]; R2[1, 3] = R2[3, 1] = ifisherz(y)
    if(any(eigen(R2)$values < 0)) next;
    
    summand = apply(resid, 1, function(resid) t(resid) %*% ginv(S[c(1, 3), c(1, 3)] %*% R2[c(1, 3), c(1, 3)] %*% S[c(1, 3), c(1, 3)]) %*% (resid))
    
    ratio = exp(fdelt(n / 2, R = R2, j = sum(summand)))*(1 / (1 - holdR[1, 3, sim - 1] ^ 2))
    
    R2 = holdR[, , sim - 1]; 
    if(any(eigen(R2)$values < 0)) next;
    
    summand = apply(resid, 1, function(resid) t(resid) %*% ginv(S[c(1, 3), c(1, 3)] %*% R2[c(1, 3), c(1, 3)] %*% S[c(1, 3), c(1, 3)]) %*% (resid))
    
    ratio2 = exp(fdelt(n / 2, R = R2, j = sum(summand)))*(1 / (1 - ifisherz(y) ^ 2))
    
    prob = max(0, min(1, (ratio / ratio2)))
    if(is.na(prob)) next
    
    z = rbern(1, prob)
    
    R = holdR[, , sim - 1]
    r13 = ifisherz(y) * z + (1 - z) * holdR[1, 3, sim - 1]
    r13 = r13[1]
    R[1, 3] = R[3, 1] = r13
    
    # estimate correlation parameters r11
    resid = cbind(tmp2, tmp4)
    phat = cor(tmp2, tmp4)
    
    y = rnorm(1, 0.5 * log((1 + holdR[2, 4, sim - 1]) / (1 - holdR[2, 4, sim - 1])), sd = sqrt(1 / (n - 3)))
    
    R2 = holdR[, , sim - 1]; R2[2, 4] = R2[4, 2] = ifisherz(y)
    
    if(any(eigen(R2)$values < 0)) next;
    
    summand = apply(resid, 1, function(resid) t(resid) %*% ginv(S[c(2, 4), c(2, 4)] %*% R2[c(2, 4), c(2, 4)] %*% S[c(2, 4), c(2, 4)]) %*% (resid))
    
    ratio = exp(fdelt(n / 2, R = R2, j = sum(summand)))*(1 / (1 - holdR[2, 4, sim - 1] ^ 2))
    
    R2 = holdR[, , sim - 1]; R2[1, 3] = R2[3, 1] = r13
    if(any(eigen(R2)$values < 0)) next;
    
    summand = apply(resid, 1, function(resid) t(resid) %*% ginv(S[c(2, 4), c(2, 4)] %*% R2[c(2, 4), c(2, 4)] %*% S[c(2, 4), c(2, 4)]) %*% (resid))
    
    ratio2 = exp(fdelt(n / 2, R = R2, j = sum(summand)))*(1 / (1 - ifisherz(y) ^ 2))
    
    prob = max(0, min(1, (ratio / ratio2)))
    if(is.na(prob)) next
    
    z = rbern(1, prob)
    
    R = holdR[, , sim - 1]
    r24 = ifisherz(y) * z + (1 - z) * holdR[2, 4, sim - 1]
    r24 = r24[1]
    R[2, 4] = R[4, 2] = r24
    
    R[1, 3] = R[3, 1] = r13
    
    # non - identified correlations
    
    SurICA <- ICA.ContCont(T0S0 = R[1, 3], T1S1 = R[2, 4], T0T0 = S[3, 3] ^ 2, T1T1 = S[4, 4] ^ 2, S0S0 = S[1, 1] ^ 2, S1S1 = S[2, 2] ^ 2, 
                            T0T1 = seq(min, max, by = grid), T0S1 = seq(min, max, by = grid), T1S0 = seq(min, max, by = grid), 
                            S0S1 = seq(min, max, by = grid)) # may be able to put in the variance or correlations as a range and do this once
    
    matrices = SurICA$Pos.Def
    m = matrices[sample(nrow(matrices), 1), ]
    R[1, 2] = R[2, 1] = c(m$S0S1)
    R[1, 4] = R[4, 1] = c(m$T1S0)
    R[2, 3] = R[3, 2] = c(m$T0S1)
    R[3, 4] = R[4, 3] = c(m$T0T1)
    
    if(condindfit){
      R[1, 2] = R[2, 1] = R[1, 4] / R[2, 4]
      R[3, 2] = R[2, 3] = R[1, 2] * R[1, 3]
    }
    
    if(any(eigen(R)$values < 0)) next; if(any(abs(R)>1)) next
    
    
    holdR[, , sim] = R
    
    
    slope[sim] = (holdR[2, 4, sim] * holdS[2, 2, sim] * holdS[4, 4, sim] - holdR[1, 4, sim] * holdS[1, 1, sim] * holdS[4, 4, sim] - holdR[2, 3, sim] * holdS[2, 2, sim] * holdS[3, 3, sim] + 
                    holdR[1, 3, sim] * holdS[1, 1, sim] * holdS[3, 3, sim]) / (holdS[1, 1, sim] ^ 2 + holdS[2, 2, sim] ^ 2 - 2*holdR[1, 2, sim] * holdS[1, 1, sim] * holdS[2, 2, sim])
    
    int[sim] =(holdmu[4, sim] - holdmu[3, sim]) - ((holdR[2, 4, sim] * holdS[2, 2, sim] * holdS[4, 4, sim] - holdR[1, 4, sim] * holdS[1, 1, sim] * holdS[4, 4, sim] - holdR[2, 3, sim] * holdS[2, 2, sim] * holdS[3, 3, sim] + 
                                                      holdR[1, 3, sim] * holdS[1, 1, sim] * holdS[3, 3, sim]) / (holdS[1, 1, sim] ^ 2 + holdS[2, 2, sim] ^ 2 - 2*holdR[1, 2, sim] * holdS[1, 1, sim] * holdS[2, 2, sim]))*(holdmu[2, sim] - holdmu[1, sim])
    
    if(sim %% 20 == 0) print(sim)
    
    sim = sim + 1
    
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
  
  colMeans(params_matrix[2 : sim, ], na.rm = T)
  
  result = list(params = params, params_matrix = params_matrix, 
                args = list(SIM = SIM, burnin = burnin, n = n))
  
  
  return(result)
  
  
}
