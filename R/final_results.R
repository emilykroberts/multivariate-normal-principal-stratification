#' results
#'
#' @description creates output of final results
#'
#' @param params_matrix results from mcmc
#' @param res results
#' @param holdR holdR correlations
#' @param holdS holdS variances
#' @param write boolean of writing results to csv
#' @param array_id identifier related to random seed for simulation replication
#'
#' @return final results
#'
#' @examples
#' example(final_results(params_matrix = params, holdR = holdR, holdS = holdS,
#' res = results, write = F))
#'
final_results = function(params_matrix, write, holdR, holdS, res, array_id){
  ## save results
  param = params_matrix
  
  n = res$args$n
  burnin = res$args$burnin
  sim = res$args$SIM
  
  params= data.frame(mus0 = numeric(1), SEmus0 = numeric(1), mus1 = numeric(1), SEmus1 = numeric(1), 
                     mut0 = numeric(1), SEmut0 = numeric(1), mut1 = numeric(1), SEmut1 = numeric(1), 
                     sigs0 = numeric(1), SEsigs0 = numeric(1), sigs1 = numeric(1), SEsigs1 = numeric(1), 
                     sigt0 = numeric(1), SEsigt0 = numeric(1), sigt1 = numeric(1), SEsigt1 = numeric(1), 
                     ps= numeric(1), SEps= numeric(1), p00 = numeric(1), SEp00 = numeric(1), 
                     p01 = numeric(1), SEp01 = numeric(1), p10 = numeric(1), SEp10 = numeric(1), p11 = numeric(1), 
                     SEp11 = numeric(1), pt = numeric(1), SEpt = numeric(1), S0l = numeric(1), S0u = numeric(1), 
                     S1l = numeric(1), S1u = numeric(1), T0l = numeric(1), T0u = numeric(1), T1l = numeric(1), 
                     T1u = numeric(1), psl = numeric(1), PSu = numeric(1), p00l = numeric(1), p00u = numeric(1), 
                     p01l = numeric(1), p01u = numeric(1), p10l = numeric(1), p10u = numeric(1), p11l = numeric(1), 
                     p11u = numeric(1), ptl = numeric(1), ptu = numeric(1))
  
  PS= data.frame(dat_int = numeric(1), dat_intSE= numeric(1), dat_sl = numeric(1), dat_slSE= numeric(1), 
                 mean_int = numeric(1), SEmean_int = numeric(1), L_int = numeric(1), U_int = numeric(1), mean_sl = numeric(1), 
                 SEmean_sl = numeric(1), L_sl = numeric(1), U_sl = numeric(1), 
                 postdat_int = numeric(1), postdat_intSE= numeric(1), postdat_sl = numeric(1), postdat_slSE= numeric(1), 
                 mean_int = numeric(1), SEmean_int = numeric(1), L_int = numeric(1), U_int = numeric(1), covslint = numeric(1), covslint1 = numeric(1), 
                 int_coverage= numeric(1), slope_coverage= numeric(1), int1_coverage= numeric(1))
  
  covs= data.frame(psl = numeric(1), psu = numeric(1), p00l = numeric(1), p00u = numeric(1), 
                   p01l = numeric(1), p01u = numeric(1), p10l = numeric(1), p10u = numeric(1), p11l = numeric(1), 
                   p11u = numeric(1), ptl = numeric(1), ptu = numeric(1), s0l = numeric(1), s0u = numeric(1), 
                   s1l = numeric(1), s1u = numeric(1), t0l = numeric(1), t0u = numeric(1), t1l = numeric(1), t1u = numeric(1), 
                   psind = numeric(1), p00ind = numeric(1), p01ind = numeric(1), p10ind = numeric(1), p11ind = numeric(1), ptind = numeric(1), 
                   s0ind = numeric(1), s1ind = numeric(1), t0ind = numeric(1), t1ind = numeric(1))
  
  
  covs[2] = quantile(holdR[1, 2, burnin : sim - 1], probs = 0.975, na.rm = T)
  covs[3] = quantile(holdR[1, 3, burnin : sim - 1], probs = 0.025, na.rm = T)
  covs[4] = quantile(holdR[1, 3, burnin : sim - 1], probs = 0.975, na.rm = T)
  covs[5] = quantile(holdR[1, 4, burnin : sim - 1], probs = 0.025, na.rm = T)
  covs[6] = quantile(holdR[1, 4, burnin : sim - 1], probs = 0.975, na.rm = T)
  covs[7] = quantile(holdR[2, 3, burnin : sim - 1], probs = 0.025, na.rm = T)
  covs[8] = quantile(holdR[2, 3, burnin : sim - 1], probs = 0.975, na.rm = T)
  covs[9] = quantile(holdR[2, 4, burnin : sim - 1], probs = 0.025, na.rm = T)
  covs[10] = quantile(holdR[2, 4, burnin : sim - 1], probs = 0.975, na.rm = T)
  covs[11] = quantile(holdR[3, 4, burnin : sim - 1], probs = 0.025, na.rm = T)
  covs[12] = quantile(holdR[3, 4, burnin : sim - 1], probs = 0.975, na.rm = T)
  
  covs[13] = quantile(holdS[1, 1, burnin : sim - 1], probs = 0.025, na.rm = T)
  covs[14] = quantile(holdS[1, 1, burnin : sim - 1], probs = 0.975, na.rm = T)
  covs[15] = quantile(holdS[2, 2, burnin : sim - 1], probs = 0.025, na.rm = T)
  covs[16] = quantile(holdS[2, 2, burnin : sim - 1], probs = 0.975, na.rm = T)
  covs[17] = quantile(holdS[3, 3, burnin : sim - 1], probs = 0.025, na.rm = T)
  covs[18] = quantile(holdS[3, 3, burnin : sim - 1], probs = 0.975, na.rm = T)
  covs[19] = quantile(holdS[4, 4, burnin : sim - 1], probs = 0.025, na.rm = T)
  covs[20] = quantile(holdS[4, 4, burnin : sim - 1], probs = 0.975, na.rm = T)
  
  covs[21] = as.numeric((holdR[1, 2, 1]>covs[1])&(holdR[1, 2, 1]<covs[2]))
  covs[22] = as.numeric((holdR[1, 3, 1]>covs[3])&(holdR[1, 2, 1]<covs[4]))
  covs[23] = as.numeric((holdR[1, 4, 1]>covs[5])&(holdR[1, 4, 1]<covs[6]))
  covs[24] = as.numeric((holdR[2, 3, 1]>covs[7])&(holdR[2, 3, 1]<covs[8]))
  covs[25] = as.numeric((holdR[2, 4, 1]>covs[9])&(holdR[2, 4, 1]<covs[10]))
  covs[26] = as.numeric((holdR[3, 4, 1]>covs[11])&(holdR[3, 4, 1]<covs[12]))
  covs[27] = as.numeric((holdS[1, 1, 1]>covs[13])&(holdS[1, 1, 1]<covs[14]))
  covs[28] = as.numeric((holdS[2, 2, 1]>covs[15])&(holdS[2, 2, 1]<covs[16]))
  covs[29] = as.numeric((holdS[3, 3, 1]>covs[17])&(holdS[3, 3, 1]<covs[18]))
  covs[30] = as.numeric((holdS[4, 4, 1]>covs[19])&(holdS[4, 4, 1]<covs[20]))
  
  params[9] = mean(holdS[1, 1, burnin:(sim - 1)])
  params[10] = sqrt(var(holdS[1, 1, burnin:(sim - 1)]))
  params[11] = mean(holdS[2, 2, burnin:(sim - 1)])
  params[12] = sqrt(var(holdS[2, 2, burnin:(sim - 1)]))
  params[13] = mean(holdS[3, 3, burnin:(sim - 1)])
  params[14] = sqrt(var(holdS[3, 3, burnin:(sim - 1)]))
  params[15] = mean(holdS[4, 4, burnin:(sim - 1)])
  params[16] = sqrt(var(holdS[4, 4, burnin:(sim - 1)]))
  
  params[17] = mean(holdR[1, 2, burnin:(sim - 1)])
  params[18] = sqrt(var(holdR[1, 2, burnin:(sim - 1)]))
  params[19] = mean(holdR[1, 3, burnin:(sim - 1)])
  params[20] = sqrt(var(holdR[1, 3, burnin:(sim - 1)]))
  params[21] = mean(holdR[1, 4, burnin:(sim - 1)])
  params[22] = sqrt(var(holdR[1, 4, burnin:(sim - 1)]))
  params[23] = mean(holdR[2, 3, burnin:(sim - 1)])
  params[24] = sqrt(var(holdR[2, 3, burnin:(sim - 1)]))
  params[25] = mean(holdR[2, 4, burnin:(sim - 1)])
  params[26] = sqrt(var(holdR[2, 4, burnin:(sim - 1)]))
  params[27] = mean(holdR[3, 4, burnin:(sim - 1)])
  params[28] = sqrt(var(holdR[3, 4, burnin:(sim - 1)]))
  
  
  params[29] = quantile(holdS[1, 1, burnin : sim - 1], probs = 0.025, na.rm = T)
  params[30] = quantile(holdS[1, 1, burnin : sim - 1], probs = 0.975, na.rm = T)
  params[31] = quantile(holdS[2, 2, burnin : sim - 1], probs = 0.025, na.rm = T)
  params[32] = quantile(holdS[2, 2, burnin : sim - 1], probs = 0.975, na.rm = T)
  params[33] = quantile(holdS[3, 3, burnin : sim - 1], probs = 0.025, na.rm = T)
  params[34] = quantile(holdS[3, 3, burnin : sim - 1], probs = 0.975, na.rm = T)
  params[35] = quantile(holdS[4, 4, burnin : sim - 1], probs = 0.025, na.rm = T)
  params[36] = quantile(holdS[4, 4, burnin : sim - 1], probs = 0.975, na.rm = T)
  
  params[37] = quantile(holdR[1, 2, burnin : sim - 1], probs = 0.025, na.rm = T)
  params[38] = quantile(holdR[1, 2, burnin : sim - 1], probs = 0.975, na.rm = T)
  params[39] = quantile(holdR[1, 3, burnin : sim - 1], probs = 0.025, na.rm = T)
  params[40] = quantile(holdR[1, 3, burnin : sim - 1], probs = 0.975, na.rm = T)
  params[41] = quantile(holdR[1, 4, burnin : sim - 1], probs = 0.025, na.rm = T)
  params[42] = quantile(holdR[1, 4, burnin : sim - 1], probs = 0.975, na.rm = T)
  params[43] = quantile(holdR[2, 3, burnin : sim - 1], probs = 0.025, na.rm = T)
  params[44] = quantile(holdR[2, 3, burnin : sim - 1], probs = 0.975, na.rm = T)
  params[45] = quantile(holdR[2, 4, burnin : sim - 1], probs = 0.025, na.rm = T)
  params[46] = quantile(holdR[2, 4, burnin : sim - 1], probs = 0.975, na.rm = T)
  params[47] = quantile(holdR[3, 4, burnin : sim - 1], probs = 0.025, na.rm = T)
  params[48] = quantile(holdR[3, 4, burnin : sim - 1], probs = 0.975, na.rm = T)
  
  print(params)
  
  PS[5] = mean(param$int[burnin : sim - 1], na.rm = T)
  PS[6] = sqrt(var(param$int[burnin : sim - 1], na.rm = T))
  PS[7] = quantile(param$int[burnin : sim - 1], probs = 0.025, na.rm = T)
  PS[8] = quantile(param$int[burnin : sim - 1], probs = 0.975, na.rm = T)
  PS[9] = mean(param$slope[burnin : sim - 1], na.rm = T)
  PS[10] = sqrt(var(param$slope[burnin : sim - 1], na.rm = T))
  PS[11] = quantile(param$slope[burnin : sim - 1], probs = 0.025, na.rm = T)
  PS[12] = quantile(param$slope[burnin : sim - 1], probs = 0.975, na.rm = T)
  
  print(PS)
  
  estcoef= data.frame(beta0mean= numeric(1), beta01mean= numeric(1), omega1mean= numeric(1), alpha0mean= numeric(1), 
                      alpha01mean= numeric(1), psi1mean= numeric(1), psi12mean= numeric(1), omega12mean= numeric(1), 
                      beta0l = numeric(1), beta0u = numeric(1), beta01l = numeric(1), beta01u = numeric(1), omega1l = numeric(1), 
                      omega1u = numeric(1), alpha0l = numeric(1), alpha0u = numeric(1), 
                      alpha01l = numeric(1), alpha01u = numeric(1), psi1l = numeric(1), psi1u = numeric(1), 
                      psi12l = numeric(1), psi12u = numeric(1), omega12l = numeric(1), omega12u = numeric(1))
  
  estcoef[1] = mean(param$holdbeta0[burnin : sim - 1], na.rm = T)
  estcoef[2] = mean(param$holdbeta01[burnin : sim - 1], na.rm = T)
  estcoef[3] = mean(param$holdomega1[burnin : sim - 1], na.rm = T)
  estcoef[4] = mean(param$holdalpha0[burnin : sim - 1], na.rm = T)
  estcoef[5] = mean(param$holdalpha01[burnin : sim - 1], na.rm = T)
  estcoef[6] = mean(param$holdpsi1[burnin : sim - 1], na.rm = T)
  estcoef[7] = mean(param$holdpsi2[burnin : sim - 1], na.rm = T)
  estcoef[8] = mean(param$holdomega2[burnin : sim - 1], na.rm = T)
  estcoef[9] = quantile(param$holdbeta0[burnin : sim - 1], probs = 0.025, na.rm = T)
  estcoef[10] = quantile(param$holdbeta0[burnin : sim - 1], probs = 0.975, na.rm = T)
  estcoef[11] = quantile(param$holdbeta01[burnin : sim - 1], probs = 0.025, na.rm = T)
  estcoef[12] = quantile(param$holdbeta01[burnin : sim - 1], probs = 0.975, na.rm = T)
  estcoef[13] = quantile(param$holdomega1[burnin : sim - 1], probs = 0.025, na.rm = T)
  estcoef[14] = quantile(param$holdomega1[burnin : sim - 1], probs = 0.975, na.rm = T)
  estcoef[15] = quantile(param$holdalpha0[burnin : sim - 1], probs = 0.025, na.rm = T)
  estcoef[16] = quantile(param$holdalpha0[burnin : sim - 1], probs = 0.975, na.rm = T)
  estcoef[17] = quantile(param$holdalpha01[burnin : sim - 1], probs = 0.025, na.rm = T)
  estcoef[18] = quantile(param$holdalpha01[burnin : sim - 1], probs = 0.975, na.rm = T)
  estcoef[19] = quantile(param$holdpsi1[burnin : sim - 1], probs = 0.025, na.rm = T)
  estcoef[20] = quantile(param$holdpsi1[burnin : sim - 1], probs = 0.975, na.rm = T)
  estcoef[21] = quantile(param$holdpsi2[burnin : sim - 1], probs = 0.025, na.rm = T)
  estcoef[22] = quantile(param$holdpsi2[burnin : sim - 1], probs = 0.975, na.rm = T)
  estcoef[23] = quantile(param$holdomega2[burnin : sim - 1], probs = 0.025, na.rm = T)
  estcoef[24] = quantile(param$holdomega2[burnin : sim - 1], probs = 0.975, na.rm = T)
  
  print(estcoef)
  
  if(write){
    fname = paste('params', array_id, '.txt', sep="")
    write.table(params, file=fname, sep="\t", row.names=F, col.names= T)
    fname2 = paste('PS', array_id, '.txt', sep="")
    write.table(PS, file=fname2, sep="\t", row.names=F, col.names= T)
    fname3 = paste('prentice', array_id, '.txt', sep="")
    #write.table(prentice, file=fname3, sep="\t", row.names=F, col.names= T)
    fname4 = paste('postpren', array_id, '.txt', sep="")
    #write.table(pren.post, file=fname4, sep="\t", row.names=F, col.names= T)
    fname5 = paste('naivemodels', array_id, '.txt', sep="")
    #write.table(naiveresults, file=fname5, sep="\t", row.names=F, col.names= T)
    fname6 = paste('estimatedcoef', array_id, '.txt', sep="")
    write.table(estcoef, file=fname6, sep="\t", row.names=F, col.names= T)
    fname7 = paste('covs', array_id, '.txt', sep="")
    #write.table(covs, file=fname7, sep="\t", col.names= T)
  }
  
}
