library(corpcor)
library(bayesSurv)
library(MASS)
library(coda)
library(mvtnorm)
library(MCMCpack)
library(Surrogate)
library(HardyWeinberg)
library(extraDistr)

if(T){
 source("~/OneDrive - University of Iowa/Rpackage code/multivariate-normal-principal-stratification/R/generatedata.R")
 source("~/OneDrive - University of Iowa/Rpackage code/multivariate-normal-principal-stratification/R/run_sim.R")
 source("~/OneDrive - University of Iowa/Rpackage code/multivariate-normal-principal-stratification/R/run_sim_unequalarms.R")
 source("~/OneDrive - University of Iowa/Rpackage code/multivariate-normal-principal-stratification/R/fdelt.R")
 source("~/OneDrive - University of Iowa/Rpackage code/multivariate-normal-principal-stratification/R/fdeltBeta.R")
 source("~/OneDrive - University of Iowa/Rpackage code/multivariate-normal-principal-stratification/R/nearest.R")
}

##simulate data
array_id = 1; array_id = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
if(is.na(array_id)) array_id = 1
n = 600
SIM = 100

mu = c(2, 3.4, 4, 5)
sig = 0.5
psi2= 0; psi1 = 0; omega1 = 0; omega2 = 0 #main effects of x

set.seed(323)
ST = generatedata(n = n, mu = mu, psi2 = psi2, psi1 = psi1, omega1 = omega1, omega2 = omega2, sig = sig)
X = ST[,5]

run_sim_obsdata = function(SIM, ST, X, n){
 
burnin = 0.3 * SIM
trt = c(rep(0, n/2), rep(1, n/2))

{holdmu= matrix(rep(0,4*SIM),4,SIM)
holdmu1 = matrix(rep(0,4*SIM),4,SIM)
holdS = array(rep(0,4*4*SIM),dim = c(4,4,SIM))
holdR = array(rep(0,4*4*SIM),dim = c(4,4,SIM))
holdST = array(rep(0,4*n*SIM),dim = c(4,n,SIM))
holdR[,,1] = R = (diag(c(1,1,1,1))); holdS[,,1] = S = (diag(c(1,1,1,1)))
# identified parameters
holdR[1,3,1] = holdR[3,1,1] = cor(ST[trt== 0,1], ST[trt== 0,3])
holdR[2,4,1] = holdR[4,2,1] = cor(ST[trt==1,2], ST[trt==1,4])

holdS[1,1,1] = sd(ST[trt== 0, 1])
holdS[2,2,1] = sd(ST[trt==1, 2])
holdS[3,3,1] = sd(ST[trt== 0, 3])
holdS[4,4,1] = sd(ST[trt==1, 4])

holdpsi1 = (rep(0,1*SIM));holdpsi2= (rep(0,1*SIM))
holdomega1 = (rep(0,1*SIM)); holdomega2= (rep(0,1*SIM))
holdalpha0 = (rep(0,1*SIM));holdalpha01 = (rep(0,1*SIM))
holdbeta0 = (rep(0,1*SIM)); holdbeta01 = (rep(0,1*SIM))

slope= array(0,c((SIM),1))
int= array(0,c((SIM),1))

SIG0Smat=matrix(c(0.1, 0.1),nrow=2, ncol = 2)
theta0S= matrix(c(1,0),nrow=2)
tauSS0 =tauSS1 =tauST0 =tauST1 = c(1)

#prelim values
R[1,2] = 0.3; R[1,3] = 0.7; R[1,4] = 0.15; R[2,3] = 0.15
R[2,4] = 0.7; R[3,4] = 0.3
for(i in 2:4){ for(j in 1:(i-1)){ R[i,j] = R[j,i]}}
 Xmat= cbind(rep(1,n),X)
sim = 2

holdalpha0[1] = coef(lm(ST[trt== 0,1] ~ X[trt== 0]))[1]
holdalpha01[1] = coef(lm(ST[trt==1,2] ~ X[trt==1]))[1]
holdbeta0[1] = coef(lm(ST[trt== 0,3] ~ X[trt== 0]))[1]
holdbeta01[1] = coef(lm(ST[trt==1,4] ~ X[trt==1]))[1]

ST = ST[,1:4]
}


holdSmatrix = F
holdRmatrix = T


while(sim<=SIM){

S = holdS[,,sim-1]

#estimate coefficients 
Lambda0t = diag(c(rep(0.1, 2))); tauST0 = holdS[1,1,sim-1]^2
v = ginv(Lambda0t + as.numeric(tauST0) * (t(Xmat[trt==0,]) %*% Xmat[trt==0,]))
m = v %*% (tauST0*t(Xmat[trt==0,]) %*% as.matrix(ST[trt==0,1])) 
betaS= c(rmvnorm(1,m,v/n))
holdalpha0[sim] = betaS[1]
holdpsi1[sim] = betaS[2]
tmp1 = (ST[trt==0,1]) - Xmat[trt==0,]%*%betaS

Lambda0t = diag(c(rep(0.1, 2))); tauST0 = holdS[2,2,sim-1]^2
v = ginv(Lambda0t + as.numeric(tauST0) * (t(Xmat[trt==1,]) %*% Xmat[trt==1,]))
m = v %*% (tauST0*t(Xmat[trt==1,]) %*% as.matrix(ST[trt==1,2])) 
betaS= c(rmvnorm(1,m,v/n))
holdalpha01[sim] = betaS[1]
holdpsi2[sim] = betaS[2]
tmp2= (ST[trt==1,2]) - Xmat[trt==1,]%*%betaS

Lambda0t = diag(c(rep(0.1, 2))); tauST0 = holdS[3,3,sim-1]^2
v = ginv(Lambda0t + as.numeric(tauST0) * (t(Xmat[trt==0,]) %*% Xmat[trt==0,]))
m = v %*% (tauST0*t(Xmat[trt==0,]) %*% as.matrix(ST[trt==0,3])) 
betaT= c(rmvnorm(1,m,v/n))
holdbeta0[sim] = betaT[1]
holdomega1[sim] = betaT[2]
tmp3= (ST[trt==0,3]) - Xmat[trt==0,]%*%betaT

Lambda0t = diag(c(rep(0.1, 2))); tauST0 = holdS[4,4,sim-1]^2
v = ginv(Lambda0t + as.numeric(tauST0) * (t(Xmat[trt==1,]) %*% Xmat[trt==1,]))
m = v %*% (tauST0*t(Xmat[trt==1,]) %*% as.matrix(ST[trt==1,4])) 
betaT= c(rmvnorm(1,m,v/n))
holdbeta01[sim] = betaT[1]
holdomega2[sim] = betaT[2]
tmp4= (ST[trt==1,4]) - Xmat[trt==1,]%*%betaT

#update entire sigma
mu= cbind(rep(holdalpha0[sim],n) + holdpsi1[sim]*X, 
 rep(holdalpha01[sim],n) + holdpsi2[sim]*X,
 rep(holdbeta0[sim],n) + holdomega1[sim]*X,
 rep(holdbeta01[sim],n) + holdomega2[sim]*X)

a = b = 0.1
s1 = MCMCpack::rinvgamma(1, shape = a + n/8, scale = (sum(tmp1^2)/2 + b))
s2 = MCMCpack::rinvgamma(1, shape = a + n/8, scale = (sum(tmp2^2)/2 + b))
s3 = MCMCpack::rinvgamma(1, shape = a + n/8, scale = (sum(tmp3^2)/2 + b))
s4 = MCMCpack::rinvgamma(1, shape = a + n/8, scale = (sum(tmp4^2)/2 + b))

S[1,1] = s1
S[2,2] = s2
S[3,3] = s3
S[4,4] = s4

if(holdSmatrix){S = holdS[,,1]}

holdmu[,sim] = c(holdalpha0[sim],holdalpha01[sim],holdbeta0[sim],holdbeta01[sim])
holdmu1[,sim] = c(holdalpha0[sim]+holdpsi1[sim],holdalpha01[sim]+holdpsi2[sim],holdbeta0[sim]+holdomega1[sim],holdbeta01[sim]+holdomega2[sim])

holdS[,,sim] = S

r13 = holdR[1,3,sim]; r24 = holdR[2,4,sim-1]

# estimate correlation parameters r00
resid = cbind(tmp1, tmp3)
phat = cor(tmp1, tmp3)

y = rnorm(1, 0.5 * log((1+holdR[1,3,sim-1]) / (1- holdR[1,3,sim-1])), sd = sqrt(1/(n-3)))

R2 = holdR[,,sim-1]; R2[1,3] = R2[3,1] = ifisherz(y)

summand = apply(resid, 1, function(resid) t(resid) %*% ginv(S[c(1,3), c(1,3)] %*% R2[c(1,3), c(1,3)] %*% S[c(1,3), c(1,3)]) %*% (resid) )

ratio = exp(fdelt(n/2, R = R2, j = sum(summand)))*(1/(1- holdR[1,3,sim-1]^2))

R2 = holdR[,,sim-1]; R2[1,3] = R2[3,1] = holdR[1,3,sim-1]

summand = apply(resid, 1, function(resid) t(resid) %*% ginv(S[c(1,3), c(1,3)] %*% R2[c(1,3), c(1,3)] %*% S[c(1,3), c(1,3)]) %*% (resid) )

ratio2 = exp(fdelt(n/2, R = R2, j = sum(summand)))*(1/(1-ifisherz(y)^2))

prob = max(0, min(1, (ratio / ratio2)))
if(is.na(prob)) next

z = rbern(1, prob)

R = holdR[,,sim-1]
r13 = ifisherz(y)*z + (1-z) * holdR[1,3,sim-1]
r13 = r13[1]
R[1,3] = R[3,1] = r13

# estimate correlation parameters r11
resid = cbind(tmp2, tmp4)
phat = cor(tmp2, tmp4)

y = rnorm(1, 0.5 * log((1+holdR[2,4,sim-1]) / (1- holdR[2,4,sim-1])), sd = sqrt(1/(n-3)))

R2 = holdR[,,sim-1]; R2[2,4] = R2[4,2] = ifisherz(y)

summand = apply(resid, 1, function(resid) t(resid) %*% ginv(S[c(2,4), c(2,4)] %*% R2[c(1,3), c(1,3)] %*% S[c(2,4), c(2,4)]) %*% (resid) )

ratio = exp(fdelt(n/2, R = R2, j = sum(summand)))*(1/(1- holdR[2,4,sim-1]^2))

R2 = holdR[,,sim-1]; R2[2,4] = R2[3,1] = holdR[2,4,sim-1]

summand = apply(resid, 1, function(resid) t(resid) %*% ginv(S[c(2,4), c(2,4)] %*% R2[c(2,4), c(2,4)] %*% S[c(2,4), c(2,4)]) %*% (resid) )

ratio2 = exp(fdelt(n/2, R = R2, j = sum(summand)))*(1/(1-ifisherz(y)^2))

prob = max(0, min(1, (ratio / ratio2)))
if(is.na(prob)) next

z = rbern(1, prob)

R = holdR[,,sim-1]
r24 = ifisherz(y)*z + (1-z) * holdR[2,4,sim-1]
r24 = r24[1]
R[2,4] = R[4,2] = r24

R[1,3] = R[3,1] = r13

# non-identified correlations

SurICA <- ICA.ContCont(T0S0=R[1,3], T1S1=R[2,4], T0T0=S[3,3]^2, T1T1=S[4,4]^2, S0S0=S[1,1]^2, S1S1=S[2,2]^2,
                       T0T1=seq(0, 1, by=.2), T0S1=seq(0, 1, by=.2), T1S0=seq(0, 1, by=.2),
                       S0S1=seq(0, 1, by=.2))

matrices = SurICA$Pos.Def
m = matrices[sample(nrow(matrices), 1), ]
R[1,2] = R[2,1] = c(m$S0S1)
R[1,4] = R[4,1] = c(m$T1S0)
R[2,3] = R[3,2] = c(m$T0S1)
R[3,4] = R[4,3] = c(m$T0T1)

if(any(eigen(R)$values<0)) next; if(any(abs(R)>1)) next


holdR[,,sim] = R
  
  
slope[sim] = (holdR[2,4,sim]*holdS[2,2,sim]*holdS[4,4,sim]-holdR[1,4,sim]*holdS[1,1,sim]*holdS[4,4,sim]-holdR[2,3,sim]*holdS[2,2,sim]*holdS[3,3,sim]+
 holdR[1,3,sim]*holdS[1,1,sim]*holdS[3,3,sim])/(holdS[1,1,sim]^2+holdS[2,2,sim]^2-2*holdR[1,2,sim]*holdS[1,1,sim]*holdS[2,2,sim])

int[sim] =(holdmu[4,sim]-holdmu[3,sim]) -((holdR[2,4,sim]*holdS[2,2,sim]*holdS[4,4,sim]-holdR[1,4,sim]*holdS[1,1,sim]*holdS[4,4,sim]-holdR[2,3,sim]*holdS[2,2,sim]*holdS[3,3,sim]+
    holdR[1,3,sim]*holdS[1,1,sim]*holdS[3,3,sim])/(holdS[1,1,sim]^2+holdS[2,2,sim]^2-2*holdR[1,2,sim]*holdS[1,1,sim]*holdS[2,2,sim]))*(holdmu[2,sim]-holdmu[1,sim])

if(sim %% 10 == 0) print(sim)

sim = sim + 1

} 

params_matrix = data.frame(holdpsi1 = holdpsi1, holdpsi2 = holdpsi2, holdomega1 = holdomega1, holdomega2 = holdomega2,
   holdalpha0 = holdalpha0, holdalpha01 = holdalpha01, holdbeta0 = holdbeta0, holdbeta01 = holdbeta01,
   int = int, slope = slope,
   r12 = holdR[1,2,], r13 = holdR[1,3,], r14 = holdR[1,4,],
   r23 = holdR[2,3,], r24 = holdR[2,4,], r34 = holdR[3,4,],
   s1 = holdS[1,1,], s2 = holdS[2,2,], s3 = holdS[3,3,], s4 = holdS[4,4,])

params = list(holdS = holdS, holdR = holdR)

colMeans(params_matrix[2:sim,], na.rm = T)

result = list(params = params, params_matrix = params_matrix,
  args = list(SIM = SIM, burnin = burnin, n = n))


return(result)


}

res = run_sim_obsdata(SIM = SIM, ST = ST, X = X, n = n)

params_matrix = res$params_matrix

holdR = res$params$holdR
holdS = res$params$holdS


 
plot_traceplots = function(params_matrix, variable){
 
 plot(eval(parse(text = paste0("params_matrix$", variable))), ylab = "Parameter Draw",
 xlab = "MCMC Iteration", main = paste("Traceplot of Parameter", variable))
 
 
}

plot_traceplots(params_matrix = params_matrix, variable = "int")
plot_traceplots(params_matrix = params_matrix, variable = "r13")
plot_traceplots(params_matrix = params_matrix, variable = "r24")
plot_traceplots(params_matrix = params_matrix, variable = "r34")
plot_traceplots(params_matrix = params_matrix, variable = "s1")
plot_traceplots(params_matrix = params_matrix, variable = "holdalpha0")
plot_traceplots(params_matrix = params_matrix, variable = "holdalpha01")
plot_traceplots(params_matrix = params_matrix, variable = "holdbeta0")
plot_traceplots(params_matrix = params_matrix, variable = "holdbeta01")

other_methods = function(){
 gam = array(0,c((sim-1),1))
 for (i in 1:(sim-1)){
 gam[i] = 0.5*((holdS[3,3,i]/holdS[1,1,i])*holdR[1,3,i]+(holdS[4,4,i]/holdS[2,2,i])*holdR[2,4,i])
 }
 alph= array(0,c((sim-1),1))
 for (i in 1:(sim-1)){
 alph[i] =holdmu[4,i]-holdmu[3,i]
 }
 
 b1 = array(0,c((sim-1),1))
 for (i in 1:(sim-1)){
 b1[i] =(holdmu[4,i]-holdmu[3,i]) -((holdR[2,4,i]*holdS[4,4,i]/holdS[2,2,i])*holdmu[2,i]-(holdR[1,3,i]*holdS[3,3,i]/holdS[1,1,i])*holdmu[1,i])
 }
 b2= array(0,c((sim-1),1))
 for (i in 1:(sim-1)){
 b2[i] =(holdR[1,3,i]*holdS[3,3,i]/holdS[1,1,i])
 }
 
 S10m = (holdST[2,,1:(sim-1)]-holdST[1,,1:(sim-1)])%*%rep(1,(sim-1))/(sim-1)
 Y10m = (holdST[4,,1:(sim-1)]-holdST[3,,1:(sim-1)])%*%rep(1,(sim-1))/(sim-1)
 mod5= lm(Y10m~S10m)
 
 mu= array(0,c((sim-1),1))
 for (i in 1:(sim-1)){
 mu[i] =holdmu[2,i]-holdmu[1,i]
 }
 
 
 prentice= data.frame(bet_s=numeric(1), SEbs=numeric(1),bet_trtt=numeric(1), SEbtt=numeric(1), 
   bet_trts=numeric(1), SEbts=numeric(1), gam_trt=numeric(1), SEgt=numeric(1),
   gam_s=numeric(1), SEgs=numeric(1))
 
 pren.post= data.frame(bet_sp=numeric(1), SEbet_sp=numeric(1),bet_trttp=numeric(1), SEbet_trttp=numeric(1), 
  bet_trtsp=numeric(1), SEbet_trtsp=numeric(1), gam_trtp=numeric(1), SEgam_trtp=numeric(1), 
  gam_sp=numeric(1),SEgam_sp=numeric(1))
 
 
 prentice[1] = summary(mod)$coefficient[2,1] ##effect of s on t
 prentice[2] = summary(mod)$coefficient[2,2] ## SE 
 prentice[3] = summary(mod1)$coefficient[2,1] ##effect of trt on t
 prentice[4] = summary(mod1)$coefficient[2,2] ## SE 
 prentice[5] = summary(mod2)$coefficient[2,1] ##effect of trt on s
 prentice[6] = summary(mod2)$coefficient[2,2] ## SE
 prentice[7] = summary(mod3)$coefficient[2,1] ##effect of trt on t given s
 prentice[8] = summary(mod3)$coefficient[2,2] ## SE 
 prentice[9] = summary(mod3)$coefficient[3,1] ##effect of s on t
 prentice[10] = summary(mod3)$coefficient[3,2] ## SE 
 
 
 PS[1] = summary(mod4)$coefficient[1,1] ##effect of s1-s0 on t1-t0 int
 PS[2] = summary(mod4)$coefficient[1,2] ## SE
 PS[3] = summary(mod4)$coefficient[2,1] ##effect of s1-s0 on t1-t0 slope
 PS[4] = summary(mod4)$coefficient[2,2] ## SE 
 
 naive6= lm(Y10m~S10m+X)
 
 naiveint= summary(naive6)$coefficient[1,1]
 naivedeltas= summary(naive6)$coefficient[2,1]
 naivex= summary(naive6)$coefficient[3,1]
 naivedeltase= summary(naive6)$coefficient[2,2]
 naivexse= summary(naive6)$coefficient[3,2]
 
 naive7= lm(Y10m~S10m+X+S10m*X)
 
 naiveint2= summary(naive7)$coefficient[1,1]
 naivedeltas2= summary(naive7)$coefficient[2,1]
 naivex2= summary(naive7)$coefficient[3,1]
 naiveinteract2= summary(naive7)$coefficient[4,1]
 naivedeltase2= summary(naive7)$coefficient[2,2]
 naivexse2= summary(naive7)$coefficient[3,2]
 naiveinteractse2= summary(naive7)$coefficient[4,2]
 
 naiveresults= cbind(naiveint,naivedeltas,naivex,naivedeltase,
  naivedeltase,naivexse,naiveint2,naivedeltas2,
  naivex2,naiveinteract2,naivedeltase2,naivexse2,naiveinteractse2)
 
}

final_results = function(params_matrix, write, holdR, holdS, res){
 ## save results
 param = params_matrix
 
 n = res$args$n
 burnin = res$args$burnin
 sim = res$args$SIM
 
params= data.frame(mus0 =numeric(1), SEmus0 =numeric(1), mus1 =numeric(1), SEmus1 =numeric(1),
  mut0 =numeric(1), SEmut0 =numeric(1), mut1 =numeric(1), SEmut1 =numeric(1),
  sigs0 =numeric(1), SEsigs0 =numeric(1), sigs1 =numeric(1), SEsigs1 =numeric(1),
  sigt0 =numeric(1), SEsigt0 =numeric(1),sigt1 =numeric(1), SEsigt1 =numeric(1), 
  ps=numeric(1), SEps=numeric(1),p00 =numeric(1), SEp00 =numeric(1),
  p01 =numeric(1),SEp01 =numeric(1),p10 =numeric(1),SEp10 =numeric(1),p11 =numeric(1),
  SEp11 =numeric(1),pt=numeric(1),SEpt=numeric(1), S0l=numeric(1), S0u=numeric(1),
  S1l=numeric(1), S1u=numeric(1), T0l=numeric(1),T0u=numeric(1), T1l=numeric(1),
  T1u=numeric(1), psl=numeric(1), PSu=numeric(1),p00l=numeric(1), p00u=numeric(1),
  p01l=numeric(1),p01u=numeric(1),p10l=numeric(1),p10u=numeric(1),p11l=numeric(1),
  p11u=numeric(1),ptl=numeric(1),ptu=numeric(1))

PS= data.frame(dat_int=numeric(1), dat_intSE=numeric(1),dat_sl=numeric(1), dat_slSE=numeric(1), 
  mean_int=numeric(1), SEmean_int=numeric(1),L_int=numeric(1), U_int=numeric(1),mean_sl=numeric(1)
  ,SEmean_sl=numeric(1),L_sl=numeric(1),U_sl=numeric(1),
  postdat_int=numeric(1), postdat_intSE=numeric(1), postdat_sl=numeric(1), postdat_slSE=numeric(1),
  mean_int=numeric(1), SEmean_int=numeric(1),L_int=numeric(1), U_int=numeric(1), covslint=numeric(1),covslint1 =numeric(1),
  int_coverage=numeric(1), slope_coverage=numeric(1), int1_coverage=numeric(1))

covs= data.frame(psl=numeric(1), psu=numeric(1),p00l=numeric(1), p00u=numeric(1),
  p01l=numeric(1),p01u=numeric(1),p10l=numeric(1),p10u=numeric(1),p11l=numeric(1),
  p11u=numeric(1),ptl=numeric(1),ptu=numeric(1),s0l=numeric(1),s0u=numeric(1)
  ,s1l=numeric(1),s1u=numeric(1),t0l=numeric(1),t0u=numeric(1),t1l=numeric(1),t1u=numeric(1),
  psind=numeric(1),p00ind=numeric(1),p01ind=numeric(1),p10ind=numeric(1),p11ind=numeric(1),ptind=numeric(1),
  s0ind=numeric(1),s1ind=numeric(1),t0ind=numeric(1),t1ind=numeric(1))


covs[2] = quantile(holdR[1,2,burnin:sim-1], probs = 0.975,na.rm =T)
covs[3] = quantile(holdR[1,3,burnin:sim-1], probs = 0.025,na.rm =T)
covs[4] = quantile(holdR[1,3,burnin:sim-1], probs = 0.975,na.rm =T)
covs[5] = quantile(holdR[1,4,burnin:sim-1], probs = 0.025,na.rm =T)
covs[6] = quantile(holdR[1,4,burnin:sim-1], probs = 0.975,na.rm =T)
covs[7] = quantile(holdR[2,3,burnin:sim-1], probs = 0.025,na.rm =T)
covs[8] = quantile(holdR[2,3,burnin:sim-1], probs = 0.975,na.rm =T)
covs[9] = quantile(holdR[2,4,burnin:sim-1], probs = 0.025,na.rm =T)
covs[10] = quantile(holdR[2,4,burnin:sim-1], probs = 0.975,na.rm =T)
covs[11] = quantile(holdR[3,4,burnin:sim-1], probs = 0.025,na.rm =T)
covs[12] = quantile(holdR[3,4,burnin:sim-1], probs = 0.975,na.rm =T)

covs[13] = quantile(holdS[1,1,burnin:sim-1], probs = 0.025,na.rm =T)
covs[14] = quantile(holdS[1,1,burnin:sim-1], probs = 0.975,na.rm =T)
covs[15] = quantile(holdS[2,2,burnin:sim-1], probs = 0.025,na.rm =T)
covs[16] = quantile(holdS[2,2,burnin:sim-1], probs = 0.975,na.rm =T)
covs[17] = quantile(holdS[3,3,burnin:sim-1], probs = 0.025,na.rm =T)
covs[18] = quantile(holdS[3,3,burnin:sim-1], probs = 0.975,na.rm =T)
covs[19] = quantile(holdS[4,4,burnin:sim-1], probs = 0.025,na.rm =T)
covs[20] = quantile(holdS[4,4,burnin:sim-1], probs = 0.975,na.rm =T)

covs[21] = as.numeric((holdR[1,2,1]>covs[1])&(holdR[1,2,1]<covs[2]))
covs[22] = as.numeric((holdR[1,3,1]>covs[3])&(holdR[1,2,1]<covs[4]))
covs[23] = as.numeric((holdR[1,4,1]>covs[5])&(holdR[1,4,1]<covs[6]))
covs[24] = as.numeric((holdR[2,3,1]>covs[7])&(holdR[2,3,1]<covs[8]))
covs[25] = as.numeric((holdR[2,4,1]>covs[9])&(holdR[2,4,1]<covs[10]))
covs[26] = as.numeric((holdR[3,4,1]>covs[11])&(holdR[3,4,1]<covs[12]))
covs[27] = as.numeric((holdS[1,1,1]>covs[13])&(holdS[1,1,1]<covs[14]))
covs[28] = as.numeric((holdS[2,2,1]>covs[15])&(holdS[2,2,1]<covs[16]))
covs[29] = as.numeric((holdS[3,3,1]>covs[17])&(holdS[3,3,1]<covs[18]))
covs[30] = as.numeric((holdS[4,4,1]>covs[19])&(holdS[4,4,1]<covs[20]))

params[9] = mean(holdS[1,1,burnin:(sim-1)])
params[10] = sqrt(var(holdS[1,1,burnin:(sim-1)]))
params[11] = mean(holdS[2,2,burnin:(sim-1)])
params[12] = sqrt(var(holdS[2,2,burnin:(sim-1)]))
params[13] = mean(holdS[3,3,burnin:(sim-1)])
params[14] = sqrt(var(holdS[3,3,burnin:(sim-1)]))
params[15] = mean(holdS[4,4,burnin:(sim-1)])
params[16] = sqrt(var(holdS[4,4,burnin:(sim-1)]))

params[17] = mean(holdR[1,2,burnin:(sim-1)])
params[18] = sqrt(var(holdR[1,2,burnin:(sim-1)]))
params[19] = mean(holdR[1,3,burnin:(sim-1)])
params[20] = sqrt(var(holdR[1,3,burnin:(sim-1)]))
params[21] = mean(holdR[1,4,burnin:(sim-1)])
params[22] = sqrt(var(holdR[1,4,burnin:(sim-1)]))
params[23] = mean(holdR[2,3,burnin:(sim-1)])
params[24] = sqrt(var(holdR[2,3,burnin:(sim-1)]))
params[25] = mean(holdR[2,4,burnin:(sim-1)])
params[26] = sqrt(var(holdR[2,4,burnin:(sim-1)]))
params[27] = mean(holdR[3,4,burnin:(sim-1)])
params[28] = sqrt(var(holdR[3,4,burnin:(sim-1)]))


params[29] = quantile(holdS[1,1,burnin:sim-1], probs = 0.025,na.rm =T)
params[30] = quantile(holdS[1,1,burnin:sim-1], probs = 0.975,na.rm =T)
params[31] = quantile(holdS[2,2,burnin:sim-1], probs = 0.025,na.rm =T)
params[32] = quantile(holdS[2,2,burnin:sim-1], probs = 0.975,na.rm =T)
params[33] = quantile(holdS[3,3,burnin:sim-1], probs = 0.025,na.rm =T)
params[34] = quantile(holdS[3,3,burnin:sim-1], probs = 0.975,na.rm =T)
params[35] = quantile(holdS[4,4,burnin:sim-1], probs = 0.025,na.rm =T)
params[36] = quantile(holdS[4,4,burnin:sim-1], probs = 0.975,na.rm =T)

params[37] = quantile(holdR[1,2,burnin:sim-1], probs = 0.025,na.rm =T)
params[38] = quantile(holdR[1,2,burnin:sim-1], probs = 0.975,na.rm =T)
params[39] = quantile(holdR[1,3,burnin:sim-1], probs = 0.025,na.rm =T)
params[40] = quantile(holdR[1,3,burnin:sim-1], probs = 0.975,na.rm =T)
params[41] = quantile(holdR[1,4,burnin:sim-1], probs = 0.025,na.rm =T)
params[42] = quantile(holdR[1,4,burnin:sim-1], probs = 0.975,na.rm =T)
params[43] = quantile(holdR[2,3,burnin:sim-1], probs = 0.025,na.rm =T)
params[44] = quantile(holdR[2,3,burnin:sim-1], probs = 0.975,na.rm =T)
params[45] = quantile(holdR[2,4,burnin:sim-1], probs = 0.025,na.rm =T)
params[46] = quantile(holdR[2,4,burnin:sim-1], probs = 0.975,na.rm =T)
params[47] = quantile(holdR[3,4,burnin:sim-1], probs = 0.025,na.rm =T)
params[48] = quantile(holdR[3,4,burnin:sim-1], probs = 0.975,na.rm =T)

print(params)

PS[5] = mean(param$int[burnin:sim-1],na.rm =T)
PS[6] = sqrt(var(param$int[burnin:sim-1],na.rm =T))
PS[7] = quantile(param$int[burnin:sim-1], probs = 0.025,na.rm =T)
PS[8] = quantile(param$int[burnin:sim-1], probs = 0.975,na.rm =T)
PS[9] = mean(param$slope[burnin:sim-1],na.rm =T)
PS[10] = sqrt(var(param$slope[burnin:sim-1],na.rm =T))
PS[11] = quantile(param$slope[burnin:sim-1], probs = 0.025,na.rm =T)
PS[12] = quantile(param$slope[burnin:sim-1], probs = 0.975,na.rm =T)

print(PS)

estcoef= data.frame(beta0mean=numeric(1), beta01mean=numeric(1),omega1mean=numeric(1), alpha0mean=numeric(1),
  alpha01mean=numeric(1),psi1mean=numeric(1),psi12mean=numeric(1),omega12mean=numeric(1),
  beta0l=numeric(1),beta0u=numeric(1), beta01l=numeric(1),beta01u=numeric(1),omega1l=numeric(1),
  omega1u=numeric(1), alpha0l=numeric(1),alpha0u=numeric(1),
  alpha01l=numeric(1),alpha01u=numeric(1),psi1l=numeric(1),psi1u=numeric(1),
  psi12l=numeric(1),psi12u=numeric(1),omega12l=numeric(1),omega12u=numeric(1))
  
estcoef[1] = mean(param$holdbeta0[burnin:sim-1],na.rm =T)
estcoef[2] = mean(param$holdbeta01[burnin:sim-1],na.rm =T)
estcoef[3] = mean(param$holdomega1[burnin:sim-1],na.rm =T)
estcoef[4] = mean(param$holdalpha0[burnin:sim-1],na.rm =T)
estcoef[5] = mean(param$holdalpha01[burnin:sim-1],na.rm =T)
estcoef[6] = mean(param$holdpsi1[burnin:sim-1],na.rm =T)
estcoef[7] = mean(param$holdpsi2[burnin:sim-1],na.rm =T)
estcoef[8] = mean(param$holdomega2[burnin:sim-1],na.rm =T)
estcoef[9] = quantile(param$holdbeta0[burnin:sim-1], probs = 0.025,na.rm =T)
estcoef[10] = quantile(param$holdbeta0[burnin:sim-1], probs = 0.975,na.rm =T)
estcoef[11] = quantile(param$holdbeta01[burnin:sim-1], probs = 0.025,na.rm =T)
estcoef[12] = quantile(param$holdbeta01[burnin:sim-1], probs = 0.975,na.rm =T)
estcoef[13] = quantile(param$holdomega1[burnin:sim-1], probs = 0.025,na.rm =T)
estcoef[14] = quantile(param$holdomega1[burnin:sim-1], probs = 0.975,na.rm =T)
estcoef[15] = quantile(param$holdalpha0[burnin:sim-1], probs = 0.025,na.rm =T)
estcoef[16] = quantile(param$holdalpha0[burnin:sim-1], probs = 0.975,na.rm =T)
estcoef[17] = quantile(param$holdalpha01[burnin:sim-1], probs = 0.025,na.rm =T)
estcoef[18] = quantile(param$holdalpha01[burnin:sim-1], probs = 0.975,na.rm =T)
estcoef[19] = quantile(param$holdpsi1[burnin:sim-1], probs = 0.025,na.rm =T)
estcoef[20] = quantile(param$holdpsi1[burnin:sim-1], probs = 0.975,na.rm =T)
estcoef[21] = quantile(param$holdpsi2[burnin:sim-1], probs = 0.025,na.rm =T)
estcoef[22] = quantile(param$holdpsi2[burnin:sim-1], probs = 0.975,na.rm =T)
estcoef[23] = quantile(param$holdomega2[burnin:sim-1], probs = 0.025,na.rm =T)
estcoef[24] = quantile(param$holdomega2[burnin:sim-1], probs = 0.975,na.rm =T)

print(estcoef)

if(write){
fname = paste('params',array_id,'.txt',sep="")
write.table(params, file=fname, sep="\t", row.names=F, col.names=T)
fname2 = paste('PS',array_id,'.txt',sep="")
write.table(PS, file=fname2, sep="\t", row.names=F, col.names=T)
fname3 = paste('prentice',array_id,'.txt',sep="")
#write.table(prentice, file=fname3, sep="\t", row.names=F, col.names=T)
fname4 = paste('postpren',array_id,'.txt',sep="")
#write.table(pren.post, file=fname4, sep="\t", row.names=F, col.names=T)
fname5 = paste('naivemodels',array_id,'.txt',sep="")
#write.table(naiveresults, file=fname5, sep="\t", row.names=F, col.names=T)
fname6 = paste('estimatedcoef',array_id,'.txt',sep="")
write.table(estcoef, file=fname6, sep="\t", row.names=F, col.names=T)
fname7 = paste('covs',array_id,'.txt',sep="")
#write.table(covs, file=fname7, sep="\t", col.names=T)
}

}

final_results(params_matrix = params_matrix, write = F, holdR = res$params$holdR,
  holdS = res$params$holdS, res = res)
