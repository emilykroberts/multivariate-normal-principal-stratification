library(corpcor)
library(bayesSurv)
library(MASS)
library(coda)
library(mvtnorm)
library(MCMCpack)

holdSmatrix = F
holdRmatrix = F

##simulate data
array_id = 1; array_id = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
if(is.na(array_id)) array_id = 1
n = 300
SIM = 100

mu= c(2, 3.4, 4, 5)
sig = 0.5
psi2=0; psi1 =0; omega1 =0; omega2=0 #main effects of x

set.seed(323)
ST = generatedata(n = n, mu = mu, psi2 = psi2, psi1 = psi1, omega1 = omega1, omega2 = omega2, sig = sig)
X = ST[,5]

run_sim = function(SIM, ST, X, n){
 
burnin = 0.3 * SIM
trt = c(rep(0, n/2), rep(1, n/2))

{holdmu= matrix(rep(0,4*SIM),4,SIM)
holdmu1 = matrix(rep(0,4*SIM),4,SIM)

holdS= array(rep(0,4*4*SIM),dim=c(4,4,SIM))
holdR= array(rep(0,4*4*SIM),dim=c(4,4,SIM))
holdST= array(rep(0,4*n*SIM),dim=c(4,n,SIM))
holdR[,,1] = R = (diag(c(1,1,1,1))); holdS[,,1] = S = (diag(c(1,1,1,1)))
# identified parameters
holdR[1,3,1] = holdR[3,1,1] = cor(ST[trt==0,1], ST[trt==0,3])
holdR[2,4,1] = holdR[4,2,1] = cor(ST[trt==1,2], ST[trt==1,4])

holdS[1,1,1] = sd(ST[trt==0, 1])
holdS[2,2,1] = sd(ST[trt==1, 2])
holdS[3,3,1] = sd(ST[trt==0, 3])
holdS[4,4,1] = sd(ST[trt==1, 4])

holdpsi1 = (rep(0,1*SIM));holdpsi2= (rep(0,1*SIM))
holdomega1 = (rep(0,1*SIM)); holdomega2= (rep(0,1*SIM))
holdalpha0 = (rep(0,1*SIM));holdalpha01 = (rep(0,1*SIM))
holdbeta0 = (rep(0,1*SIM)); holdbeta01 = (rep(0,1*SIM))

slope= array(0,c((SIM),1))
int= array(0,c((SIM),1))

SIG0Smat=matrix(c(0.1, 0.1),nrow=2, ncol = 2)
theta0S= matrix(c(1,0),nrow=2)
tauSS0 =tauSS1 =tauST0 =tauST1 =c(1)

#prelim values
R[1,2] = 0.3; R[1,3] = 0.7; R[1,4] = 0.15; R[2,3] = 0.15
R[2,4] = 0.7; R[3,4] = 0.3
for(i in 2:4){ for(j in 1:(i-1)){ R[i,j] = R[j,i]}}
XmatS=cbind(rep(1,n),X)
sim=2

holdalpha0[1] = coef(lm(ST[trt==0,1] ~ X[trt==0]))[1]
holdalpha01[1] = coef(lm(ST[trt==1,2] ~ X[trt==1]))[1]
holdbeta0[1] = coef(lm(ST[trt==0,3] ~ X[trt==0]))[1]
holdbeta01[1] = coef(lm(ST[trt==1,4] ~ X[trt==1]))[1]

ST = ST[,1:4]
}

while(sim<=SIM){

S = holdS[,,sim-1]
R = holdR[,,sim-1]
SIG= S%*%R%*%S

mu=cbind(rep(holdalpha0[sim-1],n)+holdpsi1[sim-1]*X,
  rep(holdalpha01[sim-1],n)+holdpsi2[sim-1]*X,
  rep(holdbeta0[sim-1])+holdomega1[sim-1]*X,
  rep(holdbeta01[sim-1],n)+holdomega2[sim-1]*X)


for(i in 1:n){
if(trt[i] ==0){
ST[i,c(2,4)] = c(mu[i,c(2,4)]+(SIG[c(2,4),c(1,3)]%*%ginv(SIG[c(1,3),c(1,3)]))%*%(ST[i,c(1,3)]-mu[i,c(1,3)]))+
mvrnorm(1,c(0,0),SIG[c(2,4),c(2,4)]-SIG[c(2,4),c(1,3)]%*%ginv(SIG[c(1,3),c(1,3)])%*%SIG[c(1,3),c(2,4)])
}
if(trt[i] ==1){
ST[i,c(1,3)] = c(mu[i,c(1,3)]+(SIG[c(1,3),c(2,4)]%*%ginv(SIG[c(2,4),c(2,4)]))%*%(ST[i,c(2,4)]-mu[i,c(2,4)]))+
mvrnorm(1, c(0,0),SIG[c(1,3),c(1,3)]-SIG[c(1,3),c(2,4)]%*%ginv(SIG[c(2,4),c(2,4)])%*%SIG[c(2,4),c(1,3)])
}
 }

colMeans(ST)
cor(ST)

#estimate coefficients #Estimate S0 model
Xmat = XmatS = cbind(rep(1,n),X)

v= ginv(SIG0Smat+as.numeric(tauSS0)*(t(XmatS)%*%XmatS))
m= v%*%((diag(SIG0Smat)*theta0S)+tauSS0*t(XmatS)%*%ST[,1]) 

Lambda0t = diag(c(rep(0.1, 2))); tauST0 = holdS[1,1,sim-1]^2
v = ginv(Lambda0t + as.numeric(tauST0) * (t(Xmat) %*% Xmat))
m = v %*% (tauST0*t(Xmat) %*% as.matrix(ST[,1])) 
betaS=c(rmvnorm(1,m,v/n))
holdalpha0[sim] =betaS[1]
holdpsi1[sim] =betaS[2]
tmp1 = (ST[,1])-XmatS%*%betaS

Lambda0t = diag(c(rep(0.1, 2))); tauST0 = holdS[2,2,sim-1]^2
v = ginv(Lambda0t + as.numeric(tauST0) * (t(Xmat) %*% Xmat))
m = v %*% (tauST0*t(Xmat) %*% as.matrix(ST[,2])) 
betaS=c(rmvnorm(1,m,v/n))
holdalpha01[sim] =betaS[1]
holdpsi2[sim] =betaS[2]
tmp2= (ST[,2])-XmatS%*%betaS

Lambda0t = diag(c(rep(0.1, 2))); tauST0 = holdS[3,3,sim-1]^2
v = ginv(Lambda0t + as.numeric(tauST0) * (t(Xmat) %*% Xmat))
m = v %*% (tauST0*t(Xmat) %*% as.matrix(ST[,3])) 
betaT=c(rmvnorm(1,m,v/n))
holdbeta0[sim] =betaT[1]
holdomega1[sim] =betaT[2]
tmp3= (ST[,3])-XmatS%*%betaT

Lambda0t = diag(c(rep(0.1, 2))); tauST0 = holdS[4,4,sim-1]^2
v = ginv(Lambda0t + as.numeric(tauST0) * (t(Xmat) %*% Xmat))
m = v %*% (tauST0*t(Xmat) %*% as.matrix(ST[,4])) 
betaT=c(rmvnorm(1,m,v/n))
holdbeta01[sim] =betaT[1]
holdomega2[sim] =betaT[2]
tmp4= (ST[,4])-XmatS%*%betaT

#update entire sigma
tmp = rbind(t(tmp1),t(tmp2),t(tmp3),t(tmp4))

#cor(tmp1[trt==0], tmp3[trt==0])
#cor(tmp1[trt==1], tmp3[trt==1])

mu=cbind(rep(holdalpha0[sim],n)+holdpsi1[sim]*X,  
  rep(holdalpha01[sim],n)+holdpsi2[sim]*X,
  rep(holdbeta0[sim])+holdomega1[sim]*X,
  rep(holdbeta01[sim],n)+holdomega2[sim]*X)

resid= ST-t(matrix((mu),4,n,byrow=T))
cor(resid)
var(resid)

a = b = 0.1
s1 = rinvgamma(1, shape = a + n/2, scale = (sum(tmp1^2)/2 + b))
s2 = rinvgamma(1, shape = a + n/2, scale = (sum(tmp2^2)/2 + b))
s3 = rinvgamma(1, shape = a + n/2, scale = (sum(tmp3^2)/2 + b))
s4 = rinvgamma(1, shape = a + n/2, scale = (sum(tmp4^2)/2 + b))

s1 = rinvgamma(1, shape = a + n/4, scale = (sum(tmp1^2)/2 + b))
s2 = rinvgamma(1, shape = a + n/4, scale = (sum(tmp2^2)/2 + b))
s3 = rinvgamma(1, shape = a + n/4, scale = (sum(tmp3^2)/2 + b))
s4 = rinvgamma(1, shape = a + n/4, scale = (sum(tmp4^2)/2 + b))


s1 = sd(tmp1)
s2 = sd(tmp2)
s3 = sd(tmp3)
s4 = sd(tmp4)

S[1,1] = s1
S[2,2] = s2
S[3,3] = s3
S[4,4] = s4

#print(S)

#if(any(S>1)) break

#print(S)
if(holdSmatrix) S = holdS[,,1]


if(!holdRmatrix){
##r12
a12=(R[3,4]^2-1)
b12=2*R[1,4]*R[2,4]-2*R[1,3]*R[2,4]*R[3,4]-2*R[1,4]*R[2,3]*R[3,4]+2*R[1,3]*R[2,3]
c12=1-R[3,4]^2-R[2,3]^2-R[2,4]^2-R[1,3]^2-R[1,4]^2+R[1,3]^2*R[2,4]^2+R[1,4]^2*R[2,3]^2+2*R[2,3]*R[2,4]*R[3,4]-2*R[1,3]*R[1,4]*R[2,3]*R[2,4]+2*R[1,3]*R[1,4]*R[3,4]

L12=(-b12+sqrt((b12^2)-4*a12*c12))/(2*a12)
U12=(-b12-sqrt((b12^2)-4*a12*c12))/(2*a12)

low12=ceiling(100*max(0,min(L12,U12)))
up12=floor(100*max(min(1,L12),min(1,U12)))

d12=up12-low12+1

fr12= matrix(rep(0,d12*4),d12,4)
for (k in low12:up12){
r= k/100
Rho= R
Rho[1,2] = r
Rho[2,1] = r
fr12[(k-low12+1),1] = r

summand = apply(resid, 1, function(resid) resid%*%ginv(S%*%Rho%*%S)%*%resid)

fr12[(k-low12+1),2] = fdelt(n = n, R = Rho, j = sum(summand))
}

fr12=matrix(fr12,ncol=4)
fr12= fr12[!is.na(fr12[,2]),] 
fr12=matrix(fr12,ncol=4)
fr12= fr12[!is.infinite(fr12[,2]),] 
fr12=matrix(fr12,ncol=4)
fr12= fr12[complete.cases(fr12),]
fr12=matrix(fr12,ncol=4)
fr12[,2] =fr12[,2]-median(fr12[,2])
fr12=matrix(fr12,ncol=4)
fr12[,2] =exp(fr12[,2])
fr12=matrix(fr12,ncol=4)
fr12= fr12[!is.infinite(fr12[,2]),] 
fr12=matrix(fr12,ncol=4)

m= which(fr12[,2] ==max(fr12[,2]))[1]
x1 =max(1, (m-10)) #why 10?
x2=min(length(fr12[,1]), (m+10))
low12=500*fr12[x1,1]
up12=500*fr12[x2,1]
d12=up12-low12+1

fr12= matrix(rep(0,d12*4),d12,4)
for (k in low12:up12){
r= k/500
Rho= R
Rho[1,2] = r
Rho[2,1] = r
fr12[(k-low12+1),1] = r

summand = apply(resid, 1, function(resid) resid%*%ginv(S%*%Rho%*%S)%*%resid)

fr12[(k-low12+1),2] = fdelt(n = n, R = Rho, j = sum(summand))
}
fr12=matrix(fr12,ncol=4)

fr12=matrix(fr12,ncol=4)
fr12= fr12[!is.na(fr12[,2]),] 
fr12=matrix(fr12,ncol=4)
fr12= fr12[!is.infinite(fr12[,2]),] 
fr12=matrix(fr12,ncol=4)
fr12= fr12[complete.cases(fr12),]
fr12=matrix(fr12,ncol=4)
fr12[,2] =fr12[,2]-median(fr12[,2])
fr12=matrix(fr12,ncol=4)
fr12[,2] =exp(fr12[,2])
fr12=matrix(fr12,ncol=4)
fr12= fr12[!is.infinite(fr12[,2]),] 
fr12=matrix(fr12,ncol=4)

for (k in 1: length(fr12[,1])){
fr12[k,3] = fr12[k,2]/(sum(fr12[,2]))
}

for (k in 1:length(fr12[,1])){
fr12[k,4] = sum(fr12[1:k,3])
}

u= runif(1,0,1)

if (u<fr12[1,4]){
r12=fr12[1,1]
}
if (u>fr12[length(fr12[,1]),4]){
r12=fr12[length(fr12[,1]),1]
}
if (u>=fr12[1,4] & u<=fr12[length(fr12[,1]),4]){
nr= nearest(fr12[,4], u)
r12=fr12[nr,1]
}
r12[r12>1] =1
r12 = r12[1]

R[1,2] = r12
R[2,1] = r12


###r13
a13=(R[2,4]^2-1)
b13=2*R[1,2]*R[2,3]+2*R[1,4]*R[3,4]-2*R[1,2]*R[2,4]*R[3,4]-2*R[1,4]*R[2,3]*R[2,4]
c13=1-R[3,4]^2-R[2,3]^2-R[2,4]^2-R[1,2]^2-R[1,4]^2+R[1,2]^2*R[3,4]^2+R[1,4]^2*R[2,3]^2+2*R[2,3]*R[2,4]*R[3,4]+2*R[1,2]*R[1,4]*R[2,4]-2*R[1,2]*R[2,3]*R[3,4]*R[1,4]

L13=(-b13+sqrt((b13^2)-4*a13*c13))/(2*a13)
U13=(-b13-sqrt((b13^2)-4*a13*c13))/(2*a13)

a13=-R[3,4]*R[2,4]*R[1,4]-R[2,4]^2*R[2,3]*R[2,4]^2*R[1,4]-R[2,4]*R[1,2]+R[1,2]*R[3,4]
b13=-R[2,4]*R[1,4]*R[2,3]*R[1,2]*R[3,4]-R[2,4]^2*R[1,4]*R[1,2]-R[3,4]*R[2,4]^2*R[2,3]+
 R[1,2]*R[3,4]*R[1,4]*R[2,3]*R[2,4]+R[3,4]^2*R[2,4]-R[1,2]*R[3,4]^2*R[1,4]
c13=-(R[1,4]+R[1,4]*R[2,3]^2)*(R[2,3]+R[2,3]*R[1,4]^2)*
 (R[3,4]*R[1,2]*R[2,4]^2-R[1,2]^2*R[2,4]*R[1,4]*R[2,4]-R[1,4]^2*R[2,4]*R[2,3]*R[1,2]+R[1,2]^2*R[3,4]^2*R[2,3]*R[1,4])

low13=ceiling(100*max(0,min(L13,U13)))
up13=floor(100*max(L13,U13))
d13=up13-low13+1

fr13= matrix(rep(0,d13*4),d13,4)
for (k in low13:up13){
r= k/100
Rho= R
Rho[1,3] = r
Rho[3,1] = r
fr13[(k-low13+1),1] = r

summand = apply(resid, 1, function(resid) resid%*%ginv(S%*%Rho%*%S)%*%resid)
#summand = apply(resid[trt==0,], 1, function(resid) t(resid)%*%ginv(S%*%Rho%*%S)%*%resid)

fr13[(k-low13+1),2] = fdelt(n = n, R = Rho, j = sum(summand))
#fr13[(k-low13+1),2] = fdelt(n = n/2, R = Rho, j = sum(summand))

}

#plot(fr13[,1], fr13[,2])

fr13= fr13[!is.na(fr13[,2]),]
fr13= fr13[!is.infinite(fr13[,2]),] 
fr13= fr13[(fr13[,2]!=0),] 

fr13= fr13[complete.cases(fr13),]
fr13[,2] =fr13[,2]-median(fr13[,2])
fr13[,2] =exp(fr13[,2])
fr13= fr13[!is.infinite(fr13[,2]),] 
m= which(fr13[,2] ==max(fr13[,2]))[1]
x1 =max(1, (m-10))
x2=min(length(fr13[,1]), (m+10))
fr13[m,]

low13=500*fr13[x1,1]
up13=500*fr13[x2,1]

d13=up13-low13+1

fr13= matrix(rep(0,d13*4),d13,4)
for (k in low13:up13){
r= k/500
Rho= R
Rho[1,3] = r
Rho[3,1] = r
fr13[(k-low13+1),1] = r
summand = apply(resid, 1, function(resid) resid%*%ginv(S%*%Rho%*%S)%*%resid)

fr13[(k-low13+1),2] = fdelt(n = n, R = Rho, j = sum(summand))
}

#plot(fr13[,1], fr13[,2])

fr13= fr13[!is.na(fr13[,2]),]
fr13= fr13[!is.infinite(fr13[,2]),] 
fr13= fr13[complete.cases(fr13),]
fr13[,2] =fr13[,2]-median(fr13[,2])
fr13[,2] =exp(fr13[,2])
fr13= fr13[!is.infinite(fr13[,2]),] 
for (k in 1: length(fr13[,1])){
fr13[k,3] = fr13[k,2]/(sum(fr13[,2]))
}

for (k in 1:length(fr13[,1])){
fr13[k,4] = sum(fr13[1:k,3])
}


u= runif(1,0,1)
if (u<fr13[1,4]){
r13=fr13[1,1]
}
if (u>fr13[length(fr13[,1]),4]){
r13=fr13[length(fr13[,1]),1]
}
if (u>=fr13[1,4] & u<=fr13[length(fr13[,1]),4]){
nr= nearest(fr13[,4], u)
r13=fr13[nr,1]
}

r13 = r13[1]
R[1,3] = r13
R[3,1] = r13

###r24

a24=(R[1,3]^2-1)
b24=2*R[2,3]*R[3,4]+2*R[1,2]*R[1,4]-2*R[1,2]*R[1,3]*R[3,4]-2*R[1,3]*R[1,4]*R[2,3]
c24=1-R[3,4]^2-R[2,3]^2-R[1,2]^2-R[1,3]^2-R[1,4]^2+R[1,2]^2*R[3,4]^2+R[1,4]^2*R[2,3]^2+2*R[1,2]*R[1,3]*R[2,3]+2*R[1,3]*R[1,4]*R[3,4]-2*R[1,2]*R[2,3]*R[3,4]*R[1,4]
L24=(-b24+sqrt((b24^2)-4*a24*c24))/(2*a24)
U24=(-b24-sqrt((b24^2)-4*a24*c24))/(2*a24)

#low24=ceiling(100*max(0,R[2,3],R[1,4],min(L24,U24)))
low24=ceiling(100*max(0,min(L24,U24)))
up24=floor(100*max(min(1,L24),min(1,U24)))
d24=up24-low24+1
fr24= matrix(rep(0,d24*4),d24,4)
for (k in low24:up24){
r= k/100
Rho= R
Rho[2,4] = r
Rho[4,2] = r
fr24[(k-low24+1),1] = r
summand = apply(resid, 1, function(resid) resid%*%ginv(S%*%Rho%*%S)%*%resid)

fr24[(k-low24+1),2] = fdelt(n = n, R = Rho, j = sum(summand))
}
fr24=matrix(fr24,ncol=4)

fr24= fr24[!is.na(fr24[,2]),] 
fr24=matrix(fr24,ncol=4)
fr24= fr24[!is.infinite(fr24[,2]),]
fr24=matrix(fr24,ncol=4)
fr24= fr24[complete.cases(fr24),] 
fr24=matrix(fr24,ncol=4)
fr24[,2] =fr24[,2]-median(fr24[,2])
fr24=matrix(fr24,ncol=4)
fr24[,2] =exp(fr24[,2])
fr24=matrix(fr24,ncol=4)
fr24= fr24[!is.infinite(fr24[,2]),] 
fr24=matrix(fr24,ncol=4)

m= which(fr24[,2] ==max(fr24[,2]))[1]
x1 =max(1, (m-10))
x2=min(length(fr24[,1]), (m+10))
low24=500*fr24[x1,1]
up24=500*fr24[x2,1]
fr24[m,]

d24=up24-low24+1

fr24= matrix(rep(0,d24*4),d24,4)
for (k in low24:up24){
r= k/500
Rho= R
Rho[2,4] = r
Rho[4,2] = r
fr24[(k-low24+1),1] = r
summand = apply(resid, 1, function(resid) resid%*%ginv(S%*%Rho%*%S)%*%resid)

fr24[(k-low24+1),2] = fdelt(n = n, R = Rho, j = sum(summand))
}
fr24=matrix(fr24,ncol=4)
fr24= fr24[!is.na(fr24[,2]),] 
fr24=matrix(fr24,ncol=4)
fr24= fr24[!is.infinite(fr24[,2]),]
fr24=matrix(fr24,ncol=4)
fr24= fr24[complete.cases(fr24),] 
fr24=matrix(fr24,ncol=4)
fr24[,2] =fr24[,2]-median(fr24[,2])
fr24=matrix(fr24,ncol=4)
fr24[,2] =exp(fr24[,2])
fr24=matrix(fr24,ncol=4)
fr24= fr24[!is.infinite(fr24[,2]),] 
fr24=matrix(fr24,ncol=4)

for (k in 1: length(fr24[,1])){
fr24[k,3] = fr24[k,2]/(sum(fr24[,2]))
}

for (k in 1:length(fr24[,1])){
fr24[k,4] = sum(fr24[1:k,3])
}

u= runif(1,0,1)

if (u<fr24[1,4]){
r24=fr24[1,1]
}
if (u>fr24[length(fr24[,1]),4]){
r24=fr24[length(fr24[,1]),1]
}
if (u>=fr24[1,4] & u<=fr24[length(fr24[,1]),4]){
nr= nearest(fr24[,4], u)
r24=fr24[nr,1]
}
r24 = r24[1]

R[2,4] = r24
R[4,2] = r24

###r14
a14=(R[2,3]^2-1)
b14=2*R[1,2]*R[2,4]+2*R[1,3]*R[3,4]-2*R[1,2]*R[2,3]*R[3,4]-2*R[1,3]*R[2,3]*R[2,4]
c14=1-R[1,2]^2-R[1,3]^2-R[2,3]^2-R[2,4]^2-R[3,4]^2+R[1,2]^2*R[3,4]^2+R[1,3]^2*R[2,4]^2+2*R[2,3]*R[2,4]*R[3,4]+2*R[1,2]*R[1,3]*R[2,3]-2*R[1,2]*R[1,3]*R[2,4]*R[3,4]

L14=(-b14+sqrt((b14^2)-4*a14*c14))/(2*a14)
U14=(-b14-sqrt((b14^2)-4*a14*c14))/(2*a14)
L14[is.na(L14)] =-1
U14[is.na(U14)] =1

low14=ceiling(100*max(0,min(L14,U14)))
up14=floor(100*max(min(1,L14),min(1,U14)))
#up14=floor(100*min(R[1,2],R[1,3],R[2,4],R[3,4],max(L14,U14)))
d14=up14-low14+1
d14[is.na(d14)] =1
fr14= matrix(rep(0,d14*4),d14,4)
for (k in low14:up14){
r= k/100
Rho= R
Rho[1,4] = r
Rho[4,1] = r
fr14[(k-low14+1),1] = r
summand = apply(resid, 1, function(resid) resid%*%ginv(S%*%Rho%*%S)%*%resid)

fr14[(k-low14+1),2] = fdelt(n = n, R = Rho, j = sum(summand))
}
fr14=matrix(fr14,ncol=4)
fr14= fr14[!is.na(fr14[,2]),] 
fr14=matrix(fr14,ncol=4)
fr14= fr14[!is.infinite(fr14[,2]),] 
fr14=matrix(fr14,ncol=4)
fr14= fr14[complete.cases(fr14),] 
fr14=matrix(fr14,ncol=4)
fr14[,2] =fr14[,2]-median(fr14[,2])
fr14=matrix(fr14,ncol=4)

fr14[,2] =exp(fr14[,2])
fr14=matrix(fr14,ncol=4)

fr14= fr14[!is.infinite(fr14[,2]),] 
fr14=matrix(fr14,ncol=4)

m= which(fr14[,2] ==max(fr14[,2]))[1]
x1 =max(1, (m-10))
x2=min(length(fr14[,1]), (m+10))
low14=500*fr14[x1,1]
up14=500*fr14[x2,1]

d14=up14-low14+1
d14[is.na(d14)] =1
fr14= matrix(rep(0,d14*4),d14,4)
for (k in low14:up14){
r= k/500
Rho= R
Rho[1,4] = r
Rho[4,1] = r
fr14[(k-low14+1),1] = r

summand = apply(resid, 1, function(resid) resid%*%ginv(S%*%Rho%*%S)%*%resid)

fr14[(k-low14+1),2] = fdelt(n = n, R = Rho, j = sum(summand))
}
fr14=matrix(fr14,ncol=4)
fr14= fr14[!is.na(fr14[,2]),] 
fr14=matrix(fr14,ncol=4)
fr14= fr14[!is.infinite(fr14[,2]),] 
fr14=matrix(fr14,ncol=4)
fr14= fr14[complete.cases(fr14),] 
fr14=matrix(fr14,ncol=4)
fr14[,2] =fr14[,2]-median(fr14[,2])
fr14=matrix(fr14,ncol=4)
fr14[,2] =exp(fr14[,2])
fr14=matrix(fr14,ncol=4)
fr14= fr14[!is.infinite(fr14[,2]),] 
fr14=matrix(fr14,ncol=4)


for (k in 1: length(fr14[,1])){
fr14[k,3] = fr14[k,2]/(sum(fr14[,2]))
}

for (k in 1:length(fr14[,1])){
fr14[k,4] = sum(fr14[1:k,3])
}

u= runif(1,0,1)

if (u<fr14[1,4]){
r14=fr14[1,1]
}
if (u>fr14[length(fr14[,1]),4]){
r14=fr14[length(fr14[,1]),1]
}
if (u>=fr14[1,4] & u<=fr14[length(fr14[,1]),4]){
nr= nearest(fr14[,4], u)
r14=fr14[nr,1]
}

r14 = r14[1]

R[1,4] = r14
R[4,1] = r14

###r23
a23=(R[1,4]^2-1)
b23=2*R[2,4]*R[3,4]+2*R[1,2]*R[1,3]-2*R[1,2]*R[3,4]*R[1,4]-2*R[1,3]*R[1,4]*R[2,4]
c23=1-R[3,4]^2-R[2,4]^2-R[1,2]^2-R[1,3]^2-R[1,4]^2+R[1,2]^2*R[3,4]^2+R[1,3]^2*R[2,4]^2+2*R[1,2]*R[1,4]*R[2,4]+2*R[1,3]*R[1,4]*R[3,4]-2*R[1,2]*R[1,3]*R[2,4]*R[3,4]

L23=(-b23+sqrt((b23^2)-4*a23*c23))/(2*a23)
U23=(-b23-sqrt((b23^2)-4*a23*c23))/(2*a23)

low23=ceiling(100*max(0,min(L23,U23)))
#up23=floor(100*min(R[1,2],R[1,3],R[2,4],R[3,4],max(L23,U23)))
low23=ceiling(100*max(0,min(L23,U23)))
up23=floor(100*max(min(1,L23),min(1,U23)))

d23=up23-low23+1
d23[is.na(d23)] =1
fr23= matrix(rep(0,d23*4),d23,4)

for (k in low23:up23){
r= k/100
Rho= R
Rho[2,3] = r
Rho[3,2] = r
fr23[(k-low23+1),1] = r

summand = apply(resid, 1, function(resid) resid%*%ginv(S%*%Rho%*%S)%*%resid)

fr23[(k-low23+1),2] = fdelt(n = n, R = Rho, j = sum(summand))
}

fr23=matrix(fr23,ncol=4)

fr23= fr23[!is.na(fr23[,2]),] 
fr23=matrix(fr23,ncol=4)

fr23= fr23[!is.infinite(fr23[,2]),]
fr23=matrix(fr23,ncol=4)

fr23= fr23[complete.cases(fr23),] 

fr23=matrix(fr23,ncol=4)
fr23[,2] =fr23[,2]-median(fr23[,2])
fr23=matrix(fr23,ncol=4)

fr23[,2] =exp(fr23[,2])
fr23=matrix(fr23,ncol=4)

fr23= fr23[!is.infinite(fr23[,2]),] 
fr23=matrix(fr23,ncol=4)

m= which(fr23[,2] ==max(fr23[,2]))[1]
x1 =max(1, (m-10))
x2=min(length(fr23[,1]), (m+10))
low23=500*fr23[x1,1]
up23=500*fr23[x2,1]

d23=up23-low23+1

fr23= matrix(rep(0,d23*4),d23,4)

for (k in low23:up23){
r= k/500
Rho= R
Rho[2,3] = r
Rho[3,2] = r
fr23[(k-low23+1),1] = r

summand = apply(resid, 1, function(resid) resid%*%ginv(S%*%Rho%*%S)%*%resid)

fr23[(k-low23+1),2] = fdelt(n = n, R = Rho, j = sum(summand))
}
fr23=matrix(fr23,ncol=4)

fr23= fr23[!is.na(fr23[,2]),] 
fr23=matrix(fr23,ncol=4)

fr23= fr23[!is.infinite(fr23[,2]),]
fr23=matrix(fr23,ncol=4)

fr23= fr23[complete.cases(fr23),] 
fr23=matrix(fr23,ncol=4)

fr23[,2] =fr23[,2]-median(fr23[,2])
fr23=matrix(fr23,ncol=4)

fr23[,2] =exp(fr23[,2])
fr23=matrix(fr23,ncol=4)

fr23= fr23[!is.infinite(fr23[,2]),] 
fr23=matrix(fr23,ncol=4)


for (k in 1: length(fr23[,1])){
fr23[k,3] = fr23[k,2]/(sum(fr23[,2]))
}

for (k in 1:length(fr23[,1])){
fr23[k,4] = sum(fr23[1:k,3])
}

u= runif(1,0,1)

if (u<fr23[1,4]){
r23=fr23[1,1]
}
if (u>fr23[length(fr23[,1]),4]){
r23=fr23[length(fr23[,1]),1]
}
if (u>=fr23[1,4] & u<=fr23[length(fr23[,1]),4]){
nr= nearest(fr23[,4], u)
r23=fr23[nr,1]
}

r23 = r23[1]

R[2,3] = r23
R[3,2] = r23


###r34
a34=(R[1,2]^2-1)
b34=2*R[2,3]*R[2,4]+2*R[1,3]*R[1,4]-2*R[1,2]*R[1,4]*R[2,3]-2*R[1,2]*R[1,3]*R[2,4]
c34=1-R[1,2]^2-R[1,3]^2-R[1,4]^2-R[2,3]^2-R[2,4]^2+R[1,3]^2*R[2,4]^2+R[2,3]^2*R[1,4]^2+2*R[1,2]*R[1,3]*R[2,3]+2*R[1,4]*R[1,2]*R[2,4]-2*R[1,3]*R[2,4]*R[1,4]*R[2,3]

L34=(-b34+sqrt((b34^2)-4*a34*c34))/(2*a34)
U34=(-b34-sqrt((b34^2)-4*a34*c34))/(2*a34)

low34=ceiling(100*max(0,min(L34,U34)))
up34=floor(100*max(min(1,L34),min(1,U34)))

d34=up34-low34+1

fr34= matrix(rep(0,d34*4),d34,4)
for (k in low34:up34){
 r= k/100
 Rho= R
 Rho[3,4] = r
 Rho[4,3] = r
 fr34[(k-low34+1),1] = r
 
 summand = apply(resid, 1, function(resid) resid%*%ginv(S%*%Rho%*%S)%*%resid)

 
 fr34[(k-low34+1),2] = fdelt(n, Rho, j = sum(summand))
}

fr34= fr34[!is.na(fr34[,2]),] 
fr34= fr34[!is.infinite(fr34[,2]),] 
fr34= fr34[complete.cases(fr34),] 
fr34[,2] =fr34[,2]-median(fr34[,2])
fr34[,2] =exp(fr34[,2])
fr34= fr34[!is.infinite(fr34[,2]),] 
if (length(fr34)==4){r34=fr34[1]}
if (length(fr34)!=4){
 m= which(fr34[,2] ==max(fr34[,2]))[1]
 x1 =max(1, (m-10))
 x2=min(length(fr34[,1]), (m+10))
 low34=500*fr34[x1,1]
 up34=500*fr34[x2,1]
 
 d34=up34-low34+1
 
 fr34= matrix(rep(0,d34*4),d34,4)
 for (k in low34:up34){
 r= k/500
 Rho= R
 Rho[3,4] = r
 Rho[4,3] = r
 fr34[(k-low34+1),1] = r
 
 summand = apply(resid, 1, function(resid) resid%*%ginv(S%*%Rho%*%S)%*%resid)

 fr34[(k-low34+1),2] = fdelt(n = n, R = Rho, j = sum(summand))
 }
 
 fr34= fr34[!is.na(fr34[,2]),] 
 fr34= fr34[!is.infinite(fr34[,2]),] 
 fr34= fr34[complete.cases(fr34),] 
 fr34[,2] =fr34[,2]-median(fr34[,2])
 fr34[,2] =exp(fr34[,2])
 fr34= fr34[!is.infinite(fr34[,2]),] 
 for (k in 1: length(fr34[,1])){
 fr34[k,3] = fr34[k,2]/(sum(fr34[,2]))
 }
 
 for (k in 1:length(fr34[,1])){
 fr34[k,4] = sum(fr34[1:k,3])
 }
 
 u= runif(1,0,1)
 
 if (u<fr34[1,4]){
 r34=fr34[1,1]
 }
 if (u>fr34[length(fr34[,1]),4]){
 r34=fr34[length(fr34[,1]),1]
 }
 if (u>=fr34[1,4] & u<=fr34[length(fr34[,1]),4]){
 nr= nearest(fr34[,4], u)
 r34=fr34[nr,1]
 }
 r34 = r34[1]
}
R[3,4] = r34
R[4,3] = r34

if(any(eigen(R)$values<0)){ next}

}

holdmu[,sim] = c(holdalpha0[sim],holdalpha01[sim],holdbeta0[sim],holdbeta01[sim])
holdmu1[,sim] = c(holdalpha0[sim]+holdpsi1[sim],holdalpha01[sim]+holdpsi2[sim],holdbeta0[sim]+holdomega1[sim],holdbeta01[sim]+holdomega2[sim])

holdS[,,sim] = S
holdR[,,sim] = R
holdST[,,sim] = t(ST)

if(sim %% 10 == 0) print(sim)

slope[sim] = (holdR[2,4,sim]*holdS[2,2,sim]*holdS[4,4,sim]-holdR[1,4,sim]*holdS[1,1,sim]*holdS[4,4,sim]-holdR[2,3,sim]*holdS[2,2,sim]*holdS[3,3,sim]+
  holdR[1,3,sim]*holdS[1,1,sim]*holdS[3,3,sim])/(holdS[1,1,sim]^2+holdS[2,2,sim]^2-2*holdR[1,2,sim]*holdS[1,1,sim]*holdS[2,2,sim])

int[sim] =(holdmu[4,sim]-holdmu[3,sim])-((holdR[2,4,sim]*holdS[2,2,sim]*holdS[4,4,sim]-holdR[1,4,sim]*holdS[1,1,sim]*holdS[4,4,sim]-holdR[2,3,sim]*holdS[2,2,sim]*holdS[3,3,sim]+
       holdR[1,3,sim]*holdS[1,1,sim]*holdS[3,3,sim])/(holdS[1,1,sim]^2+holdS[2,2,sim]^2-2*holdR[1,2,sim]*holdS[1,1,sim]*holdS[2,2,sim]))*(holdmu[2,sim]-holdmu[1,sim])

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

res = run_sim(SIM = SIM, ST = ST, X = X, n = n)

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
 gam= array(0,c((sim-1),1))
 for (i in 1:(sim-1)){
 gam[i] =0.5*((holdS[3,3,i]/holdS[1,1,i])*holdR[1,3,i]+(holdS[4,4,i]/holdS[2,2,i])*holdR[2,4,i])
 }
 alph= array(0,c((sim-1),1))
 for (i in 1:(sim-1)){
 alph[i] =holdmu[4,i]-holdmu[3,i]
 }
 
 b1 = array(0,c((sim-1),1))
 for (i in 1:(sim-1)){
 b1[i] =(holdmu[4,i]-holdmu[3,i])-((holdR[2,4,i]*holdS[4,4,i]/holdS[2,2,i])*holdmu[2,i]-(holdR[1,3,i]*holdS[3,3,i]/holdS[1,1,i])*holdmu[1,i])
 }
 b2= array(0,c((sim-1),1))
 for (i in 1:(sim-1)){
 b2[i] =(holdR[1,3,i]*holdS[3,3,i]/holdS[1,1,i])
 }
 
 S10m= (holdST[2,,1:(sim-1)]-holdST[1,,1:(sim-1)])%*%rep(1,(sim-1))/(sim-1)
 Y10m= (holdST[4,,1:(sim-1)]-holdST[3,,1:(sim-1)])%*%rep(1,(sim-1))/(sim-1)
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
 
 naiveresults=cbind(naiveint,naivedeltas,naivex,naivedeltase,
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


covs[2] = quantile(holdR[1,2,burnin:sim-1], probs = 0.975,na.rm=T)
covs[3] = quantile(holdR[1,3,burnin:sim-1], probs = 0.025,na.rm=T)
covs[4] = quantile(holdR[1,3,burnin:sim-1], probs = 0.975,na.rm=T)
covs[5] = quantile(holdR[1,4,burnin:sim-1], probs = 0.025,na.rm=T)
covs[6] = quantile(holdR[1,4,burnin:sim-1], probs = 0.975,na.rm=T)
covs[7] = quantile(holdR[2,3,burnin:sim-1], probs = 0.025,na.rm=T)
covs[8] = quantile(holdR[2,3,burnin:sim-1], probs = 0.975,na.rm=T)
covs[9] = quantile(holdR[2,4,burnin:sim-1], probs = 0.025,na.rm=T)
covs[10] = quantile(holdR[2,4,burnin:sim-1], probs = 0.975,na.rm=T)
covs[11] = quantile(holdR[3,4,burnin:sim-1], probs = 0.025,na.rm=T)
covs[12] = quantile(holdR[3,4,burnin:sim-1], probs = 0.975,na.rm=T)

covs[13] = quantile(holdS[1,1,burnin:sim-1], probs = 0.025,na.rm=T)
covs[14] = quantile(holdS[1,1,burnin:sim-1], probs = 0.975,na.rm=T)
covs[15] = quantile(holdS[2,2,burnin:sim-1], probs = 0.025,na.rm=T)
covs[16] = quantile(holdS[2,2,burnin:sim-1], probs = 0.975,na.rm=T)
covs[17] = quantile(holdS[3,3,burnin:sim-1], probs = 0.025,na.rm=T)
covs[18] = quantile(holdS[3,3,burnin:sim-1], probs = 0.975,na.rm=T)
covs[19] = quantile(holdS[4,4,burnin:sim-1], probs = 0.025,na.rm=T)
covs[20] = quantile(holdS[4,4,burnin:sim-1], probs = 0.975,na.rm=T)

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


params[29] = quantile(holdS[1,1,burnin:sim-1], probs = 0.025,na.rm=T)
params[30] = quantile(holdS[1,1,burnin:sim-1], probs = 0.975,na.rm=T)
params[31] = quantile(holdS[2,2,burnin:sim-1], probs = 0.025,na.rm=T)
params[32] = quantile(holdS[2,2,burnin:sim-1], probs = 0.975,na.rm=T)
params[33] = quantile(holdS[3,3,burnin:sim-1], probs = 0.025,na.rm=T)
params[34] = quantile(holdS[3,3,burnin:sim-1], probs = 0.975,na.rm=T)
params[35] = quantile(holdS[4,4,burnin:sim-1], probs = 0.025,na.rm=T)
params[36] = quantile(holdS[4,4,burnin:sim-1], probs = 0.975,na.rm=T)

params[37] = quantile(holdR[1,2,burnin:sim-1], probs = 0.025,na.rm=T)
params[38] = quantile(holdR[1,2,burnin:sim-1], probs = 0.975,na.rm=T)
params[39] = quantile(holdR[1,3,burnin:sim-1], probs = 0.025,na.rm=T)
params[40] = quantile(holdR[1,3,burnin:sim-1], probs = 0.975,na.rm=T)
params[41] = quantile(holdR[1,4,burnin:sim-1], probs = 0.025,na.rm=T)
params[42] = quantile(holdR[1,4,burnin:sim-1], probs = 0.975,na.rm=T)
params[43] = quantile(holdR[2,3,burnin:sim-1], probs = 0.025,na.rm=T)
params[44] = quantile(holdR[2,3,burnin:sim-1], probs = 0.975,na.rm=T)
params[45] = quantile(holdR[2,4,burnin:sim-1], probs = 0.025,na.rm=T)
params[46] = quantile(holdR[2,4,burnin:sim-1], probs = 0.975,na.rm=T)
params[47] = quantile(holdR[3,4,burnin:sim-1], probs = 0.025,na.rm=T)
params[48] = quantile(holdR[3,4,burnin:sim-1], probs = 0.975,na.rm=T)

print(params)

PS[5] = mean(param$int[burnin:sim-1],na.rm=T)
PS[6] = sqrt(var(param$int[burnin:sim-1],na.rm=T))
PS[7] = quantile(param$int[burnin:sim-1], probs = 0.025,na.rm=T)
PS[8] = quantile(param$int[burnin:sim-1], probs = 0.975,na.rm=T)
PS[9] = mean(param$slope[burnin:sim-1],na.rm=T)
PS[10] = sqrt(var(param$slope[burnin:sim-1],na.rm=T))
PS[11] = quantile(param$slope[burnin:sim-1], probs = 0.025,na.rm=T)
PS[12] = quantile(param$slope[burnin:sim-1], probs = 0.975,na.rm=T)

print(PS)

estcoef= data.frame(beta0mean=numeric(1), beta01mean=numeric(1),omega1mean=numeric(1), alpha0mean=numeric(1),
    alpha01mean=numeric(1),psi1mean=numeric(1),psi12mean=numeric(1),omega12mean=numeric(1),
    beta0l=numeric(1),beta0u=numeric(1), beta01l=numeric(1),beta01u=numeric(1),omega1l=numeric(1),
    omega1u=numeric(1), alpha0l=numeric(1),alpha0u=numeric(1),
    alpha01l=numeric(1),alpha01u=numeric(1),psi1l=numeric(1),psi1u=numeric(1),
    psi12l=numeric(1),psi12u=numeric(1),omega12l=numeric(1),omega12u=numeric(1))
    
estcoef[1] = mean(param$holdbeta0[burnin:sim-1],na.rm=T)
estcoef[2] = mean(param$holdbeta01[burnin:sim-1],na.rm=T)
estcoef[3] = mean(param$holdomega1[burnin:sim-1],na.rm=T)
estcoef[4] = mean(param$holdalpha0[burnin:sim-1],na.rm=T)
estcoef[5] = mean(param$holdalpha01[burnin:sim-1],na.rm=T)
estcoef[6] = mean(param$holdpsi1[burnin:sim-1],na.rm=T)
estcoef[7] = mean(param$holdpsi2[burnin:sim-1],na.rm=T)
estcoef[8] = mean(param$holdomega2[burnin:sim-1],na.rm=T)
estcoef[9] = quantile(param$holdbeta0[burnin:sim-1], probs = 0.025,na.rm=T)
estcoef[10] = quantile(param$holdbeta0[burnin:sim-1], probs = 0.975,na.rm=T)
estcoef[11] = quantile(param$holdbeta01[burnin:sim-1], probs = 0.025,na.rm=T)
estcoef[12] = quantile(param$holdbeta01[burnin:sim-1], probs = 0.975,na.rm=T)
estcoef[13] = quantile(param$holdomega1[burnin:sim-1], probs = 0.025,na.rm=T)
estcoef[14] = quantile(param$holdomega1[burnin:sim-1], probs = 0.975,na.rm=T)
estcoef[15] = quantile(param$holdalpha0[burnin:sim-1], probs = 0.025,na.rm=T)
estcoef[16] = quantile(param$holdalpha0[burnin:sim-1], probs = 0.975,na.rm=T)
estcoef[17] = quantile(param$holdalpha01[burnin:sim-1], probs = 0.025,na.rm=T)
estcoef[18] = quantile(param$holdalpha01[burnin:sim-1], probs = 0.975,na.rm=T)
estcoef[19] = quantile(param$holdpsi1[burnin:sim-1], probs = 0.025,na.rm=T)
estcoef[20] = quantile(param$holdpsi1[burnin:sim-1], probs = 0.975,na.rm=T)
estcoef[21] = quantile(param$holdpsi2[burnin:sim-1], probs = 0.025,na.rm=T)
estcoef[22] = quantile(param$holdpsi2[burnin:sim-1], probs = 0.975,na.rm=T)
estcoef[23] = quantile(param$holdomega2[burnin:sim-1], probs = 0.025,na.rm=T)
estcoef[24] = quantile(param$holdomega2[burnin:sim-1], probs = 0.975,na.rm=T)

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
