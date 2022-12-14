library(corpcor)
library(bayesSurv)
library(MASS)
library(coda)
library(mvtnorm)
library(MCMCpack)
library(Surrogate)
library(HardyWeinberg)
library(extraDistr)

#setwd("./R")
setwd("~/Library/CloudStorage/OneDrive-UniversityofIowa/Rpackage code/multivariate-normal-principal-stratification/R")

# keep source files in a folder called R directory
for(i in 1:(length(list.files(pattern = "\\.R$")))){
  source(list.files(pattern = "\\.R$")[i]
  )
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

res = run_sim_obsdata(SIM = SIM, ST = ST, X = X, trt = c(rep(0, n/2), rep(1, n/2)))

params_matrix = res$params_matrix

plot_traceplots(params = res, variable = "int")
plot_traceplots(params = res, variable = "slope")
plot_traceplots(params = res, variable = "r13")
plot_traceplots(params = res, variable = "r24")
plot_traceplots(params = res, variable = "r34")
plot_traceplots(params = res, variable = "s1")
plot_traceplots(params = res, variable = "holdalpha0")
plot_traceplots(params = res, variable = "holdalpha01")
plot_traceplots(params = res, variable = "holdbeta0")
plot_traceplots(params = res, variable = "holdbeta01")


final_results(params_matrix = params_matrix, write = F, holdR = res$params$holdR,
  holdS = res$params$holdS, res = res)
