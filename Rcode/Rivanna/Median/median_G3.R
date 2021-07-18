library(loo)
library(rjags)
library(mvtnorm)
library(nimble)
library(cubature)
library(parallel)
library(mvnfast)

cmdArgs <- commandArgs(trailingOnly = TRUE)
numCores <- as.integer(cmdArgs[1])
options(mc.cores=numCores)

filename1 <- cmdArgs[2]  # "gmm_g2_${SLURM_ARRAY_TASK_ID}.rda"
filename2 <- cmdArgs[3]  # "est_${SLURM_ARRAY_TASK_ID}.rda"

Nsmp <- as.numeric(cmdArgs[4])   #"300/500"
MD <- cmdArgs[5]     #"MD0 / MD1"
PP <- cmdArgs[6]     #"p5/p10/p15"

# Load Data
load(file=paste0("../data/N", Nsmp, "/", MD,"/", PP, "/", filename1))


# settings
# number of iterations & seed for initial values
Nadapt = 5000
Niter = 5000
seed = as.numeric(substr(filename1, 8, 8))
Nstart = 50
a = 10 # prior for lambda


# ----------------------------------------
# gmm (G=3)
# ----------------------------------------
G = 3
dat.de3 <- list("N" = Nsmp, "y" = resp, "Time" = ncol(resp), "alpha" = rep(a,3))
for(i in 1:Nstart){
  seed.jags = seed*i
  initial = list(".RNG.name" = "base::Wichmann-Hill", ".RNG.seed" = seed.jags)
  jags.de3 <- jags.model( file = "../JAGS/model_gmm_de_G3.txt", 
                          data=dat.de3, n.chains=1, n.adapt=Nadapt, inits = initial)
  params <- c("par", "loglik", "LS")
  samps.de3 <- coda.samples(jags.de3, params, n.iter = Niter)
  name.list = names(samps.de3[[1]][1,])
  list.par.de3 = which(grepl("par", name.list)==TRUE)
  geweke.de3 = apply(samps.de3[[1]][,list.par.de3], 2, function(x) geweke.diag(x)$z)
  if (sum(abs(geweke.de3)<2)==13) break
}

total_iter <-i

name.list = names(samps.de3[[1]][1,])
list.par.de3 = which(grepl("par", name.list)==TRUE)
list.par.ll3 = which(grepl("loglik", name.list)==TRUE)
list.LS = which(grepl("LS", name.list)==TRUE)
geweke.de3 = apply(samps.de3[[1]][,list.par.de3], 2, function(x) geweke.diag(x)$z)
ll.de3 = as.matrix(samps.de3[[1]][,list.par.ll3])

source('../sources/de_WAIC_LOOCV_DIC_gmm_sameSig_G3_RV2.R')

samps.par = samps.de3[[1]][,list.par.de3]
save(samps.par, waic_c, loo_c, dic_c, w1_c, l1_c, w1_m, l1_m,
     waic_m, loo_m, dic_m, geweke.de3, total_iter, del.list, 
     file=paste0("../results/N", Nsmp, "/", MD,"/", PP, "/g3_de_", filename2))



