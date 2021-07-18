library(loo)
library(rjags)
library(mvtnorm)
library(metRology)
library(parallel)
# library(cubature)


cmdArgs <- commandArgs(trailingOnly = TRUE)
numCores <- as.integer(cmdArgs[1])
options(mc.cores=numCores)

filename1 <- cmdArgs[2]  # "gmm_g2_${SLURM_ARRAY_TASK_ID}.rda"
filename2 <- cmdArgs[3]  # "est_${SLURM_ARRAY_TASK_ID}.rda"

Nsmp <- as.numeric(cmdArgs[4])   #"300/500"
MD <- cmdArgs[5]     #"MD1 / MD2"
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
dat.td3 <- list("N" = Nsmp, "y" = resp, "Time" = ncol(resp), "alpha" = rep(a,3))
for(i in 1:Nstart){
  seed.jags = seed*i
  initial = list(".RNG.name" = "base::Wichmann-Hill", ".RNG.seed" = seed.jags)
  jags.td3 <- jags.model( file = "../JAGS/model_gmm_td_G3.txt", 
                          data=dat.td3, n.chains=1, n.adapt=Nadapt, inits = initial)
  params <- c("par", "loglik", "LS")
  samps.td3 <- coda.samples(jags.td3, params, n.iter = Niter)
  name.list = names(samps.td3[[1]][1,])
  list.par.td3 = which(grepl("par", name.list)==TRUE)
  geweke.td3 = apply(samps.td3[[1]][,list.par.td3], 2, function(x) geweke.diag(x)$z)
  if (sum(abs(geweke.td3)<2)==14) break
}

total_iter <-i

name.list = names(samps.td3[[1]][1,])
list.par.td3 = which(grepl("par", name.list)==TRUE)
list.par.ll3 = which(grepl("loglik", name.list)==TRUE)
list.LS = which(grepl("LS", name.list)==TRUE)
geweke.td3 = apply(samps.td3[[1]][,list.par.td3], 2, function(x) geweke.diag(x)$z)
ll.td3 = as.matrix(samps.td3[[1]][,list.par.ll3])

source('../sources/td_WAIC_LOOCV_DIC_gmm_sameSig_G3_RV.R')

samps.par = samps.td3[[1]][,list.par.td3]
save(samps.par, waic_c, loo_c, dic_c, w1, l1, w1_ver2, l1_ver2,
     waic_m, loo_m, dic_m, geweke.td3, total_iter, del.list, 
     file=paste0("../results/N", Nsmp, "/", MD,"/", PP, "/g3_td_", filename2))


