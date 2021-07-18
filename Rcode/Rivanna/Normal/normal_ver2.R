library(loo)
library(rjags)
library(mvtnorm)


cmdArgs <- commandArgs(trailingOnly = TRUE)
numCores <- as.integer(cmdArgs[1])
options(mc.cores=numCores)

filename1 <- cmdArgs[2]  # "gmm_g2_${SLURM_ARRAY_TASK_ID}.rda"
filename2 <- cmdArgs[3]  # "est_${SLURM_ARRAY_TASK_ID}.rda"

Nsmp <- as.numeric(cmdArgs[4])   #"300/500"
MD <- cmdArgs[5]     # "MD0 / MD1"
PP <- cmdArgs[6]     # "p5/p10/p15"

# Load Data
load(file=paste0("../data/N", Nsmp, "/", MD,"/", PP, "/", filename1))


# settings
# number of iterations & seed for initial values
seed = as.numeric(substr(filename1, 8, 8))
Nstart = 50
a = 10 # prior for lambda

# ----------------------------------------
# gcm (G=1)
# ----------------------------------------
Niter = 20000
Nadapt = 10000
burnIn = Niter/2

dat <- list("N" = Nsmp, "y" = resp, "Time" = ncol(resp))
for(i in 1:Nstart){
  seed.jags = seed*i
  initial = list(".RNG.name" = "base::Wichmann-Hill", ".RNG.seed" = seed.jags)
  jags.m1 <- jags.model( file = "../JAGS/model_gcm_nm.txt", 
                         data=dat, n.chains=1, n.adapt=Nadapt, inits = initial)
  params <- c("par", "loglik", "sig2", "LS")
  samps.nm1 <- coda.samples(jags.m1, params, n.iter = Niter)
  name.list = names(samps.nm1[[1]][1,])
  list.par.nm1 = which(grepl("par", name.list)==TRUE)
  geweke.nm1 = apply(samps.nm1[[1]][,list.par.nm1], 2, function(x) geweke.diag(x)$z)
  if (sum(abs(geweke.nm1)<2)==6) break
}

total_iter <-i

name.list = names(samps.nm1[[1]][1,])
list.par.nm1 = which(grepl("par", name.list)==TRUE)
list.par.ll1 = which(grepl("loglik", name.list)==TRUE)
list.s = which(grepl("sig2", name.list)==TRUE)
list.LS = which(grepl("LS", name.list)==TRUE)
geweke.nm1 = apply(samps.nm1[[1]][,list.par.nm1], 2, function(x) geweke.diag(x)$z)
ll.nm1 = as.matrix(samps.nm1[[1]][(burnIn+1):Niter,list.par.ll1])


source('../sources/nm_WAIC_LOOCV_DIC_gcm.R')

samps.par = samps.nm1[[1]][,list.par.nm1]
save(samps.par, waic_c, loo_c, dic_c,
     waic_m, loo_m, dic_m, geweke.nm1, total_iter, 
     file=paste0("../results/N", Nsmp, "/", MD,"/", PP, "/g1_nm_", filename2))
rm(samps.nm1, samps.par, geweke.nm1)

# ----------------------------------------
# gmm (G=2)
# ----------------------------------------
Nadapt = 10000
Niter = 10000

G = 2
dat <- list("N" = Nsmp, "y" = resp, "Time" = ncol(resp), "alpha" = rep(a,2))
for(i in 1:Nstart){
  seed.jags = seed*i
  initial = list(".RNG.name" = "base::Wichmann-Hill", ".RNG.seed" = seed.jags)
  jags.nm2 <- jags.model( file = "../JAGS/model_gmm_nm_G2.txt", 
                          data=dat, n.chains=1, n.adapt=Nadapt, inits = initial)
  params <- c("par", "loglik", "sig2", "mem", "LS")
  samps.nm2 <- coda.samples(jags.nm2, params, n.iter = Niter)
  name.list = names(samps.nm2[[1]][1,])
  list.par.nm2 = which(grepl("par", name.list)==TRUE)
  geweke.nm2 = apply(samps.nm2[[1]][,list.par.nm2], 2, function(x) geweke.diag(x)$z)
  if (sum(abs(geweke.nm2)<2)==10) break
}

total_iter <-i

name.list = names(samps.nm2[[1]][1,])
list.par.nm2 = which(grepl("par", name.list)==TRUE)
list.par.ll = which(grepl("loglik", name.list)==TRUE)
list.s = which(grepl("sig2", name.list)==TRUE)
list.LS = which(grepl("LS", name.list)==TRUE)
list.m = which(grepl("mem", name.list)==TRUE)
geweke.nm2 = apply(samps.nm2[[1]][,list.par.nm2], 2, function(x) geweke.diag(x)$z)
ll.nm = as.matrix(samps.nm2[[1]][,list.par.ll])

source('../sources/nm_WAIC_LOOCV_DIC_gmm_sameSig_G2.R')

samps.par = samps.nm2[[1]][,list.par.nm2]
save(samps.par, waic_c, loo_c, dic_c,
     waic_m, loo_m, dic_m, geweke.nm2, total_iter, 
     file=paste0("../results/N", Nsmp, "/", MD,"/", PP, "/g2_nm_", filename2))
rm(samps.nm2, samps.par, geweke.nm2, ll.nm)

# ----------------------------------------
# gmm (G=3)
# ----------------------------------------
G = 3
dat <- list("N" = Nsmp, "y" = resp, "Time" = ncol(resp), "alpha" = rep(a,3))
for(i in 1:Nstart){
  seed.jags = seed*i
  initial = list(".RNG.name" = "base::Wichmann-Hill", ".RNG.seed" = seed.jags)
  jags.nm3 <- jags.model( file = "../JAGS/model_gmm_nm_G3.txt", 
                          data=dat, n.chains=1, n.adapt=Nadapt, inits = initial)
  params <- c("par", "loglik", "sig2", "mem", "LS")
  samps.nm3 <- coda.samples(jags.nm3, params, n.iter = Niter)
  name.list = names(samps.nm3[[1]][1,])
  list.par.nm3 = which(grepl("par", name.list)==TRUE)
  geweke.nm3 = apply(samps.nm3[[1]][,list.par.nm3], 2, function(x) geweke.diag(x)$z)
  if (sum(abs(geweke.nm3)<2)==13) break
}

total_iter <-i

name.list = names(samps.nm3[[1]][1,])
list.par.nm3 = which(grepl("par", name.list)==TRUE)
list.par.ll3 = which(grepl("loglik", name.list)==TRUE)
list.s = which(grepl("sig2", name.list)==TRUE)
list.LS = which(grepl("LS", name.list)==TRUE)
list.m = which(grepl("mem", name.list)==TRUE)
geweke.nm3 = apply(samps.nm3[[1]][,list.par.nm3], 2, function(x) geweke.diag(x)$z)

ll.nm = as.matrix(samps.nm3[[1]][,list.par.ll3])

source('../sources/nm_WAIC_LOOCV_DIC_gmm_sameSig_G3.R')

samps.par = samps.nm3[[1]][,list.par.nm3]
save(samps.par, waic_c, loo_c, dic_c,
     waic_m, loo_m, dic_m, geweke.nm3, total_iter,
     file=paste0("../results/N", Nsmp, "/", MD,"/", PP, "/g3_nm_", filename2))
rm(samps.nm3, samps.par, geweke.nm3, ll.nm)

