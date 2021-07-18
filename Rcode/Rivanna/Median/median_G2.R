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
# gmm (G=2)
# ----------------------------------------
G = 2
dat.de2 <- list("N" = Nsmp, "y" = resp, "Time" = ncol(resp), "alpha" = rep(a,2))
for(i in 1:Nstart){
  seed.jags = seed*i
  initial = list(".RNG.name" = "base::Wichmann-Hill", ".RNG.seed" = seed.jags)
  jags.de2 <- jags.model( file = "../JAGS/model_gmm_de_G2.txt", 
                          data=dat.de2, n.chains=1, n.adapt=Nadapt, inits = initial)
  params <- c("par", "loglik", "LS")
  samps.de2 <- coda.samples(jags.de2, params, n.iter = Niter)
  name.list = names(samps.de2[[1]][1,])
  list.par.de2 = which(grepl("par", name.list)==TRUE)
  geweke.de2 = apply(samps.de2[[1]][,list.par.de2], 2, function(x) geweke.diag(x)$z)
  if (sum(abs(geweke.de2)<2)==10) break
}

total_iter <-i

name.list = names(samps.de2[[1]][1,])
list.par.de2 = which(grepl("par", name.list)==TRUE)
list.par.ll2 = which(grepl("loglik", name.list)==TRUE)
list.LS = which(grepl("LS", name.list)==TRUE)
geweke.de2 = apply(samps.de2[[1]][,list.par.de2], 2, function(x) geweke.diag(x)$z)
ll.de2 = as.matrix(samps.de2[[1]][,list.par.ll2])

source('../sources/de_WAIC_LOOCV_DIC_gmm_sameSig_G2_RV2.R')

samps.par = samps.de2[[1]][,list.par.de2]
save(samps.par, waic_c, loo_c, dic_c, w1_c, l1_c, w1_m, l1_m,
     waic_m, loo_m, dic_m, geweke.de2, total_iter, del.list,
     file=paste0("../results/N", Nsmp, "/", MD,"/", PP, "/g2_de_", filename2))


