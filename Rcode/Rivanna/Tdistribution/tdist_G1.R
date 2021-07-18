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
# gcm (G=1)
# ----------------------------------------
dat.td1 <- list("N" = Nsmp, "y" = resp, "Time" = ncol(resp))
for(i in 1:Nstart){
  seed.jags = seed*i
  initial = list(".RNG.name" = "base::Wichmann-Hill", ".RNG.seed" = seed.jags)
  jags.td1 <- jags.model( file = "../JAGS/model_gcm_td.txt",
                          data=dat.td1, n.chains=1, n.adapt=Nadapt, inits = initial)
  params <- c("par", "loglik", "LS")
  samps.td1 <- coda.samples(jags.td1, params, n.iter = Niter)
  name.list = names(samps.td1[[1]][1,])
  list.par.td1 = which(grepl("par", name.list)==TRUE)
  geweke.td1 = apply(samps.td1[[1]][,list.par.td1], 2, function(x) geweke.diag(x)$z)
  if (sum(abs(geweke.td1)<2)==7) break
}

total_iter <-i

name.list = names(samps.td1[[1]][1,])
list.par.td1 = which(grepl("par", name.list)==TRUE)
list.par.ll1 = which(grepl("loglik", name.list)==TRUE)
list.LS = which(grepl("LS", name.list)==TRUE)
geweke.td1 = apply(samps.td1[[1]][,list.par.td1], 2, function(x) geweke.diag(x)$z)
ll.td1 = as.matrix(samps.td1[[1]][,list.par.ll1])

source('../sources/td_WAIC_LOOCV_DIC_gcm_ver3_RV.R')

samps.par = samps.td1[[1]][,list.par.td1]
save(samps.par, waic_c, loo_c, dic_c, w1, l1, w1_ver2, l1_ver2,
     waic_m, loo_m, dic_m, geweke.td1, total_iter,  del.list,
     file=paste0("../results/N", Nsmp, "/", MD,"/", PP, "/g1_td_", filename2))
