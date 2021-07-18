# ---------------------------------------------
#  WAIC, LOO-CV, DIC calculation for normal gmm
#  G =2
#  Traditional (normal)
# ---------------------------------------------
# library(mvtnorm)

S = Niter
y = resp
Nsmp = nrow(y)
J = ncol(y)
Lambda = matrix(c(rep(1,J), (c(1:J)-1)), nrow = J)

# Estimates for DIC
tmp = summary(samps.nm2)
stat.samp = (tmp$statistics)
sig.h = stat.samp[list.s,1]
LS.h = as.numeric(stat.samp[list.LS,1])
LS.h.m = matrix(LS.h, nrow=Nsmp, ncol = 2 )

# conditional likelihood from JAGS
loglik = array(NA, dim=c(Nsmp, J, S))
for(s in 1:S){
  loglik[,,s] = matrix(ll.nm[s,], nrow=Nsmp, ncol=J)
}


# Likelihood for DIC
log.lik.h = array(NA, dim=c(Nsmp, J))
for(i in 1:Nsmp){
  for(j in 1:J){
    log.lik.h[i,j] = dnorm(as.vector(y[i,j]), 
                           mean = Lambda[j,]%*%(LS.h.m[i, ]), 
                           sd = sqrt(sig.h), log=TRUE)
  }
}

# -------------------------------
# Conditional likelihood from JAGS
# -------------------------------
# WAIC
w1 = waic(ll.nm)
waic_c = w1$estimates[3,1]
print(paste("waic_c:", round(waic_c,2)))

# LOO
l1 = loo(ll.nm)
loo_c = l1$estimates[3,1]
print(paste("loo_c:", round(loo_c,2)))

# DIC
logp.h = sum(log.lik.h)
p.dic = 2*(logp.h - mean(apply(loglik, 3, sum)))
dic_c = -2*logp.h + 2*p.dic
print(paste("dic_c:", round(dic_c,2)))



# ---------------------------------------- #
# Marginal likelihood
# ---------------------------------------- #

# parameter estimates
beta = matrix(NA, nrow = 2, ncol = 2)
beta[,1] = as.vector(stat.samp[list.par.nm2[1:2],1])
beta[,2] = as.vector(stat.samp[list.par.nm2[3:4],1])
Phi = matrix(NA, nrow=2, ncol =2)
Phi[lower.tri(Phi, diag = T)] = as.vector(stat.samp[list.par.nm2[5:7],1])
Phi[1,2] = Phi[2,1]
pi_g = as.vector(stat.samp[list.par.nm2[8:9],1])
mem.est = as.numeric(tmp$quantiles[list.m, 3])

# for each iteration
beta.s.elements = samps.nm2[[1]][,list.par.nm2[1:4]]
Phi.s.elements = samps.nm2[[1]][,list.par.nm2[5:7]]
Phi.s = beta.s = vector("list", S)
LS.s = array(NA, dim=c(Nsmp, 2, S))
for(s in 1:S){
  beta.s[[s]] = matrix(NA, nrow = 2, ncol =2)
  beta.s[[s]][,1] = as.vector(beta.s.elements[s,c(1:2)])
  beta.s[[s]][,2] = as.vector(beta.s.elements[s,c(3:4)])
  Phi.s[[s]] = matrix(NA, nrow = 2, ncol =2)
  Phi.s[[s]][lower.tri(Phi.s[[s]], diag = T)] = as.vector(Phi.s.elements[s,])
  Phi.s[[s]][1,2] = Phi.s[[s]][2,1]
  LS.s[,,s] = as.numeric(samps.nm2[s,list.LS][[1]])
}
pi.s = samps.nm2[[1]][,list.par.nm2[8:9]]
sig.s = samps.nm2[[1]][,list.s]
mem.s = samps.nm2[[1]][,list.m]

log.lik.i.s2 = array(NA, dim=c(S, Nsmp))
for(s in 1:S){
  for(i in 1:Nsmp){
    log.lik.i.s2[s,i] = log(pi.s[s,1]*dmvnorm(y[i,],
                                mean = Lambda %*% beta.s[[s]][,1],
                                sigma = Lambda %*% Phi.s[[s]] %*% t(Lambda) + sig.s[s]*diag(J), log=FALSE) +
      pi.s[s,2]*dmvnorm(y[i,],
                        mean = Lambda %*% beta.s[[s]][,2],
                        sigma = Lambda %*% Phi.s[[s]] %*% t(Lambda) + sig.s[s]*diag(J), log=FALSE))
  }
}


# For DIC
log.lik.i2 = rep(NA, Nsmp)
for(i in 1:Nsmp){
  log.lik.i2[i] = log(pi_g[1]*dmvnorm(y[i,],
                          mean = Lambda %*% beta[,1],
                          sigma = Lambda %*% Phi %*% t(Lambda) + sig.h*diag(J), log=FALSE) +
    pi_g[2]*dmvnorm(y[i,],
                    mean = Lambda %*% beta[,2],
                    sigma = Lambda %*% Phi %*% t(Lambda) + sig.h*diag(J), log=FALSE))
}



# WAIC
w2_m = waic(log.lik.i.s2)
waic_m = w2_m$estimates[3,1]
print(paste("waic_m:", round(waic_m,2)))

# LOOCV
l2_m = loo(log.lik.i.s2)
loo_m = l2_m$estimates[3,1]
print(paste("loo_m:", round(loo_m,2)))

# DIC
logp.h = sum(log.lik.i2)
p.dic = 2*(logp.h - mean(apply(log.lik.i.s2, 1, sum)))
dic_m = -2*logp.h + 2*p.dic
print(paste("dic_m:", round(dic_m,2)))




