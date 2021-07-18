# ---------------------------------------------
#  WAIC, LOO-CV, DIC calculation for normal gcm
#  G =1
#  Traditional (normal)
# ---------------------------------------------
library(mvtnorm)

S = Niter - burnIn
y = resp
Nsmp = nrow(y)
J = ncol(y)
Lambda = matrix(c(rep(1,J), (c(1:J)-1)), nrow = J)

tmp = summary(window(samps.nm1, start = (burnIn+1)))
stat.samp = (tmp$statistics)

# mean sig2 for DIC
sig.h = stat.samp[list.s,1]

# mean LS for DIC
LS.h = as.numeric(stat.samp[list.LS,1])
LS.h.m = matrix(LS.h, nrow=Nsmp, ncol = 2 )


# conditional likelihood from JAGS
loglik = array(NA, dim=c(Nsmp, J, S))
for(s in 1:S){
  loglik[,,s] = matrix(ll.nm1[s,], nrow=Nsmp, ncol=J)
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
w1 = waic(ll.nm1)
waic_c = w1$estimates[3,1]
print(paste("waic_c:", round(waic_c,2)))

# LOO
l1 = loo(ll.nm1)
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
beta = as.vector(stat.samp[list.par.nm1[1:2],1])
Phi = matrix(NA, nrow=2, ncol =2)
Phi[lower.tri(Phi, diag = T)] = as.vector(stat.samp[list.par.nm1[3:5],1])
Phi[1,2] = Phi[2,1]

# for each iteration
beta.s = samps.nm1[[1]][((burnIn+1):Niter),list.par.nm1[1:2]]
Phi.s.elements = samps.nm1[[1]][((burnIn+1):Niter),list.par.nm1[3:5]]
Phi.s = vector("list", S)
LS.s = array(NA, dim=c(Nsmp, 2, S))
for(s in 1:S){
  Phi.s[[s]] = matrix(NA, nrow = 2, ncol =2)
  Phi.s[[s]][lower.tri(Phi.s[[s]], diag = T)] = as.vector(Phi.s.elements[s,])
  Phi.s[[s]][1,2] = Phi.s[[s]][2,1]
  LS.s[,,s] = as.numeric(samps.nm1[(burnIn+s),list.LS][[1]])
}
sig.s = samps.nm1[[1]][((burnIn+1):Niter),list.s]

log.lik.i.s = array(NA, dim=c(S, Nsmp))
for(s in 1:S){
  for(i in 1:Nsmp){
    log.lik.i.s[s,i] = dmvnorm(y[i,], 
                               mean = Lambda %*% beta.s[s,], 
                               sigma = Lambda %*% Phi.s[[s]] %*% t(Lambda) + sig.s[s]*diag(J), log=TRUE)
  }
}

# For DIC
log.lik.i = rep(NA, Nsmp)
for(i in 1:Nsmp){
  log.lik.i[i] = dmvnorm(y[i,], 
                         mean = Lambda %*% beta, 
                         sigma = Lambda %*% Phi %*% t(Lambda) + sig.h*diag(J), log=TRUE)
}





# # WAIC
w1_m = waic(log.lik.i.s)
waic_m = w1_m$estimates[3,1]
print(paste("waic_m:", round(waic_m,2)))

# LOOCV
l1_m = loo(log.lik.i.s)
loo_m = l1_m$estimates[3,1]
print(paste("loo_m:", round(loo_m,2)))

# DIC
logp.h = sum(log.lik.i)
p.dic = 2*(logp.h - mean(apply(log.lik.i.s, 1, sum)))
dic_m = -2*logp.h + 2*p.dic
print(paste("dic_m:", round(dic_m,2)))


