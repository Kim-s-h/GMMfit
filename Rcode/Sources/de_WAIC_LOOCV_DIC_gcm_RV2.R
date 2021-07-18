# ---------------------------------------------
#  WAIC, LOO-CV, DIC calculation for gcm based on conditional medians
#  G = 1
# ---------------------------------------------

S = Niter
Nsmp = nrow(resp)
J = ncol(resp)
Lambda = matrix(c(rep(1,J), (c(1:J)-1)), nrow = J)

m.tol = 0.0001   # maximum tol
m.eval = 1024L  # The maximum number of function evaluations needed

stat.samp = summary(samps.de1)$statistics

# mean statistics 
tau.m = (stat.samp[list.par.de1[6],1])
LS.h = as.numeric(stat.samp[list.LS,1])
LS.h.m = matrix(LS.h, nrow=Nsmp, ncol = 2 )

beta = as.vector(stat.samp[list.par.de1[1:2],1])
Phi = matrix(NA, nrow=2, ncol =2)
Phi[lower.tri(Phi, diag = T)] = as.vector(stat.samp[list.par.de1[3:5],1])
Phi[1,2] = Phi[2,1]


# estimates from each iteration
beta.s = samps.de1[[1]][,list.par.de1[1:2]]
Phi.s.elements = samps.de1[[1]][,list.par.de1[3:5]]
Phi.s = vector("list", S)
for(s in 1:S){
  Phi.s[[s]] = matrix(NA, nrow = 2, ncol =2)
  Phi.s[[s]][lower.tri(Phi.s[[s]], diag = T)] = as.vector(Phi.s.elements[s,])
  Phi.s[[s]][1,2] = Phi.s[[s]][2,1]
}
tau.s = as.numeric(samps.de1[[1]][,list.par.de1[6]])



# conditional likelihood from JAGS
loglik = array(NA, dim=c(Nsmp, J, S))
for(s in 1:S){
  loglik[,,s] = matrix(ll.de1[s,], nrow=Nsmp, ncol=J)
}


# Likelihood for DIC
log.lik.h = array(NA, dim=c(Nsmp, J))
for(i in 1:Nsmp){
  for(j in 1:J){
    log.lik.h[i,j] = ddexp(as.vector(resp[i,j]), 
                           location = Lambda[j,]%*%(LS.h.m[i, ]), 
                           scale = (1/tau.m), log=TRUE)
  }
}

# -------------------------------
# Conditional likelihood from JAGS
# -------------------------------
# WAIC
w1_c = waic(ll.de1)
waic_c = w1_c$estimates[3,1]
print(paste("waic_c:", round(waic_c,2)))

# LOO
l1_c = loo(ll.de1)
loo_c = l1_c$estimates[3,1]
print(paste("loo_c:", round(loo_c,2)))

# DIC
logp.h = sum(log.lik.h)
p.dic = 2*(logp.h - mean(apply(loglik, 3, sum)))
dic_c = -2*logp.h + 2*p.dic
print(paste("dic_c:", round(dic_c,2)))






# -----------------------------------
# Marginalized likelihood based method
# -----------------------------------


# INTEGRATE - vectorized (function)
f_yb <- function(x, Phi, beta, tau, resp){
  abs_sum <- apply(x, 2, function(z) sum(tau*abs(resp - Lambda %*% z)))
  exp(matrix(J*log(tau) - J*log(2) - abs_sum, ncol = ncol(x)))*matrix(dmvn(t(x), mu=beta, sigma=Phi), ncol = ncol(x))
}

# marginalized likelihood ----------------------------------
log.lik.m.s = array(NA, dim=c(S, Nsmp))
ptm <- proc.time()
for(s in 1:S){
  log.lik.m.s[s,] = unlist(mclapply(
    lapply(1:Nsmp, function(x) resp[x,]),
    function(a)
      log(hcubature(f_yb, lowerLimit = rep(-Inf,2), upperLimit = rep(Inf,2), vectorInterface = T, tol = m.tol, maxEval = m.eval, 
                    Phi = Phi.s[[s]], beta = beta.s[s,], tau = tau.s[s], resp = a)$integral)
    , mc.cores = numCores))
  if(s %% 20 == 0) print(paste0("s =", s))
}
print(proc.time() - ptm)



# marginalized likelihood for DIC (d-hat)

log.lik.m2 = unlist(mclapply(
  lapply(1:Nsmp, function(x) resp[x,]),
  function(a)
    log(hcubature(f_yb, lowerLimit = rep(-Inf,2), upperLimit = rep(Inf,2), vectorInterface = T, tol = m.tol, maxEval = m.eval, 
                  Phi = Phi, beta = beta, tau = tau.m, resp = a)$integral)
  , mc.cores = numCores))



# -------------------------------
# Marginalized likelihood
# -------------------------------
if(length(which(is.nan(log.lik.m.s)))>0){
  nan.loc = unique(floor(which(is.nan(log.lik.m.s))/S) + 1)
  del.list = c()
  for(w in nan.loc){
    del.list = c(del.list, which(is.nan(log.lik.m.s[,w])==TRUE))
  }
  del.list = unique(del.list)
  log.lik.m.s = log.lik.m.s[-del.list,]
} else {
  del.list = NA
}

# # WAIC
w1_m = waic(log.lik.m.s)
waic_m = w1_m$estimates[3,1]
print(paste("waic_m:", round(waic_m,2)))

# LOOCV
l1_m = loo(log.lik.m.s)
loo_m = l1_m$estimates[3,1]
print(paste("loo_m:", round(loo_m,2)))

# DIC
logp.h = sum(log.lik.m2)
p.dic = 2*(logp.h - mean(apply(log.lik.m.s, 1, sum)))
dic_m = -2*logp.h + 2*p.dic
print(paste("dic_m:", round(dic_m,2)))


