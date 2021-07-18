# ---------------------------------------------
#  WAIC, LOO-CV, DIC calculation for t-distributed gcm
#  G = 1
# ---------------------------------------------
# library(mvtnorm)
# library(metRology)

S = Niter
y = resp
Nsmp = nrow(y)
J = ncol(y)
Lambda = matrix(c(rep(1,J), (c(1:J)-1)), nrow = J)

stat.samp = summary(samps.td1)$statistics

# mean statistics for DIC
sig.h = stat.samp[list.par.td1[6],1]
df.h = stat.samp[list.par.td1[7],1]
LS.h = as.numeric(stat.samp[list.LS,1])
LS.h.m = matrix(LS.h, nrow=Nsmp, ncol = 2 )
beta = as.vector(stat.samp[list.par.td1[1:2],1])
Phi = matrix(NA, nrow=2, ncol =2)
Phi[lower.tri(Phi, diag = T)] = as.vector(stat.samp[list.par.td1[3:5],1])
Phi[1,2] = Phi[2,1]


# conditional likelihood from JAGS
loglik = array(NA, dim=c(Nsmp, J, S))
for(s in 1:S){
  loglik[,,s] = matrix(ll.td1[s,], nrow=Nsmp, ncol=J)
}

# conditonal likelihood for DIC
log.lik.h = array(NA, dim=c(Nsmp, J))
for(i in 1:Nsmp){
  for(j in 1:J){
    log.lik.h[i,j] = dt.scaled(as.vector(y[i,j]),
                               df = df.h,
                               mean = Lambda[j,]%*%(LS.h.m[i, ]), 
                               sd = sqrt(sig.h), log=TRUE)
  }
}


# -------------------------------
# Conditional likelihood from JAGS
# -------------------------------

# WAIC
w1 = waic(ll.td1)
waic_c = w1$estimates[3,1]
print(paste("waic_c:", round(waic_c,2)))

# LOO
l1 = loo(ll.td1)
loo_c = l1$estimates[3,1]
print(paste("loo_c:", round(loo_c,2)))

# DIC
logp.h = sum(log.lik.h)
p.dic = 2*(logp.h - mean(apply(loglik, 3, sum)))
dic_c = -2*logp.h + 2*p.dic
print(paste("dic_c:", round(dic_c,2)))






# -----------------------------------
# Marginalized likelihood based method
# -----------------------------------

# estimates from each iteration
beta.s = samps.td1[[1]][,list.par.td1[1:2]]
LS.s = array(NA, dim=c(Nsmp, 2, S))
Phi.s.elements = samps.td1[[1]][,list.par.td1[3:5]]
Phi.s = vector("list", S)
for(s in 1:S){
  LS.s = matrix(samps.td1[[1]][s,list.LS], nrow=Nsmp, ncol = 2)
  Phi.s[[s]] = matrix(NA, nrow = 2, ncol =2)
  Phi.s[[s]][lower.tri(Phi.s[[s]], diag = T)] = as.vector(Phi.s.elements[s,])
  Phi.s[[s]][1,2] = Phi.s[[s]][2,1]
}
sig.s = as.numeric(samps.td1[[1]][,list.par.td1[6]])
d.s = as.numeric(samps.td1[[1]][,list.par.td1[7]])


# INTEGRATE - vectorized (function)
f_yw <- function(x, resp, beta, Phi, sig, dof){
  d1 <- 0.5*dof
  distval <- sapply(x, function(z) stats::mahalanobis(resp, center=Lambda%*%beta, cov = Lambda %*% Phi %*% t(Lambda) + (sig/z)*diag(J)))
  logdet <- sapply(x, function(z) sum(log(eigen(Lambda %*% Phi %*% t(Lambda) + (sig/z)*diag(J), symmetric = TRUE, only.value = TRUE)$value)))
  exp(-(J*log(2*pi) + logdet +distval)/2 + d1*log(d1) - log(gamma(d1)) + (d1-1)*log(x) - d1*x)
}


# marginalized likelihood ----------------------------------
log.lik.m.s = array(NA, dim=c(S, Nsmp))
ptm <- proc.time()
for(s in 1:S){
  log.lik.m.s[s,] = unlist(mclapply(
    lapply(1:Nsmp, function(x) resp[x,]),
    function(a)
      log(integrate(f_yw, lower = 0, upper = Inf, stop.on.error = FALSE, resp = a, beta = beta.s[s,], Phi = Phi.s[[s]], sig = sig.s[s], dof = d.s[s])$value)
    , mc.cores = numCores))
  if(s %% 50 == 0) print(paste0("s =", s))
}
print(proc.time() - ptm)


# marglinalized likelihood for DIC (d-hat) -------------------------
log.lik.m2 = rep(NA, Nsmp)
log.lik.m2 = unlist(mclapply(
  lapply(1:Nsmp, function(x) resp[x,]),
  function(a)
    log(integrate(f_yw, lower = 0, upper = Inf, stop.on.error = FALSE, resp = a, beta = beta, Phi = Phi, sig = sig.h, dof = df.h)$value)
  , mc.cores = numCores))


# -------------------------------
# WAIC, LOO, DIC, AICM using marginalized likelihood
# -------------------------------

if(length(which(is.nan(log.lik.m.s)))>0){
  nan.loc = unique(floor(which(is.nan(log.lik.m.s))/S) + 1)
  del.list = c()
  for(w in nan.loc){
    del.list = c(del.list, which(is.nan(log.lik.m.s[,w])==T))
  }
  del.list = unique(del.list)
  log.lik.m.s = log.lik.m.s[-del.list,]
} else {
  del.list = NA
}

# # WAIC
w1_ver2 = waic(log.lik.m.s)
waic_m = w1_ver2$estimates[3,1]
print(paste("waic_m:", round(waic_m,2)))

# LOOCV
l1_ver2 = loo(log.lik.m.s)
loo_m = l1_ver2$estimates[3,1]
print(paste("loo_m:", round(loo_m,2)))

# DIC
logp.h = sum(log.lik.m2)
p.dic = 2*(logp.h - mean(apply(log.lik.m.s, 1, sum)))
dic_m = -2*logp.h + 2*p.dic
print(paste("dic_m:", round(dic_m,2)))

