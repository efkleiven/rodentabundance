model{
  ## DETECTION PROCESS ##
  # y_mtk = number of pictures at day k
  for(t in 1:T){
    for(m in 1:M){
      for(k in 1:K){
        y[m,t,k] ~ dpois(rho[m,t] * n[t])
        y.sim[m,t,k] ~ dpois(rho[m,t] * n[t])
      }
    }
  }
  
  ## POPULATION PROCESS ##
  n0 ~ dpois(lam)
  n[1] <- n0
  for(t in 2:T){
    n[t] <- S[t-1] + G[t-1]
    S[t-1] ~ dbin(tau[t-1], n[t-1])
    G[t-1] ~ dpois(gam[t-1])
  }
  
  ## LINEAR PREDICTORS ## 
  for(t in 1:T){
    for(m in 1:M){
      log(rho[m,t])  <- a_RE_stat[m] + inprod(a, rho_covs[t,])
    }
    log(gam[t]) <- inprod(b, gam_covs[t,])
    logit(tau[t]) <- inprod(d, tau_covs[t,])
  }

  ## PRIORS ##
  for(cov in 1:Ncov_rho){
    a[cov] ~ dlogis(0,1)
  }
  for(cov in 1:Ncov_gam){
    b[cov] ~ dlogis(0,1)
  }
  for(cov in 1:Ncov_tau){
    d[cov] ~ dlogis(0,1)
  }
  for(m in 1:M){
    a_RE_stat[m] ~ dlogis(0,1)
  }
  lam ~ dunif(1,100)
  
  for(m in 1:M){
    for(t in 1:T){
      for(k in 1:K){
        chi2[m,t,k] <- pow(y[m,t,k] - rho[m,t] * n[t], 2)/(rho[m,t] * n[t] + .001)
        chi2.sim[m,t,k] <- pow(y.sim[m,t,k] - rho[m,t] * n[t], 2)/(rho[m,t] * n[t] + .001)
      }
    }
  }
  Chi2.obs <- sum(chi2)
  Chi2.sim <- sum(chi2.sim)
}