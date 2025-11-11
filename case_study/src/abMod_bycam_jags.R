model{
  ## DETECTION PROCESS ##
  # y_mtk is a binary variable 1 = at least one observation at day k
  for(t in 1:T){
    for(m in 1:M){
      p[m,t] <- 1 - (1 - rho[m,t])**n[m, t]
      for(k in 1:K){
        y[m,t,k] ~ dcat(c(1-p[m,t], p[m,t]))
      }
    }
  }
  
  ## POPULATION PROCESS ##
  for(m in 1:M){
    n[m, 1] ~ dpois(lam)
    for(t in 2:T){
      n[m, t] <- S[m, t-1] + G[m, t-1]
      S[m, t-1] ~ dbin(tau[t-1], n[m, t-1])
      G[m, t-1] ~ dpois(gam[t-1])
    }
  }

  ## LINEAR PREDICTORS ## 
  for(t in 1:T){
    for(m in 1:M){
      # logit(rho[m,t])  <- a_RE_stat[m] + inprod(a, rho_covs[t,])
      logit(rho[m,t])  <- inprod(a, rho_covs[t,])
    }
    gam[t] <- exp(inprod(b, gam_covs[t,]))
    tau[t] <- ilogit(inprod(d, tau_covs[t,]))
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
  # for(m in 1:M){
  #   a_RE_stat[m] ~ dlogis(0,1)
  # }
  lam ~ dunif(1,100)
  
}