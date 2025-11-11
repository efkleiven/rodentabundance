model{
  ## DETECTION PROCESS ##
  # y_mtk is a binary variable 1 = at least one observation at day k
  for(t in 1:T){
    for(m in 1:M){
      p[m,t] <- 1 - (1 - rho[m,t])**n[blocks[m], t]
      for(k in 1:K){
        y[m,t,k] ~ dcat(c(1-p[m,t], p[m,t]))
      }
    }
  }
  
  ## POPULATION PROCESS ##
  for(b in 1:B){
    n[b, 1] ~ dpois(lam)
    for(t in 2:T){
      n[b, t] <- S[b, t-1] + G[b, t-1]
      S[b, t-1] ~ dbin(tau[t-1], n[b, t-1])
      G[b, t-1] ~ dpois(gam[t-1])
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