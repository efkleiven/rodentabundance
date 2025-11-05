model{
  ## DETECTION PROCESS ##
  # y_mtk is a binary variable 1 = at least one observation at day k
  for(t in 1:T){
    for(m in 1:M){
      p[m,t] <- 1 - (1 - rho[m,t])**n[t]
      for(k in 1:K){
        y[m,t,k] ~ dcat(c(1-p[m,t], p[m,t]))
        y.sim[m,t,k] ~ dcat(c(1-p[m,t], p[m,t]))
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
      logit(rho[m,t])  <- a_RE_stat[m] + inprod(a, rho_covs[t,])
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
  for(m in 1:M){
    a_RE_stat[m] ~ dlogis(0,1)
  }
  lam ~ dunif(1,100)
  

  ## GET DISCREPANCIES ##
  
  for(m in 1:M){
    for(t in 1:T){
      chi2[m,t] <- (sum(y[m,t,]) - K*p[m,t])**2/sqrt((K*p[m,t]+.001)*K*(1-p[m,t]+0.001))
      chi2.sim[m,t] <- (sum(y.sim[m,t,]) - K*p[m,t])**2/sqrt((K*p[m,t]+0.001)*K*(1-p[m,t]+0.001))
    }
  }
  Chi2.obs <- sum(chi2)
  Chi2.sim <- sum(chi2.sim)
}