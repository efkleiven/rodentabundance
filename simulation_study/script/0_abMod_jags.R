model{
  ## DETECTION PROCESS ##
  # y_mtk is a binary variable 1 = at least one observation at day k
  for(s in 1:R){
    for(t in 1:T){
      p[s, t] <- 1 - (1 - rho)**n[s, t]
      for(m in 1:M){
        for(k in 1:K){
          y[s,t,m,k] ~ dbern(p[s, t])
        }
      }
    }
  }

  
  ## POPULATION PROCESS ##
  for(s in 1:R){
    n[s, 1] ~ dpois(lam)
    for(t in 2:T){
      n[s, t] <- S[s, t-1] + G[s, t-1]
      S[s, t-1] ~ dbin(omega, n[s, t-1])
      G[s, t-1] ~ dpois(gamma[t-1])
    }
  }

  for(t in 1:T){
    gamma[t] <- exp(gam) * breeding[t]
  }
  
  ## PRIORS ##
  lam ~ dunif(1,100)
  gam ~ dnorm(0,1)
  omega ~ dlogis(0,1)
  rho ~ dlogis(0,1)
}
