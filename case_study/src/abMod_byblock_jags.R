model{
  ## DETECTION PROCESS ##
  for(t in 1:T){
    for(m in 1:M){
      p[m,t] <- 1 - (1 - theta)**n[region[m], t]
      for(k in 1:K){
        y[m,t,k] ~ dcat(c(1-p[m,t], p[m,t]))
      }
    }
  }
  
  ## POPULATION PROCESS ##
  for(r in 1:R){
    n[r, 1] ~ dpois(lambda)
    for(t in 2:T){
      n[r, t] <- S[r, t-1] + G[r, t-1]
      S[r, t-1] ~ dbin(omega[t-1], n[r, t-1])
      G[r, t-1] ~ dpois(gamma[t-1])
    }
  }

  ## LINEAR PREDICTORS ## 
  logit(theta)  <- mu_theta
  for(t in 1:T){
    logit(omega[t]) <- mu_omega + beta_omega_year[year[t]] + beta_omega_season[season[t]]
    log(gamma[t]) <- mu_gamma + beta_gamma_year[year[t]] + beta_gamma_season[season[t]]
  }

  ## PRIORS ##
  mu_gamma ~ dlogis(0,1)
  mu_omega ~ dlogis(0,1)
  mu_theta ~ dlogis(0,1)
  lambda ~ dunif(1,100)
  
  beta_gamma_year[1] ~ dnorm(0,1e09)
  beta_omega_year[1] ~ dnorm(0,1e09)
  for(year in 2:Nyear){
    beta_gamma_year[year] ~ dlogis(0,1)
    beta_omega_year[year] ~ dlogis(0,1)
  }
  
  beta_gamma_season[1] ~ dnorm(0,1e09)
  beta_gamma_season[2] ~ dlogis(0,1)
  beta_omega_season[1] ~ dnorm(0,1e09)
  beta_omega_season[2] ~ dlogis(0,1)
  
}