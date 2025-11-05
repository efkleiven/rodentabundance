library(tidyverse)
library(runjags)
library(rjags)
library(coda)

rm(list = ls())

nYears <- 5
nWeeks <- 52 * nYears
nRegions <- 10
nSites <- 10
nOccas <- 7

winterStart <- 41
winterLength <- 26

# "True" values
lam <- 5
omega <- .8
gamma <- 1

#Simulate true abundances, N, for each location

N <- matrix(NA, nRegions, nWeeks)
S <- G <- matrix(NA, nRegions, nWeeks-1)

#First year of sampling follows a Poisson distribution
N[,1] <- rpois(nRegions, lam)

#subsequent years follow the birth-death-immigration process
for(t in 2:nWeeks) {
  S[,t-1] <- rbinom(nSites, N[,t-1], omega)
  G[,t-1] <- rpois(nSites, gamma)
  N[,t] <- S[,t-1] + G[,t-1] 
}

# Generate data vector y for the counts

nReps <- 7
p <- 0.3

y <- array(NA, c(nRegions, nWeeks, nSites, nOccas))
for(t in 1:nWeeks) {
  for(j in 1:nSites) {
    for(k in 1:nOccas) {
      y[,t,j, k] <- sign(rbinom(nSites, N[,t], p))
    }
  }
}

plot(apply(y[ , , , 1], 2,mean))

N.df <- data.frame(t(N))
regions <- colnames(N.df)

N.df <- N.df %>%
  mutate(t = 1:n())%>%
  pivot_longer(cols = all_of(regions),
               names_to = "region", values_to = "N")

ggplot(N.df) + 
  geom_line(aes(x = t, y = N, col = region, group = region), show.legend = F)

data.list <- list(y = y,
                  T = nWeeks,
                  R = nRegions,
                  M = nSites,
                  K = nOccas)

Ni <- y[,,1,1] + 5
Ni[, -1] <- NA
Si <- S
Si[] <- 4
Gi<-matrix(rep(5, 10),nrow=nRegions,ncol=nWeeks-1, byrow = T)

inits <- map(1:2, ~ list(gamma = .8, omega= .8, rho = .3, lam = 5))

Mod_base <- run.jags(model = "abMod_jags_base.R",
                monitor = c( "gamma", "omega", "rho", "lam", "n"),
                data = data.list,
                n.chains = 4,
                inits = inits,
                adapt = 1000,
                burnin = 1000,
                sample = 1000,
                thin = 2,
                method = "parallel")

mod_base <- as.matrix(as.mcmc.list(Mod_base))

save(mod_base, file="sim_output_base.RData")

#boxplot(exp(mod[, "gam"]))
#boxplot(mod[, "omega"])
#boxplot(mod[, "lam"])
#boxplot(mod[, "rho"])

#-end of script