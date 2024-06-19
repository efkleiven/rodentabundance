library(tidyverse)
library(foreach)
library(lubridate)

library(LaplacesDemon)
library(runjags)
library(rjags)
library(coda)

rm(list = ls())

modname = ""
K = 7
SITE = "Valdres"

### 1. Load Data ---------------------------------------------------------------
  
ct_data.filename <- "data/camData.rds"
snow_data.filename <- "data/snowData.rds"

ct_data <- readRDS(ct_data.filename) %>%
  filter(site == SITE) %>%
  left_join(readRDS(snow_data.filename), by = c("station", "date"))%>%
  filter(year >= 2017)

day.orig <- min(ct_data$date)

ct_data <- ct_data %>%
  mutate(jday = julian(date, origin = day.orig) + 1,
         t = (jday - 1) %/% K + 1)

detHist <- ct_data %>%
  group_by(station, t, date, jday) %>% 
  summarize(Microtus = length(which(species == "microtus")))

### 2. Get Data Matrix ---------------------------------------------------------

y.df <- detHist%>%
  group_by(station, t) %>%
  mutate(k = 1:n()) %>% 
  pivot_wider(names_from = station, values_from = Microtus) 

M <- length(unique(detHist$station))
T <- max(detHist$t)

y.mat <- array(NA, dim = c(M, T, K))
for(i in 1:nrow(y.df)){
  y.mat[, y.df$t[i], y.df$k[i]] <- t(y.df[i,unique(detHist$station)])
}

### 3. Get covariates ----------------------------------------------------------

covs <- ct_data %>% 
  group_by(t) %>%
  summarize(year = year[1],
            date = date[1],
            yday = yday(date),
            cos1 = cos(2*pi*yday/365),
            cos2 = cos(4*pi*yday/365),
            sin1 = sin(2*pi*yday/365),
            sin2 = sin(4*pi*yday/365),
            int = 1,
            zeros = 0,
            SNOW = as.numeric(any(SNOW,  na.rm = T)))%>%
  mutate(breeding = as.numeric(month(date) %in% 5:9)) %>%
  left_join(data.frame(t = 1:T), .)

years <- covs$year
Nyears <- length(unique(years))

### 4. Prepare to run the model ------------------------------------------------

# individual detection
rho_covs <- as.matrix(covs[, c("SNOW")])

# natality
gam_covs <- as.matrix(covs[, c("int", "breeding")])

#survival rate
tau_covs <- as.matrix(covs[, c("int", "SNOW")])

rho_covs[is.na(rho_covs)] <- 0 
gam_covs[is.na(gam_covs)] <- 0 
tau_covs[is.na(tau_covs)] <- 0 

data_list <- list(rho_covs = rho_covs,
                  gam_covs = gam_covs,
                  tau_covs = tau_covs,
                  Ncov_rho = ncol(rho_covs),
                  Ncov_gam = ncol(gam_covs),
                  Ncov_tau = ncol(tau_covs),
                  T = T, M = M, K = K,
                  y = y.mat)


### 5. Run Model ---------------------------------------------------------------

NCHAINS = 2

inits <- foreach(ch = 1:NCHAINS) %do%{
  list(a = rnorm(ncol(rho_covs), 0, 0),
       b = rnorm(ncol(gam_covs), 0, 0),
       d = rnorm(ncol(tau_covs), 0, 0), 
       lam = runif(1, 1, 10),
       a_RE_stat = rnorm(M, 0, 0))
}

Mod <- run.jags(model = "src/abModPoiss_jags.R",
              monitor = c("a", "b", "d", "lam",
                          "a_RE_stat",
                          "n", "G", "S",
                          "Chi2.obs", "Chi2.sim"),
              data = data_list,
              n.chains = NCHAINS,
              inits = inits,
              adapt = 500,
              burnin = 2000,
              sample = 2000,
              thin = 1,
              method = "parallel")

# saveRDS(Mod, file = paste0("outputs/Models/abModPoiss_", SITE,"_",modname, ".rds"))
# saveRDS(covs, file = paste0("outputs/covsMat/", SITE,"_K", K, ".rds"))

### 6. Model Analysis ----------------------------------------------------------

M.mat_ <- as.matrix(as.mcmc.list(Mod))

M.mat <- M.mat_[, c(-grep(c("Chi2.sim"), colnames(M.mat_)),
                    -grep(c("Chi2.obs"), colnames(M.mat_)))]

P.V <- M.mat_[, c(grep(c("Chi2.sim"), colnames(M.mat_)),
                  grep(c("Chi2.obs"), colnames(M.mat_)))] %>% 
  apply(1, function(x){as.numeric(x[2]>x[1])})%>%
  mean

L.O.F <- M.mat_[, c(grep(c("Chi2.sim"), colnames(M.mat_)),
                    grep(c("Chi2.obs"), colnames(M.mat_)))] %>% 
  apply(1, function(x){as.numeric(x[2]/x[1])})%>%
  mean


M.mat <- M.mat_[, c(-grep(c("n"), colnames(M.mat_)),
                    -grep(c("G"), colnames(M.mat_)),
                    -grep(c("S"), colnames(M.mat_)))]

n.qu <- M.mat_[, grep("n", colnames(M.mat_))] %>%
  apply(2, function(x){quantile(x, c(0.025,0.5,0.975))})

n.df <- data.frame(date = covs$date,
                          Nest.inf = rep(n.qu[1,]),
                          Nest.med = rep(n.qu[2,]),
                          Nest.sup = rep(n.qu[3,]))

ggplot(n.df)+
  geom_line(data = detHist, aes(x=date, y = Microtus, col = station))+
  geom_line(aes(x=date, y = Nest.med), col = "black")+
  geom_ribbon(aes(x = date, ymin = Nest.inf, ymax = Nest.sup),
              fill = "black", alpha = 0.5)+
  geom_rect(data = covs, aes(xmin = date, xmax = lag(date, 1), 
                             ymin = -3, ymax = -2, fill = factor(SNOW)))+
  scale_fill_manual(values = c("orange", "lightblue"), name = "Snow", labels = c("NO", "YES"))+
  scale_x_date(breaks = "6 months")+
  theme_bw()