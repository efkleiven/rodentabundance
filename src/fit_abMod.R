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
SITE = "Hakoya"

### 1. Load Data ---------------------------------------------------------------

ct_data.filename <- "data/camData.rds"
snow_data.filename <- "data/snowData.rds"

ct_data <- readRDS(ct_data.filename) %>%
  filter(site == SITE) %>%
  left_join(readRDS(snow_data.filename), by = c("station", "date")) %>%
  filter(year >= 2017)

day.orig <- min(ct_data$date)

ct_data <- ct_data %>%
  mutate(jday = julian(date, origin = day.orig) + 1,
         t = (jday - 1) %/% K + 1)

detHist <- ct_data %>%
  group_by(station, t, date, jday) %>% 
  summarize(Microtus = as.numeric(any(species == "microtus")))

### 2. Get Data Matrix ---------------------------------------------------------

y.df <- detHist%>%
  group_by(station, t) %>%
  mutate(k = 1:n()) %>% 
  pivot_wider(names_from = station, values_from = Microtus) 

M <- length(unique(detHist$station))
T <- max(detHist$t)

y.mat <- array(NA, dim = c(M, T, K))
for(i in 1:nrow(y.df)){
  y.mat[, y.df$t[i], y.df$k[i]] <- t(y.df[i,unique(detHist$station)]) + 1
}

### 3. Get covariates ----------------------------------------------------------

covs <- ct_data %>% 
  group_by(t) %>%
  summarize(year = year[1],
            date = date[1],
            yday = yday(date),
            yday2 = yday**2,
            cos1 = cos(2*pi*yday/365),
            cos2 = cos(4*pi*yday/365),
            sin1 = sin(2*pi*yday/365),
            sin2 = sin(4*pi*yday/365),
            int = 1, 
            zeros = 0,
            SNOW = as.numeric(any(SNOW,  na.rm = T))) %>%
  mutate(breeding = as.numeric(month(date) %in% 5:9)) %>%
  left_join(data.frame(t = 1:T), .)

years <- covs$year
Nyears <- length(unique(years))

### 4. Prepare to run the model ------------------------------------------------

# individual detection
rho_covs <- as.matrix(covs[, c("SNOW")])

# natality
gam_covs <- as.matrix(covs[, c("int")])

#survival rate
tau_covs <- as.matrix(covs[, c("int", "yday", "yday2")])

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

Mod <- run.jags(model = "src/abMod_jags.R",
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

# saveRDS(Mod, file = paste0("outputs/Models/abMod_", SITE,"_",modname, ".rds"))

### 6. Model Analysis ----------------------------------------------------------

M.mat_ <- as.matrix(as.mcmc.list(Mod))

source("src/model_summary.R")

# M.mat <- M.mat_[, c(-grep(c("Chi2.sim"), colnames(M.mat_)),
#                     -grep(c("Chi2.obs"), colnames(M.mat_)))]
# 
# P.V <- M.mat_[, c(grep(c("Chi2.sim"), colnames(M.mat_)),
#                     grep(c("Chi2.obs"), colnames(M.mat_)))] %>% 
#   apply(1, function(x){as.numeric(x[2]>x[1])})%>%
#   mean
# 
# L.O.F <- M.mat_[, c(grep(c("Chi2.sim"), colnames(M.mat_)),
#                   grep(c("Chi2.obs"), colnames(M.mat_)))] %>% 
#   apply(1, function(x){as.numeric(x[2]/x[1])})%>%
#   mean

M.mat <- M.mat_[, c(-grep(c("Chi2.obs"), colnames(M.mat_)),
                    -grep(c("Chi2.sim"), colnames(M.mat_)))]

sumStats <- mod.sumStat(mm = M.mat, data_list = data_list)

M.mat <- M.mat[, c(-grep(c("n"), colnames(M.mat_)),
                    -grep(c("G"), colnames(M.mat_)),
                    -grep(c("S"), colnames(M.mat_)))]

n.qu <- M.mat_[, grep("n", colnames(M.mat_))] %>%
  apply(2, function(x){quantile(x, c(0.025,0.5,0.975), na.rm = T)})

n.df <- data.frame(date = covs$date,
                   Nest.inf = rep(n.qu[1,]),
                   Nest.med = rep(n.qu[2,]),
                   Nest.sup = rep(n.qu[3,])) %>%
  mutate(year = year(date),
         month = month(date))

ggplot(n.df)+
  geom_line(data = detHist %>% group_by(t, station) %>% summarize(Microtus = sum(Microtus), date = date[1]),
            aes(x=date, y = Microtus, col = station))+
  geom_line(aes(x=date, y = Nest.med), col = "black")+
  geom_ribbon(aes(x = date, ymin = Nest.inf, ymax = Nest.sup),
              fill = "black", alpha = 0.5)+
  geom_rect(data = covs, aes(xmin = date, xmax = lag(date, 1), 
                             ymin = -1.5, ymax = -.5, fill = factor(SNOW)))+
  scale_fill_manual(values = c("orange", "lightblue"), name = "Snow", labels = c("NO", "YES"))+
  scale_x_date(breaks = "6 months")+
  theme_bw()


date.df <- data.frame(date = seq(as.Date("2019/1/1"), as.Date("2019/12/31"), by = "day"),
                      int = 1,
                      yday = 1:365) %>%
  mutate(yday2 = yday ** 2)

date.df$omega <- M.mat[, grep("d\\[", colnames(M.mat_))]%>%
  apply(2, median) %*% t(as.matrix(date.df[,c("int", "yday", "yday2")])) %>%
  c %>%
  invlogit

date.df %>% 
ggplot()+
  geom_line(aes(x = date, y = omega)) +
  ggtitle("survival rate") +
  theme_bw()

a_int <- M.mat[, grep("a_RE_stat\\[", colnames(M.mat_))]%>%
  apply(1, mean)
a_SNOW <- M.mat[,"a"]

a.df <- data.frame(NOSNOW = invlogit(a_int),
                   SNOW = invlogit(a_int+a_SNOW))

ggplot(a.df)+
  geom_violin(aes(x = "NOSNOW", y = NOSNOW), fill = "lightgreen")+
  geom_violin(aes(x = "SNOW", y = SNOW), fill = "lightblue")+
  geom_point(aes(x = "NOSNOW", y = NOSNOW), col = "red", size = .5)+
  geom_point(aes(x = "SNOW", y = SNOW), col = "red", size = .5)+
  geom_segment(aes(x = "NOSNOW", xend = "SNOW", y = NOSNOW, yend = SNOW), size = .1) +
  ggtitle("individual detection probability") +
  theme_bw()

d_int <- M.mat[, "d[1]"]
d_SNOW <- M.mat[,"d[2]"]

d.df <- data.frame(NOSNOW = invlogit(d_int),
                   SNOW = invlogit(d_int+d_SNOW))

ggplot(d.df)+
  geom_violin(aes(x = "NOSNOW", y = NOSNOW), fill = "lightgreen")+
  geom_violin(aes(x = "SNOW", y = SNOW), fill = "lightblue")+
  geom_point(aes(x = "NOSNOW", y = NOSNOW), col = "red", size = .5)+
  geom_point(aes(x = "SNOW", y = SNOW), col = "red", size = .5)+
  geom_segment(aes(x = "NOSNOW", xend = "SNOW", y = NOSNOW, yend = SNOW), size = .1) +
  ggtitle("survival rate") +
  theme_bw()