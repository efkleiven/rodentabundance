library(tidyverse)
library(lubridate)

library(LaplacesDemon)
library(runjags)
library(rjags)
library(coda)
library(patchwork)

rm(list = ls())

# POOLING = "all"
POOLING = "blocks"
# POOLING = "none"

K = 5

NCHAINS = 2

### 1. Load Data ---------------------------------------------------------------

ct_data.filename <- "case_study/data/camData_pors.rds"
regions.filename <- "case_study/data/Porsanger/camera_blocks.csv"

ct_data <- readRDS(ct_data.filename) %>%
  filter(year >= 2017)

ct_data <- ct_data %>%
  mutate(jday = julian(date, origin = min(ct_data$date)) + 1,
         t = (jday - 1) %/% K + 1)

regions <- read.csv(regions.filename, sep = ";") %>%
  arrange(site)

if(POOLING == "all"){
  regions$block <- 1
}
if(POOLING == "none"){
  regions$block <- 1:nrow(regions)
}

### 2. Get detection history ---------------------------------------------------

detHist <- ct_data %>%
  mutate(m = as.numeric(as.factor(site)),
         season = ifelse(yday(date) < 274 & yday(date) > 152 , "spring", "winter")) %>% # spring between June 1st and October 1st
  group_by(site) %>%
  mutate(year = (data.table::rleid(season) - 1) %/% 2) %>%
  filter(year < 6)

detHist <- detHist %>%
  group_by(site, m, t, date, season, year) %>% 
  summarize(vole_presence = as.numeric(any(species == "vole")),
            vole_count = sum(as.numeric(species == "vole"))) %>% 
  group_by(site, m, t) %>%
  mutate(k = 1:n()) %>%
  ungroup

### 3. Get data matrix ---------------------------------------------------------

M <- length(unique(detHist$site))
T <- max(detHist$t)

y <- expand.grid(t = 1:T, m = 1:M, k = 1:K) %>%
  left_join(detHist) %>%
  select(m,t,k,vole_presence) %>%
  mutate(vole_presence = replace_na(vole_presence, 0))

y <- array(y$vole_presence, dim = c(T, M, K)) %>%
  aperm(c(2,1,3))

y <- y

### 3. Get covariates ----------------------------------------------------------

time_covariates <- detHist %>% 
  group_by(t) %>%
  summarize(season = season[1],
            year = year[1]) %>%
  distinct

### 4. Set up model ------------------------------------------------------------

data_list <- list(y = y,
                  region = regions$block,
                  year = as.numeric(as_factor(time_covariates$year)),
                  season = as.numeric(as_factor(time_covariates$season)),
                  T = T, M = M, K = K, 
                  R = length(unique(regions$block)),
                  Nyear = length(unique(time_covariates$year)))

inits <- map(1:NCHAINS, 
             ~ list(mu_gamma = rnorm(1, 0, 0),
                    mu_omega = rnorm(1, 0, 0),
                    mu_theta = rnorm(1, 0, 0),
                    beta_gamma_year = rnorm(data_list$Nyear, 0, 0),
                    beta_omega_year = rnorm(data_list$Nyear, 0, 0),
                    beta_gamma_season = c(0, rnorm(1, 0, 0)),
                    beta_omega_season = c(0, rnorm(1, 0, 0)),
                    lambda = runif(1, 1, 10)))

### 5. Run Model ---------------------------------------------------------------

Mod <- run.jags(model = "case_study/src/abMod_byblock_jags.R",
              monitor = c("mu_gamma", "mu_omega", "mu_theta",
                          "beta_gamma_year", "beta_gamma_season", 
                          "beta_omega_year", "beta_omega_season", 
                          "lambda", "n"),
              data = data_list,
              n.chains = NCHAINS,
              inits = inits,
              adapt = 500,
              burnin = 2000,
              sample = 2000,
              thin = 1,
              method = "parallel")

Model_matrix1 <- as.matrix(as.mcmc.list(Mod))

Model_matrix2 <- Model_matrix1[, -grep(c("^n\\["), colnames(Model_matrix1))]

# saveRDS(Mod, file = paste0("outputs/Models/abMod_", SITE,"_",modname, ".rds"))

### 6. Plot estimated abundance ------------------------------------------------

### a. get other abundance metrics ---

source("case_study/src/get_cmr_data.R")

CT_abundances <- detHist %>%
  left_join(regions) %>%
  group_by(block, t) %>%
    summarize(date = date[1],
              vole_count = sum(vole_count, ))

all_dates <- detHist %>% 
  group_by(t) %>%
  summarize(date = min(date))

### b. get estimated abundance by the model ---

RN_abundances <- Model_matrix1[, grep("^n\\[", colnames(Model_matrix1))] %>%
  apply(2, function(x){quantile(x, c(0.025,0.5,0.975), na.rm = T)}) %>%
  t %>%
  data.frame() 

colnames(RN_abundances) <- c("inf", "med", "sup")

RN_abundances <- RN_abundances %>%
  mutate(date = as_date(rep(all_dates$date,
                            each = length(unique(regions$block)))),
         block = rep(unique(regions$block),
                     dim(y)[2]))

### c. retrieve site information ---

if(POOLING == "none") {
  RN_abundances <- left_join(RN_abundances, regions) %>%
    mutate(block = site) %>%
    select(-site)
  cmr_data <- left_join(cmr_data, regions) %>%
    mutate(block = site) %>%
    select(-site)
  CT_abundances <- left_join(CT_abundances, regions) %>%
    mutate(block = site) %>%
    select(-site)
}

## d. plot

p1 <- ggplot(RN_abundances)+
  geom_point(data = cmr_data, aes(x = date, y = cmr_estimate/3)) +
  geom_line(data = cmr_data, aes(x = date, y = cmr_estimate/3)) +
  geom_line(aes(x=date, y = med, col = factor(block)), linewidth = 1, show.legend = F) +
  geom_ribbon(aes(x = date, ymin = inf, ymax = sup, fill = factor(block)),
              alpha = 0.50, show.legend = F)+
  scale_y_continuous(name = "RN estimate",
                     sec.axis = sec_axis( transform=~.*3, name="CMR estimate"))+
  facet_wrap(~block)+
  theme_bw()

p2 <- ggplot(CT_abundances)+
  geom_line(aes(x=date, y = vole_count), col = "blue")+
  geom_point(data = cmr_data, aes(x = date, y = cmr_estimate))+
  geom_line(data = cmr_data, aes(x = date, y = cmr_estimate))+
  scale_x_date(breaks = "6 months")+
  scale_y_continuous(name = "CT estimate",
                     sec.axis = sec_axis( transform=~.*5, name="CMR estimate"))+
  facet_wrap(~block)+
  theme_bw()

(p1 | p2)

ggsave("case_study/plots/CMR_vs_RN_vs_CT_byblock.png",  width = 29.7, height = 29.7, unit = "cm")

### 7. Plot parameter values ---------------------------------------------------

mu <- Model_matrix2[, grep("mu_", colnames(Model_matrix2))]

beta_gamma <- Model_matrix2[, grep("beta_gamma", colnames(Model_matrix2))]
beta_omega <- Model_matrix2[, grep("beta_omega", colnames(Model_matrix2))]

gamma <- exp(apply(beta_gamma, 2, function(x){x + mu[,1]}))
omega <- invlogit(apply(beta_omega, 2, function(x){x + mu[,2]}))
theta <- invlogit(mu[,3])

colnames(gamma) <- c(2018:2023, "spring", "winter")
colnames(omega) <- c(2018:2023, "spring", "winter")

par(mfrow = (c(2,2)))
boxplot(gamma[,1:6], main = "gamma - year")
boxplot(gamma[,7:8], main = "gamma - season")
boxplot(omega[,1:6], main = "omega - year")
boxplot(omega[,7:8], main = "omega - season")

# ggsave("case_study/plots/estimates_byblock.png",  width = 29.7, height = 15, unit = "cm")
