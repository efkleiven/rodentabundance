library(tidyverse)
library(foreach)
library(lubridate)

library(LaplacesDemon)
library(runjags)
library(rjags)
library(coda)
library(patchwork)

rm(list = ls())

modname = ""
K = 5

### 1. Load Data ---------------------------------------------------------------

ct_data.filename <- "case_study/data/camData_pors.rds"
cmr_data.filename <- "case_study/data/Porsanger/cmr/porsanger_mus_reg_2025.csv"

ct_data <- readRDS(ct_data.filename) %>%
  filter(year >= 2017)

day.orig <- min(ct_data$date)

ct_data <- ct_data %>%
  mutate(jday = julian(date, origin = day.orig) + 1,
         t = (jday - 1) %/% K + 1)

detHist <- ct_data %>%
  group_by(site, t, date, jday) %>% 
  summarize(vole = as.numeric(any(species == "vole")))

countHist <- ct_data %>%
  mutate(period = jday %/% 5) %>%
  group_by(period, jday) %>% 
  summarize(t = t[1], date = date[1],
            volecount = sum(as.numeric(species == "vole")))

### 2. Get Data Matrix ---------------------------------------------------------

y.df <- detHist%>%
  group_by(site, t) %>%
  mutate(k = 1:n()) %>% 
  pivot_wider(names_from = site, values_from = vole) 

M <- length(unique(detHist$site))
T <- max(detHist$t)

y.mat <- array(NA, dim = c(M, T, K))
for(i in 1:nrow(y.df)){
  y.mat[, y.df$t[i], y.df$k[i]] <- t(y.df[i,unique(detHist$site)]) + 1
}

### 3. Get covariates ----------------------------------------------------------

covs <- ct_data %>% 
  group_by(t) %>%
  summarize(year = year[1],
            date = date[1],
            int = 1) %>%
  mutate(breeding = as.numeric(month(date) %in% 5:9)) %>%
  left_join(data.frame(t = 1:T), .)

years <- covs$year
Nyears <- length(unique(years))

### 4. Prepare to run the model ------------------------------------------------

# individual detection
rho_covs <- as.matrix(covs[, c("int")])

# natality
gam_covs <- as.matrix(covs[, c("int")])

#survival rate
tau_covs <- as.matrix(covs[, c("int")])

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

Mod <- run.jags(model = "case_study/src/abMod_bycam_jags.R",
              monitor = c("a", "b", "d", "lam", "n"),
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

M.mat <- M.mat_[, -grep(c("n"), colnames(M.mat_))]

n.qu <- M.mat_[, grep("n", colnames(M.mat_))] %>%
  apply(2, function(x){quantile(x, c(0.025,0.5,0.975), na.rm = T)})

dates <- detHist %>% 
  group_by(t) %>%
  summarize(date = min(date))

n.df <- data.frame(Nest.inf = n.qu[1,],
                   Nest.med = n.qu[2,],
                   Nest.sup = n.qu[3,],
                   date = as_date(rep(dates$date, each = M)),
                   site = rep(unique(detHist$site), T))

cmr_data.por.filename <- "case_study/data/Porsanger/cmr/karma_grid.csv"
cmr_data.kar.filename <- "case_study/data/Porsanger/cmr/pors_grid.csv"

cmr_data.por <- read.csv(cmr_data.por.filename) %>%
  filter(year > 2017) %>%
  mutate(date = ifelse(seas == "FALL", "-09-15", "-06-15")) %>%
  mutate(date = as_date(paste0(year, date)))

cmr_data.por <- cmr_data.por %>%
  group_by(year, seas, date) %>%   
  summarise(
    across(all_of(colnames(cmr_data.por)[1:11]), \(x) sum(x, na.rm = TRUE)),
    .groups = "drop"
  )

cmr_data.kar <- read.csv(cmr_data.kar.filename) %>%
  filter(year > 2017) %>%
  mutate(date = ifelse(seas == "FALL", "-09-15", "-06-15")) %>%
  mutate(date = as_date(paste0(year, date))) %>%
  select(-X)

cmr_data.kar <- cmr_data.kar %>%
  group_by(year, seas, date) %>%   
  summarise(
    across(all_of(colnames(cmr_data.kar)[1:8]), \(x) sum(x, na.rm = TRUE)),
    .groups = "drop"
  )

cmr_data <- left_join(cmr_data.por, cmr_data.kar)
cmr_data <- cmr_data %>%
  pivot_longer(cols = all_of(colnames(cmr_data)[4:ncol(cmr_data)]),
               values_to = "cmr_estimate", names_to = "site") %>%
  mutate(site = gsub("\\.", "-", site)) %>%
  filter(site %in% unique(ct_data$site))

ggplot(n.df)+
  geom_point(data = cmr_data, aes(x = date, y = cmr_estimate/2)) +
  geom_line(data = cmr_data, aes(x = date, y = cmr_estimate/2)) +
  geom_line(aes(x=date, y = Nest.med, col = site), linewidth = 1) +
  geom_ribbon(aes(x = date, ymin = Nest.inf, ymax = Nest.sup, fill = site),
              alpha = 0.50)+
  scale_y_continuous(name = "RN estimate",
                     sec.axis = sec_axis( transform=~.*2, name="CMR estimate"))+
  facet_wrap(~site)+
  theme_bw()

ggsave("case_study/plots/CMR_vs_RN_bycam.png",  width = 29.7, height = 29.7, unit = "cm")
