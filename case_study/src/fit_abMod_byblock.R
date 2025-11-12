library(tidyverse)
library(lubridate)

library(LaplacesDemon)
library(runjags)
library(rjags)
library(coda)
library(patchwork)

rm(list = ls())

K = 5

### 1. Load Data ---------------------------------------------------------------

ct_data.filename <- "case_study/data/camData_pors.rds"
cmr_data.por.filename <- "case_study/data/Porsanger/cmr/karma_grid.csv"
cmr_data.kar.filename <- "case_study/data/Porsanger/cmr/pors_grid.csv"
regions.filename <- "case_study/data/Porsanger/camera_blocks.csv"

ct_data <- readRDS(ct_data.filename) %>%
  filter(year >= 2017)

ct_data <- ct_data %>%
  mutate(jday = julian(date, origin = min(ct_data$date)) + 1,
         t = (jday - 1) %/% K + 1)

detHist <- ct_data %>%
  group_by(site, t, date, jday) %>% 
  summarize(vole = as.numeric(any(species == "vole")))

regions <- read.csv(regions.filename, sep = ";")

### 2. Get Data Matrix ---------------------------------------------------------

y.df <- detHist%>%
  group_by(site, t) %>%
  mutate(k = 1:n()) %>% 
  pivot_wider(names_from = site, values_from = vole) 

M <- length(unique(detHist$site))
T <- max(detHist$t)

y <- array(NA, dim = c(M, T, K))
for(i in 1:nrow(y.df)){
  y[, y.df$t[i], y.df$k[i]] <- t(y.df[i,unique(detHist$site)]) + 1
}

### 3. Get covariates ----------------------------------------------------------

covs <- ct_data %>% 
  group_by(t) %>%
  summarize(date = date[1],
            int = 1) %>%
  mutate(season = ifelse(yday(date) < 274 & yday(date) > 152 , "spring", "winter")) %>% # spring between June 1st and October 1st
  mutate(year = (data.table::rleid(season) - 1) %/% 2) %>% # change year every begining of spring
  left_join(data.frame(t = 1:T), .)

y <- y[, covs$year != 6,]
covs <- covs[covs$year != 6, ]

### 4. Prepare to run the model ------------------------------------------------

data_list <- list(y = y,
                  region = regions$block,
                  year = as.numeric(as_factor(covs$year)),
                  season = as.numeric(as_factor(covs$season)),
                  T = dim(y)[2], M = M, K = K, 
                  R = length(unique(regions$block)),
                  Nyear = length(unique(covs$year))
                  )

### 5. Run Model ---------------------------------------------------------------

NCHAINS = 2

inits <- map(1:NCHAINS, 
             ~ list(mu_gamma = rnorm(1, 0, 0),
                    mu_omega = rnorm(1, 0, 0),
                    mu_theta = rnorm(1, 0, 0),
                    beta_gamma_year = rnorm(data_list$Nyear, 0, 0),
                    beta_omega_year = rnorm(data_list$Nyear, 0, 0),
                    beta_gamma_season = c(0, rnorm(1, 0, 0)),
                    beta_omega_season = c(0, rnorm(1, 0, 0)),
                    lambda = runif(1, 1, 10)))

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

# saveRDS(Mod, file = paste0("outputs/Models/abMod_", SITE,"_",modname, ".rds"))

### 6. Model Analysis ----------------------------------------------------------

## a. get mcmcs matrices

M.mat_ <- as.matrix(as.mcmc.list(Mod))

M.mat <- M.mat_[, -grep(c("^n\\["), colnames(M.mat_))]

n.qu <- M.mat_[, grep("^n\\[", colnames(M.mat_))] %>%
  apply(2, function(x){quantile(x, c(0.025,0.5,0.975), na.rm = T)})

dates <- detHist %>% 
  filter(t <= dim(y)[2]) %>%
  group_by(t) %>%
  summarize(date = min(date))

## b. summarise estimated abundances

n.df <- data.frame(Nest.inf = n.qu[1,],
                   Nest.med = n.qu[2,],
                   Nest.sup = n.qu[3,],
                   date = as_date(rep(dates$date, each = length(unique(regions$block)))),
                   block = rep(unique(regions$block), dim(y)[2]))

## c. get CMR data 

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

cmr_data <- cmr_data %>%
  left_join(regions) %>%
  group_by(year, seas, date, block) %>% 
  summarize(cmr_estimate = sum(cmr_estimate))

countHist <- ct_data %>%
  mutate(period = jday %/% 5) %>%
  group_by(period, site, jday) %>% 
  summarize(t = t[1], date = date[1],
            volecount = sum(as.numeric(species == "vole"))) %>% 
  left_join(regions) %>%
  group_by(date, block) %>% 
  summarize(volecount = sum(volecount))

## d. plot time series of CT, RN, and CMR abundance estimates

p1 <- ggplot(n.df)+
  geom_point(data = cmr_data, aes(x = date, y = cmr_estimate/3)) +
  geom_line(data = cmr_data, aes(x = date, y = cmr_estimate/3)) +
  geom_line(aes(x=date, y = Nest.med, col = factor(block)), linewidth = 1, show.legend = F) +
  geom_ribbon(aes(x = date, ymin = Nest.inf, ymax = Nest.sup, fill = factor(block)),
              alpha = 0.50, show.legend = F)+
  scale_y_continuous(name = "RN estimate",
                     sec.axis = sec_axis( transform=~.*3, name="CMR estimate"))+
  facet_wrap(~block)+
  theme_bw()

p2 <- ggplot(countHist)+
  geom_line(aes(x=date, y = volecount), col = "grey")+
  geom_point(data = cmr_data, aes(x = date, y = cmr_estimate/5))+
  geom_line(data = cmr_data, aes(x = date, y = cmr_estimate/5))+
  scale_x_date(breaks = "6 months")+
  scale_y_continuous(name = "CT estimate",
                     sec.axis = sec_axis( transform=~.*5, name="CMR estimate"))+
  facet_wrap(~block)+
  ylim(c(0,100))+
  theme_bw()

(p1 | p2)

ggsave("case_study/plots/CMR_vs_RN_vs_CT_byblock.png",  width = 29.7, height = 29.7, unit = "cm")

## e. plot parameter values

gam <- exp(M.mat[, grep("b\\[", colnames(M.mat_))])
colnames(gam) <- colnames(gam_covs)
gam <- pivot_longer(data.frame(gam), cols = colnames(gam), names_to = "level") %>% 
  mutate(param = "gamma")

omega <- invlogit(M.mat[, grep("d\\[", colnames(M.mat_))])
colnames(omega) <- colnames(tau_covs)
omega <- pivot_longer(data.frame(omega), cols = colnames(omega), names_to = "level") %>% 
  mutate(param = "omega")

demo_pars <- rbind(gam, omega) %>%
  mutate(variable = ifelse(grepl("yr", level), "year", "season")) %>%
  filter(level != "yr.6")

ggplot(demo_pars) +
  geom_boxplot(aes(x = level, y = value)) + 
  facet_wrap( ~ paste(param, "-", variable), scales = "free")

# data.frame(M.mat) %>%
#   mutate(rho = invlogit(a),
#          gamma = exp(b),
#          omega = invlogit(d)) %>% 
#   pivot_longer(cols = c("rho", "gamma", "omega"), names_to = "param", values_to = "value") %>%
#   ggplot() + 
#   geom_histogram(aes(x = value)) + 
#   facet_wrap( ~ param, scales = "free") + 
#   theme_bw()
#   
# ggsave("case_study/plots/estimates_byblock.png",  width = 29.7, height = 15, unit = "cm")
