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

n.df <- data.frame(Nest.inf = n.qu[1,],
                   Nest.med = n.qu[2,],
                   Nest.sup = n.qu[3,],
                   t = rep(1:T, each = M),
                   site = rep(unique(detHist$site), T))

ggplot(n.df)+
  geom_line(aes(x=t, y = Nest.med, col = site), linewidth = 1)+
  geom_ribbon(aes(x = t, ymin = Nest.inf, ymax = Nest.sup, fill = site),
              alpha = 0.50)+
  facet_wrap(~site)+
  theme_bw()

cmr_data <- read.csv(cmr_data.filename) %>%
  filter(year > 2017) %>%
  mutate(date = ifelse(seas == "FALL", "-09-15", "-06-15")) %>%
  mutate(date = as_date(paste0(year, date)))

cmr_data.por.filename <- "case_study/data/Porsanger/cmr/karma_grid.csv"
cmr_data.kar.filename <- "case_study/data/Porsanger/cmr/pors_grid.csv"

cmr_data.por <- read.csv(cmr_data.por.filename) %>%
  filter(year > 2017) %>%
  mutate(date = ifelse(seas == "FALL", "-09-15", "-06-15")) %>%
  mutate(date = as_date(paste0(year, date)))%>%
  group_by(year, seas) %>%   
  summarise(
    across(G1:G20, \(x) sum(x, na.rm = TRUE)),
    .groups = "drop"
  )

cmr_data.kar <- read.csv(cmr_data.kar.filename) %>%
  filter(year > 2017) %>%
  mutate(date = ifelse(seas == "FALL", "-09-15", "-06-15")) %>%
  mutate(date = as_date(paste0(year, date))) %>%
  group_by(year, seas) %>%   
  summarise(
    across(T12:T52, \(x) sum(x, na.rm = TRUE)),
    .groups = "drop"
  )


p1 <- ggplot(n.df)+
  geom_point(data = cmr_data, aes(x = date, y = PorsSum/7), col = "red")+
  geom_line(data = cmr_data, aes(x = date, y = PorsSum/7), col = "red")+
  geom_line(aes(x=date, y = Nest.med), col = "blue")+
  geom_ribbon(aes(x = date, ymin = Nest.inf, ymax = Nest.sup),
              fill = "lightblue", alpha = 0.5)+
  scale_x_date(breaks = "6 months")+
  scale_y_continuous(name = "RN estimate",
                     sec.axis = sec_axis( transform=~.*7, name="CMR estimate"))+
  theme_bw()

p2 <- ggplot(countHist)+
  geom_point(data = cmr_data, aes(x = date, y = PorsSum/2), col = "red")+
  geom_line(data = cmr_data, aes(x = date, y = PorsSum/2), col = "red")+
  geom_line(aes(x=date, y = volecount))+
  scale_x_date(breaks = "6 months")+
  scale_y_continuous(name = "CT estimate",
                     sec.axis = sec_axis( transform=~.*2, name="CMR estimate"))+
  theme_bw()

p1/p2

ggsave("case_study/plots/CMR_vs_RN_vs_CT.png",  width = 29.7, height = 29.7, unit = "cm")

# date.df <- data.frame(date = seq(as.Date("2019/1/1"), as.Date("2019/12/31"), by = "day"),
#                       int = 1,
#                       yday = 1:365) %>%
#   mutate(yday2 = yday ** 2)
# 
# date.df$omega <- M.mat[, grep("d\\[", colnames(M.mat_))]%>%
#   apply(2, median) %*% t(as.matrix(date.df[,c("int", "yday", "yday2")])) %>%
#   c %>%
#   invlogit
# 
# date.df %>% 
# ggplot()+
#   geom_line(aes(x = date, y = omega)) +
#   ggtitle("survival rate") +
#   theme_bw()
# 
# a_int <- M.mat[, grep("a_RE_stat\\[", colnames(M.mat_))]%>%
#   apply(1, mean)
# a_SNOW <- M.mat[,"a"]
# 
# a.df <- data.frame(NOSNOW = invlogit(a_int),
#                    SNOW = invlogit(a_int+a_SNOW))
# 
# ggplot(a.df)+
#   geom_violin(aes(x = "NOSNOW", y = NOSNOW), fill = "lightgreen")+
#   geom_violin(aes(x = "SNOW", y = SNOW), fill = "lightblue")+
#   geom_point(aes(x = "NOSNOW", y = NOSNOW), col = "red", size = .5)+
#   geom_point(aes(x = "SNOW", y = SNOW), col = "red", size = .5)+
#   geom_segment(aes(x = "NOSNOW", xend = "SNOW", y = NOSNOW, yend = SNOW), size = .1) +
#   ggtitle("individual detection probability") +
#   theme_bw()
# 
# d_int <- M.mat[, "d[1]"]
# d_SNOW <- M.mat[,"d[2]"]
# 
# d.df <- data.frame(NOSNOW = invlogit(d_int),
#                    SNOW = invlogit(d_int+d_SNOW))
# 
# ggplot(d.df)+
#   geom_violin(aes(x = "NOSNOW", y = NOSNOW), fill = "lightgreen")+
#   geom_violin(aes(x = "SNOW", y = SNOW), fill = "lightblue")+
#   geom_point(aes(x = "NOSNOW", y = NOSNOW), col = "red", size = .5)+
#   geom_point(aes(x = "SNOW", y = SNOW), col = "red", size = .5)+
#   geom_segment(aes(x = "NOSNOW", xend = "SNOW", y = NOSNOW, yend = SNOW), size = .1) +
#   ggtitle("survival rate") +
#   theme_bw()