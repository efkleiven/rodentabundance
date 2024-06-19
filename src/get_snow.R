library(tidyverse)
library(foreach)
library(lubridate)

rm(list = ls())

ct_data.filename <- "data/camData.rds"

ct_data <- readRDS(ct_data.filename)

snow_data <- ct_data %>%
  filter(trigger == "T")%>%
  group_by(station) %>%
  mutate(cond1.L = zoo::rollapply(temp, 6,  sd,
                                  partial = T, align = "left") < 0.7,
         cond1.R = zoo::rollapply(temp, 6,  sd,
                                  partial = T, align = "right") < 0.7,
         cond2.L = zoo::rollapply(temp, 6,  FUN = function(x){all(abs(x) <= 2)},
                                  partial = T, align = "left"),
         cond2.R = zoo::rollapply(temp, 6,  FUN = function(x){all(abs(x) <= 2)},
                                  partial = T, align = "right"),
         SNOW = (cond1.L & cond2.L)|(cond1.R & cond2.R)) %>%
  mutate(SNOW = ifelse(is.na(SNOW), FALSE, SNOW))

snow_data <- snow_data %>%
  group_by(station, date) %>%
  summarize(SNOW = any(SNOW))
saveRDS(snow_data, "data/snowData.rds")

