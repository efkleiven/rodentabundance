library(tidyverse)

ct_data.filename <- "case_study/data/camData_pors.rds"

ct_data <- readRDS(ct_data.filename) %>%
  # filter(site == SITE) %>%
  # left_join(readRDS(snow_data.filename), by = c("station", "date")) %>%
  filter(year >= 2017)

day.orig <- min(ct_data$date)
K <- 5

ct_data <- ct_data %>%
  mutate(jday = julian(date, origin = day.orig) + 1,
         t = (jday - 1) %/% K + 1)

table(ct_data$species)

cam_function <- ct_data %>%
  group_by(site, t, date, jday) %>% 
  #summarize(bad_q = as.numeric(any(species == "bad_quality")))
  summarize(bad_q = any(species == "bad.quality"))

table(cam_function$bad_q)
summary(cam_function$bad_q)

cam_weekly <- cam_function %>%
  group_by(site, t) %>%                          # group by site, t, and week
  summarise(bad_q_sum = sum(bad_q, na.rm = TRUE),      # sum bad_q within each group
            n_days = n(),
            first_date = min(date), 
            .groups = "drop")

cam_weekly
save(cam_weekly, file="missing_data.rda")

# get one date per t
t_labels <- cam_weekly %>%
  group_by(t) %>%
  summarise(first_date = min(first_date)) %>%
  mutate(year_month = format(first_date, "%Y-%m")) %>%
  group_by(year_month) %>%
  slice(1) %>%   # first t of each month
  ungroup() %>%
  slice(seq(1, n(), by = 2))  # keep every 2nd month

#plot heatmap
p1 <- ggplot(cam_weekly, aes(x = t, y = site, fill = bad_q_sum)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma") +
  scale_x_continuous(
    breaks = t_labels$t,
    labels = t_labels$year_month
  ) +
  #facet_wrap(~ t, ncol = 1, scales = "free_y") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "", y = "", fill = "bad_quality")

p1

ggsave(p1, file="case_study/plots/NA_heatmap_porsanger.png", width=24, height = 12, units="cm")
