cmr_data.paths <- "case_study/data/Porsanger/cmr/"

site_labs <- gsub("-", "\\.", unique(ct_data$site))
                     
cmr_data <- map(list.files(cmr_data.paths, full.names = TRUE, pattern = "grid"),
                read.csv) %>%
  reduce(right_join) %>%
  select(sp, seas, year, all_of(site_labs))

cmr_data <- cmr_data %>%
  filter(year > 2017) %>%
  mutate(date = ifelse(seas == "FALL", "-09-15", "-06-15")) %>%
  mutate(date = as_date(paste0(year, date)))

cmr_data <- cmr_data %>%
  group_by(year, seas, date) %>%   
  summarise(
    across(all_of(site_labs), \(x) sum(x, na.rm = TRUE)),
    .groups = "drop"
  )

cmr_data <- cmr_data %>%
  pivot_longer(cols = all_of(colnames(cmr_data)[4:ncol(cmr_data)]),
               values_to = "cmr_estimate", names_to = "site") %>%
  mutate(site = gsub("\\.", "-", site)) %>%
  filter(site %in% unique(ct_data$site))

cmr_data <- cmr_data %>%
  left_join(regions) %>%
  group_by(year, seas, date, block) %>% 
  summarize(cmr_estimate = sum(cmr_estimate))
