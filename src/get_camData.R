library(tidyverse)
library(foreach)
library(lubridate)

rm(list = ls())

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#################################### Hakoya ####################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
years.hakoy <- 2018:2020

answ_to_sp_hakoy <- data.frame(answer = 0:7,
                               species = c("bad.quality","empty","bird","microtus","least_weasel",
                                     "lemming","shrew","stoat"))

hakoy_data.path <- "data/Hakoya"
hakoy_data.filenames <- list.files(hakoy_data.path, full.names = TRUE, pattern = "classifications")
hakoy_metadata.filenames <- list.files(hakoy_data.path, full.names = TRUE, pattern = "metadata")

hakoy_data <- foreach(file = 1:length(hakoy_data.filenames), .combine = rbind)%do%{
  read.csv(hakoy_data.filenames[file], stringsAsFactors = FALSE, header=TRUE) %>%
    mutate(site = "Hakoya") %>%
    mutate(filename = str_split_i(fileName, pattern = "\\\\", i = 3),
           filename = gsub("'", "", filename),
           station = str_split_i(filename, pattern = "_", i = 1),
           year = years.hakoy[file])
}

hakoy_metadata <- foreach(file = hakoy_metadata.filenames, .combine = rbind)%do%{
  read.table(file, stringsAsFactors = FALSE, header=TRUE)%>%
    select(NewFileName, DateTimeOriginal, AmbientTemperature, TriggerMode)%>%
    rename(filename = NewFileName,
           datetime = DateTimeOriginal,
           temp = AmbientTemperature,
           trigger = TriggerMode)
}

hakoy_data <- left_join(hakoy_data, hakoy_metadata, by = "filename")%>%
  arrange(station, datetime)

hakoy_data <- hakoy_data %>%
  mutate(answer = as.numeric(confidence1 > 0.95) * guess1)%>%
  left_join(answ_to_sp_hakoy, by = "answer")

hakoy_data <- hakoy_data %>%
  mutate(datetime = as_datetime(datetime),
         date = date(datetime), 
         time = format(datetime, format="%H:%M:%S")) 

pics_to_keep <- hakoy_data %>%
  group_by(station)%>%
  mutate(dt = as.numeric(difftime(datetime, lag(datetime, 1), units = "secs")),
         event = cumsum(dt>10|is.na(dt)))%>%
  group_by(station, event)%>% 
  summarise(filename = filename[which.max(confidence1)])%>%
  .$filename

hakoy_data <- hakoy_data %>%
  filter(filename %in% pics_to_keep)

hakoy_data <- hakoy_data %>%
  select(site, station, year, trigger, datetime, date, time, temp, species)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
################################### Valdres ####################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

answ_to_sp_valdres <- data.frame(Species = c("M", "S", "DS", "VS", "K",
                                            "L", "SK", "R", "W"),
                                 species = c("microtus", "shrew", "shrew", "shrew", "bank_vole",
                                             "lemming", "lemming", "stoat", "least_weasel"))

Valdres_data.path <- "data/Valdres"
Valdres_data.filenames <- list.files(Valdres_data.path, full.names = TRUE, pattern = ".csv")

valdres_data <- foreach(file = Valdres_data.filenames, .combine = rbind) %do% {
  read.csv(file, sep = ";", header=TRUE)%>%
    mutate(Temp_C = ifelse(Temp_C > Temp_F, Temp_F, Temp_C)) %>%
    rename(site = Site,
           year = Year,
           date = Date,
           time = Time,
           temp = Temp_C,
           trigger = Trigger)
}

valdres_data <- valdres_data %>%
  mutate(station = paste0("V", substr(Location_season, 4,4)))

valdres_data <- left_join(valdres_data, answ_to_sp_valdres, by = "Species") %>%
  filter(!is.na(species) | trigger == "T") %>%
  mutate(species = ifelse(is.na(species), "empty", species))

valdres_data <- valdres_data %>% 
  mutate(date = gsub("/",".", date),
         datetime = as.POSIXct(paste(date, time, sep = "-"), format="%d.%m.%Y-%H:%M:%S", tz = "UTC"),
         date = date(datetime)) %>%
  arrange(station, datetime) %>%
  select(site, station, year, trigger, datetime, date, time, temp, species)

ct_data_allsites <- rbind(valdres_data, hakoy_data)

# saveRDS(ct_data_allsites, "data/camData.rds")
