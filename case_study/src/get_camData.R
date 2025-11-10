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

saveRDS(ct_data_allsites, "data/camData.rds")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
################################### Porsanger ###################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

getwd()

porsanger_years <- dir("case_study/data/Porsanger/automatic_classification")
porsanger_path <- "case_study/data/Porsanger/automatic_classification"
porsanger_meta_path <- "case_study/data/Porsanger/image_metadata"

cat_to_species <- data.frame(answer = 0:7,
                               species = c("bad.quality","empty","bird","vole","least_weasel",
                                           "lemming","shrew","stoat"))

#load camera data
# Get all .txt file paths from all subfolders
cam_files <- list.files(
  path = file.path(porsanger_path, porsanger_years),
  pattern = "\\.txt$",        # only .txt files
  full.names = TRUE,          # return full paths
  recursive = TRUE            # include files in nested folders
)

# Check what files were found
print(cam_files)

# Read all .txt files into a list
library(readr)
cam_list <- lapply(cam_files, read_delim)

# (optional) Combine into one big data frame (if structure is identical)
cam_data <- dplyr::bind_rows(cam_list, .id = "source_file")

#load metadata
# Get all .txt file paths from all subfolders
meta_files <- list.files(
  path = file.path(porsanger_meta_path, porsanger_years),
  pattern = "\\.txt$",        # only .txt files
  full.names = TRUE,          # return full paths
  recursive = TRUE            # include files in nested folders
)

# Check what files were found
print(meta_files)

# Read all .txt files into a list
meta_list <- lapply(meta_files, read_delim, delim = "\t")

# (optional) Combine into one big data frame (if structure is identical)
meta_data <- dplyr::bind_rows(meta_list, .id = "source_file") %>%
  mutate(filename = NewFileName,
         datetime = DateTimeOriginal,
         temp = AmbientTemperature,
         trigger = TriggerMode) %>%
  select(filename, datetime, temp, trigger)


pors_data <- left_join(cam_data, meta_data, by = "filename")%>%
  arrange(site, datetime)

pors_data <- pors_data %>%
  mutate(answer = as.numeric(confidence1 > 0.95) * guess1)%>%
  left_join(cat_to_species, by = "answer")

pors_data <- pors_data %>%
  mutate(datetime = as_datetime(datetime),
         date = date(datetime), 
         time = format(datetime, format="%H:%M:%S")) 

pics_to_keep <- pors_data %>%
  group_by(site)%>%
  mutate(dt = as.numeric(difftime(datetime, lag(datetime, 1), units = "secs")),
         event = cumsum(dt>10|is.na(dt)))%>%
  group_by(site, event)%>% 
  summarise(filename = filename[which.max(confidence1)])%>%
  .$filename

pors_data <- pors_data %>%
  filter(filename %in% pics_to_keep)

pors_data <- pors_data %>%
  select(site, year, trigger, datetime, date, time, temp, species)

#
saveRDS(pors_data, "case_study/data/camData_pors.rds")



