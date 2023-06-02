library(daymetr)
library(lubridate)
library(dplyr)
rm(list=ls())

DataFolder <- "~/Documents/My Files/USA-NPN/Data/Analysis/R_default/npn_analyses/TimetoRestore/Data"

#https://tmieno2.github.io/R-as-GIS-for-Economists/daymet-with-daymetr-and-feddata.html

get_daymet <- function(i){
  
  temp_lat <- stations[i, ] %>% pull(latitude)
  temp_lon <- stations[i, ] %>% pull(longitude)
  temp_site <-stations[i, ] %>% pull(site_id)
  
  temp_daymet <- download_daymet(
    lat = temp_lat,
    lon = temp_lon,
    start = 2021,
    end = 2022
  ) %>% 
    #--- just get the data part ---#
    .$data %>% 
    #--- convert to tibble (not strictly necessary) ---#
    as_tibble() %>% 
    #--- assign site_id so you know which record is for which site_id ---#
    mutate(site_id = temp_site) %>% 
    #--- get date from day of the year ---#
    mutate(date = as.Date(paste(year, yday, sep = "-"), "%Y-%j"))
  
  return(temp_daymet)
}  

#check out the results of the first row
get_daymet(1)

#run for all stations
(
  daymet_df <- lapply(1:nrow(stations), get_daymet) %>% 
    #--- need to combine the list of data.frames into a single data.frame ---#
    bind_rows()
)

#https://carpentries-incubator.github.io/data-harvesting-for-agriculture/09-WeatherData/index.html
daymet_df$month <- lubridate::month(daymet_df$date, label = TRUE)
head(daymet_df$month)

str(daymet_df)

write.csv(daymet_df, file="daymet_stations_16spp.csv")
                              
fall_ref_df <- daymet_df %>%
  group_by(site_id)  %>%
  filter(year == 2021, month %in% c("Sep","Oct","Nov")) %>%
  mutate(prcp_fall = (sum(prcp..mm.day.))) %>%
  mutate(tmin_fall = (mean(tmin..deg.c.)))%>%
  mutate(tmax_fall = (mean(tmax..deg.c.)))%>%
  distinct(site_id, .keep_all = TRUE)


winter_ref_df <- daymet_df %>%
  group_by(site_id)  %>%
  filter(year == 2021 & month %in% c("Dec") | year == 2022 & month %in% c("Jan","Feb")) %>%
  mutate(prcp_winter = (sum(prcp..mm.day.))) %>%
  mutate(tmin_winter = (mean(tmin..deg.c.)))%>%
  mutate(tmax_winter = (mean(tmax..deg.c.)))%>%
  distinct(site_id, .keep_all = TRUE)

spring_ref_df <- daymet_df %>%
  group_by(site_id)  %>%
  filter(year == 2022, month %in% c("Mar","Apr","May")) %>%
  mutate(prcp_spring = (sum(prcp..mm.day.))) %>%
  mutate(tmin_spring = (mean(tmin..deg.c.)))%>%
  mutate(tmax_spring = (mean(tmax..deg.c.)))%>%
  distinct(site_id, .keep_all = TRUE)


summer_ref_df <- daymet_df %>%
  group_by(site_id)  %>%
  filter(year == 2022, month %in% c("Jun","Jul","Aug")) %>%
  mutate(prcp_summer = (sum(prcp..mm.day.))) %>%
  mutate(tmin_summer = (mean(tmin..deg.c.)))%>%
  mutate(tmax_summer = (mean(tmax..deg.c.)))%>%
  distinct(site_id, .keep_all = TRUE)


spring_fall <- merge(fall_ref_df, spring_ref_df, by = "site_id")

sum_win <- merge(summer_ref_df, winter_ref_df, by = "site_id")

all_ref <- merge(sum_win, spring_fall, by = "site_id")
all_ref1 <- all_ref %>%
  select("site_id","prcp_summer","prcp_fall","prcp_winter", "prcp_spring",
         "tmin_summer","tmin_fall","tmin_winter","tmin_spring",
         "tmax_summer","tmax_fall","tmax_winter","tmax_spring")


hist(all_ref1$prcp_summer)
hist(all_ref1$tmax_summer)
hist(all_ref1$tmin_summer)

write.csv(all_ref1, file="daymet_ref_file_2022_16spp.csv")

setwd(DataFolder)
daymet2022 <- (read.csv("daymet_ref_file_2022.csv"))

#code to use if you just want one point
download_daymet(
  site = "Daymet",
  lat = 36.0133,
  lon = -84.2625,
  start = 2021,
  end = 2022,
  internal = FALSE,
  path = DataFolder,
  silent = FALSE,
  force = FALSE,
  simplify = TRUE
)
