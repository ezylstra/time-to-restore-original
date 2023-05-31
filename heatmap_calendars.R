library(rnpn)
library(stringr)
library(tidyverse)
library(dplyr)
library(httr)
library(readr)
library(lubridate)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(ggbreak)
rm(list= ls())

#Notes
#- One viz to ask their pref on peak flowering vs open flowers, vs both (assuming fl/buds not relevant)

#- One  viz to ask about uncertainty - will be calendar with all the uncertainty combined with three sub vizes that show interannual, spatial and observer separately

#- One viz to see prefs on community heatmap vs species heatmap + open_flowers_estimate curve 


#this script is to develop the Time to Restore prototype calendars of type, heatmap (Alyssa, Sept 2022)
#code significantly updated by Travis Matlock(May 2023)
setwd("~/Documents/My Files/USA-NPN/Data/Analysis/R_default/npn_analyses/TimetoRestore/Data")


species_list <- npn_species()
#Find Time to Restore - priority species - https://docs.google.com/spreadsheets/d/1WyxsCqZAAATIgZvR5ewNcZkaN2mumVi2d1vFw8I0HJo/edit#gid=0

species <- c(186,197,200,201,202,203,204,207,224,781,845,916,931,1163,1167,1334)
phenophases = c(501, 390) # open flowers and ripe fruits

df <- npn_download_status_data(request_source="Alyssa",
                               years=c(2009:2023),
                               species_ids = c(paste0(species)), # species codes
                               state = c("TX", "LA", "OK", "NM"),
                               additional_fields = c("observedby_person_id"),
                               phenophase_ids= c(paste0(phenophases)), 
                               climate_data = FALSE)


write.csv(df, file="status_16spp_2009-2023.csv")
df <- (read.csv("status_16spp_2009-2023.csv"))

#format dates as needed and remove uncertain (?) records
df <- df %>%
  mutate(year = lubridate::year(observation_date)) %>%
  mutate(month = lubridate::month(observation_date)) %>%
  mutate(week = lubridate::week(observation_date))  %>%
  filter(phenophase_status != -1) # -1 is when people reported ?
  #filter(!species_id %in% c(1167,1163,197,931,202))

str(df)
colnames(df)

df$common_name <- as.factor(df$common_name)
df$phenophase_status <- as.numeric(df$phenophase_status)

#Ignore for now - filters used to get buttonbush observer variability plot
  filter(individual_id == 141920) %>%
  filter(observedby_person_id < 43389, observedby_person_id > 38000, observedby_person_id != 39173) %>%
  filter(year !=2017, year !=2022)

#calculate the proportion of records in to the total records by week for heatmaps
df1 <- df %>%
  group_by(week, common_name, phenophase_id) %>%
  mutate(num_records = n()) %>%
  mutate(sum_pp = sum(phenophase_status))%>%
  mutate(proportion = sum_pp/num_records) %>%
  ungroup()

#Wrangle peak dataset
#to get all species/states etc load the file 8priority_spp_2013-2022_w_open_flower_estimate.csv
df_peak <- (read.csv("buttonbush2013-2022_w_open_flower_estimate.csv"))
df_peak <- df_peak %>%
  mutate(week = lubridate::week(observation_date))  %>%
  filter(state == "LA")

#Ignore for now, this code creates the option "Heatmaps for Both Measures" on the survey, that no one chose
df_peak <- df_peak %>%
  group_by(week) %>%
  mutate(num_records = n()) %>%
  mutate(sum_peak = sum(peak))%>%
  mutate(proportion = sum_peak/num_records) %>%
  filter(year !=2017, year !=2022) %>%
  filter(!is.na(proportion)) %>%
  ungroup()

#We should check this df_peak for NAs in the open flowers estimate column and decide what to do with them

#COMMUNITY HEATMAPS 
#Creates dataframe for which species should use which colors - not used currently
species_color <- data.frame(species_id = species, color=c("lightgoldenrod1", "palevioletred1", "palegreen1", "burlywood1",
                                                          "plum1", "gold", "orchid2", "papayawhip", "rosybrown1", "indianred1",
                                                          "purple1", "firebrick1", "thistle2", "darkolivegreen1", "yellow", "lightcyan2"))

# Create data frames based on states' priority species
df_LA <- df1 %>%
  # Filter based on states relevant to selected species for this state
  filter(state == "LA") %>%
  # Filter for species
  filter(species_id == 201 | species_id == 1334) %>%
  # Assign different colors to different species
  mutate(color = case_when(species_id == 201 ~ "burlywood1",
                           species_id == 1334 ~ "lightcyan2"))
# Repeat for all states
df_OK <- df1 %>%
  filter(state != "NM" ) %>%
  filter(species_id == 781 | species_id == 201 | species_id == 203) %>%
  mutate(color = case_when(species_id == 201 ~ "burlywood1",
                           species_id == 781 ~ "indianred1",
                           species_id == 203 ~ "gold"))
df_NM <- df1 %>%
  filter(state == "NM" | state == "TX") %>%
  filter(species_id == 224 | species_id == 1163 | species_id == 1167 | species_id == 203) %>%
  mutate(color = case_when(species_id == 224 ~ "rosybrown1",
                           species_id == 1163 ~ "darkolivegreen1",
                           species_id == 1167 ~ "yellow",
                           species_id == 203 ~ "gold"))
df_TX <- df1 %>%
  filter(state == "TX") %>%
  filter(species_id == 201 | species_id == 202 | species_id == 203) %>%
  mutate(color = case_when(species_id == 201 ~ "burlywood1",
                           species_id == 202 ~ "plum1",
                           species_id == 203 ~ "gold"))



#Iterate to make plots for each species-state combo for the community heatmap calendars
plots <- list()
# Loop over each species
for (i in unique(df_TX$common_name)) {
  # Create temporary dataset filtered by individual species
  tmp <- df_TX %>%
    filter(common_name == i) 
  # Add the graphs for each individual species to the plot list
  plots[[i]] <- tmp %>%
    ggplot() +
    geom_tile(aes(x=week, y=phenophase_description, fill=proportion)) +
    theme(axis.text.x = element_text(hjust=-0.6), axis.text.y = element_text(), axis.ticks.y = element_blank(), legend.position="none", panel.background = element_blank()) +
    labs(x = "", y ="") +
    scale_fill_gradient(low="grey", high=tmp$color[1], (name="Proportion")) +
    facet_wrap(tmp$common_name, strip.position="top") + 
    scale_x_continuous(breaks=seq(1, 52, by=4.25), limits=c(1,52),
                       labels=c("Jan","Feb","Mar", "Apr", "May", "June", "July", "Aug", "Sep", "Oct", "Nov", "Dec", "")) #+
    #coord_fixed(ratio=2) # Can be used to adjust aspect ratio if needed
}

# Stack all plots
purrr::reduce(plots, `/`)
# For a 2x2 grid
(plots[[1]] + plots[[2]]) / (plots[[3]] + plots[[4]])


#SPECIES HEATMAP + OPEN_FLOWERS_ESTIMATE CURVE
#create plot
 p2 = ggplot() +
  geom_line(data = df_peak, aes(x=week, y=open_flower_estimate),color="orange", size = 2) +
  labs(y = "# Open Flowers ", x = " ") +
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank())
plot(p2) + theme_minimal()

#use the copy plot, and resizing to get it to come out long and narrow
#maybe in the future explore the reshape library to make it circular like the year?

