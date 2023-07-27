library(rnpn)
library(stringr)
library(tidyverse)
library(dplyr)
library(httr)
library(readr)
library(lubridate)
library(ggplot2)
library(patchwork)
library(ggbreak)
library(ggtext)
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
phenophases = c(390, 501) # open flowers and ripe fruits

df <- npn_download_status_data(request_source="Alyssa",
                               years=c(2009:2023),
                               species_ids = c(paste0(species)), # species codes
                               state = c("TX", "LA", "OK", "NM"),
                               additional_fields = c("observedby_person_id", "Site_Name", "Network_Name"),
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
  mutate(sum_pp = sum(phenophase_status)) %>%
  mutate(proportion = sum_pp/num_records) %>%
  ungroup() %>%
  # Format data table for plots requested by Erin
  mutate(phenophase_description = ifelse(phenophase_description == "Ripe fruits", "Ripe seeds", "Open flowers")) %>%
  rowwise() %>%
  mutate(both_names = paste(common_name, " (", "*",genus, " ", species, "*", ")", sep="")) %>%
  ungroup()

df1$both_names <- as.factor(df1$both_names)

#Wrangle peak dataset
#to get all species/states etc load the file 8priority_spp_2013-2022_w_open_flower_estimate.csv
df_peak <- read.csv("16priority_spp_2013-2022_w_open_flower_estimate.csv")
df_peak <- df_peak %>%
  mutate(week = lubridate::week(observation_date)) # %>%
  #filter(state == "LA")

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
species_color <- data.frame(species_id = species, color=c("lightgoldenrod1", "palevioletred1", "palegreen1", "sienna1",
                                                          "darkorchid1", "gold", "orchid2", "papayawhip", "hotpink", "indianred1",
                                                          "purple", "firebrick1", "thistle2", "seagreen3", "yellow", "chartreuse"))

# Create data frames based on states' priority species
df_LA <- df1 %>%
  # Filter based on states relevant to selected species for this state
  filter(state == "LA") %>%
  # Filter for species
  filter(species_id == 201 | species_id == 1334) %>%
  # Assign different colors to different species
  mutate(color = case_when(species_id == 201 ~ "sienna1",
                           species_id == 1334 ~ "chartreuse1")) %>%
  arrange(common_name)

# Repeat for all states
df_OK <- df1 %>%
  filter(state != "NM" ) %>%
  filter(species_id == 781 | species_id == 201 | species_id == 203) %>%
  mutate(color = case_when(species_id == 201 ~ "sienna1",
                           species_id == 781 ~ "indianred1",
                           species_id == 203 ~ "gold")) %>%
  arrange(common_name)

df_NM <- df1 %>%
  filter(state == "NM" | state == "TX") %>%
  filter(species_id == 224 | species_id == 1163 | species_id == 1167 | species_id == 203) %>%
  mutate(color = case_when(species_id == 224 ~ "hotpink",
                           species_id == 1163 ~ "seagreen3",
                           species_id == 1167 ~ "yellow",
                           species_id == 203 ~ "gold")) %>%
  arrange(common_name)

df_TX <- df1 %>%
  filter(state == "TX") %>%
  filter(species_id == 201 | species_id == 202 | species_id == 203) %>%
  mutate(color = case_when(species_id == 201 ~ "sienna1",
                           species_id == 202 ~ "darkorchid1",
                           species_id == 203 ~ "gold")) %>%
  arrange(common_name)


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
    facet_wrap(tmp$both_names, strip.position="top") + 
    theme(axis.text.x = element_text(hjust=-0.6), axis.text.y = element_text(), strip.text = ggtext::element_markdown(),
          axis.ticks.y = element_blank(), legend.position="left", panel.background = element_rect(fill="white"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.key.size = unit(0.5, "cm"), strip.background = element_rect("white")) +
    labs(x = "", y ="") +
    scale_fill_gradient(low="gray90", high=tmp$color[1], (name="Proportion")) +
    scale_x_continuous(breaks=seq(1, 52, by=4.25), limits=c(1,52),
                       labels=c("Jan","Feb","Mar", "Apr", "May", "June", "July", "Aug", "Sep", "Oct", "Nov", "Dec", "")) +
    scale_y_discrete(limits=rev) 
    #coord_fixed(ratio=2) # Can be used to adjust aspect ratio if needed
}

# Stack all plots
purrr::reduce(plots, `/`)

(plots[[1]]+plots[[2]])/(plots[[3]]+plots[[4]])


 #SPECIES HEATMAP + OPEN_FLOWERS_ESTIMATE CURVE

# NM showy milkweed, LA buttonbush, OK silver maple, TX sunflower

df_peak_TX <- df_peak %>%
  filter(state == "TX") %>%
  filter(species_id == 203) %>%
  filter(open_flower_estimate != 0.00)

inat_TX <- read.csv("inat_sunflower.csv") %>%
  rowwise() %>%
  mutate(date = strsplit(time_observed_at, " ")[[1]][1]) %>%
  mutate(year = lubridate::year(date)) %>%
  mutate(month = lubridate::month(date)) %>%
  mutate(week = lubridate::week(date)) %>%
  ungroup() %>%
  group_by(week) %>%
  mutate(num_records = n()) %>%
  ungroup()

df_peak_NM <- df_peak %>%
  filter(state == "NM") %>%
  filter(species_id == 224) %>%
  filter(open_flower_estimate != 0.00)

inat_NM <- read.csv("inat_showy_milkweed.csv") %>%
  rowwise() %>%
  mutate(date = strsplit(time_observed_at, " ")[[1]][1]) %>%
  mutate(year = lubridate::year(date)) %>%
  mutate(month = lubridate::month(date)) %>%
  mutate(week = lubridate::week(date)) %>%
  ungroup() %>%
  group_by(week) %>%
  mutate(num_records = n()) %>%
  ungroup()

df_peak_LA <- df_peak %>%
  filter(state == "LA") %>%
  filter(species_id == 201) %>%
  filter(open_flower_estimate != 0.00)

inat_LA <- read.csv("inat_buttonbush.csv") %>%
  rowwise() %>%
  mutate(date = strsplit(time_observed_at, " ")[[1]][1]) %>%
  mutate(year = lubridate::year(date)) %>%
  mutate(month = lubridate::month(date)) %>%
  mutate(week = lubridate::week(date)) %>%
  ungroup() %>%
  group_by(week) %>%
  mutate(num_records = n()) %>%
  ungroup()

df_peak_OK <- df_peak %>%
  filter(state == "OK") %>%
  filter(species_id == 781) %>%
  filter(open_flower_estimate != 0.00)

inat_OK <- read.csv("inat_silver_maple.csv") %>%
  rowwise() %>%
  mutate(date = strsplit(time_observed_at, " ")[[1]][1]) %>%
  mutate(year = lubridate::year(date)) %>%
  mutate(month = lubridate::month(date)) %>%
  mutate(week = lubridate::week(date)) %>%
  ungroup() %>%
  group_by(week) %>%
  mutate(num_records = n()) %>%
  ungroup()

df_LA_single <- filter(df_LA, phenophase_id == 501, species_id == 201, state == "LA")
df_TX_single <- filter(df_TX, phenophase_id == 501, species_id == 203, state == "TX")
df_NM_single <- filter(df_NM, phenophase_id == 501, species_id == 224, state == "NM")
df_OK_single <- filter(df_OK, phenophase_id == 501, species_id == 781, state == "OK")


#create plot
# These are hard coded and the script must be edited for each state.
# If editing, change any postal abbreviations as well as years in each labs arg.
p1 = ggplot() +
  geom_tile(data = df_TX_single, aes(x=week, y=phenophase_description, fill=proportion)) +
  theme(axis.text.x = element_text(hjust=-.7), axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), legend.position="left", panel.background = element_blank()) +
  scale_fill_gradient(low="grey90", high=df_TX_single$color[1]) +
  labs(x="", y="Proportion of\nOpen Flowers\n2011-2023") +
  scale_x_continuous(breaks=seq(1, 52, by=4.25), limits=c(1,52),
                     labels=c("Jan","Feb","Mar", "Apr", "May", "June", "July", "Aug", "Sep", "Oct", "Nov", "Dec", "")) 
  

p2 = ggplot() +
  geom_line(data = df_peak_TX, aes(x=week, y=open_flower_estimate), color="#3CBB75FF", linewidth = 1) +
  theme(axis.text.y = element_text(), axis.text.x = element_text(hjust=-.7), legend.position = "none",
        panel.background = element_rect("white"), panel.grid = element_blank()) +
  labs(x="", y="Number of\nOpen Flowers\n2014-2022") +
  scale_x_continuous(breaks=seq(1, 52, by=4.25), limits=c(1,52),
                     labels=c("Jan","Feb","Mar", "Apr", "May", "June", "July", "Aug", "Sep", "Oct", "Nov", "Dec", "")) 
  

p3 = ggplot() +
  geom_line(data = inat_TX, aes(x=week, y=num_records), color="#440154FF", linewidth = 1) +
  theme(axis.text.y = element_text(), axis.text.x = element_text(hjust=-.7), legend.position="none",
        panel.background = element_rect("white"), panel.grid = element_blank()) +
  labs(x="", y="Number of\niNat Records\n2011-2023", title=df_peak_TX$common_name[1]) +
  scale_x_continuous(breaks=seq(1, 52, by=4.25), limits=c(1,52),
                     labels=c("Jan","Feb","Mar", "Apr", "May", "June", "July", "Aug", "Sep", "Oct", "Nov", "Dec", ""))
  

wrap_elements(p3 / p2 / p1 & theme(plot.margin = margin(l=1)))
#use the copy plot, and resizing to get it to come out long and narrow
#maybe in the future explore the reshape library to make it circular like the year?


# SITE SPECIFIC CALENDAR
# Site 27474 in LA is good for buttonbush, site 35069 would be only sufficient data in OK,
# but it already constitutes the vast majority of data for silver maple
df_LA_site <- df_LA %>%
  filter(site_id == 27474)
df_peak_LA_site <- df_peak_LA %>%
  filter(site_id == 27474)
  

p4 = ggplot() +
  geom_tile(data = df_LA_site,
            aes(x=week, y=phenophase_description, fill=proportion)) +
  theme(axis.text.x = element_text(hjust=-.7), axis.text.y = element_text(),
        axis.ticks.y = element_blank(), legend.position="left", panel.background = element_blank()) +
  scale_fill_gradient(low="grey90", high=df_LA_site$color[1], (name="")) +
  labs(x="", y="Proportion") +
  scale_x_continuous(breaks=seq(1, 52, by=4.25), limits=c(1,52),
                     labels=c("Jan","Feb","Mar", "Apr", "May", "June", "July", "Aug", "Sep", "Oct", "Nov", "Dec", "")) +
  scale_y_discrete(limits=rev) 

p5 = ggplot() +
  geom_line(data = df_peak_LA_site, aes(x=week, y=open_flower_estimate),color="sienna1", linewidth = 1) +
  theme(axis.text.y = element_text(), axis.text.x = element_text(hjust=-.7)) +
  labs(x="", y="# Open Flowers", title="Jean Lafitte NHPP: Barataria Preserve\nVisitor Center Trail, Louisiana") +
  scale_x_continuous(breaks=seq(1, 52, by=4.25), limits=c(1,52),
                     labels=c("Jan","Feb","Mar", "Apr", "May", "June", "July", "Aug", "Sep", "Oct", "Nov", "Dec", "")) 
p5/p4

# NM SITE SPECIFIC CALENDAR

df_NM_site <- df_NM %>%
  filter(site_id == 16786)
df_peak_NM_site <- df_peak_NM %>%
  filter(site_id == 16786)

p7 = ggplot() +
  geom_tile(data = df_NM_site,
            aes(x=week, y=phenophase_description, fill=proportion)) +
  theme(axis.text.x = element_text(hjust=-.7), axis.text.y = element_text(angle=90, hjust=.5),
        axis.ticks.y = element_blank(), legend.position="left", panel.background = element_blank()) +
  scale_fill_gradient(low="grey90", high=df_NM_site$color[1], (name="")) +
  labs(x="", y="Proportion") +
  scale_x_continuous(breaks=seq(1, 52, by=4.25), limits=c(1,52),
                     labels=c("Jan","Feb","Mar", "Apr", "May", "June", "July", "Aug", "Sep", "Oct", "Nov", "Dec", "")) 

p8 = ggplot() +
  geom_line(data = df_peak_NM_site, aes(x=week, y=open_flower_estimate),color="hotpink", linewidth = 1) +
  theme(axis.text.y = element_text(), axis.text.x = element_text(hjust=-.7)) +
  labs(x="", y="# Open Flowers", title="Cottonwood Gallery 3") +
  scale_x_continuous(breaks=seq(1, 52, by=4.25), limits=c(1,52),
                     labels=c("Jan","Feb","Mar", "Apr", "May", "June", "July", "Aug", "Sep", "Oct", "Nov", "Dec", "")) 
p8/p7



# INTERANNUAL VARIATION
p6 = ggplot() +
  geom_tile(data = filter(df_LA, year != 2014, year != 2017, year != 2023, species_id == 201),
            aes(x=week, y=phenophase_description, fill=proportion)) +
  theme(axis.text.x = element_text(hjust=-.7), axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), legend.position="left", panel.background = element_blank()) +
  scale_fill_gradient(low="grey90", high=df_LA_site$color[1], (name="")) +
  labs(x="", y="Proportion", title='LA common buttonbush') +
  facet_wrap(~year,ncol=1, strip.position="right") +
  scale_x_continuous(breaks=seq(1, 52, by=4.25), limits=c(1,52),
                     labels=c("Jan","Feb","Mar", "Apr", "May", "June", "July", "Aug", "Sep", "Oct", "Nov", "Dec", "")) 
p6