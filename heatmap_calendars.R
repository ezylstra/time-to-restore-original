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
rm(list= ls())

#Notes
#- One viz to ask their pref on peak flowering vs open flowers, vs both (assuming fl/buds not relevant)

#- One  viz to ask about uncertainty - will be calendar with all the uncertainty combined with three sub vizes that show interannual, spatial and observer separately

#- One viz to see prefs on heatmap vs activity curve (including one that looks like our proposed one, a bar for open flowers and a curve with average number of flowers, right?)


#this script is to develop the Time to Restore prototype calendars of type, heatmap (Alyssa, Sept 2022)
setwd("~/Documents/My Files/USA-NPN/Data/Analysis/R_default/npn_analyses/TimetoRestore/Data")

#get buttonbush ind pm
#get buttonbush intensity pm
species_list <- npn_species()
#Find Time to Restore - priority species - https://docs.google.com/spreadsheets/d/1WyxsCqZAAATIgZvR5ewNcZkaN2mumVi2d1vFw8I0HJo/edit#gid=0


species <- c(931,916,201,202,203,224,200,204,845,207,781,197,1334,1163,186,1167)

df <- npn_download_status_data(request_source="Alyssa",
                               years=c(2009:2023),
                               species_ids = c(paste0(species)), # species codes
                               state = c("TX", "LA", "OK", "NM"),
                               additional_fields = c("observedby_person_id"),
                               phenophase_ids= c(501), # open flowers
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
df$phenophase_description <- as.factor(df$phenophase_description)

#Ignore for now - filters used to get buttonbush observer variability plot
  filter(individual_id == 141920) %>%
  filter(observedby_person_id < 43389, observedby_person_id > 38000, observedby_person_id != 39173) %>%
  filter(year !=2017, year !=2022)

#calculate the proportion of records in to the total records by week
df1 <- df %>%
  group_by(week, common_name) %>%
  mutate(num_records = n()) %>%
  mutate(sum_pp = sum(phenophase_status))%>%
  mutate(proportion = sum_pp/num_records) %>%
  ungroup()

#Ignore for now - wrangle peak dataset
df_peak <- (read.csv("buttonbush2013-2022_w_peak.csv"))
df_peak <- df_peak %>%
  mutate(week = lubridate::week(observation_date))  %>%
  filter(state == "LA")

df_peak <- df_peak %>%
  group_by(week) %>%
  mutate(num_records = n()) %>%
  mutate(sum_peak = sum(peak))%>%
  mutate(proportion = sum_peak/num_records) %>%
  filter(year !=2017, year !=2022) %>%
  filter(!is.na(proportion)) %>%
  ungroup()

#ignore for now - hm, there are some NAs unexpected in both estimate of open flowers and in proportions - why?
str(df_peak)

# Create data frames based on states' priority species
df_LA <- df1 %>%
  # Filter based on states relevant to selected species for this state
  filter(state == "LA") %>%
  # Filter for species
  filter(species_id == 201 | species_id == 1334) %>%
  # Assign different colors to different species
  mutate(color = case_when(species_id == 201 ~ "green",
                           species_id == 1334 ~ "pink"))
# Repeat for all states
df_OK <- df1 %>%
  filter(state != "NM" ) %>%
  filter(species_id == 781 | species_id == 201 | species_id == 203) %>%
  mutate(color = case_when(species_id == 201 ~ "green",
                           species_id == 781 ~ "red",
                           species_id == 203 ~ "yellow"))
df_NM <- df1 %>%
  filter(state == "NM" | state == "TX") %>%
  filter(species_id == 224 | species_id == 1163 | species_id == 1167 | species_id == 203) %>%
  mutate(color = case_when(species_id == 224 ~ "pink",
                           species_id == 1163 ~ "green",
                           species_id == 1167 ~ "purple",
                           species_id == 203 ~ "yellow"))
df_TX <- df1 %>%
  filter(state == "TX") %>%
  filter(species_id == 201 | species_id == 202 | species_id == 203) %>%
  mutate(color = case_when(species_id == 201 ~ "green",
                           species_id == 202 ~ "purple",
                           species_id == 203 ~ "yellow"))

# Plotting for peak
p = ggplot() +
  geom_tile(data = df_peak, aes(x = week, y = phenophase_description, fill = proportion)) +
  scale_fill_gradient(low = "grey", high = "yellow", (name = "Proportion")) +
  #facet_grid(df$latitude, as.table = FALSE) +
  #facet_grid(df$observedby_person_id) +
  labs(y = "  ", x = " ") +
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank(), legend.position="left")
plot(p)

#community calendar May 2023

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
    theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank(), legend.position="left") +
    labs(x = "", y ="") +
    scale_fill_gradient(low="grey", high=tmp$color[1], (name="Proportion")) +
    # Uncomment following line for color palette, same for all species
    # scale_fill_gradientn(colors=hcl.colors(35, "Emrld"), (name="Proportion")) +
    facet_wrap(tmp$common_name, nrow=1, strip.position="top")
}
# Stack all plots
purrr::reduce(plots, `+`)

# Code for desired plots where species are outlined in different colors, but gradient is uniform
# Delete color input in geom_tile for smooth gradient (no outline)
p <- ggplot(df_TX) +
  geom_tile(aes(x = week, y = phenophase_description, color=common_name, fill = proportion)) +
  # Code for potential coloration fix - not working
  #scale_fill_manual(values = brewer.pal(name="Set1", n = 3), guide="none") +
  #scale_alpha_manual(values=rep(c("red","purple","green", "blue", "yellow"), 7)) +
  scale_fill_gradient(low = "grey", high = "purple", (name = "Proportion")) +
  facet_wrap(df_TX$common_name, ncol=1, strip.position="top") +
  labs(x = "", y ="") +
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank(), legend.position="left") +
  theme(axis.text = element_text(size = 13))
plot(p)

p2 = ggplot() +
  geom_line(data = df_peak, aes(x=week, y=open_flower_estimate),color="orange", size = 2) +
  labs(y = "# Open Flowers ", x = " ") +
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank())
plot(p2) + theme_minimal()

#use the copy plot, and resizing to get it to come out long and narrow
#maybe in the future explore the reshape library to make it circular like the year

#stuff that didn't work to get months labeled
df$week1 <- as.Date(df$week, origin = "2017-01-01")
month_numeric <- as.numeric(format(unique(df$month, format = "%U")))
labels.minor <- c("Jan", "Feb", "Mar", "Apr", "May", "June", "July", "August", "Sept", "Oct", "Nov", "Dec")
label_df = tibble(mon = month(df$month, origin = "2017-01-01", label = T), week = week(df$week))

#stuff that didn't work for month, inside  ggplot
#scale_x_date(date_breaks = "1 month", date_labels = "%b") +
#scale_x_date(breaks = as.Date(df$month, origin="2022-01-01"), labels = df$month) +
