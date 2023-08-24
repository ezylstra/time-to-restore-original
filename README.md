# TimeToRestore

Time to Restore is a SC Climate Adaptation Science Center-funded project aimed at better understanding how the plants that pollinators depend on are responding to climate change. Project results will ultimately help shape species selection for restoration. For more information see www.usanpn.org/TimetoRestore

##Peak Phenometrics
peak_phenometrics_flower.R

peak_phenometrics_fruits.R

peak_phenometrics_leaf.R

These three scripts take as their input Nature's Notebook status data, and format and process the intensity data to create peak phenometrics - estimates of start and end date of the period when the highest number of flowers, fruits or leaves were present on the plant. 

##Cues Analysis

cues_flower_fruits.R

daymet_download.R

cues_lmer_model-by-model.R

The first script is the primary script that compiles, formats and process Nature's Notebook individual phenometric and peak phenometric data to determine temperature and precipitation cues for flowering and fruiting onset, peak and duration. daymet_download.R brings in the missing 2022 daymet climate data to use in the analysis. cues_lmer_model-by-model.R is the workflow tracing, shows how I built and tested each linear mixed effect model for each species and phenophase combo. 

##Calendars and Graphs

heatmap_calendars.R

This script creates Community Calendars, using Nature's Notebook status data to estimate active periods for multiple species and phenophases, with ggplot heat maps. It also creates Species and State Calendars which combine the heatmap calendars with line plots of the estimate of the number of open flowers from Nature's Notebook, as well as the a count of records per week from iNaturalist. 



