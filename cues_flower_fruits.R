library(rnpn)
library(dplyr)
library(lubridate)
library(ggplot2)
library(leaflet)
library(ggpmisc)
library(ggpubr)
library(ggforce)
library(lme4)
library(mblm)
rm(list=ls())

setwd("~/Documents/My Files/USA-NPN/Data/Analysis/R_default/npn_analyses/TimetoRestore/Data")

#This script performs preliminary analyses on the cues for blooming in pollinator resource
#plants for the Time to Restore project. 

#helper call for species IDs
species<-npn_species()
species <- c(931,916,201,202,203,224,200,204,845,207,781,197,1334,1163,186,1167) # top 16


#download individual phenometric data, with supporting fields, including Daymet climate data (only thru 2021 currently - 2022 coming soon)
#just 2022
df2022 <- npn_download_individual_phenometrics(
  request_source = 'Alyssa', 
  years = c(2022), 
  species_ids = c(paste0(species)),
  additional_fields = c("Site_Name", "Network_Name", "Phenophase_Category", "Observed_Status_Conflict_Flag"),
  phenophase_ids = c(501,390), #open flowers and ripe fruits 
  climate_data = TRUE
)

write.csv(df2022, file="16priority_spp_OF_RF_2022.csv")
df2022 <- (read.csv("16priority_spp_OF_RF_2022.csv"))

#create the list of stations needed to run daymet_download.R
stations <- df2022 %>% distinct(site_id, latitude, longitude)

#pre 2022 other years
df <- npn_download_individual_phenometrics(
  request_source = 'Alyssa', 
  years = c(2009:2021), 
  species_ids = c(paste0(species)),
  additional_fields = c("Site_Name", "Network_Name", "Phenophase_Category", "Observed_Status_Conflict_Flag"),
  phenophase_ids = c(501,390), #open flowers and ripe fruits 
  climate_data = TRUE
)

write.csv(df, file="16priority_spp_OF_RF_2009-2021.csv")
df <- (read.csv("16priority_spp_OF_RF_2009-2021.csv"))

#load daymet data by station reference file (made using daymet_download.R)
daymet2022 <- (read.csv("daymet_ref_file_2022_16spp.csv"))

#merge daymet ref file back with 2022 observational data
df2022_daymet <- merge(df2022, daymet2022, by = c("site_id"), all=TRUE)
colnames(df2022_daymet)

#remove unneeded columns (note - these are for the version from excel with the extra index added)
df2022_daymet <- select(df2022_daymet, -c(2,30:51))
#remove the .y from the merged daymet file col names
colnames(df2022_daymet) <- sub("\\.y","",colnames(df2022_daymet))

#helper for figuring out mismatch between columns in two files
cols_intersection <- intersect(df2022_daymet, df)

#remove unneeded cols in df (gdd, acc_prcp, daylength, that we aren't using, to facilitate the bind)
colnames(df)
df <- select(df, -c(29,30,35,36,41,42,47:49)) #web service
#df <- select(df, -c(1,30,31,36,37,42,43,48:50)) #excel

#row bind, FINALLY make a dataset with all the priority species data and daymet data for 2022
df_complete <- rbind(df2022_daymet, df)

write.csv(df_complete, file="16priority_spp_OF_RF_2009-2022.csv")
df_complete <- (read.csv("16priority_spp_OF_RF_2009-2022.csv"))

#take only the first yes in each year
df_complete = df_complete %>% 
  group_by(individual_id,first_yes_year, phenophase_id) %>%
  filter(first_yes_doy == min(first_yes_doy))

#add a column to facilitate merging with peak data
df_complete$ind_pp_year <- paste0(df_complete$individual_id, "_", df_complete$phenophase_id, "_", df_complete$first_yes_year)

#read in the flower and fruit intensity data created using peak_phenometrics_flower.R and peak_phenometrics_fruit.R
df_of_peak <- (read.csv("intensity_phenometrics_flower_16priority_spp_2013-2022.csv"))
df_rf_peak <- (read.csv("intensity_phenometrics_fruit_16priority_spp_2013-2022.csv"))
df_peak <- rbind(df_of_peak,df_rf_peak)

df_peak$ind_pp_year <- paste0(df_peak$individual_id, "_", df_peak$phenophase_id, "_", df_peak$year)

df_complete_peak <- merge(df_complete, df_peak, by = c("ind_pp_year"), all=TRUE)

colnames(df_complete_peak)
df_complete_peak <- select(df_complete_peak, -c(2,43:46,52,53))
colnames(df_complete_peak) <- sub("\\.x","",colnames(df_complete_peak))

write.csv(df_complete_peak, file="16priority_spp_flower_fruit_w_peak_2009-2022.csv")

#YOU CAN START HERE
#if you want to reuse the same file with daymet integrated, TTR 16 spp, peak and ind pmetrics
df <- (read.csv("16priority_spp_flower_fruit_w_peak_2009-2022.csv"))

#data cleaning
df = df %>% 
  subset(state != -9999) %>% #dropping no-state records - removes Midway Atoll Verbesina records w no climate data
  subset(numdays_since_prior_no != -9999) %>% #dropping no prior no
  subset(numdays_since_prior_no < 14) #dropping if prior no > 14 days prior

#mapping data distribution by species
pal <- colorFactor(palette = c('red','grey','yellow','purple','pink', 'orange','green','blue'), domain = df$common_name)

leaflet(df) %>% 
  addTiles() %>% 
  addLegend("bottomright", pal = pal, values = ~common_name, title = "Species", opacity = 1) %>% 
  addCircles(lng = ~df$longitude, lat = ~df$latitude, label = df$common_name, weight = 8, color = ~pal(common_name)) 

#determining onsets by state to decide re making a table with mean and SE and range by state for TTR states
onsets_by_state <- df %>%
  subset(state %in% c("NM", "TX", "OK", "LA"))  %>%
  group_by(common_name, phenophase_description, state) %>%
  summarise(n_first_yes = sum(!is.na(first_yes_doy))) 

#now that we've dropped all records with a prior no greater than 14 days, how does rest of distribution look?
hist(df$numdays_since_prior_no)

#what species are we working with here
length(unique(df$common_name))
print(unique(df$common_name))
print(unique(df$species))

#set site id to factor so it can be used as a random factor in models
df$site_id_factor <- as.factor(df$site_id)

#I did not use Lat bands in the end, but leaving the code in case
#Create five lat bands named by major mid western cities for looking at temp sensitivity by lat
df$lat_bans <- as.factor(cut(df$latitude, c(-Inf,30,35,40,45,Inf), c("Orlando", "Jackson", "OKCity", "Madison", "Bismark")))

#explore records with observer conflicts
#likely better to do this by spp/phenophase, bc some error-prone combos can be missed
conflict_summary <- df %>%
  count(common_name, phenophase_category, observed_status_conflict_flag) %>%
  group_by(common_name, phenophase_category) %>%
  mutate(observed_status_conflict_flag=recode(
    observed_status_conflict_flag,'MultiObserver-StatusConflict'='Multi', 'OneObserver-StatusConflict'='One')) %>%
  mutate(Percent_Conflict = n / sum(n))

#decide not to remove any conflicting records, because represent less than 4% of records, thus Yes will override a No in all cases (for now)

#determine number of onsets for each spp-phenophase (if you also group by site_id, can calc site_years, only dif if multiple onsets/year, but we already chose first onset only above)
onsets_n_records <- df %>%
  group_by(common_name, phenophase_description) %>% 
  summarise(n_first_yes =  sum(!is.na(first_yes_doy))) 

#determine number of peak records (for some reason cannot make this and the above work as one table)
#also use of sum here is unintuitive to me, but length doesn't work, it counts the NAs even if you put the !is.na
peak_n_records <- df %>%
  group_by(common_name, phenophase_description) %>% 
  summarize(n_peak = sum(!is.na(peak_onset_doy))) 

#drop species with insufficient data
df <- df  %>% subset(species_id %in% c(931,201:203,197,781,1163,916))

#determine mean onset day of year by species and phenophase
mean_onsets <- df  %>%
  group_by(common_name, phenophase_description) %>%
  summarize(mean_onset = mean(first_yes_doy))

write.csv(mean_onsets, file="cues_analysis_mean_onsets.csv")

#determine mean peak onset day of year by species and phenophase
mean_peak_onsets <- df  %>%
  group_by(common_name, phenophase_description) %>%
  filter(!is.na(peak_onset_doy)) %>%
  summarize(mean_peak_onset = mean(peak_onset_doy)) 

write.csv(mean_peak_onsets, file="cues_analysis_mean_peak_onsets.csv")

#split the dataset to facilitate histograms 
#creates this list with 16 elements (one for each spp phenophase combo)
s <- split(df, list(df$phenophase_description, df$common_name))

#histograms for each species-phenophase combo 
#why is xlab doubled?

setwd("~/Documents/My Files/USA-NPN/Data/Analysis/R_default/npn_analyses/TimetoRestore/Data/output")

#look at distribution of records over the years and latitude (where and when data collected)
pdf("latitude_histograms.pdf")
par(mfrow = c(2,2))
for (i in c(1:16)) {
  hist(s[[i]]$latitude, xlab = paste(s[[i]]$common_name, " ",s[[i]]$phenophase_description))
}
dev.off()

pdf("year_histograms.pdf")
par(mfrow = c(2,2))
for (i in c(1:16)) {
  histograms1 <- hist(s[[i]]$first_yes_year, xlab = paste(s[[i]]$common_name, " ",s[[i]]$phenophase_description))
}
dev.off()

#look at distributions of DOY for phenophase onset and peak onset
pdf("onset_peak_histograms.pdf")
par(mfcol = c(2,2))
for (i in c(1:16)) {
  hist(s[[i]]$first_yes_doy, xlab = paste(s[[i]]$common_name, " ",s[[i]]$phenophase_description))
  hist(s[[i]]$peak_onset_doy, xlab = paste(s[[i]]$common_name, " ",s[[i]]$phenophase_description))
}
dev.off()

#look at distribution of the length of peak duration (days)
pdf("duration_histograms.pdf")
par(mfcol = c(2,2))
for (i in c(1:16)) {
  hist(s[[i]]$peak_duration, xlab = paste(s[[i]]$common_name, " ",s[[i]]$phenophase_description))
}
dev.off()

#look at distribution of predictor variables for each species-phenophase combo
pdf("climate_variables_histograms.pdf")
par(mfrow = c(2,2))
for (i in c(1:16)) {
  hist(s[[i]]$tmin_fall, xlab = paste(s[[i]]$common_name, " ",s[[i]]$phenophase_description))
  hist(s[[i]]$tmin_winter, xlab = paste(s[[i]]$common_name, " ",s[[i]]$phenophase_description))
  hist(s[[i]]$tmin_spring, xlab = paste(s[[i]]$common_name, " ",s[[i]]$phenophase_description))
  hist(s[[i]]$tmin_summer, xlab = paste(s[[i]]$common_name, " ",s[[i]]$phenophase_description))
  hist(s[[i]]$tmax_fall, xlab = paste(s[[i]]$common_name, " ",s[[i]]$phenophase_description))
  hist(s[[i]]$tmax_winter, xlab = paste(s[[i]]$common_name, " ",s[[i]]$phenophase_description))
  hist(s[[i]]$tmax_spring, xlab = paste(s[[i]]$common_name, " ",s[[i]]$phenophase_description))
  hist(s[[i]]$tmax_summer, xlab = paste(s[[i]]$common_name, " ",s[[i]]$phenophase_description))
  hist(s[[i]]$prcp_fall, xlab = paste(s[[i]]$common_name, " ",s[[i]]$phenophase_description))
  hist(s[[i]]$prcp_winter, xlab = paste(s[[i]]$common_name, " ",s[[i]]$phenophase_description))
  hist(s[[i]]$prcp_spring, xlab = paste(s[[i]]$common_name, " ",s[[i]]$phenophase_description))
  hist(s[[i]]$prcp_summer, xlab = paste(s[[i]]$common_name, " ",s[[i]]$phenophase_description))
}
dev.off()

#outlier removal
#figure out which number is which spp-pp combo
summary(s)
#removing ripe fruits before day 200 for wild bergamot (no spring fruits)
s[[16]] <- s[[16]]  %>% subset(s[[16]]$first_yes_doy > 200)
#removing open flowers in silver maple after day 250 (no fall flowers)
s[[11]] <- s[[11]]  %>% subset(s[[11]]$first_yes_doy < 250)
#removing ripe fruits in silver maple after day 250 (no fall fruits)
s[[12]] <- s[[12]]  %>% subset(s[[12]]$first_yes_doy < 200)
#removing open flowers in Gooding's willow after day 250 (no fall flowers)
s[[9]] <- s[[9]]  %>% subset(s[[9]]$first_yes_doy < 250)
#removing ripe fruit onset and ripe fruit peak in common buttonbush before day 150 (no spring fruits)
s[[4]] <- s[[4]]  %>% subset(s[[4]]$first_yes_doy > 150)
s[[4]] <- s[[4]]  %>% subset(s[[4]]$peak_onset_doy > 150)

#Simple Linear Regression
#plot a linear model of first day that open flowers or ripe fruits were observed against climate variables


################################################################################
# Create multiple ggplots for each species with a variety of combinations of 
# variables (Jeff Oliver's code)

#' ggplot object of a linear model with a single predictor and single response
#' 
#' @param species_df data frame with values for a single species
#' @param x character indicating predictor variable
#' @param y character indicating response variable
plot_lm <- function(species_df, x, y) {
  species_name <- species_df$common_name[1]
  phenophase <- species_df$phenophase_description[1]
  # Have to use esoteric ggplot syntax first converting x to a symbol, then 
  # unquoting it with !!
  species_plot <- ggplot(data = species_df,
                         mapping = aes(x = !!sym(x), y = !!sym(y))) +
    stat_cor() +
    geom_point() +
    stat_smooth(method = "lm", formula = y ~ x, linewidth = 1) +
    labs(title = paste(species_name, " ", phenophase))
  return(species_plot)
}
# Test of above function
cardinal_flower_plot <- plot_lm(species_df = s[[1]],
                                x = "tmin_fall",
                                y = "first_yes_doy")
print(cardinal_flower_plot)

# Now we want this over:
#    All species-phenophase combos (N = 16)
#    All predictors (N = 12)
#    All responses (N = 3)
# Total number of plots: 16 x 12 x 3 = 576
# Make a vector of the 12 predictors
pvars <- c("tmin", "tmax", "prcp")
pseasons <- c("fall", "winter", "spring", "summer")
predictors <- paste(rep(pvars, each = length(pseasons)), 
                    pseasons, 
                    sep = "_")
# Make a vector of the three response variables
responses <- c("first_yes_doy", "peak_onset_doy", "peak_duration")

# Make sure the output folder exists; if not, create it
if (!dir.exists("output/simple_linear_models")) {
  dir.create("output/simple_linear_models")
}

# We iterate over all species in the list s, creating one plot at a time, 
# saving each ggplot object to a list. After all pairwise predictor/response 
# plots are made, we write these plots to a multi-page pdf - one plot per page
# and one file for each species/phenophase combination
for (species_i in 1:length(s)) {
  # A list to hold each plot for this species/phenophase
  species_plots <- vector(mode = "list")
  # Pull out the data frame for this one species/phenophase
  species_df <- s[[species_i]]
  # Iterate over each response
  for (response in responses) {
    # Iterate over each predictor
    for (predictor in predictors) {
      one_plot <- plot_lm(species_df = species_df,
                          x = predictor,
                          y = response)
      # Add this new plot to the end of the list
      species_plots[[(length(species_plots) + 1)]] <- one_plot
    }
  }
  # Done creating the plots, now save to pdf; we make a "nice name" variable 
  # that doesn't contain spaces or apostrophes (I'm looking at you S. gooddinii)
  nice_name <- gsub(pattern = " ",
                    replacement = "_",
                    x = species_df$common_name[1])
  nice_name <- gsub(pattern = "'",
                    replacement = "",
                    x = nice_name)
  nice_phenophase <- tolower(gsub(pattern = " ",
                             replacement = "_",
                             x = species_df$phenophase_description[1]))
  message("Saving: ", nice_name, ", ", nice_phenophase)
  filename <- paste0("output/simple_linear_models", nice_name, "-", nice_phenophase, "-linear-models.pdf")
  pdf(filename)
  # Note this produces warnings if rows are removed because they contain 
  # missing (NA) values.
  print(species_plots)
  dev.off()
}

# End Jeff solves my ggplot issues
################################################################################

#data format for building linear mixed effect models
#subset to the species / phenophase combo you are working with below 
df1 <- subset(df, species_id == 201  & phenophase_id == 390 & peak_onset_doy > 150) #remember to repeat any cutting of fall flowering etc as above with S object


#relevant predictors tracking - https://docs.google.com/spreadsheets/d/1vknYKsH1cqDSJGtwZIaRp1I3iFCjL55084kmbvJ_noI/edit#gid=1084297830

#Linear Mixed Effect Regression
#Individual plants at the same site are likely to behave similarly - site should probably be a random effect
#Plants are probably responding to multiple climatic conditions, daylength, and these responses may differ
#by latitude based on local adaptation or other factors.
#Compare models based on lowest REML or ML criterion at convergence, and using anova compare more to less complex models

#help for this section https://ourcodingclub.github.io/tutorials/mixed-models/#what and
#https://webinar-free-download.s3.amazonaws.com/COSA-Fixed_and_Random_Factors_in_Mixed_Models-Slides-Handout.pdf
# and https://mspeekenbrink.github.io/sdam-r-companion/linear-mixed-effects-models.html

#look at how variables range at sites - helps to decide if site should indeed be a random effect
par(mfrow = c(4,4))
for (i in c(1:16)) {
  boxplot(s[[i]]$first_yes_doy ~ s[[i]]$site_id, data = s[[i]], xlab = paste(s[[i]]$common_name, " ",s[[i]]$phenophase_description))
  boxplot(s[[i]]$latitude ~ s[[i]]$site_id, data = s[[i]], xlab = paste(s[[i]]$common_name, " ",s[[i]]$phenophase_description))
  boxplot(s[[i]]$prcp_fall ~ s[[i]]$site_id, data = s[[i]], xlab = paste(s[[i]]$common_name, " ",s[[i]]$phenophase_description))
  boxplot(s[[i]]$tmin_winter ~ s[[i]]$site_id, data = s[[i]], xlab = paste(s[[i]]$common_name, " ",s[[i]]$phenophase_description))
}

#object df1 created above
#and manually building these linear models
#as needed - review boxplots for effect of site
boxplot(peak_duration~site_id, data = df1)

#see cues_lmer_model-by-model.R for the code as run for each spp*pp combo
#model 1, first yes predicted by X, with site as a random effect
lmer1 <- lmer(peak_duration~tmax_spring + (1| site_id_factor), data = df1)
summary(lmer1)

#REML criterion at convergence: 114
#(compare these across models - looking for lowest)
#proportion of variance from sites- 1/2
#(you get this from the Random Effects section of the model summary - Random Effects - the site_id variance vs the site_id + residual variance
#(ie, of all the variance how much is attributable to site)

par(mfrow = c(1,1))
#Q-Q plots to check residuals
qqnorm(resid(lmer1))
qqline(resid(lmer1))  

#model 2
lmer2 <- lmer(peak_duration~tmax_spring + tmin_fall + (1| site_id_factor), data = df1)
summary(lmer2)

#REML criterion at convergence: 110
#proportion of variance from sites - 40%

#Q-Q plots to check residuals
qqnorm(resid(lmer2))
qqline(resid(lmer2)) 

#is the more complex model worth the extra parameters?
#make the more complex model the first argument
anova(lmer2, lmer1)

#adding fall tmin does not improve the model (P ChiSq 0.44) 

#model 3
lmer3 <- lmer(peak_duration~tmax_spring + prcp_summer +(1| site_id_factor), data = df1)
summary(lmer3)

#REML criterion at convergence: 119
#proportion of variance from sites 40%

#Q-Q plots to check residuals
qqnorm(resid(lmer3))
qqline(resid(lmer3))


anova(lmer3, lmer1)

#adding summer precip does not improve the model - p is 0.6

#model 4
lmer4 <- lmer(peak_duration~tmax_spring + latitude +(1| site_id_factor), data = df1)
summary(lmer4)

#REML criterion at convergence: 109
#proportion of variance from sites - 40%

#Q-Q plots to check residuals
qqnorm(resid(lmer4))
qqline(resid(lmer4))

anova(lmer4, lmer1)

#adding latitude does improve P of 0.3

#model 5
lmer5 <- lmer(first_yes_doy~tmax_spring+ tmin_winter*latitude  +(1| site_id_factor), data = df1)
summary(lmer5)

#REML criterion at convergence: 892
#proportion of variance from sites 1/2

#Q-Q plots to check residuals
qqnorm(resid(lmer5))
qqline(resid(lmer5))

anova(lmer5, lmer4)

#adding latitude as interactive w winter tmin does not improve p 0.3



