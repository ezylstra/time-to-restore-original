library(rnpn)
library(dplyr)
library(lubridate)
library(ggplot2)
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

#download individual phenometric data, with supporting fields, including Daymet climate data (only thru 2021 currently - 2022 coming soon)
#just 2022
df2022 <- npn_download_individual_phenometrics(
  request_source = 'Alyssa', 
  years = c(2022), 
  species_ids = c(931,916,201,202,203,224,200,204),
  additional_fields = c("Site_Name", "Network_Name", "Phenophase_Category", "Observed_Status_Conflict_Flag"),
  climate_data = TRUE
)

write.csv(df2022, file="8priority_spp_2022.csv")
df2022 <- (read.csv("8priority_spp_2022.csv"))

#create the list of stations needed to run daymet_download.R
stations <- df2022 %>% distinct(site_id, latitude, longitude)

#pre 2022 other years
df <- npn_download_individual_phenometrics(
  request_source = 'Alyssa', 
  years = c(2009:2021), 
  species_ids = c(931,916,201,202,203,224,200,204),
  additional_fields = c("Site_Name", "Network_Name", "Phenophase_Category", "Observed_Status_Conflict_Flag"),
  climate_data = TRUE
)

write.csv(df, file="8priority_spp_2009-2021.csv")
df <- (read.csv("8priority_spp_2009-2021.csv"))

#load daymet data by station reference file (made using daymet_download.R)
daymet2022 <- (read.csv("daymet_ref_file_2022.csv"))

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
df <- select(df, -c(1,30,31,36,37,42,43,48:50))

#row bind, FINALLY make a dataset with all the priority species data and daymet data for 2022
df_complete <- rbind(df2022_daymet, df)

write.csv(df_complete, file="8priority_spp_2009-2022.csv")

df_complete <- (read.csv("8priority_spp_2009-2022.csv"))
df_complete = df_complete %>% 
  subset(phenophase_description == c('Open flowers')) %>% 
  group_by(individual_id,first_yes_year) %>%
  filter(first_yes_doy == min(first_yes_doy))

df_complete$ind_year <- paste0(df_complete$individual_id, "_", df_complete$first_yes_year)

#read in the peak data created using peak_phenometrics_flower.R
df_peak <- (read.csv("intensity_phenometrics_8priority_spp_2013-2022.csv"))
df_peak$ind_year <- paste0(df_peak$individual_id, "_", df_peak$year)

df_more_complete <- merge(df_complete, df_peak, by = c("ind_year"), all=TRUE)

colnames(df_more_complete)
df_more_complete <- select(df_more_complete, -c(2,43:46))
colnames(df_more_complete) <- sub("\\.x","",colnames(df_more_complete))

write.csv(df_more_complete, file="8priority_spp_w_peak_2009-2022.csv")
df_more_complete <- (read.csv("8priority_spp_w_peak_2009-2022.csv"))

df_more_complete = df_more_complete %>% 
  subset(numdays_since_prior_no != -9999) %>%
  subset(numdays_since_prior_no < 14)


#weird to me that this gives you a row count, rather than a sum of peak dates.. but it does what i need
sum(!is.na(df_more_complete$peak_onset_doy))

hist(df_more_complete$numdays_since_prior_no)

df <- df_more_complete %>% subset(species_id == 201) 

#Calls that can be used to see how much data in the project region/lower lats:
#subset(state == c("TX", "LA", "OK", "NM")) 
#subset(latitude < 37)

#subsetting to only first yes in the year - didn't impact model fit/change outliers, so leaving all first yeses
#group_by(first_yes_year, individual_id, phenophase_description) %>%
#filter(first_yes_doy == min(first_yes_doy))

#viz amt of data by year
hist(df$first_yes_year)
hist(df$first_yes_doy)
hist(df$peak_onset_doy)

#set site id to factor so it can be used as a random factor
df$site_id_factor <- as.factor(df$site_id)

#create five lat bands named by major mid western cities for looking at temp sensitivity by lat
df$lat_bans <- as.factor(cut(df$latitude, c(-Inf,30,35,40,45,Inf), c("Orlando", "Jackson", "OKCity", "Madison", "Bismark")))
  
#explore records with observer conflicts
conflict_summary <- df %>%
  count(phenophase_category, observed_status_conflict_flag) %>%
  group_by(phenophase_category) %>%
  mutate(observed_status_conflict_flag=recode(
    observed_status_conflict_flag,'MultiObserver-StatusConflict'='Multi', 'OneObserver-StatusConflict'='One')) %>%
  mutate(Percent_Conflict = n / sum(n))

#decide not to remove any conflicting records, because represent less than 4% of records, thus Yes will override a No

# dropping very warm spring outlier, didn't effect the model, didn't keep
#NoOutlier <- df %>% subset(tmin_spring < 15)

#determine how many site-years there are - this is only dif from # records where there are multiple onsets/year
onsets_summary <- df %>%
  group_by(site_id, phenophase_description) %>%
  summarize(years = n_distinct(first_yes_year)) 

sum(onsets_summary$years)

#weird to me that this gives you a row count, rather than a sum of peak dates.. but it does what i need
sum(!is.na(df$peak_onset_doy))

#look at distribution of data

par(mfrow = c(2,2))
hist(df$first_yes_doy, main = "Day of Year")
hist(df$peak_onset_doy, main = "Peak Onset")
hist(df$peak_duration, main = "Peak Duration")
hist(df$latitude, main = "Latitude")
par(mfrow = c(2, 2))
hist(df$tmax_spring, main = "Spring Tmax")
hist(df$tmax_fall, main = "Fall Tmax")
hist(df$tmax_summer, main = "Summer Tmax")
hist(df$tmax_winter, main = "Winter Tmax")
par(mfrow = c(2, 2))
hist(df$tmin_spring, main = "Spring Tmin")
hist(df$tmin_fall, main = "Fall Tmin")
hist(df$tmin_summer, main = "Summer Tmin")
hist(df$tmin_winter, main = "Winter Tmin")
par(mfrow = c(2, 2))
hist(df$prcp_spring, main = "Spring Precip")
hist(df$prcp_fall, main = "Fall Precip")
hist(df$prcp_summer, main = "Summer Precip")
hist(df$prcp_winter, main = "Winter Precip")

#Simple Linear Regression
#plot a linear model of first day that "open flowers" were observed against climate variables

ggplot(data = df, aes(x = tmin_spring, y = first_yes_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = df, aes(x = tmax_spring, y = first_yes_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = df, aes(x = tmin_winter, y = first_yes_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = df, aes(x = tmax_winter, y = first_yes_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = df, aes(x = tmin_fall, y = first_yes_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = df, aes(x = tmax_fall, y = first_yes_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = df, aes(x = tmin_summer, y = first_yes_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = df, aes(x = tmax_summer, y = first_yes_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = df, aes(x = prcp_spring, y = first_yes_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = df, aes(x = prcp_winter, y = first_yes_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = df, aes(x = prcp_fall, y = first_yes_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = df, aes(x = prcp_summer, y = first_yes_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)


#FOR PEAK
ggplot(data = df, aes(x = tmin_spring, y = peak_onset_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = df, aes(x = tmax_spring, y = peak_onset_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = df, aes(x = tmax_winter, y = peak_onset_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = df, aes(x = tmin_winter, y = peak_onset_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = df, aes(x = tmax_fall, y = peak_onset_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = df, aes(x = tmin_fall, y = peak_onset_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = df, aes(x = tmax_summer, y = peak_onset_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = df, aes(x = tmin_summer, y = peak_onset_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = df, aes(x = prcp_winter, y = peak_onset_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = df, aes(x = prcp_spring, y = peak_onset_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = df, aes(x = prcp_fall, y = peak_onset_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = df, aes(x = prcp_summer, y = peak_onset_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)


#FOR PEAK DURATION 
ggplot(data = df, aes(x = tmin_spring, y = peak_duration)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = df, aes(x = tmax_spring, y = peak_duration)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = df, aes(x = tmax_winter, y = peak_duration)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = df, aes(x = tmin_winter, y = peak_duration)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = df, aes(x = tmax_fall, y = peak_duration)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = df, aes(x = tmin_fall, y = peak_duration)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = df, aes(x = tmax_summer, y = peak_duration)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = df, aes(x = tmin_summer, y = peak_duration)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = df, aes(x = prcp_winter, y = peak_duration)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = df, aes(x = prcp_spring, y = peak_duration)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = df, aes(x = prcp_fall, y = peak_duration)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = df, aes(x = prcp_summer, y = peak_duration)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)


#relevant predictors - https://docs.google.com/spreadsheets/d/1vknYKsH1cqDSJGtwZIaRp1I3iFCjL55084kmbvJ_noI/edit#gid=1084297830

#Linear Mixed Effect Regression
#Individual plants at the same site are likely to behave similarly - site should probably be a random effect
#Plants are probably responding to multiple climatic conditions, daylength, and these responses may differ
#by latitude based on local adaptation or other factors.
#Compare models based on lowest REML or ML criterion at convergence, and using anova compare more to less complex models

#help for this section https://ourcodingclub.github.io/tutorials/mixed-models/#what and
#https://webinar-free-download.s3.amazonaws.com/COSA-Fixed_and_Random_Factors_in_Mixed_Models-Slides-Handout.pdf
# and https://mspeekenbrink.github.io/sdam-r-companion/linear-mixed-effects-models.html

#look at how variables range at sites - helps to decide if site should indeed be a random effect
boxplot(first_yes_doy ~ site_id, data = df)
boxplot(latitude ~ site_id, data = df)
boxplot(tmin_spring ~ site_id, data = df)
boxplot(first_yes_doy ~ state, data = df)


#model 1, first yes predicted by tmin_spring, with site as a random effect
lmer1 <- lmer(peak_onset_doy~tmin_spring + (1| site_id_factor), data = df)
summary(lmer1)

#REML criterion at convergence: 302
#(compare these across models - looking for lowest)
#proportion of variance from sites ~1/4
#(you get this from the Random Effects section of the model summary - Random Effects - the site_id variance vs the site_id + residual variance
#(ie, of all the variance how much is attributable to site)

#Q-Q plots to check residuals
qqnorm(resid(lmer1))
qqline(resid(lmer1))  

#model 2
lmer2 <- lmer(peak_onset_doy~tmin_spring + prcp_summer + (1| site_id_factor), data = df)
summary(lmer2)

#REML criterion at convergence: 300
#proportion of variance from sites ~1/4

#Q-Q plots to check residuals
qqnorm(resid(lmer2))
qqline(resid(lmer2)) 

#is the more complex model worth the extra parameters?
#make the more complex model the first argument
anova(lmer2, lmer1)

#adding summer prcp does  improve the model (P ChiSq 0.95)

#model 3
lmer3 <- lmer(peak_onset_doy~tmin_spring + tmin_winter + (1| site_id_factor), data = df)
summary(lmer3)

#REML criterion at convergence: 295
#proportion of variance from sites ~1/5

#Q-Q plots to check residuals
qqnorm(resid(lmer3))
qqline(resid(lmer3))

anova(lmer3, lmer1)

#adding winter tmin does improve the model (0.8)

#model 4
lmer4 <- lmer(peak_onset_doy~tmin_spring + latitude + (1| site_id_factor), data = df)
summary(lmer4)

#REML criterion at convergence: 295
#proportion of variance from sites ~1/5

#Q-Q plots to check residuals
qqnorm(resid(lmer4))
qqline(resid(lmer4))

anova(lmer4, lmer1)

#adding latitude as additive doesn't improve P of 0.4

#model 5
lmer5 <- lmer(first_yes_doy~tmin_spring + prcp_summer + (1| site_id_factor), data = df)
summary(lmer5)

#REML criterion at convergence: 309
#proportion of variance from sites ~1/5

#Q-Q plots to check residuals
qqnorm(resid(lmer5))
qqline(resid(lmer5))

anova(lmer5, lmer3)

#adding summer precip doesn't help



#visually explore a model with latitude
p <- ggplot(df, aes(x=tmin_spring, y=first_yes_doy, color=lat_bans)) +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1) 
p 
