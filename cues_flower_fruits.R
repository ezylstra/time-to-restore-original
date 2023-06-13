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
df <- (read.csv("16priority_spp_flower_fruit_w_peak_2009-2022.csv"))

df = df %>% 
  subset(state != -9999) %>% #dropping no-state records - removes Midway Atoll Verbesina records w no climate data
  subset(numdays_since_prior_no != -9999) %>% #dropping no prior no
  subset(numdays_since_prior_no < 14) #dropping if prior no > 14 days prior

#now that we've dropped all records with a prior no greater than 14 days, how does rest of distribution look?
hist(df$numdays_since_prior_no)

#what species are we working with here
length(unique(df$common_name))
print(unique(df$common_name))
print(unique(df$species))

#set site id to factor so it can be used as a random factor in models
df$site_id_factor <- as.factor(df$site_id)

#create five lat bands named by major mid western cities for looking at temp sensitivity by lat
df$lat_bans <- as.factor(cut(df$latitude, c(-Inf,30,35,40,45,Inf), c("Orlando", "Jackson", "OKCity", "Madison", "Bismark")))

#explore records with observer conflicts - likely better to do this by spp/phenophase, bc some error-prone combos can be missed
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

#split the dataset to facilitate histograms and passing data to linear models
#creates this list with 28 elements (one for each spp phenophase combo)
s <- split(df, list(df$phenophase_description, df$common_name))

#histograms for each species-phenophase combo 
#why is xlab doubled?

#look at distribution of records over the years and latitude (where and when data collected)
par(mfrow = c(2,2))
for (i in c(1:16)) {
  hist(s[[i]]$first_yes_year, xlab = paste(s[[i]]$common_name, " ",s[[i]]$phenophase_description))
  hist(s[[i]]$latitude, xlab = paste(s[[i]]$common_name, " ",s[[i]]$phenophase_description))
}

#look at distributions of DOY for phenophase onset and peak onset
par(mfcol = c(2,2))
for (i in c(1:16)) {
  hist(s[[i]]$first_yes_doy, xlab = paste(s[[i]]$common_name, " ",s[[i]]$phenophase_description))
  hist(s[[i]]$peak_onset_doy, xlab = paste(s[[i]]$common_name, " ",s[[i]]$phenophase_description))
}



#look at distribution of the length of peak duration (days)
par(mfcol = c(2,2))
for (i in c(1:16)) {
  hist(s[[i]]$peak_duration, xlab = paste(s[[i]]$common_name, " ",s[[i]]$phenophase_description))
}

#look at distribution of predictor variables for each species-phenophase combo
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

#outliers
#removing peaks before day 200 for wild bergamot ripe fruit
s[[16]] <- s[[16]]  %>% subset(s[[16]]$first_yes_doy > 200)


#Simple Linear Regression
#plot a linear model of first day that open flowers or ripe fruits were observed against climate variables

#function for x and y
fun <- function(x,y) {
  p <-  ggplot(data = s[[i]], aes(x = x, y = y)) +
    stat_cor() +
    geom_point() +
    stat_smooth(method = "lm", formula = y~x , linewidth = 1) +
    labs(title = paste(s[[i]]$common_name, " ",s[[i]]$phenophase_description)) 
  plot(p)
}

#works - but only allows evaluating one predictor across the whole dataset, really
#i want to look sequentially at each species-phenophase combo, first at all predictors by first_yes_doy
#then all predictors by peak_onset_doy, then all predictors by peak_duration
for (i in c(1:16)) {
  y = s[[i]]$first_yes_doy
  x = s[[i]]$tmax_fall
  fun(x,y)
}

#doesn't work (using c(1:2) so it runs/fails faster)
predictors <- as.list(s[[i]]$tmax_fall, s[[i]]$tmax_winter, s[[i]]$tmax_spring, s[[i]]$tmax_summer)
for (i in c(1:2)) {
  for (x in c(predictors)) {
    y = s[[i]]$first_yes_doy
    fun(x,y)
  }
}

#doesn't work
for (i in c(1:2)) { 
  for (x in c(s[[i]]$tmax_fall, s[[i]]$tmax_winter, s[[i]]$tmax_spring)) {
    y = s[[i]]$first_yes_doy
      fun(x,y)
  }
}

#this doesn't work bc there are no longer colnames in the s object
for (i in c(1:28)) { 
  for (x in colnames(s[[i]])[31:42]) {
    for (y in colnames(s[[i]])[c(22,41,46)]) {
      fun(x,y)
    }
  }
}

#This makes a 4K page PDF with just a single point in the middle of each graph
#Create PDF
pdf(paste("SimpleLinearModels",".pdf",sep=""))
for (i in c(1:2)) { 
  for (y in c(s[[i]]$first_yes_doy,s[[i]]$peak_onset_doy)) {
   for (x in c(s[[i]]$tmax_fall,s[[i]]$tmax_winter)) {
    fun(x,y)
   }
  }
}
dev.off()

df1 <- subset(df, species_id == 202 & phenophase_id == 390)


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
par(mfrow = c(4,4))
for (i in c(1:16)) {
  boxplot(s[[i]]$first_yes_doy ~ s[[i]]$site_id, data = s[[i]], xlab = paste(s[[i]]$common_name, " ",s[[i]]$phenophase_description))
  boxplot(s[[i]]$latitude ~ s[[i]]$site_id, data = s[[i]], xlab = paste(s[[i]]$common_name, " ",s[[i]]$phenophase_description))
  boxplot(s[[i]]$prcp_fall ~ s[[i]]$site_id, data = s[[i]], xlab = paste(s[[i]]$common_name, " ",s[[i]]$phenophase_description))
  boxplot(s[[i]]$tmin_winter ~ s[[i]]$site_id, data = s[[i]], xlab = paste(s[[i]]$common_name, " ",s[[i]]$phenophase_description))
}

#object df1 is created below - still manually ggplotting tor review
#and manually building these linear models
#as needed - review boxplots for effect of site
boxplot(peak_duration~site_id, data = df1)

#model 1, first yes predicted by X, with site as a random effect
lmer1 <- lmer(first_yes_doy~tmax_spring + (1| site_id_factor), data = df1)
summary(lmer1)

#REML criterion at convergence: 680
#(compare these across models - looking for lowest)
#proportion of variance from sites- 1/5
#(you get this from the Random Effects section of the model summary - Random Effects - the site_id variance vs the site_id + residual variance
#(ie, of all the variance how much is attributable to site)

par(mfrow = c(1,1))
#Q-Q plots to check residuals
qqnorm(resid(lmer1))
qqline(resid(lmer1))  

#model 2
lmer2 <- lmer(first_yes_doy~tmax_spring + tmax_winter + (1| site_id_factor), data = df1)
summary(lmer2)

#REML criterion at convergence: 678
#proportion of variance from sites - 1/5

#Q-Q plots to check residuals
qqnorm(resid(lmer2))
qqline(resid(lmer2)) 

#is the more complex model worth the extra parameters?
#make the more complex model the first argument
anova(lmer2, lmer1)

#adding winter tmax does not improve the model (P ChiSq 0.2) marginal

#model 3
lmer3 <- lmer(first_yes_doy~tmax_spring + prcp_winter +(1| site_id_factor), data = df1)
summary(lmer3)

#REML criterion at convergence: 686
#proportion of variance from sites 1/5

#Q-Q plots to check residuals
qqnorm(resid(lmer3))
qqline(resid(lmer3))

anova(lmer3, lmer1)

#adding prcp winter (nor summer tmax) does not improve the model - p is 0.9

#model 4
lmer4 <- lmer(first_yes_doy~tmax_spring + latitude +(1| site_id_factor), data = df1)
summary(lmer4)

#REML criterion at convergence: 678
#proportion of variance from sites - 1/5

#Q-Q plots to check residuals
qqnorm(resid(lmer4))
qqline(resid(lmer4))

anova(lmer4, lmer1)

#adding latitude does improve P of 0.7

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
#also tried making lat interactive with spring tmax, also did not improve

#visually explore a model with latitude
p <- ggplot(df1, aes(x=tmin_summer, y=peak_duration, color=lat_bans)) +
  geom_point() 
  #stat_smooth(method = "lm", formula = y~x , linewidth = 1) 
p 

##GGPLOTS for simple linear regression

df1 <- subset(df, species_id == 781 & phenophase_id == 390 & first_yes_doy < 200)

peak_n_records <- df1 %>%
  summarize(n_peak = sum(!is.na(peak_onset_doy))) 


#FOR PHENOPHASE ONSET
ggplot(data = df1, aes(x = tmin_spring, y = first_yes_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = df1, aes(x = tmax_spring, y = first_yes_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = df1, aes(x = tmin_winter, y = first_yes_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = df1, aes(x = tmax_winter, y = first_yes_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = df1, aes(x = tmin_fall, y = first_yes_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = df1, aes(x = tmax_fall, y = first_yes_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = df1, aes(x = tmin_summer, y = first_yes_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = df1, aes(x = tmax_summer, y = first_yes_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = df1, aes(x = prcp_spring, y = first_yes_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = df1, aes(x = prcp_winter, y = first_yes_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = df1, aes(x = prcp_fall, y = first_yes_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = df1, aes(x = prcp_summer, y = first_yes_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)


#FOR PEAK ONSET
ggplot(data = df1, aes(x = tmin_spring, y = peak_onset_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = df1, aes(x = tmax_spring, y = peak_onset_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = df1, aes(x = tmax_winter, y = peak_onset_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = df1, aes(x = tmin_winter, y = peak_onset_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = df1, aes(x = tmax_fall, y = peak_onset_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = df1, aes(x = tmin_fall, y = peak_onset_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = df1, aes(x = tmax_summer, y = peak_onset_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = df1, aes(x = tmin_summer, y = peak_onset_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = df1, aes(x = prcp_winter, y = peak_onset_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = df1, aes(x = prcp_spring, y = peak_onset_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = df1, aes(x = prcp_fall, y = peak_onset_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = df1, aes(x = prcp_summer, y = peak_onset_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)


#FOR PEAK DURATION 
ggplot(data = df1, aes(x = tmin_spring, y = peak_duration)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = df1, aes(x = tmax_spring, y = peak_duration)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = df1, aes(x = tmax_winter, y = peak_duration)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = df1, aes(x = tmin_winter, y = peak_duration)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = df1, aes(x = tmax_fall, y = peak_duration)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = df1, aes(x = tmin_fall, y = peak_duration)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = df1, aes(x = tmax_summer, y = peak_duration)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = df1, aes(x = tmin_summer, y = peak_duration)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = df1, aes(x = prcp_winter, y = peak_duration)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = df1, aes(x = prcp_spring, y = peak_duration)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = df1, aes(x = prcp_fall, y = peak_duration)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = df1, aes(x = prcp_summer, y = peak_duration)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

