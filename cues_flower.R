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

EPC = df_complete %>% subset(phenophase_category == 'Flowers' & species_id == 202) 

#Calls that can be used to see how much data in the project region/lower lats:
#subset(state == c("TX", "LA", "OK", "NM")) 
#subset(latitude < 37)

#subsetting to only first yes in the year - didn't impact model fit/change outliers, so leaving all first yeses
#group_by(first_yes_year, individual_id, phenophase_description) %>%
#filter(first_yes_doy == min(first_yes_doy))

#viz amt of data by year
hist(EPC$first_yes_year)

#set site id to factor so it can be used as a random factor
EPC$site_id_factor <- as.factor(EPC$site_id)

#create five lat bands named by major mid western cities for looking at temp sensitivity by lat
EPC$lat_bans <- as.factor(cut(EPC$latitude, c(-Inf,30,35,40,45,Inf), c("Orlando", "Jackson", "OKCity", "Madison", "Bismark")))
  
#explore records with observer conflicts
conflict_summary <- EPC %>%
  count(phenophase_category, observed_status_conflict_flag) %>%
  group_by(phenophase_category) %>%
  mutate(observed_status_conflict_flag=recode(
    observed_status_conflict_flag,'MultiObserver-StatusConflict'='Multi', 'OneObserver-StatusConflict'='One')) %>%
  mutate(Percent_Conflict = n / sum(n))

#decide not to remove any conflicting records, because represent less than 4% of records, thus Yes will override a No

#limit to open flowers
OpenFlowers  <- EPC %>% subset(phenophase_description == c('Open flowers')) 

# dropping very warm spring outlier, didn't effect the model, didn't keep
#OpenFlowers_NoOutlier <- OpenFlowers %>% subset(tmin_spring < 15)

#determine how many site-years there are - this is only dif from # records where there are multiple onsets/year
onsets_summary <- OpenFlowers %>%
  group_by(site_id, phenophase_description) %>%
  summarize(years = n_distinct(first_yes_year))

sum(onsets_summary$years)

#look at distribution of data

par(mfrow = c(1, 2))
hist(OpenFlowers$first_yes_doy, main = "Day of Year")
hist(OpenFlowers$latitude, main = "Latitude")
par(mfrow = c(2, 2))
hist(OpenFlowers$tmax_spring, main = "Spring Tmax")
hist(OpenFlowers$tmax_fall, main = "Fall Tmax")
hist(OpenFlowers$tmax_summer, main = "Summer Tmax")
hist(OpenFlowers$tmax_winter, main = "Winter Tmax")
par(mfrow = c(2, 2))
hist(OpenFlowers$tmin_spring, main = "Spring Tmin")
hist(OpenFlowers$tmin_fall, main = "Fall Tmin")
hist(OpenFlowers$tmin_summer, main = "Summer Tmin")
hist(OpenFlowers$tmin_winter, main = "Winter Tmin")
par(mfrow = c(2, 2))
hist(OpenFlowers$prcp_spring, main = "Spring Precip")
hist(OpenFlowers$prcp_fall, main = "Fall Precip")
hist(OpenFlowers$prcp_summer, main = "Summer Precip")
hist(OpenFlowers$prcp_winter, main = "Winter Precip")

#Simple Linear Regression
#plot a linear model of first day that "open flowers" were observed against climate variables

ggplot(data = OpenFlowers, aes(x = tmin_spring, y = first_yes_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = OpenFlowers, aes(x = tmax_spring, y = first_yes_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = OpenFlowers, aes(x = tmax_winter, y = first_yes_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = OpenFlowers, aes(x = tmin_winter, y = first_yes_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = OpenFlowers, aes(x = tmax_fall, y = first_yes_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = OpenFlowers, aes(x = tmin_fall, y = first_yes_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = OpenFlowers, aes(x = tmax_summer, y = first_yes_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = OpenFlowers, aes(x = tmin_summer, y = first_yes_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = OpenFlowers, aes(x = prcp_winter, y = first_yes_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = OpenFlowers, aes(x = prcp_spring, y = first_yes_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = OpenFlowers, aes(x = prcp_fall, y = first_yes_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = OpenFlowers, aes(x = prcp_summer, y = first_yes_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

#relevant predictors - tmin_spring, tmax_spring, tmax_winter, tmin_winter, tmax_summer, tmin_summer, tmax_fall, tmin_fall, prcp_fall,prcp_summer

#create a model object for this linear model and then summarize it 
linear_model <- lm(first_yes_doy~tmin_spring, data = OpenFlowers)
summary(linear_model)

plot(linear_model, which = 2)

# check for normally distributed residuals (to meet assumptions, should be normal)
hist(model$residuals, main = "Residuals Histogram")

#check for normality in distribution of residuals with Shapiro Wilks, pvalue should be >0.05 if the distribution is normal
shapiro.test(model$residuals) 
##doesn't pass this test - points to non-parametric tests

#check for patterns in residuals (to meet assumptions, shd be random)
plot(model$residuals, pch = 16, col = "pink", main = "Residuals Plot")


#Kendall–Theil Sen Siegel nonparametric linear regression
#see https://rcompanion.org/handbook/F_12.html - uses Siegel estimator by default

siegel_model = mblm(first_yes_doy~tmin_spring, data=OpenFlowers)
summary(siegel_model)

plot(first_yes_doy~tmax_spring, data=OpenFlowers, pch  = 16)
abline(siegel_model, col="blue",lwd=2)

#median absolute deviation is pretty large, otherwise seems a decent model

#Linear Mixed Effect Regression
#Individual plants at the same site are likely to behave similarly - site should probably be a random effect
#Plants are probably responding to multiple climatic conditions, daylength, and these responses may differ
#by latitude based on local adaptation or other factors.
#Compare models based on lowest REML or ML criterion at convergence, and using anova compare more to less complex models

#help for this section https://ourcodingclub.github.io/tutorials/mixed-models/#what and
#https://webinar-free-download.s3.amazonaws.com/COSA-Fixed_and_Random_Factors_in_Mixed_Models-Slides-Handout.pdf
# and https://mspeekenbrink.github.io/sdam-r-companion/linear-mixed-effects-models.html

#look at how variables range at sites - helps to decide if site should indeed be a random effect
boxplot(first_yes_doy ~ site_id, data = OpenFlowers)
boxplot(latitude ~ site_id, data = OpenFlowers)
boxplot(tmin_spring ~ site_id, data = OpenFlowers)
boxplot(first_yes_doy ~ state, data = OpenFlowers)


#model 1, first yes predicted by tmin_spring, with site as a random effect
lmer1 <- lmer(first_yes_doy~tmin_spring + (1| site_id_factor), data = OpenFlowers)
summary(lmer1)

#REML criterion at convergence: 2040 
#(compare these arcoss models - looking for lowest)
#proportion of variance from sites ~1/3 
#(you get this from the Random Effects section of the model summary - Random Effects - the site_id variance vs the site_id + residual variance
#(ie, of all the variance how much is attributable to site)

#Q-Q plots to check residuals
qqnorm(resid(lmer1))
qqline(resid(lmer1))  

#model 2
lmer2 <- lmer(first_yes_doy~tmin_spring + prcp_winter + (1| site_id_factor), data = OpenFlowers)
summary(lmer2)

#REML criterion at convergence: 2045
#proportion of variance from sites ~1/3

#Q-Q plots to check residuals
qqnorm(resid(lmer2))
qqline(resid(lmer2)) 

#is the more complex model worth the extra parameters?
#make the more complex model the first argument
anova(lmer2, lmer1)

#adding winter prcp does not improve the model (P ChiSq 0.74)

#model 3
lmer3 <- lmer(first_yes_doy~tmin_spring + tmin_winter + (1| site_id_factor), data = OpenFlowers)
summary(lmer3)

#REML criterion at convergence: 2038
#proportion of variance from sites ~1/3

#Q-Q plots to check residuals
qqnorm(resid(lmer3))
qqline(resid(lmer3))

anova(lmer3, lmer1)

#adding winter temp also does not improve, P chisq 0.53

#model 4
lmer4 <- lmer(first_yes_doy~tmin_spring*latitude + (1| site_id_factor), data = OpenFlowers)
summary(lmer4)

#REML criterion at convergence: 2031
#proportion of variance from sites ~1/3

#Q-Q plots to check residuals
qqnorm(resid(lmer4))
qqline(resid(lmer4))

anova(lmer4, lmer1)

#interaction of spring temp and latitude does improve model - P chiSq 0.02
#19 day delay is very unlikely

#model 5 - is using lat_bands defensible?)
lmer5 <- lmer(first_yes_doy~tmin_spring*lat_bans + (1| site_id_factor), data = OpenFlowers)
summary(lmer5)

#REML criterion at convergence: 1984 (lowest of this set of 5)
#proportion of variance from sites ~1/5

#Q-Q plots to check residuals
qqnorm(resid(lmer5))
qqline(resid(lmer5))

anova(lmer5, lmer4)

#interaction of spring temp and latitude bands doesn't improve model, tho not far - P chiSq 0.066
#effect size of 9 delay very unlikely

#visually explore model 4
p <- ggplot(OpenFlowers, aes(x=tmin_spring, y=first_yes_doy, color=lat_bans)) +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1) 
p 


#And, then maybe it would be easier to get the effects, like spring min temps advance/delay by X days, the interaction
#effect is like: after accounting for latitude, spring temps delay bloom by X days?
