library(rnpn)
library(dplyr)
library(lubridate)
library(ggplot2)
library(ggpmisc)
library(ggpubr)
library(lme4)
library(mblm)
rm(list=ls())

setwd("~/Documents/My Files/USA-NPN/Data/Analysis/R_default/npn_analyses/TimetoRestore/Data")

#This script performs preliminary analyses on the cues for blooming in pollinator resource
#plants for the Time to Restore project. 

#helper call for species IDs
species<-npn_species()

#download individual phenometric data, with supporting fields, including Daymet climate data (only thru 2021 currently - 2022 coming soon)
df <- npn_download_individual_phenometrics(
  request_source = 'Alyssa', 
  years = c(2009:2023), 
  species_ids = c(202), 
  additional_fields = c("Site_Name", "Network_Name", "Phenophase_Category", "Observed_Status_Conflict_Flag"),
  climate_data = TRUE
)

write.csv(df, file="EPurpConeflower_2009-2023.csv")

df <- (read.csv("EPurpConeflower_2009-2023.csv"))

#subset to flowering data (flowers/flower buds and open flowers)
df = df %>% 
  subset(phenophase_category == 'Flowers') 

#set site id to factor so it can be used as a random factor
df$site_id_factor <- as.factor(df$site_id) 
  
#remedy NAs in climate data
df = df %>% 
  mutate(tmax_spring = na_if(tmax_spring, "-9999")) %>%
  mutate(tmax_fall = na_if(tmax_fall, "-9999")) %>%
  mutate(tmax_winter = na_if(tmax_winter, "-9999"))%>%
  mutate(tmax_summer = na_if(tmax_summer, "-9999"))%>%
  mutate(prcp_spring = na_if(prcp_spring, "-9999")) %>%
  mutate(prcp_summer = na_if(prcp_summer, "-9999")) %>%
  mutate(prcp_winter = na_if(prcp_winter, "-9999")) %>%
  mutate(prcp_fall = na_if(prcp_fall, "-9999"))  %>%
  mutate(acc_prcp = na_if(acc_prcp, "-9999"))
  
#explore records with observer conflicts
conflict_summary <- df %>%
  count(phenophase_category, observed_status_conflict_flag) %>%
  group_by(phenophase_category) %>%
  mutate(observed_status_conflict_flag=recode(
    observed_status_conflict_flag,'MultiObserver-StatusConflict'='Multi', 'OneObserver-StatusConflict'='One')) %>%
  mutate(Percent_Conflict = n / sum(n))

#decide not to remove any conflicting records, because represent less than 1% of records, thus Yes will override a No

#limit to open flowers, and drop NAs
OpenFlowers  <- df %>%
  subset(phenophase_description == c('Open flowers')) %>%
  subset(!is.na(tmax_spring)) 

#determine how many site-years there are - this is only dif from # records where there are multiple onsets/year
#for now, using all onsets, not limiting to first onset of the year
onsets_summary <- OpenFlowers %>%
  group_by(site_id, phenophase_description) %>%
  summarize(years = n_distinct(first_yes_year))

sum(onsets_summary$years)

#look at distribution of variables
hist(OpenFlowers$first_yes_doy)
hist(OpenFlowers$tmax_spring)
hist(OpenFlowers$tmin_spring)
hist(OpenFlowers$acc_prcp)

hist(OpenFlowers$elevation_in_meters)
hist(OpenFlowers$latitude)


#Simple Linear Regression
#plot a linear model of first day that "open flowers" were observed against spring tmin

ggplot(data = OpenFlowers, aes(x = tmin_spring, y = first_yes_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , size = 1)

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

#help for this section https://ourcodingclub.github.io/tutorials/mixed-models/#what and
#https://webinar-free-download.s3.amazonaws.com/COSA-Fixed_and_Random_Factors_in_Mixed_Models-Slides-Handout.pdf
# and https://mspeekenbrink.github.io/sdam-r-companion/linear-mixed-effects-models.html

#look at how variables range at sites
boxplot(first_yes_doy ~ site_id, data = OpenFlowers)
boxplot(latitude ~ site_id, data = OpenFlowers)
boxplot(tmin_spring ~ site_id, data = OpenFlowers)
boxplot(first_yes_doy ~ state, data = OpenFlowers)

#model 1
lmer1 <- lmer(first_yes_doy~(tmin_spring + acc_prcp)*latitude + (1| site_id_factor), data = OpenFlowers)
summary(lmer1)

#REML criterion at convergence: 1617.2
#consider rescaling
#proportion of variance from sites not overwhelming

#Q-Q plots to check residuals
qqnorm(resid(lmer1))
qqline(resid(lmer1))  

#model 2
lmer2 <- lmer(first_yes_doy~tmin_spring + (1| site_id_factor), data = OpenFlowers)
summary(lmer2)

#REML criterion at convergence: 1741.5
#proportion of variance from sites low

#Q-Q plots to check residuals
qqnorm(resid(lmer2))
qqline(resid(lmer2))  


#model 3
lmer3 <- lmer(first_yes_doy~tmin_spring + acc_prcp + (1| site_id_factor), data = OpenFlowers)
summary(lmer3)

#REML criterion at convergence: 1633.8
#proportion of variance from sites ~1/3

#Q-Q plots to check residuals
qqnorm(resid(lmer3))
qqline(resid(lmer3)) 

#model 4
lmer4 <- lmer(first_yes_doy~tmin_spring + acc_prcp + elevation_in_meters + (1| site_id_factor), data = OpenFlowers)
summary(lmer4)

#REML criterion at convergence: 1636.5
#proportion of variance from sites ~1/3

#Q-Q plots to check residuals
qqnorm(resid(lmer4))
qqline(resid(lmer4)) 

#Ok, so to interpret this, am I looking for the model with the lowest criterion? Which is model 1
#So, I should rescale some variables
#And, then maybe it would be easier to get the effects, like spring min temps advance/delay by X days, the interaction
#effect is like: after accounting for latitude, spring temps delay bloom by X days?
