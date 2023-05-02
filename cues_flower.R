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
  subset(phenophase_category == 'Flowers') %>% 


#not used:
#subsetting to the region yields only 3 records and to lower lats 35(?) 
subset(state == c("TX", "LA", "OK", "NM")) 
subset(latitude < 37)
#subsetting to only first yes in the year didn't impact model fit/change outliers, so leaving all first yeses
group_by(first_yes_year, individual_id, phenophase_description) %>%
filter(first_yes_doy == min(first_yes_doy))

#set site id to factor so it can be used as a random factor
df$site_id_factor <- as.factor(df$site_id)

#create five lat bands named by major mid western cities for looking at temp sensitivity by lat
df$lat_bans <- as.factor(cut(df$latitude, c(-Inf,30,35,40,45,Inf), c("Orlando", "Jackson", "OKCity", "Madison", "Bismark")))

#remedy NAs in climate data (not using gdd or acc_prcp bc correlated w DOY)
df <- df %>% 
  mutate(tmax_spring = as.numeric(tmax_spring)) %>%
  mutate(tmax_fall = as.numeric(tmax_fall)) %>%
  mutate(tmax_summer = as.numeric(tmax_summer)) %>%
  mutate(tmax_winter = as.numeric(tmax_winter)) %>%
  mutate(tmin_spring = as.numeric(tmin_spring)) %>%
  mutate(tmin_fall = as.numeric(tmin_fall)) %>%
  mutate(tmin_summer = as.numeric(tmin_summer)) %>%
  mutate(tmin_winter = as.numeric(tmin_winter)) %>%
  mutate(prcp_spring = as.numeric(prcp_spring)) %>%
  mutate(prcp_summer = as.numeric(prcp_summer)) %>%
  mutate(prcp_winter = as.numeric(prcp_winter)) %>%
  mutate(prcp_fall = as.numeric(prcp_fall)) 

df <- df %>% 
  mutate(tmax_spring = na_if(tmax_spring, -9999)) %>%
  mutate(tmax_fall = na_if(tmax_fall, -9999)) %>%
  mutate(tmax_winter = na_if(tmax_winter, -9999))%>%
  mutate(tmax_summer = na_if(tmax_summer, -9999))%>%
  mutate(tmin_spring = na_if(tmin_spring, -9999)) %>%
  mutate(tmin_fall = na_if(tmin_fall, -9999)) %>%
  mutate(tmin_winter = na_if(tmin_winter, -9999))%>%
  mutate(tmin_summer = na_if(tmin_summer, -9999))%>%
  mutate(prcp_spring = na_if(prcp_spring, -9999)) %>%
  mutate(prcp_summer = na_if(prcp_summer, -9999)) %>%
  mutate(prcp_winter = na_if(prcp_winter, -9999)) %>%
  mutate(prcp_fall = na_if(prcp_fall, -9999))  
  
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

OpenFlowers_NoOutlier <- OpenFlowers %>%
  subset(tmin_spring < 15)

#only 14 records
OpenFlowers_Jackson <- OpenFlowers %>%
  subset(lat_bans == c("Jackson"))

OpenFlowers_OKCity <- OpenFlowers %>%
  subset(lat_bans == c("OKCity"))

#determine how many site-years there are - this is only dif from # records where there are multiple onsets/year
onsets_summary <- OpenFlowers %>%
  group_by(site_id, phenophase_description) %>%
  summarize(years = n_distinct(first_yes_year))

sum(onsets_summary$years)

#look at distribution of data
hist(OpenFlowers$first_yes_doy)

hist(OpenFlowers$tmax_spring)
hist(OpenFlowers$tmax_fall)
hist(OpenFlowers$tmax_summer)
hist(OpenFlowers$tmax_winter)

hist(OpenFlowers$tmin_spring)
hist(OpenFlowers$tmin_fall)
hist(OpenFlowers$tmin_summer)
hist(OpenFlowers$tmin_winter)

hist(OpenFlowers$prcp_spring)
hist(OpenFlowers$prcp_fall)
hist(OpenFlowers$prcp_summer)
hist(OpenFlowers$prcp_winter)

hist(OpenFlowers$latitude)


#Simple Linear Regression
#plot a linear model of first day that "open flowers" were observed against climate variables

ggplot(data = OpenFlowers, aes(x = tmin_spring, y = first_yes_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , size = 1) 

ggplot(data = OpenFlowers, aes(x = tmax_spring, y = first_yes_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , size = 1)

ggplot(data = OpenFlowers, aes(x = tmax_winter, y = first_yes_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , size = 1)

ggplot(data = OpenFlowers, aes(x = tmin_winter, y = first_yes_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , size = 1)

ggplot(data = OpenFlowers, aes(x = tmax_fall, y = first_yes_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , size = 1)

ggplot(data = OpenFlowers, aes(x = tmin_fall, y = first_yes_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , size = 1)

ggplot(data = OpenFlowers, aes(x = tmax_summer, y = first_yes_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , size = 1)

ggplot(data = OpenFlowers, aes(x = tmin_summer, y = first_yes_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , size = 1)

ggplot(data = OpenFlowers, aes(x = prcp_winter, y = first_yes_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , size = 1)

ggplot(data = OpenFlowers, aes(x = prcp_spring, y = first_yes_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , size = 1)

ggplot(data = OpenFlowers, aes(x = prcp_fall, y = first_yes_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , size = 1)

ggplot(data = OpenFlowers, aes(x = prcp_summer, y = first_yes_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , size = 1)

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

#look at how variables range at sites
boxplot(first_yes_doy ~ site_id, data = OpenFlowers)
boxplot(latitude ~ site_id, data = OpenFlowers)
boxplot(tmin_spring ~ site_id, data = OpenFlowers)
boxplot(first_yes_doy ~ state, data = OpenFlowers)


#model 1
lmer1 <- lmer(first_yes_doy~tmin_spring + (1| site_id_factor), data = OpenFlowers)
summary(lmer1)

#REML criterion at convergence: 1432.4
#proportion of variance from sites ~1/3

#Q-Q plots to check residuals
qqnorm(resid(lmer1))
qqline(resid(lmer1))  

#model 2
lmer2 <- lmer(first_yes_doy~tmin_spring + prcp_summer + (1| site_id_factor), data = OpenFlowers)
summary(lmer2)

#REML criterion at convergence: 1437.31763.2
#proportion of variance from sites ~1/3

#Q-Q plots to check residuals
qqnorm(resid(lmer2))
qqline(resid(lmer2)) 

#is the more complex model worth the extra parameters?
#make the more complex model the first argument
anova(lmer2, lmer1)

#adding summer prcp does not improve the model (P ChiSq 0.36)

#model 3
lmer3 <- lmer(first_yes_doy~tmin_spring + tmin_summer + (1| site_id_factor), data = OpenFlowers)
summary(lmer3)

#REML criterion at convergence: 1427.8
#proportion of variance from sites ~1/3

#Q-Q plots to check residuals
qqnorm(resid(lmer3))
qqline(resid(lmer3))

anova(lmer3, lmer1)

#adding summer temp also does not improve, P chisq 0.36

#model 4
lmer4 <- lmer(first_yes_doy~tmin_spring*latitude + (1| site_id_factor), data = OpenFlowers)
summary(lmer4)

#REML criterion at convergence: 1415.7
#proportion of variance from sites ~1/3

#Q-Q plots to check residuals
qqnorm(resid(lmer4))
qqline(resid(lmer4))

anova(lmer4, lmer1)

#interaction of spring temp and latitude does improve model - P chiSq 0.0005


#model 5
lmer5 <- lmer(first_yes_doy~tmin_spring*lat_bans + (1| site_id_factor), data = OpenFlowers)
summary(lmer5)

#REML criterion at convergence: 1749.9
#proportion of variance from sites ~1/5

#Q-Q plots to check residuals
qqnorm(resid(lmer5))
qqline(resid(lmer5))

anova(lmer5, lmer1)

#interaction of spring temp and latitude bands does improve model - P chiSq 0.000004

#what about this very warm spring 



#visually explore this model
p <- ggplot(OpenFlowers, aes(x=tmin_spring, y=first_yes_doy, color=lat_bans)) +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , size = 1) 
p  

#seems to be driven by low latitude very warm spring, what if remove the outlier in the Jackson bin where tmin_spring >15

p <- ggplot(OpenFlowers_NoOutlier, aes(x=tmin_spring, y=first_yes_doy, color=lat_bans)) +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , size = 1) 
p 

#Ok, so to interpret this, am I looking for the model with the lowest criterion? Which is model 1

#And, then maybe it would be easier to get the effects, like spring min temps advance/delay by X days, the interaction
#effect is like: after accounting for latitude, spring temps delay bloom by X days?

ggplot(data = OpenFlowers_Jackson, aes(x = tmin_spring, y = first_yes_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , size = 1) 

ggplot(data = OpenFlowers_OKCity, aes(x = tmin_spring, y = first_yes_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , size = 1) 
