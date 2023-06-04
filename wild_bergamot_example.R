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

df <- (read.csv("16priority_spp_flower_fruit_w_peak_2009-2022.csv"))

df <- subset(df, species_id == 931) # wild bergamot

#filter down to onsets with a prior no within 14 days
df = df %>% 
  subset(numdays_since_prior_no != -9999) %>%
  subset(numdays_since_prior_no < 14)

#how does <14 day of distribution look of prior Nos look?
hist(df$numdays_since_prior_no)

#set site id to factor so it can be used as a random factor in models
df$site_id_factor <- as.factor(df$site_id)

#prolly not using - create five lat bands named by major mid western cities for looking at temp sensitivity by lat
df$lat_bans <- as.factor(cut(df$latitude, c(-Inf,30,35,40,45,Inf), c("Orlando", "Jackson", "OKCity", "Madison", "Bismark")))

#explore records with observer conflicts for this species
conflict_summary <- df %>%
  count(common_name, phenophase_category, observed_status_conflict_flag) %>%
  group_by(common_name, phenophase_category) %>%
  mutate(observed_status_conflict_flag=recode(
    observed_status_conflict_flag,'MultiObserver-StatusConflict'='Multi', 'OneObserver-StatusConflict'='One')) %>%
  mutate(Percent_Conflict = n / sum(n))

#decide not to remove any conflicting records, let yes override no - it is 8% of for flowers, maybe revisit

#determine number of onsets for each spp-phenophase (if you also group by site_id, can calc site_years, only dif if multiple onsets/year, but we already chose first onset only when making the base file)
onsets_n_records <- df %>%
  group_by(common_name, phenophase_description) %>% 
  summarize(n_first_yes = length(first_yes_doy)) 

#determine number of peak records (for some reason cannot make this and the above work as one table)
#use of sum here is unintuitive to me, but length doesn't work, it counts the NAs as rows
peak_n_records <- df %>%
  group_by(common_name, phenophase_description) %>% 
  summarize(n_peak = sum(!is.na(peak_onset_doy))) 


#histograms for response variables species-phenophase combo 
hist(subset(df, phenophase_description == 'Open flowers')$first_yes_year)
hist(subset(df, phenophase_description == 'Open flowers')$first_yes_doy)
hist(subset(df, phenophase_description == 'Open flowers')$peak_onset_doy)
hist(subset(df, phenophase_description == 'Open flowers')$peak_duration)
hist(subset(df, phenophase_description == 'Ripe fruits')$first_yes_doy)
hist(subset(df, phenophase_description == 'Ripe fruits')$peak_onset_doy)
hist(subset(df, phenophase_description == 'Ripe fruits')$peak_duration)

#would be great to get the Open Flowers peak onsets in December captured in the model correctly
#as is, they show up as late records, when they are really early/continuous with spring peak
#To address, would create water years by hand - 
#Set day 300 as the start date for the following year (2020 starts day 300 of 2019)
#just for this spp, just for peak flowering

#look at distribution of predictors
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
#look at each potential predictor of a phenophase onset (or peak or duration)
ggplot(data = subset(df, phenophase_description == 'Open flowers'), 
                     aes(x = tmin_spring, y = first_yes_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = subset(df, phenophase_description == 'Open flowers'), 
       aes(x = tmax_spring, y = first_yes_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = subset(df, phenophase_description == 'Open flowers'), 
       aes(x = tmin_winter, y = first_yes_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = subset(df, phenophase_description == 'Open flowers'), 
       aes(x = tmax_winter, y = first_yes_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = subset(df, phenophase_description == 'Open flowers'), 
       aes(x = tmin_fall, y = first_yes_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = subset(df, phenophase_description == 'Open flowers'), 
       aes(x = tmax_fall, y = first_yes_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = subset(df, phenophase_description == 'Open flowers'), 
       aes(x = tmin_summer, y = first_yes_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = subset(df, phenophase_description == 'Open flowers'), 
       aes(x = tmax_summer, y = first_yes_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = subset(df, phenophase_description == 'Open flowers'), 
       aes(x = prcp_spring, y = first_yes_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = subset(df, phenophase_description == 'Open flowers'), 
       aes(x = prcp_winter, y = first_yes_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = subset(df, phenophase_description == 'Open flowers'), 
       aes(x = prcp_fall, y = first_yes_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)

ggplot(data = subset(df, phenophase_description == 'Open flowers'), 
       aes(x = prcp_summer, y = first_yes_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , linewidth = 1)


#looking at all the above ggplots - for Wild Bergamot, Open Flowers, the only simple linear model that has a 
#decent R2 and P value is fall precip, so I note that in col E of the second tab 
#of the results spreadshet- https://docs.google.com/spreadsheets/d/1vknYKsH1cqDSJGtwZIaRp1I3iFCjL55084kmbvJ_noI/edit#gid=1084297830

#LINEAR MIXED EFFECT REGRESSIONS
#Individual plants at the same site are likely to behave similarly - site should probably be a random effect
#Plants are probably responding to multiple climatic conditions, daylength, and these responses may differ
#by latitude based on local adaptation or other factors.
#Compare models based on lowest REML or ML criterion at convergence, and using anova compare more to less complex models

#help pages for this section https://ourcodingclub.github.io/tutorials/mixed-models/#what and
#https://webinar-free-download.s3.amazonaws.com/COSA-Fixed_and_Random_Factors_in_Mixed_Models-Slides-Handout.pdf
# and https://mspeekenbrink.github.io/sdam-r-companion/linear-mixed-effects-models.html

#look at how variables range at sites - helps to decide if site should indeed be a random effect
boxplot(first_yes_doy ~ site_id, data = subset(df, phenophase_description == 'Open flowers'))
boxplot(latitude ~ site_id, data = subset(df, phenophase_description == 'Open flowers'))
boxplot(tmin_spring ~ site_id, data = subset(df, phenophase_description == 'Open flowers'))
boxplot(first_yes_doy ~ state, data = subset(df, phenophase_description == 'Open flowers'))

#model 1, first yes predicted by tmin_spring, with site as a random effect
lmer1 <- lmer(first_yes_doy~prcp_fall + (1| site_id_factor), data = subset(df, phenophase_description == 'Open flowers'))
summary(lmer1)

#REML criterion at convergence: 554.9
#(compare these across models - looking for lowest)
#proportion of variance from sites ~1/2
#(you get this from the Random Effects section of the model summary - Random Effects - the site_id variance vs the site_id + residual variance
#to see - of all the variance how much is attributable to site)

#Q-Q plots to check residuals 
qqnorm(resid(lmer1))
qqline(resid(lmer1))  

#model 2
lmer2 <- lmer(first_yes_doy~prcp_fall + latitude + (1| site_id_factor), data = subset(df, phenophase_description == 'Open flowers'))
summary(lmer2)

#REML criterion at convergence:  551.7
#proportion of variance from sites ~1/2

#Q-Q plots to check residuals
qqnorm(resid(lmer2))
qqline(resid(lmer2)) 

#is the more complex model worth the extra parameters?
#make the more complex model the first argument
anova(lmer2, lmer1)

#adding latitude does not improve the model (P ChiSq 0.82)

#there are no other candidate cues, so I use the effect sizes from lmer1
#to populate the first tab (cells 2H and 2I)
#https://docs.google.com/spreadsheets/d/1vknYKsH1cqDSJGtwZIaRp1I3iFCjL55084kmbvJ_noI/edit#gid=0

#if there were additional variables that popped out in the simple models, 
#i would be testing them in increasingly complex linear mixed effect models
#only adding interactive effects if I a priori think very likely, and already
#tried additive
