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

#code to explore climate drivers for the SC CASC Time to Restore project

species<-npn_species()

df <- npn_download_individual_phenometrics(
  request_source = 'Alyssa', 
  years = c(2009:2023), 
  species_ids = c(202), 
  additional_fields = c("Site_Name", "Network_Name", "Phenophase_Category", "Observed_Status_Conflict_Flag"),
  climate_data = TRUE
)

write.csv(df, file="EPurpConeflower_2009-2023.csv")

df <- (read.csv("EPurpConeflower_2009-2023.csv"))

#Keep only Leaves Category (BLB, Leaves, Colored Leaves, Falling Leaves) and 2010 bc seems they started in Sept, only colored leaves would work, and some unneeded columns
df = df %>% 
  subset(phenophase_category == 'Flowers') 

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
  

#not working - get number of site years
onsets_summary <- df %>%
  group_by(site_id, phenophase_description) %>%
  summarize(distinct_years = n_distinct(df$first_yes_year))
  

#Phenology Data Cleaning
#Conflicting records (see vignette) - go back to the full dataset - df, with 0/1s
conflict_summary <- df %>%
  count(phenophase_category, observed_status_conflict_flag) %>%
  group_by(phenophase_category) %>%
  mutate(observed_status_conflict_flag=recode(
    observed_status_conflict_flag,'MultiObserver-StatusConflict'='Multi', 'OneObserver-StatusConflict'='One')) %>%
  mutate(Percent_Conflict = n / sum(n))

#Decide not to remove any conflicting records, because represent less than 1% of records, thus Yes will override


#Create a dataset of just the open flowers phenophase across species
OpenFlowers  <- df %>%
  subset(phenophase_description == c('Open flowers')) %>%
  subset(!is.na(tmax_spring)) 

hist(OpenFlowers$first_yes_doy)
hist(OpenFlowers$tmax_spring)

#Plot a linear model of first day that "open flowers" were observed against spring tmax
#Note this method is not great, because doesn't account for species and 
#sites, which are likely to behave similarly, random effects model likely better

ggplot(data = OpenFlowers, aes(x = tmax_spring, y = first_yes_doy)) +
  stat_cor() +
  geom_point() +
  stat_smooth(method = "lm", formula = y~x , size = 1)

#create a model object for this linear model and then summarize it 
model <- lm(first_yes_doy~tmax_spring, data = OpenFlowers)
summary(model)

plot(model, which = 2)

boxplot(first_yes_doy ~ site_id, data = OpenFlowers)
boxplot(latitude ~ site_id, data = OpenFlowers)
boxplot(tmax_spring ~ site_id, data = OpenFlowers)
boxplot(first_yes_doy ~ state, data = OpenFlowers)

# look for normally distributed residuals
hist(model$residuals, main = "Residuals Histogram")

#look for pattern, for model assumptions to be met, it should look random
plot(model$residuals, pch = 16, col = "pink", main = "Residuals Plot")

#check for normality in distribution of residuals with Shapiro Wilks, pvalue should be >0.05 if the distribution is normal
shapiro.test(model$residuals) 
#doesn't pass this test, need to ask about this

#help for this section https://ourcodingclub.github.io/tutorials/mixed-models/#what and
#https://webinar-free-download.s3.amazonaws.com/COSA-Fixed_and_Random_Factors_in_Mixed_Models-Slides-Handout.pdf

model2 <- lmer(first_yes_doy~tmax_spring + acc_prcp + (1| site_id_factor), data = OpenFlowers)
summary(model2)

qqnorm(resid(model2))
qqline(resid(model2))  # points fall nicely onto the line - good!

#https://rcompanion.org/handbook/F_12.html - uses Siegals estimator by default


model3 = mblm(first_yes_doy~tmax_spring, data=OpenFlowers)
summary(model3)

plot(first_yes_doy~tmax_spring, data=OpenFlowers, pch  = 16)
abline(model3, col="blue",lwd=2)

#median absolute deviation overlaps 0 - seems like might be tricky
Pvalue    = as.numeric(summary(model3)$coefficients[2,4])
Intercept = as.numeric(summary(model3)$coefficients[1,1])
Slope     = as.numeric(summary(model3)$coefficients[2,1])
R2        = NULL

t1     = paste0("p-value: ", signif(Pvalue, digits=3))
t2     = paste0("R-squared: ", "NULL")
t3     = paste0("Intercept: ", signif(Intercept, digits=3))
t4     = paste0("Slope: ", signif(Slope, digits=3))

