#MODELS FOR CARDINALFLOWER FLOWER PEAK DURATION
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

#MODELS FOR CARDINALFLOWER FLOWER PEAK ONSET
#model 1, first yes predicted by X, with site as a random effect
lmer1 <- lmer(peak_onset_doy~prcp_fall + (1| site_id_factor), data = df1)
summary(lmer1)

#REML criterion at convergence: 117
#(compare these across models - looking for lowest)
#proportion of variance from sites- 6x
#(you get this from the Random Effects section of the model summary - Random Effects - the site_id variance vs the site_id + residual variance
#(ie, of all the variance how much is attributable to site)

par(mfrow = c(1,1))
#Q-Q plots to check residuals
qqnorm(resid(lmer1))
qqline(resid(lmer1))  

#model 2
lmer2 <- lmer(peak_onset_doy~prcp_fall + prcp_winter + (1| site_id_factor), data = df1)
summary(lmer2)

#REML criterion at convergence: 120
#proportion of variance from sites - 6x

#Q-Q plots to check residuals
qqnorm(resid(lmer2))
qqline(resid(lmer2)) 

#is the more complex model worth the extra parameters?
#make the more complex model the first argument
anova(lmer2, lmer1)

#adding winter precip does not improve the model (P ChiSq 0.16) 

#model 3
lmer3 <- lmer(peak_onset_doy~prcp_fall + tmin_spring +(1| site_id_factor), data = df1)
summary(lmer3)

#REML criterion at convergence: 113
#proportion of variance from sites 6x

#Q-Q plots to check residuals
qqnorm(resid(lmer3))
qqline(resid(lmer3))


anova(lmer3, lmer1)

#adding spring tmin does not improve the model - p is 0.2

#model 4
lmer4 <- lmer(peak_onset_doy~prcp_fall + latitude +(1| site_id_factor), data = df1)
summary(lmer4)

#REML criterion at convergence: 112
#proportion of variance from sites - 6x

#Q-Q plots to check residuals
qqnorm(resid(lmer4))
qqline(resid(lmer4))

anova(lmer4, lmer1)

#adding latitude does improve P of 0.08 (marginal)

#MODELS FOR CARDINALFLOWER FLOWER ONSET
#model 1, first yes predicted by X, with site as a random effect
lmer1 <- lmer(first_yes_doy~tmin_summer + (1| site_id_factor), data = df1)
summary(lmer1)

#REML criterion at convergence: 211
#(compare these across models - looking for lowest)
#proportion of variance from sites- 20 x
#(you get this from the Random Effects section of the model summary - Random Effects - the site_id variance vs the site_id + residual variance
#(ie, of all the variance how much is attributable to site)

par(mfrow = c(1,1))
#Q-Q plots to check residuals
qqnorm(resid(lmer1))
qqline(resid(lmer1))  

#model 2
lmer2 <- lmer(first_yes_doy~tmin_summer + tmin_fall + (1| site_id_factor), data = df1)
summary(lmer2)

#REML criterion at convergence: 207
#proportion of variance from sites - 20x

#Q-Q plots to check residuals
qqnorm(resid(lmer2))
qqline(resid(lmer2)) 

#is the more complex model worth the extra parameters?
#make the more complex model the first argument
anova(lmer2, lmer1)

#adding fall tmin does not improve the model (P ChiSq 0.5) 

#model 3
lmer3 <- lmer(first_yes_doy~tmin_summer + prcp_fall +(1| site_id_factor), data = df1)
summary(lmer3)

#REML criterion at convergence: 216
#proportion of variance from sites 20x

#Q-Q plots to check residuals
qqnorm(resid(lmer3))
qqline(resid(lmer3))

#maybe doesn't pass the qq test - pretty messed up!

anova(lmer3, lmer1)

#adding fall precip does not improve the model - p is 0.6

#model 4
lmer4 <- lmer(first_yes_doy~tmin_summer + latitude +(1| site_id_factor), data = df1)
summary(lmer4)

#REML criterion at convergence: 207
#proportion of variance from sites - 20x

#Q-Q plots to check residuals
qqnorm(resid(lmer4))
qqline(resid(lmer4))

#also messed up qq

anova(lmer4, lmer1)

#adding latitude does improve P of 0.5

#MODELS FOR SILVER MAPLE PEAK FRUIT ONSET
#model 1, first yes predicted by X, with site as a random effect
lmer1 <- lmer(peak_onset_doy~tmax_fall + (1| site_id_factor), data = df1)
summary(lmer1)

#REML criterion at convergence: 562
#(compare these across models - looking for lowest)
#proportion of variance from sites- 3/4
#(you get this from the Random Effects section of the model summary - Random Effects - the site_id variance vs the site_id + residual variance
#(ie, of all the variance how much is attributable to site)

par(mfrow = c(1,1))
#Q-Q plots to check residuals
qqnorm(resid(lmer1))
qqline(resid(lmer1))  

#model 2
lmer2 <- lmer(peak_onset_doy~tmax_fall + tmin_spring + (1| site_id_factor), data = df1)
summary(lmer2)

#REML criterion at convergence: 557
#proportion of variance from sites - 3/4

#Q-Q plots to check residuals
qqnorm(resid(lmer2))
qqline(resid(lmer2)) 

#is the more complex model worth the extra parameters?
#make the more complex model the first argument
anova(lmer2, lmer1)

#adding spring tmin does not improve the model (P ChiSq 0.13) 

#model 3
lmer3 <- lmer(peak_onset_doy~tmax_fall + latitude +(1| site_id_factor), data = df1)
summary(lmer3)

#REML criterion at convergence: 558
#proportion of variance from sites 3/4

#Q-Q plots to check residuals
qqnorm(resid(lmer3))
qqline(resid(lmer3))

anova(lmer3, lmer1)

#adding latitude does not improve the model - p is 0.9

#MODELS FOR SILVER MAPLE FRUIT ONSET
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



#MODELS FOR SILVER MAPLE PEAK FLOWER ONSET
#model 1, first yes predicted by X, with site as a random effect
lmer1 <- lmer(peak_onset_doy~tmax_spring + (1| site_id_factor), data = df1)
summary(lmer1)

#REML criterion at convergence: 764
#(compare these across models - looking for lowest)
#proportion of variance from sites-1/4
#(you get this from the Random Effects section of the model summary - Random Effects - the site_id variance vs the site_id + residual variance
#(ie, of all the variance how much is attributable to site)

par(mfrow = c(1,1))
#Q-Q plots to check residuals
qqnorm(resid(lmer1))
qqline(resid(lmer1))  

#model 2
lmer2 <- lmer(peak_onset_doy~tmax_spring + tmin_winter + (1| site_id_factor), data = df1)
summary(lmer2)

#REML criterion at convergence: 760
#proportion of variance from sites - 1/4

#Q-Q plots to check residuals
qqnorm(resid(lmer2))
qqline(resid(lmer2)) 

#is the more complex model worth the extra parameters?
#make the more complex model the first argument
anova(lmer2, lmer1)

#adding winter tmin does not improve the model (P ChiSq 0.6) marginal

#model 3
lmer3 <- lmer(peak_onset_doy~tmax_spring + prcp_winter +(1| site_id_factor), data = df1)
summary(lmer3)

#REML criterion at convergence: 770
#proportion of variance from sites 1/4

#Q-Q plots to check residuals
qqnorm(resid(lmer3))
qqline(resid(lmer3))

anova(lmer3, lmer1)

#adding winter precip does not improve the model - p is 0.6

#model 4
lmer4 <- lmer(peak_onset_doy~tmax_spring + latitude +(1| site_id_factor), data = df1)
summary(lmer4)

#REML criterion at convergence: 761 
#proportion of variance from sites - 1/4

#Q-Q plots to check residuals
qqnorm(resid(lmer4))
qqline(resid(lmer4))

anova(lmer4, lmer1)

#adding latitude does improve P of 0.3

#MODELS FOR SILVER MAPLE FLOWER ONSET

#model 1, first yes predicted by X, with site as a random effect
lmer1 <- lmer(first_yes_doy~tmax_spring + (1| site_id_factor), data = df1)
summary(lmer1)

#REML criterion at convergence: 909
#(compare these across models - looking for lowest)
#proportion of variance from sites-1/2
#(you get this from the Random Effects section of the model summary - Random Effects - the site_id variance vs the site_id + residual variance
#(ie, of all the variance how much is attributable to site)

par(mfrow = c(1,1))
#Q-Q plots to check residuals
qqnorm(resid(lmer1))
qqline(resid(lmer1))  

#model 2
lmer2 <- lmer(first_yes_doy~tmax_spring + tmin_winter + (1| site_id_factor), data = df1)
summary(lmer2)

#REML criterion at convergence: 904
#proportion of variance from sites - 1/2

#Q-Q plots to check residuals
qqnorm(resid(lmer2))
qqline(resid(lmer2)) 

#is the more complex model worth the extra parameters?
#make the more complex model the first argument
anova(lmer2, lmer1)

#adding winter tmin does  improve the model (P ChiSq 0.02)

#model 3
lmer3 <- lmer(first_yes_doy~tmax_spring + tmin_winter + prcp_winter +(1| site_id_factor), data = df1)
summary(lmer3)

#REML criterion at convergence: 910
#proportion of variance from sites 1/2

#Q-Q plots to check residuals
qqnorm(resid(lmer3))
qqline(resid(lmer3))

anova(lmer3, lmer2)

#adding winter precip does not improve the model - p is 0.7

#model 4
lmer4 <- lmer(first_yes_doy~tmax_spring + tmin_winter + latitude +(1| site_id_factor), data = df1)
summary(lmer4)

#REML criterion at convergence: 891
#proportion of variance from sites - 1/2

#Q-Q plots to check residuals
qqnorm(resid(lmer4))
qqline(resid(lmer4))

anova(lmer4, lmer2)

#adding latitude does improve P of 0.0009

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


#MODELS FOR SWAMP MILKWEED RIPE FRUIT PEAK DURATION
boxplot(peak_duration~site_id, data = df1)

#model 1, first yes predicted by X, with site as a random effect
lmer1 <- lmer(peak_duration~prcp_winter + (1| site_id_factor), data = df1)
summary(lmer1)

#REML criterion at convergence: 160
#(compare these across models - looking for lowest)
#proportion of variance from sites- singular

#MODELS FOR SWAMP MILKWEED RIPE FRUIT PEAK ONSET

#model 1, first yes predicted by X, with site as a random effect
lmer1 <- lmer(peak_onset_doy~tmin_spring + (1| site_id_factor), data = df1)
summary(lmer1)

#REML criterion at convergence: 166
#(compare these across models - looking for lowest)
#proportion of variance from sites- 6x
#(you get this from the Random Effects section of the model summary - Random Effects - the site_id variance vs the site_id + residual variance
#(ie, of all the variance how much is attributable to site)

par(mfrow = c(1,1))
#Q-Q plots to check residuals
qqnorm(resid(lmer1))
qqline(resid(lmer1))  

#model 2
lmer2 <- lmer(peak_onset_doy~tmin_spring + tmin_fall + (1| site_id_factor), data = df1)
summary(lmer2)

#REML criterion at convergence: 157
#proportion of variance from sites - 12x

#Q-Q plots to check residuals
qqnorm(resid(lmer2))
qqline(resid(lmer2)) 

#is the more complex model worth the extra parameters?
#make the more complex model the first argument
anova(lmer2, lmer1)

#adding fall tmin does not improve the model (P ChiSq 0.06) - marginal

#model 3
lmer3 <- lmer(peak_onset_doy~tmin_spring + prcp_summer + (1| site_id_factor), data = df1)
summary(lmer3)

#REML criterion at convergence: 169
#proportion of variance from sites -5x

#Q-Q plots to check residuals
qqnorm(resid(lmer3))
qqline(resid(lmer3))

anova(lmer3, lmer1)

#adding summer precip does not improve the model - p is 0.1

#model 4
lmer4 <- lmer(peak_onset_doy~tmin_spring + tmin_summer +(1| site_id_factor), data = df1)
summary(lmer4)

#REML criterion at convergence: 161
#proportion of variance from sites - 6x

#Q-Q plots to check residuals
qqnorm(resid(lmer4))
qqline(resid(lmer4))

anova(lmer4, lmer1)

#adding summer tmin  doesn't improve P of 0.6

#model 5
lmer5 <- lmer(peak_onset_doy~tmin_spring + latitude + (1| site_id_factor), data = df1)
summary(lmer5)

#REML criterion at convergence: 160
#proportion of variance from sites ~10x

#Q-Q plots to check residuals
qqnorm(resid(lmer5))
qqline(resid(lmer5))

anova(lmer5, lmer1)

#adding latitude does improve - 0.4

#MODELS FOR SWAMP MILKWEED RIPE FRUIT ONSET

#model 1, first yes predicted by X, with site as a random effect
lmer1 <- lmer(first_yes_doy~prcp_summer + (1| site_id_factor), data = df1)
summary(lmer1)

#REML criterion at convergence: 303
#(compare these across models - looking for lowest)
#proportion of variance from sites- 1/3
#(you get this from the Random Effects section of the model summary - Random Effects - the site_id variance vs the site_id + residual variance
#(ie, of all the variance how much is attributable to site)

par(mfrow = c(1,1))
#Q-Q plots to check residuals
qqnorm(resid(lmer1))
qqline(resid(lmer1))  

#model 2
lmer2 <- lmer(first_yes_doy~prcp_summer + tmin_winter + (1| site_id_factor), data = df1)
summary(lmer2)

#REML criterion at convergence: 298
#proportion of variance from sites - 1/3

#Q-Q plots to check residuals
qqnorm(resid(lmer2))
qqline(resid(lmer2)) 

#is the more complex model worth the extra parameters?
#make the more complex model the first argument
anova(lmer2, lmer1)

#adding winter tmin does not improve the model (P ChiSq 0.1)

#model 3
lmer3 <- lmer(first_yes_doy~prcp_summer + tmin_spring + (1| site_id_factor), data = df1)
summary(lmer3)

#REML criterion at convergence: 299
#proportion of variance from sites - 1/3

#Q-Q plots to check residuals
qqnorm(resid(lmer3))
qqline(resid(lmer3))

anova(lmer3, lmer1)

#adding spring tmin does not improve the model - p is 0.28

#model 4
lmer4 <- lmer(first_yes_doy~prcp_summer + latitude + (1| site_id_factor), data = df1)
summary(lmer4)

#REML criterion at convergence: 298
#proportion of variance from sites - 40%

#Q-Q plots to check residuals
qqnorm(resid(lmer4))
qqline(resid(lmer4))

anova(lmer4, lmer1)

#adding latitude  doesn't improve P of 0.4


#MODELS FOR SWAMP MILKWEED FLOWER PEAK DURATION

#model 1, first yes predicted by X, with site as a random effect
lmer1 <- lmer(peak_duration~tmin_winter + (1| site_id_factor), data = df1)
summary(lmer1)

#REML criterion at convergence: 232
#(compare these across models - looking for lowest)
#proportion of variance from sites- 1/3
#(you get this from the Random Effects section of the model summary - Random Effects - the site_id variance vs the site_id + residual variance
#(ie, of all the variance how much is attributable to site)

par(mfrow = c(1,1))
#Q-Q plots to check residuals
qqnorm(resid(lmer1))
qqline(resid(lmer1))  

#model 2
lmer2 <- lmer(peak_duration~tmin_winter + tmin_fall + (1| site_id_factor), data = df1)
summary(lmer2)

#REML criterion at convergence: 227
#proportion of variance from sites - 1/3

#Q-Q plots to check residuals
qqnorm(resid(lmer2))
qqline(resid(lmer2)) 

#is the more complex model worth the extra parameters?
#make the more complex model the first argument
anova(lmer2, lmer1)

#adding tmin fall  does not improve the model (P ChiSq 0.4)

#model 3
lmer3 <- lmer(peak_duration~tmin_winter + prcp_summer +  (1| site_id_factor), data = df1)
summary(lmer3)

#REML criterion at convergence: 236
#proportion of variance from sites - 1/6

#Q-Q plots to check residuals
qqnorm(resid(lmer3))
qqline(resid(lmer3))

anova(lmer3, lmer1)

#adding summper precip does not improve the model - p is 0.4

#model 4
lmer4 <- lmer(peak_duration~tmin_winter + latitude +  (1| site_id_factor), data = df1)
summary(lmer4)

#REML criterion at convergence: 227
#proportion of variance from sites - 1/3

#Q-Q plots to check residuals
qqnorm(resid(lmer4))
qqline(resid(lmer4))

anova(lmer4, lmer1)

#adding latitude  doesn't improve P of 0.5

#MODELS FOR SWAMP MILKWEED FLOWER PEAK ONSET

#model 1, first yes predicted by X, with site as a random effect
lmer1 <- lmer(peak_onset_doy~prcp_winter + (1| site_id_factor), data = df1)
summary(lmer1)

#REML criterion at convergence: 236.6
#(compare these across models - looking for lowest)
#proportion of variance from sites- 5x
#(you get this from the Random Effects section of the model summary - Random Effects - the site_id variance vs the site_id + residual variance
#(ie, of all the variance how much is attributable to site)

par(mfrow = c(1,1))
#Q-Q plots to check residuals
qqnorm(resid(lmer1))
qqline(resid(lmer1))  

#model 2
lmer2 <- lmer(peak_onset_doy~prcp_winter + tmax_spring + (1| site_id_factor), data = df1)
summary(lmer2)

#REML criterion at convergence: 234
#proportion of variance from sites - 5x

#Q-Q plots to check residuals
qqnorm(resid(lmer2))
qqline(resid(lmer2)) 

#is the more complex model worth the extra parameters?
#make the more complex model the first argument
anova(lmer2, lmer1)

#adding spring tmax does not improve the model (P ChiSq 0.7)

#model 3
lmer3 <- lmer(peak_onset_doy~prcp_winter + tmin_fall +  (1| site_id_factor), data = df1)
summary(lmer3)

#REML criterion at convergence: 233
#proportion of variance from sites - 5x

#Q-Q plots to check residuals
qqnorm(resid(lmer3))
qqline(resid(lmer3))

anova(lmer3, lmer1)

#adding fall tmin does not improve the model - p is 0.9

#model 4
lmer4 <- lmer(peak_onset_doy~prcp_winter + latitude +  (1| site_id_factor), data = df1)
summary(lmer4)

#REML criterion at convergence: 232
#proportion of variance from sites - 6x

#Q-Q plots to check residuals
qqnorm(resid(lmer4))
qqline(resid(lmer4))

anova(lmer4, lmer1)

#adding latitude  doesn't improve P of 0.5


#MODELS FOR SWAMP MILKWEED FLOWER ONSET

#model 1, first yes predicted by X, with site as a random effect
lmer1 <- lmer(first_yes_doy~tmin_fall + (1| site_id_factor), data = df1)
summary(lmer1)

#REML criterion at convergence: 303
#(compare these across models - looking for lowest)
#proportion of variance from sites - 2x
#(you get this from the Random Effects section of the model summary - Random Effects - the site_id variance vs the site_id + residual variance
#(ie, of all the variance how much is attributable to site)

par(mfrow = c(1,1))
#Q-Q plots to check residuals
qqnorm(resid(lmer1))
qqline(resid(lmer1))  

#model 2
lmer2 <- lmer(first_yes_doy~tmin_fall + tmin_spring + (1| site_id_factor), data = df1)
summary(lmer2)

#REML criterion at convergence: 293
#proportion of variance from sites 2x

#Q-Q plots to check residuals
qqnorm(resid(lmer2))
qqline(resid(lmer2)) 

#is the more complex model worth the extra parameters?
#make the more complex model the first argument
anova(lmer2, lmer1)

#adding spring tmin does improve the model (P ChiSq 0.02)

#model 3 - disregarded
lmer3 <- lmer(first_yes_doy~tmin_fall + tmin_spring + prcp_summer + (1| site_id_factor), data = df1)
summary(lmer3)

#REML criterion at convergence: 290
#proportion of variance from sites 1.3x

#Q-Q plots to check residuals
qqnorm(resid(lmer3))
qqline(resid(lmer3))

anova(lmer3, lmer2)

#adding summer precip does improve the model - p is 0.002
#however upon further review, flowering onset is before summer, so we should not include this predictor

#model 4
lmer4 <- lmer(first_yes_doy~tmin_fall + tmin_spring + latitude + (1| site_id_factor), data = df1)
summary(lmer4)

#REML criterion at convergence: 286
#proportion of variance from sites 2x

#Q-Q plots to check residuals
qqnorm(resid(lmer4))
qqline(resid(lmer4))

anova(lmer4, lmer2)

#adding latitude doesn't improve P of 0.9


#MODELS FOR SUNFLOWER PEAK FRUIT ONSET
#None ID-ed bc all these models are singular and n is 12
#model 1, first yes predicted by X, with site as a random effect
lmer1 <- lmer(peak_onset_doy~prcp_summer + (1| site_id_factor), data = df1)
summary(lmer1)

#REML criterion at convergence: 127


#MODELS FOR SUNFLOWER RIPE FRUIT ONSET
#None ID-ed bc all these models are singular and n is 17
#model 1, first yes predicted by X, with site as a random effect
lmer1 <- lmer(first_yes_doy~tmin_spring + (1| site_id_factor), data = df1)
summary(lmer1)

#REML criterion at convergence: 177
#(compare these across models - looking for lowest)
#proportion of variance from sites -0 - is singular
#(you get this from the Random Effects section of the model summary - Random Effects - the site_id variance vs the site_id + residual variance
#(ie, of all the variance how much is attributable to site)

par(mfrow = c(1,1))
#Q-Q plots to check residuals
qqnorm(resid(lmer1))
qqline(resid(lmer1))  

#model 2
lmer2 <- lmer(first_yes_doy~prcp_summer + (1| site_id_factor), data = df1)
summary(lmer2)

#REML criterion at convergence: 183
#proportion of variance from sites 2x


#MODELS FOR SUNFLOWER PEAK FLOWERING ONSET

#model 1, first yes predicted by X, with site as a random effect
lmer1 <- lmer(peak_onset_doy~tmin_fall + (1| site_id_factor), data = df1)
summary(lmer1)

#REML criterion at convergence: 168
#(compare these across models - looking for lowest)
#proportion of variance from sites -2x
#(you get this from the Random Effects section of the model summary - Random Effects - the site_id variance vs the site_id + residual variance
#(ie, of all the variance how much is attributable to site)

par(mfrow = c(1,1))
#Q-Q plots to check residuals
qqnorm(resid(lmer1))
qqline(resid(lmer1))  

#model 2
lmer2 <- lmer(peak_onset_doy~tmin_fall + tmin_summer + (1| site_id_factor), data = df1)
summary(lmer2)

#REML criterion at convergence: 163
#proportion of variance from sites 2x

#Q-Q plots to check residuals
qqnorm(resid(lmer2))
qqline(resid(lmer2)) 

#is the more complex model worth the extra parameters?
#make the more complex model the first argument
anova(lmer2, lmer1)

#adding summer tmin does not improve the model (P ChiSq 0.6)

#model 3
lmer3 <- lmer(peak_onset_doy~tmin_fall + tmin_spring + (1| site_id_factor), data = df1)
summary(lmer3)

#REML criterion at convergence: 162
#proportion of variance from sites 1/3

#Q-Q plots to check residuals
qqnorm(resid(lmer3))
qqline(resid(lmer3))

anova(lmer3, lmer1)

#adding spring tmin  does not improve the model (0.09) - also get the "is singular warning", maybe over fit?

#model 4
lmer4 <- lmer(peak_onset_doy~tmin_fall + prcp_spring + (1| site_id_factor), data = df1)
summary(lmer4)

#REML criterion at convergence: 172
#proportion of variance from sites 2x

#Q-Q plots to check residuals
qqnorm(resid(lmer4))
qqline(resid(lmer4))

anova(lmer4, lmer1)

#adding spring prcp  doesn't improve P of 0.6

#model 5
lmer5 <- lmer(peak_onset_doy~tmin_fall + latitude + (1| site_id_factor), data = df1)
summary(lmer5)

#REML criterion at convergence: 162
#proportion of variance from sites ~1/2

#Q-Q plots to check residuals
qqnorm(resid(lmer5))
qqline(resid(lmer5))

anova(lmer5, lmer1)

#adding latitude does improve - 0.03

#MODELS FOR SUNFLOWER FLOWERING ONSET

#model 1, first yes predicted by X, with site as a random effect
lmer1 <- lmer(first_yes_doy~tmin_fall + (1| site_id_factor), data = df1)
summary(lmer1)

#REML criterion at convergence: 201
#(compare these across models - looking for lowest)
#proportion of variance from sites -40%
#(you get this from the Random Effects section of the model summary - Random Effects - the site_id variance vs the site_id + residual variance
#(ie, of all the variance how much is attributable to site)

par(mfrow = c(1,1))
#Q-Q plots to check residuals
qqnorm(resid(lmer1))
qqline(resid(lmer1))  

#model 2
lmer2 <- lmer(first_yes_doy~tmin_fall + tmin_summer + (1| site_id_factor), data = df1)
summary(lmer2)

#REML criterion at convergence: 195
#proportion of variance from sites 40%

#Q-Q plots to check residuals
qqnorm(resid(lmer2))
qqline(resid(lmer2)) 

#is the more complex model worth the extra parameters?
#make the more complex model the first argument
anova(lmer2, lmer1)

#adding summer tmin does not improve the model (P ChiSq 0.9)

#model 3
lmer3 <- lmer(first_yes_doy~tmin_fall + tmin_winter + (1| site_id_factor), data = df1)
summary(lmer3)

#REML criterion at convergence: 196
#proportion of variance from sites 1/3

#Q-Q plots to check residuals
qqnorm(resid(lmer3))
qqline(resid(lmer3))

anova(lmer3, lmer1)

#adding winter tmin  does not improve the model (0.3) 

#model 4
lmer4 <- lmer(first_yes_doy~tmin_fall + latitude + (1| site_id_factor), data = df1)
summary(lmer4)

#REML criterion at convergence: 195
#proportion of variance from sites ~40%

#Q-Q plots to check residuals
qqnorm(resid(lmer4))
qqline(resid(lmer4))

anova(lmer4, lmer1)

#adding latitude  doesn't improve P of 0.14

#MODELS FOR BUTTON BUSH FRUIT PEAK DURATION

#model 1, first yes predicted by tmin_spring, with site as a random effect
lmer1 <- lmer(peak_duration~tmin_summer + (1| site_id_factor), data = df1)
summary(lmer1)

#REML criterion at convergence: 244
#(compare these across models - looking for lowest)
#proportion of variance from sites ~1/2
#(you get this from the Random Effects section of the model summary - Random Effects - the site_id variance vs the site_id + residual variance
#(ie, of all the variance how much is attributable to site)

par(mfrow = c(1,1))
#Q-Q plots to check residuals
qqnorm(resid(lmer1))
qqline(resid(lmer1))  

#model 2
lmer2 <- lmer(peak_duration~tmin_summer + latitude + (1| site_id_factor), data = df1)
summary(lmer2)

#REML criterion at convergence: 237
#proportion of variance from sites 1/3

#Q-Q plots to check residuals
qqnorm(resid(lmer2))
qqline(resid(lmer2)) 

#is the more complex model worth the extra parameters?
#make the more complex model the first argument
anova(lmer2, lmer1)

#adding latitude does not improve the model (P ChiSq 0.09) - marginal


#MODELS FOR BUTONBUSH FRUIT PEAK ONSET
#model 1, first yes predicted by tmin_spring, with site as a random effect
lmer1 <- lmer(peak_onset_doy~tmin_summer + (1| site_id_factor), data = df1)
summary(lmer1)

#REML criterion at convergence: 170
#(compare these across models - looking for lowest)
#proportion of variance from sites 50%
#(you get this from the Random Effects section of the model summary - Random Effects - the site_id variance vs the site_id + residual variance
#(ie, of all the variance how much is attributable to site)

par(mfrow = c(1,1))
#Q-Q plots to check residuals
qqnorm(resid(lmer1))
qqline(resid(lmer1))  

#model 2
lmer2 <- lmer(peak_onset_doy~tmin_summer + tmin_spring + (1| site_id_factor), data = df1)
summary(lmer2)

#REML criterion at convergence: 157
#proportion of variance from sites 50%

#Q-Q plots to check residuals
qqnorm(resid(lmer2))
qqline(resid(lmer2)) 

#is the more complex model worth the extra parameters?
#make the more complex model the first argument
anova(lmer2, lmer1)

#adding spring tmin does improve the model (P ChiSq 0.006)
#however, this model has nonsensical effect sizes, and the N is small, so excluding results at this stage

#model 3
lmer3 <- lmer(peak_onset_doy~tmin_summer + tmin_spring + tmin_winter + (1| site_id_factor), data = df1)
summary(lmer3)

#REML criterion at convergence: 152
#proportion of variance from sites 55%

#Q-Q plots to check residuals
qqnorm(resid(lmer3))
qqline(resid(lmer3))

anova(lmer3, lmer2)

#adding winter tmin  does not  improve the model (0.9) 

#model 4
lmer4 <- lmer(peak_onset_doy~tmin_summer + tmin_spring + latitude + (1| site_id_factor), data = df1)
summary(lmer4)

#REML criterion at convergence: 151
#proportion of variance from sites 55%

#Q-Q plots to check residuals
qqnorm(resid(lmer4))
qqline(resid(lmer4))

anova(lmer4, lmer2)

#adding latitide does not improve the model (0.6) 


#MODELS FOR COMMON BUTTONBUSH FRUIT ONSET

lmer1 <- lmer(first_yes_doy~prcp_summer + (1| site_id_factor), data = df1)
summary(lmer1)

#REML criterion at convergence: 296
#(compare these across models - looking for lowest)
#proportion of variance from sites ~95%
#(you get this from the Random Effects section of the model summary - Random Effects - the site_id variance vs the site_id + residual variance
#(ie, of all the variance how much is attributable to site)

par(mfrow = c(1,1))
#Q-Q plots to check residuals
qqnorm(resid(lmer1))
qqline(resid(lmer1))  

#model 2
lmer2 <- lmer(first_yes_doy~prcp_summer + latitude + (1| site_id_factor), data = df1)
summary(lmer2)

#REML criterion at convergence: 287
#proportion of variance from sites 95%

#Q-Q plots to check residuals
qqnorm(resid(lmer2))
qqline(resid(lmer2)) 

#is the more complex model worth the extra parameters?
#make the more complex model the first argument
anova(lmer2, lmer1)

#adding latitude does  improve the model (P ChiSq 0.03)

#model 3
lmer3 <- lmer(first_yes_doy~prcp_summer*latitude + (1| site_id_factor), data = df1)
summary(lmer3)

#REML criterion at convergence: 287
#proportion of variance from sites 95%

#Q-Q plots to check residuals
qqnorm(resid(lmer3))
qqline(resid(lmer3))

anova(lmer3, lmer2)

#adding latitude as interactive does improve the model (0.03) 


#MODELS FOR EPC FRUIT PEAK DURATION

#model 1, first yes predicted by tmin_spring, with site as a random effect
lmer1 <- lmer(peak_duration~tmin_summer + (1| site_id_factor), data = df1)
summary(lmer1)

#REML criterion at convergence: 245
#(compare these across models - looking for lowest)
#proportion of variance from sites ~1/5
#(you get this from the Random Effects section of the model summary - Random Effects - the site_id variance vs the site_id + residual variance
#(ie, of all the variance how much is attributable to site)

par(mfrow = c(1,1))
#Q-Q plots to check residuals
qqnorm(resid(lmer1))
qqline(resid(lmer1))  

#model 2
lmer2 <- lmer(peak_duration~tmin_summer + tmin_fall + (1| site_id_factor), data = df1)
summary(lmer2)

#REML criterion at convergence: 239
#proportion of variance from sites ~1/2

#Q-Q plots to check residuals
qqnorm(resid(lmer2))
qqline(resid(lmer2)) 

#is the more complex model worth the extra parameters?
#make the more complex model the first argument
anova(lmer2, lmer1)

#adding fall tmin does not improve the model (P ChiSq 0.6)

#model 3
lmer3 <- lmer(peak_duration~tmin_summer + latitude + (1| site_id_factor), data = df1)
summary(lmer3)

#REML criterion at convergence: 237
#proportion of variance from sites ~1/3

#Q-Q plots to check residuals
qqnorm(resid(lmer3))
qqline(resid(lmer3))

anova(lmer3, lmer1)

#adding latitude does not improve the model (0.09) - borderline tho



#MODELS FOR EPC FRUIT PEAK ONSET

#model 1, first yes predicted by tmin_spring, with site as a random effect
lmer1 <- lmer(peak_onset_doy~tmax_spring + (1| site_id_factor), data = df1)
summary(lmer1)

#REML criterion at convergence: 229
#(compare these across models - looking for lowest)
#proportion of variance from sites ~2% - so remove site factor I think
#(you get this from the Random Effects section of the model summary - Random Effects - the site_id variance vs the site_id + residual variance
#(ie, of all the variance how much is attributable to site)

par(mfrow = c(1,1))
#Q-Q plots to check residuals
qqnorm(resid(lmer1))
qqline(resid(lmer1))  

#model 2
lm2 <- lm(peak_onset_doy~tmax_spring, data = df1)
summary(lm2)

#Seems to me that a simple linear model on spring tmax makes most sense
#since sites don't make up a lot of the variance, and adding fall_tmin didn't help this model

#MODELS FOR EPC FRUITING ONSET

#model 1, first yes predicted by tmin_spring, with site as a random effect
lmer1 <- lmer(first_yes_doy~tmax_winter + (1| site_id_factor), data = df1)
summary(lmer1)

#REML criterion at convergence: 459
#(compare these across models - looking for lowest)
#proportion of variance from sites ~1/3
#(you get this from the Random Effects section of the model summary - Random Effects - the site_id variance vs the site_id + residual variance
#(ie, of all the variance how much is attributable to site)

par(mfrow = c(1,1))
#Q-Q plots to check residuals
qqnorm(resid(lmer1))
qqline(resid(lmer1))  

#model 2
lmer2 <- lmer(first_yes_doy~tmax_winter + tmin_fall + (1| site_id_factor), data = df1)
summary(lmer2)

#REML criterion at convergence: 451
#proportion of variance from sites ~2/3

#Q-Q plots to check residuals
qqnorm(resid(lmer2))
qqline(resid(lmer2)) 

#is the more complex model worth the extra parameters?
#make the more complex model the first argument
anova(lmer2, lmer1)

#adding fall tmin does  improve the model (P ChiSq 0.032)

#model 3
lmer3 <- lmer(first_yes_doy~tmax_winter + tmin_fall + prcp_summer + (1| site_id_factor), data = df1)
summary(lmer3)

#REML criterion at convergence: 456
#proportion of variance from sites ~1/3

#Q-Q plots to check residuals
qqnorm(resid(lmer3))
qqline(resid(lmer3))

anova(lmer3, lmer2)

#adding winter tmin does not improve the model (0.24)

#model 4
lmer4 <- lmer(first_yes_doy~tmax_winter + tmin_fall + prcp_winter + (1| site_id_factor), data = df1)
summary(lmer4)

#REML criterion at convergence: 457
#proportion of variance from sites ~2/3

#Q-Q plots to check residuals
qqnorm(resid(lmer4))
qqline(resid(lmer4))

anova(lmer4, lmer2)

#adding winter precip instead doesn't improve P of 0.8

#model 5
lmer5 <- lmer(first_yes_doy~tmax_winter + tmin_fall + latitude + (1| site_id_factor), data = df1)
summary(lmer5)

#REML criterion at convergence: 447
#proportion of variance from sites ~2/3

#Q-Q plots to check residuals
qqnorm(resid(lmer5))
qqline(resid(lmer5))

anova(lmer5, lmer2)

#adding latitude doesn't improve 0.3

#MODELS FOR EPC FLOWERING ONSET

#model 1, first yes predicted by tmin_spring, with site as a random effect
#disregarding tmin summer bc happens before
lmer1 <- lmer(first_yes_doy~tmin_spring + (1| site_id_factor), data = df1)
summary(lmer1)

#REML criterion at convergence: 681
#(compare these across models - looking for lowest)
#proportion of variance from sites ~1/9
#(you get this from the Random Effects section of the model summary - Random Effects - the site_id variance vs the site_id + residual variance
#(ie, of all the variance how much is attributable to site)

#Q-Q plots to check residuals
qqnorm(resid(lmer1))
qqline(resid(lmer1))  

#model 2
lmer2 <- lmer(first_yes_doy~tmin_spring + tmin_fall + (1| site_id_factor), data = df1)
summary(lmer2)

#REML criterion at convergence: 679
#proportion of variance from sites ~1/9

#Q-Q plots to check residuals
qqnorm(resid(lmer2))
qqline(resid(lmer2)) 

#is the more complex model worth the extra parameters?
#make the more complex model the first argument
anova(lmer2, lmer1)

#adding fall tmin does not improve the model (P ChiSq 0.8)

#model 3
lmer3 <- lmer(first_yes_doy~tmin_spring + prcp_winter + (1| site_id_factor), data = df1)
summary(lmer3)

#REML criterion at convergence: 688
#proportion of variance from sites ~1/9

#Q-Q plots to check residuals
qqnorm(resid(lmer3))
qqline(resid(lmer3))

anova(lmer3, lmer1)

#adding winter prcp doesn't improve the model p of 0.7

#model 4
lmer4 <- lmer(first_yes_doy~tmin_spring + latitude + (1| site_id_factor), data = df1)
summary(lmer4)

#REML criterion at convergence: 675
#proportion of variance from sites ~1/8

#Q-Q plots to check residuals
qqnorm(resid(lmer4))
qqline(resid(lmer4))

anova(lmer4, lmer1)

#adding latitude doesn't improve, but it is close - Pr(>Chisq) of 0.06


#EPC PEAK FLOWRING

#model 1, first yes predicted by tmin_spring, with site as a random effect
lmer1 <- lmer(peak_onset_doy~tmax_summer + (1| site_id_factor), data = df)
summary(lmer1)

#REML criterion at convergence: 440.6
#(compare these across models - looking for lowest)
#proportion of variance from sites ~1/8 
#(you get this from the Random Effects section of the model summary - Random Effects - the site_id variance vs the site_id + residual variance
#(ie, of all the variance how much is attributable to site)

#Q-Q plots to check residuals
qqnorm(resid(lmer1))
qqline(resid(lmer1))  

#model 2
lmer2 <- lmer(peak_onset_doy~tmax_summer + tmax_spring + (1| site_id_factor), data = df)
summary(lmer2)

#REML criterion at convergence: 674
#proportion of variance from sites ~1/8

#Q-Q plots to check residuals
qqnorm(resid(lmer2))
qqline(resid(lmer2)) 

#is the more complex model worth the extra parameters?
#make the more complex model the first argument
anova(lmer2, lmer1)

#adding spring tmax does  improve the model (P ChiSq 0.01)

#model 3
lmer3 <- lmer(peak_onset_doy~tmax_summer + tmax_spring + prcp_spring + (1| site_id_factor), data = df)
summary(lmer3)

#REML criterion at convergence: 671.1
#proportion of variance from sites ~1/8

#Q-Q plots to check residuals
qqnorm(resid(lmer3))
qqline(resid(lmer3))

anova(lmer3, lmer2)

#adding spring prcp doesn't improve the model

#model 4
lmer4 <- lmer(peak_onset_doy~tmin_summer + tmin_spring + latitude + (1| site_id_factor), data = df)
summary(lmer4)

#REML criterion at convergence: 430.1
#proportion of variance from sites ~1/4

#Q-Q plots to check residuals
qqnorm(resid(lmer4))
qqline(resid(lmer4))

anova(lmer4, lmer2)

#adding latitude doesn't improve


#Wild Bergamont FLOWERING ONSET
#model 1, first yes predicted by prcp fall, with site as a random effect
lmer1 <- lmer(first_yes_doy~prcp_fall + (1| site_id_factor), data = df)
summary(lmer1)

#REML criterion at convergence: 554.9
#(compare these across models - looking for lowest)
#proportion of variance from sites ~1/8 
#(you get this from the Random Effects section of the model summary - Random Effects - the site_id variance vs the site_id + residual variance
#(ie, of all the variance how much is attributable to site)

#Q-Q plots to check residuals
qqnorm(resid(lmer1))
qqline(resid(lmer1))  

#model 2
lmer2 <- lmer(first_yes_doy~prcp_fall + latitude + (1| site_id_factor), data = df)
summary(lmer2)

#REML criterion at convergence: 551.7
#proportion of variance from sites ~1/8

#Q-Q plots to check residuals
qqnorm(resid(lmer2))
qqline(resid(lmer2)) 

#is the more complex model worth the extra parameters?
#make the more complex model the first argument
anova(lmer2, lmer1)

#adding latitude does not improve the model (P ChiSq 0.82)


#WILD BERG DROPPING PRE DAY 200 FOR PEAK ONSET and seeing impact
df1 <- subset(df, species_id == 931 & phenophase_id == 390 & first_yes_doy > 200)

par(mfrow = c(1,1))
hg_onset <- hist(df1$first_yes_doy)
hg_peak_onset <- hist(df1$peak_onset_doy, col = 'red')

plot(hg_onset)
plot(hg_peak_onset, col = 'red', add = TRUE)


#WILD BERG FRUIT PEAK DURATION

#model 1, first yes predicted by X, with site as a random effect
lmer1 <- lmer(peak_duration~prcp_fall + (1| site_id_factor), data = df1)
summary(lmer1)

#REML criterion at convergence: 229
#(compare these across models - looking for lowest)
#proportion of variance from sites - 1/10
#(you get this from the Random Effects section of the model summary - Random Effects - the site_id variance vs the site_id + residual variance
#(ie, of all the variance how much is attributable to site)

par(mfrow = c(1,1))
#Q-Q plots to check residuals
qqnorm(resid(lmer1))
qqline(resid(lmer1))  

#model 2
lmer2 <- lmer(peak_duration~prcp_fall + latitude + (1| site_id_factor), data = df1)
summary(lmer2)

#REML criterion at convergence: 225
#proportion of variance from sites 1/10

#Q-Q plots to check residuals
qqnorm(resid(lmer2))
qqline(resid(lmer2)) 

#is the more complex model worth the extra parameters?
#make the more complex model the first argument
anova(lmer2, lmer1)

#adding latitude does not improve the model (P ChiSq 0.6)


#WILD BERG FRUIT PEAK ONSET

#model 1, first yes predicted by X, with site as a random effect
lmer1 <- lmer(peak_onset_doy~prcp_fall + (1| site_id_factor), data = df1)
summary(lmer1)

#REML criterion at convergence: 218
#(compare these across models - looking for lowest)
#proportion of variance from sites - 4x
#(you get this from the Random Effects section of the model summary - Random Effects - the site_id variance vs the site_id + residual variance
#(ie, of all the variance how much is attributable to site)

par(mfrow = c(1,1))
#Q-Q plots to check residuals
qqnorm(resid(lmer1))
qqline(resid(lmer1))  

#model 2
lmer2 <- lmer(peak_onset_doy~prcp_fall + latitude + (1| site_id_factor), data = df1)
summary(lmer2)

#REML criterion at convergence: 214
#proportion of variance from sites 5x

#Q-Q plots to check residuals
qqnorm(resid(lmer2))
qqline(resid(lmer2)) 

#is the more complex model worth the extra parameters?
#make the more complex model the first argument
anova(lmer2, lmer1)

#adding latitude does not improve the model (P ChiSq 0.9)

#WILD BERG - MODELING RIPE FRUIT ONSET

#model 1, first yes predicted by tmin_spring, with site as a random effect
lmer1 <- lmer(first_yes_doy~tmin_summer + (1| site_id_factor), data = df1)
summary(lmer1)

#REML criterion at convergence: 296
#(compare these across models - looking for lowest)
#proportion of variance from sites ~1/3
#(you get this from the Random Effects section of the model summary - Random Effects - the site_id variance vs the site_id + residual variance
#(ie, of all the variance how much is attributable to site)

par(mfrow = c(1,1))
#Q-Q plots to check residuals
qqnorm(resid(lmer1))
qqline(resid(lmer1))  

#model 2
lmer2 <- lmer(first_yes_doy~tmin_summer + prcp_fall + (1| site_id_factor), data = df1)
summary(lmer2)

#REML criterion at convergence: 298
#proportion of variance from sites 1/3

#Q-Q plots to check residuals
qqnorm(resid(lmer2))
qqline(resid(lmer2)) 

#is the more complex model worth the extra parameters?
#make the more complex model the first argument
anova(lmer2, lmer1)

#adding fall precip does not improve the model (P ChiSq 0.07) - marginal

#model 3
lmer3 <- lmer(first_yes_doy~tmin_summer + latitude + (1| site_id_factor), data = df1)
summary(lmer3)

#REML criterion at convergence: 292
#proportion of variance from sites 1/3

#Q-Q plots to check residuals
qqnorm(resid(lmer3))
qqline(resid(lmer3))

anova(lmer3, lmer1)

#adding latitude  does not improve the model (0.8) 

##WILD BERG FLOWER PEAK DURATION

#model 1, first yes predicted by tmin_spring, with site as a random effect
lmer1 <- lmer(peak_duration~prcp_fall + (1| site_id_factor), data = df)
summary(lmer1)

#REML criterion at convergence: 279.1
#(compare these across models - looking for lowest)
#proportion of variance from sites ~1/8 
#(you get this from the Random Effects section of the model summary - Random Effects - the site_id variance vs the site_id + residual variance
#(ie, of all the variance how much is attributable to site)

#Q-Q plots to check residuals
qqnorm(resid(lmer1))
qqline(resid(lmer1))  

#model 2
lmer2 <- lmer(peak_duration~prcp_fall + tmin_spring + (1| site_id_factor), data = df)
summary(lmer2)

#REML criterion at convergence: 274.1
#proportion of variance from sites ~1/8

#Q-Q plots to check residuals
qqnorm(resid(lmer2))
qqline(resid(lmer2)) 

#is the more complex model worth the extra parameters?
#make the more complex model the first argument
anova(lmer2, lmer1)

#adding spring tmin does  improve the model (P ChiSq 0.02)

#model 3
lmer3 <- lmer(peak_duration~prcp_fall + tmin_spring + latitude + (1| site_id_factor), data = df)
summary(lmer3)

#REML criterion at convergence: 273.1
#proportion of variance from sites ~1/8

#Q-Q plots to check residuals
qqnorm(resid(lmer3))
qqline(resid(lmer3))

anova(lmer3, lmer2)

#adding latitude doesn't improve the model

##Common ButtonBush FLowering Onset

#model 1, first yes predicted by tmin_spring, with site as a random effect
lmer1 <- lmer(first_yes_doy~tmin_spring + (1| site_id_factor), data = df)
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
lmer2 <- lmer(first_yes_doy~tmin_spring + tmin_winter + (1| site_id_factor), data = df)
summary(lmer2)

#REML criterion at convergence: 300
#proportion of variance from sites ~1/4

#Q-Q plots to check residuals
qqnorm(resid(lmer2))
qqline(resid(lmer2)) 

#is the more complex model worth the extra parameters?
#make the more complex model the first argument
anova(lmer2, lmer1)

#adding winter tmin does  improve the model (P ChiSq 0.45)

#model 3
lmer3 <- lmer(first_yes_doy~tmin_spring + latitude + (1| site_id_factor), data = df)
summary(lmer3)

#REML criterion at convergence: 295
#proportion of variance from sites ~1/5

#Q-Q plots to check residuals
qqnorm(resid(lmer3))
qqline(resid(lmer3))

anova(lmer3, lmer1)

#adding latitude does improve the model (0.03)

#model 4
lmer4 <- lmer(first_yes_doy~tmin_spring*latitude + (1| site_id_factor), data = df)
summary(lmer4)

#REML criterion at convergence: 295
#proportion of variance from sites ~1/5

#Q-Q plots to check residuals
qqnorm(resid(lmer4))
qqline(resid(lmer4))

anova(lmer4, lmer3)

#adding latitude as interaction doesn't improve

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

