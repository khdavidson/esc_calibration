# calibration work 
# june 2020
# k davidson based on early work on EStu calibration by S. Decker 
# relates to previous script: Early_Stuart_calibration_expanded.R written by S. Decker 
# there were also some scripts produced by M. Hawkshaw

# Reading in data is one of the first issues with the previous script. It called to 'esc_data_expanded.txt', the source of which is not known, which
# has issues for reproducibility/transparency. 

# In general, there are a lot of R scripts floating around relating to this work. I will be recreating all from scratch, beginning with the base 
# spreadsheet (originally compiled by P. Welch and K. Benner) and filtering/cleaning as needed. 


#####################################################################################################################################################

#                                                                          PREAMBLE

# Objective: to calibrate escapement estimates for Fraser sockeye using high precision and low precision paired data. From approx. 2010-2018, a
# calibration program was funded whereby streams would have both high-precision (fence or M-R, or in later years, SONAR) projects and low-precision
# (visual counts/carcass recoveries) surveys conducted independently. The point was to obtain estimates of the "real" (HP) and "realized" (LP)
# escapement estimates to get an index (or calibration factor). 
# The index would be used to convert LP estimates to HP estimates in years where only LP projects happened. 

# Historical methods were simply to take the ratio of the  HP estimate/LP estimate each year, and average those ratios into one index for the creek.

# More advanced and in-depth work by S. Decker and M. Hawkshaw showed that this was not adequate, and that a more rigorous regression model 
# accounting for more factors was needed. (see 'Early_Stuart_esc_estimation_issudes.docx' for details). 

# This script picks up at the multiple regression model approach. It focuses on the Early Stuart sub-sample to start, as that has the best sample
# size and data quality. 

# Eventually this should include environmental data, possibly a non-quantitative grouping association (eg, NMDS) or PCA to cluster like-streams to
# apply calibration factors to similar systems across the watershed (i.e., those without HP and LP paired data).

#####################################################################################################################################################

setwd("~/Documents/ANALYSIS/Data/calibration")

library(tidyverse)
library(AICcmodavg)   # for AICc()
library(bbmle)        # for AICctab()
library(MuMIn)        # for model.avg(), model.sel(), get.models()
#library(nlme)        # if GLMM was fit (was not)
#library(lme4)        # if GLMM was fit (was not)
library(car)          # for vif()

# will eventually need core calibration data
# note: The Excel spreadsheet is too big to read in using xlsx packages, so the main sheet will need to be exported to a csv first
# This will have to be repeated whenever new years of data or edits are made

# for now load estu subset created by S. Decker 
raw_estu_data <- read.csv("esc_data_expanded.csv")

#####################################################################################################################################################

#                                                         DATA CLEANING

# clean and replicate calculations that gave last 2 columns in original csv
# remove cases where there is no fence estimate (=NA)
# cubed variables for linear models as per Decker previous analysis 
estu_data <- raw_estu_data %>% 
  rename(stream = Stream, 
    year = Year,
    total_carc = TotalCarc,
    hp_est = FenceEsc,
    max_live = MaxLive,
    lp_est1 = MaxLiveDead1, 
    lp_est2 = MaxLiveDead2, 
    ratio_dead_live = DeadvsLive,
    ratio_dead_livedead = DeadvsLiveAndDead) %>% 
  mutate(ratio_dead_live = total_carc/max_live,
    ratio_dead_livedead = total_carc/lp_est1) %>%
  mutate(total_carc_cube = total_carc^(1/3),
    hp_est_cube = hp_est^(1/3), 
    max_live_cube = max_live^(1/3),
    lp_est1_cube = lp_est1^(1/3),
    lp_est2_cube = lp_est2^(1/3)) %>%
  mutate(usid = paste0(substr(stream, 1, 3),"-", year)) %>% 
  filter(!is.na(hp_est)) %>%
  print()

#####################################################################################################################################################

#                                                         DATA EXPLORATION 

# check out histograms for each creek and all creeks together 
ggplot(estu_data, aes(x=hp_est)) +
  geom_histogram(bins=20, fill="light blue", colour="blue") +
  labs(x="Fence count", y="Frequency") +
  facet_wrap(~stream) +
  theme_bw()

ggplot(estu_data, aes(x=hp_est)) +
  geom_histogram(bins=20, fill="light blue", colour="blue") +
  labs(x="Fence count", y="Frequency") +
  theme_bw()

# fence counts are pretty skewed. there are a lot of cases with escapement < 5000, and very few cases over 10,000
# this may affect the model structure, as seen from scott's work where cube-root transformation was required to
# normalize. however, we can avoid normality requirements using GLMMs too.

#####################################################################################################################################################

#                                                          TEST ASSUMPTIONS 

####################
# DATA EXPLORATION #
####################

hist(estu_data$hp_est)

# plot the data and simulated distributions
# observed data: 
op <- par(mfrow = c(2, 2))
hist(estu_data$hp_est, nclass=10, xlab="HP estimate", main="Observed data")
hist((estu_data$hp_est)^(1/3), nclass=10, xlab="HP estimate cubed", main="Observed data cubed")

# simulated normal data:
Y <- rnorm(1281, mean=mean(estu_data$hp_est), sd=sd(estu_data$hp_est))
hist(Y, nclass=10, main="Simulated normal data", xlab="HP estimate")
#X <- seq(from=0, to=30, length=200)
#Y <- dnorm(X, mean=mean(estu_data$hp_est), sd=sd(estu_data$hp_est))
#plot(X, Y, type="l", xlab="HP estimate", ylab="Probabilities", ylim=c(0,0.25), xlim=c(0, 30), main="Normal density curve")

# poisson
x1 <- 0:10
y1 = dpois(x1, lambda=1)
plot(x1, y1, type="h", main="Poisson with u = 1")
par(op)

#### Raw data appear poisson-distributed, so will have to do detailed tests for linear model assumptions prior to proceeding with linear models


##########################
# MODEL ASSUMPTION TESTS #
##########################

#--------- Linearity
ggplot(estu_data, aes(x=lp_est1, y=hp_est)) +
  geom_point() +
  geom_smooth(method="lm")
# linear

ggplot(estu_data, aes(x=total_carc, y=hp_est)) +
  geom_point() +
  geom_smooth(method="lm")
# a little bit of a curvi-linear relationship here 




#--------- Normal distribution
lm1 <- lm(hp_est ~ lp_est1, data=estu_data)
plot(lm1)
  # residuals vs fitted: obvious cone-shaped spread pattern to residuals 
  # qqplot: very mid-peaked
  # scale-location: cone-shaped spread again
  # residuals vs leverage: clustered, one point outside 0.5 line
r1<-resid(lm1)
hist(r1)
plot(r1)
qqnorm(r1)
qqline(r1)


lm2 <- lm(hp_est ~ lp_est1 + total_carc, data=estu_data)
plot(lm2)
  # residuals vs fitted: obvious cone-shaped spread pattern to residuals 
  # qqplot: mid-peaked less intense than lm1 but still high
  # scale-location: cone-shaped spread again
  # residuals vs leverage: clustered, one point outside 0.5 line
r2<-resid(lm2)
hist(r2)
plot(r2)
qqnorm(r2)
qqline(r2)

# shows clear violation of normal distribution in both cases


###
# LOG RESPONSE    *does not improve lm1 or lm2*
###
lm1.log <- lm(log(hp_est) ~ lp_est1, data=estu_data)
plot(lm1.log)
  # residuals vs fitted: obvious peak spread pattern to residuals 
  # qqplot: skewed
  # scale-location: not bad
  # residuals vs leverage: clustered, point outside 0.5 and 1 lines 
r1.log<-resid(lm1.log)
hist(r1.log)
plot(r1.log)
qqnorm(r1.log)
qqline(r1.log)

lm2.log <- lm(log(hp_est) ~ lp_est1 + total_carc, data=estu_data)
plot(lm2.log)
  # residuals vs fitted: obvious peak spread pattern to residuals 
  # qqplot: skewed
  # scale-location: not bad
  # residuals vs leverage: clustered, point outside 0.5 and 1 lines 
r2.log<-resid(lm2.log)
hist(r2.log)
plot(r2.log)
qqnorm(r2.log)
qqline(r2.log)

###
# CUBE-ROOT RESPONSE   *better than above* 
###
lm1.cube <- lm(hp_est_cube ~ lp_est1, data=estu_data)
plot(lm1.cube)
  # residuals vs fitted: obvious peak spread pattern to residuals 
  # qqplot: little skewed, perhaps best so far
  # scale-location: not bad
  # residuals vs leverage: clustered, point outside 0.5 and 1 lines 
r1.cube<-resid(lm1.cube)
hist(r1.cube)
plot(r1.cube)
qqnorm(r1.cube)
qqline(r1.cube)

lm2.cube <- lm(hp_est_cube ~ lp_est1 + total_carc, data=estu_data)
plot(lm2.cube)
  # residuals vs fitted: obvious peak spread pattern to residuals 
  # qqplot: little skewed, perhaps best so far
  # scale-location: not bad
  # residuals vs leverage: clustered, point outside 0.5 and 1 lines 
r2.cube<-resid(lm2.cube)
hist(r2.cube)
plot(r2.cube)
qqnorm(r2.cube)
qqline(r2.cube)


# if linear model is used, cube-root transformed response should be used (Decker observation), but this may change the relationship b/w response 
# and predictors.



#--------- Collinearity of predictors 

# pairs plot
z <- cbind(estu_data$hp_est, estu_data$lp_est1, estu_data$total_carc)
colnames(z) <- c("HP estimate","LP estimate", "total carcs")
cor(z)
pairs(z)

panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
{
  usr = par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r = cor(x, y,use="na.or.complete")
  txt = format(c(r, 0.123456789), digits=digits)[1]
  txt = paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex.cor = 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = 1.5)#cex.cor * r)
}
pairs(z,
      upper.panel = panel.cor,
      cex=1,
      pch=16)

# VIF
vif(lm2)


# SUMMARY: LP_estimate and total_carcasses are highly correlated (r2=0.95), and vif = 9.7  -->  too high

# NOTE: With interaction lp_est*total_carc, VIF scores worsened for total_carc which makes sense as it adds structural multicollinearity 
lm. <- lm(hp_est ~ lp_est1 + total_carc + lp_est1*total_carc, data=estu_data)
plot(lm3)
vif(lm3)


###
# CENTER - typically for equalizing effect sizes between vastly different predictors, or for structural collinearity, but could help
###
estu_data <- estu_data %>% 
  mutate(lp_est1_crs = lp_est1-mean(lp_est1),
    total_carc_crs = total_carc-mean(total_carc)) %>% 
  print()

# pairs plot
z.c <- cbind(estu_data$hp_est, estu_data$lp_est1_crs, estu_data$total_carc_crs)
colnames(z.c) <- c("HP estimate","LP estimate centered", "total carcs centered")
cor(z.c)
pairs(z.c)

panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
{
  usr = par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r = cor(x, y,use="na.or.complete")
  txt = format(c(r, 0.123456789), digits=digits)[1]
  txt = paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex.cor = 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = 1.5)#cex.cor * r)
}
pairs(z.c,
      upper.panel = panel.cor,
      cex=1,
      pch=16)

# VIF
lm.center <- lm(hp_est ~ lp_est1_crs + total_carc_crs, data=estu_data)
vif(lm.center)

# No change. Still correlated.  



###
# Assess models separately 
###
lm1 <- lm(hp_est ~ lp_est1, data=estu_data)
lm2 <- lm(hp_est ~ lp_est1 + total_carc, data=estu_data)
lm3 <- lm(hp_est ~ total_carc, data=estu_data)
summary(lm4)   # adj r2 = 0.95
summary(lm1)   # adj r2 = 0.92


#--------- Summary of model exploration:
# - Relationships are linear 
# - Data are non-normal: appear Poisson-distributed, or require cube-root transformation to the response variable 
# - Predictors are highly correlated: although Decker showed improved model fit when total_carc was added, this is likely because over-fitting
# with redundant highly correlated variables occurred (VIF=9.7 and R2=0.95). In fact, preliminary model fitting to re-create analysis had model fit
# issues when came to fitting hp_est~lp_est1+total_carcs; model weight=1.0 and <0.00001 for other models with only 1 predictor is suspicious and
# suggests poor model fitting/over-fitting to higher paramaterized, multicollinear model. Centering did nothing to improve VIF or R2 as multicollinearity
# is data not structural. 
# - Univariate analysis shows that total_carcs is a slightly better predictor than the lp_est1 
##--------- 

##--------- Conclusions moving forward: 
# Some options to move forward with this analysis exist: 
# - Explore Poisson-distributed GLM (GLMM?) 
# - Ensure predictive model predicts within the range of values (i.e., model data set encompasses a range of potential escapements)


#####################################################################################################################################################

#                                                             MODEL TRAINING/TESTING

# Here we select 25% of the data to hold aside in order to train and test the model fits. 

#############
# SUBSAMPLE #
#############
# Subsample and training dataframes 

# Select a random 75% of the data (47 observations*0.75 = 35)
set.seed(100)
estu_sub_data <- estu_data[sample(1:nrow(estu_data), 35, replace=F), ]

# Extract the 25% that aren't having the model fitted to them - guinea pig dataframe 
estu_train_data <- anti_join(estu_data, estu_sub_data, by="usid")



#######################
# FIT TRAINING MODELS #
#######################
# We know non-normal data so just skip lm()s as they won't be appropriate

glm0 <- glm(hp_est ~ 1, data=estu_sub_data, family = poisson(link = "log"))

glm1 <- glm(hp_est ~ lp_est1, data=estu_sub_data, family = poisson(link = "log"))
summary(glm1)

glm2 <- glm(hp_est ~ lp_est1 + total_carc, data=estu_sub_data, family = poisson(link = "log"))
summary(glm2)

glm3 <- glm(hp_est ~ total_carc, data=estu_sub_data, family = poisson(link = "log"))
summary(glm3)


###########
# PREDICT #
###########

# glm1
estu_train_data$predict.hpe.1 <- predict.glm(glm1, estu_train_data, type="response")

ggplot(estu_train_data) +
  geom_bar(aes(x=year, y=hp_est), stat="identity", fill="gray40", colour="black", alpha=0.5) + 
  geom_bar(aes(x=year, y=predict.hpe.1), stat="identity", fill="light blue", colour="blue", alpha=0.5) +
  scale_y_continuous(limits=c(0,23000)) +
  theme_bw() +
  facet_wrap(.~stream)


# glm2
estu_train_data$predict.hpe.2 <- predict.glm(glm2, estu_train_data, type="response")

ggplot(estu_train_data) +
  geom_bar(aes(x=year, y=hp_est), stat="identity", fill="gray40", colour="black", alpha=0.5) + 
  geom_bar(aes(x=year, y=predict.hpe.2), stat="identity", fill="light blue", colour="blue", alpha=0.5) +
  scale_y_continuous(limits=c(0,23000)) +
  theme_bw() +
  facet_wrap(.~stream)

# glm3
estu_train_data$predict.hpe.3 <- predict.glm(glm3, estu_train_data, type="response")

ggplot(estu_train_data) +
  geom_bar(aes(x=year, y=hp_est), position="dodge", stat="identity", fill="gray40", colour="black", alpha=0.5) + 
  geom_bar(aes(x=year, y=predict.hpe.3), position="dodge", stat="identity", fill="light blue", colour="blue", alpha=0.5) +
  scale_y_continuous(limits=c(0,23000)) +
  theme_bw() +
  facet_wrap(.~stream)


#---------- DBEs
# reshape for plotting and calculate DBEs
estu_train_data <- estu_train_data %>% 
  gather(predicted.glm, predicted.vals, c(predict.hpe.1:predict.hpe.3)) %>% 
  mutate(diff = hp_est-predicted.vals) %>%
  print()

ggplot(estu_train_data, aes(x=year, y=diff, fill=predicted.glm, colour=predicted.glm)) +
  geom_bar(stat="identity", alpha=0.5, position="dodge", width=0.5) +
  geom_hline(yintercept=0) +
  geom_text(data=estu_train_data %>% filter(predicted.glm=="predict.hpe.3"), aes(label=hp_est), colour="black", vjust=2) +
  scale_x_continuous(breaks=seq(1990,2010,2)) +
  labs(x="") +
  theme_bw() +
  facet_wrap(.~stream)

# mean DBEs - average of predicted estimates from 3 model-fitted predictions compared to known escapements 
train_sum <- estu_train_data %>%
  group_by(stream, hp_est) %>% 
  summarize(mean.pred=mean(predicted.vals), sd.pred=sd(predicted.vals)) %>%
  print()

ggplot(train_sum) +
  geom_bar(aes(x=stream, y=hp_est, group=hp_est), fill="white", colour="gray80", stat="identity", position="dodge") +
  geom_bar(aes(x=stream, y=mean.pred, group=hp_est), fill="gray60", colour="black", stat="identity", position="dodge", alpha=0.5) +
  geom_errorbar(aes(x=stream, ymin=mean.pred-sd.pred, ymax=mean.pred+sd.pred, group=hp_est), width=0, position=position_dodge(width=0.9)) +
  #scale_fill_distiller(palette="Spectral") +
  #geom_text(aes(x=stream, y=mean.pred, group=hp_est, label=hp_est), colour="black", vjust=-3, position=position_dodge(width=0.9)) +
  theme_bw()


##################
# COMPARE MODELS #
##################
AICc(glm0, glm1, glm2, glm3)

# AIC comparison table 
cand.models <- list(glm0, glm1, glm2, glm3)
cand.names <- c("glm0_null", "glm1", "glm2", "glm3")

# Table of all AIC vals w weights 
t <- AICctab(glm0, glm1, glm2, glm3, nobs=35, logLik=T, base=T, weights=T, delta=T, sort=T)
print(t)

# can't model average because 100% of the weights are in one model 

###############
## WHAT HAVE WE LEARNED FROM THIS? 
# These models have terrible predictive power.

#####################################################################################################################################################

#                                                                MODELS

#################
# LINEAR MODELS #
#################

#-------------Model structure and fit

# cube-root transformed linear models - these should be made from all systems pooled as it will go towards creating models for the entire system,
# so it should include all the variability
# while Decker cube-root transformed predictors and response, linear model assumptions only apply to the response variable so i don't think the
# predictors should have been cube-root transformed???? 

# structure: 

# high precision fence count ~ (max live + cumulative dead up to max live date) + total # carcasses

lm0 <- lm(hp_est_cube ~ 1, data=estu_data)
lm1 <- lm(hp_est_cube ~ lp_est1_cube, data=estu_data)
lm2 <- lm(hp_est_cube ~ lp_est1_cube + total_carc_cube, data=estu_data)
lm3 <- lm(hp_est_cube ~ lp_est1, data=estu_data)
lm4 <- lm(hp_est_cube ~ lp_est1 + total_carc, data=estu_data)

#-------------Compete models 

AICc(lm0, lm1, lm2, lm3, lm4)
anova(lm0, lm1, lm2, lm3, lm4)

# AIC comparison table 
cand.models <- list(lm0, lm1, lm2, lm3, lm4)
cand.names <- c("lm0_null", "lm1", "lm2", "lm3", "lm4")

# Table of all AIC vals w weights 
t <- AICctab(lm0, lm1, lm2, lm3, lm4, nobs=83, logLik=T, base=T, weights=T, delta=T, sort=T)
t_2 <- AICctab(lm0, lm3, lm4, nobs=83, logLik=T, base=T, weights=T, delta=T, sort=T)
print(t)

# export table as csv
#class(t) <- "data.frame"
#write.csv(t, "AIC_table_dec72019_propnmods.csv", row.names = T)

# Select top models 95% model weight
top.set <- model.sel(cand.models)
top.comp.models.95 <- get.models(top.set, cumsum(weight)<=0.95)
modavg.95 <- model.avg(top.comp.models.95)
summary(modavg.95)





################
# POISSON GLMS #
################

#-------------Model structure and fit

# poisson GLM - Decker model formulae but non-linear so no cube-root transform
glm0 <- glm(hp_est ~ 1, data=estu_data, family = poisson(link = "log"))
glm1 <- glm(hp_est ~ lp_est1, data=estu_data, family = poisson(link = "log"))
glm2 <- glm(hp_est ~ lp_est1 + total_carc, data=estu_data, family = poisson(link = "log"))
glm3 <- glm(hp_est ~ lp_est1 + total_carc + lp_est1*total_carc, data=estu_data, family = poisson(link = "log"))

# should this be GLMM with random effect for year and/or stream ? 



#-------------Compete models 

AICc(glm0, glm1, glm2)
anova(glm0, glm1, glm2)

# AIC comparison table 
cand.models <- list(glm0, glm1, glm2, glm3)
cand.names <- c("glm0_null", "glm1", "glm2", "glm3")

# Table of all AIC vals w weights 
t <- AICctab(glm0, glm1, glm2, glm3, nobs=83, logLik=T, base=T, weights=T, delta=T, sort=T)
print(t)

# export table as csv
#class(t) <- "data.frame"
#write.csv(t, "AIC_table_dec72019_propnmods.csv", row.names = T)

# Select top models 95% model weight
top.set <- model.sel(cand.models)
top.comp.models.95 <- get.models(top.set, cumsum(weight)<=0.95)
modavg.95 <- model.avg(top.comp.models.95)
summary(modavg.95)


######################
# COMPETE ALL MODELS #   Bad mojo
######################

# AIC comparison table 
cand.models <- list(lm0, lm3, lm4, glm0, glm1, glm2, glm3)
cand.names <- c("lm0_null", "lm3", "lm4", "glm0_null", "glm1", "glm2", "glm3")

# Table of all AIC vals w weights 
t <- AICctab(lm0, lm3, lm4, glm0, glm1, glm2, glm3, nobs=83, logLik=T, base=T, weights=T, delta=T, sort=T)
print(t)











