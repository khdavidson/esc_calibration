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
library(AICcmodavg)
library(bbmle)
library(MuMIn)
library(nlme)
library(lme4)

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
  filter(!is.na(hp_est)) %>%
  print()

#####################################################################################################################################################

#                                                         DATA EXPLORATION 

# check out histograms for each creek and all creeks together 
ggplot(estu_data, aes(x=hp_est)) +
  geom_histogram(bins=6, fill="light blue", colour="blue") +
  labs(x="Fence count", y="Frequency") +
  facet_wrap(~stream) +
  theme_bw()

ggplot(estu_data, aes(x=hp_est)) +
  geom_histogram(bins=6, fill="light blue", colour="blue") +
  labs(x="Fence count", y="Frequency") +
  theme_bw()

# fence counts are pretty skewed. there are a lot of cases with escapement < 5000, and very few cases over 10,000
# this may affect the model structure, as seen from scott's work where cube-root transformation was required to
# normalize. however, we can avoid normality requirements using GLMMs 

#####################################################################################################################################################

# TEST ASSUMPTIONS 

lm_global <- lm(hp_est ~ lp_est1, data=estu_data) 
lm_global <- lm(hp_est ~ total_carc, data=estu_data)
hist(estu_data$hp_est)

plot(lm_global)
r<-resid(lm_global)
hist(r)
plot(r)
qqnorm(r)
qqline(r)

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

x1 <- 0:10
y1 = dpois(x1, lambda=1)
plot(x1, y1, type="h", main="Poisson with u = 1")
par(op)

# although cube-root transformation improves normality perhaps a bit, it still isn't optimal as it creats a bimodal distribution. 
# i think it is likely best to 






#####################################################################################################################################################

#                                                                MODELS


#-------------Model formulae

# cube-root transformed linear models - these should be made from all systems pooled as it will go towards creating models for the entire system,
# so it should include all the variability
lm0 <- lm(hp_est_cube ~ 1, data=estu_data)
lm1 <- lm(hp_est_cube ~ lp_est1_cube, data=estu_data)
lm2 <- lm(hp_est_cube ~ lp_est1_cube + total_carc_cube, data=estu_data)

# poisson GLM - Decker model formulae but non-linear so no cube-root transform
glm0 <- glm(hp_est ~ 1, data=estu_data, family = poisson(link = "log"))
glm1 <- glm(hp_est ~ lp_est1, data=estu_data, family = poisson(link = "log"))
glm2 <- glm(hp_est ~ lp_est1 + total_carc, data=estu_data, family = poisson(link = "log"))



#-------------Compete models 

AICc(lm0, lm1, lm2, glm0, glm1, glm2)
anova(lm0, lm1, lm2, glm0, glm1, glm2)

# AIC comparison table 
cand.models <- list(lm0, lm1, lm2, glm0, glm1, glm2)
cand.names <- c("lm0_null", "lm1", "lm2", "glm0_null", "glm1", "glm2")

# Table of all AIC vals w weights 
t <- AICctab(lm0, lm1, lm2, glm0, glm1, glm2, nobs=83, logLik=T, base=T, weights=T, delta=T, sort=T)
print(t)

# export table as csv
#class(t) <- "data.frame"
#write.csv(t, "AIC_table_dec72019_propnmods.csv", row.names = T)

# Select top models 95% model weight
top.set <- model.sel(cand.models)
top.comp.models.95 <- get.models(top.set, cumsum(weight)<=0.95)
modavg.95 <- model.avg(top.comp.models.95)
summary(modavg.95)








