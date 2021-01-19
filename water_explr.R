# discharge data exploration 
# cot 2020
# NOT ON GITHUB YET 

#####################################################################################################################################################

library(tidyverse)

setwd("~/Documents/ANALYSIS/data/calibration")

water.dat <- read.csv("water_data.csv")

#####################################################################################################################################################

#                                                                          EXPLORATION


#######################
# Level vs. discharge #
#######################

ggplot(water.dat %>% filter(gazetted_stream_name != "Harrison River"), aes(x=water_discharge, y=water_level)) + 
  geom_point()

# there is a surprisingly weak relationship between water level and discharge. Harrison River was removed from this visual analysis as it only had 
# one discharge data point which was obtained from an unreliable proxy station; it is much higher in both discharge and level and was significantly 
# driving the relationship, creating an unrealistic relationship. 



















