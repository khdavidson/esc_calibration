
# stellako timing assessment 
# feb 17,2021


#####################################################################################################################################################

setwd("~/Documents/ANALYSIS/data")

library(tidyverse)
library(xlsx)
library(openxlsx)
library(egg)

stellako.nadina.raw <- read.xlsx("stellako-nadina data overview.xlsx", sheet=2)


#####################################################################################################################################################

#                                                                         CLEAN

data <- stellako.nadina.raw %>% 
  mutate(date=as.Date(date, origin = "1899-12-30")) %>% 
  mutate(yday=lubridate::yday(date)) %>%
  mutate(year=lubridate::year(date)) %>% 
  mutate_at("year", as.factor) %>%
  print()

#####################################################################################################################################################

# SCALES
scales<-ggplot() +
  geom_line(data=data%>%filter(!is.na(daily_perc_ND_scales)), 
    aes(x=as.Date(yday, origin="2005-01-01"), y=daily_perc_ND_scales, group=year), 
      colour="gray70", size=1, alpha=0.8) +
  geom_point(data=data%>%filter(!is.na(daily_perc_ND_scales)), 
    aes(x=as.Date(yday, origin="2005-01-01"), y=daily_perc_ND_scales, group=year), 
      fill="gray70", colour="gray70", stroke=1.5, shape=21, size=3, alpha=0.8) +
  
  geom_line(data=data%>%filter(!is.na(daily_perc_ND_scales)), 
    aes(x=as.Date(yday, origin="2005-01-01"), y=daily_perc_ND_applied, group=year), 
      colour="black", size=1, alpha=0.8) +
  geom_point(data=data%>%filter(!is.na(daily_perc_ND_scales)), 
    aes(x=as.Date(yday, origin="2005-01-01"), y=daily_perc_ND_applied, group=year), 
      colour="black", fill="#00b8ff", stroke=1.5, shape=21, size=3, alpha=0.8) +
  labs(x="", y="Daily % Nadina (at Stellako) - SCALES") +
  theme_bw()   +
  facet_wrap(~year)

# DNA 
dna<-ggplot() +
  geom_line(data=data%>%filter(!is.na(daily_perc_DNA)), 
    aes(x=as.Date(yday, origin="2005-01-01"), y=daily_perc_DNA, group=year, colour=year), 
      size=1, alpha=0.8) +
  geom_point(data=data%>%filter(!is.na(daily_perc_DNA)), 
    aes(x=as.Date(yday, origin="2005-01-01"), y=daily_perc_DNA, group=year, fill=year), 
      shape=21, stroke=1.5, size=3, alpha=0.8) +
  scale_x_date(date_labels="%b %d", limits=c(as.Date(215, origin="2005-01-01"), as.Date(270, origin="2005-01-01"))) +
  labs(x="", y="Daily % Nadina (at Stellako) - DNA") +
  theme_bw()   

# TAG %
tags<-ggplot() +
  geom_line(data=data%>%filter(!is.na(daily_perc_tags)), 
    aes(x=as.Date(yday, origin="2005-01-01"), y=daily_perc_tags, group=year, colour=year), 
      size=1) +
  geom_point(data=data%>%filter(!is.na(daily_perc_tags)), 
    aes(x=as.Date(yday, origin="2005-01-01"), y=daily_perc_tags, group=year, fill=year), 
      shape=21, size=3, stroke=1.5) +
  scale_x_date(date_labels="%b %d", limits=c(as.Date(215, origin="2005-01-01"), as.Date(270, origin="2005-01-01"))) +
  labs(x="", y="Daily % Nadina (at Stellako) - RADIO TAGS") +
  theme_bw()   


# SUBTRACTION
subs<-ggplot() +
  geom_line(data=data%>%filter(!is.na(daily_perc_SN_sub)), 
    aes(x=as.Date(yday, origin="2005-01-01"), y=daily_perc_SN_sub, group=year, colour=year), 
      size=1, alpha=0.8) +
  geom_point(data=data%>%filter(!is.na(daily_perc_SN_sub)), 
    aes(x=as.Date(yday, origin="2005-01-01"), y=daily_perc_SN_sub, group=year, fill=year), 
      shape=21, size=3, stroke=1.5, alpha=0.8) +
  scale_x_date(date_labels="%b %d", limits=c(as.Date(215, origin="2005-01-01"), as.Date(270, origin="2005-01-01"))) +
  labs(x="", y="Daily % Nadina (at Stellako) - SUBTRACTION") +
  theme_bw()   

ggarrange(scales, dna, tags, subs, ncol=2, nrow=2)









