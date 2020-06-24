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
library(openxlsx)
library(xlsx)

# Core calibration data. The initial Excel spreadsheet is too big to read in using xlsx packages, so the main sheet was exported to a csv first
# This will have to be repeated whenever new years of data or edits are made
raw_data <- 
  
  
  
  
  
  
  
  
  
  
  









