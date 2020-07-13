# 1 - Download and wrangle Neon Data

library(neonUtilities)
library(tidyverse)

source("scripts/custom_functions.R")

# load Length Weight coefficient table (used in part C below)
coeff <- read.csv("data/LW_coeffs.csv")


# A) Download data ---------------------------------------------------------

# these are large files, and may take a while to fully download dependent on web traffic and connection speed. 

# download all macroinvertebrate colllection data from January 2017 to december 2019
macro <- loadByProduct(dpID = "DP1.20120.001",
                       startdate = "2017-01",
                       enddate = "2019-12",
                       check.size = FALSE,
                       nCores = 3)
#saveRDS(macro, "data/Benth_macro_download.RDS")
# macro <- readRDS("data/Benth_macro_download.RDS")

# download all surface water temperature data from January 2017 to December 2019
# note - as of 7/7, this is no longer downloading all of the same data. Used to have ~90 site:season:year samples with > 200 daily temperature recordings, but now only have ~8. 
# I emailed NEON on 7/9, have yet to hear back. 
water.t <- loadByProduct(dpID = "DP1.20053.001",
                         startdate = "2017-01",
                         enddate = "2019-12",
                         avg = "30",  
                         check.size = FALSE,
                         nCores = 3)
#saveRDS(water.t, "data/water_T_download.RDS")
#water.t <- readRDS("data/water_T_download.RDS")
# B) Wrangle Temperature data ------------------------------------------

# calculate cumulative degree days for each site and year
degree.days <- deg_days(x = water.t$TSW_30min)

# add a column with yearly observations by site for model sensitivity analysis later
degree.days <- degree.days %>%
  group_by(siteID, year) %>%
  mutate(yearly.n = n()) %>%
  ungroup()

degree.days %>%
  select(siteID, year, yearly.n) %>%
  filter(yearly.n >200, year != 2019) %>%
  unique()
# this should result in ~ 19 sites and 60 unique collections in 2017/2018

degree.days %>%
  select(siteID, year, yearly.n) %>%
  filter(yearly.n >200, year == 2019) %>%
  unique()
# this should result in ~ 31 unique collections in 2019

# C) Wrangle macroinvertebrate data ------------------------------------

# add length weight coefficients by taxon
MN.lw <- LW_coef(x = macro$inv_taxonomyProcessed,
              lw_coef = coeff,
              percent = TRUE)
# 95.8% of observations have length weight equations

# # estimate dry weights based on length weight coefficients
# MN <- est_dw(MN.lw, fieldData = macro$inv_fieldData)

# questionable measurements ####
# filter out individuals that were "damaged" and measurement was affected
MN.no.damage <- MN.lw %>%
  filter(!str_detect(sampleCondition, "measurement")) %>%
  est_dw( fieldData = macro$inv_fieldData)

# what percent of data is damaged?
message(paste(
  round(nrow(MN.no.damage) / nrow(MN.lw), 4)*100,
  "percent of individuals were not damaged"))


# need to redo once temp data is fixed???? ####
# remove individuals smaller than 0.00195 mg dry mass
# this corresonds to individuals in bins smaller than 2^-9 

# MN.no.damage <- MN.no.damage %>%
#   filter(dw >= 0.00195)

# total number of measured individuals
length(na.omit(MN.no.damage$sizeClass))
#87.4k individual measurements

# clean up MN.no.damage object
MN.no.damage <- MN.no.damage %>%
  # calculate total count / per m2 for each size class
  group_by(siteID, collectDate, dw) %>%
  summarise(count = sum(no_m2)) %>%
  ungroup() %>%
  # remove counts that are NA
  na.omit() %>%
  # duplicate number of rows based on count
  #i.e. if count = 10, duplicate that row 10 times
  uncount(count)

nrow(MN.no.damage)
# 21 million rows

# bin and center macroinvertebrate dry weight estimates

# first, calculate breaks that are inclusive of data range
breaks <-  2^seq(floor(range(log2(MN.no.damage$dw))[1]),
                 ceiling(range(log2(MN.no.damage$dw))[2]))
# log2 bins from -9 to 14
length(breaks)
# 24 log2 width bins

MN.no.damage <- MN.no.damage %>%
  # group combos of site and collectdata
  group_by(siteID, collectDate) %>% 
  # create a list column of data
  nest() %>% 
  # Create another list column with binned dw for each site and collection independently
  mutate(bins = map(data, 
                    bin_and_center,
                    "dw",
                    breaks = breaks)) 


# remove binned data from list column to "normal" data frame format and get rid of "data" list column
MN.no.damage <- MN.no.damage %>%
  unnest(cols = bins) %>%
  select(-data) %>%
  ungroup()

# D) join binned and temperature data -------------------------------------

# no damage data set ####
MN.no.damage <- MN.no.damage %>%
  separate(collectDate,
           c("year", "month", "day", "hour", "min", "sec")) %>%
  select(-hour, -min, -sec)
# join tables
dat.filtered <- left_join(MN.no.damage, degree.days)

# save dat for modeling in next script
saveRDS(dat.filtered, "data/binned_no_damage_MN_degree_days_2017_2019.RDS")
