# old stuff

# binning full data set, including damaged individuals ####
# clean up MN object
MN <- MN %>%
  # calculate total count / per m2 for each size class
  group_by(siteID, collectDate, dw) %>%
  summarise(count = sum(no_m2)) %>%
  ungroup() %>%
  # remove counts that are NA
  na.omit() %>%
  # duplicate number of rows based on count
  #i.e. if count = 10, duplicate that row 10 times
  uncount(count)

# bin and center macroinvertebrate dry weight estimates

# first, calculate breaks that are inclusive of data range
breaks <-  2^seq(floor(range(log2(MN$dw))[1]),
                 ceiling(range(log2(MN$dw))[2]))
MN <- MN %>%
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
MN <- MN %>%
  unnest(cols = bins) %>%
  select(-data) %>%
  ungroup()

# full data set ####
# first, separate "collectDate" in MN so it matches degree.days
MN <- MN %>%
  separate(collectDate,
           c("year", "month", "day", "hour", "min", "sec")) %>%
  select(-hour, -min, -sec)
# join tables
dat <- left_join(MN, degree.days)

# save dat for modeling in next script
saveRDS(dat, "data/binned_MN_degree_days_2017_2019.RDS")