# custom functions for analysis 
# written by Justin Pomeranz
# jfpomeranz@gmail.com
# started June 2020

# function to add published length weight regression coefficients to macroinvertebrate taxonomic data from NEON sites. Only returns observations which has a corresponding LW coefficient. i.e. output may be smaller than input
LW_coef <- function(x, lw_coef, percent = FALSE){
  # x = inv_taxonomyProcessed table from NEON data product "DP1.20120.001"
  # lw_coef = table of taxon specific length_weight coefficients including higher taxonomic classification. Note that this table must have a column named "taxon" which is used to match up with the different taxonomic levels in x
  # percent = Report percent of rows in x which have a taxon-specific LW equation? default is FALSE
  
  
  # record number of rows in x to calculate percent of observations with a LW equation
  x.nrow <- nrow(x)
  # make indexes for each taxonomic level and extract those rows
  # remove rows from x that are accounted for in taxononmic tables
  
  # genus
  genus_index <- which(x$genus %in% lw_coef$taxon)
  # make table of observations which have a genus-specific LW equation
  genus_table <- x[genus_index,]
  # remove these from global taxonomic data
  x <- x[-genus_index,]
  
  # repeat process for other taxonomic levels
  
  # family
  family_index <- which(x$family %in% lw_coef$taxon)
  family_table <- x[family_index,]
  x <- x[-family_index,]
  
  # order
  order_index <- which(x$order %in% lw_coef$taxon)
  order_table <- x[order_index,]
  x <- x[-order_index,]
  
  # subclass
  subclass_index <- which(x$subclass %in% lw_coef$taxon)
  subclass_table <- x[subclass_index,]
  x <- x[-subclass_index,]
  
  # class
  class_index <- which(x$class %in% lw_coef$taxon)
  class_table <- x[class_index,]
  x <- x[-class_index,]
  
  # combine coefficients to taxonomic tables
  lw_cols <- c("taxon",
               "a",
               "b",
               "L_units",
               "dw_units",
               "formula_type",
               "formula")
  genus_table <- merge(genus_table, lw_coef[,lw_cols], 
                       by.x = "genus",
                       by.y = "taxon", all.x = TRUE)
  family_table <- merge(family_table, lw_coef[,lw_cols], 
                        by.x = "family",
                        by.y = "taxon", all.x = TRUE)
  order_table <- merge(order_table, lw_coef[,lw_cols], 
                       by.x = "order",
                       by.y = "taxon", all.x = TRUE)
  subclass_table <- merge(subclass_table, lw_coef[,lw_cols], 
                          by.x = "subclass",
                          by.y = "taxon", all.x = TRUE)
  class_table <- merge(class_table, lw_coef[,lw_cols], 
                       by.x = "class",
                       by.y = "taxon", all.x = TRUE)
  datout <- rbind(genus_table, family_table, order_table, subclass_table, class_table)
  
  if (percent == TRUE){
    message(paste(nrow(datout)/x.nrow * 100,
                "percent of input data had Length-Weight equations"))
  }
  return(datout)
}

# estimate dry weight (dw) values from length-weight regression coefficients
est_dw <- function(x, fieldData){
  # x = inv_taxonomyProcessed table from NEON data product "DP1.20120.001" with LW coefficients added using the LW_coeff() function
  # fieldData = inv_fieldData table from NEON data product "DP1.20120.001"
  
  # functions uses tidyverse functions
  require(tidyverse)
  
  # simplify fieldData to three columns
  field = fieldData %>%
    select(siteID, sampleID, benthicArea)
  # add benthicArea column from fieldData to x. This is necessary to calculate number per m2 below
  
  # join by siteID and sampleID
  x <- left_join(x, field, by = c("siteID", "sampleID"))
  
  # correct counts to per meter squared. Round to nearest integer
  x <- x %>%
    mutate(no_m2 = round(estimatedTotalCount / benthicArea))
  
  # calculate dw based on different formula types
  # sizeClass = length of individual in mm
  x <- x %>% 
    mutate(dw = case_when(formula_type == 1 ~ a * sizeClass^b,
                          formula_type == 2 ~ exp(a + b * log(sizeClass))))
  return(x)
}

# function modified from Pomeranz, Warburton, and Harding. Anthropogenic mining alters macroinvertebrate size spectra in streams. Freshwater Biology. 2019. 
bin_and_center <- function(data, var, breaks, ...){
  # data is a data frame
  # var is a string, and is the name of a column in data which you want to bin
  # breaks controls the number of bins as defined in hist(). Can be a vector giving breakpoints of bins [i.e. Log10 breaks = 10^seq(min, max)], a single number designating number of bins, etc. See ?hist for details. If not supplied by the user a default of log2 width bins are calculated 
  
  # are breaks supplied by the user?
  if(exists("breaks") == FALSE){
    
    # calculate log2 width bins which are inclusive of the range of data supplied
    breaks = 2^seq(floor(range(log2(data[[var]]))[1]),
                   ceiling(range(log2(data[[var]]))[2]))
    message("breaks not supplied, using log2 width bins")
  }
  
  # bin values using hist()
  binned_hist = hist(data[[var]], 
                     breaks = breaks,
                     include.lowest = TRUE, plot = FALSE)
  # calculate "left" and "right" edge of bins
  breaks_orig = binned_hist$breaks[1:(length(breaks)-1)]
  breaks_offset = binned_hist$breaks[2:length(breaks)]
  # total bin width = right edge - left edge
  break_width = breaks_offset - breaks_orig
  count = binned_hist$counts 
  dataout = data.frame(
    # normalize counts =count / width (White et al 1997)
    log_count_corrected = log10(count / break_width),
    # original midpoint of bin log10 transformed
    log_mids = log10(binned_hist$mids),
    log_mids_center = NA)
  # remove bins with 0 counts
  # -Inf comes from log10(count / break_width) above
  dataout = dataout[dataout$log_count_corrected !=-Inf,]
  # recenter data at x=0
  mid_row = ceiling(nrow(dataout)/2)
  # subtract value of mid row from all mids
  dataout$log_mids_center = 
    dataout[,"log_mids"] - dataout[mid_row,"log_mids"]
  dataout
}

deg_days <- function(x){
  # x = TSW_30min object from NEON surface water temperature data product dpID = "DP1.20053.001"
  
  require(tidyverse)
  # separate startDateTime into components
  x <- x %>%
    as_tibble() %>%
    # first, need to filter out missing data
    filter(!is.na(surfWaterTempMean)) %>%
    # select only columns that we will need below
    select(siteID, startDateTime, surfWaterTempMean) %>%
    separate(startDateTime,
             c("year", "month", "day", "hour", "min", "sec")) %>%
    # group by individual days
    group_by(siteID, year, month, day) %>%
    # caluculate daily mean temperature
    summarize(daily = mean(surfWaterTempMean)) %>%
    ungroup() %>% # ungroup for next step
    # only want the positive degree days
    # i.e. a negative temperature day doesn't "remove" degree days
    filter(daily >=0) %>%
      # group by site ID and year
    group_by(siteID, year) %>%
      # calculate cumulutaive sum of mean daily temperatures
    mutate(degree_days = cumsum(daily))
    return(x)
}

# plot prior and posterior distributions
# this function is currently only for Intercept and b prioirs
plot_prior_vs_posterior <- function(model,
                                    prior_intercept_mu,
                                    prior_intercept_sigma,
                                    prior_b_mu,
                                    prior_b_sigma){
  require(dplyr)
  message("this function is currently programed for 3 beta coefficients")
  # int_coef <- fixef(model)[1,1:2]
  # slope_coef <- fixef(model)[2,1:2]
  
  # get distributions of coefficients
  fixed_effects <- row.names(fixef(model))
  prior.list = NULL
  for (i in seq_along(fixed_effects)){
    prior.list[[i]] =
      data.frame(value =
                   rnorm(1000, fixef(model)[i,1], fixef(model)[i,2]),
                 sample = fixed_effects[i])
  }
  # sample prior distributions
  post.dat <- data.frame(
       value = c(rnorm(1000, prior_intercept_mu, prior_intercept_sigma),
                 rnorm(1000, prior_b_mu, prior_b_sigma)),
       sample = rep(c("Intercept prior", "b priors"), each = 1000))
  # combine prior and posts for plot
  plot.dat <- bind_rows(prior.list, post.dat)
  # plot
  ggplot(plot.dat, aes(x = value, fill = sample))+
    geom_density(alpha = 0.5) +
    theme_bw() +
    scale_fill_manual(values = c("mediumorchid4", "mediumorchid1",
                                 "turquoise4", "steelblue1", "cadetblue1",
                                 "skyblue1"),
                      breaks = c("Intercept prior", "Intercept",
                                 "b priors", "log_mids_center",
                                 "log10degree_days",
                                 "log_mids_center:log10degree_days"))
}
