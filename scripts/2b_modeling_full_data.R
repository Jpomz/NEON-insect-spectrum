# 2b modeling full binned data

# libraries
library(tidyverse)
library(brms)

# load data
dat <- readRDS("data/binned_MN_degree_days_2017_2019.RDS")
dat <- ungroup(dat)

# filter out 2017-2018 data
dat <- dat %>%
  filter(year != 2019)

# combine date columns into one collect date factor
dat <- dat %>%
  unite(collectDate, year:day, remove = FALSE)

# Subset data which has degree_days values
dat.with.dd <- dat %>%
  filter(!is.na(degree_days))

# run models with / without yearly.n < X
# Subset data with more than 200 daily observations per year
dat.200.days <- dat.with.dd %>%
  filter(yearly.n >= 200)

ggplot(dat.200.days,
       aes(y = log_count_corrected, 
           x = log_mids_center,
           color = as.factor(log10(degree_days)))) +
  geom_point() +
  facet_wrap(.~siteID) +
  theme(legend.position = "none")

ggplot(dat.200.days,
       aes(y = log_count_corrected, 
           x = log_mids_center,
           color = siteID)) +
  geom_point() +
  facet_wrap(.~year) +
  theme(legend.position = "none")

# Bayesian models ---------------------------------------------------------

# define my priors
my_priors <- c(prior(normal(-1, 0.5), class = "b"),
               prior(normal(4.5, 1), class = "Intercept"),
               prior(cauchy(0, 0.1), class = "sigma"))
# For models that also include sd priors
my_priors_sd <- c(my_priors, prior(cauchy(0,0.1), class = "sd"))


# random intercept by site models
# model log_count_corrected ~ log_mids_center * log10(degree_days) + (1|siteID) + (1|year)

ri3 <- brm(data = dat.with.dd,
           log_count_corrected ~
             log_mids_center * log10(degree_days) + (1|siteID), 
           family = gaussian(),
           prior = my_priors_sd,
           chains = 1, 
           iter = 1000)
ri3_200 <- update(ri3,
                  newdata = dat.200.days)

ri4 <- update(ri3,
              formula. = ~. + (1|year),
              newdata = dat.with.dd,
              iter = 4000)
ri5 <- update(ri4,
              newdata = dat.200.days)
