# Sensitivity analysis of priors

# justin pomeranz July 6 2020

# libraries
library(tidyverse)
library(brms)
library(viridis)
library(ggridges)
library(gridExtra)
source("scripts/custom_functions.R")

# load data
# this data set has been filtered to remove the individuals which had damage that affected their measurements
dat <- readRDS("data/binned_no_damage_MN_degree_days_2017_2019.RDS")

# # combine date columns into one collect date factor
# dat <- dat %>%
#   unite(collectDate, year:day, remove = FALSE)

# Subset data which has degree_days values and more than 200 daily observations per year
dat <- dat %>%
  filter(!is.na(degree_days), 
         yearly.n >= 200)

# filter out 2019 data for prediction later
dat.2019 <- dat %>%
  filter(year == 2019)

# filter out 2017-2018 data for bayesian model
dat <- dat %>%
  filter(year != 2019)


my_priors <- c()

# Bayesian Model
mod <- brm(data = dat,
           log_count_corrected ~
             log_mids_center * log10(degree_days) +
             (1+ log10(degree_days)*log_mids_center|siteID) + (1|year), 
           family = gaussian(),
           prior = c(prior(normal(-1, 1),
                           class = "b", coef = "log_mids_center"),
                     prior(normal(0,1), class = "b"),
                     prior(normal(4.5, 1.5), class = "Intercept"),
                     prior(cauchy(0, 0.1), class = "sigma"),
                     prior(cauchy(0,0.1), class = "sd")),
           chains = 4, 
           iter = 4000,
           cores = 4,
           control = list(adapt_delta = 0.99,
                          max_treedepth = 12))
mod_sdx2 <- brm(data = dat,
           log_count_corrected ~
             log_mids_center * log10(degree_days) +
             (1+ log10(degree_days)*log_mids_center|siteID) + (1|year), 
           family = gaussian(),
           prior = c(prior(normal(-1, 2), class = "b",
                         coef = "log_mids_center"),
           prior(normal(0,2), class = "b"),
           prior(normal(4.5, 3.0), class = "Intercept"),
           prior(cauchy(0, 0.2), class = "sigma"),
           prior(cauchy(0,0.2), class = "sd")),
           chains = 4, 
           iter = 4000,
           cores = 4,
           control = list(adapt_delta = 0.99,
                          max_treedepth = 12))

mod_mu_p_1 <- brm(data = dat,
                log_count_corrected ~
                  log_mids_center * log10(degree_days) +
                  (1+ log10(degree_days)*log_mids_center|siteID) +
                  (1|year), 
                family = gaussian(),
                prior = c(prior(normal(0, 1), class = "b",
                                coef = "log_mids_center"),
                          prior(normal(1,1), class = "b"),
                          prior(normal(5.5, 1.5), class = "Intercept"),
                          prior(cauchy(1, 0.1), class = "sigma"),
                          prior(cauchy(1,0.1), class = "sd")),
                chains = 4, 
                iter = 4000,
                cores = 4,
                control = list(adapt_delta = 0.99,
                               max_treedepth = 12))

mod_mu_m_1 <- brm(data = dat,
                  log_count_corrected ~
                    log_mids_center * log10(degree_days) +
                    (1+ log10(degree_days)*log_mids_center|siteID) +
                    (1|year), 
                  family = gaussian(),
                  prior = c(prior(normal(-2, 1), class = "b",
                                  coef = "log_mids_center"),
                            prior(normal(-1,1), class = "b"),
                            prior(normal(3.5, 1.5), class = "Intercept"),
                            prior(cauchy(-1, 0.1), class = "sigma"),
                            prior(cauchy(-1,0.1), class = "sd")),
                  chains = 4, 
                  iter = 4000,
                  cores = 4,
                  control = list(adapt_delta = 0.99,
                                 max_treedepth = 12))

data.frame(value = c(
  # Int_post = 
  rnorm(1000, fixef(mod)[1,1], fixef(mod)[1,2]),
  # Int_post = 
  rnorm(1000, fixef(mod_sdx2)[1,1], fixef(mod_sdx2)[1,2]),
  # Int_post = 
  rnorm(1000, fixef(mod_mu_p_1)[1,1], fixef(mod_mu_p_1)[1,2]),
  # Int_post = 
  rnorm(1000, fixef(mod_mu_m_1)[1,1], fixef(mod_mu_m_1)[1,2]),
  # M_post =
  rnorm(1000, fixef(mod)[2,1], fixef(mod)[2,2]),
  # M_post =
  rnorm(1000, fixef(mod_sdx2)[2,1], fixef(mod_sdx2)[2,2]),
  # M_post =
  rnorm(1000, fixef(mod_mu_p_1)[2,1], fixef(mod_mu_p_1)[2,2]),
  # M_post =
  rnorm(1000, fixef(mod_mu_m_1)[2,1], fixef(mod_mu_m_1)[2,2]),
  # DD_post = 
  rnorm(1000, fixef(mod)[3,1], fixef(mod)[3,2]),
  # DD_post = 
  rnorm(1000, fixef(mod_sdx2)[3,1], fixef(mod_sdx2)[3,2]),
  # DD_post = 
  rnorm(1000, fixef(mod_mu_p_1)[3,1], fixef(mod_mu_p_1)[3,2]),
  # DD_post = 
  rnorm(1000, fixef(mod_mu_m_1)[3,1], fixef(mod_mu_m_1)[3,2]),
  # M_DD_post =
  rnorm(1000, fixef(mod)[4,1], fixef(mod)[4,2]),
  # M_DD_post =
  rnorm(1000, fixef(mod_sdx2)[4,1], fixef(mod_sdx2)[4,2]),
  # M_DD_post =
  rnorm(1000, fixef(mod_mu_p_1)[4,1], fixef(mod_mu_p_1)[4,2]),
  # M_DD_post =
  rnorm(1000, fixef(mod_mu_m_1)[4,1], fixef(mod_mu_m_1)[4,2])),
  b_name = rep(c("Intercept", "Intercept","Intercept", "Intercept",
                 "M", "M","M", "M",
                 "Degree_days", "Degree_days","Degree_days", "Degree_days",
                 "M:DD", "M:DD", "M:DD", "M:DD"),
               each = 1000),
  sample = rep(c("model","sd x 2","mu + 1", "mu - 1",
                 "model","sd x 2","mu + 1", "mu - 1",
                 "model","sd x 2","mu + 1", "mu - 1",
                 "model","sd x 2","mu + 1", "mu - 1"), each = 1000)) %>%
  ggplot() +
  geom_density(aes(x = value,
                   fill = sample),
               alpha = 0.5) +
  facet_wrap(.~b_name,
             #cols = 2,
             scales = "free") +
  theme_bw() +
  labs(title =
         "Sensitivity analysis of Prior distributions on posterior",
       subtitle = "Posterior parameter distributions")
ggsave("ms/SI_prior_sensitivity.png")
