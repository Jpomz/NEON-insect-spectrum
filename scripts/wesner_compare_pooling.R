# libraries
library(tidyverse)
library(brms)
library(viridis)
library(ggridges)
library(gridExtra)
library(janitor)
library(tidyselect)
library(tidybayes)
library(nlme)
library(ggdist)
library(cowplot)

# load all NEON data
dat_all <- readRDS(file = "data/dat_all.rds")
data_conditioned_on <- dat_all %>% filter(year != 2019) %>%
  filter(!is.na(degree_days), 
         yearly.n >= 200)
data_new <- dat_all %>% filter(year == 2019)
dat <- data_conditioned_on

# load model
mod <- readRDS("data/dd_as_r_slope_6-30-20.RDS") # r slopes after talking with JSW 6/30

# load posterior summaries of partial pooling model (i.e. the main model in the manuscript)
sitedd.all <- readRDS(file = "data/sitedd.all.rds")

# Compare pooling ---------------------------------------------------------
# Model with no pooling (separate regressions at each site)
dat_group <- dat %>% mutate(group = paste0(siteID,"_", round(degree_days, 0),"_", year))
dat_lists <- split(dat_group, f = dat_group$group)

# mod_nopool <- brm_multiple(log_count_corrected ~ log_mids_center,
#                   family = gaussian(),
#                   data = dat_lists,
#                   prior = c(prior(normal(-1, 1),
#                                   class = "b", coef = "log_mids_center"),
#                             prior(normal(0,1), class = "b"),
#                             # prior(normal(4.5, 1.5), class = "Intercept"),
#                             prior(cauchy(0,1), class = "sigma")),
#                   chains = 2,
#                   iter = 2000)

# saveRDS(mod_nopool, file = "data/mod_nopool.rds")
mod_nopool <- readRDS(file = "data/mod_nopool.rds")
# mod_nopool <- update(mod_nopool, newdata = dat_lists, iter = 2000) #ignore warnings of "sampling not done" and of "rhats". Seems like a bug with brm_multiple

#sites to add to no pooling model
sites_rep <- dat_group %>% select(group) %>% distinct(group) %>% slice(rep(1:n(), each=2000))

# summarize no pooling estimates
post_het <- posterior_samples(mod_nopool, summary = F) %>% as_tibble() %>% clean_names() %>% mutate(iter = 1:nrow(.),
                                                                                                  group = sites_rep$group) %>%
  separate(group, c("siteID", "degree_days", "year")) %>% 
  mutate(degree_days = as.numeric(degree_days)) %>% 
  select(-lp) %>%
  mutate(Slope = b_log_mids_center,
         y_fit = b_intercept + 
           b_log_mids_center*mean(dat$log_mids_center),
         Intercept = rnorm(nrow(.), y_fit, sigma)) %>%
  mutate(prediction_level = "Sites in original model",
         method = "No pooling") %>% 
  pivot_longer(cols = c(Slope, Intercept), names_to = "Parameter") %>% 
  glimpse()

post_het_new <- posterior_samples(mod_nopool, summary = F) %>% as_tibble() %>% clean_names() %>% mutate(iter = 1:nrow(.),
                                                                                                        group = sites_rep$group) %>%
  separate(group, c("siteID", "degree_days", "year")) %>% 
  mutate(degree_days = as.numeric(degree_days)) %>% 
  select(-lp) %>%
  mutate(Slope = b_log_mids_center,
         y_fit = b_intercept + 
           b_log_mids_center*mean(dat_group$log_mids_center),
         Intercept = rnorm(nrow(.), y_fit, sigma)) %>%
  mutate(prediction_level = "Sites not in original model",
         method = "No pooling") %>% 
  pivot_longer(cols = c(Slope, Intercept), names_to = "Parameter") %>% 
  glimpse()


post_het_all <- bind_rows(post_het, post_het_new) 

pool_nopool <- sitedd.all %>% 
  mutate(method = "Partial Pooling") %>% 
  bind_rows(post_het_all) %>% 
  filter(prediction_level == "Sites in original model") %>% 
  group_by(method, degree_days, Parameter) %>% 
  median_qi(value) %>% 
  ggplot(aes(x = degree_days, y = value, ymin = .lower, ymax = .upper)) +
  geom_pointrange(position = position_dodge(width = 0.2)) +
  facet_grid(Parameter~method, scales = "free", switch = "y") +
  scale_x_log10() +
  labs(y = "Parameter value") +
  guides(color = F) +
  theme_bw()

pool_nopool
ggsave(pool_nopool, file = "plots/pool_nopool.png", dpi = 600, width = 7, height = 5)





# Estimate 2019 size spectra (no pooling model) ---------------------------

dat_2019 <- dat_all %>% filter(year == 2019)
dat2019_lists <- split(dat_2019, f = dat_2019$siteID)
# 
# mod_nopool2019 <- brm_multiple(log_count_corrected ~ log_mids_center,
#                                              family = gaussian(),
#                                              data = dat2019_lists,
#                                              prior = c(prior(normal(-2, 1),
#                                                              class = "b", coef = "log_mids_center"),
#                                                        prior(normal(0,1), class = "b"),
#                                                        # prior(normal(4.5, 1.5), class = "Intercept"),
#                                                        prior(cauchy(0, 0.1), class = "sigma")),
#                                              chains = 1,
#                                              iter = 2000)
# 
# mod_nopool2019
# 
#saveRDS(mod_nopool2019, file = "data/mod_nopool2019.rds")
mod_nopool2019 <- readRDS(file = "data/mod_nopool2019.rds")

#sites to add to no pooling model
sites_rep2019 <- dat_2019 %>% select(siteID) %>% distinct(siteID) %>% slice(rep(1:n(), each=1000))

# degree days to add
dd2019 <- dat_2019 %>% select(degree_days) %>% distinct()

# summarize no pooling estimates for 2019
post_het2019 <- posterior_samples(mod_nopool2019, summary = F) %>% as_tibble() %>% 
  clean_names() %>% mutate(iter = 1:nrow(.),
                           siteID = sites_rep2019$siteID) %>% 
  select(-sigma, -lp) %>%
  left_join(dat_2019 %>% select(siteID, degree_days) %>% distinct(siteID, degree_days)) %>% 
  mutate(Slope = b_log_mids_center,
         Intercept = b_intercept + 
           b_log_mids_center*mean(dat$log_mids_center)) %>%
  mutate(prediction_level = "New data from 2019",
         method = "No pooling 2019") %>% 
  pivot_longer(cols = c(Slope, Intercept), names_to = "Parameter") %>% 
  glimpse()


partial_predict <- posterior_samples(mod) %>% clean_names() %>% mutate(iter = 1:nrow(.)) %>% 
  select(!contains("cor_site")) %>% 
  pivot_longer(contains("r_site"), names_to = "ranef", values_to = "offset") %>%
  mutate(siteID = toupper(str_sub(ranef, 11,14)),
         ranef = paste0("r_",str_sub(ranef, 16))) %>%
  pivot_wider(names_from = ranef, values_from = offset) %>%
  right_join(dat_2019 %>% select(siteID, degree_days) %>% distinct(siteID, degree_days)) %>%
  mutate(Slope = b_log_mids_center + r_log_mids_center + (b_log_mids_center_log10degree_days + 
                                                            r_log10degree_days_log_mids_center)*log10(degree_days),
         y_fit = b_intercept + r_intercept + rnorm(nrow(.), 0, sd_site_id_intercept) + 
           (b_log_mids_center + r_log_mids_center)*mean(dat$log_mids_center) +
           (b_log10degree_days + r_log10degree_days)*log10(degree_days) + 
           (b_log_mids_center_log10degree_days + r_log10degree_days_log_mids_center)*log10(degree_days)*mean(dat$log_mids_center),
         Intercept = rnorm(nrow(.), y_fit, sigma)) %>%
  mutate(prediction_level = "Sites in original model",
         method = "Partial Pooling") %>% 
  pivot_longer(cols = c(Slope, Intercept), names_to = "Parameter") %>% 
  filter(!is.na(value)) %>% 
  glimpse()

nopool_predict <- posterior_samples(mod_nopool, summary = F) %>% as_tibble() %>% clean_names() %>% mutate(iter = 1:nrow(.),
                                                                                        group = sites_rep$group) %>%
  separate(group, c("siteID", "degree_days", "year")) %>% 
  select(-degree_days) %>% 
  select(-lp) %>%
  right_join(data_new %>% select(siteID, degree_days) %>% distinct(siteID, degree_days)) %>%
  mutate(Slope = b_log_mids_center,
         y_fit = b_intercept + 
           b_log_mids_center*mean(dat$log_mids_center),
         Intercept = rnorm(nrow(.), y_fit, sigma)) %>%
  mutate(prediction_level = "Sites in original model",
         method = "No pooling") %>% 
  pivot_longer(cols = c(Slope, Intercept), names_to = "Parameter") %>% 
  glimpse()

post_hetall2019 <- bind_rows(nopool_predict, post_het2019, partial_predict)

post_hetall2019_medians <- post_hetall2019 %>% 
  filter(!is.na(value)) %>% 
  # filter(prediction_level != "Sites in original model") %>%
  group_by(method, degree_days, Parameter, prediction_level) %>% 
  median_qi(value)

post_hetall_error <- post_hetall2019 %>%
  select(method, degree_days, Parameter, method, value, siteID) %>% 
  group_by(method, degree_days, Parameter) %>% 
  summarize(median = median(value)) %>%
  pivot_wider(names_from = c(Parameter, method), values_from = median) %>%
  clean_names() %>%
  mutate(slope_pp_error = abs(slope_partial_pooling - slope_no_pooling_2019),
         slope_np_error = abs(slope_no_pooling - slope_no_pooling_2019),
         intercept_pp_error = abs(intercept_partial_pooling - intercept_no_pooling_2019),
         intercept_np_error = abs(intercept_no_pooling - intercept_no_pooling_2019)) %>% 
  pivot_longer(cols = contains("error")) %>%
  separate(name, c("Parameter", "model", "delete")) %>% 
  rename(abs_error = value) %>% 
  glimpse()
  
post_hetall2019_medians %>% 
ggplot(aes(x = degree_days, y = value, ymin = .lower, ymax = .upper, color = method)) +
  # geom_pointrange(position = position_dodge(width = 0.2)) +
  geom_pointrange() +
  facet_grid(Parameter~., scales = "free", switch = "y") +
  scale_x_log10() +
  labs(y = "Parameter value") +
  # guides(color = F) +
  theme_bw()

post_hetall_error %>%  
  group_by(Parameter, model) %>% 
  arrange(abs_error) %>% 
  mutate(rank = 1:31,
         rank_nudge = as.numeric(case_when(model == "pp" ~ rank + .5,
                                TRUE ~ rank - .5))) %>%
  filter(!is.na(abs_error)) %>% 
  ggplot(aes(x = rank_nudge, y = abs_error, fill = model)) + 
  geom_bar(stat = "identity", position = "dodge") +
  # geom_segment(aes(yend = abs_error, xend = rank_nudge, y = 0)) +
  facet_grid(Parameter~., scales = "free") +
  coord_flip() +
  scale_fill_brewer(type = "qual")




